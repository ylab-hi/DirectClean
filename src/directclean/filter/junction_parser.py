"""
Junction parser for chimeric (split-mapped) reads.

Parses BAM records to reconstruct the linear layout of supplementary
alignment segments along a read, then identifies inter-segment
junctions that represent potential fusion or artifact breakpoints.

Key concepts
------------
* **Segment**: one continuous alignment of part of a read to the genome,
  described by a CIGAR string.  A chimeric read has ≥ 2 segments.
* **Junction**: the boundary between two adjacent segments on the read.
  This is where the read "jumps" from one genomic locus to another.

The module intentionally does NO filtering or classification — it only
extracts structural information.  Downstream modules (homopolymer.py,
artifact_classifier.py) decide whether a junction is an artifact.
"""

from __future__ import annotations

import re
import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import pysam

logger = logging.getLogger(__name__)

# Pre-compiled regex for parsing CIGAR strings from the SA tag
# Matches pairs like ('562', 'S'), ('28', 'M'), etc.
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SegmentInfo:
    """One alignment segment of a chimeric read.

    Attributes:
        chrom:        Reference chromosome / contig name.
        ref_start:    0-based leftmost mapping position on the reference.
        strand:       '+' or '-'.
        mapq:         Mapping quality.
        read_start:   0-based start on the *original* read sequence.
        read_end:     0-based exclusive end on the original read sequence.
        cigar_string: Raw CIGAR string (for debugging / downstream use).
    """
    chrom: str
    ref_start: int
    strand: str
    mapq: int
    read_start: int
    read_end: int
    cigar_string: str


@dataclass(frozen=True)
class JunctionInfo:
    """A junction between two adjacent segments on a read.

    Attributes:
        read_position:  Position on the read where the junction occurs.
        upstream_seq:   Sequence upstream of the junction (within window).
        downstream_seq: Sequence downstream of the junction (within window).
        left_segment:   The segment on the 5' side (lower read coord).
        right_segment:  The segment on the 3' side (higher read coord).
    """
    read_position: int
    upstream_seq: str
    downstream_seq: str
    left_segment: SegmentInfo
    right_segment: SegmentInfo


@dataclass
class ChimericRead:
    """All structural information for one chimeric read.

    Attributes:
        read_id:   Read name.
        sequence:  Full original read sequence.
        segments:  Ordered list of alignment segments (by read coordinate).
        junctions: Derived list of inter-segment junctions.
    """
    read_id: str
    sequence: str
    segments: List[SegmentInfo] = field(default_factory=list)
    junctions: List[JunctionInfo] = field(default_factory=list)


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

# pysam CIGAR operation codes
_CONSUMES_QUERY = {
    0,  # M  (alignment match)
    1,  # I  (insertion to reference)
    4,  # S  (soft clip)
    7,  # =  (sequence match)
    8,  # X  (sequence mismatch)
}

_CIGAR_OP_MAP = {
    "M": 0, "I": 1, "D": 2, "N": 3,
    "S": 4, "H": 5, "P": 6, "=": 7, "X": 8,
}


def _cigar_to_tuples(cigar_str: str) -> List[Tuple[int, int]]:
    """Convert a CIGAR string to a list of (operation, length) tuples.

    Matches the pysam convention: operation is an int code.

    Examples:
        >>> _cigar_to_tuples("562S28M1D4M")
        [(4, 562), (0, 28), (2, 1), (0, 4)]
    """
    return [
        (_CIGAR_OP_MAP[op], int(length))
        for length, op in _CIGAR_RE.findall(cigar_str)
    ]


def cigar_read_span(
    cigar_tuples: List[Tuple[int, int]],
    is_reverse: bool = False,
) -> Tuple[int, int]:
    """Compute the read-coordinate span [start, end) from CIGAR tuples.

    **Strand handling**: when a supplementary alignment maps to the
    reverse strand, minimap2 writes the CIGAR against the reverse-
    complemented read.  The leading soft-clip in the CIGAR therefore
    corresponds to the *trailing* end of the original (forward) read.
    Setting ``is_reverse=True`` flips the coordinates so that the
    returned interval is always on the **original forward-strand read**.

    Args:
        cigar_tuples: List of (op_code, length) pairs.
        is_reverse:   True if this segment maps to the reverse strand.

    Returns:
        (read_start, read_end) — 0-based half-open interval on the
        original (forward) read sequence.

    Examples:
        >>> # Plus strand: 562S 28M 1D 4M 1069S
        >>> tuples = [(4,562), (0,28), (2,1), (0,4), (4,1069)]
        >>> cigar_read_span(tuples, is_reverse=False)
        (562, 594)

        >>> # Minus strand: 2170S 534M 25S  (total query = 2729)
        >>> # Forward coords: start = 25, end = 25+534 = 559
        >>> tuples = [(4,2170), (0,534), (4,25)]
        >>> cigar_read_span(tuples, is_reverse=True)
        (25, 559)
    """
    # --- compute leading clip, aligned bases, trailing clip ---
    leading_clip = 0
    if cigar_tuples and cigar_tuples[0][0] == 4:  # S
        leading_clip = cigar_tuples[0][1]

    trailing_clip = 0
    if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 4:  # S
        trailing_clip = cigar_tuples[-1][1]

    consumed = 0
    for i, (op, length) in enumerate(cigar_tuples):
        if i == 0 and op == 4:
            continue
        if i == len(cigar_tuples) - 1 and op == 4:
            continue
        if op == 5:  # H — hard clip, no query bases
            continue
        if op in _CONSUMES_QUERY:
            consumed += length

    if not is_reverse:
        # Forward strand: leading clip is the read-start offset
        return leading_clip, leading_clip + consumed
    else:
        # Reverse strand: the CIGAR is written against the RC read.
        # On the original forward read the "leading" clip in the CIGAR
        # is actually at the 3' end, so the true start = trailing_clip.
        forward_start = trailing_clip
        return forward_start, forward_start + consumed


# ---------------------------------------------------------------------------
# SA tag parser
# ---------------------------------------------------------------------------

def _parse_sa_tag(sa_string: str) -> List[dict]:
    """Parse the SA:Z auxiliary tag into a list of segment dicts.

    SA format: ``rname,pos,strand,CIGAR,mapQ,NM;``  (semicolon-separated,
    with a trailing semicolon).

    Returns:
        List of dicts with keys: chrom, pos (1-based), strand, cigar, mapq.
    """
    segments = []
    for entry in sa_string.rstrip(";").split(";"):
        parts = entry.split(",")
        if len(parts) < 5:
            continue
        segments.append({
            "chrom": parts[0],
            "pos": int(parts[1]),      # 1-based in SA tag
            "strand": parts[2],
            "cigar": parts[3],
            "mapq": int(parts[4]),
        })
    return segments


# ---------------------------------------------------------------------------
# Main public API
# ---------------------------------------------------------------------------

def parse_chimeric_read(
    alignment: pysam.AlignedSegment,
    window_size: int = 30,
) -> Optional[ChimericRead]:
    """Parse a BAM alignment into a ChimericRead with junction info.

    Only processes reads that carry an SA (Supplementary Alignment)
    tag — i.e. chimeric / split-mapped reads.  Returns ``None`` for
    non-chimeric reads.

    Steps:
        1. Collect the primary segment and all SA-tag segments.
        2. Compute each segment's read-coordinate span from its CIGAR.
        3. Sort segments by read_start to recover linear order along
           the read.
        4. Identify junctions between adjacent segments and extract
           the surrounding sequence window.

    Args:
        alignment:   A pysam.AlignedSegment (primary or supplementary).
        window_size: Bases to extract on each side of every junction.

    Returns:
        ChimericRead object, or None if the read is not chimeric.
    """
    # --- guard: skip unmapped or non-chimeric reads ---
    if alignment.is_unmapped:
        return None

    sa_value = alignment.get_tag("SA") if alignment.has_tag("SA") else None
    if sa_value is None:
        return None

    # We only want to process each read once.  Use the primary alignment
    # as the canonical entry point; skip supplementary records.
    if alignment.is_supplementary:
        return None

    read_id = alignment.query_name
    sequence = alignment.get_forward_sequence()
    if sequence is None:
        return None

    print(f"[DEBUG] read_id={read_id}")
    print(f"[DEBUG] is_reverse={alignment.is_reverse}")
    print(f"[DEBUG] seq_head={sequence[:80]}")
    print(f"[DEBUG] seq_tail={sequence[-80:]}")

    # ---- 1. Build segment list ----
    segments: List[SegmentInfo] = []

    # Primary segment
    primary_cigar = alignment.cigartuples
    if primary_cigar is None:
        return None

    p_strand = "-" if alignment.is_reverse else "+"
    p_start, p_end = cigar_read_span(primary_cigar, is_reverse=alignment.is_reverse)
    segments.append(SegmentInfo(
        chrom=alignment.reference_name,
        ref_start=alignment.reference_start,
        strand=p_strand,
        mapq=alignment.mapping_quality,
        read_start=p_start,
        read_end=p_end,
        cigar_string=alignment.cigarstring,
    ))

    # Supplementary segments from SA tag
    for sa_seg in _parse_sa_tag(sa_value):
        cigar_tuples = _cigar_to_tuples(sa_seg["cigar"])
        sa_is_reverse = sa_seg["strand"] == "-"
        s_start, s_end = cigar_read_span(cigar_tuples, is_reverse=sa_is_reverse)
        segments.append(SegmentInfo(
            chrom=sa_seg["chrom"],
            ref_start=sa_seg["pos"] - 1,   # SA is 1-based → 0-based
            strand=sa_seg["strand"],
            mapq=sa_seg["mapq"],
            read_start=s_start,
            read_end=s_end,
            cigar_string=sa_seg["cigar"],
        ))

    # Need at least 2 segments for a junction
    if len(segments) < 2:
        return None

    # ---- 2. Sort by read coordinate ----
    segments.sort(key=lambda s: s.read_start)

    # ---- 3. Identify junctions ----
    junctions: List[JunctionInfo] = []
    for i in range(len(segments) - 1):
        left = segments[i]
        right = segments[i + 1]

        # Junction position: midpoint of the gap / overlap between segments
        junction_pos = (left.read_end + right.read_start) // 2

        # Clamp to valid range
        junction_pos = max(0, min(junction_pos, len(sequence)))

        # Extract context window
        up_start = max(0, junction_pos - window_size)
        dn_end = min(len(sequence), junction_pos + window_size)
        upstream_seq = sequence[up_start:junction_pos]
        downstream_seq = sequence[junction_pos:dn_end]

        junctions.append(JunctionInfo(
            read_position=junction_pos,
            upstream_seq=upstream_seq,
            downstream_seq=downstream_seq,
            left_segment=left,
            right_segment=right,
        ))

    return ChimericRead(
        read_id=read_id,
        sequence=sequence,
        segments=segments,
        junctions=junctions,
    )


def iter_chimeric_reads(
    bam_path: str,
    window_size: int = 30,
    min_mapq: int = 0,
) -> "Iterator[ChimericRead]":
    """Iterate over chimeric reads in a BAM file.

    Yields one ChimericRead per read that has an SA tag.  Reads are
    deduplicated by processing only primary alignments.

    Args:
        bam_path:    Path to a coordinate-sorted, indexed BAM file.
        window_size: Context window for junction extraction.
        min_mapq:    Skip reads with mapping quality below this.

    Yields:
        ChimericRead objects.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch():
            # Skip secondary; supplementary handled inside parse_chimeric_read
            if aln.is_secondary:
                continue
            if aln.mapping_quality < min_mapq:
                continue

            chimeric = parse_chimeric_read(aln, window_size=window_size)
            if chimeric is not None:
                yield chimeric
