"""
Unknowns rescuer — rescue reads from Restrander's unknowns output.

Restrander discards reads it cannot orient (unknowns, RTP-RTP and
TSO-TSO artefacts) into a separate FASTQ.  Some of these reads are
actually two cDNA molecules ligated together through internal adapters.

This module:
    1. Scans unknowns FASTQ for internal adapter junctions (reusing
       Stage 3's AdapterFinder).
    2. Chops reads at detected junctions to produce sub-reads.
    3. Determines strand orientation of each sub-read using polyA/T
       tails and TSO/RTP primer signals.
    4. Reverse-complements reverse-strand sub-reads so all output
       reads are in 5'→3' mRNA orientation.
    5. Writes oriented sub-reads to output FASTQ for merging into
       the main pipeline.

Only sub-reads that can be confidently oriented are kept; reads that
remain ambiguous are discarded.

Typical usage (within the pipeline)::

    rescuer = UnknownsRescuer(
        unknowns_fastq=restrander_unknowns,
        output_fastq=oriented_output,
        config=adapter_config,
    )
    report = rescuer.run()
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import edlib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from directclean.rescuer.adaptor_seq import (
    AdapterConfig,
    TSO_SEQUENCE,
    RTP_SEQUENCE,
)
from directclean.rescuer.adapter_finder import AdapterFinder, InternalJunction
from directclean.utils.io import read_fastq, write_fastq
from directclean.utils.sequence_operator import reverse_complement

logger = logging.getLogger(__name__)

# Primer sequences for orientation detection
_TSO = TSO_SEQUENCE.upper()
_RTP = RTP_SEQUENCE.upper()
_TSO_RC = reverse_complement(_TSO)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


@dataclass
class UnknownsRescueReport:
    """Summary statistics for the unknowns rescue operation.

    Attributes:
        total_unknowns:     Total reads in the unknowns FASTQ.
        reads_with_adapter: Reads with ≥1 internal adapter junction.
        reads_without:      Reads with no internal adapter (discarded).
        segments_produced:  Total sub-reads after chopping.
        segments_discarded_short: Sub-reads too short (<min_segment_length).
        oriented_forward:   Sub-reads oriented as forward (no flip).
        oriented_reverse:   Sub-reads oriented as reverse (flipped).
        oriented_unknown:   Sub-reads that could not be oriented (discarded).
        output_reads:       Final reads written to output.
    """

    total_unknowns: int = 0
    reads_with_adapter: int = 0
    reads_without: int = 0
    segments_produced: int = 0
    segments_discarded_short: int = 0
    oriented_forward: int = 0
    oriented_reverse: int = 0
    oriented_unknown: int = 0
    output_reads: int = 0

    def __str__(self) -> str:
        pct_adapter = (
            f"{self.reads_with_adapter / self.total_unknowns * 100:.1f}%"
            if self.total_unknowns > 0
            else "N/A"
        )
        return (
            "=== DirectClean Unknowns Rescue Report ===\n"
            f"  Total unknowns scanned  : {self.total_unknowns:,}\n"
            f"  Reads with internal adpt: {self.reads_with_adapter:,} ({pct_adapter})\n"
            f"  Reads without (skipped) : {self.reads_without:,}\n"
            f"  ---\n"
            f"  Sub-reads produced      : {self.segments_produced:,}\n"
            f"  Discarded (too short)   : {self.segments_discarded_short:,}\n"
            f"  Oriented forward        : {self.oriented_forward:,}\n"
            f"  Oriented reverse (RC)   : {self.oriented_reverse:,}\n"
            f"  Orientation unknown     : {self.oriented_unknown:,}\n"
            f"  ---\n"
            f"  Output reads            : {self.output_reads:,}\n"
            "==========================================="
        )


# ---------------------------------------------------------------------------
# Orientation logic
# ---------------------------------------------------------------------------


def _fuzzy_match(seq: str, query: str, max_ed: int = 3) -> bool:
    """Check if query exists in seq within edit distance."""
    result = edlib.align(query, seq.upper(), mode="HW", k=max_ed)
    return result["editDistance"] != -1


def _has_polya(seq: str, min_run: int = 10) -> bool:
    """Check for polyA run in the last 200bp of sequence."""
    tail = seq[-200:].upper()
    max_run = 0
    current = 0
    for ch in tail:
        if ch == "A":
            current += 1
            if current > max_run:
                max_run = current
        else:
            current = 0
    return max_run >= min_run


def _has_polyt(seq: str, min_run: int = 10) -> bool:
    """Check for polyT run in the first 200bp of sequence."""
    head = seq[:200].upper()
    max_run = 0
    current = 0
    for ch in head:
        if ch == "T":
            current += 1
            if current > max_run:
                max_run = current
        else:
            current = 0
    return max_run >= min_run


def orient_subread(seq: str) -> str | None:
    """Determine strand orientation and return correctly oriented sequence.

    Detection logic (in priority order):
        1. TSO at 5' end OR polyA at 3' end → forward → return as-is.
        2. RTP at 5' end OR polyT at 5' end → reverse → return RC.
        3. TSO_rc at 3' end → reverse → return RC.
        4. None of the above → unknown → return None.

    Args:
        seq: Raw sub-read sequence.

    Returns:
        Oriented sequence (5'→3'), or None if orientation unknown.
    """
    head = seq[:200]
    tail = seq[-200:]

    # Forward signals: TSO at 5' or polyA at 3'
    if _fuzzy_match(head, _TSO) or _has_polya(seq):
        return seq

    # Reverse signals: RTP or polyT at 5'
    if _fuzzy_match(head, _RTP) or _has_polyt(seq):
        return reverse_complement(seq)

    # Reverse signal: TSO_rc at 3'
    if _fuzzy_match(tail, _TSO_RC):
        return reverse_complement(seq)

    return None


# ---------------------------------------------------------------------------
# Main rescuer class
# ---------------------------------------------------------------------------


class UnknownsRescuer:
    """Rescue oriented sub-reads from Restrander's unknowns FASTQ.

    Workflow:
        1. Scan each read for internal adapter junctions.
        2. Chop reads with detected junctions.
        3. Orient each sub-read using primer/tail signals.
        4. Write successfully oriented sub-reads to output.

    Args:
        unknowns_fastq: Path to Restrander's unknowns FASTQ.
        output_fastq:   Path for oriented rescued reads.
        config:         AdapterConfig for internal adapter detection.
        min_confidence: Minimum confidence to chop (default 2).
        min_segment_length: Minimum sub-read length to keep (default 50).
    """

    def __init__(
        self,
        unknowns_fastq: str | Path,
        output_fastq: str | Path,
        config: AdapterConfig | None = None,
        min_confidence: int = 2,
        min_segment_length: int = 50,
    ) -> None:
        self.unknowns_fastq = Path(unknowns_fastq)
        self.output_fastq = Path(output_fastq)
        self.config = config or AdapterConfig()
        self.min_confidence = min_confidence
        self.min_segment_length = min_segment_length

        self.output_fastq.parent.mkdir(parents=True, exist_ok=True)

    def run(self) -> UnknownsRescueReport:
        """Execute the unknowns rescue pipeline.

        Returns:
            UnknownsRescueReport with statistics.
        """
        logger.info("Unknowns Rescue: scanning for internal adapters ...")

        finder = AdapterFinder(self.config)
        report = UnknownsRescueReport()
        output_records: list[SeqRecord] = []

        for record in read_fastq(self.unknowns_fastq):
            report.total_unknowns += 1
            seq = str(record.seq)
            quals = record.letter_annotations.get("phred_quality", [])

            # Step 1: detect internal adapters
            finder_result = finder.find(
                read_id=record.id,
                sequence=seq,
            )
            qualified = [
                j
                for j in finder_result.junctions
                if j.confidence >= self.min_confidence
            ]

            if not qualified:
                report.reads_without += 1
                continue

            report.reads_with_adapter += 1

            # Step 2: chop at junction positions
            boundaries = [0]
            for junc in qualified:
                pos = junc.chop_position
                if 0 < pos < len(seq):
                    boundaries.append(pos)
            boundaries.append(len(seq))

            for i in range(len(boundaries) - 1):
                start = boundaries[i]
                end = boundaries[i + 1]
                seg_len = end - start

                if seg_len < self.min_segment_length:
                    report.segments_discarded_short += 1
                    continue

                report.segments_produced += 1
                seg_seq = seq[start:end]

                # Step 3: orient the sub-read
                oriented_seq = orient_subread(seg_seq)

                if oriented_seq is None:
                    report.oriented_unknown += 1
                    continue

                if oriented_seq == seg_seq:
                    report.oriented_forward += 1
                else:
                    report.oriented_reverse += 1

                # Build output record
                seg_id = f"{record.id}_rescued{i + 1}"
                seg_record = SeqRecord(
                    seq=Seq(oriented_seq),
                    id=seg_id,
                    name=seg_id,
                    description=(
                        f"unknowns_rescued_from={record.id} start={start} end={end}"
                    ),
                )

                # Preserve quality scores (reverse if RC'd)
                if quals:
                    seg_quals = quals[start:end]
                    if oriented_seq != seg_seq:
                        # Quality scores must be reversed for RC reads
                        seg_quals = seg_quals[::-1]
                    seg_record.letter_annotations["phred_quality"] = seg_quals

                output_records.append(seg_record)
                report.output_reads += 1

        # Write output
        write_fastq(output_records, self.output_fastq)

        logger.info(
            f"Unknowns Rescue complete: {report.total_unknowns:,} scanned, "
            f"{report.output_reads:,} oriented reads rescued"
        )
        logger.info(f"\n{report}")
        return report
