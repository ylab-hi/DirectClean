"""
Internal adapter finder using anchor-and-extend strategy.

Detects internal ``polyA → RTP_rc → TSO`` signatures in reads that
indicate two cDNA molecules were ligated together.  Uses a tiered
confidence system:

* **High confidence** (3/3): polyA + RTP_rc + TSO all found.
* **Medium confidence** (2/3): any two of the three signals found.
* **Low confidence** (1/3): only TSO found alone in the middle of the
  read, with very low edit distance — optional, controlled by
  ``allow_single_tso``.

Algorithm overview
------------------
1. Scan the read for poly-A anchors (long A-runs in the interior).
2. For each poly-A anchor, search downstream for RTP_rc (fuzzy).
3. For each RTP_rc hit, search further downstream for TSO (fuzzy).
4. Also do a standalone TSO scan to catch cases where polyA/RTP_rc
   are degraded but TSO is clear.
5. Merge and deduplicate hits, assign confidence tiers, and return
   candidate chop sites.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import edlib

from directclean.rescuer.adaptor_seq import AdapterConfig
from directclean.utils.sequence_operator import at_density

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AdapterHit:
    """A single fuzzy match of an adapter sequence in a read.

    Attributes:
        start:         0-based start position on the read.
        end:           0-based exclusive end position on the read.
        edit_distance: Number of edits (substitutions + indels).
        label:         Which adapter: 'TSO', 'RTP_rc', or 'polyA'.
    """
    start: int
    end: int
    edit_distance: int
    label: str


@dataclass(frozen=True)
class InternalJunction:
    """A detected internal adapter junction — a candidate chop site.

    Attributes:
        chop_position: Where to cut the read (0-based).  This is the
                       start of the TSO if found, otherwise the end
                       of RTP_rc, otherwise the end of polyA.
        confidence:    Number of signals found (1, 2, or 3).
        polya_hit:     The poly-A anchor hit, or None.
        rtp_rc_hit:    The RTP_rc hit, or None.
        tso_hit:       The TSO hit, or None.
    """
    chop_position: int
    confidence: int
    polya_hit: Optional[AdapterHit]
    rtp_rc_hit: Optional[AdapterHit]
    tso_hit: Optional[AdapterHit]


@dataclass
class FinderResult:
    """All internal junctions found in a single read.

    Attributes:
        read_id:    Read name.
        read_len:   Length of the read.
        junctions:  List of InternalJunction, sorted by position.
        n_chops:    Number of chop sites (= number of junctions).
    """
    read_id: str
    read_len: int
    junctions: List[InternalJunction] = field(default_factory=list)

    @property
    def n_chops(self) -> int:
        return len(self.junctions)

    @property
    def has_internal_adapter(self) -> bool:
        return self.n_chops > 0


# ---------------------------------------------------------------------------
# Core search functions
# ---------------------------------------------------------------------------

def _find_polya_anchors(
    sequence: str,
    min_run: int,
    density_window: int,
    density_threshold: float,
    five_prime_tol: int,
    three_prime_tol: int,
) -> List[AdapterHit]:
    """Find interior poly-A regions in a read.

    Scans for runs of consecutive A's, then validates each run with
    a density check in a surrounding window.  Hits near the 5' or
    3' ends are excluded (those are normal poly-A tails).

    Returns:
        List of AdapterHit with label='polyA'.
    """
    seq = sequence.upper()
    n = len(seq)
    hits: List[AdapterHit] = []

    i = 0
    while i < n:
        if seq[i] == "A":
            start = i
            while i < n and seq[i] == "A":
                i += 1
            run_len = i - start

            if run_len < min_run:
                continue

            # Check density in a window around the run
            win_start = max(0, start - density_window // 2)
            win_end = min(n, i + density_window // 2)
            window_seq = seq[win_start:win_end]
            if at_density(window_seq) < density_threshold:
                continue

            # Exclude hits at read ends (those are normal poly-A tails)
            if start < five_prime_tol:
                continue
            if i > n - three_prime_tol:
                continue

            hits.append(AdapterHit(
                start=start, end=i, edit_distance=0, label="polyA"
            ))
        else:
            i += 1

    return hits


def _fuzzy_search(
    sequence: str,
    query: str,
    region_start: int,
    region_end: int,
    max_edit_distance: int,
) -> Optional[AdapterHit]:
    """Fuzzy search for a short query in a region of the read.

    Uses edlib in 'HW' (semi-global) mode: the query is fully
    aligned but the read region can have unaligned flanks.

    Args:
        sequence:          Full read sequence.
        query:             Adapter sequence to find.
        region_start:      Start of the search region (0-based).
        region_end:        End of the search region (exclusive).
        max_edit_distance: Maximum allowed edits.

    Returns:
        AdapterHit if found within edit distance, else None.
    """
    seq = sequence.upper()
    qry = query.upper()

    # Clamp to valid range
    region_start = max(0, region_start)
    region_end = min(len(seq), region_end)

    if region_end - region_start < len(qry) // 2:
        return None

    region = seq[region_start:region_end]

    result = edlib.align(
        qry, region,
        mode="HW",                    # semi-global
        task="locations",             # we need start/end positions
        k=max_edit_distance,          # max edit distance cutoff
    )

    if result["editDistance"] == -1:
        # No match within edit distance
        return None

    # edlib returns list of (start, end) tuples; take the best (first)
    locations = result["locations"]
    if not locations:
        return None

    loc_start, loc_end = locations[0]

    return AdapterHit(
        start=region_start + loc_start,
        end=region_start + loc_end + 1,  # edlib end is inclusive
        edit_distance=result["editDistance"],
        label="",  # caller sets this
    )


def _search_all_tso(
    sequence: str,
    tso_seq: str,
    max_edit_distance: int,
    five_prime_tol: int,
    three_prime_tol: int,
) -> List[AdapterHit]:
    """Scan the entire read for all TSO occurrences.

    Used for the standalone TSO search path (catches cases where
    polyA/RTP_rc are degraded).  Excludes hits near read ends.

    Returns:
        List of AdapterHit with label='TSO'.
    """
    seq = sequence.upper()
    qry = tso_seq.upper()
    n = len(seq)
    qlen = len(qry)
    hits: List[AdapterHit] = []

    # Slide a window across the interior of the read
    search_start = five_prime_tol
    search_end = n - three_prime_tol

    if search_end - search_start < qlen:
        return hits

    # Use edlib on the full interior region
    region = seq[search_start:search_end]
    result = edlib.align(
        qry, region,
        mode="HW",
        task="locations",
        k=max_edit_distance,
    )

    if result["editDistance"] == -1 or not result.get("locations"):
        return hits

    for loc_start, loc_end in result["locations"]:
        abs_start = search_start + loc_start
        abs_end = search_start + loc_end + 1

        hits.append(AdapterHit(
            start=abs_start,
            end=abs_end,
            edit_distance=result["editDistance"],
            label="TSO",
        ))

    return hits


# ---------------------------------------------------------------------------
# Main finder class
# ---------------------------------------------------------------------------

class AdapterFinder:
    """Detect internal adapter junctions in Direct-cDNA reads.

    Implements the anchor-and-extend strategy with tiered confidence:

    1. Find poly-A anchors in the interior of the read.
    2. Extend downstream to find RTP_rc.
    3. Extend further downstream to find TSO.
    4. Also do a standalone TSO scan for degraded cases.
    5. Merge results and assign confidence levels.

    Usage::

        finder = AdapterFinder(AdapterConfig())
        result = finder.find(read_id="read1", sequence="ATCG...")
        for junc in result.junctions:
            print(f"Chop at {junc.chop_position}, confidence={junc.confidence}")
    """

    def __init__(self, config: AdapterConfig | None = None) -> None:
        self.config = config or AdapterConfig()

    def find(self, read_id: str, sequence: str) -> FinderResult:
        """Find all internal adapter junctions in a read.

        Args:
            read_id:  Read name.
            sequence: Full read sequence (5'→3', already restranded).

        Returns:
            FinderResult with detected junctions.
        """
        cfg = self.config
        seq = sequence.upper()
        n = len(seq)

        result = FinderResult(read_id=read_id, read_len=n)

        if n < cfg.five_prime_tolerance + cfg.three_prime_tolerance:
            return result  # read too short to have internal adapters

        # ---- Step 1: anchor-and-extend from polyA ----
        polya_hits = _find_polya_anchors(
            seq,
            min_run=cfg.polya_min_run,
            density_window=cfg.polya_density_window,
            density_threshold=cfg.polya_density_threshold,
            five_prime_tol=cfg.five_prime_tolerance,
            three_prime_tol=cfg.three_prime_tolerance,
        )

        used_tso_positions: set[int] = set()

        for pa in polya_hits:
            # ---- Step 2: search for RTP_rc downstream of polyA ----
            rtp_rc_hit = _fuzzy_search(
                seq,
                query=cfg.rtp_rc_seq,
                region_start=pa.end,
                region_end=pa.end + cfg.rtp_rc_search_range,
                max_edit_distance=cfg.max_edit_distance,
            )
            if rtp_rc_hit is not None:
                rtp_rc_hit = AdapterHit(
                    rtp_rc_hit.start, rtp_rc_hit.end,
                    rtp_rc_hit.edit_distance, "RTP_rc"
                )

            # ---- Step 3: search for TSO downstream of RTP_rc (or polyA) ----
            tso_search_start = (
                rtp_rc_hit.end if rtp_rc_hit is not None else pa.end
            )
            tso_hit = _fuzzy_search(
                seq,
                query=cfg.tso_seq,
                region_start=tso_search_start,
                region_end=tso_search_start + cfg.tso_search_range,
                max_edit_distance=cfg.max_edit_distance,
            )
            if tso_hit is not None:
                tso_hit = AdapterHit(
                    tso_hit.start, tso_hit.end,
                    tso_hit.edit_distance, "TSO"
                )
                used_tso_positions.add(tso_hit.start)

            # ---- Determine confidence and chop position ----
            signals = [pa, rtp_rc_hit, tso_hit]
            confidence = sum(1 for s in signals if s is not None)

            if confidence < 2:
                # polyA alone is not enough — skip
                continue

            # Chop position priority: TSO start > RTP_rc end > polyA end
            if tso_hit is not None:
                chop_pos = tso_hit.start
            elif rtp_rc_hit is not None:
                chop_pos = rtp_rc_hit.end
            else:
                chop_pos = pa.end

            result.junctions.append(InternalJunction(
                chop_position=chop_pos,
                confidence=confidence,
                polya_hit=pa,
                rtp_rc_hit=rtp_rc_hit,
                tso_hit=tso_hit,
            ))

        # ---- Step 4: standalone TSO scan for degraded cases ----
        standalone_tso_hits = _search_all_tso(
            seq,
            tso_seq=cfg.tso_seq,
            max_edit_distance=min(1, cfg.max_edit_distance),  # stricter
            five_prime_tol=cfg.five_prime_tolerance,
            three_prime_tol=cfg.three_prime_tolerance,
        )

        for tso in standalone_tso_hits:
            # Skip if already captured by anchor-and-extend
            if tso.start in used_tso_positions:
                continue

            # Check if any existing junction is close (within 150bp)
            too_close = any(
                abs(tso.start - j.chop_position) < 150
                for j in result.junctions
            )
            if too_close:
                continue

            # Standalone TSO: confidence = 1
            result.junctions.append(InternalJunction(
                chop_position=tso.start,
                confidence=1,
                polya_hit=None,
                rtp_rc_hit=None,
                tso_hit=tso,
            ))

        # Deduplicate junctions at the same (or very close) chop position.
        # Keep the one with highest confidence.
        result.junctions.sort(key=lambda j: (-j.confidence, j.chop_position))
        deduped: list[InternalJunction] = []
        used_positions: set[int] = set()
        for j in result.junctions:
            # Consider positions within 20bp as the same junction
            if any(abs(j.chop_position - p) < 20 for p in used_positions):
                continue
            deduped.append(j)
            used_positions.add(j.chop_position)
        result.junctions = sorted(deduped, key=lambda j: j.chop_position)

        if result.has_internal_adapter:
            logger.debug(
                f"Read {read_id}: found {result.n_chops} internal junction(s)"
            )

        return result