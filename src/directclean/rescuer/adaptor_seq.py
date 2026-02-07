"""
Adapter sequence definitions for Oxford Nanopore Direct-cDNA kits.

Contains TSO, RTP, and derived sequences used by the Rescuer module
to detect internal adapter junctions in reads that have already been
processed by Breakinator + Restrander.

In a correctly oriented 5'→3' read the expected layout is::

    TSO --- gene body --- polyA --- RTP_rc

When two cDNA molecules are ligated together the read looks like::

    TSO --- geneA --- polyA --- RTP_rc --- [junk] --- TSO --- geneB --- polyA --- RTP_rc

The Rescuer detects the internal ``polyA → RTP_rc → TSO`` signature
and chops the read at the TSO boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from directclean.utils.sequence_operator import reverse_complement


# ---------------------------------------------------------------------------
# Raw primer sequences (as ordered / from kit spec)
# ---------------------------------------------------------------------------

# Strand-Switching Primer (SSP / TSO)
# Original: 5'-TTTCTGTTGGTGCTGATATTGCTmGmGmG-3'
# We drop the 2'-O-methyl G modifications for sequence matching.
TSO_SEQUENCE = "TTTCTGTTGGTGCTGATATTGCT"

# VN Primer (RTP)
# Original: 5'-/5phos/ACTTGCCTGTCGCTCTATCTTCTTTTTTTTTTTTTTTTTTTTVN-3'
# We use the core sequence without the 5' phosphate and 3' VN anchor.
RTP_SEQUENCE = "ACTTGCCTGTCGCTCTATCTTC"

# Reverse complement of RTP — this is what appears in a 5'→3' read
# downstream of the poly-A tail.
RTP_RC_SEQUENCE = reverse_complement(RTP_SEQUENCE)


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AdapterConfig:
    """Tuneable parameters for internal adapter detection.

    Attributes:
        tso_seq:            TSO sequence to search for.
        rtp_rc_seq:         Reverse-complement of RTP to search for.

        max_edit_distance:  Maximum edit distance for fuzzy adapter
                            matching (default 3).

        polya_min_run:      Minimum consecutive A's to call a poly-A
                            anchor (default 10).
        polya_density_window: Window size for A-density check around a
                            poly-A candidate (default 15).
        polya_density_threshold: Minimum A/(A+other) density within
                            the window (default 0.8).

        rtp_rc_search_range: How far downstream of a poly-A anchor to
                            look for RTP_rc (default 80 bp).
        tso_search_range:   How far downstream of an RTP_rc hit to
                            look for TSO (default 80 bp).

        five_prime_tolerance:  Hits within this distance from the 5'
                            end are considered normal 5' adapters and
                            are ignored (default 100 bp).
        three_prime_tolerance: Hits within this distance from the 3'
                            end are considered normal 3' adapters and
                            are ignored (default 50 bp).

        min_segment_length: Minimum length of a rescued sub-read to
                            keep (default 50 bp).  Very short fragments
                            are likely junk.
    """
    tso_seq: str = TSO_SEQUENCE
    rtp_rc_seq: str = RTP_RC_SEQUENCE

    max_edit_distance: int = 3

    polya_min_run: int = 10
    polya_density_window: int = 15
    polya_density_threshold: float = 0.8

    rtp_rc_search_range: int = 80
    tso_search_range: int = 80

    five_prime_tolerance: int = 100
    three_prime_tolerance: int = 50

    min_segment_length: int = 50