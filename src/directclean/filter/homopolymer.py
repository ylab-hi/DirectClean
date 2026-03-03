"""
Homopolymer-mediated RT artifact detector.

Wraps the dual-criteria scanning logic (A/T density + longest run)
into a configurable detector class that operates on JunctionInfo and
ChimericRead objects produced by junction_parser.

Biological background
---------------------
During reverse transcription of Direct-cDNA libraries, the RT enzyme
can dissociate at A/T-rich (especially poly-A tail) regions and
re-prime on a different mRNA molecule.  This "template switching"
produces chimeric reads that look like gene fusions but are artifacts.

The detector flags junctions whose flanking sequences are dominated
by A/T homopolymers — the hallmark of this mechanism.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List

from directclean.filter.junction_parser import ChimericRead, JunctionInfo
from directclean.utils.sequence_operator import scan_homopolymer, HomopolymerHit

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HomopolymerConfig:
    """Tuneable parameters for the homopolymer filter.

    Attributes:
        scan_window:        Sliding window size used by scan_homopolymer
                            on each flanking sequence (default 10 bp).
        density_threshold:  Minimum A/T fraction within a scanning window
                            to count as a hit (default 0.8).
        min_run:            Minimum consecutive A or T to count as a hit
                            (default 5).
        context_window:     Already applied in junction_parser when
                            extracting upstream/downstream — stored here
                            for bookkeeping (default 30 bp).
        require_both_sides: If True, both upstream AND downstream must
                            hit to call an artifact.  If False, either
                            side hitting is sufficient (default True).
    """
    scan_window: int = 10
    density_threshold: float = 0.8
    min_run: int = 5
    context_window: int = 30
    require_both_sides: bool = True


# ---------------------------------------------------------------------------
# Per-junction verdict
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class JunctionVerdict:
    """Artifact verdict for a single junction.

    Attributes:
        junction:       The JunctionInfo that was evaluated.
        is_artifact:    True if this junction is flagged as an RT artifact.
        upstream_hit:   HomopolymerHit result for the upstream flank.
        downstream_hit: HomopolymerHit result for the downstream flank.
    """
    junction: JunctionInfo
    is_artifact: bool
    upstream_hit: HomopolymerHit
    downstream_hit: HomopolymerHit


# ---------------------------------------------------------------------------
# Per-read verdict
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReadVerdict:
    """Artifact verdict for an entire read.

    A read is called artifact if **any** of its junctions is flagged.

    Attributes:
        read_id:            Read name.
        is_artifact:        True → goes to removed.fastq.
        junction_verdicts:  Per-junction details.
        n_junctions:        Total number of inter-segment junctions.
        n_artifact_junctions: How many were flagged.
    """
    read_id: str
    is_artifact: bool
    junction_verdicts: List[JunctionVerdict]
    n_junctions: int
    n_artifact_junctions: int


# ---------------------------------------------------------------------------
# Detector class
# ---------------------------------------------------------------------------

class HomopolymerDetector:
    """Configurable detector for homopolymer-mediated RT artifacts.

    Usage::

        cfg = HomopolymerConfig(scan_window=10, density_threshold=0.8, min_run=3)
        detector = HomopolymerDetector(cfg)

        # From a ChimericRead produced by junction_parser
        verdict = detector.judge_read(chimeric_read)
        if verdict.is_artifact:
            print(f"{verdict.read_id} is an artifact")
    """

    def __init__(self, config: HomopolymerConfig | None = None) -> None:
        self.config = config or HomopolymerConfig()

    # ----- single junction -----

    def judge_junction(self, junction: JunctionInfo) -> JunctionVerdict:
        """Evaluate one junction for homopolymer artifact signal.

        Scans the upstream and downstream flanking sequences extracted
        by junction_parser, using the dual-criteria sliding window.

        Args:
            junction: JunctionInfo with upstream_seq / downstream_seq.

        Returns:
            JunctionVerdict with detailed hit information.
        """
        cfg = self.config

        upstream_hit = scan_homopolymer(
            junction.upstream_seq,
            window_size=cfg.scan_window,
            density_threshold=cfg.density_threshold,
            min_run=cfg.min_run,
        )
        downstream_hit = scan_homopolymer(
            junction.downstream_seq,
            window_size=cfg.scan_window,
            density_threshold=cfg.density_threshold,
            min_run=cfg.min_run,
        )

        if cfg.require_both_sides:
            is_artifact = upstream_hit.is_hit and downstream_hit.is_hit
        else:
            is_artifact = upstream_hit.is_hit or downstream_hit.is_hit

        return JunctionVerdict(
            junction=junction,
            is_artifact=is_artifact,
            upstream_hit=upstream_hit,
            downstream_hit=downstream_hit,
        )

    # ----- whole read -----

    def judge_read(self, chimeric_read: ChimericRead) -> ReadVerdict:
        """Evaluate all junctions of a chimeric read.

        A read is flagged as artifact if **any** junction is positive.

        Args:
            chimeric_read: ChimericRead from junction_parser.

        Returns:
            ReadVerdict summarising the outcome.
        """
        verdicts: List[JunctionVerdict] = []
        n_artifact = 0

        for junc in chimeric_read.junctions:
            v = self.judge_junction(junc)
            verdicts.append(v)
            if v.is_artifact:
                n_artifact += 1

        return ReadVerdict(
            read_id=chimeric_read.read_id,
            is_artifact=n_artifact > 0,
            junction_verdicts=verdicts,
            n_junctions=len(verdicts),
            n_artifact_junctions=n_artifact,
        )