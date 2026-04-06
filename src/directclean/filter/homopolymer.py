"""
Homopolymer-mediated RT artifact detector.

Wraps the dual-criteria scanning logic (A/T density + longest run)
into a configurable detector class that operates on JunctionInfo and
ChimericRead objects produced by junction_parser.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from directclean.filter.junction_parser import ChimericRead, JunctionInfo
from directclean.utils.sequence_operator import HomopolymerHit, scan_homopolymer

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HomopolymerConfig:
    """Tuneable parameters for the homopolymer filter."""

    scan_window: int = 10
    density_threshold: float = 0.85
    min_run: int = 5
    context_window: int = 50
    require_both_sides: bool = False


# ---------------------------------------------------------------------------
# Per-junction verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JunctionVerdict:
    """Artifact verdict for a single junction."""

    junction: JunctionInfo
    is_artifact: bool
    upstream_hit: HomopolymerHit
    downstream_hit: HomopolymerHit
    hit_source: str = "none"  # upstream | downstream | combined_cross_boundary | combined_upstream | combined_downstream | none


# ---------------------------------------------------------------------------
# Per-read verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReadVerdict:
    """Artifact verdict for an entire read."""

    read_id: str
    is_artifact: bool
    junction_verdicts: list[JunctionVerdict]
    n_junctions: int
    n_artifact_junctions: int


# ---------------------------------------------------------------------------
# Detector class
# ---------------------------------------------------------------------------


class HomopolymerDetector:
    """Configurable detector for homopolymer-mediated RT artifacts."""

    def __init__(self, config: HomopolymerConfig | None = None) -> None:
        self.config = config or HomopolymerConfig()

    def judge_junction(self, junction: JunctionInfo) -> JunctionVerdict:
        """Evaluate one junction for homopolymer artifact signal."""
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

        hit_source = "none"

        if cfg.require_both_sides:
            is_artifact = upstream_hit.is_hit and downstream_hit.is_hit
            if is_artifact:
                hit_source = "both_sides"
        else:
            is_artifact = upstream_hit.is_hit or downstream_hit.is_hit
            if upstream_hit.is_hit and downstream_hit.is_hit:
                hit_source = "both_sides"
            elif upstream_hit.is_hit:
                hit_source = "upstream"
            elif downstream_hit.is_hit:
                hit_source = "downstream"

        # Combined scan: rescue motifs that straddle the junction boundary.
        if not is_artifact:
            combined_seq = junction.upstream_seq + junction.downstream_seq
            combined_hit = scan_homopolymer(
                combined_seq,
                window_size=cfg.scan_window,
                density_threshold=cfg.density_threshold,
                min_run=cfg.min_run,
            )

            if combined_hit.is_hit:
                is_artifact = True
                boundary = len(junction.upstream_seq)

                start = combined_hit.window_start
                end = combined_hit.window_end
                center = (start + end) / 2.0

                crosses_boundary = (
                    (start < boundary < end) or (start == boundary) or (end == boundary)
                )

                if crosses_boundary:
                    upstream_hit = combined_hit
                    hit_source = "combined_cross_boundary"
                elif center < boundary:
                    upstream_hit = combined_hit
                    hit_source = "combined_upstream"
                else:
                    downstream_hit = combined_hit
                    hit_source = "combined_downstream"

        return JunctionVerdict(
            junction=junction,
            is_artifact=is_artifact,
            upstream_hit=upstream_hit,
            downstream_hit=downstream_hit,
            hit_source=hit_source,
        )

    def judge_read(self, chimeric_read: ChimericRead) -> ReadVerdict:
        """Evaluate all junctions of a chimeric read."""
        verdicts: list[JunctionVerdict] = []
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
