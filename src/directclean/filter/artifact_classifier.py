"""
Artifact classifier — main orchestrator for the homopolymer filter.

Reads a BAM file, identifies chimeric reads whose inter-segment
junctions show homopolymer-mediated RT template switching, then
splits the original FASTQ into clean and artifact fractions.

Typical usage::

    classifier = ArtifactClassifier(
        bam_path="aligned.bam",
        input_fastq="restranded.fastq",
        output_dir="results/",
    )
    report = classifier.run()
    print(report)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set

from directclean.filter.junction_parser import iter_chimeric_reads
from directclean.filter.homopolymer import (
    HomopolymerConfig,
    HomopolymerDetector,
    ReadVerdict,
)
from directclean.utils.io import split_fastq_by_ids

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Report data class
# ---------------------------------------------------------------------------

@dataclass
class FilterReport:
    """Summary statistics produced after filtering.

    Attributes:
        total_chimeric_reads:   Reads with SA tag (candidates examined).
        artifact_reads:         Reads flagged as homopolymer artifacts.
        clean_chimeric_reads:   Chimeric reads that passed the filter.
        total_junctions:        Total inter-segment junctions examined.
        artifact_junctions:     Junctions that triggered the filter.
        total_reads_in_fastq:   Reads in the input FASTQ (chimeric + non-chimeric).
        kept_reads:             Reads written to cleaned output.
        removed_reads:          Reads written to removed output.
    """
    total_chimeric_reads: int = 0
    artifact_reads: int = 0
    clean_chimeric_reads: int = 0
    total_junctions: int = 0
    artifact_junctions: int = 0
    total_reads_in_fastq: int = 0
    kept_reads: int = 0
    removed_reads: int = 0
    verdicts: List[ReadVerdict] = field(default_factory=list, repr=False)

    def __str__(self) -> str:
        pct = (
            f"{self.artifact_reads / self.total_chimeric_reads * 100:.1f}%"
            if self.total_chimeric_reads > 0 else "N/A"
        )
        return (
            "=== DirectClean Homopolymer Filter Report ===\n"
            f"  Chimeric reads examined : {self.total_chimeric_reads:,}\n"
            f"  Artifact reads          : {self.artifact_reads:,} ({pct})\n"
            f"  Clean chimeric reads    : {self.clean_chimeric_reads:,}\n"
            f"  Total junctions         : {self.total_junctions:,}\n"
            f"  Artifact junctions      : {self.artifact_junctions:,}\n"
            f"  ---\n"
            f"  Input FASTQ reads       : {self.total_reads_in_fastq:,}\n"
            f"  Kept reads              : {self.kept_reads:,}\n"
            f"  Removed reads           : {self.removed_reads:,}\n"
            "=============================================="
        )


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class ArtifactClassifier:
    """End-to-end homopolymer artifact filter.

    Workflow:
        1. Iterate over chimeric reads in the BAM.
        2. For each read, extract junctions and run the homopolymer
           detector on the flanking sequences.
        3. Collect the set of artifact read IDs.
        4. Stream through the original FASTQ and split into
           ``cleaned.fastq`` and ``removed.fastq``.

    Args:
        bam_path:      Path to minimap2 BAM (sorted & indexed).
        input_fastq:   Original FASTQ that was aligned.
        output_dir:    Directory for output files.
        config:        HomopolymerConfig with detection parameters.
        context_window: Bases to extract on each side of a junction
                        (passed to junction_parser).
        min_mapq:      Minimum mapping quality for a segment to be
                        considered.
        prefix:        Filename prefix for outputs (default "directclean").
    """

    def __init__(
        self,
        bam_path: str | Path,
        input_fastq: str | Path,
        output_dir: str | Path,
        config: HomopolymerConfig | None = None,
        context_window: int = 30,
        min_mapq: int = 0,
        prefix: str = "directclean",
    ) -> None:
        self.bam_path = str(bam_path)
        self.input_fastq = Path(input_fastq)
        self.output_dir = Path(output_dir)
        self.context_window = context_window
        self.min_mapq = min_mapq
        self.prefix = prefix

        self.config = config or HomopolymerConfig()
        self.detector = HomopolymerDetector(self.config)

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ---- output paths ----

    @property
    def cleaned_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.cleaned.fastq"

    @property
    def removed_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.removed.fastq"

    @property
    def report_path(self) -> Path:
        return self.output_dir / f"{self.prefix}.homopolymer_report.tsv"

    # ---- Phase 1: scan BAM ----

    def scan_bam(self) -> tuple[Set[str], FilterReport]:
        """Scan BAM for chimeric reads and classify junctions.

        Returns:
            (artifact_ids, report) — set of artifact read IDs and
            a partially filled FilterReport.
        """
        logger.info("Phase 1: scanning BAM for homopolymer artifacts ...")

        artifact_ids: Set[str] = set()
        report = FilterReport()

        for chimeric in iter_chimeric_reads(
            self.bam_path,
            window_size=self.context_window,
            min_mapq=self.min_mapq,
        ):
            report.total_chimeric_reads += 1
            verdict = self.detector.judge_read(chimeric)
            report.verdicts.append(verdict)
            report.total_junctions += verdict.n_junctions
            report.artifact_junctions += verdict.n_artifact_junctions

            if verdict.is_artifact:
                report.artifact_reads += 1
                artifact_ids.add(verdict.read_id)

        report.clean_chimeric_reads = (
            report.total_chimeric_reads - report.artifact_reads
        )

        logger.info(
            f"BAM scan complete: {report.total_chimeric_reads:,} chimeric reads, "
            f"{report.artifact_reads:,} artifacts"
        )
        return artifact_ids, report

    # ---- Phase 2: split FASTQ ----

    def split_fastq(
        self,
        artifact_ids: Set[str],
        report: FilterReport,
    ) -> FilterReport:
        """Split FASTQ into clean and artifact fractions.

        Args:
            artifact_ids: Read IDs flagged in Phase 1.
            report:       FilterReport to fill with FASTQ counts.

        Returns:
            Updated FilterReport.
        """
        logger.info("Phase 2: splitting FASTQ ...")

        kept, removed = split_fastq_by_ids(
            input_fastq=self.input_fastq,
            remove_ids=artifact_ids,
            output_kept=self.cleaned_fastq,
            output_removed=self.removed_fastq,
        )

        report.kept_reads = kept
        report.removed_reads = removed
        report.total_reads_in_fastq = kept + removed
        return report

    # ---- Phase 3: write TSV report ----

    def write_report(self, report: FilterReport) -> None:
        """Write a per-read TSV report for downstream inspection.

        Columns: read_id, verdict, n_junctions, n_artifact_junctions,
                 junction_details (semicolon-separated).
        """
        logger.info(f"Writing report to {self.report_path}")

        with open(self.report_path, "w") as fh:
            header = (
                "read_id\tverdict\tn_junctions\tn_artifact_junctions"
                "\tjunction_details\n"
            )
            fh.write(header)

            for rv in report.verdicts:
                details_parts = []
                for jv in rv.junction_verdicts:
                    left = jv.junction.left_segment
                    right = jv.junction.right_segment
                    detail = (
                        f"{left.chrom}:{left.ref_start}({left.strand})->"
                        f"{right.chrom}:{right.ref_start}({right.strand})"
                        f"|pos={jv.junction.read_position}"
                        f"|up_hit={jv.upstream_hit.is_hit}"
                        f"(d={jv.upstream_hit.density:.2f},"
                        f"run={jv.upstream_hit.longest_run_length})"
                        f"|dn_hit={jv.downstream_hit.is_hit}"
                        f"(d={jv.downstream_hit.density:.2f},"
                        f"run={jv.downstream_hit.longest_run_length})"
                        f"|artifact={jv.is_artifact}"
                    )
                    details_parts.append(detail)

                line = (
                    f"{rv.read_id}\t"
                    f"{'ARTIFACT' if rv.is_artifact else 'CLEAN'}\t"
                    f"{rv.n_junctions}\t"
                    f"{rv.n_artifact_junctions}\t"
                    f"{';'.join(details_parts)}\n"
                )
                fh.write(line)

    # ---- Public entry point ----

    def run(self) -> FilterReport:
        """Execute the full homopolymer filter pipeline.

        Returns:
            FilterReport with all statistics.
        """
        artifact_ids, report = self.scan_bam()
        report = self.split_fastq(artifact_ids, report)
        self.write_report(report)

        logger.info("Homopolymer filter complete.")
        logger.info(f"\n{report}")
        return report