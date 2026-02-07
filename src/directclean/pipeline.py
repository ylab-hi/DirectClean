"""
DirectClean pipeline — end-to-end orchestration.

Chains the three processing stages in order::

    Input FASTQ (from Breakinator + Restrander)
        │
        ▼
    Stage 1: Rescuer
        Detect internal TSO/RTP adapters, chop reads
        │
        ▼
    Stage 2: Minimap2
        Align rescued reads to reference genome
        │
        ▼
    Stage 3: Homopolymer Filter
        Identify and remove RT template switching artifacts
        │
        ▼
    Final outputs:
        cleaned.fastq   — reads passing all filters
        removed.fastq    — artifact reads
        reports/         — per-stage TSV reports
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from directclean.rescuer import ReadChopper, AdapterConfig, RescueReport
from directclean.external.minimap2 import Minimap2Aligner
from directclean.filter import (
    ArtifactClassifier,
    HomopolymerConfig,
    FilterReport,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline configuration
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    """All tuneable parameters for the full pipeline.

    Attributes:
        adapter_config:     Parameters for internal adapter detection.
        homopolymer_config: Parameters for homopolymer filtering.
        min_confidence:     Minimum confidence to chop (1-3, default 2).
        context_window:     Bases flanking each junction for homo check.
        min_mapq:           Minimum mapping quality for filter stage.
        threads:            Threads for minimap2 and samtools.
        skip_rescue:        Skip the Rescuer stage entirely.
        skip_filter:        Skip the Homopolymer Filter stage entirely.
    """
    adapter_config: AdapterConfig = None
    homopolymer_config: HomopolymerConfig = None
    min_confidence: int = 2
    context_window: int = 30
    min_mapq: int = 0
    threads: int = 4
    skip_rescue: bool = False
    skip_filter: bool = False

    def __post_init__(self):
        if self.adapter_config is None:
            self.adapter_config = AdapterConfig()
        if self.homopolymer_config is None:
            self.homopolymer_config = HomopolymerConfig()


# ---------------------------------------------------------------------------
# Pipeline report
# ---------------------------------------------------------------------------

@dataclass
class PipelineReport:
    """Combined report from all pipeline stages.

    Attributes:
        rescue_report: Report from the Rescuer stage (or None if skipped).
        filter_report: Report from the Homopolymer Filter (or None).
        elapsed_seconds: Total wall-clock time.
    """
    rescue_report: Optional[RescueReport] = None
    filter_report: Optional[FilterReport] = None
    elapsed_seconds: float = 0.0

    def __str__(self) -> str:
        parts = [
            "╔══════════════════════════════════════════════╗",
            "║       DirectClean Pipeline Report            ║",
            "╚══════════════════════════════════════════════╝",
        ]

        if self.rescue_report is not None:
            parts.append("")
            parts.append(str(self.rescue_report))

        if self.filter_report is not None:
            parts.append("")
            parts.append(str(self.filter_report))

        parts.append("")
        parts.append(f"  Total elapsed time: {self.elapsed_seconds:.1f}s")
        return "\n".join(parts)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class DirectCleanPipeline:
    """End-to-end DirectClean pipeline.

    Usage::

        pipeline = DirectCleanPipeline(
            input_fastq="restranded.fastq",
            reference="genome.fa",
            output_dir="results/",
            threads=8,
        )
        report = pipeline.run()
        print(report)

    Args:
        input_fastq: Input FASTQ from Breakinator + Restrander.
        reference:   Reference genome FASTA for minimap2.
        output_dir:  Output directory for all results.
        config:      PipelineConfig with all parameters.
        threads:     Shortcut for config.threads (overrides if set).
        prefix:      Filename prefix for outputs.
    """

    def __init__(
        self,
        input_fastq: str | Path,
        reference: str | Path,
        output_dir: str | Path,
        config: PipelineConfig | None = None,
        threads: int | None = None,
        prefix: str = "directclean",
    ) -> None:
        self.input_fastq = Path(input_fastq)
        self.reference = Path(reference)
        self.output_dir = Path(output_dir)
        self.prefix = prefix

        self.config = config or PipelineConfig()
        if threads is not None:
            self.config.threads = threads

        # Validate inputs
        if not self.input_fastq.exists():
            raise FileNotFoundError(f"Input FASTQ not found: {self.input_fastq}")
        if not self.reference.exists():
            raise FileNotFoundError(f"Reference not found: {self.reference}")

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._reports_dir.mkdir(parents=True, exist_ok=True)

    # ---- output paths ----

    @property
    def _reports_dir(self) -> Path:
        return self.output_dir / "reports"

    @property
    def _rescued_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.rescued.fastq"

    @property
    def _aligned_bam(self) -> Path:
        return self.output_dir / f"{self.prefix}.aligned.sorted.bam"

    @property
    def cleaned_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.cleaned.fastq"

    @property
    def removed_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.removed.fastq"

    # ---- Stage 1: Rescue ----

    def _run_rescue(self) -> tuple[Path, Optional[RescueReport]]:
        """Stage 1: detect and chop internal TSO/RTP adapters.

        Returns:
            (output_fastq, report) — the FASTQ to feed into stage 2.
        """
        if self.config.skip_rescue:
            logger.info("Stage 1 (Rescue): SKIPPED")
            return self.input_fastq, None

        logger.info("=" * 50)
        logger.info("Stage 1: Rescue — internal adapter detection")
        logger.info("=" * 50)

        chopper = ReadChopper(
            input_fastq=self.input_fastq,
            output_fastq=self._rescued_fastq,
            config=self.config.adapter_config,
            min_confidence=self.config.min_confidence,
            report_path=self._reports_dir / f"{self.prefix}.rescue_report.tsv",
        )
        report = chopper.run()
        return self._rescued_fastq, report

    # ---- Stage 2: Minimap2 ----

    def _run_alignment(self, fastq: Path) -> Path:
        """Stage 2: align reads to reference genome.

        Args:
            fastq: Input FASTQ (from rescue stage or original).

        Returns:
            Path to sorted, indexed BAM.
        """
        logger.info("=" * 50)
        logger.info("Stage 2: Minimap2 — splice-aware alignment")
        logger.info("=" * 50)

        aligner = Minimap2Aligner(
            reference=self.reference,
            threads=self.config.threads,
            sample_id=self.prefix,
        )
        bam = aligner.align(
            fastq=fastq,
            output_bam=self._aligned_bam,
        )
        return bam

    # ---- Stage 3: Homopolymer Filter ----

    def _run_filter(
        self,
        bam: Path,
        fastq: Path,
    ) -> Optional[FilterReport]:
        """Stage 3: detect and remove homopolymer-mediated artifacts.

        Args:
            bam:   Aligned BAM from stage 2.
            fastq: The FASTQ that was aligned (for splitting).

        Returns:
            FilterReport, or None if skipped.
        """
        if self.config.skip_filter:
            logger.info("Stage 3 (Homopolymer Filter): SKIPPED")
            # If filter is skipped, just copy the input FASTQ as "cleaned"
            import shutil
            shutil.copy2(fastq, self.cleaned_fastq)
            return None

        logger.info("=" * 50)
        logger.info("Stage 3: Homopolymer Filter — RT artifact removal")
        logger.info("=" * 50)

        classifier = ArtifactClassifier(
            bam_path=bam,
            input_fastq=fastq,
            output_dir=self.output_dir,
            config=self.config.homopolymer_config,
            context_window=self.config.context_window,
            min_mapq=self.config.min_mapq,
            prefix=self.prefix,
        )
        report = classifier.run()
        return report

    # ---- Public entry point ----

    def run(self) -> PipelineReport:
        """Execute the full DirectClean pipeline.

        Stages:
            1. Rescue: chop reads with internal TSO/RTP adapters.
            2. Minimap2: align rescued reads to reference.
            3. Homopolymer Filter: remove RT template switching artifacts.

        Returns:
            PipelineReport with combined statistics.
        """
        start_time = time.time()
        pipeline_report = PipelineReport()

        logger.info(f"DirectClean v{self._version} starting")
        logger.info(f"  Input:     {self.input_fastq}")
        logger.info(f"  Reference: {self.reference}")
        logger.info(f"  Output:    {self.output_dir}")
        logger.info(f"  Threads:   {self.config.threads}")

        # Stage 1
        stage1_fastq, rescue_report = self._run_rescue()
        pipeline_report.rescue_report = rescue_report

        # Stage 2
        bam = self._run_alignment(stage1_fastq)

        # Stage 3
        filter_report = self._run_filter(bam, stage1_fastq)
        pipeline_report.filter_report = filter_report

        pipeline_report.elapsed_seconds = time.time() - start_time

        logger.info("")
        logger.info(str(pipeline_report))

        return pipeline_report

    @property
    def _version(self) -> str:
        try:
            from directclean import __version__
            return __version__
        except ImportError:
            return "unknown"