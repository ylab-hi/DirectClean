"""
DirectClean pipeline — end-to-end orchestration.

Chains all five processing stages::

    Raw FASTQ (user input)
        │
        ▼
    Stage 1: Breakinator
        minimap2 → SAM → Breakinator → remove foldback reads
        │
        ▼
    Stage 2: Restrander
        Correct strand orientation → trim primers → 5'→3' reads
        │
        ▼
    Stage 3: Rescuer
        Detect internal TSO/RTP adapters → chop chimeric reads
        │
        ▼
    Stage 4: Minimap2
        Splice-aware alignment to reference genome
        │
        ▼
    Stage 5: Homopolymer Filter
        Identify and remove RT template switching artifacts
        │
        ▼
    Final outputs:
        cleaned.fastq   — reads passing all stages
        removed.fastq    — artifact reads from homopolymer filter
        reports/         — per-stage reports
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from directclean.external.breakinator import BreakinatorRunner, BreakReport
from directclean.external.restrander import RestranderRunner, RestranderReport
from directclean.external.minimap2 import Minimap2Aligner
from directclean.external.dependencies import check_all_dependencies
from directclean.rescuer import ReadChopper, AdapterConfig, RescueReport
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
        threads:            Threads for minimap2, samtools, breakinator.
        junc_bed:           Optional junction BED for Breakinator alignment.
    """
    adapter_config: AdapterConfig = None
    homopolymer_config: HomopolymerConfig = None
    min_confidence: int = 2
    context_window: int = 30
    min_mapq: int = 0
    threads: int = 4
    junc_bed: Optional[Path] = None

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
        break_report:      Report from Breakinator (Stage 1).
        restrander_report: Report from Restrander (Stage 2).
        rescue_report:     Report from Rescuer (Stage 3).
        filter_report:     Report from Homopolymer Filter (Stage 5).
        elapsed_seconds:   Total wall-clock time.
    """
    break_report: Optional[BreakReport] = None
    restrander_report: Optional[RestranderReport] = None
    rescue_report: Optional[RescueReport] = None
    filter_report: Optional[FilterReport] = None
    elapsed_seconds: float = 0.0

    def __str__(self) -> str:
        parts = [
            "",
            "╔══════════════════════════════════════════════╗",
            "║       DirectClean Pipeline Report            ║",
            "╚══════════════════════════════════════════════╝",
        ]

        if self.break_report is not None:
            parts.append("")
            parts.append(str(self.break_report))

        if self.restrander_report is not None:
            parts.append("")
            parts.append(str(self.restrander_report))

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

    Processes raw Direct-cDNA FASTQ through all five stages:
    Breakinator → Restrander → Rescuer → Minimap2 → Homopolymer Filter.

    Usage::

        pipeline = DirectCleanPipeline(
            input_fastq="raw.fastq",
            reference="genome.fa",
            output_dir="results/",
        )
        report = pipeline.run()
        print(report)

    Args:
        input_fastq: Raw input FASTQ from sequencing.
        reference:   Reference genome FASTA.
        output_dir:  Output directory for all results.
        config:      PipelineConfig with all parameters.
        prefix:      Filename prefix for outputs.
    """

    def __init__(
        self,
        input_fastq: str | Path,
        reference: str | Path,
        output_dir: str | Path,
        config: PipelineConfig | None = None,
        prefix: str = "directclean",
    ) -> None:
        self.input_fastq = Path(input_fastq)
        self.reference = Path(reference)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.config = config or PipelineConfig()

        # Validate inputs
        if not self.input_fastq.exists():
            raise FileNotFoundError(
                f"Input FASTQ not found: {self.input_fastq}"
            )
        if not self.reference.exists():
            raise FileNotFoundError(
                f"Reference not found: {self.reference}"
            )

        # Validate external dependencies up front
        check_all_dependencies()

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._reports_dir.mkdir(parents=True, exist_ok=True)
        self._intermediates_dir.mkdir(parents=True, exist_ok=True)

    # ---- output paths ----

    @property
    def _reports_dir(self) -> Path:
        return self.output_dir / "reports"

    @property
    def _intermediates_dir(self) -> Path:
        return self.output_dir / "intermediates"

    @property
    def _no_foldback_fastq(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.no_foldback.fastq"

    @property
    def _restranded_fastq(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.restranded.fastq"

    @property
    def _rescued_fastq(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.rescued.fastq"

    @property
    def _aligned_bam(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.aligned.sorted.bam"

    @property
    def cleaned_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.cleaned.fastq"

    @property
    def removed_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.removed.fastq"

    # ---- Stage 1: Breakinator ----

    def _run_breakinator(self) -> tuple[Path, BreakReport]:
        """Stage 1: remove foldback inversion artifacts.

        Returns:
            (output_fastq, report)
        """
        logger.info("=" * 55)
        logger.info("Stage 1/5: Breakinator — foldback artifact removal")
        logger.info("=" * 55)

        runner = BreakinatorRunner(
            reference=self.reference,
            threads=self.config.threads,
            junc_bed=self.config.junc_bed,
        )
        report = runner.run(
            input_fastq=self.input_fastq,
            output_fastq=self._no_foldback_fastq,
            work_dir=self._intermediates_dir,
            prefix=self.prefix,
        )
        return self._no_foldback_fastq, report

    # ---- Stage 2: Restrander ----

    def _run_restrander(self, fastq: Path) -> tuple[Path, RestranderReport]:
        """Stage 2: correct strand orientation and trim primers.

        Args:
            fastq: Foldback-free FASTQ from Stage 1.

        Returns:
            (output_fastq, report)
        """
        logger.info("=" * 55)
        logger.info("Stage 2/5: Restrander — strand orientation correction")
        logger.info("=" * 55)

        runner = RestranderRunner()
        report = runner.run(
            input_fastq=fastq,
            output_fastq=self._restranded_fastq,
        )
        return self._restranded_fastq, report

    # ---- Stage 3: Rescue ----

    def _run_rescue(self, fastq: Path) -> tuple[Path, RescueReport]:
        """Stage 3: detect and chop internal TSO/RTP adapters.

        Args:
            fastq: Restranded FASTQ from Stage 2.

        Returns:
            (output_fastq, report)
        """
        logger.info("=" * 55)
        logger.info("Stage 3/5: Rescuer — internal adapter detection")
        logger.info("=" * 55)

        chopper = ReadChopper(
            input_fastq=fastq,
            output_fastq=self._rescued_fastq,
            config=self.config.adapter_config,
            min_confidence=self.config.min_confidence,
            report_path=self._reports_dir / f"{self.prefix}.rescue_report.tsv",
            threads=self.config.threads,
        )
        report = chopper.run()
        return self._rescued_fastq, report

    # ---- Stage 4: Minimap2 ----

    def _run_alignment(self, fastq: Path) -> Path:
        """Stage 4: align reads to reference genome.

        Args:
            fastq: Rescued FASTQ from Stage 3.

        Returns:
            Path to sorted, indexed BAM.
        """
        logger.info("=" * 55)
        logger.info("Stage 4/5: Minimap2 — splice-aware alignment")
        logger.info("=" * 55)

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

    # ---- Stage 5: Homopolymer Filter ----

    def _run_filter(self, bam: Path, fastq: Path) -> FilterReport:
        """Stage 5: detect and remove homopolymer-mediated artifacts.

        Args:
            bam:   Aligned BAM from Stage 4.
            fastq: The FASTQ that was aligned (for splitting).

        Returns:
            FilterReport with statistics.
        """
        logger.info("=" * 55)
        logger.info("Stage 5/5: Homopolymer Filter — RT artifact removal")
        logger.info("=" * 55)

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
            1. Breakinator: remove foldback inversion artifacts.
            2. Restrander: correct strand orientation, trim primers.
            3. Rescuer: chop reads with internal TSO/RTP adapters.
            4. Minimap2: splice-aware alignment to reference.
            5. Homopolymer Filter: remove RT template switching artifacts.

        Returns:
            PipelineReport with combined statistics from all stages.
        """
        start_time = time.time()
        pipeline_report = PipelineReport()

        logger.info(f"DirectClean v{self._version} starting")
        logger.info(f"  Input:     {self.input_fastq}")
        logger.info(f"  Reference: {self.reference}")
        logger.info(f"  Output:    {self.output_dir}")
        logger.info(f"  Threads:   {self.config.threads}")

        # Stage 1: Breakinator
        stage1_fastq, break_report = self._run_breakinator()
        pipeline_report.break_report = break_report

        # Stage 2: Restrander
        stage2_fastq, restrander_report = self._run_restrander(stage1_fastq)
        pipeline_report.restrander_report = restrander_report

        # Stage 3: Rescuer
        stage3_fastq, rescue_report = self._run_rescue(stage2_fastq)
        pipeline_report.rescue_report = rescue_report

        # Stage 4: Minimap2
        bam = self._run_alignment(stage3_fastq)

        # Stage 5: Homopolymer Filter
        filter_report = self._run_filter(bam, stage3_fastq)
        pipeline_report.filter_report = filter_report

        pipeline_report.elapsed_seconds = time.time() - start_time

        logger.info(str(pipeline_report))

        return pipeline_report

    @property
    def _version(self) -> str:
        try:
            from directclean import __version__
            return __version__
        except ImportError:
            return "unknown"
