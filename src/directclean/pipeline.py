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
    Stage 5: Homopolymer Rescue
        Detect RT template switching junctions → chop and rescue
        │
        ▼
    Final outputs:
        cleaned.fastq   — all reads passing stages 1-2, plus rescued
                           sub-reads from stages 3 and 5
        rescued.fastq    — sub-reads rescued by homopolymer chopping
        reports/         — per-stage reports
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path

from directclean.external.breakinator import BreakinatorRunner, BreakReport
from directclean.external.restrander import RestranderRunner, RestranderReport
from directclean.external.minimap2 import (
    Minimap2Aligner,
    get_or_build_minimap2_index,
)
from directclean.external.dependencies import check_all_dependencies
from directclean.rescuer import ReadChopper, AdapterConfig, RescueReport
from directclean.rescuer.unknowns_rescuer import UnknownsRescuer, UnknownsRescueReport
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
    context_window: int = 50
    min_mapq: int = 0
    threads: int = 4
    junc_bed: Path | None = None

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
        filter_report:     Report from Homopolymer Rescue (Stage 5).
        elapsed_seconds:   Total wall-clock time.
    """

    break_report: BreakReport | None = None
    restrander_report: RestranderReport | None = None
    unknowns_rescue_report: UnknownsRescueReport | None = None
    rescue_report: RescueReport | None = None
    filter_report: FilterReport | None = None
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

        if self.unknowns_rescue_report is not None:
            parts.append("")
            parts.append(str(self.unknowns_rescue_report))

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
    Breakinator → Restrander → Rescuer → Minimap2 → Homopolymer Rescue.

    Stages 1-2 (Breakinator, Restrander) may **remove** reads that are
    definitively artifactual (foldback inversions, invalid primer configs).

    Stages 3 and 5 (Rescuer, Homopolymer Rescue) never remove reads —
    they **chop** chimeric reads at artifact junctions and rescue the
    flanking sub-reads, increasing the effective read count.

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
        self._minimap2_index: Path | None = None

        # Validate inputs
        if not self.input_fastq.exists():
            raise FileNotFoundError(f"Input FASTQ not found: {self.input_fastq}")
        if not self.reference.exists():
            raise FileNotFoundError(f"Reference not found: {self.reference}")

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
    def _unknowns_fastq(self) -> Path:
        """Restrander's unknowns output (auto-generated by Restrander)."""
        return self._intermediates_dir / f"{self.prefix}-unknowns.restranded.fastq"

    @property
    def _unknowns_rescued_fastq(self) -> Path:
        """Oriented sub-reads rescued from unknowns."""
        return self._intermediates_dir / f"{self.prefix}.unknowns_rescued.fastq"

    @property
    def _merged_fastq(self) -> Path:
        """Merged FASTQ: Restrander output + unknowns rescued reads."""
        return self._intermediates_dir / f"{self.prefix}.merged.fastq"

    @property
    def _rescued_fastq(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.rescued.fastq"

    @property
    def _aligned_bam(self) -> Path:
        return self._intermediates_dir / f"{self.prefix}.aligned.bam"

    @property
    def cleaned_fastq(self) -> Path:
        return self.output_dir / f"{self.prefix}.cleaned.fastq"

    @property
    def rescued_fastq(self) -> Path:
        """Sub-reads rescued by homopolymer chopping (Stage 5)."""
        return self.output_dir / f"{self.prefix}.rescued.fastq"

    # ---- Shared persistent minimap2 index ----

    def _prepare_minimap2_index(self) -> Path:
        """Build once when needed and reuse across stages and samples."""
        index_path = get_or_build_minimap2_index(
            reference=self.reference,
            threads=self.config.threads,
        )
        self._minimap2_index = index_path
        return index_path

    @property
    def _mapping_reference(self) -> Path:
        """Reference index used by both minimap2 alignment stages."""
        if self._minimap2_index is None:
            raise RuntimeError(
                "Minimap2 index has not been prepared before alignment."
            )
        return self._minimap2_index

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
            reference=self._mapping_reference,
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

    # ---- Stage 2.5: Unknowns Rescue ----

    def _run_unknowns_rescue(
        self, restranded_fastq: Path
    ) -> tuple[Path, UnknownsRescueReport]:
        """Stage 2.5: rescue oriented reads from Restrander unknowns.

        Scans the unknowns FASTQ for internal adapter junctions, chops
        reads at those junctions, orients each sub-read using polyA/T
        and primer signals.

        Args:
            restranded_fastq: Normal Restrander output (for reference).

        Returns:
            (rescued_fastq, report) — FASTQ containing only the
            rescued and oriented sub-reads from unknowns.
        """
        logger.info("=" * 55)
        logger.info("Stage 2.5/6: Unknowns Rescue — recovering discarded reads")
        logger.info("=" * 55)

        unknowns_report = UnknownsRescueReport()

        if not self._unknowns_fastq.exists():
            logger.warning(
                f"Unknowns FASTQ not found: {self._unknowns_fastq}. "
                f"Skipping unknowns rescue."
            )
            return self._unknowns_rescued_fastq, unknowns_report

        rescuer = UnknownsRescuer(
            unknowns_fastq=self._unknowns_fastq,
            output_fastq=self._unknowns_rescued_fastq,
            config=self.config.adapter_config,
            min_confidence=self.config.min_confidence,
        )
        unknowns_report = rescuer.run()
        return self._unknowns_rescued_fastq, unknowns_report

    # ---- Merge helper ----

    def _merge_for_alignment(
        self,
        stage3_fastq: Path,
        unknowns_rescued_fastq: Path,
        unknowns_report: UnknownsRescueReport,
    ) -> Path:
        """Merge Stage 3 output with Stage 2.5 rescued reads.

        If no unknowns were rescued, returns Stage 3 output directly
        without creating a merged file.

        Args:
            stage3_fastq:          Output from Stage 3 (adapter rescue).
            unknowns_rescued_fastq: Output from Stage 2.5.
            unknowns_report:       Report to check if any reads rescued.

        Returns:
            Path to FASTQ for Stage 4 input.
        """
        if unknowns_report.output_reads == 0:
            logger.info("No unknowns rescued, skipping merge.")
            return stage3_fastq

        logger.info(
            f"Merging Stage 3 output with {unknowns_report.output_reads:,} "
            f"rescued unknowns reads ..."
        )
        with open(self._merged_fastq, "w") as f_out:
            with open(stage3_fastq) as f_in:
                for line in f_in:
                    f_out.write(line)
            if unknowns_rescued_fastq.exists():
                with open(unknowns_rescued_fastq) as f_in:
                    for line in f_in:
                        f_out.write(line)

        logger.info(f"Merged FASTQ written: {self._merged_fastq}")
        return self._merged_fastq

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
            threads=1,
        )
        report = chopper.run()
        return self._rescued_fastq, report

    # ---- Stage 4: Minimap2 ----

    def _run_alignment(self, fastq: Path) -> Path:
        """Stage 4: align reads to reference genome.

        Args:
            fastq: Rescued FASTQ from Stage 3.

        Returns:
            Path to the alignment-order BAM.
        """
        logger.info("=" * 55)
        logger.info("Stage 4/5: Minimap2 — splice-aware alignment")
        logger.info("=" * 55)

        aligner = Minimap2Aligner(
            reference=self._mapping_reference,
            threads=self.config.threads,
            sample_id=self.prefix,
        )
        bam = aligner.align(
            fastq=fastq,
            output_bam=self._aligned_bam,
        )
        return bam

    # ---- Stage 5: Homopolymer Rescue ----

    def _run_filter(self, bam: Path, fastq: Path) -> FilterReport:
        """Stage 5: detect and rescue homopolymer-mediated artifacts.

        Artifact reads are chopped at homopolymer junctions; flanking
        sub-reads ≥100 bp are rescued into the cleaned output.

        Args:
            bam:   Aligned BAM from Stage 4.
            fastq: The FASTQ that was aligned (for chopping).

        Returns:
            FilterReport with statistics.
        """
        logger.info("=" * 55)
        logger.info("Stage 5/5: Homopolymer Rescue — RT artifact chopping")
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
            2.5. Unknowns Rescue: recover oriented reads from unknowns.
            3. Rescuer: chop reads with internal TSO/RTP adapters.
            4. Minimap2: splice-aware alignment to reference.
            5. Homopolymer Rescue: chop RT template switching artifacts.

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

        # Build once if absent; reuse for both alignment stages and future runs.
        index_path = self._prepare_minimap2_index()
        logger.info(f"  MM2 index: {index_path}")

        # Stage 1: Breakinator
        stage1_fastq, break_report = self._run_breakinator()
        pipeline_report.break_report = break_report

        # Stage 2: Restrander
        stage2_fastq, restrander_report = self._run_restrander(stage1_fastq)
        pipeline_report.restrander_report = restrander_report

        # Stage 2.5: Unknowns Rescue
        unknowns_rescued_fastq, unknowns_report = self._run_unknowns_rescue(
            stage2_fastq
        )
        pipeline_report.unknowns_rescue_report = unknowns_report

        # Stage 3: Rescuer (only on Restrander normal output, not unknowns)
        stage3_fastq, rescue_report = self._run_rescue(stage2_fastq)
        pipeline_report.rescue_report = rescue_report

        # Merge Stage 3 output + Stage 2.5 rescued reads → Stage 4 input
        stage4_input = self._merge_for_alignment(
            stage3_fastq, unknowns_rescued_fastq, unknowns_report
        )

        # Stage 4: Minimap2
        bam = self._run_alignment(stage4_input)

        # Stage 5: Homopolymer Rescue
        filter_report = self._run_filter(bam, stage4_input)
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
