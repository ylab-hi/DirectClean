"""
Breakinator wrapper — remove foldback inversion artifacts.

Breakinator detects two types of artifacts in Direct-cDNA reads:
*foldback inversions* (hairpin structures from sequencing) and
*chimeric reads* (ligation artifacts).  DirectClean only removes
**foldback** reads because chimeric reads are handled downstream
by our own Rescuer and Homopolymer Filter modules.

The workflow mirrors the manual steps::

    1. minimap2  → raw name-sorted SAM  (Breakinator needs SAM input)
    2. breakinator → artifacts.txt       (tabular classification)
    3. Parse artifacts.txt → collect Foldback read IDs
    4. split_fastq_by_ids → keep non-foldback reads

Note: This minimap2 run is independent from Stage 4 (the alignment
used by the Homopolymer Filter).  Breakinator requires a SAM file,
not a coordinate-sorted BAM.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Set

from directclean.external.dependencies import check_binary
from directclean.utils.io import split_fastq_by_ids

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


@dataclass
class BreakReport:
    """Summary statistics from the Breakinator stage.

    Attributes:
        total_breakpoints:  Total breakpoints classified by Breakinator.
        foldback_count:     Breakpoints classified as Foldback.
        chimeric_count:     Breakpoints classified as Chimeric.
        pass_count:         Breakpoints classified as Pass.
        foldback_read_ids:  Number of unique reads with Foldback.
        input_reads:        Total reads in input FASTQ.
        kept_reads:         Reads after removing foldback.
        removed_reads:      Foldback reads removed.
    """

    total_breakpoints: int = 0
    foldback_count: int = 0
    chimeric_count: int = 0
    pass_count: int = 0
    foldback_read_ids: int = 0
    input_reads: int = 0
    kept_reads: int = 0
    removed_reads: int = 0

    def __str__(self) -> str:
        pct = (
            f"{self.removed_reads / self.input_reads * 100:.1f}%"
            if self.input_reads > 0
            else "N/A"
        )
        return (
            "=== Breakinator Report ===\n"
            f"  Breakpoints total       : {self.total_breakpoints:,}\n"
            f"    Foldback              : {self.foldback_count:,}\n"
            f"    Chimeric              : {self.chimeric_count:,}\n"
            f"    Pass                  : {self.pass_count:,}\n"
            f"  ---\n"
            f"  Input reads             : {self.input_reads:,}\n"
            f"  Foldback reads removed  : {self.removed_reads:,} ({pct})\n"
            f"  Reads kept              : {self.kept_reads:,}\n"
            "=========================="
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _run_minimap2_for_breakinator(
    reference: Path,
    input_fastq: Path,
    output_sam: Path,
    threads: int,
    junc_bed: Optional[Path] = None,
) -> None:
    """Run minimap2 to produce a SAM file for Breakinator.

    This is NOT coordinate-sorted — Breakinator reads SAM directly.
    Uses the same Direct-cDNA splice parameters as the main alignment
    but outputs SAM instead of piping to samtools sort.

    Args:
        reference:    Reference genome FASTA.
        input_fastq:  Raw input FASTQ.
        output_sam:   Output SAM path.
        threads:      Number of threads.
        junc_bed:     Optional junction BED for guided alignment.
    """
    minimap2_bin = check_binary("minimap2")

    cmd = [
        minimap2_bin,
        "-Y",  # soft-clip with original sequence
        "-t",
        str(threads),
        "-ax",
        "splice",  # splice-aware
        "-uf",  # forward strand for Direct-cDNA
        "-k14",  # k-mer size
        "--secondary=no",  # no secondary alignments
    ]

    if junc_bed is not None:
        cmd.extend(["--junc-bed", str(junc_bed)])

    cmd.extend([str(reference), str(input_fastq)])

    logger.info(f"Running minimap2 for Breakinator: {' '.join(cmd[:6])}...")

    with open(output_sam, "w") as sam_fh:
        proc = subprocess.run(
            cmd,
            stdout=sam_fh,
            stderr=subprocess.PIPE,
            text=True,
        )

    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed (exit {proc.returncode}):\n{proc.stderr}")

    logger.info(f"SAM written: {output_sam}")


def _run_breakinator(
    input_sam: Path,
    output_artifacts: Path,
    threads: int,
) -> None:
    """Run Breakinator on a SAM file.

    Args:
        input_sam:         Input SAM from minimap2.
        output_artifacts:  Output tabular artifacts file.
        threads:           Number of threads.
    """
    breakinator_bin = check_binary("breakinator")

    cmd = [
        breakinator_bin,
        "-i",
        str(input_sam),
        "-o",
        str(output_artifacts),
        "--tabular",
    ]

    logger.info(f"Running Breakinator: {' '.join(cmd[:4])}...")

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Breakinator failed (exit {proc.returncode}):\n{proc.stderr}"
        )

    logger.info(f"Artifacts written: {output_artifacts}")


def _parse_foldback_ids(artifacts_file: Path) -> tuple[Set[str], dict]:
    """Parse Breakinator tabular output and extract Foldback read IDs.

    Only Foldback reads are removed — Chimeric reads are left for
    DirectClean's own Rescuer and Homopolymer Filter to handle.

    Args:
        artifacts_file: Path to Breakinator --tabular output.

    Returns:
        (foldback_ids, breakpoint_stats) where breakpoint_stats has
        counts for each classification.
    """
    foldback_ids: Set[str] = set()
    stats = {"total": 0, "Foldback": 0, "Chimeric": 0, "Pass": 0}

    with open(artifacts_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            read_id = fields[6]
            classification = fields[7]

            stats["total"] += 1
            if classification in stats:
                stats[classification] += 1

            if classification == "Foldback":
                foldback_ids.add(read_id)

    logger.info(
        f"Breakinator: {stats['total']:,} breakpoints, "
        f"{len(foldback_ids):,} unique foldback reads"
    )
    return foldback_ids, stats


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


class BreakinatorRunner:
    """End-to-end Breakinator wrapper for foldback removal.

    Runs minimap2 → breakinator → foldback filtering as a single
    stage in the DirectClean pipeline.

    Usage::

        runner = BreakinatorRunner(
            reference=Path("genome.fa"),
            threads=8,
        )
        report = runner.run(
            input_fastq=Path("raw.fastq"),
            output_fastq=Path("no_foldback.fastq"),
            work_dir=Path("output/"),
        )

    Args:
        reference:  Reference genome FASTA.
        threads:    Number of threads for minimap2 and breakinator.
        junc_bed:   Optional junction BED file for guided alignment.
    """

    def __init__(
        self,
        reference: Path,
        threads: int = 4,
        junc_bed: Optional[Path] = None,
    ) -> None:
        self.reference = Path(reference)
        self.threads = threads
        self.junc_bed = Path(junc_bed) if junc_bed is not None else None

    def run(
        self,
        input_fastq: Path,
        output_fastq: Path,
        work_dir: Path,
        prefix: str = "directclean",
    ) -> BreakReport:
        """Execute the full Breakinator stage.

        Steps:
            1. minimap2 alignment → SAM
            2. Breakinator → artifacts.txt
            3. Parse foldback IDs
            4. Split FASTQ → kept + removed

        Args:
            input_fastq:  Raw input FASTQ.
            output_fastq: Output FASTQ with foldback reads removed.
            work_dir:     Working directory for intermediate files.
            prefix:       Filename prefix for intermediates.

        Returns:
            BreakReport with statistics.
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        report = BreakReport()

        # Intermediate file paths
        sam_path = work_dir / f"{prefix}.breakinator.sam"
        artifacts_path = work_dir / f"{prefix}.breakinator_artifacts.txt"
        removed_path = work_dir / f"{prefix}.foldback_removed.fastq"

        # Step 1: minimap2 → SAM
        logger.info("Step 1/3: Aligning reads for Breakinator...")
        _run_minimap2_for_breakinator(
            reference=self.reference,
            input_fastq=input_fastq,
            output_sam=sam_path,
            threads=self.threads,
            junc_bed=self.junc_bed,
        )

        # Step 2: Breakinator → artifacts.txt
        logger.info("Step 2/3: Running Breakinator...")
        _run_breakinator(
            input_sam=sam_path,
            output_artifacts=artifacts_path,
            threads=self.threads,
        )

        # Step 3: Parse foldback IDs and filter FASTQ
        logger.info("Step 3/3: Filtering foldback reads...")
        foldback_ids, bp_stats = _parse_foldback_ids(artifacts_path)

        report.total_breakpoints = bp_stats["total"]
        report.foldback_count = bp_stats["Foldback"]
        report.chimeric_count = bp_stats["Chimeric"]
        report.pass_count = bp_stats["Pass"]
        report.foldback_read_ids = len(foldback_ids)

        # Use our own split_fastq_by_ids — no seqkit dependency
        kept, removed = split_fastq_by_ids(
            input_fastq=input_fastq,
            remove_ids=foldback_ids,
            output_kept=output_fastq,
            output_removed=removed_path,
        )

        report.input_reads = kept + removed
        report.kept_reads = kept
        report.removed_reads = removed

        # Clean up large SAM file (user still has artifacts.txt for reference)
        if sam_path.exists():
            sam_path.unlink()
            logger.debug(f"Cleaned up intermediate SAM: {sam_path}")

        logger.info(f"Breakinator stage complete.\n{report}")
        return report
