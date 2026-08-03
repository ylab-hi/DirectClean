"""
Minimap2 alignment wrapper for DirectClean.

Runs minimap2 with Direct-cDNA-optimised parameters and converts the
SAM stream directly to an alignment-order BAM with samtools view.
Coordinate sorting and BAM indexing are intentionally skipped because
the downstream homopolymer classifier scans the complete BAM
sequentially and does not perform genomic-region queries.

Typical usage::

    from directclean.external.minimap2 import Minimap2Aligner

    aligner = Minimap2Aligner(
        reference="genome.fa",
        threads=8,
    )
    bam_path = aligner.align(
        fastq="reads.fastq",
        output_bam="results/aligned.bam",
    )
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check_binary(name: str) -> str:
    """Verify that an external binary is on PATH.

    Args:
        name: Binary name (e.g. "minimap2", "samtools").

    Returns:
        Absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary cannot be found.
    """
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(
            f"'{name}' not found on PATH. "
            f"Please install {name} and ensure it is accessible."
        )
    return path


# ---------------------------------------------------------------------------
# Aligner class
# ---------------------------------------------------------------------------


class Minimap2Aligner:
    """Minimap2 wrapper with Direct-cDNA defaults.

    Parameters are baked in for the Oxford Nanopore Direct-cDNA
    protocol but can be overridden where needed.

    Args:
        reference:    Path to the reference genome FASTA.
        threads:      Number of threads for minimap2 and samtools.
        extra_args:   Additional minimap2 arguments (list of strings).
        sample_id:    Sample name for the @RG read-group tag.
    """

    # Default minimap2 flags for Direct-cDNA
    _DEFAULTS = [
        "-Y",  # soft-clip with original seq (critical for CIGAR parsing)
        "-ax",
        "splice",  # spliced alignment for cDNA
        "-uf",  # transcript on forward strand (Direct-cDNA protocol)
        "-k14",  # shorter kmer, better sensitivity for Nanopore error profile
        "--secondary=no",  # drop secondary alignments (we only need primary + suppl)
        "--cs",  # cs tag for debugging / variant calling
    ]

    def __init__(
        self,
        reference: str | Path,
        threads: int = 4,
        extra_args: list[str] | None = None,
        sample_id: str = "directclean",
    ) -> None:
        self.reference = Path(reference)
        self.threads = threads
        self.extra_args = extra_args or []
        self.sample_id = sample_id

        # Validate inputs early
        if not self.reference.exists():
            raise FileNotFoundError(f"Reference not found: {self.reference}")
        _check_binary("minimap2")
        _check_binary("samtools")

    # ---- public API ----

    def align(
        self,
        fastq: str | Path,
        output_bam: str | Path,
    ) -> Path:
        """Align FASTQ to reference and produce an alignment-order BAM.

        Pipeline::

            minimap2 ... ref.fa reads.fq | samtools view -b → aligned.bam

        The BAM is intentionally not coordinate-sorted or indexed.
        DirectClean scans every alignment sequentially and does not perform
        genomic-region queries, so sorting and indexing add unnecessary
        memory, runtime, temporary-file, and I/O costs.

        Args:
            fastq:      Input FASTQ file (plain or gzipped).
            output_bam: Path for the output BAM.

        Returns:
            Path to the alignment-order BAM file.

        Raises:
            FileNotFoundError:              If input files are missing.
            subprocess.CalledProcessError:  If minimap2 or samtools fails.
        """
        fastq = Path(fastq)
        output_bam = Path(output_bam)

        if not fastq.exists():
            raise FileNotFoundError(f"FASTQ not found: {fastq}")

        output_bam.parent.mkdir(parents=True, exist_ok=True)

        self._align_to_bam(fastq, output_bam)

        logger.info(f"Alignment complete: {output_bam}")
        return output_bam

    # ---- internal steps ----

    def _build_minimap2_cmd(self, fastq: Path) -> list[str]:
        """Construct the minimap2 command."""
        rg_string = (
            f"@RG\\tID:{self.sample_id}_direct_cDNA"
            f"\\tSM:{self.sample_id}\\tLB:lib\\tPL:ONT"
        )

        cmd = [
            "minimap2",
            *self._DEFAULTS,
            "-t",
            str(self.threads),
            "-R",
            rg_string,
            *self.extra_args,
            str(self.reference),
            str(fastq),
        ]
        return cmd

    def _build_samtools_view_cmd(self, output_bam: Path) -> list[str]:
        """Construct the samtools view command for SAM-to-BAM conversion."""
        return [
            "samtools",
            "view",
            "-@",
            str(self.threads),
            "-b",
            "-o",
            str(output_bam),
            "-",  # read SAM from stdin
        ]

    def _align_to_bam(self, fastq: Path, output_bam: Path) -> None:
        """Run minimap2 piped directly into samtools view.

        The output preserves minimap2 alignment order. It is not
        coordinate-sorted and does not require a BAM index for DirectClean's
        sequential downstream scan.
        """
        mm2_cmd = self._build_minimap2_cmd(fastq)
        view_cmd = self._build_samtools_view_cmd(output_bam)

        logger.info(f"[minimap2] {' '.join(mm2_cmd)}")
        logger.info(f"[samtools view] {' '.join(view_cmd)}")

        # minimap2 SAM stdout → samtools view BAM conversion
        mm2_proc = subprocess.Popen(
            mm2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        view_proc = subprocess.Popen(
            view_cmd,
            stdin=mm2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if mm2_proc.stdout is None or mm2_proc.stderr is None:
            raise RuntimeError("Failed to create minimap2 subprocess pipes")

        # Allow minimap2 to receive SIGPIPE if samtools exits early.
        mm2_proc.stdout.close()

        _, view_stderr = view_proc.communicate()
        mm2_stderr = mm2_proc.stderr.read()
        mm2_proc.wait()

        if mm2_stderr:
            for line in mm2_stderr.decode().strip().split("\n"):
                logger.debug(f"  [minimap2] {line}")

        if mm2_proc.returncode != 0:
            logger.error(f"minimap2 failed (exit {mm2_proc.returncode})")
            logger.error(mm2_stderr.decode())
            raise subprocess.CalledProcessError(
                mm2_proc.returncode,
                mm2_cmd,
                stderr=mm2_stderr,
            )

        if view_proc.returncode != 0:
            logger.error(f"samtools view failed (exit {view_proc.returncode})")
            logger.error(view_stderr.decode())
            raise subprocess.CalledProcessError(
                view_proc.returncode,
                view_cmd,
                stderr=view_stderr,
            )

