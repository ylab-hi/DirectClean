"""
Minimap2 alignment wrapper for DirectClean.

Runs minimap2 with Direct-cDNA-optimised parameters, pipes output
through samtools sort, and builds a BAM index.  All external calls
go through a single helper so that logging, error handling, and
dry-run support live in one place.

Typical usage::

    from directclean.external.minimap2 import Minimap2Aligner

    aligner = Minimap2Aligner(
        reference="genome.fa",
        threads=8,
    )
    bam_path = aligner.align(
        fastq="reads.fastq",
        output_bam="results/aligned.sorted.bam",
    )
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

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


def _run(cmd: list[str], description: str, **kwargs) -> subprocess.CompletedProcess:
    """Run a command with logging and error handling.

    Args:
        cmd:         Command as a list of strings.
        description: Human-readable description for log messages.
        **kwargs:    Forwarded to subprocess.run().

    Returns:
        CompletedProcess instance.

    Raises:
        subprocess.CalledProcessError: If the command exits non-zero.
    """
    cmd_str = " ".join(cmd)
    logger.info(f"[{description}] {cmd_str}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        **kwargs,
    )

    if result.returncode != 0:
        logger.error(f"[{description}] failed (exit {result.returncode})")
        logger.error(f"  stderr: {result.stderr.strip()}")
        result.check_returncode()  # raises CalledProcessError

    return result


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
        "-Y",              # soft-clip with original seq (critical for CIGAR parsing)
        "-ax", "splice",   # spliced alignment for cDNA
        "-uf",             # transcript on forward strand (Direct-cDNA protocol)
        "-k14",            # shorter kmer, better sensitivity for Nanopore error profile
        "--secondary=no",  # drop secondary alignments (we only need primary + suppl)
        "--cs",            # cs tag for debugging / variant calling
    ]

    def __init__(
        self,
        reference: str | Path,
        threads: int = 4,
        extra_args: Optional[list[str]] = None,
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
        """Align FASTQ to reference and produce a sorted, indexed BAM.

        Pipeline::

            minimap2 ... ref.fa reads.fq | samtools sort → sorted.bam
            samtools index sorted.bam

        Args:
            fastq:      Input FASTQ file (plain or gzipped).
            output_bam: Path for the output sorted BAM.

        Returns:
            Path to the sorted, indexed BAM file.

        Raises:
            FileNotFoundError:              If input files are missing.
            subprocess.CalledProcessError:  If minimap2 or samtools fails.
        """
        fastq = Path(fastq)
        output_bam = Path(output_bam)

        if not fastq.exists():
            raise FileNotFoundError(f"FASTQ not found: {fastq}")

        # Create output directory
        output_bam.parent.mkdir(parents=True, exist_ok=True)

        # --- Step 1: minimap2 | samtools sort ---
        self._align_and_sort(fastq, output_bam)

        # --- Step 2: samtools index ---
        self._index(output_bam)

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
            "-t", str(self.threads),
            "-R", rg_string,
            *self.extra_args,
            str(self.reference),
            str(fastq),
        ]
        return cmd

    def _build_samtools_sort_cmd(self, output_bam: Path) -> list[str]:
        """Construct the samtools sort command."""
        return [
            "samtools", "sort",
            "-@", str(self.threads),
            "-O", "BAM",
            "-o", str(output_bam),
            "-",  # read from stdin
        ]

    def _align_and_sort(self, fastq: Path, output_bam: Path) -> None:
        """Run minimap2 piped into samtools sort.

        Uses two subprocesses connected by a pipe, avoiding a
        temporary SAM file that could be hundreds of GB.
        """
        mm2_cmd = self._build_minimap2_cmd(fastq)
        sort_cmd = self._build_samtools_sort_cmd(output_bam)

        logger.info(f"[minimap2] {' '.join(mm2_cmd)}")
        logger.info(f"[samtools sort] {' '.join(sort_cmd)}")

        # minimap2 stdout → samtools sort stdin
        mm2_proc = subprocess.Popen(
            mm2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        sort_proc = subprocess.Popen(
            sort_cmd,
            stdin=mm2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # Allow mm2 to receive SIGPIPE if sort exits early
        mm2_proc.stdout.close()

        # Wait for both processes
        sort_stdout, sort_stderr = sort_proc.communicate()
        mm2_stderr = mm2_proc.stderr.read()
        mm2_proc.wait()

        # Log minimap2 stats (written to stderr)
        if mm2_stderr:
            for line in mm2_stderr.decode().strip().split("\n"):
                logger.debug(f"  [minimap2] {line}")

        # Check exit codes
        if mm2_proc.returncode != 0:
            logger.error(f"minimap2 failed (exit {mm2_proc.returncode})")
            logger.error(mm2_stderr.decode())
            raise subprocess.CalledProcessError(
                mm2_proc.returncode, mm2_cmd, stderr=mm2_stderr
            )

        if sort_proc.returncode != 0:
            logger.error(f"samtools sort failed (exit {sort_proc.returncode})")
            logger.error(sort_stderr.decode())
            raise subprocess.CalledProcessError(
                sort_proc.returncode, sort_cmd, stderr=sort_stderr
            )

    def _index(self, bam_path: Path) -> None:
        """Run samtools index."""
        _run(
            ["samtools", "index", str(bam_path)],
            description="samtools index",
        )