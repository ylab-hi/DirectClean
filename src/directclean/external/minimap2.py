"""
Minimap2 alignment wrapper for DirectClean.

Builds or reuses a persistent splice-aware minimap2 index, then runs
minimap2 with Direct-cDNA-optimised parameters and converts the SAM
stream directly to an alignment-order BAM with samtools view.

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

import hashlib
import logging
import os
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
# Persistent reference index
# ---------------------------------------------------------------------------


_INDEX_SIGNATURE = "minimap2-splice-k14-w5-v1"


def _default_index_cache_dir() -> Path:
    """Return the persistent DirectClean minimap2 cache directory.

    The location can be overridden with ``DIRECTCLEAN_CACHE_DIR``.
    """
    configured = os.environ.get("DIRECTCLEAN_CACHE_DIR")
    if configured:
        return Path(configured).expanduser() / "minimap2"
    return Path.home() / ".cache" / "directclean" / "minimap2"


def _reference_cache_key(reference: Path) -> str:
    """Build a stable cache key for one reference file and index recipe."""
    resolved = reference.resolve()
    stat = resolved.stat()
    payload = (
        f"{resolved}\0{stat.st_size}\0{stat.st_mtime_ns}\0{_INDEX_SIGNATURE}"
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def get_or_build_minimap2_index(
    reference: str | Path,
    threads: int = 1,
    cache_dir: str | Path | None = None,
) -> Path:
    """Return a persistent splice-aware ``k=14`` minimap2 index.

    A supplied ``.mmi`` is returned unchanged.  For a FASTA reference,
    DirectClean stores the index in a user-writable persistent cache so
    all samples and both alignment stages can reuse it.

    The cache key includes the resolved reference path, file size,
    modification time, and index recipe.  A changed reference therefore
    receives a new index automatically.

    Args:
        reference: Reference FASTA or an existing ``.mmi`` index.
        threads: Threads used while building a missing index.
        cache_dir: Optional cache root.  Defaults to
            ``$DIRECTCLEAN_CACHE_DIR/minimap2`` when configured, otherwise
            ``~/.cache/directclean/minimap2``.

    Returns:
        Path to the existing or newly built ``.mmi`` file.

    Raises:
        FileNotFoundError: If the reference does not exist.
        subprocess.CalledProcessError: If minimap2 index construction fails.
        RuntimeError: If minimap2 reports success without producing an index.
    """
    reference = Path(reference).expanduser()
    if not reference.exists():
        raise FileNotFoundError(f"Reference not found: {reference}")

    if reference.suffix == ".mmi":
        if reference.stat().st_size == 0:
            raise RuntimeError(f"Minimap2 index is empty: {reference}")
        logger.info(f"Using supplied minimap2 index: {reference}")
        return reference

    cache_root = (
        Path(cache_dir).expanduser()
        if cache_dir is not None
        else _default_index_cache_dir()
    )
    cache_root.mkdir(parents=True, exist_ok=True)

    key = _reference_cache_key(reference)
    safe_stem = reference.name.replace(os.sep, "_")
    output_index = cache_root / f"{safe_stem}.{key}.splice_k14.mmi"

    if output_index.exists() and output_index.stat().st_size > 0:
        logger.info(f"Reusing persistent minimap2 index: {output_index}")
        return output_index

    minimap2_bin = _check_binary("minimap2")
    temporary_index = output_index.with_name(
        f".{output_index.name}.{os.getpid()}.tmp"
    )
    if temporary_index.exists():
        temporary_index.unlink()

    # The current mapping recipe is `-x splice -k14`.  The splice preset
    # supplies w=5; k=14 overrides the preset's default k while preserving
    # the same index-affecting settings used by both pipeline stages.
    cmd = [
        minimap2_bin,
        "-x",
        "splice",
        "-k14",
        "-t",
        str(max(1, threads)),
        "-d",
        str(temporary_index),
        str(reference),
    ]

    logger.info(f"Building persistent minimap2 index: {output_index}")
    logger.info(f"[minimap2 index] {' '.join(cmd)}")

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if proc.stderr:
        for line in proc.stderr.strip().splitlines():
            logger.debug(f"  [minimap2 index] {line}")

    if proc.returncode != 0:
        if temporary_index.exists():
            temporary_index.unlink()
        raise subprocess.CalledProcessError(
            proc.returncode,
            cmd,
            output=proc.stdout,
            stderr=proc.stderr,
        )

    if not temporary_index.exists() or temporary_index.stat().st_size == 0:
        raise RuntimeError(
            "minimap2 reported a successful index build but no non-empty "
            f"index was produced: {temporary_index}"
        )

    # Atomic within the same filesystem.  If another job completed the same
    # cache entry first, replacing it with an equivalent index is harmless.
    temporary_index.replace(output_index)
    logger.info(f"Persistent minimap2 index ready: {output_index}")
    return output_index


# ---------------------------------------------------------------------------
# Aligner class
# ---------------------------------------------------------------------------


class Minimap2Aligner:
    """Minimap2 wrapper with Direct-cDNA defaults.

    Parameters are baked in for the Oxford Nanopore Direct-cDNA
    protocol but can be overridden where needed.

    Args:
        reference:    Path to the reference genome FASTA or prebuilt .mmi.
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
        """Construct the lightweight SAM-to-BAM conversion command.

        Minimap2 is the main CPU consumer.  Samtools therefore uses at
        most two compression threads and BAM compression level 1 for this
        temporary alignment-order file.
        """
        samtools_threads = min(2, self.threads)

        return [
            "samtools",
            "view",
            "-@",
            str(samtools_threads),
            "-b",
            "-1",
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

