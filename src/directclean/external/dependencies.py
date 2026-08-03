"""
External dependency checker for DirectClean.

Validates that all required external binaries (minimap2, samtools,
breakinator, restrander) are available on ``$PATH`` before the
pipeline starts.

For Breakinator, DirectClean also validates the required command-line
capabilities. This is necessary because different Bioconda builds may
share the same package version while providing incompatible interfaces.

Fail-fast design: raises immediately with a clear message telling the
user what is missing or incompatible and how to install it.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Registry of required tools
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExternalTool:
    """Metadata for a required external binary.

    Attributes:
        name:            Binary name to search on PATH.
        description:     One-line description shown in error messages.
        install_hint:    How the user can install it.
        required:        Whether the pipeline cannot run without it.
    """

    name: str
    description: str
    install_hint: str
    required: bool = True


# All tools DirectClean depends on
REQUIRED_TOOLS: list[ExternalTool] = [
    ExternalTool(
        name="minimap2",
        description="Long-read splice-aware aligner",
        install_hint="conda install -c bioconda minimap2",
    ),
    ExternalTool(
        name="samtools",
        description="SAM/BAM manipulation toolkit",
        install_hint="conda install -c bioconda samtools",
    ),
    ExternalTool(
        name="breakinator",
        description="Foldback / inversion artifact detector",
        install_hint=(
            "mamba install -c bioconda "
            "'breakinator=1.1.1=h067a5f5_1'"
        ),
    ),
    ExternalTool(
        name="restrander",
        description="Direct-cDNA strand orientation corrector",
        install_hint="conda install -c genomedk restrander",
    ),
]


# ---------------------------------------------------------------------------
# Breakinator compatibility
# ---------------------------------------------------------------------------


def check_breakinator_compatibility(binary: str) -> None:
    """Verify that Breakinator provides the interface DirectClean requires.

    DirectClean requires the Rust Breakinator implementation, which:

    - accepts SAM/BAM/CRAM input;
    - supports ``--threads``;
    - supports ``--tabular``.

    Compatibility is checked using the actual command-line capabilities,
    rather than the reported version string. Different Conda builds may
    share the same package version while exposing incompatible interfaces,
    and the compatible binary may report a version string that differs
    from the Conda package version.

    Args:
        binary: Absolute path to the Breakinator executable.

    Raises:
        RuntimeError: If Breakinator cannot be inspected or lacks one or
            more capabilities required by DirectClean.
    """

    install_command = (
        "mamba install -c bioconda "
        "'breakinator=1.1.1=h067a5f5_1'"
    )

    try:
        proc = subprocess.run(
            [binary, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            "Timed out while checking the Breakinator installation.\n"
            f"Binary: {binary}\n\n"
            "DirectClean requires the compatible Rust Breakinator build.\n"
            "Install it with:\n"
            f"    {install_command}"
        ) from exc
    except OSError as exc:
        raise RuntimeError(
            "Unable to inspect the Breakinator installation.\n"
            f"Binary: {binary}\n"
            f"Error: {exc}\n\n"
            "DirectClean requires the compatible Rust Breakinator build.\n"
            "Install it with:\n"
            f"    {install_command}"
        ) from exc

    help_text = f"{proc.stdout}\n{proc.stderr}".lower()

    required_capabilities = {
        "SAM/BAM/CRAM input support": "sam/bam/cram",
        "--threads option": "--threads",
        "--tabular option": "--tabular",
    }

    missing_capabilities = [
        description
        for description, marker in required_capabilities.items()
        if marker not in help_text
    ]

    if proc.returncode != 0 or missing_capabilities:
        missing_text = (
            ", ".join(missing_capabilities)
            if missing_capabilities
            else "valid --help output"
        )

        raise RuntimeError(
            "An incompatible Breakinator installation was found.\n"
            f"Binary: {binary}\n"
            f"Missing capability: {missing_text}\n\n"
            "DirectClean requires the Rust Breakinator implementation "
            "with SAM/BAM/CRAM input support and the --threads and "
            "--tabular options.\n\n"
            "Install the compatible Bioconda build with:\n"
            f"    {install_command}"
        )

    logger.debug(f"  ✓ Breakinator interface is compatible: {binary}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def check_binary(name: str) -> str:
    """Check that a single binary is available on PATH.

    Args:
        name: Binary name (e.g. ``"minimap2"``).

    Returns:
        Absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary is not found.
    """

    path = shutil.which(name)

    if path is None:
        raise FileNotFoundError(
            f"Required tool '{name}' not found on PATH. "
            "Please install it first."
        )

    return path


def check_all_dependencies() -> dict[str, str]:
    """Validate all required external tools are available and compatible.

    Returns:
        Dictionary mapping tool name to absolute path.

    Raises:
        FileNotFoundError: Lists all missing tools.
        RuntimeError: If Breakinator is installed but incompatible.
    """

    found: dict[str, str] = {}
    missing: list[ExternalTool] = []

    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool.name)

        if path is not None:
            found[tool.name] = path
            logger.debug(f"  ✓ {tool.name}: {path}")
        elif tool.required:
            missing.append(tool)

    if missing:
        lines = ["The following required tools are missing:\n"]

        for tool in missing:
            lines.append(f"  • {tool.name} — {tool.description}")
            lines.append(f"    Install: {tool.install_hint}\n")

        lines.append(
            "Tip: use the provided environment.yml to install everything:\n"
            "    mamba env create -f environment.yml"
        )

        raise FileNotFoundError("\n".join(lines))

    # Different Breakinator builds may share the same package version while
    # exposing incompatible command-line interfaces. Check the real binary
    # before any expensive pipeline stage begins.
    check_breakinator_compatibility(found["breakinator"])

    logger.info("All external dependencies found.")
    return found