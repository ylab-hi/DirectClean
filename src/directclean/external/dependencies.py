"""
External dependency checker for DirectClean.

Checks required external command-line tools before running the pipeline.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExternalTool:
    name: str
    description: str
    install_hint: str
    required: bool = True


REQUIRED_TOOLS: list[ExternalTool] = [
    ExternalTool(
        name="minimap2",
        description="Long-read splice-aware aligner",
        install_hint=(
            "mamba install -c conda-forge -c bioconda minimap2"
        ),
    ),
    ExternalTool(
        name="samtools",
        description="SAM/BAM manipulation toolkit",
        install_hint=(
            "mamba install -c conda-forge -c bioconda samtools"
        ),
    ),
    ExternalTool(
        name="breakinator",
        description="Foldback / inversion artifact detector",
        install_hint=(
            "mamba install -c conda-forge -c bioconda "
            "'breakinator=1.1.1=*_2'"
        ),
    ),
    ExternalTool(
        name="restrander",
        description="Direct-cDNA strand orientation corrector",
        install_hint=(
            "mamba install -c conda-forge -c bioconda "
            "'restrander=1.1.3=*_1'"
        ),
    ),
]


def check_breakinator_compatibility(binary: str) -> None:
    """
    Check that Breakinator provides the interface required by DirectClean.

    Required:
      - SAM/BAM/CRAM input support
      - --threads
      - --tabular
    """

    install_command = (
        "mamba install -c conda-forge -c bioconda "
        "'breakinator=1.1.1=*_2'"
    )

    try:
        proc = subprocess.run(
            [binary, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
        )
    except Exception as exc:
        raise RuntimeError(
            "Unable to inspect Breakinator installation.\n"
            f"Binary: {binary}\n\n"
            "Install the compatible Bioconda build with:\n"
            f"    {install_command}"
        ) from exc

    help_text = f"{proc.stdout}\n{proc.stderr}".lower()

    required = {
        "SAM/BAM/CRAM input support": "sam/bam/cram",
        "--threads option": "--threads",
        "--tabular option": "--tabular",
    }

    missing = [
        name for name, marker in required.items()
        if marker not in help_text
    ]

    if proc.returncode != 0 or missing:
        missing_text = ", ".join(missing) if missing else "valid help output"

        raise RuntimeError(
            "An incompatible Breakinator installation was found.\n"
            f"Binary: {binary}\n"
            f"Missing capability: {missing_text}\n\n"
            "DirectClean requires the Rust Breakinator implementation "
            "with SAM/BAM/CRAM input support and --threads/--tabular.\n\n"
            "Install the compatible build with:\n"
            f"    {install_command}"
        )

    logger.debug(
        "Breakinator interface is compatible: %s",
        binary,
    )


def check_binary(name: str) -> str:
    """Check that a binary exists on PATH."""

    path = shutil.which(name)

    if path is None:
        raise FileNotFoundError(
            f"Required tool '{name}' not found on PATH."
        )

    return path


def check_all_dependencies() -> dict[str, str]:
    """
    Validate all DirectClean external dependencies.

    Returns:
        Mapping from tool name to executable path.
    """

    found: dict[str, str] = {}
    missing: list[ExternalTool] = []

    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool.name)

        if path:
            found[tool.name] = path
            logger.debug("%s: %s", tool.name, path)
        elif tool.required:
            missing.append(tool)

    if missing:
        lines = [
            "The following required tools are missing:",
            "",
        ]

        for tool in missing:
            lines.append(
                f"  • {tool.name} — {tool.description}"
            )
            lines.append(
                f"    Install: {tool.install_hint}"
            )
            lines.append("")

        lines.extend(
            [
                "For source installation, install dependencies with:",
                "    mamba env create -f environment.yml",
            ]
        )

        raise FileNotFoundError("\n".join(lines))

    check_breakinator_compatibility(
        found["breakinator"]
    )

    logger.info(
        "All external dependencies found."
    )

    return found