"""
External dependency checker for DirectClean.

Validates that all required external binaries (minimap2, samtools,
breakinator, restrander) are available on ``$PATH`` before the
pipeline starts.  Fail-fast design: raises immediately with a
clear message telling the user what is missing and how to install it.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from typing import Dict, List, Optional

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
REQUIRED_TOOLS: List[ExternalTool] = [
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
        install_hint="conda install -c bioconda breakinator",
    ),
    ExternalTool(
        name="restrander",
        description="Direct-cDNA strand orientation corrector",
        install_hint="conda install -c genomedk restrander",
    ),
]


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
            f"Please install it first."
        )
    return path


def check_all_dependencies() -> Dict[str, str]:
    """Validate all required external tools are available.

    Returns:
        Dictionary mapping tool name → absolute path.

    Raises:
        FileNotFoundError: Lists *all* missing tools (not just the
            first one) so the user can fix everything at once.
    """
    found: Dict[str, str] = {}
    missing: List[ExternalTool] = []

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

    logger.info("All external dependencies found.")
    return found
