"""
Restrander wrapper — correct strand orientation of Direct-cDNA reads.

Restrander classifies each read as forward or reverse based on polyA
tail position and TSO/RTP primer sequences, then reverse-complements
reverse reads so all output reads are in 5'→3' orientation.  It also
trims primer sequences and removes reads that cannot be classified
(unknowns) or have aberrant primer configurations (RTP-RTP / TSO-TSO
artefacts).

DirectClean bundles the PCB109 configuration file (for SQK-LSK114 kit)
and calls Restrander as a subprocess.

Restrander is MIT-licensed: https://github.com/Oshlack/restrander
"""

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from directclean.external.dependencies import check_binary

logger = logging.getLogger(__name__)

# Path to the bundled PCB109 config file
_CONFIGS_DIR = Path(__file__).parent / "configs"
DEFAULT_CONFIG = _CONFIGS_DIR / "PCB109.json"


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

@dataclass
class RestranderReport:
    """Summary statistics from the Restrander stage.

    Attributes:
        total_input:      Total reads in input FASTQ.
        forward:          Reads classified as forward (already 5'→3').
        reverse:          Reads classified as reverse (flipped to 5'→3').
        unknown:          Reads that could not be classified (excluded).
        artefacts:        Reads with aberrant primer configs (excluded).
        output_reads:     Reads in the output FASTQ (forward + reverse).
    """
    total_input: int = 0
    forward: int = 0
    reverse: int = 0
    unknown: int = 0
    artefacts: int = 0
    output_reads: int = 0

    def __str__(self) -> str:
        pct_kept = (
            f"{self.output_reads / self.total_input * 100:.1f}%"
            if self.total_input > 0 else "N/A"
        )
        pct_artefact = (
            f"{self.artefacts / self.total_input * 100:.1f}%"
            if self.total_input > 0 else "N/A"
        )
        return (
            "=== Restrander Report ===\n"
            f"  Total input reads       : {self.total_input:,}\n"
            f"  Forward (5'→3')         : {self.forward:,}\n"
            f"  Reverse (flipped)       : {self.reverse:,}\n"
            f"  Unknown (excluded)      : {self.unknown:,}\n"
            f"  Artefacts (excluded)    : {self.artefacts:,} ({pct_artefact})\n"
            f"  ---\n"
            f"  Output reads            : {self.output_reads:,} ({pct_kept})\n"
            "========================="
        )


# ---------------------------------------------------------------------------
# Stats parser
# ---------------------------------------------------------------------------

def _parse_restrander_output(raw_output: str) -> RestranderReport:
    """Parse Restrander's JSON statistics output.

    Restrander writes a JSON object to stdout with classification
    counts.  The format varies slightly between versions, so we
    parse defensively.

    Args:
        raw_output: Combined stdout+stderr from Restrander.

    Returns:
        RestranderReport with parsed statistics.
    """
    report = RestranderReport()

    # Try to find and parse JSON in the output
    # Restrander may print log messages before the JSON
    json_str = None
    brace_depth = 0
    json_start = -1

    for i, ch in enumerate(raw_output):
        if ch == "{":
            if brace_depth == 0:
                json_start = i
            brace_depth += 1
        elif ch == "}":
            brace_depth -= 1
            if brace_depth == 0 and json_start >= 0:
                json_str = raw_output[json_start:i + 1]
                break

    if json_str is None:
        logger.warning(
            "Could not parse Restrander JSON output. "
            "Statistics will be unavailable."
        )
        return report

    try:
        data = json.loads(json_str)
    except json.JSONDecodeError:
        logger.warning(f"Invalid JSON from Restrander: {json_str[:200]}")
        return report

    # Restrander output keys (may vary by version)
    # Common keys: "forward", "reverse", "unknown", "artefact"/"artefacts"
    report.forward = data.get("forward", data.get("Forward", 0))
    report.reverse = data.get("reverse", data.get("Reverse", 0))
    report.unknown = data.get("unknown", data.get("Unknown", 0))
    report.artefacts = (
        data.get("artefacts", 0)
        or data.get("artefact", 0)
        or data.get("Artefacts", 0)
        or data.get("Artefact", 0)
    )
    report.total_input = (
        report.forward + report.reverse
        + report.unknown + report.artefacts
    )
    report.output_reads = report.forward + report.reverse

    return report


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class RestranderRunner:
    """Restrander wrapper for strand orientation correction.

    Calls the ``restrander`` binary with the bundled PCB109 config.
    All output reads are oriented 5'→3' with primers trimmed.

    Usage::

        runner = RestranderRunner()
        report = runner.run(
            input_fastq=Path("no_foldback.fastq"),
            output_fastq=Path("restranded.fastq"),
        )

    Args:
        config_json: Path to Restrander config JSON.
                     Defaults to bundled PCB109.json.
    """

    def __init__(
        self,
        config_json: Optional[Path] = None,
    ) -> None:
        self.config_json = Path(config_json) if config_json else DEFAULT_CONFIG

        if not self.config_json.exists():
            raise FileNotFoundError(
                f"Restrander config not found: {self.config_json}"
            )

    def run(
        self,
        input_fastq: Path,
        output_fastq: Path,
    ) -> RestranderReport:
        """Run Restrander on a FASTQ file.

        Args:
            input_fastq:  Input FASTQ (foldback-free from Breakinator).
            output_fastq: Output FASTQ with reads oriented 5'→3'.

        Returns:
            RestranderReport with classification statistics.
        """
        restrander_bin = check_binary("restrander")

        # Ensure output directory exists
        output_fastq = Path(output_fastq)
        output_fastq.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            restrander_bin,
            str(input_fastq),
            str(output_fastq),
            str(self.config_json),
        ]

        logger.info(f"Running Restrander: {cmd[0]} ...")
        logger.debug(f"  Config: {self.config_json}")

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Restrander mixes stdout/stderr
            text=True,
        )

        if proc.returncode != 0:
            raise RuntimeError(
                f"Restrander failed (exit {proc.returncode}):\n{proc.stdout}"
            )

        # Parse statistics
        report = _parse_restrander_output(proc.stdout)

        logger.info(f"Restrander stage complete.\n{report}")
        return report
