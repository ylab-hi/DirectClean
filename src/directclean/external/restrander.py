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
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

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
        rtp_rtp:          RTP-RTP artefacts (excluded).
        tso_tso:          TSO-TSO artefacts (excluded).
        output_reads:     Reads in the output FASTQ (forward + reverse).
    """

    total_input: int = 0
    forward: int = 0
    reverse: int = 0
    unknown: int = 0
    rtp_rtp: int = 0
    tso_tso: int = 0
    output_reads: int = 0

    @property
    def total_artefacts(self) -> int:
        return self.rtp_rtp + self.tso_tso

    def __str__(self) -> str:
        pct_kept = (
            f"{self.output_reads / self.total_input * 100:.1f}%"
            if self.total_input > 0
            else "N/A"
        )
        pct_artefact = (
            f"{self.total_artefacts / self.total_input * 100:.1f}%"
            if self.total_input > 0
            else "N/A"
        )
        pct_unknown = (
            f"{self.unknown / self.total_input * 100:.1f}%"
            if self.total_input > 0
            else "N/A"
        )
        return (
            "=== Restrander Report ===\n"
            f"  Total input reads       : {self.total_input:,}\n"
            f"  Forward (+)             : {self.forward:,}\n"
            f"  Reverse (-)             : {self.reverse:,}\n"
            f"  Unknown (?)             : {self.unknown:,} ({pct_unknown})\n"
            f"  Artefacts               : {self.total_artefacts:,} ({pct_artefact})\n"
            f"    RTP-RTP               : {self.rtp_rtp:,}\n"
            f"    TSO-TSO               : {self.tso_tso:,}\n"
            f"  ---\n"
            f"  Output reads            : {self.output_reads:,} ({pct_kept})\n"
            "========================="
        )


# ---------------------------------------------------------------------------
# Stats parser
# ---------------------------------------------------------------------------


def _strip_ansi(text: str) -> str:
    """Remove ANSI escape codes from text.

    Restrander outputs coloured log messages mixed with JSON.
    We need to strip these before parsing.
    """
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


def _parse_restrander_output(raw_output: str) -> RestranderReport:
    """Parse Restrander's JSON statistics from its stdout.

    Restrander output format (after stripping ANSI codes)::

        Restrander initialised.
        Input file  : ...
        ...
        {
            "stats": {
                "artefactStats": {
                    "RTP-RTP": 992122,
                    "TSO-TSO": 22696,
                    "no artefact": 3899743
                },
                "strandStats": {
                    "+": 1286044,
                    "-": 2102154,
                    "?": 1526363
                },
                "totalReads": 4914561
            }
        }

    Args:
        raw_output: Combined stdout+stderr from Restrander.

    Returns:
        RestranderReport with parsed statistics.
    """
    report = RestranderReport()

    # Step 1: strip ANSI colour codes
    cleaned = _strip_ansi(raw_output)

    # Step 2: find the LAST JSON object by scanning from the end
    # The JSON stats block is always at the tail of the output
    brace_depth = 0
    json_start = -1
    json_end = -1

    for i in range(len(cleaned) - 1, -1, -1):
        if cleaned[i] == "}":
            if brace_depth == 0:
                json_end = i
            brace_depth += 1
        elif cleaned[i] == "{":
            brace_depth -= 1
            if brace_depth == 0:
                json_start = i
                break

    if json_start < 0 or json_end < 0:
        logger.warning(
            "Could not find JSON in Restrander output. Statistics will be unavailable."
        )
        return report

    json_str = cleaned[json_start : json_end + 1]

    try:
        data = json.loads(json_str)
    except json.JSONDecodeError as e:
        logger.warning(f"Invalid JSON from Restrander: {e}")
        logger.debug(f"JSON string was: {json_str[:500]}")
        return report

    # Step 3: extract stats from the actual Restrander format
    stats = data.get("stats", data)

    # Strand stats: {"+": N, "-": N, "?": N}
    strand_stats = stats.get("strandStats", {})
    report.forward = strand_stats.get("+", 0)
    report.reverse = strand_stats.get("-", 0)
    report.unknown = strand_stats.get("?", 0)

    # Artefact stats: {"RTP-RTP": N, "TSO-TSO": N, "no artefact": N}
    artefact_stats = stats.get("artefactStats", {})
    report.rtp_rtp = artefact_stats.get("RTP-RTP", 0)
    report.tso_tso = artefact_stats.get("TSO-TSO", 0)

    # Total reads
    report.total_input = stats.get("totalReads", 0)

    # Output = forward + reverse (unknowns and artefacts are excluded)
    report.output_reads = report.forward + report.reverse

    logger.debug(f"Parsed Restrander stats: {data}")

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
        config_json: Path | None = None,
    ) -> None:
        self.config_json = Path(config_json) if config_json else DEFAULT_CONFIG

        if not self.config_json.exists():
            raise FileNotFoundError(f"Restrander config not found: {self.config_json}")

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
