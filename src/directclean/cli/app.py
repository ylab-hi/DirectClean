"""
DirectClean command-line interface.

Single-command design — no subcommands::

    directclean -i reads.fastq -r genome.fa -o results/ -t 8
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import typer
from rich.console import Console
from rich.logging import RichHandler

from directclean import __version__
from directclean.pipeline import DirectCleanPipeline, PipelineConfig
from directclean.rescuer.adaptor_seq import AdapterConfig
from directclean.filter.homopolymer import HomopolymerConfig

console = Console()

app = typer.Typer(
    name="directclean",
    help=(
        "DirectClean — Remove RT artifacts from Oxford Nanopore "
        "Direct-cDNA sequencing data.\n\n"
        "Processes raw FASTQ through five stages: "
        "Breakinator → Restrander → Rescuer → Minimap2 → Homopolymer Filter."
    ),
    add_completion=False,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _setup_logging(verbose: bool) -> None:
    """Configure logging with Rich handler."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[
            RichHandler(
                console=console,
                show_path=False,
                rich_tracebacks=True,
            )
        ],
    )


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"DirectClean v{__version__}")
        raise typer.Exit()


# ---------------------------------------------------------------------------
# Main command
# ---------------------------------------------------------------------------


@app.callback(invoke_without_command=True)
def main(
    # ---- Required ----
    input_fastq: Path = typer.Option(
        ...,
        "--input",
        "-i",
        help="Raw input FASTQ file from Oxford Nanopore Direct-cDNA sequencing.",
        exists=True,
        dir_okay=False,
        readable=True,
    ),
    reference: Path = typer.Option(
        ...,
        "--reference",
        "-r",
        help="Reference genome FASTA file.",
        exists=True,
        dir_okay=False,
        readable=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output directory for results.",
    ),
    # ---- General ----
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of threads for minimap2, samtools, and breakinator.",
        min=1,
    ),
    prefix: str = typer.Option(
        "directclean",
        "--prefix",
        "-p",
        help="Filename prefix for output files.",
    ),
    # ---- Breakinator parameters ----
    junc_bed: Path | None = typer.Option(
        None,
        "--junc-bed",
        "-j",
        help=(
            "Junction BED12 file for guided minimap2 alignment "
            "(used by Breakinator stage). Recommended: GENCODE annotation."
        ),
        exists=True,
        dir_okay=False,
        readable=True,
    ),
    # ---- Rescuer parameters ----
    max_edit_distance: int = typer.Option(
        3,
        "--max-edit-dist",
        help="Maximum edit distance for adapter fuzzy matching.",
        min=0,
        max=5,
    ),
    min_confidence: int = typer.Option(
        2,
        "--min-confidence",
        help=(
            "Minimum signals (1-3) to chop: "
            "1=TSO only, 2=two of polyA/RTP/TSO, 3=all three."
        ),
        min=1,
        max=3,
    ),
    min_segment_length: int = typer.Option(
        50,
        "--min-segment-len",
        help="Minimum sub-read length to keep after chopping (bp).",
        min=10,
    ),
    # ---- Homopolymer filter parameters ----
    scan_window: int = typer.Option(
        10,
        "--scan-window",
        help="Sliding window size for A/T density scan (bp).",
        min=5,
    ),
    density_threshold: float = typer.Option(
        0.85,
        "--density-threshold",
        help="Minimum A/T fraction in scanning window to flag.",
        min=0.5,
        max=1.0,
    ),
    min_run: int = typer.Option(
        5,
        "--min-run",
        help="Minimum consecutive A or T to flag.",
        min=2,
    ),
    context_window: int = typer.Option(
        50,
        "--context-window",
        help="Bases to extract on each side of a junction.",
        min=10,
    ),
    # ---- Flags ----
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable debug logging.",
    ),
    html_report: bool = typer.Option(
        False,
        "--html-report",
        help="Generate an interactive HTML summary report with charts.",
    ),
    version: bool | None = typer.Option(
        None,
        "--version",
        "-V",
        help="Show version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """
    Remove RT artifacts from Oxford Nanopore Direct-cDNA sequencing data.

    DirectClean processes raw FASTQ files through five stages:

    \b
    1. BREAKINATOR:  Remove foldback inversion artifacts.
    2. RESTRANDER:   Correct strand orientation, trim primers.
    3. RESCUER:      Detect internal TSO/RTP adapters, chop chimeric reads.
    4. MINIMAP2:     Splice-aware alignment to reference genome.
    5. FILTER:       Remove chimeric junctions caused by RT template
                     switching at poly-A/T homopolymer regions.
    """
    _setup_logging(verbose)

    # Build configuration from CLI options
    adapter_cfg = AdapterConfig(
        max_edit_distance=max_edit_distance,
        min_segment_length=min_segment_length,
    )
    homopolymer_cfg = HomopolymerConfig(
        scan_window=scan_window,
        density_threshold=density_threshold,
        min_run=min_run,
        context_window=context_window,
    )
    pipeline_cfg = PipelineConfig(
        adapter_config=adapter_cfg,
        homopolymer_config=homopolymer_cfg,
        min_confidence=min_confidence,
        context_window=context_window,
        threads=threads,
        junc_bed=junc_bed,
    )

    # Run pipeline
    try:
        pipeline = DirectCleanPipeline(
            input_fastq=input_fastq,
            reference=reference,
            output_dir=output_dir,
            config=pipeline_cfg,
            prefix=prefix,
        )
        report = pipeline.run()

        console.print()
        console.print("[bold green]✓ DirectClean completed successfully.[/]")
        console.print(f"  Cleaned reads: [cyan]{pipeline.cleaned_fastq}[/]")
        console.print(f"  Rescued reads: [cyan]{pipeline.rescued_fastq}[/]")

        if html_report:
            from directclean.report import HtmlReportGenerator

            report_path = output_dir / f"{prefix}.report.html"
            generator = HtmlReportGenerator(
                report=report,
                config=pipeline_cfg,
                input_fastq=input_fastq,
                output_dir=output_dir,
                prefix=prefix,
            )
            generator.write(report_path)
            console.print(f"  HTML report: [cyan]{report_path}[/]")

    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/] {e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]Error:[/] {e}")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1)
