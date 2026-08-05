"""Read chopper — resolve adapter-associated junctions in Direct-cDNA reads.

AdapterFinder supplies candidate chop sites.  Stage 3 keeps the existing
junction detection and chop positions, but it does not automatically treat
every length-qualified fragment as an independent cDNA molecule.

The final fragment is trimmed as an unsupported terminal residual only when:

* it extends from the final qualified junction to the read end;
* that junction has no downstream TSO signal; and
* the final fragment is shorter than the immediately preceding fragment.

All other fragments continue to use the configured minimum-length rule.  A
single shared segment planner is used by serial and bounded multi-process
execution so both modes preserve identical read-level decisions and ordering.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from directclean.rescuer.adaptor_seq import AdapterConfig
from directclean.rescuer.adapter_finder import (
    AdapterFinder,
    FinderResult,
    InternalJunction,
)
from directclean.utils.io import read_fastq
from directclean.utils.parallel import (
    SerializedFastqRecord,
    bounded_ordered_process_map,
    deserialize_fastq_record,
    iter_serialized_fastq_chunks,
    temporary_output_path,
)

logger = logging.getLogger(__name__)

_DEFAULT_CHUNK_SIZE = 1000
_DEFAULT_CHUNK_BASES = 5_000_000

_KEEP = "kept"
_TOO_SHORT = "too_short"
_TERMINAL_RESIDUAL = "terminal_residual"
_SPLIT_TSO_SUPPORTED = "split_tso_supported"
_SPLIT_NO_TSO_RETAINED = "split_no_tso_retained_by_length_guard"
_TRIM_TERMINAL_RESIDUAL = "trim_terminal_residual"

_REPORT_HEADER = (
    "read_id\tn_junctions\tconfidences\tchop_positions\tdetails"
    "\tn_qualified_junctions\tqualified_actions"
    "\tterminal_residual_bases\n"
)


@dataclass(frozen=True)
class SegmentDecision:
    """One segment produced by a set of qualified chop sites."""

    index: int
    start: int
    end: int
    keep: bool
    reason: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class ChopPlan:
    """Deterministic segment-level decisions for one read."""

    junctions: tuple[InternalJunction, ...]
    segments: tuple[SegmentDecision, ...]
    junction_actions: tuple[str, ...]

    @property
    def terminal_residual_bases(self) -> int:
        return sum(
            segment.length
            for segment in self.segments
            if segment.reason == _TERMINAL_RESIDUAL
        )


@dataclass
class RescueReport:
    """Summary statistics for Stage 3 adapter-structure resolution."""

    total_reads: int = 0
    reads_with_internal: int = 0
    reads_without: int = 0
    total_segments: int = 0
    segments_rescued: int = 0
    segments_discarded: int = 0
    terminal_residuals_trimmed: int = 0
    tso_supported_junctions: int = 0
    no_tso_junctions_retained: int = 0
    input_bases: int = 0
    output_bases: int = 0
    segments_discarded_bases: int = 0
    terminal_residual_bases: int = 0
    # Retained for API compatibility. Streaming paths leave it empty.
    details: list[FinderResult] = field(default_factory=list, repr=False)

    @property
    def accounted_bases(self) -> int:
        """Bases explicitly assigned to output or one discard category."""
        return (
            self.output_bases
            + self.segments_discarded_bases
            + self.terminal_residual_bases
        )

    @property
    def base_accounting_delta(self) -> int:
        """Input minus all explicitly accounted output/discarded bases."""
        return self.input_bases - self.accounted_bases

    def __str__(self) -> str:
        pct = (
            f"{self.reads_with_internal / self.total_reads * 100:.1f}%"
            if self.total_reads > 0
            else "N/A"
        )
        base_pct = (
            f"{self.output_bases / self.input_bases * 100:.1f}%"
            if self.input_bases > 0
            else "N/A"
        )
        return (
            "=== DirectClean Adapter Resolution Report ===\n"
            f"  Total reads processed        : {self.total_reads:,}\n"
            f"  Reads with adapter structure : "
            f"{self.reads_with_internal:,} ({pct})\n"
            f"  Reads passed unchanged       : {self.reads_without:,}\n"
            f"  TSO-supported junctions      : "
            f"{self.tso_supported_junctions:,}\n"
            f"  No-TSO junctions retained: "
            f"{self.no_tso_junctions_retained:,}\n"
            f"  ---\n"
            f"  Segments retained            : {self.segments_rescued:,}\n"
            f"  Short segments discarded     : {self.segments_discarded:,}\n"
            f"  Terminal residuals trimmed   : "
            f"{self.terminal_residuals_trimmed:,}\n"
            f"  Total output records         : {self.total_segments:,}\n"
            f"  ---\n"
            f"  Input bases                  : {self.input_bases:,}\n"
            f"  Output bases                 : "
            f"{self.output_bases:,} ({base_pct})\n"
            f"  Short-segment bases removed  : "
            f"{self.segments_discarded_bases:,}\n"
            f"  Terminal-residual bases trim : "
            f"{self.terminal_residual_bases:,}\n"
            f"  Base accounting delta        : "
            f"{self.base_accounting_delta:,}\n"
            "=============================================="
        )

    def merge(self, other: "RescueReport") -> None:
        """Merge counters from one processed chunk."""
        self.total_reads += other.total_reads
        self.reads_with_internal += other.reads_with_internal
        self.reads_without += other.reads_without
        self.total_segments += other.total_segments
        self.segments_rescued += other.segments_rescued
        self.segments_discarded += other.segments_discarded
        self.terminal_residuals_trimmed += other.terminal_residuals_trimmed
        self.tso_supported_junctions += other.tso_supported_junctions
        self.no_tso_junctions_retained += (
            other.no_tso_junctions_retained
        )
        self.input_bases += other.input_bases
        self.output_bases += other.output_bases
        self.segments_discarded_bases += other.segments_discarded_bases
        self.terminal_residual_bases += other.terminal_residual_bases


def _valid_sorted_junctions(
    read_length: int,
    junctions: list[InternalJunction],
) -> tuple[InternalJunction, ...]:
    """Return strictly increasing in-range junctions."""
    valid: list[InternalJunction] = []
    last_position = -1

    for junction in sorted(junctions, key=lambda item: item.chop_position):
        position = junction.chop_position
        if position <= 0 or position >= read_length:
            continue
        if position <= last_position:
            continue
        valid.append(junction)
        last_position = position

    return tuple(valid)


def _build_chop_plan(
    read_length: int,
    junctions: list[InternalJunction],
    min_segment_length: int,
) -> ChopPlan:
    """Build the one authoritative set of Stage 3 segment decisions.

    The only new biological gate applies to the final read-end segment.
    It is treated as an unsupported terminal residual when its preceding
    junction lacks TSO and the segment is shorter than its immediate upstream
    neighbor.  The length comparison prevents an early no-TSO chop from
    discarding the long downstream body of a read.
    """
    valid_junctions = _valid_sorted_junctions(read_length, junctions)

    if not valid_junctions:
        return ChopPlan(
            junctions=(),
            segments=(
                SegmentDecision(
                    index=1,
                    start=0,
                    end=read_length,
                    keep=True,
                    reason=_KEEP,
                ),
            ),
            junction_actions=(),
        )

    boundaries = [0]
    boundaries.extend(junction.chop_position for junction in valid_junctions)
    boundaries.append(read_length)

    raw_segments = [
        (index + 1, boundaries[index], boundaries[index + 1])
        for index in range(len(boundaries) - 1)
    ]

    terminal_index, terminal_start, terminal_end = raw_segments[-1]
    _, previous_start, previous_end = raw_segments[-2]
    terminal_length = terminal_end - terminal_start
    previous_length = previous_end - previous_start
    last_junction = valid_junctions[-1]

    trim_terminal_residual = (
        last_junction.tso_hit is None
        and terminal_length < previous_length
    )

    decisions: list[SegmentDecision] = []
    for index, start, end in raw_segments:
        is_terminal = index == terminal_index
        length = end - start

        if is_terminal and trim_terminal_residual:
            decisions.append(
                SegmentDecision(
                    index=index,
                    start=start,
                    end=end,
                    keep=False,
                    reason=_TERMINAL_RESIDUAL,
                )
            )
        elif length < min_segment_length:
            decisions.append(
                SegmentDecision(
                    index=index,
                    start=start,
                    end=end,
                    keep=False,
                    reason=_TOO_SHORT,
                )
            )
        else:
            decisions.append(
                SegmentDecision(
                    index=index,
                    start=start,
                    end=end,
                    keep=True,
                    reason=_KEEP,
                )
            )

    actions: list[str] = []
    for index, junction in enumerate(valid_junctions):
        is_last = index == len(valid_junctions) - 1
        if is_last and trim_terminal_residual:
            actions.append(_TRIM_TERMINAL_RESIDUAL)
        elif junction.tso_hit is not None:
            actions.append(_SPLIT_TSO_SUPPORTED)
        else:
            actions.append(_SPLIT_NO_TSO_RETAINED)

    return ChopPlan(
        junctions=valid_junctions,
        segments=tuple(decisions),
        junction_actions=tuple(actions),
    )


def _update_report_from_plan(report: RescueReport, plan: ChopPlan) -> None:
    """Add one chopped read's junction and segment outcomes to a report."""
    report.tso_supported_junctions += sum(
        1 for junction in plan.junctions if junction.tso_hit is not None
    )
    report.no_tso_junctions_retained += sum(
        1
        for action in plan.junction_actions
        if action == _SPLIT_NO_TSO_RETAINED
    )

    for segment in plan.segments:
        if segment.keep:
            report.segments_rescued += 1
            report.total_segments += 1
            report.output_bases += segment.length
        elif segment.reason == _TERMINAL_RESIDUAL:
            report.terminal_residuals_trimmed += 1
            report.terminal_residual_bases += segment.length
        elif segment.reason == _TOO_SHORT:
            report.segments_discarded += 1
            report.segments_discarded_bases += segment.length
        else:
            raise RuntimeError(f"Unknown segment decision: {segment.reason}")


def _chop_record(
    record: SeqRecord,
    junctions: list[InternalJunction],
    min_segment_length: int,
) -> tuple[list[SeqRecord], ChopPlan]:
    """Apply the shared plan and build kept SeqRecord fragments."""
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])
    read_id = record.id
    plan = _build_chop_plan(
        read_length=len(seq),
        junctions=junctions,
        min_segment_length=min_segment_length,
    )

    segments: list[SeqRecord] = []
    for decision in plan.segments:
        if not decision.keep:
            continue

        segment_id = f"{read_id}_part{decision.index}"
        segment = SeqRecord(
            seq=Seq(seq[decision.start:decision.end]),
            id=segment_id,
            name=segment_id,
            description=(
                f"rescued_from={read_id} "
                f"start={decision.start} end={decision.end}"
            ),
        )
        if quals:
            segment.letter_annotations["phred_quality"] = quals[
                decision.start:decision.end
            ]
        segments.append(segment)

    return segments, plan


def _format_finder_result_line(
    result: FinderResult,
    qualified_junctions: tuple[InternalJunction, ...] = (),
    junction_actions: tuple[str, ...] = (),
    terminal_residual_bases: int = 0,
) -> str:
    """Format one per-read TSV row while preserving the original columns."""
    confidences = ",".join(str(j.confidence) for j in result.junctions)
    positions = ",".join(str(j.chop_position) for j in result.junctions)

    detail_parts: list[str] = []
    for junction in result.junctions:
        parts: list[str] = []
        if junction.polya_hit:
            parts.append(
                f"polyA[{junction.polya_hit.start}-{junction.polya_hit.end}]"
            )
        if junction.rtp_rc_hit:
            parts.append(
                f"RTP_rc[{junction.rtp_rc_hit.start}-{junction.rtp_rc_hit.end}]"
                f"(ed={junction.rtp_rc_hit.edit_distance})"
            )
        if junction.tso_hit:
            parts.append(
                f"TSO[{junction.tso_hit.start}-{junction.tso_hit.end}]"
                f"(ed={junction.tso_hit.edit_distance})"
            )
        detail_parts.append("+".join(parts))

    qualified_actions = ";".join(
        f"{junction.chop_position}:{action}"
        for junction, action in zip(qualified_junctions, junction_actions)
    )

    return (
        f"{result.read_id}\t{result.n_chops}\t{confidences}\t"
        f"{positions}\t{';'.join(detail_parts)}\t"
        f"{len(qualified_junctions)}\t{qualified_actions}\t"
        f"{terminal_residual_bases}\n"
    )


_CHOPPER_FINDER: AdapterFinder | None = None
_CHOPPER_MIN_CONFIDENCE = 2
_CHOPPER_MIN_SEGMENT_LENGTH = 50


def _init_chopper_worker(
    config_values: dict,
    min_confidence: int,
    min_segment_length: int,
) -> None:
    """Initialize one long-lived Stage 3 worker process."""
    global _CHOPPER_FINDER
    global _CHOPPER_MIN_CONFIDENCE
    global _CHOPPER_MIN_SEGMENT_LENGTH

    config = AdapterConfig(**config_values)
    _CHOPPER_FINDER = AdapterFinder(config)
    _CHOPPER_MIN_CONFIDENCE = min_confidence
    _CHOPPER_MIN_SEGMENT_LENGTH = min_segment_length


def _process_chopper_chunk(
    chunk: list[SerializedFastqRecord],
) -> tuple[list[SerializedFastqRecord], list[str], RescueReport]:
    """Process one Stage 3 chunk using the shared segment planner."""
    if _CHOPPER_FINDER is None:
        raise RuntimeError("Stage 3 worker was not initialized")

    output_records: list[SerializedFastqRecord] = []
    report_lines: list[str] = []
    report = RescueReport()

    for read_id, name, description, sequence, qualities in chunk:
        report.total_reads += 1
        report.input_bases += len(sequence)
        finder_result = _CHOPPER_FINDER.find(
            read_id=read_id,
            sequence=sequence,
        )
        qualified = [
            junction
            for junction in finder_result.junctions
            if junction.confidence >= _CHOPPER_MIN_CONFIDENCE
        ]

        if not qualified:
            report.reads_without += 1
            report.total_segments += 1
            report.output_bases += len(sequence)
            output_records.append(
                (read_id, name, description, sequence, qualities)
            )
            continue

        plan = _build_chop_plan(
            read_length=len(sequence),
            junctions=qualified,
            min_segment_length=_CHOPPER_MIN_SEGMENT_LENGTH,
        )
        if not plan.junctions:
            report.reads_without += 1
            report.total_segments += 1
            report.output_bases += len(sequence)
            output_records.append(
                (read_id, name, description, sequence, qualities)
            )
            continue

        report.reads_with_internal += 1
        _update_report_from_plan(report, plan)
        report_lines.append(
            _format_finder_result_line(
                finder_result,
                qualified_junctions=plan.junctions,
                junction_actions=plan.junction_actions,
                terminal_residual_bases=plan.terminal_residual_bases,
            )
        )

        for decision in plan.segments:
            if not decision.keep:
                continue

            segment_id = f"{read_id}_part{decision.index}"
            output_records.append(
                (
                    segment_id,
                    segment_id,
                    f"rescued_from={read_id} "
                    f"start={decision.start} end={decision.end}",
                    sequence[decision.start:decision.end],
                    qualities[decision.start:decision.end] if qualities else b"",
                )
            )

    return output_records, report_lines, report


class ReadChopper:
    """Resolve adapter-associated reads and write usable FASTQ records."""

    def __init__(
        self,
        input_fastq: str | Path,
        output_fastq: str | Path,
        config: AdapterConfig | None = None,
        min_confidence: int = 2,
        report_path: str | Path | None = None,
        threads: int = 1,
        chunk_size: int = _DEFAULT_CHUNK_SIZE,
        chunk_bases: int = _DEFAULT_CHUNK_BASES,
    ) -> None:
        self.input_fastq = Path(input_fastq)
        self.output_fastq = Path(output_fastq)
        self.config = config or AdapterConfig()
        self.min_confidence = min_confidence
        self.report_path = Path(report_path) if report_path else None
        self.threads = max(1, threads)
        self.chunk_size = max(1, chunk_size)
        self.chunk_bases = max(1, chunk_bases)

        self.output_fastq.parent.mkdir(parents=True, exist_ok=True)

    def _run_single_thread(self) -> RescueReport:
        """Process and write records immediately in input order."""
        report = RescueReport()
        finder = AdapterFinder(self.config)

        report_handle = None
        if self.report_path is not None:
            self.report_path.parent.mkdir(parents=True, exist_ok=True)
            report_handle = open(self.report_path, "w")
            report_handle.write(_REPORT_HEADER)

        try:
            with open(self.output_fastq, "w") as output_handle:
                for record in read_fastq(self.input_fastq):
                    report.total_reads += 1
                    sequence = str(record.seq)
                    report.input_bases += len(sequence)
                    finder_result = finder.find(
                        read_id=record.id,
                        sequence=sequence,
                    )
                    qualified = [
                        junction
                        for junction in finder_result.junctions
                        if junction.confidence >= self.min_confidence
                    ]

                    if not qualified:
                        report.reads_without += 1
                        report.total_segments += 1
                        report.output_bases += len(sequence)
                        SeqIO.write(record, output_handle, "fastq")
                        continue

                    segments, plan = _chop_record(
                        record,
                        qualified,
                        self.config.min_segment_length,
                    )
                    if not plan.junctions:
                        report.reads_without += 1
                        report.total_segments += 1
                        report.output_bases += len(sequence)
                        SeqIO.write(record, output_handle, "fastq")
                        continue

                    report.reads_with_internal += 1
                    _update_report_from_plan(report, plan)
                    if report_handle is not None:
                        report_handle.write(
                            _format_finder_result_line(
                                finder_result,
                                qualified_junctions=plan.junctions,
                                junction_actions=plan.junction_actions,
                                terminal_residual_bases=(
                                    plan.terminal_residual_bases
                                ),
                            )
                        )

                    for segment in segments:
                        SeqIO.write(segment, output_handle, "fastq")
        finally:
            if report_handle is not None:
                report_handle.close()

        return report

    def _run_parallel(self) -> RescueReport:
        """Run bounded, ordered Stage 3 multiprocessing."""
        workers = self.threads
        max_in_flight = max(workers, workers * 2)
        logger.info(
            "Rescuer: bounded parallel mode with %d workers, "
            "max_in_flight=%d, chunk_reads=%d, chunk_bases=%d",
            workers,
            max_in_flight,
            self.chunk_size,
            self.chunk_bases,
        )

        output_tmp = temporary_output_path(self.output_fastq)
        report_tmp = (
            temporary_output_path(self.report_path)
            if self.report_path is not None
            else None
        )
        output_tmp.unlink(missing_ok=True)
        if report_tmp is not None:
            report_tmp.parent.mkdir(parents=True, exist_ok=True)
            report_tmp.unlink(missing_ok=True)

        merged_report = RescueReport()
        chunks = iter_serialized_fastq_chunks(
            self.input_fastq,
            max_reads=self.chunk_size,
            max_bases=self.chunk_bases,
        )
        config_values = dict(vars(self.config))

        try:
            report_handle = open(report_tmp, "w") if report_tmp else None
            try:
                if report_handle is not None:
                    report_handle.write(_REPORT_HEADER)

                with open(output_tmp, "w") as output_handle:
                    for _, result in bounded_ordered_process_map(
                        _process_chopper_chunk,
                        chunks,
                        max_workers=workers,
                        max_in_flight=max_in_flight,
                        initializer=_init_chopper_worker,
                        initargs=(
                            config_values,
                            self.min_confidence,
                            self.config.min_segment_length,
                        ),
                    ):
                        output_records, report_lines, partial_report = result
                        for serialized in output_records:
                            SeqIO.write(
                                deserialize_fastq_record(serialized),
                                output_handle,
                                "fastq",
                            )
                        if report_handle is not None:
                            report_handle.writelines(report_lines)
                        merged_report.merge(partial_report)
            finally:
                if report_handle is not None:
                    report_handle.close()

            output_tmp.replace(self.output_fastq)
            if report_tmp is not None and self.report_path is not None:
                report_tmp.replace(self.report_path)
        except BaseException:
            output_tmp.unlink(missing_ok=True)
            if report_tmp is not None:
                report_tmp.unlink(missing_ok=True)
            raise

        return merged_report

    def run(self) -> RescueReport:
        """Process the FASTQ and verify exact Stage 3 base accounting."""
        report = (
            self._run_single_thread()
            if self.threads <= 1
            else self._run_parallel()
        )

        if report.base_accounting_delta != 0:
            raise RuntimeError(
                "Stage 3 base accounting failed: "
                f"input={report.input_bases}, output={report.output_bases}, "
                f"short_discarded={report.segments_discarded_bases}, "
                f"terminal_trimmed={report.terminal_residual_bases}, "
                f"delta={report.base_accounting_delta}"
            )

        logger.info("Adapter resolution complete.")
        logger.info("\n%s", report)
        return report

    @staticmethod
    def _format_report_line(result: FinderResult) -> str:
        """Compatibility wrapper for external callers."""
        return _format_finder_result_line(result)
