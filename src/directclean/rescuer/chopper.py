"""Read chopper — split reads at internal adapter junctions.

Takes the junction sites detected by AdapterFinder and produces rescued
sub-reads as independent FASTQ records.  Short fragments below a configurable
threshold are discarded.

Single-thread and bounded multi-process modes preserve the same FASTQ and TSV
ordering.
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


@dataclass
class RescueReport:
    """Summary statistics for the rescue operation."""

    total_reads: int = 0
    reads_with_internal: int = 0
    reads_without: int = 0
    total_segments: int = 0
    segments_rescued: int = 0
    segments_discarded: int = 0
    # Retained for API compatibility. Streaming paths leave it empty.
    details: list[FinderResult] = field(default_factory=list, repr=False)

    def __str__(self) -> str:
        pct = (
            f"{self.reads_with_internal / self.total_reads * 100:.1f}%"
            if self.total_reads > 0
            else "N/A"
        )
        return (
            "=== DirectClean Rescue Report ===\n"
            f"  Total reads processed   : {self.total_reads:,}\n"
            f"  Reads with internal adpt: {self.reads_with_internal:,} ({pct})\n"
            f"  Reads passed unchanged  : {self.reads_without:,}\n"
            f"  ---\n"
            f"  Segments rescued        : {self.segments_rescued:,}\n"
            f"  Segments discarded      : {self.segments_discarded:,}\n"
            f"  Total output reads      : {self.total_segments:,}\n"
            "================================="
        )

    def merge(self, other: "RescueReport") -> None:
        """Merge integer counters from one processed chunk."""
        self.total_reads += other.total_reads
        self.reads_with_internal += other.reads_with_internal
        self.reads_without += other.reads_without
        self.total_segments += other.total_segments
        self.segments_rescued += other.segments_rescued
        self.segments_discarded += other.segments_discarded


def _chop_record(
    record: SeqRecord,
    junctions: list[InternalJunction],
    min_segment_length: int,
) -> tuple[list[SeqRecord], int]:
    """Split a SeqRecord at the supplied junction positions."""
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])
    read_id = record.id

    boundaries = [0]
    for junction in junctions:
        boundaries.append(junction.chop_position)
    boundaries.append(len(seq))

    segments: list[SeqRecord] = []
    discarded = 0

    for index in range(len(boundaries) - 1):
        start = boundaries[index]
        end = boundaries[index + 1]

        if end - start < min_segment_length:
            discarded += 1
            continue

        segment_id = f"{read_id}_part{index + 1}"
        segment = SeqRecord(
            seq=Seq(seq[start:end]),
            id=segment_id,
            name=segment_id,
            description=f"rescued_from={read_id} start={start} end={end}",
        )
        if quals:
            segment.letter_annotations["phred_quality"] = quals[start:end]
        segments.append(segment)

    return segments, discarded


def _format_finder_result_line(result: FinderResult) -> str:
    """Format one per-read TSV row exactly as the serial implementation."""
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

    return (
        f"{result.read_id}\t{result.n_chops}\t{confidences}\t"
        f"{positions}\t{';'.join(detail_parts)}\n"
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
    """Process one Stage 3 chunk without changing read-level logic."""
    if _CHOPPER_FINDER is None:
        raise RuntimeError("Stage 3 worker was not initialized")

    output_records: list[SerializedFastqRecord] = []
    report_lines: list[str] = []
    report = RescueReport()

    for read_id, name, description, sequence, qualities in chunk:
        report.total_reads += 1
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
            output_records.append(
                (read_id, name, description, sequence, qualities)
            )
            continue

        report.reads_with_internal += 1
        report_lines.append(_format_finder_result_line(finder_result))

        boundaries = [0]
        for junction in qualified:
            boundaries.append(junction.chop_position)
        boundaries.append(len(sequence))

        for index in range(len(boundaries) - 1):
            start = boundaries[index]
            end = boundaries[index + 1]

            if end - start < _CHOPPER_MIN_SEGMENT_LENGTH:
                report.segments_discarded += 1
                continue

            segment_id = f"{read_id}_part{index + 1}"
            output_records.append(
                (
                    segment_id,
                    segment_id,
                    f"rescued_from={read_id} start={start} end={end}",
                    sequence[start:end],
                    qualities[start:end] if qualities else b"",
                )
            )
            report.segments_rescued += 1
            report.total_segments += 1

    return output_records, report_lines, report


class ReadChopper:
    """Chop reads at internal adapter junctions and write rescued FASTQ."""

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
            report_handle.write(
                "read_id\tn_junctions\tconfidences\tchop_positions"
                "\tdetails\n"
            )

        try:
            with open(self.output_fastq, "w") as output_handle:
                for record in read_fastq(self.input_fastq):
                    report.total_reads += 1
                    sequence = str(record.seq)
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
                        SeqIO.write(record, output_handle, "fastq")
                        continue

                    report.reads_with_internal += 1
                    if report_handle is not None:
                        report_handle.write(
                            _format_finder_result_line(finder_result)
                        )

                    segments, discarded = _chop_record(
                        record,
                        qualified,
                        self.config.min_segment_length,
                    )
                    report.segments_rescued += len(segments)
                    report.segments_discarded += discarded
                    report.total_segments += len(segments)

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
                    report_handle.write(
                        "read_id\tn_junctions\tconfidences\tchop_positions"
                        "\tdetails\n"
                    )

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
        """Process the entire FASTQ and write rescued output."""
        report = (
            self._run_single_thread()
            if self.threads <= 1
            else self._run_parallel()
        )
        logger.info("Rescue complete.")
        logger.info("\n%s", report)
        return report

    @staticmethod
    def _format_report_line(result: FinderResult) -> str:
        """Compatibility wrapper for external callers and tests."""
        return _format_finder_result_line(result)
