"""Unknowns rescuer — rescue reads from Restrander's unknowns output.

The parallel implementation is bounded and writes results in original chunk
order, preserving the serial algorithm's FASTQ content and statistics.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import edlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from directclean.rescuer.adaptor_seq import (
    AdapterConfig,
    RTP_SEQUENCE,
    TSO_SEQUENCE,
)
from directclean.rescuer.adapter_finder import AdapterFinder
from directclean.utils.io import read_fastq
from directclean.utils.parallel import (
    SerializedFastqRecord,
    bounded_ordered_process_map,
    deserialize_fastq_record,
    iter_serialized_fastq_chunks,
    temporary_output_path,
)
from directclean.utils.sequence_operator import reverse_complement

logger = logging.getLogger(__name__)

_TSO = TSO_SEQUENCE.upper()
_RTP = RTP_SEQUENCE.upper()
_TSO_RC = reverse_complement(_TSO)

_DEFAULT_CHUNK_SIZE = 1000
_DEFAULT_CHUNK_BASES = 5_000_000


@dataclass
class UnknownsRescueReport:
    """Summary statistics for the unknowns rescue operation."""

    total_unknowns: int = 0
    reads_with_adapter: int = 0
    reads_without: int = 0
    segments_produced: int = 0
    segments_discarded_short: int = 0
    oriented_forward: int = 0
    oriented_reverse: int = 0
    oriented_unknown: int = 0
    output_reads: int = 0

    def __str__(self) -> str:
        pct_adapter = (
            f"{self.reads_with_adapter / self.total_unknowns * 100:.1f}%"
            if self.total_unknowns > 0
            else "N/A"
        )
        return (
            "=== DirectClean Unknowns Rescue Report ===\n"
            f"  Total unknowns scanned  : {self.total_unknowns:,}\n"
            f"  Reads with internal adpt: {self.reads_with_adapter:,} ({pct_adapter})\n"
            f"  Reads without (skipped) : {self.reads_without:,}\n"
            f"  ---\n"
            f"  Sub-reads produced      : {self.segments_produced:,}\n"
            f"  Discarded (too short)   : {self.segments_discarded_short:,}\n"
            f"  Oriented forward        : {self.oriented_forward:,}\n"
            f"  Oriented reverse (RC)   : {self.oriented_reverse:,}\n"
            f"  Orientation unknown     : {self.oriented_unknown:,}\n"
            f"  ---\n"
            f"  Output reads            : {self.output_reads:,}\n"
            "==========================================="
        )

    def merge(self, other: "UnknownsRescueReport") -> None:
        """Merge integer counters from one processed chunk."""
        self.total_unknowns += other.total_unknowns
        self.reads_with_adapter += other.reads_with_adapter
        self.reads_without += other.reads_without
        self.segments_produced += other.segments_produced
        self.segments_discarded_short += other.segments_discarded_short
        self.oriented_forward += other.oriented_forward
        self.oriented_reverse += other.oriented_reverse
        self.oriented_unknown += other.oriented_unknown
        self.output_reads += other.output_reads


def _fuzzy_match(seq: str, query: str, max_ed: int = 3) -> bool:
    """Check if query exists in seq within edit distance."""
    result = edlib.align(query, seq.upper(), mode="HW", k=max_ed)
    return result["editDistance"] != -1


def _has_polya(seq: str, min_run: int = 10) -> bool:
    """Check for a polyA run in the final 200 bases."""
    tail = seq[-200:].upper()
    max_run = 0
    current = 0
    for base in tail:
        if base == "A":
            current += 1
            if current > max_run:
                max_run = current
        else:
            current = 0
    return max_run >= min_run


def _has_polyt(seq: str, min_run: int = 10) -> bool:
    """Check for a polyT run in the first 200 bases."""
    head = seq[:200].upper()
    max_run = 0
    current = 0
    for base in head:
        if base == "T":
            current += 1
            if current > max_run:
                max_run = current
        else:
            current = 0
    return max_run >= min_run


def orient_subread(seq: str) -> str | None:
    """Return the existing serial orientation decision for one segment."""
    head = seq[:200]
    tail = seq[-200:]

    if _fuzzy_match(head, _TSO) or _has_polya(seq):
        return seq
    if _fuzzy_match(head, _RTP) or _has_polyt(seq):
        return reverse_complement(seq)
    if _fuzzy_match(tail, _TSO_RC):
        return reverse_complement(seq)
    return None


_UNKNOWNS_FINDER: AdapterFinder | None = None
_UNKNOWNS_MIN_CONFIDENCE = 2
_UNKNOWNS_MIN_SEGMENT_LENGTH = 50


def _init_unknowns_worker(
    config_values: dict,
    min_confidence: int,
    min_segment_length: int,
) -> None:
    """Initialize one long-lived Unknowns Rescue worker process."""
    global _UNKNOWNS_FINDER
    global _UNKNOWNS_MIN_CONFIDENCE
    global _UNKNOWNS_MIN_SEGMENT_LENGTH

    config = AdapterConfig(**config_values)
    _UNKNOWNS_FINDER = AdapterFinder(config)
    _UNKNOWNS_MIN_CONFIDENCE = min_confidence
    _UNKNOWNS_MIN_SEGMENT_LENGTH = min_segment_length


def _process_unknowns_chunk(
    chunk: list[SerializedFastqRecord],
) -> tuple[list[SerializedFastqRecord], UnknownsRescueReport]:
    """Process one Unknowns Rescue chunk with unchanged read logic."""
    if _UNKNOWNS_FINDER is None:
        raise RuntimeError("Unknowns Rescue worker was not initialized")

    output_records: list[SerializedFastqRecord] = []
    report = UnknownsRescueReport()

    for read_id, _name, _description, sequence, qualities in chunk:
        report.total_unknowns += 1
        finder_result = _UNKNOWNS_FINDER.find(
            read_id=read_id,
            sequence=sequence,
        )
        qualified = [
            junction
            for junction in finder_result.junctions
            if junction.confidence >= _UNKNOWNS_MIN_CONFIDENCE
        ]

        if not qualified:
            report.reads_without += 1
            continue

        report.reads_with_adapter += 1
        boundaries = [0]
        for junction in qualified:
            position = junction.chop_position
            if 0 < position < len(sequence):
                boundaries.append(position)
        boundaries.append(len(sequence))

        for index in range(len(boundaries) - 1):
            start = boundaries[index]
            end = boundaries[index + 1]

            if end - start < _UNKNOWNS_MIN_SEGMENT_LENGTH:
                report.segments_discarded_short += 1
                continue

            report.segments_produced += 1
            segment_sequence = sequence[start:end]
            oriented_sequence = orient_subread(segment_sequence)

            if oriented_sequence is None:
                report.oriented_unknown += 1
                continue

            reverse_oriented = oriented_sequence != segment_sequence
            if reverse_oriented:
                report.oriented_reverse += 1
            else:
                report.oriented_forward += 1

            segment_id = f"{read_id}_rescued{index + 1}"
            segment_qualities = qualities[start:end] if qualities else b""
            if reverse_oriented:
                segment_qualities = segment_qualities[::-1]

            output_records.append(
                (
                    segment_id,
                    segment_id,
                    f"unknowns_rescued_from={read_id} start={start} end={end}",
                    oriented_sequence,
                    segment_qualities,
                )
            )
            report.output_reads += 1

    return output_records, report


class UnknownsRescuer:
    """Rescue oriented sub-reads from Restrander's unknowns FASTQ."""

    def __init__(
        self,
        unknowns_fastq: str | Path,
        output_fastq: str | Path,
        config: AdapterConfig | None = None,
        min_confidence: int = 2,
        min_segment_length: int = 50,
        threads: int = 1,
        chunk_size: int = _DEFAULT_CHUNK_SIZE,
        chunk_bases: int = _DEFAULT_CHUNK_BASES,
    ) -> None:
        self.unknowns_fastq = Path(unknowns_fastq)
        self.output_fastq = Path(output_fastq)
        self.config = config or AdapterConfig()
        self.min_confidence = min_confidence
        self.min_segment_length = min_segment_length
        self.threads = max(1, threads)
        self.chunk_size = max(1, chunk_size)
        self.chunk_bases = max(1, chunk_bases)

        self.output_fastq.parent.mkdir(parents=True, exist_ok=True)

    def _run_single_thread(self) -> UnknownsRescueReport:
        """Execute the unchanged serial algorithm with streaming output."""
        finder = AdapterFinder(self.config)
        report = UnknownsRescueReport()

        with open(self.output_fastq, "w") as output_handle:
            for record in read_fastq(self.unknowns_fastq):
                report.total_unknowns += 1
                sequence = str(record.seq)
                qualities = record.letter_annotations.get("phred_quality", [])

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
                    continue

                report.reads_with_adapter += 1
                boundaries = [0]
                for junction in qualified:
                    position = junction.chop_position
                    if 0 < position < len(sequence):
                        boundaries.append(position)
                boundaries.append(len(sequence))

                for index in range(len(boundaries) - 1):
                    start = boundaries[index]
                    end = boundaries[index + 1]

                    if end - start < self.min_segment_length:
                        report.segments_discarded_short += 1
                        continue

                    report.segments_produced += 1
                    segment_sequence = sequence[start:end]
                    oriented_sequence = orient_subread(segment_sequence)

                    if oriented_sequence is None:
                        report.oriented_unknown += 1
                        continue

                    reverse_oriented = oriented_sequence != segment_sequence
                    if reverse_oriented:
                        report.oriented_reverse += 1
                    else:
                        report.oriented_forward += 1

                    segment_id = f"{record.id}_rescued{index + 1}"
                    segment = SeqRecord(
                        seq=Seq(oriented_sequence),
                        id=segment_id,
                        name=segment_id,
                        description=(
                            f"unknowns_rescued_from={record.id} "
                            f"start={start} end={end}"
                        ),
                    )
                    if qualities:
                        segment_qualities = qualities[start:end]
                        if reverse_oriented:
                            segment_qualities = segment_qualities[::-1]
                        segment.letter_annotations["phred_quality"] = (
                            segment_qualities
                        )

                    SeqIO.write(segment, output_handle, "fastq")
                    report.output_reads += 1

        return report

    def _run_parallel(self) -> UnknownsRescueReport:
        """Run bounded, ordered Unknowns Rescue multiprocessing."""
        workers = self.threads
        max_in_flight = max(workers, workers * 2)
        logger.info(
            "Unknowns Rescue: bounded parallel mode with %d workers, "
            "max_in_flight=%d, chunk_reads=%d, chunk_bases=%d",
            workers,
            max_in_flight,
            self.chunk_size,
            self.chunk_bases,
        )

        output_tmp = temporary_output_path(self.output_fastq)
        output_tmp.unlink(missing_ok=True)
        merged_report = UnknownsRescueReport()
        chunks = iter_serialized_fastq_chunks(
            self.unknowns_fastq,
            max_reads=self.chunk_size,
            max_bases=self.chunk_bases,
        )
        config_values = dict(vars(self.config))

        try:
            with open(output_tmp, "w") as output_handle:
                for _, result in bounded_ordered_process_map(
                    _process_unknowns_chunk,
                    chunks,
                    max_workers=workers,
                    max_in_flight=max_in_flight,
                    initializer=_init_unknowns_worker,
                    initargs=(
                        config_values,
                        self.min_confidence,
                        self.min_segment_length,
                    ),
                ):
                    output_records, partial_report = result
                    for serialized in output_records:
                        SeqIO.write(
                            deserialize_fastq_record(serialized),
                            output_handle,
                            "fastq",
                        )
                    merged_report.merge(partial_report)

            output_tmp.replace(self.output_fastq)
        except BaseException:
            output_tmp.unlink(missing_ok=True)
            raise

        return merged_report

    def run(self) -> UnknownsRescueReport:
        """Execute the unknowns rescue pipeline."""
        logger.info("Unknowns Rescue: scanning for internal adapters ...")
        report = (
            self._run_single_thread()
            if self.threads <= 1
            else self._run_parallel()
        )
        logger.info(
            "Unknowns Rescue complete: %s scanned, %s oriented reads rescued",
            f"{report.total_unknowns:,}",
            f"{report.output_reads:,}",
        )
        logger.info("\n%s", report)
        return report
