"""Bounded, ordered process-pool helpers for DirectClean.

The utilities in this module keep multiprocessing memory bounded while
preserving the exact input order of FASTQ-derived output records.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from pathlib import Path
from typing import TypeVar

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from directclean.utils.io import read_fastq


InputT = TypeVar("InputT")
OutputT = TypeVar("OutputT")

# id, name, description, sequence, phred qualities as compact bytes
SerializedFastqRecord = tuple[str, str, str, str, bytes]


def serialize_fastq_record(record: SeqRecord) -> SerializedFastqRecord:
    """Convert a SeqRecord into a compact, pickle-friendly tuple."""
    qualities = record.letter_annotations.get("phred_quality", [])
    return (
        record.id,
        record.name,
        record.description,
        str(record.seq),
        bytes(qualities),
    )


def deserialize_fastq_record(data: SerializedFastqRecord) -> SeqRecord:
    """Reconstruct a SeqRecord without changing its FASTQ header or data."""
    read_id, name, description, sequence, qualities = data
    record = SeqRecord(
        seq=Seq(sequence),
        id=read_id,
        name=name,
        description=description,
    )
    if qualities:
        record.letter_annotations["phred_quality"] = list(qualities)
    return record


def iter_serialized_fastq_chunks(
    fastq_path: str | Path,
    max_reads: int,
    max_bases: int,
) -> Iterator[list[SerializedFastqRecord]]:
    """Yield bounded FASTQ chunks in original input order.

    A chunk closes when adding another record would exceed either the read
    count or total-base limit.  A single read longer than ``max_bases`` is
    emitted alone.
    """
    if max_reads < 1:
        raise ValueError("max_reads must be at least 1")
    if max_bases < 1:
        raise ValueError("max_bases must be at least 1")

    chunk: list[SerializedFastqRecord] = []
    chunk_bases = 0

    for record in read_fastq(fastq_path):
        serialized = serialize_fastq_record(record)
        record_bases = len(serialized[3])

        if chunk and (
            len(chunk) >= max_reads or chunk_bases + record_bases > max_bases
        ):
            yield chunk
            chunk = []
            chunk_bases = 0

        chunk.append(serialized)
        chunk_bases += record_bases

    if chunk:
        yield chunk


def temporary_output_path(path: str | Path) -> Path:
    """Return a hidden sibling path suitable for atomic output replacement."""
    path = Path(path)
    return path.with_name(f".{path.name}.directclean.tmp")


def bounded_ordered_process_map(
    function: Callable[[InputT], OutputT],
    items: Iterable[InputT],
    *,
    max_workers: int,
    max_in_flight: int,
    initializer: Callable[..., None] | None = None,
    initargs: tuple = (),
) -> Iterator[tuple[int, OutputT]]:
    """Process items concurrently and yield results in submission order.

    At most ``max_in_flight`` submitted-but-not-yielded items are retained.
    Completed results count toward that bound, so a slow early chunk cannot
    cause an unbounded reorder buffer.
    """
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    if max_in_flight < max_workers:
        raise ValueError("max_in_flight must be at least max_workers")

    indexed_items = enumerate(items)
    pending: dict[Future[OutputT], int] = {}
    ready: dict[int, OutputT] = {}
    next_to_yield = 0
    exhausted = False

    def fill(executor: ProcessPoolExecutor) -> None:
        nonlocal exhausted
        while not exhausted and len(pending) + len(ready) < max_in_flight:
            try:
                index, item = next(indexed_items)
            except StopIteration:
                exhausted = True
                return
            future = executor.submit(function, item)
            pending[future] = index

    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=initializer,
        initargs=initargs,
    ) as executor:
        fill(executor)

        try:
            while pending or ready:
                while next_to_yield in ready:
                    result = ready.pop(next_to_yield)
                    yield next_to_yield, result
                    next_to_yield += 1
                    fill(executor)

                if not pending:
                    if ready:
                        raise RuntimeError(
                            "Ordered process map ended with a missing chunk result"
                        )
                    break

                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    index = pending.pop(future)
                    ready[index] = future.result()
        except BaseException:
            for future in pending:
                future.cancel()
            raise
