"""
File I/O utilities for DirectClean.

Handles reading/writing FASTQ files with memory-efficient streaming.
"""

from pathlib import Path
from typing import Iterator, List, Union, Set
import gzip
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def read_fastq(
    fastq_file: Union[str, Path],
    compressed: bool = None
) -> Iterator[SeqRecord]:
    """
    Read FASTQ file and yield SeqRecord objects.
    
    Memory-efficient: uses generator to yield records one at a time
    instead of loading entire file into memory.
    
    Args:
        fastq_file: Path to FASTQ file
        compressed: If None, auto-detect from .gz extension
    
    Yields:
        SeqRecord objects with sequence, ID, and quality
    
    Examples:
        >>> for record in read_fastq("reads.fastq"):
        ...     print(record.id, len(record.seq))
        
        >>> # Works with compressed files too
        >>> for record in read_fastq("reads.fastq.gz"):
        ...     process_read(record)
    
    Notes:
        - Handles both plain and gzip-compressed files
        - Generator pattern: memory usage stays constant
        - Suitable for very large files (10+ GB)
    """
    fastq_file = Path(fastq_file)
    
    # Auto-detect compression from extension
    if compressed is None:
        compressed = fastq_file.suffix == '.gz'
    
    # Open with appropriate handler
    if compressed:
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')
    
    try:
        for record in SeqIO.parse(handle, "fastq"):
            yield record
    finally:
        handle.close()


def write_fastq(
    records: Union[List[SeqRecord], Iterator[SeqRecord]],
    output_file: Union[str, Path],
    compressed: bool = False
) -> int:
    """
    Write SeqRecord objects to FASTQ file.
    
    Args:
        records: List or iterator of SeqRecord objects
        output_file: Path to output file
        compressed: Whether to gzip compress output
    
    Returns:
        Number of records written
    
    Examples:
        >>> # Write from list
        >>> records = list(read_fastq("input.fastq"))
        >>> count = write_fastq(records, "output.fastq")
        
        >>> # Write from generator (memory-efficient)
        >>> filtered = (r for r in read_fastq("in.fq") if len(r.seq) > 1000)
        >>> count = write_fastq(filtered, "long_reads.fastq")
    
    Notes:
        - Creates parent directories automatically
        - Accepts both lists and generators
        - Can write gzip-compressed output
    """
    output_file = Path(output_file)
    
    # Create parent directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Determine compression
    if compressed or output_file.suffix == '.gz':
        handle = gzip.open(output_file, 'wt')
    else:
        handle = open(output_file, 'w')
    
    try:
        count = SeqIO.write(records, handle, "fastq")
        logger.info(f"Wrote {count} records to {output_file}")
        return count
    finally:
        handle.close()


def count_reads(fastq_file: Union[str, Path]) -> int:
    """
    Fast read counting without loading sequences into memory.
    
    Args:
        fastq_file: Path to FASTQ file
    
    Returns:
        Number of reads
    
    Examples:
        >>> count = count_reads("reads.fastq")
        >>> print(f"File contains {count:,} reads")
    
    Notes:
        - Very fast: only counts lines
        - Does not validate FASTQ format
        - Works with gzipped files
        - FASTQ format: 4 lines per read
    """
    fastq_file = Path(fastq_file)
    
    # Open with appropriate handler
    if fastq_file.suffix == '.gz':
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')
    
    try:
        # Count lines (FASTQ: 4 lines per record)
        line_count = sum(1 for _ in handle)
        read_count = line_count // 4
        
        logger.info(f"Counted {read_count:,} reads in {fastq_file}")
        return read_count
    
    finally:
        handle.close()


def split_fastq_by_ids(
    input_fastq: Union[str, Path],
    remove_ids: Set[str],
    output_kept: Union[str, Path],
    output_removed: Union[str, Path]
) -> tuple:
    """
    Split FASTQ into clean vs artifact reads using streaming I/O.
    
    Memory-efficient implementation: reads one record at a time,
    immediately writes to appropriate output file. Suitable for
    very large files (10+ GB) that won't fit in memory.
    
    Args:
        input_fastq: Input FASTQ file
        remove_ids: Set of read IDs to mark as artifacts
        output_kept: Output file for clean reads
        output_removed: Output file for artifact reads
    
    Returns:
        Tuple of (kept_count, removed_count)
    
    Examples:
        >>> artifact_ids = {'read_001', 'read_005', 'read_010'}
        >>> kept, removed = split_fastq_by_ids(
        ...     "all_reads.fastq",
        ...     artifact_ids,
        ...     "clean.fastq",
        ...     "artifacts.fastq"
        ... )
        >>> print(f"Kept: {kept:,}, Removed: {removed:,}")
    
    Notes:
        - Uses streaming I/O to minimize memory usage
        - Keeps file handles open during entire operation
        - Creates output directories automatically
        - Handles gzip compression based on file extension
    """
    # Setup paths
    kept_path = Path(output_kept)
    removed_path = Path(output_removed)
    
    # Create output directories
    kept_path.parent.mkdir(parents=True, exist_ok=True)
    removed_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Determine compression for each output
    def open_output(path: Path):
        """Helper to open output with correct compression."""
        if path.suffix == '.gz':
            return gzip.open(path, 'wt')
        else:
            return open(path, 'w')
    
    kept_count = 0
    removed_count = 0
    
    # Stream through input, write directly to outputs
    with open_output(kept_path) as f_kept, \
         open_output(removed_path) as f_removed:
        
        for record in read_fastq(input_fastq):
            if record.id in remove_ids:
                # Write to removed file immediately
                SeqIO.write(record, f_removed, "fastq")
                removed_count += 1
            else:
                # Write to kept file immediately
                SeqIO.write(record, f_kept, "fastq")
                kept_count += 1
    
    logger.info(
        f"Split {input_fastq}: "
        f"{kept_count:,} kept, {removed_count:,} removed"
    )
    
    return kept_count, removed_count


def filter_fastq_by_function(
    input_fastq: Union[str, Path],
    output_fastq: Union[str, Path],
    filter_func,
    keep_filtered: bool = False
) -> int:
    """
    Filter FASTQ using a custom function with streaming I/O.
    
    Generic filtering function that can be used for various
    filtering operations while maintaining memory efficiency.
    
    Args:
        input_fastq: Input FASTQ file
        output_fastq: Output FASTQ file
        filter_func: Function that takes SeqRecord and returns bool
                    True = keep record, False = discard
        keep_filtered: If True, output discarded reads instead
    
    Returns:
        Number of records written to output
    
    Examples:
        >>> # Keep only reads longer than 1000bp
        >>> count = filter_fastq_by_function(
        ...     "input.fastq",
        ...     "long_reads.fastq",
        ...     lambda r: len(r.seq) > 1000
        ... )
        
        >>> # Keep reads without homopolymers
        >>> from directclean.utils.sequence import has_homopolymer
        >>> count = filter_fastq_by_function(
        ...     "input.fastq",
        ...     "clean.fastq",
        ...     lambda r: not has_homopolymer(str(r.seq))[0]
        ... )
    
    Notes:
        - Memory-efficient: processes one record at a time
        - Flexible: works with any filtering function
        - Can be used to collect filtered-out reads
    """
    output_path = Path(output_fastq)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Determine compression
    if output_path.suffix == '.gz':
        handle = gzip.open(output_path, 'wt')
    else:
        handle = open(output_path, 'w')
    
    count = 0
    
    try:
        for record in read_fastq(input_fastq):
            # Apply filter function
            should_keep = filter_func(record)
            
            # Invert logic if we want filtered reads
            if keep_filtered:
                should_keep = not should_keep
            
            if should_keep:
                SeqIO.write(record, handle, "fastq")
                count += 1
        
        logger.info(f"Wrote {count:,} filtered records to {output_path}")
        return count
    
    finally:
        handle.close()


def get_read_length_stats(fastq_file: Union[str, Path]) -> dict:
    """
    Calculate read length statistics.
    
    Useful for quality control reporting.
    
    Args:
        fastq_file: Path to FASTQ file
    
    Returns:
        Dictionary with min, max, mean, median lengths
    
    Examples:
        >>> stats = get_read_length_stats("reads.fastq")
        >>> print(f"Mean: {stats['mean']:.1f} bp")
        >>> print(f"Median: {stats['median']} bp")
    
    Notes:
        - Loads all lengths into memory (only integers, not sequences)
        - For very large files (>100M reads), consider sampling
    """
    lengths = []
    
    for record in read_fastq(fastq_file):
        lengths.append(len(record.seq))
    
    if not lengths:
        return {
            'count': 0,
            'min': 0,
            'max': 0,
            'mean': 0,
            'median': 0
        }
    
    lengths.sort()
    
    return {
        'count': len(lengths),
        'min': min(lengths),
        'max': max(lengths),
        'mean': sum(lengths) / len(lengths),
        'median': lengths[len(lengths) // 2]
    }
