"""
Statistics utilities for DirectClean reporting.

Calculate and format filtering statistics.
"""

from typing import Dict


def calculate_filtering_stats(
    total_reads: int,
    homopolymer_removed: int,
    adapter_removed: int
) -> Dict[str, any]:
    """
    Calculate filtering statistics.
    
    Args:
        total_reads: Total input reads
        homopolymer_removed: Reads removed by homopolymer filter
        adapter_removed: Reads removed by adapter filter
    
    Returns:
        Dictionary with counts and percentages
    
    Examples:
        >>> stats = calculate_filtering_stats(
        ...     total_reads=10000,
        ...     homopolymer_removed=500,
        ...     adapter_removed=300
        ... )
        >>> print(f"Clean: {stats['clean_reads']} ({stats['clean_pct']:.1f}%)")
    """
    # Note: some reads may be removed by both filters
    total_removed = homopolymer_removed + adapter_removed
    clean_reads = total_reads - total_removed
    
    return {
        'total_reads': total_reads,
        'homopolymer_removed': homopolymer_removed,
        'homopolymer_pct': (homopolymer_removed / total_reads * 100) if total_reads > 0 else 0,
        'adapter_removed': adapter_removed,
        'adapter_pct': (adapter_removed / total_reads * 100) if total_reads > 0 else 0,
        'total_removed': total_removed,
        'removed_pct': (total_removed / total_reads * 100) if total_reads > 0 else 0,
        'clean_reads': clean_reads,
        'clean_pct': (clean_reads / total_reads * 100) if total_reads > 0 else 0,
    }


def format_stats_table(stats: Dict) -> str:
    """
    Format statistics as a text table.
    
    Args:
        stats: Statistics dictionary from calculate_filtering_stats
    
    Returns:
        Formatted string table
    
    Examples:
        >>> stats = calculate_filtering_stats(10000, 500, 300)
        >>> print(format_stats_table(stats))
    """
    lines = [
        "DirectClean Filtering Statistics",
        "=" * 50,
        f"Total input reads:        {stats['total_reads']:>10,}",
        f"Homopolymer artifacts:    {stats['homopolymer_removed']:>10,} ({stats['homopolymer_pct']:>5.1f}%)",
        f"Internal adapter reads:   {stats['adapter_removed']:>10,} ({stats['adapter_pct']:>5.1f}%)",
        "-" * 50,
        f"Total removed:            {stats['total_removed']:>10,} ({stats['removed_pct']:>5.1f}%)",
        f"Clean reads:              {stats['clean_reads']:>10,} ({stats['clean_pct']:>5.1f}%)",
        "=" * 50,
    ]
    
    return "\n".join(lines)
