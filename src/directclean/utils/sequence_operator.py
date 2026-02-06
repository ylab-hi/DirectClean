"""
Sequence manipulation utilities for DirectClean.

This module provides basic DNA sequence operations for RT artifact detection.
"""

from typing import Tuple, Optional


def reverse_complement(sequence: str) -> str:
    """
    Generate reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence string (e.g., "ATCG")
    
    Returns:
        Reverse complement of the input sequence (e.g., "CGAT")
    
    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("AAAAA")
        'TTTTT'
    
    Notes:
        - Handles standard DNA bases: A, T, C, G
        - Non-standard bases (N) are preserved
        - Returns uppercase sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'a': 't', 't': 'a',
        'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }
    
    rc = ''.join(complement_map.get(base, base) for base in reversed(sequence))
    return rc.upper()


def has_homopolymer(
    sequence: str,
    min_length: int = 5,
    bases: str = 'AT'
) -> Tuple[bool, Optional[str]]:
    """
    Detect homopolymer runs in a DNA sequence.
    
    Homopolymers (e.g., AAAAA or TTTTT) are prone to RT template
    switching artifacts in Direct-cDNA sequencing.
    
    Args:
        sequence: DNA sequence string to check
        min_length: Minimum consecutive bases (default: 5)
        bases: Bases to check (default: 'AT')
    
    Returns:
        Tuple of (has_homopolymer, base_type):
            - has_homopolymer (bool): True if found
            - base_type (str or None): Which base ('A' or 'T')
    
    Examples:
        >>> has_homopolymer("ATCGAAAAA")
        (True, 'A')
        >>> has_homopolymer("ATCGTTTTT")
        (True, 'T')
        >>> has_homopolymer("ATCGATCG")
        (False, None)
    """
    sequence = sequence.upper()
    bases = bases.upper()
    
    for base in bases:
        pattern = base * min_length
        if pattern in sequence:
            return True, base
    
    return False, None


def extract_junction_context(
    sequence: str,
    position: int,
    window_size: int = 30
) -> Tuple[str, str]:
    """
    Extract upstream and downstream sequence around a junction position.
    
    Used to check for homopolymers near alignment breakpoints.
    
    Args:
        sequence: Full read sequence
        position: Junction position (0-indexed)
        window_size: Bases to extract on each side (default: 30)
    
    Returns:
        Tuple of (upstream_seq, downstream_seq)
    
    Examples:
        >>> seq = "ATCGATCGATCGATCG"
        >>> extract_junction_context(seq, 8, window_size=4)
        ('ATCG', 'ATCG')
    
    Notes:
        - Handles edge cases at sequence boundaries
        - Returns shorter sequences if window exceeds bounds
    """
    upstream_start = max(0, position - window_size)
    downstream_end = min(len(sequence), position + window_size)
    
    upstream = sequence[upstream_start:position]
    downstream = sequence[position:downstream_end]
    
    return upstream, downstream


def find_all_homopolymers(
    sequence: str,
    min_length: int = 5,
    bases: str = 'AT'
) -> list:
    """
    Find all homopolymer runs with their positions.
    
    Useful for detailed artifact analysis and reporting.
    
    Args:
        sequence: DNA sequence string
        min_length: Minimum length to report (default: 5)
        bases: Bases to check (default: 'AT')
    
    Returns:
        List of tuples: (base, start, end, length)
    
    Examples:
        >>> find_all_homopolymers("ATCGAAAAACCTTTTT")
        [('A', 4, 9, 5), ('T', 11, 16, 5)]
    """
    sequence = sequence.upper()
    bases = bases.upper()
    homopolymers = []
    
    for base in bases:
        i = 0
        while i < len(sequence):
            if sequence[i] == base:
                start = i
                while i < len(sequence) and sequence[i] == base:
                    i += 1
                length = i - start
                
                if length >= min_length:
                    homopolymers.append((base, start, i, length))
            else:
                i += 1
    
    return sorted(homopolymers, key=lambda x: x[1])
