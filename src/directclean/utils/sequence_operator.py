"""
Sequence manipulation utilities for DirectClean.

Provides pure functions for DNA sequence operations, including
robust homopolymer detection with dual-criteria scoring
(A/T density + longest consecutive run) designed to catch
imperfect homopolymer regions common in Nanopore sequencing.
"""

from typing import Tuple, Optional, List
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HomopolymerHit:
    """Result of a homopolymer scan on a single sequence window."""
    is_hit: bool
    density: float          # A/T fraction in the window (0.0 - 1.0)
    longest_run_length: int # longest consecutive A or T stretch
    longest_run_base: str   # which base formed the longest run ('A' or 'T')
    window_seq: str         # the actual window sequence examined


@dataclass(frozen=True)
class HomopolymerRun:
    """A single homopolymer run found in a sequence."""
    base: str
    start: int
    end: int      # exclusive
    length: int


# ---------------------------------------------------------------------------
# Core: reverse complement
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Uses str.maketrans for speed over dictionary lookups.

    Args:
        sequence: DNA string (A/T/C/G/N, case-insensitive).

    Returns:
        Upper-case reverse complement.

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("AAAAA")
        'TTTTT'
    """
    return sequence.translate(_COMPLEMENT)[::-1].upper()


# ---------------------------------------------------------------------------
# Core: homopolymer detection — pure functions
# ---------------------------------------------------------------------------

def longest_consecutive_run(sequence: str, bases: str = "AT") -> Tuple[int, str]:
    """
    Find the longest consecutive run of any base in *bases*.

    Scans the sequence once (O(n)) and tracks the longest stretch of
    each target base independently, then returns the overall winner.

    Args:
        sequence: DNA string (will be upper-cased internally).
        bases:    Characters to consider (default "AT").

    Returns:
        (run_length, base) — length of the longest run and which base.
        Returns (0, '') if no target base is found at all.

    Examples:
        >>> longest_consecutive_run("AAAATAAA")
        (4, 'A')
        >>> longest_consecutive_run("TTTTCTTTTA")
        (5, 'T')
        >>> longest_consecutive_run("GCGCGCGC")
        (0, '')
    """
    seq = sequence.upper()
    target = set(bases.upper())

    best_len = 0
    best_base = ""
    cur_len = 0
    cur_base = ""

    for ch in seq:
        if ch in target:
            if ch == cur_base:
                cur_len += 1
            else:
                cur_base = ch
                cur_len = 1
        else:
            if cur_len > best_len:
                best_len = cur_len
                best_base = cur_base
            cur_len = 0
            cur_base = ""

    # handle run that extends to the end of the string
    if cur_len > best_len:
        best_len = cur_len
        best_base = cur_base

    return best_len, best_base


def at_density(sequence: str) -> float:
    """
    Fraction of A + T in a sequence.

    Args:
        sequence: DNA string.

    Returns:
        Float in [0.0, 1.0].  Returns 0.0 for empty strings.

    Examples:
        >>> at_density("AAAATAAA")   # 8 A/T out of 8
        1.0
        >>> at_density("ATCGATCG")   # 4 A/T out of 8
        0.5
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    return (seq.count("A") + seq.count("T")) / len(seq)


def scan_homopolymer(
    sequence: str,
    window_size: int = 10,
    density_threshold: float = 0.8,
    min_run: int = 3,
) -> HomopolymerHit:
    """
    Sliding-window scan for A/T-rich homopolymer regions.

    Uses **dual criteria** so that imperfect homopolymers with
    1-2 Nanopore sequencing errors (e.g. AAAATAAA, TTTTCTTTTA)
    are still caught, while true dinucleotide repeats like
    ATATATATAT are correctly rejected.

    Algorithm
    ---------
    A sliding window of *window_size* bp moves across *sequence*.
    A window is a "hit" when **both** conditions are met:

    1. A/T density  ≥  *density_threshold*
    2. Longest consecutive A or T run  ≥  *min_run*

    The function returns the **worst** (most homopolymer-like) window
    found, i.e. the one with the highest density among all windows
    that satisfy both criteria.

    Args:
        sequence:           DNA string to scan.
        window_size:        Sliding window width in bp (default 10).
        density_threshold:  Minimum A/T fraction to call a hit (default 0.8).
        min_run:            Minimum consecutive A or T run (default 3).

    Returns:
        HomopolymerHit with is_hit=True if any window passes both criteria.

    Examples:
        >>> scan_homopolymer("AAAATAAA").is_hit          # imperfect poly-A
        True
        >>> scan_homopolymer("TTTTCTTTTA").is_hit         # imperfect poly-T
        True
        >>> scan_homopolymer("ATATATATAT").is_hit         # dinuc repeat
        False
        >>> scan_homopolymer("AGTCAGTCAG").is_hit         # normal seq
        False
    """
    seq = sequence.upper()
    n = len(seq)

    # If sequence is shorter than window, treat the whole thing as one window
    effective_window = min(window_size, n)
    if effective_window == 0:
        return HomopolymerHit(False, 0.0, 0, "", "")

    best_hit: Optional[HomopolymerHit] = None
    best_density = -1.0

    for start in range(n - effective_window + 1):
        win = seq[start : start + effective_window]

        d = at_density(win)
        if d < density_threshold:
            continue

        run_len, run_base = longest_consecutive_run(win)
        if run_len < min_run:
            continue

        # Both criteria met — track the most extreme window
        if d > best_density:
            best_density = d
            best_hit = HomopolymerHit(
                is_hit=True,
                density=d,
                longest_run_length=run_len,
                longest_run_base=run_base,
                window_seq=win,
            )

    if best_hit is not None:
        return best_hit

    # No window passed both criteria — return the worst-case stats
    # across the whole sequence so callers can still inspect values.
    overall_run, overall_base = longest_consecutive_run(seq)
    return HomopolymerHit(
        is_hit=False,
        density=at_density(seq),
        longest_run_length=overall_run,
        longest_run_base=overall_base,
        window_seq=seq[:effective_window],
    )


# ---------------------------------------------------------------------------
# Context extraction around a junction
# ---------------------------------------------------------------------------

def extract_junction_context(
    sequence: str,
    position: int,
    window_size: int = 30,
) -> Tuple[str, str]:
    """
    Extract upstream and downstream sequence around a junction position.

    Args:
        sequence:    Full read sequence.
        position:    Junction position (0-indexed, between two bases).
        window_size: Bases to extract on each side (default 30).

    Returns:
        (upstream_seq, downstream_seq).  May be shorter than *window_size*
        at sequence boundaries.

    Examples:
        >>> extract_junction_context("ATCGATCGATCGATCG", 8, window_size=4)
        ('ATCG', 'ATCG')
    """
    up_start = max(0, position - window_size)
    dn_end = min(len(sequence), position + window_size)
    return sequence[up_start:position], sequence[position:dn_end]


# ---------------------------------------------------------------------------
# Utility: enumerate all runs (for reporting / debugging)
# ---------------------------------------------------------------------------

def find_all_homopolymers(
    sequence: str,
    min_length: int = 5,
    bases: str = "AT",
) -> List[HomopolymerRun]:
    """
    Enumerate every homopolymer run of length ≥ *min_length*.

    Useful for detailed artifact reports and debugging.

    Args:
        sequence:   DNA string.
        min_length: Minimum run length to report (default 5).
        bases:      Which bases to scan (default "AT").

    Returns:
        Sorted list of HomopolymerRun objects.

    Examples:
        >>> find_all_homopolymers("ATCGAAAAACCTTTTT")
        [HomopolymerRun(base='A', start=4, end=9, length=5),
        HomopolymerRun(base='T', start=11, end=16, length=5)]
    """
    seq = sequence.upper()
    target = set(bases.upper())
    runs: List[HomopolymerRun] = []

    i = 0
    while i < len(seq):
        ch = seq[i]
        if ch in target:
            start = i
            while i < len(seq) and seq[i] == ch:
                i += 1
            length = i - start
            if length >= min_length:
                runs.append(HomopolymerRun(ch, start, i, length))
        else:
            i += 1

    runs.sort(key=lambda r: r.start)
    return runs