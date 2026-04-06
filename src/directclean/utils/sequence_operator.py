"""
Sequence manipulation utilities for DirectClean.

Provides pure functions for DNA sequence operations, including
robust homopolymer detection with dual-criteria scoring
(A/T density + longest consecutive run) designed to catch
imperfect homopolymer regions common in Nanopore sequencing.
"""

from __future__ import annotations

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HomopolymerHit:
    """Result of a homopolymer scan on a single sequence window."""

    is_hit: bool
    density: float  # A/T fraction in the selected window (0.0 - 1.0)
    longest_run_length: int  # longest consecutive A or T stretch
    longest_run_base: str  # which base formed the longest run ('A' or 'T')
    window_seq: str  # the selected window sequence
    window_start: int  # 0-based start on the scanned sequence
    window_end: int  # 0-based exclusive end on the scanned sequence


@dataclass(frozen=True)
class HomopolymerRun:
    """A single homopolymer run found in a sequence."""

    base: str
    start: int
    end: int  # exclusive
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


def longest_consecutive_run(sequence: str, bases: str = "AT") -> tuple[int, str]:
    """
    Find the longest consecutive run of any base in *bases*.

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
        Float in [0.0, 1.0]. Returns 0.0 for empty strings.

    Examples:
        >>> at_density("AAAATAAA")
        1.0
        >>> at_density("ATCGATCG")
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

    Uses dual criteria so that imperfect homopolymers with
    1-2 Nanopore sequencing errors are still caught, while
    true dinucleotide repeats like ATATATATAT are rejected.

    Among all passing windows, returns the most relevant one:
    - highest density first
    - then longer consecutive run
    - then earlier window position

    Args:
        sequence:           DNA string to scan.
        window_size:        Sliding window width in bp.
        density_threshold:  Minimum A/T fraction to call a hit.
        min_run:            Minimum consecutive A or T run.

    Returns:
        HomopolymerHit with exact window coordinates.
    """
    seq = sequence.upper()
    n = len(seq)

    effective_window = min(window_size, n)
    if effective_window == 0:
        return HomopolymerHit(
            is_hit=False,
            density=0.0,
            longest_run_length=0,
            longest_run_base="",
            window_seq="",
            window_start=0,
            window_end=0,
        )

    best_hit: HomopolymerHit | None = None
    best_key = None  # (density, run_len, -start)

    for start in range(n - effective_window + 1):
        end = start + effective_window
        win = seq[start:end]

        d = at_density(win)
        if d < density_threshold:
            continue

        run_len, run_base = longest_consecutive_run(win)
        if run_len < min_run:
            continue

        key = (d, run_len, -start)
        if best_key is None or key > best_key:
            best_key = key
            best_hit = HomopolymerHit(
                is_hit=True,
                density=d,
                longest_run_length=run_len,
                longest_run_base=run_base,
                window_seq=win,
                window_start=start,
                window_end=end,
            )

    if best_hit is not None:
        return best_hit

    overall_run, overall_base = longest_consecutive_run(seq)
    return HomopolymerHit(
        is_hit=False,
        density=at_density(seq),
        longest_run_length=overall_run,
        longest_run_base=overall_base,
        window_seq=seq[:effective_window],
        window_start=0,
        window_end=effective_window,
    )


# ---------------------------------------------------------------------------
# Context extraction around a junction
# ---------------------------------------------------------------------------


def extract_junction_context(
    sequence: str,
    position: int,
    window_size: int = 30,
) -> tuple[str, str]:
    """
    Extract upstream and downstream sequence around a junction position.

    Args:
        sequence:    Full read sequence.
        position:    Junction position (0-indexed, between two bases).
        window_size: Bases to extract on each side.

    Returns:
        (upstream_seq, downstream_seq). May be shorter at boundaries.

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
) -> list[HomopolymerRun]:
    """
    Enumerate every homopolymer run of length ≥ *min_length*.

    Useful for detailed artifact reports and debugging.

    Args:
        sequence:   DNA string.
        min_length: Minimum run length to report.
        bases:      Which bases to scan.

    Returns:
        Sorted list of HomopolymerRun objects.

    Examples:
        >>> find_all_homopolymers("ATCGAAAAACCTTTTT")
        [HomopolymerRun(base='A', start=4, end=9, length=5),
         HomopolymerRun(base='T', start=11, end=16, length=5)]
    """
    seq = sequence.upper()
    target = set(bases.upper())
    runs: list[HomopolymerRun] = []

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
