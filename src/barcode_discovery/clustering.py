"""Sequence-cluster ratio (CD-HIT ``-c 1`` equivalent: exact de-duplication).

Each region is cut from the alignment with the same coordinates for every
sequence (equal length, gaps kept), so ``c = 1`` means the aligned slices are
100% identical. ratio = number of clusters / number of sequences * 100.
"""

from .fasta import Alignment


def extract_region(aln: Alignment, start: int, end: int) -> list[str]:
    """Aligned sub-sequences (gaps kept) for 1-based inclusive [start, end]."""
    return [s[start - 1:end] for s in aln.seqs]


def extract_ranges(aln: Alignment, ranges: list[tuple[int, int]]) -> list[str]:
    """Concatenate the aligned slices of several 1-based inclusive ranges."""
    return ["".join(s[a - 1:b] for a, b in ranges) for s in aln.seqs]


def cluster_ratio(seqs: list[str]) -> tuple[int, int, float]:
    """Return (n_clusters, n_sequences, ratio_percent)."""
    if not seqs:
        return 0, 0, float("nan")
    n_clusters = len(set(seqs))
    return n_clusters, len(seqs), 100.0 * n_clusters / len(seqs)
