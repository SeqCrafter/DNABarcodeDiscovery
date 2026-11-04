"""Combine Pi and cluster ratio into a candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Candidate:
    name: str
    aln_ranges: list[tuple[int, int]]            # merged, alignment columns (1-based)
    ref_ranges: list[tuple[int, int] | None]     # same ranges on the reference species
    n_intervals: int                             # user-specified intervals fused
    n_windows: int
    length: int                                  # reference bases, gaps excluded
    valid_sites: int                             # columns used for Pi
    pi: float
    mean_window_pi: float
    n_clusters: int
    n_seqs: int
    cluster_ratio: float                         # percent
    passed: bool
    rank: int = field(default=0)


def rank(cands: list[Candidate]) -> list[Candidate]:
    """Order candidates: passing the ratio filter first, then by the sum of the
    Pi rank and the cluster-ratio rank (lower is better); ties are broken by
    cluster ratio, then Pi.
    """
    pi_rank = {id(c): i for i, c in enumerate(sorted(cands, key=lambda c: -c.pi), 1)}
    ratio_rank = {id(c): i for i, c in
                  enumerate(sorted(cands, key=lambda c: -c.cluster_ratio), 1)}
    ordered = sorted(
        cands,
        key=lambda c: (not c.passed, pi_rank[id(c)] + ratio_rank[id(c)],
                       -c.cluster_ratio, -c.pi),
    )
    for i, c in enumerate(ordered, 1):
        c.rank = i
    return ordered
