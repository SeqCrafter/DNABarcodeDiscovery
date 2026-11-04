"""Nucleotide diversity (Pi), following DnaSP's polymorphism module.

DnaSP (``analyse_region`` in the reference dnasp.py) computes, per region:

* complete deletion: columns where any sequence has a gap / non-ACGT are dropped;
  ``L_net`` is the number of remaining sites;
* ``k`` = mean number of pairwise differences over all n(n-1)/2 sequence pairs;
* ``Pi = k / L_net``.

Here ``k`` is obtained per column from allele counts, which is exactly equal to
the pairwise-difference average but O(n) instead of O(n^2) per site, so that
sliding windows can be evaluated with cumulative sums.
"""

from __future__ import annotations

import numpy as np

from .fasta import Alignment

_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(aln: Alignment) -> np.ndarray:
    """(n, L) int8 matrix: A,C,G,T -> 0..3, anything else -> -1."""
    lut = np.full(256, -1, dtype=np.int8)
    for ch, code in _CODE.items():
        lut[ord(ch)] = code
    raw = np.frombuffer("".join(aln.seqs).encode("ascii", "replace"), dtype=np.uint8)
    return lut[raw].reshape(aln.n, aln.length)


def site_profile(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-column (valid_mask, mean pairwise difference).

    valid_mask  : True where all sequences have an unambiguous base.
    pair_diff   : fraction-free count of differing pairs / number of pairs,
                  i.e. each column's contribution to k (0 for invalid columns).
    """
    n, _ = matrix.shape
    if n < 2:
        raise ValueError("At least 2 sequences are required")
    valid = (matrix >= 0).all(axis=0)
    n_pairs = n * (n - 1) / 2
    same_pairs = np.zeros(matrix.shape[1], dtype=np.float64)
    for code in range(4):
        c = (matrix == code).sum(axis=0).astype(np.float64)
        same_pairs += c * (c - 1) / 2
    diff = (n_pairs - same_pairs) / n_pairs
    diff[~valid] = 0.0
    return valid, diff


def pi_of_region(matrix: np.ndarray) -> float:
    """Pi of an arbitrary alignment slice (complete deletion), NaN if no site."""
    valid, diff = site_profile(matrix)
    l_net = int(valid.sum())
    return float(diff.sum() / l_net) if l_net else float("nan")
