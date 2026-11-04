"""Sliding-window Pi scan, replicating DnaSP's window iteration.

DnaSP loop (dnasp.py ``run_analysis``): the window starts at 0, covers
``[pos, min(pos + size, L))``, and the scan stops after the window that
reaches the alignment end. Gap columns inside a window are removed
(complete deletion) and the window midpoint is the middle *valid* site.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Window:
    start: int      # 1-based, inclusive (alignment coordinates)
    end: int        # 1-based, inclusive
    midpoint: int   # 1-based
    n_sites: int    # valid sites used (L_net)
    pi: float       # NaN when n_sites == 0


def scan(
    valid: np.ndarray,
    diff: np.ndarray,
    size: int,
    step: int,
    gaps_excluded: bool = True,
) -> list[Window]:
    """Slide a window over the alignment.

    gaps_excluded=True  (DnaSP "Gaps in Sliding Window: not considered"):
        size/step are counted in *valid* sites (columns without gap / non-ACGT
        in any sequence), so each window holds exactly ``size`` valid sites and
        spans a variable number of alignment columns. start/end are the
        alignment columns of the first/last valid site of the window.
    gaps_excluded=False ("considered"): size/step are alignment columns.
    """
    if size <= 0 or step <= 0:
        raise ValueError("window size and step must be positive")
    if not gaps_excluded:
        return _scan_columns(valid, diff, size, step)

    valid_idx = np.flatnonzero(valid)
    n_valid = len(valid_idx)
    if n_valid == 0:
        return []
    d = diff[valid_idx]
    cum = np.concatenate(([0.0], np.cumsum(d)))

    windows: list[Window] = []
    pos = 0
    while pos < n_valid:
        end = min(pos + size, n_valid)  # exclusive, in valid-site index
        n_sites = end - pos
        pi = float((cum[end] - cum[pos]) / n_sites)
        mid = int(valid_idx[pos + (n_sites - 1) // 2]) + 1
        windows.append(
            Window(int(valid_idx[pos]) + 1, int(valid_idx[end - 1]) + 1, mid, n_sites, pi)
        )
        if end == n_valid:
            break
        pos = min(pos + step, n_valid - 1)
    return windows


def _scan_columns(valid: np.ndarray, diff: np.ndarray, size: int, step: int) -> list[Window]:
    length = len(valid)
    cum = np.concatenate(([0.0], np.cumsum(diff)))
    cum_valid = np.concatenate(([0], np.cumsum(valid, dtype=np.int64)))
    valid_idx = np.flatnonzero(valid)

    windows: list[Window] = []
    pos = 0
    while pos < length:
        end = min(pos + size, length)
        n_sites = int(cum_valid[end] - cum_valid[pos])
        if n_sites:
            pi = float((cum[end] - cum[pos]) / n_sites)
            first = cum_valid[pos]
            mid = int(valid_idx[first + (n_sites - 1) // 2]) + 1
        else:
            pi = float("nan")
            mid = pos + 1
        windows.append(Window(pos + 1, end, mid, n_sites, pi))
        if end == length:
            break
        pos = min(pos + step, length - 1)
    return windows


def write_tsv(windows: list[Window], path: str | Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["start", "end", "midpoint", "valid_sites", "pi"])
        for x in windows:
            w.writerow([x.start, x.end, x.midpoint, x.n_sites, f"{x.pi:.6f}"])


def read_tsv(path: str | Path) -> list[Window]:
    with open(path, newline="") as fh:
        return [
            Window(int(r["start"]), int(r["end"]), int(r["midpoint"]),
                   int(r["valid_sites"]), float(r["pi"]))
            for r in csv.DictReader(fh, delimiter="\t")
        ]
