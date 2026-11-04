"""Step 2: fuse user-selected windows, compute Pi / cluster ratio, rank."""

from __future__ import annotations

import csv

import numpy as np

from . import clustering, diversity, ranking, sliding_window
from .config import CandidateSpec, Config
from .fasta import Alignment, find_sequence, read_alignment


def _merge(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or touching column ranges."""
    out: list[tuple[int, int]] = []
    for a, b in sorted(ranges):
        if out and a <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _fmt(ranges) -> str:
    return ";".join("NA" if r is None else f"{r[0]}-{r[1]}" for r in ranges)


def _ref_range(cum: np.ndarray, a: int, b: int) -> tuple[int, int] | None:
    """Reference (gap-free) coordinates of alignment columns a..b (1-based)."""
    if cum[b] - cum[a - 1] == 0:
        return None
    return int(cum[a - 1]) + 1, int(cum[b])


def _build(spec: CandidateSpec, windows, aln: Alignment, valid, diff, ref_cum,
           ratio_threshold: float) -> ranking.Candidate:
    by_start = {w.start: w for w in windows}
    cols: list[tuple[int, int]] = []
    sel = set()
    for s1, s2 in spec.intervals:
        for s in (s1, s2):
            if s not in by_start:
                raise ValueError(
                    f"candidate '{spec.name}': {s} is not a window start in window_pi.tsv")
        first, last = by_start[s1], by_start[s2]
        cols.append((first.start, last.end))
        sel.update(w.start for w in windows if first.start <= w.start <= last.start)

    merged = _merge(cols)
    mask = np.zeros(aln.length, dtype=bool)
    for a, b in merged:
        mask[a - 1:b] = True
    n_valid = int(valid[mask].sum())
    pi = float(diff[mask].sum() / n_valid) if n_valid else float("nan")
    wpis = [by_start[s].pi for s in sel if by_start[s].pi == by_start[s].pi]

    ref_ranges = [_ref_range(ref_cum, a, b) for a, b in merged]
    length = sum(r[1] - r[0] + 1 for r in ref_ranges if r)
    k, n, ratio = clustering.cluster_ratio(clustering.extract_ranges(aln, merged))
    return ranking.Candidate(
        name=spec.name, aln_ranges=merged, ref_ranges=ref_ranges,
        n_intervals=len(spec.intervals), n_windows=len(sel), length=length,
        valid_sites=n_valid, pi=pi,
        mean_window_pi=float(np.mean(wpis)) if wpis else float("nan"),
        n_clusters=k, n_seqs=n, cluster_ratio=ratio,
        passed=ratio > ratio_threshold,
    )


def run(cfg: Config) -> list[ranking.Candidate]:
    cfg.out.mkdir(parents=True, exist_ok=True)
    aln = read_alignment(cfg.fasta)
    ref_idx = find_sequence(aln, cfg.reference_id)
    ref_seq = aln.seqs[ref_idx]
    ref_cum = np.concatenate(([0], np.cumsum([c != "-" for c in ref_seq])))

    windows = sliding_window.read_tsv(cfg.window_pi)
    valid, diff = diversity.site_profile(diversity.encode(aln))

    cands = ranking.rank([
        _build(s, windows, aln, valid, diff, ref_cum, cfg.ratio_threshold)
        for s in cfg.candidates
    ])

    with open(cfg.out / "barcode_candidates.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "name", "aln_region", f"{cfg.reference_id}_region",
                    "pi", "mean_window_pi", "cluster_ratio_pct", "n_clusters",
                    "n_seqs", "n_intervals", "n_windows", "length_bp",
                    "valid_sites", "pass"])
        for c in cands:
            w.writerow([c.rank, c.name, _fmt(c.aln_ranges), _fmt(c.ref_ranges),
                        f"{c.pi:.6f}", f"{c.mean_window_pi:.6f}",
                        f"{c.cluster_ratio:.2f}", c.n_clusters, c.n_seqs,
                        c.n_intervals, c.n_windows, c.length, c.valid_sites,
                        "PASS" if c.passed else "FAIL"])

    for c in cands:  # aligned (gaps kept) fused region of every sequence
        sub = clustering.extract_ranges(aln, c.aln_ranges)
        with open(cfg.out / f"{c.name}.fasta", "w") as fh:
            for name, s in zip(aln.names, sub):
                fh.write(f">{name}\n{s}\n")
    return cands
