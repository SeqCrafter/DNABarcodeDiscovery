"""Step 1: sliding-window Pi scan -> window_pi.tsv + interactive plot."""

from __future__ import annotations

from pathlib import Path

from . import diversity, plotting, sliding_window
from .fasta import read_alignment


def run(
    fasta: str | Path,
    out_dir: str | Path,
    window: int = 200,
    step: int = 10,
    pi_threshold: float = 0.05,
    gaps_excluded: bool = True,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    aln = read_alignment(fasta)
    valid, diff = diversity.site_profile(diversity.encode(aln))
    windows = sliding_window.scan(valid, diff, window, step, gaps_excluded)
    sliding_window.write_tsv(windows, out / "window_pi.tsv")
    plotting.plot_pi(windows, pi_threshold, out / "pi_sliding_window.html")
    print(f"{len(windows)} windows -> {out}/window_pi.tsv, {out}/pi_sliding_window.html")
