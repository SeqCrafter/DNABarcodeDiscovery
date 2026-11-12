"""Step 1: sliding-window Pi scan -> window_pi.tsv + interactive plot."""

from __future__ import annotations

from pathlib import Path

from . import annotation, diversity, plotting, sliding_window
from .fasta import find_sequence, read_alignment


def run(
    fasta: str | Path,
    out_dir: str | Path,
    window: int = 200,
    step: int = 10,
    pi_threshold: float = 0.05,
    gaps_excluded: bool = True,
    gff3: str | Path | None = None,
    reference_id: str | None = None,
    gff3_seqid: str | None = None,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    aln = read_alignment(fasta)
    genes = None
    if gff3:
        if not reference_id:
            raise ValueError("--gff3 requires --reference-id (the genome the GFF3 describes)")
        ref_seq = aln.seqs[find_sequence(aln, reference_id)]
        genes = annotation.to_alignment_columns(
            annotation.load_genes(gff3, gff3_seqid), ref_seq)
    valid, diff = diversity.site_profile(diversity.encode(aln))
    windows = sliding_window.scan(valid, diff, window, step, gaps_excluded)
    sliding_window.write_tsv(windows, out / "window_pi.tsv")
    plotting.plot_pi(windows, pi_threshold, out / "pi_sliding_window.html", genes)
    print(f"{len(windows)} windows -> {out}/window_pi.tsv, {out}/pi_sliding_window.html")
