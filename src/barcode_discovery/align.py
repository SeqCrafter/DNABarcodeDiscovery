"""Step 0 (optional): multiple sequence alignment with pymafft (Rust MAFFT)."""

from __future__ import annotations

from pathlib import Path

import pymafft


def run(
    fasta: str | Path,
    out: str | Path,
    strategy: str = "auto",
    maxiterate: int = 0,
    threads: int = 0,
    adjust_direction: bool = False,
    seq_type: str = "nuc",
) -> None:
    """Align an unaligned multi-FASTA and write an upper-case aligned FASTA."""
    result = pymafft.align_file(
        str(fasta),
        strategy=strategy,
        maxiterate=maxiterate,
        seq_type=seq_type,
        adjust_direction=adjust_direction,
        threads=threads,
        progress=lambda line: print(line, flush=True),
    )
    with open(out, "w") as fh:
        for seq in result:
            fh.write(f">{seq.name}\n{seq.sequence.upper()}\n")
    print(f"Aligned {result.nseq} sequences, width {result.width} -> {out}")
