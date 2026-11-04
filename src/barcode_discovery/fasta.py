"""Minimal FASTA reader for aligned sequences."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Alignment:
    names: list[str]
    seqs: list[str]  # upper-case, equal length, '-' for gaps

    @property
    def n(self) -> int:
        return len(self.seqs)

    @property
    def length(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0


def read_alignment(path: str | Path) -> Alignment:
    names: list[str] = []
    chunks: list[list[str]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                names.append(line[1:].strip())
                chunks.append([])
            elif names:
                chunks[-1].append(line)
    if not names:
        raise ValueError(f"No sequences found in {path}")
    seqs = ["".join(c).upper().replace(".", "-") for c in chunks]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(
            f"Sequences are not aligned (found {len(lengths)} different lengths)"
        )
    return Alignment(names=names, seqs=seqs)


def find_sequence(aln: Alignment, seq_id: str) -> int:
    """Index of the sequence whose FASTA header equals ``seq_id`` (full header
    or its first whitespace-delimited token)."""
    hits = [i for i, n in enumerate(aln.names)
            if n == seq_id or n.split()[0] == seq_id]
    if not hits:
        raise ValueError(f"Species ID '{seq_id}' not found in the FASTA headers")
    if len(hits) > 1:
        raise ValueError(f"Species ID '{seq_id}' matches {len(hits)} sequences")
    return hits[0]
