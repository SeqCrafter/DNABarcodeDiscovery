"""Annotate a reference-genome interval with genes and intergenic spacers.

GFF3 is parsed into memory (a chloroplast annotation is small). Only ``gene``
features are used and named by their ``Name`` attribute (fallbacks: ``gene``,
``ID``).

An interval is split into segments:

* covered by one gene               -> ``geneA``
* covered by several overlapping    -> ``geneA/geneB``
* covered by no gene (spacer)       -> ``{previous gene}-{next gene}``
  The genome is treated as circular, so a spacer before the first gene or
  after the last one is named ``{last gene}-{first gene}``.

Formatted as ``geneA(10bp)+geneA-geneB(180bp)+geneB(20bp)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


@dataclass(frozen=True)
class Gene:
    name: str
    start: int  # 1-based inclusive
    end: int


@dataclass(frozen=True)
class Segment:
    label: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def _attrs(field: str) -> dict[str, str]:
    out = {}
    for item in field.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def load_genes(path: str | Path, seqid: str | None = None) -> list[Gene]:
    """Read ``gene`` features from a GFF3 file (single sequence)."""
    genes: list[Gene] = []
    seqids: set[str] = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            if seqid is not None and cols[0] != seqid:
                continue
            a = _attrs(cols[8])
            name = a.get("Name") or a.get("gene") or a.get("ID")
            if not name:
                raise ValueError(f"{path}: gene at {cols[0]}:{cols[3]}-{cols[4]} has no Name")
            s, e = sorted((int(cols[3]), int(cols[4])))
            genes.append(Gene(name, s, e))
            seqids.add(cols[0])
    if len(seqids) > 1:
        raise ValueError(
            f"{path}: genes on several sequences {sorted(seqids)}; set 'gff3_seqid' in the config"
        )
    if not genes:
        raise ValueError(f"{path}: no 'gene' features found"
                         + (f" for seqid '{seqid}'" if seqid else ""))
    return sorted(genes, key=lambda g: (g.start, g.end))


def annotate(genes: list[Gene], start: int, end: int) -> list[Segment]:
    """Split 1-based inclusive [start, end] into gene / spacer segments."""
    points = {start, end + 1}
    for g in genes:
        for p in (g.start, g.end + 1):
            if start < p <= end:
                points.add(p)
    pts = sorted(points)

    first = min(genes, key=lambda g: (g.start, g.end))
    last = max(genes, key=lambda g: (g.end, g.start))

    segs: list[Segment] = []
    for a, b in zip(pts, pts[1:]):
        b -= 1  # elementary segment [a, b]
        covering = [g.name for g in genes if g.start <= a and g.end >= b]
        if covering:
            label = "/".join(dict.fromkeys(covering))
        else:
            before = [g for g in genes if g.end < a]
            after = [g for g in genes if g.start > b]
            prev = max(before, key=lambda g: (g.end, g.start)) if before else last
            nxt = min(after, key=lambda g: (g.start, g.end)) if after else first
            label = f"{prev.name}-{nxt.name}"
        if segs and segs[-1].label == label:
            segs[-1] = Segment(label, segs[-1].start, b)
        else:
            segs.append(Segment(label, a, b))
    return segs


def format_segments(segs: list[Segment]) -> str:
    return "+".join(f"{s.label}({s.length}bp)" for s in segs)
