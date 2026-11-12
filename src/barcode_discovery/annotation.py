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
    strand: str = "+"
    exons: tuple[tuple[int, int], ...] = ()  # sub-intervals (exon, else CDS)


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
    gene_ids: list[str | None] = []
    seqids: set[str] = set()
    parent_of: dict[str, str] = {}   # feature ID -> Parent (to climb exon -> RNA -> gene)
    children: list[tuple[str, str, int, int]] = []  # (type, parent, start, end)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or (seqid is not None and cols[0] != seqid):
                continue
            a = _attrs(cols[8])
            s, e = sorted((int(cols[3]), int(cols[4])))
            if cols[2] in ("exon", "CDS") and a.get("Parent"):
                children.append((cols[2], a["Parent"].split(",")[0], s, e))
                continue
            if cols[2] != "gene":
                if a.get("ID") and a.get("Parent"):
                    parent_of[a["ID"]] = a["Parent"].split(",")[0]
                continue
            name = a.get("Name") or a.get("gene") or a.get("ID")
            if not name:
                raise ValueError(f"{path}: gene at {cols[0]}:{cols[3]}-{cols[4]} has no Name")
            genes.append(Gene(name, s, e, cols[6] if cols[6] in ("+", "-") else "+"))
            gene_ids.append(a.get("ID"))
            seqids.add(cols[0])
    index = {gid: i for i, gid in enumerate(gene_ids) if gid}
    parts: dict[int, dict[str, set]] = {}
    for kind, parent, s, e in children:
        for _ in range(3):  # exon -> RNA -> gene
            if parent in index or parent not in parent_of:
                break
            parent = parent_of[parent]
        if parent in index:
            parts.setdefault(index[parent], {"exon": set(), "CDS": set()})[kind].add((s, e))
    for i, p in parts.items():
        chosen = p["exon"] or p["CDS"]
        genes[i] = Gene(genes[i].name, genes[i].start, genes[i].end, genes[i].strand,
                        tuple(sorted(chosen)))
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


def to_alignment_columns(genes: list[Gene], ref_seq: str) -> list[Gene]:
    """Map genes from reference (gap-free) coordinates to 1-based alignment
    columns, using the gap pattern of the reference's aligned sequence."""
    cols = [i + 1 for i, c in enumerate(ref_seq) if c != "-"]
    far = max(g.end for g in genes)
    if far > len(cols):
        raise ValueError(
            f"GFF3 gene ends at {far} but the reference has only {len(cols)} bases: "
            "the GFF3 does not match the reference sequence")
    return [Gene(g.name, cols[g.start - 1], cols[g.end - 1], g.strand,
                 tuple((cols[a - 1], cols[b - 1]) for a, b in g.exons if b <= len(cols)))
            for g in genes]
