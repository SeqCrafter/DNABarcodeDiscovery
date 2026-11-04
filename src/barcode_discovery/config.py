"""TOML config for step 2 (cluster evaluation)."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

TEMPLATE = """\
# barcode-discovery 第二步配置文件
# 运行: barcode-discovery evaluate config.toml

# 第一步使用的 MAFFT 比对 FASTA
fasta = "aligned.fasta"

# 第一步输出的窗口 Pi 表
window_pi = "results/window_pi.tsv"

# 用于注释真实坐标的物种 ID (与 FASTA 头 '>' 后的 ID 完全一致, 或其第一个空白前的部分)
reference_id = ""

# (可选) 参考基因组的 GFF3 注释; 必须与 reference_id 对应同一个基因组。
# 填写后 evaluate 会把每个区间注释为 geneA(10bp)+geneA-geneB(180bp)+geneB(20bp)。
# 不需要注释时留空。
gff3 = ""
# (可选) GFF3 含多条序列时, 指定与参考基因组对应的 seqid
gff3_seqid = ""

# 输出目录
out = "results"

# 序列聚类比阈值 (%), 大于该值判定为 PASS
ratio_threshold = 80.0

# ---- 候选条形码 ----
# 每个 [[candidate]] 是一个候选条形码, 可以写任意多个。
# intervals: 一个或多个区间, 每个区间 = [起始窗口的起点坐标, 结束窗口的起点坐标]
#   起点坐标取自 window_pi.tsv 的 start 列 (或交互图悬停显示的 window start)。
#   区间内的窗口会被融合为 [起始窗口 start, 结束窗口 end]。
#   多个区间会拼接为同一个条形码 (重叠/相邻的区间自动合并)。
[[candidate]]
name = "candidate_1"
intervals = [
    # [101, 301],
    # [901, 1001],
]

# [[candidate]]
# name = "candidate_2"
# intervals = [
#     [2001, 2101],
# ]
"""


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    intervals: list[tuple[int, int]]


@dataclass(frozen=True)
class Config:
    fasta: Path
    window_pi: Path
    reference_id: str
    out: Path
    ratio_threshold: float
    gff3: Path | None
    gff3_seqid: str | None
    candidates: list[CandidateSpec]


def write_template(path: str | Path, force: bool = False) -> None:
    p = Path(path)
    if p.exists() and not force:
        raise FileExistsError(f"{p} already exists (use --force to overwrite)")
    p.write_text(TEMPLATE, encoding="utf-8")


def load(path: str | Path) -> Config:
    path = Path(path)
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    base = path.parent

    def need(key: str):
        if key not in raw or raw[key] in ("", None):
            raise ValueError(f"config: '{key}' is required")
        return raw[key]

    def resolve(p: str) -> Path:
        q = Path(p)
        return q if q.is_absolute() else base / q

    specs: list[CandidateSpec] = []
    for i, c in enumerate(raw.get("candidate", []), 1):
        name = c.get("name") or f"candidate_{i}"
        ivs = c.get("intervals", [])
        if not ivs:
            raise ValueError(f"config: candidate '{name}' has no intervals")
        parsed = []
        for iv in ivs:
            if len(iv) != 2 or not all(isinstance(x, int) for x in iv):
                raise ValueError(
                    f"config: candidate '{name}': each interval must be two integers, got {iv}")
            if iv[0] > iv[1]:
                raise ValueError(
                    f"config: candidate '{name}': interval {iv} start > end")
            parsed.append((iv[0], iv[1]))
        specs.append(CandidateSpec(name, parsed))
    if not specs:
        raise ValueError("config: no [[candidate]] defined")
    if len({s.name for s in specs}) != len(specs):
        raise ValueError("config: candidate names must be unique")

    return Config(
        fasta=resolve(need("fasta")),
        window_pi=resolve(need("window_pi")),
        reference_id=need("reference_id"),
        out=resolve(raw.get("out", "results")),
        ratio_threshold=float(raw.get("ratio_threshold", 80.0)),
        gff3=resolve(raw["gff3"]) if raw.get("gff3") else None,
        gff3_seqid=raw.get("gff3_seqid") or None,
        candidates=specs,
    )
