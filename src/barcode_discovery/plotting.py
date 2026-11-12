"""Interactive Pi sliding-window dot plot (plotly)."""

from __future__ import annotations

from pathlib import Path

import json

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .annotation import Gene
from .sliding_window import Window

_FILL = {"+": "#9ecae1", "-": "#fdae6b"}
_LINE = {"+": "#3182bd", "-": "#e6550d"}

# Show a gene label only when its drawn width fits the text (re-evaluated on zoom).
_LABEL_JS = '''
var gd = document.getElementById('{plot_id}');
var NAMES = %s, LEN = %s, IDX = %d;
function updateLabels() {
  var ax = gd._fullLayout.xaxis2 || gd._fullLayout.xaxis;
  var span = ax.range[1] - ax.range[0], px = ax._length;
  var t = NAMES.map(function (n, i) { return LEN[i] / span * px > n.length * 6.5 + 8 ? n : ''; });
  Plotly.restyle(gd, {text: [t]}, [IDX]);
}
gd.on('plotly_relayout', updateLabels);
updateLabels();
'''


def _lanes(genes: list[Gene]) -> list[int]:
    """Greedy lane packing so overlapping genes do not collide."""
    order = sorted(range(len(genes)), key=lambda i: (genes[i].start, genes[i].end))
    ends: list[int] = []
    lane = [0] * len(genes)
    for i in order:
        g = genes[i]
        for k, e in enumerate(ends):
            if g.start > e:
                lane[i], ends[k] = k, g.end
                break
        else:
            lane[i] = len(ends)
            ends.append(g.end)
    return lane


def _arrow(g: Gene, y: float) -> tuple[list[float], list[float]]:
    body, head = 0.18, 0.34
    hl = min((g.end - g.start + 1) * 0.4, 150)
    if g.strand == "-":
        x = [g.end, g.start + hl, g.start + hl, g.start, g.start + hl, g.start + hl, g.end, g.end]
    else:
        x = [g.start, g.end - hl, g.end - hl, g.end, g.end - hl, g.end - hl, g.start, g.start]
    return x, [y - body, y - body, y - head, y, y + head, y + body, y + body, y - body]


def _add_gene_track(fig: go.Figure, genes: list[Gene]) -> tuple[list[str], list[int]]:
    lanes = _lanes(genes)
    for g, k in zip(genes, lanes):
        x, y = _arrow(g, -k)
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines", fill="toself", hoveron="fills",
            fillcolor=_FILL[g.strand], line=dict(color=_LINE[g.strand], width=1),
            text=f"{g.name} ({g.strand})<br>alignment columns {g.start}-{g.end}",
            hoverinfo="text", showlegend=False,
        ), row=2, col=1)
    fig.add_trace(go.Scatter(
        x=[(g.start + g.end) / 2 for g in genes], y=[-k for k in lanes],
        mode="text", text=[""] * len(genes), textfont=dict(size=10, color="#222"),
        hoverinfo="skip", showlegend=False,
    ), row=2, col=1)
    fig.update_yaxes(visible=False, fixedrange=True,
                     range=[-max(lanes) - 0.6, 0.6], row=2, col=1)
    return [g.name for g in genes], [g.end - g.start + 1 for g in genes]


def plot_pi(
    windows: list[Window],
    threshold: float,
    out_path: str | Path,
    genes: list[Gene] | None = None,
    title: str = "Sliding-window nucleotide diversity (Pi)",
) -> None:
    pts = [w for w in windows if w.pi == w.pi]
    fig = make_subplots(
        rows=2 if genes else 1, cols=1, shared_xaxes=True, vertical_spacing=0.03,
        row_heights=[0.78, 0.22] if genes else None,
    )
    fig.add_trace(
        go.Scattergl(
            x=[w.midpoint for w in pts],
            y=[w.pi for w in pts],
            mode="markers",
            marker=dict(size=4, color="#1f77b4"),
            customdata=[[w.start, w.end, w.n_sites] for w in pts],
            hovertemplate=(
                "window start: %{customdata[0]}<br>"
                "window end: %{customdata[1]}<br>"
                "midpoint: %{x}<br>"
                "valid sites: %{customdata[2]}<br>"
                "Pi: %{y:.5f}<extra></extra>"
            ),
            name="window Pi",
        ),
        row=1, col=1,
    )
    fig.add_hline(
        y=threshold, line=dict(color="red", dash="dash", width=1),
        annotation_text=f"Pi = {threshold:g}", annotation_position="top right",
        row=1, col=1,
    )
    post = None
    if genes:
        names, lens = _add_gene_track(fig, genes)
        post = _LABEL_JS % (json.dumps(names), json.dumps(lens), len(fig.data) - 1)
    fig.update_xaxes(title_text="Alignment position (bp)", row=2 if genes else 1, col=1)
    fig.update_yaxes(title_text="Pi", row=1, col=1)
    fig.update_layout(title=title, template="plotly_white",
                      height=640 if genes else 520, hovermode="closest")
    fig.write_html(str(out_path), include_plotlyjs=True, post_script=post)
