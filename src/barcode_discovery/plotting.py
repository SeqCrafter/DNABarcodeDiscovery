"""Interactive Pi sliding-window dot plot (plotly)."""

from __future__ import annotations

from pathlib import Path

import json

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .annotation import Gene
from .sliding_window import Window

_COLOR = {"+": "#2b6cb0", "-": "#dd6b20"}  # by strand
_BLOCK = 0.08   # half height of an exon block (y units)
_PITCH = 1.0    # vertical distance between gene lanes

# Gene names below the genes; shown only when a gene is wide enough for its
# name (re-evaluated on every zoom).
_TRACK_JS = '''
var gd = document.getElementById('{plot_id}');
var NAMES = %s, LEN = %s, LAB = %d;
function updateTrack() {
  var ax = gd._fullLayout.xaxis2 || gd._fullLayout.xaxis;
  var span = ax.range[1] - ax.range[0], px = ax._length;
  var t = NAMES.map(function (n, i) { return LEN[i] / span * px > n.length * 6.5 + 8 ? n : ''; });
  Plotly.restyle(gd, {text: [t]}, [LAB]);
}
gd.on('plotly_relayout', updateTrack);
updateTrack();
'''



# Dashed vertical guide that follows the cursor over the plot area and spans
# both panels, so a Pi point can be read against the gene track.
_SPIKE_JS = '''
var gd = document.getElementById('{plot_id}');
gd.style.position = 'relative';
var guide = document.createElement('div');
guide.style.cssText = 'position:absolute;width:0;border-left:1px dashed #555;pointer-events:none;display:none;z-index:10';
var tag = document.createElement('div');
tag.style.cssText = 'position:absolute;top:-16px;left:3px;font:11px sans-serif;color:#333;background:rgba(255,255,255,.85);padding:0 3px;white-space:nowrap';
guide.appendChild(tag);
gd.appendChild(guide);
gd.addEventListener('mousemove', function (ev) {
  var r = gd.getBoundingClientRect(), s = gd._fullLayout._size;
  var x = ev.clientX - r.left, y = ev.clientY - r.top;
  if (x < s.l || x > s.l + s.w || y < s.t || y > s.t + s.h) { guide.style.display = 'none'; return; }
  guide.style.display = 'block';
  guide.style.left = x + 'px'; guide.style.top = s.t + 'px'; guide.style.height = s.h + 'px';
  tag.textContent = Math.round(gd._fullLayout.xaxis.p2d(x - s.l));
});
gd.addEventListener('mouseleave', function () { guide.style.display = 'none'; });
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


def _rect(a: float, b: float, y: float, h: float) -> tuple[list, list]:
    return [a, b, b, a, a, None], [y - h, y - h, y + h, y + h, y - h, None]


def _arrow_rect(a: float, b: float, y: float, h: float, strand: str) -> tuple[list, list]:
    """Rectangle whose end in the strand direction is a pointed arrow head."""
    hl = min((b - a) * 0.5, 120)
    if strand == "-":
        x = [b, a + hl, a, a + hl, b, b, None]
    else:
        x = [a, b - hl, b, b - hl, a, a, None]
    return x, [y + h, y + h, y, y - h, y - h, y + h, None]


def _add_gene_track(fig: go.Figure, genes: list[Gene]) -> int:
    """Gene models: ``|-----|`` span line plus directional exon arrows,
    coloured by strand. Returns the index of the label trace."""
    lanes = _lanes(genes)
    ys = [-k * _PITCH for k in lanes]
    for strand in "+-":
        sel = [(g, y) for g, y in zip(genes, ys) if g.strand == strand]
        if not sel:
            continue
        lx: list = []
        ly: list = []
        bx: list = []
        by: list = []
        for g, y in sel:
            t = _BLOCK * 1.5  # end ticks
            lx += [g.start, g.end, None, g.start, g.start, None, g.end, g.end, None]
            ly += [y, y, None, y - t, y + t, None, y - t, y + t, None]
            for a, b in (g.exons or ((g.start, g.end),)):
                x, yy = _arrow_rect(a, b, y, _BLOCK, strand)
                bx += x
                by += yy
        fig.add_trace(go.Scatter(x=lx, y=ly, mode="lines", hoverinfo="skip",
                                 showlegend=False,
                                 line=dict(color=_COLOR[strand], width=1.5)), row=2, col=1)
        fig.add_trace(go.Scatter(x=bx, y=by, mode="lines", fill="toself",
                                 fillcolor=_COLOR[strand], hoverinfo="skip",
                                 showlegend=False,
                                 line=dict(color=_COLOR[strand], width=1)), row=2, col=1)
    for g, y in zip(genes, ys):  # transparent hover area over the whole gene
        x, yy = _rect(g.start, g.end, y, _BLOCK * 1.6)
        fig.add_trace(go.Scatter(
            x=x, y=yy, mode="lines", fill="toself", fillcolor="rgba(0,0,0,0)",
            line=dict(width=0), hoveron="fills", hoverinfo="text", showlegend=False,
            text=f"{g.name} ({g.strand})<br>alignment columns {g.start}-{g.end}"
                 + (f"<br>{len(g.exons)} exons" if len(g.exons) > 1 else ""),
        ), row=2, col=1)
    lab = len(fig.data)
    fig.add_trace(go.Scatter(
        x=[(g.start + g.end) / 2 for g in genes], y=[y - _BLOCK - 0.04 for y in ys],
        text=[""] * len(genes), textposition="bottom center", mode="text",
        textfont=dict(size=10, color="#333"), hoverinfo="skip", showlegend=False,
    ), row=2, col=1)
    fig.update_yaxes(visible=False, fixedrange=True,
                     range=[min(ys) - 0.55, 0.3], row=2, col=1)
    return lab


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
        row_heights=[0.72, 0.28] if genes else None,
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
    post = _SPIKE_JS
    if genes:
        lab = _add_gene_track(fig, genes)
        post += _TRACK_JS % (json.dumps([g.name for g in genes]),
                             json.dumps([g.end - g.start + 1 for g in genes]), lab)
    fig.update_xaxes(title_text="Alignment position (bp)", row=2 if genes else 1, col=1)
    fig.update_yaxes(title_text="Pi", row=1, col=1)
    fig.update_layout(title=title, template="plotly_white",
                      height=760 if genes else 520, hovermode="closest")
    fig.write_html(str(out_path), include_plotlyjs=True, post_script=post)
