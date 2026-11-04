"""Interactive Pi sliding-window dot plot (plotly)."""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go

from .sliding_window import Window


def plot_pi(
    windows: list[Window],
    threshold: float,
    out_path: str | Path,
    title: str = "Sliding-window nucleotide diversity (Pi)",
) -> None:
    pts = [w for w in windows if w.pi == w.pi]
    fig = go.Figure(
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
        )
    )
    fig.add_hline(
        y=threshold, line=dict(color="red", dash="dash", width=1),
        annotation_text=f"Pi = {threshold:g}", annotation_position="top right",
    )
    fig.update_layout(
        title=title, xaxis_title="Alignment position (bp)", yaxis_title="Pi",
        template="plotly_white", height=520,
    )
    fig.write_html(str(out_path), include_plotlyjs=True)
