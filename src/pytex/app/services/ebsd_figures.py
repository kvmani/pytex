"""Figures of the data an EBSD analysis starts from.

Purpose
-------
Every EBSD result — grains, boundaries, textures — is computed from the points
that passed an indexing-quality threshold, and a threshold chosen without
looking at the distribution it cuts is a guess. This figure shows each quality
channel's distribution with the threshold on it, and the grain-size
distribution the segmentation produced, so the reader sees what was kept and
what was discarded before reading anything computed from it.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from pytex.app.figures import COLORS, render_figure
from pytex.app.results import ResultFigure

__all__ = ["scan_quality_figure"]


def scan_quality_figure(
    channels: Mapping[str, tuple[str, np.ndarray]],
    *,
    confidence_threshold: float,
    diameters_um: np.ndarray,
) -> ResultFigure:
    """Histograms of the indexing-quality channels and of the grain sizes."""

    panels = [(key, label, values) for key, (label, values) in channels.items()]
    count = len(panels) + (1 if diameters_um.size else 0)
    columns = min(2, max(1, count))
    rows = int(np.ceil(count / columns))

    def draw(figure: Any) -> None:
        grid = figure.subplots(rows, columns, squeeze=False)
        axes_list = list(grid.flat)
        for axes, (key, label, values) in zip(axes_list, panels, strict=False):
            finite = np.asarray(values, dtype=float)
            finite = finite[np.isfinite(finite)]
            axes.hist(
                finite, bins=50, color=COLORS["band_2"], edgecolor=COLORS["difference"], lw=0.4
            )
            if key == "confidence_index":
                kept = float(np.mean(finite >= confidence_threshold)) if finite.size else 0.0
                axes.axvline(
                    confidence_threshold,
                    color=COLORS["warning"],
                    lw=1.0,
                    ls="--",
                    label=f"Threshold {confidence_threshold:g}: {100 * kept:.1f} % kept",
                )
                axes.legend(loc="best", frameon=False, fontsize=7)
            axes.set_xlabel(label)
            axes.set_ylabel("Points")
        if diameters_um.size:
            axes = axes_list[len(panels)]
            axes.hist(
                diameters_um,
                bins=min(40, max(5, diameters_um.size // 3)),
                color=COLORS["band_3"],
                edgecolor=COLORS["model"],
                lw=0.4,
            )
            axes.axvline(
                float(np.mean(diameters_um)),
                color=COLORS["model"],
                lw=1.0,
                label=f"Mean {float(np.mean(diameters_um)):.3g} µm",
            )
            axes.set_xlabel("Equivalent-circle diameter (µm)")
            axes.set_ylabel("Grains")
            axes.legend(loc="best", frameon=False, fontsize=7)
        for axes in axes_list[count:]:
            axes.set_visible(False)

    return render_figure(
        draw,
        key="scan_quality",
        title="Indexing quality and grain sizes",
        caption=(
            "Distribution of each indexing-quality channel over every point of the scan, with "
            "the confidence-index threshold dashed, and the equivalent-circle diameters of the "
            "grains the segmentation found."
        ),
        interpretation=(
            "A threshold that falls inside the main peak of the confidence index discards "
            "well-indexed points; one below a separate low-CI population removes the "
            "misindexed ones. Grains of one or two points are usually noise, and dominate the "
            "count though not the area."
        ),
        height_in=2.4 * rows + 0.4,
    )
