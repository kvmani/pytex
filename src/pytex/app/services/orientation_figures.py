"""Figures of an orientation relationship determined from measured grains.

Purpose
-------
A relationship fitted to measured parent-child pairs is reported as a
rotation and a name. Two pictures decide whether to believe it: how far each
measured pair sits from the fitted relationship (one pair far off is a
misassigned parent or a different relationship), and how far the fit sits from
every catalogued relationship (a name is earned only by a clear margin). Both
are read from the operation's own payload.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.figures import COLORS, render_figure
from pytex.app.results import ResultFigure

__all__ = ["catalog_distance_figure", "pair_residual_figure"]


def catalog_distance_figure(
    catalog: Sequence[Mapping[str, Any]], *, tolerance_deg: float
) -> ResultFigure:
    """Distance from the fitted relationship to every catalogued one, nearest first."""

    rows = sorted(catalog, key=lambda row: float(row["deviation_deg"]))
    names = [str(row["relationship"]) for row in rows]
    distances = np.array([float(row["deviation_deg"]) for row in rows])

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(names))
        colors = [
            COLORS["accent"] if value <= tolerance_deg else COLORS["guide"] for value in distances
        ]
        axes.barh(positions, np.maximum(distances, 1e-3), color=colors, height=0.6)
        axes.axvline(
            tolerance_deg,
            color=COLORS["warning"],
            lw=0.9,
            ls="--",
            label=f"Naming tolerance {tolerance_deg:g}°",
        )
        for y, value in zip(positions, distances, strict=True):
            axes.text(max(value, 1e-3) * 1.15, y, f"{value:.2f}°", va="center", fontsize=7)
        axes.set_yticks(positions, labels=names, fontsize=7)
        axes.invert_yaxis()
        axes.set_xscale("log")
        axes.set_xlabel("Misorientation from the fitted relationship (°, log scale)")
        axes.legend(loc="lower right", frameon=False, fontsize=7)

    margin = float(distances[1] - distances[0]) if distances.size > 1 else float("nan")
    return render_figure(
        draw,
        key="or_catalog_distances",
        title="How close the fit is to each named relationship",
        caption=(
            "The misorientation between the fitted relationship and every catalogued "
            "relationship for these crystal systems, nearest first (logarithmic scale); green "
            "bars are within the naming tolerance (dashed)."
        ),
        interpretation=(
            f"The nearest relationship leads the next by {margin:.2f}°. A name is only "
            "justified when one relationship is inside the tolerance and the runner-up is well "
            "outside it; two within the tolerance means these measurements cannot tell the two "
            "relationships apart."
        ),
        height_in=0.8 + 0.28 * len(names),
    )


def pair_residual_figure(
    pairs: Sequence[Mapping[str, Any]],
    *,
    mean_deg: float,
    weighted_mean_deg: float | None,
    tolerance_deg: float,
) -> ResultFigure:
    """How far each measured pair sits from the fitted relationship, and how much it counted."""

    labels = [str(row["pair"]) for row in pairs]
    residual = np.array([float(row["residual_deg"]) for row in pairs])
    weights = np.array([float(row.get("normalized_weight") or 0.0) for row in pairs])

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(labels))
        shade = weights / float(weights.max()) if weights.size and weights.max() > 0 else weights
        colors = [
            COLORS["guide"]
            if weight == 0
            else (COLORS["warning"] if value > tolerance_deg else COLORS["data"])
            for value, weight in zip(residual, shade, strict=True)
        ]
        axes.bar(positions, residual, color=colors, width=0.7)
        axes.axhline(
            mean_deg, color=COLORS["model"], lw=1.0, ls="--", label=f"Mean {mean_deg:.3f}°"
        )
        if weighted_mean_deg is not None:
            axes.axhline(
                weighted_mean_deg,
                color=COLORS["difference"],
                lw=1.0,
                ls=":",
                label=f"Weighted mean {weighted_mean_deg:.3f}°",
            )
        axes.axhline(
            tolerance_deg,
            color=COLORS["warning"],
            lw=0.8,
            ls="-.",
            label=f"Naming tolerance {tolerance_deg:g}°",
        )
        axes.set_xticks(positions, labels=labels, fontsize=7)
        axes.set_xlabel("Measured pair")
        axes.set_ylabel("Residual to the fitted relationship (°)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    excluded = int(np.count_nonzero(weights == 0))
    return render_figure(
        draw,
        key="or_pair_residuals",
        title="Each measured pair against the fitted relationship",
        caption=(
            "For every parent-child pair, the misorientation between its measured relationship "
            "and the fitted one (through the nearest variant). Grey pairs carried no weight in "
            f"the fit ({excluded} excluded); red pairs exceed the naming tolerance."
        ),
        interpretation=(
            "Residuals of similar size are the orientation-measurement noise. One pair far "
            "above the rest is a misassigned parent grain or a boundary of a different kind, "
            "and it drags an unweighted fit towards itself; exclude it or give it weight zero "
            "and refit."
        ),
        height_in=3.0,
    )
