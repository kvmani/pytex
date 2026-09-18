# ruff: noqa: RUF001
"""Figures that test an ODF against the pole figures it was inverted from.

Purpose
-------
Pole-figure inversion is ill-posed, so an ODF is only as credible as its
recalculated pole figures are close to the measured ones. The panel draws the
figures themselves; these plots put the comparison in numbers a reader can
judge at a glance — every measured intensity against its recalculated value,
the misfit against tilt (where defocusing and a truncated measurement show up),
and the component volume fractions against what a random specimen would hold.
They read the operation's own arrays and compute nothing new.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import pairwise
from typing import Any

import numpy as np

from pytex.app.figures import COLORS, render_figure
from pytex.app.results import ResultFigure

__all__ = ["fraction_figure", "parity_figure", "tilt_misfit_figure"]

_PALETTE = ("#2563eb", "#d97706", "#0f766e", "#7c3aed", "#db2777", "#65a30d")


def parity_figure(figures: Sequence[Mapping[str, Any]]) -> ResultFigure:
    """Recalculated against measured intensity, one panel per pole figure."""

    count = len(figures)
    columns = min(3, max(1, count))
    rows = int(np.ceil(count / columns))

    def draw(figure: Any) -> None:
        grid = figure.subplots(rows, columns, squeeze=False)
        for slot, axes in enumerate(grid.flat):
            if slot >= count:
                axes.set_visible(False)
                continue
            entry = figures[slot]
            measured = np.array([point["measured"] for point in entry["points"]], dtype=float)
            recalculated = np.array(
                [point["recalculated"] for point in entry["points"]], dtype=float
            )
            top = float(max(measured.max(initial=1.0), recalculated.max(initial=1.0))) * 1.05
            axes.plot([0.0, top], [0.0, top], color=COLORS["guide"], lw=0.8, ls="--")
            axes.plot(
                measured, recalculated, "o", ms=1.8, alpha=0.6, color=_PALETTE[slot % len(_PALETTE)]
            )
            axes.set_xlim(0.0, top)
            axes.set_ylim(0.0, top)
            axes.set_aspect("equal")
            axes.set_title(str(entry["label"]), fontsize=8)
            rp = entry["residual"]["rp_percent"]
            fit = entry["fit_residual"]["rp_percent"]
            axes.text(
                0.04,
                0.96,
                f"RP {rp:.1f} %\nRP vs inverted {fit:.1f} %",
                transform=axes.transAxes,
                va="top",
                fontsize=6.5,
            )
            axes.tick_params(labelsize=6.5)
        figure.supxlabel("Measured intensity (m.r.d.)", fontsize=8)
        figure.supylabel("Recalculated from the ODF (m.r.d.)", fontsize=8)

    return render_figure(
        draw,
        key="odf_parity",
        title="Recalculated against measured pole figures",
        caption=(
            "Every measured point of each pole figure: its measured intensity against the "
            "intensity the ODF predicts in the same direction, with the 1:1 line. RP is the mean "
            "relative misfit against the measured figure, and against the symmetrized figure "
            "that was inverted."
        ),
        interpretation=(
            "Points along the dashed line are reproduced. A cloud flattened below the line at "
            "high intensity means the ODF is too smooth to reach the peaks (a wide kernel or "
            "strong regularisation); a gap between the two RP values is the cost of the sample "
            "symmetry that was imposed."
        ),
        height_in=2.3 * rows + 0.5,
    )


def tilt_misfit_figure(figures: Sequence[Mapping[str, Any]]) -> ResultFigure:
    """Mean recalculated-minus-measured misfit against tilt, per pole figure."""

    edges = np.arange(0.0, 95.0, 5.0)

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        for index, entry in enumerate(figures):
            polar = np.array([point["polar_deg"] for point in entry["points"]], dtype=float)
            difference = np.array([point["difference"] for point in entry["points"]], dtype=float)
            centres, means, spreads = [], [], []
            for low, high in pairwise(edges):
                inside = (polar >= low) & (polar < high)
                if np.count_nonzero(inside) >= 3:
                    centres.append(0.5 * (low + high))
                    means.append(float(np.mean(difference[inside])))
                    spreads.append(float(np.std(difference[inside])))
            color = _PALETTE[index % len(_PALETTE)]
            axes.errorbar(
                centres,
                means,
                yerr=spreads,
                fmt="o-",
                ms=3.5,
                lw=1.0,
                capsize=2,
                color=color,
                label=str(entry["label"]),
            )
        axes.axhline(0.0, color=COLORS["guide"], lw=0.8)
        axes.set_xlabel("Tilt from the specimen normal (°)")
        axes.set_ylabel("Recalculated − measured (m.r.d.)")
        axes.set_xlim(0.0, 90.0)
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="odf_tilt_misfit",
        title="Misfit against tilt",
        caption=(
            "Mean of recalculated minus measured intensity in 5° rings of tilt, with ±1 "
            "standard deviation of the points in each ring as error bars."
        ),
        interpretation=(
            "A misfit that grows towards the rim, with the recalculated figure above the "
            "measured one, is the signature of uncorrected defocusing and absorption in the "
            "measurement rather than of the ODF. A misfit confined to one ring is a component "
            "the ODF did not resolve."
        ),
        height_in=3.0,
    )


def fraction_figure(rows: Sequence[Mapping[str, Any]], *, tolerance_deg: float) -> ResultFigure:
    """Volume fraction of each ideal orientation against the random reference."""

    names = [str(row["component"]) for row in rows]
    measured = np.array([float(row["percent"]) for row in rows])
    random = np.array([float(row["random_percent"]) for row in rows])

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(names))
        axes.bar(positions - 0.2, measured, 0.4, color=COLORS["difference"], label="This ODF")
        axes.bar(positions + 0.2, random, 0.4, color=COLORS["guide"], label="Random specimen")
        for x, value, reference in zip(positions, measured, random, strict=True):
            if reference > 0:
                axes.text(
                    x - 0.2,
                    value,
                    f"×{value / reference:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                )
        axes.set_xticks(positions, labels=names, rotation=25, ha="right")
        axes.set_ylabel(f"Volume within {tolerance_deg:g}° (%)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="component_fractions",
        title="Ideal-orientation volume fractions",
        caption=(
            f"Volume of the ODF within {tolerance_deg:g}° of each ideal orientation and its "
            "symmetry equivalents (blue), against the volume a random specimen holds in the "
            "same ball (grey); the label is the ratio. The fractions carry a sampling "
            "uncertainty of a few percent of their value."
        ),
        interpretation=(
            "Only a ratio well above 1 is a texture component; a fraction near its random "
            "value is no evidence for the component, however large the percentage looks. "
            "Neighbouring balls can overlap, so the fractions need not sum to 100 %."
        ),
        height_in=3.0,
    )
