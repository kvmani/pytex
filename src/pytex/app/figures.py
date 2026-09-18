"""Publication-quality figures for application results, drawn once on the server.

Purpose
-------
Every analysis in the workbench reaches its answer through intermediate results
— a scan, the peaks found in it, a fit, its residuals — and a reader can only
judge the answer by seeing them. This module is the one place those figures are
drawn. A service describes *what* to plot with the helpers below and gets back a
:class:`~pytex.app.results.ResultFigure`; the generic renderer shows it, offers
it as SVG and PNG, and the Markdown report embeds it. No panel draws its own
copy, so the figure on screen, in the download and in the report cannot differ.

Conventions
-----------
- Vector SVG with text drawn as paths, so a figure renders identically on a
  machine without the font, in an offline report viewer, and when rasterized in
  the browser. No external resource is ever referenced.
- The house style is the journal theme of :mod:`pytex.plotting.figure`, with a
  sans-serif face that ships with matplotlib.
- Measured data is drawn in near-black, a model in orange, a difference in blue,
  and guide lines in grey. The pairing survives greyscale print because data are
  points or thin lines and the model is a thicker line.
- Error bars are one standard uncertainty unless a caption says otherwise.
- Figures are drawn on a bare :class:`matplotlib.figure.Figure`, never through
  ``pyplot``, so no global figure registry is touched: the server is threaded,
  and nothing here can leak an open figure into a test.
"""

from __future__ import annotations

import io
import re
from collections.abc import Callable, Sequence
from contextlib import ExitStack
from typing import Any

import numpy as np

from pytex.app.results import ResultFigure

__all__ = [
    "COLORS",
    "draw_correlation",
    "draw_image",
    "draw_normalized_residuals",
    "draw_observed_model",
    "draw_residuals",
    "draw_ticks",
    "render_figure",
]

#: The palette every application figure uses, by role rather than by hue.
COLORS: dict[str, str] = {
    "data": "#111827",
    "model": "#d97706",
    "difference": "#2563eb",
    "background": "#64748b",
    "guide": "#94a3b8",
    "accent": "#0f766e",
    "warning": "#dc2626",
    "band_2": "#dbeafe",
    "band_3": "#fef3c7",
    "unassigned": "#dc2626",
}

#: Width, in inches, of a figure that spans the result column.
FULL_WIDTH_IN = 6.4

_RC_OVERRIDES: dict[str, Any] = {
    "font.family": "DejaVu Sans",
    "font.size": 9.0,
    "axes.titlesize": 10.0,
    "axes.labelsize": 9.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 8.0,
    "lines.linewidth": 1.2,
    "svg.fonttype": "path",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "savefig.bbox": "standard",
}


def render_figure(
    draw: Callable[[Any], None],
    *,
    key: str,
    title: str,
    caption: str,
    interpretation: str | None = None,
    width_in: float = FULL_WIDTH_IN,
    height_in: float = 3.6,
) -> ResultFigure:
    """Draw one figure and return it as a :class:`ResultFigure`.

    Parameters
    ----------
    draw : callable
        Receives an empty, constrained-layout :class:`matplotlib.figure.Figure`
        and adds its axes and artists. Called inside the house style, so every
        axis it creates inherits the fonts and line weights.
    key, title, caption, interpretation
        As for :class:`ResultFigure`.
    width_in, height_in : float
        Physical size of the figure in inches.

    Returns
    -------
    ResultFigure
        The figure as self-contained SVG, with its caption and reading.
    """

    import matplotlib
    from matplotlib.figure import Figure

    from pytex.plotting.figure import publication_style

    with ExitStack() as stack:
        stack.enter_context(publication_style(theme="journal"))
        stack.enter_context(
            matplotlib.rc_context({**_RC_OVERRIDES, "svg.hashsalt": f"pytex-{key}"})
        )
        figure = Figure(figsize=(width_in, height_in), layout="constrained")
        draw(figure)
        buffer = io.StringIO()
        figure.savefig(buffer, format="svg", metadata={"Date": None})
    markup = buffer.getvalue()
    # Only the <svg> element: the XML declaration and DOCTYPE are invalid inside
    # an HTML document and inside a data URL wrapped by a Markdown viewer.
    markup = markup[markup.index("<svg") :]
    markup = re.sub(r"\n\s*\n", "\n", markup)
    return ResultFigure(
        key=key,
        title=title,
        svg=markup,
        caption=caption,
        interpretation=interpretation,
        width_in=float(width_in),
        height_in=float(height_in),
    )


def draw_observed_model(
    axes: Any,
    x: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    model: Sequence[float] | np.ndarray | None = None,
    *,
    difference_axes: Any | None = None,
    observed_label: str = "Observed",
    model_label: str = "Model",
    observed_as_points: bool = False,
    observed_sigma: Sequence[float] | np.ndarray | None = None,
) -> None:
    """Overlay measured data on its model, with the difference below.

    The difference is drawn on its own axes, sharing the abscissa and at the
    same intensity scale, because the Rietveld convention of a difference curve
    under the pattern is the one plot that shows *where* a model fails.
    """

    x_values = np.asarray(x, dtype=float)
    y_observed = np.asarray(observed, dtype=float)
    if observed_as_points:
        if observed_sigma is not None:
            axes.errorbar(
                x_values,
                y_observed,
                yerr=np.asarray(observed_sigma, dtype=float),
                fmt="o",
                ms=2.5,
                lw=0.6,
                color=COLORS["data"],
                label=observed_label,
            )
        else:
            axes.plot(x_values, y_observed, "o", ms=2.0, color=COLORS["data"], label=observed_label)
    else:
        axes.plot(x_values, y_observed, lw=0.7, color=COLORS["data"], label=observed_label)
    if model is not None:
        y_model = np.asarray(model, dtype=float)
        axes.plot(x_values, y_model, lw=1.3, color=COLORS["model"], label=model_label)
        if difference_axes is not None:
            difference = y_observed - y_model
            difference_axes.plot(x_values, difference, lw=0.7, color=COLORS["difference"])
            difference_axes.axhline(0.0, color=COLORS["guide"], lw=0.6)
            difference_axes.set_ylabel("Obs − model")
    axes.legend(loc="best", frameon=False)


def draw_residuals(
    axes: Any,
    x: Sequence[float] | np.ndarray,
    residual: Sequence[float] | np.ndarray,
    sigma: Sequence[float] | np.ndarray | None = None,
    *,
    labels: Sequence[str] | None = None,
    xlabel: str = "",
    ylabel: str = "",
) -> None:
    """Residuals with one-sigma error bars about a zero line.

    A residual is only large or small relative to its own uncertainty, so the
    bars are drawn whenever they are known: a point whose bar crosses zero is
    consistent with the model.
    """

    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(residual, dtype=float)
    axes.axhline(0.0, color=COLORS["guide"], lw=0.8)
    axes.errorbar(
        x_values,
        y_values,
        yerr=None if sigma is None else np.asarray(sigma, dtype=float),
        fmt="o",
        ms=4.0,
        lw=0.9,
        capsize=2.5,
        color=COLORS["data"],
        ecolor=COLORS["background"],
    )
    if labels is not None:
        _label_points(axes, x_values, y_values, labels)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)


def draw_normalized_residuals(
    axes: Any,
    x: Sequence[float] | np.ndarray,
    normalized: Sequence[float] | np.ndarray,
    *,
    labels: Sequence[str] | None = None,
    xlabel: str = "",
    ylabel: str = "",
) -> None:
    """Residuals in units of their own standard uncertainty, with 2σ and 3σ guides.

    Under a correct model with correct uncertainties about 95 % of the points lie
    inside ±2 and almost all inside ±3; points outside ±3 are outliers worth
    naming. The shaded bands make that judgement at a glance.
    """

    x_values = np.asarray(x, dtype=float)
    z_values = np.asarray(normalized, dtype=float)
    finite = z_values[np.isfinite(z_values)]
    span = max(3.8, float(np.max(np.abs(finite))) * 1.15 if finite.size else 3.8)
    axes.axhspan(-3.0, 3.0, color=COLORS["band_3"], lw=0, zorder=0)
    axes.axhspan(-2.0, 2.0, color=COLORS["band_2"], lw=0, zorder=0)
    for level, style in ((2.0, "--"), (3.0, ":")):
        axes.axhline(level, color=COLORS["guide"], lw=0.8, ls=style, zorder=1)
        axes.axhline(-level, color=COLORS["guide"], lw=0.8, ls=style, zorder=1)
    axes.axhline(0.0, color=COLORS["guide"], lw=0.8, zorder=1)
    outside = np.abs(z_values) > 3.0
    axes.plot(x_values[~outside], z_values[~outside], "o", ms=4.5, color=COLORS["data"], zorder=3)
    if np.any(outside):
        axes.plot(
            x_values[outside],
            z_values[outside],
            "D",
            ms=5.0,
            color=COLORS["warning"],
            zorder=3,
            label="beyond ±3σ",
        )
        axes.legend(loc="best", frameon=False)
    if labels is not None:
        _label_points(axes, x_values, z_values, labels)
    axes.set_ylim(-span, span)
    right = axes.get_xlim()[1]
    axes.text(right, 2.0, " ±2σ", va="center", ha="left", fontsize=7, color=COLORS["background"])
    axes.text(right, 3.0, " ±3σ", va="center", ha="left", fontsize=7, color=COLORS["background"])
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)


def draw_correlation(
    figure: Any,
    axes: Any,
    matrix: np.ndarray,
    names: Sequence[str],
) -> None:
    """A parameter-correlation matrix as an annotated diverging heat map.

    Fixed to the full ``[-1, 1]`` scale so a pale cell always means a weak
    correlation, whatever the other entries are.
    """

    values = np.asarray(matrix, dtype=float)
    image = axes.imshow(values, cmap="RdBu_r", vmin=-1.0, vmax=1.0)
    size = values.shape[0]
    axes.set_xticks(range(size), labels=list(names))
    axes.set_yticks(range(size), labels=list(names))
    axes.tick_params(length=0)
    for spine in axes.spines.values():
        spine.set_visible(False)
    for row in range(size):
        for column in range(size):
            value = values[row, column]
            axes.text(
                column,
                row,
                f"{value:+.3f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if abs(value) > 0.6 else COLORS["data"],
            )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.85)
    colorbar.set_label("Correlation coefficient")


def draw_image(
    figure: Any,
    axes: Any,
    image: np.ndarray,
    *,
    extent: Sequence[float] | None = None,
    cmap: str = "gray",
    colorbar_label: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    origin: str = "lower",
    xlabel: str = "",
    ylabel: str = "",
) -> None:
    """A two-dimensional map in physical coordinates, with a labelled colour bar."""

    artist = axes.imshow(
        np.asarray(image, dtype=float),
        extent=None if extent is None else tuple(extent),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin=origin,
        interpolation="nearest",
        aspect="equal",
    )
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if colorbar_label is not None:
        colorbar = figure.colorbar(artist, ax=axes, shrink=0.85)
        colorbar.set_label(colorbar_label)


def draw_ticks(
    axes: Any,
    positions: Sequence[float] | np.ndarray,
    *,
    labels: Sequence[str] | None = None,
    color: str = COLORS["accent"],
    height: float = 0.08,
    label: str | None = None,
) -> None:
    """Short vertical ticks along the bottom of an axis, one per predicted line.

    Drawn in axis-fraction height so they stay the same size whatever the
    intensity scale, which is how reflection markers are conventionally shown
    under a powder pattern.
    """

    transform = axes.get_xaxis_transform()
    values = np.asarray(positions, dtype=float)
    for index, value in enumerate(values):
        axes.plot(
            [value, value],
            [0.0, height],
            color=color,
            lw=0.9,
            transform=transform,
            label=label if index == 0 else None,
        )
        if labels is not None and index < len(labels) and labels[index]:
            axes.text(
                value,
                height + 0.01,
                labels[index],
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=color,
                transform=transform,
            )


def _label_points(axes: Any, x: np.ndarray, y: np.ndarray, labels: Sequence[str]) -> None:
    for x_value, y_value, text in zip(x, y, labels, strict=False):
        if not text or not np.isfinite(y_value):
            continue
        axes.annotate(
            text,
            (x_value, y_value),
            xytext=(5, 3),
            textcoords="offset points",
            ha="left",
            fontsize=6.5,
            color=COLORS["background"],
        )
