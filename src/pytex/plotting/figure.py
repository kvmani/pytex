"""Publication-figure foundation: rc bridging, panel grids, and export.

This module turns the YAML theme system (`pytex.plotting.styles`) into
matplotlib-wide behavior and provides the reusable composition primitives that
publication figures need:

- `rc_params_from_style` / `publication_style`: derive a complete rcParams set
  from a resolved theme (fonts, DPI, faces, grids, journal-compliant font
  embedding: text stays editable text in SVG, TrueType in PDF/PS).
- `PanelGrid`: a multi-panel figure composer with uniform per-panel sizing,
  automatic (a), (b), (c) panel labels, and shared colorbars.
- `add_scale_bar`: anchored micrograph-style scale bar for map axes.
- `export_figure`: one-call multi-format export (SVG/PDF/PNG) for manuscripts.

Every PyTex plotting entry point can build on these so figures across the
library share one publication-grade visual system.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from string import ascii_lowercase
from typing import Any

import numpy as np

from pytex.plotting.styles import resolve_style


def _require_matplotlib() -> Any:
    try:
        import matplotlib
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return matplotlib


def rc_params_from_style(style: dict[str, Any]) -> dict[str, Any]:
    """Translate a resolved PyTex style mapping into matplotlib rcParams.

    Covers fonts, figure/axes faces, DPI, grid and line weights, legend
    framing, and publication-grade export defaults: ``svg.fonttype: none``
    keeps SVG text editable in vector editors, and TrueType (type 42) fonts in
    PDF/PS satisfy journal submission systems.
    """

    common = style["common"]
    figure = common["figure"]
    font = common["font"]
    size = float(font["size"])
    return {
        "figure.dpi": float(figure["dpi"]),
        "figure.facecolor": figure["facecolor"],
        "axes.facecolor": figure["axes_facecolor"],
        "font.family": font["family"],
        "font.size": size,
        "axes.titlesize": float(font["title_size"]),
        "axes.labelsize": size,
        "xtick.labelsize": size - 1.0,
        "ytick.labelsize": size - 1.0,
        "legend.fontsize": size - 1.0,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "grid.alpha": float(figure["grid_alpha"]),
        "grid.linewidth": 0.6,
        "lines.linewidth": float(common["line"]["width"]),
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.85",
        "savefig.dpi": 300.0,
        "savefig.bbox": "tight",
        "savefig.facecolor": figure["facecolor"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }


@contextmanager
def publication_style(
    *,
    theme: str = "journal",
    style_path: str | Path | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> Iterator[dict[str, Any]]:
    """Context manager applying a theme's rcParams; yields the resolved style.

    Everything plotted inside the context inherits the publication rc set, so
    ad-hoc user plots and PyTex entry points compose into one visual system::

        with publication_style(theme="journal") as style:
            fig, ax = plt.subplots()
            ...
    """

    matplotlib = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    with matplotlib.rc_context(rc_params_from_style(style)):
        yield style


def label_panels(
    axes: Sequence[Any],
    *,
    fmt: str = "({})",
    fontsize: float | None = None,
    fontweight: str = "bold",
    x_offset: float = -0.02,
    y_offset: float = 1.02,
    color: str | None = None,
) -> list[Any]:
    """Add sequential (a), (b), (c) labels above the given axes.

    Labels sit just outside the top-left corner in axes coordinates, the
    standard multi-panel convention of Nature/Science/Acta journals. Returns
    the created text artists.
    """

    if len(axes) > len(ascii_lowercase):
        raise ValueError("label_panels supports at most 26 panels.")
    artists = []
    for letter, axis in zip(ascii_lowercase, axes, strict=False):
        artists.append(
            axis.text(
                x_offset,
                y_offset,
                fmt.format(letter),
                transform=axis.transAxes,
                fontsize=fontsize,
                fontweight=fontweight,
                color=color,
                ha="right",
                va="bottom",
            )
        )
    return artists


class PanelGrid:
    """A publication multi-panel figure with uniform panels and shared color.

    Creates a ``rows x cols`` grid sized per panel (not per figure), under the
    active or requested theme. Panels are exposed as a 2-D array (`axes`) and
    a flat list (`axes_flat`); `label` stamps (a), (b), (c) labels and
    `shared_colorbar` attaches one colorbar to any subset of panels.
    """

    def __init__(
        self,
        rows: int,
        cols: int,
        *,
        panel_size: tuple[float, float] = (3.4, 3.2),
        theme: str = "journal",
        style_path: str | Path | None = None,
        style_overrides: dict[str, Any] | None = None,
        sharex: bool = False,
        sharey: bool = False,
        constrained: bool = True,
    ) -> None:
        if rows < 1 or cols < 1:
            raise ValueError("PanelGrid needs at least one row and one column.")
        matplotlib = _require_matplotlib()
        import matplotlib.pyplot as plt

        self.style = resolve_style(
            theme=theme, style_path=style_path, overrides=style_overrides
        )
        rc = rc_params_from_style(self.style)
        with matplotlib.rc_context(rc):
            self.figure, axes = plt.subplots(
                rows,
                cols,
                figsize=(panel_size[0] * cols, panel_size[1] * rows),
                dpi=rc["figure.dpi"],
                facecolor=rc["figure.facecolor"],
                sharex=sharex,
                sharey=sharey,
                squeeze=False,
                layout="constrained" if constrained else None,
            )
        self.axes = np.asarray(axes, dtype=object)
        self._rc = rc

    @property
    def axes_flat(self) -> list[Any]:
        return list(self.axes.ravel())

    def label(self, **kwargs: Any) -> list[Any]:
        """Stamp (a), (b), (c) panel labels; forwards to `label_panels`."""

        return label_panels(self.axes_flat, **kwargs)

    def shared_colorbar(
        self,
        mappable: Any,
        *,
        axes: Sequence[Any] | None = None,
        label: str | None = None,
        location: str = "right",
        shrink: float = 0.85,
    ) -> Any:
        """Attach one colorbar for `mappable` spanning the given (or all) panels."""

        targets = list(axes) if axes is not None else self.axes_flat
        colorbar = self.figure.colorbar(
            mappable, ax=targets, location=location, shrink=shrink
        )
        if label is not None:
            colorbar.set_label(label)
        return colorbar

    def hide_unused(self, used: int) -> None:
        """Hide trailing panels beyond the first ``used`` (ragged last row)."""

        for axis in self.axes_flat[used:]:
            axis.set_visible(False)

    def export(self, path: str | Path, **kwargs: Any) -> tuple[Path, ...]:
        """Export the figure; forwards to `export_figure`."""

        return export_figure(self.figure, path, **kwargs)


def add_scale_bar(
    axis: Any,
    length: float,
    *,
    label: str | None = None,
    units: str = "µm",
    location: str = "lower right",
    color: str = "#111111",
    bar_height_fraction: float = 0.012,
    fontsize: float | None = None,
) -> Any:
    """Add an anchored micrograph-style scale bar to a map axis.

    ``length`` is in data units (the map's coordinate units, typically
    micrometres). The default label is ``"<length> <units>"``. Returns the
    anchored artist.
    """

    _require_matplotlib()
    import matplotlib.font_manager as font_manager
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    if length <= 0.0:
        raise ValueError("scale-bar length must be positive.")
    text = label if label is not None else f"{length:g} {units}"
    y_min, y_max = axis.get_ylim()
    bar = AnchoredSizeBar(
        axis.transData,
        length,
        text,
        location,
        pad=0.4,
        borderpad=0.6,
        sep=3.0,
        frameon=True,
        color=color,
        size_vertical=abs(y_max - y_min) * bar_height_fraction,
        fontproperties=font_manager.FontProperties(size=fontsize),
    )
    bar.patch.set_alpha(0.75)
    axis.add_artist(bar)
    return bar


_EXPORT_FORMATS = ("svg", "png", "pdf")


def export_figure(
    figure: Any,
    path: str | Path,
    *,
    formats: Sequence[str] = ("svg", "png"),
    dpi: float = 300.0,
    transparent: bool = False,
) -> tuple[Path, ...]:
    """Export a figure to one or more publication formats.

    ``path`` may be a bare stem (each requested format is written next to it)
    or carry an explicit ``.svg``/``.png``/``.pdf`` suffix (only that format
    is written). Parent directories are created. Returns the written paths.
    """

    target = Path(path)
    if target.suffix.lstrip(".").lower() in _EXPORT_FORMATS:
        stem = target.with_suffix("")
        requested: tuple[str, ...] = (target.suffix.lstrip(".").lower(),)
    else:
        stem = target
        requested = tuple(fmt.lower() for fmt in formats)
    unknown = set(requested) - set(_EXPORT_FORMATS)
    if unknown:
        raise ValueError(f"Unsupported export formats: {sorted(unknown)!r}")
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in requested:
        output = stem.with_suffix(f".{fmt}")
        figure.savefig(output, format=fmt, dpi=dpi, transparent=transparent)
        written.append(output)
    return tuple(written)


__all__ = [
    "PanelGrid",
    "add_scale_bar",
    "export_figure",
    "label_panels",
    "publication_style",
    "rc_params_from_style",
]
