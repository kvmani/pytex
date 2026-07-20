"""Publication-grade rendering for composite OR SAED patterns.

Purpose: a typed, layered, highly configurable matplotlib renderer for
:class:`pytex.diffraction.composite.CompositeSAEDPattern`. Every visual
aspect — marker shape, fill, size mapping, color, opacity, edge, z-order,
variant subset, axes units, legend, transmitted beam — is controlled by an
explicit frozen configuration object with defaults tuned for publication
figures (open indigo parent circles over filled colored variant markers on a
white background, equal-aspect detector axes).

Rendering is structural: each sub-pattern becomes one scatter collection
tagged with a machine-readable ``gid`` (``pytex-composite:<label>``), which
is also how the regression tests verify the figure semantically (per repo
policy, no image-byte baselines).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, Literal

import numpy as np

from pytex.diffraction.composite import CompositeSAEDPattern, VariantZonePattern
from pytex.diffraction.kinematic import SpotTable

SizeModeName = Literal["intensity_area", "intensity_radius", "constant"]
AxesUnitsName = Literal["mm", "inv_angstrom"]

GID_PREFIX = "pytex-composite"

#: Colorblind-aware categorical palette for variant sub-patterns (12 hues:
#: Paul Tol's muted scheme extended with two high-contrast extras).
VARIANT_COLOR_PALETTE: tuple[str, ...] = (
    "#cc6677",
    "#332288",
    "#ddcc77",
    "#117733",
    "#88ccee",
    "#882255",
    "#44aa99",
    "#999933",
    "#aa4499",
    "#dd7788",
    "#6699cc",
    "#661100",
)

#: Marker cycle for variant sub-patterns; combined with the color palette the
#: default styling distinguishes up to 96 variants before repeating.
VARIANT_MARKER_CYCLE: tuple[str, ...] = ("s", "^", "D", "v", "P", "X", "h", "*")


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt


@dataclass(frozen=True, slots=True)
class SpotStyle:
    """Appearance of one sub-pattern's diffraction spots.

    ``size_mode`` maps normalized intensity to the matplotlib scatter size
    (points^2): ``"intensity_area"`` makes marker **area** proportional to
    intensity (the perceptually honest default), ``"intensity_radius"`` makes
    marker radius proportional to intensity (dramatizes strong spots), and
    ``"constant"`` ignores intensity. ``min_size_pt2`` keeps weak reflections
    visible. ``filled=False`` renders hollow markers (edge in ``color``),
    the classic convention for parent-phase overlays.
    """

    marker: str = "o"
    color: str = "#3f51b5"
    filled: bool = True
    size_scale: float = 90.0
    size_mode: SizeModeName = "intensity_area"
    min_size_pt2: float = 6.0
    alpha: float = 0.9
    edge_color: str = "#111111"
    edge_width: float = 0.6
    zorder: float = 3.0

    def __post_init__(self) -> None:
        if not self.marker:
            raise ValueError("SpotStyle.marker must be non-empty.")
        if not np.isfinite(self.size_scale) or self.size_scale <= 0.0:
            raise ValueError("SpotStyle.size_scale must be finite and strictly positive.")
        if self.size_mode not in {"intensity_area", "intensity_radius", "constant"}:
            raise ValueError(
                "SpotStyle.size_mode must be 'intensity_area', 'intensity_radius' or "
                "'constant'."
            )
        if not np.isfinite(self.min_size_pt2) or self.min_size_pt2 <= 0.0:
            raise ValueError("SpotStyle.min_size_pt2 must be finite and strictly positive.")
        if not 0.0 < self.alpha <= 1.0:
            raise ValueError("SpotStyle.alpha must lie in the interval (0, 1].")
        if not np.isfinite(self.edge_width) or self.edge_width < 0.0:
            raise ValueError("SpotStyle.edge_width must be finite and non-negative.")

    def marker_sizes_pt2(self, intensity: np.ndarray) -> np.ndarray:
        """Scatter sizes (points^2) for max-normalized intensities."""

        values = np.asarray(intensity, dtype=np.float64)
        if self.size_mode == "constant":
            sizes = np.full(values.shape, self.size_scale)
        elif self.size_mode == "intensity_area":
            sizes = self.size_scale * values
        else:
            sizes = self.size_scale * values**2
        return np.asarray(np.maximum(sizes, self.min_size_pt2), dtype=np.float64)


_DEFAULT_PARENT_STYLE = SpotStyle(
    marker="o",
    color="#3f51b5",
    filled=False,
    size_scale=140.0,
    edge_width=1.4,
    alpha=1.0,
    zorder=5.0,
)
_DEFAULT_CHILD_STYLE = SpotStyle(
    marker="s",
    color="#f57c00",
    filled=True,
    size_scale=80.0,
    edge_width=0.5,
    alpha=0.85,
    zorder=3.0,
)


@dataclass(frozen=True, slots=True)
class CompositeSAEDPlotConfig:
    """Full rendering configuration for composite SAED figures.

    ``variant_styles`` overrides the palette/marker cycling per 1-based
    variant index; unlisted variants fall back to ``child_base_style`` with
    the color palette and marker cycle applied in rendering order.
    ``variant_indices`` restricts rendering to a subset without re-simulating.
    ``axes_units`` switches between calibrated detector millimetres and
    reciprocal angstrom. The transmitted beam is drawn as a distinct central
    marker (never part of any spot table).
    """

    parent_style: SpotStyle = field(default_factory=lambda: _DEFAULT_PARENT_STYLE)
    child_base_style: SpotStyle = field(default_factory=lambda: _DEFAULT_CHILD_STYLE)
    variant_styles: Mapping[int, SpotStyle] | None = None
    variant_color_palette: tuple[str, ...] = VARIANT_COLOR_PALETTE
    variant_marker_cycle: tuple[str, ...] = VARIANT_MARKER_CYCLE
    variant_indices: tuple[int, ...] | None = None
    show_parent: bool = True
    show_transmitted_beam: bool = True
    transmitted_beam_size_pt2: float = 240.0
    transmitted_beam_color: str = "#111111"
    axes_units: AxesUnitsName = "mm"
    show_legend: bool = True
    legend_max_entries: int = 13
    legend_outside: bool = True
    show_title: bool = True
    title: str | None = None
    background: str = "#ffffff"
    figsize: tuple[float, float] = (7.0, 7.0)
    dpi: int = 200
    limit_padding_fraction: float = 0.08

    def __post_init__(self) -> None:
        if not self.variant_color_palette:
            raise ValueError("variant_color_palette must contain at least one color.")
        if not self.variant_marker_cycle:
            raise ValueError("variant_marker_cycle must contain at least one marker.")
        if (
            not np.isfinite(self.transmitted_beam_size_pt2)
            or self.transmitted_beam_size_pt2 <= 0.0
        ):
            raise ValueError("transmitted_beam_size_pt2 must be finite and strictly positive.")
        if self.axes_units not in {"mm", "inv_angstrom"}:
            raise ValueError("axes_units must be 'mm' or 'inv_angstrom'.")
        if self.legend_max_entries < 1:
            raise ValueError("legend_max_entries must be at least 1.")
        if not 0.0 <= self.limit_padding_fraction < 1.0:
            raise ValueError("limit_padding_fraction must lie in the interval [0, 1).")
        if self.dpi < 1:
            raise ValueError("dpi must be at least 1.")

    def style_for_variant(self, variant_index: int, render_position: int) -> SpotStyle:
        """Resolved style for a variant: explicit override or cycled default."""

        if self.variant_styles is not None and variant_index in self.variant_styles:
            return self.variant_styles[variant_index]
        color = self.variant_color_palette[
            render_position % len(self.variant_color_palette)
        ]
        marker = self.variant_marker_cycle[
            render_position % len(self.variant_marker_cycle)
        ]
        return replace(self.child_base_style, color=color, marker=marker)


def _spot_coordinates(table: SpotTable, axes_units: AxesUnitsName) -> np.ndarray:
    if axes_units == "mm":
        return table.detector_mm
    return table.g_detector_inv_angstrom


def _scatter_kwargs(style: SpotStyle) -> dict[str, Any]:
    if style.filled:
        return {
            "facecolors": style.color,
            "edgecolors": style.edge_color,
            "linewidths": style.edge_width,
        }
    return {
        "facecolors": "none",
        "edgecolors": style.color,
        "linewidths": max(style.edge_width, 0.8),
    }


def _scatter_sub_pattern(
    ax: Any,
    table: SpotTable,
    style: SpotStyle,
    *,
    label: str,
    gid: str,
    axes_units: AxesUnitsName,
) -> None:
    coordinates = _spot_coordinates(table, axes_units)
    collection = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        s=style.marker_sizes_pt2(table.intensity),
        marker=style.marker,
        alpha=style.alpha,
        zorder=style.zorder,
        label=label,
        **_scatter_kwargs(style),
    )
    collection.set_gid(gid)


def render_composite_saed(
    pattern: CompositeSAEDPattern,
    *,
    config: CompositeSAEDPlotConfig | None = None,
    ax: Any | None = None,
) -> Any:
    """Render a composite OR SAED pattern to a matplotlib figure.

    Purpose: the presentation layer over
    :func:`pytex.diffraction.composite.simulate_composite_saed`. Children are
    drawn first (palette/marker cycling or explicit per-variant styles), the
    parent on top as hollow markers, then the transmitted beam at the origin;
    axes are equal-aspect with units chosen by the configuration.

    Inputs: ``config`` — a :class:`CompositeSAEDPlotConfig` (defaults are
    publication-ready); ``ax`` — optional existing matplotlib axes (its
    figure is reused), otherwise a new figure is created with the configured
    size, dpi and background.

    Output: the matplotlib ``Figure``. Every sub-pattern scatter carries the
    gid ``pytex-composite:<label>`` for structural inspection and testing.
    """

    plot_config = CompositeSAEDPlotConfig() if config is None else config
    plt = _require_matplotlib()

    rendered = pattern
    if plot_config.variant_indices is not None:
        rendered = pattern.select_variants(list(plot_config.variant_indices))

    if ax is None:
        fig, axes = plt.subplots(
            figsize=plot_config.figsize,
            dpi=plot_config.dpi,
            facecolor=plot_config.background,
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(plot_config.background)

    child_phase_name = rendered.relationship.child_phase.name
    variant_patterns: tuple[VariantZonePattern, ...] = rendered.variant_patterns
    for position, variant_pattern in enumerate(variant_patterns):
        style = plot_config.style_for_variant(variant_pattern.variant_index, position)
        label = f"{child_phase_name} {variant_pattern.label()}"
        _scatter_sub_pattern(
            axes,
            variant_pattern.spots,
            style,
            label=label,
            gid=f"{GID_PREFIX}:variant:{variant_pattern.variant_index}",
            axes_units=plot_config.axes_units,
        )

    parent_table = rendered.parent_spots
    if plot_config.show_parent and parent_table is not None:
        parent_zone = "[" + " ".join(
            str(int(v)) for v in rendered.parent_zone_axis.indices
        ) + "]"
        parent_label = f"{rendered.relationship.parent_phase.name} {parent_zone}"
        _scatter_sub_pattern(
            axes,
            parent_table,
            plot_config.parent_style,
            label=parent_label,
            gid=f"{GID_PREFIX}:parent",
            axes_units=plot_config.axes_units,
        )

    if plot_config.show_transmitted_beam:
        beam = axes.scatter(
            [0.0],
            [0.0],
            s=plot_config.transmitted_beam_size_pt2,
            marker="o",
            facecolors=plot_config.transmitted_beam_color,
            edgecolors="none",
            zorder=6.0,
            label="transmitted beam",
        )
        beam.set_gid(f"{GID_PREFIX}:transmitted-beam")

    coordinates = rendered.all_detector_coordinates()
    if plot_config.axes_units == "inv_angstrom":
        coordinates = (
            coordinates / rendered.config.camera_constant_mm_angstrom
            if coordinates.size
            else coordinates
        )
    if coordinates.size:
        extent = float(np.max(np.abs(coordinates)))
    else:
        extent = 1.0
    extent = extent * (1.0 + plot_config.limit_padding_fraction) or 1.0
    axes.set_xlim(-extent, extent)
    axes.set_ylim(-extent, extent)
    axes.set_aspect("equal", adjustable="box")

    if plot_config.axes_units == "mm":
        axes.set_xlabel("detector u (mm)")
        axes.set_ylabel("detector v (mm)")
    else:
        axes.set_xlabel(r"$g \cdot u$ ($\mathrm{\AA}^{-1}$)")
        axes.set_ylabel(r"$g \cdot v$ ($\mathrm{\AA}^{-1}$)")

    if plot_config.show_title:
        if plot_config.title is not None:
            axes.set_title(plot_config.title)
        else:
            parent_zone = "[" + " ".join(
                str(int(v)) for v in rendered.parent_zone_axis.indices
            ) + "]"
            axes.set_title(
                f"{rendered.relationship.name}: composite SAED, parent {parent_zone} zone"
            )

    if plot_config.show_legend:
        handles, labels = axes.get_legend_handles_labels()
        if len(handles) > plot_config.legend_max_entries:
            handles = handles[: plot_config.legend_max_entries]
            labels = labels[: plot_config.legend_max_entries]
        if handles:
            legend_kwargs: dict[str, Any] = {"fontsize": 8, "framealpha": 0.9}
            if plot_config.legend_outside:
                legend_kwargs.update({"loc": "upper left", "bbox_to_anchor": (1.02, 1.0)})
            axes.legend(handles, labels, **legend_kwargs)

    fig.tight_layout()
    return fig


__all__ = [
    "GID_PREFIX",
    "VARIANT_COLOR_PALETTE",
    "VARIANT_MARKER_CYCLE",
    "CompositeSAEDPlotConfig",
    "SpotStyle",
    "render_composite_saed",
]
