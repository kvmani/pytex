from __future__ import annotations

from typing import Any

import numpy as np

from pytex.ebsd.models import (
    CrystalMap,
    GrainBoundaryNetwork,
    GrainSegmentation,
    _specimen_direction_vector,
)
from pytex.plotting.colormaps import categorical_colors, register_pytex_colormaps
from pytex.plotting.figure import add_scale_bar
from pytex.plotting.ipf import IPFColorKey
from pytex.plotting.styles import resolve_style


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, object


def _boundary_network_for_overlay(
    crystal_map: CrystalMap,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None,
) -> GrainBoundaryNetwork | None:
    if boundary_overlay is None:
        return None
    if isinstance(boundary_overlay, GrainSegmentation):
        if boundary_overlay.crystal_map is not crystal_map:
            raise ValueError("boundary_overlay.crystal_map must be this CrystalMap instance.")
        return boundary_overlay.boundary_network()
    if boundary_overlay.segmentation.crystal_map is not crystal_map:
        raise ValueError("boundary_overlay.segmentation.crystal_map must be this CrystalMap.")
    return boundary_overlay


def _regular_grid_extent(crystal_map: CrystalMap) -> tuple[float, float, float, float]:
    rows, cols = crystal_map._require_regular_2d_grid()
    if crystal_map.step_sizes is not None:
        dx, dy = crystal_map.step_sizes
    else:
        x_values = np.unique(crystal_map.coordinates[:, 0])
        y_values = np.unique(crystal_map.coordinates[:, 1])
        dx = float(np.min(np.diff(x_values))) if x_values.size > 1 else 1.0
        dy = float(np.min(np.diff(y_values))) if y_values.size > 1 else 1.0
    origin = crystal_map.coordinates[0, :2]
    xmin = float(origin[0] - 0.5 * dx)
    xmax = float(origin[0] + (cols - 0.5) * dx)
    ymin = float(origin[1] - 0.5 * dy)
    ymax = float(origin[1] + (rows - 0.5) * dy)
    return xmin, xmax, ymin, ymax


def _overlay_boundaries(
    ax: Any,
    crystal_map: CrystalMap,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None,
    *,
    color: str,
    linewidth: float,
) -> None:
    network = _boundary_network_for_overlay(crystal_map, boundary_overlay)
    if network is None:
        return
    for segment in network.segments:
        pair = crystal_map.coordinates[[segment.left_index, segment.right_index], :2]
        ax.plot(pair[:, 0], pair[:, 1], color=color, linewidth=linewidth, alpha=0.9)


def plot_ipf_map(
    crystal_map: CrystalMap,
    *,
    direction: str | np.ndarray = "z",
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    scale_bar: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    plt, _ = _require_matplotlib()
    if crystal_map.orientations.symmetry is None:
        raise ValueError("plot_ipf_map() requires crystal-map orientations with crystal symmetry.")
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    specimen_direction = _specimen_direction_vector(
        direction,
        crystal_map.orientations.specimen_frame,
    )
    colors = IPFColorKey(
        crystal_symmetry=crystal_map.orientations.symmetry,
        specimen_direction=specimen_direction,
    ).colors_from_orientations(crystal_map.orientations)
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    try:
        rows, cols = crystal_map._require_regular_2d_grid()
        image = colors.reshape((rows, cols, 3))
        axes.imshow(
            image,
            origin="lower",
            extent=_regular_grid_extent(crystal_map),
            interpolation="nearest",
        )
    except ValueError:
        axes.scatter(
            crystal_map.coordinates[:, 0],
            crystal_map.coordinates[:, 1],
            c=colors,
            s=float(common["marker"]["size"]) * 1.2,
            edgecolors="none",
        )
    _overlay_boundaries(
        axes,
        crystal_map,
        boundary_overlay,
        color=boundary_color,
        linewidth=boundary_linewidth,
    )
    if scale_bar is not None:
        add_scale_bar(axes, scale_bar)
    direction_label = direction if isinstance(direction, str) else "custom"
    axes.set_xlabel(crystal_map.map_frame.axes[0])
    axes.set_ylabel(crystal_map.map_frame.axes[1])
    axes.set_title(f"IPF Map ({direction_label})")
    axes.set_aspect("equal", adjustable="box")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig


def plot_kam_map(
    crystal_map: CrystalMap,
    *,
    symmetry_aware: bool = True,
    connectivity: int = 4,
    order: int = 1,
    threshold_deg: float | None = None,
    statistic: str = "mean",
    segmentation: GrainSegmentation | None = None,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    cmap: str = "pytex.misorientation",
    scale_bar: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    plt, _ = _require_matplotlib()
    register_pytex_colormaps()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    values = crystal_map.kernel_average_misorientation_deg(
        symmetry_aware=symmetry_aware,
        connectivity=connectivity,
        order=order,
        threshold_deg=threshold_deg,
        statistic=statistic,
        segmentation=segmentation,
    )
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    image = axes.imshow(
        values,
        origin="lower",
        extent=_regular_grid_extent(crystal_map),
        interpolation="nearest",
        cmap=cmap,
    )
    colorbar = fig.colorbar(image, ax=axes)
    colorbar.set_label("KAM (deg)")
    _overlay_boundaries(
        axes,
        crystal_map,
        boundary_overlay,
        color=boundary_color,
        linewidth=boundary_linewidth,
    )
    if scale_bar is not None:
        add_scale_bar(axes, scale_bar)
    axes.set_xlabel(crystal_map.map_frame.axes[0])
    axes.set_ylabel(crystal_map.map_frame.axes[1])
    axes.set_title("Kernel Average Misorientation")
    axes.set_aspect("equal", adjustable="box")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig


def plot_property_map(
    crystal_map: CrystalMap,
    name: str,
    *,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar_label: str | None = None,
    scale_bar: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Render a per-point scalar property channel (IQ, CI, BC, MAD, ...)."""

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    values = crystal_map.get_property(name)
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    try:
        image = axes.imshow(
            crystal_map.property_map(name),
            origin="lower",
            extent=_regular_grid_extent(crystal_map),
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    except ValueError:
        image = axes.scatter(
            crystal_map.coordinates[:, 0],
            crystal_map.coordinates[:, 1],
            c=values,
            s=float(common["marker"]["size"]) * 1.2,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            edgecolors="none",
        )
    colorbar = fig.colorbar(image, ax=axes)
    colorbar.set_label(colorbar_label if colorbar_label is not None else name)
    _overlay_boundaries(
        axes,
        crystal_map,
        boundary_overlay,
        color=boundary_color,
        linewidth=boundary_linewidth,
    )
    if scale_bar is not None:
        add_scale_bar(axes, scale_bar)
    axes.set_xlabel(crystal_map.map_frame.axes[0])
    axes.set_ylabel(crystal_map.map_frame.axes[1])
    axes.set_title(f"Property Map ({name})")
    axes.set_aspect("equal", adjustable="box")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig


def plot_phase_map(
    crystal_map: CrystalMap,
    *,
    colors: dict[str, str] | None = None,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    scale_bar: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Render a phase map, coloring each point by its phase assignment."""

    plt, _ = _require_matplotlib()
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    phase_ids = crystal_map.phase_id_array
    if phase_ids is None:
        raise ValueError("plot_phase_map() requires a CrystalMap with phase assignments.")
    entries = crystal_map.resolved_phase_entries
    # fixed-order CVD-safe identity palette; a phase keeps its color as
    # maps gain or lose phases
    default_cycle = categorical_colors(len(entries))
    ordered_ids = [entry.phase_id for entry in entries]
    id_to_index = {phase_id: index for index, phase_id in enumerate(ordered_ids)}
    palette: list[str] = []
    for index, entry in enumerate(entries):
        if colors is not None and entry.name in colors:
            palette.append(colors[entry.name])
        else:
            palette.append(default_cycle[index])
    indexed = np.array([id_to_index[int(value)] for value in phase_ids], dtype=np.float64)
    colormap = ListedColormap(palette)
    norm = BoundaryNorm(np.arange(-0.5, len(entries) + 0.5, 1.0), colormap.N)
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    try:
        rows, cols = crystal_map._require_regular_2d_grid()
        axes.imshow(
            indexed.reshape((rows, cols)),
            origin="lower",
            extent=_regular_grid_extent(crystal_map),
            interpolation="nearest",
            cmap=colormap,
            norm=norm,
        )
    except ValueError:
        axes.scatter(
            crystal_map.coordinates[:, 0],
            crystal_map.coordinates[:, 1],
            c=indexed,
            s=float(common["marker"]["size"]) * 1.2,
            cmap=colormap,
            norm=norm,
            edgecolors="none",
        )
    legend_handles = [
        Patch(facecolor=palette[index], edgecolor="none", label=entry.name)
        for index, entry in enumerate(entries)
    ]
    axes.legend(handles=legend_handles, loc="upper right", framealpha=0.9)
    _overlay_boundaries(
        axes,
        crystal_map,
        boundary_overlay,
        color=boundary_color,
        linewidth=boundary_linewidth,
    )
    if scale_bar is not None:
        add_scale_bar(axes, scale_bar)
    axes.set_xlabel(crystal_map.map_frame.axes[0])
    axes.set_ylabel(crystal_map.map_frame.axes[1])
    axes.set_title("Phase Map")
    axes.set_aspect("equal", adjustable="box")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig


def plot_ipf_xyz_maps(
    crystal_map: CrystalMap,
    *,
    directions: tuple[str, ...] = ("x", "y", "z"),
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> Any:
    """Render a side-by-side triptych of IPF maps for the given sample directions."""

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    width, height = common["figure"]["figsize"]
    fig, axes_row = plt.subplots(
        1,
        len(directions),
        figsize=(float(width) * len(directions), float(height)),
        dpi=int(common["figure"]["dpi"]),
        facecolor=common["figure"]["facecolor"],
    )
    axes_list = np.atleast_1d(axes_row).ravel().tolist()
    for axis, direction in zip(axes_list, directions, strict=True):
        plot_ipf_map(
            crystal_map,
            direction=direction,
            boundary_overlay=boundary_overlay,
            boundary_color=boundary_color,
            boundary_linewidth=boundary_linewidth,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
            ax=axis,
        )
    fig.tight_layout()
    return fig
