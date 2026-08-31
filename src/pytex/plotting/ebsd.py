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


def _matplotlib() -> tuple[Any, Any]:
    import matplotlib.pyplot as plt

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
    if not network.segments:
        return
    import matplotlib
    from matplotlib.collections import LineCollection

    # One artist for the whole network, not one per segment. A real scan
    # carries tens of thousands of boundary faces, and a Line2D apiece costs
    # both the construction and every redraw the figure ever does.
    endpoints = np.array(
        [[segment.left_index, segment.right_index] for segment in network.segments],
        dtype=np.int64,
    )
    coordinates = np.asarray(crystal_map.coordinates, dtype=np.float64)[:, :2]
    ax.add_collection(
        LineCollection(
            list(coordinates[endpoints]),
            colors=color,
            linewidths=linewidth,
            alpha=0.9,
            # A collection defaults to butt caps and a line to whatever
            # `lines.solid_capstyle` says, which is `projecting`. Taking the
            # line's default keeps the drawn boundary the length it has always
            # been — the difference is a whole pixel at each end of a face.
            capstyle=matplotlib.rcParams["lines.solid_capstyle"],
        )
    )


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
    """Plot an inverse-pole-figure coloured orientation map.

    Purpose
    -------
    The standard EBSD map view: every point is coloured by which crystal
    direction lies along a chosen specimen axis, folded into the symmetry
    fundamental sector so the colour encodes crystallography rather than the
    arbitrary symmetry branch a measurement reported.

    Parameters
    ----------
    crystal_map : CrystalMap
        Must be on a regular 2-D grid.
    Remaining parameters select the reference direction, the colour key, and
    the Matplotlib target axes.

    Returns
    -------
    Any
        The Matplotlib axes.
    """

    plt, _ = _matplotlib()
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
    connectivity: int | None = None,
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
    """Plot a kernel-average-misorientation map.

    Purpose
    -------
    Display local misorientation in degrees — the standard visualization of
    stored deformation and subgrain structure. See
    :meth:`~pytex.ebsd.CrystalMap.kernel_average_misorientation_deg` for the
    threshold and neighbourhood choices, which materially change what the map
    shows.

    Returns
    -------
    Any
        The Matplotlib axes, with a labelled colorbar in degrees.
    """

    plt, _ = _matplotlib()
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
    if values.ndim == 2:
        image = axes.imshow(
            values,
            origin="lower",
            extent=_regular_grid_extent(crystal_map),
            interpolation="nearest",
            cmap=cmap,
        )
    else:
        image = axes.scatter(
            crystal_map.coordinates[:, 0],
            crystal_map.coordinates[:, 1],
            c=values,
            s=float(common["marker"]["size"]) * 1.2,
            cmap=cmap,
            edgecolors="none",
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


def plot_gnd_density_map(
    crystal_map: CrystalMap,
    *,
    burgers_vector_nm: float,
    method: str = "curvature",
    step_scale_m: float = 1e-6,
    kam_threshold_deg: float | None = 5.0,
    log_scale: bool = True,
    boundary_overlay: GrainSegmentation | GrainBoundaryNetwork | None = None,
    boundary_color: str = "#111111",
    boundary_linewidth: float = 0.85,
    cmap: str = "inferno",
    scale_bar: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot a geometrically necessary dislocation density map.

    Purpose
    -------
    Display the dislocation content implied by the measured lattice curvature,
    in dislocations per square metre — the map that connects an orientation
    measurement to stored energy, work hardening, and recrystallization driving
    force.

    Parameters
    ----------
    crystal_map : CrystalMap
        Regular 2-D grid carrying ``step_sizes``.
    burgers_vector_nm : float
        Burgers vector magnitude in nanometres.
    method : str
        ``"curvature"`` (Nye route, default) or ``"kam"``; see
        :func:`~pytex.ebsd.gnd.geometrically_necessary_dislocation_density`.
    step_scale_m : float
        Metres per map coordinate unit; the default treats them as micrometres.
    kam_threshold_deg : float, optional
        Grain-boundary exclusion threshold for the ``"kam"`` method.
    log_scale : bool
        Plot ``log10(density)`` (default). GND densities span orders of
        magnitude, so a linear scale is dominated by boundary artefacts and
        shows nothing of the grain interiors.
    boundary_overlay, boundary_color, boundary_linewidth :
        Optional grain-boundary overlay, as for :func:`plot_kam_map`. Strongly
        recommended: the curvature method does not exclude boundaries, so
        overlaying them shows the reader which features are artefacts.
    cmap, scale_bar, theme, style_path, style_overrides, ax :
        Styling and target axes.

    Returns
    -------
    Any
        The figure, with a colorbar labelled in the units actually plotted.

    Notes
    -----
    The plotted density is a **lower bound** and is **resolution dependent**;
    see the module documentation of :mod:`pytex.ebsd.gnd`. The step size should
    be quoted in any figure caption.
    """

    plt, _ = _matplotlib()
    from pytex.ebsd.gnd import geometrically_necessary_dislocation_density

    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    density = np.asarray(
        geometrically_necessary_dislocation_density(
            crystal_map,
            burgers_vector_nm=burgers_vector_nm,
            method=method,  # type: ignore[arg-type]
            step_scale_m=step_scale_m,
            kam_threshold_deg=kam_threshold_deg,
        ),
        dtype=np.float64,
    )
    if log_scale:
        # Zero density is a real, meaningful value (an unbent lattice) but has no
        # logarithm, so those points are left blank rather than clamped to an
        # arbitrary floor that would read as a measurement. ``out`` is required
        # alongside ``where``: without it the untouched entries are uninitialized.
        values = np.full_like(density, np.nan)
        positive = density > 0.0
        np.log10(density, out=values, where=positive)
        label = "log$_{10}$ GND density (m$^{-2}$)"
    else:
        values = density
        label = "GND density (m$^{-2}$)"
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
    colorbar.set_label(label)
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
    axes.set_title("Geometrically Necessary Dislocation Density")
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

    plt, _ = _matplotlib()
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

    plt, _ = _matplotlib()
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

    plt, _ = _matplotlib()
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
