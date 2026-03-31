from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from pytex.core.lattice import CrystalDirection, CrystalPlane
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.stereonets import (
    generate_stereonet_grid,
    project_great_circle_trace,
    projection_boundary_radius,
)
from pytex.plotting._render import (
    FigureSpec2D,
    LineLayer2D,
    MarkerLayer2D,
    TextLayer2D,
    render_figure_spec,
)
from pytex.plotting.styles import resolve_style
from pytex.texture.projections import project_directions

_SYMMETRY_MARKERS = {
    2: "D",
    3: "^",
    4: "s",
    6: "h",
}


def _as_tuple(items: Any) -> tuple[Any, ...]:
    if isinstance(items, Sequence) and not isinstance(items, (str, bytes)):
        return tuple(items)
    return (items,)


def _style_bundle(
    *,
    theme: str,
    style_path: str | None,
    style_overrides: dict[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    return style["common"], style["spherical"]


def _palette(section: dict[str, Any], key: str, fallback: tuple[str, ...]) -> tuple[str, ...]:
    values = section.get(key, fallback)
    if not isinstance(values, (list, tuple)) or not values:
        return fallback
    return tuple(str(value) for value in values)


def _radial_label_position(point: np.ndarray, *, offset: float) -> np.ndarray:
    norm = float(np.linalg.norm(point))
    if norm <= 1e-9:
        direction = np.array([1.0, 1.0], dtype=np.float64) / np.sqrt(2.0)
    else:
        direction = point / norm
    return np.asarray(point + offset * direction, dtype=np.float64)


def _axis_to_indices(axis: np.ndarray, *, max_index: int = 6) -> tuple[int, int, int] | None:
    abs_axis = np.abs(axis)
    nonzero = abs_axis > 1e-8
    if not np.any(nonzero):
        return None
    scale = 1.0 / float(np.min(abs_axis[nonzero]))
    scaled = np.rint(axis * scale).astype(np.int64)
    gcd_value = int(np.gcd.reduce(np.abs(scaled[np.abs(scaled) > 0])))
    gcd_value = max(gcd_value, 1)
    scaled //= gcd_value
    if np.max(np.abs(scaled)) > max_index:
        return None
    normalized = scaled / np.linalg.norm(scaled)
    if not np.allclose(normalized, axis / np.linalg.norm(axis), atol=1e-6):
        return None
    return (int(scaled[0]), int(scaled[1]), int(scaled[2]))


def _direction_label(direction: CrystalDirection, label: str | Sequence[int] | None) -> str | None:
    if isinstance(label, str):
        return label
    if label is not None:
        return format_direction_indices(tuple(int(value) for value in label))
    rounded = np.rint(direction.coordinates).astype(np.int64)
    if np.allclose(rounded.astype(np.float64), direction.coordinates, atol=1e-8):
        return format_direction_indices(tuple(int(value) for value in rounded))
    return None


def _plane_label(plane: CrystalPlane, label: str | Sequence[int] | None) -> str | None:
    if isinstance(label, str):
        return label
    if label is not None:
        return format_plane_indices(tuple(int(value) for value in label))
    return format_plane_indices(tuple(int(value) for value in plane.miller.indices))


def _wulff_net_layers(
    *,
    method: str,
    spherical_style: dict[str, Any],
) -> tuple[LineLayer2D, ...]:
    minor_step = None
    if spherical_style.get("show_minor_grid", True):
        minor_step = float(spherical_style.get("minor_step_deg", 2.0))
    grid = generate_stereonet_grid(
        method=method,
        major_step_deg=float(spherical_style.get("major_step_deg", 10.0)),
        minor_step_deg=minor_step,
    )
    minor_layers = tuple(
        LineLayer2D(
            points=line,
            color=str(spherical_style.get("net_minor_color", "#cbd5e1")),
            linewidth=float(spherical_style.get("net_minor_linewidth", 0.45)),
            alpha=float(spherical_style.get("net_alpha", 0.85)),
        )
        for line in grid.minor_lines
    )
    major_layers = tuple(
        LineLayer2D(
            points=line,
            color=str(spherical_style.get("net_major_color", "#94a3b8")),
            linewidth=float(spherical_style.get("net_major_linewidth", 0.8)),
            alpha=float(spherical_style.get("net_alpha", 0.9)),
        )
        for line in grid.major_lines
    )
    return minor_layers + major_layers


def build_wulff_net_figure_spec(
    *,
    method: str = "stereographic",
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    _, spherical_style = _style_bundle(
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    return FigureSpec2D(
        title=title or "Wulff Net",
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-projection_boundary_radius(method), projection_boundary_radius(method)),
        ylim=(-projection_boundary_radius(method), projection_boundary_radius(method)),
        boundary_circle_radius=projection_boundary_radius(method),
        boundary_circle_color=str(spherical_style.get("boundary_color", "#0f172a")),
        boundary_circle_linewidth=float(spherical_style.get("boundary_linewidth", 1.15)),
        boundary_circle_linestyle="-",
        equal_aspect=True,
        grid=False,
        show_axes=False,
        line_layers=_wulff_net_layers(method=method, spherical_style=spherical_style),
    )


def build_crystal_direction_figure_spec(
    directions: CrystalDirection | Sequence[CrystalDirection],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    direction_items = _as_tuple(directions)
    label_items = tuple(labels) if labels is not None else (None,) * len(direction_items)
    if len(label_items) != len(direction_items):
        raise ValueError("labels must match the number of crystal directions.")
    common, spherical_style = _style_bundle(
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    palette = _palette(
        spherical_style,
        "direction_colors",
        ("#1f3a5f", "#bc6c25", "#4c956c", "#7c3aed"),
    )
    direction_vectors = np.stack([direction.unit_vector for direction in direction_items], axis=0)
    points = project_directions(direction_vectors, method=method, antipodal=True)
    text_layers = []
    label_offset = float(spherical_style.get("label_offset", 0.04))
    for point, direction, label in zip(points, direction_items, label_items, strict=True):
        formatted = _direction_label(direction, label)
        if formatted is None:
            continue
        text_layers.append(
            TextLayer2D(
                position=_radial_label_position(point, offset=label_offset),
                text=formatted,
                color=str(spherical_style.get("label_color", "#111111")),
                fontsize=float(spherical_style.get("label_fontsize", common["font"]["size"])),
                bbox_facecolor=str(spherical_style.get("label_bbox_color", "#ffffff")),
                bbox_alpha=float(spherical_style.get("label_bbox_alpha", 0.82)),
            )
        )
    marker_layer = MarkerLayer2D(
        points=points,
        marker=str(spherical_style.get("direction_marker", "o")),
        facecolors=[palette[index % len(palette)] for index in range(len(direction_items))],
        edgecolors=str(spherical_style.get("direction_edgecolor", "#ffffff")),
        sizes=float(spherical_style.get("direction_size", 88.0)),
        linewidths=float(spherical_style.get("direction_linewidth", 1.0)),
        label="directions",
    )
    line_layers = (
        _wulff_net_layers(method=method, spherical_style=spherical_style)
        if include_wulff_net
        else ()
    )
    radius = projection_boundary_radius(method)
    return FigureSpec2D(
        title=title or "Crystal Directions",
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        boundary_circle_radius=radius,
        boundary_circle_color=str(spherical_style.get("boundary_color", "#0f172a")),
        boundary_circle_linewidth=float(spherical_style.get("boundary_linewidth", 1.15)),
        boundary_circle_linestyle="-",
        equal_aspect=True,
        grid=False,
        show_axes=False,
        marker_layers=(marker_layer,),
        line_layers=line_layers,
        text_layers=tuple(text_layers),
    )


def build_crystal_plane_figure_spec(
    planes: CrystalPlane | Sequence[CrystalPlane],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    render: str = "trace",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    if render not in {"pole", "trace", "both"}:
        raise ValueError("render must be one of 'pole', 'trace', or 'both'.")
    plane_items = _as_tuple(planes)
    label_items = tuple(labels) if labels is not None else (None,) * len(plane_items)
    if len(label_items) != len(plane_items):
        raise ValueError("labels must match the number of crystal planes.")
    common, spherical_style = _style_bundle(
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    palette = _palette(
        spherical_style,
        "plane_colors",
        ("#2563eb", "#dc2626", "#0f766e", "#b45309"),
    )
    line_layers = list(
        _wulff_net_layers(method=method, spherical_style=spherical_style)
        if include_wulff_net
        else ()
    )
    marker_layers: list[MarkerLayer2D] = []
    text_layers: list[TextLayer2D] = []
    label_offset = float(spherical_style.get("label_offset", 0.04))
    poles = np.stack([plane.normal for plane in plane_items], axis=0)
    pole_points = project_directions(poles, method=method, antipodal=True)
    if render in {"trace", "both"}:
        for index, plane in enumerate(plane_items):
            line_layers.append(
                LineLayer2D(
                    points=project_great_circle_trace(plane.normal, method=method),
                    color=palette[index % len(palette)],
                    linewidth=float(spherical_style.get("plane_trace_linewidth", 1.45)),
                    alpha=float(spherical_style.get("plane_trace_alpha", 0.95)),
                    label=None,
                )
            )
    if render in {"pole", "both"}:
        marker_layers.append(
            MarkerLayer2D(
                points=pole_points,
                marker=str(spherical_style.get("plane_pole_marker", "o")),
                facecolors=str(spherical_style.get("plane_pole_facecolor", "#ffffff")),
                edgecolors=[palette[index % len(palette)] for index in range(len(plane_items))],
                sizes=float(spherical_style.get("plane_pole_size", 72.0)),
                linewidths=float(spherical_style.get("plane_pole_linewidth", 1.3)),
                label="plane poles",
            )
        )
    for point, plane, label in zip(pole_points, plane_items, label_items, strict=True):
        formatted = _plane_label(plane, label)
        if formatted is None:
            continue
        text_layers.append(
            TextLayer2D(
                position=_radial_label_position(point, offset=label_offset),
                text=formatted,
                color=str(spherical_style.get("label_color", "#111111")),
                fontsize=float(spherical_style.get("label_fontsize", common["font"]["size"])),
                bbox_facecolor=str(spherical_style.get("label_bbox_color", "#ffffff")),
                bbox_alpha=float(spherical_style.get("label_bbox_alpha", 0.82)),
            )
        )
    radius = projection_boundary_radius(method)
    return FigureSpec2D(
        title=title or "Crystal Planes",
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        boundary_circle_radius=radius,
        boundary_circle_color=str(spherical_style.get("boundary_color", "#0f172a")),
        boundary_circle_linewidth=float(spherical_style.get("boundary_linewidth", 1.15)),
        boundary_circle_linestyle="-",
        equal_aspect=True,
        grid=False,
        show_axes=False,
        marker_layers=tuple(marker_layers),
        line_layers=tuple(line_layers),
        text_layers=tuple(text_layers),
    )


def _canonical_axis(axis: np.ndarray) -> np.ndarray:
    canonical = np.array(axis, copy=True)
    nonzero = np.flatnonzero(np.abs(canonical) > 1e-10)
    if nonzero.size > 0 and canonical[int(nonzero[0])] < 0.0:
        canonical *= -1.0
    canonical = np.asarray(canonical / np.linalg.norm(canonical), dtype=np.float64)
    return canonical


def _symmetry_axes_by_order(symmetry: SymmetrySpec) -> dict[int, tuple[np.ndarray, ...]]:
    order_by_axis: dict[tuple[float, float, float], int] = {}
    vector_by_axis: dict[tuple[float, float, float], np.ndarray] = {}
    for operator in symmetry.operators:
        rotation = Rotation.from_matrix(operator)
        if np.isclose(rotation.angle_deg, 0.0, atol=1e-8):
            continue
        axis = _canonical_axis(rotation.axis)
        order = max(2, round(360.0 / rotation.angle_deg))
        key = (
            float(np.round(axis[0], 8)),
            float(np.round(axis[1], 8)),
            float(np.round(axis[2], 8)),
        )
        if key not in order_by_axis or order > order_by_axis[key]:
            order_by_axis[key] = order
            vector_by_axis[key] = axis
    grouped: dict[int, list[np.ndarray]] = {}
    for key, order in order_by_axis.items():
        grouped.setdefault(order, []).append(vector_by_axis[key])
    return {order: tuple(vectors) for order, vectors in grouped.items()}


def build_symmetry_elements_figure_spec(
    symmetry: SymmetrySpec,
    *,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    annotate_axes: bool = False,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    common, spherical_style = _style_bundle(
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    line_layers = (
        _wulff_net_layers(method=method, spherical_style=spherical_style)
        if include_wulff_net
        else ()
    )
    marker_layers: list[MarkerLayer2D] = []
    text_layers: list[TextLayer2D] = []
    grouped = _symmetry_axes_by_order(symmetry)
    color_map = dict(spherical_style.get("symmetry_colors", {}))
    size_map = dict(spherical_style.get("symmetry_size", {}))
    label_offset = float(spherical_style.get("label_offset", 0.04))
    for order in sorted(grouped):
        axes = np.stack(grouped[order], axis=0)
        points = project_directions(axes, method=method, antipodal=True)
        marker_layers.append(
            MarkerLayer2D(
                points=points,
                marker=_SYMMETRY_MARKERS.get(order, "o"),
                facecolors=str(color_map.get(str(order), "#ffffff")),
                edgecolors=str(spherical_style.get("symmetry_edgecolor", "#0f172a")),
                sizes=float(size_map.get(str(order), 96.0 + 10.0 * order)),
                linewidths=float(spherical_style.get("symmetry_linewidth", 1.1)),
                label=f"{order}-fold axes",
            )
        )
        if annotate_axes:
            for axis, point in zip(axes, points, strict=True):
                indices = _axis_to_indices(axis)
                if indices is None:
                    continue
                text_layers.append(
                    TextLayer2D(
                        position=_radial_label_position(point, offset=label_offset),
                        text=format_direction_indices(indices),
                        color=str(spherical_style.get("label_color", "#111111")),
                        fontsize=float(
                            spherical_style.get("label_fontsize", common["font"]["size"])
                        ),
                        bbox_facecolor=str(spherical_style.get("label_bbox_color", "#ffffff")),
                        bbox_alpha=float(spherical_style.get("label_bbox_alpha", 0.82)),
                    )
                )
    radius = projection_boundary_radius(method)
    return FigureSpec2D(
        title=title or f"Symmetry Elements: {symmetry.point_group}",
        xlabel="projection x",
        ylabel="projection y",
        xlim=(-radius, radius),
        ylim=(-radius, radius),
        boundary_circle_radius=radius,
        boundary_circle_color=str(spherical_style.get("boundary_color", "#0f172a")),
        boundary_circle_linewidth=float(spherical_style.get("boundary_linewidth", 1.15)),
        boundary_circle_linestyle="-",
        equal_aspect=True,
        grid=False,
        show_axes=False,
        marker_layers=tuple(marker_layers),
        line_layers=line_layers,
        text_layers=tuple(text_layers),
    )


def plot_wulff_net(
    *,
    method: str = "stereographic",
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_wulff_net_figure_spec(
            method=method,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
    )


def plot_crystal_directions(
    directions: CrystalDirection | Sequence[CrystalDirection],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_crystal_direction_figure_spec(
            directions,
            labels=labels,
            method=method,
            include_wulff_net=include_wulff_net,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
    )


def plot_crystal_planes(
    planes: CrystalPlane | Sequence[CrystalPlane],
    *,
    labels: Sequence[str | Sequence[int] | None] | None = None,
    method: str = "stereographic",
    render: str = "trace",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_crystal_plane_figure_spec(
            planes,
            labels=labels,
            method=method,
            render=render,
            include_wulff_net=include_wulff_net,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
    )


def plot_symmetry_elements(
    symmetry: SymmetrySpec,
    *,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    annotate_axes: bool = False,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    return render_figure_spec(
        build_symmetry_elements_figure_spec(
            symmetry,
            method=method,
            include_wulff_net=include_wulff_net,
            annotate_axes=annotate_axes,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
    )


__all__ = [
    "build_crystal_direction_figure_spec",
    "build_crystal_plane_figure_spec",
    "build_symmetry_elements_figure_spec",
    "build_wulff_net_figure_spec",
    "plot_crystal_directions",
    "plot_crystal_planes",
    "plot_symmetry_elements",
    "plot_wulff_net",
]
