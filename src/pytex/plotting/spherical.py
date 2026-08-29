from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
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
    if isinstance(items, Sequence) and not isinstance(items, str | bytes):
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
    if not isinstance(values, list | tuple) or not values:
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


def build_vector_stereogram_figure_spec(
    vectors: Any,
    *,
    labels: Sequence[str | None] | None = None,
    colors: Sequence[str] | None = None,
    method: str = "stereographic",
    render: str = "pole",
    antipodal: bool = True,
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    """Build a stereogram spec for arbitrary Cartesian direction vectors.

    The frame-agnostic projection primitive behind `plot_stereographic_vectors`:
    it plots any ``(n, 3)`` set of world-frame directions as poles and/or
    great-circle traces on a Wulff net, without requiring a `CrystalDirection`.
    Because the inputs are plain Cartesian vectors, directions from two crystals
    that were placed in a common frame with a `Transform3D` can be overlaid on
    one stereogram — the projection analog of a composite 3D scene.
    """

    if render not in {"pole", "trace", "both"}:
        raise ValueError("render must be one of 'pole', 'trace', or 'both'.")
    array = np.asarray(vectors, dtype=np.float64)
    if array.ndim == 1:
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("vectors must have shape (n, 3) or (3,).")
    norms = np.linalg.norm(array, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError("vectors must not contain a zero vector.")
    unit_vectors = array / norms[:, None]
    label_items = tuple(labels) if labels is not None else (None,) * unit_vectors.shape[0]
    if len(label_items) != unit_vectors.shape[0]:
        raise ValueError("labels must match the number of vectors.")
    common, spherical_style = _style_bundle(
        theme=theme,
        style_path=style_path,
        style_overrides=style_overrides,
    )
    palette = tuple(colors) if colors is not None else _palette(
        spherical_style,
        "direction_colors",
        ("#1f3a5f", "#bc6c25", "#4c956c", "#7c3aed"),
    )
    points = project_directions(unit_vectors, method=method, antipodal=antipodal)
    line_layers = list(
        _wulff_net_layers(method=method, spherical_style=spherical_style)
        if include_wulff_net
        else ()
    )
    if render in {"trace", "both"}:
        for index in range(unit_vectors.shape[0]):
            line_layers.append(
                LineLayer2D(
                    points=project_great_circle_trace(unit_vectors[index], method=method),
                    color=palette[index % len(palette)],
                    linewidth=float(spherical_style.get("plane_trace_linewidth", 1.45)),
                    alpha=float(spherical_style.get("plane_trace_alpha", 0.95)),
                )
            )
    marker_layers: tuple[MarkerLayer2D, ...] = ()
    if render in {"pole", "both"}:
        marker_layers = (
            MarkerLayer2D(
                points=points,
                marker=str(spherical_style.get("direction_marker", "o")),
                facecolors=[
                    palette[index % len(palette)] for index in range(unit_vectors.shape[0])
                ],
                edgecolors=str(spherical_style.get("direction_edgecolor", "#ffffff")),
                sizes=float(spherical_style.get("direction_size", 88.0)),
                linewidths=float(spherical_style.get("direction_linewidth", 1.0)),
                label="directions",
            ),
        )
    text_layers = []
    label_offset = float(spherical_style.get("label_offset", 0.04))
    for point, label in zip(points, label_items, strict=True):
        if label is None:
            continue
        text_layers.append(
            TextLayer2D(
                position=_radial_label_position(point, offset=label_offset),
                text=str(label),
                color=str(spherical_style.get("label_color", "#111111")),
                fontsize=float(spherical_style.get("label_fontsize", common["font"]["size"])),
                bbox_facecolor=str(spherical_style.get("label_bbox_color", "#ffffff")),
                bbox_alpha=float(spherical_style.get("label_bbox_alpha", 0.82)),
            )
        )
    radius = projection_boundary_radius(method)
    return FigureSpec2D(
        title=title or "Direction Stereogram",
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
        marker_layers=marker_layers,
        line_layers=tuple(line_layers),
        text_layers=tuple(text_layers),
    )


def plot_stereographic_vectors(
    vectors: Any,
    *,
    labels: Sequence[str | None] | None = None,
    colors: Sequence[str] | None = None,
    method: str = "stereographic",
    render: str = "pole",
    antipodal: bool = True,
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Plot arbitrary Cartesian direction vectors on a stereographic net.

    The general, frame-agnostic stereographic primitive: give it any ``(n, 3)``
    world-frame directions (a bare 3D vector, plane poles, or crystal directions
    already rotated into a common frame) and it renders them as poles and/or
    great-circle traces on a Wulff net. Overlay two crystals by concatenating
    their world-frame directions and passing per-vector ``colors`` / ``labels``.
    """

    return render_figure_spec(
        build_vector_stereogram_figure_spec(
            vectors,
            labels=labels,
            colors=colors,
            method=method,
            render=render,
            antipodal=antipodal,
            include_wulff_net=include_wulff_net,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
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


# --------------------------------------------------------------------------- #
# The orientation-relationship stereogram (F18)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORStereogramPair:
    """One OR-parallel (or near-parallel) pair placed on a common stereogram.

    Both vectors are Cartesian unit vectors **in the parent crystal frame**:
    the child-side object has already been carried into the parent frame by the
    variant rotation, which is what makes the pair comparable on one net. The
    child vector's sign is aligned to the parent's, so ``deviation_deg`` is the
    acute angle between the two poles and the tie-line drawn between them is
    the short arc.
    """

    kind: str
    label: str
    parent_vector: np.ndarray
    child_vector: np.ndarray
    deviation_deg: float

    def __post_init__(self) -> None:
        if self.kind not in {"plane", "direction"}:
            raise ValueError("ORStereogramPair.kind must be 'plane' or 'direction'.")
        vectors = []
        for name in ("parent_vector", "child_vector"):
            vector = np.asarray(getattr(self, name), dtype=np.float64).reshape(-1)
            if vector.shape != (3,):
                raise ValueError(f"{name} must have shape (3,).")
            norm = float(np.linalg.norm(vector))
            if not np.isfinite(norm) or norm <= 1e-12:
                raise ValueError(f"{name} must be a non-zero finite vector.")
            unit = np.ascontiguousarray(vector / norm)
            unit.setflags(write=False)
            vectors.append(unit)
        if not np.isfinite(self.deviation_deg) or self.deviation_deg < 0.0:
            raise ValueError("deviation_deg must be finite and non-negative.")
        object.__setattr__(self, "parent_vector", vectors[0])
        object.__setattr__(self, "child_vector", vectors[1])

    def __eq__(self, other: object) -> bool:
        # The generated __eq__ compares the vectors with ``==``, which yields an
        # array and raises "truth value is ambiguous" for every pair of equal
        # objects. Equality here is the pair's geometry and its label.
        if not isinstance(other, ORStereogramPair):
            return NotImplemented
        return (
            self.kind == other.kind
            and self.label == other.label
            and bool(np.array_equal(self.parent_vector, other.parent_vector))
            and bool(np.array_equal(self.child_vector, other.child_vector))
            and self.deviation_deg == other.deviation_deg
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.kind,
                self.label,
                self.parent_vector.tobytes(),
                self.child_vector.tobytes(),
                self.deviation_deg,
            )
        )


#: Below this, a pole component is numerical noise rather than geometry. It is
#: the scale at which a symmetry image of an equatorial pole comes back with
#: ``z = -8e-16`` instead of ``z = 0``.
_POLE_NOISE = 1e-9


def _fold_pair_upper(
    parent_unit: np.ndarray, child_unit: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Fold an already sign-aligned pair onto the upper hemisphere *together*.

    A stereogram folds antipodal directions onto one point, and the fold rule
    has to break ties for poles lying on the equator. Folding the two ends of a
    tie-line independently is then unsafe: a defining parallelism whose pole
    lies in the equatorial plane comes back from the variant rotation with
    ``z = -8e-16`` on one side and ``z = 0`` on the other, and the two ends of a
    zero-deviation tie-line land on **opposite rims** of the net — a diameter-long
    line across a figure whose whole claim is that the two poles coincide.

    The fold decision is therefore taken once, from the parent pole, and applied
    to both; the child's remaining sub-noise dip below the equator is then
    flattened onto it. A pair that genuinely straddles the equator (a real
    non-zero deviation) is left alone, and its tie-line is split by
    :func:`_tie_line_segments` as it should be.
    """

    below = float(parent_unit[2]) < -_POLE_NOISE
    on_equator = abs(float(parent_unit[2])) <= _POLE_NOISE
    equator_negative = on_equator and (
        float(parent_unit[0]) < -_POLE_NOISE
        or (abs(float(parent_unit[0])) <= _POLE_NOISE and float(parent_unit[1]) < 0.0)
    )
    if below or equator_negative:
        parent_unit, child_unit = -parent_unit, -child_unit
    if -_POLE_NOISE < float(child_unit[2]) < 0.0:
        child_unit = child_unit.copy()
        child_unit[2] = 0.0
        child_unit = child_unit / np.linalg.norm(child_unit)
    return parent_unit, child_unit


def _acute_pair(
    parent_vector: np.ndarray, child_vector: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Sign-align the child to the parent, fold the pair up, return the angle."""

    parent_unit = parent_vector / np.linalg.norm(parent_vector)
    child_unit = child_vector / np.linalg.norm(child_vector)
    if float(parent_unit @ child_unit) < 0.0:
        child_unit = -child_unit
    parent_unit, child_unit = _fold_pair_upper(parent_unit, child_unit)
    # arccos of a dot product loses half the significant digits exactly where
    # these angles live -- near zero, where cos is flat -- and reports 1e-6 deg
    # for a parallelism that is exact to machine precision. The half-angle
    # tangent form stays conditioned there.
    angle = 2.0 * np.arctan2(
        float(np.linalg.norm(parent_unit - child_unit)),
        float(np.linalg.norm(parent_unit + child_unit)),
    )
    return parent_unit, child_unit, float(np.rad2deg(angle))


def _canonical_index_sign(indices: tuple[int, ...]) -> tuple[int, ...]:
    """Flip an index triple so its first non-zero component is positive.

    A plane has no sign — ``(111)`` and ``(-1-1-1)`` name the same plane, and a
    stereogram folds them onto one point. Without this, the 24 Kurdjumov-Sachs
    variants appear to name eight distinct parent planes where they name four.
    """

    for value in indices:
        if value > 0:
            return indices
        if value < 0:
            return tuple(-value for value in indices)
    return indices


def _or_pair_label(parent_indices: Any, child_indices: Any, *, kind: str) -> str:
    formatter = format_plane_indices if kind == "plane" else format_direction_indices
    parent = tuple(round(float(value)) for value in np.asarray(parent_indices).reshape(-1))
    child = tuple(round(float(value)) for value in np.asarray(child_indices).reshape(-1))
    if kind == "plane":
        parent = _canonical_index_sign(parent)
        child = _canonical_index_sign(child)
    return f"{formatter(parent, style='plain')} ∥ {formatter(child, style='plain')}"


def _defining_stereogram_pairs(variant: Any, rotation: np.ndarray) -> list[ORStereogramPair]:
    pairs: list[ORStereogramPair] = []
    for parent_plane, child_plane in variant.parallel_planes:
        parent_unit, child_unit, deviation = _acute_pair(
            parent_plane.normal, rotation.T @ child_plane.normal
        )
        pairs.append(
            ORStereogramPair(
                kind="plane",
                label=_or_pair_label(
                    parent_plane.miller.indices, child_plane.miller.indices, kind="plane"
                ),
                parent_vector=parent_unit,
                child_vector=child_unit,
                deviation_deg=deviation,
            )
        )
    for parent_direction, child_direction in variant.parallel_directions:
        parent_unit, child_unit, deviation = _acute_pair(
            parent_direction.unit_vector, rotation.T @ child_direction.unit_vector
        )
        pairs.append(
            ORStereogramPair(
                kind="direction",
                label=_or_pair_label(
                    parent_direction.coordinates,
                    child_direction.coordinates,
                    kind="direction",
                ),
                parent_vector=parent_unit,
                child_vector=child_unit,
                deviation_deg=deviation,
            )
        )
    return pairs


def _nominated_stereogram_pairs(
    relationship: Any,
    variant: Any,
    rotation: np.ndarray,
    *,
    planes: Sequence[Any],
    directions: Sequence[Any],
    tolerance_deg: float,
) -> list[ORStereogramPair]:
    from pytex.core.lattice import MillerIndex
    from pytex.core.transformation import find_parallel_directions, find_parallel_planes

    parent_phase = relationship.parent_phase
    child_phase = relationship.child_phase
    pairs: list[ORStereogramPair] = []
    for kind, nominated in (("plane", planes), ("direction", directions)):
        for item in nominated:
            report = (
                find_parallel_planes(
                    relationship, item, tolerance_deg=tolerance_deg, variants=(variant,)
                )
                if kind == "plane"
                else find_parallel_directions(
                    relationship, item, tolerance_deg=tolerance_deg, variants=(variant,)
                )
            )
            for match in report.matches:
                if kind == "plane":
                    parent_cartesian = CrystalPlane(
                        MillerIndex(match.parent_indices, phase=parent_phase),
                        phase=parent_phase,
                    ).normal
                    child_cartesian = CrystalPlane(
                        MillerIndex(match.child_indices, phase=child_phase),
                        phase=child_phase,
                    ).normal
                else:
                    parent_cartesian = CrystalDirection(
                        match.parent_indices.astype(np.float64), phase=parent_phase
                    ).unit_vector
                    child_cartesian = CrystalDirection(
                        match.child_indices.astype(np.float64), phase=child_phase
                    ).unit_vector
                parent_unit, child_unit, _ = _acute_pair(
                    parent_cartesian, rotation.T @ child_cartesian
                )
                pairs.append(
                    ORStereogramPair(
                        kind=kind,
                        label=_or_pair_label(
                            match.parent_indices, match.child_indices, kind=kind
                        ),
                        parent_vector=parent_unit,
                        child_vector=child_unit,
                        # the report's own deviation, not a re-derived one:
                        # one number, one owner
                        deviation_deg=float(match.angular_deviation_deg),
                    )
                )
    return pairs


def or_stereogram_pairs(
    relationship: Any,
    *,
    variant: int | Any | None = None,
    parent_planes: Sequence[Any] = (),
    parent_directions: Sequence[Any] = (),
    tolerance_deg: float = 3.0,
) -> tuple[ORStereogramPair, ...]:
    """The OR-parallel pairs of one variant, both sides in the parent frame.

    Purpose
    -------
    The data behind the OR stereogram: what is parallel to what, where each
    object lands on a net drawn in the parent crystal frame, and how exactly
    the parallelism holds.

    Parameters
    ----------
    relationship : OrientationRelationship
    variant : int or TransformationVariant, optional
        One-based index or a variant object. ``None`` means variant 1, which
        *is* the relationship as stated: its parent symmetry operator is the
        identity, so its rotation and its parallelisms are the nominal ones.
    parent_planes, parent_directions : sequence of CrystalPlane / CrystalDirection
        Extra parent families to nominate. Their child partners are found with
        ``find_parallel_planes`` / ``find_parallel_directions`` within
        ``tolerance_deg``, and those functions own the deviation each pair is
        labelled with. Read that number precisely: the *exact* child image of
        any parent plane is parallel to it by construction, so what the search
        reports — and what this figure draws and labels — is the
        **rationalization residual**, the angle by which the nearest low-index
        child index misses the exact image. A small tolerance therefore keeps
        the pairs for which a low-index child object really is parallel, and
        drops the parent members for which none is.
    tolerance_deg : float
        Angular tolerance for the nominated families, in the sense just stated.

    Returns
    -------
    tuple of ORStereogramPair
        The variant's own defining pairs first, then any nominated matches.

    Notes
    -----
    The defining pairs come from ``TransformationVariant.parallel_planes`` and
    ``.parallel_directions``, so they are *this variant's* symmetry images and
    not the relationship's nominal pair. Their deviation is zero by
    construction; it is still reported rather than asserted, because a figure
    that prints the number it claims is worth more than one that asserts
    perfection.
    """

    from pytex.plotting.scene3d import resolve_transformation_variant

    resolved = resolve_transformation_variant(relationship, variant)
    if resolved is None:
        resolved = relationship.generate_variants()[0]
    rotation = np.asarray(resolved.parent_to_child_rotation.as_matrix(), dtype=np.float64)
    pairs = _defining_stereogram_pairs(resolved, rotation)
    pairs.extend(
        _nominated_stereogram_pairs(
            relationship,
            resolved,
            rotation,
            planes=tuple(parent_planes),
            directions=tuple(parent_directions),
            tolerance_deg=tolerance_deg,
        )
    )
    # A nominated family contains the defining member, so the defining pair
    # comes back a second time; drawn twice it is two tie-lines and two labels
    # on one another. Deduplicate on the geometry actually plotted.
    seen: set[tuple[str, tuple[float, ...], tuple[float, ...]]] = set()
    unique: list[ORStereogramPair] = []
    for pair in pairs:
        key = (
            pair.kind,
            tuple(np.round(pair.parent_vector, 9) + 0.0),
            tuple(np.round(pair.child_vector, 9) + 0.0),
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(pair)
    return tuple(unique)


def _slerp(start: np.ndarray, end: np.ndarray, *, samples: int) -> np.ndarray:
    """Great-circle interpolation between two unit vectors, endpoints included."""

    cosine = float(np.clip(start @ end, -1.0, 1.0))
    angle = float(np.arccos(cosine))
    if angle < 1e-9:
        return np.repeat(start[None, :], samples, axis=0)
    fractions = np.linspace(0.0, 1.0, samples)[:, None]
    sin_angle = np.sin(angle)
    weights_start = np.sin((1.0 - fractions) * angle) / sin_angle
    weights_end = np.sin(fractions * angle) / sin_angle
    arc = weights_start * start[None, :] + weights_end * end[None, :]
    return np.asarray(arc, dtype=np.float64)


def _tie_line_segments(
    start: np.ndarray,
    end: np.ndarray,
    *,
    method: str,
    samples: int = 33,
) -> list[np.ndarray]:
    """Projected tie-line polylines, split where antipodal folding jumps.

    A pair straddling the equator folds onto opposite sides of the disc, and a
    single polyline drawn through that fold would cross the whole net along a
    chord no crystallography justifies. The polyline is therefore cut wherever
    consecutive projected points jump by more than the disc radius.
    """

    arc = _slerp(start, end, samples=samples)
    projected = np.asarray(project_directions(arc, method=method, antipodal=True))
    return _split_on_fold_jumps(projected, method=method)


def _split_on_fold_jumps(points: np.ndarray, *, method: str) -> list[np.ndarray]:
    """Cut a projected polyline wherever the antipodal fold jumps the rim.

    A curve that leaves the projected hemisphere re-enters at the antipodal rim
    point. Drawn as one polyline it is a chord straight across the net that no
    crystallography justifies, so the polyline is cut wherever consecutive
    points jump by more than the disc radius.
    """

    if points.shape[0] < 2:
        return []
    jumps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    breaks = np.flatnonzero(jumps > projection_boundary_radius(method)) + 1
    return [segment for segment in np.split(points, breaks) if segment.shape[0] >= 2]


def build_or_stereogram_figure_spec(
    relationship: Any,
    *,
    variant: int | Any | None = None,
    parent_planes: Sequence[Any] = (),
    parent_directions: Sequence[Any] = (),
    tolerance_deg: float = 3.0,
    show_great_circles: bool = True,
    show_tie_lines: bool = True,
    label_pairs: bool = True,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> FigureSpec2D:
    """Build the orientation-relationship stereogram: the parallelism, drawn.

    Purpose
    -------
    The figure by which orientation relationships are read in the literature,
    and the natural teaching object for a parallelism statement. One net, drawn
    in the **parent crystal frame**, carrying:

    - the parent pole of each pair as an **open** symbol and the child pole,
      carried into the parent frame by the variant rotation, as a **filled**
      one, so a parallelism reads as two symbols on top of each other;
    - a **tie-line** joining each pair along its great circle, labelled with
      the deviation in degrees;
    - for plane pairs, the **great circles** of both planes, so a plane
      parallelism reads as two coincident circles rather than two coincident
      points -- the parent circle dashed, the child solid.

    Parameters
    ----------
    relationship : OrientationRelationship
    variant : int or TransformationVariant, optional
        One-based index or a variant object; ``None`` means variant 1, which is
        the relationship as stated. Each variant carries its **own** parallel
        pair (see ``TransformationVariant.parallel_planes``), so the labels move
        with the variant instead of repeating variant 1's indices.
    parent_planes, parent_directions : sequence of CrystalPlane / CrystalDirection
        Extra parent families to nominate; their partners within
        ``tolerance_deg`` are added as further tie-lines, and this is where the
        deviation labels stop being zero. The number is the rationalization
        residual of the child index, not a departure from parallelism — see
        `or_stereogram_pairs`, which states it in full. The tie-line drawn is
        exactly that gap, so the figure and its label agree.
    tolerance_deg : float
    show_great_circles, show_tie_lines, label_pairs : bool
    method : str
        ``"stereographic"`` (default, conformal -- the right net for reading
        angles) or ``"equal_area"``.
    include_wulff_net : bool
    title : str, optional
    theme, style_path, style_overrides
        Standard style resolution.

    Returns
    -------
    FigureSpec2D
        Render it with `render_figure_spec`, or call `plot_or_stereogram`.

    Raises
    ------
    ValueError
        If the relationship states no parallelisms and none were nominated,
        leaving nothing to draw.

    See Also
    --------
    or_stereogram_pairs : the pairs and deviations this figure draws.
    """

    pairs = or_stereogram_pairs(
        relationship,
        variant=variant,
        parent_planes=parent_planes,
        parent_directions=parent_directions,
        tolerance_deg=tolerance_deg,
    )
    if not pairs:
        raise ValueError(
            "The orientation relationship states no parallel planes or directions "
            "and none were nominated; there is nothing to draw."
        )
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
    line_layers = list(
        _wulff_net_layers(method=method, spherical_style=spherical_style)
        if include_wulff_net
        else ()
    )
    parent_points = project_directions(
        np.stack([pair.parent_vector for pair in pairs]), method=method, antipodal=True
    )
    child_points = project_directions(
        np.stack([pair.child_vector for pair in pairs]), method=method, antipodal=True
    )
    colors = [palette[index % len(palette)] for index in range(len(pairs))]
    trace_linewidth = float(spherical_style.get("plane_trace_linewidth", 1.45))
    for index, pair in enumerate(pairs):
        color = colors[index]
        if show_great_circles and pair.kind == "plane":
            # The child circle is drawn first as a wide pale halo and the
            # parent's dashes over it. Two coincident circles in the same
            # colour and weight are indistinguishable from one -- the figure
            # would show one circle where it means two -- so the difference is
            # carried by weight and opacity, not only by dash pattern.
            for pole, linestyle, alpha, width in (
                (pair.child_vector, "-", 0.35, trace_linewidth + 2.2),
                (pair.parent_vector, "--", 1.0, trace_linewidth),
            ):
                trace = project_great_circle_trace(pole, method=method)
                for segment in _split_on_fold_jumps(np.asarray(trace), method=method):
                    line_layers.append(
                        LineLayer2D(
                            points=segment,
                            color=color,
                            linewidth=width,
                            linestyle=linestyle,
                            alpha=alpha,
                        )
                    )
        if show_tie_lines:
            for segment in _tie_line_segments(
                pair.parent_vector, pair.child_vector, method=method
            ):
                line_layers.append(
                    LineLayer2D(
                        points=segment,
                        color=color,
                        linewidth=trace_linewidth + 0.9,
                        linestyle="-",
                    )
                )
    marker_size = float(spherical_style.get("direction_size", 88.0))
    marker_layers = (
        MarkerLayer2D(
            points=parent_points,
            marker="o",
            facecolors=["none"] * len(pairs),
            edgecolors=colors,
            sizes=marker_size * 1.6,
            linewidths=1.6,
            label="parent",
        ),
        MarkerLayer2D(
            points=child_points,
            marker="o",
            facecolors=colors,
            edgecolors=str(spherical_style.get("direction_edgecolor", "#ffffff")),
            sizes=marker_size * 0.8,
            linewidths=float(spherical_style.get("direction_linewidth", 1.0)),
            label="child",
        ),
    )
    text_layers: list[TextLayer2D] = []
    if label_pairs:
        label_offset = float(spherical_style.get("label_offset", 0.04))
        for index, pair in enumerate(pairs):
            midpoint = 0.5 * (parent_points[index] + child_points[index])
            text_layers.append(
                TextLayer2D(
                    position=_radial_label_position(midpoint, offset=2.0 * label_offset),
                    text=f"{pair.label}  {pair.deviation_deg:.2f} deg",
                    color=colors[index],
                    fontsize=float(
                        spherical_style.get("label_fontsize", common["font"]["size"])
                    ),
                    bbox_facecolor=str(spherical_style.get("label_bbox_color", "#ffffff")),
                    bbox_alpha=float(spherical_style.get("label_bbox_alpha", 0.82)),
                )
            )
    radius = projection_boundary_radius(method)
    return FigureSpec2D(
        title=title or f"{relationship.name}: OR stereogram (parent frame)",
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
        marker_layers=marker_layers,
        line_layers=tuple(line_layers),
        text_layers=tuple(text_layers),
    )


def plot_or_stereogram(
    relationship: Any,
    *,
    variant: int | Any | None = None,
    parent_planes: Sequence[Any] = (),
    parent_directions: Sequence[Any] = (),
    tolerance_deg: float = 3.0,
    show_great_circles: bool = True,
    show_tie_lines: bool = True,
    label_pairs: bool = True,
    method: str = "stereographic",
    include_wulff_net: bool = True,
    title: str | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    """Render the orientation-relationship stereogram.

    See `build_or_stereogram_figure_spec` for what the figure contains and why;
    this is the one-call form, returning the Matplotlib figure.
    """

    return render_figure_spec(
        build_or_stereogram_figure_spec(
            relationship,
            variant=variant,
            parent_planes=parent_planes,
            parent_directions=parent_directions,
            tolerance_deg=tolerance_deg,
            show_great_circles=show_great_circles,
            show_tie_lines=show_tie_lines,
            label_pairs=label_pairs,
            method=method,
            include_wulff_net=include_wulff_net,
            title=title,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        ),
        ax=ax,
    )


__all__ = [
    "ORStereogramPair",
    "build_crystal_direction_figure_spec",
    "build_crystal_plane_figure_spec",
    "build_or_stereogram_figure_spec",
    "build_symmetry_elements_figure_spec",
    "build_vector_stereogram_figure_spec",
    "build_wulff_net_figure_spec",
    "or_stereogram_pairs",
    "plot_crystal_directions",
    "plot_crystal_planes",
    "plot_or_stereogram",
    "plot_stereographic_vectors",
    "plot_symmetry_elements",
    "plot_wulff_net",
]
