"""Renderer-independent 3D geometric primitives and their composition.

This module is the *core vocabulary* every composite PyTex visualization is
built from. A primitive is an immutable, backend-agnostic description of one
geometric element in a shared **world Cartesian frame** (angstrom):

- `Arrow3D` — a 3D vector (the atom of the whole system: a bare vector, a
  Miller direction, a plane normal, a triad axis are all arrows).
- `PolyLine3D` — an open or closed polyline (cell edges, projected traces).
- `PlanePatch3D` — a planar polygon with an outward normal (a Miller-plane
  patch, a slip plane, a habit plane).
- `PointCloud3D` — a set of point markers (lattice nodes, atom sites, poles).
- `Label3D` — a free text label anchored at a world point.
- `AxisTriad3D` — an orthonormal axis triad (a reference-frame gizmo:
  specimen RD/TD/ND, crystal a/b/c), expanded to three arrows plus labels.

`Transform3D` places any primitive (or a whole `PrimitiveScene3D`) in the world
frame from a `Rotation`/`Orientation` plus a translation, so the same primitive
description can be reused across crystals held in different orientations — the
mechanism that lets *two crystals in an orientation relationship* be assembled
with a handful of primitives.

Two layers are kept strictly separate so the model stays portable:

1. **Model** — the dataclasses above and `PrimitiveScene3D`, carrying only
   geometry, colors, and annotation intent (nothing matplotlib-specific).
2. **Renderer** — `render_primitive_scene_3d` (and the shared `_draw_primitive_scene`
   helper reused by the composite world-scene renderer) rasterizes a scene with
   matplotlib.

Builders (`vector_arrow`, `direction_arrow`, `plane_normal_arrow`,
`crystal_plane_patch`, `reference_frame_triad`, `unit_cell_polylines`,
`lattice_point_cloud`) turn the canonical crystallographic objects
(`CrystalDirection`, `CrystalPlane`, `ReferenceFrame`, `Phase`) into primitives
so the crystallographic surface and the drawing surface share one language.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, is_rotation_matrix
from pytex.core.lattice import CrystalDirection, CrystalPlane, Lattice, Phase
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.plotting.styles import resolve_style

# Okabe-Ito derived triad palette (a=blue, b=green, c=vermillion): chosen so
# axis identity survives grayscale printing and the common color-vision
# deficiencies, matching the repository categorical palette.
_DEFAULT_TRIAD_COLORS: tuple[str, str, str] = ("#1d4ed8", "#059669", "#dc2626")

#: The canonical per-axis triad colors, in first/second/third axis order.
#:
#: Exported so every renderer that draws a reference frame — 3D scene triads,
#: the 2D corner gizmo in `pytex.plotting.frames`, and documentation SVG — uses
#: one palette, and a frame therefore looks the same wherever it appears.
TRIAD_AXIS_COLORS: tuple[str, str, str] = _DEFAULT_TRIAD_COLORS


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:  # pragma: no cover - environment-dependent branch
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, Poly3DCollection


@dataclass(frozen=True, slots=True)
class Transform3D:
    """A rigid-or-affine placement of geometry into the shared world frame.

    Applies as ``x_world = matrix @ x_local + translation`` (equivalently
    ``x_local @ matrix.T + translation`` for row-stacked point arrays). The
    ``matrix`` is the linear part (a proper rotation for crystal placement, any
    invertible linear map for general primitives); ``translation`` shifts the
    placed geometry so several crystals can sit side by side in one scene.

    When to use it: to hold one crystal (or primitive group) at an
    `Orientation`/`Rotation` relative to another — for example placing the
    child crystal of an orientation relationship so its parallel planes and
    directions coincide with the parent's.
    """

    matrix: np.ndarray = field(default_factory=lambda: np.eye(3, dtype=np.float64))
    translation: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", as_float_array(self.matrix, shape=(3, 3)))
        object.__setattr__(self, "translation", as_float_array(self.translation, shape=(3,)))
        if abs(float(np.linalg.det(self.matrix))) < 1e-12:
            raise ValueError("Transform3D.matrix must be invertible (non-zero determinant).")

    @classmethod
    def identity(cls) -> Transform3D:
        """The identity placement (no rotation, no translation)."""

        return cls()

    @classmethod
    def from_matrix(cls, matrix: Any, *, translation: Any = (0.0, 0.0, 0.0)) -> Transform3D:
        """Build from an explicit ``(3, 3)`` linear part and translation."""

        return cls(matrix=as_float_array(matrix, shape=(3, 3)), translation=translation)

    @classmethod
    def from_rotation(cls, rotation: Any, *, translation: Any = (0.0, 0.0, 0.0)) -> Transform3D:
        """Build from a `pytex.core.orientation.Rotation` and a translation.

        The rotation acts on world geometry directly: ``x_world = R x_local``.
        """

        return cls(matrix=rotation.as_matrix(), translation=translation)

    @classmethod
    def from_orientation(
        cls,
        orientation: Any,
        *,
        translation: Any = (0.0, 0.0, 0.0),
        sense: str = "crystal_to_sample",
    ) -> Transform3D:
        """Build a crystal placement from an `Orientation`.

        With ``sense="crystal_to_sample"`` (the default) the transform maps
        crystal-frame Cartesian geometry into the specimen/world frame, matching
        `Orientation.map_crystal_vector`; ``sense="sample_to_crystal"`` uses the
        inverse. Use the default when drawing a grain's structure in the sample
        frame it was measured in.
        """

        matrix = np.asarray(orientation.as_matrix(), dtype=np.float64)
        if sense == "crystal_to_sample":
            linear = matrix
        elif sense == "sample_to_crystal":
            linear = matrix.T
        else:
            raise ValueError("sense must be 'crystal_to_sample' or 'sample_to_crystal'.")
        return cls(matrix=linear, translation=translation)

    @property
    def is_rigid(self) -> bool:
        """Whether the linear part is a proper rotation (lengths/angles preserved)."""

        return is_rotation_matrix(self.matrix)

    def apply_points(self, points: Any) -> np.ndarray:
        """Map an ``(n, 3)`` (or ``(3,)``) point array into the world frame."""

        array = np.asarray(points, dtype=np.float64)
        single = array.ndim == 1
        stacked = array.reshape(1, 3) if single else array
        mapped = stacked @ self.matrix.T + self.translation[None, :]
        return mapped[0] if single else np.ascontiguousarray(mapped)

    def apply_vector(self, vectors: Any) -> np.ndarray:
        """Map direction vectors (no translation) into the world frame."""

        array = np.asarray(vectors, dtype=np.float64)
        single = array.ndim == 1
        stacked = array.reshape(1, 3) if single else array
        mapped = stacked @ self.matrix.T
        return mapped[0] if single else np.ascontiguousarray(mapped)

    def apply_normal(self, normals: Any) -> np.ndarray:
        """Map plane normals covariantly (via the inverse-transpose).

        For a proper rotation this coincides with `apply_vector`; for a scaled
        or sheared transform it keeps a normal perpendicular to its plane.
        Returned vectors are re-normalized to unit length.
        """

        array = np.asarray(normals, dtype=np.float64)
        single = array.ndim == 1
        stacked = array.reshape(1, 3) if single else array
        cotransform = np.linalg.inv(self.matrix)
        mapped = stacked @ cotransform
        lengths = np.linalg.norm(mapped, axis=1, keepdims=True)
        mapped = mapped / np.where(lengths == 0.0, 1.0, lengths)
        return mapped[0] if single else np.ascontiguousarray(mapped)

    def compose(self, other: Transform3D) -> Transform3D:
        """Return ``self ∘ other`` (apply ``other`` first, then ``self``)."""

        return Transform3D(
            matrix=self.matrix @ other.matrix,
            translation=self.matrix @ other.translation + self.translation,
        )

    def inverse(self) -> Transform3D:
        """Return the inverse placement."""

        inverse_matrix = np.linalg.inv(self.matrix)
        return Transform3D(
            matrix=inverse_matrix,
            translation=-inverse_matrix @ self.translation,
        )


@dataclass(frozen=True, slots=True)
class Arrow3D:
    """A 3D vector drawn tail -> head, the base primitive of the vocabulary."""

    tail: np.ndarray
    head: np.ndarray
    color: str = "#2563eb"
    label: str | None = None
    linewidth: float = 2.2
    arrow_ratio: float = 0.14
    alpha: float = 0.97
    label_color: str | None = None
    fontsize: float = 11.0
    label_offset_fraction: float = 0.03

    def __post_init__(self) -> None:
        object.__setattr__(self, "tail", as_float_array(self.tail, shape=(3,)))
        object.__setattr__(self, "head", as_float_array(self.head, shape=(3,)))
        if np.allclose(self.tail, self.head):
            raise ValueError("Arrow3D requires distinct tail and head points.")

    @property
    def vector(self) -> np.ndarray:
        return as_float_array(self.head - self.tail, shape=(3,))

    def transformed(self, transform: Transform3D) -> Arrow3D:
        return Arrow3D(
            tail=transform.apply_points(self.tail),
            head=transform.apply_points(self.head),
            color=self.color,
            label=self.label,
            linewidth=self.linewidth,
            arrow_ratio=self.arrow_ratio,
            alpha=self.alpha,
            label_color=self.label_color,
            fontsize=self.fontsize,
            label_offset_fraction=self.label_offset_fraction,
        )


@dataclass(frozen=True, slots=True)
class PolyLine3D:
    """An open or closed polyline through an ordered ``(n, 3)`` point set."""

    points: np.ndarray
    color: str = "#334155"
    linewidth: float = 1.4
    linestyle: str = "-"
    alpha: float = 1.0
    closed: bool = False
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "points", as_float_array(self.points, shape=(None, 3)))
        if self.points.shape[0] < 2:
            raise ValueError("PolyLine3D requires at least two points.")

    def transformed(self, transform: Transform3D) -> PolyLine3D:
        return PolyLine3D(
            points=transform.apply_points(self.points),
            color=self.color,
            linewidth=self.linewidth,
            linestyle=self.linestyle,
            alpha=self.alpha,
            closed=self.closed,
            label=self.label,
        )


@dataclass(frozen=True, slots=True)
class PlanePatch3D:
    """A planar polygon with an outward unit normal and optional label."""

    vertices: np.ndarray
    normal: np.ndarray
    color: str = "#0f766e"
    alpha: float = 0.28
    edge_color: str | None = None
    edge_width: float = 0.8
    label: str | None = None
    label_color: str | None = None
    fontsize: float = 11.0
    label_offset_fraction: float = 0.035

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", as_float_array(self.vertices, shape=(None, 3)))
        if self.vertices.shape[0] < 3:
            raise ValueError("PlanePatch3D requires at least three vertices.")
        normal = np.asarray(self.normal, dtype=np.float64)
        norm = float(np.linalg.norm(normal))
        if np.isclose(norm, 0.0):
            raise ValueError("PlanePatch3D.normal must be non-zero.")
        object.__setattr__(self, "normal", as_float_array(normal / norm, shape=(3,)))

    def transformed(self, transform: Transform3D) -> PlanePatch3D:
        return PlanePatch3D(
            vertices=transform.apply_points(self.vertices),
            normal=transform.apply_normal(self.normal),
            color=self.color,
            alpha=self.alpha,
            edge_color=self.edge_color,
            edge_width=self.edge_width,
            label=self.label,
            label_color=self.label_color,
            fontsize=self.fontsize,
            label_offset_fraction=self.label_offset_fraction,
        )


@dataclass(frozen=True, slots=True)
class PointCloud3D:
    """A set of point markers (lattice nodes, atom sites, projected poles)."""

    points: np.ndarray
    color: str | Sequence[str] | np.ndarray = "#1f3a5f"
    size: float = 40.0
    marker: str = "o"
    alpha: float = 0.95
    edgecolor: str | None = "#f8fafc"
    linewidth: float = 0.3
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "points", as_float_array(self.points, shape=(None, 3)))
        if isinstance(self.color, np.ndarray):
            object.__setattr__(self, "color", as_float_array(self.color, shape=(None, 3)))

    def transformed(self, transform: Transform3D) -> PointCloud3D:
        return PointCloud3D(
            points=transform.apply_points(self.points),
            color=self.color,
            size=self.size,
            marker=self.marker,
            alpha=self.alpha,
            edgecolor=self.edgecolor,
            linewidth=self.linewidth,
            label=self.label,
        )


@dataclass(frozen=True, slots=True)
class Label3D:
    """A free text label anchored at a world-frame point."""

    position: np.ndarray
    text: str
    color: str = "#111111"
    fontsize: float = 11.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", as_float_array(self.position, shape=(3,)))

    def transformed(self, transform: Transform3D) -> Label3D:
        return Label3D(
            position=transform.apply_points(self.position),
            text=self.text,
            color=self.color,
            fontsize=self.fontsize,
        )


@dataclass(frozen=True, slots=True)
class AxisTriad3D:
    """An orthonormal-looking axis triad (reference-frame gizmo).

    ``axes`` holds the three axis vectors as columns (already scaled to the
    desired display length); ``origin`` anchors them. Expanding a triad yields
    one `Arrow3D` per axis plus optional `Label3D` tip labels, so triads reuse
    exactly the same rendering path as any other arrow.
    """

    origin: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    axes: np.ndarray = field(default_factory=lambda: np.eye(3, dtype=np.float64))
    colors: tuple[str, str, str] = _DEFAULT_TRIAD_COLORS
    labels: tuple[str, str, str] | None = None
    linewidth: float = 2.0
    arrow_ratio: float = 0.16
    fontsize: float = 11.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin", as_float_array(self.origin, shape=(3,)))
        object.__setattr__(self, "axes", as_float_array(self.axes, shape=(3, 3)))
        if len(self.colors) != 3:
            raise ValueError("AxisTriad3D.colors must contain exactly three colors.")
        if self.labels is not None and len(self.labels) != 3:
            raise ValueError("AxisTriad3D.labels must contain exactly three labels when provided.")

    def arrows(self) -> tuple[Arrow3D, ...]:
        return tuple(
            Arrow3D(
                tail=self.origin,
                head=self.origin + self.axes[:, index],
                color=self.colors[index],
                label=None,
                linewidth=self.linewidth,
                arrow_ratio=self.arrow_ratio,
                fontsize=self.fontsize,
            )
            for index in range(3)
        )

    def tip_labels(self) -> tuple[Label3D, ...]:
        if self.labels is None:
            return ()
        return tuple(
            Label3D(
                position=self.origin + 1.08 * self.axes[:, index],
                text=self.labels[index],
                color=self.colors[index],
                fontsize=self.fontsize,
            )
            for index in range(3)
        )

    def transformed(self, transform: Transform3D) -> AxisTriad3D:
        return AxisTriad3D(
            origin=transform.apply_points(self.origin),
            axes=transform.apply_vector(self.axes.T).T,
            colors=self.colors,
            labels=self.labels,
            linewidth=self.linewidth,
            arrow_ratio=self.arrow_ratio,
            fontsize=self.fontsize,
        )


@dataclass(frozen=True, slots=True)
class PrimitiveScene3D:
    """An immutable, backend-agnostic bag of world-frame primitives.

    Compose a scene by construction or with `merge`; reposition an entire scene
    with `transformed`. `render_primitive_scene_3d` draws it standalone, and the
    composite world-scene renderer draws it alongside placed crystals — so the
    primitives that annotate a lone stereographic vector and those that annotate
    a two-crystal orientation-relationship figure are described identically.
    """

    arrows: tuple[Arrow3D, ...] = ()
    polylines: tuple[PolyLine3D, ...] = ()
    patches: tuple[PlanePatch3D, ...] = ()
    point_clouds: tuple[PointCloud3D, ...] = ()
    labels: tuple[Label3D, ...] = ()
    triads: tuple[AxisTriad3D, ...] = ()

    def is_empty(self) -> bool:
        return not (
            self.arrows
            or self.polylines
            or self.patches
            or self.point_clouds
            or self.labels
            or self.triads
        )

    def merge(self, other: PrimitiveScene3D) -> PrimitiveScene3D:
        return PrimitiveScene3D(
            arrows=self.arrows + other.arrows,
            polylines=self.polylines + other.polylines,
            patches=self.patches + other.patches,
            point_clouds=self.point_clouds + other.point_clouds,
            labels=self.labels + other.labels,
            triads=self.triads + other.triads,
        )

    def transformed(self, transform: Transform3D) -> PrimitiveScene3D:
        return PrimitiveScene3D(
            arrows=tuple(arrow.transformed(transform) for arrow in self.arrows),
            polylines=tuple(line.transformed(transform) for line in self.polylines),
            patches=tuple(patch.transformed(transform) for patch in self.patches),
            point_clouds=tuple(cloud.transformed(transform) for cloud in self.point_clouds),
            labels=tuple(label.transformed(transform) for label in self.labels),
            triads=tuple(triad.transformed(transform) for triad in self.triads),
        )

    def bounds(self) -> np.ndarray:
        """Axis-aligned ``(2, 3)`` [min; max] bounds of every primitive vertex."""

        points: list[np.ndarray] = []
        for arrow in self.arrows:
            points.append(np.vstack([arrow.tail, arrow.head]))
        for line in self.polylines:
            points.append(line.points)
        for patch in self.patches:
            points.append(patch.vertices)
        for cloud in self.point_clouds:
            points.append(cloud.points)
        for label in self.labels:
            points.append(label.position[None, :])
        for triad in self.triads:
            corners = triad.origin[None, :] + np.vstack([np.zeros(3), triad.axes.T])
            points.append(corners)
        if not points:
            return np.zeros((2, 3), dtype=np.float64)
        stacked = np.vstack(points)
        return np.vstack([np.min(stacked, axis=0), np.max(stacked, axis=0)])


# --------------------------------------------------------------------------- #
# Builders: canonical crystallographic objects -> primitives
# --------------------------------------------------------------------------- #


def vector_arrow(
    vector: Any,
    *,
    origin: Any = (0.0, 0.0, 0.0),
    scale: float = 1.0,
    color: str = "#2563eb",
    label: str | None = None,
    **arrow_kwargs: Any,
) -> Arrow3D:
    """Build an `Arrow3D` for a bare Cartesian 3D vector.

    The vector is drawn from ``origin`` to ``origin + scale * vector``. Use this
    for any world-frame vector that is not tied to a crystal (a loading axis, a
    displacement, a diffraction ``g``-vector already in Cartesian form).
    """

    tail = as_float_array(origin, shape=(3,))
    displacement = float(scale) * as_float_array(vector, shape=(3,))
    return Arrow3D(tail=tail, head=tail + displacement, color=color, label=label, **arrow_kwargs)


def direction_arrow(
    direction: CrystalDirection,
    *,
    origin: Any = (0.0, 0.0, 0.0),
    length: float | None = None,
    color: str = "#2563eb",
    label: str | Sequence[int] | None = None,
    **arrow_kwargs: Any,
) -> Arrow3D:
    """Build an `Arrow3D` for a `CrystalDirection` in its crystal Cartesian frame.

    The arrow points along the direction's `unit_vector`. ``length`` sets the
    arrow length in angstrom (default: the longest direct-cell edge, so the
    arrow reads at the scale of the cell); when ``label`` is omitted an integer
    direction is auto-labelled as ``[uvw]``.
    """

    unit = direction.unit_vector
    if length is None:
        basis = direction.phase.lattice.direct_basis().matrix
        length = float(np.max(np.linalg.norm(basis, axis=0)))
    text = _direction_label_text(direction, label)
    return vector_arrow(
        unit,
        origin=origin,
        scale=float(length),
        color=color,
        label=text,
        **arrow_kwargs,
    )


def plane_normal_arrow(
    plane: CrystalPlane,
    *,
    origin: Any = (0.0, 0.0, 0.0),
    length: float | None = None,
    color: str = "#0f766e",
    label: str | Sequence[int] | None = None,
    **arrow_kwargs: Any,
) -> Arrow3D:
    """Build an `Arrow3D` along a `CrystalPlane` normal (the plane pole)."""

    unit = plane.normal
    if length is None:
        basis = plane.phase.lattice.direct_basis().matrix
        length = float(np.max(np.linalg.norm(basis, axis=0)))
    text = _plane_label_text(plane, label)
    return vector_arrow(
        unit,
        origin=origin,
        scale=float(length),
        color=color,
        label=text,
        **arrow_kwargs,
    )


def crystal_plane_patch(
    plane: CrystalPlane,
    *,
    center: Any = (0.0, 0.0, 0.0),
    extent: float | None = None,
    offset: float = 0.0,
    color: str = "#0f766e",
    alpha: float = 0.28,
    label: str | Sequence[int] | None = None,
    **patch_kwargs: Any,
) -> PlanePatch3D:
    """Build a square `PlanePatch3D` representing a `CrystalPlane`.

    The patch is centered at ``center + offset * n`` (``n`` the unit normal) and
    spans ``2 * extent`` angstrom on each in-plane axis. Use it to show a lattice
    plane, a slip plane, or a habit plane as a translucent sheet in a 3D scene.
    """

    normal = plane.normal
    if extent is None:
        basis = plane.phase.lattice.direct_basis().matrix
        extent = float(np.max(np.linalg.norm(basis, axis=0)))
    u_axis, v_axis = _in_plane_axes(normal)
    anchor = as_float_array(center, shape=(3,)) + float(offset) * normal
    vertices = np.vstack(
        [
            anchor + extent * (u_axis + v_axis),
            anchor + extent * (-u_axis + v_axis),
            anchor + extent * (-u_axis - v_axis),
            anchor + extent * (u_axis - v_axis),
        ]
    )
    text = _plane_label_text(plane, label)
    return PlanePatch3D(
        vertices=vertices,
        normal=normal,
        color=color,
        alpha=alpha,
        label=text,
        **patch_kwargs,
    )


def reference_frame_triad(
    frame: Any = None,
    *,
    basis: Any = None,
    origin: Any = (0.0, 0.0, 0.0),
    length: float = 1.0,
    colors: tuple[str, str, str] = _DEFAULT_TRIAD_COLORS,
    labels: tuple[str, str, str] | None = None,
    orthonormalize: bool = True,
    **triad_kwargs: Any,
) -> AxisTriad3D:
    """Build an `AxisTriad3D` for a reference frame or an explicit basis.

    Provide either a `ReferenceFrame` — whose axis labels name the triad and
    whose `ReferenceFrame.axis_vectors` set where the axes point, so a frame
    recorded as rotated draws rotated — or an explicit ``(3, 3)`` ``basis``
    whose columns are the axis vectors (e.g. a `Lattice.direct_basis().matrix`
    for the crystal ``a/b/c`` axes). An explicit ``basis`` wins over the frame's
    own geometry. ``length`` scales the drawn axes. With ``orthonormalize`` a
    non-orthonormal basis is drawn as unit direction arrows (so an oblique cell
    still yields a legible gizmo); set it ``False`` to draw the true edge
    vectors.
    """

    if basis is None and frame is not None and getattr(frame, "basis_matrix", None) is not None:
        basis = frame.basis_matrix
    if basis is not None:
        matrix = as_float_array(basis, shape=(3, 3))
        if orthonormalize:
            norms = np.linalg.norm(matrix, axis=0, keepdims=True)
            matrix = matrix / np.where(norms == 0.0, 1.0, norms)
        axes = float(length) * matrix
    else:
        axes = float(length) * np.eye(3, dtype=np.float64)
    if labels is None and frame is not None and getattr(frame, "axes", None) is not None:
        frame_axes = tuple(str(name) for name in frame.axes)
        if len(frame_axes) == 3:
            labels = frame_axes
    return AxisTriad3D(
        origin=as_float_array(origin, shape=(3,)),
        axes=axes,
        colors=colors,
        labels=labels,
        **triad_kwargs,
    )


def unit_cell_polylines(
    source: Phase | Lattice,
    *,
    origin: Any = (0.0, 0.0, 0.0),
    repeats: tuple[int, int, int] = (1, 1, 1),
    color: str = "#334155",
    linewidth: float = 1.2,
    alpha: float = 0.9,
) -> tuple[PolyLine3D, ...]:
    """Build the twelve edge polylines of a (super)cell parallelepiped.

    ``source`` is a `Phase` or a `Lattice`; ``repeats`` extends the box to a
    supercell block. Every edge is one `PolyLine3D`, so the cell frame composes
    with any other primitives and can be placed by a `Transform3D`.
    """

    lattice = source.lattice if isinstance(source, Phase) else source
    basis = lattice.direct_basis().matrix
    anchor = as_float_array(origin, shape=(3,))
    span = np.array(repeats, dtype=np.float64)
    unit_offsets = np.array(
        [
            [0.0, 0.0, 0.0],
            [span[0], 0.0, 0.0],
            [0.0, span[1], 0.0],
            [0.0, 0.0, span[2]],
            [span[0], span[1], 0.0],
            [span[0], 0.0, span[2]],
            [0.0, span[1], span[2]],
            [span[0], span[1], span[2]],
        ],
        dtype=np.float64,
    )
    corners = anchor[None, :] + (basis @ unit_offsets.T).T
    edge_pairs = (
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
        (2, 4),
        (2, 6),
        (3, 5),
        (3, 6),
        (4, 7),
        (5, 7),
        (6, 7),
    )
    return tuple(
        PolyLine3D(
            points=np.vstack([corners[a], corners[b]]),
            color=color,
            linewidth=linewidth,
            alpha=alpha,
        )
        for a, b in edge_pairs
    )


def lattice_point_cloud(
    source: Phase | Lattice,
    *,
    origin: Any = (0.0, 0.0, 0.0),
    repeats: tuple[int, int, int] = (1, 1, 1),
    color: str = "#1f3a5f",
    size: float = 42.0,
    **cloud_kwargs: Any,
) -> PointCloud3D:
    """Build a `PointCloud3D` of Bravais lattice nodes over a supercell block."""

    if any(value < 1 for value in repeats):
        raise ValueError("lattice_point_cloud repeats must contain positive integers.")
    lattice = source.lattice if isinstance(source, Phase) else source
    basis = lattice.direct_basis().matrix
    anchor = as_float_array(origin, shape=(3,))
    grid = np.array(
        [
            [i, j, k]
            for i in range(repeats[0] + 1)
            for j in range(repeats[1] + 1)
            for k in range(repeats[2] + 1)
        ],
        dtype=np.float64,
    )
    points = anchor[None, :] + (basis @ grid.T).T
    return PointCloud3D(points=points, color=color, size=size, **cloud_kwargs)


def _in_plane_axes(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unit = np.asarray(normal, dtype=np.float64)
    unit = unit / np.linalg.norm(unit)
    trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(trial, unit))), 1.0, atol=1e-8):
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u_axis = np.cross(unit, trial)
    u_axis /= np.linalg.norm(u_axis)
    v_axis = np.cross(unit, u_axis)
    return u_axis, v_axis


def _direction_label_text(
    direction: CrystalDirection, label: str | Sequence[int] | None
) -> str | None:
    if isinstance(label, str):
        return label
    if label is not None:
        return format_direction_indices(tuple(int(value) for value in label))
    rounded = np.rint(direction.coordinates).astype(np.int64)
    if np.allclose(rounded.astype(np.float64), direction.coordinates, atol=1e-8):
        return format_direction_indices(tuple(int(value) for value in rounded))
    return None


def _plane_label_text(plane: CrystalPlane, label: str | Sequence[int] | None) -> str | None:
    if isinstance(label, str):
        return label
    if label is not None:
        return format_plane_indices(tuple(int(value) for value in label))
    return format_plane_indices(tuple(int(value) for value in plane.miller.indices))


# --------------------------------------------------------------------------- #
# Renderer
# --------------------------------------------------------------------------- #


def _draw_primitive_scene(
    axes: Any,
    scene: PrimitiveScene3D,
    *,
    scene_span: float,
) -> None:
    """Draw every primitive of ``scene`` onto an existing 3D matplotlib axes.

    Shared by `render_primitive_scene_3d` and the composite world-scene renderer
    so standalone primitive figures and multi-crystal figures use identical
    drawing logic. ``scene_span`` scales label offsets to the figure extent.
    """

    _, poly3d_collection = _require_matplotlib()
    arrows = list(scene.arrows)
    labels = list(scene.labels)
    for triad in scene.triads:
        arrows.extend(triad.arrows())
        labels.extend(triad.tip_labels())
    for patch in scene.patches:
        axes.add_collection3d(
            poly3d_collection(
                [patch.vertices],
                facecolors=patch.color,
                edgecolors=patch.edge_color or patch.color,
                linewidths=patch.edge_width,
                alpha=patch.alpha,
            )
        )
        if patch.label:
            center = np.mean(patch.vertices, axis=0)
            offset = patch.label_offset_fraction * scene_span * patch.normal
            axes.text(
                center[0] + offset[0],
                center[1] + offset[1],
                center[2] + offset[2],
                patch.label,
                color=patch.label_color or patch.color,
                fontsize=patch.fontsize,
            )
    for line in scene.polylines:
        points = line.points
        if line.closed:
            points = np.vstack([points, points[0]])
        axes.plot(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            color=line.color,
            linewidth=line.linewidth,
            linestyle=line.linestyle,
            alpha=line.alpha,
        )
    for cloud in scene.point_clouds:
        color_arg = cloud.color
        axes.scatter(
            cloud.points[:, 0],
            cloud.points[:, 1],
            cloud.points[:, 2],
            s=cloud.size,
            c=[color_arg] if isinstance(color_arg, str) else color_arg,
            marker=cloud.marker,
            alpha=cloud.alpha,
            edgecolors=cloud.edgecolor,
            linewidths=cloud.linewidth,
        )
    for arrow in arrows:
        vector = arrow.vector
        axes.quiver(
            arrow.tail[0],
            arrow.tail[1],
            arrow.tail[2],
            vector[0],
            vector[1],
            vector[2],
            color=arrow.color,
            alpha=arrow.alpha,
            arrow_length_ratio=arrow.arrow_ratio,
            linewidth=arrow.linewidth,
        )
        if arrow.label:
            unit = vector / np.linalg.norm(vector)
            label_point = arrow.head + arrow.label_offset_fraction * scene_span * unit
            axes.text(
                label_point[0],
                label_point[1],
                label_point[2],
                arrow.label,
                color=arrow.label_color or arrow.color,
                fontsize=arrow.fontsize,
            )
    for label in labels:
        axes.text(
            label.position[0],
            label.position[1],
            label.position[2],
            label.text,
            color=label.color,
            fontsize=label.fontsize,
        )


def scene_span(bounds: np.ndarray) -> float:
    """The largest side of an axis-aligned ``(2, 3)`` [min; max] bounding box."""

    extent = np.asarray(bounds, dtype=np.float64)
    return float(np.max(extent[1] - extent[0]) + 1e-6)


def render_primitive_scene_3d(
    scene: PrimitiveScene3D,
    *,
    ax: Any | None = None,
    title: str | None = None,
    show_axes: bool = False,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    elev_deg: float = 22.0,
    azim_deg: float = 34.0,
    projection: str = "persp",
) -> Any:
    """Render a `PrimitiveScene3D` as a standalone 3D figure.

    Draws every primitive (arrows, polylines, plane patches, point clouds,
    labels, triads) with equal box aspect so vectors and angles read true. Pass
    ``ax`` to draw into an existing 3D axes (for panels or overlays). Use this to
    visualize a lone 3D vector, a reference-frame triad, a unit-cell wireframe,
    or any hand-built composition that does not need the full atomistic crystal
    renderer.
    """

    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    crystal_style = style["crystal"]
    background = crystal_style["background"]
    if ax is None:
        fig = plt.figure(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=background,
        )
        axes = fig.add_subplot(111, projection="3d", proj_type=projection)
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(background)
    bounds = scene.bounds()
    span = scene_span(bounds)
    _draw_primitive_scene(axes, scene, scene_span=span)
    center = 0.5 * (bounds[0] + bounds[1])
    radius = 0.55 * span
    axes.set_xlim(center[0] - radius, center[0] + radius)
    axes.set_ylim(center[1] - radius, center[1] + radius)
    axes.set_zlim(center[2] - radius, center[2] + radius)
    axes.set_box_aspect((1.0, 1.0, 1.0))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    if bool(crystal_style.get("hide_grid", True)):
        axes.grid(False)
    if not show_axes:
        axes.set_axis_off()
    if title is not None:
        axes.set_title(title)
    fig.tight_layout()
    return fig


__all__ = [
    "Arrow3D",
    "AxisTriad3D",
    "Label3D",
    "PlanePatch3D",
    "PointCloud3D",
    "PolyLine3D",
    "PrimitiveScene3D",
    "Transform3D",
    "crystal_plane_patch",
    "direction_arrow",
    "lattice_point_cloud",
    "plane_normal_arrow",
    "reference_frame_triad",
    "render_primitive_scene_3d",
    "unit_cell_polylines",
    "vector_arrow",
]
