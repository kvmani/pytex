"""VESTA-class 3D crystal-structure scenes and their matplotlib renderer.

Architecture (kept strictly two-layer so the scene model stays portable):

1. **Scene model** — `build_crystal_scene` turns a `Phase` (lattice + unit
   cell + symmetry) into a `CrystalScene`: an immutable, renderer-independent
   scene graph of typed glyphs (`CrystalAtomGlyph`, `CrystalBondGlyph`,
   `CrystalCellGlyph`, `CrystalPlaneGlyph`, `CrystalDirectionGlyph`). Glyphs
   carry only geometry (angstrom coordinates), colors, and annotation intent;
   nothing matplotlib-specific. Any future backend (OpenGL/GUI, ray tracer,
   web) can consume a `CrystalScene` unchanged.
2. **Renderer** — `plot_crystal_structure_3d` rasterizes a scene with
   matplotlib. Atoms and bonds become quad meshes lit by a Blinn-Phong-style
   model (`_lit_face_colors`) and are drawn as ONE depth-sorted
   `Poly3DCollection`, so atoms and bonds occlude each other correctly from
   every viewing angle (per-artist painter's-order artifacts are avoided by
   construction). Bonds render two-tone (each half in its atom's color, the
   VESTA convention) whenever their glyphs carry per-end colors.

Styling comes from the YAML theme system (`crystal` section); every knob used
here (lighting strengths, mesh resolutions, render modes, bond color mode)
has a theme default and can be overridden per call via ``style_overrides``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core._chemistry import (
    covalent_radius_angstrom,
    cpk_color,
    display_radius_angstrom,
)
from pytex.core.lattice import AtomicSite, CrystalDirection, CrystalPlane, MillerIndex, Phase
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.plotting.frames import add_frame_indicator
from pytex.plotting.primitives import (
    CellRegion,
    Transform3D,
    hexagonal_prism_region,
    lattice_plane_polygon,
)
from pytex.plotting.styles import _deep_merge, resolve_style


def _to_hex(color: Any) -> str:
    """``matplotlib.colors.to_hex``, imported on demand.

    matplotlib is a required dependency, but a heavy one, and the repository
    forbids import-time coupling to heavy scientific stacks: a module-level
    ``from matplotlib.colors import to_hex`` would make every ``import pytex``
    pay for it. Importing inside the call costs one cached dict lookup.
    """

    from matplotlib.colors import to_hex

    return str(to_hex(color))


def _to_rgb(color: Any) -> tuple[float, float, float]:
    """``matplotlib.colors.to_rgb``, imported on demand. See :func:`_to_hex`."""

    from matplotlib.colors import to_rgb
    red, green, blue = to_rgb(color)
    return (float(red), float(green), float(blue))


# VESTA-style render presets: each maps to style overrides (and build behavior)
# so one keyword switches the whole visual system. User style_overrides win.
_RENDER_STYLE_OVERRIDES: dict[str, dict[str, Any]] = {
    "ball_and_stick": {},
    "space_filling": {
        "atom_radius_scale": 1.0,
        "atom_radius_kind": "atomic",
    },
    "stick": {
        "atom_radius_scale": 0.24,
        "bond_radius_scale": 1.0,
    },
    "wireframe": {
        # VESTA wireframe: the bond network alone, no atom bodies
        "bond_render_mode": "line",
        "atom_render_mode": "none",
        "atom_radius_scale": 0.3,
    },
    "polyhedral": {
        "atom_radius_scale": 0.4,
    },
}

# Minimum Euclidean sRGB separation between two species colours in one scene.
# Tuned so Jmol's Ni/Cl greens separate while chemically distinct pairs keep
# their exact CPK colours.
_SPECIES_COLOR_MIN_DISTANCE = 0.30


def _merged_style_overrides(
    render_style: str,
    style_overrides: dict[str, Any] | None,
) -> dict[str, Any] | None:
    """Merge a render-style preset under any explicit user overrides."""

    if render_style not in _RENDER_STYLE_OVERRIDES:
        raise ValueError(
            f"render_style must be one of {sorted(_RENDER_STYLE_OVERRIDES)!r}, "
            f"received {render_style!r}."
        )
    preset = _RENDER_STYLE_OVERRIDES[render_style]
    if not preset:
        return style_overrides
    merged: dict[str, Any] = {"crystal": dict(preset)}
    if style_overrides is not None:
        merged = _deep_merge(merged, style_overrides)
    return merged


def _matplotlib() -> tuple[Any, Any]:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    return plt, Poly3DCollection


@dataclass(frozen=True, slots=True)
class CrystalAtomGlyph:
    """A rendered atom: position and display radius in angstrom, CPK color.

    Partially occupied and shared (mixed-species) sites carry VESTA-style
    sector metadata: the sphere is divided azimuthally, this glyph paints the
    fraction ``[sector_start, sector_start + occupancy)`` of one full turn in
    its species color, and ``vacancy_fraction`` (set on the last glyph of a
    shared site) paints the unoccupied remainder in the theme vacancy color.
    Fully occupied sites keep the defaults and render as plain spheres.
    """

    position_angstrom: np.ndarray
    species: str
    radius_angstrom: float
    color: str
    alpha: float
    occupancy: float = 1.0
    sector_start: float = 0.0
    vacancy_fraction: float = 0.0
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "position_angstrom", as_float_array(self.position_angstrom, shape=(3,))
        )
        if not 0.0 < self.occupancy <= 1.0:
            raise ValueError("CrystalAtomGlyph.occupancy must lie in (0, 1].")
        if self.vacancy_fraction < 0.0:
            raise ValueError("CrystalAtomGlyph.vacancy_fraction must be non-negative.")

    @property
    def is_full_sphere(self) -> bool:
        """Whether this atom is drawn as a whole sphere rather than a partial one.

        Atoms on a cell boundary are drawn clipped when periodic images are
        shown, matching the VESTA convention.
        """

        return (
            self.occupancy >= 1.0 - 1e-9
            and self.sector_start <= 1e-9
            and self.vacancy_fraction <= 1e-9
        )


@dataclass(frozen=True, slots=True)
class CrystalBondGlyph:
    """A rendered bond cylinder between two atom centers.

    ``color`` is the uniform bond color; when ``start_color`` / ``end_color``
    are set (the two-tone VESTA convention: each half in its atom's color),
    renderers split the cylinder at the midpoint and color the halves
    independently. One glyph always represents one physical bond.
    """

    start_angstrom: np.ndarray
    end_angstrom: np.ndarray
    color: str
    alpha: float
    radius_angstrom: float
    start_color: str | None = None
    end_color: str | None = None
    start_species: str | None = None
    end_species: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_angstrom", as_float_array(self.start_angstrom, shape=(3,)))
        object.__setattr__(self, "end_angstrom", as_float_array(self.end_angstrom, shape=(3,)))

    @property
    def length_angstrom(self) -> float:
        """Length of the bond, in angstroms.
        """

        return float(np.linalg.norm(self.end_angstrom - self.start_angstrom))

    def half_segments(self) -> tuple[tuple[np.ndarray, np.ndarray, str], ...]:
        """Return renderable (start, end, color) segments honoring two-tone."""

        if self.start_color is None and self.end_color is None:
            return ((self.start_angstrom, self.end_angstrom, self.color),)
        midpoint = 0.5 * (self.start_angstrom + self.end_angstrom)
        return (
            (self.start_angstrom, midpoint, self.start_color or self.color),
            (midpoint, self.end_angstrom, self.end_color or self.color),
        )


@dataclass(frozen=True, slots=True)
class CrystalPolyhedronGlyph:
    """A coordination polyhedron around a central atom.

    ``triangles_angstrom`` holds the convex-hull triangles of the bonded
    neighbor positions as an ``(n_faces, 3, 3)`` array, with matching outward
    unit ``face_normals`` ``(n_faces, 3)``. Renderers draw the faces
    translucent over the ball-and-stick model (the classic VESTA polyhedral
    view) and may outline the edges.
    """

    center_angstrom: np.ndarray
    center_species: str
    triangles_angstrom: np.ndarray
    face_normals: np.ndarray
    color: str
    alpha: float
    edge_color: str
    edge_width: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "center_angstrom", as_float_array(self.center_angstrom, shape=(3,))
        )
        triangles = np.ascontiguousarray(
            np.asarray(self.triangles_angstrom, dtype=np.float64)
        )
        normals = np.ascontiguousarray(np.asarray(self.face_normals, dtype=np.float64))
        if triangles.ndim != 3 or triangles.shape[1:] != (3, 3):
            raise ValueError("triangles_angstrom must have shape (n_faces, 3, 3).")
        if normals.shape != (triangles.shape[0], 3):
            raise ValueError("face_normals must have shape (n_faces, 3).")
        triangles.setflags(write=False)
        normals.setflags(write=False)
        object.__setattr__(self, "triangles_angstrom", triangles)
        object.__setattr__(self, "face_normals", normals)


@dataclass(frozen=True, slots=True)
class CrystalCellOverlay:
    """A request to draw unit-cell outlines in a crystal-structure view.

    A declarative overlay: it states *what* to draw, and the renderer decides
    how. Kept separate from the resulting :class:`CrystalCellGlyph` so scene
    description stays independent of rendering.

    Attributes
    ----------
    kind : str
        Cell outline style.
    anchor_fractional : np.ndarray
        Where the outline starts, in fractional cell coordinates.
    span_cells : tuple of int
        How many cells to span along each axis.
    color, alpha, linewidth, show_faces, face_alpha : optional
        Appearance; ``None`` defers to the active style theme.
    """

    kind: str = "parallelepiped"
    anchor_fractional: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    span_cells: tuple[int, int, int] = (1, 1, 1)
    color: str | None = None
    alpha: float | None = None
    linewidth: float | None = None
    show_faces: bool = False
    face_alpha: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "anchor_fractional", as_float_array(self.anchor_fractional, shape=(3,))
        )
        if self.kind not in {"parallelepiped", "hexagonal_prism"}:
            raise ValueError(
                "CrystalCellOverlay.kind must be either 'parallelepiped' or 'hexagonal_prism'."
            )
        if any(value <= 0 for value in self.span_cells):
            raise ValueError(
                "CrystalCellOverlay.span_cells must contain strictly positive integers."
            )


@dataclass(frozen=True, slots=True)
class PlaneAnnotationStyle:
    """Label appearance for a crystal plane drawn in a 3-D view.

    Attributes
    ----------
    color : str, optional
        ``None`` defers to the active style theme.
    fontsize : float
    offset_fraction : float
        Label offset from the plane, as a fraction of the scene extent, so
        labels scale with the figure rather than being fixed in points.
    """

    color: str | None = None
    fontsize: float = 11.0
    offset_fraction: float = 0.035


@dataclass(frozen=True, slots=True)
class DirectionAnnotationStyle:
    """Label appearance for a crystal direction drawn in a 3-D view.

    Attributes
    ----------
    color : str, optional
        ``None`` defers to the active style theme.
    fontsize : float
    offset_fraction : float
        Label offset from the arrow, as a fraction of the scene extent.
    """

    color: str | None = None
    fontsize: float = 11.0
    offset_fraction: float = 0.03


@dataclass(frozen=True, slots=True)
class CrystalPlaneOverlay:
    """A request to draw a crystal plane in a structure view.

    Attributes
    ----------
    plane : CrystalPlane
        The plane to draw; carries its own phase, so the patch is placed
        correctly for the lattice rather than assuming cubic geometry.
    offset : float, optional
        Displacement along the normal, for drawing a plane away from the
        origin.
    color, alpha : optional
        Appearance; ``None`` defers to the style theme.
    label : str, optional
        Explicit label text; when omitted, ``label_indices`` is formatted in
        proper crystallographic notation.
    label_indices : tuple of int, optional
    annotation_style : PlaneAnnotationStyle, optional
    """

    plane: CrystalPlane
    offset: float | None = None
    color: str | None = None
    alpha: float | None = None
    label: str | None = None
    label_indices: tuple[int, ...] | None = None
    annotation_style: PlaneAnnotationStyle | None = None


@dataclass(frozen=True, slots=True)
class CrystalDirectionOverlay:
    """A request to draw a crystal direction arrow in a structure view.

    Attributes
    ----------
    direction : CrystalDirection
        The direction to draw; resolved through the direct basis, so it is
        correct in non-cubic lattices.
    anchor_fractional : np.ndarray
        Arrow origin in fractional cell coordinates.
    color, alpha : optional
    label : str, optional
        Explicit label text; when omitted, ``label_indices`` is formatted in
        proper crystallographic notation.
    label_indices : tuple of int, optional
    annotation_style : DirectionAnnotationStyle, optional
    arrow_length_scale : float
        Arrow length as a fraction of the lattice repeat, so the arrowhead
        stays inside the cell it belongs to.
    """

    direction: CrystalDirection
    anchor_fractional: np.ndarray
    color: str | None = None
    alpha: float | None = None
    label: str | None = None
    label_indices: tuple[int, ...] | None = None
    annotation_style: DirectionAnnotationStyle | None = None
    arrow_length_scale: float = 0.92

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "anchor_fractional", as_float_array(self.anchor_fractional, shape=(3,))
        )
        if not 0.0 < self.arrow_length_scale <= 1.0:
            raise ValueError(
                "CrystalDirectionOverlay.arrow_length_scale must lie in the interval (0, 1]."
            )


@dataclass(frozen=True, slots=True)
class CrystalPlaneGlyph:
    """A crystal plane resolved to drawable geometry in Cartesian coordinates.

    The rendered counterpart of :class:`CrystalPlaneOverlay`, with the
    lattice geometry already applied.

    Attributes
    ----------
    vertices_angstrom : np.ndarray
        Polygon vertices of the plane patch.
    normal_angstrom : np.ndarray
        Plane normal, for label placement and face orientation.
    color, alpha : str, float
        Resolved appearance, no longer deferred to the theme.
    label : str
    annotation_style : PlaneAnnotationStyle
    """

    vertices_angstrom: np.ndarray
    normal_angstrom: np.ndarray
    color: str
    alpha: float
    label: str
    annotation_style: PlaneAnnotationStyle

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "vertices_angstrom", as_float_array(self.vertices_angstrom, shape=(None, 3))
        )
        object.__setattr__(
            self, "normal_angstrom", as_float_array(self.normal_angstrom, shape=(3,))
        )


@dataclass(frozen=True, slots=True)
class CrystalDirectionGlyph:
    """A crystal direction resolved to a drawable arrow in Cartesian coordinates.

    The rendered counterpart of :class:`CrystalDirectionOverlay`.

    Attributes
    ----------
    start_angstrom, end_angstrom : np.ndarray
        Arrow endpoints.
    color, alpha : str, float
        Resolved appearance.
    label : str
    annotation_style : DirectionAnnotationStyle
    """

    start_angstrom: np.ndarray
    end_angstrom: np.ndarray
    color: str
    alpha: float
    label: str
    annotation_style: DirectionAnnotationStyle

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_angstrom", as_float_array(self.start_angstrom, shape=(3,)))
        object.__setattr__(self, "end_angstrom", as_float_array(self.end_angstrom, shape=(3,)))


@dataclass(frozen=True, slots=True)
class CrystalCellGlyph:
    """Unit-cell outlines resolved to drawable geometry.

    The rendered counterpart of :class:`CrystalCellOverlay`.

    Attributes
    ----------
    kind : str
    edges_angstrom : tuple of np.ndarray
        Cell edge segments.
    faces_angstrom : tuple of np.ndarray
        Cell face polygons, drawn only when faces were requested.
    color, alpha, face_alpha, linewidth : str, float
        Resolved appearance.
    """

    kind: str
    edges_angstrom: tuple[np.ndarray, ...]
    faces_angstrom: tuple[np.ndarray, ...]
    color: str
    alpha: float
    face_alpha: float
    linewidth: float

    def __post_init__(self) -> None:
        normalized_edges = tuple(as_float_array(edge, shape=(2, 3)) for edge in self.edges_angstrom)
        normalized_faces = tuple(
            as_float_array(face, shape=(None, 3)) for face in self.faces_angstrom
        )
        object.__setattr__(self, "edges_angstrom", normalized_edges)
        object.__setattr__(self, "faces_angstrom", normalized_faces)


@dataclass(frozen=True, slots=True)
class CrystalScene:
    """Immutable, renderer-independent scene graph for one crystal structure.

    Produced by `build_crystal_scene`; consumed by `plot_crystal_structure_3d`
    (matplotlib) and by any future backend. All geometry is Cartesian
    angstrom in the crystal frame; glyph tuples are ordered but carry no
    drawing order semantics (renderers depth-sort).
    """

    phase: Phase
    atoms: tuple[CrystalAtomGlyph, ...]
    bonds: tuple[CrystalBondGlyph, ...]
    cells: tuple[CrystalCellGlyph, ...]
    planes: tuple[CrystalPlaneGlyph, ...]
    directions: tuple[CrystalDirectionGlyph, ...]
    lattice_edges: tuple[np.ndarray, ...]
    repeats: tuple[int, int, int]
    polyhedra: tuple[CrystalPolyhedronGlyph, ...] = ()

    def bounds(self) -> np.ndarray:
        """Axis-aligned bounding box of the scene, as min and max corners.

        Used to set equal-aspect 3-D axis limits so the cell is not visually
        distorted.
        """

        points: list[np.ndarray] = []
        points.extend(edge for edge in self.lattice_edges)
        points.extend(edge for cell in self.cells for edge in cell.edges_angstrom)
        # atoms contribute their full sphere extent (center +/- radius) so
        # space-filling spheres are never clipped by the axis limits
        points.extend(
            np.vstack(
                [
                    atom.position_angstrom - atom.radius_angstrom,
                    atom.position_angstrom + atom.radius_angstrom,
                ]
            )
            for atom in self.atoms
        )
        points.extend(
            np.vstack([direction.start_angstrom, direction.end_angstrom])
            for direction in self.directions
        )
        if self.planes:
            points.extend(plane.vertices_angstrom for plane in self.planes)
        stacked = np.vstack(points) if points else np.zeros((1, 3), dtype=np.float64)
        return np.vstack([np.min(stacked, axis=0), np.max(stacked, axis=0)])

    def bond_lengths_angstrom(self) -> np.ndarray:
        """Length of every bond glyph in angstrom, aligned with `bonds`.

        The programmatic analog of VESTA's interactive distance readout: every
        detected bond becomes a measurable number, so bond statistics can be
        scripted, tested, and tabulated instead of clicked one at a time.
        """

        if not self.bonds:
            return as_float_array(np.zeros(0, dtype=np.float64), shape=(0,))
        return as_float_array(
            np.array([bond.length_angstrom for bond in self.bonds], dtype=np.float64),
            shape=(len(self.bonds),),
        )

    def bond_length_summary(self) -> dict[tuple[str, str], dict[str, float]]:
        """Per species-pair bond statistics: count, min, mean, max (angstrom).

        Pairs are alphabetically ordered tuples such as ``("Cl", "Na")``. Bonds
        whose glyphs carry no species metadata are grouped under ``("?", "?")``.
        """

        grouped: dict[tuple[str, str], list[float]] = {}
        for bond in self.bonds:
            pair = tuple(sorted((bond.start_species or "?", bond.end_species or "?")))
            grouped.setdefault((pair[0], pair[1]), []).append(bond.length_angstrom)
        return {
            pair: {
                "count": float(len(lengths)),
                "min": float(np.min(lengths)),
                "mean": float(np.mean(lengths)),
                "max": float(np.max(lengths)),
            }
            for pair, lengths in sorted(grouped.items())
        }

    def transformed(self, transform: Transform3D) -> CrystalScene:
        """Return a copy with all geometry placed by ``transform``.

        Maps every glyph's world coordinates (atoms, bonds, cells, planes,
        directions, lattice box, polyhedra) through the placement while keeping
        colors, radii, and labels. Intended for **rigid** transforms (a
        `Rotation`/`Orientation` plus a translation); with a rigid transform the
        lit sphere/cylinder meshes stay geometrically exact, so a placed crystal
        renders identically to one built in its own frame. This is what lets two
        crystals in an orientation relationship share one scene.
        """

        if not transform.is_rigid:
            raise ValueError(
                "CrystalScene.transformed requires a rigid (rotation + translation) "
                "Transform3D so atom spheres and bond cylinders stay undistorted."
            )
        return CrystalScene(
            phase=self.phase,
            atoms=tuple(
                CrystalAtomGlyph(
                    position_angstrom=transform.apply_points(atom.position_angstrom),
                    species=atom.species,
                    radius_angstrom=atom.radius_angstrom,
                    color=atom.color,
                    alpha=atom.alpha,
                    occupancy=atom.occupancy,
                    sector_start=atom.sector_start,
                    vacancy_fraction=atom.vacancy_fraction,
                    label=atom.label,
                )
                for atom in self.atoms
            ),
            bonds=tuple(
                CrystalBondGlyph(
                    start_angstrom=transform.apply_points(bond.start_angstrom),
                    end_angstrom=transform.apply_points(bond.end_angstrom),
                    color=bond.color,
                    alpha=bond.alpha,
                    radius_angstrom=bond.radius_angstrom,
                    start_color=bond.start_color,
                    end_color=bond.end_color,
                    start_species=bond.start_species,
                    end_species=bond.end_species,
                )
                for bond in self.bonds
            ),
            cells=tuple(
                CrystalCellGlyph(
                    kind=cell.kind,
                    edges_angstrom=tuple(
                        transform.apply_points(edge) for edge in cell.edges_angstrom
                    ),
                    faces_angstrom=tuple(
                        transform.apply_points(face) for face in cell.faces_angstrom
                    ),
                    color=cell.color,
                    alpha=cell.alpha,
                    face_alpha=cell.face_alpha,
                    linewidth=cell.linewidth,
                )
                for cell in self.cells
            ),
            planes=tuple(
                CrystalPlaneGlyph(
                    vertices_angstrom=transform.apply_points(plane.vertices_angstrom),
                    normal_angstrom=transform.apply_normal(plane.normal_angstrom),
                    color=plane.color,
                    alpha=plane.alpha,
                    label=plane.label,
                    annotation_style=plane.annotation_style,
                )
                for plane in self.planes
            ),
            directions=tuple(
                CrystalDirectionGlyph(
                    start_angstrom=transform.apply_points(direction.start_angstrom),
                    end_angstrom=transform.apply_points(direction.end_angstrom),
                    color=direction.color,
                    alpha=direction.alpha,
                    label=direction.label,
                    annotation_style=direction.annotation_style,
                )
                for direction in self.directions
            ),
            lattice_edges=tuple(transform.apply_points(edge) for edge in self.lattice_edges),
            repeats=self.repeats,
            polyhedra=tuple(
                CrystalPolyhedronGlyph(
                    center_angstrom=transform.apply_points(polyhedron.center_angstrom),
                    center_species=polyhedron.center_species,
                    triangles_angstrom=transform.apply_points(
                        polyhedron.triangles_angstrom.reshape(-1, 3)
                    ).reshape(polyhedron.triangles_angstrom.shape),
                    face_normals=transform.apply_normal(polyhedron.face_normals),
                    color=polyhedron.color,
                    alpha=polyhedron.alpha,
                    edge_color=polyhedron.edge_color,
                    edge_width=polyhedron.edge_width,
                )
                for polyhedron in self.polyhedra
            ),
        )


def has_hexagonal_axes(phase: Phase) -> bool:
    """Whether this lattice is on hexagonal axes: ``a = b``, 90/90/120.

    The test the prism needs, and the reason it is a function: "hexagonal" as a
    *crystal system* is a symmetry statement, while the prism is a statement
    about the axes the cell is written on. Trigonal phases on hexagonal axes
    get the prism too, and correctly so — the drawing follows the cell.
    """

    lattice = phase.lattice
    return bool(
        np.isclose(lattice.a, lattice.b, atol=1e-6)
        and np.isclose(lattice.alpha_deg, 90.0, atol=1e-6)
        and np.isclose(lattice.beta_deg, 90.0, atol=1e-6)
        and np.isclose(lattice.gamma_deg, 120.0, atol=1e-6)
    )


def prism_axis_anchor(phase: Phase) -> np.ndarray:
    """Where to put the prism's axis: through a column of atoms.

    The prism's placement is a free choice — it is a drawing of the lattice,
    not a cell of it — and the choice that gives the familiar picture is an
    axis through an atomic column, so the six corner columns carry atoms too.
    A site already at the cell origin keeps the origin; otherwise the first
    site's basal coordinates are used, which for hcp written with sites at
    ``(1/3, 2/3)`` and ``(2/3, 1/3)`` is the column that makes the corners.
    """

    if phase.unit_cell is None or not phase.unit_cell.sites:
        return np.zeros(2, dtype=np.float64)
    for site in phase.unit_cell.sites:
        basal = np.asarray(site.fractional_coordinates, dtype=np.float64)[:2]
        if np.allclose(basal % 1.0, 0.0, atol=1e-9):
            return np.zeros(2, dtype=np.float64)
    first = np.asarray(phase.unit_cell.sites[0].fractional_coordinates, dtype=np.float64)
    return np.asarray(first[:2], dtype=np.float64)


def prism_region_for(phase: Phase, repeats: tuple[int, int, int]) -> CellRegion:
    """The prism this phase is drawn as, at these repeats."""

    return hexagonal_prism_region(
        scale=max(repeats[0], repeats[1]),
        height=repeats[2],
        anchor=prism_axis_anchor(phase),
    )


def _hexagonal_prism_atom_positions(
    phase: Phase,
    repeats: tuple[int, int, int],
    *,
    include_boundary_atoms: bool = True,
) -> list[tuple[AtomicSite, np.ndarray]]:
    """Every atom of the hexagonal prism, which is three rhombic cells.

    The prism is not a supercell of the rhombic cell, so its atoms cannot come
    from a block of translations: they come from a *larger* block, filtered by
    the prism itself. Drawn without the filter the picture is a parallelogram
    of atoms with a hexagon outlined over part of it, which is worse than
    either alone.

    Boundary atoms are kept, as everywhere else in this module: an atom on the
    prism wall belongs to the drawing of the prism, and dropping it would leave
    the six corner columns — the columns that make it *look* hexagonal — half
    empty.
    """

    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError("Crystal visualization requires phase.unit_cell with atomic sites.")
    scale = int(max(repeats[0], repeats[1]))
    height = int(repeats[2])
    region = prism_region_for(phase, repeats)
    direct_basis = phase.lattice.direct_basis().matrix
    tolerance = 1e-9 if include_boundary_atoms else -1e-9
    atoms: list[tuple[AtomicSite, np.ndarray]] = []
    seen: set[tuple[str, tuple[float, float, float]]] = set()
    for site in phase.unit_cell.sites:
        base = np.asarray(site.fractional_coordinates, dtype=np.float64)
        for i in range(-scale - 1, scale + 2):
            for j in range(-scale - 1, scale + 2):
                for k in range(-1, height + 2):
                    frac = base + np.array([i, j, k], dtype=np.float64)
                    if not region.contains(frac, tolerance=tolerance):
                        continue
                    position = direct_basis @ frac
                    key = (
                        site.species,
                        (
                            round(float(position[0]), 6),
                            round(float(position[1]), 6),
                            round(float(position[2]), 6),
                        ),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    atoms.append((site, position))
    return atoms


def _supercell_atom_positions(
    phase: Phase,
    repeats: tuple[int, int, int],
    *,
    include_boundary_atoms: bool = True,
) -> list[tuple[AtomicSite, np.ndarray]]:
    """Atom (site, Cartesian position) list for the repeated cell block.

    With ``include_boundary_atoms`` (the VESTA convention), sites lying on
    the supercell boundary also appear as translated copies on the opposite
    faces/edges/corners, so the displayed block looks complete: a corner atom
    of a single cell renders eight times. Coincident duplicates (e.g. a
    boundary copy landing on an explicitly listed far-face site) are removed.
    Distinct sites sharing one position (mixed occupancy) both survive, so the
    renderer can draw them as occupancy sectors of one sphere.
    """

    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError("Crystal visualization requires phase.unit_cell with atomic sites.")
    direct_basis = phase.lattice.direct_basis().matrix
    tolerance = 1e-9
    atoms: list[tuple[AtomicSite, np.ndarray]] = []
    seen: set[tuple[str, tuple[float, float, float]]] = set()
    for site in phase.unit_cell.sites:
        base = np.asarray(site.fractional_coordinates, dtype=np.float64)
        axis_translations = []
        for axis in range(3):
            if include_boundary_atoms:
                low = int(np.ceil(-base[axis] - tolerance))
                high = int(np.floor(repeats[axis] - base[axis] + tolerance))
                axis_translations.append(range(low, high + 1))
            else:
                axis_translations.append(range(repeats[axis]))
        for i in axis_translations[0]:
            for j in axis_translations[1]:
                for k in axis_translations[2]:
                    frac = base + np.array([i, j, k], dtype=np.float64)
                    position = direct_basis @ frac
                    # dedup by (species, position): boundary copies collapse onto
                    # explicitly listed far-face sites of the same species, while
                    # mixed-species shared sites survive for occupancy sectors
                    key = (
                        site.species,
                        (round(position[0], 6), round(position[1], 6), round(position[2], 6)),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    atoms.append((site, position))
    return atoms


def _supercell_box_edges(phase: Phase, repeats: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    corners = _supercell_corners_cartesian(phase, repeats)
    return _polyhedron_edges_from_corners(corners)


def _polyhedron_edges_from_corners(corners: np.ndarray) -> tuple[np.ndarray, ...]:
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
    return tuple(np.vstack([corners[a], corners[b]]) for a, b in edge_pairs)


def _polyhedron_faces_from_corners(corners: np.ndarray) -> tuple[np.ndarray, ...]:
    face_indices = (
        (0, 1, 4, 2),
        (0, 1, 5, 3),
        (0, 2, 6, 3),
        (7, 4, 1, 5),
        (7, 4, 2, 6),
        (7, 5, 3, 6),
    )
    return tuple(np.vstack([corners[index] for index in face]) for face in face_indices)


def _supercell_corners_fractional(repeats: tuple[int, int, int]) -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [repeats[0], 0.0, 0.0],
            [0.0, repeats[1], 0.0],
            [0.0, 0.0, repeats[2]],
            [repeats[0], repeats[1], 0.0],
            [repeats[0], 0.0, repeats[2]],
            [0.0, repeats[1], repeats[2]],
            [repeats[0], repeats[1], repeats[2]],
        ],
        dtype=np.float64,
    )


def _supercell_corners_cartesian(phase: Phase, repeats: tuple[int, int, int]) -> np.ndarray:
    basis = phase.lattice.direct_basis().matrix
    corners_frac = _supercell_corners_fractional(repeats)
    return np.asarray((basis @ corners_frac.T).T, dtype=np.float64)


def _parallelepiped_corners_fractional(
    anchor_fractional: np.ndarray,
    span_cells: tuple[int, int, int],
) -> np.ndarray:
    anchor = as_float_array(anchor_fractional, shape=(3,))
    span = np.array(span_cells, dtype=np.float64)
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
    return anchor + unit_offsets


def _cell_overlay_cartesian(
    phase: Phase,
    overlay: CrystalCellOverlay,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    basis = phase.lattice.direct_basis().matrix
    if overlay.kind == "parallelepiped":
        corners = (
            basis
            @ _parallelepiped_corners_fractional(overlay.anchor_fractional, overlay.span_cells).T
        ).T
        return _polyhedron_edges_from_corners(corners), _polyhedron_faces_from_corners(corners)
    return _hexagonal_prism_geometry(phase, overlay)


def _hexagonal_prism_geometry(
    phase: Phase,
    overlay: CrystalCellOverlay,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    lattice = phase.lattice
    if not (
        np.isclose(lattice.alpha_deg, 90.0, atol=1e-6)
        and np.isclose(lattice.beta_deg, 90.0, atol=1e-6)
        and np.isclose(lattice.gamma_deg, 120.0, atol=1e-6)
        and np.isclose(lattice.a, lattice.b, atol=1e-6)
    ):
        raise ValueError(
            "hexagonal_prism cell overlays require a lattice with a=b, "
            "alpha=beta=90 deg, gamma=120 deg."
        )
    basis = lattice.direct_basis().matrix
    # The hexagon's circumradius in cells. One is the single prism — three
    # rhombic cells, the figure every hexagonal crystal is drawn as — and
    # larger values scale it, so "cells along each axis" still means something
    # in prism mode instead of being quietly ignored.
    scale = float(max(overlay.span_cells[0], overlay.span_cells[1]))
    a_vec = scale * basis[:, 0]
    b_vec = scale * basis[:, 1]
    c_vec = float(overlay.span_cells[2]) * basis[:, 2]
    anchor_cart = basis @ overlay.anchor_fractional
    basal_vectors = (
        a_vec,
        a_vec + b_vec,
        b_vec,
        -a_vec,
        -(a_vec + b_vec),
        -b_vec,
    )
    planar_u = a_vec / np.linalg.norm(a_vec)
    basal_normal = np.cross(a_vec, b_vec)
    basal_normal /= np.linalg.norm(basal_normal)
    planar_v = np.cross(basal_normal, planar_u)
    ordered_basal = sorted(
        (anchor_cart + vector for vector in basal_vectors),
        key=lambda point: np.arctan2(
            float(np.dot(point - anchor_cart, planar_v)),
            float(np.dot(point - anchor_cart, planar_u)),
        ),
    )
    bottom = tuple(ordered_basal)
    top = tuple(point + c_vec for point in bottom)
    edges = (
        tuple(np.vstack([bottom[index], bottom[(index + 1) % 6]]) for index in range(6))
        + tuple(np.vstack([top[index], top[(index + 1) % 6]]) for index in range(6))
        + tuple(np.vstack([bottom[index], top[index]]) for index in range(6))
    )
    faces = (
        np.vstack(bottom),
        np.vstack(top),
        *[
            np.vstack(
                [
                    bottom[index],
                    bottom[(index + 1) % 6],
                    top[(index + 1) % 6],
                    top[index],
                ]
            )
            for index in range(6)
        ],
    )
    return edges, faces


def _plane_polygon_for_box(
    phase: Phase,
    hkl: tuple[int, ...],
    repeats: tuple[int, int, int],
    offset: float | None,
    *,
    region: CellRegion | None = None,
) -> np.ndarray | None:
    """The lattice plane clipped to the supercell box.

    One line, because the geometry belongs to `pytex.plotting.primitives` where
    every other overlay can reach it. It used to live here, which is why the
    orientation-relationship overlays — built in `scene3d` — drew origin-centred
    squares instead: the correct construction was in a module they did not use.

    ``offset`` of ``None`` takes the policy the shared function documents: the
    member of the family with the largest cross-section through the box.
    """

    return lattice_plane_polygon(phase, hkl, repeats=repeats, offset=offset, region=region)


def _coerce_plane_overlay(
    overlay: CrystalPlane | CrystalPlaneOverlay, *, phase: Phase
) -> CrystalPlaneOverlay:
    if isinstance(overlay, CrystalPlaneOverlay):
        if overlay.plane.phase != phase:
            raise ValueError("CrystalPlaneOverlay.plane.phase must match the scene phase.")
        return overlay
    if overlay.phase != phase:
        raise ValueError("CrystalPlane.phase must match the scene phase.")
    return CrystalPlaneOverlay(plane=overlay)


def _coerce_direction_overlay(
    overlay: CrystalDirection | CrystalDirectionOverlay,
    *,
    phase: Phase,
) -> CrystalDirectionOverlay:
    if isinstance(overlay, CrystalDirectionOverlay):
        if overlay.direction.phase != phase:
            raise ValueError("CrystalDirectionOverlay.direction.phase must match the scene phase.")
        return overlay
    if overlay.phase != phase:
        raise ValueError("CrystalDirection.phase must match the scene phase.")
    return CrystalDirectionOverlay(
        direction=overlay, anchor_fractional=np.zeros(3, dtype=np.float64)
    )


def _coerce_cell_overlay(overlay: CrystalCellOverlay, *, phase: Phase) -> CrystalCellOverlay:
    del phase
    if np.any(overlay.anchor_fractional < -1e-9):
        raise ValueError("CrystalCellOverlay.anchor_fractional must be non-negative.")
    return overlay


def _direction_fractional_vector(direction: CrystalDirection) -> np.ndarray:
    direct_basis = direction.phase.lattice.direct_basis().matrix
    cartesian = direction.unit_vector
    fractional = np.linalg.solve(direct_basis, cartesian)
    if np.allclose(fractional, 0.0):
        raise ValueError("CrystalDirection overlay produced a degenerate fractional vector.")
    return fractional


def _direction_endpoint_fractional(
    direction: CrystalDirection,
    *,
    anchor_fractional: np.ndarray,
    repeats: tuple[int, int, int],
    arrow_length_scale: float,
) -> np.ndarray:
    fractional_vector = _direction_fractional_vector(direction)
    candidates: list[float] = []
    for anchor_value, vector_value, repeat in zip(
        anchor_fractional, fractional_vector, repeats, strict=True
    ):
        if vector_value > 1e-12:
            candidates.append((repeat - anchor_value) / vector_value)
        elif vector_value < -1e-12:
            candidates.append((0.0 - anchor_value) / vector_value)
    positive_candidates = [value for value in candidates if value > 1e-12]
    if not positive_candidates:
        raise ValueError("CrystalDirection overlay does not intersect the repeated cell volume.")
    return np.asarray(
        anchor_fractional + arrow_length_scale * min(positive_candidates) * fractional_vector,
        dtype=np.float64,
    )


def _default_direction_indices(direction: CrystalDirection) -> tuple[int, ...] | None:
    rounded = np.rint(direction.coordinates).astype(np.int64)
    if np.allclose(direction.coordinates, rounded.astype(np.float64), atol=1e-8):
        return tuple(int(value) for value in rounded)
    return None


def _default_plane_indices(plane: CrystalPlane) -> tuple[int, ...]:
    return tuple(int(value) for value in plane.miller.indices)


def _coordination_polyhedra(
    atoms: tuple[CrystalAtomGlyph, ...],
    *,
    species: tuple[str, ...],
    bond_tolerance_angstrom: float,
    crystal_style: dict[str, Any],
) -> tuple[CrystalPolyhedronGlyph, ...]:
    """Convex-hull coordination polyhedra around every atom of the species.

    Neighbors are the atoms within the chemical bond cutoff of the center
    (covalent-radius sum plus tolerance, the same rule bonds use). Centers
    with fewer than four neighbors, or with degenerate (coplanar) neighbor
    sets, produce no polyhedron.
    """

    from scipy.spatial import ConvexHull, QhullError

    if not species:
        return ()
    positions = np.vstack([atom.position_angstrom for atom in atoms])
    radii = np.array([covalent_radius_angstrom(atom.species) for atom in atoms])
    polyhedra: list[CrystalPolyhedronGlyph] = []
    for index, atom in enumerate(atoms):
        if atom.species not in species:
            continue
        distances = np.linalg.norm(positions - positions[index], axis=1)
        cutoffs = radii[index] + radii + bond_tolerance_angstrom
        neighbor_mask = (distances <= cutoffs) & (distances > 1e-9)
        neighbors = positions[neighbor_mask]
        if neighbors.shape[0] < 4:
            continue
        try:
            hull = ConvexHull(neighbors)
        except QhullError:
            continue
        triangles = neighbors[hull.simplices]
        centers = triangles.mean(axis=1)
        normals = np.asarray(hull.equations[:, :3], dtype=np.float64)
        lengths = np.linalg.norm(normals, axis=1, keepdims=True)
        normals = normals / np.where(lengths == 0.0, 1.0, lengths)
        # orient outward from the central atom (hull equations point outward
        # from the hull interior, which contains the center)
        flip = np.einsum("ni,ni->n", centers - positions[index], normals) < 0.0
        normals[flip] *= -1.0
        polyhedra.append(
            CrystalPolyhedronGlyph(
                center_angstrom=positions[index],
                center_species=atom.species,
                triangles_angstrom=triangles,
                face_normals=normals,
                color=str(crystal_style.get("polyhedron_color") or atom.color),
                alpha=float(crystal_style["polyhedron_alpha"]),
                edge_color=str(crystal_style["polyhedron_edge_color"]),
                edge_width=float(crystal_style["polyhedron_edge_width"]),
            )
        )
    return tuple(polyhedra)


def _relative_luminance(color: str) -> float:
    """Perceived brightness of a colour on the 0-1 sRGB luminance scale."""

    red, green, blue = _to_rgb(color)
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def _rgb_tuple(values: np.ndarray) -> tuple[float, float, float]:
    """Matplotlib-facing RGB triple from a length-3 float array."""

    red, green, blue = (float(value) for value in values)
    return (red, green, blue)


def _separated_species_colors(species: Sequence[str]) -> dict[str, str]:
    """CPK colours for one scene's species, pulled apart where they collide.

    The CPK/Jmol palette assigns colours per element with no regard for which
    elements share a structure, so a real phase can put two near-identical
    colours side by side — NiCl2 draws Jmol's Ni (#50d050) against its Cl
    (#1ff01f) and the two species become indistinguishable, which defeats the
    purpose of colouring by element at all. Every species keeps its CPK identity
    when there is no clash; when two do clash the later one (in a stable
    alphabetical order, so a figure never changes between runs) is darkened or
    lightened until it separates.
    """

    ordered = sorted({str(name) for name in species})
    assigned: dict[str, str] = {}
    placed: list[np.ndarray] = []
    for name in ordered:
        base = np.asarray(_to_rgb(cpk_color(name)), dtype=np.float64)
        candidate = base
        for _ in range(6):
            if all(
                float(np.linalg.norm(candidate - existing)) >= _SPECIES_COLOR_MIN_DISTANCE
                for existing in placed
            ):
                break
            # Darken bright colours and lighten dark ones, so the adjustment
            # always moves into the available contrast range.
            if _relative_luminance(_to_hex(_rgb_tuple(candidate))) > 0.45:
                candidate = candidate * 0.68
            else:
                candidate = candidate + (1.0 - candidate) * 0.34
            candidate = np.clip(candidate, 0.0, 1.0)
        placed.append(candidate)
        assigned[name] = _to_hex(_rgb_tuple(candidate))
    return assigned


def _atom_glyphs_from_sites(
    atom_data: list[tuple[AtomicSite, np.ndarray]],
    crystal_style: dict[str, Any],
    *,
    atom_label_mode: str,
) -> tuple[CrystalAtomGlyph, ...]:
    """Turn (site, position) pairs into atom glyphs with occupancy sectors.

    Sites sharing one Cartesian position (mixed occupancy, the VESTA pie-sphere
    case) become consecutive azimuthal sectors of one shared-radius sphere; a
    total occupancy below one leaves a vacancy sector on the last glyph. Fully
    occupied lone sites keep the plain full-sphere defaults.
    """

    if atom_label_mode not in {"none", "species", "site"}:
        raise ValueError("atom_label_mode must be 'none', 'species', or 'site'.")
    radius_scale = float(crystal_style["atom_radius_scale"])
    radius_kind = str(crystal_style.get("atom_radius_kind", "covalent"))
    alpha = float(crystal_style["atom_alpha"])
    species_colors = _separated_species_colors([site.species for site, _ in atom_data])
    configured_colors = crystal_style.get("species_colors", {})
    if isinstance(configured_colors, Mapping):
        for species, color in configured_colors.items():
            if species in species_colors:
                species_colors[species] = _to_hex(color)
    groups: dict[tuple[float, float, float], list[tuple[AtomicSite, np.ndarray]]] = {}
    for site, position in atom_data:
        key = (round(position[0], 6), round(position[1], 6), round(position[2], 6))
        groups.setdefault(key, []).append((site, position))
    glyphs: list[CrystalAtomGlyph] = []
    for members in groups.values():
        shared_radius = radius_scale * max(
            display_radius_angstrom(site.species, kind=radius_kind) for site, _ in members
        )
        total_occupancy = sum(site.occupancy for site, _ in members)
        plain = len(members) == 1 and total_occupancy >= 1.0 - 1e-9
        cumulative = 0.0
        for index, (site, position) in enumerate(members):
            label: str | None = None
            if atom_label_mode == "species":
                label = site.species
            elif atom_label_mode == "site":
                label = site.label
            is_last = index == len(members) - 1
            vacancy = max(0.0, 1.0 - total_occupancy) if is_last and not plain else 0.0
            glyphs.append(
                CrystalAtomGlyph(
                    position_angstrom=position,
                    species=site.species,
                    radius_angstrom=shared_radius,
                    color=species_colors[site.species],
                    alpha=alpha,
                    occupancy=1.0 if plain else site.occupancy,
                    sector_start=0.0 if plain else cumulative,
                    vacancy_fraction=vacancy,
                    label=label,
                )
            )
            cumulative += site.occupancy
    return tuple(glyphs)


def build_crystal_scene(
    phase: Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    render_style: str = "ball_and_stick",
    show_bonds: bool = True,
    bond_tolerance_angstrom: float = 0.45,
    include_boundary_atoms: bool = True,
    polyhedra_species: tuple[str, ...] = (),
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
    plane_overlays: tuple[CrystalPlane | CrystalPlaneOverlay, ...] = (),
    direction_overlays: tuple[CrystalDirection | CrystalDirectionOverlay, ...] = (),
    show_unit_cells: bool = False,
    hexagonal_prism: bool = False,
    cell_overlays: tuple[CrystalCellOverlay, ...] = (),
    slab_hkl: tuple[int, int, int] | None = None,
    slab_thickness_angstrom: float | None = None,
    atom_label_mode: str = "none",
    site_vectors: Mapping[str, Any] | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> CrystalScene:
    """Build a renderer-independent `CrystalScene` for one phase.

    ``render_style`` selects a VESTA-style visual system in one keyword:
    ``"ball_and_stick"`` (default), ``"space_filling"`` (atomic-radius spheres,
    bonds suppressed), ``"stick"`` (uniform thin cylinders), ``"wireframe"``
    (line bonds, marker atoms), or ``"polyhedral"`` (coordination polyhedra for
    every eligible species unless ``polyhedra_species`` narrows it). Partially
    occupied or mixed-species sites automatically render as VESTA-style
    occupancy sectors. ``atom_label_mode`` stamps per-atom text labels
    (``"species"`` or ``"site"``), and ``site_vectors`` maps site labels to
    crystal-Cartesian vectors (angstrom) drawn as arrows on every periodic copy
    of that site — the magnetic-moment / displacement-vector convention.
    """

    if any(value <= 0 for value in repeats):
        raise ValueError("repeats must contain strictly positive integers.")
    effective_overrides = _merged_style_overrides(render_style, style_overrides)
    style = resolve_style(theme=theme, style_path=style_path, overrides=effective_overrides)
    crystal_style = style["crystal"]
    if render_style == "space_filling":
        show_bonds = False
    # The prism is the conventional drawing of a hexagonal cell, and it is a
    # different volume from the rhombic block: three cells rather than one, so
    # its atoms are filtered out of a larger block rather than tiled.
    as_prism = bool(hexagonal_prism) and has_hexagonal_axes(phase)
    atom_data = (
        _hexagonal_prism_atom_positions(
            phase, repeats, include_boundary_atoms=include_boundary_atoms
        )
        if as_prism
        else _supercell_atom_positions(
            phase, repeats, include_boundary_atoms=include_boundary_atoms
        )
    )
    if slab_hkl is not None and slab_thickness_angstrom is not None:
        normal = phase.lattice.reciprocal_basis().matrix @ np.array(slab_hkl, dtype=np.float64)
        normal = normal / np.linalg.norm(normal)
        filtered: list[tuple[AtomicSite, np.ndarray]] = []
        for site, position in atom_data:
            distance = abs(float(np.dot(position, normal)) - 1.0 / np.linalg.norm(normal))
            if distance <= slab_thickness_angstrom:
                filtered.append((site, position))
        atom_data = filtered
    atoms = _atom_glyphs_from_sites(atom_data, crystal_style, atom_label_mode=atom_label_mode)
    bonds: list[CrystalBondGlyph] = []
    if show_bonds:
        two_tone = str(crystal_style.get("bond_color_mode", "two_tone")).lower() == "two_tone"
        for i, atom_i in enumerate(atoms):
            for atom_j in atoms[i + 1 :]:
                cutoff = (
                    covalent_radius_angstrom(atom_i.species)
                    + covalent_radius_angstrom(atom_j.species)
                    + bond_tolerance_angstrom
                )
                pair_distance = float(
                    np.linalg.norm(atom_i.position_angstrom - atom_j.position_angstrom)
                )
                if 1e-6 < pair_distance <= cutoff:
                    bonds.append(
                        CrystalBondGlyph(
                            start_angstrom=atom_i.position_angstrom,
                            end_angstrom=atom_j.position_angstrom,
                            color=crystal_style["bond_color"],
                            alpha=float(crystal_style["bond_alpha"]),
                            radius_angstrom=float(crystal_style["bond_radius_scale"])
                            * min(atom_i.radius_angstrom, atom_j.radius_angstrom),
                            start_color=atom_i.color if two_tone else None,
                            end_color=atom_j.color if two_tone else None,
                            start_species=atom_i.species,
                            end_species=atom_j.species,
                        )
                    )
    if not show_unit_cells:
        default_cells: tuple[CrystalCellOverlay, ...] = ()
    elif as_prism:
        # One outline, not three: the prism *is* the cell being shown, and
        # drawing the rhombi inside it would say the crystal has boundaries
        # where the picture does not.
        anchor = prism_axis_anchor(phase)
        default_cells = (
            CrystalCellOverlay(
                kind="hexagonal_prism",
                anchor_fractional=np.array([anchor[0], anchor[1], 0.0], dtype=np.float64),
                span_cells=(max(repeats[0], repeats[1]), max(repeats[0], repeats[1]), repeats[2]),
            ),
        )
    else:
        default_cells = tuple(
            CrystalCellOverlay(anchor_fractional=np.array([i, j, k], dtype=np.float64))
            for i in range(repeats[0])
            for j in range(repeats[1])
            for k in range(repeats[2])
        )
    merged_cell_overlays = default_cells + tuple(
        _coerce_cell_overlay(overlay, phase=phase) for overlay in cell_overlays
    )
    cells: list[CrystalCellGlyph] = []
    for cell_overlay in merged_cell_overlays:
        edges, faces = _cell_overlay_cartesian(phase, cell_overlay)
        cells.append(
            CrystalCellGlyph(
                kind=cell_overlay.kind,
                edges_angstrom=edges,
                faces_angstrom=faces if cell_overlay.show_faces else (),
                color=cell_overlay.color or crystal_style["cell_color"],
                alpha=float(
                    crystal_style["cell_alpha"]
                    if cell_overlay.alpha is None
                    else cell_overlay.alpha
                ),
                face_alpha=float(
                    crystal_style["cell_face_alpha"]
                    if cell_overlay.face_alpha is None
                    else cell_overlay.face_alpha
                ),
                linewidth=float(
                    crystal_style["cell_linewidth"]
                    if cell_overlay.linewidth is None
                    else cell_overlay.linewidth
                ),
            )
        )
    merged_plane_overlays = tuple(
        CrystalPlaneOverlay(
            plane=CrystalPlane(MillerIndex(np.array(hkl), phase=phase), phase=phase),
            label_indices=tuple(int(value) for value in hkl),
        )
        for hkl in plane_hkls
    ) + tuple(_coerce_plane_overlay(overlay, phase=phase) for overlay in plane_overlays)
    planes: list[CrystalPlaneGlyph] = []
    for plane_overlay in merged_plane_overlays:
        display_indices = plane_overlay.label_indices or _default_plane_indices(plane_overlay.plane)
        polygon = _plane_polygon_for_box(
            phase,
            tuple(int(value) for value in plane_overlay.plane.miller.indices),
            repeats,
            region=prism_region_for(phase, repeats) if as_prism else None,
            # `None` asks for the family member with the largest cross-section
            # through the box, which is what a reader means by "the (110) plane
            # of this cell". An explicit overlay offset still wins.
            offset=(
                None if plane_overlay.offset is None else float(plane_overlay.offset)
            ),
        )
        if polygon is None:
            continue
        planes.append(
            CrystalPlaneGlyph(
                vertices_angstrom=polygon,
                normal_angstrom=plane_overlay.plane.normal,
                color=plane_overlay.color or crystal_style["plane_color"],
                alpha=float(
                    crystal_style["plane_alpha"]
                    if plane_overlay.alpha is None
                    else plane_overlay.alpha
                ),
                label=plane_overlay.label or format_plane_indices(display_indices),
                annotation_style=plane_overlay.annotation_style
                or PlaneAnnotationStyle(
                    color=plane_overlay.color or crystal_style["plane_color"],
                    fontsize=float(crystal_style["plane_label_fontsize"]),
                    offset_fraction=float(crystal_style["plane_label_offset_fraction"]),
                ),
            )
        )
    directions: list[CrystalDirectionGlyph] = []
    direct_basis = phase.lattice.direct_basis().matrix
    for direction_overlay in tuple(
        _coerce_direction_overlay(item, phase=phase) for item in direction_overlays
    ):
        max_repeats = np.array(repeats, dtype=np.float64)
        if np.any(direction_overlay.anchor_fractional < -1e-9) or np.any(
            direction_overlay.anchor_fractional > max_repeats + 1e-9
        ):
            raise ValueError(
                "CrystalDirectionOverlay.anchor_fractional must lie within the "
                "repeated cell volume."
            )
        endpoint_fractional = _direction_endpoint_fractional(
            direction_overlay.direction,
            anchor_fractional=direction_overlay.anchor_fractional,
            repeats=repeats,
            arrow_length_scale=direction_overlay.arrow_length_scale,
        )
        direction_display_indices: tuple[int, ...] | None = (
            direction_overlay.label_indices
            or _default_direction_indices(direction_overlay.direction)
        )
        directions.append(
            CrystalDirectionGlyph(
                start_angstrom=direct_basis @ direction_overlay.anchor_fractional,
                end_angstrom=direct_basis @ endpoint_fractional,
                color=direction_overlay.color or crystal_style["direction_color"],
                alpha=float(
                    crystal_style["direction_alpha"]
                    if direction_overlay.alpha is None
                    else direction_overlay.alpha
                ),
                label=direction_overlay.label
                or (
                    format_direction_indices(direction_display_indices)
                    if direction_display_indices is not None
                    else ""
                ),
                annotation_style=direction_overlay.annotation_style
                or DirectionAnnotationStyle(
                    color=direction_overlay.color or crystal_style["direction_color"],
                    fontsize=float(crystal_style["direction_label_fontsize"]),
                    offset_fraction=float(crystal_style["direction_label_offset_fraction"]),
                ),
            )
        )
    if site_vectors:
        # VESTA vector convention: every periodic copy of a labelled site
        # carries the same arrow (magnetic moment, displacement, force),
        # drawn from the atom center in crystal Cartesian angstrom.
        vector_color = str(crystal_style.get("site_vector_color", "#b91c1c"))
        vector_style = DirectionAnnotationStyle(
            color=vector_color,
            fontsize=float(crystal_style["direction_label_fontsize"]),
            offset_fraction=float(crystal_style["direction_label_offset_fraction"]),
        )
        vectors_by_label = {
            str(label): as_float_array(vector, shape=(3,))
            for label, vector in site_vectors.items()
        }
        unknown_labels = set(vectors_by_label) - {site.label for site, _ in atom_data}
        if unknown_labels:
            raise ValueError(
                f"site_vectors labels not present in the displayed cell: {sorted(unknown_labels)!r}"
            )
        for site, position in atom_data:
            vector = vectors_by_label.get(site.label)
            if vector is None:
                continue
            directions.append(
                CrystalDirectionGlyph(
                    start_angstrom=position,
                    end_angstrom=position + vector,
                    color=vector_color,
                    alpha=float(crystal_style["direction_alpha"]),
                    label="",
                    annotation_style=vector_style,
                )
            )
    effective_polyhedra_species = polyhedra_species
    if render_style == "polyhedral" and not effective_polyhedra_species:
        # VESTA polyhedral style: every species that achieves >= 4-coordination
        # gets its coordination polyhedron unless the caller narrows the set.
        effective_polyhedra_species = tuple(
            dict.fromkeys(site.species for site, _ in atom_data)
        )
    return CrystalScene(
        phase=phase,
        atoms=atoms,
        bonds=tuple(bonds),
        cells=tuple(cells),
        planes=tuple(planes),
        directions=tuple(directions),
        # In prism mode the prism *is* the cell being shown, so the rhombic
        # supercell box is not drawn around it: two outlines would say the
        # crystal has boundaries in two different places.
        lattice_edges=() if as_prism else _supercell_box_edges(phase, repeats),
        repeats=repeats,
        polyhedra=_coordination_polyhedra(
            atoms,
            species=effective_polyhedra_species,
            bond_tolerance_angstrom=bond_tolerance_angstrom,
            crystal_style=crystal_style,
        ),
    )


def _view_angles_from_direction(direction: np.ndarray) -> tuple[float, float]:
    vector = np.asarray(direction, dtype=np.float64)
    vector = vector / np.linalg.norm(vector)
    azim = float(np.rad2deg(np.arctan2(vector[1], vector[0])))
    elev = float(np.rad2deg(np.arcsin(vector[2])))
    return elev, azim


def _normalize_light_direction(direction: Any) -> np.ndarray:
    vector = as_float_array(direction, shape=(3,))
    norm = np.linalg.norm(vector)
    if np.isclose(norm, 0.0):
        raise ValueError("crystal.light_direction must be non-zero.")
    return vector / norm


def _view_vector_from_angles(elev_deg: float, azim_deg: float) -> np.ndarray:
    """Unit vector from the scene toward the matplotlib camera."""

    elev = np.deg2rad(elev_deg)
    azim = np.deg2rad(azim_deg)
    return np.array(
        [np.cos(elev) * np.cos(azim), np.cos(elev) * np.sin(azim), np.sin(elev)],
        dtype=np.float64,
    )


@lru_cache(maxsize=8)
def _unit_sphere_quads(resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """Quad faces and outward face normals of a unit sphere at the origin.

    Cached per resolution: every atom reuses the same unit mesh, scaled and
    translated by the renderer. Faces have shape ``(n_faces, 4, 3)``, normals
    ``(n_faces, 3)``.
    """

    if resolution < 4:
        raise ValueError("sphere mesh resolution must be at least 4.")
    u = np.linspace(0.0, 2.0 * np.pi, resolution)
    v = np.linspace(0.0, np.pi, resolution)
    uu, vv = np.meshgrid(u, v, indexing="xy")
    vertices = np.stack(
        [np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv)], axis=-1
    )
    quads = np.stack(
        [vertices[:-1, :-1], vertices[:-1, 1:], vertices[1:, 1:], vertices[1:, :-1]],
        axis=2,
    ).reshape(-1, 4, 3)
    centers = quads.mean(axis=1)
    lengths = np.linalg.norm(centers, axis=1, keepdims=True)
    normals = centers / np.where(lengths == 0.0, 1.0, lengths)
    quads.setflags(write=False)
    normals.setflags(write=False)
    return quads, normals


def _sector_quad_mask(
    quads: np.ndarray,
    start_fraction: float,
    span_fraction: float,
) -> np.ndarray:
    """Boolean mask of unit-sphere quads inside an azimuthal occupancy sector.

    The sector covers ``[start, start + span)`` as fractions of one full turn
    about +z (VESTA's pie-slice convention for partially occupied sites), with
    wrap-around handled. Quads are selected by their center azimuth.
    """

    if span_fraction <= 0.0:
        return np.zeros(quads.shape[0], dtype=bool)
    if span_fraction >= 1.0 - 1e-12:
        return np.ones(quads.shape[0], dtype=bool)
    centers = quads.mean(axis=1)
    azimuth_fraction = np.mod(np.arctan2(centers[:, 1], centers[:, 0]), 2.0 * np.pi) / (
        2.0 * np.pi
    )
    offset = np.mod(azimuth_fraction - float(start_fraction), 1.0)
    return np.asarray(offset < float(span_fraction))


def _apply_depth_cue(
    faces: np.ndarray,
    colors: np.ndarray,
    *,
    elev_deg: float,
    azim_deg: float,
    strength: float,
    background: str,
) -> np.ndarray:
    """Fade mesh face colors toward the background with distance from the viewer.

    VESTA-style depth cueing ("fog") computed for the initial view direction:
    faces farther along the line of sight blend toward the background color by
    up to ``strength``. Static per render — an interactively rotated axes keeps
    the fade of the export view, which is the publication use case.
    """

    if strength <= 0.0 or faces.shape[0] == 0:
        return colors
    view = _view_vector_from_angles(elev_deg, azim_deg)
    depth = faces.mean(axis=1) @ view
    span = float(np.max(depth) - np.min(depth))
    if span <= 1e-12:
        return colors
    nearness = (depth - float(np.min(depth))) / span
    fade = float(np.clip(strength, 0.0, 1.0)) * (1.0 - nearness)
    background_rgb = np.asarray(_to_rgb(background), dtype=np.float64)
    faded = colors.copy()
    faded[:, :3] = colors[:, :3] * (1.0 - fade[:, None]) + background_rgb[None, :] * fade[:, None]
    return faded


def _cylinder_quads(
    start: np.ndarray,
    end: np.ndarray,
    radius: float,
    *,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Quad faces and outward face normals of an open cylinder start -> end."""

    axis = end - start
    length = np.linalg.norm(axis)
    if np.isclose(length, 0.0):
        raise ValueError("Bond cylinder requires distinct start and end points.")
    direction = axis / length
    trial = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(trial, direction))), 1.0, atol=1e-8):
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u_axis = np.cross(direction, trial)
    u_axis /= np.linalg.norm(u_axis)
    v_axis = np.cross(direction, u_axis)
    theta = np.linspace(0.0, 2.0 * np.pi, resolution)
    ring = np.cos(theta)[:, None] * u_axis[None, :] + np.sin(theta)[:, None] * v_axis[None, :]
    bottom = start[None, :] + radius * ring
    top = end[None, :] + radius * ring
    faces = np.stack([bottom[:-1], bottom[1:], top[1:], top[:-1]], axis=1)
    normals = ring[:-1] + ring[1:]
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / np.where(lengths == 0.0, 1.0, lengths)
    return faces, normals


def _lit_face_colors(
    color: str,
    normals: np.ndarray,
    *,
    alpha: float,
    light_direction: np.ndarray,
    view_direction: np.ndarray,
    ambient: float,
    diffuse: float,
    specular: float,
    shininess: float,
) -> np.ndarray:
    """Blinn-Phong-style RGBA per face for ``(n, 3)`` outward unit normals.

    Lambert diffuse plus a camera-aware specular highlight, matching the
    ball-and-stick look of dedicated crystal viewers. Returns ``(n, 4)``
    facecolors.
    """

    base_rgb = np.asarray(_to_rgb(color), dtype=np.float64)
    lambert = np.clip(normals @ light_direction, 0.0, 1.0)
    reflected = 2.0 * lambert[:, None] * normals - light_direction
    reflected_norm = np.linalg.norm(reflected, axis=-1, keepdims=True)
    reflected = np.divide(reflected, np.where(reflected_norm == 0.0, 1.0, reflected_norm))
    specular_term = np.clip(reflected @ view_direction, 0.0, 1.0) ** shininess
    intensity = np.clip(ambient + diffuse * lambert, 0.0, 1.35)
    rgb = np.clip(base_rgb[None, :] * intensity[:, None], 0.0, 1.0)
    rgb = np.clip(rgb + specular * specular_term[:, None], 0.0, 1.0)
    alpha_column = np.full((normals.shape[0], 1), alpha, dtype=np.float64)
    return np.concatenate([rgb, alpha_column], axis=-1)


def _zoom_to_fit(
    axes: Any,
    *,
    spans: np.ndarray,
    center: np.ndarray,
    crystal_style: dict[str, Any],
) -> None:
    """Scale the axes so the structure fills the frame it is drawn in.

    A 3D axes sizes its bounding cube to stay inside the frame under *any*
    rotation, so at a fixed publication view most of the canvas is dead space
    and the structure reads as a small object floating in white. The corners of
    the data box are projected through the live camera and the axes is zoomed by
    whatever slack that projection leaves, capped so a strongly anisotropic cell
    (a layered structure with c >> a) can never be zoomed past the frame edge.
    """

    from mpl_toolkits.mplot3d import proj3d

    limit = float(crystal_style.get("view_zoom_limit", 1.6))
    fill = float(crystal_style.get("view_fill_fraction", 0.94))
    figure = axes.figure
    figure.canvas.draw()
    half = 0.5 * np.asarray(spans, dtype=np.float64)
    offsets = np.array(list(product((-1.0, 1.0), repeat=3)), dtype=np.float64)
    corners = np.asarray(center, dtype=np.float64)[None, :] + offsets * half[None, :]
    projected_x, projected_y, _ = proj3d.proj_transform(
        corners[:, 0], corners[:, 1], corners[:, 2], axes.get_proj()
    )
    display = axes.transData.transform(np.column_stack([projected_x, projected_y]))
    width = float(display[:, 0].max() - display[:, 0].min())
    height = float(display[:, 1].max() - display[:, 1].min())
    if width <= 1e-9 or height <= 1e-9:
        return
    window = axes.get_window_extent()
    zoom = min(fill * window.width / width, fill * window.height / height)
    # Below 1.0 the structure is overflowing the frame and must shrink; above,
    # it is floating in dead space and may grow up to the configured cap.
    zoom = float(np.clip(zoom, 0.4, limit))
    if abs(zoom - 1.0) > 1e-3:
        axes.set_box_aspect(axes.get_box_aspect(), zoom=zoom)


def _draw_crystal_frame(axes: Any, scene: CrystalScene, crystal_style: dict[str, Any]) -> None:
    """Draw the lattice box edges and any unit-cell overlays for one scene."""

    _, poly3d_collection = _matplotlib()
    lattice_linewidth = float(crystal_style.get("lattice_linewidth", 1.2))
    for edge in scene.lattice_edges:
        axes.plot(
            edge[:, 0],
            edge[:, 1],
            edge[:, 2],
            color=crystal_style["lattice_color"],
            linewidth=lattice_linewidth,
        )
    for cell in scene.cells:
        if cell.faces_angstrom:
            axes.add_collection3d(
                poly3d_collection(
                    list(cell.faces_angstrom),
                    facecolors=cell.color,
                    edgecolors="none",
                    linewidths=0.0,
                    alpha=cell.face_alpha,
                )
            )
        for edge in cell.edges_angstrom:
            axes.plot(
                edge[:, 0],
                edge[:, 1],
                edge[:, 2],
                color=cell.color,
                alpha=cell.alpha,
                linewidth=cell.linewidth,
            )


def _accumulate_crystal_mesh(
    axes: Any,
    scene: CrystalScene,
    crystal_style: dict[str, Any],
    *,
    light_direction: np.ndarray,
    view_direction: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Accumulate lit atom/bond/polyhedron faces for one scene.

    Returns ``(mesh_faces, mesh_colors)`` lists that the caller concatenates into
    a single depth-sorted `Poly3DCollection`; sharing one collection across every
    scene in a composite figure is what keeps global occlusion correct. Line-mode
    bonds, scatter-mode atoms, and polyhedron edges are per-scene artists and are
    drawn directly onto ``axes`` here.
    """

    mesh_faces: list[np.ndarray] = []
    mesh_colors: list[np.ndarray] = []
    ambient = float(crystal_style["light_ambient"])
    diffuse = float(crystal_style["light_diffuse"])
    specular = float(crystal_style["light_specular"])
    bond_render_mode = str(crystal_style.get("bond_render_mode", "cylinder")).lower()
    if bond_render_mode == "line":
        for bond in scene.bonds:
            for segment_start, segment_end, segment_color in bond.half_segments():
                axes.plot(
                    [segment_start[0], segment_end[0]],
                    [segment_start[1], segment_end[1]],
                    [segment_start[2], segment_end[2]],
                    color=segment_color,
                    alpha=bond.alpha,
                    linewidth=float(crystal_style["bond_radius"]),
                )
    else:
        bond_resolution = int(crystal_style["bond_surface_resolution"])
        bond_shininess = float(crystal_style["bond_shininess"])
        bond_specular = float(crystal_style["bond_specular_strength"]) * specular
        for bond in scene.bonds:
            for segment_start, segment_end, segment_color in bond.half_segments():
                faces, normals = _cylinder_quads(
                    segment_start,
                    segment_end,
                    bond.radius_angstrom,
                    resolution=bond_resolution,
                )
                mesh_faces.append(faces)
                mesh_colors.append(
                    _lit_face_colors(
                        segment_color,
                        normals,
                        alpha=bond.alpha,
                        light_direction=light_direction,
                        view_direction=view_direction,
                        ambient=ambient,
                        diffuse=diffuse,
                        specular=bond_specular,
                        shininess=bond_shininess,
                    )
                )
    if scene.atoms:
        atom_render_mode = str(crystal_style.get("atom_render_mode", "sphere")).lower()
        if atom_render_mode == "none":
            pass  # wireframe convention: bond network only, no atom bodies
        elif atom_render_mode == "scatter":
            positions = np.vstack([atom.position_angstrom for atom in scene.atoms])
            sizes = np.array(
                [(atom.radius_angstrom * 175.0) ** 2 for atom in scene.atoms], dtype=np.float64
            )
            colors = [atom.color for atom in scene.atoms]
            axes.scatter(
                positions[:, 0],
                positions[:, 1],
                positions[:, 2],
                s=sizes,
                c=colors,
                alpha=float(crystal_style["atom_alpha"]),
                edgecolors=crystal_style["atom_edgecolor"],
                linewidths=float(crystal_style["atom_edgewidth"]),
            )
        else:
            resolution = int(crystal_style["atom_surface_resolution"])
            shininess = float(crystal_style["atom_shininess"])
            atom_specular = float(crystal_style["atom_specular_strength"]) * specular
            vacancy_color = str(crystal_style.get("vacancy_color", "#ffffff"))
            unit_quads, unit_normals = _unit_sphere_quads(resolution)

            def _append_atom_sector(
                atom: CrystalAtomGlyph, quads: np.ndarray, normals: np.ndarray, color: str
            ) -> None:
                mesh_faces.append(
                    atom.position_angstrom[None, None, :] + atom.radius_angstrom * quads
                )
                mesh_colors.append(
                    _lit_face_colors(
                        color,
                        normals,
                        alpha=float(crystal_style["atom_alpha"]),
                        light_direction=light_direction,
                        view_direction=view_direction,
                        ambient=ambient,
                        diffuse=diffuse,
                        specular=atom_specular,
                        shininess=shininess,
                    )
                )

            for atom in scene.atoms:
                if atom.is_full_sphere:
                    _append_atom_sector(atom, unit_quads, unit_normals, atom.color)
                    continue
                # VESTA occupancy pie: this glyph's species sector, plus the
                # vacancy remainder when this glyph closes an underfilled site
                occupied = _sector_quad_mask(unit_quads, atom.sector_start, atom.occupancy)
                if np.any(occupied):
                    _append_atom_sector(
                        atom, unit_quads[occupied], unit_normals[occupied], atom.color
                    )
                if atom.vacancy_fraction > 1e-9:
                    vacant = _sector_quad_mask(
                        unit_quads,
                        atom.sector_start + atom.occupancy,
                        atom.vacancy_fraction,
                    )
                    if np.any(vacant):
                        _append_atom_sector(
                            atom, unit_quads[vacant], unit_normals[vacant], vacancy_color
                        )
        label_fontsize = float(crystal_style.get("atom_label_fontsize", 10.0))
        label_color = str(crystal_style.get("atom_label_color", "#111111"))
        for atom in scene.atoms:
            if atom.label:
                axes.text(
                    atom.position_angstrom[0],
                    atom.position_angstrom[1],
                    atom.position_angstrom[2] + 1.15 * atom.radius_angstrom,
                    atom.label,
                    color=label_color,
                    fontsize=label_fontsize,
                    ha="center",
                )
    for polyhedron in scene.polyhedra:
        # triangles ride in the same depth-sorted mesh as degenerate quads
        triangles = polyhedron.triangles_angstrom
        quads = np.concatenate([triangles, triangles[:, 2:3, :]], axis=1)
        mesh_faces.append(quads)
        mesh_colors.append(
            _lit_face_colors(
                polyhedron.color,
                polyhedron.face_normals,
                alpha=polyhedron.alpha,
                light_direction=light_direction,
                view_direction=view_direction,
                ambient=ambient,
                diffuse=diffuse,
                specular=0.0,
                shininess=1.0,
            )
        )
        for triangle in triangles:
            closed = np.vstack([triangle, triangle[0]])
            axes.plot(
                closed[:, 0],
                closed[:, 1],
                closed[:, 2],
                color=polyhedron.edge_color,
                linewidth=polyhedron.edge_width,
                alpha=min(1.0, polyhedron.alpha + 0.25),
            )
    return mesh_faces, mesh_colors


def _draw_crystal_planes_and_directions(
    axes: Any,
    scene: CrystalScene,
    crystal_style: dict[str, Any],
    *,
    scene_span: float,
) -> None:
    """Draw plane patches (with labels) and direction quivers (with labels)."""

    _, poly3d_collection = _matplotlib()
    for plane in scene.planes:
        axes.add_collection3d(
            poly3d_collection(
                [plane.vertices_angstrom],
                facecolors=plane.color,
                edgecolors=plane.color,
                linewidths=0.8,
                alpha=plane.alpha,
            )
        )
        if plane.label:
            # Anchor the label on the plane's rim rather than its centroid. A
            # centroid label lands in the middle of the atom cloud, where it is
            # hidden behind spheres and collides with the label of any second
            # plane cutting through the same region.
            vertices = plane.vertices_angstrom
            center = np.mean(vertices, axis=0)
            rim = vertices[int(np.argmax(np.linalg.norm(vertices - center, axis=1)))]
            outward = rim - center
            outward_norm = float(np.linalg.norm(outward))
            outward = outward / outward_norm if outward_norm > 1e-9 else np.zeros(3)
            offset = plane.annotation_style.offset_fraction * scene_span * (
                plane.normal_angstrom + 1.6 * outward
            )
            axes.text(
                rim[0] + offset[0],
                rim[1] + offset[1],
                rim[2] + offset[2],
                plane.label,
                color=plane.annotation_style.color or plane.color,
                fontsize=plane.annotation_style.fontsize,
                ha="center",
                va="center",
                zorder=10,
            )
    direction_arrow_ratio = float(crystal_style["direction_arrow_ratio"])
    direction_linewidth = float(crystal_style["direction_linewidth"])
    for direction in scene.directions:
        vector = direction.end_angstrom - direction.start_angstrom
        axes.quiver(
            direction.start_angstrom[0],
            direction.start_angstrom[1],
            direction.start_angstrom[2],
            vector[0],
            vector[1],
            vector[2],
            color=direction.color,
            alpha=direction.alpha,
            arrow_length_ratio=direction_arrow_ratio,
            linewidth=direction_linewidth,
        )
        if direction.label:
            unit = vector / np.linalg.norm(vector)
            offset = direction.annotation_style.offset_fraction * scene_span * unit
            label_point = direction.end_angstrom + offset
            axes.text(
                label_point[0],
                label_point[1],
                label_point[2],
                direction.label,
                color=direction.annotation_style.color or direction.color,
                fontsize=direction.annotation_style.fontsize,
                ha="center",
                va="center",
                zorder=10,
            )


def plot_crystal_structure_3d(
    scene_or_phase: CrystalScene | Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    render_style: str = "ball_and_stick",
    show_bonds: bool = True,
    bond_tolerance_angstrom: float = 0.45,
    include_boundary_atoms: bool = True,
    polyhedra_species: tuple[str, ...] = (),
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
    plane_overlays: tuple[CrystalPlane | CrystalPlaneOverlay, ...] = (),
    direction_overlays: tuple[CrystalDirection | CrystalDirectionOverlay, ...] = (),
    show_unit_cells: bool = False,
    hexagonal_prism: bool = False,
    cell_overlays: tuple[CrystalCellOverlay, ...] = (),
    slab_hkl: tuple[int, int, int] | None = None,
    slab_thickness_angstrom: float | None = None,
    atom_label_mode: str = "none",
    site_vectors: Mapping[str, Any] | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    elev_deg: float = 22.0,
    azim_deg: float = 34.0,
    projection: str = "persp",
    view_direction: CrystalDirection | np.ndarray | None = None,
    view_preset: str | None = None,
    show_legend: bool = False,
    show_frame_indicator: bool = False,
    frame_indicator_loc: str = "lower left",
    ax: Any | None = None,
) -> Any:
    """Render a crystal scene (or phase) as a VESTA-class 3D figure.

    ``render_style`` switches the whole visual system in one keyword:
    ``"ball_and_stick"``, ``"space_filling"``, ``"stick"``, ``"wireframe"``, or
    ``"polyhedral"`` (see `build_crystal_scene`). ``view_preset`` selects a
    crystallographic viewing direction by name ("a", "b", or "c": look along
    that lattice vector); ``view_direction`` accepts an arbitrary vector or
    `CrystalDirection` and wins over the preset. ``show_legend`` adds a
    per-species color legend. Set the ``depth_cue_strength`` style key above
    zero for VESTA-style distance fog computed for the initial view.

    ``show_frame_indicator`` adds a small `pytex.plotting.frames` gizmo in the
    chosen corner showing where the phase's crystal axes point *in the rendered
    view*, so a reader can tell the orientation at a glance even when the camera
    is set by a crystallographic direction rather than by angles. The gizmo is
    drawn at the figure's own ``elev_deg``/``azim_deg``, so it always agrees with
    the scene. Off by default, so existing figures are unchanged.
    """

    plt, poly3d_collection = _matplotlib()
    effective_overrides = _merged_style_overrides(render_style, style_overrides)
    style = resolve_style(theme=theme, style_path=style_path, overrides=effective_overrides)
    common = style["common"]
    crystal_style = style["crystal"]
    if isinstance(scene_or_phase, CrystalScene):
        scene = scene_or_phase
    else:
        scene = build_crystal_scene(
            scene_or_phase,
            repeats=repeats,
            render_style=render_style,
            show_bonds=show_bonds,
            bond_tolerance_angstrom=bond_tolerance_angstrom,
            include_boundary_atoms=include_boundary_atoms,
            polyhedra_species=polyhedra_species,
            plane_hkls=plane_hkls,
            plane_overlays=plane_overlays,
            direction_overlays=direction_overlays,
            show_unit_cells=show_unit_cells,
            hexagonal_prism=hexagonal_prism,
            cell_overlays=cell_overlays,
            slab_hkl=slab_hkl,
            slab_thickness_angstrom=slab_thickness_angstrom,
            atom_label_mode=atom_label_mode,
            site_vectors=site_vectors,
            theme=theme,
            style_path=style_path,
            style_overrides=style_overrides,
        )
    if ax is None:
        fig = plt.figure(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=crystal_style["background"],
        )
        axes = fig.add_subplot(111, projection="3d", proj_type=projection)
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(crystal_style["background"])
    light_direction = _normalize_light_direction(crystal_style["light_direction"])
    # resolve the viewing direction FIRST so depth cueing can fade along it
    if view_direction is not None:
        if isinstance(view_direction, CrystalDirection):
            elev_deg, azim_deg = _view_angles_from_direction(view_direction.unit_vector)
        else:
            elev_deg, azim_deg = _view_angles_from_direction(
                np.asarray(view_direction, dtype=np.float64)
            )
    elif view_preset is not None:
        axis_index = {"a": 0, "b": 1, "c": 2}.get(view_preset.lower())
        if axis_index is None:
            raise ValueError("view_preset must be one of 'a', 'b', or 'c'.")
        lattice_vector = scene.phase.lattice.direct_basis().matrix[:, axis_index]
        elev_deg, azim_deg = _view_angles_from_direction(
            lattice_vector / np.linalg.norm(lattice_vector)
        )
    _draw_crystal_frame(axes, scene, crystal_style)
    # Atoms and bonds accumulate into ONE Poly3DCollection so matplotlib
    # depth-sorts every face globally: bonds correctly disappear behind
    # atoms (and vice versa) from any viewing angle, which per-artist
    # painter's ordering cannot guarantee.
    camera_direction = _view_vector_from_angles(elev_deg, azim_deg)
    mesh_faces, mesh_colors = _accumulate_crystal_mesh(
        axes,
        scene,
        crystal_style,
        light_direction=light_direction,
        view_direction=camera_direction,
    )
    if mesh_faces:
        all_faces = np.concatenate(mesh_faces, axis=0)
        all_colors = _apply_depth_cue(
            all_faces,
            np.concatenate(mesh_colors, axis=0),
            elev_deg=elev_deg,
            azim_deg=azim_deg,
            strength=float(crystal_style.get("depth_cue_strength", 0.0)),
            background=str(crystal_style["background"]),
        )
        # Each quad is stroked in its own face colour. With `edgecolors="none"`
        # matplotlib antialiases every quad against the background, and the
        # hairline gaps between neighbouring quads read as a wireframe grid
        # drawn over each sphere. Stroking the seam closed is what makes the
        # surfaces render as solid bodies.
        mesh = poly3d_collection(
            all_faces,
            facecolors=all_colors,
            edgecolors=all_colors,
            linewidths=float(crystal_style.get("mesh_seam_linewidth", 0.3)),
        )
        mesh.set_zsort("average")
        axes.add_collection3d(mesh)
    scene_span = float(np.max(scene.bounds()[1] - scene.bounds()[0]) + 1e-6)
    _draw_crystal_planes_and_directions(axes, scene, crystal_style, scene_span=scene_span)
    bounds = scene.bounds()
    center = 0.5 * (bounds[0] + bounds[1])
    extent = np.asarray(bounds[1] - bounds[0], dtype=np.float64)
    pad = 1.0 + 2.0 * float(crystal_style.get("view_padding_fraction", 0.05))
    # Per-axis limits sized to the data, with a matching box aspect. A single
    # cubic bounding box wastes most of the canvas on an elongated cell (a
    # layered structure with c >> a shrinks to a sliver); scaling the box the
    # same way as the limits keeps the units-per-inch identical on all three
    # axes, so spheres still render round.
    spans = np.maximum(extent * pad, float(np.max(extent)) * 0.15 + 1e-6)
    for setter, axis in (
        (axes.set_xlim, 0),
        (axes.set_ylim, 1),
        (axes.set_zlim, 2),
    ):
        half = 0.5 * float(spans[axis])
        setter(center[axis] - half, center[axis] + half)
    axes.set_box_aspect(tuple(spans / float(np.max(spans))))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    _zoom_to_fit(axes, spans=spans, center=center, crystal_style=crystal_style)
    if show_legend and scene.atoms:
        from matplotlib.lines import Line2D

        species_colors: dict[str, str] = {}
        for atom in scene.atoms:
            species_colors.setdefault(atom.species, atom.color)
        handles = [
            Line2D(
                [0.0],
                [0.0],
                marker="o",
                linestyle="none",
                markersize=9,
                markerfacecolor=color,
                markeredgecolor="#334155",
                markeredgewidth=0.4,
                label=species,
            )
            for species, color in sorted(species_colors.items())
        ]
        axes.legend(handles=handles, loc="upper right", framealpha=0.85)
    if bool(crystal_style.get("hide_grid", True)):
        axes.grid(False)
    pane_rgba = (*_to_rgb(crystal_style["background"]), float(crystal_style["pane_alpha"]))
    axes.xaxis.set_pane_color(pane_rgba)
    axes.yaxis.set_pane_color(pane_rgba)
    axes.zaxis.set_pane_color(pane_rgba)
    if bool(crystal_style.get("show_axes", False)):
        axes.set_xlabel(scene.phase.crystal_frame.axes[0], color=crystal_style["axis_label_color"])
        axes.set_ylabel(scene.phase.crystal_frame.axes[1], color=crystal_style["axis_label_color"])
        axes.set_zlabel(scene.phase.crystal_frame.axes[2], color=crystal_style["axis_label_color"])
        axes.tick_params(colors=crystal_style["axis_label_color"])
    else:
        axes.set_axis_off()
    if bool(crystal_style.get("show_title", True)):
        axes.set_title(f"{scene.phase.name} Crystal Structure")
    if show_frame_indicator:
        # The lattice basis columns are the a/b/c edge vectors in crystal
        # Cartesian coordinates, which is the space the scene is drawn in, so an
        # oblique cell's gizmo leans the way the cell does.
        add_frame_indicator(
            axes,
            scene.phase.crystal_frame,
            loc=frame_indicator_loc,
            basis=scene.phase.lattice.direct_basis().matrix,
            elev_deg=elev_deg,
            azim_deg=azim_deg,
            label_frame=True,
        )
    fig.tight_layout()
    return fig
