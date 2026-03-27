from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from matplotlib.colors import to_rgb, to_rgba

from pytex.core._arrays import as_float_array
from pytex.core._chemistry import covalent_radius_angstrom, cpk_color
from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex, Phase
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.plotting.styles import resolve_style


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, Poly3DCollection


@dataclass(frozen=True, slots=True)
class CrystalAtomGlyph:
    position_angstrom: np.ndarray
    species: str
    radius_angstrom: float
    color: str
    alpha: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "position_angstrom", as_float_array(self.position_angstrom, shape=(3,)))


@dataclass(frozen=True, slots=True)
class CrystalBondGlyph:
    start_angstrom: np.ndarray
    end_angstrom: np.ndarray
    color: str
    alpha: float
    radius_angstrom: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_angstrom", as_float_array(self.start_angstrom, shape=(3,)))
        object.__setattr__(self, "end_angstrom", as_float_array(self.end_angstrom, shape=(3,)))


@dataclass(frozen=True, slots=True)
class CrystalCellOverlay:
    kind: str = "parallelepiped"
    anchor_fractional: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    span_cells: tuple[int, int, int] = (1, 1, 1)
    color: str | None = None
    alpha: float | None = None
    linewidth: float | None = None
    show_faces: bool = False
    face_alpha: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "anchor_fractional", as_float_array(self.anchor_fractional, shape=(3,)))
        if self.kind not in {"parallelepiped", "hexagonal_prism"}:
            raise ValueError("CrystalCellOverlay.kind must be either 'parallelepiped' or 'hexagonal_prism'.")
        if any(value <= 0 for value in self.span_cells):
            raise ValueError("CrystalCellOverlay.span_cells must contain strictly positive integers.")


@dataclass(frozen=True, slots=True)
class PlaneAnnotationStyle:
    color: str | None = None
    fontsize: float = 11.0
    offset_fraction: float = 0.035


@dataclass(frozen=True, slots=True)
class DirectionAnnotationStyle:
    color: str | None = None
    fontsize: float = 11.0
    offset_fraction: float = 0.03


@dataclass(frozen=True, slots=True)
class CrystalPlaneOverlay:
    plane: CrystalPlane
    offset: float | None = None
    color: str | None = None
    alpha: float | None = None
    label: str | None = None
    label_indices: tuple[int, ...] | None = None
    annotation_style: PlaneAnnotationStyle | None = None


@dataclass(frozen=True, slots=True)
class CrystalDirectionOverlay:
    direction: CrystalDirection
    anchor_fractional: np.ndarray
    color: str | None = None
    alpha: float | None = None
    label: str | None = None
    label_indices: tuple[int, ...] | None = None
    annotation_style: DirectionAnnotationStyle | None = None
    arrow_length_scale: float = 0.92

    def __post_init__(self) -> None:
        object.__setattr__(self, "anchor_fractional", as_float_array(self.anchor_fractional, shape=(3,)))
        if not 0.0 < self.arrow_length_scale <= 1.0:
            raise ValueError("CrystalDirectionOverlay.arrow_length_scale must lie in the interval (0, 1].")


@dataclass(frozen=True, slots=True)
class CrystalPlaneGlyph:
    vertices_angstrom: np.ndarray
    normal_angstrom: np.ndarray
    color: str
    alpha: float
    label: str
    annotation_style: PlaneAnnotationStyle

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices_angstrom", as_float_array(self.vertices_angstrom, shape=(None, 3)))
        object.__setattr__(self, "normal_angstrom", as_float_array(self.normal_angstrom, shape=(3,)))


@dataclass(frozen=True, slots=True)
class CrystalDirectionGlyph:
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
    kind: str
    edges_angstrom: tuple[np.ndarray, ...]
    faces_angstrom: tuple[np.ndarray, ...]
    color: str
    alpha: float
    face_alpha: float
    linewidth: float

    def __post_init__(self) -> None:
        normalized_edges = tuple(as_float_array(edge, shape=(2, 3)) for edge in self.edges_angstrom)
        normalized_faces = tuple(as_float_array(face, shape=(None, 3)) for face in self.faces_angstrom)
        object.__setattr__(self, "edges_angstrom", normalized_edges)
        object.__setattr__(self, "faces_angstrom", normalized_faces)


@dataclass(frozen=True, slots=True)
class CrystalScene:
    phase: Phase
    atoms: tuple[CrystalAtomGlyph, ...]
    bonds: tuple[CrystalBondGlyph, ...]
    cells: tuple[CrystalCellGlyph, ...]
    planes: tuple[CrystalPlaneGlyph, ...]
    directions: tuple[CrystalDirectionGlyph, ...]
    lattice_edges: tuple[np.ndarray, ...]
    repeats: tuple[int, int, int]

    def bounds(self) -> np.ndarray:
        points: list[np.ndarray] = []
        points.extend(edge for edge in self.lattice_edges)
        points.extend(edge for cell in self.cells for edge in cell.edges_angstrom)
        points.extend(atom.position_angstrom[None, :] for atom in self.atoms)
        points.extend(np.vstack([direction.start_angstrom, direction.end_angstrom]) for direction in self.directions)
        if self.planes:
            points.extend(plane.vertices_angstrom for plane in self.planes)
        stacked = np.vstack(points) if points else np.zeros((1, 3), dtype=np.float64)
        return np.vstack([np.min(stacked, axis=0), np.max(stacked, axis=0)])


def _supercell_atom_positions(phase: Phase, repeats: tuple[int, int, int]) -> list[tuple[str, np.ndarray]]:
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError("Crystal visualization requires phase.unit_cell with atomic sites.")
    direct_basis = phase.lattice.direct_basis().matrix
    atoms: list[tuple[str, np.ndarray]] = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                translation = np.array([i, j, k], dtype=np.float64)
                for site in phase.unit_cell.sites:
                    frac = site.fractional_coordinates + translation
                    position = direct_basis @ frac
                    atoms.append((site.species, position))
    return atoms


def _supercell_box_edges(phase: Phase, repeats: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    corners = _supercell_corners_cartesian(phase, repeats)
    return _polyhedron_edges_from_corners(corners)


def _polyhedron_edges_from_corners(corners: np.ndarray) -> tuple[np.ndarray, ...]:
    edge_pairs = (
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 4), (2, 6),
        (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7),
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
    return (basis @ corners_frac.T).T


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
        corners = (basis @ _parallelepiped_corners_fractional(overlay.anchor_fractional, overlay.span_cells).T).T
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
        raise ValueError("hexagonal_prism cell overlays require a lattice with a=b, alpha=beta=90 deg, gamma=120 deg.")
    basis = lattice.direct_basis().matrix
    a_vec = basis[:, 0]
    b_vec = basis[:, 1]
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
    edges = tuple(
        np.vstack([bottom[index], bottom[(index + 1) % 6]])
        for index in range(6)
    ) + tuple(
        np.vstack([top[index], top[(index + 1) % 6]])
        for index in range(6)
    ) + tuple(
        np.vstack([bottom[index], top[index]])
        for index in range(6)
    )
    faces = (
        np.vstack(bottom),
        np.vstack(top),
    ) + tuple(
        np.vstack([bottom[index], bottom[(index + 1) % 6], top[(index + 1) % 6], top[index]])
        for index in range(6)
    )
    return edges, faces


def _default_plane_offset(indices: np.ndarray) -> float:
    return 1.0 if np.any(indices > 0) else -1.0


def _plane_polygon_for_box(
    phase: Phase,
    hkl: tuple[int, int, int],
    repeats: tuple[int, int, int],
    offset: float,
) -> np.ndarray | None:
    basis = phase.lattice.direct_basis().matrix
    corners_frac = _supercell_corners_fractional(repeats)
    edge_pairs = (
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 4), (2, 6),
        (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7),
    )
    normal_frac = np.array(hkl, dtype=np.float64)
    if np.allclose(normal_frac, 0.0):
        return None
    intersections: list[np.ndarray] = []
    for a, b in edge_pairs:
        start = corners_frac[a]
        end = corners_frac[b]
        direction = end - start
        denominator = float(np.dot(normal_frac, direction))
        if np.isclose(denominator, 0.0):
            continue
        t = (offset - float(np.dot(normal_frac, start))) / denominator
        if -1e-10 <= t <= 1.0 + 1e-10:
            point_frac = start + t * direction
            max_frac = np.array(repeats, dtype=np.float64) + 1e-9
            if np.all(point_frac >= -1e-9) and np.all(point_frac <= max_frac):
                intersections.append(basis @ point_frac)
    if len(intersections) < 3:
        return None
    unique_points: list[np.ndarray] = []
    for point in intersections:
        if not any(np.allclose(point, existing, atol=1e-8) for existing in unique_points):
            unique_points.append(point)
    points = np.vstack(unique_points)
    normal_cart = phase.lattice.reciprocal_basis().matrix @ normal_frac
    normal_cart = normal_cart / np.linalg.norm(normal_cart)
    center = np.mean(points, axis=0)
    trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(trial, normal_cart))), 1.0, atol=1e-8):
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u_axis = np.cross(normal_cart, trial)
    u_axis /= np.linalg.norm(u_axis)
    v_axis = np.cross(normal_cart, u_axis)
    local = points - center
    angles = np.arctan2(local @ v_axis, local @ u_axis)
    polygon = points[np.argsort(angles)]
    polygon.setflags(write=False)
    return polygon


def _coerce_plane_overlay(overlay: CrystalPlane | CrystalPlaneOverlay, *, phase: Phase) -> CrystalPlaneOverlay:
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
    return CrystalDirectionOverlay(direction=overlay, anchor_fractional=np.zeros(3, dtype=np.float64))


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
    for anchor_value, vector_value, repeat in zip(anchor_fractional, fractional_vector, repeats, strict=True):
        if vector_value > 1e-12:
            candidates.append((repeat - anchor_value) / vector_value)
        elif vector_value < -1e-12:
            candidates.append((0.0 - anchor_value) / vector_value)
    positive_candidates = [value for value in candidates if value > 1e-12]
    if not positive_candidates:
        raise ValueError("CrystalDirection overlay does not intersect the repeated cell volume.")
    return anchor_fractional + arrow_length_scale * min(positive_candidates) * fractional_vector


def _default_direction_indices(direction: CrystalDirection) -> tuple[int, ...] | None:
    rounded = np.rint(direction.coordinates).astype(np.int64)
    if np.allclose(direction.coordinates, rounded.astype(np.float64), atol=1e-8):
        return tuple(int(value) for value in rounded)
    return None


def _default_plane_indices(plane: CrystalPlane) -> tuple[int, ...]:
    return tuple(int(value) for value in plane.miller.indices)


def build_crystal_scene(
    phase: Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    show_bonds: bool = True,
    bond_tolerance_angstrom: float = 0.45,
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
    plane_overlays: tuple[CrystalPlane | CrystalPlaneOverlay, ...] = (),
    direction_overlays: tuple[CrystalDirection | CrystalDirectionOverlay, ...] = (),
    show_unit_cells: bool = False,
    cell_overlays: tuple[CrystalCellOverlay, ...] = (),
    slab_hkl: tuple[int, int, int] | None = None,
    slab_thickness_angstrom: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
) -> CrystalScene:
    if any(value <= 0 for value in repeats):
        raise ValueError("repeats must contain strictly positive integers.")
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    crystal_style = style["crystal"]
    atom_data = _supercell_atom_positions(phase, repeats)
    if slab_hkl is not None and slab_thickness_angstrom is not None:
        normal = phase.lattice.reciprocal_basis().matrix @ np.array(slab_hkl, dtype=np.float64)
        normal = normal / np.linalg.norm(normal)
        filtered: list[tuple[str, np.ndarray]] = []
        for species, position in atom_data:
            distance = abs(float(np.dot(position, normal)) - 1.0 / np.linalg.norm(normal))
            if distance <= slab_thickness_angstrom:
                filtered.append((species, position))
        atom_data = filtered
    atoms = tuple(
        CrystalAtomGlyph(
            position_angstrom=position,
            species=species,
            radius_angstrom=covalent_radius_angstrom(species) * float(crystal_style["atom_radius_scale"]),
            color=cpk_color(species),
            alpha=float(crystal_style["atom_alpha"]),
        )
        for species, position in atom_data
    )
    bonds: list[CrystalBondGlyph] = []
    if show_bonds:
        for i, atom_i in enumerate(atoms):
            for atom_j in atoms[i + 1 :]:
                cutoff = (
                    covalent_radius_angstrom(atom_i.species)
                    + covalent_radius_angstrom(atom_j.species)
                    + bond_tolerance_angstrom
                )
                if np.linalg.norm(atom_i.position_angstrom - atom_j.position_angstrom) <= cutoff:
                    bonds.append(
                        CrystalBondGlyph(
                            start_angstrom=atom_i.position_angstrom,
                            end_angstrom=atom_j.position_angstrom,
                            color=crystal_style["bond_color"],
                            alpha=float(crystal_style["bond_alpha"]),
                            radius_angstrom=float(crystal_style["bond_radius_scale"])
                            * min(atom_i.radius_angstrom, atom_j.radius_angstrom),
                        )
                    )
    merged_cell_overlays = (
        tuple(
            CrystalCellOverlay(anchor_fractional=np.array([i, j, k], dtype=np.float64))
            for i in range(repeats[0])
            for j in range(repeats[1])
            for k in range(repeats[2])
        )
        if show_unit_cells
        else ()
    ) + tuple(_coerce_cell_overlay(overlay, phase=phase) for overlay in cell_overlays)
    cells: list[CrystalCellGlyph] = []
    for overlay in merged_cell_overlays:
        edges, faces = _cell_overlay_cartesian(phase, overlay)
        cells.append(
            CrystalCellGlyph(
                kind=overlay.kind,
                edges_angstrom=edges,
                faces_angstrom=faces if overlay.show_faces else (),
                color=overlay.color or crystal_style["cell_color"],
                alpha=float(crystal_style["cell_alpha"] if overlay.alpha is None else overlay.alpha),
                face_alpha=float(crystal_style["cell_face_alpha"] if overlay.face_alpha is None else overlay.face_alpha),
                linewidth=float(crystal_style["cell_linewidth"] if overlay.linewidth is None else overlay.linewidth),
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
    for overlay in merged_plane_overlays:
        display_indices = overlay.label_indices or _default_plane_indices(overlay.plane)
        polygon = _plane_polygon_for_box(
            phase,
            tuple(int(value) for value in overlay.plane.miller.indices),
            repeats,
            offset=_default_plane_offset(overlay.plane.miller.indices) if overlay.offset is None else float(overlay.offset),
        )
        if polygon is None:
            continue
        planes.append(
            CrystalPlaneGlyph(
                vertices_angstrom=polygon,
                normal_angstrom=overlay.plane.normal,
                color=overlay.color or crystal_style["plane_color"],
                alpha=float(crystal_style["plane_alpha"] if overlay.alpha is None else overlay.alpha),
                label=overlay.label or format_plane_indices(display_indices),
                annotation_style=overlay.annotation_style or PlaneAnnotationStyle(
                    color=overlay.color or crystal_style["plane_color"],
                    fontsize=float(crystal_style["plane_label_fontsize"]),
                    offset_fraction=float(crystal_style["plane_label_offset_fraction"]),
                ),
            )
        )
    directions: list[CrystalDirectionGlyph] = []
    direct_basis = phase.lattice.direct_basis().matrix
    for overlay in tuple(_coerce_direction_overlay(item, phase=phase) for item in direction_overlays):
        max_repeats = np.array(repeats, dtype=np.float64)
        if np.any(overlay.anchor_fractional < -1e-9) or np.any(overlay.anchor_fractional > max_repeats + 1e-9):
            raise ValueError("CrystalDirectionOverlay.anchor_fractional must lie within the repeated cell volume.")
        endpoint_fractional = _direction_endpoint_fractional(
            overlay.direction,
            anchor_fractional=overlay.anchor_fractional,
            repeats=repeats,
            arrow_length_scale=overlay.arrow_length_scale,
        )
        display_indices = overlay.label_indices or _default_direction_indices(overlay.direction)
        directions.append(
            CrystalDirectionGlyph(
                start_angstrom=direct_basis @ overlay.anchor_fractional,
                end_angstrom=direct_basis @ endpoint_fractional,
                color=overlay.color or crystal_style["direction_color"],
                alpha=float(crystal_style["direction_alpha"] if overlay.alpha is None else overlay.alpha),
                label=overlay.label or (
                    format_direction_indices(display_indices) if display_indices is not None else ""
                ),
                annotation_style=overlay.annotation_style or DirectionAnnotationStyle(
                    color=overlay.color or crystal_style["direction_color"],
                    fontsize=float(crystal_style["direction_label_fontsize"]),
                    offset_fraction=float(crystal_style["direction_label_offset_fraction"]),
                ),
            )
        )
    return CrystalScene(
        phase=phase,
        atoms=atoms,
        bonds=tuple(bonds),
        cells=tuple(cells),
        planes=tuple(planes),
        directions=tuple(directions),
        lattice_edges=_supercell_box_edges(phase, repeats),
        repeats=repeats,
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


def _atom_surface_mesh(center: np.ndarray, radius: float, *, resolution: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u = np.linspace(0.0, 2.0 * np.pi, resolution)
    v = np.linspace(0.0, np.pi, resolution)
    uu, vv = np.meshgrid(u, v, indexing="xy")
    x = center[0] + radius * np.cos(uu) * np.sin(vv)
    y = center[1] + radius * np.sin(uu) * np.sin(vv)
    z = center[2] + radius * np.cos(vv)
    normals = np.stack(
        [
            np.cos(uu) * np.sin(vv),
            np.sin(uu) * np.sin(vv),
            np.cos(vv),
        ],
        axis=-1,
    )
    return x, y, z, normals


def _cylinder_surface_mesh(
    start: np.ndarray,
    end: np.ndarray,
    radius: float,
    *,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    t = np.linspace(0.0, length, 2)
    theta_grid, t_grid = np.meshgrid(theta, t, indexing="xy")
    circle = radius * (
        np.cos(theta_grid)[..., None] * u_axis[None, None, :]
        + np.sin(theta_grid)[..., None] * v_axis[None, None, :]
    )
    points = start[None, None, :] + t_grid[..., None] * direction[None, None, :] + circle
    normals = circle / radius
    return points[..., 0], points[..., 1], points[..., 2], normals


def _lit_surface_facecolors(
    color: str,
    normals: np.ndarray,
    *,
    alpha: float,
    light_direction: np.ndarray,
    ambient: float,
    diffuse: float,
    specular: float,
    shininess: float,
) -> np.ndarray:
    base_rgb = np.asarray(to_rgb(color), dtype=np.float64)
    lambert = np.clip(normals @ light_direction, 0.0, 1.0)
    view_direction = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    reflected = 2.0 * lambert[..., None] * normals - light_direction
    reflected_norm = np.linalg.norm(reflected, axis=-1, keepdims=True)
    reflected = np.divide(reflected, np.where(reflected_norm == 0.0, 1.0, reflected_norm))
    specular_term = np.clip(reflected @ view_direction, 0.0, 1.0) ** shininess
    intensity = np.clip(ambient + diffuse * lambert, 0.0, 1.35)
    color_grid = np.clip(base_rgb[None, None, :] * intensity[..., None], 0.0, 1.0)
    color_grid = np.clip(color_grid + specular * specular_term[..., None], 0.0, 1.0)
    alpha_grid = np.full((*normals.shape[:2], 1), alpha, dtype=np.float64)
    return np.concatenate([color_grid, alpha_grid], axis=-1)


def plot_crystal_structure_3d(
    scene_or_phase: CrystalScene | Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    show_bonds: bool = True,
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
    plane_overlays: tuple[CrystalPlane | CrystalPlaneOverlay, ...] = (),
    direction_overlays: tuple[CrystalDirection | CrystalDirectionOverlay, ...] = (),
    show_unit_cells: bool = False,
    cell_overlays: tuple[CrystalCellOverlay, ...] = (),
    slab_hkl: tuple[int, int, int] | None = None,
    slab_thickness_angstrom: float | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    elev_deg: float = 22.0,
    azim_deg: float = 34.0,
    projection: str = "persp",
    view_direction: CrystalDirection | np.ndarray | None = None,
    ax: Any | None = None,
) -> Any:
    plt, poly3d_collection = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    crystal_style = style["crystal"]
    if isinstance(scene_or_phase, CrystalScene):
        scene = scene_or_phase
    else:
        scene = build_crystal_scene(
            scene_or_phase,
            repeats=repeats,
            show_bonds=show_bonds,
            plane_hkls=plane_hkls,
            plane_overlays=plane_overlays,
            direction_overlays=direction_overlays,
            show_unit_cells=show_unit_cells,
            cell_overlays=cell_overlays,
            slab_hkl=slab_hkl,
            slab_thickness_angstrom=slab_thickness_angstrom,
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
    for edge in scene.lattice_edges:
        axes.plot(edge[:, 0], edge[:, 1], edge[:, 2], color=crystal_style["lattice_color"], linewidth=1.2)
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
    bond_render_mode = str(crystal_style.get("bond_render_mode", "cylinder")).lower()
    if bond_render_mode == "line":
        for bond in scene.bonds:
            axes.plot(
                [bond.start_angstrom[0], bond.end_angstrom[0]],
                [bond.start_angstrom[1], bond.end_angstrom[1]],
                [bond.start_angstrom[2], bond.end_angstrom[2]],
                color=bond.color,
                alpha=bond.alpha,
                linewidth=float(crystal_style["bond_radius"]),
            )
    else:
        bond_resolution = int(crystal_style["bond_surface_resolution"])
        bond_shininess = float(crystal_style["bond_shininess"])
        bond_specular = float(crystal_style["bond_specular_strength"]) * float(crystal_style["light_specular"])
        for bond in scene.bonds:
            x, y, z, normals = _cylinder_surface_mesh(
                bond.start_angstrom,
                bond.end_angstrom,
                bond.radius_angstrom,
                resolution=bond_resolution,
            )
            facecolors = _lit_surface_facecolors(
                bond.color,
                normals,
                alpha=bond.alpha,
                light_direction=light_direction,
                ambient=float(crystal_style["light_ambient"]),
                diffuse=float(crystal_style["light_diffuse"]),
                specular=bond_specular,
                shininess=bond_shininess,
            )
            axes.plot_surface(
                x,
                y,
                z,
                rstride=1,
                cstride=1,
                facecolors=facecolors,
                linewidth=0.0,
                antialiased=True,
                shade=False,
            )
    if scene.atoms:
        atom_render_mode = str(crystal_style.get("atom_render_mode", "sphere")).lower()
        if atom_render_mode == "scatter":
            positions = np.vstack([atom.position_angstrom for atom in scene.atoms])
            sizes = np.array([(atom.radius_angstrom * 175.0) ** 2 for atom in scene.atoms], dtype=np.float64)
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
            ambient = float(crystal_style["light_ambient"])
            diffuse = float(crystal_style["light_diffuse"])
            specular = float(crystal_style["light_specular"])
            shininess = float(crystal_style["atom_shininess"])
            edge_rgba = to_rgba(crystal_style["atom_edgecolor"], alpha=float(crystal_style["atom_alpha"]))
            for atom in scene.atoms:
                x, y, z, normals = _atom_surface_mesh(atom.position_angstrom, atom.radius_angstrom, resolution=resolution)
                facecolors = _lit_surface_facecolors(
                    atom.color,
                    normals,
                    alpha=float(crystal_style["atom_alpha"]),
                    light_direction=light_direction,
                    ambient=ambient,
                    diffuse=diffuse,
                    specular=float(crystal_style["atom_specular_strength"]) * specular,
                    shininess=shininess,
                )
                axes.plot_surface(
                    x,
                    y,
                    z,
                    rstride=1,
                    cstride=1,
                    facecolors=facecolors,
                    linewidth=float(crystal_style["atom_edgewidth"]),
                    edgecolor=edge_rgba,
                    antialiased=True,
                    shade=False,
                )
    scene_span = float(np.max(scene.bounds()[1] - scene.bounds()[0]) + 1e-6)
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
            center = np.mean(plane.vertices_angstrom, axis=0)
            offset = plane.annotation_style.offset_fraction * scene_span * plane.normal_angstrom
            axes.text(
                center[0] + offset[0],
                center[1] + offset[1],
                center[2] + offset[2],
                plane.label,
                color=plane.annotation_style.color or plane.color,
                fontsize=plane.annotation_style.fontsize,
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
            )
    bounds = scene.bounds()
    center = 0.5 * (bounds[0] + bounds[1])
    radius = 0.55 * float(np.max(bounds[1] - bounds[0]) + 1e-6)
    axes.set_xlim(center[0] - radius, center[0] + radius)
    axes.set_ylim(center[1] - radius, center[1] + radius)
    axes.set_zlim(center[2] - radius, center[2] + radius)
    if view_direction is not None:
        if isinstance(view_direction, CrystalDirection):
            elev_deg, azim_deg = _view_angles_from_direction(view_direction.unit_vector)
        else:
            elev_deg, azim_deg = _view_angles_from_direction(np.asarray(view_direction, dtype=np.float64))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    if bool(crystal_style.get("hide_grid", True)):
        axes.grid(False)
    pane_rgba = (*to_rgb(crystal_style["background"]), float(crystal_style["pane_alpha"]))
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
    fig.tight_layout()
    return fig
