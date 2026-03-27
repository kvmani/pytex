from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core._chemistry import covalent_radius_angstrom, cpk_color
from pytex.core.lattice import CrystalDirection, Phase
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
        object.__setattr__(
            self,
            "position_angstrom",
            as_float_array(self.position_angstrom, shape=(3,)),
        )


@dataclass(frozen=True, slots=True)
class CrystalBondGlyph:
    start_angstrom: np.ndarray
    end_angstrom: np.ndarray
    color: str
    alpha: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_angstrom", as_float_array(self.start_angstrom, shape=(3,)))
        object.__setattr__(self, "end_angstrom", as_float_array(self.end_angstrom, shape=(3,)))


@dataclass(frozen=True, slots=True)
class CrystalPlaneGlyph:
    vertices_angstrom: np.ndarray
    color: str
    alpha: float
    label: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "vertices_angstrom",
            as_float_array(self.vertices_angstrom, shape=(None, 3)),
        )


@dataclass(frozen=True, slots=True)
class CrystalScene:
    phase: Phase
    atoms: tuple[CrystalAtomGlyph, ...]
    bonds: tuple[CrystalBondGlyph, ...]
    planes: tuple[CrystalPlaneGlyph, ...]
    lattice_edges: tuple[np.ndarray, ...]
    repeats: tuple[int, int, int]

    def bounds(self) -> np.ndarray:
        points: list[np.ndarray] = []
        points.extend(edge for edge in self.lattice_edges)
        points.extend(atom.position_angstrom[None, :] for atom in self.atoms)
        if self.planes:
            points.extend(plane.vertices_angstrom for plane in self.planes)
        stacked = np.vstack(points) if points else np.zeros((1, 3), dtype=np.float64)
        return np.vstack([np.min(stacked, axis=0), np.max(stacked, axis=0)])


def _supercell_atom_positions(
    phase: Phase,
    repeats: tuple[int, int, int],
) -> list[tuple[str, np.ndarray]]:
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
    basis = phase.lattice.direct_basis().matrix
    corners_frac = np.array(
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
    corners = (basis @ corners_frac.T).T
    edge_pairs = (
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 4), (2, 6),
        (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7),
    )
    return tuple(np.vstack([corners[a], corners[b]]) for a, b in edge_pairs)


def _plane_polygon_for_box(
    phase: Phase,
    hkl: tuple[int, int, int],
    repeats: tuple[int, int, int],
    offset: float = 1.0,
) -> np.ndarray | None:
    basis = phase.lattice.direct_basis().matrix
    corners_frac = np.array(
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
                point_cart = basis @ point_frac
                intersections.append(point_cart)
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
    order = np.argsort(angles)
    polygon = points[order]
    polygon.setflags(write=False)
    return polygon


def build_crystal_scene(
    phase: Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    show_bonds: bool = True,
    bond_tolerance_angstrom: float = 0.45,
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
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
            radius_angstrom=covalent_radius_angstrom(species)
            * float(crystal_style["atom_radius_scale"]),
            color=cpk_color(species),
            alpha=float(crystal_style["atom_alpha"]),
        )
        for species, position in atom_data
    )
    bonds: list[CrystalBondGlyph] = []
    if show_bonds:
        for i, atom_i in enumerate(atoms):
            for atom_j in atoms[i + 1 :]:
                cutoff = atom_i.radius_angstrom + atom_j.radius_angstrom + bond_tolerance_angstrom
                if np.linalg.norm(atom_i.position_angstrom - atom_j.position_angstrom) <= cutoff:
                    bonds.append(
                        CrystalBondGlyph(
                            start_angstrom=atom_i.position_angstrom,
                            end_angstrom=atom_j.position_angstrom,
                            color=crystal_style["bond_color"],
                            alpha=float(crystal_style["bond_alpha"]),
                        )
                    )
    planes = tuple(
        CrystalPlaneGlyph(
            vertices_angstrom=polygon,
            color=crystal_style["plane_color"],
            alpha=float(crystal_style["plane_alpha"]),
            label="(" + " ".join(str(value) for value in hkl) + ")",
        )
        for hkl in plane_hkls
        for polygon in [_plane_polygon_for_box(phase, hkl, repeats)]
        if polygon is not None
    )
    return CrystalScene(
        phase=phase,
        atoms=atoms,
        bonds=tuple(bonds),
        planes=planes,
        lattice_edges=_supercell_box_edges(phase, repeats),
        repeats=repeats,
    )


def _view_angles_from_direction(direction: np.ndarray) -> tuple[float, float]:
    vector = np.asarray(direction, dtype=np.float64)
    vector = vector / np.linalg.norm(vector)
    azim = float(np.rad2deg(np.arctan2(vector[1], vector[0])))
    elev = float(np.rad2deg(np.arcsin(vector[2])))
    return elev, azim


def plot_crystal_structure_3d(
    scene_or_phase: CrystalScene | Phase,
    *,
    repeats: tuple[int, int, int] = (1, 1, 1),
    show_bonds: bool = True,
    plane_hkls: tuple[tuple[int, int, int], ...] = (),
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
    for edge in scene.lattice_edges:
        axes.plot(
            edge[:, 0],
            edge[:, 1],
            edge[:, 2],
            color=crystal_style["lattice_color"],
            linewidth=1.2,
        )
    for bond in scene.bonds:
        axes.plot(
            [bond.start_angstrom[0], bond.end_angstrom[0]],
            [bond.start_angstrom[1], bond.end_angstrom[1]],
            [bond.start_angstrom[2], bond.end_angstrom[2]],
            color=bond.color,
            alpha=bond.alpha,
            linewidth=float(crystal_style["bond_radius"]),
        )
    if scene.atoms:
        positions = np.vstack([atom.position_angstrom for atom in scene.atoms])
        sizes = np.array(
            [(atom.radius_angstrom * 175.0) ** 2 for atom in scene.atoms],
            dtype=np.float64,
        )
        colors = [atom.color for atom in scene.atoms]
        axes.scatter(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            s=sizes,
            c=colors,
            alpha=float(crystal_style["atom_alpha"]),
            edgecolors="#f8fafc",
            linewidths=0.35,
        )
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
        center = np.mean(plane.vertices_angstrom, axis=0)
        axes.text(center[0], center[1], center[2], plane.label, color=plane.color)
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
            elev_deg, azim_deg = _view_angles_from_direction(
                np.asarray(view_direction, dtype=np.float64)
            )
    axes.view_init(elev=elev_deg, azim=azim_deg)
    axes.set_xlabel(scene.phase.crystal_frame.axes[0])
    axes.set_ylabel(scene.phase.crystal_frame.axes[1])
    axes.set_zlabel(scene.phase.crystal_frame.axes[2])
    axes.set_title(f"{scene.phase.name} Crystal Structure")
    fig.tight_layout()
    return fig
