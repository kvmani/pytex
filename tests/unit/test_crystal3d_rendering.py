from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    build_crystal_scene,
    plot_crystal_structure_3d,
)
from pytex.core._chemistry import cpk_color
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.plotting.crystal3d import _cylinder_quads, _unit_sphere_quads


def make_two_species_phase() -> Phase:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Na1", species="Na", fractional_coordinates=np.zeros(3)),
            AtomicSite(
                label="Cl1", species="Cl", fractional_coordinates=np.array([0.5, 0.0, 0.0])
            ),
        ),
    )
    return Phase(
        "halite_pair",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )


def test_bond_glyphs_default_to_two_tone_atom_colors() -> None:
    scene = build_crystal_scene(make_two_species_phase())
    assert len(scene.bonds) == 1
    bond = scene.bonds[0]
    assert bond.start_color == cpk_color("Na")
    assert bond.end_color == cpk_color("Cl")
    segments = bond.half_segments()
    assert len(segments) == 2
    midpoint = 0.5 * (bond.start_angstrom + bond.end_angstrom)
    np.testing.assert_allclose(segments[0][1], midpoint)
    np.testing.assert_allclose(segments[1][0], midpoint)
    assert segments[0][2] == cpk_color("Na")
    assert segments[1][2] == cpk_color("Cl")


def test_uniform_bond_mode_keeps_single_segment() -> None:
    scene = build_crystal_scene(
        make_two_species_phase(),
        style_overrides={"crystal": {"bond_color_mode": "uniform"}},
    )
    bond = scene.bonds[0]
    assert bond.start_color is None and bond.end_color is None
    segments = bond.half_segments()
    assert len(segments) == 1
    assert segments[0][2] == bond.color


def test_renderer_merges_atoms_and_bonds_into_one_depth_sorted_mesh() -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    phase = make_two_species_phase()
    atom_resolution = 16
    bond_resolution = 12
    figure = plot_crystal_structure_3d(
        phase,
        style_overrides={
            "crystal": {
                "atom_surface_resolution": atom_resolution,
                "bond_surface_resolution": bond_resolution,
            }
        },
    )
    figure.canvas.draw()
    axes = figure.axes[0]
    meshes = [
        artist for artist in axes.collections if isinstance(artist, Poly3DCollection)
    ]
    # one unified collection carrying every atom and bond face
    assert len(meshes) == 1
    expected_faces = 2 * (atom_resolution - 1) ** 2 + 2 * (bond_resolution - 1)
    assert len(meshes[0].get_paths()) == expected_faces


def test_unit_sphere_quads_are_cached_unit_normals() -> None:
    quads_a, normals_a = _unit_sphere_quads(16)
    quads_b, normals_b = _unit_sphere_quads(16)
    assert quads_a is quads_b and normals_a is normals_b
    np.testing.assert_allclose(np.linalg.norm(normals_a, axis=1), 1.0, atol=1e-12)
    assert quads_a.shape == ((16 - 1) ** 2, 4, 3)
    with pytest.raises(ValueError, match="at least 4"):
        _unit_sphere_quads(3)


def test_cylinder_quads_geometry_and_degenerate_rejection() -> None:
    start = np.zeros(3)
    end = np.array([0.0, 0.0, 2.0])
    faces, normals = _cylinder_quads(start, end, 0.3, resolution=12)
    assert faces.shape == (11, 4, 3)
    # outward normals are perpendicular to the cylinder axis
    np.testing.assert_allclose(normals @ np.array([0.0, 0.0, 1.0]), 0.0, atol=1e-12)
    radii = np.linalg.norm(faces[..., :2], axis=-1)
    np.testing.assert_allclose(radii, 0.3, atol=1e-12)
    with pytest.raises(ValueError, match="distinct"):
        _cylinder_quads(start, start, 0.3, resolution=12)
