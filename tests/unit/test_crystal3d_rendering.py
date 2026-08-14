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
from pytex.plotting.crystal3d import _cylinder_quads, _lit_face_colors, _unit_sphere_quads


def test_specular_highlight_tracks_the_actual_camera_direction() -> None:
    normals = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    light = np.array([1.0, 0.0, 1.0]) / np.sqrt(2.0)
    from_z = _lit_face_colors(
        "#000000",
        normals,
        alpha=1.0,
        light_direction=light,
        view_direction=np.array([0.0, 0.0, 1.0]),
        ambient=0.0,
        diffuse=0.0,
        specular=1.0,
        shininess=2.0,
    )
    from_x = _lit_face_colors(
        "#000000",
        normals,
        alpha=1.0,
        light_direction=light,
        view_direction=np.array([1.0, 0.0, 0.0]),
        ambient=0.0,
        diffuse=0.0,
        specular=1.0,
        shininess=2.0,
    )

    assert from_z[0, 0] > from_z[1, 0]
    assert from_x[1, 0] > from_x[0, 0]


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
    scene = build_crystal_scene(make_two_species_phase(), include_boundary_atoms=False)
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
        include_boundary_atoms=False,
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
        include_boundary_atoms=False,
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


def make_octahedral_phase() -> Phase:
    """One Ti center with six O face-center neighbors: ideal octahedron."""

    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    oxygen_fractionals = (
        (0.5, 0.5, 0.0),
        (0.5, 0.5, 1.0),
        (0.5, 0.0, 0.5),
        (0.5, 1.0, 0.5),
        (0.0, 0.5, 0.5),
        (1.0, 0.5, 0.5),
    )
    sites = [
        AtomicSite(label="Ti1", species="Ti", fractional_coordinates=np.array([0.5, 0.5, 0.5]))
    ]
    sites.extend(
        AtomicSite(label=f"O{i}", species="O", fractional_coordinates=np.array(frac))
        for i, frac in enumerate(oxygen_fractionals)
    )
    return Phase(
        "octahedron_demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=UnitCell(lattice=lattice, sites=tuple(sites)),
    )


def test_coordination_polyhedron_is_octahedron_with_outward_normals() -> None:
    scene = build_crystal_scene(make_octahedral_phase(), polyhedra_species=("Ti",))
    assert len(scene.polyhedra) == 1
    polyhedron = scene.polyhedra[0]
    assert polyhedron.center_species == "Ti"
    # an octahedron's convex hull has exactly eight triangular faces
    assert polyhedron.triangles_angstrom.shape == (8, 3, 3)
    centers = polyhedron.triangles_angstrom.mean(axis=1)
    outward = np.einsum(
        "ni,ni->n", centers - polyhedron.center_angstrom, polyhedron.face_normals
    )
    assert np.all(outward > 0.0)
    assert polyhedron.color == cpk_color("Ti")


def test_polyhedra_skip_low_coordination_centers() -> None:
    scene = build_crystal_scene(make_two_species_phase(), polyhedra_species=("Na",))
    assert scene.polyhedra == ()


def test_renderer_draws_polyhedron_faces_in_unified_mesh() -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    figure = plot_crystal_structure_3d(
        make_octahedral_phase(),
        polyhedra_species=("Ti",),
        style_overrides={
            "crystal": {"atom_surface_resolution": 8, "bond_surface_resolution": 6}
        },
    )
    figure.canvas.draw()
    axes = figure.axes[0]
    meshes = [
        artist for artist in axes.collections if isinstance(artist, Poly3DCollection)
    ]
    assert len(meshes) == 1
    atom_faces = 7 * (8 - 1) ** 2
    bond_faces = 6 * 2 * (6 - 1)
    assert len(meshes[0].get_paths()) == atom_faces + bond_faces + 8


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


def make_corner_site_phase() -> Phase:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(AtomicSite(label="Na1", species="Na", fractional_coordinates=np.zeros(3)),),
    )
    return Phase(
        "corner_demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )


def test_boundary_atoms_complete_the_cell_like_vesta() -> None:
    phase = make_corner_site_phase()
    # a corner atom renders at all eight corners of a single cell
    scene = build_crystal_scene(phase, show_bonds=False)
    assert len(scene.atoms) == 8
    # and at 3 x 2 x 2 = 12 lattice points of a 2x1x1 block
    block = build_crystal_scene(phase, repeats=(2, 1, 1), show_bonds=False)
    assert len(block.atoms) == 12
    # the historical behavior remains available
    interior = build_crystal_scene(phase, show_bonds=False, include_boundary_atoms=False)
    assert len(interior.atoms) == 1
    # boundary copies landing on explicitly listed far-face sites deduplicate
    octahedron = build_crystal_scene(make_octahedral_phase(), show_bonds=False)
    assert len(octahedron.atoms) == 7


def test_view_preset_sets_crystallographic_view() -> None:
    phase = make_corner_site_phase()
    figure = plot_crystal_structure_3d(phase, view_preset="c")
    axes = figure.axes[0]
    # looking along +c of a cubic lattice: elevation 90 deg
    assert axes.elev == pytest.approx(90.0, abs=1e-6)
    figure = plot_crystal_structure_3d(phase, view_preset="a")
    assert figure.axes[0].elev == pytest.approx(0.0, abs=1e-6)
    with pytest.raises(ValueError, match="view_preset"):
        plot_crystal_structure_3d(phase, view_preset="d")


def test_show_legend_lists_species() -> None:
    figure = plot_crystal_structure_3d(make_octahedral_phase(), show_legend=True)
    legend = figure.axes[0].get_legend()
    assert legend is not None
    labels = [text.get_text() for text in legend.get_texts()]
    assert labels == ["O", "Ti"]
