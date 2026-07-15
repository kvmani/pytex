"""VESTA-parity feature tests: render styles, occupancy sectors, labels,
site vectors, depth cueing, and programmatic distance measurement.

These lock in the visual-system semantics documented in
docs/testing/vesta_parity_matrix.md.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pytex import (
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    Transform3D,
    build_crystal_scene,
    plot_crystal_structure_3d,
)
from pytex.core._chemistry import (
    atomic_radius_angstrom,
    covalent_radius_angstrom,
    display_radius_angstrom,
    van_der_waals_radius_angstrom,
)
from pytex.core.conventions import FrameDomain
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.plotting.crystal3d import _apply_depth_cue, _sector_quad_mask, _unit_sphere_quads


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def _phase_with_sites(sites: tuple[AtomicSite, ...], a: float = 4.0) -> Phase:
    crystal = _crystal_frame()
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    return Phase(
        "demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
    )


def _halite_pair_phase() -> Phase:
    return _phase_with_sites(
        (
            AtomicSite(label="Na1", species="Na", fractional_coordinates=np.zeros(3)),
            AtomicSite(
                label="Cl1", species="Cl", fractional_coordinates=np.array([0.5, 0.0, 0.0])
            ),
        )
    )


# --------------------------------------------------------------------------- #
# Display radii
# --------------------------------------------------------------------------- #


def test_display_radius_kinds_are_ordered_and_positive() -> None:
    for species in ("Na", "Fe", "O"):
        covalent = covalent_radius_angstrom(species)
        atomic = atomic_radius_angstrom(species)
        vdw = van_der_waals_radius_angstrom(species)
        assert covalent > 0.0 and atomic > 0.0 and vdw > 0.0
        # van der Waals radii are the largest of the three systems
        assert vdw > covalent
        assert display_radius_angstrom(species, kind="covalent") == covalent
        assert display_radius_angstrom(species, kind="atomic") == atomic
        assert display_radius_angstrom(species, kind="van_der_waals") == vdw


def test_display_radius_rejects_unknown_kind() -> None:
    with pytest.raises(ValueError, match="kind"):
        display_radius_angstrom("Fe", kind="metallic")


# --------------------------------------------------------------------------- #
# Render styles
# --------------------------------------------------------------------------- #


def test_space_filling_uses_atomic_radii_and_suppresses_bonds() -> None:
    scene = build_crystal_scene(
        _halite_pair_phase(), render_style="space_filling", include_boundary_atoms=False
    )
    assert scene.bonds == ()
    by_species = {atom.species: atom for atom in scene.atoms}
    assert by_species["Na"].radius_angstrom == pytest.approx(atomic_radius_angstrom("Na"))
    assert by_species["Cl"].radius_angstrom == pytest.approx(atomic_radius_angstrom("Cl"))


def test_stick_style_matches_bond_and_atom_radii() -> None:
    scene = build_crystal_scene(
        _halite_pair_phase(), render_style="stick", include_boundary_atoms=False
    )
    assert len(scene.bonds) == 1
    smallest_atom = min(atom.radius_angstrom for atom in scene.atoms)
    assert scene.bonds[0].radius_angstrom == pytest.approx(smallest_atom)


def test_wireframe_style_renders_bond_lines_only() -> None:
    figure = plot_crystal_structure_3d(
        _halite_pair_phase(), render_style="wireframe", include_boundary_atoms=False
    )
    figure.canvas.draw()
    axes = figure.axes[0]
    meshes = [artist for artist in axes.collections if isinstance(artist, Poly3DCollection)]
    # the VESTA wireframe convention: bond lines only -- no lit quad mesh
    # and no atom bodies (scatter) at all
    assert meshes == []
    assert list(axes.collections) == []
    assert len(axes.lines) > 0


def test_polyhedral_style_autoselects_species() -> None:
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
    scene = build_crystal_scene(_phase_with_sites(tuple(sites)), render_style="polyhedral")
    assert len(scene.polyhedra) >= 1
    assert any(polyhedron.center_species == "Ti" for polyhedron in scene.polyhedra)


def test_unknown_render_style_is_rejected() -> None:
    with pytest.raises(ValueError, match="render_style"):
        build_crystal_scene(_halite_pair_phase(), render_style="cartoon")


def test_user_style_overrides_beat_render_style_preset() -> None:
    scene = build_crystal_scene(
        _halite_pair_phase(),
        render_style="space_filling",
        include_boundary_atoms=False,
        style_overrides={"crystal": {"atom_radius_scale": 0.5}},
    )
    by_species = {atom.species: atom for atom in scene.atoms}
    assert by_species["Na"].radius_angstrom == pytest.approx(0.5 * atomic_radius_angstrom("Na"))


# --------------------------------------------------------------------------- #
# Occupancy pie-spheres
# --------------------------------------------------------------------------- #


def _mixed_occupancy_phase() -> Phase:
    return _phase_with_sites(
        (
            AtomicSite(
                label="Fe1",
                species="Fe",
                fractional_coordinates=np.array([0.5, 0.5, 0.5]),
                occupancy=0.5,
            ),
            AtomicSite(
                label="Ni1",
                species="Ni",
                fractional_coordinates=np.array([0.5, 0.5, 0.5]),
                occupancy=0.5,
            ),
            AtomicSite(
                label="Na1", species="Na", fractional_coordinates=np.zeros(3), occupancy=0.7
            ),
        )
    )


def test_shared_site_becomes_consecutive_sectors_with_shared_radius() -> None:
    scene = build_crystal_scene(
        _mixed_occupancy_phase(), include_boundary_atoms=False, show_bonds=False
    )
    shared = [atom for atom in scene.atoms if atom.species in {"Fe", "Ni"}]
    assert len(shared) == 2
    by_species = {atom.species: atom for atom in shared}
    assert by_species["Fe"].sector_start == pytest.approx(0.0)
    assert by_species["Fe"].occupancy == pytest.approx(0.5)
    assert by_species["Ni"].sector_start == pytest.approx(0.5)
    assert by_species["Ni"].occupancy == pytest.approx(0.5)
    # one shared sphere radius for a coherent pie
    assert shared[0].radius_angstrom == pytest.approx(shared[1].radius_angstrom)
    assert not shared[0].is_full_sphere


def test_partial_occupancy_carries_vacancy_remainder() -> None:
    scene = build_crystal_scene(
        _mixed_occupancy_phase(), include_boundary_atoms=False, show_bonds=False
    )
    partial = next(atom for atom in scene.atoms if atom.species == "Na")
    assert partial.occupancy == pytest.approx(0.7)
    assert partial.vacancy_fraction == pytest.approx(0.3)


def test_full_occupancy_site_stays_plain_sphere() -> None:
    scene = build_crystal_scene(
        _halite_pair_phase(), include_boundary_atoms=False, show_bonds=False
    )
    assert all(atom.is_full_sphere for atom in scene.atoms)


def test_sector_quad_mask_selects_expected_fraction() -> None:
    resolution = 17  # 16 azimuth columns: 0.5 selects exactly half
    quads, _ = _unit_sphere_quads(resolution)
    half = _sector_quad_mask(quads, 0.0, 0.5)
    assert int(half.sum()) == quads.shape[0] // 2
    everything = _sector_quad_mask(quads, 0.25, 1.0)
    assert bool(everything.all())
    nothing = _sector_quad_mask(quads, 0.25, 0.0)
    assert not bool(nothing.any())
    # complementary sectors partition the sphere
    other_half = _sector_quad_mask(quads, 0.5, 0.5)
    assert not np.any(half & other_half)
    assert bool((half | other_half).all())


def test_occupancy_render_smoke() -> None:
    figure = plot_crystal_structure_3d(_mixed_occupancy_phase(), show_bonds=False)
    assert len(figure.axes) == 1


def test_transformed_scene_preserves_occupancy_sectors() -> None:
    scene = build_crystal_scene(
        _mixed_occupancy_phase(), include_boundary_atoms=False, show_bonds=False
    )
    moved = scene.transformed(Transform3D.from_matrix(np.eye(3), translation=[5.0, 0.0, 0.0]))
    original = {(atom.species, atom.sector_start, atom.occupancy) for atom in scene.atoms}
    preserved = {(atom.species, atom.sector_start, atom.occupancy) for atom in moved.atoms}
    assert original == preserved


# --------------------------------------------------------------------------- #
# Atom labels and site vectors
# --------------------------------------------------------------------------- #


def test_atom_label_modes() -> None:
    phase = _halite_pair_phase()
    species_scene = build_crystal_scene(
        phase, atom_label_mode="species", include_boundary_atoms=False
    )
    assert {atom.label for atom in species_scene.atoms} == {"Na", "Cl"}
    site_scene = build_crystal_scene(phase, atom_label_mode="site", include_boundary_atoms=False)
    assert {atom.label for atom in site_scene.atoms} == {"Na1", "Cl1"}
    none_scene = build_crystal_scene(phase, include_boundary_atoms=False)
    assert all(atom.label is None for atom in none_scene.atoms)
    with pytest.raises(ValueError, match="atom_label_mode"):
        build_crystal_scene(phase, atom_label_mode="index")


def test_atom_labels_are_drawn_as_text() -> None:
    figure = plot_crystal_structure_3d(
        _halite_pair_phase(), atom_label_mode="species", include_boundary_atoms=False
    )
    texts = {text.get_text() for text in figure.axes[0].texts}
    assert {"Na", "Cl"} <= texts


def test_site_vectors_attach_to_every_periodic_copy() -> None:
    phase = _halite_pair_phase()
    moment = np.array([0.0, 0.0, 1.6])
    # boundary atoms: the corner Na site renders 8 times, each with an arrow
    scene = build_crystal_scene(phase, site_vectors={"Na1": moment})
    na_count = sum(1 for atom in scene.atoms if atom.species == "Na")
    assert len(scene.directions) == na_count == 8
    for glyph in scene.directions:
        np.testing.assert_allclose(glyph.end_angstrom - glyph.start_angstrom, moment, atol=1e-12)


def test_site_vectors_reject_unknown_labels() -> None:
    with pytest.raises(ValueError, match="site_vectors"):
        build_crystal_scene(_halite_pair_phase(), site_vectors={"Xx9": np.array([0.0, 0.0, 1.0])})


# --------------------------------------------------------------------------- #
# Depth cueing
# --------------------------------------------------------------------------- #


def test_depth_cue_fades_far_faces_and_zero_is_identity() -> None:
    faces = np.array(
        [
            [[0.0, 0.0, 0.0]] * 4,
            [[10.0, 0.0, 0.0]] * 4,
        ]
    )
    colors = np.array([[1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
    unchanged = _apply_depth_cue(
        faces, colors, elev_deg=0.0, azim_deg=0.0, strength=0.0, background="#ffffff"
    )
    np.testing.assert_allclose(unchanged, colors)
    # viewer along +x: the face at x=10 is NEAR, the face at x=0 is FAR
    faded = _apply_depth_cue(
        faces, colors, elev_deg=0.0, azim_deg=0.0, strength=0.8, background="#ffffff"
    )
    np.testing.assert_allclose(faded[1], colors[1])  # near face untouched
    assert faded[0][1] > 0.5  # far red face gained white
    np.testing.assert_allclose(faded[:, 3], 1.0)  # alpha untouched


def test_depth_cue_render_smoke() -> None:
    figure = plot_crystal_structure_3d(
        _halite_pair_phase(),
        style_overrides={"crystal": {"depth_cue_strength": 0.6}},
    )
    assert len(figure.axes) == 1


# --------------------------------------------------------------------------- #
# Measurement
# --------------------------------------------------------------------------- #


def test_bond_lengths_match_halite_geometry() -> None:
    # NaCl-type nearest neighbour distance is a/2 = 2.0 angstrom for a = 4
    scene = build_crystal_scene(_halite_pair_phase(), include_boundary_atoms=False)
    np.testing.assert_allclose(scene.bond_lengths_angstrom(), 2.0, atol=1e-12)
    summary = scene.bond_length_summary()
    assert set(summary) == {("Cl", "Na")}
    stats = summary[("Cl", "Na")]
    assert stats["count"] == 1.0
    assert stats["min"] == stats["mean"] == stats["max"] == pytest.approx(2.0)


def test_bond_length_summary_empty_scene() -> None:
    scene = build_crystal_scene(
        _halite_pair_phase(), show_bonds=False, include_boundary_atoms=False
    )
    assert scene.bond_lengths_angstrom().shape == (0,)
    assert scene.bond_length_summary() == {}
