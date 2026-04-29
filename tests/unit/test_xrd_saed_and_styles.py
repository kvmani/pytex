from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.figure import Figure

from pytex import (
    CrystalCellOverlay,
    CrystalDirection,
    CrystalDirectionOverlay,
    CrystalPlane,
    CrystalPlaneOverlay,
    DirectionAnnotationStyle,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    PlaneAnnotationStyle,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
    build_crystal_scene,
    format_direction_indices,
    format_plane_indices,
    generate_saed_pattern,
    generate_xrd_pattern,
    list_style_themes,
    plot_crystal_structure_3d,
    plot_saed_pattern,
    plot_xrd_pattern,
    read_style_yaml,
    resolve_style,
)
from pytex.core.lattice import AtomicSite, MillerIndex, UnitCell


def make_phase() -> Phase:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Si1", species="Si", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
        ),
    )
    return Phase(
        "demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )


def test_style_system_lists_themes_and_merges_yaml(tmp_path: Path) -> None:
    assert {"base", "journal", "presentation", "dark"} <= set(list_style_themes())
    override_path = tmp_path / "custom.yaml"
    override_path.write_text(
        "xrd:\n  annotate_peaks: false\ncommon:\n  font:\n    size: 14\n",
        encoding="utf-8",
    )
    loaded = read_style_yaml(override_path)
    assert loaded["xrd"]["annotate_peaks"] is False
    resolved = resolve_style(theme="journal", style_path=override_path)
    assert resolved["common"]["font"]["size"] == 14
    assert resolved["xrd"]["annotate_peaks"] is False


def test_style_system_rejects_unknown_top_level_sections(tmp_path: Path) -> None:
    override_path = tmp_path / "invalid.yaml"
    override_path.write_text("unknown:\n  color: '#000000'\n", encoding="utf-8")
    try:
        read_style_yaml(override_path)
    except ValueError as exc:
        assert "Unknown top-level style sections" in str(exc)
    else:  # pragma: no cover - defensive branch
        raise AssertionError("Expected invalid style YAML to raise ValueError.")


def test_generate_xrd_pattern_contains_expected_reflection() -> None:
    phase = make_phase()
    pattern = generate_xrd_pattern(
        phase,
        radiation=RadiationSpec.cu_ka(),
        two_theta_range_deg=(20.0, 80.0),
        resolution_deg=0.02,
        max_index=4,
    )
    reflection_111 = next(
        reflection
        for reflection in pattern.reflections
        if tuple(reflection.miller_indices) == (1, 1, 1)
    )
    expected_two_theta = np.rad2deg(2.0 * np.arcsin(1.5406 / (2.0 * (3.0 / np.sqrt(3.0)))))
    assert np.isclose(reflection_111.two_theta_deg, expected_two_theta, atol=1e-3)
    assert 0.0 <= np.max(pattern.intensity_grid) <= 1.0 + 1e-12


def test_generate_saed_pattern_respects_zone_axis_geometry() -> None:
    phase = make_phase()
    zone_axis = ZoneAxis(np.array([0, 0, 1]), phase=phase)
    pattern = generate_saed_pattern(
        phase,
        zone_axis,
        camera_constant_mm_angstrom=150.0,
        max_index=4,
    )
    assert pattern.spots
    for spot in pattern.spots:
        assert np.isclose(spot.reciprocal_vector_detector[2], 0.0, atol=1e-6)


def test_xrd_saed_and_crystal_plotters_return_figures() -> None:
    phase = make_phase()
    xrd_pattern = generate_xrd_pattern(phase, two_theta_range_deg=(20.0, 80.0), max_index=4)
    saed_pattern = generate_saed_pattern(
        phase,
        ZoneAxis(np.array([0, 0, 1]), phase=phase),
        max_index=4,
    )
    scene = build_crystal_scene(phase, repeats=(1, 1, 1), plane_hkls=((1, 0, 0),))
    figures = [
        plot_xrd_pattern(xrd_pattern),
        plot_saed_pattern(saed_pattern),
        plot_crystal_structure_3d(scene),
    ]
    for figure in figures:
        assert isinstance(figure, Figure)


def test_crystal_scene_contains_atoms_and_plane() -> None:
    phase = make_phase()
    scene = build_crystal_scene(phase, repeats=(2, 1, 1), plane_hkls=((1, 0, 0),))
    assert len(scene.atoms) == 2
    assert len(scene.planes) == 1


def test_crystal_scene_can_overlay_repeated_unit_cells() -> None:
    phase = make_phase()
    scene = build_crystal_scene(phase, repeats=(2, 1, 1), show_unit_cells=True)
    assert len(scene.cells) == 2
    assert all(cell.kind == "parallelepiped" for cell in scene.cells)
    assert all(len(cell.edges_angstrom) == 12 for cell in scene.cells)


def test_crystal_scene_slab_filter_reduces_visible_atoms() -> None:
    phase = make_phase()
    scene = build_crystal_scene(
        phase,
        repeats=(3, 1, 1),
        slab_hkl=(1, 0, 0),
        slab_thickness_angstrom=1.1,
    )
    assert len(scene.atoms) < 3


def test_crystal_scene_preserves_species_dependent_colors_and_sizes() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(5.6402, 5.6402, 5.6402, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Na1", species="Na", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
            AtomicSite(label="Cl1", species="Cl", fractional_coordinates=np.array([0.5, 0.5, 0.5])),
        ),
    )
    phase = Phase(
        "NaCl",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )
    scene = build_crystal_scene(phase, repeats=(1, 1, 1), show_bonds=False)
    assert len(scene.atoms) == 2
    assert scene.atoms[0].color != scene.atoms[1].color
    assert scene.atoms[0].radius_angstrom != scene.atoms[1].radius_angstrom


def test_crystal_scene_bonds_use_chemical_cutoffs_not_render_scale() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(2.4, 2.4, 2.4, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Si1", species="Si", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
            AtomicSite(label="Si2", species="Si", fractional_coordinates=np.array([0.5, 0.0, 0.0])),
        ),
    )
    phase = Phase(
        "bond_demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )
    scene = build_crystal_scene(
        phase,
        repeats=(1, 1, 1),
        show_bonds=True,
        style_overrides={"crystal": {"atom_radius_scale": 0.05}},
    )
    assert len(scene.bonds) == 1
    assert scene.bonds[0].radius_angstrom > 0.0


def test_format_miller_indices_supports_negative_components() -> None:
    assert format_direction_indices((1, 1, -2, 0)) == "$[11\\bar{2}0]$"
    assert format_plane_indices((1, 1, -2, 1)) == "$(11\\bar{2}1)$"
    assert format_plane_indices((0, 0, 0, -1), style="plain") == "(000-1)"


def test_hexagonal_miller_bravais_helpers_construct_primitives() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal)
    phase = Phase(
        "zr_hcp",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
        crystal_frame=crystal,
    )
    direction = CrystalDirection.from_miller_bravais((2, -1, -1, 0), phase=phase)
    plane = CrystalPlane.from_miller_bravais((1, 1, -2, 1), phase=phase)
    zone = ZoneAxis.from_miller_bravais((0, 0, 0, 1), phase=phase)
    assert np.allclose(direction.coordinates, np.array([1.0, 0.0, 0.0]))
    assert tuple(plane.miller.indices) == (1, 1, 1)
    assert tuple(zone.indices) == (0, 0, 1)


def test_hexagonal_prism_cell_overlay_is_constructed_for_hexagonal_axes() -> None:
    pytest.importorskip(
        "pymatgen.core",
        reason=(
            "CIF-backed crystal-scene fixture loading requires the optional "
            "pymatgen dependency."
        ),
    )
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    phase = Phase.from_cif("fixtures/phases/zr_hcp/phase.cif", crystal_frame=crystal)
    scene = build_crystal_scene(
        phase,
        repeats=(2, 2, 1),
        cell_overlays=(
            CrystalCellOverlay(
                kind="hexagonal_prism",
                anchor_fractional=np.array([1.0, 1.0, 0.0]),
                show_faces=True,
                alpha=0.75,
                face_alpha=0.12,
            ),
        ),
    )
    assert len(scene.cells) == 1
    assert scene.cells[0].kind == "hexagonal_prism"
    assert len(scene.cells[0].edges_angstrom) == 18
    assert len(scene.cells[0].faces_angstrom) == 8


def test_crystal_scene_supports_planes_and_directions_with_annotations() -> None:
    phase = make_phase()
    scene = build_crystal_scene(
        phase,
        repeats=(2, 2, 2),
        plane_overlays=(
            CrystalPlaneOverlay(
                plane=CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
                label_indices=(1, 1, 1),
                color="#ef4444",
                alpha=0.22,
                annotation_style=PlaneAnnotationStyle(fontsize=10.0),
            ),
            CrystalPlaneOverlay(
                plane=CrystalPlane(MillerIndex([0, 0, -1], phase=phase), phase=phase),
                offset=-1.0,
                label_indices=(0, 0, 0, -1),
                color="#0f766e",
                alpha=0.18,
            ),
        ),
        direction_overlays=(
            CrystalDirectionOverlay(
                direction=CrystalDirection(np.array([1.0, 1.0, 1.0]), phase=phase),
                anchor_fractional=np.array([0.0, 0.0, 0.0]),
                label_indices=(1, 1, 1),
                color="#2563eb",
                annotation_style=DirectionAnnotationStyle(fontsize=10.0),
            ),
        ),
    )
    assert len(scene.planes) == 2
    assert len(scene.directions) == 1
    assert scene.planes[0].label == "$(111)$"
    assert scene.directions[0].label == "$[111]$"


def test_crystal_plotter_supports_cylinder_bond_rendering() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(2.4, 2.4, 2.4, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite(label="Si1", species="Si", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
            AtomicSite(label="Si2", species="Si", fractional_coordinates=np.array([0.5, 0.0, 0.0])),
        ),
    )
    phase = Phase(
        "bond_demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
    )
    figure = plot_crystal_structure_3d(
        phase,
        show_bonds=True,
        style_overrides={"crystal": {"bond_render_mode": "cylinder", "atom_radius_scale": 0.05}},
    )
    assert isinstance(figure, Figure)
