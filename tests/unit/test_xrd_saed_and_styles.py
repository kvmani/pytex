from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
from matplotlib.figure import Figure

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
    build_crystal_scene,
    generate_saed_pattern,
    generate_xrd_pattern,
    list_style_themes,
    plot_crystal_structure_3d,
    plot_saed_pattern,
    plot_xrd_pattern,
    read_style_yaml,
    resolve_style,
)
from pytex.core.lattice import AtomicSite, UnitCell


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


def test_crystal_scene_slab_filter_reduces_visible_atoms() -> None:
    phase = make_phase()
    scene = build_crystal_scene(
        phase,
        repeats=(3, 1, 1),
        slab_hkl=(1, 0, 0),
        slab_thickness_angstrom=1.1,
    )
    assert len(scene.atoms) < 3
