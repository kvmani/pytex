"""Regression tests for the vectorized kinematic zone-axis engine (CD1).

Pinned references:

- Relativistic electron wavelengths: De Graef, Introduction to Conventional
  Transmission Electron Microscopy (2003), Table 2.2.
- Ni lattice parameter 3.52387 angstrom (repo fixture corpus; ICDD 04-0850),
  d(111) = a / sqrt(3) = 2.03451 angstrom.
- Legacy parity target: pytex.diffraction.saed.generate_saed_pattern.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable

import numpy as np
import pytest

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import (
    AtomicSite,
    CrystalDirection,
    Lattice,
    Phase,
    UnitCell,
    ZoneAxis,
)
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.kinematic import (
    KinematicSimulationConfig,
    SpotTable,
    centering_allowed_mask,
    double_diffraction_sums,
    electron_structure_factors,
    electron_wavelength_angstrom,
    simulate_zone_axis_spots,
    zone_basis_from_axis,
)
from pytex.diffraction.physics import ReflectionCondition
from pytex.diffraction.saed import _choose_zone_basis, generate_saed_pattern
from pytex.diffraction.shape_factors import FiniteThicknessShapeFactor


def _crystal_frame(name: str = "crystal") -> ReferenceFrame:
    return ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def make_nickel_phase() -> Phase:
    crystal = _crystal_frame()
    lattice = Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal)
    sites = tuple(
        AtomicSite(label=f"Ni{i}", species="Ni", fractional_coordinates=np.array(coords))
        for i, coords in enumerate(
            [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
        )
    )
    return Phase(
        "nickel-fcc",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group_symbol="Fm-3m",
    )


def make_iron_bcc_phase() -> Phase:
    crystal = _crystal_frame()
    lattice = Lattice(2.8665, 2.8665, 2.8665, 90.0, 90.0, 90.0, crystal_frame=crystal)
    sites = (
        AtomicSite(label="Fe1", species="Fe", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
        AtomicSite(label="Fe2", species="Fe", fractional_coordinates=np.array([0.5, 0.5, 0.5])),
    )
    return Phase(
        "iron-bcc",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group_symbol="Im-3m",
    )


class TestElectronWavelength:
    @pytest.mark.parametrize(
        ("kev", "expected_angstrom"),
        [(100.0, 0.037014), (200.0, 0.025079), (300.0, 0.019687)],
    )
    def test_pinned_literature_values(self, kev: float, expected_angstrom: float) -> None:
        assert electron_wavelength_angstrom(kev) == pytest.approx(expected_angstrom, abs=5e-7)

    def test_monotonically_decreasing_with_voltage(self) -> None:
        values = [electron_wavelength_angstrom(kv) for kv in (80.0, 120.0, 200.0, 300.0)]
        assert all(a > b for a, b in itertools.pairwise(values))

    @pytest.mark.parametrize("bad", [0.0, -100.0, float("nan"), float("inf")])
    def test_invalid_voltage_raises(self, bad: float) -> None:
        with pytest.raises(ValueError):
            electron_wavelength_angstrom(bad)


class TestZoneBasis:
    @pytest.mark.parametrize(
        "zone",
        [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.123, -0.456, 0.789],
            [1.0, 0.0, 0.0],
        ],
    )
    def test_orthonormal_right_handed(self, zone: list[float]) -> None:
        basis = zone_basis_from_axis(np.array(zone))
        assert np.allclose(basis.T @ basis, np.eye(3), atol=1e-12)
        assert np.allclose(np.cross(basis[:, 0], basis[:, 1]), basis[:, 2], atol=1e-12)
        assert np.allclose(basis[:, 2], np.array(zone) / np.linalg.norm(zone), atol=1e-12)

    def test_matches_legacy_construction(self) -> None:
        zone = np.array([0.0, 1.0, 1.0]) / np.sqrt(2.0)
        assert np.allclose(zone_basis_from_axis(zone), _choose_zone_basis(zone), atol=1e-12)

    def test_align_g_puts_vector_on_positive_u(self) -> None:
        zone = np.array([0.0, 0.0, 1.0])
        g_target = np.array([1.0, 1.0, 0.4])
        basis = zone_basis_from_axis(zone, align_g_cartesian=g_target)
        in_plane = g_target - np.dot(g_target, zone) * zone
        assert float(in_plane @ basis[:, 0]) == pytest.approx(float(np.linalg.norm(in_plane)))
        assert float(in_plane @ basis[:, 1]) == pytest.approx(0.0, abs=1e-12)

    def test_align_g_parallel_to_zone_raises(self) -> None:
        with pytest.raises(ValueError, match="parallel"):
            zone_basis_from_axis(np.array([0.0, 0.0, 1.0]), align_g_cartesian=np.array([0, 0, 2.0]))

    def test_in_plane_rotation_rotates_pattern_counterclockwise(self) -> None:
        zone = np.array([0.0, 0.0, 1.0])
        base = zone_basis_from_axis(zone)
        rotated = zone_basis_from_axis(zone, in_plane_rotation_deg=90.0)
        g = np.array([0.7, 0.2, 0.0])
        coords = g @ base[:, :2]
        coords_rotated = g @ rotated[:, :2]
        expected = np.array([-coords[1], coords[0]])
        assert np.allclose(coords_rotated, expected, atol=1e-12)

    def test_basis_is_read_only(self) -> None:
        basis = zone_basis_from_axis(np.array([0.0, 0.0, 1.0]))
        with pytest.raises(ValueError):
            basis[0, 0] = 5.0


class TestCenteringMask:
    @pytest.mark.parametrize("centering", ["P", "I", "F", "A", "B", "C", "R"])
    def test_matches_scalar_reference(self, centering: str) -> None:
        condition = ReflectionCondition(centering=centering)
        rng = np.random.default_rng(7)
        hkl = rng.integers(-4, 5, size=(200, 3))
        hkl = hkl[np.any(hkl != 0, axis=1)]
        mask = centering_allowed_mask(hkl, condition)
        expected = np.array([condition.is_allowed(row) for row in hkl])
        assert np.array_equal(mask, expected)

    def test_bad_shape_raises(self) -> None:
        with pytest.raises(ValueError, match="shape"):
            centering_allowed_mask(np.array([1, 1, 1]), ReflectionCondition())


class TestStructureFactors:
    def test_fcc_all_even_all_odd_rule(self) -> None:
        phase = make_nickel_phase()
        hkl = np.array([[1, 1, 1], [2, 0, 0], [1, 1, 0], [2, 1, 0]])
        g = np.linalg.norm(hkl.astype(float) @ phase.lattice.reciprocal_basis().matrix.T, axis=1)
        factors = np.abs(electron_structure_factors(phase, hkl, g))
        assert factors[0] == pytest.approx(4.0 * 28.0)
        assert factors[1] == pytest.approx(4.0 * 28.0)
        assert factors[2] == pytest.approx(0.0, abs=1e-10)
        assert factors[3] == pytest.approx(0.0, abs=1e-10)

    def test_phase_without_sites_gives_unit_factors(self) -> None:
        crystal = _crystal_frame()
        lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
        phase = Phase(
            "bare",
            lattice=lattice,
            symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
            crystal_frame=crystal,
        )
        hkl = np.array([[1, 0, 0], [2, 2, 0]])
        factors = electron_structure_factors(phase, hkl, np.array([0.3, 0.9]))
        assert np.allclose(factors, 1.0)

    def test_debye_waller_damps_high_g(self) -> None:
        crystal = _crystal_frame()
        lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
        site = AtomicSite(
            label="Ni1",
            species="Ni",
            fractional_coordinates=np.array([0.0, 0.0, 0.0]),
            b_iso=1.0,
        )
        phase = Phase(
            "damped",
            lattice=lattice,
            symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
            crystal_frame=crystal,
            unit_cell=UnitCell(lattice=lattice, sites=(site,)),
        )
        hkl = np.array([[1, 0, 0], [3, 0, 0]])
        g = np.array([1.0 / 3.0, 1.0])
        factors = np.abs(electron_structure_factors(phase, hkl, g))
        assert factors[1] < factors[0]


class TestSimulateZoneAxisSpots:
    def test_parity_with_legacy_saed_ni_011(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        legacy = generate_saed_pattern(phase, zone, camera_constant_mm_angstrom=180.0, max_index=3)
        config = KinematicSimulationConfig(
            beam_energy_kev=200.0,
            camera_constant_mm_angstrom=180.0,
            max_index=3,
            g_max_inv_angstrom=None,
        )
        table = simulate_zone_axis_spots(phase, zone, config=config)
        legacy_map = {
            tuple(int(v) for v in spot.miller_indices): spot.detector_coordinates
            for spot in legacy.spots
            if spot.intensity > 1e-8
        }
        new_map = {
            tuple(int(v) for v in row): table.detector_mm[i] for i, row in enumerate(table.hkl)
        }
        assert set(new_map) == set(legacy_map)
        for key, coords in new_map.items():
            assert np.allclose(coords, legacy_map[key], atol=1e-9)

    def test_fcc_forbidden_reflections_absent(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        present = {tuple(int(v) for v in row) for row in table.hkl}
        assert (2, 0, 0) in present
        assert (2, 2, 0) in present
        assert (1, 0, 0) not in present
        assert (1, 1, 0) not in present

    def test_bcc_forbidden_reflections_absent(self) -> None:
        phase = make_iron_bcc_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        present = {tuple(int(v) for v in row) for row in table.hkl}
        assert (1, 1, 0) in present
        assert (2, 0, 0) in present
        assert (1, 0, 0) not in present

    def test_d_spacing_pinned_ni_111(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, -1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        rows = {tuple(int(v) for v in row): i for i, row in enumerate(table.hkl)}
        assert (1, 1, 1) in rows
        d_111 = float(table.d_spacing_angstrom[rows[(1, 1, 1)]])
        assert d_111 == pytest.approx(3.52387 / np.sqrt(3.0), abs=1e-9)
        assert d_111 == pytest.approx(2.03451, abs=5e-6)

    def test_zolz_excitation_error_identity(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        config = KinematicSimulationConfig()
        table = simulate_zone_axis_spots(phase, zone, config=config)
        g_magnitude = np.linalg.norm(table.g_crystal, axis=1)
        expected = -0.5 * config.wavelength_angstrom * g_magnitude**2
        assert np.allclose(table.excitation_error_inv_angstrom, expected, atol=1e-12)
        assert np.all(table.excitation_error_inv_angstrom <= 0.0)

    def test_holz_excluded_by_excitation_error(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        zone_law = table.hkl @ np.array([0, 1, 1])
        assert np.all(zone_law == 0)

    def test_intensity_sorted_and_normalized(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        assert len(table) > 0
        assert float(np.max(table.intensity)) == pytest.approx(1.0)
        assert np.all(np.diff(table.intensity) <= 1e-12)

    def test_unit_intensity_model(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=phase)
        config = KinematicSimulationConfig(intensity_model="unit")
        table = simulate_zone_axis_spots(phase, zone, config=config)
        assert np.allclose(table.intensity, 1.0)

    def test_relrod_damping_reduces_curved_ewald_spots(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        base = simulate_zone_axis_spots(
            phase, zone, config=KinematicSimulationConfig(intensity_model="unit")
        )
        damped = simulate_zone_axis_spots(
            phase,
            zone,
            config=KinematicSimulationConfig(
                intensity_model="unit", relrod_sigma_inv_angstrom=0.01
            ),
        )
        base_map = dict(zip(base.hkl_labels(), base.intensity, strict=False))
        damped_map = dict(zip(damped.hkl_labels(), damped.intensity, strict=False))
        outer = max(
            base_map,
            key=lambda label: float(np.linalg.norm(base.g_crystal[base.hkl_labels().index(label)])),
        )
        assert damped_map[outer] < base_map[outer]

    def test_finite_thickness_uses_exact_normalized_shape_factor(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        thickness = 100.0
        table = simulate_zone_axis_spots(
            phase,
            zone,
            config=KinematicSimulationConfig(
                intensity_model="unit",
                foil_thickness_angstrom=thickness,
                min_relative_intensity=0.0,
            ),
        )
        expected = np.asarray(
            FiniteThicknessShapeFactor(thickness).intensity_factor(
                table.excitation_error_inv_angstrom
            )
        )
        expected /= np.max(expected)

        assert table.intensity == pytest.approx(expected)
        assert "sinc^2(t s_g)" in table.describe()

    def test_g_max_filter(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        config = KinematicSimulationConfig(g_max_inv_angstrom=0.6)
        table = simulate_zone_axis_spots(phase, zone, config=config)
        g_magnitude = np.linalg.norm(table.g_crystal, axis=1)
        assert len(table) > 0
        assert np.all(g_magnitude <= 0.6 + 1e-12)

    def test_irrational_zone_axis_supported(self) -> None:
        phase = make_nickel_phase()
        direction = CrystalDirection(np.array([0.0, 1.0, 1.02]), phase=phase)
        table = simulate_zone_axis_spots(phase, direction)
        assert np.all(
            np.abs(table.excitation_error_inv_angstrom)
            <= table.config.max_excitation_error_inv_angstrom
        )

    def test_deterministic_output(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([1, 1, 1]), phase=phase)
        first = simulate_zone_axis_spots(phase, zone)
        second = simulate_zone_axis_spots(phase, zone)
        assert np.array_equal(first.hkl, second.hkl)
        assert np.array_equal(first.detector_mm, second.detector_mm)
        assert np.array_equal(first.intensity, second.intensity)

    def test_explicit_shared_basis_round_trip(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        basis = zone_basis_from_axis(zone.unit_vector, in_plane_rotation_deg=30.0)
        table = simulate_zone_axis_spots(phase, zone, basis=basis)
        default = simulate_zone_axis_spots(phase, zone)
        angle = np.deg2rad(30.0)
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        default_map = dict(zip(default.hkl_labels(), default.detector_mm, strict=False))
        for label, coords in zip(table.hkl_labels(), table.detector_mm, strict=False):
            assert np.allclose(coords, rotation @ default_map[label], atol=1e-9)

    def test_wrong_phase_zone_axis_raises(self) -> None:
        nickel = make_nickel_phase()
        iron = make_iron_bcc_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=iron)
        with pytest.raises(ValueError, match=r"zone_axis\.phase"):
            simulate_zone_axis_spots(nickel, zone)

    def test_invalid_basis_raises(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=phase)
        with pytest.raises(ValueError, match="orthonormal"):
            simulate_zone_axis_spots(phase, zone, basis=np.eye(3) * 2.0)
        wrong_zone = zone_basis_from_axis(np.array([1.0, 0.0, 0.0]))
        with pytest.raises(ValueError, match="zone-axis"):
            simulate_zone_axis_spots(phase, zone, basis=wrong_zone)

    def test_empty_pattern_allowed(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 0, 1]), phase=phase)
        config = KinematicSimulationConfig(g_max_inv_angstrom=0.05)
        table = simulate_zone_axis_spots(phase, zone, config=config)
        assert len(table) == 0
        assert "empty" in table.describe()

    def test_spot_table_arrays_read_only(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        with pytest.raises(ValueError):
            table.intensity[0] = 2.0
        with pytest.raises(ValueError):
            table.hkl[0, 0] = 9

    def test_describe_content(self) -> None:
        phase = make_nickel_phase()
        zone = ZoneAxis(np.array([0, 1, 1]), phase=phase)
        table = simulate_zone_axis_spots(phase, zone)
        text = table.describe()
        assert "nickel-fcc" in text
        assert "[011]" in text
        assert "200 kV" in text
        assert "0.025079" in text
        assert "s_g" in text
        assert str(len(table)) in text


class TestConfigValidation:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"beam_energy_kev": 0.0},
            {"beam_energy_kev": -5.0},
            {"camera_constant_mm_angstrom": 0.0},
            {"max_index": 0},
            {"max_index": -2},
            {"g_max_inv_angstrom": 0.0},
            {"max_excitation_error_inv_angstrom": 0.0},
            {"intensity_model": "dynamical"},
            {"relrod_sigma_inv_angstrom": -1.0},
            {"foil_thickness_angstrom": 0.0},
            {"min_relative_intensity": 1.0},
            {"min_relative_intensity": -0.1},
            {"double_diffraction_coupling": 0.0},
            {"double_diffraction_coupling": 1.5},
            {"double_diffraction_coupling": -0.2},
        ],
    )
    def test_invalid_config_raises(self, kwargs: dict[str, object]) -> None:
        with pytest.raises(ValueError):
            KinematicSimulationConfig(**kwargs)  # type: ignore[arg-type]

    def test_wavelength_property(self) -> None:
        config = KinematicSimulationConfig(beam_energy_kev=200.0)
        assert config.wavelength_angstrom == pytest.approx(0.025079, abs=5e-7)


def make_silicon_phase() -> Phase:
    """Diamond-cubic silicon, a = 5.4309 angstrom (ICDD 27-1402)."""

    crystal = _crystal_frame()
    lattice = Lattice(5.4309, 5.4309, 5.4309, 90.0, 90.0, 90.0, crystal_frame=crystal)
    face_centred = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
    coordinates = [
        tuple(np.asarray(base) + offset)
        for base in face_centred
        for offset in (np.zeros(3), np.full(3, 0.25))
    ]
    sites = tuple(
        AtomicSite(label=f"Si{i}", species="Si", fractional_coordinates=np.asarray(coords))
        for i, coords in enumerate(coordinates)
    )
    return Phase(
        "silicon",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group_symbol="Fd-3m",
    )


class TestDoubleDiffractionSums:
    def test_sums_and_weights_are_the_pairwise_products(self) -> None:
        hkl = np.array([[1, 1, 1], [1, 1, -1]], dtype=np.int64)
        intensity = np.array([0.8, 0.5])
        sums, weight, parents = double_diffraction_sums(hkl, intensity)
        lookup = {
            tuple(int(v) for v in row): (float(w), p)
            for row, w, p in zip(sums, weight, parents, strict=True)
        }
        assert lookup[(2, 2, 2)][0] == pytest.approx(0.64)
        assert lookup[(2, 2, 0)][0] == pytest.approx(0.40)
        assert lookup[(2, 2, -2)][0] == pytest.approx(0.25)
        assert np.array_equal(lookup[(2, 2, 0)][1].sum(axis=0), np.array([2, 2, 0]))

    def test_transmitted_beam_is_excluded(self) -> None:
        hkl = np.array([[1, 1, 1], [-1, -1, -1]], dtype=np.int64)
        sums, _, _ = double_diffraction_sums(hkl, np.array([1.0, 1.0]))
        assert not np.any(np.all(sums == 0, axis=1))

    def test_repeated_paths_accumulate(self) -> None:
        # (2 0 0) is reachable both as (1 1 1) + (1 -1 -1) and as (1 -1 1) +
        # (1 1 -1), so its weight is the sum of the two path products.
        hkl = np.array([[1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, -1]], dtype=np.int64)
        sums, weight, _ = double_diffraction_sums(hkl, np.full(4, 0.5))
        row = int(np.flatnonzero(np.all(sums == np.array([2, 0, 0]), axis=1))[0])
        assert float(weight[row]) == pytest.approx(0.5)

    def test_empty_input(self) -> None:
        sums, weight, parents = double_diffraction_sums(
            np.zeros((0, 3), dtype=np.int64), np.zeros(0)
        )
        assert sums.shape == (0, 3)
        assert weight.shape == (0,)
        assert parents.shape == (0, 2, 3)

    def test_parents_always_sum_to_their_reflection(self) -> None:
        rng = np.random.default_rng(7)
        hkl = np.unique(rng.integers(-3, 4, size=(40, 3)), axis=0)
        sums, _, parents = double_diffraction_sums(hkl, rng.random(hkl.shape[0]))
        assert np.array_equal(parents.sum(axis=1), sums)

    @pytest.mark.parametrize(
        "bad",
        [np.zeros((4, 2), dtype=np.int64), np.zeros(3, dtype=np.int64)],
    )
    def test_malformed_indices_raise(self, bad: np.ndarray) -> None:
        with pytest.raises(ValueError, match=r"shape \(M, 3\)"):
            double_diffraction_sums(bad, np.zeros(bad.shape[0]))

    def test_mismatched_intensity_length_raises(self) -> None:
        with pytest.raises(ValueError, match=r"shape \(M,\)"):
            double_diffraction_sums(np.zeros((4, 3), dtype=np.int64), np.zeros(3))


class TestDoubleDiffractionInPatterns:
    """Silicon [110] is the textbook case.

    In diamond cubic, ``F`` vanishes unless ``h, k, l`` are all odd or all even
    with ``h + k + l = 4n``, so 002 is forbidden. It is nevertheless present in
    every recorded Si [110] pattern, produced by the two-step path
    ``(1 1 1) + (-1 -1 1)``; see Williams and Carter, *Transmission Electron
    Microscopy*, 2nd ed. (2009), ch. 16.
    """

    @staticmethod
    def _silicon_table(*, include: bool) -> SpotTable:
        phase = make_silicon_phase()
        config = KinematicSimulationConfig(
            max_index=4, g_max_inv_angstrom=1.2, include_double_diffraction=include
        )
        return simulate_zone_axis_spots(phase, ZoneAxis([1, 1, 0], phase=phase), config=config)

    def test_forbidden_002_absent_by_default(self) -> None:
        table = self._silicon_table(include=False)
        assert not np.any(np.all(np.abs(table.hkl) == np.array([0, 0, 2]), axis=1))
        assert not np.any(table.is_double_diffraction)

    def test_forbidden_002_appears_and_is_designated(self) -> None:
        table = self._silicon_table(include=True)
        rows = np.flatnonzero(np.all(np.abs(table.hkl) == np.array([0, 0, 2]), axis=1))
        assert rows.size == 2
        for row in rows:
            assert bool(table.is_double_diffraction[row])
            assert float(table.structure_factor_amplitude[row]) == pytest.approx(0.0, abs=1e-9)
            parents = table.double_diffraction_parents[row]
            assert np.array_equal(parents.sum(axis=0), table.hkl[row])
            assert np.all(np.abs(parents) == 1)

    def test_added_spots_are_weaker_than_the_genuine_ones(self) -> None:
        table = self._silicon_table(include=True)
        forbidden = table.is_double_diffraction
        assert float(np.max(table.intensity[forbidden])) < float(
            np.min(table.intensity[~forbidden])
        )
        assert float(np.max(table.intensity)) == pytest.approx(1.0)

    def test_enabling_only_adds_rows(self) -> None:
        plain = self._silicon_table(include=False)
        extended = self._silicon_table(include=True)
        plain_rows = {tuple(int(v) for v in row) for row in plain.hkl}
        extended_rows = {tuple(int(v) for v in row) for row in extended.hkl}
        assert plain_rows < extended_rows
        marked = {
            tuple(int(v) for v in row) for row in extended.hkl[extended.is_double_diffraction]
        }
        assert extended_rows - plain_rows == marked

    def test_coupling_scales_the_added_intensity_linearly(self) -> None:
        phase = make_silicon_phase()
        zone = ZoneAxis([1, 1, 0], phase=phase)
        intensities = []
        for coupling in (0.02, 0.04):
            table = simulate_zone_axis_spots(
                phase,
                zone,
                config=KinematicSimulationConfig(
                    max_index=4,
                    g_max_inv_angstrom=1.2,
                    include_double_diffraction=True,
                    double_diffraction_coupling=coupling,
                ),
            )
            row = int(np.flatnonzero(np.all(table.hkl == np.array([0, 0, 2]), axis=1))[0])
            intensities.append(float(table.intensity[row]))
        assert intensities[1] == pytest.approx(2.0 * intensities[0])

    def test_describe_reports_the_designation(self) -> None:
        text = self._silicon_table(include=True).describe()
        assert "kinematically forbidden" in text
        assert "double diffraction" in text
        assert "Double diffraction is disabled" in self._silicon_table(include=False).describe()

    @pytest.mark.parametrize(
        ("factory", "rule"),
        [(make_iron_bcc_phase, "I"), (make_nickel_phase, "F")],
    )
    def test_centring_absences_are_never_revived(
        self, factory: Callable[[], Phase], rule: str
    ) -> None:
        # Centring conditions define a sublattice of reciprocal space, and a
        # sublattice is closed under addition, so no sum of two allowed
        # reflections can land on a centring absence. Double diffraction can
        # only revive a basis absence.
        phase = factory()
        table = simulate_zone_axis_spots(
            phase,
            ZoneAxis([1, 1, 0], phase=phase),
            config=KinematicSimulationConfig(
                max_index=5, include_double_diffraction=True, double_diffraction_coupling=1.0
            ),
        )
        allowed = centering_allowed_mask(table.hkl, ReflectionCondition(centering=rule))
        assert bool(np.all(allowed))

    def test_origin_label_is_empty_for_genuine_reflections(self) -> None:
        table = self._silicon_table(include=True)
        assert table.double_diffraction_origin_label(0) == ""
        row = int(np.flatnonzero(table.is_double_diffraction)[0])
        assert "=" in table.double_diffraction_origin_label(row)

    def test_designation_arrays_are_read_only(self) -> None:
        table = self._silicon_table(include=True)
        with pytest.raises(ValueError):
            table.is_double_diffraction[0] = True
        with pytest.raises(ValueError):
            table.double_diffraction_parents[0, 0, 0] = 9


class TestSpotTableDesignationValidation:
    @staticmethod
    def _minimal_kwargs(phase: Phase) -> dict[str, object]:
        table = simulate_zone_axis_spots(
            phase,
            ZoneAxis([0, 0, 1], phase=phase),
            config=KinematicSimulationConfig(max_index=2),
        )
        return {
            "phase": table.phase,
            "zone_axis": table.zone_axis,
            "basis": table.basis,
            "config": table.config,
            "hkl": table.hkl,
            "g_crystal": table.g_crystal,
            "g_detector_inv_angstrom": table.g_detector_inv_angstrom,
            "detector_mm": table.detector_mm,
            "d_spacing_angstrom": table.d_spacing_angstrom,
            "structure_factor_amplitude": table.structure_factor_amplitude,
            "intensity": table.intensity,
            "excitation_error_inv_angstrom": table.excitation_error_inv_angstrom,
        }

    def test_wrong_flag_length_raises(self) -> None:
        kwargs = self._minimal_kwargs(make_nickel_phase())
        with pytest.raises(ValueError, match="is_double_diffraction"):
            SpotTable(**kwargs, is_double_diffraction=np.zeros(3, dtype=bool))  # type: ignore[arg-type]

    def test_wrong_parent_shape_raises(self) -> None:
        kwargs = self._minimal_kwargs(make_nickel_phase())
        with pytest.raises(ValueError, match="double_diffraction_parents"):
            SpotTable(**kwargs, double_diffraction_parents=np.zeros((2, 2, 3), dtype=np.int64))  # type: ignore[arg-type]

    def test_parents_that_do_not_sum_to_their_reflection_raise(self) -> None:
        kwargs = self._minimal_kwargs(make_nickel_phase())
        count = int(np.asarray(kwargs["hkl"]).shape[0])  # type: ignore[arg-type]
        flags = np.zeros(count, dtype=bool)
        flags[0] = True
        with pytest.raises(ValueError, match="must sum to the reflection"):
            SpotTable(  # type: ignore[arg-type]
                **kwargs,
                is_double_diffraction=flags,
                double_diffraction_parents=np.zeros((count, 2, 3), dtype=np.int64),
            )
