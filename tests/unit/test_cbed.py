"""Tests for `pytex.diffraction.cbed`.

The absolute scale is the thing that can silently be wrong here: a structure
factor with a missing relativistic factor, or an extinction distance off by
``8 pi^2 a0``, still produces plausible-looking fringes. So the extinction
distances are checked against **published values for aluminium** (Williams and
Carter, *Transmission Electron Microscopy*, 2nd ed., Table 23.1), and the
thickness analysis is checked by a round trip: simulate fringes at a known
thickness, read the minima back, and recover it.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import crystal_frame
from pytex.core.lattice import AtomicSite, Lattice, Phase, SpaceGroupSpec, UnitCell, ZoneAxis
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.cbed import (
    CBED_PATTERN_SCHEMA,
    CBED_THICKNESS_SCHEMA,
    ConvergentBeamConfig,
    electron_structure_factor_angstrom,
    extinction_distance_angstrom,
    fringe_minimum_excitation_errors,
    holz_ring_radii_inv_angstrom,
    simulate_cbed_pattern,
    thickness_from_fringe_minima,
    two_beam_rocking_curve,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom
from pytex.diffraction.scattering import electron_scattering_factors, xray_form_factors

#: Aluminium extinction distances at 100 kV, in angstrom, from Williams and
#: Carter, *Transmission Electron Microscopy* (2nd ed.), Table 23.1. Aluminium is
#: the calibration case because the fitted scattering-factor parametrization is
#: most accurate for light elements.
ALUMINIUM_EXTINCTION_DISTANCES_100KV = {(1, 1, 1): 556.0, (2, 0, 0): 673.0, (2, 2, 0): 1057.0}


def _cubic_fcc_phase(name: str, species: str, parameter_angstrom: float) -> Phase:
    frame = crystal_frame()
    lattice = Lattice(
        parameter_angstrom,
        parameter_angstrom,
        parameter_angstrom,
        90.0,
        90.0,
        90.0,
        crystal_frame=frame,
    )
    sites = tuple(
        AtomicSite(
            label=f"{species}{index}",
            species=species,
            fractional_coordinates=np.asarray(position, dtype=np.float64),
        )
        for index, position in enumerate(
            [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
        )
    )
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=frame),
    )


@pytest.fixture(scope="module")
def aluminium() -> Phase:
    """FCC aluminium, built inline so the benchmark needs no optional dependency."""

    return _cubic_fcc_phase("aluminium-fcc", "Al", 4.0495)


@pytest.fixture(scope="module")
def nickel() -> Phase:
    """FCC nickel with the fixture lattice parameter, for the FCC worked case."""

    return _cubic_fcc_phase("nickel-fcc", "Ni", 3.5239)


# --------------------------------------------------------------------------- #
# Electron scattering factors and the absolute scale
# --------------------------------------------------------------------------- #


def test_mott_bethe_inverts_the_tabulated_x_ray_parametrization() -> None:
    """``f_e = (Z - f_x) / (8 pi^2 a0 s^2)`` must reproduce the fitted coefficients.

    The X-ray table is stored as ``Z - 41.78214 s^2 sum a_i exp(-b_i s^2)``, so
    the Mott-Bethe inversion should return exactly ``sum a_i exp(-b_i s^2)``.
    Checking that identity pins the constant: a wrong prefactor would scale
    every extinction distance without changing the shape of anything.
    """

    s_values = np.array([0.05, 0.2, 0.5, 1.0])
    for species in ("Al", "Ni", "Zr"):
        from pytex.core._chemistry import atomic_number

        expected = (float(atomic_number(species)) - xray_form_factors(species, s_values)) / (
            41.78214 * s_values * s_values
        )
        assert electron_scattering_factors(species, s_values) == pytest.approx(
            expected, rel=1e-12
        )


def test_electron_scattering_factor_is_finite_at_zero_angle() -> None:
    """The Mott-Bethe singularity at ``s = 0`` is removable, not a division by zero."""

    value = electron_scattering_factors("Al", np.array([0.0]))
    assert np.all(np.isfinite(value))
    assert float(value[0]) > 0.0


@pytest.mark.parametrize("hkl,reference", sorted(ALUMINIUM_EXTINCTION_DISTANCES_100KV.items()))
def test_extinction_distances_match_published_aluminium_values(
    aluminium: Phase, hkl: tuple[int, int, int], reference: float
) -> None:
    """The absolute scale, against a citation rather than a prior program output."""

    computed = float(extinction_distance_angstrom(aluminium, hkl, beam_energy_kev=100.0)[0])
    assert computed == pytest.approx(reference, rel=0.02)


def test_a_forbidden_reflection_has_an_infinite_extinction_distance(aluminium: Phase) -> None:
    """``F_g = 0`` means no coupling and no oscillation; ``inf`` is the answer."""

    assert np.isinf(extinction_distance_angstrom(aluminium, (1, 0, 0))[0])
    assert abs(electron_structure_factor_angstrom(aluminium, (1, 0, 0))[0]) < 1e-9


def test_the_relativistic_correction_is_present(aluminium: Phase) -> None:
    """Omitting ``gamma`` would make every extinction distance 39 percent too long.

    Checked structurally: the ratio of structure factors at two voltages must be
    the ratio of ``1 + E / m0 c^2``, since nothing else in ``F_g`` depends on the
    accelerating voltage.
    """

    low = abs(electron_structure_factor_angstrom(aluminium, (1, 1, 1), beam_energy_kev=100.0)[0])
    high = abs(electron_structure_factor_angstrom(aluminium, (1, 1, 1), beam_energy_kev=300.0)[0])
    expected = (1.0 + 300.0 / 510.99895) / (1.0 + 100.0 / 510.99895)
    assert high / low == pytest.approx(expected, rel=1e-12)


def test_structure_factors_require_a_unit_cell() -> None:
    frame = crystal_frame()
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=frame)
    bare = Phase(
        "bare",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
    )
    with pytest.raises(ValueError, match="carries no unit cell"):
        electron_structure_factor_angstrom(bare, (1, 1, 1))


# --------------------------------------------------------------------------- #
# The two-beam rocking curve
# --------------------------------------------------------------------------- #


def test_rocking_curve_at_exact_bragg_is_the_pendelloesung_sine() -> None:
    """``I_g(0) = sin^2(pi t / xi)``: complete oscillation with thickness."""

    xi = 500.0
    thicknesses = np.array([100.0, 250.0, 500.0, 750.0])
    values = np.array(
        [
            float(
                two_beam_rocking_curve(
                    0.0, thickness_angstrom=t, extinction_distance_angstrom=xi
                )
            )
            for t in thicknesses
        ]
    )
    assert values == pytest.approx(np.square(np.sin(np.pi * thicknesses / xi)), abs=1e-12)
    # At t = xi / 2 the beam is fully diffracted; at t = xi it is fully back.
    assert values[1] == pytest.approx(1.0, abs=1e-12)
    assert values[2] == pytest.approx(0.0, abs=1e-12)


def test_rocking_curve_minima_sit_where_the_theory_says(aluminium: Phase) -> None:
    """Zeros at ``t s_eff = n``: the relation the thickness method inverts."""

    xi = float(extinction_distance_angstrom(aluminium, (1, 1, 1))[0])
    thickness = 1500.0
    s_values = np.linspace(1e-6, 0.02, 200_001)
    intensity = two_beam_rocking_curve(
        s_values, thickness_angstrom=thickness, extinction_distance_angstrom=xi
    )
    minima = fringe_minimum_excitation_errors(s_values, intensity)
    orders = np.round(thickness * np.sqrt(minima**2 + xi**-2))
    assert np.all(orders >= 1.0)
    predicted = np.sqrt(np.maximum((orders / thickness) ** 2 - xi**-2, 0.0))
    assert minima == pytest.approx(predicted, abs=2e-6)


def test_rocking_curve_is_bounded_and_never_negative() -> None:
    s_values = np.linspace(-0.05, 0.05, 5001)
    intensity = two_beam_rocking_curve(
        s_values, thickness_angstrom=800.0, extinction_distance_angstrom=350.0
    )
    assert float(intensity.min()) >= 0.0
    assert float(intensity.max()) <= 1.0 + 1e-12


def test_a_forbidden_reflection_diffracts_nothing() -> None:
    values = two_beam_rocking_curve(
        np.linspace(-0.01, 0.01, 11),
        thickness_angstrom=1000.0,
        extinction_distance_angstrom=float("inf"),
    )
    assert values == pytest.approx(np.zeros(11))


def test_rocking_curve_rejects_a_non_physical_thickness() -> None:
    with pytest.raises(ValueError, match="thickness_angstrom"):
        two_beam_rocking_curve(0.0, thickness_angstrom=0.0, extinction_distance_angstrom=300.0)


# --------------------------------------------------------------------------- #
# Thickness determination: the round trip
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("thickness", [600.0, 1200.0, 2500.0])
def test_thickness_round_trips_through_the_fringe_analysis(
    aluminium: Phase, thickness: float
) -> None:
    """Simulate fringes at a known thickness; read them back; recover it.

    This is the end-to-end proof of the method, and it recovers the extinction
    distance too — which is what makes CBED thickness determination independent
    of a tabulated constant.
    """

    xi = float(extinction_distance_angstrom(aluminium, (1, 1, 1))[0])
    s_values = np.linspace(1e-6, 0.02, 200_001)
    intensity = two_beam_rocking_curve(
        s_values, thickness_angstrom=thickness, extinction_distance_angstrom=xi
    )
    minima = fringe_minimum_excitation_errors(s_values, intensity)
    report = thickness_from_fringe_minima(minima)

    assert report.thickness_angstrom == pytest.approx(thickness, rel=1e-3)
    assert report.extinction_distance_angstrom == pytest.approx(xi, rel=1e-2)
    assert report.r_squared > 0.999
    assert report.to_json_dict()["schema"] == CBED_THICKNESS_SCHEMA
    assert "does not depend on a tabulated extinction distance" in report.describe()


def test_the_order_search_recovers_the_true_first_order(aluminium: Phase) -> None:
    """The innermost visible minimum is usually not ``n = 1``, and assuming so is wrong."""

    xi = float(extinction_distance_angstrom(aluminium, (1, 1, 1))[0])
    thickness = 1800.0
    s_values = np.linspace(1e-6, 0.02, 200_001)
    intensity = two_beam_rocking_curve(
        s_values, thickness_angstrom=thickness, extinction_distance_angstrom=xi
    )
    minima = fringe_minimum_excitation_errors(s_values, intensity)
    true_first_order = round(thickness * float(np.sqrt(minima[0] ** 2 + xi**-2)))

    searched = thickness_from_fringe_minima(minima)
    assert searched.first_order == true_first_order
    assert true_first_order > 1, "the case is only interesting when n = 1 is wrong"

    # Forcing the classic wrong assumption does not merely bias the answer here:
    # it tilts the fitted line the wrong way, so the implied 1/xi^2 is negative
    # and there is no thickness to report. Failing loudly beats a plausible number.
    with pytest.raises(ValueError, match="No order assignment"):
        thickness_from_fringe_minima(minima, first_order=1)


def test_thickness_analysis_needs_two_minima() -> None:
    with pytest.raises(ValueError, match="At least two fringe minima"):
        thickness_from_fringe_minima([0.004])


def test_thickness_analysis_refuses_data_that_are_not_two_beam_fringes() -> None:
    """A plausible number from unphysical data would be worse than an error.

    The two-beam law ``s_n = sqrt((n/t)^2 - 1/xi^2)`` grows towards ``n/t``, so
    ``(s_n/n)^2`` must *increase* with ``n``. Minima that crowd together instead
    of spreading out violate that for every order assignment, and there is no
    thickness to report.
    """

    crowding = [0.00100, 0.00110, 0.00115, 0.00118, 0.00120]
    with pytest.raises(ValueError, match="No order assignment"):
        thickness_from_fringe_minima(crowding, max_first_order=6)


# --------------------------------------------------------------------------- #
# HOLZ geometry
# --------------------------------------------------------------------------- #


def test_holz_radius_follows_the_layer_spacing(nickel: Phase) -> None:
    """``G_n = sqrt(2 n H / lambda)`` with ``H = 1 / |r_uvw|``, for FCC [001]."""

    zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
    orders, radii = holz_ring_radii_inv_angstrom(nickel, zone, beam_energy_kev=200.0)
    wavelength = electron_wavelength_angstrom(200.0)
    layer_spacing = 1.0 / 3.5239
    expected = np.sqrt(2.0 * orders * layer_spacing / wavelength)
    assert radii == pytest.approx(expected, rel=1e-12)
    assert list(orders) == [1, 2]


def test_holz_rings_are_reported_only_for_layers_that_can_diffract(nickel: Phase) -> None:
    """A centred lattice can extinguish a whole Laue zone; reporting it would mislead.

    For FCC down ``[111]`` the layer index of an allowed reflection is
    ``h + k + l``, which for all-odd or all-even indices is never an odd number
    that is not a multiple of... — rather than assert the arithmetic in prose,
    the test simply requires that every reported order does have an allowed
    reflection, and that a layer with none is absent from the result.
    """

    zone = ZoneAxis(np.array([1, 1, 1]), phase=nickel)
    orders, radii = holz_ring_radii_inv_angstrom(nickel, zone, beam_energy_kev=200.0, orders=3)
    assert orders.size == radii.size

    values = np.arange(-6, 7)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid != 0, axis=1)]
    all_odd = np.all(grid % 2 != 0, axis=1)
    all_even = np.all(grid % 2 == 0, axis=1)
    allowed = grid[all_odd | all_even]
    layers = set(int(value) for value in np.abs(allowed @ np.array([1, 1, 1])))
    for order in orders:
        assert int(order) in layers


def test_holz_rejects_a_zone_axis_of_another_phase(nickel: Phase, aluminium: Phase) -> None:
    with pytest.raises(ValueError, match=r"zone_axis\.phase must match phase"):
        holz_ring_radii_inv_angstrom(nickel, ZoneAxis(np.array([0, 0, 1]), phase=aluminium))


# --------------------------------------------------------------------------- #
# Pattern simulation
# --------------------------------------------------------------------------- #


def test_disc_radius_is_the_convergence_angle_times_the_camera_length(nickel: Phase) -> None:
    """``R = (L lambda) alpha / lambda = L alpha``: the wavelength cancels."""

    config = ConvergentBeamConfig(
        beam_energy_kev=200.0, convergence_semi_angle_mrad=5.0,
        camera_constant_mm_angstrom=180.0,
    )
    wavelength = electron_wavelength_angstrom(200.0)
    assert config.disc_radius_inv_angstrom == pytest.approx(5.0e-3 / wavelength)
    assert config.disc_radius_mm == pytest.approx(180.0 * 5.0e-3 / wavelength)

    pattern = simulate_cbed_pattern(
        nickel, ZoneAxis(np.array([0, 0, 1]), phase=nickel), config=config
    )
    for disc in pattern.discs:
        assert disc.radius_mm == pytest.approx(config.disc_radius_mm)


def test_convergence_angle_selects_the_regime(nickel: Phase) -> None:
    """Small alpha: separated discs. Large alpha: overlap. The threshold is geometric."""

    zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
    tight = simulate_cbed_pattern(
        nickel, zone, config=ConvergentBeamConfig(convergence_semi_angle_mrad=2.0)
    )
    assert tight.is_kossel_moellenstedt
    assert tight.regime == "kossel-moellenstedt"

    # Push the convergence past half the closest disc separation.
    threshold_mrad = (
        tight.nearest_disc_separation_mm
        / (2.0 * tight.config.camera_constant_mm_angstrom)
        * tight.config.wavelength_angstrom
        * 1e3
    )
    wide = simulate_cbed_pattern(
        nickel,
        zone,
        config=ConvergentBeamConfig(convergence_semi_angle_mrad=threshold_mrad * 1.2),
    )
    assert not wide.is_kossel_moellenstedt
    assert wide.regime == "kossel"
    assert "overlap" in wide.describe()


def test_excitation_error_varies_linearly_along_g_and_not_across_it(nickel: Phase) -> None:
    """The reason CBED fringes are straight lines perpendicular to ``g``."""

    zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
    pattern = simulate_cbed_pattern(
        nickel, zone, config=ConvergentBeamConfig(disc_samples=41)
    )
    disc = next(d for d in pattern.discs if not d.is_transmitted)
    direction = disc.g_detector_inv_angstrom
    s_map = disc.excitation_error_inv_angstrom
    middle = s_map.shape[0] // 2

    # Along the axis carrying the larger component of g, s must change; along a
    # line of constant projection onto g it must not.
    if abs(direction[0]) >= abs(direction[1]):
        varying, constant = s_map[:, middle], s_map[middle, :]
    else:
        varying, constant = s_map[middle, :], s_map[:, middle]
    varying = varying[np.isfinite(varying)]
    constant = constant[np.isfinite(constant)]
    assert float(np.ptp(varying)) > 1e-6
    assert float(np.ptp(constant)) < 1e-12
    # Linear: second differences vanish.
    assert np.allclose(np.diff(varying, n=2), 0.0, atol=1e-12)


def test_transmitted_disc_is_the_two_beam_complement(nickel: Phase) -> None:
    """Two beams exchange intensity; nothing is absorbed in this model."""

    pattern = simulate_cbed_pattern(
        nickel, ZoneAxis(np.array([0, 0, 1]), phase=nickel),
        config=ConvergentBeamConfig(disc_samples=31),
    )
    transmitted = pattern.transmitted_disc
    strongest = pattern.discs[1]
    inside = np.isfinite(transmitted.intensity)
    assert transmitted.is_transmitted
    assert (transmitted.intensity[inside] + strongest.intensity[inside]) == pytest.approx(
        np.ones(int(inside.sum())), abs=1e-12
    )


def test_simulated_disc_fringes_recover_the_input_thickness(nickel: Phase) -> None:
    """The full round trip: pattern in, thickness out.

    This is the capability claim of the whole module, so it is tested through
    the public simulation path rather than through the rocking-curve function
    alone.
    """

    thickness = 1400.0
    pattern = simulate_cbed_pattern(
        nickel,
        ZoneAxis(np.array([0, 0, 1]), phase=nickel),
        config=ConvergentBeamConfig(
            thickness_angstrom=thickness,
            convergence_semi_angle_mrad=12.0,
            disc_samples=1501,
        ),
    )
    disc = pattern.disc_for((2, 0, 0))
    s_values, intensity = disc.radial_profile()
    minima = fringe_minimum_excitation_errors(s_values, intensity)

    # The disc is centred at s = -lambda g^2 / 2, so the exact Bragg condition
    # sits off-centre and the two branches of the rocking curve are unequal.
    # Read the richer branch, as an experimenter does.
    negative, positive = minima[minima < 0.0], minima[minima > 0.0]
    branch = negative if negative.size >= positive.size else positive
    assert branch.size >= 3

    report = thickness_from_fringe_minima(branch)
    assert report.thickness_angstrom == pytest.approx(thickness, rel=0.05)
    assert report.extinction_distance_angstrom == pytest.approx(
        disc.extinction_distance_angstrom, rel=0.1
    )


def test_pattern_describe_and_json_stay_in_lockstep(nickel: Phase) -> None:
    pattern = simulate_cbed_pattern(
        nickel, ZoneAxis(np.array([0, 0, 1]), phase=nickel),
        config=ConvergentBeamConfig(disc_samples=21),
    )
    payload = pattern.to_json_dict()
    prose = pattern.describe()

    assert payload["schema"] == CBED_PATTERN_SCHEMA
    assert payload["regime"] == pattern.regime
    assert payload["zone_axis"] == [0, 0, 1]
    assert len(payload["discs"]) == len(pattern.discs)
    assert payload["disc_radius_mm"] == pytest.approx(pattern.config.disc_radius_mm)
    # The limits must be stated, not implied.
    assert "two-beam" in prose
    assert "diffraction-group" in prose


def test_pattern_reports_a_missing_reflection_rather_than_a_wrong_one(nickel: Phase) -> None:
    pattern = simulate_cbed_pattern(
        nickel, ZoneAxis(np.array([0, 0, 1]), phase=nickel),
        config=ConvergentBeamConfig(disc_samples=11),
    )
    with pytest.raises(KeyError, match="not in this pattern"):
        pattern.disc_for((1, 0, 0))


def test_a_zero_convergence_angle_is_rejected_as_saed() -> None:
    with pytest.raises(ValueError, match="parallel beam"):
        ConvergentBeamConfig(convergence_semi_angle_mrad=0.0)


def test_simulation_requires_a_unit_cell() -> None:
    frame = crystal_frame()
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=frame)
    bare = Phase(
        "bare",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
    )
    with pytest.raises(ValueError, match="carries no unit cell"):
        simulate_cbed_pattern(bare, ZoneAxis(np.array([0, 0, 1]), phase=bare))
