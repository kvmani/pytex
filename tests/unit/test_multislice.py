"""Multislice engine: potentials, propagation, imaging and series against independent answers.

Every expected value here comes from an analytic identity, a closed form, or a
different method (Bloch waves) - never from a previous run of the engine.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest
import scipy.fft
import scipy.special

from pytex.app.phases import phase_from_request
from pytex.core.lattice import ZoneAxis
from pytex.diffraction.dynamical import beam_set_for_zone, solve_bloch_waves
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    MicroscopeAberrations,
    relativistic_interaction_parameter_inv_v_angstrom,
)
from pytex.diffraction.multislice import (
    POTENTIAL_PREFACTOR_V_ANGSTROM2,
    MultisliceGrid,
    PotentialParametrization,
    SlicedPotential,
    TemporalCoherence,
    antialias_aperture,
    fresnel_propagator,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    simulate_multislice_hrem,
    slice_potential,
    zone_axis_cell,
)


def _phase(name: str):  # type: ignore[no-untyped-def]
    return phase_from_request({"builtin": name})[1]


def _lens(**overrides: float) -> MicroscopeAberrations:
    base = MicroscopeAberrations(
        energy_kev=200.0,
        defocus_angstrom=0.0,
        cs_mm=0.0,
        focal_spread_angstrom=0.0,
        convergence_semiangle_mrad=0.0,
    )
    return replace(base, **overrides)  # type: ignore[arg-type]


# --- Potentials -------------------------------------------------------------


def test_prefactor_is_two_pi_bohr_radius_times_electron_charge() -> None:
    # Kirkland (2010) eq. 5.9: h^2 / (2 pi m0 e) = 2 pi a0 e with a0 = 0.529177 Å
    # and e^2/(4 pi eps0) = 14.3996 eV Å.
    assert POTENTIAL_PREFACTOR_V_ANGSTROM2 == pytest.approx(
        2.0 * math.pi * 0.529177210903 * 14.39964548, rel=1e-8
    )


def test_three_parametrizations_agree_for_silicon() -> None:
    g = np.linspace(0.0, 2.0, 9)
    lobato = parametrized_electron_scattering_factor("Si", g, "lobato")
    kirkland = parametrized_electron_scattering_factor("Si", g, "kirkland")
    mott = parametrized_electron_scattering_factor("Si", g, "mott_bethe")
    np.testing.assert_allclose(kirkland, lobato, rtol=0.01)
    np.testing.assert_allclose(mott, lobato, rtol=0.03)
    # Lobato's zero-angle limit is 2 sum(a_i), from the closed form.
    assert float(lobato[0]) == pytest.approx(5.836, abs=1e-3)


def test_unknown_element_is_refused() -> None:
    with pytest.raises(ValueError, match="No lobato"):
        parametrized_electron_scattering_factor("Xx", 0.0)


def test_projected_potential_integrates_to_the_zero_angle_scattering_factor() -> None:
    grid = MultisliceGrid((60, 50), (7.0, 8.0))
    species = ["Si", "O", "O"]
    xy = np.array([[1.0, 2.0], [3.3, 4.1], [6.9, 7.95]])
    v = slice_potential(species, xy, grid)
    dx, dy = grid.sampling_angstrom
    expected = POTENTIAL_PREFACTOR_V_ANGSTROM2 * sum(
        float(parametrized_electron_scattering_factor(s, 0.0)) for s in species
    )
    assert float(np.sum(v) * dx * dy) == pytest.approx(expected, rel=1e-12)


def test_kirkland_potential_matches_its_real_space_closed_form() -> None:
    # Kirkland (2010) eq. C.19: the infinite projection of the Kirkland fit is
    # v(r) = 4 pi^2 a0 e sum a_i K0(2 pi r sqrt(b_i)) + 2 pi^2 a0 e sum c_i/d_i
    # exp(-pi^2 r^2 / d_i). Away from the core the grid synthesis must agree.
    from pytex.diffraction.multislice import _parametrization_table

    a, b, c, d = _parametrization_table()["kirkland"]["Au"]
    a0e = POTENTIAL_PREFACTOR_V_ANGSTROM2 / (2.0 * math.pi)
    n, length = 1200, 12.0
    grid = MultisliceGrid((n, n), (length, length))
    v = slice_potential(["Au"], np.array([[6.0, 6.0]]), grid, "kirkland")
    columns = np.arange(n)[int(6.5 / 0.01) : int(8.0 / 0.01)]  # 0.5 Å to 2 Å from the core
    r = columns * 0.01 - 6.0
    closed = 4 * math.pi**2 * a0e * np.sum(
        a[:, None] * scipy.special.k0(2 * math.pi * r[None] * np.sqrt(b[:, None])), axis=0
    ) + 2 * math.pi**2 * a0e * np.sum(
        (c / d)[:, None] * np.exp(-(math.pi**2) * r[None] ** 2 / d[:, None]), axis=0
    )
    # The grid synthesis truncates f_e at the Nyquist frequency, which rings about
    # the closed form (Gibbs); at 0.01 Å sampling the ringing is under 1 % of the
    # potential 0.5 Å from the nucleus.
    np.testing.assert_allclose(v[n // 2, columns], closed, atol=0.01 * float(closed[0]))


def test_static_debye_waller_damping_preserves_the_mean_and_lowers_the_peak() -> None:
    grid = MultisliceGrid((64, 64), (6.0, 6.0))
    xy = np.array([[3.0, 3.0]])
    sharp = slice_potential(["Cu"], xy, grid)
    blurred = slice_potential(["Cu"], xy, grid, debye_waller_sigma_angstrom=0.1)
    assert float(np.mean(blurred)) == pytest.approx(float(np.mean(sharp)), rel=1e-12)
    assert float(blurred.max()) < float(sharp.max())


# --- Periodic zone-axis cells -----------------------------------------------


@pytest.mark.parametrize(
    ("name", "zone", "vectors", "cells"),
    [
        ("si_diamond", (1, 1, 0), [[0, 0, 1], [1, -1, 0], [1, 1, 0]], 2),
        ("si_diamond", (2, 2, 0), [[0, 0, 1], [1, -1, 0], [1, 1, 0]], 2),
        ("al_fcc", (0, 0, 1), [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1),
        ("ti_hcp", (0, 0, 1), [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], 2),
    ],
)
def test_zone_axis_cell_uses_perpendicular_lattice_vectors(
    name: str, zone: tuple[int, int, int], vectors: list[list[int]], cells: int
) -> None:
    cell = zone_axis_cell(_phase(name), zone)
    np.testing.assert_array_equal(cell.vectors_uvw, vectors)
    gram = cell.vectors_cartesian @ cell.vectors_cartesian.T
    np.testing.assert_allclose(gram - np.diag(np.diag(gram)), 0.0, atol=1e-9)
    assert cell.cells_per_box == cells
    assert np.linalg.det(cell.vectors_cartesian) > 0.0


def test_periodic_slab_is_exactly_periodic_and_complete() -> None:
    phase = _phase("si_diamond")
    snap, cell, repeats = periodic_slab(phase, (1, 1, 0), (2, 1), thickness_angstrom=15.0)
    assert repeats == (2, 1, 2)
    assert snap.natoms == 8 * cell.cells_per_box * 4
    lx, ly, _ = np.diag(snap.cell)
    assert lx == pytest.approx(2 * 5.43102) and ly == pytest.approx(5.43102 * math.sqrt(2))
    # Nearest-neighbour distance across the periodic boundary is the bond length.
    delta = snap.positions[:, None, :] - snap.positions[None, :, :]
    delta[..., 0] -= lx * np.round(delta[..., 0] / lx)
    delta[..., 1] -= ly * np.round(delta[..., 1] / ly)
    distance = np.linalg.norm(delta, axis=-1) + np.eye(snap.natoms) * 99.0
    np.testing.assert_allclose(distance.min(axis=1), 5.43102 * math.sqrt(3) / 4, rtol=1e-9)


def test_zero_zone_axis_is_refused() -> None:
    with pytest.raises(ValueError, match="not all zero"):
        zone_axis_cell(_phase("al_fcc"), (0, 0, 0))


# --- Propagation --------------------------------------------------------------


def test_propagators_compose() -> None:
    grid = MultisliceGrid((32, 48), (5.0, 7.0))
    inside = antialias_aperture(grid) == 1.0
    p1 = fresnel_propagator(grid, 300.0, 3.0, (1.0, -2.0))
    p2 = fresnel_propagator(grid, 300.0, 4.5, (1.0, -2.0))
    p12 = fresnel_propagator(grid, 300.0, 7.5, (1.0, -2.0))
    np.testing.assert_allclose((p1 * p2)[inside], p12[inside], atol=1e-12)


def test_vacuum_leaves_the_plane_wave_unchanged() -> None:
    empty = AtomicSnapshot(species=(), positions=np.zeros((0, 3)), cell=np.diag([6.0, 6.0, 20.0]))
    wave = multislice(empty, 200.0, sampling_angstrom=0.2, slice_thickness_angstrom=2.0)
    np.testing.assert_allclose(wave.wave(), 1.0, atol=1e-12)
    assert wave.retained_intensity[-1] == pytest.approx(1.0, abs=1e-12)


def test_one_thin_slice_is_the_band_limited_phase_object() -> None:
    phase = _phase("al_fcc")
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (2, 2))
    exit_wave = multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=100.0)
    grid = exit_wave.grid
    sliced = SlicedPotential.from_snapshot(snap, grid, 100.0)
    sigma = relativistic_interaction_parameter_inv_v_angstrom(200.0)
    t = np.exp(1j * sigma * sliced.total_projected_potential())
    t = scipy.fft.ifft2(scipy.fft.fft2(t) * antialias_aperture(grid))
    expected = scipy.fft.ifft2(
        scipy.fft.fft2(t) * fresnel_propagator(grid, 200.0, exit_wave.slice_thickness_angstrom)
    )
    np.testing.assert_allclose(exit_wave.wave(), expected, atol=1e-10)


def test_exit_depths_in_one_pass_equal_separate_thinner_runs() -> None:
    phase = _phase("si_diamond")
    thick, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=6)
    thin, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=2)
    a = 5.43102
    series = multislice(
        thick, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=a / 4,
        exit_depths_angstrom=[2 * a, 6 * a],
    )
    single = multislice(thin, 200.0, grid=series.grid, slice_thickness_angstrom=a / 4)
    np.testing.assert_allclose(series.depths_angstrom, [2 * a, 6 * a])
    np.testing.assert_allclose(series.wave(0), single.wave(), atol=1e-10)


def test_repeated_crystal_slices_are_built_once() -> None:
    snap, _, _ = periodic_slab(_phase("si_diamond"), (0, 0, 1), (1, 1), beam_repeats=10)
    grid = MultisliceGrid.from_sampling((5.43102, 5.43102), 0.2)
    sliced = SlicedPotential.from_snapshot(snap, grid, 5.43102 / 4)
    assert sliced.num_slices == 40
    assert sliced.unique_slice_count == 4


def _projected_multislice(phase, grid, periods, splits, depths):  # type: ignore[no-untyped-def]
    """ZOLZ multislice: one period's projected potential split into ``splits`` slices."""
    a = phase.lattice.a
    cell, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1))
    v = slice_potential(cell.species, cell.positions[:, :2], grid, "mott_bethe")
    sigma = relativistic_interaction_parameter_inv_v_angstrom(200.0)
    phase_grating = np.exp(1j * sigma * v / splits)
    tau = scipy.fft.ifft2(scipy.fft.fft2(phase_grating) * antialias_aperture(grid))
    propagator = fresnel_propagator(grid, 200.0, a / splits)
    psi = np.ones(grid.shape, dtype=complex)
    stored = []
    for step in range(periods * splits):
        psi = scipy.fft.ifft2(scipy.fft.fft2(psi * tau) * propagator)
        if (step + 1) % splits == 0 and (step + 1) // splits in depths:
            stored.append(psi)
    return stored


def test_one_slice_per_period_is_the_projected_multislice() -> None:
    phase = _phase("si_diamond")
    a = phase.lattice.a
    grid = MultisliceGrid((64, 64), (a, a))
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=10)
    wave = multislice(snap, 200.0, grid=grid, slice_thickness_angstrom=a,
                      parametrization=PotentialParametrization.MOTT_BETHE)
    (manual,) = _projected_multislice(phase, grid, 10, 1, {10})
    np.testing.assert_allclose(wave.wave(), manual, atol=1e-12)


def test_multislice_converges_to_bloch_waves_on_identical_potential() -> None:
    # Two solutions of the same Schrodinger equation from the same Mott-Bethe
    # potential, zero-order Laue zone only. One slice per period carries the
    # splitting error of a 5.4 Å slice (about 1 % in the transmitted beam at
    # 108 Å); split into eight, the multislice agrees with Bloch waves.
    phase = _phase("si_diamond")
    a = phase.lattice.a
    grid = MultisliceGrid((64, 64), (a, a))
    depths = (10, 20, 37)
    waves = _projected_multislice(phase, grid, 37, 8, set(depths))
    beams = beam_set_for_zone(
        phase, ZoneAxis(indices=(0, 0, 1), phase=phase), beam_energy_kev=200.0,
        max_index=24, g_max_inv_angstrom=3.9, max_excitation_error_inv_angstrom=0.6,
    )
    hkl = [(0, 0, 0), (2, 2, 0), (4, 0, 0)]
    for periods, psi in zip(depths, waves, strict=True):
        spectrum = scipy.fft.fft2(psi, norm="forward")
        ours = [abs(spectrum[k % 64, h % 64]) ** 2 for h, k, _ in hkl]
        bloch = solve_bloch_waves(beams, [[0.0, 0.0]], thickness_angstrom=periods * a)
        expected = [float(bloch.intensity_of(h)[0]) for h in hkl]
        np.testing.assert_allclose(ours, expected, atol=1.5e-3)


def test_forbidden_reflection_stays_forbidden_in_a_periodic_cell() -> None:
    phase = _phase("si_diamond")
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=4)
    a = phase.lattice.a
    g = np.array([[2 / a, 0.0], [2 / a, 2 / a]])
    # (200) is kinematically forbidden in diamond. With one slice per period only
    # the zero-order Laue zone couples, and no sum of allowed (hk0) reflections
    # reaches (200): it stays exactly dark. Finer slices admit the first-order
    # Laue zone, (111) + (1-1-1) = (200), and it appears by double diffraction.
    zolz = multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=a)
    holz = multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=a / 4)
    forbidden_zolz, allowed = zolz.beam_intensities(g, -1)
    forbidden_holz, _ = holz.beam_intensities(g, -1)
    assert forbidden_zolz < 1e-15
    assert 1e-9 < forbidden_holz < 1e-3 * allowed


def test_tilt_mirrors_the_pattern_of_a_mirror_symmetric_crystal() -> None:
    phase = _phase("al_fcc")
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=20)
    kwargs = {"sampling_angstrom": 0.1, "slice_thickness_angstrom": 2.02475}
    plus = multislice(snap, 200.0, tilt_mrad=(8.0, 0.0), **kwargs)
    minus = multislice(snap, 200.0, tilt_mrad=(-8.0, 0.0), **kwargs)
    level = multislice(snap, 200.0, **kwargs)
    a = phase.lattice.a
    g = np.array([[2 / a, 0.0], [-2 / a, 0.0]])
    ip, im, i0 = (w.beam_intensities(g, -1) for w in (plus, minus, level))
    assert ip[0] == pytest.approx(im[1], rel=1e-9) and ip[1] == pytest.approx(im[0], rel=1e-9)
    assert i0[0] == pytest.approx(i0[1], rel=1e-9)
    assert abs(ip[0] - ip[1]) > 1e-3


# --- Imaging and series -----------------------------------------------------


def _exit_wave():  # type: ignore[no-untyped-def]
    snap, _, _ = periodic_slab(_phase("si_diamond"), (1, 1, 0), (1, 1), thickness_angstrom=40.0)
    return multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=1.0)


def test_defocus_equals_free_space_propagation_of_the_exit_wave() -> None:
    exit_wave = _exit_wave()
    image = exit_wave.image(_lens(defocus_angstrom=-150.0))
    propagated = scipy.fft.ifft2(
        scipy.fft.fft2(exit_wave.wave())
        * np.exp(-1j * math.pi * exit_wave.wavelength_angstrom * -150.0
                 * exit_wave.grid.frequency_magnitude() ** 2)
    )
    np.testing.assert_allclose(image, np.abs(propagated) ** 2, atol=1e-10)


def test_in_focus_aberration_free_image_is_the_exit_wave_intensity() -> None:
    exit_wave = _exit_wave()
    np.testing.assert_allclose(exit_wave.image(_lens()), np.abs(exit_wave.wave()) ** 2,
                               atol=1e-10)


def test_focal_series_images_are_the_single_images() -> None:
    exit_wave = _exit_wave()
    lens = _lens(cs_mm=0.001, focal_spread_angstrom=20.0)
    series = exit_wave.focal_series(lens, [-100.0, 0.0, 100.0])
    np.testing.assert_allclose(
        series.images[2], exit_wave.image(replace(lens, defocus_angstrom=100.0)), atol=1e-12
    )
    assert series.images.shape[0] == 3
    assert "Focal series" in series.describe()


def test_focal_integration_reduces_to_the_coherent_image_without_spread() -> None:
    exit_wave = _exit_wave()
    lens = _lens(defocus_angstrom=-60.0, cs_mm=0.002)
    np.testing.assert_allclose(
        exit_wave.image(lens, temporal_coherence=TemporalCoherence.FOCAL_INTEGRATION),
        exit_wave.image(lens),
        atol=1e-12,
    )


def test_focal_integration_matches_the_envelope_only_for_a_weak_phase_object() -> None:
    # For a weak phase object the image is linear in the potential, where
    # Frank's envelope is exact, so the two treatments of focal spread agree; a
    # heavy atom is a strong object and the envelope then misstates the image.
    lens = _lens(energy_kev=300.0, defocus_angstrom=-80.0, cs_mm=0.01, focal_spread_angstrom=30.0)
    gap = {}
    for species in ("H", "Au"):
        snap = AtomicSnapshot(
            species=(species,), positions=np.array([[3.0, 3.0, 0.5]]),
            cell=np.diag([6.0, 6.0, 1.0]),
        )
        wave = multislice(snap, 300.0, sampling_angstrom=0.05, slice_thickness_angstrom=1.0)
        quasi = wave.image(lens)
        exact = wave.image(lens, temporal_coherence="focal_integration")
        gap[species] = float(np.max(np.abs(quasi - exact)) / np.ptp(quasi))
    assert gap["H"] < 0.01
    assert gap["Au"] > 0.05


def test_frozen_phonons_average_intensities_and_add_diffuse_background() -> None:
    phase = _phase("al_fcc")
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (2, 2), beam_repeats=10)
    common = {"sampling_angstrom": 0.15, "slice_thickness_angstrom": 2.02475}
    static = multislice(snap, 200.0, **common)
    thermal = multislice(
        snap, 200.0, frozen_phonon_sigma_angstrom=0.1, frozen_phonon_configurations=4,
        seed=7, **common,
    )
    again = multislice(
        snap, 200.0, frozen_phonon_sigma_angstrom=0.1, frozen_phonon_configurations=4,
        seed=7, **common,
    )
    assert thermal.configurations == 4
    np.testing.assert_array_equal(thermal.waves, again.waves)
    # Between the Bragg beams (a half-integer reciprocal-lattice position of
    # the 2 x 2 cell) the static pattern is empty and the thermal one is not.
    a = phase.lattice.a
    between = np.array([[0.5 / a, 0.0]])
    assert float(static.beam_intensities(between, -1)[0]) < 1e-20
    assert float(thermal.beam_intensities(between, -1)[0]) > 1e-8
    assert "frozen-phonon" in thermal.describe()


def test_defocus_thickness_map_rows_are_focal_series_at_each_depth() -> None:
    snap, _, _ = periodic_slab(_phase("si_diamond"), (1, 1, 0), (1, 1), thickness_angstrom=40.0)
    wave = multislice(
        snap, 200.0, sampling_angstrom=0.15, slice_thickness_angstrom=1.92,
        exit_depths_angstrom=[15.0, 40.0],
    )
    lens = _lens(cs_mm=0.001)
    table = wave.defocus_thickness_map(lens, [-50.0, 50.0])
    assert table.images.shape[:2] == (2, 2)
    np.testing.assert_allclose(
        table.images[0, 1], wave.image(replace(lens, defocus_angstrom=50.0), depth_index=0),
        atol=1e-12,
    )
    assert "Defocus-thickness map" in table.describe()


def test_hrtem_image_of_vacuum_is_uniform_unity() -> None:
    grid = MultisliceGrid((32, 32), (4.0, 4.0))
    lens = _lens(defocus_angstrom=-300.0, cs_mm=1.0, focal_spread_angstrom=40.0,
                 convergence_semiangle_mrad=0.5)
    np.testing.assert_allclose(hrtem_image(np.ones((32, 32)), grid, lens), 1.0, atol=1e-12)


def test_simulate_multislice_hrem_returns_an_explained_result() -> None:
    snap, _, _ = periodic_slab(_phase("si_diamond"), (1, 1, 0), (2, 1), thickness_angstrom=30.0)
    lens = MicroscopeAberrations.cs_corrected(energy_kev=300.0)
    result = simulate_multislice_hrem(snap, lens, sampling_angstrom=0.1)
    assert result.image.shape == result.power_spectrum.shape
    assert result.extent_angstrom == pytest.approx((2 * 5.43102, 5.43102 * math.sqrt(2)))
    assert result.exit_wave is not None
    assert "Simulated High-Resolution TEM image" in result.describe()


def test_describe_reports_the_calculation() -> None:
    text = _exit_wave().describe()
    for phrase in ("Multislice exit wave", "band limit", "slices of", "lobato"):
        assert phrase in text


def test_grid_rejects_nonsense() -> None:
    with pytest.raises(ValueError):
        MultisliceGrid((1, 10), (5.0, 5.0))
    with pytest.raises(ValueError):
        MultisliceGrid((10, 10), (0.0, 5.0))
    with pytest.raises(ValueError):
        MultisliceGrid.from_sampling((5.0, 5.0), 0.0)


def test_focal_integration_equals_a_dense_defocus_average() -> None:
    # The adaptive Gauss-Hermite rule against brute force: the image averaged
    # over a dense uniform grid of defoci with Gaussian weights.
    from pytex.diffraction.multislice import _lens_transfer

    film = AtomicSnapshot.amorphous_sample(
        species="C", density_g_cm3=2.0, dimensions_angstrom=(6.0, 6.0, 10.0), seed=5
    )
    wave = multislice(film, 300.0, sampling_angstrom=0.1, slice_thickness_angstrom=2.0)
    lens = _lens(energy_kev=300.0, defocus_angstrom=-40.0, cs_mm=0.005, focal_spread_angstrom=40.0)
    exact = wave.image(lens, temporal_coherence="focal_integration")
    spectrum = scipy.fft.fft2(wave.wave())
    total = np.zeros(wave.grid.shape)
    weights = 0.0
    for delta in np.linspace(-320.0, 320.0, 4001):
        weight = math.exp(-0.5 * (delta / 40.0) ** 2)
        focused = replace(lens, defocus_angstrom=-40.0 + delta, focal_spread_angstrom=0.0)
        image = scipy.fft.ifft2(spectrum * _lens_transfer(wave.grid, focused))
        total += weight * np.abs(image) ** 2
        weights += weight
    np.testing.assert_allclose(exact, total / weights, atol=1e-8)
