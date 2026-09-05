"""Peak detection and single-peak fitting on patterns whose answer is known.

Every position assertion here is checked against the reflection list that
generated the synthetic pattern, not against a stored output of this code, so a
regression in the fit shows up as a disagreement with Bragg's law rather than
with a previous run.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import (
    PEAK_FIT_SCHEMA,
    PEAK_TABLE_SCHEMA,
    PeakFit,
    PeakTable,
    detect_and_fit_peaks,
    detect_peaks,
    fit_peaks,
    pseudo_voigt_area,
    pseudo_voigt_profile,
    split_pseudo_voigt_profile,
)

TEST_FWHM_DEG = 0.12


def _synthetic_scan(
    *,
    identifier: str = "ni_fcc",
    radiation: RadiationSpec | None = None,
    two_theta_range_deg: tuple[float, float] = (30.0, 150.0),
    peak_counts: float = 20000.0,
    background_counts: float = 200.0,
    seed: int = 7,
    resolution_deg: float = 0.01,
) -> tuple[MeasuredPowderPattern, tuple[float, ...]]:
    """Return a Poisson-noised synthetic scan and its true K-alpha1 positions."""

    phase = builtin_phase(identifier).to_phase()
    spec = radiation if radiation is not None else RadiationSpec.cu_ka_doublet()
    pattern = generate_xrd_pattern(
        phase,
        radiation=spec,
        two_theta_range_deg=two_theta_range_deg,
        resolution_deg=resolution_deg,
        broadening_fwhm_deg=TEST_FWHM_DEG,
        profile="pseudo_voigt",
        max_index=6,
    )
    axis = np.asarray(pattern.two_theta_grid_deg, dtype=float)
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    expected = profile / profile.max() * peak_counts + background_counts
    counts = np.random.default_rng(seed).poisson(expected).astype(float)
    measured = MeasuredPowderPattern(
        name=f"{identifier} synthetic",
        two_theta_deg=axis,
        intensity=counts,
        radiation=spec,
        synthetic=True,
    )
    return measured, tuple(float(item.two_theta_deg) for item in pattern.reflections)


# ---------------------------------------------------------------------------
# Profile shapes
# ---------------------------------------------------------------------------


def test_pseudo_voigt_is_unit_height_and_has_the_requested_fwhm() -> None:
    axis = np.linspace(-2.0, 2.0, 40001)
    for eta in (0.0, 0.37, 1.0):
        profile = pseudo_voigt_profile(axis, centre_deg=0.0, fwhm_deg=0.5, eta=eta)
        assert profile.max() == pytest.approx(1.0, abs=1e-12)
        # The half-maximum crossings are one FWHM apart, by construction of
        # both components sharing the width.
        above = axis[profile >= 0.5]
        assert float(above[-1] - above[0]) == pytest.approx(0.5, abs=2e-4)


@pytest.mark.parametrize("eta", [0.0, 0.5, 1.0])
def test_pseudo_voigt_area_matches_numerical_quadrature(eta: float) -> None:
    axis = np.linspace(-600.0, 600.0, 4_000_001)
    profile = pseudo_voigt_profile(axis, centre_deg=0.0, fwhm_deg=0.4, eta=eta)
    numerical = float(np.trapezoid(profile, axis))
    analytic = pseudo_voigt_area(height=1.0, fwhm_deg=0.4, eta=eta)
    assert numerical == pytest.approx(analytic, rel=2e-3)


def test_split_pseudo_voigt_is_continuous_and_asymmetric() -> None:
    # Pure Gaussian halves, so the tails are enclosed by a finite window and
    # the area ratio is exactly the width ratio with no truncation bias.
    axis = np.linspace(-5.0, 5.0, 200001)
    profile = split_pseudo_voigt_profile(
        axis, centre_deg=0.0, fwhm_left_deg=0.2, fwhm_right_deg=0.4, eta=0.0
    )
    assert profile.max() == pytest.approx(1.0, abs=1e-12)
    left_area = float(np.trapezoid(profile[axis < 0.0], axis[axis < 0.0]))
    right_area = float(np.trapezoid(profile[axis > 0.0], axis[axis > 0.0]))
    # The tolerance covers the half-cell each half-integral drops at the join.
    assert right_area == pytest.approx(2.0 * left_area, rel=1e-3)

    # Continuity at the join: both halves equal one at the centre.
    join = split_pseudo_voigt_profile(
        np.array([-1e-9, 0.0, 1e-9]),
        centre_deg=0.0,
        fwhm_left_deg=0.2,
        fwhm_right_deg=0.4,
        eta=0.4,
    )
    assert np.allclose(join, 1.0, atol=1e-12)


def test_profile_rejects_invalid_shape_parameters() -> None:
    with pytest.raises(ValueError, match="positive fwhm_deg"):
        pseudo_voigt_profile(np.zeros(3), centre_deg=0.0, fwhm_deg=0.0, eta=0.5)
    with pytest.raises(ValueError, match=r"eta in \[0, 1\]"):
        pseudo_voigt_profile(np.zeros(3), centre_deg=0.0, fwhm_deg=0.1, eta=1.5)


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def test_detection_finds_every_reflection_and_no_others() -> None:
    measured, expected = _synthetic_scan()
    found = detect_peaks(
        measured, instrument=InstrumentBroadening.ideal(TEST_FWHM_DEG), prominence_sigma=5.0
    )
    assert len(found) == len(expected)
    for candidate, truth in zip(found, expected, strict=True):
        assert candidate == pytest.approx(truth, abs=0.05)


def test_detection_suppresses_the_resolved_kalpha2_partners() -> None:
    """Above ~90 degrees the doublet resolves; alpha2 must not become a peak."""

    measured, expected = _synthetic_scan()
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    suppressed = detect_peaks(measured, instrument=instrument, suppress_kalpha2=True)
    raw = detect_peaks(measured, instrument=instrument, suppress_kalpha2=False)
    assert len(raw) > len(expected), "the alpha2 lines should be detectable at all"
    assert len(suppressed) == len(expected)

    # Each extra raw candidate is at a K-alpha2 position of a real reflection.
    ratio = 1.544390 / 1.540562
    partners = [
        float(np.rad2deg(2.0 * np.arcsin(ratio * np.sin(np.deg2rad(0.5 * position)))))
        for position in expected
    ]
    for candidate in raw:
        near_real = min(abs(candidate - value) for value in expected)
        near_partner = min(abs(candidate - value) for value in partners)
        assert min(near_real, near_partner) < 0.05


def test_detection_threshold_is_expressed_in_noise_sigmas() -> None:
    measured, expected = _synthetic_scan(peak_counts=400.0, background_counts=120.0, seed=3)
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    strict = detect_peaks(measured, instrument=instrument, prominence_sigma=30.0)
    permissive = detect_peaks(measured, instrument=instrument, prominence_sigma=4.0)
    assert len(strict) <= len(permissive)
    assert len(permissive) >= len(expected) - 2


def test_detection_honours_the_angular_window_and_the_peak_cap() -> None:
    measured, _ = _synthetic_scan()
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    windowed = detect_peaks(
        measured, instrument=instrument, two_theta_range_deg=(60.0, 100.0)
    )
    assert windowed
    assert all(60.0 <= value <= 100.0 for value in windowed)
    capped = detect_peaks(measured, instrument=instrument, max_peaks=2)
    assert len(capped) == 2


def test_detection_rejects_invalid_settings() -> None:
    measured, _ = _synthetic_scan()
    with pytest.raises(ValueError, match="positive prominence_sigma"):
        detect_peaks(measured, prominence_sigma=0.0)
    with pytest.raises(ValueError, match="max_peaks"):
        detect_peaks(measured, max_peaks=0)
    with pytest.raises(ValueError, match="increasing"):
        detect_peaks(measured, two_theta_range_deg=(90.0, 40.0))
    with pytest.raises(ValueError, match="sixteen points"):
        detect_peaks(measured, two_theta_range_deg=(60.0, 60.05))


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------


def test_fitted_positions_recover_the_generating_reflections() -> None:
    measured, expected = _synthetic_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(TEST_FWHM_DEG))
    assert len(table) == len(expected)
    assert table.converged_count == len(table)
    for peak, truth in zip(table, expected, strict=True):
        # One millidegree is roughly five standard uncertainties here; the fit
        # is limited by counting statistics, not by the model.
        assert peak.two_theta_deg == pytest.approx(truth, abs=1.0e-3)
        assert 0.0 < peak.two_theta_standard_uncertainty_deg < 1.0e-3


def test_fit_uncertainties_shrink_as_counting_statistics_improve() -> None:
    """Ten times the counts should roughly halve the position uncertainty."""

    weak, _ = _synthetic_scan(peak_counts=2_000.0, seed=11)
    strong, _ = _synthetic_scan(peak_counts=200_000.0, seed=11)
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    weak_table = detect_and_fit_peaks(weak, instrument=instrument)
    strong_table = detect_and_fit_peaks(strong, instrument=instrument)
    weak_median = float(np.median(weak_table.standard_uncertainty_deg))
    strong_median = float(np.median(strong_table.standard_uncertainty_deg))
    assert strong_median < weak_median
    # sigma scales as 1/sqrt(N); a hundredfold count increase is a factor ten.
    assert strong_median == pytest.approx(weak_median / 10.0, rel=0.6)


def test_reduced_chi_squared_is_near_unity_for_a_correct_model() -> None:
    measured, _ = _synthetic_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(TEST_FWHM_DEG))
    values = np.array([peak.reduced_chi_squared for peak in table])
    assert np.all(values > 0.5)
    assert np.all(values < 2.0)


def test_modelling_the_doublet_beats_ignoring_it() -> None:
    """A single line fitted to an unresolved doublet sits at the wrong angle."""

    measured, expected = _synthetic_scan(two_theta_range_deg=(40.0, 60.0))
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    candidates = detect_peaks(measured, instrument=instrument)
    modelled = fit_peaks(measured, candidates, instrument=instrument, model_doublet=True)
    ignored = fit_peaks(measured, candidates, instrument=instrument, model_doublet=False)
    for peak, truth in zip(modelled, expected, strict=False):
        assert peak.two_theta_deg == pytest.approx(truth, abs=1.0e-3)
    # Ignoring the alpha2 line pulls every centre to higher angle, because the
    # partner always sits at higher 2*theta.
    shifts = [
        wrong.two_theta_deg - right.two_theta_deg
        for wrong, right in zip(ignored, modelled, strict=True)
    ]
    assert all(shift > 5.0e-3 for shift in shifts)


def test_split_profile_fits_and_reports_its_asymmetry() -> None:
    measured, expected = _synthetic_scan()
    instrument = InstrumentBroadening.ideal(TEST_FWHM_DEG)
    candidates = detect_peaks(measured, instrument=instrument)
    table = fit_peaks(
        measured, candidates, instrument=instrument, shape="split_pseudo_voigt"
    )
    assert all(peak.shape == "split_pseudo_voigt" for peak in table)
    # The generator is symmetric, so the fitted asymmetry must stay near one.
    for peak, truth in zip(table, expected, strict=True):
        assert peak.asymmetry == pytest.approx(1.0, abs=0.25)
        assert peak.two_theta_deg == pytest.approx(truth, abs=3.0e-3)


def test_bragg_conversion_round_trips_through_the_peak_position() -> None:
    measured, expected = _synthetic_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(TEST_FWHM_DEG))
    wavelength = 1.540562
    spacings = table.d_spacing_angstrom()
    for spacing, truth in zip(spacings, expected, strict=True):
        angle = 2.0 * np.rad2deg(np.arcsin(wavelength / (2.0 * spacing)))
        assert angle == pytest.approx(truth, abs=1.0e-3)


def test_fitting_rejects_impossible_requests() -> None:
    measured, _ = _synthetic_scan()
    with pytest.raises(ValueError, match="shape in"):
        fit_peaks(measured, [45.0], shape="voigt")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="positive window_fwhm"):
        fit_peaks(measured, [45.0], window_fwhm=0.0)
    with pytest.raises(ValueError, match="no candidate positions"):
        fit_peaks(measured, [])
    with pytest.raises(ValueError, match="enough measured points"):
        fit_peaks(measured, [400.0], expected_fwhm_deg=0.1)


# ---------------------------------------------------------------------------
# Result objects
# ---------------------------------------------------------------------------


def _one_fit(**overrides: object) -> PeakFit:
    values: dict[str, object] = {
        "two_theta_deg": 44.5,
        "two_theta_standard_uncertainty_deg": 2.0e-4,
        "height": 1000.0,
        "integrated_intensity": 130.0,
        "fwhm_deg": 0.12,
        "fwhm_left_deg": 0.12,
        "fwhm_right_deg": 0.12,
        "eta": 0.4,
        "shape": "pseudo_voigt",
        "doublet_modelled": True,
        "background_intercept": 200.0,
        "background_slope": 0.0,
        "reduced_chi_squared": 1.05,
        "point_count": 96,
        "window_deg": (44.0, 45.0),
        "converged": True,
    }
    values.update(overrides)
    return PeakFit(**values)  # type: ignore[arg-type]


def test_peak_fit_validates_its_invariants() -> None:
    with pytest.raises(ValueError, match="strictly inside"):
        _one_fit(two_theta_deg=181.0)
    with pytest.raises(ValueError, match="standard_uncertainty"):
        _one_fit(two_theta_standard_uncertainty_deg=0.0)
    with pytest.raises(ValueError, match="half-widths"):
        _one_fit(fwhm_left_deg=0.0)
    with pytest.raises(ValueError, match=r"eta must lie in \[0, 1\]"):
        _one_fit(eta=2.0)
    with pytest.raises(ValueError, match="point_count"):
        _one_fit(point_count=0)
    with pytest.raises(ValueError, match="shape must be one of"):
        _one_fit(shape="voigt")


def test_peak_table_sorts_by_angle_and_reports_convergence() -> None:
    table = PeakTable(
        name="table",
        peaks=(_one_fit(two_theta_deg=90.0, converged=False), _one_fit(two_theta_deg=44.5)),
        radiation=RadiationSpec.cu_ka_doublet(),
        source_name="scan",
    )
    assert [peak.two_theta_deg for peak in table] == [44.5, 90.0]
    assert len(table) == 2
    assert table.converged_count == 1
    assert len(table.filter_converged()) == 1
    with pytest.raises(ValueError, match="No peak in this table converged"):
        PeakTable(name="none", peaks=(_one_fit(converged=False),)).filter_converged()


def test_peak_table_requires_radiation_to_produce_spacings() -> None:
    table = PeakTable(name="table", peaks=(_one_fit(),))
    with pytest.raises(ValueError, match="requires a radiation"):
        table.d_spacing_angstrom()


def test_json_contracts_carry_their_schema_and_match_describe() -> None:
    measured, _ = _synthetic_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(TEST_FWHM_DEG))
    payload = table.to_json()
    assert payload["schema"] == PEAK_TABLE_SCHEMA
    assert len(payload["peaks"]) == len(table)
    assert payload["peaks"][0]["schema"] == PEAK_FIT_SCHEMA
    assert payload["radiation"]["kalpha2_wavelength_angstrom"] == pytest.approx(1.544390)

    prose = table.describe()
    assert f"{len(table)} fitted reflections" in prose
    assert "Cullity" in prose
    assert table[0].describe().startswith("Reflection fitted at")


def test_describe_names_the_doublet_treatment_it_actually_used() -> None:
    modelled = _one_fit(doublet_modelled=True).describe()
    ignored = _one_fit(doublet_modelled=False).describe()
    assert "modelled jointly" in modelled
    assert "must not be used for precise parameters" in ignored


def test_empty_detection_is_reported_rather_than_faked() -> None:
    axis = np.linspace(20.0, 80.0, 3000)
    flat = MeasuredPowderPattern(
        name="amorphous",
        two_theta_deg=axis,
        intensity=np.full_like(axis, 500.0),
        radiation=RadiationSpec.cu_ka_doublet(),
        synthetic=True,
    )
    assert detect_peaks(flat, expected_fwhm_deg=0.12, prominence_sigma=8.0) == ()
    with pytest.raises(ValueError, match="No reflection was detected"):
        detect_and_fit_peaks(flat, expected_fwhm_deg=0.12, prominence_sigma=8.0)
