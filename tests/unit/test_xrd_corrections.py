"""Corrections and display transforms, checked against their closed forms.

The aberration tests assert the *angular signature* of each correction, not a
stored magnitude, because the signature is what makes the aberrations separable
and is the whole basis of extrapolation methods for precise lattice parameters.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_corrections import (
    PROFILE_VIEW_SCHEMA,
    ProfileView,
    correct_peak_positions,
    kalpha2_partner_two_theta_deg,
    monochromator_polarization_factor,
    position_correction_deg,
    profile_view,
    refraction_decrement,
    refraction_shift_deg,
    smooth_savitzky_golay,
    specimen_displacement_shift_deg,
    specimen_transparency_shift_deg,
    strip_kalpha2,
    variable_to_fixed_slit,
    zero_shift_deg,
)
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import detect_and_fit_peaks, detect_peaks

CU_KA1 = 1.540562
CU_KA2 = 1.544390


def _ni_scan(
    *, noisy: bool = True, seed: int = 7, peak_counts: float = 20000.0
) -> tuple[MeasuredPowderPattern, tuple[float, ...]]:
    """Return a synthetic Cu K-alpha doublet scan of nickel and its alpha1 lines."""

    pattern = generate_xrd_pattern(
        builtin_phase("ni_fcc").to_phase(),
        radiation=RadiationSpec.cu_ka_doublet(),
        two_theta_range_deg=(30.0, 150.0),
        resolution_deg=0.01,
        broadening_fwhm_deg=0.12,
        profile="pseudo_voigt",
        max_index=6,
    )
    axis = np.asarray(pattern.two_theta_grid_deg, dtype=float)
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    expected = profile / profile.max() * peak_counts + 200.0
    counts = (
        np.random.default_rng(seed).poisson(expected).astype(float) if noisy else expected
    )
    measured = MeasuredPowderPattern(
        name="Ni synthetic",
        two_theta_deg=axis,
        intensity=counts,
        radiation=RadiationSpec.cu_ka_doublet(),
        synthetic=True,
    )
    return measured, tuple(float(item.two_theta_deg) for item in pattern.reflections)


# ---------------------------------------------------------------------------
# Doublet geometry and stripping
# ---------------------------------------------------------------------------


def test_kalpha2_partner_obeys_braggs_law_at_fixed_spacing() -> None:
    angles = np.array([30.0, 60.0, 90.0, 120.0, 150.0])
    partners = kalpha2_partner_two_theta_deg(angles, wavelength_ratio=CU_KA2 / CU_KA1)
    # Same d for both lines is the defining property.
    d_one = CU_KA1 / (2.0 * np.sin(np.deg2rad(0.5 * angles)))
    d_two = CU_KA2 / (2.0 * np.sin(np.deg2rad(0.5 * partners)))
    assert np.allclose(d_one, d_two)
    # The separation grows as tan(theta), so it is negligible low and large high.
    separations = partners - angles
    assert np.all(np.diff(separations) > 0.0)
    assert separations[0] < 0.1 < separations[-1]


def test_kalpha2_partner_rejects_impossible_geometry() -> None:
    with pytest.raises(ValueError, match="wavelength_ratio"):
        kalpha2_partner_two_theta_deg([45.0], wavelength_ratio=0.0)
    with pytest.raises(ValueError, match="Ewald limit"):
        kalpha2_partner_two_theta_deg([179.9], wavelength_ratio=CU_KA2 / CU_KA1)


def test_stripping_removes_most_of_each_kalpha2_line() -> None:
    measured, expected = _ni_scan()
    stripped = strip_kalpha2(measured)
    axis = np.asarray(measured.two_theta_deg)
    before = np.asarray(measured.intensity)
    after = np.asarray(stripped.intensity)
    partners = kalpha2_partner_two_theta_deg(
        np.asarray(expected), wavelength_ratio=CU_KA2 / CU_KA1
    )
    # Only the well-separated partners are a meaningful test; below about 90
    # degrees the pair is unresolved and "the alpha2 position" has no peak.
    for partner in partners[partners > 95.0]:
        index = int(np.argmin(np.abs(axis - partner)))
        removed = 1.0 - after[index] / before[index]
        assert removed > 0.8


def test_stripping_clears_the_alpha2_lines_from_peak_detection() -> None:
    measured, expected = _ni_scan()
    stripped = strip_kalpha2(measured)
    instrument = InstrumentBroadening.ideal(0.12)
    raw = detect_peaks(measured, instrument=instrument, suppress_kalpha2=False)
    after = detect_peaks(stripped, instrument=instrument, suppress_kalpha2=False)
    partners = kalpha2_partner_two_theta_deg(
        np.asarray(expected), wavelength_ratio=CU_KA2 / CU_KA1
    )
    resolved = [float(value) for value in partners if 95.0 < value < 140.0]
    assert resolved, "the fixture must contain resolvable doublets to make this test mean anything"
    for partner in resolved:
        assert any(abs(candidate - partner) < 0.1 for candidate in raw)
        assert not any(abs(candidate - partner) < 0.1 for candidate in after)


def test_stripping_preserves_the_alpha1_positions() -> None:
    measured, expected = _ni_scan()
    stripped = strip_kalpha2(measured)
    table = detect_and_fit_peaks(
        stripped, instrument=InstrumentBroadening.ideal(0.12), model_doublet=False
    )
    for truth in expected:
        assert any(
            abs(peak.two_theta_deg - truth) < 0.01 for peak in table
        ), f"no stripped peak near the alpha1 line at {truth}"


def test_stripping_relabels_the_radiation_and_records_itself() -> None:
    measured, _ = _ni_scan()
    stripped = strip_kalpha2(measured)
    assert stripped.radiation is not None
    assert stripped.radiation.kalpha2_wavelength_angstrom is None
    assert stripped.radiation.wavelength_angstrom == pytest.approx(CU_KA1)
    assert "rachinger" in stripped.metadata["corrections"]


def test_stripping_propagates_and_inflates_the_uncertainty() -> None:
    """Each subtraction feeds the next, so the stripped pattern is noisier."""

    measured, _ = _ni_scan()
    counts = np.asarray(measured.intensity)
    with_sigma = MeasuredPowderPattern(
        name=measured.name,
        two_theta_deg=measured.two_theta_deg,
        intensity=counts,
        standard_uncertainty=np.sqrt(np.maximum(counts, 1.0)),
        radiation=measured.radiation,
        synthetic=True,
    )
    stripped = strip_kalpha2(with_sigma)
    assert stripped.standard_uncertainty is not None
    original = np.asarray(with_sigma.standard_uncertainty)
    propagated = np.asarray(stripped.standard_uncertainty)
    # Nowhere smaller, and materially larger by the high-angle end.
    assert np.all(propagated >= original - 1e-9)
    assert float(np.mean(propagated[-2000:] / original[-2000:])) > 1.05


def test_stripping_rejects_patterns_it_cannot_strip() -> None:
    measured, _ = _ni_scan()
    single = MeasuredPowderPattern(
        name="single line",
        two_theta_deg=measured.two_theta_deg,
        intensity=measured.intensity,
        radiation=RadiationSpec.cu_ka(),
        synthetic=True,
    )
    with pytest.raises(ValueError, match="declares no K-alpha2 line"):
        strip_kalpha2(single)
    bare = MeasuredPowderPattern(
        name="no radiation",
        two_theta_deg=measured.two_theta_deg,
        intensity=measured.intensity,
        synthetic=True,
    )
    with pytest.raises(ValueError, match="needs a radiation"):
        strip_kalpha2(bare)
    with pytest.raises(ValueError, match=r"intensity ratio in \(0, 1\]"):
        strip_kalpha2(measured, intensity_ratio=1.5)


# ---------------------------------------------------------------------------
# Position aberrations: each is identified by its angular signature
# ---------------------------------------------------------------------------


def test_zero_shift_is_constant_in_angle() -> None:
    angles = np.linspace(20.0, 160.0, 15)
    shift = zero_shift_deg(angles, zero_deg=0.02)
    assert np.allclose(shift, 0.02)


def test_displacement_shift_follows_cos_theta_and_vanishes_at_backscatter() -> None:
    angles = np.linspace(10.0, 179.0, 40)
    shift = specimen_displacement_shift_deg(
        angles, displacement_mm=0.05, goniometer_radius_mm=240.0
    )
    theta = np.deg2rad(0.5 * angles)
    predicted = np.rad2deg(-2.0 * 0.05 * np.cos(theta) / 240.0)
    assert np.allclose(shift, predicted)
    # This is the reason extrapolation to theta = 90 degrees works at all.
    assert abs(shift[-1]) < abs(shift[0]) / 50.0
    assert np.all(np.diff(np.abs(shift)) < 0.0)


def test_displacement_shift_is_linear_in_the_displacement() -> None:
    angles = np.linspace(20.0, 160.0, 20)
    one = specimen_displacement_shift_deg(
        angles, displacement_mm=0.05, goniometer_radius_mm=240.0
    )
    two = specimen_displacement_shift_deg(
        angles, displacement_mm=0.10, goniometer_radius_mm=240.0
    )
    assert np.allclose(two, 2.0 * one)


def test_transparency_shift_peaks_at_forty_five_degrees_theta() -> None:
    angles = np.linspace(5.0, 175.0, 341)
    shift = specimen_transparency_shift_deg(
        angles, linear_absorption_coefficient_inv_mm=2.0, goniometer_radius_mm=240.0
    )
    # sin(2 theta) is extremal at 2 theta = 90 degrees.
    assert angles[int(np.argmax(np.abs(shift)))] == pytest.approx(90.0, abs=1.0)
    assert abs(shift[0]) < abs(shift[int(np.argmax(np.abs(shift)))])


def test_transparency_shift_scales_inversely_with_absorption() -> None:
    angles = np.linspace(20.0, 160.0, 20)
    weak = specimen_transparency_shift_deg(
        angles, linear_absorption_coefficient_inv_mm=1.0, goniometer_radius_mm=240.0
    )
    strong = specimen_transparency_shift_deg(
        angles, linear_absorption_coefficient_inv_mm=100.0, goniometer_radius_mm=240.0
    )
    assert np.allclose(strong, weak / 100.0)


def test_refraction_decrement_has_the_expected_order_of_magnitude() -> None:
    # Nickel at Cu K-alpha: delta is a few parts in 1e5, which is exactly the
    # precision level a Cohen determination reaches, so it is not negligible.
    delta = refraction_decrement(density_g_cm3=8.9, wavelength_angstrom=CU_KA1)
    assert 1.0e-5 < delta < 1.0e-4
    # delta scales as rho lambda^2.
    doubled = refraction_decrement(density_g_cm3=17.8, wavelength_angstrom=CU_KA1)
    assert doubled == pytest.approx(2.0 * delta)
    quadrupled = refraction_decrement(density_g_cm3=8.9, wavelength_angstrom=2.0 * CU_KA1)
    assert quadrupled == pytest.approx(4.0 * delta)


def test_refraction_shift_grows_towards_low_angle() -> None:
    angles = np.linspace(20.0, 160.0, 30)
    shift = refraction_shift_deg(angles, decrement=3.0e-5)
    assert np.all(shift > 0.0)
    # 1 / sin(2 theta) is symmetric about 90 degrees and minimal there.
    assert shift[int(np.argmin(shift))] == pytest.approx(shift.min())
    assert angles[int(np.argmin(shift))] == pytest.approx(90.0, abs=5.0)


def test_aberrations_compose_additively() -> None:
    angles = np.linspace(20.0, 160.0, 25)
    total = position_correction_deg(
        angles,
        zero_deg=0.01,
        displacement_mm=0.05,
        goniometer_radius_mm=240.0,
        linear_absorption_coefficient_inv_mm=5.0,
        refraction_decrement_value=3.0e-5,
    )
    parts = (
        zero_shift_deg(angles, zero_deg=0.01)
        + specimen_displacement_shift_deg(
            angles, displacement_mm=0.05, goniometer_radius_mm=240.0
        )
        + specimen_transparency_shift_deg(
            angles, linear_absorption_coefficient_inv_mm=5.0, goniometer_radius_mm=240.0
        )
        + refraction_shift_deg(angles, decrement=3.0e-5)
    )
    assert np.allclose(total, parts)


def test_correcting_a_peak_table_undoes_a_known_displacement() -> None:
    """Inject a displacement into the positions, then remove it exactly."""

    measured, expected = _ni_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(0.12))
    displaced_axis = np.asarray(measured.two_theta_deg) - specimen_displacement_shift_deg(
        measured.two_theta_deg, displacement_mm=0.08, goniometer_radius_mm=240.0
    )
    displaced = MeasuredPowderPattern(
        name="displaced",
        two_theta_deg=displaced_axis,
        intensity=measured.intensity,
        radiation=measured.radiation,
        synthetic=True,
    )
    shifted_table = detect_and_fit_peaks(displaced, instrument=InstrumentBroadening.ideal(0.12))
    # The displacement is large enough to matter before correction.
    assert np.max(np.abs(shifted_table.two_theta_deg - table.two_theta_deg)) > 0.01
    restored = correct_peak_positions(
        shifted_table, displacement_mm=-0.08, goniometer_radius_mm=240.0
    )
    for peak, truth in zip(restored, expected, strict=True):
        assert peak.two_theta_deg == pytest.approx(truth, abs=2.0e-3)


def test_correcting_positions_leaves_the_uncertainties_alone() -> None:
    """A known aberration is a bias, not a random error."""

    measured, _ = _ni_scan()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(0.12))
    corrected = correct_peak_positions(table, zero_deg=0.03)
    assert np.allclose(corrected.standard_uncertainty_deg, table.standard_uncertainty_deg)
    assert np.allclose(corrected.two_theta_deg, table.two_theta_deg - 0.03)
    assert corrected.settings["position_correction_zero_deg"] == pytest.approx(0.03)


def test_aberrations_reject_impossible_geometry() -> None:
    with pytest.raises(ValueError, match="goniometer_radius_mm"):
        specimen_displacement_shift_deg([45.0], displacement_mm=0.05, goniometer_radius_mm=0.0)
    with pytest.raises(ValueError, match="linear_absorption_coefficient_inv_mm"):
        specimen_transparency_shift_deg(
            [45.0], linear_absorption_coefficient_inv_mm=0.0, goniometer_radius_mm=240.0
        )
    with pytest.raises(ValueError, match="density_g_cm3"):
        refraction_decrement(density_g_cm3=0.0, wavelength_angstrom=CU_KA1)
    with pytest.raises(ValueError, match="non-negative decrement"):
        refraction_shift_deg([45.0], decrement=-1.0)
    with pytest.raises(ValueError, match="diverges"):
        refraction_shift_deg([0.0], decrement=1e-5)


# ---------------------------------------------------------------------------
# Intensity corrections
# ---------------------------------------------------------------------------


def test_monochromator_polarization_matches_its_limiting_cases() -> None:
    angles = np.linspace(0.0, 180.0, 19)
    # A monochromator at 2 theta_M = 0 does not polarize: recover the
    # unpolarized (1 + cos^2 2theta) / 2 factor.
    unpolarized = monochromator_polarization_factor(angles, monochromator_two_theta_deg=0.0)
    assert np.allclose(unpolarized, 0.5 * (1.0 + np.square(np.cos(np.deg2rad(angles)))))
    # An ideal polarizer at 90 degrees passes only the component the specimen
    # scatters without angular dependence, so in the parallel setting the
    # polarization correction disappears altogether.
    ideal_parallel = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=90.0
    )
    assert np.allclose(ideal_parallel, 1.0)
    # In the perpendicular setting the surviving component is the one the
    # specimen modulates fully, so the factor is cos^2(2 theta).
    ideal_perpendicular = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=90.0, perpendicular=True
    )
    assert np.allclose(ideal_perpendicular, np.square(np.cos(np.deg2rad(angles))))
    # Normalized to one in the forward direction in every configuration.
    for setting in (0.0, 26.6, 90.0):
        assert monochromator_polarization_factor(
            [0.0], monochromator_two_theta_deg=setting
        )[0] == pytest.approx(1.0)


def test_perpendicular_geometry_exchanges_the_two_cosines() -> None:
    angles = np.array([0.0, 45.0, 90.0])
    parallel = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=30.0, perpendicular=False
    )
    perpendicular = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=30.0, perpendicular=True
    )
    assert not np.allclose(parallel, perpendicular)
    # (1 + m c) - (m + c) = (1 - m)(1 - c) >= 0 for m, c in [0, 1], so the
    # parallel factor is never below the perpendicular one, and they meet only
    # where one of the two cosines squared reaches one.
    assert np.all(parallel >= perpendicular - 1e-12)
    assert parallel[0] == pytest.approx(perpendicular[0])  # 2 theta = 0, c = 1
    assert parallel[1] > perpendicular[1]
    unpolarized_parallel = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=0.0
    )
    unpolarized_perpendicular = monochromator_polarization_factor(
        angles, monochromator_two_theta_deg=0.0, perpendicular=True
    )
    assert np.allclose(unpolarized_parallel, unpolarized_perpendicular)  # m = 1


def test_monochromator_angle_is_range_checked() -> None:
    with pytest.raises(ValueError, match=r"\[0, 180\]"):
        monochromator_polarization_factor([45.0], monochromator_two_theta_deg=200.0)


def test_variable_to_fixed_slit_suppresses_the_low_angle_end() -> None:
    measured, _ = _ni_scan()
    converted = variable_to_fixed_slit(measured)
    axis = np.asarray(measured.two_theta_deg)
    ratio = np.asarray(converted.intensity) / np.maximum(np.asarray(measured.intensity), 1.0)
    assert np.allclose(ratio, np.sin(np.deg2rad(0.5 * axis)), atol=1e-9)
    assert converted.intensity_unit == "arbitrary"
    assert "variable_to_fixed_slit" in converted.metadata["corrections"]


def test_smoothing_broadens_peaks_and_drops_the_uncertainties() -> None:
    """The docstring's warning is a testable claim, so it is tested."""

    measured, _ = _ni_scan(peak_counts=800.0, seed=5)
    instrument = InstrumentBroadening.ideal(0.12)
    raw_widths = detect_and_fit_peaks(measured, instrument=instrument).peaks
    smoothed = smooth_savitzky_golay(measured, window_points=31, polynomial_order=3)
    smooth_widths = detect_and_fit_peaks(smoothed, instrument=instrument).peaks
    assert smoothed.standard_uncertainty is None
    raw_median = float(np.median([peak.fwhm_deg for peak in raw_widths]))
    smooth_median = float(np.median([peak.fwhm_deg for peak in smooth_widths]))
    assert smooth_median > raw_median


def test_smoothing_rejects_invalid_windows() -> None:
    measured, _ = _ni_scan()
    with pytest.raises(ValueError, match="odd window_points"):
        smooth_savitzky_golay(measured, window_points=10)
    with pytest.raises(ValueError, match="window_points > polynomial_order"):
        smooth_savitzky_golay(measured, window_points=3, polynomial_order=5)
    too_long = len(measured) + 1 + (len(measured) % 2)
    assert too_long % 2 == 1 and too_long > len(measured)
    with pytest.raises(ValueError, match="longer than the measured pattern"):
        smooth_savitzky_golay(measured, window_points=too_long)


# ---------------------------------------------------------------------------
# Display transforms
# ---------------------------------------------------------------------------


def test_abscissa_transforms_are_the_textbook_identities() -> None:
    measured, _ = _ni_scan()
    angles = np.asarray(measured.two_theta_deg)
    sine = np.sin(np.deg2rad(0.5 * angles))

    assert np.allclose(profile_view(measured, abscissa="two_theta_deg").abscissa, angles)
    assert np.allclose(
        profile_view(measured, abscissa="sin_squared_theta").abscissa, np.square(sine)
    )
    spacing = profile_view(measured, abscissa="d_angstrom").abscissa
    assert np.allclose(spacing, CU_KA1 / (2.0 * sine))
    scattering = profile_view(measured, abscissa="q_inv_angstrom").abscissa
    assert np.allclose(scattering, 4.0 * np.pi * sine / CU_KA1)
    # Q = 2 pi / d ties the two wavelength-dependent abscissae together.
    assert np.allclose(scattering, 2.0 * np.pi / spacing)
    # d runs backwards relative to angle.
    assert np.all(np.diff(spacing) < 0.0)


def test_square_root_scale_equalises_the_counting_noise() -> None:
    """The point of the sqrt ordinate is statistical, not aesthetic."""

    rng = np.random.default_rng(3)
    axis = np.linspace(20.0, 120.0, 20001)
    level = np.where(axis < 70.0, 100.0, 10000.0)
    pattern = MeasuredPowderPattern(
        name="two count levels",
        two_theta_deg=axis,
        intensity=rng.poisson(level).astype(float),
        radiation=RadiationSpec.cu_ka(),
        synthetic=True,
    )
    linear = profile_view(pattern, scale="linear").ordinate
    root = profile_view(pattern, scale="sqrt").ordinate
    low = axis < 70.0
    # On a linear ordinate the noise amplitude differs by a factor of ten.
    assert float(np.std(linear[~low]) / np.std(linear[low])) > 5.0
    # On the variance-stabilized ordinate it is the same everywhere.
    assert float(np.std(root[~low]) / np.std(root[low])) == pytest.approx(1.0, abs=0.2)


def test_logarithmic_scale_clips_rather_than_producing_negative_infinity() -> None:
    axis = np.linspace(20.0, 60.0, 501)
    intensity = np.zeros_like(axis)
    intensity[250] = 1000.0
    pattern = MeasuredPowderPattern(
        name="one spike",
        two_theta_deg=axis,
        intensity=intensity,
        radiation=RadiationSpec.cu_ka(),
        synthetic=True,
    )
    view = profile_view(pattern, scale="log10", log_floor_fraction=1.0e-3)
    assert np.all(np.isfinite(view.ordinate))
    assert float(view.ordinate.min()) == pytest.approx(np.log10(1.0))
    assert float(view.ordinate.max()) == pytest.approx(3.0)
    assert "clipped at" in view.ordinate_label


def test_normalizations_do_what_they_say() -> None:
    measured, _ = _ni_scan()
    maximum = profile_view(measured, normalization="maximum")
    assert float(maximum.ordinate.max()) == pytest.approx(1.0)
    integral = profile_view(measured, normalization="integral")
    area = float(np.trapezoid(integral.ordinate, integral.abscissa))
    assert area == pytest.approx(1.0, rel=1e-9)


def test_profile_view_contract_and_prose_agree() -> None:
    measured, _ = _ni_scan()
    view = profile_view(
        measured, abscissa="q_inv_angstrom", scale="sqrt", normalization="maximum"
    )
    payload = view.to_json()
    assert payload["schema"] == PROFILE_VIEW_SCHEMA
    assert payload["abscissa_kind"] == "q_inv_angstrom"
    assert len(payload["abscissa"]) == len(view) == len(payload["ordinate"])
    prose = view.describe()
    assert "wavelength-free abscissa" in prose
    assert "variance-stabilizing" in prose
    assert "strongest point is one" in prose


def test_profile_view_validates_its_inputs() -> None:
    measured, _ = _ni_scan()
    with pytest.raises(ValueError, match="abscissa in"):
        profile_view(measured, abscissa="theta")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="scale in"):
        profile_view(measured, scale="ln")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="normalization in"):
        profile_view(measured, normalization="unit")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="log_floor_fraction"):
        profile_view(measured, scale="log10", log_floor_fraction=2.0)
    bare = MeasuredPowderPattern(
        name="no radiation",
        two_theta_deg=measured.two_theta_deg,
        intensity=measured.intensity,
        synthetic=True,
    )
    with pytest.raises(ValueError, match="needs a wavelength"):
        profile_view(bare, abscissa="d_angstrom")


def test_profile_view_rejects_mismatched_arrays() -> None:
    with pytest.raises(ValueError, match="same shape"):
        ProfileView(
            name="bad",
            abscissa=np.array([1.0, 2.0, 3.0]),
            abscissa_kind="two_theta_deg",
            abscissa_label="2*theta",
            ordinate=np.array([1.0, 2.0]),
            ordinate_label="intensity",
            scale="linear",
            normalization="none",
            source_name="source",
        )
    with pytest.raises(ValueError, match="abscissa_kind"):
        ProfileView(
            name="bad",
            abscissa=np.array([1.0, 2.0]),
            abscissa_kind="theta",  # type: ignore[arg-type]
            abscissa_label="theta",
            ordinate=np.array([1.0, 2.0]),
            ordinate_label="intensity",
            scale="linear",
            normalization="none",
            source_name="source",
        )
