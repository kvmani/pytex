"""Background estimation, instrumental broadening, and whole-profile refinement.

The load-bearing tests here are the *recovery* tests: a pattern is built from a
known cell, zero shift, width and texture, and the refinement is required to
find those values back. A refinement that merely converges to a low ``R_wp``
proves nothing, because a wrong model with enough background coefficients can
do that too.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    InstrumentBroadening,
    MeasuredPowderPattern,
    PowderBackground,
    RadiationSpec,
    calibrate_instrument_broadening,
    deconvolve_instrument_width,
    estimate_background,
    generate_xrd_pattern,
    refine_rietveld,
    scherrer_size_nm,
    williamson_hall,
)
from pytex.app.phases import phase_from_request
from pytex.diffraction.rietveld import _scaled_phase

CU_KA = RadiationSpec.cu_ka()


def _silicon():
    _, phase = phase_from_request({"builtin": "si_diamond"})
    return phase


def _synthetic_measurement(
    *,
    lattice_scale: float = 1.0,
    zero_shift_deg: float = 0.0,
    fwhm_deg: float = 0.14,
    background_level: float = 150.0,
    peak_counts: float = 20000.0,
    seed: int | None = 7,
    preferred_orientation=None,
) -> MeasuredPowderPattern:
    """Build a measured-looking pattern from a known, deliberately perturbed model."""

    phase = _scaled_phase(_silicon(), lattice_scale)
    pattern = generate_xrd_pattern(
        phase,
        radiation=CU_KA,
        two_theta_range_deg=(20.0, 120.0),
        resolution_deg=0.02,
        broadening_fwhm_deg=fwhm_deg,
        intensity_model="xray_tabulated",
        preferred_orientation=preferred_orientation,
    )
    angles = pattern.two_theta_grid_deg + zero_shift_deg
    counts = peak_counts * pattern.intensity_grid + background_level
    counts = counts + 60.0 * np.exp(-0.5 * ((angles - 25.0) / 9.0) ** 2)
    if seed is not None:
        counts = np.random.default_rng(seed).poisson(counts).astype(float)
    return MeasuredPowderPattern(
        name="synthetic silicon",
        two_theta_deg=angles,
        intensity=counts,
        radiation=CU_KA,
        synthetic=True,
    )


class TestBackgroundEstimation:
    def test_snip_recovers_a_flat_background_under_sharp_peaks(self) -> None:
        angles = np.linspace(20.0, 80.0, 601)
        peaks = sum(
            120.0 * np.exp(-0.5 * ((angles - centre) / 0.15) ** 2)
            for centre in (33.0, 45.0, 58.0, 72.0)
        )
        measured = MeasuredPowderPattern(
            name="flat", two_theta_deg=angles, intensity=peaks + 10.0, synthetic=True
        )
        estimate = estimate_background(measured, half_window_deg=2.0)
        assert np.isclose(float(np.median(estimate.background)), 10.0, atol=0.5)

    @pytest.mark.parametrize(
        ("method", "kwargs"),
        (("snip", {"half_window_deg": 2.0}), ("chebyshev", {"degree": 6})),
    )
    def test_both_methods_follow_a_curved_background(self, method, kwargs) -> None:
        angles = np.linspace(20.0, 80.0, 601)
        truth = 10.0 + 40.0 * np.exp(-0.5 * ((angles - 30.0) / 8.0) ** 2)
        peaks = sum(
            140.0 * np.exp(-0.5 * ((angles - centre) / 0.15) ** 2)
            for centre in (33.0, 45.0, 58.0, 72.0)
        )
        measured = MeasuredPowderPattern(
            name="curved", two_theta_deg=angles, intensity=truth + peaks, synthetic=True
        )
        estimate = estimate_background(measured, method=method, **kwargs)
        # Within a tenth of the background's own amplitude: close enough to be
        # useful, loose enough not to pin an estimator to its current tuning.
        assert float(np.mean(np.abs(estimate.background - truth))) < 4.0
        assert estimate.method == method

    def test_subtraction_is_non_negative_and_records_itself(self) -> None:
        measured = _synthetic_measurement(seed=None)
        estimate = estimate_background(measured, half_window_deg=2.0)
        subtracted = estimate.subtract(measured)
        assert np.all(subtracted.intensity >= 0.0)
        assert subtracted.metadata["background_subtracted"] == "snip"
        assert 0.0 < estimate.background_fraction < 1.0

    def test_a_background_cannot_be_subtracted_from_a_different_pattern(self) -> None:
        measured = _synthetic_measurement(seed=None)
        estimate = estimate_background(measured, half_window_deg=2.0)
        other = MeasuredPowderPattern(
            name="other",
            two_theta_deg=measured.two_theta_deg + 5.0,
            intensity=measured.intensity,
            synthetic=True,
        )
        with pytest.raises(ValueError, match="estimated on"):
            estimate.subtract(other)

    def test_a_window_wider_than_the_scan_is_refused(self) -> None:
        measured = _synthetic_measurement(seed=None)
        with pytest.raises(ValueError, match="half the measured range"):
            estimate_background(measured, half_window_deg=90.0)

    def test_describe_names_the_method_its_settings_and_its_limits(self) -> None:
        measured = _synthetic_measurement(seed=None)
        text = estimate_background(measured, half_window_deg=2.0).describe()
        assert "snip" in text
        assert "half_window_deg" in text
        assert "doi:10.1016/0168-583X(88)90063-8" in text
        assert "modelling choice" in text

    def test_a_negative_background_cannot_be_constructed(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            PowderBackground(
                two_theta_deg=np.array([10.0, 11.0]),
                background=np.array([-1.0, 1.0]),
                observed_intensity=np.array([5.0, 5.0]),
                method="snip",
                parameters={},
                source_name="x",
            )


class TestInstrumentBroadening:
    def test_calibration_recovers_the_coefficients_it_was_built_from(self) -> None:
        truth = InstrumentBroadening(caglioti_u=0.004, caglioti_v=-0.002, caglioti_w=0.006)
        angles = np.array([21.4, 30.4, 37.5, 43.6, 53.9, 67.6, 79.3, 95.2, 115.4])
        fitted = calibrate_instrument_broadening(angles, truth.gaussian_fwhm_deg(angles))
        assert np.isclose(fitted.caglioti_u, 0.004, atol=1e-9)
        assert np.isclose(fitted.caglioti_v, -0.002, atol=1e-9)
        assert np.isclose(fitted.caglioti_w, 0.006, atol=1e-9)

    def test_calibration_needs_three_peaks(self) -> None:
        with pytest.raises(ValueError, match="at least three"):
            calibrate_instrument_broadening(np.array([30.0, 60.0]), np.array([0.1, 0.12]))

    def test_the_pseudo_voigt_width_lies_between_its_two_components(self) -> None:
        instrument = InstrumentBroadening.laboratory_bragg_brentano()
        angles = np.array([20.0, 60.0, 120.0])
        gaussian = instrument.gaussian_fwhm_deg(angles)
        lorentzian = instrument.lorentzian_fwhm_deg(angles)
        combined = instrument.fwhm_deg(angles)
        assert np.all(combined >= np.maximum(gaussian, lorentzian) - 1e-12)
        assert np.all(combined <= gaussian + lorentzian + 1e-12)
        assert np.all((instrument.eta(angles) >= 0.0) & (instrument.eta(angles) <= 1.0))

    def test_a_pure_gaussian_instrument_has_zero_lorentzian_fraction(self) -> None:
        instrument = InstrumentBroadening.ideal(0.1)
        assert np.allclose(instrument.eta(np.array([30.0, 90.0])), 0.0)
        assert np.allclose(instrument.fwhm_deg(np.array([30.0, 90.0])), 0.1)

    def test_a_peak_narrower_than_the_instrument_is_refused(self) -> None:
        instrument = InstrumentBroadening.laboratory_bragg_brentano()
        with pytest.raises(ValueError, match="cannot be sharper than the"):
            deconvolve_instrument_width(np.array([0.01]), instrument, np.array([30.0]))

    def test_williamson_hall_recovers_a_known_size_and_strain(self) -> None:
        instrument = InstrumentBroadening.laboratory_bragg_brentano()
        angles = np.array([21.4, 30.4, 37.5, 43.6, 53.9, 67.6, 79.3, 95.2, 115.4])
        theta = np.deg2rad(0.5 * angles)
        size_nm, strain = 25.0, 0.002
        sample = np.rad2deg(
            0.9 * CU_KA.wavelength_angstrom / (10.0 * size_nm) / np.cos(theta)
            + 4.0 * strain * np.tan(theta)
        )
        observed = np.sqrt(sample**2 + instrument.fwhm_deg(angles) ** 2)
        recovered = deconvolve_instrument_width(observed, instrument, angles, mode="gaussian")
        analysis = williamson_hall(
            angles, recovered, wavelength_angstrom=CU_KA.wavelength_angstrom
        )
        assert np.isclose(analysis.crystallite_size_nm, size_nm, rtol=1e-6)
        assert np.isclose(analysis.microstrain, strain, rtol=1e-6)
        assert analysis.r_squared > 0.999

    def test_scherrer_is_a_lower_bound_on_size_when_strain_is_present(self) -> None:
        angles = np.array([30.0, 60.0, 90.0, 120.0])
        theta = np.deg2rad(0.5 * angles)
        widths = np.rad2deg(
            0.9 * CU_KA.wavelength_angstrom / 300.0 / np.cos(theta)
            + 4.0 * 0.003 * np.tan(theta)
        )
        sizes = scherrer_size_nm(widths, angles, wavelength_angstrom=CU_KA.wavelength_angstrom)
        analysis = williamson_hall(
            angles, widths, wavelength_angstrom=CU_KA.wavelength_angstrom
        )
        assert np.all(sizes < analysis.crystallite_size_nm)

    def test_describe_states_the_convention_and_the_fit_quality(self) -> None:
        angles = np.array([30.0, 60.0, 90.0, 120.0])
        theta = np.deg2rad(0.5 * angles)
        widths = np.rad2deg(0.9 * CU_KA.wavelength_angstrom / 300.0 / np.cos(theta))
        text = williamson_hall(
            angles, widths, wavelength_angstrom=CU_KA.wavelength_angstrom
        ).describe()
        assert "beta cos(theta) = K lambda / D + 4 epsilon sin(theta)" in text
        assert "doi:10.1016/0001-6160(53)90006-6" in text
        assert "anisotropic" in text


class TestRietveldRefinement:
    def test_a_known_cell_zero_and_width_are_recovered(self) -> None:
        measured = _synthetic_measurement(
            lattice_scale=1.003, zero_shift_deg=0.05, fwhm_deg=0.14
        )
        result = refine_rietveld(
            measured,
            _silicon(),
            radiation=CU_KA,
            instrument=InstrumentBroadening.ideal(0.10),
            background_degree=4,
        )
        assert result.converged
        assert np.isclose(result.parameter("lattice_scale").value, 1.003, atol=2e-5)
        assert np.isclose(result.parameter("zero_shift_deg").value, 0.05, atol=5e-3)
        # W is the squared FWHM for a constant-width model.
        assert np.isclose(result.parameter("caglioti_w").value, 0.14**2, rtol=0.05)
        # Poisson noise and a correct model give a goodness of fit near one.
        assert 0.8 < result.goodness_of_fit < 1.5

    def test_the_refined_phase_carries_the_refined_cell(self) -> None:
        measured = _synthetic_measurement(lattice_scale=1.002)
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        expected = _silicon().lattice.a * result.parameter("lattice_scale").value
        assert np.isclose(result.phase.lattice.a, expected, rtol=1e-12)
        # The phase and its unit cell must not disagree about the metric.
        assert result.phase.unit_cell is not None
        assert np.isclose(result.phase.unit_cell.lattice.a, expected, rtol=1e-12)

    def test_a_noise_free_pattern_fits_essentially_exactly(self) -> None:
        measured = _synthetic_measurement(seed=None, background_level=100.0)
        result = refine_rietveld(
            measured,
            _silicon(),
            radiation=CU_KA,
            instrument=InstrumentBroadening.ideal(0.10),
            background_degree=6,
        )
        assert result.weighted_profile_r_factor < 0.02
        assert result.bragg_r_factor < 0.02

    def test_r_factors_and_the_residual_are_self_consistent(self) -> None:
        measured = _synthetic_measurement()
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        assert np.allclose(
            result.residual_intensity,
            result.observed_intensity - result.calculated_intensity,
        )
        expected_rp = float(
            np.sum(np.abs(result.residual_intensity)) / np.sum(result.observed_intensity)
        )
        assert np.isclose(result.profile_r_factor, expected_rp)
        assert np.isclose(
            result.goodness_of_fit,
            result.weighted_profile_r_factor / result.expected_r_factor,
        )
        assert np.isclose(result.reduced_chi_squared, result.goodness_of_fit**2)

    def test_texture_strength_is_recovered_when_an_axis_is_stated(self) -> None:
        from pytex.core.miller import MillerPlane
        from pytex.diffraction.preferred_orientation import MarchDollaseModel

        phase = _silicon()
        model = MarchDollaseModel(
            preferred_orientation=MillerPlane(
                indices=np.array([1, 1, 1], dtype=np.int64), phase=phase
            ),
            march_coefficient=0.75,
        )
        measured = _synthetic_measurement(seed=None, preferred_orientation=model)
        result = refine_rietveld(
            measured,
            phase,
            radiation=CU_KA,
            instrument=InstrumentBroadening.ideal(0.10),
            refine=("scale", "caglioti_w", "march_coefficient"),
            preferred_orientation_plane=(1, 1, 1),
        )
        assert np.isclose(result.parameter("march_coefficient").value, 0.75, atol=0.02)

    def test_refining_texture_without_an_axis_is_refused(self) -> None:
        measured = _synthetic_measurement(seed=None)
        with pytest.raises(ValueError, match="preferred_orientation_plane"):
            refine_rietveld(
                measured,
                _silicon(),
                radiation=CU_KA,
                refine=("scale", "march_coefficient"),
            )

    def test_an_unknown_parameter_name_is_refused(self) -> None:
        measured = _synthetic_measurement(seed=None)
        with pytest.raises(ValueError, match="Unknown refinement parameters"):
            refine_rietveld(measured, _silicon(), radiation=CU_KA, refine=("scale", "u_iso"))

    def test_a_phase_with_no_reflections_in_the_window_stops_the_refinement(self) -> None:
        measured = _synthetic_measurement(seed=None)
        with pytest.raises(ValueError, match="do not correspond"):
            refine_rietveld(
                measured,
                _silicon(),
                radiation=CU_KA,
                two_theta_range_deg=(20.0, 24.0),
            )

    def test_scale_is_always_refined_even_if_it_is_not_asked_for(self) -> None:
        measured = _synthetic_measurement(seed=None)
        result = refine_rietveld(
            measured,
            _silicon(),
            radiation=CU_KA,
            instrument=InstrumentBroadening.ideal(0.10),
            refine=("caglioti_w",),
        )
        assert result.parameter("scale").refined

    def test_a_fixed_parameter_carries_no_uncertainty(self) -> None:
        measured = _synthetic_measurement(seed=None)
        result = refine_rietveld(
            measured,
            _silicon(),
            radiation=CU_KA,
            instrument=InstrumentBroadening.ideal(0.10),
            refine=("scale", "caglioti_w"),
        )
        fixed = result.parameter("lattice_scale")
        assert not fixed.refined
        assert fixed.standard_uncertainty is None
        assert result.parameter("caglioti_w").standard_uncertainty is not None

    def test_counting_statistics_choose_poisson_weights(self) -> None:
        measured = _synthetic_measurement(seed=None)
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        assert result.weight_model == "poisson"

    def test_stated_uncertainties_choose_inverse_variance_weights(self) -> None:
        base = _synthetic_measurement(seed=None)
        measured = MeasuredPowderPattern(
            name=base.name,
            two_theta_deg=base.two_theta_deg,
            intensity=base.intensity,
            standard_uncertainty=np.sqrt(np.clip(base.intensity, 1.0, None)),
            radiation=CU_KA,
            synthetic=True,
        )
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        assert result.weight_model == "inverse_variance"

    def test_uncertainties_are_formatted_in_crystallographic_notation(self) -> None:
        measured = _synthetic_measurement()
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        formatted = result.parameter("lattice_scale").format()
        assert formatted.endswith(")") and "(" in formatted

    def test_as_pattern_returns_the_bragg_profile_without_the_background(self) -> None:
        measured = _synthetic_measurement(seed=None)
        result = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        )
        pattern = result.as_pattern()
        assert pattern.phase is result.phase
        assert np.isclose(float(np.max(pattern.intensity_grid)), 1.0)
        assert pattern.reflections

    def test_describe_reports_the_indices_the_scope_and_the_verdict(self) -> None:
        measured = _synthetic_measurement()
        text = refine_rietveld(
            measured, _silicon(), radiation=CU_KA, instrument=InstrumentBroadening.ideal(0.10)
        ).describe()
        assert "R_wp" in text and "R_exp" in text and "R_Bragg" in text
        assert "doi:10.1107/S0021889869006558" in text
        assert "does not refine atomic coordinates" in text
        assert "Poisson counting weights" in text
