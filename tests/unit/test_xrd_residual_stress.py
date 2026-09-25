"""Tests for residual-stress determination by the sin^2(psi) method.

Every expected value here has independent provenance: a closed-form identity
(the strain projection, the isotropic constants, the Reuss constants of a cubic
crystal, Kroener's cubic equation, the transversely isotropic constants of a
hexagonal basal reflection), or a stress that a noise-free synthetic data set
was generated from and must be recovered exactly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_residual_stress import (
    RESIDUAL_STRESS_SCHEMA,
    STRESS_COMPONENTS,
    DiffractionElasticConstants,
    Sin2PsiMeasurement,
    StressPeak,
    StressScan,
    _window_centroid_weights,
    determine_residual_stress,
    fit_sin2psi_lines,
    kroener_shear_modulus_cubic,
    locate_stress_peak,
    lpa_factor,
    measurement_direction,
    parse_stress_peak_positions,
    parse_stress_scans,
    residual_stress_pipeline,
    simulate_sin2psi_measurement,
    single_crystal_stiffness,
    strain_design_matrix,
)
from pytex.properties.tensors import StiffnessTensor

WAVELENGTH = RadiationSpec.cr_ka().wavelength_angstrom
D0_FE_211 = 2.8665 / math.sqrt(6.0)
PHIS = (0.0, 45.0, 90.0)
PSIS = (-45.0, -35.0, -25.0, -15.0, 0.0, 15.0, 25.0, 35.0, 45.0)
BIAXIAL = {"sigma_11": -350.0, "sigma_22": -150.0, "sigma_12": 60.0}


def _fe_decs(model: str = "kroener", **kwargs: float) -> DiffractionElasticConstants:
    return DiffractionElasticConstants.from_single_crystal(
        single_crystal_stiffness("fe_bcc"), [2, 1, 1], model=model, hkl=(2, 1, 1), **kwargs
    )


def _exact_peaks(
    stress: np.ndarray,
    dec: DiffractionElasticConstants,
    *,
    d0: float = D0_FE_211,
    phis: tuple[float, ...] = PHIS,
    psis: tuple[float, ...] = PSIS,
    u_two_theta: float = 1e-3,
) -> tuple[StressPeak, ...]:
    """Noise-free peak positions of a known stress, from the defining equation."""

    peaks = []
    for phi in phis:
        for psi in psis:
            m = measurement_direction(phi, psi)
            strain = 1e-6 * (
                dec.half_s2_per_tpa * float(m @ stress @ m)
                + dec.s1_per_tpa * float(np.trace(stress))
            )
            d = d0 * (1.0 + strain)
            two_theta = math.degrees(2.0 * math.asin(WAVELENGTH / (2.0 * d)))
            peaks.append(
                StressPeak(
                    phi_deg=phi,
                    psi_deg=psi,
                    two_theta_deg=two_theta,
                    two_theta_uncertainty_deg=u_two_theta,
                    method="given",
                )
            )
    return tuple(peaks)


def _tensor(values: dict[str, float]) -> np.ndarray:
    tensor = np.zeros((3, 3))
    index = {
        "sigma_11": (0, 0),
        "sigma_22": (1, 1),
        "sigma_12": (0, 1),
        "sigma_13": (0, 2),
        "sigma_23": (1, 2),
        "sigma_33": (2, 2),
    }
    for name, value in values.items():
        i, j = index[name]
        tensor[i, j] = tensor[j, i] = value
    return tensor


# ---------------------------------------------------------------------------
# Geometry and the fundamental equation
# ---------------------------------------------------------------------------


def test_measurement_direction_convention() -> None:
    assert np.allclose(measurement_direction(0.0, 0.0), [0.0, 0.0, 1.0])
    assert np.allclose(measurement_direction(0.0, 90.0), [1.0, 0.0, 0.0])
    assert np.allclose(measurement_direction(90.0, 90.0), [0.0, 1.0, 0.0])
    directions = measurement_direction(np.array(PHIS)[:, None], np.array(PSIS)[None, :])
    assert directions.shape == (3, len(PSIS), 3)
    assert np.allclose(np.linalg.norm(directions, axis=-1), 1.0)


def test_design_matrix_reproduces_the_tensor_contraction() -> None:
    generator = np.random.default_rng(7)
    raw = generator.normal(size=(3, 3)) * 200.0
    stress = 0.5 * (raw + raw.T)
    vector = np.array(
        [stress[0, 0], stress[1, 1], stress[0, 1], stress[0, 2], stress[1, 2], stress[2, 2]]
    )
    phi = generator.uniform(0.0, 360.0, 20)
    psi = generator.uniform(-70.0, 70.0, 20)
    design = strain_design_matrix(phi, psi, s1_per_tpa=-1.25, half_s2_per_tpa=5.8)
    m = measurement_direction(phi, psi)
    expected = 1e-6 * (5.8 * np.einsum("ni,ij,nj->n", m, stress, m) - 1.25 * np.trace(stress))
    assert np.allclose(design @ vector, expected, rtol=1e-12, atol=1e-18)


def test_design_matrix_column_selection_and_unknown_component() -> None:
    full = strain_design_matrix([0.0, 30.0], [10.0, 20.0], s1_per_tpa=-1.0, half_s2_per_tpa=5.0)
    part = strain_design_matrix(
        [0.0, 30.0],
        [10.0, 20.0],
        s1_per_tpa=-1.0,
        half_s2_per_tpa=5.0,
        components=("sigma_22", "sigma_11"),
    )
    assert np.allclose(part, full[:, [1, 0]])
    with pytest.raises(ValueError, match="Unknown stress component"):
        strain_design_matrix(
            [0.0], [0.0], s1_per_tpa=-1.0, half_s2_per_tpa=5.0, components=("tau",)
        )


# ---------------------------------------------------------------------------
# Diffraction elastic constants
# ---------------------------------------------------------------------------


def test_isotropic_constants_from_youngs_modulus_and_poisson_ratio() -> None:
    dec = DiffractionElasticConstants.isotropic(210.0, 0.28)
    assert dec.s1_per_tpa == pytest.approx(-0.28 / 210.0 * 1000.0)
    assert dec.half_s2_per_tpa == pytest.approx(1.28 / 210.0 * 1000.0)
    assert dec.youngs_modulus_gpa == pytest.approx(210.0)
    assert dec.poisson_ratio == pytest.approx(0.28)


@pytest.mark.parametrize("hkl", [(2, 0, 0), (2, 1, 1), (2, 2, 2), (3, 1, 0)])
def test_reuss_constants_match_the_cubic_closed_form(hkl: tuple[int, int, int]) -> None:
    # S1 = S12 + S0 Gamma and 1/2 S2 = S11 - S12 - 3 S0 Gamma, with
    # S0 = S11 - S12 - S44/2 and Gamma = (h2k2 + k2l2 + l2h2)/(h2 + k2 + l2)^2
    # (Voigt compliance, S44 carrying its factor of four).
    stiffness = single_crystal_stiffness("fe_bcc")
    compliance = np.linalg.inv(np.asarray(stiffness.voigt_matrix())) * 1000.0
    s11, s12, s44 = compliance[0, 0], compliance[0, 1], compliance[3, 3]
    h, k, l = (float(value) for value in hkl)  # noqa: E741
    gamma = (h * h * k * k + k * k * l * l + l * l * h * h) / (h * h + k * k + l * l) ** 2
    s0 = s11 - s12 - 0.5 * s44
    dec = DiffractionElasticConstants.from_single_crystal(stiffness, hkl, model="reuss")
    assert dec.s1_per_tpa == pytest.approx(s12 + s0 * gamma, rel=1e-10)
    assert dec.half_s2_per_tpa == pytest.approx(s11 - s12 - 3.0 * s0 * gamma, rel=1e-10)


def test_voigt_constants_are_reflection_independent_and_hill_is_the_mean() -> None:
    stiffness = single_crystal_stiffness("fe_bcc")
    voigt = [
        DiffractionElasticConstants.from_single_crystal(stiffness, hkl, model="voigt")
        for hkl in ((2, 0, 0), (2, 1, 1), (2, 2, 2))
    ]
    assert all(dec.half_s2_per_tpa == pytest.approx(voigt[0].half_s2_per_tpa) for dec in voigt)
    # Voigt shear modulus of a cubic crystal: (C11 - C12 + 3 C44) / 5.
    shear = (231.4 - 134.7 + 3.0 * 116.4) / 5.0
    assert voigt[0].half_s2_per_tpa == pytest.approx(1000.0 / (2.0 * shear), rel=1e-10)
    reuss = DiffractionElasticConstants.from_single_crystal(stiffness, (2, 1, 1), model="reuss")
    hill = DiffractionElasticConstants.from_single_crystal(stiffness, (2, 1, 1), model="hill")
    assert hill.half_s2_per_tpa == pytest.approx(
        0.5 * (reuss.half_s2_per_tpa + voigt[1].half_s2_per_tpa)
    )
    assert hill.s1_per_tpa == pytest.approx(0.5 * (reuss.s1_per_tpa + voigt[1].s1_per_tpa))


def test_kroener_iteration_reaches_the_root_of_kroeners_cubic() -> None:
    from pytex.diffraction.xrd_residual_stress import _kroener_grain_compliance, _to_mandel

    for c11, c12, c44 in ((231.4, 134.7, 116.4), (168.4, 121.4, 75.4), (107.3, 60.9, 28.3)):
        stiffness = StiffnessTensor.cubic(c11, c12, c44)
        _, bulk, shear = _kroener_grain_compliance(_to_mandel(np.asarray(stiffness.tensor)))
        assert shear == pytest.approx(kroener_shear_modulus_cubic(c11, c12, c44), rel=1e-10)
        # The bulk modulus of a cubic aggregate is exact, (C11 + 2 C12)/3.
        assert bulk == pytest.approx((c11 + 2.0 * c12) / 3.0, rel=1e-10)


@pytest.mark.parametrize("hkl", [(2, 0, 0), (2, 1, 1), (2, 2, 2), (3, 1, 0)])
def test_kroener_constants_lie_between_the_reuss_and_voigt_bounds(
    hkl: tuple[int, int, int],
) -> None:
    stiffness = single_crystal_stiffness("fe_bcc")
    values = {
        model: DiffractionElasticConstants.from_single_crystal(
            stiffness, hkl, model=model
        ).half_s2_per_tpa
        for model in ("reuss", "voigt", "kroener")
    }
    low, high = sorted((values["reuss"], values["voigt"]))
    assert low - 1e-12 <= values["kroener"] <= high + 1e-12


def test_every_model_agrees_for_an_isotropic_crystal() -> None:
    stiffness = StiffnessTensor.isotropic(youngs_modulus=200.0, poisson_ratio=0.3)
    reference = DiffractionElasticConstants.isotropic(200.0, 0.3)
    for model in ("reuss", "voigt", "hill", "kroener"):
        dec = DiffractionElasticConstants.from_single_crystal(stiffness, (1, 2, 3), model=model)
        assert dec.s1_per_tpa == pytest.approx(reference.s1_per_tpa, rel=1e-9)
        assert dec.half_s2_per_tpa == pytest.approx(reference.half_s2_per_tpa, rel=1e-9)


def test_hexagonal_basal_reflection_reuss_constants() -> None:
    # With the plane normal along c, the Reuss grain strain along c is S13 per
    # unit stress across c and S33 per unit stress along c, for any rotation
    # about c: S1 = S13 and S1 + 1/2 S2 = S33.
    stiffness = single_crystal_stiffness("ti_hcp")
    compliance = np.linalg.inv(np.asarray(stiffness.voigt_matrix())) * 1000.0
    dec = DiffractionElasticConstants.from_single_crystal(stiffness, [0, 0, 1], model="reuss")
    assert dec.s1_per_tpa == pytest.approx(compliance[0, 2], rel=1e-10)
    assert dec.s1_per_tpa + dec.half_s2_per_tpa == pytest.approx(compliance[2, 2], rel=1e-10)


def test_for_reflection_resolves_the_normal_through_the_reciprocal_basis() -> None:
    from pytex.app.phases import builtin_phase

    phase = builtin_phase("fe_bcc").to_phase()
    stiffness = single_crystal_stiffness("fe_bcc")
    via_phase = DiffractionElasticConstants.for_reflection(
        phase, (2, 1, 1), stiffness, model="reuss"
    )
    direct = DiffractionElasticConstants.from_single_crystal(stiffness, [2, 1, 1], model="reuss")
    assert via_phase.half_s2_per_tpa == pytest.approx(direct.half_s2_per_tpa, rel=1e-12)
    assert via_phase.hkl == (2, 1, 1)


def test_constants_are_validated_at_construction() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        DiffractionElasticConstants(s1_per_tpa=-1.0, half_s2_per_tpa=0.0)
    with pytest.raises(ValueError, match="compressibility"):
        DiffractionElasticConstants(s1_per_tpa=-3.0, half_s2_per_tpa=5.0)
    with pytest.raises(ValueError, match="model"):
        DiffractionElasticConstants(s1_per_tpa=-1.0, half_s2_per_tpa=5.0, model="bogus")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Poisson"):
        DiffractionElasticConstants.isotropic(200.0, 0.6)
    with pytest.raises(KeyError, match="known"):
        single_crystal_stiffness("unobtainium")


def test_constants_describe_and_serialize() -> None:
    dec = _fe_decs(relative_standard_uncertainty=0.05)
    text = dec.describe()
    assert "(211)" in text and "Kroener" in text and "5.0 %" in text
    payload = dec.to_json()
    assert payload["model"] == "kroener" and payload["hkl"] == [2, 1, 1]
    assert dec.with_uncertainty(0.0).relative_standard_uncertainty == 0.0


# ---------------------------------------------------------------------------
# The tensor fit on exact data
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("state", "stress"),
    [
        ("biaxial", BIAXIAL),
        ("biaxial_shear", {**BIAXIAL, "sigma_13": 25.0, "sigma_23": -15.0}),
        (
            "triaxial",
            {**BIAXIAL, "sigma_13": 25.0, "sigma_23": -15.0, "sigma_33": 40.0},
        ),
    ],
)
def test_exact_data_recover_the_generating_stress(state: str, stress: dict[str, float]) -> None:
    dec = _fe_decs()
    peaks = _exact_peaks(_tensor(stress), dec)
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        stress_state=state,  # type: ignore[arg-type]
    )
    assert result.tensor is not None
    for name in STRESS_COMPONENTS:
        value, _ = result.tensor.component(name)
        assert value == pytest.approx(stress.get(name, 0.0), abs=1e-5)
    assert result.tensor.reduced_chi_squared < 1e-6
    assert np.allclose(result.tensor.tensor_mpa, _tensor(stress), atol=1e-5)


def test_refining_d0_under_plane_stress_recovers_both() -> None:
    dec = _fe_decs()
    peaks = _exact_peaks(_tensor(BIAXIAL), dec)
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211 * 1.002,
        dec=dec,
        refine_d0=True,
    )
    assert result.d0_refined
    assert result.d0_angstrom == pytest.approx(D0_FE_211, rel=1e-10)
    assert result.tensor is not None
    assert np.allclose(result.tensor.values_mpa, [-350.0, -150.0, 60.0], atol=1e-4)
    with pytest.raises(ValueError, match="triaxial"):
        determine_residual_stress(
            peaks,
            wavelength_angstrom=WAVELENGTH,
            d0_angstrom=D0_FE_211,
            dec=dec,
            stress_state="triaxial",
            refine_d0=True,
        )


def test_per_azimuth_regressions_read_sigma_phi_and_tau_phi() -> None:
    dec = _fe_decs()
    stress = {**BIAXIAL, "sigma_13": 25.0, "sigma_23": -15.0}
    peaks = _exact_peaks(_tensor(stress), dec)
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        stress_state="biaxial_shear",
    )
    assert [line.phi_deg for line in result.regressions] == list(PHIS)
    for line in result.regressions:
        phi = math.radians(line.phi_deg)
        sigma_phi = (
            -350.0 * math.cos(phi) ** 2 + 60.0 * math.sin(2 * phi) - 150.0 * math.sin(phi) ** 2
        )
        tau_phi = 25.0 * math.cos(phi) - 15.0 * math.sin(phi)
        # The slope is d0 1/2 S2 sigma_phi exactly; the regression divides by
        # d0 1/2 S2, so the only departure is the second-order strain in d.
        assert line.sigma_phi_mpa == pytest.approx(sigma_phi, abs=0.5)
        assert line.has_splitting_term
        assert line.tau_phi_mpa == pytest.approx(tau_phi, abs=0.2)
        s2psi, a1, _, a2, _ = line.branch_averages()
        assert s2psi.size == 4
        # a2 = d0 1/2 S2 tau_phi sin|2 psi|, to first order in the strain.
        expected = (
            D0_FE_211
            * dec.half_s2_per_tpa
            * 1e-6
            * tau_phi
            * np.sin(2.0 * np.arcsin(np.sqrt(s2psi)))
        )
        assert np.allclose(a2, expected, rtol=2e-3, atol=1e-9)
        assert np.all(np.isfinite(a1))


def test_curvature_is_detected() -> None:
    psi = np.array(PSIS)
    sin2 = np.sin(np.deg2rad(psi)) ** 2
    d = 1.17 - 2e-3 * sin2 + 2e-3 * sin2**2
    lines = fit_sin2psi_lines(
        np.zeros_like(psi),
        psi,
        d,
        np.full_like(psi, 2e-6),
        d0_angstrom=1.17,
        dec=_fe_decs(),
    )
    assert abs(lines[0].curvature_t) > 3.0
    assert "curved" in lines[0].describe()


def test_single_azimuth_reports_sigma_phi_without_a_tensor() -> None:
    dec = _fe_decs()
    peaks = _exact_peaks(_tensor(BIAXIAL), dec, phis=(0.0,))
    result = determine_residual_stress(
        peaks, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211, dec=dec
    )
    assert result.tensor is None
    assert "three azimuths" in (result.tensor_unavailable_reason or "")
    assert result.regressions[0].sigma_phi_mpa == pytest.approx(-350.0, abs=0.5)
    assert "not determined" in result.describe()


def test_shear_from_one_sign_of_psi_is_possible_and_much_less_certain() -> None:
    # With d0 known, sin^2(psi) and sin(2 psi) are independent functions even
    # over psi > 0 alone, so the shear components are identifiable -- but they
    # lean on the curvature of one branch rather than on the difference of two,
    # and their uncertainty shows it.
    dec = _fe_decs()
    stress = _tensor({**BIAXIAL, "sigma_13": 25.0})
    one_sided = determine_residual_stress(
        _exact_peaks(stress, dec, psis=(0.0, 15.0, 25.0, 35.0, 45.0)),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        stress_state="biaxial_shear",
    )
    two_sided = determine_residual_stress(
        _exact_peaks(stress, dec, psis=(-45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0)),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        stress_state="biaxial_shear",
    )
    assert one_sided.tensor is not None and two_sided.tensor is not None
    assert one_sided.tensor.component("sigma_13")[0] == pytest.approx(25.0, abs=1e-4)
    index = one_sided.tensor.component_names.index("sigma_13")
    assert (
        one_sided.tensor.covariance_internal_mpa2[index, index]
        > 2.0 * two_sided.tensor.covariance_internal_mpa2[index, index]
    )
    # Only two tilts cannot carry a slope, a shear term and an intercept.
    starved = determine_residual_stress(
        _exact_peaks(stress, dec, psis=(0.0, 40.0)),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        stress_state="biaxial_shear",
    )
    assert starved.tensor is None and starved.tensor_unavailable_reason


# ---------------------------------------------------------------------------
# The uncertainty budget
# ---------------------------------------------------------------------------


def test_d0_budget_is_the_exact_sensitivity() -> None:
    dec = _fe_decs()
    peaks = _exact_peaks(_tensor(BIAXIAL), dec)
    u_d0 = 1e-4
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        d0_uncertainty_angstrom=u_d0,
    )
    shifted = determine_residual_stress(
        peaks, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211 + u_d0, dec=dec
    )
    assert result.tensor is not None and shifted.tensor is not None
    expected = np.abs(shifted.tensor.values_mpa - result.tensor.values_mpa)
    assert np.allclose(result.tensor.budget_mpa["d0"], expected, rtol=1e-3)
    # A d0 error is a near-uniform strain, and a uniform strain is what S1 tr(sigma)
    # produces, so it lands on the in-plane normal stresses almost equally.
    d0_part = result.tensor.budget_mpa["d0"]
    assert d0_part[0] == pytest.approx(d0_part[1], rel=0.05)


def test_elastic_constant_budget_and_monte_carlo_agree_with_linear_propagation() -> None:
    dec = _fe_decs(relative_standard_uncertainty=0.05)
    # Elastic constants alone: a Monte Carlo draw of S1 and 1/2 S2 is an
    # independent check on the linearized sensitivity.
    alone = determine_residual_stress(
        _exact_peaks(_tensor(BIAXIAL), dec, u_two_theta=1e-7),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        monte_carlo_draws=6000,
        seed=5,
    ).tensor
    assert alone is not None and alone.monte_carlo_uncertainty_mpa is not None
    assert np.allclose(
        alone.budget_mpa["elastic constants"], alone.monte_carlo_uncertainty_mpa, rtol=0.06
    )
    # And everything together.
    peaks = _exact_peaks(_tensor(BIAXIAL), dec, u_two_theta=5e-3)
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        d0_uncertainty_angstrom=5e-5,
        monte_carlo_draws=4000,
        seed=11,
    )
    fit = result.tensor
    assert fit is not None
    assert set(fit.budget_mpa) == {"statistical", "d0", "elastic constants"}
    combined = np.sqrt(sum(value**2 for value in fit.budget_mpa.values()))
    assert np.allclose(fit.combined_uncertainty_mpa, combined, rtol=1e-9)
    assert fit.monte_carlo_uncertainty_mpa is not None
    assert np.allclose(fit.monte_carlo_uncertainty_mpa, fit.combined_uncertainty_mpa, rtol=0.1)


def test_nominal_uncertainties_take_the_error_from_the_scatter() -> None:
    dec = _fe_decs()
    generator = np.random.default_rng(3)
    exact = _exact_peaks(_tensor(BIAXIAL), dec)
    noisy = [
        f"{peak.phi_deg} {peak.psi_deg} {peak.two_theta_deg + generator.normal(0.0, 0.004):.6f}"
        for peak in exact
    ]
    peaks = parse_stress_peak_positions("phi psi 2theta\n" + "\n".join(noisy))
    assert all(peak.uncertainty_is_nominal for peak in peaks)
    result = determine_residual_stress(
        peaks, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211, dec=dec
    )
    stated = parse_stress_peak_positions("\n".join(f"{line} 0.004" for line in noisy))
    assert not any(peak.uncertainty_is_nominal for peak in stated)
    reference = determine_residual_stress(
        stated, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211, dec=dec
    )
    assert result.tensor is not None and reference.tensor is not None
    # The scatter is what it is whatever placeholder the positions carried.
    ratio = result.tensor.budget_mpa["statistical"] / reference.tensor.budget_mpa["statistical"]
    assert np.allclose(ratio, ratio[0])
    assert 0.5 < ratio[0] < 2.0


def test_principal_stresses_and_von_mises() -> None:
    dec = _fe_decs()
    result = determine_residual_stress(
        _exact_peaks(_tensor(BIAXIAL), dec),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
    )
    fit = result.tensor
    assert fit is not None
    principal = fit.in_plane_principal()
    eigen = np.linalg.eigvalsh(_tensor(BIAXIAL)[:2, :2])
    assert principal["sigma_I_mpa"] == pytest.approx(eigen[1], abs=1e-4)
    assert principal["sigma_II_mpa"] == pytest.approx(eigen[0], abs=1e-4)
    angle = math.radians(principal["angle_deg"])
    direction = np.array([math.cos(angle), math.sin(angle)])
    assert direction @ _tensor(BIAXIAL)[:2, :2] @ direction == pytest.approx(eigen[1], abs=1e-3)
    equivalent, _ = fit.von_mises()
    s11, s22, s12 = -350.0, -150.0, 60.0
    assert equivalent == pytest.approx(
        math.sqrt(s11**2 - s11 * s22 + s22**2 + 3.0 * s12**2), abs=1e-4
    )


def test_result_serializes_and_describes_itself() -> None:
    dec = _fe_decs(relative_standard_uncertainty=0.05)
    result = determine_residual_stress(
        _exact_peaks(_tensor(BIAXIAL), dec),
        wavelength_angstrom=WAVELENGTH,
        d0_angstrom=D0_FE_211,
        dec=dec,
        reflection_label="(211)",
        phase_name="ferrite",
    )
    payload = result.to_json()
    assert payload["schema"] == RESIDUAL_STRESS_SCHEMA
    assert payload["tensor"]["component_names"] == ["sigma_11", "sigma_22", "sigma_12"]
    assert len(payload["peaks"]) == len(PHIS) * len(PSIS)
    text = result.describe()
    for fragment in ("(211)", "ferrite", "S1 towards S2", "tensile stress is positive", "doi:"):
        assert fragment in text


# ---------------------------------------------------------------------------
# Peak location, corrections and synthetic data
# ---------------------------------------------------------------------------


def test_lpa_factor_forms() -> None:
    theta = math.radians(78.0)
    lp = (1.0 + math.cos(2 * theta) ** 2) / math.sin(theta) ** 2
    absorption = 1.0 - math.tan(math.radians(30.0)) / math.tan(theta)
    omega = lpa_factor([156.0], psi_deg=30.0, geometry="omega")
    chi = lpa_factor([156.0], psi_deg=30.0, geometry="chi")
    assert omega[0] == pytest.approx(lp * absorption, rel=1e-12)
    assert chi[0] == pytest.approx(lp, rel=1e-12)
    with pytest.raises(ValueError, match="below the"):
        lpa_factor([40.0], psi_deg=60.0, geometry="omega")


def test_window_centroid_weights_integrate_exactly() -> None:
    axis = np.linspace(0.0, 10.0, 101)
    weights, _ = _window_centroid_weights(axis, 4.37, 1.23)
    constant = np.ones_like(axis)
    assert weights[1] @ constant == pytest.approx(2.46, rel=1e-12)
    assert weights[0] @ constant / (weights[1] @ constant) == pytest.approx(4.37, rel=1e-12)
    linear = axis.copy()
    # Integral of x over [a, b] is (b^2 - a^2)/2, exact for a linear profile.
    lower, upper = 4.37 - 1.23, 4.37 + 1.23
    assert weights[1] @ linear == pytest.approx(0.5 * (upper**2 - lower**2), rel=1e-12)


def test_synthetic_measurement_is_what_it_says() -> None:
    dec = _fe_decs()
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=dec, seed=2
    )
    assert measurement.synthetic
    assert len(measurement) == 3 * 7
    assert measurement.true_stress_mpa is not None
    assert np.allclose(measurement.true_stress_mpa, _tensor(BIAXIAL))
    assert "3 azimuth(s)" in measurement.describe()


@pytest.mark.parametrize(
    ("method", "tolerance_mpa"), [("pseudo_voigt", 8.0), ("parabola", 25.0), ("centroid", 25.0)]
)
def test_pipeline_recovers_a_known_stress_from_noisy_scans(
    method: str, tolerance_mpa: float
) -> None:
    dec = _fe_decs()
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=dec, seed=5
    )
    result = residual_stress_pipeline(
        measurement,
        d0_angstrom=D0_FE_211,
        dec=dec,
        window_deg=8.0,
        peak_method=method,  # type: ignore[arg-type]
    )
    assert result.tensor is not None
    assert np.allclose(result.tensor.values_mpa, [-350.0, -150.0, 60.0], atol=tolerance_mpa)
    if method == "pseudo_voigt":
        assert 0.3 < result.tensor.reduced_chi_squared < 3.0
        assert all(peak.peak_fit is not None for peak in result.peaks)
    else:
        assert all(peak.kalpha2_stripped for peak in result.peaks)


def test_an_unstripped_doublet_biases_the_parabola() -> None:
    # The reason the stripping exists: at Cr K-alpha and 2theta ~ 156 degrees
    # the doublet is ~0.9 degrees apart and blended into a profile whose width
    # grows with tilt, so the parabola's bias changes from tilt to tilt.
    dec = _fe_decs()
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=dec, seed=5
    )
    result = residual_stress_pipeline(
        measurement,
        d0_angstrom=D0_FE_211,
        dec=dec,
        window_deg=8.0,
        peak_method="parabola",
        strip_doublet=False,
    )
    assert result.tensor is not None
    assert abs(result.tensor.values_mpa[0] + 350.0) > 50.0


def test_scan_table_round_trip() -> None:
    dec = _fe_decs()
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211,
        stress_mpa=BIAXIAL,
        dec=dec,
        seed=4,
        phi_deg=(0.0, 60.0, 120.0),
        psi_deg=(0.0, 20.0, 30.0, 40.0),
    )
    lines = ["phi psi two_theta intensity"]
    for scan in reversed(measurement.scans):
        for angle, count in zip(scan.pattern.two_theta_deg, scan.pattern.intensity, strict=True):
            lines.append(f"{scan.phi_deg:g}, {scan.psi_deg:g}, {angle:.6f}, {count:g}")
    parsed = parse_stress_scans("\n".join(lines), radiation=measurement.radiation)
    assert len(parsed) == len(measurement)
    original = {(scan.phi_deg, scan.psi_deg): scan for scan in measurement.scans}
    for scan in parsed.scans:
        source = original[(scan.phi_deg, scan.psi_deg)]
        assert np.allclose(scan.pattern.intensity, source.pattern.intensity)
        located = locate_stress_peak(scan, expected_two_theta_deg=156.0, window_deg=8.0)
        reference = locate_stress_peak(source, expected_two_theta_deg=156.0, window_deg=8.0)
        assert located.two_theta_deg == pytest.approx(reference.two_theta_deg, abs=1e-5)
    with pytest.raises(ValueError, match="not numeric"):
        parse_stress_scans("0 0 150 10\n0 0 abc 11", radiation=measurement.radiation)


def test_containers_validate_their_invariants() -> None:
    pattern = (
        simulate_sin2psi_measurement(
            d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=_fe_decs(), seed=1
        )
        .scans[0]
        .pattern
    )
    with pytest.raises(ValueError, match="psi_deg"):
        StressScan(phi_deg=0.0, psi_deg=90.0, pattern=pattern)
    scans = tuple(StressScan(phi_deg=0.0, psi_deg=10.0, pattern=pattern) for _ in range(3))
    with pytest.raises(ValueError, match="two distinct tilts"):
        Sin2PsiMeasurement(name="x", scans=scans, radiation=RadiationSpec.cr_ka())
    with pytest.raises(ValueError, match="positive"):
        StressPeak(
            phi_deg=0.0,
            psi_deg=0.0,
            two_theta_deg=150.0,
            two_theta_uncertainty_deg=0.0,
            method="given",
        )
    with pytest.raises(ValueError, match="three peak"):
        determine_residual_stress(
            _exact_peaks(_tensor(BIAXIAL), _fe_decs())[:2],
            wavelength_angstrom=WAVELENGTH,
            d0_angstrom=D0_FE_211,
            dec=_fe_decs(),
        )


# ---------------------------------------------------------------------------
# Excluding measurements, and finding the ones to exclude
# ---------------------------------------------------------------------------


def _with_one_bad_point(shift_deg: float = 0.3) -> tuple[tuple[StressPeak, ...], int]:
    dec = _fe_decs()
    peaks = list(_exact_peaks(_tensor(BIAXIAL), dec))
    bad = next(
        index for index, peak in enumerate(peaks) if (peak.phi_deg, peak.psi_deg) == (45.0, 25.0)
    )
    original = peaks[bad]
    peaks[bad] = StressPeak(
        phi_deg=original.phi_deg,
        psi_deg=original.psi_deg,
        two_theta_deg=original.two_theta_deg + shift_deg,
        two_theta_uncertainty_deg=original.two_theta_uncertainty_deg,
        method="given",
    )
    return tuple(peaks), bad


def test_excluding_a_bad_point_restores_the_exact_stress() -> None:
    dec = _fe_decs()
    peaks, bad = _with_one_bad_point()
    kwargs = {"wavelength_angstrom": WAVELENGTH, "d0_angstrom": D0_FE_211, "dec": dec}
    spoiled = determine_residual_stress(peaks, **kwargs)
    assert spoiled.tensor is not None
    assert abs(spoiled.tensor.values_mpa[2] - 60.0) > 5.0
    cleaned = determine_residual_stress(peaks, excluded=(bad,), **kwargs)
    assert cleaned.tensor is not None
    assert np.allclose(cleaned.tensor.values_mpa, [-350.0, -150.0, 60.0], atol=1e-4)
    assert cleaned.excluded_indices == (bad,)
    assert not cleaned.included_mask[bad] and cleaned.included_mask.sum() == len(peaks) - 1
    # The excluded point stays in the result, with the prediction of the clean
    # fit and a residual that shows how far off it was.
    assert cleaned.d_angstrom.size == len(peaks)
    assert abs(cleaned.tensor.residual_strain[bad]) > 1e-4
    assert all(bad not in np.flatnonzero(line.psi_deg == 25.0) for line in cleaned.regressions)
    payload = cleaned.to_json()
    assert payload["excluded_indices"] == [bad] and payload["included"][bad] is False
    assert "excluded from every fit by the analyst" in cleaned.describe()


def test_deleted_residuals_point_at_the_bad_measurement_only() -> None:
    dec = _fe_decs()
    peaks, bad = _with_one_bad_point()
    noisy = [
        StressPeak(
            phi_deg=peak.phi_deg,
            psi_deg=peak.psi_deg,
            two_theta_deg=peak.two_theta_deg + noise,
            two_theta_uncertainty_deg=2e-3,
            method="given",
        )
        for peak, noise in zip(
            peaks, np.random.default_rng(4).normal(0.0, 2e-3, len(peaks)), strict=True
        )
    ]
    result = determine_residual_stress(
        noisy, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211, dec=dec
    )
    assert result.suggested_outliers() == (bad,)
    assert result.deleted_residuals is not None
    good = np.delete(result.deleted_residuals, bad)
    assert np.max(np.abs(good)) < 3.5 < abs(result.deleted_residuals[bad])
    # Once excluded it is no longer suggested, and its deleted residual, now
    # against the fit of every other point, still shows why it was excluded.
    cleaned = determine_residual_stress(
        noisy, wavelength_angstrom=WAVELENGTH, d0_angstrom=D0_FE_211, dec=dec, excluded=[bad]
    )
    assert cleaned.suggested_outliers() == ()
    assert cleaned.deleted_residuals is not None and abs(cleaned.deleted_residuals[bad]) > 10.0


def test_orientations_name_measurements_and_bad_names_are_refused() -> None:
    from pytex.diffraction.xrd_residual_stress import indices_of_orientations

    peaks = _exact_peaks(_tensor(BIAXIAL), _fe_decs())
    index = indices_of_orientations(peaks, [(405.0, -15.0)])
    assert len(index) == 1
    assert (peaks[index[0]].phi_deg, peaks[index[0]].psi_deg) == (45.0, -15.0)
    with pytest.raises(ValueError, match="No measurement"):
        indices_of_orientations(peaks, [(44.0, -15.0)])
    with pytest.raises(ValueError, match="outside"):
        determine_residual_stress(
            peaks,
            wavelength_angstrom=WAVELENGTH,
            d0_angstrom=D0_FE_211,
            dec=_fe_decs(),
            excluded=[len(peaks)],
        )
    with pytest.raises(ValueError, match="remain after the exclusions"):
        determine_residual_stress(
            peaks[:4],
            wavelength_angstrom=WAVELENGTH,
            d0_angstrom=D0_FE_211,
            dec=_fe_decs(),
            excluded=[0, 1],
        )


def test_the_simulator_plants_a_bad_measurement_where_asked() -> None:
    dec = _fe_decs()
    clean = simulate_sin2psi_measurement(d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=dec, seed=3)
    planted = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211,
        stress_mpa=BIAXIAL,
        dec=dec,
        seed=3,
        corrupted=[(45.0, 35.3)],
        corruption_shift_deg=0.3,
    )
    shifted = [
        (a.phi_deg, a.psi_deg)
        for a, b in zip(clean.scans, planted.scans, strict=True)
        if not np.allclose(a.pattern.two_theta_deg, b.pattern.two_theta_deg)
    ]
    assert shifted == [(45.0, 35.3)]


# ---------------------------------------------------------------------------
# A polarized (synchrotron) beam
# ---------------------------------------------------------------------------


def test_lpa_factor_follows_the_beam_polarization() -> None:
    theta = math.radians(40.0)
    unpolarized = lpa_factor([80.0], psi_deg=0.0, geometry="chi")[0]
    assert unpolarized == pytest.approx((1 + math.cos(2 * theta) ** 2) / math.sin(theta) ** 2)
    assert lpa_factor([80.0], psi_deg=0.0, geometry="chi", perpendicular_fraction=0.5)[0] == (
        pytest.approx(unpolarized)
    )
    vertical = lpa_factor([80.0], psi_deg=0.0, geometry="chi", perpendicular_fraction=1.0)[0]
    assert vertical == pytest.approx(2.0 / math.sin(theta) ** 2)
    horizontal = lpa_factor([80.0], psi_deg=0.0, geometry="chi", perpendicular_fraction=0.0)[0]
    assert horizontal == pytest.approx(2.0 * math.cos(2 * theta) ** 2 / math.sin(theta) ** 2)
    with pytest.raises(ValueError, match="perpendicular_fraction"):
        lpa_factor([80.0], psi_deg=0.0, perpendicular_fraction=1.5)


def test_a_synchrotron_measurement_recovers_the_stress_with_chi_tilting() -> None:
    dec = _fe_decs()
    beam = RadiationSpec.synchrotron(0.5)
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=D0_FE_211,
        stress_mpa=BIAXIAL,
        dec=dec,
        radiation=beam,
        geometry="chi",
        fwhm_deg=0.3,
        window_half_width_deg=1.5,
        step_deg=0.01,
        seed=6,
    )
    result = residual_stress_pipeline(
        measurement, d0_angstrom=D0_FE_211, dec=dec, window_deg=3.0, expected_fwhm_deg=0.3
    )
    assert result.tensor is not None
    assert all(
        peak.peak_fit is not None and not peak.peak_fit.doublet_modelled for peak in result.peaks
    )
    assert np.allclose(result.tensor.values_mpa, [-350.0, -150.0, 60.0], atol=25.0)
    with pytest.raises(ValueError, match="below the"):
        simulate_sin2psi_measurement(
            d0_angstrom=D0_FE_211, stress_mpa=BIAXIAL, dec=dec, radiation=beam, geometry="omega"
        )
