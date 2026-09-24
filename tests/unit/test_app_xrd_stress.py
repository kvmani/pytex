# ruff: noqa: RUF001
"""The residual-stress operation and its report.

The science is tested in ``test_xrd_residual_stress.py``. Here: the operation
recovers the stress its own demonstration was generated with; every quantity
the report derives is the rearrangement it claims to be; each input route and
each refusal behaves; the report carries every stage and figure it promises;
and the downloadable bundle holds the report, the figures and the result.
"""

from __future__ import annotations

import io
import math
import zipfile

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.export import result_to_bundle, result_to_markdown
from pytex.app.services.xrd_stress_report import tensor_sigma_phi
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants,
    determine_residual_stress,
    measurement_direction,
    parse_stress_peak_positions,
    simulate_sin2psi_measurement,
    single_crystal_stiffness,
)

_OPERATION = "xrd.residual_stress"
_TRUE = {"sigma_11": -350.0, "sigma_22": -150.0, "sigma_12": 60.0}
_STAGES = (
    "stress_tensor",
    "uncertainty_budget",
    "measurement",
    "peak_positions",
    "sin2psi_lines",
    "fit_residuals",
    "linearity",
    "theory",
    "elastic_constants",
    "algorithm",
)


def _run(**overrides: object) -> dict:
    request: dict[str, object] = {"phase": {"builtin": "fe_bcc"}, "reflection": [2, 1, 1]}
    request.update(overrides)
    return REGISTRY.call(_OPERATION, request)


@pytest.fixture(scope="module")
def default() -> dict:
    return _run()


def _stage(result: dict, key: str) -> dict:
    return next(stage for stage in result["stages"] if stage["key"] == key)


def _positions_text(*, phis: tuple[float, ...] = (0.0, 45.0, 90.0), with_u: bool = True) -> str:
    """Exact peak positions of the known stress, as a pasted table."""

    dec = DiffractionElasticConstants.for_reflection(
        _ferrite(), (2, 1, 1), single_crystal_stiffness("fe_bcc"), model="kroener"
    )
    stress = np.array([[-350.0, 60.0, 0.0], [60.0, -150.0, 0.0], [0.0, 0.0, 0.0]])
    d0 = 2.8665 / math.sqrt(6.0)
    wavelength = RadiationSpec.cr_ka().wavelength_angstrom
    lines = ["# phi psi two_theta u"]
    for phi in phis:
        for psi in (-45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0):
            m = measurement_direction(phi, psi)
            strain = 1e-6 * (
                dec.half_s2_per_tpa * float(m @ stress @ m)
                + dec.s1_per_tpa * float(np.trace(stress))
            )
            angle = math.degrees(2.0 * math.asin(wavelength / (2.0 * d0 * (1.0 + strain))))
            lines.append(f"{phi} {psi} {angle:.7f}" + (" 0.002" if with_u else ""))
    return "\n".join(lines)


def _ferrite():  # type: ignore[no-untyped-def]
    from pytex.app.phases import builtin_phase

    return builtin_phase("fe_bcc").to_phase()


# ---------------------------------------------------------------------------
# The answer
# ---------------------------------------------------------------------------


def test_the_demonstration_recovers_its_generating_stress(default: dict) -> None:
    tensor = default["data"]["tensor"]
    for name, value in _TRUE.items():
        entry = tensor[name]
        # Combined u is dominated by u(d0), which is systematic; the
        # statistical part is what the noise draws from, so compare with it.
        statistical = next(
            row["statistical_mpa"]
            for row in _stage(default, "stress_tensor")["table"]["rows"]
            if row["true_mpa"] == value
        )
        assert abs(entry["value_mpa"] - value) < max(4.0 * statistical, 3.0)
        assert entry["uncertainty_mpa"] > 0.0


def test_the_report_has_every_stage_and_figure(default: dict) -> None:
    assert [stage["key"] for stage in default["stages"]] == list(_STAGES)
    figures = {figure["key"] for stage in default["stages"] for figure in stage.get("figures", [])}
    assert figures == {
        "sigma_phi",
        "mohr_circle",
        "uncertainty_budget",
        "measurement_directions",
        "peak_shift",
        "peak_fits",
        "peak_quality",
        "d_vs_sin2psi",
        "strain_vs_sin2psi",
        "strain_residuals",
        "stress_correlation",
        "psi_splitting",
        "elastic_constants",
    }
    sections = {stage["key"]: stage["section"] for stage in default["stages"]}
    assert sections["stress_tensor"] == "result"
    assert sections["sin2psi_lines"] == "evidence"
    assert sections["linearity"] == "diagnostics"
    assert sections["theory"] == "method"
    for stage in default["stages"]:
        for figure in stage.get("figures", []):
            assert figure["svg"].startswith("<svg")
            assert figure["caption"] and figure["interpretation"]


def test_the_table_has_one_row_per_measurement(default: dict) -> None:
    rows = default["table"]["rows"]
    assert len(rows) == 27
    for row in rows:
        assert row["sin2psi"] == pytest.approx(math.sin(math.radians(row["psi_deg"])) ** 2)
        # u(d) = d cot(theta) u(theta), in the table's own units.
        theta = math.radians(0.5 * row["two_theta_deg"])
        expected = (
            row["d_angstrom"]
            / math.tan(theta)
            * math.radians(0.5e-3 * row["two_theta_uncertainty_mdeg"])
        )
        assert row["d_uncertainty_angstrom"] == pytest.approx(expected, rel=1e-9)


def test_normalized_residuals_reproduce_the_strain_fit_chi_squared(default: dict) -> None:
    rows = default["table"]["rows"]
    z = np.array([row["normalized_residual"] for row in rows])
    dof = int(
        next(
            metric["value"]
            for metric in _stage(default, "fit_residuals")["metrics"]
            if metric["label"] == "Degrees of freedom"
        )
    )
    assert float(np.sum(z**2) / dof) == pytest.approx(default["data"]["reduced_chi_squared"])


def test_the_tensor_curve_is_sigma_phi_of_the_components() -> None:
    peaks = parse_stress_peak_positions(_positions_text())
    dec = DiffractionElasticConstants.for_reflection(
        _ferrite(), (2, 1, 1), single_crystal_stiffness("fe_bcc"), model="kroener"
    )
    result = determine_residual_stress(
        peaks,
        wavelength_angstrom=RadiationSpec.cr_ka().wavelength_angstrom,
        d0_angstrom=2.8665 / math.sqrt(6.0),
        dec=dec,
    )
    assert result.tensor is not None
    phi = np.array([0.0, 30.0, 45.0, 90.0, 135.0])
    values, uncertainty = tensor_sigma_phi(result.tensor, phi)
    radians = np.deg2rad(phi)
    expected = (
        -350.0 * np.cos(radians) ** 2 + 60.0 * np.sin(2 * radians) - 150.0 * np.sin(radians) ** 2
    )
    assert np.allclose(values, expected, atol=1e-3)
    assert np.all(uncertainty >= 0.0)


def test_budget_and_monte_carlo_are_reported(default: dict) -> None:
    rows = _stage(default, "stress_tensor")["table"]["rows"]
    for row in rows:
        combined = math.sqrt(row["statistical_mpa"] ** 2 + row["d0_mpa"] ** 2 + row["dec_mpa"] ** 2)
        assert row["combined_mpa"] == pytest.approx(combined, rel=1e-9)
        assert row["monte_carlo_mpa"] == pytest.approx(row["combined_mpa"], rel=0.15)
    assert "Monte Carlo" in _stage(default, "uncertainty_budget")["summary"]


# ---------------------------------------------------------------------------
# Input routes
# ---------------------------------------------------------------------------


def test_pasted_peak_positions() -> None:
    result = _run(data_source="positions", measurement=_positions_text())
    tensor = result["data"]["tensor"]
    for name, value in _TRUE.items():
        assert tensor[name]["value_mpa"] == pytest.approx(value, abs=0.05)
    # No scans: nothing to draw per profile.
    figures = {f["key"] for s in result["stages"] for f in s.get("figures", [])}
    assert "peak_fits" not in figures and "d_vs_sin2psi" in figures
    assert "measurement" not in result["inputs"]
    assert result["inputs"]["measurement_lines"] == 22


def test_pasted_scans_match_the_library_pipeline() -> None:
    dec = DiffractionElasticConstants.for_reflection(
        _ferrite(), (2, 1, 1), single_crystal_stiffness("fe_bcc"), model="kroener"
    )
    measurement = simulate_sin2psi_measurement(
        d0_angstrom=2.8665 / math.sqrt(6.0),
        stress_mpa=_TRUE,
        dec=dec,
        phi_deg=(0.0, 60.0, 120.0),
        psi_deg=(0.0, 21.1, 30.0, 37.8, 45.0),
        seed=3,
    )
    text = "\n".join(
        f"{scan.phi_deg} {scan.psi_deg} {angle:.5f} {count:g}"
        for scan in measurement.scans
        for angle, count in zip(scan.pattern.two_theta_deg, scan.pattern.intensity, strict=True)
    )
    result = _run(data_source="scans", measurement=text)
    assert len(result["table"]["rows"]) == 15
    tensor = result["data"]["tensor"]
    for name, value in _TRUE.items():
        assert tensor[name]["value_mpa"] == pytest.approx(value, abs=25.0)


def test_one_azimuth_gives_sigma_phi_and_says_why_there_is_no_tensor() -> None:
    result = _run(data_source="positions", measurement=_positions_text(phis=(0.0,)))
    assert "tensor" not in result["data"]
    assert any("three azimuths" in warning for warning in result["warnings"])
    assert _stage(result, "stress_tensor")["status"] == "warning"
    assert "σφ = -350" in result["summary"].replace("−", "-")


def test_nominal_positions_are_flagged() -> None:
    result = _run(data_source="positions", measurement=_positions_text(with_u=False))
    assert "scatter alone" in _stage(result, "measurement")["summary"]


def test_refined_d0_and_the_elastic_constant_models() -> None:
    refined = _run(refine_d0=True, a0_angstrom=2.870)
    assert refined["data"]["library"]["d0_refined"] is True
    tabulated = 2.8665 / math.sqrt(6.0)
    assert refined["data"]["library"]["d0_angstrom"] == pytest.approx(tabulated, rel=2e-4)
    isotropic = _run(dec_model="isotropic", youngs_modulus_gpa=210.0, poisson_ratio=0.28)
    assert isotropic["data"]["library"]["dec"]["half_s2_per_tpa"] == pytest.approx(
        1.28 / 210.0 * 1000.0
    )
    user = _run(dec_model="user", s1_per_tpa=-1.2, half_s2_per_tpa=5.6)
    assert user["data"]["library"]["dec"]["model"] == "user"
    custom = _run(stiffness_source="custom", c11_gpa=231.4, c12_gpa=134.7, c44_gpa=116.4)
    assert custom["data"]["library"]["dec"]["half_s2_per_tpa"] == pytest.approx(
        _run()["data"]["library"]["dec"]["half_s2_per_tpa"]
    )


def test_a_wrong_d0_moves_the_normal_stresses_not_the_slopes(default: dict) -> None:
    shifted = _run(a0_angstrom=2.8665 * 1.0005)
    for before, after in zip(default["data"]["series"], shifted["data"]["series"], strict=True):
        # A 5e-4 change of d0 rescales a slope-derived stress by 5e-4; the rest
        # of the difference is the fit window following the new Bragg angle,
        # which is noise, and stays inside the statistical uncertainty.
        assert (
            abs(after["sigma_phi_mpa"] - before["sigma_phi_mpa"])
            < before["sigma_phi_uncertainty_mpa"]
        )
    moved = (
        shifted["data"]["tensor"]["sigma_11"]["value_mpa"]
        - default["data"]["tensor"]["sigma_11"]["value_mpa"]
    )
    assert abs(moved) > 50.0


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def test_refusals_name_the_field_to_fix() -> None:
    with pytest.raises(InvalidInputError) as caught:
        _run(data_source="scans", measurement="0 0 abc 1")
    assert caught.value.details["field"] == "measurement"
    with pytest.raises(InvalidInputError) as caught:
        _run(reflection=[1, 1, 0], radiation="cr_ka", phase={"builtin": "al_fcc"}, a0_angstrom=1.5)
    assert caught.value.details["field"] == "reflection"
    with pytest.raises(InvalidInputError) as caught:
        _run(phase={"builtin": "zr_hcp"}, reflection=[1, 0, 3])
    assert caught.value.details["field"] == "stiffness_source"
    with pytest.raises(InvalidInputError) as caught:
        _run(stress_state="triaxial", refine_d0=True)
    assert caught.value.details["field"] == "refine_d0"


# ---------------------------------------------------------------------------
# The downloadable report
# ---------------------------------------------------------------------------


def test_the_bundle_holds_the_report_every_figure_and_the_result(default: dict) -> None:
    archive = zipfile.ZipFile(io.BytesIO(result_to_bundle(default)))
    names = set(archive.namelist())
    figures = [f["key"] for s in default["stages"] for f in s.get("figures", [])]
    assert {"report.md", "report.html", "result.json"} <= names
    assert {f"figures/{key.replace('_', '-')}.svg" for key in figures} <= names | {
        f"figures/{key}.svg" for key in figures
    }
    report = archive.read("report.md").decode("utf-8")
    for heading in (
        "## Result and reliability",
        "The stress tensor",
        "Uncertainty budget",
        "Peak positions",
        "d against sin²ψ at each azimuth",
        "Theory: from peak shift to stress",
        "Diffraction elastic constants",
        "Algorithm",
        "## Sources",
    ):
        assert heading in report
    single = result_to_markdown(default).decode("utf-8")
    assert single.count("data:image/svg+xml;base64,") == len(figures)
    printable = archive.read("report.html").decode("utf-8")
    assert printable.count("data:image/svg+xml;base64,") == len(figures)
    for heading in ("The stress tensor", "Theory: from peak shift to stress", "Sources"):
        assert heading in printable
