# ruff: noqa: RUF001, RUF002
"""The lattice-parameter report: unchanged numbers, correct derived quantities, honest words.

The report was rewritten for readability, and the rewrite had one hard
constraint: it must not change a single validated number. The first class pins
the numbers the operation produced before the rewrite. The rest check that the
quantities the report *derives* - normalized residuals, the correction curve,
the fitted peak profiles - are the algebraic rearrangements they claim to be,
that each warning fires on the condition it names and not otherwise, that the
figures are exportable, and that the report keeps apart the things it promises
to keep apart: peak-fit and lattice-fit chi-squared, precision and accuracy, a
change from the database cell and an elastic strain.
"""

from __future__ import annotations

import dataclasses
import io
import re
import zipfile
from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.export import result_to_bundle, result_to_markdown
from pytex.app.phases import phase_from_request
from pytex.app.services.xrd import _RADIATION, _measured_from_request
from pytex.app.services.xrd_lattice_report import (
    lattice_warnings,
    normalized_residuals,
    systematic_correction_curve,
)
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_lattice_parameter import LatticeParameterResult
from pytex.diffraction.xrd_peaks import detect_and_fit_peaks, kalpha_doublet_parameters

_LATTICE = "xrd.lattice_parameters"


def _run(**overrides: object) -> dict:
    request: dict[str, object] = {"phase": {"builtin": "ni_fcc"}}
    request.update(overrides)
    return REGISTRY.call(_LATTICE, request)


@pytest.fixture(scope="module")
def cohen() -> dict:
    return _run()


@pytest.fixture(scope="module")
def le_bail() -> dict:
    return _run(method="le_bail")


# Values produced by the operation on its deterministic demonstration scan. A
# regression pin, not a reference value: the independent checks of these
# methods live in test_xrd_lattice_parameter.py and the worked examples. The
# numbers were first taken at commit 87b1c38, before the report was rewritten,
# and re-taken when the peak fit gained its closed-form Jacobian and
# Gauss-Newton polish - a change that moved the cell by two parts in 1e12,
# five orders of magnitude below the uncertainty the same fit reports.
#
# The tolerances differ by quantity because the quantities differ in how much
# of the floating-point world they can see. A cell parameter is a smooth
# function of six fitted peak centres, and those centres are now bit-identical
# between Windows and Linux, so it pins to a part in 1e9 with room to spare
# (the measured spread is a part in 1e12). The uncertainties and the
# chi-squared are not: a standard uncertainty is the square root of an element
# of an inverted Hessian, and the Le Bail refinement builds that Hessian from
# finite differences, which are accurate to about the cube root of the machine
# epsilon however the fit is driven. Those pin to a part in 1e7, against a
# measured Windows-Linux spread of about a part in 1e9.
_CELL_TOLERANCE = 1.0e-9
_DERIVED_TOLERANCE = 1.0e-7
_CELL_KEYS = frozenset({"a", "b", "c"})

_BASELINE = {
    "cohen_nelson_riley": (
        {},
        {
            "a": 3.533881136172797,
            "a_standard_uncertainty": 1.513066770421227e-05,
            "reduced_chi_squared": 3.599242749919796,
            "drift_coefficient": 0.000814800518152223,
            "drift_standard_uncertainty": 7.84395463103818e-06,
            "reflection_count": 6,
        },
    ),
    "cohen_no_correction": (
        {"extrapolation": "none"},
        {
            "a": 3.5326029143107553,
            "a_standard_uncertainty": 0.00040820827328433366,
            "reduced_chi_squared": 7770.24202037487,
            "reflection_count": 6,
        },
    ),
    "average": (
        {"method": "average"},
        {
            "a": 3.532339281503157,
            "a_standard_uncertainty": 0.0004635606795724213,
            "reduced_chi_squared": 8418.572488274502,
        },
    ),
    "le_bail": (
        {"method": "le_bail"},
        {
            "a": 3.533931746475637,
            "a_standard_uncertainty": 1.9702216161025456e-05,
            "reduced_chi_squared": 6.129212634128142,
            "weighted_profile_r": 0.4404597011600123,
        },
    ),
    "hexagonal": (
        {"phase": {"builtin": "ti_hcp"}},
        {
            "a": 2.959212856068566,
            "c": 4.698849481692253,
            "c_standard_uncertainty": 3.097333238893004e-05,
            "reduced_chi_squared": 4.89737090541353,
            "reflection_count": 20,
        },
    ),
}

# A per-reflection residual is a difference of two angles near 2 theta = 100
# that comes out near 1e-4 degrees, so it carries eight fewer significant
# figures than the angles it is formed from. It holds to a part in 1e10
# between Windows and Linux; it is pinned at a part in 1e7.
_BASELINE_RESIDUALS_MDEG = [
    -0.07633598695015076,
    -0.2742989105833083,
    0.24275835211917496,
    0.5217337916794396,
    -0.05815062388274163,
    -1.8941820684288422,
]

class TestNumbersAreUnchanged:
    @pytest.mark.parametrize("case", sorted(_BASELINE))
    def test_the_rewrite_changed_no_validated_number(self, case: str) -> None:
        overrides, expected = _BASELINE[case]
        data = _run(**overrides)["data"]
        for key, value in expected.items():
            tolerance = _CELL_TOLERANCE if key in _CELL_KEYS else _DERIVED_TOLERANCE
            assert data[key] == pytest.approx(value, rel=tolerance), key

    def test_the_per_reflection_residuals_are_unchanged(self, cohen: dict) -> None:
        residuals = [row["residual_mdeg"] for row in cohen["table"]["rows"]]
        assert residuals == pytest.approx(
            _BASELINE_RESIDUALS_MDEG, rel=_DERIVED_TOLERANCE, abs=1e-12
        )


class TestDerivedQuantities:
    def test_normalized_residuals_sum_to_the_lattice_fit_chi_squared(self, cohen: dict) -> None:
        """Σz²/ν is the reported reduced χ²: the residual and σ convert by one derivative."""

        rows = cohen["table"]["rows"]
        z = np.array([row["normalized_residual"] for row in rows])
        degrees = next(
            metric["value"]
            for metric in next(s for s in cohen["stages"] if s["key"] == "lattice_fit")["metrics"]
            if metric["label"] == "Degrees of freedom"
        )
        assert degrees == len(rows) - 2  # a*² and the drift coefficient D
        assert float(np.sum(z**2) / degrees) == pytest.approx(
            cohen["data"]["reduced_chi_squared"], rel=1e-9
        )

    @pytest.mark.parametrize("overrides", [{"extrapolation": "none"}, {"method": "average"}])
    def test_the_identity_holds_for_the_other_position_methods(
        self, overrides: dict[str, str]
    ) -> None:
        result = _run(**overrides)
        z = np.array([row["normalized_residual"] for row in result["table"]["rows"]])
        # One refined parameter in both: a*² for Cohen without D, the mean for the average.
        assert float(np.sum(z**2) / (z.size - 1)) == pytest.approx(
            result["data"]["reduced_chi_squared"], rel=1e-9
        )

    def test_each_normalized_residual_uses_that_peaks_own_sigma(self, cohen: dict) -> None:
        assignment = next(s for s in cohen["stages"] if s["key"] == "assignment")
        sigma_by_angle = {
            round(row["two_theta_observed_deg"], 9): row["sigma_mdeg"]
            for row in assignment["table"]["rows"]
        }
        for row in cohen["table"]["rows"]:
            sigma = sigma_by_angle[round(row["two_theta_observed_deg"], 9)]
            assert row["standard_uncertainty_mdeg"] == pytest.approx(sigma)
            assert row["normalized_residual"] == pytest.approx(row["residual_mdeg"] / sigma)

    def test_normalized_residuals_are_undefined_without_an_uncertainty(self) -> None:
        result = _synthetic(residual_mdeg=(0.5, -0.5, 1.0))
        values = normalized_residuals(result, np.array([0.001, 0.0, 0.002]))
        assert values[0] == pytest.approx(0.5)
        assert np.isnan(values[1])
        assert values[2] == pytest.approx(0.5)

    def test_the_correction_curve_passes_through_every_reflection(self) -> None:
        result = _synthetic(drift=8.0e-4, drift_sigma=1.0e-5)
        curve, band = systematic_correction_curve(result, result.two_theta_deg)
        assert curve == pytest.approx(result.systematic_shift_deg, rel=1e-12)
        assert band == pytest.approx(np.abs(curve) * 1.0e-5 / 8.0e-4, rel=1e-12)

    def test_no_correction_draws_a_flat_zero(self) -> None:
        result = _synthetic(extrapolation="none")
        curve, band = systematic_correction_curve(result, np.linspace(20.0, 140.0, 9))
        assert np.all(curve == 0.0)
        assert np.all(band == 0.0)

    def test_a_fitted_profile_reproduces_its_own_chi_squared(self) -> None:
        """PeakFit.evaluate is the model the fit minimized, doublet and background included."""

        radiation = _RADIATION["cu_ka_doublet"]()
        _, phase = phase_from_request({"builtin": "ni_fcc"})
        request = REGISTRY.get(_LATTICE).bind({"phase": {"builtin": "ni_fcc"}})
        measured, _ = _measured_from_request(request, phase, radiation)
        table = detect_and_fit_peaks(measured, radiation=radiation)
        doublet = kalpha_doublet_parameters(radiation)
        assert doublet is not None
        axis = np.asarray(measured.two_theta_deg)
        counts = np.asarray(measured.intensity)
        for peak in table:
            assert peak.doublet_modelled
            inside = (axis >= peak.window_deg[0]) & (axis <= peak.window_deg[1])
            model = peak.evaluate(axis[inside], doublet=doublet)
            weights = 1.0 / np.sqrt(np.maximum(counts[inside], 1.0))
            chi = float(np.sum(((model - counts[inside]) * weights) ** 2)) / (inside.sum() - 6)
            assert chi == pytest.approx(peak.reduced_chi_squared, rel=1e-6)
        with pytest.raises(ValueError, match="doublet"):
            table[0].evaluate(axis[:5])

    def test_a_single_line_radiation_has_no_doublet(self) -> None:
        assert kalpha_doublet_parameters(None) is None
        single = RadiationSpec.cu_ka()
        assert kalpha_doublet_parameters(single) is None


def _synthetic(
    *,
    residual_mdeg: tuple[float, ...] = (0.2, -0.1, 0.3),
    chi: float = 1.0,
    drift: float = 0.0,
    drift_sigma: float = 0.0,
    extrapolation: str = "nelson_riley",
    method: str = "cohen",
    correlation: float = 0.5,
    count: int | None = None,
) -> LatticeParameterResult:
    """A hand-made cubic result, so a warning can be provoked by exactly one condition."""

    residuals = np.asarray(residual_mdeg, dtype=float) / 1000.0
    size = residuals.size if count is None else count
    angles = np.linspace(44.0, 120.0, residuals.size)
    names = ("a*^2", "D") if extrapolation != "none" else ("a*^2",)
    matrix = (
        np.array([[1.0, correlation], [correlation, 1.0]]) if len(names) == 2 else np.array([[1.0]])
    )
    return LatticeParameterResult(
        method=method,  # type: ignore[arg-type]
        phase_name="test",
        crystal_system="cubic",
        a=3.52,
        b=3.52,
        c=3.52,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        a_standard_uncertainty=1e-5,
        b_standard_uncertainty=1e-5,
        c_standard_uncertainty=1e-5,
        free_parameter_names=("a*^2",),
        extrapolation=extrapolation,  # type: ignore[arg-type]
        drift_coefficient=drift,
        drift_standard_uncertainty=drift_sigma,
        reflection_count=size,
        reduced_chi_squared=chi,
        residual_two_theta_deg=residuals,
        miller_indices=tuple((1, 1, index) for index in range(residuals.size)),
        two_theta_deg=angles,
        reference_lattice=phase_from_request({"builtin": "ni_fcc"})[1].lattice,
        parameter_correlation=matrix if method == "cohen" else None,
        correlation_parameter_names=names if method == "cohen" else (),
    )


class TestWarnings:
    def _warnings(self, **kwargs: object) -> tuple[str, ...]:
        normalized = kwargs.pop("normalized", None)
        labels = kwargs.pop("labels", ())
        result = _synthetic(**kwargs)  # type: ignore[arg-type]
        return lattice_warnings(result, normalized=normalized, labels=labels)  # type: ignore[arg-type]

    def test_a_clean_fit_raises_nothing(self) -> None:
        clean = self._warnings(
            residual_mdeg=(0.1, -0.1, 0.1, -0.1, 0.1, -0.1), drift=8e-4, drift_sigma=1e-5
        )
        assert clean == ()

    def test_poor_fit(self) -> None:
        warnings = self._warnings(residual_mdeg=(0.1,) * 6, chi=9.0, drift=8e-4, drift_sigma=1e-5)
        assert any(item.startswith("Poor lattice fit: reduced χ² = 9.00") for item in warnings)
        assert any("3.0 times" in item for item in warnings)

    def test_overstated_uncertainties(self) -> None:
        warnings = self._warnings(residual_mdeg=(0.1,) * 6, chi=0.1, drift=8e-4, drift_sigma=1e-5)
        assert any("well below 1" in item for item in warnings)

    def test_exactly_determined(self) -> None:
        warnings = self._warnings(residual_mdeg=(0.0, 0.0), drift=8e-4, drift_sigma=1e-5)
        assert any("exactly determined" in item for item in warnings)

    def test_few_degrees_of_freedom(self) -> None:
        warnings = self._warnings(residual_mdeg=(0.1, 0.1, 0.1, 0.1), drift=8e-4, drift_sigma=1e-5)
        assert any("Only 2 degrees of freedom" in item for item in warnings)

    def test_outliers_are_named(self) -> None:
        warnings = self._warnings(
            residual_mdeg=(0.1,) * 6,
            drift=8e-4,
            drift_sigma=1e-5,
            normalized=np.array([0.1, -3.4, 0.2, 2.9, 0.0, 3.1]),
            labels=["(111)", "(200)", "(220)", "(311)", "(222)", "(400)"],
        )
        outlier = next(item for item in warnings if "beyond ±3σ" in item)
        assert "(200) (-3.4σ)" in outlier
        assert "(400) (+3.1σ)" in outlier
        assert "(311)" not in outlier

    def test_strong_correlation(self) -> None:
        warnings = self._warnings(
            residual_mdeg=(0.1,) * 6, drift=8e-4, drift_sigma=1e-5, correlation=-0.97
        )
        assert any("correlated at -0.970" in item for item in warnings)

    def test_insignificant_correction(self) -> None:
        warnings = self._warnings(residual_mdeg=(0.1,) * 6, drift=1e-5, drift_sigma=1e-5)
        assert any("not significant (1.0σ)" in item for item in warnings)

    def test_no_correction_and_average(self) -> None:
        none = self._warnings(residual_mdeg=(0.1,) * 6, extrapolation="none")
        assert any("No systematic correction was refined" in item for item in none)
        average = self._warnings(residual_mdeg=(0.1,) * 6, extrapolation="none", method="average")
        assert any("average method cannot remove" in item for item in average)

    def test_the_demonstration_scan_is_not_a_defect(self) -> None:
        result = _synthetic(residual_mdeg=(0.1,) * 6, drift=8e-4, drift_sigma=1e-5)
        assert lattice_warnings(result, generated=True) == ()

    def test_a_large_change_from_the_reference_cell(self) -> None:
        result = dataclasses.replace(
            _synthetic(residual_mdeg=(0.1,) * 6, drift=8e-4, drift_sigma=1e-5), a=3.60
        )
        warnings = lattice_warnings(result)
        assert any("far more than any elastic strain" in item for item in warnings)

    def test_the_operation_reports_its_poor_fit(self, cohen: dict) -> None:
        # The demonstration scan carries a constant zero error, which the
        # Nelson-Riley form absorbs only approximately: chi-squared 3.6.
        assert any(item.startswith("Poor lattice fit") for item in cohen["warnings"])
        assert cohen["data"]["warnings"] == cohen["warnings"]


class TestReportLayout:
    def test_the_headline_leads_with_the_cell_and_its_reliability(self, cohen: dict) -> None:
        labels = [metric["label"] for metric in cohen["highlights"]]
        assert labels[0] == "a"
        assert "±" in cohen["highlights"][0]["value"]
        for required in (
            "Reflections used",
            "Degrees of freedom",
            "Lattice-fit reduced χ²",
            "Strongest parameter correlation",
            "Systematic correction",
            "Change from the reference cell (a)",
        ):
            assert required in labels
        assert cohen["summary"].startswith("a = 3.533881 ± 0.000015 Å")
        assert "precision" in cohen["summary"] and "not its accuracy" in cohen["summary"]

    def test_every_requested_figure_is_drawn(self, cohen: dict) -> None:
        keys = [figure["key"] for stage in cohen["stages"] for figure in stage.get("figures", [])]
        for required in (
            "scan",
            "fitted_peaks",
            "peak_windows",
            "indexing",
            "peak_quality",
            "residuals",
            "normalized_residuals",
            "systematic_correction",
            "correlation",
            "extrapolation",
            "cross_check",
        ):
            assert required in keys
        assert len(keys) == len(set(keys))

    def test_every_figure_is_self_contained_svg_with_a_reading(self, cohen: dict) -> None:
        for stage in cohen["stages"]:
            for figure in stage.get("figures", []):
                ElementTree.fromstring(figure["svg"])
                assert "<text" not in figure["svg"]
                assert figure["caption"] and figure["interpretation"]

    def test_the_two_chi_squared_statistics_are_never_confused(self, cohen: dict) -> None:
        peaks = next(s for s in cohen["stages"] if s["key"] == "peaks")
        assert "not the lattice fit" in peaks["table"]["caption"]
        quality = next(s for s in cohen["stages"] if s["key"] == "peak_quality")
        assert "unrelated to the lattice-fit χ²" in quality["explanation"]
        fit = next(s for s in cohen["stages"] if s["key"] == "lattice_fit")
        assert any(m["label"] == "Lattice-fit reduced χ²" for m in fit["metrics"])

    def test_the_correction_is_not_called_a_displacement(self, cohen: dict) -> None:
        figure = next(
            figure
            for stage in cohen["stages"]
            for figure in stage.get("figures", [])
            if figure["key"] == "systematic_correction"
        )
        assert "not as a measured displacement" in figure["interpretation"]
        text = " ".join([cohen["summary"], *(stage["summary"] for stage in cohen["stages"])])
        assert "displacement of" not in text.lower()

    def test_database_change_is_not_called_a_strain(self, cohen: dict) -> None:
        cell = next(s for s in cohen["stages"] if s["key"] == "cell")
        labels = [metric["label"] for metric in cell["metrics"]]
        assert "Lattice strain along a" not in labels
        assert "Change from the reference cell (a)" in labels
        assert "not a measured elastic strain" in cell["summary"]
        assert cohen["data"]["relative_change_from_reference"] == pytest.approx(
            cohen["data"]["strain_relative_to_reference"]
        )

    def test_precision_is_distinguished_from_accuracy(self, cohen: dict) -> None:
        cell = next(s for s in cohen["stages"] if s["key"] == "cell")
        assert "precision" in cell["explanation"] and "not accuracy" in cell["explanation"]
        assert any("SRM 640" in note for note in cohen["notes"])

    def test_markdown_reads_result_evidence_diagnostics_method_audit(self, cohen: dict) -> None:
        text = result_to_markdown(cohen).decode("utf-8")
        headings = [
            "## Result and reliability",
            "## Warnings",
            "## Result in detail",
            "## Evidence",
            "## Diagnostics",
            "## Method",
            "## Audit details",
        ]
        positions = [text.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "Poor lattice fit" in text
        assert len(re.findall(r"!\[[^\]]+\]\(data:image/svg\+xml;base64,", text)) == 11
        assert "**What it shows.**" in text

    def test_the_bundle_carries_every_figure_as_a_file(self, cohen: dict) -> None:
        archive = zipfile.ZipFile(io.BytesIO(result_to_bundle(cohen)))
        figures = sorted(name for name in archive.namelist() if name.startswith("figures/"))
        assert len(figures) == 11
        assert "figures/normalized-residuals.svg" in figures


class TestLeBail:
    def test_le_bail_shows_profiles_not_invented_peak_residuals(self, le_bail: dict) -> None:
        keys = [figure["key"] for stage in le_bail["stages"] for figure in stage.get("figures", [])]
        assert "le_bail_profile" in keys
        for absent in ("residuals", "normalized_residuals", "indexing", "peak_windows"):
            assert absent not in keys
        columns = [column["key"] for column in le_bail["table"]["columns"]]
        assert "residual_mdeg" not in columns
        assert "normalized_residual" not in columns

    def test_le_bail_headline_names_the_profile_statistics(self, le_bail: dict) -> None:
        labels = [metric["label"] for metric in le_bail["highlights"]]
        assert "Profile-fit reduced χ²" in labels
        assert "R_wp (background removed)" in labels
        assert "Lattice-fit reduced χ²" not in labels

    def test_le_bail_warns_on_a_poor_profile_fit(self, le_bail: dict) -> None:
        assert le_bail["data"]["reduced_chi_squared"] > 3.0
        assert any("whole-pattern fit" in item for item in le_bail["warnings"])


class TestPrecisionIsNotAccuracy:
    """The demonstration scan carries a known detector zero error, so accuracy is checkable.

    These are the statements the documentation makes about it (section 8 of
    docs/site/algorithms/precise_lattice_parameter_determination.md), and why the
    report must never present its standard uncertainty as an accuracy.
    """

    @staticmethod
    def _true_a() -> float:
        from pytex.app.services.xrd import _DEMO_LATTICE_SCALE

        return float(phase_from_request({"builtin": "ni_fcc"})[1].lattice.a) * _DEMO_LATTICE_SCALE

    def test_nelson_riley_is_precise_but_biased_and_says_so(self, cohen: dict) -> None:
        data = cohen["data"]
        error = data["a"] - self._true_a()
        assert abs(error) / self._true_a() > 100e-6
        assert abs(error) > 10.0 * data["a_standard_uncertainty"]
        assert data["reduced_chi_squared"] > 3.0
        assert any(item.startswith("Poor lattice fit") for item in cohen["warnings"])

    def test_the_exact_form_for_a_zero_error_recovers_the_true_cell(self) -> None:
        from pytex.app.services.xrd import _DEMO_ZERO_SHIFT_DEG

        data = _run(extrapolation="cot_theta")["data"]
        assert abs(data["a"] - self._true_a()) < 3.0 * data["a_standard_uncertainty"]
        assert 0.3 < data["reduced_chi_squared"] < 3.0
        # For f = cot(theta) the fitted D is the zero error itself, in radians.
        assert np.rad2deg(data["drift_coefficient"]) == pytest.approx(
            _DEMO_ZERO_SHIFT_DEG, rel=0.02
        )

    def test_an_acceptable_chi_squared_does_not_prove_accuracy(self) -> None:
        data = _run(extrapolation="cos_squared_over_sin")["data"]
        assert data["reduced_chi_squared"] < 3.0
        assert abs(data["a"] - self._true_a()) > 10.0 * data["a_standard_uncertainty"]
