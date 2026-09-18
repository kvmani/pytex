# ruff: noqa: RUF002
"""Figures of the other XRD operations, and the statistics they state.

Each operation gains figures of its intermediate results. Two of them state a
derived statistic that the operation itself did not report - the weighted
Rietveld residual and the ordinary least-squares uncertainty of the
Williamson-Hall line - and those are checked against an independent
computation here. The rest are checked for being present, parseable and
self-contained.
"""

from __future__ import annotations

from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.services.xrd_figures import williamson_hall_uncertainties

_NICKEL = {"phase": {"builtin": "ni_fcc"}}


def _figure_keys(result: dict) -> list[str]:
    keys = [figure["key"] for figure in result.get("figures", [])]
    for figure in result.get("figures", []):
        ElementTree.fromstring(figure["svg"])
        assert figure["caption"] and figure["interpretation"]
    return keys


def test_background_draws_the_estimate_under_the_scan() -> None:
    result = REGISTRY.call("xrd.background", _NICKEL)
    assert _figure_keys(result) == ["background"]


def test_the_simulated_pattern_shows_what_sets_each_intensity() -> None:
    result = REGISTRY.call("xrd.powder_pattern", _NICKEL)
    assert _figure_keys(result) == ["powder_pattern", "intensity_factors"]


def test_phase_identification_draws_every_candidates_lines() -> None:
    request = {
        **_NICKEL,
        "candidates": {
            "phases": [
                {"phase": {"builtin": identifier}} for identifier in ("ni_fcc", "cu_fcc", "fe_bcc")
            ]
        },
    }
    result = REGISTRY.call("xrd.phase_identification", request)
    assert _figure_keys(result) == ["identification_lines", "identification_scores"]


class TestRietveld:
    @pytest.fixture(scope="class")
    def refined(self) -> dict:
        return REGISTRY.call("xrd.rietveld", _NICKEL)

    def test_the_fit_is_drawn_with_its_residuals_and_shifts(self, refined: dict) -> None:
        assert _figure_keys(refined) == [
            "rietveld_profile",
            "rietveld_weighted_residuals",
            "rietveld_parameter_shifts",
        ]

    def test_weighted_residuals_reproduce_the_goodness_of_fit(self) -> None:
        """Σ[(y_obs − y_calc)√w]² / (N − P) is the square of the reported goodness of fit."""

        from pytex.app.phases import phase_from_request
        from pytex.app.services.xrd import _RADIATION, _measured_from_request
        from pytex.app.services.xrd_figures import rietveld_weighted_residuals
        from pytex.diffraction.rietveld import refine_rietveld

        request = REGISTRY.get("xrd.rietveld").bind(_NICKEL)
        radiation = _RADIATION[str(request["radiation"])]()
        _, phase = phase_from_request(request["phase"])
        measured, _ = _measured_from_request(request, phase, radiation)
        result = refine_rietveld(measured, phase, radiation=radiation)
        weighted = rietveld_weighted_residuals(result)
        assert weighted is not None
        refined_count = sum(1 for item in result.parameters if item.refined)
        degrees = weighted.size - refined_count
        assert float(np.sum(weighted**2) / degrees) == pytest.approx(
            result.goodness_of_fit**2, rel=1e-6
        )


class TestSizeStrain:
    def test_the_line_uncertainties_match_numpy_covariance(self) -> None:
        rng = np.random.default_rng(3)
        x = np.linspace(0.8, 2.8, 9)
        y = 0.004 + 0.0017 * x + rng.normal(0.0, 5e-5, x.size)
        stats = williamson_hall_uncertainties(x, y, wavelength_angstrom=1.5406, shape_factor=0.9)
        (slope, intercept), covariance = np.polyfit(x, y, 1, cov=True)
        assert stats["slope"] == pytest.approx(slope)
        assert stats["intercept"] == pytest.approx(intercept)
        assert stats["sigma_slope"] == pytest.approx(np.sqrt(covariance[0, 0]), rel=1e-9)
        assert stats["sigma_intercept"] == pytest.approx(np.sqrt(covariance[1, 1]), rel=1e-9)
        assert stats["covariance"] == pytest.approx(covariance[0, 1], rel=1e-9)
        # D = K lambda / intercept, in nm from angstrom, and its relative error.
        assert stats["size_nm"] == pytest.approx(0.9 * 1.5406 / intercept / 10.0)
        assert stats["sigma_size_nm"] / stats["size_nm"] == pytest.approx(
            stats["sigma_intercept"] / intercept
        )

    def test_two_points_have_no_scatter_to_estimate(self) -> None:
        stats = williamson_hall_uncertainties(
            np.array([1.0, 2.0]),
            np.array([0.005, 0.006]),
            wavelength_angstrom=1.54,
            shape_factor=0.9,
        )
        assert np.isnan(stats["sigma_slope"]) and np.isnan(stats["sigma_size_nm"])

    def test_the_operation_reports_size_and_strain_with_uncertainties(self) -> None:
        result = REGISTRY.call("xrd.size_strain", {})
        assert _figure_keys(result) == [
            "instrument_calibration",
            "width_decomposition",
            "williamson_hall",
            "scherrer_sizes",
        ]
        labels = [metric["label"] for metric in result["highlights"]]
        assert labels == ["Crystallite size D", "Microstrain ε", "R² of the line"]
        size = float(result["highlights"][0]["value"].split(" ± ")[0])
        assert size == pytest.approx(result["data"]["crystallite_size_nm"], rel=1e-3)
