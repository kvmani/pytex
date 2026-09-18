"""Every analysis operation shows the intermediate results its answer rests on.

One test per operation outside XRD and Kearns (which have their own modules):
the figures exist, parse as self-contained SVG, and say what they plot and what
it implies. Where a figure is attached to a stage, it is attached to the stage
whose evidence it is.
"""

from __future__ import annotations

from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.app import REGISTRY


def _figures_by_stage(result: dict) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {"": [figure["key"] for figure in result.get("figures", [])]}
    for stage in result.get("stages", []):
        found[stage["key"]] = [figure["key"] for figure in stage.get("figures", [])]
    for figures in (
        result.get("figures", []),
        *(s.get("figures", []) for s in result.get("stages", [])),
    ):
        for figure in figures:
            ElementTree.fromstring(figure["svg"])
            assert "<text" not in figure["svg"]
            assert figure["caption"] and figure["interpretation"]
    return found


def test_texture_analysis_tests_its_odf_against_the_measurement() -> None:
    found = _figures_by_stage(REGISTRY.call("texture.analysis", {}))
    assert found["recalculated"] == ["odf_parity", "odf_tilt_misfit"]
    assert found["fractions"] == ["component_fractions"]


def _plate() -> dict:
    return REGISTRY.call("tem.gallery_pattern", {"pattern": "fcc_al_001"})


def test_a_lattice_fit_shows_each_pick_against_its_node() -> None:
    result = REGISTRY.call("tem.fit_lattice", {"picks": _plate()["data"]["suggested_picks"]})
    assert _figures_by_stage(result)[""] == ["lattice_fit_residuals"]


def test_an_indexed_pattern_shows_its_spacings_and_its_rivals() -> None:
    opened = _plate()
    calibration = opened["data"]["calibration"]
    result = REGISTRY.call(
        "tem.solve_pattern",
        {
            "phase": calibration["phase"],
            "picks": opened["data"]["suggested_picks"],
            "units": "px",
            "camera_constant_mm_angstrom": calibration["camera_constant_mm_angstrom"],
            "pixel_size_mm": calibration["pixel_size_mm"],
        },
    )
    assert _figures_by_stage(result)[""] == ["solve_spacings", "solve_candidates"]


def test_an_hrtem_simulation_shows_what_the_lens_passed() -> None:
    from pytex.app.services.tem_figures import radial_average

    result = REGISTRY.call("tem.simulate_hrem", {})
    assert _figures_by_stage(result)[""] == ["hrem_spectrum_vs_ctf"]
    # The rotational average of a radially symmetric image is its radial profile:
    # for f = |q| each ring's mean is its centre, to within the ring's width.
    size, step = 64, 0.5
    y, x = np.indices((size, size)) - size // 2
    q = np.hypot(x / (size * step), y / (size * step))
    centres, average = radial_average(q, pixel_size_angstrom=step, bins=16)
    width = float(centres[1] - centres[0])
    assert np.allclose(average, centres, atol=0.3 * width)


def test_a_cbed_thickness_fit_states_its_uncertainty_when_it_has_one() -> None:
    from pytex.app.fitstats import fit_line

    three = REGISTRY.call("cbed.thickness_from_fringes", {"s1": 0.0071, "s2": 0.0128, "s3": 0.0165})
    assert _figures_by_stage(three)[""] == ["thickness_fit"]
    orders = np.asarray(three["data"]["orders"], dtype=float)
    values = np.asarray(three["data"]["excitation_errors_inv_angstrom"], dtype=float)
    line = fit_line(1.0 / orders**2, (values / orders) ** 2)
    thickness = three["data"]["thickness_angstrom"]
    expected = thickness * line.sigma_intercept / (2.0 * line.intercept) / 10.0
    stated = three["highlights"][0]["value"]
    assert float(stated.split(" ± ")[1]) == pytest.approx(expected, abs=0.006)

    two = REGISTRY.call("cbed.thickness_from_fringes", {})
    assert "no uncertainty from two minima" in two["highlights"][0]["value"]


def test_orientation_relationships_show_pairs_and_catalogue_margins() -> None:
    for operation, expected in (
        ("ebsd.or_from_grains", ["or_catalog_distances", "or_pair_residuals"]),
        ("variants.or_from_grains", ["or_catalog_distances"]),
    ):
        example = next(item for item in REGISTRY.examples() if item.operation == operation)
        result = REGISTRY.call(operation, example.request)
        assert _figures_by_stage(result)[""] == expected


def test_an_ebsd_scan_summary_shows_what_its_threshold_keeps() -> None:
    result = REGISTRY.call("ebsd.scan_summary", {})
    assert _figures_by_stage(result)[""] == ["scan_quality"]
    assert "confidence-index threshold" in result["figures"][0]["caption"]
    # One histogram per quality channel present in the scan, plus the grain sizes.
    assert set(result["data"]["channels"]) == {"confidence_index", "fit", "image_quality"}
