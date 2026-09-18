"""Figures of the Kearns routes: the running sum ends at f, and every route draws its steps."""

from __future__ import annotations

from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.services.kearns_figures import running_kearns_sum


def _figures(result: dict) -> list[dict]:
    figures = list(result.get("figures", []))
    for stage in result.get("stages", []):
        figures.extend(stage.get("figures", []))
    for figure in figures:
        ElementTree.fromstring(figure["svg"])
        assert figure["caption"] and figure["interpretation"]
    return figures


def test_the_running_sum_ends_at_kearns_published_value() -> None:
    """Kearns' Table 3 longitudinal profile integrates to his tabulated 0.488."""

    tilts = np.arange(5.0, 90.0, 10.0)
    density = [3.27, 2.71, 1.69, 1.35, 1.17, 0.97, 0.73, 0.62, 0.55]
    running = running_kearns_sum(tilts, density)
    assert running[-1] == pytest.approx(0.4879, abs=5e-5)
    assert np.all(np.diff(running) >= 0.0)


def test_the_tilt_profile_route_draws_the_profile_it_integrated() -> None:
    result = REGISTRY.call("kearns.from_tilt_profile", {})
    figures = _figures(result)
    assert [figure["key"] for figure in figures] == ["kearns_profile"]
    rows = result["table"]["rows"]
    running = running_kearns_sum(
        [row["polar_deg"] for row in rows], [row["intensity"] for row in rows]
    )
    assert running[-1] == pytest.approx(result["data"]["directions"][0]["f"], rel=1e-12)


def test_the_diffractogram_route_draws_its_profile() -> None:
    result = REGISTRY.call("kearns.from_diffractogram", {})
    assert [figure["key"] for figure in _figures(result)] == ["kearns_profile"]
    profile = result["data"]["profile"]
    running = running_kearns_sum(profile["polar_deg"], profile["intensity"])
    assert running[-1] == pytest.approx(result["data"]["directions"][0]["f"], rel=1e-9)


def test_the_orientation_route_quotes_a_sampling_standard_error() -> None:
    result = REGISTRY.call("kearns.from_orientations", {})
    keys = [figure["key"] for figure in _figures(result)]
    assert keys == ["kearns_triad", "basal_tilt_distribution"]
    for metric, value in zip(
        result["highlights"], [item["f"] for item in result["data"]["directions"]], strict=True
    ):
        mean, error = (float(part) for part in metric["value"].split(" ± "))
        assert mean == pytest.approx(value, abs=5e-5)
        assert 0.0 < error < 0.05


def test_three_sections_draw_every_section_and_the_triad() -> None:
    result = REGISTRY.call("kearns.from_three_sections", {})
    keys = {figure["key"] for figure in _figures(result)}
    for section in ("axial", "radial", "transverse"):
        assert {f"{section}_scan", f"{section}_profile"} <= keys
    assert "kearns_triad" in keys
    for section in result["data"]["sections"]:
        profile = section["profile"]
        running = running_kearns_sum(profile["polar_deg"], profile["intensity"])
        assert running[-1] == pytest.approx(section["f"], rel=1e-9)
