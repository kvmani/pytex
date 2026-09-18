"""Every analysis operation shows the intermediate results its answer rests on.

One test per operation outside XRD and Kearns (which have their own modules):
the figures exist, parse as self-contained SVG, and say what they plot and what
it implies. Where a figure is attached to a stage, it is attached to the stage
whose evidence it is.
"""

from __future__ import annotations

from xml.etree import ElementTree

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
