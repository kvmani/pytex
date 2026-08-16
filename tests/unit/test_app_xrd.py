"""Scientific and wire-contract tests for the powder-XRD workbench module."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError


def simulate(**request: object) -> dict:
    return REGISTRY.call("xrd.powder_pattern", request)


def test_default_nickel_pattern_has_the_textbook_fcc_sequence() -> None:
    result = simulate()
    rows = result["table"]["rows"]
    assert [row["hkl_label"] for row in rows[:4]] == ["(111)", "(200)", "(220)", "(311)"]
    assert rows[0]["two_theta_deg"] == pytest.approx(44.495, abs=0.02)
    assert rows[0]["relative_intensity"] == 1.0


def test_profile_and_reflection_rows_share_one_physical_interval() -> None:
    result = simulate(two_theta_min_deg=40.0, two_theta_max_deg=100.0, resolution_deg=0.05)
    angles = result["data"]["two_theta_deg"]
    intensity = result["data"]["intensity"]
    assert angles[0] == 40.0
    assert angles[-1] == pytest.approx(100.0)
    assert len(angles) == len(intensity)
    assert max(intensity) == pytest.approx(1.0)
    assert all(40.0 <= row["two_theta_deg"] <= 100.0 for row in result["table"]["rows"])


def test_hover_columns_are_the_export_table_columns() -> None:
    result = simulate()
    assert result["data"]["columns"] == result["table"]["columns"]
    declared = {column["key"] for column in result["data"]["columns"]}
    assert all(declared <= set(row) for row in result["data"]["reflections"])


def test_doublet_changes_the_profile_without_duplicating_primary_rows() -> None:
    common = {
        "two_theta_min_deg": 90.0,
        "two_theta_max_deg": 100.0,
        "resolution_deg": 0.01,
        "fwhm_deg": 0.05,
    }
    single = simulate(radiation="cu_ka", **common)
    doublet = simulate(radiation="cu_ka_doublet", **common)
    assert single["data"]["doublet"] is False
    assert doublet["data"]["doublet"] is True
    assert len(single["table"]["rows"]) == len(doublet["table"]["rows"])
    assert single["data"]["intensity"] != doublet["data"]["intensity"]


def test_pseudo_voigt_and_gaussian_are_distinct_runtime_choices() -> None:
    common = {
        "two_theta_min_deg": 40.0,
        "two_theta_max_deg": 50.0,
        "resolution_deg": 0.02,
        "fwhm_deg": 0.2,
    }
    gaussian = simulate(profile="gaussian", **common)
    pseudo_voigt = simulate(profile="pseudo_voigt", pseudo_voigt_eta=0.7, **common)
    assert gaussian["data"]["intensity"] != pseudo_voigt["data"]["intensity"]


def test_hexagonal_labels_use_miller_bravais_notation() -> None:
    result = simulate(phase={"builtin": "zr_hcp"}, two_theta_min_deg=25.0)
    assert any(len(row["hkl_label"].strip("() ").split()) == 4 for row in result["table"]["rows"])


def test_reversed_scan_range_is_a_field_error() -> None:
    with pytest.raises(InvalidInputError) as caught:
        simulate(two_theta_min_deg=80.0, two_theta_max_deg=20.0)
    assert caught.value.details["field"] == "two_theta_max_deg"


def test_every_xrd_example_is_canonical_and_runnable() -> None:
    examples = [example for example in REGISTRY.examples() if example.panel == "xrd"]
    assert {example.id for example in examples} == {
        "xrd.example.nickel_doublet",
        "xrd.example.silicon",
        "xrd.example.molybdenum_nickel",
        "xrd.example.zirconium",
    }
    for example in examples:
        result = REGISTRY.call(example.operation, example.request)
        assert result["table"]["rows"], example.id
        assert result["data"]["two_theta_deg"], example.id


def test_shared_frontend_exposes_xrd_plot_and_live_appearance_controls() -> None:
    static = Path("src/pytex/app/static")
    main = (static / "js" / "main.js").read_text(encoding="utf-8")
    panel = (static / "js" / "panels" / "xrd.js").read_text(encoding="utf-8")
    assert "import * as xrd" in main
    # Mounted in the tab bar. Asserted as membership rather than as a substring
    # of the whole list: pinning the neighbours makes inserting an unrelated
    # workspace fail an XRD test, which is what happened when EBSD was added
    # between XRD and Variants.
    panels = re.search(r"const PANELS = \[([^\]]*)\]", main)
    assert panels is not None, "main.js must declare a PANELS list"
    assert "xrd" in [name.strip() for name in panels.group(1).split(",")]
    assert "xrd.powder_pattern" in panel
    assert "frame.hoverable(hit, reflection, data.columns)" in panel
    assert "Display controls redraw the existing profile" in panel
    for control in (
        "Profile colour",
        "Reflection-stick colour",
        "Profile line width",
        "Vertical display scale",
        "Fill below profile",
        "Show reflection sticks",
        "Label strong peaks",
        "Peak-label threshold",
        "Reset XRD appearance",
    ):
        assert control in panel
    assert panel.count("call(operation.id") == 1
