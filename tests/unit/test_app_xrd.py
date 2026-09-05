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


#: Which key in `data` carries the abscissa each view plots against. The
#: size/strain view has no profile at all -- it plots Williamson-Hall
#: coordinates -- so asserting one shape for every example would either miss
#: that view or force it to fake a profile it does not have.
_EXAMPLE_ABSCISSA = {
    "xrd.powder_pattern": "two_theta_deg",
    "xrd.background": "two_theta_deg",
    "xrd.rietveld": "two_theta_deg",
    "xrd.size_strain": "abscissa",
}


def test_every_xrd_example_is_canonical_and_runnable() -> None:
    examples = [example for example in REGISTRY.examples() if example.panel == "xrd"]
    assert {example.id for example in examples} == {
        "xrd.example.nickel_doublet",
        "xrd.example.silicon",
        "xrd.example.molybdenum_nickel",
        "xrd.example.zirconium",
        "xrd.example.background_snip",
        "xrd.example.rietveld_zero_and_cell",
        "xrd.example.rietveld_width_wrong",
        "xrd.example.size_strain",
    }
    for example in examples:
        result = REGISTRY.call(example.operation, example.request)
        assert result["table"]["rows"], example.id
        assert result["data"][_EXAMPLE_ABSCISSA[example.operation]], example.id


def test_a_generated_scan_is_always_declared_synthetic() -> None:
    """A demonstration scan must never be presentable as a measurement."""

    for operation in ("xrd.background", "xrd.rietveld"):
        result = REGISTRY.call(
            operation, {"phase": {"builtin": "ni_fcc"}, "data_source": "demonstration"}
        )
        assert result["data"]["synthetic"] is True, operation
        assert any("generated, not measured" in note for note in result["notes"]), operation


def test_refinement_recovers_the_cell_and_zero_the_demonstration_scan_was_built_with() -> None:
    """The demonstration scan states its own answer, so the view can be checked."""

    result = REGISTRY.call(
        "xrd.rietveld", {"phase": {"builtin": "ni_fcc"}, "data_source": "demonstration"}
    )
    parameters = {row["parameter"]: row for row in result["table"]["rows"]}
    assert parameters["lattice_scale"]["value"].startswith("1.003")
    assert parameters["zero_shift_deg"]["value"].startswith("0.050")
    assert result["data"]["converged"] is True
    # Poisson noise and a model that is right give a goodness of fit near one.
    assert 0.8 < result["data"]["goodness_of_fit"] < 1.6


def test_holding_the_zero_pushes_the_error_into_the_cell() -> None:
    """The failure the zero-shift control exists to prevent, asserted rather than described.

    A detector zero error displaces every peak by the same angle; a cell dilation
    displaces them in proportion to ``tan(theta)``. Over a wide scan those are
    genuinely different signatures, so the cell cannot absorb the zero error
    *silently* -- it takes a wrong value and the fit degrades as well. Both halves
    are asserted, because the second is what makes the first diagnosable.
    """

    common = {"phase": {"builtin": "ni_fcc"}, "data_source": "demonstration"}
    honest = REGISTRY.call("xrd.rietveld", {**common, "refine_zero_shift": True})
    absorbed = REGISTRY.call("xrd.rietveld", {**common, "refine_zero_shift": False})
    assert abs(absorbed["data"]["refined_lattice_a"] - honest["data"]["refined_lattice_a"]) > 1e-4
    assert (
        absorbed["data"]["weighted_profile_r_factor"]
        > 2.0 * honest["data"]["weighted_profile_r_factor"]
    )
    # And the residual stops looking like noise, which is the general signature
    # of a parameter absorbing something that is not its own.
    assert absorbed["data"]["durbin_watson"] < honest["data"]["durbin_watson"]


def test_a_pasted_scan_is_read_and_a_malformed_one_is_refused() -> None:
    generated = REGISTRY.call(
        "xrd.background", {"phase": {"builtin": "ni_fcc"}, "data_source": "demonstration"}
    )
    angles = generated["data"]["two_theta_deg"]
    observed = generated["data"]["observed"]
    pasted = "\n".join(
        f"{angle} {intensity}" for angle, intensity in zip(angles, observed, strict=True)
    )
    result = REGISTRY.call(
        "xrd.background",
        {
            "phase": {"builtin": "ni_fcc"},
            "data_source": "paste",
            "scan": "# a header line, which the reader must skip\n" + pasted,
        },
    )
    assert result["data"]["synthetic"] is False
    assert len(result["data"]["two_theta_deg"]) == len(angles)

    with pytest.raises(InvalidInputError) as caught:
        REGISTRY.call(
            "xrd.background",
            {
                "phase": {"builtin": "ni_fcc"},
                "data_source": "paste",
                "scan": "30.0\n30.02\n",
            },
        )
    assert caught.value.details["field"] == "scan"


def test_size_and_strain_defaults_reproduce_their_stated_answer() -> None:
    """The untouched form returns 25 nm and 0.2%, which is what its help text claims."""

    result = REGISTRY.call("xrd.size_strain", {})
    assert result["data"]["crystallite_size_nm"] == pytest.approx(25.0, abs=1e-3)
    assert result["data"]["microstrain"] == pytest.approx(0.002, abs=1e-6)
    assert result["data"]["r_squared"] > 0.999


def test_a_specimen_sharper_than_the_instrument_is_refused() -> None:
    with pytest.raises(InvalidInputError) as caught:
        REGISTRY.call(
            "xrd.size_strain",
            {"sample_peaks": "21.36  0.0100\n53.99  0.0100\n87.79  0.0100"},
        )
    assert caught.value.details["field"] == "sample_peaks"


def test_shared_frontend_exposes_xrd_plot_and_live_appearance_controls() -> None:
    static = Path("src/pytex/app/static")
    main = (static / "js" / "main.js").read_text(encoding="utf-8")
    panel = (static / "js" / "panels" / "xrd.js").read_text(encoding="utf-8")
    assert "import * as xrd" in main
    # Mounted in the tab bar. Asserted as membership rather than as a substring
    # of the whole list: pinning the neighbours makes inserting an unrelated
    # workspace fail an XRD test, which is what happened when EBSD was added
    # between XRD and Variants.
    workspaces = re.search(r"const WORKSPACES = \[([^\]]*)\]", main)
    assert workspaces is not None, "main.js must declare a WORKSPACES list"
    assert "solo(xrd)" in [entry.strip() for entry in workspaces.group(1).split(",")]
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
    # Still exactly one call site. More than one would let two views drift apart
    # in how they handle errors.
    assert panel.count("call(launched.id") == 1
    # And it dispatches on the view the run was *launched* for, not on whichever
    # view happens to be selected when the answer arrives. A request in flight
    # across a view change would otherwise be drawn by a renderer expecting
    # different keys, which is a type error rather than a wrong picture.
    assert "const launched = state.operation;" in panel
    assert panel.count("if (state.operation !== launched) return;") == 2
    # A view the panel does not list is a view no user can reach, so every
    # registered xrd operation must appear in the panel's own VIEWS list.
    views = re.search(r"const VIEWS = \[([^\]]*)\]", panel)
    assert views is not None, "xrd.js must declare a VIEWS list"
    listed = {entry.strip().strip("',\"") for entry in views.group(1).split("\n") if entry.strip()}
    registered = {operation.id for operation in REGISTRY.operations() if operation.panel == "xrd"}
    assert registered <= listed, registered - listed
