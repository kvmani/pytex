"""Scientific and wire-contract tests for the powder-XRD workbench module."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.core.fixtures import get_phase_fixture


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
    "xrd.lattice_parameters": "abscissa",
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
        "xrd.example.lattice_average_fails",
        "xrd.example.lattice_cohen_extrapolation",
        "xrd.example.lattice_hexagonal_le_bail",
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


# ---------------------------------------------------------------------------
# Determine lattice parameters
# ---------------------------------------------------------------------------

#: The demonstration scan is built with the cell dilated by this factor, so the
#: answer every method below is chasing is known before the call is made.
DEMONSTRATION_DILATION = 1.003


def determine(**request: object) -> dict:
    base: dict[str, object] = {"phase": {"builtin": "ni_fcc"}, "radiation": "cu_ka_doublet"}
    base.update(request)
    return REGISTRY.call("xrd.lattice_parameters", base)


def nickel_target() -> float:
    from pytex.app.phases import builtin_phase

    return builtin_phase("ni_fcc").to_phase().lattice.a * DEMONSTRATION_DILATION


def test_matching_the_extrapolation_to_the_aberration_beats_averaging() -> None:
    """The argument the view exists to make, on a scan whose answer is known.

    The demonstration scan carries a 0.05 degree detector zero error. A constant
    angular offset gives Delta(sin^2 theta) = sin^2(theta) cot(theta) Delta(2
    theta), so cot(theta) is the extrapolation function that matches it; the
    others are the wrong shape and leave most of the error behind.
    """

    target = nickel_target()

    def error(**request: object) -> float:
        return abs(determine(**request)["data"]["a"] - target) / target

    naive = error(method="average", extrapolation="none")
    uncorrected = error(method="cohen", extrapolation="none")
    mismatched = error(method="cohen", extrapolation="nelson_riley")
    matched = error(method="cohen", extrapolation="cot_theta")

    assert naive > 3.0e-4
    assert uncorrected < naive
    assert mismatched < uncorrected
    assert matched < 1.0e-5
    assert matched < mismatched / 20.0


def test_the_goodness_of_fit_tracks_the_accuracy() -> None:
    """A user cannot check against a known answer on real data, so chi-squared
    has to be the thing that tells them which method worked."""

    target = nickel_target()
    rows = [
        determine(method="cohen", extrapolation=function)["data"]
        for function in ("none", "nelson_riley", "cot_theta")
    ]
    errors = [abs(row["a"] - target) for row in rows]
    chi_squared = [row["reduced_chi_squared"] for row in rows]
    assert errors == sorted(errors, reverse=True)
    assert chi_squared == sorted(chi_squared, reverse=True)


def test_re_indexing_recovers_the_high_angle_reflections() -> None:
    """A starting cell 0.3 per cent wrong misplaces a 121 degree line by 0.6
    degrees, well outside the indexing tolerance. Without a second indexing
    pass the tool would silently drop exactly the reflections that carry the
    precision, and would therefore require the answer in advance."""

    result = determine(method="cohen", extrapolation="cot_theta")
    rows = result["table"]["rows"]
    assert len(rows) == 6
    assert [row["hkl_label"] for row in rows[:3]] == ["(111)", "(200)", "(220)"]
    assert rows[-1]["two_theta_observed_deg"] > 115.0


def test_the_average_method_refuses_a_hexagonal_cell_as_a_field_error() -> None:
    with pytest.raises(InvalidInputError) as raised:
        determine(phase={"builtin": "ti_hcp"}, method="average")
    assert raised.value.details["field"] == "method"
    assert "cubic" in str(raised.value)


def test_a_hexagonal_determination_reports_both_parameters_and_their_ratio() -> None:
    from pytex.app.phases import builtin_phase

    lattice = builtin_phase("ti_hcp").to_phase().lattice
    data = determine(phase={"builtin": "ti_hcp"}, method="le_bail", systematic="zero")["data"]
    assert data["crystal_system"] == "hexagonal"
    assert data["a"] == pytest.approx(lattice.a * DEMONSTRATION_DILATION, rel=5.0e-5)
    assert data["c"] == pytest.approx(lattice.c * DEMONSTRATION_DILATION, rel=5.0e-5)
    assert data["axial_ratio"] == pytest.approx(lattice.c / lattice.a, rel=1.0e-4)
    # No per-reflection lattice parameter exists here, so the view must not
    # draw an extrapolation plot that pretends one does. A whole-pattern fit
    # measures no individual position either, so its picture is the difference
    # curve rather than a residual scatter.
    assert data["plot_kind"] == "profile"
    assert len(data["abscissa"]) == len(data["calculated"]) == len(data["difference"])
    assert data["weighted_profile_r"] is not None


def test_le_bail_does_not_report_positions_it_never_measured() -> None:
    """A whole-pattern fit measures no individual peak position, so the table
    must not carry a residual or an uncertainty column full of zeros."""

    result = determine(method="le_bail", systematic="zero")
    keys = {column["key"] for column in result["table"]["columns"]}
    assert "residual_mdeg" not in keys
    assert "standard_uncertainty_mdeg" not in keys
    assert "two_theta_calculated_deg" in keys
    assert result["data"]["columns"] == result["table"]["columns"]


def test_the_drift_column_recovers_the_zero_error_that_was_injected() -> None:
    """The strongest available check that the drift machinery is right.

    A detector zero error is a *constant* Delta(2 theta), and the cot(theta)
    extrapolation function is derived to be exactly the one that represents it.
    So if the algebra is right, the angular shift the drift term removes must
    come back constant across the pattern and equal to the 0.05 degree zero the
    demonstration scan was built with -- neither of which was fitted for
    directly.
    """

    rows = determine(method="cohen", extrapolation="cot_theta")["table"]["rows"]
    shifts = [row["systematic_shift_mdeg"] for row in rows]
    sigmas = [row["standard_uncertainty_mdeg"] for row in rows]
    assert all(shift == pytest.approx(shifts[0], rel=1.0e-6) for shift in shifts)
    assert abs(shifts[0]) == pytest.approx(50.0, abs=2.0)
    # And it must dwarf the position uncertainties, or it would not be worth
    # refining in the first place.
    assert max(sigmas) > 0.0
    assert abs(shifts[0]) > 50.0 * max(sigmas)


def test_a_displacement_correction_falls_towards_back_reflection() -> None:
    """The other drift shape, and the reason extrapolation works at all."""

    rows = determine(
        specimen_displacement_mm=0.2, method="cohen", extrapolation="cos_squared_over_sin"
    )["table"]["rows"]
    shifts = [abs(row["systematic_shift_mdeg"]) for row in rows]
    assert shifts[0] > shifts[-1]
    assert shifts[-1] < 0.6 * shifts[0]


def test_no_correction_leaves_the_drift_column_empty_and_says_so() -> None:
    result = determine(method="cohen", extrapolation="none")
    assert all(row["systematic_shift_mdeg"] == 0.0 for row in result["table"]["rows"])
    assert "inside the quoted cell" in result["data"]["describe"]


def test_an_injected_displacement_is_declared_in_the_notes() -> None:
    result = determine(specimen_displacement_mm=0.1)
    joined = " ".join(result["notes"])
    assert "100 \u00b5m specimen displacement" in joined
    assert "two aberrations of different angular form" in joined
    assert result["data"]["synthetic"] is True


def test_the_result_refuses_to_call_a_lattice_parameter_a_stress() -> None:
    result = determine()
    joined = " ".join(result["notes"]) + result["data"]["describe"]
    assert "not a stress" in joined
    assert "X-ray elastic constants" in joined


def test_the_plot_payload_carries_what_the_view_draws() -> None:
    data = determine(method="cohen", extrapolation="cot_theta")["data"]
    assert data["plot_kind"] == "extrapolation"
    assert len(data["abscissa"]) == len(data["ordinate"]) == data["reflection_count"]
    assert data["abscissa_label"]
    assert data["ordinate_label"]
    # The intercept of the drawn line and the determined value must agree, or
    # the picture would contradict the number printed beside it.
    assert data["line_intercept"] == pytest.approx(data["a"], rel=5.0e-5)


def test_a_pasted_scan_is_determined_without_a_demonstration_claim() -> None:
    simulated = REGISTRY.call(
        "xrd.powder_pattern",
        {
            "phase": {"builtin": "ni_fcc"},
            "radiation": "cu_ka_doublet",
            "two_theta_min_deg": 40.0,
            "two_theta_max_deg": 145.0,
            "profile": "pseudo_voigt",
            "fwhm_deg": 0.12,
            "resolution_deg": 0.01,
        },
    )["data"]
    lines = [
        f"{angle} {20000.0 * value + 150.0}"
        for angle, value in zip(simulated["two_theta_deg"], simulated["intensity"], strict=True)
    ]
    result = determine(
        data_source="paste",
        scan="\n".join(lines),
        method="cohen",
        extrapolation="nelson_riley",
        expected_fwhm_deg=0.12,
    )
    from pytex.app.phases import builtin_phase

    assert result["data"]["synthetic"] is False
    assert result["data"]["a"] == pytest.approx(
        builtin_phase("ni_fcc").to_phase().lattice.a, rel=1.0e-4
    )


def test_the_panel_opens_on_the_view_it_lists_first() -> None:
    """`xrd.js` opens on `examples[0]`, so registration order decides the
    landing view. Adding examples for a new view put them at the front of the
    manifest and silently changed the XRD panel from opening on the simulation
    to opening on a determination -- the same failure the Variants panel had
    when its example list was reordered. The order is a contract, so it is
    pinned here rather than remembered.
    """

    panel_source = (
        Path("src/pytex/app/static/js/panels/xrd.js").read_text(encoding="utf-8")
    )
    assert "loadExample(examples[0])" in panel_source, (
        "the panel no longer opens on the first example; update this guard to match"
    )
    views = re.search(r"const VIEWS = \[(.*?)\];", panel_source, re.S)
    assert views is not None
    first_view = re.findall(r"'([^']+)'", views.group(1))[0]
    examples = [example for example in REGISTRY.examples() if example.panel == "xrd"]
    assert examples[0].operation == first_view


def test_uncertainties_stay_with_their_reflections_under_an_angular_floor() -> None:
    """An angular floor filters the determination but not the indexing behind it.

    The table looks up each position uncertainty by angle rather than by row
    position for exactly this case: a positional lookup would attach every
    uncertainty to the wrong reflection the moment an operator restricts the
    range -- which is the case they would restrict it for.
    """

    everything = determine(method="cohen", extrapolation="cot_theta")
    restricted = determine(
        method="cohen", extrapolation="cot_theta", minimum_two_theta_deg=90.0
    )
    assert len(restricted["table"]["rows"]) < len(everything["table"]["rows"])
    assert all(
        row["two_theta_observed_deg"] >= 90.0 for row in restricted["table"]["rows"]
    )

    # Every surviving row must keep the uncertainty it had in the full table.
    full = {
        round(row["two_theta_observed_deg"], 6): row["standard_uncertainty_mdeg"]
        for row in everything["table"]["rows"]
    }
    for row in restricted["table"]["rows"]:
        assert row["standard_uncertainty_mdeg"] > 0.0
        assert row["standard_uncertainty_mdeg"] == pytest.approx(
            full[round(row["two_theta_observed_deg"], 6)]
        )


ROOT = Path(__file__).resolve().parents[2]
NI_FCC_XY = ROOT / "fixtures" / "diffraction" / "experimental_ni_fcc_pattern.xy"
NI_FCC_XRDML = ROOT / "fixtures" / "diffraction" / "experimental_ni_fcc_pattern.xrdml"


def test_operations_accept_experimental_pattern_file_xy() -> None:
    assert NI_FCC_XY.exists()
    text = NI_FCC_XY.read_text(encoding="utf-8")
    scan_file = {"name": "experimental_ni_fcc_pattern.xy", "text": text}

    # 1. Background estimation
    bg = REGISTRY.call(
        "xrd.background", {"phase": {"builtin": "ni_fcc"}, "scan_file": scan_file}
    )
    assert bg["data"]["synthetic"] is False
    assert len(bg["data"]["two_theta_deg"]) == 4001

    # 2. Lattice parameters
    lp = REGISTRY.call(
        "xrd.lattice_parameters",
        {"phase": {"builtin": "ni_fcc"}, "scan_file": scan_file},
    )
    assert lp["data"]["synthetic"] is False
    assert lp["data"]["a"] == pytest.approx(3.52387, rel=1e-3)

    # 3. Rietveld refinement
    riet = REGISTRY.call(
        "xrd.rietveld",
        {"phase": {"builtin": "ni_fcc"}, "scan_file": scan_file},
    )
    assert riet["data"]["synthetic"] is False
    assert riet["data"]["weighted_profile_r_factor"] < 0.45


def test_operations_accept_experimental_pattern_file_xrdml() -> None:
    assert NI_FCC_XRDML.exists()
    text = NI_FCC_XRDML.read_text(encoding="utf-8")
    scan_file = {"name": "experimental_ni_fcc_pattern.xrdml", "text": text}

    bg = REGISTRY.call(
        "xrd.background", {"phase": {"builtin": "ni_fcc"}, "scan_file": scan_file}
    )
    assert bg["data"]["synthetic"] is False
    assert len(bg["data"]["two_theta_deg"]) == 4001


def test_operation_rejects_missing_scan_file_when_source_is_file() -> None:
    with pytest.raises(InvalidInputError) as caught:
        REGISTRY.call(
            "xrd.background", {"phase": {"builtin": "ni_fcc"}, "data_source": "file"}
        )
    assert caught.value.details["field"] == "scan_file"


def test_pattern_controls_are_wired_in_panel_js() -> None:
    static = Path("src/pytex/app/static/js")
    panel = (static / "panels" / "xrd.js").read_text(encoding="utf-8")
    xrdscan = (static / "core" / "xrdscan.js").read_text(encoding="utf-8")

    assert "patternControls" in panel
    assert "adoptForm" in panel
    assert "withPattern" in panel
    assert "PATTERN_OPERATIONS" in panel
    assert "Experimental pattern" in xrdscan
    assert ".xy,.xrdml,.csv,.dat,.txt" in xrdscan



# ---------------------------------------------------------------------------
# Phase identification: the one operation told several phases rather than one
# ---------------------------------------------------------------------------


def identify(**request: object) -> dict:
    return REGISTRY.call("xrd.phase_identification", request)


def _builtin(*identifiers: str) -> dict:
    return {"phases": [{"phase": {"builtin": identifier}} for identifier in identifiers]}


def _from_cif(identifier: str, name: str) -> dict:
    """Return a candidate entry carrying CIF text, as the browser sends one."""

    return {
        "label": name,
        "phase": {"cif": {"name": name, "text": get_phase_fixture(identifier).read_cif_text()}},
    }


def test_the_demonstration_specimen_is_recovered_from_a_list_of_candidates() -> None:
    """The synthetic scan is made of nickel, so the ranking must return nickel.

    The demonstration scan carries a cell dilated by 1.003 on top of a detector
    zero error, which is the whole point: a tabulated CIF never matches a real
    specimen exactly, and an identification that only worked on an undilated
    cell would not survive first contact with a laboratory.
    """

    result = identify(candidates=_builtin("ni_fcc", "cu_fcc", "fe_bcc", "nacl"))
    assert result["data"]["best_phase_name"] == "Nickel (fcc)"
    assert result["data"]["is_conclusive"] is True
    assert result["data"]["is_decisive"] is True
    assert result["table"]["rows"][0]["phase_name"] == "Nickel (fcc)"


def test_candidates_may_be_uploaded_as_cif_files() -> None:
    """The point of the feature: the candidates are the user's own structure files."""

    result = identify(
        candidates={
            "phases": [
                _from_cif("ni_fcc", "nickel.cif"),
                _from_cif("fe_bcc", "ferrite.cif"),
            ]
        }
    )
    assert result["data"]["best_phase_name"] == "nickel.cif"
    rows = {row["phase_name"]: row for row in result["table"]["rows"]}
    assert set(rows) == {"nickel.cif", "ferrite.cif"}
    # Provenance survives into the table, so a ranking can be traced to the
    # files that produced it rather than only to the labels on screen.
    assert "nickel.cif" in rows["nickel.cif"]["source"]


def test_built_in_and_uploaded_candidates_compete_in_one_ranking() -> None:
    result = identify(
        candidates={
            "phases": [
                {"phase": {"builtin": "fe_bcc"}},
                _from_cif("ni_fcc", "unknown_powder.cif"),
                {"phase": {"builtin": "nacl"}},
            ]
        }
    )
    assert result["data"]["best_phase_name"] == "unknown_powder.cif"


def test_the_table_reports_every_criterion_and_not_only_the_total() -> None:
    """Which criterion a candidate fails is the diagnosis; the total is not."""

    result = identify(candidates=_builtin("ni_fcc", "fe_bcc"))
    keys = {column["key"] for column in result["table"]["columns"]}
    assert {"score", "explained", "completeness", "position", "intensity"} <= keys
    assert result["data"]["columns"] == result["table"]["columns"]
    for row in result["table"]["rows"]:
        assert set(keys) <= set(row)


def test_a_centring_is_separated_by_lines_seen_rather_than_by_position() -> None:
    """A bcc candidate on an fcc pattern must lose on absent lines, not on angles."""

    result = identify(candidates=_builtin("ni_fcc", "fe_bcc"))
    rows = {row["phase_name"]: row for row in result["table"]["rows"]}
    assert rows["Nickel (fcc)"]["completeness"] > rows["Ferrite (bcc Fe)"]["completeness"]


def test_the_refined_cell_dilation_is_reported_for_every_candidate() -> None:
    """A candidate stretched to the edge of the search range has to be visible.

    The dilation is what makes a tabulated CIF usable on a real specimen, and
    it is also the thing that could flatter a wrong candidate, so it is a
    reported column rather than a hidden step.
    """

    result = identify(candidates=_builtin("ni_fcc", "cu_fcc"), cell_scale_range=0.02)
    rows = {row["phase_name"]: row for row in result["table"]["rows"]}
    # The demonstration scan dilates the cell by 1.003, so the true phase needs
    # about +0.3 per cent and no more.
    assert rows["Nickel (fcc)"]["cell_dilation_percent"] == pytest.approx(0.3, abs=0.1)
    # Copper's cell is 2.6 per cent larger than nickel's, so the search pins it
    # at the edge of the two per cent range: it is being stretched as far as it
    # is allowed and still does not fit. That is the reading the column exists
    # to make possible, and it is why the dilation is reported rather than
    # silently applied.
    assert rows["Copper (fcc)"]["cell_dilation_percent"] == pytest.approx(-2.0, abs=0.01)


def test_every_indexed_candidate_travels_as_an_overlay_for_the_plot() -> None:
    """The drawing has to let a reader check the ranking, not only read it."""

    result = identify(candidates=_builtin("ni_fcc", "fe_bcc", "nacl"))
    overlays = result["data"]["overlays"]
    assert overlays[0]["phase_name"] == "Nickel (fcc)"
    for entry in overlays:
        assert len(entry["two_theta_deg"]) == len(entry["labels"])
        assert len(entry["two_theta_deg"]) == len(entry["relative_intensity"])


def test_the_fitted_peaks_travel_with_the_result() -> None:
    """Peak detection is where an identification most often goes wrong."""

    result = identify(candidates=_builtin("ni_fcc"))
    peaks = result["data"]["peaks"]
    assert peaks and len(peaks) == result["data"]["peak_count"]
    assert all(peak["two_theta_deg"] > 0.0 for peak in peaks)


def test_one_candidate_is_a_check_rather_than_an_identification_and_says_so() -> None:
    result = identify(candidates=_builtin("ni_fcc"))
    assert "a check on one phase" in result["data"]["describe"]


def test_a_pattern_no_candidate_explains_is_reported_as_inconclusive() -> None:
    result = identify(candidates=_builtin("quartz_alpha", "alpha_u"))
    assert result["data"]["is_conclusive"] is False
    assert "No candidate accounts for this pattern" in result["summary"]


def test_the_textured_weighting_discards_intensity_evidence() -> None:
    """A rolled sheet's intensities are moved by orientation, not by structure."""

    balanced = identify(candidates=_builtin("ni_fcc", "fe_bcc"), weighting="standard")
    textured = identify(candidates=_builtin("ni_fcc", "fe_bcc"), weighting="textured")
    assert balanced["data"]["best_score"] != pytest.approx(textured["data"]["best_score"])
    assert textured["data"]["best_phase_name"] == "Nickel (fcc)"


def test_an_empty_candidate_list_is_refused_beside_the_right_control() -> None:
    with pytest.raises(InvalidInputError) as error:
        identify(candidates={"phases": []})
    assert error.value.details["field"] == "candidates"


def test_an_unreadable_candidate_names_itself_in_the_error() -> None:
    """A user with five CIFs open must be told which one failed, not that one did."""

    with pytest.raises(InvalidInputError) as error:
        identify(
            candidates={
                "phases": [
                    {"phase": {"builtin": "ni_fcc"}},
                    {"phase": {"cif": {"name": "broken.cif", "text": "not a cif"}}},
                ]
            }
        )
    assert error.value.details["field"] == "candidates"
    assert "Candidate 2" in error.value.message


def test_the_operation_declares_its_search_match_citations() -> None:
    spec = REGISTRY.get("xrd.phase_identification")
    joined = " ".join(spec.citations)
    assert "10.1021/ac50125a001" in joined  # Hanawalt, Rinn & Frevel
    assert "10.1107/S0021889886089458" in joined  # Dollase, on why intensities are weak evidence
