"""The SAED simulator: the forward pattern, and the orientation it was built from.

Every assertion here has a source outside the code under test.

* The radial identity ``r = (camera constant) / d`` is the calibration relation the
  whole indexing chain rests on, so a simulated spot's distance from the beam is
  checked against its own d-spacing rather than against a stored number.
* The orientation matrix is checked by *using* it: multiplying a reflection's
  reciprocal vector by it must land on the pixel the pattern drew, and the zone
  axis must map onto the beam. That is what makes the matrix usable by the
  Kikuchi overlay, which consumes it and draws nothing the spots can check.
* The hexagonal basal net has ``|g|`` ratios of ``sqrt(3)`` between the first two
  rings whatever ``c/a`` is — a property of the hexagonal net, not of zirconium.
* Stating an orientation in Bunge angles and reading the angles back must return
  the same zone axis and the same roll, because the two ways of pointing the
  crystal are the same statement.
"""

from __future__ import annotations

import math
import re

import numpy as np
import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase
from pytex.app.registry import REGISTRY

pytest.importorskip("matplotlib", reason="the diffraction stack pulls in the plotting layer")

BASE_REQUEST = {
    "phase": {"builtin": "zr_hcp"},
    "orientation_source": "uvw",
    "zone_axis": [0, 0, 1],
    "in_plane_rotation_deg": 0.0,
    "camera_length_mm": 400.0,
    "beam_energy_kev": 200.0,
    "detector_px": 1024,
    "pixel_size_mm": 0.024,
    "show_kikuchi": False,
    "max_index": 6,
    "max_zone_index": 3,
}


def simulate(**overrides):
    request = dict(BASE_REQUEST)
    request.update(overrides)
    return REGISTRY.call("tem.simulate_saed", request)


def test_every_spot_sits_at_the_camera_constant_over_its_own_spacing() -> None:
    result = simulate()
    pattern = result["data"]["pattern"]
    centre = np.asarray(pattern["centre_px"], dtype=float)
    camera_constant = float(pattern["camera_constant_mm_angstrom"])
    pixel_size = float(pattern["pixel_size_mm"])

    assert pattern["spots"], "a basal zirconium pattern has reflections on a 1024 px detector"
    for spot in pattern["spots"]:
        radius_mm = math.hypot(spot["x"] - centre[0], spot["y"] - centre[1]) * pixel_size
        assert radius_mm == pytest.approx(camera_constant / spot["d_angstrom"], rel=1e-9, abs=1e-9)


def test_the_basal_net_puts_the_first_two_rings_a_root_three_apart() -> None:
    """A property of the hexagonal net rather than of the material.

    Down [0001] the reflections are the basal reciprocal net alone: the six
    prism {10-10} at one radius and the six {11-20} at sqrt(3) times it.
    """

    result = simulate()
    radii = sorted(
        {round(spot["g_inv_angstrom"], 6) for spot in result["data"]["pattern"]["spots"]}
    )
    assert len(radii) >= 2
    # The tolerance is a part in a hundred thousand rather than machine
    # precision because the pinned zirconium cell is a measured one, whose a and
    # b agree to the digits the CIF states rather than exactly.
    assert radii[1] / radii[0] == pytest.approx(math.sqrt(3.0), rel=1e-5)


def test_the_reported_orientation_places_every_spot_it_drew() -> None:
    """The matrix is checked by using it, not by comparing it with itself."""

    result = simulate(in_plane_rotation_deg=23.0)
    pattern = result["data"]["pattern"]
    matrix = np.asarray(result["data"]["orientation"]["crystal_to_pattern"], dtype=float).reshape(
        3, 3
    )
    assert float(np.linalg.det(matrix)) == pytest.approx(1.0, abs=1e-12)

    phase = builtin_phase("zr_hcp").to_phase()
    reciprocal = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=float)
    scale = float(pattern["camera_constant_mm_angstrom"]) / float(pattern["pixel_size_mm"])
    centre = np.asarray(pattern["centre_px"], dtype=float)

    for spot in pattern["spots"]:
        g_crystal = reciprocal @ np.asarray(spot["hkl"], dtype=float)
        projected = matrix @ g_crystal
        # In the zone, so the third component — the excitation error — is zero.
        assert float(projected[2]) == pytest.approx(0.0, abs=1e-9)
        placed = centre + projected[:2] * scale
        assert placed[0] == pytest.approx(spot["x"], abs=1e-6)
        assert placed[1] == pytest.approx(spot["y"], abs=1e-6)


def test_the_beam_direction_is_the_zone_axis_that_was_asked_for() -> None:
    result = simulate(zone_axis=[1, 1, 0])
    matrix = np.asarray(result["data"]["orientation"]["crystal_to_pattern"], dtype=float).reshape(
        3, 3
    )
    phase = builtin_phase("zr_hcp").to_phase()
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    axis_cartesian = direct @ np.array([1.0, 1.0, 0.0])
    axis_cartesian /= np.linalg.norm(axis_cartesian)
    assert np.allclose(matrix @ axis_cartesian, [0.0, 0.0, 1.0], atol=1e-12)


def test_hexagonal_three_index_input_is_reported_as_an_explicit_equivalence() -> None:
    """A notation conversion must never look like a changed physical direction."""

    result = simulate(zone_axis=[2, -1, 0])

    assert result["data"]["zone_axis"] == [2, -1, 0]
    assert result["data"]["zone_axis_three_index_label"] == "[2 -1 0]"
    assert result["data"]["zone_axis_label"] == "[5 -4 -1 0]"
    assert result["data"]["zone_axis_notation_label"] == "Zone axis [uvtw] (Miller-Bravais)"
    assert result["data"]["zone_axis_conversion_note"] == (
        "[2 -1 0] in the three-index hexagonal basis and [5 -4 -1 0] in four-index "
        "Miller-Bravais notation are the same crystal direction."
    )
    assert result["data"]["zone_axis_conversion_note"] in result["summary"]


def test_an_orientation_round_trips_through_its_own_euler_angles() -> None:
    """The two ways of pointing the crystal are one statement, so they agree."""

    from_axis = simulate(zone_axis=[1, 0, 0], in_plane_rotation_deg=17.0)
    phi1, phi, phi2 = from_axis["data"]["orientation"]["euler_bunge_deg"]

    from_euler = simulate(orientation_source="bunge", phi1_deg=phi1, Phi_deg=phi, phi2_deg=phi2)
    assert from_euler["data"]["zone_axis"] == from_axis["data"]["zone_axis"]
    assert from_euler["data"]["orientation"]["deviation_deg"] == pytest.approx(0.0, abs=1e-6)
    assert from_euler["data"]["orientation"]["in_plane_rotation_deg"] == pytest.approx(
        17.0, abs=1e-6
    )
    assert np.allclose(
        from_euler["data"]["orientation"]["crystal_to_pattern"],
        from_axis["data"]["orientation"]["crystal_to_pattern"],
        atol=1e-9,
    )


def test_an_orientation_between_two_zones_says_how_far_off_it_is() -> None:
    """An orientation that is not a zone axis must not be presented as one."""

    on_axis = simulate(
        phase={"builtin": "fe_bcc"},
        orientation_source="bunge",
        phi1_deg=0.0,
        Phi_deg=0.0,
        phi2_deg=0.0,
    )
    assert on_axis["data"]["zone_axis"] == [0, 0, 1]
    assert on_axis["data"]["orientation"]["deviation_deg"] == pytest.approx(0.0, abs=1e-9)

    tilted = simulate(
        phase={"builtin": "fe_bcc"},
        orientation_source="bunge",
        phi1_deg=0.0,
        Phi_deg=8.0,
        phi2_deg=0.0,
    )
    # Eight degrees off [001] about the specimen X axis, and [001] is still the
    # nearest low-index axis, so the deviation is the tilt itself.
    assert tilted["data"]["zone_axis"] == [0, 0, 1]
    assert tilted["data"]["orientation"]["deviation_deg"] == pytest.approx(8.0, abs=1e-6)
    assert "nearest zone axis" in tilted["summary"]
    assert any("nearest zone axis" in note for note in tilted["notes"])


def test_a_camera_length_that_puts_every_reflection_off_the_plate_is_refused() -> None:
    with pytest.raises(InvalidInputError) as excinfo:
        simulate(camera_length_mm=4000.0, detector_px=128)
    assert excinfo.value.details["field"] == "camera_length_mm"
    assert "camera length" in str(excinfo.value.hint)


def test_the_calibration_travels_with_the_pattern() -> None:
    """The solver beside it must never need a number retyped by hand."""

    result = simulate()
    calibration = result["data"]["calibration"]
    pattern = result["data"]["pattern"]
    assert calibration["units"] == "px"
    assert calibration["camera_constant_mm_angstrom"] == pytest.approx(
        pattern["camera_constant_mm_angstrom"]
    )
    assert calibration["pixel_size_mm"] == pytest.approx(pattern["pixel_size_mm"])
    assert calibration["scale_px_per_inv_angstrom"] == pytest.approx(
        pattern["camera_constant_mm_angstrom"] / pattern["pixel_size_mm"]
    )
    assert calibration["phase"] == {"builtin": "zr_hcp"}


def test_the_simulated_pattern_indexes_back_to_the_zone_it_was_built_from() -> None:
    """The forward and inverse operations must share their conventions."""

    result = simulate(zone_axis=[1, 1, 0], in_plane_rotation_deg=11.0)
    pattern = result["data"]["pattern"]
    centre = pattern["centre_px"]
    # Three non-collinear reflections, which is what the indexer needs.
    picks = []
    directions: list[np.ndarray] = []
    for spot in pattern["spots"]:
        offset = np.array([spot["x"] - centre[0], spot["y"] - centre[1]], dtype=float)
        norm = float(np.linalg.norm(offset))
        if norm < 1e-9:
            continue
        unit = offset / norm
        if any(abs(float(unit @ other)) > 0.999 for other in directions):
            continue
        directions.append(unit)
        picks.append({"x": spot["x"], "y": spot["y"]})
        if len(picks) == 3:
            break

    solved = REGISTRY.call(
        "tem.solve_pattern",
        {
            "phase": result["data"]["calibration"]["phase"],
            "picks": {"centre": list(centre), "spots": picks},
            "units": "px",
            "camera_constant_mm_angstrom": pattern["camera_constant_mm_angstrom"],
            "pixel_size_mm": pattern["pixel_size_mm"],
        },
    )
    assert solved["table"]["rows"], "the simulated pattern must index to something"

    # The answer is compared as a *family*, not as a triple. One pattern cannot
    # distinguish the members of a zone-axis family — they are the same picture —
    # so the solver is free to name any of them, and here it returns
    # [-1 2 -1 0] for a plate built down [1 1 -2 0]. Both are ⟨1 1 -2 0⟩.
    def family(label: str) -> list[int]:
        return sorted(abs(int(value)) for value in re.findall(r"-?\d+", label))

    assert family(solved["data"]["zone_axis_label"]) == family(result["data"]["zone_axis_label"])
    # Every spot indexed to the spacing it was drawn at.
    for row in solved["table"]["rows"]:
        assert abs(float(row["d_deviation_percent"])) < 1e-6


def test_double_diffraction_is_off_unless_asked_for() -> None:
    """The default pattern is purely kinematic, and its limits say so.

    Default-off is what makes the option safe to add: no existing pattern,
    gallery entry or pinned test changes meaning because the capability exists.
    """

    result = simulate()
    assert result["data"]["include_double_diffraction"] is False
    assert result["data"]["double_diffraction_count"] == 0
    assert all(row["origin"] == "Kinematic" for row in result["table"]["rows"])
    assert any("not modelled" in note for note in result["notes"])


def test_the_toggle_adds_the_forbidden_reflections_and_names_their_origin() -> None:
    """Silicon down [011] gains {200}, the textbook double-diffraction case."""

    result = simulate(
        phase={"builtin": "si_diamond"},
        zone_axis=[0, 1, -1],
        include_double_diffraction=True,
    )
    assert result["data"]["include_double_diffraction"] is True
    forbidden = [row for row in result["table"]["rows"] if row["origin"] != "Kinematic"]
    assert result["data"]["double_diffraction_count"] == len(forbidden)
    assert forbidden
    indices = {row["hkl"] for row in forbidden}
    assert "(200)" in indices
    origin = next(row["origin"] for row in forbidden if row["hkl"] == "(200)")
    assert origin == "(1 -1 -1) + (111)"


def test_a_pattern_that_includes_double_diffraction_does_not_claim_otherwise() -> None:
    """A stale limit is worse than no limit, so the note must be replaced.

    A reader who checks one stated limit and finds it false has no reason to
    trust the others, so the standing "double diffraction is not modelled" line
    must not survive onto a pattern that models it.
    """

    result = simulate(
        phase={"builtin": "si_diamond"},
        zone_axis=[0, 1, -1],
        include_double_diffraction=True,
    )
    assert not any("not modelled" in note for note in result["notes"])
    assert any("Double diffraction is included" in note for note in result["notes"])
    assert "Double diffraction is on" in result["summary"]


def test_a_hexagonal_row_names_its_parents_in_the_same_notation_as_itself() -> None:
    """One row must not mix four-index and three-index Miller notation.

    The repository's notation standard is not a formatting preference: a reader
    who sees (0001) explained by three-index parents has to work out which
    convention each half is in before the row means anything.
    """

    result = simulate(
        phase={"builtin": "ti_hcp"},
        zone_axis=[1, 1, 0],
        include_double_diffraction=True,
    )
    forbidden = [row for row in result["table"]["rows"] if row["origin"] != "Kinematic"]
    assert forbidden
    row = next(row for row in forbidden if row["hkl"] == "(0001)")
    # Four indices on both sides of the sum, as the phase is hexagonal.
    assert row["origin"] == "(-1 1 0 1) + (1 -1 0 0)"
    for other in forbidden:
        assert other["origin"].count(" + ") == 1


def test_a_centred_lattice_gains_nothing_and_the_summary_explains_why() -> None:
    """Turning it on for ferrite must not silently look like a broken feature."""

    result = simulate(
        phase={"builtin": "fe_bcc"},
        zone_axis=[0, 0, 1],
        include_double_diffraction=True,
    )
    assert result["data"]["double_diffraction_count"] == 0
    assert "adds nothing to this pattern" in result["summary"]
    assert "closed under addition" in result["summary"]


def test_the_marked_spots_reach_the_browser_with_their_origin() -> None:
    """The plate is drawn from this payload, so the marking has to be in it."""

    result = simulate(
        phase={"builtin": "si_diamond"},
        zone_axis=[0, 1, -1],
        include_double_diffraction=True,
    )
    spots = result["data"]["pattern"]["spots"]
    forbidden = [spot for spot in spots if spot["double_diffraction"]]
    assert forbidden
    for spot in forbidden:
        assert spot["double_diffraction_origin"]
        parents = spot["double_diffraction_parents"]
        assert parents is not None
        assert [sum(pair) for pair in zip(*parents, strict=True)] == spot["hkl"]
    for spot in spots:
        if not spot["double_diffraction"]:
            assert spot["double_diffraction_origin"] == ""
            assert spot["double_diffraction_parents"] is None
