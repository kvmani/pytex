"""The TEM solver, driven with patterns whose answer is known by construction.

The indexing tests build a pattern from a phase's own geometry — spots placed at
r = (camera constant)/d in the directions a chosen zone axis dictates — and then
require the solver to recover that zone axis and those indices. That is a round
trip rather than a pinned output: if the calibration, the geometry, or the
indexer disagree, the recovered answer stops matching the construction.

The tilt tests check the properties that must hold whatever the holder is, and
one that is easy to assume and false: reachability depends on the rotation about
the beam, which a single pattern cannot supply.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase

CAMERA_CONSTANT = 180.0
PIXEL_SIZE = 0.05
CENTRE = (512.0, 512.0)


def picks_for_cubic_001(phase_id: str, reflections: list[tuple[int, int]]) -> dict:
    """Place spots of a cubic [001] pattern at their exact pixel positions."""

    a = builtin_phase(phase_id).a
    spots = []
    for h, k in reflections:
        d = a / math.hypot(h, k)
        radius_mm = CAMERA_CONSTANT / d
        angle = math.atan2(k, h)
        spots.append(
            {
                "x": CENTRE[0] + radius_mm * math.cos(angle) / PIXEL_SIZE,
                "y": CENTRE[1] + radius_mm * math.sin(angle) / PIXEL_SIZE,
            }
        )
    return {"centre": list(CENTRE), "spots": spots}


def solve(phase_id: str = "ni_fcc", **overrides: object) -> dict:
    request: dict = {
        "phase": {"builtin": phase_id},
        "picks": picks_for_cubic_001(
            phase_id, [(2, 0), (0, 2), (-2, 0), (0, -2), (2, 2), (-2, 2), (2, -2), (-2, -2)]
        ),
        "units": "px",
        "camera_constant_mm_angstrom": CAMERA_CONSTANT,
        "pixel_size_mm": PIXEL_SIZE,
    }
    request.update(overrides)
    return REGISTRY.call("tem.solve_pattern", request)


class TestIndexing:
    def test_a_constructed_pattern_indexes_to_the_axis_it_was_built_from(self) -> None:
        result = solve()
        assert result["data"]["zone_axis_label"] == "[001]"
        assert result["data"]["phase_name"] == "Nickel (fcc)"

    def test_every_spot_is_indexed_with_zero_residual(self) -> None:
        result = solve()
        assert result["data"]["matched_fraction"] == pytest.approx(1.0)
        assert result["data"]["max_residual_inv_angstrom"] < 1e-9

    def test_observed_and_calculated_spacings_agree(self) -> None:
        result = solve()
        for row in result["table"]["rows"]:
            assert row["d_observed"] == pytest.approx(row["d_calculated"], rel=1e-9)

    def test_the_indexed_spacing_is_the_one_the_pattern_was_built_with(self) -> None:
        result = solve()
        a = builtin_phase("ni_fcc").a
        for row in result["table"]["rows"]:
            assert row["d_observed"] == pytest.approx(a / 2.0, rel=1e-6) or row[
                "d_observed"
            ] == pytest.approx(a / (2.0 * math.sqrt(2.0)), rel=1e-6)

    def test_the_answer_is_conclusive_for_one_candidate_phase(self) -> None:
        assert solve()["data"]["conclusive"] is True

    def test_a_wrong_camera_constant_fails_with_the_observed_spacings(self) -> None:
        # The most common real mistake: a camera constant from a different
        # camera length. The error must name the spacings so the user can see
        # at a glance that they are wrong for the phase.
        with pytest.raises(InvalidInputError) as excinfo:
            solve(camera_constant_mm_angstrom=25.0)
        assert "Å" in (excinfo.value.hint or "")

    def test_pixels_need_a_pixel_size(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            solve(pixel_size_mm=0.0)
        assert excinfo.value.details["field"] == "pixel_size_mm"

    def test_millimetres_need_a_camera_constant(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            solve(units="mm", camera_constant_mm_angstrom=0.0)
        assert excinfo.value.details["field"] == "camera_constant_mm_angstrom"

    def test_one_spot_is_not_enough(self) -> None:
        picks = picks_for_cubic_001("ni_fcc", [(2, 0)])
        with pytest.raises(InvalidInputError, match="At least two spots"):
            solve(picks=picks)

    def test_an_unmarked_beam_is_refused(self) -> None:
        picks = picks_for_cubic_001("ni_fcc", [(2, 0), (0, 2)])
        del picks["centre"]
        with pytest.raises(InvalidInputError, match="transmitted beam"):
            solve(picks=picks)

    def test_bcc_indexes_to_its_own_reflections(self) -> None:
        # 110-type spots: the first allowed family of a body-centred lattice.
        a = builtin_phase("fe_bcc").a
        spots = []
        for h, k in [(1, 1), (-1, 1), (1, -1), (-1, -1), (2, 0), (0, 2)]:
            d = a / math.hypot(h, k)
            radius_mm = CAMERA_CONSTANT / d
            angle = math.atan2(k, h)
            spots.append(
                {
                    "x": CENTRE[0] + radius_mm * math.cos(angle) / PIXEL_SIZE,
                    "y": CENTRE[1] + radius_mm * math.sin(angle) / PIXEL_SIZE,
                }
            )
        result = solve("fe_bcc", picks={"centre": list(CENTRE), "spots": spots})
        assert result["data"]["zone_axis_label"] == "[001]"


def plan(**overrides: object) -> dict:
    request: dict = {
        "phase": {"builtin": "austenite_fcc"},
        "current_zone_axis": [0, 0, 1],
        "target_zone_axis": [0, 1, 1],
    }
    request.update(overrides)
    return REGISTRY.call("tem.plan_tilt", request)


class TestTiltPlanning:
    def test_a_45_degree_move_is_out_of_reach_of_a_narrow_holder(self) -> None:
        result = plan(alpha_limit_deg=30.0, beta_limit_deg=20.0)
        assert result["data"]["exact"] is False
        assert result["data"]["reachable_orbit_size"] == 0

    def test_an_unreachable_target_still_reports_the_nearest_approach(self) -> None:
        result = plan(alpha_limit_deg=30.0, beta_limit_deg=20.0)
        assert result["table"]["rows"]
        assert result["data"]["nearest_approach_deg"] > 0.0

    def test_reachability_depends_on_the_rotation_about_the_beam(self) -> None:
        """The claim worth pinning, because it is easy to assume otherwise.

        One indexed pattern fixes the zone axis and leaves the roll about the
        beam free. Every roll gives identical spot positions — and yet with a
        ±40° holder, roll 0 puts every member of ⟨011⟩ out of reach while roll
        45° brings eight of them within it.
        """

        flat = plan(alpha_limit_deg=40.0, beta_limit_deg=40.0, beam_rotation_deg=0.0)
        rolled = plan(alpha_limit_deg=40.0, beta_limit_deg=40.0, beam_rotation_deg=45.0)
        assert flat["data"]["reachable_orbit_size"] == 0
        assert rolled["data"]["reachable_orbit_size"] > 0

    def test_the_orbit_size_is_the_symmetry_multiplicity(self) -> None:
        result = plan()
        # <011> has 12 members under m-3m with antipodal equivalence off, which
        # is the count the planner must consider.
        assert result["data"]["orbit_size"] == 12

    def test_a_reachable_move_reports_tilts_inside_the_envelope(self) -> None:
        result = plan(alpha_limit_deg=40.0, beta_limit_deg=40.0, beam_rotation_deg=45.0)
        assert result["data"]["exact"] is True
        for row in result["table"]["rows"]:
            if row["margin_deg"] >= 0.0:
                assert abs(row["alpha_deg"]) <= 40.0 + 1e-9
                assert abs(row["beta_deg"]) <= 40.0 + 1e-9

    def test_the_delta_is_the_move_from_the_current_position(self) -> None:
        result = plan(
            alpha_deg=5.0,
            beta_deg=-3.0,
            alpha_limit_deg=60.0,
            beta_limit_deg=60.0,
            beam_rotation_deg=45.0,
        )
        for row in result["table"]["rows"]:
            assert row["delta_alpha_deg"] == pytest.approx(row["alpha_deg"] - 5.0, abs=1e-6)
            assert row["delta_beta_deg"] == pytest.approx(row["beta_deg"] + 3.0, abs=1e-6)

    def test_travel_is_the_crystal_rotation_not_the_sum_of_the_tilts(self) -> None:
        result = plan(alpha_limit_deg=60.0, beta_limit_deg=60.0, beam_rotation_deg=45.0)
        row = result["table"]["rows"][0]
        assert row["travel_deg"] != pytest.approx(
            abs(row["delta_alpha_deg"]) + abs(row["delta_beta_deg"])
        )

    def test_uncertainty_propagates_into_the_predicted_tilts(self) -> None:
        precise = plan(
            alpha_limit_deg=60.0,
            beta_limit_deg=60.0,
            beam_rotation_deg=45.0,
            orientation_uncertainty_deg=0.1,
        )
        vague = plan(
            alpha_limit_deg=60.0,
            beta_limit_deg=60.0,
            beam_rotation_deg=45.0,
            orientation_uncertainty_deg=2.0,
        )
        assert (
            vague["table"]["rows"][0]["sigma_alpha_deg"]
            > (precise["table"]["rows"][0]["sigma_alpha_deg"])
        )

    def test_the_roll_caveat_is_always_stated(self) -> None:
        result = plan(alpha_limit_deg=60.0, beta_limit_deg=60.0, beam_rotation_deg=45.0)
        assert any("rotation about the beam" in note for note in result["notes"])

    def test_a_hexagonal_target_uses_four_index_labels(self) -> None:
        result = REGISTRY.call(
            "tem.plan_tilt",
            {
                "phase": {"builtin": "zr_hcp"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [1, 0, 0],
                "alpha_limit_deg": 60.0,
                "beta_limit_deg": 60.0,
            },
        )
        assert result["title"].startswith("[0001]")


class TestOrientationHelper:
    """The placement helper, checked against the geometry it claims."""

    def test_the_axis_ends_up_along_the_beam(self) -> None:
        from pytex.app.services.tem import _orientation_with_axis_on_beam
        from pytex.core.lattice import ZoneAxis
        from pytex.tem.stage import StagePosition, beam_direction_holder

        phase = builtin_phase("austenite_fcc").to_phase()
        axis = ZoneAxis(indices=np.array([1, 1, 0]), phase=phase)
        position = StagePosition(alpha_deg=12.0, beta_deg=-7.0)
        orientation = _orientation_with_axis_on_beam(phase, axis, position)
        matrix = np.asarray(orientation.rotation.as_matrix(), dtype=float)
        placed = matrix @ np.asarray(axis.unit_vector, dtype=float)
        beam = np.asarray(
            beam_direction_holder(position.alpha_deg, position.beta_deg), dtype=float
        ).reshape(3)
        assert np.allclose(placed, beam / np.linalg.norm(beam), atol=1e-9)

    def test_the_roll_keeps_the_axis_on_the_beam(self) -> None:
        from pytex.app.services.tem import _orientation_with_axis_on_beam
        from pytex.core.lattice import ZoneAxis
        from pytex.tem.stage import StagePosition, beam_direction_holder

        phase = builtin_phase("austenite_fcc").to_phase()
        axis = ZoneAxis(indices=np.array([0, 0, 1]), phase=phase)
        position = StagePosition(alpha_deg=0.0, beta_deg=0.0)
        orientation = _orientation_with_axis_on_beam(phase, axis, position, roll_deg=37.0)
        matrix = np.asarray(orientation.rotation.as_matrix(), dtype=float)
        placed = matrix @ np.asarray(axis.unit_vector, dtype=float)
        beam = np.asarray(beam_direction_holder(0.0, 0.0), dtype=float).reshape(3)
        # A roll is a rotation *about* the beam, so it must leave the axis on it.
        assert np.allclose(placed, beam / np.linalg.norm(beam), atol=1e-9)
