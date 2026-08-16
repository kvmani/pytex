"""The stereogram is a map, so it is tested as one: does it put poles where they are.

Every assertion here is against a number that is known independently of this
code — an interzonal angle fixed by cubic geometry, the ``tan(rho/2)`` law of the
stereographic projection, or the closed form in :mod:`pytex.tem.navigation` that
the panel's stage readings must reproduce.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.services.tem import _primitive_directions, _stage_angles_for_holder
from pytex.core.sphere import project_directions
from pytex.tem.navigation import solve_tilts_for_direction
from pytex.tem.stage import beam_direction_holder


def stereogram(**overrides: object) -> dict:
    request: dict[str, object] = {
        "phase": {"builtin": "al_fcc"},
        "zone_axis": "0 0 1",
    }
    request.update(overrides)
    return REGISTRY.call("tem.stereogram", request)


def axis_named(result: dict, label: str) -> dict:
    matches = [entry for entry in result["data"]["axes"] if entry["label"] == label]
    assert matches, f"the stereogram plotted no pole labelled {label}"
    assert len(matches) == 1, f"{label} was plotted {len(matches)} times"
    return matches[0]


def axis_indexed(result: dict, indices: list[int]) -> dict:
    """Find a pole by its indices rather than its label.

    Hexagonal directions are labelled in four-index Weber notation, so ``[100]``
    is written ``[2 1̄ 1̄ 0]`` and cannot be looked up by the triple that was
    asked for.
    """

    reversed_indices = [-value for value in indices]
    matches = [
        entry for entry in result["data"]["axes"] if entry["indices"] in (indices, reversed_indices)
    ]
    assert matches, f"the stereogram plotted no pole with indices {indices}"
    assert len(matches) == 1
    return matches[0]


def test_stage_angles_match_the_navigation_closed_form() -> None:
    """The label beside a pole must be the angle the planner would give it.

    The panel inlines the principal branch of `solve_tilts_for_direction` so it
    can label a thousand poles without building a solution object for each. That
    is only safe while the two agree, which is what this checks — including the
    gimbal case along the beta axis, where the branch is degenerate.
    """

    rng = np.random.default_rng(11)
    directions = rng.normal(size=(400, 3))
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    for direction in directions:
        expected = solve_tilts_for_direction(direction)[0]
        assert_allclose(_stage_angles_for_holder(direction), expected, atol=1e-9)

    for pole in ([0.0, 1.0, 0.0], [0.0, -1.0, 0.0]):
        assert_allclose(
            _stage_angles_for_holder(np.asarray(pole)),
            solve_tilts_for_direction(pole)[0],
            atol=1e-9,
        )


def test_the_stage_reading_really_puts_the_pole_on_the_beam() -> None:
    """A round trip through the stage model, not a comparison with itself."""

    result = stereogram()
    for entry in result["data"]["axes"][:40]:
        beam = beam_direction_holder(entry["alpha_deg"], entry["beta_deg"])
        # The pole is recovered up to sense: the pattern cannot tell the two apart.
        recovered = np.asarray(beam, dtype=float)
        pole = np.asarray(
            [
                2.0 * entry["x"] / (1.0 + entry["x"] ** 2 + entry["y"] ** 2),
                2.0 * entry["y"] / (1.0 + entry["x"] ** 2 + entry["y"] ** 2),
                (1.0 - entry["x"] ** 2 - entry["y"] ** 2)
                / (1.0 + entry["x"] ** 2 + entry["y"] ** 2),
            ]
        )
        assert abs(float(np.dot(recovered, pole))) == pytest.approx(1.0, abs=1e-9)


def test_projection_follows_the_tangent_half_angle_law() -> None:
    """``r = tan(rho / 2)``: the defining property of a stereographic net."""

    result = stereogram()
    for entry in result["data"]["axes"]:
        radius = math.hypot(entry["x"], entry["y"])
        rho = math.radians(entry["angle_from_beam_deg"])
        assert radius == pytest.approx(math.tan(rho / 2.0), abs=1e-9)
    # And the centre is the axis on the beam.
    centre = axis_named(result, "[001]")
    assert math.hypot(centre["x"], centre["y"]) == pytest.approx(0.0, abs=1e-12)


def test_cubic_interzonal_angles_are_the_textbook_ones() -> None:
    result = stereogram()
    assert axis_named(result, "[011]")["angle_from_beam_deg"] == pytest.approx(45.0, abs=1e-9)
    assert axis_named(result, "[111]")["angle_from_beam_deg"] == pytest.approx(
        54.735610317245346, abs=1e-9
    )
    # atan(1/2) and atan(sqrt(2)/2): the two standard first stops out of [001].
    assert axis_named(result, "[012]")["angle_from_beam_deg"] == pytest.approx(
        math.degrees(math.atan(0.5)), abs=1e-9
    )
    assert axis_named(result, "[112]")["angle_from_beam_deg"] == pytest.approx(
        math.degrees(math.atan(math.sqrt(2.0) / 2.0)), abs=1e-9
    )


def test_poles_are_primitive_and_appear_once() -> None:
    """A stereogram plots poles, not integer triples.

    ``[002]`` is ``[001]``; ``[uvw]`` and its reverse are the same axis once the
    hemisphere is folded. Both used to appear as separate markers, which put two
    labels on one point.
    """

    directions = _primitive_directions(3)
    magnitudes = np.abs(directions)
    divisor = np.gcd(np.gcd(magnitudes[:, 0], magnitudes[:, 1]), magnitudes[:, 2])
    assert np.all(divisor == 1)
    keys = {tuple(int(value) for value in row) for row in directions}
    for row in directions:
        assert tuple(int(-value) for value in row) not in keys

    labels = [entry["label"] for entry in stereogram()["data"]["axes"]]
    assert len(labels) == len(set(labels))
    assert "[002]" not in labels


def test_the_envelope_outline_encloses_exactly_the_reachable_poles() -> None:
    """The drawn region and the reachability flag must be the same claim."""

    result = stereogram(alpha_limit_deg=30.0, beta_limit_deg=20.0)
    for entry in result["data"]["axes"]:
        inside = abs(entry["alpha_deg"]) <= 30.0 + 1e-9 and abs(entry["beta_deg"]) <= 20.0 + 1e-9
        assert entry["reachable"] is inside, entry["label"]

    boundary = np.asarray(result["data"]["envelope"]["boundary"], dtype=float)
    assert boundary.shape[1] == 2
    # The outline is a closed loop through the four corners of the tilt range.
    corner = project_directions(beam_direction_holder(30.0, 20.0), method="stereographic")[0]
    assert np.min(np.linalg.norm(boundary - corner, axis=1)) < 1e-9


def test_the_route_runs_from_the_beam_to_the_target_through_low_index_zones() -> None:
    result = stereogram(target_zone_axis="0 1 1")
    data = result["data"]
    target = data["target"]
    assert target["angle_from_beam_deg"] == pytest.approx(45.0, abs=1e-9)
    assert target["alpha_deg"] == pytest.approx(45.0, abs=1e-9)
    assert target["reachable"] is False  # 45 degrees of alpha, on a +/-30 holder

    points = np.asarray(data["path"]["points"], dtype=float)
    assert points.shape[0] > 10
    assert_allclose(points[0], [data["beam"]["x"], data["beam"]["y"]], atol=1e-9)
    assert_allclose(points[-1], [target["x"], target["y"]], atol=1e-9)

    # [012] is the standard intermediate between [001] and [011], and it lies on
    # the drawn route rather than merely being named beside it.
    waypoints = data["path"]["waypoints"]
    assert [waypoint["label"] for waypoint in waypoints] == ["[012]"]
    offsets = np.linalg.norm(points - np.asarray([waypoints[0]["x"], waypoints[0]["y"]]), axis=1)
    assert float(np.min(offsets)) < 0.02


def test_no_target_draws_no_route() -> None:
    data = stereogram(target_zone_axis="0 0 0")["data"]
    assert data["target"] is None
    assert data["path"] is None


def test_rolling_the_crystal_turns_the_drawing_and_leaves_the_angles_alone() -> None:
    """The roll is exactly the thing one indexed pattern does not determine.

    It must move where a pole is drawn — that is how a move divides between
    alpha and beta — and must not move how far away the pole is, which is fixed
    by the lattice.
    """

    straight = stereogram()
    rolled = stereogram(beam_rotation_deg=90.0)
    for entry in straight["data"]["axes"]:
        other = axis_named(rolled, entry["label"])
        assert other["angle_from_beam_deg"] == pytest.approx(entry["angle_from_beam_deg"], abs=1e-9)
    moved = [
        entry
        for entry in straight["data"]["axes"]
        if math.hypot(entry["x"], entry["y"]) > 0.05
        and not math.isclose(axis_named(rolled, entry["label"])["x"], entry["x"], abs_tol=1e-6)
    ]
    assert moved, "a 90 degree roll must move the poles on the drawing"


def test_the_stage_position_moves_the_beam_off_the_centre() -> None:
    """The centre is the holder's zero, not wherever the crystal happens to be."""

    data = stereogram(alpha_deg=15.0, beta_deg=-10.0)["data"]
    assert data["beam"]["alpha_deg"] == 15.0
    assert math.hypot(data["beam"]["x"], data["beam"]["y"]) > 0.1
    # The axis said to be on the beam is drawn at the beam, wherever that is.
    on_beam = axis_named(stereogram(alpha_deg=15.0, beta_deg=-10.0), "[001]")
    assert on_beam["x"] == pytest.approx(data["beam"]["x"], abs=1e-9)
    assert on_beam["y"] == pytest.approx(data["beam"]["y"], abs=1e-9)
    # A thousandth of a degree, not machine epsilon: the angle is an arccos near
    # 1, where the square root turns 1e-16 of cosine error into ~1e-7 of a
    # degree. The atlas carries the same caveat for the same reason.
    assert on_beam["angle_from_beam_deg"] == pytest.approx(0.0, abs=1e-3)


def test_a_hexagonal_phase_projects_its_own_angles() -> None:
    """Not a cubic special case: zirconium's c/a decides these angles."""

    result = stereogram(phase={"builtin": "zr_hcp"}, zone_axis="0 0 1", max_index=2)
    basal = axis_indexed(result, [0, 0, 1])
    assert basal["label"] == "[0001]"  # four-index Weber notation, as hcp is written
    assert basal["angle_from_beam_deg"] == pytest.approx(0.0, abs=1e-9)
    # [100] and [110] both lie in the basal plane, 90 degrees from [0001], and
    # 60 degrees apart from each other — the hexagonal net, not a cubic one.
    prism = axis_indexed(result, [1, 0, 0])
    assert prism["angle_from_beam_deg"] == pytest.approx(90.0, abs=1e-9)
    assert math.hypot(prism["x"], prism["y"]) == pytest.approx(1.0, abs=1e-9)
    other = axis_indexed(result, [1, 1, 0])
    assert other["angle_from_beam_deg"] == pytest.approx(90.0, abs=1e-9)
    separation = math.degrees(
        math.acos(
            (prism["x"] * other["x"] + prism["y"] * other["y"])
            / (math.hypot(prism["x"], prism["y"]) * math.hypot(other["x"], other["y"]))
        )
    )
    assert separation == pytest.approx(60.0, abs=1e-6)


def test_an_impossible_index_limit_is_refused_rather_than_absorbed() -> None:
    with pytest.raises(InvalidInputError):
        stereogram(max_index=0)
