"""The Kikuchi map served to the crystal viewer, and the traces it is drawn from.

Every assertion has a source outside the code under test.

* **Where the zone axes are.** In a cubic crystal the angles from [001] are exact:
  0 degrees to itself, 45 to <011>, and arccos(1/sqrt(3)) = 54.7356 to <111>.
  Those come from the cubic metric, not from a stored result.
* **How wide a band is.** ``2 theta_B`` with ``sin theta_B = lambda / 2d``, checked
  against the spacing and wavelength the map itself reports, so the geometry and
  the numbers beside it cannot disagree.
* **Where the map's rim falls.** The stereographic radius of a direction at polar
  angle rho is ``tan(rho / 2)``, which is the definition of the projection.
* **What a trace is.** A band's centre line is the great circle perpendicular to
  the plane normal, so every sampled direction on it is perpendicular to that
  normal; an edge is the Kossel cone at the Bragg angle, so every direction on it
  makes exactly ``sin theta_B`` with it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase
from pytex.app.registry import REGISTRY
from pytex.diffraction.kikuchi_map import (
    compute_kikuchi_map,
    projected_trace_runs,
)


def call_map(**overrides):
    request = {"phase": {"builtin": "ni_fcc"}}
    request.update(overrides)
    return REGISTRY.call("crystal.kikuchi_map", request)


def test_the_map_finds_the_cubic_zone_axes_at_their_exact_angles() -> None:
    result = call_map()
    by_label = {axis["label"]: axis for axis in result["data"]["zone_axes"]}

    assert "[001]" in by_label
    assert by_label["[001]"]["polar_angle_deg"] == pytest.approx(0.0, abs=1e-9)
    # The centre of the projection is the centre of the map.
    assert by_label["[001]"]["x"] == pytest.approx(0.0, abs=1e-9)
    assert by_label["[001]"]["y"] == pytest.approx(0.0, abs=1e-9)

    hundred_and_ten = [
        axis
        for axis in result["data"]["zone_axes"]
        if sorted(abs(int(value)) for value in axis["uvw"]) == [0, 1, 1]
    ]
    assert hundred_and_ten, "the <011> axes are 45 degrees from [001] and inside a 60 degree map"
    for axis in hundred_and_ten:
        assert axis["polar_angle_deg"] == pytest.approx(45.0, abs=1e-6)

    triples = [
        axis
        for axis in result["data"]["zone_axes"]
        if sorted(abs(int(value)) for value in axis["uvw"]) == [1, 1, 1]
    ]
    for axis in triples:
        assert axis["polar_angle_deg"] == pytest.approx(
            math.degrees(math.acos(1.0 / math.sqrt(3.0))), abs=1e-6
        )


def test_every_band_is_as_wide_as_its_own_bragg_angle_says() -> None:
    result = call_map()
    wavelength = result["data"]["wavelength_angstrom"]
    assert result["data"]["bands"], "a nickel map has bands"
    for band in result["data"]["bands"]:
        expected = 2.0 * math.degrees(math.asin(wavelength / (2.0 * band["d_angstrom"])))
        assert band["width_deg"] == pytest.approx(expected, rel=1e-9)
        assert band["bragg_angle_deg"] == pytest.approx(expected / 2.0, rel=1e-9)


def test_the_map_rim_is_where_the_projection_puts_that_polar_angle() -> None:
    for polar in (30.0, 60.0, 90.0):
        result = call_map(max_polar_angle_deg=polar)
        assert result["data"]["boundary_radius"] == pytest.approx(
            math.tan(math.radians(polar) / 2.0), abs=1e-6
        )
        for axis in result["data"]["zone_axes"]:
            assert axis["polar_angle_deg"] <= polar + 1e-9


def test_the_traces_stay_inside_the_projection_disc() -> None:
    result = call_map(max_polar_angle_deg=90.0)
    for band in result["data"]["bands"]:
        runs = [*band["centre"], *band["edges"][0], *band["edges"][1]]
        assert runs, f"{band['label']} should cross the mapped region"
        for run in runs:
            assert len(run) >= 2
            for x, y in run:
                # A one-hemisphere stereographic projection lands inside the
                # unit circle by construction — a direction on the equator sits
                # exactly on it — so a point beyond it would mean the antipodal
                # fold was skipped. The slack is the four-decimal rounding the
                # payload is transported at, which can move a rim point by up to
                # about seven parts in a hundred thousand.
                assert math.hypot(x, y) <= 1.0 + 1e-4


def test_recentring_the_map_moves_the_axis_that_is_at_its_middle() -> None:
    centred = call_map(centre_direction="1 1 1")
    by_label = {axis["label"]: axis for axis in centred["data"]["zone_axes"]}
    assert by_label["[111]"]["polar_angle_deg"] == pytest.approx(0.0, abs=1e-9)
    assert math.hypot(by_label["[111]"]["x"], by_label["[111]"]["y"]) == pytest.approx(
        0.0, abs=1e-9
    )
    # [001] is now 54.74 degrees out rather than at the centre: the same crystal,
    # a different page of its atlas.
    assert by_label["[001]"]["polar_angle_deg"] == pytest.approx(
        math.degrees(math.acos(1.0 / math.sqrt(3.0))), abs=1e-6
    )


def test_the_view_matrix_puts_the_centre_direction_on_the_map_axis() -> None:
    result = call_map(centre_direction="1 1 1")
    view = np.asarray(result["data"]["view_matrix"], dtype=float).reshape(3, 3)
    assert float(np.linalg.det(view)) == pytest.approx(1.0, abs=1e-12)

    phase = builtin_phase("ni_fcc").to_phase()
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    centre = direct @ np.array([1.0, 1.0, 1.0])
    centre /= np.linalg.norm(centre)
    assert np.allclose(view @ centre, [0.0, 0.0, 1.0], atol=1e-12)


def test_a_direction_that_is_not_one_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        call_map(centre_direction="0 0 0")


def test_a_centre_parallel_to_the_horizontal_is_refused_with_a_usable_hint() -> None:
    with pytest.raises(InvalidInputError) as excinfo:
        call_map(centre_direction="1 0 0", horizontal_direction="1 0 0")
    assert "parallel" in str(excinfo.value.hint)


# ----------------------------------------------------- the library helpers


def test_a_band_centre_line_is_perpendicular_to_its_own_plane_normal() -> None:
    phase = builtin_phase("ni_fcc").to_phase()
    kikuchi_map = compute_kikuchi_map(phase, max_bands=4)
    for band in kikuchi_map.bands:
        directions = np.asarray(band.centre_directions(samples=91), dtype=float)
        assert directions.shape == (91, 3)
        assert np.allclose(directions @ np.asarray(band.normal_map, dtype=float), 0.0, atol=1e-12)


def test_a_band_edge_stands_at_the_bragg_angle_from_the_trace() -> None:
    phase = builtin_phase("ni_fcc").to_phase()
    kikuchi_map = compute_kikuchi_map(phase, max_bands=4)
    for band in kikuchi_map.bands:
        expected = math.sin(math.radians(band.bragg_angle_deg))
        for edge in band.edge_directions(samples=91):
            cosines = np.asarray(edge, dtype=float) @ np.asarray(band.normal_map, dtype=float)
            assert np.allclose(np.abs(cosines), expected, atol=1e-12)


def test_a_curve_crossing_the_equator_is_split_rather_than_chorded() -> None:
    """The defect this helper exists to prevent, stated as a test.

    A great circle through the poles crosses the equator twice. Folded onto one
    hemisphere and drawn as a single polyline it gains a chord straight across
    the middle of the figure, which is not a feature of the crystal.
    """

    angles = np.linspace(0.0, 2.0 * np.pi, 181)
    circle = np.column_stack([np.cos(angles), np.zeros_like(angles), np.sin(angles)])
    runs = projected_trace_runs(circle)
    assert len(runs) > 1
    for run in runs:
        assert run.shape[1] == 2
        assert run.shape[0] >= 2
        # Within a run the steps are small; a chord would be a single long jump.
        steps = np.linalg.norm(np.diff(run, axis=0), axis=1)
        assert float(steps.max()) < 0.2


def test_an_exaggerated_band_width_cannot_be_pushed_past_the_pole() -> None:
    phase = builtin_phase("ni_fcc").to_phase()
    band = compute_kikuchi_map(phase, max_bands=1).bands[0]
    # A scale large enough to send the edge beyond 90 degrees is clamped rather
    # than allowed to turn the band inside out.
    narrow, far = band.edge_directions(samples=31, width_scale=1000.0)
    for edge in (narrow, far):
        cosines = np.abs(np.asarray(edge, dtype=float) @ np.asarray(band.normal_map, dtype=float))
        assert float(cosines.max()) <= math.sin(math.radians(89.0)) + 1e-9
    with pytest.raises(ValueError):
        band.edge_directions(width_scale=0.0)
