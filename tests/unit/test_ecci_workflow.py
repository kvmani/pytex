"""The ECCI workflow is tested as geometry: every claim here is checkable by
hand from the stage convention stated in ``pytex.app.services.ecci``, never
from a value captured off a run of the code under test.

The central claim - a solved tilt/rotation actually reaches the two-beam
condition it was solved for - is checked the same way the notebook
demonstrates it: solve, then re-simulate at the solved state, and show the
target direction lands on the beam and the on-axis reflections' excitation
errors collapse to (numerically) zero.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase
from pytex.app.services.ecci import (
    _beam_direction_crystal,
    _beam_direction_specimen,
    _nearest_zone_axis,
    _solve_stage_for_direction,
    _stage_branches,
    _stage_to_lab_matrix,
)
from pytex.core.frame_catalog import SPECIMEN_FRAME
from pytex.core.orientation import Orientation


def _base_request(**overrides: object) -> dict:
    request: dict[str, object] = {
        "phase": {"builtin": "ni_fcc"},
        "phi1_deg": 0.0,
        "Phi_deg": 0.0,
        "phi2_deg": 0.0,
        "stage_tilt_deg": 70.0,
        "stage_rotation_deg": 0.0,
        "target_zone_axis": [1, 1, 1],
    }
    request.update(overrides)
    return request


# --------------------------------------------------------------------------
# Pure geometry: the stage matrix and the closed-form solver, independent of
# the registered operations.
# --------------------------------------------------------------------------


def test_stage_matrix_matches_for_ebsd_at_zero_rotation() -> None:
    """At rotation = 0 the extended stage matrix must reduce to for_ebsd's own.

    This is the contract the whole module leans on: an EBSD-measured
    orientation and geometry, at the stage's own zero-rotation reading, must
    pass through this module unchanged.
    """

    from pytex.diffraction.models import DiffractionGeometry

    for tilt in (0.0, 20.0, 54.7, 70.0, 89.0):
        geometry = DiffractionGeometry.for_ebsd(sample_tilt_deg=tilt)
        assert_allclose(
            _stage_to_lab_matrix(tilt, 0.0), geometry.specimen_to_lab_matrix, atol=1e-12
        )


def test_stage_matrix_is_a_proper_rotation() -> None:
    rng = np.random.default_rng(0)
    for _ in range(20):
        tilt = float(rng.uniform(0.0, 89.0))
        rotation = float(rng.uniform(-180.0, 180.0))
        matrix = _stage_to_lab_matrix(tilt, rotation)
        assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-10)
        assert math.isclose(float(np.linalg.det(matrix)), 1.0, abs_tol=1e-10)


def test_stage_branches_place_the_direction_exactly_on_beam() -> None:
    """Every branch _stage_branches returns must satisfy its own construction.

    Rx(180 - tilt) @ Rz(rotation) @ v must equal exactly (0, 0, 1) for a
    branch's own (tilt, rotation) - this is what the closed form claims to
    solve, checked directly against random unit vectors rather than against
    any value the module itself produced elsewhere.
    """

    rng = np.random.default_rng(1)
    for _ in range(30):
        v = rng.normal(size=3)
        v = v / np.linalg.norm(v)
        for tilt_deg, rotation_deg in _stage_branches(v):
            matrix = _stage_to_lab_matrix(tilt_deg, rotation_deg)
            achieved = matrix @ v
            assert_allclose(achieved, [0.0, 0.0, 1.0], atol=1e-9)


def test_beam_direction_specimen_is_minus_z_at_zero_tilt_zero_rotation() -> None:
    """for_ebsd's own documented convention: an untilted specimen normal faces
    the beam, i.e. the beam lies along specimen -z, not +z."""

    beam = _beam_direction_specimen(0.0, 0.0)
    assert_allclose(beam, [0.0, 0.0, -1.0], atol=1e-12)


def test_beam_direction_specimen_is_plus_z_at_180_tilt() -> None:
    """At the (unphysical but mathematically valid) 180 deg tilt the stage
    matrix is the identity, so the beam lies along specimen +z."""

    beam = _beam_direction_specimen(180.0, 0.0)
    assert_allclose(beam, [0.0, 0.0, 1.0], atol=1e-10)


def test_solve_stage_for_direction_is_forward_consistent() -> None:
    """Every returned solution's own reported residual must be (numerically)
    zero: a solution is only ever reported after being forward-validated."""

    phase = builtin_phase("ni_fcc").to_phase()
    orientation = Orientation.from_euler(
        12.0, 34.0, 56.0, degrees=True, specimen_frame=SPECIMEN_FRAME, phase=phase
    )
    direction = np.array([1.0, 1.0, 1.0])
    solutions = _solve_stage_for_direction(
        orientation, direction, current_tilt_deg=70.0, current_rotation_deg=0.0, allow_reverse=True
    )
    assert solutions
    for solution in solutions:
        assert solution["residual_deg"] < 1e-6
        assert 0.0 <= solution["tilt_deg"] < 89.9
    # Shortest travel first.
    travels = [solution["travel_deg"] for solution in solutions]
    assert travels == sorted(travels)


def test_nearest_zone_axis_of_an_exact_pole_is_itself() -> None:
    phase = builtin_phase("fe_bcc").to_phase()
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    exact = direct @ np.array([1.0, 0.0, 1.0])
    exact = exact / np.linalg.norm(exact)
    indices, deviation_deg = _nearest_zone_axis(phase, exact, max_index=4)
    assert deviation_deg < 1e-6
    assert sorted(abs(v) for v in indices) == [0, 1, 1]


# --------------------------------------------------------------------------
# The registered operations, exercised through REGISTRY.call, as
# test_app_tem_kikuchi.py does for the analogous TEM operations.
# --------------------------------------------------------------------------


def test_solve_workflow_forward_validates_to_zero_residual() -> None:
    result = REGISTRY.call("ecci.solve_workflow", _base_request())
    solution = result["data"]["solution"]
    assert solution["residual_deg"] == pytest.approx(0.0, abs=1e-6)


def test_solve_workflow_then_resimulate_reaches_the_target() -> None:
    """The end-to-end claim: solving, then re-simulating at the solved state,
    actually brings the target direction onto the beam and collapses the
    on-axis excitation errors - the same demonstration the tutorial notebook
    makes."""

    request = _base_request()
    solved = REGISTRY.call("ecci.solve_workflow", request)
    solution = solved["data"]["solution"]

    before = REGISTRY.call("ecci.resimulate", request)
    after = REGISTRY.call(
        "ecci.resimulate",
        {
            **request,
            "stage_tilt_deg": solution["tilt_deg"],
            "stage_rotation_deg": solution["rotation_deg"],
        },
    )

    assert before["data"]["target"]["angle_from_beam_deg"] > 1.0
    assert after["data"]["target"]["angle_from_beam_deg"] == pytest.approx(0.0, abs=1e-6)

    # The on-axis pattern's reflections belong to the zone of the (now
    # on-beam) target direction, so their excitation errors - dot(g, beam) -
    # must all be driven to (numerically) zero: the two-beam/zone-axis
    # condition, demonstrated rather than asserted by construction alone.
    assert before["data"]["on_axis"]["max_abs_excitation_error_inv_angstrom"] > 1e-3
    assert after["data"]["on_axis"]["max_abs_excitation_error_inv_angstrom"] < 1e-6

    assert after["data"]["proximity"]["deviation_deg"] == pytest.approx(0.0, abs=1e-6)


def test_resimulate_reports_the_target_direction_it_was_asked_for() -> None:
    """A direct hand check: the target direction and the actual beam direction
    are compared through the module's own (independently written) helper, and
    must agree with the reported angle."""

    request = _base_request(stage_tilt_deg=54.735610317245346, stage_rotation_deg=-135.0)
    result = REGISTRY.call("ecci.resimulate", request)

    phase = builtin_phase("ni_fcc").to_phase()
    orientation = Orientation.from_euler(
        0.0, 0.0, 0.0, degrees=True, specimen_frame=SPECIMEN_FRAME, phase=phase
    )
    beam_crystal = _beam_direction_crystal(orientation, 54.735610317245346, -135.0)
    beam_crystal = beam_crystal / np.linalg.norm(beam_crystal)
    target = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    expected_angle = math.degrees(math.acos(min(1.0, abs(float(np.dot(beam_crystal, target))))))

    assert result["data"]["target"]["angle_from_beam_deg"] == pytest.approx(
        expected_angle, abs=1e-6
    )


def test_solve_workflow_kikuchi_pattern_matches_ebsd_pattern_operation() -> None:
    """At stage_rotation_deg = 0 the embedded Kikuchi payload must be the same
    geometry ebsd.simulate_kikuchi_pattern reports, since both call the same
    underlying simulator with the same inputs."""

    request = _base_request(stage_rotation_deg=0.0, target_zone_axis=[0, 0, 1])
    ecci_result = REGISTRY.call("ecci.solve_workflow", request)
    ebsd_result = REGISTRY.call(
        "ebsd.simulate_kikuchi_pattern",
        {
            "phase": {"builtin": "ni_fcc"},
            "phi1_deg": 0.0,
            "Phi_deg": 0.0,
            "phi2_deg": 0.0,
            "sample_tilt_deg": 70.0,
            "detector_distance": 0.65,
            "beam_energy_kev": 20.0,
            "max_bands": 24,
            "max_index": 4,
            "zone_axis_max_index": 3,
        },
    )
    ecci_bands = {tuple(band["hkl"]) for band in ecci_result["data"]["kikuchi"]["bands"]}
    ebsd_bands = {tuple(band["hkl"]) for band in ebsd_result["data"]["bands"]}
    assert ecci_bands == ebsd_bands
    assert ecci_result["data"]["kikuchi"]["wavelength_angstrom"] == pytest.approx(
        ebsd_result["data"]["wavelength_angstrom"]
    )


def test_on_axis_pattern_excitation_error_zero_at_exact_zone() -> None:
    """When the beam is placed exactly along a low-index pole (not merely
    close to one), every reflection of that zone must show zero excitation
    error - the defining property of the quantity, checked at a point where
    the answer is known by construction rather than by solving anything."""

    # At tilt = 0, rotation = 0 the beam lies along specimen -z (for_ebsd's own
    # documented convention), which for the identity orientation is crystal
    # -z, i.e. exactly the [0 0 -1]/[001] pole.
    request = _base_request(
        phi1_deg=0.0, Phi_deg=0.0, phi2_deg=0.0, stage_tilt_deg=0.0, stage_rotation_deg=0.0
    )
    result = REGISTRY.call("ecci.resimulate", request)
    assert sorted(abs(v) for v in result["data"]["on_axis"]["nearest_zone_axis"]) == [0, 0, 1]
    assert result["data"]["on_axis"]["max_abs_excitation_error_inv_angstrom"] < 1e-9
    for spot in result["data"]["on_axis"]["spots"]:
        assert abs(spot["excitation_error_inv_angstrom"]) < 1e-9


# --------------------------------------------------------------------------
# Input validation.
# --------------------------------------------------------------------------


def test_unknown_phase_is_rejected() -> None:
    with pytest.raises(InvalidInputError):
        REGISTRY.call("ecci.solve_workflow", _base_request(phase={"builtin": "not_a_phase"}))


def test_non_numeric_euler_angle_is_rejected() -> None:
    with pytest.raises(InvalidInputError):
        REGISTRY.call("ecci.solve_workflow", _base_request(phi1_deg="not a number"))


def test_all_zero_target_direction_is_rejected() -> None:
    with pytest.raises(InvalidInputError):
        REGISTRY.call("ecci.solve_workflow", _base_request(target_zone_axis=[0, 0, 0]))


def test_missing_phase_is_rejected() -> None:
    request = _base_request()
    del request["phase"]
    with pytest.raises(InvalidInputError):
        REGISTRY.call("ecci.solve_workflow", request)


def test_unknown_parameter_is_rejected() -> None:
    with pytest.raises(InvalidInputError):
        REGISTRY.call("ecci.solve_workflow", _base_request(bogus_field=1.0))


def test_resimulate_matches_solve_workflow_current_state() -> None:
    """solve_workflow's own current-state view and a direct resimulate call at
    the same stage state must describe the same geometry."""

    request = _base_request()
    solved = REGISTRY.call("ecci.solve_workflow", request)
    direct = REGISTRY.call("ecci.resimulate", request)
    assert solved["data"]["proximity"] == direct["data"]["proximity"]
    assert (
        solved["data"]["kikuchi"]["wavelength_angstrom"]
        == direct["data"]["kikuchi"]["wavelength_angstrom"]
    )
