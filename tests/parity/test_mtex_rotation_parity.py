from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pytex.core.orientation import quaternion_to_matrix
from pytex import Rotation

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_json(path: str) -> dict:
    return json.loads((REPO_ROOT / path).read_text(encoding="utf-8"))


def test_euler_quaternion_cases_match_pinned_parity_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/rotation/euler_quaternion_cases.json")
    tolerance = float(fixture["tolerance"])
    for case in fixture["cases"]:
        rotation = Rotation.from_euler(
            *case["angles"],
            convention=case["convention"],
            degrees=bool(case["degrees"]),
        )
        expected = np.asarray(case["quaternion_wxyz"], dtype=np.float64)
        assert_allclose(rotation.canonicalized().quaternion, expected, atol=tolerance)
        recovered = Rotation(rotation.quaternion).to_euler(
            convention=case["convention"],
            degrees=bool(case["degrees"]),
        )
        round_trip = Rotation.from_euler(
            *recovered,
            convention=case["convention"],
            degrees=bool(case["degrees"]),
        )
        assert_allclose(round_trip.as_matrix(), rotation.as_matrix(), atol=tolerance)


def test_axis_angle_cases_match_pinned_quaternion_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/rotation/axis_angle_cases.json")
    tolerance = float(fixture["tolerance"])
    for case in fixture["cases"]:
        rotation = Rotation.from_axis_angle(case["axis"], np.deg2rad(case["angle_deg"]))
        expected = np.asarray(case["quaternion_wxyz"], dtype=np.float64)
        assert_allclose(rotation.canonicalized().quaternion, expected, atol=tolerance)


def test_quaternion_matrix_round_trip_preserves_rotation_matrix() -> None:
    fixture = _load_json("fixtures/mtex_parity/rotation/euler_quaternion_cases.json")
    tolerance = float(fixture["tolerance"])
    for case in fixture["cases"]:
        quaternion = np.asarray(case["quaternion_wxyz"], dtype=np.float64)
        matrix = quaternion_to_matrix(quaternion)
        recovered = Rotation.from_matrix(matrix)
        assert_allclose(recovered.as_matrix(), matrix, atol=tolerance)
