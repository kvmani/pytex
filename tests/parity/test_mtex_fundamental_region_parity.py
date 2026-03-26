from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pytex import (
    FrameDomain,
    Handedness,
    Orientation,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_json(path: str) -> dict:
    return json.loads((REPO_ROOT / path).read_text(encoding="utf-8"))


def _make_frames() -> tuple[ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    return crystal, specimen


def test_vector_reduction_cases_match_pinned_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/fundamental_region/vector_cases.json")
    tolerance = float(fixture["tolerance"])
    crystal, _ = _make_frames()
    for case in fixture["cases"]:
        symmetry = SymmetrySpec.from_point_group(case["point_group"], reference_frame=crystal)
        reduced = symmetry.reduce_vector_to_fundamental_sector(
            case["input_vector"],
            antipodal=bool(case["antipodal"]),
        )
        expected = np.asarray(case["expected_vector"], dtype=np.float64)
        assert_allclose(reduced, expected, atol=tolerance)
        assert symmetry.vector_in_fundamental_sector(reduced, antipodal=bool(case["antipodal"]))


def test_orientation_projection_cases_match_pinned_projected_angles() -> None:
    fixture = _load_json("fixtures/mtex_parity/fundamental_region/orientation_cases.json")
    tolerance = float(fixture["tolerance"])
    crystal, specimen = _make_frames()
    identity = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=None,
    )
    for case in fixture["cases"]:
        symmetry = SymmetrySpec.from_point_group(case["point_group"], reference_frame=crystal)
        orientation = Orientation(
            rotation=Rotation.from_euler(
                *case["angles_deg"],
                convention=case["euler_convention"],
                degrees=True,
            ),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        )
        projected = orientation.project_to_fundamental_region(reference_orientation=identity)
        assert_allclose(
            projected.rotation.angle_deg, case["expected_projected_angle_deg"], atol=tolerance
        )
