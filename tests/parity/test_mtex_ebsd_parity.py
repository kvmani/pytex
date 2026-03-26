from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    OrientationSet,
    ReferenceFrame,
    SymmetrySpec,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_json(path: str) -> dict:
    return json.loads((REPO_ROOT / path).read_text(encoding="utf-8"))


def _make_orientation_set(
    point_group: str, angles_deg: list[list[float]]
) -> tuple[OrientationSet, ReferenceFrame]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    orientations = OrientationSet.from_euler_angles(
        angles_deg,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        convention="bunge",
        degrees=True,
    )
    return orientations, specimen


def test_kam_cases_match_pinned_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/ebsd/kam_cases.json")
    case = fixture["case"]
    tolerance = float(fixture["tolerance"])
    orientations, specimen = _make_orientation_set(case["point_group"], case["angles_deg"])
    crystal_map = CrystalMap(
        coordinates=np.asarray(case["coordinates"], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=tuple(case["grid_shape"]),
        step_sizes=(1.0, 1.0),
    )
    assert_allclose(
        crystal_map.kernel_average_misorientation_deg(symmetry_aware=False, order=1),
        np.asarray(case["expected_order1_mean"], dtype=np.float64),
        atol=tolerance,
    )
    assert_allclose(
        crystal_map.kernel_average_misorientation_deg(symmetry_aware=False, order=2),
        np.asarray(case["expected_order2_mean"], dtype=np.float64),
        atol=tolerance,
    )
    assert_allclose(
        crystal_map.kernel_average_misorientation_deg(
            symmetry_aware=False,
            order=2,
            threshold_deg=3.0,
        ),
        np.asarray(case["expected_threshold3_mean"], dtype=np.float64),
        atol=tolerance,
    )
    assert_allclose(
        crystal_map.kernel_average_misorientation_deg(
            symmetry_aware=False,
            order=1,
            statistic="max",
        ),
        np.asarray(case["expected_order1_max"], dtype=np.float64),
        atol=tolerance,
    )


def test_grain_and_boundary_cases_match_pinned_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/ebsd/grains_cases.json")
    case = fixture["segmentation_case"]
    tolerance = float(fixture["tolerance"])
    orientations, specimen = _make_orientation_set(case["point_group"], case["angles_deg"])
    crystal_map = CrystalMap(
        coordinates=np.asarray(case["coordinates"], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=tuple(case["grid_shape"]),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=case["max_misorientation_deg"],
        symmetry_aware=False,
        connectivity=4,
    )
    assert_allclose(
        segmentation.label_grid, np.asarray(case["expected_labels"], dtype=np.int64), atol=tolerance
    )
    assert_allclose(
        segmentation.grod_map_deg(),
        np.asarray(case["expected_grod_deg"], dtype=np.float64),
        atol=tolerance,
    )
    network = segmentation.boundary_network(high_angle_threshold_deg=15.0)
    assert network.count == case["expected_boundary_count"]
    assert network.high_angle_count == case["expected_high_angle_count"]
    assert_allclose(network.total_length, case["expected_boundary_total_length"], atol=tolerance)


def test_cleanup_case_matches_pinned_values() -> None:
    fixture = _load_json("fixtures/mtex_parity/ebsd/grains_cases.json")
    case = fixture["cleanup_case"]
    tolerance = float(fixture["tolerance"])
    orientations, specimen = _make_orientation_set(case["point_group"], case["angles_deg"])
    crystal_map = CrystalMap(
        coordinates=np.asarray(case["coordinates"], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=tuple(case["grid_shape"]),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=case["max_misorientation_deg"],
        symmetry_aware=False,
        connectivity=4,
    )
    merged = segmentation.merge_small_grains(min_size=2)
    assert_allclose(
        merged.label_grid, np.asarray(case["expected_merged_label"], dtype=np.int64), atol=tolerance
    )
