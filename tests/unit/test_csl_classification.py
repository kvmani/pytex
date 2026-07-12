from __future__ import annotations

import math

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.ebsd import (
    CSLType,
    brandon_tolerance_deg,
    classify_misorientations,
)


def _cubic_operators() -> np.ndarray:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    return SymmetrySpec.from_point_group("m-3m", reference_frame=crystal).operators


def test_brandon_tolerance_scales_inversely_with_sqrt_sigma() -> None:
    assert brandon_tolerance_deg(1) == pytest.approx(15.0)
    assert brandon_tolerance_deg(3) == pytest.approx(15.0 / math.sqrt(3))
    assert brandon_tolerance_deg(9) == pytest.approx(5.0)


def test_classify_ideal_csl_misorientations() -> None:
    operators = _cubic_operators()
    sigma3 = CSLType(3, 60.0, (1, 1, 1)).matrix()
    sigma5 = CSLType(5, 36.86, (1, 0, 0)).matrix()
    (match3,) = classify_misorientations(sigma3, operators=operators)
    (match5,) = classify_misorientations(sigma5, operators=operators)
    assert match3 is not None and match3.sigma == 3
    assert match3.deviation_deg == pytest.approx(0.0, abs=1e-6)
    assert match5 is not None and match5.sigma == 5


def test_classify_batch_preserves_per_row_assignment() -> None:
    # A batch (n > 1) exercises the vectorized reduction; each row must be
    # classified independently (regression against a batch-axis reshape bug).
    operators = _cubic_operators()
    batch = np.stack(
        [
            CSLType(3, 60.0, (1, 1, 1)).matrix(),
            CSLType(5, 36.86, (1, 0, 0)).matrix(),
            CSLType(9, 38.94, (1, 1, 0)).matrix(),
        ]
    )
    matches = classify_misorientations(batch, operators=operators)
    assert [m.sigma for m in matches if m is not None] == [3, 5, 9]
    assert all(m is not None and m.deviation_deg == pytest.approx(0.0, abs=1e-6) for m in matches)


def test_low_angle_boundary_is_unclassified_by_default() -> None:
    operators = _cubic_operators()
    low_angle = CSLType(1, 3.0, (1, 0, 0)).matrix()
    (match,) = classify_misorientations(low_angle, operators=operators)
    assert match is None
    # opting into Sigma1 assigns the low-angle boundary
    (match1,) = classify_misorientations(low_angle, operators=operators, include_sigma1=True)
    assert match1 is not None and match1.sigma == 1


def _sigma3_twin_map() -> CrystalMap:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    axis = np.array([1.0, 1.0, 1.0])
    twin = Rotation.from_axis_angle(axis, math.radians(60.0))
    orientations = OrientationSet.from_quaternions(
        np.array([Rotation.identity().quaternion, twin.quaternion]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    return CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(1, 2),
        step_sizes=(1.0, 1.0),
    )


def test_boundary_network_classifies_sigma3_twin() -> None:
    crystal_map = _sigma3_twin_map()
    segmentation = crystal_map.segment_grains(max_misorientation_deg=15.0, symmetry_aware=True)
    assert len(segmentation.grains) == 2
    network = segmentation.boundary_network()
    matches = network.classify_csl()
    assert len(matches) == network.count
    assert any(match is not None and match.sigma == 3 for match in matches)
    assert network.csl_fraction(3) == pytest.approx(1.0)
    assert len(network.select_csl(3)) == network.count


def test_twin_merge_unions_sigma3_related_grains() -> None:
    crystal_map = _sigma3_twin_map()
    segmentation = crystal_map.segment_grains(max_misorientation_deg=15.0, symmetry_aware=True)
    assert len(segmentation.grains) == 2
    merged = segmentation.boundary_network().twin_merge()
    # the two Sigma3-related grains collapse into a single parent grain
    assert len(merged.grains) == 1
    assert merged.grains[0].size == 2


def test_classify_csl_rejects_non_cubic_maps() -> None:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal)
    tilted = Rotation.from_bunge_euler(30.0, 0.0, 0.0)
    orientations = OrientationSet.from_quaternions(
        np.array([Rotation.identity().quaternion, tilted.quaternion]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(1, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=15.0, symmetry_aware=True)
    network = segmentation.boundary_network()
    with pytest.raises(ValueError, match="cubic"):
        network.classify_csl()
