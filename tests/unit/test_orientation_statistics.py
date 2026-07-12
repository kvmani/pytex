from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    FrameDomain,
    Handedness,
    Orientation,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def make_frames() -> tuple[ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    return crystal, specimen


def cubic_symmetry(crystal: ReferenceFrame) -> SymmetrySpec:
    return SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)


def make_set(
    euler_deg: list[tuple[float, float, float]],
    *,
    symmetry: SymmetrySpec | None = None,
) -> OrientationSet:
    crystal, specimen = make_frames()
    resolved_symmetry = symmetry
    return OrientationSet.from_euler_angles(
        np.asarray(euler_deg, dtype=np.float64),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=resolved_symmetry,
        convention="bunge",
        degrees=True,
    )


def test_mean_of_identical_orientations_is_that_orientation() -> None:
    orientations = make_set([(30.0, 40.0, 50.0)] * 5)
    mean = orientations.mean_orientation()
    spread = orientations.spread_angles_deg(reference=mean)
    assert_allclose(spread, 0.0, atol=1e-8)


def test_mean_of_small_cluster_recovers_center_without_symmetry() -> None:
    center = (20.0, 30.0, 40.0)
    perturbations = [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, -2.0, 0.0), (0.0, 2.0, 0.0)]
    euler = [tuple(np.add(center, delta)) for delta in perturbations]
    orientations = make_set(euler)  # type: ignore[arg-type]
    mean = orientations.mean_orientation()
    center_orientation = make_set([center])[0]
    distance_deg = np.rad2deg(mean.distance_to(center_orientation, symmetry_aware=False))
    assert distance_deg < 0.5


def test_mean_is_symmetry_aware_for_scattered_equivalents() -> None:
    crystal, specimen = make_frames()
    symmetry = cubic_symmetry(crystal)
    base = Orientation(
        rotation=Rotation.from_euler(15.0, 25.0, 35.0, convention="bunge", degrees=True),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    equivalents = base.equivalent_orientations()
    rng = np.random.default_rng(5)
    picks = rng.integers(0, len(equivalents), size=12)
    scattered = equivalents.subset(picks)
    mean = scattered.mean_orientation()
    distance_deg = np.rad2deg(mean.distance_to(base, symmetry_aware=True))
    assert distance_deg < 1e-6
    spread = scattered.spread_angles_deg(reference=mean)
    assert_allclose(spread, 0.0, atol=1e-6)


def test_mean_with_weights_prefers_heavier_cluster() -> None:
    orientations = make_set([(0.0, 10.0, 0.0), (0.0, 50.0, 0.0)])
    weighted_mean = orientations.mean_orientation(weights=[0.999, 0.001])
    near = make_set([(0.0, 10.0, 0.0)])[0]
    distance_deg = np.rad2deg(weighted_mean.distance_to(near, symmetry_aware=False))
    assert distance_deg < 0.5


def test_mean_orientation_input_validation() -> None:
    orientations = make_set([(0.0, 10.0, 0.0), (0.0, 12.0, 0.0)])
    with pytest.raises(ValueError, match="one value per orientation"):
        orientations.mean_orientation(weights=[1.0])
    with pytest.raises(ValueError, match="non-negative"):
        orientations.mean_orientation(weights=[1.0, -1.0])
    with pytest.raises(ValueError, match="sum to zero"):
        orientations.mean_orientation(weights=[0.0, 0.0])


def test_spread_reports_expected_disorientation_scale() -> None:
    crystal, specimen = make_frames()
    symmetry = cubic_symmetry(crystal)
    euler = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (355.0, 0.0, 0.0)]
    orientations = OrientationSet.from_euler_angles(
        np.asarray(euler),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        degrees=True,
    )
    mean = orientations.mean_orientation()
    spread = orientations.spread_angles_deg(reference=mean)
    assert spread.shape == (3,)
    assert spread.max() < 6.0
    assert spread.max() > 2.0
