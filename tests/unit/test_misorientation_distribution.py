from __future__ import annotations

import math

import numpy as np
import pytest

from pytex import (
    FrameDomain,
    Handedness,
    MisorientationDistribution,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def _frames() -> tuple[ReferenceFrame, ReferenceFrame]:
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
    return crystal, specimen


def _cubic_symmetry(crystal: ReferenceFrame) -> SymmetrySpec:
    return SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)


def test_random_baseline_matches_cubic_mackenzie_shape() -> None:
    crystal, _ = _frames()
    symmetry = _cubic_symmetry(crystal)
    mdf = MisorientationDistribution.random_baseline(symmetry, count=6000, seed=3)
    # The cubic disorientation is bounded by ~62.8 deg; the MacKenzie mean is
    # near 42-45 deg and the distribution peaks close to 45 deg.
    assert mdf.max_angle_deg <= 62.8 + 1e-6
    assert mdf.max_angle_deg > 55.0
    assert 38.0 < mdf.mean_angle_deg < 46.0
    density, edges = mdf.histogram(bins=20, angle_range=(0.0, 62.8))
    peak_center = 0.5 * (edges[int(np.argmax(density))] + edges[int(np.argmax(density)) + 1])
    assert 40.0 < peak_center < 50.0


def test_random_baseline_is_reproducible() -> None:
    crystal, _ = _frames()
    symmetry = _cubic_symmetry(crystal)
    first = MisorientationDistribution.random_baseline(symmetry, count=500, seed=7)
    second = MisorientationDistribution.random_baseline(symmetry, count=500, seed=7)
    np.testing.assert_allclose(first.angles_deg, second.angles_deg)


def test_uncorrelated_mdf_of_identical_orientations_is_zero() -> None:
    crystal, specimen = _frames()
    symmetry = _cubic_symmetry(crystal)
    orientations = OrientationSet.from_quaternions(
        np.tile(Rotation.identity().quaternion, (4, 1)),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    mdf = MisorientationDistribution.from_orientations(orientations)
    assert mdf.count == 6  # C(4, 2) unique pairs
    assert mdf.max_angle_deg == pytest.approx(0.0, abs=1e-9)


def test_mdf_recovers_known_pair_misorientation() -> None:
    crystal, specimen = _frames()
    symmetry = _cubic_symmetry(crystal)
    tilt = Rotation.from_axis_angle(np.array([0.0, 0.0, 1.0]), math.radians(30.0))
    orientations = OrientationSet.from_quaternions(
        np.array([Rotation.identity().quaternion, tilt.quaternion]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    mdf = MisorientationDistribution.from_orientations(orientations)
    assert mdf.count == 1
    # 30 deg about <001> is already the disorientation for cubic symmetry.
    assert mdf.angles_deg[0] == pytest.approx(30.0, abs=1e-6)


def test_correlated_pairs_flag_and_selection() -> None:
    crystal, specimen = _frames()
    symmetry = _cubic_symmetry(crystal)
    angles = [0.0, 10.0, 25.0, 40.0]
    orientations = OrientationSet.from_quaternions(
        np.array(
            [
                Rotation.from_axis_angle(np.array([0.0, 0.0, 1.0]), math.radians(a)).quaternion
                for a in angles
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    pairs = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int64)
    mdf = MisorientationDistribution.from_orientations(orientations, pairs=pairs)
    assert mdf.correlated is True
    assert mdf.count == 3
    np.testing.assert_allclose(sorted(mdf.angles_deg), [10.0, 15.0, 15.0], atol=1e-6)
