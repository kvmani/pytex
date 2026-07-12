from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    FrameDomain,
    Handedness,
    Orientation,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.core.misorientation_distribution import _haar_uniform_quaternions


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


def _brute_force(
    left: OrientationSet, right: OrientationSet, *, symmetry_aware: bool
) -> np.ndarray:
    angles = np.empty((len(left), len(right)), dtype=np.float64)
    for i, qa in enumerate(left.quaternions):
        oa = Orientation(
            Rotation(qa),
            crystal_frame=left.crystal_frame,
            specimen_frame=left.specimen_frame,
            symmetry=left.symmetry,
        )
        for j, qb in enumerate(right.quaternions):
            ob = Orientation(
                Rotation(qb),
                crystal_frame=right.crystal_frame,
                specimen_frame=right.specimen_frame,
                symmetry=right.symmetry,
            )
            angles[i, j] = oa.distance_to(ob, symmetry_aware=symmetry_aware)
    return angles


@pytest.mark.parametrize("symmetry_aware", [True, False])
def test_vectorized_matches_per_pair_reference(symmetry_aware: bool) -> None:
    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    left = OrientationSet.from_quaternions(
        _haar_uniform_quaternions(11, np.random.default_rng(4)),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    right = OrientationSet.from_quaternions(
        _haar_uniform_quaternions(9, np.random.default_rng(5)),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    fast = left.misorientation_angles_to(right, symmetry_aware=symmetry_aware)
    reference = _brute_force(left, right, symmetry_aware=symmetry_aware)
    np.testing.assert_allclose(fast, reference, atol=1e-10)


def test_symmetry_equivalents_measure_zero() -> None:
    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    base = OrientationSet.from_quaternions(
        Rotation.identity().quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    # a 90-degree rotation about <001> is a cubic symmetry operator
    equivalent = OrientationSet.from_quaternions(
        Rotation.from_axis_angle(np.array([0.0, 0.0, 1.0]), np.pi / 2.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    angle = base.misorientation_angles_to(equivalent)
    assert angle.shape == (1, 1)
    assert angle[0, 0] == pytest.approx(0.0, abs=1e-9)


def test_frame_mismatch_is_rejected() -> None:
    crystal, specimen = _frames()
    other_specimen = ReferenceFrame(
        name="other",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    left = OrientationSet.from_quaternions(
        Rotation.identity().quaternion[None, :], crystal_frame=crystal, specimen_frame=specimen
    )
    right = OrientationSet.from_quaternions(
        Rotation.identity().quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=other_specimen,
    )
    with pytest.raises(ValueError, match="specimen frame"):
        left.misorientation_angles_to(right)
