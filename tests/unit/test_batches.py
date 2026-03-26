from __future__ import annotations

import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    EulerSet,
    FrameDomain,
    Handedness,
    QuaternionSet,
    ReferenceFrame,
    RotationSet,
    VectorSet,
)


def make_frame(name: str = "specimen") -> ReferenceFrame:
    return ReferenceFrame(
        name=name,
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )


def test_vector_set_is_read_only_and_preserves_frame() -> None:
    frame = make_frame()
    vectors = VectorSet(values=[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]], reference_frame=frame)
    assert vectors.reference_frame == frame
    assert vectors.values.flags.writeable is False
    assert_allclose(vectors.normalized().values, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def test_euler_set_round_trips_through_rotation_set() -> None:
    eulers = EulerSet(angles=[[45.0, 35.0, 15.0], [90.0, 10.0, 5.0]], convention="bunge")
    rotation_set = eulers.to_rotation_set()
    recovered = rotation_set.as_euler_set(convention="bunge", degrees=True)
    assert recovered.convention == "bunge"
    assert_allclose(
        RotationSet.from_euler_set(recovered).as_matrices(),
        rotation_set.as_matrices(),
        atol=1e-8,
    )


def test_quaternion_set_normalizes_input() -> None:
    quaternion_set = QuaternionSet(quaternions=[[2.0, 0.0, 0.0, 0.0]])
    assert_allclose(quaternion_set.quaternions, [[1.0, 0.0, 0.0, 0.0]])


def test_rotation_set_applies_pairwise_to_vector_set() -> None:
    frame = make_frame()
    rotation_set = EulerSet(
        angles=[[90.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        convention="bunge",
        degrees=True,
    ).to_rotation_set()
    vectors = VectorSet(values=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], reference_frame=frame)
    rotated = rotation_set.apply(vectors)
    assert isinstance(rotated, VectorSet)
    assert rotated.reference_frame == frame
    assert rotated.values.shape == (2, 3)


def test_rotation_set_rejects_mismatched_vector_count() -> None:
    rotation_set = EulerSet(angles=[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]).to_rotation_set()
    with pytest.raises(ValueError):
        rotation_set.apply([[1.0, 0.0, 0.0]])
