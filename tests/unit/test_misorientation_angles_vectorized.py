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
from pytex.core.orientation import (
    _disorientation_medoid_index,
    _disorientation_scalar_projection,
    _reduced_pair_disorientation_angles,
    quaternions_to_matrices,
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


def _matrix_reference_angles(
    relative_matrices: np.ndarray,
    left_operators: np.ndarray,
    right_operators: np.ndarray,
) -> np.ndarray:
    """Trace-based ``min over S_l, S_r of angle(S_l M S_r^T)``.

    Deliberately written the slow, obvious way: it is the specification the
    quaternion/GEMM kernel in ``pytex.core.orientation`` must reproduce.
    """

    angles = np.empty(relative_matrices.shape[0], dtype=np.float64)
    for index, matrix in enumerate(relative_matrices):
        best = np.inf
        for left in left_operators:
            for right in right_operators:
                trace = np.trace(left @ matrix @ right.T)
                best = min(best, float(np.arccos(np.clip((trace - 1.0) * 0.5, -1.0, 1.0))))
        angles[index] = best
    return angles


def _random_rotation_matrices(count: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return quaternions_to_matrices(_haar_uniform_quaternions(count, rng))


@pytest.mark.parametrize(
    ("left_group", "right_group"),
    [
        ("m-3m", "m-3m"),
        ("6/mmm", "6/mmm"),
        ("mmm", "mmm"),
        ("-3m", "-3m"),
        ("1", "1"),
        ("m-3m", "6/mmm"),
    ],
)
def test_reduced_disorientation_kernel_matches_trace_reference(
    left_group: str, right_group: str
) -> None:
    crystal, _ = _frames()
    left_operators = SymmetrySpec.from_point_group(left_group, reference_frame=crystal).operators
    right_operators = SymmetrySpec.from_point_group(right_group, reference_frame=crystal).operators
    matrices = _random_rotation_matrices(24, seed=17)

    fast = _reduced_pair_disorientation_angles(matrices, left_operators, right_operators)
    reference = _matrix_reference_angles(matrices, left_operators, right_operators)

    np.testing.assert_allclose(fast, reference, atol=1e-11)


def test_scalar_projection_deduplicates_redundant_operator_pairs() -> None:
    """Only ``|a . q|`` matters, so sign-equal functionals are dropped.

    For same-phase cubic symmetry the 24 x 24 operator pairs collapse to 24
    distinct rows. This is what makes the symmetry reduction cheap, so it is
    pinned rather than left as an implementation detail.
    """

    crystal, _ = _frames()
    operators = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal).operators
    projection = _disorientation_scalar_projection(operators, operators)

    assert projection.shape == (24, 4)


def test_medoid_resolves_exact_ties_to_the_lowest_index() -> None:
    """A symmetric cluster has no unique medoid; the choice must be stable.

    Two members placed symmetrically about a centre have identical total
    disorientation, so the representative may not be decided by summation
    order, BLAS build or machine.
    """

    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    axis = np.array([0.0, 0.0, 1.0])
    quaternions = np.stack(
        [
            Rotation.from_axis_angle(axis, -np.deg2rad(2.0)).quaternion,
            Rotation.identity().quaternion,
            Rotation.from_axis_angle(axis, np.deg2rad(2.0)).quaternion,
        ]
    )
    members = OrientationSet.from_quaternions(
        quaternions, crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry
    )

    # Index 1 is the unique centre here.
    assert _disorientation_medoid_index(members) == 1

    # Drop the centre: indices 0 and 1 are now exactly tied, and the lower wins.
    tied = members.subset(np.array([0, 2]))
    assert _disorientation_medoid_index(tied) == 0


def test_medoid_is_independent_of_block_size() -> None:
    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    members = OrientationSet.from_quaternions(
        _haar_uniform_quaternions(40, np.random.default_rng(9)),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )

    unblocked = _disorientation_medoid_index(members, max_pairs_per_block=10_000_000)
    blocked = _disorientation_medoid_index(members, max_pairs_per_block=41)

    assert unblocked == blocked
