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
    _align_quaternions_to_common_branch,
    _certified_single_branch_row_sums,
    _disorientation_medoid_index,
    _disorientation_scalar_projection,
    _exact_pair_angle_row_sums,
    _group_minimum_nonidentity_angle,
    _reduced_pair_disorientation_angles,
    matrices_to_quaternions,
    quaternions_from_axes_angles,
    quaternions_multiply,
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

    The full 4x4 conjugation of each surviving pair is returned alongside, and
    must be deduplicated in step with the rows: the reduction finds the winning
    pair from the scalar rows and then reads its angle off that pair's whole
    quaternion, so a row and an action that disagreed about which pair they
    describe would give the right winner and the wrong angle.
    """

    crystal, _ = _frames()
    operators = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal).operators
    projection, actions = _disorientation_scalar_projection(operators, operators)

    assert projection.shape == (24, 4)
    assert actions.shape == (24, 4, 4)
    np.testing.assert_allclose(actions[:, 0, :], projection, atol=0.0)


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


def _clustered_set(
    count: int,
    spread_deg: float,
    point_group: str,
    *,
    seed: int,
) -> OrientationSet:
    """A grain-like cluster whose members are stored on scrambled symmetry branches.

    A measurement reports whichever symmetry branch its indexing landed on, so
    a routine that folds the branches must be exercised against that, not
    against a set someone already tidied.
    """

    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    rng = np.random.default_rng(seed)
    centre = _haar_uniform_quaternions(1, rng)
    axes = rng.normal(size=(count, 3))
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    angles = np.deg2rad(rng.uniform(0.0, spread_deg, size=count))
    quaternions = quaternions_multiply(centre, quaternions_from_axes_angles(axes, angles))
    operator_quaternions = matrices_to_quaternions(symmetry.operators)
    branch = rng.integers(0, operator_quaternions.shape[0], size=count)
    quaternions = quaternions_multiply(quaternions, operator_quaternions[branch])
    quaternions *= np.where(rng.random(count) < 0.5, -1.0, 1.0)[:, None]
    return OrientationSet.from_quaternions(
        quaternions, crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry
    )


@pytest.mark.parametrize("point_group", ["m-3m", "6/mmm", "4/mmm", "mmm", "-1"])
@pytest.mark.parametrize("spread_deg", [0.0, 5.0, 20.0, 60.0, 180.0])
def test_certified_single_branch_row_sums_match_the_full_group_search(
    point_group: str, spread_deg: float
) -> None:
    """The fast medoid path must be exact, not merely close.

    Where the certificate holds the row sums must equal the full-group search
    to the accuracy ``arccos`` itself allows near zero; where it does not hold
    the fast path must decline rather than answer, and the medoid index must be
    the same either way.
    """

    members = _clustered_set(64, spread_deg, point_group, seed=len(point_group) * 101 + 7)
    quaternions = np.asarray(members.quaternions, dtype=np.float64)
    operators = members.symmetry.operators

    exact = _exact_pair_angle_row_sums(quaternions, operators, max_pairs_per_block=4_000_000)
    fast = _certified_single_branch_row_sums(quaternions, operators, max_pairs_per_block=4_000_000)

    if spread_deg == 0.0:
        # One orientation repeated: every total is zero, so there is nothing to
        # compare and nothing to choose between. The representative comes from
        # the degeneracy rule, tested separately.
        assert _disorientation_medoid_index(members) == 0
        return
    if fast is not None:
        # Each pair angle carries the arccos endpoint noise, about 2e-8 rad,
        # and the signs do not conspire across a row.
        assert np.abs(fast - exact).max() < 1e-6
    assert _disorientation_medoid_index(members) == int(
        np.flatnonzero(exact <= exact.min() + 1e-9 * max(abs(exact.min()), 1.0))[0]
    )


def test_certificate_declines_a_cluster_spanning_the_symmetry_angle() -> None:
    """A set too spread to certify must fall back, not answer from one branch."""

    crystal, _specimen = _frames()
    operators = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal).operators
    scattered = _haar_uniform_quaternions(32, np.random.default_rng(3))

    assert (
        _certified_single_branch_row_sums(scattered, operators, max_pairs_per_block=4_000_000)
        is None
    )


def test_group_minimum_nonidentity_angle_is_the_smallest_proper_rotation() -> None:
    """The certificate constant is a property of the group, and must be read off it."""

    crystal, _specimen = _frames()
    for point_group, expected_deg in (("m-3m", 90.0), ("6/mmm", 60.0), ("mmm", 180.0)):
        operators = SymmetrySpec.from_point_group(point_group, reference_frame=crystal).operators
        assert np.isclose(np.rad2deg(_group_minimum_nonidentity_angle(operators)), expected_deg)

    triclinic = SymmetrySpec.from_point_group("-1", reference_frame=crystal).operators
    assert not np.isfinite(_group_minimum_nonidentity_angle(triclinic))


def test_aligning_to_a_common_branch_preserves_the_orientations() -> None:
    """Folding branches may only change the representative, never the rotation."""

    crystal, _specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    operator_quaternions = matrices_to_quaternions(symmetry.operators)
    original = _haar_uniform_quaternions(16, np.random.default_rng(11))

    aligned = _align_quaternions_to_common_branch(original, operator_quaternions)

    relative = quaternions_multiply(original * np.array([1.0, -1.0, -1.0, -1.0]), aligned)
    equivalent = np.abs(relative @ operator_quaternions.T).max(axis=1)
    assert np.allclose(equivalent, 1.0)


def test_a_cluster_of_one_repeated_orientation_answers_with_the_lowest_index() -> None:
    """Where the totals are all zero, the choice must be a definition, not dust.

    Every member of such a set is equally central, so no arithmetic can single
    one out; letting the residue of the summation decide would make a grain's
    reference orientation a property of the machine that computed it.
    """

    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    single = _haar_uniform_quaternions(1, np.random.default_rng(2))
    operator_quaternions = matrices_to_quaternions(symmetry.operators)
    rng = np.random.default_rng(6)
    branch = rng.integers(0, operator_quaternions.shape[0], size=200)
    repeated = quaternions_multiply(single, operator_quaternions[branch])
    members = OrientationSet.from_quaternions(
        repeated, crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry
    )

    assert _disorientation_medoid_index(members) == 0
    # And it stays 0 however the members were shuffled or blocked.
    shuffled = members.subset(rng.permutation(200).astype(np.int64))
    assert _disorientation_medoid_index(shuffled, max_pairs_per_block=97) == 0


def test_a_real_separation_is_not_swallowed_by_the_tie_band() -> None:
    """The tie rule breaks ties; it may not overrule a member that is genuinely central."""

    crystal, specimen = _frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    axis = np.array([0.0, 0.0, 1.0])
    # Index 0 sits at one end, index 1 at the centre of the three.
    quaternions = np.stack(
        [
            Rotation.from_axis_angle(axis, np.deg2rad(-3.0)).quaternion,
            Rotation.from_axis_angle(axis, np.deg2rad(0.2)).quaternion,
            Rotation.from_axis_angle(axis, np.deg2rad(3.0)).quaternion,
            Rotation.from_axis_angle(axis, np.deg2rad(3.1)).quaternion,
        ]
    )
    members = OrientationSet.from_quaternions(
        quaternions, crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry
    )

    assert _disorientation_medoid_index(members) == 1
