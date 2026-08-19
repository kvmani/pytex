from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import FrameDomain, Handedness, ReferenceFrame, Rotation, SymmetrySpec, VectorSet


def make_crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def test_cubic_point_group_has_expected_order() -> None:
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=make_crystal_frame())
    assert symmetry.order == 24


def test_hexagonal_point_group_has_expected_order() -> None:
    symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=make_crystal_frame())
    assert symmetry.order == 12


def test_canonicalize_vector_folds_antipodal_direction() -> None:
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=make_crystal_frame())
    canonical = symmetry.canonicalize_vector([0.0, 0.0, -1.0], antipodal=True)
    assert_allclose(canonical, [0.0, 0.0, 1.0], atol=1e-12)


def test_apply_to_rotation_matrices_preserves_shape() -> None:
    symmetry = SymmetrySpec.from_point_group("4/mmm", reference_frame=make_crystal_frame())
    rotation = Rotation.from_bunge_euler(20.0, 30.0, 40.0).as_matrix()
    transformed = symmetry.apply_to_rotation_matrices(rotation, side="right")
    assert transformed.shape == (symmetry.order, 3, 3)


def test_cubic_sector_reduction_returns_standard_ipf_triangle_direction() -> None:
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=make_crystal_frame())
    reduced = symmetry.reduce_vector_to_fundamental_sector([0.0, 0.0, -1.0], antipodal=True)
    assert symmetry.vector_in_fundamental_sector(reduced, antipodal=True)
    assert reduced[2] >= reduced[0] >= reduced[1] >= 0.0


def test_hexagonal_sector_vertices_span_expected_wedge() -> None:
    symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=make_crystal_frame())
    sector = symmetry.fundamental_sector()
    assert sector.vertices.shape == (3, 3)
    assert_allclose(sector.vertices[0], [0.0, 0.0, 1.0], atol=1e-12)
    basal_angles = np.rad2deg(np.arctan2(sector.vertices[1:, 1], sector.vertices[1:, 0]))
    assert_allclose(basal_angles, [0.0, 30.0], atol=1e-8)


def test_reduce_vector_set_to_fundamental_sector_preserves_frame() -> None:
    frame = make_crystal_frame()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=frame)
    vectors = VectorSet(values=[[0.0, 0.0, -1.0]], reference_frame=frame)
    reduced = symmetry.reduce_vectors_to_fundamental_sector(vectors, antipodal=True)
    assert isinstance(reduced, VectorSet)
    assert reduced.reference_frame == frame
    assert symmetry.vector_in_fundamental_sector(reduced.values[0], antipodal=True)


def test_symmetry_spec_equality_is_well_defined_for_distinct_instances() -> None:
    frame = make_crystal_frame()
    left = SymmetrySpec.from_point_group("m-3m", reference_frame=frame)
    right = SymmetrySpec.from_point_group("m-3m", reference_frame=frame)
    assert left is not right
    assert left == right
    assert hash(left) == hash(right)
    assert left != SymmetrySpec.from_point_group("6/mmm", reference_frame=frame)
    assert {left: "cubic"}[right] == "cubic"


@pytest.mark.parametrize("point_group", ["m-3m", "6/mmm", "4/mmm", "2/m", "-1"])
@pytest.mark.parametrize("antipodal", [True, False])
def test_batch_sector_reduction_agrees_with_the_scalar_rule(
    point_group: str, antipodal: bool
) -> None:
    """The vectorized reduction is an optimization, so it may not answer differently.

    It skips renormalizing the orbit, which the operators being rotations makes
    redundant; the check is that the representatives still come back unit and
    identical to the one-vector-at-a-time rule, member for member.
    """

    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    rng = np.random.default_rng(4)
    directions = rng.normal(size=(400, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    reduced = np.asarray(
        symmetry.reduce_vectors_to_fundamental_sector(directions, antipodal=antipodal)
    )

    assert_allclose(np.linalg.norm(reduced, axis=1), 1.0, atol=1e-12)
    for index, direction in enumerate(directions):
        expected = symmetry.reduce_vector_to_fundamental_sector(direction, antipodal=antipodal)
        assert_allclose(reduced[index], np.asarray(expected), atol=1e-10)


def test_batch_sector_reduction_is_independent_of_the_chunk_size() -> None:
    """Blocking is bookkeeping; it may not reach the edge of the answer."""

    from pytex.core import symmetry as symmetry_module

    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    rng = np.random.default_rng(5)
    directions = rng.normal(size=(300, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    whole = np.asarray(symmetry.reduce_vectors_to_fundamental_sector(directions))
    original = symmetry_module._SECTOR_REDUCTION_CHUNK
    try:
        symmetry_module._SECTOR_REDUCTION_CHUNK = 7
        chunked = np.asarray(symmetry.reduce_vectors_to_fundamental_sector(directions))
    finally:
        symmetry_module._SECTOR_REDUCTION_CHUNK = original

    assert_allclose(whole, chunked, atol=0.0)
