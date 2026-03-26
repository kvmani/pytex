from __future__ import annotations

import numpy as np
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
