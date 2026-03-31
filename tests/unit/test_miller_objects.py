from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerDirectionSet,
    MillerIndex,
    MillerPlane,
    MillerPlaneSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
    angle_dir_dir_rad,
    angle_dir_plane_inclination_rad,
    angle_dir_plane_normal_rad,
    angle_plane_plane_rad,
    antipodal_keys,
    canonicalize_sign,
    direction_uvtw_to_uvw_array,
    direction_uvw_to_uvtw_array,
    plane_hkil_to_hkl_array,
    plane_hkl_to_hkil_array,
    project_directions_onto_planes,
    reduce_indices,
)


def _make_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def _make_cubic_phase() -> Phase:
    crystal = _make_frame()
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    return Phase("cubic_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def _make_hex_phase() -> Phase:
    crystal = _make_frame()
    lattice = Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal)
    return Phase("hcp_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def test_miller_construction_and_validation() -> None:
    phase = _make_cubic_phase()

    plane = MillerPlane.from_hkl([2, 0, 0], phase=phase)
    direction = MillerDirection.from_uvw([1, 1, 0], phase=phase)
    plane_set = MillerPlaneSet.from_hkl([[2, 0, 0], [1, 1, 1]], phase=phase)
    direction_set = MillerDirectionSet.from_uvw([[1, 0, 0], [1, 1, 0]], phase=phase)

    assert_array_equal(plane.indices, np.array([2, 0, 0], dtype=np.int64))
    assert_array_equal(direction.indices, np.array([1, 1, 0], dtype=np.int64))
    assert plane_set.indices.shape == (2, 3)
    assert direction_set.indices.shape == (2, 3)

    with pytest.raises(ValueError, match="zero triplet"):
        MillerPlane.from_hkl([0, 0, 0], phase=phase)
    with pytest.raises(ValueError, match="zero triplet"):
        MillerDirection.from_uvw([0, 0, 0], phase=phase)
    with pytest.raises(ValueError, match="shape"):
        MillerPlaneSet.from_hkl([[1, 0], [0, 1]], phase=phase)
    with pytest.raises(ValueError, match="satisfy i = -\\(h \\+ k\\)"):
        MillerPlane.from_hkil([1, 0, 0, 0], phase=phase)
    with pytest.raises(ValueError, match="satisfy U \\+ V \\+ T = 0"):
        MillerDirection.from_UVTW([1, 0, 0, 0], phase=phase)


def test_reduction_sign_canonicalization_and_antipodal_keys() -> None:
    values = np.array([[2, -4, 2], [-2, 4, -2], [0, -6, 0]], dtype=np.int64)

    reduced = reduce_indices(values)
    signs = canonicalize_sign(reduced)
    keys = antipodal_keys(values)

    assert_array_equal(reduced, np.array([[1, -2, 1], [-1, 2, -1], [0, -1, 0]], dtype=np.int64))
    assert_array_equal(signs, np.array([[1, -2, 1], [1, -2, 1], [0, 1, 0]], dtype=np.int64))
    assert_array_equal(keys, np.array([[1, -2, 1], [1, -2, 1], [0, 1, 0]], dtype=np.int64))


def test_vectorized_hex_conversions_round_trip() -> None:
    hkl = np.array([[1, 0, 0], [2, -1, 3]], dtype=np.int64)
    uvw = np.array([[1, 0, 0], [2, 1, 0]], dtype=np.int64)

    hkil = plane_hkl_to_hkil_array(hkl)
    uvtw = direction_uvw_to_uvtw_array(uvw)

    assert_array_equal(plane_hkil_to_hkl_array(hkil), hkl)
    assert_array_equal(direction_uvtw_to_uvw_array(uvtw), reduce_indices(uvw))
    assert_array_equal(hkil[0], np.array([1, 0, -1, 0], dtype=np.int64))
    assert_array_equal(uvtw[0], np.array([2, -1, -1, 0], dtype=np.int64))


def test_plane_and_direction_compatibility_helpers_round_trip() -> None:
    phase = _make_cubic_phase()
    crystal_plane = CrystalPlane(miller=MillerIndex([2, 0, 0], phase=phase), phase=phase)
    zone_axis = ZoneAxis([2, 0, 0], phase=phase)

    plane = MillerPlane.from_crystal_plane(crystal_plane)
    direction = MillerDirection.from_zone_axis(zone_axis)

    converted_plane = plane.to_crystal_plane()
    converted_axis = direction.to_zone_axis()

    assert converted_plane.phase == crystal_plane.phase
    assert converted_axis.phase == zone_axis.phase
    assert_array_equal(converted_plane.miller.indices, crystal_plane.miller.indices)
    assert_array_equal(converted_axis.indices, zone_axis.indices)


def test_plane_d_spacing_parity_for_cubic() -> None:
    phase = _make_cubic_phase()

    d_200 = MillerPlane.from_hkl([2, 0, 0], phase=phase).d_spacing_angstrom
    d_111 = MillerPlane.from_hkl([1, 1, 1], phase=phase).d_spacing_angstrom

    assert_allclose(d_200, phase.lattice.a / 2.0, atol=1e-12)
    assert_allclose(d_111, phase.lattice.a / np.sqrt(3.0), atol=1e-12)


def test_plane_d_spacing_parity_for_hexagonal_lattice() -> None:
    phase = _make_hex_phase()
    planes = MillerPlaneSet.from_hkl([[1, 0, 0], [1, 0, 1], [0, 0, 2]], phase=phase)
    d_spacings = planes.d_spacings_angstrom()

    h = planes.indices[:, 0].astype(np.float64)
    k = planes.indices[:, 1].astype(np.float64)
    ell = planes.indices[:, 2].astype(np.float64)
    inverse_d_sq = ((4.0 / 3.0) * (h**2 + h * k + k**2) / (phase.lattice.a**2)) + (
        ell**2 / (phase.lattice.c**2)
    )
    analytic = 1.0 / np.sqrt(inverse_d_sq)

    assert_allclose(d_spacings, analytic, atol=1e-12)


def test_cubic_vector_parity_for_plane_normals_and_directions() -> None:
    phase = _make_cubic_phase()
    plane = MillerPlane.from_hkl([1, 1, 1], phase=phase)
    direction = MillerDirection.from_uvw([1, 1, 0], phase=phase)

    assert_allclose(plane.normal_cartesian, np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0), atol=1e-12)
    assert_allclose(
        direction.unit_vector_cartesian,
        np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
        atol=1e-12,
    )


def test_symmetry_families_for_cubic_planes_and_directions() -> None:
    phase = _make_cubic_phase()

    plane_family = MillerPlane.from_hkl([1, 0, 0], phase=phase).symmetry_equivalents()
    direction_family = MillerDirection.from_uvw([1, 0, 0], phase=phase).symmetry_equivalents()

    assert plane_family.indices.shape == (3, 3)
    assert direction_family.indices.shape == (3, 3)
    assert_array_equal(
        plane_family.indices,
        np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.int64),
    )
    assert_array_equal(
        direction_family.indices,
        np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.int64),
    )


def test_angle_functions_and_angle_matrices() -> None:
    phase = _make_cubic_phase()
    dir_100 = MillerDirection.from_uvw([1, 0, 0], phase=phase)
    dir_110 = MillerDirection.from_uvw([1, 1, 0], phase=phase)
    plane_100 = MillerPlane.from_hkl([1, 0, 0], phase=phase)
    plane_110 = MillerPlane.from_hkl([1, 1, 0], phase=phase)

    assert_allclose(angle_dir_dir_rad(dir_100, dir_110), np.pi / 4.0, atol=1e-12)
    assert_allclose(angle_plane_plane_rad(plane_100, plane_110), np.pi / 4.0, atol=1e-12)
    assert_allclose(angle_dir_plane_normal_rad(dir_100, plane_110), np.pi / 4.0, atol=1e-12)
    assert_allclose(
        angle_dir_plane_inclination_rad(dir_100, plane_110),
        np.pi / 4.0,
        atol=1e-12,
    )

    directions = MillerDirectionSet.from_uvw([[1, 0, 0], [1, 1, 0]], phase=phase)
    matrix = directions.angle_matrix_rad()
    assert_allclose(matrix, matrix.T, atol=1e-12)
    assert_allclose(matrix[0, 1], np.pi / 4.0, atol=1e-12)


def test_projection_returns_projected_vectors_and_degenerate_mask() -> None:
    phase = _make_cubic_phase()
    directions = MillerDirectionSet.from_uvw([[0, 0, 1], [1, 1, 0]], phase=phase)
    planes = MillerPlaneSet.from_hkl([[0, 0, 1], [0, 0, 1]], phase=phase)

    projected, degenerate = project_directions_onto_planes(directions, planes)

    assert_array_equal(degenerate, np.array([True, False]))
    assert_allclose(projected[0], np.zeros(3), atol=1e-12)
    assert_allclose(projected[1], np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0), atol=1e-12)


def test_projection_requires_phase_consistency() -> None:
    cubic = _make_cubic_phase()
    hexagonal = _make_hex_phase()
    direction = MillerDirection.from_uvw([1, 0, 0], phase=cubic)
    plane = MillerPlane.from_hkl([0, 0, 1], phase=hexagonal)

    with pytest.raises(ValueError, match="phase must match"):
        project_directions_onto_planes(direction, plane)


def test_batch_symmetry_indices_and_unique_inverse_mapping() -> None:
    phase = _make_cubic_phase()
    planes = MillerPlaneSet.from_hkl([[1, 0, 0], [2, 0, 0], [0, 1, 0]], phase=phase)
    unique_planes, inverse = planes.unique()
    expanded_indices, mask = planes.symmetry_equivalent_indices()

    assert_array_equal(unique_planes.indices, np.array([[0, 1, 0], [1, 0, 0]], dtype=np.int64))
    assert_array_equal(inverse, np.array([1, 1, 0], dtype=np.int64))
    assert expanded_indices.shape[:2] == mask.shape
    assert expanded_indices.shape[0] == 3
    assert np.all(np.sum(mask, axis=1) == 3)
