from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerDirectionSet,
    MillerPlane,
    MillerPlaneSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    from_orix_miller,
    to_orix_miller_direction,
    to_orix_miller_plane,
    to_orix_phase,
    to_orix_rotation,
)


def _make_phase() -> Phase:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    return Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def test_orix_miller_adapter_round_trip_for_scalars_and_sets() -> None:
    pytest.importorskip("orix")
    phase = _make_phase()

    plane = MillerPlane.from_hkl([1, 0, 0], phase=phase)
    direction = MillerDirection.from_uvw([1, 1, 0], phase=phase)
    plane_set = MillerPlaneSet.from_hkl([[1, 0, 0], [1, 1, 1]], phase=phase)
    direction_set = MillerDirectionSet.from_uvw([[1, 0, 0], [1, 1, 0]], phase=phase)

    orix_phase = to_orix_phase(phase)
    assert getattr(orix_phase, "name", None) == phase.name

    recovered_plane = from_orix_miller(to_orix_miller_plane(plane), phase=phase)
    recovered_direction = from_orix_miller(to_orix_miller_direction(direction), phase=phase)
    recovered_plane_set = from_orix_miller(to_orix_miller_plane(plane_set), phase=phase)
    recovered_direction_set = from_orix_miller(to_orix_miller_direction(direction_set), phase=phase)

    assert isinstance(recovered_plane, MillerPlane)
    assert isinstance(recovered_direction, MillerDirection)
    assert isinstance(recovered_plane_set, MillerPlaneSet)
    assert isinstance(recovered_direction_set, MillerDirectionSet)
    assert_array_equal(recovered_plane.indices, plane.indices)
    assert_array_equal(recovered_direction.indices, direction.indices)
    assert_array_equal(recovered_plane_set.indices, plane_set.indices)
    assert_array_equal(recovered_direction_set.indices, direction_set.indices)


def test_orix_rotation_adapter_preserves_scalar_first_quaternion_ordering() -> None:
    pytest.importorskip("orix")
    rotation = Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 3.0)

    orix_rotation = to_orix_rotation(rotation)
    quaternions = np.asarray(
        getattr(orix_rotation, "data", orix_rotation),
        dtype=np.float64,
    ).reshape(-1, 4)

    assert_allclose(quaternions[0], rotation.quaternion, atol=1e-12)
