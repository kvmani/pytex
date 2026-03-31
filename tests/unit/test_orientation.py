from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalDirection,
    CrystalPlane,
    EulerSet,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    RotationSet,
    SymmetrySpec,
    VectorSet,
)


def make_orientation_frames() -> tuple[ReferenceFrame, ReferenceFrame]:
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


def make_phase_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def test_bunge_zero_yields_identity_rotation() -> None:
    rotation = Rotation.from_bunge_euler(0.0, 0.0, 0.0)
    assert_allclose(rotation.as_matrix(), np.eye(3))


def test_rotation_matrix_round_trip_preserves_orientation() -> None:
    rotation = Rotation.from_bunge_euler(45.0, 35.0, 15.0)
    recovered = Rotation.from_matrix(rotation.as_matrix())
    assert_allclose(recovered.as_matrix(), rotation.as_matrix(), atol=1e-8)


def test_bunge_euler_round_trip_preserves_orientation() -> None:
    rotation = Rotation.from_bunge_euler(270.0, 10.0, 45.0)
    recovered = Rotation.from_bunge_euler(*rotation.to_bunge_euler())
    assert_allclose(recovered.as_matrix(), rotation.as_matrix(), atol=1e-8)


def test_orientation_set_builds_from_orientations() -> None:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ),
        Orientation(
            Rotation.from_bunge_euler(45.0, 35.0, 15.0),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ),
    ]
    orientation_set = OrientationSet.from_orientations(orientations)
    assert len(orientation_set) == 2
    assert orientation_set.as_matrices().shape == (2, 3, 3)


def test_misorientation_identity_for_same_orientation() -> None:
    crystal, specimen = make_orientation_frames()
    orientation = Orientation(
        Rotation.from_bunge_euler(30.0, 10.0, 20.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=SymmetrySpec.identity(reference_frame=crystal),
    )
    misorientation = orientation.misorientation_to(orientation)
    assert_allclose(misorientation.as_matrix(), np.eye(3), atol=1e-8)


def test_orientation_set_rejects_mixed_specimen_frames() -> None:
    crystal, specimen = make_orientation_frames()
    specimen_alt = ReferenceFrame(
        name="specimen_alt",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ),
        Orientation(
            Rotation.from_bunge_euler(10.0, 20.0, 30.0),
            crystal_frame=crystal,
            specimen_frame=specimen_alt,
            symmetry=symmetry,
        ),
    ]
    with pytest.raises(ValueError):
        OrientationSet.from_orientations(orientations)


def test_orientation_set_preserves_shared_provenance() -> None:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    provenance = ProvenanceRecord.minimal("demo-import", note="shared")
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
            provenance=provenance,
        ),
        Orientation(
            Rotation.from_bunge_euler(10.0, 20.0, 30.0),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
            provenance=provenance,
        ),
    ]
    orientation_set = OrientationSet.from_orientations(orientations)
    assert orientation_set.provenance == provenance


def test_orientation_set_rejects_mixed_provenance_until_aggregate_model_exists() -> None:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
            provenance=ProvenanceRecord.minimal("vendor-a"),
        ),
        Orientation(
            Rotation.from_bunge_euler(10.0, 20.0, 30.0),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
            provenance=ProvenanceRecord.minimal("vendor-b"),
        ),
    ]
    with pytest.raises(ValueError):
        OrientationSet.from_orientations(orientations)


def test_orientation_set_accepts_euler_set() -> None:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    eulers = EulerSet(angles=[[45.0, 35.0, 15.0], [0.0, 0.0, 0.0]])
    orientation_set = OrientationSet.from_euler_angles(
        eulers,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    assert len(orientation_set) == 2
    assert orientation_set.as_euler_set().convention == "bunge"


def test_orientation_maps_vector_set_between_frames() -> None:
    crystal, specimen = make_orientation_frames()
    orientation = Orientation(
        Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=SymmetrySpec.identity(reference_frame=crystal),
    )
    vectors = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=crystal)
    mapped = orientation.map_crystal_vector(vectors)
    assert isinstance(mapped, VectorSet)
    assert mapped.reference_frame == specimen


def test_orientation_set_maps_vector_set_pairwise() -> None:
    crystal, specimen = make_orientation_frames()
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    rotations = RotationSet.from_rotations(
        [Rotation.identity(), Rotation.from_bunge_euler(0.0, 0.0, 0.0)]
    )
    orientation_set = OrientationSet(
        quaternions=rotations.quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    vectors = VectorSet(
        values=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        reference_frame=crystal,
    )
    mapped = orientation_set.map_crystal_directions(vectors)
    assert isinstance(mapped, VectorSet)
    assert mapped.reference_frame == specimen


def test_orientation_set_from_axes_angles_preserves_phase_and_symmetry() -> None:
    crystal, specimen, phase = make_phase_context()
    orientations = OrientationSet.from_axes_angles(
        axes=np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
        angles_rad=np.array([np.pi / 2.0, np.pi / 3.0]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    axes, angles = orientations.to_axes_angles()
    recovered = OrientationSet.from_axes_angles(
        axes,
        angles,
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    assert recovered.phase == phase
    assert recovered.symmetry == phase.symmetry
    assert_allclose(recovered.as_matrices(), orientations.as_matrices(), atol=1e-8)


def test_orientation_set_named_quaternion_and_matrix_constructors_match() -> None:
    crystal, specimen, phase = make_phase_context()
    base = OrientationSet.from_euler_angles(
        [[0.0, 0.0, 0.0], [35.0, 25.0, 10.0]],
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    from_quaternions = OrientationSet.from_quaternions(
        base.quaternions,
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    from_matrices = OrientationSet.from_matrices(
        base.as_matrices(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    assert_allclose(from_quaternions.as_matrices(), base.as_matrices(), atol=1e-8)
    assert_allclose(from_matrices.as_matrices(), base.as_matrices(), atol=1e-8)


def test_orientation_from_plane_direction_maps_cube_component_to_identity() -> None:
    _, specimen, phase = make_phase_context()
    plane = CrystalPlane(miller=MillerIndex([0, 0, 1], phase=phase), phase=phase)
    direction = CrystalDirection([1.0, 0.0, 0.0], phase=phase)
    orientation = Orientation.from_plane_direction(
        plane,
        direction,
        specimen_frame=specimen,
    )
    assert_allclose(orientation.as_matrix(), np.eye(3), atol=1e-8)
    assert_allclose(orientation.map_crystal_vector(plane.normal), [0.0, 0.0, 1.0], atol=1e-8)
    assert_allclose(
        orientation.map_crystal_vector(direction.unit_vector),
        [1.0, 0.0, 0.0],
        atol=1e-8,
    )


def test_orientation_from_plane_direction_projects_direction_into_plane() -> None:
    _, specimen, phase = make_phase_context()
    plane = CrystalPlane(miller=MillerIndex([0, 0, 1], phase=phase), phase=phase)
    direction = CrystalDirection([1.0, 0.0, 1.0], phase=phase)
    orientation = Orientation.from_plane_direction(
        plane,
        direction,
        specimen_frame=specimen,
    )
    projected = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    assert_allclose(orientation.map_crystal_vector(plane.normal), [0.0, 0.0, 1.0], atol=1e-8)
    assert_allclose(orientation.map_crystal_vector(projected), [1.0, 0.0, 0.0], atol=1e-8)


def test_orientation_set_from_plane_direction_vectorized_path_preserves_invariants() -> None:
    _, specimen, phase = make_phase_context()
    planes = [
        CrystalPlane(miller=MillerIndex([0, 0, 1], phase=phase), phase=phase),
        CrystalPlane(miller=MillerIndex([0, 0, 1], phase=phase), phase=phase),
    ]
    directions = [
        CrystalDirection([1.0, 0.0, 0.0], phase=phase),
        CrystalDirection([0.0, 1.0, 0.0], phase=phase),
    ]
    orientations = OrientationSet.from_plane_direction(
        planes,
        directions,
        specimen_frame=specimen,
    )
    matrices = orientations.as_matrices()
    mapped_normals = orientations.map_crystal_directions(
        np.vstack([plane.normal for plane in planes])
    )
    mapped_directions = orientations.map_crystal_directions(
        np.vstack([direction.unit_vector for direction in directions])
    )
    assert_allclose(mapped_normals, [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], atol=1e-8)
    assert_allclose(mapped_directions, [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], atol=1e-8)
    assert_allclose(np.linalg.det(matrices), [1.0, 1.0], atol=1e-8)
    assert_allclose(
        np.einsum("nij,nkj->nik", matrices, matrices),
        np.broadcast_to(np.eye(3), matrices.shape),
        atol=1e-8,
    )
    raw = OrientationSet.from_plane_direction(
        plane=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        direction=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        specimen_frame=specimen,
        phase=phase,
    )
    assert_allclose(raw.as_matrices(), matrices, atol=1e-8)
