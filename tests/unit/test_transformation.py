from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Phase,
    PhaseTransformationRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    TransformationVariant,
    VectorSet,
)


def make_phases() -> tuple[ReferenceFrame, ReferenceFrame, Phase, Phase]:
    crystal_parent = ReferenceFrame(
        name="parent_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    crystal_child = ReferenceFrame(
        name="child_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    symmetry_parent = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal_parent)
    symmetry_child = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal_child)
    lattice_parent = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal_parent)
    lattice_child = Lattice(3.2, 3.2, 3.2, 90.0, 90.0, 90.0, crystal_frame=crystal_child)
    parent = Phase(
        "austenite", lattice=lattice_parent, symmetry=symmetry_parent, crystal_frame=crystal_parent
    )
    child = Phase(
        "martensite", lattice=lattice_child, symmetry=symmetry_child, crystal_frame=crystal_child
    )
    return crystal_parent, crystal_child, parent, child


def test_orientation_relationship_maps_parent_vectors_to_child_frame() -> None:
    crystal_parent, crystal_child, parent, child = make_phases()
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 2.0),
    )
    vector_set = VectorSet(
        values=np.array([[1.0, 0.0, 0.0]]),
        reference_frame=crystal_parent,
    )
    mapped = relationship.map_parent_vector_to_child(vector_set)
    assert mapped.reference_frame == crystal_child
    assert_allclose(mapped.values[0], [0.0, 1.0, 0.0], atol=1e-8)


def test_orientation_relationship_can_be_built_from_parallel_plane_direction() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_parallel_plane_direction(
        name="plane_direction_or",
        parent_plane=CrystalPlane(
            miller=MillerIndex(indices=np.array([0, 0, 1]), phase=parent),
            phase=parent,
        ),
        child_plane=CrystalPlane(
            miller=MillerIndex(indices=np.array([0, 0, 1]), phase=child),
            phase=child,
        ),
        parent_direction=CrystalDirection([1.0, 0.0, 0.0], phase=parent),
        child_direction=CrystalDirection([0.0, 1.0, 0.0], phase=child),
    )
    assert relationship.parent_phase == parent
    assert relationship.child_phase == child
    assert len(relationship.parallel_planes) == 1
    assert len(relationship.parallel_directions) == 1
    assert_allclose(
        relationship.parent_to_child_rotation.as_matrix(),
        Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 2.0).as_matrix(),
        atol=1e-8,
    )


def test_orientation_relationship_from_parallel_plane_direction_rejects_phase_mismatch() -> None:
    _, _, parent, child = make_phases()
    with pytest.raises(
        ValueError,
        match=r"parent_plane\.phase must match parent_direction\.phase",
    ):
        OrientationRelationship.from_parallel_plane_direction(
            name="bad_or",
            parent_plane=CrystalPlane(
                miller=MillerIndex(indices=np.array([0, 0, 1]), phase=parent),
                phase=parent,
            ),
            child_plane=CrystalPlane(
                miller=MillerIndex(indices=np.array([0, 0, 1]), phase=child),
                phase=child,
            ),
            parent_direction=CrystalDirection([1.0, 0.0, 0.0], phase=child),
            child_direction=CrystalDirection([0.0, 1.0, 0.0], phase=child),
        )


def test_orientation_relationship_can_build_named_bain_correspondence() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent,
        child_phase=child,
    )
    assert relationship.name == "bain"
    assert_allclose(
        relationship.parallel_planes[0][0].miller.indices,
        np.array([0, 0, 1]),
    )
    assert_allclose(
        relationship.parallel_planes[0][1].miller.indices,
        np.array([0, 0, 1]),
    )
    mapped = relationship.map_parent_vector_to_child(
        VectorSet(values=np.array([[1.0, 1.0, 0.0]]), reference_frame=parent.crystal_frame)
    )
    assert_allclose(
        mapped.values[0] / np.linalg.norm(mapped.values[0]),
        [1.0, 0.0, 0.0],
        atol=1e-8,
    )


def test_bain_correspondence_rejects_non_cubic_phases() -> None:
    crystal_parent = ReferenceFrame(
        name="hcp_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    crystal_child = ReferenceFrame(
        name="bcc_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    parent = Phase(
        "hcp_parent",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=crystal_parent),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal_parent),
        crystal_frame=crystal_parent,
    )
    child = Phase(
        "bcc_child",
        lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=crystal_child),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal_child),
        crystal_frame=crystal_child,
    )
    with pytest.raises(ValueError, match="cubic parent phase"):
        OrientationRelationship.from_bain_correspondence(
            parent_phase=parent,
            child_phase=child,
        )


def test_orientation_relationship_can_build_named_nw_correspondence() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent,
        child_phase=child,
    )
    assert relationship.name == "nishiyama_wassermann"
    assert_allclose(
        relationship.parallel_planes[0][0].miller.indices,
        np.array([1, 1, 1]),
    )
    assert_allclose(
        relationship.parallel_planes[0][1].miller.indices,
        np.array([0, 1, 1]),
    )
    mapped = relationship.map_parent_vector_to_child(
        VectorSet(values=np.array([[1.0, -1.0, 0.0]]), reference_frame=parent.crystal_frame)
    )
    assert_allclose(
        mapped.values[0] / np.linalg.norm(mapped.values[0]),
        [1.0, 0.0, 0.0],
        atol=1e-8,
    )


def test_nw_correspondence_rejects_non_cubic_child() -> None:
    crystal_parent = ReferenceFrame(
        name="fcc_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    crystal_child = ReferenceFrame(
        name="hcp_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    parent = Phase(
        "fcc_parent",
        lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal_parent),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal_parent),
        crystal_frame=crystal_parent,
    )
    child = Phase(
        "hcp_child",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=crystal_child),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal_child),
        crystal_frame=crystal_child,
    )
    with pytest.raises(ValueError, match="cubic child phase"):
        OrientationRelationship.from_nishiyama_wassermann_correspondence(
            parent_phase=parent,
            child_phase=child,
        )


def test_orientation_relationship_generates_unique_variants() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.identity(),
    )
    variants = relationship.generate_variants()
    assert variants
    assert all(isinstance(variant, TransformationVariant) for variant in variants)
    assert len(
        {tuple(np.round(variant.parent_to_child_rotation.quaternion, 12)) for variant in variants}
    ) == len(variants)


def test_phase_transformation_record_requires_phase_alignment() -> None:
    crystal_parent, crystal_child, parent, child = make_phases()
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.identity(),
    )
    parent_orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal_parent,
        specimen_frame=specimen,
        symmetry=parent.symmetry,
        phase=parent,
    )
    child_orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal_child,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            )
        ]
    )
    record = PhaseTransformationRecord(
        name="record",
        orientation_relationship=relationship,
        parent_orientation=parent_orientation,
        child_orientations=child_orientations,
        variant_indices=np.array([1]),
    )
    assert record.variant_count == 1


def test_transformation_variant_rejects_wrong_phase_planes() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.identity(),
    )
    bad_plane = CrystalPlane(
        miller=MillerIndex(indices=np.array([1, 0, 0]), phase=parent), phase=parent
    )
    child_plane = CrystalPlane(
        miller=MillerIndex(indices=np.array([1, 0, 0]), phase=child), phase=child
    )
    with pytest.raises(ValueError):
        TransformationVariant(
            orientation_relationship=relationship,
            variant_index=1,
            parent_operator_index=0,
            child_operator_index=0,
            parent_to_child_rotation=Rotation.identity(),
            habit_plane_pairs=((child_plane, bad_plane),),
        )


def test_phase_transformation_record_predicted_orientations_follow_variant_indices() -> None:
    crystal_parent, crystal_child, parent, child = make_phases()
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 4.0),
    )
    variants = relationship.generate_variants()
    parent_orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal_parent,
        specimen_frame=specimen,
        symmetry=parent.symmetry,
        phase=parent,
    )
    child_orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal_child,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            ),
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal_child,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            ),
        ]
    )
    record = PhaseTransformationRecord(
        name="variant_record",
        orientation_relationship=relationship,
        parent_orientation=parent_orientation,
        child_orientations=child_orientations,
        variant_indices=np.array([variants[0].variant_index, variants[-1].variant_index]),
    )
    predicted = record.predicted_child_orientations()
    assert predicted.quaternions.shape == (2, 4)
    expected_first = variants[0].parent_to_child_rotation.compose(parent_orientation.rotation)
    expected_last = variants[-1].parent_to_child_rotation.compose(parent_orientation.rotation)
    assert_allclose(predicted.quaternions[0], expected_first.quaternion, atol=1e-8)
    assert_allclose(predicted.quaternions[1], expected_last.quaternion, atol=1e-8)
