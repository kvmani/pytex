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
    ORCharacterizationReport,
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Phase,
    PhaseTransformationRecord,
    RationalizedORResult,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    TransformationVariant,
    VectorSet,
    find_parallel_directions,
    find_parallel_planes,
    fit_orientation_relationship,
    map_plane_across_variants,
    or_deviation,
    phases_semantically_match,
    reconstruct_parent_orientation,
    select_variants,
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
    expected_first = parent_orientation.rotation.compose(
        variants[0].parent_to_child_rotation.inverse()
    )
    expected_last = parent_orientation.rotation.compose(
        variants[-1].parent_to_child_rotation.inverse()
    )
    # q and -q describe the same rotation; compare up to the double cover so
    # the pinned behavior is the rotation itself, not an incidental sign.
    assert abs(float(predicted.quaternions[0] @ expected_first.quaternion)) == pytest.approx(
        1.0, abs=1e-8
    )
    assert abs(float(predicted.quaternions[1] @ expected_last.quaternion)) == pytest.approx(
        1.0, abs=1e-8
    )


def make_hcp_child() -> Phase:
    crystal_child = ReferenceFrame(
        name="hcp_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    symmetry_child = SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal_child)
    lattice_child = Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=crystal_child)
    return Phase(
        "alpha_titanium",
        lattice=lattice_child,
        symmetry=symmetry_child,
        crystal_frame=crystal_child,
    )


def _symmetry_reduced_angle_deg(
    relationship_a: OrientationRelationship,
    relationship_b: OrientationRelationship,
) -> float:
    parent_ops = relationship_a.parent_phase.symmetry.operators
    child_ops = relationship_a.child_phase.symmetry.operators
    matrix_a = relationship_a.parent_to_child_rotation.as_matrix()
    matrix_b = relationship_b.parent_to_child_rotation.as_matrix()
    best = 180.0
    for child_op in child_ops:
        product_left = child_op @ matrix_a
        for parent_op in parent_ops:
            relative = product_left @ parent_op @ matrix_b.T
            cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
            best = min(best, float(np.degrees(np.arccos(cosine))))
    return best


def test_named_or_variant_counts_match_literature() -> None:
    _, _, parent, child = make_phases()
    expected_counts = {
        "bain": (OrientationRelationship.from_bain_correspondence, 3),
        "kurdjumov_sachs": (OrientationRelationship.from_kurdjumov_sachs_correspondence, 24),
        "nishiyama_wassermann": (
            OrientationRelationship.from_nishiyama_wassermann_correspondence,
            12,
        ),
        "greninger_troiano": (
            OrientationRelationship.from_greninger_troiano_correspondence,
            24,
        ),
        "pitsch": (OrientationRelationship.from_pitsch_correspondence, 12),
    }
    for name, (constructor, expected) in expected_counts.items():
        relationship = constructor(parent_phase=parent, child_phase=child)
        variants = relationship.generate_variants()
        assert len(variants) == expected, name
        assert [variant.variant_index for variant in variants] == list(
            range(1, expected + 1)
        )


def test_burgers_or_requires_hexagonal_child_and_yields_12_variants() -> None:
    _, _, parent, cubic_child = make_phases()
    hcp_child = make_hcp_child()
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=hcp_child
    )
    assert len(relationship.generate_variants()) == 12
    # the defining parallelism: {110}_bcc || (0001)_hcp and <-111> || <11-20>
    basal_normal = relationship.parallel_planes[0][1].normal
    mapped_normal = relationship.map_parent_vector_to_child(
        relationship.parallel_planes[0][0].normal
    )
    assert abs(float(np.dot(mapped_normal, basal_normal))) == pytest.approx(1.0, abs=1e-9)
    with pytest.raises(ValueError, match="hexagonal child"):
        OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=cubic_child
        )


def test_ks_gt_pitsch_geometric_constants() -> None:
    # the classic fcc->bcc OR triangle: KS-NW 5.26 deg, GT between them,
    # Pitsch mirrored 5.26 deg from KS
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    pitsch = OrientationRelationship.from_pitsch_correspondence(
        parent_phase=parent, child_phase=child
    )
    assert _symmetry_reduced_angle_deg(ks, nw) == pytest.approx(5.264, abs=0.01)
    assert _symmetry_reduced_angle_deg(gt, ks) == pytest.approx(2.404, abs=0.01)
    assert _symmetry_reduced_angle_deg(gt, nw) == pytest.approx(2.861, abs=0.01)
    assert _symmetry_reduced_angle_deg(pitsch, ks) == pytest.approx(5.264, abs=0.01)


def test_generate_variants_child_symmetry_reduction_default() -> None:
    _, _, parent, child = make_phases()
    nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    reduced = nw.generate_variants()
    raw = nw.generate_variants(reduce_by_child_symmetry=False)
    assert len(reduced) == 12
    # the historical raw enumeration counts every symmetry-equivalent
    # description separately
    assert len(raw) > len(reduced)
    # each reduced variant is crystallographically distinct: no two are
    # related by a child symmetry operator
    child_ops = child.symmetry.operators
    matrices = [variant.parent_to_child_rotation.as_matrix() for variant in reduced]
    for i in range(len(matrices)):
        for j in range(i + 1, len(matrices)):
            for child_op in child_ops:
                relative = child_op @ matrices[i] @ matrices[j].T
                cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
                assert np.degrees(np.arccos(cosine)) > 1e-6


def test_standard_catalogs_resolve_named_relationships() -> None:
    from pytex.core import standard_bcc_hcp_relationships, standard_fcc_bcc_relationships

    _, _, parent, child = make_phases()
    catalog = standard_fcc_bcc_relationships(parent_phase=parent, child_phase=child)
    assert set(catalog.names()) == {
        "bain",
        "kurdjumov_sachs",
        "nishiyama_wassermann",
        "greninger_troiano",
        "pitsch",
    }
    assert len(catalog.get("kurdjumov_sachs").generate_variants()) == 24
    hcp_catalog = standard_bcc_hcp_relationships(
        parent_phase=parent, child_phase=make_hcp_child()
    )
    assert hcp_catalog.names() == ("burgers",)


def test_ks_intervariant_angles_match_morito_table() -> None:
    from pytex.core import intervariant_misorientation_angles_deg

    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    angles = intervariant_misorientation_angles_deg(ks)
    assert angles.shape == (24, 24)
    assert_allclose(angles, angles.T, atol=1e-9)
    assert_allclose(np.diag(angles), 0.0, atol=1e-9)
    # the published KS intervariant disorientation classes (Morito et al.,
    # Acta Mater. 51 (2003) 1789, Table 3)
    expected_classes = np.array(
        [10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00]
    )
    off_diagonal = angles[~np.eye(24, dtype=bool)]
    observed_classes = np.unique(np.round(off_diagonal, 2))
    assert observed_classes.size == expected_classes.size
    assert_allclose(observed_classes, expected_classes, atol=0.01)


def test_intervariant_misorientations_axes_are_consistent() -> None:
    from pytex.core import intervariant_misorientation_angles_deg, intervariant_misorientations

    _, _, parent, child = make_phases()
    nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    pairs = intervariant_misorientations(nw)
    assert len(pairs) == 12 * 11 // 2
    matrix = intervariant_misorientation_angles_deg(nw)
    for pair in pairs:
        assert pair.angle_deg == pytest.approx(
            matrix[pair.variant_a - 1, pair.variant_b - 1], abs=1e-6
        )
        assert np.linalg.norm(pair.axis_child_frame) == pytest.approx(1.0)


def test_select_variants_recovers_planted_assignments() -> None:
    from pytex.core import select_variants

    crystal_parent, crystal_child, parent, child = make_phases()
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    variants = ks.generate_variants()
    parent_orientation = Orientation(
        rotation=Rotation.from_bunge_euler(20.0, 35.0, 50.0),
        crystal_frame=crystal_parent,
        specimen_frame=specimen,
        symmetry=parent.symmetry,
        phase=parent,
    )
    planted = [3, 7, 3, 18]
    perturbation = Rotation.from_axis_angle([0.0, 0.0, 1.0], np.deg2rad(0.5))
    quaternions = np.stack(
        [
            perturbation.compose(
                parent_orientation.rotation.compose(
                    variants[index - 1].parent_to_child_rotation.inverse()
                )
            ).quaternion
            for index in planted
        ]
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=crystal_child,
        specimen_frame=specimen,
        symmetry=child.symmetry,
        phase=child,
    )
    record = PhaseTransformationRecord(
        name="planted",
        orientation_relationship=ks,
        parent_orientation=parent_orientation,
        child_orientations=children,
    )
    report = select_variants(record)
    assert report.variant_indices.tolist() == planted
    assert_allclose(report.scores_deg, 0.5, atol=1e-6)
    frequencies = report.variant_frequencies(variant_count=24)
    assert frequencies.sum() == 4
    assert frequencies[2] == 2 and frequencies[6] == 1 and frequencies[17] == 1
    with pytest.raises(ValueError, match="variant_count"):
        report.variant_frequencies(variant_count=2)


def test_phases_semantically_match_handles_independent_construction_and_none() -> None:
    _, _, parent, child = make_phases()
    _, _, parent_again, _ = make_phases()
    assert phases_semantically_match(parent, parent_again)
    assert not phases_semantically_match(parent, child)
    assert not phases_semantically_match(parent, None)
    assert phases_semantically_match(None, None)


def test_parallel_directions_are_typed_crystal_directions_with_index_meaning() -> None:
    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    parent_direction, child_direction = relationship.parallel_directions[0]
    assert isinstance(parent_direction, CrystalDirection)
    assert isinstance(child_direction, CrystalDirection)
    assert_allclose(parent_direction.coordinates, [-1.0, 0.0, 1.0], atol=1e-12)
    assert_allclose(child_direction.coordinates, [-1.0, -1.0, 1.0], atol=1e-12)


def test_parallel_directions_accept_legacy_cartesian_vectors() -> None:
    _, _, parent, child = make_phases()
    base = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=child
    )
    legacy = OrientationRelationship(
        name="legacy_bain",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=base.parent_to_child_rotation,
        parallel_directions=((np.array([1.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])),),
    )
    parent_direction, child_direction = legacy.parallel_directions[0]
    assert isinstance(parent_direction, CrystalDirection)
    assert_allclose(
        parent_direction.unit_vector, np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0), atol=1e-12
    )
    assert_allclose(child_direction.unit_vector, [1.0, 0.0, 0.0], atol=1e-12)


def test_parallel_directions_reject_phase_mismatch() -> None:
    _, _, parent, child = make_phases()
    base = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=child
    )
    with pytest.raises(ValueError, match="parallel parent directions"):
        OrientationRelationship(
            name="mismatched",
            parent_phase=parent,
            child_phase=child,
            parent_to_child_rotation=base.parent_to_child_rotation,
            parallel_directions=(
                (
                    CrystalDirection([1.0, 1.0, 0.0], phase=child),
                    CrystalDirection([1.0, 0.0, 0.0], phase=child),
                ),
            ),
        )


def test_burgers_parallel_directions_keep_miller_bravais_meaning() -> None:
    _crystal_parent, crystal_child, parent, _ = make_phases()
    hex_symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal_child)
    hex_lattice = Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=crystal_child)
    alpha = Phase(
        "alpha", lattice=hex_lattice, symmetry=hex_symmetry, crystal_frame=crystal_child
    )
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=alpha
    )
    _, child_direction = relationship.parallel_directions[0]
    assert_allclose(child_direction.coordinates, [1.0, 1.0, 0.0], atol=1e-12)


def test_ks_defining_parallelisms_map_exactly() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    plane = ks.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    )
    assert_allclose(plane.rational_indices, [0, 1, 1])
    assert plane.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    direction = ks.map_direction_to_child(CrystalDirection([-1.0, 0.0, 1.0], phase=parent))
    assert_allclose(direction.rational_indices, [-1, -1, 1])
    assert direction.angular_residual_deg == pytest.approx(0.0, abs=1e-9)


def test_bain_direction_correspondence_matches_literature() -> None:
    _, _, parent, child = make_phases()
    bain = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=child
    )
    mapped = bain.map_direction_to_child(CrystalDirection([1.0, 1.0, 0.0], phase=parent))
    assert_allclose(mapped.rational_indices, [1, 0, 0])
    assert mapped.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    cube_axis = bain.map_direction_to_child(CrystalDirection([1.0, 0.0, 0.0], phase=parent))
    assert tuple(sorted(np.abs(cube_axis.rational_indices))) == (0, 1, 1)
    assert cube_axis.angular_residual_deg == pytest.approx(0.0, abs=1e-9)


def test_correspondence_reciprocal_is_inverse_transpose_and_preserves_zone_law() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    direct = ks.correspondence_direct()
    reciprocal = ks.correspondence_reciprocal()
    assert_allclose(reciprocal, np.linalg.inv(direct).T, atol=1e-12)
    plane_indices = np.array([1.0, 1.0, 1.0])
    direction_indices = np.array([-1.0, 0.0, 1.0])
    assert float((reciprocal @ plane_indices) @ (direct @ direction_indices)) == pytest.approx(
        0.0, abs=1e-12
    )


def test_direction_round_trip_recovers_indices() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    forward = ks.map_direction_to_child(CrystalDirection([-1.0, 0.0, 1.0], phase=parent))
    back = ks.map_direction_to_parent(forward.target)
    assert_allclose(back.rational_indices, [-1, 0, 1])
    assert back.angular_residual_deg == pytest.approx(0.0, abs=1e-9)


def test_ks_close_packed_group_membership_across_variants() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    plane = CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    results = map_plane_across_variants(ks, plane)
    assert len(results) == 24
    exact_images = [
        result
        for result in results
        if result.angular_residual_deg < 1e-6
        and tuple(sorted(np.abs(result.rational_indices))) == (0, 1, 1)
    ]
    # Exactly the six variants sharing (111) as their close-packed plane map
    # it onto a {011} child plane; the other 18 land on irrational images.
    assert len(exact_images) == 6
    assert all(result.variant_index is not None for result in results)


def test_burgers_hexagonal_correspondence_keeps_index_meaning() -> None:
    _crystal_parent, _, parent, _ = make_phases()
    hex_frame = ReferenceFrame(
        name="hex_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    alpha = Phase(
        "alpha",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=hex_frame),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=hex_frame),
        crystal_frame=hex_frame,
    )
    burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=alpha
    )
    basal = burgers.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([1, 1, 0]), phase=parent), phase=parent)
    )
    assert_allclose(basal.rational_indices, [0, 0, 1])
    assert basal.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    close_packed = burgers.map_direction_to_child(
        CrystalDirection([-1.0, 1.0, 1.0], phase=parent)
    )
    assert_allclose(close_packed.rational_indices, [1, 1, 0])
    assert close_packed.angular_residual_deg == pytest.approx(0.0, abs=1e-9)


def test_index_mapping_validates_phase_and_variant_membership() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    with pytest.raises(ValueError, match="parent phase"):
        ks.map_direction_to_child(CrystalDirection([1.0, 0.0, 0.0], phase=child))
    with pytest.raises(ValueError, match="child phase"):
        ks.map_direction_to_parent(CrystalDirection([1.0, 0.0, 0.0], phase=parent))
    foreign_variant = nw.generate_variants()[0]
    with pytest.raises(ValueError, match="must belong to this OrientationRelationship"):
        ks.map_direction_to_child(
            CrystalDirection([1.0, 0.0, 0.0], phase=parent), variant=foreign_variant
        )


def _paired_sets_for_deviation(
    parent: Phase, child: Phase
) -> tuple[OrientationSet, OrientationSet, OrientationRelationship]:
    specimen = ReferenceFrame(
        name="specimen_dev",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    parents = OrientationSet.from_orientations(
        [
            Orientation.from_euler(
                20.0, 30.0, 40.0, specimen_frame=specimen, symmetry=parent.symmetry, phase=parent
            ),
            Orientation.from_euler(
                70.0, 15.0, 5.0, specimen_frame=specimen, symmetry=parent.symmetry, phase=parent
            ),
        ]
    )
    variants = gt.generate_variants()
    picked = (variants[0], variants[7])
    quaternions = np.stack(
        [
            parents[index].rotation.compose(variant.parent_to_child_rotation.inverse()).quaternion
            for index, variant in enumerate(picked)
        ],
        axis=0,
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=child.crystal_frame,
        specimen_frame=specimen,
        symmetry=child.symmetry,
        phase=child,
    )
    return parents, children, gt


def test_named_or_misorientation_representations_match_literature() -> None:
    _, _, parent, child = make_phases()
    cases = {
        "ks": (
            OrientationRelationship.from_kurdjumov_sachs_correspondence(
                parent_phase=parent, child_phase=child
            ),
            42.85,
            (0.9679, 0.1776, 0.1776),
        ),
        "nw": (
            OrientationRelationship.from_nishiyama_wassermann_correspondence(
                parent_phase=parent, child_phase=child
            ),
            45.99,
            (0.9761, 0.2007, 0.0831),
        ),
        "bain": (
            OrientationRelationship.from_bain_correspondence(
                parent_phase=parent, child_phase=child
            ),
            45.0,
            (1.0, 0.0, 0.0),
        ),
    }
    for _, (relationship, expected_angle, expected_axis) in cases.items():
        misorientation = relationship.misorientation()
        assert misorientation.angle_deg == pytest.approx(expected_angle, abs=0.01)
        axis = np.sort(np.abs(misorientation.rotation.axis))[::-1]
        assert_allclose(axis, expected_axis, atol=1e-3)


def test_or_deviation_zero_for_exact_children_and_recovers_variants() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _paired_sets_for_deviation(parent, child)
    report = or_deviation(parents, children, gt)
    assert_allclose(report.deviations_deg, [0.0, 0.0], atol=1e-8)
    assert tuple(report.best_variant_indices) == (1, 8)
    assert report.mean_deviation_deg == pytest.approx(0.0, abs=1e-8)
    assert report.relationship_name == "greninger_troiano"


def test_or_deviation_measures_known_or_separations() -> None:
    _, _, parent, child = make_phases()
    parents, children, _ = _paired_sets_for_deviation(parent, child)
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    ks_report = or_deviation(parents, children, ks)
    nw_report = or_deviation(parents, children, nw)
    # Greninger-Troiano sits 2.40 deg from Kurdjumov-Sachs and 2.86 deg from
    # Nishiyama-Wassermann; children generated with GT must reproduce those
    # separations as their minimal deviations.
    assert_allclose(ks_report.deviations_deg, [2.404, 2.404], atol=5e-3)
    assert_allclose(nw_report.deviations_deg, [2.861, 2.861], atol=5e-3)


def test_or_deviation_validates_pairing_and_phase_semantics() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _paired_sets_for_deviation(parent, child)
    single_child = OrientationSet(
        quaternions=children.quaternions[:1],
        crystal_frame=children.crystal_frame,
        specimen_frame=children.specimen_frame,
        symmetry=children.symmetry,
        phase=children.phase,
    )
    with pytest.raises(ValueError, match="paired"):
        or_deviation(parents, single_child, gt)
    with pytest.raises(ValueError, match="parent phase"):
        or_deviation(children, children, gt)


def test_find_parallel_planes_recovers_ks_close_packed_pairing() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = find_parallel_planes(
        ks,
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent),
        tolerance_deg=0.01,
    )
    assert len(report.matches) == 24
    assert {tuple(sorted(np.abs(m.child_indices))) for m in report.matches} == {(0, 1, 1)}
    per_variant: dict[int, int] = {}
    for match in report.matches:
        per_variant[match.variant_index] = per_variant.get(match.variant_index, 0) + 1
    # Every variant pairs exactly one {111} member with a {011} child plane.
    assert set(per_variant.values()) == {1}
    assert len(per_variant) == 24
    assert all(m.angular_deviation_deg < 1e-6 for m in report.matches)


def test_find_parallel_directions_recovers_ks_close_packed_directions() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = find_parallel_directions(
        ks, CrystalDirection([1.0, 1.0, 0.0], phase=parent), tolerance_deg=0.01
    )
    assert len(report.matches) == 24
    assert {tuple(sorted(np.abs(m.child_indices))) for m in report.matches} == {(1, 1, 1)}


def test_describe_surfaces_state_conventions_and_key_numbers() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    text = ks.describe()
    assert "kurdjumov_sachs" in text
    assert "(111) parent || (011) child" in text
    # Negative components force a separator: "[-101]" would be ambiguous
    # between [-1, 0, 1] and [-10, 1].
    assert "[-1 0 1] parent || [-1 -1 1] child" in text
    assert "42.85 deg" in text
    assert "24 crystallographically distinct" in text

    plane = ks.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    )
    plane_text = plane.describe()
    assert "(111)" in plane_text and "(011)" in plane_text
    assert "residual of 0.0000 deg" in plane_text

    direction = ks.map_direction_to_child(CrystalDirection([-1.0, 0.0, 1.0], phase=parent))
    direction_text = direction.describe()
    assert "[-1 0 1]" in direction_text and "[-1 -1 1]" in direction_text

    parents, children, gt = _paired_sets_for_deviation(parent, child)
    deviation_text = or_deviation(parents, children, gt).describe()
    assert "greninger_troiano" in deviation_text
    assert "mean 0.000 deg" in deviation_text

    parallels_text = find_parallel_planes(
        ks,
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent),
        tolerance_deg=0.01,
    ).describe()
    assert "24 match(es)" in parallels_text
    assert "variant 1:" in parallels_text


def test_variant_selection_and_reconstruction_describe() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _paired_sets_for_deviation(parent, child)
    record = PhaseTransformationRecord(
        name="gt_record",
        orientation_relationship=gt,
        parent_orientation=parents[0],
        child_orientations=children,
    )
    selection = select_variants(record)
    text = selection.describe()
    assert "2 child orientation(s)" in text
    assert "deg" in text
    reconstruction = reconstruct_parent_orientation(record, parents)
    recon_text = reconstruction.describe()
    assert "gt_record" in recon_text
    assert "best candidate index 0" in recon_text


def _random_gt_pairs(
    parent: Phase, child: Phase, *, count: int = 30, seed: int = 42
) -> tuple[OrientationSet, OrientationSet, OrientationRelationship]:
    specimen = ReferenceFrame(
        name="specimen_fit",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    rng = np.random.default_rng(seed)
    eulers = rng.uniform(0.0, 60.0, size=(count, 3))
    parents = OrientationSet.from_orientations(
        [
            Orientation.from_euler(
                *euler, specimen_frame=specimen, symmetry=parent.symmetry, phase=parent
            )
            for euler in eulers
        ]
    )
    variants = gt.generate_variants()
    picks = rng.integers(0, len(variants), size=count)
    quaternions = np.stack(
        [
            parents[index]
            .rotation.compose(variants[int(picks[index])].parent_to_child_rotation.inverse())
            .quaternion
            for index in range(count)
        ],
        axis=0,
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=child.crystal_frame,
        specimen_frame=specimen,
        symmetry=child.symmetry,
        phase=child,
    )
    return parents, children, gt


def _symmetry_reduced_distance_deg(
    left: OrientationRelationship, right: OrientationRelationship
) -> float:
    from pytex.core.transformation import _symmetry_reduced_angle_between_deg

    return _symmetry_reduced_angle_between_deg(
        left.parent_to_child_rotation.as_matrix(),
        right.parent_to_child_rotation.as_matrix(),
        child_operators=left.child_phase.symmetry.operators,
        parent_operators=left.parent_phase.symmetry.operators,
    )


def test_fit_recovers_exact_relationship_from_matching_nominal() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _random_gt_pairs(parent, child)
    report = fit_orientation_relationship(parents, children, gt)
    assert report.converged
    assert report.deviation_from_nominal_deg == pytest.approx(0.0, abs=1e-8)
    assert report.mean_residual_deg == pytest.approx(0.0, abs=1e-8)


def test_fit_recovers_gt_from_ks_nominal_start() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _random_gt_pairs(parent, child)
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = fit_orientation_relationship(parents, children, ks)
    assert report.converged
    # The fit must land on the true (GT) relationship even though the
    # starting nominal was KS, and must report the documented KS-GT distance.
    assert _symmetry_reduced_distance_deg(report.relationship, gt) == pytest.approx(
        0.0, abs=1e-6
    )
    assert report.deviation_from_nominal_deg == pytest.approx(2.404, abs=5e-3)
    assert report.mean_residual_deg == pytest.approx(0.0, abs=1e-8)
    assert report.relationship.name == "kurdjumov_sachs_fitted"


def test_fit_averages_noise_toward_true_relationship() -> None:
    _, _, parent, child = make_phases()
    parents, children, gt = _random_gt_pairs(parent, child, count=40)
    rng = np.random.default_rng(7)
    noisy = []
    for index in range(len(children)):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        perturbation = Rotation.from_axis_angle(axis, np.deg2rad(rng.normal(0.0, 0.5)))
        noisy.append(
            perturbation.compose(Rotation(quaternion=children.quaternions[index])).quaternion
        )
    children_noisy = OrientationSet(
        quaternions=np.stack(noisy, axis=0),
        crystal_frame=child.crystal_frame,
        specimen_frame=children.specimen_frame,
        symmetry=child.symmetry,
        phase=child,
    )
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = fit_orientation_relationship(parents, children_noisy, ks)
    # Averaging must pull the estimate well inside the noise level.
    assert _symmetry_reduced_distance_deg(report.relationship, gt) < 0.15
    assert 0.1 < report.mean_residual_deg < 1.0
    text = report.describe()
    assert "kurdjumov_sachs_fitted" in text
    assert "Per-pair residuals" in text


def test_fit_validates_pairing_and_phases() -> None:
    _, _, parent, child = make_phases()
    _parents, children, gt = _random_gt_pairs(parent, child, count=4)
    with pytest.raises(ValueError, match="parent phase"):
        fit_orientation_relationship(children, children, gt)


def test_child_composition_follows_canonical_crystal_to_specimen_convention() -> None:
    """Regression for the V @ P convention bug (development-guide Phase 12).

    Under the canonical crystal->specimen orientation convention, the child
    orientation is g_child = g_parent o V^T: corresponding parent and child
    crystal directions must then map to the SAME specimen direction, and
    variant selection must recover planted variants from canonically built
    children.
    """

    _, _, parent, child = make_phases()
    specimen = ReferenceFrame(
        name="specimen_conv",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    variants = ks.generate_variants()
    parent_orientation = Orientation.from_euler(
        20.0, 30.0, 40.0, specimen_frame=specimen, symmetry=parent.symmetry, phase=parent
    )
    # Physics: the defining parallel directions coincide in specimen space.
    child_rotation = parent_orientation.rotation.compose(
        variants[0].parent_to_child_rotation.inverse()
    )
    parent_direction = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    child_direction = variants[0].parent_to_child_rotation.as_matrix() @ parent_direction
    specimen_from_parent = parent_orientation.rotation.as_matrix() @ parent_direction
    specimen_from_child = child_rotation.as_matrix() @ child_direction
    assert_allclose(specimen_from_parent, specimen_from_child, atol=1e-12)
    # Workflow: planted variants are recovered from canonically built children.
    picks = (1, 6, 12)
    quaternions = np.stack(
        [
            parent_orientation.rotation.compose(
                variants[index - 1].parent_to_child_rotation.inverse()
            ).quaternion
            for index in picks
        ],
        axis=0,
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=child.crystal_frame,
        specimen_frame=specimen,
        symmetry=child.symmetry,
        phase=child,
    )
    record = PhaseTransformationRecord(
        name="convention_record",
        orientation_relationship=ks,
        parent_orientation=parent_orientation,
        child_orientations=children,
    )
    report = select_variants(record)
    assert tuple(report.variant_indices) == picks
    assert report.scores_deg.max() == pytest.approx(0.0, abs=1e-8)


def test_orientation_set_slicing_returns_metadata_preserving_subset() -> None:
    _, _, parent, child = make_phases()
    _parents, children, _ = _paired_sets_for_deviation(parent, child)
    sliced = children[:1]
    assert isinstance(sliced, OrientationSet)
    assert len(sliced) == 1
    assert sliced.phase is children.phase
    assert sliced.symmetry == children.symmetry
    assert sliced.specimen_frame == children.specimen_frame
    assert_allclose(sliced.quaternions, children.quaternions[:1])
    single = children[-1]
    assert isinstance(single, Orientation)
    assert_allclose(
        np.abs(float(single.rotation.quaternion @ children.quaternions[-1])), 1.0, atol=1e-12
    )


def test_bain_deformation_gradient_matches_textbook_strains() -> None:
    _, _, parent, child = make_phases()
    bain = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = bain.deformation_gradient()
    ratio = 3.2 / 3.6
    expected = sorted([np.sqrt(2.0) * ratio, np.sqrt(2.0) * ratio, ratio], reverse=True)
    assert_allclose(report.principal_stretches, expected, atol=1e-10)
    assert report.volume_ratio == pytest.approx(2.0 * ratio**3, abs=1e-10)
    # Bain IS the pure correspondence distortion: no residual rigid rotation.
    assert report.polar_rotation_deg == pytest.approx(0.0, abs=1e-8)
    assert_allclose(report.stretch_tensor, report.stretch_tensor.T, atol=1e-12)


def test_ks_class_relationships_share_bain_stretches_with_literature_rotations() -> None:
    _, _, parent, child = make_phases()
    bain = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=child
    ).deformation_gradient()
    cases = {
        "ks": (
            OrientationRelationship.from_kurdjumov_sachs_correspondence(
                parent_phase=parent, child_phase=child
            ),
            11.06,
        ),
        "nw": (
            OrientationRelationship.from_nishiyama_wassermann_correspondence(
                parent_phase=parent, child_phase=child
            ),
            9.74,
        ),
    }
    for _, (relationship, expected_rotation) in cases.items():
        report = relationship.deformation_gradient()
        # Same lattice distortion as Bain: KS-class ORs differ only rigidly.
        assert_allclose(report.principal_stretches, bain.principal_stretches, atol=1e-10)
        assert report.volume_ratio == pytest.approx(bain.volume_ratio, abs=1e-10)
        # The polar rotation is the literature rigid-body rotation from Bain.
        assert report.polar_rotation_deg == pytest.approx(expected_rotation, abs=0.01)
        # F = R_polar U reconstructs the gradient.
        polar = report.deformation_gradient @ np.linalg.inv(report.stretch_tensor)
        assert_allclose(
            polar @ report.stretch_tensor, report.deformation_gradient, atol=1e-10
        )
        assert_allclose(polar @ polar.T, np.eye(3), atol=1e-10)
    text = cases["ks"][0].deformation_gradient().describe()
    assert "+25.71%" in text
    assert "-11.11%" in text
    assert "11.06 deg" in text


def test_shoji_nishiyama_correspondence_and_catalog() -> None:
    from pytex.core import standard_fcc_hcp_relationships, variant_close_packed_groups

    _crystal_parent, _, parent, _ = make_phases()
    hex_frame = ReferenceFrame(
        name="sn_hex_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    epsilon = Phase(
        "epsilon",
        lattice=Lattice(2.53, 2.53, 4.13, 90.0, 90.0, 120.0, crystal_frame=hex_frame),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=hex_frame),
        crystal_frame=hex_frame,
    )
    sn = OrientationRelationship.from_shoji_nishiyama_correspondence(
        parent_phase=parent, child_phase=epsilon
    )
    # Literature variant count: one variant per {111} close-packed plane.
    variants = sn.generate_variants()
    assert len(variants) == 4
    # Defining parallelisms map exactly: (111)_fcc -> (0001)_hcp, i.e. (001).
    basal = sn.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    )
    assert_allclose(basal.rational_indices, [0, 0, 1])
    assert basal.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    # Each variant descends from its own {111} plane: four packets of one.
    packets = variant_close_packed_groups(
        sn, CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    )
    assert np.bincount(packets).tolist() == [1, 1, 1, 1]
    catalog = standard_fcc_hcp_relationships(parent_phase=parent, child_phase=epsilon)
    assert catalog.names() == ("shoji_nishiyama",)
    with pytest.raises(ValueError, match="hexagonal child phase"):
        OrientationRelationship.from_shoji_nishiyama_correspondence(
            parent_phase=parent, child_phase=parent
        )


def test_pitsch_schrader_correspondence_and_hcp_bcc_catalog() -> None:
    from pytex.core import standard_hcp_bcc_relationships
    from pytex.core.transformation import _symmetry_reduced_angle_between_deg

    _, _, _, cubic_child = make_phases()
    hex_frame = ReferenceFrame(
        name="ps_hex_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    alpha = Phase(
        "alpha_ti",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=hex_frame),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=hex_frame),
        crystal_frame=hex_frame,
    )
    ps = OrientationRelationship.from_pitsch_schrader_correspondence(
        parent_phase=alpha, child_phase=cubic_child
    )
    # Three distinct variants: one per <11-20> axis of the hexagonal parent
    # (internally derived orbit count; 12 proper parent operators over the
    # order-4 stabilizer of the pairing).
    assert len(ps.generate_variants()) == 3
    # Defining parallelism maps exactly: basal (001)_hcp -> (110)_bcc.
    basal = CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=alpha)
    mapped = ps.map_plane_to_child(basal)
    assert tuple(sorted(np.abs(mapped.rational_indices))) == (0, 1, 1)
    assert mapped.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    # Literature separation: PS sits 5.26 deg from the inverse Burgers OR,
    # the hexagonal analogue of the KS-Pitsch separation.
    inverse_burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=cubic_child, child_phase=alpha
    ).inverse()
    separation = _symmetry_reduced_angle_between_deg(
        ps.parent_to_child_rotation.as_matrix(),
        inverse_burgers.parent_to_child_rotation.as_matrix(),
        child_operators=cubic_child.symmetry.operators,
        parent_operators=alpha.symmetry.operators,
    )
    assert separation == pytest.approx(5.26, abs=0.01)
    catalog = standard_hcp_bcc_relationships(parent_phase=alpha, child_phase=cubic_child)
    assert catalog.names() == ("pitsch_schrader", "burgers_inverse", "potter")
    with pytest.raises(ValueError, match="hexagonal parent phase"):
        OrientationRelationship.from_pitsch_schrader_correspondence(
            parent_phase=cubic_child, child_phase=cubic_child
        )


def test_potter_correspondence_pins_pyramidal_parallelism_and_burgers_proximity() -> None:
    from pytex.core.transformation import _symmetry_reduced_angle_between_deg

    _, _, _, cubic_child = make_phases()
    hex_frame = ReferenceFrame(
        name="potter_hex_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    alpha = Phase(
        "alpha_ti",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=hex_frame),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=hex_frame),
        crystal_frame=hex_frame,
    )
    potter = OrientationRelationship.from_potter_correspondence(
        parent_phase=alpha, child_phase=cubic_child
    )
    # Twelve distinct variants (internally derived orbit count: trivial
    # stabilizer of the pyramidal-plane/close-packed-direction pairing).
    assert len(potter.generate_variants()) == 12
    # Defining parallelism maps exactly: (01-11)_hcp -> {110}_bcc.
    pyramidal = CrystalPlane.from_miller_bravais((0, 1, -1, 1), phase=alpha)
    mapped = potter.map_plane_to_child(pyramidal)
    assert tuple(sorted(np.abs(mapped.rational_indices))) == (0, 1, 1)
    assert mapped.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    # The basal plane is NOT exactly parallel to {110}: it sits the small
    # Potter rotation away from its Burgers partner. The residual of the
    # basal-plane image equals the symmetry-reduced separation from the
    # inverse Burgers OR (the "rotation of Burgers about the shared
    # close-packed direction" structure; literature quotes ~2 deg, the exact
    # value is c/a-dependent — 1.370 deg at c/a = 4.68/2.95).
    basal = CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=alpha)
    basal_mapped = potter.map_plane_to_child(basal)
    inverse_burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=cubic_child, child_phase=alpha
    ).inverse()
    separation = _symmetry_reduced_angle_between_deg(
        potter.parent_to_child_rotation.as_matrix(),
        inverse_burgers.parent_to_child_rotation.as_matrix(),
        child_operators=cubic_child.symmetry.operators,
        parent_operators=alpha.symmetry.operators,
    )
    assert separation == pytest.approx(1.370, abs=0.005)
    assert basal_mapped.angular_residual_deg == pytest.approx(separation, abs=1e-6)
    assert 0.5 < separation < 3.0
    # The shared close-packed direction maps exactly, as in Burgers.
    direction = CrystalDirection.from_miller_bravais((2, -1, -1, 0), phase=alpha)
    mapped_direction = potter.map_direction_to_child(direction)
    assert tuple(sorted(np.abs(mapped_direction.rational_indices))) == (1, 1, 1)
    assert mapped_direction.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError, match="hexagonal parent phase"):
        OrientationRelationship.from_potter_correspondence(
            parent_phase=cubic_child, child_phase=cubic_child
        )


def _make_ferrite_and_cementite() -> tuple[Phase, Phase]:
    ferrite_frame = ReferenceFrame(
        name="ferrite_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    cementite_frame = ReferenceFrame(
        name="cementite_crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    # Pnma (Lipson-Petch) cementite setting, b > a > c; parameters as used by
    # Bhadeshia, Mater. Sci. Technol. 34 (2018) 1666.
    ferrite = Phase(
        "ferrite",
        lattice=Lattice(2.8662, 2.8662, 2.8662, 90.0, 90.0, 90.0, crystal_frame=ferrite_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=ferrite_frame),
        crystal_frame=ferrite_frame,
    )
    cementite = Phase(
        "cementite",
        lattice=Lattice(5.0883, 6.7416, 4.5241, 90.0, 90.0, 90.0, crystal_frame=cementite_frame),
        symmetry=SymmetrySpec.from_point_group("mmm", reference_frame=cementite_frame),
        crystal_frame=cementite_frame,
    )
    return ferrite, cementite


def test_bagaryatsky_correspondence_pins_all_three_axis_parallelisms() -> None:
    from pytex.core.transformation import _crystal_direction

    ferrite, cementite = _make_ferrite_and_cementite()
    bag = OrientationRelationship.from_bagaryatsky_correspondence(
        parent_phase=ferrite, child_phase=cementite
    )
    # Twelve distinct variants (internally derived orbit count: the parent
    # 180-deg rotation about [0-11] pairs with the child 180-deg rotation
    # about c to stabilize the correspondence).
    assert len(bag.generate_variants()) == 12
    # All three Bagaryatsky axis parallelisms map exactly (Pnma setting):
    # [1-1-1]_a -> [100]_th, [211]_a -> [010]_th, (0-11)_a -> (001)_th.
    mapped_a = bag.map_direction_to_child(_crystal_direction((1, -1, -1), phase=ferrite))
    assert_allclose(mapped_a.rational_indices, [1, 0, 0])
    assert mapped_a.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    mapped_b = bag.map_direction_to_child(_crystal_direction((2, 1, 1), phase=ferrite))
    assert_allclose(mapped_b.rational_indices, [0, 1, 0])
    assert mapped_b.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    mapped_c = bag.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([0, -1, 1]), phase=ferrite), phase=ferrite)
    )
    assert_allclose(mapped_c.rational_indices, [0, 0, 1])
    assert mapped_c.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError, match="orthorhombic child phase"):
        OrientationRelationship.from_bagaryatsky_correspondence(
            parent_phase=ferrite, child_phase=ferrite
        )


def test_isaichev_correspondence_and_ferrite_cementite_catalog() -> None:
    from pytex.core import standard_ferrite_cementite_relationships
    from pytex.core.transformation import (
        _crystal_direction,
        _symmetry_reduced_angle_between_deg,
    )

    ferrite, cementite = _make_ferrite_and_cementite()
    isa = OrientationRelationship.from_isaichev_correspondence(
        parent_phase=ferrite, child_phase=cementite
    )
    # The irrational (031) pairing breaks the Bagaryatsky stabilizer:
    # twenty-four distinct variants (internally derived orbit count).
    assert len(isa.generate_variants()) == 24
    # Defining parallelisms exact: (101)_a -> (031)_th, [1-1-1]_a -> [100]_th.
    mapped_plane = isa.map_plane_to_child(
        CrystalPlane(MillerIndex(np.array([1, 0, 1]), phase=ferrite), phase=ferrite)
    )
    assert_allclose(mapped_plane.rational_indices, [0, 3, 1])
    assert mapped_plane.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    mapped_direction = isa.map_direction_to_child(
        _crystal_direction((1, -1, -1), phase=ferrite)
    )
    assert_allclose(mapped_direction.rational_indices, [1, 0, 0])
    assert mapped_direction.angular_residual_deg == pytest.approx(0.0, abs=1e-9)
    # Isaichev is the Bagaryatsky orientation rotated about the cementite
    # a-axis; the magnitude depends on the cementite axial ratios (3.586 deg
    # for these lattice parameters; the modern literature quotes ~3.8 deg).
    bag = OrientationRelationship.from_bagaryatsky_correspondence(
        parent_phase=ferrite, child_phase=cementite
    )
    separation = _symmetry_reduced_angle_between_deg(
        bag.parent_to_child_rotation.as_matrix(),
        isa.parent_to_child_rotation.as_matrix(),
        child_operators=cementite.symmetry.operators,
        parent_operators=ferrite.symmetry.operators,
    )
    assert separation == pytest.approx(3.586, abs=0.005)
    relative = (
        isa.parent_to_child_rotation.as_matrix() @ bag.parent_to_child_rotation.as_matrix().T
    )
    eigenvalues, eigenvectors = np.linalg.eig(relative)
    axis = np.real(eigenvectors[:, np.argmin(np.abs(eigenvalues - 1.0))])
    assert_allclose(np.abs(axis), [1.0, 0.0, 0.0], atol=1e-9)
    catalog = standard_ferrite_cementite_relationships(
        parent_phase=ferrite, child_phase=cementite
    )
    assert catalog.names() == ("bagaryatsky", "isaichev")
    with pytest.raises(ValueError, match="cubic parent phase"):
        OrientationRelationship.from_isaichev_correspondence(
            parent_phase=cementite, child_phase=cementite
        )


# ---------------------------------------------------------------------------
# Reconstructive transformations: the Burgers shuffle
# ---------------------------------------------------------------------------


def _zirconium_allotropes() -> tuple[Phase, Phase]:
    """Beta (bcc) and alpha (hcp) zirconium, the classic Burgers pair.

    Lattice parameters: beta-Zr a = 3.574 A (high-temperature bcc allotrope) and
    alpha-Zr a = 3.232 A, c = 5.147 A, matching the ``zr_hcp`` fixture.
    """

    # The alpha phase comes from its CIF fixture, which pymatgen parses.

    from pytex import crystal_frame, get_phase_fixture

    alpha_frame = crystal_frame("alpha_zr")
    beta_frame = crystal_frame("beta_zr")
    alpha = get_phase_fixture("zr_hcp").load_phase(crystal_frame=alpha_frame)
    beta = Phase(
        "beta-zirconium",
        lattice=Lattice(3.574, 3.574, 3.574, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
        crystal_frame=beta_frame,
    )
    return beta, alpha


def test_burgers_relationship_has_twelve_variants_and_the_literature_misorientation() -> None:
    beta, alpha = _zirconium_allotropes()
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    )
    assert len(relationship.generate_variants()) == 12
    assert relationship.misorientation().angle_deg == pytest.approx(45.29, abs=0.02)


def test_burgers_correspondence_is_half_integer_because_of_the_shuffle() -> None:
    """bcc-to-hcp has no integer Bravais-lattice correspondence.

    The bcc primitive cell has half the volume of the hcp one, so two bcc lattice
    points map onto one hcp lattice point plus one motif atom. The missing half is
    the Burgers shuffle, and it shows up as a half-integer correspondence. A
    denominator of 2 is therefore a scientific statement: this transformation
    cannot be a pure lattice strain.
    """

    beta, alpha = _zirconium_allotropes()
    report = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    ).deformation_gradient()

    assert report.correspondence_denominator == 2
    doubled = np.asarray(report.correspondence) * 2.0
    assert np.allclose(doubled, np.rint(doubled), atol=1e-9)
    assert not np.allclose(
        report.correspondence, np.rint(report.correspondence), atol=1e-9
    )
    # A unit-determinant correspondence: one parent cell of atoms per child cell.
    assert float(np.linalg.det(report.correspondence)) == pytest.approx(1.0, abs=1e-9)


def test_burgers_transformation_strain_matches_an_independent_volume_calculation() -> None:
    """The volume change is checked against the cell volumes, not a prior output.

    beta-Zr conventional cell (2 atoms) has volume a^3; alpha-Zr conventional cell
    (2 atoms) has volume (sqrt(3)/2) a^2 c. Their ratio is an independent
    derivation of ``det F``, so agreement validates the whole correspondence and
    deformation-gradient path rather than restating it.
    """

    beta, alpha = _zirconium_allotropes()
    report = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    ).deformation_gradient()

    beta_volume = beta.lattice.a**3
    alpha_volume = (np.sqrt(3.0) / 2.0) * alpha.lattice.a**2 * alpha.lattice.c
    assert report.volume_ratio == pytest.approx(alpha_volume / beta_volume, rel=1e-9)
    # About +2% for zirconium.
    assert report.volume_ratio == pytest.approx(1.0199, abs=1e-3)


def test_burgers_principal_strains_are_the_classic_ten_percent_pair() -> None:
    """One near +10% extension, one near -10% contraction, one small.

    This is the shape of the Burgers lattice strain reported for the beta-to-alpha
    transformation in titanium and zirconium.
    """

    beta, alpha = _zirconium_allotropes()
    report = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    ).deformation_gradient()

    strains = np.sort((np.asarray(report.principal_stretches) - 1.0) * 100.0)[::-1]
    assert strains[0] == pytest.approx(10.76, abs=0.05)
    assert strains[1] == pytest.approx(1.83, abs=0.05)
    assert strains[2] == pytest.approx(-9.57, abs=0.05)


def test_cubic_relationships_still_use_an_integer_correspondence() -> None:
    """Widening the denominator search must not change the cubic cases."""

    # Both phases come from CIF fixtures, which is pymatgen's job.

    from pytex import crystal_frame, get_phase_fixture

    parent = get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal_frame("fcc"))
    child = get_phase_fixture("fe_bcc").load_phase(crystal_frame=crystal_frame("bcc"))
    for builder in (
        OrientationRelationship.from_bain_correspondence,
        OrientationRelationship.from_kurdjumov_sachs_correspondence,
        OrientationRelationship.from_nishiyama_wassermann_correspondence,
    ):
        report = builder(parent_phase=parent, child_phase=child).deformation_gradient()
        assert report.correspondence_denominator == 1
        assert np.allclose(report.correspondence, np.rint(report.correspondence), atol=1e-9)
        # The textbook Bain stretches, shared by every KS/NW variant.
        assert np.sort(report.principal_stretches)[::-1] == pytest.approx(
            [1.1504, 1.1504, 0.8135], abs=1e-3
        )


def test_deformation_gradient_accepts_an_explicit_correspondence() -> None:
    """A literature correspondence can be pinned instead of searched for."""

    beta, alpha = _zirconium_allotropes()
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    )
    searched = relationship.deformation_gradient()
    pinned = relationship.deformation_gradient(correspondence=searched.correspondence)
    assert pinned.correspondence_denominator == 2
    assert np.allclose(pinned.deformation_gradient, searched.deformation_gradient, atol=1e-12)

    with pytest.raises(ValueError, match="singular"):
        relationship.deformation_gradient(correspondence=np.zeros((3, 3)))
    with pytest.raises(ValueError, match=r"shape \(3, 3\)"):
        relationship.deformation_gradient(correspondence=np.eye(2))


def _reverse_burgers() -> OrientationRelationship:
    """The hcp-to-bcc relationship, built explicitly.

    ``from_burgers_correspondence`` cannot be reused with the arguments swapped:
    it requires a cubic parent (432) and hexagonal child (622). The reverse is a
    different object, with the phases exchanged.
    """

    beta, alpha = _zirconium_allotropes()
    return OrientationRelationship.from_parallel_plane_direction(
        name="burgers_reverse",
        parent_plane=CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=alpha),
        child_plane=CrystalPlane(MillerIndex(np.array([1, 1, 0]), phase=beta), phase=beta),
        parent_direction=CrystalDirection.from_miller_bravais((1, 1, -2, 0), phase=alpha),
        child_direction=CrystalDirection([-1.0, 1.0, 1.0], phase=beta),
    )


def test_reverse_burgers_has_six_variants_and_the_same_misorientation() -> None:
    """The reverse is not the inverse rotation: the symmetry reduction differs.

    Forward bcc-to-hcp gives 12 variants because the cubic parent's 24 proper
    operators reduce against the hexagonal child's 12. Reversed, the roles swap
    and the orbit reduces to 6. The misorientation angle is unchanged, because it
    is the same relationship read in the other direction.
    """

    reverse = _reverse_burgers()
    assert len(reverse.generate_variants()) == 6
    assert reverse.misorientation().angle_deg == pytest.approx(45.29, abs=0.02)


def test_reverse_burgers_strain_undoes_the_forward_volume_change() -> None:
    """hcp-to-bcc must contract by what bcc-to-hcp expanded."""

    beta, alpha = _zirconium_allotropes()
    forward = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    ).deformation_gradient()
    reverse = _reverse_burgers().deformation_gradient()

    assert reverse.correspondence_denominator == 2
    assert reverse.volume_ratio == pytest.approx(1.0 / forward.volume_ratio, rel=1e-6)
    # Independently: the beta cell volume over the alpha cell volume.
    beta_volume = beta.lattice.a**3
    alpha_volume = (np.sqrt(3.0) / 2.0) * alpha.lattice.a**2 * alpha.lattice.c
    assert reverse.volume_ratio == pytest.approx(beta_volume / alpha_volume, rel=1e-6)


def test_correspondence_search_rejects_an_invertible_but_poor_integer_fit() -> None:
    """Fit quality, not mere invertibility, picks the denominator.

    For the reverse relationship the denominator-1 rounding *is* invertible, but
    it is a poor fit with determinant 2 — taking it would report a doubled cell
    and a nonsensical +96% volume change. The search requires the fit to be within
    the rationalization tolerance before accepting a denominator.
    """

    reverse = _reverse_burgers()
    report = reverse.deformation_gradient()
    assert report.correspondence_denominator == 2
    assert report.correspondence_max_component_error <= 0.25
    assert float(np.linalg.det(report.correspondence)) == pytest.approx(1.0, abs=1e-9)
    # The volume change is a contraction of about 2%, not an expansion of 96%.
    assert -0.03 < report.volume_ratio - 1.0 < 0.0


def test_a_large_contraction_does_not_tip_a_cubic_pair_to_a_finer_grid() -> None:
    """The coarsest grid wins unless a finer one fits substantially better.

    For a 3.60 to 2.87 Bain pair the denominator-2 rounding fits marginally
    better than denominator 1 — by 0.009 — but with determinant 3 instead of 2,
    which is not the physical correspondence. Preferring the coarser grid unless
    the improvement is substantial keeps this case integer.
    """

    _, _, parent, child = make_phases()
    contracted = Phase(
        "martensite_contracted",
        lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child.crystal_frame),
        symmetry=child.symmetry,
        crystal_frame=child.crystal_frame,
    )
    report = OrientationRelationship.from_bain_correspondence(
        parent_phase=parent, child_phase=contracted
    ).deformation_gradient()

    assert report.correspondence_denominator == 1
    assert float(np.linalg.det(report.correspondence)) == pytest.approx(2.0, abs=1e-9)
    # A contraction of a=3.60 to a=2.87 is a small volume change, not +52%.
    assert abs(report.volume_ratio - 1.0) < 0.1


def _random_rotation_matrices(count: int, rng: np.random.Generator) -> np.ndarray:
    quaternions = rng.normal(size=(count, 4))
    quaternions /= np.linalg.norm(quaternions, axis=1)[:, None]
    return np.stack(
        [Rotation(quaternion=quaternion).as_matrix() for quaternion in quaternions], axis=0
    )


def test_boundary_fingerprint_contains_the_identity_and_is_symmetry_closed() -> None:
    """Group-theoretic identities of the same-parent boundary set.

    Two children of one parent formed through the *same* variant have zero
    misorientation, so the identity must belong to the set; and because each
    child orientation is only defined up to its own crystal symmetry, the set
    must be invariant under left and right multiplication by child operators.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
    )

    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)
    assert fingerprint.ndim == 3
    assert fingerprint.shape[1:] == (3, 3)
    assert_allclose(np.linalg.det(fingerprint), 1.0, atol=1e-9)
    # The identity (same-variant pair) is present: trace 3.
    assert float(np.einsum("kii->k", fingerprint).max()) == pytest.approx(3.0, abs=1e-9)

    # Closure is asserted as a set-distance statement rather than by comparing
    # rounded quaternion keys: the keys are brittle at rounding boundaries, and
    # the mathematical claim is that the transformed set lands back on the
    # original set, which is exactly "every transformed element is at zero
    # distance from the fingerprint".
    operators = child.symmetry.operators
    for operator in (operators[1], operators[5]):
        for transformed in (
            np.einsum("ij,kjl->kil", operator, fingerprint),
            np.einsum("kij,jl->kil", fingerprint, operator),
        ):
            distances = boundary_fingerprint_distances_deg(transformed, fingerprint)
            # 1e-5 deg, not 0: arccos near zero and the quaternion/matrix round
            # trip both floor at ~1e-6 deg (see the index-correspondence notes).
            assert float(distances.max()) < 1e-5


def test_boundary_fingerprint_is_exact_on_same_parent_boundaries() -> None:
    """Every variant pair of one parent scores zero against the fingerprint.

    This is the defining property: the admissible set is exactly the
    ``V_i V_j^T`` family, so a boundary built from any two variants of a
    common parent has zero distance to it, up to the quaternion/matrix
    round-trip floor.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
    )

    _, _, parent, child = make_phases()
    specimen = ReferenceFrame(
        name="fingerprint_specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    variants = relationship.generate_variants()
    parent_orientation = Orientation.from_euler(
        20.0, 30.0, 40.0, specimen_frame=specimen, symmetry=parent.symmetry, phase=parent
    )
    children = np.stack(
        [
            parent_orientation.rotation.compose(
                variant.parent_to_child_rotation.inverse()
            ).as_matrix()
            for variant in variants
        ],
        axis=0,
    )
    left, right = np.triu_indices(len(variants), k=1)
    relative = np.einsum("nji,njk->nik", children[left], children[right], optimize=True)
    distances = boundary_fingerprint_distances_deg(
        relative, intervariant_boundary_fingerprint(relationship)
    )
    assert distances.shape == (left.size,)
    assert float(distances.max()) < 1e-5


def test_boundary_fingerprint_rejects_unrelated_boundaries_far_better_than_angles() -> None:
    """The axis carries most of the discriminating power of the fingerprint.

    Matching only the misorientation *angle* against the intervariant spectrum
    is very permissive: for a cubic-cubic relationship those angles are spread
    densely enough that a few-degree window admits a large fraction of
    entirely unrelated boundaries. Matching the full rotation removes most of
    them. The rates below are measured properties of the method on uniformly
    random rotations with a pinned seed, not literature values, so they are
    asserted with wide margins that pin the qualitative separation.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
        intervariant_misorientation_angles_deg,
    )
    from pytex.core.orientation import _reduced_pair_disorientation_angles

    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    rng = np.random.default_rng(20260730)
    left = _random_rotation_matrices(4000, rng)
    right = _random_rotation_matrices(4000, rng)
    relative = np.einsum("eji,ejk->eik", left, right, optimize=True)

    table = intervariant_misorientation_angles_deg(relationship)
    spectrum = np.unique(
        np.round(np.concatenate([[0.0], table[np.triu_indices(table.shape[0], k=1)]]), 6)
    )
    operators = child.symmetry.operators
    angles = np.degrees(_reduced_pair_disorientation_angles(relative, operators, operators))
    angle_only = np.min(np.abs(angles[:, None] - spectrum[None, :]), axis=1)
    full = boundary_fingerprint_distances_deg(
        relative, intervariant_boundary_fingerprint(relationship)
    )

    angle_rate = float(np.mean(angle_only <= 3.0))
    full_rate = float(np.mean(full <= 3.0))
    assert angle_rate > 0.40, f"angle-only false-accept fell to {angle_rate:.3f}"
    assert full_rate < 0.15, f"fingerprint false-accept rose to {full_rate:.3f}"
    assert angle_rate > 4.0 * full_rate
    # Tightening the tolerance sharpens the fingerprint test quickly.
    assert float(np.mean(full <= 1.0)) < 0.02


def _hexagonal_child(reference: Phase, *, c_over_a: float = 1.587) -> Phase:
    """An hcp child sharing a cubic phase's crystal frame, for group-theory tests."""

    lattice = Lattice(
        2.95,
        2.95,
        2.95 * c_over_a,
        90.0,
        90.0,
        120.0,
        crystal_frame=reference.crystal_frame,
    )
    return Phase(
        "alpha",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(
            "6/mmm", reference_frame=reference.crystal_frame
        ),
        crystal_frame=reference.crystal_frame,
    )


def test_intervariant_fingerprint_sizes_are_pinned() -> None:
    """The admissible same-parent set is much smaller for Burgers than for KS.

    The size of $G_c (R G_p R^T) G_c$ governs how much of orientation space a
    chance boundary can land in, and therefore how reconstructable a
    relationship is from boundary evidence alone. It is pure group theory — it
    depends on the point groups and the relationship, never on lattice
    parameters — so it is an exact integer, not a measurement, and the
    reconstruction robustness study reasons from it in prose.

    Both numbers this test pins were wrong in the documentation before it
    existed: the study quoted "about 2 800" for Burgers against "about 10 700"
    for Kurdjumov-Sachs. The Burgers estimate had treated the hexagonal child as
    contributing 24 proper operators rather than 12, and the Kurdjumov-Sachs
    figure came from a deduplication that double-counted 81 elements. A number a
    document argues from is a number a test should hold.
    """

    from pytex.core import intervariant_boundary_fingerprint

    _, _, cubic_parent, cubic_child = make_phases()
    kurdjumov_sachs = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=cubic_parent, child_phase=cubic_child
    )
    assert len(kurdjumov_sachs.generate_variants()) == 24
    assert intervariant_boundary_fingerprint(kurdjumov_sachs).shape[0] == 10584

    burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=cubic_parent, child_phase=_hexagonal_child(cubic_child)
    )
    assert len(burgers.generate_variants()) == 12
    assert intervariant_boundary_fingerprint(burgers).shape[0] == 684


def test_intervariant_fingerprint_size_is_independent_of_lattice_parameters() -> None:
    """A group-theoretic count must not move with the cell it was computed from.

    This is the property the previous deduplication broke. It canonicalized a
    quaternion's sign on its largest-magnitude component, and for the 90 and 180
    degree elements of a crystal point group two components tie in magnitude, so
    ``argmax`` broke the tie arbitrarily: two numerically identical rotations
    canonicalized to ``q`` and ``-q`` and were counted twice. It also rounded to
    a fixed number of decimals, so a value near a rounding boundary rounded
    either way depending on floating-point noise. Together those made the count
    of a fixed set depend on the lattice parameters that entered the rotation.
    """

    from pytex.core import intervariant_boundary_fingerprint

    _, _, parent, child = make_phases()
    cubic_counts = set()
    for parameter in (2.87, 3.2, 3.61):
        scaled = Phase(
            "child",
            lattice=Lattice(
                parameter,
                parameter,
                parameter,
                90.0,
                90.0,
                90.0,
                crystal_frame=child.crystal_frame,
            ),
            symmetry=child.symmetry,
            crystal_frame=child.crystal_frame,
        )
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=scaled
        )
        cubic_counts.add(intervariant_boundary_fingerprint(relationship).shape[0])
    assert cubic_counts == {10584}

    hexagonal_counts = set()
    for ratio in (1.587, 1.60, 1.633):
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=_hexagonal_child(child, c_over_a=ratio)
        )
        hexagonal_counts.add(intervariant_boundary_fingerprint(relationship).shape[0])
    assert hexagonal_counts == {684}


def test_intervariant_fingerprint_holds_no_duplicates() -> None:
    """No two elements of the returned set may name the same rotation.

    Duplicates were never a correctness problem downstream, because the distance
    kernel takes a maximum over the set. They mattered because the set's size is
    quoted as a scientific quantity, and because 81 wasted elements is 0.8% of
    every distance evaluation.
    """

    from pytex.core import intervariant_boundary_fingerprint

    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)
    traces = np.einsum("aij,bij->ab", fingerprint, fingerprint, optimize=True)
    np.fill_diagonal(traces, -3.0)
    # trace(A B^T) = 3 only when A and B are the same rotation.
    assert float(np.max(traces)) < 3.0 - 1e-6


def test_boundary_fingerprint_distance_is_blocked_but_exact() -> None:
    """The blocked kernel must equal the direct all-pairs evaluation.

    Edges are processed in fixed-size blocks so the transient stays bounded at
    map scale; that must not change the result, including when the edge count
    is not a multiple of the block size.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
    )
    from pytex.core.transformation import _FINGERPRINT_BLOCK_SIZE

    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent, child_phase=child
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)
    rng = np.random.default_rng(11)
    count = 2 * _FINGERPRINT_BLOCK_SIZE + 7
    relative = _random_rotation_matrices(count, rng)
    blocked = boundary_fingerprint_distances_deg(relative, fingerprint)
    traces = np.einsum("eij,kij->ek", relative, fingerprint, optimize=True).max(axis=1)
    direct = np.degrees(np.arccos(np.clip((traces - 1.0) * 0.5, -1.0, 1.0)))
    assert blocked.shape == (count,)
    assert_allclose(blocked, direct, atol=1e-9)


def test_boundary_fingerprint_distance_validates_shapes() -> None:
    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
    )

    _, _, parent, child = make_phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)
    with pytest.raises(ValueError, match="relative_matrices"):
        boundary_fingerprint_distances_deg(np.eye(3), fingerprint)
    with pytest.raises(ValueError, match="fingerprint must have shape"):
        boundary_fingerprint_distances_deg(np.eye(3)[None, :, :], np.eye(3))
    with pytest.raises(ValueError, match="at least one rotation"):
        boundary_fingerprint_distances_deg(np.eye(3)[None, :, :], np.empty((0, 3, 3)))


# --------------------------------------------------------------------------- #
# Per-variant parallelisms
# --------------------------------------------------------------------------- #


def test_variant_parallelisms_map_exactly_for_every_variant() -> None:
    # V = S_c R S_p^T maps THIS variant's parent objects onto its child objects;
    # the relationship's nominal pair is only variant 1's.
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    for variant in ks.generate_variants():
        rotation = variant.parent_to_child_rotation.as_matrix()
        for parent_plane, child_plane in variant.parallel_planes:
            assert_allclose(rotation @ parent_plane.normal, child_plane.normal, atol=1e-12)
        for parent_direction, child_direction in variant.parallel_directions:
            assert_allclose(
                rotation @ parent_direction.unit_vector,
                child_direction.unit_vector,
                atol=1e-12,
            )


def test_variant_parallel_planes_stay_in_the_defining_family() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    variants = ks.generate_variants()
    parent_members = {
        tuple(sorted(abs(int(value)) for value in variant.parallel_planes[0][0].miller.indices))
        for variant in variants
    }
    child_members = {
        tuple(sorted(abs(int(value)) for value in variant.parallel_planes[0][1].miller.indices))
        for variant in variants
    }
    assert parent_members == {(1, 1, 1)}  # every variant carries a {111}_gamma
    assert child_members == {(0, 1, 1)}  # onto a {011}_alpha
    # and the parent members are not all the same specific plane
    specific = {
        tuple(int(value) for value in variant.parallel_planes[0][0].miller.indices)
        for variant in variants
    }
    assert len(specific) > 1


def test_variant_parallelisms_agree_with_close_packed_grouping() -> None:
    # the packet a variant belongs to is exactly the parent family member its
    # own parallelism names, so the two independent routes must agree
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    variants = ks.generate_variants()
    from pytex.core import variant_close_packed_groups

    labels = variant_close_packed_groups(
        ks, CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)
    )
    grouped: dict[int, set[tuple[int, ...]]] = {}
    for variant, label in zip(variants, labels, strict=True):
        indices = np.asarray(variant.parallel_planes[0][0].miller.indices, dtype=np.int64)
        canonical = tuple(int(value) for value in (indices if indices[0] >= 0 else -indices))
        grouped.setdefault(int(label), set()).add(canonical)
    assert len(grouped) == 4
    assert all(len(members) == 1 for members in grouped.values())
    assert sum(int(np.sum(labels == label)) for label in set(labels.tolist())) == 24


def test_variant_parallelisms_hold_for_a_hexagonal_child() -> None:
    _, _, parent, _ = make_phases()
    hcp = make_hcp_child()
    burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=hcp
    )
    variants = burgers.generate_variants()
    assert len(variants) == 12
    for variant in variants:
        rotation = variant.parent_to_child_rotation.as_matrix()
        for parent_plane, child_plane in variant.parallel_planes:
            assert_allclose(rotation @ parent_plane.normal, child_plane.normal, atol=1e-12)


def test_variant_symmetry_operators_reproduce_the_variant_rotation() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    base = ks.parent_to_child_rotation.as_matrix()
    for variant in ks.generate_variants():
        expected = (
            variant.child_symmetry_operator @ base @ variant.parent_symmetry_operator.T
        )
        assert_allclose(variant.parent_to_child_rotation.as_matrix(), expected, atol=1e-12)


# --------------------------------------------------------------------------- #
# Rationalizing a measured relationship, and pricing the idealization
# --------------------------------------------------------------------------- #


def _measured_report(relationship, picks=(0, 4, 8, 13, 17, 22)):  # type: ignore[no-untyped-def]
    """Characterize exact children of one parent under `relationship`.

    Exact rather than noisy on purpose: with zero scatter the fit residual is
    zero, so every number the rationalizer reports is the cost of the
    idealization alone and nothing else.
    """

    from pytex.core import characterize_orientation_relationship, specimen_frame

    variants = relationship.generate_variants()
    chosen = [index for index in picks if index < len(variants)]
    parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
    # Canonical crystal->specimen convention: C = P V^T.
    children = np.stack(
        [parent_matrix @ variants[index].parent_to_child_rotation.as_matrix().T
         for index in chosen]
    )
    frame = specimen_frame()
    parents = OrientationSet.from_matrices(
        np.stack([parent_matrix] * len(chosen)),
        specimen_frame=frame,
        phase=relationship.parent_phase,
    )
    return characterize_orientation_relationship(
        parents,
        OrientationSet.from_matrices(
            children, specimen_frame=frame, phase=relationship.child_phase
        ),
    )


def test_rationalizing_exact_ks_costs_nothing_and_recovers_its_statement() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    result = _measured_report(ks).as_rational_relationship()
    assert result.residual_rotation_deg == pytest.approx(0.0, abs=1e-6)
    assert result.plane_statement.deviation_deg == pytest.approx(0.0, abs=1e-6)
    assert result.direction_statement.deviation_deg == pytest.approx(0.0, abs=1e-6)
    plane = tuple(sorted(abs(int(v)) for v in result.plane_statement.parent_indices))
    child_plane = tuple(sorted(abs(int(v)) for v in result.plane_statement.child_indices))
    direction = tuple(sorted(abs(int(v)) for v in result.direction_statement.parent_indices))
    child_direction = tuple(
        sorted(abs(int(v)) for v in result.direction_statement.child_indices)
    )
    assert plane == (1, 1, 1)  # {111}_gamma
    assert child_plane == (0, 1, 1)  # {011}_alpha
    assert direction == (0, 1, 1)  # <110>_gamma
    assert child_direction == (1, 1, 1)  # <111>_alpha


def test_rationalizing_greninger_troiano_to_low_indices_costs_the_ks_separation() -> None:
    """The number this whole surface exists to report.

    Greninger-Troiano has no low-index direction pair. Held to |index| <= 2 the
    tidiest integer statement available *is* the Kurdjumov-Sachs one, and the
    price of writing it is the published KS-GT separation of 2.40 deg — not
    zero, and not hidden. An idealization returned without that number would
    look like a measurement of Kurdjumov-Sachs.
    """

    _, _, parent, child = make_phases()
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = _measured_report(gt)
    assert report.mean_residual_deg == pytest.approx(0.0, abs=1e-6)
    result = report.as_rational_relationship(max_index=2)
    plane = tuple(sorted(abs(int(v)) for v in result.plane_statement.parent_indices))
    direction = tuple(sorted(abs(int(v)) for v in result.direction_statement.parent_indices))
    assert plane == (1, 1, 1)
    assert direction == (0, 1, 1)
    assert result.residual_rotation_deg == pytest.approx(2.404, abs=0.01)


def test_the_cost_of_the_idealization_falls_as_the_indices_are_allowed_to_grow() -> None:
    """Tidier integers cost more; the trade is the user's to make, so it has to
    be visible."""

    _, _, parent, child = make_phases()
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = _measured_report(gt)
    residuals = [
        report.as_rational_relationship(max_index=bound).residual_rotation_deg
        for bound in (2, 3, 4)
    ]
    assert residuals == sorted(residuals, reverse=True)
    assert residuals[0] > residuals[-1] + 1.0


def test_the_reported_direction_lies_in_the_reported_plane() -> None:
    """Without the zone-law filter the constructed relationship is not the one
    the two printed labels describe: `from_parallel_plane_direction` drops the
    direction's normal component."""

    _, _, parent, child = make_phases()
    for constructor in (
        OrientationRelationship.from_kurdjumov_sachs_correspondence,
        OrientationRelationship.from_nishiyama_wassermann_correspondence,
        OrientationRelationship.from_greninger_troiano_correspondence,
        OrientationRelationship.from_pitsch_correspondence,
    ):
        relationship = constructor(parent_phase=parent, child_phase=child)
        result = _measured_report(relationship).as_rational_relationship()
        plane = np.asarray(result.plane_statement.parent_indices, dtype=np.int64)
        direction = np.asarray(result.direction_statement.parent_indices, dtype=np.int64)
        assert int(np.dot(plane, direction)) == 0, relationship.name
        assert result.zone_law_deviation_deg == pytest.approx(0.0, abs=1e-9)


def test_the_idealized_relationship_realizes_its_own_statement_exactly() -> None:
    """The result is a genuine relationship, so its own parallelisms are exact.
    The deviations it carries are distances from the *measurement*, not from
    being a relationship."""

    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    idealized = _measured_report(ks).as_rational_relationship().relationship
    rotation = idealized.parent_to_child_rotation.as_matrix()
    for parent_plane, child_plane in idealized.parallel_planes:
        assert_allclose(rotation @ parent_plane.normal, child_plane.normal, atol=1e-12)
    for parent_direction, child_direction in idealized.parallel_directions:
        assert_allclose(
            rotation @ parent_direction.unit_vector, child_direction.unit_vector, atol=1e-12
        )


def test_the_idealization_is_named_so_it_cannot_be_mistaken_for_the_measurement() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = _measured_report(ks)
    result = report.as_rational_relationship()
    assert result.relationship.name.endswith("_rationalized")
    assert result.source_relationship_name == report.relationship.name
    assert report.as_rational_relationship(name="ks_ideal").relationship.name == "ks_ideal"


def test_describe_says_it_is_an_idealization_and_what_it_cost() -> None:
    _, _, parent, child = make_phases()
    gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent, child_phase=child
    )
    text = _measured_report(gt).as_rational_relationship(max_index=2).describe()
    assert "idealization" in text
    assert "2.404" in text
    assert "|index| <= 2" in text
    assert "Compare the residual against the scatter" in text


def test_the_rationalizer_rejects_a_search_it_cannot_run() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    with pytest.raises(ValueError, match="max_index must be at least 1"):
        _measured_report(ks).as_rational_relationship(max_index=0)


def test_exact_kurdjumov_sachs_rationalizes_at_index_one_with_no_tolerance() -> None:
    """Zero tolerance is the strongest thing a caller can ask for, and K-S meets it.

    Both clauses of K-S are index-1 families -- {111} against {110}, <110>
    against <111> -- so an exactly K-S rotation satisfies them to rounding and
    the search has to say so. This used to raise: the deviations were recovered
    with ``arccos`` of a dot product, which cannot report better than about
    ``1e-06`` degrees for a pair that is exactly parallel, and the test that
    stood here asserted the refusal that produced.
    """

    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    result = _measured_report(ks).as_rational_relationship(tolerance_deg=0.0, max_index=1)
    assert result.plane_statement.deviation_deg == pytest.approx(0.0, abs=1e-9)
    assert result.direction_statement.deviation_deg == pytest.approx(0.0, abs=1e-9)


def test_a_plane_without_a_direction_in_it_is_refused_by_name() -> None:
    """The second refusal: a clause found, but nothing to complete it with.

    A rotation about [100] leaves ``(100)`` exactly parallel to itself, so the
    plane clause is there at ``|index| <= 1``; the directions *in* that plane
    are turned by the rotation angle, which is not a low-index parallelism at
    any tolerance worth the name. The error must say which plane it could not
    complete rather than reporting the plane alone as though it were the
    statement.
    """

    from pytex.core import characterize_orientation_relationship, specimen_frame

    _, _, parent, child = make_phases()
    parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
    in_plane_turn = Rotation.from_axis_angle([1.0, 0.0, 0.0], 0.37).as_matrix()
    frame = specimen_frame()
    parents = OrientationSet.from_matrices(
        parent_matrix[None, :, :], specimen_frame=frame, phase=parent
    )
    children = OrientationSet.from_matrices(
        (parent_matrix @ in_plane_turn.T)[None, :, :], specimen_frame=frame, phase=child
    )
    report = characterize_orientation_relationship(parents, children)
    with pytest.raises(ValueError, match="No direction clause lies in"):
        report.as_rational_relationship(tolerance_deg=0.0, max_index=1)


def test_a_rotation_with_no_low_index_plane_is_refused_outright() -> None:
    """The first refusal: nothing to rationalize at all, which asks for a wider
    search rather than for a different completion."""

    from pytex.core import characterize_orientation_relationship, specimen_frame

    _, _, parent, child = make_phases()
    parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
    awkward = Rotation.from_axis_angle([3.0, 5.0, 7.0], 0.37).as_matrix()
    frame = specimen_frame()
    parents = OrientationSet.from_matrices(
        parent_matrix[None, :, :], specimen_frame=frame, phase=parent
    )
    children = OrientationSet.from_matrices(
        (parent_matrix @ awkward.T)[None, :, :], specimen_frame=frame, phase=child
    )
    report = characterize_orientation_relationship(parents, children)
    with pytest.raises(ValueError, match="nothing to rationalize"):
        report.as_rational_relationship(max_index=2, tolerance_deg=0.05)


def test_the_search_is_rerun_rather_than_reusing_the_reports_clause_list() -> None:
    """`max_index` has to mean something. A report carrying no direction clauses
    at all still rationalizes, because the bounds passed here drive a fresh
    search; the report's own clauses only steer which family is preferred."""

    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    report = _measured_report(ks)
    statements = report.plane_statements
    stripped = ORCharacterizationReport(
        relationship=report.relationship,
        pair_count=report.pair_count,
        residuals_deg=report.residuals_deg,
        iterations=report.iterations,
        converged=report.converged,
        catalog_names=report.catalog_names,
        catalog_deviations_deg=report.catalog_deviations_deg,
        best_catalog_name=report.best_catalog_name,
        best_catalog_deviation_deg=report.best_catalog_deviation_deg,
        margin_deg=report.margin_deg,
        catalog_tolerance_deg=report.catalog_tolerance_deg,
        plane_statements=statements,
        direction_statements=(),
        provenance=report.provenance,
    )
    # The search is re-run, so an empty statement list does not by itself starve
    # it; the tolerance is what has to be closed down.
    assert stripped.as_rational_relationship().residual_rotation_deg == pytest.approx(
        0.0, abs=1e-6
    )


def test_the_result_rejects_a_swapped_clause_pair() -> None:
    _, _, parent, child = make_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    result = _measured_report(ks).as_rational_relationship()
    with pytest.raises(ValueError, match="plane_statement must be a plane clause"):
        RationalizedORResult(
            relationship=result.relationship,
            source_relationship_name=result.source_relationship_name,
            plane_statement=result.direction_statement,
            direction_statement=result.direction_statement,
            residual_rotation_deg=0.0,
            zone_law_deviation_deg=0.0,
            max_index=4,
            tolerance_deg=3.0,
        )
