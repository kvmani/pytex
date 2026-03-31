from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Phase,
    PhaseTransformationRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.experimental import score_parent_orientations


def _make_frames_and_phases(
) -> tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame, Phase, Phase]:
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
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    parent_symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal_parent)
    child_symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal_child)
    parent = Phase(
        "parent",
        lattice=Lattice(3.6, 3.8, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal_parent),
        symmetry=parent_symmetry,
        crystal_frame=crystal_parent,
    )
    child = Phase(
        "child",
        lattice=Lattice(3.2, 3.3, 3.5, 90.0, 90.0, 90.0, crystal_frame=crystal_child),
        symmetry=child_symmetry,
        crystal_frame=crystal_child,
    )
    return crystal_parent, crystal_child, specimen, parent, child


def test_experimental_parent_scoring_ranks_true_parent_best() -> None:
    crystal_parent, crystal_child, specimen, parent, child = _make_frames_and_phases()
    relationship = OrientationRelationship(
        name="demo_or",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 6.0),
    )
    variants = relationship.generate_variants()
    true_parent = Orientation(
        rotation=Rotation.from_axis_angle([0.0, 1.0, 0.0], np.pi / 12.0),
        crystal_frame=crystal_parent,
        specimen_frame=specimen,
        symmetry=parent.symmetry,
        phase=parent,
    )
    variant_indices = np.array(
        [variants[0].variant_index, variants[-1].variant_index],
        dtype=np.int64,
    )
    observed_children = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=variants[0].parent_to_child_rotation.compose(true_parent.rotation),
                crystal_frame=crystal_child,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            ),
            Orientation(
                rotation=variants[-1].parent_to_child_rotation.compose(true_parent.rotation),
                crystal_frame=crystal_child,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            ),
        ]
    )
    record = PhaseTransformationRecord(
        name="experimental_record",
        orientation_relationship=relationship,
        parent_orientation=true_parent,
        child_orientations=observed_children,
        variant_indices=variant_indices,
    )
    candidates = OrientationSet.from_orientations(
        [
            true_parent,
            Orientation(
                rotation=Rotation.from_axis_angle([1.0, 0.0, 0.0], np.pi / 9.0),
                crystal_frame=crystal_parent,
                specimen_frame=specimen,
                symmetry=parent.symmetry,
                phase=parent,
            ),
            Orientation(
                rotation=Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 8.0),
                crystal_frame=crystal_parent,
                specimen_frame=specimen,
                symmetry=parent.symmetry,
                phase=parent,
            ),
        ]
    )

    result = score_parent_orientations(record, candidates, symmetry_aware=True, reduction="mean")

    assert result.best_index == 0
    assert result.best_score_rad <= 1e-10
    assert result.scores_rad.shape == (3,)
    predicted = result.predicted_child_orientations()
    assert_allclose(predicted.as_matrices(), observed_children.as_matrices(), atol=1e-8)


def test_experimental_parent_scoring_requires_parent_semantic_match() -> None:
    crystal_parent, _, specimen, parent, child = _make_frames_and_phases()
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
                crystal_frame=child.crystal_frame,
                specimen_frame=specimen,
                symmetry=child.symmetry,
                phase=child,
            )
        ]
    )
    record = PhaseTransformationRecord(
        name="experimental_record",
        orientation_relationship=relationship,
        parent_orientation=parent_orientation,
        child_orientations=child_orientations,
    )
    wrong_phase = Phase(
        "wrong",
        lattice=parent.lattice,
        symmetry=parent.symmetry,
        crystal_frame=crystal_parent,
    )
    candidates = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal_parent,
                specimen_frame=specimen,
                symmetry=parent.symmetry,
                phase=wrong_phase,
            )
        ]
    )

    with pytest.raises(ValueError, match="parent phase semantics"):
        score_parent_orientations(record, candidates)
