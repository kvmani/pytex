from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from pytex.core import (
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Rotation,
    standard_fcc_bcc_relationships,
)
from pytex.experimental import identify_orientation_relationship
from tests.unit.test_parent_grain_reconstruction import _phases


def _microstructure(
    builder, *, noise_deg: float = 0.0, seed: int = 3
) -> tuple[OrientationSet, np.ndarray]:
    """Three parents, five children each from `builder`'s OR, same-parent edges."""

    parent_phase, child_phase, specimen = _phases()
    relationship = builder(parent_phase=parent_phase, child_phase=child_phase)
    variants = relationship.generate_variants()
    rng = np.random.default_rng(seed)
    quaternions: list[np.ndarray] = []
    planted: list[int] = []
    parent_eulers = [(20.0, 30.0, 40.0), (10.0, 50.0, 20.0), (65.0, 20.0, 50.0)]
    for parent_index, euler in enumerate(parent_eulers):
        parent = Orientation.from_euler(
            *euler, specimen_frame=specimen, symmetry=parent_phase.symmetry, phase=parent_phase
        )
        for pick in rng.choice(len(variants), size=5, replace=False):
            rotation = parent.rotation.compose(
                variants[int(pick)].parent_to_child_rotation.inverse()
            )
            if noise_deg > 0.0:
                axis = rng.normal(size=3)
                axis /= np.linalg.norm(axis)
                rotation = Rotation.from_axis_angle(
                    axis, np.deg2rad(rng.normal(0.0, noise_deg))
                ).compose(rotation)
            quaternions.append(rotation.quaternion)
            planted.append(parent_index)
    children = OrientationSet(
        quaternions=np.stack(quaternions, axis=0),
        crystal_frame=child_phase.crystal_frame,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    labels = np.asarray(planted)
    edges = []
    for parent_index in range(3):
        members = np.flatnonzero(labels == parent_index)
        edges.extend((int(a), int(b)) for a, b in pairwise(members))
    return children, np.asarray(edges, dtype=np.int64)


def _catalog():
    parent_phase, child_phase, _ = _phases()
    return standard_fcc_bcc_relationships(
        parent_phase=parent_phase, child_phase=child_phase
    ).relationships


@pytest.mark.parametrize(
    "builder_name",
    ["from_kurdjumov_sachs_correspondence", "from_greninger_troiano_correspondence"],
)
def test_identification_recovers_the_generating_relationship(builder_name: str) -> None:
    builder = getattr(OrientationRelationship, builder_name)
    children, edges = _microstructure(builder)
    report = identify_orientation_relationship(children, edges, _catalog())
    generating = builder(
        parent_phase=_phases()[0], child_phase=_phases()[1]
    ).name
    assert report.best_name == generating
    assert float(np.min(report.mean_distances_deg)) == pytest.approx(0.0, abs=1e-8)
    # Every other candidate sits at a clearly positive distance.
    others = [
        distance
        for name, distance in zip(report.candidate_names, report.mean_distances_deg, strict=False)
        if name != generating
    ]
    assert min(others) > 1.0


def test_identification_survives_orientation_noise() -> None:
    children, edges = _microstructure(
        OrientationRelationship.from_kurdjumov_sachs_correspondence, noise_deg=0.3
    )
    report = identify_orientation_relationship(children, edges, _catalog())
    assert report.best_name == "kurdjumov_sachs"
    assert float(np.min(report.mean_distances_deg)) < 1.0
    text = report.describe()
    assert "no parent orientations used" in text
    assert "kurdjumov_sachs" in text
    assert "refinement of the rotation itself is not implemented" in text


def test_identification_validates_inputs() -> None:
    children, edges = _microstructure(
        OrientationRelationship.from_kurdjumov_sachs_correspondence
    )
    with pytest.raises(ValueError, match="requires candidates"):
        identify_orientation_relationship(children, edges, ())
    with pytest.raises(ValueError, match="self-edges"):
        identify_orientation_relationship(children, np.array([[1, 1]]), _catalog())
    with pytest.raises(ValueError, match="at least one edge"):
        identify_orientation_relationship(
            children, np.empty((0, 2), dtype=np.int64), _catalog()
        )
