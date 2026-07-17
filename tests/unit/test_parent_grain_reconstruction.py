from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from pytex.core import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.core.transformation import _symmetry_reduced_angle_between_deg
from pytex.experimental import reconstruct_parent_grains


def _phases() -> tuple[Phase, Phase, ReferenceFrame]:
    parent_frame = ReferenceFrame(
        name="recon_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    child_frame = ReferenceFrame(
        name="recon_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="recon_specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    parent = Phase(
        "austenite",
        lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
        crystal_frame=parent_frame,
    )
    child = Phase(
        "ferrite",
        lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
        crystal_frame=child_frame,
    )
    return parent, child, specimen


def _planted_microstructure(
    *, noise_deg: float = 0.0, seed: int = 3
) -> tuple[OrientationSet, np.ndarray, np.ndarray, list[Orientation], OrientationRelationship]:
    """Three planted parents, five KS children each, chain adjacency plus cross edges."""

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variants = ks.generate_variants()
    parent_orientations = [
        Orientation.from_euler(
            *euler, specimen_frame=specimen, symmetry=parent_phase.symmetry, phase=parent_phase
        )
        for euler in [(0.0, 0.0, 0.0), (10.0, 50.0, 20.0), (75.0, 33.0, 44.0)]
    ]
    rng = np.random.default_rng(seed)
    quaternions: list[np.ndarray] = []
    planted: list[int] = []
    for parent_index, parent_orientation in enumerate(parent_orientations):
        picks = rng.choice(len(variants), size=5, replace=False)
        for pick in picks:
            rotation = parent_orientation.rotation.compose(
                variants[int(pick)].parent_to_child_rotation.inverse()
            )
            if noise_deg > 0.0:
                axis = rng.normal(size=3)
                axis /= np.linalg.norm(axis)
                perturbation = Rotation.from_axis_angle(
                    axis, np.deg2rad(rng.normal(0.0, noise_deg))
                )
                rotation = perturbation.compose(rotation)
            quaternions.append(rotation.quaternion)
            planted.append(parent_index)
    children = OrientationSet(
        quaternions=np.stack(quaternions, axis=0),
        crystal_frame=child_phase.crystal_frame,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    planted_labels = np.asarray(planted, dtype=np.int64)
    edges: list[tuple[int, int]] = []
    for parent_index in range(3):
        members = np.flatnonzero(planted_labels == parent_index)
        edges.extend(
            (int(a), int(b)) for a, b in pairwise(members)
        )
    # Two cross-parent edges that must be rejected by the fingerprint test.
    adjacency = np.asarray([*edges, (0, 5), (7, 12)], dtype=np.int64)
    return children, adjacency, planted_labels, parent_orientations, ks


def _partitions_equal(left: np.ndarray, right: np.ndarray) -> bool:
    count = left.size
    return all(
        (left[i] == left[j]) == (right[i] == right[j])
        for i in range(count)
        for j in range(count)
    )


def test_reconstruction_recovers_planted_partition_and_parents_exactly() -> None:
    children, adjacency, planted, parents, ks = _planted_microstructure()
    result = reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=2.0)
    assert result.parent_count == 3
    assert _partitions_equal(result.parent_labels, planted)
    # The two cross-parent edges must have been rejected.
    assert result.edges_tested == adjacency.shape[0]
    assert result.edges_linked == adjacency.shape[0] - 2
    assert result.max_deviation_deg == pytest.approx(0.0, abs=1e-8)
    parent_symmetry = ks.parent_phase.symmetry
    for parent_index, parent_orientation in enumerate(parents):
        cluster = int(result.parent_labels[np.flatnonzero(planted == parent_index)[0]])
        # Canonical crystal->specimen parents are equivalent up to RIGHT
        # multiplication by parent crystal operators: P' = P S_p.
        distance = _symmetry_reduced_angle_between_deg(
            result.parent_orientations.as_matrices()[cluster],
            parent_orientation.rotation.as_matrix(),
            child_operators=np.eye(3, dtype=np.float64)[None, :, :],
            parent_operators=parent_symmetry.operators,
        )
        assert distance == pytest.approx(0.0, abs=1e-6)


def test_reconstruction_tolerates_orientation_noise() -> None:
    children, adjacency, planted, parents, ks = _planted_microstructure(noise_deg=0.3)
    result = reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=3.0)
    assert result.parent_count == 3
    assert _partitions_equal(result.parent_labels, planted)
    # Residuals sit at the noise level, and the quaternion-averaged parent
    # estimate must beat any single member's noise (sigma/sqrt(n) behavior).
    assert result.mean_deviation_deg < 1.0
    parent_symmetry = ks.parent_phase.symmetry
    for parent_index, parent_orientation in enumerate(parents):
        cluster = int(result.parent_labels[np.flatnonzero(planted == parent_index)[0]])
        # Canonical crystal->specimen parents are equivalent up to RIGHT
        # multiplication by parent crystal operators: P' = P S_p.
        distance = _symmetry_reduced_angle_between_deg(
            result.parent_orientations.as_matrices()[cluster],
            parent_orientation.rotation.as_matrix(),
            child_operators=np.eye(3, dtype=np.float64)[None, :, :],
            parent_operators=parent_symmetry.operators,
        )
        assert distance < 0.5


def test_singleton_grains_are_reported_as_ambiguous() -> None:
    children, adjacency, planted, _, ks = _planted_microstructure()
    # Drop every edge touching grain 0 so it becomes a singleton cluster.
    kept = adjacency[(adjacency[:, 0] != 0) & (adjacency[:, 1] != 0)]
    result = reconstruct_parent_grains(children, kept, ks, tolerance_deg=2.0)
    assert result.singleton_count >= 1
    text = result.describe()
    assert "singleton" in text
    assert "kurdjumov_sachs" in text


def test_reconstruction_validates_inputs() -> None:
    children, adjacency, _, _, ks = _planted_microstructure()
    with pytest.raises(ValueError, match="child phase"):
        parent_phase, _, specimen = _phases()
        wrong = OrientationSet(
            quaternions=children.quaternions,
            crystal_frame=parent_phase.crystal_frame,
            specimen_frame=specimen,
            symmetry=parent_phase.symmetry,
            phase=parent_phase,
        )
        reconstruct_parent_grains(wrong, adjacency, ks)
    with pytest.raises(ValueError, match="self-edges"):
        reconstruct_parent_grains(children, np.array([[1, 1]]), ks)
    with pytest.raises(ValueError, match="reference child_orientations"):
        reconstruct_parent_grains(children, np.array([[0, 99]]), ks)
    with pytest.raises(ValueError, match="tolerance_deg"):
        reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=0.0)


def test_reconstruction_from_ebsd_grain_graph() -> None:
    from pytex.ebsd import CrystalMap
    from pytex.experimental import reconstruct_parent_grains_from_graph

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variants = ks.generate_variants()
    parent_orientations = [
        Orientation.from_euler(
            *euler, specimen_frame=specimen, symmetry=parent_phase.symmetry, phase=parent_phase
        )
        # General orientations: a cube-oriented parent makes some variant
        # pairs symmetry-degenerate as child orientations, and the second
        # parent is chosen so the cross-parent boundary disorientation sits
        # ~5 deg from every intervariant fingerprint angle.
        for euler in [(20.0, 30.0, 40.0), (65.0, 20.0, 50.0)]
    ]
    # 2x4 pixel map: columns 0-1 are two variants of parent 0, columns 2-3 two
    # variants of parent 1; each column is one grain (two identical pixels).
    column_rotations = [
        parent_orientations[0].rotation.compose(variants[0].parent_to_child_rotation.inverse()),
        parent_orientations[0].rotation.compose(variants[7].parent_to_child_rotation.inverse()),
        parent_orientations[1].rotation.compose(variants[2].parent_to_child_rotation.inverse()),
        parent_orientations[1].rotation.compose(variants[11].parent_to_child_rotation.inverse()),
    ]
    pixel_orientations = [
        Orientation(
            rotation=column_rotations[x],
            crystal_frame=child_phase.crystal_frame,
            specimen_frame=specimen,
            symmetry=child_phase.symmetry,
            phase=child_phase,
        )
        for _y in range(2)
        for x in range(4)
    ]
    coordinates = np.array(
        [[float(x), float(y)] for y in range(2) for x in range(4)], dtype=np.float64
    )
    crystal_map = CrystalMap(
        coordinates=coordinates,
        orientations=OrientationSet.from_orientations(pixel_orientations),
        map_frame=specimen,
        grid_shape=(2, 4),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=2.0, symmetry_aware=True, connectivity=4
    )
    graph = segmentation.grain_graph()
    assert len(graph.node_grain_ids) == 4
    result = reconstruct_parent_grains_from_graph(graph, ks, tolerance_deg=2.0)
    assert result.grain_ids is not None
    assert result.grain_ids.shape == result.parent_labels.shape
    assert result.parent_count == 2
    assert result.cluster_sizes.tolist() == [2, 2]
    assert result.max_deviation_deg == pytest.approx(0.0, abs=1e-6)
