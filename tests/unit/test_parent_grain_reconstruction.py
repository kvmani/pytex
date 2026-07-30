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


def test_lath_martensite_structure_full_24_variant_parent() -> None:
    """Literature-structure fixture: one austenite grain, all 24 KS variants.

    Reproduces the lath-martensite hierarchy of Morito et al.: 24 variants
    from one parent form four close-packed packets of six variants; the
    reconstruction must gather all 24 children into one parent, recover the
    parent orientation exactly, and variant selection must recover every
    planted index.
    """

    from pytex.core import (
        CrystalPlane,
        MillerIndex,
        PhaseTransformationRecord,
        select_variants,
        variant_close_packed_groups,
    )

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variants = ks.generate_variants()
    assert len(variants) == 24
    parent_orientation = Orientation.from_euler(
        20.0, 30.0, 40.0, specimen_frame=specimen, symmetry=parent_phase.symmetry,
        phase=parent_phase,
    )
    quaternions = np.stack(
        [
            parent_orientation.rotation.compose(
                variant.parent_to_child_rotation.inverse()
            ).quaternion
            for variant in variants
        ],
        axis=0,
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=child_phase.crystal_frame,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    # Packet structure: four close-packed packets of six variants each.
    packets = variant_close_packed_groups(
        ks,
        CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent_phase), phase=parent_phase),
    )
    assert packets.shape == (24,)
    assert np.bincount(packets).tolist() == [6, 6, 6, 6]
    # Reconstruction: chain adjacency over all 24 children -> one parent.
    adjacency = np.asarray(list(pairwise(range(24))), dtype=np.int64)
    result = reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=2.0)
    assert result.parent_count == 1
    assert result.cluster_sizes.tolist() == [24]
    assert result.edges_linked == adjacency.shape[0]
    assert result.max_deviation_deg == pytest.approx(0.0, abs=1e-8)
    distance = _symmetry_reduced_angle_between_deg(
        result.parent_orientations.as_matrices()[0],
        parent_orientation.rotation.as_matrix(),
        child_operators=np.eye(3, dtype=np.float64)[None, :, :],
        parent_operators=parent_phase.symmetry.operators,
    )
    assert distance == pytest.approx(0.0, abs=1e-6)
    # Variant selection recovers every planted variant index.
    record = PhaseTransformationRecord(
        name="lath_block",
        orientation_relationship=ks,
        parent_orientation=parent_orientation,
        child_orientations=children,
    )
    report = select_variants(record)
    assert report.variant_indices.tolist() == [variant.variant_index for variant in variants]
    assert report.scores_deg.max() == pytest.approx(0.0, abs=1e-8)


def test_burgers_variant_groups_form_six_pairs() -> None:
    from pytex.core import CrystalPlane, MillerIndex, variant_close_packed_groups

    parent_phase, _, _ = _phases()
    hex_frame = ReferenceFrame(
        name="recon_hex_child",
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
    burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent_phase, child_phase=alpha
    )
    groups = variant_close_packed_groups(
        burgers,
        CrystalPlane(MillerIndex(np.array([1, 1, 0]), phase=parent_phase), phase=parent_phase),
    )
    # Twelve Burgers variants inherit their basal plane from the six {110}
    # parent planes: six groups of two.
    assert groups.shape == (12,)
    assert np.bincount(groups).tolist() == [2, 2, 2, 2, 2, 2]


def test_variant_pole_figure_predicts_packet_plane_coincidence() -> None:
    """Each variant's child {011} pole set contains its packet's parent {111} pole.

    The defining KS parallelism ties every variant's (011) child plane to its
    close-packed parent {111} member, so the predicted specimen-frame poles
    must coincide — the textbook signature read off variant pole figures.
    """

    from pytex.core import CrystalPlane, MillerIndex, variant_pole_figure
    from pytex.core.transformation import _integer_index_orbit

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    parent_orientation = Orientation.from_euler(
        20.0, 30.0, 40.0, specimen_frame=specimen, symmetry=parent_phase.symmetry,
        phase=parent_phase,
    )
    child_plane = CrystalPlane(
        MillerIndex(np.array([0, 1, 1]), phase=child_phase), phase=child_phase
    )
    prediction = variant_pole_figure(parent_orientation, ks, child_plane)
    assert prediction.poles.values.shape == (24 * 6, 3)
    assert np.allclose(np.linalg.norm(prediction.poles.values, axis=1), 1.0)
    assert prediction.poles.reference_frame == specimen
    assert prediction.variant_count == 24
    members = _integer_index_orbit(
        np.array([1, 1, 1]), phase=parent_phase, reciprocal=True
    )
    reciprocal = parent_phase.lattice.reciprocal_basis().matrix
    parent_normals = members.astype(np.float64) @ reciprocal.T
    parent_normals /= np.linalg.norm(parent_normals, axis=1)[:, None]
    specimen_normals = parent_normals @ parent_orientation.rotation.as_matrix().T
    for variant in ks.generate_variants():
        rows = prediction.poles.values[
            np.flatnonzero(prediction.variant_indices == variant.variant_index)
        ]
        residuals = [
            ks.map_plane_to_child(
                CrystalPlane(MillerIndex(member, phase=parent_phase), phase=parent_phase),
                variant=variant,
            ).angular_residual_deg
            for member in members
        ]
        target = specimen_normals[int(np.argmin(residuals))]
        assert float(np.max(np.abs(rows @ target))) == pytest.approx(1.0, abs=1e-9)
    text = prediction.describe()
    assert "24 variant(s)" in text
    # A pole figure plots the whole symmetry orbit, so the family brackets
    # {hkl} are the correct notation for the quantity, not (hkl).
    assert "{011}" in text


def test_variant_pole_figure_validates_phases_and_plots() -> None:
    from pytex.core import CrystalPlane, MillerIndex, variant_pole_figure
    from pytex.plotting import plot_variant_pole_figure

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    parent_orientation = Orientation.from_euler(
        10.0, 20.0, 30.0, specimen_frame=specimen, symmetry=parent_phase.symmetry,
        phase=parent_phase,
    )
    with pytest.raises(ValueError, match="child phase"):
        variant_pole_figure(
            parent_orientation,
            ks,
            CrystalPlane(MillerIndex(np.array([0, 1, 1]), phase=parent_phase),
                         phase=parent_phase),
        )
    prediction = variant_pole_figure(
        parent_orientation,
        ks,
        CrystalPlane(MillerIndex(np.array([0, 1, 1]), phase=child_phase), phase=child_phase),
    )
    axes = plot_variant_pole_figure(prediction)
    assert axes is not None


def _many_parent_microstructure(
    *, seed: int = 187, parent_count: int = 6, per_parent: int = 5
) -> tuple[OrientationSet, np.ndarray, np.ndarray, OrientationRelationship]:
    """Randomly oriented parents, each transformed through distinct KS variants.

    Adjacency chains the children of each parent and adds one contact edge
    between consecutive parents, so the fixture contains genuine cross-parent
    boundaries of exactly the kind the edge test must reject. Seed 187 is
    pinned because its cross-parent boundaries are all far from the same-parent
    fingerprint, which the test asserts rather than assumes.
    """

    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variants = ks.generate_variants()
    rng = np.random.default_rng(seed)
    quaternions: list[np.ndarray] = []
    planted: list[int] = []
    for parent_index in range(parent_count):
        quaternion = rng.normal(size=4)
        quaternion /= np.linalg.norm(quaternion)
        parent_rotation = Rotation(quaternion=quaternion)
        for pick in rng.choice(len(variants), size=per_parent, replace=False):
            quaternions.append(
                parent_rotation.compose(
                    variants[int(pick)].parent_to_child_rotation.inverse()
                ).quaternion
            )
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
    for parent_index in range(parent_count):
        members = np.flatnonzero(planted_labels == parent_index)
        edges.extend((int(a), int(b)) for a, b in pairwise(members))
    for parent_index in range(parent_count - 1):
        edges.append(
            (
                int(np.flatnonzero(planted_labels == parent_index)[-1]),
                int(np.flatnonzero(planted_labels == parent_index + 1)[0]),
            )
        )
    return children, np.asarray(edges, dtype=np.int64), planted_labels, ks


def test_edge_test_matches_the_full_rotation_not_only_the_angle() -> None:
    """Regression: an angle-only edge test merges unrelated parent grains.

    The same-parent fingerprint is the double coset ``G_c (R G_p R^T) G_c``. A
    test that compares only the misorientation *angle* against the intervariant
    spectrum discards the axis, and for a cubic-cubic relationship that spectrum
    is dense enough that a 3 deg window admits most unrelated boundaries.

    This fixture makes the consequence concrete. Its cross-parent boundaries are
    first asserted to be unambiguously far from the fingerprint, so the planted
    partition is the only defensible answer; the angle-only rule nonetheless
    links most of them and collapses distinct parents together, while the
    shipped full-rotation rule recovers the planted partition exactly.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
        intervariant_misorientation_angles_deg,
    )
    from pytex.core.orientation import _reduced_pair_disorientation_angles

    children, adjacency, planted, ks = _many_parent_microstructure()
    matrices = children.as_matrices()
    relative = np.einsum(
        "eji,ejk->eik",
        matrices[adjacency[:, 0]],
        matrices[adjacency[:, 1]],
        optimize=True,
    )
    cross = planted[adjacency[:, 0]] != planted[adjacency[:, 1]]
    assert int(np.count_nonzero(cross)) == 5

    fingerprint = intervariant_boundary_fingerprint(ks)
    distances = boundary_fingerprint_distances_deg(relative, fingerprint)
    # Ground truth is unambiguous: same-parent edges sit on the fingerprint and
    # cross-parent edges sit far off it, well beyond any sane tolerance.
    assert float(distances[~cross].max()) < 1e-5
    assert float(distances[cross].min()) > 8.0

    # The discarded-axis rule links most cross-parent boundaries anyway.
    operators = ks.child_phase.symmetry.operators
    table = intervariant_misorientation_angles_deg(ks)
    spectrum = np.unique(
        np.round(np.concatenate([[0.0], table[np.triu_indices(table.shape[0], k=1)]]), 6)
    )
    angles = np.degrees(_reduced_pair_disorientation_angles(relative, operators, operators))
    angle_only_linked = (
        np.min(np.abs(angles[:, None] - spectrum[None, :]), axis=1) <= 3.0
    )
    assert int(np.count_nonzero(angle_only_linked & cross)) >= 3

    # The shipped rule rejects every one of them and recovers the partition.
    result = reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=3.0)
    assert result.edges_linked == int(np.count_nonzero(~cross))
    assert result.parent_count == 6
    assert _partitions_equal(result.parent_labels, planted)
    assert result.max_deviation_deg == pytest.approx(0.0, abs=1e-6)


def test_reconstruction_scales_to_many_parents_without_merging_them() -> None:
    """Partition recovery must hold across independent random microstructures.

    A single pinned fixture can flatter an edge test, so this sweeps several
    seeds and requires that the reconstruction never merges two planted parents
    whose boundary is genuinely separable, and never splits a parent whose
    children are all linked.
    """

    from pytex.core import (
        boundary_fingerprint_distances_deg,
        intervariant_boundary_fingerprint,
    )

    for seed in (11, 23, 57, 91):
        children, adjacency, planted, ks = _many_parent_microstructure(seed=seed)
        matrices = children.as_matrices()
        relative = np.einsum(
            "eji,ejk->eik",
            matrices[adjacency[:, 0]],
            matrices[adjacency[:, 1]],
            optimize=True,
        )
        cross = planted[adjacency[:, 0]] != planted[adjacency[:, 1]]
        distances = boundary_fingerprint_distances_deg(
            relative, intervariant_boundary_fingerprint(ks)
        )
        # Only judge seeds whose ground truth is unambiguous at this tolerance;
        # a cross-parent boundary that genuinely lands on the fingerprint is a
        # real physical ambiguity, not an algorithmic defect.
        if float(distances[cross].min()) <= 3.0:
            continue
        result = reconstruct_parent_grains(children, adjacency, ks, tolerance_deg=3.0)
        assert _partitions_equal(result.parent_labels, planted), f"seed {seed}"
        assert result.parent_count == 6, f"seed {seed}"
        assert result.max_deviation_deg == pytest.approx(0.0, abs=1e-6), f"seed {seed}"
