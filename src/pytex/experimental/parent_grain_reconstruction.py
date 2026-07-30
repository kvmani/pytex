"""Experimental map-scale parent-grain reconstruction from variant graphs.

Given child-grain mean orientations, a grain-adjacency list, and a nominal
orientation relationship, this module groups child grains that are mutually
consistent with a common parent (their boundary disorientations match the
relationship's intervariant fingerprint), then estimates one parent
orientation per group. This is the bounded v1 of the parent-grain
reconstruction flagship (development-guide finding 4 / OR-foundation F8):
grain-mean level, explicit adjacency input, staged under
``pytex.experimental`` until literature-fixture validation broadens.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core.lattice import phases_semantically_match
from pytex.core.orientation import (
    OrientationSet,
    Rotation,
    _reduced_pair_disorientation_angles,
    matrices_to_quaternions,
)
from pytex.core.point_groups import normalize_point_group_symbol
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import (
    OrientationRelationship,
    boundary_fingerprint_distances_deg,
    intervariant_boundary_fingerprint,
)
from pytex.ebsd.models import GrainGraph


@dataclass(frozen=True, slots=True)
class ParentGrainReconstructionResult:
    """Grouping of child grains into parents with estimated parent orientations.

    ``parent_labels[i]`` is the 0-based parent-cluster index of child grain
    ``i``; ``parent_orientations[k]`` is the estimated orientation of cluster
    ``k``; ``per_grain_deviation_deg[i]`` is the child-symmetry-reduced angle
    between observed child ``i`` and its best variant prediction from the
    assigned parent. Singleton clusters (no orientation-consistent neighbor)
    are 24-fold ambiguous for cubic-cubic relationships: their parent is the
    exact-fit candidate from the only child and must be treated as one of
    several equivalent possibilities.

    ``chance_link_probability`` is the probability that two *unrelated* grains
    would be linked at this tolerance — the fraction of uniformly random
    misorientations lying within ``tolerance_deg`` of the same-parent
    fingerprint. It is the controlling reliability number for the clustering:
    a densely connected grain graph gives many chances for such a coincidence,
    and one is enough to merge two parents irreversibly. See
    ``describe()``, which warns when the expected number of chance links is
    not small.
    """

    relationship_name: str
    parent_labels: np.ndarray
    parent_orientations: OrientationSet
    per_grain_deviation_deg: np.ndarray
    cluster_sizes: np.ndarray
    edges_tested: int
    edges_linked: int
    tolerance_deg: float
    chance_link_probability: float = 0.0
    grain_ids: np.ndarray | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        labels = np.asarray(self.parent_labels, dtype=np.int64).reshape(-1)
        deviations = np.asarray(self.per_grain_deviation_deg, dtype=np.float64).reshape(-1)
        sizes = np.asarray(self.cluster_sizes, dtype=np.int64).reshape(-1)
        if labels.shape != deviations.shape:
            raise ValueError("parent_labels and per_grain_deviation_deg must align.")
        if labels.size == 0:
            raise ValueError("ParentGrainReconstructionResult requires at least one grain.")
        if sizes.shape != (len(self.parent_orientations),):
            raise ValueError("cluster_sizes must have one entry per parent orientation.")
        if np.any(labels < 0) or np.any(labels >= sizes.size):
            raise ValueError("parent_labels must index into the parent clusters.")
        if np.any(~np.isfinite(deviations)) or np.any(deviations < 0.0):
            raise ValueError("per_grain_deviation_deg must be finite and non-negative.")
        for array in (labels, deviations, sizes):
            array.setflags(write=False)
        object.__setattr__(self, "parent_labels", labels)
        object.__setattr__(self, "per_grain_deviation_deg", deviations)
        object.__setattr__(self, "cluster_sizes", sizes)
        if self.grain_ids is not None:
            grain_ids = np.asarray(self.grain_ids, dtype=np.int64).reshape(-1)
            if grain_ids.shape != labels.shape:
                raise ValueError("grain_ids must align with parent_labels.")
            grain_ids.setflags(write=False)
            object.__setattr__(self, "grain_ids", grain_ids)

    @property
    def parent_count(self) -> int:
        return int(self.cluster_sizes.size)

    @property
    def singleton_count(self) -> int:
        return int(np.count_nonzero(self.cluster_sizes == 1))

    @property
    def mean_deviation_deg(self) -> float:
        return float(np.mean(self.per_grain_deviation_deg))

    @property
    def max_deviation_deg(self) -> float:
        return float(np.max(self.per_grain_deviation_deg))

    def describe(self) -> str:
        """Prose summary: clusters found, fit quality, ambiguity caveats."""

        singleton_note = (
            f" {self.singleton_count} cluster(s) are singletons whose parent is "
            "symmetry-ambiguous (no orientation-consistent neighbor)."
            if self.singleton_count
            else ""
        )
        # The controlling reliability number: how many of the tested edges would
        # be expected to link by chance alone. One chance link merges two
        # parents irreversibly, so the warning triggers on the expected count,
        # not the rate — a small rate over a dense graph is still unreliable.
        expected_chance_links = self.chance_link_probability * self.edges_tested
        if expected_chance_links >= 1.0:
            reliability_note = (
                f" WARNING: at this tolerance {self.chance_link_probability:.1%} of "
                f"unrelated grain pairs fall within the same-parent fingerprint, so "
                f"about {expected_chance_links:.0f} of the {self.edges_tested} tested "
                "edge(s) are expected to link by chance alone. One such link merges "
                "two parents irreversibly, so these groupings are unreliable — reduce "
                "tolerance_deg, or treat merged clusters as unresolved."
            )
        else:
            reliability_note = (
                f" At this tolerance {self.chance_link_probability:.1%} of unrelated "
                f"grain pairs would fall within the fingerprint, so fewer than one of "
                f"the {self.edges_tested} tested edge(s) is expected to link by chance."
            )
        return (
            f"Parent-grain reconstruction under orientation relationship "
            f"'{self.relationship_name}': {self.parent_labels.size} child grain(s) grouped "
            f"into {self.parent_count} parent(s) using {self.edges_linked} of "
            f"{self.edges_tested} adjacency edge(s) whose boundary misorientations matched "
            f"the same-parent fingerprint (full rotation, axis and angle) within "
            f"{self.tolerance_deg:.2f} deg. Residual of "
            f"each child to its assigned parent (best variant, child-symmetry-reduced): mean "
            f"{self.mean_deviation_deg:.3f} deg, max {self.max_deviation_deg:.3f} deg."
            f"{singleton_note}{reliability_note} Experimental surface: measured-data "
            "validation is still ahead; treat groupings as hypotheses at map scale."
        )


#: Sample count for the chance-link estimate. Fixed and seeded so the reported
#: probability is deterministic; 4096 resolves a ~0.3% rate to roughly a third
#: of its own value, which is enough to drive the reliability warning.
_CHANCE_SAMPLES = 4096
_CHANCE_SEED = 20260730


def _chance_link_probability(fingerprint: np.ndarray, tolerance_deg: float) -> float:
    """Probability that two unrelated grains fall within tolerance of the set.

    Estimated by scoring uniformly random misorientations against the
    fingerprint. This is a property of the relationship and the tolerance
    alone — not of the data — so it can be reported before any conclusion is
    drawn from the clustering.
    """

    generator = np.random.default_rng(_CHANCE_SEED)
    quaternions = generator.normal(size=(_CHANCE_SAMPLES, 4))
    quaternions /= np.linalg.norm(quaternions, axis=1)[:, None]
    w, x, y, z = quaternions.T
    matrices = np.empty((_CHANCE_SAMPLES, 3, 3), dtype=np.float64)
    matrices[:, 0, 0] = 1.0 - 2.0 * (y * y + z * z)
    matrices[:, 0, 1] = 2.0 * (x * y - z * w)
    matrices[:, 0, 2] = 2.0 * (x * z + y * w)
    matrices[:, 1, 0] = 2.0 * (x * y + z * w)
    matrices[:, 1, 1] = 1.0 - 2.0 * (x * x + z * z)
    matrices[:, 1, 2] = 2.0 * (y * z - x * w)
    matrices[:, 2, 0] = 2.0 * (x * z - y * w)
    matrices[:, 2, 1] = 2.0 * (y * z + x * w)
    matrices[:, 2, 2] = 1.0 - 2.0 * (x * x + y * y)
    distances = boundary_fingerprint_distances_deg(matrices, fingerprint)
    return float(np.count_nonzero(distances <= tolerance_deg) / _CHANCE_SAMPLES)


def _variant_descriptions(
    variant_matrices: np.ndarray, child_operators: np.ndarray
) -> np.ndarray:
    r"""Every description of a parent-to-child variant under child symmetry.

    A child of parent :math:`\mathbf{P}` through variant :math:`\mathbf{V}_k` is
    :math:`\mathbf{C} = \mathbf{P}\mathbf{V}_k^{\mathsf{T}}`, so
    :math:`\mathbf{C}^{\mathsf{T}}\mathbf{P} = \mathbf{V}_k`. The measured child
    orientation is only defined up to its own crystal symmetry, so the testable
    statement "grain ``j`` descends from parent ``P``" is that
    :math:`\mathbf{C}_j^{\mathsf{T}}\mathbf{P}` lies near the set
    :math:`\{S_c \mathbf{V}_k\}`, returned here deduplicated.
    """

    products = np.einsum(
        "aij,vjk->avik", child_operators, variant_matrices, optimize=True
    ).reshape(-1, 3, 3)
    quaternions = matrices_to_quaternions(products)
    pivot = np.argmax(np.abs(quaternions), axis=1)
    signs = np.sign(quaternions[np.arange(quaternions.shape[0]), pivot])
    keys = np.round(quaternions * signs[:, None], 8)
    _, unique_indices = np.unique(keys, axis=0, return_index=True)
    descriptions = np.ascontiguousarray(products[unique_indices])
    descriptions.setflags(write=False)
    return descriptions


def _vote_partition_cluster(
    cluster_children: np.ndarray,
    *,
    descriptions: np.ndarray,
    variant_matrices: np.ndarray,
    tolerance_deg: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split one connected cluster into groups a single parent can explain.

    Connectivity alone is not evidence of a shared parent. Two unrelated parent
    grains that happen to share one fingerprint-consistent boundary are joined
    by union-find into a cluster no single parent orientation explains, and on
    a dense grain graph that coincidence is common enough to dominate the error
    (see the reconstruction robustness study).

    The split is decided by agreement rather than by connectivity: every member
    proposes the parent it implies, each proposal is scored by how many members
    of the cluster it explains — ``C_j^T P`` near the variant-description set —
    and the best-supported proposal claims its supporters. Unexplained members
    repeat the vote, so a cluster spanning two parents separates into two
    groups while a genuine single-parent cluster is returned untouched.

    Only one proposal per member is needed: the candidates ``C_i V_k`` for
    different ``k`` differ by a parent symmetry operation, and the support set
    is invariant under that, so they are the same hypothesis.

    Returns ``(local_member_indices, seed_parent_matrix)`` per group.
    """

    member_count = cluster_children.shape[0]
    remaining = np.arange(member_count, dtype=np.int64)
    groups: list[tuple[np.ndarray, np.ndarray]] = []
    while remaining.size:
        children = cluster_children[remaining]
        proposals = np.einsum("nij,jk->nik", children, variant_matrices[0], optimize=True)
        # relative[c, j] = C_j^T P_c, for every proposal against every member.
        relative = np.einsum("nji,pjk->pnik", children, proposals, optimize=True)
        distances = boundary_fingerprint_distances_deg(
            relative.reshape(-1, 3, 3), descriptions
        ).reshape(remaining.size, remaining.size)
        supported = distances <= tolerance_deg
        # A proposal always explains its own proposer, so support is never
        # empty and the loop always makes progress.
        best = int(np.argmax(supported.sum(axis=1)))
        claimed = supported[best]
        groups.append((remaining[claimed], proposals[best]))
        remaining = remaining[~claimed]
    return groups


def _child_operators(relationship: OrientationRelationship) -> np.ndarray:
    symmetry = relationship.child_phase.symmetry
    if symmetry is None:
        return np.eye(3, dtype=np.float64)[None, :, :]
    return np.asarray(symmetry.operators, dtype=np.float64)


def _member_deviations_deg(
    cluster_children: np.ndarray,
    parent_matrix: np.ndarray,
    *,
    variant_matrices: np.ndarray,
    operators: np.ndarray,
) -> np.ndarray:
    """Min-over-variants disorientation of each member to one parent estimate."""

    predicted = np.einsum("ij,lkj->lik", parent_matrix, variant_matrices, optimize=True)
    relative = np.einsum("nji,ljk->nlik", cluster_children, predicted, optimize=True)
    member_count, variant_count = relative.shape[:2]
    angles = np.degrees(
        _reduced_pair_disorientation_angles(relative.reshape(-1, 3, 3), operators, operators)
    ).reshape(member_count, variant_count)
    return np.asarray(np.min(angles, axis=1), dtype=np.float64)


def _refine_cluster_parent(
    cluster_children: np.ndarray,
    seed_parent: np.ndarray,
    *,
    variant_matrices: np.ndarray,
    parent_operators: np.ndarray,
) -> np.ndarray:
    """Average every member's parent estimate around the seed (noise reduction).

    Each member contributes the candidate descriptions ``C_i V_k S_p`` over
    all variants and parent symmetry operators; the description with the
    largest trace against the seed is that member's aligned parent estimate,
    and the quaternion eigen-mean of the aligned estimates is the refined
    parent. On exact data this reproduces the seed; on noisy data it averages
    the per-member noise instead of inheriting the first member's.
    """

    estimates = np.einsum(
        "nij,vjk->nvik", cluster_children, variant_matrices, optimize=True
    )
    # Crystal-symmetry equivalents of a crystal->specimen orientation act by
    # right multiplication: P' = P S_p.
    described = np.einsum(
        "nvij,ajk->navik", estimates, parent_operators, optimize=True
    )
    member_count = cluster_children.shape[0]
    flat = described.reshape(member_count, -1, 3, 3)
    scores = np.einsum("ncij,ij->nc", flat, seed_parent, optimize=True)
    aligned = flat[np.arange(member_count), np.argmax(scores, axis=1)]
    quaternions = matrices_to_quaternions(aligned)
    scatter = quaternions.T @ quaternions
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    mean_quaternion = eigenvectors[:, int(np.argmax(eigenvalues))]
    return Rotation(quaternion=mean_quaternion).as_matrix()


def reconstruct_parent_grains(
    child_orientations: OrientationSet,
    adjacency: ArrayLike,
    relationship: OrientationRelationship,
    *,
    tolerance_deg: float = 3.0,
    provenance: ProvenanceRecord | None = None,
) -> ParentGrainReconstructionResult:
    """Group child grains into parents and estimate each parent orientation.

    Purpose: the map-scale reconstruction step — decides which neighboring
    product grains descend from one parent grain and what that parent's
    orientation was, using only the orientation relationship (no parent phase
    retained in the microstructure required).

    Algorithm, in two stages: connectivity proposes, consistency disposes.

    First, an adjacency edge links two child grains when their boundary
    misorientation ``M = C_i^T C_j`` lies within ``tolerance_deg`` of the
    same-parent fingerprint ``G_c (R G_p R^T) G_c`` (see
    ``intervariant_boundary_fingerprint``) — the full rotation is matched,
    axis and angle, since a same-variant pair contributes the identity and a
    distinct-variant pair contributes ``V_i V_j^T``. Union-find over linked
    edges yields candidate clusters.

    Connectivity alone is not sufficient evidence, because two unrelated
    parents can share a boundary that genuinely lies inside the fingerprint —
    no edge test can reject those, and on a dense grain graph they are common
    enough to dominate the error (``chance_link_probability`` reports how
    common). Each candidate cluster is therefore split by agreement: every
    member proposes the parent it implies, each proposal is scored by how many
    members it explains (``C_j^T P`` near the variant-description set), and the
    best-supported proposal claims its supporters; unexplained members repeat
    the vote. A cluster spanning two parents separates, while a genuine
    single-parent cluster is returned whole.

    Each resulting group's winning proposal then seeds a quaternion-eigen-mean
    refinement over its members, and per-grain residuals to the refined parent
    are reported.

    Inputs: child grain-mean orientations (phase must match the
    relationship's child phase), an ``(m, 2)`` integer adjacency array of
    child indices, the nominal relationship, and the edge tolerance.

    Output: a ``ParentGrainReconstructionResult`` (see its ``describe()``).
    """

    if not _child_semantics_match(child_orientations, relationship):
        raise ValueError(
            "child_orientations must match the relationship child phase "
            "(or, for phase-less map data, carry the child point-group symmetry)."
        )
    grain_count = len(child_orientations)
    if grain_count == 0:
        raise ValueError("reconstruct_parent_grains requires at least one child grain.")
    edges = as_int_array(np.asarray(adjacency, dtype=np.int64).reshape(-1, 2), shape=(None, 2))
    if edges.size and (edges.min() < 0 or edges.max() >= grain_count):
        raise ValueError("adjacency indices must reference child_orientations entries.")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("adjacency must not contain self-edges.")
    if not np.isfinite(tolerance_deg) or tolerance_deg <= 0.0:
        raise ValueError("tolerance_deg must be positive and finite.")

    child_matrices = child_orientations.as_matrices()
    operators = _child_operators(relationship)
    parent_symmetry = relationship.parent_phase.symmetry
    parent_operators = (
        parent_symmetry.operators
        if parent_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)

    # Edge test: the boundary misorientation must match the same-parent
    # fingerprint as a *rotation* — axis and angle. Matching the angle alone
    # is far too permissive: for a cubic-cubic relationship the intervariant
    # angle spectrum covers so much of [0, 62] deg that a 3 deg angle-only
    # window admits over half of all entirely unrelated boundaries, which
    # merges distinct parents at map scale.
    linked = np.zeros(edges.shape[0], dtype=bool)
    if edges.shape[0]:
        # Crystal-frame boundary relative under the canonical crystal->specimen
        # orientation convention: M = C_i^T C_j.
        relative = np.einsum(
            "eji,ejk->eik",
            child_matrices[edges[:, 0]],
            child_matrices[edges[:, 1]],
            optimize=True,
        )
        # The fingerprint is closed under child symmetry on both sides, so the
        # raw relative is scored directly (no prior disorientation reduction).
        distance = boundary_fingerprint_distances_deg(relative, fingerprint)
        linked = distance <= tolerance_deg

    parent = np.arange(grain_count, dtype=np.int64)

    def find(index: int) -> int:
        root = index
        while parent[root] != root:
            root = int(parent[root])
        while parent[index] != root:
            parent[index], index = root, int(parent[index])
        return root

    for edge_index in np.flatnonzero(linked):
        left_root = find(int(edges[edge_index, 0]))
        right_root = find(int(edges[edge_index, 1]))
        if left_root != right_root:
            parent[right_root] = left_root

    roots = np.array([find(int(index)) for index in range(grain_count)], dtype=np.int64)
    unique_roots, labels = np.unique(roots, return_inverse=True)

    variants = relationship.generate_variants()
    variant_matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in variants], axis=0
    )
    descriptions = _variant_descriptions(variant_matrices, operators)
    final_labels = np.empty(grain_count, dtype=np.int64)
    deviations = np.empty(grain_count, dtype=np.float64)
    parent_list: list[np.ndarray] = []
    size_list: list[int] = []
    split_clusters = 0
    for cluster in range(unique_roots.size):
        members = np.flatnonzero(labels == cluster)
        # Connectivity proposes the cluster; single-parent consistency disposes.
        # A cluster joined through a coincidental boundary is separated here.
        groups = _vote_partition_cluster(
            child_matrices[members],
            descriptions=descriptions,
            variant_matrices=variant_matrices,
            tolerance_deg=tolerance_deg,
        )
        if len(groups) > 1:
            split_clusters += 1
        for local_indices, seed_parent in groups:
            group_members = members[local_indices]
            group_children = child_matrices[group_members]
            refined = _refine_cluster_parent(
                group_children,
                seed_parent,
                variant_matrices=variant_matrices,
                parent_operators=parent_operators,
            )
            final_labels[group_members] = len(parent_list)
            parent_list.append(refined)
            size_list.append(int(local_indices.size))
            deviations[group_members] = _member_deviations_deg(
                group_children,
                refined,
                variant_matrices=variant_matrices,
                operators=operators,
            )
    parent_matrices = np.stack(parent_list, axis=0)
    sizes = np.asarray(size_list, dtype=np.int64)
    labels = final_labels

    parent_phase = relationship.parent_phase
    parent_set = OrientationSet.from_matrices(
        parent_matrices,
        crystal_frame=parent_phase.crystal_frame,
        specimen_frame=child_orientations.specimen_frame,
        symmetry=parent_phase.symmetry,
        phase=parent_phase,
        provenance=provenance or relationship.provenance,
    )
    return ParentGrainReconstructionResult(
        relationship_name=relationship.name,
        parent_labels=labels.astype(np.int64),
        parent_orientations=parent_set,
        per_grain_deviation_deg=deviations,
        cluster_sizes=sizes,
        edges_tested=int(edges.shape[0]),
        edges_linked=int(np.count_nonzero(linked)),
        tolerance_deg=float(tolerance_deg),
        chance_link_probability=_chance_link_probability(fingerprint, float(tolerance_deg)),
        provenance=provenance or relationship.provenance,
    )


def _child_semantics_match(
    child_orientations: OrientationSet, relationship: OrientationRelationship
) -> bool:
    if phases_semantically_match(child_orientations.phase, relationship.child_phase):
        return True
    child_symmetry = relationship.child_phase.symmetry
    return (
        child_orientations.phase is None
        and child_orientations.symmetry is not None
        and child_symmetry is not None
        and normalize_point_group_symbol(child_orientations.symmetry.point_group)
        == normalize_point_group_symbol(child_symmetry.point_group)
    )


def reconstruct_parent_grains_from_graph(
    graph: GrainGraph,
    relationship: OrientationRelationship,
    *,
    tolerance_deg: float = 3.0,
    provenance: ProvenanceRecord | None = None,
) -> ParentGrainReconstructionResult:
    """Reconstruct parent grains directly from an EBSD grain graph.

    Purpose: the map-facing entry point — takes the ``GrainGraph`` produced by
    ``GrainSegmentation.grain_graph()``, uses each grain's mean orientation as
    the child orientation and each graph edge as adjacency, and runs
    ``reconstruct_parent_grains``. Result rows follow ``graph.node_grain_ids``
    order, and the returned ``grain_ids`` field records that mapping.

    Inputs: a grain graph over the transformed (child) phase and the nominal
    orientation relationship; the map's orientations must carry the child
    phase, or a symmetry whose point group matches it.

    Output: a ``ParentGrainReconstructionResult`` with ``grain_ids`` set.
    """

    segmentation = graph.segmentation
    means = segmentation.grain_mean_orientations()
    grain_ids = np.asarray(graph.node_grain_ids, dtype=np.int64)
    children = OrientationSet.from_orientations(
        [means[int(grain_id)] for grain_id in grain_ids]
    )
    node_index = {int(grain_id): index for index, grain_id in enumerate(grain_ids)}
    adjacency = np.asarray(
        [
            [node_index[int(edge.left_grain_id)], node_index[int(edge.right_grain_id)]]
            for edge in graph.edges
        ],
        dtype=np.int64,
    ).reshape(-1, 2)
    result = reconstruct_parent_grains(
        children,
        adjacency,
        relationship,
        tolerance_deg=tolerance_deg,
        provenance=provenance,
    )
    return ParentGrainReconstructionResult(
        relationship_name=result.relationship_name,
        parent_labels=result.parent_labels,
        parent_orientations=result.parent_orientations,
        per_grain_deviation_deg=result.per_grain_deviation_deg,
        cluster_sizes=result.cluster_sizes,
        edges_tested=result.edges_tested,
        edges_linked=result.edges_linked,
        tolerance_deg=result.tolerance_deg,
        chance_link_probability=result.chance_link_probability,
        grain_ids=grain_ids,
        provenance=result.provenance,
    )
