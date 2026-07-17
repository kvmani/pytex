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
    _reduced_pair_disorientation_angles,
)
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import (
    OrientationRelationship,
    intervariant_misorientation_angles_deg,
)


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
    """

    relationship_name: str
    parent_labels: np.ndarray
    parent_orientations: OrientationSet
    per_grain_deviation_deg: np.ndarray
    cluster_sizes: np.ndarray
    edges_tested: int
    edges_linked: int
    tolerance_deg: float
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
        return (
            f"Parent-grain reconstruction under orientation relationship "
            f"'{self.relationship_name}': {self.parent_labels.size} child grain(s) grouped "
            f"into {self.parent_count} parent(s) using {self.edges_linked} of "
            f"{self.edges_tested} adjacency edge(s) whose boundary disorientations matched "
            f"the intervariant fingerprint within {self.tolerance_deg:.2f} deg. Residual of "
            f"each child to its assigned parent (best variant, child-symmetry-reduced): mean "
            f"{self.mean_deviation_deg:.3f} deg, max {self.max_deviation_deg:.3f} deg."
            f"{singleton_note} Experimental surface: literature-fixture validation is still "
            "ahead; treat groupings as hypotheses at map scale."
        )


def _child_operators(relationship: OrientationRelationship) -> np.ndarray:
    symmetry = relationship.child_phase.symmetry
    if symmetry is None:
        return np.eye(3, dtype=np.float64)[None, :, :]
    return np.asarray(symmetry.operators, dtype=np.float64)


def _reference_angles_deg(relationship: OrientationRelationship) -> np.ndarray:
    """Intervariant fingerprint: 0 (same variant) plus all distinct pair angles."""

    table = intervariant_misorientation_angles_deg(relationship)
    upper = table[np.triu_indices(table.shape[0], k=1)]
    reference = np.unique(np.round(np.concatenate([[0.0], upper]), 6))
    return np.asarray(reference, dtype=np.float64)


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

    Algorithm: an adjacency edge links two child grains when their
    child-symmetry-reduced disorientation angle matches the relationship's
    intervariant fingerprint (including 0 deg for same-variant neighbors)
    within ``tolerance_deg``; union-find over linked edges yields parent
    clusters; per cluster, the candidate parents ``V_k^T C_first`` generated
    from its first child are scored against every member (minimum-over-
    variants disorientation), and the best-scoring candidate becomes the
    parent estimate with per-grain residuals reported.

    Inputs: child grain-mean orientations (phase must match the
    relationship's child phase), an ``(m, 2)`` integer adjacency array of
    child indices, the nominal relationship, and the edge tolerance.

    Output: a ``ParentGrainReconstructionResult`` (see its ``describe()``).
    """

    if not phases_semantically_match(child_orientations.phase, relationship.child_phase):
        raise ValueError("child_orientations.phase must match the relationship child phase.")
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
    reference = _reference_angles_deg(relationship)

    # Edge test: boundary disorientation matches the intervariant fingerprint.
    linked = np.zeros(edges.shape[0], dtype=bool)
    if edges.shape[0]:
        relative = np.einsum(
            "eij,ekj->eik",
            child_matrices[edges[:, 0]],
            child_matrices[edges[:, 1]],
            optimize=True,
        )
        angles = np.degrees(
            _reduced_pair_disorientation_angles(relative, operators, operators)
        )
        distance = np.min(np.abs(angles[:, None] - reference[None, :]), axis=1)
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
    parent_matrices = np.empty((unique_roots.size, 3, 3), dtype=np.float64)
    deviations = np.empty(grain_count, dtype=np.float64)
    sizes = np.empty(unique_roots.size, dtype=np.int64)
    for cluster in range(unique_roots.size):
        members = np.flatnonzero(labels == cluster)
        sizes[cluster] = members.size
        cluster_children = child_matrices[members]
        # Candidate parents from the first member: P_k = V_k^T C_first.
        candidates = np.einsum(
            "vji,jk->vik", variant_matrices, cluster_children[0], optimize=True
        )
        # Predicted children for every (candidate, variant): V_l P_k.
        predicted = np.einsum("lij,vjk->vlik", variant_matrices, candidates, optimize=True)
        relative = np.einsum(
            "nij,vlkj->nvlik", cluster_children, predicted, optimize=True
        )
        member_count, candidate_count, variant_count = relative.shape[:3]
        angles = np.degrees(
            _reduced_pair_disorientation_angles(
                relative.reshape(-1, 3, 3), operators, operators
            )
        ).reshape(member_count, candidate_count, variant_count)
        per_member = np.min(angles, axis=2)
        best_candidate = int(np.argmin(np.mean(per_member, axis=0)))
        parent_matrices[cluster] = candidates[best_candidate]
        deviations[members] = per_member[:, best_candidate]

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
        provenance=provenance or relationship.provenance,
    )
