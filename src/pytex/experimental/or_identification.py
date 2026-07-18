"""Experimental OR identification from child-child boundary misorientations.

Given a fully transformed microstructure (no parent phase retained), the
boundary misorientations between child grains that share a parent are drawn
from the relationship's intervariant set — the double coset
``G_c (R S_p R^T) G_c`` of the parent symmetry conjugated by the OR rotation.
Scoring measured boundaries against each candidate's fingerprint set answers
"which orientation relationship is operative" without any parent orientation
(OR-foundation F7, first stage). The second stage — boundary-based
*refinement* of the winning rotation — lives in
``pytex.experimental.or_refinement``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core.orientation import OrientationSet, matrices_to_quaternions
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import OrientationRelationship


@dataclass(frozen=True, slots=True)
class ORIdentificationReport:
    """Ranking of candidate orientation relationships against boundary data.

    ``mean_distances_deg[i]`` is the mean angular distance of every measured
    child-child boundary misorientation to the nearest element of candidate
    ``candidate_names[i]``'s intervariant fingerprint; ``best_name`` is the
    minimizer and ``best_edge_distances_deg`` its per-edge residuals. A clear
    winner near zero (relative to orientation noise) identifies the operative
    relationship; several close candidates mean the boundary data cannot
    discriminate them.
    """

    candidate_names: tuple[str, ...]
    mean_distances_deg: np.ndarray
    median_distances_deg: np.ndarray
    best_name: str
    best_edge_distances_deg: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        names = tuple(str(name) for name in self.candidate_names)
        means = np.asarray(self.mean_distances_deg, dtype=np.float64).reshape(-1)
        medians = np.asarray(self.median_distances_deg, dtype=np.float64).reshape(-1)
        edges = np.asarray(self.best_edge_distances_deg, dtype=np.float64).reshape(-1)
        if len(names) == 0:
            raise ValueError("ORIdentificationReport requires at least one candidate.")
        if means.shape != (len(names),) or medians.shape != (len(names),):
            raise ValueError("Distance arrays must have one entry per candidate.")
        if self.best_name not in names:
            raise ValueError("best_name must be one of candidate_names.")
        if edges.size == 0:
            raise ValueError("best_edge_distances_deg must contain at least one edge.")
        for array in (means, medians, edges):
            if np.any(~np.isfinite(array)) or np.any(array < 0.0):
                raise ValueError("Distances must be finite and non-negative.")
            array.setflags(write=False)
        object.__setattr__(self, "candidate_names", names)
        object.__setattr__(self, "mean_distances_deg", means)
        object.__setattr__(self, "median_distances_deg", medians)
        object.__setattr__(self, "best_edge_distances_deg", edges)

    def describe(self) -> str:
        """Prose summary: ranking, winner margin, honest ambiguity statement."""

        order = np.argsort(self.mean_distances_deg)
        ranking = "; ".join(
            f"{self.candidate_names[int(index)]}: "
            f"{float(self.mean_distances_deg[int(index)]):.3f} deg"
            for index in order[:5]
        )
        margin = (
            float(
                self.mean_distances_deg[order[1]] - self.mean_distances_deg[order[0]]
            )
            if order.size > 1
            else float("nan")
        )
        margin_text = (
            f" The winner leads the runner-up by {margin:.3f} deg; margins comparable "
            "to the orientation noise mean the boundary data cannot discriminate the "
            "candidates."
            if order.size > 1
            else ""
        )
        return (
            f"Orientation-relationship identification from "
            f"{self.best_edge_distances_deg.size} child-child boundary "
            f"misorientation(s), no parent orientations used: best candidate "
            f"'{self.best_name}' with mean fingerprint distance "
            f"{float(np.min(self.mean_distances_deg)):.3f} deg "
            f"(per-edge max {float(np.max(self.best_edge_distances_deg)):.3f} deg). "
            f"Ranking: {ranking}.{margin_text} Experimental surface: identification "
            "among supplied candidates only — for boundary-based refinement of the "
            "rotation itself, pass the winner to "
            "refine_orientation_relationship_from_boundaries."
        )


def _fingerprint_set(relationship: OrientationRelationship) -> np.ndarray:
    """Deduplicated double coset G_c (R S_p R^T) G_c of same-parent boundaries."""

    rotation = relationship.parent_to_child_rotation.as_matrix()
    parent_symmetry = relationship.parent_phase.symmetry
    child_symmetry = relationship.child_phase.symmetry
    parent_ops = (
        parent_symmetry.operators
        if parent_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    child_ops = (
        child_symmetry.operators
        if child_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    conjugated = np.einsum("ij,pjk,lk->pil", rotation, parent_ops, rotation, optimize=True)
    products = np.einsum(
        "aij,pjk,bkl->apbil", child_ops, conjugated, child_ops, optimize=True
    ).reshape(-1, 3, 3)
    quaternions = matrices_to_quaternions(products)
    # q and -q are one rotation: canonicalize sign on the largest component.
    pivot = np.argmax(np.abs(quaternions), axis=1)
    signs = np.sign(quaternions[np.arange(quaternions.shape[0]), pivot])
    keys = np.round(quaternions * signs[:, None], 8)
    _, unique_indices = np.unique(keys, axis=0, return_index=True)
    fingerprint = np.ascontiguousarray(products[unique_indices])
    fingerprint.setflags(write=False)
    return fingerprint


def identify_orientation_relationship(
    child_orientations: OrientationSet,
    adjacency: ArrayLike,
    candidates: tuple[OrientationRelationship, ...],
    *,
    provenance: ProvenanceRecord | None = None,
) -> ORIdentificationReport:
    """Rank candidate ORs against measured child-child boundary misorientations.

    Purpose: identifies the operative orientation relationship of a fully
    transformed microstructure without any parent orientation — the boundary
    misorientations of same-parent child grains must populate the candidate's
    intervariant fingerprint, so the candidate with the smallest mean
    distance of every boundary to its fingerprint set wins.

    Inputs: child grain-mean orientations (phase must match every candidate's
    child phase), an ``(m, 2)`` adjacency array of same-parent child pairs,
    and the candidate relationships (e.g. from
    ``standard_fcc_bcc_relationships(...).relationships``).

    Output: an ``ORIdentificationReport`` (see its ``describe()``).
    """

    if len(candidates) == 0:
        raise ValueError("identify_orientation_relationship requires candidates.")
    from pytex.core.lattice import phases_semantically_match

    for candidate in candidates:
        if not phases_semantically_match(
            child_orientations.phase, candidate.child_phase
        ):
            raise ValueError(
                f"child_orientations.phase must match candidate '{candidate.name}' "
                "child phase."
            )
    edges = as_int_array(np.asarray(adjacency, dtype=np.int64).reshape(-1, 2), shape=(None, 2))
    if edges.shape[0] == 0:
        raise ValueError("identify_orientation_relationship requires at least one edge.")
    grain_count = len(child_orientations)
    if edges.min() < 0 or edges.max() >= grain_count:
        raise ValueError("adjacency indices must reference child_orientations entries.")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("adjacency must not contain self-edges.")
    child_matrices = child_orientations.as_matrices()
    relatives = np.einsum(
        "eji,ejk->eik",
        child_matrices[edges[:, 0]],
        child_matrices[edges[:, 1]],
        optimize=True,
    )
    means: list[float] = []
    medians: list[float] = []
    per_edge_best: dict[str, np.ndarray] = {}
    for candidate in candidates:
        fingerprint = _fingerprint_set(candidate)
        # angle(M F^T) needs only trace(M F^T) = sum(M * F), elementwise.
        traces = np.einsum("eij,kij->ek", relatives, fingerprint, optimize=True)
        cosines = np.clip((traces.max(axis=1) - 1.0) * 0.5, -1.0, 1.0)
        distances = np.degrees(np.arccos(cosines))
        means.append(float(np.mean(distances)))
        medians.append(float(np.median(distances)))
        per_edge_best[candidate.name] = distances
    best_index = int(np.argmin(means))
    best_name = candidates[best_index].name
    return ORIdentificationReport(
        candidate_names=tuple(candidate.name for candidate in candidates),
        mean_distances_deg=np.asarray(means, dtype=np.float64),
        median_distances_deg=np.asarray(medians, dtype=np.float64),
        best_name=best_name,
        best_edge_distances_deg=per_edge_best[best_name],
        provenance=provenance,
    )
