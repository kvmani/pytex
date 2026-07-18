"""Experimental boundary-based refinement of an orientation-relationship rotation.

Second stage of OR-foundation F7: given child-child boundary misorientations
of a fully transformed microstructure (no parent orientations) and a nominal
orientation relationship, recover the operative OR rotation. Same-parent
boundary misorientations populate the double coset ``G_c (R S_p R^T) G_c``;
the refinement alternates (1) assigning each measured boundary to its nearest
coset element at the current rotation estimate with (2) a least-squares
rotation update with the assignments held fixed, until the assignments are
stable and the update step is below tolerance.

The rotation is only identifiable up to the symmetry ambiguity of the coset
structure, so distances between rotations are reported symmetry-reduced.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import least_squares

from pytex.core._arrays import as_int_array
from pytex.core.orientation import OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import (
    OrientationRelationship,
    _symmetry_reduced_angle_between_deg,
)


@dataclass(frozen=True, slots=True)
class ORRefinementReport:
    """Boundary-refined orientation relationship with fit diagnostics.

    ``relationship`` carries the refined rotation (named ``<nominal>_refined``;
    defining parallelisms are deliberately not carried over — they describe
    the nominal relationship, not the refined rotation). ``edge_distances_deg``
    are the per-boundary residuals to the refined fingerprint;
    ``rotation_update_deg`` is the symmetry-reduced angle between the nominal
    and refined rotations.
    """

    relationship: OrientationRelationship
    nominal_name: str
    initial_mean_distance_deg: float
    refined_mean_distance_deg: float
    edge_distances_deg: np.ndarray
    iterations: int
    converged: bool
    rotation_update_deg: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        edges = np.ascontiguousarray(
            np.asarray(self.edge_distances_deg, dtype=np.float64).reshape(-1)
        )
        if edges.size == 0:
            raise ValueError("edge_distances_deg must contain at least one edge.")
        if np.any(~np.isfinite(edges)) or np.any(edges < 0.0):
            raise ValueError("edge_distances_deg must be finite and non-negative.")
        edges.setflags(write=False)
        object.__setattr__(self, "edge_distances_deg", edges)

    def describe(self) -> str:
        """Prose summary: refinement result, residuals, honest caveats."""

        return (
            f"Boundary-based refinement of orientation relationship "
            f"'{self.nominal_name}' from {self.edge_distances_deg.size} "
            f"child-child boundary misorientation(s), no parent orientations "
            f"used: mean fingerprint distance improved from "
            f"{self.initial_mean_distance_deg:.3f} deg to "
            f"{self.refined_mean_distance_deg:.3f} deg "
            f"(per-edge max {float(np.max(self.edge_distances_deg)):.3f} deg) "
            f"after {self.iterations} iteration(s) "
            f"({'converged' if self.converged else 'NOT converged'}). The "
            f"refined rotation sits {self.rotation_update_deg:.3f} deg "
            f"(symmetry-reduced) from the nominal. Experimental surface: the "
            "rotation is identifiable only up to the coset symmetry ambiguity, "
            "and a nominal relationship far from the operative one may "
            "converge to a wrong local assignment."
        )


def _symmetry_operator_stacks(
    relationship: OrientationRelationship,
) -> tuple[np.ndarray, np.ndarray]:
    identity = np.eye(3, dtype=np.float64)[None, :, :]
    parent_symmetry = relationship.parent_phase.symmetry
    child_symmetry = relationship.child_phase.symmetry
    parent_ops = parent_symmetry.operators if parent_symmetry is not None else identity
    child_ops = child_symmetry.operators if child_symmetry is not None else identity
    return np.asarray(parent_ops, dtype=np.float64), np.asarray(child_ops, dtype=np.float64)


def _assign_coset_elements(
    relatives: np.ndarray,
    rotation: np.ndarray,
    parent_ops: np.ndarray,
    child_ops: np.ndarray,
    *,
    block_size: int = 128,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest coset element A (R S R^T) B per edge -> (a, s, b) indices, distances."""

    conjugated = np.einsum("ij,pjk,lk->pil", rotation, parent_ops, rotation, optimize=True)
    elements = np.einsum(
        "aij,pjk,bkl->apbil", child_ops, conjugated, child_ops, optimize=True
    )
    flat = elements.reshape(-1, 3, 3)
    edge_count = relatives.shape[0]
    assignments = np.empty((edge_count, 3), dtype=np.int64)
    distances = np.empty(edge_count, dtype=np.float64)
    shape = (child_ops.shape[0], parent_ops.shape[0], child_ops.shape[0])
    for start in range(0, edge_count, block_size):
        stop = min(start + block_size, edge_count)
        traces = np.einsum(
            "eij,kij->ek", relatives[start:stop], flat, optimize=True
        )
        best = np.argmax(traces, axis=1)
        cosines = np.clip((traces[np.arange(best.size), best] - 1.0) * 0.5, -1.0, 1.0)
        distances[start:stop] = np.degrees(np.arccos(cosines))
        assignments[start:stop] = np.stack(np.unravel_index(best, shape), axis=1)
    return assignments, distances


def _chordal_residuals(
    x: np.ndarray,
    estimate: np.ndarray,
    symmetry_parent: np.ndarray,
    collapsed: np.ndarray,
) -> np.ndarray:
    """Per-edge 2 sin(theta/2) residuals of the left-perturbed rotation."""

    angle = float(np.linalg.norm(x))
    if angle < 1e-16:
        delta = np.eye(3, dtype=np.float64)
    else:
        axis = x / angle
        skew = np.array(
            [
                [0.0, -axis[2], axis[1]],
                [axis[2], 0.0, -axis[0]],
                [-axis[1], axis[0], 0.0],
            ],
            dtype=np.float64,
        )
        delta = np.eye(3) + np.sin(angle) * skew + (1.0 - np.cos(angle)) * (skew @ skew)
    rotation = delta @ estimate
    conjugated = np.einsum(
        "eij,jk->eik",
        np.einsum("ij,ejk->eik", rotation, symmetry_parent, optimize=True),
        rotation.T,
        optimize=True,
    )
    cosines = np.clip(
        (np.einsum("eij,eij->e", collapsed, conjugated, optimize=True) - 1.0) * 0.5,
        -1.0,
        1.0,
    )
    # 2 sin(theta/2): monotone in the angle and smooth at zero.
    return np.asarray(np.sqrt(2.0 * np.maximum(0.0, 1.0 - cosines)), dtype=np.float64)


def refine_orientation_relationship_from_boundaries(
    child_orientations: OrientationSet,
    adjacency: ArrayLike,
    nominal: OrientationRelationship,
    *,
    max_iterations: int = 30,
    convergence_tol_deg: float = 1e-6,
    provenance: ProvenanceRecord | None = None,
) -> ORRefinementReport:
    """Refine an OR rotation against child-child boundary misorientations.

    Purpose: recovers the operative orientation-relationship rotation of a
    fully transformed microstructure from boundary misorientations alone —
    the F7 second stage that ``identify_orientation_relationship`` (the
    identification-only first stage) points to. Same-parent boundary
    misorientations populate the double coset ``G_c (R S_p R^T) G_c`` of the
    operative rotation ``R``; starting from the nominal rotation, the method
    alternates nearest-coset-element assignment with a bounded least-squares
    rotation update (smooth ``2 sin(theta/2)`` chordal residuals, so exact
    data is a regular point) until the assignments are stable.

    Inputs: child grain-mean orientations (phase must match the nominal's
    child phase), an ``(m, 2)`` adjacency array of same-parent child pairs,
    and the nominal relationship (typically the
    ``identify_orientation_relationship`` winner).

    Output: an ``ORRefinementReport`` (see its ``describe()``).
    """

    from pytex.core.lattice import phases_semantically_match

    if not phases_semantically_match(child_orientations.phase, nominal.child_phase):
        raise ValueError("child_orientations.phase must match nominal.child_phase.")
    edges = as_int_array(
        np.asarray(adjacency, dtype=np.int64).reshape(-1, 2), shape=(None, 2)
    )
    if edges.shape[0] == 0:
        raise ValueError(
            "refine_orientation_relationship_from_boundaries requires at least one edge."
        )
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
    parent_ops, child_ops = _symmetry_operator_stacks(nominal)
    estimate = nominal.parent_to_child_rotation.as_matrix()
    _, initial_distances = _assign_coset_elements(
        relatives, estimate, parent_ops, child_ops
    )
    initial_mean = float(np.mean(initial_distances))

    previous_assignments: np.ndarray | None = None
    converged = False
    iterations = 0
    distances = initial_distances
    for _ in range(max_iterations):
        iterations += 1
        assignments, distances = _assign_coset_elements(
            relatives, estimate, parent_ops, child_ops
        )
        symmetry_left = child_ops[assignments[:, 0]]
        symmetry_parent = parent_ops[assignments[:, 1]]
        symmetry_right = child_ops[assignments[:, 2]]
        # trace(M (A R S R^T B)^T) = sum((A^T M B^T) * (R S R^T)) elementwise.
        collapsed = np.einsum(
            "eji,ejk,elk->eil", symmetry_left, relatives, symmetry_right, optimize=True
        )
        solution = least_squares(
            _chordal_residuals,
            np.zeros(3),
            args=(estimate, symmetry_parent, collapsed),
            method="lm",
        )
        step_angle_deg = float(np.degrees(np.linalg.norm(solution.x)))
        angle = float(np.linalg.norm(solution.x))
        estimate = Rotation.from_axis_angle(
            solution.x / angle if angle > 1e-16 else np.array([0.0, 0.0, 1.0]),
            angle,
        ).as_matrix() @ estimate
        stable = previous_assignments is not None and bool(
            np.array_equal(assignments, previous_assignments)
        )
        previous_assignments = assignments
        if stable and step_angle_deg <= convergence_tol_deg:
            converged = True
            break
    _, distances = _assign_coset_elements(relatives, estimate, parent_ops, child_ops)

    refined_rotation = Rotation.from_matrix(estimate).canonicalized()
    refined = OrientationRelationship(
        name=f"{nominal.name}_refined",
        parent_phase=nominal.parent_phase,
        child_phase=nominal.child_phase,
        parent_to_child_rotation=refined_rotation,
        provenance=provenance or nominal.provenance,
    )
    update = _symmetry_reduced_angle_between_deg(
        estimate,
        nominal.parent_to_child_rotation.as_matrix(),
        child_operators=child_ops,
        parent_operators=parent_ops,
    )
    return ORRefinementReport(
        relationship=refined,
        nominal_name=nominal.name,
        initial_mean_distance_deg=initial_mean,
        refined_mean_distance_deg=float(np.mean(distances)),
        edge_distances_deg=distances,
        iterations=iterations,
        converged=converged,
        rotation_update_deg=update,
        provenance=provenance or nominal.provenance,
    )
