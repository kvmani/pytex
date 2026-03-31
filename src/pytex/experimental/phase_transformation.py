from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core.lattice import Phase
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.core.transformation import PhaseTransformationRecord

ReductionMode = Literal["mean", "median", "max"]


def _phase_semantically_matches(left: Phase | None, right: Phase | None) -> bool:
    if left is None or right is None:
        return left is right
    return (
        left.name == right.name
        and left.crystal_frame == right.crystal_frame
        and left.lattice == right.lattice
        and left.symmetry == right.symmetry
    )


def _rotation_angles_from_matrices(matrices: np.ndarray) -> np.ndarray:
    traces = np.trace(matrices, axis1=-2, axis2=-1)
    cos_theta = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    return np.asarray(np.arccos(cos_theta), dtype=np.float64)


def _disorientation_angles(
    relative_matrices: np.ndarray,
    *,
    left_symmetry: SymmetrySpec | None,
    right_symmetry: SymmetrySpec | None,
) -> np.ndarray:
    left_ops = (
        left_symmetry.operators
        if left_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    right_ops = (
        right_symmetry.operators
        if right_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    left_applied = np.einsum("aij,mnjk->mnaik", left_ops, relative_matrices, optimize=True)
    candidates = np.einsum(
        "mnaij,bkj->mnabik",
        left_applied,
        right_ops,
        optimize=True,
    )
    angles = _rotation_angles_from_matrices(candidates.reshape(-1, 3, 3)).reshape(
        relative_matrices.shape[0],
        relative_matrices.shape[1],
        left_ops.shape[0] * right_ops.shape[0],
    )
    return np.asarray(np.min(angles, axis=2), dtype=np.float64)


def _resolved_child_symmetry(record: PhaseTransformationRecord) -> SymmetrySpec | None:
    if record.child_orientations.phase is not None:
        return record.child_orientations.phase.symmetry
    return record.child_orientations.symmetry


def _variant_rotation_matrices(record: PhaseTransformationRecord) -> np.ndarray:
    child_count = len(record.child_orientations)
    if record.variant_indices is None:
        matrix = record.orientation_relationship.parent_to_child_rotation.as_matrix()
        return np.repeat(matrix[None, :, :], child_count, axis=0)
    variant_lookup = {
        variant.variant_index: variant.parent_to_child_rotation.as_matrix()
        for variant in record.orientation_relationship.generate_variants()
    }
    variant_indices = np.asarray(record.variant_indices, dtype=np.int64)
    unique_indices, inverse = np.unique(variant_indices, return_inverse=True)
    missing = [
        int(index) for index in unique_indices.tolist() if int(index) not in variant_lookup
    ]
    if missing:
        raise ValueError(
            "PhaseTransformationRecord.variant_indices contain values not produced by "
            "OrientationRelationship.generate_variants(): "
            + ", ".join(str(value) for value in missing)
        )
    unique_matrices = np.stack(
        [variant_lookup[int(index)] for index in unique_indices.tolist()],
        axis=0,
    )
    return np.asarray(unique_matrices[inverse], dtype=np.float64)


def _validate_candidate_parent_set(
    record: PhaseTransformationRecord,
    candidate_parents: OrientationSet,
) -> None:
    reference_parent = record.parent_orientation
    if candidate_parents.crystal_frame != reference_parent.crystal_frame:
        raise ValueError("Candidate parent orientations must use the parent crystal frame.")
    if candidate_parents.specimen_frame != reference_parent.specimen_frame:
        raise ValueError("Candidate parent orientations must use the parent specimen frame.")
    if not _phase_semantically_matches(candidate_parents.phase, reference_parent.phase):
        raise ValueError("Candidate parent orientations must match the parent phase semantics.")
    if candidate_parents.symmetry != reference_parent.symmetry:
        raise ValueError("Candidate parent orientations must match the parent symmetry.")


def _reduce_candidate_scores(
    angles_rad: np.ndarray,
    *,
    reduction: ReductionMode,
) -> np.ndarray:
    if reduction == "mean":
        return np.asarray(np.mean(angles_rad, axis=1), dtype=np.float64)
    if reduction == "median":
        return np.asarray(np.median(angles_rad, axis=1), dtype=np.float64)
    return np.asarray(np.max(angles_rad, axis=1), dtype=np.float64)


def _predicted_child_matrices(
    record: PhaseTransformationRecord,
    candidate_parent_matrices: np.ndarray,
) -> np.ndarray:
    variant_matrices = _variant_rotation_matrices(record)
    return np.asarray(
        np.einsum("nij,mjk->mnik", variant_matrices, candidate_parent_matrices, optimize=True),
        dtype=np.float64,
    )


@dataclass(frozen=True, slots=True)
class ParentReconstructionResult:
    record: PhaseTransformationRecord
    candidate_parents: OrientationSet
    scores_rad: np.ndarray
    reduction: ReductionMode
    symmetry_aware: bool
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.reduction not in {"mean", "median", "max"}:
            raise ValueError("reduction must be 'mean', 'median', or 'max'.")
        scores = np.asarray(self.scores_rad, dtype=np.float64).reshape(-1)
        if scores.shape != (len(self.candidate_parents),):
            raise ValueError("scores_rad must contain one score per candidate parent orientation.")
        scores = np.ascontiguousarray(scores)
        scores.setflags(write=False)
        object.__setattr__(self, "scores_rad", scores)

    @property
    def scores_deg(self) -> np.ndarray:
        values = np.asarray(np.rad2deg(self.scores_rad), dtype=np.float64)
        values.setflags(write=False)
        return values

    @property
    def best_index(self) -> int:
        return int(np.argmin(self.scores_rad))

    @property
    def best_score_rad(self) -> float:
        return float(self.scores_rad[self.best_index])

    @property
    def best_score_deg(self) -> float:
        return float(np.rad2deg(self.best_score_rad))

    def best_parent_orientation(self) -> Orientation:
        return self.candidate_parents[self.best_index]

    def predicted_child_orientations(self, *, index: int | None = None) -> OrientationSet:
        candidate_index = self.best_index if index is None else int(index)
        candidate_matrix = self.candidate_parents.as_matrices()[
            candidate_index : candidate_index + 1
        ]
        predicted = _predicted_child_matrices(self.record, candidate_matrix)[0]
        return OrientationSet.from_matrices(
            predicted,
            crystal_frame=self.record.child_orientations.crystal_frame,
            specimen_frame=self.record.child_orientations.specimen_frame,
            symmetry=self.record.child_orientations.symmetry,
            phase=self.record.child_orientations.phase,
            provenance=self.provenance or self.record.provenance,
        )


def score_parent_orientations(
    record: PhaseTransformationRecord,
    candidate_parents: OrientationSet,
    *,
    symmetry_aware: bool = True,
    reduction: ReductionMode = "mean",
    provenance: ProvenanceRecord | None = None,
) -> ParentReconstructionResult:
    """Score candidate parent orientations against observed child orientations.

    This function is intentionally experimental. It provides a bounded parent-reconstruction
    scoring primitive without promoting full parent reconstruction into the stable public API.

    Parameters
    ----------
    record:
        Transformation record containing the observed child orientations and any selected variant
        assignments.
    candidate_parents:
        Parent-orientation hypotheses to score. The set must match the parent frame, phase, and
        symmetry semantics recorded in ``record``.
    symmetry_aware:
        When ``True``, child-child comparison uses disorientation under the child symmetry.
        When ``False``, raw misorientation angles are used.
    reduction:
        Reduction applied across the per-child angular residuals for each candidate parent.
        Supported values are ``"mean"``, ``"median"``, and ``"max"``.
    provenance:
        Optional provenance attached to the result object and to any predicted-child orientations
        materialized from it.
    """

    if len(candidate_parents) == 0:
        raise ValueError("candidate_parents must contain at least one orientation.")
    _validate_candidate_parent_set(record, candidate_parents)
    parent_matrices = candidate_parents.as_matrices()
    predicted_child_matrices = _predicted_child_matrices(record, parent_matrices)
    observed_child_matrices = record.child_orientations.as_matrices()
    relative = np.asarray(
        np.einsum(
            "nij,mnjk->mnik",
            observed_child_matrices,
            np.swapaxes(predicted_child_matrices, -1, -2),
            optimize=True,
        ),
        dtype=np.float64,
    )
    if symmetry_aware:
        child_symmetry = _resolved_child_symmetry(record)
        residuals = _disorientation_angles(
            relative,
            left_symmetry=child_symmetry,
            right_symmetry=child_symmetry,
        )
    else:
        residuals = _rotation_angles_from_matrices(relative)
    scores = _reduce_candidate_scores(residuals, reduction=reduction)
    return ParentReconstructionResult(
        record=record,
        candidate_parents=candidate_parents,
        scores_rad=scores,
        reduction=reduction,
        symmetry_aware=symmetry_aware,
        provenance=provenance or record.provenance or candidate_parents.provenance,
    )
