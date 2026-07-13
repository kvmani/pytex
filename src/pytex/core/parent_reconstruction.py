from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core.lattice import Phase
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import OrientationRelationship, PhaseTransformationRecord
from pytex.experimental.phase_transformation import score_parent_orientations

ReductionMode = Literal["mean", "median", "max"]


@dataclass(frozen=True, slots=True)
class ParentReconstructionConfig:
    reduction: ReductionMode = "mean"
    symmetry_aware: bool = True
    ambiguity_tolerance_deg: float = 1.0
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.reduction not in {"mean", "median", "max"}:
            raise ValueError("reduction must be 'mean', 'median', or 'max'.")
        if not np.isfinite(self.ambiguity_tolerance_deg) or self.ambiguity_tolerance_deg < 0.0:
            raise ValueError("ambiguity_tolerance_deg must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class VariantSelectionReport:
    variant_indices: np.ndarray
    scores_deg: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        indices = np.asarray(self.variant_indices, dtype=np.int64).reshape(-1)
        scores = np.asarray(self.scores_deg, dtype=np.float64).reshape(-1)
        if indices.shape != scores.shape:
            raise ValueError("variant_indices and scores_deg must have the same shape.")
        if np.any(indices <= 0):
            raise ValueError("variant_indices must be strictly positive.")
        if np.any(~np.isfinite(scores)) or np.any(scores < 0.0):
            raise ValueError("scores_deg must be finite and non-negative.")
        indices = np.ascontiguousarray(indices)
        scores = np.ascontiguousarray(scores)
        indices.setflags(write=False)
        scores.setflags(write=False)
        object.__setattr__(self, "variant_indices", indices)
        object.__setattr__(self, "scores_deg", scores)


@dataclass(frozen=True, slots=True)
class ParentReconstructionReport:
    record: PhaseTransformationRecord
    candidate_parents: OrientationSet
    scores_deg: np.ndarray
    best_index: int
    best_score_deg: float
    ambiguous_indices: tuple[int, ...]
    reduction: ReductionMode
    symmetry_aware: bool
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        scores = np.asarray(self.scores_deg, dtype=np.float64).reshape(-1)
        if scores.shape != (len(self.candidate_parents),):
            raise ValueError("scores_deg must contain one score per candidate parent.")
        if np.any(~np.isfinite(scores)) or np.any(scores < 0.0):
            raise ValueError("scores_deg must be finite and non-negative.")
        if not 0 <= self.best_index < len(self.candidate_parents):
            raise ValueError("best_index is out of range.")
        if not np.isclose(float(scores[self.best_index]), self.best_score_deg, atol=1e-10):
            raise ValueError("best_score_deg must match scores_deg[best_index].")
        scores = np.ascontiguousarray(scores, dtype=np.float64)
        scores.setflags(write=False)
        object.__setattr__(self, "scores_deg", scores)
        object.__setattr__(self, "ambiguous_indices", tuple(int(i) for i in self.ambiguous_indices))

    @property
    def is_ambiguous(self) -> bool:
        return len(self.ambiguous_indices) > 1

    def best_parent_orientation(self) -> Orientation:
        return self.candidate_parents[self.best_index]


@dataclass(frozen=True, slots=True)
class OrientationRelationshipCatalog:
    relationships: tuple[OrientationRelationship, ...]
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "relationships", tuple(self.relationships))
        names = [relationship.name for relationship in self.relationships]
        if len(set(names)) != len(names):
            raise ValueError("OrientationRelationshipCatalog relationship names must be unique.")

    def names(self) -> tuple[str, ...]:
        return tuple(relationship.name for relationship in self.relationships)

    def get(self, name: str) -> OrientationRelationship:
        for relationship in self.relationships:
            if relationship.name == name:
                return relationship
        raise KeyError(name)


def reconstruct_parent_orientation(
    record: PhaseTransformationRecord,
    candidate_parents: OrientationSet,
    *,
    config: ParentReconstructionConfig | None = None,
    provenance: ProvenanceRecord | None = None,
) -> ParentReconstructionReport:
    if len(candidate_parents) == 0:
        raise ValueError("candidate_parents must contain at least one orientation.")
    reconstruction_config = ParentReconstructionConfig() if config is None else config
    result = score_parent_orientations(
        record,
        candidate_parents,
        symmetry_aware=reconstruction_config.symmetry_aware,
        reduction=reconstruction_config.reduction,
        provenance=provenance or reconstruction_config.provenance,
    )
    scores = result.scores_deg
    best = int(np.argmin(scores))
    cutoff = float(scores[best] + reconstruction_config.ambiguity_tolerance_deg)
    ambiguous = tuple(int(index) for index in np.flatnonzero(scores <= cutoff))
    return ParentReconstructionReport(
        record=record,
        candidate_parents=candidate_parents,
        scores_deg=scores,
        best_index=best,
        best_score_deg=float(scores[best]),
        ambiguous_indices=ambiguous,
        reduction=reconstruction_config.reduction,
        symmetry_aware=reconstruction_config.symmetry_aware,
        provenance=provenance or reconstruction_config.provenance or result.provenance,
    )


def standard_fcc_bcc_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The standard named fcc->bcc orientation-relationship catalog.

    Returns Bain, Kurdjumov-Sachs, Nishiyama-Wassermann, Greninger-Troiano,
    and Pitsch bound to the given cubic parent (austenite-like) and child
    (ferrite/martensite-like) phases, resolvable by their default names.
    """

    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_bain_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_kurdjumov_sachs_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_nishiyama_wassermann_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_greninger_troiano_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_pitsch_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
        ),
        provenance=provenance,
    )


def standard_bcc_hcp_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The standard named bcc->hcp orientation-relationship catalog.

    Returns the Burgers relationship (beta->alpha titanium/zirconium class)
    bound to the given cubic parent and hexagonal child phases.
    """

    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_burgers_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
        ),
        provenance=provenance,
    )


__all__ = [
    "OrientationRelationshipCatalog",
    "ParentReconstructionConfig",
    "ParentReconstructionReport",
    "VariantSelectionReport",
    "reconstruct_parent_orientation",
    "standard_bcc_hcp_relationships",
    "standard_fcc_bcc_relationships",
]
