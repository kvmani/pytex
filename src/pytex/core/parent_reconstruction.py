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

    def variant_frequencies(self, *, variant_count: int | None = None) -> np.ndarray:
        """Occurrence count of each variant index (1-based).

        Returns an array of length ``variant_count`` (defaults to the largest
        selected index) whose entry ``k`` counts selections of variant
        ``k + 1`` — the variant-selection histogram.
        """

        highest = int(self.variant_indices.max()) if self.variant_indices.size else 0
        total = highest if variant_count is None else int(variant_count)
        if total < highest:
            raise ValueError("variant_count must cover the largest selected variant index.")
        counts = np.bincount(self.variant_indices, minlength=total + 1)[1:]
        return np.ascontiguousarray(counts)

    def describe(self) -> str:
        """Prose summary: children assigned, variants used, residual statistics."""

        frequencies = self.variant_frequencies()
        order = np.argsort(frequencies)[::-1]
        top = [
            f"V{int(index) + 1} ({int(frequencies[index])}x)"
            for index in order[:5]
            if frequencies[index] > 0
        ]
        distinct = int(np.count_nonzero(frequencies))
        return (
            f"Variant selection over {self.variant_indices.size} child orientation(s): "
            f"{distinct} distinct variant(s) selected; most frequent: {', '.join(top)}. "
            f"Residual misorientation to the assigned variant: mean "
            f"{float(np.mean(self.scores_deg)):.3f} deg, max "
            f"{float(np.max(self.scores_deg)):.3f} deg (child-symmetry-reduced angles). "
            "Uniform frequencies indicate no variant selection; strong imbalance "
            "indicates selection relative to the equal-probability baseline."
        )


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

    def describe(self) -> str:
        """Prose summary: best candidate, score, reduction rule, ambiguity."""

        ambiguity = (
            f"AMBIGUOUS: {len(self.ambiguous_indices)} candidates lie within the "
            "ambiguity tolerance of the best score; treat the selection as unresolved."
            if self.is_ambiguous
            else "The selection is unambiguous within the ambiguity tolerance."
        )
        symmetry_text = (
            "child-symmetry-reduced" if self.symmetry_aware else "raw (not symmetry-reduced)"
        )
        return (
            f"Parent reconstruction over {len(self.candidate_parents)} candidate parent "
            f"orientation(s) for record '{self.record.name}': best candidate index "
            f"{self.best_index} with a {self.reduction} residual of "
            f"{self.best_score_deg:.3f} deg across the observed children "
            f"({symmetry_text} angles). {ambiguity}"
        )


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


def select_variants(
    record: PhaseTransformationRecord,
    *,
    symmetry_aware: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> VariantSelectionReport:
    """Assign each observed child orientation to its nearest transformation variant.

    For every child orientation in ``record``, predicts the child orientation
    of each variant of the record's orientation relationship applied to the
    parent orientation, and selects the variant with the smallest
    (child-symmetry-reduced) misorientation. Returns a
    `VariantSelectionReport` with the 1-based variant index and the residual
    angle per child; ``report.variant_frequencies()`` gives the
    variant-selection histogram.
    """

    variants = record.orientation_relationship.generate_variants()
    children = record.child_orientations
    if len(children) == 0:
        raise ValueError("select_variants requires at least one child orientation.")
    variant_matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in variants], axis=0
    )
    parent_matrix = record.parent_orientation.rotation.as_matrix()
    # Canonical crystal->specimen convention: C = P @ V^T per variant.
    predicted_matrices = np.einsum(
        "ij,vkj->vik", parent_matrix, variant_matrices, optimize=True
    )
    predicted = OrientationSet.from_matrices(
        predicted_matrices,
        crystal_frame=children.crystal_frame,
        specimen_frame=children.specimen_frame,
        symmetry=children.symmetry,
        phase=children.phase,
    )
    angles_deg = np.degrees(
        children.misorientation_angles_to(predicted, symmetry_aware=symmetry_aware)
    )
    best_columns = np.argmin(angles_deg, axis=1)
    indices = np.array(
        [variants[int(column)].variant_index for column in best_columns], dtype=np.int64
    )
    scores = angles_deg[np.arange(len(children)), best_columns]
    return VariantSelectionReport(
        variant_indices=indices, scores_deg=scores, provenance=provenance
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


def standard_fcc_hcp_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The standard named fcc->hcp orientation-relationship catalog.

    Returns the Shoji-Nishiyama relationship (austenite -> epsilon-martensite
    class) bound to the given cubic parent and hexagonal child phases.
    """

    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_shoji_nishiyama_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
        ),
        provenance=provenance,
    )


def standard_hcp_bcc_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The standard named hcp->bcc orientation-relationship catalog.

    Returns the Pitsch-Schrader relationship and the inverse Burgers
    relationship bound to the given hexagonal parent and cubic child phases.
    """

    inverse_burgers = OrientationRelationship.from_burgers_correspondence(
        parent_phase=child_phase, child_phase=parent_phase, provenance=provenance
    ).inverse(name="burgers_inverse")
    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_pitsch_schrader_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            inverse_burgers,
        ),
        provenance=provenance,
    )


__all__ = [
    "OrientationRelationshipCatalog",
    "ParentReconstructionConfig",
    "ParentReconstructionReport",
    "VariantSelectionReport",
    "reconstruct_parent_orientation",
    "select_variants",
    "standard_bcc_hcp_relationships",
    "standard_fcc_bcc_relationships",
    "standard_fcc_hcp_relationships",
    "standard_hcp_bcc_relationships",
]
