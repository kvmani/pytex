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
    """Settings governing how candidate parent orientations are scored.

    Attributes
    ----------
    reduction : str
        How per-child residuals are reduced to one score: ``"mean"``
        (default), ``"median"`` (robust to a few misindexed children), or
        ``"max"`` (worst case).
    symmetry_aware : bool
        Use symmetry-reduced disorientation angles (default).
    ambiguity_tolerance_deg : float
        Candidates scoring within this of the best are reported as
        ambiguous. Setting it to zero does not remove ambiguity, it hides it.
    provenance : ProvenanceRecord, optional
    """

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
    """Which transformation variant best explains each child orientation.

    Purpose
    -------
    Variant selection is the observable signature of the transformation
    mechanism: a random variant distribution indicates no selection, while a
    skewed one points to stress, boundary, or prior-deformation effects.

    Attributes
    ----------
    variant_indices : np.ndarray
        One-based variant index per child; strictly positive.
    scores_deg : np.ndarray
        The residual angle of the assignment, per child. A large residual
        means the child is not well explained by *any* variant, which is a
        different conclusion from "assigned to variant k".
    provenance : ProvenanceRecord, optional
    """

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
    """The outcome of scoring candidate parent orientations, with its ambiguity.

    Purpose
    -------
    Reports which candidate parent best explains a set of measured child
    orientations, and — as importantly — whether the answer is determined at
    all. With few children, several parents typically explain the data
    equally well; ``ambiguous_indices`` records that rather than letting the
    top score pass as a determination.

    Attributes
    ----------
    record : PhaseTransformationRecord
        The parent/child data and the orientation relationship used.
    candidate_parents : OrientationSet
        The candidates that were scored.
    scores_deg : np.ndarray
        One residual score per candidate.
    best_index : int
    best_score_deg : float
    ambiguous_indices : tuple of int
        Candidates within the configured tolerance of the best. More than one
        entry means the reconstruction is unresolved.
    reduction : str
    symmetry_aware : bool
        The settings used, recorded so the result is reproducible.
    provenance : ProvenanceRecord, optional
    """

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
        """Whether more than one candidate parent scores within the tolerance.

        An ambiguous reconstruction must not be reported as a determination: the
        data admit several parents, which for a single child orientation is the
        normal case rather than a failure.
        """

        return len(self.ambiguous_indices) > 1

    def best_parent_orientation(self) -> Orientation:
        """The highest-scoring candidate parent orientation.

        Check :attr:`is_ambiguous` before relying on it — the best candidate is
        only meaningful when the selection is unambiguous.
        """

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
    """A named collection of orientation relationships.

    Purpose
    -------
    The lookup for standard relationships — Kurdjumov-Sachs,
    Nishiyama-Wassermann, Bain, Burgers, Pitsch — so a fitted relationship
    can be compared against the literature by name rather than by
    hand-entered angles. Names are required to be unique.

    Attributes
    ----------
    relationships : tuple of OrientationRelationship
    provenance : ProvenanceRecord, optional
    """

    relationships: tuple[OrientationRelationship, ...]
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "relationships", tuple(self.relationships))
        names = [relationship.name for relationship in self.relationships]
        if len(set(names)) != len(names):
            raise ValueError("OrientationRelationshipCatalog relationship names must be unique.")

    def names(self) -> tuple[str, ...]:
        """Names of the catalogued orientation relationships, in catalog order.
        """

        return tuple(relationship.name for relationship in self.relationships)

    def get(self, name: str) -> OrientationRelationship:
        """One catalogued orientation relationship by name.

        Raises ``KeyError`` for an unknown name.
        """

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
    """Choose the parent orientation best explaining a transformation record.

    Purpose
    -------
    Given measured child orientations and a set of candidate parents, score
    each candidate by how well it reproduces the children under the record's
    orientation relationship, and report the best with an explicit ambiguity
    verdict.

    Method and limits
    -----------------
    Scores *given* candidates; it does not search orientation space, so the
    candidate set bounds what can be found. With few child orientations
    several parents typically explain the data equally well, which is why the
    report carries :attr:`~ParentReconstructionReport.is_ambiguous` rather
    than presenting the top score as an answer.

    Parameters
    ----------
    record : PhaseTransformationRecord
        The parent/child orientations and the orientation relationship.
    candidate_parents : OrientationSet
        Candidates to score; must be non-empty.
    config : ParentReconstructionConfig, optional
        Scoring and ambiguity-tolerance settings.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    ParentReconstructionReport
        Scores, the best index, the set of ambiguous indices, and a
        :meth:`~ParentReconstructionReport.describe` prose summary.
    """

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


def standard_cubic_cubic_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The cubic-to-cubic relationships that are not transformation ORs.

    Purpose
    -------
    Two cubic phases are related by more than the austenite-to-ferrite family.
    The two that a cubic-cubic ranking must be able to name are:

    ``cube_on_cube``
        Parallel axes, the identity rotation. A coherent cubic precipitate on a
        cubic matrix -- gamma-prime in a superalloy, TiN in ferrite, an
        epitaxial film on a cubic substrate.
    ``fcc_twin``
        The coherent ``{111}`` twin, 60 degrees about ``<111>`` after symmetry
        reduction: the Sigma 3 boundary of an annealing or deformation twin.

    Why they belong in the default cubic-cubic catalog
    --------------------------------------------------
    Without them a measured Sigma 3 boundary is reported as some number of
    degrees away from Bain, which is true and useless. Neither is a
    transformation relationship, and neither displaces one: they are added to
    the fcc->bcc family rather than replacing it, and a genuine
    martensite relationship still ranks against Bain, Kurdjumov-Sachs,
    Nishiyama-Wassermann, Greninger-Troiano and Pitsch exactly as before.

    Note on the twin's two phases. Matrix and twin are the same material, and
    :class:`~pytex.core.transformation.OrientationRelationship` rejects a
    relationship between a phase and itself, so the twin is built from the two
    phase objects handed in. Passing the same phase twice raises; pass a
    distinguishable child, e.g. ``Phase("nickel (twin)", ...)``.
    """

    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_cube_on_cube_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_fcc_twin_correspondence(
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

    Returns the Pitsch-Schrader relationship, the inverse Burgers
    relationship, and the Potter relationship bound to the given hexagonal
    parent and cubic child phases.
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
            OrientationRelationship.from_potter_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
        ),
        provenance=provenance,
    )


def standard_ferrite_cementite_relationships(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog:
    """The standard named ferrite->cementite orientation-relationship catalog.

    Returns the Bagaryatsky and Isaichev relationships (tempered-martensite /
    pearlite carbide precipitation class, Pnma cementite setting) bound to
    the given cubic parent and orthorhombic child phases.
    """

    return OrientationRelationshipCatalog(
        relationships=(
            OrientationRelationship.from_bagaryatsky_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
            OrientationRelationship.from_isaichev_correspondence(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ),
        ),
        provenance=provenance,
    )


#: Which standard catalog serves which (parent, child) crystal-system pair.
#:
#: One table rather than a chain of conditionals, so the dispatch is auditable
#: and extending it is a one-line change. Cubic-to-cubic cannot distinguish an
#: fcc parent from a bcc one by point group alone, so it resolves to the full
#: fcc->bcc family (Bain, Kurdjumov-Sachs, Nishiyama-Wassermann,
#: Greninger-Troiano, Pitsch); that is the intended reading of a cubic-cubic
#: transformation catalog and is stated in `default_relationship_catalog`. It
#: also carries the two cubic-cubic relationships that are not transformations
#: at all -- cube-on-cube and the coherent {111} twin -- because a measured
#: Sigma 3 boundary reported as a deviation from Bain is a named answer the
#: ranking could have given and did not.
_CATALOG_DISPATCH: dict[tuple[str, str], tuple[str, ...]] = {
    ("cubic", "cubic"): (
        "standard_fcc_bcc_relationships",
        "standard_cubic_cubic_relationships",
    ),
    ("cubic", "hexagonal"): (
        "standard_bcc_hcp_relationships",
        "standard_fcc_hcp_relationships",
    ),
    ("hexagonal", "cubic"): ("standard_hcp_bcc_relationships",),
    ("cubic", "orthorhombic"): ("standard_ferrite_cementite_relationships",),
}


def default_relationship_catalog(
    *,
    parent_phase: Phase,
    child_phase: Phase,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipCatalog | None:
    """The standard OR catalog for a parent/child phase pair, chosen by crystal system.

    Purpose: lets callers such as
    `pytex.core.transformation.characterize_orientation_relationship` rank a
    measured relationship against the named relationships that are actually
    plausible for the two phases, without the caller having to know which
    builder applies.

    Inputs: the parent and child phases. The choice depends only on their
    crystal systems, per the `_CATALOG_DISPATCH` table.

    Output: an `OrientationRelationshipCatalog`, or ``None`` when no standard
    catalog covers the pair — in which case a fitted relationship is reported
    without a name rather than being forced onto an inapplicable list.

    Note the cubic-to-cubic entry returns the fcc->bcc family. Point-group
    symmetry cannot tell an fcc phase from a bcc one, so a cubic-cubic
    transformation is assumed to be of the austenite->ferrite class; supply an
    explicit catalog when it is not.
    """

    parent_system = parent_phase.symmetry.to_point_group().crystal_system
    child_system = child_phase.symmetry.to_point_group().crystal_system
    builder_names = _CATALOG_DISPATCH.get((parent_system, child_system))
    if builder_names is None:
        return None
    relationships: list[OrientationRelationship] = []
    for builder_name in builder_names:
        builder = globals()[builder_name]
        relationships.extend(
            builder(
                parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
            ).relationships
        )
    return OrientationRelationshipCatalog(
        relationships=tuple(relationships), provenance=provenance
    )


__all__ = [
    "OrientationRelationshipCatalog",
    "ParentReconstructionConfig",
    "ParentReconstructionReport",
    "VariantSelectionReport",
    "default_relationship_catalog",
    "reconstruct_parent_orientation",
    "select_variants",
    "standard_bcc_hcp_relationships",
    "standard_cubic_cubic_relationships",
    "standard_fcc_bcc_relationships",
    "standard_fcc_hcp_relationships",
    "standard_ferrite_cementite_relationships",
    "standard_hcp_bcc_relationships",
]
