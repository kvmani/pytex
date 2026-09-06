"""Deciding which of several candidate phases a measured powder pattern is.

Indexing (:mod:`pytex.diffraction.xrd_indexing`) answers a question about one
phase: *how well does this cell account for these lines?* Identification asks a
different one: *given these several candidates, which accounts for them best,
and is the winner actually convincing?* The second question is not the first
asked repeatedly, for two reasons.

**A ranking needs criteria that a wrong phase can fail.** The figures of merit
of indexing -- de Wolff's ``M_N``, Smith and Snyder's ``F_N`` -- reward a cell
that puts lines where lines were seen. They say nothing about a strong measured
peak the candidate does not explain at all, nor about a strong line the
candidate insists on that the scan does not show, and those two are the
evidence that actually separates phases. A face-centred and a body-centred cell
of similar size both land *some* lines correctly; what distinguishes them is
which lines are absent.

**A ranking needs to be able to say "none of these".** A comparison always has
a winner. A user who uploads four CIFs and gets back "the best match is the
third" has learned nothing unless the report also says whether the third is a
good match in absolute terms, and whether it is meaningfully better than the
fourth. Both statements are made here, explicitly, as
:attr:`PhaseIdentification.is_conclusive` and
:attr:`PhaseIdentification.is_decisive`.

The four criteria
-----------------

Each candidate is indexed against the measured peaks, and then scored on four
quantities that are independent of one another and each bounded to ``[0, 1]``.
They are reported individually as well as combined, because which one a
candidate fails is diagnostic in a way its total is not.

``explained_intensity_fraction``
    The share of the measured *integrated intensity* carried by peaks the
    candidate indexed. Intensity-weighted rather than counted, because a strong
    unindexed peak is decisive evidence against a candidate while a weak one is
    routinely a trace impurity or an artefact. This is the criterion that
    detects "the sample contains something this phase is not".

``completeness``
    The share of the candidate's own strong reflections, inside the measured
    angular range, that were actually observed. This is the criterion that
    detects "this phase would have shown a line here, and the scan is flat".
    It is what separates two cells that differ by a centring, since a centring
    is a statement about which lines are *absent*.

``position_score``
    ``1 - <|Delta 2 theta|> / tolerance``, clipped at zero: how far inside the
    matching window the assigned lines actually landed, rather than merely
    whether they landed inside it. A candidate whose lines sit at the edge of
    the tolerance throughout has the wrong cell dimensions even though every
    peak was formally indexed.

``intensity_agreement``
    ``1 - (1/2) sum |I_obs - I_calc|`` over the indexed reflections with each
    set normalized to unit sum, which is one minus the Bray-Curtis
    dissimilarity and so is scale-free and bounded. This is the criterion that
    separates two phases with nearly the same cell but different atomic bases
    -- an ordered and a disordered structure, say, or the same framework with a
    different cation.

    It is deliberately given the least weight. Measured powder intensities are
    the least trustworthy quantity in the pattern: preferred orientation,
    microabsorption, extinction and a coarse powder all move them, sometimes by
    a factor of several, without moving a single peak position. A candidate
    should not be rejected on intensity alone, and this weighting is why it is
    not.

The combined score is their weighted mean, over whichever criteria are defined
for that candidate. The weights are a declared, overridable parameter rather
than a hidden constant, because they encode a judgement about evidence rather
than a law of diffraction, and a laboratory whose specimens are strongly
textured is right to lower ``intensity_agreement`` further.

What this is not
----------------

This is not a search of a powder-diffraction database. Hanawalt, Rinn and
Frevel's method -- and every search-match system descended from it -- indexes
hundreds of thousands of patterns by their three strongest lines so that
candidates can be *retrieved*. Here the candidates are already in hand, because
the user chose them, and the problem is only to rank them. Retrieval and
ranking are separate problems, and the retrieval half needs a licensed database
this library does not ship.

Nor does this quantify a mixture. When several phases each explain part of the
pattern, the honest answer is the ranked list plus the unexplained peaks, and
the quantitative next step is a Rietveld refinement of the phases together
(:mod:`pytex.diffraction.rietveld`), not a bigger score.

References
----------
Hanawalt, J. D., Rinn, H. W. & Frevel, L. K., *Ind. Eng. Chem. Anal. Ed.* **10**
(1938) 457-512, doi:10.1021/ac50125a001 -- the search-match method by
characteristic ``d`` spacings, of which the criteria here are the
already-have-the-candidates case.

Smith, G. S. & Snyder, R. L., *J. Appl. Crystallogr.* **12** (1979) 60-65,
doi:10.1107/S002188987901178X -- ``F_N``, reported alongside each candidate.

de Wolff, P. M., *J. Appl. Crystallogr.* **1** (1968) 108-113,
doi:10.1107/S002188986800508X -- ``M_20``, likewise.

Dollase, W. A., *J. Appl. Crystallogr.* **19** (1986) 267-272,
doi:10.1107/S0021889886089458 -- the March model of preferred orientation, and
the reason measured intensities carry the least weight of the four criteria.

Gates-Rector, S. & Blanton, T., *Powder Diffr.* **34** (2019) 352-360,
doi:10.1017/S0885715619000812 -- the Powder Diffraction File, the database this
module deliberately does not attempt to replace.

Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Ch. 14 -- chemical analysis by diffraction.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np

from pytex.core.lattice import Phase
from pytex.diffraction.xrd import RadiationSpec, generate_powder_reflections
from pytex.diffraction.xrd_indexing import PeakIndexing, index_peaks
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import PeakTable, detect_and_fit_peaks

PHASE_CANDIDATE_SCORE_SCHEMA = "pytex.diffraction.phase_candidate_score"
PHASE_IDENTIFICATION_SCHEMA = "pytex.diffraction.phase_identification"

#: The four criteria, in the order they are reported everywhere.
CRITERION_NAMES: tuple[str, ...] = (
    "explained_intensity_fraction",
    "completeness",
    "position_score",
    "intensity_agreement",
)

#: Default weight of each criterion in the combined score. Position and
#: explained intensity dominate because peak positions are the trustworthy part
#: of a powder pattern; intensity agreement is lowest because preferred
#: orientation moves measured intensities without moving a single position.
DEFAULT_CRITERION_WEIGHTS: Mapping[str, float] = MappingProxyType(
    {
        "explained_intensity_fraction": 0.40,
        "completeness": 0.25,
        "position_score": 0.20,
        "intensity_agreement": 0.15,
    }
)

_CITATION_HANAWALT = (
    "Hanawalt, Rinn & Frevel, Ind. Eng. Chem. Anal. Ed. 10 (1938) 457, "
    "doi:10.1021/ac50125a001."
)
_CITATION_SMITH_SNYDER = (
    "Smith & Snyder, J. Appl. Crystallogr. 12 (1979) 60, doi:10.1107/S002188987901178X."
)
_CITATION_DOLLASE = (
    "Dollase, J. Appl. Crystallogr. 19 (1986) 267, doi:10.1107/S0021889886089458."
)
_CITATION_CULLITY = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 14."
)


def _normalized_weights(weights: Mapping[str, float] | None) -> Mapping[str, float]:
    """Validate a criterion weighting and return it as a read-only mapping."""

    supplied = DEFAULT_CRITERION_WEIGHTS if weights is None else weights
    unknown = sorted(set(supplied) - set(CRITERION_NAMES))
    if unknown:
        raise ValueError(
            "Unknown scoring criterion(s): "
            + ", ".join(unknown)
            + ". Known criteria are: "
            + ", ".join(CRITERION_NAMES)
            + "."
        )
    resolved = {name: float(supplied.get(name, 0.0)) for name in CRITERION_NAMES}
    if any(value < 0.0 or not np.isfinite(value) for value in resolved.values()):
        raise ValueError("Criterion weights must be finite and non-negative.")
    if sum(resolved.values()) <= 0.0:
        raise ValueError("At least one criterion weight must be strictly positive.")
    return MappingProxyType(resolved)


@dataclass(frozen=True, slots=True)
class PhaseCandidateScore:
    """One candidate phase, indexed against the measurement and scored.

    Purpose
    -------
    Hold the evidence for and against a single candidate in a form that can be
    read on its own. The combined :attr:`score` orders the candidates, but the
    four criteria are what say *why* a candidate lost, and they are carried
    individually for exactly that reason.

    Attributes
    ----------
    phase_name : str
        The candidate's name, as the user knows it.
    source : str
        Where the candidate came from -- a CIF file name, a built-in catalogue
        entry -- so a ranking can be traced back to the files that produced it.
    indexing : PeakIndexing | None
        The assignment of measured peaks to this candidate's reflections.
        ``None`` when the candidate could not be indexed at all, in which case
        :attr:`rejection` says why and every criterion is zero.
    rejection : str
        Empty for a candidate that was scored; otherwise the reason indexing
        was impossible, which is itself a decisive result.
    explained_intensity_fraction : float
        Share of the measured integrated intensity carried by indexed peaks.
    completeness : float
        Share of this candidate's strong reflections, inside the measured
        range, that were observed.
    position_score : float
        ``1 - <|Delta 2 theta|> / tolerance``, clipped at zero. ``nan`` when
        nothing was indexed.
    intensity_agreement : float
        One minus the Bray-Curtis dissimilarity between observed and calculated
        relative intensities over the indexed reflections. ``nan`` when fewer
        than two reflections were indexed, since a similarity between two
        one-element distributions is identically one and would be misleading.
    strongest_unexplained_fraction : float
        Integrated intensity of the strongest *unindexed* peak, as a fraction
        of the strongest measured peak. A value near one means the candidate
        fails to explain the most prominent feature of the pattern.
    strongest_unobserved_relative_intensity : float
        Calculated relative intensity of the strongest reflection this
        candidate predicts inside the measured range but that was not observed.
    strong_line_threshold : float
        Relative-intensity threshold above which a calculated line counted
        towards :attr:`completeness`.
    weights : Mapping[str, float]
        The criterion weighting the combined score was formed with.
    """

    phase_name: str
    source: str = ""
    indexing: PeakIndexing | None = None
    rejection: str = ""
    explained_intensity_fraction: float = 0.0
    completeness: float = 0.0
    position_score: float = float("nan")
    intensity_agreement: float = float("nan")
    strongest_unexplained_fraction: float = 0.0
    strongest_unobserved_relative_intensity: float = 0.0
    strong_line_threshold: float = 0.05
    weights: Mapping[str, float] = field(default_factory=lambda: DEFAULT_CRITERION_WEIGHTS)

    def __post_init__(self) -> None:
        if not self.phase_name.strip():
            raise ValueError("PhaseCandidateScore.phase_name must be non-empty.")
        if self.indexing is None and not self.rejection.strip():
            raise ValueError(
                "PhaseCandidateScore with no indexing must state a rejection reason."
            )
        object.__setattr__(self, "weights", _normalized_weights(self.weights))

    @property
    def criteria(self) -> Mapping[str, float]:
        """Return the four criterion values keyed by name, in report order."""

        return MappingProxyType(
            {
                "explained_intensity_fraction": float(self.explained_intensity_fraction),
                "completeness": float(self.completeness),
                "position_score": float(self.position_score),
                "intensity_agreement": float(self.intensity_agreement),
            }
        )

    @property
    def score(self) -> float:
        """Return the combined score in ``[0, 1]``, higher being a better match.

        The weighted mean of whichever criteria are defined, with the weights
        renormalized over those. A criterion that is undefined for this
        candidate -- intensity agreement on a single indexed line, say -- is
        omitted rather than treated as zero, because "not measurable here" and
        "measured and bad" are different findings and only the second is
        evidence against the candidate.
        """

        if self.indexing is None:
            return 0.0
        total_weight = 0.0
        total = 0.0
        for name, value in self.criteria.items():
            if not np.isfinite(value):
                continue
            weight = float(self.weights[name])
            total += weight * float(value)
            total_weight += weight
        return 0.0 if total_weight <= 0.0 else float(total / total_weight)

    @property
    def indexed_count(self) -> int:
        """Return the number of measured peaks this candidate indexed."""

        return 0 if self.indexing is None else self.indexing.indexed_count

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this candidate."""

        merit_m = merit_f = None
        if self.indexing is not None:
            value_m, count_m = self.indexing.figure_of_merit_m()
            value_f, count_f = self.indexing.figure_of_merit_f()
            merit_m = {"value": float(value_m), "count": int(count_m)}
            merit_f = {"value": float(value_f), "count": int(count_f)}
        return {
            "schema": PHASE_CANDIDATE_SCORE_SCHEMA,
            "phase_name": self.phase_name,
            "source": self.source,
            "rejection": self.rejection,
            "score": self.score,
            "criteria": dict(self.criteria),
            "weights": dict(self.weights),
            "indexed_count": self.indexed_count,
            "unindexed_count": (
                0 if self.indexing is None else len(self.indexing.unindexed_peaks)
            ),
            "indexed_fraction": (
                0.0 if self.indexing is None else self.indexing.indexed_fraction
            ),
            "mean_absolute_delta_two_theta_deg": (
                float("nan")
                if self.indexing is None
                else self.indexing.mean_absolute_delta_two_theta_deg
            ),
            "strongest_unexplained_fraction": float(self.strongest_unexplained_fraction),
            "strongest_unobserved_relative_intensity": float(
                self.strongest_unobserved_relative_intensity
            ),
            "strong_line_threshold": float(self.strong_line_threshold),
            "figure_of_merit_m": merit_m,
            "figure_of_merit_f": merit_f,
            "indexing": None if self.indexing is None else self.indexing.to_json(),
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this candidate."""

        origin = f" (from {self.source})" if self.source else ""
        if self.indexing is None:
            return (
                f"Candidate '{self.phase_name}'{origin} could not be indexed against this "
                f"pattern and scores 0.000: {self.rejection}"
            )
        merit_m, count_m = self.indexing.figure_of_merit_m()
        merit_f, count_f = self.indexing.figure_of_merit_f()
        unindexed = len(self.indexing.unindexed_peaks)
        parts = [
            f"Candidate '{self.phase_name}'{origin} scores {self.score:.3f}. It indexed "
            f"{self.indexed_count} of {self.indexed_count + unindexed} measured peaks, "
            f"accounting for {100.0 * self.explained_intensity_fraction:.1f} per cent of the "
            f"measured integrated intensity, with a mean absolute position discrepancy of "
            f"{self.indexing.mean_absolute_delta_two_theta_deg:.4f} degrees 2*theta against a "
            f"{self.indexing.tolerance_deg:.3f} degree tolerance "
            f"(position score {self.position_score:.3f}). "
            f"{100.0 * self.completeness:.1f} per cent of the reflections it predicts above "
            f"{100.0 * self.strong_line_threshold:.0f} per cent relative intensity inside the "
            f"measured range were observed."
        ]
        if np.isfinite(self.intensity_agreement):
            parts.append(
                f" Observed and calculated relative intensities agree to "
                f"{self.intensity_agreement:.3f} on the bounded similarity used here."
            )
        else:
            parts.append(
                " Too few reflections were indexed for the intensity comparison to mean "
                "anything, so it was left out of the score rather than counted as a failure."
            )
        if self.strongest_unexplained_fraction >= 0.2:
            parts.append(
                f" A measured peak carrying "
                f"{100.0 * self.strongest_unexplained_fraction:.0f} per cent of the strongest "
                "peak's intensity is left unexplained, which is the signature of a second "
                "phase rather than of a poor fit."
            )
        if self.strongest_unobserved_relative_intensity >= 0.5:
            parts.append(
                f" A reflection this candidate predicts at "
                f"{100.0 * self.strongest_unobserved_relative_intensity:.0f} per cent relative "
                "intensity was not observed; absent strong lines argue against a candidate as "
                "forcefully as extra ones do, and are how a centring is recognised."
            )
        parts.append(
            f" de Wolff M_{count_m} = {merit_m:.1f} and Smith-Snyder F_{count_f} = "
            f"{merit_f:.1f} for this assignment."
        )
        return "".join(parts)


@dataclass(frozen=True, slots=True)
class PhaseIdentification:
    """A ranking of candidate phases against one measured pattern.

    Purpose
    -------
    Answer the operator's actual question -- *which of these is it?* -- and
    qualify the answer with the two statements that make it usable: whether the
    winner is good in absolute terms (:attr:`is_conclusive`) and whether it is
    better than the runner-up by enough to matter (:attr:`is_decisive`).

    Attributes
    ----------
    name : str
        A human name for the result.
    candidates : tuple[PhaseCandidateScore, ...]
        Every candidate, in descending order of score.
    minimum_score : float
        The score at or above which the best match is called conclusive.
    decisive_margin : float
        The score gap over the runner-up at or above which the ranking is
        called decisive.
    peak_count : int
        Number of measured peaks the candidates were compared against.
    settings : Mapping[str, float | str]
        The settings that produced the result.
    """

    name: str
    candidates: tuple[PhaseCandidateScore, ...]
    minimum_score: float = 0.55
    decisive_margin: float = 0.05
    peak_count: int = 0
    settings: Mapping[str, float | str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("PhaseIdentification.name must be non-empty.")
        if not self.candidates:
            raise ValueError("PhaseIdentification needs at least one candidate.")
        # Ties are broken by explained intensity and then by name, so that the
        # same inputs always produce the same ranking regardless of the order
        # the candidates were supplied in.
        ordered = tuple(
            sorted(
                self.candidates,
                key=lambda item: (
                    -item.score,
                    -item.explained_intensity_fraction,
                    item.phase_name,
                ),
            )
        )
        object.__setattr__(self, "candidates", ordered)
        object.__setattr__(self, "settings", MappingProxyType(dict(self.settings)))

    def __len__(self) -> int:
        return len(self.candidates)

    def __iter__(self) -> Any:
        return iter(self.candidates)

    @property
    def best(self) -> PhaseCandidateScore:
        """Return the highest-scoring candidate."""

        return self.candidates[0]

    @property
    def runner_up(self) -> PhaseCandidateScore | None:
        """Return the second-placed candidate, or ``None`` if there is only one."""

        return self.candidates[1] if len(self.candidates) > 1 else None

    @property
    def margin(self) -> float:
        """Return the best candidate's score lead over the runner-up.

        ``nan`` when only one candidate was offered, because a comparison of
        one has no margin and reporting zero would read as a tie.
        """

        second = self.runner_up
        return float("nan") if second is None else float(self.best.score - second.score)

    @property
    def is_conclusive(self) -> bool:
        """Return whether the best candidate explains the pattern well enough to accept."""

        return bool(self.best.score >= self.minimum_score)

    @property
    def is_decisive(self) -> bool:
        """Return whether the best candidate beats the runner-up by enough to matter.

        True by convention when only one candidate was offered: nothing was
        chosen between, so nothing is ambiguous about the choice.
        """

        margin = self.margin
        return True if not np.isfinite(margin) else bool(margin >= self.decisive_margin)

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this identification."""

        return {
            "schema": PHASE_IDENTIFICATION_SCHEMA,
            "name": self.name,
            "peak_count": int(self.peak_count),
            "candidate_count": len(self.candidates),
            "best_phase_name": self.best.phase_name,
            "best_score": self.best.score,
            "margin": self.margin,
            "minimum_score": float(self.minimum_score),
            "decisive_margin": float(self.decisive_margin),
            "is_conclusive": self.is_conclusive,
            "is_decisive": self.is_decisive,
            "candidates": [item.to_json() for item in self.candidates],
            "settings": {key: value for key, value in self.settings.items()},
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this identification."""

        best = self.best
        second = self.runner_up
        head = (
            f"Identification '{self.name}' compared {len(self.candidates)} candidate phases "
            f"against {self.peak_count} measured peaks. The best match is "
            f"'{best.phase_name}' at a score of {best.score:.3f}"
        )
        if second is None:
            head += ", the only candidate offered, so this is a check on one phase rather "
            head += "than a choice between several."
        else:
            head += (
                f", ahead of '{second.phase_name}' at {second.score:.3f}, a margin of "
                f"{self.margin:.3f}."
            )

        if not self.is_conclusive:
            verdict = (
                f" That is below the {self.minimum_score:.2f} threshold for an accepted "
                "identification, so the honest reading is that none of the candidates offered "
                "accounts for this pattern. Either the phase present is not among them, the "
                "specimen is a mixture, or the matching tolerance is too tight for the "
                "instrument's zero-point and displacement errors."
            )
        elif not self.is_decisive:
            verdict = (
                f" The margin over the runner-up is below {self.decisive_margin:.2f}, so the "
                "top two are not distinguished by this scan. Distinguishing them needs the "
                "evidence a comparison of positions cannot supply: a longer count at high "
                "angle where their calculated lines diverge, a different wavelength, or "
                "chemistry."
            )
        else:
            verdict = (
                " The match is both above the acceptance threshold and clear of the "
                "runner-up, so the identification stands on this scan alone."
            )

        leftover = ""
        if best.indexing is not None and best.indexing.unindexed_peaks:
            leftover = (
                f" {len(best.indexing.unindexed_peaks)} measured peaks remain unexplained by "
                "the best match; if any is strong, the specimen holds a second phase and the "
                "quantitative step is a multi-phase Rietveld refinement rather than a further "
                "search."
            )

        caution = (
            " Ranking rests mainly on peak positions: preferred orientation, microabsorption "
            "and a coarse powder all move measured intensities substantially without moving a "
            "position, which is why intensity agreement carries the least weight of the four "
            "criteria."
        )

        return (
            head
            + verdict
            + leftover
            + caution
            + f" {_CITATION_HANAWALT} {_CITATION_SMITH_SNYDER} {_CITATION_DOLLASE} "
            + _CITATION_CULLITY
        )


def _as_named_candidates(
    candidates: Mapping[str, Phase] | Sequence[Phase] | Sequence[tuple[str, Phase]],
) -> tuple[tuple[str, Phase], ...]:
    """Normalize the several accepted candidate spellings to ``(name, phase)`` pairs."""

    pairs: list[tuple[str, Any]] = []
    if isinstance(candidates, Mapping):
        pairs = [(str(key), value) for key, value in candidates.items()]
    else:
        for index, entry in enumerate(candidates):
            if isinstance(entry, tuple):
                if len(entry) != 2:
                    raise ValueError(
                        "A candidate given as a tuple must be a (name, phase) pair; "
                        f"entry {index + 1} has {len(entry)} elements."
                    )
                pairs.append((str(entry[0]), entry[1]))
            else:
                default = f"candidate {index + 1}"
                pairs.append((str(getattr(entry, "name", None) or default), entry))
    if not pairs:
        raise ValueError("identify_phase needs at least one candidate phase.")
    for name, phase in pairs:
        if not isinstance(phase, Phase):
            raise TypeError(
                f"Candidate '{name}' is not a pytex.core.lattice.Phase. Load a CIF with "
                "Phase.from_cif or Phase.from_cif_string first."
            )
    seen: dict[str, int] = {}
    unique: list[tuple[str, Phase]] = []
    for name, phase in pairs:
        # Two uploaded CIFs can carry the same chemical formula as their name.
        # Silently ranking two rows both called "Fe2O3" would make the result
        # unreadable, so the duplicates are numbered instead.
        count = seen.get(name, 0)
        seen[name] = count + 1
        unique.append((name if count == 0 else f"{name} ({count + 1})", phase))
    return tuple(unique)


def _bray_curtis_similarity(observed: np.ndarray, calculated: np.ndarray) -> float:
    """Return ``1 - `` the Bray-Curtis dissimilarity of two intensity vectors.

    Both are normalized to unit sum first, so the result is invariant to the
    arbitrary scale of a measured count rate and of a kinematic relative
    intensity, and is bounded to ``[0, 1]``.
    """

    observed_sum = float(np.sum(observed))
    calculated_sum = float(np.sum(calculated))
    if observed_sum <= 0.0 or calculated_sum <= 0.0:
        return float("nan")
    return float(
        1.0 - 0.5 * np.sum(np.abs(observed / observed_sum - calculated / calculated_sum))
    )


def _score_candidate(
    *,
    phase_name: str,
    phase: Phase,
    source: str,
    table: PeakTable,
    radiation: RadiationSpec,
    tolerance_deg: float,
    max_index: int,
    minimum_relative_intensity: float,
    strong_line_threshold: float,
    weights: Mapping[str, float],
    total_measured_intensity: float,
    strongest_measured_intensity: float,
) -> PhaseCandidateScore:
    """Index one candidate against the peak table and score it on the four criteria."""

    try:
        indexing = index_peaks(
            table,
            phase,
            radiation=radiation,
            tolerance_deg=tolerance_deg,
            max_index=max_index,
            minimum_relative_intensity=minimum_relative_intensity,
            phase_name=phase_name,
            name=f"{table.name} vs {phase_name}",
        )
    except ValueError as error:
        # A candidate that predicts nothing in the measured range, or nothing
        # above the intensity floor, is not an error in the comparison -- it is
        # a decisive result about that candidate. Recording it keeps the other
        # candidates rankable, which is the whole point of offering several.
        return PhaseCandidateScore(
            phase_name=phase_name,
            source=source,
            indexing=None,
            rejection=str(error),
            weights=weights,
            strong_line_threshold=strong_line_threshold,
        )

    explained = sum(item.peak.integrated_intensity for item in indexing.reflections)
    if total_measured_intensity > 0.0:
        explained_fraction = float(explained / total_measured_intensity)
    else:  # pragma: no cover - a table of zero-area peaks is pathological
        explained_fraction = indexing.indexed_fraction
    explained_fraction = float(np.clip(explained_fraction, 0.0, 1.0))

    if indexing.unindexed_peaks and strongest_measured_intensity > 0.0:
        strongest_unexplained = float(
            max(peak.integrated_intensity for peak in indexing.unindexed_peaks)
            / strongest_measured_intensity
        )
    else:
        strongest_unexplained = 0.0

    # Completeness is judged strictly inside the measured span. A strong line
    # predicted beyond the last measured point was not looked for, and counting
    # it as missing would penalize a candidate for the operator's scan range.
    observed_angles = table.two_theta_deg
    low = max(float(observed_angles[0]), 1.0e-3)
    high = min(float(observed_angles[-1]), 179.999)
    predicted = generate_powder_reflections(
        phase,
        radiation=radiation,
        two_theta_range_deg=(low, high),
        max_index=max_index,
    )
    completeness = 0.0
    strongest_unobserved = 0.0
    if predicted:
        strongest_predicted = max(item.intensity for item in predicted)
        if strongest_predicted > 0.0:
            strong = tuple(
                item
                for item in predicted
                if item.intensity >= strong_line_threshold * strongest_predicted
            )
            matched = {tuple(int(value) for value in item.miller_indices) for item in indexing}
            if strong:
                seen = [
                    item
                    for item in strong
                    if tuple(int(value) for value in item.miller_indices) in matched
                ]
                completeness = float(len(seen) / len(strong))
                missing = [
                    item.intensity / strongest_predicted
                    for item in strong
                    if tuple(int(value) for value in item.miller_indices) not in matched
                ]
                strongest_unobserved = float(max(missing)) if missing else 0.0

    if indexing.reflections:
        position_score = float(
            np.clip(
                1.0 - indexing.mean_absolute_delta_two_theta_deg / float(tolerance_deg),
                0.0,
                1.0,
            )
        )
    else:
        position_score = float("nan")

    if len(indexing.reflections) >= 2:
        intensity_agreement = _bray_curtis_similarity(
            np.array([item.peak.integrated_intensity for item in indexing.reflections]),
            np.array(
                [item.relative_intensity_calculated for item in indexing.reflections]
            ),
        )
    else:
        intensity_agreement = float("nan")

    return PhaseCandidateScore(
        phase_name=phase_name,
        source=source,
        indexing=indexing,
        explained_intensity_fraction=explained_fraction,
        completeness=completeness,
        position_score=position_score,
        intensity_agreement=intensity_agreement,
        strongest_unexplained_fraction=strongest_unexplained,
        strongest_unobserved_relative_intensity=strongest_unobserved,
        strong_line_threshold=float(strong_line_threshold),
        weights=weights,
    )


def identify_phase(
    table: PeakTable,
    candidates: Mapping[str, Phase] | Sequence[Phase] | Sequence[tuple[str, Phase]],
    *,
    radiation: RadiationSpec | None = None,
    sources: Mapping[str, str] | None = None,
    tolerance_deg: float = 0.3,
    max_index: int = 6,
    minimum_relative_intensity: float = 0.001,
    strong_line_threshold: float = 0.05,
    weights: Mapping[str, float] | None = None,
    minimum_score: float = 0.55,
    decisive_margin: float = 0.05,
    name: str | None = None,
) -> PhaseIdentification:
    """Rank candidate phases against a table of fitted peaks.

    Purpose
    -------
    Decide which of several candidate structures a measured powder pattern
    belongs to, and say how confident that decision is. This is the step
    between fitting the peaks of an unknown specimen and doing anything
    quantitative with it: a lattice parameter, a texture, a stress and a
    Rietveld refinement all require the phase to be settled first.

    When and where to use it
    ------------------------
    On a laboratory diffractogram of a specimen whose phase is suspected but
    not established -- an as-received alloy, a corrosion product, a powder from
    a synthesis that may or may not have gone as intended. Supply the plausible
    structures as CIF files (:meth:`pytex.core.lattice.Phase.from_cif`) or from
    the built-in catalogue, and read the ranked report.

    Do not use it as a database search. The candidates must be supplied; this
    ranks what it is given and cannot propose a phase nobody thought of.

    Method
    ------
    Each candidate is indexed against the measured peaks by
    :func:`~pytex.diffraction.xrd_indexing.index_peaks` -- a global one-to-one
    assignment, not a greedy nearest-line pass -- and then scored on the four
    criteria described in this module's documentation: explained intensity,
    completeness, position agreement and intensity agreement. The combined
    score is their weighted mean over whichever are defined, and the candidates
    are ranked by it.

    Two qualifications are then applied to the winner and reported separately
    from the ranking, because a ranking always has a winner and that is not the
    same as having an answer:

    - :attr:`PhaseIdentification.is_conclusive` -- the best score reaches
      ``minimum_score``, so the winner explains the pattern in absolute terms.
    - :attr:`PhaseIdentification.is_decisive` -- the winner leads the runner-up
      by at least ``decisive_margin``, so the ranking distinguishes them.

    Parameters
    ----------
    table
        Fitted peaks from the measured pattern. Use
        :meth:`~pytex.diffraction.xrd_peaks.PeakTable.filter_converged` first
        if any fit failed: a peak whose position did not converge carries a
        position that is not a measurement, and every candidate is penalized by
        it equally and wrongly.
    candidates
        The candidate phases: a mapping of name to phase, a sequence of phases
        (named from each phase's own ``name``), or a sequence of
        ``(name, phase)`` pairs. Duplicate names are numbered rather than
        merged, since two uploaded CIFs can legitimately share a formula.
    radiation
        Radiation to calculate positions with. Falls back to the table's, which
        is the normal case.
    sources
        Optional provenance per candidate name -- typically the CIF file name
        -- carried into the report so a ranking can be traced to its files.
    tolerance_deg
        Matching tolerance for indexing, in degrees ``2*theta``. It must exceed
        the instrument's uncorrected zero-point and displacement errors and
        stay below the spacing between neighbouring calculated lines. It also
        sets the scale of ``position_score``, so widening it to rescue a
        candidate lowers that candidate's score rather than raising it.
    max_index
        Largest ``|h|``, ``|k|``, ``|l|`` enumerated for every candidate.
    minimum_relative_intensity
        Calculated families weaker than this fraction of the strongest are not
        offered for matching at all.
    strong_line_threshold
        Relative intensity above which a predicted line counts towards
        ``completeness``. The default of 0.05 asks only about lines an operator
        would expect to see on the plot.
    weights
        Criterion weighting, overriding :data:`DEFAULT_CRITERION_WEIGHTS`.
        Missing criteria default to zero weight. Lower
        ``intensity_agreement`` further for a specimen known to be textured.
    minimum_score, decisive_margin
        Thresholds for the two qualifications above.
    name
        A name for the report.

    Returns
    -------
    PhaseIdentification
        The ranked candidates, the best match, and the two qualifications.

    Raises
    ------
    ValueError
        If the table is empty, no candidate is offered, the table carries no
        radiation and none was supplied, the tolerance is not positive, or a
        weighting is invalid.
    TypeError
        If a candidate is not a :class:`~pytex.core.lattice.Phase`.

    See Also
    --------
    identify_phase_from_pattern : the same, starting from a raw diffractogram.
    pytex.diffraction.xrd_indexing.index_peaks : the per-candidate assignment.
    pytex.diffraction.rietveld : the quantitative step once the phases are settled.

    Notes
    -----
    A high score is evidence, not proof. Two structures related by a small
    distortion, or the same framework with a different cation, can be
    indistinguishable on a laboratory scan of limited angular range; the
    result says so through ``is_decisive`` rather than leaving the reader to
    notice. Read the runner-up's criteria as well as the winner's.

    Examples
    --------
    >>> from pytex.diffraction.xrd_phase_identification import identify_phase
    >>> report = identify_phase(table, {"fcc Ni": nickel, "bcc Fe": iron})  # doctest: +SKIP
    >>> report.best.phase_name  # doctest: +SKIP
    'fcc Ni'
    """

    if len(table) == 0:
        raise ValueError("identify_phase was given an empty peak table.")
    if not np.isfinite(tolerance_deg) or tolerance_deg <= 0.0:
        raise ValueError("identify_phase requires a finite, positive tolerance_deg.")
    if not 0.0 < float(strong_line_threshold) <= 1.0:
        raise ValueError("strong_line_threshold must lie in (0, 1].")
    spec = radiation if radiation is not None else table.radiation
    if spec is None:
        raise ValueError(
            "identify_phase needs a radiation: the peak table declares none and none was "
            "passed."
        )

    resolved_weights = _normalized_weights(weights)
    pairs = _as_named_candidates(candidates)
    provenance = dict(sources or {})

    total_intensity = float(sum(peak.integrated_intensity for peak in table))
    strongest_intensity = float(max(peak.integrated_intensity for peak in table))

    scored = tuple(
        _score_candidate(
            phase_name=phase_name,
            phase=phase,
            source=provenance.get(phase_name, ""),
            table=table,
            radiation=spec,
            tolerance_deg=float(tolerance_deg),
            max_index=int(max_index),
            minimum_relative_intensity=float(minimum_relative_intensity),
            strong_line_threshold=float(strong_line_threshold),
            weights=resolved_weights,
            total_measured_intensity=total_intensity,
            strongest_measured_intensity=strongest_intensity,
        )
        for phase_name, phase in pairs
    )

    return PhaseIdentification(
        name=name or f"{table.name} phase identification",
        candidates=scored,
        minimum_score=float(minimum_score),
        decisive_margin=float(decisive_margin),
        peak_count=len(table),
        settings={
            "tolerance_deg": float(tolerance_deg),
            "max_index": float(max_index),
            "minimum_relative_intensity": float(minimum_relative_intensity),
            "strong_line_threshold": float(strong_line_threshold),
            "source_table": table.name,
            "radiation": spec.name,
        },
    )


def identify_phase_from_pattern(
    measured: MeasuredPowderPattern,
    candidates: Mapping[str, Phase] | Sequence[Phase] | Sequence[tuple[str, Phase]],
    *,
    radiation: RadiationSpec | None = None,
    sources: Mapping[str, str] | None = None,
    prominence_sigma: float = 5.0,
    expected_fwhm_deg: float | None = None,
    minimum_two_theta_deg: float | None = None,
    max_peaks: int = 128,
    tolerance_deg: float = 0.3,
    max_index: int = 6,
    minimum_relative_intensity: float = 0.001,
    strong_line_threshold: float = 0.05,
    weights: Mapping[str, float] | None = None,
    minimum_score: float = 0.55,
    decisive_margin: float = 0.05,
    name: str | None = None,
) -> tuple[PhaseIdentification, PeakTable]:
    """Detect and fit the peaks of a measured pattern, then rank candidate phases.

    Purpose
    -------
    Be the single call an operator makes on an unknown specimen: raw scan and a
    handful of candidate CIFs in, a ranked and qualified identification out.

    Method
    ------
    :func:`~pytex.diffraction.xrd_peaks.detect_and_fit_peaks` locates the
    reflections by a scale-matched Ricker filter on the variance-stabilized
    profile and fits each with a pseudo-Voigt, then :func:`identify_phase`
    ranks the candidates against the resulting table.

    The peak table is returned alongside the identification, not consumed
    inside it. Peak detection is the step where an identification most often
    goes wrong -- a threshold that swallows the weak lines, or one that
    promotes shoulders of the background into peaks -- and the table is where
    that is visible.

    Parameters
    ----------
    measured
        The measured diffractogram.
    candidates, radiation, sources, tolerance_deg, max_index,
    minimum_relative_intensity, strong_line_threshold, weights, minimum_score,
    decisive_margin, name
        As for :func:`identify_phase`.
    prominence_sigma
        Detection threshold in units of the local noise. Lower it to admit weak
        lines, at the cost of admitting background structure with them.
    expected_fwhm_deg
        Expected peak width, which sets the filter scale. Estimated from the
        pattern when omitted.
    minimum_two_theta_deg
        Ignore everything below this angle. Useful when the low-angle end
        carries a beam-stop or air-scatter rise that is not diffraction.
    max_peaks
        Largest number of peaks fitted.

    Returns
    -------
    tuple[PhaseIdentification, PeakTable]
        The ranked identification, and the fitted peaks it was formed from.

    Raises
    ------
    ValueError
        If no peak survives detection, which for a real diffractogram means the
        threshold is too high rather than that the specimen is amorphous, or
        for any reason :func:`identify_phase` raises.

    See Also
    --------
    identify_phase : the ranking step alone, on an existing peak table.
    pytex.diffraction.xrd_peaks.detect_and_fit_peaks : the detection step alone.
    """

    lower = (
        None
        if minimum_two_theta_deg is None
        else (float(minimum_two_theta_deg), float(measured.two_theta_deg[-1]))
    )
    table = detect_and_fit_peaks(
        measured,
        radiation=radiation,
        expected_fwhm_deg=expected_fwhm_deg,
        prominence_sigma=prominence_sigma,
        two_theta_range_deg=lower,
        max_peaks=max_peaks,
    )
    converged = table.filter_converged()
    if len(converged) > 0:
        table = converged
    identification = identify_phase(
        table,
        candidates,
        radiation=radiation,
        sources=sources,
        tolerance_deg=tolerance_deg,
        max_index=max_index,
        minimum_relative_intensity=minimum_relative_intensity,
        strong_line_threshold=strong_line_threshold,
        weights=weights,
        minimum_score=minimum_score,
        decisive_margin=decisive_margin,
        name=name,
    )
    return identification, table


__all__ = [
    "CRITERION_NAMES",
    "DEFAULT_CRITERION_WEIGHTS",
    "PHASE_CANDIDATE_SCORE_SCHEMA",
    "PHASE_IDENTIFICATION_SCHEMA",
    "PhaseCandidateScore",
    "PhaseIdentification",
    "identify_phase",
    "identify_phase_from_pattern",
]
