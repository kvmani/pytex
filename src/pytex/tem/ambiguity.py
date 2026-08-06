"""What an indexed diffraction pattern does *not* determine, stated exactly.

The "180 degree ambiguity" of SAED indexing is three unrelated things with
different causes, different consequences and different remedies. Conflating them
is why the failure mode has a reputation for being intractable, so this module
keeps them apart and names them:

**Layer 1 — crystallographic (Friedel/Laue).** Kinematic intensities obey
``|F(g)| = |F(-g)|``, so the recorded pattern is centrosymmetric even when the
crystal is not. The reconstructed orientation is therefore determined only up to
the rotations of the **Laue class** that map the zone plane to itself. The
decisive and often-missed consequence: for a centrosymmetric crystal the Laue
rotation group *equals* the crystal's own proper group, so every "ambiguous"
solution is a genuine symmetry equivalent and nothing is lost. Enumerating the
32 point groups, the Laue rotation group is strictly larger for exactly ten of
them — those with improper operations other than inversion. A chiral crystal such
as quartz (32) is non-centrosymmetric yet ambiguity-free, so "non-centrosymmetric"
is the wrong test and would warn on eleven groups without cause.

**Layer 2 — instrumental rotation.** An error in the diffraction rotation is a
rotation about the beam axis, is *not* absorbed by crystal symmetry, and produces
a residual ``2 asin(sin(dphi/2) sin theta)`` where ``theta`` is the angle between
current and target zones. At ``dphi = 180 deg`` this negates both alpha and beta:
the operator tilts exactly the wrong way while the calculation reports a clean
zero residual. This — not Friedel — is the ambiguity that wastes a session.

**Layer 3 — parity.** A mirrored pattern rendering or a reversed readout sign
reflects the trajectory rather than rotating it. Detected by the reconstructed
orientation coming out improper, or by the two-zone consistency test.

The engine must never silently pick one member of an ambiguity class. It reports
the count of genuinely distinct families, emits each as its own solution, and
states the experiment that discriminates them — including the predicted outcome
for each alternative, so the observation is a decisive test rather than a hint.

See ``docs/architecture/tem_tilt_navigation_foundation.md`` section 8.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.lattice import Phase
from pytex.core.point_groups import laue_class_symbol_for, normalize_point_group_symbol
from pytex.core.symmetry import SymmetrySpec

__all__ = [
    "AMBIGUOUS_POINT_GROUPS",
    "AmbiguityFamily",
    "AmbiguityLayer",
    "AmbiguityReport",
    "DiscriminatingExperiment",
    "analyze_ambiguity",
    "laue_rotation_operators",
    "observation_stabilizer",
]

#: The ten point groups for which Friedel's law adds rotations the crystal does
#: not have, so a single kinematic pattern leaves a genuine layer-1 ambiguity.
#:
#: These are exactly the non-centrosymmetric groups containing improper
#: operations other than inversion. The eleven centrosymmetric groups and the
#: eleven purely rotational (enantiomorphic) groups are unaffected: for those the
#: Laue rotation group equals the crystal's proper point group. Verified by
#: enumeration against ``pytex.core.point_groups`` in the TN test suite.
AMBIGUOUS_POINT_GROUPS: frozenset[str] = frozenset(
    {"m", "mm2", "-4", "4mm", "-42m", "3m", "-6", "6mm", "-6m2", "-43m"}
)

_OPERATOR_TOLERANCE = 1e-8


class AmbiguityLayer(StrEnum):
    """Which of the three independent ambiguity mechanisms is in play."""

    #: Friedel/Laue: intrinsic to a single kinematic pattern.
    CRYSTALLOGRAPHIC = "crystallographic"
    #: Diffraction-rotation error: rotation about the beam axis.
    INSTRUMENTAL_ROTATION = "instrumental_rotation"
    #: Mirrored rendering or reversed readout sign.
    PARITY = "parity"


@dataclass(frozen=True, slots=True)
class DiscriminatingExperiment:
    """An observation that tells two ambiguity families apart.

    Purpose
    -------
    Makes an ambiguity *actionable*. Reporting that two answers exist is only
    half an answer; the operator needs to know which observation separates them
    and what each alternative predicts, so that the observation is a decisive
    test rather than a hint.

    Attributes
    ----------
    name : str
        Short identifier, e.g. ``"tilt_excursion"``.
    procedure : str
        What the operator does, concretely enough to follow at the microscope.
    predicted_outcomes : tuple of str
        One prediction per family, in the same order as
        :attr:`AmbiguityReport.families`. A family whose prediction matches the
        observation is the true one.
    cost : str
        Rough experimental cost, so the operator can choose the cheapest
        decisive test.
    resolves : tuple of AmbiguityLayer
    """

    name: str
    procedure: str
    predicted_outcomes: tuple[str, ...]
    cost: str
    resolves: tuple[AmbiguityLayer, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "predicted_outcomes", tuple(self.predicted_outcomes))
        object.__setattr__(self, "resolves", tuple(self.resolves))
        if not self.predicted_outcomes:
            raise ValueError("A discriminating experiment needs at least one outcome.")

    def describe(self) -> str:
        outcomes = "; ".join(
            f"family {index + 1}: {outcome}"
            for index, outcome in enumerate(self.predicted_outcomes)
        )
        return f"{self.procedure} ({self.cost}). Predicted — {outcomes}."


@dataclass(frozen=True, slots=True)
class AmbiguityFamily:
    """One competing hypothesis about the true crystal-to-holder orientation.

    Purpose
    -------
    Distinguishes a *hypothesis* from a *choice*. A symmetry-equivalent target is
    a choice the operator may make freely — every equivalent gives the same
    pattern. A family is a statement about reality of which only one member is
    true, and acting on the wrong one sends the specimen somewhere else. The two
    are generated by similar group operations and must never be presented alike.

    Attributes
    ----------
    index : int
        1-based, matching the order in :attr:`AmbiguityReport.families`.
    operator : np.ndarray
        The 3x3 crystal-frame rotation relating this family to family 1. The
        identity for the first family.
    is_symmetry_equivalent : bool
        Whether ``operator`` is a proper crystal symmetry operation. When true
        the family is not really a competing hypothesis at all and is reported
        for completeness only.
    layer : AmbiguityLayer
    rationale : str
        Why this family exists, in one sentence.
    """

    index: int
    operator: np.ndarray
    is_symmetry_equivalent: bool
    layer: AmbiguityLayer
    rationale: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "operator", as_float_array(self.operator, shape=(3, 3)))
        if self.index < 1:
            raise ValueError("AmbiguityFamily.index is 1-based.")

    @property
    def rotation_angle_deg(self) -> float:
        """Angle of :attr:`operator`, in degrees; zero for the identity family."""

        trace = float(np.trace(self.operator))
        return float(math.degrees(math.acos(max(-1.0, min(1.0, (trace - 1.0) / 2.0)))))

    def describe(self) -> str:
        if self.index == 1:
            return "Family 1: the reference reconstruction."
        status = (
            "a genuine crystal symmetry, so this family is the same physical "
            "orientation described differently"
            if self.is_symmetry_equivalent
            else "NOT a crystal symmetry, so this is a physically different orientation"
        )
        return (
            f"Family {self.index}: related to family 1 by a "
            f"{self.rotation_angle_deg:.1f} deg rotation, which is {status}. "
            f"{self.rationale}"
        )


@dataclass(frozen=True, slots=True)
class AmbiguityReport:
    """What the available observations leave undetermined, and how to fix it.

    Purpose
    -------
    The honest accounting demanded by the foundation document: how many
    genuinely distinct orientation hypotheses the data admit, which mechanism
    produced them, and what to measure to collapse them. When the answer is
    "none", it says so plainly — warning fatigue is itself a failure mode, and a
    report that always warns is a report nobody reads.

    Attributes
    ----------
    families : tuple of AmbiguityFamily
        At least one. Length greater than one means the data do not determine
        the orientation uniquely.
    stabilizer_order : int
        Order of the observation stabilizer: the rotations of the Laue class
        that map the zone plane to itself.
    symmetry_stabilizer_order : int
        Order of its intersection with the crystal's proper point group.
    layers : tuple of AmbiguityLayer
        Which mechanisms contribute.
    experiments : tuple of DiscriminatingExperiment
    point_group : str
    zone_axis_sense_determined : bool
        Whether ``[uvw]`` can be told from its reverse.
    notes : tuple of str
    """

    families: tuple[AmbiguityFamily, ...]
    stabilizer_order: int
    symmetry_stabilizer_order: int
    layers: tuple[AmbiguityLayer, ...]
    experiments: tuple[DiscriminatingExperiment, ...]
    point_group: str
    zone_axis_sense_determined: bool
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "families", tuple(self.families))
        object.__setattr__(self, "layers", tuple(self.layers))
        object.__setattr__(self, "experiments", tuple(self.experiments))
        object.__setattr__(self, "notes", tuple(self.notes))
        if not self.families:
            raise ValueError("An AmbiguityReport must contain at least one family.")

    @property
    def is_unique(self) -> bool:
        """Whether the observations determine one physical orientation."""

        return len(self.families) == 1

    @property
    def distinct_family_count(self) -> int:
        """Number of physically distinct orientation hypotheses."""

        return sum(1 for family in self.families if not family.is_symmetry_equivalent) or 1

    def describe(self) -> str:
        """Convention-explicit prose: what is determined, what is not, what to do."""

        if self.is_unique:
            head = (
                f"The observations determine the crystal-to-holder orientation "
                f"uniquely for point group {self.point_group}. The Laue rotation group "
                f"adds nothing the crystal does not already have, so every alternative "
                f"indexing of this pattern is a symmetry-equivalent description of the "
                f"same physical orientation (observation stabilizer of order "
                f"{self.stabilizer_order}, entirely contained in the proper point group)."
            )
            sense = (
                "The sense of the zone axis is determined."
                if self.zone_axis_sense_determined
                else (
                    "The sense of the zone axis is not determined — a two-fold "
                    "perpendicular to it maps [uvw] to its reverse — but that operator is "
                    "a crystal symmetry, so both descriptions denote the same orientation."
                )
            )
            return f"{head} {sense}"

        head = (
            f"The observations leave {len(self.families)} distinct orientation "
            f"hypotheses for point group {self.point_group}. Layers involved: "
            f"{', '.join(layer.value for layer in self.layers)}. The observation "
            f"stabilizer has order {self.stabilizer_order}, of which "
            f"{self.symmetry_stabilizer_order} operators are crystal symmetries; the "
            f"quotient is what remains undetermined."
        )
        families = " ".join(family.describe() for family in self.families)
        experiments = " ".join(
            f"To discriminate: {experiment.describe()}" for experiment in self.experiments
        )
        warning = (
            "Each family is emitted as its own ranked solution with its own tilts; "
            "none is presented as the answer."
        )
        return f"{head} {families} {experiments} {warning}"

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "point_group": self.point_group,
            "is_unique": self.is_unique,
            "family_count": len(self.families),
            "distinct_family_count": self.distinct_family_count,
            "stabilizer_order": self.stabilizer_order,
            "symmetry_stabilizer_order": self.symmetry_stabilizer_order,
            "layers": [layer.value for layer in self.layers],
            "zone_axis_sense_determined": self.zone_axis_sense_determined,
            "families": [
                {
                    "index": family.index,
                    "operator": family.operator.tolist(),
                    "is_symmetry_equivalent": family.is_symmetry_equivalent,
                    "rotation_angle_deg": family.rotation_angle_deg,
                    "layer": family.layer.value,
                    "rationale": family.rationale,
                }
                for family in self.families
            ],
            "discriminating_experiments": [
                {
                    "name": experiment.name,
                    "procedure": experiment.procedure,
                    "predicted_outcomes": list(experiment.predicted_outcomes),
                    "cost": experiment.cost,
                    "resolves": [layer.value for layer in experiment.resolves],
                }
                for experiment in self.experiments
            ],
            "notes": list(self.notes),
        }


def laue_rotation_operators(point_group: str) -> np.ndarray:
    """Proper rotations of the Laue class of ``point_group``.

    Purpose
    -------
    The group that governs what a *single* kinematic diffraction pattern can
    determine. Friedel's law symmetrizes the observation to the Laue class, so
    two orientations differing by one of these rotations produce the same spot
    pattern.

    This is deliberately **not** the crystal's own proper point group. For ten of
    the thirty-two point groups the Laue rotation group is strictly larger — by a
    factor of exactly two — and the difference is precisely the extra ambiguity
    Friedel's law introduces.

    Parameters
    ----------
    point_group : str
        Hermann-Mauguin symbol.

    Returns
    -------
    np.ndarray
        Shape ``(n, 3, 3)`` rotation matrices in the crystal Cartesian frame.

    Examples
    --------
    A centrosymmetric group gains nothing::

        >>> len(laue_rotation_operators("m-3m"))
        24

    A group with improper operations other than inversion gains a factor of two::

        >>> len(laue_rotation_operators("-43m"))
        24
    """

    symbol = normalize_point_group_symbol(point_group)
    laue_symbol = laue_class_symbol_for(symbol)
    return np.asarray(
        SymmetrySpec.from_point_group(laue_symbol).operators, dtype=np.float64
    )


def observation_stabilizer(
    point_group: str, zone_axis_cartesian: ArrayLike
) -> np.ndarray:
    """Rotations of the Laue class that map the observed zone plane to itself.

    Purpose
    -------
    The exact set of orientations a single indexed pattern cannot tell apart. An
    operator qualifies when it maps the zone axis to plus or minus itself: the
    plus case is an in-plane rotation of the pattern, the minus case reverses the
    zone axis while preserving the observed spot positions, which is the classic
    "flip" ambiguity.

    Parameters
    ----------
    point_group : str
    zone_axis_cartesian : array_like
        The zone axis as a Cartesian crystal-frame vector; need not be
        normalized.

    Returns
    -------
    np.ndarray
        Shape ``(n, 3, 3)``. Always contains the identity.
    """

    axis = normalize_vector(zone_axis_cartesian)
    operators = laue_rotation_operators(point_group)
    images = np.einsum("nij,j->ni", operators, axis)
    keeps = np.all(np.abs(images - axis) < _OPERATOR_TOLERANCE, axis=1)
    reverses = np.all(np.abs(images + axis) < _OPERATOR_TOLERANCE, axis=1)
    return np.ascontiguousarray(operators[keeps | reverses])


def _contains_operator(group: np.ndarray, operator: np.ndarray) -> bool:
    if group.size == 0:
        return False
    return bool(np.any(np.all(np.abs(group - operator) < _OPERATOR_TOLERANCE, axis=(1, 2))))


def _coset_representatives(
    stabilizer: np.ndarray, proper_group: np.ndarray
) -> list[np.ndarray]:
    """Representatives of ``stabilizer / (stabilizer & proper_group)``.

    The quotient counts the genuinely distinct orientation hypotheses: operators
    inside the proper group merely relabel the same physical orientation, while
    each remaining coset is a different one.
    """

    intersection = np.asarray(
        [op for op in stabilizer if _contains_operator(proper_group, op)],
        dtype=np.float64,
    )
    if intersection.size == 0:
        intersection = np.eye(3, dtype=np.float64)[None, :, :]
    representatives: list[np.ndarray] = [np.eye(3, dtype=np.float64)]
    for candidate in stabilizer:
        covered = False
        for representative in representatives:
            # candidate lies in representative * intersection?
            relative = representative.T @ candidate
            if _contains_operator(intersection, relative):
                covered = True
                break
        if not covered:
            representatives.append(np.asarray(candidate, dtype=np.float64))
    return representatives


def _sense_is_determined(stabilizer: np.ndarray, axis: np.ndarray) -> bool:
    """Whether no stabilizer operator reverses the zone axis."""

    images = np.einsum("nij,j->ni", stabilizer, axis)
    return not bool(np.any(np.all(np.abs(images + axis) < _OPERATOR_TOLERANCE, axis=1)))


def _standard_experiments(family_count: int) -> tuple[DiscriminatingExperiment, ...]:
    """The discriminating observations, ordered cheapest-decisive first."""

    def outcomes(template: str, opposite: str) -> tuple[str, ...]:
        values = [template, opposite]
        while len(values) < family_count:
            values.append("a distinct, computable pattern motion")
        return tuple(values[:family_count])

    return (
        DiscriminatingExperiment(
            name="tilt_excursion",
            procedure=(
                "Apply a known small positive alpha tilt (5-10 deg) and record how the "
                "Kikuchi pattern moves"
            ),
            predicted_outcomes=outcomes(
                "the pattern translates along the predicted azimuth",
                "the pattern translates along the opposite azimuth",
            ),
            cost="two exposures, about one minute",
            resolves=(
                AmbiguityLayer.INSTRUMENTAL_ROTATION,
                AmbiguityLayer.PARITY,
            ),
        ),
        DiscriminatingExperiment(
            name="second_zone_axis",
            procedure=(
                "Index a second zone axis at a different stage position and rebuild the "
                "orientation from the two-zone path, which needs no rotation calibration"
            ),
            predicted_outcomes=outcomes(
                "the interzonal angle matches the crystallographic value",
                "the interzonal angle disagrees, indicting this family",
            ),
            cost="free when a second zone was already visited",
            resolves=(
                AmbiguityLayer.INSTRUMENTAL_ROTATION,
                AmbiguityLayer.PARITY,
            ),
        ),
        DiscriminatingExperiment(
            name="cbed_holz_symmetry",
            procedure=(
                "Record a convergent-beam pattern and read the HOLZ ring symmetry, which "
                "is sensitive to the polarity Friedel's law hides"
            ),
            predicted_outcomes=outcomes(
                "the HOLZ symmetry matches this polarity",
                "the HOLZ symmetry matches the reversed polarity",
            ),
            cost="convergent probe on a thin, clean area",
            resolves=(AmbiguityLayer.CRYSTALLOGRAPHIC,),
        ),
    )


def analyze_ambiguity(
    phase: Phase,
    zone_axis_cartesian: ArrayLike,
    *,
    rotation_calibrated: bool = True,
    reconstruction_note: str = "",
) -> AmbiguityReport:
    """Classify what the current observations leave undetermined.

    Purpose
    -------
    Answers, for a specific phase and a specific indexed zone axis, the question
    the foundation document insists must never be answered silently: *is this
    orientation actually determined, and if not, in how many ways, and what
    should I measure?*

    When to use
    -----------
    Called automatically by :func:`pytex.tem.navigation.plan_tilt_to_zone_axis`.
    Call it directly when auditing a reconstruction before trusting it, or when
    teaching why a particular crystal class is or is not affected.

    Parameters
    ----------
    phase : Phase
        Supplies the point group. Its Laue class decides layer 1.
    zone_axis_cartesian : array_like
        The indexed zone axis, as a Cartesian crystal-frame vector.
    rotation_calibrated : bool, default True
        Whether a diffraction rotation is available. ``False`` adds a layer-2
        family, because without the calibration the pattern azimuth — and hence
        the tilt direction — is unknown.
    reconstruction_note : str, optional
        Carried into the report notes; used to record which reconstruction mode
        produced the estimate.

    Returns
    -------
    AmbiguityReport

    Examples
    --------
    A centrosymmetric cubic crystal down ``[001]`` is unambiguous, and the
    report says so rather than warning. The observation stabilizer has order 8 —
    the four-fold about ``[001]`` and the two-folds about the in-plane
    ``<100>`` and ``<110>`` axes, that is the group 422 — and every one of those
    operators is a crystal symmetry of ``m-3m``, so nothing is left undetermined.
    """

    axis = normalize_vector(zone_axis_cartesian)
    point_group = normalize_point_group_symbol(phase.symmetry.point_group)
    stabilizer = observation_stabilizer(point_group, axis)
    proper_operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
    symmetry_stabilizer = np.asarray(
        [op for op in stabilizer if _contains_operator(proper_operators, op)],
        dtype=np.float64,
    )
    representatives = _coset_representatives(stabilizer, proper_operators)

    layers: list[AmbiguityLayer] = []
    families: list[AmbiguityFamily] = [
        AmbiguityFamily(
            index=1,
            operator=np.eye(3, dtype=np.float64),
            is_symmetry_equivalent=True,
            layer=AmbiguityLayer.CRYSTALLOGRAPHIC,
            rationale="Reference reconstruction.",
        )
    ]
    for representative in representatives[1:]:
        layers.append(AmbiguityLayer.CRYSTALLOGRAPHIC)
        families.append(
            AmbiguityFamily(
                index=len(families) + 1,
                operator=representative,
                is_symmetry_equivalent=False,
                layer=AmbiguityLayer.CRYSTALLOGRAPHIC,
                rationale=(
                    "Friedel's law makes the recorded pattern centrosymmetric, so this "
                    f"Laue-class rotation is indistinguishable from the identity for "
                    f"point group {point_group}, which lacks it. The destination is a "
                    "Friedel partner: the diffraction pattern there is identical, but "
                    "the physical orientation and any polarity-sensitive measurement "
                    "differ."
                ),
            )
        )

    if not rotation_calibrated:
        layers.append(AmbiguityLayer.INSTRUMENTAL_ROTATION)
        families.append(
            AmbiguityFamily(
                index=len(families) + 1,
                operator=np.diag([-1.0, -1.0, 1.0]).astype(np.float64),
                is_symmetry_equivalent=_contains_operator(
                    proper_operators, np.diag([-1.0, -1.0, 1.0]).astype(np.float64)
                ),
                layer=AmbiguityLayer.INSTRUMENTAL_ROTATION,
                rationale=(
                    "The diffraction rotation is not calibrated, so the pattern azimuth "
                    "is known only up to a 180 degree turn about the beam. That error is "
                    "not absorbed by crystal symmetry: it negates both alpha and beta, "
                    "so the specimen tilts exactly the wrong way while the calculation "
                    "still reports a zero residual."
                ),
            )
        )

    notes: list[str] = []
    if reconstruction_note:
        notes.append(reconstruction_note)
    if point_group in AMBIGUOUS_POINT_GROUPS:
        notes.append(
            f"Point group {point_group} is one of the ten for which Friedel's law adds "
            "rotations the crystal does not have; a layer-1 ambiguity is possible for "
            "some zones."
        )
    else:
        notes.append(
            f"Point group {point_group} is centrosymmetric or purely rotational, so "
            "Friedel's law adds no rotation the crystal lacks and layer 1 contributes "
            "no ambiguity."
        )

    return AmbiguityReport(
        families=tuple(families),
        stabilizer_order=int(stabilizer.shape[0]),
        symmetry_stabilizer_order=int(
            symmetry_stabilizer.shape[0] if symmetry_stabilizer.size else 0
        ),
        layers=tuple(dict.fromkeys(layers)),
        experiments=_standard_experiments(len(families)) if len(families) > 1 else (),
        point_group=point_group,
        zone_axis_sense_determined=_sense_is_determined(stabilizer, axis),
        notes=tuple(notes),
    )
