"""TEM specimen-tilt navigation: reaching a target zone axis on a real holder.

The problem this subpackage solves is the one asked in front of a microscope.
An operator has reached some zone axis by chasing Kikuchi bands, the pattern is
indexed, the stage reads some ``(alpha, beta)`` — and now a *different* zone axis
is wanted. Which way to tilt, how far, is it even reachable, and which of the
symmetry-equivalent versions of the target will actually be found there?

`pytex.tem` answers that end to end: it reconstructs the crystal-to-holder
orientation from indexed diffraction data, enumerates every symmetry-equivalent
and ambiguity-related version of the requested target, solves the holder angles
in closed form, validates each solution by forward simulation through the
calibrated stage model, ranks the survivors, and plans a path to the best one.

Why a separate subpackage
-------------------------
This is a new *domain of use* — instrument operation — rather than new
crystallography. It composes `pytex.core` symmetry and orientation primitives
with `pytex.diffraction` pattern solving; nothing in the repository depends on
it in turn.

The scientific formulation, including the derivations, the three-layer ambiguity
analysis and the validation strategy, is
``docs/architecture/tem_tilt_navigation_foundation.md``. Read that before
changing conventions here.

Where to start
--------------
`plan_tilt_to_zone_axis` is the single entry point most callers need; build its
``current`` argument with one of the :class:`CurrentState` constructors, whose
choice determines what can and cannot be resolved:

- :meth:`CurrentState.from_two_zone_axes` — **recommended**. Two indexed zones at
  two stage positions determine the orientation with no pattern-rotation
  calibration at all.
- :meth:`CurrentState.from_pattern_solution` — one indexed pattern plus a
  calibrated diffraction rotation. Convenient, but it carries the instrumental
  ambiguity that causes an operator to tilt in exactly the wrong direction.
- :meth:`CurrentState.from_orientation` — a known crystal-to-holder orientation,
  for synthetic work and teaching.
"""

from __future__ import annotations

from pytex.tem.ambiguity import (
    AmbiguityFamily,
    AmbiguityLayer,
    AmbiguityReport,
    DiscriminatingExperiment,
    analyze_ambiguity,
    observation_stabilizer,
)
from pytex.tem.calibration import (
    TiltExcursionObservation,
    calibrate_from_tilt_excursions,
    fit_stage_and_orientation,
)
from pytex.tem.navigation import (
    DEFAULT_RANKING,
    RankingWeights,
    Reachability,
    TiltPlanReport,
    TiltSolution,
    plan_tilt_to_zone_axis,
    solve_tilts_for_direction,
)
from pytex.tem.path import (
    PathStrategy,
    TiltPath,
    TiltPathSample,
    connecting_band,
    plan_path,
)
from pytex.tem.reconstruction import (
    CurrentState,
    ReconstructionMode,
    ZoneAxisObservation,
)
from pytex.tem.stage import (
    BEAM_AXIS_LABORATORY,
    DoubleTiltStage,
    EllipticalEnvelope,
    GeneralStageAxes,
    HolderKind,
    MaskedEnvelope,
    PolygonEnvelope,
    RectangularEnvelope,
    SingleTiltStage,
    StageCalibration,
    StageModel,
    StagePosition,
    TiltEnvelope,
    TiltRotateStage,
    beam_direction_holder,
)

__all__ = [
    "BEAM_AXIS_LABORATORY",
    "DEFAULT_RANKING",
    "AmbiguityFamily",
    "AmbiguityLayer",
    "AmbiguityReport",
    "CurrentState",
    "DiscriminatingExperiment",
    "DoubleTiltStage",
    "EllipticalEnvelope",
    "GeneralStageAxes",
    "HolderKind",
    "MaskedEnvelope",
    "PathStrategy",
    "PolygonEnvelope",
    "RankingWeights",
    "Reachability",
    "ReconstructionMode",
    "RectangularEnvelope",
    "SingleTiltStage",
    "StageCalibration",
    "StageModel",
    "StagePosition",
    "TiltEnvelope",
    "TiltExcursionObservation",
    "TiltPath",
    "TiltPathSample",
    "TiltPlanReport",
    "TiltRotateStage",
    "TiltSolution",
    "ZoneAxisObservation",
    "analyze_ambiguity",
    "beam_direction_holder",
    "calibrate_from_tilt_excursions",
    "connecting_band",
    "fit_stage_and_orientation",
    "observation_stabilizer",
    "plan_path",
    "plan_tilt_to_zone_axis",
    "solve_tilts_for_direction",
]
