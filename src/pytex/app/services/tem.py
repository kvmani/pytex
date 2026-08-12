"""The TEM pattern solver: from a picked pattern to the next zone axis.

The workflow this implements is the one a microscopist actually performs, in
order:

1. **Calibrate.** Say what the picked coordinates mean — pixels, millimetres, or
   already-reciprocal — and supply the camera constant that scales them.
2. **Pick.** Mark the transmitted beam and the spots. This happens in the
   browser; what arrives here is a list of coordinates.
3. **Index.** Decide which phase, which zone axis, and which reflection each
   spot is, with residuals and the alternatives that were rejected.
4. **Navigate.** From the orientation that indexing established, work out the
   alpha and beta tilts that reach the next zone axis, which axes lie on the
   way, and whether the stage can get there at all.

Steps 3 and 4 are the ones worth automating. Indexing by hand means measuring
ratios and angles against tables; tilt planning by hand means composing two
non-orthogonal holder rotations in your head and then discovering, ten minutes
in, that the target was outside the beta range all along. Both are exactly the
kind of thing a computer should do while the microscope time is running.

Everything here is a thin surface over :mod:`pytex.diffraction.solving` and
:mod:`pytex.tem`; no crystallography is implemented in this module.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import direction_label, phase_parameter, plane_label

__all__ = ["measured_pattern_from_picks"]

_CITATION_WILLIAMS = "Williams & Carter, Transmission Electron Microscopy, 2nd ed., chapters 16-18."
_CITATION_EDINGTON = "Edington, Practical Electron Microscopy in Materials Science (1975)."


_UNIT_OPTIONS = (
    (
        "px",
        "Pixels",
        "Coordinates read off the image. Needs the pixel size as well as the camera constant.",
    ),
    ("mm", "Millimetres", "Coordinates measured on the detector or on a printed plate."),
    (
        "reciprocal_angstrom",
        "Å⁻¹ (already calibrated)",
        "Coordinates already converted to reciprocal space; no camera constant is used.",
    ),
)


def _picks(payload: Mapping[str, Any] | None) -> tuple[tuple[float, float], list[dict[str, Any]]]:
    """Split a picker payload into the beam centre and the spots."""

    if not isinstance(payload, Mapping):
        raise InvalidInputError(
            "Pick the transmitted beam and at least two spots.",
            field="picks",
            hint="Click the direct beam first, then two or more reflections.",
        )
    centre_raw = payload.get("centre")
    if not isinstance(centre_raw, Sequence) or len(centre_raw) != 2:
        raise InvalidInputError(
            "The transmitted beam has not been marked.",
            field="picks",
            hint=(
                "The direct beam is the origin every spot is measured from; it is not itself a "
                "reflection."
            ),
        )
    spots_raw = payload.get("spots") or []
    if not isinstance(spots_raw, Sequence) or len(spots_raw) < 2:
        raise InvalidInputError(
            "At least two spots are needed to index a pattern.",
            field="picks",
            hint=(
                "Two non-collinear reflections are the minimum that fixes a zone axis. Pick the "
                "two shortest ones you trust."
            ),
        )
    spots: list[dict[str, Any]] = []
    for index, entry in enumerate(spots_raw):
        if not isinstance(entry, Mapping):
            raise InvalidInputError(f"Spot {index + 1} is malformed.", field="picks")
        spots.append(
            {
                "x": float(entry["x"]),
                "y": float(entry["y"]),
                "intensity": (
                    None if entry.get("intensity") is None else float(entry["intensity"])
                ),
                "label": entry.get("label"),
            }
        )
    return (float(centre_raw[0]), float(centre_raw[1])), spots


def measured_pattern_from_picks(
    request: Mapping[str, Any],
) -> tuple[Any, tuple[float, float], list[dict[str, Any]]]:
    """Build a :class:`MeasuredSAEDPattern` from a picker payload and calibration.

    Kept as a named function rather than inlined because the solver and any
    future re-solve, refine, or export path must all interpret a set of picks the
    same way; two readings of one calibration is exactly how a pattern comes to
    be indexed against a camera constant nobody set.
    """

    from pytex.diffraction.solving import MeasuredSAEDPattern, MeasuredSpot, PatternCalibration

    centre, spots = _picks(request.get("picks"))
    units = str(request["units"])
    camera_constant = request.get("camera_constant_mm_angstrom")
    pixel_size = request.get("pixel_size_mm")

    if units == "px" and not pixel_size:
        raise InvalidInputError(
            "Pixel coordinates need a pixel size to become millimetres.",
            field="pixel_size_mm",
            hint="Read it from the detector specification, or calibrate against a known ring.",
        )
    if units in {"px", "mm"} and not camera_constant:
        raise InvalidInputError(
            "A camera constant is needed to turn detector distance into 1/d.",
            field="camera_constant_mm_angstrom",
            hint=(
                "Use the microscope's calibrated value for this camera length, or derive one "
                "from a spot whose d-spacing you already know."
            ),
        )

    calibration = PatternCalibration(
        units=units,
        centre=centre,
        camera_constant_mm_angstrom=(
            float(camera_constant) if units != "reciprocal_angstrom" and camera_constant else None
        ),
        pixel_size_mm=float(pixel_size) if units == "px" and pixel_size else None,
    )
    pattern = MeasuredSAEDPattern(
        name=str(request.get("pattern_name") or "picked pattern"),
        spots=tuple(
            MeasuredSpot(
                position=(spot["x"], spot["y"]),
                intensity=spot["intensity"],
                label=spot["label"],
            )
            for spot in spots
        ),
        calibration=calibration,
    )
    return pattern, centre, spots


_CALIBRATION_PARAMETERS = (
    ChoiceParameter(
        name="units",
        label="Coordinate units",
        help_text="What the picked coordinates are measured in.",
        options=_UNIT_OPTIONS,
        default="px",
        group="Calibration",
    ),
    NumberParameter(
        name="camera_constant_mm_angstrom",
        label="Camera constant",
        help_text=(
            "The instrument constant L·λ relating a spot's distance from the beam to 1/d: "
            "r = (camera constant) / d. This is the single number that turns a picture into a "
            "measurement, and using a value from the wrong camera length is the most common way "
            "to index a pattern to the wrong phase."
        ),
        units="mm·Å",
        default=180.0,
        minimum=0.0,
        required=False,
        group="Calibration",
    ),
    NumberParameter(
        name="pixel_size_mm",
        label="Pixel size",
        help_text="Detector pixel pitch, needed only when the coordinates are in pixels.",
        units="mm",
        default=0.05,
        minimum=0.0,
        required=False,
        group="Calibration",
    ),
)


@REGISTRY.operation(
    "tem.solve_pattern",
    title="Index a diffraction pattern",
    summary="From picked spots to phase, zone axis, and an index for every spot.",
    help_text=(
        "Takes the spots you picked and works out which phase they belong to, which zone axis the "
        "beam is down, and the indices of every reflection — with the residual of each, the "
        "fraction of spots that were indexed, and the alternative solutions that were "
        "considered.\n\n"
        "**Pick the beam first.** The transmitted beam is not a reflection; it is the origin "
        "every spot is measured from, and an error in its position biases every d-spacing in the "
        "pattern.\n\n"
        "**The method is geometric.** Indexing uses the two shortest non-collinear spots as a "
        "seed, requires both their lengths to match a calculated pair within the length "
        "tolerance and their angle within the angle tolerance, and then projects every remaining "
        "spot. Intensities are never used, which is deliberate: kinematic intensities are "
        "unreliable in a real pattern and geometry alone is enough.\n\n"
        "**Read the verdict, not just the answer.** A pattern with two solutions of similar score "
        "is genuinely ambiguous, and the result says so. That usually means the zone axis is a "
        "high-symmetry one where two phases give the same spot arrangement, and the way out is a "
        "second pattern at a different tilt rather than a tighter tolerance."
    ),
    parameters=(
        phase_parameter(help_text="The phase to test the pattern against."),
        ObjectParameter(
            name="picks",
            label="Picked spots",
            help_text=(
                "The transmitted beam and the reflections, marked on the image. Click the direct "
                "beam first, then the spots."
            ),
            editor="spot-picker",
        ),
        *_CALIBRATION_PARAMETERS,
        NumberParameter(
            name="length_tolerance",
            label="Length tolerance",
            help_text=(
                "Fractional agreement required between a measured |g| and a calculated one. "
                "0.03 means 3 percent, which is realistic for a well-calibrated instrument; "
                "loosening it admits more solutions rather than better ones."
            ),
            default=0.03,
            minimum=0.001,
            maximum=0.3,
            advanced=True,
        ),
        NumberParameter(
            name="angle_tolerance_deg",
            label="Angle tolerance",
            help_text="Agreement required between the measured and calculated interspot angle.",
            units="°",
            default=2.0,
            minimum=0.1,
            maximum=15.0,
            advanced=True,
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest |h|, |k| or |l| considered when matching.",
            default=4,
            minimum=1,
            maximum=8,
            advanced=True,
        ),
        ObjectParameter(
            name="second_phase",
            label="Second candidate phase",
            help_text=(
                "Optional. Give a second phase to have the solver choose between them — which is "
                "how a pattern is used to identify a precipitate rather than to confirm one."
            ),
            editor="phase",
            required=False,
        ),
    ),
    returns="One row per spot with its index and residual; the ranked solutions under `data`.",
    panel="tem",
    citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON),
    tags=("TEM", "SAED", "indexing", "solve", "zone axis", "pattern", "calibration"),
)
def _solve_pattern(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.diffraction.solving import solve_saed_pattern

    spec, phase = phase_from_request(request["phase"])
    phases = [phase]
    specs = [spec]
    if request.get("second_phase"):
        second_spec, second_phase = phase_from_request(request["second_phase"])
        phases.append(second_phase)
        specs.append(second_spec)

    pattern, centre, picks = measured_pattern_from_picks(request)
    report = solve_saed_pattern(
        pattern,
        phases,
        max_index=int(request["max_index"]),
        length_tolerance_relative=float(request["length_tolerance"]),
        angle_tolerance_deg=float(request["angle_tolerance_deg"]),
    )
    if not report.solutions:
        observed = pattern.d_spacings_angstrom()
        raise InvalidInputError(
            "No solution matches these spots within the tolerances.",
            field="picks",
            hint=(
                "Observed spacings are "
                + ", ".join(f"{value:.3f}" for value in sorted(observed, reverse=True)[:6])
                + " Å. If those look wrong for this phase the camera constant is the first thing "
                "to check; if they look right, loosen the length tolerance or raise the index "
                "limit."
            ),
        )
    best = report.best()

    spec_by_name = {item.name: item for item in specs}
    best_spec = spec_by_name.get(best.phase_name, spec)
    observed_d = np.asarray(pattern.d_spacings_angstrom(), dtype=float)
    rows: list[dict[str, Any]] = []
    for spot in best.solved_spots:
        index = int(spot.measured_index)
        indices = tuple(int(value) for value in np.asarray(spot.hkl, dtype=int))
        # `predicted_g_inv_angstrom` is the in-plane vector, not its length.
        predicted = float(np.linalg.norm(np.asarray(spot.predicted_g_inv_angstrom, dtype=float)))
        rows.append(
            {
                "spot": index + 1,
                "hkl": plane_label(indices, spec=best_spec) if any(indices) else "—",
                "d_observed": float(observed_d[index]),
                "d_calculated": 1.0 / predicted if predicted > 0.0 else float("nan"),
                "residual": float(spot.residual_inv_angstrom),
                "x": picks[index]["x"],
                "y": picks[index]["y"],
            }
        )
    rows.sort(key=lambda row: int(row["spot"]))

    zone_indices = tuple(int(value) for value in np.asarray(best.zone_axis.indices, dtype=int))
    zone_text = direction_label(zone_indices, spec=best_spec)
    alternatives = [
        {
            "phase": solution.phase_name,
            "zone_axis": solution.zone_axis_label,
            # `score` is a sort key (matched fraction, then residual), not a
            # scalar quality: reporting it as a number would invite a comparison
            # it does not support.
            "matched_spots": round(solution.matched_fraction * solution.measured_spot_count),
            "matched_fraction": float(solution.matched_fraction),
            "mean_residual_inv_angstrom": float(solution.mean_residual_inv_angstrom),
        }
        for solution in report.solutions
    ]
    conclusive = bool(report.is_conclusive)
    notes: list[str] = []
    if not conclusive:
        notes.append(
            "Two or more solutions score similarly, so this pattern does not identify the phase "
            "or the axis on its own. Record a second pattern at a different tilt: two zone axes "
            "settle what one cannot."
        )
    if best.unindexed_spot_indices:
        notes.append(
            f"{len(best.unindexed_spot_indices)} spot(s) were not indexed. Double diffraction, a "
            "second phase, or a misplaced pick are the usual causes, in that order."
        )

    result = AppResult(
        title=f"{best.phase_name} down {zone_text}",
        summary=(
            f"The pattern indexes as {best.phase_name} with the beam along {zone_text}. "
            f"{len(rows)} of {best.measured_spot_count} picked spots were indexed "
            f"({best.matched_fraction * 100:.0f}%), with a mean residual of "
            f"{best.mean_residual_inv_angstrom:.4f} Å⁻¹. "
            + (
                "No competing solution comes close, so the answer is unambiguous."
                if conclusive
                else "The answer is not unambiguous — see the note below."
            )
        ),
        table=ResultTable(
            columns=(
                Column("spot", "Spot", numeric=True),
                Column("hkl", "Index"),
                Column("d_observed", "d measured", units="Å", numeric=True, digits=4),
                Column("d_calculated", "d calculated", units="Å", numeric=True, digits=4),
                Column(
                    "residual",
                    "Residual",
                    units="Å⁻¹",
                    numeric=True,
                    digits=5,
                    help_text="Distance in reciprocal space between the measured and indexed spot.",
                ),
                Column("x", "x", numeric=True, digits=2),
                Column("y", "y", numeric=True, digits=2),
            ),
            rows=tuple(rows),
            caption=f"Indexed spots of {best.phase_name} down {zone_text}.",
        ),
        data={
            "phase_name": best.phase_name,
            "zone_axis": list(zone_indices),
            "zone_axis_label": zone_text,
            "crystal_to_pattern": np.asarray(best.orientation.as_matrix(), dtype=float)
            .reshape(-1)
            .tolist(),
            "matched_fraction": float(best.matched_fraction),
            "mean_residual_inv_angstrom": float(best.mean_residual_inv_angstrom),
            "max_residual_inv_angstrom": float(best.max_residual_inv_angstrom),
            "unindexed_spots": [int(value) + 1 for value in best.unindexed_spot_indices],
            "seed_spots": [int(value) + 1 for value in best.seed_spot_indices],
            "alternatives": alternatives,
            "conclusive": conclusive,
            "centre": list(centre),
            "describe": report.describe(),
        },
        inputs={
            "phase": spec.to_json(),
            "second_phase": specs[1].to_json() if len(specs) > 1 else None,
            "picks": {"centre": list(centre), "spots": picks},
            "units": request["units"],
            "camera_constant_mm_angstrom": request.get("camera_constant_mm_angstrom"),
            "pixel_size_mm": request.get("pixel_size_mm"),
            "length_tolerance": float(request["length_tolerance"]),
            "angle_tolerance_deg": float(request["angle_tolerance_deg"]),
            "max_index": int(request["max_index"]),
        },
        notes=notes,
        citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON),
    )
    return result.to_json()


@REGISTRY.operation(
    "tem.plan_tilt",
    title="Plan the tilt to another zone axis",
    summary="Alpha and beta tilts to reach a target axis, and whether the stage can.",
    help_text=(
        "Given where the crystal is now — the orientation indexing established, and the stage "
        "angles it was recorded at — this works out the tilts that bring a chosen zone axis onto "
        "the beam.\n\n"
        "It answers three questions a microscopist has to answer before touching the controls. "
        "**Which member?** A cubic ⟨011⟩ has twelve symmetry-equivalent members and they are not "
        "equally far away; the plan ranks them by how much tilting each costs. **Can the stage "
        "get there?** Every solution reports its margin against the tilt envelope, so a target "
        "outside the holder's beta range is stated as unreachable before ten minutes are spent "
        "discovering it. **What will I see on the way?** Each solution carries the path the "
        "crystal traverses, which is what makes it possible to recognise the intermediate axes "
        "and Kikuchi bands as they arrive rather than getting lost between them.\n\n"
        "The uncertainty columns matter. An orientation known to half a degree does not place a "
        "zone axis to half a degree: the holder's two axes are not orthogonal in the crystal, and "
        "the conditioning number reports how badly that amplifies the error for this particular "
        "move."
    ),
    parameters=(
        phase_parameter(),
        IndicesParameter(
            name="current_zone_axis",
            label="Current zone axis [uvw]",
            help_text="The axis on the beam now — normally the one indexing just found.",
            default=(0, 0, 1),
        ),
        IndicesParameter(
            name="target_zone_axis",
            label="Target zone axis [uvw]",
            help_text=(
                "Where you want to be. Any symmetry-equivalent member is acceptable, and the "
                "plan picks the cheapest reachable one."
            ),
            default=(0, 1, 1),
        ),
        NumberParameter(
            name="alpha_deg",
            label="Current alpha",
            help_text="The holder's primary tilt, as the stage reads now.",
            units="°",
            default=0.0,
            group="Stage",
        ),
        NumberParameter(
            name="beta_deg",
            label="Current beta",
            help_text="The holder's secondary tilt, as the stage reads now.",
            units="°",
            default=0.0,
            group="Stage",
        ),
        NumberParameter(
            name="alpha_limit_deg",
            label="Alpha limit",
            help_text=(
                "The holder's alpha range, plus and minus. A double-tilt holder is often ±30°."
            ),
            units="°",
            default=30.0,
            minimum=1.0,
            maximum=90.0,
            group="Stage",
        ),
        NumberParameter(
            name="beta_limit_deg",
            label="Beta limit",
            help_text=(
                "The holder's beta range. Usually smaller than alpha, and usually what makes a "
                "target unreachable."
            ),
            units="°",
            default=20.0,
            minimum=1.0,
            maximum=90.0,
            group="Stage",
        ),
        NumberParameter(
            name="beam_rotation_deg",
            label="Rotation about the beam",
            help_text=(
                "How the crystal is rolled about the beam direction. One indexed pattern cannot "
                "give this — every roll produces the same spot positions — so it comes from a "
                "second pattern at a different tilt, from a Kikuchi pattern, or from the known "
                "orientation of a feature in the image.\n\n"
                "It matters: which member of the target family is nearest and how far the crystal "
                "must turn do not depend on the roll, but how that turn divides between alpha and "
                "beta does, and so therefore does whether the holder can make the move at all."
            ),
            units="°",
            default=0.0,
            group="Stage",
        ),
        NumberParameter(
            name="orientation_uncertainty_deg",
            label="Orientation uncertainty",
            help_text=(
                "How well the current orientation is known. It propagates into the predicted "
                "tilts through the stage geometry, which is what the sigma columns report."
            ),
            units="°",
            default=0.5,
            minimum=0.0,
            maximum=10.0,
            advanced=True,
        ),
        NumberParameter(
            name="tolerance_deg",
            label="On-axis tolerance",
            help_text="How close to exactly on-axis counts as arrived.",
            units="°",
            default=0.5,
            minimum=0.05,
            maximum=5.0,
            advanced=True,
        ),
        BooleanParameter(
            name="allow_reverse",
            label="Allow the opposite sense",
            help_text=(
                "Treat [uvw] and [-u-v-w] as the same destination. Correct for diffraction, "
                "where the two give the same pattern, and usually halves the travel."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns="One row per candidate move; the full report and paths under `data`.",
    panel="tem",
    citations=(
        _CITATION_WILLIAMS,
        "Liu, J. Appl. Crystallogr. 27 (1994) 755 (double-tilt holder geometry).",
    ),
    tags=("TEM", "tilt", "navigation", "zone axis", "stage", "alpha", "beta", "holder"),
)
def _plan_tilt(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import CrystalDirection, ZoneAxis
    from pytex.tem.navigation import plan_tilt_to_zone_axis
    from pytex.tem.reconstruction import CurrentState
    from pytex.tem.stage import DoubleTiltStage, RectangularEnvelope, StagePosition

    spec, phase = phase_from_request(request["phase"])
    current_indices = tuple(request["current_zone_axis"])
    target_indices = tuple(request["target_zone_axis"])

    alpha_limit = float(request["alpha_limit_deg"])
    beta_limit = float(request["beta_limit_deg"])
    stage = DoubleTiltStage(
        envelope=RectangularEnvelope(
            alpha_min_deg=-alpha_limit,
            alpha_max_deg=alpha_limit,
            beta_min_deg=-beta_limit,
            beta_max_deg=beta_limit,
        )
    )
    position = StagePosition(
        alpha_deg=float(request["alpha_deg"]), beta_deg=float(request["beta_deg"])
    )
    current_axis = ZoneAxis(indices=np.asarray(current_indices, dtype=int), phase=phase)
    state = CurrentState.from_orientation(
        _orientation_with_axis_on_beam(
            phase, current_axis, position, roll_deg=float(request["beam_rotation_deg"])
        ),
        position,
        current_zone_axis=current_axis,
        orientation_uncertainty_deg=float(request["orientation_uncertainty_deg"]),
    )
    report = plan_tilt_to_zone_axis(
        state,
        CrystalDirection(coordinates=np.asarray(target_indices, dtype=float), phase=phase),
        stage,
        tolerance_deg=float(request["tolerance_deg"]),
        allow_reverse=bool(request["allow_reverse"]),
    )

    # An unreachable target is a result, not an error. "You can get within 15
    # degrees, at alpha -30" is what a microscopist needs to hear; refusing to
    # answer would send them to tilt blindly and find out the same thing slowly.
    solutions = list(report.solutions)
    exact = bool(solutions)
    if not exact and report.nearest_approach is not None:
        solutions = [report.nearest_approach]

    rows = [
        {
            "member": direction_label(
                tuple(
                    int(value)
                    for value in np.asarray(solution.orbit_member_indices, dtype=int).reshape(-1)
                ),
                spec=spec,
            ),
            "verdict": str(solution.verdict),
            "alpha_deg": float(solution.position.alpha_deg),
            "beta_deg": float(solution.position.beta_deg),
            "delta_alpha_deg": float(solution.delta_alpha_deg),
            "delta_beta_deg": float(solution.delta_beta_deg),
            "travel_deg": float(solution.travel_deg),
            "margin_deg": float(solution.envelope_margin_deg),
            "sigma_alpha_deg": float(solution.sigma_alpha_deg),
            "sigma_beta_deg": float(solution.sigma_beta_deg),
            "conditioning": float(solution.conditioning),
            "residual_deg": float(solution.residual_deg),
        }
        for solution in solutions
    ]
    if not rows:
        raise InvalidInputError(
            f"No member of this target family can be approached within ±{alpha_limit:g}° alpha "
            f"and ±{beta_limit:g}° beta.",
            field="target_zone_axis",
            hint=(
                "Widen the holder limits if the holder really allows it, or choose an "
                "intermediate axis, tilt there, re-index, and plan again."
            ),
        )

    best = rows[0]
    current_text = direction_label(current_indices, spec=spec)
    target_text = direction_label(target_indices, spec=spec)
    route = (
        f"{best['member']}: alpha {best['delta_alpha_deg']:+.2f}° to {best['alpha_deg']:.2f}°, "
        f"beta {best['delta_beta_deg']:+.2f}° to {best['beta_deg']:.2f}°, a crystal rotation of "
        f"{best['travel_deg']:.2f}°"
    )
    result = AppResult(
        title=f"{current_text} → {target_text}",
        summary=(
            (
                f"{report.reachable_orbit_size} of {report.orbit_size} symmetry-equivalent "
                f"members of {target_text} are reachable from {current_text} within the holder's "
                f"range. The cheapest is {route}. Predicted to ±"
                f"{best['sigma_alpha_deg']:.2f}° in alpha and ±{best['sigma_beta_deg']:.2f}° in "
                "beta from the stated orientation uncertainty."
            )
            if exact
            else (
                f"No member of {target_text} is reachable from {current_text} within ±"
                f"{alpha_limit:g}° alpha and ±{beta_limit:g}° beta. The closest the holder can "
                f"come is {route}, leaving the target {best['residual_deg']:.2f}° off axis — "
                "often still enough to see and work with, and the usual answer is to tilt to an "
                "intermediate axis, re-index there, and plan the rest from the new orientation."
            )
        ),
        table=ResultTable(
            columns=(
                Column("member", "Target member"),
                Column("verdict", "Reachable"),
                Column("delta_alpha_deg", "Δα", units="°", numeric=True, digits=2),
                Column("delta_beta_deg", "Δβ", units="°", numeric=True, digits=2),
                Column("alpha_deg", "α", units="°", numeric=True, digits=2),
                Column("beta_deg", "β", units="°", numeric=True, digits=2),
                Column(
                    "travel_deg",
                    "Crystal rotation",
                    units="°",
                    numeric=True,
                    digits=2,
                    help_text="How far the crystal turns, which is not the sum of the tilts.",
                ),
                Column(
                    "margin_deg",
                    "Envelope margin",
                    units="°",
                    numeric=True,
                    digits=2,
                    help_text="How much tilt range is left over. Negative means out of reach.",
                ),
                Column(
                    "residual_deg",
                    "Off axis",
                    units="°",
                    numeric=True,
                    digits=2,
                    help_text="How far the target still is from the beam after this move.",
                ),
                Column("sigma_alpha_deg", "σα", units="°", numeric=True, digits=2),
                Column("sigma_beta_deg", "σβ", units="°", numeric=True, digits=2),
                Column(
                    "conditioning",
                    "Conditioning",
                    numeric=True,
                    digits=2,
                    help_text=(
                        "How much this move amplifies orientation error into tilt error. Large "
                        "means the two holder axes are nearly degenerate for this target."
                    ),
                ),
            ),
            rows=tuple(rows),
            caption=f"Routes from {current_text} to {target_text}.",
        ),
        data={
            "current_zone_axis": list(current_indices),
            "target_zone_axis": list(target_indices),
            "orbit_size": int(report.orbit_size),
            "reachable_orbit_size": int(report.reachable_orbit_size),
            "exact": exact,
            "nearest_approach_deg": (
                float(report.nearest_approach.residual_deg)
                if report.nearest_approach is not None
                else None
            ),
            "waypoints": [
                {
                    "indices": [
                        round(value) for value in np.asarray(waypoint.coordinates, dtype=float)
                    ],
                    "label": direction_label(
                        tuple(
                            round(value) for value in np.asarray(waypoint.coordinates, dtype=float)
                        ),
                        spec=spec,
                    ),
                }
                for waypoint in report.waypoints
            ]
            if report.waypoints
            else [],
            "envelope": {
                "alpha_limit_deg": alpha_limit,
                "beta_limit_deg": beta_limit,
            },
            "start": {
                "alpha_deg": float(request["alpha_deg"]),
                "beta_deg": float(request["beta_deg"]),
            },
            "describe": solutions[0].describe() if solutions else "",
        },
        inputs={
            "phase": spec.to_json(),
            "current_zone_axis": list(current_indices),
            "target_zone_axis": list(target_indices),
            "alpha_deg": float(request["alpha_deg"]),
            "beta_deg": float(request["beta_deg"]),
            "alpha_limit_deg": alpha_limit,
            "beta_limit_deg": beta_limit,
            "beam_rotation_deg": float(request["beam_rotation_deg"]),
            "orientation_uncertainty_deg": float(request["orientation_uncertainty_deg"]),
            "tolerance_deg": float(request["tolerance_deg"]),
            "allow_reverse": bool(request["allow_reverse"]),
        },
        notes=(
            *report.notes,
            "The tilts above hold for the stated rotation about the beam. Which member is "
            "nearest and how far the crystal must turn do not depend on that rotation; how the "
            "turn divides between alpha and beta does, and so does reachability. If the roll is "
            "not known, index a second pattern at a different tilt to fix it.",
        ),
        citations=(_CITATION_WILLIAMS,),
    )
    return result.to_json()


def _orientation_with_axis_on_beam(
    phase: Any, axis: Any, position: Any, *, roll_deg: float = 0.0
) -> Any:
    """Place the crystal with ``axis`` along the beam, rolled by ``roll_deg``.

    Purpose
    -------
    Tilt planning needs a crystal-to-holder orientation, and one indexed pattern
    supplies less than that. It fixes which axis is along the beam and leaves the
    rotation *about* the beam undetermined, because every roll produces the same
    spot positions. The caller supplies the roll — from a second pattern, from a
    Kikuchi pattern, or as an assumption — and this places the crystal
    accordingly.

    What depends on the roll, and what does not, is worth being precise about.
    Which symmetry-equivalent member of a target family is nearest, and how far
    the crystal must turn to reach it, depend only on angles between crystal
    directions and are therefore roll-independent. How that turn divides between
    the holder's alpha and beta is not, and neither, consequently, is whether the
    holder can make the move.

    Parameters
    ----------
    phase : Phase
    axis : ZoneAxis
        The axis to place along the beam.
    position : StagePosition
        Where the stage reads, which fixes the beam direction in the holder.
    roll_deg : float
        Rotation about the beam, positive anticlockwise looking along it.

    Returns
    -------
    Orientation
        Crystal-to-holder, in the holder frame the navigation code expects.
    """

    from pytex.core.orientation import Orientation
    from pytex.tem.reconstruction import HOLDER_FRAME
    from pytex.tem.stage import beam_direction_holder

    forward = np.asarray(axis.unit_vector, dtype=float)
    forward = forward / (np.linalg.norm(forward) or 1.0)
    beam = np.asarray(
        beam_direction_holder(position.alpha_deg, position.beta_deg), dtype=float
    ).reshape(3)
    beam = beam / (np.linalg.norm(beam) or 1.0)
    matrix = _rotation_between(forward, beam)
    if roll_deg:
        matrix = _rotation_about(beam, math.radians(roll_deg)) @ matrix
    return Orientation.from_matrix(
        matrix,
        crystal_frame=phase.crystal_frame,
        specimen_frame=HOLDER_FRAME,
        symmetry=phase.symmetry,
        phase=phase,
    )


def _rotation_about(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation of ``angle_rad`` about a unit ``axis``."""

    unit = np.asarray(axis, dtype=float)
    unit = unit / (np.linalg.norm(unit) or 1.0)
    cross_matrix = np.array(
        [
            [0.0, -unit[2], unit[1]],
            [unit[2], 0.0, -unit[0]],
            [-unit[1], unit[0], 0.0],
        ]
    )
    rotation = (
        np.eye(3)
        + math.sin(angle_rad) * cross_matrix
        + (1.0 - math.cos(angle_rad)) * (cross_matrix @ cross_matrix)
    )
    return np.asarray(rotation, dtype=float)


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """The minimal rotation carrying one unit vector onto another.

    Rodrigues' formula on the cross product, with the antiparallel case handled
    separately because there the axis is undetermined and any perpendicular will
    do.
    """

    axis = np.cross(source, target)
    norm = float(np.linalg.norm(axis))
    cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
    if norm < 1e-12:
        if cosine > 0.0:
            return np.eye(3)
        # Antiparallel: a half turn about any perpendicular axis.
        perpendicular = np.cross(source, np.array([1.0, 0.0, 0.0]))
        if float(np.linalg.norm(perpendicular)) < 1e-9:
            perpendicular = np.cross(source, np.array([0.0, 1.0, 0.0]))
        perpendicular = perpendicular / float(np.linalg.norm(perpendicular))
        return 2.0 * np.outer(perpendicular, perpendicular) - np.eye(3)
    axis = axis / norm
    angle = math.atan2(norm, cosine)
    cross_matrix = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    rotation = (
        np.eye(3)
        + math.sin(angle) * cross_matrix
        + (1.0 - math.cos(angle)) * (cross_matrix @ cross_matrix)
    )
    return np.asarray(rotation, dtype=float)


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="tem.example.tilt_001_to_011",
            title="A 45 degree move a standard holder cannot make",
            panel="tem",
            summary="Austenite [001] to [011] on a ±30°/±20° double-tilt holder.",
            teaches=(
                "⟨011⟩ is 45° from [001], and no combination of ±30° alpha and ±20° beta reaches "
                "45°. The planner says so and reports the closest approach instead of refusing — "
                "which is the useful answer, because a pole 15° off axis is often still workable. "
                "This is the case worth recognising before the tilting starts, not after."
            ),
            operation="tem.plan_tilt",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [0, 1, 1],
                "alpha_limit_deg": 30.0,
                "beta_limit_deg": 20.0,
            },
        ),
        ExampleScenario(
            id="tem.example.tilt_roll_matters",
            title="The same move, reachable — because of the roll",
            panel="tem",
            summary="The identical target on a ±40°/±40° holder, rolled 45° about the beam.",
            teaches=(
                "Change only the rotation about the beam and eight members of ⟨011⟩ become "
                "reachable where none were. That is the point of the roll field: one indexed "
                "pattern cannot give it, every roll produces identical spot positions, and yet "
                "reachability depends on it entirely. Set it from a second pattern rather than "
                "guessing, and read the waypoint list for the axes on the way."
            ),
            operation="tem.plan_tilt",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [0, 1, 1],
                "alpha_limit_deg": 40.0,
                "beta_limit_deg": 40.0,
                "beam_rotation_deg": 45.0,
            },
        ),
        ExampleScenario(
            id="tem.example.tilt_111",
            title="Reaching [111] from [001]",
            panel="tem",
            summary="A 54.7° move that needs a wide holder, on a ±60° stage.",
            teaches=(
                "A given zone axis corresponds to essentially one (alpha, beta) pair, not to a "
                "choice of routes that can be traded off — so reachability is that single pair "
                "against the envelope, and there is no combining a large alpha with a large beta "
                "to reach further. Here the move is 54.7° and lands almost entirely on one axis, "
                "which is why a ±30° holder cannot make it however the tilts are divided."
            ),
            operation="tem.plan_tilt",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [1, 1, 1],
                "alpha_limit_deg": 60.0,
                "beta_limit_deg": 60.0,
                "beam_rotation_deg": 45.0,
            },
        ),
        ExampleScenario(
            id="tem.example.tilt_hcp",
            title="Basal to prism in zirconium",
            panel="tem",
            summary="A 90° move in a hexagonal crystal.",
            teaches=(
                "Basal and prism axes are exactly 90° apart whatever c/a is, so this is beyond "
                "any conventional holder in one step and the planner reports the nearest "
                "approach. The honest route is the waypoint list: tilt to an intermediate axis, "
                "re-index there, and plan the rest from the new orientation — which also resets "
                "the accumulated calibration error rather than compounding it."
            ),
            operation="tem.plan_tilt",
            request={
                "phase": {"builtin": "zr_hcp"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [1, 0, 0],
                "alpha_limit_deg": 60.0,
                "beta_limit_deg": 60.0,
            },
        ),
    )
)
