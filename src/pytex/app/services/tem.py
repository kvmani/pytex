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
from pytex.app.services.calculator import (
    direction_label,
    family_label,
    phase_parameter,
    plane_label,
)
from pytex.app.tem_gallery import GALLERY, gallery_entry, gallery_options

__all__ = ["measured_pattern_from_picks"]

_CITATION_HIRSCH = (
    "Hirsch et al., Electron Microscopy of Thin Crystals (1965), appendix on zone-axis patterns."
)

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
        "**Every candidate carries its calculated pattern.** Superimposing what a solution "
        "*predicts* on what was *measured* turns accepting it into a judgement made by looking: a "
        "calculated pattern uniformly too large is a camera constant, one turned is a roll, one "
        "with rows the plate does not show is the wrong phase. The prediction is bounded by the "
        "index limit below, so a plate spot with no calculated node beside it means "
        "*check the index limit* before doubting the solution.\n\n"
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
        NumberParameter(
            name="score_length_weight",
            label="Weight: d-spacing agreement",
            help_text=(
                "How much a d-spacing disagreement counts in the fused accuracy score. Sensitive "
                "to the camera constant, which is the one calibration that can be wrong while "
                "everything else stays self-consistent."
            ),
            default=1.0,
            minimum=0.0,
            maximum=10.0,
            group="Scoring",
            advanced=True,
        ),
        NumberParameter(
            name="score_angle_weight",
            label="Weight: angle agreement",
            help_text=(
                "How much an interspot-angle disagreement counts. Weighted above d-spacings by "
                "default because angles are calibration-free: a wrong camera constant scales "
                "every length and leaves every angle alone, so an angular disagreement is "
                "evidence about the crystallography rather than about the instrument."
            ),
            default=1.5,
            minimum=0.0,
            maximum=10.0,
            group="Scoring",
            advanced=True,
        ),
        NumberParameter(
            name="score_coverage_weight",
            label="Weight: spots explained",
            help_text=(
                "How much the fraction of picked spots indexed counts. Weighted highest by "
                "default: an unindexed spot is unexplained evidence, and precision on the others "
                "does not answer it."
            ),
            default=2.0,
            minimum=0.0,
            maximum=10.0,
            group="Scoring",
            advanced=True,
        ),
        NumberParameter(
            name="score_length_tolerance",
            label="Half-score d-spacing deviation",
            help_text=(
                "The relative d-spacing deviation that scores one half on its term. 0.02 is two "
                "percent, about what a well-calibrated instrument achieves."
            ),
            default=0.02,
            minimum=0.0001,
            maximum=0.5,
            group="Scoring",
            advanced=True,
        ),
        NumberParameter(
            name="score_angle_tolerance_deg",
            label="Half-score angle deviation",
            help_text=(
                "The angular deviation that scores one half on its term. One degree is roughly "
                "the precision of picking two spots and measuring the angle between them."
            ),
            units="deg",
            default=1.0,
            minimum=0.01,
            maximum=30.0,
            group="Scoring",
            advanced=True,
        ),
        IndicesParameter(
            name="expected_zone_axis",
            label="Check against axis [uvw]",
            help_text=(
                "Optional. If you already know which axis this pattern is down — a practice "
                "plate, a reference specimen, a pattern you have indexed before — give it here "
                "and the result states whether the indexing agrees, measuring the disagreement in "
                "degrees. The comparison is made up to symmetry, because a bcc [110] pattern is "
                "indistinguishable from a [101] one and calling that a mismatch would be wrong."
            ),
            required=False,
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
    returns=(
        "One row per spot with its index and its d-spacing deviation; the scored, ranked "
        "solutions and their calculated patterns under `data`."
    ),
    panel="tem",
    citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON),
    tags=("TEM", "SAED", "indexing", "solve", "zone axis", "pattern", "calibration", "score"),
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
    from pytex.diffraction.solution_scoring import ScoringWeights, score_solution

    weights = ScoringWeights(
        length=float(request["score_length_weight"]),
        angle=float(request["score_angle_weight"]),
        coverage=float(request["score_coverage_weight"]),
        length_tolerance=float(request["score_length_tolerance"]),
        angle_tolerance_deg=float(request["score_angle_tolerance_deg"]),
    )
    measured_g = np.asarray(pattern.g_vectors_inv_angstrom(), dtype=float)
    scored = [
        (solution, score_solution(solution, measured_g, weights=weights))
        for solution in report.solutions
    ]
    # Rank by the fused score rather than by the solver own sort key. That key
    # orders by matched fraction then residual and is explicitly not a quality;
    # the score is one, it is the thing the user configured, and a list sorted
    # by something other than the number printed beside it would be a trap.
    ranked = sorted(scored, key=lambda item: -item[1].score)
    best, best_score = ranked[0]
    reordered = best is not report.solutions[0]

    spec_by_name = {item.name: item for item in specs}
    best_spec = spec_by_name.get(best.phase_name, spec)
    observed_d = np.asarray(pattern.d_spacings_angstrom(), dtype=float)
    rows: list[dict[str, Any]] = []
    for spot in best.solved_spots:
        index = int(spot.measured_index)
        indices = tuple(int(value) for value in np.asarray(spot.hkl, dtype=int))
        # `predicted_g_inv_angstrom` is the in-plane vector, not its length.
        predicted = float(np.linalg.norm(np.asarray(spot.predicted_g_inv_angstrom, dtype=float)))
        calculated_d = 1.0 / predicted if predicted > 0.0 else float("nan")
        rows.append(
            {
                "spot": index + 1,
                "hkl": plane_label(indices, spec=best_spec) if any(indices) else "—",
                "d_observed": float(observed_d[index]),
                "d_calculated": calculated_d,
                "d_deviation_percent": (
                    100.0 * (observed_d[index] - calculated_d) / calculated_d
                    if calculated_d > 0.0
                    else float("nan")
                ),
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
            # `PatternSolution.zone_axis_label` is the solver's own three-index
            # rendering. Everything the user sees elsewhere goes through
            # `direction_label`, which writes four indices for a hexagonal phase,
            # and a candidate list that says [010] beside a title reading
            # [1̄21̄0] is describing the same axis in two notations.
            "zone_axis": direction_label(
                tuple(int(value) for value in np.asarray(solution.zone_axis.indices, dtype=int)),
                spec=spec_by_name.get(solution.phase_name, spec),
            ),
            "zone_axis_indices": [
                int(value) for value in np.asarray(solution.zone_axis.indices, dtype=int)
            ],
            "matched_spots": round(solution.matched_fraction * solution.measured_spot_count),
            "matched_fraction": float(solution.matched_fraction),
            "mean_residual_inv_angstrom": float(solution.mean_residual_inv_angstrom),
            "score": float(item.score),
            "length_agreement": float(item.length_agreement),
            "angle_agreement": float(item.angle_agreement),
            "coverage_agreement": float(item.coverage_agreement),
            "rms_relative_length_deviation": float(item.rms_relative_length_deviation),
            "rms_angle_deviation_deg": float(item.rms_angle_deviation_deg),
            # Everything the browser needs to draw this candidate over the
            # measured pattern, so accepting a solution can be a judgement made
            # by looking rather than by reading a residual column.
            "overlay": _calculated_overlay(
                spec_by_name.get(solution.phase_name, spec),
                solution,
                request=request,
                centre=centre,
            ),
            "describe": item.describe(),
        }
        for solution, item in ranked
    ]
    conclusive = bool(report.is_conclusive)
    notes: list[str] = []
    check: dict[str, Any] | None = None
    if request.get("expected_zone_axis"):
        expected = tuple(int(value) for value in request["expected_zone_axis"])
        deviation = _symmetry_angle_deg(phase, expected, zone_indices)
        expected_text = direction_label(expected, spec=best_spec)
        # Half a degree is the same tolerance the tilt planner treats as
        # on-axis; anything inside it is the same zone, differently rounded.
        matched = deviation <= 0.5
        check = {
            "expected": list(expected),
            "expected_label": expected_text,
            "deviation_deg": deviation,
            "correct": matched,
        }
        notes.append(
            (
                f"This agrees with the expected axis {expected_text}, up to symmetry — the two "
                f"differ by {deviation:.3f}°."
            )
            if matched
            else (
                f"This does *not* agree with the expected axis {expected_text}: the indexed axis "
                f"is {deviation:.2f}° away from it, which is more than a rounding difference. "
                "Check the beam pick first — an error there biases every d-spacing — then the "
                "camera constant, then whether every picked spot really is a reflection."
            )
        )
    if reordered:
        notes.append(
            "The accuracy score puts a different solution first than the solver ordering does. "
            "The two ask different questions - the solver ranks by how many spots were indexed "
            "and how tightly, the score weighs d-spacings, angles and coverage by the policy you "
            "set - and a disagreement between them is a sign the pattern does not settle the "
            "answer on its own."
        )
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
                    "d_deviation_percent",
                    "Δd",
                    units="%",
                    numeric=True,
                    digits=2,
                    help_text=(
                        "Measured minus calculated, as a percentage. The same sign on every spot "
                        "is the signature of a wrong camera constant rather than a wrong indexing."
                    ),
                ),
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
            "score": best_score.to_json(),
            "reordered_by_score": reordered,
            "conclusive": conclusive,
            "check": check,
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
            "score_length_weight": float(request["score_length_weight"]),
            "score_angle_weight": float(request["score_angle_weight"]),
            "score_coverage_weight": float(request["score_coverage_weight"]),
            "score_length_tolerance": float(request["score_length_tolerance"]),
            "score_angle_tolerance_deg": float(request["score_angle_tolerance_deg"]),
            "expected_zone_axis": (
                list(request["expected_zone_axis"]) if request.get("expected_zone_axis") else None
            ),
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
            "member": _member_label(solution, target_indices, spec)[0],
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
    # The repository notation standard is explicit: a specific direction is
    # [uvw] and a symmetry family is <uvw>. "Members of [012]" is a category
    # error — [012] has no members — so the sentence that counts members writes
    # the family form even though the title keeps the direction the user typed.
    target_family = family_label(target_indices, spec=spec, family="direction")
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
                f"members of {target_family} are reachable from {current_text} within the "
                f"holder's "
                f"range. The cheapest is {route}. Predicted to ±"
                f"{best['sigma_alpha_deg']:.2f}° in alpha and ±{best['sigma_beta_deg']:.2f}° in "
                "beta from the stated orientation uncertainty."
            )
            if exact
            else (
                f"No member of {target_family} is reachable from {current_text} within ±"
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


@REGISTRY.operation(
    "tem.fit_lattice",
    title="Fit a lattice to the picks",
    summary="Refine the beam centre from the spots, and show which picks do not belong.",
    help_text=(
        "A zone-axis pattern is a plane of the reciprocal lattice, so its spots lie on a "
        "two-dimensional lattice. Imposing that — before any phase, camera constant or zone axis "
        "enters — buys two things.\n\n"
        "**The beam centre stops being a guess.** Picking the transmitted beam by eye is the "
        "largest avoidable error in the whole workflow: it biases every d-spacing at once, and it "
        "does so while leaving the pattern self-consistent, so the result is a plausible answer "
        "for the wrong material rather than an obvious failure. But four or more spots give more "
        "equations than unknowns, so the centre can be *solved for* instead of clicked.\n\n"
        "**A mis-picked spot becomes visible.** The fitted lattice is drawn over the pattern; a "
        "spot clicked one node out, or clicked on a dust particle, stops matching the grid. That "
        "is a judgement anyone makes instantly from a picture, and slowly from residuals.\n\n"
        "**This is geometry, not indexing.** A good fit says the picks are consistent with some "
        "lattice, which is necessary for a correct indexing and nowhere near sufficient. It does "
        "not know what phase this is and cannot tell you.\n\n"
        "Two limits are worth knowing. A centre wrong by an exact lattice vector is undetectable "
        "here — every spot is still an exact node — and what identifies the transmitted beam is "
        "that it is the brightest thing on the plate. And a centre more than half a spacing out "
        "cannot be refined without changing which node the origin is, so the fit is held there and "
        "says so."
    ),
    parameters=(
        ObjectParameter(
            name="picks",
            label="Picked spots",
            help_text=(
                "The transmitted beam and the reflections. The beam is the starting point for the "
                "refinement, not the answer."
            ),
            editor="spot-picker",
        ),
        BooleanParameter(
            name="refine_centre",
            label="Refine the beam centre",
            help_text=(
                "Solve for the centre along with the lattice. Turn it off to hold a centre "
                "established another way — a measured beam-stop position, for instance."
            ),
            default=True,
        ),
        NumberParameter(
            name="inlier_fraction",
            label="Pick tolerance",
            help_text=(
                "How close a click must land to a lattice node to count as that node, as a "
                "fraction of the shortest separation between two picked spots. It stands in for "
                "picking precision: 0.04 is a few pixels on a typical plate."
            ),
            default=0.04,
            minimum=0.001,
            maximum=0.5,
            advanced=True,
        ),
        IntegerParameter(
            name="node_limit",
            label="Overlay extent",
            help_text="Largest lattice index drawn in the overlay, in each direction.",
            default=6,
            minimum=1,
            maximum=20,
            advanced=True,
        ),
        NumberParameter(
            name="frame_width",
            label="Frame width",
            help_text="Image width, so the overlay carries no nodes outside the picture.",
            default=1024.0,
            minimum=1.0,
            required=False,
            advanced=True,
        ),
        NumberParameter(
            name="frame_height",
            label="Frame height",
            help_text="Image height, so the overlay carries no nodes outside the picture.",
            default=1024.0,
            minimum=1.0,
            required=False,
            advanced=True,
        ),
    ),
    returns="One row per picked spot with its lattice node and residual; the overlay under `data`.",
    panel="tem",
    citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON),
    tags=("TEM", "SAED", "lattice", "centre", "refine", "picks", "overlay"),
)
def _fit_lattice(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.diffraction.lattice_fit import fit_planar_lattice

    centre, picks = _picks(request.get("picks"))
    positions = np.asarray([[spot["x"], spot["y"]] for spot in picks], dtype=float)
    try:
        fit = fit_planar_lattice(
            positions,
            centre,
            refine_centre=bool(request["refine_centre"]),
            inlier_fraction=float(request["inlier_fraction"]),
        )
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            field="picks",
            hint=(
                "A lattice needs two directions. Pick at least one spot off the row you already "
                "have, and check that the transmitted beam was marked first."
            ),
        ) from error

    rows = [
        {
            "spot": spot.index + 1,
            "node": f"({spot.lattice_indices[0]}, {spot.lattice_indices[1]})",
            "x": float(spot.position[0]),
            "y": float(spot.position[1]),
            "predicted_x": float(spot.predicted[0]),
            "predicted_y": float(spot.predicted[1]),
            "residual": float(spot.residual),
            "verdict": "on the lattice" if spot.inlier else "off the lattice",
        }
        for spot in fit.spots
    ]
    # The two basis vectors, aimed at the picks that actually generate them.
    #
    # Pointing the arrow at the *picked spot* rather than at the ideal node is
    # the whole value of drawing it: the arrow then moves with the spot as it is
    # adjusted, and the gap between where it lands and where the lattice says it
    # should be is visible directly. When no pick sits on a unit node — the user
    # picked only second-order reflections, say — the ideal node is used instead
    # and the row says so, because an arrow to a spot that is not there would be
    # a lie about which pick is driving the fit.
    by_node = {spot.lattice_indices: spot for spot in fit.spots}
    basis_vectors: list[dict[str, Any]] = []
    for label, node, opposite in (("a", (1, 0), (-1, 0)), ("b", (0, 1), (0, -1))):
        spot = by_node.get(node)
        sense = 1
        if spot is None:
            spot = by_node.get(opposite)
            sense = -1
        row = 0 if label == "a" else 1
        ideal = np.asarray(fit.centre, dtype=float) + sense * np.asarray(fit.basis[row])
        tip = np.asarray(spot.position, dtype=float) if spot is not None else ideal
        basis_vectors.append(
            {
                "label": label if sense > 0 else f"-{label}",
                "from": [float(value) for value in fit.centre],
                "to": [float(tip[0]), float(tip[1])],
                "length": float(np.linalg.norm(tip - np.asarray(fit.centre, dtype=float))),
                "spot": (spot.index + 1) if spot is not None else None,
                "on_a_pick": spot is not None,
            }
        )

    bounds = (
        float(request["frame_width"]) if request.get("frame_width") else None,
        float(request["frame_height"]) if request.get("frame_height") else None,
    )
    nodes = fit.node_positions(
        max_index=int(request["node_limit"]),
        bounds=(bounds[0], bounds[1]) if bounds[0] and bounds[1] else None,
    )
    short, long_ = fit.basis_lengths

    result = AppResult(
        title="Lattice fitted to the picks",
        summary=(
            f"The {len(fit.spots)} picked spots fit a lattice with sides of {short:.1f} and "
            f"{long_:.1f} picking units at {fit.basis_angle_deg:.2f} degrees, explaining "
            f"{fit.inlier_count} of them to an r.m.s. of {fit.rms_residual:.2f} units. "
            + (
                f"The transmitted beam is held where it was picked, at ({fit.centre[0]:.1f}, "
                f"{fit.centre[1]:.1f})."
                if not fit.centre_refined
                else (
                    f"The transmitted beam refines to ({fit.centre[0]:.1f}, {fit.centre[1]:.1f}), "
                    f"{fit.centre_shift:.1f} units from where it was picked."
                    if fit.centre_shift > 0.05
                    else "The transmitted beam is already where the spots say it should be."
                )
            )
        ),
        table=ResultTable(
            columns=(
                Column("spot", "Spot", numeric=True),
                Column(
                    "node",
                    "Lattice node",
                    help_text="Which (m, n) node of the fitted lattice this pick was assigned to.",
                ),
                Column("x", "x", numeric=True, digits=1),
                Column("y", "y", numeric=True, digits=1),
                Column("predicted_x", "Node x", numeric=True, digits=1),
                Column("predicted_y", "Node y", numeric=True, digits=1),
                Column(
                    "residual",
                    "Residual",
                    numeric=True,
                    digits=2,
                    help_text="Distance from the pick to its node, in picking units.",
                ),
                Column("verdict", "Verdict"),
            ),
            rows=tuple(rows),
            caption="Every pick, and the lattice node it was assigned to.",
        ),
        data={
            "fit": fit.to_json(),
            "centre": [float(value) for value in fit.centre],
            "supplied_centre": [float(value) for value in fit.supplied_centre],
            "centre_shift": float(fit.centre_shift),
            # Whether `centre` is a measurement or an echo of the click. The
            # overlay draws the grid from it and the panel offers to adopt it, so
            # a caller that cannot tell the two apart will present a clicked
            # position as a refined one — which is the misreading this whole
            # operation exists to prevent.
            "centre_refined": bool(fit.centre_refined),
            "nodes": [
                {"x": float(row[0]), "y": float(row[1]), "m": int(row[2]), "n": int(row[3])}
                for row in nodes
            ],
            "basis_vectors": basis_vectors,
            "outliers": [spot.index + 1 for spot in fit.outliers],
            "describe": fit.describe(),
        },
        inputs={
            "picks": {"centre": list(centre), "spots": picks},
            "refine_centre": bool(request["refine_centre"]),
            "inlier_fraction": float(request["inlier_fraction"]),
        },
        notes=(
            *fit.notes,
            "This is geometry, not indexing: a lattice that fits says the picks are mutually "
            "consistent, which is necessary for a correct indexing and far from sufficient.",
        ),
        citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON),
    )
    return result.to_json()


@REGISTRY.operation(
    "tem.gallery_pattern",
    title="Open a practice pattern",
    summary="A simulated SAED plate with a known answer, ready to pick and index.",
    help_text=(
        "Loads one of three simulated diffraction patterns onto the picking canvas, so the whole "
        "indexing workflow can be run without a micrograph — and, because the pattern was built "
        "from a known zone axis, checked afterwards against the answer.\n\n"
        "**These are calculations, not pictures.** Every spot sits exactly where the lattice, the "
        "zone axis and the camera constant put it: radial position is the camera constant divided "
        "by the d-spacing, which is the one identity a real pattern is calibrated by. Relative "
        "brightness is kinematic and therefore indicative; double diffraction is not modelled, so "
        "a reflection that a real plate shows through double diffraction is absent here.\n\n"
        "**The realism is deliberate.** The beam is not always at the centre of the frame, the "
        "pattern is rolled about the beam so it does not line up with the detector axes, and each "
        "spot carries a small centroiding scatter. A workflow that only works on an idealised "
        "pattern has not been tested.\n\n"
        "**The calibration is computed, not typed.** The camera constant is the camera length "
        "times the relativistic electron wavelength at the stated accelerating voltage — change "
        "either and the pattern rescales, which is exactly what happens at the microscope."
    ),
    parameters=(
        ChoiceParameter(
            name="pattern",
            label="Practice pattern",
            help_text="Which of the gallery plates to open.",
            options=gallery_options(),
            default=GALLERY[0].identifier,
        ),
        NumberParameter(
            name="camera_length_mm",
            label="Camera length",
            help_text=(
                "The projected camera length L. With the accelerating voltage it fixes the camera "
                "constant L·λ, and therefore how much of reciprocal space fits on the plate: a "
                "longer camera spreads the pattern out and loses the outer reflections."
            ),
            units="mm",
            default=400.0,
            minimum=50.0,
            maximum=4000.0,
            group="Instrument",
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text=(
                "Sets the relativistic electron wavelength. A higher voltage shortens λ and so "
                "shrinks the pattern at fixed camera length."
            ),
            units="kV",
            default=200.0,
            minimum=20.0,
            maximum=1000.0,
            group="Instrument",
        ),
        NumberParameter(
            name="extra_rotation_deg",
            label="Extra roll about the beam",
            help_text=(
                "Added to the entry's own roll. Turn it to get a fresh exercise from the same "
                "pattern: the indexed answer must not change, because one pattern cannot fix the "
                "rotation about the beam in the first place."
            ),
            units="°",
            default=0.0,
            group="Instrument",
        ),
        BooleanParameter(
            name="realistic_scatter",
            label="Include centroiding scatter",
            help_text=(
                "Displace each spot by a fraction of a pixel, as a measurement of a real plate "
                "would be. Switch it off for an exactly constructed pattern, where indexing "
                "residuals fall to machine precision and stop being informative."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns="One row per simulated spot; the pattern, the answer and the calibration under `data`.",
    panel="tem",
    citations=(_CITATION_WILLIAMS, _CITATION_HIRSCH),
    tags=("TEM", "SAED", "gallery", "practice", "simulation", "pattern", "teaching"),
)
def _gallery_pattern(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.kinematic import electron_wavelength_angstrom
    from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

    entry = gallery_entry(str(request["pattern"]))
    spec = entry.phase_spec()
    phase = spec.to_phase()
    camera_length = float(request["camera_length_mm"])
    beam_energy = float(request["beam_energy_kev"])
    camera_constant = float(camera_length * electron_wavelength_angstrom(beam_energy))
    roll = float(entry.in_plane_rotation_deg) + float(request["extra_rotation_deg"])
    image = synthesize_saed_image(
        phase,
        ZoneAxis(indices=np.asarray(entry.zone_axis, dtype=int), phase=phase),
        camera_constant_mm_angstrom=camera_constant,
        raster=DetectorRaster(
            width_px=entry.width_px,
            height_px=entry.height_px,
            pixel_size_mm=entry.pixel_size_mm,
            centre_px=entry.centre_px,
        ),
        in_plane_rotation_deg=roll,
        position_jitter_px=entry.jitter_px if request["realistic_scatter"] else 0.0,
        rng_seed=entry.rng_seed,
    )
    if not image.spots:
        raise InvalidInputError(
            "No reflection of this zone falls on the detector at that camera length.",
            field="camera_length_mm",
            hint=(
                f"At {camera_length:g} mm the innermost reflection is already off the plate. "
                "Shorten the camera length until the pattern fits."
            ),
        )

    zone_text = direction_label(tuple(entry.zone_axis), spec=spec)
    rows = [
        {
            "spot": index + 1,
            "hkl": plane_label(tuple(int(value) for value in spot.miller_indices), spec=spec),
            "d": float(spot.d_spacing_angstrom),
            "g": float(spot.g_inv_angstrom),
            "intensity": float(spot.relative_intensity),
            "x": float(spot.position_px[0]),
            "y": float(spot.position_px[1]),
        }
        for index, spot in enumerate(image.spots)
    ]

    payload = image.to_json()
    # The browser draws from the same numbers the table shows, and labels a
    # hexagonal reflection in four indices as the literature does; the library
    # deals only in the three-index basis it computes with.
    for spot_payload, row in zip(payload["spots"], rows, strict=True):
        spot_payload["label"] = row["hkl"]

    suggested = image.independent_seed_spots(6)
    targets = [
        {
            **target.to_json(),
            "label": direction_label(tuple(target.indices), spec=spec),
        }
        for target in entry.targets
    ]

    result = AppResult(
        title=f"{entry.title}",
        summary=(
            f"A simulated {spec.name} pattern with the beam along {zone_text}, recorded at "
            f"{beam_energy:g} kV with a {camera_length:g} mm camera length — a camera constant of "
            f"{camera_constant:.4f} mm·Å. {len(rows)} reflections fall on the "
            f"{entry.width_px}×{entry.height_px} detector. Click the transmitted beam, then the "
            "spots, and index it; the answer is known, so the panel can tell you whether you got "
            "it right."
        ),
        table=ResultTable(
            columns=(
                Column("spot", "Spot", numeric=True),
                Column("hkl", "Index"),
                Column("d", "d", units="Å", numeric=True, digits=4),
                Column("g", "|g|", units="Å⁻¹", numeric=True, digits=4),
                Column(
                    "intensity",
                    "Relative intensity",
                    numeric=True,
                    digits=3,
                    help_text=(
                        "Kinematic, normalized to the strongest spot. Indicative rather than "
                        "quantitative — a real plate redistributes intensity dynamically."
                    ),
                ),
                Column("x", "x", units="px", numeric=True, digits=1),
                Column("y", "y", units="px", numeric=True, digits=1),
            ),
            rows=tuple(rows),
            caption=f"Simulated reflections of {spec.name} down {zone_text}.",
        ),
        data={
            "entry": {
                "id": entry.identifier,
                "title": entry.title,
                "summary": entry.summary,
                "teaches": entry.teaches,
                "phase_key": entry.phase_key,
            },
            "pattern": payload,
            "zone_axis_label": zone_text,
            "phase_name": spec.name,
            # Everything the solver needs, so a user never has to transcribe a
            # calibration from this panel into the one beside it.
            "calibration": {
                "units": "px",
                "camera_constant_mm_angstrom": camera_constant,
                "pixel_size_mm": float(entry.pixel_size_mm),
                # A catalogue *reference*, not an expanded description. The phase
                # picker treats a full description as a user-edited phase and
                # renames it "(edited)", which is both untrue and alarming on a
                # pattern the user has not touched.
                "phase": {"builtin": entry.phase_key},
            },
            "suggested_picks": {
                "centre": [float(value) for value in image.centre_px],
                "spots": [
                    {"x": float(spot.position_px[0]), "y": float(spot.position_px[1])}
                    for spot in suggested
                ],
            },
            "targets": targets,
            "describe": image.describe(),
        },
        inputs={
            "pattern": entry.identifier,
            "camera_length_mm": camera_length,
            "beam_energy_kev": beam_energy,
            "extra_rotation_deg": float(request["extra_rotation_deg"]),
            "realistic_scatter": bool(request["realistic_scatter"]),
        },
        notes=(
            "Spot positions are exact for this lattice, zone axis and camera constant. Relative "
            "intensities are kinematic, and double diffraction is not modelled — so a forbidden "
            "reflection that a real plate shows will be missing here.",
        ),
        citations=(_CITATION_WILLIAMS, _CITATION_HIRSCH),
    )
    return result.to_json()


@REGISTRY.operation(
    "tem.zone_axis_atlas",
    title="Find the zone axes worth going to",
    summary="Nearby axes of this phase: pattern, angle, and whether the holder can get there.",
    help_text=(
        "The other half of tilt planning. Planning a tilt answers *can I reach the axis I named*; "
        "this answers *which axis should I name* — which is the question actually asked at the "
        "column, with the beam already down something and the specimen still uncharacterised.\n\n"
        "Each row is one symmetry-distinct zone-axis family, with four things that decide the "
        "choice. **How far** it is, measured to the nearest member of the family rather than to "
        "the one you typed, because they are all the same pattern. **How many members** it has: a "
        "family of twelve offers twelve chances that one of them is inside the holder's range. "
        "**How rich** its pattern is, as the number of reflections inside a fixed cut-off — a "
        "trip to a two-spot zone buys very little. **What it looks like**, as the n-fold symmetry "
        "you will recognise on the screen when you arrive, which is the first confirmation that "
        "you arrived where you meant to.\n\n"
        "Reachability is computed by the same planner the tilt panel uses, against the same holder "
        "envelope and rotation about the beam, so a row marked reachable here is reachable there."
    ),
    parameters=(
        phase_parameter(help_text="The phase whose zone axes to enumerate."),
        IndicesParameter(
            name="current_zone_axis",
            label="Current zone axis [uvw]",
            help_text="The axis on the beam now. Angles in the table are measured from it.",
            default=(0, 0, 1),
        ),
        NumberParameter(
            name="max_angle_deg",
            label="Search within",
            help_text=(
                "Ignore families farther than this from the current axis. 60° reaches every "
                "low-index cubic axis — ⟨110⟩ at 45° and ⟨111⟩ at 54.74° — while staying short "
                "of the far side of the stereographic triangle."
            ),
            units="°",
            default=60.0,
            minimum=1.0,
            maximum=90.0,
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
            help_text="The holder's alpha range, plus and minus.",
            units="°",
            default=30.0,
            minimum=1.0,
            maximum=90.0,
            group="Stage",
        ),
        NumberParameter(
            name="beta_limit_deg",
            label="Beta limit",
            help_text="The holder's beta range. Usually what makes a target unreachable.",
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
                "How the crystal is rolled about the beam. It does not change which axes are "
                "near, but it does change how the move divides between alpha and beta, and "
                "therefore which rows come back reachable."
            ),
            units="°",
            default=0.0,
            group="Stage",
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text=(
                "Largest |u|, |v| or |w| considered. 2 gives the axes a standard stereogram "
                "labels — ⟨100⟩, ⟨110⟩, ⟨111⟩, ⟨210⟩, ⟨211⟩, ⟨221⟩ — which is what a session "
                "normally uses. Raising it admits axes such as ⟨320⟩ and ⟨331⟩: nearer, but with "
                "sparse patterns that are harder to recognise and buy less information."
            ),
            default=2,
            minimum=1,
            maximum=5,
            advanced=True,
        ),
        IntegerParameter(
            name="limit",
            label="Rows",
            help_text="How many families to return, after ranking by distance.",
            default=12,
            minimum=1,
            maximum=40,
            advanced=True,
        ),
    ),
    returns="One row per zone-axis family, nearest first; the full atlas under `data`.",
    panel="tem",
    citations=(_CITATION_WILLIAMS, _CITATION_EDINGTON, _CITATION_HIRSCH),
    tags=("TEM", "zone axis", "atlas", "navigation", "tilt", "stage", "planning"),
)
def _zone_axis_atlas(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import CrystalDirection, ZoneAxis
    from pytex.tem.atlas import zone_axis_atlas
    from pytex.tem.navigation import plan_tilt_to_zone_axis
    from pytex.tem.reconstruction import CurrentState
    from pytex.tem.stage import DoubleTiltStage, RectangularEnvelope, StagePosition

    spec, phase = phase_from_request(request["phase"])
    current_indices = tuple(request["current_zone_axis"])
    current_axis = ZoneAxis(indices=np.asarray(current_indices, dtype=int), phase=phase)
    limit = int(request["limit"])

    atlas = zone_axis_atlas(
        phase,
        current_zone_axis=current_axis,
        max_index=int(request["max_index"]),
        max_angle_deg=float(request["max_angle_deg"]),
        limit=limit,
    )

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
    state = CurrentState.from_orientation(
        _orientation_with_axis_on_beam(
            phase, current_axis, position, roll_deg=float(request["beam_rotation_deg"])
        ),
        position,
        current_zone_axis=current_axis,
    )

    rows: list[dict[str, Any]] = []
    reachable_count = 0
    for entry in atlas.entries:
        indices = tuple(int(value) for value in entry.indices)
        row: dict[str, Any] = {
            "family": family_label(indices, spec=spec, family="direction"),
            "target": direction_label(indices, spec=spec),
            "angle_deg": float(entry.angle_from_current_deg),
            "family_size": int(entry.family_size),
            "reflections": int(entry.reflection_count),
            "symmetry": f"{entry.rotational_order}-fold",
            "verdict": "current axis",
            "reachable": False,
            "delta_alpha_deg": 0.0,
            "delta_beta_deg": 0.0,
            "travel_deg": 0.0,
            "margin_deg": float("nan"),
        }
        # A thousandth of a degree, not machine epsilon: the angle is an arccos
        # near 1, where the square-root behaviour turns 1e-16 of cosine error
        # into ~1e-6 of a degree — enough for the axis already on the beam to
        # miss a tighter test and be reported as somewhere worth tilting to.
        if entry.angle_from_current_deg > 1e-3:
            report = plan_tilt_to_zone_axis(
                state,
                CrystalDirection(coordinates=np.asarray(indices, dtype=float), phase=phase),
                stage,
            )
            solution = report.solutions[0] if report.solutions else report.nearest_approach
            if solution is None:
                row["verdict"] = "no solution"
            else:
                label, member = _member_label(solution, indices, spec)
                row["target"] = label
                # The same four verdicts the tilt panel reports, so a row here
                # and the plan it produces there cannot describe the same move
                # in two different vocabularies.
                row["verdict"] = str(solution.verdict)
                row["delta_alpha_deg"] = float(solution.delta_alpha_deg)
                row["delta_beta_deg"] = float(solution.delta_beta_deg)
                row["travel_deg"] = float(solution.travel_deg)
                row["margin_deg"] = float(solution.envelope_margin_deg)
                row["indices"] = member
                row["reachable"] = bool(report.solutions)
                if report.solutions:
                    reachable_count += 1
        row.setdefault("indices", list(indices))
        rows.append(row)

    if not rows:
        raise InvalidInputError(
            "No zone-axis family of this phase lies within the search angle.",
            field="max_angle_deg",
            hint="Widen the search angle, or raise the index limit to admit higher-index axes.",
        )

    current_text = direction_label(current_indices, spec=spec)
    candidates = [row for row in rows if row["verdict"] != "current axis"]
    nearest_reachable = next((row for row in candidates if row["reachable"]), None)
    result = AppResult(
        title=f"Zone axes near {current_text}",
        summary=(
            f"{len(candidates)} zone-axis families of {spec.name} lie within "
            f"{float(request['max_angle_deg']):g}° of {current_text}, of which {reachable_count} "
            f"are reachable within ±{alpha_limit:g}° alpha and ±{beta_limit:g}° beta. "
            + (
                f"The nearest reachable one is {nearest_reachable['target']} at "
                f"{nearest_reachable['angle_deg']:.2f}°, showing "
                f"{nearest_reachable['reflections']} reflections with "
                f"{nearest_reachable['symmetry']} symmetry."
                if nearest_reachable is not None
                else (
                    "None of them is reachable in one move from this position: tilt to the "
                    "closest partial approach, re-index there, and plan the rest from the new "
                    "orientation."
                )
            )
        ),
        table=ResultTable(
            columns=(
                Column("family", "Family"),
                Column(
                    "target",
                    "Nearest member",
                    help_text="The member of the family the holder would actually go to.",
                ),
                Column("angle_deg", "Angle", units="°", numeric=True, digits=2),
                Column(
                    "family_size",
                    "Members",
                    numeric=True,
                    help_text=(
                        "Symmetry-equivalent axes giving the same pattern. More members means "
                        "more chances one of them is inside the holder's range."
                    ),
                ),
                Column(
                    "reflections",
                    "Reflections",
                    numeric=True,
                    help_text=(
                        "Kinematic reflections inside 1.5 Å⁻¹. How much the pattern has to say."
                    ),
                ),
                Column(
                    "symmetry",
                    "Pattern",
                    help_text="Apparent rotational symmetry — what you will recognise on arrival.",
                ),
                Column("verdict", "Reachable"),
                Column("delta_alpha_deg", "Δα", units="°", numeric=True, digits=2),
                Column("delta_beta_deg", "Δβ", units="°", numeric=True, digits=2),
                Column("travel_deg", "Crystal rotation", units="°", numeric=True, digits=2),
                Column(
                    "margin_deg",
                    "Envelope margin",
                    units="°",
                    numeric=True,
                    digits=2,
                    help_text="Tilt range left over. Negative means out of reach.",
                ),
            ),
            rows=tuple(rows),
            caption=f"Zone-axis families of {spec.name} within reach of {current_text}.",
        ),
        data={
            "current_zone_axis": list(current_indices),
            "current_zone_axis_label": current_text,
            "entries": [entry.to_json() for entry in atlas.entries],
            "reachable_count": reachable_count,
            "envelope": {"alpha_limit_deg": alpha_limit, "beta_limit_deg": beta_limit},
            "start": {
                "alpha_deg": float(request["alpha_deg"]),
                "beta_deg": float(request["beta_deg"]),
            },
            "describe": atlas.describe(),
        },
        inputs={
            "phase": spec.to_json(),
            "current_zone_axis": list(current_indices),
            "max_angle_deg": float(request["max_angle_deg"]),
            "alpha_deg": float(request["alpha_deg"]),
            "beta_deg": float(request["beta_deg"]),
            "alpha_limit_deg": alpha_limit,
            "beta_limit_deg": beta_limit,
            "beam_rotation_deg": float(request["beam_rotation_deg"]),
            "max_index": int(request["max_index"]),
            "limit": limit,
        },
        notes=(
            "Reflection counts are kinematic and exclude double diffraction, so a real plate of a "
            "diamond-structure or hcp phase will show a few more spots than the count states.",
            "Angles between zone axes are fixed by the lattice and are independent of the holder, "
            "the rotation about the beam, and the calibration. Only the reachability columns "
            "depend on those.",
        ),
        citations=(_CITATION_WILLIAMS, _CITATION_HIRSCH),
    )
    return result.to_json()


def _member_label(solution: Any, requested: Sequence[int], spec: Any) -> tuple[str, list[int]]:
    """How to name the orbit member a tilt solution places on the beam.

    ``TiltSolution.orbit_member_indices`` is ``None`` whenever the member does
    not rationalize to a low-index triple within the navigation module's bound,
    which happens routinely for a high-index hexagonal family: the member is a
    perfectly good lattice direction, it simply has no tidy integer form at that
    bound. Reading it as an integer array in that case raised, so asking the
    planner for a target such as [4 3̄ 1] in zirconium returned a 500 rather than
    a plan — and the plan itself was fine; only its label was missing.

    Falling back to the *family* form is the honest answer. The move is to some
    member of that family; every number in the row — the tilts, the travel, the
    margin — belongs to the member the planner chose, and only the name of that
    member is unavailable.
    """

    indices = getattr(solution, "orbit_member_indices", None)
    if indices is not None:
        member = [int(value) for value in np.asarray(indices, dtype=int).reshape(-1)]
        return direction_label(tuple(member), spec=spec), member
    requested_indices = [int(value) for value in requested]
    return (
        family_label(tuple(requested_indices), spec=spec, family="direction"),
        requested_indices,
    )


def _picking_scale(request: Mapping[str, Any]) -> float:
    """Picking units per inverse angstrom, for the calibration in ``request``.

    The inverse of what :meth:`PatternCalibration.to_reciprocal_angstrom` does.
    Kept as one function because a calculated pattern drawn at a different scale
    from the measured one it is superimposed on would look like a disagreement
    the crystallography never had.
    """

    units = str(request["units"])
    if units == "reciprocal_angstrom":
        return 1.0
    camera_constant = float(request.get("camera_constant_mm_angstrom") or 0.0)
    if camera_constant <= 0.0:
        return 0.0
    if units == "px":
        pixel_size = float(request.get("pixel_size_mm") or 0.0)
        if pixel_size <= 0.0:
            return 0.0
        return camera_constant / pixel_size
    return camera_constant


def _calculated_overlay(
    spec: Any,
    solution: Any,
    *,
    request: Mapping[str, Any],
    centre: tuple[float, float],
    limit: int = 240,
) -> list[dict[str, Any]]:
    """Where this solution says every reflection of its zone should appear.

    Purpose
    -------
    Accepting an indexing is a judgement, and the honest way to make it is to
    look: draw the pattern the candidate *predicts* on top of the pattern that
    was *measured* and see whether the two coincide. A residual column can hide a
    systematic error that this makes obvious at a glance — a calculated pattern
    uniformly too large is a camera constant, one rotated is a roll, one with
    extra rows is the wrong phase.

    Every reflection of the zone is returned, not only the ones a spot was picked
    for. The unmatched predictions are the informative part: they are where a
    user should look for the spots they have not picked yet, and their absence
    from the plate is evidence against the candidate.

    Notes
    -----
    Positions come from the solution's own orientation, through the same
    projection the indexer used, and are converted to picking coordinates by the
    inverse of the supplied calibration — so the overlay is in the coordinates
    the user clicks in and needs no separate scale.
    """

    from pytex.diffraction.solving import _allowed_reflections

    scale = _picking_scale(request)
    if scale <= 0.0:
        return []
    phase = spec.to_phase()
    hkl, g_crystal = _allowed_reflections(phase, int(request["max_index"]))
    rotation = np.asarray(solution.orientation.as_matrix(), dtype=float)
    projected = g_crystal @ rotation.T
    in_plane = np.linalg.norm(projected[:, :2], axis=1)
    # The zone tolerance is relative to the pattern's own reach, so a plate
    # recorded at a long camera length does not admit reflections a short one
    # would reject.
    reach = float(in_plane.max()) if in_plane.size else 0.0
    on_zone = np.abs(projected[:, 2]) <= max(1e-3 * reach, 1e-9)
    visible = on_zone & (in_plane > 0.0)
    order = np.argsort(np.where(visible, in_plane, np.inf))[: int(limit)]
    overlay: list[dict[str, Any]] = []
    for index in order:
        if not bool(visible[index]):
            continue
        indices = tuple(int(value) for value in hkl[index])
        magnitude = float(in_plane[index])
        overlay.append(
            {
                "hkl": list(indices),
                "label": plane_label(indices, spec=spec),
                "x": float(centre[0] + projected[index, 0] * scale),
                "y": float(centre[1] + projected[index, 1] * scale),
                "g": magnitude,
                "d": 1.0 / magnitude,
            }
        )
    return overlay


def _symmetry_angle_deg(phase: Any, first: Any, second: Any) -> float:
    """Smallest angle between two lattice directions, over the symmetry orbit.

    Purpose
    -------
    Comparing an indexed zone axis with an expected one is not a comparison of
    index triples. A bcc [110] pattern is indistinguishable from a [101] one —
    the crystal symmetry maps one onto the other — so the honest question is how
    far apart the two directions are once every symmetry-equivalent version of
    the first has been tried. That is what this returns.

    Parameters
    ----------
    phase : Phase
    first, second : array_like
        Direction indices in the phase's three-index basis.

    Returns
    -------
    float
        Degrees, in ``[0, 90]``. Zero when the two are the same direction up to
        symmetry and sense.
    """

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    left = direct @ np.asarray(first, dtype=float)
    left = left / float(np.linalg.norm(left))
    right = direct @ np.asarray(second, dtype=float)
    right = right / float(np.linalg.norm(right))
    cosines = np.abs(np.einsum("nij,j->ni", operators, left) @ right)
    return float(math.degrees(math.acos(float(np.clip(cosines.max(), -1.0, 1.0)))))


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
            id="tem.example.gallery_fcc_001",
            title="Practice: index an fcc [001] pattern",
            panel="tem",
            summary="A simulated aluminium plate down [001], with the answer known.",
            teaches=(
                "Start here. The four innermost spots are 200-type, at 90° to each other and all "
                "the same length; the next four are 220-type, at 45° and longer by √2. That "
                "ratio-and-angle signature identifies a cubic ⟨001⟩ zone without using the camera "
                "constant at all — the calibration only enters when you want the lattice "
                "parameter, which is exactly why a wrong camera constant produces a "
                "self-consistent pattern of the wrong material."
            ),
            operation="tem.gallery_pattern",
            request={"pattern": "fcc_al_001"},
        ),
        ExampleScenario(
            id="tem.example.gallery_bcc_110",
            title="Practice: the bcc [110] rectangle",
            panel="tem",
            summary="A simulated ferrite plate down [110], off-centre and rolled.",
            teaches=(
                "The two shortest vectors here are perpendicular but unequal, in the ratio √2 for "
                "any bcc metal — so the rectangle is a lattice-independent check on your beam "
                "centre and your calibration. This is the pattern most often misread as a cubic "
                "⟨001⟩ square, which then indexes to a plausible and entirely wrong lattice "
                "parameter. Note also that the beam is not at the middle of the frame, as it "
                "generally is not on a real plate."
            ),
            operation="tem.gallery_pattern",
            request={"pattern": "bcc_fe_110"},
        ),
        ExampleScenario(
            id="tem.example.gallery_hcp_prism",
            title="Practice: measuring c/a from one hcp pattern",
            panel="tem",
            summary="A simulated zirconium plate down [2̄110], the prism zone.",
            teaches=(
                "The rectangle's aspect ratio here is √3·a/c — 1.088 for zirconium, 1.091 for "
                "titanium, 1.067 for magnesium — so this single pattern measures the axial ratio "
                "and separates the hcp metals from one another with no calibration whatsoever, "
                "because a ratio of two lengths on the same plate does not care what the camera "
                "constant is. Note that 0001 is absent while 0002 is present: the hcp basis "
                "extinguishes odd l on the 000l row, so a 0001 spot on a real plate is double "
                "diffraction rather than a lattice reflection."
            ),
            operation="tem.gallery_pattern",
            request={"pattern": "hcp_zr_2-1-10"},
        ),
        ExampleScenario(
            id="tem.example.atlas_from_bcc_110",
            title="Where can I go from bcc [110]?",
            panel="tem",
            summary="Every zone axis within 60° of ⟨110⟩ in ferrite, ranked and checked.",
            teaches=(
                "The answer is not the axis with the smallest angle. ⟨210⟩ is nearest at 18.43° "
                "but shows twelve reflections; ⟨111⟩ is twice as far at 35.26° and shows "
                "thirty-six, with six-fold symmetry that is unmistakable the moment it arrives. "
                "The cost of a tilt is a few minutes and the risk of losing the grain; the value "
                "is the information the new pattern carries, and this table is what lets you "
                "weigh one against the other before touching the controls. Raise the index limit "
                "to 3 and nearer families appear — ⟨320⟩ at 11.31° with eight reflections — which "
                "is the trade-off made explicit."
            ),
            operation="tem.zone_axis_atlas",
            request={
                "phase": {"builtin": "fe_bcc"},
                "current_zone_axis": [1, 1, 0],
                "alpha_limit_deg": 30.0,
                "beta_limit_deg": 20.0,
            },
        ),
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
