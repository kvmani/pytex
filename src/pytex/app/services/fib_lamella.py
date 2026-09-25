"""FIB lamella planning in the workbench: which grain, where, at what azimuth.

One operation, ``fib.lamella_plan``, in the EBSD workspace. It reads the scan the
workspace has open (or a practice map), segments it, and asks
:func:`pytex.fib.rank_grains` for the grains that make the best lamella for a
target TEM zone axis --- or, when a grain is named, plans that grain alone. The
answer is the ranked table, the seven figures of the foundation document, the
prose of :meth:`pytex.fib.LamellaPlanReport.describe`, and a printable one-page
work order for the operator at the FIB.

The form follows the input groups of section 10 of
``docs/architecture/fib_lamella_planning_foundation.md``: file and phase, target,
frames and registration (decisions D2 and D3), chamber (D4), lamella dimensions,
holder envelope (D7), uncertainty inputs, ranking weights. Everything mandatory
has a default, so the button works on first press; everything instrument-
specific is a declared input rather than an assumption.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.ebsd_gallery import build_map, entry_ids, get_entry
from pytex.app.errors import InvalidInputError
from pytex.app.figures import render_figure
from pytex.app.logbook import APP_LOG
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    DocumentationLink,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    TextParameter,
)
from pytex.app.results import (
    AppResult,
    Column,
    ResultFigure,
    ResultMetric,
    ResultStage,
    ResultTable,
)
from pytex.app.services.ebsd import _imported_map, _source_parameters
from pytex.core.symbols import symbol_text

__all__: tuple[str, ...] = ()

_THETA_S = symbol_text("lamella_azimuth")
_THETA_I = symbol_text("lamella_image_azimuth")
_THETA_ION = symbol_text("ion_azimuth")

_DOCUMENTATION = DocumentationLink("FIB lamella planning", "workflows/fib_lamella_planning")

_CITATION_ITA = (
    "International Tables for Crystallography, Vol. A (point-group orbits), "
    "doi:10.1107/97809553602060000100."
)
_CITATION_DE_GRAEF = (
    "De Graef, Introduction to Conventional Transmission Electron Microscopy, CUP 2003, "
    "doi:10.1017/CBO9780511615092 (double-tilt geometry)."
)
_CITATION_GIANNUZZI = (
    "Giannuzzi & Stevie, A review of focused ion beam milling techniques for TEM specimen "
    "preparation, Micron 30 (1999) 197, doi:10.1016/S0968-4328(99)00005-0."
)
_CITATION_GUM = "JCGM 100:2008 (GUM), evaluation of measurement data, uncertainty in quadrature."

_CANDIDATE_COLUMNS = (
    Column("rank", "Rank", numeric=True),
    Column("grain_id", "Grain", numeric=True),
    Column(
        "member",
        "Axis on beam",
        help_text="The symmetry-equivalent of the target that lies closest to the surface plane.",
    ),
    Column(
        "eps_deg",
        "Residual tilt",
        units="°",
        numeric=True,
        digits=2,
        help_text="eps*: how far the chosen axis lies out of the surface, the tilt the holder "
        "must supply.",
    ),
    Column(
        "expanded_deg",
        "± U",
        units="°",
        numeric=True,
        digits=2,
        help_text="Expanded uncertainty k u(eps*) of the residual tilt.",
    ),
    Column("theta_sample_deg", f"Azimuth {_THETA_S}", units="°", numeric=True, digits=2),
    Column("theta_image_deg", f"Azimuth {_THETA_I}", units="°", numeric=True, digits=2),
    Column(
        "theta_ion_deg",
        f"FIB rotation {_THETA_ION}",
        units="°",
        numeric=True,
        digits=2,
        help_text="Pattern rotation for the lamella normal. Only as good as the chamber "
        "calibration; see the warnings.",
    ),
    Column("feasibility", "Feasibility"),
    Column(
        "mount_fraction",
        "Mounts that work",
        numeric=True,
        digits=2,
        help_text="Fraction of the unknown in-plane mounting rotations for which the holder "
        "reaches the axis.",
    ),
    Column("fits", "Fits in grain"),
    Column("margin_um", "Clearance", units="µm", numeric=True, digits=2),
    Column("gos_deg", "GOS", units="°", numeric=True, digits=2),
    Column("score", "Score", numeric=True, digits=3),
    Column("x_um", "Site x", units="µm", numeric=True, digits=2),
    Column("y_um", "Site y", units="µm", numeric=True, digits=2),
)


def _parameters() -> tuple[Any, ...]:
    return (
        *_source_parameters(),
        NumberParameter(
            name="practice_step_um",
            label="Practice-map step",
            help_text=(
                "Step of the practice maps, in micrometres. The EBSD tabs draw them at 0.5 µm, "
                "which makes grains too small for a 15 µm lamella; planning draws the same "
                "construction coarser so real lamellae fit. Ignored for a scan of your own."
            ),
            units="µm",
            default=1.5,
            minimum=0.1,
            maximum=20.0,
            advanced=True,
        ),
        NumberParameter(
            name="grain_threshold_deg",
            label="Grain threshold",
            help_text=(
                "Misorientation above which neighbouring points belong to different grains. "
                "The lamella must sit inside one grain, so this decides where the grain "
                "boundaries the rectangle keeps clear of are."
            ),
            units="°",
            default=5.0,
            minimum=0.5,
            maximum=30.0,
            group="Scan and phase",
        ),
        TextParameter(
            name="phase_name",
            label="Phase",
            help_text=(
                "Which phase to cut, by the name the scan file gives it. Leave empty for the "
                "scan's primary phase; points of other phases are left out of the ranking."
            ),
            required=False,
            max_length=80,
            group="Scan and phase",
        ),
        IndicesParameter(
            name="zone_axis",
            label="Target zone axis [uvw]",
            help_text=(
                "The TEM zone axis you want to look down. Any symmetry equivalent, in either "
                "sense, is accepted (decision D8), so [0 1 1] and [1 -1 0] of a cubic crystal "
                "are the same request."
            ),
            default=[0, 1, 1],
            group="Target",
        ),
        IntegerParameter(
            name="grain_id",
            label="Plan one grain",
            help_text=(
                "Grain id to plan on its own, as the EBSD grain table numbers it. Leave empty "
                "to rank every grain of the map and recommend the best."
            ),
            required=False,
            minimum=0,
            maximum=1_000_000,
            group="Target",
        ),
        BooleanParameter(
            name="both_senses",
            label="Either beam sense is acceptable",
            help_text=(
                "On (the default) accepts the target and its reverse, which a kinematic pattern "
                "cannot tell apart. Turn it off only for CBED or polarity work where the sense "
                "has been determined; the lamella's front/back flip then decides success."
            ),
            default=True,
            group="Target",
            advanced=True,
        ),
        ChoiceParameter(
            name="surface_normal",
            label="Outward surface normal",
            help_text=(
                "Which sense of the EBSD specimen z axis points out of the polished surface. "
                "Vendors differ and no file header says; a wrong choice mirrors every azimuth. "
                "Check it once per instrument against a fiducial."
            ),
            options=(
                ("+z", "+z of the scan's specimen frame", "The specimen z axis points out."),
                ("-z", "-z of the scan's specimen frame", "The specimen z axis points in."),
            ),
            default="+z",
            group="Frames and registration",
        ),
        ChoiceParameter(
            name="scan_rows",
            label="Scan rows run along",
            help_text=(
                "Which sense of the specimen y axis the scan's y coordinate increases along."
            ),
            options=(
                ("+y", "+y of the specimen frame", "Scan y increases along specimen +y."),
                ("-y", "-y of the specimen frame", "Scan y increases along specimen -y."),
            ),
            default="+y",
            group="Frames and registration",
        ),
        ChoiceParameter(
            name="registration",
            label="SEM image registration",
            help_text=(
                "How EBSD scan coordinates map onto the SEM image the rectangle is milled "
                "against. *Identity* is decision D3: image axes along the scan axes, one unit "
                "per micrometre. *Similarity* applies the rotation, scale and flip below. "
                "*Affine* fits the control points below and reports their residual, which "
                "enters the uncertainty budget."
            ),
            options=(
                ("identity", "Identity (scan = image)", "The default of decision D3."),
                ("similarity", "Rotation, scale, flip", "A declared similarity transform."),
                ("affine", "Affine from control points", "Least squares on matched features."),
            ),
            default="identity",
            group="Frames and registration",
        ),
        NumberParameter(
            name="image_rotation_deg",
            label="Image rotation",
            help_text="Rotation of the scan axes into the image axes, for a similarity.",
            units="°",
            default=0.0,
            minimum=-360.0,
            maximum=360.0,
            group="Frames and registration",
            row="similarity",
        ),
        NumberParameter(
            name="image_scale",
            label="Image units per µm",
            help_text="Scale of the similarity: image pixels (or micrometres) per scan micrometre.",
            default=1.0,
            minimum=1e-6,
            maximum=1e6,
            exclusive_minimum=True,
            field_width="short",
            group="Frames and registration",
            row="similarity",
        ),
        BooleanParameter(
            name="image_flip",
            label="Mirror the image (flip y)",
            help_text="The image is a mirror of the scan: reverse scan y before rotating.",
            default=False,
            group="Frames and registration",
        ),
        TextParameter(
            name="control_points",
            label="Control points",
            help_text=(
                "For the affine registration: one feature per line, 'x_scan y_scan x_image "
                "y_image' in scan micrometres and image units. Three points fit exactly; four "
                "or more also measure the misfit."
            ),
            required=False,
            multiline=True,
            placeholder="0 0 12 8\n60 0 132 14\n0 40 6 88\n60 40 128 94",
            group="Frames and registration",
            advanced=True,
        ),
        ChoiceParameter(
            name="column_angle",
            label="Ion-column angle",
            help_text=(
                "Angle between the electron and ion columns (decision D4). The stage is tilted "
                "to it for milling, so the ion beam is normal to the surface."
            ),
            options=(
                ("54", "54°", "The default."),
                ("52", "52°", "The other common column geometry."),
            ),
            default="54",
            symbol="column_angle",
            field_width="short",
            group="Chamber",
        ),
        ChoiceParameter(
            name="rotation_sense",
            label="Rotation sense",
            help_text=(
                "s_R: +1 when a positive FIB pattern rotation turns the box the same way as a "
                "positive sample azimuth. Measured by the fiducial calibration; never assumed."
            ),
            options=(("+1", "+1", "Same sense."), ("-1", "-1", "Opposite sense.")),
            default="+1",
            symbol="rotation_sense",
            field_width="tiny",
            group="Chamber",
            row="calibration",
        ),
        NumberParameter(
            name="rotation_offset_deg",
            label="Rotation offset",
            help_text="R0: the pattern rotation that aligns the box with sample azimuth zero.",
            units="°",
            default=0.0,
            minimum=-360.0,
            maximum=360.0,
            symbol="rotation_offset",
            group="Chamber",
            row="calibration",
        ),
        NumberParameter(
            name="stage_rotation_deg",
            label="Stage rotation",
            help_text="The stage rotation reading at which the lamella will be milled.",
            units="°",
            default=0.0,
            minimum=-360.0,
            maximum=360.0,
            symbol="stage_rotation",
            group="Chamber",
            row="calibration",
        ),
        BooleanParameter(
            name="chamber_calibrated",
            label="Sense and offset come from the fiducial calibration",
            help_text=(
                "Tick only after running the calibration protocol on this instrument. Until "
                "then every plan and the work order carry an UNCALIBRATED warning on the FIB "
                "rotation, because a wrong sign destroys a real lamella."
            ),
            default=False,
            group="Chamber",
        ),
        NumberParameter(
            name="length_um",
            label="Length",
            help_text="L, along the lamella's long axis.",
            units="µm",
            default=15.0,
            minimum=1.0,
            maximum=200.0,
            symbol="lamella_length",
            group="Lamella",
            row="lamella",
        ),
        NumberParameter(
            name="width_um",
            label="Width",
            help_text="W, trench to trench across the lamella, including the protective strap.",
            units="µm",
            default=2.0,
            minimum=0.2,
            maximum=50.0,
            symbol="lamella_width",
            group="Lamella",
            row="lamella",
        ),
        NumberParameter(
            name="depth_um",
            label="Depth",
            help_text="D, the milling depth below the surface.",
            units="µm",
            default=8.0,
            minimum=0.5,
            maximum=100.0,
            symbol="mill_depth",
            group="Lamella",
            row="lamella",
        ),
        NumberParameter(
            name="thickness_nm",
            label="Final thickness",
            help_text="t, the electron-transparent thickness after thinning.",
            units="nm",
            default=80.0,
            minimum=5.0,
            maximum=1000.0,
            symbol="foil_thickness",
            group="Lamella",
            row="lamella",
        ),
        NumberParameter(
            name="alpha_limit_deg",
            label="Alpha limit",
            help_text="The holder's alpha range is ± this (decision D7).",
            units="°",
            default=30.0,
            minimum=1.0,
            maximum=80.0,
            symbol="alpha_tilt",
            group="TEM holder",
            row="envelope",
        ),
        NumberParameter(
            name="beta_limit_deg",
            label="Beta limit",
            help_text="The holder's beta range is ± this (decision D7).",
            units="°",
            default=30.0,
            minimum=1.0,
            maximum=80.0,
            symbol="beta_tilt",
            group="TEM holder",
            row="envelope",
        ),
        NumberParameter(
            name="margin_deg",
            label="Safety margin",
            help_text=(
                "How far inside the envelope the residual must stay, at its upper bound, for "
                "the plan to count as guaranteed for every mounting rotation."
            ),
            units="°",
            default=5.0,
            minimum=0.0,
            maximum=20.0,
            group="TEM holder",
            row="envelope",
        ),
        NumberParameter(
            name="ebsd_accuracy_deg",
            label="EBSD accuracy",
            help_text="Standard uncertainty of an indexed orientation.",
            units="°",
            default=0.5,
            minimum=0.0,
            maximum=10.0,
            group="Uncertainty",
            row="uncertainty",
        ),
        NumberParameter(
            name="mount_repeatability_deg",
            label="Mount repeatability",
            help_text="Standard uncertainty added by lift-out, welding and seating the grid.",
            units="°",
            default=2.0,
            minimum=0.0,
            maximum=20.0,
            group="Uncertainty",
            row="uncertainty",
        ),
        NumberParameter(
            name="coverage_factor",
            label="Coverage factor",
            help_text="k of the expanded uncertainty U = k u; 2 is about 95 %.",
            default=2.0,
            minimum=0.5,
            maximum=5.0,
            symbol="coverage_factor",
            group="Uncertainty",
            row="uncertainty",
        ),
        NumberParameter(
            name="weight_eps",
            label="Residual weight",
            help_text="Weight of the residual-tilt term, the primary one.",
            default=0.5,
            minimum=0.0,
            maximum=10.0,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_a",
        ),
        NumberParameter(
            name="weight_fit",
            label="Fit weight",
            help_text="Weight of the clearance between the rectangle and the grain boundary.",
            default=0.2,
            minimum=0.0,
            maximum=10.0,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_a",
        ),
        NumberParameter(
            name="weight_reliability",
            label="Reliability weight",
            help_text="Weight of exp(-GOS), how well one orientation describes the grain.",
            default=0.1,
            minimum=0.0,
            maximum=10.0,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_a",
        ),
        NumberParameter(
            name="weight_edge",
            label="Edge weight",
            help_text="Weight of the site's distance from the edge of the scan.",
            default=0.1,
            minimum=0.0,
            maximum=10.0,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_b",
        ),
        NumberParameter(
            name="weight_neighbour",
            label="Neighbour weight",
            help_text="Weight of the grain lying wholly inside the scan rather than cut by it.",
            default=0.1,
            minimum=0.0,
            maximum=10.0,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_b",
        ),
        IntegerParameter(
            name="shortlist",
            label="Grains planned in full",
            help_text=(
                "Every grain is screened; this many of the best receive the full mounting "
                "sweep and footprint fit."
            ),
            default=6,
            minimum=1,
            maximum=40,
            group="Ranking weights",
            group_collapsed=True,
            row="weights_b",
        ),
    )


def _crystal_map(request: dict[str, Any]) -> tuple[Any, str]:
    scan_file = request.get("scan_file")
    if scan_file:
        crystal_map, entry = _imported_map(scan_file)
        return crystal_map, str(entry.title)
    entry_id = str(request["dataset"])
    try:
        entry = get_entry(entry_id)
    except KeyError as error:
        raise InvalidInputError(
            f"Unknown dataset {entry_id!r}.",
            field="dataset",
            hint=f"Choose one of: {', '.join(entry_ids())}.",
        ) from error
    return (
        build_map(
            entry_id,
            grid=int(request["grid_points"]),
            step_um=float(request["practice_step_um"]),
        ),
        f"{entry.title} (practice map)",
    )


def _phase(crystal_map: Any, name: str | None) -> Any:
    if not name:
        phase = crystal_map.primary_phase
        if phase is None:
            raise InvalidInputError(
                "The scan declares no phase to plan for.",
                field="phase_name",
                hint="Name the phase, or open a scan whose header declares one.",
            )
        return phase
    candidates = [entry.phase for entry in crystal_map.resolved_phase_entries]
    if crystal_map.primary_phase is not None:
        candidates.append(crystal_map.primary_phase)
    for phase in candidates:
        if phase is not None and str(phase.name).strip().lower() == name.strip().lower():
            return phase
    known = ", ".join(sorted({str(phase.name) for phase in candidates if phase is not None}))
    raise InvalidInputError(
        f"The scan has no phase called {name!r}.",
        field="phase_name",
        hint=f"Phases in this scan: {known or 'none declared'}.",
    )


def _registration(request: dict[str, Any]) -> Any:
    from pytex.fib import ImageRegistration

    method = request["registration"]
    if method == "identity":
        return ImageRegistration.identity()
    if method == "similarity":
        return ImageRegistration.similarity(
            rotation_deg=float(request["image_rotation_deg"]),
            scale=float(request["image_scale"]),
            flip_y=bool(request["image_flip"]),
        )
    text = str(request.get("control_points") or "").strip()
    rows: list[list[float]] = []
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        try:
            values = [float(item) for item in line.replace(",", " ").split()]
        except ValueError as error:
            raise InvalidInputError(
                f"Control point line {number} is not four numbers: {line!r}.",
                field="control_points",
            ) from error
        if len(values) != 4:
            raise InvalidInputError(
                f"Control point line {number} needs four numbers, got {len(values)}.",
                field="control_points",
                hint="Write x_scan y_scan x_image y_image on each line.",
            )
        rows.append(values)
    if len(rows) < 3:
        raise InvalidInputError(
            "The affine registration needs at least three control points.",
            field="control_points",
            hint="Locate three or more features in both the EBSD map and the SEM image.",
        )
    array = np.asarray(rows)
    try:
        return ImageRegistration.from_control_points(array[:, :2], array[:, 2:])
    except ValueError as error:
        raise InvalidInputError(str(error), field="control_points") from error


def _settings(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.fib import (
        ChamberGeometry,
        LamellaRankingWeights,
        LamellaSpec,
        MountModel,
        SurfaceGeometry,
        UncertaintyInputs,
    )
    from pytex.tem.stage import DoubleTiltStage, RectangularEnvelope

    alpha, beta = float(request["alpha_limit_deg"]), float(request["beta_limit_deg"])
    try:
        spec = LamellaSpec(
            length_um=float(request["length_um"]),
            width_um=float(request["width_um"]),
            depth_um=float(request["depth_um"]),
            final_thickness_nm=float(request["thickness_nm"]),
        )
        weights = LamellaRankingWeights(
            eps=float(request["weight_eps"]),
            fit=float(request["weight_fit"]),
            reliability=float(request["weight_reliability"]),
            edge=float(request["weight_edge"]),
            neighbour=float(request["weight_neighbour"]),
        )
    except ValueError as error:
        raise InvalidInputError(str(error), field="length_um") from error
    return {
        "surface": SurfaceGeometry(
            normal_sign=1 if request["surface_normal"] == "+z" else -1,
            scan_y_sign=1 if request["scan_rows"] == "+y" else -1,
        ),
        "registration": _registration(request),
        "chamber": ChamberGeometry(
            column_angle_deg=float(request["column_angle"]),
            rotation_sense=1 if request["rotation_sense"] == "+1" else -1,
            rotation_offset_deg=float(request["rotation_offset_deg"]),
            stage_rotation_deg=float(request["stage_rotation_deg"]),
            calibrated=bool(request["chamber_calibrated"]),
        ),
        "spec": spec,
        "stage": DoubleTiltStage(
            envelope=RectangularEnvelope(-alpha, alpha, -beta, beta),
            name=f"double tilt, +/-{alpha:g} deg x +/-{beta:g} deg",
        ),
        "mount": MountModel(margin_deg=float(request["margin_deg"]), solver="closed_form"),
        "uncertainty": UncertaintyInputs(
            ebsd_accuracy_deg=float(request["ebsd_accuracy_deg"]),
            mount_repeatability_deg=float(request["mount_repeatability_deg"]),
            coverage_factor=float(request["coverage_factor"]),
        ),
        "weights": weights,
    }


def _single_grain(segmentation: Any, grain_id: int, report: Any) -> Any:
    """Restrict a ranking to one named grain, or explain why it cannot be."""

    for plan in report.plans:
        if plan.grain_id == grain_id:
            return plan
    raise InvalidInputError(
        f"Grain {grain_id} was not planned.",
        field="grain_id",
        hint=(
            f"The segmentation has {len(segmentation.grains)} grains numbered from 0; grains "
            "smaller than ten points or of another phase are not planned."
        ),
    )


@REGISTRY.operation(
    "fib.lamella_plan",
    title="FIB lamella for a target zone axis",
    summary=(
        "Which grain to cut, where, at what azimuth in every frame, and what TEM tilt remains."
    ),
    help_text=(
        "Plans a site-specific TEM lamella from the open EBSD scan so that the finished "
        "lamella shows a chosen zone axis.\n\n"
        "A lamella milled **strictly vertically** has its faces perpendicular to the surface, "
        "so its normal - the TEM beam direction - lies in the surface plane. A grain is usable "
        "for a zone axis exactly when some symmetry equivalent of it lies close to that plane, "
        "and the angle by which it misses, **eps***, is the tilt the TEM holder must supply. "
        "Every grain is screened for that angle; the best are then planned in full: the lamella "
        "azimuth in the sample frame, in the SEM image and as the FIB pattern rotation, the "
        "rectangle placed as far inside the grain as it goes, and the holder tilts for every "
        "possible way the lamella can land on the grid.\n\n"
        "**Read the feasibility class, not only eps*.** The in-plane rotation of a lamella on "
        "the TEM grid is not controlled, so a plan is *guaranteed* only when the residual, at "
        "its upper uncertainty bound and with a safety margin, is reachable for every mounting "
        "rotation; *probabilistic* plans work for the stated fraction of them.\n\n"
        "**The FIB rotation is only as good as its calibration.** Its sign and offset are "
        "instrument-specific: run the fiducial protocol once and tick *Sense and offset come "
        "from the fiducial calibration*, or every result carries an UNCALIBRATED warning.\n\n"
        "The grain is assumed columnar below the surface, and the chain has not yet been "
        "validated against a lamella cut on a real instrument; both are stated on every result. "
        "The *Print work order* button gives the one page to take to the FIB."
    ),
    parameters=_parameters(),
    returns=(
        "The ranked candidates, the recommended plan with its figures (plan view, stereogram, "
        "mounting-rotation feasibility, section, preparability map, predicted SAED), the "
        "uncertainty budget, the prose explanation and a printable work order."
    ),
    panel="fib_lamella",
    citations=(_CITATION_ITA, _CITATION_DE_GRAEF, _CITATION_GIANNUZZI, _CITATION_GUM),
    tags=("FIB", "lamella", "lift-out", "TEM preparation", "zone axis", "site-specific"),
    documentation=_DOCUMENTATION,
)
def _lamella_plan(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.fib import rank_grains, work_order_html

    crystal_map, source = _crystal_map(request)
    phase = _phase(crystal_map, request.get("phase_name"))
    settings = _settings(request)
    target = tuple(int(value) for value in request["zone_axis"])
    try:
        segmentation = crystal_map.segment_grains(
            max_misorientation_deg=float(request["grain_threshold_deg"])
        )
        grain_id = request.get("grain_id")
        report = rank_grains(
            segmentation,
            target,
            phase=phase,
            surface=settings["surface"],
            registration=settings["registration"],
            chamber=settings["chamber"],
            spec=settings["spec"],
            stage=settings["stage"],
            mount=settings["mount"],
            uncertainty=settings["uncertainty"],
            weights=settings["weights"],
            both_senses=bool(request["both_senses"]),
            shortlist=(
                len(segmentation.grains) if grain_id is not None else int(request["shortlist"])
            ),
        )
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            field="zone_axis",
            hint="Check the zone axis, the phase and that the scan is a square grid with a step.",
        ) from error
    if not report.plans:
        raise InvalidInputError(
            " ".join(report.notes) or "No grain could be planned.",
            field="grain_threshold_deg",
            hint="Lower the grain threshold, or open a scan with larger grains.",
        )
    if grain_id is not None:
        best = _single_grain(segmentation, int(grain_id), report)
        plans = (best,)
    else:
        best = report.best()
        plans = report.plans
    APP_LOG.info(
        f"Planned {len(report.plans)} lamellae for {report.orbit.family_text}.",
        source="fib.lamella_plan",
        detail={"grains": report.grain_count, "best": best.grain_id},
    )

    rows, cols = crystal_map.grid_shape
    background_label = "image quality"
    if "image_quality" in crystal_map.property_names:
        background = np.asarray(crystal_map.get_property("image_quality"), dtype=float)
    elif "band_contrast" in crystal_map.property_names:
        background = np.asarray(crystal_map.get_property("band_contrast"), dtype=float)
        background_label = "band contrast"
    else:
        background = np.asarray(segmentation.labels, dtype=float) % 7.0
        background_label = "grain"
    background = background.reshape(rows, cols)
    mask = segmentation.label_grid == best.grain_id
    figures = _figures(best, report, background, mask, background_label, settings)
    order_html = work_order_html(
        best,
        figures=tuple(
            (figure.caption, figure.svg)
            for figure in figures.values()
            if figure.key in {"fib_plan_view", "fib_section"}
        ),
    )
    result = _result(
        request=request,
        source=source,
        report=report,
        best=best,
        plans=plans,
        figures=figures,
        order_html=order_html,
    )
    return result.to_json()


def _figures(
    plan: Any,
    report: Any,
    background: np.ndarray,
    mask: np.ndarray,
    background_label: str,
    settings: dict[str, Any],
) -> dict[str, ResultFigure]:
    from pytex.plotting.fib_figures import (
        draw_phi_feasibility,
        draw_plan_view,
        draw_predicted_saed,
        draw_preparability_map,
        draw_section_schematic,
        draw_stereogram,
    )

    extent = report.raster.extent_um
    figures = {
        "fib_plan_view": render_figure(
            lambda figure: draw_plan_view(
                figure,
                plan,
                background=background,
                extent=extent,
                grain_mask=mask,
                background_label=background_label,
            ),
            key="fib_plan_view",
            title="Plan view of the site",
            caption=(
                f"The {background_label} map around the grain (outlined), with the lamella "
                "rectangle, its normal and azimuths, a scale bar and the sample axes. Scan y runs "
                "downwards, as the image is displayed; the diamond marks the +t_L end, a "
                "suggested lift-out end."
            ),
            interpretation=plan.placement.describe() if plan.placement is not None else None,
            height_in=4.4,
        ),
        "fib_stereogram": render_figure(
            lambda figure: draw_stereogram(
                figure, plan, guaranteed_radius_deg=report.guaranteed_radius_deg
            ),
            key="fib_stereogram",
            title="Why this azimuth: the orbit against the surface plane",
            caption=(
                "Upper-hemisphere stereogram of the sample frame, Z_s at the centre. The outer "
                "circle is the surface plane, the only directions a vertical lamella can put on "
                "the beam; every orbit member is a dot, the chosen one ringed, its distance from "
                "the circle is eps*. The shaded band is reachable for every mounting rotation."
            ),
            interpretation=(
                f"The chosen member lies {plan.eps_deg:.2f} degrees from the surface plane; the "
                f"band reaches {report.guaranteed_radius_deg:.1f} degrees."
            ),
            width_in=5.2,
            height_in=5.2,
        ),
        "fib_phi_feasibility": render_figure(
            lambda figure: draw_phi_feasibility(figure, plan, envelope=settings["stage"].envelope),
            key="fib_phi_feasibility",
            title="Feasibility over the unknown mounting rotation",
            caption=(
                "Polar plot of the in-plane mounting rotation phi. Shaded: the largest residual "
                "the holder removes at each phi (exact solution); dashed: the guaranteed disc; "
                "solid and dotted circles: eps* and its upper bound; thick arcs: the rotations "
                "for which the front and back branches reach the axis (back drawn just outside)."
            ),
            interpretation=(
                f"Reachable for {100 * plan.sweep.fraction:.0f}% of mounting rotations; class "
                f"{plan.feasibility.value}."
            ),
            height_in=4.2,
        ),
        "fib_section": render_figure(
            lambda figure: draw_section_schematic(figure, plan),
            key="fib_section",
            title="Section across the lamella",
            caption=(
                "Viewed along the long axis: the two trenches, the slab of width W and depth D, "
                "its thinned window, the ion beam milling along -Z_s, and the target axis at eps* "
                "from the lamella normal, in the plane of the drawing."
            ),
            height_in=3.8,
        ),
        "fib_preparability": render_figure(
            lambda figure: draw_preparability_map(
                figure,
                report.raster,
                plans=report.plans,
                guaranteed_radius_deg=report.guaranteed_radius_deg,
            ),
            key="fib_preparability",
            title="Preparability map",
            caption=(
                "eps* at every point of the map (viridis, low is good), with the dashed contour "
                "at the guaranteed radius and the planned rectangles numbered by rank."
            ),
            height_in=4.4,
        ),
        "fib_saed": render_figure(
            lambda figure: draw_predicted_saed(figure, plan),
            key="fib_saed",
            title="Predicted diffraction pattern",
            caption=(
                "Kinematic SAED pattern down the chosen member, the axis the lamella shows once "
                "the holder has removed the residual tilt."
            ),
            width_in=4.8,
            height_in=4.4,
        ),
    }
    return figures


def _number(value: float | None) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return float(value)


def _candidate_rows(plans: tuple[Any, ...]) -> tuple[dict[str, Any], ...]:
    rows = []
    for rank, plan in enumerate(plans, start=1):
        placement = plan.placement
        centre = placement.center_scan_um if placement is not None else plan.location_um
        rows.append(
            {
                "rank": rank,
                "grain_id": plan.grain_id,
                "member": plan.member_text,
                "eps_deg": plan.eps_deg,
                "expanded_deg": plan.budget.expanded_deg,
                "theta_sample_deg": plan.geometry.theta_sample_deg,
                "theta_image_deg": plan.theta_image_deg,
                "theta_ion_deg": plan.theta_ion_deg,
                "feasibility": plan.feasibility.value,
                "mount_fraction": plan.sweep.fraction,
                "fits": "yes" if placement is not None and placement.fits else "no",
                "margin_um": placement.margin_um if placement is not None else None,
                "gos_deg": plan.grain_spread_deg,
                "score": plan.score,
                "x_um": None if centre is None else float(centre[0]),
                "y_um": None if centre is None else float(centre[1]),
            }
        )
    return tuple(rows)


def _result(
    *,
    request: dict[str, Any],
    source: str,
    report: Any,
    best: Any,
    plans: tuple[Any, ...],
    figures: dict[str, ResultFigure],
    order_html: str,
) -> AppResult:
    placement = best.placement
    budget = best.budget
    grain_text = f"grain {best.grain_id}" if best.grain_id is not None else "the site"
    summary = (
        f"Cut {grain_text} of {best.phase_name} for {report.orbit.family_text}: the member "
        f"{best.member_text} lies eps* = {best.eps_deg:.2f} ± {budget.expanded_deg:.2f}° out of "
        f"the surface, so the TEM holder must tilt that far. Mill the lamella with its normal at "
        f"{best.geometry.theta_sample_deg:.2f}° in the sample frame "
        f"({best.theta_image_deg:.2f}° in the SEM image, FIB pattern rotation "
        f"{best.theta_ion_deg:.2f}°). Feasibility: {best.feasibility.value}, reachable for "
        f"{100 * best.sweep.fraction:.0f}% of the unknown mounting rotations."
    )
    if placement is not None:
        summary += f" The {best.spec.length_um:g} x {best.spec.width_um:g} µm rectangle " + (
            f"fits with {placement.margin_um:.2f} µm clearance."
            if placement.fits
            else f"does not fit; the longest that fits is {placement.largest_length_um:.1f} µm."
        )
    summary += f" Scan: {source}."
    highlights = [
        ResultMetric("Recommended grain", best.grain_id if best.grain_id is not None else "-"),
        ResultMetric(
            "Residual tilt eps*",
            f"{best.eps_deg:.2f} ± {budget.expanded_deg:.2f}",
            "°",
            f"Expanded uncertainty with k = {budget.coverage_factor:g}; the verdict uses the "
            f"upper bound {budget.eps_upper_deg:.2f}°.",
        ),
        ResultMetric("Axis on the beam", best.member_text),
        ResultMetric(
            f"Lamella normal azimuth {_THETA_S}", round(best.geometry.theta_sample_deg, 2), "°"
        ),
        ResultMetric(f"SEM image azimuth {_THETA_I}", round(best.theta_image_deg, 2), "°"),
        ResultMetric(
            f"FIB pattern rotation {_THETA_ION}",
            round(best.theta_ion_deg, 2),
            "°",
            best.chamber.calibration_caveat(),
        ),
        ResultMetric(
            "Feasibility",
            best.feasibility.value,
            None,
            f"Reachable for {100 * best.sweep.fraction:.0f}% of mounting rotations.",
        ),
    ]
    if placement is not None:
        highlights.append(
            ResultMetric(
                "Clearance to the grain boundary",
                round(placement.margin_um, 2) if placement.fits else "does not fit",
                "µm" if placement.fits else None,
            )
        )
    warnings = [flag.message for flag in best.risk_flags if flag.severity != "info"]
    notes = [flag.message for flag in best.risk_flags if flag.severity == "info"]
    notes.extend(report.notes)

    budget_rows = (
        *(
            {"term": item.name, "value_deg": item.value_deg, "source": item.source}
            for item in budget.components
        ),
        {"term": "combined u(eps*)", "value_deg": budget.combined_deg, "source": "quadrature"},
        {
            "term": f"expanded U (k = {budget.coverage_factor:g})",
            "value_deg": budget.expanded_deg,
            "source": "k u",
        },
    )
    work_rows = tuple({"item": item, "value": value} for item, value in best.work_order_lines())
    option_rows = tuple(
        {
            "member": "[" + " ".join(str(value) for value in option.member_indices) + "]",
            "eps_deg": option.eps_deg,
            "theta_sample_deg": option.theta_sample_deg,
        }
        for option in best.options
    )
    stages = (
        ResultStage(
            key="recommended",
            title="1. The recommended lamella",
            summary=(
                f"{best.member_text} at eps* = {best.eps_deg:.2f}°; normal at "
                f"{_THETA_S} = {best.geometry.theta_sample_deg:.2f}°."
            ),
            metrics=tuple(highlights[1:6]),
            explanation=(
                "The plan view shows where to mill; the stereogram shows why this azimuth: of "
                "every symmetry equivalent of the target, the ringed one is nearest the circle "
                "of directions a vertical lamella can put on the beam."
            ),
            section="result",
            figures=(figures["fib_plan_view"], figures["fib_stereogram"]),
        ),
        ResultStage(
            key="mounting",
            title="2. The unknown mounting rotation",
            summary=(
                f"{best.feasibility.value}: front branch {100 * best.sweep.fraction_front:.0f}%, "
                f"back branch {100 * best.sweep.fraction_back:.0f}% of rotations; "
                f"{100 * best.sweep.conservative_fraction:.0f}% at the upper bound."
            ),
            metrics=(
                ResultMetric("Guaranteed radius", round(report.guaranteed_radius_deg, 2), "°"),
                ResultMetric(
                    "Solver", "TEM navigation" if best.sweep.solver == "navigation" else "exact"
                ),
            ),
            explanation=(
                "A lamella's in-plane rotation on the grid and its front/back flip are unknown "
                "until it is in the TEM, so the residual is a known angle whose split into alpha "
                "and beta is not. The plan therefore never quotes one (alpha, beta) pair: it "
                "reports for which rotations the holder can close the gap."
            ),
            section="result",
            figures=(figures["fib_phi_feasibility"],),
        ),
        ResultStage(
            key="candidates",
            title="3. The other candidates",
            summary=(
                f"{report.candidate_count} of {report.grain_count} grains screened; "
                f"{len(report.plans)} planned in full."
            ),
            explanation=report.weights.describe(),
            section="evidence",
            figures=(figures["fib_preparability"],),
        ),
        ResultStage(
            key="section",
            title="4. The lamella in section, and what it will show",
            summary=best.spec.describe(),
            explanation=(
                "The section shows the geometry of section 2 of the foundation document in one "
                "picture: the target axis lies in the plane of the drawing, at eps* from the "
                "lamella normal. The pattern is what the lamella shows down that axis."
            ),
            section="evidence",
            figures=(figures["fib_section"], figures["fib_saed"]),
        ),
        ResultStage(
            key="uncertainty",
            title="5. Uncertainty budget",
            summary=budget.describe(),
            table=ResultTable(
                columns=(
                    Column("term", "Term"),
                    Column("value_deg", "Standard uncertainty", units="°", numeric=True, digits=3),
                    Column("source", "Source"),
                ),
                rows=budget_rows,
            ),
            section="diagnostics",
        ),
        ResultStage(
            key="options",
            title="6. Alternative azimuths in the same grain",
            summary="Distinct lamella azimuths for the same grain, best first.",
            table=ResultTable(
                columns=(
                    Column("member", "Axis"),
                    Column("eps_deg", "Residual tilt", units="°", numeric=True, digits=2),
                    Column(
                        "theta_sample_deg", f"Azimuth {_THETA_S}", units="°", numeric=True, digits=2
                    ),
                ),
                rows=option_rows,
            ),
            explanation=(
                "If the best azimuth is blocked - a scratch, a neighbouring grain, the edge of the "
                "specimen - these are the next best cuts through the same grain."
            ),
            section="diagnostics",
        ),
        ResultStage(
            key="frames",
            title="7. Frames and calibration",
            summary=best.chamber.calibration_caveat(),
            explanation=" ".join(
                (best.surface.describe(), best.registration.describe(), best.chamber.describe())
            ),
            status="ok" if best.chamber.calibrated else "warning",
            section="method",
        ),
        ResultStage(
            key="work_order",
            title="8. Work order",
            summary="Everything the FIB operator needs, in the order it is needed.",
            table=ResultTable(
                columns=(Column("item", "Item"), Column("value", "Value")),
                rows=work_rows,
            ),
            section="method",
        ),
        ResultStage(
            key="explanation",
            title="9. The whole answer in words",
            summary=report.describe(),
            section="audit",
            status="info",
        ),
    )
    data = {
        "report": report.to_json_dict(),
        "work_order_html": order_html,
        "describe": report.describe(),
        "best_grain_id": best.grain_id,
        "source": source,
    }
    return AppResult(
        title=f"FIB lamella for {report.orbit.family_text} in {best.phase_name}",
        summary=summary,
        highlights=tuple(highlights),
        warnings=tuple(warnings),
        table=ResultTable(
            columns=_CANDIDATE_COLUMNS,
            rows=_candidate_rows(plans),
            caption=(
                "Candidates in order: feasibility class first, then the weighted score. "
                "Azimuths are of the lamella normal."
            ),
        ),
        data=data,
        inputs={key: value for key, value in request.items() if key != "scan_file"},
        notes=tuple(notes),
        citations=(_CITATION_ITA, _CITATION_DE_GRAEF, _CITATION_GIANNUZZI, _CITATION_GUM),
        stages=stages,
    )


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="fib.example.polycrystal_011",
            title="Rank a polycrystal for a <011> lamella",
            panel="fib_lamella",
            summary=(
                "Twelve nickel grains, a <011> target and a +/-30 degree holder: which grain, "
                "where, and at what azimuth."
            ),
            teaches=(
                "Almost every cubic grain has some <011> equivalent near the surface plane, so "
                "the ranking is decided by clearance and spread as much as by eps*. Read the "
                "stereogram: the ringed member sits just inside the outer circle, and that small "
                "gap is exactly the tilt the holder must supply."
            ),
            operation="fib.lamella_plan",
            request={"dataset": "equiaxed_polycrystal", "zone_axis": [0, 1, 1]},
        ),
        ExampleScenario(
            id="fib.example.polycrystal_111",
            title="A harder target: <111>, and the unknown mount",
            panel="fib_lamella",
            summary=(
                "The same map for <111>, which has only four upper-hemisphere members, so "
                "several grains are only probabilistically reachable."
            ),
            teaches=(
                "With fewer equivalents the residual grows and the mounting rotation starts to "
                "matter: compare the polar plot of a guaranteed grain with one whose arcs are "
                "broken - that grain works only if the lamella happens to land the right way "
                "round on the grid."
            ),
            operation="fib.lamella_plan",
            request={"dataset": "equiaxed_polycrystal", "zone_axis": [1, 1, 1]},
        ),
        ExampleScenario(
            id="fib.example.calibrated_52",
            title="One chosen grain on a calibrated 52-degree instrument",
            panel="fib_lamella",
            summary=(
                "Plan grain 0 for <001> on a chamber whose fiducial calibration gave a reversed "
                "sense and a 90 degree offset."
            ),
            teaches=(
                "The sample azimuth is a property of the crystal; the FIB rotation is a property "
                "of the instrument. Here the pattern rotation differs from the sample azimuth by "
                "the calibrated sense and offset, and the UNCALIBRATED warning is gone because "
                "the calibration was declared."
            ),
            operation="fib.lamella_plan",
            request={
                "dataset": "equiaxed_polycrystal",
                "zone_axis": [0, 0, 1],
                "grain_id": 0,
                "column_angle": "52",
                "rotation_sense": "-1",
                "rotation_offset_deg": 90.0,
                "chamber_calibrated": True,
            },
        ),
        ExampleScenario(
            id="fib.example.affine_registration",
            title="Registering the SEM image with control points",
            panel="fib_lamella",
            summary=(
                "The twinned map for <011>, with the SEM image registered to the scan by four "
                "control points that do not fit perfectly."
            ),
            teaches=(
                "The rectangle and the azimuth are now reported in SEM pixels, and the affine "
                "misfit becomes a term of the uncertainty budget: a sloppy registration widens "
                "eps* +/- U and can demote a plan from guaranteed."
            ),
            operation="fib.lamella_plan",
            request={
                "dataset": "sigma3_twin",
                "zone_axis": [0, 1, 1],
                "registration": "affine",
                "control_points": (
                    "0 0 20 10\n60 0 139 25\n0 60 5 130\n60 60 125 146\n30 30 72 79"
                ),
            },
        ),
    )
)
