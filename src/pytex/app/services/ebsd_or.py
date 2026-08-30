"""The orientation relationship of a real microstructure, from measured grains.

What it does
    Takes the orientations of grains that were *measured* — a parent grain and
    a product grain, one pair per line, as they come out of an indexing run —
    and answers the question those numbers were collected to answer: what
    orientation relationship is this material transforming through? The answer
    is a fitted rotation, the name of the catalogued relationship it matches,
    the parallel planes and directions that state it the way a paper states it,
    and a ranked list of the coincident directions the fit admits.

Why it lives in the EBSD workspace
    Its input is a column of Euler angles from a scan, which is the EBSD
    subject. The catalogue side of the same question — "what does Burgers look
    like, and what variants does it give?" — lives in the Variants panel, and
    the two are joined by the handoff: identify the relationship here, then see
    it drawn there.

Why several pairs and not one
    One pair fits any rotation exactly. Its residual is zero by construction,
    so a single pair can never be contradicted by its own data and gives no
    evidence that the relationship is real. Several pairs give a *scatter*,
    which is the only number in the answer that says whether the fit means
    anything — and pairs from different variants are welcome, because the
    symmetry reduction absorbs the parent operator that distinguishes them.

The four angles, which are not interchangeable
    The *scatter* is how far the measured pairs sit from one fitted rotation.
    The *catalogue distance* is how far that fit sits from a named
    relationship, and is what identifies it. The *rationalization cost* is what
    writing the fit in integers costs. The *clause deviation* is how far one
    index pair sits from the exact image. Each is labelled where it appears.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    ExampleScenario,
    IntegerParameter,
    NumberParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import (
    direction_label,
    phase_parameter,
    plane_label,
    relationship_name,
)
from pytex.app.services.crystal import _EULER_CONVENTIONS, _euler_convention
from pytex.app.services.variants import (
    _ANGLE_MEANINGS,
    _CANONICAL_CHILD,
    _CANONICAL_PARENT,
    _axis_in_basis,
)

__all__: tuple[str, ...] = ()

_CITATION_BUNGE = "Bunge, Texture Analysis in Materials Science (1982), chapter 2."
_CITATION_BURGERS = (
    "Burgers, Physica 1 (1934) 561 (the bcc-to-hcp relationship in zirconium)."
)

#: The panel this operation belongs to. Its own sub-tab of the EBSD workspace.
_PANEL = "ebsd_or"

#: The default measurement: a beta grain and three alpha grains descended from
#: it, one per variant, computed exactly through Burgers.
#:
#: Exact rather than noisy, and stated as such in the help, so the panel opens
#: on an answer that can be checked rather than on a plausible-looking number.
#: Three pairs rather than one because one pair fits any rotation exactly: the
#: opening screen has to show a scatter, or the panel teaches that a single
#: measurement is conclusive, which is the error it exists to prevent.
_DEFAULT_PAIRS = """# parent phi1 Phi phi2   |   child phi1 Phi phi2   (degrees, Bunge)
# A beta-zirconium grain at (30, 40, 10) and three alpha grains grown from it
# through Burgers variants 1, 5 and 9. Exact, so the answer is checkable.
30 40 10    167.5709  58.2280   0.9653
30 40 10    338.2303  62.4354  19.6967
30 40 10     91.7913 111.6347 189.7671"""


def _parse_pairs(text: str) -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    """Read the grain-pair table: six numbers a line, parent then child.

    The format is what a user actually has — six columns pasted out of a
    spreadsheet or an indexing export — rather than a form of six boxes, which
    cannot express more than one pair and is the reason the single-pair view
    could not answer this question. Blank lines and ``#`` comments are ignored
    so the default can carry its own header, and commas count as separators
    because a CSV pasted whole should work.

    Errors name the line, because "row 4 has five numbers" is actionable and
    "invalid input" is not.
    """

    pairs: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = []
    for number, raw in enumerate(str(text).splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = [field for field in line.replace(",", " ").replace("\t", " ").split() if field]
        if len(fields) != 6:
            raise InvalidInputError(
                f"Line {number} has {len(fields)} numbers; a grain pair needs six.",
                field="pairs",
                hint=(
                    "One pair per line: the parent's three Euler angles, then the child's "
                    "three, in degrees. Blank lines and lines starting with # are ignored."
                ),
            )
        try:
            values = [float(field) for field in fields]
        except ValueError as error:
            raise InvalidInputError(
                f"Line {number} contains something that is not a number.",
                field="pairs",
                hint="Euler angles are decimal degrees, e.g. `30 40 10 167.57 58.23 0.97`.",
            ) from error
        if not all(math.isfinite(value) for value in values):
            raise InvalidInputError(
                f"Line {number} contains a value that is not finite.", field="pairs"
            )
        pairs.append(
            ((values[0], values[1], values[2]), (values[3], values[4], values[5]))
        )
    if not pairs:
        raise InvalidInputError(
            "No grain pairs were given, so there is nothing to fit.",
            field="pairs",
            hint="Enter at least one line of six Euler angles: parent, then child.",
        )
    return pairs


def _statement_rows(
    statements: Any, *, kind: str, parent_spec: PhaseSpec, child_spec: PhaseSpec
) -> list[dict[str, Any]]:
    """Every candidate parallelism of one kind, relabelled in each phase's notation.

    The core returns three-index labels for both sides. A hexagonal child
    indexed three ways beside a four-index plane in the next column would be
    two notations for one crystal in one table, so both sides are relabelled
    here through the same helpers the rest of the application uses.

    The list is the *ranking*, not the winner: the question "what are the best
    coincident directions" is a search over families, and showing only the
    chosen pair hides whether it won clearly or was picked from a tie.
    """

    label = plane_label if kind == "plane" else direction_label
    rows: list[dict[str, Any]] = []
    for statement in statements:
        parent_indices = tuple(int(value) for value in statement.parent_indices)
        child_indices = tuple(int(value) for value in statement.child_indices)
        rows.append(
            {
                "kind": kind,
                "parent": label(parent_indices, spec=parent_spec),
                "child": label(child_indices, spec=child_spec),
                "parent_indices": list(parent_indices),
                "child_indices": list(child_indices),
                "deviation_deg": float(statement.deviation_deg),
            }
        )
    return rows


def _variant_assignments(
    report: Any,
    *,
    parent_orientations: list[Any],
    child_orientations: list[Any],
) -> list[dict[str, Any]]:
    """Which variant of the fitted relationship each measured child sits on.

    This is the question that follows immediately from a named relationship and
    is not answered by naming it: given the parent, *which* of the twelve
    admissible children is this grain? Answered by predicting every variant's
    child orientation from the measured parent and taking the nearest, with the
    distance reported — a large distance means the pair does not belong to this
    parent at all, and is worth seeing rather than rounding to a variant index.
    """

    variants = report.relationship.generate_variants()
    count = len(variants)
    rows: list[dict[str, Any]] = []
    for index, (parent, child) in enumerate(
        zip(parent_orientations, child_orientations, strict=True), start=1
    ):
        parent_matrix = np.asarray(parent.rotation.as_matrix(), dtype=float)
        best_index = 0
        best_angle = math.inf
        for variant in variants:
            rotation = np.asarray(variant.parent_to_child_rotation.as_matrix(), dtype=float)
            predicted = _orientation_like(child, parent_matrix @ rotation.T)
            angle = float(child.misorientation_to(predicted).angle_deg)
            if angle < best_angle:
                best_angle = angle
                best_index = int(variant.variant_index)
        rows.append(
            {
                "pair": index,
                "variant": best_index,
                "variant_count": count,
                "distance_deg": best_angle,
                "residual_deg": float(report.residuals_deg[index - 1]),
            }
        )
    return rows


def _orientation_like(template: Any, matrix: np.ndarray) -> Any:
    """An orientation with the given matrix, carrying the template's phase and frame."""

    from pytex.core.orientation import Orientation, Rotation

    return Orientation(
        rotation=Rotation.from_matrix(matrix),
        specimen_frame=template.specimen_frame,
        crystal_frame=template.crystal_frame,
        symmetry=template.symmetry,
        phase=template.phase,
    )


_PARAMETERS: tuple[Any, ...] = (
    phase_parameter(
        label="Parent phase",
        help_text="The phase of the grains that transformed — beta zirconium in the default.",
        builtin=_CANONICAL_PARENT,
    ),
    phase_parameter(
        name="child_phase",
        label="Product phase",
        help_text="The phase of the product grains — alpha zirconium in the default.",
        builtin=_CANONICAL_CHILD,
    ),
    TextParameter(
        name="pairs",
        label="Measured grain pairs",
        help_text=(
            "One pair per line: the parent grain's three Euler angles, then the product "
            "grain's three, in degrees. Blank lines and `#` comments are ignored, and commas "
            "count as separators so a pasted CSV works.\n\n"
            "Several pairs are worth far more than one. A single pair fits any rotation "
            "exactly, so its scatter is zero by construction and proves nothing; several "
            "pairs give a scatter that can contradict the fit. Pairs from *different* "
            "variants are welcome — the symmetry reduction absorbs the parent operator that "
            "tells variants apart, so they average rather than fight."
        ),
        multiline=True,
        default=_DEFAULT_PAIRS,
        placeholder="30 40 10   167.5709 58.2280 0.9653",
    ),
    ChoiceParameter(
        name="euler_convention",
        label="Euler convention",
        help_text=(
            "Which axis sequence the angles name. Every grain is read in the same "
            "convention, because they came from one indexing run."
        ),
        options=_EULER_CONVENTIONS,
        default="bunge",
    ),
    NumberParameter(
        name="catalog_tolerance_deg",
        label="Naming tolerance",
        help_text=(
            "How close the fit must sit to a catalogued relationship before it is named. "
            "Three degrees is the working figure: above the orientation noise of a "
            "well-calibrated map, below the separation between the catalogued relationships."
        ),
        units="deg",
        default=3.0,
        minimum=0.1,
        maximum=15.0,
    ),
    IntegerParameter(
        name="max_index",
        label="Largest index in the statement",
        help_text=(
            "Bound on the integers the rational statement may use. Two gives the tidiest "
            "statement and the largest cost; raising it buys a closer one with untidier "
            "indices."
        ),
        default=3,
        minimum=1,
        maximum=6,
        advanced=True,
    ),
    IntegerParameter(
        name="max_statements",
        label="How many parallelisms to rank",
        help_text=(
            "How many candidate plane pairs and direction pairs to report. The extras are "
            "the runners-up: they say whether the chosen pair won clearly or was picked out "
            "of a near-tie, which the winner alone cannot."
        ),
        default=5,
        minimum=1,
        maximum=12,
        advanced=True,
    ),
)


@REGISTRY.operation(
    "ebsd.or_from_grains",
    title="Orientation relationship from measured grains",
    summary="Measured parent and product orientations in; the relationship they show out.",
    help_text=(
        "The everyday question of a partially transformed microstructure, answered end to "
        "end. Index a few parent grains and the product grains inside them, paste the Euler "
        "angles in, and get the rotation fitted to the whole set, its distance from every "
        "catalogued relationship, a conclusive-or-not verdict, the statement in integers, and "
        "the coincident directions the fit admits — ranked, with the runners-up visible.\n\n"
        "**The defaults are an exact Burgers measurement.** A beta-zirconium grain at "
        "(30, 40, 10) and three alpha grains grown from it through variants 1, 5 and 9. They "
        "are computed rather than measured, so the panel opens on an answer that can be "
        "checked: Burgers, at essentially zero, with zero scatter. Replace them with your own "
        "and the scatter becomes the number that matters.\n\n"
        "**Four different angles appear here and they are not interchangeable.** The "
        "*scatter* is how far the measured pairs sit from one fitted rotation — zero for a "
        "single pair, by construction, which is why one pair can never be contradicted by its "
        "own residual. The *catalogue distance* is how far the fit sits from a named "
        "relationship, and it is what identifies it. The *rationalization cost* is what "
        "writing the fit in integers costs. The *clause deviation* is how far one index pair "
        "sits from the exact image. Each is labelled where it appears.\n\n"
        "**Naming is deliberately conservative.** A relationship is called conclusive only "
        "when it fits within the tolerance *and* leads the runner-up by more than the data's "
        "own scatter. Two relationships a degree apart cannot be told apart by data that "
        "scatters by two, and the honest answer is to say so."
    ),
    parameters=_PARAMETERS,
    returns=(
        "The catalogue ranking as the table; the fit, the verdict, the rational statement, "
        "the ranked coincident planes and directions, and one row per measured pair with the "
        "variant it sits on, under `data`."
    ),
    panel=_PANEL,
    citations=(_CITATION_BUNGE, _CITATION_BURGERS),
    tags=(
        "orientation relationship",
        "OR",
        "EBSD",
        "measured",
        "grain",
        "Euler",
        "Burgers",
        "variant",
        "rationalization",
    ),
)
def _or_from_grains(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.frame_catalog import specimen_frame
    from pytex.core.orientation import Orientation, OrientationSet
    from pytex.core.transformation import characterize_orientation_relationship

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    convention = _euler_convention(request["euler_convention"])
    frame = specimen_frame()
    pairs = _parse_pairs(str(request["pairs"]))
    tolerance = float(request["catalog_tolerance_deg"])
    max_index = int(request["max_index"])
    max_statements = int(request["max_statements"])

    def _orientation(angles: tuple[float, float, float], phase: Any) -> Any:
        return Orientation.from_euler(
            angles[0],
            angles[1],
            angles[2],
            specimen_frame=frame,
            symmetry=phase.symmetry,
            phase=phase,
            convention=convention,
            degrees=True,
        )

    parent_orientations = [_orientation(parent, parent_phase) for parent, _ in pairs]
    child_orientations = [_orientation(child, child_phase) for _, child in pairs]

    try:
        report = characterize_orientation_relationship(
            OrientationSet.from_orientations(parent_orientations),
            OrientationSet.from_orientations(child_orientations),
            catalog_tolerance_deg=tolerance,
            max_index=max_index,
            max_statements=max_statements,
        )
    except ValueError as error:
        raise InvalidInputError(
            f"These grains cannot be characterized: {error}",
            field="phase",
            hint=(
                "A relationship is defined between two distinct phases; choose a parent and "
                "a product that are not the same phase, and give at least one pair."
            ),
        ) from error

    order = np.argsort(np.asarray(report.catalog_deviations_deg, dtype=float))
    catalog_rows = [
        {
            "relationship": relationship_name(report.catalog_names[int(index)]),
            "deviation_deg": float(report.catalog_deviations_deg[int(index)]),
            "within_tolerance": (
                "yes" if report.catalog_deviations_deg[int(index)] <= tolerance else "no"
            ),
        }
        for index in order
    ]

    misorientation = report.relationship.misorientation()
    axis = np.asarray(misorientation.rotation.axis, dtype=float)
    fit_euler = report.relationship.parent_to_child_rotation.to_bunge_euler(degrees=True)

    statement: dict[str, Any] | None = None
    statement_note: str | None = None
    try:
        rationalized = report.as_rational_relationship(
            max_index=max_index, tolerance_deg=tolerance
        )
    except ValueError as error:
        # Not an input error: a rotation with no low-index statement is a real
        # answer, and the panel says so rather than failing the whole run.
        statement_note = str(error)
    else:
        statement = {
            "text": rationalized.statement,
            "plane": {
                "parent": plane_label(
                    tuple(int(v) for v in rationalized.plane_statement.parent_indices),
                    spec=parent_spec,
                ),
                "child": plane_label(
                    tuple(int(v) for v in rationalized.plane_statement.child_indices),
                    spec=child_spec,
                ),
                "deviation_deg": float(rationalized.plane_statement.deviation_deg),
            },
            "direction": {
                "parent": direction_label(
                    tuple(int(v) for v in rationalized.direction_statement.parent_indices),
                    spec=parent_spec,
                ),
                "child": direction_label(
                    tuple(int(v) for v in rationalized.direction_statement.child_indices),
                    spec=child_spec,
                ),
                "deviation_deg": float(rationalized.direction_statement.deviation_deg),
            },
            "cost_deg": float(rationalized.residual_rotation_deg),
            "zone_law_deviation_deg": float(rationalized.zone_law_deviation_deg),
            "max_index": int(rationalized.max_index),
        }

    planes = _statement_rows(
        report.plane_statements, kind="plane", parent_spec=parent_spec, child_spec=child_spec
    )
    directions = _statement_rows(
        report.direction_statements,
        kind="direction",
        parent_spec=parent_spec,
        child_spec=child_spec,
    )
    assignments = _variant_assignments(
        report,
        parent_orientations=parent_orientations,
        child_orientations=child_orientations,
    )

    best_label = (
        None if report.best_catalog_name is None else relationship_name(report.best_catalog_name)
    )
    verdict = (
        f"{best_label} at {report.best_catalog_deviation_deg:.2f} degrees"
        if report.is_conclusive
        else "no relationship can be named conclusively"
    )
    result = AppResult(
        title=(
            f"{parent_spec.name} to {child_spec.name}: "
            f"{report.pair_count} measured pair(s)"
        ),
        summary=(
            f"Fitted to {report.pair_count} pair(s) with a scatter of "
            f"{report.mean_residual_deg:.3f} degrees; the fit is "
            f"{report.best_catalog_deviation_deg:.3f} degrees from "
            f"{best_label or 'no catalogued relationship'} and leads the runner-up by "
            f"{report.margin_deg:.3f} degrees, so the verdict is {verdict}. The scatter, the "
            "catalogue distance, the rationalization cost and the clause deviations are four "
            "different angles and are labelled separately wherever they appear."
        ),
        table=ResultTable(
            columns=(
                Column("relationship", "Relationship"),
                Column(
                    "deviation_deg",
                    "Catalogue distance",
                    units="°",
                    numeric=True,
                    digits=3,
                    help_text=_ANGLE_MEANINGS["catalog"],
                ),
                Column(
                    "within_tolerance",
                    "Within tolerance",
                    help_text="Whether this relationship is close enough to be named.",
                ),
            ),
            rows=tuple(catalog_rows),
            caption=(
                "Every catalogued relationship for these two crystal systems, ordered by how "
                "far the fitted rotation sits from it. A conclusive naming needs the winner to "
                "lead the runner-up by more than the scatter and its own misfit."
            ),
        ),
        data={
            "fit": {
                "angle_deg": float(misorientation.angle_deg),
                "axis": [float(value) for value in axis],
                "axis_parent": _axis_in_basis(axis, parent_phase, parent_spec),
                "axis_child": _axis_in_basis(axis, child_phase, child_spec),
                "euler_deg": [float(value) for value in fit_euler],
                "matrix": [
                    [float(value) for value in row]
                    for row in report.relationship.parent_to_child_rotation.as_matrix()
                ],
                "mean_residual_deg": float(report.mean_residual_deg),
                "max_residual_deg": float(report.max_residual_deg),
                "pair_count": int(report.pair_count),
                "converged": bool(report.converged),
            },
            "naming": {
                "best": report.best_catalog_name,
                "best_label": best_label,
                "best_deviation_deg": float(report.best_catalog_deviation_deg),
                "margin_deg": float(report.margin_deg),
                "is_conclusive": bool(report.is_conclusive),
                "tolerance_deg": tolerance,
            },
            "catalog": catalog_rows,
            "statement": statement,
            "statement_note": statement_note,
            "coincidences": {"planes": planes, "directions": directions},
            "pairs": [
                {
                    "pair": row["pair"],
                    "variant": row["variant"],
                    "variant_count": row["variant_count"],
                    "distance_deg": row["distance_deg"],
                    "residual_deg": row["residual_deg"],
                    "parent_euler": list(pairs[row["pair"] - 1][0]),
                    "child_euler": list(pairs[row["pair"] - 1][1]),
                }
                for row in assignments
            ],
            "phases": {"parent": parent_spec.name, "child": child_spec.name},
            "variant_count": (
                assignments[0]["variant_count"] if assignments else 0
            ),
            "angle_meanings": dict(_ANGLE_MEANINGS),
            "euler_convention": convention,
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "pairs": str(request["pairs"]),
            "euler_convention": convention,
            "catalog_tolerance_deg": tolerance,
            "max_index": max_index,
            "max_statements": max_statements,
        },
        citations=(_CITATION_BUNGE, _CITATION_BURGERS),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="ebsd_or.example.burgers_three_variants",
            title="Three alpha grains from one beta grain",
            panel=_PANEL,
            summary="The canonical case: exact Burgers pairs, named at zero with zero scatter.",
            teaches=(
                "What a clean answer looks like, so a real one can be judged against it. "
                "Three product grains from one parent, on three different variants, all "
                "reduce to the same rotation: the scatter is zero, the catalogue distance is "
                "zero, and Burgers leads the runner-up by tens of degrees.\n\n"
                "Note that the three pairs sit on different variants and still average. That "
                "is the symmetry reduction doing its work — the parent symmetry operator is "
                "exactly what distinguishes one variant from another, and it is absorbed "
                "before the pairs are compared."
            ),
            operation="ebsd.or_from_grains",
            request={
                "phase": {"builtin": _CANONICAL_PARENT},
                "child_phase": {"builtin": _CANONICAL_CHILD},
                "pairs": _DEFAULT_PAIRS,
                "euler_convention": "bunge",
                "catalog_tolerance_deg": 3.0,
                "max_index": 3,
                "max_statements": 5,
            },
        ),
        ExampleScenario(
            id="ebsd_or.example.wrong_parent",
            title="One grain that does not belong to this parent",
            panel=_PANEL,
            summary="Two good pairs and one misattributed grain: the panel refuses to name it.",
            teaches=(
                "The third product grain was grown from a *different* beta grain and reported "
                "here against this one, which is exactly what a mis-segmented map delivers. "
                "The fit is dragged 7.7 degrees off Burgers and the scatter rises to about 10 "
                "degrees, so the verdict becomes **inconclusive** — not because the "
                "calculation failed, but because a set of pairs that scatters by ten degrees "
                "cannot distinguish relationships that sit closer together than that.\n\n"
                "Read the pair table to find the culprit: two rows sit within a hundredth of "
                "a degree of a variant and the third is tens of degrees from every one of "
                "them. That column exists for this case. A panel that answered 'Burgers' here "
                "would be reporting the two grains that agreed and hiding the one that did "
                "not."
            ),
            operation="ebsd.or_from_grains",
            request={
                "phase": {"builtin": _CANONICAL_PARENT},
                "child_phase": {"builtin": _CANONICAL_CHILD},
                "pairs": (
                    "30 40 10   167.5709  58.2280   0.9653\n"
                    "30 40 10   338.2303  62.4354  19.6967\n"
                    "30 40 10    45.0789  65.6736  47.9195"
                ),
                "euler_convention": "bunge",
                "catalog_tolerance_deg": 3.0,
                "max_index": 3,
                "max_statements": 5,
            },
        ),
        ExampleScenario(
            id="ebsd_or.example.noisy_pairs",
            title="The same grains, measured with noise",
            panel=_PANEL,
            summary="Half a degree of orientation noise: still Burgers, and now with a scatter.",
            teaches=(
                "The same three pairs with the product angles perturbed by a few tenths of a "
                "degree, which is what a well-calibrated map delivers. The relationship is "
                "still named, and the number that changed is the scatter — which is the point "
                "of collecting more than one pair.\n\n"
                "Compare the scatter with the gap between the first and second rows of the "
                "catalogue table. While the gap is much larger than the scatter, the naming "
                "is safe; when they become comparable, the panel stops calling it conclusive, "
                "and that is the honest answer rather than a failure."
            ),
            operation="ebsd.or_from_grains",
            request={
                "phase": {"builtin": _CANONICAL_PARENT},
                "child_phase": {"builtin": _CANONICAL_CHILD},
                "pairs": (
                    "30 40 10   167.94  58.51   0.62\n"
                    "30 40 10   338.02  62.13  20.05\n"
                    "30 40 10    91.51 111.94 189.41"
                ),
                "euler_convention": "bunge",
                "catalog_tolerance_deg": 3.0,
                "max_index": 3,
                "max_statements": 5,
            },
        ),
    )
)
