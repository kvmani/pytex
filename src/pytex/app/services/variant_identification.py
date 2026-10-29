# ruff: noqa: RUF001
"""Which variant is this child grain? The calculation, shown step by step.

What it does
    Takes the measured Euler angles of a parent grain and of a neighbouring
    child grain, and an orientation relationship, and identifies which variant
    of that relationship the child is. The answer comes with its derivation:
    every equation used, the measured values substituted into it, and the
    number each step produced, typeset as a paper would print them. The
    symmetry operators that turned out to match -- the parent operator that
    generated the variant and the child operator that reconciles the measured
    child's crystal setting with the prediction -- are written out as matrices.

When to use it
    For a few hand-picked grain pairs, when the reasoning has to be followed or
    documented: teaching, checking a reconstruction, writing a methods section.
    For a whole map use the EBSD OR panel, which assigns the same variants
    without the derivation.

How the steps are ordered
    The answer first (variant, misfit, margin, verdict), then the derivation in
    the order a reader verifies it: the two orientation matrices, the
    relationship's rotation, the measured rotation, the variants, the misfit of
    each, the winner substituted in full, the physical parallelism check, and
    the decision. Equations are numbered once across the whole derivation so
    the prose can cite them.

Every number shown is read from
:class:`pytex.core.variant_identification.VariantIdentificationReport` or is an
algebraic rearrangement of its fields (a product of two of its matrices, a
trace); the tests hold each displayed quantity to the report.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pytex.app import mathml as mm
from pytex.app.errors import InvalidInputError
from pytex.app.figures import COLORS, render_figure
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    ExampleScenario,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import (
    AppResult,
    Column,
    ResultEquation,
    ResultFigure,
    ResultMetric,
    ResultStage,
    ResultTable,
)
from pytex.app.services.calculator import (
    _RELATIONSHIPS,
    _four_index,
    custom_relationship_parameters,
    direction_label,
    plane_label,
    relationship_name,
    resolve_relationship,
)
from pytex.app.services.crystal import _euler_convention
from pytex.app.services.variants import (
    _CANONICAL_RELATIONSHIP,
    _CITATION_BUNGE,
    _CITATION_MORITO,
    _MEASURED_PAIR_PARAMETERS,
    _euler_angles,
    _nearest_axis,
    _orientation_from_euler,
)
from pytex.core._angles import rotation_angle_from_matrix_rad
from pytex.core.miller import canonicalize_sign
from pytex.core.symbols import symbol_text
from pytex.core.variant_identification import (
    DEFAULT_ASSIGNMENT_TOLERANCE_DEG,
    VariantIdentificationReport,
    identify_variant,
    operator_axis_angle,
)

__all__: tuple[str, ...] = ()

_CITATION_HUMBERT = (
    "Humbert, Gey, Muller & Esling, J. Appl. Cryst. 29 (1996) 662 (variant identification "
    "from measured orientations)."
)
_CITATION_GEY = (
    "Gey & Humbert, Acta Mater. 50 (2002) 277 (variant identification in Ti and Zr, Burgers)."
)

#: Euler conventions: the axis of each factor and the registered angle symbols.
_EULER_FORMS: dict[str, tuple[tuple[str, str, str], tuple[str, str, str]]] = {
    "bunge": (("z", "x", "z"), ("phi_1", "Phi", "phi_2")),
    "matthies": (("z", "y", "z"), ("euler_alpha", "euler_beta", "euler_gamma")),
}

_VERDICT_TEXT = {
    "assigned": "assigned (unambiguous)",
    "ambiguous": "ambiguous",
    "no_fit": "no variant fits",
}


# -- small typesetting helpers ------------------------------------------------------


def _s(name: str) -> mm.Node:
    """A registered symbol."""

    return mm.symbol(name)


def _sub(name: str, script: str | mm.Node) -> mm.Node:
    """A registered symbol with an extra subscript (``P`` → ``P`` sub ``p``)."""

    script_node = mm.mi(script) if isinstance(script, str) else script
    return mm.sub(_s(name), script_node)


def _transpose(node: mm.Node) -> mm.Node:
    return mm.sup(node, mm.mi("T", font="sans"))


def _degrees(value: float, digits: int) -> mm.Node:
    return mm.row(mm.number(value, digits), mm.mo("°"))


def _index_node(k: int | str) -> mm.Node:
    return mm.mn(str(k)) if isinstance(k, int) else mm.mi(k)


def _base(name: str) -> mm.Node:
    """The base glyph of a registered subscripted symbol (``V`` of ``V_k``)."""

    node = _s(name)
    return node["children"][0] if node["tag"] in {"msub", "msup", "msubsup"} else node


def _child_operator_star() -> mm.Node:
    """``S_c*``: the registered child operator, starred as the matching one."""

    return mm.subsup(_base("child_operator"), mm.mi("c"), mm.mo("*"))


def _parent_operator_k(k: int | str) -> mm.Node:
    """``S_p,k``: the registered parent operator, indexed by variant."""

    return mm.sub(_base("parent_operator"), mm.row(mm.mi("p"), mm.mo(","), _index_node(k)))


def _variant_k(k: int | str) -> mm.Node:
    return mm.sub(_base("variant_rotation"), _index_node(k))


def _predicted_child(k: int | str) -> mm.Node:
    return mm.sub(_s("orientation_child"), _index_node(k))


def _misfit_k(k: int | str) -> mm.Node:
    return mm.sub(_base("variant_misfit"), _index_node(k))


def _angle_of(argument: mm.Node) -> mm.Node:
    """``ω(X)``: the rotation angle of a rotation matrix."""

    return mm.row(_s("omega"), mm.fenced(argument))


def _miller_node(indices: Sequence[int], *, spec: PhaseSpec, family: str = "plane") -> mm.Node:
    """Indices in the notation natural to the phase, four-index for hexagonal."""

    from pytex.core.miller import direction_uvw_to_uvtw_array, plane_hkl_to_hkil_array

    values = tuple(int(value) for value in indices)
    if _four_index(spec):
        converter = plane_hkl_to_hkil_array if family == "plane" else direction_uvw_to_uvtw_array
        values = tuple(int(value) for value in converter(np.asarray(values))[0])
    return mm.miller(values, family=family)


def _direction_indices(member: Any) -> tuple[int, int, int]:
    coordinates = np.asarray(member.coordinates, dtype=float)
    return (
        round(float(coordinates[0])),
        round(float(coordinates[1])),
        round(float(coordinates[2])),
    )


def _plane_indices(member: Any) -> tuple[int, int, int]:
    values = np.asarray(member.miller.indices, dtype=int)
    return (int(values[0]), int(values[1]), int(values[2]))


@dataclass
class _Equations:
    """Numbers the derivation's equations continuously across stages."""

    count: int = 0
    collected: list[ResultEquation] = field(default_factory=list)

    def add(self, *items: mm.Node | str, lead: str | None = None, note: str | None = None) -> str:
        self.count += 1
        tag = str(self.count)
        self.collected.append(
            ResultEquation(math=mm.equation(*items), lead=lead, tag=tag, note=note)
        )
        return tag

    def take(self) -> tuple[ResultEquation, ...]:
        taken = tuple(self.collected)
        self.collected = []
        return taken


# -- the operation ------------------------------------------------------------------

_PARAMETERS = (
    *_MEASURED_PAIR_PARAMETERS[:2],
    ChoiceParameter(
        name="relationship",
        label="Orientation relationship",
        help_text=(
            "The relationship the two grains are tested against; its variants are the "
            "candidates. **Custom** takes the parallelisms from the boxes below, for a "
            "relationship from a paper or one fitted to your own data."
        ),
        options=_RELATIONSHIPS,
        default=_CANONICAL_RELATIONSHIP,
    ),
    *_MEASURED_PAIR_PARAMETERS[2:9],
    NumberParameter(
        name="tolerance_deg",
        label="Assignment tolerance",
        help_text=(
            "The misfit inside which a variant is accepted as the child's. Five degrees sits "
            "above EBSD noise and the natural scatter of a transformation about its nominal "
            "relationship, and below the smallest angle between two variants (10.53 degrees "
            "for Kurdjumov-Sachs), so two variants cannot both be inside it."
        ),
        units="deg",
        default=DEFAULT_ASSIGNMENT_TOLERANCE_DEG,
        minimum=0.1,
        maximum=30.0,
        symbol="assignment_tolerance",
    ),
    IntegerParameter(
        name="digits",
        label="Decimals in the derivation",
        help_text=(
            "How many decimals the matrices and angles of the derivation are printed with. "
            "Every computation is carried at full precision whatever is chosen; this changes "
            "only what is printed."
        ),
        default=4,
        minimum=2,
        maximum=8,
        advanced=True,
    ),
    *custom_relationship_parameters(),
)


@REGISTRY.operation(
    "variants.identify_variant",
    title="Identify the variant of a measured child grain",
    summary="Two measured orientations in; which variant the child is, derived step by step.",
    help_text=(
        "Enter the Euler angles of a parent grain and of a child grain next to it, as your EBSD "
        "software exported them, and choose the orientation relationship the transformation "
        "follows. The result says which variant of that relationship the child is — and shows "
        "how: every equation, the measured values substituted into it, and what each step "
        "produced.\n\n"
        "**The derivation, in nine steps.** The Euler angles become orientation matrices "
        "**P** and **C**. The relationship becomes its rotation **R**, built from the plane and "
        "direction it holds parallel. The measured rotation is **V**ₘ = **C**ᵀ**P**. The "
        "variants are **V**ₖ = **R** **S**ₚ,ₖᵀ for the parent symmetry operators **S**ₚ,ₖ. "
        "Each variant predicts the child **C**ₖ = **P** **V**ₖᵀ, and its misfit θₖ is the "
        "smallest rotation from that prediction to the measured child, over the child's own "
        "symmetry operators **S**c — because the measured child's axes could have been "
        "labelled in any of its equivalent settings. The smallest misfit identifies the "
        "variant; the winning operators are written out in full; the variant's parallel plane "
        "and direction are checked on the two measured grains; and the margin to the "
        "runner-up decides whether the answer is unambiguous.\n\n"
        "**Variant numbers** are PyTex's generation order, not a published table. The "
        "parallelism each variant realizes is stated beside its number, and that statement "
        "is the convention-free identity of the variant."
    ),
    parameters=_PARAMETERS,
    returns=(
        "The misfit of every variant as the table; the derivation as numbered stages with "
        "typeset equations; the identification under `data`."
    ),
    panel="variants",
    citations=(_CITATION_BUNGE, _CITATION_HUMBERT, _CITATION_GEY, _CITATION_MORITO),
    tags=(
        "variant",
        "variant identification",
        "orientation relationship",
        "EBSD",
        "measured",
        "Euler",
        "derivation",
        "step by step",
    ),
)
def _identify_variant(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.frame_catalog import specimen_frame

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    relationship_key = str(request["relationship"])
    relationship = resolve_relationship(request, parent_phase, child_phase)
    label = relationship_name(relationship_key)
    convention = _euler_convention(request["euler_convention"])
    if convention not in _EULER_FORMS:  # pragma: no cover - the choice list is the guard
        raise InvalidInputError(
            f"Unknown Euler convention {convention!r}.", field="euler_convention"
        )
    frame = specimen_frame()
    parent_angles = _euler_angles(request, "parent")
    child_angles = _euler_angles(request, "child")
    parent = _orientation_from_euler(
        parent_angles, phase=parent_phase, convention=convention, frame=frame
    )
    child = _orientation_from_euler(
        child_angles, phase=child_phase, convention=convention, frame=frame
    )
    tolerance = float(request["tolerance_deg"])
    digits = int(request["digits"])
    report = identify_variant(parent, child, relationship, tolerance_deg=tolerance)

    context = _Context(
        report=report,
        parent_spec=parent_spec,
        child_spec=child_spec,
        relationship_label=label,
        convention=convention,
        parent_angles=parent_angles,
        child_angles=child_angles,
        digits=digits,
    )
    stages = context.stages()
    misfit_rows = context.misfit_rows()
    figure = _misfit_figure(context)

    winner = report.variant
    statement = context.variant_statement(winner)
    runner_up = report.runner_up
    margin = report.margin_deg
    verdict = report.verdict
    highlights = [
        ResultMetric(
            "Identified variant",
            f"V{report.variant_index}",
            help_text="Index in PyTex generation order; the parallelism below is its identity.",
        ),
        ResultMetric("Parallelism it realizes", statement),
        ResultMetric(
            f"Misfit {symbol_text('identified_misfit')}",
            report.misfit_deg,
            units="°",
            help_text="Smallest rotation from the predicted to the measured child orientation.",
        ),
    ]
    if runner_up is not None and margin is not None:
        highlights += [
            ResultMetric(
                "Runner-up",
                f"V{runner_up.variant_index} at {report.runner_up_misfit_deg:.2f}°",
            ),
            ResultMetric(
                "Margin to the runner-up",
                margin,
                units="°",
                help_text="How decisively the measurement picks one variant over the rest.",
            ),
        ]
    highlights.append(ResultMetric("Verdict", _VERDICT_TEXT[verdict]))

    warnings: list[str] = []
    if verdict == "no_fit":
        warnings.append(
            f"No variant of {label} lies within the {tolerance:g}° tolerance; the best misfit is "
            f"{report.misfit_deg:.2f}°. Either the transformation does not follow this "
            "relationship, the two grains are not parent and child, or an Euler triple was "
            "entered in a different convention or crystal setting from the one assumed."
        )
    elif verdict == "ambiguous":
        warnings.append(
            f"Variants V{report.variant_index} and "
            f"V{runner_up.variant_index if runner_up else '?'} both lie within the "
            f"{tolerance:g}° tolerance, so this pair cannot tell them apart. Tighten the "
            "tolerance only if the measurement justifies it."
        )

    summary = (
        f"Under {label}, the child grain is variant V{report.variant_index} of "
        f"{len(report.variants)}: it realizes {statement}, and the orientation that variant "
        f"predicts from the measured parent lies {report.misfit_deg:.3f}° from the measured "
        "child."
        + (
            f" The runner-up, V{runner_up.variant_index}, lies "
            f"{report.runner_up_misfit_deg:.2f}° away, {margin:.2f}° behind."
            if runner_up is not None and margin is not None
            else ""
        )
        + f" Verdict at a {tolerance:g}° tolerance: {_VERDICT_TEXT[verdict]}. Orientations are "
        f"{_convention_label(convention)} Euler angles read as crystal-to-specimen matrices; "
        "the derivation below shows every equation with its values."
    )

    result = AppResult(
        title=(f"Variant identification: {parent_spec.name} to {child_spec.name} under {label}"),
        summary=summary,
        table=ResultTable(
            columns=_MISFIT_COLUMNS,
            rows=tuple(misfit_rows),
            caption=(
                "Every variant of the relationship, with its misfit to the measured child and "
                "the child symmetry operator that achieves it."
            ),
        ),
        highlights=tuple(highlights),
        warnings=tuple(warnings),
        stages=stages,
        figures=(figure,),
        data={
            "identification": report.to_json_dict(),
            "variant_index": report.variant_index,
            "misfit_deg": report.misfit_deg,
            "verdict": verdict,
            "statement": statement,
            "euler_convention": convention,
            "relationship": relationship_key,
            "parent_operator_index": int(winner.parent_operator_index),
            "child_operator_index": report.child_operator_index,
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": relationship_key,
            "euler_convention": convention,
            "parent_angle1": parent_angles[0],
            "parent_angle2": parent_angles[1],
            "parent_angle3": parent_angles[2],
            "child_angle1": child_angles[0],
            "child_angle2": child_angles[1],
            "child_angle3": child_angles[2],
            "tolerance_deg": tolerance,
            "digits": digits,
        },
        citations=(_CITATION_BUNGE, _CITATION_HUMBERT, _CITATION_GEY, _CITATION_MORITO),
    )
    return result.to_json()


def _angle_text(value: float) -> str:
    """An entered angle as typed: up to four decimals, no trailing zeros."""

    return f"{value:.4f}".rstrip("0").rstrip(".")


def _convention_label(convention: str) -> str:
    return {"bunge": "Bunge (ZXZ)", "matthies": "Matthies (ZYZ)"}[convention]


_MISFIT_COLUMNS: tuple[Column, ...] = (
    Column("variant", "Variant", numeric=True),
    Column(
        "packet",
        "Packet",
        numeric=True,
        help_text="Variants sharing the parent plane they hold parallel.",
    ),
    Column("parent_plane", "Parent plane held parallel"),
    Column(
        "parent_operator",
        "Generating operator",
        help_text="The parent point-group rotation S_p,k with V_k = R S_p,k^T.",
    ),
    Column(
        "misfit_deg",
        "Misfit",
        units="°",
        numeric=True,
        digits=3,
        help_text="Smallest rotation from this variant's predicted child to the measured child.",
    ),
    Column(
        "child_operator",
        "Matching child operator",
        help_text="The child point-group rotation S_c achieving the misfit.",
    ),
    Column("rank", "Rank", numeric=True),
)


@dataclass
class _Context:
    """Everything the stages need, and the one numbering of their equations."""

    report: VariantIdentificationReport
    parent_spec: PhaseSpec
    child_spec: PhaseSpec
    relationship_label: str
    convention: str
    parent_angles: tuple[float, float, float]
    child_angles: tuple[float, float, float]
    digits: int
    equations: _Equations = field(default_factory=_Equations)

    # -- descriptions ------------------------------------------------------------

    def variant_statement(self, variant: Any) -> str:
        """The variant's own parallelism in plain text, in the measured settings."""

        parts: list[str] = []
        planes = variant.parallel_planes
        directions = variant.parallel_directions
        if planes:
            parent_plane, child_plane = planes[0]
            parts.append(
                f"{plane_label(_plane_indices(parent_plane), spec=self.parent_spec)} ∥ "
                f"{plane_label(_plane_indices(child_plane), spec=self.child_spec)}"
            )
        if directions:
            parent_direction, child_direction = directions[0]
            parts.append(
                f"{direction_label(_direction_indices(parent_direction), spec=self.parent_spec)}"
                " ∥ "
                f"{direction_label(_direction_indices(child_direction), spec=self.child_spec)}"
            )
        return ", ".join(parts) if parts else "no stated parallelism"

    def operator_text(self, operator: np.ndarray, *, parent: bool) -> str:
        """``120° about [111]``: a point-group rotation as a crystallographer says it."""

        axis, angle = operator_axis_angle(operator)
        if angle == 0.0:
            return "identity"
        spec = self.parent_spec if parent else self.child_spec
        phase = (
            self.report.relationship.parent_phase
            if parent
            else self.report.relationship.child_phase
        )
        axis_label, _residual = _nearest_axis(axis, phase, spec)
        return f"{angle:.0f}° about {axis_label}"

    def packets(self) -> list[int]:
        """Packet number of each variant: one per distinct parent plane held parallel."""

        keys: list[tuple[int, ...]] = []
        numbers: list[int] = []
        for variant in self.report.variants:
            planes = variant.parallel_planes
            if not planes:
                numbers.append(1)
                continue
            key = tuple(int(value) for value in canonicalize_sign(_plane_indices(planes[0][0]))[0])
            if key not in keys:
                keys.append(key)
            numbers.append(keys.index(key) + 1)
        return numbers

    def misfit_rows(self) -> list[dict[str, Any]]:
        report = self.report
        order = np.argsort(report.misfits_deg, kind="stable")
        ranks = np.empty(len(order), dtype=int)
        ranks[order] = np.arange(1, len(order) + 1)
        packets = self.packets()
        rows = []
        for position, variant in enumerate(report.variants):
            planes = variant.parallel_planes
            rows.append(
                {
                    "variant": int(variant.variant_index),
                    "packet": packets[position],
                    "parent_plane": (
                        plane_label(_plane_indices(planes[0][0]), spec=self.parent_spec)
                        if planes
                        else "—"
                    ),
                    "parent_operator": self.operator_text(
                        variant.parent_symmetry_operator, parent=True
                    ),
                    "misfit_deg": float(report.misfits_deg[position]),
                    "child_operator": self.operator_text(
                        report.child_operators[int(report.child_operator_indices[position])],
                        parent=False,
                    ),
                    "rank": int(ranks[position]),
                }
            )
        return rows

    # -- the stages --------------------------------------------------------------

    def stages(self) -> tuple[ResultStage, ...]:
        return (
            self._orientation_stage("parent"),
            self._orientation_stage("child"),
            self._relationship_stage(),
            self._measured_stage(),
            self._variants_stage(),
            self._misfit_stage(),
            self._winner_stage(),
            self._parallelism_stage(),
            self._decision_stage(),
        )

    def _orientation_stage(self, role: str) -> ResultStage:
        digits = self.digits
        axes, names = _EULER_FORMS[self.convention]
        angles = self.parent_angles if role == "parent" else self.child_angles
        matrix = self.report.parent_matrix if role == "parent" else self.report.child_matrix
        symbol = "orientation_parent" if role == "parent" else "orientation_child"
        spec = self.parent_spec if role == "parent" else self.child_spec
        step = 1 if role == "parent" else 2

        def elementary(axis: str, angle: mm.Node) -> mm.Node:
            return mm.row(_sub("elementary_rotation", axis), mm.fenced(angle))

        lead = (
            f"The {_convention_label(self.convention)} angles of the {role} grain compose "
            "three elementary rotations into the orientation matrix, which takes "
            f"{role}-crystal components to specimen components:"
        )
        if role == "parent":
            definition = self.equations.add(
                _s(symbol),
                "=",
                *(elementary(axis, _s(name)) for axis, name in zip(axes, names, strict=True)),
                lead=lead,
            )
            omega = _s("omega")
            cos = mm.function("cos", omega)
            sin = mm.function("sin", omega)
            minus_sin = mm.row(mm.mo("−"), sin)
            zero, one = mm.mn("0"), mm.mn("1")
            rz = mm.matrix_of([[cos, minus_sin, zero], [sin, cos, zero], [zero, zero, one]])
            second_axis = axes[1]
            if second_axis == "x":
                r2 = mm.matrix_of([[one, zero, zero], [zero, cos, minus_sin], [zero, sin, cos]])
            else:
                r2 = mm.matrix_of([[cos, zero, sin], [zero, one, zero], [minus_sin, zero, cos]])
            self.equations.add(
                elementary("z", omega),
                "=",
                rz,
                mm.mo(","),
                mm.space(1.5),
                elementary(second_axis, omega),
                "=",
                r2,
                lead=(
                    "with the elementary rotations, active and right-handed, by an angle ω "
                    "about each axis:"
                ),
            )
            cite = f"Eq. ({definition})"
        else:
            cite = "Eq. (1)"
        substituted = self.equations.add(
            _s(symbol),
            "=",
            *(
                elementary(axis, _degrees(value, digits))
                for axis, value in zip(axes, angles, strict=True)
            ),
            "=",
            mm.matrix(matrix, digits=digits),
            lead=(f"Substituting the measured angles of the {role} grain into {cite}:"),
        )
        orthogonality = float(np.max(np.abs(matrix.T @ matrix - np.eye(3))))
        return ResultStage(
            key=f"{role}_orientation",
            title=f"Step {step}. The {role} orientation matrix",
            summary=(
                f"The {role} grain ({spec.name}) was measured at "
                + ", ".join(
                    f"{symbol_text(name)} = {_angle_text(value)}°"
                    for name, value in zip(names, angles, strict=True)
                )
                + f". Its orientation matrix {symbol_text(symbol)} is Eq. ({substituted})."
            ),
            equations=self.equations.take(),
            metrics=(
                ResultMetric(f"det {symbol_text(symbol)}", float(np.linalg.det(matrix))),
                ResultMetric(
                    "Largest entry of PᵀP − I" if role == "parent" else "Largest entry of CᵀC − I",
                    orthogonality,
                    help_text="Zero to rounding for a proper rotation.",
                ),
            ),
            explanation=(
                "Each column of the matrix is one crystal axis written in specimen "
                "coordinates (RD, TD, ND). A determinant of one and an orthogonality error at "
                "rounding level confirm the angles describe a proper rotation."
            ),
            section="derivation",
        )

    def _relationship_stage(self) -> ResultStage:
        digits = self.digits
        relationship = self.report.relationship
        rotation = np.asarray(relationship.parent_to_child_rotation.as_matrix(), dtype=float)
        statement = self.variant_statement_nominal()
        triads = self._triads()
        if triads is not None:
            (n_p, d_p, f_p), (n_c, d_c, f_c) = triads
            (parent_plane, child_plane) = relationship.parallel_planes[0]
            (parent_direction, child_direction) = relationship.parallel_directions[0]
            x_hat, y_hat, n_hat = _s("triad_x"), _s("triad_y"), _s("unit_normal")
            d_hat = _s("unit_direction")
            dot = mm.row(d_hat, mm.mo("⋅"), n_hat)
            numerator = mm.row(d_hat, mm.mo("−"), mm.fenced(dot), n_hat)
            self.equations.add(
                x_hat,
                "=",
                mm.frac(
                    numerator,
                    mm.row(mm.mo("‖", stretchy="false"), numerator, mm.mo("‖", stretchy="false")),
                ),
                mm.mo(","),
                mm.space(1),
                y_hat,
                "=",
                n_hat,
                mm.mo("×"),
                x_hat,
                mm.mo(","),
                mm.space(1),
                _s("statement_triad"),
                "=",
                mm.fenced(mm.row(x_hat, mm.space(0.4), y_hat, mm.space(0.4), n_hat), "[", "]"),
                lead=(
                    "A plane-and-direction statement fixes a rotation through one orthonormal "
                    "triad per phase: the plane normal, the stated direction made exactly normal "
                    "to it, and their cross product."
                ),
            )
            for role, n_vec, d_vec, f_mat, plane, direction, spec in (
                ("p", n_p, d_p, f_p, parent_plane, parent_direction, self.parent_spec),
                ("c", n_c, d_c, f_c, child_plane, child_direction, self.child_spec),
            ):
                self.equations.add(
                    mm.sub(n_hat, mm.mi(role)),
                    "∥",
                    _miller_node(_plane_indices(plane), spec=spec),
                    mm.mo(":"),
                    mm.space(0.4),
                    mm.sub(n_hat, mm.mi(role)),
                    "=",
                    mm.vector(n_vec, digits=digits),
                    mm.mo(","),
                    mm.space(1),
                    mm.sub(d_hat, mm.mi(role)),
                    "∥",
                    _miller_node(_direction_indices(direction), spec=spec, family="direction"),
                    mm.mo(":"),
                    mm.space(0.4),
                    mm.sub(d_hat, mm.mi(role)),
                    "=",
                    mm.vector(d_vec, digits=digits),
                    lead=(
                        f"In the {'parent' if role == 'p' else 'child'} crystal's Cartesian "
                        "frame the stated plane normal and direction are"
                    ),
                )
                self.equations.add(
                    _sub("statement_triad", role),
                    "=",
                    mm.matrix(f_mat, digits=digits),
                )
            tag = self.equations.add(
                _s("or_rotation"),
                "=",
                _sub("statement_triad", "c"),
                _transpose(_sub("statement_triad", "p")),
                "=",
                mm.matrix(rotation, digits=digits),
                lead=(
                    "The rotation that carries the parent triad onto the child triad is the "
                    "relationship itself:"
                ),
                note=(
                    "R maps parent-crystal Cartesian components to child-crystal ones: the "
                    "parent plane normal goes onto the child plane normal, and the parent "
                    "direction onto the child direction."
                ),
            )
        else:
            tag = self.equations.add(
                _s("or_rotation"),
                "=",
                mm.matrix(rotation, digits=digits),
                lead=(
                    "This relationship is not defined by a single plane-and-direction pair, so "
                    "its rotation is taken as constructed:"
                ),
            )
        misorientation = relationship.misorientation()
        axis = np.asarray(misorientation.rotation.axis, dtype=float)
        return ResultStage(
            key="relationship",
            title="Step 3. The orientation relationship as a rotation",
            summary=(
                f"{self.relationship_label} states {statement}. As a rotation from parent to "
                f"child crystal axes it is Eq. ({tag})."
            ),
            equations=self.equations.take(),
            metrics=(
                ResultMetric(
                    "Disorientation angle of R", float(misorientation.angle_deg), units="°"
                ),
                ResultMetric(
                    "Disorientation axis",
                    "[" + ", ".join(f"{value:.4f}" for value in axis) + "]",
                    help_text="Cartesian components; the same in both crystal frames.",
                ),
            ),
            explanation=(
                "The triad construction is how PyTex builds a relationship from its statement; "
                "the matrix shown is the relationship's own rotation and agrees with "
                "the product of the two triads to rounding."
            ),
            section="derivation",
        )

    def variant_statement_nominal(self) -> str:
        relationship = self.report.relationship
        parts = []
        if relationship.parallel_planes:
            parent_plane, child_plane = relationship.parallel_planes[0]
            parts.append(
                f"{plane_label(_plane_indices(parent_plane), spec=self.parent_spec)} ∥ "
                f"{plane_label(_plane_indices(child_plane), spec=self.child_spec)}"
            )
        if relationship.parallel_directions:
            parent_direction, child_direction = relationship.parallel_directions[0]
            parts.append(
                f"{direction_label(_direction_indices(parent_direction), spec=self.parent_spec)}"
                " ∥ "
                f"{direction_label(_direction_indices(child_direction), spec=self.child_spec)}"
            )
        return " with ".join(parts) if parts else "a rotation with no stated parallelism"

    def _triads(self) -> Any:
        """The two triads of the statement, when they reproduce R; otherwise ``None``."""

        relationship = self.report.relationship
        if len(relationship.parallel_planes) != 1 or len(relationship.parallel_directions) != 1:
            return None
        (parent_plane, child_plane) = relationship.parallel_planes[0]
        (parent_direction, child_direction) = relationship.parallel_directions[0]

        def triad(normal: Any, direction: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            n = np.asarray(normal, dtype=float)
            n = n / np.linalg.norm(n)
            d = np.asarray(direction, dtype=float)
            d = d / np.linalg.norm(d)
            x = d - float(d @ n) * n
            x = x / np.linalg.norm(x)
            y = np.cross(n, x)
            return n, d, np.stack([x, y, n], axis=1)

        parent_triad = triad(parent_plane.normal, parent_direction.unit_vector)
        child_triad = triad(child_plane.normal, child_direction.unit_vector)
        rotation = np.asarray(relationship.parent_to_child_rotation.as_matrix(), dtype=float)
        if not np.allclose(child_triad[2] @ parent_triad[2].T, rotation, atol=1e-9):
            return None
        return parent_triad, child_triad

    def _measured_stage(self) -> ResultStage:
        digits = self.digits
        measured = self.report.measured_rotation
        angle = float(np.degrees(rotation_angle_from_matrix_rad(measured)))
        trace = float(np.trace(measured))
        tag = self.equations.add(
            _s("measured_rotation"),
            "=",
            _transpose(_s("orientation_child")),
            _s("orientation_parent"),
            "=",
            mm.matrix(measured, digits=digits),
            lead=(
                "A parent-crystal vector u reaches the specimen as P u; the child reads that "
                "specimen vector back in its own axes as Cᵀ P u. The measured parent-to-child "
                "rotation is therefore"
            ),
        )
        self.equations.add(
            _angle_of(_s("measured_rotation")),
            "=",
            mm.function(
                "arccos",
                mm.frac(
                    mm.row(mm.function("tr", _s("measured_rotation")), mm.mo("−"), mm.mn("1")),
                    mm.mn("2"),
                ),
            ),
            "=",
            mm.function(
                "arccos",
                mm.frac(mm.row(mm.number(trace, digits), mm.mo("−"), mm.mn("1")), mm.mn("2")),
            ),
            "=",
            _degrees(angle, digits),
            lead="Its rotation angle, read from the trace, is",
            note=(
                "This angle is not yet comparable with the relationship: the axes of each "
                "measured crystal could have been labelled in any of its symmetry-equivalent "
                "settings, and each choice changes this matrix. The comparison in Steps 5 and "
                "6 removes that arbitrariness."
            ),
        )
        return ResultStage(
            key="measured_rotation",
            title="Step 4. The measured parent-to-child rotation",
            summary=(
                f"The measured rotation {symbol_text('measured_rotation')} = CᵀP is Eq. ({tag}); "
                f"as it stands it turns by {angle:.2f}°."
            ),
            equations=self.equations.take(),
            section="derivation",
        )

    def _variants_stage(self) -> ResultStage:
        report = self.report
        count = len(report.variants)
        gp = len(report.parent_operators)
        gc = len(report.child_operators)
        k = mm.mi("k")
        tag = self.equations.add(
            _variant_k("k"),
            "=",
            _s("or_rotation"),
            _transpose(_parent_operator_k("k")),
            mm.mo(","),
            mm.space(1),
            _parent_operator_k("k"),
            "∈",
            _s("parent_point_group"),
            mm.mo(","),
            mm.space(1),
            k,
            "=",
            mm.mn("1"),
            mm.mo(","),
            mm.mo("…"),
            mm.mo(","),
            mm.mn(str(count)),
            lead=(
                "Every rotation of the parent point group gives the parent an equivalent "
                "description, and carries the relationship onto a different but equivalent "
                "child orientation — a variant:"
            ),
            note=(
                f"Of the {gp} operators of {symbol_text('parent_point_group')}, those whose "
                f"variants differ only by an operator of {symbol_text('child_point_group')} "
                f"({gc} operators) describe the same child crystal and are counted once, "
                f"leaving {count} variants."
            ),
        )
        packets = self.packets()
        rows = []
        for position, variant in enumerate(report.variants):
            rotation = variant.parent_to_child_rotation
            rows.append(
                {
                    "variant": int(variant.variant_index),
                    "operator_index": int(variant.parent_operator_index) + 1,
                    "operator": self.operator_text(variant.parent_symmetry_operator, parent=True),
                    "packet": packets[position],
                    "statement": self.variant_statement(variant),
                    "angle_deg": float(rotation.angle_deg),
                }
            )
        return ResultStage(
            key="variants",
            title="Step 5. The variants of the relationship",
            summary=(
                f"{self.relationship_label} has {count} variants (Eq. ({tag})), in "
                f"{max(packets)} packet{'s' if max(packets) != 1 else ''} of variants sharing a "
                "parent plane."
            ),
            equations=self.equations.take(),
            metrics=(
                ResultMetric(f"Operators in {symbol_text('parent_point_group')}", gp),
                ResultMetric(f"Operators in {symbol_text('child_point_group')}", gc),
                ResultMetric("Distinct variants", count),
            ),
            table=ResultTable(
                columns=(
                    Column("variant", "Variant k", numeric=True),
                    Column(
                        "operator_index",
                        "Operator in Gp",
                        numeric=True,
                        help_text="One-based position of S_p,k in the parent point group.",
                    ),
                    Column("operator", "S_p,k as a rotation"),
                    Column("packet", "Packet", numeric=True),
                    Column(
                        "statement",
                        "Parallelism realized",
                        help_text="The relationship's plane and direction pair, carried by S_p,k.",
                    ),
                    Column(
                        "angle_deg",
                        "Rotation angle of V_k",
                        units="°",
                        numeric=True,
                        digits=3,
                    ),
                ),
                rows=tuple(rows),
                caption=(
                    "Each variant, the parent operator that generates it, and what it holds "
                    "parallel."
                ),
            ),
            section="derivation",
        )

    def _misfit_stage(self) -> ResultStage:
        report = self.report
        k = "k"
        sc = _s("child_operator")
        self.equations.add(
            _predicted_child(k),
            "=",
            _s("orientation_parent"),
            _transpose(_variant_k(k)),
            lead=(
                "Variant k predicts the child orientation from the measured parent (inverting "
                "C = P Vᵀ):"
            ),
        )
        self.equations.add(
            _transpose(_predicted_child(k)),
            _s("orientation_child"),
            sc,
            "=",
            _variant_k(k),
            _transpose(_s("orientation_parent")),
            _s("orientation_child"),
            sc,
            "=",
            _variant_k(k),
            _transpose(_s("measured_rotation")),
            sc,
            lead=(
                "The rotation from the prediction to the measured child, whose axes may be in "
                "any of its equivalent settings C Sc, is"
            ),
        )
        tag = self.equations.add(
            _s("variant_misfit"),
            "=",
            mm.under(
                mm.mi("min", font="upright"),
                mm.row(sc, mm.mo("∈"), _s("child_point_group")),
            ),
            mm.space(0.2),
            mm.function(
                "arccos",
                mm.frac(
                    mm.row(
                        mm.function(
                            "tr", mm.row(_variant_k(k), _transpose(_s("measured_rotation")), sc)
                        ),
                        mm.mo("−"),
                        mm.mn("1"),
                    ),
                    mm.mn("2"),
                ),
            ),
            lead="and the misfit of variant k is the smallest such rotation angle:",
            note=(
                "Evaluated for every variant against every operator of the child point group; "
                "the table lists the minimum for each variant and the operator that attains it."
            ),
        )
        order = np.argsort(report.misfits_deg, kind="stable")
        top = ", ".join(
            f"V{report.variants[int(item)].variant_index} {report.misfits_deg[int(item)]:.2f}°"
            for item in order[: min(3, len(order))]
        )
        return ResultStage(
            key="misfits",
            title="Step 6. The misfit of every variant",
            summary=(
                f"Eq. ({tag}) evaluated for all {len(report.variants)} variants. The three "
                f"smallest misfits are {top}."
            ),
            equations=self.equations.take(),
            table=ResultTable(
                columns=_MISFIT_COLUMNS,
                rows=tuple(self.misfit_rows()),
                caption=(
                    "Misfit of each variant to the measured child; rank 1 is the identification."
                ),
            ),
            section="derivation",
        )

    def _winner_stage(self) -> ResultStage:
        report = self.report
        digits = self.digits
        index = report.variant_index
        parent_operator = report.parent_operator
        child_operator = report.child_operator
        residual = report.residual_matrix
        trace = float(np.trace(residual))
        cosine = 0.5 * (trace - 1.0)
        skew = 0.5 * np.array(
            [
                residual[2, 1] - residual[1, 2],
                residual[0, 2] - residual[2, 0],
                residual[1, 0] - residual[0, 1],
            ]
        )
        sine = float(np.linalg.norm(skew))
        theta = report.misfit_deg
        position = report.best_position
        self.equations.add(
            _parent_operator_k(index),
            "=",
            mm.matrix(parent_operator, digits=digits),
            lead=(
                f"The parent operator that generates V{index} — "
                f"{self.operator_text(parent_operator, parent=True)}, entry "
                f"{int(report.variant.parent_operator_index) + 1} of "
                f"{symbol_text('parent_point_group')} — is"
            ),
        )
        self.equations.add(
            _variant_k(index),
            "=",
            _s("or_rotation"),
            _transpose(_parent_operator_k(index)),
            "=",
            mm.matrix(report.variant_matrix, digits=digits),
            lead="so that the variant rotation is",
        )
        self.equations.add(
            _predicted_child(index),
            "=",
            _s("orientation_parent"),
            _transpose(_variant_k(index)),
            "=",
            mm.matrix(report.predicted_child_matrix, digits=digits),
            lead="and the child orientation it predicts is",
        )
        self.equations.add(
            _child_operator_star(),
            "=",
            mm.matrix(child_operator, digits=digits),
            lead=(
                "The child operator that brings the measured child into the variant's setting — "
                f"{self.operator_text(child_operator, parent=False)}, entry "
                f"{report.child_operator_index + 1} of {symbol_text('child_point_group')} — is"
            ),
        )
        self.equations.add(
            _s("residual_rotation"),
            "=",
            _variant_k(index),
            _transpose(_s("measured_rotation")),
            _child_operator_star(),
            "=",
            mm.matrix(residual, digits=digits),
            lead="The residual rotation between prediction and measurement is then",
        )
        tag = self.equations.add(
            _s("identified_misfit"),
            "=",
            mm.function(
                "arccos",
                mm.frac(
                    mm.row(mm.number(trace, digits + 2), mm.mo("−"), mm.mn("1")),
                    mm.mn("2"),
                ),
            ),
            "=",
            mm.function(
                "atan2",
                mm.row(
                    mm.number(sine, digits + 2),
                    mm.mo(","),
                    mm.space(0.3),
                    mm.number(cosine, digits + 2),
                ),
            ),
            "=",
            _degrees(theta, digits),
            lead=f"Its trace gives the misfit of V{index}:",
            note=(
                "PyTex evaluates the angle as atan2(½‖vec(Δg − Δgᵀ)‖, ½(tr Δg − 1)): the "
                "same angle as the arccos form, but without losing digits when the misfit is "
                "small, where the trace is within θ² of 3."
            ),
        )
        metrics = [
            ResultMetric("Trace of Δg", trace),
            ResultMetric(f"{symbol_text('identified_misfit')}", theta, units="°"),
        ]
        if theta > 1e-6:
            axis = skew / sine if sine > 0 else np.array([0.0, 0.0, 1.0])
            metrics.append(
                ResultMetric(
                    "Residual axis (child crystal frame)",
                    "[" + ", ".join(f"{value:.4f}" for value in axis) + "]",
                    help_text="Axis about which the measured child is turned from the prediction.",
                )
            )
        return ResultStage(
            key="winner",
            title=f"Step 7. Variant V{index}, substituted in full",
            summary=(
                f"V{index} has the smallest misfit, {theta:.4f}° (Eq. ({tag})). The matrices "
                "below are the operators and rotations that produce it."
            ),
            equations=self.equations.take(),
            metrics=tuple(metrics),
            section="derivation",
            status="ok" if report.misfits_deg[position] <= report.tolerance_deg else "warning",
        )

    def _parallelism_stage(self) -> ResultStage:
        report = self.report
        digits = self.digits
        rows = []
        delta = _s("parallelism_deviation")
        for check in report.parallelism_checks:
            is_plane = check.kind == "plane"
            parent_indices = (
                _plane_indices(check.parent) if is_plane else _direction_indices(check.parent)
            )
            child_indices = (
                _plane_indices(check.child) if is_plane else _direction_indices(check.child)
            )
            family = "plane" if is_plane else "direction"
            parent_label = (
                plane_label(parent_indices, spec=self.parent_spec)
                if is_plane
                else direction_label(parent_indices, spec=self.parent_spec)
            )
            child_label = (
                plane_label(child_indices, spec=self.child_spec)
                if is_plane
                else direction_label(child_indices, spec=self.child_spec)
            )
            vector = _s("unit_normal") if is_plane else _s("unit_direction")
            self.equations.add(
                _miller_node(parent_indices, spec=self.parent_spec, family=family),
                mm.mo(":"),
                mm.space(0.4),
                _s("orientation_parent"),
                mm.sub(vector, mm.mi("p")),
                "=",
                mm.vector(check.parent_specimen, digits=digits),
                mm.mo(","),
                mm.space(1.2),
                _miller_node(child_indices, spec=self.child_spec, family=family),
                mm.mo(":"),
                mm.space(0.4),
                _s("orientation_child"),
                mm.sub(vector, mm.mi("c")),
                "=",
                mm.vector(check.child_specimen, digits=digits),
                lead=(
                    f"The variant holds the parent {check.kind} {parent_label} parallel to the "
                    f"child {check.kind} {child_label} (indices of the measured child). Carried "
                    "into the specimen frame by the measured orientations:"
                ),
            )
            dot = float(np.clip(check.parent_specimen @ check.child_specimen, -1.0, 1.0))
            self.equations.add(
                delta,
                "=",
                mm.function("arccos", mm.number(dot, digits + 2)),
                "=",
                _degrees(check.angle_deg, digits),
            )
            rows.append(
                {
                    "kind": check.kind,
                    "parent": parent_label,
                    "child": child_label,
                    "angle_deg": float(check.angle_deg),
                }
            )
        worst = max((float(check.angle_deg) for check in report.parallelism_checks), default=0.0)
        return ResultStage(
            key="parallelism",
            title="Step 8. Check: the variant's parallelism on the measured grains",
            summary=(
                "The identified variant's own plane and direction pair, carried into the "
                "specimen frame through the two measured orientations, are parallel to within "
                f"{worst:.3f}°. This is the metallurgist's check: no symmetry bookkeeping, just "
                "the two measured crystals."
            ),
            equations=self.equations.take(),
            table=ResultTable(
                columns=(
                    Column("kind", "Kind"),
                    Column("parent", "Parent member"),
                    Column("child", "Child member (measured indices)"),
                    Column("angle_deg", "Deviation δ", units="°", numeric=True, digits=3),
                ),
                rows=tuple(rows),
            ),
            explanation=(
                "Each deviation can be no larger than the misfit θ* of Step 7: a rotation by θ* "
                "moves no vector by more than θ*. Equality holds when the vector is normal to "
                "the residual axis."
            ),
            section="derivation",
        )

    def _decision_stage(self) -> ResultStage:
        report = self.report
        digits = self.digits
        tau = _s("assignment_tolerance")
        theta_star = _s("identified_misfit")
        verdict = report.verdict
        relation = "≤" if report.misfit_deg <= report.tolerance_deg else ">"
        items: list[mm.Node | str] = [
            theta_star,
            "=",
            _degrees(report.misfit_deg, digits),
            relation,
            tau,
            "=",
            _degrees(report.tolerance_deg, digits),
        ]
        runner_up = report.runner_up
        if runner_up is not None and report.runner_up_misfit_deg is not None:
            runner_relation = "≤" if report.runner_up_misfit_deg <= report.tolerance_deg else ">"
            items += [
                mm.mo(","),
                mm.space(1.5),
                _misfit_k(runner_up.variant_index),
                "=",
                _degrees(report.runner_up_misfit_deg, digits),
                runner_relation,
                tau,
            ]
        tag = self.equations.add(
            *items,
            lead=(
                "A variant is assigned when its misfit is inside the tolerance and no other "
                "variant's is:"
            ),
        )
        conclusion = {
            "assigned": (
                f"V{report.variant_index} is inside the tolerance and the runner-up is not, so "
                f"the child grain is assigned to V{report.variant_index}."
            ),
            "ambiguous": (
                "More than one variant is inside the tolerance; this pair cannot decide "
                "between them."
            ),
            "no_fit": (
                "No variant is inside the tolerance; the pair is not explained by this "
                "relationship."
            ),
        }[verdict]
        return ResultStage(
            key="decision",
            title="Step 9. Decision",
            summary=f"Eq. ({tag}). {conclusion}",
            equations=self.equations.take(),
            metrics=(
                ResultMetric("Verdict", _VERDICT_TEXT[verdict]),
                ResultMetric(
                    "Margin to the runner-up",
                    report.margin_deg if report.margin_deg is not None else "—",
                    units="°" if report.margin_deg is not None else None,
                ),
            ),
            status="ok" if verdict == "assigned" else "warning",
            section="derivation",
        )


def _misfit_figure(context: _Context) -> ResultFigure:
    report = context.report
    labels = [f"V{variant.variant_index}" for variant in report.variants]
    misfits = np.asarray(report.misfits_deg, dtype=float)
    best = report.best_position
    tolerance = report.tolerance_deg

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(labels))
        colors = [COLORS["accent"] if index == best else COLORS["guide"] for index in positions]
        axes.bar(positions, np.maximum(misfits, 0.05), color=colors, width=0.7)
        axes.axhline(
            tolerance,
            color=COLORS["warning"],
            lw=0.9,
            ls="--",
            label=f"Assignment tolerance {tolerance:g}°",
        )
        axes.text(
            best,
            max(float(misfits[best]), 0.05) + 1.0,
            f"{misfits[best]:.2f}°",
            ha="center",
            va="bottom",
            fontsize=7,
            color=COLORS["accent"],
        )
        axes.set_xticks(
            positions, labels=labels, fontsize=7, rotation=90 if len(labels) > 12 else 0
        )
        axes.set_xlabel("Variant")
        axes.set_ylabel("Misfit to the measured child (°)")
        axes.set_ylim(0.0, max(float(np.max(misfits)) * 1.08, tolerance * 1.5))
        axes.legend(loc="upper right", frameon=False, fontsize=7)

    margin = report.margin_deg
    return render_figure(
        draw,
        key="variant_misfits",
        title="Misfit of every variant",
        caption=(
            "The smallest rotation from each variant's predicted child orientation to the "
            "measured one, minimized over the child point group (Step 6). The identified "
            "variant is highlighted; the dashed line is the assignment tolerance."
        ),
        interpretation=(
            f"V{report.variant_index} stands alone below the tolerance"
            + (f", {margin:.2f}° ahead of the next variant." if margin is not None else ".")
            if report.verdict == "assigned"
            else "No single variant stands alone below the tolerance; read the warning above."
        ),
        height_in=3.2,
    )


# -- examples -------------------------------------------------------------------------

REGISTRY.add_examples(
    (
        ExampleScenario(
            id="variants.example.identify_burgers",
            title="Which Burgers variant is this alpha grain?",
            panel="variants",
            summary="A beta-zirconium parent and an alpha child, identified step by step.",
            teaches=(
                "The two Euler triples are an exact Burgers pair, so the derivation should end at "
                "a misfit of zero, and it does: variant 1, with the runner-up more than ten "
                "degrees behind. Read Step 7 for the two operators that produce it, and Step 8 "
                "for the check a metallurgist trusts most — the (110) of the parent and the "
                "basal plane of the child, measured parallel on the two grains."
            ),
            operation="variants.identify_variant",
            request={
                "phase": {"builtin": "zr_bcc_beta"},
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
                "euler_convention": "bunge",
                "parent_angle1": 30.0,
                "parent_angle2": 40.0,
                "parent_angle3": 10.0,
                "child_angle1": 167.5709,
                "child_angle2": 58.2280,
                "child_angle3": 0.9653,
                "tolerance_deg": 5.0,
            },
        ),
        ExampleScenario(
            id="variants.example.identify_ks",
            title="A martensite lath and its parent austenite",
            panel="variants",
            summary="Kurdjumov-Sachs: which of 24 variants a ferrite grain is.",
            teaches=(
                "Twenty-four candidates, and the misfit figure shows one standing alone. The "
                "pair is exact Kurdjumov-Sachs variant 1, so the misfit is zero and the "
                "runner-up sits at 10.53 degrees — the smallest angle between two "
                "Kurdjumov-Sachs variants (Morito et al.). That gap is why a five-degree "
                "tolerance can never admit two variants at once."
            ),
            operation="variants.identify_variant",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "euler_convention": "bunge",
                "parent_angle1": 30.0,
                "parent_angle2": 40.0,
                "parent_angle3": 10.0,
                "child_angle1": 45.2774,
                "child_angle2": 34.9979,
                "child_angle3": 316.2482,
                "tolerance_deg": 5.0,
            },
        ),
        ExampleScenario(
            id="variants.example.identify_gt_as_ks",
            title="A real transformation is not exactly the nominal one",
            panel="variants",
            summary="A Greninger-Troiano child tested against Kurdjumov-Sachs.",
            teaches=(
                "The child here obeys Greninger-Troiano exactly, but is tested against "
                "Kurdjumov-Sachs, which is how measured martensite usually behaves: close to "
                "the nominal relationship, not on it. The variant is still identified without "
                "ambiguity, but now at a misfit of about 2.4 degrees — the separation of the "
                "two relationships — and Step 8 shows the parallel planes are no longer exactly "
                "parallel on the grains. A nonzero misfit with a clear margin is the normal, "
                "healthy result on real data."
            ),
            operation="variants.identify_variant",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "euler_convention": "bunge",
                "parent_angle1": 30.0,
                "parent_angle2": 40.0,
                "parent_angle3": 10.0,
                "child_angle1": 126.5536,
                "child_angle2": 73.6774,
                "child_angle3": 352.9432,
                "tolerance_deg": 5.0,
            },
        ),
    )
)
