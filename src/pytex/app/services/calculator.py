"""The crystallographic calculator: the questions a researcher asks daily.

What it does
    Interplanar and interdirection angles, plane-direction geometry, symmetry
    families, zone axes, d-spacing tables, and the same geometry *between* two
    phases of arbitrary space group once their relative orientation is stated.

When to use it
    This is the "look it up quickly" surface. Every operation takes a phase
    described by cell parameters and a point group — not a file — so it works
    for any material a user can write down, and every result comes back with
    the numbers in a table so it can leave the application as CSV or XLSX.

Conventions, stated once
    Plane angles are angles between plane *normals*, computed through the
    reciprocal metric, with antipodal equivalence, so they lie in [0°, 90°].
    Direction angles are computed through the direct metric and are reported
    with antipodal equivalence by default for the same reason: ``[110]`` and
    ``[-1-10]`` are the same line. Indices are rendered by
    :mod:`pytex.core.notation`, so overbars appear in labels and plain minus
    signs in machine-readable columns.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import PhaseSpec, list_builtin_phases, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesListParameter,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    Parameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.core.hexagonal import (
    direction_uvtw_to_uvw,
    direction_uvw_to_uvtw,
    plane_hkil_to_hkl,
    plane_hkl_to_hkil,
)
from pytex.core.lattice import Phase
from pytex.core.miller import (
    MillerDirectionSet,
    MillerPlaneSet,
    angle_dir_plane_inclination_rad,
    angle_dir_plane_normal_rad,
    reduce_indices,
)
from pytex.core.notation import format_miller_indices
from pytex.core.orientation import Rotation
from pytex.core.point_groups import PointGroup, all_point_group_symbols

__all__ = ["direction_label", "family_label", "phase_parameter", "plane_label"]

_CITATION_ITA = "International Tables for Crystallography, Volume A (Space-Group Symmetry), 6th ed."
_CITATION_CULLITY = "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3."
_CITATION_FRANK = (
    "Frank, Acta Crystallogr. 18 (1965) 862 (Miller-Bravais indices as four-dimensional vectors)."
)
_CITATION_OTTE_CROCKER = (
    "Otte & Crocker, Phys. Status Solidi 9 (1965) 441 "
    "(crystallographic formulae for hexagonal lattices)."
)


# --------------------------------------------------------------------------
# Shared parameters and labelling
# --------------------------------------------------------------------------


def phase_parameter(
    name: str = "phase",
    *,
    label: str = "Phase",
    help_text: str | None = None,
    builtin: str | None = None,
) -> ObjectParameter:
    """Return the standard phase-selection parameter.

    Every operation that needs a material takes exactly this parameter, so the
    frontend renders one phase picker everywhere and a phase chosen in one tab
    can be carried into another without translation.

    ``builtin`` names the phase the picker should start on. Give it wherever the
    operation relates two phases: the picker otherwise starts both on the first
    entry in the catalogue, and an operation that needs two *different* phases
    then refuses the very first press of its own button.
    """

    phase_help = help_text or (
        "Choose a built-in phase or enter six cell parameters and a point group. "
        "The point group sets the symmetry; the space-group symbol, if given, sets "
        "which reflections are systematically absent."
    )
    phase_help += (
        " You may instead load a .cif file; PyTex imports its lattice, symmetry, space group "
        "and atomic sites through the canonical CIF phase constructor."
    )
    return ObjectParameter(
        name=name,
        label=label,
        help_text=phase_help,
        editor="phase",
        default=None if builtin is None else {"builtin": builtin},
    )


def _second_phase_parameter(builtin: str | None = None) -> ObjectParameter:
    return ObjectParameter(
        name="other_phase",
        label="Second phase",
        help_text=(
            "The other crystal. Angles between two phases are only defined once their "
            "relative orientation is stated, which is what the rotation below does."
        ),
        editor="phase",
        required=False,
        default=None if builtin is None else {"builtin": builtin},
    )


_ANGLE_DIGITS = 4


def _four_index(spec: PhaseSpec) -> bool:
    return spec.uses_miller_bravais


def plane_label(indices: Sequence[int], *, spec: PhaseSpec, style: str = "plain") -> str:
    """Render one plane in the notation natural to its phase.

    Hexagonal and trigonal phases get four-index Miller-Bravais labels, because
    that is what the literature and every hexagonal-symmetry argument uses; all
    other systems get three.
    """

    from pytex.core.miller import plane_hkl_to_hkil_array

    values = tuple(int(value) for value in indices)
    if _four_index(spec):
        values = tuple(int(value) for value in plane_hkl_to_hkil_array(np.asarray(values))[0])
    return format_miller_indices(values, family="plane", style=style, scope="specific")


def direction_label(indices: Sequence[int], *, spec: PhaseSpec, style: str = "plain") -> str:
    """Render one direction in the notation natural to its phase."""

    from pytex.core.miller import direction_uvw_to_uvtw_array

    values = tuple(int(value) for value in indices)
    if _four_index(spec):
        values = tuple(int(value) for value in direction_uvw_to_uvtw_array(np.asarray(values))[0])
    return format_miller_indices(values, family="direction", style=style, scope="specific")


def family_label(indices: Sequence[int], *, spec: PhaseSpec, family: str) -> str:
    """Render one symmetry family in the notation natural to its phase.

    ``family`` is ``"plane"`` for ``{hkl}`` or ``"direction"`` for
    ``<uvw>``. Hexagonal and trigonal phases get four-index labels, matching
    :func:`plane_label` and :func:`direction_label` for the specific forms.
    """

    from pytex.core.miller import direction_uvw_to_uvtw_array, plane_hkl_to_hkil_array

    values = tuple(int(value) for value in indices)
    if _four_index(spec):
        converter = plane_hkl_to_hkil_array if family == "plane" else direction_uvw_to_uvtw_array
        values = tuple(int(value) for value in converter(np.asarray(values))[0])
    return format_miller_indices(values, family=family, style="plain", scope="family")


def _require_rows(
    rows: tuple[tuple[int, ...], ...] | None, *, field: str, what: str
) -> tuple[tuple[int, ...], ...]:
    if not rows:
        raise InvalidInputError(
            f"Give at least one {what}.",
            field=field,
            hint=f"Enter one {what} per row, for example '1 1 1'.",
        )
    return rows


def _plane_set(phase: Phase, rows: Sequence[Sequence[int]] | np.ndarray) -> MillerPlaneSet:
    return MillerPlaneSet.from_hkl(np.asarray(rows, dtype=int), phase=phase)


def _direction_set(phase: Phase, rows: Sequence[Sequence[int]] | np.ndarray) -> MillerDirectionSet:
    return MillerDirectionSet.from_uvw(np.asarray(rows, dtype=int), phase=phase)


def _plane_family_representatives(phase: Phase, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Group index rows into symmetry families, vectorized over the whole grid.

    Returns one representative row per input and the size of that input's
    antipodal orbit. Two subtleties are handled here rather than left to the
    caller:

    - The library's orbit enumeration reduces indices to lowest terms, so the
      orbit of (222) is reported as the orbit of (111). A d-spacing table must
      keep those apart, because their spacings differ by a factor of two, so the
      common divisor is stripped before grouping and multiplied back afterwards.
    - The representative is the orbit member with the largest packed key, which
      is a deterministic choice that lands on the conventional form — (311)
      rather than (1 -1 3) — without a special case per crystal system.
    """

    indices = np.asarray(rows, dtype=int)
    divisors = np.gcd.reduce(np.abs(indices), axis=1)
    divisors[divisors == 0] = 1
    reduced = indices // divisors[:, np.newaxis]
    members, mask = _plane_set(phase, reduced).symmetry_equivalent_indices()
    members = np.asarray(members, dtype=int)
    valid = np.asarray(mask, dtype=bool)
    span = int(np.abs(members).max(initial=1)) + 1
    packed = ((members[..., 0] + span) * (2 * span) + (members[..., 1] + span)) * (2 * span) + (
        members[..., 2] + span
    )
    packed = np.where(valid, packed, -1)
    chosen = np.argmax(packed, axis=1)
    representatives = members[np.arange(members.shape[0]), chosen]
    return representatives * divisors[:, np.newaxis], valid.sum(axis=1)


def _pair_rows(
    left_rows: Sequence[Sequence[int]],
    right_rows: Sequence[Sequence[int]],
    angles_deg: np.ndarray,
    *,
    spec: PhaseSpec,
    other_spec: PhaseSpec | None,
    family: str,
    left_extra: np.ndarray | None = None,
    right_extra: np.ndarray | None = None,
) -> tuple[dict[str, Any], ...]:
    """Flatten an angle matrix into one row per pair, skipping the diagonal.

    A matrix is what a user wants to *look* at; one row per pair is what a user
    wants to *export*. Both come from the same array, so they cannot disagree.
    """

    labeller = plane_label if family == "plane" else direction_label
    right_spec = other_spec or spec
    rows: list[dict[str, Any]] = []
    for i, left in enumerate(left_rows):
        for j, right in enumerate(right_rows):
            if other_spec is None and j <= i:
                continue
            entry: dict[str, Any] = {
                "left": labeller(left, spec=spec),
                "right": labeller(right, spec=right_spec),
                "angle_deg": float(angles_deg[i, j]),
            }
            if left_extra is not None:
                entry["left_d_angstrom"] = float(left_extra[i])
            if right_extra is not None:
                entry["right_d_angstrom"] = float(right_extra[j])
            rows.append(entry)
    return tuple(rows)


def _angle_columns(*, family: str, with_spacing: bool) -> tuple[Column, ...]:
    what = "Plane" if family == "plane" else "Direction"
    columns = [
        Column("left", f"{what} 1"),
        Column("right", f"{what} 2"),
        Column("angle_deg", "Angle", units="°", numeric=True, digits=_ANGLE_DIGITS),
    ]
    if with_spacing:
        columns.extend(
            [
                Column(
                    "left_d_angstrom",
                    "d₁",
                    units="Å",
                    numeric=True,
                    digits=5,
                    help_text="Interplanar spacing of the first plane.",
                ),
                Column("right_d_angstrom", "d₂", units="Å", numeric=True, digits=5),
            ]
        )
    return tuple(columns)


def _matrix_payload(
    left_rows: Sequence[Sequence[int]],
    right_rows: Sequence[Sequence[int]],
    angles_deg: np.ndarray,
    *,
    spec: PhaseSpec,
    other_spec: PhaseSpec | None,
    family: str,
) -> dict[str, Any]:
    labeller = plane_label if family == "plane" else direction_label
    right_spec = other_spec or spec
    return {
        "row_labels": [labeller(row, spec=spec) for row in left_rows],
        "column_labels": [labeller(row, spec=right_spec) for row in right_rows],
        "values_deg": [[float(value) for value in row] for row in angles_deg],
    }


def _rotation_from_request(
    convention: str, values: tuple[float, float, float], angle_deg: float
) -> Rotation:
    """Build the rotation that carries the second crystal into the first.

    Two conventions are offered because researchers hold relative orientations
    in two different ways: Bunge Euler angles, which is how an OR arrives from
    EBSD software, and an axis with an angle, which is how one arrives from a
    misorientation table.
    """

    if convention == "bunge":
        return Rotation.from_bunge_euler(values[0], values[1], values[2], degrees=True)
    axis = np.asarray(values, dtype=float)
    if not np.any(axis):
        raise InvalidInputError(
            "The rotation axis must not be the zero vector.",
            field="rotation_values",
            hint="Give a direction such as 1 1 1.",
        )
    return Rotation.from_axis_angle(axis, math.radians(angle_deg))


#: The first Kurdjumov-Sachs variant, as the rotation this operation consumes.
#:
#: It is the opening state of the cross-phase angle operation, so that the first
#: press of the button answers the question the help text poses — does this
#: relationship really put a {111} of the parent on a {110} of the product —
#: instead of comparing a phase with itself through a null rotation. It does:
#: the top row of that first table is austenite (111) against ferrite (011) at 0°.
#:
#: Two conventions are folded in and both matter. The operation carries the
#: *second* crystal into the first, so the rotation wanted here is the inverse of
#: the variant's parent-to-child rotation; the inverse is written by negating the
#: axis rather than the angle, because a negative angle in a form field reads as
#: a mistake. And a variant rotation is used rather than the disorientation,
#: because the disorientation is a symmetry-reduced representative that need not
#: map (111) onto (011) — it maps *some* member of each family onto the other,
#: which is the right answer to a different question.
#:
#: The numbers are not transcribed from the literature: they are what
#: :meth:`OrientationRelationship.from_kurdjumov_sachs_correspondence` computes
#: for austenite and ferrite, and
#: ``test_app_calculator.py::test_kurdjumov_sachs_default_matches_the_computed_relationship``
#: recomputes them and fails if the two ever part company.
KS_VARIANT_AXIS: tuple[float, float, float] = (-0.177620, 0.177620, -0.967937)
KS_VARIANT_ANGLE_DEG: float = 42.847760


def _rotation_parameters(
    *,
    convention: str = "bunge",
    components: tuple[float, float, float] = (0.0, 0.0, 0.0),
    angle_deg: float = 0.0,
) -> tuple[Parameter, ...]:
    """Return the relative-orientation parameter group.

    The defaults are per-operation because "no rotation" is a sensible opening
    state for some cross-phase questions and a meaningless one for others.
    """

    return (
        ChoiceParameter(
            name="rotation_convention",
            label="Relative orientation",
            help_text=(
                "How the rotation carrying the second crystal into the first is stated. "
                "Bunge Euler angles are what EBSD software reports; axis and angle is what "
                "a misorientation or orientation-relationship table reports."
            ),
            options=(
                ("bunge", "Bunge Euler (φ₁, Φ, φ₂)", "Three ZXZ Euler angles in degrees."),
                ("axis_angle", "Axis and angle", "A rotation axis in crystal 1, plus an angle."),
            ),
            default=convention,
            group="Relative orientation",
        ),
        NumberParameter(
            name="rotation_1",
            label="φ₁ or axis x",
            help_text="First Euler angle in degrees, or the x component of the rotation axis.",
            default=components[0],
            group="Relative orientation",
        ),
        NumberParameter(
            name="rotation_2",
            label="Φ or axis y",
            help_text="Second Euler angle in degrees, or the y component of the rotation axis.",
            default=components[1],
            group="Relative orientation",
        ),
        NumberParameter(
            name="rotation_3",
            label="φ₂ or axis z",
            help_text="Third Euler angle in degrees, or the z component of the rotation axis.",
            default=components[2],
            group="Relative orientation",
        ),
        NumberParameter(
            name="rotation_angle_deg",
            label="Rotation angle",
            help_text=(
                "Rotation angle for the axis-and-angle convention. Ignored for Euler angles."
            ),
            units="°",
            default=angle_deg,
            group="Relative orientation",
        ),
    )


# --------------------------------------------------------------------------
# Catalogue and phase summary
# --------------------------------------------------------------------------


@REGISTRY.operation(
    "calc.catalog",
    title="Phase catalogue",
    summary="Built-in phases and the 32 crystallographic point groups.",
    help_text=(
        "Everything the phase picker needs before a phase is chosen. Built-in phases carry "
        "literal, cited cell parameters and a full atomic basis, so they work with no optional "
        "dependency installed. Point groups are listed with their crystal system and the metric "
        "constraints that system imposes, which is what the phase editor validates against."
    ),
    returns="A table of built-in phases, plus the point-group list under `data.point_groups`.",
    panel="calculator",
    citations=(_CITATION_ITA,),
    tags=("phase", "catalogue", "point group", "space group", "material"),
)
def _catalog(_: dict[str, Any]) -> dict[str, Any]:
    entries = list_builtin_phases()
    point_groups = []
    for symbol in all_point_group_symbols():
        group = PointGroup.from_symbol(symbol)
        point_groups.append(
            {
                "symbol": symbol,
                "crystal_system": str(group.crystal_system),
                "order": int(group.order),
                "laue_class": str(group.laue_class_symbol),
                "centrosymmetric": bool(group.is_centrosymmetric),
            }
        )
    table = ResultTable(
        columns=(
            Column("id", "Identifier"),
            Column("name", "Phase"),
            Column("crystal_system", "System"),
            Column("point_group", "Point group"),
            Column("space_group_symbol", "Space group"),
            Column("a", "a", units="Å", numeric=True, digits=5),
            Column("c", "c", units="Å", numeric=True, digits=5),
            Column("source", "Source"),
        ),
        rows=tuple(
            {
                "id": entry["id"],
                "name": entry["name"],
                "crystal_system": entry["crystal_system"],
                "point_group": entry["point_group"],
                "space_group_symbol": entry.get("space_group_symbol", "—"),
                "a": entry["a"],
                "c": entry["c"],
                "source": entry.get("source", ""),
            }
            for entry in entries
        ),
    )
    result = AppResult(
        title="Phase catalogue",
        summary=(
            f"{len(entries)} built-in phases are available, spanning cubic, hexagonal, "
            "trigonal and orthorhombic symmetry. Any other material can be entered directly "
            "as six cell parameters and one of the 32 crystallographic point groups."
        ),
        table=table,
        data={"phases": list(entries), "point_groups": point_groups},
        citations=(_CITATION_ITA,),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.phase_summary",
    title="Phase summary",
    summary="Cell metric, reciprocal cell, symmetry order, and centring absences.",
    help_text=(
        "The derived quantities that follow from the six cell parameters alone: cell volume, "
        "the direct and reciprocal metric tensors, the reciprocal cell parameters, and the "
        "order of the point group.\n\n"
        "Read the reciprocal parameters when you need to reason about a diffraction pattern's "
        "scale, and the metric tensor when you want to see *why* an angle came out as it did — "
        "every interplanar angle in this calculator is a quadratic form in one of these two "
        "matrices."
    ),
    parameters=(phase_parameter(),),
    returns="A one-row table of derived cell quantities; matrices under `data`.",
    panel="calculator",
    citations=(_CITATION_ITA,),
    tags=("cell", "volume", "reciprocal", "metric tensor", "symmetry"),
)
def _phase_summary(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    lattice = phase.lattice
    metric = np.asarray(lattice.metric_tensor(), dtype=float)
    reciprocal_metric = np.asarray(lattice.reciprocal_metric_tensor(), dtype=float)
    volume = lattice.volume_angstrom3()
    # Basis vectors are the columns of the basis matrix, so the reciprocal axes
    # are its columns and the lengths are column norms.
    reciprocal_basis = np.asarray(lattice.reciprocal_basis().matrix, dtype=float).T
    reciprocal_lengths = np.linalg.norm(reciprocal_basis, axis=1)
    reciprocal_angles = [
        math.degrees(
            math.acos(
                float(
                    np.clip(
                        np.dot(reciprocal_basis[i], reciprocal_basis[j])
                        / (reciprocal_lengths[i] * reciprocal_lengths[j]),
                        -1.0,
                        1.0,
                    )
                )
            )
        )
        for i, j in ((1, 2), (0, 2), (0, 1))
    ]
    group = PointGroup.from_symbol(spec.point_group)
    notes: list[str] = []
    if spec.space_group_symbol is None:
        notes.append(
            "No space group was given, so systematic absences cannot be applied. Reflection "
            "tables will list every index, including ones a real pattern does not show."
        )
    if not spec.has_structure:
        notes.append(
            "No atomic basis was given, so structure factors and intensities are unavailable; "
            "lattice geometry is unaffected."
        )
    table = ResultTable(
        columns=(
            Column("quantity", "Quantity"),
            Column("value", "Value", numeric=True, digits=6),
            Column("units", "Units"),
        ),
        rows=(
            {"quantity": "a", "value": spec.a, "units": "Å"},
            {"quantity": "b", "value": spec.b, "units": "Å"},
            {"quantity": "c", "value": spec.c, "units": "Å"},
            {"quantity": "alpha", "value": spec.alpha, "units": "°"},
            {"quantity": "beta", "value": spec.beta, "units": "°"},
            {"quantity": "gamma", "value": spec.gamma, "units": "°"},
            {"quantity": "cell volume", "value": volume, "units": "Å³"},
            {"quantity": "a*", "value": float(reciprocal_lengths[0]), "units": "Å⁻¹"},
            {"quantity": "b*", "value": float(reciprocal_lengths[1]), "units": "Å⁻¹"},
            {"quantity": "c*", "value": float(reciprocal_lengths[2]), "units": "Å⁻¹"},
            {"quantity": "alpha*", "value": reciprocal_angles[0], "units": "°"},
            {"quantity": "beta*", "value": reciprocal_angles[1], "units": "°"},
            {"quantity": "gamma*", "value": reciprocal_angles[2], "units": "°"},
            {"quantity": "point-group order", "value": float(group.order), "units": ""},
        ),
    )
    result = AppResult(
        title=f"{spec.name}: cell summary",
        summary=(
            f"{spec.name} is {spec.crystal_system} with point group {spec.point_group} "
            f"(order {group.order}, Laue class {group.laue_class_symbol}). The unit cell has "
            f"volume {volume:.4f} Å³. Reciprocal parameters are stated with the "
            "crystallographic convention a*·a = 1, so |a*| is in Å⁻¹ and reciprocal-lattice "
            "vector magnitudes equal 1/d."
        ),
        table=table,
        data={
            "phase": spec.to_json(),
            "cell_volume_angstrom3": volume,
            "metric_tensor": metric.tolist(),
            "reciprocal_metric_tensor": reciprocal_metric.tolist(),
            "reciprocal_basis_cartesian": reciprocal_basis.tolist(),
            "point_group": {
                "symbol": spec.point_group,
                "order": int(group.order),
                "laue_class": str(group.laue_class_symbol),
                "crystal_system": str(group.crystal_system),
                "centrosymmetric": bool(group.is_centrosymmetric),
            },
        },
        inputs={"phase": spec.to_json()},
        notes=notes,
        citations=(_CITATION_ITA,),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Angles
# --------------------------------------------------------------------------


@REGISTRY.operation(
    "calc.plane_angles",
    title="Interplanar angles",
    summary="Angles between plane normals, with interplanar spacings.",
    help_text=(
        "Enter any number of planes; the result is the full angle matrix between them, plus "
        "one row per pair for export.\n\n"
        "The angle reported is the angle between plane **normals**, evaluated through the "
        "reciprocal metric tensor, with antipodal equivalence — a normal and its opposite "
        "describe the same plane — so every value lies in [0°, 90°]. For a cubic phase the "
        "(100)/(110) angle is exactly 45° whatever the lattice parameter, which makes it a "
        "quick check that the phase is set up as intended."
    ),
    parameters=(
        phase_parameter(),
        IndicesListParameter(
            name="planes",
            label="Planes (hkl)",
            help_text=(
                "One plane per row, for example '1 1 1'. Use plain minus signs for negative "
                "indices; the labels come back with overbars."
            ),
            default=((1, 0, 0), (1, 1, 0), (1, 1, 1)),
        ),
        IndicesListParameter(
            name="against",
            label="Against planes (optional)",
            help_text=(
                "Leave empty to get the angles among the planes above. Fill it in to get a "
                "rectangular matrix of the first set against this one."
            ),
            required=False,
        ),
    ),
    returns="Pairwise angle rows; the full matrix under `data.matrix`.",
    panel="calculator",
    citations=(_CITATION_ITA, _CITATION_CULLITY),
    tags=("angle", "interplanar", "plane", "normal", "hkl"),
)
def _plane_angles(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    left_rows = _require_rows(request["planes"], field="planes", what="plane")
    right_rows = request.get("against") or left_rows
    left = _plane_set(phase, left_rows)
    right = _plane_set(phase, right_rows)
    angles_deg = np.degrees(np.asarray(left.angle_matrix_rad(right), dtype=float))
    left_d = np.asarray(left.d_spacings_angstrom(), dtype=float)
    right_d = np.asarray(right.d_spacings_angstrom(), dtype=float)
    symmetric = request.get("against") is None
    rows = _pair_rows(
        left_rows,
        right_rows,
        angles_deg,
        spec=spec,
        other_spec=None if symmetric else spec,
        family="plane",
        left_extra=left_d,
        right_extra=right_d,
    )
    if not rows:
        raise InvalidInputError(
            "Give at least two planes, or a second set to compare against.",
            field="planes",
            hint="One plane on its own has no angle to report.",
        )
    smallest = min(rows, key=lambda row: row["angle_deg"])
    result = AppResult(
        title=f"Interplanar angles in {spec.name}",
        summary=(
            f"{len(rows)} plane pairs in {spec.name} ({spec.crystal_system}, "
            f"{spec.point_group}). Angles are between plane normals through the reciprocal "
            "metric, with antipodal equivalence, so they lie in [0°, 90°]. The closest pair is "
            f"{smallest['left']} to {smallest['right']} at {smallest['angle_deg']:.4f}°."
        ),
        table=ResultTable(
            columns=_angle_columns(family="plane", with_spacing=True),
            rows=rows,
            caption=f"Interplanar angles in {spec.name}.",
        ),
        data={
            "matrix": _matrix_payload(
                left_rows,
                right_rows,
                angles_deg,
                spec=spec,
                other_spec=None,
                family="plane",
            ),
            "d_spacings_angstrom": {
                "left": left_d.tolist(),
                "right": right_d.tolist(),
            },
        },
        inputs={
            "phase": spec.to_json(),
            "planes": [list(row) for row in left_rows],
            "against": None if symmetric else [list(row) for row in right_rows],
        },
        citations=(_CITATION_ITA, _CITATION_CULLITY),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.direction_angles",
    title="Interdirection angles",
    summary="Angles between lattice directions, through the direct metric.",
    help_text=(
        "Enter any number of directions; the result is the full angle matrix between them.\n\n"
        "Angles are evaluated through the **direct** metric tensor, which is what makes this "
        "different from the interplanar case: in a non-cubic crystal the angle between [uvw] "
        "and [u'v'w'] is not the angle between the planes with the same indices. Antipodal "
        "equivalence is on by default, because [110] and its reverse are the same line; turn it "
        "off when the sense of the direction matters, as it does for a Burgers vector."
    ),
    parameters=(
        phase_parameter(),
        IndicesListParameter(
            name="directions",
            label="Directions [uvw]",
            help_text="One direction per row, for example '1 1 1'.",
            default=((1, 0, 0), (1, 1, 0), (1, 1, 1)),
        ),
        IndicesListParameter(
            name="against",
            label="Against directions (optional)",
            help_text="Leave empty for angles among the directions above.",
            required=False,
        ),
        BooleanParameter(
            name="antipodal",
            label="Treat reversed directions as equal",
            help_text=(
                "On: angles lie in [0°, 90°] and [uvw] equals [-u-v-w]. Off: angles lie in "
                "[0°, 180°] and the sense of the direction is kept."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns="Pairwise angle rows; the full matrix under `data.matrix`.",
    panel="calculator",
    citations=(_CITATION_ITA,),
    tags=("angle", "direction", "uvw", "zone"),
)
def _direction_angles(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    left_rows = _require_rows(request["directions"], field="directions", what="direction")
    right_rows = request.get("against") or left_rows
    antipodal = bool(request["antipodal"])
    left = _direction_set(phase, left_rows)
    right = _direction_set(phase, right_rows)
    angles_deg = np.degrees(
        np.asarray(left.angle_matrix_rad(right, antipodal=antipodal), dtype=float)
    )
    symmetric = request.get("against") is None
    rows = _pair_rows(
        left_rows,
        right_rows,
        angles_deg,
        spec=spec,
        other_spec=None if symmetric else spec,
        family="direction",
    )
    if not rows:
        raise InvalidInputError(
            "Give at least two directions, or a second set to compare against.",
            field="directions",
            hint="One direction on its own has no angle to report.",
        )
    result = AppResult(
        title=f"Interdirection angles in {spec.name}",
        summary=(
            f"{len(rows)} direction pairs in {spec.name} ({spec.crystal_system}). Angles are "
            "evaluated through the direct metric tensor"
            + (
                ", with antipodal equivalence, so they lie in [0°, 90°]."
                if antipodal
                else ", keeping the sense of each direction, so they lie in [0°, 180°]."
            )
        ),
        table=ResultTable(
            columns=_angle_columns(family="direction", with_spacing=False),
            rows=rows,
            caption=f"Interdirection angles in {spec.name}.",
        ),
        data={
            "matrix": _matrix_payload(
                left_rows,
                right_rows,
                angles_deg,
                spec=spec,
                other_spec=None,
                family="direction",
            ),
            "antipodal": antipodal,
        },
        inputs={
            "phase": spec.to_json(),
            "directions": [list(row) for row in left_rows],
            "against": None if symmetric else [list(row) for row in right_rows],
            "antipodal": antipodal,
        },
        citations=(_CITATION_ITA,),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.plane_direction_angles",
    title="Plane and direction geometry",
    summary="Angle to the plane normal, inclination to the plane, and zone membership.",
    help_text=(
        "For every direction against every plane, three related numbers:\n\n"
        "- the angle to the plane **normal**;\n"
        "- the **inclination** to the plane itself, which is 90° minus the first;\n"
        "- the zone-law value h·u + k·v + l·w, which is zero exactly when the direction lies "
        "in the plane.\n\n"
        "The zone-law column is the one to read when checking whether a direction is a valid "
        "zone axis for a set of reflections: an integer zero is exact, not approximate, so it "
        "settles the question that a rounded angle cannot."
    ),
    parameters=(
        phase_parameter(),
        IndicesListParameter(
            name="planes",
            label="Planes (hkl)",
            help_text="One plane per row.",
            default=((1, 1, 1),),
        ),
        IndicesListParameter(
            name="directions",
            label="Directions [uvw]",
            help_text="One direction per row.",
            default=((1, 1, 0), (1, -1, 0)),
        ),
    ),
    returns="One row per plane-direction pair.",
    panel="calculator",
    citations=(_CITATION_ITA,),
    tags=("zone law", "inclination", "plane", "direction", "Weiss"),
)
def _plane_direction_angles(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    plane_rows = _require_rows(request["planes"], field="planes", what="plane")
    direction_rows = _require_rows(request["directions"], field="directions", what="direction")
    rows: list[dict[str, Any]] = []
    in_zone = 0
    for plane_indices in plane_rows:
        planes = _plane_set(phase, [plane_indices] * len(direction_rows))
        directions = _direction_set(phase, direction_rows)
        normal_angles = np.degrees(
            np.asarray(angle_dir_plane_normal_rad(directions, planes), dtype=float)
        )
        inclinations = np.degrees(
            np.asarray(angle_dir_plane_inclination_rad(directions, planes), dtype=float)
        )
        for index, direction_indices in enumerate(direction_rows):
            zone_value = int(np.dot(np.asarray(plane_indices), np.asarray(direction_indices)))
            in_zone += int(zone_value == 0)
            rows.append(
                {
                    "plane": plane_label(plane_indices, spec=spec),
                    "direction": direction_label(direction_indices, spec=spec),
                    "normal_angle_deg": float(normal_angles[index]),
                    "inclination_deg": float(inclinations[index]),
                    "zone_law": zone_value,
                    "in_zone": zone_value == 0,
                }
            )
    result = AppResult(
        title=f"Plane and direction geometry in {spec.name}",
        summary=(
            f"{len(rows)} plane-direction pairs in {spec.name}. {in_zone} of them satisfy the "
            "zone law h·u + k·v + l·w = 0 exactly, meaning the direction lies in the plane and "
            "the plane belongs to that zone. The normal angle and the inclination are "
            "complementary by construction."
        ),
        table=ResultTable(
            columns=(
                Column("plane", "Plane"),
                Column("direction", "Direction"),
                Column(
                    "normal_angle_deg",
                    "Angle to normal",
                    units="°",
                    numeric=True,
                    digits=_ANGLE_DIGITS,
                ),
                Column(
                    "inclination_deg",
                    "Inclination to plane",
                    units="°",
                    numeric=True,
                    digits=_ANGLE_DIGITS,
                    help_text="90° minus the angle to the normal.",
                ),
                Column(
                    "zone_law",
                    "h·u + k·v + l·w",
                    numeric=True,
                    help_text="Exactly zero when the direction lies in the plane.",
                ),
                Column("in_zone", "In zone"),
            ),
            rows=tuple(rows),
        ),
        data={"in_zone_count": in_zone},
        inputs={
            "phase": spec.to_json(),
            "planes": [list(row) for row in plane_rows],
            "directions": [list(row) for row in direction_rows],
        },
        citations=(_CITATION_ITA,),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Symmetry families, zone axes, spacings
# --------------------------------------------------------------------------


@REGISTRY.operation(
    "calc.symmetry_family",
    title="Symmetry-equivalent set",
    summary="The full {hkl} or ⟨uvw⟩ family, with its multiplicity.",
    help_text=(
        "Applies the point group of the phase to one plane or direction and lists every "
        "distinct member of the resulting family, together with the multiplicity.\n\n"
        "The multiplicity is *not* simply the order of the point group: whenever an operation "
        "maps the plane onto itself the orbit is shorter, which is why {100} has 3 members in "
        "m-3m and not 48.\n\n"
        "Planes are always treated with antipodal equivalence — (hkl) and (-h-k-l) are the same "
        "plane — so {100} is listed as (100), (010), (001). Powder multiplicity counts the two "
        "senses separately and is therefore twice this number, 6 for {100}; the d-spacing table "
        "reports that convention, since it is the one that weights a diffraction peak. "
        "Directions can be treated either way."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="family",
            label="Kind",
            help_text="Whether the indices describe a plane or a direction.",
            options=(
                (
                    "plane",
                    "Plane {hkl}",
                    "Symmetry family of a plane, through the reciprocal basis.",
                ),
                (
                    "direction",
                    "Direction ⟨uvw⟩",
                    "Symmetry family of a lattice direction, through the direct basis.",
                ),
            ),
            default="plane",
        ),
        IndicesParameter(
            name="indices",
            label="Indices",
            help_text="Three indices, for example '1 1 1'.",
            default=(1, 1, 1),
        ),
        BooleanParameter(
            name="antipodal",
            label="Treat reversed directions as equal",
            help_text=(
                "Applies to directions only; planes are always antipodal. Off doubles most "
                "direction families, listing [uvw] and [-u-v-w] separately."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns="One row per family member; multiplicity under `data.multiplicity`.",
    panel="calculator",
    citations=(_CITATION_ITA, _CITATION_CULLITY),
    tags=("family", "multiplicity", "orbit", "symmetry", "equivalent"),
)
def _symmetry_family(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    indices = tuple(request["indices"])
    family = str(request["family"])
    antipodal = bool(request["antipodal"])
    if family == "plane":
        members, mask = _plane_set(phase, [indices]).symmetry_equivalent_indices()
    else:
        members, mask = _direction_set(phase, [indices]).symmetry_equivalent_indices(
            antipodal=antipodal
        )
    valid = np.asarray(members[0])[np.asarray(mask[0], dtype=bool)]
    labeller = plane_label if family == "plane" else direction_label
    rows = tuple(
        {
            "index": position + 1,
            "member": labeller(tuple(int(value) for value in row), spec=spec),
            "h": int(row[0]),
            "k": int(row[1]),
            "l": int(row[2]),
        }
        for position, row in enumerate(valid)
    )
    multiplicity = len(rows)
    group = PointGroup.from_symbol(spec.point_group)
    label = family_label(indices, spec=spec, family=family)
    stabilizer = int(group.order) // max(multiplicity, 1)
    extras: dict[str, Any] = {}
    if family == "plane":
        d_spacing = float(_plane_set(phase, [indices]).d_spacings_angstrom()[0])
        extras["d_spacing_angstrom"] = d_spacing
    result = AppResult(
        title=f"{label} in {spec.name}",
        summary=(
            f"The family {label} has {multiplicity} distinct members under point group "
            f"{spec.point_group} (order {group.order}). The ratio {group.order}/{multiplicity} "
            f"= {stabilizer} is the number of operations that map "
            + ("the plane" if family == "plane" else "the direction")
            + " onto itself, which is why the family is shorter than the group. "
            + (
                "Every member shares the interplanar spacing "
                f"d = {extras['d_spacing_angstrom']:.5f} Å."
                if family == "plane"
                else "Every member is crystallographically indistinguishable from the input."
            )
        ),
        table=ResultTable(
            columns=(
                Column("index", "#", numeric=True),
                Column("member", "Member"),
                Column("h", "h" if family == "plane" else "u", numeric=True),
                Column("k", "k" if family == "plane" else "v", numeric=True),
                Column("l", "l" if family == "plane" else "w", numeric=True),
            ),
            rows=rows,
            caption=f"Members of {label}.",
        ),
        data={
            "family_label": label,
            "multiplicity": multiplicity,
            "point_group_order": int(group.order),
            "self_mapping_operations": stabilizer,
            "members": [[int(value) for value in row] for row in valid],
            **extras,
        },
        inputs={
            "phase": spec.to_json(),
            "family": family,
            "indices": list(indices),
            "antipodal": antipodal,
        },
        citations=(_CITATION_ITA, _CITATION_CULLITY),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.hexagonal_indices",
    title="Hexagonal index converter (3 ↔ 4)",
    summary="Miller (hkl) or [uvw] to Miller-Bravais (hkil) or [UVTW], and back.",
    help_text=(
        "Converts one plane or one direction between the three-index Miller form and the "
        "four-index Miller-Bravais form used for hexagonal and trigonal crystals, in whichever "
        "of the two directions is wanted.\n\n"
        "**Why the four-index form exists.** The hexagonal basal plane has three equivalent "
        "⟨a⟩ axes at 120°, but a three-index basis names only two of them. The three "
        "crystallographically identical close-packed directions therefore get the unrelated "
        "labels [100], [010] and [-1-10], and nothing about those symbols says they are the same "
        "family. Adding the redundant third basal index makes them permutations of each other — "
        "[2-1-10], [-12-10], [-1-120] — so a family can be read off the indices.\n\n"
        "**Planes take the redundant index; directions take a change of basis.** For a plane, "
        "`i = -(h + k)` and nothing else changes, because the reciprocal basis vectors are "
        "already at 120° to one another. For a direction the transformation is real: "
        "`U = (2u - v)/3`, `V = (2v - u)/3`, `T = -(U + V)`, `W = w`, and back through "
        "`u = 2U + V`, `v = 2V + U`, `w = W`. This is why (100) and [100] convert differently, "
        "and why converting a direction as though it were a plane is a silent error rather than "
        "a loud one.\n\n"
        "**Redundancy is checked, not assumed.** A four-index plane must satisfy "
        "`i = -(h + k)` and a four-index direction `U + V + T = 0`; a quadruple that does not is "
        "refused rather than quietly reinterpreted. The three-index direction is carried through "
        "exact rational arithmetic and the denominators cleared, so the answer is always an "
        "integer quadruple describing the same line — [111] becomes [11-23], not a decimal.\n\n"
        "The table is the symmetry family of whichever entity was given, written in both "
        "notations side by side, because seeing the four-index members as permutations of one "
        "another is the point of the notation."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "A hexagonal or trigonal phase. The conversion itself needs only the crystal "
                "system, but the family, the spacing and the angle reported alongside it are "
                "properties of this cell and this point group."
            ),
            builtin="ti_hcp",
        ),
        ChoiceParameter(
            name="kind",
            label="Kind",
            help_text=(
                "Whether the indices describe a plane or a direction. The two convert by "
                "different rules, so this is not a labelling choice."
            ),
            options=(
                (
                    "plane",
                    "Plane (hkl) ↔ (hkil)",
                    "Insert or drop the redundant index i = -(h + k).",
                ),
                (
                    "direction",
                    "Direction [uvw] ↔ [UVTW]",
                    "Change basis between the two- and three-axis basal descriptions.",
                ),
            ),
            default="direction",
        ),
        ChoiceParameter(
            name="input_notation",
            label="Given in",
            help_text=(
                "Which of the two rows below is read. The other is ignored, and the answer is "
                "reported in both notations either way."
            ),
            options=(
                ("three", "Three indices", "Convert the three-index row to four indices."),
                ("four", "Four indices", "Convert the four-index row to three indices."),
            ),
            default="four",
        ),
        IndicesParameter(
            name="three_index",
            label="Three-index row (hkl)",
            help_text=(
                "Read when *Given in* is Three indices. A plane (hkl) or a direction [uvw] — "
                "which one is set by *Kind*, not by this box."
            ),
            default=(1, 1, 1),
            group="Indices",
        ),
        IndicesParameter(
            name="four_index",
            label="Four-index row (hkil)",
            help_text=(
                "Read when *Given in* is Four indices. A plane (hkil) with i = -(h + k), or a "
                "direction [UVTW] with U + V + T = 0."
            ),
            default=(1, 1, -2, 0),
            width=4,
            group="Indices",
        ),
        BooleanParameter(
            name="antipodal",
            label="Treat reversed directions as equal",
            help_text=(
                "Applies to directions only; planes are always antipodal. On, the family lists "
                "one representative per reversed pair, which is the shorter list. Off lists both "
                "senses, and is the setting under which the four-index members are literally "
                "permutations of one another rather than permutations up to sign."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns=(
        "One row per member of the family, in both notations; the converted row itself under "
        "`data.three_index` and `data.four_index`."
    ),
    panel="calculator",
    citations=(
        _CITATION_ITA,
        _CITATION_FRANK,
        _CITATION_OTTE_CROCKER,
    ),
    tags=(
        "hexagonal",
        "Miller-Bravais",
        "four-index",
        "hkil",
        "UVTW",
        "conversion",
        "notation",
        "hcp",
    ),
)
def _hexagonal_indices(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    if not spec.uses_miller_bravais:
        raise InvalidInputError(
            f"{spec.name} is {spec.point_group}, which is neither hexagonal nor trigonal, so "
            "the four-index Miller-Bravais notation does not apply to it.",
            field="phase",
            hint=(
                "Miller-Bravais indices exist to make the three equivalent basal ⟨a⟩ axes "
                "permutations of one another, which only happens under hexagonal symmetry. "
                "Choose a hexagonal phase such as titanium, zirconium or magnesium, or load a "
                "hexagonal CIF."
            ),
        )
    kind = str(request["kind"])
    given = str(request["input_notation"])
    is_plane = kind == "plane"

    widen = plane_hkl_to_hkil if is_plane else direction_uvw_to_uvtw
    if given == "three":
        three = tuple(int(value) for value in request["three_index"])
        four = tuple(int(value) for value in widen(three))
    else:
        four = tuple(int(value) for value in request["four_index"])
        narrow = plane_hkil_to_hkl if is_plane else direction_uvtw_to_uvw
        try:
            three = tuple(int(value) for value in narrow(four))
        except ValueError as error:
            constraint = "i = -(h + k)" if is_plane else "U + V + T = 0"
            raise InvalidInputError(
                f"That four-index {kind} is not self-consistent: {error}",
                field="four_index",
                hint=(
                    f"The third index is redundant, not free: it is fixed by {constraint}. For "
                    + (
                        "the first-order prism plane the row is 1 0 -1 0, not 1 0 1 0."
                        if is_plane
                        else "the close-packed direction the row is 2 -1 -1 0, not 2 -1 1 0."
                    )
                ),
            ) from error

    three_text = format_miller_indices(three, family=kind, style="plain", scope="specific")
    four_text = format_miller_indices(four, family=kind, style="plain", scope="specific")
    family_text = format_miller_indices(four, family=kind, style="plain", scope="family")

    antipodal = bool(request["antipodal"])
    if is_plane:
        members, mask = _plane_set(phase, [three]).symmetry_equivalent_indices()
    else:
        members, mask = _direction_set(phase, [three]).symmetry_equivalent_indices(
            antipodal=antipodal
        )
    valid = np.asarray(members[0])[np.asarray(mask[0], dtype=bool)]
    rows: list[dict[str, Any]] = []
    for position, row in enumerate(valid):
        triple = tuple(int(value) for value in row)
        quadruple = tuple(int(value) for value in widen(triple))
        rows.append(
            {
                "index": position + 1,
                "three_index": format_miller_indices(
                    triple, family=kind, style="plain", scope="specific"
                ),
                "four_index": format_miller_indices(
                    quadruple, family=kind, style="plain", scope="specific"
                ),
                "i1": quadruple[0],
                "i2": quadruple[1],
                "i3": quadruple[2],
                "i4": quadruple[3],
            }
        )

    # One geometric fact about the converted entity, so the answer is checkable
    # against something other than the arithmetic that produced it: a spacing
    # for a plane, and for a direction its inclination from the c axis, which is
    # 90 degrees for every basal direction whatever its indices.
    extras: dict[str, Any] = {}
    if is_plane:
        extras["d_spacing_angstrom"] = float(_plane_set(phase, [three]).d_spacings_angstrom()[0])
        fact = f"Every member has the interplanar spacing d = {extras['d_spacing_angstrom']:.5f} Å."
    else:
        c_axis = _direction_set(phase, [(0, 0, 1)])
        angle_deg = float(np.degrees(_direction_set(phase, [three]).angle_matrix_rad(c_axis)[0, 0]))
        extras["angle_to_c_axis_deg"] = angle_deg
        fact = f"It lies {angle_deg:.3f}° from the c axis [0001]" + (
            ", so it is a basal direction." if abs(angle_deg - 90.0) < 1e-6 else "."
        )

    symbols = "hkl → hkil" if is_plane else "uvw → UVTW"
    rule = (
        "the redundant index i = -(h + k) is inserted, because the reciprocal basis vectors "
        "already lie at 120° to one another"
        if is_plane
        else "the basal components are re-expressed on three coplanar axes through "
        "U = (2u - v)/3, V = (2v - u)/3, T = -(U + V), W = w, and the denominators cleared"
    )
    result = AppResult(
        title=f"{three_text} = {four_text} in {spec.name}",
        summary=(
            f"The {kind} {three_text} is written {four_text} in Miller-Bravais notation "
            f"({symbols}): {rule}. {fact} Its family {family_text} has {len(rows)} distinct "
            "members. In the four-index form those members are permutations of the first three "
            "indices, up to the sign convention in force; in the three-index form they are not, "
            "which is what the notation exists to fix."
        ),
        table=ResultTable(
            columns=(
                Column("index", "#", numeric=True),
                Column("three_index", "Three-index"),
                Column("four_index", "Four-index"),
                Column("i1", "h" if is_plane else "U", numeric=True),
                Column("i2", "k" if is_plane else "V", numeric=True),
                Column("i3", "i" if is_plane else "T", numeric=True),
                Column("i4", "l" if is_plane else "W", numeric=True),
            ),
            rows=tuple(rows),
            caption=f"Members of {family_text}, in both notations.",
        ),
        data={
            "kind": kind,
            "three_index": list(three),
            "four_index": list(four),
            "three_index_label": three_text,
            "four_index_label": four_text,
            "family_label": family_text,
            "multiplicity": len(rows),
            "members": [
                {"three_index": row["three_index"], "four_index": row["four_index"]} for row in rows
            ],
            **extras,
        },
        inputs={
            "phase": spec.to_json(),
            "kind": kind,
            "input_notation": given,
            "three_index": list(three),
            "four_index": list(four),
            "antipodal": antipodal,
        },
        notes=(
            "Planes and directions convert by different rules. A four-index plane is the same "
            "three numbers with i = -(h + k) inserted; a four-index direction is a genuine change "
            "of basis. The two therefore disagree in general — (100) is (10-10) while [100] is "
            "[2-1-10] — and where they do coincide, as (110) and [110] both do, that is a "
            "property of those particular indices rather than a rule.",
            "Directions are reduced to lowest terms, so [11-23] rather than [22-46].",
        ),
        citations=(_CITATION_ITA, _CITATION_FRANK, _CITATION_OTTE_CROCKER),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.zone_axis",
    title="Zone axis from two planes",
    summary="The common direction of two planes, and the reflections in that zone.",
    help_text=(
        "The zone axis of two planes is the direction lying in both, obtained as the cross "
        "product of their index triples and then reduced to the smallest integers.\n\n"
        "This is the operation behind indexing a diffraction pattern: two indexed spots fix "
        "the zone axis, and every other spot in the pattern must satisfy the zone law against "
        "it. The table lists the reflections that do, up to the index limit, which is the set "
        "of spots the pattern can contain."
    ),
    parameters=(
        phase_parameter(),
        IndicesParameter(
            name="first",
            label="First plane (hkl)",
            help_text="One of the two planes defining the zone.",
            default=(1, 1, 1),
        ),
        IndicesParameter(
            name="second",
            label="Second plane (hkl)",
            help_text="The other plane. It must not be parallel to the first.",
            default=(1, -1, 0),
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit for the zone list",
            help_text=(
                "Largest absolute value of h, k or l when listing reflections in the zone. "
                "Raise it to reach higher-order spots; the count grows roughly as its cube."
            ),
            default=3,
            minimum=0,
            maximum=8,
            advanced=True,
        ),
    ),
    returns="The zone axis under `data.zone_axis`; in-zone reflections as rows.",
    panel="calculator",
    citations=(_CITATION_ITA, _CITATION_CULLITY),
    tags=("zone axis", "cross product", "indexing", "Weiss zone law"),
)
def _zone_axis(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    first = np.asarray(request["first"], dtype=int)
    second = np.asarray(request["second"], dtype=int)
    axis = np.cross(first, second)
    if not np.any(axis):
        raise InvalidInputError(
            "The two planes are parallel, so they do not define a zone axis.",
            field="second",
            hint="Choose a second plane that is not a multiple of the first.",
        )
    axis = np.asarray(reduce_indices(axis[np.newaxis, :])[0], dtype=int)
    max_index = int(request["max_index"])
    rows: list[dict[str, Any]] = []
    if max_index > 0:
        grid = np.arange(-max_index, max_index + 1)
        candidates = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), axis=-1).reshape(-1, 3)
        candidates = candidates[np.any(candidates != 0, axis=1)]
        in_zone = candidates[candidates @ axis == 0]
        if in_zone.size:
            planes = _plane_set(phase, in_zone)
            spacings = np.asarray(planes.d_spacings_angstrom(), dtype=float)
            order = np.argsort(-spacings)
            for position in order:
                indices = tuple(int(value) for value in in_zone[position])
                rows.append(
                    {
                        "plane": plane_label(indices, spec=spec),
                        "h": indices[0],
                        "k": indices[1],
                        "l": indices[2],
                        "d_angstrom": float(spacings[position]),
                        "g_inv_angstrom": float(1.0 / spacings[position]),
                    }
                )
    axis_text = direction_label(tuple(int(value) for value in axis), spec=spec)
    result = AppResult(
        title=f"Zone axis {axis_text} in {spec.name}",
        summary=(
            f"{plane_label(tuple(int(v) for v in first), spec=spec)} and "
            f"{plane_label(tuple(int(v) for v in second), spec=spec)} intersect along "
            f"{axis_text}. Within |h|, |k|, |l| ≤ {max_index}, {len(rows)} reflections satisfy "
            "the zone law against it and could therefore appear in a pattern taken down this "
            "axis. Whether they actually appear also depends on the structure factor and on "
            "systematic absences, which this table does not apply."
        ),
        table=ResultTable(
            columns=(
                Column("plane", "Reflection"),
                Column("h", "h", numeric=True),
                Column("k", "k", numeric=True),
                Column("l", "l", numeric=True),
                Column("d_angstrom", "d", units="Å", numeric=True, digits=5),
                Column("g_inv_angstrom", "|g| = 1/d", units="Å⁻¹", numeric=True, digits=5),
            ),
            rows=tuple(rows),
            caption=f"Reflections in the {axis_text} zone of {spec.name}.",
        ),
        data={
            "zone_axis": [int(value) for value in axis],
            "zone_axis_label": axis_text,
            "reflection_count": len(rows),
        },
        inputs={
            "phase": spec.to_json(),
            "first": [int(value) for value in first],
            "second": [int(value) for value in second],
            "max_index": max_index,
        },
        notes=(
            "Systematic absences are not applied here; use the d-spacing table for a list "
            "filtered by the centring condition.",
        ),
        citations=(_CITATION_ITA, _CITATION_CULLITY),
    )
    return result.to_json()


@REGISTRY.operation(
    "calc.d_spacings",
    title="d-spacing table",
    summary="Reflections sorted by spacing, with multiplicity and centring absences.",
    help_text=(
        "Enumerates reflections up to an index limit, groups them into symmetry families, and "
        "reports the spacing, the reciprocal-lattice magnitude 1/d, and the multiplicity of "
        "each family.\n\n"
        "When the phase carries a space-group symbol, the centring letter is applied and "
        "systematically absent families are marked. Centring absences are the integral "
        "conditions only — glide and screw absences depend on the full space group and are not "
        "applied, so a family marked allowed here may still be extinct in a real pattern."
    ),
    parameters=(
        phase_parameter(),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest absolute value of h, k or l to enumerate.",
            default=3,
            minimum=1,
            maximum=8,
        ),
        NumberParameter(
            name="min_d_angstrom",
            label="Smallest spacing to report",
            help_text="Families with d below this are dropped. Set to 0 to keep everything.",
            units="Å",
            default=0.5,
            minimum=0.0,
        ),
        BooleanParameter(
            name="allowed_only",
            label="Hide systematically absent families",
            help_text=(
                "Applies the centring condition of the space group. Has no effect when the "
                "phase carries no space-group symbol."
            ),
            default=True,
        ),
    ),
    returns="One row per symmetry family, sorted by decreasing d.",
    panel="calculator",
    citations=(_CITATION_ITA, _CITATION_CULLITY),
    tags=("d-spacing", "reflection", "multiplicity", "absence", "powder", "XRD"),
)
def _d_spacings(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.diffraction import ReflectionCondition

    spec, phase = phase_from_request(request["phase"])
    max_index = int(request["max_index"])
    min_d = float(request["min_d_angstrom"])
    allowed_only = bool(request["allowed_only"])
    grid = np.arange(-max_index, max_index + 1)
    candidates = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), axis=-1).reshape(-1, 3)
    candidates = candidates[np.any(candidates != 0, axis=1)]
    planes = _plane_set(phase, candidates)
    spacings = np.asarray(planes.d_spacings_angstrom(), dtype=float)
    keep = spacings >= min_d
    candidates = candidates[keep]
    spacings = spacings[keep]
    condition = ReflectionCondition.from_phase(phase) if spec.space_group_symbol else None
    if not candidates.size:
        raise InvalidInputError(
            f"No reflection within |h|, |k|, |l| ≤ {max_index} has d ≥ {min_d:g} Å.",
            field="min_d_angstrom",
            hint="Lower the smallest spacing, or raise the index limit.",
        )
    # One vectorized pass for the whole grid: grouping family by family would
    # apply the point group thousands of times for the same answer.
    representatives, orbit_sizes = _plane_family_representatives(phase, candidates)
    seen: dict[tuple[int, ...], dict[str, Any]] = {}
    for indices, spacing, representative, orbit in zip(
        candidates, spacings, representatives, orbit_sizes, strict=True
    ):
        family_key = tuple(int(value) for value in representative)
        if family_key in seen:
            continue
        row = tuple(int(value) for value in indices)
        allowed = bool(condition.is_allowed(row)) if condition is not None else True
        seen[family_key] = {
            "family": family_label(family_key, spec=spec, family="plane"),
            "h": family_key[0],
            "k": family_key[1],
            "l": family_key[2],
            "d_angstrom": float(spacing),
            "g_inv_angstrom": float(1.0 / spacing),
            # Powder convention: (hkl) and (-h-k-l) diffract into different
            # points of the pattern, so both senses count, and the multiplicity
            # is twice the antipodal orbit size reported by the family operation.
            "multiplicity": 2 * int(orbit),
            "allowed": allowed,
        }
    rows = sorted(seen.values(), key=lambda entry: -float(entry["d_angstrom"]))
    absent = [row for row in rows if not row["allowed"]]
    if allowed_only and condition is not None:
        rows = [row for row in rows if row["allowed"]]
    notes: list[str] = []
    if condition is None:
        notes.append(
            "No space group was given, so no systematic absences were applied and every "
            "family is listed as allowed."
        )
    else:
        notes.append(
            f"Centring letter {condition.centering!r} applied; glide and screw absences are "
            "not included."
        )
    result = AppResult(
        title=f"d-spacings of {spec.name}",
        summary=(
            f"{len(rows)} reflection families with |h|, |k|, |l| ≤ {max_index} and "
            f"d ≥ {min_d:g} Å"
            + (
                f", after removing {len(absent)} families forbidden by the centring condition."
                if condition is not None and allowed_only
                else "."
            )
            + " Multiplicity is the powder convention — both senses of each family member — so "
            "it is what weights the corresponding powder peak: 8 for {111} and 6 for {200} in a "
            "cubic crystal."
        ),
        table=ResultTable(
            columns=(
                Column("family", "Family"),
                Column("h", "h", numeric=True),
                Column("k", "k", numeric=True),
                Column("l", "l", numeric=True),
                Column("d_angstrom", "d", units="Å", numeric=True, digits=5),
                Column("g_inv_angstrom", "1/d", units="Å⁻¹", numeric=True, digits=5),
                Column(
                    "multiplicity",
                    "Multiplicity",
                    numeric=True,
                    help_text=(
                        "Powder multiplicity: both senses of each family member counted, "
                        "so 8 for {111} in a cubic crystal."
                    ),
                ),
                Column("allowed", "Allowed"),
            ),
            rows=tuple(rows),
            caption=f"Reflection families of {spec.name}.",
        ),
        data={
            "family_count": len(rows),
            "forbidden_count": len(absent),
            "centering": str(condition.centering) if condition is not None else None,
        },
        inputs={
            "phase": spec.to_json(),
            "max_index": max_index,
            "min_d_angstrom": min_d,
            "allowed_only": allowed_only,
        },
        notes=notes,
        citations=(_CITATION_ITA, _CITATION_CULLITY),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Between two phases
# --------------------------------------------------------------------------


@REGISTRY.operation(
    "calc.interphase_angles",
    title="Angles between two phases",
    summary="Plane or direction angles across a stated orientation relationship.",
    help_text=(
        "Angles *between* two crystals are only defined once their relative orientation is "
        "stated — the indices alone say nothing, because they are components in two different "
        "bases. State the relationship as Bunge Euler angles or as an axis and angle, and this "
        "operation carries the second crystal's vectors into the first crystal's Cartesian "
        "frame before measuring.\n\n"
        "The natural use is checking a proposed orientation relationship: a Kurdjumov-Sachs "
        "relationship, for instance, should put a {111} of the parent within a fraction of a "
        "degree of a {110} of the product. Rows near 0° are the parallelisms the relationship "
        "asserts."
    ),
    parameters=(
        phase_parameter(
            help_text="The reference crystal. Its Cartesian frame is the one used.",
            builtin="austenite_fcc",
        ),
        _second_phase_parameter(builtin="fe_bcc"),
        ChoiceParameter(
            name="family",
            label="Compare",
            help_text="Whether to compare plane normals or lattice directions.",
            options=(
                ("plane", "Planes", "Angles between plane normals across the two crystals."),
                ("direction", "Directions", "Angles between lattice directions."),
            ),
            default="plane",
        ),
        IndicesListParameter(
            name="first_indices",
            label="Indices in phase 1",
            help_text="One row per plane or direction of the first crystal.",
            default=((1, 1, 1), (2, 0, 0)),
        ),
        IndicesListParameter(
            name="second_indices",
            label="Indices in phase 2",
            help_text="One row per plane or direction of the second crystal.",
            default=((0, 1, 1), (0, 0, 2)),
        ),
        *_rotation_parameters(
            convention="axis_angle",
            components=KS_VARIANT_AXIS,
            angle_deg=KS_VARIANT_ANGLE_DEG,
        ),
    ),
    returns="One row per cross-phase pair, sorted by increasing angle.",
    panel="calculator",
    citations=(
        _CITATION_ITA,
        "Bunge, Texture Analysis in Materials Science (1982), chapter 2.",
    ),
    tags=("orientation relationship", "interphase", "parallelism", "misorientation", "OR"),
)
def _interphase_angles(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    other_payload = request.get("other_phase")
    other_spec, other_phase = (
        phase_from_request(other_payload) if other_payload is not None else (spec, phase)
    )
    family = str(request["family"])
    first_rows = _require_rows(
        request["first_indices"], field="first_indices", what="plane or direction"
    )
    second_rows = _require_rows(
        request["second_indices"], field="second_indices", what="plane or direction"
    )
    rotation = _rotation_from_request(
        str(request["rotation_convention"]),
        (
            float(request["rotation_1"]),
            float(request["rotation_2"]),
            float(request["rotation_3"]),
        ),
        float(request["rotation_angle_deg"]),
    )
    if family == "plane":
        first_vectors = np.asarray(_plane_set(phase, first_rows).normals_cartesian(), dtype=float)
        second_vectors = np.asarray(
            _plane_set(other_phase, second_rows).normals_cartesian(), dtype=float
        )
    else:
        first_vectors = np.asarray(
            _direction_set(phase, first_rows).unit_vectors_cartesian(), dtype=float
        )
        second_vectors = np.asarray(
            _direction_set(other_phase, second_rows).unit_vectors_cartesian(), dtype=float
        )
    rotated = np.asarray(rotation.apply(second_vectors), dtype=float)
    cosines = np.clip(first_vectors @ rotated.T, -1.0, 1.0)
    angles_deg = np.degrees(np.arccos(np.abs(cosines)))
    labeller = plane_label if family == "plane" else direction_label
    entries: list[dict[str, Any]] = [
        {
            "first": labeller(first, spec=spec),
            "second": labeller(second, spec=other_spec),
            "angle_deg": float(angles_deg[i, j]),
        }
        for i, first in enumerate(first_rows)
        for j, second in enumerate(second_rows)
    ]
    rows = sorted(entries, key=lambda entry: float(entry["angle_deg"]))
    closest = rows[0]
    closest_angle = float(closest["angle_deg"])
    result = AppResult(
        title=f"{spec.name} against {other_spec.name}",
        summary=(
            f"{len(rows)} cross-phase pairs, measured after carrying {other_spec.name} into the "
            f"Cartesian crystal frame of {spec.name} through a rotation of "
            f"{rotation.angle_deg:.4f}°. Angles use antipodal equivalence and lie in [0°, 90°]. "
            f"The closest pair is {closest['first']} to {closest['second']} at "
            f"{closest_angle:.4f}°"
            + (
                ", which is a parallelism to within the tolerance normally quoted for an "
                "orientation relationship."
                if closest_angle < 1.0
                else "."
            )
        ),
        table=ResultTable(
            columns=(
                Column("first", f"In {spec.name}"),
                Column("second", f"In {other_spec.name}"),
                Column("angle_deg", "Angle", units="°", numeric=True, digits=_ANGLE_DIGITS),
            ),
            rows=tuple(rows),
            caption=f"Cross-phase angles, {spec.name} to {other_spec.name}.",
        ),
        data={
            "rotation_matrix": np.asarray(rotation.as_matrix(), dtype=float).tolist(),
            "rotation_angle_deg": float(rotation.angle_deg),
            "rotation_axis": [float(value) for value in np.asarray(rotation.axis, dtype=float)],
            "matrix": {
                "row_labels": [labeller(row, spec=spec) for row in first_rows],
                "column_labels": [labeller(row, spec=other_spec) for row in second_rows],
                "values_deg": [[float(value) for value in row] for row in angles_deg],
            },
        },
        inputs={
            "phase": spec.to_json(),
            "other_phase": other_spec.to_json(),
            "family": family,
            "first_indices": [list(row) for row in first_rows],
            "second_indices": [list(row) for row in second_rows],
            "rotation_convention": request["rotation_convention"],
            "rotation_values": [
                float(request["rotation_1"]),
                float(request["rotation_2"]),
                float(request["rotation_3"]),
            ],
            "rotation_angle_deg": float(request["rotation_angle_deg"]),
        },
        citations=(_CITATION_ITA,),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Named orientation relationships
# --------------------------------------------------------------------------

#: The named relationships the calculator offers, as choice options.
#:
#: Names and defining parallelisms follow the literature they are attributed to;
#: the constructors live in :mod:`pytex.core.transformation`, so the application
#: never re-derives a relationship it can look up.
#: The identifier of the user-defined relationship.
#:
#: Not a constructor on :class:`OrientationRelationship` but a *statement* the
#: user makes -- a parent plane parallel to a child plane, a parent direction
#: parallel to a child direction -- which is what every named relationship in
#: the table below also is. It is resolved by :func:`resolve_relationship`
#: through the same ``from_parallel_plane_direction`` constructor the named ones
#: use, so a custom relationship is not a second class of object: variants,
#: packets, intervariant boundaries and composite diffraction all work on it
#: unchanged.
CUSTOM_RELATIONSHIP = "custom"

_RELATIONSHIPS: tuple[tuple[str, str, str], ...] = (
    (
        "kurdjumov_sachs",
        "Kurdjumov-Sachs (fcc to bcc)",
        "{111}γ ∥ {110}α with ⟨110⟩γ ∥ ⟨111⟩α; 24 variants.",
    ),
    (
        "nishiyama_wassermann",
        "Nishiyama-Wassermann (fcc to bcc)",
        "{111}γ ∥ {110}α with ⟨112⟩γ ∥ ⟨110⟩α; 12 variants.",
    ),
    ("bain", "Bain (fcc to bcc)", "The pure strain path; 3 variants, no shear."),
    (
        "greninger_troiano",
        "Greninger-Troiano (fcc to bcc)",
        "An irrational relationship lying between Kurdjumov-Sachs and Nishiyama-Wassermann.",
    ),
    ("pitsch", "Pitsch (fcc to bcc)", "{100}γ ∥ {110}α with ⟨110⟩γ ∥ ⟨111⟩α."),
    (
        "burgers",
        "Burgers (bcc to hcp)",
        "{110}β ∥ {0001}α with ⟨111⟩β ∥ ⟨11-20⟩α; the titanium and zirconium relationship.",
    ),
    (
        "cube_on_cube",
        "Cube-on-cube (cubic to cubic)",
        "(001) ∥ (001) with [100] ∥ [100]; parallel axes, the identity rotation; 1 variant.",
    ),
    (
        "fcc_twin",
        "Coherent twin (fcc, Σ3)",
        "(111) ∥ (111) with [1-10] ∥ [-110]; 60° about ⟨111⟩ after reduction; 4 variants.",
    ),
    (
        CUSTOM_RELATIONSHIP,
        "Custom — your own parallelisms",
        "State a plane pair and a direction pair yourself; everything else is derived from them.",
    ),
)

#: How each relationship is written in prose and in a title.
#:
#: The identifier is a slug (``kurdjumov_sachs``) and the choice label carries a
#: parenthetical the picker needs but a sentence does not
#: (``Kurdjumov-Sachs (fcc to bcc)``). Titles and prose take the name alone, from
#: here, so that no surface has to reconstruct it — and so that a heading never
#: reads ``kurdjumov-sachs``, which is a variable name, not a pair of surnames.
RELATIONSHIP_NAMES: dict[str, str] = {
    identifier: label.split(" (")[0] for identifier, label, _ in _RELATIONSHIPS
} | {
    # A sentence says "the custom relationship", not "the Custom - your own
    # parallelisms relationship". The picker needs the longer label to explain
    # itself; prose needs the noun.
    CUSTOM_RELATIONSHIP: "custom",
    # Catalogued for comparison but not offered as a construction, so it has no
    # entry among the options above. It still reaches the screen: every
    # characterization of a bcc-to-hcp pair ranks it as the runner-up to
    # Burgers, and without a display name it printed as "shoji-nishiyama"
    # beside properly cased names.
    "shoji_nishiyama": "Shoji-Nishiyama",
}


def relationship_name(identifier: str) -> str:
    """Return the human-readable name of a relationship identifier."""

    return RELATIONSHIP_NAMES.get(identifier, identifier.replace("_", "-"))


#: Constructor name on :class:`OrientationRelationship` for each option.
_RELATIONSHIP_CONSTRUCTORS = {
    "kurdjumov_sachs": "from_kurdjumov_sachs_correspondence",
    "nishiyama_wassermann": "from_nishiyama_wassermann_correspondence",
    "bain": "from_bain_correspondence",
    "greninger_troiano": "from_greninger_troiano_correspondence",
    "pitsch": "from_pitsch_correspondence",
    "burgers": "from_burgers_correspondence",
    "cube_on_cube": "from_cube_on_cube_correspondence",
    "fcc_twin": "from_fcc_twin_correspondence",
}


def custom_relationship_parameters(
    *, group: str = "Custom relationship", collapsed: bool = True
) -> tuple[Parameter, ...]:
    """The four index rows that define a relationship the user states themselves.

    Purpose
    -------
    A published orientation relationship is a pair of parallelisms, and the
    named entries in :data:`_RELATIONSHIPS` are nothing more than those pairs
    written down once. These parameters let a user write down a pair that is not
    on the list -- one from a paper, one they have just fitted, or one they want
    to test a hypothesis with -- and get the same treatment: the misorientation,
    every crystallographically distinct variant, the packet grouping, the
    intervariant boundaries, and the composite diffraction pattern.

    Declared once and shared by every operation that offers the relationship
    picker, so a custom relationship means the same thing in the variant table
    as it does in the composite pattern.

    When to use
    -----------
    On an operation whose ``relationship`` picker offers *Custom*, where they
    apply only to that choice: they are ignored otherwise, and grouped behind
    their own disclosure so they cost nothing on screen while the relationship
    is a named one. On an operation that is *about* a user-stated relationship
    -- ``variants.custom_relationship`` -- they are the subject rather than a
    fallback, so it passes its own heading and ``collapsed=False`` and they open
    with the panel.

    Parameters
    ----------
    group : str
        Control-panel section the four rows are gathered under.
    collapsed : bool
        Whether that section starts closed.

    Returns
    -------
    tuple of Parameter
        Parent plane, child plane, parent direction, child direction.
    """

    return (
        IndicesParameter(
            name="custom_parent_plane",
            label="Parent plane (hkl)",
            help_text=(
                "The parent plane held parallel to the child plane. With the direction pair "
                "below, this is the whole statement of the relationship: Kurdjumov-Sachs is "
                "(111) here, (011) as the child plane, [-101] and [-1-11] as the directions."
            ),
            default=(1, 1, 1),
            group=group,
            group_collapsed=collapsed,
        ),
        IndicesParameter(
            name="custom_child_plane",
            label="Child plane (hkl)",
            help_text="The child plane brought parallel to the parent plane above.",
            default=(0, 1, 1),
            group=group,
            group_collapsed=collapsed,
        ),
        IndicesParameter(
            name="custom_parent_direction",
            label="Parent direction [uvw]",
            help_text=(
                "The parent direction held parallel to the child direction, fixing the rotation "
                "that remains about the common plane normal. It need not lie exactly in the "
                "plane: the normal component is removed, so a literature statement that is only "
                "approximately consistent still yields a proper rotation."
            ),
            default=(-1, 0, 1),
            group=group,
            group_collapsed=collapsed,
        ),
        IndicesParameter(
            name="custom_child_direction",
            label="Child direction [uvw]",
            help_text="The child direction brought parallel to the parent direction above.",
            default=(-1, -1, 1),
            group=group,
            group_collapsed=collapsed,
        ),
    )


def custom_relationship_request(request: dict[str, Any]) -> dict[str, Any]:
    """The ``custom_*`` fields of a request, ready to forward to a sub-request.

    Purpose
    -------
    Some operations compute by calling another operation's handler with a
    hand-built request -- the variant figure draws from the variant pole figure
    rather than reimplementing it -- and a hand-built request lists its keys.
    A custom relationship adds four more, and a forwarder that does not know
    about them turns "Custom" into a missing field at the far end.

    Returning them as a mapping to splat in means a forwarder cannot list three
    of the four, and a fifth added later reaches every forwarder at once.
    """

    return {
        parameter.name: request[parameter.name]
        for parameter in custom_relationship_parameters()
        if parameter.name in request
    }


def resolve_relationship(
    request: dict[str, Any],
    parent_phase: Any,
    child_phase: Any,
    *,
    name: str | None = None,
    custom_name: str | None = None,
) -> Any:
    """Build the relationship an operation was asked for, named or custom.

    Purpose
    -------
    The one place that turns the ``relationship`` choice into an
    :class:`~pytex.core.transformation.OrientationRelationship`. Three services
    offer the same picker -- the calculator, the variant tools, and composite
    diffraction -- and each used to carry its own copy of the lookup, which is
    three places for an option to go missing from.

    Parameters
    ----------
    request : dict
        The validated request. Read for ``relationship`` unless ``name`` says
        otherwise, and for the ``custom_*`` index rows when the choice is
        ``"custom"``.
    parent_phase, child_phase : Phase
        The two phases the relationship is built between.
    name : str, optional
        Use this identifier instead of the one in the request.
    custom_name : str, optional
        The name a user-stated relationship is built under, and therefore the
        name it carries into ``describe()``, reports and figures. Ignored for a
        named relationship, whose name is its own. Defaults to ``"custom"``.

    Returns
    -------
    OrientationRelationship

    Raises
    ------
    InvalidInputError
        When a named relationship does not apply to these phases, or when a
        custom statement does not define a rotation. Both name ``relationship``
        so the message lands on the control the user chose from.
    """

    from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex
    from pytex.core.transformation import OrientationRelationship

    identifier = str(name if name is not None else request["relationship"])
    if identifier != CUSTOM_RELATIONSHIP:
        constructor = getattr(OrientationRelationship, _RELATIONSHIP_CONSTRUCTORS[identifier])
        try:
            return constructor(parent_phase=parent_phase, child_phase=child_phase)
        except (ValueError, TypeError) as error:
            raise InvalidInputError(
                f"The {relationship_name(identifier)} relationship does not apply to these "
                f"phases: {error}",
                field="relationship",
                hint=(
                    "The fcc-to-bcc relationships, cube-on-cube and the coherent twin need a "
                    "cubic parent and a cubic child; Burgers needs a cubic parent and a "
                    "hexagonal child. Choose Custom to state a relationship these two phases "
                    "can actually carry."
                ),
            ) from error

    def _row(field: str) -> list[float]:
        row = request.get(field)
        if row is None:
            raise InvalidInputError(
                "A custom relationship needs all four index rows.",
                field=field,
                hint=(
                    "Choose a named relationship, or fill in the parent and child plane and "
                    "the parent and child direction."
                ),
            )
        return [float(value) for value in row]

    def _plane(field: str, phase: Any) -> Any:
        indices = np.asarray([int(value) for value in _row(field)], dtype=np.int64)
        return CrystalPlane(MillerIndex(indices=indices, phase=phase), phase=phase)

    def _direction(field: str, phase: Any) -> Any:
        return CrystalDirection(np.asarray(_row(field), dtype=np.float64), phase=phase)

    try:
        return OrientationRelationship.from_parallel_plane_direction(
            name=custom_name or CUSTOM_RELATIONSHIP,
            parent_plane=_plane("custom_parent_plane", parent_phase),
            child_plane=_plane("custom_child_plane", child_phase),
            parent_direction=_direction("custom_parent_direction", parent_phase),
            child_direction=_direction("custom_child_direction", child_phase),
        )
    except (ValueError, TypeError, np.linalg.LinAlgError) as error:
        raise InvalidInputError(
            f"Those parallelisms do not define a rotation: {error}",
            field="relationship",
            hint=(
                "The direction must not be parallel to its own plane normal — with nothing "
                "left lying in the plane there is no rotation about the normal to fix. "
                "Kurdjumov-Sachs is the worked shape: parent (111) with [-101], child (011) "
                "with [-1-11]."
            ),
        ) from error


@REGISTRY.operation(
    "calc.orientation_relationship",
    title="Named orientation relationship",
    summary="Parallelisms, misorientation, and variants of a classical OR.",
    help_text=(
        "Builds one of the classical orientation relationships between the two phases and "
        "reports what it actually asserts: the defining plane and direction parallelisms, the "
        "misorientation angle and axis, and every crystallographically distinct variant.\n\n"
        "The variant table is the reason to run this rather than to quote a number. One parent "
        "grain transforming under Kurdjumov-Sachs produces 24 child orientations, not one, and "
        "that multiplicity is what makes packets and blocks visible in an EBSD map.\n\n"
        "Bain is included for contrast: 3 variants and no shear, so it is the reference against "
        "which the shear-carrying relationships are understood."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text=(
                "The phase that transforms — austenite for the fcc-to-bcc relationships, beta "
                "for Burgers."
            ),
            builtin="austenite_fcc",
        ),
        ObjectParameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase — ferrite or martensite, or alpha for Burgers.",
            editor="phase",
            default={"builtin": "fe_bcc"},
        ),
        ChoiceParameter(
            name="relationship",
            label="Relationship",
            help_text=(
                "Which relationship to build. **Custom** takes the parallelisms from the boxes "
                "below instead of from this list, and everything else — the misorientation, the "
                "variants, the parallelism table — is derived from them exactly as it is for a "
                "named relationship."
            ),
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        *custom_relationship_parameters(),
        BooleanParameter(
            name="reduce_by_child_symmetry",
            label="Reduce variants by child symmetry",
            help_text=(
                "On, the list holds one entry per crystallographically distinct child "
                "orientation, which is the count quoted in the literature. Off lists every "
                "symmetry image, which is longer and rarely what is wanted."
            ),
            default=True,
            advanced=True,
        ),
    ),
    returns="One row per variant, with the parallelisms and the OR prose under `data`.",
    panel="calculator",
    citations=(
        "Kurdjumov & Sachs, Z. Phys. 64 (1930) 325.",
        "Nishiyama, Sci. Rep. Tohoku Univ. 23 (1934) 637.",
        "Burgers, Physica 1 (1934) 561.",
        "Morito et al., Acta Mater. 51 (2003) 1789 (variant numbering).",
    ),
    tags=("orientation relationship", "variant", "martensite", "Bain", "Burgers", "OR"),
)
def _orientation_relationship(request: dict[str, Any]) -> dict[str, Any]:
    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    name = str(request["relationship"])
    relationship = resolve_relationship(request, parent_phase, child_phase, name=name)
    variants = relationship.generate_variants(
        reduce_by_child_symmetry=bool(request["reduce_by_child_symmetry"])
    )
    representative = relationship.misorientation()
    rows: list[dict[str, Any]] = []
    for variant in variants:
        rotation = variant.parent_to_child_rotation
        variant_axis = np.asarray(rotation.axis, dtype=float)
        rows.append(
            {
                "variant": int(variant.variant_index),
                "angle_deg": float(rotation.angle_deg),
                "axis_x": float(variant_axis[0]),
                "axis_y": float(variant_axis[1]),
                "axis_z": float(variant_axis[2]),
            }
        )
    plane_pairs = [
        {
            "parent": plane_label(
                tuple(int(value) for value in parent.miller.as_array()), spec=parent_spec
            ),
            "child": plane_label(
                tuple(int(value) for value in child.miller.as_array()), spec=child_spec
            ),
        }
        for parent, child in relationship.parallel_planes
    ]
    direction_pairs = [
        {
            "parent": direction_label(
                tuple(round(value) for value in parent.coordinates), spec=parent_spec
            ),
            "child": direction_label(
                tuple(round(value) for value in child.coordinates), spec=child_spec
            ),
        }
        for parent, child in relationship.parallel_directions
    ]
    axis = np.asarray(representative.rotation.axis, dtype=float)
    label = relationship_name(name)
    result = AppResult(
        title=f"{label}: {parent_spec.name} to {child_spec.name}",
        summary=(
            f"The {label} relationship carries {parent_spec.name} into {child_spec.name} with "
            + (
                f"{plane_pairs[0]['parent']} ∥ {plane_pairs[0]['child']}"
                if plane_pairs
                else "no defining plane parallelism"
            )
            + (
                f" and {direction_pairs[0]['parent']} ∥ {direction_pairs[0]['child']}"
                if direction_pairs
                else ""
            )
            + f". The representative misorientation is {representative.angle_deg:.2f}° about "
            f"[{axis[0]:.3f} {axis[1]:.3f} {axis[2]:.3f}], and one parent grain produces "
            f"{len(rows)} crystallographically distinct child orientations."
        ),
        table=ResultTable(
            columns=(
                Column("variant", "Variant", numeric=True),
                Column("angle_deg", "Angle", units="°", numeric=True, digits=4),
                Column("axis_x", "Axis x", numeric=True, digits=5),
                Column("axis_y", "Axis y", numeric=True, digits=5),
                Column("axis_z", "Axis z", numeric=True, digits=5),
            ),
            rows=tuple(rows),
            caption=f"Variants of the {label} relationship.",
        ),
        data={
            "relationship": name,
            "variant_count": len(rows),
            "misorientation_angle_deg": float(representative.angle_deg),
            "misorientation_axis": [float(value) for value in axis],
            "parallel_planes": plane_pairs,
            "parallel_directions": direction_pairs,
            "description": relationship.describe(),
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": name,
            "reduce_by_child_symmetry": bool(request["reduce_by_child_symmetry"]),
        },
        citations=("Morito et al., Acta Mater. 51 (2003) 1789.",),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Canonical examples
# --------------------------------------------------------------------------

REGISTRY.add_examples(
    (
        ExampleScenario(
            id="calc.example.fcc_powder",
            title="Why an fcc powder pattern starts at {111}",
            panel="calculator",
            summary="The reflection families of austenite, sorted by spacing.",
            teaches=(
                "The first four allowed families are 111, 200, 220 and 311 — the mixed-parity "
                "families such as 100 and 110 are absent because the F centring cancels them "
                "exactly. The multiplicity column, 8/6/12/24, is what makes 311 a strong peak "
                "despite its small spacing."
            ),
            operation="calc.d_spacings",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "max_index": 4,
                "min_d_angstrom": 0.8,
            },
        ),
        ExampleScenario(
            id="calc.example.bcc_vs_fcc",
            title="Ferrite: the other centring condition",
            panel="calculator",
            summary="The reflection families of bcc iron, for comparison with austenite.",
            teaches=(
                "Body centring keeps h+k+l even, so the sequence is 110, 200, 211, 220 — a "
                "different fingerprint from the fcc one. Comparing the two tables is how a bcc "
                "pattern is told from an fcc pattern without indexing a single spot."
            ),
            operation="calc.d_spacings",
            request={"phase": {"builtin": "fe_bcc"}, "max_index": 4, "min_d_angstrom": 0.8},
        ),
        ExampleScenario(
            id="calc.example.hcp_angles",
            title="Basal, prism and pyramidal planes in zirconium",
            panel="calculator",
            summary="Angles among the three plane families that carry hcp slip.",
            teaches=(
                "The basal plane is exactly perpendicular to every prism plane whatever c/a is, "
                "but the pyramidal angle is not fixed by symmetry — it depends on c/a, which is "
                "why the same table for magnesium or titanium gives a different number."
            ),
            operation="calc.plane_angles",
            request={
                "phase": {"builtin": "zr_hcp"},
                "planes": [[0, 0, 1], [1, 0, 0], [1, 0, 1]],
            },
        ),
        ExampleScenario(
            id="calc.example.basal_slip_indices",
            title="Why the hcp slip direction is written [2-1-10] and not [100]",
            panel="calculator",
            summary="The close-packed basal direction of titanium, in both notations.",
            teaches=(
                "The three close-packed basal directions are crystallographically identical, but "
                "in three-index form they read [100], [010] and [-1-10], which share no visible "
                "relationship. Converted, they are [2-1-10], [-12-10] and [-1-120]: permutations "
                "of one another, and obviously one family. That is the whole reason the "
                "four-index notation exists, and why the hcp slip system is quoted as "
                "(0001)<11-20>. Reversed senses are listed here, so the six members are the "
                "permutations of the first three indices and their reverses. Note also that a "
                "direction converts by a change of basis, not by inserting an index: the "
                "plane (100) becomes (1 0 -1 0) while the direction [100] becomes [2 -1 -1 0]."
            ),
            operation="calc.hexagonal_indices",
            request={
                "phase": {"builtin": "ti_hcp"},
                "kind": "direction",
                "input_notation": "three",
                "three_index": [1, 0, 0],
                "four_index": [1, 1, -2, 0],
                "antipodal": False,
            },
        ),
        ExampleScenario(
            id="calc.example.nacl_family",
            title="Why {100} has three members, not forty-eight",
            panel="calculator",
            summary="The symmetry family of the cube plane in halite.",
            teaches=(
                "The point group m-3m has 48 operations but the family has 3 members, because "
                "16 of those operations map (100) onto itself or onto its opposite. That ratio "
                "— not the group order — is what sets the weight of the corresponding peak."
            ),
            operation="calc.symmetry_family",
            request={"phase": {"builtin": "nacl"}, "family": "plane", "indices": [1, 0, 0]},
        ),
        ExampleScenario(
            id="calc.example.ks_variants",
            title="Why one austenite grain gives 24 martensite orientations",
            panel="calculator",
            summary="The Kurdjumov-Sachs relationship between austenite and ferrite.",
            teaches=(
                "K-S asserts (111)γ ∥ (011)α and ⟨110⟩γ ∥ ⟨111⟩α, and the misorientation that "
                "satisfies both is 42.85° about ⟨0.968 0.178 0.178⟩. The 24 rows are the reason "
                "a prior-austenite grain shows a patchwork rather than one orientation in an "
                "EBSD map."
            ),
            operation="calc.orientation_relationship",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
            },
        ),
        ExampleScenario(
            id="calc.example.burgers",
            title="The Burgers path into hexagonal zirconium",
            panel="calculator",
            summary="The bcc-to-hcp relationship that governs zirconium and titanium.",
            teaches=(
                "{110}β ∥ {0001}α with ⟨111⟩β ∥ ⟨11-20⟩α gives 12 variants rather than 24, "
                "because the hexagonal child has lower symmetry to absorb. This is the same "
                "machinery as the steel example, applied across a change of crystal system."
            ),
            operation="calc.orientation_relationship",
            request={
                "phase": {"builtin": "zr_bcc_beta"},
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
            },
        ),
    )
)
