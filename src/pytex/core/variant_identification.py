"""Which variant is this child grain? One measured pair, every step kept.

Purpose
-------
A parent grain transforming under an orientation relationship produces a
finite set of child orientations, the variants. Given the measured orientation
of one parent grain and of one child grain next to it, the question an EBSD
analyst asks is *which* of those variants the child is. This module answers it
and keeps every intermediate quantity of the answer, so that a report can show
the calculation rather than only its outcome.

The calculation, in the canonical crystal-to-specimen convention
(``docs/standards/terminology_and_symbol_registry.md``):

1. The measured orientations are the matrices :math:`\\mathbf{P}` and
   :math:`\\mathbf{C}`, each taking crystal-frame components to specimen-frame
   components.
2. The measured parent-to-child rotation is
   :math:`\\mathbf{V}_{\\mathrm{m}} = \\mathbf{C}^{\\mathsf{T}}\\mathbf{P}`.
3. The variants of the relationship :math:`\\mathbf{R}` are
   :math:`\\mathbf{V}_k = \\mathbf{R}\\,\\mathbf{S}_{p,k}^{\\mathsf{T}}` for parent
   point-group operators :math:`\\mathbf{S}_{p,k} \\in G_p`, kept when distinct
   modulo the child point group :math:`G_c`
   (:meth:`~pytex.core.transformation.OrientationRelationship.generate_variants`).
4. Variant :math:`k` predicts the child orientation
   :math:`\\mathbf{C}_k = \\mathbf{P}\\,\\mathbf{V}_k^{\\mathsf{T}}`. The measured
   child is only defined up to its own symmetry, :math:`\\mathbf{C}\\mathbf{S}_c`,
   so the misfit of variant :math:`k` is the smallest rotation angle of
   :math:`\\Delta\\mathbf{g}_k = \\mathbf{C}_k^{\\mathsf{T}}\\mathbf{C}\\,\\mathbf{S}_c
   = \\mathbf{V}_k\\mathbf{V}_{\\mathrm{m}}^{\\mathsf{T}}\\mathbf{S}_c` over
   :math:`\\mathbf{S}_c \\in G_c`.
5. The variant with the smallest misfit is the identification; the gap to the
   runner-up is what makes it unambiguous.

As a physical cross-check, the identified variant's own parallelism — the
parent plane and direction it carries onto the child's — is placed in the
specimen frame through both *measured* orientations, and the angle between the
two members is reported. For a correct identification it equals the misfit in
size, and it is the quantity a metallurgist recognizes.

When to use it
--------------
For a handful of hand-picked grain pairs, when the reasoning matters as much as
the label: teaching, checking a reconstruction, or documenting a variant in a
paper. For whole maps use :func:`pytex.core.transformation.or_deviation` or
:func:`pytex.core.parent_reconstruction.select_variants`, which assign the same
variants without retaining the intermediates; the tests hold the three to the
same answer.

Variant numbering
-----------------
Variant indices are those of
:meth:`~pytex.core.transformation.OrientationRelationship.generate_variants`,
which is PyTex's enumeration order, not a published table (only V1 of
Kurdjumov-Sachs is pinned to Morito et al.). The report therefore also states
the variant by the parallelism it realizes, which is convention-free.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._angles import angle_between_unit_vectors_rad, rotation_angle_from_matrix_rad
from pytex.core.lattice import CrystalDirection, CrystalPlane, phases_semantically_match
from pytex.core.orientation import Orientation, Rotation
from pytex.core.transformation import (
    OrientationRelationship,
    TransformationVariant,
    _symmetry_image_direction,
    _symmetry_image_plane,
)

__all__ = [
    "DEFAULT_ASSIGNMENT_TOLERANCE_DEG",
    "ParallelismCheck",
    "VariantIdentificationReport",
    "VariantVerdict",
    "identify_variant",
    "operator_axis_angle",
]

#: Misfit below which a variant is accepted as the child's, in degrees.
#:
#: Five degrees: well above the angular resolution of a good EBSD map (about
#: half a degree) and the scatter of a real transformation about its nominal
#: relationship (Greninger-Troiano sits 2.4 degrees from Kurdjumov-Sachs), and
#: well below the smallest Kurdjumov-Sachs intervariant disorientation, 10.53
#: degrees (Morito et al. 2003), so two variants cannot both be inside it.
DEFAULT_ASSIGNMENT_TOLERANCE_DEG = 5.0

VariantVerdict = Literal["assigned", "ambiguous", "no_fit"]


@dataclass(frozen=True, slots=True)
class ParallelismCheck:
    """One defining parallelism of the identified variant, tested on the measurement.

    Attributes
    ----------
    kind : {"plane", "direction"}
    parent : CrystalPlane or CrystalDirection
        The parent member, in the measured parent's own indices.
    child : CrystalPlane or CrystalDirection
        The child member, in the measured child's own indices: the nominal
        child member carried by the matching child operator.
    parent_specimen, child_specimen : np.ndarray
        Unit vectors of the two members in the specimen frame, computed with
        the *measured* orientations.
    angle_deg : float
        Angle between them. Zero for an exact relationship and error-free data.
    """

    kind: Literal["plane", "direction"]
    parent: CrystalPlane | CrystalDirection
    child: CrystalPlane | CrystalDirection
    parent_specimen: np.ndarray
    child_specimen: np.ndarray
    angle_deg: float


@dataclass(frozen=True, slots=True)
class VariantIdentificationReport:
    """Every quantity behind the identification of one child grain's variant.

    Attributes
    ----------
    relationship : OrientationRelationship
    parent_matrix, child_matrix : np.ndarray
        :math:`\\mathbf{P}` and :math:`\\mathbf{C}`, crystal to specimen.
    measured_rotation : np.ndarray
        :math:`\\mathbf{V}_{\\mathrm{m}} = \\mathbf{C}^{\\mathsf{T}}\\mathbf{P}`.
    variants : tuple of TransformationVariant
        The variants compared, in generation order.
    misfits_deg : np.ndarray
        :math:`\\theta_k` for each variant, in degrees.
    child_operator_indices : np.ndarray
        Index into the child point group of the operator achieving each
        :math:`\\theta_k`.
    best_position : int
        Position (zero-based) in ``variants`` of the smallest misfit.
    runner_up_position : int or None
        Position of the second smallest misfit; ``None`` for a one-variant
        relationship.
    tolerance_deg : float
        The assignment tolerance :math:`\\tau`.
    parallelism_checks : tuple of ParallelismCheck
        The identified variant's defining parallelisms, measured.
    """

    relationship: OrientationRelationship
    parent_matrix: np.ndarray
    child_matrix: np.ndarray
    measured_rotation: np.ndarray
    variants: tuple[TransformationVariant, ...]
    misfits_deg: np.ndarray
    child_operator_indices: np.ndarray
    best_position: int
    runner_up_position: int | None
    tolerance_deg: float
    parallelism_checks: tuple[ParallelismCheck, ...]

    def __post_init__(self) -> None:
        count = len(self.variants)
        if count == 0:
            raise ValueError("A variant identification needs at least one variant.")
        misfits = np.asarray(self.misfits_deg, dtype=np.float64).reshape(-1)
        operators = np.asarray(self.child_operator_indices, dtype=np.int64).reshape(-1)
        if misfits.shape != (count,) or operators.shape != (count,):
            raise ValueError("misfits_deg and child_operator_indices need one entry per variant.")
        if not 0 <= self.best_position < count:
            raise ValueError("best_position must index a variant.")
        if not np.isfinite(self.tolerance_deg) or self.tolerance_deg <= 0.0:
            raise ValueError("tolerance_deg must be a positive finite angle.")
        for name in ("parent_matrix", "child_matrix", "measured_rotation"):
            array = np.array(getattr(self, name), dtype=np.float64)
            if array.shape != (3, 3):
                raise ValueError(f"{name} must be a 3x3 matrix.")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        misfits = misfits.copy()
        operators = operators.copy()
        misfits.setflags(write=False)
        operators.setflags(write=False)
        object.__setattr__(self, "misfits_deg", misfits)
        object.__setattr__(self, "child_operator_indices", operators)

    # -- the answer ---------------------------------------------------------

    @property
    def variant(self) -> TransformationVariant:
        """The identified variant."""

        return self.variants[self.best_position]

    @property
    def variant_index(self) -> int:
        """One-based index of the identified variant."""

        return int(self.variant.variant_index)

    @property
    def misfit_deg(self) -> float:
        """:math:`\\theta^{*}`, the misfit of the identified variant."""

        return float(self.misfits_deg[self.best_position])

    @property
    def runner_up(self) -> TransformationVariant | None:
        """The second-best variant, or ``None`` for a one-variant relationship."""

        if self.runner_up_position is None:
            return None
        return self.variants[self.runner_up_position]

    @property
    def runner_up_misfit_deg(self) -> float | None:
        """Misfit of the second-best variant."""

        if self.runner_up_position is None:
            return None
        return float(self.misfits_deg[self.runner_up_position])

    @property
    def margin_deg(self) -> float | None:
        """How far the runner-up trails the identified variant, in degrees."""

        runner_up = self.runner_up_misfit_deg
        return None if runner_up is None else runner_up - self.misfit_deg

    @property
    def verdict(self) -> VariantVerdict:
        """``assigned``, ``ambiguous`` or ``no_fit`` against the tolerance.

        ``no_fit`` when even the best variant lies outside the tolerance: the
        pair does not obey this relationship, or the two grains are not parent
        and child. ``ambiguous`` when the runner-up is inside it too.
        """

        if self.misfit_deg > self.tolerance_deg:
            return "no_fit"
        runner_up = self.runner_up_misfit_deg
        if runner_up is not None and runner_up <= self.tolerance_deg:
            return "ambiguous"
        return "assigned"

    # -- the operators that produced it ---------------------------------------

    @property
    def parent_operators(self) -> np.ndarray:
        """The parent point group :math:`G_p` as an ``(n, 3, 3)`` array."""

        return _operators(self.relationship.parent_phase.symmetry)

    @property
    def child_operators(self) -> np.ndarray:
        """The child point group :math:`G_c` as an ``(n, 3, 3)`` array."""

        return _operators(self.relationship.child_phase.symmetry)

    @property
    def parent_operator(self) -> np.ndarray:
        """:math:`\\mathbf{S}_{p,k^{*}}`, which generated the identified variant."""

        return np.asarray(self.variant.parent_symmetry_operator, dtype=np.float64)

    @property
    def child_operator_index(self) -> int:
        """Index in :math:`G_c` of the operator matching the identified variant."""

        return int(self.child_operator_indices[self.best_position])

    @property
    def child_operator(self) -> np.ndarray:
        """:math:`\\mathbf{S}_{c}^{*}`, the child operator achieving :math:`\\theta^{*}`."""

        return np.asarray(self.child_operators[self.child_operator_index], dtype=np.float64)

    @property
    def variant_matrix(self) -> np.ndarray:
        """:math:`\\mathbf{V}_{k^{*}}`."""

        return np.asarray(self.variant.parent_to_child_rotation.as_matrix(), dtype=np.float64)

    @property
    def predicted_child_matrix(self) -> np.ndarray:
        """:math:`\\mathbf{C}_{k^{*}} = \\mathbf{P}\\,\\mathbf{V}_{k^{*}}^{\\mathsf{T}}`."""

        predicted: np.ndarray = self.parent_matrix @ self.variant_matrix.T
        return predicted

    @property
    def residual_matrix(self) -> np.ndarray:
        """The residual rotation of the identified variant.

        :math:`\\Delta\\mathbf{g} =
        \\mathbf{V}_{k^{*}}\\mathbf{V}_{\\mathrm{m}}^{\\mathsf{T}}\\mathbf{S}_{c}^{*}`.
        """

        residual: np.ndarray = self.variant_matrix @ self.measured_rotation.T @ self.child_operator
        return residual

    def describe(self) -> str:
        """The identification as prose, with its conventions and its margin."""

        relationship = self.relationship
        parent = relationship.parent_phase.name
        child = relationship.child_phase.name
        statement = "; ".join(
            f"{check.kind} {_indices_text(check.parent)} parallel to "
            f"{_indices_text(check.child)} within {check.angle_deg:.2f} deg"
            for check in self.parallelism_checks
        )
        margin = self.margin_deg
        runner_text = (
            ""
            if self.runner_up is None or margin is None
            else (
                f" The runner-up is variant {self.runner_up.variant_index} at "
                f"{self.runner_up_misfit_deg:.2f} deg, {margin:.2f} deg behind."
            )
        )
        verdict_text = {
            "assigned": "The assignment is unambiguous at this tolerance.",
            "ambiguous": (
                "More than one variant lies inside the tolerance, so the assignment is ambiguous."
            ),
            "no_fit": (
                "No variant lies inside the tolerance: the pair does not obey this "
                "relationship, or the grains are not parent and child."
            ),
        }[self.verdict]
        return (
            f"Variant identification under '{relationship.name}' ({parent} to {child}): of "
            f"{len(self.variants)} variants, variant {self.variant_index} predicts the measured "
            f"child to within {self.misfit_deg:.3f} deg (the smallest rotation taking the "
            "predicted child orientation onto the measured one, minimized over the child point "
            f"group).{runner_text} {verdict_text} Measured on the two grains: {statement}. "
            "Orientations are crystal-to-specimen matrices, the measured rotation is "
            "V_m = C^T P, and variants are V_k = R S_p^T in PyTex generation order."
        )

    def to_json_dict(self) -> dict[str, object]:
        """A JSON-ready record of the identification, in lockstep with :meth:`describe`."""

        return {
            "relationship": self.relationship.name,
            "parent_phase": self.relationship.parent_phase.name,
            "child_phase": self.relationship.child_phase.name,
            "variant_index": self.variant_index,
            "misfit_deg": self.misfit_deg,
            "runner_up_variant_index": (
                None if self.runner_up is None else int(self.runner_up.variant_index)
            ),
            "runner_up_misfit_deg": self.runner_up_misfit_deg,
            "margin_deg": self.margin_deg,
            "tolerance_deg": float(self.tolerance_deg),
            "verdict": self.verdict,
            "parent_matrix": self.parent_matrix.tolist(),
            "child_matrix": self.child_matrix.tolist(),
            "measured_rotation": self.measured_rotation.tolist(),
            "variants": [
                {
                    "variant_index": int(variant.variant_index),
                    "parent_operator_index": int(variant.parent_operator_index),
                    "misfit_deg": float(self.misfits_deg[position]),
                    "child_operator_index": int(self.child_operator_indices[position]),
                }
                for position, variant in enumerate(self.variants)
            ],
            "parallelism_checks": [
                {
                    "kind": check.kind,
                    "parent": _indices_text(check.parent),
                    "child": _indices_text(check.child),
                    "angle_deg": float(check.angle_deg),
                }
                for check in self.parallelism_checks
            ],
            "describe": self.describe(),
        }


def _operators(symmetry: object) -> np.ndarray:
    if symmetry is None:
        return np.eye(3, dtype=np.float64)[None, :, :]
    return np.asarray(symmetry.operators, dtype=np.float64)  # type: ignore[attr-defined]


def _indices_text(member: CrystalPlane | CrystalDirection) -> str:
    if isinstance(member, CrystalPlane):
        values = member.miller.indices
        return "(" + " ".join(str(int(value)) for value in values) + ")"
    coordinates = np.asarray(member.coordinates, dtype=np.float64)
    return "[" + " ".join(f"{value:g}" for value in coordinates) + "]"


def _require_orientation(
    orientation: Orientation, relationship: OrientationRelationship, *, role: str
) -> None:
    phase = relationship.parent_phase if role == "parent" else relationship.child_phase
    if orientation.phase is not None and not phases_semantically_match(orientation.phase, phase):
        raise ValueError(
            f"The {role} orientation's phase does not match the relationship's {role} phase."
        )


def identify_variant(
    parent: Orientation,
    child: Orientation,
    relationship: OrientationRelationship,
    *,
    tolerance_deg: float = DEFAULT_ASSIGNMENT_TOLERANCE_DEG,
) -> VariantIdentificationReport:
    """Identify which variant of ``relationship`` the measured child grain is.

    Purpose
    -------
    The single-pair, fully explained form of variant assignment: every
    intermediate a report needs to show the calculation is kept.

    Parameters
    ----------
    parent, child : Orientation
        Measured orientations of a parent grain and of a child grain, in one
        specimen frame. Their phases, when set, must be the relationship's.
    relationship : OrientationRelationship
        The relationship the pair is tested against.
    tolerance_deg : float
        :math:`\\tau`: the misfit inside which a variant is accepted.

    Returns
    -------
    VariantIdentificationReport

    Raises
    ------
    ValueError
        If the specimen frames differ, a phase does not match, or the
        tolerance is not a positive finite angle.

    Examples
    --------
    A child built as variant 4 of Kurdjumov-Sachs, with the child's
    description changed by one of its own symmetry operators, is identified as
    variant 4 at zero misfit; see ``tests/unit/test_variant_identification.py``.
    """

    if parent.specimen_frame != child.specimen_frame:
        raise ValueError("The parent and child orientations must share a specimen frame.")
    _require_orientation(parent, relationship, role="parent")
    _require_orientation(child, relationship, role="child")
    if not np.isfinite(tolerance_deg) or tolerance_deg <= 0.0:
        raise ValueError("tolerance_deg must be a positive finite angle.")

    parent_matrix = np.asarray(parent.as_matrix(), dtype=np.float64)
    child_matrix = np.asarray(child.as_matrix(), dtype=np.float64)
    measured = child_matrix.T @ parent_matrix

    variants = relationship.generate_variants()
    variant_matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in variants], axis=0
    )
    child_operators = _operators(relationship.child_phase.symmetry)
    # Delta[k, s] = V_k V_m^T S_c[s]: the residual of variant k under child operator s.
    residuals = np.einsum(
        "kij,lj,slm->ksim", variant_matrices, measured, child_operators, optimize=True
    )
    angles = np.degrees(rotation_angle_from_matrix_rad(residuals))
    best_operators = np.argmin(angles, axis=1)
    misfits = angles[np.arange(len(variants)), best_operators]
    order = np.argsort(misfits, kind="stable")
    best_position = int(order[0])
    runner_up_position = int(order[1]) if len(variants) > 1 else None

    variant = variants[best_position]
    child_operator = child_operators[int(best_operators[best_position])]
    checks: list[ParallelismCheck] = []
    for parent_plane, child_plane in variant.parallel_planes:
        measured_child_plane = _symmetry_image_plane(child_plane, child_operator)
        checks.append(
            _check("plane", parent_plane, measured_child_plane, parent_matrix, child_matrix)
        )
    for parent_direction, child_direction in variant.parallel_directions:
        measured_child_direction = _symmetry_image_direction(child_direction, child_operator)
        checks.append(
            _check(
                "direction",
                parent_direction,
                measured_child_direction,
                parent_matrix,
                child_matrix,
            )
        )

    return VariantIdentificationReport(
        relationship=relationship,
        parent_matrix=parent_matrix,
        child_matrix=child_matrix,
        measured_rotation=measured,
        variants=variants,
        misfits_deg=misfits,
        child_operator_indices=best_operators,
        best_position=best_position,
        runner_up_position=runner_up_position,
        tolerance_deg=float(tolerance_deg),
        parallelism_checks=tuple(checks),
    )


def _unit(member: CrystalPlane | CrystalDirection) -> np.ndarray:
    vector = member.normal if isinstance(member, CrystalPlane) else member.unit_vector
    array = np.asarray(vector, dtype=np.float64)
    return array / np.linalg.norm(array)


def _check(
    kind: Literal["plane", "direction"],
    parent_member: CrystalPlane | CrystalDirection,
    child_member: CrystalPlane | CrystalDirection,
    parent_matrix: np.ndarray,
    child_matrix: np.ndarray,
) -> ParallelismCheck:
    parent_specimen = parent_matrix @ _unit(parent_member)
    child_specimen = child_matrix @ _unit(child_member)
    angle = float(np.degrees(angle_between_unit_vectors_rad(parent_specimen, child_specimen)))
    for array in (parent_specimen, child_specimen):
        array.setflags(write=False)
    return ParallelismCheck(
        kind=kind,
        parent=parent_member,
        child=child_member,
        parent_specimen=parent_specimen,
        child_specimen=child_specimen,
        angle_deg=angle,
    )


def operator_axis_angle(operator: np.ndarray) -> tuple[np.ndarray, float]:
    """Axis (unit, crystal frame) and angle in degrees of a point-group operator.

    The identity has no axis; it is returned as ``[0, 0, 1]`` with angle 0 so a
    caller can print it without a special case.
    """

    rotation = Rotation.from_matrix(np.asarray(operator, dtype=np.float64))
    angle = float(rotation.angle_deg)
    if angle < 1e-9:
        return np.array([0.0, 0.0, 1.0]), 0.0
    return np.asarray(rotation.axis, dtype=np.float64), angle
