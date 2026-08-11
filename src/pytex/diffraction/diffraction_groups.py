"""Diffraction groups: the crystal point group from a CBED pattern.

Kinematic diffraction cannot see a centre of symmetry. Friedel's law makes
:math:`I_{g} = I_{-g}` whatever the structure, so a selected-area pattern
determines the Laue class and stops there — 11 possibilities where there are 32.
Convergent-beam diffraction breaks that limit, and doing so is its most
celebrated capability: from the symmetry of one zone-axis exposure the full point
group follows, centre included, and with it the answer to whether a crystal can
be piezoelectric, ferroelectric, or optically active.

This module implements the group theory of that determination, following the
formulation of Buxton, Eades, Steeds and Rackham (1976). The 31 diffraction
groups are **derived here rather than transcribed**, as is the table relating
them to the 32 point groups, so there is no copied list to fall out of step with
the operators PyTex actually stores.

The construction
----------------

Fix a beam direction :math:`\\mathbf{b}`. Each operator :math:`S` of the crystal
point group falls into one of three cases:

- :math:`S\\mathbf{b} = \\mathbf{b}`. The operator is a genuine symmetry of the
  experiment. It contributes its restriction :math:`T = S|_{\\perp}` to the
  plane normal to the beam, **untagged**.
- :math:`S\\mathbf{b} = -\\mathbf{b}`. The operator reverses the beam, so it is
  not a symmetry on its own; combined with the reciprocity theorem it becomes
  one. It contributes :math:`T = S|_{\\perp}` **tagged**, written with the
  subscript :math:`R`.
- Neither. The operator relates different patterns and contributes nothing.

The map :math:`S \\mapsto (S|_{\\perp}, \\text{tag})` is a homomorphism onto a
subgroup of :math:`G_{2}\\times\\mathbb{Z}_{2}`, with :math:`G_{2}` one of the
ten two-dimensional crystallographic point groups. Enumerating those subgroups
gives 10 with no tagged element, 10 of the form :math:`G_{2}\\times\\mathbb{Z}_2`
(written with the suffix :math:`1_{R}`), and 11 that are graphs of a surjection
onto :math:`\\mathbb{Z}_{2}` — 31 in all, which is Buxton's count, obtained by
construction.

What each group predicts
------------------------

Two of the observable symmetries follow directly from the element list.

**Whole pattern.** A tagged element requires reciprocity, which relates a point
in one disc to a point in another at an incident direction outside the
illumination cone. Only the untagged elements are rigid symmetries of the
recorded pattern, so

.. math::  \\text{WP} = \\{\\,T : (T, \\text{untagged}) \\in D\\,\\}.

**Bright field.** Inside the direct disc the reciprocity displacement vanishes,
because it is proportional to :math:`\\mathbf{g}_{\\perp}` and
:math:`\\mathbf{g} = 0` there. What remains is reciprocity's own inversion of the
incident direction, so a tagged element acts on the bright-field disc as
:math:`-T` rather than :math:`T`:

.. math::  \\text{BF} = \\varphi(D), \\qquad
   \\varphi(T, \\text{untagged}) = T, \\quad
   \\varphi(T, \\text{tagged}) = -T .

:math:`\\varphi` is a homomorphism because :math:`-1` is central in two
dimensions, so BF is again one of the ten plane groups. The consequences are the
familiar ones: a two-fold axis perpendicular to the beam gives a mirror in the
bright-field disc that the whole pattern does not have (:math:`m_{R}`), and a
horizontal mirror gives a two-fold in the bright-field disc alone
(:math:`1_{R}`).

The centre of symmetry
----------------------

The element :math:`2_{R}` needs an operator with :math:`S\\mathbf{b} =
-\\mathbf{b}` and :math:`S|_{\\perp} = -\\mathbf{1}` — that is, an operator
acting as :math:`-1` on the beam direction and on the plane, which is the
inversion and nothing else. Therefore

    :math:`2_{R} \\in D` at **any** beam direction
    :math:`\\iff` the crystal is centrosymmetric,

and the observation that decides it is whether the :math:`+\\mathbf{g}` and
:math:`-\\mathbf{g}` discs are related by a two-fold rotation of the pattern.
That is the same statement `pytex.diffraction.dynamical` reaches from the other
end, where it appears as the symmetry of the propagator; the two must agree, and
the tests check that they do on a controlled pair of structures.

Note that :math:`2_{R}` is invisible in both BF and WP — :math:`\\varphi(2,
\\text{tagged}) = -2 = 1` — which is precisely why the centrosymmetry question
requires the :math:`\\pm\\mathbf{g}` observation and cannot be settled from disc
symmetry alone.

What is not implemented
-----------------------

Buxton's table also lists the symmetries of dark-field discs and
:math:`\\pm\\mathbf{g}` pairs for reflections lying on symmetry lines, recorded
with that reflection at its own Bragg condition. Those *special* observations
sharpen a determination that BF, WP and the general :math:`\\pm\\mathbf{g}`
relation leave ambiguous. They are not implemented, and
:meth:`PointGroupDetermination.describe` says so instead of presenting a
narrower candidate list than the evidence supports.

See Also
--------
`pytex.diffraction.dynamical` : the forward model that produces the symmetry.
`pytex.core.point_groups` : the 32 crystallographic point groups and operators.
`docs/site/theory/dynamical_cbed_and_symmetry_determination.md` : derivations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.point_groups import PointGroup, all_point_group_symbols
from pytex.core.provenance import ProvenanceRecord

__all__ = [
    "PLANE_POINT_GROUP_SYMBOLS",
    "POINT_GROUP_DETERMINATION_SCHEMA",
    "DiffractionGroup",
    "PointGroupDetermination",
    "SymmetryObservations",
    "ZoneAxisDiffractionGroup",
    "determine_point_group",
    "diffraction_group_for",
    "diffraction_group_for_zone_axis",
    "diffraction_group_symbols",
    "diffraction_group_table",
    "plane_point_group_symbol",
]

#: Schema identifier of the point-group determination payload.
POINT_GROUP_DETERMINATION_SCHEMA = "pytex.cbed_point_group_determination/1"

#: The ten two-dimensional crystallographic point groups, in the order used for
#: reporting: rotation order ascending, mirror-free before mirror-bearing.
PLANE_POINT_GROUP_SYMBOLS = ("1", "2", "3", "4", "6", "m", "2mm", "3m", "4mm", "6mm")

_AXIS_TOLERANCE = 1e-8
_MATRIX_DECIMALS = 6


# --------------------------------------------------------------------------- #
# Two-dimensional group bookkeeping
# --------------------------------------------------------------------------- #


def _matrix_key(matrix: np.ndarray) -> tuple[float, ...]:
    rounded = np.round(np.asarray(matrix, dtype=np.float64).reshape(-1), _MATRIX_DECIMALS)
    # ``+ 0.0`` normalizes negative zero, which would otherwise split a key.
    return tuple(rounded + 0.0)


def _plane_group_symbol(matrices: tuple[np.ndarray, ...]) -> str:
    """The Hermann-Mauguin symbol of a set of 2x2 orthogonal operators.

    Repeated matrices are collapsed first. They occur legitimately: a diffraction
    group of the form ``G x Z_2`` carries every plane operation twice, once
    tagged and once not, and its *spatial* projection is still ``G``.
    """

    distinct: dict[tuple[float, ...], np.ndarray] = {}
    for matrix in matrices:
        distinct.setdefault(_matrix_key(matrix), matrix)
    rotations = [matrix for matrix in distinct.values() if np.linalg.det(matrix) > 0.0]
    mirrors = [matrix for matrix in distinct.values() if np.linalg.det(matrix) < 0.0]
    order = len(rotations)
    if order not in {1, 2, 3, 4, 6}:
        raise ValueError(
            f"{order} rotations do not form a two-dimensional crystallographic point group; "
            "only orders 1, 2, 3, 4 and 6 are possible."
        )
    if not mirrors:
        return str(order)
    if order == 1:
        return "m"
    if order == 3:
        return "3m"
    return f"{order}mm"


def _mirror_line_angle(matrix: np.ndarray) -> float:
    """Angle of a 2D mirror's *line* (not its normal), in ``[0, pi)``."""

    angle = 0.5 * math.atan2(float(matrix[0, 1]), float(matrix[0, 0]))
    return angle % math.pi


def _rotation_angle(matrix: np.ndarray) -> float:
    """Rotation angle of a 2D proper operator, in ``[0, 2 pi)``."""

    return math.atan2(float(matrix[1, 0]), float(matrix[0, 0])) % (2.0 * math.pi)


def plane_point_group_symbol(operators: ArrayLike) -> str:
    """Name the two-dimensional point group formed by a set of plane operations.

    What it does
        Collapses duplicates, counts the proper rotations and asks whether any
        mirror is present, and returns one of the ten symbols in
        :data:`PLANE_POINT_GROUP_SYMBOLS`.

    When to use it
        When a symmetry has been *measured* rather than derived — for example by
        testing which operations leave a simulated CBED disc unchanged — and the
        surviving operations must be named before they can be compared with a
        diffraction group's prediction. This is the bridge
        `pytex.diffraction.cbed.CBEDPattern.symmetry_observations` crosses.

    Parameters
    ----------
    operators:
        ``(n, 2, 2)`` orthogonal matrices. They must form a group; a set that
        does not will produce a rotation count that is not 1, 2, 3, 4 or 6 and
        will raise, which is the intended failure for a measurement that has
        found an inconsistent set.

    Returns
    -------
    str
        The Hermann-Mauguin symbol.

    Raises
    ------
    ValueError
        If the operations do not form a crystallographic plane point group.

    Examples
    --------
    >>> import numpy as np
    >>> plane_point_group_symbol([np.eye(2), -np.eye(2)])
    '2'
    """

    matrices = as_float_array(operators, shape=(None, 2, 2))
    return _plane_group_symbol(tuple(matrices))


# --------------------------------------------------------------------------- #
# The diffraction group
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class DiffractionGroup:
    """One of the 31 diffraction groups, with its elements and its predictions.

    Purpose
    -------
    The symmetry of a convergent-beam pattern, expressed as a group of plane
    operations each carrying a flag saying whether it needs the reciprocity
    theorem. It is the object that stands between a crystal point group and an
    observation: the forward direction predicts what a pattern will look like,
    and the inverse direction is the point-group determination.

    Attributes
    ----------
    symbol : str
        The Buxton symbol, for example ``"4mm1_R"``, ``"2_Rmm_R"``, ``"m_R"``.
        Built from the element list rather than looked up.
    operators : np.ndarray
        ``(order, 2, 2)`` plane operations, in the plane normal to the beam.
    reciprocity_flags : np.ndarray
        ``(order,)`` booleans: ``True`` where the operation needs reciprocity,
        which is where the crystal operator reversed the beam direction.
    projection_symbol : str
        The plane point group of the spatial parts alone, ignoring the flags.
    whole_pattern_symbol : str
        Symmetry of the recorded pattern as a whole: the untagged elements only.
    bright_field_symbol : str
        Symmetry of the direct disc: ``phi(D)`` with tagged elements contributing
        ``-T``.
    """

    symbol: str
    operators: np.ndarray
    reciprocity_flags: np.ndarray
    projection_symbol: str
    whole_pattern_symbol: str
    bright_field_symbol: str

    def __eq__(self, other: object) -> bool:
        # The generated __eq__ would compare ndarrays with ``==`` and raise; the
        # symbol is a complete invariant of the group by construction.
        if not isinstance(other, DiffractionGroup):
            return NotImplemented
        return self.symbol == other.symbol

    def __hash__(self) -> int:
        return hash(self.symbol)

    @property
    def order(self) -> int:
        """Number of elements, tagged and untagged together."""

        return int(self.operators.shape[0])

    @property
    def has_friedel_symmetry(self) -> bool:
        """Whether ``2_R`` is present: the ``+-g`` discs related by a two-fold.

        Purpose
        -------
        The centrosymmetry test. ``2_R`` requires a crystal operator that
        reverses the beam and acts as ``-1`` on the transverse plane, and the
        inversion is the only operator that does both. So this property is true
        for every beam direction of a centrosymmetric crystal and for none of a
        non-centrosymmetric one.
        """

        target = -np.eye(2)
        return any(
            bool(flag) and np.allclose(operator, target, atol=1e-9)
            for operator, flag in zip(self.operators, self.reciprocity_flags, strict=True)
        )

    @property
    def has_projection_reciprocity(self) -> bool:
        """Whether ``1_R`` is present: reciprocity with no spatial operation.

        It arises from a mirror perpendicular to the beam, and its visible
        consequence is a two-fold in the bright-field disc that the whole pattern
        does not share.
        """

        return any(
            bool(flag) and np.allclose(operator, np.eye(2), atol=1e-9)
            for operator, flag in zip(self.operators, self.reciprocity_flags, strict=True)
        )

    def describe(self) -> str:
        """Convention-explicit prose: the group, what it predicts, and the centre."""

        tagged = int(np.count_nonzero(self.reciprocity_flags))
        centre = (
            "It contains 2_R, so the +g and -g discs are related by a two-fold rotation of "
            "the pattern and the crystal is centrosymmetric."
            if self.has_friedel_symmetry
            else (
                "It does not contain 2_R, so the +g and -g discs are *not* related by a "
                "two-fold and the crystal has no centre of symmetry. This is the observation "
                "that kinematic diffraction cannot make."
            )
        )
        return (
            f"Diffraction group {self.symbol}, of order {self.order}: {self.order - tagged} "
            f"operations that are symmetries of the experiment outright and {tagged} that "
            "become symmetries only through the reciprocity theorem, because the crystal "
            "operator they come from reverses the beam direction. It predicts a bright-field "
            f"disc of symmetry {self.bright_field_symbol} and a whole-pattern symmetry of "
            f"{self.whole_pattern_symbol}; these differ whenever a tagged element is present, "
            f"since reciprocity contributes its own inversion inside the direct disc. {centre}"
        )


def _diffraction_group_from_elements(
    operators: list[np.ndarray], flags: list[bool]
) -> DiffractionGroup:
    """Name and characterize a diffraction group from its element list."""

    matrices = np.stack(operators)
    tags = np.asarray(flags, dtype=bool)

    projection = _plane_group_symbol(tuple(matrices))
    untagged = tuple(matrices[~tags])
    whole_pattern = _plane_group_symbol(untagged) if untagged else "1"
    bright = tuple(
        -matrix if tag else matrix for matrix, tag in zip(matrices, tags, strict=True)
    )
    bright_unique: dict[tuple[float, ...], np.ndarray] = {}
    for matrix in bright:
        bright_unique.setdefault(_matrix_key(matrix), matrix)
    bright_field = _plane_group_symbol(tuple(bright_unique.values()))

    symbol = _diffraction_group_symbol(matrices, tags, projection)
    return DiffractionGroup(
        symbol=symbol,
        operators=np.ascontiguousarray(matrices),
        reciprocity_flags=np.ascontiguousarray(tags),
        projection_symbol=projection,
        whole_pattern_symbol=whole_pattern,
        bright_field_symbol=bright_field,
    )


def _diffraction_group_symbol(
    matrices: np.ndarray, tags: np.ndarray, projection: str
) -> str:
    """Assemble the Buxton symbol from the tagged/untagged element structure.

    Three cases exhaust the possibilities, because the tag is a homomorphism onto
    ``Z_2``: no tagged element (the symbol is the plane group), the identity
    tagged (the group is a direct product, suffix ``1_R``), or a proper
    surjection (the generators carry ``_R`` individually).
    """

    if not tags.any():
        return projection
    identity_tagged = any(
        bool(tag) and np.allclose(matrix, np.eye(2), atol=1e-9)
        for matrix, tag in zip(matrices, tags, strict=True)
    )
    if identity_tagged:
        # The trivial plane group is written "1_R", not "11_R": the leading "1"
        # of the projection is implicit in the reciprocity symbol itself.
        return "1_R" if projection == "1" else f"{projection}1_R"

    rotations = [
        (matrix, tag)
        for matrix, tag in zip(matrices, tags, strict=True)
        if np.linalg.det(matrix) > 0.0
    ]
    mirrors = [
        (matrix, tag)
        for matrix, tag in zip(matrices, tags, strict=True)
        if np.linalg.det(matrix) < 0.0
    ]
    order = len(rotations)

    generator_angle = 2.0 * math.pi / order if order > 1 else 0.0
    rotation_tagged = False
    for matrix, tag in rotations:
        if order > 1 and abs(_rotation_angle(matrix) - generator_angle) < 1e-6:
            rotation_tagged = bool(tag)
            break
    rotation_part = ("1" if order == 1 else str(order)) + ("_R" if rotation_tagged else "")

    if not mirrors:
        return rotation_part
    if order == 1:
        return "m" + ("_R" if bool(mirrors[0][1]) else "")

    reference = _mirror_line_angle(mirrors[0][0])
    classes: dict[int, bool] = {}
    for matrix, tag in mirrors:
        offset = (_mirror_line_angle(matrix) - reference) % math.pi
        index = round(offset * order / math.pi) % 2 if order % 2 == 0 else 0
        classes.setdefault(index, bool(tag))

    if order == 3 or len(classes) == 1:
        return rotation_part + "m" + ("_R" if classes[next(iter(classes))] else "")

    ordered = sorted(classes.items(), key=lambda item: (item[1], item[0]))
    suffix = "".join("m" + ("_R" if tag else "") for _, tag in ordered)
    return rotation_part + suffix


def diffraction_group_for(
    point_group: PointGroup | str, beam_direction: ArrayLike
) -> DiffractionGroup:
    """The diffraction group of a point group viewed along a beam direction.

    What it does
        Partitions the point-group operators by their action on the beam
        direction — fixing it, reversing it, or neither — restricts the first two
        classes to the plane normal to the beam, tags the reversing ones with the
        reciprocity flag, and names the resulting group.

    When to use it
        To predict what symmetry a CBED exposure of a known structure will show,
        which is how a simulation is checked and how an experiment is planned;
        and, through :func:`diffraction_group_table`, to find the zone axis that
        separates two candidate point groups.

    Parameters
    ----------
    point_group:
        A :class:`~pytex.core.point_groups.PointGroup` or its Hermann-Mauguin
        symbol. This must be the **true** point group, not the Laue class:
        determining which of them the crystal has is the entire exercise.
    beam_direction:
        ``(3,)`` beam direction in the point group's Cartesian setting. Only the
        direction matters; the sign does not, because the construction treats
        ``+b`` and ``-b`` together by design.

    Returns
    -------
    DiffractionGroup

    Raises
    ------
    ValueError
        If the beam direction is the zero vector.

    Examples
    --------
    A centrosymmetric cubic crystal down a four-fold axis:

    >>> diffraction_group_for("m-3m", [0, 0, 1]).symbol
    '4mm1_R'

    Its non-centrosymmetric relative, along the same direction, loses the
    inversion and with it the two-fold relation between ``+g`` and ``-g``:

    >>> diffraction_group_for("-43m", [0, 0, 1]).has_friedel_symmetry
    False
    """

    group = (
        point_group
        if isinstance(point_group, PointGroup)
        else PointGroup.from_symbol(point_group)
    )
    axis = normalize_vector(as_float_array(beam_direction, shape=(3,)))
    plane = _transverse_basis(axis)

    seen: dict[tuple[float, ...], tuple[np.ndarray, bool]] = {}
    for operator in group.operators:
        image = operator @ axis
        if np.allclose(image, axis, atol=_AXIS_TOLERANCE):
            tagged = False
        elif np.allclose(image, -axis, atol=_AXIS_TOLERANCE):
            tagged = True
        else:
            continue
        transverse = plane.T @ operator @ plane
        seen.setdefault((*_matrix_key(transverse), float(tagged)), (transverse, tagged))

    operators = [entry[0] for entry in seen.values()]
    flags = [entry[1] for entry in seen.values()]
    return _diffraction_group_from_elements(operators, flags)


def _transverse_basis(axis: np.ndarray) -> np.ndarray:
    """``(3, 2)`` orthonormal basis of the plane normal to ``axis``."""

    reference = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    first = normalize_vector(np.cross(reference, axis))
    second = np.cross(axis, first)
    return np.stack([first, second], axis=1)


def diffraction_group_for_zone_axis(
    phase: Phase, zone_axis: ZoneAxis
) -> DiffractionGroup:
    """The diffraction group of a phase viewed down one of its zone axes.

    What it does
        Reads the phase's declared point group and the zone axis's Cartesian
        direction and calls :func:`diffraction_group_for`.

    When to use it
        The practical entry point: it keeps the beam direction in the same
        crystal frame as the rest of `pytex.diffraction`, so the answer refers to
        the same zone axis a pattern would be simulated down.

    Parameters
    ----------
    phase:
        Its ``symmetry.point_group`` must name the crystal's true point group.
        If it names a Laue class, the answer will describe a centrosymmetric
        crystal because that is what was declared, and no calculation can
        recover a centre that was thrown away in the input.
    zone_axis:
        Must belong to the phase.

    Returns
    -------
    DiffractionGroup

    Raises
    ------
    ValueError
        If the zone axis belongs to a different phase.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    return diffraction_group_for(
        phase.symmetry.point_group, np.asarray(zone_axis.unit_vector, dtype=np.float64)
    )


# --------------------------------------------------------------------------- #
# The forward table: which diffraction groups each point group can show
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ZoneAxisDiffractionGroup:
    """One diffraction group a point group can show, and where to look for it.

    Attributes
    ----------
    diffraction_group : DiffractionGroup
    beam_direction : np.ndarray
        ``(3,)`` unit vector, in the point group's Cartesian setting, of a beam
        direction that produces this diffraction group.
    is_special : bool
        Whether the direction lies on a symmetry element of the point group. A
        non-special direction is the generic case, which every crystal shows
        somewhere and which carries the least information.
    """

    diffraction_group: DiffractionGroup
    beam_direction: np.ndarray
    is_special: bool


def _characteristic_directions(group: PointGroup) -> list[np.ndarray]:
    """Beam directions that can produce a distinct diffraction group.

    The special directions of a point group are its symmetry axes and mirror
    normals; a direction *inside* a mirror plane but off every axis is reached by
    summing two of those; and one generic direction covers the rest. Enumerating
    from the group itself avoids the trap of a fixed index list, which misses the
    in-plane two-folds of a hexagonal group entirely because they do not sit at
    integer Cartesian coordinates.
    """

    axes: list[np.ndarray] = [np.array([0.0, 0.0, 1.0])]
    for operator in group.operators:
        eigenvalues, eigenvectors = np.linalg.eig(np.asarray(operator, dtype=np.float64))
        for value, vector in zip(eigenvalues, eigenvectors.T, strict=True):
            if abs(value.imag) > 1e-9 or abs(abs(value.real) - 1.0) > 1e-9:
                continue
            candidate = np.real(vector)
            if float(np.linalg.norm(candidate)) < 1e-9:
                continue
            axes.append(normalize_vector(candidate))

    combined = list(axes)
    for first_index, first in enumerate(axes):
        for second in axes[first_index + 1 :]:
            for signed in (first + second, first - second):
                if float(np.linalg.norm(signed)) > 1e-6:
                    combined.append(normalize_vector(signed))
    combined.append(normalize_vector(np.array([0.1234, 0.2345, 0.4567])))

    unique: dict[tuple[float, ...], np.ndarray] = {}
    for vector in combined:
        signed = vector if _first_significant_component(vector) > 0.0 else -vector
        unique.setdefault(tuple(np.round(signed, 6) + 0.0), signed)
    return list(unique.values())


def _first_significant_component(vector: np.ndarray) -> float:
    for component in vector:
        if abs(component) > 1e-9:
            return float(component)
    return 1.0


@cache
def _table_for_symbol(symbol: str) -> tuple[ZoneAxisDiffractionGroup, ...]:
    group = PointGroup.from_symbol(symbol)
    generic = diffraction_group_for(group, normalize_vector(np.array([0.1234, 0.2345, 0.4567])))
    found: dict[str, ZoneAxisDiffractionGroup] = {}
    for direction in _characteristic_directions(group):
        diffraction = diffraction_group_for(group, direction)
        entry = ZoneAxisDiffractionGroup(
            diffraction_group=diffraction,
            beam_direction=np.ascontiguousarray(direction),
            is_special=diffraction.order > generic.order,
        )
        found.setdefault(diffraction.symbol, entry)
    return tuple(sorted(found.values(), key=lambda item: -item.diffraction_group.order))


def diffraction_group_table(point_group: PointGroup | str) -> tuple[ZoneAxisDiffractionGroup, ...]:
    """Every diffraction group a point group can show, with a beam direction for each.

    What it does
        Scans the characteristic directions of the point group — its rotation
        axes, its mirror normals, directions lying inside mirror planes, and one
        generic direction — and collects the distinct diffraction groups, each
        with a direction that produces it.

    When to use it
        Twice in a determination. Forwards, to choose the zone axis worth
        recording: two candidate point groups that share a diffraction group at
        one zone may differ at another, and this is how that zone is found.
        Backwards, it is the table :func:`determine_point_group` inverts, which
        is why the point-group-to-diffraction-group correspondence in PyTex is
        computed rather than copied.

    Parameters
    ----------
    point_group:
        A point group or its Hermann-Mauguin symbol.

    Returns
    -------
    tuple of ZoneAxisDiffractionGroup
        Ordered by decreasing diffraction-group order, so the most informative
        zone axis comes first.

    Notes
    -----
    The scan is over symmetry-derived directions rather than a list of low
    indices, because a fixed index list misses the in-plane two-fold axes of a
    hexagonal group: they lie at 0, 60 and 120 degrees in the Cartesian setting
    and are not integer triples.
    """

    group = (
        point_group
        if isinstance(point_group, PointGroup)
        else PointGroup.from_symbol(point_group)
    )
    return _table_for_symbol(group.hermann_mauguin)


@cache
def diffraction_group_symbols() -> tuple[str, ...]:
    """The 31 diffraction groups, obtained by scanning all 32 point groups.

    Purpose
    -------
    The count is the check. Buxton, Eades, Steeds and Rackham derived 31
    diffraction groups; this function reaches the same 31 by construction from
    PyTex's own operator tables, so a discrepancy would mean the construction or
    the operators are wrong, not that a transcribed list has a typo.

    Returns
    -------
    tuple of str
        Sorted symbols.
    """

    symbols: set[str] = set()
    for symbol in all_point_group_symbols():
        for entry in diffraction_group_table(symbol):
            symbols.add(entry.diffraction_group.symbol)
    return tuple(sorted(symbols))


# --------------------------------------------------------------------------- #
# The inverse problem
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class SymmetryObservations:
    """What an experimenter reads off a convergent-beam exposure.

    Purpose
    -------
    The input to a point-group determination, stated in the three observations
    that PyTex can predict from first principles. Each may be left unknown, and
    leaving one unknown widens the answer rather than silently assuming a value.

    Attributes
    ----------
    bright_field : str or None
        Symmetry of the direct disc, as a plane point group: one of
        :data:`PLANE_POINT_GROUP_SYMBOLS`.
    whole_pattern : str or None
        Symmetry of the whole recorded pattern, as a plane point group. It is
        always a subgroup of the bright-field symmetry.
    friedel_pair_two_fold : bool or None
        Whether the ``+g`` and ``-g`` discs are related by a two-fold rotation of
        the pattern. This is the ``2_R`` observation, and it alone decides
        centrosymmetry.

    Notes
    -----
    The third observation must be made with higher-order Laue zone detail
    visible. A zeroth-Laue-zone pattern samples the *projected* potential, whose
    symmetry is often higher than the crystal's, and reports a two-fold that the
    crystal does not have — see `pytex.diffraction.dynamical`, where the same
    trap appears as a propagator that is accidentally symmetric.
    """

    bright_field: str | None = None
    whole_pattern: str | None = None
    friedel_pair_two_fold: bool | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("bright_field", self.bright_field),
            ("whole_pattern", self.whole_pattern),
        ):
            if value is not None and value not in PLANE_POINT_GROUP_SYMBOLS:
                supported = ", ".join(PLANE_POINT_GROUP_SYMBOLS)
                raise ValueError(
                    f"{name} must be one of the ten two-dimensional crystallographic point "
                    f"groups ({supported}); got '{value}'."
                )

    @property
    def is_empty(self) -> bool:
        """Whether nothing at all was observed."""

        return (
            self.bright_field is None
            and self.whole_pattern is None
            and self.friedel_pair_two_fold is None
        )

    def matches(self, group: DiffractionGroup) -> bool:
        """Whether a diffraction group is consistent with these observations."""

        if self.bright_field is not None and group.bright_field_symbol != self.bright_field:
            return False
        if self.whole_pattern is not None and group.whole_pattern_symbol != self.whole_pattern:
            return False
        return not (
            self.friedel_pair_two_fold is not None
            and group.has_friedel_symmetry != self.friedel_pair_two_fold
        )

    def describe(self) -> str:
        """Convention-explicit prose: what was recorded and what was left open."""

        parts = []
        parts.append(
            f"bright-field disc symmetry {self.bright_field}"
            if self.bright_field is not None
            else "bright-field disc symmetry not recorded"
        )
        parts.append(
            f"whole-pattern symmetry {self.whole_pattern}"
            if self.whole_pattern is not None
            else "whole-pattern symmetry not recorded"
        )
        if self.friedel_pair_two_fold is None:
            parts.append("the +-g relation not recorded")
        elif self.friedel_pair_two_fold:
            parts.append("the +g and -g discs related by a two-fold")
        else:
            parts.append("the +g and -g discs NOT related by a two-fold")
        return "CBED symmetry observations: " + ", ".join(parts) + "."


@dataclass(frozen=True, slots=True)
class PointGroupDetermination:
    """The point groups consistent with a set of CBED symmetry observations.

    Purpose
    -------
    The result of the determination, carrying the candidates, the diffraction
    groups behind them, the centrosymmetry verdict, and — when the answer is not
    unique — what would narrow it.

    Attributes
    ----------
    observations : SymmetryObservations
    diffraction_groups : tuple of str
        The diffraction-group symbols consistent with the observations.
    point_groups : tuple of str
        The Hermann-Mauguin symbols of the crystal point groups that can produce
        at least one of those diffraction groups at some beam direction.
    is_centrosymmetric : bool or None
        ``True`` or ``False`` when every candidate agrees, ``None`` when they do
        not. This is the answer CBED exists to give.
    provenance : ProvenanceRecord or None
    """

    observations: SymmetryObservations
    diffraction_groups: tuple[str, ...]
    point_groups: tuple[str, ...]
    is_centrosymmetric: bool | None
    provenance: ProvenanceRecord | None = None

    @property
    def is_unique(self) -> bool:
        """Whether exactly one point group survives."""

        return len(self.point_groups) == 1

    def describe(self) -> str:
        """Convention-explicit prose: the verdict, its basis, and what is missing."""

        if not self.point_groups:
            return (
                f"{self.observations.describe()} No crystallographic point group can produce "
                "these observations at any beam direction. Either a symmetry was misread, or "
                "the pattern was recorded off the zone axis, or higher-order Laue zone detail "
                "was absent and the projected symmetry was read as the true one."
            )
        centre = {
            True: (
                "The crystal is centrosymmetric: every candidate contains the inversion, "
                "which is what the two-fold relation between the +g and -g discs reports."
            ),
            False: (
                "The crystal is NOT centrosymmetric: no candidate contains the inversion. "
                "This is the determination kinematic diffraction cannot make, because "
                "Friedel's law imposes the two-fold whatever the structure."
            ),
            None: (
                "The presence of a centre of symmetry is not decided by these observations. "
                "The +-g relation is the observation that decides it, and it must be read "
                "with higher-order Laue zone detail present."
            ),
        }[self.is_centrosymmetric]
        narrowing = (
            ""
            if self.is_unique
            else (
                " To narrow this further, record a second zone axis: two point groups sharing "
                "a diffraction group at one beam direction generally differ at another, and "
                "diffraction_group_table names the directions. Buxton's dark-field and +-g "
                "observations for reflections on symmetry lines would also discriminate, and "
                "are not implemented here."
            )
        )
        return (
            f"{self.observations.describe()} Consistent diffraction groups: "
            f"{', '.join(self.diffraction_groups)}. Consistent crystal point groups: "
            f"{', '.join(self.point_groups)}. {centre}{narrowing}"
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": POINT_GROUP_DETERMINATION_SCHEMA,
            "observations": {
                "bright_field": self.observations.bright_field,
                "whole_pattern": self.observations.whole_pattern,
                "friedel_pair_two_fold": self.observations.friedel_pair_two_fold,
            },
            "diffraction_groups": list(self.diffraction_groups),
            "point_groups": list(self.point_groups),
            "is_centrosymmetric": self.is_centrosymmetric,
            "is_unique": self.is_unique,
        }


def determine_point_group(
    observations: SymmetryObservations,
    *,
    candidate_point_groups: tuple[str, ...] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> PointGroupDetermination:
    """Determine the crystal point group from CBED symmetry observations.

    What it does
        Collects every diffraction group consistent with the observations, then
        every crystal point group that can produce one of them at some beam
        direction, and reports the centrosymmetry verdict on which all the
        survivors agree — or ``None`` when they do not.

    When to use it
        After recording a zone-axis exposure and reading its symmetry: the disc
        symmetry, the pattern symmetry, and whether the ``+g`` and ``-g`` discs
        are related by a two-fold. This is the step that turns those three
        observations into a point group, and it is what distinguishes CBED from
        selected-area diffraction, which cannot get past the Laue class.

    Parameters
    ----------
    observations:
        What was read off the pattern; any field may be left unknown.
    candidate_point_groups:
        Restrict the search, for example to the point groups compatible with an
        already-known Laue class or crystal system. Supplying prior knowledge
        here is usually what makes the answer unique.
    provenance:
        Optional record.

    Returns
    -------
    PointGroupDetermination

    Raises
    ------
    ValueError
        If the observations are entirely empty, since every point group would
        then be consistent and the answer would be no answer at all.

    Notes
    -----
    **Algorithm.**

    1. Enumerate the diffraction groups of every candidate point group at every
       characteristic beam direction (:func:`diffraction_group_table`, cached).
    2. Keep those matching the observations.
    3. Keep the point groups that produced at least one survivor.
    4. Report ``is_centrosymmetric`` when every survivor agrees.

    **The centre of symmetry needs the third observation.** ``2_R`` maps to the
    identity in the bright-field homomorphism and is absent from the whole
    pattern, so disc and pattern symmetry alone can never decide it. Leaving
    ``friedel_pair_two_fold`` unknown therefore leaves the verdict open by
    construction rather than by accident.
    """

    if observations.is_empty:
        raise ValueError(
            "At least one symmetry observation is needed: with none, every point group is "
            "consistent and the determination has no content. The observation that decides "
            "centrosymmetry is friedel_pair_two_fold."
        )

    symbols = candidate_point_groups or all_point_group_symbols()
    matching_diffraction: set[str] = set()
    matching_points: list[str] = []
    for symbol in symbols:
        matched = False
        for entry in diffraction_group_table(symbol):
            if observations.matches(entry.diffraction_group):
                matching_diffraction.add(entry.diffraction_group.symbol)
                matched = True
        if matched:
            matching_points.append(PointGroup.from_symbol(symbol).hermann_mauguin)

    centrosymmetric_flags = {
        PointGroup.from_symbol(symbol).is_centrosymmetric for symbol in matching_points
    }
    verdict = centrosymmetric_flags.pop() if len(centrosymmetric_flags) == 1 else None

    return PointGroupDetermination(
        observations=observations,
        diffraction_groups=tuple(sorted(matching_diffraction)),
        point_groups=tuple(dict.fromkeys(matching_points)),
        is_centrosymmetric=verdict,
        provenance=provenance,
    )
