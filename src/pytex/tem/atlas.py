"""The zone-axis atlas: which axes a phase has, and what each one shows.

Why this module exists
----------------------
:func:`pytex.tem.navigation.plan_tilt_to_zone_axis` answers "can I get to the
axis I named?". At the column the question usually comes the other way round: the
beam is down some axis, the specimen must be characterised, and what is wanted is
the *list* of axes worth going to — which are close, which give a rich pattern,
which are so nearly degenerate with the one already recorded that going there
settles nothing.

That list is a property of the phase and its symmetry, and it is the same list
every time, so it is computed once here rather than assembled by hand from a
textbook table. Each entry carries what decides the choice:

- **the family** it represents, and how many symmetry-equivalent members it has,
  because a family with 12 members offers 12 chances of one being reachable;
- **the angle** from the axis currently on the beam, which is the cost of going;
- **how many reflections** the pattern shows inside a stated cut-off, which is
  how much information the trip buys;
- **the rotational symmetry of the pattern**, which is what the microscopist
  recognises on the screen when they arrive, and the first check that they
  arrived where they intended.

Everything here is kinematic and geometric. See
:doc:`/theory/tem_specimen_tilt_navigation` for the tilt geometry that consumes
this, and :doc:`/theory/reciprocal_space_and_kinematic_spots` for the zone law
and the structure factors.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core._arrays import as_int_array
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.notation import format_direction_family_indices, format_direction_indices

__all__ = [
    "ZoneAxisAtlas",
    "ZoneAxisEntry",
    "pattern_rotational_order",
    "zone_axis_atlas",
]

#: Reciprocal-space cut-off used when counting the reflections of a zone.
#:
#: 1.5 Å⁻¹ corresponds to d = 0.667 Å, which is comfortably beyond the useful
#: range of a conventional SAED plate and therefore counts every reflection a
#: microscopist would index while excluding the arbitrarily large ones that a
#: bare index bound would admit.
DEFAULT_CUTOFF_INV_ANGSTROM = 1.5

#: Relative-intensity threshold for a reflection to count as visible.
DEFAULT_INTENSITY_FLOOR = 0.002

#: Angular tolerance when testing a pattern for n-fold rotational symmetry.
_SYMMETRY_TOLERANCE = 1e-6


def _reduce(indices: np.ndarray) -> np.ndarray:
    """Divide out the common factor, so [002] and [001] are one direction."""

    values = np.asarray(indices, dtype=np.int64).reshape(3)
    divisor = int(np.gcd.reduce(np.abs(values)))
    if divisor > 1:
        values = values // divisor
    return values


def _preference_key(indices: np.ndarray) -> tuple[int, int, int, tuple[int, ...], tuple[int, ...]]:
    """Sort key selecting the conventional representative of a direction family.

    The literature writes the ⟨110⟩ family as [110], not as [1̄01̄]: small
    indices, positive where the symmetry allows, and the larger components
    first. This orders an orbit so that ``min`` picks that member. The final
    component breaks the remaining ties towards positive leading indices, so a
    hexagonal family is reported as [310] rather than the equally valid [3̄10].
    """

    values = [int(value) for value in np.asarray(indices).reshape(3)]
    negatives = sum(1 for value in values if value < 0)
    return (
        int(sum(abs(value) for value in values)),
        int(max(abs(value) for value in values)),
        negatives,
        tuple(-abs(value) for value in values),
        tuple(-value for value in values),
    )


def _orbit_indices(phase: Phase, indices: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """The symmetry orbit of a lattice direction, in index space, both senses.

    Symmetry operators act on Cartesian crystal vectors and map the lattice onto
    itself, so pushing the direction through the direct basis, applying every
    operator, and pulling back gives integers exactly. Rounding therefore repairs
    floating-point dust rather than approximating anything.
    """

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    inverse = np.linalg.inv(direct)
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    cartesian = direct @ np.asarray(indices, dtype=float).reshape(3)
    images = np.einsum("nij,j->ni", operators, cartesian)
    both = np.vstack([images, -images])
    lattice_images = np.rint(both @ inverse.T).astype(np.int64)
    unique = {tuple(int(value) for value in _reduce(row)) for row in lattice_images}
    return tuple(sorted(unique))


def pattern_rotational_order(positions: np.ndarray, intensities: np.ndarray) -> int:
    """The n-fold rotational symmetry of a projected spot pattern.

    Purpose
    -------
    What a microscopist reads off the screen the instant a zone axis arrives: an
    fcc [001] pattern is square, a [111] is hexagonal, a bcc [110] is
    rectangular. It is the fastest confirmation that the intended axis, and not
    a neighbour, is on the beam.

    Method
    ------
    Rotate the spot set about the transmitted beam by ``360/n`` for
    ``n = 6, 4, 3, 2`` in turn and return the largest ``n`` that maps the set
    onto itself, matching intensity as well as position. Measured on the pattern
    rather than deduced from the point group, so it reports the symmetry actually
    present — which for a kinematic pattern includes the Friedel centre of
    symmetry whether or not the crystal has one.

    Parameters
    ----------
    positions : np.ndarray
        ``(n, 2)`` in-plane spot coordinates, relative to the transmitted beam.
    intensities : np.ndarray
        ``(n,)`` relative intensities, matched to within one part in a hundred.

    Returns
    -------
    int
        6, 4, 3, 2, or 1 when no rotation above the identity maps the set onto
        itself.
    """

    coordinates = np.asarray(positions, dtype=float).reshape(-1, 2)
    weights = np.asarray(intensities, dtype=float).reshape(-1)
    if coordinates.shape[0] == 0:
        return 1
    scale = float(np.max(np.linalg.norm(coordinates, axis=1)))
    if scale <= 0.0:
        return 1
    tolerance = max(scale * 1e-4, _SYMMETRY_TOLERANCE)
    for order in (6, 4, 3, 2):
        angle = 2.0 * math.pi / order
        rotation = np.array(
            [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
        )
        rotated = coordinates @ rotation.T
        distances = np.linalg.norm(rotated[:, None, :] - coordinates[None, :, :], axis=2)
        nearest = np.argmin(distances, axis=1)
        if not np.all(distances[np.arange(rotated.shape[0]), nearest] <= tolerance):
            continue
        if np.allclose(weights[nearest], weights, rtol=1e-2, atol=1e-6):
            return order
    return 1


@dataclass(frozen=True, slots=True)
class ZoneAxisEntry:
    """One zone-axis family of a phase, described for a navigation decision.

    Attributes
    ----------
    indices : np.ndarray
        The conventional representative of the family, in the phase's
        three-index basis.
    label : str
        The representative written as a specific direction ``[uvw]``.
    family_label : str
        The same axis written as a symmetry family ``⟨uvw⟩``, which is what the
        entry actually stands for: every member gives the same pattern.
    family_size : int
        Number of symmetry-equivalent members, with ``[uvw]`` and ``[ūv̄w̄]``
        counted once because they give the same diffraction pattern.
    angle_from_current_deg : float
        Smallest angle between the axis on the beam and any member of this
        family. ``0`` when this *is* the current family, and ``nan`` when no
        current axis was supplied.
    nearest_member : np.ndarray
        The member achieving that angle — the one worth tilting to.
    reflection_count : int
        Reflections in the zone inside the cut-off, above the intensity floor.
        A measure of how much the pattern has to say.
    largest_d_angstrom : float
        The innermost reflection's spacing: the first ring a microscopist
        measures, and the one most sensitive to a wrong camera constant.
    rotational_order : int
        The pattern's apparent n-fold symmetry. See
        :func:`pattern_rotational_order`.
    """

    indices: np.ndarray
    label: str
    family_label: str
    family_size: int
    angle_from_current_deg: float
    nearest_member: np.ndarray
    reflection_count: int
    largest_d_angstrom: float
    rotational_order: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        object.__setattr__(self, "nearest_member", as_int_array(self.nearest_member, shape=(3,)))

    def to_json(self) -> dict[str, Any]:
        """The entry as JSON-ready data."""

        return {
            "indices": [int(value) for value in self.indices],
            "label": self.label,
            "family_label": self.family_label,
            "family_size": int(self.family_size),
            "angle_from_current_deg": float(self.angle_from_current_deg),
            "nearest_member": [int(value) for value in self.nearest_member],
            "reflection_count": int(self.reflection_count),
            "largest_d_angstrom": float(self.largest_d_angstrom),
            "rotational_order": int(self.rotational_order),
        }


@dataclass(frozen=True, slots=True)
class ZoneAxisAtlas:
    """The zone axes of one phase, ranked for a navigation decision.

    Attributes
    ----------
    phase : Phase
    current_zone_axis : ZoneAxis, optional
        The axis the angles are measured from.
    entries : tuple of ZoneAxisEntry
        Ordered by angle from the current axis when there is one, and by
        reflection count otherwise.
    cutoff_inv_angstrom : float
        The reciprocal-space limit the reflection counts were taken inside.
    max_index : int
        The index bound the families were enumerated within.
    """

    phase: Phase
    current_zone_axis: ZoneAxis | None
    entries: tuple[ZoneAxisEntry, ...]
    cutoff_inv_angstrom: float
    max_index: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "entries", tuple(self.entries))

    def describe(self) -> str:
        """A prose account of the atlas, stating its bounds and its ranking."""

        if not self.entries:
            return (
                f"No zone-axis family of {self.phase.name} within an index bound of "
                f"{self.max_index} shows a reflection inside {self.cutoff_inv_angstrom:g} Å⁻¹."
            )
        richest = max(self.entries, key=lambda entry: entry.reflection_count)
        head = (
            f"{len(self.entries)} zone-axis families of {self.phase.name} within an index bound "
            f"of {self.max_index}. Reflection counts are kinematic and taken inside "
            f"{self.cutoff_inv_angstrom:g} Å⁻¹; the richest is {richest.label} with "
            f"{richest.reflection_count} reflections and {richest.rotational_order}-fold pattern "
            "symmetry."
        )
        if self.current_zone_axis is None:
            return head + " No current axis was given, so the list is ordered by pattern richness."
        current = format_direction_indices(
            tuple(int(value) for value in self.current_zone_axis.indices), style="plain"
        )
        nearest = min(
            (entry for entry in self.entries if entry.angle_from_current_deg > 1e-6),
            key=lambda entry: entry.angle_from_current_deg,
            default=None,
        )
        tail = (
            f" Angles are measured from {current}."
            if nearest is None
            else (
                f" Angles are measured from {current}; the nearest distinct family is "
                f"{nearest.label} at {nearest.angle_from_current_deg:.2f}°."
            )
        )
        return head + tail


def _zone_reflection_statistics(
    phase: Phase,
    axis_cartesian: np.ndarray,
    hkl: np.ndarray,
    g_cartesian: np.ndarray,
    g_magnitude: np.ndarray,
    amplitude: np.ndarray,
    *,
    cutoff_inv_angstrom: float,
    intensity_floor: float,
) -> tuple[int, float, int]:
    """Count, innermost spacing, and rotational order of one zone's pattern.

    The reflection table is enumerated once by the caller and shared across every
    candidate axis; this selects the rows obeying the zone law for one axis. The
    selection is a dot product and a comparison, so scanning a few hundred axes
    costs one pass over the table each rather than a re-enumeration.
    """

    projection = g_cartesian @ np.asarray(axis_cartesian, dtype=float)
    # The zone law is exact in index space, so the tolerance only absorbs
    # floating-point error in the basis product, not a physical relrod width.
    on_zone = np.abs(projection) <= 1e-8 * np.maximum(g_magnitude, 1.0)
    inside = on_zone & (g_magnitude <= cutoff_inv_angstrom) & (g_magnitude > 0.0)
    if not np.any(inside):
        return 0, float("nan"), 1
    selected_intensity = amplitude[inside] ** 2 / (1.0 + g_magnitude[inside] ** 2)
    peak = float(selected_intensity.max())
    if peak <= 0.0:
        return 0, float("nan"), 1
    relative = selected_intensity / peak
    visible = relative >= intensity_floor
    if not np.any(visible):
        return 0, float("nan"), 1
    magnitudes = g_magnitude[inside][visible]
    vectors = g_cartesian[inside][visible]
    del hkl  # kept in the signature so the caller's table stays one object

    # Project onto any orthonormal pair spanning the zone plane; rotational
    # order is a property of the spot arrangement, not of the basis chosen.
    axis = np.asarray(axis_cartesian, dtype=float)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(trial, axis))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    first = np.cross(axis, trial)
    first = first / float(np.linalg.norm(first))
    second = np.cross(axis, first)
    positions = np.column_stack([vectors @ first, vectors @ second])
    order = pattern_rotational_order(positions, relative[visible])
    return int(magnitudes.size), float(1.0 / magnitudes.min()), order


def zone_axis_atlas(
    phase: Phase,
    *,
    current_zone_axis: ZoneAxis | None = None,
    max_index: int = 3,
    cutoff_inv_angstrom: float = DEFAULT_CUTOFF_INV_ANGSTROM,
    intensity_floor: float = DEFAULT_INTENSITY_FLOOR,
    max_angle_deg: float | None = None,
    limit: int = 32,
    reflection_max_index: int = 6,
) -> ZoneAxisAtlas:
    """Enumerate the zone-axis families of a phase, ranked for navigation.

    Purpose
    -------
    Answer "where should I go next?" rather than "can I get to X?". Returns one
    entry per symmetry-distinct zone axis, with the angle from the axis currently
    on the beam, how rich its pattern is, and what its symmetry looks like.

    When and where to use it
    ------------------------
    At the start of a tilting session, to choose a target; after indexing, to
    find a second axis that will resolve an ambiguous pattern; and in teaching,
    where the table replaces a hand-copied list of interplanar angles. Feed the
    chosen entry's ``nearest_member`` to
    :func:`pytex.tem.navigation.plan_tilt_to_zone_axis` to find out whether the
    holder can actually make the move.

    Parameters
    ----------
    phase : Phase
    current_zone_axis : ZoneAxis, optional
        The axis on the beam. When given, entries are ordered by angle from it
        and the family containing it appears first, at zero.
    max_index : int
        Index bound for enumerating axes. 3 covers every axis a conventional
        session uses; raising it adds increasingly obscure high-index axes.
    cutoff_inv_angstrom : float
        Reciprocal-space limit for counting reflections. See
        :data:`DEFAULT_CUTOFF_INV_ANGSTROM`.
    intensity_floor : float
        Relative intensity below which a reflection is not counted as visible.
    max_angle_deg : float, optional
        Drop families farther than this from the current axis. Requires
        ``current_zone_axis``.
    limit : int
        Largest number of entries returned, applied after ranking.
    reflection_max_index : int
        Index bound for the reflection table the counts are taken from.

    Returns
    -------
    ZoneAxisAtlas

    Raises
    ------
    ValueError
        If ``current_zone_axis`` belongs to another phase, ``max_angle_deg`` is
        given without one, or a bound is not positive.

    Notes
    -----
    Reflection counts are kinematic and ignore double diffraction, so a
    diamond-structure zone will be reported with fewer spots than a real plate
    shows. That understates richness; it never invents it.
    """

    from pytex.diffraction.kinematic import electron_structure_factors

    if max_index <= 0 or reflection_max_index <= 0:
        raise ValueError("max_index and reflection_max_index must be strictly positive.")
    if limit <= 0:
        raise ValueError("limit must be strictly positive.")
    if cutoff_inv_angstrom <= 0.0:
        raise ValueError("cutoff_inv_angstrom must be strictly positive.")
    if current_zone_axis is not None and current_zone_axis.phase != phase:
        raise ValueError("current_zone_axis.phase must match phase.")
    if max_angle_deg is not None and current_zone_axis is None:
        raise ValueError("max_angle_deg requires a current_zone_axis to measure from.")

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    reciprocal = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=float)

    values = np.arange(-reflection_max_index, reflection_max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    hkl = grid[np.any(grid != 0, axis=1)]
    g_cartesian = hkl.astype(float) @ reciprocal.T
    g_magnitude = np.linalg.norm(g_cartesian, axis=1)
    keep = g_magnitude <= cutoff_inv_angstrom
    hkl = hkl[keep]
    g_cartesian = g_cartesian[keep]
    g_magnitude = g_magnitude[keep]
    amplitude = np.abs(electron_structure_factors(phase, hkl, g_magnitude))

    axis_values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    axis_grid = np.stack(
        np.meshgrid(axis_values, axis_values, axis_values, indexing="ij"), axis=-1
    ).reshape(-1, 3)
    candidates = axis_grid[np.any(axis_grid != 0, axis=1)]

    current_unit: np.ndarray | None = None
    if current_zone_axis is not None:
        current_unit = np.asarray(current_zone_axis.unit_vector, dtype=float)

    seen: set[tuple[int, ...]] = set()
    entries: list[ZoneAxisEntry] = []
    for candidate in candidates:
        reduced = _reduce(candidate)
        key = tuple(int(value) for value in reduced)
        if key in seen:
            continue
        orbit = _orbit_indices(phase, reduced)
        seen.update(orbit)
        representative = min((np.asarray(member) for member in orbit), key=_preference_key)
        if current_zone_axis is not None:
            # "You are here" should be labelled with the indices the user has in
            # hand. [001] and [100] are one cubic family, and telling someone
            # looking down [001] that they are at [100] is a needless puzzle.
            current_reduced = _reduce(np.asarray(current_zone_axis.indices, dtype=np.int64))
            if tuple(int(value) for value in current_reduced) in orbit:
                representative = current_reduced

        members = np.asarray(orbit, dtype=float) @ direct.T
        norms = np.linalg.norm(members, axis=1)
        members = members[norms > 0.0] / norms[norms > 0.0, None]

        angle = float("nan")
        nearest = representative
        if current_unit is not None:
            cosines = np.abs(members @ current_unit)
            best = int(np.argmax(np.clip(cosines, -1.0, 1.0)))
            angle = float(math.degrees(math.acos(float(np.clip(cosines[best], -1.0, 1.0)))))
            nearest = np.asarray(orbit[best], dtype=np.int64)
            # The tolerance is not cosmetic. ⟨110⟩ is at exactly 45° from
            # ⟨001⟩, and computing that through a basis product lands a few
            # ulps either side of it, so a bare comparison drops the single
            # most-wanted target from a 45° search about half the time.
            if max_angle_deg is not None and angle > float(max_angle_deg) + 1e-9:
                continue

        axis_cartesian = direct @ representative.astype(float)
        axis_cartesian = axis_cartesian / float(np.linalg.norm(axis_cartesian))
        count, largest_d, order = _zone_reflection_statistics(
            phase,
            axis_cartesian,
            hkl,
            g_cartesian,
            g_magnitude,
            amplitude,
            cutoff_inv_angstrom=cutoff_inv_angstrom,
            intensity_floor=intensity_floor,
        )
        if count == 0:
            continue
        entries.append(
            ZoneAxisEntry(
                indices=representative,
                label=format_direction_indices(
                    tuple(int(value) for value in representative), style="plain"
                ),
                family_label=format_direction_family_indices(
                    tuple(int(value) for value in representative), style="plain"
                ),
                # Both senses give the same pattern, so the orbit — which holds
                # them both — is counted in pairs.
                family_size=max(len(orbit) // 2, 1),
                angle_from_current_deg=angle,
                nearest_member=nearest,
                reflection_count=count,
                largest_d_angstrom=largest_d,
                rotational_order=order,
            )
        )

    if current_unit is None:
        entries.sort(key=lambda entry: (-entry.reflection_count, _preference_key(entry.indices)))
    else:
        entries.sort(
            key=lambda entry: (
                entry.angle_from_current_deg,
                -entry.reflection_count,
                _preference_key(entry.indices),
            )
        )
    return ZoneAxisAtlas(
        phase=phase,
        current_zone_axis=current_zone_axis,
        entries=tuple(entries[:limit]),
        cutoff_inv_angstrom=float(cutoff_inv_angstrom),
        max_index=int(max_index),
    )
