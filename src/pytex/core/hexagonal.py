from __future__ import annotations

from collections.abc import Iterable
from fractions import Fraction
from math import gcd
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import IntArray, as_int_array

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from pytex.core.lattice import Phase


def _reduce_integer_tuple(values: Iterable[int]) -> IntArray:
    values_tuple = tuple(int(value) for value in values)
    non_zero = [abs(value) for value in values_tuple if value != 0]
    if not non_zero:
        return as_int_array(values_tuple, shape=(len(values_tuple),))
    divisor = non_zero[0]
    for value in non_zero[1:]:
        divisor = gcd(divisor, value)
    reduced = tuple(value // divisor for value in values_tuple)
    return as_int_array(reduced, shape=(len(values_tuple),))


def direction_uvw_to_uvtw(uvw: ArrayLike) -> IntArray:
    """One three-index direction ``[uvw]`` to four-index ``[UVTW]``.

    Applies ``U = (2u - v)/3``, ``V = (2v - u)/3``, ``T = -(U + V)``,
    ``W = w`` in exact rational arithmetic and clears denominators, so the
    result is always an integer quadruple describing the same direction.
    The four-index form exists so that the three symmetry-equivalent
    ``<11-20>`` directions of a hexagonal crystal receive indices that are
    permutations of one another, which the three-index form cannot do.

    See :func:`~pytex.core.miller.direction_uvw_to_uvtw_array` for the
    vectorized form.
    """

    u, v, w = (int(value) for value in as_int_array(uvw, shape=(3,)))
    fractions = (
        Fraction(2 * u - v, 3),
        Fraction(2 * v - u, 3),
        Fraction(-(u + v), 3),
        Fraction(w, 1),
    )
    common_multiple = 1
    for value in fractions:
        common_multiple = int(np.lcm(common_multiple, value.denominator))
    integer_form = tuple(int(value * common_multiple) for value in fractions)
    return _reduce_integer_tuple(integer_form)


def direction_uvtw_to_uvw(uvtw: ArrayLike) -> IntArray:
    """One four-index direction ``[UVTW]`` to three-index ``[uvw]``.

    Applies ``u = 2U + V``, ``v = 2V + U``, ``w = W`` and reduces to lowest
    terms. The redundancy constraint ``U + V + T = 0`` is checked, not
    assumed, so a mistyped quadruple raises instead of silently converting.
    """

    u4, v4, t4, w4 = (int(value) for value in as_int_array(uvtw, shape=(4,)))
    if u4 + v4 + t4 != 0:
        raise ValueError("Hexagonal four-index directions must satisfy U + V + T = 0.")
    three_index = (2 * u4 + v4, 2 * v4 + u4, w4)
    return _reduce_integer_tuple(three_index)


def plane_hkl_to_hkil(hkl: ArrayLike) -> IntArray:
    """One three-index plane ``(hkl)`` to four-index ``(hkil)``.

    Inserts the redundant index ``i = -(h + k)``.
    """

    h, k, ell = (int(value) for value in as_int_array(hkl, shape=(3,)))
    return as_int_array((h, k, -(h + k), ell), shape=(4,))


def plane_hkil_to_hkl(hkil: ArrayLike) -> IntArray:
    """One four-index plane ``(hkil)`` to three-index ``(hkl)``.

    Drops the redundant ``i`` after checking ``i = -(h + k)``.
    """

    h, k, i, ell = (int(value) for value in as_int_array(hkil, shape=(4,)))
    if i != -(h + k):
        raise ValueError("Hexagonal four-index planes must satisfy i = -(h + k).")
    return as_int_array((h, k, ell), shape=(3,))


def is_hexagonal_phase(phase: Phase) -> bool:
    """Whether ``phase`` belongs to the hexagonal crystal system.

    Purpose: decide when to present directions and reflections in four-index
    Miller-Bravais notation (the convention for hexagonal crystals such as
    the alpha-hcp product of the Burgers relationship). Trigonal phases are
    intentionally excluded: PyTex stores them in the hexagonal setting but
    three-index labels remain conventional there.

    This is the single definition used across the library; the diffraction
    package re-exports it rather than defining its own.
    """

    return phase.symmetry.to_point_group().crystal_system == "hexagonal"
