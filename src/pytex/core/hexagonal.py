from __future__ import annotations

from collections.abc import Iterable
from fractions import Fraction
from math import gcd

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import IntArray, as_int_array


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
    u4, v4, t4, w4 = (int(value) for value in as_int_array(uvtw, shape=(4,)))
    if u4 + v4 + t4 != 0:
        raise ValueError("Hexagonal four-index directions must satisfy U + V + T = 0.")
    three_index = (2 * u4 + v4, 2 * v4 + u4, w4)
    return _reduce_integer_tuple(three_index)


def plane_hkl_to_hkil(hkl: ArrayLike) -> IntArray:
    h, k, ell = (int(value) for value in as_int_array(hkl, shape=(3,)))
    return as_int_array((h, k, -(h + k), ell), shape=(4,))


def plane_hkil_to_hkl(hkil: ArrayLike) -> IntArray:
    h, k, i, ell = (int(value) for value in as_int_array(hkil, shape=(4,)))
    if i != -(h + k):
        raise ValueError("Hexagonal four-index planes must satisfy i = -(h + k).")
    return as_int_array((h, k, ell), shape=(3,))
