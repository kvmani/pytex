from __future__ import annotations

from numpy.testing import assert_array_equal

from pytex.core import (
    direction_uvtw_to_uvw,
    direction_uvw_to_uvtw,
    plane_hkil_to_hkl,
    plane_hkl_to_hkil,
)


def test_direction_three_index_to_four_index_for_basal_a_direction() -> None:
    assert_array_equal(direction_uvw_to_uvtw([1, 0, 0]), [2, -1, -1, 0])


def test_direction_four_index_round_trip() -> None:
    four_index = [2, -1, -1, 3]
    three_index = direction_uvtw_to_uvw(four_index)
    assert_array_equal(three_index, [1, 0, 1])
    assert_array_equal(direction_uvw_to_uvtw(three_index), four_index)


def test_plane_three_index_to_four_index_round_trip() -> None:
    plane = [1, 0, 2]
    four_index = plane_hkl_to_hkil(plane)
    assert_array_equal(four_index, [1, 0, -1, 2])
    assert_array_equal(plane_hkil_to_hkl(four_index), plane)
