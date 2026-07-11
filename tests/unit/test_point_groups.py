from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    PointGroup,
    SymmetrySpec,
    all_point_group_symbols,
    laue_class_symbols,
    normalize_point_group_symbol,
)

EXPECTED_ORDERS = {
    "1": 1,
    "-1": 2,
    "2": 2,
    "m": 2,
    "2/m": 4,
    "222": 4,
    "mm2": 4,
    "mmm": 8,
    "4": 4,
    "-4": 4,
    "4/m": 8,
    "422": 8,
    "4mm": 8,
    "-42m": 8,
    "4/mmm": 16,
    "3": 3,
    "-3": 6,
    "32": 6,
    "3m": 6,
    "-3m": 12,
    "6": 6,
    "-6": 6,
    "6/m": 12,
    "622": 12,
    "6mm": 12,
    "-6m2": 12,
    "6/mmm": 24,
    "23": 12,
    "m-3": 24,
    "432": 24,
    "-43m": 24,
    "m-3m": 48,
}

EXPECTED_LAUE = {
    "1": "-1",
    "-1": "-1",
    "2": "2/m",
    "m": "2/m",
    "2/m": "2/m",
    "222": "mmm",
    "mm2": "mmm",
    "mmm": "mmm",
    "4": "4/m",
    "-4": "4/m",
    "4/m": "4/m",
    "422": "4/mmm",
    "4mm": "4/mmm",
    "-42m": "4/mmm",
    "4/mmm": "4/mmm",
    "3": "-3",
    "-3": "-3",
    "32": "-3m",
    "3m": "-3m",
    "-3m": "-3m",
    "6": "6/m",
    "-6": "6/m",
    "6/m": "6/m",
    "622": "6/mmm",
    "6mm": "6/mmm",
    "-6m2": "6/mmm",
    "6/mmm": "6/mmm",
    "23": "m-3",
    "m-3": "m-3",
    "432": "m-3m",
    "-43m": "m-3m",
    "m-3m": "m-3m",
}

EXPECTED_PROPER_SUBGROUP = {
    "1": "1",
    "-1": "1",
    "2": "2",
    "m": "1",
    "2/m": "2",
    "222": "222",
    "mm2": "2",
    "mmm": "222",
    "4": "4",
    "-4": "2",
    "4/m": "4",
    "422": "422",
    "4mm": "4",
    "-42m": "222",
    "4/mmm": "422",
    "3": "3",
    "-3": "3",
    "32": "32",
    "3m": "3",
    "-3m": "32",
    "6": "6",
    "-6": "3",
    "6/m": "6",
    "622": "622",
    "6mm": "6",
    "-6m2": "32",
    "6/mmm": "622",
    "23": "23",
    "m-3": "23",
    "432": "432",
    "-43m": "23",
    "m-3m": "432",
}

EXPECTED_SCHOENFLIES = {
    "1": "C1",
    "-1": "Ci",
    "2": "C2",
    "m": "Cs",
    "2/m": "C2h",
    "222": "D2",
    "mm2": "C2v",
    "mmm": "D2h",
    "4": "C4",
    "-4": "S4",
    "4/m": "C4h",
    "422": "D4",
    "4mm": "C4v",
    "-42m": "D2d",
    "4/mmm": "D4h",
    "3": "C3",
    "-3": "C3i",
    "32": "D3",
    "3m": "C3v",
    "-3m": "D3d",
    "6": "C6",
    "-6": "C3h",
    "6/m": "C6h",
    "622": "D6",
    "6mm": "C6v",
    "-6m2": "D3h",
    "6/mmm": "D6h",
    "23": "T",
    "m-3": "Th",
    "432": "O",
    "-43m": "Td",
    "m-3m": "Oh",
}


def test_all_32_point_groups_are_defined() -> None:
    symbols = all_point_group_symbols()
    assert len(symbols) == 32
    assert set(symbols) == set(EXPECTED_ORDERS)


@pytest.mark.parametrize("symbol", sorted(EXPECTED_ORDERS))
def test_group_order_matches_international_tables(symbol: str) -> None:
    group = PointGroup.from_symbol(symbol)
    assert group.order == EXPECTED_ORDERS[symbol]


@pytest.mark.parametrize("symbol", sorted(EXPECTED_ORDERS))
def test_group_operators_are_closed_orthogonal_set(symbol: str) -> None:
    group = PointGroup.from_symbol(symbol)
    operators = group.operators
    keys = {tuple(np.round(operator, decimals=8).ravel()) for operator in operators}
    assert len(keys) == group.order
    for left in operators:
        assert_allclose(left.T @ left, np.eye(3), atol=1e-10)
        for right in operators:
            product_key = tuple(np.round(left @ right, decimals=8).ravel())
            assert product_key in keys


@pytest.mark.parametrize("symbol", sorted(EXPECTED_LAUE))
def test_laue_class_assignment(symbol: str) -> None:
    group = PointGroup.from_symbol(symbol)
    assert group.laue_class_symbol == EXPECTED_LAUE[symbol]
    laue = group.laue_class()
    assert laue.is_laue
    assert laue.is_centrosymmetric
    assert laue.order == 2 * laue.proper_subgroup().order


@pytest.mark.parametrize("symbol", sorted(EXPECTED_PROPER_SUBGROUP))
def test_proper_subgroup_contains_exactly_the_rotations(symbol: str) -> None:
    group = PointGroup.from_symbol(symbol)
    assert group.proper_subgroup_symbol == EXPECTED_PROPER_SUBGROUP[symbol]
    subgroup = group.proper_subgroup()
    assert subgroup.is_proper
    assert subgroup.order == group.rotations.shape[0]
    subgroup_keys = {
        tuple(np.round(operator, decimals=8).ravel()) for operator in subgroup.operators
    }
    rotation_keys = {
        tuple(np.round(operator, decimals=8).ravel()) for operator in group.rotations
    }
    assert subgroup_keys == rotation_keys


@pytest.mark.parametrize("symbol", sorted(EXPECTED_SCHOENFLIES))
def test_schoenflies_names_round_trip(symbol: str) -> None:
    group = PointGroup.from_symbol(symbol)
    assert group.schoenflies == EXPECTED_SCHOENFLIES[symbol]
    assert normalize_point_group_symbol(group.schoenflies) == symbol
    assert normalize_point_group_symbol(group.schoenflies.lower()) == symbol


def test_symbol_aliases_and_rejection() -> None:
    assert normalize_point_group_symbol("m3m") == "m-3m"
    assert normalize_point_group_symbol("m3") == "m-3"
    assert normalize_point_group_symbol("43m") == "-43m"
    assert normalize_point_group_symbol("s6") == "-3"
    assert normalize_point_group_symbol("3/m") == "-6"
    assert normalize_point_group_symbol(" M-3M ") == "m-3m"
    with pytest.raises(ValueError):
        normalize_point_group_symbol("5")


def test_laue_class_symbols_are_the_11_laue_groups() -> None:
    laue = laue_class_symbols()
    assert len(laue) == 11
    assert set(laue) == {
        "-1",
        "2/m",
        "mmm",
        "4/m",
        "4/mmm",
        "-3",
        "-3m",
        "6/m",
        "6/mmm",
        "m-3",
        "m-3m",
    }


def test_centrosymmetric_groups_are_exactly_the_laue_groups() -> None:
    for symbol in all_point_group_symbols():
        group = PointGroup.from_symbol(symbol)
        assert group.is_centrosymmetric == group.is_laue


def test_mirror_counts_for_reference_groups() -> None:
    assert PointGroup.from_symbol("mmm").mirror_normals().shape[0] == 3
    assert PointGroup.from_symbol("4mm").mirror_normals().shape[0] == 4
    assert PointGroup.from_symbol("6mm").mirror_normals().shape[0] == 6
    assert PointGroup.from_symbol("-6m2").mirror_normals().shape[0] == 4
    assert PointGroup.from_symbol("-43m").mirror_normals().shape[0] == 6
    assert PointGroup.from_symbol("m-3m").mirror_normals().shape[0] == 9
    assert PointGroup.from_symbol("432").mirror_normals().shape[0] == 0


def test_minus_6_contains_horizontal_mirror() -> None:
    group = PointGroup.from_symbol("-6")
    normals = group.mirror_normals()
    assert normals.shape[0] == 1
    assert_allclose(np.abs(normals[0]), [0.0, 0.0, 1.0], atol=1e-10)


def test_equivalent_directions_distinguish_improper_coverage() -> None:
    assert PointGroup.from_symbol("m-3m").equivalent_directions([1, 1, 1]).shape[0] == 8
    assert PointGroup.from_symbol("-43m").equivalent_directions([1, 1, 1]).shape[0] == 4
    assert PointGroup.from_symbol("23").equivalent_directions([1, 1, 1]).shape[0] == 4
    assert PointGroup.from_symbol("m-3").equivalent_directions([1, 1, 1]).shape[0] == 8
    assert PointGroup.from_symbol("m-3m").equivalent_directions([1, 0, 0]).shape[0] == 6
    assert (
        PointGroup.from_symbol("23").equivalent_directions([1, 1, 1], antipodal=True).shape[0]
        == 8
    )


def test_symmetry_spec_accepts_all_32_symbols_with_proper_subgroup_operators() -> None:
    for symbol in all_point_group_symbols():
        spec = SymmetrySpec.from_point_group(symbol)
        group = PointGroup.from_symbol(symbol)
        assert spec.order == group.proper_subgroup().order
        assert spec.proper_point_group == EXPECTED_PROPER_SUBGROUP[symbol]
        assert spec.laue_group_symbol == EXPECTED_LAUE[symbol]


def test_symmetry_spec_point_group_bridges() -> None:
    spec = SymmetrySpec.from_point_group("-42m")
    assert spec.order == 4
    assert spec.is_laue is False
    group = spec.to_point_group()
    assert group.hermann_mauguin == "-42m"
    assert group.order == 8
    laue_spec = spec.laue_symmetry()
    assert laue_spec.point_group == "4/mmm"
    assert laue_spec.is_laue is True
    assert laue_spec.order == 8
    round_trip = group.to_symmetry_spec()
    assert isinstance(round_trip, SymmetrySpec)
    assert round_trip.proper_point_group == "222"


def test_specimen_symmetry_constructor() -> None:
    triclinic = SymmetrySpec.specimen("triclinic")
    assert triclinic.order == 1
    monoclinic = SymmetrySpec.specimen("monoclinic")
    assert monoclinic.order == 2
    orthorhombic = SymmetrySpec.specimen("orthorhombic")
    assert orthorhombic.order == 4
    orthotropic = SymmetrySpec.specimen("orthotropic")
    assert orthotropic.order == 4
    assert orthotropic.specimen_symmetry == "orthotropic"
    with pytest.raises(ValueError):
        SymmetrySpec.specimen("axial")
