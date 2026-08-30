"""The chemistry fallbacks, which are what runs where pymatgen is absent.

pymatgen is an optional extra. It is not installed in the base test lane and it
is deliberately not deployed to the office server, so every lookup in this
module has to work without it. That is not a degraded mode used by nobody: it
is the mode the shipped application runs in.

The failure this guards against had exactly that shape. The fallback atomic
numbers held Ga and Ge but not As, so simulating a CBED pattern of gallium
arsenide -- the textbook non-centrosymmetric case -- raised `ValueError` on
every machine without pymatgen while passing on every machine with it.
"""

from __future__ import annotations

import pytest

from pytex.core._chemistry import (
    _FALLBACK_ATOMIC_NUMBERS,
    atomic_number,
    covalent_radius_angstrom,
    normalize_species_symbol,
)


def test_the_fallback_covers_every_element() -> None:
    """All 118, because an atomic number is definitional.

    A partial table is a list of the elements someone happened to need, and the
    next element anybody needs is the one that breaks. Completeness is checked
    rather than trusted: the values must be exactly 1 through 118, which fails
    on a duplicate, a gap, or a typo in any single entry.
    """

    assert len(_FALLBACK_ATOMIC_NUMBERS) == 118
    assert sorted(_FALLBACK_ATOMIC_NUMBERS.values()) == list(range(1, 119))
    # Spot values across the table, including the one whose absence was the bug.
    assert _FALLBACK_ATOMIC_NUMBERS["H"] == 1
    assert _FALLBACK_ATOMIC_NUMBERS["As"] == 33
    assert _FALLBACK_ATOMIC_NUMBERS["U"] == 92
    assert _FALLBACK_ATOMIC_NUMBERS["Og"] == 118


def test_every_element_resolves_without_pymatgen(monkeypatch: pytest.MonkeyPatch) -> None:
    """The lookup itself, with the optional dependency made unavailable.

    Checking the table alone would not catch a lookup that stopped consulting
    it, so the import is broken on purpose and every element asked for.
    """

    import builtins

    real_import = builtins.__import__

    def refuse_pymatgen(name: str, *args: object, **kwargs: object) -> object:
        if name.startswith("pymatgen"):
            raise ImportError("pymatgen is not installed in this lane")
        return real_import(name, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(builtins, "__import__", refuse_pymatgen)

    for symbol, expected in _FALLBACK_ATOMIC_NUMBERS.items():
        assert atomic_number(symbol) == expected
    # And a covalent radius still comes back for an element the radius table
    # does not name, because that table carries a documented default.
    assert covalent_radius_angstrom("As") > 0.0


def test_the_fallback_agrees_with_pymatgen_where_it_is_installed() -> None:
    """Two sources of one definitional quantity must not disagree.

    Skipped where the optional extra is absent, which is the lane this file
    otherwise exists for.
    """

    element = pytest.importorskip(
        "pymatgen.core.periodic_table", reason="the cross-check needs the 'adapters' extra"
    ).Element
    for symbol, expected in _FALLBACK_ATOMIC_NUMBERS.items():
        try:
            reference = int(element(symbol).Z)
        except Exception:  # pragma: no cover - pymatgen not knowing an element
            continue
        assert reference == expected, symbol


def test_species_labels_are_normalised_before_lookup() -> None:
    """Site labels carry charges and suffixes; the table is keyed by element."""

    assert normalize_species_symbol("As3-") == "As"
    assert atomic_number("As3-") == 33
