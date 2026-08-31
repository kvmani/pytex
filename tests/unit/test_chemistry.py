"""Element data, which pymatgen is the single authority for.

pymatgen is a required dependency, so these lookups answer the same way on
every machine PyTex is installed on. That is the point of this module, and it
is what the arrangement it replaced could not promise: a hand-maintained
fallback table was consulted wherever the optional extra was absent, it held Ga
and Ge but not As, and so simulating a CBED pattern of gallium arsenide -- the
textbook non-centrosymmetric case -- raised `ValueError` on the deployed office
server while passing on every developer machine.

There is no longer a second table to disagree with the first.
"""

from __future__ import annotations

import pytest

from pytex.core._chemistry import (
    _UNTABULATED_COVALENT_RADIUS_ANGSTROM,
    atomic_number,
    atomic_radius_angstrom,
    covalent_radius_angstrom,
    cpk_color,
    normalize_species_symbol,
    van_der_waals_radius_angstrom,
)

#: Symbol and Z for every element, checked end to end rather than spot-checked.
#: A partial expectation is a list of the elements someone happened to need,
#: and the next element anybody needs is the one that breaks.
_PERIODIC_TABLE = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def test_every_element_has_its_atomic_number() -> None:
    """All 118, in order, because an atomic number is definitional."""

    assert len(_PERIODIC_TABLE) == 118
    assert [atomic_number(symbol) for symbol in _PERIODIC_TABLE] == list(range(1, 119))


@pytest.mark.parametrize("species", ["As", "As3-", "As5+", "As1", "As_2"])
def test_arsenic_resolves_however_the_site_labels_it(species: str) -> None:
    """The regression case, in the label forms a CIF actually writes.

    Site labels carry charges and site suffixes; the lookup is keyed by
    element, so it must reduce the label first.
    """

    assert normalize_species_symbol(species) == "As"
    assert atomic_number(species) == 33


def test_gallium_arsenide_resolves_as_a_pair() -> None:
    """Both species of the structure whose failure motivated this change."""

    assert atomic_number("Ga") == 31
    assert atomic_number("As") == 33
    assert cpk_color("Ga") != cpk_color("As")


def test_an_unparseable_species_is_rejected() -> None:
    with pytest.raises(ValueError, match="element symbol"):
        normalize_species_symbol("3+")


def test_a_symbol_that_is_not_an_element_is_rejected() -> None:
    """`Xx` parses as a symbol and is not an element; the error says so."""

    with pytest.raises(ValueError, match="No atomic number"):
        atomic_number("Xx")


def test_covalent_radii_cover_the_structural_elements() -> None:
    """pymatgen's Cordero table, not a sixteen-entry list of favourites.

    The elements below span the table well past where the old hand-maintained
    list stopped; every one of them must carry a tabulated value rather than
    the trans-curium default.
    """

    for species in ("H", "C", "O", "Na", "Cl", "Fe", "Ga", "As", "Zr", "W", "Pb", "U"):
        radius = covalent_radius_angstrom(species)
        assert radius > 0.0
        assert radius != _UNTABULATED_COVALENT_RADIUS_ANGSTROM, species


def test_the_trans_curium_elements_fall_back_to_a_documented_default() -> None:
    """pymatgen's covalent-radius table stops at curium; the sizes stay sane."""

    assert covalent_radius_angstrom("Og") == _UNTABULATED_COVALENT_RADIUS_ANGSTROM
    assert atomic_radius_angstrom("Og") > 0.0
    assert van_der_waals_radius_angstrom("Og") > 0.0
