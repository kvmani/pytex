from __future__ import annotations

import re
from collections.abc import Mapping
from functools import lru_cache
from typing import Any

# Jmol/CPK element colors, the palette VESTA and every other crystal viewer
# renders from. This is a *display* convention, not element data: pymatgen has
# no equivalent table, so it is carried here rather than looked up.
#
# Two species must never come out the same colour on a figure, so the table is
# kept complete through the common structural elements rather than partial: a
# missing element used to fall back to one shared grey, and near-duplicates
# (Ni/Cl both green) made two-species structures unreadable.
_FALLBACK_CPK_COLORS = {
    "H": "#ffffff",
    "He": "#d9ffff",
    "Li": "#cc80ff",
    "Be": "#c2ff00",
    "B": "#ffb5b5",
    "C": "#909090",
    "N": "#3050f8",
    "O": "#ff0d0d",
    "F": "#90e050",
    "Ne": "#b3e3f5",
    "Na": "#ab5cf2",
    "Mg": "#8aff00",
    "Al": "#bfa6a6",
    "Si": "#f0c8a0",
    "P": "#ff8000",
    "S": "#ffff30",
    "Cl": "#1ff01f",
    "Ar": "#80d1e3",
    "K": "#8f40d4",
    "Ca": "#3dff00",
    "Sc": "#e6e6e6",
    "Ti": "#bfc2c7",
    "V": "#a6a6ab",
    "Cr": "#8a99c7",
    "Mn": "#9c7ac7",
    "Fe": "#e06633",
    "Co": "#f090a0",
    "Ni": "#50d050",
    "Cu": "#c88033",
    "Zn": "#7d80b0",
    "Ga": "#c28f8f",
    "Ge": "#668f8f",
    "As": "#bd80e3",
    "Se": "#ffa100",
    "Br": "#a62929",
    "Kr": "#5cb8d1",
    "Rb": "#702eb0",
    "Sr": "#00ff00",
    "Y": "#94ffff",
    "Zr": "#94e0e0",
    "Nb": "#73c2c9",
    "Mo": "#54b5b5",
    "Tc": "#3b9e9e",
    "Ru": "#248f8f",
    "Rh": "#0a7d8c",
    "Pd": "#006985",
    "Ag": "#c0c0c0",
    "Cd": "#ffd98f",
    "In": "#a67573",
    "Sn": "#668080",
    "Sb": "#9e63b5",
    "Te": "#d47a00",
    "I": "#940094",
    "Xe": "#429eb0",
    "Cs": "#57178f",
    "Ba": "#00c900",
    "La": "#70d4ff",
    "Ce": "#ffffc7",
    "Nd": "#c7ffc7",
    "Sm": "#8fffc7",
    "Eu": "#61ffc7",
    "Gd": "#45ffc7",
    "Tb": "#30ffc7",
    "Dy": "#1fffc7",
    "Er": "#00e675",
    "Yb": "#00bf38",
    "Lu": "#00ab24",
    "Hf": "#4dc2ff",
    "Ta": "#4da6ff",
    "W": "#2194d6",
    "Re": "#267dab",
    "Os": "#266696",
    "Ir": "#175487",
    "Pt": "#d0d0e0",
    "Au": "#ffd123",
    "Hg": "#b8b8d0",
    "Tl": "#a6544d",
    "Pb": "#575961",
    "Bi": "#9e4fb5",
    "Th": "#00baff",
    "U": "#008fff",
}


def _element(symbol: str) -> Any:
    """pymatgen's element record for a bare element symbol.

    Imported at call time rather than at module import: pymatgen loads a large
    data table on first use, and `pytex.core` must stay cheap to import. It is
    a required dependency, so the import cannot fail for a reason worth
    recovering from.
    """

    from pymatgen.core.periodic_table import Element

    return Element(symbol)


def normalize_species_symbol(species: str) -> str:
    match = re.match(r"([A-Z][a-z]?)", species.strip())
    if match is None:
        raise ValueError(f"Could not parse element symbol from species {species!r}.")
    return match.group(1)


def atomic_number(species: str) -> int:
    """Atomic number ``Z`` for a species label, from pymatgen's periodic table.

    ``species`` may carry a charge or a site suffix (``"As3-"``); it is reduced
    to its element symbol first. pymatgen is a required dependency, so there is
    one answer for every element on every machine -- the arrangement that
    replaced a hand-maintained partial table which held Ga and Ge but not As,
    and so raised on gallium arsenide wherever the extra was not installed.

    Raises
    ------
    ValueError
        If the symbol is not an element.
    """

    symbol = normalize_species_symbol(species)
    try:
        return int(_element(symbol).Z)
    except Exception as exc:
        raise ValueError(f"No atomic number available for species {species!r}.") from exc


#: Covalent radius used for elements pymatgen tabulates none for.
#:
#: Its table stops at curium, so this covers the trans-curium elements only. A
#: gap in the reference data, not a missing dependency: a mid-range value chosen
#: so bond detection stays plausible rather than collapsing to zero-length or
#: bonding everything to everything.
_UNTABULATED_COVALENT_RADIUS_ANGSTROM = 1.15


@lru_cache(maxsize=1)
def _covalent_radii() -> Mapping[str, float]:
    """pymatgen's Cordero covalent-radius table, keyed by element symbol.

    Read through both module paths because pymatgen moved this table from
    ``pymatgen.analysis`` to ``pymatgen.core`` and left a deprecation stub
    behind; importing the stub emits a ``DeprecationWarning`` that the test
    suite turns into an error. This is a version shim for one relocated import,
    not an optional-dependency guard -- pymatgen itself is required.

    These radii are *not* interchangeable with `Element.atomic_radius`: the
    covalent radius is what bond detection sums, and the atomic radius is a
    display size.
    """

    try:
        from pymatgen.core.molecule_structure_comparator import CovalentRadius
    except ImportError:  # pragma: no cover - pymatgen older than the move
        from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
    return {str(symbol): float(radius) for symbol, radius in CovalentRadius.radius.items()}


def covalent_radius_angstrom(species: str) -> float:
    """Covalent radius in angstrom, from pymatgen's table, for bond detection.

    This is the length bond detection sums; `display_radius_angstrom` chooses
    what an atom is *drawn* as, which is a separate question. Falls back to
    `_UNTABULATED_COVALENT_RADIUS_ANGSTROM` for the trans-curium elements the
    table does not reach.
    """

    symbol = normalize_species_symbol(species)
    return _covalent_radii().get(symbol, float(_UNTABULATED_COVALENT_RADIUS_ANGSTROM))


# Distinct hues for species outside the CPK table. A single shared grey fallback
# made two unlisted species indistinguishable in a figure, which is the one thing
# element colouring exists to prevent.
_UNLISTED_SPECIES_COLORS = (
    "#8c564b",
    "#17becf",
    "#bcbd22",
    "#e377c2",
    "#7f7f7f",
    "#9467bd",
)


def cpk_color(species: str) -> str:
    """Jmol/CPK display colour for an element, as used by VESTA-class viewers."""

    symbol = normalize_species_symbol(species)
    color = _FALLBACK_CPK_COLORS.get(symbol)
    if color is not None:
        return color
    index = sum(ord(character) for character in symbol) % len(_UNLISTED_SPECIES_COLORS)
    return _UNLISTED_SPECIES_COLORS[index]


def _element_data_radius(species: str, key: str) -> float | None:
    """Read a radius entry from pymatgen's element data table, if it has one."""

    symbol = normalize_species_symbol(species)
    try:
        value = _element(symbol).data.get(key)
    except Exception:
        return None
    if value is None:
        return None
    return float(value)


def atomic_radius_angstrom(species: str) -> float:
    """Empirical atomic radius in angstrom (Slater), for space-filling display.

    Falls back to 1.25x the covalent radius when pymatgen (or the element's
    tabulated value) is unavailable, which keeps relative atom sizes sensible.
    """

    value = _element_data_radius(species, "Atomic radius")
    if value is not None:
        return value
    return 1.25 * covalent_radius_angstrom(species)


def van_der_waals_radius_angstrom(species: str) -> float:
    """Van der Waals radius in angstrom, for molecular space-filling display.

    Falls back to 1.8x the covalent radius when no tabulated value exists.
    """

    value = _element_data_radius(species, "Van der waals radius")
    if value is not None:
        return value
    return 1.8 * covalent_radius_angstrom(species)


_DISPLAY_RADIUS_KINDS = ("covalent", "atomic", "van_der_waals")


def display_radius_angstrom(species: str, *, kind: str = "covalent") -> float:
    """Radius used for *display* sizing of atoms (never for bond chemistry).

    ``kind`` selects the radius system: ``"covalent"`` (ball-and-stick),
    ``"atomic"`` (space-filling of crystals), or ``"van_der_waals"``
    (space-filling of molecules). Bond detection always uses
    `covalent_radius_angstrom` regardless of the display choice.
    """

    if kind == "covalent":
        return covalent_radius_angstrom(species)
    if kind == "atomic":
        return atomic_radius_angstrom(species)
    if kind == "van_der_waals":
        return van_der_waals_radius_angstrom(species)
    raise ValueError(f"display radius kind must be one of {_DISPLAY_RADIUS_KINDS!r}.")
