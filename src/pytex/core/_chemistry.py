from __future__ import annotations

import re

_FALLBACK_ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Cl": 17,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Ti": 22,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "Zr": 40,
    "Mo": 42,
    "Ag": 47,
    "Sn": 50,
    "W": 74,
    "Au": 79,
    "Pb": 82,
}

_FALLBACK_COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ti": 1.60,
    "Fe": 1.24,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
}

_FALLBACK_CPK_COLORS = {
    "H": "#f8fafc",
    "C": "#334155",
    "N": "#2563eb",
    "O": "#dc2626",
    "Na": "#7c3aed",
    "Mg": "#22c55e",
    "Al": "#94a3b8",
    "Cl": "#22c55e",
    "Si": "#c084fc",
    "P": "#d97706",
    "S": "#facc15",
    "Ti": "#64748b",
    "Fe": "#b45309",
    "Ni": "#10b981",
    "Cu": "#b45309",
    "Zn": "#94a3b8",
}


def normalize_species_symbol(species: str) -> str:
    match = re.match(r"([A-Z][a-z]?)", species.strip())
    if match is None:
        raise ValueError(f"Could not parse element symbol from species {species!r}.")
    return match.group(1)


def atomic_number(species: str) -> int:
    symbol = normalize_species_symbol(species)
    try:
        from pymatgen.core.periodic_table import Element

        return int(Element(symbol).Z)
    except Exception as exc:
        if symbol not in _FALLBACK_ATOMIC_NUMBERS:
            raise ValueError(f"No atomic number available for species {species!r}.") from exc
        return _FALLBACK_ATOMIC_NUMBERS[symbol]


def covalent_radius_angstrom(species: str) -> float:
    symbol = normalize_species_symbol(species)
    try:
        from pymatgen.core.periodic_table import Element

        radius = Element(symbol).covalent_radius
        if radius is None:
            raise ValueError
        return float(radius)
    except Exception:
        return float(_FALLBACK_COVALENT_RADII.get(symbol, 1.15))


def cpk_color(species: str) -> str:
    return _FALLBACK_CPK_COLORS.get(normalize_species_symbol(species), "#64748b")


def _element_data_radius(species: str, key: str) -> float | None:
    """Read a radius entry from pymatgen's element data table, if available."""

    symbol = normalize_species_symbol(species)
    try:
        from pymatgen.core.periodic_table import Element

        value = Element(symbol).data.get(key)
        if value is None:
            return None
        return float(value)
    except Exception:
        return None


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
