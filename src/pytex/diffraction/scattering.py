"""Tabulated X-ray atomic form factors for structure-factor calculations.

Implements the De Graef & McHenry parametrization used by pymatgen's XRD
calculator::

    f(s) = Z - 41.78214 * s^2 * sum_i a_i * exp(-b_i * s^2),   s = sin(theta)/lambda = |g| / 2

with per-element fitted coefficients ``(a_i, b_i)`` stored in
``_data/xray_scattering_factors.json``. The table is generated (with
provenance) by ``scripts/generate_xray_scattering_factor_table.py`` from
pymatgen's bundled ``ATOMIC_SCATTERING_PARAMS``, so PyTex intensities computed
with ``intensity_model="xray_tabulated"`` are directly comparable against the
pinned pymatgen external baselines.

At ``s = 0`` the form factor equals the atomic number ``Z`` (all electrons
scatter in phase), recovering the historical ``f = Z`` proxy as the
zero-angle limit; elements missing from the table fall back to ``Z`` so
exotic species keep working.
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from typing import Any

import numpy as np

from pytex.core._chemistry import atomic_number

_DATA_PACKAGE = "pytex.diffraction._data"
_TABLE_FILENAME = "xray_scattering_factors.json"
_DE_GRAEF_MCHENRY_SCALE = 41.78214


@lru_cache(maxsize=1)
def _scattering_table() -> dict[str, Any]:
    payload = (
        resources.files(_DATA_PACKAGE).joinpath(_TABLE_FILENAME).read_text(encoding="utf-8")
    )
    table = json.loads(payload)
    coefficients: dict[str, Any] = table["coefficients"]
    return coefficients


def tabulated_species() -> tuple[str, ...]:
    """Element symbols with fitted form-factor coefficients."""

    return tuple(sorted(_scattering_table()))


def xray_form_factors(species: str, s_values: np.ndarray) -> np.ndarray:
    """X-ray atomic form factor ``f(s)`` for one element at ``s = sin(theta)/lambda``.

    ``s_values`` is in 1/angstrom; returns an array of the same shape.
    Elements without tabulated coefficients fall back to the constant ``Z``
    (the zero-angle limit).
    """

    s = np.asarray(s_values, dtype=np.float64)
    z = float(atomic_number(species))
    entry = _scattering_table().get(species)
    if entry is None:
        return np.full_like(s, z)
    coefficients = np.asarray(entry, dtype=np.float64)  # (4, 2) rows of (a, b)
    s_squared = np.atleast_1d(s * s)
    summation = np.sum(
        coefficients[:, 0][:, None] * np.exp(-coefficients[:, 1][:, None] * s_squared[None, :]),
        axis=0,
    )
    result: np.ndarray = (z - _DE_GRAEF_MCHENRY_SCALE * s_squared * summation).reshape(s.shape)
    return result


def xray_form_factor_matrix(species: tuple[str, ...], s_values: np.ndarray) -> np.ndarray:
    """Form factors for several species at once: shape ``(len(s), len(species))``."""

    s = np.asarray(s_values, dtype=np.float64).reshape(-1)
    columns = [xray_form_factors(symbol, s) for symbol in species]
    return np.stack(columns, axis=1)


__all__ = [
    "tabulated_species",
    "xray_form_factor_matrix",
    "xray_form_factors",
]
