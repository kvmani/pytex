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


def electron_scattering_factors(species: str, s_values: np.ndarray) -> np.ndarray:
    """Electron atomic scattering factor ``f_e(s)`` in angstrom, by Mott-Bethe.

    What it does
        Converts the tabulated X-ray form factor into the electron scattering
        factor through the Mott-Bethe relation

        .. math::

           f_{e}(s) = \\frac{Z - f_{x}(s)}{8\\pi^{2}a_{0}\\,s^{2}}
                    = \\frac{Z - f_{x}(s)}{41.78214\\,s^{2}},
           \\qquad s = \\frac{\\sin\\theta}{\\lambda} = \\frac{|g|}{2}.

        Electrons are scattered by the *electrostatic potential*, so they see
        the nuclear charge screened by the electron cloud — which is exactly
        what :math:`Z - f_{x}` measures. The result is a length, in angstrom,
        not a dimensionless electron count.

    When to use it
        For any electron-diffraction quantity that needs an absolute scale
        rather than relative intensities: structure factors in angstrom,
        extinction distances, and therefore CBED thickness determination. The
        kinematic SAED path deliberately does **not** use it, because relative
        spot intensities within one zone do not need the absolute scale.

    Parameters
    ----------
    species:
        Element symbol.
    s_values:
        ``sin(theta)/lambda`` in 1/angstrom, any shape. The Mott-Bethe relation
        has a removable singularity at ``s = 0``; the zero-angle limit
        ``sum(a_i)`` is returned there.

    Returns
    -------
    np.ndarray
        ``f_e(s)`` in angstrom, the shape of ``s_values``.

    Notes
    -----
    Because the tabulated coefficients are of the form
    ``f_x = Z - 41.78214 s^2 sum_i a_i exp(-b_i s^2)``, this relation returns
    exactly ``sum_i a_i exp(-b_i s^2)``: the same coefficients *are* the
    electron scattering factor parameters, and the 41.78214 in the X-ray
    parametrization is precisely :math:`8\\pi^{2}a_{0}`. No additional
    constant is introduced here.

    The parametrization is fitted, and its accuracy at small ``s`` degrades for
    heavy elements. Absolute extinction distances computed from it should be
    treated as good to a few percent for light elements and as indicative for
    heavy ones — which is one reason CBED practice *measures* the extinction
    distance from the fringe positions rather than taking it from a table.

    See Also
    --------
    `pytex.diffraction.cbed.extinction_distance_angstrom` : the main consumer.
    """

    s = np.asarray(s_values, dtype=np.float64)
    entry = _scattering_table().get(species)
    if entry is None:
        # No fitted coefficients: the X-ray fallback is the constant Z, so the
        # Mott-Bethe numerator vanishes identically and there is nothing to
        # report rather than a value to invent.
        return np.zeros_like(s)
    coefficients = np.asarray(entry, dtype=np.float64)
    zero_angle_limit = float(np.sum(coefficients[:, 0]))

    # Evaluate Mott-Bethe away from the removable singularity, and substitute
    # the analytic limit at s = 0 rather than dividing by zero there.
    safe_s = np.where(np.abs(s) < 1e-12, 1.0, s)
    numerator = float(atomic_number(species)) - xray_form_factors(species, safe_s)
    factors = numerator / (_DE_GRAEF_MCHENRY_SCALE * safe_s * safe_s)
    return np.asarray(np.where(np.abs(s) < 1e-12, zero_angle_limit, factors), dtype=np.float64)


__all__ = [
    "electron_scattering_factors",
    "tabulated_species",
    "xray_form_factor_matrix",
    "xray_form_factors",
]
