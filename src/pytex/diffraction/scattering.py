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
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import Phase

_DATA_PACKAGE = "pytex.diffraction._data"
_TABLE_FILENAME = "xray_scattering_factors.json"
_DE_GRAEF_MCHENRY_SCALE = 41.78214

#: Electron rest energy in keV, for the relativistic correction ``gamma``.
_ELECTRON_REST_ENERGY_KEV = 510.99895


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


def electron_structure_factor_angstrom(
    phase: Phase,
    hkl: ArrayLike,
    *,
    beam_energy_kev: float = 200.0,
) -> np.ndarray:
    """Electron structure factors ``F_g`` in angstrom, relativistically corrected.

    What it does
        Sums the Mott-Bethe electron scattering factors over the unit cell,

        .. math::

           F_{g} = \\gamma \\sum_{j} o_{j}\\, f_{e}^{j}(s)\\,
                   e^{-B_{j}s^{2}}\\, e^{2\\pi i\\,\\mathbf{g}\\cdot\\mathbf{r}_{j}},
           \\qquad s = |\\mathbf{g}|/2,

        with the relativistic factor
        :math:`\\gamma = 1 + E/m_{0}c^{2}` applied because the incident electron
        is not slow: at 200 kV it is 1.39, and omitting it would make every
        extinction distance 39 percent too long.

    When to use it
        Whenever the *absolute* scale matters — extinction distances, dynamical
        rocking curves, thickness determination. For relative spot intensities
        within one zone, `pytex.diffraction.kinematic.electron_structure_factors`
        is the cheaper atomic-number proxy and is sufficient.

    Parameters
    ----------
    phase:
        Must carry a unit cell with sites; a bare lattice has no structure
        factor and raises rather than returning a fabricated one.
    hkl:
        ``(n, 3)`` integer Miller indices, or one triple.
    beam_energy_kev:
        Accelerating voltage, for the wavelength and for ``gamma``.

    Returns
    -------
    np.ndarray
        ``(n,)`` complex structure factors in angstrom.

    Raises
    ------
    ValueError
        If the phase has no unit-cell sites.
    """

    indices = np.atleast_2d(np.asarray(hkl, dtype=np.int64))
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (3,) or (n, 3).")
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(
            "Electron structure factors on an absolute scale need the atom positions: "
            f"phase '{phase.name}' carries no unit cell. Load the phase from a CIF, or "
            "use the relative-intensity path in pytex.diffraction.kinematic."
        )
    if not np.isfinite(beam_energy_kev) or beam_energy_kev <= 0.0:
        raise ValueError("beam_energy_kev must be finite and strictly positive.")

    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_cartesian = indices.astype(np.float64) @ reciprocal.T
    s_values = np.linalg.norm(g_cartesian, axis=1) / 2.0

    gamma = 1.0 + float(beam_energy_kev) / _ELECTRON_REST_ENERGY_KEV
    total = np.zeros(indices.shape[0], dtype=np.complex128)
    for site in phase.unit_cell.sites:
        factors = electron_scattering_factors(site.species, s_values)
        b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
        damping = np.exp(-b_iso * s_values * s_values)
        phase_argument = indices.astype(np.float64) @ site.fractional_coordinates
        total += (
            float(site.occupancy)
            * factors
            * damping
            * np.exp(2.0j * np.pi * phase_argument)
        )
    return np.asarray(gamma * total, dtype=np.complex128)


__all__ = [
    "electron_scattering_factors",
    "electron_structure_factor_angstrom",
    "tabulated_species",
    "xray_form_factor_matrix",
    "xray_form_factors",
]
