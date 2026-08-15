"""Finite-size shape factors for kinematic electron diffraction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

FINITE_THICKNESS_SHAPE_FACTOR_SCHEMA = "pytex.diffraction.finite_thickness_shape_factor"


@dataclass(frozen=True, slots=True)
class FiniteThicknessShapeFactor:
    """Normalized kinematic shape factor for a plane-parallel thin foil.

    Purpose
    -------
    Convert excitation error into the finite-thickness relrod modulation for
    selected-area electron diffraction (SAED). For a uniform slab of thickness
    ``t`` normal to the incident beam, the normalized amplitude is
    ``sin(pi t s_g) / (pi t s_g)`` and the normalized intensity is its square.

    When to use
    -----------
    Use this model when foil thickness is known and the single-scattering,
    plane-parallel approximation is adequate. It replaces an empirical
    Lorentzian relrod-width proxy; it does not model dynamical scattering,
    absorption, bending, mosaic spread, or a thickness distribution.

    Parameters
    ----------
    thickness_angstrom : float
        Positive specimen thickness normal to the beam, in angstrom.

    Returns
    -------
    FiniteThicknessShapeFactor
        Immutable model exposing amplitude and intensity factors, its first
        relrod zero, and a convention-explicit :meth:`describe` summary.
    """

    thickness_angstrom: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.thickness_angstrom) or self.thickness_angstrom <= 0.0:
            raise ValueError("thickness_angstrom must be finite and strictly positive.")
        object.__setattr__(self, "thickness_angstrom", float(self.thickness_angstrom))

    @property
    def first_zero_inv_angstrom(self) -> float:
        """Absolute excitation error of the first relrod zero, ``1 / t``."""

        return 1.0 / self.thickness_angstrom

    def amplitude_factor(
        self, excitation_error_inv_angstrom: float | np.ndarray
    ) -> float | np.ndarray:
        """Return ``sinc(t s_g)`` for scalar or array excitation error.

        NumPy's normalized ``sinc(x) = sin(pi x)/(pi x)`` supplies the
        analytic limiting value one at the exact Bragg condition.
        """

        excitation = np.asarray(excitation_error_inv_angstrom, dtype=np.float64)
        if np.any(~np.isfinite(excitation)):
            raise ValueError("excitation_error_inv_angstrom must be finite.")
        factor = np.sinc(self.thickness_angstrom * excitation)
        if excitation.ndim == 0:
            return float(factor)
        return np.asarray(factor, dtype=np.float64)

    def intensity_factor(
        self, excitation_error_inv_angstrom: float | np.ndarray
    ) -> float | np.ndarray:
        """Return normalized ``sinc^2(t s_g)`` intensity modulation."""

        amplitude = self.amplitude_factor(excitation_error_inv_angstrom)
        return amplitude * amplitude

    def describe(self) -> str:
        """Explain the physical convention, scale, and limits of this model."""

        return (
            "Plane-parallel finite-thickness SAED shape factor for a "
            f"{self.thickness_angstrom:g} angstrom foil normal to the beam: normalized "
            "amplitude sinc(t s_g) = sin(pi t s_g)/(pi t s_g), normalized intensity "
            f"sinc^2(t s_g), and first zero at |s_g| = "
            f"{self.first_zero_inv_angstrom:g} 1/angstrom. This single-scattering "
            "rectangular-slab model excludes dynamical scattering, absorption, bending, "
            "mosaic spread, and thickness variation."
        )
