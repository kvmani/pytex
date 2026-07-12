"""SO(3) radial kernel functions for ODF estimation.

The de la Vallee Poussin kernel is the MTEX workhorse kernel:
``psi(omega) = C(kappa) * cos^(2*kappa)(omega / 2)`` with the normalization
``C = B(3/2, 1/2) / B(3/2, kappa + 1/2)`` chosen so the kernel integrates to
one over SO(3) with the normalized Haar measure
``d mu = (2 / pi) * sin^2(omega / 2) d omega``.

Chebyshev (harmonic character) coefficients are computed by deterministic
quadrature of ``A_l = integral psi(omega) * chi_l(omega) d mu`` with the
SO(3) characters ``chi_l(omega) = sin((2l + 1) omega / 2) / sin(omega / 2)``,
so ``A_0 = 1`` exactly for a normalized kernel.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike

_QUADRATURE_SAMPLES = 20001


def _beta(a: float, b: float) -> float:
    # Log-gamma form: the direct gamma product overflows for the large kappa
    # values produced by sharp kernels (halfwidth below ~10 degrees).
    return math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))


def _kappa_from_halfwidth(halfwidth_rad: float) -> float:
    cosine = math.cos(halfwidth_rad / 2.0)
    if cosine <= 0.0:
        raise ValueError("Kernel halfwidth must be smaller than 180 degrees.")
    return math.log(0.5) / (2.0 * math.log(cosine))


@dataclass(frozen=True)
class DeLaValleePoussinKernel:
    halfwidth_deg: float

    def __post_init__(self) -> None:
        halfwidth = float(self.halfwidth_deg)
        if not 0.0 < halfwidth < 180.0 or not np.isfinite(halfwidth):
            raise ValueError(
                "DeLaValleePoussinKernel.halfwidth_deg must lie in (0, 180) degrees."
            )
        object.__setattr__(self, "halfwidth_deg", halfwidth)

    @property
    def halfwidth_rad(self) -> float:
        return math.radians(self.halfwidth_deg)

    @cached_property
    def kappa(self) -> float:
        return _kappa_from_halfwidth(self.halfwidth_rad)

    @cached_property
    def normalization(self) -> float:
        return _beta(1.5, 0.5) / _beta(1.5, self.kappa + 0.5)

    def evaluate(self, omega_rad: ArrayLike) -> np.ndarray:
        omega = np.asarray(omega_rad, dtype=np.float64)
        if np.any(np.abs(omega) > np.pi + 1e-12):
            raise ValueError("Kernel angles must lie in [-pi, pi] radians.")
        values = self.normalization * np.cos(omega / 2.0) ** (2.0 * self.kappa)
        values = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        values.setflags(write=False)
        return values

    def evaluate_deg(self, omega_deg: ArrayLike) -> np.ndarray:
        return self.evaluate(np.deg2rad(np.asarray(omega_deg, dtype=np.float64)))

    def chebyshev_coefficients(self, bandwidth: int) -> np.ndarray:
        """Return the character coefficients A_0 .. A_bandwidth by quadrature."""

        if bandwidth < 0:
            raise ValueError("bandwidth must be non-negative.")
        omega = np.linspace(0.0, np.pi, _QUADRATURE_SAMPLES)
        half = omega / 2.0
        psi = self.normalization * np.cos(half) ** (2.0 * self.kappa)
        measure = (2.0 / np.pi) * np.sin(half) ** 2
        orders = np.arange(bandwidth + 1, dtype=np.float64)
        # chi_l(omega) * sin(omega/2) = sin((2l+1) omega/2); folding the
        # character's sin(omega/2) denominator into the measure keeps the
        # integrand finite at omega = 0.
        numerators = np.sin((2.0 * orders[:, None] + 1.0) * half[None, :])
        integrand = psi[None, :] * (2.0 / np.pi) * np.sin(half)[None, :] * numerators
        del measure
        # Manual trapezoid rule on the uniform grid: portable across the
        # numpy 1.x (trapz) / 2.x (trapezoid) rename.
        spacing = float(omega[1] - omega[0])
        coefficients = spacing * (
            0.5 * integrand[:, 0]
            + integrand[:, 1:-1].sum(axis=1)
            + 0.5 * integrand[:, -1]
        )
        coefficients = np.ascontiguousarray(coefficients)
        coefficients.setflags(write=False)
        return coefficients

    def bandwidth(self, *, threshold: float = 1e-3, max_bandwidth: int = 512) -> int:
        """Smallest order whose coefficient magnitude drops below `threshold`."""

        if not 0.0 < threshold < 1.0:
            raise ValueError("threshold must lie in (0, 1).")
        coefficients = self.chebyshev_coefficients(max_bandwidth)
        below = np.nonzero(np.abs(coefficients) < threshold)[0]
        if below.size == 0:
            return max_bandwidth
        return int(below[0])


__all__ = ["DeLaValleePoussinKernel"]
