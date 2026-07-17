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
from collections.abc import Callable
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




def _character_values(orders: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """SO(3) characters chi_l(omega) = sin((2l+1) omega/2) / sin(omega/2)."""

    half = omega / 2.0
    sin_half = np.sin(half)
    numerators = np.sin((2.0 * orders[:, None] + 1.0) * half[None, :])
    limits = np.broadcast_to((2.0 * orders[:, None] + 1.0), numerators.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = numerators / sin_half[None, :]
    return np.asarray(np.where(np.abs(sin_half)[None, :] > 1e-12, ratio, limits))


def _series_evaluate(coefficients: np.ndarray, omega: np.ndarray) -> np.ndarray:
    orders = np.arange(coefficients.shape[0], dtype=np.float64)
    characters = _character_values(orders, omega)
    return np.asarray(coefficients @ characters, dtype=np.float64)


def _truncated_coefficients(
    coefficient_of_order: Callable[[int], float],
    *,
    tail_threshold: float = 1e-10,
    max_order: int = 2048,
) -> np.ndarray:
    values = []
    for order in range(max_order + 1):
        value = float(coefficient_of_order(order))
        values.append(value)
        if order > 4 and abs(value) < tail_threshold:
            break
    return np.asarray(values, dtype=np.float64)


def _solve_halfwidth_parameter(
    halfwidth_rad: float,
    coefficients_for: Callable[[float], np.ndarray],
    *,
    low: float,
    high: float,
    broader_is_larger: bool,
) -> float:
    """Bisection on the defining property psi(halfwidth) = psi(0) / 2."""

    target = np.asarray([0.0, halfwidth_rad], dtype=np.float64)

    def ratio(parameter: float) -> float:
        coefficients = coefficients_for(parameter)
        values = _series_evaluate(coefficients, target)
        return float(values[1] / values[0])

    low_value, high_value = ratio(low), ratio(high)
    if not min(low_value, high_value) <= 0.5 <= max(low_value, high_value):
        raise ValueError("Kernel halfwidth is outside the solvable parameter range.")
    for _ in range(200):
        middle = 0.5 * (low + high)
        middle_value = ratio(middle)
        if (middle_value < 0.5) == (low_value < 0.5):
            low, low_value = middle, middle_value
        else:
            high, high_value = middle, middle_value
        if abs(high - low) < 1e-14 * max(1.0, abs(high)):
            break
    return 0.5 * (low + high)


@dataclass(frozen=True)
class GaussianSO3Kernel:
    """Gauss-Weierstrass kernel on SO(3), defined spectrally.

    Purpose: the "Gaussian" texture kernel — Chebyshev (character)
    coefficients ``A_l = (2l + 1) exp(-l (l + 1) epsilon)`` with the spread
    ``epsilon`` solved from the halfwidth via the defining property
    ``psi(halfwidth) = psi(0) / 2``. Being the heat-kernel spectrum, it is
    the smoothest kernel for a given halfwidth; ``A_0 = 1`` exactly, so the
    kernel is normalized over SO(3).

    Inputs: ``halfwidth_deg`` in (0, 180).

    Output surface mirrors ``DeLaValleePoussinKernel``: ``evaluate`` /
    ``evaluate_deg``, closed-form ``chebyshev_coefficients``, ``bandwidth``.
    """

    halfwidth_deg: float

    def __post_init__(self) -> None:
        halfwidth = float(self.halfwidth_deg)
        if not 0.0 < halfwidth < 180.0 or not np.isfinite(halfwidth):
            raise ValueError("GaussianSO3Kernel.halfwidth_deg must lie in (0, 180) degrees.")
        object.__setattr__(self, "halfwidth_deg", halfwidth)

    @property
    def halfwidth_rad(self) -> float:
        return math.radians(self.halfwidth_deg)

    @cached_property
    def epsilon(self) -> float:
        return _solve_halfwidth_parameter(
            self.halfwidth_rad,
            lambda parameter: _truncated_coefficients(
                lambda order: (2 * order + 1) * math.exp(-order * (order + 1) * parameter)
            ),
            low=1e-8,
            high=5.0,
            broader_is_larger=True,
        )

    @cached_property
    def _coefficients(self) -> np.ndarray:
        epsilon = self.epsilon
        values = _truncated_coefficients(
            lambda order: (2 * order + 1) * math.exp(-order * (order + 1) * epsilon)
        )
        values.setflags(write=False)
        return values

    def evaluate(self, omega_rad: ArrayLike) -> np.ndarray:
        omega = np.asarray(omega_rad, dtype=np.float64)
        if np.any(np.abs(omega) > np.pi + 1e-12):
            raise ValueError("Kernel angles must lie in [-pi, pi] radians.")
        values = _series_evaluate(self._coefficients, np.atleast_1d(omega)).reshape(omega.shape)
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def evaluate_deg(self, omega_deg: ArrayLike) -> np.ndarray:
        return self.evaluate(np.deg2rad(np.asarray(omega_deg, dtype=np.float64)))

    def chebyshev_coefficients(self, bandwidth: int) -> np.ndarray:
        if bandwidth < 0:
            raise ValueError("bandwidth must be non-negative.")
        orders = np.arange(bandwidth + 1, dtype=np.float64)
        coefficients = (2.0 * orders + 1.0) * np.exp(-orders * (orders + 1.0) * self.epsilon)
        coefficients = np.ascontiguousarray(coefficients)
        coefficients.setflags(write=False)
        return coefficients

    def bandwidth(self, *, threshold: float = 1e-3, max_bandwidth: int = 512) -> int:
        if not 0.0 < threshold < 1.0:
            raise ValueError("threshold must lie in (0, 1).")
        coefficients = self.chebyshev_coefficients(max_bandwidth)
        below = np.nonzero(np.abs(coefficients) < threshold)[0]
        if below.size == 0:
            return max_bandwidth
        return int(below[0])


@dataclass(frozen=True)
class AbelPoissonKernel:
    """Abel-Poisson kernel on SO(3), defined spectrally.

    Purpose: the classical Abel-Poisson texture kernel — Chebyshev
    coefficients ``A_l = (2l + 1) kappa^(2l)`` with ``kappa`` in (0, 1)
    solved from the halfwidth via ``psi(halfwidth) = psi(0) / 2``. Its
    geometric coefficient decay makes it broader-tailed than the Gaussian at
    equal halfwidth; ``A_0 = 1`` exactly (normalized over SO(3)).

    Inputs: ``halfwidth_deg`` in (0, 180). Surface mirrors
    ``DeLaValleePoussinKernel``.
    """

    halfwidth_deg: float

    def __post_init__(self) -> None:
        halfwidth = float(self.halfwidth_deg)
        if not 0.0 < halfwidth < 180.0 or not np.isfinite(halfwidth):
            raise ValueError("AbelPoissonKernel.halfwidth_deg must lie in (0, 180) degrees.")
        object.__setattr__(self, "halfwidth_deg", halfwidth)

    @property
    def halfwidth_rad(self) -> float:
        return math.radians(self.halfwidth_deg)

    @cached_property
    def kappa(self) -> float:
        return _solve_halfwidth_parameter(
            self.halfwidth_rad,
            lambda parameter: _truncated_coefficients(
                lambda order: (2 * order + 1) * parameter ** (2 * order)
            ),
            low=1e-6,
            high=1.0 - 1e-9,
            broader_is_larger=False,
        )

    @cached_property
    def _coefficients(self) -> np.ndarray:
        kappa = self.kappa
        values = _truncated_coefficients(
            lambda order: (2 * order + 1) * kappa ** (2 * order)
        )
        values.setflags(write=False)
        return values

    def evaluate(self, omega_rad: ArrayLike) -> np.ndarray:
        omega = np.asarray(omega_rad, dtype=np.float64)
        if np.any(np.abs(omega) > np.pi + 1e-12):
            raise ValueError("Kernel angles must lie in [-pi, pi] radians.")
        values = _series_evaluate(self._coefficients, np.atleast_1d(omega)).reshape(omega.shape)
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def evaluate_deg(self, omega_deg: ArrayLike) -> np.ndarray:
        return self.evaluate(np.deg2rad(np.asarray(omega_deg, dtype=np.float64)))

    def chebyshev_coefficients(self, bandwidth: int) -> np.ndarray:
        if bandwidth < 0:
            raise ValueError("bandwidth must be non-negative.")
        orders = np.arange(bandwidth + 1, dtype=np.float64)
        coefficients = (2.0 * orders + 1.0) * self.kappa ** (2.0 * orders)
        coefficients = np.ascontiguousarray(coefficients)
        coefficients.setflags(write=False)
        return coefficients

    def bandwidth(self, *, threshold: float = 1e-3, max_bandwidth: int = 512) -> int:
        if not 0.0 < threshold < 1.0:
            raise ValueError("threshold must lie in (0, 1).")
        coefficients = self.chebyshev_coefficients(max_bandwidth)
        below = np.nonzero(np.abs(coefficients) < threshold)[0]
        if below.size == 0:
            return max_bandwidth
        return int(below[0])


__all__ = ["AbelPoissonKernel", "DeLaValleePoussinKernel", "GaussianSO3Kernel"]
