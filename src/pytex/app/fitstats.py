# ruff: noqa: RUF002
"""Standard errors of an unweighted straight-line fit, for the figures that state them.

Several workbench analyses end in an unweighted least-squares line whose
intercept or slope *is* the answer — the Williamson–Hall size and strain, the
CBED thickness and extinction distance. Their operations report the line and
not its uncertainty. The figures state the uncertainty the scatter about the
line implies, computed here once by the textbook formulas so every figure
quotes the same statistic:

``s² = Σr²/(n − 2)``, ``σ(slope)² = s²/Sxx``,
``σ(intercept)² = s²(1/n + x̄²/Sxx)`` and ``cov = −x̄ s²/Sxx``.

With two points the line passes through both, there is no scatter to estimate
from, and every uncertainty is ``nan`` — which is reported as such rather than
as zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["LineFit", "fit_line"]


@dataclass(frozen=True)
class LineFit:
    """An unweighted least-squares line and the standard errors of its coefficients."""

    slope: float
    intercept: float
    sigma_slope: float
    sigma_intercept: float
    covariance: float
    count: int

    @property
    def has_uncertainty(self) -> bool:
        """Whether there were more points than coefficients to estimate the scatter from."""

        return bool(np.isfinite(self.sigma_slope))

    def band(self, x: np.ndarray) -> np.ndarray:
        """One-standard-error band of the line's value at ``x``."""

        values = np.asarray(x, dtype=float)
        variance = (
            self.sigma_intercept**2
            + values**2 * self.sigma_slope**2
            + 2.0 * values * self.covariance
        )
        band: np.ndarray = np.sqrt(np.maximum(variance, 0.0))
        return band


def fit_line(x: np.ndarray, y: np.ndarray) -> LineFit:
    """Fit ``y = slope·x + intercept`` without weights and return its standard errors."""

    abscissa = np.asarray(x, dtype=float)
    ordinate = np.asarray(y, dtype=float)
    count = int(abscissa.size)
    slope, intercept = (float(value) for value in np.polyfit(abscissa, ordinate, 1))
    if count <= 2:
        nan = float("nan")
        return LineFit(slope, intercept, nan, nan, nan, count)
    residual = ordinate - (slope * abscissa + intercept)
    variance = float(np.sum(residual**2) / (count - 2))
    mean = float(np.mean(abscissa))
    sxx = float(np.sum((abscissa - mean) ** 2))
    return LineFit(
        slope=slope,
        intercept=intercept,
        sigma_slope=float(np.sqrt(variance / sxx)),
        sigma_intercept=float(np.sqrt(variance * (1.0 / count + mean**2 / sxx))),
        covariance=-mean * variance / sxx,
        count=count,
    )
