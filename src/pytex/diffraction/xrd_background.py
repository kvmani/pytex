"""Background estimation for measured powder-diffraction profiles.

A raw diffractogram is a sum of Bragg intensity and a slowly varying background
built from air scatter, sample fluorescence, incoherent scattering, the sample
holder and detector electronic noise. Every quantitative use of the pattern --
integrated intensities, peak widths, whole-profile refinement -- depends on
separating the two, and the separation is a *modelling choice* rather than a
measurement, so this module makes the choice explicit, records it, and explains
it.

Two estimators are offered because they answer different questions:

``"snip"``
    Statistics-sensitive Non-linear Iterative Peak-clipping (Ryan et al. 1988).
    Non-parametric: it assumes only that the background varies more slowly than
    the peaks, so it copes with curved, structured backgrounds -- amorphous
    humps, fluorescence steps -- that no low-order polynomial follows. Use it on
    unfamiliar data, and to *see* what the background is doing.

``"chebyshev"``
    A Chebyshev polynomial of the first kind fitted by iteratively reweighted
    least squares with asymmetric peak clipping. Parametric: it yields a handful
    of coefficients that a refinement can carry and refine, which is why Rietveld
    programs use this family. Use it when the background feeds a refinement.

Neither estimator knows where the peaks are. That is deliberate: a background
that was told where the peaks are cannot be used to *find* them, and the failure
would be silent.

References
----------
Ryan, C. G. et al., *Nucl. Instrum. Methods Phys. Res. B* **34** (1988) 396-402,
doi:10.1016/0168-583X(88)90063-8 -- the SNIP algorithm.

Morhac, M. et al., *Nucl. Instrum. Methods Phys. Res. A* **401** (1997) 113-132,
doi:10.1016/S0168-9002(97)01023-1 -- the decreasing-window refinement of SNIP.

Young, R. A. (ed.), *The Rietveld Method*, IUCr/OUP (1993), Ch. 1 -- background
treatment in whole-profile refinement.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

POWDER_BACKGROUND_SCHEMA = "pytex.diffraction.powder_background_estimate"

BACKGROUND_METHODS = ("snip", "chebyshev")

BackgroundMethod = Literal["snip", "chebyshev"]

_CITATION_SNIP = (
    "Ryan et al., Nucl. Instrum. Methods B 34 (1988) 396, doi:10.1016/0168-583X(88)90063-8."
)
_CITATION_MORHAC = (
    "Morhac et al., Nucl. Instrum. Methods A 401 (1997) 113, doi:10.1016/S0168-9002(97)01023-1."
)
_CITATION_RIETVELD_BOOK = "Young (ed.), The Rietveld Method, IUCr/OUP (1993), Ch. 1."


@dataclass(frozen=True, slots=True)
class PowderBackground:
    """An estimated background curve and the choices that produced it.

    Purpose
    -------
    Carry the background *and* its provenance together, so a subtracted pattern
    can always answer "what was removed, by what rule, with what settings".

    Attributes
    ----------
    two_theta_deg : np.ndarray
        The measured angular support the background was estimated on.
    background : np.ndarray
        The estimated background at each angle, in the measured intensity unit.
    observed_intensity : np.ndarray
        The raw intensity the estimate was made from, retained so the estimate
        is self-contained and re-checkable.
    method : str
        ``"snip"`` or ``"chebyshev"``.
    parameters : Mapping[str, float]
        The settings the estimator ran with, by name.
    source_name : str
        Name of the measured pattern the background belongs to.
    """

    two_theta_deg: np.ndarray
    background: np.ndarray
    observed_intensity: np.ndarray
    method: BackgroundMethod
    parameters: Mapping[str, float]
    source_name: str

    def __post_init__(self) -> None:
        axis = as_float_array(self.two_theta_deg, shape=(None,))
        background = as_float_array(self.background, shape=(None,))
        observed = as_float_array(self.observed_intensity, shape=(None,))
        if axis.size < 2:
            raise ValueError("A background estimate needs at least two angular points.")
        if axis.shape != background.shape or axis.shape != observed.shape:
            raise ValueError("Background arrays must align with the angular support.")
        if np.any(~np.isfinite(axis)) or np.any(np.diff(axis) <= 0.0):
            raise ValueError("Background angles must be finite and strictly increasing.")
        if np.any(~np.isfinite(background)) or np.any(background < 0.0):
            raise ValueError("Estimated background must be finite and non-negative.")
        if np.any(~np.isfinite(observed)) or np.any(observed < 0.0):
            raise ValueError("Observed intensity must be finite and non-negative.")
        if self.method not in BACKGROUND_METHODS:
            raise ValueError(f"Background method must be one of {BACKGROUND_METHODS}.")
        if not self.source_name.strip():
            raise ValueError("PowderBackground.source_name must be non-empty.")
        object.__setattr__(self, "two_theta_deg", axis)
        object.__setattr__(self, "background", background)
        object.__setattr__(self, "observed_intensity", observed)
        object.__setattr__(self, "parameters", MappingProxyType(dict(self.parameters)))

    @property
    def point_count(self) -> int:
        """Return the number of angular points the background was estimated on."""

        return int(self.two_theta_deg.size)

    @property
    def background_fraction(self) -> float:
        """Return the share of the total measured signal assigned to background.

        A value near one says the pattern is mostly background, which is the
        usual signature of a weak or badly aligned measurement; a value near
        zero on a laboratory pattern usually means the estimator has clipped
        into the peaks.
        """

        total = float(np.sum(self.observed_intensity))
        if total <= 0.0:
            return 0.0
        return float(np.sum(self.background) / total)

    def subtracted_intensity(self) -> np.ndarray:
        """Return the observed intensity less the background, clipped at zero.

        The clip is not cosmetic. A background estimate that crosses above the
        data in a low-count region would otherwise produce negative intensity,
        which every downstream invariant in PyTex rejects; the clip states that
        the model, not the measurement, is what went slightly wrong there.
        """

        subtracted = np.clip(self.observed_intensity - self.background, 0.0, None)
        subtracted = np.ascontiguousarray(subtracted, dtype=np.float64)
        subtracted.setflags(write=False)
        return subtracted

    def subtract(self, measured: MeasuredPowderPattern) -> MeasuredPowderPattern:
        """Return ``measured`` with this background removed.

        The result keeps the source radiation, unit and uncertainties, and
        records the subtraction in its metadata so a background-subtracted
        pattern can never be mistaken for a raw one.

        Parameters
        ----------
        measured
            The pattern the background was estimated from. Its angular support
            must match the estimate exactly, because a background interpolated
            onto a different grid is a second modelling step that this method
            deliberately refuses to take silently.
        """

        if measured.two_theta_deg.shape != self.two_theta_deg.shape or not np.allclose(
            measured.two_theta_deg, self.two_theta_deg, rtol=0.0, atol=1e-12
        ):
            raise ValueError(
                "A background can only be subtracted from the pattern it was estimated on."
            )
        metadata = dict(measured.metadata)
        metadata["background_subtracted"] = self.method
        metadata["background_fraction"] = f"{self.background_fraction:.6g}"
        return MeasuredPowderPattern(
            name=f"{measured.name} (background subtracted)",
            two_theta_deg=measured.two_theta_deg,
            intensity=self.subtracted_intensity(),
            standard_uncertainty=measured.standard_uncertainty,
            intensity_unit=measured.intensity_unit,
            radiation=measured.radiation,
            synthetic=measured.synthetic,
            metadata=metadata,
            provenance=measured.provenance,
        )

    def describe(self) -> str:
        """Return the method, its settings, what it removed, and its limits."""

        settings = ", ".join(
            f"{name} = {value:g}" for name, value in sorted(self.parameters.items())
        )
        if self.method == "snip":
            order = (
                f"widest window inward ({_CITATION_MORHAC})"
                if self.parameters.get("decreasing_window", 0.0)
                else "narrowest window outward, the original order"
            )
            basis = (
                "Statistics-sensitive Non-linear Iterative Peak-clipping in the "
                "log-log-square-root "
                f"domain ({_CITATION_SNIP}), clipped from the {order}, which assumes only that the "
                "background varies more slowly with angle than the Bragg peaks do."
            )
        else:
            basis = (
                "A Chebyshev polynomial of the first kind fitted by iteratively reweighted least "
                "squares with asymmetric clipping, the background family used in whole-profile "
                f"refinement ({_CITATION_RIETVELD_BOOK})."
            )
        return (
            f"Estimated the background of '{self.source_name}' over {self.point_count} points from "
            f"{self.two_theta_deg[0]:.4f} to {self.two_theta_deg[-1]:.4f} degrees 2*theta by the "
            f"{self.method} method with {settings}. Basis: {basis} "
            f"It accounts for {100.0 * self.background_fraction:.2f}% of the total measured "
            "signal, "
            f"with a mean level of {float(np.mean(self.background)):.6g} and a range of "
            f"{float(np.min(self.background)):.6g} to {float(np.max(self.background)):.6g} in the "
            "measured intensity unit. The estimate is a modelling choice, not a measurement: it "
            "does not identify peaks, and a background that clips into weak reflections removes "
            "intensity that a refinement will then be unable to account for."
        )


def _lls(values: np.ndarray) -> np.ndarray:
    """Log-log-square-root transform: compresses the dynamic range of counts.

    SNIP clips with a *linear* mean, so it needs a domain in which a Poisson
    peak looks the same at 10 counts and at 10,000. This is that domain.
    """

    return np.log(np.log(np.sqrt(values + 1.0) + 1.0) + 1.0)


def _inverse_lls(values: np.ndarray) -> np.ndarray:
    inverted = np.square(np.exp(np.exp(values) - 1.0) - 1.0) - 1.0
    return np.asarray(inverted, dtype=np.float64)


def _snip_background(
    intensity: np.ndarray, *, half_window: int, decreasing_window: bool
) -> np.ndarray:
    working = _lls(intensity)
    windows = range(half_window, 0, -1) if decreasing_window else range(1, half_window + 1)
    for width in windows:
        shifted_low = np.roll(working, width)
        shifted_high = np.roll(working, -width)
        shifted_low[:width] = working[:width]
        shifted_high[-width:] = working[-width:]
        clipped = 0.5 * (shifted_low + shifted_high)
        working = np.minimum(working, clipped)
    return np.asarray(np.clip(_inverse_lls(working), 0.0, None), dtype=np.float64)


def _chebyshev_background(
    axis: np.ndarray,
    intensity: np.ndarray,
    *,
    degree: int,
    iterations: int,
    clip_sigma: float,
) -> np.ndarray:
    # Chebyshev polynomials are orthogonal on [-1, 1], so the angular support is
    # mapped there; without the mapping the fit is numerically the same as a
    # power-basis fit and loses the conditioning that motivates the basis.
    span = float(axis[-1] - axis[0])
    reduced = 2.0 * (axis - axis[0]) / span - 1.0
    weights = np.ones_like(intensity)
    background = np.full_like(intensity, float(np.median(intensity)))
    for _ in range(iterations):
        coefficients = np.polynomial.chebyshev.chebfit(reduced, intensity, degree, w=weights)
        background = np.polynomial.chebyshev.chebval(reduced, coefficients)
        residual = intensity - background
        scale = float(np.std(residual))
        if scale <= 0.0:
            break
        # Asymmetric clipping: points far *above* the current curve are peaks and
        # are down-weighted; points below it are background scatter and are kept.
        # A symmetric rejection would pull the curve up into the peak feet.
        weights = np.where(residual > clip_sigma * scale, 0.0, 1.0)
        if not np.any(weights > 0.0):
            break
    return np.asarray(np.clip(background, 0.0, None), dtype=np.float64)


def estimate_background(
    measured: MeasuredPowderPattern,
    *,
    method: BackgroundMethod = "snip",
    half_window_deg: float = 2.0,
    decreasing_window: bool = True,
    degree: int = 6,
    iterations: int = 12,
    clip_sigma: float = 1.0,
) -> PowderBackground:
    """Estimate the slowly varying background of a measured powder pattern.

    Purpose
    -------
    Separate Bragg intensity from the instrumental and specimen background, so
    integrated intensities, peak widths and whole-profile refinement operate on
    the diffracted signal rather than on the signal plus the room.

    Method
    ------
    ``"snip"`` clips the pattern iteratively against the mean of its own
    neighbours at a growing (or shrinking) separation, in the log-log-square-root
    domain that makes the clip insensitive to count level. Anything narrower than
    the window is removed as a peak; anything broader survives as background.

    ``"chebyshev"`` fits a Chebyshev polynomial by iteratively reweighted least
    squares, discarding points more than ``clip_sigma`` residual standard
    deviations *above* the current curve on each pass. The asymmetry is the
    point: peaks are one-sided excursions, so a symmetric rejection would drag
    the background up into the peak feet.

    Parameters
    ----------
    measured
        The raw profile. Background estimation on an already-subtracted pattern
        is meaningless, so pass the instrument's own output.
    method
        ``"snip"`` (default, non-parametric) or ``"chebyshev"`` (parametric).
    half_window_deg
        SNIP clipping half-window in degrees 2*theta. Set it comfortably wider
        than the broadest peak half-width and comfortably narrower than the
        curvature of the background; too small leaves peak feet in the
        background, too large flattens genuine background structure.
    decreasing_window
        Run SNIP from the widest window inward (Morhac et al. 1997), which
        preserves narrow background structure better than the original
        increasing-window order. Ignored by the Chebyshev method.
    degree
        Chebyshev polynomial degree. Ignored by SNIP.
    iterations
        Number of reweighting passes for the Chebyshev fit. Ignored by SNIP.
    clip_sigma
        Chebyshev rejection threshold in residual standard deviations. Ignored
        by SNIP.

    Returns
    -------
    PowderBackground
        The background curve, the settings that produced it, and a
        :meth:`~PowderBackground.describe` summary.

    Raises
    ------
    ValueError
        If the method is unknown, or a setting is outside its valid range, or
        the SNIP window is too wide for the measured angular support.

    Examples
    --------
    >>> import numpy as np
    >>> from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
    >>> from pytex.diffraction.xrd_background import estimate_background
    >>> angles = np.linspace(20.0, 80.0, 601)
    >>> peak = 100.0 * np.exp(-0.5 * ((angles - 45.0) / 0.2) ** 2)
    >>> pattern = MeasuredPowderPattern(
    ...     name="demo", two_theta_deg=angles, intensity=peak + 10.0, synthetic=True
    ... )
    >>> background = estimate_background(pattern, half_window_deg=2.0)
    >>> bool(abs(float(np.median(background.background)) - 10.0) < 0.5)
    True

    See Also
    --------
    pytex.diffraction.rietveld.refine_rietveld : Refines a Chebyshev background
        jointly with the structural and profile parameters, which is preferable
        to subtracting a fixed background first.
    """

    if method not in BACKGROUND_METHODS:
        raise ValueError(f"method must be one of {BACKGROUND_METHODS}.")
    axis = measured.two_theta_deg
    intensity = measured.intensity
    parameters: dict[str, float]
    if method == "snip":
        if half_window_deg <= 0.0:
            raise ValueError("half_window_deg must be strictly positive.")
        step = float(np.median(np.diff(axis)))
        half_window = round(half_window_deg / step)
        if half_window < 1:
            raise ValueError(
                "half_window_deg is narrower than one measured step; SNIP would do nothing."
            )
        if half_window >= axis.size // 2:
            raise ValueError(
                "half_window_deg spans half the measured range or more; there is no background "
                "structure left for SNIP to separate from the peaks."
            )
        background = _snip_background(
            intensity, half_window=half_window, decreasing_window=decreasing_window
        )
        parameters = {
            "half_window_deg": float(half_window_deg),
            "half_window_points": float(half_window),
            "decreasing_window": float(decreasing_window),
        }
    else:
        if degree < 0:
            raise ValueError("degree must be non-negative.")
        if degree >= axis.size:
            raise ValueError("degree must be smaller than the number of measured points.")
        if iterations < 1:
            raise ValueError("iterations must be at least one.")
        if clip_sigma <= 0.0:
            raise ValueError("clip_sigma must be strictly positive.")
        background = _chebyshev_background(
            axis,
            intensity,
            degree=degree,
            iterations=iterations,
            clip_sigma=clip_sigma,
        )
        parameters = {
            "degree": float(degree),
            "iterations": float(iterations),
            "clip_sigma": float(clip_sigma),
        }
    return PowderBackground(
        two_theta_deg=axis,
        background=background,
        observed_intensity=intensity,
        method=method,
        parameters=parameters,
        source_name=measured.name,
    )


__all__ = [
    "BACKGROUND_METHODS",
    "POWDER_BACKGROUND_SCHEMA",
    "BackgroundMethod",
    "PowderBackground",
    "estimate_background",
]
