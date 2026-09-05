"""Peak detection and single-peak profile fitting for measured powder patterns.

Everything downstream of a measured diffractogram -- indexing, precise lattice
parameters, size-strain analysis, residual stress -- begins with the same two
numbers per reflection: **where the peak is** and **how well that position is
known**. This module produces both, and treats the second as a first-class
result rather than an afterthought, because a lattice parameter quoted without
an uncertainty cannot be compared with anything.

Three ideas keep the module honest:

1. **Detection uses the instrument's own resolution function.** A peak is not
   "a local maximum above a threshold": it is a feature whose *width* matches
   what the diffractometer can produce at that angle. Filtering with a
   scale-matched Ricker (Mexican-hat) kernel, whose scale is set by the
   calibrated Caglioti curve at each angle, rejects noise spikes that are too
   narrow and background undulations that are too broad, without a single
   hand-tuned threshold in angle units.

2. **K-alpha2 is modelled, never stripped, on the fitting path.** Rachinger
   stripping (see :mod:`pytex.diffraction.xrd_corrections`) is a display and
   peak-picking convenience that assumes an exact 1:2 intensity ratio and an
   identical profile shape, and it amplifies counting noise multiplicatively
   down the scan. When the doublet is *fitted* instead, the alpha2 partner is
   constrained to the position Bragg's law puts it at, and no noise is created.
   :func:`fit_peaks` therefore reports the K-alpha1 position directly.

3. **The reported position uncertainty comes from the fit's own covariance.**
   It is the square root of the leading diagonal element of
   ``(J^T W J)^-1`` scaled by the reduced chi-squared of the window, which is
   the goodness-of-fit-scaled convention used throughout crystallography.

References
----------
Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Chs. 6 and 11 -- peak position measurement and the
precision it can reach.

Marple, S. L., *Digital Spectral Analysis*; and Du, P., Kibbe, W. A. & Lin, S. M.,
*Bioinformatics* **22** (2006) 2059-2065, doi:10.1093/bioinformatics/btl355 --
continuous-wavelet peak detection with a Ricker kernel, the scheme adapted here
to a known, angle-dependent instrumental width.

Thompson, P., Cox, D. E. & Hastings, J. B., *J. Appl. Crystallogr.* **20** (1987)
79-83, doi:10.1107/S0021889887087090 -- the pseudo-Voigt profile fitted here.

Toraya, H., *J. Appl. Crystallogr.* **19** (1986) 440-447,
doi:10.1107/S0021889886088982 -- the split pseudo-Voigt used for low-angle
axial-divergence asymmetry.

Bearden, J. A., *Rev. Mod. Phys.* **39** (1967) 78-124,
doi:10.1103/RevModPhys.39.78 -- the K-alpha1/K-alpha2 wavelengths that fix the
doublet separation.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
from scipy.optimize import least_squares

from pytex.core._arrays import as_float_array
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_background import estimate_background
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

PEAK_FIT_SCHEMA = "pytex.diffraction.powder_peak_fit"
PEAK_TABLE_SCHEMA = "pytex.diffraction.powder_peak_table"

PeakShape = Literal["pseudo_voigt", "split_pseudo_voigt"]

PEAK_SHAPES: tuple[PeakShape, ...] = ("pseudo_voigt", "split_pseudo_voigt")

#: Gaussian integrated-area coefficient: a unit-height Gaussian of unit FWHM
#: encloses ``sqrt(pi / (4 ln 2))``.
_GAUSSIAN_AREA = float(np.sqrt(np.pi / (4.0 * np.log(2.0))))

#: Lorentzian integrated-area coefficient: a unit-height Lorentzian of unit
#: FWHM encloses ``pi / 2``.
_LORENTZIAN_AREA = float(np.pi / 2.0)

_CITATION_CULLITY = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Chs. 6, 11."
)
_CITATION_TORAYA = (
    "Toraya, J. Appl. Crystallogr. 19 (1986) 440, doi:10.1107/S0021889886088982."
)
_CITATION_RICKER = (
    "Du, Kibbe & Lin, Bioinformatics 22 (2006) 2059, doi:10.1093/bioinformatics/btl355."
)


# ---------------------------------------------------------------------------
# Profile shapes
# ---------------------------------------------------------------------------


def pseudo_voigt_profile(
    two_theta_deg: Any,
    *,
    centre_deg: float,
    fwhm_deg: float,
    eta: float,
) -> np.ndarray:
    """Return the unit-height pseudo-Voigt profile on ``two_theta_deg``.

    Purpose
    -------
    Provide the single peak shape that laboratory powder-diffraction profile
    fitting has converged on, in the height-normalized form that makes a fitted
    amplitude read directly as a peak height.

    Method
    ------
    The pseudo-Voigt is the linear mixture
    ``eta * L(x) + (1 - eta) * G(x)`` of a Lorentzian and a Gaussian of the
    *same* full width at half maximum, so the mixture also has that FWHM. Both
    components are normalized to unit height at the centre:

    ``G(x) = exp(-4 ln2 (x - c)^2 / w^2)`` and
    ``L(x) = 1 / (1 + 4 (x - c)^2 / w^2)``.

    Parameters
    ----------
    two_theta_deg
        Angles at which to evaluate the profile, in degrees.
    centre_deg
        Peak centre in degrees ``2*theta``.
    fwhm_deg
        Full width at half maximum in degrees ``2*theta``. Must be positive.
    eta
        Lorentzian fraction in ``[0, 1]``. ``0`` is pure Gaussian, ``1`` pure
        Lorentzian.

    Returns
    -------
    np.ndarray
        Profile values with the same shape as ``two_theta_deg``, equal to one
        at ``centre_deg``.

    Raises
    ------
    ValueError
        If ``fwhm_deg`` is not positive or ``eta`` lies outside ``[0, 1]``.

    See Also
    --------
    split_pseudo_voigt_profile : the asymmetric generalization.
    pseudo_voigt_area : the integrated area of this profile.
    """

    if not np.isfinite(fwhm_deg) or fwhm_deg <= 0.0:
        raise ValueError("pseudo_voigt_profile requires a finite, positive fwhm_deg.")
    if not np.isfinite(eta) or not 0.0 <= eta <= 1.0:
        raise ValueError("pseudo_voigt_profile requires eta in [0, 1].")
    axis = np.asarray(two_theta_deg, dtype=np.float64)
    reduced = 2.0 * (axis - centre_deg) / fwhm_deg
    gaussian = np.exp(-np.log(2.0) * np.square(reduced))
    lorentzian = 1.0 / (1.0 + np.square(reduced))
    mixture: np.ndarray = eta * lorentzian + (1.0 - eta) * gaussian
    return mixture


def split_pseudo_voigt_profile(
    two_theta_deg: Any,
    *,
    centre_deg: float,
    fwhm_left_deg: float,
    fwhm_right_deg: float,
    eta: float,
) -> np.ndarray:
    """Return the unit-height split pseudo-Voigt profile on ``two_theta_deg``.

    Purpose
    -------
    Model the asymmetry that axial divergence and flat-specimen aberration
    impose on low-angle reflections. A symmetric shape fitted to an asymmetric
    peak returns a centroid displaced towards the tail, and that displacement
    propagates straight into the lattice parameter.

    Method
    ------
    Two half-profiles of :func:`pseudo_voigt_profile` share a centre and a
    Lorentzian fraction but carry independent widths below and above it
    (Toraya 1986). The result is continuous at the centre, where both halves
    equal one, though its derivative is not.

    Parameters
    ----------
    two_theta_deg
        Angles at which to evaluate the profile, in degrees.
    centre_deg
        Peak centre in degrees ``2*theta``. This is the mode, not the centroid.
    fwhm_left_deg, fwhm_right_deg
        Full widths at half maximum of the low- and high-angle halves.
    eta
        Lorentzian fraction in ``[0, 1]``, shared by both halves.

    Returns
    -------
    np.ndarray
        Profile values with the same shape as ``two_theta_deg``.

    Raises
    ------
    ValueError
        If either width is not positive or ``eta`` lies outside ``[0, 1]``.
    """

    axis = np.asarray(two_theta_deg, dtype=np.float64)
    left = pseudo_voigt_profile(axis, centre_deg=centre_deg, fwhm_deg=fwhm_left_deg, eta=eta)
    right = pseudo_voigt_profile(axis, centre_deg=centre_deg, fwhm_deg=fwhm_right_deg, eta=eta)
    return np.where(axis < centre_deg, left, right)


def pseudo_voigt_area(*, height: float, fwhm_deg: float, eta: float) -> float:
    """Return the integrated area of a height-normalized pseudo-Voigt peak.

    Purpose
    -------
    Convert a fitted height and width into the integrated intensity that
    structure-factor and quantitative work needs, without numerical
    quadrature.

    Method
    ------
    A unit-height Gaussian of FWHM ``w`` encloses ``w * sqrt(pi / (4 ln 2))``
    and a unit-height Lorentzian of the same FWHM encloses ``w * pi / 2``; the
    pseudo-Voigt is their linear mixture, so its area is the same mixture of
    the two.

    Parameters
    ----------
    height
        Peak height above background.
    fwhm_deg
        Full width at half maximum in degrees ``2*theta``.
    eta
        Lorentzian fraction in ``[0, 1]``.

    Returns
    -------
    float
        Integrated area in (intensity unit) x degrees.
    """

    shape_area = eta * _LORENTZIAN_AREA + (1.0 - eta) * _GAUSSIAN_AREA
    return float(height) * float(fwhm_deg) * shape_area


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PeakFit:
    """One fitted powder reflection, with the uncertainty on its position.

    Purpose
    -------
    Carry a measured peak position together with everything needed to judge it:
    how well it is determined, what shape was assumed, whether the K-alpha2
    partner was modelled, and how well the model actually described the data.
    Lattice-parameter determination weights reflections by
    ``two_theta_standard_uncertainty_deg``, so this is a load-bearing number,
    not decoration.

    Attributes
    ----------
    two_theta_deg : float
        Fitted K-alpha1 peak position in degrees ``2*theta``. When the doublet
        was modelled this is the alpha1 component alone; when it was not, it is
        the position of the single fitted profile, which for an unresolved
        doublet sits between the two lines.
    two_theta_standard_uncertainty_deg : float
        Estimated standard uncertainty of the position, from the fit
        covariance scaled by the reduced chi-squared of the window.
    height : float
        Peak height above the local background, in the measured intensity unit.
    integrated_intensity : float
        Area of the K-alpha1 component in (intensity unit) x degrees.
    fwhm_deg : float
        Full width at half maximum. For a split profile this is the mean of
        the two half-widths, which is the width a symmetric fit would report.
    fwhm_left_deg, fwhm_right_deg : float
        The two half-widths. Equal for a symmetric ``"pseudo_voigt"`` fit.
    eta : float
        Fitted Lorentzian fraction in ``[0, 1]``.
    shape : str
        ``"pseudo_voigt"`` or ``"split_pseudo_voigt"``.
    doublet_modelled : bool
        Whether a K-alpha2 partner was included at the Bragg-law position.
    background_intercept, background_slope : float
        The straight line subtracted under the peak, evaluated as
        ``intercept + slope * (2*theta - centre)``.
    reduced_chi_squared : float
        Weighted residual sum of squares per degree of freedom in the fit
        window. Values far from one mean the shape, the weights, or the
        background is wrong.
    point_count : int
        Number of measured points inside the fit window.
    window_deg : tuple[float, float]
        Angular limits of the fit window.
    converged : bool
        Whether the optimizer reported success. A non-converged fit is
        retained, not discarded, so the caller can see what happened.
    """

    two_theta_deg: float
    two_theta_standard_uncertainty_deg: float
    height: float
    integrated_intensity: float
    fwhm_deg: float
    fwhm_left_deg: float
    fwhm_right_deg: float
    eta: float
    shape: PeakShape
    doublet_modelled: bool
    background_intercept: float
    background_slope: float
    reduced_chi_squared: float
    point_count: int
    window_deg: tuple[float, float]
    converged: bool

    def __post_init__(self) -> None:
        if self.shape not in PEAK_SHAPES:
            raise ValueError(f"PeakFit.shape must be one of {PEAK_SHAPES}.")
        if not np.isfinite(self.two_theta_deg):
            raise ValueError("PeakFit.two_theta_deg must be finite.")
        if not 0.0 < self.two_theta_deg < 180.0:
            raise ValueError("PeakFit.two_theta_deg must lie strictly inside (0, 180) degrees.")
        if not np.isfinite(self.two_theta_standard_uncertainty_deg):
            raise ValueError("PeakFit.two_theta_standard_uncertainty_deg must be finite.")
        if self.two_theta_standard_uncertainty_deg <= 0.0:
            raise ValueError("PeakFit.two_theta_standard_uncertainty_deg must be positive.")
        if self.fwhm_left_deg <= 0.0 or self.fwhm_right_deg <= 0.0:
            raise ValueError("PeakFit half-widths must be strictly positive.")
        if not 0.0 <= self.eta <= 1.0:
            raise ValueError("PeakFit.eta must lie in [0, 1].")
        if self.point_count < 1:
            raise ValueError("PeakFit.point_count must be at least one.")

    @property
    def asymmetry(self) -> float:
        """Return the width ratio ``fwhm_right / fwhm_left``.

        One means symmetric. Values below one are the usual low-angle case,
        where axial divergence draws out the low-angle side.
        """

        return float(self.fwhm_right_deg / self.fwhm_left_deg)

    def d_spacing_angstrom(self, wavelength_angstrom: float) -> float:
        """Return the interplanar spacing implied by this position.

        Purpose
        -------
        Apply Bragg's law to the fitted position, for the many workflows that
        want ``d`` rather than an angle.

        Parameters
        ----------
        wavelength_angstrom
            The K-alpha1 wavelength the position was fitted against.

        Returns
        -------
        float
            ``d = lambda / (2 sin(theta))`` in angstrom.

        Raises
        ------
        ValueError
            If the wavelength is not positive.
        """

        if not np.isfinite(wavelength_angstrom) or wavelength_angstrom <= 0.0:
            raise ValueError("d_spacing_angstrom requires a positive wavelength.")
        theta = np.deg2rad(0.5 * self.two_theta_deg)
        return float(wavelength_angstrom / (2.0 * np.sin(theta)))

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this fit."""

        return {
            "schema": PEAK_FIT_SCHEMA,
            "two_theta_deg": float(self.two_theta_deg),
            "two_theta_standard_uncertainty_deg": float(
                self.two_theta_standard_uncertainty_deg
            ),
            "height": float(self.height),
            "integrated_intensity": float(self.integrated_intensity),
            "fwhm_deg": float(self.fwhm_deg),
            "fwhm_left_deg": float(self.fwhm_left_deg),
            "fwhm_right_deg": float(self.fwhm_right_deg),
            "asymmetry": self.asymmetry,
            "eta": float(self.eta),
            "shape": self.shape,
            "doublet_modelled": bool(self.doublet_modelled),
            "background_intercept": float(self.background_intercept),
            "background_slope": float(self.background_slope),
            "reduced_chi_squared": float(self.reduced_chi_squared),
            "point_count": int(self.point_count),
            "window_deg": [float(self.window_deg[0]), float(self.window_deg[1])],
            "converged": bool(self.converged),
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this fit."""

        doublet = (
            "The K-alpha1/K-alpha2 pair was modelled jointly, so the quoted position is the "
            "alpha1 line"
            if self.doublet_modelled
            else "A single line was fitted; if the radiation is a doublet this position lies "
            "between alpha1 and alpha2 and must not be used for precise parameters"
        )
        shape = (
            "a symmetric pseudo-Voigt"
            if self.shape == "pseudo_voigt"
            else (
                f"a split pseudo-Voigt with width ratio {self.asymmetry:.3f}, which absorbs the "
                "low-angle axial-divergence asymmetry that would otherwise displace the centre"
            )
        )
        quality = (
            "The reduced chi-squared is near unity, so the model describes the window within "
            "counting statistics"
            if 0.5 <= self.reduced_chi_squared <= 2.0
            else (
                "The reduced chi-squared is far from unity, so the profile shape, the weights or "
                "the local background does not describe this window; treat the position "
                "uncertainty as a lower bound"
            )
        )
        convergence = (
            ""
            if self.converged
            else " The optimizer did not report convergence; this fit is reported for inspection "
            "and should not be used quantitatively."
        )
        return (
            f"Reflection fitted at {self.two_theta_deg:.5f} +/- "
            f"{self.two_theta_standard_uncertainty_deg:.5f} degrees 2*theta using {shape} over "
            f"{self.point_count} points spanning {self.window_deg[0]:.3f} to "
            f"{self.window_deg[1]:.3f} degrees. {doublet}. Height is {self.height:.4g} and "
            f"integrated intensity {self.integrated_intensity:.4g} (intensity unit x degrees), "
            f"with FWHM {self.fwhm_deg:.5f} degrees and Lorentzian fraction {self.eta:.3f}. "
            f"{quality} (reduced chi-squared {self.reduced_chi_squared:.3f}). "
            f"Position uncertainty is the fit covariance scaled by that reduced chi-squared. "
            f"{_CITATION_TORAYA}{convergence}"
        )


@dataclass(frozen=True, slots=True)
class PeakTable:
    """An ordered set of fitted reflections from one measured pattern.

    Purpose
    -------
    Be the single input type of every position-based analysis in the library:
    indexing, lattice-parameter determination, size-strain separation. It
    carries the radiation the positions were measured with, because a position
    without its wavelength is not a measurement of anything.

    Attributes
    ----------
    name : str
        A human name, normally derived from the source pattern.
    peaks : tuple[PeakFit, ...]
        The fits, sorted by ascending ``2*theta``.
    radiation : RadiationSpec | None
        The radiation used. Required by every method that converts positions
        to ``d`` spacings.
    source_name : str
        Name of the :class:`~pytex.diffraction.xrd_measurement.MeasuredPowderPattern`
        the peaks were fitted to.
    settings : Mapping[str, float | str | bool]
        The detection and fitting settings that produced the table, retained so
        a result is re-checkable.
    """

    name: str
    peaks: tuple[PeakFit, ...]
    radiation: RadiationSpec | None = None
    source_name: str = ""
    settings: Mapping[str, float | str | bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("PeakTable.name must be non-empty.")
        peaks = tuple(self.peaks)
        if any(not isinstance(peak, PeakFit) for peak in peaks):
            raise TypeError("PeakTable.peaks must contain PeakFit instances.")
        positions = [peak.two_theta_deg for peak in peaks]
        if positions != sorted(positions):
            peaks = tuple(sorted(peaks, key=lambda peak: peak.two_theta_deg))
        object.__setattr__(self, "peaks", peaks)
        object.__setattr__(self, "settings", MappingProxyType(dict(self.settings)))

    def __len__(self) -> int:
        return len(self.peaks)

    def __iter__(self) -> Any:
        return iter(self.peaks)

    def __getitem__(self, index: int) -> PeakFit:
        return self.peaks[index]

    @property
    def two_theta_deg(self) -> np.ndarray:
        """Return the fitted positions as a read-only array, in degrees."""

        return as_float_array([peak.two_theta_deg for peak in self.peaks], shape=(None,))

    @property
    def standard_uncertainty_deg(self) -> np.ndarray:
        """Return the per-position standard uncertainties, in degrees."""

        return as_float_array(
            [peak.two_theta_standard_uncertainty_deg for peak in self.peaks], shape=(None,)
        )

    @property
    def integrated_intensity(self) -> np.ndarray:
        """Return the integrated intensities as a read-only array."""

        return as_float_array([peak.integrated_intensity for peak in self.peaks], shape=(None,))

    @property
    def converged_count(self) -> int:
        """Return how many fits the optimizer reported as converged."""

        return int(sum(1 for peak in self.peaks if peak.converged))

    def d_spacing_angstrom(self) -> np.ndarray:
        """Return the Bragg ``d`` spacings of every peak, in angstrom.

        Raises
        ------
        ValueError
            If the table carries no radiation, since without a wavelength a
            position cannot be converted.
        """

        if self.radiation is None:
            raise ValueError(
                "PeakTable.d_spacing_angstrom requires a radiation; the table was built without "
                "one, so its positions cannot be converted to spacings."
            )
        wavelength = self.radiation.wavelength_angstrom
        return as_float_array(
            [peak.d_spacing_angstrom(wavelength) for peak in self.peaks], shape=(None,)
        )

    def filter_converged(self) -> PeakTable:
        """Return a copy holding only the converged fits.

        Purpose
        -------
        Provide the one-line guard that quantitative workflows should apply
        before consuming a table, without silently dropping information at the
        point of fitting.
        """

        kept = tuple(peak for peak in self.peaks if peak.converged)
        if not kept:
            raise ValueError(
                "No peak in this table converged, so there is nothing to analyse. Check the "
                "fit window, the starting width, and whether the detected positions are real."
            )
        return PeakTable(
            name=f"{self.name} (converged only)",
            peaks=kept,
            radiation=self.radiation,
            source_name=self.source_name,
            settings=self.settings,
        )

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this table."""

        return {
            "schema": PEAK_TABLE_SCHEMA,
            "name": self.name,
            "source_name": self.source_name,
            "radiation": None
            if self.radiation is None
            else {
                "name": self.radiation.name,
                "wavelength_angstrom": float(self.radiation.wavelength_angstrom),
                "kalpha2_wavelength_angstrom": (
                    None
                    if self.radiation.kalpha2_wavelength_angstrom is None
                    else float(self.radiation.kalpha2_wavelength_angstrom)
                ),
            },
            "settings": {key: value for key, value in self.settings.items()},
            "peaks": [peak.to_json() for peak in self.peaks],
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this table."""

        if not self.peaks:
            return (
                f"Peak table '{self.name}' is empty: no reflection was detected or fitted in "
                f"'{self.source_name}'."
            )
        positions = self.two_theta_deg
        uncertainties = self.standard_uncertainty_deg
        radiation = (
            f"Positions are referred to {self.radiation.name} at "
            f"{self.radiation.wavelength_angstrom:.6f} angstrom."
            if self.radiation is not None
            else "No radiation was declared, so these positions cannot be converted to spacings."
        )
        failures = len(self.peaks) - self.converged_count
        convergence = (
            "Every fit converged."
            if failures == 0
            else f"{failures} of {len(self.peaks)} fits did not converge and are flagged."
        )
        return (
            f"Peak table '{self.name}' holds {len(self.peaks)} fitted reflections from "
            f"'{self.source_name}', spanning {positions[0]:.4f} to {positions[-1]:.4f} degrees "
            f"2*theta. Median position uncertainty is "
            f"{float(np.median(uncertainties)):.5f} degrees. {radiation} {convergence} "
            f"Precision of a lattice parameter derived from these positions improves as "
            f"cot(theta) falls, so the highest-angle reflections carry the most weight. "
            f"{_CITATION_CULLITY}"
        )


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def _ricker_kernel(sigma_points: float, *, truncate: float = 4.0) -> np.ndarray:
    """Return a unit-L2-norm, zero-mean Ricker kernel of the given scale.

    Zero mean makes the response blind to a constant offset; the even symmetry
    makes it nearly blind to a linear ramp. Unit L2 norm means that filtering
    white noise of unit variance produces a response of unit variance, which is
    what lets the detection threshold be quoted in noise sigmas.
    """

    half = int(max(3, np.ceil(truncate * sigma_points)))
    offsets = np.arange(-half, half + 1, dtype=np.float64)
    reduced = offsets / sigma_points
    kernel = (1.0 - np.square(reduced)) * np.exp(-0.5 * np.square(reduced))
    kernel = kernel - kernel.mean()
    norm = float(np.linalg.norm(kernel))
    if norm <= 0.0:  # pragma: no cover - unreachable for sigma_points >= 1
        raise ValueError("Degenerate Ricker kernel.")
    normalized: np.ndarray = kernel / norm
    return normalized


def _anscombe(intensity: np.ndarray) -> np.ndarray:
    """Return the Anscombe variance-stabilizing transform of Poisson counts.

    ``2 sqrt(I + 3/8)`` has a variance very close to one for any mean above a
    few counts, which turns heteroscedastic counting noise into homoscedastic
    noise the matched filter can threshold uniformly.
    """

    stabilized: np.ndarray = 2.0 * np.sqrt(np.maximum(intensity, 0.0) + 0.375)
    return stabilized


def _expected_fwhm_deg(
    two_theta_deg: np.ndarray,
    *,
    instrument: InstrumentBroadening | None,
    expected_fwhm_deg: float | None,
) -> np.ndarray:
    """Return the expected peak FWHM at every angle, in degrees."""

    if expected_fwhm_deg is not None:
        if not np.isfinite(expected_fwhm_deg) or expected_fwhm_deg <= 0.0:
            raise ValueError("expected_fwhm_deg must be finite and positive when given.")
        return np.full_like(two_theta_deg, float(expected_fwhm_deg))
    if instrument is not None:
        widths = np.asarray(instrument.fwhm_deg(two_theta_deg), dtype=np.float64)
        return np.broadcast_to(widths, two_theta_deg.shape).copy()
    # With neither a calibrated instrument nor a stated width, assume a
    # laboratory Bragg-Brentano default. It is deliberately a single number
    # rather than a fitted guess: detection only needs the scale to within a
    # factor of about two, and a stated default is auditable.
    return np.full_like(two_theta_deg, 0.12)


def detect_peaks(
    measured: MeasuredPowderPattern,
    *,
    instrument: InstrumentBroadening | None = None,
    expected_fwhm_deg: float | None = None,
    prominence_sigma: float = 5.0,
    subtract_background: bool = True,
    background_half_window_deg: float = 2.0,
    two_theta_range_deg: tuple[float, float] | None = None,
    max_peaks: int = 128,
    radiation: RadiationSpec | None = None,
    suppress_kalpha2: bool = True,
) -> tuple[float, ...]:
    """Locate candidate reflections by matching the instrument's own width.

    Purpose
    -------
    Turn a measured diffractogram into a list of angles worth fitting, with a
    threshold expressed in units of the noise rather than in counts, and with a
    width criterion taken from the diffractometer rather than guessed.

    Method
    ------
    1. The slowly varying background is removed with SNIP
       (:func:`~pytex.diffraction.xrd_background.estimate_background`), because
       a matched filter run on top of a rising background reports the
       background's curvature.
    2. Counts are Anscombe-transformed, ``y = 2 sqrt(I + 3/8)``, so that Poisson
       noise has unit variance at every count level and one threshold applies
       across a pattern whose strongest and weakest peaks differ by decades.
    3. The transformed profile is convolved with zero-mean, unit-L2-norm Ricker
       (negative second derivative of a Gaussian) kernels over a grid of scales,
       and the response at each point is interpolated to *that point's* expected
       width, taken from the calibrated Caglioti curve. A feature narrower than
       the instrument can produce -- a cosmic-ray spike, a dead channel -- or
       broader than it -- background structure, an amorphous halo -- gives a
       weak response and is rejected without any threshold in angle units.
    4. Local maxima of the response above ``prominence_sigma`` times the robust
       (median-absolute-deviation) noise level of the response are kept, refined
       to sub-step precision by parabolic interpolation, and thinned so that no
       two survivors sit closer than one expected FWHM.
    5. When the radiation is a doublet, each accepted candidate also suppresses
       the angle its own K-alpha2 partner would sit at. Above roughly 90 degrees
       ``2*theta`` the pair separates far enough for the filter to report the
       alpha2 line as a peak in its own right, and an alpha2 line admitted to
       the candidate list becomes a reflection at the wrong ``d`` spacing --
       which is a bias in the lattice parameter, not merely a spurious row.

    Detection is deliberately separate from fitting: it answers "where should I
    look", and :func:`fit_peaks` answers "what is there, and how well do I know
    it".

    Parameters
    ----------
    measured
        The measured pattern. Pass the raw profile; background removal is part
        of the method.
    instrument
        A calibrated resolution function. When given, the filter scale tracks
        the true width at every angle, which is the configuration this method
        is designed for.
    expected_fwhm_deg
        A single width to use instead of ``instrument``. Overrides it when
        both are given.
    prominence_sigma
        Detection threshold in robust noise standard deviations of the filter
        response. Five is conservative; lower it to chase weak reflections and
        expect to reject more of them at the fitting stage.
    subtract_background
        Remove a SNIP background before filtering. Turn this off only for an
        already-subtracted pattern.
    background_half_window_deg
        SNIP half-window, in degrees. Set it wider than the broadest peak.
    two_theta_range_deg
        Restrict detection to this angular window.
    max_peaks
        Keep at most this many candidates, strongest response first. The
        returned tuple is still sorted by angle.
    radiation
        Radiation whose K-alpha2 line drives the doublet suppression of step 5.
        Falls back to ``measured.radiation``.
    suppress_kalpha2
        Discard candidates that coincide with the K-alpha2 partner of a
        stronger candidate. Turn it off only to inspect the raw detection.

    Returns
    -------
    tuple[float, ...]
        Candidate peak positions in degrees ``2*theta``, ascending. May be
        empty, which is a legitimate answer for a featureless scan.

    Raises
    ------
    ValueError
        If a setting is outside its valid range, or the requested angular
        window contains too few points to filter.

    See Also
    --------
    fit_peaks : refine these candidates into positions with uncertainties.
    detect_and_fit_peaks : both steps in one call.

    Notes
    -----
    The scale-matched Ricker filter is the continuous-wavelet peak detector of
    Du et al. (2006), specialized to the case where the correct scale is known
    in advance from the instrument rather than searched over. That knowledge is
    what removes the free parameters.
    """

    if not np.isfinite(prominence_sigma) or prominence_sigma <= 0.0:
        raise ValueError("detect_peaks requires a finite, positive prominence_sigma.")
    if max_peaks < 1:
        raise ValueError("detect_peaks requires max_peaks >= 1.")

    axis = np.asarray(measured.two_theta_deg, dtype=np.float64)
    intensity = np.asarray(measured.intensity, dtype=np.float64)

    if subtract_background:
        background = estimate_background(
            measured, method="snip", half_window_deg=background_half_window_deg
        )
        intensity = np.maximum(intensity - background.background, 0.0)

    if two_theta_range_deg is not None:
        low, high = float(two_theta_range_deg[0]), float(two_theta_range_deg[1])
        if not low < high:
            raise ValueError("two_theta_range_deg must be an increasing (low, high) pair.")
        window = (axis >= low) & (axis <= high)
        axis = axis[window]
        intensity = intensity[window]

    if axis.size < 16:
        raise ValueError(
            "detect_peaks needs at least sixteen points inside the analysed window; the scan or "
            "the requested range is too short to filter."
        )

    step_deg = float(np.median(np.diff(axis)))
    widths_deg = _expected_fwhm_deg(
        axis, instrument=instrument, expected_fwhm_deg=expected_fwhm_deg
    )
    # A Gaussian of FWHM w has sigma = w / (2 sqrt(2 ln 2)); the matched Ricker
    # scale is that sigma expressed in data points.
    sigma_points = widths_deg / (step_deg * 2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_points = np.clip(sigma_points, 1.0, max(1.0, axis.size / 8.0))

    transformed = (
        _anscombe(intensity) if measured.intensity_unit == "counts" else intensity.copy()
    )
    transformed = transformed - float(np.median(transformed))

    scale_grid = np.unique(
        np.geomspace(float(sigma_points.min()), float(sigma_points.max()) * 1.0 + 1e-9, num=8)
    )
    responses = np.empty((scale_grid.size, axis.size), dtype=np.float64)
    for row, scale in enumerate(scale_grid):
        kernel = _ricker_kernel(float(scale))
        responses[row] = np.convolve(transformed, kernel, mode="same")

    if scale_grid.size == 1:
        response = responses[0]
    else:
        # Interpolate the scale axis at each point's own expected scale.
        upper = np.searchsorted(scale_grid, sigma_points, side="left")
        upper = np.clip(upper, 1, scale_grid.size - 1)
        lower = upper - 1
        span = scale_grid[upper] - scale_grid[lower]
        blend = np.where(span > 0.0, (sigma_points - scale_grid[lower]) / span, 0.0)
        columns = np.arange(axis.size)
        response = (1.0 - blend) * responses[lower, columns] + blend * responses[upper, columns]

    noise = 1.4826 * float(np.median(np.abs(response - np.median(response))))
    if noise <= 0.0:
        noise = float(np.std(response)) or 1.0
    threshold = prominence_sigma * noise

    interior = response[1:-1]
    is_maximum = (interior > response[:-2]) & (interior >= response[2:]) & (interior > threshold)
    candidates = np.flatnonzero(is_maximum) + 1
    if candidates.size == 0:
        return ()

    # Parabolic refinement on the response, which is smooth by construction.
    left = response[candidates - 1]
    centre = response[candidates]
    right = response[candidates + 1]
    curvature = left - 2.0 * centre + right
    offset = np.where(np.abs(curvature) > 0.0, 0.5 * (left - right) / curvature, 0.0)
    offset = np.clip(offset, -0.5, 0.5)
    positions = axis[candidates] + offset * step_deg
    strengths = centre

    spec = radiation if radiation is not None else measured.radiation
    wavelength_ratio: float | None = None
    if (
        suppress_kalpha2
        and spec is not None
        and spec.kalpha2_wavelength_angstrom is not None
    ):
        wavelength_ratio = float(spec.kalpha2_wavelength_angstrom / spec.wavelength_angstrom)

    # Thin by strength: a strong peak suppresses weaker maxima within one FWHM,
    # and, for doublet radiation, also the angle its own K-alpha2 line occupies.
    order = np.argsort(strengths)[::-1]
    kept: list[float] = []
    forbidden: list[tuple[float, float]] = []
    for index in order:
        position = float(positions[index])
        separation = float(widths_deg[candidates[index]])
        if any(abs(position - previous) < separation for previous in kept):
            continue
        if any(abs(position - centre) < tolerance for centre, tolerance in forbidden):
            continue
        kept.append(position)
        if wavelength_ratio is not None:
            try:
                partner = _kalpha2_position_deg(position, wavelength_ratio)
            except ValueError:
                partner = None
            if partner is not None:
                forbidden.append((partner, separation))
        if len(kept) >= max_peaks:
            break
    return tuple(sorted(kept))


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------


def _kalpha2_position_deg(two_theta_one_deg: float, ratio: float) -> float:
    """Return the K-alpha2 position for a K-alpha1 position at the same ``d``.

    Bragg's law at fixed ``d`` gives ``sin(theta2) = (lambda2 / lambda1)
    sin(theta1)``, so the doublet separation is not a constant offset: it grows
    as ``tan(theta)`` and is what resolves the pair at high angle.
    """

    argument = ratio * np.sin(np.deg2rad(0.5 * two_theta_one_deg))
    if not -1.0 <= argument <= 1.0:
        raise ValueError(
            "The K-alpha2 partner of this reflection is beyond the Ewald limit; the position or "
            "the wavelength pair is inconsistent."
        )
    return float(np.rad2deg(2.0 * np.arcsin(argument)))


def _weights(
    intensity: np.ndarray,
    uncertainty: np.ndarray | None,
    unit: str,
) -> np.ndarray:
    """Return least-squares weights ``1/sigma`` for a fit window."""

    if uncertainty is not None:
        stated: np.ndarray = 1.0 / uncertainty
        return stated
    if unit == "counts":
        poisson: np.ndarray = 1.0 / np.sqrt(np.maximum(intensity, 1.0))
        return poisson
    return np.ones_like(intensity)


def _fit_one_peak(
    axis: np.ndarray,
    intensity: np.ndarray,
    uncertainty: np.ndarray | None,
    *,
    unit: str,
    centre_guess: float,
    fwhm_guess: float,
    shape: PeakShape,
    doublet: tuple[float, float] | None,
) -> PeakFit:
    """Fit a single window and return its :class:`PeakFit`.

    ``doublet`` is ``(lambda2 / lambda1, I2 / I1)`` when the K-alpha2 partner is
    to be modelled, and ``None`` for a single line.
    """

    weights = _weights(intensity, uncertainty, unit)
    baseline = float(np.min(intensity))
    height_guess = max(float(np.max(intensity)) - baseline, 1e-9)

    split = shape == "split_pseudo_voigt"
    # [centre, height, fwhm_left, (fwhm_right,) eta, bg_intercept, bg_slope]
    start = [centre_guess, height_guess, fwhm_guess]
    lower = [axis[0], 0.0, 1e-4]
    upper = [axis[-1], np.inf, float(axis[-1] - axis[0])]
    if split:
        start.append(fwhm_guess)
        lower.append(1e-4)
        upper.append(float(axis[-1] - axis[0]))
    start.extend([0.5, baseline, 0.0])
    lower.extend([0.0, -np.inf, -np.inf])
    upper.extend([1.0, np.inf, np.inf])

    def model(parameters: np.ndarray) -> np.ndarray:
        centre = float(parameters[0])
        height = float(parameters[1])
        fwhm_left = float(parameters[2])
        if split:
            fwhm_right = float(parameters[3])
            eta = float(parameters[4])
            background = float(parameters[5]) + float(parameters[6]) * (axis - centre)
        else:
            fwhm_right = fwhm_left
            eta = float(parameters[3])
            background = float(parameters[4]) + float(parameters[5]) * (axis - centre)
        profile = split_pseudo_voigt_profile(
            axis,
            centre_deg=centre,
            fwhm_left_deg=fwhm_left,
            fwhm_right_deg=fwhm_right,
            eta=eta,
        )
        total = height * profile
        if doublet is not None:
            wavelength_ratio, intensity_ratio = doublet
            partner = _kalpha2_position_deg(centre, wavelength_ratio)
            total = total + height * intensity_ratio * split_pseudo_voigt_profile(
                axis,
                centre_deg=partner,
                fwhm_left_deg=fwhm_left,
                fwhm_right_deg=fwhm_right,
                eta=eta,
            )
        return total + background

    def residual(parameters: np.ndarray) -> np.ndarray:
        weighted: np.ndarray = (model(parameters) - intensity) * weights
        return weighted

    solution = least_squares(
        residual,
        x0=np.asarray(start, dtype=np.float64),
        bounds=(np.asarray(lower, dtype=np.float64), np.asarray(upper, dtype=np.float64)),
        method="trf",
        max_nfev=4000,
    )

    parameters = solution.x
    degrees_of_freedom = max(axis.size - parameters.size, 1)
    residual_sum = float(np.sum(np.square(solution.fun)))
    reduced_chi_squared = residual_sum / degrees_of_freedom

    # Goodness-of-fit-scaled covariance: (J^T J)^-1 times the reduced chi-squared.
    jacobian = np.asarray(solution.jac, dtype=np.float64)
    normal = jacobian.T @ jacobian
    try:
        covariance = np.linalg.inv(normal) * reduced_chi_squared
        position_variance = float(covariance[0, 0])
    except np.linalg.LinAlgError:  # pragma: no cover - singular only for degenerate windows
        position_variance = np.nan
    if not np.isfinite(position_variance) or position_variance <= 0.0:
        # A degenerate window still deserves a number: fall back to the
        # width-over-signal-to-noise scaling that a single-peak position
        # uncertainty obeys asymptotically.
        position_variance = float(np.square(fwhm_guess / max(np.sqrt(axis.size), 1.0)))

    centre = float(parameters[0])
    height = float(parameters[1])
    fwhm_left = float(parameters[2])
    if split:
        fwhm_right = float(parameters[3])
        eta = float(parameters[4])
        intercept = float(parameters[5])
        slope = float(parameters[6])
    else:
        fwhm_right = fwhm_left
        eta = float(parameters[3])
        intercept = float(parameters[4])
        slope = float(parameters[5])

    mean_fwhm = 0.5 * (fwhm_left + fwhm_right)
    area = 0.5 * (
        pseudo_voigt_area(height=height, fwhm_deg=fwhm_left, eta=eta)
        + pseudo_voigt_area(height=height, fwhm_deg=fwhm_right, eta=eta)
    )

    return PeakFit(
        two_theta_deg=centre,
        two_theta_standard_uncertainty_deg=float(np.sqrt(position_variance)),
        height=height,
        integrated_intensity=area,
        fwhm_deg=mean_fwhm,
        fwhm_left_deg=fwhm_left,
        fwhm_right_deg=fwhm_right,
        eta=eta,
        shape=shape,
        doublet_modelled=doublet is not None,
        background_intercept=intercept,
        background_slope=slope,
        reduced_chi_squared=reduced_chi_squared,
        point_count=int(axis.size),
        window_deg=(float(axis[0]), float(axis[-1])),
        converged=bool(solution.success),
    )


def _merge_coincident_fits(
    fits: list[PeakFit], *, tolerance_fraction: float
) -> list[PeakFit]:
    """Collapse fits that converged to the same angle into one reflection.

    Two windows centred on different candidates can slide onto the same strong
    line. The duplicate is not merely redundant: an indexer matching one-to-one
    would give the two copies different ``(hkl)``, inventing a reflection at a
    spacing nothing diffracted from. The fit with the lower reduced
    chi-squared wins, because it is the one whose window actually described its
    data.
    """

    if tolerance_fraction <= 0.0 or len(fits) < 2:
        return fits
    ordered = sorted(fits, key=lambda item: item.two_theta_deg)
    kept: list[PeakFit] = [ordered[0]]
    for candidate in ordered[1:]:
        previous = kept[-1]
        separation = abs(candidate.two_theta_deg - previous.two_theta_deg)
        width = 0.5 * (candidate.fwhm_deg + previous.fwhm_deg)
        if separation < tolerance_fraction * width:
            if candidate.reduced_chi_squared < previous.reduced_chi_squared:
                kept[-1] = candidate
        else:
            kept.append(candidate)
    return kept


def fit_peaks(
    measured: MeasuredPowderPattern,
    centres_deg: Iterable[float],
    *,
    radiation: RadiationSpec | None = None,
    instrument: InstrumentBroadening | None = None,
    expected_fwhm_deg: float | None = None,
    shape: PeakShape = "pseudo_voigt",
    model_doublet: bool = True,
    window_fwhm: float = 4.0,
    merge_tolerance_fwhm: float = 0.25,
    name: str | None = None,
) -> PeakTable:
    """Fit a profile to each candidate position and return positions with ESDs.

    Purpose
    -------
    Convert candidate angles into measured reflection positions whose
    uncertainties are good enough to weight a precise lattice-parameter
    determination.

    Method
    ------
    Each candidate gets its own window of ``window_fwhm`` expected widths on
    either side, fitted independently by bounded Levenberg-Marquardt
    (``scipy.optimize.least_squares``, trust-region reflective) with:

    * a pseudo-Voigt or split pseudo-Voigt profile,
    * a straight local background under the window,
    * and, when the radiation declares a K-alpha2 line and ``model_doublet`` is
      set, a partner peak whose position is *not* free: it is fixed by Bragg's
      law at the same ``d`` spacing, ``sin(theta2) = (lambda2 / lambda1)
      sin(theta1)``, with the tabulated intensity ratio and shared width and
      mixing. This costs no extra parameter and removes the systematic shift
      that fitting a single symmetric line to an unresolved doublet produces.

    Weights are ``1/sigma`` from the pattern's own standard uncertainties when
    present, from Poisson counting statistics when the intensity unit is
    counts, and unity otherwise.

    Finally, fits that converged to the *same* angle are merged. Two candidates
    can be far enough apart to survive detection -- for instance where a weak
    reflection is predicted beside a strong one -- and then both windows slide
    onto the strong line. Two fits at one angle are one reflection, and leaving
    both in place would let an indexer pair each of them with a different
    ``(hkl)``, which puts a reflection at a spacing nothing diffracted from.
    The survivor is the fit with the lower reduced chi-squared.

    Parameters
    ----------
    measured
        The measured pattern, background included: the fit removes a local
        straight background itself.
    centres_deg
        Candidate positions, normally from :func:`detect_peaks`.
    radiation
        Radiation for the table, and the source of the K-alpha2 wavelength.
        Falls back to ``measured.radiation``.
    instrument
        Calibrated resolution function, used for the starting width and the
        window size.
    expected_fwhm_deg
        A single starting width, overriding ``instrument``.
    shape
        ``"pseudo_voigt"`` or ``"split_pseudo_voigt"``. Prefer the split shape
        below about 40 degrees ``2*theta``, where axial divergence is
        asymmetric enough to move a symmetric centroid.
    model_doublet
        Model the K-alpha2 partner. Has no effect when the radiation declares
        no second line.
    window_fwhm
        Half-width of each fit window in expected FWHM. Too small starves the
        background; too large invites the neighbouring reflection in.
    merge_tolerance_fwhm
        Fits whose centres agree to within this fraction of their mean FWHM are
        treated as one reflection. Set it to zero to keep every fit, which is
        useful only for inspecting the raw behaviour.
    name
        Name for the returned table.

    Returns
    -------
    PeakTable
        The fits, sorted by angle, with the settings that produced them.

    Raises
    ------
    ValueError
        If ``shape`` is unknown, ``window_fwhm`` is not positive, or no
        candidate lies inside the measured range.

    See Also
    --------
    detect_peaks : produce the candidate positions.
    pytex.diffraction.xrd_lattice_parameter.determine_lattice_parameters :
        the main consumer of the result.

    Notes
    -----
    Fitting the doublet is strictly better than stripping it. Rachinger
    stripping assumes an exact 1:2 ratio and identical shapes, and propagates
    each subtraction into the next, so it raises the noise of every point above
    the first. Modelling makes the same physical assumptions but leaves the
    data untouched, so the uncertainty reported here is the uncertainty of the
    measurement rather than of the measurement plus an arithmetic operation.
    """

    if shape not in PEAK_SHAPES:
        raise ValueError(f"fit_peaks requires shape in {PEAK_SHAPES}.")
    if not np.isfinite(window_fwhm) or window_fwhm <= 0.0:
        raise ValueError("fit_peaks requires a finite, positive window_fwhm.")

    axis = np.asarray(measured.two_theta_deg, dtype=np.float64)
    intensity = np.asarray(measured.intensity, dtype=np.float64)
    uncertainty = (
        None
        if measured.standard_uncertainty is None
        else np.asarray(measured.standard_uncertainty, dtype=np.float64)
    )
    spec = radiation if radiation is not None else measured.radiation

    doublet: tuple[float, float] | None = None
    if model_doublet and spec is not None and spec.kalpha2_wavelength_angstrom is not None:
        doublet = (
            float(spec.kalpha2_wavelength_angstrom / spec.wavelength_angstrom),
            float(spec.kalpha2_relative_intensity),
        )

    candidates = sorted(float(value) for value in centres_deg)
    if not candidates:
        raise ValueError("fit_peaks was given no candidate positions to fit.")

    widths = _expected_fwhm_deg(
        np.asarray(candidates, dtype=np.float64),
        instrument=instrument,
        expected_fwhm_deg=expected_fwhm_deg,
    )

    fits: list[PeakFit] = []
    for candidate, width in zip(candidates, widths, strict=True):
        half_window = window_fwhm * float(width)
        window = (axis >= candidate - half_window) & (axis <= candidate + half_window)
        minimum_points = 8 if shape == "pseudo_voigt" else 10
        if int(np.count_nonzero(window)) < minimum_points:
            continue
        fits.append(
            _fit_one_peak(
                axis[window],
                intensity[window],
                None if uncertainty is None else uncertainty[window],
                unit=measured.intensity_unit,
                centre_guess=candidate,
                fwhm_guess=float(width),
                shape=shape,
                doublet=doublet,
            )
        )

    if not fits:
        raise ValueError(
            "No candidate position had enough measured points inside its fit window. Reduce "
            "window_fwhm, widen the scan, or check that the candidates lie inside the measured "
            "angular range."
        )

    fits = _merge_coincident_fits(fits, tolerance_fraction=merge_tolerance_fwhm)

    return PeakTable(
        name=name or f"{measured.name} peaks",
        peaks=tuple(fits),
        radiation=spec,
        source_name=measured.name,
        settings={
            "shape": shape,
            "model_doublet": doublet is not None,
            "window_fwhm": float(window_fwhm),
            "candidate_count": float(len(candidates)),
        },
    )


def detect_and_fit_peaks(
    measured: MeasuredPowderPattern,
    *,
    radiation: RadiationSpec | None = None,
    instrument: InstrumentBroadening | None = None,
    expected_fwhm_deg: float | None = None,
    prominence_sigma: float = 5.0,
    shape: PeakShape = "pseudo_voigt",
    model_doublet: bool = True,
    window_fwhm: float = 4.0,
    merge_tolerance_fwhm: float = 0.25,
    two_theta_range_deg: tuple[float, float] | None = None,
    max_peaks: int = 128,
    name: str | None = None,
) -> PeakTable:
    """Detect and then fit every reflection in a measured pattern.

    Purpose
    -------
    Provide the one call that a routine analysis makes, without hiding either
    half: the two stages remain available separately for the cases where the
    operator wants to edit the candidate list before fitting.

    Parameters
    ----------
    measured
        The measured pattern.
    radiation
        Radiation for both stages: it suppresses K-alpha2 candidates during
        detection and fixes the modelled partner position during fitting.
    instrument, expected_fwhm_deg, prominence_sigma, two_theta_range_deg, max_peaks
        Passed to :func:`detect_peaks`.
    shape, model_doublet, window_fwhm, merge_tolerance_fwhm, name
        Passed to :func:`fit_peaks`.

    Returns
    -------
    PeakTable
        Fitted reflections with position uncertainties.

    Raises
    ------
    ValueError
        If no candidate survives detection, which for a real diffractogram
        means the threshold is too high rather than that the specimen is
        amorphous.
    """

    candidates = detect_peaks(
        measured,
        instrument=instrument,
        expected_fwhm_deg=expected_fwhm_deg,
        prominence_sigma=prominence_sigma,
        two_theta_range_deg=two_theta_range_deg,
        max_peaks=max_peaks,
        radiation=radiation,
        suppress_kalpha2=model_doublet,
    )
    if not candidates:
        raise ValueError(
            "No reflection was detected above the threshold. Lower prominence_sigma, widen the "
            "background half-window, or confirm the scan contains Bragg peaks."
        )
    return fit_peaks(
        measured,
        candidates,
        radiation=radiation,
        instrument=instrument,
        expected_fwhm_deg=expected_fwhm_deg,
        shape=shape,
        model_doublet=model_doublet,
        window_fwhm=window_fwhm,
        merge_tolerance_fwhm=merge_tolerance_fwhm,
        name=name,
    )


__all__ = [
    "PEAK_FIT_SCHEMA",
    "PEAK_SHAPES",
    "PEAK_TABLE_SCHEMA",
    "PeakFit",
    "PeakShape",
    "PeakTable",
    "detect_and_fit_peaks",
    "detect_peaks",
    "fit_peaks",
    "pseudo_voigt_area",
    "pseudo_voigt_profile",
    "split_pseudo_voigt_profile",
]

