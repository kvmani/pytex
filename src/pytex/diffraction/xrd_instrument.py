"""Instrumental broadening: the width a perfect crystal would still show.

A measured diffraction peak is wider than the sample alone makes it. The
divergence slits, the source size, the flat-specimen and axial-divergence
aberrations, the receiving optics and the detector each add breadth, and their
sum -- the *instrumental resolution function* -- varies systematically with
angle. Nothing about crystallite size or microstrain can be read out of a
measured width until that contribution is removed, and a size derived from an
uncorrected width is not an underestimate of the truth: it is a measurement of
the diffractometer.

This module holds three things and keeps them apart on purpose:

1. :class:`InstrumentBroadening` -- the resolution function itself, as the
   Caglioti Gaussian width plus a Lorentzian term, combined into a pseudo-Voigt
   by the Thompson-Cox-Hastings construction that Rietveld programs use.
2. :func:`calibrate_instrument_broadening` -- fitting that function to the peak
   widths of a *standard* (LaB6, Si, or any specimen whose own broadening is
   negligible). This is a measurement of the instrument and belongs to the
   instrument, not to the sample under study.
3. :func:`deconvolve_instrument_width` and :func:`williamson_hall` -- removing
   the calibrated instrument width from measured widths, and turning what
   remains into a crystallite size and a microstrain.

References
----------
Caglioti, G., Paoletti, A. & Ricci, F. P., *Nucl. Instrum.* **3** (1958) 223-228,
doi:10.1016/0369-643X(58)90029-X -- the ``U tan^2(theta) + V tan(theta) + W``
resolution function.

Thompson, P., Cox, D. E. & Hastings, J. B., *J. Appl. Crystallogr.* **20** (1987)
79-83, doi:10.1107/S0021889887087090 -- the pseudo-Voigt approximation to the
Voigt convolution used here.

Williamson, G. K. & Hall, W. H., *Acta Metall.* **1** (1953) 22-31,
doi:10.1016/0001-6160(53)90006-6 -- size-strain separation.

Scherrer, P., *Nachr. Ges. Wiss. Goettingen* **26** (1918) 98-100; see Langford &
Wilson, *J. Appl. Crystallogr.* **11** (1978) 102-113,
doi:10.1107/S0021889878012844 for the shape-factor discussion.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array

INSTRUMENT_BROADENING_SCHEMA = "pytex.diffraction.instrument_broadening"
WILLIAMSON_HALL_SCHEMA = "pytex.diffraction.williamson_hall_analysis"

DeconvolutionMode = Literal["gaussian", "lorentzian", "pseudo_voigt"]

DECONVOLUTION_MODES = ("gaussian", "lorentzian", "pseudo_voigt")

# Thompson, Cox & Hastings (1987), equations 3 and 4.
_TCH_A = 2.69269
_TCH_B = 2.42843
_TCH_C = 4.47163
_TCH_D = 0.07842
_TCH_ETA = (1.36603, -0.47719, 0.11116)

_CITATION_CAGLIOTI = (
    "Caglioti, Paoletti & Ricci, Nucl. Instrum. 3 (1958) 223, doi:10.1016/0369-643X(58)90029-X."
)
_CITATION_TCH = (
    "Thompson, Cox & Hastings, J. Appl. Crystallogr. 20 (1987) 79, doi:10.1107/S0021889887087090."
)
_CITATION_WILLIAMSON_HALL = (
    "Williamson & Hall, Acta Metall. 1 (1953) 22, doi:10.1016/0001-6160(53)90006-6."
)
_CITATION_LANGFORD_WILSON = (
    "Langford & Wilson, J. Appl. Crystallogr. 11 (1978) 102, doi:10.1107/S0021889878012844."
)

#: Scherrer shape factor for a spherical crystallite measured by FWHM.
#: Langford & Wilson (1978) tabulate this as shape-dependent; 0.9 is the
#: conventional value and the one every reported "Scherrer size" assumes unless
#: it says otherwise.
SCHERRER_SHAPE_FACTOR = 0.9


@dataclass(frozen=True, slots=True)
class InstrumentBroadening:
    """An angular resolution function: the width the instrument itself adds.

    Purpose
    -------
    State, as a function of angle, how wide a peak the diffractometer produces
    from a specimen that contributes no broadening of its own. Everything in
    this library that separates sample physics from instrument physics goes
    through this object.

    Method
    ------
    The Gaussian half of the width follows Caglioti et al. (1958):
    ``FWHM_G^2 = U tan^2(theta) + V tan(theta) + W``. The Lorentzian half
    follows the Rietveld convention ``FWHM_L = X tan(theta) + Y / cos(theta)``,
    where ``X`` carries strain-like and ``Y`` size-like angular dependence. The
    two are combined into one pseudo-Voigt width and mixing parameter by the
    Thompson-Cox-Hastings construction, which approximates the true Voigt
    convolution to better than one percent.

    Attributes
    ----------
    caglioti_u, caglioti_v, caglioti_w : float
        Gaussian resolution coefficients in degrees squared. ``W`` is the
        low-angle floor and must keep the squared width positive across the
        working range.
    lorentzian_x, lorentzian_y : float
        Lorentzian coefficients in degrees.
    zero_shift_deg : float
        Detector zero-point error, added to every calculated ``2*theta``. It is
        carried here rather than with the sample because it is an instrument
        misalignment, and treating it as one is what stops a refinement from
        absorbing it into the lattice parameter.
    name : str
        A human name for the configuration, so a pattern can say which
        instrument setting it was modelled with.
    """

    caglioti_u: float = 0.0
    caglioti_v: float = 0.0
    caglioti_w: float = 0.01
    lorentzian_x: float = 0.0
    lorentzian_y: float = 0.0
    zero_shift_deg: float = 0.0
    name: str = "instrument"

    def __post_init__(self) -> None:
        for field_name in (
            "caglioti_u",
            "caglioti_v",
            "caglioti_w",
            "lorentzian_x",
            "lorentzian_y",
            "zero_shift_deg",
        ):
            value = float(getattr(self, field_name))
            if not np.isfinite(value):
                raise ValueError(f"InstrumentBroadening.{field_name} must be finite.")
            object.__setattr__(self, field_name, value)
        if self.lorentzian_x < 0.0 or self.lorentzian_y < 0.0:
            raise ValueError("Lorentzian resolution coefficients must be non-negative.")
        if not self.name.strip():
            raise ValueError("InstrumentBroadening.name must be non-empty.")

    @classmethod
    def laboratory_bragg_brentano(cls) -> InstrumentBroadening:
        """Return a representative sealed-tube Bragg-Brentano resolution function.

        These are *plausible* coefficients for a well-aligned laboratory
        powder diffractometer with medium slits, giving about 0.08 degrees FWHM
        at 30 degrees and about 0.16 degrees at 120 degrees. They are a starting
        point for teaching and for method development, and are not a substitute
        for calibrating an actual instrument with
        :func:`calibrate_instrument_broadening`.
        """

        return cls(
            caglioti_u=0.0035,
            caglioti_v=-0.0021,
            caglioti_w=0.0060,
            lorentzian_x=0.0,
            lorentzian_y=0.012,
            name="Laboratory Bragg-Brentano (representative)",
        )

    @classmethod
    def ideal(cls, fwhm_deg: float = 0.1) -> InstrumentBroadening:
        """Return a constant-width resolution function of ``fwhm_deg``.

        Angle-independent broadening is physically unrealistic, but it is the
        right model for a synthetic teaching pattern where the point is the
        structure and not the diffractometer.
        """

        if fwhm_deg <= 0.0:
            raise ValueError("fwhm_deg must be strictly positive.")
        return cls(caglioti_w=float(fwhm_deg) ** 2, name=f"Constant {fwhm_deg:g} deg FWHM")

    def gaussian_fwhm_deg(self, two_theta_deg: np.ndarray | float) -> np.ndarray:
        """Return the Caglioti Gaussian FWHM in degrees at each ``2*theta``."""

        tangent = np.tan(np.deg2rad(0.5 * np.asarray(two_theta_deg, dtype=np.float64)))
        squared = self.caglioti_u * tangent**2 + self.caglioti_v * tangent + self.caglioti_w
        if np.any(squared <= 0.0):
            raise ValueError(
                "The Caglioti coefficients give a non-positive squared width in this angular "
                "range; U, V and W must keep U tan^2(theta) + V tan(theta) + W positive."
            )
        return np.sqrt(squared)

    def lorentzian_fwhm_deg(self, two_theta_deg: np.ndarray | float) -> np.ndarray:
        """Return the Lorentzian FWHM in degrees at each ``2*theta``."""

        theta = np.deg2rad(0.5 * np.asarray(two_theta_deg, dtype=np.float64))
        return self.lorentzian_x * np.tan(theta) + self.lorentzian_y / np.cos(theta)

    def fwhm_deg(self, two_theta_deg: np.ndarray | float) -> np.ndarray:
        """Return the combined pseudo-Voigt FWHM in degrees at each ``2*theta``.

        The Gaussian and Lorentzian halves are combined by the
        Thompson-Cox-Hastings fifth-order relation rather than added, because
        the convolution of a Gaussian with a Lorentzian is neither.
        """

        gaussian = self.gaussian_fwhm_deg(two_theta_deg)
        lorentzian = self.lorentzian_fwhm_deg(two_theta_deg)
        return _tch_width(gaussian, lorentzian)

    def eta(self, two_theta_deg: np.ndarray | float) -> np.ndarray:
        """Return the pseudo-Voigt Lorentzian fraction at each ``2*theta``.

        Zero is pure Gaussian and one pure Lorentzian. It varies with angle
        because the two halves of the resolution function have different
        angular dependences, which is exactly why a single fixed ``eta`` cannot
        describe a real instrument across a wide scan.
        """

        gaussian = self.gaussian_fwhm_deg(two_theta_deg)
        lorentzian = self.lorentzian_fwhm_deg(two_theta_deg)
        return _tch_eta(gaussian, lorentzian)

    def describe(self) -> str:
        """Return the coefficients, the widths they imply, and their standing."""

        probe = np.array([20.0, 60.0, 120.0])
        widths = self.fwhm_deg(probe)
        mixings = self.eta(probe)
        quoted = ", ".join(
            f"{angle:g} deg: {width:.4f} deg FWHM (eta = {mixing:.3f})"
            for angle, width, mixing in zip(probe, widths, mixings, strict=True)
        )
        return (
            f"Instrument resolution function '{self.name}'. Gaussian FWHM^2 = "
            f"{self.caglioti_u:.6g} tan^2(theta) + {self.caglioti_v:.6g} tan(theta) + "
            f"{self.caglioti_w:.6g} deg^2 ({_CITATION_CAGLIOTI}); Lorentzian FWHM = "
            f"{self.lorentzian_x:.6g} tan(theta) + {self.lorentzian_y:.6g} / cos(theta) deg, "
            f"combined as a pseudo-Voigt by {_CITATION_TCH} Zero shift = "
            f"{self.zero_shift_deg:.6g} deg 2*theta. Implied widths -- {quoted}. This is the width "
            "the diffractometer contributes; any width measured beyond it belongs to the specimen."
        )


@dataclass(frozen=True, slots=True)
class WilliamsonHallAnalysis:
    """Crystallite size and microstrain from the angular trend of peak widths.

    Purpose
    -------
    A single peak width cannot separate size from strain: both widen the peak.
    Their *angular dependences* differ -- size broadening scales as
    ``1 / cos(theta)`` and strain broadening as ``tan(theta)`` -- so the trend
    across several reflections separates them where one peak cannot.

    Attributes
    ----------
    two_theta_deg : np.ndarray
        The reflections used, in degrees.
    sample_fwhm_deg : np.ndarray
        Sample-only widths, after instrumental deconvolution.
    abscissa, ordinate : np.ndarray
        The plotted quantities ``4 sin(theta)`` and ``beta cos(theta)``, in
        radians, retained so the fit can be replotted and judged by eye.
    crystallite_size_nm : float
        Volume-averaged column length from the intercept.
    microstrain : float
        Dimensionless strain ``epsilon`` from the slope; multiply by 100 for
        percent. Reported as ``4 epsilon`` in some sources, so the convention is
        stated in :meth:`describe`.
    intercept, slope : float
    r_squared : float
        Coefficient of determination of the straight-line fit. A low value says
        the uniform-deformation model does not fit, most often because the
        broadening is anisotropic; the size and strain are then not meaningful.
    shape_factor : float
    wavelength_angstrom : float
    """

    two_theta_deg: np.ndarray
    sample_fwhm_deg: np.ndarray
    abscissa: np.ndarray
    ordinate: np.ndarray
    crystallite_size_nm: float
    microstrain: float
    intercept: float
    slope: float
    r_squared: float
    shape_factor: float
    wavelength_angstrom: float

    def __post_init__(self) -> None:
        arrays = tuple(
            as_float_array(value, shape=(None,))
            for value in (
                self.two_theta_deg,
                self.sample_fwhm_deg,
                self.abscissa,
                self.ordinate,
            )
        )
        if arrays[0].size < 2 or len({array.shape for array in arrays}) != 1:
            raise ValueError(
                "A Williamson-Hall analysis needs at least two reflections and aligned arrays."
            )
        if np.any(~np.isfinite(np.concatenate(arrays))):
            raise ValueError("Williamson-Hall arrays must be finite.")
        if np.any(arrays[1] <= 0.0):
            raise ValueError("Sample-only widths must be strictly positive.")
        for name in ("crystallite_size_nm", "microstrain", "intercept", "slope", "r_squared"):
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f"WilliamsonHallAnalysis.{name} must be finite.")
        if self.crystallite_size_nm <= 0.0:
            raise ValueError("WilliamsonHallAnalysis.crystallite_size_nm must be positive.")
        if self.shape_factor <= 0.0 or self.wavelength_angstrom <= 0.0:
            raise ValueError("Shape factor and wavelength must be strictly positive.")
        for index, name in enumerate(("two_theta_deg", "sample_fwhm_deg", "abscissa", "ordinate")):
            object.__setattr__(self, name, arrays[index])

    @property
    def reflection_count(self) -> int:
        """Return the number of reflections the straight line was fitted to."""

        return int(self.two_theta_deg.size)

    def describe(self) -> str:
        """Return the model, the numbers, the fit quality, and the caveats."""

        strain_note = (
            "The strain is negative, which the uniform-deformation model cannot mean physically; "
            "it usually indicates that the instrumental width has been over-subtracted or that "
            "size broadening dominates and the slope is fitting noise. "
            if self.microstrain < 0.0
            else ""
        )
        quality = (
            "The straight-line fit is poor (R^2 below 0.8), so the uniform-deformation assumption "
            "is not supported by these reflections and the separation into size and strain should "
            "not be quoted. "
            if self.r_squared < 0.8
            else ""
        )
        return (
            f"Williamson-Hall analysis of {self.reflection_count} reflections between "
            f"{self.two_theta_deg[0]:.3f} and {self.two_theta_deg[-1]:.3f} degrees 2*theta at "
            f"lambda = {self.wavelength_angstrom:.6g} A ({_CITATION_WILLIAMSON_HALL}). Fitting "
            "beta cos(theta) = K lambda / D + 4 epsilon sin(theta) with beta the sample-only FWHM "
            f"in radians gives an intercept of {self.intercept:.6g} and a slope of "
            f"{self.slope:.6g}, hence a volume-averaged crystallite size D = "
            f"{self.crystallite_size_nm:.4g} nm with shape factor K = {self.shape_factor:g} "
            f"({_CITATION_LANGFORD_WILSON}) and a microstrain epsilon = {self.microstrain:.6g} "
            f"({100.0 * self.microstrain:.4g}%). R^2 = {self.r_squared:.6g}. {quality}{strain_note}"
            "The model assumes isotropic size and strain broadening and a common peak shape; "
            "anisotropic broadening, stacking faults and a compositional gradient all violate it "
            "and show up as scatter about the line rather than as an error bar."
        )


def _tch_width(gaussian: np.ndarray, lorentzian: np.ndarray) -> np.ndarray:
    """Thompson-Cox-Hastings pseudo-Voigt FWHM from Gaussian and Lorentzian FWHM."""

    g = np.asarray(gaussian, dtype=np.float64)
    l_ = np.asarray(lorentzian, dtype=np.float64)
    combined = (
        g**5
        + _TCH_A * g**4 * l_
        + _TCH_B * g**3 * l_**2
        + _TCH_C * g**2 * l_**3
        + _TCH_D * g * l_**4
        + l_**5
    )
    return np.power(combined, 0.2)


def _tch_eta(gaussian: np.ndarray, lorentzian: np.ndarray) -> np.ndarray:
    """Thompson-Cox-Hastings Lorentzian fraction of the combined pseudo-Voigt."""

    width = _tch_width(gaussian, lorentzian)
    ratio = np.divide(
        np.asarray(lorentzian, dtype=np.float64),
        width,
        out=np.zeros_like(width),
        where=width > 0.0,
    )
    c1, c2, c3 = _TCH_ETA
    mixing = np.clip(c1 * ratio + c2 * ratio**2 + c3 * ratio**3, 0.0, 1.0)
    return np.asarray(mixing, dtype=np.float64)


def calibrate_instrument_broadening(
    two_theta_deg: np.ndarray,
    fwhm_deg: np.ndarray,
    *,
    name: str = "calibrated",
    zero_shift_deg: float = 0.0,
) -> InstrumentBroadening:
    """Fit Caglioti ``U``, ``V``, ``W`` to the peak widths of a line-profile standard.

    Purpose
    -------
    Turn a measurement of a standard -- NIST SRM 660 lanthanum hexaboride, SRM
    640 silicon, or any specimen whose own broadening is negligible against the
    instrument's -- into a resolution function that can be subtracted from the
    widths of a real sample.

    Method
    ------
    ``FWHM^2 = U tan^2(theta) + V tan(theta) + W`` is *linear* in ``U``, ``V``
    and ``W``, so the fit is an exact linear least squares in the squared
    widths, with no starting values and no convergence question. The result is
    purely Gaussian: a Lorentzian component cannot be recovered from widths
    alone and must come from a shape fit.

    Parameters
    ----------
    two_theta_deg
        Peak positions of the standard, in degrees. At least three, and they
        must span a useful angular range -- three peaks within ten degrees of
        each other determine the parabola no better than one does.
    fwhm_deg
        Measured full widths at half maximum of those peaks, in degrees.
    name
        A label for the resulting configuration.
    zero_shift_deg
        Detector zero error to record alongside the widths, if it was
        determined from the same standard.

    Returns
    -------
    InstrumentBroadening
        A Gaussian-only resolution function fitted to the standard.

    Raises
    ------
    ValueError
        If fewer than three peaks are given, if the arrays disagree, if any
        width is non-positive, or if the fitted coefficients would give a
        non-positive width somewhere inside the calibrated range -- which means
        the fit has extrapolated itself into nonsense and must not be used.

    Examples
    --------
    >>> import numpy as np
    >>> from pytex.diffraction.xrd_instrument import (
    ...     InstrumentBroadening, calibrate_instrument_broadening
    ... )
    >>> truth = InstrumentBroadening(caglioti_u=0.004, caglioti_v=-0.002, caglioti_w=0.006)
    >>> angles = np.array([21.4, 30.4, 37.5, 43.6, 53.9, 67.6, 79.3, 95.2, 115.4])
    >>> widths = truth.gaussian_fwhm_deg(angles)
    >>> fitted = calibrate_instrument_broadening(angles, widths)
    >>> bool(np.isclose(fitted.caglioti_u, 0.004, atol=1e-9))
    True
    """

    angles = as_float_array(two_theta_deg, shape=(None,))
    widths = as_float_array(fwhm_deg, shape=(None,))
    if angles.shape != widths.shape:
        raise ValueError("Standard peak positions and widths must align.")
    if angles.size < 3:
        raise ValueError(
            "Fitting U, V and W needs at least three standard peaks; with fewer the "
            "resolution function is under-determined."
        )
    if np.any(~np.isfinite(angles)) or np.any(angles <= 0.0) or np.any(angles >= 180.0):
        raise ValueError("Standard peak positions must be finite and inside (0, 180) degrees.")
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("Standard peak widths must be finite and strictly positive.")
    tangent = np.tan(np.deg2rad(0.5 * angles))
    design = np.column_stack((tangent**2, tangent, np.ones_like(tangent)))
    coefficients, *_ = np.linalg.lstsq(design, np.square(widths), rcond=None)
    u, v, w = (float(value) for value in coefficients)
    candidate = InstrumentBroadening(
        caglioti_u=u,
        caglioti_v=v,
        caglioti_w=w,
        zero_shift_deg=float(zero_shift_deg),
        name=name,
    )
    # A parabola through three or four widths can dip below zero between them,
    # and the failure is silent -- the object constructs, and only a later width
    # evaluation raises. Check the calibrated range now, where the fix (more
    # standard peaks) is still obvious.
    candidate.gaussian_fwhm_deg(np.linspace(float(angles[0]), float(angles[-1]), 128))
    return candidate


def deconvolve_instrument_width(
    observed_fwhm_deg: np.ndarray | float,
    instrument: InstrumentBroadening,
    two_theta_deg: np.ndarray | float,
    *,
    mode: DeconvolutionMode = "gaussian",
) -> np.ndarray:
    """Remove the instrumental width from measured widths, leaving the sample's.

    Purpose
    -------
    Answer "how wide would this peak be on a perfect instrument?", which is the
    only width from which a crystallite size or a microstrain may be quoted.

    Method
    ------
    Deconvolution of widths depends on peak shape, and the three available modes
    are the three standard answers:

    - ``"gaussian"``: ``beta^2 = observed^2 - instrument^2``. Correct when both
      profiles are Gaussian; the usual choice for a well-behaved laboratory
      instrument.
    - ``"lorentzian"``: ``beta = observed - instrument``. Correct when both are
      Lorentzian; appropriate when size broadening dominates, since size
      broadening is Lorentzian-like.
    - ``"pseudo_voigt"``: the Gaussian and Lorentzian components are separated
      by the Thompson-Cox-Hastings relation and deconvolved in their own way
      before being recombined. The most defensible when the shapes are mixed,
      and the most sensitive to a poor instrument calibration.

    Parameters
    ----------
    observed_fwhm_deg
        Measured widths, in degrees.
    instrument
        The calibrated resolution function.
    two_theta_deg
        Angles at which the widths were measured, in degrees.
    mode
        Deconvolution model, as above.

    Returns
    -------
    np.ndarray
        Sample-only widths in degrees.

    Raises
    ------
    ValueError
        If any measured width is narrower than the instrument's own width at
        that angle. That is not a very small crystallite: it is an impossible
        measurement, and it means the calibration does not apply to this data.
    """

    if mode not in DECONVOLUTION_MODES:
        raise ValueError(f"mode must be one of {DECONVOLUTION_MODES}.")
    observed = np.atleast_1d(np.asarray(observed_fwhm_deg, dtype=np.float64))
    angles = np.atleast_1d(np.asarray(two_theta_deg, dtype=np.float64))
    if observed.shape != angles.shape:
        raise ValueError("Observed widths and their angles must align.")
    if np.any(~np.isfinite(observed)) or np.any(observed <= 0.0):
        raise ValueError("Observed widths must be finite and strictly positive.")
    total = instrument.fwhm_deg(angles)
    if np.any(observed <= total):
        narrowest = int(np.argmin(observed - total))
        raise ValueError(
            "A measured width is not wider than the instrumental width at "
            f"{angles[narrowest]:.4f} degrees ({observed[narrowest]:.5f} deg measured against "
            f"{total[narrowest]:.5f} deg instrumental). The specimen cannot be sharper than the "
            "instrument, so either the calibration belongs to a different configuration or the "
            "widths were measured differently from the standard's."
        )
    if mode == "gaussian":
        sample = np.sqrt(observed**2 - total**2)
    elif mode == "lorentzian":
        sample = observed - total
    else:
        instrument_gaussian = instrument.gaussian_fwhm_deg(angles)
        instrument_lorentzian = instrument.lorentzian_fwhm_deg(angles)
        eta = _tch_eta(instrument_gaussian, instrument_lorentzian)
        # Split the observed width into components with the instrument's own
        # mixing, deconvolve each in its own algebra, then recombine. Using the
        # instrument's eta is an assumption -- the sample may be mixed
        # differently -- and it is the assumption this mode is named for.
        observed_lorentzian = eta * observed
        observed_gaussian = np.sqrt(np.clip(observed**2 - observed_lorentzian**2, 0.0, None))
        sample_gaussian = np.sqrt(np.clip(observed_gaussian**2 - instrument_gaussian**2, 0.0, None))
        sample_lorentzian = np.clip(observed_lorentzian - instrument_lorentzian, 0.0, None)
        sample = _tch_width(sample_gaussian, sample_lorentzian)
    sample = np.ascontiguousarray(sample, dtype=np.float64)
    sample.setflags(write=False)
    return sample


def scherrer_size_nm(
    fwhm_deg: np.ndarray | float,
    two_theta_deg: np.ndarray | float,
    *,
    wavelength_angstrom: float,
    shape_factor: float = SCHERRER_SHAPE_FACTOR,
) -> np.ndarray:
    """Return the Scherrer crystallite size in nanometres from sample-only widths.

    ``D = K lambda / (beta cos(theta))`` with ``beta`` the sample-only FWHM in
    radians. The width must already have had the instrumental contribution
    removed by :func:`deconvolve_instrument_width`; passing a raw measured width
    yields the size of the diffractometer's resolution, not of the crystallites.

    Scherrer's relation attributes *all* remaining broadening to size, so it is
    an upper bound on broadening and therefore a **lower bound on size** whenever
    strain is present. Use :func:`williamson_hall` when several reflections are
    available.
    """

    if wavelength_angstrom <= 0.0:
        raise ValueError("wavelength_angstrom must be strictly positive.")
    if shape_factor <= 0.0:
        raise ValueError("shape_factor must be strictly positive.")
    widths = np.atleast_1d(np.asarray(fwhm_deg, dtype=np.float64))
    angles = np.atleast_1d(np.asarray(two_theta_deg, dtype=np.float64))
    if widths.shape != angles.shape:
        raise ValueError("Widths and their angles must align.")
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("Widths must be finite and strictly positive.")
    beta = np.deg2rad(widths)
    cos_theta = np.cos(np.deg2rad(0.5 * angles))
    size_angstrom = shape_factor * wavelength_angstrom / (beta * cos_theta)
    size_nm = np.ascontiguousarray(size_angstrom / 10.0, dtype=np.float64)
    size_nm.setflags(write=False)
    return size_nm


def williamson_hall(
    two_theta_deg: np.ndarray,
    sample_fwhm_deg: np.ndarray,
    *,
    wavelength_angstrom: float,
    shape_factor: float = SCHERRER_SHAPE_FACTOR,
) -> WilliamsonHallAnalysis:
    """Separate crystallite size from microstrain across several reflections.

    Purpose
    -------
    Size broadening and strain broadening cannot be told apart in one peak. They
    can be told apart across peaks, because they depend on angle differently:
    fitting ``beta cos(theta) = K lambda / D + 4 epsilon sin(theta)`` puts the
    size in the intercept and the strain in the slope.

    Parameters
    ----------
    two_theta_deg
        Reflection positions in degrees. At least two, and in practice at least
        four spanning a wide angular range, because the intercept of a line
        fitted over a narrow range is almost unconstrained.
    sample_fwhm_deg
        Sample-only widths in degrees, already deconvolved from the instrument
        by :func:`deconvolve_instrument_width`.
    wavelength_angstrom
        The radiation wavelength the widths were measured at.
    shape_factor
        Scherrer constant ``K``; 0.9 by convention for spherical crystallites.

    Returns
    -------
    WilliamsonHallAnalysis
        Size, strain, the fitted line, its ``R^2``, and an explainable summary
        that states when the model does not hold.

    Raises
    ------
    ValueError
        If fewer than two reflections are given, if the arrays disagree, or if
        the fitted intercept is non-positive -- which corresponds to an
        infinite or negative crystallite size and means the data do not support
        this model rather than that the crystallites are very large.
    """

    angles = as_float_array(two_theta_deg, shape=(None,))
    widths = as_float_array(sample_fwhm_deg, shape=(None,))
    if angles.shape != widths.shape:
        raise ValueError("Reflection angles and sample widths must align.")
    if angles.size < 2:
        raise ValueError("A Williamson-Hall fit needs at least two reflections.")
    if wavelength_angstrom <= 0.0:
        raise ValueError("wavelength_angstrom must be strictly positive.")
    if shape_factor <= 0.0:
        raise ValueError("shape_factor must be strictly positive.")
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("Sample-only widths must be finite and strictly positive.")
    order = np.argsort(angles)
    angles = angles[order]
    widths = widths[order]
    theta = np.deg2rad(0.5 * angles)
    abscissa = 4.0 * np.sin(theta)
    ordinate = np.deg2rad(widths) * np.cos(theta)
    slope, intercept = (float(value) for value in np.polyfit(abscissa, ordinate, 1))
    if intercept <= 0.0:
        raise ValueError(
            "The Williamson-Hall intercept is not positive, so no finite crystallite size follows "
            "from these widths. This is a statement about the fit, not about the specimen: with "
            "few reflections or a narrow angular range the intercept is poorly constrained."
        )
    predicted = slope * abscissa + intercept
    total_variance = float(np.sum((ordinate - np.mean(ordinate)) ** 2))
    residual_variance = float(np.sum((ordinate - predicted) ** 2))
    r_squared = 1.0 if total_variance <= 0.0 else 1.0 - residual_variance / total_variance
    size_nm = shape_factor * wavelength_angstrom / intercept / 10.0
    return WilliamsonHallAnalysis(
        two_theta_deg=angles,
        sample_fwhm_deg=widths,
        abscissa=abscissa,
        ordinate=ordinate,
        crystallite_size_nm=float(size_nm),
        # The abscissa is already 4 sin(theta), so the fitted slope is epsilon
        # itself rather than 4 epsilon.
        microstrain=float(slope),
        intercept=intercept,
        slope=slope,
        r_squared=float(r_squared),
        shape_factor=float(shape_factor),
        wavelength_angstrom=float(wavelength_angstrom),
    )


__all__ = [
    "DECONVOLUTION_MODES",
    "INSTRUMENT_BROADENING_SCHEMA",
    "SCHERRER_SHAPE_FACTOR",
    "WILLIAMSON_HALL_SCHEMA",
    "DeconvolutionMode",
    "InstrumentBroadening",
    "WilliamsonHallAnalysis",
    "calibrate_instrument_broadening",
    "deconvolve_instrument_width",
    "scherrer_size_nm",
    "williamson_hall",
]
