"""Day-to-day corrections and display transforms for measured powder patterns.

A laboratory diffractogram is not a measurement of a crystal. It is a
measurement of a crystal *plus* a goniometer whose zero is not quite right, a
specimen surface that is not quite on the axis, a beam that penetrates before
it diffracts, a doublet source, and an optical train with its own polarization
and slit behaviour. This module holds the corrections that separate the two,
and keeps them as small, named, individually testable operations rather than a
single opaque "prepare data" step -- because a corrected pattern must always be
able to say what was done to it.

The corrections divide into two kinds, and confusing them is the classic error:

**Position aberrations** move peaks and therefore change every lattice
parameter derived from them. The three that matter on a Bragg-Brentano
diffractometer have *different* angular signatures, which is the only reason
they can be separated at all:

===========================  ==============================  ===================
Aberration                   ``Delta(2 theta)``              Typical magnitude
===========================  ==============================  ===================
Detector zero                constant                        0.01 to 0.05 deg
Specimen displacement ``s``  ``-2 s cos(theta) / R``         50 um gives 0.03 deg
Specimen transparency        ``-sin(2 theta) / (2 mu R)``    large for low ``mu``
Refraction                   ``2 delta tan(theta)``          about 1 part in 1e5
===========================  ==============================  ===================

**Intensity corrections** change relative peak areas and therefore quantitative
and structural work, but leave positions alone: monochromator polarization,
variable-to-fixed divergence slit conversion.

Separately, this module holds the **display** transforms -- abscissa in
``2 theta``, ``d``, ``Q`` or ``sin^2(theta)``, ordinate on a linear, square-root
or logarithmic scale -- because reading a pattern is part of analysing it, and
a square-root ordinate is the honest way to show weak reflections next to
strong ones.

One warning is repeated in the code because it is repeatedly ignored:
:func:`strip_kalpha2` is for looking and for picking candidates. Every *fitting*
path in PyTex models the K-alpha2 line instead (see
:mod:`pytex.diffraction.xrd_peaks`), because stripping assumes an exact
intensity ratio and an identical profile shape, and it propagates each
subtraction into the next, so the noise of the stripped pattern grows with
angle.

References
----------
Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Ch. 11 and Appendix -- the diffractometer aberrations and
their angular dependence.

Rachinger, W. A., *J. Sci. Instrum.* **25** (1948) 254-255,
doi:10.1088/0950-7671/25/7/125 -- the K-alpha2 stripping algorithm.

Wilson, A. J. C., *Mathematical Theory of X-ray Powder Diffractometry*, Philips
Technical Library (1963) -- specimen displacement, transparency and flat
specimen aberrations derived from the geometry.

Savitzky, A. & Golay, M. J. E., *Anal. Chem.* **36** (1964) 1627-1639,
doi:10.1021/ac60214a047 -- the smoothing filter, and the width bias it causes.

Klug, H. P. & Alexander, L. E., *X-ray Diffraction Procedures*, 2nd ed., Wiley
(1974), Chs. 3 and 5 -- monochromator polarization and divergence-slit geometry.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
from scipy.signal import savgol_filter

from pytex.core._arrays import as_float_array
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import PeakTable

PROFILE_VIEW_SCHEMA = "pytex.diffraction.powder_profile_view"

ProfileAbscissa = Literal["two_theta_deg", "d_angstrom", "q_inv_angstrom", "sin_squared_theta"]
IntensityScale = Literal["linear", "sqrt", "log10"]
IntensityNormalization = Literal["none", "maximum", "integral"]

PROFILE_ABSCISSAE: tuple[ProfileAbscissa, ...] = (
    "two_theta_deg",
    "d_angstrom",
    "q_inv_angstrom",
    "sin_squared_theta",
)
INTENSITY_SCALES: tuple[IntensityScale, ...] = ("linear", "sqrt", "log10")
INTENSITY_NORMALIZATIONS: tuple[IntensityNormalization, ...] = ("none", "maximum", "integral")

#: Classical electron radius in angstrom, for the refraction decrement.
_ELECTRON_RADIUS_ANGSTROM = 2.8179403262e-5

#: Avogadro's number per mole.
_AVOGADRO = 6.02214076e23

_CITATION_RACHINGER = (
    "Rachinger, J. Sci. Instrum. 25 (1948) 254, doi:10.1088/0950-7671/25/7/125."
)


def _record(pattern: MeasuredPowderPattern, step: str) -> dict[str, str]:
    """Return ``pattern``'s metadata with ``step`` appended to its audit trail.

    Corrections compose, and the order they were applied in changes the result,
    so the record is a list rather than a set of flags.
    """

    metadata = dict(pattern.metadata)
    existing = metadata.get("corrections", "")
    metadata["corrections"] = f"{existing}; {step}" if existing else step
    return metadata


# ---------------------------------------------------------------------------
# K-alpha2 stripping
# ---------------------------------------------------------------------------


def kalpha2_partner_two_theta_deg(two_theta_deg: Any, *, wavelength_ratio: float) -> np.ndarray:
    """Return where the K-alpha2 line of each K-alpha1 reflection falls.

    Purpose
    -------
    Give every routine that reasons about the doublet -- stripping, detection,
    fitting -- one place where the geometry is written down.

    Method
    ------
    Both lines diffract from the same planes, so ``d`` is common and Bragg's law
    gives ``sin(theta_2) = (lambda_2 / lambda_1) sin(theta_1)``. The separation
    is therefore *not* a constant offset: differentiating at fixed ``d`` gives
    ``Delta(2 theta) = 2 (Delta lambda / lambda) tan(theta)``, so the pair is
    unresolved at low angle and cleanly split above roughly 90 degrees.

    Parameters
    ----------
    two_theta_deg
        K-alpha1 positions in degrees.
    wavelength_ratio
        ``lambda_2 / lambda_1``, slightly greater than one for a K-alpha
        doublet.

    Returns
    -------
    np.ndarray
        K-alpha2 positions in degrees, same shape as the input.

    Raises
    ------
    ValueError
        If the ratio is not positive, or a partner would fall beyond the Ewald
        limit (``sin(theta) > 1``), which means the input angle is too high for
        this wavelength pair.
    """

    if not np.isfinite(wavelength_ratio) or wavelength_ratio <= 0.0:
        raise ValueError("wavelength_ratio must be finite and positive.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    argument = wavelength_ratio * np.sin(np.deg2rad(0.5 * angles))
    if np.any(np.abs(argument) > 1.0):
        raise ValueError(
            "A K-alpha2 partner falls beyond the Ewald limit: the K-alpha1 angle is too high "
            "for this wavelength pair."
        )
    partner: np.ndarray = np.rad2deg(2.0 * np.arcsin(argument))
    return partner


def strip_kalpha2(
    measured: MeasuredPowderPattern,
    *,
    radiation: RadiationSpec | None = None,
    intensity_ratio: float | None = None,
) -> MeasuredPowderPattern:
    """Remove the K-alpha2 contribution by Rachinger subtraction (display use).

    Purpose
    -------
    Produce the single-wavelength pattern that a human reads, and that a
    simple peak-picker can work on, from a doublet-source measurement.

    **Do not use this ahead of profile fitting or lattice-parameter work.** The
    library's fitting paths model the doublet instead
    (:func:`pytex.diffraction.xrd_peaks.fit_peaks`), which makes the same
    physical assumptions without touching the data. This function exists
    because a stripped pattern is genuinely easier to look at, and because
    comparing the two is instructive -- not because stripping is the better
    analysis.

    Method
    ------
    Rachinger's (1948) recursion. Assume the alpha2 profile is the alpha1
    profile scaled by ``r = I_2 / I_1`` and displaced to the angle Bragg's law
    puts it at. Then, sweeping from low angle upward,

    ``I_1(2 theta) = I(2 theta) - r I_1(2 theta')``

    where ``2 theta'`` is the alpha1 angle whose alpha2 partner lands on
    ``2 theta``, namely ``sin(theta') = (lambda_1 / lambda_2) sin(theta)``.
    Because ``lambda_1 < lambda_2``, ``2 theta'`` always lies *below*
    ``2 theta``, so the recursion only ever needs values it has already
    computed. ``I_1`` is interpolated linearly at ``2 theta'``.

    The three assumptions are exactly the three weaknesses: the ratio is taken
    as exact, the two profiles are taken as identically shaped, and each
    subtraction feeds the next. That last point is quantitative -- the variance
    propagates as
    ``var(I_1(2 theta)) = var(I(2 theta)) + r^2 var(I_1(2 theta'))`` -- so the
    stripped pattern is noisier than the measurement, increasingly so with
    angle. When the input carries standard uncertainties, the output carries
    the propagated ones, so this cost is visible rather than assumed away.

    Parameters
    ----------
    measured
        The measured doublet pattern.
    radiation
        Radiation declaring both wavelengths. Falls back to
        ``measured.radiation``.
    intensity_ratio
        ``I_2 / I_1``, overriding the radiation's tabulated value. The
        conventional 0.5 follows from the statistical weights of the
        ``2p_{3/2}`` and ``2p_{1/2}`` initial states.

    Returns
    -------
    MeasuredPowderPattern
        The stripped pattern, with radiation replaced by the single K-alpha1
        line (because the result no longer contains a doublet), propagated
        uncertainties where available, and the operation recorded in metadata.

    Raises
    ------
    ValueError
        If no radiation is available, the radiation declares no K-alpha2 line,
        or the ratio is outside ``(0, 1]``.

    See Also
    --------
    pytex.diffraction.xrd_peaks.fit_peaks : models the doublet instead, and is
        what quantitative work should use.

    Notes
    -----
    Negative excursions are clipped to zero. They arise where the assumed ratio
    or shape is slightly wrong, and a negative intensity is not a physical
    result; the clip is stated here rather than hidden because it means the
    stripped pattern's background is very slightly biased upward.

    **Ringing is expected and is not a bug.** Because the recursion subtracts a
    value it computed earlier, whatever it fails to remove at one angle is
    re-subtracted at the alpha2 partner of *that* angle, and again beyond it.
    The result is a decaying train of small alternating residuals above each
    strong reflection, spaced by the growing doublet separation. On a synthetic
    Ni pattern this method removes 86 to 94 per cent of each alpha2 line, and
    what survives above the strongest reflection is still large enough for a
    peak detector to report. The residual has two causes that no choice of
    ratio fixes: the true alpha2 profile is slightly *wider* in ``2 theta``
    than the alpha1 profile, because ``d(2 theta_2) / d(2 theta_1)`` is not
    one, and the tabulated intensity ratio is not exactly one half. The
    lowest few points of the scan are also unreliable, since the recursion has
    no measured values below the first point to draw on.

    All of this is the argument for modelling rather than stripping, stated
    quantitatively. A fit that includes the alpha2 line makes the same two
    assumptions, but makes them *inside* the model where they are visible in
    the residual, instead of writing their failure into the data.
    """

    spec = radiation if radiation is not None else measured.radiation
    if spec is None:
        raise ValueError(
            "strip_kalpha2 needs a radiation: the pattern declares none and none was passed."
        )
    if spec.kalpha2_wavelength_angstrom is None:
        raise ValueError(
            f"Radiation '{spec.name}' declares no K-alpha2 line, so there is nothing to strip. "
            "Use a doublet specification such as RadiationSpec.cu_ka_doublet()."
        )
    ratio = (
        float(spec.kalpha2_relative_intensity)
        if intensity_ratio is None
        else float(intensity_ratio)
    )
    if not np.isfinite(ratio) or not 0.0 < ratio <= 1.0:
        raise ValueError("strip_kalpha2 requires an intensity ratio in (0, 1].")

    axis = np.asarray(measured.two_theta_deg, dtype=np.float64)
    observed = np.asarray(measured.intensity, dtype=np.float64)
    inverse_ratio = float(spec.wavelength_angstrom / spec.kalpha2_wavelength_angstrom)

    # The alpha1 angle whose alpha2 partner lands on each measured angle.
    argument = inverse_ratio * np.sin(np.deg2rad(0.5 * axis))
    source = np.rad2deg(2.0 * np.arcsin(np.clip(argument, -1.0, 1.0)))

    stripped = np.empty_like(observed)
    variance = np.zeros_like(observed)
    observed_variance = (
        np.square(np.asarray(measured.standard_uncertainty, dtype=np.float64))
        if measured.standard_uncertainty is not None
        else None
    )
    # Sweep upward: every value the recursion needs lies at a lower angle and
    # has already been written.
    for index in range(axis.size):
        partner_angle = float(source[index])
        if partner_angle <= axis[0]:
            contribution = 0.0
            contribution_variance = 0.0
        else:
            contribution = float(np.interp(partner_angle, axis[:index], stripped[:index]))
            contribution_variance = float(np.interp(partner_angle, axis[:index], variance[:index]))
        stripped[index] = observed[index] - ratio * contribution
        variance[index] = (
            0.0 if observed_variance is None else float(observed_variance[index])
        ) + ratio * ratio * contribution_variance
    stripped = np.maximum(stripped, 0.0)
    if not np.any(stripped > 0.0):
        raise ValueError(
            "Rachinger stripping removed all intensity, which means the assumed ratio or the "
            "wavelength pair does not match this pattern."
        )

    single = RadiationSpec(
        name=f"{spec.name} (K-alpha1, alpha2 stripped)",
        wavelength_angstrom=spec.wavelength_angstrom,
        anode=spec.anode,
        kind=spec.kind,
    )
    return replace(
        measured,
        name=f"{measured.name} (Ka2 stripped)",
        intensity=stripped,
        standard_uncertainty=(
            None if observed_variance is None else np.sqrt(np.maximum(variance, 1e-12))
        ),
        radiation=single,
        metadata=_record(measured, f"rachinger_kalpha2_strip(ratio={ratio:.4f})"),
    )


# ---------------------------------------------------------------------------
# Position aberrations
# ---------------------------------------------------------------------------


def zero_shift_deg(two_theta_deg: Any, *, zero_deg: float) -> np.ndarray:
    """Return the constant detector zero-point error at each angle.

    Purpose
    -------
    Make the simplest aberration explicit, so that it can be told apart from
    the others by its angular signature rather than by assertion.

    Method
    ------
    A misaligned ``2 theta`` origin adds the same amount everywhere:
    ``Delta(2 theta) = zero``. Being constant is what distinguishes it from
    specimen displacement (``cos(theta)``) and transparency (``sin(2 theta)``),
    and what makes it, in principle, separable from them. In practice the
    separation is poor over a single specimen's angular range, which is why
    zero belongs to a calibrated instrument and displacement to the specimen.

    Parameters
    ----------
    two_theta_deg
        Angles in degrees, used only for shape.
    zero_deg
        The zero-point error in degrees ``2 theta``.

    Returns
    -------
    np.ndarray
        The additive error at each angle, in degrees.
    """

    angles = np.asarray(two_theta_deg, dtype=np.float64)
    return np.full_like(angles, float(zero_deg))


def specimen_displacement_shift_deg(
    two_theta_deg: Any,
    *,
    displacement_mm: float,
    goniometer_radius_mm: float,
) -> np.ndarray:
    """Return the peak shift caused by a specimen off the diffractometer axis.

    Purpose
    -------
    Quantify the aberration that dominates precise lattice-parameter work on a
    laboratory Bragg-Brentano instrument, and that no amount of averaging over
    reflections can remove.

    Method
    ------
    For a flat specimen whose surface sits a distance ``s`` from the focusing
    circle -- positive when displaced away from the source -- the diffracted
    beam converges at the wrong place on the receiving circle, and to first
    order in ``s / R`` (Wilson 1963)

    ``Delta(2 theta) = -2 s cos(theta) / R``   [radians]

    The ``cos(theta)`` dependence is the important part: the shift vanishes at
    ``2 theta = 180`` degrees, which is precisely why extrapolating a lattice
    parameter to ``theta = 90`` degrees works, and why the Nelson-Riley
    function contains a ``cos^2(theta) / sin(theta)`` term.

    Parameters
    ----------
    two_theta_deg
        Angles in degrees.
    displacement_mm
        Specimen surface displacement ``s`` in millimetres. Typical
        preparation errors are tens of micrometres.
    goniometer_radius_mm
        Goniometer radius ``R`` in millimetres. Laboratory instruments are
        usually 150 to 300 mm.

    Returns
    -------
    np.ndarray
        The additive shift at each angle, in degrees ``2 theta``. Negative for
        a positive displacement.

    Raises
    ------
    ValueError
        If the goniometer radius is not positive.
    """

    if not np.isfinite(goniometer_radius_mm) or goniometer_radius_mm <= 0.0:
        raise ValueError("goniometer_radius_mm must be finite and positive.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    theta = np.deg2rad(0.5 * angles)
    shift: np.ndarray = np.rad2deg(
        -2.0 * float(displacement_mm) * np.cos(theta) / float(goniometer_radius_mm)
    )
    return shift


def specimen_transparency_shift_deg(
    two_theta_deg: Any,
    *,
    linear_absorption_coefficient_inv_mm: float,
    goniometer_radius_mm: float,
) -> np.ndarray:
    """Return the peak shift caused by beam penetration into the specimen.

    Purpose
    -------
    Account for the fact that diffraction happens at a mean depth below the
    surface rather than at the surface, which for a weakly absorbing specimen
    displaces peaks by more than the specimen was displaced.

    Method
    ------
    For a specimen thick enough to be effectively infinite, the centroid of
    the diffracting volume lies at a depth ``sin(theta) / (2 mu)``, and the
    resulting shift is (Wilson 1963)

    ``Delta(2 theta) = -sin(2 theta) / (2 mu R)``   [radians]

    The ``sin(2 theta)`` dependence peaks at 45 degrees ``theta`` and vanishes
    at both ends, which is a different signature from displacement and zero --
    the reason all three can be modelled rather than lumped together.

    Parameters
    ----------
    two_theta_deg
        Angles in degrees.
    linear_absorption_coefficient_inv_mm
        ``mu`` in inverse millimetres. Steels and other dense metals for Cu
        K-alpha are of order 100 to 500 per mm, where this term is negligible;
        light-element powders, polymers and organics are of order 1 to 10 per
        mm, where it is not.
    goniometer_radius_mm
        Goniometer radius ``R`` in millimetres.

    Returns
    -------
    np.ndarray
        The additive shift at each angle, in degrees ``2 theta``.

    Raises
    ------
    ValueError
        If ``mu`` or ``R`` is not positive.
    """

    if (
        not np.isfinite(linear_absorption_coefficient_inv_mm)
        or linear_absorption_coefficient_inv_mm <= 0.0
    ):
        raise ValueError("linear_absorption_coefficient_inv_mm must be finite and positive.")
    if not np.isfinite(goniometer_radius_mm) or goniometer_radius_mm <= 0.0:
        raise ValueError("goniometer_radius_mm must be finite and positive.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    shift: np.ndarray = np.rad2deg(
        -np.sin(np.deg2rad(angles))
        / (
            2.0
            * float(linear_absorption_coefficient_inv_mm)
            * float(goniometer_radius_mm)
        )
    )
    return shift


def refraction_decrement(
    *,
    density_g_cm3: float,
    wavelength_angstrom: float,
    electrons_per_gram: float | None = None,
) -> float:
    """Return the refractive-index decrement ``delta = 1 - n`` of a specimen.

    Purpose
    -------
    Supply the one material constant needed for the refraction correction, at
    the level of accuracy that correction deserves.

    Method
    ------
    Well away from an absorption edge the classical Drude result holds:
    ``delta = r_e lambda^2 n_e / (2 pi)``, with ``r_e`` the classical electron
    radius and ``n_e`` the electron number density. ``n_e`` is taken as
    ``rho N_A (Z / A)``, and ``Z / A`` is close to 0.5 for every element except
    hydrogen, which is why a density and a wavelength are enough.

    Parameters
    ----------
    density_g_cm3
        Specimen density in grams per cubic centimetre.
    wavelength_angstrom
        Radiation wavelength in angstrom.
    electrons_per_gram
        Electrons per gram divided by Avogadro's number, i.e. the mean
        ``Z / A``. Defaults to 0.5.

    Returns
    -------
    float
        ``delta``, a dimensionless number of order 1e-5 for a metal at Cu
        K-alpha.

    Raises
    ------
    ValueError
        If the density or wavelength is not positive.
    """

    if not np.isfinite(density_g_cm3) or density_g_cm3 <= 0.0:
        raise ValueError("density_g_cm3 must be finite and positive.")
    if not np.isfinite(wavelength_angstrom) or wavelength_angstrom <= 0.0:
        raise ValueError("wavelength_angstrom must be finite and positive.")
    z_over_a = 0.5 if electrons_per_gram is None else float(electrons_per_gram)
    # Electron density in electrons per cubic angstrom: g/cm^3 -> g/A^3 is 1e-24.
    electron_density = float(density_g_cm3) * _AVOGADRO * z_over_a * 1.0e-24
    return float(
        _ELECTRON_RADIUS_ANGSTROM
        * float(wavelength_angstrom) ** 2
        * electron_density
        / (2.0 * np.pi)
    )


def refraction_shift_deg(two_theta_deg: Any, *, decrement: float) -> np.ndarray:
    """Return the apparent peak shift caused by refraction at the surface.

    Purpose
    -------
    Include the last systematic worth naming at the 1e-5 level in ``a``, which
    is exactly the level precise lattice-parameter work operates at.

    Method
    ------
    Inside the specimen the wavelength is longer by the factor ``1 / n``, so
    Bragg's law becomes ``lambda = 2 d sin(theta) (1 - delta / sin^2(theta))``.
    Solving for the shift in the *apparent* angle to first order in ``delta``,

    ``Delta(2 theta) = 2 delta / (sin(theta) cos(theta)) = 4 delta / sin(2 theta)``
    [radians]

    Unlike the geometric aberrations this one grows towards low angle, so it is
    small where the others are large and vice versa.

    Parameters
    ----------
    two_theta_deg
        Angles in degrees.
    decrement
        ``delta = 1 - n``, from :func:`refraction_decrement`.

    Returns
    -------
    np.ndarray
        The additive shift at each angle, in degrees ``2 theta``, positive
        because the observed angle is larger than the ideal one.

    Raises
    ------
    ValueError
        If the decrement is negative, or any angle is at 0 or 180 degrees where
        the expression diverges.
    """

    if not np.isfinite(decrement) or decrement < 0.0:
        raise ValueError("refraction_shift_deg requires a finite, non-negative decrement.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    sine = np.sin(np.deg2rad(angles))
    if np.any(np.abs(sine) < 1.0e-9):
        raise ValueError(
            "The refraction shift diverges at 0 and 180 degrees 2*theta; restrict the range."
        )
    shift: np.ndarray = np.rad2deg(4.0 * float(decrement) / sine)
    return shift


def position_correction_deg(
    two_theta_deg: Any,
    *,
    zero_deg: float = 0.0,
    displacement_mm: float = 0.0,
    goniometer_radius_mm: float = 240.0,
    linear_absorption_coefficient_inv_mm: float | None = None,
    refraction_decrement_value: float = 0.0,
) -> np.ndarray:
    """Return the total additive position aberration at each angle.

    Purpose
    -------
    Compose the individual aberrations into the one quantity a corrected
    position needs, while keeping each contribution separately inspectable.

    Method
    ------
    The four terms are additive to first order and are simply summed. Their
    different angular dependences -- constant, ``cos(theta)``,
    ``sin(2 theta)``, ``1 / sin(2 theta)`` -- are what makes a simultaneous fit
    of more than one of them possible in principle, and what makes it
    ill-conditioned in practice over a single scan.

    Parameters
    ----------
    two_theta_deg
        Observed angles in degrees.
    zero_deg
        Detector zero error, from a standard.
    displacement_mm
        Specimen surface displacement.
    goniometer_radius_mm
        Goniometer radius.
    linear_absorption_coefficient_inv_mm
        Specimen ``mu``. ``None`` omits the transparency term, which is the
        right choice for a strongly absorbing metal.
    refraction_decrement_value
        ``delta``, from :func:`refraction_decrement`. Zero omits refraction.

    Returns
    -------
    np.ndarray
        The total shift, to be *subtracted* from an observed angle to obtain
        the ideal one.
    """

    angles = np.asarray(two_theta_deg, dtype=np.float64)
    total = zero_shift_deg(angles, zero_deg=zero_deg)
    if displacement_mm != 0.0:
        total = total + specimen_displacement_shift_deg(
            angles,
            displacement_mm=displacement_mm,
            goniometer_radius_mm=goniometer_radius_mm,
        )
    if linear_absorption_coefficient_inv_mm is not None:
        total = total + specimen_transparency_shift_deg(
            angles,
            linear_absorption_coefficient_inv_mm=linear_absorption_coefficient_inv_mm,
            goniometer_radius_mm=goniometer_radius_mm,
        )
    if refraction_decrement_value != 0.0:
        total = total + refraction_shift_deg(angles, decrement=refraction_decrement_value)
    return total


def correct_peak_positions(
    table: PeakTable,
    *,
    zero_deg: float = 0.0,
    displacement_mm: float = 0.0,
    goniometer_radius_mm: float = 240.0,
    linear_absorption_coefficient_inv_mm: float | None = None,
    refraction_decrement_value: float = 0.0,
) -> PeakTable:
    """Return a peak table with the instrumental aberrations removed.

    Purpose
    -------
    Apply the position corrections where they belong -- to a handful of fitted
    positions -- rather than by resampling a whole profile, so that the
    uncertainties attached to those positions travel with them unchanged.

    Method
    ------
    Each corrected position is ``2 theta_obs - Delta(2 theta)`` with
    ``Delta`` from :func:`position_correction_deg`. The standard uncertainties
    are carried through untouched: a *known* aberration is a bias, not a
    random error, so removing it does not change how well the position is
    known. If the aberration parameters are themselves uncertain, that
    uncertainty belongs to them and is propagated by the lattice-parameter
    determination that refines them, not here.

    Parameters
    ----------
    table
        The fitted peaks.
    zero_deg, displacement_mm, goniometer_radius_mm
    linear_absorption_coefficient_inv_mm, refraction_decrement_value
        As for :func:`position_correction_deg`.

    Returns
    -------
    PeakTable
        A new table with corrected positions and the correction recorded in
        its settings.

    See Also
    --------
    pytex.diffraction.xrd_lattice_parameter.determine_lattice_parameters :
        refines a systematic-error coefficient instead of requiring these
        values to be known in advance.
    """

    positions = table.two_theta_deg
    correction = position_correction_deg(
        positions,
        zero_deg=zero_deg,
        displacement_mm=displacement_mm,
        goniometer_radius_mm=goniometer_radius_mm,
        linear_absorption_coefficient_inv_mm=linear_absorption_coefficient_inv_mm,
        refraction_decrement_value=refraction_decrement_value,
    )
    corrected = tuple(
        replace(peak, two_theta_deg=float(position - shift))
        for peak, position, shift in zip(table.peaks, positions, correction, strict=True)
    )
    settings = dict(table.settings)
    settings["position_correction_zero_deg"] = float(zero_deg)
    settings["position_correction_displacement_mm"] = float(displacement_mm)
    settings["position_correction_goniometer_radius_mm"] = float(goniometer_radius_mm)
    return PeakTable(
        name=f"{table.name} (aberration corrected)",
        peaks=corrected,
        radiation=table.radiation,
        source_name=table.source_name,
        settings=settings,
    )


# ---------------------------------------------------------------------------
# Intensity corrections
# ---------------------------------------------------------------------------


def monochromator_polarization_factor(
    two_theta_deg: Any,
    *,
    monochromator_two_theta_deg: float,
    perpendicular: bool = False,
) -> np.ndarray:
    """Return the polarization factor for a beam conditioned by a monochromator.

    Purpose
    -------
    Replace the unpolarized ``(1 + cos^2(2 theta)) / 2`` factor when the beam
    has already been polarized by Bragg reflection from a crystal
    monochromator, which is the normal laboratory configuration and is
    frequently left uncorrected.

    Method
    ------
    A monochromator set at ``2 theta_M`` partially polarizes the beam. For the
    usual *parallel* (diffracted-beam, same plane) geometry

    ``P = (1 + cos^2(2 theta_M) cos^2(2 theta)) / (1 + cos^2(2 theta_M))``

    and for the perpendicular geometry the two cosines exchange roles:

    ``P = (cos^2(2 theta_M) + cos^2(2 theta)) / (1 + cos^2(2 theta_M))``.

    The two forms come from tracking the two polarization components
    separately. A monochromator reflects the component perpendicular to *its*
    diffraction plane with efficiency one and the parallel component with
    ``cos^2(2 theta_M)``; the specimen then does the same with
    ``cos^2(2 theta)``. In the parallel setting the two perpendicular
    directions coincide, so the products pair as ``1 x 1`` and
    ``cos^2(2 theta_M) x cos^2(2 theta)``; in the perpendicular setting the
    monochromator's perpendicular component is the specimen's parallel one, and
    the products cross over.

    The limits follow. ``2 theta_M = 0`` leaves the beam unpolarized and both
    forms collapse to ``(1 + cos^2(2 theta)) / 2``. ``2 theta_M = 90`` degrees
    is an ideal polarizer: in the *parallel* setting the surviving component is
    the one the specimen scatters without angular dependence, so ``P = 1`` and
    the polarization correction disappears entirely; in the *perpendicular*
    setting the surviving component is the one the specimen modulates fully, so
    ``P = cos^2(2 theta)``. Because ``(1 + m c) - (m + c) = (1 - m)(1 - c)``
    with ``m, c`` in ``[0, 1]``, the parallel factor is never smaller than the
    perpendicular one.

    Parameters
    ----------
    two_theta_deg
        Specimen scattering angles in degrees.
    monochromator_two_theta_deg
        The monochromator's own ``2 theta``. Graphite (0002) with Cu K-alpha
        is about 26.6 degrees, where the effect is small but not zero.
    perpendicular
        Use the perpendicular geometry.

    Returns
    -------
    np.ndarray
        The polarization factor, normalized to one at ``2 theta = 0``.

    Raises
    ------
    ValueError
        If the monochromator angle is outside ``[0, 180]`` degrees.
    """

    if not 0.0 <= float(monochromator_two_theta_deg) <= 180.0:
        raise ValueError("monochromator_two_theta_deg must lie in [0, 180] degrees.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    specimen = np.square(np.cos(np.deg2rad(angles)))
    monochromator = float(np.cos(np.deg2rad(monochromator_two_theta_deg))) ** 2
    numerator = (
        monochromator + specimen if perpendicular else 1.0 + monochromator * specimen
    )
    factor: np.ndarray = numerator / (1.0 + monochromator)
    return factor


def variable_to_fixed_slit(measured: MeasuredPowderPattern) -> MeasuredPowderPattern:
    """Convert an automatic-divergence-slit scan to fixed-slit intensities.

    Purpose
    -------
    Make a scan collected with a variable (constant irradiated length) slit
    comparable with the fixed-slit intensities that structure factors,
    reference databases and Rietveld programs assume.

    Method
    ------
    A variable slit opens as ``theta`` increases so that the irradiated length
    on the specimen stays constant; a fixed slit illuminates a constant
    *angular* aperture, so the irradiated length grows as ``1 / sin(theta)``
    and a progressively larger fraction of the beam spills off the specimen.
    Comparing the two illuminated areas gives (Klug & Alexander 1974)

    ``I_fixed = I_variable * sin(theta)``.

    Positions are untouched; only relative intensities change, strongly
    suppressing the low-angle end.

    Parameters
    ----------
    measured
        A scan collected with a variable divergence slit.

    Returns
    -------
    MeasuredPowderPattern
        The converted pattern, in ``"arbitrary"`` intensity units because the
        result is no longer a count, with the conversion recorded in metadata.

    Notes
    -----
    The conversion is exact only for an ideal, infinitely thick, flat specimen
    fully covered at every angle. Near the low-angle end of a real scan the
    beam may already have overflowed the specimen, and no analytic conversion
    repairs that.
    """

    axis = np.asarray(measured.two_theta_deg, dtype=np.float64)
    converted = np.asarray(measured.intensity, dtype=np.float64) * np.sin(
        np.deg2rad(0.5 * axis)
    )
    uncertainty = (
        None
        if measured.standard_uncertainty is None
        else np.asarray(measured.standard_uncertainty, dtype=np.float64)
        * np.sin(np.deg2rad(0.5 * axis))
    )
    return replace(
        measured,
        name=f"{measured.name} (fixed-slit equivalent)",
        intensity=converted,
        standard_uncertainty=uncertainty,
        intensity_unit="arbitrary",
        metadata=_record(measured, "variable_to_fixed_slit"),
    )


def smooth_savitzky_golay(
    measured: MeasuredPowderPattern,
    *,
    window_points: int = 11,
    polynomial_order: int = 3,
) -> MeasuredPowderPattern:
    """Smooth a noisy pattern, for display only.

    Purpose
    -------
    Make a weak, noisy scan legible.

    .. warning::

       **Never smooth before a width analysis, and prefer not to smooth before
       a position analysis.** A Savitzky-Golay filter is a local polynomial
       least-squares fit, and convolving the data with it convolves the peak
       shape too: measured widths increase, so Scherrer sizes come out too
       small and Williamson-Hall strains too large. It also correlates
       neighbouring points, which invalidates the independent-observation
       assumption behind every uncertainty this library reports. Profile
       *fitting* already performs the noise averaging that smoothing attempts,
       without either side effect.

    Method
    ------
    Savitzky & Golay (1964): fit a polynomial of ``polynomial_order`` to each
    sliding window of ``window_points`` measured values by least squares and
    take the fitted value at the window centre. The uncertainty array, if
    present, is discarded rather than propagated, because after smoothing the
    points are no longer independent and a per-point sigma would be
    misleading.

    Parameters
    ----------
    measured
        The pattern to smooth.
    window_points
        Window length in data points. Must be odd and greater than
        ``polynomial_order``.
    polynomial_order
        Order of the local polynomial. Three preserves peak height and
        curvature far better than one.

    Returns
    -------
    MeasuredPowderPattern
        The smoothed pattern, in ``"arbitrary"`` units, with uncertainties
        dropped and the operation recorded in metadata.

    Raises
    ------
    ValueError
        If the window is even, too short, or longer than the pattern.
    """

    if window_points % 2 == 0:
        raise ValueError("smooth_savitzky_golay requires an odd window_points.")
    if window_points <= polynomial_order:
        raise ValueError("smooth_savitzky_golay requires window_points > polynomial_order.")
    if window_points > len(measured):
        raise ValueError("smooth_savitzky_golay window is longer than the measured pattern.")
    smoothed = savgol_filter(
        np.asarray(measured.intensity, dtype=np.float64),
        window_length=int(window_points),
        polyorder=int(polynomial_order),
    )
    return replace(
        measured,
        name=f"{measured.name} (smoothed)",
        intensity=np.maximum(smoothed, 0.0),
        standard_uncertainty=None,
        intensity_unit="arbitrary",
        metadata=_record(
            measured,
            f"savitzky_golay(window={int(window_points)}, order={int(polynomial_order)})",
        ),
    )


# ---------------------------------------------------------------------------
# Display transforms
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ProfileView:
    """A powder profile expressed on a chosen abscissa and intensity scale.

    Purpose
    -------
    Carry a plottable version of a pattern together with the axis meaning and
    the transform that produced it, so a figure can label itself correctly and
    a reader can tell a square-root ordinate from a linear one without being
    told.

    This is deliberately *not* a
    :class:`~pytex.diffraction.xrd_measurement.MeasuredPowderPattern`: once the
    abscissa is ``d`` or ``Q`` the array is no longer degrees ``2 theta``, and
    once the ordinate is logarithmic it is no longer an intensity. Keeping the
    types apart stops a transformed view being fed back into an analysis that
    assumes measured units.

    Attributes
    ----------
    name : str
        A human name for the view.
    abscissa : np.ndarray
        The transformed horizontal axis, ascending or descending depending on
        the choice (``d`` decreases as ``2 theta`` increases).
    abscissa_kind : str
        One of :data:`PROFILE_ABSCISSAE`.
    abscissa_label : str
        A display label, with units.
    ordinate : np.ndarray
        The transformed intensity.
    ordinate_label : str
        A display label naming the scale and any normalization.
    scale : str
        One of :data:`INTENSITY_SCALES`.
    normalization : str
        One of :data:`INTENSITY_NORMALIZATIONS`.
    source_name : str
        Name of the pattern this view was made from.
    metadata : Mapping[str, str]
        Provenance carried from the source pattern.
    """

    name: str
    abscissa: np.ndarray
    abscissa_kind: ProfileAbscissa
    abscissa_label: str
    ordinate: np.ndarray
    ordinate_label: str
    scale: IntensityScale
    normalization: IntensityNormalization
    source_name: str
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        abscissa = as_float_array(self.abscissa, shape=(None,))
        ordinate = as_float_array(self.ordinate, shape=(None,))
        if abscissa.shape != ordinate.shape:
            raise ValueError("ProfileView abscissa and ordinate must have the same shape.")
        if abscissa.size < 2:
            raise ValueError("A ProfileView needs at least two points.")
        if np.any(~np.isfinite(abscissa)) or np.any(~np.isfinite(ordinate)):
            raise ValueError("ProfileView arrays must be finite.")
        if self.abscissa_kind not in PROFILE_ABSCISSAE:
            raise ValueError(f"ProfileView.abscissa_kind must be one of {PROFILE_ABSCISSAE}.")
        if self.scale not in INTENSITY_SCALES:
            raise ValueError(f"ProfileView.scale must be one of {INTENSITY_SCALES}.")
        if self.normalization not in INTENSITY_NORMALIZATIONS:
            raise ValueError(
                f"ProfileView.normalization must be one of {INTENSITY_NORMALIZATIONS}."
            )
        object.__setattr__(self, "abscissa", abscissa)
        object.__setattr__(self, "ordinate", ordinate)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def __len__(self) -> int:
        return int(self.abscissa.size)

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this view."""

        return {
            "schema": PROFILE_VIEW_SCHEMA,
            "name": self.name,
            "source_name": self.source_name,
            "abscissa_kind": self.abscissa_kind,
            "abscissa_label": self.abscissa_label,
            "ordinate_label": self.ordinate_label,
            "scale": self.scale,
            "normalization": self.normalization,
            "abscissa": [float(value) for value in self.abscissa],
            "ordinate": [float(value) for value in self.ordinate],
        }

    def describe(self) -> str:
        """Return convention-explicit prose about what this view shows."""

        reading = {
            "two_theta_deg": (
                "The measured abscissa, and the only one on which instrumental aberrations have "
                "their tabulated angular form"
            ),
            "d_angstrom": (
                "Interplanar spacing, which runs backwards relative to angle and makes the same "
                "phase look identical whatever wavelength measured it"
            ),
            "q_inv_angstrom": (
                "Scattering vector magnitude Q = 4 pi sin(theta) / lambda, the wavelength-free "
                "abscissa that makes X-ray and neutron patterns directly comparable"
            ),
            "sin_squared_theta": (
                "The abscissa in which Bragg's law is linear in the reciprocal metric tensor, "
                "which is why lattice-parameter least squares is done here"
            ),
        }[self.abscissa_kind]
        scale = {
            "linear": (
                "A linear ordinate, on which a weak reflection next to a strong one is invisible"
            ),
            "sqrt": (
                "A square-root ordinate, which is the variance-stabilizing scale for counting "
                "statistics: noise has the same visual amplitude everywhere, so a feature that "
                "looks significant is significant"
            ),
            "log10": (
                "A base-ten logarithmic ordinate, which shows several decades at once but "
                "exaggerates background structure and cannot display a zero count"
            ),
        }[self.scale]
        normalization = {
            "none": "Intensities are unnormalized.",
            "maximum": "Intensities are scaled so the strongest point is one.",
            "integral": "Intensities are scaled to unit integrated area on this abscissa.",
        }[self.normalization]
        return (
            f"Profile view '{self.name}' of '{self.source_name}' plots {self.ordinate_label} "
            f"against {self.abscissa_label} over {len(self)} points from "
            f"{float(self.abscissa[0]):.5g} to {float(self.abscissa[-1]):.5g}. {reading}. "
            f"{scale}. {normalization}"
        )


def profile_view(
    measured: MeasuredPowderPattern,
    *,
    abscissa: ProfileAbscissa = "two_theta_deg",
    scale: IntensityScale = "linear",
    normalization: IntensityNormalization = "none",
    radiation: RadiationSpec | None = None,
    log_floor_fraction: float = 1.0e-4,
) -> ProfileView:
    """Return a plottable view of a pattern on a chosen abscissa and scale.

    Purpose
    -------
    Give every display surface -- the workbench panel, a notebook figure, an
    export -- one transform with one set of conventions, instead of each
    reimplementing ``4 pi sin(theta) / lambda`` slightly differently.

    Method
    ------
    Abscissa transforms are the standard identities, all at fixed wavelength:

    * ``d = lambda / (2 sin(theta))``,
    * ``Q = 4 pi sin(theta) / lambda = 2 pi / d``,
    * ``sin^2(theta)``, in which Bragg's law becomes
      ``sin^2(theta) = (lambda^2 / 4) h^T G* h``, linear in the reciprocal
      metric tensor. This is the abscissa lattice-parameter least squares works
      in, which is why it is offered for looking at as well.

    Ordinate transforms:

    * ``"sqrt"`` is the variance-stabilizing scale for Poisson counts. Its
      virtue is not aesthetic: after the transform the noise amplitude is the
      same everywhere, so visual significance and statistical significance
      agree. It should be the default for judging whether a weak feature is
      real.
    * ``"log10"`` shows decades but cannot represent a zero count, so values
      below ``log_floor_fraction`` of the maximum are clipped to that floor;
      the clip is stated in the ordinate label.

    Parameters
    ----------
    measured
        The pattern to transform.
    abscissa
        One of :data:`PROFILE_ABSCISSAE`.
    scale
        One of :data:`INTENSITY_SCALES`.
    normalization
        One of :data:`INTENSITY_NORMALIZATIONS`. Normalization is applied
        before the scale transform, so a logarithm of a normalized pattern is
        negative rather than clipped away.
    radiation
        Radiation supplying the wavelength for ``d`` and ``Q``. Falls back to
        ``measured.radiation``.
    log_floor_fraction
        Clip floor for the logarithmic scale, as a fraction of the maximum.

    Returns
    -------
    ProfileView
        The transformed view, labelled and self-describing.

    Raises
    ------
    ValueError
        If a choice is unknown, a wavelength is needed but unavailable, or the
        clip floor is outside ``(0, 1)``.
    """

    if abscissa not in PROFILE_ABSCISSAE:
        raise ValueError(f"profile_view requires abscissa in {PROFILE_ABSCISSAE}.")
    if scale not in INTENSITY_SCALES:
        raise ValueError(f"profile_view requires scale in {INTENSITY_SCALES}.")
    if normalization not in INTENSITY_NORMALIZATIONS:
        raise ValueError(f"profile_view requires normalization in {INTENSITY_NORMALIZATIONS}.")
    if not 0.0 < log_floor_fraction < 1.0:
        raise ValueError("log_floor_fraction must lie strictly inside (0, 1).")

    angles = np.asarray(measured.two_theta_deg, dtype=np.float64)
    intensity = np.asarray(measured.intensity, dtype=np.float64)
    sine = np.sin(np.deg2rad(0.5 * angles))

    if abscissa == "two_theta_deg":
        values = angles
        label = "2*theta (degrees)"
    elif abscissa == "sin_squared_theta":
        values = np.square(sine)
        label = "sin^2(theta)"
    else:
        spec = radiation if radiation is not None else measured.radiation
        if spec is None:
            raise ValueError(
                f"The '{abscissa}' abscissa needs a wavelength, and neither the pattern nor the "
                "call supplied a radiation."
            )
        wavelength = float(spec.wavelength_angstrom)
        if np.any(sine <= 0.0):
            raise ValueError(
                "A spacing or scattering-vector abscissa is undefined at 2*theta = 0; restrict "
                "the measured range."
            )
        if abscissa == "d_angstrom":
            values = wavelength / (2.0 * sine)
            label = "d (angstrom)"
        else:
            values = 4.0 * np.pi * sine / wavelength
            label = "Q (1/angstrom)"

    ordinate = intensity
    if normalization == "maximum":
        peak = float(np.max(ordinate))
        if peak <= 0.0:
            raise ValueError("Cannot normalize to the maximum of an all-zero pattern.")
        ordinate = ordinate / peak
        normalization_label = "normalized"
    elif normalization == "integral":
        order = np.argsort(values)
        area = float(np.trapezoid(ordinate[order], values[order]))
        if area == 0.0:
            raise ValueError("Cannot normalize to a zero integrated area.")
        ordinate = ordinate / abs(area)
        normalization_label = "area-normalized"
    else:
        normalization_label = ""

    unit = measured.intensity_unit.replace("_", " ")
    if scale == "linear":
        ordinate_label = f"intensity ({normalization_label or unit})"
    elif scale == "sqrt":
        ordinate = np.sqrt(np.maximum(ordinate, 0.0))
        ordinate_label = f"sqrt(intensity) ({normalization_label or unit})"
    else:
        peak = float(np.max(ordinate))
        if peak <= 0.0:
            raise ValueError("Cannot take the logarithm of an all-zero pattern.")
        floor = peak * float(log_floor_fraction)
        ordinate = np.log10(np.maximum(ordinate, floor))
        ordinate_label = (
            f"log10(intensity) ({normalization_label or unit}), "
            f"clipped at {log_floor_fraction:g} of maximum"
        )

    return ProfileView(
        name=f"{measured.name} [{abscissa}, {scale}]",
        abscissa=values,
        abscissa_kind=abscissa,
        abscissa_label=label,
        ordinate=ordinate,
        ordinate_label=ordinate_label,
        scale=scale,
        normalization=normalization,
        source_name=measured.name,
        metadata=dict(measured.metadata),
    )


__all__ = [
    "INTENSITY_NORMALIZATIONS",
    "INTENSITY_SCALES",
    "PROFILE_ABSCISSAE",
    "PROFILE_VIEW_SCHEMA",
    "IntensityNormalization",
    "IntensityScale",
    "ProfileAbscissa",
    "ProfileView",
    "correct_peak_positions",
    "kalpha2_partner_two_theta_deg",
    "monochromator_polarization_factor",
    "position_correction_deg",
    "profile_view",
    "refraction_decrement",
    "refraction_shift_deg",
    "smooth_savitzky_golay",
    "specimen_displacement_shift_deg",
    "specimen_transparency_shift_deg",
    "strip_kalpha2",
    "variable_to_fixed_slit",
    "zero_shift_deg",
]
