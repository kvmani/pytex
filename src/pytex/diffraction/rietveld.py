"""Whole-profile (Rietveld) refinement of a powder pattern against a structure.

Rietveld's insight was that a powder pattern should not be reduced to a list of
integrated intensities before it is used. Overlapping reflections cannot be
separated reliably, but they can be *calculated* from a structural model and
compared point by point with the measurement, so the fit is done on the profile
itself and the overlap never has to be resolved.

What is refined here
--------------------
A curated set of parameters, each of which is a statement about a distinct piece
of physics, so a user can tell which part of the model a refinement is actually
adjusting:

============================ ==========================================================
``scale``                    Overall intensity scale. Always refined; without it every
                             other parameter fits the wrong amplitude.
``zero_shift_deg``           Detector zero error. An *instrument* misalignment, so
                             refining it is what stops the lattice parameter from
                             silently absorbing it.
``lattice_scale``            Isotropic dilation of the unit cell. Symmetry-preserving
                             in every crystal system, which is why it is the safe
                             general cell parameter.
``caglioti_u/v/w``           The Gaussian angular width, as in
                             :class:`~pytex.diffraction.xrd_instrument.InstrumentBroadening`.
``lorentzian_y``             Size-like Lorentzian width, ``Y / cos(theta)``.
``b_iso_overall``            One isotropic displacement parameter applied to every
                             site. It correlates strongly with scale and with the
                             background, which is why it is not refined by default.
``march_coefficient``        March-Dollase preferred orientation about a stated axis.
``background_0 .. N``        Chebyshev background coefficients, refined jointly with
                             everything else rather than subtracted beforehand.
============================ ==========================================================

What is deliberately *not* refined
----------------------------------
Individual atomic coordinates, site occupancies, anisotropic displacement
parameters, and independent cell edges of a low-symmetry cell. Those are the
parameters that make a refinement a structure determination, and they need
constraints, restraints, and a rigid-body or symmetry-mode apparatus that this
module does not have. Refining them without that apparatus produces a lower
``R_wp`` and a structure nobody should publish, which is a worse outcome than
not offering them.

This is therefore a *profile and phase-quantification* refinement, honest about
its scope, and the right tool for lattice parameters, peak shape, texture
strength, background, and goodness of fit against a known structure.

References
----------
Rietveld, H. M., *J. Appl. Crystallogr.* **2** (1969) 65-71,
doi:10.1107/S0021889869006558 -- the method.

Young, R. A. (ed.), *The Rietveld Method*, IUCr/OUP (1993) -- ``R_p``, ``R_wp``,
``R_exp``, and the goodness of fit as used here.

Toby, B. H., *Powder Diffr.* **21** (2006) 67-70, doi:10.1154/1.2179804 --
"R factors in Rietveld analysis: how good is good enough?", the source for the
interpretation given in :meth:`RietveldResult.describe`.

Hill, R. J. & Flack, H. D., *J. Appl. Crystallogr.* **20** (1987) 356-361,
doi:10.1107/S0021889887086485 -- the Durbin-Watson statistic for serial
correlation in powder-profile residuals.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Literal

import numpy as np
from scipy.optimize import least_squares

from pytex.core._arrays import as_float_array
from pytex.core.lattice import Phase
from pytex.core.miller import MillerPlane
from pytex.diffraction.preferred_orientation import MarchDollaseModel
from pytex.diffraction.xrd import (
    PowderPattern,
    PowderReflection,
    RadiationSpec,
    generate_powder_reflections,
)
from pytex.diffraction.xrd_instrument import InstrumentBroadening, _tch_eta, _tch_width
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

RIETVELD_REFINEMENT_SCHEMA = "pytex.diffraction.rietveld_refinement"

WeightModel = Literal["unit", "poisson", "inverse_variance"]

#: Every parameter name :func:`refine_rietveld` understands, other than the
#: background coefficients, which are named ``background_0`` upward.
REFINABLE_PARAMETERS = (
    "scale",
    "zero_shift_deg",
    "lattice_scale",
    "caglioti_u",
    "caglioti_v",
    "caglioti_w",
    "lorentzian_y",
    "b_iso_overall",
    "march_coefficient",
)

#: The default refinement set: scale, the cell, the instrument zero, one width,
#: and the background. These are the parameters that are almost always worth
#: refining and are least likely to run away when refined together.
DEFAULT_REFINEMENT_SET = (
    "scale",
    "zero_shift_deg",
    "lattice_scale",
    "caglioti_w",
)

_PARAMETER_BOUNDS: dict[str, tuple[float, float]] = {
    "scale": (1e-12, np.inf),
    # A zero error beyond a degree is a broken diffractometer, not a refinable
    # offset, and letting it wander further lets it swap places with the cell.
    "zero_shift_deg": (-1.0, 1.0),
    # Ten percent is far outside any real thermal or compositional cell change;
    # a refinement that wants more has misidentified the phase.
    "lattice_scale": (0.9, 1.1),
    "caglioti_u": (-1.0, 1.0),
    "caglioti_v": (-1.0, 1.0),
    "caglioti_w": (1e-8, 1.0),
    "lorentzian_y": (0.0, 1.0),
    "b_iso_overall": (0.0, 20.0),
    "march_coefficient": (0.2, 5.0),
}

_PARAMETER_UNITS: dict[str, str] = {
    "scale": "",
    "zero_shift_deg": "deg 2theta",
    "lattice_scale": "",
    "caglioti_u": "deg^2",
    "caglioti_v": "deg^2",
    "caglioti_w": "deg^2",
    "lorentzian_y": "deg",
    "b_iso_overall": "A^2",
    "march_coefficient": "",
}

_PARAMETER_DESCRIPTIONS: dict[str, str] = {
    "scale": "Overall intensity scale relating calculated to measured intensity.",
    "zero_shift_deg": "Detector zero-point error added to every calculated 2*theta.",
    "lattice_scale": "Isotropic dilation of the unit cell about its published parameters.",
    "caglioti_u": "Caglioti U: the tan^2(theta) term of the squared Gaussian width.",
    "caglioti_v": "Caglioti V: the tan(theta) term of the squared Gaussian width.",
    "caglioti_w": "Caglioti W: the angle-independent term of the squared Gaussian width.",
    "lorentzian_y": "Lorentzian Y: size-like broadening varying as 1 / cos(theta).",
    "b_iso_overall": "One isotropic atomic displacement parameter shared by every site.",
    "march_coefficient": "March-Dollase texture strength about the stated preferred axis.",
}

_CITATION_RIETVELD = "Rietveld, J. Appl. Crystallogr. 2 (1969) 65, doi:10.1107/S0021889869006558."
_CITATION_TOBY = "Toby, Powder Diffr. 21 (2006) 67, doi:10.1154/1.2179804."
_CITATION_HILL_FLACK = (
    "Hill & Flack, J. Appl. Crystallogr. 20 (1987) 356, doi:10.1107/S0021889887086485."
)


@dataclass(frozen=True, slots=True)
class RefinedParameter:
    """One refinement parameter: what it means, where it started, where it went.

    Attributes
    ----------
    name : str
    value : float
        The value after refinement, or the fixed value when ``refined`` is false.
    initial_value : float
    standard_uncertainty : float or None
        One estimated standard deviation from the least-squares covariance,
        scaled by the goodness of fit. ``None`` for a fixed parameter. These are
        *precision* estimates conditional on the model being right; they say
        nothing about whether the model is right, and are routinely optimistic
        by a factor of two or three because neighbouring profile points are
        correlated.
    refined : bool
    units : str
    description : str
    """

    name: str
    value: float
    initial_value: float
    standard_uncertainty: float | None
    refined: bool
    units: str
    description: str

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("RefinedParameter.name must be non-empty.")
        for field_name in ("value", "initial_value"):
            if not np.isfinite(getattr(self, field_name)):
                raise ValueError(f"RefinedParameter.{field_name} must be finite.")
        if self.standard_uncertainty is not None:
            if not np.isfinite(self.standard_uncertainty) or self.standard_uncertainty < 0.0:
                raise ValueError(
                    "RefinedParameter.standard_uncertainty must be finite and non-negative."
                )
        if not self.refined and self.standard_uncertainty is not None:
            raise ValueError("A fixed parameter cannot carry a standard uncertainty.")

    @property
    def shift(self) -> float:
        """Return how far the parameter moved from its starting value."""

        return float(self.value - self.initial_value)

    def format(self) -> str:
        """Return the value with its uncertainty in crystallographic notation.

        The uncertainty is written in parentheses on the last quoted digits,
        which is the IUCr convention: ``3.5241(6)`` means ``3.5241 +/- 0.0006``.
        """

        if self.standard_uncertainty is None or self.standard_uncertainty <= 0.0:
            return f"{self.value:.6g}"
        magnitude = int(np.floor(np.log10(self.standard_uncertainty)))
        digits = max(0, -magnitude + 1)
        scaled = round(self.standard_uncertainty * 10.0**digits)
        return f"{self.value:.{digits}f}({scaled})"


@dataclass(frozen=True, slots=True)
class RietveldResult:
    """The outcome of a whole-profile refinement, with everything needed to judge it.

    Purpose
    -------
    Carry the fitted profile, the refined parameters with uncertainties, the
    agreement indices, and enough context that the fit can be replotted,
    re-examined, and argued about without re-running it.

    Attributes
    ----------
    measured : MeasuredPowderPattern
    phase : Phase
        The phase *as refined* -- its lattice already carries the refined
        dilation, so the returned phase and the returned profile agree.
    radiation : RadiationSpec
    two_theta_deg, observed_intensity, calculated_intensity : np.ndarray
        The fitted region: measured angles, measured intensity, and the full
        calculated profile including background.
    background_intensity : np.ndarray
        The refined background alone, so the Bragg contribution can be seen.
    residual_intensity : np.ndarray
        ``observed - calculated``. Structure in this curve is the most
        informative single output of a refinement.
    parameters : tuple of RefinedParameter
    profile_r_factor, weighted_profile_r_factor, expected_r_factor : float
        ``R_p``, ``R_wp`` and ``R_exp`` as fractions, not percentages.
    goodness_of_fit : float
        ``R_wp / R_exp``. Its square is the reduced chi-squared.
    bragg_r_factor : float
        ``R_B`` from Rietveld-partitioned integrated intensities.
    durbin_watson : float
        Serial correlation of the weighted residuals; near 2 is uncorrelated,
        and values well below 1 mean the misfit is systematic.
    weight_model : str
    converged : bool
    function_evaluations : int
    """

    measured: MeasuredPowderPattern
    phase: Phase
    radiation: RadiationSpec
    two_theta_deg: np.ndarray
    observed_intensity: np.ndarray
    calculated_intensity: np.ndarray
    background_intensity: np.ndarray
    residual_intensity: np.ndarray
    reflections: tuple[PowderReflection, ...]
    parameters: tuple[RefinedParameter, ...]
    profile_r_factor: float
    weighted_profile_r_factor: float
    expected_r_factor: float
    goodness_of_fit: float
    bragg_r_factor: float
    durbin_watson: float
    weight_model: WeightModel
    converged: bool
    function_evaluations: int

    def __post_init__(self) -> None:
        arrays = tuple(
            as_float_array(value, shape=(None,))
            for value in (
                self.two_theta_deg,
                self.observed_intensity,
                self.calculated_intensity,
                self.background_intensity,
                self.residual_intensity,
            )
        )
        if arrays[0].size < 2 or len({array.shape for array in arrays}) != 1:
            raise ValueError("Rietveld profile arrays must align and hold at least two points.")
        if np.any(~np.isfinite(np.concatenate(arrays))):
            raise ValueError("Rietveld profile arrays must be finite.")
        if np.any(np.diff(arrays[0]) <= 0.0):
            raise ValueError("Rietveld profile angles must be strictly increasing.")
        if not np.allclose(arrays[4], arrays[1] - arrays[2], rtol=1e-10, atol=1e-10):
            raise ValueError("Residual intensity must equal observed minus calculated intensity.")
        for name in (
            "profile_r_factor",
            "weighted_profile_r_factor",
            "expected_r_factor",
            "goodness_of_fit",
            "bragg_r_factor",
            "durbin_watson",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"RietveldResult.{name} must be finite and non-negative.")
        if not self.parameters:
            raise ValueError("A refinement must report at least one parameter.")
        if self.weight_model not in {"unit", "poisson", "inverse_variance"}:
            raise ValueError("RietveldResult.weight_model is not recognized.")
        object.__setattr__(self, "parameters", tuple(self.parameters))
        object.__setattr__(self, "reflections", tuple(self.reflections))
        for index, name in enumerate(
            (
                "two_theta_deg",
                "observed_intensity",
                "calculated_intensity",
                "background_intensity",
                "residual_intensity",
            )
        ):
            object.__setattr__(self, name, arrays[index])

    @property
    def point_count(self) -> int:
        """Return the number of profile points included in the fit."""

        return int(self.two_theta_deg.size)

    @property
    def refined_parameter_count(self) -> int:
        """Return how many parameters were actually varied."""

        return sum(1 for parameter in self.parameters if parameter.refined)

    @property
    def reduced_chi_squared(self) -> float:
        """Return the reduced chi-squared, the square of the goodness of fit."""

        return float(self.goodness_of_fit**2)

    def parameter(self, name: str) -> RefinedParameter:
        """Return one refinement parameter by name.

        Raises
        ------
        KeyError
            If no parameter of that name took part in the refinement.
        """

        for parameter in self.parameters:
            if parameter.name == name:
                return parameter
        raise KeyError(
            f"'{name}' was not part of this refinement. Present: "
            f"{', '.join(item.name for item in self.parameters)}."
        )

    def as_pattern(self) -> PowderPattern:
        """Return the calculated Bragg profile as a :class:`PowderPattern`.

        The background is *removed*, because ``PowderPattern`` means a
        simulated diffraction profile and a fitted background is not part of
        one. The result is normalized like any simulated pattern, so it plots
        against other simulations rather than against the measured counts.
        """

        bragg = np.clip(self.calculated_intensity - self.background_intensity, 0.0, None)
        peak = float(np.max(bragg))
        if peak > 0.0:
            bragg = bragg / peak
        return PowderPattern(
            phase=self.phase,
            radiation=self.radiation,
            reflections=self.reflections,
            two_theta_grid_deg=self.two_theta_deg,
            intensity_grid=bragg,
            provenance=self.measured.provenance,
        )

    def describe(self) -> str:
        """Return the model, the refined values, the indices, and how to read them."""

        refined = [parameter for parameter in self.parameters if parameter.refined]
        fixed = [parameter for parameter in self.parameters if not parameter.refined]
        refined_text = (
            "; ".join(
                f"{parameter.name} = {parameter.format()}"
                + (f" {parameter.units}" if parameter.units else "")
                for parameter in refined
            )
            or "none"
        )
        fixed_text = (
            ", ".join(f"{parameter.name} = {parameter.value:.6g}" for parameter in fixed) or "none"
        )
        # Toby (2006): R_wp is bounded below by R_exp, so their ratio is the
        # statistic worth quoting; a GoF near 1 means the model explains the
        # data down to the counting noise and no further.
        if self.goodness_of_fit < 1.0:
            verdict = (
                "The goodness of fit is below 1, which means the model is following noise: "
                "either the weights are too pessimistic or too many parameters are being refined."
            )
        elif self.goodness_of_fit <= 2.0:
            verdict = (
                "A goodness of fit at or below about 2 is the conventional sign of an acceptable "
                "profile fit for laboratory data."
            )
        elif self.goodness_of_fit <= 5.0:
            verdict = (
                "The goodness of fit indicates real, unmodelled misfit -- most often peak shape, "
                "background curvature, texture, or a second phase."
            )
        else:
            verdict = (
                "The goodness of fit is far above the counting noise, so the calculated pattern "
                "does not describe these data; check the phase identity before reading any "
                "refined value."
            )
        correlation = (
            f"The Durbin-Watson statistic is {self.durbin_watson:.4g}; values well below 1 mean "
            "neighbouring residuals share sign, so the misfit is systematic rather than random "
            f"({_CITATION_HILL_FLACK})"
            if self.durbin_watson < 1.0
            else (
                f"The Durbin-Watson statistic is {self.durbin_watson:.4g}, consistent with "
                f"largely uncorrelated residuals ({_CITATION_HILL_FLACK})"
            )
        )
        convergence = (
            "The optimizer converged"
            if self.converged
            else (
                "The optimizer stopped without converging, so the reported values are wherever "
                "it happened to halt and must not be quoted"
            )
        )
        weights = {
            "poisson": "Poisson counting weights w = 1 / max(y_obs, 1)",
            "inverse_variance": "inverse-variance weights from the measured standard uncertainties",
            "unit": "unit weights, because the intensities carry no counting statistics",
        }[self.weight_model]
        return (
            f"Rietveld refinement of {self.phase.name} against '{self.measured.name}' over "
            f"{self.point_count} points from {self.two_theta_deg[0]:.4f} to "
            f"{self.two_theta_deg[-1]:.4f} degrees 2*theta with {self.radiation.name} radiation "
            f"({_CITATION_RIETVELD}). {convergence} after {self.function_evaluations} function "
            f"evaluations, varying {self.refined_parameter_count} parameters against "
            f"{self.point_count} observations using {weights}. "
            f"Refined: {refined_text}. Fixed: {fixed_text}. "
            f"R_p = {100.0 * self.profile_r_factor:.3f}%, R_wp = "
            f"{100.0 * self.weighted_profile_r_factor:.3f}%, R_exp = "
            f"{100.0 * self.expected_r_factor:.3f}%, R_Bragg = "
            f"{100.0 * self.bragg_r_factor:.3f}%, goodness of fit = "
            f"{self.goodness_of_fit:.4g} (reduced chi-squared "
            f"{self.reduced_chi_squared:.4g}). {verdict} {correlation}. "
            f"R factors follow the definitions in {_CITATION_TOBY} "
            "This refinement varies scale, cell dilation, zero point, profile width, background "
            "and texture strength only. It does not refine atomic coordinates, occupancies or "
            "anisotropic displacement parameters, so it tests and adjusts a structural model "
            "rather than determining a structure."
        )


def _weights(measured: MeasuredPowderPattern, mask: np.ndarray) -> tuple[np.ndarray, WeightModel]:
    observed = measured.intensity[mask]
    if measured.standard_uncertainty is not None:
        sigma = measured.standard_uncertainty[mask]
        return 1.0 / np.square(np.clip(sigma, 1e-12, None)), "inverse_variance"
    if measured.intensity_unit in {"counts", "counts_per_second"}:
        # Poisson counting statistics: sigma^2 = y. The clip at one keeps an
        # empty channel from carrying infinite weight, which is the standard
        # treatment and the reason a Rietveld fit is dominated by peak tops
        # rather than by the background.
        return 1.0 / np.clip(observed, 1.0, None), "poisson"
    return np.ones_like(observed), "unit"


def _scaled_phase(phase: Phase, lattice_scale: float) -> Phase:
    """Return ``phase`` with its cell dilated isotropically by ``lattice_scale``.

    The unit cell carries its own copy of the lattice and ``Phase`` requires the
    two to agree, so both are replaced together; dilating only the phase lattice
    raises rather than producing a phase whose sites sit in a different cell
    from its metric.
    """

    if lattice_scale == 1.0:
        return phase
    lattice = replace(
        phase.lattice,
        a=phase.lattice.a * lattice_scale,
        b=phase.lattice.b * lattice_scale,
        c=phase.lattice.c * lattice_scale,
    )
    unit_cell = None if phase.unit_cell is None else replace(phase.unit_cell, lattice=lattice)
    return replace(phase, lattice=lattice, unit_cell=unit_cell)


@dataclass(frozen=True, slots=True)
class _ForwardModel:
    """The calculation a refinement step performs, held apart from the optimizer.

    Keeping the physics in one callable object is what lets the same code make
    the starting profile, every trial profile, and the final reported profile,
    so none of the three can drift from the others.
    """

    phase: Phase
    radiation: RadiationSpec
    axis: np.ndarray
    reduced_axis: np.ndarray
    hkls: np.ndarray
    multiplicities: np.ndarray
    intensity_model: str
    preferred_plane: tuple[int, int, int] | None
    background_degree: int

    def reflection_intensities(
        self, values: dict[str, float]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return per-family 2*theta, d-spacing and integrated intensity."""

        phase = _scaled_phase(self.phase, values["lattice_scale"])
        reciprocal = phase.lattice.reciprocal_basis().matrix
        g_cartesian = self.hkls.astype(np.float64) @ reciprocal.T
        g_magnitude = np.linalg.norm(g_cartesian, axis=1)
        d_spacing = 1.0 / g_magnitude
        argument = self.radiation.wavelength_angstrom / (2.0 * d_spacing)
        visible = argument <= 1.0
        two_theta = np.full_like(argument, np.nan)
        two_theta[visible] = np.rad2deg(2.0 * np.arcsin(argument[visible]))
        from pytex.diffraction.xrd import _lorentz_polarization, _structure_factors_xray

        structure_factors = (
            np.ones(self.hkls.shape[0], dtype=np.complex128)
            if self.intensity_model == "unit"
            else _structure_factors_xray(
                phase, self.hkls, tabulated=self.intensity_model == "xray_tabulated"
            )
        )
        amplitude_squared = np.abs(structure_factors) ** 2
        # One overall isotropic displacement parameter multiplies every
        # structure factor by exp(-B sin^2(theta) / lambda^2) = exp(-B |g|^2 / 4),
        # so it is applied here rather than inside the structure-factor sum.
        b_iso = values.get("b_iso_overall", 0.0)
        if b_iso > 0.0:
            amplitude_squared = amplitude_squared * np.exp(-0.5 * b_iso * g_magnitude**2)
        lorentz = np.array(
            [
                _lorentz_polarization(np.deg2rad(angle)) if np.isfinite(angle) else 0.0
                for angle in two_theta
            ],
            dtype=np.float64,
        )
        intensities = self.multiplicities * amplitude_squared * lorentz
        march = values.get("march_coefficient")
        if march is not None and self.preferred_plane is not None:
            model = MarchDollaseModel(
                preferred_orientation=MillerPlane(
                    indices=np.asarray(self.preferred_plane, dtype=np.int64), phase=phase
                ),
                march_coefficient=march,
            )
            planes = [MillerPlane(indices=hkl, phase=phase) for hkl in self.hkls]
            intensities = intensities * np.asarray(model.factors(planes), dtype=np.float64)
        intensities = np.where(visible, intensities, 0.0)
        return two_theta, d_spacing, intensities

    def peak_matrix(self, values: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the per-reflection profile matrix, positions and intensities.

        The matrix has one row per reflection family and one column per measured
        angle. It is what makes the Rietveld intensity partitioning behind
        ``R_Bragg`` possible, so it is computed rather than the summed profile
        alone.
        """

        two_theta, _d_spacing, intensities = self.reflection_intensities(values)
        instrument = InstrumentBroadening(
            caglioti_u=values["caglioti_u"],
            caglioti_v=values["caglioti_v"],
            caglioti_w=values["caglioti_w"],
            lorentzian_y=values["lorentzian_y"],
        )
        centres = two_theta + values["zero_shift_deg"]
        matrix = np.zeros((self.hkls.shape[0], self.axis.size), dtype=np.float64)
        lines: list[tuple[np.ndarray, float]] = [(centres, 1.0)]
        if self.radiation.kalpha2_wavelength_angstrom is not None:
            # The Ka2 line diffracts from the same planes, so its position
            # follows from the same d-spacings and only the wavelength changes.
            ratio = self.radiation.kalpha2_wavelength_angstrom / (
                self.radiation.wavelength_angstrom
            )
            argument = np.clip(np.sin(np.deg2rad(0.5 * two_theta)) * ratio, -1.0, 1.0)
            kalpha2 = np.rad2deg(2.0 * np.arcsin(argument)) + values["zero_shift_deg"]
            lines.append((kalpha2, self.radiation.kalpha2_relative_intensity))
        for index in range(self.hkls.shape[0]):
            if not np.isfinite(two_theta[index]) or intensities[index] <= 0.0:
                continue
            gaussian = float(instrument.gaussian_fwhm_deg(two_theta[index]))
            lorentzian = float(instrument.lorentzian_fwhm_deg(two_theta[index]))
            width = float(_tch_width(np.array(gaussian), np.array(lorentzian)))
            eta = float(_tch_eta(np.array(gaussian), np.array(lorentzian)))
            for centre, weight in lines:
                matrix[index] += (
                    weight
                    * intensities[index]
                    * _pseudo_voigt(self.axis, float(centre[index]), width, eta)
                )
        return matrix, centres, intensities

    def background(self, coefficients: np.ndarray) -> np.ndarray:
        """Return the Chebyshev background on the fitted angular support."""

        if coefficients.size == 0:
            return np.zeros_like(self.axis)
        values = np.polynomial.chebyshev.chebval(self.reduced_axis, coefficients)
        return np.asarray(values, dtype=np.float64)

    def profile(self, values: dict[str, float], coefficients: np.ndarray) -> np.ndarray:
        matrix, _, _ = self.peak_matrix(values)
        profile = values["scale"] * matrix.sum(axis=0) + self.background(coefficients)
        return np.asarray(profile, dtype=np.float64)


def _pseudo_voigt(axis: np.ndarray, centre: float, fwhm: float, eta: float) -> np.ndarray:
    """Area-normalized pseudo-Voigt of full width ``fwhm`` centred at ``centre``."""

    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gaussian = np.exp(-0.5 * ((axis - centre) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))
    gamma = 0.5 * fwhm
    lorentzian = gamma / (np.pi * ((axis - centre) ** 2 + gamma * gamma))
    return np.asarray(eta * lorentzian + (1.0 - eta) * gaussian, dtype=np.float64)


def refine_rietveld(
    measured: MeasuredPowderPattern,
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    instrument: InstrumentBroadening | None = None,
    refine: Sequence[str] = DEFAULT_REFINEMENT_SET,
    background_degree: int = 4,
    refine_background: bool = True,
    preferred_orientation_plane: tuple[int, int, int] | None = None,
    march_coefficient: float = 1.0,
    b_iso_overall: float = 0.0,
    two_theta_range_deg: tuple[float, float] | None = None,
    max_index: int = 6,
    intensity_model: Literal["xray_atomic_number", "xray_tabulated", "unit"] = "xray_tabulated",
    max_function_evaluations: int = 400,
) -> RietveldResult:
    """Refine a structural and profile model against a measured powder pattern.

    Purpose
    -------
    Fit the whole measured profile -- not extracted peak intensities -- with a
    calculated pattern from ``phase``, adjusting the scale, unit-cell dilation,
    instrument zero, peak widths, texture strength and background until the
    calculated and measured profiles agree, and report how well they do.

    Method
    ------
    Reflection families are enumerated once from the starting cell, over an
    angular window padded so that a shifting cell cannot move a reflection into
    or out of the fit and change the model discontinuously. Each least-squares
    step recomputes positions from the dilated cell, structure factors with the
    current displacement parameter, March-Dollase texture factors, and
    Thompson-Cox-Hastings pseudo-Voigt profiles, and adds a Chebyshev
    background. Trust-region least squares (``scipy.optimize.least_squares``)
    minimizes the weighted residual with every parameter bounded, since an
    unbounded profile refinement will happily reach a lower ``R_wp`` at a
    physically impossible cell.

    Parameters
    ----------
    measured
        The observed profile. Do **not** subtract a background first: the
        background is refined jointly here, and subtracting it beforehand
        removes the correlation between background and scale that the
        uncertainties depend on.
    phase
        The structural model. It needs an atomic basis for structure-sensitive
        intensities; without one only geometry and multiplicity contribute.
    radiation
        Defaults to the measured pattern's own radiation when it carries one,
        and to Cu K-alpha otherwise.
    instrument
        Starting resolution function. Its coefficients are the starting values
        for any width parameter in ``refine``. Defaults to a constant
        0.1-degree width.
    refine
        Names from :data:`REFINABLE_PARAMETERS` to vary. Everything else is
        held at its starting value and reported as fixed.
    background_degree
        Chebyshev degree for the background. Degree 4 to 8 covers most
        laboratory patterns; a high degree will absorb genuine broad features,
        including an amorphous halo you may want to see.
    refine_background
        Vary the background coefficients. Fixing them at a flat median is
        occasionally useful for teaching what the background does, and is
        almost never right for real data.
    preferred_orientation_plane
        The ``(hkl)`` whose normals cluster along the specimen axis. Required
        before ``march_coefficient`` may be refined.
    march_coefficient
        Starting texture strength; 1 is a random powder.
    b_iso_overall
        Starting overall isotropic displacement parameter in square angstroms.
    two_theta_range_deg
        Restrict the fit to this window. Defaults to the whole measured range.
    max_index
        Largest absolute Miller index enumerated.
    intensity_model
        Structure-factor model, as in
        :func:`~pytex.diffraction.xrd.generate_powder_reflections`.
    max_function_evaluations
        Optimizer budget.

    Returns
    -------
    RietveldResult
        Fitted profile, refined parameters with uncertainties, ``R_p``,
        ``R_wp``, ``R_exp``, ``R_Bragg``, goodness of fit, Durbin-Watson, and
        an explainable summary.

    Raises
    ------
    ValueError
        If a name in ``refine`` is not refinable, if ``march_coefficient`` is
        refined without a preferred-orientation plane, if the fitted window
        holds fewer points than parameters, or if no reflection of the phase
        falls inside the window -- which means the phase and the pattern do not
        belong together, and is a result worth stopping on rather than fitting
        through.

    See Also
    --------
    pytex.diffraction.xrd_measurement.compare_powder_patterns : A scaled
        comparison with no refinement, for phase identification.
    pytex.diffraction.xrd_background.estimate_background : A standalone
        background estimate, for inspection rather than for refinement.
    """

    unknown = tuple(name for name in refine if name not in REFINABLE_PARAMETERS)
    if unknown:
        raise ValueError(
            f"Unknown refinement parameters {unknown}. Refinable: {REFINABLE_PARAMETERS}."
        )
    refine_set = set(refine) | {"scale"}
    if "march_coefficient" in refine_set and preferred_orientation_plane is None:
        raise ValueError(
            "Refining march_coefficient needs preferred_orientation_plane: a texture strength "
            "about an unstated axis has no meaning."
        )
    if background_degree < 0:
        raise ValueError("background_degree must be non-negative.")
    radiation_spec = radiation or measured.radiation or RadiationSpec.cu_ka()
    resolution = instrument or InstrumentBroadening.ideal(0.1)

    axis_full = measured.two_theta_deg
    if two_theta_range_deg is None:
        mask = np.ones_like(axis_full, dtype=bool)
    else:
        lower, upper = (float(value) for value in two_theta_range_deg)
        if lower >= upper:
            raise ValueError("two_theta_range_deg must satisfy min < max.")
        mask = (axis_full >= lower) & (axis_full <= upper)
    axis = axis_full[mask]
    observed = measured.intensity[mask]
    if axis.size < 10:
        raise ValueError("A refinement needs at least ten measured points in the fitted window.")

    # Enumerate once, over a padded window: the cell may dilate by up to the
    # lattice_scale bound, and a reflection that enters or leaves the model
    # mid-refinement makes the residual discontinuous and the optimizer's
    # derivative estimates wrong.
    padding = 3.0
    enumeration_range = (
        max(0.0, float(axis[0]) - padding),
        min(180.0, float(axis[-1]) + padding),
    )
    families = generate_powder_reflections(
        phase,
        radiation=radiation_spec,
        two_theta_range_deg=enumeration_range,
        max_index=max_index,
        intensity_model=intensity_model,
    )
    if not families:
        raise ValueError(
            f"No reflection of {phase.name} falls between {enumeration_range[0]:.2f} and "
            f"{enumeration_range[1]:.2f} degrees 2*theta at this wavelength. The phase and the "
            "measurement do not correspond; refining would fit background to noise."
        )
    hkls = np.array([reflection.miller_indices for reflection in families], dtype=np.int64)
    multiplicities = np.array(
        [reflection.multiplicity for reflection in families], dtype=np.float64
    )
    span = float(axis[-1] - axis[0])
    model = _ForwardModel(
        phase=phase,
        radiation=radiation_spec,
        axis=axis,
        reduced_axis=2.0 * (axis - axis[0]) / span - 1.0,
        hkls=hkls,
        multiplicities=multiplicities,
        intensity_model=intensity_model,
        preferred_plane=preferred_orientation_plane,
        background_degree=background_degree,
    )

    starting: dict[str, float] = {
        "scale": 1.0,
        "zero_shift_deg": resolution.zero_shift_deg,
        "lattice_scale": 1.0,
        "caglioti_u": resolution.caglioti_u,
        "caglioti_v": resolution.caglioti_v,
        "caglioti_w": resolution.caglioti_w,
        "lorentzian_y": resolution.lorentzian_y,
        "b_iso_overall": float(b_iso_overall),
        "march_coefficient": float(march_coefficient),
    }
    # The starting scale is set from the data rather than left at one, because a
    # trust-region step sized for a scale of 1 is meaningless when the true
    # scale is 10^4 counts, and the refinement would spend its budget crawling
    # towards the right order of magnitude.
    trial = model.profile(starting, np.zeros(0))
    trial_peak = float(np.max(trial))
    if trial_peak <= 0.0:
        raise ValueError(
            "The starting model calculates no intensity in the fitted window; check the phase's "
            "atomic basis and the radiation."
        )
    starting["scale"] = float(np.max(observed)) / trial_peak

    background_names = tuple(f"background_{index}" for index in range(background_degree + 1))
    background_start = np.zeros(background_degree + 1)
    # Chebyshev T_0 is the constant, so the median observed level is the honest
    # starting background and every higher coefficient starts at zero.
    background_start[0] = float(np.median(observed))

    if "march_coefficient" in refine_set and march_coefficient == 1.0:
        # r = 1 is the untextured point, and it is a bad place to start from:
        # plate-like (r < 1) and needle-like (r > 1) textures lie on opposite
        # sides of it, so a local optimizer commits to whichever side its first
        # derivative points and can settle in the wrong one. A coarse scan costs
        # a handful of profile evaluations and picks the side from the data
        # instead. It runs only when the caller left the default, because a
        # stated starting value is information the scan would discard.
        starting["march_coefficient"] = _best_march_start(
            model, starting, background_start, observed, np.sqrt(_weights(measured, mask)[0])
        )

    varied = [name for name in REFINABLE_PARAMETERS if name in refine_set]
    weights, weight_model = _weights(measured, mask)
    root_weight = np.sqrt(weights)

    def unpack(vector: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
        values = dict(starting)
        for index, name in enumerate(varied):
            values[name] = float(vector[index])
        coefficients = vector[len(varied) :] if refine_background else background_start
        return values, np.asarray(coefficients, dtype=np.float64)

    def residual(vector: np.ndarray) -> np.ndarray:
        values, coefficients = unpack(vector)
        calculated = model.profile(values, coefficients)
        return np.asarray(root_weight * (observed - calculated), dtype=np.float64)

    initial_vector = np.concatenate(
        (
            np.array([starting[name] for name in varied], dtype=np.float64),
            background_start if refine_background else np.zeros(0),
        )
    )
    lower_bounds = np.concatenate(
        (
            np.array([_PARAMETER_BOUNDS[name][0] for name in varied], dtype=np.float64),
            np.full(background_start.size, -np.inf) if refine_background else np.zeros(0),
        )
    )
    upper_bounds = np.concatenate(
        (
            np.array([_PARAMETER_BOUNDS[name][1] for name in varied], dtype=np.float64),
            np.full(background_start.size, np.inf) if refine_background else np.zeros(0),
        )
    )
    parameter_count = int(initial_vector.size)
    if axis.size <= parameter_count:
        raise ValueError(
            f"The fitted window holds {axis.size} points but the model varies {parameter_count} "
            "parameters. Widen the window or refine fewer parameters."
        )

    solution = least_squares(
        residual,
        initial_vector,
        bounds=(lower_bounds, upper_bounds),
        max_nfev=max_function_evaluations,
        xtol=1e-12,
        ftol=1e-12,
    )
    values, coefficients = unpack(solution.x)
    matrix, centres, integrated = model.peak_matrix(values)
    background = model.background(coefficients)
    calculated = values["scale"] * matrix.sum(axis=0) + background
    residual_profile = observed - calculated

    denominator = float(np.sum(observed))
    rp = float(np.sum(np.abs(residual_profile)) / denominator) if denominator > 0.0 else np.inf
    weighted_denominator = float(np.sum(weights * np.square(observed)))
    rwp = float(np.sqrt(np.sum(weights * np.square(residual_profile)) / weighted_denominator))
    degrees_of_freedom = max(axis.size - parameter_count, 1)
    rexp = float(np.sqrt(degrees_of_freedom / weighted_denominator))
    goodness = rwp / rexp if rexp > 0.0 else np.inf

    # Rietveld intensity partitioning: the observed intensity at each point is
    # shared between reflections in proportion to their calculated contribution
    # there. It is a model-dependent quantity by construction, which is why
    # R_Bragg is quoted alongside R_wp rather than instead of it.
    scaled_matrix = values["scale"] * matrix
    total = scaled_matrix.sum(axis=0)
    share = np.divide(scaled_matrix, total, out=np.zeros_like(scaled_matrix), where=total > 0.0)
    observed_bragg = share @ np.clip(observed - background, 0.0, None)
    calculated_bragg = scaled_matrix.sum(axis=1)
    bragg_denominator = float(np.sum(observed_bragg))
    bragg_r = (
        float(np.sum(np.abs(observed_bragg - calculated_bragg)) / bragg_denominator)
        if bragg_denominator > 0.0
        else 0.0
    )

    weighted_residual = root_weight * residual_profile
    residual_energy = float(np.sum(np.square(weighted_residual)))
    durbin_watson = (
        float(np.sum(np.square(np.diff(weighted_residual))) / residual_energy)
        if residual_energy > 0.0
        else 2.0
    )

    uncertainties = _standard_uncertainties(solution.jac, residual_energy, degrees_of_freedom)

    parameters: list[RefinedParameter] = []
    for name in REFINABLE_PARAMETERS:
        if name == "march_coefficient" and preferred_orientation_plane is None:
            continue
        is_refined = name in refine_set
        index = varied.index(name) if is_refined else None
        parameters.append(
            RefinedParameter(
                name=name,
                value=values[name],
                initial_value=starting[name] if name != "scale" else starting["scale"],
                standard_uncertainty=(None if index is None else float(uncertainties[index])),
                refined=is_refined,
                units=_PARAMETER_UNITS[name],
                description=_PARAMETER_DESCRIPTIONS[name],
            )
        )
    for offset, name in enumerate(background_names):
        varied_index = len(varied) + offset if refine_background else None
        parameters.append(
            RefinedParameter(
                name=name,
                value=float(coefficients[offset]),
                initial_value=float(background_start[offset]),
                standard_uncertainty=(
                    None if varied_index is None else float(uncertainties[varied_index])
                ),
                refined=refine_background,
                units="intensity" if offset == 0 else "",
                description=(
                    f"Chebyshev background coefficient of order {offset} on the fitted window."
                ),
            )
        )

    refined_phase = _scaled_phase(phase, values["lattice_scale"])
    reflections = _result_reflections(
        refined_phase,
        model,
        values,
        centres,
        integrated,
        window=(float(axis[0]), float(axis[-1])),
        intensity_model=intensity_model,
    )
    return RietveldResult(
        measured=measured,
        phase=refined_phase,
        radiation=radiation_spec,
        two_theta_deg=axis,
        observed_intensity=observed,
        calculated_intensity=calculated,
        background_intensity=background,
        residual_intensity=residual_profile,
        reflections=reflections,
        parameters=tuple(parameters),
        profile_r_factor=rp,
        weighted_profile_r_factor=rwp,
        expected_r_factor=rexp,
        goodness_of_fit=float(goodness),
        bragg_r_factor=bragg_r,
        durbin_watson=durbin_watson,
        weight_model=weight_model,
        converged=bool(solution.success),
        function_evaluations=int(solution.nfev),
    )


def _best_march_start(
    model: _ForwardModel,
    starting: dict[str, float],
    background_start: np.ndarray,
    observed: np.ndarray,
    root_weight: np.ndarray,
) -> float:
    """Return the March coefficient whose profile best matches before refining.

    The scale is re-optimized analytically at each trial coefficient, because
    otherwise the scan would compare texture strengths at the wrong amplitude
    and simply pick whichever one happens to suit the fixed scale.
    """

    background = model.background(background_start)
    residual_target = observed - background
    best_coefficient, best_cost = starting["march_coefficient"], np.inf
    for coefficient in (0.4, 0.55, 0.7, 0.85, 1.0, 1.2, 1.5, 1.9, 2.4):
        trial = dict(starting)
        trial["march_coefficient"] = coefficient
        matrix, _, _ = model.peak_matrix(trial)
        calculated = matrix.sum(axis=0)
        denominator = float(np.sum(root_weight**2 * calculated**2))
        if denominator <= 0.0:
            continue
        scale = float(np.sum(root_weight**2 * calculated * residual_target) / denominator)
        cost = float(np.sum((root_weight * (residual_target - scale * calculated)) ** 2))
        if cost < best_cost:
            best_coefficient, best_cost = coefficient, cost
    return best_coefficient


def _standard_uncertainties(
    jacobian: np.ndarray, residual_energy: float, degrees_of_freedom: int
) -> np.ndarray:
    """Return one estimated standard deviation per varied parameter.

    The covariance is ``(J^T J)^-1`` scaled by the reduced chi-squared, which is
    the standard treatment when the weights are only known up to a constant. A
    rank-deficient Jacobian -- two parameters the data cannot tell apart -- gives
    an infinite uncertainty rather than a small one, so a correlated pair is
    visible instead of being quietly reported as well determined.
    """

    try:
        _, singular_values, right = np.linalg.svd(jacobian, full_matrices=False)
    except np.linalg.LinAlgError:
        return np.full(jacobian.shape[1], np.inf)
    threshold = (
        np.finfo(float).eps
        * max(jacobian.shape)
        * (singular_values[0] if singular_values.size else 0.0)
    )
    usable = singular_values > threshold
    inverse = np.zeros_like(singular_values)
    inverse[usable] = 1.0 / singular_values[usable] ** 2
    covariance_diagonal = np.sum((right.T**2) * inverse, axis=1)
    covariance_diagonal[~np.isfinite(covariance_diagonal)] = np.inf
    if not np.all(usable):
        covariance_diagonal = np.where(
            np.isclose(covariance_diagonal, 0.0), np.inf, covariance_diagonal
        )
    scale = residual_energy / max(degrees_of_freedom, 1)
    sigmas = np.sqrt(np.clip(covariance_diagonal * scale, 0.0, None))
    return np.asarray(sigmas, dtype=np.float64)


def _result_reflections(
    phase: Phase,
    model: _ForwardModel,
    values: dict[str, float],
    centres: np.ndarray,
    integrated: np.ndarray,
    *,
    window: tuple[float, float],
    intensity_model: str,
) -> tuple[PowderReflection, ...]:
    """Build the reported reflection list from the refined model."""

    _, d_spacing, _ = model.reflection_intensities(values)
    from pytex.diffraction.xrd import _lorentz_polarization, _structure_factors_xray

    structure_factors = (
        np.ones(model.hkls.shape[0], dtype=np.complex128)
        if intensity_model == "unit"
        else _structure_factors_xray(
            phase, model.hkls, tabulated=intensity_model == "xray_tabulated"
        )
    )
    reflections: list[PowderReflection] = []
    for index in range(model.hkls.shape[0]):
        centre = float(centres[index])
        if not np.isfinite(centre) or not window[0] <= centre <= window[1]:
            continue
        if integrated[index] <= 0.0:
            continue
        structure_factor = complex(structure_factors[index])
        reflections.append(
            PowderReflection(
                miller_indices=model.hkls[index],
                d_spacing_angstrom=float(d_spacing[index]),
                two_theta_deg=centre,
                intensity=float(integrated[index]),
                structure_factor_amplitude=float(abs(structure_factor)),
                multiplicity=int(model.multiplicities[index]),
                structure_factor_real=float(structure_factor.real),
                structure_factor_imag=float(structure_factor.imag),
                lorentz_polarization_factor=_lorentz_polarization(np.deg2rad(centre)),
                intensity_model=intensity_model,
            )
        )
    reflections.sort(key=lambda reflection: reflection.two_theta_deg)
    return tuple(reflections)


__all__ = [
    "DEFAULT_REFINEMENT_SET",
    "REFINABLE_PARAMETERS",
    "RIETVELD_REFINEMENT_SCHEMA",
    "RefinedParameter",
    "RietveldResult",
    "WeightModel",
    "refine_rietveld",
]
