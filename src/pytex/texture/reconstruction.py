from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core.lattice import CrystalPlane
from pytex.core.provenance import ProvenanceRecord
from pytex.core.sphere import directions_to_spherical_angles
from pytex.texture.harmonics import HarmonicODF, HarmonicODFReconstructionReport
from pytex.texture.models import (
    ODF,
    KernelSpec,
    ODFInversionReport,
    PoleFigure,
    PoleFigureDifference,
    random_pole_density,
)

CorrectionPolicy = Literal["clip_zero", "raise"]
ReconstructionAlgorithm = Literal["discrete", "harmonic"]
DefocusReducer = Literal["mean", "median"]

POLE_FIGURE_DEFOCUS_CALIBRATION_SCHEMA = "pytex.texture.pole_figure_defocus_calibration"


@dataclass(frozen=True, slots=True)
class PoleFigureDefocusCalibration:
    """A radial defocusing curve measured from a same-reflection random standard.

    Use this result to evaluate auditable correction factors on a specimen pole
    figure before ODF inversion. Each factor is the background-subtracted random
    standard intensity at one tilt ring divided by the corresponding intensity
    at the lowest calibrated tilt. The curve is reflection-specific and never
    extrapolates beyond its measured range.
    """

    pole: CrystalPlane
    tilt_deg: np.ndarray
    defocus_factors: np.ndarray
    ring_intensities: np.ndarray
    azimuthal_relative_std: np.ndarray
    ring_counts: np.ndarray
    reference_tilt_deg: float
    reference_intensity: float
    background: float = 0.0
    reducer: DefocusReducer = "mean"
    synthetic: bool = False
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        tilt = as_float_array(self.tilt_deg, shape=(None,))
        factors = as_float_array(self.defocus_factors, shape=tilt.shape)
        intensities = as_float_array(self.ring_intensities, shape=tilt.shape)
        spread = as_float_array(self.azimuthal_relative_std, shape=tilt.shape)
        counts = np.asarray(self.ring_counts, dtype=np.int64).reshape(-1)
        if tilt.size < 2:
            raise ValueError("A defocusing calibration requires at least two tilt rings.")
        if np.any(~np.isfinite(tilt)) or np.any(np.diff(tilt) <= 0.0):
            raise ValueError("Calibration tilt angles must be finite and strictly increasing.")
        if tilt[0] < 0.0 or tilt[-1] > 90.0:
            raise ValueError("Reflection-geometry calibration tilts must lie in [0, 90] degrees.")
        if np.any(~np.isfinite(factors)) or np.any(factors <= 0.0):
            raise ValueError("Defocusing factors must be positive and finite.")
        if not np.isclose(factors[0], 1.0, rtol=1e-12, atol=1e-12):
            raise ValueError("The lowest-tilt defocusing factor must be normalized to one.")
        if np.any(~np.isfinite(intensities)) or np.any(intensities <= 0.0):
            raise ValueError("Background-subtracted ring intensities must be positive and finite.")
        if np.any(~np.isfinite(spread)) or np.any(spread < 0.0):
            raise ValueError("Azimuthal relative standard deviations must be non-negative.")
        if counts.shape != tilt.shape or np.any(counts <= 0):
            raise ValueError("ring_counts must contain one positive count per tilt ring.")
        if not np.isclose(self.reference_tilt_deg, tilt[0], rtol=0.0, atol=1e-12):
            raise ValueError("reference_tilt_deg must equal the lowest calibrated tilt.")
        if not np.isclose(self.reference_intensity, intensities[0], rtol=1e-12, atol=1e-12):
            raise ValueError("reference_intensity must equal the lowest-tilt ring intensity.")
        if not np.allclose(
            factors, intensities / self.reference_intensity, rtol=1e-12, atol=1e-12
        ):
            raise ValueError(
                "Defocusing factors must equal ring intensities divided by reference intensity."
            )
        if not np.isfinite(self.background) or self.background < 0.0:
            raise ValueError("Calibration background must be finite and non-negative.")
        if self.reducer not in {"mean", "median"}:
            raise ValueError("Calibration reducer must be 'mean' or 'median'.")
        counts = np.ascontiguousarray(counts)
        counts.setflags(write=False)
        object.__setattr__(self, "tilt_deg", tilt)
        object.__setattr__(self, "defocus_factors", factors)
        object.__setattr__(self, "ring_intensities", intensities)
        object.__setattr__(self, "azimuthal_relative_std", spread)
        object.__setattr__(self, "ring_counts", counts)

    @property
    def max_azimuthal_relative_std(self) -> float:
        """Return the largest within-ring standard deviation divided by its mean."""

        return float(np.max(self.azimuthal_relative_std))

    def factors_for(self, pole_figure: PoleFigure) -> np.ndarray:
        """Interpolate factors onto a same-reflection pole figure without extrapolation.

        Parameters
        ----------
        pole_figure
            Measured specimen figure. Its pole must match the random standard,
            and every effective polar angle must lie inside the calibrated range.

        Returns
        -------
        np.ndarray
            One positive read-only defocusing factor per measured direction.
        """

        if pole_figure.pole != self.pole:
            raise ValueError(
                "A defocusing calibration is reflection-specific; target and standard poles differ."
            )
        polar_deg, _ = directions_to_spherical_angles(
            pole_figure.sample_directions, antipodal=pole_figure.antipodal
        )
        tolerance = 1e-9
        if np.any(polar_deg < self.tilt_deg[0] - tolerance) or np.any(
            polar_deg > self.tilt_deg[-1] + tolerance
        ):
            raise ValueError(
                "Target pole-figure tilt lies outside the calibrated range; extrapolation is "
                "not permitted."
            )
        factors = np.interp(polar_deg, self.tilt_deg, self.defocus_factors)
        factors = np.ascontiguousarray(factors, dtype=np.float64)
        factors.setflags(write=False)
        return factors

    def correction_spec(
        self,
        pole_figure: PoleFigure,
        *,
        scale: float = 1.0,
        background: float = 0.0,
        missing_intensity_policy: CorrectionPolicy = "clip_zero",
    ) -> PoleFigureCorrectionSpec:
        """Build a correction spec evaluated on one target pole-figure support.

        ``background`` is the target specimen scan's background in its own
        intensity units; it is deliberately independent of the background used
        to calibrate the random standard. The returned spec subtracts that value
        before dividing by the interpolated factors.
        """

        return PoleFigureCorrectionSpec(
            scale=scale,
            background=background,
            defocus_factors=self.factors_for(pole_figure),
            missing_intensity_policy=missing_intensity_policy,
            provenance=self.provenance or pole_figure.provenance,
        )

    def describe(self) -> str:
        """Return the calibration source, normalization, diagnostics, and limits."""

        origin = "synthetic validation standard" if self.synthetic else "experimental standard"
        return (
            f"Random-standard pole-figure defocusing calibration from an explicitly labelled "
            f"{origin}: {len(self.tilt_deg)} tilt rings span {self.tilt_deg[0]:.3f} to "
            f"{self.tilt_deg[-1]:.3f} degrees. Background {self.background:.6g} was subtracted "
            f"before {self.reducer} ring reduction; factors are normalized to 1 at "
            f"{self.reference_tilt_deg:.3f} degrees and range from "
            f"{np.min(self.defocus_factors):.6g} to {np.max(self.defocus_factors):.6g}. "
            f"The largest azimuthal relative standard deviation is "
            f"{self.max_azimuthal_relative_std:.6g}. The curve is valid only for the same "
            "reflection and measured tilt interval; factors_for() refuses extrapolation. "
            "Target correction subtracts its declared background before dividing by this curve."
        )


def defocus_from_random_standard(
    random_standard: PoleFigure,
    *,
    background: float = 0.0,
    reducer: DefocusReducer = "mean",
    ring_tolerance_deg: float = 1e-6,
    synthetic: bool | None = None,
    provenance: ProvenanceRecord | None = None,
) -> PoleFigureDefocusCalibration:
    """Calibrate a radial defocusing curve from an untextured reference specimen.

    Use this after importing a random-standard scan for the same reflection and
    instrument configuration as the specimen scan. Intensities are first
    background-subtracted, then grouped by polar-angle ring and reduced over
    azimuth. Ring values are divided by the lowest-tilt ring value, matching the
    established experimental correction ``(I - background) / defocus``.

    Parameters
    ----------
    random_standard
        A ``sampling='sampled_density'`` pole figure whose ideal texture signal
        is azimuthally and radially constant before instrument losses.
    background
        Constant background in the random-standard intensity units.
    reducer
        Arithmetic mean (default) or median across each azimuthal ring.
    ring_tolerance_deg
        Maximum adjacent tilt difference joined into one nominal ring.
    synthetic
        Explicit data-origin label. When omitted, a ``synthetic=true``
        provenance metadata entry is honored; otherwise the data are treated as
        experimental.
    provenance
        Override for the result; the standard provenance is retained by default.

    Returns
    -------
    PoleFigureDefocusCalibration
        Radial factors, ring intensities/counts, azimuthal-scatter diagnostics,
        and the interpolation/correction helpers.
    """

    if random_standard.sampling != "sampled_density":
        raise ValueError(
            "Random-standard calibration requires sampled_density intensities, not pole weights."
        )
    if not np.isfinite(background) or background < 0.0:
        raise ValueError("Random-standard background must be finite and non-negative.")
    if reducer not in {"mean", "median"}:
        raise ValueError("Random-standard reducer must be 'mean' or 'median'.")
    if not np.isfinite(ring_tolerance_deg) or ring_tolerance_deg <= 0.0:
        raise ValueError("ring_tolerance_deg must be positive and finite.")
    polar_deg, _ = directions_to_spherical_angles(
        random_standard.sample_directions, antipodal=random_standard.antipodal
    )
    if np.any(polar_deg > 90.0 + ring_tolerance_deg):
        raise ValueError("Random-standard directions must lie on the reflection upper hemisphere.")
    corrected = np.asarray(random_standard.intensities, dtype=np.float64) - background
    if np.any(corrected <= 0.0):
        raise ValueError(
            "Random-standard intensities must remain positive after background subtraction."
        )
    order = np.argsort(polar_deg, kind="stable")
    sorted_tilt = polar_deg[order]
    sorted_intensity = corrected[order]
    group_starts = np.concatenate(
        ([0], np.flatnonzero(np.diff(sorted_tilt) > ring_tolerance_deg) + 1)
    )
    group_ends = np.concatenate((group_starts[1:], [sorted_tilt.size]))
    if group_starts.size < 2:
        raise ValueError("Random-standard calibration requires at least two distinct tilt rings.")
    tilt = np.array(
        [
            np.mean(sorted_tilt[start:end])
            for start, end in zip(group_starts, group_ends, strict=True)
        ]
    )
    ring_values = [
        sorted_intensity[start:end] for start, end in zip(group_starts, group_ends, strict=True)
    ]
    reduce = np.mean if reducer == "mean" else np.median
    ring_intensities = np.array([reduce(values) for values in ring_values], dtype=np.float64)
    ring_spread = np.array(
        [
            np.std(values, ddof=0) / mean
            for values, mean in zip(ring_values, ring_intensities, strict=True)
        ],
        dtype=np.float64,
    )
    ring_counts = np.array([values.size for values in ring_values], dtype=np.int64)
    factors = ring_intensities / ring_intensities[0]
    effective_provenance = provenance or random_standard.provenance
    if synthetic is None:
        declared = (
            effective_provenance.metadata.get("synthetic", "false")
            if effective_provenance is not None
            else "false"
        )
        synthetic = str(declared).lower() in {"1", "true", "yes"}
    return PoleFigureDefocusCalibration(
        pole=random_standard.pole,
        tilt_deg=tilt,
        defocus_factors=factors,
        ring_intensities=ring_intensities,
        azimuthal_relative_std=ring_spread,
        ring_counts=ring_counts,
        reference_tilt_deg=float(tilt[0]),
        reference_intensity=float(ring_intensities[0]),
        background=float(background),
        reducer=reducer,
        synthetic=synthetic,
        provenance=effective_provenance,
    )


@dataclass(frozen=True, slots=True)
class PoleFigureCorrectionSpec:
    """Deterministic correction metadata for pole-figure intensities."""

    scale: float = 1.0
    background: float = 0.0
    defocus_factors: np.ndarray | None = None
    missing_intensity_policy: CorrectionPolicy = "clip_zero"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.scale) or self.scale <= 0.0:
            raise ValueError("PoleFigureCorrectionSpec.scale must be positive and finite.")
        if not np.isfinite(self.background) or self.background < 0.0:
            raise ValueError("PoleFigureCorrectionSpec.background must be finite and non-negative.")
        if self.missing_intensity_policy not in {"clip_zero", "raise"}:
            raise ValueError("missing_intensity_policy must be either 'clip_zero' or 'raise'.")
        if self.defocus_factors is not None:
            factors = np.asarray(self.defocus_factors, dtype=np.float64).reshape(-1)
            if np.any(~np.isfinite(factors)) or np.any(factors <= 0.0):
                raise ValueError("defocus_factors must contain positive finite values.")
            factors = np.ascontiguousarray(factors, dtype=np.float64)
            factors.setflags(write=False)
            object.__setattr__(self, "defocus_factors", factors)

    def apply(self, pole_figure: PoleFigure) -> PoleFigure:
        """Apply background, defocusing, and scale corrections to a pole figure.

        Purpose
        -------
        Measured X-ray pole figures need instrumental correction before
        inversion: a background is subtracted, the tilt-dependent defocusing loss
        is divided out, and a scale is applied. Doing this in a declared spec —
        rather than ad hoc at the call site — keeps the correction auditable and
        reproducible.

        Order of operations
        -------------------
        Background subtraction first, then defocus division, then scaling.
        This is the experimental random-standard convention: the background
        is not part of the specimen signal and must not be amplified by the
        defocusing correction.

        Parameters
        ----------
        pole_figure : PoleFigure
            The measured figure. Defocus factors, when given, must match its
            intensity shape.

        Returns
        -------
        PoleFigure
            The corrected figure. Negative intensities are physically
            meaningless; the spec's ``missing_intensity_policy`` decides whether
            they raise or are clipped to zero, so the choice is explicit rather
            than hidden.
        """

        intensities = np.asarray(pole_figure.intensities, dtype=np.float64) - self.background
        if self.defocus_factors is not None:
            if self.defocus_factors.shape != intensities.shape:
                raise ValueError("defocus_factors must match the pole-figure intensity shape.")
            intensities = intensities / self.defocus_factors
        corrected = self.scale * intensities
        if np.any(corrected < 0.0):
            if self.missing_intensity_policy == "raise":
                raise ValueError("Pole-figure correction produced negative intensities.")
            corrected = np.maximum(corrected, 0.0)
        corrected = np.ascontiguousarray(corrected, dtype=np.float64)
        corrected.setflags(write=False)
        return PoleFigure(
            pole=pole_figure.pole,
            sample_directions=pole_figure.sample_directions,
            intensities=corrected,
            specimen_frame=pole_figure.specimen_frame,
            antipodal=pole_figure.antipodal,
            sample_symmetry=pole_figure.sample_symmetry,
            provenance=self.provenance or pole_figure.provenance,
            includes_symmetry_family=pole_figure.includes_symmetry_family,
            # Correcting intensities does not change what they mean.
            sampling=pole_figure.sampling,
        )


@dataclass(frozen=True, slots=True)
class PoleFigureResidualReport:
    """How well an ODF reproduces one measured pole figure.

    Purpose
    -------
    The goodness-of-fit check on a reconstruction. An inversion that cannot
    reproduce its own input data is not usable, whatever its internal
    convergence reported, so this comparison is the meaningful acceptance
    test.

    Attributes
    ----------
    pole_figure : PoleFigure
        The measurement compared against.
    predicted_intensities : np.ndarray
        Pole densities the ODF implies at the measured directions.
    residuals : np.ndarray
        Predicted minus measured, per direction.
    Remaining attributes summarize the residual magnitudes.
    provenance : ProvenanceRecord, optional
    """

    pole_figure: PoleFigure
    predicted_intensities: np.ndarray
    residuals: np.ndarray
    residual_norm: float
    relative_residual_norm: float
    mean_absolute_error: float
    max_absolute_error: float
    observation_count: int
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        predicted = as_float_array(self.predicted_intensities, shape=(self.observation_count,))
        residuals = as_float_array(self.residuals, shape=(self.observation_count,))
        if self.observation_count <= 0:
            raise ValueError("PoleFigureResidualReport.observation_count must be positive.")
        for name, value in (
            ("residual_norm", self.residual_norm),
            ("relative_residual_norm", self.relative_residual_norm),
            ("mean_absolute_error", self.mean_absolute_error),
            ("max_absolute_error", self.max_absolute_error),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(
                    f"PoleFigureResidualReport.{name} must be non-negative and finite."
                )
        object.__setattr__(self, "predicted_intensities", predicted)
        object.__setattr__(self, "residuals", residuals)

    def difference_figure(self) -> PoleFigureDifference:
        """The residual as a figure that can be projected and plotted.

        Purpose
        -------
        Turns the goodness-of-fit numbers into a diagnosis. ``residual_norm``
        says how badly the ODF misses; only the residual *figure* says
        **where**, and a systematic miss concentrated in one region of the
        specimen sphere means something quite different from noise spread over
        all of it — an unmodelled component, an uncorrected defocusing loss at
        high tilt, or a bandwidth too low to carry a sharp peak.

        Returns
        -------
        PoleFigureDifference
            Recalculated minus measured, on the measured directions. Positive
            means the ODF over-predicts. Read its ``describe()``, or hand it to
            ``pytex.plotting.plot_pole_figure_difference``.
        """

        return PoleFigureDifference(
            pole=self.pole_figure.pole,
            sample_directions=self.pole_figure.sample_directions,
            values=self.residuals,
            specimen_frame=self.pole_figure.specimen_frame,
            antipodal=self.pole_figure.antipodal,
            left_label="recalculated",
            right_label="measured",
            includes_symmetry_family=self.pole_figure.includes_symmetry_family,
            provenance=self.provenance,
        )

    def describe(self) -> str:
        """Prose summary: how well the ODF reproduces this figure, and where not."""

        difference = self.difference_figure()
        quality = (
            "The reconstruction reproduces this figure well"
            if self.relative_residual_norm <= 0.05
            else "The reconstruction reproduces this figure only approximately"
            if self.relative_residual_norm <= 0.2
            else "The reconstruction does not reproduce this figure"
        )
        return (
            f"Pole-figure residual over {self.observation_count} measured directions: "
            f"relative residual norm {self.relative_residual_norm:.4f}, mean absolute error "
            f"{self.mean_absolute_error:.4f}, maximum absolute error "
            f"{self.max_absolute_error:.4f}. {quality}. "
            f"{difference.describe()} An inversion that cannot reproduce its own input data is "
            "not usable whatever its internal convergence reported, so this comparison — not "
            "the solver's iteration count — is the acceptance test."
        )

    @classmethod
    def from_odf(
        cls,
        odf: ODF | HarmonicODF,
        pole_figure: PoleFigure,
        *,
        include_symmetry_family: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigureResidualReport:
        """Compare an ODF's predicted pole densities against a measured figure.

        Purpose
        -------
        The goodness-of-fit check on a PF-to-ODF inversion: reconstruct the pole
        density the ODF implies at the measured directions and report the
        residuals. An inversion that cannot reproduce its own input data is not
        usable, whatever its internal convergence said.

        Parameters
        ----------
        odf : ODF or HarmonicODF
            Either ODF representation is accepted.
        pole_figure : PoleFigure
            The measured figure to compare against.
        include_symmetry_family : bool
            Include the whole ``{hkl}`` family (default), matching the
            measurement.
        provenance : ProvenanceRecord, optional
        """

        predicted = np.asarray(
            odf.evaluate_pole_density(
                pole_figure.pole,
                pole_figure.sample_directions,
                include_symmetry_family=include_symmetry_family,
            ),
            dtype=np.float64,
        )
        if isinstance(odf, ODF):
            # A discrete ODF's pole density is a kernel-weighted *response*, not
            # a value in multiples of random: the kernel peaks at 1 rather than
            # integrating to 1, so a random texture returns the kernel's
            # spherical mean (about 0.016 at a 12 degree halfwidth). Differencing
            # that against a measured figure in m.r.d. would report a ~100%
            # residual for a perfect fit — a scale error, not a misfit. A
            # HarmonicODF needs no correction; its densities are already m.r.d.
            predicted = predicted / random_pole_density(odf.kernel)
        residuals = np.ascontiguousarray(predicted - pole_figure.intensities, dtype=np.float64)
        residuals.setflags(write=False)
        residual_norm = float(np.linalg.norm(residuals))
        observation_norm = max(float(np.linalg.norm(pole_figure.intensities)), 1e-12)
        return cls(
            pole_figure=pole_figure,
            predicted_intensities=predicted,
            residuals=residuals,
            residual_norm=residual_norm,
            relative_residual_norm=float(residual_norm / observation_norm),
            mean_absolute_error=float(np.mean(np.abs(residuals))),
            max_absolute_error=float(np.max(np.abs(residuals))),
            observation_count=int(residuals.shape[0]),
            provenance=provenance or pole_figure.provenance,
        )


@dataclass(frozen=True, slots=True)
class ODFReconstructionConfig:
    """A declared, reproducible PF-to-ODF reconstruction procedure.

    Purpose
    -------
    Bundles the correction spec, the algorithm choice, the kernel, and the
    regularization into one object, so a study's inversion settings live in a
    single declared configuration rather than scattered across call sites —
    and so the correction step can never be accidentally skipped.

    Attributes
    ----------
    algorithm : str
        ``"discrete"`` (weighted support; needs an orientation dictionary) or
        the harmonic series method.
    correction : PoleFigureCorrectionSpec, optional
        Applied before inversion; ``None`` means the figures are used as
        given.
    kernel : KernelSpec
        Smoothing kernel for the response model.
    regularization : float
        Tikhonov weight; larger is smoother, more stable, less detailed.
    Remaining attributes carry algorithm-specific settings.
    """

    algorithm: ReconstructionAlgorithm = "harmonic"
    kernel: KernelSpec = field(default_factory=KernelSpec)
    correction: PoleFigureCorrectionSpec | None = None
    regularization: float = 1e-6
    include_symmetry_family: bool = True
    degree_bandlimit: int = 6
    even_degrees_only: bool | None = None
    max_iterations: int = 500
    tolerance: float = 1e-8
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.algorithm not in {"discrete", "harmonic"}:
            raise ValueError("algorithm must be either 'discrete' or 'harmonic'.")
        if self.regularization < 0.0:
            raise ValueError("regularization must be non-negative.")
        if self.degree_bandlimit < 0:
            raise ValueError("degree_bandlimit must be non-negative.")
        if self.max_iterations <= 0:
            raise ValueError("max_iterations must be positive.")
        if self.tolerance <= 0.0:
            raise ValueError("tolerance must be positive.")

    def corrected_pole_figures(self, pole_figures: Sequence[PoleFigure]) -> tuple[PoleFigure, ...]:
        """Apply this configuration's correction spec to a set of pole figures.

        Returns the input unchanged when no correction is configured, so the
        call is always safe to make.
        """

        if self.correction is None:
            return tuple(pole_figures)
        return tuple(self.correction.apply(pole_figure) for pole_figure in pole_figures)

    def reconstruct(
        self,
        pole_figures: Sequence[PoleFigure],
        *,
        orientation_dictionary: object | None = None,
    ) -> ODFInversionReport | HarmonicODFReconstructionReport:
        """Run the configured PF-to-ODF reconstruction.

        Purpose
        -------
        One entry point that applies the correction spec and then dispatches to
        the chosen inversion algorithm, so the correction can never be
        accidentally skipped between the two steps.

        Parameters
        ----------
        pole_figures : sequence of PoleFigure
            The measured figures.
        orientation_dictionary : OrientationSet, optional
            Required by the ``"discrete"`` algorithm, which needs a support to
            place weights on; unused by the harmonic algorithm.

        Returns
        -------
        ODFInversionReport or HarmonicODFReconstructionReport
            Depending on the configured algorithm.
        """

        corrected = self.corrected_pole_figures(pole_figures)
        if self.algorithm == "discrete":
            if orientation_dictionary is None:
                raise ValueError("Discrete ODF reconstruction requires orientation_dictionary.")
            return ODF.invert_pole_figures(
                corrected,
                orientation_dictionary=orientation_dictionary,  # type: ignore[arg-type]
                kernel=self.kernel,
                regularization=self.regularization,
                include_symmetry_family=self.include_symmetry_family,
                max_iterations=self.max_iterations,
                tolerance=self.tolerance,
                provenance=self.provenance,
            )
        return HarmonicODF.invert_pole_figures(
            corrected,
            degree_bandlimit=self.degree_bandlimit,
            regularization=self.regularization,
            include_symmetry_family=self.include_symmetry_family,
            even_degrees_only=self.even_degrees_only,
            pole_kernel=self.kernel,
            provenance=self.provenance,
        )


def residual_reports_for_pole_figures(
    odf: ODF | HarmonicODF,
    pole_figures: Sequence[PoleFigure],
    *,
    include_symmetry_family: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> tuple[PoleFigureResidualReport, ...]:
    """Residual reports for an ODF against several measured pole figures.

    The batch form of :meth:`PoleFigureResidualReport.from_odf`, and the
    standard way to check a reconstruction against every figure that went
    into it.
    """

    return tuple(
        PoleFigureResidualReport.from_odf(
            odf,
            pole_figure,
            include_symmetry_family=include_symmetry_family,
            provenance=provenance,
        )
        for pole_figure in pole_figures
    )


__all__ = [
    "POLE_FIGURE_DEFOCUS_CALIBRATION_SCHEMA",
    "ODFReconstructionConfig",
    "PoleFigureCorrectionSpec",
    "PoleFigureDefocusCalibration",
    "PoleFigureResidualReport",
    "defocus_from_random_standard",
    "residual_reports_for_pole_figures",
]
