from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core.provenance import ProvenanceRecord
from pytex.texture.harmonics import HarmonicODF, HarmonicODFReconstructionReport
from pytex.texture.models import ODF, KernelSpec, ODFInversionReport, PoleFigure

CorrectionPolicy = Literal["clip_zero", "raise"]
ReconstructionAlgorithm = Literal["discrete", "harmonic"]


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
        Defocus division first, then background subtraction, then scaling.

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

        intensities = np.asarray(pole_figure.intensities, dtype=np.float64)
        if self.defocus_factors is not None:
            if self.defocus_factors.shape != intensities.shape:
                raise ValueError("defocus_factors must match the pole-figure intensity shape.")
            intensities = intensities / self.defocus_factors
        corrected = self.scale * (intensities - self.background)
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

        predicted = odf.evaluate_pole_density(
            pole_figure.pole,
            pole_figure.sample_directions,
            include_symmetry_family=include_symmetry_family,
        )
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
    "ODFReconstructionConfig",
    "PoleFigureCorrectionSpec",
    "PoleFigureResidualReport",
    "residual_reports_for_pole_figures",
]
