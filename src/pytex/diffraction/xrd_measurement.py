"""Measured powder-XRD profiles and convention-explicit profile comparison."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.xrd import PowderPattern, RadiationSpec

MEASURED_POWDER_PATTERN_SCHEMA = "pytex.diffraction.measured_powder_pattern"
POWDER_PATTERN_COMPARISON_SCHEMA = "pytex.diffraction.powder_pattern_comparison"
INTENSITY_UNITS = ("counts", "counts_per_second", "arbitrary")

IntensityUnit = Literal["counts", "counts_per_second", "arbitrary"]


@dataclass(frozen=True, slots=True)
class MeasuredPowderPattern:
    """A measured one-dimensional powder diffractogram.

    Use this type at the boundary between an instrument/export file and PyTex.
    It preserves the measured ``2*theta`` support, intensity meaning, optional
    standard uncertainties, radiation, provenance, and the explicit distinction
    between experimental and synthetic data. It does not imply background
    subtraction, normalization, or a structural model.
    """

    name: str
    two_theta_deg: np.ndarray
    intensity: np.ndarray
    standard_uncertainty: np.ndarray | None = None
    intensity_unit: IntensityUnit = "counts"
    radiation: RadiationSpec | None = None
    synthetic: bool = False
    metadata: Mapping[str, str] = field(default_factory=dict)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        axis = as_float_array(self.two_theta_deg, shape=(None,))
        intensity = as_float_array(self.intensity, shape=(None,))
        if not self.name.strip():
            raise ValueError("MeasuredPowderPattern.name must be non-empty.")
        if axis.size < 2:
            raise ValueError("A measured powder pattern needs at least two data points.")
        if axis.shape != intensity.shape:
            raise ValueError("Measured powder-pattern axis and intensity arrays must align.")
        if np.any(~np.isfinite(axis)) or np.any(np.diff(axis) <= 0.0):
            raise ValueError("Measured 2*theta values must be finite and strictly increasing.")
        if np.any(~np.isfinite(intensity)) or np.any(intensity < 0.0):
            raise ValueError("Measured intensities must be finite and non-negative.")
        if not np.any(intensity > 0.0):
            raise ValueError("A measured powder pattern cannot contain only zero intensity.")
        if self.intensity_unit not in INTENSITY_UNITS:
            raise ValueError(
                "MeasuredPowderPattern.intensity_unit must be one of "
                + ", ".join(INTENSITY_UNITS)
                + "."
            )
        uncertainty = self.standard_uncertainty
        if uncertainty is not None:
            uncertainty = as_float_array(uncertainty, shape=(None,))
            if uncertainty.shape != axis.shape:
                raise ValueError("Standard uncertainties must align with measured intensities.")
            if np.any(~np.isfinite(uncertainty)) or np.any(uncertainty <= 0.0):
                raise ValueError("Standard uncertainties must be finite and strictly positive.")
        object.__setattr__(self, "two_theta_deg", axis)
        object.__setattr__(self, "intensity", intensity)
        object.__setattr__(self, "standard_uncertainty", uncertainty)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def __len__(self) -> int:
        return int(self.two_theta_deg.size)

    def describe(self) -> str:
        """Return convention- and provenance-explicit scientific prose."""

        origin = "synthetic validation profile" if self.synthetic else "experimental profile"
        uncertainty = (
            " Per-point standard uncertainties are available for inverse-variance weighting."
            if self.standard_uncertainty is not None
            else " No per-point uncertainties are present; comparisons use unit weights."
        )
        radiation = (
            f" Radiation is {self.radiation.name} at "
            f"{self.radiation.wavelength_angstrom:.6f} angstrom."
            if self.radiation is not None
            else " Radiation was not declared."
        )
        source = (
            f" Source: {self.provenance.source_system}"
            + (
                f" ({self.provenance.source_path})."
                if self.provenance.source_path is not None
                else "."
            )
            if self.provenance is not None
            else " Source provenance was not supplied."
        )
        return (
            f"Measured powder-XRD pattern '{self.name}' is an explicitly labelled {origin} with "
            f"{len(self)} points from {self.two_theta_deg[0]:.4f} to "
            f"{self.two_theta_deg[-1]:.4f} degrees 2*theta; intensity is in "
            f"{self.intensity_unit}.{radiation}{uncertainty}{source}"
        )


@dataclass(frozen=True, slots=True)
class PowderPatternComparison:
    """An explainable measured-versus-simulated whole-profile comparison."""

    measured: MeasuredPowderPattern
    simulated: PowderPattern
    two_theta_deg: np.ndarray
    observed_intensity: np.ndarray
    calculated_intensity: np.ndarray
    residual_intensity: np.ndarray
    scale_factor: float
    background_offset: float
    profile_r_factor: float
    weighted_profile_r_factor: float
    correlation_coefficient: float
    weight_model: Literal["unit", "inverse_variance"]
    fitted_background: bool

    def __post_init__(self) -> None:
        arrays = tuple(
            as_float_array(value, shape=(None,))
            for value in (
                self.two_theta_deg,
                self.observed_intensity,
                self.calculated_intensity,
                self.residual_intensity,
            )
        )
        if arrays[0].size < 2 or len({array.shape for array in arrays}) != 1:
            raise ValueError("Powder-pattern comparison arrays must align and contain two points.")
        if np.any(~np.isfinite(np.concatenate(arrays))):
            raise ValueError("Powder-pattern comparison arrays must be finite.")
        if np.any(np.diff(arrays[0]) <= 0.0):
            raise ValueError("Powder-pattern comparison angles must be strictly increasing.")
        if np.any(arrays[1] < 0.0):
            raise ValueError("Observed powder-pattern intensities must be non-negative.")
        if not np.allclose(arrays[3], arrays[1] - arrays[2], rtol=1e-12, atol=1e-12):
            raise ValueError("Residual intensity must equal observed minus calculated intensity.")
        for name in (
            "scale_factor",
            "background_offset",
            "profile_r_factor",
            "weighted_profile_r_factor",
            "correlation_coefficient",
        ):
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f"PowderPatternComparison.{name} must be finite.")
        if self.scale_factor < 0.0:
            raise ValueError("PowderPatternComparison.scale_factor must be non-negative.")
        if self.profile_r_factor < 0.0 or self.weighted_profile_r_factor < 0.0:
            raise ValueError("Powder-pattern profile R factors must be non-negative.")
        if not -1.0 <= self.correlation_coefficient <= 1.0:
            raise ValueError("Correlation coefficient must lie in [-1, 1].")
        if self.weight_model not in {"unit", "inverse_variance"}:
            raise ValueError("Powder-pattern comparison weight model is not recognized.")
        object.__setattr__(self, "two_theta_deg", arrays[0])
        object.__setattr__(self, "observed_intensity", arrays[1])
        object.__setattr__(self, "calculated_intensity", arrays[2])
        object.__setattr__(self, "residual_intensity", arrays[3])

    @property
    def point_count(self) -> int:
        """Return the number of measured angles retained in the shared interval."""

        return int(self.two_theta_deg.size)

    def describe(self) -> str:
        """Return fit choices, metrics, limits, and their normative definition."""

        background = (
            f"and a constant background of {self.background_offset:.6g}"
            if self.fitted_background
            else "with no fitted background"
        )
        return (
            f"Compared {self.point_count} overlapping profile points from "
            f"{self.two_theta_deg[0]:.4f} to {self.two_theta_deg[-1]:.4f} degrees 2*theta. "
            f"Weighted least squares used {self.weight_model.replace('_', ' ')} weights, a "
            f"non-negative scale of {self.scale_factor:.6g}, {background}. "
            f"R_p = {self.profile_r_factor:.6g}, R_wp = "
            f"{self.weighted_profile_r_factor:.6g}, and Pearson correlation = "
            f"{self.correlation_coefficient:.6g}. R_p and R_wp follow the IUCr pdCIF "
            "definitions (_pd_proc_ls.prof_R_factor and _pd_proc_ls.prof_wR_factor). "
            "The calculation interpolates the simulated profile onto measured angles but does "
            "not shift peaks, refine structure, or model specimen/instrument physics."
        )


def _header_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("#"):
            break
        key, separator, value = stripped[1:].partition(":")
        if separator:
            metadata[key.strip().lower()] = value.strip()
    return metadata


def read_powder_xy(
    path: str | Path,
    *,
    delimiter: str | None = None,
    skip_rows: int = 0,
    two_theta_column: int = 0,
    intensity_column: int = 1,
    uncertainty_column: int | None = None,
    name: str | None = None,
    intensity_unit: IntensityUnit | None = None,
    radiation: RadiationSpec | None = None,
    synthetic: bool | None = None,
) -> MeasuredPowderPattern:
    """Read a commentable whitespace/CSV powder profile into the canonical model.

    Use this for conventional two- or three-column instrument exports before
    comparison with :class:`~pytex.diffraction.xrd.PowderPattern`. Numeric rows
    hold ``2*theta`` and intensity, with an optional standard-uncertainty
    column. Comment lines may declare ``name``, ``intensity_unit``,
    ``wavelength_angstrom``, ``radiation_name``, and ``synthetic`` as
    ``# key: value`` metadata. Explicit keyword arguments take precedence.

    Parameters
    ----------
    path
        Whitespace numeric file, or comma-separated file when its suffix is
        ``.csv`` (an explicit ``delimiter`` overrides detection).
    skip_rows, two_theta_column, intensity_column, uncertainty_column
        Layout controls for vendor exports. Column numbers are zero-based.
    name, intensity_unit, radiation, synthetic
        Semantic overrides for missing or incorrect header metadata.

    Returns
    -------
    MeasuredPowderPattern
        Validated, immutable measured arrays with file provenance.
    """

    source = Path(path)
    metadata = _header_metadata(source)
    effective_delimiter = (
        "," if delimiter is None and source.suffix.lower() == ".csv" else delimiter
    )
    values = np.loadtxt(
        source,
        comments="#",
        delimiter=effective_delimiter,
        skiprows=skip_rows,
        ndmin=2,
    )
    columns = [two_theta_column, intensity_column]
    if uncertainty_column is not None:
        columns.append(uncertainty_column)
    if min(columns) < 0 or max(columns) >= values.shape[1]:
        raise ValueError(
            f"Requested powder-profile column outside the {values.shape[1]} available columns."
        )
    if radiation is None and "wavelength_angstrom" in metadata:
        radiation = RadiationSpec(
            name=metadata.get("radiation_name", "declared radiation"),
            wavelength_angstrom=float(metadata["wavelength_angstrom"]),
        )
    declared_synthetic = metadata.get("synthetic", "false").lower() in {"1", "true", "yes"}
    is_synthetic = declared_synthetic if synthetic is None else synthetic
    unit = intensity_unit or metadata.get("intensity_unit", "counts")
    if unit not in INTENSITY_UNITS:
        raise ValueError(f"Unsupported measured intensity unit '{unit}'.")
    provenance = ProvenanceRecord(
        source_system="powder_xy",
        source_path=str(source),
        metadata={**metadata, "synthetic": str(is_synthetic).lower()},
        notes=("Synthetic validation data; not an experimental measurement.",)
        if is_synthetic
        else (),
    )
    return MeasuredPowderPattern(
        name=name or metadata.get("name", source.stem),
        two_theta_deg=values[:, two_theta_column],
        intensity=values[:, intensity_column],
        standard_uncertainty=(
            None if uncertainty_column is None else values[:, uncertainty_column]
        ),
        intensity_unit=unit,  # type: ignore[arg-type]
        radiation=radiation,
        synthetic=is_synthetic,
        metadata=metadata,
        provenance=provenance,
    )


def write_powder_xy(pattern: MeasuredPowderPattern, path: str | Path) -> Path:
    """Write the canonical measured profile as deterministic whitespace columns.

    Use this at a workflow boundary when another tool needs a plain numeric
    profile rather than the richer JSON contract. The returned path contains a
    comment metadata header followed by ``2*theta``, intensity, and optional
    standard-uncertainty columns; no normalization or rounding beyond the
    declared 12-significant-digit text representation is applied.
    """

    output = Path(path)
    header = [
        "pytex-powder-pattern: 1",
        f"name: {pattern.name}",
        f"synthetic: {str(pattern.synthetic).lower()}",
        f"intensity_unit: {pattern.intensity_unit}",
    ]
    if pattern.radiation is not None:
        header.extend(
            (
                f"radiation_name: {pattern.radiation.name}",
                f"wavelength_angstrom: {pattern.radiation.wavelength_angstrom:.12g}",
            )
        )
    columns = [pattern.two_theta_deg, pattern.intensity]
    labels = "two_theta_deg intensity"
    if pattern.standard_uncertainty is not None:
        columns.append(pattern.standard_uncertainty)
        labels += " standard_uncertainty"
    header.append(f"columns: {labels}")
    np.savetxt(output, np.column_stack(columns), fmt="%.12g", header="\n".join(header))
    return output


def compare_powder_patterns(
    measured: MeasuredPowderPattern,
    simulated: PowderPattern,
    *,
    fit_background: bool = True,
) -> PowderPatternComparison:
    """Scale and compare a simulated profile on the measured angular support.

    Use this for phase-identification diagnostics and simulation validation,
    not for structural refinement. The simulated profile is linearly
    interpolated only inside the shared ``2*theta`` interval. Scale and optional
    constant background are fitted by weighted least squares. Per-point
    measured standard uncertainties imply inverse-variance weights; otherwise
    all points receive unit weight.

    Parameters
    ----------
    measured
        Observed profile and its uncertainty/provenance semantics.
    simulated
        Kinematic profile to interpolate onto the measured angular support.
    fit_background
        Fit one constant additive offset in addition to a non-negative scale.

    Returns
    -------
    PowderPatternComparison
        Observed, calculated and residual profiles plus scale/background,
        IUCr ``R_p``/``R_wp``, correlation, and an explainable summary.
    """

    lower = max(float(measured.two_theta_deg[0]), float(simulated.two_theta_grid_deg[0]))
    upper = min(float(measured.two_theta_deg[-1]), float(simulated.two_theta_grid_deg[-1]))
    mask = (measured.two_theta_deg >= lower) & (measured.two_theta_deg <= upper)
    if np.count_nonzero(mask) < 2:
        raise ValueError("Measured and simulated powder patterns need at least two overlap points.")
    axis = measured.two_theta_deg[mask]
    observed = measured.intensity[mask]
    if not np.any(observed > 0.0):
        raise ValueError("The measured profile has no positive intensity in the overlap interval.")
    simulated_on_measured = np.interp(
        axis, simulated.two_theta_grid_deg, simulated.intensity_grid
    )
    if not np.any(simulated_on_measured > 0.0):
        raise ValueError("The simulated profile has no positive intensity in the overlap interval.")
    if measured.standard_uncertainty is None:
        weights = np.ones_like(observed)
        weight_model: Literal["unit", "inverse_variance"] = "unit"
    else:
        weights = 1.0 / np.square(measured.standard_uncertainty[mask])
        weight_model = "inverse_variance"
    if fit_background:
        design = np.column_stack((simulated_on_measured, np.ones_like(simulated_on_measured)))
    else:
        design = simulated_on_measured[:, None]
    root_weight = np.sqrt(weights)
    coefficients, *_ = np.linalg.lstsq(
        design * root_weight[:, None], observed * root_weight, rcond=None
    )
    scale = float(coefficients[0])
    if scale < 0.0:
        scale = 0.0
        background = (
            float(np.average(observed, weights=weights)) if fit_background else 0.0
        )
    else:
        background = float(coefficients[1]) if fit_background else 0.0
    calculated = scale * simulated_on_measured + background
    residual = observed - calculated
    rp = float(np.sum(np.abs(residual)) / np.sum(observed))
    rwp = float(
        np.sqrt(np.sum(weights * np.square(residual)) / np.sum(weights * np.square(observed)))
    )
    if np.isclose(np.std(observed), 0.0) or np.isclose(np.std(calculated), 0.0):
        correlation = 0.0
    else:
        correlation = float(np.corrcoef(observed, calculated)[0, 1])
    return PowderPatternComparison(
        measured=measured,
        simulated=simulated,
        two_theta_deg=axis,
        observed_intensity=observed,
        calculated_intensity=calculated,
        residual_intensity=residual,
        scale_factor=scale,
        background_offset=background,
        profile_r_factor=rp,
        weighted_profile_r_factor=rwp,
        correlation_coefficient=correlation,
        weight_model=weight_model,
        fitted_background=fit_background,
    )


__all__ = [
    "INTENSITY_UNITS",
    "MEASURED_POWDER_PATTERN_SCHEMA",
    "POWDER_PATTERN_COMPARISON_SCHEMA",
    "MeasuredPowderPattern",
    "PowderPatternComparison",
    "compare_powder_patterns",
    "read_powder_xy",
    "write_powder_xy",
]
