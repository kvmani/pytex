from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, MillerIndex, Phase
from pytex.core.orientation import OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.sphere import raster_solid_angle_weights
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.texture.models import ODF, KernelSpec, ODFInversionReport, PoleFigure

_INTENSITY_NORMALIZATION_MODES = {"none", "max", "sum", "mrd"}
_FORMAT_EXTENSIONS = {".ppf": "PPF", ".epf": "EPF"}


def _normalize_intensity_grid(
    intensities: np.ndarray,
    *,
    mode: str,
    polar_deg: np.ndarray | None = None,
) -> np.ndarray:
    if mode not in _INTENSITY_NORMALIZATION_MODES:
        raise ValueError(
            "intensity_normalization must be one of 'none', 'max', 'sum', or 'mrd'."
        )
    normalized = np.array(intensities, copy=True, dtype=np.float64)
    if mode == "none":
        normalized = np.ascontiguousarray(normalized, dtype=np.float64)
        normalized.setflags(write=False)
        return normalized
    if mode == "mrd":
        if polar_deg is None:  # pragma: no cover - callers always supply the raster
            raise ValueError("The 'mrd' normalization requires the raster polar angles.")
        weights = raster_solid_angle_weights(np.asarray(polar_deg, dtype=np.float64).reshape(-1))
        scale = float(np.sum(weights * normalized.reshape(-1)))
    else:
        scale = float(np.nanmax(normalized)) if mode == "max" else float(np.sum(normalized))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(
            "Cannot normalize LaboTex intensities because the selected scale is non-positive."
        )
    normalized /= scale
    normalized = np.ascontiguousarray(normalized, dtype=np.float64)
    normalized.setflags(write=False)
    return normalized


def _split_numeric_tokens(lines: Sequence[str]) -> np.ndarray:
    numeric_tokens: list[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        numeric_tokens.extend(float(token) for token in stripped.split())
    values = np.asarray(numeric_tokens, dtype=np.float64)
    values = np.ascontiguousarray(values, dtype=np.float64)
    values.setflags(write=False)
    return values


@dataclass(frozen=True, slots=True)
class LaboTexPoleFigureDescriptor:
    """The scan geometry of one pole figure in a LaboTex file.

    Attributes
    ----------
    two_theta_deg : float
        Diffraction angle of the measured reflection.
    alpha_start_deg, alpha_end_deg, alpha_step_deg : float
        Tilt-axis scan range and step.
    beta_start_deg, beta_end_deg, beta_step_deg : float
        Azimuth-axis scan range and step.
    scan_index : int
    hkl : tuple of int
        The measured reflection.
    background_flag : int
        Whether a background correction was recorded by the source system.
    """

    two_theta_deg: float
    alpha_start_deg: float
    alpha_end_deg: float
    alpha_step_deg: float
    beta_start_deg: float
    beta_end_deg: float
    beta_step_deg: float
    scan_index: int
    hkl: tuple[int, int, int]
    background_flag: int
    figure_flag: int

    @property
    def alpha_values_deg(self) -> np.ndarray:
        """The pole-figure tilt (alpha) sample positions, in degrees.
        """

        count = round((self.alpha_end_deg - self.alpha_start_deg) / self.alpha_step_deg) + 1
        values = np.linspace(
            self.alpha_start_deg,
            self.alpha_end_deg,
            count,
            dtype=np.float64,
        )
        values = np.ascontiguousarray(values, dtype=np.float64)
        values.setflags(write=False)
        return values

    @property
    def beta_values_deg(self) -> np.ndarray:
        """The pole-figure azimuth (beta) sample positions, in degrees.
        """

        count = round((self.beta_end_deg - self.beta_start_deg) / self.beta_step_deg) + 1
        values = np.linspace(
            self.beta_start_deg,
            self.beta_end_deg,
            count,
            dtype=np.float64,
        )
        values = np.ascontiguousarray(values, dtype=np.float64)
        values.setflags(write=False)
        return values

    @property
    def shape(self) -> tuple[int, int]:
        """Grid shape of the measurement as ``(n_alpha, n_beta)``.
        """

        return (int(self.alpha_values_deg.size), int(self.beta_values_deg.size))

    @property
    def sample_directions(self) -> np.ndarray:
        """Specimen directions of every grid point, as unit vectors.

        Converts the LaboTex ``(alpha, beta)`` angle grid into the Cartesian
        specimen directions a PyTex pole figure is defined on.
        """

        alpha_grid = np.repeat(self.alpha_values_deg[:, None], self.beta_values_deg.size, axis=1)
        beta_grid = np.repeat(self.beta_values_deg[None, :], self.alpha_values_deg.size, axis=0)
        directions = spherical_angles_to_directions(alpha_grid, beta_grid).reshape(-1, 3)
        directions = np.ascontiguousarray(directions, dtype=np.float64)
        directions.setflags(write=False)
        return directions


@dataclass(frozen=True, slots=True)
class LaboTexPoleFigureMeasurement:
    """A parsed LaboTex pole-figure file: geometry, intensities, and metadata.

    Purpose
    -------
    The intermediate between the file and PyTex pole figures. It deliberately
    holds the measurement *as recorded* — raw intensity grids and scan
    descriptors — without crystallographic interpretation; attaching a phase
    and poles is a separate, explicit step via :meth:`to_pole_figures`.

    Attributes
    ----------
    title : str
    format_kind : str
        Which LaboTex format variant was parsed.
    lattice_parameters : tuple of float
        The six cell parameters as recorded in the file.
    descriptors : tuple of LaboTexPoleFigureDescriptor
        Scan geometry, one per pole figure.
    intensity_grids : tuple of np.ndarray
        Raw intensity grids, one per pole figure.
    source_path : str, optional
    metadata : Mapping[str, str]
    comments : tuple of str
        File comments, retained for traceability.
    """

    title: str
    format_kind: str
    lattice_parameters: tuple[float, float, float, float, float, float]
    descriptors: tuple[LaboTexPoleFigureDescriptor, ...]
    intensity_grids: tuple[np.ndarray, ...]
    source_path: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    comments: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.format_kind not in {"PPF", "EPF"}:
            raise ValueError("format_kind must be either 'PPF' or 'EPF'.")
        if len(self.descriptors) == 0:
            raise ValueError("LaboTexPoleFigureMeasurement requires at least one pole figure.")
        if len(self.descriptors) != len(self.intensity_grids):
            raise ValueError("descriptors and intensity_grids must have matching lengths.")
        normalized_grids: list[np.ndarray] = []
        for descriptor, grid in zip(self.descriptors, self.intensity_grids, strict=True):
            grid_array = as_float_array(grid, shape=descriptor.shape)
            if np.any(~np.isfinite(grid_array)) or np.any(grid_array < 0.0):
                raise ValueError("LaboTex intensity grids must be finite and non-negative.")
            normalized_grids.append(grid_array)
        object.__setattr__(self, "intensity_grids", tuple(normalized_grids))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
        object.__setattr__(self, "comments", tuple(self.comments))

    def normalized_intensity_grids(self, *, mode: str = "none") -> tuple[np.ndarray, ...]:
        """The measured intensity grids under a chosen normalization.

        Parameters
        ----------
        mode : str
            Normalization to apply. ``"none"`` (default) leaves the raw counts,
            which is the honest starting point; the other modes rescale so
            figures from different measurements become comparable. Because the
            choice changes every downstream density, it is explicit rather than
            implied.

        Returns
        -------
        tuple of np.ndarray
            One grid per pole figure in the file.
        """

        return tuple(
            _normalize_intensity_grid(
                grid,
                mode=mode,
                # LaboTex alpha is the tilt from the specimen normal, so it is
                # the polar angle whose bands carry the solid angle.
                polar_deg=np.repeat(
                    descriptor.alpha_values_deg[:, None],
                    descriptor.beta_values_deg.size,
                    axis=1,
                ),
            )
            for descriptor, grid in zip(self.descriptors, self.intensity_grids, strict=True)
        )

    def to_pole_figures(
        self,
        phase: Phase,
        *,
        specimen_frame: ReferenceFrame,
        antipodal: bool = True,
        sample_symmetry: SymmetrySpec | None = None,
        intensity_normalization: str = "none",
        provenance: ProvenanceRecord | None = None,
    ) -> tuple[PoleFigure, ...]:
        """Convert the measurement into PyTex pole figures.

        Attaches the phase, the poles, and the specimen frame to the measured
        grids, producing figures that can be inverted to an ODF or compared
        against a computed one.
        """

        normalized_grids = self.normalized_intensity_grids(mode=intensity_normalization)
        figures: list[PoleFigure] = []
        for descriptor, grid in zip(self.descriptors, normalized_grids, strict=True):
            pole = CrystalPlane(
                miller=MillerIndex(np.asarray(descriptor.hkl, dtype=np.int64), phase=phase),
                phase=phase,
            )
            record = provenance or ProvenanceRecord(
                source_system="labotex",
                source_path=self.source_path,
                metadata={
                    "format_kind": self.format_kind,
                    "title": self.title,
                    "intensity_normalization": intensity_normalization,
                    "two_theta_deg": f"{descriptor.two_theta_deg:.6f}",
                    "scan_index": str(descriptor.scan_index),
                },
                notes=self.comments,
            )
            figures.append(
                PoleFigure(
                    pole=pole,
                    sample_directions=descriptor.sample_directions,
                    intensities=np.ascontiguousarray(grid.reshape(-1), dtype=np.float64),
                    specimen_frame=specimen_frame,
                    antipodal=antipodal,
                    sample_symmetry=sample_symmetry,
                    # A LaboTex export is a measured density raster, not a
                    # cloud of individually weighted poles.
                    sampling="sampled_density",
                    provenance=record,
                )
            )
        return tuple(figures)


def read_labotex_pole_figures(path: str | Path) -> LaboTexPoleFigureMeasurement:
    """Parse a LaboTex pole-figure file into a measurement object.

    Reads the angle grid and intensity blocks without interpreting them
    crystallographically; use :meth:`LaboTexPoleFigureMeasurement.to_pole_figures`
    to attach phase and pole meaning.
    """

    source = Path(path)
    suffix = source.suffix.lower()
    if suffix not in _FORMAT_EXTENSIONS:
        raise ValueError("LaboTex reader supports .ppf and .epf files only.")
    try:
        lines = source.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        lines = source.read_text(encoding="latin-1").splitlines()
    nonempty_lines = [line.rstrip() for line in lines if line.strip()]
    if len(nonempty_lines) < 7:
        raise ValueError("LaboTex pole figure file is too short to be valid.")
    title = nonempty_lines[0].strip()
    structure_header_index = next(
        (index for index, line in enumerate(nonempty_lines) if "Structure Code" in line),
        None,
    )
    if structure_header_index is None or structure_header_index + 1 >= len(nonempty_lines):
        raise ValueError("LaboTex pole figure file is missing the structure header block.")
    lattice_tokens = nonempty_lines[structure_header_index + 1].split()
    if len(lattice_tokens) < 7:
        raise ValueError("LaboTex lattice line must contain structure code and six cell values.")
    lattice_parameters = tuple(float(token) for token in lattice_tokens[1:7])
    count_line = nonempty_lines[structure_header_index + 2]
    count_tokens = count_line.split()
    if not count_tokens:
        raise ValueError("LaboTex pole figure count line is missing.")
    pole_figure_count = int(count_tokens[0])
    descriptor_header_index = structure_header_index + 3
    descriptor_lines = nonempty_lines[
        descriptor_header_index + 1 : descriptor_header_index + 1 + pole_figure_count
    ]
    descriptors: list[LaboTexPoleFigureDescriptor] = []
    for line in descriptor_lines:
        tokens = line.split()
        if len(tokens) < 12:
            raise ValueError("LaboTex descriptor lines must contain at least 12 tokens.")
        background_flag = int(tokens[11]) if len(tokens) >= 13 else 0
        figure_flag = int(tokens[12]) if len(tokens) >= 13 else int(tokens[11])
        descriptors.append(
            LaboTexPoleFigureDescriptor(
                two_theta_deg=float(tokens[0]),
                alpha_start_deg=float(tokens[1]),
                alpha_end_deg=float(tokens[2]),
                alpha_step_deg=float(tokens[3]),
                beta_start_deg=float(tokens[4]),
                beta_end_deg=float(tokens[5]),
                beta_step_deg=float(tokens[6]),
                scan_index=int(tokens[7]),
                hkl=(int(tokens[8]), int(tokens[9]), int(tokens[10])),
                background_flag=background_flag,
                figure_flag=figure_flag,
            )
        )
    numeric_values = _split_numeric_tokens(
        nonempty_lines[descriptor_header_index + 1 + pole_figure_count :]
    )
    intensity_grids: list[np.ndarray] = []
    offset = 0
    for descriptor in descriptors:
        size = descriptor.shape[0] * descriptor.shape[1]
        if offset + size > numeric_values.size:
            raise ValueError("LaboTex pole figure file ended before all grids were populated.")
        grid = np.ascontiguousarray(
            numeric_values[offset : offset + size].reshape(descriptor.shape),
            dtype=np.float64,
        )
        grid.setflags(write=False)
        intensity_grids.append(grid)
        offset += size
    if offset != numeric_values.size:
        raise ValueError("LaboTex pole figure file contains trailing numeric data after the grids.")
    return LaboTexPoleFigureMeasurement(
        title=title,
        format_kind=_FORMAT_EXTENSIONS[suffix],
        lattice_parameters=(
            float(lattice_parameters[0]),
            float(lattice_parameters[1]),
            float(lattice_parameters[2]),
            float(lattice_parameters[3]),
            float(lattice_parameters[4]),
            float(lattice_parameters[5]),
        ),
        descriptors=tuple(descriptors),
        intensity_grids=tuple(intensity_grids),
        source_path=str(source),
        metadata={
            "format_kind": _FORMAT_EXTENSIONS[suffix],
            "structure_code": lattice_tokens[0],
            "pole_figure_count": str(pole_figure_count),
        },
        comments=tuple(nonempty_lines[:structure_header_index]),
    )


def load_labotex_pole_figures(
    path: str | Path,
    *,
    phase: Phase,
    specimen_frame: ReferenceFrame,
    antipodal: bool = True,
    sample_symmetry: SymmetrySpec | None = None,
    intensity_normalization: str = "none",
) -> tuple[PoleFigure, ...]:
    """Read a LaboTex file and return PyTex pole figures directly.

    Convenience composition of :func:`read_labotex_pole_figures` with
    :meth:`LaboTexPoleFigureMeasurement.to_pole_figures`.
    """

    measurement = read_labotex_pole_figures(path)
    return measurement.to_pole_figures(
        phase,
        specimen_frame=specimen_frame,
        antipodal=antipodal,
        sample_symmetry=sample_symmetry,
        intensity_normalization=intensity_normalization,
    )


def invert_labotex_pole_figures(
    measurements: Sequence[LaboTexPoleFigureMeasurement | str | Path],
    *,
    phase: Phase,
    specimen_frame: ReferenceFrame,
    orientation_dictionary: OrientationSet,
    kernel: KernelSpec | None = None,
    regularization: float = 1e-6,
    include_symmetry_family: bool = True,
    antipodal: bool = True,
    sample_symmetry: SymmetrySpec | None = None,
    intensity_normalization: str = "none",
    max_iterations: int = 500,
    tolerance: float = 1e-8,
    provenance: ProvenanceRecord | None = None,
) -> ODFInversionReport:
    """Read LaboTex pole figures and reconstruct an ODF from them.

    The end-to-end path from a measured LaboTex file to an orientation
    distribution. The inversion is ill-posed — see
    :meth:`~pytex.texture.ODF.invert_pole_figures` — so check the returned
    residuals before using the result.
    """

    pole_figures: list[PoleFigure] = []
    for measurement in measurements:
        parsed = (
            measurement
            if isinstance(measurement, LaboTexPoleFigureMeasurement)
            else read_labotex_pole_figures(measurement)
        )
        pole_figures.extend(
            parsed.to_pole_figures(
                phase,
                specimen_frame=specimen_frame,
                antipodal=antipodal,
                sample_symmetry=sample_symmetry,
                intensity_normalization=intensity_normalization,
            )
        )
    return ODF.invert_pole_figures(
        pole_figures,
        orientation_dictionary=orientation_dictionary,
        kernel=kernel,
        regularization=regularization,
        include_symmetry_family=include_symmetry_family,
        max_iterations=max_iterations,
        tolerance=tolerance,
        provenance=provenance,
    )


__all__ = [
    "LaboTexPoleFigureDescriptor",
    "LaboTexPoleFigureMeasurement",
    "invert_labotex_pole_figures",
    "load_labotex_pole_figures",
    "read_labotex_pole_figures",
]
