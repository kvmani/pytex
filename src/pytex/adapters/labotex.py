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
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.texture.models import ODF, KernelSpec, ODFInversionReport, PoleFigure

_INTENSITY_NORMALIZATION_MODES = {"none", "max", "sum"}
_FORMAT_EXTENSIONS = {".ppf": "PPF", ".epf": "EPF"}


def _normalize_intensity_grid(intensities: np.ndarray, *, mode: str) -> np.ndarray:
    if mode not in _INTENSITY_NORMALIZATION_MODES:
        raise ValueError("intensity_normalization must be one of 'none', 'max', or 'sum'.")
    normalized = np.array(intensities, copy=True, dtype=np.float64)
    if mode == "none":
        normalized = np.ascontiguousarray(normalized, dtype=np.float64)
        normalized.setflags(write=False)
        return normalized
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
        return (int(self.alpha_values_deg.size), int(self.beta_values_deg.size))

    @property
    def sample_directions(self) -> np.ndarray:
        alpha_grid = np.repeat(self.alpha_values_deg[:, None], self.beta_values_deg.size, axis=1)
        beta_grid = np.repeat(self.beta_values_deg[None, :], self.alpha_values_deg.size, axis=0)
        directions = spherical_angles_to_directions(alpha_grid, beta_grid).reshape(-1, 3)
        directions = np.ascontiguousarray(directions, dtype=np.float64)
        directions.setflags(write=False)
        return directions


@dataclass(frozen=True, slots=True)
class LaboTexPoleFigureMeasurement:
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
        return tuple(_normalize_intensity_grid(grid, mode=mode) for grid in self.intensity_grids)

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
                    provenance=record,
                )
            )
        return tuple(figures)


def read_labotex_pole_figures(path: str | Path) -> LaboTexPoleFigureMeasurement:
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
