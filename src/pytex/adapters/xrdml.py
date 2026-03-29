from __future__ import annotations

import bz2
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from xml.etree import ElementTree

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.texture.models import KernelSpec, ODF, ODFInversionReport, PoleFigure

XRDML_NAMESPACE = "http://www.xrdml.com/XRDMeasurement/1.3"


def _open_xrdml_text(path: str | Path) -> str:
    source = Path(path)
    if source.suffix.lower() == ".bz2":
        return bz2.open(source, mode="rt", encoding="utf-8").read()
    return source.read_text(encoding="utf-8")


def _namespace(root: ElementTree.Element) -> str:
    if root.tag.startswith("{") and "}" in root.tag:
        return root.tag[: root.tag.index("}") + 1]
    return ""


def _find_text(element: ElementTree.Element, query: str) -> str | None:
    child = element.find(query)
    if child is None or child.text is None:
        return None
    value = child.text.strip()
    return value or None


def _scan_position_array(
    position_element: ElementTree.Element,
    *,
    point_count: int,
    namespace: str,
) -> np.ndarray:
    list_positions = position_element.find(namespace + "listPositions")
    if list_positions is not None and list_positions.text:
        values = np.fromstring(list_positions.text, sep=" ", dtype=np.float64)
    else:
        start_position = position_element.find(namespace + "startPosition")
        end_position = position_element.find(namespace + "endPosition")
        common_position = position_element.find(namespace + "commonPosition")
        if start_position is not None and end_position is not None:
            values = np.linspace(
                float(start_position.text),
                float(end_position.text),
                point_count,
                dtype=np.float64,
            )
        elif common_position is not None and common_position.text is not None:
            values = np.full(point_count, float(common_position.text), dtype=np.float64)
        else:
            raise ValueError("Unsupported XRDML positions encoding.")
    if values.shape != (point_count,):
        raise ValueError("XRDML positions arrays must align with the scan point count.")
    values = np.ascontiguousarray(values, dtype=np.float64)
    values.setflags(write=False)
    return values


@dataclass(frozen=True, slots=True)
class XRDMLPoleFigureMeasurement:
    phi_deg: np.ndarray
    psi_deg: np.ndarray
    intensity_grid: np.ndarray
    two_theta_deg: np.ndarray | None = None
    omega_deg: np.ndarray | None = None
    wavelength_angstrom: float | None = None
    source_path: str | None = None
    sample_name: str | None = None
    sample_mode: str | None = None
    measurement_axis: str | None = None
    scan_axis: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    comment_entries: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        phi = as_float_array(self.phi_deg, shape=(None, None))
        psi = as_float_array(self.psi_deg, shape=phi.shape)
        intensities = as_float_array(self.intensity_grid, shape=phi.shape)
        if np.any(~np.isfinite(intensities)) or np.any(intensities < 0.0):
            raise ValueError("XRDML intensities must be finite and non-negative.")
        object.__setattr__(self, "phi_deg", phi)
        object.__setattr__(self, "psi_deg", psi)
        object.__setattr__(self, "intensity_grid", intensities)
        if self.two_theta_deg is not None:
            object.__setattr__(
                self,
                "two_theta_deg",
                as_float_array(self.two_theta_deg, shape=phi.shape),
            )
        if self.omega_deg is not None:
            object.__setattr__(self, "omega_deg", as_float_array(self.omega_deg, shape=phi.shape))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
        object.__setattr__(self, "comment_entries", tuple(self.comment_entries))

    @property
    def sample_directions(self) -> np.ndarray:
        directions = spherical_angles_to_directions(self.psi_deg, self.phi_deg).reshape(-1, 3)
        directions = np.ascontiguousarray(directions, dtype=np.float64)
        directions.setflags(write=False)
        return directions

    @property
    def flattened_intensities(self) -> np.ndarray:
        flattened = np.ascontiguousarray(self.intensity_grid.reshape(-1), dtype=np.float64)
        flattened.setflags(write=False)
        return flattened

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(int(value) for value in self.intensity_grid.shape)

    def to_pole_figure(
        self,
        pole: CrystalPlane,
        *,
        specimen_frame: ReferenceFrame,
        antipodal: bool = True,
        sample_symmetry: SymmetrySpec | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigure:
        record = provenance or ProvenanceRecord(
            source_system="xrdml",
            source_path=self.source_path,
            metadata={
                "measurement_axis": self.measurement_axis or "",
                "scan_axis": self.scan_axis or "",
                "sample_mode": self.sample_mode or "",
                "reader": "pytex.adapters.read_xrdml_pole_figure",
            },
            notes=self.comment_entries,
        )
        return PoleFigure(
            pole=pole,
            sample_directions=self.sample_directions,
            intensities=self.flattened_intensities,
            specimen_frame=specimen_frame,
            antipodal=antipodal,
            sample_symmetry=sample_symmetry,
            provenance=record,
        )


def read_xrdml_pole_figure(path: str | Path) -> XRDMLPoleFigureMeasurement:
    source = Path(path)
    root = ElementTree.fromstring(_open_xrdml_text(source))
    namespace = _namespace(root)
    measurement = root.find(namespace + "xrdMeasurement")
    if measurement is None:
        raise ValueError("XRDML file does not contain an xrdMeasurement element.")
    scans = measurement.findall(namespace + "scan")
    if not scans:
        raise ValueError("XRDML pole-figure import requires at least one scan.")

    phi_rows: list[np.ndarray] = []
    psi_rows: list[np.ndarray] = []
    intensity_rows: list[np.ndarray] = []
    two_theta_rows: list[np.ndarray] = []
    omega_rows: list[np.ndarray] = []
    comments = tuple(
        entry.text.strip()
        for entry in root.findall(f".//{namespace}comment/{namespace}entry")
        if entry.text and entry.text.strip()
    )

    for scan in scans:
        points = scan.find(namespace + "dataPoints")
        if points is None:
            raise ValueError("XRDML scan is missing a dataPoints block.")
        intensities_node = points.find(namespace + "intensities")
        count_time_text = _find_text(points, namespace + "commonCountingTime") or "1.0"
        count_time = max(float(count_time_text), 1e-12)
        if intensities_node is not None and intensities_node.text:
            intensities = np.fromstring(intensities_node.text, sep=" ", dtype=np.float64)
        else:
            counts_node = points.find(namespace + "counts")
            if counts_node is None or counts_node.text is None:
                raise ValueError("XRDML scan is missing intensities or counts data.")
            intensities = np.fromstring(counts_node.text, sep=" ", dtype=np.float64) / count_time
        point_count = int(intensities.size)
        positions = {
            position.attrib["axis"]: _scan_position_array(
                position,
                point_count=point_count,
                namespace=namespace,
            )
            for position in points.findall(namespace + "positions")
            if "axis" in position.attrib
        }
        if "Phi" not in positions or "Psi" not in positions:
            raise ValueError("XRDML pole-figure import requires Phi and Psi positions.")
        phi_rows.append(positions["Phi"])
        psi_rows.append(positions["Psi"])
        intensity_rows.append(np.ascontiguousarray(intensities, dtype=np.float64))
        two_theta_rows.append(
            positions.get("2Theta", np.full(point_count, np.nan, dtype=np.float64))
        )
        omega_rows.append(positions.get("Omega", np.full(point_count, np.nan, dtype=np.float64)))

    used_wavelength = measurement.find(namespace + "usedWavelength")
    wavelength_text = None
    if used_wavelength is not None:
        wavelength_text = _find_text(used_wavelength, namespace + "kAlpha1")
    sample = root.find(namespace + "sample")
    sample_name = _find_text(sample, namespace + "name") if sample is not None else None
    return XRDMLPoleFigureMeasurement(
        phi_deg=np.vstack(phi_rows),
        psi_deg=np.vstack(psi_rows),
        intensity_grid=np.vstack(intensity_rows),
        two_theta_deg=np.vstack(two_theta_rows),
        omega_deg=np.vstack(omega_rows),
        wavelength_angstrom=None if wavelength_text is None else float(wavelength_text),
        source_path=str(source),
        sample_name=sample_name,
        sample_mode=measurement.attrib.get("sampleMode"),
        measurement_axis=measurement.attrib.get("measurementStepAxis"),
        scan_axis=scans[0].attrib.get("scanAxis"),
        metadata={"namespace": namespace.strip("{}")},
        comment_entries=comments,
    )


def load_xrdml_pole_figure(
    path: str | Path,
    *,
    pole: CrystalPlane,
    specimen_frame: ReferenceFrame,
    antipodal: bool = True,
    sample_symmetry: SymmetrySpec | None = None,
) -> PoleFigure:
    measurement = read_xrdml_pole_figure(path)
    return measurement.to_pole_figure(
        pole,
        specimen_frame=specimen_frame,
        antipodal=antipodal,
        sample_symmetry=sample_symmetry,
    )


def invert_xrdml_pole_figures(
    measurements: Sequence[XRDMLPoleFigureMeasurement | str | Path],
    *,
    poles: Sequence[CrystalPlane],
    specimen_frame: ReferenceFrame,
    orientation_dictionary,
    kernel: KernelSpec | None = None,
    regularization: float = 1e-6,
    include_symmetry_family: bool = True,
    antipodal: bool = True,
    sample_symmetry: SymmetrySpec | None = None,
    max_iterations: int = 500,
    tolerance: float = 1e-8,
    provenance: ProvenanceRecord | None = None,
) -> ODFInversionReport:
    if len(measurements) != len(poles):
        raise ValueError("measurements and poles must have matching lengths.")
    pole_figures = []
    for measurement, pole in zip(measurements, poles, strict=True):
        parsed = (
            measurement
            if isinstance(measurement, XRDMLPoleFigureMeasurement)
            else read_xrdml_pole_figure(measurement)
        )
        pole_figures.append(
            parsed.to_pole_figure(
                pole,
                specimen_frame=specimen_frame,
                antipodal=antipodal,
                sample_symmetry=sample_symmetry,
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
    "XRDMLPoleFigureMeasurement",
    "XRDML_NAMESPACE",
    "invert_xrdml_pole_figures",
    "load_xrdml_pole_figure",
    "read_xrdml_pole_figure",
]
