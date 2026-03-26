from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

from pytex.core._arrays import normalize_vector
from pytex.core.conventions import FrameDomain
from pytex.core.frames import FrameTransform, ReferenceFrame
from pytex.core.provenance import ProvenanceRecord

_CALIBRATION_STATUSES = {"nominal", "calibrated", "refined", "verified"}
_ACQUISITION_MODALITIES = {"generic", "ebsd", "xrd", "neutron", "tem"}
_RADIATION_TYPES = {"electron", "xray", "neutron", "generic"}


def _freeze_float_mapping(values: Mapping[str, float]) -> Mapping[str, float]:
    frozen = MappingProxyType({str(key): float(value) for key, value in values.items()})
    if any(not np.isfinite(value) or value < 0.0 for value in frozen.values()):
        raise ValueError("Mapping values must be finite and non-negative.")
    return frozen


@dataclass(frozen=True, slots=True)
class MeasurementQuality:
    confidence: float | None = None
    valid_fraction: float | None = None
    masked_fraction: float | None = None
    uncertainty: Mapping[str, float] = field(default_factory=dict)
    flags: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("confidence", self.confidence),
            ("valid_fraction", self.valid_fraction),
            ("masked_fraction", self.masked_fraction),
        ):
            if value is not None and not 0.0 <= value <= 1.0:
                raise ValueError(f"MeasurementQuality.{name} must lie in [0, 1] when provided.")
        object.__setattr__(self, "uncertainty", _freeze_float_mapping(self.uncertainty))
        object.__setattr__(self, "flags", tuple(str(flag) for flag in self.flags))


@dataclass(frozen=True, slots=True)
class CalibrationRecord:
    source: str
    status: str = "nominal"
    residual_error: float | None = None
    parameter_uncertainty: Mapping[str, float] = field(default_factory=dict)
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        normalized_source = self.source.strip()
        if not normalized_source:
            raise ValueError("CalibrationRecord.source must be non-empty.")
        normalized_status = self.status.strip().lower()
        if normalized_status not in _CALIBRATION_STATUSES:
            supported = ", ".join(sorted(_CALIBRATION_STATUSES))
            raise ValueError(
                f"Unsupported calibration status '{self.status}'. Supported statuses: {supported}"
            )
        if self.residual_error is not None and (
            not np.isfinite(self.residual_error) or self.residual_error < 0.0
        ):
            raise ValueError("CalibrationRecord.residual_error must be finite and non-negative.")
        object.__setattr__(self, "source", normalized_source)
        object.__setattr__(self, "status", normalized_status)
        object.__setattr__(
            self, "parameter_uncertainty", _freeze_float_mapping(self.parameter_uncertainty)
        )
        object.__setattr__(self, "notes", tuple(str(note) for note in self.notes))

    @property
    def is_calibrated(self) -> bool:
        return self.status != "nominal"


@dataclass(frozen=True, slots=True)
class ScatteringSetup:
    laboratory_frame: ReferenceFrame
    radiation_type: str = "electron"
    incident_beam_direction: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 0.0, 1.0], dtype=np.float64)
    )
    wavelength_angstrom: float | None = None
    beam_energy_kev: float | None = None
    scattering_mode: str = "elastic"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.laboratory_frame.domain is not FrameDomain.LABORATORY:
            raise ValueError(
                "ScatteringSetup.laboratory_frame must belong to the laboratory domain."
            )
        radiation_type = self.radiation_type.strip().lower()
        if radiation_type not in _RADIATION_TYPES:
            supported = ", ".join(sorted(_RADIATION_TYPES))
            raise ValueError(
                f"Unsupported radiation_type '{self.radiation_type}'. Supported types: {supported}"
            )
        if self.wavelength_angstrom is None and self.beam_energy_kev is None:
            raise ValueError(
                "ScatteringSetup requires either wavelength_angstrom or beam_energy_kev."
            )
        if self.wavelength_angstrom is not None and self.wavelength_angstrom <= 0.0:
            raise ValueError("ScatteringSetup.wavelength_angstrom must be strictly positive.")
        if self.beam_energy_kev is not None and self.beam_energy_kev <= 0.0:
            raise ValueError("ScatteringSetup.beam_energy_kev must be strictly positive.")
        if not self.scattering_mode.strip():
            raise ValueError("ScatteringSetup.scattering_mode must be non-empty.")
        object.__setattr__(self, "radiation_type", radiation_type)
        object.__setattr__(
            self,
            "incident_beam_direction",
            normalize_vector(self.incident_beam_direction),
        )
        object.__setattr__(self, "scattering_mode", self.scattering_mode.strip().lower())

    @property
    def effective_wavelength_angstrom(self) -> float:
        if self.wavelength_angstrom is not None:
            return float(self.wavelength_angstrom)
        if self.beam_energy_kev is None:
            raise ValueError(
                "ScatteringSetup.beam_energy_kev must be provided when wavelength is absent."
            )
        voltage = float(self.beam_energy_kev) * 1000.0
        numerator = 12.2639
        denominator = np.sqrt(voltage * (1.0 + 0.97845e-6 * voltage))
        return float(numerator / denominator)


@dataclass(frozen=True, slots=True)
class AcquisitionGeometry:
    specimen_frame: ReferenceFrame
    modality: str = "generic"
    map_frame: ReferenceFrame | None = None
    detector_frame: ReferenceFrame | None = None
    laboratory_frame: ReferenceFrame | None = None
    reciprocal_frame: ReferenceFrame | None = None
    specimen_to_map: FrameTransform | None = None
    specimen_to_detector: FrameTransform | None = None
    specimen_to_laboratory: FrameTransform | None = None
    laboratory_to_reciprocal: FrameTransform | None = None
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError(
                "AcquisitionGeometry.specimen_frame must belong to the specimen domain."
            )
        modality = self.modality.strip().lower()
        if modality not in _ACQUISITION_MODALITIES:
            supported = ", ".join(sorted(_ACQUISITION_MODALITIES))
            raise ValueError(
                "Unsupported acquisition modality "
                f"'{self.modality}'. Supported modalities: {supported}"
            )
        if self.map_frame is not None and self.map_frame.domain is not FrameDomain.MAP:
            raise ValueError("AcquisitionGeometry.map_frame must belong to the map domain.")
        if (
            self.detector_frame is not None
            and self.detector_frame.domain is not FrameDomain.DETECTOR
        ):
            raise ValueError(
                "AcquisitionGeometry.detector_frame must belong to the detector domain."
            )
        if (
            self.laboratory_frame is not None
            and self.laboratory_frame.domain is not FrameDomain.LABORATORY
        ):
            raise ValueError(
                "AcquisitionGeometry.laboratory_frame must belong to the laboratory domain."
            )
        if (
            self.reciprocal_frame is not None
            and self.reciprocal_frame.domain is not FrameDomain.RECIPROCAL
        ):
            raise ValueError(
                "AcquisitionGeometry.reciprocal_frame must belong to the reciprocal domain."
            )
        self._validate_transform(
            self.specimen_to_map,
            target_frame=self.map_frame,
            expected_source=self.specimen_frame,
            label="specimen_to_map",
        )
        self._validate_transform(
            self.specimen_to_detector,
            target_frame=self.detector_frame,
            expected_source=self.specimen_frame,
            label="specimen_to_detector",
        )
        self._validate_transform(
            self.specimen_to_laboratory,
            target_frame=self.laboratory_frame,
            expected_source=self.specimen_frame,
            label="specimen_to_laboratory",
        )
        self._validate_transform(
            self.laboratory_to_reciprocal,
            target_frame=self.reciprocal_frame,
            expected_source=self.laboratory_frame,
            label="laboratory_to_reciprocal",
        )
        object.__setattr__(
            self,
            "metadata",
            MappingProxyType({str(key): str(value) for key, value in self.metadata.items()}),
        )
        object.__setattr__(self, "modality", modality)

    def _validate_transform(
        self,
        transform: FrameTransform | None,
        *,
        target_frame: ReferenceFrame | None,
        expected_source: ReferenceFrame | None,
        label: str,
    ) -> None:
        if transform is None:
            if (
                target_frame is not None
                and expected_source is not None
                and target_frame != expected_source
            ):
                raise ValueError(
                    f"AcquisitionGeometry.{label} must be provided when the "
                    "corresponding frames differ."
                )
            return
        if expected_source is None:
            raise ValueError(
                f"AcquisitionGeometry.{label} cannot be provided without its source frame."
            )
        if transform.source != expected_source:
            raise ValueError(
                f"AcquisitionGeometry.{label}.source must match the expected source frame."
            )
        if target_frame is None:
            raise ValueError(
                f"AcquisitionGeometry.{label} requires the corresponding target frame."
            )
        if transform.target != target_frame:
            raise ValueError(
                f"AcquisitionGeometry.{label}.target must match the corresponding target frame."
            )

    @property
    def supports_mapping(self) -> bool:
        return self.map_frame is not None

    @property
    def supports_detection(self) -> bool:
        return self.detector_frame is not None

    @property
    def supports_laboratory_context(self) -> bool:
        return self.laboratory_frame is not None

    def specimen_to_reciprocal_transform(self) -> FrameTransform | None:
        if self.specimen_to_laboratory is None or self.laboratory_to_reciprocal is None:
            return None
        return self.laboratory_to_reciprocal.compose(self.specimen_to_laboratory)
