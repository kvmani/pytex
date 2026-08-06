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
    """Quality and completeness statistics of a measurement.

    Purpose
    -------
    Carries what is known about how good the data are — confidence, the
    fraction successfully indexed, the fraction masked out, and named
    warning flags — so that quality can be reported and used for weighting
    rather than being assumed uniform.

    Attributes
    ----------
    confidence : float, optional
        Overall confidence in ``[0, 1]``.
    valid_fraction : float, optional
        Fraction of points successfully measured.
    masked_fraction : float, optional
        Fraction excluded by masking.
    uncertainty : Mapping[str, float]
        Named uncertainty estimates.
    flags : tuple of str
        Named quality warnings.
    provenance : ProvenanceRecord, optional
    """

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
    """The calibration state of an instrument geometry.

    Purpose
    -------
    Distinguishes a measured calibration from an assumed nominal one. This
    matters for honest reporting: results derived from nominal geometry must
    not be presented as if the geometry had been calibrated, and
    :attr:`is_calibrated` is the flag that keeps the two apart.

    Attributes
    ----------
    status : str
        ``"nominal"`` for assumed instrument defaults, anything else for a
        real calibration.
    Remaining attributes record the method, date, and parameters of the
    calibration when there was one.
    """

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
        """Whether the record describes a real calibration rather than a nominal one.

        ``False`` for status ``"nominal"``, meaning the geometry was assumed from
        instrument defaults. Reports must not present nominal geometry as
        measured, so this distinction is carried explicitly rather than inferred.
        """

        return self.status != "nominal"


@dataclass(frozen=True, slots=True)
class PatternCenter:
    """The pattern centre of a diffraction detector, in fractional coordinates.

    Purpose
    -------
    The projection origin of an EBSD or TEM pattern. Vendors differ in how
    they define it — which corner is the origin, and how the detector
    distance is normalized — so the convention name is stored with the
    numbers rather than assumed.

    Attributes
    ----------
    x_fraction, y_fraction : float
        In-plane position as fractions of the detector extent.
    detector_distance_fraction : float
        Specimen-to-detector distance, as a fraction of the detector width.
    convention : str
        Which vendor convention the fractions follow.
    provenance : ProvenanceRecord, optional
    """

    x_fraction: float
    y_fraction: float
    detector_distance_fraction: float
    convention: str = "fractional_detector"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("x_fraction", self.x_fraction),
            ("y_fraction", self.y_fraction),
        ):
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"PatternCenter.{name} must lie in [0, 1].")
        if (
            not np.isfinite(self.detector_distance_fraction)
            or self.detector_distance_fraction <= 0.0
        ):
            raise ValueError("PatternCenter.detector_distance_fraction must be positive.")
        convention = self.convention.strip().lower()
        if not convention:
            raise ValueError("PatternCenter.convention must be non-empty.")
        object.__setattr__(self, "x_fraction", float(self.x_fraction))
        object.__setattr__(self, "y_fraction", float(self.y_fraction))
        object.__setattr__(
            self,
            "detector_distance_fraction",
            float(self.detector_distance_fraction),
        )
        object.__setattr__(self, "convention", convention)

    def as_array(self) -> np.ndarray:
        """The pattern centre as ``[x_fraction, y_fraction, distance_fraction]``.

        Fractional detector coordinates; read-only. Vendors differ in their
        pattern-centre conventions, which is why the convention name is stored
        alongside the numbers.
        """

        values = np.array(
            [self.x_fraction, self.y_fraction, self.detector_distance_fraction],
            dtype=np.float64,
        )
        values.setflags(write=False)
        return values


@dataclass(frozen=True, slots=True)
class EBSDDetectorGeometry:
    """The detector geometry of an EBSD acquisition.

    Purpose
    -------
    Detector frame, pattern centre, shape, and pixel size together — the
    geometry that turns a Kikuchi pattern into orientation information.

    Attributes
    ----------
    detector_frame : ReferenceFrame
        Must belong to the detector domain.
    pattern_center : PatternCenter
    Remaining attributes record the detector shape, pixel size, and
    specimen-to-detector distance.
    """

    detector_frame: ReferenceFrame
    pattern_center: PatternCenter
    detector_distance_mm: float
    pixel_size_um: tuple[float, float]
    detector_shape: tuple[int, int]
    tilt_degrees: tuple[float, float, float] = (0.0, 0.0, 0.0)
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.detector_frame.domain is not FrameDomain.DETECTOR:
            raise ValueError("EBSDDetectorGeometry.detector_frame must belong to detector domain.")
        if not np.isfinite(self.detector_distance_mm) or self.detector_distance_mm <= 0.0:
            raise ValueError("EBSDDetectorGeometry.detector_distance_mm must be positive.")
        if len(self.pixel_size_um) != 2 or any(value <= 0.0 for value in self.pixel_size_um):
            raise ValueError("EBSDDetectorGeometry.pixel_size_um must contain two positive values.")
        if len(self.detector_shape) != 2 or any(value <= 0 for value in self.detector_shape):
            raise ValueError("EBSDDetectorGeometry.detector_shape must contain two positive ints.")
        if len(self.tilt_degrees) != 3:
            raise ValueError("EBSDDetectorGeometry.tilt_degrees must contain three values.")
        object.__setattr__(self, "detector_distance_mm", float(self.detector_distance_mm))
        object.__setattr__(
            self,
            "pixel_size_um",
            tuple(float(value) for value in self.pixel_size_um),
        )
        object.__setattr__(
            self,
            "detector_shape",
            tuple(int(value) for value in self.detector_shape),
        )
        object.__setattr__(self, "tilt_degrees", tuple(float(value) for value in self.tilt_degrees))

    @property
    def pattern_center_array(self) -> np.ndarray:
        """The pattern centre as a fractional 3-vector; see
        :meth:`PatternCenter.as_array`.
        """

        return self.pattern_center.as_array()

    @property
    def pattern_center_px(self) -> np.ndarray:
        """The in-plane pattern centre in detector pixels, as ``(x, y)``.

        Converted from fractional coordinates using the detector shape.
        Read-only.
        """

        height, width = self.detector_shape
        values = np.array(
            [
                self.pattern_center.x_fraction * (width - 1),
                self.pattern_center.y_fraction * (height - 1),
            ],
            dtype=np.float64,
        )
        values.setflags(write=False)
        return values


@dataclass(frozen=True, slots=True)
class EBSDCalibrationGeometry:
    """An EBSD detector geometry together with its calibration record.

    Purpose
    -------
    Pairs the geometry with the evidence for it, so that a downstream
    consumer can tell measured geometry from assumed geometry without
    looking elsewhere.
    """

    acquisition_geometry: AcquisitionGeometry
    detector_geometry: EBSDDetectorGeometry
    map_to_specimen: FrameTransform | None = None
    calibration_record: CalibrationRecord | None = None
    measurement_quality: MeasurementQuality | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.acquisition_geometry.modality != "ebsd":
            raise ValueError("EBSDCalibrationGeometry.acquisition_geometry.modality must be ebsd.")
        if self.acquisition_geometry.detector_frame != self.detector_geometry.detector_frame:
            raise ValueError(
                "EBSDCalibrationGeometry detector frame must match acquisition geometry."
            )
        if self.map_to_specimen is not None:
            if self.acquisition_geometry.map_frame is None:
                raise ValueError("map_to_specimen requires acquisition_geometry.map_frame.")
            if self.map_to_specimen.source != self.acquisition_geometry.map_frame:
                raise ValueError("map_to_specimen.source must match acquisition map_frame.")
            if self.map_to_specimen.target != self.acquisition_geometry.specimen_frame:
                raise ValueError("map_to_specimen.target must match acquisition specimen_frame.")
        if (
            self.calibration_record is not None
            and self.detector_geometry.calibration_record is not None
            and self.calibration_record != self.detector_geometry.calibration_record
        ):
            raise ValueError("calibration_record must match detector geometry when both exist.")

    @property
    def pattern_center(self) -> PatternCenter:
        """The pattern centre of the underlying detector geometry.
        """

        return self.detector_geometry.pattern_center

    @property
    def detector_distance_mm(self) -> float:
        """Specimen-to-detector distance in millimetres.
        """

        return self.detector_geometry.detector_distance_mm


@dataclass(frozen=True, slots=True)
class ScatteringSetup:
    """The radiation and beam configuration of a scattering experiment.

    Purpose
    -------
    Declares what is doing the scattering: the radiation type, the incident
    beam direction in the laboratory frame, and the wavelength — either
    directly or as a beam energy from which the relativistically corrected
    electron wavelength is derived.

    Attributes
    ----------
    laboratory_frame : ReferenceFrame
        Must belong to the laboratory domain.
    radiation_type : str
        ``"electron"`` by default.
    incident_beam_direction : np.ndarray
        Beam direction in the laboratory frame.
    wavelength_angstrom : float, optional
        Explicit wavelength; takes precedence when given.
    beam_energy_kev : float, optional
        Accelerating voltage, from which the electron wavelength is derived.
        At least one of these two must be present.
    provenance : ProvenanceRecord, optional
    """

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
        """The wavelength to use, in angstroms.

        Returns the explicitly stored wavelength when there is one. Otherwise it
        derives the relativistically corrected electron wavelength from the beam
        energy, ``lambda = 12.2639 / sqrt(V (1 + 0.97845e-6 V))`` with ``V`` in
        volts. Raises when neither is available, rather than silently assuming a
        default that would misplace every simulated reflection.
        """

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
    """The complete frame and transform context of a measurement.

    Purpose
    -------
    The bridge across the canonical frame chain
    ``crystal -> specimen -> map -> detector -> laboratory -> reciprocal``.
    Any workflow that crosses a tool boundary must carry this, so the
    receiving side inherits the frames and transforms rather than assuming
    them.

    Not every workflow instantiates every frame; the ``supports_*``
    properties report which parts of the chain are actually declared, and
    composed transforms return ``None`` rather than inventing an identity
    when a leg is missing.

    Attributes
    ----------
    specimen_frame : ReferenceFrame
        The one required frame.
    modality : str
        ``"ebsd"``, ``"tem"``, ``"xrd"``, or ``"generic"``.
    map_frame, detector_frame, laboratory_frame : ReferenceFrame, optional
        The remaining domains, when the workflow reaches them.
    specimen_to_detector, specimen_to_laboratory : FrameTransform, optional
    laboratory_to_reciprocal : FrameTransform, optional
        The declared transforms between the frames above.
    calibration_record : CalibrationRecord, optional
    measurement_quality : MeasurementQuality, optional
    provenance : ProvenanceRecord, optional
    """

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
        """Whether a map-domain frame is declared, so mapped workflows are usable.
        """

        return self.map_frame is not None

    @property
    def supports_detection(self) -> bool:
        """Whether a detector-domain frame is declared.
        """

        return self.detector_frame is not None

    @property
    def supports_laboratory_context(self) -> bool:
        """Whether a laboratory-domain frame is declared.
        """

        return self.laboratory_frame is not None

    def specimen_to_reciprocal_transform(self) -> FrameTransform | None:
        """The composed specimen-to-reciprocal transform, when it can be formed.

        Composes specimen-to-laboratory with laboratory-to-reciprocal. Returns
        ``None`` when either leg is missing, rather than inventing an identity
        that would silently assert an unmeasured alignment.
        """

        if self.specimen_to_laboratory is None or self.laboratory_to_reciprocal is None:
            return None
        return self.laboratory_to_reciprocal.compose(self.specimen_to_laboratory)
