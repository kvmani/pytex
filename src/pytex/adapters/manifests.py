from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from pytex.core import (
    AcquisitionGeometry,
    CalibrationRecord,
    MeasurementQuality,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    ScatteringSetup,
)
from pytex.core.frames import FrameTransform

EXPERIMENT_MANIFEST_SCHEMA_ID = "pytex.experiment_manifest"
EXPERIMENT_MANIFEST_SCHEMA_VERSION = "1.0.0"
BENCHMARK_MANIFEST_SCHEMA_ID = "pytex.benchmark_manifest"
BENCHMARK_MANIFEST_SCHEMA_VERSION = "1.0.0"
VALIDATION_MANIFEST_SCHEMA_ID = "pytex.validation_manifest"
VALIDATION_MANIFEST_SCHEMA_VERSION = "1.0.0"
WORKFLOW_RESULT_MANIFEST_SCHEMA_ID = "pytex.workflow_result_manifest"
WORKFLOW_RESULT_MANIFEST_SCHEMA_VERSION = "1.0.0"
TRANSFORMATION_MANIFEST_SCHEMA_ID = "pytex.transformation_manifest"
TRANSFORMATION_MANIFEST_SCHEMA_VERSION = "1.0.0"
_PYTEX_VERSION = "0.1.0.dev0"


def _iso8601_now() -> str:
    return datetime.now(tz=UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _validate_non_empty_string(name: str, value: str) -> str:
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"{name} must be a non-empty string.")
    return normalized


def _string_mapping(values: Mapping[str, Any], *, name: str) -> dict[str, str]:
    normalized: dict[str, str] = {}
    for key, value in values.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(f"{name} keys and values must both be strings.")
        normalized[key] = value
    return normalized


def _string_tuple(values: Sequence[str], *, name: str) -> tuple[str, ...]:
    normalized = tuple(_validate_non_empty_string(name, str(value)) for value in values)
    return normalized


def _frame_declaration(frame: ReferenceFrame) -> dict[str, Any]:
    return {
        "name": frame.name,
        "domain": frame.domain.value,
        "axis_1": frame.axes[0],
        "axis_2": frame.axes[1],
        "axis_3": frame.axes[2],
        "handedness": frame.handedness.value,
        "convention": frame.convention.name,
    }


def _transform_declaration(transform: FrameTransform) -> dict[str, Any]:
    return {
        "source": transform.source.name,
        "target": transform.target.name,
        "rotation_matrix": transform.rotation_matrix.tolist(),
        "translation_vector": transform.translation_vector.tolist(),
    }


def _provenance_declaration(provenance: ProvenanceRecord | None) -> dict[str, Any] | None:
    if provenance is None:
        return None
    return {
        "source_system": provenance.source_system,
        "source_identifier": provenance.source_identifier,
        "source_path": provenance.source_path,
        "source_version": provenance.source_version,
        "imported_at": provenance.imported_at,
        "metadata": dict(provenance.metadata),
        "notes": list(provenance.notes),
    }


def _calibration_declaration(record: CalibrationRecord | None) -> dict[str, Any] | None:
    if record is None:
        return None
    return {
        "source": record.source,
        "status": record.status,
        "residual_error": record.residual_error,
        "parameter_uncertainty": dict(record.parameter_uncertainty),
        "notes": list(record.notes),
        "provenance": _provenance_declaration(record.provenance),
    }


def _measurement_quality_declaration(record: MeasurementQuality | None) -> dict[str, Any] | None:
    if record is None:
        return None
    return {
        "confidence": record.confidence,
        "valid_fraction": record.valid_fraction,
        "masked_fraction": record.masked_fraction,
        "uncertainty": dict(record.uncertainty),
        "flags": list(record.flags),
        "provenance": _provenance_declaration(record.provenance),
    }


def _scattering_setup_declaration(setup: ScatteringSetup | None) -> dict[str, Any] | None:
    if setup is None:
        return None
    return {
        "laboratory_frame": _frame_declaration(setup.laboratory_frame),
        "radiation_type": setup.radiation_type,
        "incident_beam_direction": setup.incident_beam_direction.tolist(),
        "wavelength_angstrom": setup.wavelength_angstrom,
        "beam_energy_kev": setup.beam_energy_kev,
        "scattering_mode": setup.scattering_mode,
        "provenance": _provenance_declaration(setup.provenance),
    }


def _phase_declaration(phase: Phase | None) -> dict[str, Any] | None:
    if phase is None:
        return None
    return {
        "name": phase.name,
        "point_group": phase.symmetry.point_group,
        "space_group_symbol": phase.space_group_symbol,
        "space_group_number": phase.space_group_number,
        "chemical_formula": phase.chemical_formula,
        "crystal_frame": _frame_declaration(phase.crystal_frame),
        "provenance": _provenance_declaration(phase.provenance),
    }


def _mapping_declaration(value: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if value is None:
        return None
    return dict(value)


def _mapping_sequence_declaration(
    values: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    return tuple(dict(value) for value in values)


@dataclass(frozen=True, slots=True)
class ExperimentManifest:
    source_system: str
    modality: str
    specimen_frame: Mapping[str, Any]
    map_frame: Mapping[str, Any] | None = None
    detector_frame: Mapping[str, Any] | None = None
    laboratory_frame: Mapping[str, Any] | None = None
    reciprocal_frame: Mapping[str, Any] | None = None
    specimen_to_map: Mapping[str, Any] | None = None
    specimen_to_detector: Mapping[str, Any] | None = None
    specimen_to_laboratory: Mapping[str, Any] | None = None
    laboratory_to_reciprocal: Mapping[str, Any] | None = None
    calibration_record: Mapping[str, Any] | None = None
    measurement_quality: Mapping[str, Any] | None = None
    scattering_setup: Mapping[str, Any] | None = None
    phase: Mapping[str, Any] | None = None
    phases: tuple[Mapping[str, Any], ...] = ()
    referenced_files: tuple[str, ...] = ()
    schema_id: str = EXPERIMENT_MANIFEST_SCHEMA_ID
    schema_version: str = EXPERIMENT_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    canonical_convention_set: str = "pytex_canonical"
    created_at: str = field(default_factory=_iso8601_now)
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "source_system", _validate_non_empty_string("source_system", self.source_system)
        )
        object.__setattr__(self, "modality", _validate_non_empty_string("modality", self.modality))
        if self.schema_id != EXPERIMENT_MANIFEST_SCHEMA_ID:
            raise ValueError(
                "ExperimentManifest.schema_id must match the stable schema identifier."
            )
        if self.schema_version != EXPERIMENT_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "ExperimentManifest.schema_version must match the stable schema version."
            )
        object.__setattr__(
            self, "referenced_files", _string_tuple(self.referenced_files, name="referenced_files")
        )
        object.__setattr__(
            self, "metadata", _string_mapping(self.metadata, name="ExperimentManifest.metadata")
        )
        if self.phase is not None:
            object.__setattr__(self, "phase", _mapping_declaration(self.phase))
        object.__setattr__(self, "phases", _mapping_sequence_declaration(self.phases))

    @classmethod
    def from_acquisition_geometry(
        cls,
        acquisition_geometry: AcquisitionGeometry,
        *,
        source_system: str,
        phase: Phase | None = None,
        phases: Sequence[Phase] = (),
        scattering_setup: ScatteringSetup | None = None,
        referenced_files: Sequence[str] = (),
        metadata: Mapping[str, str] | None = None,
    ) -> ExperimentManifest:
        phase_declarations = tuple(
            declaration
            for declaration in (_phase_declaration(candidate) for candidate in phases)
            if declaration is not None
        )
        return cls(
            source_system=source_system,
            modality=acquisition_geometry.modality,
            specimen_frame=_frame_declaration(acquisition_geometry.specimen_frame),
            map_frame=None
            if acquisition_geometry.map_frame is None
            else _frame_declaration(acquisition_geometry.map_frame),
            detector_frame=None
            if acquisition_geometry.detector_frame is None
            else _frame_declaration(acquisition_geometry.detector_frame),
            laboratory_frame=None
            if acquisition_geometry.laboratory_frame is None
            else _frame_declaration(acquisition_geometry.laboratory_frame),
            reciprocal_frame=None
            if acquisition_geometry.reciprocal_frame is None
            else _frame_declaration(acquisition_geometry.reciprocal_frame),
            specimen_to_map=None
            if acquisition_geometry.specimen_to_map is None
            else _transform_declaration(acquisition_geometry.specimen_to_map),
            specimen_to_detector=None
            if acquisition_geometry.specimen_to_detector is None
            else _transform_declaration(acquisition_geometry.specimen_to_detector),
            specimen_to_laboratory=None
            if acquisition_geometry.specimen_to_laboratory is None
            else _transform_declaration(acquisition_geometry.specimen_to_laboratory),
            laboratory_to_reciprocal=None
            if acquisition_geometry.laboratory_to_reciprocal is None
            else _transform_declaration(acquisition_geometry.laboratory_to_reciprocal),
            calibration_record=_calibration_declaration(acquisition_geometry.calibration_record),
            measurement_quality=_measurement_quality_declaration(
                acquisition_geometry.measurement_quality
            ),
            scattering_setup=_scattering_setup_declaration(scattering_setup),
            phase=_phase_declaration(phase),
            phases=phase_declarations,
            referenced_files=tuple(referenced_files),
            metadata={} if metadata is None else dict(metadata),
        )

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "canonical_convention_set": self.canonical_convention_set,
            "created_at": self.created_at,
            "source_system": self.source_system,
            "modality": self.modality,
            "specimen_frame": dict(self.specimen_frame),
            "referenced_files": list(self.referenced_files),
            "metadata": dict(self.metadata),
        }
        for key, value in (
            ("map_frame", self.map_frame),
            ("detector_frame", self.detector_frame),
            ("laboratory_frame", self.laboratory_frame),
            ("reciprocal_frame", self.reciprocal_frame),
            ("specimen_to_map", self.specimen_to_map),
            ("specimen_to_detector", self.specimen_to_detector),
            ("specimen_to_laboratory", self.specimen_to_laboratory),
            ("laboratory_to_reciprocal", self.laboratory_to_reciprocal),
            ("calibration_record", self.calibration_record),
            ("measurement_quality", self.measurement_quality),
            ("scattering_setup", self.scattering_setup),
            ("phase", self.phase),
        ):
            declaration = _mapping_declaration(value)
            if declaration is not None:
                payload[key] = declaration
        if self.phases:
            payload["phases"] = [dict(item) for item in self.phases]
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ExperimentManifest:
        validate_experiment_manifest(payload)
        return cls(**payload)

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path


@dataclass(frozen=True, slots=True)
class BenchmarkManifest:
    benchmark_id: str
    subsystem: str
    baseline_kind: str
    workflows: tuple[str, ...]
    fixture_ids: tuple[str, ...] = ()
    artifact_paths: tuple[str, ...] = ()
    tolerances: Mapping[str, float] = field(default_factory=dict)
    status: str = "active"
    notes: tuple[str, ...] = ()
    schema_id: str = BENCHMARK_MANIFEST_SCHEMA_ID
    schema_version: str = BENCHMARK_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    created_at: str = field(default_factory=_iso8601_now)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "benchmark_id", _validate_non_empty_string("benchmark_id", self.benchmark_id)
        )
        object.__setattr__(
            self, "subsystem", _validate_non_empty_string("subsystem", self.subsystem)
        )
        object.__setattr__(
            self, "baseline_kind", _validate_non_empty_string("baseline_kind", self.baseline_kind)
        )
        if self.schema_id != BENCHMARK_MANIFEST_SCHEMA_ID:
            raise ValueError("BenchmarkManifest.schema_id must match the stable schema identifier.")
        if self.schema_version != BENCHMARK_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "BenchmarkManifest.schema_version must match the stable schema version."
            )
        object.__setattr__(self, "workflows", _string_tuple(self.workflows, name="workflows"))
        object.__setattr__(self, "fixture_ids", _string_tuple(self.fixture_ids, name="fixture_ids"))
        object.__setattr__(
            self, "artifact_paths", _string_tuple(self.artifact_paths, name="artifact_paths")
        )
        object.__setattr__(
            self, "notes", _string_tuple(self.notes, name="notes") if self.notes else ()
        )
        tolerances = {str(key): float(value) for key, value in self.tolerances.items()}
        if any(value < 0.0 for value in tolerances.values()):
            raise ValueError("BenchmarkManifest tolerances must be non-negative.")
        object.__setattr__(self, "tolerances", tolerances)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "created_at": self.created_at,
            "benchmark_id": self.benchmark_id,
            "subsystem": self.subsystem,
            "baseline_kind": self.baseline_kind,
            "workflows": list(self.workflows),
            "fixture_ids": list(self.fixture_ids),
            "artifact_paths": list(self.artifact_paths),
            "tolerances": dict(self.tolerances),
            "status": self.status,
            "notes": list(self.notes),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> BenchmarkManifest:
        validate_benchmark_manifest(payload)
        return cls(**payload)

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path


@dataclass(frozen=True, slots=True)
class ValidationManifest:
    campaign_name: str
    subsystem: str
    baseline_kind: str
    status: str
    reference_ids: tuple[str, ...]
    fixture_ids: tuple[str, ...] = ()
    artifact_paths: tuple[str, ...] = ()
    linked_benchmark_ids: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    schema_id: str = VALIDATION_MANIFEST_SCHEMA_ID
    schema_version: str = VALIDATION_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    created_at: str = field(default_factory=_iso8601_now)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "campaign_name", _validate_non_empty_string("campaign_name", self.campaign_name)
        )
        object.__setattr__(
            self, "subsystem", _validate_non_empty_string("subsystem", self.subsystem)
        )
        object.__setattr__(
            self, "baseline_kind", _validate_non_empty_string("baseline_kind", self.baseline_kind)
        )
        object.__setattr__(self, "status", _validate_non_empty_string("status", self.status))
        if self.schema_id != VALIDATION_MANIFEST_SCHEMA_ID:
            raise ValueError(
                "ValidationManifest.schema_id must match the stable schema identifier."
            )
        if self.schema_version != VALIDATION_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "ValidationManifest.schema_version must match the stable schema version."
            )
        object.__setattr__(
            self, "reference_ids", _string_tuple(self.reference_ids, name="reference_ids")
        )
        object.__setattr__(self, "fixture_ids", _string_tuple(self.fixture_ids, name="fixture_ids"))
        object.__setattr__(
            self, "artifact_paths", _string_tuple(self.artifact_paths, name="artifact_paths")
        )
        object.__setattr__(
            self,
            "linked_benchmark_ids",
            _string_tuple(self.linked_benchmark_ids, name="linked_benchmark_ids"),
        )
        object.__setattr__(
            self, "notes", _string_tuple(self.notes, name="notes") if self.notes else ()
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "created_at": self.created_at,
            "campaign_name": self.campaign_name,
            "subsystem": self.subsystem,
            "baseline_kind": self.baseline_kind,
            "status": self.status,
            "reference_ids": list(self.reference_ids),
            "fixture_ids": list(self.fixture_ids),
            "artifact_paths": list(self.artifact_paths),
            "linked_benchmark_ids": list(self.linked_benchmark_ids),
            "notes": list(self.notes),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ValidationManifest:
        validate_validation_manifest(payload)
        return cls(**payload)

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path


@dataclass(frozen=True, slots=True)
class WorkflowResultManifest:
    result_id: str
    workflow_name: str
    modality: str
    produced_by: str
    input_manifest_ids: tuple[str, ...]
    artifact_paths: tuple[str, ...]
    notes: tuple[str, ...] = ()
    metadata: Mapping[str, str] = field(default_factory=dict)
    schema_id: str = WORKFLOW_RESULT_MANIFEST_SCHEMA_ID
    schema_version: str = WORKFLOW_RESULT_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    created_at: str = field(default_factory=_iso8601_now)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "result_id", _validate_non_empty_string("result_id", self.result_id)
        )
        object.__setattr__(
            self, "workflow_name", _validate_non_empty_string("workflow_name", self.workflow_name)
        )
        object.__setattr__(self, "modality", _validate_non_empty_string("modality", self.modality))
        object.__setattr__(
            self, "produced_by", _validate_non_empty_string("produced_by", self.produced_by)
        )
        if self.schema_id != WORKFLOW_RESULT_MANIFEST_SCHEMA_ID:
            raise ValueError(
                "WorkflowResultManifest.schema_id must match the stable schema identifier."
            )
        if self.schema_version != WORKFLOW_RESULT_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "WorkflowResultManifest.schema_version must match the stable schema version."
            )
        object.__setattr__(
            self,
            "input_manifest_ids",
            _string_tuple(self.input_manifest_ids, name="input_manifest_ids"),
        )
        object.__setattr__(
            self, "artifact_paths", _string_tuple(self.artifact_paths, name="artifact_paths")
        )
        object.__setattr__(
            self, "notes", _string_tuple(self.notes, name="notes") if self.notes else ()
        )
        object.__setattr__(
            self, "metadata", _string_mapping(self.metadata, name="WorkflowResultManifest.metadata")
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "created_at": self.created_at,
            "result_id": self.result_id,
            "workflow_name": self.workflow_name,
            "modality": self.modality,
            "produced_by": self.produced_by,
            "input_manifest_ids": list(self.input_manifest_ids),
            "artifact_paths": list(self.artifact_paths),
            "notes": list(self.notes),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> WorkflowResultManifest:
        validate_workflow_result_manifest(payload)
        return cls(**payload)

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path


@dataclass(frozen=True, slots=True)
class TransformationManifest:
    transformation_id: str
    orientation_relationship_name: str
    parent_phase: Mapping[str, Any]
    child_phase: Mapping[str, Any]
    specimen_frame: Mapping[str, Any]
    variant_count: int
    observed_child_count: int
    has_variant_assignments: bool
    referenced_files: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    metadata: Mapping[str, str] = field(default_factory=dict)
    schema_id: str = TRANSFORMATION_MANIFEST_SCHEMA_ID
    schema_version: str = TRANSFORMATION_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    created_at: str = field(default_factory=_iso8601_now)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "transformation_id",
            _validate_non_empty_string("transformation_id", self.transformation_id),
        )
        object.__setattr__(
            self,
            "orientation_relationship_name",
            _validate_non_empty_string(
                "orientation_relationship_name", self.orientation_relationship_name
            ),
        )
        if self.schema_id != TRANSFORMATION_MANIFEST_SCHEMA_ID:
            raise ValueError(
                "TransformationManifest.schema_id must match the stable schema identifier."
            )
        if self.schema_version != TRANSFORMATION_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "TransformationManifest.schema_version must match the stable schema version."
            )
        if self.variant_count < 0:
            raise ValueError("TransformationManifest.variant_count must be non-negative.")
        if self.observed_child_count < 0:
            raise ValueError("TransformationManifest.observed_child_count must be non-negative.")
        object.__setattr__(
            self,
            "referenced_files",
            _string_tuple(self.referenced_files, name="referenced_files"),
        )
        object.__setattr__(
            self,
            "notes",
            _string_tuple(self.notes, name="notes") if self.notes else (),
        )
        object.__setattr__(
            self,
            "metadata",
            _string_mapping(self.metadata, name="TransformationManifest.metadata"),
        )

    @classmethod
    def from_phase_transformation_record(
        cls,
        record: Any,
        *,
        referenced_files: Sequence[str] = (),
        metadata: Mapping[str, str] | None = None,
    ) -> TransformationManifest:
        return cls(
            transformation_id=record.name,
            orientation_relationship_name=record.orientation_relationship.name,
            parent_phase=_phase_declaration(record.orientation_relationship.parent_phase) or {},
            child_phase=_phase_declaration(record.orientation_relationship.child_phase) or {},
            specimen_frame=_frame_declaration(record.parent_orientation.specimen_frame),
            variant_count=record.variant_count,
            observed_child_count=len(record.child_orientations),
            has_variant_assignments=record.variant_indices is not None,
            referenced_files=tuple(referenced_files),
            notes=tuple(record.notes),
            metadata={} if metadata is None else dict(metadata),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "created_at": self.created_at,
            "transformation_id": self.transformation_id,
            "orientation_relationship_name": self.orientation_relationship_name,
            "parent_phase": dict(self.parent_phase),
            "child_phase": dict(self.child_phase),
            "specimen_frame": dict(self.specimen_frame),
            "variant_count": self.variant_count,
            "observed_child_count": self.observed_child_count,
            "has_variant_assignments": self.has_variant_assignments,
            "referenced_files": list(self.referenced_files),
            "notes": list(self.notes),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> TransformationManifest:
        validate_transformation_manifest(payload)
        return cls(**payload)

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path


def _require_fields(payload: Mapping[str, Any], required: set[str], *, name: str) -> None:
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"{name} is missing required fields: {', '.join(missing)}")


def validate_experiment_manifest(payload: dict[str, Any]) -> None:
    _require_fields(
        payload,
        {
            "schema_id",
            "schema_version",
            "pytex_version",
            "canonical_convention_set",
            "created_at",
            "source_system",
            "modality",
            "specimen_frame",
            "referenced_files",
            "metadata",
        },
        name="Experiment manifest",
    )
    ExperimentManifest(**payload)


def validate_benchmark_manifest(payload: dict[str, Any]) -> None:
    _require_fields(
        payload,
        {
            "schema_id",
            "schema_version",
            "pytex_version",
            "created_at",
            "benchmark_id",
            "subsystem",
            "baseline_kind",
            "workflows",
            "fixture_ids",
            "artifact_paths",
            "tolerances",
            "status",
            "notes",
        },
        name="Benchmark manifest",
    )
    BenchmarkManifest(**payload)


def validate_validation_manifest(payload: dict[str, Any]) -> None:
    _require_fields(
        payload,
        {
            "schema_id",
            "schema_version",
            "pytex_version",
            "created_at",
            "campaign_name",
            "subsystem",
            "baseline_kind",
            "status",
            "reference_ids",
            "fixture_ids",
            "artifact_paths",
            "linked_benchmark_ids",
            "notes",
        },
        name="Validation manifest",
    )
    ValidationManifest(**payload)


def validate_workflow_result_manifest(payload: dict[str, Any]) -> None:
    _require_fields(
        payload,
        {
            "schema_id",
            "schema_version",
            "pytex_version",
            "created_at",
            "result_id",
            "workflow_name",
            "modality",
            "produced_by",
            "input_manifest_ids",
            "artifact_paths",
            "notes",
            "metadata",
        },
        name="Workflow result manifest",
    )
    WorkflowResultManifest(**payload)


def validate_transformation_manifest(payload: dict[str, Any]) -> None:
    _require_fields(
        payload,
        {
            "schema_id",
            "schema_version",
            "pytex_version",
            "created_at",
            "transformation_id",
            "orientation_relationship_name",
            "parent_phase",
            "child_phase",
            "specimen_frame",
            "variant_count",
            "observed_child_count",
            "has_variant_assignments",
            "referenced_files",
            "notes",
            "metadata",
        },
        name="Transformation manifest",
    )
    TransformationManifest(**payload)


def _schema_path(filename: str) -> Path:
    return Path(__file__).resolve().parents[3] / "schemas" / filename


def experiment_manifest_schema_path() -> Path:
    return _schema_path("experiment_manifest.schema.json")


def benchmark_manifest_schema_path() -> Path:
    return _schema_path("benchmark_manifest.schema.json")


def validation_manifest_schema_path() -> Path:
    return _schema_path("validation_manifest.schema.json")


def workflow_result_manifest_schema_path() -> Path:
    return _schema_path("workflow_result_manifest.schema.json")


def transformation_manifest_schema_path() -> Path:
    return _schema_path("transformation_manifest.schema.json")


def _read_manifest(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Manifest JSON must decode to an object.")
    return payload


def read_experiment_manifest(path: str | Path) -> ExperimentManifest:
    return ExperimentManifest.from_dict(_read_manifest(path))


def read_benchmark_manifest(path: str | Path) -> BenchmarkManifest:
    return BenchmarkManifest.from_dict(_read_manifest(path))


def read_validation_manifest(path: str | Path) -> ValidationManifest:
    return ValidationManifest.from_dict(_read_manifest(path))


def read_workflow_result_manifest(path: str | Path) -> WorkflowResultManifest:
    return WorkflowResultManifest.from_dict(_read_manifest(path))


def read_transformation_manifest(path: str | Path) -> TransformationManifest:
    return TransformationManifest.from_dict(_read_manifest(path))
