from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from pytex.core.conventions import PYTEX_CANONICAL_CONVENTIONS
from pytex.core.frames import ReferenceFrame
from pytex.core.orientation import OrientationSet
from pytex.core.symmetry import SymmetrySpec
from pytex.ebsd import CrystalMap

EBSD_IMPORT_MANIFEST_SCHEMA_ID = "pytex.ebsd_import_manifest"
EBSD_IMPORT_MANIFEST_SCHEMA_VERSION = "1.0.0"
_PYTEX_VERSION = "0.1.0.dev0"
_MISSING = object()


def _iso8601_now() -> str:
    return datetime.now(tz=UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _validate_mapping_string(
    mapping: dict[str, Any], required_keys: tuple[str, ...], name: str
) -> None:
    missing = [key for key in required_keys if key not in mapping]
    if missing:
        raise ValueError(f"{name} is missing required keys: {', '.join(missing)}")
    for key in required_keys:
        value = mapping[key]
        if not isinstance(value, str) or not value:
            raise ValueError(f"{name}.{key} must be a non-empty string.")


@dataclass(frozen=True, slots=True)
class EBSDImportManifest:
    source_system: str
    source_file: str
    phase_name: str
    point_group: str
    crystal_frame: dict[str, str]
    specimen_frame: dict[str, str]
    orientation_convention: str = "bunge"
    angle_unit: str = "degree"
    schema_id: str = EBSD_IMPORT_MANIFEST_SCHEMA_ID
    schema_version: str = EBSD_IMPORT_MANIFEST_SCHEMA_VERSION
    pytex_version: str = _PYTEX_VERSION
    canonical_convention_set: str = PYTEX_CANONICAL_CONVENTIONS.name
    created_at: str = field(default_factory=_iso8601_now)
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("source_system", "source_file", "phase_name", "point_group"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"EBSDImportManifest.{name} must be a non-empty string.")
        if self.schema_id != EBSD_IMPORT_MANIFEST_SCHEMA_ID:
            raise ValueError(
                "EBSDImportManifest.schema_id must match the stable schema identifier."
            )
        if self.schema_version != EBSD_IMPORT_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                "EBSDImportManifest.schema_version must match the stable schema version."
            )
        if self.angle_unit not in {"degree", "radian"}:
            raise ValueError("EBSDImportManifest.angle_unit must be either 'degree' or 'radian'.")
        _validate_mapping_string(
            self.crystal_frame,
            required_keys=("name", "axis_1", "axis_2", "axis_3"),
            name="crystal_frame",
        )
        _validate_mapping_string(
            self.specimen_frame,
            required_keys=("name", "axis_1", "axis_2", "axis_3"),
            name="specimen_frame",
        )
        if not isinstance(self.metadata, dict):
            raise ValueError("EBSDImportManifest.metadata must be a string-to-string mapping.")
        for key, value in self.metadata.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise ValueError("EBSDImportManifest.metadata keys and values must be strings.")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "pytex_version": self.pytex_version,
            "canonical_convention_set": self.canonical_convention_set,
            "created_at": self.created_at,
            "source_system": self.source_system,
            "source_file": self.source_file,
            "phase_name": self.phase_name,
            "point_group": self.point_group,
            "orientation_convention": self.orientation_convention,
            "angle_unit": self.angle_unit,
            "crystal_frame": dict(self.crystal_frame),
            "specimen_frame": dict(self.specimen_frame),
            "metadata": dict(self.metadata),
        }

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8"
        )
        return output_path

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> EBSDImportManifest:
        validate_ebsd_import_manifest(payload)
        return cls(
            schema_id=payload["schema_id"],
            schema_version=payload["schema_version"],
            pytex_version=payload["pytex_version"],
            canonical_convention_set=payload["canonical_convention_set"],
            created_at=payload["created_at"],
            source_system=payload["source_system"],
            source_file=payload["source_file"],
            phase_name=payload["phase_name"],
            point_group=payload["point_group"],
            orientation_convention=payload["orientation_convention"],
            angle_unit=payload["angle_unit"],
            crystal_frame=dict(payload["crystal_frame"]),
            specimen_frame=dict(payload["specimen_frame"]),
            metadata=dict(payload["metadata"]),
        )


@dataclass(frozen=True, slots=True)
class NormalizedEBSDDataset:
    crystal_map: CrystalMap
    manifest: EBSDImportManifest
    source_phase_aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        map_phase = self.crystal_map.orientations.phase
        if map_phase is not None and map_phase.name != self.manifest.phase_name:
            raise ValueError(
                "NormalizedEBSDDataset manifest phase_name must match "
                "crystal_map.orientations.phase.name."
            )
        object.__setattr__(self, "source_phase_aliases", tuple(self.source_phase_aliases))

    def write_manifest_json(self, path: str | Path) -> Path:
        return self.manifest.write_json(path)


def validate_ebsd_import_manifest(payload: dict[str, Any]) -> None:
    required = {
        "schema_id",
        "schema_version",
        "pytex_version",
        "canonical_convention_set",
        "created_at",
        "source_system",
        "source_file",
        "phase_name",
        "point_group",
        "orientation_convention",
        "angle_unit",
        "crystal_frame",
        "specimen_frame",
        "metadata",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"EBSD import manifest is missing required fields: {', '.join(missing)}")
    EBSDImportManifest(
        schema_id=payload["schema_id"],
        schema_version=payload["schema_version"],
        pytex_version=payload["pytex_version"],
        canonical_convention_set=payload["canonical_convention_set"],
        created_at=payload["created_at"],
        source_system=payload["source_system"],
        source_file=payload["source_file"],
        phase_name=payload["phase_name"],
        point_group=payload["point_group"],
        orientation_convention=payload["orientation_convention"],
        angle_unit=payload["angle_unit"],
        crystal_frame=dict(payload["crystal_frame"]),
        specimen_frame=dict(payload["specimen_frame"]),
        metadata=dict(payload["metadata"]),
    )


def manifest_schema_path() -> Path:
    return Path(__file__).resolve().parents[3] / "schemas" / "ebsd_import_manifest.schema.json"


def read_ebsd_import_manifest(path: str | Path) -> EBSDImportManifest:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("EBSD import manifest JSON must decode to an object.")
    return EBSDImportManifest.from_dict(payload)


def _extract_from_object_paths(obj: Any, paths: tuple[str, ...], *, required: bool = True) -> Any:
    for path in paths:
        current = obj
        found = True
        for part in path.split("."):
            if isinstance(current, dict):
                if part not in current:
                    found = False
                    break
                current = current[part]
            elif hasattr(current, part):
                current = getattr(current, part)
            else:
                found = False
                break
        if found and current is not None:
            return current
    if required:
        raise ValueError(f"Unable to extract any of the required paths: {', '.join(paths)}")
    return _MISSING


def _payload_from_object(
    obj: Any,
    *,
    coordinate_paths: tuple[str, ...],
    angle_paths: tuple[str, ...],
    point_group_paths: tuple[str, ...],
    phase_name_paths: tuple[str, ...],
    grid_shape_paths: tuple[str, ...] = (),
    step_size_paths: tuple[str, ...] = (),
    source_file_paths: tuple[str, ...] = (),
    convention_paths: tuple[str, ...] = (),
    angle_unit_paths: tuple[str, ...] = (),
    metadata_paths: tuple[str, ...] = (),
    phase_alias_paths: tuple[str, ...] = (),
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "coordinates": _extract_from_object_paths(obj, coordinate_paths),
        "point_group": _extract_from_object_paths(obj, point_group_paths),
        "phase_name": _extract_from_object_paths(obj, phase_name_paths),
    }
    angle_value = _extract_from_object_paths(obj, angle_paths)
    angle_key = (
        "euler_angles_deg" if any("deg" in path.lower() for path in angle_paths) else "euler_angles"
    )
    payload[angle_key] = angle_value
    optional_mappings = (
        ("grid_shape", grid_shape_paths),
        ("step_sizes", step_size_paths),
        ("source_file", source_file_paths),
        ("orientation_convention", convention_paths),
        ("angle_unit", angle_unit_paths),
        ("metadata", metadata_paths),
        ("phase_aliases", phase_alias_paths),
    )
    for key, paths in optional_mappings:
        if not paths:
            continue
        value = _extract_from_object_paths(obj, paths, required=False)
        if value is not _MISSING:
            payload[key] = value
    return payload


def _normalize_vendor_payload(
    payload: dict[str, Any],
    *,
    source_system: str,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
    angle_key_candidates: tuple[str, ...],
) -> NormalizedEBSDDataset:
    angle_key = next(
        (candidate for candidate in angle_key_candidates if candidate in payload), None
    )
    if angle_key is None:
        raise ValueError(
            f"{source_system} payload must contain one of: {', '.join(angle_key_candidates)}"
        )
    if "coordinates" not in payload:
        raise ValueError(f"{source_system} payload must contain 'coordinates'.")
    if "point_group" not in payload or "phase_name" not in payload:
        raise ValueError(f"{source_system} payload must contain 'point_group' and 'phase_name'.")
    convention = str(payload.get("orientation_convention", "bunge"))
    degrees = str(payload.get("angle_unit", "degree")) == "degree"
    symmetry = SymmetrySpec.from_point_group(
        str(payload["point_group"]), reference_frame=crystal_frame
    )
    orientations = OrientationSet.from_euler_angles(
        payload[angle_key],
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        symmetry=symmetry,
        convention=convention,
        degrees=degrees,
    )
    crystal_map = CrystalMap(
        coordinates=payload["coordinates"],
        orientations=orientations,
        map_frame=map_frame,
        grid_shape=tuple(payload["grid_shape"]) if "grid_shape" in payload else None,
        step_sizes=tuple(payload["step_sizes"]) if "step_sizes" in payload else None,
    )
    manifest = EBSDImportManifest(
        source_system=source_system,
        source_file=str(payload.get("source_file", f"{source_system}_in_memory")),
        phase_name=str(payload["phase_name"]),
        point_group=str(payload["point_group"]),
        orientation_convention=convention,
        angle_unit="degree" if degrees else "radian",
        crystal_frame={
            "name": crystal_frame.name,
            "axis_1": crystal_frame.axes[0],
            "axis_2": crystal_frame.axes[1],
            "axis_3": crystal_frame.axes[2],
        },
        specimen_frame={
            "name": specimen_frame.name,
            "axis_1": specimen_frame.axes[0],
            "axis_2": specimen_frame.axes[1],
            "axis_3": specimen_frame.axes[2],
        },
        metadata={str(key): str(value) for key, value in dict(payload.get("metadata", {})).items()},
    )
    aliases = tuple(str(alias) for alias in payload.get("phase_aliases", ()))
    return NormalizedEBSDDataset(
        crystal_map=crystal_map,
        manifest=manifest,
        source_phase_aliases=aliases,
    )


def normalize_kikuchipy_payload(
    payload: dict[str, Any],
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
) -> NormalizedEBSDDataset:
    return _normalize_vendor_payload(
        payload,
        source_system="kikuchipy",
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
        angle_key_candidates=("euler_angles_deg", "euler_angles"),
    )


def normalize_kikuchipy_dataset(
    dataset: Any,
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
) -> NormalizedEBSDDataset:
    payload = _payload_from_object(
        dataset,
        coordinate_paths=("coordinates", "xmap.coordinates", "xmap.xy"),
        angle_paths=(
            "euler_angles_deg",
            "euler_angles",
            "xmap.euler_angles_deg",
            "xmap.euler_angles",
            "xmap.rotations.euler_angles_deg",
            "xmap.rotations.euler_angles",
        ),
        point_group_paths=("point_group", "phase.point_group", "xmap.phase.point_group"),
        phase_name_paths=("phase_name", "phase.name", "xmap.phase.name"),
        grid_shape_paths=("grid_shape", "xmap.grid_shape"),
        step_size_paths=("step_sizes", "xmap.step_sizes"),
        source_file_paths=("source_file", "path", "xmap.source_file"),
        convention_paths=("orientation_convention", "xmap.orientation_convention"),
        angle_unit_paths=("angle_unit", "xmap.angle_unit"),
        metadata_paths=("metadata", "xmap.metadata"),
        phase_alias_paths=("phase_aliases",),
    )
    return normalize_kikuchipy_payload(
        payload,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
    )


def normalize_pyebsdindex_payload(
    payload: dict[str, Any],
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
) -> NormalizedEBSDDataset:
    return _normalize_vendor_payload(
        payload,
        source_system="pyebsdindex",
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
        angle_key_candidates=("phi1_Phi_phi2_deg", "euler_angles_deg", "euler_angles"),
    )


def normalize_pyebsdindex_result(
    result: Any,
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
) -> NormalizedEBSDDataset:
    payload = _payload_from_object(
        result,
        coordinate_paths=("coordinates", "scan.coordinates", "xy"),
        angle_paths=(
            "phi1_Phi_phi2_deg",
            "euler_angles_deg",
            "euler_angles",
            "scan.phi1_Phi_phi2_deg",
            "scan.euler_angles_deg",
            "scan.euler_angles",
        ),
        point_group_paths=("point_group", "phase.point_group", "scan.phase.point_group"),
        phase_name_paths=("phase_name", "phase.name", "scan.phase.name"),
        grid_shape_paths=("grid_shape", "scan.grid_shape"),
        step_size_paths=("step_sizes", "scan.step_sizes"),
        source_file_paths=("source_file", "path", "scan.source_file"),
        convention_paths=("orientation_convention", "scan.orientation_convention"),
        angle_unit_paths=("angle_unit", "scan.angle_unit"),
        metadata_paths=("metadata", "scan.metadata"),
        phase_alias_paths=("phase_aliases",),
    )
    return normalize_pyebsdindex_payload(
        payload,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
    )
