from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np

from pytex.core import (
    AcquisitionGeometry,
    AtomicSite,
    Basis,
    BasisKind,
    CalibrationRecord,
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    FrameTransform,
    Handedness,
    Lattice,
    MeasurementQuality,
    MillerIndex,
    Misorientation,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReciprocalLatticeVector,
    ReferenceFrame,
    Rotation,
    ScatteringSetup,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
)
from pytex.diffraction import (
    DiffractionGeometry,
    DiffractionPattern,
    PowderPattern,
    PowderReflection,
    RadiationSpec,
    SAEDPattern,
    SAEDSpot,
)

JSON_CONTRACT_SCHEMA_VERSION = "1.0.0"


def _as_float_list(array: np.ndarray) -> list[float] | list[list[float]] | list[list[list[float]]]:
    return cast(
        list[float] | list[list[float]] | list[list[list[float]]],
        np.asarray(array, dtype=np.float64).tolist(),
    )


def _as_int_list(array: np.ndarray) -> list[int] | list[list[int]]:
    return cast(list[int] | list[list[int]], np.asarray(array, dtype=np.int64).tolist())


def _serialize_provenance_record(provenance: ProvenanceRecord) -> dict[str, Any]:
    payload = _serialize_provenance(provenance)
    if payload is None:
        raise TypeError("ProvenanceRecord serialization unexpectedly returned None.")
    return payload


def _require_schema(payload: dict[str, Any], schema_id: str) -> None:
    if payload.get("schema_id") != schema_id:
        raise ValueError(
            f"Expected schema_id '{schema_id}', received '{payload.get('schema_id')}'."
        )
    if payload.get("schema_version") != JSON_CONTRACT_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported schema_version "
            f"'{payload.get('schema_version')}'. Expected '{JSON_CONTRACT_SCHEMA_VERSION}'."
        )


def _base_payload(schema_id: str) -> dict[str, Any]:
    return {"schema_id": schema_id, "schema_version": JSON_CONTRACT_SCHEMA_VERSION}


def _serialize_provenance(provenance: ProvenanceRecord | None) -> dict[str, Any] | None:
    if provenance is None:
        return None
    return {
        **_base_payload("pytex.core.provenance_record"),
        "source_system": provenance.source_system,
        "source_identifier": provenance.source_identifier,
        "source_path": provenance.source_path,
        "source_version": provenance.source_version,
        "imported_at": provenance.imported_at,
        "metadata": dict(provenance.metadata),
        "notes": list(provenance.notes),
    }


def _deserialize_provenance(payload: dict[str, Any] | None) -> ProvenanceRecord | None:
    if payload is None:
        return None
    _require_schema(payload, "pytex.core.provenance_record")
    return ProvenanceRecord(
        source_system=payload["source_system"],
        source_identifier=payload.get("source_identifier"),
        source_path=payload.get("source_path"),
        source_version=payload.get("source_version"),
        imported_at=payload.get("imported_at"),
        metadata=payload.get("metadata", {}),
        notes=tuple(payload.get("notes", [])),
    )


def _serialize_reference_frame(frame: ReferenceFrame) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.reference_frame"),
        "name": frame.name,
        "domain": frame.domain.value,
        "axes": list(frame.axes),
        "handedness": frame.handedness.value,
        "convention": frame.convention.name,
        "description": frame.description,
        "provenance": _serialize_provenance(frame.provenance),
        "metadata": dict(frame.metadata),
    }


def _deserialize_reference_frame(payload: dict[str, Any]) -> ReferenceFrame:
    _require_schema(payload, "pytex.core.reference_frame")
    return ReferenceFrame(
        name=payload["name"],
        domain=FrameDomain(payload["domain"]),
        axes=tuple(payload["axes"]),
        handedness=Handedness(payload["handedness"]),
        description=payload.get("description", ""),
        provenance=_deserialize_provenance(payload.get("provenance")),
        metadata=payload.get("metadata", {}),
    )


def _serialize_frame_transform(transform: FrameTransform) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.frame_transform"),
        "source": _serialize_reference_frame(transform.source),
        "target": _serialize_reference_frame(transform.target),
        "rotation_matrix": _as_float_list(transform.rotation_matrix),
        "translation_vector": _as_float_list(transform.translation_vector),
        "provenance": _serialize_provenance(transform.provenance),
    }


def _deserialize_frame_transform(payload: dict[str, Any]) -> FrameTransform:
    _require_schema(payload, "pytex.core.frame_transform")
    return FrameTransform(
        source=_deserialize_reference_frame(payload["source"]),
        target=_deserialize_reference_frame(payload["target"]),
        rotation_matrix=np.asarray(payload["rotation_matrix"], dtype=np.float64),
        translation_vector=np.asarray(payload["translation_vector"], dtype=np.float64),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_measurement_quality(record: MeasurementQuality) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.measurement_quality"),
        "confidence": record.confidence,
        "valid_fraction": record.valid_fraction,
        "masked_fraction": record.masked_fraction,
        "uncertainty": dict(record.uncertainty),
        "flags": list(record.flags),
        "provenance": _serialize_provenance(record.provenance),
    }


def _deserialize_measurement_quality(payload: dict[str, Any]) -> MeasurementQuality:
    _require_schema(payload, "pytex.core.measurement_quality")
    return MeasurementQuality(
        confidence=payload.get("confidence"),
        valid_fraction=payload.get("valid_fraction"),
        masked_fraction=payload.get("masked_fraction"),
        uncertainty=payload.get("uncertainty", {}),
        flags=tuple(payload.get("flags", [])),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_calibration_record(record: CalibrationRecord) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.calibration_record"),
        "source": record.source,
        "status": record.status,
        "residual_error": record.residual_error,
        "parameter_uncertainty": dict(record.parameter_uncertainty),
        "notes": list(record.notes),
        "provenance": _serialize_provenance(record.provenance),
    }


def _deserialize_calibration_record(payload: dict[str, Any]) -> CalibrationRecord:
    _require_schema(payload, "pytex.core.calibration_record")
    return CalibrationRecord(
        source=payload["source"],
        status=payload.get("status", "nominal"),
        residual_error=payload.get("residual_error"),
        parameter_uncertainty=payload.get("parameter_uncertainty", {}),
        notes=tuple(payload.get("notes", [])),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_scattering_setup(setup: ScatteringSetup) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.scattering_setup"),
        "laboratory_frame": _serialize_reference_frame(setup.laboratory_frame),
        "radiation_type": setup.radiation_type,
        "incident_beam_direction": _as_float_list(setup.incident_beam_direction),
        "wavelength_angstrom": setup.wavelength_angstrom,
        "beam_energy_kev": setup.beam_energy_kev,
        "scattering_mode": setup.scattering_mode,
        "provenance": _serialize_provenance(setup.provenance),
    }


def _deserialize_scattering_setup(payload: dict[str, Any]) -> ScatteringSetup:
    _require_schema(payload, "pytex.core.scattering_setup")
    return ScatteringSetup(
        laboratory_frame=_deserialize_reference_frame(payload["laboratory_frame"]),
        radiation_type=payload.get("radiation_type", "electron"),
        incident_beam_direction=np.asarray(payload["incident_beam_direction"], dtype=np.float64),
        wavelength_angstrom=payload.get("wavelength_angstrom"),
        beam_energy_kev=payload.get("beam_energy_kev"),
        scattering_mode=payload.get("scattering_mode", "elastic"),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_acquisition_geometry(geometry: AcquisitionGeometry) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.acquisition_geometry"),
        "specimen_frame": _serialize_reference_frame(geometry.specimen_frame),
        "modality": geometry.modality,
        "map_frame": None
        if geometry.map_frame is None
        else _serialize_reference_frame(geometry.map_frame),
        "detector_frame": None
        if geometry.detector_frame is None
        else _serialize_reference_frame(geometry.detector_frame),
        "laboratory_frame": None
        if geometry.laboratory_frame is None
        else _serialize_reference_frame(geometry.laboratory_frame),
        "reciprocal_frame": None
        if geometry.reciprocal_frame is None
        else _serialize_reference_frame(geometry.reciprocal_frame),
        "specimen_to_map": None
        if geometry.specimen_to_map is None
        else _serialize_frame_transform(geometry.specimen_to_map),
        "specimen_to_detector": None
        if geometry.specimen_to_detector is None
        else _serialize_frame_transform(geometry.specimen_to_detector),
        "specimen_to_laboratory": None
        if geometry.specimen_to_laboratory is None
        else _serialize_frame_transform(geometry.specimen_to_laboratory),
        "laboratory_to_reciprocal": None
        if geometry.laboratory_to_reciprocal is None
        else _serialize_frame_transform(geometry.laboratory_to_reciprocal),
        "calibration_record": None
        if geometry.calibration_record is None
        else _serialize_calibration_record(geometry.calibration_record),
        "measurement_quality": None
        if geometry.measurement_quality is None
        else _serialize_measurement_quality(geometry.measurement_quality),
        "metadata": dict(geometry.metadata),
        "provenance": _serialize_provenance(geometry.provenance),
    }


def _deserialize_acquisition_geometry(payload: dict[str, Any]) -> AcquisitionGeometry:
    _require_schema(payload, "pytex.core.acquisition_geometry")
    return AcquisitionGeometry(
        specimen_frame=_deserialize_reference_frame(payload["specimen_frame"]),
        modality=payload.get("modality", "generic"),
        map_frame=None
        if payload.get("map_frame") is None
        else _deserialize_reference_frame(payload["map_frame"]),
        detector_frame=None
        if payload.get("detector_frame") is None
        else _deserialize_reference_frame(payload["detector_frame"]),
        laboratory_frame=None
        if payload.get("laboratory_frame") is None
        else _deserialize_reference_frame(payload["laboratory_frame"]),
        reciprocal_frame=None
        if payload.get("reciprocal_frame") is None
        else _deserialize_reference_frame(payload["reciprocal_frame"]),
        specimen_to_map=None
        if payload.get("specimen_to_map") is None
        else _deserialize_frame_transform(payload["specimen_to_map"]),
        specimen_to_detector=None
        if payload.get("specimen_to_detector") is None
        else _deserialize_frame_transform(payload["specimen_to_detector"]),
        specimen_to_laboratory=None
        if payload.get("specimen_to_laboratory") is None
        else _deserialize_frame_transform(payload["specimen_to_laboratory"]),
        laboratory_to_reciprocal=None
        if payload.get("laboratory_to_reciprocal") is None
        else _deserialize_frame_transform(payload["laboratory_to_reciprocal"]),
        calibration_record=None
        if payload.get("calibration_record") is None
        else _deserialize_calibration_record(payload["calibration_record"]),
        measurement_quality=None
        if payload.get("measurement_quality") is None
        else _deserialize_measurement_quality(payload["measurement_quality"]),
        metadata=payload.get("metadata", {}),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_symmetry(symmetry: SymmetrySpec) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.symmetry_spec"),
        "name": symmetry.name,
        "point_group": symmetry.point_group,
        "specimen_symmetry": symmetry.specimen_symmetry,
        "reference_frame": None
        if symmetry.reference_frame is None
        else _serialize_reference_frame(symmetry.reference_frame),
        "operators": _as_float_list(symmetry.operators),
        "provenance": _serialize_provenance(symmetry.provenance),
    }


def _deserialize_symmetry(payload: dict[str, Any]) -> SymmetrySpec:
    _require_schema(payload, "pytex.core.symmetry_spec")
    return SymmetrySpec(
        name=payload["name"],
        point_group=payload["point_group"],
        operators=np.asarray(payload["operators"], dtype=np.float64),
        specimen_symmetry=payload.get("specimen_symmetry"),
        reference_frame=None
        if payload.get("reference_frame") is None
        else _deserialize_reference_frame(payload["reference_frame"]),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_space_group(space_group: SpaceGroupSpec) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.space_group_spec"),
        "symbol": space_group.symbol,
        "number": space_group.number,
        "reference_frame": _serialize_reference_frame(space_group.reference_frame),
        "setting": space_group.setting,
        "crystal_system": space_group.crystal_system,
        "provenance": _serialize_provenance(space_group.provenance),
    }


def _deserialize_space_group(payload: dict[str, Any]) -> SpaceGroupSpec:
    _require_schema(payload, "pytex.core.space_group_spec")
    return SpaceGroupSpec(
        symbol=payload["symbol"],
        number=int(payload["number"]),
        reference_frame=_deserialize_reference_frame(payload["reference_frame"]),
        setting=payload.get("setting"),
        crystal_system=payload.get("crystal_system"),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_basis(basis: Basis) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.basis"),
        "frame": _serialize_reference_frame(basis.frame),
        "kind": basis.kind.value,
        "matrix": _as_float_list(basis.matrix),
        "unit": basis.unit,
    }


def _deserialize_basis(payload: dict[str, Any]) -> Basis:
    _require_schema(payload, "pytex.core.basis")
    return Basis(
        frame=_deserialize_reference_frame(payload["frame"]),
        kind=BasisKind(payload["kind"]),
        matrix=np.asarray(payload["matrix"], dtype=np.float64),
        unit=payload.get("unit", "angstrom"),
    )


def _serialize_lattice(lattice: Lattice) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.lattice"),
        "a": lattice.a,
        "b": lattice.b,
        "c": lattice.c,
        "alpha_deg": lattice.alpha_deg,
        "beta_deg": lattice.beta_deg,
        "gamma_deg": lattice.gamma_deg,
        "crystal_frame": _serialize_reference_frame(lattice.crystal_frame),
        "provenance": _serialize_provenance(lattice.provenance),
    }


def _deserialize_lattice(payload: dict[str, Any]) -> Lattice:
    _require_schema(payload, "pytex.core.lattice")
    return Lattice(
        a=float(payload["a"]),
        b=float(payload["b"]),
        c=float(payload["c"]),
        alpha_deg=float(payload["alpha_deg"]),
        beta_deg=float(payload["beta_deg"]),
        gamma_deg=float(payload["gamma_deg"]),
        crystal_frame=_deserialize_reference_frame(payload["crystal_frame"]),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_atomic_site(site: AtomicSite) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.atomic_site"),
        "label": site.label,
        "species": site.species,
        "fractional_coordinates": _as_float_list(site.fractional_coordinates),
        "occupancy": site.occupancy,
        "b_iso": site.b_iso,
    }


def _deserialize_atomic_site(payload: dict[str, Any]) -> AtomicSite:
    _require_schema(payload, "pytex.core.atomic_site")
    return AtomicSite(
        label=payload["label"],
        species=payload["species"],
        fractional_coordinates=np.asarray(payload["fractional_coordinates"], dtype=np.float64),
        occupancy=float(payload.get("occupancy", 1.0)),
        b_iso=payload.get("b_iso"),
    )


def _serialize_unit_cell(cell: UnitCell) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.unit_cell"),
        "lattice": _serialize_lattice(cell.lattice),
        "sites": [_serialize_atomic_site(site) for site in cell.sites],
        "provenance": _serialize_provenance(cell.provenance),
    }


def _deserialize_unit_cell(payload: dict[str, Any]) -> UnitCell:
    _require_schema(payload, "pytex.core.unit_cell")
    return UnitCell(
        lattice=_deserialize_lattice(payload["lattice"]),
        sites=tuple(
            _deserialize_atomic_site(site_payload) for site_payload in payload.get("sites", [])
        ),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_phase(phase: Phase) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.phase"),
        "name": phase.name,
        "lattice": _serialize_lattice(phase.lattice),
        "symmetry": _serialize_symmetry(phase.symmetry),
        "crystal_frame": _serialize_reference_frame(phase.crystal_frame),
        "unit_cell": None if phase.unit_cell is None else _serialize_unit_cell(phase.unit_cell),
        "space_group": None
        if phase.space_group is None
        else _serialize_space_group(phase.space_group),
        "space_group_symbol": phase.space_group_symbol,
        "space_group_number": phase.space_group_number,
        "chemical_formula": phase.chemical_formula,
        "aliases": list(phase.aliases),
        "provenance": _serialize_provenance(phase.provenance),
    }


def _deserialize_phase(payload: dict[str, Any]) -> Phase:
    _require_schema(payload, "pytex.core.phase")
    return Phase(
        name=payload["name"],
        lattice=_deserialize_lattice(payload["lattice"]),
        symmetry=_deserialize_symmetry(payload["symmetry"]),
        crystal_frame=_deserialize_reference_frame(payload["crystal_frame"]),
        unit_cell=None
        if payload.get("unit_cell") is None
        else _deserialize_unit_cell(payload["unit_cell"]),
        space_group=None
        if payload.get("space_group") is None
        else _deserialize_space_group(payload["space_group"]),
        space_group_symbol=payload.get("space_group_symbol"),
        space_group_number=payload.get("space_group_number"),
        chemical_formula=payload.get("chemical_formula"),
        aliases=tuple(payload.get("aliases", [])),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_miller_index(index: MillerIndex) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.miller_index"),
        "indices": _as_int_list(index.indices),
        "phase": _serialize_phase(index.phase),
        "basis_kind": index.basis_kind.value,
    }


def _deserialize_miller_index(payload: dict[str, Any]) -> MillerIndex:
    _require_schema(payload, "pytex.core.miller_index")
    return MillerIndex(
        indices=np.asarray(payload["indices"], dtype=np.int64),
        phase=_deserialize_phase(payload["phase"]),
        basis_kind=BasisKind(payload.get("basis_kind", "reciprocal")),
    )


def _serialize_crystal_direction(direction: CrystalDirection) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.crystal_direction"),
        "coordinates": _as_float_list(direction.coordinates),
        "phase": _serialize_phase(direction.phase),
    }


def _deserialize_crystal_direction(payload: dict[str, Any]) -> CrystalDirection:
    _require_schema(payload, "pytex.core.crystal_direction")
    return CrystalDirection(
        coordinates=np.asarray(payload["coordinates"], dtype=np.float64),
        phase=_deserialize_phase(payload["phase"]),
    )


def _serialize_crystal_plane(plane: CrystalPlane) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.crystal_plane"),
        "miller": _serialize_miller_index(plane.miller),
        "phase": _serialize_phase(plane.phase),
    }


def _deserialize_crystal_plane(payload: dict[str, Any]) -> CrystalPlane:
    _require_schema(payload, "pytex.core.crystal_plane")
    miller = _deserialize_miller_index(payload["miller"])
    phase = miller.phase
    if payload.get("phase") is not None:
        phase_payload = payload["phase"]
        if phase_payload != payload["miller"]["phase"]:
            phase = _deserialize_phase(phase_payload)
    return CrystalPlane(miller=miller, phase=phase)


def _serialize_reciprocal_vector(vector: ReciprocalLatticeVector) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.reciprocal_lattice_vector"),
        "coordinates": _as_float_list(vector.coordinates),
        "phase": _serialize_phase(vector.phase),
    }


def _deserialize_reciprocal_vector(payload: dict[str, Any]) -> ReciprocalLatticeVector:
    _require_schema(payload, "pytex.core.reciprocal_lattice_vector")
    return ReciprocalLatticeVector(
        coordinates=np.asarray(payload["coordinates"], dtype=np.float64),
        phase=_deserialize_phase(payload["phase"]),
    )


def _serialize_zone_axis(axis: ZoneAxis) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.zone_axis"),
        "indices": _as_int_list(axis.indices),
        "phase": _serialize_phase(axis.phase),
    }


def _deserialize_zone_axis(payload: dict[str, Any]) -> ZoneAxis:
    _require_schema(payload, "pytex.core.zone_axis")
    return ZoneAxis(
        indices=np.asarray(payload["indices"], dtype=np.int64),
        phase=_deserialize_phase(payload["phase"]),
    )


def _serialize_rotation(rotation: Rotation) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.rotation"),
        "quaternion_wxyz": _as_float_list(rotation.quaternion),
        "provenance": _serialize_provenance(rotation.provenance),
    }


def _deserialize_rotation(payload: dict[str, Any]) -> Rotation:
    _require_schema(payload, "pytex.core.rotation")
    return Rotation(
        quaternion=np.asarray(payload["quaternion_wxyz"], dtype=np.float64),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_orientation(orientation: Orientation) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.orientation"),
        "rotation": _serialize_rotation(orientation.rotation),
        "crystal_frame": _serialize_reference_frame(orientation.crystal_frame),
        "specimen_frame": _serialize_reference_frame(orientation.specimen_frame),
        "symmetry": None
        if orientation.symmetry is None
        else _serialize_symmetry(orientation.symmetry),
        "phase": None if orientation.phase is None else _serialize_phase(orientation.phase),
        "provenance": _serialize_provenance(orientation.provenance),
    }


def _deserialize_orientation(payload: dict[str, Any]) -> Orientation:
    _require_schema(payload, "pytex.core.orientation")
    phase = None if payload.get("phase") is None else _deserialize_phase(payload["phase"])
    symmetry = None
    if payload.get("symmetry") is not None:
        if phase is not None and payload["phase"]["symmetry"] == payload["symmetry"]:
            symmetry = phase.symmetry
        else:
            symmetry = _deserialize_symmetry(payload["symmetry"])
    return Orientation(
        rotation=_deserialize_rotation(payload["rotation"]),
        crystal_frame=_deserialize_reference_frame(payload["crystal_frame"]),
        specimen_frame=_deserialize_reference_frame(payload["specimen_frame"]),
        symmetry=symmetry,
        phase=phase,
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_misorientation(misorientation: Misorientation) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.misorientation"),
        "rotation": _serialize_rotation(misorientation.rotation),
        "left_symmetry": None
        if misorientation.left_symmetry is None
        else _serialize_symmetry(misorientation.left_symmetry),
        "right_symmetry": None
        if misorientation.right_symmetry is None
        else _serialize_symmetry(misorientation.right_symmetry),
        "provenance": _serialize_provenance(misorientation.provenance),
    }


def _deserialize_misorientation(payload: dict[str, Any]) -> Misorientation:
    _require_schema(payload, "pytex.core.misorientation")
    return Misorientation(
        rotation=_deserialize_rotation(payload["rotation"]),
        left_symmetry=None
        if payload.get("left_symmetry") is None
        else _deserialize_symmetry(payload["left_symmetry"]),
        right_symmetry=None
        if payload.get("right_symmetry") is None
        else _deserialize_symmetry(payload["right_symmetry"]),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_orientation_set(orientations: OrientationSet) -> dict[str, Any]:
    return {
        **_base_payload("pytex.core.orientation_set"),
        "quaternions_wxyz": _as_float_list(orientations.quaternions),
        "crystal_frame": _serialize_reference_frame(orientations.crystal_frame),
        "specimen_frame": _serialize_reference_frame(orientations.specimen_frame),
        "symmetry": None
        if orientations.symmetry is None
        else _serialize_symmetry(orientations.symmetry),
        "phase": None if orientations.phase is None else _serialize_phase(orientations.phase),
        "provenance": _serialize_provenance(orientations.provenance),
    }


def _deserialize_orientation_set(payload: dict[str, Any]) -> OrientationSet:
    _require_schema(payload, "pytex.core.orientation_set")
    phase = None if payload.get("phase") is None else _deserialize_phase(payload["phase"])
    symmetry = None
    if payload.get("symmetry") is not None:
        if phase is not None and payload["phase"]["symmetry"] == payload["symmetry"]:
            symmetry = phase.symmetry
        else:
            symmetry = _deserialize_symmetry(payload["symmetry"])
    return OrientationSet(
        quaternions=np.asarray(payload["quaternions_wxyz"], dtype=np.float64),
        crystal_frame=_deserialize_reference_frame(payload["crystal_frame"]),
        specimen_frame=_deserialize_reference_frame(payload["specimen_frame"]),
        symmetry=symmetry,
        phase=phase,
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_radiation(radiation: RadiationSpec) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.radiation_spec"),
        "name": radiation.name,
        "wavelength_angstrom": radiation.wavelength_angstrom,
    }


def _deserialize_radiation(payload: dict[str, Any]) -> RadiationSpec:
    _require_schema(payload, "pytex.diffraction.radiation_spec")
    return RadiationSpec(
        name=payload["name"], wavelength_angstrom=float(payload["wavelength_angstrom"])
    )


def _serialize_powder_reflection(reflection: PowderReflection) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.powder_reflection"),
        "miller_indices": _as_int_list(reflection.miller_indices),
        "d_spacing_angstrom": reflection.d_spacing_angstrom,
        "two_theta_deg": reflection.two_theta_deg,
        "intensity": reflection.intensity,
        "structure_factor_amplitude": reflection.structure_factor_amplitude,
        "multiplicity": reflection.multiplicity,
    }


def _deserialize_powder_reflection(payload: dict[str, Any]) -> PowderReflection:
    _require_schema(payload, "pytex.diffraction.powder_reflection")
    return PowderReflection(
        miller_indices=np.asarray(payload["miller_indices"], dtype=np.int64),
        d_spacing_angstrom=float(payload["d_spacing_angstrom"]),
        two_theta_deg=float(payload["two_theta_deg"]),
        intensity=float(payload["intensity"]),
        structure_factor_amplitude=float(payload["structure_factor_amplitude"]),
        multiplicity=int(payload["multiplicity"]),
    )


def _serialize_powder_pattern(pattern: PowderPattern) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.powder_pattern"),
        "phase": _serialize_phase(pattern.phase),
        "radiation": _serialize_radiation(pattern.radiation),
        "reflections": [
            _serialize_powder_reflection(reflection) for reflection in pattern.reflections
        ],
        "two_theta_grid_deg": _as_float_list(pattern.two_theta_grid_deg),
        "intensity_grid": _as_float_list(pattern.intensity_grid),
        "broadening_fwhm_deg": pattern.broadening_fwhm_deg,
        "provenance": _serialize_provenance(pattern.provenance),
    }


def _deserialize_powder_pattern(payload: dict[str, Any]) -> PowderPattern:
    _require_schema(payload, "pytex.diffraction.powder_pattern")
    return PowderPattern(
        phase=_deserialize_phase(payload["phase"]),
        radiation=_deserialize_radiation(payload["radiation"]),
        reflections=tuple(
            _deserialize_powder_reflection(item) for item in payload.get("reflections", [])
        ),
        two_theta_grid_deg=np.asarray(payload["two_theta_grid_deg"], dtype=np.float64),
        intensity_grid=np.asarray(payload["intensity_grid"], dtype=np.float64),
        broadening_fwhm_deg=payload.get("broadening_fwhm_deg"),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_saed_spot(spot: SAEDSpot) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.saed_spot"),
        "miller_indices": _as_int_list(spot.miller_indices),
        "reciprocal_vector_crystal": _as_float_list(spot.reciprocal_vector_crystal),
        "reciprocal_vector_detector": _as_float_list(spot.reciprocal_vector_detector),
        "detector_coordinates": _as_float_list(spot.detector_coordinates),
        "intensity": spot.intensity,
        "excitation_error_inv_angstrom": spot.excitation_error_inv_angstrom,
        "label": spot.label,
    }


def _deserialize_saed_spot(payload: dict[str, Any]) -> SAEDSpot:
    _require_schema(payload, "pytex.diffraction.saed_spot")
    return SAEDSpot(
        miller_indices=np.asarray(payload["miller_indices"], dtype=np.int64),
        reciprocal_vector_crystal=np.asarray(
            payload["reciprocal_vector_crystal"], dtype=np.float64
        ),
        reciprocal_vector_detector=np.asarray(
            payload["reciprocal_vector_detector"], dtype=np.float64
        ),
        detector_coordinates=np.asarray(payload["detector_coordinates"], dtype=np.float64),
        intensity=float(payload["intensity"]),
        excitation_error_inv_angstrom=float(payload["excitation_error_inv_angstrom"]),
        label=payload.get("label"),
    )


def _serialize_saed_pattern(pattern: SAEDPattern) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.saed_pattern"),
        "phase": _serialize_phase(pattern.phase),
        "zone_axis": _serialize_zone_axis(pattern.zone_axis),
        "detector_frame": _serialize_reference_frame(pattern.detector_frame),
        "reciprocal_frame": _serialize_reference_frame(pattern.reciprocal_frame),
        "camera_constant_mm_angstrom": pattern.camera_constant_mm_angstrom,
        "spots": [_serialize_saed_spot(spot) for spot in pattern.spots],
        "zone_basis_crystal": _as_float_list(pattern.zone_basis_crystal),
        "provenance": _serialize_provenance(pattern.provenance),
    }


def _deserialize_saed_pattern(payload: dict[str, Any]) -> SAEDPattern:
    _require_schema(payload, "pytex.diffraction.saed_pattern")
    return SAEDPattern(
        phase=_deserialize_phase(payload["phase"]),
        zone_axis=_deserialize_zone_axis(payload["zone_axis"]),
        detector_frame=_deserialize_reference_frame(payload["detector_frame"]),
        reciprocal_frame=_deserialize_reference_frame(payload["reciprocal_frame"]),
        camera_constant_mm_angstrom=float(payload["camera_constant_mm_angstrom"]),
        spots=tuple(_deserialize_saed_spot(item) for item in payload.get("spots", [])),
        zone_basis_crystal=np.asarray(payload["zone_basis_crystal"], dtype=np.float64),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_diffraction_geometry(geometry: DiffractionGeometry) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.diffraction_geometry"),
        "detector_frame": _serialize_reference_frame(geometry.detector_frame),
        "specimen_frame": _serialize_reference_frame(geometry.specimen_frame),
        "laboratory_frame": _serialize_reference_frame(geometry.laboratory_frame),
        "beam_energy_kev": geometry.beam_energy_kev,
        "camera_length_mm": geometry.camera_length_mm,
        "pattern_center": _as_float_list(geometry.pattern_center),
        "detector_pixel_size_um": list(geometry.detector_pixel_size_um),
        "detector_shape": list(geometry.detector_shape),
        "beam_direction_lab": _as_float_list(geometry.beam_direction_lab),
        "specimen_to_lab_matrix": _as_float_list(geometry.specimen_to_lab_matrix),
        "tilt_degrees": list(geometry.tilt_degrees),
        "acquisition_geometry": None
        if geometry.acquisition_geometry is None
        else _serialize_acquisition_geometry(geometry.acquisition_geometry),
        "calibration_record": None
        if geometry.calibration_record is None
        else _serialize_calibration_record(geometry.calibration_record),
        "measurement_quality": None
        if geometry.measurement_quality is None
        else _serialize_measurement_quality(geometry.measurement_quality),
        "scattering_setup": None
        if geometry.scattering_setup is None
        else _serialize_scattering_setup(geometry.scattering_setup),
        "provenance": _serialize_provenance(geometry.provenance),
    }


def _deserialize_diffraction_geometry(payload: dict[str, Any]) -> DiffractionGeometry:
    _require_schema(payload, "pytex.diffraction.diffraction_geometry")
    return DiffractionGeometry(
        detector_frame=_deserialize_reference_frame(payload["detector_frame"]),
        specimen_frame=_deserialize_reference_frame(payload["specimen_frame"]),
        laboratory_frame=_deserialize_reference_frame(payload["laboratory_frame"]),
        beam_energy_kev=float(payload["beam_energy_kev"]),
        camera_length_mm=float(payload["camera_length_mm"]),
        pattern_center=np.asarray(payload["pattern_center"], dtype=np.float64),
        detector_pixel_size_um=(
            float(payload["detector_pixel_size_um"][0]),
            float(payload["detector_pixel_size_um"][1]),
        ),
        detector_shape=(
            int(payload["detector_shape"][0]),
            int(payload["detector_shape"][1]),
        ),
        beam_direction_lab=np.asarray(payload["beam_direction_lab"], dtype=np.float64),
        specimen_to_lab_matrix=np.asarray(payload["specimen_to_lab_matrix"], dtype=np.float64),
        tilt_degrees=(
            float(payload.get("tilt_degrees", (0.0, 0.0, 0.0))[0]),
            float(payload.get("tilt_degrees", (0.0, 0.0, 0.0))[1]),
            float(payload.get("tilt_degrees", (0.0, 0.0, 0.0))[2]),
        ),
        acquisition_geometry=None
        if payload.get("acquisition_geometry") is None
        else _deserialize_acquisition_geometry(payload["acquisition_geometry"]),
        calibration_record=None
        if payload.get("calibration_record") is None
        else _deserialize_calibration_record(payload["calibration_record"]),
        measurement_quality=None
        if payload.get("measurement_quality") is None
        else _deserialize_measurement_quality(payload["measurement_quality"]),
        scattering_setup=None
        if payload.get("scattering_setup") is None
        else _deserialize_scattering_setup(payload["scattering_setup"]),
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


def _serialize_diffraction_pattern(pattern: DiffractionPattern) -> dict[str, Any]:
    return {
        **_base_payload("pytex.diffraction.diffraction_pattern"),
        "coordinates_px": _as_float_list(pattern.coordinates_px),
        "intensities": _as_float_list(pattern.intensities),
        "geometry": _serialize_diffraction_geometry(pattern.geometry),
        "phase": _serialize_phase(pattern.phase),
        "orientation": None
        if pattern.orientation is None
        else _serialize_orientation(pattern.orientation),
        "provenance": _serialize_provenance(pattern.provenance),
    }


def _deserialize_diffraction_pattern(payload: dict[str, Any]) -> DiffractionPattern:
    _require_schema(payload, "pytex.diffraction.diffraction_pattern")
    phase = _deserialize_phase(payload["phase"])
    orientation = None
    if payload.get("orientation") is not None:
        orientation_payload = payload["orientation"]
        if orientation_payload.get("phase") == payload["phase"]:
            symmetry = None
            if orientation_payload.get("symmetry") is not None:
                if orientation_payload["symmetry"] == payload["phase"]["symmetry"]:
                    symmetry = phase.symmetry
                else:
                    symmetry = _deserialize_symmetry(orientation_payload["symmetry"])
            orientation = Orientation(
                rotation=_deserialize_rotation(orientation_payload["rotation"]),
                crystal_frame=_deserialize_reference_frame(orientation_payload["crystal_frame"]),
                specimen_frame=_deserialize_reference_frame(orientation_payload["specimen_frame"]),
                symmetry=symmetry,
                phase=phase,
                provenance=_deserialize_provenance(orientation_payload.get("provenance")),
            )
        else:
            orientation = _deserialize_orientation(orientation_payload)
    return DiffractionPattern(
        coordinates_px=np.asarray(payload["coordinates_px"], dtype=np.float64),
        intensities=np.asarray(payload["intensities"], dtype=np.float64),
        geometry=_deserialize_diffraction_geometry(payload["geometry"]),
        phase=phase,
        orientation=orientation,
        provenance=_deserialize_provenance(payload.get("provenance")),
    )


_SERIALIZERS: dict[type[Any], Callable[[Any], dict[str, Any]]] = {
    ProvenanceRecord: _serialize_provenance_record,
    ReferenceFrame: _serialize_reference_frame,
    FrameTransform: _serialize_frame_transform,
    MeasurementQuality: _serialize_measurement_quality,
    CalibrationRecord: _serialize_calibration_record,
    ScatteringSetup: _serialize_scattering_setup,
    AcquisitionGeometry: _serialize_acquisition_geometry,
    SymmetrySpec: _serialize_symmetry,
    SpaceGroupSpec: _serialize_space_group,
    Basis: _serialize_basis,
    Lattice: _serialize_lattice,
    AtomicSite: _serialize_atomic_site,
    UnitCell: _serialize_unit_cell,
    Phase: _serialize_phase,
    MillerIndex: _serialize_miller_index,
    CrystalDirection: _serialize_crystal_direction,
    CrystalPlane: _serialize_crystal_plane,
    ReciprocalLatticeVector: _serialize_reciprocal_vector,
    ZoneAxis: _serialize_zone_axis,
    Rotation: _serialize_rotation,
    Orientation: _serialize_orientation,
    Misorientation: _serialize_misorientation,
    OrientationSet: _serialize_orientation_set,
    RadiationSpec: _serialize_radiation,
    PowderReflection: _serialize_powder_reflection,
    PowderPattern: _serialize_powder_pattern,
    SAEDSpot: _serialize_saed_spot,
    SAEDPattern: _serialize_saed_pattern,
    DiffractionGeometry: _serialize_diffraction_geometry,
    DiffractionPattern: _serialize_diffraction_pattern,
}

_DESERIALIZERS: dict[str, Callable[[dict[str, Any]], Any]] = {
    "pytex.core.provenance_record": _deserialize_provenance,
    "pytex.core.reference_frame": _deserialize_reference_frame,
    "pytex.core.frame_transform": _deserialize_frame_transform,
    "pytex.core.measurement_quality": _deserialize_measurement_quality,
    "pytex.core.calibration_record": _deserialize_calibration_record,
    "pytex.core.scattering_setup": _deserialize_scattering_setup,
    "pytex.core.acquisition_geometry": _deserialize_acquisition_geometry,
    "pytex.core.symmetry_spec": _deserialize_symmetry,
    "pytex.core.space_group_spec": _deserialize_space_group,
    "pytex.core.basis": _deserialize_basis,
    "pytex.core.lattice": _deserialize_lattice,
    "pytex.core.atomic_site": _deserialize_atomic_site,
    "pytex.core.unit_cell": _deserialize_unit_cell,
    "pytex.core.phase": _deserialize_phase,
    "pytex.core.miller_index": _deserialize_miller_index,
    "pytex.core.crystal_direction": _deserialize_crystal_direction,
    "pytex.core.crystal_plane": _deserialize_crystal_plane,
    "pytex.core.reciprocal_lattice_vector": _deserialize_reciprocal_vector,
    "pytex.core.zone_axis": _deserialize_zone_axis,
    "pytex.core.rotation": _deserialize_rotation,
    "pytex.core.orientation": _deserialize_orientation,
    "pytex.core.misorientation": _deserialize_misorientation,
    "pytex.core.orientation_set": _deserialize_orientation_set,
    "pytex.diffraction.radiation_spec": _deserialize_radiation,
    "pytex.diffraction.powder_reflection": _deserialize_powder_reflection,
    "pytex.diffraction.powder_pattern": _deserialize_powder_pattern,
    "pytex.diffraction.saed_spot": _deserialize_saed_spot,
    "pytex.diffraction.saed_pattern": _deserialize_saed_pattern,
    "pytex.diffraction.diffraction_geometry": _deserialize_diffraction_geometry,
    "pytex.diffraction.diffraction_pattern": _deserialize_diffraction_pattern,
}


def to_json_contract(obj: Any) -> dict[str, Any]:
    for cls, serializer in _SERIALIZERS.items():
        if isinstance(obj, cls):
            payload = serializer(obj)
            if payload is None:
                raise TypeError(f"Serializer for {cls.__name__} returned None.")
            return payload
    raise TypeError(f"No JSON contract serializer registered for {type(obj).__name__}.")


def from_json_contract(payload: dict[str, Any]) -> Any:
    schema_id = payload.get("schema_id")
    if not isinstance(schema_id, str):
        raise ValueError("JSON contract payload must contain a string schema_id.")
    deserializer = _DESERIALIZERS.get(schema_id)
    if deserializer is None:
        raise ValueError(f"No JSON contract deserializer registered for schema_id '{schema_id}'.")
    return deserializer(payload)


def write_json_contract(obj: Any, path: str | Path) -> Path:
    output_path = Path(path)
    output_path.write_text(
        json.dumps(to_json_contract(obj), indent=2, sort_keys=True), encoding="utf-8"
    )
    return output_path


def read_json_contract(path: str | Path) -> Any:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("JSON contract file must decode to an object.")
    return from_json_contract(payload)
