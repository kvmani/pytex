from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pytex import (
    AcquisitionGeometry,
    BenchmarkManifest,
    CalibrationRecord,
    CrystalMap,
    DiffractionGeometry,
    ExperimentManifest,
    FrameDomain,
    FrameTransform,
    Handedness,
    Lattice,
    MeasurementQuality,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    ScatteringSetup,
    SymmetrySpec,
    ValidationManifest,
    WorkflowResultManifest,
    benchmark_manifest_schema_path,
    experiment_manifest_schema_path,
    read_benchmark_manifest,
    read_experiment_manifest,
    read_validation_manifest,
    read_workflow_result_manifest,
    validation_manifest_schema_path,
    validate_benchmark_manifest,
    validate_validation_manifest,
    validate_workflow_result_manifest,
    workflow_result_manifest_schema_path,
)


def make_foundation() -> (
    tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame, ReferenceFrame, Phase]
):
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab",
        domain=FrameDomain.LABORATORY,
        axes=("X", "Y", "Z"),
        handedness=Handedness.RIGHT,
    )
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    phase = Phase(name="fcc", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, detector, lab, phase


def test_experiment_manifest_round_trip_for_acquisition_geometry(tmp_path: Path) -> None:
    crystal, specimen, detector, lab, phase = make_foundation()
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="ebsd",
        detector_frame=detector,
        laboratory_frame=lab,
        specimen_to_detector=FrameTransform(
            source=specimen, target=detector, rotation_matrix=np.eye(3)
        ),
        specimen_to_laboratory=FrameTransform(
            source=specimen, target=lab, rotation_matrix=np.eye(3)
        ),
        calibration_record=CalibrationRecord(
            source="vendor-fit", status="calibrated", residual_error=0.1
        ),
        measurement_quality=MeasurementQuality(confidence=0.95, valid_fraction=1.0),
    )
    manifest = ExperimentManifest.from_acquisition_geometry(
        acquisition,
        source_system="pytex",
        phase=phase,
        referenced_files=("scan.ang",),
        metadata={"operator": "test"},
    )
    payload = manifest.to_dict()
    assert payload["schema_id"] == "pytex.experiment_manifest"
    path = tmp_path / "experiment.json"
    manifest.write_json(path)
    restored = read_experiment_manifest(path)
    assert restored.to_dict() == payload


def test_crystal_map_can_emit_experiment_manifest() -> None:
    crystal, specimen, _, _, phase = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            )
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0]]),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(1, 1),
        step_sizes=(1.0, 1.0),
        calibration_record=CalibrationRecord(source="import", status="nominal"),
        measurement_quality=MeasurementQuality(confidence=0.9, valid_fraction=1.0),
    )
    manifest = crystal_map.to_experiment_manifest(source_system="pytex")
    assert manifest.modality == "ebsd"
    assert manifest.phase is not None
    assert manifest.phase["name"] == "fcc"
    assert manifest.metadata["grid_shape"] == "1x1"


def test_diffraction_geometry_can_emit_experiment_manifest() -> None:
    _, specimen, detector, lab, phase = make_foundation()
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(512, 512),
        calibration_record=CalibrationRecord(source="pattern-fit", status="refined"),
        measurement_quality=MeasurementQuality(confidence=0.92, valid_fraction=0.98),
        scattering_setup=ScatteringSetup(laboratory_frame=lab, beam_energy_kev=200.0),
    )
    manifest = geometry.to_experiment_manifest(source_system="pytex", phase=phase)
    assert manifest.detector_frame is not None
    assert manifest.laboratory_frame is not None
    assert manifest.scattering_setup is not None
    assert manifest.phase is not None


def test_benchmark_validation_and_workflow_manifests_round_trip(tmp_path: Path) -> None:
    benchmark = BenchmarkManifest(
        benchmark_id="demo_benchmark",
        subsystem="diffraction",
        baseline_kind="internal",
        workflows=("spots", "families"),
        tolerances={"projection_atol": 1e-8},
    )
    validation = ValidationManifest(
        campaign_name="demo_validation",
        subsystem="structure_import",
        baseline_kind="iucr",
        status="foundational",
        reference_ids=("IUCr_A",),
        linked_benchmark_ids=("demo_benchmark",),
    )
    workflow = WorkflowResultManifest(
        result_id="demo_result",
        workflow_name="diffraction_geometry",
        modality="tem",
        produced_by="pytex",
        input_manifest_ids=("pytex.experiment_manifest",),
        artifact_paths=("docs/figures/diffraction_geometry.svg",),
    )
    benchmark_path = tmp_path / "benchmark.json"
    validation_path = tmp_path / "validation.json"
    workflow_path = tmp_path / "workflow.json"
    benchmark.write_json(benchmark_path)
    validation.write_json(validation_path)
    workflow.write_json(workflow_path)
    assert read_benchmark_manifest(benchmark_path).to_dict() == benchmark.to_dict()
    assert read_validation_manifest(validation_path).to_dict() == validation.to_dict()
    assert read_workflow_result_manifest(workflow_path).to_dict() == workflow.to_dict()


def test_manifest_schema_files_match_stable_identifiers() -> None:
    schema_paths = {
        experiment_manifest_schema_path(): "pytex.experiment_manifest",
        benchmark_manifest_schema_path(): "pytex.benchmark_manifest",
        validation_manifest_schema_path(): "pytex.validation_manifest",
        workflow_result_manifest_schema_path(): "pytex.workflow_result_manifest",
    }
    for path, schema_id in schema_paths.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["$id"] == schema_id
        assert payload["properties"]["schema_version"]["const"] == "1.0.0"


def test_benchmark_manifest_rejects_negative_tolerance() -> None:
    with pytest.raises(ValueError):
        BenchmarkManifest(
            benchmark_id="bad",
            subsystem="diffraction",
            baseline_kind="internal",
            workflows=("spots",),
            tolerances={"projection_atol": -1.0},
        )


def test_workflow_result_manifest_requires_non_empty_artifact_paths() -> None:
    with pytest.raises(ValueError):
        WorkflowResultManifest(
            result_id="bad",
            workflow_name="demo",
            modality="ebsd",
            produced_by="pytex",
            input_manifest_ids=("input",),
            artifact_paths=("",),
        )


def test_repo_manifests_are_schema_valid_and_readable() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    for manifest_path in sorted((repo_root / "benchmarks").glob("**/*.json")):
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        schema_id = payload["schema_id"]
        if schema_id == "pytex.benchmark_manifest":
            validate_benchmark_manifest(payload)
            assert read_benchmark_manifest(manifest_path).to_dict() == payload
        elif schema_id == "pytex.validation_manifest":
            validate_validation_manifest(payload)
            assert read_validation_manifest(manifest_path).to_dict() == payload
        elif schema_id == "pytex.workflow_result_manifest":
            validate_workflow_result_manifest(payload)
            assert read_workflow_result_manifest(manifest_path).to_dict() == payload
        else:
            pytest.fail(f"Unexpected schema_id in repo manifest {manifest_path}: {schema_id}")
