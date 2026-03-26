from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    EBSDImportManifest,
    FrameDomain,
    Handedness,
    IPFColorKey,
    Lattice,
    NormalizedEBSDDataset,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    manifest_schema_path,
    normalize_kikuchipy_dataset,
    normalize_kikuchipy_payload,
    normalize_pyebsdindex_result,
    normalize_pyebsdindex_payload,
    read_ebsd_import_manifest,
    validate_ebsd_import_manifest,
)


def make_foundation() -> tuple[ReferenceFrame, ReferenceFrame, SymmetrySpec]:
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
    symmetry = SymmetrySpec.from_point_group("432", reference_frame=crystal)
    return crystal, specimen, symmetry


def test_ipf_color_key_maps_sector_vertices_to_rgb_primaries() -> None:
    crystal, _, symmetry = make_foundation()
    key = IPFColorKey(
        crystal_symmetry=symmetry,
        specimen_direction=np.array([0.0, 0.0, 1.0]),
    )
    colors = key.colors_from_crystal_directions(key.sector_vertices)
    assert np.allclose(colors, np.eye(3), atol=1e-8)


def test_ipf_color_key_generates_finite_colors_from_orientations() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(35.0, 25.0, 10.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction=np.array([0.0, 0.0, 1.0]))
    colors = key.colors_from_orientations(orientations)
    assert colors.shape == (2, 3)
    assert np.all(np.isfinite(colors))
    assert np.all((colors >= 0.0) & (colors <= 1.0))


def test_ebsd_import_manifest_round_trip_and_validation() -> None:
    manifest = EBSDImportManifest(
        source_system="kikuchipy",
        source_file="demo.h5",
        phase_name="fcc_demo",
        point_group="m-3m",
        crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
        specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
        metadata={"vendor": "demo"},
    )
    payload = manifest.to_dict()
    validate_ebsd_import_manifest(payload)
    restored = EBSDImportManifest.from_dict(payload)
    assert restored.to_dict() == payload


def test_ebsd_import_manifest_requires_schema_fields() -> None:
    manifest = EBSDImportManifest(
        source_system="pyebsdindex",
        source_file="demo.h5",
        phase_name="fcc_demo",
        point_group="m-3m",
        crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
        specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
    ).to_dict()
    manifest.pop("phase_name")
    with pytest.raises(ValueError):
        validate_ebsd_import_manifest(manifest)


def test_normalized_ebsd_dataset_requires_phase_alignment() -> None:
    crystal, specimen, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
                phase=phase,
            )
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0]]),
        orientations=orientations,
        map_frame=specimen,
    )
    manifest = EBSDImportManifest(
        source_system="kikuchipy",
        source_file="demo.h5",
        phase_name="fcc_demo",
        point_group="m-3m",
        crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
        specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
    )
    dataset = NormalizedEBSDDataset(crystal_map=crystal_map, manifest=manifest)
    assert dataset.manifest.phase_name == "fcc_demo"


def test_manifest_schema_file_matches_stable_identifiers() -> None:
    schema_path = manifest_schema_path()
    assert schema_path == Path("schemas/ebsd_import_manifest.schema.json").resolve()
    payload = json.loads(schema_path.read_text(encoding="utf-8"))
    assert payload["$id"] == "pytex.ebsd_import_manifest"
    assert payload["properties"]["schema_version"]["const"] == "1.0.0"


def test_ebsd_import_manifest_can_round_trip_through_json_file(tmp_path: Path) -> None:
    manifest = EBSDImportManifest(
        source_system="kikuchipy",
        source_file="demo.h5",
        phase_name="fcc_demo",
        point_group="m-3m",
        crystal_frame={"name": "crystal", "axis_1": "a", "axis_2": "b", "axis_3": "c"},
        specimen_frame={"name": "specimen", "axis_1": "x", "axis_2": "y", "axis_3": "z"},
    )
    manifest_path = tmp_path / "manifest.json"
    manifest.write_json(manifest_path)
    restored = read_ebsd_import_manifest(manifest_path)
    assert restored.to_dict() == manifest.to_dict()


def test_normalize_kikuchipy_payload_builds_normalized_dataset() -> None:
    crystal, specimen, symmetry = make_foundation()
    payload = {
        "coordinates": np.array([[0.0, 0.0], [1.0, 0.0]]),
        "euler_angles_deg": np.array([[0.0, 0.0, 0.0], [10.0, 20.0, 30.0]]),
        "point_group": symmetry.point_group,
        "phase_name": "fcc_demo",
        "grid_shape": (1, 2),
        "step_sizes": (1.0, 1.0),
    }
    dataset = normalize_kikuchipy_payload(
        payload,
        crystal_frame=crystal,
        specimen_frame=specimen,
        map_frame=specimen,
    )
    assert isinstance(dataset, NormalizedEBSDDataset)
    assert dataset.manifest.source_system == "kikuchipy"
    assert len(dataset.crystal_map.orientations) == 2


def test_normalize_pyebsdindex_payload_accepts_vendor_angle_key() -> None:
    crystal, specimen, symmetry = make_foundation()
    payload = {
        "coordinates": np.array([[0.0, 0.0]]),
        "phi1_Phi_phi2_deg": np.array([[5.0, 10.0, 15.0]]),
        "point_group": symmetry.point_group,
        "phase_name": "fcc_demo",
    }
    dataset = normalize_pyebsdindex_payload(
        payload,
        crystal_frame=crystal,
        specimen_frame=specimen,
        map_frame=specimen,
    )
    assert dataset.manifest.source_system == "pyebsdindex"
    assert len(dataset.crystal_map.orientations) == 1


def test_normalize_kikuchipy_dataset_extracts_from_object_bridge() -> None:
    crystal, specimen, symmetry = make_foundation()

    class Rotations:
        def __init__(self) -> None:
            self.euler_angles_deg = np.array([[0.0, 0.0, 0.0], [5.0, 10.0, 15.0]])

    class PhaseStub:
        def __init__(self) -> None:
            self.name = "fcc_demo"
            self.point_group = symmetry.point_group

    class XMapStub:
        def __init__(self) -> None:
            self.coordinates = np.array([[0.0, 0.0], [1.0, 0.0]])
            self.rotations = Rotations()
            self.phase = PhaseStub()
            self.grid_shape = (1, 2)
            self.step_sizes = (1.0, 1.0)
            self.source_file = "demo.h5"
            self.metadata = {"vendor": "stub"}

    class KikuchiPyDatasetStub:
        def __init__(self) -> None:
            self.xmap = XMapStub()

    dataset = normalize_kikuchipy_dataset(
        KikuchiPyDatasetStub(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        map_frame=specimen,
    )
    assert dataset.manifest.source_system == "kikuchipy"
    assert dataset.manifest.source_file == "demo.h5"
    assert len(dataset.crystal_map.orientations) == 2


def test_normalize_pyebsdindex_result_extracts_from_object_bridge() -> None:
    crystal, specimen, symmetry = make_foundation()

    class PhaseStub:
        def __init__(self) -> None:
            self.name = "fcc_demo"
            self.point_group = symmetry.point_group

    class ScanStub:
        def __init__(self) -> None:
            self.coordinates = np.array([[0.0, 0.0]])
            self.phi1_Phi_phi2_deg = np.array([[10.0, 20.0, 30.0]])
            self.phase = PhaseStub()
            self.source_file = "scan.ang"

    class PyEBSDIndexResultStub:
        def __init__(self) -> None:
            self.scan = ScanStub()

    dataset = normalize_pyebsdindex_result(
        PyEBSDIndexResultStub(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        map_frame=specimen,
    )
    assert dataset.manifest.source_system == "pyebsdindex"
    assert dataset.manifest.source_file == "scan.ang"
    assert len(dataset.crystal_map.orientations) == 1
