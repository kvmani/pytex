from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pytex import (
    CrystalDirection,
    CrystalMap,
    EBSDImportManifest,
    FrameDomain,
    Handedness,
    IPFColorKey,
    KikuchiPyWorkflowResult,
    Lattice,
    NormalizedEBSDDataset,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    RotationSet,
    SymmetrySpec,
    TextureReport,
    direction_from_orix_miller,
    from_orix_orientation,
    from_orix_rotation,
    from_orix_symmetry,
    index_hough,
    manifest_schema_path,
    normalize_ebsd,
    normalize_kikuchipy_dataset,
    normalize_kikuchipy_payload,
    normalize_pyebsdindex_payload,
    normalize_pyebsdindex_result,
    plane_from_orix_miller,
    plot_ipf_map,
    plot_kam_map,
    read_ebsd_import_manifest,
    refine_orientations,
    to_orix_direction,
    to_orix_orientation,
    to_orix_plane,
    to_orix_rotation,
    to_orix_symmetry,
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


def make_phase() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal, specimen, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def make_crystal_map() -> tuple[CrystalMap, Phase]:
    crystal, specimen, phase = make_phase()
    orientations = OrientationSet.from_euler_angles(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [45.0, 0.0, 0.0],
                [0.0, 45.0, 0.0],
                [45.0, 45.0, 0.0],
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
            ]
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    return crystal_map, phase


def test_ipf_color_key_maps_sector_vertices_to_rgb_primaries() -> None:
    _, _, symmetry = make_foundation()
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


def test_normalize_kikuchipy_payload_preserves_multiphase_assignments() -> None:
    crystal, specimen, _ = make_foundation()
    lattice_a = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    lattice_b = Lattice(2.9, 2.9, 2.9, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase_a = Phase(
        "fcc_demo",
        lattice=lattice_a,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
    )
    phase_b = Phase(
        "bcc_demo",
        lattice=lattice_b,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
        crystal_frame=crystal,
    )
    payload = {
        "coordinates": np.array([[0.0, 0.0], [1.0, 0.0]]),
        "euler_angles_deg": np.array([[0.0, 0.0, 0.0], [10.0, 20.0, 30.0]]),
        "phases": (
            {"phase_id": "0", "name": phase_a.name, "point_group": phase_a.symmetry.point_group},
            {"phase_id": "1", "name": phase_b.name, "point_group": phase_b.symmetry.point_group},
        ),
        "phase_ids": np.array([0, 1], dtype=np.int64),
        "grid_shape": (1, 2),
        "step_sizes": (1.0, 1.0),
    }
    dataset = normalize_kikuchipy_payload(
        payload,
        crystal_frame=crystal,
        specimen_frame=specimen,
        map_frame=specimen,
        phases={0: phase_a, 1: phase_b},
    )
    assert dataset.crystal_map.is_multiphase
    assert dataset.crystal_map.phase_summary() == {"fcc_demo": 1, "bcc_demo": 1}
    assert {phase["name"] for phase in dataset.manifest.phases} == {"fcc_demo", "bcc_demo"}


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


def test_crystal_map_texture_methods_build_odf_pf_and_ipf() -> None:
    crystal_map, phase = make_crystal_map()
    odf = crystal_map.to_odf(weights=[4.0, 3.0, 2.0, 1.0])
    pole_figure = crystal_map.pole_figure(np.array([1, 0, 0]))
    inverse_pole_figure = crystal_map.inverse_pole_figure("z")
    report = crystal_map.texture_report(
        poles=[np.array([1, 0, 0])],
        sample_directions=("x", "z"),
        plot=True,
    )
    assert odf.orientations.phase == phase
    assert pole_figure.pole.phase == phase
    assert inverse_pole_figure.crystal_frame == phase.crystal_frame
    assert isinstance(report, TextureReport)
    assert len(report.pole_figures) == 1
    assert len(report.inverse_pole_figures) == 2
    assert report.odf_figure is not None
    for figure in (*report.pole_figure_figures, *report.inverse_pole_figure_figures):
        assert figure is not None


def test_plot_ipf_map_and_kam_map_return_matplotlib_figures() -> None:
    crystal_map, _ = make_crystal_map()
    ipf_figure = plot_ipf_map(crystal_map, direction="z")
    kam_figure = plot_kam_map(crystal_map)
    assert ipf_figure.axes[0].get_title() == "IPF Map (z)"
    assert kam_figure.axes[0].get_title() == "Kernel Average Misorientation"


def test_orix_rotation_and_orientation_adapters_round_trip() -> None:
    pytest.importorskip("orix")
    crystal_map, phase = make_crystal_map()
    orix_rotation = to_orix_rotation(crystal_map.orientations.as_rotation_set())
    recovered_rotation = from_orix_rotation(orix_rotation)
    assert isinstance(recovered_rotation, RotationSet)
    assert np.allclose(recovered_rotation.as_matrices(), crystal_map.orientations.as_matrices())
    orix_orientation = to_orix_orientation(crystal_map.orientations)
    recovered_orientation = from_orix_orientation(
        orix_orientation,
        crystal_frame=phase.crystal_frame,
        specimen_frame=crystal_map.orientations.specimen_frame,
        phase=phase,
    )
    assert isinstance(recovered_orientation, OrientationSet)
    assert np.allclose(recovered_orientation.as_matrices(), crystal_map.orientations.as_matrices())
    recovered_symmetry = from_orix_symmetry(
        to_orix_symmetry(phase.symmetry),
        reference_frame=phase.crystal_frame,
    )
    assert recovered_symmetry.proper_point_group == phase.symmetry.proper_point_group
    plane = crystal_map.pole_figure(np.array([1, 0, 0])).pole
    orix_plane = to_orix_plane(plane)
    recovered_plane = plane_from_orix_miller(orix_plane, phase=phase)
    direction = CrystalDirection([1.0, 0.0, 0.0], phase=phase)
    orix_direction = to_orix_direction(direction)
    recovered_direction = direction_from_orix_miller(orix_direction, phase=phase)
    assert np.array_equal(recovered_plane.miller.indices, plane.miller.indices)
    assert np.allclose(recovered_direction.coordinates, direction.coordinates)


def test_kikuchipy_workflow_wrappers_normalize_hough_and_refinement_results() -> None:
    crystal, specimen, phase = make_phase()

    class XMapStub:
        def __init__(self, quaternions: np.ndarray, *, source_file: str) -> None:
            self.rotations = quaternions
            self.x = np.array([0.0, 1.0], dtype=np.float64)
            self.y = np.array([0.0, 0.0], dtype=np.float64)
            self.shape = (1, 2)
            self.dx = 1.0
            self.dy = 1.0
            self.source_file = source_file

    indexed_xmap = XMapStub(
        np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.92387953, 0.0, 0.0, 0.38268343],
            ]
        ),
        source_file="indexed.h5",
    )
    refined_xmap = XMapStub(
        np.array(
            [
                [0.98078528, 0.0, 0.0, 0.19509032],
                [0.92387953, 0.0, 0.0, 0.38268343],
            ]
        ),
        source_file="refined.h5",
    )

    class SignalStub:
        def __init__(self) -> None:
            self.xmap = indexed_xmap

        def hough_indexing(
            self,
            phase_list: object,
            indexer: object,
            **kwargs: object,
        ) -> tuple[object, np.ndarray, np.ndarray]:
            del phase_list, indexer, kwargs
            return indexed_xmap, np.array([1.0]), np.array([2.0])

        def refine_orientation(self, xmap: object, **kwargs: object) -> object:
            del xmap, kwargs
            return refined_xmap

    signal = SignalStub()
    dataset = normalize_ebsd(
        signal,
        frames=(crystal, specimen, specimen),
        phase=phase,
    )
    assert dataset.manifest.phase_name == phase.name
    result = index_hough(
        signal,
        phase_list=object(),
        indexer=object(),
        frames=(crystal, specimen, specimen),
        phase=phase,
    )
    assert isinstance(result, KikuchiPyWorkflowResult)
    assert result.index_data is not None
    assert result.band_data is not None
    refined = refine_orientations(
        signal,
        result,
        frames=(crystal, specimen, specimen),
        phase=phase,
        detector=object(),
        master_pattern=object(),
        energy=20.0,
    )
    assert refined.dataset.manifest.metadata["workflow"] == "orientation_refinement"
    assert not np.allclose(
        refined.dataset.crystal_map.orientations.quaternions,
        result.dataset.crystal_map.orientations.quaternions,
    )
