from __future__ import annotations

import numpy as np
import pytest

from pytex.core import (
    CrystalPlane,
    ReciprocalLatticeVector,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.diffraction import DiffractionGeometry, DiffractionPattern, KinematicSimulation
from pytex.diffraction import (
    DetectorAcceptanceMask,
    DetectedSpotCluster,
    FamilyIndexingReport,
    IndexingCandidate,
    OrientationIndexingCandidate,
    OrientationRefinementResult,
    ReflectionFamily,
    SpotAssignment,
)
from pytex.ebsd import CrystalMap, GrainGraph
from pytex.texture import ODF


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
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    return crystal, specimen, symmetry


def make_orientation_set() -> OrientationSet:
    crystal, specimen, symmetry = make_foundation()
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ),
        Orientation(
            Rotation.from_bunge_euler(15.0, 20.0, 25.0),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ),
    ]
    return OrientationSet.from_orientations(orientations)


def test_provenance_metadata_is_immutable() -> None:
    provenance = ProvenanceRecord(source_system="test", metadata={"vendor": "demo"})
    with pytest.raises(TypeError):
        provenance.metadata["vendor"] = "other"  # type: ignore[index]


def test_crystal_map_validates_lengths() -> None:
    orientation_set = make_orientation_set()
    specimen = orientation_set.specimen_frame
    with pytest.raises(ValueError):
        CrystalMap(
            coordinates=np.array([[0.0, 0.0]]),
            orientations=orientation_set,
            map_frame=specimen,
        )


def test_crystal_map_validates_positive_step_sizes() -> None:
    orientation_set = make_orientation_set()
    specimen = orientation_set.specimen_frame
    with pytest.raises(ValueError):
        CrystalMap(
            coordinates=np.array([[0.0, 0.0], [1.0, 0.0]]),
            orientations=orientation_set,
            map_frame=specimen,
            step_sizes=(1.0, 0.0),
        )


def test_crystal_map_neighbor_pairs_cover_regular_grid_edges() -> None:
    orientation_set = make_orientation_set()
    specimen = orientation_set.specimen_frame
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0]]),
        orientations=orientation_set,
        map_frame=specimen,
        grid_shape=(1, 2),
        step_sizes=(1.0, 1.0),
    )
    assert np.array_equal(crystal_map.neighbor_pairs(), np.array([[0, 1]]))


def test_crystal_map_kernel_average_misorientation_returns_grid() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(5.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(10.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(15.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    kam = crystal_map.kernel_average_misorientation_deg(symmetry_aware=False)
    assert kam.shape == (2, 2)
    assert np.all(kam >= 0.0)
    assert np.all(np.isfinite(kam))


def test_crystal_map_segment_grains_splits_components_by_threshold() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(1.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(25.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(26.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    assert segmentation.label_grid.shape == (2, 2)
    assert len(segmentation.grains) == 2
    assert np.array_equal(segmentation.label_grid[0], np.array([0, 0]))
    assert np.array_equal(segmentation.label_grid[1], np.array([1, 1]))


def test_grain_segmentation_grod_map_is_zero_at_reference_orientations() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(2.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(25.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(27.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    grod = segmentation.grod_map_deg()
    assert grod.shape == (2, 2)
    assert np.all(grod >= 0.0)
    assert np.all(np.isfinite(grod))
    for grain in segmentation.grains:
        row, col = divmod(grain.reference_orientation_index, 2)
        assert grod[row, col] == pytest.approx(0.0, abs=1e-12)


def test_grain_boundary_network_extracts_cross_grain_edges() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(1.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(25.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(26.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    network = segmentation.boundary_network()
    assert network.count == 2
    assert network.mean_misorientation_deg > 0.0
    for segment in network.segments:
        assert segment.left_grain_id != segment.right_grain_id
        assert segment.misorientation_deg >= 0.0


def test_grain_graph_aggregates_boundary_connectivity() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(1.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(25.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(26.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    graph = segmentation.grain_graph()
    assert isinstance(graph, GrainGraph)
    assert graph.edge_count == 1
    assert np.array_equal(graph.adjacency_matrix, np.array([[0, 1], [1, 0]]))
    assert np.array_equal(graph.neighbors(0), np.array([1]))
    assert graph.edges[0].total_length == pytest.approx(2.0)


def test_merge_small_grains_absorbs_isolated_single_pixel_grain() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
            Orientation(
                Rotation.from_bunge_euler(25.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            ),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    merged = segmentation.merge_small_grains(min_size=2)
    assert len(segmentation.grains) == 2
    assert len(merged.grains) == 1
    assert np.all(merged.label_grid == 0)


def test_grain_segmentation_majority_smoothed_removes_isolated_label_noise() -> None:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.from_bunge_euler(25.0, 0.0, 0.0), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
            Orientation(Rotation.identity(), crystal_frame=crystal, specimen_frame=specimen, symmetry=symmetry),
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 1.0], [0.0, 2.0], [1.0, 2.0], [2.0, 2.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(3, 3),
        step_sizes=(1.0, 1.0),
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    smoothed = segmentation.majority_smoothed(iterations=1, min_neighbor_votes=3)
    assert len(smoothed.grains) == 1
    assert np.all(smoothed.label_grid == 0)


def test_odf_requires_non_negative_weights() -> None:
    orientation_set = make_orientation_set()
    with pytest.raises(ValueError):
        ODF(orientations=orientation_set, weights=np.array([0.6, -0.1]))


def test_diffraction_geometry_produces_positive_wavelength() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(60.0, 60.0),
        detector_shape=(480, 480),
    )
    assert geometry.electron_wavelength_angstrom > 0.0


def test_crystal_plane_exposes_expected_cubic_d_spacing() -> None:
    crystal, _, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    plane = CrystalPlane(miller=MillerIndex(np.array([2, 0, 0]), phase=phase), phase=phase)
    assert np.isclose(plane.d_spacing_angstrom, 1.5)


def test_reciprocal_lattice_vector_matches_plane_spacing_magnitude() -> None:
    crystal, _, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    reciprocal_vector = ReciprocalLatticeVector.from_miller_index(
        MillerIndex(np.array([2, 0, 0]), phase=phase)
    )
    assert np.isclose(reciprocal_vector.magnitude_inv_angstrom, 2.0 / 3.0)


def test_zone_axis_exposes_unit_direction() -> None:
    crystal, _, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    zone_axis = ZoneAxis(indices=np.array([1, 1, 0]), phase=phase)
    assert np.allclose(zone_axis.unit_vector, np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0))


def test_zone_axis_zone_law_helpers_match_expected_plane_membership() -> None:
    crystal, _, symmetry = make_foundation()
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    zone_axis = ZoneAxis(indices=np.array([0, 0, 1]), phase=phase)
    in_zone = MillerIndex(np.array([1, 0, 0]), phase=phase)
    out_of_zone = MillerIndex(np.array([0, 0, 1]), phase=phase)
    assert zone_axis.zone_law_value(in_zone) == 0
    assert zone_axis.contains_miller_index(in_zone)
    assert zone_axis.zone_law_value(out_of_zone) == 1
    assert not zone_axis.contains_miller_index(out_of_zone)


def test_diffraction_geometry_rejects_invalid_detector_domain() -> None:
    crystal, specimen, _ = make_foundation()
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    with pytest.raises(ValueError):
        DiffractionGeometry(
            detector_frame=crystal,
            specimen_frame=specimen,
            laboratory_frame=lab,
            beam_energy_kev=20.0,
            camera_length_mm=100.0,
            pattern_center=np.array([0.5, 0.5, 0.7]),
            detector_pixel_size_um=(60.0, 60.0),
            detector_shape=(480, 480),
        )


def test_diffraction_pattern_validates_lengths() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(60.0, 60.0),
        detector_shape=(480, 480),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    with pytest.raises(ValueError):
        DiffractionPattern(
            coordinates_px=np.array([[1.0, 2.0], [3.0, 4.0]]),
            intensities=np.array([1.0]),
            geometry=geometry,
            phase=phase,
        )


def test_diffraction_geometry_center_pixel_has_zero_scattering_angle() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
    )
    center_pixel = np.array([[50.0, 50.0]])
    assert np.allclose(geometry.detector_coordinates_mm(center_pixel), np.array([[0.0, 0.0]]))
    assert np.allclose(geometry.two_theta_rad(center_pixel), np.array([0.0]), atol=1e-10)


def test_diffraction_geometry_accepts_specimen_to_lab_rotation_matrix() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    ninety_deg_about_z = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
        specimen_to_lab_matrix=ninety_deg_about_z,
    )
    mapped = geometry.specimen_vectors_to_lab(np.array([[1.0, 0.0, 0.0]]))
    assert np.allclose(mapped, np.array([[0.0, 1.0, 0.0]]))


def test_diffraction_geometry_rejects_non_positive_pattern_center_depth() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    with pytest.raises(ValueError):
        DiffractionGeometry(
            detector_frame=detector,
            specimen_frame=specimen,
            laboratory_frame=lab,
            beam_energy_kev=20.0,
            camera_length_mm=100.0,
            pattern_center=np.array([0.5, 0.5, 0.0]),
            detector_pixel_size_um=(100.0, 100.0),
            detector_shape=(101, 101),
        )


def test_diffraction_geometry_rejects_pattern_center_outside_detector_extent() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    with pytest.raises(ValueError):
        DiffractionGeometry(
            detector_frame=detector,
            specimen_frame=specimen,
            laboratory_frame=lab,
            beam_energy_kev=20.0,
            camera_length_mm=100.0,
            pattern_center=np.array([1.1, 0.5, 0.7]),
            detector_pixel_size_um=(100.0, 100.0),
            detector_shape=(101, 101),
        )


def test_diffraction_geometry_reports_expected_azimuth_and_positive_two_theta() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
    )
    right_pixel = np.array([[60.0, 50.0]])
    up_pixel = np.array([[50.0, 40.0]])
    assert geometry.two_theta_rad(right_pixel)[0] > 0.0
    assert np.isclose(geometry.azimuth_rad(right_pixel)[0], 0.0, atol=1e-8)
    assert np.isclose(abs(geometry.azimuth_rad(up_pixel)[0]), np.pi / 2.0, atol=1e-8)


def test_diffraction_geometry_rejects_projection_behind_detector_plane() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
    )
    coordinates_px, valid = geometry.project_directions_to_detector_px(np.array([[0.0, 0.0, -1.0]]))
    assert not valid[0]
    assert np.isnan(coordinates_px[0, 0])
    assert np.isnan(coordinates_px[0, 1])


def test_diffraction_geometry_predicts_bragg_ring_radius_from_plane_spacing() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(512, 512),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    plane = CrystalPlane(miller=MillerIndex(np.array([2, 0, 0]), phase=phase), phase=phase)
    two_theta = geometry.bragg_two_theta_rad(plane.d_spacing_angstrom)
    assert np.isclose(
        geometry.ring_radius_mm_for_plane(plane),
        geometry.camera_length_mm * np.tan(two_theta),
    )
    assert geometry.ring_radius_mm_for_plane(plane) > 0.0


def test_diffraction_geometry_clips_bragg_argument_at_numerical_boundary() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
    )
    boundary_spacing = geometry.electron_wavelength_angstrom / (2.0 * (1.0 + 5e-13))
    assert np.isclose(geometry.bragg_two_theta_rad(boundary_spacing), np.pi)


def test_diffraction_pattern_exposes_geometry_derived_arrays() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        camera_length_mm=100.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(100.0, 100.0),
        detector_shape=(101, 101),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    pattern = DiffractionPattern(
        coordinates_px=np.array([[50.0, 50.0], [60.0, 50.0]]),
        intensities=np.array([1.0, 2.0]),
        geometry=geometry,
        phase=phase,
    )
    assert pattern.detector_coordinates_mm().shape == (2, 2)
    assert pattern.outgoing_directions_lab().shape == (2, 3)
    assert pattern.scattering_vectors_lab().shape == (2, 3)
    assert pattern.two_theta_rad().shape == (2,)
    assert pattern.azimuth_rad().shape == (2,)


def test_kinematic_simulation_projects_zone_axis_spots() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 1]]),
        orientation=orientation,
        zone_axis=ZoneAxis(indices=np.array([0, 0, 1]), phase=phase),
        max_excitation_error_inv_angstrom=0.2,
    )
    assert len(simulation.spots) >= 2
    for spot in simulation.spots:
        assert spot.on_detector
        assert spot.two_theta_rad >= 0.0


def test_kinematic_simulation_rejects_mismatched_orientation_frame() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    other_specimen = ReferenceFrame(
        name="other_specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=other_specimen,
        symmetry=symmetry,
        phase=phase,
    )
    with pytest.raises(ValueError):
        KinematicSimulation.simulate_spots(
            geometry,
            phase,
            np.array([[1, 0, 0]]),
            orientation=orientation,
        )


def test_kinematic_simulation_rejects_non_integer_miller_indices() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    with pytest.raises(ValueError):
        KinematicSimulation.simulate_spots(
            geometry,
            phase,
            np.array([[1.0, 0.5, 0.0]]),
        )


def test_kinematic_simulation_rejects_negative_excitation_threshold() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    with pytest.raises(ValueError):
        KinematicSimulation.simulate_spots(
            geometry,
            phase,
            np.array([[1, 0, 0]]),
            max_excitation_error_inv_angstrom=-0.1,
        )


def test_kinematic_simulation_marks_off_detector_spots_without_invalid_angles() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(32, 32),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [-1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.2,
    )
    assert simulation.spots
    assert any(not spot.on_detector for spot in simulation.spots)
    for spot in simulation.spots:
        assert np.isfinite(spot.two_theta_rad)
        assert np.isfinite(spot.azimuth_rad)


def test_detector_acceptance_mask_filters_by_radius_and_inset() -> None:
    crystal, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    mask = DetectorAcceptanceMask(inset_px=(20.0, 20.0), max_radius_px=40.0)
    accepted = mask.contains(
        geometry,
        np.array(
            [
                geometry.pattern_center_px,
                geometry.pattern_center_px + np.array([50.0, 0.0]),
                np.array([5.0, 5.0]),
            ]
        ),
    )
    assert np.array_equal(accepted, np.array([True, False, False]))


def test_kinematic_simulation_groups_symmetry_equivalent_reflections_into_families() -> None:
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
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    symmetry = SymmetrySpec.from_point_group("432", reference_frame=crystal)
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.2,
    )
    assert len(simulation.spots) == 3
    assert len(simulation.reflection_families) == 1
    family = simulation.reflection_families[0]
    assert isinstance(family, ReflectionFamily)
    assert family.multiplicity == 3
    assert np.array_equal(family.spot_indices, np.array([0, 1, 2]))
    assert all(spot.family_id == 0 for spot in simulation.spots)


def test_kinematic_simulation_can_deduplicate_reflection_families() -> None:
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
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    symmetry = SymmetrySpec.from_point_group("432", reference_frame=crystal)
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.2,
        deduplicate_families=True,
    )
    assert len(simulation.spots) == 1
    assert len(simulation.reflection_families) == 1
    assert simulation.spots[0].family_id == 0


def test_kinematic_simulation_applies_acceptance_mask_without_dropping_spots() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[0, 0, 1], [1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.4,
        acceptance_mask=DetectorAcceptanceMask(max_radius_px=10.0),
    )
    assert len(simulation.spots) == 2
    assert simulation.spots[0].accepted_by_mask
    assert not simulation.spots[1].accepted_by_mask
    assert len(simulation.accepted_spots()) == 1


def test_kinematic_simulation_proxy_intensity_penalizes_higher_order_reflections() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [2, 0, 0]]),
        max_excitation_error_inv_angstrom=0.4,
        intensity_model="kinematic_proxy",
        excitation_sigma_inv_angstrom=0.2,
    )
    intensities = {tuple(spot.miller_indices.tolist()): spot.intensity for spot in simulation.spots}
    assert intensities[(1, 0, 0)] > intensities[(2, 0, 0)]


def test_diffraction_pattern_clusters_nearby_observations() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    pattern = DiffractionPattern(
        coordinates_px=np.array([[100.0, 100.0], [103.0, 100.0], [200.0, 200.0]]),
        intensities=np.array([2.0, 1.0, 5.0]),
        geometry=geometry,
        phase=phase,
    )
    clusters = pattern.cluster_observations(max_distance_px=5.0)
    assert len(clusters) == 2
    assert isinstance(clusters[0], DetectedSpotCluster)
    assert clusters[0].member_indices.shape[0] == 2


def test_kinematic_simulation_can_associate_to_pattern() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [-1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.2,
    )
    coordinates = np.vstack([spot.detector_coordinates_px for spot in simulation.accepted_spots()])
    pattern = DiffractionPattern(
        coordinates_px=coordinates + np.array([[1.0, 0.0], [-1.0, 0.5]]),
        intensities=np.array([10.0, 8.0]),
        geometry=geometry,
        phase=phase,
    )
    indexing = simulation.associate_to_pattern(pattern, max_distance_px=5.0, cluster_radius_px=2.0)
    assert isinstance(indexing, IndexingCandidate)
    assert len(indexing.matches) == 2
    assert isinstance(indexing.matches[0], SpotAssignment)
    assert indexing.match_fraction == pytest.approx(1.0)
    assert indexing.mean_residual_px < 2.0
    assert indexing.score > 0.0


def test_kinematic_simulation_can_rank_orientation_candidates() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    true_orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    wrong_orientation = Orientation(
        rotation=Rotation.from_bunge_euler(0.0, 30.0, 0.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
        orientation=true_orientation,
        max_excitation_error_inv_angstrom=0.2,
    )
    pattern = DiffractionPattern(
        coordinates_px=np.vstack([spot.detector_coordinates_px for spot in simulation.accepted_spots()]),
        intensities=np.array([10.0, 9.0, 8.0, 7.0]),
        geometry=geometry,
        phase=phase,
    )
    ranked = KinematicSimulation.rank_orientation_candidates(
        geometry,
        phase,
        pattern,
        np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
        [wrong_orientation, true_orientation],
        max_excitation_error_inv_angstrom=0.2,
        max_distance_px=5.0,
        cluster_radius_px=2.0,
    )
    assert isinstance(ranked[0], OrientationIndexingCandidate)
    assert ranked[0].orientation_index == 1
    assert ranked[0].score >= ranked[1].score


def test_indexing_candidate_reports_family_level_aggregation() -> None:
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
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    symmetry = SymmetrySpec.from_point_group("432", reference_frame=crystal)
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0]]),
        max_excitation_error_inv_angstrom=0.2,
    )
    pattern = DiffractionPattern(
        coordinates_px=np.vstack([spot.detector_coordinates_px for spot in simulation.accepted_spots()]),
        intensities=np.array([10.0, 9.0, 8.0]),
        geometry=geometry,
        phase=phase,
    )
    indexing = simulation.associate_to_pattern(pattern, max_distance_px=5.0, cluster_radius_px=2.0)
    reports = indexing.family_reports()
    assert len(reports) == 1
    assert isinstance(reports[0], FamilyIndexingReport)
    assert reports[0].family_id == 0
    assert reports[0].multiplicity == 3
    assert reports[0].simulated_spot_count == 3
    assert reports[0].matched_spot_count == 3
    assert reports[0].matched_fraction == pytest.approx(1.0)


def test_refine_orientation_candidate_improves_or_preserves_candidate_score() -> None:
    crystal, specimen, symmetry = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    lab = ReferenceFrame(
        name="lab", domain=FrameDomain.LABORATORY, axes=("X", "Y", "Z"), handedness=Handedness.RIGHT
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 0.7]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
    )
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="demo",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
    )
    true_orientation = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    seed_orientation = Orientation(
        rotation=Rotation.from_bunge_euler(0.0, 6.0, 0.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    miller = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])
    pattern_simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        miller,
        orientation=true_orientation,
        max_excitation_error_inv_angstrom=0.2,
    )
    pattern = DiffractionPattern(
        coordinates_px=np.vstack([spot.detector_coordinates_px for spot in pattern_simulation.accepted_spots()]),
        intensities=np.array([10.0, 9.0, 8.0, 7.0]),
        geometry=geometry,
        phase=phase,
    )
    seed_simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        miller,
        orientation=seed_orientation,
        max_excitation_error_inv_angstrom=0.2,
    )
    seed_indexing = seed_simulation.associate_to_pattern(
        pattern,
        max_distance_px=5.0,
        cluster_radius_px=2.0,
    )
    refinement = KinematicSimulation.refine_orientation_candidate(
        geometry,
        phase,
        pattern,
        miller,
        seed_orientation,
        max_excitation_error_inv_angstrom=0.2,
        max_distance_px=5.0,
        cluster_radius_px=2.0,
        search_half_width_deg=6.0,
        step_deg=3.0,
        iterations=2,
    )
    assert isinstance(refinement, OrientationRefinementResult)
    assert refinement.evaluated_candidates > 1
    assert refinement.refined_candidate.score >= seed_indexing.score
    assert refinement.refined_candidate.indexing.mean_residual_px <= seed_indexing.mean_residual_px


def test_fundamental_region_key_matches_for_symmetry_equivalent_orientations() -> None:
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
    base = Orientation(
        rotation=Rotation.from_bunge_euler(40.0, 30.0, 20.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    equivalent = Orientation(
        rotation=Rotation.from_matrix(base.as_matrix() @ symmetry.operators[1]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    assert np.allclose(base.fundamental_region_key(), equivalent.fundamental_region_key())
    projected = base.project_to_fundamental_region()
    reduced_axis = symmetry.reduce_vector_to_fundamental_sector(projected.rotation.axis, antipodal=True)
    assert symmetry.vector_in_fundamental_sector(reduced_axis, antipodal=True)


def test_exact_fundamental_region_projection_selects_minimum_angle_representative() -> None:
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
    orientation = Orientation(
        rotation=Rotation.from_bunge_euler(80.0, 35.0, 25.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
    )
    projected = orientation.project_to_exact_fundamental_region()
    equivalent_angles = [
        Orientation(
            rotation=Rotation(quaternion),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=symmetry,
        ).rotation.angle_rad
        for quaternion in orientation.equivalent_orientations().quaternions
    ]
    assert projected.rotation.angle_rad == pytest.approx(min(equivalent_angles))
    assert projected.is_in_fundamental_region()
