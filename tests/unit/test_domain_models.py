from __future__ import annotations

import numpy as np
import pytest

from pytex.core import (
    AcquisitionGeometry,
    CalibrationRecord,
    CrystalPlane,
    FrameDomain,
    FrameTransform,
    Handedness,
    Lattice,
    MeasurementQuality,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReciprocalLatticeVector,
    ReferenceFrame,
    Rotation,
    ScatteringSetup,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.diffraction import (
    DetectedSpotCluster,
    DetectorAcceptanceMask,
    DiffractionGeometry,
    DiffractionPattern,
    FamilyIndexingReport,
    IndexingCandidate,
    KinematicSimulation,
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


def test_crystal_map_accepts_multimodal_acquisition_context() -> None:
    orientation_set = make_orientation_set()
    specimen = orientation_set.specimen_frame
    map_frame = ReferenceFrame(
        name="map",
        domain=FrameDomain.MAP,
        axes=("i", "j", "k"),
        handedness=Handedness.RIGHT,
    )
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="ebsd",
        map_frame=map_frame,
        specimen_to_map=FrameTransform(
            source=specimen,
            target=map_frame,
            rotation_matrix=np.eye(3),
        ),
        calibration_record=CalibrationRecord(source="stage-fit", status="calibrated"),
        measurement_quality=MeasurementQuality(confidence=0.95, valid_fraction=1.0),
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0]]),
        orientations=orientation_set,
        map_frame=map_frame,
        grid_shape=(1, 2),
        step_sizes=(1.0, 1.0),
        acquisition_geometry=acquisition,
        calibration_record=acquisition.calibration_record,
        measurement_quality=acquisition.measurement_quality,
    )
    assert crystal_map.acquisition_geometry == acquisition


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


def _two_grain_map(euler_deg: list[list[float]]) -> CrystalMap:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.from_bunge_euler(*angles),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            )
            for angles in euler_deg
        ]
    )
    return CrystalMap(
        coordinates=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            dtype=np.float64,
        ),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(2.0, 2.0),
    )


def test_grain_scalar_metrics_gos_gam_and_equivalent_diameter() -> None:
    # Two grains (top row / bottom row). Grain 0 has a 2 deg internal spread,
    # grain 1 is perfectly uniform.
    crystal_map = _two_grain_map(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [25.0, 0.0, 0.0],
            [25.0, 0.0, 0.0],
        ]
    )
    segmentation = crystal_map.segment_grains(
        max_misorientation_deg=5.0,
        symmetry_aware=False,
        connectivity=4,
    )
    assert len(segmentation.grains) == 2

    gos = segmentation.grain_orientation_spread_deg()
    # Grain 0 mean sits ~1 deg from each member; grain 1 has zero spread.
    assert gos[0] == pytest.approx(1.0, abs=0.05)
    assert gos[1] == pytest.approx(0.0, abs=1e-6)

    gam = segmentation.grain_average_misorientation_deg()
    # Grain 0 neighbors differ by 2 deg; grain 1 by 0.
    assert gam[0] == pytest.approx(2.0, abs=0.05)
    assert gam[1] == pytest.approx(0.0, abs=1e-6)

    # Each grain has two 2x2 um pixels -> area 8 um^2 -> d = 2 sqrt(8/pi).
    diameters = segmentation.grain_equivalent_diameters()
    expected = 2.0 * np.sqrt(8.0 / np.pi)
    assert diameters[0] == pytest.approx(expected)
    assert diameters[1] == pytest.approx(expected)

    gos_map = segmentation.gos_map_deg()
    assert gos_map.shape == (2, 2)
    assert np.allclose(gos_map[1], 0.0)

    means = segmentation.grain_mean_orientations()
    assert set(means) == {0, 1}
    # Grain 1's mean recovers its single orientation (25, 0, 0).
    assert means[1].rotation.to_bunge_euler(degrees=True)[0] == pytest.approx(25.0, abs=1e-3)


def test_grain_equivalent_diameters_requires_step_sizes() -> None:
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
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0, symmetry_aware=False)
    with pytest.raises(ValueError, match="step_sizes"):
        segmentation.grain_equivalent_diameters()


def _single_grain_map(coordinates: np.ndarray) -> CrystalMap:
    crystal, specimen, symmetry = make_foundation()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            )
            for _ in range(len(coordinates))
        ]
    )
    return CrystalMap(
        coordinates=np.asarray(coordinates, dtype=np.float64),
        orientations=orientations,
        map_frame=specimen,
        step_sizes=(1.0, 1.0),
    )


def test_grain_fitted_ellipse_captures_elongation_and_orientation() -> None:
    # A horizontal 1x4 strip forms one grain elongated along x.
    crystal_map = _single_grain_map(
        np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0, symmetry_aware=False)
    assert len(segmentation.grains) == 1
    ellipse = segmentation.grain_fitted_ellipse(segmentation.grains[0])
    assert ellipse.semi_axes[0] > ellipse.semi_axes[1]
    assert ellipse.semi_axes[1] == pytest.approx(0.0, abs=1e-9)
    assert ellipse.aspect_ratio == float("inf")
    assert ellipse.angle_deg == pytest.approx(0.0, abs=1e-6)
    assert segmentation.grain_bounding_boxes()[0] == (3.0, 0.0)


def test_grain_fitted_ellipse_is_isotropic_for_square_grain() -> None:
    crystal_map = _single_grain_map(
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0, symmetry_aware=False)
    ellipse = segmentation.grain_fitted_ellipse(segmentation.grains[0])
    assert ellipse.aspect_ratio == pytest.approx(1.0, abs=1e-9)


def test_grain_shape_orientation_tracks_vertical_elongation() -> None:
    # A vertical 4x1 strip is elongated along y, so the major axis is at 90 deg.
    crystal_map = _single_grain_map(
        np.array([[0.0, 0.0], [0.0, 1.0], [0.0, 2.0], [0.0, 3.0]])
    )
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0, symmetry_aware=False)
    angle = segmentation.grain_shape_orientations_deg()[0]
    assert angle == pytest.approx(90.0, abs=1e-6)


def test_remove_wild_spikes_corrects_isolated_pixel() -> None:
    crystal, specimen, symmetry = make_foundation()
    # 3x3 grid: uniform orientation except a single wild spike in the center.
    euler = [[0.0, 0.0, 0.0]] * 9
    euler[4] = [40.0, 0.0, 0.0]
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.from_bunge_euler(*angles),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
            )
            for angles in euler
        ]
    )
    coordinates = np.array(
        [[float(c), float(r)] for r in range(3) for c in range(3)],
        dtype=np.float64,
    )
    crystal_map = CrystalMap(
        coordinates=coordinates,
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(3, 3),
        step_sizes=(1.0, 1.0),
        properties={"iq": np.arange(9, dtype=np.float64)},
    )
    cleaned = crystal_map.remove_wild_spikes(threshold_deg=10.0, symmetry_aware=False)
    # the spike is pulled back to the uniform neighborhood orientation
    spike_after = cleaned.orientations[4].rotation.to_bunge_euler(degrees=True)
    assert spike_after[0] == pytest.approx(0.0, abs=1e-6)
    # geometry and property channels are preserved
    assert cleaned.grid_shape == (3, 3)
    assert np.array_equal(cleaned.get_property("iq"), np.arange(9))
    # a spike-free map is returned unchanged (same instance, no rebuild)
    assert cleaned.remove_wild_spikes(threshold_deg=10.0, symmetry_aware=False) is cleaned


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
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, 1.0],
                [0.0, 2.0],
                [1.0, 2.0],
                [2.0, 2.0],
            ],
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
    _, specimen, _ = make_foundation()
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


def test_diffraction_geometry_accepts_multimodal_context() -> None:
    _, specimen, _ = make_foundation()
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
    calibration = CalibrationRecord(source="pattern-fit", status="refined", residual_error=0.2)
    quality = MeasurementQuality(confidence=0.9, valid_fraction=0.98, masked_fraction=0.02)
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="tem",
        detector_frame=detector,
        laboratory_frame=lab,
        specimen_to_detector=FrameTransform(
            source=specimen,
            target=detector,
            rotation_matrix=np.eye(3),
        ),
        specimen_to_laboratory=FrameTransform(
            source=specimen,
            target=lab,
            rotation_matrix=np.eye(3),
        ),
        calibration_record=calibration,
        measurement_quality=quality,
    )
    scattering = ScatteringSetup(
        laboratory_frame=lab,
        beam_energy_kev=20.0,
        incident_beam_direction=np.array([0.0, 0.0, 1.0]),
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
        acquisition_geometry=acquisition,
        calibration_record=calibration,
        measurement_quality=quality,
        scattering_setup=scattering,
    )
    assert geometry.acquisition_geometry == acquisition
    assert geometry.scattering_setup == scattering


def test_diffraction_geometry_rejects_inconsistent_acquisition_context() -> None:
    _, specimen, _ = make_foundation()
    detector = ReferenceFrame(
        name="detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
    )
    other_detector = ReferenceFrame(
        name="other_detector",
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
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="tem",
        detector_frame=other_detector,
        laboratory_frame=lab,
        specimen_to_detector=FrameTransform(
            source=specimen,
            target=other_detector,
            rotation_matrix=np.eye(3),
        ),
        specimen_to_laboratory=FrameTransform(
            source=specimen,
            target=lab,
            rotation_matrix=np.eye(3),
        ),
    )
    with pytest.raises(ValueError):
        DiffractionGeometry(
            detector_frame=detector,
            specimen_frame=specimen,
            laboratory_frame=lab,
            beam_energy_kev=20.0,
            camera_length_mm=100.0,
            pattern_center=np.array([0.5, 0.5, 0.7]),
            detector_pixel_size_um=(60.0, 60.0),
            detector_shape=(480, 480),
            acquisition_geometry=acquisition,
        )


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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
    _, specimen, _ = make_foundation()
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
        coordinates_px=np.vstack(
            [spot.detector_coordinates_px for spot in simulation.accepted_spots()]
        ),
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
        coordinates_px=np.vstack(
            [spot.detector_coordinates_px for spot in simulation.accepted_spots()]
        ),
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
        coordinates_px=np.vstack(
            [spot.detector_coordinates_px for spot in pattern_simulation.accepted_spots()]
        ),
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
    reduced_axis = symmetry.reduce_vector_to_fundamental_sector(
        projected.rotation.axis, antipodal=True
    )
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
