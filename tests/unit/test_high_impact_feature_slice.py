from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pytex import (
    AcquisitionGeometry,
    CalibrationRecord,
    DiffractionGeometry,
    DiffractionPattern,
    EBSDCalibrationGeometry,
    EBSDDetectorGeometry,
    FrameDomain,
    FrameTransform,
    Handedness,
    Lattice,
    MeasurementQuality,
    MillerBravaisDirection,
    MillerBravaisPlane,
    PatternCenter,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    ScatteringSetup,
    SymmetrySpec,
    ZoneAxis,
    direct_to_reciprocal_components,
    estimate_zone_axis,
    from_json_contract,
    index_saed_pattern,
    metric_tensor,
    reciprocal_metric_tensor,
    reciprocal_to_direct_components,
    to_json_contract,
)
from pytex.diffraction import KinematicSimulation


def _frame(name: str, domain: FrameDomain, axes: tuple[str, str, str]) -> ReferenceFrame:
    return ReferenceFrame(name=name, domain=domain, axes=axes, handedness=Handedness.RIGHT)


def _phase(
    *,
    a: float = 3.0,
    b: float = 3.0,
    c: float = 3.0,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
    point_group: str = "m-3m",
) -> Phase:
    crystal = _frame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
    lattice = Lattice(a, b, c, alpha, beta, gamma, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=crystal)
    return Phase("demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def test_miller_bravais_plane_direction_round_trip_and_zone_law() -> None:
    phase = _phase(a=2.95, b=2.95, c=4.68, gamma=120.0, point_group="6/mmm")

    plane = MillerBravaisPlane.from_hkil([1, 0, -1, 0], phase=phase)
    direction = MillerBravaisDirection.from_UVTW([0, 0, 0, 1], phase=phase)
    equivalent_direction = MillerBravaisDirection.from_uvw([1, 0, 0], phase=phase)

    assert_array_equal(plane.hkl, np.array([1, 0, 0], dtype=np.int64))
    assert_array_equal(direction.uvw, np.array([0, 0, 1], dtype=np.int64))
    assert_array_equal(equivalent_direction.indices, np.array([2, -1, -1, 0], dtype=np.int64))
    assert plane.contains_direction(direction)
    assert not plane.contains_direction(equivalent_direction)
    assert_array_equal(plane.to_miller_plane().hkil, plane.indices)
    assert_array_equal(direction.to_miller_direction().UVTW, direction.indices)


def test_miller_bravais_rejects_invalid_four_index_constraints() -> None:
    phase = _phase(a=2.95, b=2.95, c=4.68, gamma=120.0, point_group="6/mmm")

    with pytest.raises(ValueError, match="i = -\\(h \\+ k\\)"):
        MillerBravaisPlane.from_hkil([1, 0, 0, 0], phase=phase)
    with pytest.raises(ValueError, match="U \\+ V \\+ T = 0"):
        MillerBravaisDirection.from_UVTW([1, 0, 0, 0], phase=phase)


def test_metric_tensor_public_helpers_for_common_lattices() -> None:
    cubic = _phase(a=3.0)
    tetragonal = _phase(a=2.0, b=2.0, c=5.0, point_group="4/mmm")
    orthorhombic = _phase(a=2.0, b=3.0, c=4.0, point_group="mmm")
    hexagonal = _phase(a=2.0, b=2.0, c=5.0, gamma=120.0, point_group="6/mmm")

    assert_allclose(metric_tensor(cubic.lattice), np.diag([9.0, 9.0, 9.0]), atol=1e-12)
    assert_allclose(
        reciprocal_metric_tensor(cubic.lattice),
        np.diag([1 / 9, 1 / 9, 1 / 9]),
        atol=1e-12,
    )
    assert_allclose(metric_tensor(tetragonal.lattice), np.diag([4.0, 4.0, 25.0]), atol=1e-12)
    assert_allclose(
        metric_tensor(orthorhombic.lattice),
        np.diag([4.0, 9.0, 16.0]),
        atol=1e-12,
    )
    assert_allclose(
        metric_tensor(hexagonal.lattice),
        np.array([[4.0, -2.0, 0.0], [-2.0, 4.0, 0.0], [0.0, 0.0, 25.0]]),
        atol=1e-12,
    )
    direct = np.array([1.0, 0.0, 0.0])
    reciprocal = direct_to_reciprocal_components(direct, cubic.lattice)
    assert_allclose(reciprocal, np.array([9.0, 0.0, 0.0]), atol=1e-12)
    assert_allclose(reciprocal_to_direct_components(reciprocal, cubic.lattice), direct, atol=1e-12)


def test_ebsd_calibration_geometry_contract_and_frame_checks() -> None:
    specimen = _frame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
    detector = _frame("detector", FrameDomain.DETECTOR, ("u", "v", "n"))
    map_frame = _frame("map", FrameDomain.MAP, ("x", "y", "z"))
    calibration = CalibrationRecord(source="pattern-center-fit", status="calibrated")
    quality = MeasurementQuality(confidence=0.9, valid_fraction=0.95)
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="ebsd",
        map_frame=map_frame,
        detector_frame=detector,
        specimen_to_map=FrameTransform(
            source=specimen,
            target=map_frame,
            rotation_matrix=np.eye(3),
        ),
        specimen_to_detector=FrameTransform(
            source=specimen,
            target=detector,
            rotation_matrix=np.eye(3),
        ),
        calibration_record=calibration,
        measurement_quality=quality,
    )
    detector_geometry = EBSDDetectorGeometry(
        detector_frame=detector,
        pattern_center=(center := PatternCenter(0.48, 0.52, 0.62)),
        detector_distance_mm=18.0,
        pixel_size_um=(50.0, 50.0),
        detector_shape=(480, 640),
        calibration_record=calibration,
        measurement_quality=quality,
        provenance=ProvenanceRecord.minimal("test"),
    )
    geometry = EBSDCalibrationGeometry(
        acquisition_geometry=acquisition,
        detector_geometry=detector_geometry,
        map_to_specimen=FrameTransform(
            source=map_frame,
            target=specimen,
            rotation_matrix=np.eye(3),
        ),
        calibration_record=calibration,
        measurement_quality=quality,
    )

    assert geometry.pattern_center == center
    assert_allclose(detector_geometry.pattern_center_px, np.array([306.72, 249.08]))
    payload = to_json_contract(geometry)
    restored = from_json_contract(payload)
    assert to_json_contract(restored) == payload

    with pytest.raises(ValueError, match="x_fraction"):
        PatternCenter(-0.1, 0.5, 0.6)


def test_synthetic_saed_indexing_and_zone_axis_estimation() -> None:
    phase = _phase()
    specimen = _frame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
    detector = _frame("detector", FrameDomain.DETECTOR, ("u", "v", "n"))
    lab = _frame("lab", FrameDomain.LABORATORY, ("X", "Y", "Z"))
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=200.0,
        pattern_center=np.array([0.5, 0.5, 1.0]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(1024, 1024),
        scattering_setup=ScatteringSetup(laboratory_frame=lab, beam_energy_kev=200.0),
    )
    zone_axis = ZoneAxis([0, 0, 1], phase=phase)
    hkl = np.array(
        [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 1, 0], [-1, -1, 0]],
        dtype=np.int64,
    )
    simulation = KinematicSimulation.simulate_spots(
        geometry,
        phase,
        hkl,
        zone_axis=zone_axis,
        max_excitation_error_inv_angstrom=1.0,
    )
    accepted = simulation.accepted_spots()
    assert len(accepted) >= 4
    pattern = DiffractionPattern(
        coordinates_px=np.stack([spot.detector_coordinates_px for spot in accepted], axis=0),
        intensities=np.array([spot.intensity for spot in accepted], dtype=np.float64),
        geometry=geometry,
        phase=phase,
    )

    estimated = estimate_zone_axis(pattern, max_index=2)
    indexing = index_saed_pattern(
        pattern,
        hkl,
        zone_axis=estimated,
        max_excitation_error_inv_angstrom=1.0,
        max_distance_px=1e-6,
        cluster_radius_px=1e-6,
    )

    assert_array_equal(estimated.indices, np.array([0, 0, 1], dtype=np.int64))
    assert indexing.match_fraction == 1.0
    assert_allclose(indexing.mean_residual_px, 0.0, atol=1e-12)
