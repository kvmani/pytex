from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

from pytex import (
    AcquisitionGeometry,
    AtomicSite,
    CalibrationRecord,
    CrystalDirection,
    CrystalPlane,
    DiffractionGeometry,
    DiffractionPattern,
    FrameDomain,
    FrameTransform,
    Handedness,
    JSON_CONTRACT_SCHEMA_VERSION,
    Lattice,
    MeasurementQuality,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    ScatteringSetup,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    from_json_contract,
    generate_saed_pattern,
    generate_xrd_pattern,
    get_phase_fixture,
    read_json_contract,
    to_json_contract,
    write_json_contract,
)


def _core_context() -> tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
        metadata={"role": "crystal"},
        provenance=ProvenanceRecord.minimal("test", note="frame fixture"),
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
    lattice = Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group(
        "m-3m",
        reference_frame=crystal,
        provenance=ProvenanceRecord.minimal("test", note="symmetry fixture"),
    )
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite("Ni1", "Ni", np.array([0.0, 0.0, 0.0])),
            AtomicSite("Ni2", "Ni", np.array([0.0, 0.5, 0.5])),
            AtomicSite("Ni3", "Ni", np.array([0.5, 0.0, 0.5])),
            AtomicSite("Ni4", "Ni", np.array([0.5, 0.5, 0.0])),
        ),
        provenance=ProvenanceRecord.minimal("test", note="unit cell fixture"),
    )
    phase = Phase(
        name="nickel-fcc",
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
        chemical_formula="Ni",
        aliases=("Ni", "fcc-nickel"),
        provenance=ProvenanceRecord.minimal("test", note="phase fixture"),
    )
    return crystal, specimen, detector, lab, phase


def _sample_diffraction_geometry() -> tuple[DiffractionGeometry, Phase]:
    crystal, specimen, detector, lab, phase = _core_context()
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
        calibration_record=CalibrationRecord(source="fit", status="calibrated", residual_error=0.1),
        measurement_quality=MeasurementQuality(confidence=0.95, valid_fraction=0.99),
        metadata={"kind": "tem"},
    )
    geometry = DiffractionGeometry(
        detector_frame=detector,
        specimen_frame=specimen,
        laboratory_frame=lab,
        beam_energy_kev=200.0,
        camera_length_mm=150.0,
        pattern_center=np.array([0.5, 0.5, 1.0]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(512, 512),
        acquisition_geometry=acquisition,
        calibration_record=acquisition.calibration_record,
        measurement_quality=acquisition.measurement_quality,
        scattering_setup=ScatteringSetup(laboratory_frame=lab, beam_energy_kev=200.0),
        provenance=ProvenanceRecord.minimal("test", note="geometry fixture"),
    )
    return geometry, phase


def test_core_json_contract_round_trip() -> None:
    crystal, specimen, _, _, phase = _core_context()
    orientation = Orientation(
        rotation=Rotation.from_bunge_euler(35.0, 45.0, 10.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
        provenance=ProvenanceRecord.minimal("test", note="orientation fixture"),
    )
    orientation_set = OrientationSet.from_orientations(
        [
            orientation,
            Orientation(
                rotation=Rotation.from_bunge_euler(15.0, 25.0, 70.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
                provenance=orientation.provenance,
            ),
        ]
    )
    objects = [
        crystal,
        FrameTransform(source=crystal, target=crystal, rotation_matrix=np.eye(3)),
        phase.lattice,
        phase.unit_cell.sites[0],
        phase.unit_cell,
        phase,
        MillerIndex([1, 1, 1], phase=phase),
        CrystalDirection([1.0, 1.0, 0.0], phase=phase),
        CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
        ZoneAxis([1, 1, 1], phase=phase),
        Rotation.from_axis_angle([0.0, 0.0, 1.0], np.deg2rad(30.0)),
        orientation,
        orientation.misorientation_to(orientation, reduce_by_symmetry=False),
        orientation_set,
    ]

    for obj in objects:
        payload = to_json_contract(obj)
        assert payload["schema_version"] == JSON_CONTRACT_SCHEMA_VERSION
        restored = from_json_contract(payload)
        assert to_json_contract(restored) == payload


def test_diffraction_json_contract_round_trip(tmp_path: Path) -> None:
    geometry, phase = _sample_diffraction_geometry()
    zone_axis = ZoneAxis([0, 0, 1], phase=phase)
    saed = generate_saed_pattern(phase, zone_axis, camera_constant_mm_angstrom=180.0, max_index=3)
    xrd = generate_xrd_pattern(phase, max_index=4)
    pattern = DiffractionPattern(
        coordinates_px=np.array([[200.0, 220.0], [256.0, 260.0], [310.0, 240.0]]),
        intensities=np.array([10.0, 4.0, 7.0]),
        geometry=geometry,
        phase=phase,
        orientation=Orientation(
            rotation=Rotation.identity(),
            crystal_frame=phase.crystal_frame,
            specimen_frame=geometry.specimen_frame,
            symmetry=phase.symmetry,
            phase=phase,
        ),
    )
    for obj in (geometry, xrd.radiation, xrd.reflections[0], xrd, saed.spots[0], saed, pattern):
        payload = to_json_contract(obj)
        restored = from_json_contract(payload)
        assert to_json_contract(restored) == payload

    output = tmp_path / "saed_pattern.json"
    write_json_contract(saed, output)
    restored = read_json_contract(output)
    assert to_json_contract(restored) == to_json_contract(saed)


def test_fixture_loaded_phase_round_trips_via_json_contract() -> None:
    crystal = ReferenceFrame(
        "crystal",
        FrameDomain.CRYSTAL,
        ("a", "b", "c"),
        Handedness.RIGHT,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz.*")
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF:.*")
        warnings.filterwarnings("ignore", message="dict interface is deprecated.*")
        phase = get_phase_fixture("zr_hcp").load_phase(crystal_frame=crystal)
    payload = to_json_contract(phase)
    restored = from_json_contract(payload)
    assert to_json_contract(restored) == payload
