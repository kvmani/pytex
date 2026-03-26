from __future__ import annotations

import numpy as np
import pytest

from pytex.core import (
    AcquisitionGeometry,
    CalibrationRecord,
    FrameDomain,
    FrameTransform,
    Handedness,
    MeasurementQuality,
    ReferenceFrame,
    ScatteringSetup,
)


def make_frames() -> tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame, ReferenceFrame]:
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    map_frame = ReferenceFrame(
        name="map",
        domain=FrameDomain.MAP,
        axes=("i", "j", "k"),
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
    return specimen, map_frame, detector, lab


def test_measurement_quality_validates_unit_interval_fields() -> None:
    with pytest.raises(ValueError):
        MeasurementQuality(confidence=1.5)


def test_calibration_record_normalizes_status() -> None:
    record = CalibrationRecord(source="vendor-fit", status="Calibrated", residual_error=0.1)
    assert record.status == "calibrated"
    assert record.is_calibrated


def test_scattering_setup_derives_effective_wavelength_from_energy() -> None:
    _, _, _, lab = make_frames()
    setup = ScatteringSetup(laboratory_frame=lab, beam_energy_kev=20.0)
    assert setup.effective_wavelength_angstrom > 0.0


def test_acquisition_geometry_requires_transforms_for_distinct_frames() -> None:
    specimen, map_frame, _, _ = make_frames()
    with pytest.raises(ValueError):
        AcquisitionGeometry(specimen_frame=specimen, modality="ebsd", map_frame=map_frame)


def test_acquisition_geometry_composes_specimen_to_reciprocal_transform() -> None:
    specimen, _, _, lab = make_frames()
    reciprocal = ReferenceFrame(
        name="reciprocal",
        domain=FrameDomain.RECIPROCAL,
        axes=("a*", "b*", "c*"),
        handedness=Handedness.RIGHT,
    )
    specimen_to_lab = FrameTransform(
        source=specimen,
        target=lab,
        rotation_matrix=np.eye(3),
    )
    lab_to_reciprocal = FrameTransform(
        source=lab,
        target=reciprocal,
        rotation_matrix=np.eye(3),
    )
    acquisition = AcquisitionGeometry(
        specimen_frame=specimen,
        modality="xrd",
        laboratory_frame=lab,
        reciprocal_frame=reciprocal,
        specimen_to_laboratory=specimen_to_lab,
        laboratory_to_reciprocal=lab_to_reciprocal,
    )
    transform = acquisition.specimen_to_reciprocal_transform()
    assert transform is not None
    assert transform.source == specimen
    assert transform.target == reciprocal
