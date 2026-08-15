from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from pytex import (
    CRYSTAL_FRAME,
    SPECIMEN_FRAME,
    Lattice,
    Phase,
    PoleFigure,
    PoleFigureCorrectionSpec,
    PoleFigureDefocusCalibration,
    SymmetrySpec,
    defocus_from_random_standard,
    from_json_contract,
    load_xrdml_pole_figure,
    spherical_angles_to_directions,
    to_json_contract,
)
from pytex.core.lattice import CrystalPlane, MillerIndex

ROOT = Path(__file__).resolve().parents[2]
SYNTHETIC_STANDARD = ROOT / "fixtures" / "xrdml" / "synthetic_random_standard.xrdml"


def _phase() -> Phase:
    lattice = Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=CRYSTAL_FRAME)
    return Phase(
        "synthetic nickel reference",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
    )


def _pole(indices: tuple[int, int, int] = (1, 1, 1)) -> CrystalPlane:
    phase = _phase()
    return CrystalPlane(MillerIndex(indices, phase=phase), phase=phase)


def _standard() -> PoleFigure:
    return load_xrdml_pole_figure(
        SYNTHETIC_STANDARD,
        pole=_pole(),
        specimen_frame=SPECIMEN_FRAME,
        antipodal=True,
        intensity_normalization="none",
    )


def test_random_standard_recovers_exact_background_subtracted_curve() -> None:
    """(I_ring - 10) / (I_0 - 10) gives exactly [1, 0.8, 0.5]."""

    calibration = defocus_from_random_standard(_standard(), background=10.0, synthetic=True)

    assert calibration.tilt_deg.tolist() == pytest.approx([0.0, 30.0, 60.0])
    assert calibration.ring_intensities.tolist() == pytest.approx([100.0, 80.0, 50.0])
    assert calibration.defocus_factors.tolist() == pytest.approx([1.0, 0.8, 0.5])
    assert calibration.ring_counts.tolist() == [4, 4, 4]
    assert calibration.max_azimuthal_relative_std == pytest.approx(0.0)
    assert calibration.synthetic is True
    assert "synthetic validation standard" in calibration.describe()
    assert "refuses extrapolation" in calibration.describe()


def test_calibrated_correction_subtracts_background_before_defocus_division() -> None:
    calibration = defocus_from_random_standard(_standard(), background=10.0, synthetic=True)
    tilt = np.array([0.0, 15.0, 30.0, 45.0, 60.0])
    directions = spherical_angles_to_directions(tilt, np.zeros_like(tilt))
    factors = np.interp(tilt, [0.0, 30.0, 60.0], [1.0, 0.8, 0.5])
    true_density = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
    target_background = 2.0
    observed = target_background + factors * true_density
    target = PoleFigure(
        pole=_pole(),
        sample_directions=directions,
        intensities=observed,
        specimen_frame=SPECIMEN_FRAME,
        antipodal=True,
        sampling="sampled_density",
    )

    correction = calibration.correction_spec(
        target, background=target_background, missing_intensity_policy="raise"
    )
    corrected = correction.apply(target)

    assert corrected.intensities.tolist() == pytest.approx(true_density, abs=1e-12)
    # This one point distinguishes (I - b) / d from the old, incorrect I / d - b order.
    direct = PoleFigureCorrectionSpec(
        background=2.0,
        defocus_factors=np.full(target.intensities.shape, 0.5),
        missing_intensity_policy="raise",
    ).apply(
        PoleFigure(
            pole=target.pole,
            sample_directions=target.sample_directions,
            intensities=np.full(target.intensities.shape, 12.0),
            specimen_frame=target.specimen_frame,
            antipodal=True,
            sampling="sampled_density",
        )
    )
    assert direct.intensities.tolist() == pytest.approx([20.0] * 5)


def test_defocus_calibration_json_contract_round_trips() -> None:
    calibration = defocus_from_random_standard(_standard(), background=10.0, synthetic=True)

    payload = to_json_contract(calibration)
    restored = from_json_contract(payload)

    assert isinstance(restored, PoleFigureDefocusCalibration)
    assert to_json_contract(restored) == payload

    with pytest.raises(ValueError, match="ring intensities divided"):
        replace(calibration, defocus_factors=np.array([1.0, 0.7, 0.5]))


def test_calibration_refuses_wrong_reflection_and_extrapolation() -> None:
    calibration = defocus_from_random_standard(_standard(), background=10.0, synthetic=True)
    outside = PoleFigure(
        pole=_pole(),
        sample_directions=spherical_angles_to_directions([75.0], [0.0]),
        intensities=np.ones(1),
        specimen_frame=SPECIMEN_FRAME,
        antipodal=True,
        sampling="sampled_density",
    )
    wrong_pole = PoleFigure(
        pole=_pole((2, 0, 0)),
        sample_directions=spherical_angles_to_directions([30.0], [0.0]),
        intensities=np.ones(1),
        specimen_frame=SPECIMEN_FRAME,
        antipodal=True,
        sampling="sampled_density",
    )

    with pytest.raises(ValueError, match="extrapolation"):
        calibration.factors_for(outside)
    with pytest.raises(ValueError, match="reflection-specific"):
        calibration.factors_for(wrong_pole)


def test_calibration_rejects_pole_weights_and_nonpositive_corrected_standard() -> None:
    standard = _standard()
    scattered = PoleFigure(
        pole=standard.pole,
        sample_directions=standard.sample_directions,
        intensities=standard.intensities,
        specimen_frame=standard.specimen_frame,
        antipodal=standard.antipodal,
        sampling="scattered_poles",
    )

    with pytest.raises(ValueError, match="sampled_density"):
        defocus_from_random_standard(scattered)
    with pytest.raises(ValueError, match="remain positive"):
        defocus_from_random_standard(standard, background=110.0)
