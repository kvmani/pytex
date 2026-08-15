from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pytex import (
    CRYSTAL_FRAME,
    Lattice,
    MeasuredPowderPattern,
    Phase,
    PowderPattern,
    PowderPatternComparison,
    RadiationSpec,
    SymmetrySpec,
    compare_powder_patterns,
    from_json_contract,
    read_powder_xy,
    to_json_contract,
    write_powder_xy,
)

ROOT = Path(__file__).resolve().parents[2]
SYNTHETIC_PROFILE = ROOT / "fixtures" / "diffraction" / "synthetic_powder_profile.xy"


def _phase() -> Phase:
    lattice = Lattice(3.5, 3.5, 3.5, 90.0, 90.0, 90.0, crystal_frame=CRYSTAL_FRAME)
    return Phase(
        "analytical test phase",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
    )


def _simulated_profile(intensity: np.ndarray | None = None) -> PowderPattern:
    return PowderPattern(
        phase=_phase(),
        radiation=RadiationSpec.cu_ka(),
        reflections=(),
        two_theta_grid_deg=np.arange(20.0, 25.0),
        intensity_grid=(
            np.arange(1.0, 6.0) if intensity is None else np.asarray(intensity, dtype=float)
        ),
    )


def test_reader_preserves_synthetic_label_units_uncertainty_and_provenance() -> None:
    measured = read_powder_xy(SYNTHETIC_PROFILE, uncertainty_column=2)

    assert measured.name == "synthetic-affine-powder-reference"
    assert measured.synthetic is True
    assert measured.intensity_unit == "counts"
    assert measured.radiation is not None
    assert measured.radiation.name == "Cu Ka"
    assert measured.radiation.wavelength_angstrom == pytest.approx(1.5406)
    assert measured.standard_uncertainty is not None
    assert measured.standard_uncertainty.tolist() == [1.0] * 5
    assert measured.provenance is not None
    assert measured.provenance.metadata["synthetic"] == "true"
    assert "not an experimental measurement" in measured.provenance.notes[0].lower()
    assert "synthetic validation profile" in measured.describe()


def test_affine_comparison_recovers_independently_known_exact_answer() -> None:
    """For y_obs = 5*y_sim + 5, weighted least squares has an exact solution."""

    measured = read_powder_xy(SYNTHETIC_PROFILE, uncertainty_column=2)
    comparison = compare_powder_patterns(measured, _simulated_profile())

    assert comparison.scale_factor == pytest.approx(5.0, abs=1e-12)
    assert comparison.background_offset == pytest.approx(5.0, abs=1e-12)
    assert comparison.profile_r_factor == pytest.approx(0.0, abs=1e-12)
    assert comparison.weighted_profile_r_factor == pytest.approx(0.0, abs=1e-12)
    assert comparison.correlation_coefficient == pytest.approx(1.0, abs=1e-12)
    assert comparison.weight_model == "inverse_variance"
    assert "IUCr pdCIF" in comparison.describe()
    assert "does not shift peaks" in comparison.describe()


def test_profile_r_factors_follow_the_iucr_definitions() -> None:
    measured = MeasuredPowderPattern(
        name="two-point analytical residual",
        two_theta_deg=np.array([20.0, 21.0]),
        intensity=np.array([10.0, 20.0]),
        standard_uncertainty=np.array([1.0, 2.0]),
        synthetic=True,
    )
    simulated = PowderPattern(
        phase=_phase(),
        radiation=RadiationSpec.cu_ka(),
        reflections=(),
        two_theta_grid_deg=np.array([20.0, 21.0]),
        intensity_grid=np.array([9.0, 18.0]),
    )
    comparison = compare_powder_patterns(measured, simulated, fit_background=False)
    calculated = comparison.calculated_intensity
    residual = measured.intensity - calculated
    weights = 1.0 / np.square(measured.standard_uncertainty)
    expected_rp = np.sum(np.abs(residual)) / np.sum(measured.intensity)
    expected_rwp = np.sqrt(
        np.sum(weights * np.square(residual))
        / np.sum(weights * np.square(measured.intensity))
    )

    assert comparison.profile_r_factor == pytest.approx(expected_rp)
    assert comparison.weighted_profile_r_factor == pytest.approx(expected_rwp)


def test_measurement_and_comparison_json_contracts_round_trip() -> None:
    measured = read_powder_xy(SYNTHETIC_PROFILE, uncertainty_column=2)
    comparison = compare_powder_patterns(measured, _simulated_profile())

    for value in (measured, comparison):
        restored = from_json_contract(to_json_contract(value))
        assert isinstance(restored, type(value))
        assert to_json_contract(restored) == to_json_contract(value)
    restored_comparison = from_json_contract(to_json_contract(comparison))
    assert isinstance(restored_comparison, PowderPatternComparison)
    assert restored_comparison.simulated.radiation == RadiationSpec.cu_ka()


def test_writer_round_trips_the_canonical_columns(tmp_path: Path) -> None:
    measured = read_powder_xy(SYNTHETIC_PROFILE, uncertainty_column=2)
    output = write_powder_xy(measured, tmp_path / "round_trip.xy")
    restored = read_powder_xy(output, uncertainty_column=2)

    assert restored.synthetic is True
    assert restored.intensity_unit == measured.intensity_unit
    assert np.array_equal(restored.two_theta_deg, measured.two_theta_deg)
    assert np.array_equal(restored.intensity, measured.intensity)
    assert np.array_equal(restored.standard_uncertainty, measured.standard_uncertainty)


def test_csv_reader_detects_comma_delimiter(tmp_path: Path) -> None:
    source = tmp_path / "instrument_export.csv"
    source.write_text(
        "# name: comma profile\n# intensity_unit: counts_per_second\n20,4\n21,9\n",
        encoding="utf-8",
    )

    measured = read_powder_xy(source)

    assert measured.name == "comma profile"
    assert measured.intensity_unit == "counts_per_second"
    assert measured.intensity.tolist() == [4.0, 9.0]


@pytest.mark.parametrize(
    ("axis", "intensity", "match"),
    [
        ([20.0], [1.0], "at least two"),
        ([20.0, 20.0], [1.0, 2.0], "strictly increasing"),
        ([20.0, 21.0], [0.0, 0.0], "only zero"),
        ([20.0, 21.0], [1.0, -1.0], "non-negative"),
    ],
)
def test_measurement_rejects_ambiguous_or_nonphysical_profiles(
    axis: list[float], intensity: list[float], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        MeasuredPowderPattern("invalid", np.array(axis), np.array(intensity))


def test_comparison_requires_an_overlap_with_positive_simulated_intensity() -> None:
    measured = MeasuredPowderPattern(
        "measured",
        np.array([30.0, 31.0]),
        np.array([1.0, 2.0]),
        synthetic=True,
    )
    with pytest.raises(ValueError, match="overlap"):
        compare_powder_patterns(measured, _simulated_profile())


def test_comparison_rejects_an_inconsistent_portable_result() -> None:
    measured = MeasuredPowderPattern(
        "measured",
        np.array([20.0, 21.0]),
        np.array([1.0, 2.0]),
        synthetic=True,
    )
    with pytest.raises(ValueError, match="Residual intensity"):
        PowderPatternComparison(
            measured=measured,
            simulated=_simulated_profile(np.array([1.0, 2.0, 3.0, 4.0, 5.0])),
            two_theta_deg=np.array([20.0, 21.0]),
            observed_intensity=np.array([1.0, 2.0]),
            calculated_intensity=np.array([1.0, 2.0]),
            residual_intensity=np.array([1.0, 1.0]),
            scale_factor=1.0,
            background_offset=0.0,
            profile_r_factor=0.0,
            weighted_profile_r_factor=0.0,
            correlation_coefficient=1.0,
            weight_model="unit",
            fitted_background=False,
        )
