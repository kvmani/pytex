"""Unit tests for canonical HREM simulation physics and models."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.core.fixtures import get_phase_fixture, phase_fixtures_available
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_interaction_parameter_inv_v_angstrom,
    relativistic_wavelength_angstrom,
)


def test_relativistic_wavelength_against_codata_values() -> None:
    # Check 200 kV, 300 kV, 80 kV against standard electron microscopy reference values
    lam_200 = relativistic_wavelength_angstrom(200.0)
    assert math.isclose(lam_200, 0.025079, rel_tol=1e-4)

    lam_300 = relativistic_wavelength_angstrom(300.0)
    assert math.isclose(lam_300, 0.019687, rel_tol=1e-4)

    lam_80 = relativistic_wavelength_angstrom(80.0)
    assert math.isclose(lam_80, 0.041756, rel_tol=1e-4)

    with pytest.raises(ValueError, match="Accelerating energy must be positive"):
        relativistic_wavelength_angstrom(0.0)
    with pytest.raises(ValueError, match="Accelerating energy must be positive"):
        relativistic_wavelength_angstrom(-200.0)


def test_relativistic_interaction_parameter() -> None:
    sigma_200 = relativistic_interaction_parameter_inv_v_angstrom(200.0)
    assert sigma_200 > 0.0
    # At 200 kV sigma is approx 0.000729 V^-1 Å^-1 (0.729 V^-1 nm^-1)
    assert 0.0005 < sigma_200 < 0.0010


def test_microscope_aberrations_presets_and_scherzer() -> None:
    # Conventional TEM at 200 kV, Cs = 1.0 mm
    conv = MicroscopeAberrations.conventional_tem(energy_kev=200.0, cs_mm=1.0)
    assert conv.mode == DoubleCorrectionMode.UNCORRECTED
    assert conv.cs_mm == 1.0
    assert conv.cs_um == 1000.0
    # Scherzer defocus: Delta f_Sch = -1.2 * sqrt(Cs * lambda)
    # Cs = 1e7 Å, lambda = 0.02508 Å -> sqrt(250790) ~ 500.79 Å -> Delta f_Sch ~ -600.9 Å
    assert math.isclose(conv.scherzer_defocus_angstrom, conv.defocus_angstrom, rel_tol=1e-4)
    assert -650.0 < conv.scherzer_defocus_angstrom < -550.0
    # Scherzer resolution: 0.64 * (Cs * lambda^3)^0.25 ~ 2.4 Å
    assert 2.0 < conv.scherzer_resolution_angstrom < 3.0

    # Cs-corrected preset
    cs_corr = MicroscopeAberrations.cs_corrected(energy_kev=200.0, cs_um=5.0)
    assert cs_corr.mode == DoubleCorrectionMode.CS_CORRECTED
    assert math.isclose(cs_corr.cs_um, 5.0)
    assert cs_corr.cs_mm == 0.005

    # Double-corrected preset
    dbl = MicroscopeAberrations.double_corrected(energy_kev=300.0, cs_um=0.0)
    assert dbl.mode == DoubleCorrectionMode.DOUBLE_CORRECTED
    assert dbl.cs_mm == 0.0
    assert dbl.focal_spread_angstrom == 5.0
    assert "Double-corrected" in dbl.describe()

    # NCSI preset
    ncsi = MicroscopeAberrations.ncsi(energy_kev=200.0, cs_um=-15.0, defocus_angstrom=50.0)
    assert ncsi.mode == DoubleCorrectionMode.NCSI
    assert ncsi.cs_um == -15.0
    assert ncsi.defocus_angstrom == 50.0
    assert "Negative Cs Imaging" in ncsi.describe()


def test_wave_aberration_and_envelopes() -> None:
    aberr = MicroscopeAberrations(
        energy_kev=200.0,
        defocus_angstrom=-400.0,
        cs_mm=1.0,
        astigmatism_angstrom=50.0,
        astigmatism_angle_deg=30.0,
        focal_spread_angstrom=30.0,
        convergence_semiangle_mrad=0.3,
        aperture_cutoff_mrad=20.0,
    )
    q = np.linspace(0.0, 1.5, 100)
    chi = aberr.wave_aberration(q)
    assert chi.shape == (100,)
    assert chi[0] == 0.0  # zero aberration at zero frequency

    # With azimuthal angle
    theta = np.full_like(q, math.radians(30.0))
    chi_theta = aberr.wave_aberration(q, theta=theta)
    assert not np.allclose(chi, chi_theta)

    # Coherence envelopes: must be in (0, 1] and decay monotonically with q
    ec = aberr.temporal_envelope(q)
    es = aberr.spatial_envelope(q)
    assert np.all(ec <= 1.0) and np.all(ec > 0.0)
    assert np.all(es <= 1.0) and np.all(es > 0.0)
    assert ec[0] == 1.0
    assert ec[-1] < 0.1  # strongly damped at 1.5 Å^-1

    # Aperture mask
    ap = aberr.aperture_mask(q)
    assert ap[0] == 1.0
    q_cutoff = (20.0 * 1e-3) / aberr.wavelength_angstrom
    assert np.all(ap[q <= q_cutoff] == 1.0)
    assert np.all(ap[q > q_cutoff] == 0.0)


def test_ctf_1d_evaluation() -> None:
    conv = MicroscopeAberrations.conventional_tem(energy_kev=200.0, cs_mm=1.0)
    ctf = conv.evaluate_ctf_1d(max_q_inv_angstrom=1.5, num_points=300)
    assert len(ctf.spatial_frequencies_inv_angstrom) == 300
    assert not math.isnan(ctf.first_zero_q_inv_angstrom)
    assert 0.3 < ctf.first_zero_q_inv_angstrom < 0.6
    assert 1.8 < ctf.point_resolution_angstrom < 3.0
    assert not math.isnan(ctf.information_limit_angstrom)

    desc = ctf.describe()
    assert "Contrast Transfer Function" in desc
    assert "point resolution" in desc


@pytest.mark.skipif(not phase_fixtures_available(), reason="Phase fixture corpus required")
def test_atomic_snapshot_from_phase_and_defects() -> None:
    from pytex.core import crystal_frame

    phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal_frame("fcc"))
    snap = AtomicSnapshot.from_phase(phase, supercell=(2, 2, 3), zone_axis=(0, 0, 1))
    assert snap.natoms == 2 * 2 * 3 * 4  # 4 atoms per FCC unit cell
    assert snap.positions.shape == (snap.natoms, 3)
    assert snap.cell.shape == (3, 3)
    assert snap.thickness_angstrom > 5.0
    assert "Ni" in snap.species

    # Vacancy creation
    vac = AtomicSnapshot.crystalline_with_vacancy(phase, supercell=(2, 2, 2), vacancy_count=2)
    assert vac.natoms == (2 * 2 * 2 * 4) - 2
    assert "vacancy" in vac.label

    # Dislocation creation
    disl = AtomicSnapshot.crystalline_with_dislocation(
        phase, supercell=(3, 3, 2), dislocation_type="edge"
    )
    assert disl.natoms == (3 * 3 * 2 * 4)
    assert "dislocation" in disl.label


def test_atomic_snapshot_amorphous_and_xyz_roundtrip() -> None:
    amorphous = AtomicSnapshot.amorphous_sample(
        species="C",
        density_g_cm3=2.0,
        dimensions_angstrom=(15.0, 15.0, 10.0),
        seed=123,
    )
    assert amorphous.natoms > 50
    assert all(s == "C" for s in amorphous.species)
    dims = amorphous.dimensions_angstrom
    assert math.isclose(dims[0], 15.0) and math.isclose(dims[1], 15.0)

    # XYZ roundtrip
    xyz_str = amorphous.to_xyz()
    loaded = AtomicSnapshot.from_xyz(xyz_str)
    assert loaded.natoms == amorphous.natoms
    assert np.allclose(loaded.positions, amorphous.positions, atol=1e-4)

def test_atomic_snapshot_ase_roundtrip() -> None:
    pytest.importorskip("ase", reason="ASE is an optional interoperability dependency")
    amorphous = AtomicSnapshot.amorphous_sample(
        species="C", density_g_cm3=2.0, dimensions_angstrom=(15.0, 15.0, 10.0), seed=123,
    )
    # ASE roundtrip
    ase_atoms = amorphous.to_ase()
    assert len(ase_atoms) == amorphous.natoms
    from_ase_snap = AtomicSnapshot.from_ase(ase_atoms)
    assert from_ase_snap.natoms == amorphous.natoms


def test_pure_python_phase_object_simulation_end_to_end() -> None:
    # Build small amorphous carbon sample
    snap = AtomicSnapshot.amorphous_sample(
        species="C",
        density_g_cm3=1.8,
        dimensions_angstrom=(16.0, 16.0, 8.0),
        seed=42,
    )
    aberr = MicroscopeAberrations.double_corrected(
        energy_kev=200.0,
        cs_um=0.0,
        defocus_angstrom=-20.0,
        focal_spread_angstrom=10.0,
    )
    result = pure_python_phase_object_simulation(snap, aberr, sampling_angstrom=0.2)

    assert result.image.ndim == 2
    ny, nx = result.image.shape
    assert nx >= 32 and ny >= 32
    assert result.power_spectrum.shape == (ny, nx)
    assert result.pixel_size_angstrom > 0.0
    assert result.michelson_contrast > 0.0

    # Line profile
    dists, intensities = result.line_profile(num_points=50)
    assert len(dists) == 50
    assert len(intensities) == 50
    assert dists[0] == 0.0
    assert dists[-1] > 0.0

    # Base64 image and power spectrum
    png_b64 = result.to_png_base64()
    assert png_b64.startswith("data:image/png;base64,")
    ps_b64 = result.to_power_spectrum_base64()
    assert ps_b64.startswith("data:image/png;base64,")

    # Explainable description
    desc = result.describe()
    assert "Simulated High-Resolution TEM image" in desc
    assert "Michelson contrast" in desc
