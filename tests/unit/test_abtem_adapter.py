"""Unit tests for abTEM adapter layer in pytex.adapters.abtem."""

from __future__ import annotations

import math
from unittest.mock import patch

import numpy as np
import pytest

from pytex.adapters.abtem import (
    from_ase_atoms,
    is_abtem_available,
    simulate_hrem,
    simulate_hrem_multislice,
    to_abtem_ctf,
    to_ase_atoms,
)
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
)


def test_is_abtem_available() -> None:
    # Environment has abtem and ase installed
    assert is_abtem_available() is True


def test_snapshot_and_ase_conversion() -> None:
    snap = AtomicSnapshot.amorphous_sample(
        species="C",
        density_g_cm3=2.0,
        dimensions_angstrom=(12.0, 12.0, 10.0),
        seed=42,
    )
    atoms = to_ase_atoms(snap)
    assert len(atoms) == snap.natoms
    assert np.allclose(atoms.positions, snap.positions)
    assert np.allclose(atoms.cell[:], snap.cell)

    reconstructed = from_ase_atoms(atoms, label="Test ASE")
    assert reconstructed.natoms == snap.natoms
    assert reconstructed.label == "Test ASE"
    assert np.allclose(reconstructed.positions, snap.positions)


@pytest.mark.skipif(not is_abtem_available(), reason="abTEM required")
def test_to_abtem_ctf_parameters() -> None:
    aberr = MicroscopeAberrations(
        energy_kev=300.0,
        mode=DoubleCorrectionMode.DOUBLE_CORRECTED,
        defocus_angstrom=-150.0,
        cs_mm=0.01,
        c5_mm=2.0,
        astigmatism_angstrom=25.0,
        astigmatism_angle_deg=45.0,
        aperture_cutoff_mrad=18.0,
        focal_spread_angstrom=15.0,
        convergence_semiangle_mrad=0.2,
    )
    ctf = to_abtem_ctf(aberr)
    assert ctf.energy == 300000.0
    assert math.isclose(ctf.semiangle_cutoff, 18.0)
    assert math.isclose(ctf.focal_spread, 15.0)


@pytest.mark.skipif(not is_abtem_available(), reason="abTEM required")
def test_simulate_hrem_multislice_end_to_end() -> None:
    snap = AtomicSnapshot.amorphous_sample(
        species="C",
        density_g_cm3=2.0,
        dimensions_angstrom=(10.0, 10.0, 6.0),
        seed=99,
    )
    aberr = MicroscopeAberrations.double_corrected(
        energy_kev=200.0,
        cs_um=0.0,
        defocus_angstrom=-30.0,
        focal_spread_angstrom=10.0,
    )
    result = simulate_hrem_multislice(
        snapshot=snap,
        aberrations=aberr,
        sampling_angstrom=0.2,
        slice_thickness_angstrom=2.0,
    )
    assert result.image.ndim == 2
    ny, nx = result.image.shape
    assert ny >= 20 and nx >= 20
    assert result.exit_wave is not None
    assert result.exit_wave.shape == (ny, nx)
    assert result.power_spectrum.shape == (ny, nx)
    assert result.pixel_size_angstrom > 0.0
    assert "double_corrected" in result.describe()
    assert "Michelson contrast" in result.describe()


def test_simulate_hrem_fallback_dispatch() -> None:
    snap = AtomicSnapshot.amorphous_sample(
        species="C",
        density_g_cm3=1.8,
        dimensions_angstrom=(10.0, 10.0, 6.0),
        seed=10,
    )
    aberr = MicroscopeAberrations.cs_corrected(energy_kev=200.0, cs_um=5.0)

    # Test explicit fallback with prefer_abtem=False
    result_fallback = simulate_hrem(
        snapshot=snap,
        aberrations=aberr,
        sampling_angstrom=0.25,
        prefer_abtem=False,
    )
    assert result_fallback.image.ndim == 2
    assert result_fallback.exit_wave is not None

    # Test mocked abtem unavailable
    with patch("pytex.adapters.abtem.is_abtem_available", return_value=False):
        result_mocked = simulate_hrem(
            snapshot=snap,
            aberrations=aberr,
            sampling_angstrom=0.25,
            prefer_abtem=True,
        )
        assert result_mocked.image.ndim == 2
