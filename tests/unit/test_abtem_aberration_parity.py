"""Optical phase parity with abTEM's documented polar coefficient convention."""

import math

import numpy as np
import pytest

from pytex.adapters.abtem import to_abtem_ctf
from pytex.diffraction.hrem import MicroscopeAberrations


@pytest.mark.parametrize("coefficients", [
    {"coma_angstrom": 150.0, "coma_angle_deg": 37.0},
    {"trefoil_angstrom": -230.0, "trefoil_angle_deg": -21.0},
    {"astigmatism_angstrom": 1e-8, "astigmatism_angle_deg": 45.0},
    {"defocus_angstrom": -150.0, "cs_mm": 0.01, "c5_mm": 2.0,
     "astigmatism_angstrom": 25.0, "astigmatism_angle_deg": 45.0,
     "coma_angstrom": 150.0, "coma_angle_deg": 37.0,
     "trefoil_angstrom": -230.0, "trefoil_angle_deg": -21.0},
])
def test_polar_coefficients_and_phase_transfer(coefficients: dict[str, float]) -> None:
    pytest.importorskip("abtem")
    aberrations = MicroscopeAberrations(energy_kev=300.0, **coefficients)
    ctf = to_abtem_ctf(aberrations)
    assert ctf.C10 == aberrations.defocus_angstrom
    assert ctf.C12 == aberrations.astigmatism_angstrom
    assert ctf.C21 == aberrations.coma_angstrom
    assert ctf.C23 == aberrations.trefoil_angstrom
    assert ctf.C50 == aberrations.c5_angstrom
    assert ctf.phi21 == math.radians(aberrations.coma_angle_deg)
    assert ctf.phi23 == math.radians(aberrations.trefoil_angle_deg)
    q = np.linspace(0.0, 1.0, 15)[:, None]
    theta = np.linspace(-math.pi, math.pi, 17)[None, :]
    q, theta = np.broadcast_arrays(q, theta)
    expected = np.exp(-1j * aberrations.wave_aberration(q, theta))
    # abTEM accepts scattering angles in radians; PyTex accepts reciprocal angstroms.
    actual = ctf._aberrations._evaluate_from_angular_grid(
        q * aberrations.wavelength_angstrom, theta
    )
    np.testing.assert_allclose(actual, expected, atol=2e-5, rtol=2e-5)
