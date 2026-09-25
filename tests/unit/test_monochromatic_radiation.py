"""A monochromatic beam of any wavelength, and the polarization of a synchrotron.

Expected values are closed forms: hc = 12.398419843 keV Å (CODATA 2018), and
the polarization factor f + (1 - f) cos²2θ of a beam whose fraction f is
polarized perpendicular to the scattering plane.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.contracts import from_json_contract, to_json_contract
from pytex.diffraction.physics import lorentz_polarization_factor
from pytex.diffraction.xrd import (
    HC_KEV_ANGSTROM,
    RadiationSpec,
    _lorentz_polarization,
    generate_powder_reflections,
)


def test_energy_and_wavelength_convert_through_hc() -> None:
    beam = RadiationSpec.from_energy_kev(30.0)
    assert beam.wavelength_angstrom == pytest.approx(12.398419843 / 30.0, rel=1e-10)
    assert beam.energy_kev == pytest.approx(30.0)
    assert HC_KEV_ANGSTROM == pytest.approx(12.398419843320026)
    assert RadiationSpec.cu_ka().energy_kev == pytest.approx(8.04, abs=0.01)


def test_a_monochromatic_beam_has_one_line_and_says_what_it_is() -> None:
    beam = RadiationSpec.monochromatic(0.4133)
    assert beam.is_monochromatic and beam.kalpha2_wavelength_angstrom is None
    assert beam.polarization_perpendicular_fraction == 0.5
    assert "0.4133" in beam.name and "keV" in beam.name
    synchrotron = RadiationSpec.synchrotron(0.2066)
    assert synchrotron.polarization_perpendicular_fraction == 0.95
    assert not RadiationSpec.cu_ka_doublet().is_monochromatic
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError, match="wavelength"):
            RadiationSpec.monochromatic(bad)
    with pytest.raises(ValueError, match="polarization"):
        RadiationSpec.monochromatic(0.5, polarization_perpendicular_fraction=1.2)


@pytest.mark.parametrize("two_theta", [20.0, 60.0, 90.0, 130.0])
def test_the_polarization_factor_follows_the_perpendicular_fraction(two_theta: float) -> None:
    angle = math.radians(two_theta)
    theta = angle / 2
    lorentz = 1.0 / (math.sin(theta) ** 2 * math.cos(theta))
    cos2 = math.cos(angle) ** 2
    # f = 1/2 is the textbook unpolarized factor, exactly.
    assert _lorentz_polarization(angle) == pytest.approx((1 + cos2) * lorentz, rel=1e-12)
    assert lorentz_polarization_factor(angle) == pytest.approx((1 + cos2) * lorentz, rel=1e-12)
    for fraction in (0.0, 0.05, 0.95, 1.0):
        expected = 2.0 * (fraction + (1 - fraction) * cos2) * lorentz
        assert _lorentz_polarization(angle, fraction) == pytest.approx(expected, rel=1e-12)
        assert lorentz_polarization_factor(angle, perpendicular_fraction=fraction) == (
            pytest.approx(expected, rel=1e-12)
        )


def test_a_horizontal_scattering_plane_extinguishes_the_ninety_degree_reflection() -> None:
    from pytex.app.phases import builtin_phase

    nickel = builtin_phase("ni_fcc").to_phase()
    # A wavelength that puts (220) at 2θ = 90° exactly: λ = 2 d sin 45°.
    d220 = 3.52387 / math.sqrt(8.0)
    wavelength = 2.0 * d220 * math.sin(math.radians(45.0))
    reflections = {}
    for fraction in (0.5, 0.0, 1.0):
        beam = RadiationSpec.monochromatic(wavelength, polarization_perpendicular_fraction=fraction)
        found = generate_powder_reflections(nickel, radiation=beam, two_theta_range_deg=(10, 170))
        reflections[fraction] = {tuple(sorted(np.abs(r.miller_indices))): r for r in found}
    key = (0, 2, 2)
    assert reflections[0.5][key].two_theta_deg == pytest.approx(90.0, abs=1e-9)
    # P(90°) = f, so f = 0 removes it, f = 1 doubles it relative to unpolarized.
    assert key not in reflections[0.0] or reflections[0.0][key].intensity < 1e-9
    assert reflections[1.0][key].intensity == pytest.approx(
        2.0 * reflections[0.5][key].intensity, rel=1e-9
    )


def test_the_polarization_survives_the_json_contract() -> None:
    beam = RadiationSpec.synchrotron(0.3)
    restored = from_json_contract(to_json_contract(beam))
    assert restored == beam
    legacy = to_json_contract(RadiationSpec.cu_ka())
    legacy.pop("polarization_perpendicular_fraction")
    assert from_json_contract(legacy).polarization_perpendicular_fraction == 0.5
