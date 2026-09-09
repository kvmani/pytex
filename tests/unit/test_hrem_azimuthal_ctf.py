"""Analytic tests for azimuth-resolved contrast transfer.

Every expected value here has provenance independent of the implementation: it is
either an algebraic identity of the wave aberration function or a symmetry the
aberration order imposes. None is a copied prior output.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.diffraction.hrem import AzimuthalCTF, MicroscopeAberrations

_ENERGY_KEV = 300.0
_DEFOCUS = -50.0
_CS_MM = 0.001
_MAX_Q = 2.5
_RADIAL = 4000


def _round_lens(defocus_angstrom: float) -> MicroscopeAberrations:
    return MicroscopeAberrations(
        energy_kev=_ENERGY_KEV, defocus_angstrom=defocus_angstrom, cs_mm=_CS_MM
    )


def test_round_lens_reports_no_azimuthal_aberration() -> None:
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV, defocus_angstrom=_DEFOCUS, cs_mm=_CS_MM, c5_mm=0.5
    )
    assert lens.has_azimuthal_aberrations is False
    assert lens.residual_aberration_terms() == ()

    band = lens.evaluate_ctf_azimuthal(_MAX_Q, 400, 36)
    # A round lens has the same profile at every azimuth, so the band has no width.
    assert np.allclose(band.transfer_min, band.transfer_max, atol=1e-12)
    assert band.resolution_anisotropy_angstrom == pytest.approx(0.0, abs=1e-12)


def test_two_fold_astigmatism_is_a_defocus_offset_at_its_own_azimuth() -> None:
    """chi gains pi*lambda*q^2*C12*cos(2(theta-phi12)).

    At theta = phi12 the cosine is +1, so the term is algebraically identical to
    adding C12 to the defocus; at theta = phi12 + 90 deg it is -1, subtracting it.
    Both cuts must therefore reproduce a round lens at the offset defocus exactly.
    """
    c12 = 20.0
    phi12 = 30.0
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        astigmatism_angstrom=c12,
        astigmatism_angle_deg=phi12,
    )

    along = lens.evaluate_ctf_1d(_MAX_Q, _RADIAL, azimuth_deg=phi12)
    across = lens.evaluate_ctf_1d(_MAX_Q, _RADIAL, azimuth_deg=phi12 + 90.0)

    assert along.point_resolution_angstrom == pytest.approx(
        _round_lens(_DEFOCUS + c12).evaluate_ctf_1d(_MAX_Q, _RADIAL).point_resolution_angstrom,
        rel=1e-12,
    )
    assert across.point_resolution_angstrom == pytest.approx(
        _round_lens(_DEFOCUS - c12).evaluate_ctf_1d(_MAX_Q, _RADIAL).point_resolution_angstrom,
        rel=1e-12,
    )
    assert along.azimuth_deg == pytest.approx(phi12)


def test_ctf_1d_no_longer_discards_non_round_terms() -> None:
    """A radial cut must respond to astigmatism; before this surface it did not."""
    round_lens = _round_lens(_DEFOCUS)
    astigmatic = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        astigmatism_angstrom=20.0,
        astigmatism_angle_deg=0.0,
    )
    assert not np.allclose(
        round_lens.evaluate_ctf_1d(_MAX_Q, 800).phase_shift_rad,
        astigmatic.evaluate_ctf_1d(_MAX_Q, 800).phase_shift_rad,
    )


def test_astigmatic_resolution_extremes_lie_ninety_degrees_apart() -> None:
    """Two-fold astigmatism has period pi, so best and worst transfer are orthogonal."""
    phi12 = 30.0
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        astigmatism_angstrom=20.0,
        astigmatism_angle_deg=phi12,
    )
    band = lens.evaluate_ctf_azimuthal(_MAX_Q, _RADIAL, 180)
    resolutions = band.point_resolution_angstrom
    best_azimuth = float(band.azimuths_deg[int(np.nanargmin(resolutions))])
    separation = abs(band.worst_azimuth_deg - best_azimuth) % 180.0

    assert separation == pytest.approx(90.0, abs=2.0)
    assert band.resolution_anisotropy_angstrom > 0.0


@pytest.mark.parametrize(
    ("field", "angle_field", "period_deg"),
    [
        ("astigmatism_angstrom", "astigmatism_angle_deg", 180.0),
        ("coma_angstrom", "coma_angle_deg", 360.0),
        ("trefoil_angstrom", "trefoil_angle_deg", 120.0),
    ],
)
def test_each_non_round_term_has_its_own_azimuthal_period(
    field: str, angle_field: str, period_deg: float
) -> None:
    """C12 varies as cos(2 theta), C21 as cos(theta) and C23 as cos(3 theta).

    The period follows from the multiplicity alone, so chi at theta and at
    theta + period must agree exactly, while a half-period apart it must not.
    """
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        **{field: 25.0, angle_field: 0.0},
    )
    q = np.linspace(0.05, _MAX_Q, 200)

    def chi_at(azimuth_deg: float) -> np.ndarray:
        return lens.wave_aberration(q, np.full_like(q, math.radians(azimuth_deg)))

    assert np.allclose(chi_at(17.0), chi_at(17.0 + period_deg), atol=1e-9)
    assert not np.allclose(chi_at(17.0), chi_at(17.0 + period_deg / 2.0), atol=1e-6)


def test_azimuthal_band_bounds_every_sampled_cut() -> None:
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        coma_angstrom=30.0,
        coma_angle_deg=15.0,
        trefoil_angstrom=12.0,
        trefoil_angle_deg=70.0,
    )
    band = lens.evaluate_ctf_azimuthal(_MAX_Q, 600, 24)

    for index, azimuth in enumerate(band.azimuths_deg):
        cut = lens.evaluate_ctf_1d(_MAX_Q, 600, azimuth_deg=float(azimuth))
        assert np.allclose(cut.transfer_function, band.transfer_function[index], atol=1e-12)
        assert np.all(band.transfer_min <= band.transfer_function[index] + 1e-12)
        assert np.all(band.transfer_function[index] <= band.transfer_max + 1e-12)


def test_azimuthal_describe_states_the_terms_and_the_envelope_limit() -> None:
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        astigmatism_angstrom=20.0,
        astigmatism_angle_deg=30.0,
    )
    band = lens.evaluate_ctf_azimuthal(_MAX_Q, 800, 72)
    prose = band.describe()

    assert "C12" in prose
    assert "anisotropic" in prose
    assert "Frank" in prose
    assert "not of the damping" in prose

    lens_prose = lens.describe()
    assert "C12" in lens_prose
    assert "depends on azimuth" in lens_prose


def test_round_lens_describe_says_transfer_is_isotropic() -> None:
    band = _round_lens(_DEFOCUS).evaluate_ctf_azimuthal(_MAX_Q, 400, 12)
    assert "isotropic" in band.describe()


def test_azimuthal_sampling_must_be_positive() -> None:
    with pytest.raises(ValueError, match="num_azimuthal must be at least 1"):
        _round_lens(_DEFOCUS).evaluate_ctf_azimuthal(_MAX_Q, 100, 0)


def test_non_oscillating_lens_reports_no_resolution_rather_than_a_number() -> None:
    """A perfectly corrected lens at zero defocus never crosses zero: NaN, not a value."""
    lens = MicroscopeAberrations(energy_kev=_ENERGY_KEV, defocus_angstrom=0.0, cs_mm=0.0)
    band = lens.evaluate_ctf_azimuthal(_MAX_Q, 400, 8)
    best, worst = band.point_resolution_range_angstrom

    assert math.isnan(best) and math.isnan(worst)
    assert math.isnan(band.resolution_anisotropy_angstrom)
    assert math.isnan(band.worst_azimuth_deg)
    assert "No zero crossing" in band.describe()


def test_azimuthal_ctf_shapes_are_consistent() -> None:
    lens = MicroscopeAberrations(
        energy_kev=_ENERGY_KEV,
        defocus_angstrom=_DEFOCUS,
        cs_mm=_CS_MM,
        astigmatism_angstrom=15.0,
    )
    band = lens.evaluate_ctf_azimuthal(_MAX_Q, 256, 36)

    assert isinstance(band, AzimuthalCTF)
    assert band.spatial_frequencies_inv_angstrom.shape == (256,)
    assert band.azimuths_rad.shape == (36,)
    assert band.phase_shift_rad.shape == (36, 256)
    assert band.transfer_function.shape == (36, 256)
    assert band.total_envelope.shape == (256,)
    assert band.first_zero_q_inv_angstrom.shape == (36,)
    assert float(band.azimuths_deg[0]) == pytest.approx(0.0)
    assert float(band.azimuths_deg[-1]) < 360.0
