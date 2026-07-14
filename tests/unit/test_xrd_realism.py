from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
    generate_powder_reflections,
    generate_xrd_pattern,
)
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import AtomicSite, UnitCell
from pytex.diffraction.scattering import (
    tabulated_species,
    xray_form_factor_matrix,
    xray_form_factors,
)


def make_nickel() -> Phase:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.524, 3.524, 3.524, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    fcc_fractionals = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
    unit_cell = UnitCell(
        lattice=lattice,
        sites=tuple(
            AtomicSite(label=f"Ni{i}", species="Ni", fractional_coordinates=np.array(frac))
            for i, frac in enumerate(fcc_fractionals)
        ),
    )
    return Phase(
        "nickel", lattice=lattice, symmetry=symmetry, crystal_frame=crystal, unit_cell=unit_cell
    )


def test_form_factors_zero_angle_limit_is_atomic_number() -> None:
    for species in ("Ni", "Fe", "Si", "O", "Zr"):
        assert species in tabulated_species()
        value = xray_form_factors(species, np.array([0.0]))[0]
        assert value == pytest.approx(float(atomic_number(species)), abs=1e-9)


def test_form_factors_decrease_with_scattering_vector() -> None:
    s_values = np.linspace(0.0, 1.0, 40)
    for species in ("Ni", "Fe", "C"):
        factors = xray_form_factors(species, s_values)
        assert np.all(np.diff(factors) < 0.0)
        assert factors[-1] > 0.0


def test_form_factor_matrix_shape_and_fallback() -> None:
    s_values = np.array([0.0, 0.25, 0.5])
    matrix = xray_form_factor_matrix(("Ni", "O"), s_values)
    assert matrix.shape == (3, 2)
    # an element missing from the table falls back to constant Z
    unknown = xray_form_factors("Es", s_values) if "Es" not in tabulated_species() else None
    if unknown is not None:
        assert np.all(unknown == unknown[0])


def test_tabulated_intensities_weaken_high_angle_reflections() -> None:
    nickel = make_nickel()
    z_proxy = generate_powder_reflections(nickel, intensity_model="xray_atomic_number")
    tabulated = generate_powder_reflections(nickel, intensity_model="xray_tabulated")
    assert len(z_proxy) == len(tabulated)
    # same reflection set, same geometry
    for a, b in zip(z_proxy, tabulated, strict=True):
        assert tuple(a.miller_indices) == tuple(b.miller_indices)
        assert a.two_theta_deg == pytest.approx(b.two_theta_deg)
    # relative to the first peak, angle-dependent form factors suppress
    # high-angle intensities compared with the constant-Z proxy
    ratio_z = z_proxy[-1].intensity / z_proxy[0].intensity
    ratio_tab = tabulated[-1].intensity / tabulated[0].intensity
    assert ratio_tab < ratio_z


def test_kalpha_doublet_positions_follow_braggs_law() -> None:
    nickel = make_nickel()
    radiation = RadiationSpec.cu_ka_doublet()
    pattern = generate_xrd_pattern(
        nickel,
        radiation=radiation,
        two_theta_range_deg=(40.0, 55.0),
        resolution_deg=0.005,
        broadening_fwhm_deg=0.05,
    )
    (reflection,) = (r for r in pattern.reflections if tuple(r.miller_indices) == (1, 1, 1))
    # analytic Ka2 position from Bragg's law
    d = reflection.d_spacing_angstrom
    two_theta_ka2 = np.degrees(
        2.0 * np.arcsin(radiation.kalpha2_wavelength_angstrom / (2.0 * d))
    )
    assert two_theta_ka2 > reflection.two_theta_deg
    window = (pattern.two_theta_grid_deg > reflection.two_theta_deg + 0.05) & (
        pattern.two_theta_grid_deg < two_theta_ka2 + 0.2
    )
    ka2_peak_position = pattern.two_theta_grid_deg[window][
        int(np.argmax(pattern.intensity_grid[window]))
    ]
    assert ka2_peak_position == pytest.approx(two_theta_ka2, abs=0.02)
    # Ka2 peak height ~ half the Ka1 height for a resolved doublet
    ka1_height = pattern.intensity_grid[
        int(np.argmin(np.abs(pattern.two_theta_grid_deg - reflection.two_theta_deg)))
    ]
    ka2_height = pattern.intensity_grid[
        int(np.argmin(np.abs(pattern.two_theta_grid_deg - ka2_peak_position)))
    ]
    assert ka2_height / ka1_height == pytest.approx(0.5, abs=0.08)


def test_pseudo_voigt_limits_and_tails() -> None:
    nickel = make_nickel()
    common = {
        "two_theta_range_deg": (40.0, 50.0),
        "resolution_deg": 0.01,
        "broadening_fwhm_deg": 0.2,
    }
    gaussian = generate_xrd_pattern(nickel, profile="gaussian", **common)
    pv_zero = generate_xrd_pattern(nickel, profile="pseudo_voigt", pseudo_voigt_eta=0.0, **common)
    np.testing.assert_allclose(pv_zero.intensity_grid, gaussian.intensity_grid, atol=1e-12)
    lorentzian = generate_xrd_pattern(
        nickel, profile="pseudo_voigt", pseudo_voigt_eta=1.0, **common
    )
    # Lorentzian tails carry far more weight than Gaussian at 1 deg off-peak
    (reflection,) = (r for r in gaussian.reflections if tuple(r.miller_indices) == (1, 1, 1))
    tail_index = int(
        np.argmin(np.abs(gaussian.two_theta_grid_deg - (reflection.two_theta_deg + 1.0)))
    )
    assert lorentzian.intensity_grid[tail_index] > 100.0 * gaussian.intensity_grid[tail_index]


def test_caglioti_widths_grow_with_angle() -> None:
    nickel = make_nickel()
    pattern = generate_xrd_pattern(
        nickel,
        two_theta_range_deg=(40.0, 100.0),
        resolution_deg=0.01,
        broadening_fwhm_deg=None,
        caglioti_uvw=(0.6, 0.0, 0.01),
    )

    def peak_fwhm(two_theta: float) -> float:
        center = int(np.argmin(np.abs(pattern.two_theta_grid_deg - two_theta)))
        height = pattern.intensity_grid[center]
        above = pattern.intensity_grid >= 0.5 * height
        # walk outward from the center to the half-maximum crossings
        left = center
        while left > 0 and above[left - 1]:
            left -= 1
        right = center
        while right < above.size - 1 and above[right + 1]:
            right += 1
        return float(
            pattern.two_theta_grid_deg[right] - pattern.two_theta_grid_deg[left]
        )

    ordered = sorted(pattern.reflections, key=lambda r: r.two_theta_deg)
    low, high = ordered[0], ordered[-1]
    assert peak_fwhm(high.two_theta_deg) > peak_fwhm(low.two_theta_deg)


def test_radiation_spec_validation_and_presets() -> None:
    for preset in (
        RadiationSpec.cu_ka_doublet(),
        RadiationSpec.mo_ka_doublet(),
        RadiationSpec.co_ka(),
        RadiationSpec.cr_ka(),
        RadiationSpec.fe_ka(),
    ):
        assert preset.kalpha2_wavelength_angstrom > preset.wavelength_angstrom
        assert preset.kalpha2_relative_intensity == 0.5
        assert preset.anode is not None
    neutron = RadiationSpec.neutron(1.798)
    assert neutron.kind == "neutron"
    with pytest.raises(ValueError, match="kind"):
        RadiationSpec(name="x", wavelength_angstrom=1.0, kind="gamma")
    with pytest.raises(ValueError, match="relative_intensity"):
        RadiationSpec(
            name="x",
            wavelength_angstrom=1.0,
            kalpha2_wavelength_angstrom=1.1,
            kalpha2_relative_intensity=0.0,
        )
    # legacy presets unchanged (pinned external baselines depend on them)
    assert RadiationSpec.cu_ka().wavelength_angstrom == 1.5406
    assert RadiationSpec.cu_ka().kalpha2_wavelength_angstrom is None
