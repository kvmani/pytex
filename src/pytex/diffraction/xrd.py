from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord

_TWO_THETA_TOLERANCE_DEG = 1e-8


def _gaussian_profile(x: np.ndarray, center: float, sigma: float) -> np.ndarray:
    exponent = -0.5 * ((x - center) / sigma) ** 2
    profile = np.exp(exponent) / (sigma * np.sqrt(2.0 * np.pi))
    profile = np.ascontiguousarray(profile, dtype=np.float64)
    profile.setflags(write=False)
    return profile


def _lorentz_polarization(two_theta_rad: float) -> float:
    theta = 0.5 * two_theta_rad
    sin_theta = max(float(np.sin(theta)), 1e-8)
    cos_theta = max(float(np.cos(theta)), 1e-8)
    cos_two_theta = float(np.cos(two_theta_rad))
    return float((1.0 + cos_two_theta * cos_two_theta) / (sin_theta * sin_theta * cos_theta))


def _structure_factor_xray(phase: Phase, hkl: np.ndarray) -> complex:
    if phase.unit_cell is None or not phase.unit_cell.sites:
        return complex(1.0, 0.0)
    reciprocal = phase.lattice.reciprocal_basis().matrix
    g_cart = reciprocal @ hkl.astype(np.float64)
    g_sq = float(np.dot(g_cart, g_cart))
    total = 0.0j
    for site in phase.unit_cell.sites:
        z = float(atomic_number(site.species))
        occupancy = float(site.occupancy)
        b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
        debye_waller = np.exp(-(b_iso * g_sq) / max(16.0 * np.pi * np.pi, 1e-12))
        phase_factor = np.exp(
            2.0j * np.pi * float(np.dot(hkl.astype(np.float64), site.fractional_coordinates))
        )
        total += occupancy * z * debye_waller * phase_factor
    return complex(total)


def _equivalent_hkls(phase: Phase, hkl: np.ndarray) -> tuple[np.ndarray, ...]:
    reciprocal_basis = phase.lattice.reciprocal_basis().matrix
    inverse_basis = np.linalg.inv(reciprocal_basis)
    g_cart = reciprocal_basis @ hkl.astype(np.float64)
    equivalents: set[tuple[int, int, int]] = set()
    for operator in phase.symmetry.operators:
        transformed = operator @ g_cart
        transformed_hkl = inverse_basis @ transformed
        rounded = np.rint(transformed_hkl).astype(np.int64)
        if np.allclose(transformed_hkl, rounded, atol=1e-6):
            equivalents.add((int(rounded[0]), int(rounded[1]), int(rounded[2])))
            negatives = tuple([-int(value) for value in rounded])
            equivalents.add((negatives[0], negatives[1], negatives[2]))
    return tuple(np.array(item, dtype=np.int64) for item in sorted(equivalents))


def _reflection_multiplicity(phase: Phase, hkl: np.ndarray) -> int:
    return len(_equivalent_hkls(phase, hkl))


def _enumerate_hkls(max_index: int) -> np.ndarray:
    values = range(-max_index, max_index + 1)
    hkls = [
        (h, k, ell)
        for h in values
        for k in values
        for ell in values
        if not (h == 0 and k == 0 and ell == 0)
    ]
    array = np.asarray(hkls, dtype=np.int64)
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class RadiationSpec:
    name: str
    wavelength_angstrom: float

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("RadiationSpec.name must be non-empty.")
        if self.wavelength_angstrom <= 0.0:
            raise ValueError("RadiationSpec.wavelength_angstrom must be strictly positive.")

    @classmethod
    def cu_ka(cls) -> RadiationSpec:
        return cls(name="Cu Ka", wavelength_angstrom=1.5406)

    @classmethod
    def mo_ka(cls) -> RadiationSpec:
        return cls(name="Mo Ka", wavelength_angstrom=0.71073)


@dataclass(frozen=True, slots=True)
class PowderReflection:
    miller_indices: np.ndarray
    d_spacing_angstrom: float
    two_theta_deg: float
    intensity: float
    structure_factor_amplitude: float
    multiplicity: int
    structure_factor_real: float | None = None
    structure_factor_imag: float | None = None
    lorentz_polarization_factor: float | None = None
    intensity_model: str = "xray_atomic_number"

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        if self.d_spacing_angstrom <= 0.0:
            raise ValueError("PowderReflection.d_spacing_angstrom must be strictly positive.")
        if self.two_theta_deg < 0.0:
            raise ValueError("PowderReflection.two_theta_deg must be non-negative.")
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("PowderReflection.intensity must be finite and non-negative.")
        if (
            not np.isfinite(self.structure_factor_amplitude)
            or self.structure_factor_amplitude < 0.0
        ):
            raise ValueError(
                "PowderReflection.structure_factor_amplitude must be finite and non-negative."
            )
        if self.multiplicity <= 0:
            raise ValueError("PowderReflection.multiplicity must be strictly positive.")
        for name, value in (
            ("structure_factor_real", self.structure_factor_real),
            ("structure_factor_imag", self.structure_factor_imag),
            ("lorentz_polarization_factor", self.lorentz_polarization_factor),
        ):
            if value is not None and not np.isfinite(value):
                raise ValueError(f"PowderReflection.{name} must be finite when provided.")
        if not self.intensity_model.strip():
            raise ValueError("PowderReflection.intensity_model must be non-empty.")


@dataclass(frozen=True, slots=True)
class PowderPattern:
    phase: Phase
    radiation: RadiationSpec
    reflections: tuple[PowderReflection, ...]
    two_theta_grid_deg: np.ndarray
    intensity_grid: np.ndarray
    broadening_fwhm_deg: float | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "reflections", tuple(self.reflections))
        object.__setattr__(
            self,
            "two_theta_grid_deg",
            as_float_array(self.two_theta_grid_deg, shape=(None,)),
        )
        object.__setattr__(
            self,
            "intensity_grid",
            as_float_array(self.intensity_grid, shape=(None,)),
        )
        if self.two_theta_grid_deg.shape != self.intensity_grid.shape:
            raise ValueError("PowderPattern grid arrays must have the same shape.")
        if np.any(~np.isfinite(self.intensity_grid)) or np.any(self.intensity_grid < 0.0):
            raise ValueError("PowderPattern intensity_grid must be finite and non-negative.")


def generate_powder_reflections(
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    two_theta_range_deg: tuple[float, float] = (5.0, 120.0),
    max_index: int = 6,
    intensity_model: Literal["xray_atomic_number", "unit"] = "xray_atomic_number",
) -> tuple[PowderReflection, ...]:
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    radiation_spec = RadiationSpec.cu_ka() if radiation is None else radiation
    min_two_theta, max_two_theta = (float(two_theta_range_deg[0]), float(two_theta_range_deg[1]))
    if not 0.0 <= min_two_theta < max_two_theta <= 180.0:
        raise ValueError("two_theta_range_deg must satisfy 0 <= min < max <= 180.")

    reflections: list[PowderReflection] = []
    for hkl in _enumerate_hkls(max_index):
        d_spacing = phase.lattice.reciprocal_basis().matrix @ hkl.astype(np.float64)
        g_magnitude = float(np.linalg.norm(d_spacing))
        if np.isclose(g_magnitude, 0.0):
            continue
        d_spacing_angstrom = 1.0 / g_magnitude
        argument = radiation_spec.wavelength_angstrom / (2.0 * d_spacing_angstrom)
        if argument <= 0.0 or argument > 1.0:
            continue
        two_theta_deg = float(np.rad2deg(2.0 * np.arcsin(argument)))
        if not (
            min_two_theta - _TWO_THETA_TOLERANCE_DEG
            <= two_theta_deg
            <= max_two_theta + _TWO_THETA_TOLERANCE_DEG
        ):
            continue
        if intensity_model == "unit":
            structure_factor = complex(1.0, 0.0)
        else:
            structure_factor = _structure_factor_xray(phase, hkl)
        amplitude = float(abs(structure_factor))
        multiplicity = _reflection_multiplicity(phase, hkl)
        lorentz_polarization = _lorentz_polarization(np.deg2rad(two_theta_deg))
        intensity = multiplicity * amplitude * amplitude * lorentz_polarization
        if intensity <= 1e-14:
            continue
        reflections.append(
            PowderReflection(
                miller_indices=hkl,
                d_spacing_angstrom=d_spacing_angstrom,
                two_theta_deg=two_theta_deg,
                intensity=intensity,
                structure_factor_amplitude=amplitude,
                multiplicity=multiplicity,
                structure_factor_real=float(structure_factor.real),
                structure_factor_imag=float(structure_factor.imag),
                lorentz_polarization_factor=lorentz_polarization,
                intensity_model=intensity_model,
            )
        )
    reflections.sort(
        key=lambda reflection: (
            reflection.two_theta_deg,
            tuple(int(value) for value in reflection.miller_indices),
        )
    )
    return tuple(reflections)


def generate_xrd_pattern(
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    two_theta_range_deg: tuple[float, float] = (5.0, 120.0),
    resolution_deg: float = 0.02,
    max_index: int = 6,
    intensity_model: Literal["xray_atomic_number", "unit"] = "xray_atomic_number",
    broadening_fwhm_deg: float | None = 0.15,
    provenance: ProvenanceRecord | None = None,
) -> PowderPattern:
    if resolution_deg <= 0.0:
        raise ValueError("resolution_deg must be strictly positive.")
    if broadening_fwhm_deg is not None and broadening_fwhm_deg <= 0.0:
        raise ValueError("broadening_fwhm_deg must be strictly positive when provided.")
    radiation_spec = RadiationSpec.cu_ka() if radiation is None else radiation
    reflections = generate_powder_reflections(
        phase,
        radiation=radiation_spec,
        two_theta_range_deg=two_theta_range_deg,
        max_index=max_index,
        intensity_model=intensity_model,
    )
    min_two_theta, max_two_theta = map(float, two_theta_range_deg)
    grid = np.arange(min_two_theta, max_two_theta + 0.5 * resolution_deg, resolution_deg)
    intensity_grid = np.zeros_like(grid)
    if broadening_fwhm_deg is None:
        for reflection in reflections:
            index = int(np.argmin(np.abs(grid - reflection.two_theta_deg)))
            intensity_grid[index] += reflection.intensity
    else:
        sigma = broadening_fwhm_deg / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        for reflection in reflections:
            intensity_grid += reflection.intensity * _gaussian_profile(
                grid,
                reflection.two_theta_deg,
                sigma,
            )
        if intensity_grid.size > 0:
            intensity_grid *= resolution_deg
    if np.max(intensity_grid) > 0.0:
        intensity_grid /= float(np.max(intensity_grid))
    intensity_grid = np.ascontiguousarray(intensity_grid, dtype=np.float64)
    intensity_grid.setflags(write=False)
    grid = np.ascontiguousarray(grid, dtype=np.float64)
    grid.setflags(write=False)
    return PowderPattern(
        phase=phase,
        radiation=radiation_spec,
        reflections=reflections,
        two_theta_grid_deg=grid,
        intensity_grid=intensity_grid,
        broadening_fwhm_deg=broadening_fwhm_deg,
        provenance=provenance,
    )
