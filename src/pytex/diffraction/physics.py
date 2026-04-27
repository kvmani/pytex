from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.xrd import RadiationSpec

ScatteringModelName = Literal["unit", "atomic_number", "xray_gaussian_proxy"]
IntensityCorrectionName = Literal["none", "lorentz_polarization"]


def lorentz_polarization_factor(two_theta_rad: float) -> float:
    theta = 0.5 * float(two_theta_rad)
    sin_theta = max(float(np.sin(theta)), 1e-8)
    cos_theta = max(float(np.cos(theta)), 1e-8)
    cos_two_theta = float(np.cos(two_theta_rad))
    return float((1.0 + cos_two_theta * cos_two_theta) / (sin_theta * sin_theta * cos_theta))


@dataclass(frozen=True, slots=True)
class ScatteringFactorTable:
    model: ScatteringModelName = "atomic_number"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.model not in {"unit", "atomic_number", "xray_gaussian_proxy"}:
            raise ValueError(
                "ScatteringFactorTable.model must be 'unit', 'atomic_number', or "
                "'xray_gaussian_proxy'."
            )

    def scattering_factor(self, species: str, g_magnitude_inv_angstrom: float) -> float:
        if self.model == "unit":
            return 1.0
        z = float(atomic_number(species))
        if self.model == "atomic_number":
            return z
        # Smooth monotonic proxy. This is not a tabulated Cromer-Mann replacement.
        return float(z * np.exp(-0.02 * g_magnitude_inv_angstrom * g_magnitude_inv_angstrom))


@dataclass(frozen=True, slots=True)
class StructureFactor:
    miller_indices: np.ndarray
    value: complex
    amplitude: float
    phase_rad: float
    phase: Phase
    scattering_table: ScatteringFactorTable = field(default_factory=ScatteringFactorTable)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        if not np.isfinite(self.amplitude) or self.amplitude < 0.0:
            raise ValueError("StructureFactor.amplitude must be finite and non-negative.")
        if not np.isfinite(self.phase_rad):
            raise ValueError("StructureFactor.phase_rad must be finite.")

    @classmethod
    def from_phase_hkl(
        cls,
        phase: Phase,
        hkl: ArrayLike,
        *,
        scattering_table: ScatteringFactorTable | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> StructureFactor:
        indices = as_int_array(hkl, shape=(3,))
        table = ScatteringFactorTable() if scattering_table is None else scattering_table
        if phase.unit_cell is None or not phase.unit_cell.sites:
            value = 1.0 + 0.0j
        else:
            reciprocal = phase.lattice.reciprocal_basis().matrix
            g_cart = reciprocal @ indices.astype(np.float64)
            g_magnitude = float(np.linalg.norm(g_cart))
            total = 0.0j
            for site in phase.unit_cell.sites:
                factor = table.scattering_factor(site.species, g_magnitude)
                occupancy = float(site.occupancy)
                b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
                debye_waller = np.exp(
                    -(b_iso * g_magnitude * g_magnitude) / max(16.0 * np.pi * np.pi, 1e-12)
                )
                phase_argument = float(
                    np.dot(indices.astype(np.float64), site.fractional_coordinates)
                )
                phase_factor = np.exp(
                    2.0j * np.pi * phase_argument
                )
                total += occupancy * factor * debye_waller * phase_factor
            value = complex(total)
        return cls(
            miller_indices=indices,
            value=value,
            amplitude=float(abs(value)),
            phase_rad=float(np.angle(value)),
            phase=phase,
            scattering_table=table,
            provenance=provenance or table.provenance or phase.provenance,
        )


@dataclass(frozen=True, slots=True)
class ReflectionCondition:
    centering: str = "P"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        centering = self.centering.strip().upper()
        if centering not in {"P", "I", "F", "A", "B", "C", "R"}:
            raise ValueError("ReflectionCondition.centering must be P, I, F, A, B, C, or R.")
        object.__setattr__(self, "centering", centering)

    @classmethod
    def from_phase(cls, phase: Phase) -> ReflectionCondition:
        symbol = phase.space_group_symbol
        if symbol is None:
            symbol = phase.space_group.symbol if phase.space_group else "P"
        centering = symbol.strip()[:1].upper() if symbol.strip() else "P"
        return cls(centering=centering, provenance=phase.provenance)

    def is_allowed(self, hkl: ArrayLike) -> bool:
        h, k, ell = (int(value) for value in as_int_array(hkl, shape=(3,)))
        if self.centering == "P":
            return True
        if self.centering == "I":
            return (h + k + ell) % 2 == 0
        if self.centering == "F":
            return h % 2 == k % 2 == ell % 2
        if self.centering == "A":
            return (k + ell) % 2 == 0
        if self.centering == "B":
            return (h + ell) % 2 == 0
        if self.centering == "C":
            return (h + k) % 2 == 0
        return (-h + k + ell) % 3 == 0


@dataclass(frozen=True, slots=True)
class DiffractionIntensityModel:
    scattering_table: ScatteringFactorTable = field(default_factory=ScatteringFactorTable)
    correction: IntensityCorrectionName = "lorentz_polarization"
    reflection_condition: ReflectionCondition | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.correction not in {"none", "lorentz_polarization"}:
            raise ValueError("correction must be either 'none' or 'lorentz_polarization'.")

    def intensity(
        self,
        phase: Phase,
        hkl: ArrayLike,
        *,
        two_theta_rad: float,
        multiplicity: int = 1,
        radiation: RadiationSpec | None = None,
    ) -> float:
        del radiation
        condition = self.reflection_condition or ReflectionCondition.from_phase(phase)
        if not condition.is_allowed(hkl):
            return 0.0
        structure_factor = StructureFactor.from_phase_hkl(
            phase,
            hkl,
            scattering_table=self.scattering_table,
            provenance=self.provenance,
        )
        correction = (
            1.0
            if self.correction == "none"
            else lorentz_polarization_factor(two_theta_rad)
        )
        return float(max(0.0, multiplicity * structure_factor.amplitude**2 * correction))


__all__ = [
    "DiffractionIntensityModel",
    "ReflectionCondition",
    "ScatteringFactorTable",
    "StructureFactor",
    "lorentz_polarization_factor",
]
