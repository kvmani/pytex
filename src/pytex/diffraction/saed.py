from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array, normalize_vector
from pytex.core._chemistry import atomic_number
from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import MillerIndex, Phase, ReciprocalLatticeVector, ZoneAxis
from pytex.core.provenance import ProvenanceRecord


def _choose_zone_basis(zone_axis: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(trial, zone_axis))), 1.0, atol=1e-8):
        trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u_axis = normalize_vector(np.cross(zone_axis, trial))
    v_axis = normalize_vector(np.cross(zone_axis, u_axis))
    basis = np.column_stack([u_axis, v_axis, zone_axis])
    basis = np.ascontiguousarray(basis)
    basis.setflags(write=False)
    return basis


def _structure_factor_electron(phase: Phase, hkl: np.ndarray) -> complex:
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
        damping = np.exp(-(b_iso * g_sq) / max(16.0 * np.pi * np.pi, 1e-12))
        phase_factor = np.exp(
            2.0j * np.pi * float(np.dot(hkl.astype(np.float64), site.fractional_coordinates))
        )
        total += occupancy * z * damping * phase_factor
    return complex(total)


def _enumerate_zone_hkls(max_index: int) -> np.ndarray:
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
class SAEDSpot:
    miller_indices: np.ndarray
    reciprocal_vector_crystal: np.ndarray
    reciprocal_vector_detector: np.ndarray
    detector_coordinates: np.ndarray
    intensity: float
    excitation_error_inv_angstrom: float
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        object.__setattr__(
            self,
            "reciprocal_vector_crystal",
            as_float_array(self.reciprocal_vector_crystal, shape=(3,)),
        )
        object.__setattr__(
            self,
            "reciprocal_vector_detector",
            as_float_array(self.reciprocal_vector_detector, shape=(3,)),
        )
        object.__setattr__(
            self,
            "detector_coordinates",
            as_float_array(self.detector_coordinates, shape=(2,)),
        )
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("SAEDSpot.intensity must be finite and non-negative.")
        if not np.isfinite(self.excitation_error_inv_angstrom):
            raise ValueError("SAEDSpot.excitation_error_inv_angstrom must be finite.")


@dataclass(frozen=True, slots=True)
class SAEDPattern:
    phase: Phase
    zone_axis: ZoneAxis
    detector_frame: ReferenceFrame
    reciprocal_frame: ReferenceFrame
    camera_constant_mm_angstrom: float
    spots: tuple[SAEDSpot, ...]
    zone_basis_crystal: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "spots", tuple(self.spots))
        object.__setattr__(
            self, "zone_basis_crystal", as_float_array(self.zone_basis_crystal, shape=(3, 3))
        )
        if self.detector_frame.domain is not FrameDomain.DETECTOR:
            raise ValueError("SAEDPattern.detector_frame must belong to the detector domain.")
        if self.reciprocal_frame.domain is not FrameDomain.RECIPROCAL:
            raise ValueError("SAEDPattern.reciprocal_frame must belong to the reciprocal domain.")
        if self.camera_constant_mm_angstrom <= 0.0:
            raise ValueError("SAEDPattern.camera_constant_mm_angstrom must be positive.")

    def detector_extent_mm(self) -> float:
        if not self.spots:
            return 1.0
        radii = [float(np.linalg.norm(spot.detector_coordinates)) for spot in self.spots]
        return max(radii) * 1.15


def generate_saed_pattern(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    camera_constant_mm_angstrom: float = 180.0,
    max_index: int = 6,
    max_g_inv_angstrom: float | None = None,
    zone_tolerance_inv_angstrom: float = 1e-6,
    intensity_model: Literal["electron_atomic_number", "unit"] = "electron_atomic_number",
    label_limit: int = 20,
    provenance: ProvenanceRecord | None = None,
) -> SAEDPattern:
    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if camera_constant_mm_angstrom <= 0.0:
        raise ValueError("camera_constant_mm_angstrom must be strictly positive.")
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if max_g_inv_angstrom is not None and max_g_inv_angstrom <= 0.0:
        raise ValueError("max_g_inv_angstrom must be strictly positive when provided.")
    if zone_tolerance_inv_angstrom < 0.0:
        raise ValueError("zone_tolerance_inv_angstrom must be non-negative.")

    zone_vector = zone_axis.unit_vector
    zone_basis = _choose_zone_basis(zone_vector)
    reciprocal_basis = phase.lattice.reciprocal_basis()
    reciprocal_frame = reciprocal_basis.frame
    detector_frame = ReferenceFrame(
        name=f"{phase.name}_saed_detector",
        domain=FrameDomain.DETECTOR,
        axes=("u", "v", "n"),
        handedness=Handedness.RIGHT,
        description="Detector plane for kinematic SAED plotting.",
        provenance=provenance,
    )
    spots: list[SAEDSpot] = []
    for hkl in _enumerate_zone_hkls(max_index):
        reciprocal_vector = ReciprocalLatticeVector.from_miller_index(
            MillerIndex(hkl, phase=phase)
        ).cartesian_vector
        zone_projection = float(np.dot(reciprocal_vector, zone_vector))
        if abs(zone_projection) > zone_tolerance_inv_angstrom:
            continue
        g_magnitude = float(np.linalg.norm(reciprocal_vector))
        if np.isclose(g_magnitude, 0.0):
            continue
        if max_g_inv_angstrom is not None and g_magnitude > max_g_inv_angstrom:
            continue
        detector_coords_mm = camera_constant_mm_angstrom * np.array(
            [
                float(np.dot(reciprocal_vector, zone_basis[:, 0])),
                float(np.dot(reciprocal_vector, zone_basis[:, 1])),
            ],
            dtype=np.float64,
        )
        intensity = (
            1.0
            if intensity_model == "unit"
            else float(abs(_structure_factor_electron(phase, hkl)) ** 2)
            / (1.0 + g_magnitude * g_magnitude)
        )
        spots.append(
            SAEDSpot(
                miller_indices=hkl,
                reciprocal_vector_crystal=reciprocal_vector,
                reciprocal_vector_detector=zone_basis.T @ reciprocal_vector,
                detector_coordinates=detector_coords_mm,
                intensity=intensity,
                excitation_error_inv_angstrom=zone_projection,
                label=" ".join(str(int(value)) for value in hkl),
            )
        )
    spots.sort(
        key=lambda spot: (
            -spot.intensity,
            float(np.linalg.norm(spot.detector_coordinates)),
        )
    )
    limited_spots = []
    for index, spot in enumerate(spots):
        if index >= label_limit:
            limited_spots.append(
                SAEDSpot(
                    miller_indices=spot.miller_indices,
                    reciprocal_vector_crystal=spot.reciprocal_vector_crystal,
                    reciprocal_vector_detector=spot.reciprocal_vector_detector,
                    detector_coordinates=spot.detector_coordinates,
                    intensity=spot.intensity,
                    excitation_error_inv_angstrom=spot.excitation_error_inv_angstrom,
                    label=None,
                )
            )
        else:
            limited_spots.append(spot)
    return SAEDPattern(
        phase=phase,
        zone_axis=zone_axis,
        detector_frame=detector_frame,
        reciprocal_frame=reciprocal_frame,
        camera_constant_mm_angstrom=camera_constant_mm_angstrom,
        spots=tuple(limited_spots),
        zone_basis_crystal=zone_basis,
        provenance=provenance,
    )
