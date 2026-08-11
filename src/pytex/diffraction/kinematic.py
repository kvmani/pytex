"""Vectorized kinematic zone-axis SAED simulation.

Purpose: the batch-oriented kinematic diffraction engine behind composite
orientation-relationship (OR) pattern simulation. Unlike the legacy
loop-based :func:`pytex.diffraction.saed.generate_saed_pattern` (kept for
backward compatibility), this module:

- vectorizes every per-reflection computation (enumeration, systematic
  absences, structure factors, excitation errors, detector projection);
- selects reflections by **excitation error** rather than the integer zone
  law, so irrational zone axes (the generic case for OR-mapped child zones)
  are handled exactly;
- exposes the detector in-plane basis as an explicit, configurable object so
  several phases/variants can share one physically consistent detector frame.

Conventions (pinned; see
``docs/roadmap/working_notes_composite_saed_program.md``):

- The beam travels antiparallel to the zone axis; the zone-axis unit vector
  points toward the electron gun.
- The excitation error of reflection ``g`` is ``s_g = g_z - g^2 * lambda / 2``
  where ``g_z`` is the zone-axis component of ``g`` (small-angle kinematic
  approximation; ``s_g = -g^2 lambda / 2 <= 0`` for exact zero-order-Laue-zone
  reflections, reflecting Ewald-sphere curvature).
- Detector coordinates follow the standard SAED small-angle relation
  ``r_mm = (camera constant) * g_perp`` with the camera constant in mm*angstrom.
- Kinematic intensity is ``|F_hkl|^2`` from the electron structure-factor
  proxy (atomic-number scattering with isotropic Debye-Waller damping),
  optionally damped by a relrod profile ``1 / (1 + (s_g / sigma_s)^2)``;
  intensities are normalized to a maximum of 1 per pattern.
- Double diffraction is off by default. When enabled, reflections reachable as
  the integer sum of two excited reflections are added even where the structure
  factor forbids them, and are flagged as such; see
  :func:`double_diffraction_sums`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import constants

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import (
    CrystalDirection,
    Phase,
    ZoneAxis,
    phases_semantically_match,
)
from pytex.core.notation import format_direction_indices
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.physics import ReflectionCondition

IntensityModelName = Literal["electron_atomic_number", "unit"]

_BASIS_ORTHONORMALITY_ATOL = 1e-9
_BASIS_ZONE_ALIGNMENT_ATOL = 1e-9
_DEBYE_WALLER_DENOMINATOR = 16.0 * np.pi * np.pi


def electron_wavelength_angstrom(beam_energy_kev: float) -> float:
    """Relativistic electron wavelength for a given accelerating voltage.

    Purpose: convert the microscope accelerating voltage into the electron
    wavelength that fixes the Ewald-sphere radius ``k = 1/lambda``.

    When to use: any kinematic TEM computation; the standard values are
    0.037014 angstrom at 100 kV, 0.025079 angstrom at 200 kV and
    0.019687 angstrom at 300 kV (De Graef, Introduction to Conventional
    Transmission Electron Microscopy, 2003, Table 2.2).

    Inputs: ``beam_energy_kev`` — accelerating voltage in kilovolts,
    strictly positive.

    Output: wavelength in angstrom from
    ``lambda = h / sqrt(2 m0 e V (1 + e V / (2 m0 c^2)))``.
    """

    if not np.isfinite(beam_energy_kev) or beam_energy_kev <= 0.0:
        raise ValueError("beam_energy_kev must be finite and strictly positive.")
    voltage = float(beam_energy_kev) * 1e3
    rest_energy = constants.electron_mass * constants.c**2
    momentum_term = (
        2.0
        * constants.electron_mass
        * constants.elementary_charge
        * voltage
        * (1.0 + constants.elementary_charge * voltage / (2.0 * rest_energy))
    )
    wavelength_m = constants.h / np.sqrt(momentum_term)
    return float(wavelength_m * 1e10)


@dataclass(frozen=True, slots=True)
class KinematicSimulationConfig:
    """Configuration for vectorized kinematic zone-axis simulation.

    ``max_excitation_error_inv_angstrom`` is the reflection-selection
    half-width: a reflection is retained when ``|s_g|`` does not exceed it.
    The default (0.05 1/angstrom) keeps every zero-order-Laue-zone reflection
    within the default ``g`` range at common TEM voltages while excluding
    higher-order Laue zones. ``relrod_sigma_inv_angstrom`` optionally damps
    intensity with distance from the exact Bragg condition;
    ``None`` disables damping. ``min_relative_intensity`` removes
    numerically-forbidden reflections after normalization.

    ``include_double_diffraction`` adds the reflections that are kinematically
    forbidden but observable in a real pattern because a diffracted beam
    re-diffracts: a spot appears at ``g1 + g2`` for any two excited
    reflections. ``double_diffraction_coupling`` scales their indicative
    intensity; see :func:`simulate_zone_axis_spots`.
    """

    beam_energy_kev: float = 200.0
    camera_constant_mm_angstrom: float = 180.0
    max_index: int = 6
    g_max_inv_angstrom: float | None = 1.5
    max_excitation_error_inv_angstrom: float = 0.05
    intensity_model: IntensityModelName = "electron_atomic_number"
    relrod_sigma_inv_angstrom: float | None = None
    apply_centering_absences: bool = True
    min_relative_intensity: float = 1e-4
    include_double_diffraction: bool = False
    double_diffraction_coupling: float = 0.05

    def __post_init__(self) -> None:
        if not np.isfinite(self.beam_energy_kev) or self.beam_energy_kev <= 0.0:
            raise ValueError("beam_energy_kev must be finite and strictly positive.")
        if (
            not np.isfinite(self.camera_constant_mm_angstrom)
            or self.camera_constant_mm_angstrom <= 0.0
        ):
            raise ValueError("camera_constant_mm_angstrom must be finite and strictly positive.")
        if int(self.max_index) != self.max_index or self.max_index <= 0:
            raise ValueError("max_index must be a strictly positive integer.")
        if self.g_max_inv_angstrom is not None and (
            not np.isfinite(self.g_max_inv_angstrom) or self.g_max_inv_angstrom <= 0.0
        ):
            raise ValueError("g_max_inv_angstrom must be finite and strictly positive when set.")
        if (
            not np.isfinite(self.max_excitation_error_inv_angstrom)
            or self.max_excitation_error_inv_angstrom <= 0.0
        ):
            raise ValueError(
                "max_excitation_error_inv_angstrom must be finite and strictly positive."
            )
        if self.intensity_model not in {"electron_atomic_number", "unit"}:
            raise ValueError(
                "intensity_model must be 'electron_atomic_number' or 'unit'."
            )
        if self.relrod_sigma_inv_angstrom is not None and (
            not np.isfinite(self.relrod_sigma_inv_angstrom)
            or self.relrod_sigma_inv_angstrom <= 0.0
        ):
            raise ValueError(
                "relrod_sigma_inv_angstrom must be finite and strictly positive when set."
            )
        if (
            not np.isfinite(self.min_relative_intensity)
            or not 0.0 <= self.min_relative_intensity < 1.0
        ):
            raise ValueError("min_relative_intensity must lie in the interval [0, 1).")
        if (
            not np.isfinite(self.double_diffraction_coupling)
            or not 0.0 < self.double_diffraction_coupling <= 1.0
        ):
            raise ValueError("double_diffraction_coupling must lie in the interval (0, 1].")

    @property
    def wavelength_angstrom(self) -> float:
        return electron_wavelength_angstrom(self.beam_energy_kev)


def zone_basis_from_axis(
    zone_cartesian: np.ndarray,
    *,
    align_g_cartesian: np.ndarray | None = None,
    in_plane_rotation_deg: float = 0.0,
) -> np.ndarray:
    """Deterministic right-handed detector basis for a zone axis.

    Purpose: build the orthonormal ``(u, v, z)`` triad whose first two
    columns span the detector plane. Column 3 is the zone-axis unit vector.

    Inputs: ``zone_cartesian`` — nonzero zone-axis vector (crystal Cartesian
    frame). ``align_g_cartesian`` — optional reciprocal-space vector whose
    in-plane component is aligned with ``+u`` (it must not be parallel to
    the zone axis). ``in_plane_rotation_deg`` — rotates the rendered pattern
    counterclockwise by this angle (right-handed about the zone axis, which
    points toward the viewer).

    Output: read-only ``(3, 3)`` array with columns ``u``, ``v``, ``z``
    satisfying ``u x v = z``. Without ``align_g_cartesian`` the construction
    matches the legacy ``generate_saed_pattern`` basis for reproducibility.
    """

    zone_unit = normalize_vector(zone_cartesian)
    if align_g_cartesian is None:
        trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if np.isclose(abs(float(np.dot(trial, zone_unit))), 1.0, atol=1e-8):
            trial = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        u_axis = normalize_vector(np.cross(zone_unit, trial))
        v_axis = normalize_vector(np.cross(zone_unit, u_axis))
    else:
        align = as_float_array(align_g_cartesian, shape=(3,))
        in_plane = align - float(np.dot(align, zone_unit)) * zone_unit
        norm = float(np.linalg.norm(in_plane))
        if np.isclose(norm, 0.0, atol=1e-12):
            raise ValueError(
                "align_g_cartesian must not be parallel to the zone axis."
            )
        u_axis = normalize_vector(in_plane)
        v_axis = np.cross(zone_unit, u_axis)
    angle_rad = float(np.deg2rad(in_plane_rotation_deg))
    if not np.isfinite(angle_rad):
        raise ValueError("in_plane_rotation_deg must be finite.")
    cos_angle = float(np.cos(angle_rad))
    sin_angle = float(np.sin(angle_rad))
    u_rotated = cos_angle * u_axis - sin_angle * v_axis
    v_rotated = sin_angle * u_axis + cos_angle * v_axis
    basis = np.column_stack([u_rotated, v_rotated, zone_unit])
    basis = np.ascontiguousarray(basis)
    basis.setflags(write=False)
    return basis


def _validate_zone_basis(basis: np.ndarray, zone_unit: np.ndarray) -> np.ndarray:
    array = as_float_array(basis, shape=(3, 3))
    gram = array.T @ array
    if not np.allclose(gram, np.eye(3), atol=_BASIS_ORTHONORMALITY_ATOL):
        raise ValueError("basis columns must be orthonormal.")
    if not np.allclose(
        np.cross(array[:, 0], array[:, 1]), array[:, 2], atol=_BASIS_ORTHONORMALITY_ATOL
    ):
        raise ValueError("basis must be right-handed with u x v = z.")
    if not np.allclose(array[:, 2], zone_unit, atol=_BASIS_ZONE_ALIGNMENT_ATOL):
        raise ValueError("basis column 3 must equal the zone-axis unit vector.")
    return array


def centering_allowed_mask(hkl: np.ndarray, condition: ReflectionCondition) -> np.ndarray:
    """Vectorized lattice-centering systematic-absence mask.

    Purpose: batch equivalent of
    :meth:`pytex.diffraction.physics.ReflectionCondition.is_allowed`.

    Inputs: ``hkl`` — integer array of shape ``(N, 3)``; ``condition`` — the
    centering rule (P, I, F, A, B, C or R).

    Output: boolean array of shape ``(N,)``, ``True`` where the reflection is
    allowed by the centering translation group.
    """

    indices = np.asarray(hkl, dtype=np.int64)
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (N, 3).")
    h, k, ell = indices[:, 0], indices[:, 1], indices[:, 2]
    centering = condition.centering
    if centering == "P":
        return np.ones(indices.shape[0], dtype=bool)
    if centering == "I":
        return np.asarray((h + k + ell) % 2 == 0)
    if centering == "F":
        parity_h = h % 2
        return np.asarray((parity_h == k % 2) & (parity_h == ell % 2))
    if centering == "A":
        return np.asarray((k + ell) % 2 == 0)
    if centering == "B":
        return np.asarray((h + ell) % 2 == 0)
    if centering == "C":
        return np.asarray((h + k) % 2 == 0)
    return np.asarray((-h + k + ell) % 3 == 0)


def electron_structure_factors(
    phase: Phase, hkl: np.ndarray, g_magnitude_inv_angstrom: np.ndarray
) -> np.ndarray:
    """Vectorized electron structure factors for a batch of reflections.

    Purpose: batch form of the atomic-number electron scattering proxy used
    across the kinematic SAED surface — per site
    ``occupancy * Z * exp(-B_iso g^2 / 16 pi^2) * exp(2 pi i (hkl . x))``,
    summed over the unit cell.

    Inputs: ``hkl`` shape ``(N, 3)``; ``g_magnitude_inv_angstrom`` shape
    ``(N,)`` — the Cartesian reciprocal-vector magnitudes for each row.

    Output: complex array of shape ``(N,)``. A phase without unit-cell sites
    yields unit factors (geometry-only pattern).
    """

    indices = np.asarray(hkl, dtype=np.float64)
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (N, 3).")
    g_squared = np.asarray(g_magnitude_inv_angstrom, dtype=np.float64) ** 2
    if g_squared.shape != (indices.shape[0],):
        raise ValueError("g_magnitude_inv_angstrom must have shape (N,).")
    if phase.unit_cell is None or not phase.unit_cell.sites:
        return np.ones(indices.shape[0], dtype=np.complex128)
    total = np.zeros(indices.shape[0], dtype=np.complex128)
    for site in phase.unit_cell.sites:
        z_value = float(atomic_number(site.species))
        occupancy = float(site.occupancy)
        b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
        damping = np.exp(-(b_iso * g_squared) / _DEBYE_WALLER_DENOMINATOR)
        phase_argument = indices @ site.fractional_coordinates
        total += occupancy * z_value * damping * np.exp(2.0j * np.pi * phase_argument)
    return total


#: Detector-radius quantum used when ordering spots, in decimal places of a
#: millimetre. Nine places is 1 pm — orders of magnitude below any physically
#: meaningful detector distance, and orders of magnitude above the ~1e-14 mm
#: floating-point spread between symmetry-equivalent reflections.
_RADIUS_SORT_DECIMALS = 9

#: Relative-intensity quantum used when ordering spots. Intensities are
#: normalized to a maximum of 1, so twelve places is 1e-12 of full scale.
_INTENSITY_SORT_DECIMALS = 12


def _hkl_keys(hkl: np.ndarray, max_index: int) -> np.ndarray:
    """Injective int64 encoding of indices bounded by ``max_index``."""

    span = 2 * int(max_index) + 1
    shifted = np.asarray(hkl, dtype=np.int64) + int(max_index)
    return np.asarray(
        shifted[:, 0] * span * span + shifted[:, 1] * span + shifted[:, 2], dtype=np.int64
    )


def double_diffraction_sums(
    hkl: np.ndarray, intensity: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pairwise integer sums of excited reflections — the double-diffraction rule.

    Purpose: enumerate the reflections a doubly-diffracted beam can reach. A
    beam diffracted by ``g1`` is itself an incident beam inside the crystal, so
    diffracting it again by ``g2`` sends it out along ``g1 + g2``. The set of
    reachable spots is therefore the set of pairwise algebraic sums of the
    excited reflections, which is a geometric selection rule and does not
    depend on the dynamical theory that supplies the intensities.

    When to use: to decide which kinematically forbidden reflections may still
    appear in a zone-axis pattern. Note that lattice-centring conditions define
    a *sublattice* of reciprocal space and a sublattice is closed under
    addition, so this rule can never revive a centring absence — only a basis
    absence from a glide plane, a screw axis, or the motif.

    Inputs: ``hkl`` shape ``(M, 3)`` integer reflections, and ``intensity``
    shape ``(M,)`` their non-negative relative intensities. Each unordered pair
    is taken once, including a reflection with itself (the ``2g`` path).

    Output: ``(sums, weight, parents)`` — unique sums of shape ``(K, 3)``, the
    summed two-step coupling weight ``sum over paths of I1 * I2`` of shape
    ``(K,)``, and the strongest contributing parent pair of shape ``(K, 2, 3)``.
    The zero sum from a ``g``/``-g`` pair, which is the transmitted beam, is
    excluded.
    """

    indices = np.asarray(hkl, dtype=np.int64)
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (M, 3).")
    weights_in = np.asarray(intensity, dtype=np.float64)
    if weights_in.shape != (indices.shape[0],):
        raise ValueError("intensity must have shape (M,).")
    if indices.shape[0] == 0:
        return (
            np.zeros((0, 3), dtype=np.int64),
            np.zeros((0,), dtype=np.float64),
            np.zeros((0, 2, 3), dtype=np.int64),
        )

    left, right = np.triu_indices(indices.shape[0])
    sums = indices[left] + indices[right]
    weights = weights_in[left] * weights_in[right]
    nonzero = np.any(sums != 0, axis=1)
    sums, weights, left, right = sums[nonzero], weights[nonzero], left[nonzero], right[nonzero]

    unique, inverse = np.unique(sums, axis=0, return_inverse=True)
    inverse = np.asarray(inverse, dtype=np.int64).ravel()
    total = np.bincount(inverse, weights=weights, minlength=unique.shape[0])
    # Last write wins over ascending weight, so each group keeps its strongest
    # path — a stable way to pick a per-group argmax without a Python loop.
    order = np.argsort(weights, kind="stable")
    best = np.zeros(unique.shape[0], dtype=np.int64)
    best[inverse[order]] = order
    parents = np.stack([indices[left[best]], indices[right[best]]], axis=1)
    return unique, np.asarray(total, dtype=np.float64), parents


def _enumerate_hkl_cube(max_index: int) -> np.ndarray:
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid_h, grid_k, grid_l = np.meshgrid(values, values, values, indexing="ij")
    hkl = np.column_stack([grid_h.ravel(), grid_k.ravel(), grid_l.ravel()])
    return np.asarray(hkl[np.any(hkl != 0, axis=1)], dtype=np.int64)


@dataclass(frozen=True, slots=True)
class SpotTable:
    """Struct-of-arrays result of a vectorized zone-axis simulation.

    All arrays share the leading reflection dimension ``N`` and are
    read-only; rows are sorted by decreasing intensity, then increasing
    detector radius, then lexicographic ``hkl`` for full determinism. The two
    continuous keys are quantized before sorting (1 pm of radius, 1e-12 of
    full-scale intensity) so that symmetry-equivalent reflections, whose values
    are mathematically equal but differ in the last few ULPs, are ordered by
    their exact indices rather than by floating-point noise.
    ``basis`` columns are the detector ``u``, ``v`` axes and the zone-axis
    unit vector, all in the phase crystal Cartesian frame. ``detector_mm``
    stacks the ``(u, v)`` detector-plane coordinates in millimetres;
    ``g_detector_inv_angstrom`` holds the same coordinates in reciprocal
    angstrom before camera scaling.

    ``is_double_diffraction`` marks the rows that are **kinematically
    forbidden** and present only because the double-diffraction option was
    enabled; ``double_diffraction_parents`` gives the strongest contributing
    reflection pair ``g1 + g2`` for those rows and is all zeros elsewhere.
    Treat a marked row as an observability statement, not as a kinematic
    intensity: its structure factor is (near) zero, which is why it is marked.
    """

    phase: Phase
    zone_axis: ZoneAxis | CrystalDirection
    basis: np.ndarray
    config: KinematicSimulationConfig
    hkl: np.ndarray
    g_crystal: np.ndarray
    g_detector_inv_angstrom: np.ndarray
    detector_mm: np.ndarray
    d_spacing_angstrom: np.ndarray
    structure_factor_amplitude: np.ndarray
    intensity: np.ndarray
    excitation_error_inv_angstrom: np.ndarray
    provenance: ProvenanceRecord | None = None
    is_double_diffraction: np.ndarray | None = None
    double_diffraction_parents: np.ndarray | None = None

    def __post_init__(self) -> None:
        basis = as_float_array(self.basis, shape=(3, 3))
        hkl = np.ascontiguousarray(np.asarray(self.hkl, dtype=np.int64))
        if hkl.ndim != 2 or hkl.shape[1] != 3:
            raise ValueError("SpotTable.hkl must have shape (N, 3).")
        count = hkl.shape[0]
        g_crystal = as_float_array(self.g_crystal, shape=(count, 3) if count else (0, 3))
        g_detector = as_float_array(
            self.g_detector_inv_angstrom, shape=(count, 2) if count else (0, 2)
        )
        detector_mm = as_float_array(self.detector_mm, shape=(count, 2) if count else (0, 2))
        d_spacing = as_float_array(self.d_spacing_angstrom, shape=(count,))
        amplitude = as_float_array(self.structure_factor_amplitude, shape=(count,))
        intensity = as_float_array(self.intensity, shape=(count,))
        excitation = as_float_array(self.excitation_error_inv_angstrom, shape=(count,))
        if self.is_double_diffraction is None:
            forbidden = np.zeros(count, dtype=bool)
        else:
            forbidden = np.ascontiguousarray(np.asarray(self.is_double_diffraction, dtype=bool))
            if forbidden.shape != (count,):
                raise ValueError("SpotTable.is_double_diffraction must have shape (N,).")
        if self.double_diffraction_parents is None:
            parents = np.zeros((count, 2, 3), dtype=np.int64)
        else:
            parents = np.ascontiguousarray(
                np.asarray(self.double_diffraction_parents, dtype=np.int64)
            )
            if parents.shape != (count, 2, 3):
                raise ValueError(
                    "SpotTable.double_diffraction_parents must have shape (N, 2, 3)."
                )
            if count and np.any(parents.sum(axis=1)[forbidden] != hkl[forbidden]):
                raise ValueError(
                    "SpotTable.double_diffraction_parents must sum to the reflection they "
                    "produce."
                )
        if count:
            if np.any(~np.isfinite(g_crystal)) or np.any(~np.isfinite(detector_mm)):
                raise ValueError("SpotTable coordinate arrays must be finite.")
            if np.any(d_spacing <= 0.0) or np.any(~np.isfinite(d_spacing)):
                raise ValueError("SpotTable.d_spacing_angstrom must be finite and positive.")
            if np.any(intensity < 0.0) or np.any(~np.isfinite(intensity)):
                raise ValueError("SpotTable.intensity must be finite and non-negative.")
            if np.any(amplitude < 0.0) or np.any(~np.isfinite(amplitude)):
                raise ValueError(
                    "SpotTable.structure_factor_amplitude must be finite and non-negative."
                )
            if np.any(~np.isfinite(excitation)):
                raise ValueError("SpotTable.excitation_error_inv_angstrom must be finite.")
        for array in (
            basis,
            hkl,
            g_crystal,
            g_detector,
            detector_mm,
            d_spacing,
            amplitude,
            intensity,
            excitation,
            forbidden,
            parents,
        ):
            array.setflags(write=False)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "hkl", hkl)
        object.__setattr__(self, "g_crystal", g_crystal)
        object.__setattr__(self, "g_detector_inv_angstrom", g_detector)
        object.__setattr__(self, "detector_mm", detector_mm)
        object.__setattr__(self, "d_spacing_angstrom", d_spacing)
        object.__setattr__(self, "structure_factor_amplitude", amplitude)
        object.__setattr__(self, "intensity", intensity)
        object.__setattr__(self, "excitation_error_inv_angstrom", excitation)
        object.__setattr__(self, "is_double_diffraction", forbidden)
        object.__setattr__(self, "double_diffraction_parents", parents)

    def __len__(self) -> int:
        return int(self.hkl.shape[0])

    @property
    def wavelength_angstrom(self) -> float:
        return self.config.wavelength_angstrom

    def zone_axis_label(self) -> str:
        if isinstance(self.zone_axis, ZoneAxis):
            return format_direction_indices(
                tuple(int(value) for value in self.zone_axis.indices), style="plain"
            )
        # An irrational zone axis has no integer indices to bracket, so the
        # decimal components are shown inside direction brackets directly.
        coordinates = self.zone_axis.coordinates
        return "[" + " ".join(f"{value:.3f}" for value in coordinates) + "]"

    def hkl_labels(self) -> tuple[str, ...]:
        return tuple(
            " ".join(str(int(value)) for value in row) for row in self.hkl
        )

    def detector_radius_mm(self) -> np.ndarray:
        return np.asarray(np.linalg.norm(self.detector_mm, axis=1), dtype=np.float64)

    def forbidden_mask(self) -> np.ndarray:
        """Rows that are kinematically forbidden and present only by double diffraction."""

        return np.asarray(self.is_double_diffraction, dtype=bool)

    def double_diffraction_origin_label(self, row: int) -> str:
        """How row ``row`` arises by double diffraction, as ``g1 + g2 = g``.

        Returns the empty string for a row that is a genuine kinematic
        reflection, so a caller can annotate unconditionally.
        """

        index = int(row)
        if not bool(self.forbidden_mask()[index]):
            return ""
        first, second = np.asarray(self.double_diffraction_parents)[index]
        return (
            f"({' '.join(str(int(v)) for v in first)}) + "
            f"({' '.join(str(int(v)) for v in second)}) = "
            f"({' '.join(str(int(v)) for v in self.hkl[index])})"
        )

    def describe(self) -> str:
        """Prose summary of the simulated zone-axis pattern and conventions."""

        config = self.config
        header = (
            f"Kinematic zone-axis SAED pattern for phase '{self.phase.name}' along "
            f"{self.zone_axis_label()} (crystal-frame zone axis pointing toward the gun; "
            f"beam antiparallel): {len(self)} reflection(s) at "
            f"{config.beam_energy_kev:g} kV (lambda = {self.wavelength_angstrom:.6f} "
            f"angstrom), camera constant {config.camera_constant_mm_angstrom:g} mm*angstrom."
        )
        selection = (
            "Reflections satisfy |s_g| <= "
            f"{config.max_excitation_error_inv_angstrom:g} 1/angstrom with "
            "s_g = g_z - g^2 lambda / 2 (excitation-error selection; exact "
            "zero-order-Laue-zone spots have s_g = -g^2 lambda / 2)."
        )
        intensity = (
            "Intensities are |F_hkl|^2 from the atomic-number electron scattering proxy, "
            "normalized to max 1."
            if config.intensity_model == "electron_atomic_number"
            else "Intensities use the geometry-only unit model."
        )
        if len(self) == 0:
            return f"{header} {selection} {intensity} The pattern is empty."
        strongest = ", ".join(
            f"({label})" for label in self.hkl_labels()[: min(4, len(self))]
        )
        summary = (
            f"{header} {selection} {intensity} Strongest reflections: {strongest}. "
            f"Minimum d-spacing {float(np.min(self.d_spacing_angstrom)):.4f} angstrom."
        )
        forbidden_rows = np.flatnonzero(self.forbidden_mask())
        if not config.include_double_diffraction:
            return (
                f"{summary} Double diffraction is disabled, so kinematically forbidden "
                "reflections that a real pattern would show are absent."
            )
        if forbidden_rows.size == 0:
            return (
                f"{summary} Double diffraction is enabled but no forbidden reflection is "
                "reachable as a sum of two excited reflections here."
            )
        examples = "; ".join(
            self.double_diffraction_origin_label(int(row)) for row in forbidden_rows[:3]
        )
        return (
            f"{summary} {forbidden_rows.size} reflection(s) are kinematically forbidden and "
            "present only through double diffraction, at an indicative intensity scaled by a "
            f"coupling of {config.double_diffraction_coupling:g}: {examples}."
        )


def _append_double_diffraction_rows(
    phase: Phase,
    *,
    hkl: np.ndarray,
    g_crystal: np.ndarray,
    g_magnitude: np.ndarray,
    excitation_error: np.ndarray,
    amplitude: np.ndarray,
    intensity: np.ndarray,
    reciprocal_matrix: np.ndarray,
    zone_unit: np.ndarray,
    wavelength: float,
    normalization: float,
    config: KinematicSimulationConfig,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Extend a kinematic reflection set with its double-diffraction spots."""

    max_index = int(config.max_index)
    candidates, weight, parents = double_diffraction_sums(hkl, intensity)

    # Sums outside the enumerated cube were never candidates for the kinematic
    # pass either, so admitting them would smuggle in ordinary allowed
    # reflections under the double-diffraction label.
    keep = np.all(np.abs(candidates) <= max_index, axis=1)
    # A sum that already has a spot is that spot; double diffraction adds an
    # unresolvable contribution to it, not a new reflection.
    keep &= ~np.isin(_hkl_keys(candidates, max_index), _hkl_keys(hkl, max_index))
    candidates, weight, parents = candidates[keep], weight[keep], parents[keep]

    g_new = candidates.astype(np.float64) @ reciprocal_matrix.T
    g_magnitude_new = np.linalg.norm(g_new, axis=1)
    keep = g_magnitude_new > 0.0
    if config.g_max_inv_angstrom is not None:
        keep &= g_magnitude_new <= config.g_max_inv_angstrom
    excitation_new = g_new @ zone_unit - 0.5 * wavelength * g_magnitude_new**2
    keep &= np.abs(excitation_new) <= config.max_excitation_error_inv_angstrom
    candidates, weight, parents = candidates[keep], weight[keep], parents[keep]
    g_new = g_new[keep]
    g_magnitude_new, excitation_new = g_magnitude_new[keep], excitation_new[keep]

    amplitude_new = np.abs(electron_structure_factors(phase, candidates, g_magnitude_new))
    if config.intensity_model == "unit":
        base = np.ones_like(amplitude_new)
    else:
        base = amplitude_new**2
    if config.relrod_sigma_inv_angstrom is not None:
        base = base / (1.0 + (excitation_new / config.relrod_sigma_inv_angstrom) ** 2)
    residual = base / normalization if normalization > 0.0 else np.zeros_like(base)
    # The two-step amplitude scales as F(g1) F(g2), so the intensity scales as
    # the product of the two parent intensities; paths are summed because their
    # relative phases are outside what a kinematic treatment can supply. The
    # coupling constant carries everything the kinematic theory cannot — beam
    # coupling strength and specimen thickness — so the result is an
    # observability estimate, not a measured intensity. Clipping keeps the
    # per-pattern "normalized to max 1" invariant intact.
    intensity_new = np.minimum(
        1.0, config.double_diffraction_coupling * weight + residual
    )
    keep = intensity_new >= config.min_relative_intensity
    candidates, parents = candidates[keep], parents[keep]
    g_new = g_new[keep]
    g_magnitude_new, excitation_new = g_magnitude_new[keep], excitation_new[keep]
    amplitude_new, intensity_new = amplitude_new[keep], intensity_new[keep]

    count = hkl.shape[0]
    return (
        np.concatenate([hkl, candidates]),
        np.concatenate([g_crystal, g_new]),
        np.concatenate([g_magnitude, g_magnitude_new]),
        np.concatenate([excitation_error, excitation_new]),
        np.concatenate([amplitude, amplitude_new]),
        np.concatenate([intensity, intensity_new]),
        np.concatenate([np.zeros(count, dtype=bool), np.ones(candidates.shape[0], dtype=bool)]),
        np.concatenate([np.zeros((count, 2, 3), dtype=np.int64), parents]),
    )


def simulate_zone_axis_spots(
    phase: Phase,
    zone_axis: ZoneAxis | CrystalDirection,
    *,
    config: KinematicSimulationConfig | None = None,
    basis: np.ndarray | None = None,
    provenance: ProvenanceRecord | None = None,
) -> SpotTable:
    """Simulate a kinematic zone-axis SAED pattern, fully vectorized.

    Purpose: the core engine for single-phase and composite OR diffraction.
    Enumerates the reflection cube ``|h|,|k|,|l| <= max_index``, applies
    lattice-centering absences, selects reflections by excitation error
    (``|s_g| <= max_excitation_error_inv_angstrom``), computes vectorized
    electron structure factors and projects onto the detector basis.

    When to use: whenever a zone-axis spot pattern is needed — directly for a
    single phase, or via :mod:`pytex.diffraction.composite` for parent+variant
    composites. Irrational zone axes (``CrystalDirection`` with non-integer
    coordinates) are supported exactly.

    Inputs: ``zone_axis`` — a rational :class:`ZoneAxis` or a possibly
    irrational :class:`CrystalDirection` in the same phase; ``basis`` —
    optional explicit detector triad (columns ``u``, ``v``, ``z``; must be
    orthonormal, right-handed, with ``z`` equal to the zone-axis unit
    vector), e.g. a shared parent-anchored basis rotated into this phase's
    crystal frame; defaults to :func:`zone_basis_from_axis`.

    Output: a :class:`SpotTable` sorted by decreasing intensity with
    read-only struct-of-arrays reflection data.

    Double diffraction: with ``config.include_double_diffraction`` the table
    also carries the kinematically forbidden reflections a real pattern shows
    because a diffracted beam re-diffracts, selected by the integer sum rule of
    :func:`double_diffraction_sums` and flagged in
    ``SpotTable.is_double_diffraction``. Their intensity is
    ``coupling * sum over paths of I1 * I2``, an observability estimate whose
    scale is set by the ``double_diffraction_coupling`` constant standing in for
    the beam coupling and specimen thickness a kinematic theory cannot supply.
    """

    simulation_config = KinematicSimulationConfig() if config is None else config
    if not phases_semantically_match(zone_axis.phase, phase):
        raise ValueError("zone_axis.phase must match phase.")
    zone_unit = zone_axis.unit_vector
    if basis is None:
        zone_basis = zone_basis_from_axis(zone_unit)
    else:
        zone_basis = _validate_zone_basis(basis, zone_unit)

    hkl = _enumerate_hkl_cube(simulation_config.max_index)
    if simulation_config.apply_centering_absences:
        condition = ReflectionCondition.from_phase(phase)
        hkl = hkl[centering_allowed_mask(hkl, condition)]

    reciprocal_matrix = phase.lattice.reciprocal_basis().matrix
    g_crystal = hkl.astype(np.float64) @ reciprocal_matrix.T
    g_magnitude = np.linalg.norm(g_crystal, axis=1)
    keep = g_magnitude > 0.0
    if simulation_config.g_max_inv_angstrom is not None:
        keep &= g_magnitude <= simulation_config.g_max_inv_angstrom
    hkl, g_crystal, g_magnitude = hkl[keep], g_crystal[keep], g_magnitude[keep]

    wavelength = simulation_config.wavelength_angstrom
    g_along_zone = g_crystal @ zone_basis[:, 2]
    excitation_error = g_along_zone - 0.5 * wavelength * g_magnitude**2
    keep = np.abs(excitation_error) <= simulation_config.max_excitation_error_inv_angstrom
    hkl, g_crystal, g_magnitude = hkl[keep], g_crystal[keep], g_magnitude[keep]
    excitation_error = excitation_error[keep]

    structure_factors = electron_structure_factors(phase, hkl, g_magnitude)
    amplitude = np.abs(structure_factors)
    if simulation_config.intensity_model == "unit":
        intensity = np.ones_like(amplitude)
    else:
        intensity = amplitude**2
    if simulation_config.relrod_sigma_inv_angstrom is not None:
        ratio = excitation_error / simulation_config.relrod_sigma_inv_angstrom
        intensity = intensity / (1.0 + ratio**2)
    maximum = float(np.max(intensity)) if intensity.size else 0.0
    if maximum > 0.0:
        intensity = intensity / maximum
    keep = intensity >= simulation_config.min_relative_intensity
    hkl, g_crystal, g_magnitude = hkl[keep], g_crystal[keep], g_magnitude[keep]
    excitation_error, amplitude, intensity = (
        excitation_error[keep],
        amplitude[keep],
        intensity[keep],
    )
    is_double_diffraction = np.zeros(hkl.shape[0], dtype=bool)
    parents = np.zeros((hkl.shape[0], 2, 3), dtype=np.int64)

    if simulation_config.include_double_diffraction:
        (
            hkl,
            g_crystal,
            g_magnitude,
            excitation_error,
            amplitude,
            intensity,
            is_double_diffraction,
            parents,
        ) = _append_double_diffraction_rows(
            phase,
            hkl=hkl,
            g_crystal=g_crystal,
            g_magnitude=g_magnitude,
            excitation_error=excitation_error,
            amplitude=amplitude,
            intensity=intensity,
            reciprocal_matrix=reciprocal_matrix,
            zone_unit=zone_basis[:, 2],
            wavelength=wavelength,
            normalization=maximum,
            config=simulation_config,
        )

    g_detector = g_crystal @ zone_basis[:, :2]
    detector_mm = simulation_config.camera_constant_mm_angstrom * g_detector
    with np.errstate(divide="ignore"):
        d_spacing = np.where(g_magnitude > 0.0, 1.0 / g_magnitude, np.inf)
    radius = np.linalg.norm(detector_mm, axis=1)
    # Quantize the two continuous sort keys before ordering. Intensity and
    # radius are ranking keys, not measurements, and reflections related by
    # symmetry have mathematically equal values that differ in the last few
    # ULPs depending on how the detector basis was built. Comparing them raw
    # lets that noise decide the order, so the same pattern reached by two
    # equivalent routes — a parent zone axis, or the same axis recovered from a
    # child variant's zone — can come out permuted. Rounding first sends every
    # such tie to the exact lexicographic hkl tie-break instead. The quanta are
    # far below anything physical (1 pm on the detector, 1e-12 of full scale in
    # intensity) and far above the ~1e-14 noise they suppress.
    order = np.lexsort(
        (
            hkl[:, 2],
            hkl[:, 1],
            hkl[:, 0],
            np.round(radius, _RADIUS_SORT_DECIMALS),
            -np.round(intensity, _INTENSITY_SORT_DECIMALS),
        )
    )

    return SpotTable(
        phase=phase,
        zone_axis=zone_axis,
        basis=zone_basis,
        config=simulation_config,
        hkl=hkl[order],
        g_crystal=g_crystal[order],
        g_detector_inv_angstrom=g_detector[order],
        detector_mm=detector_mm[order],
        d_spacing_angstrom=d_spacing[order],
        structure_factor_amplitude=amplitude[order],
        intensity=intensity[order],
        excitation_error_inv_angstrom=excitation_error[order],
        provenance=provenance,
        is_double_diffraction=is_double_diffraction[order],
        double_diffraction_parents=parents[order],
    )


__all__ = [
    "IntensityModelName",
    "KinematicSimulationConfig",
    "SpotTable",
    "centering_allowed_mask",
    "double_diffraction_sums",
    "electron_structure_factors",
    "electron_wavelength_angstrom",
    "simulate_zone_axis_spots",
    "zone_basis_from_axis",
]
