"""Slip systems and Schmid-factor analysis.

A `SlipSystem` pairs a crystallographic slip plane ``{hkl}`` with a slip
direction ``<uvw>`` lying in that plane. `SlipSystemFamily` expands a
representative slip system into its full symmetry-distinct set using the crystal
point-group operators, and computes vectorized Schmid factors for orientation
populations and crystal maps.

The Schmid factor for a slip system under a uniaxial stress along the unit
specimen direction ``t`` is ``m = |(t . n)(t . b)|`` where ``n`` is the unit
slip-plane normal and ``b`` the unit slip direction, both expressed in the
specimen frame via the crystal orientation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.lattice import Phase
from pytex.core.miller import MillerDirection, MillerPlane
from pytex.core.orientation import Orientation, OrientationSet


@dataclass(frozen=True, slots=True)
class SlipSystem:
    """A slip plane normal paired with an in-plane slip direction."""

    plane: MillerPlane
    direction: MillerDirection
    name: str = ""

    def __post_init__(self) -> None:
        if self.plane.phase != self.direction.phase:
            raise ValueError("SlipSystem plane and direction must share a phase.")
        normal = self.plane_normal
        slip = self.slip_direction
        if abs(float(np.dot(normal, slip))) > 1e-6:
            raise ValueError(
                "SlipSystem slip direction must lie in the slip plane "
                "(plane normal and slip direction must be orthogonal)."
            )

    @property
    def phase(self) -> Phase:
        return self.plane.phase

    @property
    def plane_normal(self) -> np.ndarray:
        return np.asarray(self.plane.normal_cartesian, dtype=np.float64)

    @property
    def slip_direction(self) -> np.ndarray:
        return np.asarray(self.direction.unit_vector_cartesian, dtype=np.float64)


def _dedupe_axis_pairs(
    normals: np.ndarray,
    directions: np.ndarray,
    *,
    tol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Keep symmetry-distinct (normal, direction) pairs, treating +/- as equal."""

    kept_normals: list[np.ndarray] = []
    kept_directions: list[np.ndarray] = []
    for normal, direction in zip(normals, directions, strict=True):
        duplicate = False
        for existing_normal, existing_direction in zip(
            kept_normals, kept_directions, strict=True
        ):
            same_plane = np.allclose(np.abs(existing_normal @ normal), 1.0, atol=tol)
            same_dir = np.allclose(np.abs(existing_direction @ direction), 1.0, atol=tol)
            if same_plane and same_dir:
                duplicate = True
                break
        if not duplicate:
            kept_normals.append(normal)
            kept_directions.append(direction)
    return np.asarray(kept_normals, dtype=np.float64), np.asarray(
        kept_directions, dtype=np.float64
    )


@dataclass(frozen=True, slots=True)
class SlipSystemFamily:
    """The full symmetry-distinct set of slip systems for a representative."""

    name: str
    representative: SlipSystem
    plane_normals: np.ndarray
    slip_directions: np.ndarray

    def __post_init__(self) -> None:
        normals = np.ascontiguousarray(np.asarray(self.plane_normals, dtype=np.float64))
        directions = np.ascontiguousarray(np.asarray(self.slip_directions, dtype=np.float64))
        if normals.ndim != 2 or normals.shape[1] != 3:
            raise ValueError("SlipSystemFamily.plane_normals must have shape (m, 3).")
        if directions.shape != normals.shape:
            raise ValueError("SlipSystemFamily plane_normals and slip_directions must match.")
        normals.setflags(write=False)
        directions.setflags(write=False)
        object.__setattr__(self, "plane_normals", normals)
        object.__setattr__(self, "slip_directions", directions)

    @property
    def count(self) -> int:
        return int(self.plane_normals.shape[0])

    @classmethod
    def from_representative(cls, representative: SlipSystem, *, name: str = "") -> SlipSystemFamily:
        phase = representative.phase
        if phase.symmetry is None:
            raise ValueError("SlipSystemFamily expansion requires a phase with crystal symmetry.")
        operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
        base_normal = representative.plane_normal
        base_direction = representative.slip_direction
        normals = np.einsum("oij,j->oi", operators, base_normal, optimize=True)
        directions = np.einsum("oij,j->oi", operators, base_direction, optimize=True)
        unique_normals, unique_directions = _dedupe_axis_pairs(normals, directions)
        return cls(
            name=name or representative.name,
            representative=representative,
            plane_normals=unique_normals,
            slip_directions=unique_directions,
        )

    def schmid_factors(
        self,
        orientations: Orientation | OrientationSet,
        stress_direction: ArrayLike,
    ) -> np.ndarray:
        """Per-orientation, per-system absolute Schmid factors.

        Returns an array of shape ``(n_orientations, n_systems)`` (or
        ``(n_systems,)`` for a single `Orientation`).
        """

        stress = np.asarray(stress_direction, dtype=np.float64)
        norm = float(np.linalg.norm(stress))
        if norm == 0.0:
            raise ValueError("stress_direction must be a non-zero vector.")
        stress = stress / norm
        if isinstance(orientations, Orientation):
            scalar = True
            orientation_set = OrientationSet.from_orientations([orientations])
        else:
            scalar = False
            orientation_set = orientations
        matrices = orientation_set.as_matrices()
        # Map crystal-frame family vectors into the specimen frame: v_spec = R @ v_crystal.
        normals_spec = np.einsum("nij,mj->nmi", matrices, self.plane_normals, optimize=True)
        directions_spec = np.einsum("nij,mj->nmi", matrices, self.slip_directions, optimize=True)
        cos_phi = normals_spec @ stress
        cos_lambda = directions_spec @ stress
        factors = np.abs(cos_phi * cos_lambda)
        factors = np.ascontiguousarray(factors)
        factors.setflags(write=False)
        return factors[0] if scalar else factors

    def max_schmid_factor(
        self,
        orientations: Orientation | OrientationSet,
        stress_direction: ArrayLike,
    ) -> np.ndarray | float:
        """Maximum absolute Schmid factor over the family per orientation."""

        factors = self.schmid_factors(orientations, stress_direction)
        if factors.ndim == 1:
            return float(np.max(factors))
        maxima = np.max(factors, axis=1)
        maxima = np.ascontiguousarray(maxima)
        maxima.setflags(write=False)
        return maxima


def fcc_octahedral_slip(phase: Phase) -> SlipSystemFamily:
    """The 12 fcc octahedral slip systems ``{111}<1-10>``."""

    representative = SlipSystem(
        plane=MillerPlane.from_hkl((1, 1, 1), phase=phase),
        direction=MillerDirection.from_uvw((1, -1, 0), phase=phase),
        name="fcc_{111}<110>",
    )
    return SlipSystemFamily.from_representative(representative, name="fcc_octahedral")


def bcc_110_slip(phase: Phase) -> SlipSystemFamily:
    """The 12 bcc slip systems ``{110}<-111>``."""

    representative = SlipSystem(
        plane=MillerPlane.from_hkl((1, 1, 0), phase=phase),
        direction=MillerDirection.from_uvw((-1, 1, 1), phase=phase),
        name="bcc_{110}<111>",
    )
    return SlipSystemFamily.from_representative(representative, name="bcc_110")


__all__ = [
    "SlipSystem",
    "SlipSystemFamily",
    "bcc_110_slip",
    "fcc_octahedral_slip",
]
