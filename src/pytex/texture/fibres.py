"""Fibre textures: a crystal direction parallel to a specimen direction.

A `Fibre` is the one-parameter family of orientations that map one crystal
direction onto one specimen direction. Fibre distances are symmetry-aware and
antipodal-aware: the reported angle is the minimum over the symmetry family
of the crystal direction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase
from pytex.core.orientation import (
    OrientationSet,
    quaternion_from_axis_angle,
    quaternions_multiply,
    specimen_direction_vector,
)
from pytex.core.symmetry import SymmetrySpec


def _alignment_quaternion(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    cosine = float(np.clip(source @ target, -1.0, 1.0))
    if np.isclose(cosine, 1.0):
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    if np.isclose(cosine, -1.0):
        helper = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if np.isclose(abs(float(source @ helper)), 1.0):
            helper = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        axis = normalize_vector(np.cross(source, helper))
        return quaternion_from_axis_angle(axis, np.pi)
    axis = normalize_vector(np.cross(source, target))
    return quaternion_from_axis_angle(axis, float(np.arccos(cosine)))


@dataclass(frozen=True, slots=True)
class Fibre:
    name: str
    crystal_direction: tuple[float, float, float]
    specimen_direction: str | tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Fibre.name must be a non-empty string.")
        direction = tuple(float(value) for value in self.crystal_direction)
        if len(direction) != 3 or np.isclose(float(np.linalg.norm(direction)), 0.0):
            raise ValueError("Fibre.crystal_direction must be a non-zero 3-vector.")
        object.__setattr__(self, "crystal_direction", direction)
        # Validate the specimen direction eagerly so bad aliases fail here.
        specimen_direction_vector(self.specimen_direction)

    @classmethod
    def alpha_bcc(cls) -> Fibre:
        """BCC alpha fibre: <110> parallel to the rolling direction."""

        return cls("alpha_bcc", (1.0, 1.0, 0.0), "RD")

    @classmethod
    def gamma_bcc(cls) -> Fibre:
        """BCC gamma fibre: <111> parallel to the sheet normal."""

        return cls("gamma_bcc", (1.0, 1.0, 1.0), "ND")

    @classmethod
    def eta(cls) -> Fibre:
        """Eta fibre: <100> parallel to the rolling direction."""

        return cls("eta", (1.0, 0.0, 0.0), "RD")

    @classmethod
    def theta(cls) -> Fibre:
        """Theta fibre: <100> parallel to the sheet normal."""

        return cls("theta", (1.0, 0.0, 0.0), "ND")

    def crystal_unit_vector(self) -> np.ndarray:
        return normalize_vector(np.asarray(self.crystal_direction, dtype=np.float64))

    def specimen_unit_vector(self) -> np.ndarray:
        return specimen_direction_vector(self.specimen_direction)

    def orientations(
        self,
        count: int,
        *,
        specimen_frame: ReferenceFrame,
        crystal_frame: ReferenceFrame | None = None,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
    ) -> OrientationSet:
        """Sample the fibre with `count` evenly spaced rotations about the axis."""

        if count < 1:
            raise ValueError("Fibre sampling requires at least one orientation.")
        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        crystal_axis = self.crystal_unit_vector()
        specimen_axis = self.specimen_unit_vector()
        base = _alignment_quaternion(crystal_axis, specimen_axis)
        angles = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
        spins = np.stack(
            [quaternion_from_axis_angle(specimen_axis, float(angle)) for angle in angles],
            axis=0,
        )
        quaternions = quaternions_multiply(spins, base[None, :])
        return OrientationSet(
            quaternions=quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
        )

    def angles_to_deg(self, orientations: OrientationSet) -> np.ndarray:
        """Per-orientation fibre distance in degrees (symmetry- and antipodal-aware)."""

        specimen_axis = self.specimen_unit_vector()
        crystal_axis = self.crystal_unit_vector()
        if orientations.symmetry is not None:
            family = orientations.symmetry.equivalent_vectors(crystal_axis, antipodal=True)
        else:
            family = np.stack([crystal_axis, -crystal_axis], axis=0)
        best = np.full(len(orientations), np.inf, dtype=np.float64)
        for member in family:
            mapped = np.asarray(orientations.map_crystal_directions(member))
            cosines = np.clip(np.abs(mapped @ specimen_axis), 0.0, 1.0)
            best = np.minimum(best, np.arccos(cosines))
        angles = np.rad2deg(best)
        angles = np.ascontiguousarray(angles)
        angles.setflags(write=False)
        return angles

    def volume_fraction(
        self,
        orientations: OrientationSet,
        *,
        tolerance_deg: float = 10.0,
        weights: ArrayLike | None = None,
    ) -> float:
        if not 0.0 < float(tolerance_deg) <= 90.0:
            raise ValueError("tolerance_deg must lie in (0, 90] degrees.")
        count = len(orientations)
        if count == 0:
            raise ValueError("volume_fraction requires at least one orientation.")
        if weights is None:
            weight_values = np.full(count, 1.0 / count, dtype=np.float64)
        else:
            weight_values = np.asarray(weights, dtype=np.float64)
            if weight_values.shape != (count,):
                raise ValueError("weights must provide one value per orientation.")
            if np.any(weight_values < 0.0):
                raise ValueError("weights must be non-negative.")
            total = float(weight_values.sum())
            if np.isclose(total, 0.0):
                raise ValueError("weights must not sum to zero.")
            weight_values = weight_values / total
        angles = self.angles_to_deg(orientations)
        return float(weight_values[angles <= float(tolerance_deg)].sum())


NAMED_BCC_FIBRES: tuple[Fibre, ...] = (
    Fibre.alpha_bcc(),
    Fibre.gamma_bcc(),
    Fibre.eta(),
    Fibre.theta(),
)


def fibre_axis_alignment_quaternion(source: ArrayLike, target: ArrayLike) -> np.ndarray:
    """Quaternion mapping one unit direction onto another (shared helper)."""

    return as_float_array(
        _alignment_quaternion(
            normalize_vector(source),
            normalize_vector(target),
        ),
        shape=(4,),
    )


__all__ = [
    "NAMED_BCC_FIBRES",
    "Fibre",
    "fibre_axis_alignment_quaternion",
]
