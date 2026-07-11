"""Canonical S2 (unit-sphere) direction semantics.

Spherical-angle convention used across PyTex: the polar angle is measured from the
+Z axis of the owning reference frame, and the azimuth is measured from +X toward
+Y. Public angle arguments and returns are degrees unless a ``_rad`` suffix says
otherwise, matching the stereonet surface this module canonicalizes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, freeze_array, normalize_vectors
from pytex.core.batches import VectorSet
from pytex.core.provenance import ProvenanceRecord

if TYPE_CHECKING:
    from pytex.core.frames import ReferenceFrame
    from pytex.core.orientation import Rotation

_HEMISPHERES = ("upper", "sphere")


def spherical_angles_to_directions(
    polar_deg: ArrayLike,
    azimuth_deg: ArrayLike,
) -> np.ndarray:
    polar, azimuth = np.broadcast_arrays(
        np.asarray(polar_deg, dtype=np.float64),
        np.asarray(azimuth_deg, dtype=np.float64),
    )
    polar_rad = np.deg2rad(polar)
    azimuth_rad = np.deg2rad(azimuth)
    directions = np.stack(
        [
            np.sin(polar_rad) * np.cos(azimuth_rad),
            np.sin(polar_rad) * np.sin(azimuth_rad),
            np.cos(polar_rad),
        ],
        axis=-1,
    )
    return freeze_array(np.ascontiguousarray(directions, dtype=np.float64))


def directions_to_spherical_angles(
    directions: ArrayLike,
    *,
    antipodal: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    vectors = np.array(normalize_vectors(directions), copy=True)
    if antipodal:
        mask = vectors[..., 2] < 0.0
        vectors[mask] *= -1.0
    polar = np.rad2deg(np.arccos(np.clip(vectors[..., 2], -1.0, 1.0)))
    azimuth = np.mod(np.rad2deg(np.arctan2(vectors[..., 1], vectors[..., 0])), 360.0)
    polar = freeze_array(np.ascontiguousarray(polar, dtype=np.float64))
    azimuth = freeze_array(np.ascontiguousarray(azimuth, dtype=np.float64))
    return polar, azimuth


def _broadcast_unit_rows(
    left: np.ndarray,
    right: np.ndarray,
    *,
    operation: str,
) -> tuple[np.ndarray, np.ndarray]:
    if left.shape[0] == right.shape[0]:
        return left, right
    if left.shape[0] == 1:
        return np.broadcast_to(left, right.shape), right
    if right.shape[0] == 1:
        return left, np.broadcast_to(right, left.shape)
    raise ValueError(
        f"Cannot broadcast spherical vector sets of lengths {left.shape[0]} and "
        f"{right.shape[0]} for {operation}: lengths must match or one must be 1."
    )


@dataclass(frozen=True, slots=True)
class SphericalVectorSet:
    values: np.ndarray
    reference_frame: ReferenceFrame
    antipodal: bool = False
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", normalize_vectors(self.values))

    @classmethod
    def from_vectors(
        cls,
        vectors: ArrayLike,
        *,
        reference_frame: ReferenceFrame,
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> SphericalVectorSet:
        return cls(
            values=np.asarray(vectors, dtype=np.float64),
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )

    @classmethod
    def from_vector_set(
        cls,
        vector_set: VectorSet,
        *,
        antipodal: bool = False,
    ) -> SphericalVectorSet:
        return cls(
            values=vector_set.values,
            reference_frame=vector_set.reference_frame,
            antipodal=antipodal,
            provenance=vector_set.provenance,
        )

    @classmethod
    def from_polar(
        cls,
        polar: ArrayLike,
        azimuth: ArrayLike,
        *,
        reference_frame: ReferenceFrame,
        degrees: bool = True,
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> SphericalVectorSet:
        polar_values = np.atleast_1d(np.asarray(polar, dtype=np.float64))
        azimuth_values = np.atleast_1d(np.asarray(azimuth, dtype=np.float64))
        if not degrees:
            polar_values = np.rad2deg(polar_values)
            azimuth_values = np.rad2deg(azimuth_values)
        directions = spherical_angles_to_directions(polar_values, azimuth_values)
        return cls(
            values=directions.reshape(-1, 3),
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | SphericalVectorSet:
        selected = self.values[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return SphericalVectorSet(
            values=selected,
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        return self.values

    def to_vector_set(self) -> VectorSet:
        return VectorSet(
            values=self.values,
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def to_polar(self, *, degrees: bool = True) -> tuple[np.ndarray, np.ndarray]:
        polar_deg, azimuth_deg = directions_to_spherical_angles(
            self.values,
            antipodal=self.antipodal,
        )
        if degrees:
            return polar_deg, azimuth_deg
        polar_rad = freeze_array(np.ascontiguousarray(np.deg2rad(polar_deg)))
        azimuth_rad = freeze_array(np.ascontiguousarray(np.deg2rad(azimuth_deg)))
        return polar_rad, azimuth_rad

    def subset(self, indices: ArrayLike) -> SphericalVectorSet:
        return SphericalVectorSet(
            values=self.values[np.asarray(indices)],
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )

    def fold_upper_hemisphere(self) -> SphericalVectorSet:
        if not self.antipodal:
            raise ValueError(
                "fold_upper_hemisphere is only meaningful for antipodal direction sets; "
                "construct the set with antipodal=True if +v and -v are equivalent."
            )
        folded = np.array(self.values, copy=True)
        lower = folded[:, 2] < 0.0
        equator_negative = (np.isclose(folded[:, 2], 0.0)) & (
            (folded[:, 0] < 0.0) | (np.isclose(folded[:, 0], 0.0) & (folded[:, 1] < 0.0))
        )
        folded[lower | equator_negative] *= -1.0
        return SphericalVectorSet(
            values=folded,
            reference_frame=self.reference_frame,
            antipodal=True,
            provenance=self.provenance,
        )

    def _require_matching_frame(self, other: SphericalVectorSet, *, operation: str) -> None:
        if self.reference_frame != other.reference_frame:
            raise ValueError(
                f"SphericalVectorSet {operation} requires both operands to share one "
                "reference frame."
            )

    def dot(self, other: SphericalVectorSet) -> np.ndarray:
        self._require_matching_frame(other, operation="dot")
        left, right = _broadcast_unit_rows(self.values, other.values, operation="dot")
        return freeze_array(np.ascontiguousarray(np.einsum("ni,ni->n", left, right)))

    def cross(self, other: SphericalVectorSet) -> SphericalVectorSet:
        self._require_matching_frame(other, operation="cross")
        left, right = _broadcast_unit_rows(self.values, other.values, operation="cross")
        products = np.cross(left, right)
        norms = np.linalg.norm(products, axis=1)
        if np.any(np.isclose(norms, 0.0)):
            raise ValueError(
                "cross is undefined for parallel or antiparallel direction pairs."
            )
        return SphericalVectorSet(
            values=products,
            reference_frame=self.reference_frame,
            antipodal=self.antipodal or other.antipodal,
            provenance=self.provenance,
        )

    def angles_to_rad(self, other: SphericalVectorSet) -> np.ndarray:
        self._require_matching_frame(other, operation="angles_to_rad")
        left, right = _broadcast_unit_rows(
            self.values,
            other.values,
            operation="angles_to_rad",
        )
        cosines = np.clip(np.einsum("ni,ni->n", left, right), -1.0, 1.0)
        if self.antipodal or other.antipodal:
            cosines = np.abs(cosines)
        return freeze_array(np.ascontiguousarray(np.arccos(cosines)))

    def angles_to_deg(self, other: SphericalVectorSet) -> np.ndarray:
        return freeze_array(np.ascontiguousarray(np.rad2deg(self.angles_to_rad(other))))

    def orientation_tensor(self) -> np.ndarray:
        tensor = np.einsum("ni,nj->ij", self.values, self.values) / float(len(self))
        return freeze_array(np.ascontiguousarray(tensor))

    def mean_direction(self) -> np.ndarray:
        if self.antipodal:
            eigenvalues, eigenvectors = np.linalg.eigh(self.orientation_tensor())
            principal = eigenvectors[:, int(np.argmax(eigenvalues))]
            if principal[2] < 0.0 or (
                np.isclose(principal[2], 0.0)
                and (
                    principal[0] < 0.0
                    or (np.isclose(principal[0], 0.0) and principal[1] < 0.0)
                )
            ):
                principal = -principal
            return as_float_array(principal, shape=(3,))
        resultant = self.values.sum(axis=0)
        norm = float(np.linalg.norm(resultant))
        if np.isclose(norm, 0.0):
            raise ValueError(
                "Mean direction is undefined: the resultant vector is numerically zero."
            )
        return as_float_array(resultant / norm, shape=(3,))

    def rotated_by(self, rotation: Rotation) -> SphericalVectorSet:
        mapped = rotation.apply(self.values)
        return SphericalVectorSet(
            values=np.asarray(mapped, dtype=np.float64),
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )


def _require_hemisphere(hemisphere: str) -> str:
    if hemisphere not in _HEMISPHERES:
        supported = ", ".join(_HEMISPHERES)
        raise ValueError(
            f"Unsupported hemisphere '{hemisphere}'. Supported values: {supported}."
        )
    return hemisphere


def _require_resolution(resolution_deg: float) -> float:
    resolution = float(resolution_deg)
    if not 0.0 < resolution <= 90.0:
        raise ValueError("Grid resolution must lie in the interval (0, 90] degrees.")
    return resolution


def _ring_band_weights(
    polar_ring_deg: np.ndarray,
    counts: np.ndarray,
    *,
    polar_max_deg: float,
    half_band_deg: float,
) -> np.ndarray:
    lower = np.clip(polar_ring_deg - half_band_deg, 0.0, polar_max_deg)
    upper = np.clip(polar_ring_deg + half_band_deg, 0.0, polar_max_deg)
    band_measure = np.cos(np.deg2rad(lower)) - np.cos(np.deg2rad(upper))
    per_point = np.repeat(band_measure / counts, counts)
    return np.asarray(per_point / per_point.sum(), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class S2Grid:
    vectors: SphericalVectorSet
    weights: np.ndarray
    resolution_deg: float
    hemisphere: str
    method: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "weights", as_float_array(self.weights, shape=(None,)))
        object.__setattr__(self, "hemisphere", _require_hemisphere(self.hemisphere))
        object.__setattr__(self, "resolution_deg", _require_resolution(self.resolution_deg))
        if self.weights.shape[0] != len(self.vectors):
            raise ValueError("S2Grid.weights must have one weight per grid direction.")
        if np.any(self.weights <= 0.0):
            raise ValueError("S2Grid.weights must be strictly positive.")
        if not np.isclose(float(self.weights.sum()), 1.0, atol=1e-10):
            raise ValueError("S2Grid.weights must sum to 1 (normalized surface measure).")

    def __len__(self) -> int:
        return len(self.vectors)

    @classmethod
    def equispaced(
        cls,
        resolution_deg: float,
        *,
        reference_frame: ReferenceFrame,
        hemisphere: str = "upper",
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> S2Grid:
        resolution = _require_resolution(resolution_deg)
        hemisphere = _require_hemisphere(hemisphere)
        polar_max = 90.0 if hemisphere == "upper" else 180.0
        ring_count = int(round(polar_max / resolution))
        polar_rings = np.linspace(0.0, polar_max, ring_count + 1)

        ring_polar: list[np.ndarray] = []
        ring_azimuth: list[np.ndarray] = []
        counts = np.empty(polar_rings.shape[0], dtype=np.int64)
        for ring_index, polar_deg in enumerate(polar_rings):
            circumference_deg = 360.0 * float(np.sin(np.deg2rad(polar_deg)))
            count = max(1, int(round(circumference_deg / resolution)))
            counts[ring_index] = count
            azimuths = np.linspace(0.0, 360.0, count, endpoint=False)
            ring_polar.append(np.full(count, float(polar_deg)))
            ring_azimuth.append(azimuths)

        polar_all = np.concatenate(ring_polar)
        azimuth_all = np.concatenate(ring_azimuth)
        weights = _ring_band_weights(
            polar_rings,
            counts,
            polar_max_deg=polar_max,
            half_band_deg=resolution / 2.0,
        )
        vectors = SphericalVectorSet.from_polar(
            polar_all,
            azimuth_all,
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )
        return cls(
            vectors=vectors,
            weights=weights,
            resolution_deg=resolution,
            hemisphere=hemisphere,
            method="equispaced",
        )

    @classmethod
    def regular(
        cls,
        polar_step_deg: float,
        azimuth_step_deg: float,
        *,
        reference_frame: ReferenceFrame,
        hemisphere: str = "upper",
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> S2Grid:
        polar_step = _require_resolution(polar_step_deg)
        azimuth_step = float(azimuth_step_deg)
        if not 0.0 < azimuth_step <= 120.0 or not np.isclose(
            np.mod(360.0, azimuth_step), 0.0, atol=1e-10
        ):
            raise ValueError(
                "azimuth_step_deg must lie in (0, 120] and divide 360 degrees evenly."
            )
        hemisphere = _require_hemisphere(hemisphere)
        polar_max = 90.0 if hemisphere == "upper" else 180.0
        ring_count = int(round(polar_max / polar_step))
        polar_rings = np.linspace(0.0, polar_max, ring_count + 1)
        azimuth_count = int(round(360.0 / azimuth_step))

        ring_polar: list[np.ndarray] = []
        ring_azimuth: list[np.ndarray] = []
        counts = np.empty(polar_rings.shape[0], dtype=np.int64)
        for ring_index, polar_deg in enumerate(polar_rings):
            at_pole = np.isclose(float(polar_deg), 0.0) or np.isclose(float(polar_deg), 180.0)
            count = 1 if at_pole else azimuth_count
            counts[ring_index] = count
            azimuths = np.linspace(0.0, 360.0, count, endpoint=False)
            ring_polar.append(np.full(count, float(polar_deg)))
            ring_azimuth.append(azimuths)

        polar_all = np.concatenate(ring_polar)
        azimuth_all = np.concatenate(ring_azimuth)
        weights = _ring_band_weights(
            polar_rings,
            counts,
            polar_max_deg=polar_max,
            half_band_deg=polar_step / 2.0,
        )
        vectors = SphericalVectorSet.from_polar(
            polar_all,
            azimuth_all,
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )
        return cls(
            vectors=vectors,
            weights=weights,
            resolution_deg=polar_step,
            hemisphere=hemisphere,
            method="regular",
        )
