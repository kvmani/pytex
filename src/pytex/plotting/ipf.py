from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector, normalize_vectors
from pytex.core.conventions import FrameDomain
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.symmetry import SymmetrySpec


def _normalize_rgb_vertices(vertices: ArrayLike) -> np.ndarray:
    array = as_float_array(vertices, shape=(3, 3))
    if np.any(~np.isfinite(array)) or np.any((array < 0.0) | (array > 1.0)):
        raise ValueError("IPF color vertices must be finite RGB triplets in [0, 1].")
    return array


@dataclass(frozen=True, slots=True)
class IPFColorKey:
    crystal_symmetry: SymmetrySpec
    specimen_direction: np.ndarray
    antipodal: bool = True
    saturation_gamma: float = 0.5
    color_vertices: np.ndarray = field(
        default_factory=lambda: np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
    )

    def __post_init__(self) -> None:
        if self.crystal_symmetry.reference_frame is None:
            raise ValueError("IPFColorKey.crystal_symmetry must declare a crystal reference frame.")
        if self.crystal_symmetry.reference_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError(
                "IPFColorKey.crystal_symmetry must use a crystal-domain reference frame."
            )
        object.__setattr__(self, "specimen_direction", normalize_vector(self.specimen_direction))
        object.__setattr__(self, "color_vertices", _normalize_rgb_vertices(self.color_vertices))
        if not np.isfinite(self.saturation_gamma) or self.saturation_gamma <= 0.0:
            raise ValueError("IPFColorKey.saturation_gamma must be finite and strictly positive.")

    @property
    def sector_vertices(self) -> np.ndarray:
        return self.crystal_symmetry.fundamental_sector(antipodal=self.antipodal).vertices

    def colors_from_crystal_directions(self, crystal_directions: ArrayLike) -> np.ndarray:
        directions = normalize_vectors(crystal_directions)
        reduced = self.crystal_symmetry.reduce_vectors_to_fundamental_sector(
            directions,
            antipodal=self.antipodal,
        )
        sector_basis = np.asarray(self.sector_vertices, dtype=np.float64).T
        barycentric = np.linalg.solve(sector_basis, reduced.T).T
        barycentric = np.clip(barycentric, 0.0, None)
        sums = barycentric.sum(axis=1, keepdims=True)
        if np.any(np.isclose(sums, 0.0)):
            raise ValueError("Reduced crystal directions must lie inside the IPF sector cone.")
        normalized = barycentric / sums
        colors = normalized @ self.color_vertices
        colors = np.clip(colors, 0.0, 1.0) ** (1.0 / self.saturation_gamma)
        max_values = np.max(colors, axis=1, keepdims=True)
        colors = colors / np.where(max_values > 0.0, max_values, 1.0)
        colors = np.ascontiguousarray(colors)
        colors.setflags(write=False)
        return colors

    def colors_from_orientations(self, orientations: OrientationSet) -> np.ndarray:
        if orientations.symmetry is None or orientations.symmetry != self.crystal_symmetry:
            raise ValueError("OrientationSet symmetry must match the IPFColorKey crystal symmetry.")
        crystal_directions = orientations.map_sample_directions_to_crystal(self.specimen_direction)
        return self.colors_from_crystal_directions(crystal_directions)

    def color_from_orientation(self, orientation: Orientation) -> np.ndarray:
        if orientation.symmetry is None or orientation.symmetry != self.crystal_symmetry:
            raise ValueError("Orientation symmetry must match the IPFColorKey crystal symmetry.")
        crystal_direction = orientation.map_sample_vector_to_crystal(self.specimen_direction)
        return as_float_array(
            self.colors_from_crystal_directions(crystal_direction[None, :])[0],
            shape=(3,),
        )
