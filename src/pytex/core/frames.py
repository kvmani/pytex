from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.conventions import (
    PYTEX_CANONICAL_CONVENTIONS,
    ConventionSet,
    FrameDomain,
    Handedness,
)
from pytex.core.provenance import ProvenanceRecord


@dataclass(frozen=True, slots=True)
class ReferenceFrame:
    name: str
    domain: FrameDomain
    axes: tuple[str, str, str]
    handedness: Handedness = Handedness.RIGHT
    convention: ConventionSet = PYTEX_CANONICAL_CONVENTIONS
    description: str = ""
    provenance: ProvenanceRecord | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.axes) != 3:
            raise ValueError("ReferenceFrame.axes must contain exactly three axis labels.")
        object.__setattr__(self, "axes", tuple(self.axes))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, slots=True)
class FrameTransform:
    source: ReferenceFrame
    target: ReferenceFrame
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray = field(default_factory=lambda: np.zeros(3))
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        rotation = as_float_array(self.rotation_matrix, shape=(3, 3))
        translation = as_float_array(self.translation_vector, shape=(3,))
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8):
            raise ValueError("rotation_matrix must be orthonormal.")
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8):
            raise ValueError("rotation_matrix must have determinant +1.")
        object.__setattr__(self, "rotation_matrix", rotation)
        object.__setattr__(self, "translation_vector", translation)

    @classmethod
    def identity(cls, frame: ReferenceFrame) -> FrameTransform:
        return cls(source=frame, target=frame, rotation_matrix=np.eye(3))

    def apply_to_vectors(self, vectors: ArrayLike) -> np.ndarray:
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = array @ self.rotation_matrix.T + self.translation_vector
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def inverse(self) -> FrameTransform:
        inverse_rotation = self.rotation_matrix.T
        inverse_translation = -(self.translation_vector @ self.rotation_matrix)
        return FrameTransform(
            source=self.target,
            target=self.source,
            rotation_matrix=inverse_rotation,
            translation_vector=inverse_translation,
            provenance=self.provenance,
        )

    def compose(self, previous: FrameTransform) -> FrameTransform:
        if previous.target != self.source:
            raise ValueError(
                "Transform domains do not chain: previous.target must equal self.source."
            )
        rotation = self.rotation_matrix @ previous.rotation_matrix
        translation = previous.translation_vector @ self.rotation_matrix.T + self.translation_vector
        return FrameTransform(
            source=previous.source,
            target=self.target,
            rotation_matrix=rotation,
            translation_vector=translation,
            provenance=self.provenance,
        )
