from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_quaternions, normalize_vectors
from pytex.core.provenance import ProvenanceRecord

if TYPE_CHECKING:
    from pytex.core.frames import ReferenceFrame
    from pytex.core.orientation import Rotation


def normalize_euler_convention_name(convention: str) -> str:
    normalized = convention.strip().lower()
    aliases = {
        "bunge": "bunge",
        "bunge_zxz": "bunge",
        "zxz": "bunge",
        "matthies": "matthies",
        "matthies_zyz": "matthies",
        "abg": "abg",
        "abg_zyz": "abg",
        "zyz": "abg",
    }
    resolved = aliases.get(normalized)
    if resolved is None:
        supported = ", ".join(sorted(aliases))
        raise ValueError(
            f"Unsupported Euler convention '{convention}'. Supported conventions: {supported}"
        )
    return resolved


@dataclass(frozen=True, slots=True)
class VectorSet:
    values: np.ndarray
    reference_frame: ReferenceFrame
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", as_float_array(self.values, shape=(None, 3)))

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | VectorSet:
        selected = self.values[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return VectorSet(
            values=selected,
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        return self.values

    def normalized(self) -> VectorSet:
        return VectorSet(
            values=normalize_vectors(self.values),
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def subset(self, indices: ArrayLike) -> VectorSet:
        return VectorSet(
            values=self.values[np.asarray(indices)],
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )


@dataclass(frozen=True, slots=True)
class EulerSet:
    angles: np.ndarray
    convention: str = "bunge"
    degrees: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "angles", as_float_array(self.angles, shape=(None, 3)))
        object.__setattr__(self, "convention", normalize_euler_convention_name(self.convention))

    def __len__(self) -> int:
        return int(self.angles.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | EulerSet:
        selected = self.angles[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return EulerSet(
            angles=selected,
            convention=self.convention,
            degrees=self.degrees,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        return self.angles

    def subset(self, indices: ArrayLike) -> EulerSet:
        return EulerSet(
            angles=self.angles[np.asarray(indices)],
            convention=self.convention,
            degrees=self.degrees,
            provenance=self.provenance,
        )

    def to_rotation_set(self) -> RotationSet:
        return RotationSet.from_euler_set(self)


@dataclass(frozen=True, slots=True)
class QuaternionSet:
    quaternions: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternions", normalize_quaternions(self.quaternions))

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | QuaternionSet:
        selected = self.quaternions[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(4,))
        return QuaternionSet(quaternions=selected, provenance=self.provenance)

    def as_array(self) -> np.ndarray:
        return self.quaternions

    def subset(self, indices: ArrayLike) -> QuaternionSet:
        return QuaternionSet(
            quaternions=self.quaternions[np.asarray(indices)],
            provenance=self.provenance,
        )

    def to_rotation_set(self) -> RotationSet:
        return RotationSet(quaternions=self.quaternions, provenance=self.provenance)


@dataclass(frozen=True, slots=True)
class RotationSet:
    quaternions: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternions", normalize_quaternions(self.quaternions))

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    def __getitem__(self, index: Any) -> Rotation | RotationSet:
        selected = self.quaternions[index]
        if np.asarray(selected).ndim == 1:
            from pytex.core.orientation import Rotation

            return Rotation(quaternion=selected, provenance=self.provenance)
        return RotationSet(quaternions=selected, provenance=self.provenance)

    @classmethod
    def from_rotations(cls, rotations: list[Rotation]) -> RotationSet:
        quaternions = np.stack([rotation.quaternion for rotation in rotations], axis=0)
        provenance = rotations[0].provenance if rotations else None
        return cls(quaternions=quaternions, provenance=provenance)

    @classmethod
    def from_euler_set(cls, euler_set: EulerSet) -> RotationSet:
        from pytex.core.orientation import Rotation

        quaternions = np.stack(
            [
                Rotation.from_euler(
                    angle1,
                    angle2,
                    angle3,
                    convention=euler_set.convention,
                    degrees=euler_set.degrees,
                ).quaternion
                for angle1, angle2, angle3 in euler_set.angles
            ],
            axis=0,
        )
        return cls(quaternions=quaternions, provenance=euler_set.provenance)

    def as_quaternion_set(self) -> QuaternionSet:
        return QuaternionSet(quaternions=self.quaternions, provenance=self.provenance)

    def as_matrices(self) -> np.ndarray:
        from pytex.core.orientation import quaternion_to_matrix

        matrices = np.stack(
            [quaternion_to_matrix(quaternion) for quaternion in self.quaternions],
            axis=0,
        )
        matrices = np.ascontiguousarray(matrices)
        matrices.setflags(write=False)
        return matrices

    def as_euler_set(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> EulerSet:
        from pytex.core.orientation import Rotation

        angles = np.stack(
            [
                Rotation(quaternion=quaternion).to_euler(convention=convention, degrees=degrees)
                for quaternion in self.quaternions
            ],
            axis=0,
        )
        return EulerSet(
            angles=angles,
            convention=convention,
            degrees=degrees,
            provenance=self.provenance,
        )

    def apply(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        matrices = self.as_matrices()
        if isinstance(vectors, VectorSet):
            values = vectors.values
            if values.shape[0] != len(self):
                raise ValueError("VectorSet must have the same number of rows as the RotationSet.")
            mapped = np.einsum("nij,nj->ni", matrices, values, optimize=True)
            return VectorSet(
                values=mapped,
                reference_frame=vectors.reference_frame,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape == (3,):
            mapped = np.einsum("nij,j->ni", matrices, array, optimize=True)
        elif array.ndim == 2 and array.shape[1] == 3:
            if array.shape[0] != len(self):
                raise ValueError(
                    "Input vectors must have the same number of rows as the RotationSet."
                )
            mapped = np.einsum("nij,nj->ni", matrices, array, optimize=True)
        else:
            raise ValueError("Input vectors must have shape (3,) or (n, 3).")
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def apply_inverse(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        inverse_matrices = np.swapaxes(self.as_matrices(), -1, -2)
        if isinstance(vectors, VectorSet):
            values = vectors.values
            if values.shape[0] != len(self):
                raise ValueError("VectorSet must have the same number of rows as the RotationSet.")
            mapped = np.einsum("nij,nj->ni", inverse_matrices, values, optimize=True)
            return VectorSet(
                values=mapped,
                reference_frame=vectors.reference_frame,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape == (3,):
            mapped = np.einsum("nij,j->ni", inverse_matrices, array, optimize=True)
        elif array.ndim == 2 and array.shape[1] == 3:
            if array.shape[0] != len(self):
                raise ValueError(
                    "Input vectors must have the same number of rows as the RotationSet."
                )
            mapped = np.einsum("nij,nj->ni", inverse_matrices, array, optimize=True)
        else:
            raise ValueError("Input vectors must have shape (3,) or (n, 3).")
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def subset(self, indices: ArrayLike) -> RotationSet:
        return RotationSet(
            quaternions=self.quaternions[np.asarray(indices)],
            provenance=self.provenance,
        )
