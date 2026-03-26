from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    normalize_quaternion,
    normalize_quaternions,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec

_EULER_CONVENTION_ALIASES = {
    "bunge": "bunge",
    "bunge_zxz": "bunge",
    "zxz": "bunge",
    "matthies": "matthies",
    "matthies_zyz": "matthies",
    "abg": "abg",
    "abg_zyz": "abg",
    "zyz": "abg",
}


def _normalize_euler_convention(convention: str) -> str:
    normalized = convention.strip().lower()
    resolved = _EULER_CONVENTION_ALIASES.get(normalized)
    if resolved is None:
        supported = ", ".join(sorted(_EULER_CONVENTION_ALIASES))
        raise ValueError(
            f"Unsupported Euler convention '{convention}'. Supported conventions: {supported}"
        )
    return resolved


def _axis_angle_quaternion_for_axis(axis_name: str, angle_rad: float) -> np.ndarray:
    axis_map = {
        "x": [1.0, 0.0, 0.0],
        "y": [0.0, 1.0, 0.0],
        "z": [0.0, 0.0, 1.0],
    }
    return quaternion_from_axis_angle(axis_map[axis_name], angle_rad)


def _euler_axes_for_convention(convention: str) -> tuple[str, str, str]:
    normalized = _normalize_euler_convention(convention)
    if normalized == "bunge":
        return ("z", "x", "z")
    return ("z", "y", "z")


def _matrix_to_repeated_axis_euler(
    matrix: np.ndarray,
    *,
    convention: str,
) -> tuple[float, float, float]:
    normalized = _normalize_euler_convention(convention)
    phi_rad = float(_safe_arccos(matrix[2, 2]))
    if normalized == "bunge":
        if np.isclose(phi_rad, 0.0, atol=1e-10):
            first = float(np.arctan2(matrix[1, 0], matrix[0, 0]))
            third = 0.0
        elif np.isclose(phi_rad, np.pi, atol=1e-10):
            first = float(np.arctan2(matrix[0, 1], matrix[0, 0]))
            third = 0.0
        else:
            first = float(np.arctan2(matrix[0, 2], -matrix[1, 2]))
            third = float(np.arctan2(matrix[2, 0], matrix[2, 1]))
    else:
        if np.isclose(phi_rad, 0.0, atol=1e-10) or np.isclose(phi_rad, np.pi, atol=1e-10):
            first = float(np.arctan2(matrix[1, 0], matrix[0, 0]))
            third = 0.0
        else:
            first = float(np.arctan2(matrix[1, 2], matrix[0, 2]))
            third = float(np.arctan2(matrix[2, 1], -matrix[2, 0]))
    return (first, phi_rad, third)


def _safe_arccos(value: ArrayLike) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    return np.arccos(np.clip(array, -1.0, 1.0))


def _canonicalize_quaternion(quaternion: ArrayLike) -> np.ndarray:
    candidate = normalize_quaternion(quaternion)
    if candidate[0] < 0.0:
        candidate = -candidate
    return normalize_quaternion(candidate)


def _canonical_quaternion_index(quaternions: np.ndarray) -> int:
    candidates = np.asarray(quaternions, dtype=np.float64)
    canonical = np.stack([_canonicalize_quaternion(candidate) for candidate in candidates], axis=0)
    rounded = np.round(canonical, decimals=12)
    return int(np.lexsort((rounded[:, 3], rounded[:, 2], rounded[:, 1], rounded[:, 0]))[-1])


def _exact_fundamental_region_key_from_quaternion(
    quaternion: ArrayLike,
) -> tuple[float, float, float, float]:
    canonical = np.round(_canonicalize_quaternion(quaternion), decimals=12)
    return (
        -float(canonical[0]),
        -float(canonical[1]),
        -float(canonical[2]),
        -float(canonical[3]),
    )


def _fundamental_region_key(
    rotation: Rotation, symmetry: SymmetrySpec | None
) -> tuple[float, float, float, float]:
    del symmetry
    return _exact_fundamental_region_key_from_quaternion(rotation.quaternion)


def quaternion_multiply(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    lw, lx, ly, lz = left
    rw, rx, ry, rz = right
    return np.array(
        [
            lw * rw - lx * rx - ly * ry - lz * rz,
            lw * rx + lx * rw + ly * rz - lz * ry,
            lw * ry - lx * rz + ly * rw + lz * rx,
            lw * rz + lx * ry - ly * rx + lz * rw,
        ],
        dtype=np.float64,
    )


def quaternion_conjugate(quaternion: np.ndarray) -> np.ndarray:
    w, x, y, z = quaternion
    return np.array([w, -x, -y, -z], dtype=np.float64)


def quaternion_to_matrix(quaternion: np.ndarray) -> np.ndarray:
    w, x, y, z = normalize_quaternion(quaternion)
    return as_float_array(
        [
            [1 - 2 * (y**2 + z**2), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x**2 + z**2), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x**2 + y**2)],
        ],
        shape=(3, 3),
    )


def matrix_to_quaternion(matrix: ArrayLike) -> np.ndarray:
    array = as_float_array(matrix, shape=(3, 3))
    trace = float(np.trace(array))
    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (array[2, 1] - array[1, 2]) * s
        y = (array[0, 2] - array[2, 0]) * s
        z = (array[1, 0] - array[0, 1]) * s
    else:
        diagonal = np.diag(array)
        index = int(np.argmax(diagonal))
        if index == 0:
            s = 2.0 * np.sqrt(1.0 + array[0, 0] - array[1, 1] - array[2, 2])
            w = (array[2, 1] - array[1, 2]) / s
            x = 0.25 * s
            y = (array[0, 1] + array[1, 0]) / s
            z = (array[0, 2] + array[2, 0]) / s
        elif index == 1:
            s = 2.0 * np.sqrt(1.0 + array[1, 1] - array[0, 0] - array[2, 2])
            w = (array[0, 2] - array[2, 0]) / s
            x = (array[0, 1] + array[1, 0]) / s
            y = 0.25 * s
            z = (array[1, 2] + array[2, 1]) / s
        else:
            s = 2.0 * np.sqrt(1.0 + array[2, 2] - array[0, 0] - array[1, 1])
            w = (array[1, 0] - array[0, 1]) / s
            x = (array[0, 2] + array[2, 0]) / s
            y = (array[1, 2] + array[2, 1]) / s
            z = 0.25 * s
    return normalize_quaternion([w, x, y, z])


def quaternion_from_axis_angle(axis: ArrayLike, angle_rad: float) -> np.ndarray:
    unit_axis = normalize_vector(axis)
    half = angle_rad / 2.0
    sin_half = np.sin(half)
    return normalize_quaternion([np.cos(half), *(unit_axis * sin_half)])


@dataclass(frozen=True, slots=True)
class Rotation:
    quaternion: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternion", normalize_quaternion(self.quaternion))

    @classmethod
    def identity(cls) -> Rotation:
        return cls(quaternion=np.array([1.0, 0.0, 0.0, 0.0]))

    @classmethod
    def from_matrix(cls, matrix: ArrayLike) -> Rotation:
        return cls(quaternion=matrix_to_quaternion(matrix))

    @classmethod
    def from_axis_angle(cls, axis: ArrayLike, angle_rad: float) -> Rotation:
        return cls(quaternion=quaternion_from_axis_angle(axis, angle_rad))

    @classmethod
    def from_euler(
        cls,
        angle1: float,
        angle2: float,
        angle3: float,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> Rotation:
        angles = np.array([angle1, angle2, angle3], dtype=np.float64)
        if degrees:
            angles = np.deg2rad(angles)
        axes = _euler_axes_for_convention(convention)
        q1 = _axis_angle_quaternion_for_axis(axes[0], float(angles[0]))
        q2 = _axis_angle_quaternion_for_axis(axes[1], float(angles[1]))
        q3 = _axis_angle_quaternion_for_axis(axes[2], float(angles[2]))
        return cls(quaternion=quaternion_multiply(quaternion_multiply(q1, q2), q3))

    @classmethod
    def from_bunge_euler(
        cls,
        phi1: float,
        Phi: float,  # noqa: N803 - crystallographic notation is intentional
        phi2: float,
        *,
        degrees: bool = True,
    ) -> Rotation:
        return cls.from_euler(phi1, Phi, phi2, convention="bunge", degrees=degrees)

    @classmethod
    def from_matthies_euler(
        cls,
        alpha: float,
        beta: float,
        gamma: float,
        *,
        degrees: bool = True,
    ) -> Rotation:
        return cls.from_euler(alpha, beta, gamma, convention="matthies", degrees=degrees)

    @classmethod
    def from_abg_euler(
        cls,
        alpha: float,
        beta: float,
        gamma: float,
        *,
        degrees: bool = True,
    ) -> Rotation:
        return cls.from_euler(alpha, beta, gamma, convention="abg", degrees=degrees)

    def as_matrix(self) -> np.ndarray:
        return quaternion_to_matrix(self.quaternion)

    def to_euler(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> tuple[float, float, float]:
        matrix = self.as_matrix()
        angles = np.mod(_matrix_to_repeated_axis_euler(matrix, convention=convention), 2.0 * np.pi)
        if degrees:
            angles_deg = np.rad2deg(angles)
            return (
                float(angles_deg[0]),
                float(angles_deg[1]),
                float(angles_deg[2]),
            )
        return (float(angles[0]), float(angles[1]), float(angles[2]))

    def to_bunge_euler(self, *, degrees: bool = True) -> tuple[float, float, float]:
        return self.to_euler(convention="bunge", degrees=degrees)

    def to_matthies_euler(self, *, degrees: bool = True) -> tuple[float, float, float]:
        return self.to_euler(convention="matthies", degrees=degrees)

    def to_abg_euler(self, *, degrees: bool = True) -> tuple[float, float, float]:
        return self.to_euler(convention="abg", degrees=degrees)

    def compose(self, other: Rotation) -> Rotation:
        return Rotation(quaternion=quaternion_multiply(self.quaternion, other.quaternion))

    def inverse(self) -> Rotation:
        return Rotation(quaternion=quaternion_conjugate(self.quaternion))

    def canonicalized(self) -> Rotation:
        return Rotation(quaternion=_canonicalize_quaternion(self.quaternion))

    @property
    def angle_rad(self) -> float:
        return float(2.0 * _safe_arccos(abs(float(self.quaternion[0]))))

    @property
    def angle_deg(self) -> float:
        return float(np.rad2deg(self.angle_rad))

    @property
    def axis(self) -> np.ndarray:
        scalar = float(np.clip(self.quaternion[0], -1.0, 1.0))
        sin_half = np.sqrt(max(0.0, 1.0 - scalar * scalar))
        if np.isclose(sin_half, 0.0):
            return as_float_array([0.0, 0.0, 1.0], shape=(3,))
        return normalize_vector(self.quaternion[1:] / sin_half)

    def distance_to(self, other: Rotation) -> float:
        return other.compose(self.inverse()).angle_rad

    def apply(self, vectors: ArrayLike) -> np.ndarray:
        matrix = self.as_matrix()
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = array @ matrix.T
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def apply_inverse(self, vectors: ArrayLike) -> np.ndarray:
        return self.inverse().apply(vectors)


@dataclass(frozen=True, slots=True)
class Orientation:
    rotation: Rotation
    crystal_frame: ReferenceFrame
    specimen_frame: ReferenceFrame
    symmetry: SymmetrySpec | None = None
    phase: Phase | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.crystal_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("Orientation.crystal_frame must belong to the crystal domain.")
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError("Orientation.specimen_frame must belong to the specimen domain.")
        if self.symmetry is not None and self.symmetry.reference_frame != self.crystal_frame:
            raise ValueError(
                "Orientation.symmetry.reference_frame must match Orientation.crystal_frame."
            )
        if self.phase is not None:
            if self.phase.crystal_frame != self.crystal_frame:
                raise ValueError("Orientation.phase.crystal_frame must match crystal_frame.")
            if self.symmetry is not None and self.phase.symmetry != self.symmetry:
                raise ValueError("Orientation.phase.symmetry must match Orientation.symmetry.")

    def as_matrix(self) -> np.ndarray:
        return self.rotation.as_matrix()

    def map_crystal_vector(self, vector: ArrayLike) -> np.ndarray:
        return as_float_array(self.rotation.apply(vector), shape=(3,))

    def map_sample_vector_to_crystal(self, vector: ArrayLike) -> np.ndarray:
        return as_float_array(self.rotation.apply_inverse(vector), shape=(3,))

    def equivalent_orientations(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> OrientationSet:
        if (
            specimen_symmetry is not None
            and specimen_symmetry.reference_frame != self.specimen_frame
        ):
            raise ValueError(
                "specimen_symmetry.reference_frame must match Orientation.specimen_frame."
            )
        left_operators = (
            specimen_symmetry.operators
            if specimen_symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        right_operators = (
            self.symmetry.operators
            if self.symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        base = self.as_matrix()
        matrices = np.stack(
            [
                left_operator @ base @ right_operator
                for left_operator in left_operators
                for right_operator in right_operators
            ],
            axis=0,
        )
        quaternions = np.stack([matrix_to_quaternion(matrix) for matrix in matrices], axis=0)
        return OrientationSet(
            quaternions=quaternions,
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def canonicalize(self, specimen_symmetry: SymmetrySpec | None = None) -> Orientation:
        equivalents = self.equivalent_orientations(specimen_symmetry=specimen_symmetry)
        index = _canonical_quaternion_index(equivalents.quaternions)
        return Orientation(
            rotation=Rotation(equivalents.quaternions[index]),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def fundamental_region_key(
        self,
        *,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> tuple[float, float, float, float]:
        projected = self.project_to_exact_fundamental_region(specimen_symmetry=specimen_symmetry)
        return _fundamental_region_key(projected.rotation, projected.symmetry)

    def project_to_fundamental_region(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
        reference_orientation: Orientation | None = None,
    ) -> Orientation:
        return self.project_to_exact_fundamental_region(
            specimen_symmetry=specimen_symmetry,
            reference_orientation=reference_orientation,
        )

    def project_to_exact_fundamental_region(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
        reference_orientation: Orientation | None = None,
    ) -> Orientation:
        equivalents = self.equivalent_orientations(specimen_symmetry=specimen_symmetry)
        if reference_orientation is None:
            keys = [
                _fundamental_region_key(Rotation(quaternion), self.symmetry)
                for quaternion in equivalents.quaternions
            ]
            index = int(min(range(len(keys)), key=keys.__getitem__))
        else:
            if reference_orientation.crystal_frame != self.crystal_frame:
                raise ValueError(
                    "reference_orientation.crystal_frame must match Orientation.crystal_frame."
                )
            if reference_orientation.specimen_frame != self.specimen_frame:
                raise ValueError(
                    "reference_orientation.specimen_frame must match Orientation.specimen_frame."
                )
            index = int(
                np.argmin(
                    [
                        Orientation(
                            rotation=Rotation(quaternion),
                            crystal_frame=self.crystal_frame,
                            specimen_frame=self.specimen_frame,
                            symmetry=self.symmetry,
                            phase=self.phase,
                            provenance=self.provenance,
                        ).distance_to(reference_orientation, symmetry_aware=False)
                        for quaternion in equivalents.quaternions
                    ]
                )
            )
        return Orientation(
            rotation=Rotation(equivalents.quaternions[index]),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def is_in_fundamental_region(
        self,
        *,
        specimen_symmetry: SymmetrySpec | None = None,
        atol: float = 1e-10,
    ) -> bool:
        projected = self.project_to_exact_fundamental_region(specimen_symmetry=specimen_symmetry)
        return bool(
            np.allclose(
                _canonicalize_quaternion(self.rotation.quaternion),
                _canonicalize_quaternion(projected.rotation.quaternion),
                atol=atol,
            )
        )

    def misorientation_to(
        self,
        other: Orientation,
        *,
        reduce_by_symmetry: bool = True,
    ) -> Misorientation:
        if self.crystal_frame != other.crystal_frame:
            raise ValueError("Misorientation requires the same crystal frame.")
        if self.specimen_frame != other.specimen_frame:
            raise ValueError("Misorientation requires the same specimen frame.")
        if self.phase is not None and other.phase is not None and self.phase != other.phase:
            raise ValueError("Misorientation requires matching phases when both are specified.")
        delta = other.rotation.compose(self.rotation.inverse())
        misorientation = Misorientation(
            rotation=delta,
            left_symmetry=self.symmetry,
            right_symmetry=other.symmetry,
            provenance=self.provenance or other.provenance,
        )
        if reduce_by_symmetry:
            return misorientation.disorientation()
        return misorientation

    def distance_to(self, other: Orientation, *, symmetry_aware: bool = True) -> float:
        return self.misorientation_to(other, reduce_by_symmetry=symmetry_aware).angle_rad


@dataclass(frozen=True, slots=True)
class Misorientation:
    rotation: Rotation
    left_symmetry: SymmetrySpec | None = None
    right_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def as_matrix(self) -> np.ndarray:
        return self.rotation.as_matrix()

    @property
    def angle_rad(self) -> float:
        return self.rotation.angle_rad

    @property
    def angle_deg(self) -> float:
        return self.rotation.angle_deg

    def disorientation(self) -> Misorientation:
        left_operators = (
            self.left_symmetry.operators
            if self.left_symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        right_operators = (
            self.right_symmetry.operators
            if self.right_symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        base = self.as_matrix()
        candidates = np.stack(
            [
                left_operator @ base @ right_operator.T
                for left_operator in left_operators
                for right_operator in right_operators
            ],
            axis=0,
        )
        quaternions = np.stack(
            [matrix_to_quaternion(candidate) for candidate in candidates],
            axis=0,
        )
        keys = [
            _exact_fundamental_region_key_from_quaternion(quaternion) for quaternion in quaternions
        ]
        index = int(min(range(len(keys)), key=keys.__getitem__))
        return Misorientation(
            rotation=Rotation(quaternions[index]),
            left_symmetry=self.left_symmetry,
            right_symmetry=self.right_symmetry,
            provenance=self.provenance,
        )


@dataclass(frozen=True, slots=True)
class OrientationSet:
    quaternions: np.ndarray
    crystal_frame: ReferenceFrame
    specimen_frame: ReferenceFrame
    symmetry: SymmetrySpec | None = None
    phase: Phase | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternions", normalize_quaternions(self.quaternions))
        if self.crystal_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("OrientationSet.crystal_frame must belong to the crystal domain.")
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError("OrientationSet.specimen_frame must belong to the specimen domain.")
        if self.symmetry is not None and self.symmetry.reference_frame != self.crystal_frame:
            raise ValueError(
                "OrientationSet.symmetry.reference_frame must match OrientationSet.crystal_frame."
            )
        if self.phase is not None:
            if self.phase.crystal_frame != self.crystal_frame:
                raise ValueError("OrientationSet.phase.crystal_frame must match crystal_frame.")
            if self.symmetry is not None and self.phase.symmetry != self.symmetry:
                raise ValueError(
                    "OrientationSet.phase.symmetry must match OrientationSet.symmetry."
                )

    @classmethod
    def from_orientations(cls, orientations: list[Orientation]) -> OrientationSet:
        if not orientations:
            raise ValueError("OrientationSet requires at least one orientation.")
        crystal_frame = orientations[0].crystal_frame
        specimen_frame = orientations[0].specimen_frame
        symmetry = orientations[0].symmetry
        phase = orientations[0].phase
        provenance = orientations[0].provenance
        for orientation in orientations[1:]:
            if orientation.crystal_frame != crystal_frame:
                raise ValueError(
                    "All orientations in an OrientationSet must share a crystal frame."
                )
            if orientation.specimen_frame != specimen_frame:
                raise ValueError(
                    "All orientations in an OrientationSet must share a specimen frame."
                )
            if orientation.symmetry != symmetry:
                raise ValueError("All orientations in an OrientationSet must share symmetry.")
            if orientation.phase != phase:
                raise ValueError("All orientations in an OrientationSet must share phase.")
            if orientation.provenance != provenance:
                raise ValueError(
                    "All orientations in an OrientationSet must share provenance until "
                    "aggregate provenance records are implemented."
                )
        quaternions = np.vstack([orientation.rotation.quaternion for orientation in orientations])
        return cls(
            quaternions=quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    def __getitem__(self, index: int) -> Orientation:
        quaternion = as_float_array(self.quaternions[index], shape=(4,))
        return Orientation(
            rotation=Rotation(quaternion=quaternion),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    @classmethod
    def from_euler_angles(
        cls,
        angles: ArrayLike,
        *,
        crystal_frame: ReferenceFrame,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        convention: str = "bunge",
        degrees: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        angle_array = as_float_array(angles, shape=(None, 3))
        quaternions = np.stack(
            [
                Rotation.from_euler(
                    angle1,
                    angle2,
                    angle3,
                    convention=convention,
                    degrees=degrees,
                ).quaternion
                for angle1, angle2, angle3 in angle_array
            ],
            axis=0,
        )
        return cls(
            quaternions=quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    def as_matrices(self) -> np.ndarray:
        matrices = np.stack(
            [quaternion_to_matrix(quaternion) for quaternion in self.quaternions],
            axis=0,
        )
        matrices = np.ascontiguousarray(matrices)
        matrices.setflags(write=False)
        return matrices

    def as_euler(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> np.ndarray:
        euler = np.stack(
            [
                Rotation(quaternion).to_euler(convention=convention, degrees=degrees)
                for quaternion in self.quaternions
            ],
            axis=0,
        )
        euler = np.ascontiguousarray(euler)
        euler.setflags(write=False)
        return euler

    def as_bunge_euler(self, *, degrees: bool = True) -> np.ndarray:
        return self.as_euler(convention="bunge", degrees=degrees)

    def map_crystal_directions(self, directions: ArrayLike) -> np.ndarray:
        matrices = self.as_matrices()
        direction_array = np.asarray(directions, dtype=np.float64)
        if direction_array.shape == (3,):
            mapped = np.einsum("nij,j->ni", matrices, direction_array, optimize=True)
        elif direction_array.ndim == 2 and direction_array.shape[1] == 3:
            if direction_array.shape[0] != len(self):
                raise ValueError(
                    "Direction array must have the same number of rows as the OrientationSet."
                )
            mapped = np.einsum("nij,nj->ni", matrices, direction_array, optimize=True)
        else:
            raise ValueError("Directions must have shape (3,) or (n, 3).")
        mapped = normalize_vectors(mapped)
        return mapped

    def map_sample_directions_to_crystal(self, directions: ArrayLike) -> np.ndarray:
        inverse_matrices = np.swapaxes(self.as_matrices(), -1, -2)
        direction_array = np.asarray(directions, dtype=np.float64)
        if direction_array.shape == (3,):
            mapped = np.einsum("nij,j->ni", inverse_matrices, direction_array, optimize=True)
        elif direction_array.ndim == 2 and direction_array.shape[1] == 3:
            if direction_array.shape[0] != len(self):
                raise ValueError(
                    "Direction array must have the same number of rows as the OrientationSet."
                )
            mapped = np.einsum("nij,nj->ni", inverse_matrices, direction_array, optimize=True)
        else:
            raise ValueError("Directions must have shape (3,) or (n, 3).")
        mapped = normalize_vectors(mapped)
        return mapped

    def canonicalized(self, specimen_symmetry: SymmetrySpec | None = None) -> OrientationSet:
        canonical_quaternions = [
            Orientation(
                rotation=Rotation(quaternion),
                crystal_frame=self.crystal_frame,
                specimen_frame=self.specimen_frame,
                symmetry=self.symmetry,
                phase=self.phase,
                provenance=self.provenance,
            )
            .canonicalize(specimen_symmetry=specimen_symmetry)
            .rotation.quaternion
            for quaternion in self.quaternions
        ]
        return OrientationSet(
            quaternions=np.stack(canonical_quaternions, axis=0),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def projected_to_fundamental_region(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
        reference_orientation: Orientation | None = None,
    ) -> OrientationSet:
        projected_quaternions = [
            Orientation(
                rotation=Rotation(quaternion),
                crystal_frame=self.crystal_frame,
                specimen_frame=self.specimen_frame,
                symmetry=self.symmetry,
                phase=self.phase,
                provenance=self.provenance,
            )
            .project_to_fundamental_region(
                specimen_symmetry=specimen_symmetry,
                reference_orientation=reference_orientation,
            )
            .rotation.quaternion
            for quaternion in self.quaternions
        ]
        return OrientationSet(
            quaternions=np.stack(projected_quaternions, axis=0),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def exact_fundamental_region_keys(
        self,
        *,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> np.ndarray:
        keys = np.asarray(
            [
                Orientation(
                    rotation=Rotation(quaternion),
                    crystal_frame=self.crystal_frame,
                    specimen_frame=self.specimen_frame,
                    symmetry=self.symmetry,
                    phase=self.phase,
                    provenance=self.provenance,
                ).fundamental_region_key(specimen_symmetry=specimen_symmetry)
                for quaternion in self.quaternions
            ],
            dtype=np.float64,
        )
        keys = np.ascontiguousarray(keys)
        keys.setflags(write=False)
        return keys

    def fundamental_region_keys(
        self,
        *,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> np.ndarray:
        return self.exact_fundamental_region_keys(specimen_symmetry=specimen_symmetry)

    def misorientation_angles_to(
        self,
        other: OrientationSet,
        *,
        symmetry_aware: bool = True,
    ) -> np.ndarray:
        angles = np.empty((len(self), len(other)), dtype=np.float64)
        for row, quaternion_a in enumerate(self.quaternions):
            orientation_a = Orientation(
                rotation=Rotation(quaternion_a),
                crystal_frame=self.crystal_frame,
                specimen_frame=self.specimen_frame,
                symmetry=self.symmetry,
                phase=self.phase,
                provenance=self.provenance,
            )
            for column, quaternion_b in enumerate(other.quaternions):
                orientation_b = Orientation(
                    rotation=Rotation(quaternion_b),
                    crystal_frame=other.crystal_frame,
                    specimen_frame=other.specimen_frame,
                    symmetry=other.symmetry,
                    phase=other.phase,
                    provenance=other.provenance,
                )
                angles[row, column] = orientation_a.distance_to(
                    orientation_b,
                    symmetry_aware=symmetry_aware,
                )
        angles = np.ascontiguousarray(angles)
        angles.setflags(write=False)
        return angles

    def subset(self, indices: ArrayLike) -> OrientationSet:
        return OrientationSet(
            quaternions=self.quaternions[np.asarray(indices)],
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )
