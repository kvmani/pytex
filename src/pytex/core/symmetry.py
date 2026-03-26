from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    is_rotation_matrix,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.batches import VectorSet
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.provenance import ProvenanceRecord


def _rotation_matrix_from_axis_angle(axis: ArrayLike, angle_deg: float) -> np.ndarray:
    unit_axis = normalize_vector(axis)
    angle_rad = np.deg2rad(angle_deg)
    x, y, z = unit_axis
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    one_minus_c = 1.0 - c
    matrix = np.array(
        [
            [c + x * x * one_minus_c, x * y * one_minus_c - z * s, x * z * one_minus_c + y * s],
            [y * x * one_minus_c + z * s, c + y * y * one_minus_c, y * z * one_minus_c - x * s],
            [z * x * one_minus_c - y * s, z * y * one_minus_c + x * s, c + z * z * one_minus_c],
        ],
        dtype=np.float64,
    )
    return as_float_array(matrix, shape=(3, 3))


def _matrix_key(matrix: np.ndarray) -> tuple[float, ...]:
    rounded = np.round(matrix, decimals=8)
    return tuple(float(value) for value in rounded.ravel())


def _unique_rotation_matrices(matrices: list[np.ndarray]) -> np.ndarray:
    unique: dict[tuple[float, ...], np.ndarray] = {}
    for matrix in matrices:
        if not is_rotation_matrix(matrix):
            raise ValueError("Generated symmetry operator is not a proper rotation matrix.")
        unique[_matrix_key(matrix)] = as_float_array(matrix, shape=(3, 3))
    operators = np.stack(list(unique.values()), axis=0)
    operators = np.ascontiguousarray(operators)
    operators.setflags(write=False)
    return operators


def _group_from_generators(generators: list[np.ndarray]) -> np.ndarray:
    identity = np.eye(3, dtype=np.float64)
    known = {_matrix_key(identity): identity}
    frontier = [identity]
    all_generators = [as_float_array(generator, shape=(3, 3)) for generator in generators]
    while frontier:
        current = frontier.pop()
        for generator in all_generators:
            for candidate in (current @ generator, generator @ current):
                key = _matrix_key(candidate)
                if key not in known:
                    known[key] = candidate
                    frontier.append(candidate)
    return _unique_rotation_matrices(list(known.values()))


def _point_group_generators() -> dict[str, list[np.ndarray]]:
    return {
        "1": [],
        "2": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 180.0)],
        "222": [
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
            _rotation_matrix_from_axis_angle([0.0, 1.0, 0.0], 180.0),
        ],
        "4": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0)],
        "422": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "3": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 120.0)],
        "32": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 120.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "6": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 60.0)],
        "622": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 60.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "23": [
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
            _rotation_matrix_from_axis_angle([1.0, 1.0, 1.0], 120.0),
        ],
        "432": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0),
            _rotation_matrix_from_axis_angle([1.0, 1.0, 1.0], 120.0),
        ],
    }


_PROPER_POINT_GROUP_ALIASES = {
    "1": "1",
    "-1": "1",
    "2": "2",
    "2/m": "2",
    "222": "222",
    "mmm": "222",
    "4": "4",
    "4/m": "4",
    "422": "422",
    "4/mmm": "422",
    "3": "3",
    "-3": "3",
    "32": "32",
    "-3m": "32",
    "6": "6",
    "6/m": "6",
    "622": "622",
    "6/mmm": "622",
    "23": "23",
    "m-3": "23",
    "432": "432",
    "m-3m": "432",
}

_SECTOR_TOLERANCE = 1e-8


def _normalized_proper_point_group(point_group: str) -> str:
    normalized = point_group.replace(" ", "").lower()
    proper_group = _PROPER_POINT_GROUP_ALIASES.get(normalized)
    if proper_group is None:
        supported = ", ".join(sorted(_PROPER_POINT_GROUP_ALIASES))
        raise ValueError(f"Unsupported point group '{point_group}'. Supported groups: {supported}")
    return proper_group


def _angle_in_wedge(vector: np.ndarray, wedge_angle_rad: float) -> bool:
    angle = float(np.mod(np.arctan2(vector[1], vector[0]), 2.0 * np.pi))
    return angle <= wedge_angle_rad + _SECTOR_TOLERANCE


def _sector_vertices_for_group(proper_group: str) -> np.ndarray:
    if proper_group in {"23", "432"}:
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
            ],
            dtype=np.float64,
        )
    elif proper_group in {"4", "422"}:
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
    elif proper_group in {"3", "32"}:
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.5, np.sqrt(3.0) / 2.0, 0.0],
            ],
            dtype=np.float64,
        )
    elif proper_group in {"6", "622"}:
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [np.sqrt(3.0) / 2.0, 0.5, 0.0],
            ],
            dtype=np.float64,
        )
    else:
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
    return normalize_vectors(vertices)


def _vector_in_fundamental_sector(vector: np.ndarray, proper_group: str) -> bool:
    x, y, z = (float(value) for value in vector)
    if proper_group in {"23", "432"}:
        return (
            z >= -_SECTOR_TOLERANCE
            and y >= -_SECTOR_TOLERANCE
            and x >= y - _SECTOR_TOLERANCE
            and z >= x - _SECTOR_TOLERANCE
        )
    if proper_group in {"4", "422"}:
        return z >= -_SECTOR_TOLERANCE and y >= -_SECTOR_TOLERANCE and x >= y - _SECTOR_TOLERANCE
    if proper_group in {"3", "32"}:
        return (
            z >= -_SECTOR_TOLERANCE
            and x >= -_SECTOR_TOLERANCE
            and _angle_in_wedge(vector, np.deg2rad(60.0))
        )
    if proper_group in {"6", "622"}:
        return (
            z >= -_SECTOR_TOLERANCE
            and x >= -_SECTOR_TOLERANCE
            and _angle_in_wedge(vector, np.deg2rad(30.0))
        )
    if proper_group in {"2", "222"}:
        return x >= -_SECTOR_TOLERANCE and y >= -_SECTOR_TOLERANCE and z >= -_SECTOR_TOLERANCE
    return z >= -_SECTOR_TOLERANCE


def _sector_sort_key(vector: np.ndarray) -> tuple[float, float, float]:
    rounded = np.round(vector, decimals=12)
    return (
        float(rounded[2]),
        float(rounded[0]),
        float(rounded[1]),
    )


@cache
def _operators_for_proper_point_group(proper_group: str) -> np.ndarray:
    operators = _group_from_generators(_point_group_generators()[proper_group])
    operators.setflags(write=False)
    return operators


def _canonical_vector_index(vectors: np.ndarray) -> int:
    rounded = np.round(vectors, decimals=12)
    return int(np.lexsort((rounded[:, 1], rounded[:, 0], rounded[:, 2]))[-1])


@dataclass(frozen=True, slots=True)
class SymmetrySpec:
    name: str
    point_group: str
    operators: np.ndarray = field(default_factory=lambda: np.eye(3)[None, :, :])
    specimen_symmetry: str | None = None
    reference_frame: ReferenceFrame | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        operators = as_float_array(self.operators, shape=(None, 3, 3))
        for operator in operators:
            if not is_rotation_matrix(operator):
                raise ValueError("All symmetry operators must be proper rotation matrices.")
        if self.reference_frame is not None and self.reference_frame.domain not in {
            FrameDomain.CRYSTAL,
            FrameDomain.SPECIMEN,
        }:
            raise ValueError(
                "SymmetrySpec.reference_frame must belong to the crystal or specimen domain."
            )
        object.__setattr__(self, "operators", operators)

    @property
    def order(self) -> int:
        return int(self.operators.shape[0])

    @property
    def proper_point_group(self) -> str:
        return _normalized_proper_point_group(self.point_group)

    @classmethod
    def identity(
        cls,
        *,
        name: str = "identity",
        point_group: str = "1",
        reference_frame: ReferenceFrame | None = None,
    ) -> SymmetrySpec:
        return cls(
            name=name,
            point_group=point_group,
            operators=np.eye(3)[None, :, :],
            reference_frame=reference_frame,
        )

    @classmethod
    def from_point_group(
        cls,
        point_group: str,
        *,
        reference_frame: ReferenceFrame | None = None,
        specimen_symmetry: str | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> SymmetrySpec:
        proper_group = _normalized_proper_point_group(point_group)
        operators = _operators_for_proper_point_group(proper_group)
        return cls(
            name=proper_group,
            point_group=point_group,
            operators=operators,
            specimen_symmetry=specimen_symmetry,
            reference_frame=reference_frame,
            provenance=provenance,
        )

    def apply_to_vectors(self, vectors: ArrayLike | VectorSet) -> np.ndarray:
        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            array = vectors.values
        else:
            array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = np.einsum("oij,...j->o...i", self.operators, array, optimize=True)
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def apply_to_rotation_matrices(self, matrices: ArrayLike, *, side: str = "right") -> np.ndarray:
        rotations = np.asarray(matrices, dtype=np.float64)
        if rotations.shape[-2:] != (3, 3):
            raise ValueError("Input rotation matrices must have trailing shape (3, 3).")
        if side == "right":
            transformed = np.einsum("...ij,ojk->o...ik", rotations, self.operators, optimize=True)
        elif side == "left":
            transformed = np.einsum("oij,...jk->o...ik", self.operators, rotations, optimize=True)
        else:
            raise ValueError("side must be either 'left' or 'right'.")
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def equivalent_vectors(self, vector: ArrayLike, *, antipodal: bool = False) -> np.ndarray:
        candidates = self.apply_to_vectors(vector)
        candidates = normalize_vectors(candidates)
        if antipodal:
            combined = np.concatenate([candidates, -candidates], axis=0)
            candidates = normalize_vectors(combined)
        unique_vectors: dict[tuple[float, ...], np.ndarray] = {}
        for candidate in candidates:
            key = tuple(np.round(candidate, decimals=8))
            unique_vectors[key] = candidate
        array = np.stack(list(unique_vectors.values()), axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        return array

    def canonicalize_vector(self, vector: ArrayLike, *, antipodal: bool = False) -> np.ndarray:
        candidates = self.equivalent_vectors(vector, antipodal=antipodal)
        return as_float_array(candidates[_canonical_vector_index(candidates)], shape=(3,))

    def canonicalize_vectors(
        self,
        vectors: ArrayLike | VectorSet,
        *,
        antipodal: bool = False,
    ) -> np.ndarray | VectorSet:
        reference_frame = None
        provenance = None
        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            normalized = normalize_vectors(vectors.values)
            reference_frame = vectors.reference_frame
            provenance = vectors.provenance
        else:
            normalized = normalize_vectors(vectors)
        canonicalized = [
            self.canonicalize_vector(vector, antipodal=antipodal) for vector in normalized
        ]
        array = np.stack(canonicalized, axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        if reference_frame is not None:
            return VectorSet(
                values=array,
                reference_frame=reference_frame,
                provenance=provenance,
            )
        return array

    def fundamental_sector(self, *, antipodal: bool = True) -> FundamentalSector:
        return FundamentalSector(
            point_group=self.point_group,
            proper_point_group=self.proper_point_group,
            antipodal=antipodal,
            vertices=_sector_vertices_for_group(self.proper_point_group),
        )

    def vector_in_fundamental_sector(self, vector: ArrayLike, *, antipodal: bool = True) -> bool:
        candidate = normalize_vector(vector)
        if antipodal and candidate[2] < 0.0:
            candidate = -candidate
        return _vector_in_fundamental_sector(candidate, self.proper_point_group)

    def reduce_vector_to_fundamental_sector(
        self,
        vector: ArrayLike,
        *,
        antipodal: bool = True,
    ) -> np.ndarray:
        candidates = self.equivalent_vectors(vector, antipodal=antipodal)
        matching = [
            candidate
            for candidate in candidates
            if _vector_in_fundamental_sector(candidate, self.proper_point_group)
        ]
        if matching:
            selected = max(matching, key=_sector_sort_key)
            return as_float_array(selected, shape=(3,))
        return self.canonicalize_vector(vector, antipodal=antipodal)

    def reduce_vectors_to_fundamental_sector(
        self,
        vectors: ArrayLike | VectorSet,
        *,
        antipodal: bool = True,
    ) -> np.ndarray | VectorSet:
        reference_frame = None
        provenance = None
        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            normalized = normalize_vectors(vectors.values)
            reference_frame = vectors.reference_frame
            provenance = vectors.provenance
        else:
            normalized = normalize_vectors(vectors)
        reduced = [
            self.reduce_vector_to_fundamental_sector(vector, antipodal=antipodal)
            for vector in normalized
        ]
        array = np.stack(reduced, axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        if reference_frame is not None:
            return VectorSet(
                values=array,
                reference_frame=reference_frame,
                provenance=provenance,
            )
        return array


@dataclass(frozen=True, slots=True)
class FundamentalSector:
    point_group: str
    proper_point_group: str
    antipodal: bool
    vertices: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", normalize_vectors(self.vertices))
