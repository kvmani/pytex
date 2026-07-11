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
from pytex.core.point_groups import (
    PointGroup,
    laue_class_symbol_for,
    normalize_point_group_symbol,
    proper_subgroup_symbol_for,
)
from pytex.core.provenance import ProvenanceRecord

_SPECIMEN_SYMMETRY_POINT_GROUPS = {
    "triclinic": "1",
    "monoclinic": "2",
    "orthorhombic": "222",
    "orthotropic": "222",
}


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


_SECTOR_TOLERANCE = 1e-8


def _normalized_proper_point_group(point_group: str) -> str:
    return proper_subgroup_symbol_for(normalize_point_group_symbol(point_group))


def _sector_array(rows: list[list[float]]) -> np.ndarray:
    if not rows:
        return as_float_array(np.zeros((0, 3), dtype=np.float64), shape=(None, 3))
    return normalize_vectors(np.array(rows, dtype=np.float64))


@cache
def _sector_geometry_table() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    # Exact Laue-reduced (antipodal, upper hemisphere) sector polygons per
    # proper rotation group. Vertices are ordered loops; edge i is the great
    # circle arc from vertex i to vertex i+1 with the given inward normal.
    # Sector area is 4*pi / (2 * group order).
    half = 0.5
    root3_half = np.sqrt(3.0) / 2.0
    z_axis = [0.0, 0.0, 1.0]
    x_axis = [1.0, 0.0, 0.0]
    y_axis = [0.0, 1.0, 0.0]
    y_normal = [0.0, 1.0, 0.0]
    octant_vertices = _sector_array([z_axis, x_axis, [0.0, 1.0, 0.0]])
    octant_normals = _sector_array([y_normal, z_axis, x_axis])
    return {
        "1": (_sector_array([]), _sector_array([z_axis])),
        "2": (
            _sector_array([x_axis, [-1.0, 0.0, 0.0]]),
            _sector_array([z_axis, y_normal]),
        ),
        "222": (octant_vertices, octant_normals),
        "4": (octant_vertices, octant_normals),
        "422": (
            _sector_array([z_axis, x_axis, [1.0, 1.0, 0.0]]),
            _sector_array([y_normal, z_axis, [1.0, -1.0, 0.0]]),
        ),
        "3": (
            _sector_array([z_axis, x_axis, [-half, root3_half, 0.0]]),
            _sector_array([y_normal, z_axis, [root3_half, half, 0.0]]),
        ),
        # For D3 the antipodal-composed two-fold axes at azimuths 0/120/240
        # induce mirror lines at azimuths 30/90/150, so the fundamental wedge
        # is [30, 90] degrees rather than a wedge anchored at azimuth 0.
        "32": (
            _sector_array([z_axis, [root3_half, half, 0.0], y_axis]),
            _sector_array([[-half, root3_half, 0.0], z_axis, x_axis]),
        ),
        "6": (
            _sector_array([z_axis, x_axis, [half, root3_half, 0.0]]),
            _sector_array([y_normal, z_axis, [root3_half, -half, 0.0]]),
        ),
        "622": (
            _sector_array([z_axis, x_axis, [root3_half, half, 0.0]]),
            _sector_array([y_normal, z_axis, [half, -root3_half, 0.0]]),
        ),
        "23": (
            _sector_array([z_axis, x_axis, [1.0, 1.0, 1.0]]),
            _sector_array([y_normal, [0.0, -1.0, 1.0], [1.0, -1.0, 0.0]]),
        ),
        "432": (
            _sector_array([z_axis, [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]),
            _sector_array([y_normal, [-1.0, 0.0, 1.0], [1.0, -1.0, 0.0]]),
        ),
    }


def _sector_vertices_for_group(proper_group: str) -> np.ndarray:
    return _sector_geometry_table()[proper_group][0]


def _sector_edge_normals_for_group(proper_group: str) -> np.ndarray:
    return _sector_geometry_table()[proper_group][1]


def _vector_in_fundamental_sector(vector: np.ndarray, proper_group: str) -> bool:
    normals = _sector_edge_normals_for_group(proper_group)
    return bool(np.all(normals @ np.asarray(vector, dtype=np.float64) >= -_SECTOR_TOLERANCE))


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

    @property
    def laue_group_symbol(self) -> str:
        return laue_class_symbol_for(self.point_group)

    @property
    def is_laue(self) -> bool:
        return normalize_point_group_symbol(self.point_group) == self.laue_group_symbol

    def to_point_group(self) -> PointGroup:
        return PointGroup.from_symbol(self.point_group)

    def laue_symmetry(self) -> SymmetrySpec:
        return SymmetrySpec.from_point_group(
            self.laue_group_symbol,
            reference_frame=self.reference_frame,
            specimen_symmetry=self.specimen_symmetry,
            provenance=self.provenance,
        )

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

    @classmethod
    def specimen(
        cls,
        name: str = "triclinic",
        *,
        reference_frame: ReferenceFrame | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> SymmetrySpec:
        normalized = name.strip().lower()
        point_group = _SPECIMEN_SYMMETRY_POINT_GROUPS.get(normalized)
        if point_group is None:
            supported = ", ".join(sorted(_SPECIMEN_SYMMETRY_POINT_GROUPS))
            raise ValueError(
                f"Unsupported specimen symmetry '{name}'. Supported names: {supported}."
            )
        if reference_frame is not None and reference_frame.domain != FrameDomain.SPECIMEN:
            raise ValueError(
                "SymmetrySpec.specimen reference_frame must belong to the specimen domain."
            )
        spec = cls.from_point_group(
            point_group,
            reference_frame=reference_frame,
            specimen_symmetry=normalized,
            provenance=provenance,
        )
        return cls(
            name=f"specimen-{normalized}",
            point_group=spec.point_group,
            operators=spec.operators,
            specimen_symmetry=normalized,
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
            edge_normals=_sector_edge_normals_for_group(self.proper_point_group),
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
    edge_normals: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    def __post_init__(self) -> None:
        vertices = as_float_array(self.vertices, shape=(None, 3))
        if vertices.shape[0] > 0:
            vertices = normalize_vectors(vertices)
        edge_normals = as_float_array(self.edge_normals, shape=(None, 3))
        if edge_normals.shape[0] > 0:
            edge_normals = normalize_vectors(edge_normals)
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "edge_normals", edge_normals)

    def contains(self, vectors: ArrayLike, *, atol: float = _SECTOR_TOLERANCE) -> np.ndarray:
        array = np.asarray(vectors, dtype=np.float64)
        squeeze = array.ndim == 1
        rows = np.atleast_2d(array)
        if rows.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        if self.edge_normals.shape[0] == 0:
            mask = np.ones(rows.shape[0], dtype=bool)
        else:
            mask = np.all(rows @ self.edge_normals.T >= -atol, axis=1)
        if squeeze:
            return np.asarray(mask[0])
        mask = np.ascontiguousarray(mask)
        mask.setflags(write=False)
        return mask

    def center(self) -> np.ndarray:
        if self.vertices.shape[0] == 0:
            return as_float_array([0.0, 0.0, 1.0], shape=(3,))
        if self.vertices.shape[0] >= 3:
            return normalize_vector(self.vertices.sum(axis=0))
        trace = self.boundary_trace()
        return normalize_vector(trace.sum(axis=0))

    def boundary_trace(self, *, samples_per_edge: int = 64) -> np.ndarray:
        if samples_per_edge < 2:
            raise ValueError("boundary_trace requires at least two samples per edge.")
        if self.vertices.shape[0] == 0:
            angles = np.linspace(0.0, 2.0 * np.pi, 4 * samples_per_edge)
            trace = np.column_stack(
                [np.cos(angles), np.sin(angles), np.zeros_like(angles)]
            )
            return as_float_array(np.ascontiguousarray(trace), shape=(None, 3))
        if self.edge_normals.shape[0] != self.vertices.shape[0]:
            raise ValueError(
                "boundary_trace requires one edge normal per vertex of the sector loop."
            )
        segments: list[np.ndarray] = []
        vertex_count = self.vertices.shape[0]
        for index in range(vertex_count):
            start = self.vertices[index]
            end = self.vertices[(index + 1) % vertex_count]
            normal = self.edge_normals[index]
            tangent = np.cross(normal, start)
            tangent_norm = float(np.linalg.norm(tangent))
            if np.isclose(tangent_norm, 0.0):
                raise ValueError(
                    "Sector edge normal must be perpendicular to its starting vertex."
                )
            tangent = tangent / tangent_norm
            sweep = float(
                np.mod(np.arctan2(float(tangent @ end), float(start @ end)), 2.0 * np.pi)
            )
            angles = np.linspace(0.0, sweep, samples_per_edge, endpoint=False)
            segments.append(
                np.outer(np.cos(angles), start) + np.outer(np.sin(angles), tangent)
            )
        segments.append(self.vertices[:1])
        trace = np.concatenate(segments, axis=0)
        return as_float_array(np.ascontiguousarray(trace), shape=(None, 3))
