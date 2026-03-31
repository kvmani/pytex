from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np

from pytex.core._arrays import FloatArray, IntArray, as_float_array, as_int_array
from pytex.core.lattice import CrystalPlane, Phase, ZoneAxis

_ROUND_DECIMALS = 12


def _direct_basis_matrix(phase: Phase) -> FloatArray:
    return as_float_array(phase.lattice.direct_basis().matrix, shape=(3, 3))


def _reciprocal_basis_matrix(phase: Phase) -> FloatArray:
    return as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))


def _require_matching_phases(
    left: Phase,
    right: Phase,
    *,
    left_name: str,
    right_name: str,
) -> None:
    if left != right:
        raise ValueError(f"{left_name}.phase must match {right_name}.phase.")


def _as_index_rows(values: Any, *, name: str, columns: int) -> IntArray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (columns,):
            raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
        array = array[None, :]
    elif array.ndim == 2:
        if array.shape[1] != columns:
            raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
    else:
        raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
    return as_int_array(array, shape=(None, columns))


def _validate_nonzero_rows(indices: IntArray, *, name: str) -> None:
    if np.any(~np.any(indices != 0, axis=1)):
        raise ValueError(f"{name} must not contain the zero triplet.")


def _rowwise_gcd(values: IntArray) -> IntArray:
    if values.shape[0] == 0:
        return as_int_array(np.ones(0, dtype=np.int64), shape=(0,))
    abs_values = np.abs(values)
    divisors = np.gcd.reduce(abs_values, axis=1)
    divisors[divisors == 0] = 1
    return as_int_array(divisors, shape=(values.shape[0],))


def _first_nonzero_signs(values: IntArray) -> IntArray:
    non_zero = values != 0
    first_index = np.argmax(non_zero, axis=1)
    rows = np.arange(values.shape[0], dtype=np.int64)
    selected = values[rows, first_index]
    signs = np.sign(selected).astype(np.int64, copy=False)
    signs[signs == 0] = 1
    return as_int_array(signs, shape=(values.shape[0],))


def reduce_indices(values: Any, *, name: str = "indices") -> IntArray:
    rows = _as_index_rows(values, name=name, columns=3)
    _validate_nonzero_rows(rows, name=name)
    divisors = _rowwise_gcd(rows)
    reduced = rows // divisors[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 3))


def canonicalize_sign(values: Any, *, name: str = "indices") -> IntArray:
    rows = _as_index_rows(values, name=name, columns=3)
    _validate_nonzero_rows(rows, name=name)
    signs = _first_nonzero_signs(rows)
    canonical = rows * signs[:, None]
    return as_int_array(canonical, shape=(rows.shape[0], 3))


def canonicalize_family_indices(
    values: Any,
    *,
    antipodal: bool,
    name: str = "indices",
) -> IntArray:
    reduced = reduce_indices(values, name=name)
    if antipodal:
        return canonicalize_sign(reduced, name=name)
    return reduced


def antipodal_keys(values: Any, *, name: str = "indices") -> IntArray:
    return canonicalize_family_indices(values, antipodal=True, name=name)


def plane_hkl_to_hkil_array(hkl: Any) -> IntArray:
    rows = _as_index_rows(hkl, name="hkl", columns=3)
    h = rows[:, 0]
    k = rows[:, 1]
    ell = rows[:, 2]
    converted = np.column_stack([h, k, -(h + k), ell])
    return as_int_array(converted, shape=(rows.shape[0], 4))


def plane_hkil_to_hkl_array(hkil: Any) -> IntArray:
    array = np.ascontiguousarray(np.asarray(hkil, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (4,):
            raise ValueError("hkil must have shape (4,) or (n, 4).")
        array = array[None, :]
    elif array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("hkil must have shape (4,) or (n, 4).")
    rows = as_int_array(array, shape=(None, 4))
    if np.any(rows[:, 2] != -(rows[:, 0] + rows[:, 1])):
        raise ValueError("Hexagonal four-index planes must satisfy i = -(h + k).")
    converted = np.column_stack([rows[:, 0], rows[:, 1], rows[:, 3]])
    return as_int_array(converted, shape=(rows.shape[0], 3))


def direction_uvw_to_uvtw_array(uvw: Any) -> IntArray:
    rows = _as_index_rows(uvw, name="uvw", columns=3)
    two_u_minus_v = 2 * rows[:, 0] - rows[:, 1]
    two_v_minus_u = 2 * rows[:, 1] - rows[:, 0]
    minus_u_minus_v = -(rows[:, 0] + rows[:, 1])
    divisible_by_three = (
        (two_u_minus_v % 3 == 0)
        & (two_v_minus_u % 3 == 0)
        & (minus_u_minus_v % 3 == 0)
    )
    reduced_form = np.column_stack(
        [two_u_minus_v // 3, two_v_minus_u // 3, minus_u_minus_v // 3, rows[:, 2]]
    )
    expanded_form = np.column_stack(
        [two_u_minus_v, two_v_minus_u, minus_u_minus_v, 3 * rows[:, 2]]
    )
    scaled = np.where(divisible_by_three[:, None], reduced_form, expanded_form)
    reduced = scaled // _rowwise_gcd(scaled)[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 4))


def direction_uvtw_to_uvw_array(uvtw: Any) -> IntArray:
    array = np.ascontiguousarray(np.asarray(uvtw, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (4,):
            raise ValueError("UVTW must have shape (4,) or (n, 4).")
        array = array[None, :]
    elif array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("UVTW must have shape (4,) or (n, 4).")
    rows = as_int_array(array, shape=(None, 4))
    if np.any(np.sum(rows[:, :3], axis=1) != 0):
        raise ValueError("Hexagonal four-index directions must satisfy U + V + T = 0.")
    converted = np.column_stack(
        [2 * rows[:, 0] + rows[:, 1], 2 * rows[:, 1] + rows[:, 0], rows[:, 3]]
    )
    reduced = converted // _rowwise_gcd(converted)[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 3))


def _family_unique_rows(values: IntArray) -> tuple[IntArray, IntArray]:
    if values.shape[0] == 0:
        empty_rows = np.empty((0, values.shape[1]), dtype=np.int64)
        empty_inverse = np.empty(0, dtype=np.int64)
        return (
            as_int_array(empty_rows, shape=(0, values.shape[1])),
            as_int_array(empty_inverse, shape=(0,)),
        )
    unique_rows, inverse = np.unique(values, axis=0, return_inverse=True)
    return (
        as_int_array(unique_rows, shape=(unique_rows.shape[0], unique_rows.shape[1])),
        as_int_array(inverse, shape=(inverse.shape[0],)),
    )


def _freeze_bool_array(values: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype=bool)
    array.setflags(write=False)
    return array


def _broadcast_rows(
    left: np.ndarray,
    right: np.ndarray,
    *,
    left_name: str,
    right_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    if left.shape[0] == right.shape[0]:
        return left, right
    if left.shape[0] == 1:
        return np.broadcast_to(left, right.shape), right
    if right.shape[0] == 1:
        return left, np.broadcast_to(right, left.shape)
    raise ValueError(
        f"{left_name} and {right_name} must have matching lengths or one row for broadcasting."
    )


def _cartesian_from_direct_indices(indices: IntArray, phase: Phase) -> FloatArray:
    basis = _direct_basis_matrix(phase)
    cartesian = np.asarray(indices, dtype=np.float64) @ basis.T
    return as_float_array(cartesian, shape=(indices.shape[0], 3))


def _cartesian_from_reciprocal_indices(indices: IntArray, phase: Phase) -> FloatArray:
    basis = _reciprocal_basis_matrix(phase)
    cartesian = np.asarray(indices, dtype=np.float64) @ basis.T
    return as_float_array(cartesian, shape=(indices.shape[0], 3))


def _unit_vectors(vectors: FloatArray, *, name: str) -> FloatArray:
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError(f"{name} contains a zero cartesian vector.")
    unit = vectors / norms[:, None]
    return as_float_array(unit, shape=(vectors.shape[0], 3))


def _structured_sort_rows(values: IntArray) -> IntArray:
    if values.size == 0:
        return values
    structured = values.view(dtype=[("c0", np.int64), ("c1", np.int64), ("c2", np.int64)]).reshape(
        values.shape[0], values.shape[1]
    )
    sorted_structured = np.sort(structured, axis=1)
    sorted_values = sorted_structured.view(np.int64).reshape(values.shape[0], values.shape[1], 3)
    return as_int_array(sorted_values, shape=(values.shape[0], values.shape[1], 3))


def _pack_unique_rows(values: IntArray) -> tuple[IntArray, np.ndarray]:
    if values.shape[0] == 0:
        empty_values = np.empty((0, 0, 3), dtype=np.int64)
        empty_mask = np.empty((0, 0), dtype=bool)
        return as_int_array(empty_values, shape=(0, 0, 3)), _freeze_bool_array(empty_mask)
    sorted_values = _structured_sort_rows(values)
    unique_mask = np.ones(sorted_values.shape[:2], dtype=bool)
    if sorted_values.shape[1] > 1:
        unique_mask[:, 1:] = np.any(
            sorted_values[:, 1:, :] != sorted_values[:, :-1, :],
            axis=2,
        )
    counts = np.sum(unique_mask, axis=1)
    max_count = int(np.max(counts)) if counts.size else 0
    packed = np.zeros((sorted_values.shape[0], max_count, 3), dtype=np.int64)
    packed_mask = np.zeros((sorted_values.shape[0], max_count), dtype=bool)
    if max_count:
        row_ids = np.broadcast_to(
            np.arange(sorted_values.shape[0], dtype=np.int64)[:, None],
            unique_mask.shape,
        )
        destination = np.cumsum(unique_mask, axis=1) - 1
        packed[row_ids[unique_mask], destination[unique_mask]] = sorted_values[unique_mask]
        packed_mask[row_ids[unique_mask], destination[unique_mask]] = True
    return (
        as_int_array(packed, shape=(packed.shape[0], packed.shape[1], 3)),
        _freeze_bool_array(packed_mask),
    )


def _recover_integer_indices(
    cartesian_vectors: np.ndarray,
    *,
    phase: Phase,
    reciprocal: bool,
    tol: float,
) -> IntArray:
    basis = _reciprocal_basis_matrix(phase) if reciprocal else _direct_basis_matrix(phase)
    inverse = np.linalg.inv(np.asarray(basis, dtype=np.float64))
    coordinates = np.einsum("...j,ij->...i", cartesian_vectors, inverse, optimize=True)
    rounded = np.rint(coordinates)
    if not np.allclose(coordinates, rounded, atol=tol, rtol=0.0):
        raise ValueError(
            "Symmetry operators produced coordinates that could not be recovered as "
            "integer Miller indices."
        )
    return as_int_array(rounded.astype(np.int64), shape=coordinates.shape)


def _symmetry_equivalent_indices(
    indices: IntArray,
    *,
    phase: Phase,
    reciprocal: bool,
    unique: bool,
    antipodal: bool,
    tol: float,
) -> tuple[IntArray, np.ndarray]:
    if indices.shape[0] == 0:
        empty_values = np.empty((0, 0, 3), dtype=np.int64)
        empty_mask = np.empty((0, 0), dtype=bool)
        return as_int_array(empty_values, shape=(0, 0, 3)), _freeze_bool_array(empty_mask)
    cartesian = (
        _cartesian_from_reciprocal_indices(indices, phase)
        if reciprocal
        else _cartesian_from_direct_indices(indices, phase)
    )
    transformed = np.einsum(
        "oij,nj->noi",
        phase.symmetry.operators,
        np.asarray(cartesian, dtype=np.float64),
        optimize=True,
    )
    recovered = _recover_integer_indices(
        transformed,
        phase=phase,
        reciprocal=reciprocal,
        tol=tol,
    )
    flattened = recovered.reshape(-1, 3)
    canonical = canonicalize_family_indices(
        flattened,
        antipodal=antipodal,
        name="symmetry_equivalent_indices",
    ).reshape(recovered.shape)
    if unique:
        return _pack_unique_rows(as_int_array(canonical, shape=recovered.shape))
    mask = np.ones(canonical.shape[:2], dtype=bool)
    return (
        as_int_array(canonical, shape=canonical.shape),
        _freeze_bool_array(mask),
    )


def _pairwise_angles_from_unit_vectors(
    left: FloatArray,
    right: FloatArray,
    *,
    left_name: str,
    right_name: str,
    antipodal: bool,
) -> FloatArray:
    broadcast_left, broadcast_right = _broadcast_rows(
        left,
        right,
        left_name=left_name,
        right_name=right_name,
    )
    dots = np.einsum("ni,ni->n", broadcast_left, broadcast_right, optimize=True)
    if antipodal:
        dots = np.abs(dots)
    return as_float_array(np.arccos(np.clip(dots, -1.0, 1.0)), shape=(dots.shape[0],))


def _pairwise_matrix_angles_from_unit_vectors(
    left: FloatArray,
    right: FloatArray,
    *,
    antipodal: bool,
) -> FloatArray:
    dots = np.asarray(left, dtype=np.float64) @ np.asarray(right, dtype=np.float64).T
    if antipodal:
        dots = np.abs(dots)
    return as_float_array(
        np.arccos(np.clip(dots, -1.0, 1.0)),
        shape=(left.shape[0], right.shape[0]),
    )


def _projection_vectors(
    directions: FloatArray,
    plane_normals: FloatArray,
) -> tuple[FloatArray, np.ndarray]:
    broadcast_directions, broadcast_normals = _broadcast_rows(
        directions,
        plane_normals,
        left_name="directions",
        right_name="planes",
    )
    normal_components = np.einsum(
        "ni,ni->n",
        broadcast_directions,
        broadcast_normals,
        optimize=True,
    )
    projected = broadcast_directions - normal_components[:, None] * broadcast_normals
    degenerate_mask = np.isclose(np.linalg.norm(projected, axis=1), 0.0)
    projected[degenerate_mask] = 0.0
    return (
        as_float_array(projected, shape=(projected.shape[0], 3)),
        _freeze_bool_array(degenerate_mask),
    )


@dataclass(frozen=True, slots=True)
class MillerPlane:
    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("MillerPlane indices must not be the zero triplet.")

    @classmethod
    def from_hkl(cls, indices: Any, *, phase: Phase) -> MillerPlane:
        return cls(indices=indices, phase=phase)

    @classmethod
    def from_hkil(cls, indices: Any, *, phase: Phase) -> MillerPlane:
        return cls(indices=plane_hkil_to_hkl_array(indices)[0], phase=phase)

    @classmethod
    def from_crystal_plane(cls, plane: CrystalPlane) -> MillerPlane:
        return cls(indices=plane.miller.indices, phase=plane.phase)

    @property
    def hkil(self) -> IntArray:
        return as_int_array(plane_hkl_to_hkil_array(self.indices)[0], shape=(4,))

    @property
    def reduced_indices(self) -> IntArray:
        return as_int_array(reduce_indices(self.indices)[0], shape=(3,))

    @property
    def antipodal_key(self) -> IntArray:
        return as_int_array(antipodal_keys(self.indices)[0], shape=(3,))

    @property
    def reciprocal_vector_cartesian(self) -> FloatArray:
        return as_float_array(
            _cartesian_from_reciprocal_indices(self.indices[None, :], self.phase)[0],
            shape=(3,),
        )

    @property
    def normal_cartesian(self) -> FloatArray:
        return as_float_array(
            _unit_vectors(self.reciprocal_vector_cartesian[None, :], name="MillerPlane")[0],
            shape=(3,),
        )

    @property
    def d_spacing_angstrom(self) -> float:
        magnitude = float(np.linalg.norm(self.reciprocal_vector_cartesian))
        if np.isclose(magnitude, 0.0):
            raise ValueError("MillerPlane reciprocal vector magnitude must be non-zero.")
        return 1.0 / magnitude

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        if not antipodal:
            raise ValueError("Miller planes are always treated with antipodal equivalence.")
        return _symmetry_equivalent_indices(
            self.indices[None, :],
            phase=self.phase,
            reciprocal=True,
            unique=unique,
            antipodal=True,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> MillerPlaneSet:
        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return MillerPlaneSet(indices=equivalent_indices[0, mask[0]], phase=self.phase)

    def to_crystal_plane(self) -> CrystalPlane:
        return CrystalPlane.from_miller_bravais(self.hkil, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerDirection:
    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("MillerDirection indices must not be the zero triplet.")

    @classmethod
    def from_uvw(cls, indices: Any, *, phase: Phase) -> MillerDirection:
        return cls(indices=indices, phase=phase)

    @classmethod
    def from_UVTW(cls, indices: Any, *, phase: Phase) -> MillerDirection:  # noqa: N802
        return cls(indices=direction_uvtw_to_uvw_array(indices)[0], phase=phase)

    @classmethod
    def from_zone_axis(cls, zone_axis: ZoneAxis) -> MillerDirection:
        return cls(indices=zone_axis.indices, phase=zone_axis.phase)

    @property
    def UVTW(self) -> IntArray:  # noqa: N802
        return as_int_array(direction_uvw_to_uvtw_array(self.indices)[0], shape=(4,))

    @property
    def reduced_indices(self) -> IntArray:
        return as_int_array(reduce_indices(self.indices)[0], shape=(3,))

    @property
    def antipodal_key(self) -> IntArray:
        return as_int_array(antipodal_keys(self.indices)[0], shape=(3,))

    @property
    def direct_vector_cartesian(self) -> FloatArray:
        return as_float_array(
            _cartesian_from_direct_indices(self.indices[None, :], self.phase)[0],
            shape=(3,),
        )

    @property
    def unit_vector_cartesian(self) -> FloatArray:
        return as_float_array(
            _unit_vectors(self.direct_vector_cartesian[None, :], name="MillerDirection")[0],
            shape=(3,),
        )

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        return _symmetry_equivalent_indices(
            self.indices[None, :],
            phase=self.phase,
            reciprocal=False,
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> MillerDirectionSet:
        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return MillerDirectionSet(indices=equivalent_indices[0, mask[0]], phase=self.phase)

    def to_zone_axis(self) -> ZoneAxis:
        return ZoneAxis(indices=self.indices, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerPlaneSet:
    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(None, 3)))
        _validate_nonzero_rows(cast(IntArray, self.indices), name="MillerPlaneSet.indices")

    @classmethod
    def from_hkl(cls, indices: Any, *, phase: Phase) -> MillerPlaneSet:
        return cls(indices=_as_index_rows(indices, name="hkl", columns=3), phase=phase)

    @classmethod
    def from_hkil(cls, indices: Any, *, phase: Phase) -> MillerPlaneSet:
        return cls(indices=plane_hkil_to_hkl_array(indices), phase=phase)

    def to_hkil(self) -> IntArray:
        return plane_hkl_to_hkil_array(self.indices)

    def reduce_indices(self) -> IntArray:
        return reduce_indices(self.indices)

    def canonical_indices(self) -> IntArray:
        return antipodal_keys(self.indices)

    def unique(self) -> tuple[MillerPlaneSet, IntArray]:
        unique_indices, inverse = _family_unique_rows(self.canonical_indices())
        return MillerPlaneSet(indices=unique_indices, phase=self.phase), inverse

    def reciprocal_vectors_cartesian(self) -> FloatArray:
        return _cartesian_from_reciprocal_indices(cast(IntArray, self.indices), self.phase)

    def normals_cartesian(self) -> FloatArray:
        return _unit_vectors(self.reciprocal_vectors_cartesian(), name="MillerPlaneSet")

    def d_spacings_angstrom(self) -> FloatArray:
        magnitudes = np.linalg.norm(self.reciprocal_vectors_cartesian(), axis=1)
        if np.any(np.isclose(magnitudes, 0.0)):
            raise ValueError("MillerPlaneSet reciprocal vectors must be non-zero.")
        return as_float_array(1.0 / magnitudes, shape=(self.indices.shape[0],))

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        if not antipodal:
            raise ValueError("Miller planes are always treated with antipodal equivalence.")
        return _symmetry_equivalent_indices(
            cast(IntArray, self.indices),
            phase=self.phase,
            reciprocal=True,
            unique=unique,
            antipodal=True,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[MillerPlaneSet, ...]:
        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return tuple(
            MillerPlaneSet(indices=equivalent_indices[row, mask[row]], phase=self.phase)
            for row in range(equivalent_indices.shape[0])
        )

    def to_crystal_planes(self) -> tuple[CrystalPlane, ...]:
        return tuple(
            MillerPlane(indices=row, phase=self.phase).to_crystal_plane()
            for row in self.indices
        )

    def angle_matrix_rad(self, other: MillerPlaneSet | None = None) -> FloatArray:
        target = self if other is None else other
        _require_matching_phases(
            self.phase,
            target.phase,
            left_name="MillerPlaneSet",
            right_name="other",
        )
        return _pairwise_matrix_angles_from_unit_vectors(
            self.normals_cartesian(),
            target.normals_cartesian(),
            antipodal=True,
        )


@dataclass(frozen=True, slots=True)
class MillerDirectionSet:
    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(None, 3)))
        _validate_nonzero_rows(cast(IntArray, self.indices), name="MillerDirectionSet.indices")

    @classmethod
    def from_uvw(cls, indices: Any, *, phase: Phase) -> MillerDirectionSet:
        return cls(indices=_as_index_rows(indices, name="uvw", columns=3), phase=phase)

    @classmethod
    def from_UVTW(cls, indices: Any, *, phase: Phase) -> MillerDirectionSet:  # noqa: N802
        return cls(indices=direction_uvtw_to_uvw_array(indices), phase=phase)

    def to_UVTW(self) -> IntArray:  # noqa: N802
        return direction_uvw_to_uvtw_array(self.indices)

    def reduce_indices(self) -> IntArray:
        return reduce_indices(self.indices)

    def canonical_indices(self, *, antipodal: bool = True) -> IntArray:
        return canonicalize_family_indices(self.indices, antipodal=antipodal)

    def unique(self, *, antipodal: bool = True) -> tuple[MillerDirectionSet, IntArray]:
        unique_indices, inverse = _family_unique_rows(self.canonical_indices(antipodal=antipodal))
        return MillerDirectionSet(indices=unique_indices, phase=self.phase), inverse

    def direct_vectors_cartesian(self) -> FloatArray:
        return _cartesian_from_direct_indices(cast(IntArray, self.indices), self.phase)

    def unit_vectors_cartesian(self) -> FloatArray:
        return _unit_vectors(self.direct_vectors_cartesian(), name="MillerDirectionSet")

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        return _symmetry_equivalent_indices(
            cast(IntArray, self.indices),
            phase=self.phase,
            reciprocal=False,
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[MillerDirectionSet, ...]:
        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return tuple(
            MillerDirectionSet(indices=equivalent_indices[row, mask[row]], phase=self.phase)
            for row in range(equivalent_indices.shape[0])
        )

    def to_zone_axes(self) -> tuple[ZoneAxis, ...]:
        return tuple(
            MillerDirection(indices=row, phase=self.phase).to_zone_axis()
            for row in self.indices
        )

    def angle_matrix_rad(
        self,
        other: MillerDirectionSet | None = None,
        *,
        antipodal: bool = True,
    ) -> FloatArray:
        target = self if other is None else other
        _require_matching_phases(
            self.phase,
            target.phase,
            left_name="MillerDirectionSet",
            right_name="other",
        )
        return _pairwise_matrix_angles_from_unit_vectors(
            self.unit_vectors_cartesian(),
            target.unit_vectors_cartesian(),
            antipodal=antipodal,
        )


def _plane_normals(planes: MillerPlane | MillerPlaneSet) -> tuple[FloatArray, Phase]:
    if isinstance(planes, MillerPlane):
        return planes.normal_cartesian[None, :], planes.phase
    return planes.normals_cartesian(), planes.phase


def _direction_units(directions: MillerDirection | MillerDirectionSet) -> tuple[FloatArray, Phase]:
    if isinstance(directions, MillerDirection):
        return directions.unit_vector_cartesian[None, :], directions.phase
    return directions.unit_vectors_cartesian(), directions.phase


def angle_plane_plane_rad(
    left: MillerPlane | MillerPlaneSet,
    right: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    left_normals, left_phase = _plane_normals(left)
    right_normals, right_phase = _plane_normals(right)
    _require_matching_phases(left_phase, right_phase, left_name="left", right_name="right")
    angles = _pairwise_angles_from_unit_vectors(
        left_normals,
        right_normals,
        left_name="left",
        right_name="right",
        antipodal=True,
    )
    if isinstance(left, MillerPlane) and isinstance(right, MillerPlane):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_dir_rad(
    left: MillerDirection | MillerDirectionSet,
    right: MillerDirection | MillerDirectionSet,
    *,
    antipodal: bool = True,
) -> float | FloatArray:
    left_units, left_phase = _direction_units(left)
    right_units, right_phase = _direction_units(right)
    _require_matching_phases(left_phase, right_phase, left_name="left", right_name="right")
    angles = _pairwise_angles_from_unit_vectors(
        left_units,
        right_units,
        left_name="left",
        right_name="right",
        antipodal=antipodal,
    )
    if isinstance(left, MillerDirection) and isinstance(right, MillerDirection):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_plane_normal_rad(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    direction_units, direction_phase = _direction_units(directions)
    plane_normals, plane_phase = _plane_normals(planes)
    _require_matching_phases(
        direction_phase,
        plane_phase,
        left_name="directions",
        right_name="planes",
    )
    angles = _pairwise_angles_from_unit_vectors(
        direction_units,
        plane_normals,
        left_name="directions",
        right_name="planes",
        antipodal=True,
    )
    if isinstance(directions, MillerDirection) and isinstance(planes, MillerPlane):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_plane_inclination_rad(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    normal_angles = angle_dir_plane_normal_rad(directions, planes)
    if isinstance(normal_angles, float):
        return float((np.pi / 2.0) - normal_angles)
    return as_float_array((np.pi / 2.0) - normal_angles, shape=(normal_angles.shape[0],))


def project_directions_onto_planes(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> tuple[FloatArray, np.ndarray]:
    direction_units, direction_phase = _direction_units(directions)
    plane_normals, plane_phase = _plane_normals(planes)
    _require_matching_phases(
        direction_phase,
        plane_phase,
        left_name="directions",
        right_name="planes",
    )
    return _projection_vectors(direction_units, plane_normals)


__all__ = [
    "MillerDirection",
    "MillerDirectionSet",
    "MillerPlane",
    "MillerPlaneSet",
    "angle_dir_dir_rad",
    "angle_dir_plane_inclination_rad",
    "angle_dir_plane_normal_rad",
    "angle_plane_plane_rad",
    "antipodal_keys",
    "canonicalize_family_indices",
    "canonicalize_sign",
    "direction_uvtw_to_uvw_array",
    "direction_uvw_to_uvtw_array",
    "plane_hkil_to_hkl_array",
    "plane_hkl_to_hkil_array",
    "project_directions_onto_planes",
    "reduce_indices",
]
