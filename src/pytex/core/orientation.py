from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import overload

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    normalize_quaternion,
    normalize_quaternions,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.batches import EulerSet, RotationSet, VectorSet, normalize_euler_convention_name
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalDirection, CrystalPlane, Phase
from pytex.core.miller import (
    MillerDirection,
    MillerDirectionSet,
    MillerPlane,
    MillerPlaneSet,
)
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

_SPECIMEN_DIRECTION_ALIASES = {
    "rd": np.array([1.0, 0.0, 0.0], dtype=np.float64),
    "x": np.array([1.0, 0.0, 0.0], dtype=np.float64),
    "td": np.array([0.0, 1.0, 0.0], dtype=np.float64),
    "y": np.array([0.0, 1.0, 0.0], dtype=np.float64),
    "nd": np.array([0.0, 0.0, 1.0], dtype=np.float64),
    "z": np.array([0.0, 0.0, 1.0], dtype=np.float64),
}


def _normalize_euler_convention(convention: str) -> str:
    return normalize_euler_convention_name(convention)


def specimen_direction_vector(direction: str | ArrayLike) -> np.ndarray:
    """Return a normalized specimen direction from a named alias or vector."""

    if isinstance(direction, str):
        normalized = direction.strip().lower()
        if normalized not in _SPECIMEN_DIRECTION_ALIASES:
            raise ValueError("Specimen direction must be one of RD, TD, ND, x, y, z, or a vector.")
        return as_float_array(_SPECIMEN_DIRECTION_ALIASES[normalized], shape=(3,))
    return normalize_vector(direction)


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


def _matrices_to_repeated_axis_euler(
    matrices: np.ndarray,
    *,
    convention: str,
) -> np.ndarray:
    """Vectorised counterpart of `_matrix_to_repeated_axis_euler` (shape (n, 3)).

    Reproduces the scalar function's gimbal-lock branching exactly via masks so
    per-orientation results are identical.
    """

    normalized = _normalize_euler_convention(convention)
    phi = _safe_arccos(matrices[:, 2, 2])
    near_zero = np.isclose(phi, 0.0, atol=1e-10)
    near_pi = np.isclose(phi, np.pi, atol=1e-10)
    first = np.empty(matrices.shape[0], dtype=np.float64)
    third = np.zeros(matrices.shape[0], dtype=np.float64)
    if normalized == "bunge":
        regular = ~near_zero & ~near_pi
        first[near_zero] = np.arctan2(
            matrices[near_zero, 1, 0], matrices[near_zero, 0, 0]
        )
        first[near_pi] = np.arctan2(matrices[near_pi, 0, 1], matrices[near_pi, 0, 0])
        first[regular] = np.arctan2(matrices[regular, 0, 2], -matrices[regular, 1, 2])
        third[regular] = np.arctan2(matrices[regular, 2, 0], matrices[regular, 2, 1])
    else:
        degenerate = near_zero | near_pi
        regular = ~degenerate
        first[degenerate] = np.arctan2(
            matrices[degenerate, 1, 0], matrices[degenerate, 0, 0]
        )
        first[regular] = np.arctan2(matrices[regular, 1, 2], matrices[regular, 0, 2])
        third[regular] = np.arctan2(matrices[regular, 2, 1], -matrices[regular, 2, 0])
    return np.column_stack([first, phi, third])


def _batched_fundamental_representatives(
    quaternions: np.ndarray,
    *,
    crystal_operators: np.ndarray,
    specimen_operators: np.ndarray,
) -> np.ndarray:
    """Fundamental-region representative quaternion for each input orientation.

    For every orientation this builds the symmetry orbit ``specimen_op @ R @
    crystal_op`` and returns the orbit quaternion that is the lexicographically
    largest canonical quaternion -- exactly the representative selected by
    ``_canonical_quaternion_index`` / the minimum ``_fundamental_region_key`` in
    the per-orientation scalar path, but vectorised over the whole set.
    """

    quaternion_array = np.asarray(quaternions, dtype=np.float64)
    count = quaternion_array.shape[0]
    if count == 0:
        return np.empty((0, 4), dtype=np.float64)
    matrices = quaternions_to_matrices(quaternion_array)
    left_applied = np.einsum("aij,njk->naik", specimen_operators, matrices, optimize=True)
    orbit = np.einsum("naik,bkl->nabil", left_applied, crystal_operators, optimize=True)
    orbit = orbit.reshape(count, -1, 3, 3)
    orbit_count = orbit.shape[1]
    orbit_quaternions = matrices_to_quaternions(orbit.reshape(-1, 3, 3)).reshape(
        count, orbit_count, 4
    )
    # Canonicalise (w >= 0, unit norm, rounded) only to choose the representative.
    canonical = orbit_quaternions.copy()
    negative = canonical[..., 0] < 0.0
    canonical[negative] = -canonical[negative]
    canonical = canonical / np.linalg.norm(canonical, axis=-1, keepdims=True)
    rounded = np.round(canonical, 12)
    # Per-row lexicographic argmax over the orbit by (w, x, y, z), narrowing ties.
    valid = np.ones((count, orbit_count), dtype=bool)
    for component in range(4):
        masked = np.where(valid, rounded[:, :, component], -np.inf)
        column_max = masked.max(axis=1, keepdims=True)
        valid &= rounded[:, :, component] >= column_max
    selected = np.argmax(valid, axis=1)
    return np.ascontiguousarray(orbit_quaternions[np.arange(count), selected])


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


def quaternions_multiply(left: ArrayLike, right: ArrayLike) -> np.ndarray:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    lw, lx, ly, lz = np.moveaxis(left_array, -1, 0)
    rw, rx, ry, rz = np.moveaxis(right_array, -1, 0)
    return np.stack(
        [
            lw * rw - lx * rx - ly * ry - lz * rz,
            lw * rx + lx * rw + ly * rz - lz * ry,
            lw * ry - lx * rz + ly * rw + lz * rx,
            lw * rz + lx * ry - ly * rx + lz * rw,
        ],
        axis=-1,
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


def quaternions_to_matrices(quaternions: ArrayLike) -> np.ndarray:
    quaternion_array = normalize_quaternions(quaternions)
    w = quaternion_array[:, 0]
    x = quaternion_array[:, 1]
    y = quaternion_array[:, 2]
    z = quaternion_array[:, 3]
    matrices = np.stack(
        [
            np.stack(
                [
                    1.0 - 2.0 * (y * y + z * z),
                    2.0 * (x * y - z * w),
                    2.0 * (x * z + y * w),
                ],
                axis=1,
            ),
            np.stack(
                [
                    2.0 * (x * y + z * w),
                    1.0 - 2.0 * (x * x + z * z),
                    2.0 * (y * z - x * w),
                ],
                axis=1,
            ),
            np.stack(
                [
                    2.0 * (x * z - y * w),
                    2.0 * (y * z + x * w),
                    1.0 - 2.0 * (x * x + y * y),
                ],
                axis=1,
            ),
        ],
        axis=1,
    )
    matrices = np.ascontiguousarray(matrices, dtype=np.float64)
    matrices.setflags(write=False)
    return matrices


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


def matrices_to_quaternions(matrices: ArrayLike) -> np.ndarray:
    matrix_array = np.asarray(matrices, dtype=np.float64)
    if matrix_array.shape == (3, 3):
        matrix_array = matrix_array[None, :, :]
    if matrix_array.ndim != 3 or matrix_array.shape[1:] != (3, 3):
        raise ValueError("Rotation matrices must have shape (3, 3) or (n, 3, 3).")
    # Vectorised proper-rotation validation: orthonormal (M^T M = I) and det = 1.
    gram = np.einsum("nji,njk->nik", matrix_array, matrix_array, optimize=True)
    determinants = np.linalg.det(matrix_array)
    if not (
        np.allclose(gram, np.eye(3, dtype=np.float64)[None, :, :], atol=1e-8)
        and np.allclose(determinants, 1.0, atol=1e-8)
    ):
        raise ValueError("All matrices must be proper rotation matrices.")
    quaternions = np.empty((matrix_array.shape[0], 4), dtype=np.float64)
    trace = np.trace(matrix_array, axis1=1, axis2=2)
    positive_trace = trace > 0.0
    if np.any(positive_trace):
        s = np.sqrt(trace[positive_trace] + 1.0) * 2.0
        quaternions[positive_trace, 0] = 0.25 * s
        quaternions[positive_trace, 1] = (
            matrix_array[positive_trace, 2, 1] - matrix_array[positive_trace, 1, 2]
        ) / s
        quaternions[positive_trace, 2] = (
            matrix_array[positive_trace, 0, 2] - matrix_array[positive_trace, 2, 0]
        ) / s
        quaternions[positive_trace, 3] = (
            matrix_array[positive_trace, 1, 0] - matrix_array[positive_trace, 0, 1]
        ) / s
    non_positive = ~positive_trace
    if np.any(non_positive):
        diagonal = np.stack(
            [
                matrix_array[non_positive, 0, 0],
                matrix_array[non_positive, 1, 1],
                matrix_array[non_positive, 2, 2],
            ],
            axis=1,
        )
        dominant = np.argmax(diagonal, axis=1)
        dominant_x = non_positive.copy()
        dominant_x[non_positive] = dominant == 0
        dominant_y = non_positive.copy()
        dominant_y[non_positive] = dominant == 1
        dominant_z = non_positive.copy()
        dominant_z[non_positive] = dominant == 2
        if np.any(dominant_x):
            s = np.sqrt(
                1.0
                + matrix_array[dominant_x, 0, 0]
                - matrix_array[dominant_x, 1, 1]
                - matrix_array[dominant_x, 2, 2]
            ) * 2.0
            quaternions[dominant_x, 0] = (
                matrix_array[dominant_x, 2, 1] - matrix_array[dominant_x, 1, 2]
            ) / s
            quaternions[dominant_x, 1] = 0.25 * s
            quaternions[dominant_x, 2] = (
                matrix_array[dominant_x, 0, 1] + matrix_array[dominant_x, 1, 0]
            ) / s
            quaternions[dominant_x, 3] = (
                matrix_array[dominant_x, 0, 2] + matrix_array[dominant_x, 2, 0]
            ) / s
        if np.any(dominant_y):
            s = np.sqrt(
                1.0
                + matrix_array[dominant_y, 1, 1]
                - matrix_array[dominant_y, 0, 0]
                - matrix_array[dominant_y, 2, 2]
            ) * 2.0
            quaternions[dominant_y, 0] = (
                matrix_array[dominant_y, 0, 2] - matrix_array[dominant_y, 2, 0]
            ) / s
            quaternions[dominant_y, 1] = (
                matrix_array[dominant_y, 0, 1] + matrix_array[dominant_y, 1, 0]
            ) / s
            quaternions[dominant_y, 2] = 0.25 * s
            quaternions[dominant_y, 3] = (
                matrix_array[dominant_y, 1, 2] + matrix_array[dominant_y, 2, 1]
            ) / s
        if np.any(dominant_z):
            s = np.sqrt(
                1.0
                + matrix_array[dominant_z, 2, 2]
                - matrix_array[dominant_z, 0, 0]
                - matrix_array[dominant_z, 1, 1]
            ) * 2.0
            quaternions[dominant_z, 0] = (
                matrix_array[dominant_z, 1, 0] - matrix_array[dominant_z, 0, 1]
            ) / s
            quaternions[dominant_z, 1] = (
                matrix_array[dominant_z, 0, 2] + matrix_array[dominant_z, 2, 0]
            ) / s
            quaternions[dominant_z, 2] = (
                matrix_array[dominant_z, 1, 2] + matrix_array[dominant_z, 2, 1]
            ) / s
            quaternions[dominant_z, 3] = 0.25 * s
    quaternions = normalize_quaternions(quaternions)
    quaternions.setflags(write=False)
    return quaternions


def quaternion_from_axis_angle(axis: ArrayLike, angle_rad: float) -> np.ndarray:
    unit_axis = normalize_vector(axis)
    half = angle_rad / 2.0
    sin_half = np.sin(half)
    return normalize_quaternion([np.cos(half), *(unit_axis * sin_half)])


def _broadcast_rotation_inputs(
    vectors: ArrayLike,
    scalars: ArrayLike,
    *,
    vector_name: str,
    scalar_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    vector_array = np.asarray(vectors, dtype=np.float64)
    if vector_array.shape == (3,):
        vector_array = vector_array[None, :]
    if vector_array.ndim != 2 or vector_array.shape[1] != 3:
        raise ValueError(f"{vector_name} must have shape (3,) or (n, 3).")
    scalar_array = np.asarray(scalars, dtype=np.float64)
    if scalar_array.ndim == 0:
        scalar_array = scalar_array.reshape(1)
    if scalar_array.ndim != 1:
        raise ValueError(f"{scalar_name} must be a scalar or a 1D array.")
    vector_count = int(vector_array.shape[0])
    scalar_count = int(scalar_array.shape[0])
    if vector_count == 1 and scalar_count > 1:
        vector_array = np.broadcast_to(vector_array, (scalar_count, 3))
    elif scalar_count == 1 and vector_count > 1:
        scalar_array = np.broadcast_to(scalar_array, (vector_count,))
    elif vector_count != scalar_count:
        raise ValueError(
            f"{vector_name} and {scalar_name} must broadcast to the same leading length."
        )
    return np.ascontiguousarray(vector_array), np.ascontiguousarray(scalar_array)


def quaternions_to_axes_angles(quaternions: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    quaternion_array = normalize_quaternions(quaternions)
    canonical = np.array(quaternion_array, copy=True)
    canonical[canonical[:, 0] < 0.0] *= -1.0
    angles = 2.0 * np.arccos(np.clip(canonical[:, 0], -1.0, 1.0))
    sin_half = np.linalg.norm(canonical[:, 1:], axis=1)
    axes = np.zeros((canonical.shape[0], 3), dtype=np.float64)
    axes[:, 2] = 1.0
    non_identity = sin_half > 1e-12
    if np.any(non_identity):
        axes[non_identity] = canonical[non_identity, 1:] / sin_half[non_identity, None]
        axes[non_identity] = normalize_vectors(axes[non_identity])
    axes = np.ascontiguousarray(axes)
    angles = np.ascontiguousarray(angles)
    axes.setflags(write=False)
    angles.setflags(write=False)
    return axes, angles


def quaternions_from_axes_angles(axes: ArrayLike, angles_rad: ArrayLike) -> np.ndarray:
    axis_array, angle_array = _broadcast_rotation_inputs(
        axes,
        angles_rad,
        vector_name="axes",
        scalar_name="angles_rad",
    )
    if np.any(~np.isfinite(angle_array)):
        raise ValueError("angles_rad must contain only finite values.")
    unit_axes = normalize_vectors(axis_array)
    half_angles = angle_array / 2.0
    sin_half = np.sin(half_angles)
    quaternions = np.column_stack(
        [
            np.cos(half_angles),
            unit_axes[:, 0] * sin_half,
            unit_axes[:, 1] * sin_half,
            unit_axes[:, 2] * sin_half,
        ]
    )
    quaternions = normalize_quaternions(quaternions)
    quaternions.setflags(write=False)
    return quaternions


def quaternions_to_rodrigues(quaternions: ArrayLike, *, frank: bool = False) -> np.ndarray:
    axes, angles = quaternions_to_axes_angles(quaternions)
    tan_half = np.tan(angles / 2.0)
    if not frank:
        rodrigues = axes * tan_half[:, None]
    else:
        frank_scale = np.where(np.isclose(angles, np.pi, atol=1e-12), np.inf, tan_half)
        rodrigues = np.column_stack([axes, frank_scale])
    rodrigues = np.ascontiguousarray(rodrigues, dtype=np.float64)
    rodrigues.setflags(write=False)
    return rodrigues


def quaternions_from_rodrigues(rodrigues: ArrayLike, *, frank: bool = False) -> np.ndarray:
    rodrigues_array = np.asarray(rodrigues, dtype=np.float64)
    if rodrigues_array.shape == (4 if frank else 3,):
        rodrigues_array = rodrigues_array[None, :]
    expected_dim = 4 if frank else 3
    if rodrigues_array.ndim != 2 or rodrigues_array.shape[1] != expected_dim:
        raise ValueError(
            "Rodrigues input must have shape "
            + f"({expected_dim},) or (n, {expected_dim}) for the selected convention."
        )
    if frank:
        axes = normalize_vectors(rodrigues_array[:, :3])
        scale = rodrigues_array[:, 3]
        angles = 2.0 * np.arctan(scale)
        infinite = np.isinf(scale)
        if np.any(infinite):
            angles = np.array(angles, copy=True)
            angles[infinite] = np.sign(scale[infinite]) * np.pi
        return quaternions_from_axes_angles(axes, angles)
    scale = np.linalg.norm(rodrigues_array, axis=1)
    axes = np.zeros_like(rodrigues_array)
    axes[:, 2] = 1.0
    non_identity = scale > 1e-12
    if np.any(non_identity):
        axes[non_identity] = rodrigues_array[non_identity] / scale[non_identity, None]
        axes[non_identity] = normalize_vectors(axes[non_identity])
    angles = 2.0 * np.arctan(scale)
    return quaternions_from_axes_angles(axes, angles)


def _grid_axis_values(
    start_deg: float,
    stop_deg: float,
    step_deg: float,
    *,
    periodic: bool,
) -> np.ndarray:
    if step_deg <= 0.0:
        raise ValueError("Euler-grid step sizes must be strictly positive.")
    if stop_deg <= start_deg:
        raise ValueError("Euler-grid stop values must be larger than the starts.")
    if periodic:
        values = np.arange(start_deg, stop_deg, step_deg, dtype=np.float64)
    else:
        values = np.arange(start_deg, stop_deg + 0.5 * step_deg, step_deg, dtype=np.float64)
    if values.size == 0:
        raise ValueError("Euler-grid generation produced no support orientations.")
    values = np.ascontiguousarray(values, dtype=np.float64)
    values.setflags(write=False)
    return values


def _require_grid_spacing(spacing_deg: float) -> float:
    spacing = float(spacing_deg)
    if not np.isfinite(spacing) or spacing <= 0.0 or spacing > 360.0:
        raise ValueError("spacing_deg must be finite and satisfy 0 < spacing_deg <= 360.")
    return spacing


def _orientation_grid_provenance(
    *,
    method: str,
    spacing_deg: float,
    note: str,
    provenance: ProvenanceRecord | None,
) -> ProvenanceRecord | None:
    if provenance is not None:
        return provenance
    return ProvenanceRecord(
        source_system="pytex.orientation_grid",
        metadata={
            "method": method,
            "spacing_deg": f"{float(spacing_deg):.12g}",
        },
        notes=(note,),
    )


def _canonicalize_quaternion_rows(quaternions: np.ndarray) -> np.ndarray:
    canonical = np.asarray(quaternions, dtype=np.float64).copy()
    canonical[canonical[:, 0] < 0.0] *= -1.0
    canonical = normalize_quaternions(canonical)
    canonical = np.ascontiguousarray(canonical, dtype=np.float64)
    canonical.setflags(write=False)
    return canonical


def _quaternion_left_matrices(quaternions: np.ndarray) -> np.ndarray:
    """Matrices ``L`` with ``L @ q = p * q`` for each ``p`` in ``quaternions``."""

    w, x, y, z = np.moveaxis(np.asarray(quaternions, dtype=np.float64), -1, 0)
    return np.stack(
        [
            np.stack([w, -x, -y, -z], axis=-1),
            np.stack([x, w, -z, y], axis=-1),
            np.stack([y, z, w, -x], axis=-1),
            np.stack([z, -y, x, w], axis=-1),
        ],
        axis=-2,
    )


def _quaternion_right_matrices(quaternions: np.ndarray) -> np.ndarray:
    """Matrices ``R`` with ``R @ q = q * p`` for each ``p`` in ``quaternions``."""

    w, x, y, z = np.moveaxis(np.asarray(quaternions, dtype=np.float64), -1, 0)
    return np.stack(
        [
            np.stack([w, -x, -y, -z], axis=-1),
            np.stack([x, w, z, -y], axis=-1),
            np.stack([y, -z, w, x], axis=-1),
            np.stack([z, y, -x, w], axis=-1),
        ],
        axis=-2,
    )


@lru_cache(maxsize=32)
def _disorientation_scalar_projection_cached(
    left_key: bytes,
    right_key: bytes,
    left_count: int,
    right_count: int,
) -> np.ndarray:
    """Rows ``a`` such that ``a . q`` is the scalar part of ``p * q * r^-1``.

    The disorientation angle only needs the scalar part of the symmetry-
    conjugated relative quaternion, and that scalar part is *linear* in the
    relative quaternion. Precomputing one row per operator pair turns the
    whole symmetry reduction into a single matrix product.
    """

    left_quaternions = matrices_to_quaternions(
        np.frombuffer(left_key, dtype=np.float64).reshape(left_count, 3, 3)
    )
    right_quaternions = matrices_to_quaternions(
        np.frombuffer(right_key, dtype=np.float64).reshape(right_count, 3, 3)
    )
    # q -> p * q; then q -> q * conj(r). Composing gives R(conj(r)) @ L(p),
    # whose first row is the linear functional for the scalar part.
    left_action = _quaternion_left_matrices(left_quaternions)
    conjugates = right_quaternions * np.array([1.0, -1.0, -1.0, -1.0])
    right_action = _quaternion_right_matrices(conjugates)
    combined = np.einsum("bij,ajk->abik", right_action, left_action, optimize=True)
    projection = combined[:, :, 0, :].reshape(-1, 4)

    # Only |a . q| matters, so rows that agree up to an overall sign are
    # redundant. For same-phase cubic symmetry this collapses 24 x 24 operator
    # pairs to 24 distinct functionals, cutting the reduction's memory traffic
    # by the same factor without changing any result.
    leading = np.argmax(np.abs(projection) > 1e-12, axis=1)
    signs = np.sign(projection[np.arange(projection.shape[0]), leading])
    signs[signs == 0.0] = 1.0
    canonical = np.round(projection * signs[:, None], 12)
    _, unique_rows = np.unique(canonical, axis=0, return_index=True)
    projection = projection[np.sort(unique_rows)]
    return np.ascontiguousarray(projection)


def _disorientation_scalar_projection(
    left_operators: np.ndarray,
    right_operators: np.ndarray,
) -> np.ndarray:
    left = np.ascontiguousarray(np.asarray(left_operators, dtype=np.float64))
    right = np.ascontiguousarray(np.asarray(right_operators, dtype=np.float64))
    return _disorientation_scalar_projection_cached(
        left.tobytes(), right.tobytes(), left.shape[0], right.shape[0]
    )


def _reduced_pair_disorientation_angles_from_quaternions(
    relative_quaternions: np.ndarray,
    left_operators: np.ndarray,
    right_operators: np.ndarray,
    *,
    max_block_elements: int = 8_000_000,
) -> np.ndarray:
    """Symmetry-reduced disorientation angle of each relative quaternion.

    Equivalent to ``min over S_l, S_r of angle(S_l @ M @ S_r^T)``, but expressed
    as one dense matrix product per block instead of a chain of einsums over
    the operator groups: ``angle = 2 * arccos(max_k |a_k . q|)``.
    """

    quaternions = np.asarray(relative_quaternions, dtype=np.float64)
    total = quaternions.shape[0]
    if total == 0:
        return np.empty(0, dtype=np.float64)
    projection = _disorientation_scalar_projection(left_operators, right_operators)
    per_row = max(projection.shape[0], 1)
    block = max(1, int(max_block_elements // per_row))
    angles = np.empty(total, dtype=np.float64)
    for start in range(0, total, block):
        stop = min(start + block, total)
        scalars = quaternions[start:stop] @ projection.T
        np.abs(scalars, out=scalars)
        best = np.clip(scalars.max(axis=1), 0.0, 1.0)
        angles[start:stop] = 2.0 * np.arccos(best)
    return angles


def _reduced_pair_disorientation_angles(
    relative_matrices: np.ndarray,
    left_operators: np.ndarray,
    right_operators: np.ndarray,
    *,
    max_block_elements: int = 8_000_000,
) -> np.ndarray:
    """Symmetry-reduced disorientation angle of each relative rotation.

    Returns ``min over S_l, S_r of angle(S_l @ M @ S_r^T)`` per matrix ``M``.
    Callers that already hold quaternions should use
    :func:`_reduced_pair_disorientation_angles_from_quaternions` and skip the
    matrix round trip entirely.
    """

    matrices = np.asarray(relative_matrices, dtype=np.float64)
    if matrices.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    return _reduced_pair_disorientation_angles_from_quaternions(
        matrices_to_quaternions(matrices),
        left_operators,
        right_operators,
        max_block_elements=max_block_elements,
    )


def _disorientation_medoid_index(
    orientations: OrientationSet,
    *,
    symmetry_aware: bool = True,
    max_pairs_per_block: int = 4_000_000,
    tie_rtol: float = 1e-9,
) -> int:
    """Index of the member with the least total disorientation to the others.

    Equivalent to summing the full pairwise matrix, but the row sums are
    accumulated in blocks so an n-member set never allocates an ``(n, n)``
    array.

    Grains routinely contain members whose total disorientation agrees to the
    last few bits — a symmetric cluster has no unique medoid. Choosing by bare
    ``argmin`` would then let the summation order, the BLAS build or the
    machine decide the grain reference orientation. Members within
    ``tie_rtol`` of the minimum are therefore treated as tied and the lowest
    index wins, which is reproducible everywhere. The tolerance sits many
    orders of magnitude below any physically meaningful separation.
    """

    count = len(orientations)
    if count == 0:
        raise ValueError("A medoid requires at least one orientation.")
    if count <= 2:
        return 0
    quaternions = np.asarray(orientations.quaternions, dtype=np.float64)
    conjugates = quaternions * np.array([1.0, -1.0, -1.0, -1.0])
    identity = np.eye(3, dtype=np.float64)[None, :, :]
    operators = (
        orientations.symmetry.operators
        if symmetry_aware and orientations.symmetry is not None
        else identity
    )
    rows_per_block = max(1, int(max_pairs_per_block // count))
    totals = np.zeros(count, dtype=np.float64)
    for start in range(0, count, rows_per_block):
        stop = min(start + rows_per_block, count)
        relative = quaternions_multiply(
            conjugates[start:stop, None, :], quaternions[None, :, :]
        ).reshape(-1, 4)
        angles = _reduced_pair_disorientation_angles_from_quaternions(
            relative, operators, operators
        ).reshape(stop - start, count)
        totals[start:stop] = angles.sum(axis=1)
    minimum = float(totals.min())
    tied = totals <= minimum + tie_rtol * max(abs(minimum), 1.0)
    return int(np.flatnonzero(tied)[0])


def _deduplicate_orientation_set(orientations: OrientationSet) -> OrientationSet:
    keys = np.round(orientations.exact_fundamental_region_keys(), decimals=10)
    _, first_indices = np.unique(keys, axis=0, return_index=True)
    ordered = np.sort(first_indices)
    return orientations.subset(ordered)


def _deterministic_s3_quaternions(count: int) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive.")
    indices = np.arange(count, dtype=np.float64) + 0.5
    golden = (np.sqrt(5.0) - 1.0) / 2.0
    u1 = np.mod(indices * golden, 1.0)
    u2 = np.mod(indices * (np.sqrt(3.0) - 1.0), 1.0)
    u3 = np.mod(indices * (np.sqrt(2.0) - 1.0), 1.0)
    root_a = np.sqrt(1.0 - u1)
    root_b = np.sqrt(u1)
    quaternions = np.column_stack(
        [
            root_b * np.cos(2.0 * np.pi * u3),
            root_a * np.sin(2.0 * np.pi * u2),
            root_a * np.cos(2.0 * np.pi * u2),
            root_b * np.sin(2.0 * np.pi * u3),
        ]
    )
    return _canonicalize_quaternion_rows(quaternions)


def _resolve_phase_symmetry(
    *,
    phase: Phase | None,
    symmetry: SymmetrySpec | None,
    crystal_frame: ReferenceFrame,
) -> tuple[Phase | None, SymmetrySpec | None]:
    resolved_symmetry = symmetry
    if phase is not None:
        if phase.crystal_frame != crystal_frame:
            raise ValueError("phase.crystal_frame must match the target crystal_frame.")
        if resolved_symmetry is None:
            resolved_symmetry = phase.symmetry
    return phase, resolved_symmetry


def _coerce_direction_array(
    direction: ArrayLike,
    *,
    size: int,
    name: str,
) -> np.ndarray:
    array = np.asarray(direction, dtype=np.float64)
    if array.shape == (3,):
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{name} must have shape (3,) or (n, 3).")
    if array.shape[0] == 1 and size > 1:
        array = np.broadcast_to(array, (size, 3))
    elif array.shape[0] != size:
        raise ValueError(f"{name} must broadcast to the number of orientations.")
    return normalize_vectors(array)


def _project_directions_onto_planes(
    directions: np.ndarray,
    normals: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    projected = directions - np.sum(directions * normals, axis=1, keepdims=True) * normals
    norms = np.linalg.norm(projected, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError(
            f"{name} must not be parallel to the corresponding plane normal after projection."
        )
    projected = projected / norms[:, None]
    projected = np.ascontiguousarray(projected, dtype=np.float64)
    projected.setflags(write=False)
    return projected


def _orthonormal_frames_from_normals_and_directions(
    normals: np.ndarray,
    directions: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    x_axis = _project_directions_onto_planes(directions, normals, name=name)
    y_axis = np.cross(normals, x_axis)
    y_norm = np.linalg.norm(y_axis, axis=1)
    if np.any(np.isclose(y_norm, 0.0)):
        raise ValueError(f"{name} does not define a non-degenerate right-handed basis.")
    y_axis = y_axis / y_norm[:, None]
    frames = np.stack([x_axis, y_axis, normals], axis=2)
    frames = np.ascontiguousarray(frames, dtype=np.float64)
    frames.setflags(write=False)
    return frames


def _plane_direction_rotation_matrices(
    *,
    crystal_normals: np.ndarray,
    crystal_directions: np.ndarray,
    specimen_normals: np.ndarray,
    specimen_directions: np.ndarray,
) -> np.ndarray:
    crystal_frames = _orthonormal_frames_from_normals_and_directions(
        crystal_normals,
        crystal_directions,
        name="Crystal directions",
    )
    specimen_frames = _orthonormal_frames_from_normals_and_directions(
        specimen_normals,
        specimen_directions,
        name="Specimen directions",
    )
    matrices = np.einsum("nij,nkj->nik", specimen_frames, crystal_frames, optimize=True)
    matrices = np.ascontiguousarray(matrices, dtype=np.float64)
    matrices.setflags(write=False)
    return matrices


def _phase_from_plane_direction_objects(
    planes: tuple[CrystalPlane, ...],
    directions: tuple[CrystalDirection, ...],
    *,
    phase: Phase | None,
) -> Phase:
    if len(planes) != len(directions):
        raise ValueError("Plane and direction sequences must have the same length.")
    if not planes:
        raise ValueError("At least one plane/direction pair is required.")
    resolved_phase = phase or planes[0].phase
    for plane in planes:
        if plane.phase != resolved_phase:
            raise ValueError("All CrystalPlane inputs must share the same phase.")
    for direction in directions:
        if direction.phase != resolved_phase:
            raise ValueError("All CrystalDirection inputs must share the same phase.")
    return resolved_phase


def _phase_from_miller_objects(
    planes: tuple[MillerPlane | MillerPlaneSet, ...],
    directions: tuple[MillerDirection | MillerDirectionSet, ...],
    *,
    phase: Phase | None,
) -> Phase:
    if len(planes) != len(directions):
        raise ValueError("Plane and direction sequences must have the same length.")
    if not planes:
        raise ValueError("At least one plane/direction pair is required.")
    resolved_phase = phase or planes[0].phase
    for plane in planes:
        if plane.phase != resolved_phase:
            raise ValueError("All Miller-plane inputs must share the same phase.")
    for direction in directions:
        if direction.phase != resolved_phase:
            raise ValueError("All Miller-direction inputs must share the same phase.")
    return resolved_phase


def _coerce_plane_direction_vectors(
    plane: CrystalPlane
    | list[CrystalPlane]
    | tuple[CrystalPlane, ...]
    | MillerPlane
    | MillerPlaneSet
    | ArrayLike,
    direction: CrystalDirection
    | list[CrystalDirection]
    | tuple[CrystalDirection, ...]
    | MillerDirection
    | MillerDirectionSet
    | ArrayLike,
    *,
    phase: Phase | None,
) -> tuple[np.ndarray, np.ndarray, Phase]:
    if isinstance(plane, CrystalPlane):
        if not isinstance(direction, CrystalDirection):
            raise ValueError(
                "A scalar CrystalPlane input requires a matching scalar CrystalDirection."
            )
        resolved_phase = _phase_from_plane_direction_objects((plane,), (direction,), phase=phase)
        return plane.normal[None, :], direction.unit_vector[None, :], resolved_phase
    if (
        isinstance(plane, tuple)
        and isinstance(direction, tuple)
        and all(isinstance(item, CrystalPlane) for item in plane)
        and all(isinstance(item, CrystalDirection) for item in direction)
    ):
        planes = plane
        directions = direction
        resolved_phase = _phase_from_plane_direction_objects(planes, directions, phase=phase)
        return (
            normalize_vectors(np.vstack([item.normal for item in planes])),
            normalize_vectors(np.vstack([item.unit_vector for item in directions])),
            resolved_phase,
        )
    if isinstance(plane, MillerPlane):
        if not isinstance(direction, MillerDirection):
            raise ValueError(
                "A scalar MillerPlane input requires a matching scalar MillerDirection."
            )
        resolved_phase = _phase_from_miller_objects((plane,), (direction,), phase=phase)
        return (
            plane.normal_cartesian[None, :],
            direction.unit_vector_cartesian[None, :],
            resolved_phase,
        )
    if isinstance(plane, MillerPlaneSet):
        if not isinstance(direction, MillerDirectionSet):
            raise ValueError(
                "A MillerPlaneSet input requires a matching MillerDirectionSet."
            )
        if plane.indices.shape[0] != direction.indices.shape[0]:
            raise ValueError("MillerPlaneSet and MillerDirectionSet must have matching lengths.")
        resolved_phase = _phase_from_miller_objects((plane,), (direction,), phase=phase)
        return plane.normals_cartesian(), direction.unit_vectors_cartesian(), resolved_phase

    if phase is None:
        raise ValueError(
            "phase is required when constructing orientations from raw plane and "
            "direction index arrays."
        )
    plane_indices = np.asarray(plane, dtype=np.float64)
    direction_indices = np.asarray(direction, dtype=np.float64)
    if plane_indices.ndim == 1:
        plane_indices = plane_indices[None, :]
    if direction_indices.ndim == 1:
        direction_indices = direction_indices[None, :]
    if plane_indices.ndim != 2 or direction_indices.ndim != 2:
        raise ValueError(
            "Raw plane and direction index arrays must have shape (3,), (4,), (n, 3), or (n, 4)."
        )
    if plane_indices.shape[1] not in {3, 4}:
        raise ValueError("plane index arrays must have 3 or 4 columns.")
    if direction_indices.shape[1] not in {3, 4}:
        raise ValueError("direction index arrays must have 3 or 4 columns.")
    if plane_indices.shape[0] == 1 and direction_indices.shape[0] > 1:
        plane_indices = np.broadcast_to(plane_indices, direction_indices.shape)
    elif direction_indices.shape[0] == 1 and plane_indices.shape[0] > 1:
        direction_indices = np.broadcast_to(direction_indices, plane_indices.shape)
    elif plane_indices.shape[0] != direction_indices.shape[0]:
        raise ValueError("Raw plane and direction index arrays must broadcast to the same length.")
    if plane_indices.shape[1] == 4:
        plane_rows = MillerPlaneSet.from_hkil(plane_indices.astype(np.int64), phase=phase)
        crystal_normals = plane_rows.normals_cartesian()
    else:
        reciprocal_basis = phase.lattice.reciprocal_basis().matrix
        crystal_normals = normalize_vectors(plane_indices @ reciprocal_basis.T)
    if direction_indices.shape[1] == 4:
        direction_rows = MillerDirectionSet.from_UVTW(
            direction_indices.astype(np.int64),
            phase=phase,
        )
        crystal_directions = direction_rows.unit_vectors_cartesian()
    else:
        direct_basis = phase.lattice.direct_basis().matrix
        crystal_directions = normalize_vectors(direction_indices @ direct_basis.T)
    return crystal_normals, crystal_directions, phase


@dataclass(frozen=True, slots=True)
class Rotation:
    """A rotation in three dimensions, stored as a unit quaternion.

    Purpose
    -------
    The frame-agnostic rotation primitive of the library. It carries the
    rotation algebra — composition, inversion, axis/angle, Euler, and
    Rodrigues representations — with no crystallographic meaning attached.
    When the rotation *does* relate a crystal to a specimen, use
    :class:`Orientation`, which adds the frames, symmetry, and phase that
    make the rotation interpretable.

    Convention
    ----------
    Quaternions are stored in ``(w, x, y, z)`` order and normalized on
    construction. Rotations act actively: ``as_matrix()`` returns ``R`` with
    ``v' = R v`` in one fixed frame. Because unit quaternions double-cover
    SO(3), ``q`` and ``-q`` are the same rotation; use :meth:`canonicalized`
    when arrays must compare equal.

    Attributes
    ----------
    quaternion : np.ndarray
        Unit quaternion in ``(w, x, y, z)`` order.
    provenance : ProvenanceRecord, optional
        Where the value came from.
    """

    quaternion: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternion", normalize_quaternion(self.quaternion))

    @classmethod
    def identity(cls) -> Rotation:
        """The null rotation, quaternion ``(1, 0, 0, 0)``.

        Useful as the neutral element of :meth:`compose` and as an explicit
        starting point for accumulated rotations.
        """

        return cls(quaternion=np.array([1.0, 0.0, 0.0, 0.0]))

    @classmethod
    def from_matrix(cls, matrix: ArrayLike) -> Rotation:
        """Rotation from a 3x3 proper-orthogonal matrix.

        The matrix is interpreted actively: it is the matrix that
        :meth:`as_matrix` returns and that :meth:`apply` uses as
        ``v -> R v``. Input that is not a proper rotation (``det = +1``,
        orthonormal columns) is rejected.
        """

        return cls(quaternion=matrix_to_quaternion(matrix))

    @classmethod
    def from_axis_angle(cls, axis: ArrayLike, angle_rad: float) -> Rotation:
        """Rotation of ``angle_rad`` radians about ``axis``, right-handed.

        The axis need not be normalized; only its direction is used. This is
        the natural constructor for misorientation descriptions quoted in the
        literature as an angle about a crystal axis (e.g. a Sigma-3 twin as
        60 deg about ``<111>``).
        """

        return cls(quaternion=quaternion_from_axis_angle(axis, angle_rad))

    @classmethod
    def from_rodrigues(cls, rodrigues: ArrayLike, *, frank: bool = False) -> Rotation:
        """Rotation from Rodrigues (or Rodrigues-Frank) parameters.

        Parameters
        ----------
        rodrigues : ArrayLike
            A 3-vector ``tan(omega/2) * n`` when ``frank`` is ``False``, or a
            homogeneous 4-vector ``(n, tan(omega/2))`` when ``frank`` is
            ``True``. The Frank form keeps the axis separate from the
            magnitude, so it stays exactly invertible at ``omega = pi`` — the
            magnitude is the infinite projective coordinate there, which the
            inverse recognizes — whereas the 3-vector loses the axis in an
            overflowing product.
        frank : bool
            Select the homogeneous four-component form.
        """

        quaternions = quaternions_from_rodrigues(rodrigues, frank=frank)
        return cls(quaternion=quaternions[0])

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
        """Rotation from three Euler angles in a named convention.

        Purpose
        -------
        Build a rotation from the angle triple a texture measurement or a
        literature table quotes, without the caller having to know which axis
        sequence the convention implies.

        Parameters
        ----------
        angle1, angle2, angle3 : float
            The angle triple, in the order the convention names them.
        convention : str
            ``"bunge"`` (ZXZ, ``phi1, Phi, phi2``; the PyTex default and the
            texture-community standard), or ``"matthies"`` / ``"abg"``
            (both ZYZ, ``alpha, beta, gamma``). Aliases such as ``"zxz"`` and
            ``"zyz"`` are accepted.
        degrees : bool
            Interpret the angles as degrees (default) rather than radians.

        See Also
        --------
        from_bunge_euler, from_matthies_euler, from_abg_euler :
            Convention-specific spellings with named parameters.
        to_euler : The inverse.
        """

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
        """Rotation from Bunge ZXZ Euler angles ``(phi1, Phi, phi2)``.

        The convention of Bunge (1982) and of essentially all texture
        software; use this when reading Euler angles out of an EBSD file or a
        published ODF section.
        """

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
        """Rotation from Matthies ZYZ Euler angles ``(alpha, beta, gamma)``.

        The convention used by the Matthies/Roe school of quantitative texture
        analysis; provided so data quoted in that convention need not be
        hand-converted to Bunge angles first.
        """

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
        """Rotation from ZYZ ``(alpha, beta, gamma)`` Euler angles.

        Shares the ZYZ axis sequence with :meth:`from_matthies_euler` and is
        kept as a separate spelling so imported data can retain the name its
        source system used.
        """

        return cls.from_euler(alpha, beta, gamma, convention="abg", degrees=degrees)

    def as_matrix(self) -> np.ndarray:
        """The 3x3 proper-orthogonal matrix of this rotation.

        Active convention: the returned ``R`` maps a vector to its rotated
        image, ``v' = R v``, in a single fixed frame.
        """

        return quaternion_to_matrix(self.quaternion)

    def to_euler(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> tuple[float, float, float]:
        """Euler angles of this rotation in a named convention.

        Purpose
        -------
        Recover the angle triple for reporting, export, or comparison with a
        literature value. The angles are wrapped into ``[0, 2*pi)``; at the
        gimbal-degenerate second angles (``Phi = 0`` or ``pi`` for Bunge) only
        the sum or difference of the outer angles is determined, and PyTex
        resolves the ambiguity by setting the third angle to zero.

        Parameters
        ----------
        convention : str
            ``"bunge"``, ``"matthies"``, or ``"abg"``; see :meth:`from_euler`.
        degrees : bool
            Return degrees (default) rather than radians.

        Returns
        -------
        tuple of float
            The three angles in the order the convention names them.
        """

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
        """Bunge ZXZ Euler angles ``(phi1, Phi, phi2)`` of this rotation."""

        return self.to_euler(convention="bunge", degrees=degrees)

    def to_matthies_euler(self, *, degrees: bool = True) -> tuple[float, float, float]:
        """Matthies ZYZ Euler angles ``(alpha, beta, gamma)`` of this rotation."""

        return self.to_euler(convention="matthies", degrees=degrees)

    def to_abg_euler(self, *, degrees: bool = True) -> tuple[float, float, float]:
        """ZYZ ``(alpha, beta, gamma)`` Euler angles of this rotation."""

        return self.to_euler(convention="abg", degrees=degrees)

    def to_rodrigues(self, *, frank: bool = False) -> np.ndarray:
        """Rodrigues (3-vector) or Rodrigues-Frank (4-vector) parameters.

        Rodrigues space is the natural setting for fundamental-zone geometry,
        because symmetry fundamental zones are convex polyhedra there. Pass
        ``frank=True`` for the homogeneous form, which keeps the axis separate
        from the magnitude and so stays exactly invertible at ``omega = pi``,
        where the 3-vector overflows and loses the axis.
        """

        return as_float_array(
            quaternions_to_rodrigues(self.quaternion[None, :], frank=frank)[0],
            shape=(4 if frank else 3,),
        )

    def compose(self, other: Rotation) -> Rotation:
        """This rotation followed in composition order by ``other``.

        Quaternion multiplication ``q_self * q_other``, so the matrix of the
        result is ``R_self @ R_other``: ``other`` acts on a vector first.
        Composition is not commutative; the order matters and is fixed here.
        """

        return Rotation(quaternion=quaternion_multiply(self.quaternion, other.quaternion))

    def inverse(self) -> Rotation:
        """The inverse rotation (quaternion conjugate).

        ``r.compose(r.inverse())`` is the identity to numerical precision.
        """

        return Rotation(quaternion=quaternion_conjugate(self.quaternion))

    def canonicalized(self) -> Rotation:
        """The same rotation with a sign-canonical quaternion.

        Unit quaternions double-cover SO(3): ``q`` and ``-q`` are the same
        rotation. This picks one representative deterministically, so that
        quaternion arrays can be compared, hashed, or averaged without the
        sign ambiguity leaking into the result.
        """

        return Rotation(quaternion=_canonicalize_quaternion(self.quaternion))

    @property
    def angle_rad(self) -> float:
        """Rotation angle in radians, in ``[0, pi]``.

        Taken from ``|q_w|``, so the double-cover sign never produces an
        angle above ``pi``. This is the raw rotation angle; it is *not*
        symmetry-reduced, so for a crystal misorientation use
        :meth:`Misorientation.disorientation` instead.
        """

        return float(2.0 * _safe_arccos(abs(float(self.quaternion[0]))))

    @property
    def angle_deg(self) -> float:
        """Rotation angle in degrees, in ``[0, 180]``. See :attr:`angle_rad`."""

        return float(np.rad2deg(self.angle_rad))

    @property
    def axis(self) -> np.ndarray:
        """Unit rotation axis in the frame the rotation acts on.

        For a rotation angle of zero the axis is undefined; ``[0, 0, 1]`` is
        returned by convention so downstream code always receives a unit
        vector.
        """

        scalar = float(np.clip(self.quaternion[0], -1.0, 1.0))
        sin_half = np.sqrt(max(0.0, 1.0 - scalar * scalar))
        if np.isclose(sin_half, 0.0):
            return as_float_array([0.0, 0.0, 1.0], shape=(3,))
        return normalize_vector(self.quaternion[1:] / sin_half)

    def distance_to(self, other: Rotation) -> float:
        """Geodesic angle (radians) between two rotations on SO(3).

        The angle of the relative rotation ``other * self^-1``, in
        ``[0, pi]``. This is the bi-invariant Riemannian metric on SO(3) and
        the right notion of "how different" two rotations are. It ignores
        crystal symmetry: for symmetry-aware separation of two crystal
        orientations use :meth:`Orientation.distance_to`.
        """

        return other.compose(self.inverse()).angle_rad

    def apply(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Actively rotate vectors, ``v -> R v``.

        Accepts a single 3-vector, an ``(n, 3)`` array, or a
        :class:`~pytex.core.batches.VectorSet`. A ``VectorSet`` in, a
        ``VectorSet`` out with its reference frame preserved — a bare
        rotation does not change which frame the vectors live in, unlike
        :meth:`Orientation.map_crystal_vector`, which does.
        """

        matrix = self.as_matrix()
        if isinstance(vectors, VectorSet):
            transformed = vectors.values @ matrix.T
            return VectorSet(
                values=transformed,
                reference_frame=vectors.reference_frame,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = array @ matrix.T
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def apply_inverse(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Rotate vectors by the inverse rotation, ``v -> R^T v``.

        Equivalent to ``rotation.inverse().apply(vectors)`` and provided so
        the round trip reads clearly at call sites.
        """

        return self.inverse().apply(vectors)


@dataclass(frozen=True, slots=True)
class Orientation:
    """A crystal orientation: the rotation taking crystal axes to specimen axes.

    Purpose
    -------
    The central object of texture analysis. It is a :class:`Rotation` that
    additionally knows *which* frames it relates, *what* crystal symmetry
    applies, and *which* phase it belongs to — the information without which
    an orientation cannot be compared, reduced, or averaged correctly.

    Convention
    ----------
    Crystal-to-specimen throughout: ``as_matrix()`` returns ``g`` with
    ``v_specimen = g v_crystal``. Misorientations are formed in the crystal
    frame as ``g_1^-1 g_2``, matching MTEX.

    Attributes
    ----------
    rotation : Rotation
        The underlying rotation.
    crystal_frame : ReferenceFrame
        Must belong to the crystal domain; enforced at construction.
    specimen_frame : ReferenceFrame
        Must belong to the specimen domain; enforced at construction.
    symmetry : SymmetrySpec, optional
        Crystal symmetry. Its reference frame must match ``crystal_frame``.
        Without it, no symmetry reduction is applied and orientations that
        are physically identical will compare as different.
    phase : Phase, optional
        Supplies the lattice for Miller-index work; must be consistent with
        ``crystal_frame`` and ``symmetry``.
    provenance : ProvenanceRecord, optional

    See Also
    --------
    OrientationSet : The vectorized batch form, which every map-scale
        calculation should use.
    """

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

    @classmethod
    def from_euler(
        cls,
        angle1: float,
        angle2: float,
        angle3: float,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        convention: str = "bunge",
        degrees: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Crystal orientation from Euler angles.

        Purpose
        -------
        The standard entry point for measured orientation data: EBSD vendors,
        LaboTex, MTEX and the texture literature all quote Euler angles, and
        this constructor attaches the frame, symmetry, and phase meaning that
        a bare angle triple lacks.

        Parameters
        ----------
        angle1, angle2, angle3 : float
            The Euler triple in the order ``convention`` names them.
        crystal_frame : ReferenceFrame, optional
            The crystal-domain frame. Optional only when ``phase`` is given,
            in which case the phase's crystal frame is used.
        specimen_frame : ReferenceFrame
            The specimen-domain frame the orientation maps into (RD/TD/ND for
            rolled sheet, for example). Required, because an orientation is
            meaningless without knowing which specimen axes it refers to.
        symmetry : SymmetrySpec, optional
            Crystal symmetry; inferred from ``phase`` when omitted.
        phase : Phase, optional
            Supplies both crystal frame and symmetry, and lets downstream code
            reach the lattice for Miller-index work.
        convention : str
            ``"bunge"`` (default), ``"matthies"``, or ``"abg"``.
        degrees : bool
            Interpret the angles as degrees (default).
        provenance : ProvenanceRecord, optional
            Where the value came from; carried through derived objects.

        Returns
        -------
        Orientation
            The crystal-to-specimen orientation, not symmetry-reduced. Call
            :meth:`project_to_fundamental_region` if you need a canonical
            representative.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            rotation=Rotation.from_euler(
                angle1,
                angle2,
                angle3,
                convention=convention,
                degrees=degrees,
            ),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_axis_angle(
        cls,
        axis: str | ArrayLike,
        angle_rad: float,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Crystal orientation from a rotation axis and angle.

        ``axis`` may be a vector in the specimen frame or one of the named
        specimen directions ``"RD"``, ``"TD"``, ``"ND"``, ``"x"``, ``"y"``,
        ``"z"``, which makes rotations about the sheet axes readable at the
        call site. The remaining arguments carry the meaning given in
        :meth:`from_euler`.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        axis_vector = specimen_direction_vector(axis) if isinstance(axis, str) else axis
        return cls(
            rotation=Rotation.from_axis_angle(axis_vector, angle_rad),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_matrix(
        cls,
        matrix: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Crystal orientation from a 3x3 crystal-to-specimen rotation matrix.

        Use this when the orientation is already available as the matrix ``g`` —
        from an external tool, or after building it from two pairs of parallel
        directions. The matrix must be proper orthogonal. Remaining arguments
        carry the meaning given in :meth:`from_euler`.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            rotation=Rotation.from_matrix(matrix),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_quaternion(
        cls,
        quaternion: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Crystal orientation from a unit quaternion in ``(w, x, y, z)`` order.

        The storage order is fixed repository-wide; see
        ``docs/standards/notation_and_conventions.md``. The quaternion is
        normalized on construction, and its sign is not canonicalized — ``q``
        and ``-q`` describe the same orientation but store different arrays.
        Remaining arguments carry the meaning given in :meth:`from_euler`.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            rotation=Rotation(quaternion=np.asarray(quaternion, dtype=np.float64)),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_rodrigues(
        cls,
        rodrigues: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        frank: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Crystal orientation from Rodrigues or Rodrigues-Frank parameters.

        Pass ``frank=True`` for the homogeneous 4-vector form, which stays
        finite at a rotation angle of ``pi`` where the 3-vector diverges.
        Remaining arguments carry the meaning given in :meth:`from_euler`.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            rotation=Rotation.from_rodrigues(rodrigues, frank=frank),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_plane_direction(
        cls,
        plane: CrystalPlane,
        direction: CrystalDirection,
        *,
        specimen_frame: ReferenceFrame,
        specimen_plane_normal: ArrayLike = (0.0, 0.0, 1.0),
        specimen_direction: ArrayLike = (1.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Orientation from a crystal plane and an in-plane crystal direction.

        Purpose
        -------
        Build the orientation that a rolling-texture component is *named* by:
        ``(hkl)[uvw]`` means the plane ``(hkl)`` lies in the sheet plane and the
        direction ``[uvw]`` points along the rolling direction. This constructor
        turns that statement into a concrete orientation.

        Parameters
        ----------
        plane : CrystalPlane
            The crystal plane whose normal is aligned to ``specimen_plane_normal``.
        direction : CrystalDirection
            An in-plane crystal direction aligned to ``specimen_direction``. It
            need not be exactly perpendicular to the plane normal; the component
            along the normal is removed, so a slightly inconsistent pair still
            yields a proper rotation.
        specimen_frame : ReferenceFrame
            The specimen-domain target frame.
        specimen_plane_normal, specimen_direction : ArrayLike
            The specimen-frame axes that the plane normal and the direction are
            aligned to; the defaults give the usual sheet convention of the plane
            normal along ND and the direction along RD.

        Returns
        -------
        Orientation
            Carrying the phase, symmetry, and crystal frame of ``plane``.

        See Also
        --------
        from_miller : The same construction accepting raw index triples and
            named specimen directions.
        """

        if plane.phase != direction.phase:
            raise ValueError("plane.phase must match direction.phase.")
        matrices = _plane_direction_rotation_matrices(
            crystal_normals=plane.normal[None, :],
            crystal_directions=direction.unit_vector[None, :],
            specimen_normals=_coerce_direction_array(
                specimen_plane_normal,
                size=1,
                name="specimen_plane_normal",
            ),
            specimen_directions=_coerce_direction_array(
                specimen_direction,
                size=1,
                name="specimen_direction",
            ),
        )
        return cls(
            rotation=Rotation.from_matrix(matrices[0]),
            crystal_frame=plane.phase.crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=plane.phase.symmetry,
            phase=plane.phase,
            provenance=provenance,
        )

    @classmethod
    def from_miller(
        cls,
        plane: CrystalPlane | ArrayLike,
        direction: CrystalDirection | ArrayLike,
        *,
        specimen_frame: ReferenceFrame,
        phase: Phase | None = None,
        specimen_plane_normal: str | ArrayLike = "ND",
        specimen_direction: str | ArrayLike = "RD",
        provenance: ProvenanceRecord | None = None,
    ) -> Orientation:
        """Orientation from ``(hkl)[uvw]`` indices, the texture-component spelling.

        Purpose
        -------
        The convenience form of :meth:`from_plane_direction`: it accepts bare
        index triples together with a ``phase``, and accepts the specimen axes
        by name, so a named texture component transcribes directly. The copper
        component ``{112}<111>`` on a cubic phase, for example, is
        ``Orientation.from_miller((1, 1, 2), (1, 1, -1), phase=phase,
        specimen_frame=frame)``.

        Parameters
        ----------
        plane, direction : CrystalPlane / CrystalDirection or ArrayLike
            Either the typed objects, or index triples together with ``phase``.
        specimen_frame : ReferenceFrame
            The specimen-domain target frame.
        phase : Phase, optional
            Required when the indices are passed as raw arrays.
        specimen_plane_normal, specimen_direction : str or ArrayLike
            Specimen axes, by name (``"ND"``, ``"RD"``, ``"TD"``, ``"x"``,
            ``"y"``, ``"z"``) or as vectors. Defaults are ``ND`` and ``RD``.

        Returns
        -------
        Orientation
            The single orientation described by the index pair.
        """

        orientations = OrientationSet.from_plane_direction(
            plane,
            direction,
            specimen_frame=specimen_frame,
            phase=phase,
            specimen_plane_normal=specimen_direction_vector(specimen_plane_normal),
            specimen_direction=specimen_direction_vector(specimen_direction),
            provenance=provenance,
        )
        if len(orientations) != 1:
            raise ValueError("Orientation.from_miller requires scalar plane and direction inputs.")
        return orientations[0]

    def as_matrix(self) -> np.ndarray:
        """The 3x3 crystal-to-specimen orientation matrix ``g``.

        ``g`` maps a vector expressed in the crystal frame to the same physical
        vector expressed in the specimen frame, matching the repository-wide
        orientation-as-crystal-to-specimen convention.
        """

        return self.rotation.as_matrix()

    def map_crystal_vector(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Express a crystal-frame vector in the specimen frame.

        Purpose
        -------
        The forward action of the orientation: it answers "where does this
        crystal direction point in the specimen?" — the computation behind
        every pole figure.

        Parameters
        ----------
        vector : ArrayLike or VectorSet
            A 3-vector, or a ``VectorSet`` whose reference frame must equal
            :attr:`crystal_frame` (checked, not assumed).

        Returns
        -------
        np.ndarray or VectorSet
            The same physical vector in :attr:`specimen_frame`. A ``VectorSet``
            in gives a ``VectorSet`` out, re-framed to the specimen domain.

        See Also
        --------
        map_sample_vector_to_crystal : The inverse direction.
        """

        if isinstance(vector, VectorSet):
            if vector.reference_frame != self.crystal_frame:
                raise ValueError("VectorSet.reference_frame must match Orientation.crystal_frame.")
            return VectorSet(
                values=vector.values @ self.rotation.as_matrix().T,
                reference_frame=self.specimen_frame,
                provenance=vector.provenance,
            )
        return as_float_array(
            np.asarray(vector, dtype=np.float64) @ self.rotation.as_matrix().T,
            shape=(3,),
        )

    def map_sample_vector_to_crystal(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Express a specimen-frame vector in the crystal frame.

        The inverse of :meth:`map_crystal_vector`, and the computation behind
        every inverse pole figure: it answers "which crystal direction is
        parallel to this specimen axis?". A ``VectorSet`` argument must carry
        :attr:`specimen_frame` and is returned re-framed to the crystal domain.
        """

        if isinstance(vector, VectorSet):
            if vector.reference_frame != self.specimen_frame:
                raise ValueError("VectorSet.reference_frame must match Orientation.specimen_frame.")
            inverse = self.rotation.inverse().as_matrix()
            return VectorSet(
                values=vector.values @ inverse.T,
                reference_frame=self.crystal_frame,
                provenance=vector.provenance,
            )
        inverse = self.rotation.inverse().as_matrix()
        return as_float_array(np.asarray(vector, dtype=np.float64) @ inverse.T, shape=(3,))

    def equivalent_orientations(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> OrientationSet:
        """All symmetrically equivalent expressions of this orientation.

        Purpose
        -------
        Enumerate the orbit of ``g`` under the crystal symmetry operators on the
        right and, optionally, the specimen (statistical) symmetry operators on
        the left. Every member describes the same physical crystal; which one a
        measurement happens to report is arbitrary, which is exactly why
        orientation comparisons must be symmetry-aware.

        Parameters
        ----------
        specimen_symmetry : SymmetrySpec, optional
            Specimen symmetry to include on the left (for example orthorhombic
            sample symmetry for a rolled sheet). Its reference frame must match
            :attr:`specimen_frame`. Omitting it means triclinic specimen
            symmetry — the honest default for a single measured orientation.

        Returns
        -------
        OrientationSet
            ``n_specimen * n_crystal`` orientations, unreduced and in
            enumeration order.
        """

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
        """The sign-canonical representative of the symmetry orbit.

        Chooses one member of :meth:`equivalent_orientations` by a fixed
        quaternion-ordering rule, so equal orientations always yield equal
        arrays. Use it for hashing, grouping, or deduplication. For the
        representative that is also *geometrically* canonical — the one lying
        in the fundamental region — use :meth:`project_to_fundamental_region`.
        """

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
        """A hashable key identifying the symmetry class of this orientation.

        Two symmetrically equivalent orientations produce the same key, which is
        derived from the quaternion of the fundamental-region representative.
        Intended for grouping and set membership, not for metric comparison —
        nearby orientations do not have nearby keys.
        """

        projected = self.project_to_exact_fundamental_region(specimen_symmetry=specimen_symmetry)
        return _fundamental_region_key(projected.rotation, projected.symmetry)

    def project_to_fundamental_region(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
        reference_orientation: Orientation | None = None,
    ) -> Orientation:
        """The representative of this orientation inside the fundamental region.

        Purpose
        -------
        Reduce the symmetry-equivalent orbit to one canonical member, so that
        orientations can be compared, averaged, or plotted without the arbitrary
        choice of symmetry operator polluting the result.

        Parameters
        ----------
        specimen_symmetry : SymmetrySpec, optional
            Specimen symmetry to include in the orbit; see
            :meth:`equivalent_orientations`.
        reference_orientation : Orientation, optional
            When given, pick the orbit member *closest to this reference* rather
            than the globally canonical one. This is what makes grain reference
            orientation deviation (GROD) and grain averaging behave: every member
            of a grain is expressed in the branch nearest the grain reference, so
            no spurious symmetry jumps appear across the grain.

        Returns
        -------
        Orientation
            The chosen representative, carrying the original frames, symmetry,
            phase, and provenance.
        """

        return self.project_to_exact_fundamental_region(
            specimen_symmetry=specimen_symmetry,
            reference_orientation=reference_orientation,
        )

    def project_to_exact_fundamental_region(
        self,
        specimen_symmetry: SymmetrySpec | None = None,
        reference_orientation: Orientation | None = None,
    ) -> Orientation:
        """Exact fundamental-region projection by explicit orbit enumeration.

        Identical in contract to :meth:`project_to_fundamental_region`, which
        delegates here. The name records that the reduction is exact — every
        symmetry operator is applied and the best member selected — rather than
        approximated by inequality tests on Rodrigues coordinates.
        """

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
        """Whether this orientation already is its own canonical representative.

        Compares the sign-canonical quaternion against the one returned by
        :meth:`project_to_exact_fundamental_region`. Useful as a test-suite
        assertion and to skip redundant reductions.

        Parameters
        ----------
        specimen_symmetry : SymmetrySpec, optional
            Specimen symmetry to include in the orbit.
        atol : float
            Absolute tolerance on the quaternion comparison.
        """

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
        """The misorientation from this orientation to ``other``.

        Purpose
        -------
        The crystal-frame rotation relating two grains — the quantity behind
        grain-boundary character, twin identification, and
        orientation-relationship analysis.

        Convention
        ----------
        Computed as ``g_self^-1 . g_other`` in the crystal frame, matching MTEX.
        Both orientations must share crystal frame, specimen frame, and (when
        both declare one) phase; a mismatch raises rather than silently
        producing a meaningless rotation.

        Parameters
        ----------
        other : Orientation
            The second orientation.
        reduce_by_symmetry : bool
            When ``True`` (default) the returned object is reduced to the
            disorientation — the symmetry-equivalent representative with the
            smallest rotation angle. When ``False`` the raw relative rotation is
            kept.

        Returns
        -------
        Misorientation
            Carrying the crystal frame and the symmetry needed to interpret it.
        """

        if self.crystal_frame != other.crystal_frame:
            raise ValueError("Misorientation requires the same crystal frame.")
        if self.specimen_frame != other.specimen_frame:
            raise ValueError("Misorientation requires the same specimen frame.")
        if self.phase is not None and other.phase is not None and self.phase != other.phase:
            raise ValueError("Misorientation requires matching phases when both are specified.")
        # Crystal-frame misorientation (MTEX convention inv(o1) * o2): crystal
        # symmetry then acts as fixed left/right operator products, which is
        # what Misorientation.disorientation() enumerates.
        delta = self.rotation.inverse().compose(other.rotation)
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
        """Angular distance in radians between two crystal orientations.

        With ``symmetry_aware=True`` (the default) this is the disorientation
        angle: the smallest rotation angle over all symmetry equivalents, in
        ``[0, pi]`` and bounded by the symmetry's maximum disorientation
        (62.8 degrees for cubic-cubic). With ``symmetry_aware=False`` it is the
        plain SO(3) geodesic angle, which is only meaningful for orientations
        already reduced to a common symmetry branch.
        """

        return self.misorientation_to(other, reduce_by_symmetry=symmetry_aware).angle_rad


@dataclass(frozen=True, slots=True)
class Misorientation:
    """The rotation relating two crystal orientations, in the crystal frame.

    Purpose
    -------
    The quantity behind grain-boundary character, twin identification, CSL
    classification, and orientation-relationship analysis. It carries the
    symmetry groups of both sides, which is what allows it to be reduced to
    the physically meaningful *disorientation*.

    Convention
    ----------
    Formed as ``g_1^-1 g_2`` in the crystal frame. The stored rotation is not
    necessarily reduced; call :meth:`disorientation` for the
    symmetry-minimal representative, which is what
    :meth:`Orientation.misorientation_to` returns by default.

    Attributes
    ----------
    rotation : Rotation
        The relative rotation.
    left_symmetry, right_symmetry : SymmetrySpec, optional
        Symmetry groups of the two crystals. Both are needed for reduction;
        they differ for an interphase misorientation.
    provenance : ProvenanceRecord, optional
    """

    rotation: Rotation
    left_symmetry: SymmetrySpec | None = None
    right_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def as_matrix(self) -> np.ndarray:
        """The 3x3 matrix of the misorientation, in the crystal frame.
        """

        return self.rotation.as_matrix()

    @property
    def angle_rad(self) -> float:
        """Misorientation angle in radians, in ``[0, pi]``.

        This is the angle of the *stored* representative. Call
        :meth:`disorientation` first if you need the symmetry-reduced minimum;
        :meth:`Orientation.misorientation_to` already does so by default.
        """

        return self.rotation.angle_rad

    @property
    def angle_deg(self) -> float:
        """Misorientation angle in degrees. See :attr:`angle_rad`.
        """

        return self.rotation.angle_deg

    def disorientation(self) -> Misorientation:
        """The symmetry-reduced representative of this misorientation.

        Purpose
        -------
        A misorientation between two crystals is defined only up to the symmetry
        of both: the physically meaningful quantity is the *disorientation*, the
        orbit member with the smallest rotation angle, which lies in the
        misorientation fundamental zone. Grain-boundary character, twin
        identification, and CSL classification are all statements about the
        disorientation, not about an arbitrary representative.

        Method
        ------
        Enumerates ``S_left . m . S_right^T`` over both symmetry groups and
        selects the member with the canonical fundamental-region key, which for
        equal-angle ties picks a deterministic representative rather than one
        decided by floating-point noise.

        Returns
        -------
        Misorientation
            The reduced representative, carrying the same symmetry groups. For
            cubic-cubic symmetry its angle never exceeds 62.8 degrees.
        """

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
    """A batch of crystal orientations sharing frames, symmetry, and phase.

    Purpose
    -------
    The vectorized counterpart of :class:`Orientation`, and the form every
    map-scale calculation should use: an EBSD scan is one ``OrientationSet``,
    not a list of a million ``Orientation`` objects. Symmetry reduction,
    misorientation matrices, means, and spreads are all implemented as array
    operations over the stored quaternions.

    Because frames, symmetry, and phase are shared by the whole set rather
    than repeated per element, they are stated once and cannot drift between
    members.

    Attributes
    ----------
    quaternions : np.ndarray
        ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order.
    crystal_frame, specimen_frame : ReferenceFrame
        Shared domain-typed frames.
    symmetry : SymmetrySpec, optional
        Shared crystal symmetry.
    phase : Phase, optional
        Shared phase.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Supports ``len()``, integer indexing (yielding an ``Orientation``), and
    slicing (yielding an ``OrientationSet``).
    """

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
        """Build a set from individual :class:`Orientation` objects.

        All members must agree on crystal frame, specimen frame, symmetry, and
        phase; a mismatch raises rather than silently dropping the metadata of
        some members. Use this to collect scattered orientations into the
        vectorized batch form that the texture and EBSD layers expect.
        """

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

    @overload
    def __getitem__(self, index: int) -> Orientation: ...

    @overload
    def __getitem__(self, index: slice) -> OrientationSet: ...

    def __getitem__(self, index: int | slice) -> Orientation | OrientationSet:
        if isinstance(index, slice):
            return OrientationSet(
                quaternions=self.quaternions[index],
                crystal_frame=self.crystal_frame,
                specimen_frame=self.specimen_frame,
                symmetry=self.symmetry,
                phase=self.phase,
                provenance=self.provenance,
            )
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
        angles: ArrayLike | EulerSet,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        convention: str = "bunge",
        degrees: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from an ``(n, 3)`` array of Euler angles.

        Purpose
        -------
        The bulk counterpart of :meth:`Orientation.from_euler`, and the usual
        entry point for a whole EBSD scan or a measured orientation list.

        Parameters
        ----------
        angles : ArrayLike
            ``(n, 3)`` Euler triples in the order the convention names them.
        crystal_frame : ReferenceFrame, optional
            Required unless ``phase`` supplies it.
        specimen_frame : ReferenceFrame
            The specimen-domain frame shared by every member.
        symmetry : SymmetrySpec, optional
            Crystal symmetry; inferred from ``phase`` when omitted.
        phase : Phase, optional
            Supplies crystal frame and symmetry.
        convention : str
            ``"bunge"`` (default), ``"matthies"``, or ``"abg"``.
        degrees : bool
            Interpret the angles as degrees (default).
        provenance : ProvenanceRecord, optional
            Carried onto the set and its derivatives.

        Returns
        -------
        OrientationSet
            ``n`` orientations in input order, not symmetry-reduced.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        if isinstance(angles, EulerSet):
            angle_array = angles.angles
            convention = angles.convention
            degrees = angles.degrees
            if provenance is None:
                provenance = angles.provenance
        else:
            angle_array = as_float_array(angles, shape=(None, 3))
        quaternions = RotationSet.from_euler_set(
            EulerSet(
                angles=angle_array,
                convention=convention,
                degrees=degrees,
                provenance=provenance,
            )
        ).quaternions
        return cls(
            quaternions=np.asarray(quaternions, dtype=np.float64),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_quaternions(
        cls,
        quaternions: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from an ``(n, 4)`` array of quaternions in ``(w, x, y, z)``.

        Quaternions are normalized on construction. Signs are left as given, so
        if you need comparable arrays call :meth:`canonicalized` afterwards.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            quaternions=np.asarray(quaternions, dtype=np.float64),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_axes_angles(
        cls,
        axes: ArrayLike,
        angles_rad: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from ``(n, 3)`` rotation axes and ``(n,)`` angles in radians.

        Axes need not be normalized. Useful for constructing controlled
        orientation spreads about a fixed axis, and for reproducing
        literature statements of the form "rotated by omega about ``<uvw>``".
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            quaternions=RotationSet.from_axes_angles(axes, angles_rad).quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_matrices(
        cls,
        matrices: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from an ``(n, 3, 3)`` stack of crystal-to-specimen matrices.

        Each matrix must be proper orthogonal. This is the natural bridge from
        external code that produces orientation matrices directly.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            quaternions=RotationSet.from_matrices(matrices).quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_rodrigues(
        cls,
        rodrigues: ArrayLike,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        frank: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from Rodrigues ``(n, 3)`` or Rodrigues-Frank ``(n, 4)`` rows.

        Pass ``frank=True`` for the homogeneous form, which keeps the axis
        separate from the magnitude and so stays exactly invertible at a
        rotation angle of ``pi``.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        return cls(
            quaternions=RotationSet.from_rodrigues(
                rodrigues,
                frank=frank,
                provenance=provenance,
            ).quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=provenance,
        )

    @classmethod
    def from_miller(
        cls,
        plane: MillerPlane | MillerPlaneSet | ArrayLike,
        direction: MillerDirection | MillerDirectionSet | ArrayLike,
        *,
        specimen_frame: ReferenceFrame,
        phase: Phase | None = None,
        specimen_plane_normal: str | ArrayLike = "ND",
        specimen_direction: str | ArrayLike = "RD",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from ``(hkl)[uvw]`` index pairs, the texture-component form.

        The bulk counterpart of :meth:`Orientation.from_miller`: accepts arrays
        of plane and direction indices together with a ``phase``, and accepts
        the specimen axes by name. Use it to instantiate a list of named ideal
        components (copper, brass, cube, Goss, ...) in one call.
        """

        return cls.from_plane_direction(
            plane=plane,
            direction=direction,
            specimen_frame=specimen_frame,
            phase=phase,
            specimen_plane_normal=specimen_direction_vector(specimen_plane_normal),
            specimen_direction=specimen_direction_vector(specimen_direction),
            provenance=provenance,
        )

    @classmethod
    def from_plane_direction(
        cls,
        plane: CrystalPlane
        | MillerPlane
        | MillerPlaneSet
        | list[CrystalPlane]
        | tuple[CrystalPlane, ...]
        | ArrayLike,
        direction: CrystalDirection
        | MillerDirection
        | MillerDirectionSet
        | list[CrystalDirection]
        | tuple[CrystalDirection, ...]
        | ArrayLike,
        *,
        specimen_frame: ReferenceFrame,
        phase: Phase | None = None,
        specimen_plane_normal: ArrayLike = (0.0, 0.0, 1.0),
        specimen_direction: ArrayLike = (1.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Build a set from crystal plane / in-plane direction pairs.

        The bulk counterpart of :meth:`Orientation.from_plane_direction`. Each
        row aligns its plane normal to ``specimen_plane_normal`` and its
        direction to ``specimen_direction``; the direction's component along the
        normal is removed, so slightly inconsistent index pairs still give
        proper rotations.
        """

        if (
            isinstance(plane, list)
            and isinstance(direction, list)
            and all(isinstance(item, CrystalPlane) for item in plane)
            and all(isinstance(item, CrystalDirection) for item in direction)
        ):
            plane = tuple(plane)
            direction = tuple(direction)
        crystal_normals, crystal_directions, resolved_phase = _coerce_plane_direction_vectors(
            plane,
            direction,
            phase=phase,
        )
        specimen_count = int(crystal_normals.shape[0])
        matrices = _plane_direction_rotation_matrices(
            crystal_normals=crystal_normals,
            crystal_directions=crystal_directions,
            specimen_normals=_coerce_direction_array(
                specimen_plane_normal,
                size=specimen_count,
                name="specimen_plane_normal",
            ),
            specimen_directions=_coerce_direction_array(
                specimen_direction,
                size=specimen_count,
                name="specimen_direction",
            ),
        )
        return cls.from_matrices(
            matrices,
            crystal_frame=resolved_phase.crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=resolved_phase.symmetry,
            phase=resolved_phase,
            provenance=provenance,
        )

    @classmethod
    def from_bunge_grid(
        cls,
        *,
        crystal_frame: ReferenceFrame,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        phi1_range_deg: tuple[float, float] = (0.0, 360.0),
        big_phi_range_deg: tuple[float, float] = (0.0, 90.0),
        phi2_range_deg: tuple[float, float] = (0.0, 90.0),
        phi1_step_deg: float = 15.0,
        big_phi_step_deg: float = 15.0,
        phi2_step_deg: float = 15.0,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """A regular orientation grid in Bunge Euler space.

        Purpose
        -------
        Generate the evaluation support for a discrete ODF, for phi2 sections,
        or for any calculation that needs orientation space sampled on the same
        Euler grid the texture literature plots.

        Parameters
        ----------
        crystal_frame, specimen_frame : ReferenceFrame
            Frames shared by every grid point.
        symmetry, phase : optional
            Crystal symmetry and phase for the grid points.
        phi1_range_deg, big_phi_range_deg, phi2_range_deg : tuple of float
            Inclusive angular ranges in degrees. The defaults cover the
            orthorhombic-specimen cubic-crystal asymmetric domain
            (``0-360``, ``0-90``, ``0-90``).
        phi1_step_deg, big_phi_step_deg, phi2_step_deg : float
            Grid spacing per axis.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        OrientationSet
            The Cartesian product of the three axis grids, flattened in
            ``phi1``-major order. Note that a regular Euler grid is *not*
            equal-volume in SO(3) — cells shrink as ``sin(Phi)`` toward
            ``Phi = 0``. Use :meth:`from_equispaced_so3_grid` when uniform
            sampling density matters.
        """

        phi1_values = _grid_axis_values(
            float(phi1_range_deg[0]),
            float(phi1_range_deg[1]),
            float(phi1_step_deg),
            periodic=True,
        )
        big_phi_values = _grid_axis_values(
            float(big_phi_range_deg[0]),
            float(big_phi_range_deg[1]),
            float(big_phi_step_deg),
            periodic=False,
        )
        phi2_values = _grid_axis_values(
            float(phi2_range_deg[0]),
            float(phi2_range_deg[1]),
            float(phi2_step_deg),
            periodic=True,
        )
        phi1_mesh, big_phi_mesh, phi2_mesh = np.meshgrid(
            phi1_values,
            big_phi_values,
            phi2_values,
            indexing="ij",
        )
        angles = np.column_stack(
            [
                phi1_mesh.reshape(-1),
                big_phi_mesh.reshape(-1),
                phi2_mesh.reshape(-1),
            ]
        )
        return cls.from_euler_angles(
            angles,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            convention="bunge",
            degrees=True,
            provenance=provenance,
        )

    @classmethod
    def from_so2_grid(
        cls,
        axis: str | ArrayLike = "ND",
        spacing_deg: float = 5.0,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """Orientations generated by rotation about a single specimen axis.

        Purpose
        -------
        The one-parameter (fibre) family obtained by spinning the crystal about
        a fixed specimen direction — the support of a fibre texture, and a
        convenient controlled test set.

        Parameters
        ----------
        axis : str or ArrayLike
            Specimen axis, by name (``"ND"`` default, ``"RD"``, ``"TD"``,
            ``"x"``, ``"y"``, ``"z"``) or as a vector.
        spacing_deg : float
            Angular step in degrees over the full ``[0, 360)`` turn.

        Returns
        -------
        OrientationSet
            ``360 / spacing_deg`` orientations, provenance-tagged with the
            generation method and spacing.
        """

        spacing = _require_grid_spacing(spacing_deg)
        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        angles_deg = np.arange(0.0, 360.0, spacing, dtype=np.float64)
        if angles_deg.size == 0:
            raise ValueError("SO2 grid generation produced no support orientations.")
        axis_vector = specimen_direction_vector(axis)
        quaternions = RotationSet.from_axes_angles(
            np.broadcast_to(axis_vector, (angles_deg.shape[0], 3)),
            np.deg2rad(angles_deg),
        ).quaternions
        return cls(
            quaternions=quaternions,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=_orientation_grid_provenance(
                method="so2_axis_grid",
                spacing_deg=spacing,
                note=f"SO2 grid around specimen axis at {spacing:g} degree spacing.",
                provenance=provenance,
            ),
        )

    @classmethod
    def from_regular_so3_grid(
        cls,
        spacing_deg: float,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        phi1_range_deg: tuple[float, float] = (0.0, 360.0),
        big_phi_range_deg: tuple[float, float] = (0.0, 180.0),
        phi2_range_deg: tuple[float, float] = (0.0, 360.0),
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """A regular SO(3) grid at a single angular spacing, built in Euler space.

        Convenience wrapper over :meth:`from_bunge_grid` that applies one
        spacing to all three Euler axes and defaults to the full SO(3) ranges.
        The sampling density is non-uniform in SO(3) for the reason given in
        :meth:`from_bunge_grid`; the provenance record states that the grid is
        MTEX-inspired and is not a parity claim.
        """

        spacing = _require_grid_spacing(spacing_deg)
        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        return cls.from_bunge_grid(
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            phi1_range_deg=phi1_range_deg,
            big_phi_range_deg=big_phi_range_deg,
            phi2_range_deg=phi2_range_deg,
            phi1_step_deg=spacing,
            big_phi_step_deg=spacing,
            phi2_step_deg=spacing,
            provenance=_orientation_grid_provenance(
                method="regular_so3_bunge_grid",
                spacing_deg=spacing,
                note=(
                    "Regular SO3 grid generated in Bunge Euler coordinates; "
                    "MTEX-inspired, not a parity claim."
                ),
                provenance=provenance,
            ),
        )

    @classmethod
    def from_equispaced_so3_grid(
        cls,
        spacing_deg: float,
        *,
        crystal_frame: ReferenceFrame | None = None,
        specimen_frame: ReferenceFrame,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
        reduce_to_fundamental_region: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationSet:
        """A near-uniform SO(3) grid at the requested angular resolution.

        Purpose
        -------
        Sample orientation space with approximately constant density — what a
        Monte-Carlo integration, a dictionary-indexing search, or an unbiased
        ODF support needs, and what a regular Euler grid cannot provide.

        Method
        ------
        A deterministic low-discrepancy point set on the unit quaternion sphere
        ``S^3``, sized from the requested spacing by the SO(3) volume relation
        ``n ~ 8 * pi^2 / spacing^3``. Deterministic means the same spacing always
        yields the same grid, so results are reproducible.

        Parameters
        ----------
        spacing_deg : float
            Target angular resolution in degrees. Cost grows as the inverse cube.
        crystal_frame, specimen_frame, symmetry, phase : see :meth:`from_euler_angles`.
        reduce_to_fundamental_region : bool
            When ``True`` (default) and a symmetry is known, the grid is
            projected into the fundamental region and deduplicated, so the
            returned count is roughly ``n / |symmetry|``.

        Returns
        -------
        OrientationSet
            Provenance-tagged with the generation method and spacing. The
            construction is MTEX-inspired and is not a parity claim.
        """

        spacing = _require_grid_spacing(spacing_deg)
        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        phase, symmetry = _resolve_phase_symmetry(
            phase=phase,
            symmetry=symmetry,
            crystal_frame=crystal_frame,
        )
        spacing_rad = np.deg2rad(spacing)
        count = max(1, int(np.ceil((8.0 * np.pi * np.pi) / (spacing_rad**3))))
        orientations = cls(
            quaternions=_deterministic_s3_quaternions(count),
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            provenance=_orientation_grid_provenance(
                method="equispaced_so3_quaternion_grid",
                spacing_deg=spacing,
                note=(
                    "Deterministic low-discrepancy SO3 quaternion grid; "
                    "MTEX-inspired, not a parity claim."
                ),
                provenance=provenance,
            ),
        )
        if reduce_to_fundamental_region and symmetry is not None:
            orientations = orientations.projected_to_fundamental_region()
            orientations = _deduplicate_orientation_set(orientations)
        return orientations

    def as_matrices(self) -> np.ndarray:
        """``(n, 3, 3)`` stack of crystal-to-specimen orientation matrices, read-only.
        """

        matrices = np.ascontiguousarray(quaternions_to_matrices(self.quaternions))
        matrices.setflags(write=False)
        return matrices

    def to_axes_angles(self) -> tuple[np.ndarray, np.ndarray]:
        """Rotation axes and angles of every member.

        Returns
        -------
        tuple of np.ndarray
            ``(n, 3)`` unit axes and ``(n,)`` angles in radians, unreduced by
            symmetry.
        """

        return self.as_rotation_set().to_axes_angles()

    def as_euler(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> np.ndarray:
        """``(n, 3)`` Euler angles in a named convention.

        Angles are wrapped into ``[0, 2*pi)`` (or ``[0, 360)`` in degrees). At the
        gimbal-degenerate second angle only the sum or difference of the outer
        angles is determined; PyTex resolves this by setting the third angle to
        zero. See :meth:`Rotation.to_euler`.
        """

        raw = _matrices_to_repeated_axis_euler(self.as_matrices(), convention=convention)
        euler = np.mod(raw, 2.0 * np.pi)
        if degrees:
            euler = np.rad2deg(euler)
        euler = np.ascontiguousarray(euler)
        euler.setflags(write=False)
        return euler

    def as_bunge_euler(self, *, degrees: bool = True) -> np.ndarray:
        """``(n, 3)`` Bunge ZXZ Euler angles ``(phi1, Phi, phi2)``.

        The array form expected by EBSD exporters and by the texture literature.
        """

        return self.as_euler(convention="bunge", degrees=degrees)

    def as_euler_set(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> EulerSet:
        """The Euler angles as a typed :class:`~pytex.core.batches.EulerSet`.

        Prefer this over :meth:`as_euler` on stable surfaces: the ``EulerSet``
        carries the convention and the degrees/radians flag with the numbers, so
        the array cannot be reinterpreted under the wrong convention downstream.
        """

        return EulerSet(
            angles=self.as_euler(convention=convention, degrees=degrees),
            convention=convention,
            degrees=degrees,
            provenance=self.provenance,
        )

    def as_rotation_set(self) -> RotationSet:
        """The underlying rotations, stripped of crystal/specimen frame meaning.

        Use this only for frame-agnostic rotation algebra. The result no longer
        knows it maps crystal to specimen, so it will not raise if it is later
        applied in the wrong direction.
        """

        return RotationSet(quaternions=self.quaternions, provenance=self.provenance)

    def map_crystal_directions(self, directions: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Express crystal directions in the specimen frame, for the whole set.

        Purpose
        -------
        The vectorized action behind pole-figure construction: for each
        orientation, where does the given crystal direction point in the
        specimen?

        Parameters
        ----------
        directions : ArrayLike or VectorSet
            Either one ``(3,)`` direction applied to every orientation, or an
            ``(n, 3)`` array applied row-wise (one direction per orientation).
            A ``VectorSet`` must carry :attr:`crystal_frame`.

        Returns
        -------
        np.ndarray or VectorSet
            ``(n, 3)`` unit vectors in the specimen frame. A ``VectorSet`` in
            gives a ``VectorSet`` out, re-framed to the specimen domain.
        """

        matrices = self.as_matrices()
        if isinstance(directions, VectorSet):
            if directions.reference_frame != self.crystal_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match OrientationSet.crystal_frame."
                )
            direction_array = directions.values
        else:
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
        if isinstance(directions, VectorSet):
            return VectorSet(
                values=mapped,
                reference_frame=self.specimen_frame,
                provenance=directions.provenance,
            )
        return mapped

    def map_sample_directions_to_crystal(
        self,
        directions: ArrayLike | VectorSet,
    ) -> np.ndarray | VectorSet:
        """Express specimen directions in the crystal frame, for the whole set.

        The inverse of :meth:`map_crystal_directions` and the vectorized action
        behind inverse pole figures: for each orientation, which crystal
        direction lies along the given specimen axis? Accepts one shared
        direction or one direction per orientation, with the same shape rules.
        """

        inverse_matrices = np.swapaxes(self.as_matrices(), -1, -2)
        if isinstance(directions, VectorSet):
            if directions.reference_frame != self.specimen_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match OrientationSet.specimen_frame."
                )
            direction_array = directions.values
        else:
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
        if isinstance(directions, VectorSet):
            return VectorSet(
                values=mapped,
                reference_frame=self.crystal_frame,
                provenance=directions.provenance,
            )
        return mapped

    def _orbit_operators(
        self, specimen_symmetry: SymmetrySpec | None
    ) -> tuple[np.ndarray, np.ndarray]:
        if (
            specimen_symmetry is not None
            and specimen_symmetry.reference_frame != self.specimen_frame
        ):
            raise ValueError(
                "specimen_symmetry.reference_frame must match OrientationSet.specimen_frame."
            )
        identity = np.eye(3, dtype=np.float64)[None, :, :]
        crystal_operators = self.symmetry.operators if self.symmetry is not None else identity
        specimen_operators = (
            specimen_symmetry.operators if specimen_symmetry is not None else identity
        )
        return crystal_operators, specimen_operators

    def canonicalized(self, specimen_symmetry: SymmetrySpec | None = None) -> OrientationSet:
        """Every member replaced by its sign-canonical symmetry representative.

        Makes quaternion arrays comparable across members and across runs. For
        the geometrically canonical choice — the representative inside the
        fundamental region — use :meth:`projected_to_fundamental_region`.
        """

        crystal_operators, specimen_operators = self._orbit_operators(specimen_symmetry)
        if len(self) == 0:
            return self
        representatives = _batched_fundamental_representatives(
            self.quaternions,
            crystal_operators=crystal_operators,
            specimen_operators=specimen_operators,
        )
        return OrientationSet(
            quaternions=representatives,
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
        """Every member reduced into the orientation fundamental region.

        Purpose
        -------
        Remove the arbitrary symmetry branch from a whole set at once, so that
        averaging, clustering, and plotting see one representative per physical
        orientation. Vectorized over the set.

        Parameters
        ----------
        specimen_symmetry : SymmetrySpec, optional
            Specimen symmetry to include in the orbit, in addition to the crystal
            symmetry the set already carries.

        Returns
        -------
        OrientationSet
            Same length and order, with reduced representatives.
        """

        if len(self) == 0:
            self._orbit_operators(specimen_symmetry)
            return self
        if reference_orientation is None:
            # Same criterion as canonicalisation: the lexicographically-largest
            # canonical quaternion in each orbit, vectorised over the whole set.
            crystal_operators, specimen_operators = self._orbit_operators(specimen_symmetry)
            projected = _batched_fundamental_representatives(
                self.quaternions,
                crystal_operators=crystal_operators,
                specimen_operators=specimen_operators,
            )
        else:
            # Reference-guided selection stays per orientation (rare path).
            projected = np.stack(
                [
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
                ],
                axis=0,
            )
        return OrientationSet(
            quaternions=projected,
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
        """``(n, 4)`` hashable keys identifying each member's symmetry class.

        Symmetrically equivalent orientations produce identical rows, so the keys
        can be used to group or deduplicate a set. They are keys, not
        coordinates: nearby orientations do not have nearby keys.
        """

        crystal_operators, specimen_operators = self._orbit_operators(specimen_symmetry)
        if len(self) == 0:
            keys = np.empty((0, 4), dtype=np.float64)
            keys.setflags(write=False)
            return keys
        representatives = _batched_fundamental_representatives(
            self.quaternions,
            crystal_operators=crystal_operators,
            specimen_operators=specimen_operators,
        )
        # The key is the negated, rounded canonical of each representative, i.e.
        # _exact_fundamental_region_key_from_quaternion applied per orientation.
        canonical = representatives.copy()
        negative = canonical[:, 0] < 0.0
        canonical[negative] = -canonical[negative]
        canonical = canonical / np.linalg.norm(canonical, axis=1, keepdims=True)
        keys = -np.round(canonical, 12)
        keys = np.ascontiguousarray(keys)
        keys.setflags(write=False)
        return keys

    def fundamental_region_keys(
        self,
        *,
        specimen_symmetry: SymmetrySpec | None = None,
    ) -> np.ndarray:
        """``(n, 4)`` fundamental-region keys; alias of
        :meth:`exact_fundamental_region_keys`.
        """

        return self.exact_fundamental_region_keys(specimen_symmetry=specimen_symmetry)

    def misorientation_angles_to(
        self,
        other: OrientationSet,
        *,
        symmetry_aware: bool = True,
    ) -> np.ndarray:
        """Pairwise misorientation angles between two orientation sets.

        Purpose
        -------
        The core measurement behind misorientation distributions, grain-boundary
        statistics, orientation clustering, and variant assignment.

        Parameters
        ----------
        other : OrientationSet
            Must share crystal frame, specimen frame, and (when both declare one)
            phase.
        symmetry_aware : bool
            When ``True`` (default) each entry is the disorientation angle — the
            minimum over both symmetry groups. When ``False`` the raw relative
            rotation angle is returned.

        Returns
        -------
        np.ndarray
            ``(len(self), len(other))`` angles in **radians**. Computed entirely
            in quaternion algebra, so the full ``(n, m, 3, 3)`` matrix stack is
            never materialized.
        """

        if self.crystal_frame != other.crystal_frame:
            raise ValueError("Misorientation requires the same crystal frame.")
        if self.specimen_frame != other.specimen_frame:
            raise ValueError("Misorientation requires the same specimen frame.")
        if self.phase is not None and other.phase is not None and self.phase != other.phase:
            raise ValueError("Misorientation requires matching phases when both are specified.")
        rows, columns = len(self), len(other)
        if rows == 0 or columns == 0:
            angles = np.zeros((rows, columns), dtype=np.float64)
            angles.setflags(write=False)
            return angles
        # Crystal-frame relative rotation inv(o_i) o_j = conj(q_i) * q_j for every
        # pair. Staying in quaternions avoids materialising an (n*m, 3, 3) array.
        left_quaternions = np.asarray(self.quaternions, dtype=np.float64)
        right_quaternions = np.asarray(other.quaternions, dtype=np.float64)
        conjugates = left_quaternions * np.array([1.0, -1.0, -1.0, -1.0])
        relative = quaternions_multiply(conjugates[:, None, :], right_quaternions[None, :, :])
        identity = np.eye(3, dtype=np.float64)[None, :, :]
        left_operators = (
            self.symmetry.operators if symmetry_aware and self.symmetry is not None else identity
        )
        right_operators = (
            other.symmetry.operators if symmetry_aware and other.symmetry is not None else identity
        )
        angles = _reduced_pair_disorientation_angles_from_quaternions(
            relative.reshape(rows * columns, 4), left_operators, right_operators
        ).reshape(rows, columns)
        angles = np.ascontiguousarray(angles)
        angles.setflags(write=False)
        return angles

    def _crystal_symmetry_quaternions(self) -> np.ndarray:
        if self.symmetry is None:
            return np.array([[1.0, 0.0, 0.0, 0.0]], dtype=np.float64)
        return matrices_to_quaternions(self.symmetry.operators)

    def mean_orientation(
        self,
        *,
        weights: ArrayLike | None = None,
        max_iterations: int = 10,
    ) -> Orientation:
        """The symmetry-aware mean orientation of the set.

        Purpose
        -------
        The representative orientation of a grain or of a cluster — the quantity
        behind grain mean orientations, GROD reference orientations, and
        orientation-relationship fitting.

        Method
        ------
        Quaternion eigenvector averaging: the mean is the dominant eigenvector
        of the weighted outer-product accumulator ``sum w_i q_i q_i^T``, which is
        the maximum-likelihood mean for small dispersions and avoids the bias of
        naive component averaging. Because symmetry-equivalent quaternions would
        otherwise cancel, the iteration re-selects, for every member, the
        symmetry branch and sign closest to the current estimate before
        accumulating. Iteration stops on convergence or at ``max_iterations``.

        Parameters
        ----------
        weights : ArrayLike, optional
            One non-negative weight per orientation, normalized internally.
            Uniform when omitted. Use this to weight by pattern quality or by
            pixel area.
        max_iterations : int
            Cap on the branch-reselection iterations (default 10).

        Returns
        -------
        Orientation
            The mean, canonicalized, carrying the set's frames, symmetry, and
            phase.

        Notes
        -----
        The mean is only meaningful for a unimodal set. A set spanning several
        distinct orientations returns a number with no physical meaning; check
        the dispersion with :meth:`spread_angles_deg` before trusting it.
        """

        if len(self) == 0:
            raise ValueError("mean_orientation requires at least one orientation.")
        count = len(self)
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
        symmetry_quaternions = self._crystal_symmetry_quaternions()
        candidates = quaternions_multiply(
            self.quaternions[:, None, :],
            symmetry_quaternions[None, :, :],
        )
        reference = self.quaternions[0]
        row_indices = np.arange(count)
        for _ in range(max_iterations):
            dots = candidates @ reference
            best = np.argmax(np.abs(dots), axis=1)
            selected = candidates[row_indices, best]
            signs = np.sign(dots[row_indices, best])
            signs[signs == 0.0] = 1.0
            selected = selected * signs[:, None]
            accumulator = np.einsum("n,ni,nj->ij", weight_values, selected, selected)
            eigenvalues, eigenvectors = np.linalg.eigh(accumulator)
            updated = eigenvectors[:, int(np.argmax(eigenvalues))]
            if float(updated @ reference) < 0.0:
                updated = -updated
            converged = bool(np.isclose(abs(float(updated @ reference)), 1.0, atol=1e-12))
            reference = updated
            if converged:
                break
        return Orientation(
            rotation=Rotation(quaternion=_canonicalize_quaternion(reference)),
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )

    def spread_angles_deg(
        self,
        *,
        reference: Orientation | None = None,
        symmetry_aware: bool = True,
    ) -> np.ndarray:
        """Angular deviation of each member from a reference orientation.

        Purpose
        -------
        The dispersion measure behind grain orientation spread (GOS) and grain
        reference orientation deviation (GROD): how far each measurement lies
        from the grain's representative orientation.

        Parameters
        ----------
        reference : Orientation, optional
            The reference. Defaults to :meth:`mean_orientation` of this set, so
            the result is a spread about the set's own mean. Must share the set's
            frames.
        symmetry_aware : bool
            Reduce each deviation by symmetry (default ``True``).

        Returns
        -------
        np.ndarray
            ``(n,)`` angles in **degrees**.
        """

        resolved = reference if reference is not None else self.mean_orientation()
        if resolved.crystal_frame != self.crystal_frame:
            raise ValueError("spread reference must share the OrientationSet crystal frame.")
        if resolved.specimen_frame != self.specimen_frame:
            raise ValueError("spread reference must share the OrientationSet specimen frame.")
        reference_set = OrientationSet(
            quaternions=resolved.rotation.quaternion[None, :],
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )
        angles_rad = self.misorientation_angles_to(
            reference_set,
            symmetry_aware=symmetry_aware,
        )[:, 0]
        angles_deg = np.rad2deg(angles_rad)
        angles_deg = np.ascontiguousarray(angles_deg)
        angles_deg.setflags(write=False)
        return angles_deg

    def subset(self, indices: ArrayLike) -> OrientationSet:
        """The orientations at the given indices, as a new set.

        Accepts integer indices or a boolean mask. Frames, symmetry, phase, and
        provenance are preserved, so a masked subset stays fully interpretable.
        """

        return OrientationSet(
            quaternions=self.quaternions[np.asarray(indices)],
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            symmetry=self.symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )
