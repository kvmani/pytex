from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    is_rotation_matrix,
    normalize_quaternion,
    normalize_quaternions,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.batches import EulerSet, RotationSet, VectorSet, normalize_euler_convention_name
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalDirection, CrystalPlane, Phase
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
    for matrix in matrix_array:
        if not is_rotation_matrix(matrix):
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

    def apply(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
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
        return self.rotation.as_matrix()

    def map_crystal_vector(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
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
    def from_plane_direction(
        cls,
        plane: CrystalPlane | list[CrystalPlane] | tuple[CrystalPlane, ...] | ArrayLike,
        direction: CrystalDirection
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
        if isinstance(plane, CrystalPlane):
            if not isinstance(direction, CrystalDirection):
                raise ValueError(
                    "A scalar CrystalPlane input requires a matching scalar CrystalDirection."
                )
            resolved_phase = _phase_from_plane_direction_objects(
                (plane,),
                (direction,),
                phase=phase,
            )
            crystal_normals = plane.normal[None, :]
            crystal_directions = direction.unit_vector[None, :]
        elif (
            isinstance(plane, (list, tuple))
            and isinstance(direction, (list, tuple))
            and all(isinstance(item, CrystalPlane) for item in plane)
            and all(isinstance(item, CrystalDirection) for item in direction)
        ):
            planes = tuple(plane)
            directions = tuple(direction)
            resolved_phase = _phase_from_plane_direction_objects(planes, directions, phase=phase)
            crystal_normals = normalize_vectors(np.vstack([item.normal for item in planes]))
            crystal_directions = normalize_vectors(
                np.vstack([item.unit_vector for item in directions])
            )
        else:
            if phase is None:
                raise ValueError(
                    "phase is required when constructing orientations from raw plane and "
                    "direction index arrays."
                )
            plane_indices = np.asarray(plane, dtype=np.float64)
            direction_indices = np.asarray(direction, dtype=np.float64)
            if plane_indices.shape == (3,):
                plane_indices = plane_indices[None, :]
            if direction_indices.shape == (3,):
                direction_indices = direction_indices[None, :]
            if plane_indices.ndim != 2 or plane_indices.shape[1] != 3:
                raise ValueError("plane index arrays must have shape (3,) or (n, 3).")
            if direction_indices.ndim != 2 or direction_indices.shape[1] != 3:
                raise ValueError("direction index arrays must have shape (3,) or (n, 3).")
            if plane_indices.shape[0] == 1 and direction_indices.shape[0] > 1:
                plane_indices = np.broadcast_to(plane_indices, direction_indices.shape)
            elif direction_indices.shape[0] == 1 and plane_indices.shape[0] > 1:
                direction_indices = np.broadcast_to(direction_indices, plane_indices.shape)
            elif plane_indices.shape[0] != direction_indices.shape[0]:
                raise ValueError(
                    "Raw plane and direction index arrays must broadcast to the same length."
                )
            reciprocal_basis = phase.lattice.reciprocal_basis().matrix
            direct_basis = phase.lattice.direct_basis().matrix
            crystal_normals = normalize_vectors(plane_indices @ reciprocal_basis.T)
            crystal_directions = normalize_vectors(direction_indices @ direct_basis.T)
            resolved_phase = phase
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
        matrices = np.stack(
            [quaternion_to_matrix(quaternion) for quaternion in self.quaternions],
            axis=0,
        )
        matrices = np.ascontiguousarray(matrices)
        matrices.setflags(write=False)
        return matrices

    def to_axes_angles(self) -> tuple[np.ndarray, np.ndarray]:
        return self.as_rotation_set().to_axes_angles()

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

    def as_euler_set(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> EulerSet:
        return EulerSet(
            angles=self.as_euler(convention=convention, degrees=degrees),
            convention=convention,
            degrees=degrees,
            provenance=self.provenance,
        )

    def as_rotation_set(self) -> RotationSet:
        return RotationSet(quaternions=self.quaternions, provenance=self.provenance)

    def map_crystal_directions(self, directions: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
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
