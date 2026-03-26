from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _validate_shape(array: np.ndarray, shape: tuple[int | None, ...]) -> None:
    if array.ndim != len(shape):
        raise ValueError(f"Expected {len(shape)} dimensions, received {array.ndim}.")
    for actual, expected in zip(array.shape, shape, strict=True):
        if expected is not None and actual != expected:
            raise ValueError(f"Expected shape {shape}, received {array.shape}.")


def freeze_array(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


def as_float_array(
    values: ArrayLike,
    *,
    shape: tuple[int | None, ...] | None = None,
) -> FloatArray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    if shape is not None:
        _validate_shape(array, shape)
    return freeze_array(array)


def as_int_array(
    values: ArrayLike,
    *,
    shape: tuple[int | None, ...] | None = None,
) -> IntArray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    if shape is not None:
        _validate_shape(array, shape)
    return freeze_array(array)


def as_str_tuple(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(str(value) for value in values)


def normalize_vector(vector: ArrayLike) -> FloatArray:
    array = as_float_array(vector, shape=(3,))
    norm = float(np.linalg.norm(array))
    if np.isclose(norm, 0.0):
        raise ValueError("Cannot normalize a zero vector.")
    return as_float_array(array / norm, shape=(3,))


def normalize_vectors(vectors: ArrayLike) -> FloatArray:
    array = as_float_array(vectors, shape=(None, 3))
    norms = np.linalg.norm(array, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError("Cannot normalize a vector set containing a zero vector.")
    return as_float_array(array / norms[:, None], shape=(None, 3))


def normalize_quaternion(quaternion: ArrayLike) -> FloatArray:
    array = as_float_array(quaternion, shape=(4,))
    norm = float(np.linalg.norm(array))
    if np.isclose(norm, 0.0):
        raise ValueError("Cannot normalize a zero quaternion.")
    return as_float_array(array / norm, shape=(4,))


def normalize_quaternions(quaternions: ArrayLike) -> FloatArray:
    array = as_float_array(quaternions, shape=(None, 4))
    norms = np.linalg.norm(array, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError("Cannot normalize a quaternion set containing a zero quaternion.")
    return as_float_array(array / norms[:, None], shape=(None, 4))


def is_rotation_matrix(matrix: ArrayLike, *, atol: float = 1e-8) -> bool:
    array = np.asarray(matrix, dtype=np.float64)
    if array.shape != (3, 3):
        return False
    identity = np.eye(3)
    return np.allclose(array.T @ array, identity, atol=atol) and np.isclose(
        np.linalg.det(array), 1.0, atol=atol
    )
