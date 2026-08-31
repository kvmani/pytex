"""Angle recovery that keeps its digits near zero.

``arccos`` is the obvious way to turn a dot product into an angle and the wrong
one near the endpoint. Its derivative is unbounded there, so an argument correct
to machine precision -- an error of order ``1e-16`` -- comes back as an angle
wrong by order ``sqrt(1e-16) = 1e-8`` radians, about ``1e-6`` degrees. Half the
significant digits are gone, and they are exactly the digits that matter: the
answer to "how far is this from exact?" is the whole question that alignment,
disorientation and residual scatter are asking.

That is not hypothetical. It is why a rotation reconstructed from four exactly
consistent zone axes reported a scatter of ``1.5e-06`` degrees rather than zero,
and why it did so on Linux while passing on Windows: the noise the ``arccos``
amplifies is a difference in BLAS summation order, so the *platform* decided
whether an exact identity looked exact.

The cure is never to route the answer through the cosine. A cosine near ``1``
has already lost the information -- recovering the angle from it more cleverly
cannot put it back -- so each function here takes the geometry itself and finds
the small quantity directly:

- between unit vectors, Kahan's ``2 atan2(||a - b||, ||a + b||)``;
- from a rotation matrix, ``atan2`` of the skew part against the trace, where
  the skew entries are of order the angle rather than its square;
- from a quaternion, ``2 atan2(||q_vec||, |q_w|)``, where the vector part is of
  order the angle.

All three are exact to the last digit at zero, which is where the interesting
cases live. Everything here is private: these are numerical utilities, not
domain surfaces.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def angle_between_unit_vectors_rad(first: ArrayLike, second: ArrayLike) -> np.ndarray:
    r"""Angle in radians between rows of two arrays of **unit** vectors.

    Kahan's formulation, :math:`2\arctan_2(\lVert a-b\rVert, \lVert a+b\rVert)`,
    which is well conditioned at both ends: near zero the numerator is the small
    quantity and carries full relative precision, and near :math:`\pi` the
    denominator is.

    Parameters
    ----------
    first, second : ArrayLike
        ``(..., n)`` arrays of unit vectors, broadcast against each other. Rows
        that are not unit vectors give a meaningless answer rather than an
        error; normalize first.

    Returns
    -------
    np.ndarray
        Angles in radians, in ``[0, pi]``, with the broadcast leading shape.
    """

    left = np.asarray(first, dtype=np.float64)
    right = np.asarray(second, dtype=np.float64)
    difference = np.linalg.norm(left - right, axis=-1)
    total = np.linalg.norm(left + right, axis=-1)
    return np.asarray(2.0 * np.arctan2(difference, total), dtype=np.float64)


def acute_angle_between_unit_vectors_rad(first: ArrayLike, second: ArrayLike) -> np.ndarray:
    """As `angle_between_unit_vectors_rad`, ignoring the sense of each vector.

    The angle between the *lines* the vectors lie along, in ``[0, pi/2]``: the
    answer for an axis, a plane normal, or any other object whose sign carries
    no meaning. The second vector is flipped where the pair is obtuse, so the
    conditioning of the acute branch is what applies.
    """

    left = np.asarray(first, dtype=np.float64)
    right = np.asarray(second, dtype=np.float64)
    obtuse = np.sum(left * right, axis=-1) < 0.0
    aligned = np.where(obtuse[..., None], -right, right)
    return angle_between_unit_vectors_rad(left, aligned)


def rotation_angle_from_matrix_rad(matrices: ArrayLike) -> np.ndarray:
    r"""Rotation angle in radians of proper rotation matrices.

    The textbook :math:`\arccos((\operatorname{tr}R - 1)/2)` is the form this
    module exists to replace: for a near-identity rotation the trace is
    :math:`3 - \theta^2`, so the angle is recovered from a quantity that
    encodes it *squared* and half the digits are gone before the ``arccos``
    is even reached.

    The skew part does not have that problem. :math:`R - R^{T}` has entries of
    order :math:`\theta`, so

    .. math:: \theta = \arctan_2\!\left(
        \tfrac{1}{2}\lVert \operatorname{vec}(R - R^{T}) \rVert,\;
        \tfrac{1}{2}(\operatorname{tr}R - 1)\right)

    is accurate to the last digit at zero and remains correct through
    :math:`\pi`, where the trace takes over as the informative term.

    Parameters
    ----------
    matrices : ArrayLike
        ``(..., 3, 3)`` proper rotation matrices.

    Returns
    -------
    np.ndarray
        Angles in radians, in ``[0, pi]``.
    """

    array = np.asarray(matrices, dtype=np.float64)
    if array.shape[-2:] != (3, 3):
        raise ValueError("rotation_angle_from_matrix_rad expects (..., 3, 3) matrices.")
    axis_vector = 0.5 * np.stack(
        (
            array[..., 2, 1] - array[..., 1, 2],
            array[..., 0, 2] - array[..., 2, 0],
            array[..., 1, 0] - array[..., 0, 1],
        ),
        axis=-1,
    )
    sines = np.linalg.norm(axis_vector, axis=-1)
    cosines = 0.5 * (np.trace(array, axis1=-2, axis2=-1) - 1.0)
    return np.asarray(np.arctan2(sines, cosines), dtype=np.float64)


def rotation_angle_from_quaternion_rad(quaternions: ArrayLike) -> np.ndarray:
    """Rotation angle in radians of unit quaternions in ``(w, x, y, z)`` order.

    ``2 atan2(||q_vec||, |q_w|)``. The vector part is of order the angle, so an
    exactly-identity quaternion gives exactly zero rather than the ``3e-08``
    radians ``2 arccos(w)`` floors out at. The scalar part is taken in absolute
    value, which puts the answer on the short branch (``[0, pi]``) whichever of
    the two antipodal representations was handed in.

    Parameters
    ----------
    quaternions : ArrayLike
        ``(..., 4)`` unit quaternions, scalar part first.
    """

    array = np.asarray(quaternions, dtype=np.float64)
    if array.shape[-1] != 4:
        raise ValueError("rotation_angle_from_quaternion_rad expects (..., 4) quaternions.")
    scalar = np.abs(array[..., 0])
    vector = np.linalg.norm(array[..., 1:], axis=-1)
    return np.asarray(2.0 * np.arctan2(vector, scalar), dtype=np.float64)


def rotation_angle_between_quaternions_rad(first: ArrayLike, second: ArrayLike) -> np.ndarray:
    """Rotation angle in radians between two arrays of unit quaternions.

    The rotation angle is twice the four-dimensional angle between the
    quaternions, so Kahan's formula applies once the antipodal ambiguity is
    resolved: flip whichever pairs point apart, then ``4 atan2(||a - b||,
    ||a + b||)``. Identical quaternions give exactly zero.

    Parameters
    ----------
    first, second : ArrayLike
        ``(..., 4)`` unit quaternions, scalar part first, broadcast against each
        other.
    """

    left = np.asarray(first, dtype=np.float64)
    right = np.asarray(second, dtype=np.float64)
    opposed = np.sum(left * right, axis=-1) < 0.0
    aligned = np.where(opposed[..., None], -right, right)
    difference = np.linalg.norm(left - aligned, axis=-1)
    total = np.linalg.norm(left + aligned, axis=-1)
    return np.asarray(4.0 * np.arctan2(difference, total), dtype=np.float64)
