"""Coincidence-site-lattice (CSL) and twin-boundary classification for cubic crystals.

Boundaries are classified by comparing the measured crystal-to-crystal
misorientation against the ideal CSL misorientations (Sigma 1 through Sigma 29
for the cubic system) using the symmetry-reduced disorientation angle. A
boundary is assigned to a CSL Sigma when its deviation from that ideal is within
the Brandon criterion ``theta_max = theta0 * Sigma^(-1/2)`` (``theta0 = 15
degrees`` by convention).

Named twin laws are thin wrappers over the same registry; the coherent cubic
twin is the ``Sigma 3`` boundary ``60 degrees <111>``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

# Ideal cubic CSL misorientations as (Sigma, angle_deg, axis, variant). Values are
# the standard Grimmer/Randle axis-angle representatives for the disorientation
# (lowest-angle) form of each coincidence relation up to Sigma 29.
_CUBIC_CSL_TABLE: tuple[tuple[int, float, tuple[int, int, int], str], ...] = (
    (1, 0.0, (1, 0, 0), ""),
    (3, 60.0, (1, 1, 1), ""),
    (5, 36.86, (1, 0, 0), ""),
    (7, 38.21, (1, 1, 1), ""),
    (9, 38.94, (1, 1, 0), ""),
    (11, 50.47, (1, 1, 0), ""),
    (13, 22.62, (1, 0, 0), "a"),
    (13, 27.79, (1, 1, 1), "b"),
    (15, 48.19, (2, 1, 0), ""),
    (17, 28.07, (1, 0, 0), "a"),
    (17, 61.90, (2, 2, 1), "b"),
    (19, 26.53, (1, 1, 0), "a"),
    (19, 46.80, (1, 1, 1), "b"),
    (21, 21.79, (1, 1, 1), "a"),
    (21, 44.41, (2, 1, 1), "b"),
    (23, 40.45, (3, 1, 1), ""),
    (25, 16.26, (1, 0, 0), "a"),
    (25, 51.68, (3, 3, 1), "b"),
    (27, 31.59, (1, 1, 0), "a"),
    (27, 35.43, (2, 1, 0), "b"),
    (29, 43.60, (1, 0, 0), "a"),
    (29, 46.40, (2, 2, 1), "b"),
)


def brandon_tolerance_deg(sigma: int, *, theta0_deg: float = 15.0) -> float:
    """Maximum deviation for a Sigma boundary under the Brandon criterion."""

    if sigma < 1:
        raise ValueError("CSL Sigma must be a positive integer.")
    return theta0_deg / math.sqrt(sigma)


def _axis_angle_matrix(axis: ArrayLike, angle_rad: float) -> np.ndarray:
    unit = np.asarray(axis, dtype=np.float64)
    norm = float(np.linalg.norm(unit))
    if norm == 0.0:
        return np.eye(3, dtype=np.float64)
    unit = unit / norm
    x, y, z = unit
    cos = math.cos(angle_rad)
    sin = math.sin(angle_rad)
    one_minus = 1.0 - cos
    return np.array(
        [
            [cos + x * x * one_minus, x * y * one_minus - z * sin, x * z * one_minus + y * sin],
            [y * x * one_minus + z * sin, cos + y * y * one_minus, y * z * one_minus - x * sin],
            [z * x * one_minus - y * sin, z * y * one_minus + x * sin, cos + z * z * one_minus],
        ],
        dtype=np.float64,
    )


@dataclass(frozen=True, slots=True)
class CSLType:
    """An ideal cubic coincidence-site-lattice misorientation."""

    sigma: int
    angle_deg: float
    axis: tuple[int, int, int]
    variant: str = ""

    @property
    def name(self) -> str:
        return f"Sigma{self.sigma}{self.variant}"

    def matrix(self) -> np.ndarray:
        axis = np.asarray(self.axis, dtype=np.float64)
        return _axis_angle_matrix(axis, math.radians(self.angle_deg))


CUBIC_CSL_TYPES: tuple[CSLType, ...] = tuple(
    CSLType(sigma=sigma, angle_deg=angle, axis=axis, variant=variant)
    for sigma, angle, axis, variant in _CUBIC_CSL_TABLE
)


@dataclass(frozen=True, slots=True)
class TwinLaw:
    """A named twin boundary, defined by its ideal CSL misorientation."""

    name: str
    csl: CSLType

    def matrix(self) -> np.ndarray:
        return self.csl.matrix()


# The coherent cubic (fcc/bcc) twin is the Sigma3 60 deg <111> boundary.
CUBIC_TWIN_LAWS: tuple[TwinLaw, ...] = (
    TwinLaw(name="sigma3_twin", csl=CSLType(sigma=3, angle_deg=60.0, axis=(1, 1, 1))),
)


@dataclass(frozen=True, slots=True)
class CSLMatch:
    """Result of classifying a boundary against the CSL registry."""

    sigma: int
    deviation_deg: float
    csl: CSLType


def _reduced_deviation_deg(
    misorientation_matrices: np.ndarray,
    ideal_matrix: np.ndarray,
    operators: np.ndarray,
) -> np.ndarray:
    """Symmetry-reduced deviation angle between measured and ideal misorientations.

    Returns ``min over S_a, S_b of angle(S_a @ M @ S_b @ C^T)`` per matrix, i.e.
    the disorientation between each measured misorientation ``M`` and the ideal
    CSL misorientation ``C`` under the crystal point-group operators.
    """

    ops_times_ct = np.einsum("bij,kj->bik", operators, ideal_matrix, optimize=True)
    left = np.einsum("aij,sjk->saik", operators, misorientation_matrices, optimize=True)
    products = np.einsum("saij,bjk->sabik", left, ops_times_ct, optimize=True)
    traces = np.trace(products, axis1=3, axis2=4)
    cos_theta = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    angles = np.arccos(cos_theta)
    minimum = np.min(angles.reshape(misorientation_matrices.shape[0], -1), axis=1)
    return np.asarray(np.degrees(minimum), dtype=np.float64)


def classify_misorientations(
    misorientation_matrices: ArrayLike,
    *,
    operators: ArrayLike,
    csl_types: tuple[CSLType, ...] = CUBIC_CSL_TYPES,
    theta0_deg: float = 15.0,
    include_sigma1: bool = False,
) -> list[CSLMatch | None]:
    """Classify crystal-to-crystal misorientations against the CSL registry.

    ``misorientation_matrices`` are the relative rotations ``inv(o1) @ o2`` as
    3x3 matrices; ``operators`` are the crystal proper-rotation operators (24 for
    the cubic system). Each boundary is assigned the CSL type it deviates from
    least while staying within that type's Brandon tolerance; boundaries with no
    qualifying CSL type return ``None``. Sigma1 (low-angle) is excluded by
    default.
    """

    matrices = np.asarray(misorientation_matrices, dtype=np.float64)
    if matrices.ndim == 2:
        matrices = matrices[None, :, :]
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3):
        raise ValueError("misorientation_matrices must have shape (n, 3, 3) or (3, 3).")
    operator_array = np.asarray(operators, dtype=np.float64)
    if operator_array.ndim != 3 or operator_array.shape[1:] != (3, 3):
        raise ValueError("operators must have shape (k, 3, 3).")
    candidates = [
        csl for csl in csl_types if include_sigma1 or csl.sigma != 1
    ]
    best_sigma = np.full(matrices.shape[0], np.iinfo(np.int64).max, dtype=np.int64)
    best_deviation = np.full(matrices.shape[0], np.inf, dtype=np.float64)
    best_type: list[CSLType | None] = [None] * matrices.shape[0]
    for csl in candidates:
        tolerance = brandon_tolerance_deg(csl.sigma, theta0_deg=theta0_deg)
        deviation = _reduced_deviation_deg(matrices, csl.matrix(), operator_array)
        within = deviation <= tolerance
        improves = within & (
            (deviation < best_deviation)
            | ((deviation == best_deviation) & (csl.sigma < best_sigma))
        )
        for index in np.flatnonzero(improves):
            best_sigma[index] = csl.sigma
            best_deviation[index] = deviation[index]
            best_type[index] = csl
    return [
        None if csl is None else CSLMatch(sigma=csl.sigma, deviation_deg=float(dev), csl=csl)
        for csl, dev in zip(best_type, best_deviation, strict=True)
    ]


__all__ = [
    "CUBIC_CSL_TYPES",
    "CUBIC_TWIN_LAWS",
    "CSLMatch",
    "CSLType",
    "TwinLaw",
    "brandon_tolerance_deg",
    "classify_misorientations",
]
