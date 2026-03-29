from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector, normalize_vectors
from pytex.texture.projections import project_directions


@dataclass(frozen=True, slots=True)
class StereonetGrid:
    method: str
    major_lines: tuple[np.ndarray, ...]
    minor_lines: tuple[np.ndarray, ...]
    boundary_radius: float


def projection_boundary_radius(method: str) -> float:
    if method == "equal_area":
        return float(np.sqrt(2.0))
    if method == "stereographic":
        return 1.0
    raise ValueError("Projection method must be 'equal_area' or 'stereographic'.")


def spherical_angles_to_directions(
    polar_deg: ArrayLike,
    azimuth_deg: ArrayLike,
) -> np.ndarray:
    polar, azimuth = np.broadcast_arrays(
        np.asarray(polar_deg, dtype=np.float64),
        np.asarray(azimuth_deg, dtype=np.float64),
    )
    polar_rad = np.deg2rad(polar)
    azimuth_rad = np.deg2rad(azimuth)
    directions = np.stack(
        [
            np.sin(polar_rad) * np.cos(azimuth_rad),
            np.sin(polar_rad) * np.sin(azimuth_rad),
            np.cos(polar_rad),
        ],
        axis=-1,
    )
    directions = np.ascontiguousarray(directions, dtype=np.float64)
    directions.setflags(write=False)
    return directions


def directions_to_spherical_angles(
    directions: ArrayLike,
    *,
    antipodal: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    vectors = np.array(normalize_vectors(directions), copy=True)
    if antipodal:
        mask = vectors[..., 2] < 0.0
        vectors[mask] *= -1.0
    polar = np.rad2deg(np.arccos(np.clip(vectors[..., 2], -1.0, 1.0)))
    azimuth = np.mod(np.rad2deg(np.arctan2(vectors[..., 1], vectors[..., 0])), 360.0)
    polar = np.ascontiguousarray(polar, dtype=np.float64)
    azimuth = np.ascontiguousarray(azimuth, dtype=np.float64)
    polar.setflags(write=False)
    azimuth.setflags(write=False)
    return polar, azimuth


def sample_great_circle(
    normal: ArrayLike,
    *,
    samples: int = 361,
    half_circle: bool = True,
) -> np.ndarray:
    if samples < 2:
        raise ValueError("Great-circle sampling requires at least two points.")
    pole = normalize_vector(normal)
    reference = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(pole, reference))), 1.0, atol=1e-8):
        reference = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    basis_u = normalize_vector(np.cross(pole, reference))
    basis_v = normalize_vector(np.cross(pole, basis_u))
    max_angle = np.pi if half_circle else 2.0 * np.pi
    angles = np.linspace(0.0, max_angle, samples, dtype=np.float64)
    directions = np.column_stack(
        [
            basis_u[0] * np.cos(angles) + basis_v[0] * np.sin(angles),
            basis_u[1] * np.cos(angles) + basis_v[1] * np.sin(angles),
            basis_u[2] * np.cos(angles) + basis_v[2] * np.sin(angles),
        ]
    )
    directions = np.ascontiguousarray(directions, dtype=np.float64)
    directions.setflags(write=False)
    return directions


def sample_small_circle(
    polar_deg: float,
    *,
    samples: int = 361,
) -> np.ndarray:
    if not 0.0 <= float(polar_deg) <= 90.0:
        raise ValueError("Small-circle polar angles must lie in the interval [0, 90].")
    azimuth_deg = np.linspace(0.0, 360.0, samples, dtype=np.float64)
    directions = spherical_angles_to_directions(
        np.full(samples, float(polar_deg), dtype=np.float64),
        azimuth_deg,
    ).reshape(-1, 3)
    directions = np.ascontiguousarray(directions, dtype=np.float64)
    directions.setflags(write=False)
    return directions


def project_great_circle_trace(
    normal: ArrayLike,
    *,
    method: str = "stereographic",
    samples: int = 361,
) -> np.ndarray:
    projected = project_directions(
        sample_great_circle(normal, samples=samples, half_circle=True),
        method=method,
        antipodal=True,
    )
    projected = np.ascontiguousarray(projected, dtype=np.float64)
    projected.setflags(write=False)
    return projected


def generate_stereonet_grid(
    *,
    method: str = "stereographic",
    major_step_deg: float = 10.0,
    minor_step_deg: float | None = None,
    samples: int = 361,
) -> StereonetGrid:
    if major_step_deg <= 0.0:
        raise ValueError("major_step_deg must be strictly positive.")
    if minor_step_deg is not None and minor_step_deg <= 0.0:
        raise ValueError("minor_step_deg must be strictly positive when provided.")
    major_lines: list[np.ndarray] = []
    minor_lines: list[np.ndarray] = []

    def _meridian_line(azimuth_deg: float) -> np.ndarray:
        direction = spherical_angles_to_directions(90.0, azimuth_deg).reshape(3)
        normal = np.cross(direction, np.array([0.0, 0.0, 1.0], dtype=np.float64))
        return project_great_circle_trace(normal, method=method, samples=samples)

    def _parallel_line(polar_deg: float) -> np.ndarray:
        circle = project_directions(
            sample_small_circle(float(polar_deg), samples=samples),
            method=method,
            antipodal=False,
        )
        circle = np.ascontiguousarray(circle, dtype=np.float64)
        circle.setflags(write=False)
        return circle

    major_meridians = np.arange(0.0, 180.0, float(major_step_deg), dtype=np.float64)
    major_parallels = np.arange(float(major_step_deg), 90.0, float(major_step_deg), dtype=np.float64)
    for azimuth_deg in major_meridians:
        major_lines.append(_meridian_line(float(azimuth_deg)))
    for polar_deg in major_parallels:
        major_lines.append(_parallel_line(float(polar_deg)))

    if minor_step_deg is not None and not np.isclose(minor_step_deg, major_step_deg):
        minor_meridians = np.arange(0.0, 180.0, float(minor_step_deg), dtype=np.float64)
        minor_parallels = np.arange(float(minor_step_deg), 90.0, float(minor_step_deg), dtype=np.float64)
        for azimuth_deg in minor_meridians:
            if np.isclose(np.mod(float(azimuth_deg), float(major_step_deg)), 0.0, atol=1e-8):
                continue
            minor_lines.append(_meridian_line(float(azimuth_deg)))
        for polar_deg in minor_parallels:
            if np.isclose(np.mod(float(polar_deg), float(major_step_deg)), 0.0, atol=1e-8):
                continue
            minor_lines.append(_parallel_line(float(polar_deg)))

    return StereonetGrid(
        method=method,
        major_lines=tuple(major_lines),
        minor_lines=tuple(minor_lines),
        boundary_radius=projection_boundary_radius(method),
    )


def flatten_direction_grid(direction_grid: ArrayLike) -> np.ndarray:
    grid = as_float_array(direction_grid, shape=(None, None, 3))
    flattened = np.ascontiguousarray(grid.reshape(-1, 3), dtype=np.float64)
    flattened.setflags(write=False)
    return flattened
