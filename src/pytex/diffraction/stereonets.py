from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.sphere import (
    directions_to_spherical_angles,
    spherical_angles_to_directions,
)
from pytex.texture.projections import project_directions

__all__ = [
    "StereonetGrid",
    "directions_to_spherical_angles",
    "flatten_direction_grid",
    "generate_stereonet_grid",
    "project_great_circle_trace",
    "projection_boundary_radius",
    "sample_great_circle",
    "sample_small_circle",
    "spherical_angles_to_directions",
]


@dataclass(frozen=True, slots=True)
class StereonetGrid:
    """The projected line set of a stereonet graticule.

    Purpose
    -------
    The Wulff or Schmidt net drawn behind a pole figure, against which angles
    are read off by hand — the reference graticule of classical stereographic
    analysis.

    Attributes
    ----------
    major_lines, minor_lines : tuple of np.ndarray
        Projected polylines, heavy and light respectively.
    method : str
        ``"stereographic"`` (Wulff, angle-preserving) or ``"equal_area"``
        (Schmidt, area-preserving).
    boundary_radius : float
        Radius of the primitive circle in the chosen projection.
    """

    method: str
    major_lines: tuple[np.ndarray, ...]
    minor_lines: tuple[np.ndarray, ...]
    boundary_radius: float


def projection_boundary_radius(method: str) -> float:
    """Radius of the projection disc for a given projection method.

    ``sqrt(2)`` for equal-area, ``1`` for stereographic. Needed to draw the
    primitive circle and to set plot limits consistently with the data.
    """

    if method == "equal_area":
        return float(np.sqrt(2.0))
    if method == "stereographic":
        return 1.0
    raise ValueError("Projection method must be 'equal_area' or 'stereographic'.")


def sample_great_circle(
    normal: ArrayLike,
    *,
    samples: int = 361,
    half_circle: bool = True,
) -> np.ndarray:
    """Points along the great circle whose pole is the given normal.

    Purpose
    -------
    A great circle is the trace of a plane on the sphere — the curve a
    crystallographic plane draws on a stereographic projection, and the
    construction used to find zone axes as trace intersections.

    Parameters
    ----------
    normal : ArrayLike
        Pole of the great circle; normalized internally.
    samples : int
        Number of points; at least two.
    half_circle : bool
        Sample half the circle (default), which is what a one-hemisphere
        projection displays, or the full circle.

    Returns
    -------
    np.ndarray
        ``(samples, 3)`` unit vectors, read-only.
    """

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
    """Points along a small circle at a fixed polar angle.

    Small circles are loci of constant angle to the projection axis — the
    latitude lines of a Wulff net, and the construction used to find
    directions at a fixed angle to a known pole.

    Parameters
    ----------
    polar_deg : float
        Polar angle in ``[0, 90]``.
    samples : int
        Number of points around the full azimuthal turn.

    Returns
    -------
    np.ndarray
        ``(samples, 3)`` unit vectors, read-only.
    """

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
    """A great-circle trace projected onto the plotting plane.

    Composes :func:`sample_great_circle` with the projection, treating the
    directions as antipodal so the trace stays within one hemisphere.

    Parameters
    ----------
    normal : ArrayLike
        Pole of the great circle.
    method : str
        ``"stereographic"`` (default) or ``"equal_area"``.
    samples : int
        Number of points along the trace.

    Returns
    -------
    np.ndarray
        ``(samples, 2)`` plane coordinates, read-only.
    """

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
    """The meridian and small-circle line set of a stereonet.

    Purpose
    -------
    The Wulff (or Schmidt) net drawn behind a pole figure, against which
    angles are read off by hand — the reference graticule of classical
    stereographic analysis.

    Parameters
    ----------
    method : str
        ``"stereographic"`` (default; a Wulff net, angle-preserving) or
        ``"equal_area"`` (a Schmidt net, area-preserving).
    major_step_deg : float
        Angular spacing of the heavy lines.
    minor_step_deg : float, optional
        Spacing of the light lines; omitted means no minor graticule.
    samples : int
        Points sampled along each line.

    Returns
    -------
    StereonetGrid
        Projected major and minor line sets, ready to plot.
    """

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
    major_parallels = np.arange(
        float(major_step_deg),
        90.0,
        float(major_step_deg),
        dtype=np.float64,
    )
    for azimuth_deg in major_meridians:
        major_lines.append(_meridian_line(float(azimuth_deg)))
    for polar_deg in major_parallels:
        major_lines.append(_parallel_line(float(polar_deg)))

    if minor_step_deg is not None and not np.isclose(minor_step_deg, major_step_deg):
        minor_meridians = np.arange(0.0, 180.0, float(minor_step_deg), dtype=np.float64)
        minor_parallels = np.arange(
            float(minor_step_deg),
            90.0,
            float(minor_step_deg),
            dtype=np.float64,
        )
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
    """Flatten an ``(m, n, 3)`` direction grid to ``(m*n, 3)``.

    Bridges gridded direction fields to the flat ``(n, 3)`` form the
    projection and symmetry routines expect. Returns a read-only array.
    """

    grid = as_float_array(direction_grid, shape=(None, None, 3))
    flattened = np.ascontiguousarray(grid.reshape(-1, 3), dtype=np.float64)
    flattened.setflags(write=False)
    return flattened
