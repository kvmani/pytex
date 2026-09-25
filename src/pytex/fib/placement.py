"""Where inside the grain the lamella goes: the footprint fit of section 7.3.

Purpose
-------
Once the orientation is fixed the azimuth ``theta_S`` is fixed too, so the only
freedom left in placing the ``L x W`` rectangle is translation. This module finds
the translation that keeps the rectangle furthest from the grain boundary, and
reports that clearance in micrometres --- or, if the rectangle does not fit at
all, says so and reports the longest lamella that does.

Algorithm
---------
1. The grain mask is resampled onto a grid aligned with the lamella: columns
   along the long axis ``t_L``, rows along the normal ``n_L``, at the map's
   step. Sampling is nearest-neighbour, so a pixel is inside exactly when the
   EBSD point it lands on belongs to the grain; nothing is interpolated.
2. A summed-area table of the resampled mask gives the number of grain pixels
   under any axis-aligned window in constant time, so "does an ``a x b``
   window fit anywhere" is one vectorized comparison over the whole grid.
3. The clearance is the largest ``k`` for which the rectangle grown by ``k``
   cells on every side still fits; it is found by bisection, since fitting is
   monotone in ``k``.
4. Among the positions that achieve it, the one deepest inside the grain (by
   the Euclidean distance transform) is taken, which centres the rectangle in
   a symmetric grain.

Everything is NumPy and :mod:`scipy.ndimage`: no dependency beyond the core
stack (the air-gapped deployment constraint of the foundation document).

Units and frames
----------------
The mask lives on the EBSD scan grid; lengths are scan micrometres. The
lamella axes are sample-frame directions and are brought onto the scan through
:class:`pytex.fib.frames.SurfaceGeometry`, so a scan whose rows run against the
specimen ``y`` axis is handled rather than mirrored.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.fib.frames import ImageRegistration, SurfaceGeometry

__all__ = ["PlacementResult", "fit_footprint", "rectangle_corners"]


@dataclass(frozen=True, slots=True)
class PlacementResult:
    """The best placement of the lamella rectangle inside one grain.

    Attributes
    ----------
    fits : bool
        Whether the full ``L x W`` rectangle fits inside the grain.
    margin_um : float
        Clearance between the rectangle and the grain boundary, on every side;
        ``0`` when it does not fit. Resolved to the map step.
    center_scan_um : (x, y)
        Rectangle centre in scan micrometres.
    center_image : (x, y)
        The same point in SEM image units.
    corners_scan_um : (4, 2) array
        Corners in scan micrometres, counter-clockwise from the
        ``-t_L, -n_L`` corner.
    corners_image : (4, 2) array
        The corners in image units.
    length_um, width_um : float
        The rectangle placed (the requested ``L`` and ``W``).
    largest_length_um : float
        When the rectangle does not fit: the longest lamella of width ``W``
        that does (``0`` if not even the width fits). Equal to ``length_um``
        when it fits.
    resolution_um : float
        Cell size of the lamella-aligned grid.
    edge_distance_um : float
        Distance from the centre to the nearest edge of the map.
    theta_sample_deg : float
        Azimuth of ``n_L`` used.
    """

    fits: bool
    margin_um: float
    center_scan_um: tuple[float, float]
    center_image: tuple[float, float]
    corners_scan_um: np.ndarray
    corners_image: np.ndarray
    length_um: float
    width_um: float
    largest_length_um: float
    resolution_um: float
    edge_distance_um: float
    theta_sample_deg: float

    def __post_init__(self) -> None:
        for name in ("corners_scan_um", "corners_image"):
            array = np.ascontiguousarray(getattr(self, name), dtype=np.float64)
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    def describe(self) -> str:
        """The placement in one or two sentences."""

        cx, cy = self.center_scan_um
        if self.fits:
            return (
                f"Placement: the {self.length_um:g} x {self.width_um:g} um rectangle fits with "
                f"{self.margin_um:.2f} um clearance to the grain boundary on every side, centred "
                f"at scan ({cx:.2f}, {cy:.2f}) um, {self.edge_distance_um:.1f} um "
                "from the map edge "
                f"(resolved to {self.resolution_um:.2f} um)."
            )
        if self.largest_length_um > 0.0:
            return (
                f"Placement: the {self.length_um:g} um lamella does NOT fit inside the grain at "
                f"this azimuth; the longest that fits is {self.largest_length_um:.2f} um, centred "
                f"at scan ({cx:.2f}, {cy:.2f}) um."
            )
        return (
            f"Placement: not even the {self.width_um:g} um width fits inside the grain at this "
            "azimuth; choose another grain."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "fits": bool(self.fits),
            "margin_um": float(self.margin_um),
            "center_scan_um": [float(value) for value in self.center_scan_um],
            "center_image": [float(value) for value in self.center_image],
            "corners_scan_um": self.corners_scan_um.tolist(),
            "corners_image": self.corners_image.tolist(),
            "length_um": float(self.length_um),
            "width_um": float(self.width_um),
            "largest_length_um": float(self.largest_length_um),
            "resolution_um": float(self.resolution_um),
            "edge_distance_um": float(self.edge_distance_um),
            "theta_sample_deg": float(self.theta_sample_deg),
        }


def _lamella_axes_scan(
    theta_sample_deg: float, surface: SurfaceGeometry
) -> tuple[np.ndarray, np.ndarray]:
    """Unit ``t_L`` and ``n_L`` expressed in scan coordinates."""

    theta = math.radians(theta_sample_deg)
    normal = surface.sample_to_scan_xy(np.array([math.cos(theta), math.sin(theta)]))
    long_axis = surface.sample_to_scan_xy(np.array([-math.sin(theta), math.cos(theta)]))
    return long_axis, normal


def rectangle_corners(
    center_scan_um: ArrayLike,
    length_um: float,
    width_um: float,
    theta_sample_deg: float,
    surface: SurfaceGeometry | None = None,
) -> np.ndarray:
    """The four corners of the lamella rectangle, in scan micrometres."""

    surface = surface or SurfaceGeometry()
    long_axis, normal = _lamella_axes_scan(theta_sample_deg, surface)
    centre = np.asarray(center_scan_um, dtype=np.float64)
    half_l, half_w = 0.5 * length_um, 0.5 * width_um
    offsets = [(-half_l, -half_w), (half_l, -half_w), (half_l, half_w), (-half_l, half_w)]
    return np.array([centre + a * long_axis + b * normal for a, b in offsets])


def _window_sums(table: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """Sum of every ``rows x cols`` window, from a zero-padded summed-area table."""

    sums = table[rows:, cols:] - table[:-rows, cols:] - table[rows:, :-cols] + table[:-rows, :-cols]
    return np.asarray(sums)


def _fits(table: np.ndarray, rows: int, cols: int) -> np.ndarray | None:
    """Boolean map of window positions fully inside the mask, or None if none fit."""

    height, width = table.shape[0] - 1, table.shape[1] - 1
    if rows > height or cols > width or rows < 1 or cols < 1:
        return None
    sums = _window_sums(table, rows, cols)
    inside = sums == rows * cols
    return inside if bool(inside.any()) else None


def fit_footprint(
    grain_mask: ArrayLike,
    *,
    step_um: tuple[float, float],
    origin_um: tuple[float, float] = (0.0, 0.0),
    theta_sample_deg: float,
    length_um: float,
    width_um: float,
    surface: SurfaceGeometry | None = None,
    registration: ImageRegistration | None = None,
) -> PlacementResult:
    """Place an ``L x W`` rectangle at a fixed azimuth as far inside a grain as possible.

    Parameters
    ----------
    grain_mask : (rows, cols) bool array_like
        ``True`` on the grain's EBSD points. Row index increases with scan
        ``y``, column index with scan ``x``.
    step_um : (dx, dy)
        Scan step along ``x`` (columns) and ``y`` (rows).
    origin_um : (x0, y0), default (0, 0)
        Scan coordinates of pixel ``[0, 0]``.
    theta_sample_deg : float
        Azimuth of the lamella normal in the sample frame.
    length_um, width_um : float
        ``L`` along ``t_L`` and ``W`` along ``n_L``.
    surface : SurfaceGeometry, optional
    registration : ImageRegistration, optional
        For the image-frame corners.

    Returns
    -------
    PlacementResult

    Examples
    --------
    A 40 x 20 um grain aligned with the scan, a 15 x 2 um lamella along ``x``
    (normal along ``y``, ``theta_S = 90``): the clearance is limited by the short
    side, ``(20 - 2) / 2 = 9`` um.
    """

    from scipy import ndimage

    surface = surface or SurfaceGeometry()
    registration = registration or ImageRegistration.identity()
    mask = np.asarray(grain_mask, dtype=bool)
    if mask.ndim != 2 or not mask.any():
        raise ValueError("fit_footprint needs a non-empty 2D grain mask.")
    dx, dy = (float(value) for value in step_um)
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("step_um must be positive.")
    if length_um <= 0.0 or width_um <= 0.0:
        raise ValueError("The rectangle must have positive length and width.")
    x0, y0 = (float(value) for value in origin_um)
    rows_total, cols_total = mask.shape
    resolution = min(dx, dy)
    long_axis, normal = _lamella_axes_scan(theta_sample_deg, surface)

    # Points of the grain, in scan micrometres, to bound the rotated grid.
    iy, ix = np.nonzero(mask)
    points = np.column_stack([x0 + ix * dx, y0 + iy * dy])
    along = points @ long_axis
    across = points @ normal
    pad = 2.0 * resolution
    u = np.arange(along.min() - pad, along.max() + pad + resolution, resolution)
    v = np.arange(across.min() - pad, across.max() + pad + resolution, resolution)
    uu, vv = np.meshgrid(u, v)  # rows along v (n_L), columns along u (t_L)
    world = uu[..., None] * long_axis + vv[..., None] * normal
    col = np.rint((world[..., 0] - x0) / dx).astype(np.int64)
    row = np.rint((world[..., 1] - y0) / dy).astype(np.int64)
    valid = (col >= 0) & (col < cols_total) & (row >= 0) & (row < rows_total)
    resampled = np.zeros(uu.shape, dtype=bool)
    resampled[valid] = mask[row[valid], col[valid]]

    table = np.zeros((resampled.shape[0] + 1, resampled.shape[1] + 1), dtype=np.int64)
    table[1:, 1:] = np.cumsum(np.cumsum(resampled.astype(np.int64), axis=0), axis=1)
    distance = ndimage.distance_transform_edt(resampled)

    length_cells = max(1, math.ceil(length_um / resolution - 1e-9))
    width_cells = max(1, math.ceil(width_um / resolution - 1e-9))

    def best_position(rows: int, cols: int, inside: np.ndarray) -> tuple[float, float]:
        # Centre of each window, and the depth of that centre inside the grain.
        centre_rows = np.arange(inside.shape[0]) + (rows - 1) / 2.0
        centre_cols = np.arange(inside.shape[1]) + (cols - 1) / 2.0
        depth = ndimage.map_coordinates(
            distance,
            np.meshgrid(centre_rows, centre_cols, indexing="ij"),
            order=1,
            mode="nearest",
        )
        score = np.where(inside, depth, -np.inf)
        # Many positions can tie on depth (a long grain constrains only one
        # direction); among them take the one nearest their own centroid, so
        # the rectangle sits in the middle of the admissible run.
        deepest = score >= float(score.max()) - 1e-6
        rr, cc = np.nonzero(deepest)
        nearest = int(np.argmin((rr - rr.mean()) ** 2 + (cc - cc.mean()) ** 2))
        r, c = int(rr[nearest]), int(cc[nearest])
        return float(v[0] + centre_rows[r] * resolution), float(u[0] + centre_cols[c] * resolution)

    base = _fits(table, width_cells, length_cells)
    if base is not None:
        low, high = 0, (min(resampled.shape) // 2) + 1
        while low < high:
            middle = (low + high + 1) // 2
            if _fits(table, width_cells + 2 * middle, length_cells + 2 * middle) is None:
                high = middle - 1
            else:
                low = middle
        margin_cells = low
        inside = _fits(table, width_cells + 2 * margin_cells, length_cells + 2 * margin_cells)
        assert inside is not None
        v_centre, u_centre = best_position(
            width_cells + 2 * margin_cells, length_cells + 2 * margin_cells, inside
        )
        fits = True
        largest = float(length_um)
        margin_um = margin_cells * resolution
    else:
        fits = False
        margin_um = 0.0
        width_fit = _fits(table, width_cells, 1)
        if width_fit is None:
            largest = 0.0
            v_centre = float(np.mean(across))
            u_centre = float(np.mean(along))
        else:
            low, high = 1, length_cells - 1
            while low < high:
                middle = (low + high + 1) // 2
                if _fits(table, width_cells, middle) is None:
                    high = middle - 1
                else:
                    low = middle
            inside = _fits(table, width_cells, low)
            assert inside is not None
            v_centre, u_centre = best_position(width_cells, low, inside)
            largest = low * resolution

    centre = u_centre * long_axis + v_centre * normal
    placed_length = float(length_um) if fits else float(largest)
    corners = rectangle_corners(
        centre,
        placed_length if placed_length > 0 else length_um,
        width_um,
        theta_sample_deg,
        surface,
    )
    x_max = x0 + (cols_total - 1) * dx
    y_max = y0 + (rows_total - 1) * dy
    edge = float(min(centre[0] - x0, x_max - centre[0], centre[1] - y0, y_max - centre[1]))
    return PlacementResult(
        fits=fits,
        margin_um=float(margin_um),
        center_scan_um=(float(centre[0]), float(centre[1])),
        center_image=tuple(float(value) for value in registration.to_image(centre)),  # type: ignore[arg-type]
        corners_scan_um=corners,
        corners_image=registration.to_image(corners),
        length_um=float(length_um),
        width_um=float(width_um),
        largest_length_um=float(largest),
        resolution_um=float(resolution),
        edge_distance_um=max(0.0, edge),
        theta_sample_deg=float(theta_sample_deg),
    )
