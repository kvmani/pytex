from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, freeze_array, normalize_vectors
from pytex.core.batches import VectorSet
from pytex.core.conventions import FrameDomain
from pytex.core.orientation import Orientation, OrientationSet, specimen_direction_vector
from pytex.core.sphere import spherical_angles_to_directions
from pytex.core.symmetry import FundamentalSector, SymmetrySpec
from pytex.texture.projections import project_directions

_DEFAULT_ANCHOR_DIRECTIONS = np.array(
    [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ],
    dtype=np.float64,
)


def _normalize_rgb_vertices(vertices: ArrayLike) -> np.ndarray:
    array = as_float_array(vertices, shape=(3, 3))
    if np.any(~np.isfinite(array)) or np.any((array < 0.0) | (array > 1.0)):
        raise ValueError("IPF color vertices must be finite RGB triplets in [0, 1].")
    return array


@dataclass(frozen=True, slots=True)
class IPFColorKey:
    """The colour mapping from crystal directions to inverse-pole-figure colours.

    Purpose
    -------
    Defines the standard IPF colouring: each direction is folded into the
    symmetry fundamental sector and coloured by its barycentric position, so
    symmetrically equivalent directions receive the same colour and the
    sector corners anchor the primary colours. It also produces the legend
    mesh that makes an IPF map interpretable.

    Attributes
    ----------
    symmetry : SymmetrySpec
        Determines the fundamental sector, and hence the colour scheme. Two
        maps of different symmetries are not colour-comparable.
    Remaining attributes carry the projection and colour-assignment options.
    """

    crystal_symmetry: SymmetrySpec
    specimen_direction: str | np.ndarray
    antipodal: bool = True
    saturation_gamma: float = 0.5
    color_vertices: np.ndarray = field(
        default_factory=lambda: np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
    )

    def __post_init__(self) -> None:
        if self.crystal_symmetry.reference_frame is None:
            raise ValueError("IPFColorKey.crystal_symmetry must declare a crystal reference frame.")
        if self.crystal_symmetry.reference_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError(
                "IPFColorKey.crystal_symmetry must use a crystal-domain reference frame."
            )
        object.__setattr__(
            self,
            "specimen_direction",
            specimen_direction_vector(self.specimen_direction),
        )
        object.__setattr__(self, "color_vertices", _normalize_rgb_vertices(self.color_vertices))
        if not np.isfinite(self.saturation_gamma) or self.saturation_gamma <= 0.0:
            raise ValueError("IPFColorKey.saturation_gamma must be finite and strictly positive.")

    @property
    def sector(self) -> FundamentalSector:
        """The symmetry fundamental sector this colour key maps over.
        """

        return self.crystal_symmetry.fundamental_sector(antipodal=self.antipodal)

    @property
    def sector_vertices(self) -> np.ndarray:
        """Corner directions of the fundamental sector, as unit vectors.

        For cubic symmetry these are the ``[001]``, ``[101]``, ``[111]`` corners
        that anchor the colour key.
        """

        vertices = self.sector.vertices
        if vertices.shape == (3, 3):
            return vertices
        # Low-symmetry sectors (hemisphere or wedge) have fewer than three
        # corners; anchor the barycentric color basis on the reference octant.
        return as_float_array(_DEFAULT_ANCHOR_DIRECTIONS, shape=(3, 3))

    def colors_from_crystal_directions(self, crystal_directions: ArrayLike) -> np.ndarray:
        """RGB colours for crystal directions, by position in the sector.

        Purpose
        -------
        The core of IPF colouring: each direction is folded into the fundamental
        sector and assigned a colour from its barycentric position, so
        symmetrically equivalent directions receive the same colour and the
        sector corners anchor the primary colours.

        Parameters
        ----------
        directions : ArrayLike
            ``(n, 3)`` crystal directions.

        Returns
        -------
        np.ndarray
            ``(n, 3)`` RGB values in ``[0, 1]``.
        """

        directions = normalize_vectors(crystal_directions)
        reduced = self.crystal_symmetry.reduce_vectors_to_fundamental_sector(
            directions,
            antipodal=self.antipodal,
        )
        reduced_array = reduced.values if isinstance(reduced, VectorSet) else reduced
        sector_basis = np.asarray(self.sector_vertices, dtype=np.float64).T
        barycentric = np.linalg.solve(sector_basis, reduced_array.T).T
        barycentric = np.clip(barycentric, 0.0, None)
        sums = barycentric.sum(axis=1, keepdims=True)
        if np.any(np.isclose(sums, 0.0)):
            raise ValueError("Reduced crystal directions must lie inside the IPF sector cone.")
        normalized = barycentric / sums
        colors = normalized @ self.color_vertices
        colors = np.clip(colors, 0.0, 1.0) ** (1.0 / self.saturation_gamma)
        max_values = np.max(colors, axis=1, keepdims=True)
        colors = colors / np.where(max_values > 0.0, max_values, 1.0)
        colors = np.ascontiguousarray(colors)
        colors.setflags(write=False)
        return colors

    def boundary_points_2d(
        self,
        *,
        method: str = "stereographic",
        samples_per_edge: int = 64,
    ) -> np.ndarray:
        """Projected sector-boundary polyline, for drawing the colour-key outline.
        """

        trace = self.sector.boundary_trace(samples_per_edge=samples_per_edge)
        projected = project_directions(trace, method=method, antipodal=self.antipodal)
        return freeze_array(np.ascontiguousarray(projected))

    def legend_mesh(
        self,
        *,
        resolution_deg: float = 2.0,
        method: str = "stereographic",
    ) -> tuple[np.ndarray, np.ndarray]:
        """A filled colour-key mesh over the fundamental sector.

        The standard-triangle legend that accompanies an IPF map and makes its
        colours interpretable.

        Returns
        -------
        tuple of np.ndarray
            Projected mesh coordinates and their RGB colours.
        """

        if not 0.0 < float(resolution_deg) <= 15.0:
            raise ValueError("legend_mesh resolution_deg must lie in (0, 15] degrees.")
        polar = np.arange(0.0, 90.0 + resolution_deg, resolution_deg)
        azimuth = np.arange(0.0, 360.0, resolution_deg)
        polar_grid, azimuth_grid = np.meshgrid(polar, azimuth, indexing="ij")
        directions = spherical_angles_to_directions(
            polar_grid.ravel(),
            azimuth_grid.ravel(),
        ).reshape(-1, 3)
        sector = self.sector
        mask = np.asarray(sector.contains(directions, atol=1e-6))
        selected = directions[mask]
        if selected.shape[0] == 0:
            raise ValueError("legend_mesh produced no in-sector directions.")
        points = project_directions(selected, method=method, antipodal=self.antipodal)
        colors = self.colors_from_crystal_directions(selected)
        return (
            freeze_array(np.ascontiguousarray(points)),
            freeze_array(np.ascontiguousarray(colors)),
        )

    def colors_from_orientations(self, orientations: OrientationSet) -> np.ndarray:
        """RGB colours for orientations relative to a specimen direction.

        Maps the chosen specimen axis into each crystal frame and colours the
        resulting crystal direction — the operation behind an IPF map.
        """

        if orientations.symmetry is None or orientations.symmetry != self.crystal_symmetry:
            raise ValueError("OrientationSet symmetry must match the IPFColorKey crystal symmetry.")
        crystal_directions = orientations.map_sample_directions_to_crystal(self.specimen_direction)
        crystal_direction_array = (
            crystal_directions.values
            if isinstance(crystal_directions, VectorSet)
            else crystal_directions
        )
        return self.colors_from_crystal_directions(crystal_direction_array)

    def color_from_orientation(self, orientation: Orientation) -> np.ndarray:
        """The RGB colour of a single orientation; see
        :meth:`colors_from_orientations`.
        """

        if orientation.symmetry is None or orientation.symmetry != self.crystal_symmetry:
            raise ValueError("Orientation symmetry must match the IPFColorKey crystal symmetry.")
        crystal_direction = orientation.map_sample_vector_to_crystal(self.specimen_direction)
        crystal_direction_array = (
            crystal_direction.values[0]
            if isinstance(crystal_direction, VectorSet)
            else crystal_direction
        )
        return as_float_array(
            self.colors_from_crystal_directions(crystal_direction_array[None, :])[0],
            shape=(3,),
        )


def _key_for_orientation(
    orientation: Orientation,
    *,
    direction: str | ArrayLike,
    key: IPFColorKey | None,
) -> IPFColorKey:
    if key is not None:
        return key
    if orientation.symmetry is None:
        raise ValueError("orientation.symmetry is required when key is not provided.")
    return IPFColorKey(
        crystal_symmetry=orientation.symmetry,
        specimen_direction=specimen_direction_vector(direction),
    )


def _key_for_orientations(
    orientations: OrientationSet,
    *,
    direction: str | ArrayLike,
    key: IPFColorKey | None,
) -> IPFColorKey:
    if key is not None:
        return key
    if orientations.symmetry is None:
        raise ValueError("orientations.symmetry is required when key is not provided.")
    return IPFColorKey(
        crystal_symmetry=orientations.symmetry,
        specimen_direction=specimen_direction_vector(direction),
    )


def ipf_color(
    orientation: Orientation,
    *,
    direction: str | ArrayLike = "ND",
    key: IPFColorKey | None = None,
) -> np.ndarray:
    """Return the IPF RGB color for one orientation and specimen direction."""

    return _key_for_orientation(orientation, direction=direction, key=key).color_from_orientation(
        orientation
    )


def ipf_colors(
    orientations: OrientationSet,
    *,
    direction: str | ArrayLike = "ND",
    key: IPFColorKey | None = None,
) -> np.ndarray:
    """Return IPF RGB colors for an orientation set and specimen direction."""

    return _key_for_orientations(
        orientations,
        direction=direction,
        key=key,
    ).colors_from_orientations(orientations)


def plot_ipf_key(
    key: IPFColorKey,
    *,
    method: str = "stereographic",
    resolution_deg: float = 1.0,
    marker_size: float = 4.0,
    ax: Any | None = None,
    title: str | None = None,
) -> tuple[Any, Any]:
    """Render the IPF color key legend for one crystal symmetry and direction."""

    import matplotlib.pyplot as plt

    if ax is None:
        fig, axis = plt.subplots(figsize=(4.0, 4.0))
    else:
        axis = ax
        fig = axis.figure
    points, colors = key.legend_mesh(resolution_deg=resolution_deg, method=method)
    axis.scatter(points[:, 0], points[:, 1], c=colors, s=marker_size, linewidths=0.0)
    boundary = key.boundary_points_2d(method=method)
    axis.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0)
    axis.set_aspect("equal")
    axis.set_axis_off()
    if title is None:
        title = f"IPF key {key.crystal_symmetry.point_group}"
    axis.set_title(title)
    return fig, axis


__all__ = [
    "IPFColorKey",
    "ipf_color",
    "ipf_colors",
    "plot_ipf_key",
]
