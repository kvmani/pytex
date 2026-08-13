"""Canonical S2 (unit-sphere) direction semantics.

Spherical-angle convention used across PyTex: the polar angle is measured from the
+Z axis of the owning reference frame, and the azimuth is measured from +X toward
+Y. Public angle arguments and returns are degrees unless a ``_rad`` suffix says
otherwise, matching the stereonet surface this module canonicalizes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, freeze_array, normalize_vectors
from pytex.core.batches import VectorSet
from pytex.core.provenance import ProvenanceRecord

if TYPE_CHECKING:
    from pytex.core.frames import ReferenceFrame
    from pytex.core.orientation import Rotation

_HEMISPHERES = ("upper", "sphere")


def spherical_angles_to_directions(
    polar_deg: ArrayLike,
    azimuth_deg: ArrayLike,
) -> np.ndarray:
    """Unit direction vectors from polar and azimuthal angles in degrees.

    Convention
    ----------
    Polar angle measured from ``+z``, azimuth measured from ``+x`` towards
    ``+y``. The two inputs are broadcast against one another, so a single
    polar angle can be paired with an array of azimuths to trace a small
    circle.

    Returns
    -------
    np.ndarray
        Unit vectors with a trailing dimension of 3.
    """

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
    return freeze_array(np.ascontiguousarray(directions, dtype=np.float64))


def directions_to_spherical_angles(
    directions: ArrayLike,
    *,
    antipodal: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Polar and azimuthal angles in degrees from direction vectors.

    The inverse of :func:`spherical_angles_to_directions`. Inputs are
    normalized internally.

    Parameters
    ----------
    directions : ArrayLike
        Any array ending in dimension 3.
    antipodal : bool
        Fold vectors onto the upper hemisphere first, so that ``+v`` and
        ``-v`` give the same angles — the right choice for plane normals and
        for any direction family quoted without a sense.

    Returns
    -------
    tuple of np.ndarray
        Polar angles in ``[0, 180]`` and azimuths in ``[0, 360)``, both
        read-only.
    """

    vectors = np.array(normalize_vectors(directions), copy=True)
    if antipodal:
        mask = vectors[..., 2] < 0.0
        vectors[mask] *= -1.0
    polar = np.rad2deg(np.arccos(np.clip(vectors[..., 2], -1.0, 1.0)))
    azimuth = np.mod(np.rad2deg(np.arctan2(vectors[..., 1], vectors[..., 0])), 360.0)
    polar = freeze_array(np.ascontiguousarray(polar, dtype=np.float64))
    azimuth = freeze_array(np.ascontiguousarray(azimuth, dtype=np.float64))
    return polar, azimuth


def _broadcast_unit_rows(
    left: np.ndarray,
    right: np.ndarray,
    *,
    operation: str,
) -> tuple[np.ndarray, np.ndarray]:
    if left.shape[0] == right.shape[0]:
        return left, right
    if left.shape[0] == 1:
        return np.broadcast_to(left, right.shape), right
    if right.shape[0] == 1:
        return left, np.broadcast_to(right, left.shape)
    raise ValueError(
        f"Cannot broadcast spherical vector sets of lengths {left.shape[0]} and "
        f"{right.shape[0]} for {operation}: lengths must match or one must be 1."
    )


@dataclass(frozen=True, slots=True)
class SphericalVectorSet:
    """A batch of unit directions with a frame and an antipodal declaration.

    Purpose
    -------
    Directional data — poles, plane normals, specimen axes — with the two
    facts needed to interpret it: which frame it lives in, and whether ``+v``
    and ``-v`` mean the same thing. The antipodal flag is not cosmetic: it
    changes how angles, means, and hemisphere folding are computed, and it is
    the difference between a plane normal (antipodal) and a slip direction
    (not).

    Attributes
    ----------
    values : np.ndarray
        ``(n, 3)`` unit vectors.
    reference_frame : ReferenceFrame
    antipodal : bool
        Whether a direction and its reverse are the same datum.
    provenance : ProvenanceRecord, optional
    """

    values: np.ndarray
    reference_frame: ReferenceFrame
    antipodal: bool = False
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", normalize_vectors(self.values))

    @classmethod
    def from_vectors(
        cls,
        vectors: ArrayLike,
        *,
        reference_frame: ReferenceFrame,
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> SphericalVectorSet:
        """Build a direction set from raw vectors, normalizing them.

        Parameters
        ----------
        vectors : ArrayLike
            ``(n, 3)`` vectors; normalized internally.
        reference_frame : ReferenceFrame
            The frame the directions live in. Required, because a direction
            without a frame cannot be compared with anything.
        antipodal : bool
            Declare ``+v`` and ``-v`` equivalent. This is a property of the
            data's meaning — true for plane normals, false for slip directions —
            and it changes how angles, means, and folding behave downstream.
        provenance : ProvenanceRecord, optional
        """

        return cls(
            values=np.asarray(vectors, dtype=np.float64),
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )

    @classmethod
    def from_vector_set(
        cls,
        vector_set: VectorSet,
        *,
        antipodal: bool = False,
    ) -> SphericalVectorSet:
        """Build a direction set from a :class:`~pytex.core.batches.VectorSet`.

        The frame is inherited; the values are normalized.
        """

        return cls(
            values=vector_set.values,
            reference_frame=vector_set.reference_frame,
            antipodal=antipodal,
            provenance=vector_set.provenance,
        )

    @classmethod
    def from_polar(
        cls,
        polar: ArrayLike,
        azimuth: ArrayLike,
        *,
        reference_frame: ReferenceFrame,
        degrees: bool = True,
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> SphericalVectorSet:
        """Build a direction set from polar and azimuthal angles.

        The angle convention is that of
        :func:`spherical_angles_to_directions`. Pass ``degrees=False`` for
        radians. The two angle arrays are broadcast against one another.
        """

        polar_values = np.atleast_1d(np.asarray(polar, dtype=np.float64))
        azimuth_values = np.atleast_1d(np.asarray(azimuth, dtype=np.float64))
        if not degrees:
            polar_values = np.rad2deg(polar_values)
            azimuth_values = np.rad2deg(azimuth_values)
        directions = spherical_angles_to_directions(polar_values, azimuth_values)
        return cls(
            values=directions.reshape(-1, 3),
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | SphericalVectorSet:
        selected = self.values[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return SphericalVectorSet(
            values=selected,
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        """The underlying ``(n, 3)`` unit-vector array, without frame meaning.
        """

        return self.values

    def to_vector_set(self) -> VectorSet:
        """The general :class:`~pytex.core.batches.VectorSet` view of these
        directions.

        Drops the antipodal flag, which the general vector type does not carry.
        """

        return VectorSet(
            values=self.values,
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def to_polar(self, *, degrees: bool = True) -> tuple[np.ndarray, np.ndarray]:
        """Polar and azimuthal angles of the directions.

        Honours the set's antipodal flag, folding onto the upper hemisphere
        first when it is set. Pass ``degrees=False`` for radians.
        """

        polar_deg, azimuth_deg = directions_to_spherical_angles(
            self.values,
            antipodal=self.antipodal,
        )
        if degrees:
            return polar_deg, azimuth_deg
        polar_rad = freeze_array(np.ascontiguousarray(np.deg2rad(polar_deg)))
        azimuth_rad = freeze_array(np.ascontiguousarray(np.deg2rad(azimuth_deg)))
        return polar_rad, azimuth_rad

    def subset(self, indices: ArrayLike) -> SphericalVectorSet:
        """The directions at the given indices, as a new set.

        Frame, antipodal flag, and provenance are preserved.
        """

        return SphericalVectorSet(
            values=self.values[np.asarray(indices)],
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )

    def fold_upper_hemisphere(self) -> SphericalVectorSet:
        """Fold every direction onto the upper hemisphere.

        Purpose
        -------
        The convention behind one-hemisphere pole figures: because ``+v`` and
        ``-v`` denote the same pole, only one hemisphere need be plotted.

        Directions on the equator, where the sign is not resolved by the ``z``
        component, are disambiguated by a deterministic rule on ``x`` and then
        ``y``, so the fold is stable rather than sensitive to numerical noise.

        Raises
        ------
        ValueError
            When the set is not antipodal. Folding a set whose directions carry
            a genuine sense would destroy information, so it is refused rather
            than performed silently.
        """

        if not self.antipodal:
            raise ValueError(
                "fold_upper_hemisphere is only meaningful for antipodal direction sets; "
                "construct the set with antipodal=True if +v and -v are equivalent."
            )
        folded = np.array(self.values, copy=True)
        lower = folded[:, 2] < 0.0
        equator_negative = (np.isclose(folded[:, 2], 0.0)) & (
            (folded[:, 0] < 0.0) | (np.isclose(folded[:, 0], 0.0) & (folded[:, 1] < 0.0))
        )
        folded[lower | equator_negative] *= -1.0
        return SphericalVectorSet(
            values=folded,
            reference_frame=self.reference_frame,
            antipodal=True,
            provenance=self.provenance,
        )

    def _require_matching_frame(self, other: SphericalVectorSet, *, operation: str) -> None:
        if self.reference_frame != other.reference_frame:
            raise ValueError(
                f"SphericalVectorSet {operation} requires both operands to share one "
                "reference frame."
            )

    def dot(self, other: SphericalVectorSet) -> np.ndarray:
        """Row-wise dot products with another direction set.

        Both sets must share a reference frame; one row broadcasts against many.
        Returns an ``(n,)`` read-only array of cosines.
        """

        self._require_matching_frame(other, operation="dot")
        left, right = _broadcast_unit_rows(self.values, other.values, operation="dot")
        return freeze_array(np.ascontiguousarray(np.einsum("ni,ni->n", left, right)))

    def cross(self, other: SphericalVectorSet) -> SphericalVectorSet:
        """Row-wise cross products with another direction set.

        The result is normalized and carries the antipodal flag if either
        operand does. Parallel or antiparallel pairs raise, because their cross
        product has no defined direction — better than returning a normalized
        zero vector that would look like a valid direction.
        """

        self._require_matching_frame(other, operation="cross")
        left, right = _broadcast_unit_rows(self.values, other.values, operation="cross")
        products = np.cross(left, right)
        norms = np.linalg.norm(products, axis=1)
        if np.any(np.isclose(norms, 0.0)):
            raise ValueError(
                "cross is undefined for parallel or antiparallel direction pairs."
            )
        return SphericalVectorSet(
            values=products,
            reference_frame=self.reference_frame,
            antipodal=self.antipodal or other.antipodal,
            provenance=self.provenance,
        )

    def angles_to_rad(self, other: SphericalVectorSet) -> np.ndarray:
        """Row-wise angles to another direction set, in radians.

        When either set is antipodal the cosine is taken in absolute value, so
        angles lie in ``[0, pi/2]`` and a direction is never reported as
        180 degrees from its own reverse. Otherwise angles lie in ``[0, pi]``.
        Both sets must share a reference frame.
        """

        self._require_matching_frame(other, operation="angles_to_rad")
        left, right = _broadcast_unit_rows(
            self.values,
            other.values,
            operation="angles_to_rad",
        )
        cosines = np.clip(np.einsum("ni,ni->n", left, right), -1.0, 1.0)
        if self.antipodal or other.antipodal:
            cosines = np.abs(cosines)
        return freeze_array(np.ascontiguousarray(np.arccos(cosines)))

    def angles_to_deg(self, other: SphericalVectorSet) -> np.ndarray:
        """Row-wise angles to another direction set, in degrees.

        See :meth:`angles_to_rad` for the antipodal handling.
        """

        return freeze_array(np.ascontiguousarray(np.rad2deg(self.angles_to_rad(other))))

    def orientation_tensor(self) -> np.ndarray:
        """The normalized orientation tensor ``T = (1/n) sum v_i v_i^T``.

        Purpose
        -------
        The second-moment descriptor of a direction distribution, and the
        standard tool of directional statistics: its eigenvalues classify the
        distribution as clustered, girdle-like, or uniform, and its eigenvectors
        give the principal directions. It is the right summary for antipodal
        data, where a vector mean would cancel.

        Returns
        -------
        np.ndarray
            ``(3, 3)`` symmetric tensor of unit trace, read-only.
        """

        tensor = np.einsum("ni,nj->ij", self.values, self.values) / float(len(self))
        return freeze_array(np.ascontiguousarray(tensor))

    def mean_direction(self) -> np.ndarray:
        """The mean direction of the set, as a unit vector.

        Method
        ------
        For an antipodal set the mean is the principal eigenvector of the
        :meth:`orientation_tensor`, sign-canonicalized to the upper hemisphere —
        the correct estimator when ``+v`` and ``-v`` are the same datum, since a
        vector sum would cancel. For a non-antipodal set it is the normalized
        resultant vector.
        """

        if self.antipodal:
            eigenvalues, eigenvectors = np.linalg.eigh(self.orientation_tensor())
            principal = eigenvectors[:, int(np.argmax(eigenvalues))]
            if principal[2] < 0.0 or (
                np.isclose(principal[2], 0.0)
                and (
                    principal[0] < 0.0
                    or (np.isclose(principal[0], 0.0) and principal[1] < 0.0)
                )
            ):
                principal = -principal
            return as_float_array(principal, shape=(3,))
        resultant = self.values.sum(axis=0)
        norm = float(np.linalg.norm(resultant))
        if np.isclose(norm, 0.0):
            raise ValueError(
                "Mean direction is undefined: the resultant vector is numerically zero."
            )
        return as_float_array(resultant / norm, shape=(3,))

    def rotated_by(self, rotation: Rotation) -> SphericalVectorSet:
        """The directions rotated by a :class:`~pytex.core.orientation.Rotation`.

        Frame, antipodal flag, and provenance are preserved. Note that the frame
        label is unchanged: this rotates directions within a frame, it does not
        re-express them in another one.
        """

        mapped = rotation.apply(self.values)
        return SphericalVectorSet(
            values=np.asarray(mapped, dtype=np.float64),
            reference_frame=self.reference_frame,
            antipodal=self.antipodal,
            provenance=self.provenance,
        )


def _require_hemisphere(hemisphere: str) -> str:
    if hemisphere not in _HEMISPHERES:
        supported = ", ".join(_HEMISPHERES)
        raise ValueError(
            f"Unsupported hemisphere '{hemisphere}'. Supported values: {supported}."
        )
    return hemisphere


def _require_resolution(resolution_deg: float) -> float:
    resolution = float(resolution_deg)
    if not 0.0 < resolution <= 90.0:
        raise ValueError("Grid resolution must lie in the interval (0, 90] degrees.")
    return resolution


def _ring_band_weights(
    polar_ring_deg: np.ndarray,
    counts: np.ndarray,
    *,
    polar_max_deg: float,
    half_band_deg: float,
) -> np.ndarray:
    lower = np.clip(polar_ring_deg - half_band_deg, 0.0, polar_max_deg)
    upper = np.clip(polar_ring_deg + half_band_deg, 0.0, polar_max_deg)
    band_measure = np.cos(np.deg2rad(lower)) - np.cos(np.deg2rad(upper))
    per_point = np.repeat(band_measure / counts, counts)
    return np.asarray(per_point / per_point.sum(), dtype=np.float64)


def raster_solid_angle_weights(
    polar_deg: ArrayLike,
    *,
    polar_max_deg: float | None = None,
) -> np.ndarray:
    """Integration weights for a measured polar/azimuth raster.

    Purpose
    -------
    A diffractometer samples a pole figure on a tilt/rotation raster, whose
    points are not equally spaced on the sphere: rings near the pole are short
    and carry little solid angle, rings near the equator are long and carry
    much. Averaging such a raster without weights over-counts the pole. These
    are the weights that make the average an integral — and therefore make a
    conversion to multiples of a random distribution correct.

    When to use
    -----------
    Whenever a mean, an integral, or an m.r.d. normalization is taken over
    measured raster data, e.g. as the ``integration_weights`` argument of
    :meth:`~pytex.texture.PoleFigure.spherical_mean`. Grids built by
    :class:`S2Grid` already carry their own weights and do not need this.

    Method
    ------
    Points are grouped into rings of equal polar angle. Each ring is given the
    solid angle of the band midway to its neighbours, ``cos(lower) -
    cos(upper)``, shared equally among its points. Bands are clipped to the
    measured polar range, so the weights describe the region actually measured;
    a partial pole figure — the usual case, since defocusing limits the
    reachable tilt — is therefore averaged over its measured cap, which equals
    the true spherical mean only if the unmeasured cap has the same mean.

    Parameters
    ----------
    polar_deg : ArrayLike
        Polar angle of every sampled point, in degrees. One entry per point, in
        the same order as the intensities they weight. Points sharing a polar
        angle are treated as one ring.
    polar_max_deg : float, optional
        Hard upper bound for the outermost band. Without it the outermost ring
        is extended by its own half step, which for a raster ending exactly at
        the equator pushes its band past 90 degrees and gives that ring close to
        twice the solid angle it owns. On a 5 degree hemispherical raster that
        alone biases the spherical mean of ``cos^2`` by 4 percent; passing
        ``polar_max_deg=90.0`` reduces the error to 0.06 percent. Pass it
        whenever the raster is known to stop at a boundary of the integration
        domain rather than merely at the last angle the instrument reached.

    Returns
    -------
    np.ndarray
        Read-only weights summing to 1, one per input point, strictly positive.

    Examples
    --------
    On a raster stepping the tilt by a constant amount, the equatorial ring
    carries far more weight than the ring next to the pole — the ratio of their
    band areas, not 1.
    """

    polar = as_float_array(np.asarray(polar_deg, dtype=np.float64).reshape(-1), shape=(None,))
    if polar.size == 0:
        raise ValueError("raster_solid_angle_weights requires at least one point.")
    if np.any(polar < -1e-9) or np.any(polar > 180.0 + 1e-9):
        raise ValueError("polar_deg must lie in [0, 180].")
    if polar_max_deg is None:
        outer_bound = 180.0
    else:
        outer_bound = float(polar_max_deg)
        if not 0.0 < outer_bound <= 180.0:
            raise ValueError("polar_max_deg must lie in (0, 180].")
        if np.any(polar > outer_bound + 1e-9):
            raise ValueError("polar_deg must not exceed polar_max_deg.")
    rings, inverse, counts = np.unique(
        np.round(polar, 9), return_inverse=True, return_counts=True
    )
    if rings.size == 1:
        # A single ring carries the whole measured band; every point on it is
        # equivalent, so the weights are uniform.
        return freeze_array(np.full(polar.size, 1.0 / polar.size, dtype=np.float64))
    midpoints = 0.5 * (rings[:-1] + rings[1:])
    half_step = np.diff(rings) / 2.0
    # Interior rings are bounded by the midpoints to their neighbours. The two
    # edge rings would otherwise get only a half band, which would under-weight
    # them; extend each outwards by its own half step, clipped to the sphere.
    # The weights then describe the cap actually measured.
    lower = np.concatenate([[max(rings[0] - half_step[0], 0.0)], midpoints])
    upper = np.concatenate([midpoints, [min(rings[-1] + half_step[-1], outer_bound)]])
    band = np.cos(np.deg2rad(lower)) - np.cos(np.deg2rad(upper))
    if np.any(band <= 0.0):  # pragma: no cover - only reachable for degenerate rings
        raise ValueError("polar_deg produced a ring of zero solid angle.")
    per_point = (band / counts)[inverse]
    return freeze_array(np.asarray(per_point / per_point.sum(), dtype=np.float64))


@dataclass(frozen=True, slots=True)
class S2Grid:
    """A sampling grid on the unit sphere, with per-point integration weights.

    Purpose
    -------
    The evaluation and integration support for pole figures and spherical
    functions. The weights are the point: sphere sampling is never uniform in
    a naive latitude-longitude scheme, so summing without weights biases
    every integral towards the poles.

    Attributes
    ----------
    vectors : SphericalVectorSet
        The grid directions.
    weights : np.ndarray
        Per-point integration weights; use these, never uniform weights.
    resolution_deg : float
        Nominal angular spacing.
    hemisphere : str
        ``"upper"`` or ``"sphere"``.
    method : str
        How the grid was generated — ``"equispaced"`` (near-equal-area, for
        integration) or ``"regular"`` (latitude-longitude, for display).
    """

    vectors: SphericalVectorSet
    weights: np.ndarray
    resolution_deg: float
    hemisphere: str
    method: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "weights", as_float_array(self.weights, shape=(None,)))
        object.__setattr__(self, "hemisphere", _require_hemisphere(self.hemisphere))
        object.__setattr__(self, "resolution_deg", _require_resolution(self.resolution_deg))
        if self.weights.shape[0] != len(self.vectors):
            raise ValueError("S2Grid.weights must have one weight per grid direction.")
        if np.any(self.weights <= 0.0):
            raise ValueError("S2Grid.weights must be strictly positive.")
        if not np.isclose(float(self.weights.sum()), 1.0, atol=1e-10):
            raise ValueError("S2Grid.weights must sum to 1 (normalized surface measure).")

    def __len__(self) -> int:
        return len(self.vectors)

    @classmethod
    def equispaced(
        cls,
        resolution_deg: float,
        *,
        reference_frame: ReferenceFrame,
        hemisphere: str = "upper",
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> S2Grid:
        """A near-equal-area grid on the sphere at the requested resolution.

        Purpose
        -------
        The integration and evaluation support for pole figures and spherical
        functions. Equal-area sampling is what makes summing over grid points a
        valid approximation to a spherical integral; a naive latitude-longitude
        grid over-samples the poles badly and biases every such sum.

        Method
        ------
        Points are placed on rings of constant polar angle, with the number of
        points per ring scaled by the ring circumference ``sin(theta)``, so cell
        areas are nearly equal. Per-point weights are computed from the ring
        band areas and returned with the grid, so callers integrate with the
        weights rather than assuming uniformity.

        Parameters
        ----------
        resolution_deg : float
            Target angular spacing.
        reference_frame : ReferenceFrame
            Frame the grid directions live in.
        hemisphere : str
            ``"upper"`` (default) or ``"sphere"``.
        antipodal : bool
            Declare ``+v`` and ``-v`` equivalent on the resulting directions.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        S2Grid
            Directions, integration weights, and the generation metadata.

        See Also
        --------
        regular : A latitude-longitude grid, for display rather than integration.
        """

        resolution = _require_resolution(resolution_deg)
        hemisphere = _require_hemisphere(hemisphere)
        polar_max = 90.0 if hemisphere == "upper" else 180.0
        ring_count = round(polar_max / resolution)
        polar_rings = np.linspace(0.0, polar_max, ring_count + 1)

        ring_polar: list[np.ndarray] = []
        ring_azimuth: list[np.ndarray] = []
        counts = np.empty(polar_rings.shape[0], dtype=np.int64)
        for ring_index, polar_deg in enumerate(polar_rings):
            circumference_deg = 360.0 * float(np.sin(np.deg2rad(polar_deg)))
            count = max(1, round(circumference_deg / resolution))
            counts[ring_index] = count
            azimuths = np.linspace(0.0, 360.0, count, endpoint=False)
            ring_polar.append(np.full(count, float(polar_deg)))
            ring_azimuth.append(azimuths)

        polar_all = np.concatenate(ring_polar)
        azimuth_all = np.concatenate(ring_azimuth)
        weights = _ring_band_weights(
            polar_rings,
            counts,
            polar_max_deg=polar_max,
            half_band_deg=resolution / 2.0,
        )
        vectors = SphericalVectorSet.from_polar(
            polar_all,
            azimuth_all,
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )
        return cls(
            vectors=vectors,
            weights=weights,
            resolution_deg=resolution,
            hemisphere=hemisphere,
            method="equispaced",
        )

    @classmethod
    def regular(
        cls,
        polar_step_deg: float,
        azimuth_step_deg: float,
        *,
        reference_frame: ReferenceFrame,
        hemisphere: str = "upper",
        antipodal: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> S2Grid:
        """A regular latitude-longitude grid on the sphere.

        Purpose
        -------
        A rectangular grid in ``(polar, azimuth)``, convenient for contouring
        and for interchange with tools that expect gridded data.

        Sampling density is *not* uniform on the sphere — cells shrink towards
        the poles — so the returned weights must be used for any integration.
        Prefer :meth:`equispaced` when the grid is an integration support rather
        than a display raster.

        Parameters
        ----------
        polar_step_deg, azimuth_step_deg : float
            Grid spacing along each angular axis.
        reference_frame : ReferenceFrame
        hemisphere : str
            ``"upper"`` (default) or ``"sphere"``.
        antipodal : bool
            Declare ``+v`` and ``-v`` equivalent on the resulting directions.
        provenance : ProvenanceRecord, optional
        """

        polar_step = _require_resolution(polar_step_deg)
        azimuth_step = float(azimuth_step_deg)
        if not 0.0 < azimuth_step <= 120.0 or not np.isclose(
            np.mod(360.0, azimuth_step), 0.0, atol=1e-10
        ):
            raise ValueError(
                "azimuth_step_deg must lie in (0, 120] and divide 360 degrees evenly."
            )
        hemisphere = _require_hemisphere(hemisphere)
        polar_max = 90.0 if hemisphere == "upper" else 180.0
        ring_count = round(polar_max / polar_step)
        polar_rings = np.linspace(0.0, polar_max, ring_count + 1)
        azimuth_count = round(360.0 / azimuth_step)

        ring_polar: list[np.ndarray] = []
        ring_azimuth: list[np.ndarray] = []
        counts = np.empty(polar_rings.shape[0], dtype=np.int64)
        for ring_index, polar_deg in enumerate(polar_rings):
            at_pole = np.isclose(float(polar_deg), 0.0) or np.isclose(float(polar_deg), 180.0)
            count = 1 if at_pole else azimuth_count
            counts[ring_index] = count
            azimuths = np.linspace(0.0, 360.0, count, endpoint=False)
            ring_polar.append(np.full(count, float(polar_deg)))
            ring_azimuth.append(azimuths)

        polar_all = np.concatenate(ring_polar)
        azimuth_all = np.concatenate(ring_azimuth)
        weights = _ring_band_weights(
            polar_rings,
            counts,
            polar_max_deg=polar_max,
            half_band_deg=polar_step / 2.0,
        )
        vectors = SphericalVectorSet.from_polar(
            polar_all,
            azimuth_all,
            reference_frame=reference_frame,
            antipodal=antipodal,
            provenance=provenance,
        )
        return cls(
            vectors=vectors,
            weights=weights,
            resolution_deg=polar_step,
            hemisphere=hemisphere,
            method="regular",
        )
