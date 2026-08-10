"""Kikuchi band geometry and the gnomonic projection.

This module supplies the geometric layer shared by EBSD and TEM Kikuchi
patterns. A Kikuchi band is not a spot: it is the pair of curves in which the
two Kossel cones of a lattice plane meet the detector, and reading a pattern is
a geometric exercise on those curves and their intersections.

Physical basis
--------------
Electrons scattered incoherently inside the specimen illuminate every direction.
For a lattice plane with interplanar spacing :math:`d`, those travelling at the
Bragg angle :math:`\\theta_B` to the plane satisfy the diffraction condition,
where

.. math:: \\sin\\theta_B = \\frac{\\lambda}{2 d}.

The diffracting directions therefore form two cones of semi-angle
:math:`90^\\circ - \\theta_B` about the plane normal — the Kossel cones — one on
each side of the plane. Their intersections with the detector are the two edges
of the band, and the plane's own trace runs midway between them. The band's
angular width is :math:`2\\theta_B = 2\\arcsin(\\lambda/2d) \\approx \\lambda/d`, so a
*wide* band means a *small* d-spacing. Band widths measure the lattice directly,
and the widest bands in a pattern come from the high-index planes while the
strongest come from the low-index ones -- two orderings that are easy to
conflate.

Why gnomonic coordinates
------------------------
In the gnomonic projection — central projection from the diffraction source onto
a plane — every great circle maps to a straight line. The trace of a lattice
plane is a great circle, so **band centre lines are exactly straight in gnomonic
coordinates** whatever the detector tilt. The Kossel cones are not great circles,
so band *edges* are conics (hyperbolae for the small Bragg angles of electron
diffraction); the common textbook statement that band edges are straight is the
small-angle approximation, and PyTex does not make it. Edges here are sampled on
the exact cones.

Zone axes project to the points where the bands that share them intersect, which
is what makes a Kikuchi pattern readable as a map of orientation space.

Limits
------
Geometric and kinematic. Band positions and widths are exact for the stated
geometry, but intensities are a kinematic proxy: the excess/deficiency asymmetry
across a band, dynamical contrast, and higher-order Laue zone effects are not
modelled. Band positions — not the intensities — are what orientation
determination uses.

See Also
--------
pytex.diffraction.models.DiffractionGeometry : The detector and beam geometry
    this module builds on.
pytex.diffraction.saed : Spot (zone-axis) patterns, the coherent counterpart.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import FloatArray, as_float_array, freeze_array, normalize_vectors
from pytex.core.lattice import Phase
from pytex.core.miller import MillerPlane
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.orientation import Orientation
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.kinematic import centering_allowed_mask, electron_structure_factors
from pytex.diffraction.models import DiffractionGeometry
from pytex.diffraction.physics import ReflectionCondition

#: Directions closer than this to the detector plane are treated as not
#: projecting. The gnomonic coordinate diverges as the denominator approaches
#: zero, so a floor is required rather than optional.
_PROJECTION_EPSILON = 1e-9

#: Tolerance in pixels when testing whether a projected point lies on the
#: detector. A point on the exact detector edge accumulates a few units of
#: floating-point error passing through the projection, and must still be
#: reported as on-detector.
_EDGE_TOLERANCE_PX = 1e-6

#: An edge within this angle of grazing incidence projects to an unbounded
#: gnomonic offset, so the band width is reported as infinite rather than as a
#: very large finite number that would read as a measurement.
_GRAZING_ANGLE_ATOL = 1e-9


def _detector_components(geometry: DiffractionGeometry, directions_lab: ArrayLike) -> FloatArray:
    """Decompose laboratory directions onto the detector basis ``(u, v, n)``."""

    basis = geometry.detector_basis_lab
    directions = normalize_vectors(directions_lab)
    return as_float_array(directions @ basis, shape=(directions.shape[0], 3))


@dataclass(frozen=True, slots=True)
class GnomonicProjection:
    """Central projection from the diffraction source onto the detector plane.

    Purpose
    -------
    The coordinate system in which Kikuchi geometry is simple. A direction
    :math:`\\hat{r}` maps to the point where the ray from the source along
    :math:`\\hat{r}` meets the detector, expressed in units of the detector
    distance. Its defining property is that **great circles map to straight
    lines**, which is why lattice-plane traces are straight here and curved in
    pixel coordinates on a tilted detector.

    The origin of gnomonic coordinates is the pattern centre, and one gnomonic
    unit is one detector distance, so a coordinate of ``1.0`` lies at 45 degrees
    from the detector normal. Coordinates are therefore detector-size
    independent: the same crystal gives the same gnomonic pattern on any
    detector, which is what makes them the right frame for comparing or fitting
    patterns.

    Attributes
    ----------
    geometry : DiffractionGeometry
        Supplies the detector basis, distance, pattern centre, and pixel size.

    Notes
    -----
    Only the forward hemisphere projects: directions travelling away from the
    detector, or parallel to its plane, have no intersection and are reported as
    invalid rather than being given a spurious coordinate.
    """

    geometry: DiffractionGeometry

    def project_directions(self, directions_lab: ArrayLike) -> tuple[FloatArray, np.ndarray]:
        """Gnomonic coordinates of laboratory directions.

        Parameters
        ----------
        directions_lab : ArrayLike
            ``(n, 3)`` directions in the laboratory frame; normalized
            internally.

        Returns
        -------
        tuple of (FloatArray, np.ndarray)
            ``(n, 2)`` gnomonic coordinates and an ``(n,)`` validity mask.
            Directions that do not meet the detector plane are returned as
            ``NaN`` and flagged ``False``, so an invalid coordinate can never be
            mistaken for a position near the origin.
        """

        components = _detector_components(self.geometry, directions_lab)
        depth = components[:, 2]
        valid = depth > _PROJECTION_EPSILON
        coordinates = np.full((components.shape[0], 2), np.nan, dtype=np.float64)
        coordinates[valid, 0] = components[valid, 0] / depth[valid]
        coordinates[valid, 1] = components[valid, 1] / depth[valid]
        return freeze_array(coordinates), freeze_array(valid.astype(bool))

    def unproject(self, coordinates: ArrayLike) -> FloatArray:
        """Laboratory directions of gnomonic coordinates.

        The inverse of :meth:`project_directions` on the forward hemisphere.
        Every gnomonic coordinate corresponds to exactly one such direction, so
        this inverse is exact rather than approximate.

        Parameters
        ----------
        coordinates : ArrayLike
            ``(n, 2)`` gnomonic coordinates.

        Returns
        -------
        FloatArray
            ``(n, 3)`` unit directions in the laboratory frame.
        """

        gnomonic = as_float_array(coordinates, shape=(None, 2))
        basis = self.geometry.detector_basis_lab
        rays = (
            gnomonic[:, [0]] * basis[:, 0][None, :]
            + gnomonic[:, [1]] * basis[:, 1][None, :]
            + basis[:, 2][None, :]
        )
        return freeze_array(normalize_vectors(rays))

    def to_detector_px(self, coordinates: ArrayLike) -> FloatArray:
        """Detector pixel coordinates of gnomonic coordinates.

        Scales by the detector distance and the per-axis pixel size, then offsets
        to the pattern centre — so non-square pixels and an off-centre pattern
        centre are both handled.
        """

        gnomonic = as_float_array(coordinates, shape=(None, 2))
        pixel_size_mm = np.array(self.geometry.detector_pixel_size_um, dtype=np.float64) / 1000.0
        offsets_mm = gnomonic * self.geometry.camera_length_mm
        pixels = np.column_stack(
            [
                offsets_mm[:, 0] / pixel_size_mm[0] + self.geometry.pattern_center_px[0],
                offsets_mm[:, 1] / pixel_size_mm[1] + self.geometry.pattern_center_px[1],
            ]
        )
        return freeze_array(np.ascontiguousarray(pixels))

    def from_detector_px(self, coordinates_px: ArrayLike) -> FloatArray:
        """Gnomonic coordinates of detector pixel coordinates.

        The inverse of :meth:`to_detector_px`, and the first step in analysing a
        measured pattern: converting detected band positions into the frame where
        band centre lines are straight.
        """

        pixels = as_float_array(coordinates_px, shape=(None, 2))
        offsets_mm = self.geometry.detector_coordinates_mm(pixels)
        return freeze_array(
            np.ascontiguousarray(np.asarray(offsets_mm) / self.geometry.camera_length_mm)
        )

    def detector_corner_coordinates(self) -> FloatArray:
        """Gnomonic coordinates of the four detector corners.

        The gnomonic extent of the physical detector, which bounds where bands
        and zone axes can actually be observed. Corner order is
        ``(0, 0)``, ``(w-1, 0)``, ``(w-1, h-1)``, ``(0, h-1)`` in pixels.
        """

        height, width = self.geometry.detector_shape
        corners_px = np.array(
            [
                [0.0, 0.0],
                [float(width - 1), 0.0],
                [float(width - 1), float(height - 1)],
                [0.0, float(height - 1)],
            ],
            dtype=np.float64,
        )
        return self.from_detector_px(corners_px)

    def contains(
        self, coordinates: ArrayLike, *, atol_px: float = _EDGE_TOLERANCE_PX
    ) -> np.ndarray:
        """Whether gnomonic coordinates fall on the physical detector.

        Converts to pixels and tests against the detector shape, so the answer
        reflects the real sensor rather than an idealized infinite plane.

        Parameters
        ----------
        coordinates : ArrayLike
            ``(n, 2)`` gnomonic coordinates.
        atol_px : float
            Tolerance in pixels on the detector edges. A point on the exact
            boundary must be accepted: it round-trips through the projection
            with a few units of floating-point error, and rejecting it would
            make :meth:`detector_corner_coordinates` report its own corners as
            off-detector.

        Returns
        -------
        np.ndarray
            ``(n,)`` boolean mask, read-only. Non-finite coordinates — the
            output of a failed projection — are rejected.
        """

        if atol_px < 0.0:
            raise ValueError("atol_px must be non-negative.")
        pixels = self.to_detector_px(coordinates)
        height, width = self.geometry.detector_shape
        inside = (
            np.isfinite(pixels).all(axis=1)
            & (pixels[:, 0] >= -atol_px)
            & (pixels[:, 0] <= float(width - 1) + atol_px)
            & (pixels[:, 1] >= -atol_px)
            & (pixels[:, 1] <= float(height - 1) + atol_px)
        )
        return freeze_array(inside.astype(bool))


@dataclass(frozen=True, slots=True)
class KikuchiBand:
    """One Kikuchi band: a lattice plane seen as a pair of Kossel-cone traces.

    Purpose
    -------
    Carries everything needed to draw or measure a band: which plane produced
    it, the plane normal in the laboratory frame, its Bragg angle, and hence its
    angular width. The geometry in the detector plane is derived on demand from
    a :class:`GnomonicProjection`, so one band object serves any detector.

    Reading a band
    --------------
    The band's angular width is exactly ``2 * bragg_angle_rad``, and
    :math:`\\sin\\theta_B = \\lambda / 2d`, so band width is a direct measurement
    of the interplanar spacing -- an *inverse* one: the width grows as the spacing
    falls. Its centre line is the trace of the lattice plane itself.

    Attributes
    ----------
    plane : MillerPlane
        The diffracting plane, carrying its phase.
    plane_normal_lab : np.ndarray
        Unit plane normal in the laboratory frame, after the crystal orientation
        and the specimen-to-laboratory transform have been applied.
    bragg_angle_rad : float
        :math:`\\theta_B`, in ``(0, pi/2)``. Half the band's angular width.
    d_spacing_angstrom : float
        Interplanar spacing of the plane.
    intensity : float
        Relative band intensity from a kinematic ``|F|^2`` proxy. Indicative
        only; see the module limits.
    """

    plane: MillerPlane
    plane_normal_lab: np.ndarray
    bragg_angle_rad: float
    d_spacing_angstrom: float
    intensity: float

    def __post_init__(self) -> None:
        normal = as_float_array(self.plane_normal_lab, shape=(3,))
        norm = float(np.linalg.norm(normal))
        if not np.isfinite(norm) or np.isclose(norm, 0.0):
            raise ValueError("KikuchiBand.plane_normal_lab must be a finite non-zero vector.")
        object.__setattr__(self, "plane_normal_lab", freeze_array(normal / norm))
        if not 0.0 < self.bragg_angle_rad < 0.5 * np.pi:
            raise ValueError("KikuchiBand.bragg_angle_rad must lie in (0, pi/2).")
        if not np.isfinite(self.d_spacing_angstrom) or self.d_spacing_angstrom <= 0.0:
            raise ValueError("KikuchiBand.d_spacing_angstrom must be finite and positive.")
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("KikuchiBand.intensity must be finite and non-negative.")

    @property
    def angular_width_rad(self) -> float:
        """Full angular width of the band, ``2 * theta_B``, in radians.

        This is the quantity a band-width measurement yields, and it determines
        the d-spacing through Bragg's law.
        """

        return float(2.0 * self.bragg_angle_rad)

    def center_line_coefficients(self, projection: GnomonicProjection) -> FloatArray:
        """Coefficients ``(a, b, c)`` of the band centre line ``a x + b y + c = 0``.

        Purpose
        -------
        The band centre is the trace of the lattice plane, which is a great
        circle, and great circles are exactly straight lines in gnomonic
        coordinates. The coefficients are simply the plane normal expressed in
        the detector basis — a direct statement of that property, and the reason
        band-detection algorithms work in gnomonic coordinates.

        Returns
        -------
        FloatArray
            ``(3,)`` coefficients, normalized so that ``a^2 + b^2 + c^2 = 1``.
        """

        components = _detector_components(projection.geometry, self.plane_normal_lab[None, :])
        return freeze_array(as_float_array(components[0], shape=(3,)))

    def _cone_directions(self, *, sign: float, samples: int) -> FloatArray:
        """Sample one Kossel cone about the plane normal.

        The cone of directions at angle ``pi/2 - theta_B`` from the plane, i.e.
        with ``r . n = sign * sin(theta_B)``. Sampling the exact cone — rather
        than solving the projected conic — keeps the edges correct at any Bragg
        angle and needs no branch bookkeeping.
        """

        normal = self.plane_normal_lab
        trial = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        if abs(float(normal @ trial)) > 0.9:
            trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        first = np.cross(normal, trial)
        first = first / float(np.linalg.norm(first))
        second = np.cross(normal, first)
        cos_polar = sign * float(np.sin(self.bragg_angle_rad))
        sin_polar = float(np.cos(self.bragg_angle_rad))
        angles = np.linspace(0.0, 2.0 * np.pi, int(samples), endpoint=True)
        directions = (
            cos_polar * normal[None, :]
            + sin_polar * (np.cos(angles)[:, None] * first[None, :])
            + sin_polar * (np.sin(angles)[:, None] * second[None, :])
        )
        return as_float_array(directions, shape=(angles.shape[0], 3))

    def center_trace(
        self,
        projection: GnomonicProjection,
        *,
        samples: int = 361,
    ) -> FloatArray:
        """The band centre line, sampled in gnomonic coordinates.

        Parameters
        ----------
        projection : GnomonicProjection
        samples : int
            Points sampled around the great circle; at least two.

        Returns
        -------
        FloatArray
            ``(m, 2)`` gnomonic coordinates of the points that project, in
            order around the circle. Points on the far hemisphere are dropped,
            so ``m`` may be smaller than ``samples``.
        """

        if samples < 2:
            raise ValueError("samples must be at least two.")
        normal = self.plane_normal_lab
        trial = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        if abs(float(normal @ trial)) > 0.9:
            trial = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        first = np.cross(normal, trial)
        first = first / float(np.linalg.norm(first))
        second = np.cross(normal, first)
        angles = np.linspace(0.0, 2.0 * np.pi, int(samples), endpoint=True)
        directions = np.cos(angles)[:, None] * first[None, :] + (
            np.sin(angles)[:, None] * second[None, :]
        )
        coordinates, valid = projection.project_directions(directions)
        return freeze_array(np.ascontiguousarray(coordinates[valid]))

    def edge_traces(
        self,
        projection: GnomonicProjection,
        *,
        samples: int = 361,
    ) -> tuple[FloatArray, FloatArray]:
        """Both band edges, sampled in gnomonic coordinates.

        Purpose
        -------
        The two Kossel-cone intersections that bound the band. These are conic
        sections — hyperbolae at the small Bragg angles of electron
        diffraction — not straight lines; they are computed from the exact
        cones, so no small-angle approximation is made.

        Parameters
        ----------
        projection : GnomonicProjection
        samples : int
            Points sampled around each cone; at least two.

        Returns
        -------
        tuple of FloatArray
            Two ``(m, 2)`` arrays for the two edges. Either may be empty when
            that cone misses the detector hemisphere entirely.
        """

        if samples < 2:
            raise ValueError("samples must be at least two.")
        traces: list[FloatArray] = []
        for sign in (1.0, -1.0):
            directions = self._cone_directions(sign=sign, samples=samples)
            coordinates, valid = projection.project_directions(directions)
            traces.append(freeze_array(np.ascontiguousarray(coordinates[valid])))
        return traces[0], traces[1]

    def width_at_pattern_center(self, projection: GnomonicProjection) -> float:
        """Band width in gnomonic units, measured across the pattern centre.

        Purpose
        -------
        The observable width a band-detection routine measures. It is *not* a
        fixed multiple of the Bragg angle: a band far from the pattern centre
        appears wider in gnomonic units because the projection stretches with
        distance from the origin. Reporting it explicitly avoids the common
        error of converting an apparent width to a d-spacing without accounting
        for the band's position.

        Returns
        -------
        float
            The gnomonic distance between the two edges, measured along the
            line through the pattern centre perpendicular to the band.
            ``inf`` when the band lies edge-on to the detector, or when an edge
            grazes the detector plane, where the width is unbounded.

        Notes
        -----
        Derivation. Write the unit band normal in the detector basis as
        ``(a, b, c)``, with ``q = sqrt(a^2 + b^2)`` its in-plane magnitude and
        ``c`` its component along the detector normal. Restricted to the plane
        spanned by the detector normal and the in-plane direction, a unit ray at
        angle ``psi`` from the detector normal has gnomonic offset ``tan(psi)``
        and satisfies ``r . n = q sin(psi) + c cos(psi) = sin(psi + phi)`` with
        ``phi = atan2(c, q)``, since the normal is a unit vector. The Kossel
        condition ``r . n = +/- sin(theta_B)`` therefore gives
        ``psi = +/- theta_B - phi``, and the two edge offsets are the tangents
        of those angles.
        """

        coefficients = self.center_line_coefficients(projection)
        in_plane_norm = float(np.linalg.norm(coefficients[:2]))
        depth = float(coefficients[2])
        if np.isclose(in_plane_norm, 0.0):
            # The band normal lies along the detector normal, so the plane is
            # edge-on: its trace never meets the detector at finite distance.
            return float("inf")
        phase = float(np.arctan2(depth, in_plane_norm))
        offsets: list[float] = []
        for sign in (1.0, -1.0):
            angle = sign * self.bragg_angle_rad - phase
            if abs(abs(angle) - 0.5 * np.pi) < _GRAZING_ANGLE_ATOL:
                return float("inf")
            offsets.append(float(np.tan(angle)))
        return float(abs(offsets[0] - offsets[1]))

    def describe(self) -> str:
        """Convention-explicit prose summary of this band.

        States the plane, its d-spacing, the Bragg angle and the resulting
        angular width, and the relative intensity — with the notation and units
        made explicit, per the explainable-results doctrine.
        """

        indices = tuple(int(value) for value in self.plane.indices)
        return (
            f"Kikuchi band from {format_plane_indices(indices, style='plain')} of "
            f"{self.plane.phase.name}: d = {self.d_spacing_angstrom:.4f} A, Bragg angle "
            f"{np.rad2deg(self.bragg_angle_rad):.4f} deg, angular width "
            f"{np.rad2deg(self.angular_width_rad):.4f} deg, relative intensity "
            f"{self.intensity:.4f}. The centre line is the plane trace and is exactly "
            "straight in gnomonic coordinates; the edges are Kossel-cone conics."
        )


@dataclass(frozen=True, slots=True)
class KikuchiZoneAxis:
    """A zone axis and the gnomonic point where its bands intersect.

    Purpose
    -------
    Zone axes are the landmarks of a Kikuchi pattern. Every band whose plane
    contains the axis passes through its projected point, so a zone axis appears
    as a hub where several bands meet — the feature by which patterns are
    recognized by eye and indexed by algorithm.

    Attributes
    ----------
    indices : np.ndarray
        The zone-axis direction ``[uvw]`` in crystal indices.
    direction_lab : np.ndarray
        Unit direction in the laboratory frame.
    coordinates : np.ndarray
        ``(2,)`` gnomonic coordinates of the axis.
    band_count : int
        How many of the simulated bands contain this axis. A higher count means
        a more prominent, more easily recognized hub.
    on_detector : bool
        Whether the projected point falls on the physical detector. Zone axes
        just outside it still organize the visible bands, so they are reported
        rather than discarded.
    """

    indices: np.ndarray
    direction_lab: np.ndarray
    coordinates: np.ndarray
    band_count: int
    on_detector: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", freeze_array(np.rint(self.indices).astype(np.int64)))
        object.__setattr__(
            self, "direction_lab", freeze_array(as_float_array(self.direction_lab, shape=(3,)))
        )
        object.__setattr__(
            self, "coordinates", freeze_array(as_float_array(self.coordinates, shape=(2,)))
        )
        if self.band_count < 0:
            raise ValueError("KikuchiZoneAxis.band_count must be non-negative.")

    def describe(self) -> str:
        """Prose summary: the axis, where it projects, and how many bands meet."""

        indices = tuple(int(value) for value in self.indices)
        location = "on the detector" if self.on_detector else "outside the detector"
        return (
            f"Zone axis {format_direction_indices(indices, style='plain')} projects to gnomonic "
            f"({self.coordinates[0]:.4f}, {self.coordinates[1]:.4f}), {location}, with "
            f"{self.band_count} band(s) intersecting there."
        )


@dataclass(frozen=True, slots=True)
class KikuchiPattern:
    """A simulated Kikuchi pattern: bands, zone axes, and the geometry behind them.

    Purpose
    -------
    The geometric model of an EBSD or TEM Kikuchi pattern for a known phase,
    orientation, and detector. It is the forward model that orientation
    determination inverts, and the reference against which a measured pattern is
    compared.

    Limits
    ------
    Band positions and widths are exact for the stated geometry. Intensities are
    a kinematic proxy: excess/deficiency asymmetry across a band, dynamical
    contrast, and background are not modelled. Orientation determination uses
    band positions, so this limitation does not affect the geometric use.

    Attributes
    ----------
    geometry : DiffractionGeometry
        Detector, beam, and frame definitions.
    phase : Phase
    orientation : Orientation
        Crystal-to-specimen orientation; combined with the geometry's
        specimen-to-laboratory transform to place the bands.
    bands : tuple of KikuchiBand
        Ordered by decreasing intensity, so the most prominent bands come first.
    zone_axes : tuple of KikuchiZoneAxis
        Ordered by decreasing band count.
    wavelength_angstrom : float
        The relativistically corrected electron wavelength used.
    provenance : ProvenanceRecord, optional
    """

    geometry: DiffractionGeometry
    phase: Phase
    orientation: Orientation
    bands: tuple[KikuchiBand, ...]
    zone_axes: tuple[KikuchiZoneAxis, ...] = ()
    wavelength_angstrom: float = 0.0
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "bands", tuple(self.bands))
        object.__setattr__(self, "zone_axes", tuple(self.zone_axes))
        if self.wavelength_angstrom < 0.0:
            raise ValueError("KikuchiPattern.wavelength_angstrom must be non-negative.")

    @property
    def projection(self) -> GnomonicProjection:
        """The gnomonic projection of this pattern's detector geometry."""

        return GnomonicProjection(geometry=self.geometry)

    @property
    def band_count(self) -> int:
        """Number of simulated bands."""

        return len(self.bands)

    def band_for_plane(self, indices: ArrayLike) -> KikuchiBand | None:
        """The band produced by a given ``(hkl)``, or ``None`` if not simulated.

        Matching is on the exact index triple as simulated; a plane whose band
        was filtered out by the d-spacing or intensity cut-offs returns
        ``None`` rather than raising.
        """

        target = tuple(int(value) for value in np.rint(np.asarray(indices, dtype=np.float64)))
        for band in self.bands:
            if tuple(int(value) for value in band.plane.indices) == target:
                return band
        return None

    def widest_bands(self, count: int = 5) -> tuple[KikuchiBand, ...]:
        """The ``count`` widest bands, largest d-spacing first.

        The widest bands are those of the largest interplanar spacings, and they
        are the ones a band-detection routine finds most reliably — so this is
        the natural set to check a simulation against a measured pattern.
        """

        if count < 0:
            raise ValueError("count must be non-negative.")
        ordered = sorted(self.bands, key=lambda band: -band.d_spacing_angstrom)
        return tuple(ordered[:count])

    def describe(self) -> str:
        """Convention-explicit prose summary of the pattern.

        States the phase, the orientation in Bunge Euler angles, the beam energy
        and wavelength, the band and zone-axis counts, the widest bands with
        their d-spacings and angular widths, and the kinematic-intensity
        limitation. Angles are in degrees and indices are three-index.
        """

        phi1, big_phi, phi2 = self.orientation.rotation.to_bunge_euler(degrees=True)
        lines = [
            f"Kikuchi pattern for {self.phase.name} at Bunge Euler angles "
            f"({phi1:.2f}, {big_phi:.2f}, {phi2:.2f}) deg, beam energy "
            f"{self.geometry.beam_energy_kev:.1f} keV (electron wavelength "
            f"{self.wavelength_angstrom:.5f} A).",
            f"{self.band_count} band(s) and {len(self.zone_axes)} zone axis/axes were simulated. "
            "Band centre lines are exactly straight in gnomonic coordinates; band edges are "
            "Kossel-cone conics, and band angular width is 2*theta_B with "
            "sin(theta_B) = lambda / 2d.",
        ]
        for band in self.widest_bands(4):
            indices = tuple(int(value) for value in band.plane.indices)
            lines.append(
                f"  {format_plane_indices(indices, style='plain')}: d = "
                f"{band.d_spacing_angstrom:.4f} A, angular width "
                f"{np.rad2deg(band.angular_width_rad):.4f} deg."
            )
        lines.append(
            "Intensities are a kinematic |F|^2 proxy: excess/deficiency asymmetry across a "
            "band, dynamical contrast, and background are not modelled. Band positions and "
            "widths are exact for the stated geometry."
        )
        return "\n".join(lines)


def _antipodal_reduced_triples(max_index: int) -> np.ndarray:
    """Integer index triples up to ``max_index``, one per antipodal pair.

    Used for both planes and zone axes: ``(hkl)`` and ``(-h-k-l)`` describe the
    same plane and hence the same band, and ``[uvw]`` and ``[-u-v-w]`` describe
    the same zone. The representative is chosen by a deterministic sign rule on
    the first nonzero component, so the enumeration is reproducible.
    """

    span = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(span, span, span, indexing="ij"), axis=-1).reshape(-1, 3)
    nonzero = grid[np.any(grid != 0, axis=1)]
    # Keep one of each antipodal pair by a deterministic sign rule on the first
    # nonzero component.
    first_nonzero = np.argmax(nonzero != 0, axis=1)
    leading = nonzero[np.arange(nonzero.shape[0]), first_nonzero]
    return np.ascontiguousarray(nonzero[leading > 0])


def simulate_kikuchi_pattern(
    geometry: DiffractionGeometry,
    phase: Phase,
    orientation: Orientation,
    *,
    max_index: int = 3,
    min_d_spacing_angstrom: float | None = None,
    min_relative_intensity: float = 1e-3,
    max_bands: int | None = None,
    zone_axis_max_index: int = 2,
    provenance: ProvenanceRecord | None = None,
) -> KikuchiPattern:
    """Simulate the Kikuchi pattern of a phase at a known orientation.

    Purpose
    -------
    The forward geometric model of an EBSD or TEM Kikuchi pattern: which lattice
    planes produce visible bands, where their centre lines and edges fall, and
    which zone axes organize them. This is what orientation determination
    inverts, and what a measured pattern is compared against.

    Method
    ------
    Candidate planes are enumerated up to ``max_index``, reduced to one
    representative per antipodal pair (a plane and its opposite normal give the
    same band), and filtered by the lattice-centring reflection condition and by
    a kinematic ``|F|^2`` intensity threshold. Each surviving plane normal is
    carried into the laboratory frame through the crystal orientation and the
    specimen-to-laboratory transform, and its Bragg angle is obtained from
    ``sin(theta_B) = lambda / 2d`` with the relativistically corrected electron
    wavelength. Planes whose spacing cannot satisfy the Bragg condition are
    dropped.

    Parameters
    ----------
    geometry : DiffractionGeometry
        Detector, beam energy, and frames. Its specimen frame must match the
        orientation's.
    phase : Phase
        Lattice, symmetry, and — for intensities — the atomic basis. Without a
        unit cell the intensities degenerate to geometry only.
    orientation : Orientation
        Crystal-to-specimen orientation; its crystal frame and phase are checked
        against ``phase`` rather than assumed compatible.
    max_index : int
        Largest absolute Miller index enumerated. Raising it admits narrower,
        higher-order bands at cubic cost.
    min_d_spacing_angstrom : float, optional
        Drop planes below this spacing. Since band width grows as the spacing
        falls, this excludes the *widest* bands -- the weak high-order ones that
        clutter a pattern -- rather than the narrowest.
    min_relative_intensity : float
        Drop bands weaker than this fraction of the strongest.
    max_bands : int, optional
        Keep only the strongest ``max_bands`` bands, for legibility.
    zone_axis_max_index : int
        Largest absolute index of the zone axes enumerated.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KikuchiPattern
        Bands ordered by decreasing intensity and zone axes by decreasing band
        count, with a :meth:`KikuchiPattern.describe` prose summary.

    Raises
    ------
    ValueError
        For inconsistent frames or phases, or non-positive ``max_index``.

    See Also
    --------
    GnomonicProjection : The coordinate system the band geometry is expressed in.
    """

    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if zone_axis_max_index < 0:
        raise ValueError("zone_axis_max_index must be non-negative.")
    if not 0.0 <= min_relative_intensity <= 1.0:
        raise ValueError("min_relative_intensity must lie in [0, 1].")
    if max_bands is not None and max_bands <= 0:
        raise ValueError("max_bands must be strictly positive when provided.")
    if min_d_spacing_angstrom is not None and min_d_spacing_angstrom <= 0.0:
        raise ValueError("min_d_spacing_angstrom must be strictly positive when provided.")
    if orientation.phase is not None and orientation.phase != phase:
        raise ValueError("orientation.phase must match phase when specified.")
    if orientation.crystal_frame != phase.crystal_frame:
        raise ValueError("orientation.crystal_frame must match phase.crystal_frame.")
    if orientation.specimen_frame != geometry.specimen_frame:
        raise ValueError("orientation.specimen_frame must match geometry.specimen_frame.")

    wavelength = geometry.electron_wavelength_angstrom
    candidates = _antipodal_reduced_triples(max_index)
    allowed = centering_allowed_mask(candidates, ReflectionCondition.from_phase(phase))
    candidates = candidates[allowed]
    if candidates.shape[0] == 0:
        raise ValueError("No reflections survive the lattice centring condition.")

    reciprocal_basis = phase.lattice.reciprocal_basis().matrix
    g_crystal = candidates.astype(np.float64) @ reciprocal_basis.T
    g_magnitude = np.linalg.norm(g_crystal, axis=1)
    d_spacing = 1.0 / g_magnitude

    # A plane diffracts only if lambda / 2d <= 1; wide-spacing planes always do,
    # so in practice this drops only the very finest spacings.
    sin_bragg = wavelength / (2.0 * d_spacing)
    diffracting = sin_bragg < 1.0
    if min_d_spacing_angstrom is not None:
        diffracting &= d_spacing >= min_d_spacing_angstrom
    if not np.any(diffracting):
        raise ValueError("No lattice plane satisfies the Bragg condition at this wavelength.")
    candidates = candidates[diffracting]
    g_crystal = g_crystal[diffracting]
    g_magnitude = g_magnitude[diffracting]
    d_spacing = d_spacing[diffracting]
    sin_bragg = sin_bragg[diffracting]

    structure_factors = electron_structure_factors(
        phase, candidates.astype(np.float64), g_magnitude
    )
    intensities = np.abs(structure_factors) ** 2
    peak = float(np.max(intensities))
    if peak <= 0.0:
        raise ValueError("All candidate reflections have zero structure factor.")
    intensities = intensities / peak
    strong = intensities >= min_relative_intensity
    if not np.any(strong):
        raise ValueError("No reflection exceeds min_relative_intensity.")
    candidates = candidates[strong]
    g_crystal = g_crystal[strong]
    d_spacing = d_spacing[strong]
    sin_bragg = sin_bragg[strong]
    intensities = intensities[strong]

    normals_crystal = g_crystal / np.linalg.norm(g_crystal, axis=1)[:, None]
    normals_specimen = normals_crystal @ orientation.rotation.as_matrix().T
    normals_lab = np.asarray(geometry.specimen_vectors_to_lab(normals_specimen), dtype=np.float64)

    order = np.argsort(-intensities, kind="stable")
    if max_bands is not None:
        order = order[:max_bands]
    bands = tuple(
        KikuchiBand(
            plane=MillerPlane(indices=candidates[index], phase=phase),
            plane_normal_lab=normals_lab[index],
            bragg_angle_rad=float(np.arcsin(sin_bragg[index])),
            d_spacing_angstrom=float(d_spacing[index]),
            intensity=float(intensities[index]),
        )
        for index in order
    )

    zone_axes = _zone_axes_for_bands(
        geometry=geometry,
        phase=phase,
        orientation=orientation,
        bands=bands,
        max_index=zone_axis_max_index,
    )
    return KikuchiPattern(
        geometry=geometry,
        phase=phase,
        orientation=orientation,
        bands=bands,
        zone_axes=zone_axes,
        wavelength_angstrom=wavelength,
        provenance=provenance,
    )


def _zone_axes_for_bands(
    *,
    geometry: DiffractionGeometry,
    phase: Phase,
    orientation: Orientation,
    bands: tuple[KikuchiBand, ...],
    max_index: int,
) -> tuple[KikuchiZoneAxis, ...]:
    """Zone axes of the simulated bands, ordered by how many bands share them.

    A direction is a zone axis of a band when it lies in the band's plane, i.e.
    when the zone law ``hu + kv + lw = 0`` holds. Only axes shared by at least
    two bands are reported, since a single band does not define a hub.
    """

    if max_index == 0 or not bands:
        return ()
    axes = _antipodal_reduced_triples(max_index)
    plane_indices = np.stack([band.plane.indices for band in bands], axis=0).astype(np.int64)
    zone_values = plane_indices @ axes.T
    band_counts = np.count_nonzero(zone_values == 0, axis=0)
    shared = np.flatnonzero(band_counts >= 2)
    if shared.size == 0:
        return ()

    direct_basis = phase.lattice.direct_basis().matrix
    directions_crystal = axes[shared].astype(np.float64) @ direct_basis.T
    directions_crystal = normalize_vectors(directions_crystal)
    directions_specimen = directions_crystal @ orientation.rotation.as_matrix().T
    directions_lab = np.asarray(
        geometry.specimen_vectors_to_lab(directions_specimen), dtype=np.float64
    )
    projection = GnomonicProjection(geometry=geometry)
    coordinates, valid = projection.project_directions(directions_lab)
    on_detector = np.zeros(coordinates.shape[0], dtype=bool)
    if np.any(valid):
        on_detector[valid] = projection.contains(coordinates[valid])

    order = np.argsort(-band_counts[shared], kind="stable")
    return tuple(
        KikuchiZoneAxis(
            indices=axes[shared][index],
            direction_lab=directions_lab[index],
            coordinates=coordinates[index],
            band_count=int(band_counts[shared][index]),
            on_detector=bool(on_detector[index]),
        )
        for index in order
        if valid[index]
    )


__all__ = [
    "GnomonicProjection",
    "KikuchiBand",
    "KikuchiPattern",
    "KikuchiZoneAxis",
    "simulate_kikuchi_pattern",
]
