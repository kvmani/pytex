r"""The Kikuchi map: a stereographic atlas of bands and zone axes for tilt planning.

What this is for
----------------
A TEM operator at a zone axis wants to get to a different one. The pattern on the
screen shows a few bands and a few intersections; the crystal sphere holds the
whole network. A **Kikuchi map** is that network drawn once, stereographically, so
the route from where you are to where you want to be can be read off before the
stage is touched — which band to follow, how far to tilt, and which zone axes you
will pass through on the way.

That is the classical Kikuchi map of Levine, Bell and Thomas (1966) and of
Edington's monograph, where such maps were montaged by hand from dozens of
exposures. Here the same object is computed from the lattice.

The geometry, exactly
---------------------
Incoherently scattered electrons travel in every direction inside the foil. Those
leaving at the Bragg angle :math:`\theta_B` to a lattice plane diffract, so the
diffracting directions for one plane form two cones of semi-angle
:math:`90^\circ - \theta_B` about the plane normal :math:`\mathbf{g}`. On the unit
sphere of directions this is a **band**:

- its **centre line** is the great circle perpendicular to :math:`\mathbf{g}` —
  the trace of the plane itself;
- its **edges** are the two small circles at :math:`90^\circ \mp \theta_B` from
  :math:`\mathbf{g}`;
- its **angular width** is exactly :math:`2\theta_B`, with
  :math:`\sin\theta_B = \lambda / 2d`.

Because the width is :math:`2\arcsin(\lambda/2d) \approx \lambda/d`, a band is
*wider* the *smaller* the interplanar spacing. The widest bands on a map therefore
come from high-index planes, while the strongest come from low-index ones — two
different orderings that are easy to conflate.

Zone axes are where band centre lines cross. A direction :math:`[uvw]` lies on the
centre line of :math:`(hkl)` exactly when :math:`hu + kv + lw = 0`, the Weiss zone
law, so the zone axis at an intersection is :math:`\mathbf{u} \propto \mathbf{g}_1
\times \mathbf{g}_2` and the number of bands crossing there is the number of
reflections satisfying the zone law. That count is the practical measure of how
prominent an axis is on the screen: a four-band intersection is unmistakable, a
two-band one is a guess.

Why the projection is stereographic
-----------------------------------
The map has to cover a hemisphere, which rules out the gnomonic projection used by
:mod:`pytex.diffraction.kikuchi` for a physical detector — there a direction
90 degrees from the pattern centre is at infinity. The stereographic projection is
conformal and bounded: it maps the hemisphere to a disc, preserves angles, and
sends great circles to circular arcs. Angles read off the map are the angles the
stage must turn through, which is the whole point.

Routing
-------
Two zone axes lie on a common band exactly when both are perpendicular to one
reflection — equivalently, when the plane they span is a rational lattice plane.
Following that band is what an experienced operator does by eye, and it is also
the geodesic, since the shortest arc between two directions lies in the plane they
span. When no single band connects them, the map is searched for a route through
intermediate zone axes, which is a shortest-path problem on the graph whose nodes
are zone axes and whose edges are shared bands.

Limits
------
Geometric and kinematic, like the module it sits beside. Band positions and widths
are exact for the stated wavelength; intensities are a kinematic
:math:`|F_{\mathbf{g}}|^2` proxy, so they order the bands sensibly but do not
predict the excess-deficiency asymmetry across a band, dynamical contrast, or
higher-order Laue zone effects. The map is a map of the crystal, not a simulated
image: it carries no foil thickness, no absorption, and no projected potential.

See Also
--------
pytex.diffraction.kikuchi : Kikuchi bands on a physical detector, in gnomonic
    coordinates, for one known orientation.
pytex.tem.navigation : Solving stage tilts for a target zone axis, given a
    calibrated holder.
pytex.tem.path : Validating the trajectory the route implies against the holder
    envelope.
"""

from __future__ import annotations

import heapq
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import FloatArray, freeze_array, normalize_vector
from pytex.core.hexagonal import direction_uvw_to_uvtw, is_hexagonal_phase, plane_hkl_to_hkil
from pytex.core.lattice import Phase
from pytex.core.miller import MillerPlane
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.kinematic import centering_allowed_mask, electron_wavelength_angstrom
from pytex.diffraction.physics import ReflectionCondition
from pytex.diffraction.scattering import electron_structure_factor_angstrom
from pytex.diffraction.stereonets import sample_great_circle, sample_small_circle
from pytex.texture.projections import project_directions

__all__ = [
    "DEFAULT_ROUTE_MAX_LEG_DEG",
    "KikuchiMapBand",
    "KikuchiMapZoneAxis",
    "KikuchiRoute",
    "KikuchiRouteLeg",
    "StereographicKikuchiMap",
    "compute_kikuchi_map",
    "plan_kikuchi_route",
]

#: Largest tilt accepted for a single leg of a route, in degrees.
#:
#: A long excursion is fragile for a reason that is structural rather than
#: numerical: an error ``dphi`` in the initial orientation costs a residual of
#: about ``dphi sin(theta)`` after a hop of angular length ``theta``, and the
#: operator loses the band before arriving. Routing through intermediate zone
#: axes and re-indexing at each converts an open-loop tilt into a self-correcting
#: procedure. Thirty degrees is a comfortable single view-field excursion on a
#: double-tilt holder.
DEFAULT_ROUTE_MAX_LEG_DEG = 30.0

#: Angular tolerance for deciding that a direction lies on a band centre line.
#:
#: Exact zone-law membership is an integer condition, so this only absorbs
#: floating-point error in the cross products and normalizations that produce the
#: candidate axes.
_ZONE_LAW_ATOL_DEG = 0.05

def _format_direction(phase: Phase, indices: Sequence[int]) -> str:
    """Direction indices in the notation the phase's crystal system uses.

    A hexagonal direction is written in four-index Miller-Bravais form, because
    that is what the literature writes and because the three-index form hides the
    symmetry: the members of a family do not look like permutations of each other
    in ``[uvw]``. Everything else keeps three indices.
    """

    values = [int(value) for value in indices]
    if is_hexagonal_phase(phase):
        values = [int(value) for value in direction_uvw_to_uvtw(values)]
    return format_direction_indices(values, style="plain")


def _format_plane(phase: Phase, indices: Sequence[int]) -> str:
    """Plane indices in the notation the phase's crystal system uses.

    Four-index Bravais-Miller ``(hkil)`` for hexagonal phases, where ``i = -(h+k)``
    is redundant but conventional; three indices otherwise. Miller indices carry no
    star: they are already components in the reciprocal basis.
    """

    values = [int(value) for value in indices]
    if is_hexagonal_phase(phase):
        values = [int(value) for value in plane_hkl_to_hkil(values)]
    return format_plane_indices(values, style="plain")


def _antipodal_reduced_triples(max_index: int) -> np.ndarray:
    """Integer triples up to ``max_index``, one representative per antipodal pair.

    A plane and its opposite normal give the same band, so only half the grid is
    needed. The retained representative is the lexicographically positive one.
    """

    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1)
    triples = grid.reshape(-1, 3)
    triples = triples[np.any(triples != 0, axis=1)]
    first_nonzero = triples[np.arange(triples.shape[0]), np.argmax(triples != 0, axis=1)]
    return np.ascontiguousarray(triples[first_nonzero > 0])


def _view_matrix(phase: Phase, centre: ArrayLike, horizontal: ArrayLike) -> FloatArray:
    """Rotation carrying crystal Cartesian coordinates into the map frame.

    The map frame has the centre direction along ``+z`` — the projection pole,
    which is the direction the beam runs along when the operator is at the centre
    of the map — and the horizontal reference along ``+x``. The horizontal is
    orthogonalized against the centre rather than required to be perpendicular,
    so a caller may name any convenient direction.
    """

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=np.float64)
    axis_z = normalize_vector(np.asarray(centre, dtype=np.float64) @ direct.T)
    horizontal_cartesian = np.asarray(horizontal, dtype=np.float64) @ direct.T
    residual = horizontal_cartesian - float(np.dot(horizontal_cartesian, axis_z)) * axis_z
    if float(np.linalg.norm(residual)) < 1e-9:
        raise ValueError(
            "The horizontal reference direction is parallel to the map centre, so it "
            "cannot orient the map. Choose a direction not along the centre."
        )
    axis_x = normalize_vector(residual)
    axis_y = np.cross(axis_z, axis_x)
    return freeze_array(np.stack([axis_x, axis_y, axis_z], axis=0))


@dataclass(frozen=True, slots=True)
class KikuchiMapBand:
    r"""One band of a Kikuchi map: a lattice plane drawn on the crystal sphere.

    Purpose
    -------
    The road an operator follows. A band is the locus of beam directions that
    satisfy the Bragg condition for one lattice plane, so tilting along it keeps
    that reflection excited — which is why the band, and not a straight line in
    stage coordinates, is the natural path between two zone axes.

    When to use
    -----------
    Obtained from :func:`compute_kikuchi_map`; use :meth:`centre_trace` and
    :meth:`edge_traces` to draw it, :attr:`angular_width_deg` to predict how
    visible it will be, and :meth:`contains_direction` to test whether a beam
    direction is inside it.

    Attributes
    ----------
    plane : MillerPlane
        The diffracting plane, carrying its phase.
    normal_map : np.ndarray
        Unit plane normal in the map frame, in which ``+z`` is the map centre.
    bragg_angle_deg : float
        :math:`\theta_B`; half the band's angular width.
    d_spacing_angstrom : float
        Interplanar spacing.
    relative_intensity : float
        Kinematic :math:`|F_{\mathbf{g}}|^2` relative to the strongest band, in
        ``[0, 1]``. Orders the bands; does not predict their contrast.
    family_multiplicity : int
        Size of the symmetry family the plane belongs to, counting antipodal pairs
        once. A large family means many parallel-looking bands elsewhere on the
        map.
    """

    plane: MillerPlane
    normal_map: FloatArray
    bragg_angle_deg: float
    d_spacing_angstrom: float
    relative_intensity: float
    family_multiplicity: int

    def __post_init__(self) -> None:
        normal = np.asarray(self.normal_map, dtype=np.float64).reshape(3)
        norm = float(np.linalg.norm(normal))
        if not np.isfinite(norm) or np.isclose(norm, 0.0):
            raise ValueError("KikuchiMapBand.normal_map must be a finite non-zero vector.")
        object.__setattr__(self, "normal_map", freeze_array(normal / norm))
        if not 0.0 < self.bragg_angle_deg < 90.0:
            raise ValueError("KikuchiMapBand.bragg_angle_deg must lie in (0, 90).")
        if not np.isfinite(self.d_spacing_angstrom) or self.d_spacing_angstrom <= 0.0:
            raise ValueError("KikuchiMapBand.d_spacing_angstrom must be finite and positive.")
        if not 0.0 <= self.relative_intensity <= 1.0:
            raise ValueError("KikuchiMapBand.relative_intensity must lie in [0, 1].")
        if self.family_multiplicity <= 0:
            raise ValueError("KikuchiMapBand.family_multiplicity must be strictly positive.")

    @property
    def indices(self) -> tuple[int, ...]:
        """The plane's Miller indices as a plain tuple."""

        return tuple(int(value) for value in np.asarray(self.plane.indices).ravel())

    @property
    def angular_width_deg(self) -> float:
        r"""Full angular width :math:`2\theta_B`, in degrees.

        This is what a band-width measurement yields, and it fixes the spacing
        through Bragg's law. It *grows* as the spacing falls, so the widest bands
        on a map are the high-index ones.
        """

        return float(2.0 * self.bragg_angle_deg)

    def centre_trace(self, *, method: str = "stereographic", samples: int = 361) -> FloatArray:
        """The band's centre line, projected onto the map plane.

        The trace of the lattice plane: the great circle perpendicular to the
        plane normal. Returned as ``(samples, 2)`` plane coordinates.
        """

        return _project_trace(
            sample_great_circle(self.normal_map, samples=samples, half_circle=False),
            method=method,
        )

    def edge_traces(
        self, *, method: str = "stereographic", samples: int = 361
    ) -> tuple[FloatArray, FloatArray]:
        r"""The two band edges, projected onto the map plane.

        The Kossel cones at :math:`90^\circ \mp \theta_B` from the plane normal.
        These are small circles about the normal, not great circles, which is why
        a band has a measurable width at all.
        """

        return (
            _project_trace(
                _small_circle_about(self.normal_map, 90.0 - self.bragg_angle_deg, samples),
                method=method,
            ),
            _project_trace(
                _small_circle_about(self.normal_map, 90.0 + self.bragg_angle_deg, samples),
                method=method,
            ),
        )

    def contains_direction(self, direction: ArrayLike, *, atol_deg: float = 0.0) -> bool:
        """Whether a beam direction lies inside the band.

        True when the direction is within :math:`\\theta_B` of the plane's trace,
        i.e. when that plane's reflection is excited for that beam direction.
        """

        unit = normalize_vector(np.asarray(direction, dtype=np.float64))
        offset_deg = 90.0 - float(
            np.degrees(np.arccos(np.clip(abs(float(np.dot(unit, self.normal_map))), -1.0, 1.0)))
        )
        return bool(abs(offset_deg) <= self.bragg_angle_deg + atol_deg)

    def describe(self) -> str:
        """One-sentence prose summary of the band, with its conventions stated."""

        notation = _format_plane(self.plane.phase, self.indices)
        return (
            f"Kikuchi band {notation} of {self.plane.phase.name}: interplanar spacing "
            f"{self.d_spacing_angstrom:.4f} A, Bragg angle {self.bragg_angle_deg:.4f} deg, so "
            f"angular width {self.angular_width_deg:.4f} deg. Relative kinematic intensity "
            f"{self.relative_intensity:.4f}; the plane belongs to a family of "
            f"{self.family_multiplicity} equivalent orientations. The centre line is the trace "
            f"of the plane and the edges are the Kossel cones at 90 deg -/+ the Bragg angle "
            f"from the plane normal."
        )


@dataclass(frozen=True, slots=True)
class KikuchiMapZoneAxis:
    """One zone axis of a Kikuchi map: a crossing point of band centre lines.

    Purpose
    -------
    The destinations. A zone axis is a beam direction along which many
    reflections are simultaneously excited, so it gives the symmetric spot
    pattern used for orientation determination, and it is where an operator aims.

    Attributes
    ----------
    phase : Phase
        Carried so the axis can name itself in its crystal system's notation.
    indices : tuple of int
        The direction indices ``[uvw]``, always three-index, whatever notation
        :meth:`describe` prints.
    direction_map : np.ndarray
        Unit direction in the map frame.
    band_indices : tuple of int
        Positions, in the parent map's ``bands``, of the bands whose centre lines
        pass through this axis.
    polar_angle_deg : float
        Angle from the map centre. Zero at the centre of the projection.
    """

    phase: Phase
    indices: tuple[int, ...]
    direction_map: FloatArray
    band_indices: tuple[int, ...]
    polar_angle_deg: float

    def __post_init__(self) -> None:
        direction = np.asarray(self.direction_map, dtype=np.float64).reshape(3)
        norm = float(np.linalg.norm(direction))
        if not np.isfinite(norm) or np.isclose(norm, 0.0):
            raise ValueError("KikuchiMapZoneAxis.direction_map must be finite and non-zero.")
        object.__setattr__(self, "direction_map", freeze_array(direction / norm))
        if len(self.indices) != 3:
            raise ValueError("KikuchiMapZoneAxis.indices must have three components.")
        if len(self.band_indices) < 2:
            raise ValueError(
                "A zone axis is an intersection, so at least two bands must pass through it."
            )

    @property
    def order(self) -> int:
        """Number of bands crossing here — how conspicuous the axis is on screen."""

        return len(self.band_indices)

    def projected(self, *, method: str = "stereographic") -> FloatArray:
        """The axis as a single ``(2,)`` point on the map plane."""

        point: FloatArray = _project_trace(self.direction_map[None, :], method=method)[0]
        return point

    def describe(self) -> str:
        """One-sentence prose summary of the zone axis."""

        notation = _format_direction(self.phase, self.indices)
        return (
            f"Zone axis {notation}: {self.order} bands cross here, "
            f"{self.polar_angle_deg:.2f} deg from the map centre."
        )


@dataclass(frozen=True, slots=True)
class KikuchiRouteLeg:
    """One hop of a tilt route: follow one band from one zone axis to the next.

    Attributes
    ----------
    start_indices, end_indices : tuple of int
        The zone axes at each end.
    band_indices : tuple of int or None
        The band followed, as plane indices. ``None`` when the two axes share no
        band, in which case the leg is a bare geodesic and the operator has no
        line to track.
    tilt_deg : float
        Angle between the two zone axes: the total stage travel for this leg.
    waypoint_indices : tuple of tuple of int
        Zone axes lying on the followed band strictly between the endpoints, in
        the order they are passed. These are the landmarks that confirm the tilt
        is on track.
    start_direction, end_direction : np.ndarray
        Unit directions in the map frame, oriented so consecutive legs of a route
        join up. A zone axis is a *line*, so each end has two equally valid senses
        and a route has to commit to one of them, or the drawn path jumps across
        the projection between legs. The indices above carry the same sense as
        these vectors, which is why they can differ in sign from the map's own
        canonical representative of the same axis -- both name the same axis.
    """

    phase: Phase
    start_indices: tuple[int, ...]
    end_indices: tuple[int, ...]
    band_indices: tuple[int, ...] | None
    tilt_deg: float
    waypoint_indices: tuple[tuple[int, ...], ...]
    start_direction: FloatArray
    end_direction: FloatArray

    def describe(self) -> str:
        """Imperative instruction for this leg, in the operator's language."""

        start = _format_direction(self.phase, self.start_indices)
        end = _format_direction(self.phase, self.end_indices)
        if self.band_indices is None:
            return (
                f"From {start}, tilt {self.tilt_deg:.2f} deg to {end}. No single Kikuchi band "
                f"joins these axes, so there is no line to follow: tilt on the calculated "
                f"angles and re-index on arrival."
            )
        band = _format_plane(self.phase, self.band_indices)
        if not self.waypoint_indices:
            passing = "no intermediate zone axis lies on the way"
        else:
            names = ", ".join(
                _format_direction(self.phase, indices) for indices in self.waypoint_indices
            )
            passing = f"passing {names}"
        return (
            f"From {start}, follow the {band} Kikuchi band {self.tilt_deg:.2f} deg to {end}, "
            f"{passing}."
        )


@dataclass(frozen=True, slots=True)
class KikuchiRoute:
    """A planned tilt route across a Kikuchi map, as a sequence of band-following legs.

    Purpose
    -------
    The deliverable of the whole module: an instruction an operator can act on,
    with the band to follow, the angle to turn, and the landmarks to expect.

    Attributes
    ----------
    legs : tuple of KikuchiRouteLeg
        In order. Empty when the start and target are the same axis.
    reachable : bool
        Whether a route was found at all within the map's zone-axis network and
        the requested maximum leg.
    """

    phase: Phase
    legs: tuple[KikuchiRouteLeg, ...]
    reachable: bool
    start_indices: tuple[int, ...]
    target_indices: tuple[int, ...]
    start_direction: FloatArray
    target_direction: FloatArray
    provenance: ProvenanceRecord | None = None

    @property
    def total_tilt_deg(self) -> float:
        """Sum of the leg angles: the total stage travel."""

        return float(sum(leg.tilt_deg for leg in self.legs))

    @property
    def direct_tilt_deg(self) -> float:
        """Angle between start and target: the travel an unbroken hop would need.

        Compare with :attr:`total_tilt_deg`. Their ratio is what multi-hop
        routing costs in travel, in exchange for being able to follow a band and
        re-index on the way.
        """

        if not self.legs:
            return 0.0
        cosine = abs(float(np.dot(self.start_direction, self.target_direction)))
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    @property
    def hop_count(self) -> int:
        """Number of legs."""

        return len(self.legs)

    def describe(self) -> str:
        """Prose instructions for the whole route, with the conventions stated."""

        start = _format_direction(self.phase, self.start_indices)
        target = _format_direction(self.phase, self.target_indices)
        if not self.reachable:
            return (
                f"No route from {start} to {target} was found on this Kikuchi map. Either the "
                f"target is not a zone axis of the map, or every path needs a leg longer than "
                f"the requested maximum. Raising the map's zone-axis index bound adds "
                f"intermediate axes to route through."
            )
        if not self.legs:
            return f"Already at {start}; no tilt required."
        lines = [
            f"Route from {start} to {target}: {self.hop_count} "
            f"{'leg' if self.hop_count == 1 else 'legs'}, {self.total_tilt_deg:.2f} deg of total "
            f"travel against a direct separation of {self.direct_tilt_deg:.2f} deg."
        ]
        lines.extend(f"  {index + 1}. {leg.describe()}" for index, leg in enumerate(self.legs))
        lines.append(
            "Angles are between crystal directions, so they are the stage travel only for a "
            "holder whose tilt axes are calibrated; pass the route to "
            "pytex.tem.navigation to convert each hop into alpha and beta and to check it "
            "against the holder envelope."
        )
        return "\n".join(lines)

    def to_json_dict(self) -> dict[str, Any]:
        """JSON-ready payload, in lockstep with :meth:`describe`."""

        return {
            "schema": "pytex.kikuchi_route.v1",
            "start": list(self.start_indices),
            "target": list(self.target_indices),
            "reachable": self.reachable,
            "hop_count": self.hop_count,
            "total_tilt_deg": self.total_tilt_deg,
            "direct_tilt_deg": self.direct_tilt_deg,
            "legs": [
                {
                    "start": list(leg.start_indices),
                    "end": list(leg.end_indices),
                    "band": None if leg.band_indices is None else list(leg.band_indices),
                    "tilt_deg": leg.tilt_deg,
                    "waypoints": [list(indices) for indices in leg.waypoint_indices],
                }
                for leg in self.legs
            ],
        }


@dataclass(frozen=True, slots=True)
class StereographicKikuchiMap:
    """The Kikuchi band network of a phase, on a stereographic projection.

    Purpose
    -------
    The operator's road atlas: every band and every zone axis of one phase, in one
    picture, with the routing between them.

    When to use
    -----------
    Build it with :func:`compute_kikuchi_map` once per phase and beam energy, then
    query it repeatedly — :meth:`bands_through` for what is visible at an axis,
    :meth:`route_to` for how to get from one axis to another, and
    :func:`pytex.plotting.plot_kikuchi_map` to draw it.

    Attributes
    ----------
    phase : Phase
    beam_energy_kev : float
    wavelength_angstrom : float
        Relativistically corrected electron wavelength; the band widths depend on
        it.
    centre_indices, horizontal_indices : tuple of int
        The direction at the centre of the projection and the direction drawn
        along ``+x``. These fix the orientation of the map, and every map frame
        vector is relative to them.
    view_matrix : np.ndarray
        ``(3, 3)`` rotation from crystal Cartesian to map coordinates.
    bands : tuple of KikuchiMapBand
        Ordered by decreasing kinematic intensity.
    zone_axes : tuple of KikuchiMapZoneAxis
        Ordered by decreasing band count, then by proximity to the map centre.
    has_intensity_model : bool
        Whether the phase carried the atomic basis the structure factors need.
        ``False`` means every band has relative intensity 1 and the band ordering
        is enumeration order, not prominence -- the geometry is unaffected.
    """

    phase: Phase
    beam_energy_kev: float
    wavelength_angstrom: float
    centre_indices: tuple[int, ...]
    horizontal_indices: tuple[int, ...]
    view_matrix: FloatArray
    bands: tuple[KikuchiMapBand, ...]
    zone_axes: tuple[KikuchiMapZoneAxis, ...]
    has_intensity_model: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        matrix = np.asarray(self.view_matrix, dtype=np.float64).reshape(3, 3)
        if not np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-9):
            raise ValueError("StereographicKikuchiMap.view_matrix must be orthonormal.")
        object.__setattr__(self, "view_matrix", freeze_array(matrix))
        if not self.bands:
            raise ValueError("A Kikuchi map needs at least one band.")

    @property
    def band_count(self) -> int:
        """Number of bands on the map."""

        return len(self.bands)

    @property
    def zone_axis_count(self) -> int:
        """Number of zone axes on the map."""

        return len(self.zone_axes)

    def band_for_plane(self, indices: ArrayLike) -> KikuchiMapBand | None:
        """The band of a given plane, or ``None`` if that plane is not on the map.

        Matches a plane and its opposite normal, since both give the same band.
        """

        wanted = np.asarray(indices, dtype=np.int64).reshape(3)
        for band in self.bands:
            candidate = np.asarray(band.indices, dtype=np.int64)
            if np.array_equal(candidate, wanted) or np.array_equal(candidate, -wanted):
                return band
        return None

    def zone_axis_for_direction(self, indices: ArrayLike) -> KikuchiMapZoneAxis | None:
        """The zone axis with these indices, or ``None``. Matches both senses."""

        wanted = np.asarray(indices, dtype=np.int64).reshape(3)
        for axis in self.zone_axes:
            candidate = np.asarray(axis.indices, dtype=np.int64)
            if np.array_equal(candidate, wanted) or np.array_equal(candidate, -wanted):
                return axis
        return None

    def bands_through(self, indices: ArrayLike) -> tuple[KikuchiMapBand, ...]:
        """The bands crossing a zone axis, strongest first.

        This is the list of reflections excited at that axis — what the operator
        sees on arriving, and the answer to "which band do I follow from here?".
        """

        axis = self.zone_axis_for_direction(indices)
        if axis is None:
            return ()
        return tuple(self.bands[position] for position in axis.band_indices)

    def shared_bands(self, first: ArrayLike, second: ArrayLike) -> tuple[KikuchiMapBand, ...]:
        """Bands whose centre lines pass through both zone axes, strongest first.

        A non-empty result means the two axes are joined by a band the operator
        can track, and following the strongest of them is the recommended route.
        """

        left = self.zone_axis_for_direction(first)
        right = self.zone_axis_for_direction(second)
        if left is None or right is None:
            return ()
        common = [position for position in left.band_indices if position in right.band_indices]
        return tuple(self.bands[position] for position in common)

    def route_to(
        self,
        start: ArrayLike,
        target: ArrayLike,
        *,
        max_leg_deg: float = DEFAULT_ROUTE_MAX_LEG_DEG,
    ) -> KikuchiRoute:
        """Plan a route between two zone axes of this map.

        Convenience wrapper over :func:`plan_kikuchi_route`.
        """

        return plan_kikuchi_route(self, start, target, max_leg_deg=max_leg_deg)

    def describe(self) -> str:
        """Prose summary of the map, its conventions, and its limits."""

        centre = _format_direction(self.phase, self.centre_indices)
        horizontal = _format_direction(self.phase, self.horizontal_indices)
        widest = max(self.bands, key=lambda band: band.angular_width_deg)
        strongest = self.bands[0]
        prominent = [axis for axis in self.zone_axes if axis.order >= 3]
        lines = [
            f"Stereographic Kikuchi map of {self.phase.name} at {self.beam_energy_kev:.1f} kV "
            f"(wavelength {self.wavelength_angstrom:.5f} A), centred on {centre} with "
            f"{horizontal} along +x.",
            f"{self.band_count} bands and {self.zone_axis_count} zone axes, "
            f"{len(prominent)} of them crossed by three or more bands.",
            f"The strongest band is "
            f"{_format_plane(self.phase, strongest.indices)} at "
            f"{strongest.angular_width_deg:.3f} deg wide; the widest is "
            f"{_format_plane(self.phase, widest.indices)} at "
            f"{widest.angular_width_deg:.3f} deg. Width is 2 arcsin(lambda / 2d), so the widest "
            f"band is the one with the smallest spacing, not the strongest.",
            "Band centre lines are plane traces (great circles) and band edges are the Kossel "
            "cones at 90 deg -/+ the Bragg angle from the plane normal, sampled exactly rather "
            "than approximated as straight.",
            "Intensities are a kinematic |F|^2 proxy: they order the bands but do not predict "
            "excess-deficiency asymmetry, dynamical contrast, or higher-order Laue zone effects."
            if self.has_intensity_model
            else "This phase carries no atomic basis, so there are no structure factors: every "
            "band is reported at relative intensity 1 and the ordering is enumeration order, "
            "not prominence. Every geometric quantity -- traces, widths, zone axes, routes -- is "
            "unaffected.",
        ]
        return "\n".join(lines)

    def to_json_dict(self) -> dict[str, Any]:
        """JSON-ready payload, in lockstep with :meth:`describe`."""

        return {
            "schema": "pytex.stereographic_kikuchi_map.v1",
            "phase": self.phase.name,
            "beam_energy_kev": self.beam_energy_kev,
            "wavelength_angstrom": self.wavelength_angstrom,
            "centre": list(self.centre_indices),
            "horizontal": list(self.horizontal_indices),
            "has_intensity_model": self.has_intensity_model,
            "bands": [
                {
                    "indices": list(band.indices),
                    "d_spacing_angstrom": band.d_spacing_angstrom,
                    "bragg_angle_deg": band.bragg_angle_deg,
                    "angular_width_deg": band.angular_width_deg,
                    "relative_intensity": band.relative_intensity,
                    "family_multiplicity": band.family_multiplicity,
                }
                for band in self.bands
            ],
            "zone_axes": [
                {
                    "indices": list(axis.indices),
                    "order": axis.order,
                    "polar_angle_deg": axis.polar_angle_deg,
                    "bands": [list(self.bands[position].indices) for position in axis.band_indices],
                }
                for axis in self.zone_axes
            ],
        }


def _project_trace(directions: ArrayLike, *, method: str) -> FloatArray:
    """Project unit directions onto the map plane, folding to one hemisphere."""

    projected = project_directions(np.asarray(directions, dtype=np.float64), method=method,
                                   antipodal=True)
    return freeze_array(np.asarray(projected, dtype=np.float64))


def _small_circle_about(pole: ArrayLike, polar_deg: float, samples: int) -> FloatArray:
    """A small circle at a fixed angle to an arbitrary pole.

    :func:`pytex.diffraction.stereonets.sample_small_circle` builds the circle
    about ``+z`` and accepts only polar angles up to 90 degrees, since a stereonet
    graticule needs no more. A band edge at ``90 + theta_B`` is the reflection of
    the edge at ``90 - theta_B`` through the plane, so it is generated from the
    complementary angle and the opposite pole.
    """

    axis = normalize_vector(np.asarray(pole, dtype=np.float64))
    angle = float(polar_deg)
    if angle > 90.0:
        axis = -axis
        angle = 180.0 - angle
    circle = sample_small_circle(angle, samples=samples)
    reference = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if np.isclose(abs(float(np.dot(axis, reference))), 1.0, atol=1e-12):
        rotation = np.eye(3) if float(axis[2]) > 0.0 else np.diag([1.0, -1.0, -1.0])
    else:
        # Rodrigues rotation taking +z onto the pole.
        cross = np.cross(reference, axis)
        sine = float(np.linalg.norm(cross))
        cosine = float(np.dot(reference, axis))
        unit = cross / sine
        skew = np.array(
            [
                [0.0, -unit[2], unit[1]],
                [unit[2], 0.0, -unit[0]],
                [-unit[1], unit[0], 0.0],
            ]
        )
        rotation = np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)
    return freeze_array(np.asarray(circle @ rotation.T, dtype=np.float64))


def compute_kikuchi_map(
    phase: Phase,
    *,
    beam_energy_kev: float = 200.0,
    centre_direction: ArrayLike = (0, 0, 1),
    horizontal_direction: ArrayLike = (1, 0, 0),
    max_index: int = 4,
    min_d_spacing_angstrom: float | None = None,
    min_relative_intensity: float = 5e-3,
    max_bands: int | None = None,
    include_higher_orders: bool = False,
    zone_axis_max_index: int = 3,
    min_zone_axis_order: int = 2,
    max_polar_angle_deg: float = 90.0,
    provenance: ProvenanceRecord | None = None,
) -> StereographicKikuchiMap:
    r"""Compute the stereographic Kikuchi map of a phase.

    Purpose
    -------
    Builds the band-and-zone-axis network a TEM operator navigates by: which
    lattice planes give visible bands, where their traces run on the crystal
    sphere, which zone axes they cross at, and hence which routes exist between
    those axes.

    When to use
    -----------
    Once per phase and accelerating voltage, before an experiment, to plan a tilt
    series; or to produce the reference atlas a measured pattern is compared
    against. For the pattern on a *physical detector* at a *known* orientation,
    use :func:`pytex.diffraction.kikuchi.simulate_kikuchi_pattern` instead — that
    is a different question, and gnomonic coordinates are the right ones for it.

    Method
    ------
    1. Enumerate integer triples to ``max_index``, one per antipodal pair, and
       drop those forbidden by the lattice centring.
    2. Compute :math:`\mathbf{g}` in the crystal Cartesian frame, hence
       :math:`d = 1/|\mathbf{g}|` and :math:`\sin\theta_B = \lambda / 2d` with the
       relativistically corrected wavelength.
    3. Keep reflections above ``min_relative_intensity`` on a kinematic
       :math:`|F_{\mathbf{g}}|^2` scale, and above ``min_d_spacing_angstrom`` if
       given.
    4. Rotate the normals into the map frame, whose ``+z`` is
       ``centre_direction``.
    5. Cross every pair of band normals to get candidate zone axes, rationalize
       each to a low-index :math:`[uvw]`, merge duplicates, and record every band
       satisfying the zone law :math:`hu + kv + lw = 0` for it.
    6. Keep axes with at least ``min_zone_axis_order`` bands and within
       ``max_polar_angle_deg`` of the centre.

    Parameters
    ----------
    phase : Phase
        Lattice, symmetry, and -- for intensities -- the atomic basis. A phase
        carrying only a lattice still gives the full geometry: traces, widths, zone
        axes and routes are all determined without a structure factor. Only the
        band *ordering* is lost, and
        :attr:`StereographicKikuchiMap.has_intensity_model` records that it was.
    beam_energy_kev : float
        Accelerating voltage. It enters only through the wavelength, so it scales
        every band width and changes nothing else about the map's topology.
    centre_direction : ArrayLike
        Lattice direction at the centre of the projection, in direct-basis
        indices. The classical standard projection of a cubic crystal is centred
        on ``[001]``.
    horizontal_direction : ArrayLike
        Lattice direction drawn along ``+x``. Orthogonalized against the centre,
        so it need not be exactly perpendicular.
    max_index : int
        Largest absolute Miller index enumerated for bands. Raising it admits
        wider, weaker, higher-order bands at cubic cost.
    min_d_spacing_angstrom : float, optional
        Drop planes with a smaller spacing. Since width grows as spacing falls,
        this excludes the *widest* bands, which for a map is usually what is
        wanted: they are the weak high-index ones that clutter it.
    min_relative_intensity : float
        Drop bands weaker than this fraction of the strongest.
    max_bands : int, optional
        Keep only the strongest this many bands, for a legible map.
    include_higher_orders : bool
        Keep every order of each reflection as its own band. Off by default,
        because all orders of one reflection share a single centre line -- (222)
        traces the same great circle as (111) -- so keeping them draws coincident
        lines and multiplies every zone axis's band count. Turn it on to study
        higher-order line positions, where the differing Bragg angles matter.
    zone_axis_max_index : int
        Largest absolute index accepted when rationalizing a zone axis. Axes
        needing higher indices are not lattice landmarks an operator would use.
    min_zone_axis_order : int
        Smallest number of crossing bands for an axis to be kept. Two is the
        definition of an intersection; three or four gives only the conspicuous
        ones.
    max_polar_angle_deg : float
        Discard axes further than this from the map centre. Ninety degrees keeps
        the whole hemisphere.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    StereographicKikuchiMap
        Bands by decreasing intensity, zone axes by decreasing band count.

    Raises
    ------
    ValueError
        For non-positive bounds, an out-of-range intensity threshold, a
        horizontal direction parallel to the centre, or a reflection set that
        empties under the filters.

    Examples
    --------
    The zone axes of a cubic map obey the zone law exactly, which is the check
    worth making: every band recorded at ``[001]`` has ``l = 0``.

    See Also
    --------
    plan_kikuchi_route : Routing between two zone axes of the returned map.
    pytex.plotting.plot_kikuchi_map : Drawing it.
    """

    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if zone_axis_max_index <= 0:
        raise ValueError("zone_axis_max_index must be strictly positive.")
    if not 0.0 <= min_relative_intensity <= 1.0:
        raise ValueError("min_relative_intensity must lie in [0, 1].")
    if max_bands is not None and max_bands <= 0:
        raise ValueError("max_bands must be strictly positive when provided.")
    if min_d_spacing_angstrom is not None and min_d_spacing_angstrom <= 0.0:
        raise ValueError("min_d_spacing_angstrom must be strictly positive when provided.")
    if min_zone_axis_order < 2:
        raise ValueError("min_zone_axis_order must be at least 2: a crossing needs two bands.")
    if not 0.0 < float(max_polar_angle_deg) <= 90.0:
        raise ValueError("max_polar_angle_deg must lie in (0, 90].")

    wavelength = electron_wavelength_angstrom(beam_energy_kev)
    view = _view_matrix(phase, centre_direction, horizontal_direction)

    candidates = _antipodal_reduced_triples(max_index)
    allowed = centering_allowed_mask(candidates, ReflectionCondition.from_phase(phase))
    candidates = candidates[allowed]
    if candidates.shape[0] == 0:
        raise ValueError("No reflections survive the lattice centring condition.")

    reciprocal = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=np.float64)
    g_crystal = candidates.astype(np.float64) @ reciprocal.T
    g_magnitude = np.linalg.norm(g_crystal, axis=1)
    d_spacing = 1.0 / g_magnitude
    sine_bragg = wavelength / (2.0 * d_spacing)

    keep = sine_bragg < 1.0
    if min_d_spacing_angstrom is not None:
        keep &= d_spacing >= float(min_d_spacing_angstrom)
    if not np.any(keep):
        raise ValueError("No lattice plane satisfies the Bragg condition under these bounds.")
    candidates, g_crystal = candidates[keep], g_crystal[keep]
    g_magnitude, d_spacing, sine_bragg = g_magnitude[keep], d_spacing[keep], sine_bragg[keep]

    if not include_higher_orders:
        # A band is a plane *trace*, and every order of one reflection shares that
        # trace: (222) draws its centre line exactly where (111) does. Keeping all
        # of them would put several coincident lines on the map and would inflate
        # every zone axis's band count, which is the number an operator uses to
        # judge how conspicuous an intersection is. So each family of collinear
        # triples is represented by its lowest allowed order -- the largest
        # spacing, hence the strongest and narrowest band of the set. Note the
        # representative is not always coprime: in an fcc lattice the {100} trace
        # is drawn by (200), because (100) is extinguished.
        divisors = np.gcd.reduce(np.abs(candidates), axis=1)
        reduced = candidates // divisors[:, None]
        keys = [tuple(int(value) for value in row) for row in reduced]
        first_order: dict[tuple[int, ...], int] = {}
        for position, key in enumerate(keys):
            best = first_order.get(key)
            if best is None or divisors[position] < divisors[best]:
                first_order[key] = position
        selected = np.asarray(sorted(first_order.values()), dtype=np.int64)
        candidates, g_crystal = candidates[selected], g_crystal[selected]
        g_magnitude = g_magnitude[selected]
        d_spacing, sine_bragg = d_spacing[selected], sine_bragg[selected]

    # The Mott-Bethe electron structure factor, not the atomic-number proxy in
    # pytex.diffraction.kinematic. The proxy replaces f_e(s) by Z, which is
    # s-independent, so on a fixture without Debye-Waller factors every allowed
    # reflection of a monatomic crystal comes out with the *same* intensity and the
    # band ordering carries no information at all. A map's whole visual grammar is
    # that strong bands are prominent, so the ordering has to be real.
    #
    # A phase carrying only a lattice has no structure factor, and inventing one
    # would be worse than declining to order the bands: the geometry -- traces,
    # widths, zone axes, routes -- is completely determined without it, and that is
    # most of what a map is for. So the intensities fall back to uniform and
    # ``has_intensity_model`` records that they did, which describe() then states.
    has_intensity_model = phase.unit_cell is not None and bool(phase.unit_cell.sites)
    if has_intensity_model:
        structure_factors = electron_structure_factor_angstrom(
            phase, candidates, beam_energy_kev=beam_energy_kev
        )
        intensities = np.abs(structure_factors) ** 2
        peak = float(np.max(intensities))
        if peak <= 0.0:
            raise ValueError("All candidate reflections have zero structure factor.")
        intensities = intensities / peak
    else:
        intensities = np.ones(candidates.shape[0], dtype=np.float64)
    strong = intensities >= min_relative_intensity
    if not np.any(strong):
        raise ValueError("No reflection exceeds min_relative_intensity.")
    candidates, g_crystal = candidates[strong], g_crystal[strong]
    d_spacing, sine_bragg, intensities = d_spacing[strong], sine_bragg[strong], intensities[strong]

    normals_crystal = g_crystal / np.linalg.norm(g_crystal, axis=1)[:, None]
    normals_map = normals_crystal @ view.T

    order = np.argsort(-intensities, kind="stable")
    if max_bands is not None:
        order = order[: int(max_bands)]

    operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
    bands = tuple(
        KikuchiMapBand(
            plane=MillerPlane(indices=candidates[index], phase=phase),
            normal_map=normals_map[index],
            bragg_angle_deg=float(np.degrees(np.arcsin(sine_bragg[index]))),
            d_spacing_angstrom=float(d_spacing[index]),
            relative_intensity=float(intensities[index]),
            family_multiplicity=_family_multiplicity(operators, normals_crystal[index]),
        )
        for index in order
    )
    plane_indices = np.asarray([band.indices for band in bands], dtype=np.int64)
    zone_axes = _zone_axes_for_bands(
        phase=phase,
        bands=bands,
        plane_indices=plane_indices,
        view_matrix=view,
        max_index=zone_axis_max_index,
        min_order=min_zone_axis_order,
        max_polar_angle_deg=float(max_polar_angle_deg),
    )
    return StereographicKikuchiMap(
        phase=phase,
        beam_energy_kev=float(beam_energy_kev),
        wavelength_angstrom=wavelength,
        centre_indices=tuple(int(value) for value in np.asarray(centre_direction).ravel()),
        horizontal_indices=tuple(int(value) for value in np.asarray(horizontal_direction).ravel()),
        view_matrix=view,
        bands=bands,
        zone_axes=zone_axes,
        has_intensity_model=has_intensity_model,
        provenance=provenance,
    )


def _family_multiplicity(operators: np.ndarray, normal: np.ndarray) -> int:
    """Number of distinct orientations of a plane normal under the point group.

    Antipodal normals denote the same plane, so the two senses are merged. This
    is the number of parallel-looking bands the family contributes to the map.
    """

    images = np.einsum("nij,j->ni", operators, normal)
    kept: list[np.ndarray] = []
    for row in images:
        if not any(
            float(np.linalg.norm(row - other)) < 1e-8 or float(np.linalg.norm(row + other)) < 1e-8
            for other in kept
        ):
            kept.append(row)
    return len(kept)


def _zone_axes_for_bands(
    *,
    phase: Phase,
    bands: Sequence[KikuchiMapBand],
    plane_indices: np.ndarray,
    view_matrix: np.ndarray,
    max_index: int,
    min_order: int,
    max_polar_angle_deg: float,
) -> tuple[KikuchiMapZoneAxis, ...]:
    """Zone axes as the crossings of band centre lines, with the bands at each.

    The zone axis common to two planes is the **integer** cross product of their
    Miller indices, because the direct and reciprocal bases are dual: no metric
    tensor appears in ``h u + k v + l w``, so the cross product of two index
    triples is a direction triple. Working in integers rather than in Cartesian
    coordinates makes this exact, makes deduplication a set lookup on tuples
    instead of an O(n^2) sweep over float vectors, and removes the need to search
    a grid for the nearest low-index direction. On a map with a few hundred bands
    -- which a phase carrying no atomic basis produces, since nothing can be
    filtered on intensity -- the difference is minutes against milliseconds.
    """

    if len(bands) < 2:
        return ()
    left, right = np.triu_indices(len(bands), k=1)
    crossings = np.cross(plane_indices[left], plane_indices[right]).astype(np.int64)
    crossings = crossings[np.any(crossings != 0, axis=1)]
    if crossings.shape[0] == 0:
        return ()
    # Lowest terms: a direction is named by its simplest indices, and the multiples
    # denote the same direction.
    divisors = np.gcd.reduce(np.abs(crossings), axis=1)
    crossings //= divisors[:, None]
    # One representative per antipodal pair, canonicalized on the first non-zero
    # component so that [uvw] and its negative -- the same axis, since the beam
    # traverses the crystal either way -- collapse to one entry.
    leading = crossings[np.arange(crossings.shape[0]), np.argmax(crossings != 0, axis=1)]
    crossings = np.where((leading < 0)[:, None], -crossings, crossings)
    crossings = crossings[np.max(np.abs(crossings), axis=1) <= max_index]
    if crossings.shape[0] == 0:
        return ()
    unique_indices = np.unique(crossings, axis=0)

    # Cartesian directions in the map frame, and the zone law for every
    # (axis, band) pair at once.
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=np.float64)
    cartesian = unique_indices.astype(np.float64) @ direct.T
    cartesian = cartesian / np.linalg.norm(cartesian, axis=1)[:, None]
    directions = cartesian @ np.asarray(view_matrix, dtype=np.float64).T
    # A numerically-zero z is snapped to exactly +0.0. An axis on the equator has
    # z of order 1e-17 of either sign, and a one-hemisphere projection folds on the
    # sign of z, so without this it lands on whichever side of the disc the
    # round-off chose -- and the two senses of one axis can land on opposite sides.
    directions[np.abs(directions[:, 2]) < 1e-12, 2] = 0.0
    # Draw each axis in the upper hemisphere. On the equator, where both senses have
    # z = 0, break the tie lexicographically on x then y.
    flip = (directions[:, 2] < 0.0) | (
        (directions[:, 2] == 0.0)
        & ((directions[:, 0] < 0.0) | ((directions[:, 0] == 0.0) & (directions[:, 1] < 0.0)))
    )
    drawn = np.where(flip[:, None], -directions, directions)
    drawn = drawn / np.linalg.norm(drawn, axis=1)[:, None]
    polar_deg = np.degrees(np.arccos(np.clip(drawn[:, 2], -1.0, 1.0)))
    membership = (plane_indices @ unique_indices.T) == 0

    axes: list[KikuchiMapZoneAxis] = []
    for position in range(unique_indices.shape[0]):
        if polar_deg[position] > max_polar_angle_deg + 1e-9:
            continue
        members = tuple(int(index) for index in np.nonzero(membership[:, position])[0])
        if len(members) < min_order:
            continue
        axes.append(
            KikuchiMapZoneAxis(
                phase=phase,
                indices=tuple(int(value) for value in unique_indices[position]),
                direction_map=drawn[position],
                band_indices=members,
                polar_angle_deg=float(polar_deg[position]),
            )
        )
    axes.sort(key=lambda axis: (-axis.order, axis.polar_angle_deg))
    return tuple(axes)


def _not_a_zone_axis_message(phase: Phase, direction: ArrayLike) -> str:
    """Explain that a direction is absent from the map, and what to do about it."""

    indices = [int(value) for value in np.asarray(direction).ravel()]
    notation = _format_direction(phase, indices)
    return (
        f"{notation} is not a zone axis of this map. Raise zone_axis_max_index to admit "
        "higher-index axes, or lower min_zone_axis_order to admit axes crossed by fewer bands."
    )


def plan_kikuchi_route(
    kikuchi_map: StereographicKikuchiMap,
    start: ArrayLike,
    target: ArrayLike,
    *,
    max_leg_deg: float = DEFAULT_ROUTE_MAX_LEG_DEG,
) -> KikuchiRoute:
    r"""Plan a band-following tilt route between two zone axes of a map.

    Purpose
    -------
    Turns the map into an instruction. Two zone axes joined by a band can be
    reached in one hop with a line to follow; otherwise the zone-axis network is
    searched for the cheapest sequence of such hops.

    Method
    ------
    A shortest-path search over the map's zone axes. Nodes are zone axes; an edge
    joins two axes when they share at least one band — equivalently, when one
    reflection is perpendicular to both — and their separation does not exceed
    ``max_leg_deg``. Edge weight is the angular separation, so the result
    minimizes total stage travel among band-followable routes.

    > The reason to prefer several short hops over one long one is not tidiness.
    > An error ``dphi`` in the assumed orientation produces a residual of about
    > ``dphi sin(theta)`` after a hop of length ``theta``, so a long excursion
    > amplifies whatever the starting orientation got wrong, and the operator
    > loses the band. Re-indexing at each intermediate axis turns an open-loop
    > calculation into a closed loop.

    When no route exists, the search falls back to reporting a single bandless
    geodesic leg with ``reachable`` set to ``False``, so the caller always gets a
    describable answer rather than an exception.

    Parameters
    ----------
    kikuchi_map : StereographicKikuchiMap
    start, target : ArrayLike
        Zone-axis direction indices. Both must be axes of the map; use
        :meth:`StereographicKikuchiMap.zone_axis_for_direction` to check, or raise
        the map's ``zone_axis_max_index`` to admit them.
    max_leg_deg : float
        Largest single hop accepted, in degrees.

    Returns
    -------
    KikuchiRoute
        With ``reachable`` False when no band-followable route exists.

    Raises
    ------
    ValueError
        When either endpoint is not a zone axis of the map, or ``max_leg_deg`` is
        not positive.
    """

    if float(max_leg_deg) <= 0.0:
        raise ValueError("max_leg_deg must be strictly positive.")
    origin = kikuchi_map.zone_axis_for_direction(start)
    destination = kikuchi_map.zone_axis_for_direction(target)
    if origin is None:
        raise ValueError(_not_a_zone_axis_message(kikuchi_map.phase, start))
    if destination is None:
        raise ValueError(_not_a_zone_axis_message(kikuchi_map.phase, target))
    start_direction = np.asarray(origin.direction_map, dtype=np.float64)
    target_direction = np.asarray(destination.direction_map, dtype=np.float64)
    if origin.indices == destination.indices:
        return KikuchiRoute(
            phase=kikuchi_map.phase,
            legs=(),
            reachable=True,
            start_indices=origin.indices,
            target_indices=destination.indices,
            start_direction=freeze_array(start_direction),
            target_direction=freeze_array(target_direction),
        )

    axes = kikuchi_map.zone_axes
    position_of = {axis.indices: index for index, axis in enumerate(axes)}
    directions = np.asarray([axis.direction_map for axis in axes], dtype=np.float64)
    separations = np.degrees(
        np.arccos(np.clip(np.abs(directions @ directions.T), -1.0, 1.0))
    )
    band_sets = [set(axis.band_indices) for axis in axes]

    neighbours: list[list[tuple[int, float, int]]] = [[] for _ in axes]
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            common = band_sets[i] & band_sets[j]
            if not common:
                continue
            separation = float(separations[i, j])
            if separation <= 0.0 or separation > float(max_leg_deg):
                continue
            # Follow the strongest shared band: bands are ordered by intensity, so
            # the smallest position is the most visible one.
            best = min(common)
            neighbours[i].append((j, separation, best))
            neighbours[j].append((i, separation, best))

    source = position_of[origin.indices]
    sink = position_of[destination.indices]
    distances = np.full(len(axes), np.inf)
    distances[source] = 0.0
    previous: dict[int, tuple[int, int]] = {}
    queue: list[tuple[float, int]] = [(0.0, source)]
    visited: set[int] = set()
    while queue:
        cost, node = heapq.heappop(queue)
        if node in visited:
            continue
        visited.add(node)
        if node == sink:
            break
        for neighbour, weight, edge in neighbours[node]:
            candidate = cost + weight
            if candidate < distances[neighbour] - 1e-12:
                distances[neighbour] = candidate
                previous[neighbour] = (node, edge)
                heapq.heappush(queue, (candidate, neighbour))

    if not np.isfinite(distances[sink]):
        direct = float(
            np.degrees(
                np.arccos(np.clip(abs(float(np.dot(start_direction, target_direction))), -1.0, 1.0))
            )
        )
        return KikuchiRoute(
            phase=kikuchi_map.phase,
            legs=(
                KikuchiRouteLeg(
                    phase=kikuchi_map.phase,
                    start_indices=origin.indices,
                    end_indices=destination.indices,
                    band_indices=None,
                    tilt_deg=direct,
                    waypoint_indices=(),
                    start_direction=freeze_array(start_direction),
                    end_direction=freeze_array(
                        target_direction
                        if float(np.dot(start_direction, target_direction)) >= 0.0
                        else -target_direction
                    ),
                ),
            ),
            reachable=False,
            start_indices=origin.indices,
            target_indices=destination.indices,
            start_direction=freeze_array(start_direction),
            target_direction=freeze_array(target_direction),
        )

    chain: list[tuple[int, int, int]] = []
    node = sink
    while node != source:
        parent, edge = previous[node]
        chain.append((parent, node, edge))
        node = parent
    chain.reverse()

    # Orient the chain. Each node is a line, so its stored direction has an
    # arbitrary sense; the route must take the sense that continues from the
    # previous node, and the reported indices must carry that same sense or the
    # printed route and the drawn one disagree about which side of the projection
    # the path runs down.
    legs: list[KikuchiRouteLeg] = []
    current = np.asarray(start_direction, dtype=np.float64)
    current_indices = origin.indices
    for _parent, child, band_position in chain:
        stored = np.asarray(axes[child].direction_map, dtype=np.float64)
        sign = 1.0 if float(np.dot(current, stored)) >= 0.0 else -1.0
        # Adding zero is not a no-op here. Negating a component that is exactly
        # zero produces -0.0, and a one-hemisphere projection decides whether to
        # fold a direction onto its antipode from the *sign bit* of z -- which is
        # set for -0.0. An axis on the equator, where z is exactly zero, would
        # then be drawn on the opposite side of the disc from the arc that reaches
        # it. Adding zero collapses -0.0 to +0.0 and leaves everything else alone.
        oriented: np.ndarray = np.asarray(sign * stored, dtype=np.float64) + 0.0
        oriented_indices = tuple(int(sign * value) for value in axes[child].indices)
        legs.append(
            KikuchiRouteLeg(
                phase=kikuchi_map.phase,
                start_indices=current_indices,
                end_indices=oriented_indices,
                band_indices=kikuchi_map.bands[band_position].indices,
                tilt_deg=float(
                    np.degrees(
                        np.arccos(np.clip(float(np.dot(current, oriented)), -1.0, 1.0))
                    )
                ),
                waypoint_indices=_waypoints_on_band(
                    kikuchi_map, band_position, current, oriented
                ),
                start_direction=freeze_array(current),
                end_direction=freeze_array(oriented),
            )
        )
        current, current_indices = oriented, oriented_indices
    return KikuchiRoute(
        phase=kikuchi_map.phase,
        legs=tuple(legs),
        reachable=True,
        start_indices=origin.indices,
        target_indices=current_indices,
        start_direction=freeze_array(start_direction),
        target_direction=freeze_array(current),
    )


def _waypoints_on_band(
    kikuchi_map: StereographicKikuchiMap,
    band_position: int,
    start: np.ndarray,
    end: np.ndarray,
) -> tuple[tuple[int, ...], ...]:
    """Zone axes on one band, strictly between two oriented directions, in order.

    These are the landmarks that tell an operator the tilt is tracking the right
    band. A candidate lies on the arc exactly when the two angles it subtends at
    the ends sum to the whole span, which is the spherical statement of "between".
    Ordering is by the angle turned from the start, and the reported indices carry
    the sense in which the axis is actually passed.
    """

    position = int(band_position)
    first = np.asarray(start, dtype=np.float64)
    second = np.asarray(end, dtype=np.float64)
    span = float(np.degrees(np.arccos(np.clip(float(np.dot(first, second)), -1.0, 1.0))))
    if span <= 0.0:
        return ()
    passed: list[tuple[float, tuple[int, ...]]] = []
    for axis in kikuchi_map.zone_axes:
        if position not in axis.band_indices:
            continue
        candidate = np.asarray(axis.direction_map, dtype=np.float64)
        for sign in (1.0, -1.0):
            sense = sign * candidate
            from_start = float(
                np.degrees(np.arccos(np.clip(float(np.dot(first, sense)), -1.0, 1.0)))
            )
            to_end = float(
                np.degrees(np.arccos(np.clip(float(np.dot(sense, second)), -1.0, 1.0)))
            )
            if abs(from_start + to_end - span) < 1e-6 and from_start > 1e-6 and to_end > 1e-6:
                passed.append((from_start, tuple(int(sign * value) for value in axis.indices)))
                break
    passed.sort()
    return tuple(indices for _, indices in passed)
