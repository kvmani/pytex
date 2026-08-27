"""Synthetic selected-area diffraction patterns, on a detector raster.

Why this module exists
----------------------
:func:`pytex.diffraction.saed.generate_saed_pattern` answers the crystallographic
question — which reflections lie on a zone, where they land in detector
millimetres, and how strong they are. It stops there, and that is one step short
of the thing a microscopist works with. What arrives from a camera is a raster of
pixels: an image with a beam stop near the middle, spots of finite size, a
strongest spot that saturates, weak reflections that fade into the background,
and a field of view that simply ends.

This module supplies that last step. It projects a simulated zone-axis pattern
onto a stated detector raster and returns the spots in **pixel coordinates**,
with the apparent radius and the relative brightness a renderer needs. The result
is a pattern that can be clicked on and indexed by exactly the same path a real
plate takes — pick the beam, pick the spots, hand the coordinates to
:func:`pytex.diffraction.solving.solve_saed_pattern` — with the difference that
here the answer is known by construction, so the workflow can be taught, checked,
and regression-tested without a micrograph.

Conventions
-----------
The raster is the recorded image: column index increases along the detector's
``+X`` axis and row index increases along its ``+Y`` axis, which is how a camera
reads out and how every image viewer displays the result. **No handedness flip is
applied**, so a picked ``(column, row)`` pair is the detector ``(X, Y)`` pair
divided by the pixel pitch and offset to the beam centre, and a pattern
synthesized here is interpreted by the solver under exactly the convention it was
built with. See :doc:`/theory/reciprocal_space_and_kinematic_spots` for the
reciprocal-space geometry and
:doc:`/theory/saed_ratio_angle_indexing` for the indexing that consumes it.

Limits
------
Kinematic, inherited from :func:`generate_saed_pattern`: relative intensities are
indicative rather than quantitative and no dynamical redistribution occurs.
Double diffraction is off by default and added on request through
``include_double_diffraction``; the reflections it puts on the plate are marked,
and their brightness is an observability estimate rather than an intensity.
Apparent spot radius is a *display* model, not a
measurement of the disc size an aperture would produce. What is exact is the
geometry — every spot sits where the lattice, the zone axis, and the camera
constant put it — and geometry is what indexing uses.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.provenance import ProvenanceRecord

__all__ = [
    "DetectorRaster",
    "SyntheticSAEDImage",
    "SyntheticSpot",
    "synthesize_saed_image",
]


#: Relative intensity below which a reflection is treated as invisible.
#:
#: A kinematic zone contains reflections spanning many orders of magnitude, and a
#: plate records the ones above its noise floor. One part in five hundred of the
#: strongest spot is a deliberately generous floor: it keeps the weak-but-real
#: reflections a microscopist does use for indexing while dropping the ones no
#: exposure shows.
DEFAULT_INTENSITY_FLOOR = 0.002


@dataclass(frozen=True, slots=True)
class DetectorRaster:
    """The pixel grid a synthetic pattern is projected onto.

    Purpose
    -------
    Fixes the three numbers that turn detector millimetres into something
    clickable: how many pixels there are, how large each one is, and where the
    transmitted beam sits. The beam is deliberately allowed to be off-centre,
    because on a real plate it usually is, and an indexing workflow that only
    works when the beam is at the exact middle of the frame has not been tested.

    Parameters
    ----------
    width_px, height_px : int
        Raster size. Must be positive.
    pixel_size_mm : float
        Detector pixel pitch. The same quantity the solver's ``pixel_size_mm``
        calibration field takes, and it must agree with it.
    centre_px : tuple of float, optional
        Where the transmitted beam lands, in ``(column, row)``. Defaults to the
        geometric centre of the raster.

    Examples
    --------
    >>> raster = DetectorRaster(width_px=512, height_px=512, pixel_size_mm=0.05)
    >>> raster.centre
    (256.0, 256.0)
    >>> round(raster.half_extent_mm(), 4)
    12.8
    """

    width_px: int = 1024
    height_px: int = 1024
    pixel_size_mm: float = 0.05
    centre_px: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        if self.width_px <= 0 or self.height_px <= 0:
            raise ValueError("DetectorRaster width_px and height_px must be strictly positive.")
        if not math.isfinite(self.pixel_size_mm) or self.pixel_size_mm <= 0.0:
            raise ValueError("DetectorRaster.pixel_size_mm must be finite and strictly positive.")
        if self.centre_px is not None:
            centre = tuple(float(value) for value in self.centre_px)
            if len(centre) != 2 or not all(math.isfinite(value) for value in centre):
                raise ValueError("DetectorRaster.centre_px must be two finite numbers.")
            object.__setattr__(self, "centre_px", centre)

    @property
    def centre(self) -> tuple[float, float]:
        """The transmitted-beam position in pixels, defaulted to the raster centre."""

        if self.centre_px is not None:
            return self.centre_px
        return (self.width_px / 2.0, self.height_px / 2.0)

    def half_extent_mm(self) -> float:
        """Largest distance from the beam to a frame edge, in millimetres.

        This is the radial cut-off a recorded pattern has: a reflection farther
        out than this is off the plate, and simulating it would put a spot where
        no user can click.
        """

        column, row = self.centre
        reach_px = max(
            math.hypot(column, row),
            math.hypot(self.width_px - column, row),
            math.hypot(column, self.height_px - row),
            math.hypot(self.width_px - column, self.height_px - row),
        )
        return float(reach_px * self.pixel_size_mm)


@dataclass(frozen=True, slots=True)
class SyntheticSpot:
    """One reflection of a synthetic pattern, placed on the detector raster.

    Attributes
    ----------
    miller_indices : np.ndarray
        The reflection ``(hkl)``, in the three-index basis of the phase.
    position_px : np.ndarray
        ``(column, row)`` on the raster.
    g_inv_angstrom : float
        Length of the reciprocal-lattice vector; the measurable the radius
        encodes.
    d_spacing_angstrom : float
        ``1 / g_inv_angstrom``.
    relative_intensity : float
        Kinematic intensity, normalized so the strongest spot in the pattern is
        ``1.0`` — unless :attr:`is_double_diffraction` is set, in which case it
        is an indicative observability estimate and not a kinematic intensity.
    apparent_radius_px : float
        Display radius. A presentation quantity — see the module limits.
    label : str
        The rendered ``(hkl)`` label.
    is_double_diffraction : bool
        Whether this reflection is kinematically forbidden and present only
        because double diffraction was requested.
    double_diffraction_origin : str
        ``"(g1) + (g2)"``, the two-step path that puts a marked reflection on
        the plate, in the phase's three-index basis. Empty for an ordinary
        kinematic reflection. A caller that renders indices in another
        convention — four-index Miller-Bravais for a hexagonal phase, say —
        should format :attr:`double_diffraction_parents` itself rather than
        reuse this string, so one row does not mix two notations.
    double_diffraction_parents : np.ndarray, optional
        ``(2, 3)`` integer pair summing to :attr:`miller_indices`, present
        exactly for a marked reflection.
    """

    miller_indices: np.ndarray
    position_px: np.ndarray
    g_inv_angstrom: float
    d_spacing_angstrom: float
    relative_intensity: float
    apparent_radius_px: float
    label: str
    is_double_diffraction: bool = False
    double_diffraction_origin: str = ""
    double_diffraction_parents: np.ndarray | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        object.__setattr__(self, "position_px", as_float_array(self.position_px, shape=(2,)))
        if not 0.0 <= self.relative_intensity <= 1.0:
            raise ValueError("SyntheticSpot.relative_intensity must lie in [0, 1].")
        if not math.isfinite(self.g_inv_angstrom) or self.g_inv_angstrom <= 0.0:
            raise ValueError("SyntheticSpot.g_inv_angstrom must be finite and strictly positive.")
        if not math.isfinite(self.apparent_radius_px) or self.apparent_radius_px <= 0.0:
            raise ValueError("SyntheticSpot.apparent_radius_px must be finite and positive.")
        if bool(self.double_diffraction_origin) and not self.is_double_diffraction:
            raise ValueError(
                "SyntheticSpot.double_diffraction_origin is set on a spot that is not marked "
                "as double diffraction."
            )
        if self.double_diffraction_parents is not None:
            parents = as_int_array(self.double_diffraction_parents, shape=(2, 3))
            if not np.array_equal(parents.sum(axis=0), self.miller_indices):
                raise ValueError(
                    "SyntheticSpot.double_diffraction_parents must sum to the reflection they "
                    "produce."
                )
            object.__setattr__(self, "double_diffraction_parents", parents)


@dataclass(frozen=True, slots=True)
class SyntheticSAEDImage:
    """A zone-axis pattern rendered onto a detector raster, with its own answer.

    Purpose
    -------
    The transport object between the simulation and anything that draws or
    indexes the pattern. It carries both what a user can see — pixel positions,
    brightnesses, radii — and what a user is trying to work out — the phase, the
    zone axis, and the index of every spot. Keeping the answer attached is what
    makes the object usable as a teaching exercise and as a round-trip test:
    index the coordinates, then compare with the construction.

    Attributes
    ----------
    phase : Phase
    zone_axis : ZoneAxis
        The beam direction the pattern was built from — the answer to "which
        axis is this?".
    raster : DetectorRaster
    camera_constant_mm_angstrom : float
        The ``L·λ`` the radii were scaled by; the calibration a user must supply
        to index the pattern, and the one number that, set wrongly, indexes a
        pattern to the wrong phase.
    in_plane_rotation_deg : float
        How far the pattern was rolled about the beam before projection. Zero
        would align the pattern with the detector axes, which a real plate never
        does.
    spots : tuple of SyntheticSpot
        Ordered by decreasing intensity, then increasing radius.
    provenance : ProvenanceRecord, optional
    """

    phase: Phase
    zone_axis: ZoneAxis
    raster: DetectorRaster
    camera_constant_mm_angstrom: float
    in_plane_rotation_deg: float
    spots: tuple[SyntheticSpot, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "spots", tuple(self.spots))
        if self.camera_constant_mm_angstrom <= 0.0:
            raise ValueError(
                "SyntheticSAEDImage.camera_constant_mm_angstrom must be strictly positive."
            )

    @property
    def centre_px(self) -> tuple[float, float]:
        """Where the transmitted beam is, in pixels."""

        return self.raster.centre

    def brightest(self, count: int) -> tuple[SyntheticSpot, ...]:
        """The ``count`` strongest spots, in the pattern's own order.

        The spots a microscopist picks first, and therefore the ones an
        auto-picking convenience should offer.
        """

        if count < 0:
            raise ValueError("count must be non-negative.")
        return self.spots[:count]

    def independent_seed_spots(self, count: int) -> tuple[SyntheticSpot, ...]:
        """Strong spots chosen so the first two are not collinear with the beam.

        Purpose
        -------
        Indexing seeds on two non-collinear reflections. The strongest spots of a
        pattern very often include a Friedel pair — ``g`` and ``-g`` — which is
        collinear through the beam and cannot seed anything. Picking the top
        ``n`` by brightness alone therefore produces a set that looks generous
        and indexes nothing. This walks the same brightness order but skips a
        spot whose direction from the beam duplicates one already taken, so the
        selection spans the zone.

        Parameters
        ----------
        count : int
            How many spots to return. Fewer are returned if the pattern has
            fewer independent directions.
        """

        if count < 0:
            raise ValueError("count must be non-negative.")
        centre = np.asarray(self.centre_px, dtype=float)
        chosen: list[SyntheticSpot] = []
        directions: list[np.ndarray] = []
        for spot in self.spots:
            if len(chosen) >= count:
                break
            offset = np.asarray(spot.position_px, dtype=float) - centre
            norm = float(np.linalg.norm(offset))
            if norm < 1e-9:
                continue
            unit = offset / norm
            # `abs` folds g and -g together: a Friedel pair carries one direction
            # between them, not two.
            if any(abs(float(np.dot(unit, other))) > 0.999 for other in directions):
                continue
            directions.append(unit)
            chosen.append(spot)
        return tuple(chosen)

    def crystal_to_pattern(self) -> np.ndarray:
        """The rotation taking crystal Cartesian vectors into the pattern frame.

        Purpose
        -------
        State the orientation this pattern *was built from*, in the same form
        the indexer reports and the Kikuchi overlay consumes: a rotation matrix
        whose third row is the beam direction, so that ``R @ g`` has the
        detector coordinates of the reflection ``g`` in its first two components
        and its excitation error in the third.

        When and where to use it
        ------------------------
        Whenever something other than the spots has to be drawn on the same
        pattern — Kikuchi bands, a calculated overlay, a stereogram — and must
        agree with it. Without this the orientation would have to be recovered
        by indexing a pattern whose answer is already known, and any convention
        mismatch would show up as a silently rotated overlay.

        Returns
        -------
        numpy.ndarray
            A ``(3, 3)`` proper rotation. Multiply a reciprocal-lattice vector
            in crystal Cartesian coordinates by it to obtain
            ``(x, y, excitation)`` in the pattern frame, in the same units;
            multiply the first two components by
            ``camera_constant_mm_angstrom / pixel_size_mm`` and add
            :attr:`centre_px` to land on the drawn spot.

        Notes
        -----
        The construction mirrors :func:`synthesize_saed_image` exactly: the
        deterministic zone basis of :func:`~pytex.diffraction.kinematic.zone_basis_from_axis`
        for the beam direction, then the roll about the beam that the projection
        applied to the detector coordinates.

        Examples
        --------
        >>> from pytex.app.phases import builtin_phase
        >>> from pytex.core.lattice import ZoneAxis
        >>> phase = builtin_phase("ni_fcc").to_phase()
        >>> axis = ZoneAxis(indices=(0, 0, 1), phase=phase)
        >>> image = synthesize_saed_image(phase, axis, in_plane_rotation_deg=17.0)
        >>> matrix = image.crystal_to_pattern()
        >>> bool(abs(float(np.linalg.det(matrix)) - 1.0) < 1e-12)
        True

        The third row is the beam: the zone axis maps onto pattern ``+z``.

        >>> beam = matrix @ axis.unit_vector
        >>> bool(abs(beam[2] - 1.0) < 1e-12)
        True
        """

        from pytex.diffraction.kinematic import zone_basis_from_axis

        basis = np.asarray(zone_basis_from_axis(self.zone_axis.unit_vector), dtype=float)
        angle = math.radians(float(self.in_plane_rotation_deg))
        cosine, sine = math.cos(angle), math.sin(angle)
        roll = np.array([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]], dtype=float)
        return np.ascontiguousarray(roll @ basis.T)

    def describe(self) -> str:
        """A convention-explicit prose account of what this pattern is.

        Written for a reader who has the picture in front of them and needs to
        know what was assumed to make it: the phase, the axis, the calibration,
        and the honest statement that the intensities are kinematic.
        """

        zone = format_direction_indices(
            tuple(int(value) for value in self.zone_axis.indices), style="plain"
        )
        if not self.spots:
            return (
                f"A simulated {self.phase.name} pattern down {zone} at a camera constant of "
                f"{self.camera_constant_mm_angstrom:g} mm·Å contains no reflection inside the "
                "recorded field of view. Either the camera constant is too large for this "
                "detector or the intensity floor has removed every spot."
            )
        strongest = self.spots[0]
        widest = max(self.spots, key=lambda spot: spot.g_inv_angstrom)
        return (
            f"Simulated selected-area pattern of {self.phase.name} with the beam along {zone}, "
            f"projected onto a {self.raster.width_px}x{self.raster.height_px} raster of "
            f"{self.raster.pixel_size_mm:g} mm pixels at a camera constant of "
            f"{self.camera_constant_mm_angstrom:g} mm·Å, rolled "
            f"{self.in_plane_rotation_deg:g}° about the beam. "
            f"{len(self.spots)} reflections fall inside the field of view, the strongest being "
            f"{strongest.label} at d = {strongest.d_spacing_angstrom:.4f} Å and the farthest "
            f"{widest.label} at d = {widest.d_spacing_angstrom:.4f} Å. "
            "Spot positions are exact for this lattice and zone axis; relative intensities are "
            "kinematic and therefore indicative. " + self._double_diffraction_sentence()
        )

    def _double_diffraction_sentence(self) -> str:
        """The honest one-line statement of what double diffraction contributed.

        Two patterns of the same phase and axis differ visibly depending on
        whether the option was on, so the description has to say which one this
        is rather than leaving the reader to infer it from the spot count.
        """

        marked = [spot for spot in self.spots if spot.is_double_diffraction]
        if not marked:
            return (
                "Double diffraction is not included, so a forbidden reflection that a real "
                "plate shows through multiple scattering is absent here."
            )
        labels = ", ".join(spot.label for spot in marked[:3])
        more = "" if len(marked) <= 3 else f", and {len(marked) - 3} more"
        return (
            f"Double diffraction is included: {len(marked)} kinematically forbidden "
            f"reflections ({labels}{more}) are drawn because two excited reflections sum to "
            "them. Their brightness is an observability estimate scaled by a coupling "
            "constant, not a kinematic intensity."
        )

    def to_json(self) -> dict[str, Any]:
        """The pattern as JSON-ready data: geometry, brightness, and the answer."""

        return {
            "phase_name": self.phase.name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "camera_constant_mm_angstrom": float(self.camera_constant_mm_angstrom),
            "in_plane_rotation_deg": float(self.in_plane_rotation_deg),
            "width_px": int(self.raster.width_px),
            "height_px": int(self.raster.height_px),
            "pixel_size_mm": float(self.raster.pixel_size_mm),
            "centre_px": [float(value) for value in self.centre_px],
            "crystal_to_pattern": [float(value) for value in self.crystal_to_pattern().reshape(-1)],
            "spots": [
                {
                    "hkl": [int(value) for value in spot.miller_indices],
                    "label": spot.label,
                    "x": float(spot.position_px[0]),
                    "y": float(spot.position_px[1]),
                    "g_inv_angstrom": float(spot.g_inv_angstrom),
                    "d_angstrom": float(spot.d_spacing_angstrom),
                    "intensity": float(spot.relative_intensity),
                    "radius_px": float(spot.apparent_radius_px),
                    "double_diffraction": bool(spot.is_double_diffraction),
                    "double_diffraction_origin": spot.double_diffraction_origin,
                    "double_diffraction_parents": (
                        None
                        if spot.double_diffraction_parents is None
                        else [
                            [int(value) for value in parent]
                            for parent in spot.double_diffraction_parents
                        ]
                    ),
                }
                for spot in self.spots
            ],
            "describe": self.describe(),
        }


def _apparent_radius_px(relative_intensity: float, *, base_radius_px: float, bloom: float) -> float:
    """Display radius of a spot of the given relative intensity.

    A recorded spot's size is set by the illumination and the aperture, not by
    how strong the reflection is; what makes a bright spot *look* larger on a
    plate is that its tails clear the display threshold. That is a
    photographic effect, so it is modelled photographically: a common base
    radius with a bounded brightness-driven bloom, on the fourth root so the
    weakest visible reflection is small but still clickable.
    """

    return float(base_radius_px * (1.0 + bloom * float(relative_intensity) ** 0.25))


def synthesize_saed_image(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    camera_constant_mm_angstrom: float = 180.0,
    raster: DetectorRaster | None = None,
    in_plane_rotation_deg: float = 0.0,
    max_index: int = 6,
    intensity_floor: float = DEFAULT_INTENSITY_FLOOR,
    include_double_diffraction: bool = False,
    double_diffraction_coupling: float = 0.05,
    max_spots: int = 400,
    base_radius_px: float = 5.0,
    bloom: float = 1.4,
    position_jitter_px: float = 0.0,
    rng_seed: int = 0,
    provenance: ProvenanceRecord | None = None,
) -> SyntheticSAEDImage:
    """Project a simulated zone-axis pattern onto a detector raster.

    Purpose
    -------
    Produce a diffraction pattern that can be clicked on, indexed, and checked —
    a practice plate whose answer is known by construction. Use it to try the
    indexing workflow without a micrograph, to teach it, and to regression-test
    the whole picking-to-indexing path end to end.

    When and where to use it
    ------------------------
    Use it wherever a *pattern* is needed rather than a *spot list*: the
    workbench's TEM gallery, a tutorial, or a test that must exercise
    calibration and picking as well as indexing. When only the crystallography is
    wanted, call :func:`pytex.diffraction.saed.generate_saed_pattern` directly;
    this function is that one plus a detector.

    Parameters
    ----------
    phase : Phase
        The material. Its unit cell drives the structure factors, so a phase
        without a basis gives every allowed reflection equal weight.
    zone_axis : ZoneAxis
        Beam direction, in the phase's three-index basis. Must belong to
        ``phase``.
    camera_constant_mm_angstrom : float
        ``L·λ``. Radial position on the detector is this divided by ``d``.
    raster : DetectorRaster, optional
        The pixel grid. Defaults to a 1024x1024 raster of 0.05 mm pixels with
        the beam at its centre.
    in_plane_rotation_deg : float
        Roll about the beam applied before projection, positive from the
        detector ``+X`` axis towards ``+Y``. A real plate is never aligned with
        the detector axes, and a workflow tested only at zero has not been
        tested.
    max_index : int
        Largest ``|h|``, ``|k|``, ``|l|`` enumerated.
    intensity_floor : float
        Relative-intensity threshold below which a reflection is dropped as
        invisible. See :data:`DEFAULT_INTENSITY_FLOOR`.
    include_double_diffraction : bool
        Whether to show the kinematically forbidden reflections that a real
        plate displays because a diffracted beam re-diffracts. Default
        ``False``, which leaves the pattern purely kinematic. The reflections it
        adds are marked by ``SyntheticSpot.is_double_diffraction`` and carry the
        two-step path in ``double_diffraction_origin``; treat their brightness
        as an observability estimate, not an intensity. See
        :func:`pytex.diffraction.saed.generate_saed_pattern`.
    double_diffraction_coupling : float
        Scale applied to the two-step weight, in ``(0, 1]``. Ignored when
        ``include_double_diffraction`` is ``False``.
    max_spots : int
        Hard cap on returned spots, applied after sorting by intensity. Guards
        the renderer against a large-cell phase producing thousands of spots.
    base_radius_px, bloom : float
        Display-radius model. See :func:`_apparent_radius_px`.
    position_jitter_px : float
        Standard deviation of an isotropic Gaussian displacement applied to each
        spot, emulating the spot-centroiding error of a real measurement. Zero
        gives exact positions. Any non-zero value makes indexing residuals
        realistic rather than machine-epsilon.
    rng_seed : int
        Seed for the jitter, so a gallery entry looks the same every time it is
        opened.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    SyntheticSAEDImage
        Spots in pixel coordinates, ordered strongest first, with the
        construction's phase and zone axis attached.

    Raises
    ------
    ValueError
        If the zone axis belongs to another phase, or any numeric argument is
        outside its admissible range.

    Examples
    --------
    Every spot of an fcc [001] pattern lies in the ``l = 0`` layer, and each one
    sits at exactly ``(camera constant) / d`` from the beam — which is the
    identity the whole indexing chain rests on::

        >>> from pytex.core.fixtures import get_phase_fixture
        >>> phase = get_phase_fixture("ni_fcc").to_phase()
        >>> axis = ZoneAxis(indices=(0, 0, 1), phase=phase)
        >>> image = synthesize_saed_image(phase, axis, camera_constant_mm_angstrom=180.0)
        >>> all(int(spot.miller_indices[2]) == 0 for spot in image.spots)
        True
        >>> spot = image.spots[0]
        >>> centre = image.centre_px
        >>> radius_px = float(np.hypot(*(spot.position_px - np.asarray(centre))))
        >>> radius_mm = radius_px * image.raster.pixel_size_mm
        >>> bool(abs(radius_mm - 180.0 / spot.d_spacing_angstrom) < 1e-9)
        True

    See Also
    --------
    pytex.diffraction.saed.generate_saed_pattern : the crystallography this wraps.
    pytex.diffraction.solving.solve_saed_pattern : the indexer that consumes the picks.
    """

    from pytex.diffraction.saed import generate_saed_pattern

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if not 0.0 <= intensity_floor < 1.0:
        raise ValueError("intensity_floor must lie in [0, 1).")
    if max_spots <= 0:
        raise ValueError("max_spots must be strictly positive.")
    if base_radius_px <= 0.0:
        raise ValueError("base_radius_px must be strictly positive.")
    if bloom < 0.0:
        raise ValueError("bloom must be non-negative.")
    if position_jitter_px < 0.0:
        raise ValueError("position_jitter_px must be non-negative.")

    grid = raster if raster is not None else DetectorRaster()
    half_extent_mm = grid.half_extent_mm()
    pattern = generate_saed_pattern(
        phase,
        zone_axis,
        camera_constant_mm_angstrom=camera_constant_mm_angstrom,
        max_index=max_index,
        # A reflection beyond the plate edge cannot be clicked, so it is not
        # simulated. Converting the frame's reach into the reciprocal cut-off
        # keeps that decision in one place.
        max_g_inv_angstrom=half_extent_mm / camera_constant_mm_angstrom,
        label_limit=max_spots,
        include_double_diffraction=bool(include_double_diffraction),
        double_diffraction_coupling=float(double_diffraction_coupling),
    )
    if not pattern.spots:
        return SyntheticSAEDImage(
            phase=phase,
            zone_axis=zone_axis,
            raster=grid,
            camera_constant_mm_angstrom=float(camera_constant_mm_angstrom),
            in_plane_rotation_deg=float(in_plane_rotation_deg),
            spots=(),
            provenance=provenance,
        )

    coordinates_mm = np.asarray(
        [spot.detector_coordinates for spot in pattern.spots], dtype=float
    ).reshape(-1, 2)
    intensities = np.asarray([spot.intensity for spot in pattern.spots], dtype=float)
    peak = float(intensities.max())
    relative = intensities / peak if peak > 0.0 else np.ones_like(intensities)

    angle = math.radians(float(in_plane_rotation_deg))
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], dtype=float
    )
    rotated_mm = coordinates_mm @ rotation.T
    centre = np.asarray(grid.centre, dtype=float)
    positions = centre + rotated_mm / grid.pixel_size_mm
    if position_jitter_px > 0.0:
        generator = np.random.default_rng(int(rng_seed))
        positions = positions + generator.normal(
            0.0, float(position_jitter_px), size=positions.shape
        )

    inside = (
        (positions[:, 0] >= 0.0)
        & (positions[:, 0] <= float(grid.width_px))
        & (positions[:, 1] >= 0.0)
        & (positions[:, 1] <= float(grid.height_px))
    )
    visible = inside & (relative >= float(intensity_floor))

    order = np.argsort(-relative, kind="stable")
    spots: list[SyntheticSpot] = []
    for index in order:
        if not bool(visible[index]):
            continue
        if len(spots) >= max_spots:
            break
        source = pattern.spots[int(index)]
        reciprocal_vector = np.asarray(source.reciprocal_vector_crystal, dtype=float)
        g_magnitude = float(np.linalg.norm(reciprocal_vector))
        indices = tuple(int(value) for value in np.asarray(source.miller_indices, dtype=int))
        spots.append(
            SyntheticSpot(
                miller_indices=np.asarray(indices, dtype=np.int64),
                position_px=positions[index],
                g_inv_angstrom=g_magnitude,
                d_spacing_angstrom=1.0 / g_magnitude,
                relative_intensity=float(relative[index]),
                apparent_radius_px=_apparent_radius_px(
                    float(relative[index]), base_radius_px=base_radius_px, bloom=bloom
                ),
                label=format_plane_indices(indices, style="plain"),
                is_double_diffraction=source.is_double_diffraction,
                double_diffraction_origin=source.double_diffraction_origin_label(),
                double_diffraction_parents=source.double_diffraction_parents,
            )
        )

    return SyntheticSAEDImage(
        phase=phase,
        zone_axis=zone_axis,
        raster=grid,
        camera_constant_mm_angstrom=float(camera_constant_mm_angstrom),
        in_plane_rotation_deg=float(in_plane_rotation_deg),
        spots=tuple(spots),
        provenance=provenance,
    )
