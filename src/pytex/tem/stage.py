"""Holder tilt kinematics: the map from stage readouts to specimen orientation.

The question this module answers is the mechanical half of TEM tilt navigation:
*given holder angles, where does the specimen actually point?* Everything else in
`pytex.tem` is an inversion of, or a constraint on, the forward map defined here.

The geometry
------------
A double-tilt holder has two axes that are **not** both fixed. The ``alpha`` axis
is the holder rod, fixed in the column. The ``beta`` axis is a cradle carried
*inside* the rod, so it rotates when ``alpha`` does. Composing the beta rotation
about its instantaneous laboratory-frame axis with the alpha rotation gives

    R = [Rx(a) Ry(b) Rx(a)^T] Rx(a) = Rx(a) Ry(b)

so the moving-axis subtlety cancels exactly. That cancellation is a property of
this particular axis pair, not a general licence to ignore axis motion: it does
**not** survive a non-orthogonal or mis-set pair, which is why `GeneralStageAxes`
composes the two rotations explicitly rather than reusing the shortcut.

The beam direction in holder coordinates follows by transposition,

    b_H(a, b) = R^T z = (-cos a sin b, sin a, cos a cos b)

which is spherical coordinates whose **pole is the beta axis**. Three consequences
drive the rest of the package: constant-beta curves are great circles through the
beta pole and constant-alpha curves are small circles about it (so a reachable
region draws exactly, as circles); the Jacobian is ``cos a``, giving the accessible
solid angle in closed form; and ``|d b/d beta| = cos a``, so a degree of beta buys
only ``cos a`` degrees of crystal rotation — the double-tilt holder's gimbal lock.

Frames
------
Per `docs/standards/notation_and_conventions.md` this module introduces no new
frame domain. The **holder frame is the specimen-domain frame** for TEM work, and
the microscope frame is the laboratory-domain frame; `RotationAxis` vectors are
components in whichever of the two the attribute documentation names.

See `docs/architecture/tem_tilt_navigation_foundation.md` sections 3 and 10 for
the derivations and for the reachability geometry these envelopes express.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.provenance import ProvenanceRecord

__all__ = [
    "BEAM_AXIS_LABORATORY",
    "DoubleTiltStage",
    "EllipticalEnvelope",
    "GeneralStageAxes",
    "HolderKind",
    "MaskedEnvelope",
    "PolygonEnvelope",
    "RectangularEnvelope",
    "SingleTiltStage",
    "StageCalibration",
    "StageModel",
    "StagePosition",
    "TiltEnvelope",
    "TiltRotateStage",
    "beam_direction_holder",
    "rotation_x",
    "rotation_y",
    "rotation_z",
]

#: The electron-beam reference axis in laboratory coordinates.
#:
#: Points **up the column, toward the gun**; electrons propagate along its
#: negative. A zone axis is on-axis when its laboratory image is parallel or
#: antiparallel to this vector.
BEAM_AXIS_LABORATORY = np.array([0.0, 0.0, 1.0], dtype=np.float64)
BEAM_AXIS_LABORATORY.setflags(write=False)

#: Below this in-plane magnitude the beta angle of the closed-form solution is
#: indeterminate: the target lies along the beta axis and the solution degenerates
#: to a one-parameter family. See `foundation` section 6.2.
GIMBAL_TOLERANCE = 1e-9


def rotation_x(angle_rad: float) -> np.ndarray:
    """Right-handed rotation about the ``x`` axis, as a 3x3 matrix."""

    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=np.float64)


def rotation_y(angle_rad: float) -> np.ndarray:
    """Right-handed rotation about the ``y`` axis, as a 3x3 matrix."""

    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=np.float64)


def rotation_z(angle_rad: float) -> np.ndarray:
    """Right-handed rotation about the ``z`` axis, as a 3x3 matrix."""

    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


def _rotation_about(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation about an arbitrary unit ``axis``."""

    x, y, z = (float(component) for component in axis)
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]], dtype=np.float64)
    return (
        np.eye(3, dtype=np.float64)
        + math.sin(angle_rad) * cross
        + (1.0 - math.cos(angle_rad)) * (cross @ cross)
    )


def beam_direction_holder(alpha_deg: ArrayLike, beta_deg: ArrayLike) -> np.ndarray:
    """Beam direction in holder coordinates for an ideal double-tilt stage.

    Purpose
    -------
    The closed form of equation (B) of the TEM tilt-navigation foundation,

        b_H = (-cos a sin b, sin a, cos a cos b),

    vectorized over the inputs. This is the *ideal-stage* expression; a
    calibrated stage should be asked through :meth:`StageModel.beam_direction`,
    which routes through the actual axis model.

    When to use
    -----------
    Reachability sweeps, stereographic trajectory sampling, and any place a large
    grid of tilt positions must be turned into directions at once. Prefer it over
    building rotation matrices in a loop: the map is analytic and this evaluates
    it in one vectorized pass.

    Parameters
    ----------
    alpha_deg, beta_deg : array_like
        Tilt angles in degrees, broadcast against one another.

    Returns
    -------
    np.ndarray
        Unit vectors with shape ``broadcast_shape + (3,)``, in holder
        coordinates.

    Examples
    --------
    At zero tilt the beam lies along the holder ``z`` axis::

        >>> import numpy as np
        >>> np.round(beam_direction_holder(0.0, 0.0), 12)
        array([0., 0., 1.])
    """

    alpha = np.deg2rad(np.asarray(alpha_deg, dtype=np.float64))
    beta = np.deg2rad(np.asarray(beta_deg, dtype=np.float64))
    alpha, beta = np.broadcast_arrays(alpha, beta)
    cos_alpha = np.cos(alpha)
    return np.stack(
        [-cos_alpha * np.sin(beta), np.sin(alpha), cos_alpha * np.cos(beta)],
        axis=-1,
    )


class HolderKind(StrEnum):
    """The tilt geometry a holder provides.

    The kind determines what the second angle *means* and therefore what the
    reachable set looks like: a surface for ``DOUBLE_TILT`` and
    ``TILT_ROTATE``, but only a curve — a set of measure zero — for
    ``SINGLE_TILT``, where an exact zone axis is reachable only by coincidence.
    """

    DOUBLE_TILT = "double_tilt"
    TILT_ROTATE = "tilt_rotate"
    SINGLE_TILT = "single_tilt"


@dataclass(frozen=True, slots=True)
class StagePosition:
    """A holder position: the two angles an operator reads off the stage.

    Purpose
    -------
    Gives the pair a name and a unit so that a tuple of bare floats can never be
    passed in the wrong order or the wrong unit. Angles are **degrees**, matching
    the repository-wide rule that only names ending in ``_rad`` carry radians.

    Attributes
    ----------
    alpha_deg : float
        Rotation about the holder rod.
    beta_deg : float
        Rotation about the cradle axis carried inside the rod. For a
        tilt-rotate holder this is the rotation about the holder normal.
    """

    alpha_deg: float
    beta_deg: float

    def __post_init__(self) -> None:
        for name in ("alpha_deg", "beta_deg"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"StagePosition.{name} must be finite; got {value!r}.")
            object.__setattr__(self, name, value)

    def as_tuple(self) -> tuple[float, float]:
        """The position as a plain ``(alpha_deg, beta_deg)`` tuple."""

        return (self.alpha_deg, self.beta_deg)

    def __str__(self) -> str:
        return f"(alpha={self.alpha_deg:+.2f} deg, beta={self.beta_deg:+.2f} deg)"


# --------------------------------------------------------------------------- #
# Tilt envelopes
# --------------------------------------------------------------------------- #


class TiltEnvelope:
    """The set of stage positions a holder can physically reach.

    Purpose
    -------
    A *predicate with a margin*, not a box. Real holders reduce one range as the
    other increases, grid bars shadow at mounting-dependent azimuths, and
    detectors intrude — so membership must be a callable and the distance to the
    boundary must be a continuous number, because both the solution ranking and
    the path planner need to know *how close* to a limit a position sits, not
    merely whether it is inside.

    Subclasses implement :meth:`contains` and :meth:`margin_deg`. The default
    :meth:`classify` turns those into the reachability verdict used by the
    navigation report.
    """

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        """Whether the position lies inside the envelope."""

        raise NotImplementedError

    def margin_deg(self, alpha_deg: float, beta_deg: float) -> float:
        """Signed degrees to the nearest boundary; negative when outside.

        A continuous quantity is required rather than a boolean because the
        ranking function applies a soft barrier to it and the path planner
        enforces a floor on it.
        """

        raise NotImplementedError

    def describe(self) -> str:
        """One-line prose description of the envelope."""

        raise NotImplementedError

    def bounds(self) -> tuple[float, float, float, float]:
        """A bounding box ``(alpha_min, alpha_max, beta_min, beta_max)``.

        Used for sampling and for drawing; need not be tight, but must contain
        the envelope.
        """

        raise NotImplementedError

    def accessible_solid_angle_sr(self, samples: int = 512) -> float:
        """Solid angle of beam directions the envelope reaches, in steradians.

        Purpose
        -------
        Quantifies how much of orientation space a holder actually commands. For
        a rectangular envelope the integral of the Jacobian ``cos(alpha)`` is
        analytic,

            Omega = (beta_max - beta_min) * (sin alpha_max - sin alpha_min),

        and :class:`RectangularEnvelope` overrides this with that closed form.
        The general case integrates the same Jacobian numerically over the
        bounding box, counting only points inside the envelope.

        Returns
        -------
        float
            Steradians. Divide by ``4 * pi`` for the sphere fraction; a zone axis
            and its reverse give the same pattern, so the *usable* fraction is
            twice that.
        """

        alpha_min, alpha_max, beta_min, beta_max = self.bounds()
        alpha_edges = np.linspace(alpha_min, alpha_max, samples + 1)
        beta_edges = np.linspace(beta_min, beta_max, samples + 1)
        alpha_mid = 0.5 * (alpha_edges[:-1] + alpha_edges[1:])
        beta_mid = 0.5 * (beta_edges[:-1] + beta_edges[1:])
        grid_alpha, grid_beta = np.meshgrid(alpha_mid, beta_mid, indexing="ij")
        inside = np.fromiter(
            (
                self.contains(float(a), float(b))
                for a, b in zip(grid_alpha.ravel(), grid_beta.ravel(), strict=True)
            ),
            dtype=bool,
            count=grid_alpha.size,
        ).reshape(grid_alpha.shape)
        cell = math.radians(alpha_edges[1] - alpha_edges[0]) * math.radians(
            beta_edges[1] - beta_edges[0]
        )
        return float(np.sum(np.cos(np.deg2rad(grid_alpha))[inside]) * cell)


@dataclass(frozen=True, slots=True)
class RectangularEnvelope(TiltEnvelope):
    """Independent, possibly asymmetric limits on each axis.

    The simplest useful envelope and the right default when only a datasheet
    figure such as "+/-30 degrees, +/-25 degrees" is known. Asymmetric limits are
    supported because real holders frequently have them.

    Attributes
    ----------
    alpha_min_deg, alpha_max_deg : float
    beta_min_deg, beta_max_deg : float
    """

    alpha_min_deg: float = -30.0
    alpha_max_deg: float = 30.0
    beta_min_deg: float = -30.0
    beta_max_deg: float = 30.0

    def __post_init__(self) -> None:
        if self.alpha_min_deg >= self.alpha_max_deg:
            raise ValueError("RectangularEnvelope requires alpha_min_deg < alpha_max_deg.")
        if self.beta_min_deg >= self.beta_max_deg:
            raise ValueError("RectangularEnvelope requires beta_min_deg < beta_max_deg.")

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        return bool(
            self.alpha_min_deg <= alpha_deg <= self.alpha_max_deg
            and self.beta_min_deg <= beta_deg <= self.beta_max_deg
        )

    def margin_deg(self, alpha_deg: float, beta_deg: float) -> float:
        return float(
            min(
                alpha_deg - self.alpha_min_deg,
                self.alpha_max_deg - alpha_deg,
                beta_deg - self.beta_min_deg,
                self.beta_max_deg - beta_deg,
            )
        )

    def bounds(self) -> tuple[float, float, float, float]:
        return (
            self.alpha_min_deg,
            self.alpha_max_deg,
            self.beta_min_deg,
            self.beta_max_deg,
        )

    def accessible_solid_angle_sr(self, samples: int = 512) -> float:
        """Closed form: ``(beta_max - beta_min)(sin alpha_max - sin alpha_min)``."""

        return float(
            math.radians(self.beta_max_deg - self.beta_min_deg)
            * (
                math.sin(math.radians(self.alpha_max_deg))
                - math.sin(math.radians(self.alpha_min_deg))
            )
        )

    def describe(self) -> str:
        return (
            f"rectangular envelope alpha in [{self.alpha_min_deg:+.1f}, "
            f"{self.alpha_max_deg:+.1f}] deg, beta in [{self.beta_min_deg:+.1f}, "
            f"{self.beta_max_deg:+.1f}] deg"
        )


@dataclass(frozen=True, slots=True)
class EllipticalEnvelope(TiltEnvelope):
    """A coupled envelope in which each range shrinks as the other grows.

    Membership is ``(alpha/alpha_max)^2 + (beta/beta_max)^2 <= 1``. This is a
    good fit to many double-tilt cartridges, where reaching the full alpha limit
    leaves no beta travel at all — a coupling a rectangular model misses, and one
    that changes which solutions are reachable near the corners.

    Attributes
    ----------
    alpha_max_deg, beta_max_deg : float
        Semi-axes of the ellipse; the limits reached when the other angle is zero.
    """

    alpha_max_deg: float = 30.0
    beta_max_deg: float = 30.0

    def __post_init__(self) -> None:
        if self.alpha_max_deg <= 0.0 or self.beta_max_deg <= 0.0:
            raise ValueError("EllipticalEnvelope semi-axes must be positive.")

    def _radius(self, alpha_deg: float, beta_deg: float) -> float:
        return math.hypot(alpha_deg / self.alpha_max_deg, beta_deg / self.beta_max_deg)

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        return bool(self._radius(alpha_deg, beta_deg) <= 1.0)

    def margin_deg(self, alpha_deg: float, beta_deg: float) -> float:
        """Approximate degrees to the boundary along the radial direction.

        Scales the normalized radial deficit by the local semi-axis, which is
        exact on the axes and a close approximation between them — adequate for a
        soft barrier and a path-margin floor, and monotone everywhere, which is
        what those consumers actually require.
        """

        radius = self._radius(alpha_deg, beta_deg)
        if radius == 0.0:
            return float(min(self.alpha_max_deg, self.beta_max_deg))
        scale = math.hypot(alpha_deg, beta_deg) / radius
        return float((1.0 - radius) * scale)

    def bounds(self) -> tuple[float, float, float, float]:
        return (
            -self.alpha_max_deg,
            self.alpha_max_deg,
            -self.beta_max_deg,
            self.beta_max_deg,
        )

    def describe(self) -> str:
        return (
            f"elliptical envelope with semi-axes alpha {self.alpha_max_deg:.1f} deg, "
            f"beta {self.beta_max_deg:.1f} deg (each range shrinks as the other grows)"
        )


@dataclass(frozen=True, slots=True)
class PolygonEnvelope(TiltEnvelope):
    """An arbitrary envelope digitized from a datasheet or measured directly.

    Purpose
    -------
    Holder envelopes are frequently neither rectangles nor ellipses — a
    tomography holder clipped by a detector, or a cartridge whose beta range
    collapses only on one side. This accepts the measured boundary as a closed
    polygon in ``(alpha, beta)`` degrees.

    Attributes
    ----------
    vertices : np.ndarray
        Shape ``(n, 2)``, columns ``(alpha_deg, beta_deg)``, in order. The
        polygon is closed implicitly; repeating the first vertex is harmless.
    """

    vertices: np.ndarray

    def __post_init__(self) -> None:
        array = as_float_array(self.vertices, shape=(None, 2))
        if array.shape[0] < 3:
            raise ValueError("PolygonEnvelope needs at least three vertices.")
        object.__setattr__(self, "vertices", array)

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        """Even-odd ray casting against the closed polygon."""

        vertices = self.vertices
        inside = False
        count = vertices.shape[0]
        for index in range(count):
            ax, ay = vertices[index]
            bx, by = vertices[(index + 1) % count]
            if (ay > beta_deg) != (by > beta_deg):
                crossing = ax + (beta_deg - ay) * (bx - ax) / (by - ay)
                if alpha_deg < crossing:
                    inside = not inside
        return inside

    def margin_deg(self, alpha_deg: float, beta_deg: float) -> float:
        """Distance to the nearest polygon edge, signed by :meth:`contains`."""

        point = np.array([alpha_deg, beta_deg], dtype=np.float64)
        starts = self.vertices
        ends = np.roll(self.vertices, -1, axis=0)
        segments = ends - starts
        lengths_sq = np.einsum("ij,ij->i", segments, segments)
        lengths_sq = np.where(lengths_sq < 1e-15, 1.0, lengths_sq)
        t = np.clip(
            np.einsum("ij,ij->i", point - starts, segments) / lengths_sq, 0.0, 1.0
        )
        closest = starts + t[:, None] * segments
        distance = float(np.min(np.linalg.norm(point - closest, axis=1)))
        return distance if self.contains(alpha_deg, beta_deg) else -distance

    def bounds(self) -> tuple[float, float, float, float]:
        return (
            float(np.min(self.vertices[:, 0])),
            float(np.max(self.vertices[:, 0])),
            float(np.min(self.vertices[:, 1])),
            float(np.max(self.vertices[:, 1])),
        )

    def describe(self) -> str:
        alpha_min, alpha_max, beta_min, beta_max = self.bounds()
        return (
            f"polygonal envelope with {self.vertices.shape[0]} vertices, bounded by "
            f"alpha in [{alpha_min:+.1f}, {alpha_max:+.1f}] deg, "
            f"beta in [{beta_min:+.1f}, {beta_max:+.1f}] deg"
        )


@dataclass(frozen=True, slots=True)
class MaskedEnvelope(TiltEnvelope):
    """A base envelope minus operator-marked exclusion regions.

    Purpose
    -------
    Grid-bar shadowing, detector intrusion and known-unsafe pockets are specific
    to a mounting, not to a holder, so they are modelled as *subtraction* from a
    holder envelope rather than baked into it. The excluded regions are
    themselves envelopes, so any shape available above can be used as a mask.

    Attributes
    ----------
    base : TiltEnvelope
    excluded : tuple of TiltEnvelope
    reason : str
        What the exclusions represent, carried into :meth:`describe` so a report
        can say *why* a position was refused.
    """

    base: TiltEnvelope
    excluded: tuple[TiltEnvelope, ...] = ()
    reason: str = "operator-marked exclusion"

    def __post_init__(self) -> None:
        object.__setattr__(self, "excluded", tuple(self.excluded))

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        if not self.base.contains(alpha_deg, beta_deg):
            return False
        return not any(region.contains(alpha_deg, beta_deg) for region in self.excluded)

    def margin_deg(self, alpha_deg: float, beta_deg: float) -> float:
        margin = self.base.margin_deg(alpha_deg, beta_deg)
        for region in self.excluded:
            # Inside an exclusion the margin is negative; outside it, the
            # distance to the exclusion boundary also limits how far the
            # position sits from an unusable region.
            margin = min(margin, -region.margin_deg(alpha_deg, beta_deg))
        return float(margin)

    def bounds(self) -> tuple[float, float, float, float]:
        return self.base.bounds()

    def describe(self) -> str:
        if not self.excluded:
            return self.base.describe()
        return (
            f"{self.base.describe()}, less {len(self.excluded)} excluded region(s) "
            f"({self.reason})"
        )


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class GeneralStageAxes:
    """Measured deviations of the tilt axes from their nominal directions.

    Purpose
    -------
    Carries the geometric part of a stage calibration: where the axes really
    point, whether they are orthogonal, and how much one drags the other. The
    default instance is the ideal stage, so a caller who has measured nothing
    gets textbook geometry and a clearly-labelled assumption rather than a
    silently wrong correction.

    A note on adopting fitted values: a spuriously fitted non-orthogonality is
    worse than an assumed-orthogonal stage, because it looks like knowledge.
    Adopt a fitted ``non_orthogonality_deg`` only when it is more than two
    standard errors from zero.

    Attributes
    ----------
    alpha_axis : np.ndarray
        Unit vector of the alpha axis in **laboratory** coordinates. Nominally
        ``(1, 0, 0)``.
    beta_axis : np.ndarray
        Unit vector of the beta axis in the **holder** frame at ``alpha = 0``.
        Nominally ``(0, 1, 0)``.
    coupling : float
        Linear mechanical coupling: the beta angle actually applied is
        ``beta + coupling * alpha``. Dimensionless.
    """

    alpha_axis: np.ndarray = field(
        default_factory=lambda: np.array([1.0, 0.0, 0.0], dtype=np.float64)
    )
    beta_axis: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 1.0, 0.0], dtype=np.float64)
    )
    coupling: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "alpha_axis", as_float_array(normalize_vector(self.alpha_axis), shape=(3,))
        )
        object.__setattr__(
            self, "beta_axis", as_float_array(normalize_vector(self.beta_axis), shape=(3,))
        )
        if not math.isfinite(self.coupling):
            raise ValueError("GeneralStageAxes.coupling must be finite.")

    @property
    def non_orthogonality_deg(self) -> float:
        """Departure of the axis pair from perpendicular, in degrees.

        Zero for an ideal stage. The sign follows the sign of the dot product,
        so a report can state which way the pair leans.
        """

        dot = float(np.dot(self.alpha_axis, self.beta_axis))
        return float(math.degrees(math.asin(max(-1.0, min(1.0, dot)))))

    @property
    def is_ideal(self) -> bool:
        """Whether this is the textbook orthogonal, uncoupled, aligned stage."""

        return bool(
            np.allclose(self.alpha_axis, [1.0, 0.0, 0.0], atol=1e-12)
            and np.allclose(self.beta_axis, [0.0, 1.0, 0.0], atol=1e-12)
            and abs(self.coupling) < 1e-12
        )


@dataclass(frozen=True, slots=True)
class StageCalibration:
    """Everything measured about a stage that the ideal model does not supply.

    Purpose
    -------
    Separates *what the instrument reports* from *what the geometry does*. The
    readouts are metadata; the transformation between them and the diffraction
    pattern is not, and essentially none of it is stored by vendor software. This
    record is the persistent result of the calibration procedures in
    `pytex.tem.calibration`.

    On ``diffraction_rotation_deg``: it depends on the *history* of the projector
    and diffraction lenses, not only on the nominal camera length, so the camera
    length and voltage it was measured at travel with it and
    :meth:`check_applicable` refuses a mismatch rather than interpolating. A
    plausible interpolated value for a hysteretic quantity is precisely the
    failure that sends an operator tilting the wrong way.

    Attributes
    ----------
    axes : GeneralStageAxes
    alpha_sign, beta_sign : int
        ``+1`` or ``-1``. Vendors and software versions differ; never assume.
    alpha_zero_deg, beta_zero_deg : float
        Readout offsets, subtracted before the sign is applied.
    diffraction_rotation_deg : float or None
        Azimuth of the pattern frame relative to the laboratory frame. ``None``
        means "not calibrated", which is a legitimate and common state — the
        two-zone reconstruction path does not need it.
    pattern_is_mirrored : bool
        Whether the stored pattern array is a mirrored rendering of the physical
        pattern. A wrong value makes the reconstructed orientation improper,
        which is detected rather than absorbed.
    camera_length_mm, accelerating_voltage_kv : float or None
        The conditions ``diffraction_rotation_deg`` was measured at.
    backlash_deg : float
        Measured repeatability difference between approaching a position from
        opposite directions. Not correctable; used to size the approach
        overshoot and to floor the reported uncertainty.
    angular_uncertainty_deg : float
        One-sigma stage readout uncertainty.
    notes : tuple of str
    provenance : ProvenanceRecord or None
    """

    axes: GeneralStageAxes = field(default_factory=GeneralStageAxes)
    alpha_sign: int = 1
    beta_sign: int = 1
    alpha_zero_deg: float = 0.0
    beta_zero_deg: float = 0.0
    diffraction_rotation_deg: float | None = None
    pattern_is_mirrored: bool = False
    camera_length_mm: float | None = None
    accelerating_voltage_kv: float | None = None
    backlash_deg: float = 0.0
    angular_uncertainty_deg: float = 0.1
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.alpha_sign not in (1, -1) or self.beta_sign not in (1, -1):
            raise ValueError("StageCalibration axis signs must be +1 or -1.")
        if self.backlash_deg < 0.0:
            raise ValueError("StageCalibration.backlash_deg must be non-negative.")
        if self.angular_uncertainty_deg < 0.0:
            raise ValueError(
                "StageCalibration.angular_uncertainty_deg must be non-negative."
            )
        object.__setattr__(self, "notes", tuple(self.notes))

    @property
    def is_rotation_calibrated(self) -> bool:
        """Whether a diffraction rotation is available for single-pattern work."""

        return self.diffraction_rotation_deg is not None

    def check_applicable(
        self,
        *,
        camera_length_mm: float | None = None,
        accelerating_voltage_kv: float | None = None,
        relative_tolerance: float = 0.01,
    ) -> None:
        """Raise when this calibration does not apply to the stated conditions.

        Deliberately strict. The diffraction rotation is hysteretic in the lens
        settings, so applying a value measured at another camera length is not a
        small extrapolation but a different number.
        """

        if camera_length_mm is not None and self.camera_length_mm is not None:
            if not math.isclose(
                camera_length_mm, self.camera_length_mm, rel_tol=relative_tolerance
            ):
                raise ValueError(
                    "This StageCalibration was measured at a camera length of "
                    f"{self.camera_length_mm:g} mm and cannot be applied at "
                    f"{camera_length_mm:g} mm. The diffraction rotation is hysteretic "
                    "in the lens settings; re-run the two-excursion calibration at "
                    "the camera length in use."
                )
        if (
            accelerating_voltage_kv is not None
            and self.accelerating_voltage_kv is not None
            and not math.isclose(
                accelerating_voltage_kv,
                self.accelerating_voltage_kv,
                rel_tol=relative_tolerance,
            )
        ):
            raise ValueError(
                "This StageCalibration was measured at "
                f"{self.accelerating_voltage_kv:g} kV and cannot be applied at "
                f"{accelerating_voltage_kv:g} kV."
            )

    def describe(self) -> str:
        """Convention-explicit prose summary of the calibration state."""

        parts: list[str] = []
        if self.axes.is_ideal:
            parts.append(
                "Stage axes assumed ideal: alpha along the laboratory x axis, beta "
                "orthogonal to it in the holder, no mechanical coupling."
            )
        else:
            parts.append(
                f"Stage axes measured: non-orthogonality "
                f"{self.axes.non_orthogonality_deg:+.2f} deg, alpha-to-beta coupling "
                f"{self.axes.coupling:+.4f}."
            )
        parts.append(
            f"Readout convention alpha x{self.alpha_sign:+d} offset "
            f"{self.alpha_zero_deg:+.2f} deg, beta x{self.beta_sign:+d} offset "
            f"{self.beta_zero_deg:+.2f} deg."
        )
        if self.is_rotation_calibrated:
            condition = ""
            if self.camera_length_mm is not None:
                condition = f" at {self.camera_length_mm:g} mm camera length"
                if self.accelerating_voltage_kv is not None:
                    condition += f" and {self.accelerating_voltage_kv:g} kV"
            parts.append(
                f"Diffraction rotation {self.diffraction_rotation_deg:+.2f} deg"
                f"{condition}; stored pattern is "
                f"{'mirrored' if self.pattern_is_mirrored else 'not mirrored'}."
            )
        else:
            parts.append(
                "Diffraction rotation NOT calibrated. Single-pattern orientation "
                "reconstruction is therefore unavailable; use the two-zone path, "
                "which does not need it."
            )
        parts.append(
            f"Backlash {self.backlash_deg:.2f} deg; readout uncertainty "
            f"{self.angular_uncertainty_deg:.2f} deg (one sigma)."
        )
        return " ".join(parts)

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload for the stage-calibration contract."""

        return {
            "alpha_axis": [float(value) for value in self.axes.alpha_axis],
            "beta_axis": [float(value) for value in self.axes.beta_axis],
            "coupling": float(self.axes.coupling),
            "alpha_sign": int(self.alpha_sign),
            "beta_sign": int(self.beta_sign),
            "alpha_zero_deg": float(self.alpha_zero_deg),
            "beta_zero_deg": float(self.beta_zero_deg),
            "diffraction_rotation_deg": (
                None
                if self.diffraction_rotation_deg is None
                else float(self.diffraction_rotation_deg)
            ),
            "pattern_is_mirrored": bool(self.pattern_is_mirrored),
            "camera_length_mm": (
                None if self.camera_length_mm is None else float(self.camera_length_mm)
            ),
            "accelerating_voltage_kv": (
                None
                if self.accelerating_voltage_kv is None
                else float(self.accelerating_voltage_kv)
            ),
            "backlash_deg": float(self.backlash_deg),
            "angular_uncertainty_deg": float(self.angular_uncertainty_deg),
            "notes": list(self.notes),
        }

    @classmethod
    def from_json_dict(cls, payload: dict[str, Any]) -> StageCalibration:
        """Rebuild a calibration from :meth:`to_json_dict` output."""

        return cls(
            axes=GeneralStageAxes(
                alpha_axis=np.asarray(
                    payload.get("alpha_axis", (1.0, 0.0, 0.0)), dtype=np.float64
                ),
                beta_axis=np.asarray(
                    payload.get("beta_axis", (0.0, 1.0, 0.0)), dtype=np.float64
                ),
                coupling=float(payload.get("coupling", 0.0)),
            ),
            alpha_sign=int(payload.get("alpha_sign", 1)),
            beta_sign=int(payload.get("beta_sign", 1)),
            alpha_zero_deg=float(payload.get("alpha_zero_deg", 0.0)),
            beta_zero_deg=float(payload.get("beta_zero_deg", 0.0)),
            diffraction_rotation_deg=(
                None
                if payload.get("diffraction_rotation_deg") is None
                else float(payload["diffraction_rotation_deg"])
            ),
            pattern_is_mirrored=bool(payload.get("pattern_is_mirrored", False)),
            camera_length_mm=(
                None
                if payload.get("camera_length_mm") is None
                else float(payload["camera_length_mm"])
            ),
            accelerating_voltage_kv=(
                None
                if payload.get("accelerating_voltage_kv") is None
                else float(payload["accelerating_voltage_kv"])
            ),
            backlash_deg=float(payload.get("backlash_deg", 0.0)),
            angular_uncertainty_deg=float(payload.get("angular_uncertainty_deg", 0.1)),
            notes=tuple(payload.get("notes", ())),
        )


# --------------------------------------------------------------------------- #
# Stage models
# --------------------------------------------------------------------------- #


@runtime_checkable
class StageModel(Protocol):
    """The forward kinematics of a specimen holder.

    Purpose
    -------
    One small interface that every holder type satisfies, so that navigation,
    path planning and plotting are written once against "a stage" rather than
    against a family of ``if holder_type ==`` branches.

    Implementations must provide :meth:`rotation_matrix` — the holder-to-laboratory
    rotation at a stage position — plus the envelope, the calibration and the
    holder kind. The remaining behaviour is derived from those in
    :class:`_StageCommon`.
    """

    @property
    def kind(self) -> HolderKind: ...

    @property
    def envelope(self) -> TiltEnvelope: ...

    @property
    def calibration(self) -> StageCalibration: ...

    def rotation_matrix(self, alpha_deg: float, beta_deg: float) -> np.ndarray: ...

    def beam_direction(self, alpha_deg: float, beta_deg: float) -> np.ndarray: ...

    def describe(self) -> str: ...


@dataclass(frozen=True, slots=True)
class _StageCommon:
    """Shared state and derived behaviour for the concrete stage models."""

    envelope: TiltEnvelope = field(default_factory=RectangularEnvelope)
    calibration: StageCalibration = field(default_factory=StageCalibration)
    name: str = "stage"

    def _applied_angles_rad(self, alpha_deg: float, beta_deg: float) -> tuple[float, float]:
        """Readout angles turned into the angles the mechanism actually applies.

        Applies, in order: the readout zero offset, the vendor sign convention,
        and the linear alpha-to-beta coupling.
        """

        calibration = self.calibration
        alpha = calibration.alpha_sign * (alpha_deg - calibration.alpha_zero_deg)
        beta = calibration.beta_sign * (beta_deg - calibration.beta_zero_deg)
        beta = beta + calibration.axes.coupling * alpha
        return math.radians(alpha), math.radians(beta)

    def beam_direction(self, alpha_deg: float, beta_deg: float) -> np.ndarray:
        """Beam direction in **holder** coordinates at a stage position.

        This is the transpose action of :meth:`rotation_matrix` on the
        laboratory beam axis. For an ideal double-tilt stage it reduces to the
        closed form :func:`beam_direction_holder`; the general path is used
        whenever the calibration is non-ideal, so a caller never has to ask
        which model is in force.
        """

        matrix: np.ndarray = self.rotation_matrix(alpha_deg, beta_deg)  # type: ignore[attr-defined]
        return np.asarray(matrix.T @ BEAM_AXIS_LABORATORY, dtype=np.float64)

    def contains(self, alpha_deg: float, beta_deg: float) -> bool:
        """Whether a stage position lies inside the holder envelope."""

        return self.envelope.contains(alpha_deg, beta_deg)


@dataclass(frozen=True, slots=True)
class DoubleTiltStage(_StageCommon):
    """The standard double-tilt holder: alpha about the rod, beta about the cradle.

    Purpose
    -------
    The workhorse geometry, and the one the closed-form tilt solution is derived
    for. The rotation is ``Rx(alpha) Ry(beta)`` in the ideal case, obtained by
    composing the beta rotation about its *instantaneous* laboratory axis with
    the alpha rotation — the moving-axis factors cancel exactly.

    When to use
    -----------
    Any conventional double-tilt specimen holder. This is the default stage
    throughout `pytex.tem`.

    Expected inputs and outputs
    ---------------------------
    Angles in degrees; :meth:`rotation_matrix` returns the 3x3 holder-to-laboratory
    rotation, :meth:`beam_direction` the unit beam direction in holder
    coordinates.

    Examples
    --------
    At zero tilt the stage is the identity and the beam lies along holder ``z``::

        >>> stage = DoubleTiltStage()
        >>> bool(np.allclose(stage.rotation_matrix(0.0, 0.0), np.eye(3)))
        True

    A pure alpha tilt moves the beam direction toward the holder ``y`` axis,
    which is the sign convention the calibration procedure keys on::

        >>> float(np.round(stage.beam_direction(10.0, 0.0)[1], 6))
        0.173648

    See Also
    --------
    pytex.tem.navigation.plan_tilt_to_zone_axis : inverts this map.
    """

    @property
    def kind(self) -> HolderKind:
        return HolderKind.DOUBLE_TILT

    def rotation_matrix(self, alpha_deg: float, beta_deg: float) -> np.ndarray:
        """Holder-to-laboratory rotation at ``(alpha_deg, beta_deg)``.

        Uses the ideal closed form when the axes are ideal, and otherwise
        composes rotations about the calibrated axes explicitly — including the
        motion of the beta axis with alpha, which the ideal cancellation would
        otherwise hide.
        """

        alpha_rad, beta_rad = self._applied_angles_rad(alpha_deg, beta_deg)
        axes = self.calibration.axes
        if axes.is_ideal:
            return np.asarray(rotation_x(alpha_rad) @ rotation_y(beta_rad), dtype=np.float64)
        alpha_rotation = _rotation_about(axes.alpha_axis, alpha_rad)
        beta_rotation = _rotation_about(axes.beta_axis, beta_rad)
        return np.asarray(alpha_rotation @ beta_rotation, dtype=np.float64)

    def describe(self) -> str:
        """Prose statement of the geometry, conventions and limits."""

        return (
            f"Double-tilt holder '{self.name}': alpha rotates the rod about the "
            "laboratory x axis, beta rotates a cradle carried inside the rod, so the "
            "holder-to-laboratory rotation is Rx(alpha) Ry(beta). The beam direction "
            "in holder coordinates is (-cos a sin b, sin a, cos a cos b), which is a "
            "spherical coordinate system whose pole is the beta axis. Reachable set: "
            f"{self.envelope.describe()}, spanning "
            f"{self.envelope.accessible_solid_angle_sr():.3f} sr "
            f"({100.0 * self.envelope.accessible_solid_angle_sr() / (4.0 * math.pi):.1f}% "
            "of all beam directions, or twice that counting a zone axis and its "
            f"reverse as equivalent). {self.calibration.describe()}"
        )


@dataclass(frozen=True, slots=True)
class TiltRotateStage(_StageCommon):
    """A tilt-rotate holder: alpha about the rod, then rotation about the holder normal.

    Purpose
    -------
    The second angle is an in-plane rotation of the specimen rather than a
    second tilt, so the reachable set has a different shape: rotating about the
    holder normal sweeps a small circle about that normal, and alpha carries it.
    The navigation algebra is unchanged — only :meth:`rotation_matrix` differs —
    which is the point of routing everything through the stage interface.

    When to use
    -----------
    Rotation holders, and tomography holders whose second freedom is azimuthal.
    """

    @property
    def kind(self) -> HolderKind:
        return HolderKind.TILT_ROTATE

    def rotation_matrix(self, alpha_deg: float, beta_deg: float) -> np.ndarray:
        """Holder-to-laboratory rotation ``Rx(alpha) Rz(theta)``."""

        alpha_rad, beta_rad = self._applied_angles_rad(alpha_deg, beta_deg)
        axes = self.calibration.axes
        alpha_rotation = (
            rotation_x(alpha_rad)
            if axes.is_ideal
            else _rotation_about(axes.alpha_axis, alpha_rad)
        )
        return np.asarray(alpha_rotation @ rotation_z(beta_rad), dtype=np.float64)

    def describe(self) -> str:
        return (
            f"Tilt-rotate holder '{self.name}': alpha tilts about the laboratory x "
            "axis and the second angle rotates the specimen about its own normal, "
            "giving Rx(alpha) Rz(theta). Reachable set: "
            f"{self.envelope.describe()}. {self.calibration.describe()}"
        )


@dataclass(frozen=True, slots=True)
class SingleTiltStage(_StageCommon):
    """A single-tilt holder: one axis, and a reachable set of measure zero.

    Purpose
    -------
    Models the common case honestly. With one freedom the beam direction traces
    a **curve** on the sphere, not a region, so an exact zone axis is reachable
    only by coincidence. The correct answer for almost any target is therefore
    "nearest approach", and the navigation report says so rather than reporting a
    failure — a qualitative difference from the double-tilt case that the API
    must express, not round away.

    The ``beta`` argument is accepted and ignored, so the interface stays uniform.
    """

    @property
    def kind(self) -> HolderKind:
        return HolderKind.SINGLE_TILT

    def rotation_matrix(self, alpha_deg: float, beta_deg: float = 0.0) -> np.ndarray:
        """Holder-to-laboratory rotation ``Rx(alpha)``; ``beta_deg`` is ignored."""

        alpha_rad, _ = self._applied_angles_rad(alpha_deg, 0.0)
        axes = self.calibration.axes
        if axes.is_ideal:
            return rotation_x(alpha_rad)
        return _rotation_about(axes.alpha_axis, alpha_rad)

    def describe(self) -> str:
        return (
            f"Single-tilt holder '{self.name}': one axis only, so reachable beam "
            "directions form a great circle — a set of measure zero. An exact zone "
            "axis is reachable only by coincidence; expect nearest-approach "
            f"solutions. Reachable set: {self.envelope.describe()}. "
            f"{self.calibration.describe()}"
        )
