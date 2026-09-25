"""The five frames of FIB lamella planning, and the conversions between them.

Purpose
-------
A lamella is planned in one instrument and read in another. The orientation
comes from an EBSD scan, the rectangle is placed on an SEM image, the milling
pattern is entered in the ion-column view of a tilted stage, and the answer is
judged on a TEM holder. Each of those speaks in its own axes, and an error in
the bookkeeping between them is invisible until a real lamella is milled at the
wrong azimuth. This module states every one of those relationships once, as a
typed and user-editable object, so no other part of :mod:`pytex.fib` holds a
private convention.

The frames
----------
``C`` crystal
    The phase's crystal frame, exactly as :mod:`pytex.core.lattice` defines it.
``E`` EBSD specimen frame
    The frame the orientation file's Euler angles refer to. PyTex's readers
    import vendor angles into it without remapping
    (:mod:`pytex.adapters.scan_files`), and the file is assumed to carry the
    70 degree tilt correction already (decision D1 of the foundation document);
    nothing here applies a second one.
``S`` FIB sample frame
    ``X_s, Y_s`` in the surface, ``Z_s`` the **outward** surface normal. Built
    from ``E`` by :class:`SurfaceGeometry`, which is where a vendor's "z points
    into the sample" convention is declared rather than assumed.
``I`` SEM image
    The 2D pixel frame of the SEM or image-quality picture the rectangle is
    drawn on. Related to the EBSD scan coordinates by :class:`ImageRegistration`:
    identity by default (decision D3), an optional similarity, or an affine fit
    to control points with its residual reported.
``P`` FIB ion view
    The plan view the ion column sees when the stage is tilted to the column
    angle. Related to ``S`` by :class:`ChamberGeometry`, whose rotation sense
    and offset are **calibrated, never asserted** (section 5.1).
``L`` lamella
    ``n_L`` the lamella-plane normal (in the surface), ``t_L`` the long axis,
    ``Z_s`` the depth axis; the ion beam mills along ``-Z_s``.
``H`` TEM holder
    :data:`pytex.tem.reconstruction.HOLDER_FRAME`, reached from ``L`` by the
    unknown mounting rotation of :class:`pytex.fib.planning.MountModel`.

The rigid three-dimensional part of the chain, ``C -> E -> S -> L -> H``, is
registered in a :class:`pytex.core.frames.FrameGraph` by
:func:`lamella_frame_graph`. The image and ion-view maps are two-dimensional and
may be non-rigid (an affine registration shears), so they are carried by the
dataclasses below rather than by a :class:`~pytex.core.frames.FrameTransform`.

No new frame domain is introduced (``docs/standards/notation_and_conventions.md``):
the sample and lamella frames are specimen-domain frames, the SEM image is a
detector-domain frame and the ion view a laboratory-domain frame.

See ``docs/architecture/fib_lamella_planning_foundation.md`` section 5.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core import frame_catalog
from pytex.core.frames import FrameGraph, FrameTransform, ReferenceFrame

__all__ = [
    "EBSD_SPECIMEN_FRAME",
    "FIB_SAMPLE_FRAME",
    "ION_VIEW_FRAME",
    "LAMELLA_FRAME",
    "SEM_IMAGE_FRAME",
    "SUPPORTED_COLUMN_ANGLES_DEG",
    "ChamberGeometry",
    "FiducialObservation",
    "ImageRegistration",
    "SurfaceGeometry",
    "calibrate_chamber_from_fiducials",
    "lamella_frame_graph",
    "wrap_azimuth_deg",
]

#: The EBSD specimen frame: whatever frame the scan's Euler angles refer to.
#: Identical to the frame the PyTex scan readers import into, so an imported
#: orientation's specimen frame compares equal to it.
EBSD_SPECIMEN_FRAME: ReferenceFrame = frame_catalog.specimen_frame()

#: The FIB sample frame ``S``: ``X_s, Y_s`` in the surface, ``Z_s`` outward.
FIB_SAMPLE_FRAME: ReferenceFrame = frame_catalog.specimen_frame(
    "fib_sample",
    axes=("X_s", "Y_s", "Z_s"),
    description=(
        "FIB sample frame: X_s and Y_s lie in the polished surface, Z_s is the outward surface "
        "normal. Built from the EBSD specimen frame by SurfaceGeometry; the ion beam of a "
        "strictly vertical mill travels along -Z_s."
    ),
)

#: The lamella frame ``L``: plane normal, long axis, depth axis.
LAMELLA_FRAME: ReferenceFrame = frame_catalog.specimen_frame(
    "lamella",
    axes=("n_L", "t_L", "Z_s"),
    description=(
        "Lamella frame: n_L is the lamella-plane normal and lies in the surface, t_L = Z_s x n_L "
        "is the long axis of the slab, and Z_s is the depth axis shared with the sample frame. "
        "(n_L, t_L, Z_s) is right-handed by construction."
    ),
)

#: The SEM (or image-quality) image plane ``I``: pixel columns, pixel rows, view axis.
SEM_IMAGE_FRAME: ReferenceFrame = frame_catalog.detector_frame(
    "sem_image",
    axes=("x_I", "y_I", "z_I"),
)

#: The FIB ion-column view ``P`` of the stage tilted to the column angle.
ION_VIEW_FRAME: ReferenceFrame = frame_catalog.laboratory_frame(
    "fib_ion_view",
    axes=("x_P", "y_P", "z_P"),
    description=(
        "Plan view seen by the ion column when the stage is tilted to the column angle T. At that "
        "tilt the ion view is undistorted; its in-plane rotation relative to the sample frame is "
        "the calibrated pair (s_R, R0) of ChamberGeometry."
    ),
)

#: Ion-column angles the workbench offers (decision D4). The library accepts any
#: angle strictly between 0 and 90 degrees; these are the two instrument
#: families the maintainer named, 54 being the default.
SUPPORTED_COLUMN_ANGLES_DEG: tuple[float, ...] = (52.0, 54.0)


def wrap_azimuth_deg(angle_deg: ArrayLike, *, period: float = 360.0) -> np.ndarray:
    """Wrap azimuths into ``[0, period)``.

    Parameters
    ----------
    angle_deg : array_like
        Angles in degrees.
    period : float, default 360
        ``180`` for the azimuth of an axis (a plane normal and its reverse name
        the same plane), ``360`` for a directed azimuth.

    Returns
    -------
    numpy.ndarray
        Same shape as the input.

    Examples
    --------
    >>> float(wrap_azimuth_deg(-30.0))
    330.0
    >>> float(wrap_azimuth_deg(200.0, period=180.0))
    20.0
    """

    values = np.mod(np.asarray(angle_deg, dtype=np.float64), period)
    # np.mod can return `period` itself for tiny negative inputs.
    return np.where(np.isclose(values, period, rtol=0.0, atol=1e-12), 0.0, values)


def _sign(value: int, name: str) -> int:
    if int(value) not in (-1, 1):
        raise ValueError(f"{name} must be +1 or -1; got {value!r}.")
    return int(value)


@dataclass(frozen=True, slots=True)
class SurfaceGeometry:
    """Where the polished surface lies in the EBSD specimen frame.

    Purpose
    -------
    Declares the two conventions every EBSD vendor answers differently and no
    file header states: which sense of the specimen ``z`` axis is the
    **outward** surface normal, and which sense of the specimen ``y`` axis the
    scan's row coordinate increases along. With both declared, the sample frame
    ``S`` is fixed and right-handed, and every angle downstream is unambiguous.

    When to use
    -----------
    Always, if only to accept the defaults explicitly. Change ``normal_sign``
    when the scan software's specimen ``z`` points *into* the material, and
    ``scan_y_sign`` when the scan rows run opposite to the specimen ``y`` axis.
    Check both against a fiducial the first time a new instrument is used: an
    error here mirrors every azimuth and is invisible on a symmetric grain.

    Attributes
    ----------
    normal_sign : {+1, -1}, default +1
        ``+1``: the outward normal ``Z_s`` is ``+z`` of the EBSD specimen frame.
        ``-1``: it is ``-z``.
    scan_y_sign : {+1, -1}, default +1
        ``+1``: the scan's ``y`` coordinate increases along ``+y`` of the EBSD
        specimen frame. ``-1``: along ``-y``.

    Notes
    -----
    ``X_s`` is taken as ``+x`` of the EBSD specimen frame and
    ``Y_s = Z_s x X_s`` completes a right-handed set, so ``normal_sign = -1``
    also reverses ``Y_s``. The scan ``x`` coordinate always runs along ``X_s``.

    Examples
    --------
    >>> SurfaceGeometry().specimen_to_sample_matrix().tolist()
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    >>> SurfaceGeometry(normal_sign=-1).specimen_to_sample_matrix()[2].tolist()
    [0.0, 0.0, -1.0]
    """

    normal_sign: int = 1
    scan_y_sign: int = 1

    def __post_init__(self) -> None:
        object.__setattr__(self, "normal_sign", _sign(self.normal_sign, "normal_sign"))
        object.__setattr__(self, "scan_y_sign", _sign(self.scan_y_sign, "scan_y_sign"))

    def specimen_to_sample_matrix(self) -> np.ndarray:
        """Rotation taking EBSD-specimen components to sample components.

        Rows are ``X_s``, ``Y_s``, ``Z_s`` written in the EBSD specimen frame,
        so ``v_S = M v_E``. Proper (determinant +1) for every setting.
        """

        s = float(self.normal_sign)
        return np.array([[1.0, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]], dtype=np.float64)

    def frame_transform(self) -> FrameTransform:
        """The ``E -> S`` transform, for a :class:`~pytex.core.frames.FrameGraph`."""

        return FrameTransform(
            source=EBSD_SPECIMEN_FRAME,
            target=FIB_SAMPLE_FRAME,
            rotation_matrix=self.specimen_to_sample_matrix(),
        )

    @property
    def scan_to_sample_y(self) -> float:
        """Factor taking a scan ``y`` coordinate to a sample ``Y_s`` coordinate."""

        return float(self.normal_sign * self.scan_y_sign)

    def scan_to_sample_xy(self, points: ArrayLike) -> np.ndarray:
        """Convert ``(..., 2)`` scan coordinates to in-plane sample coordinates."""

        array = np.asarray(points, dtype=np.float64)
        out = array.copy()
        out[..., 1] = array[..., 1] * self.scan_to_sample_y
        return out

    def sample_to_scan_xy(self, points: ArrayLike) -> np.ndarray:
        """Inverse of :meth:`scan_to_sample_xy` (the map is its own inverse)."""

        return self.scan_to_sample_xy(points)

    def describe(self) -> str:
        """One sentence stating both conventions."""

        normal = "+z" if self.normal_sign > 0 else "-z"
        rows = "+y" if self.scan_y_sign > 0 else "-y"
        return (
            f"Surface: the outward normal Z_s is {normal} of the EBSD specimen frame, "
            "X_s is its +x, "
            f"Y_s = Z_s x X_s; scan rows increase along {rows} of the specimen frame, so a scan y "
            f"coordinate maps to Y_s with factor {self.scan_to_sample_y:+.0f}."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "normal_sign": self.normal_sign,
            "scan_y_sign": self.scan_y_sign,
            "scan_to_sample_y": self.scan_to_sample_y,
            "specimen_to_sample_matrix": self.specimen_to_sample_matrix().tolist(),
        }


@dataclass(frozen=True, slots=True)
class ImageRegistration:
    """How EBSD scan coordinates land on the SEM image.

    Purpose
    -------
    The rectangle is planned from the EBSD map but milled against an SEM (or
    ion) image, so the two must be registered. The default is the identity of
    decision D3 --- image ``x`` along scan ``x``, image ``y`` along scan ``y``,
    one image unit per micrometre, no flip --- and every departure from it is a
    declared parameter rather than a correction applied in someone's head.

    When to use
    -----------
    Leave the default when the SEM image is the scan's own image-quality map.
    Use :meth:`similarity` when the rotation, scale or handedness between the
    two is known, and :meth:`from_control_points` when three or more features
    have been located in both, which fits a full affine map and reports how
    well it fits.

    Attributes
    ----------
    matrix : (2, 2) array
        Linear part ``A`` of ``q = A p + t``, ``p`` in scan micrometres and
        ``q`` in image units.
    translation : (2,) array
        Offset ``t`` in image units.
    residual_um : float
        RMS misfit of the control points, mapped back to scan micrometres.
        Zero for a declared registration.
    residual_deg : float
        The angular uncertainty that misfit implies for an azimuth read
        through the registration: the RMS residual divided by the RMS radius
        of the control points about their centroid. Enters the uncertainty
        budget of :mod:`pytex.fib.planning`. Zero for a declared registration.
    control_points : int
        How many point pairs were fitted; zero for a declared registration.
    method : str
        ``"identity"``, ``"similarity"`` or ``"affine"``.
    image_units : str
        What an image unit is, for reporting (``"um"`` or ``"px"``).

    Examples
    --------
    The default is a proven no-op:

    >>> reg = ImageRegistration.identity()
    >>> reg.to_image([[3.0, 4.0]]).tolist()
    [[3.0, 4.0]]
    >>> reg.image_azimuth_deg(30.0, SurfaceGeometry())
    30.0
    """

    matrix: np.ndarray = field(default_factory=lambda: np.eye(2))
    translation: np.ndarray = field(default_factory=lambda: np.zeros(2))
    residual_um: float = 0.0
    residual_deg: float = 0.0
    control_points: int = 0
    method: str = "identity"
    image_units: str = "um"

    def __post_init__(self) -> None:
        matrix = np.asarray(self.matrix, dtype=np.float64)
        translation = np.asarray(self.translation, dtype=np.float64)
        if matrix.shape != (2, 2) or not np.all(np.isfinite(matrix)):
            raise ValueError("ImageRegistration.matrix must be a finite 2x2 array.")
        if translation.shape != (2,) or not np.all(np.isfinite(translation)):
            raise ValueError("ImageRegistration.translation must be a finite 2-vector.")
        if abs(float(np.linalg.det(matrix))) < 1e-12:
            raise ValueError("ImageRegistration.matrix is singular; it cannot be inverted.")
        if self.residual_um < 0.0 or self.residual_deg < 0.0:
            raise ValueError("Registration residuals must be non-negative.")
        if self.method not in ("identity", "similarity", "affine"):
            raise ValueError("ImageRegistration.method must be identity, similarity or affine.")
        matrix = np.ascontiguousarray(matrix)
        matrix.setflags(write=False)
        translation = np.ascontiguousarray(translation)
        translation.setflags(write=False)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "translation", translation)

    # -- constructors --------------------------------------------------------

    @classmethod
    def identity(cls) -> ImageRegistration:
        """Decision D3: image axes along the scan axes, unit scale, no flip."""

        return cls()

    @classmethod
    def similarity(
        cls,
        *,
        rotation_deg: float = 0.0,
        scale: float = 1.0,
        flip_y: bool = False,
        translation: Sequence[float] = (0.0, 0.0),
        image_units: str = "um",
    ) -> ImageRegistration:
        """A declared rotation, scale and optional handedness flip.

        Parameters
        ----------
        rotation_deg : float
            Counter-clockwise rotation ``omega`` of the scan axes into the image
            axes, applied after the flip.
        scale : float
            Image units per scan micrometre; must be positive.
        flip_y : bool
            Reverse the scan ``y`` axis before rotating (a mirror image).
        translation : pair of float
            Image coordinates of the scan origin.
        image_units : str
            ``"um"`` or ``"px"``.
        """

        if not scale > 0.0:
            raise ValueError("ImageRegistration scale must be positive.")
        c, s = math.cos(math.radians(rotation_deg)), math.sin(math.radians(rotation_deg))
        flip = np.diag([1.0, -1.0 if flip_y else 1.0])
        matrix = float(scale) * np.array([[c, -s], [s, c]]) @ flip
        method = (
            "identity"
            if np.allclose(matrix, np.eye(2)) and not any(translation)
            else ("similarity")
        )
        return cls(
            matrix=matrix,
            translation=np.asarray(translation, dtype=np.float64),
            method=method,
            image_units=image_units,
        )

    @classmethod
    def from_control_points(
        cls,
        scan_points: ArrayLike,
        image_points: ArrayLike,
        *,
        image_units: str = "px",
    ) -> ImageRegistration:
        """Least-squares affine registration from matched features.

        Purpose
        -------
        The optional, user-activated registration of decision D3. Three
        non-collinear pairs determine an affine map exactly; more are fitted by
        least squares and the misfit is reported, because a registration that
        cannot say how well it fits is the one that ruins a real lamella.

        Parameters
        ----------
        scan_points : (n, 2) array_like
            Feature positions in EBSD scan micrometres.
        image_points : (n, 2) array_like
            The same features in image units.
        image_units : str, default "px"

        Returns
        -------
        ImageRegistration
            With ``residual_um`` and ``residual_deg`` filled in. For exactly
            three points both are zero: the fit is exact, not validated, and
            :meth:`describe` says so.

        Raises
        ------
        ValueError
            Fewer than three pairs, mismatched shapes, or collinear points.
        """

        scan = np.asarray(scan_points, dtype=np.float64)
        image = np.asarray(image_points, dtype=np.float64)
        if scan.ndim != 2 or scan.shape[1] != 2 or scan.shape != image.shape:
            raise ValueError("Control points must be two matching (n, 2) arrays.")
        count = int(scan.shape[0])
        if count < 3:
            raise ValueError("An affine registration needs at least three control-point pairs.")
        design = np.column_stack([scan, np.ones(count)])
        if np.linalg.matrix_rank(design, tol=1e-9 * max(1.0, float(np.abs(scan).max()))) < 3:
            raise ValueError("The control points are collinear; an affine map is undetermined.")
        solution, *_ = np.linalg.lstsq(design, image, rcond=None)
        matrix = solution[:2].T
        translation = solution[2]
        predicted = scan @ matrix.T + translation
        misfit_image = image - predicted
        misfit_scan = misfit_image @ np.linalg.inv(matrix).T
        dof = max(count - 3, 1)
        residual_um = float(np.sqrt(np.sum(misfit_scan**2) / dof)) if count > 3 else 0.0
        radius = float(np.sqrt(np.mean(np.sum((scan - scan.mean(axis=0)) ** 2, axis=1))))
        residual_deg = math.degrees(residual_um / radius) if radius > 0.0 else 0.0
        return cls(
            matrix=matrix,
            translation=translation,
            residual_um=residual_um,
            residual_deg=residual_deg,
            control_points=count,
            method="affine",
            image_units=image_units,
        )

    # -- mapping -------------------------------------------------------------

    def to_image(self, scan_points: ArrayLike) -> np.ndarray:
        """Map ``(..., 2)`` scan micrometres to image coordinates."""

        points = np.asarray(scan_points, dtype=np.float64)
        return np.asarray(points @ self.matrix.T + self.translation, dtype=np.float64)

    def to_scan(self, image_points: ArrayLike) -> np.ndarray:
        """Map ``(..., 2)`` image coordinates back to scan micrometres."""

        points = np.asarray(image_points, dtype=np.float64)
        inverse = np.linalg.inv(self.matrix)
        return np.asarray((points - self.translation) @ inverse.T, dtype=np.float64)

    def image_direction(self, theta_sample_deg: float, surface: SurfaceGeometry) -> np.ndarray:
        """Unit image-frame direction of a sample-frame azimuth."""

        theta = math.radians(theta_sample_deg)
        sample = np.array([math.cos(theta), math.sin(theta)])
        scan = surface.sample_to_scan_xy(sample)
        image = self.matrix @ scan
        return np.asarray(image / np.linalg.norm(image), dtype=np.float64)

    def image_azimuth_deg(self, theta_sample_deg: float, surface: SurfaceGeometry) -> float:
        """Azimuth, in the image frame, of a direction at ``theta_sample_deg`` in ``S``.

        Measured from image ``+x`` towards image ``+y`` and wrapped to
        ``[0, 360)``. Note that image ``+y`` is whatever the image's row axis
        is; on a screen where rows run downwards the angle reads clockwise.
        """

        direction = self.image_direction(theta_sample_deg, surface)
        angle = math.degrees(math.atan2(float(direction[1]), float(direction[0])))
        return round(float(wrap_azimuth_deg(angle)), 12)

    @property
    def is_mirror(self) -> bool:
        """Whether the registration reverses handedness."""

        return bool(np.linalg.det(self.matrix) < 0.0)

    @property
    def is_identity(self) -> bool:
        """Whether the registration is exactly the D3 default."""

        return bool(np.array_equal(self.matrix, np.eye(2)) and not np.any(self.translation))

    def describe(self) -> str:
        """The registration in one or two sentences, residual included."""

        if self.is_identity:
            return (
                "Registration: identity (decision D3) - image x along scan x, "
                "image y along scan y, "
                "one image unit per micrometre, no flip."
            )
        scale = math.sqrt(abs(float(np.linalg.det(self.matrix))))
        rotation = math.degrees(math.atan2(float(self.matrix[1, 0]), float(self.matrix[0, 0])))
        text = (
            f"Registration: {self.method}, mean scale {scale:.4g} {self.image_units} per um, "
            f"rotation {rotation:+.2f} deg{', mirror' if self.is_mirror else ''}."
        )
        if self.method == "affine":
            if self.control_points > 3:
                text += (
                    f" Fitted to {self.control_points} control points: RMS residual "
                    f"{self.residual_um:.3g} um, i.e. {self.residual_deg:.3g} deg on an azimuth."
                )
            else:
                text += (
                    " Fitted exactly to three control points, so the fit has no redundancy and "
                    "its error is unknown; add a fourth point to measure it."
                )
        return text

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "method": self.method,
            "matrix": self.matrix.tolist(),
            "translation": self.translation.tolist(),
            "is_mirror": self.is_mirror,
            "residual_um": float(self.residual_um),
            "residual_deg": float(self.residual_deg),
            "control_points": int(self.control_points),
            "image_units": self.image_units,
        }


@dataclass(frozen=True, slots=True)
class ChamberGeometry:
    """The FIB-SEM chamber: column angle, stage rotation, and the calibrated azimuth map.

    Purpose
    -------
    Converts the lamella azimuth found in the sample frame into the rotation to
    type into the FIB pattern box,

        theta_ion = s_R * (theta_S + R_stage) + R0,

    and records the one fact that makes that number trustworthy or not: whether
    ``s_R`` and ``R0`` were calibrated on this instrument with a fiducial. PyTex
    does not know the sign convention of any vendor's pattern rotation and does
    not pretend to (section 5.1): an uncalibrated chamber keeps the defaults
    ``s_R = +1``, ``R0 = 0`` and every report says, prominently, that the
    number must not be used to mill until the protocol in
    ``docs/site/workflows/fib_azimuth_calibration.md`` has been run.

    Attributes
    ----------
    column_angle_deg : float, default 54
        Angle ``T`` between the electron and ion columns (decision D4). The
        stage is tilted to ``T`` so the ion beam is normal to the surface.
    rotation_sense : {+1, -1}, default +1
        ``s_R``: +1 when a positive pattern rotation turns the box the same way
        as a positive sample azimuth, as seen in the ion view.
    rotation_offset_deg : float, default 0
        ``R0``: the pattern rotation that aligns the box with sample azimuth
        zero at stage rotation zero.
    stage_rotation_deg : float, default 0
        ``R_stage``: the stage rotation reading at which the lamella is milled.
    tilt_axis_azimuth_deg : float, default 0
        Azimuth of the stage tilt axis in the sample frame, used only for the
        foreshortening of the SEM view at tilt.
    compucentric : bool, default True
        Whether the stage rotates about the field of view. When false, a
        stage rotation also translates the site and it must be re-found.
    calibrated : bool, default False
        Set by :func:`calibrate_chamber_from_fiducials`, or by the user after
        running the protocol by hand.
    calibration_residual_deg : float or None
        RMS misfit of the fiducial calibration, when one was fitted.

    Examples
    --------
    >>> ChamberGeometry().ion_azimuth_deg(30.0)
    30.0
    >>> ChamberGeometry(rotation_sense=-1, rotation_offset_deg=90.0).ion_azimuth_deg(30.0)
    60.0
    """

    column_angle_deg: float = 54.0
    rotation_sense: int = 1
    rotation_offset_deg: float = 0.0
    stage_rotation_deg: float = 0.0
    tilt_axis_azimuth_deg: float = 0.0
    compucentric: bool = True
    calibrated: bool = False
    calibration_residual_deg: float | None = None

    def __post_init__(self) -> None:
        if not 0.0 < float(self.column_angle_deg) < 90.0:
            raise ValueError("ChamberGeometry.column_angle_deg must lie strictly between 0 and 90.")
        object.__setattr__(self, "rotation_sense", _sign(self.rotation_sense, "rotation_sense"))
        for name in ("rotation_offset_deg", "stage_rotation_deg", "tilt_axis_azimuth_deg"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"ChamberGeometry.{name} must be finite.")
        if self.calibration_residual_deg is not None and self.calibration_residual_deg < 0.0:
            raise ValueError("ChamberGeometry.calibration_residual_deg must be non-negative.")

    def ion_azimuth_deg(self, theta_sample_deg: float) -> float:
        """The pattern rotation for a sample azimuth, wrapped to ``[0, 360)``."""

        value = (
            self.rotation_sense * (float(theta_sample_deg) + self.stage_rotation_deg)
            + self.rotation_offset_deg
        )
        return round(float(wrap_azimuth_deg(value)), 12)

    @property
    def sem_foreshortening(self) -> float:
        """``cos T``: the SEM view's compression perpendicular to the tilt axis at tilt."""

        return math.cos(math.radians(self.column_angle_deg))

    def apparent_sem_length(self, length_um: float, azimuth_deg: float) -> float:
        """Apparent length, in the SEM view at tilt, of a segment at ``azimuth_deg``.

        The component along the tilt axis is seen at full length, the
        perpendicular component foreshortened by ``cos T``.
        """

        relative = math.radians(azimuth_deg - self.tilt_axis_azimuth_deg)
        along = length_um * math.cos(relative)
        across = length_um * math.sin(relative) * self.sem_foreshortening
        return float(math.hypot(along, across))

    @property
    def is_standard_column_angle(self) -> bool:
        """Whether ``T`` is one of the column angles the workbench offers (52, 54)."""

        return any(
            math.isclose(self.column_angle_deg, value) for value in SUPPORTED_COLUMN_ANGLES_DEG
        )

    def calibration_caveat(self) -> str:
        """The sentence every work order carries about ``theta_ion``."""

        if self.calibrated:
            residual = (
                f" (fiducial RMS residual {self.calibration_residual_deg:.2f} deg)"
                if self.calibration_residual_deg is not None
                else ""
            )
            return (
                f"theta_ion uses the calibrated s_R = {self.rotation_sense:+d} and "
                f"R0 = {self.rotation_offset_deg:.2f} deg{residual}; re-check the calibration "
                "after any service of the stage or a software update."
            )
        return (
            "UNCALIBRATED: theta_ion uses the placeholder s_R = +1, R0 = 0 deg. Its sign and "
            "offset are instrument-specific and have not been measured; run the fiducial "
            "calibration protocol before milling, or the lamella may be cut at the wrong azimuth."
        )

    def describe(self) -> str:
        """The chamber settings and the calibration status in prose."""

        tilt = (
            f"Chamber: ion column at T = {self.column_angle_deg:g} deg from the electron column; "
            f"at that stage tilt the ion view is the undistorted plan view and the SEM view is "
            f"foreshortened by cos T = {self.sem_foreshortening:.3f} across the tilt axis "
            f"(at azimuth {self.tilt_axis_azimuth_deg:g} deg). Stage rotation "
            f"R_stage = {self.stage_rotation_deg:g} deg"
            + ("" if self.compucentric else ", not compucentric (re-find the site after rotating)")
            + ". theta_ion = s_R (theta_S + R_stage) + R0. "
        )
        if not self.is_standard_column_angle:
            tilt += (
                f"Note: T = {self.column_angle_deg:g} deg is not one of the standard 52/54 deg "
                "column angles; confirm it against the instrument. "
            )
        return tilt + self.calibration_caveat()

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "column_angle_deg": float(self.column_angle_deg),
            "rotation_sense": int(self.rotation_sense),
            "rotation_offset_deg": float(self.rotation_offset_deg),
            "stage_rotation_deg": float(self.stage_rotation_deg),
            "tilt_axis_azimuth_deg": float(self.tilt_axis_azimuth_deg),
            "compucentric": bool(self.compucentric),
            "calibrated": bool(self.calibrated),
            "calibration_residual_deg": (
                None
                if self.calibration_residual_deg is None
                else float(self.calibration_residual_deg)
            ),
            "sem_foreshortening": self.sem_foreshortening,
            "calibration_caveat": self.calibration_caveat(),
        }


@dataclass(frozen=True, slots=True)
class FiducialObservation:
    """One fiducial trench of the azimuth calibration.

    Attributes
    ----------
    pattern_rotation_deg : float
        The rotation typed into the FIB pattern box when the fiducial was milled.
    sample_azimuth_deg : float
        The azimuth ``theta_S`` of the trench's *normal*, measured afterwards in
        the sample frame (for example on the registered SEM image). Only its
        value modulo 180 degrees is meaningful: a trench has no front.
    stage_rotation_deg : float, default 0
        Stage rotation when the fiducial was milled.
    """

    pattern_rotation_deg: float
    sample_azimuth_deg: float
    stage_rotation_deg: float = 0.0


def _wrap_signed(angle: np.ndarray, period: float) -> np.ndarray:
    return (angle + period / 2.0) % period - period / 2.0


def calibrate_chamber_from_fiducials(
    observations: Sequence[FiducialObservation],
    *,
    base: ChamberGeometry | None = None,
) -> ChamberGeometry:
    """Fit ``s_R`` and ``R0`` from fiducial trenches milled at known pattern rotations.

    Purpose
    -------
    The one-time calibration section 5.1 of the foundation document requires.
    Mill two or more short trenches at deliberately different pattern
    rotations, measure the azimuth of each trench normal in the sample frame,
    and this returns a :class:`ChamberGeometry` whose ``theta_ion`` reproduces
    them, with its misfit.

    Algorithm
    ---------
    For each candidate sense ``s`` in ``{+1, -1}`` the offset is the circular
    mean, with period 180 degrees, of ``rho_i - s (theta_i + R_i)``; the sense
    with the smaller RMS residual wins. Period 180 because a trench normal and
    its reverse are the same trench. Two observations determine the sense only
    if their rotations differ by something other than a multiple of 90 degrees;
    use rotations 0, 30 and 60 degrees.

    Parameters
    ----------
    observations : sequence of FiducialObservation
        At least two, at distinct pattern rotations.
    base : ChamberGeometry, optional
        Settings to carry over (column angle, stage rotation of the lamella
        itself). Defaults to :class:`ChamberGeometry`.

    Returns
    -------
    ChamberGeometry
        ``calibrated=True`` and ``calibration_residual_deg`` set.

    Raises
    ------
    ValueError
        Fewer than two observations, or rotations that cannot tell the two
        senses apart.

    Examples
    --------
    An instrument with a reversed sense and a 90 degree offset:

    >>> obs = [FiducialObservation(r, (90.0 - r) % 180.0) for r in (0.0, 30.0, 60.0)]
    >>> chamber = calibrate_chamber_from_fiducials(obs)
    >>> chamber.rotation_sense, round(chamber.rotation_offset_deg % 180.0, 6)
    (-1, 90.0)
    """

    if len(observations) < 2:
        raise ValueError("The azimuth calibration needs at least two fiducial trenches.")
    rho = np.array([item.pattern_rotation_deg for item in observations], dtype=np.float64)
    theta = np.array(
        [item.sample_azimuth_deg + item.stage_rotation_deg for item in observations],
        dtype=np.float64,
    )
    if np.ptp(_wrap_signed(rho - rho[0], 180.0)) < 1e-6:
        raise ValueError("All fiducials share one pattern rotation; vary it between trenches.")
    fits: list[tuple[float, int, float]] = []
    for sense in (1, -1):
        offsets = rho - sense * theta
        doubled = np.radians(2.0 * offsets)
        mean = math.degrees(math.atan2(np.sin(doubled).mean(), np.cos(doubled).mean())) / 2.0
        residual = _wrap_signed(offsets - mean, 180.0)
        rms = float(np.sqrt(np.mean(residual**2)))
        fits.append((rms, sense, mean))
    fits.sort()
    if math.isclose(fits[0][0], fits[1][0], abs_tol=1e-6):
        raise ValueError(
            "The fiducial rotations cannot distinguish the two rotation senses (they differ by "
            "multiples of 90 degrees); mill a trench at 30 or 60 degrees as well."
        )
    rms, sense, offset = fits[0]
    template = base or ChamberGeometry()
    return replace(
        template,
        rotation_sense=sense,
        rotation_offset_deg=float(wrap_azimuth_deg(offset)),
        calibrated=True,
        calibration_residual_deg=rms,
    )


def lamella_frame_graph(
    surface: SurfaceGeometry,
    *,
    sample_to_lamella: np.ndarray | None = None,
    lamella_to_holder: np.ndarray | None = None,
) -> FrameGraph:
    """The rigid chain ``E -> S -> L -> H`` as a :class:`~pytex.core.frames.FrameGraph`.

    Parameters
    ----------
    surface : SurfaceGeometry
        Supplies ``E -> S``.
    sample_to_lamella : (3, 3) array, optional
        Rows ``n_L, t_L, Z_s`` in sample components, from a solved plan.
    lamella_to_holder : (3, 3) array, optional
        One mounting hypothesis, from :meth:`pytex.fib.planning.MountModel.lamella_to_holder`.

    Returns
    -------
    FrameGraph
        Ask it for ``transform_between("specimen", "holder")`` to get the whole
        chain composed, which is how the planning code forms the
        crystal-to-holder orientation it hands to the TEM solver.
    """

    from pytex.tem.reconstruction import HOLDER_FRAME

    graph = FrameGraph(name="fib_lamella_chain")
    graph.add_transform(surface.frame_transform())
    if sample_to_lamella is not None:
        graph.add_transform(
            FrameTransform(
                source=FIB_SAMPLE_FRAME, target=LAMELLA_FRAME, rotation_matrix=sample_to_lamella
            )
        )
        if lamella_to_holder is not None:
            graph.add_transform(
                FrameTransform(
                    source=LAMELLA_FRAME, target=HOLDER_FRAME, rotation_matrix=lamella_to_holder
                )
            )
    elif lamella_to_holder is not None:
        raise ValueError("lamella_to_holder needs sample_to_lamella to connect to the chain.")
    return graph
