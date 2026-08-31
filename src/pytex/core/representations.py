"""Every numerical form of one rotation, and the conversions between them.

A rotation is a single geometric object; the numbers that denote it are not. A
grain measured by EBSD arrives as Bunge Euler angles, the same grain in a TEM
report arrives as a zone axis and an in-plane direction, a misorientation is
quoted as an axis-angle pair, a rolling texture component is *named* by
``(hkl)[uvw]``, and an orientation-distribution sampler wants a coordinate in
which uniform sampling really is uniform. All of those describe rotations, and a
user should not have to know six conversion recipes to move between them.

This module is the single conversion hub. It adds the two representations the
rest of the library had no constructor for — **homochoric** and **cubochoric** —
provides fully vectorized conversions between every pair, and gives one call,
:func:`orientation_representations`, that takes an orientation in any form and
returns *all* of them together in standard notation.

The representations
-------------------

================ ================= ==========================================
name             shape             what it is good for
================ ================= ==========================================
matrix           ``(3, 3)``        applying the rotation to vectors
quaternion       ``(4,)``          composing rotations; numerically stable
axis-angle       ``(3,)`` + scalar reading the physics off directly
Rodrigues        ``(3,)``          fundamental zones are convex polyhedra
Rodrigues-Frank  ``(4,)``          the same, still invertible at ``omega = pi``
Euler (Bunge)    ``(3,)``          the texture community's lingua franca
Euler (ZYZ)      ``(3,)``          Matthies/Roe school, and ``abg`` imports
homochoric       ``(3,)``          equal-volume; a ball of radius ``(3 pi/4)^(1/3)``
cubochoric       ``(3,)``          equal-volume **and** a cube: uniform grids
``(hkl)[uvw]``   two index triples  naming a texture component
================ ================= ==========================================

The mathematics
---------------

Write the rotation as an axis-angle pair :math:`(\\hat{n}, \\omega)` with
:math:`\\omega \\in [0, \\pi]`. Then

.. math::

   q &= \\left(\\cos\\tfrac{\\omega}{2},\\; \\hat{n}\\sin\\tfrac{\\omega}{2}\\right), \\\\
   \\boldsymbol{\\rho} &= \\hat{n}\\tan\\tfrac{\\omega}{2}
       \\quad\\text{(Rodrigues)}, \\\\
   \\mathbf{h} &= \\hat{n}
       \\left[\\tfrac{3}{4}\\left(\\omega - \\sin\\omega\\right)\\right]^{1/3}
       \\quad\\text{(homochoric)}.

The homochoric vector is the one whose *volume element is the volume element of
the rotation group*. The invariant measure on SO(3) in axis-angle coordinates
carries the factor :math:`(1 - \\cos\\omega)`, so a uniform cloud of points in
axis-angle space is badly non-uniform in orientation. Choosing the radial
function :math:`f(\\omega)` such that :math:`f^2 \\mathrm{d}f \\propto
(1-\\cos\\omega)\\,\\mathrm{d}\\omega` gives exactly the cube-root above, and the
whole of SO(3) becomes a ball of radius

.. math::

   R_{1} = \\left(\\tfrac{3\\pi}{4}\\right)^{1/3} \\approx 1.3306 .

Cubochoric coordinates then map that ball to a **cube** of edge
:math:`a_{p} = \\pi^{2/3}` — same volume, :math:`\\pi^{2}` — by the equal-volume
map of Rosca, Morawiec and De Graef. A uniform Cartesian grid in the cube is
therefore a uniform grid in orientation space, which is what makes cubochoric
sampling the standard choice for dictionary indexing.

The equal-volume cube-to-ball map
---------------------------------

The map is built from two facts, and both are worth stating because they explain
every constant in the code.

1. **Nested surfaces.** The cube of half-edge :math:`z` (in the rescaled
   coordinates) must map onto the sphere of radius :math:`r(z)` enclosing the
   same volume: :math:`(2z)^{3} = \\tfrac{4}{3}\\pi r^{3}`, hence
   :math:`r = z\\,(6/\\pi)^{1/3}`.
2. **Each cube face maps to its own spherical sector.** The six faces must tile
   the sphere, and the map commutes with the octahedral symmetry, so the
   boundary between the ``+z`` and ``+x`` images lies on the mirror plane
   :math:`x = z` — that is, each face maps to the curvilinear square
   :math:`\\{n_{z} \\ge |n_{x}|, |n_{y}|\\}`, never to a spherical cap.

Those two conditions fix the map. Within the ``+z`` pyramid it factors through a
planar wedge and a Lambert azimuthal equal-area lift; see
:func:`homochoric_from_cubochoric` for the closed form and
``docs/site/theory/orientation_representations.md`` for the derivation.

Vectorization
-------------

Every free function here is vectorized over a leading batch axis and takes
``(n, ...)`` arrays, including the Euler conversions, which the object-level
`pytex.core.batches.RotationSet` methods still perform in a Python loop. Use
:func:`convert_orientations` for bulk conversion of a large orientation cloud
and :class:`OrientationRepresentationSet` when the batch also has to be
*reported*.

See Also
--------
`pytex.core.orientation` : the `Rotation` and `Orientation` objects themselves.
`pytex.core.notation` : the bracket, overbar, and reciprocal-star rules used here.
`docs/standards/notation_and_conventions.md` : the governing notation standard.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._angles import angle_between_unit_vectors_rad
from pytex.core._arrays import as_float_array, normalize_quaternions
from pytex.core.hexagonal import is_hexagonal_phase
from pytex.core.miller import (
    direction_uvw_to_uvtw_array,
    plane_hkl_to_hkil_array,
    reduce_indices,
)
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.orientation import (
    Orientation,
    Rotation,
    _euler_axes_for_convention,
    _matrices_to_repeated_axis_euler,
    matrices_to_quaternions,
    quaternions_from_axes_angles,
    quaternions_from_rodrigues,
    quaternions_multiply,
    quaternions_to_axes_angles,
    quaternions_to_matrices,
    quaternions_to_rodrigues,
    specimen_direction_vector,
)
from pytex.core.provenance import ProvenanceRecord

__all__ = [
    "CUBOCHORIC_CUBE_EDGE",
    "CUBOCHORIC_CUBE_HALF_EDGE",
    "HOMOCHORIC_BALL_RADIUS",
    "IDEAL_ORIENTATION_SCHEMA",
    "ORIENTATION_REPRESENTATIONS_SCHEMA",
    "IdealOrientationIndices",
    "OrientationRepresentationSet",
    "OrientationRepresentations",
    "RepresentationKind",
    "canonical_quaternions",
    "convert_orientations",
    "cubochoric_from_homochoric",
    "cubochoric_from_quaternions",
    "homochoric_from_cubochoric",
    "homochoric_from_quaternions",
    "ideal_orientation_indices",
    "orientation_representations",
    "quaternions_from_cubochoric",
    "quaternions_from_euler_angles",
    "quaternions_from_homochoric",
    "quaternions_to_euler_angles",
    "rotation_representations",
]

#: Schema identifier of the all-representations payload.
ORIENTATION_REPRESENTATIONS_SCHEMA = "pytex.orientation_representations/1"

#: Schema identifier of the ideal-orientation ``(hkl)[uvw]`` payload.
IDEAL_ORIENTATION_SCHEMA = "pytex.ideal_orientation_indices/1"

#: Radius of the homochoric ball, ``(3 pi / 4) ** (1 / 3)``.
#:
#: Every rotation of SO(3) lies inside it, and the bounding sphere is the set of
#: rotations by ``pi`` — where antipodal points denote the *same* rotation, which
#: is why the ball is a model of SO(3) rather than of the unit quaternions.
HOMOCHORIC_BALL_RADIUS = float((3.0 * np.pi / 4.0) ** (1.0 / 3.0))

#: Edge length of the cubochoric cube, ``pi ** (2 / 3)``.
#:
#: Chosen so the cube volume ``pi ** 2`` equals the homochoric ball volume; the
#: map between them is then volume-preserving point by point.
CUBOCHORIC_CUBE_EDGE = float(np.pi ** (2.0 / 3.0))

#: Half the cubochoric cube edge: the coordinate bound, ``|c_i| <= this``.
CUBOCHORIC_CUBE_HALF_EDGE = CUBOCHORIC_CUBE_EDGE / 2.0

# Scale taking a cubochoric coordinate into the working coordinates in which the
# nested-surface law reads r = sqrt(6 / pi) * z rather than (6 / pi) ** (1 / 3) * z.
_CUBE_SCALE = float((np.pi / 6.0) ** (1.0 / 6.0))

# Radial prefactor of the square-to-wedge map, fixed by requiring that the face
# edge land on the sector boundary: 2^(1/4) * sqrt(6 / pi).
_WEDGE_PREFACTOR = float(2.0 ** 0.25 * np.sqrt(6.0 / np.pi))

_SQRT2 = float(np.sqrt(2.0))

# Below this the rotation is treated as the identity, where the axis is
# undefined and every representation degenerates to zero.
_ANGLE_EPS = 1e-12


class RepresentationKind(StrEnum):
    """The numerical forms a rotation can be written in.

    Used by :func:`convert_orientations` to name the source and target of a bulk
    conversion, so a call cannot silently misread Rodrigues vectors as
    homochoric ones — they have the same shape and differ only in meaning.
    """

    #: ``(n, 3, 3)`` proper-orthogonal matrices, active convention ``v' = R v``.
    MATRIX = "matrix"
    #: ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order.
    QUATERNION = "quaternion"
    #: ``(n, 4)`` rows of ``(nx, ny, nz, omega_rad)``.
    AXIS_ANGLE = "axis_angle"
    #: ``(n, 3)`` Rodrigues vectors ``n tan(omega / 2)``.
    RODRIGUES = "rodrigues"
    #: ``(n, 4)`` Rodrigues-Frank rows ``(nx, ny, nz, tan(omega / 2))``; the
    #: magnitude is projective, so ``omega = pi`` is representable rather than
    #: overflowing as it does in the 3-vector form.
    RODRIGUES_FRANK = "rodrigues_frank"
    #: ``(n, 3)`` Bunge ZXZ angles ``(phi1, Phi, phi2)``.
    EULER_BUNGE = "euler_bunge"
    #: ``(n, 3)`` ZYZ angles ``(alpha, beta, gamma)``, Matthies/Roe school.
    EULER_MATTHIES = "euler_matthies"
    #: ``(n, 3)`` homochoric vectors, equal-volume, inside a ball.
    HOMOCHORIC = "homochoric"
    #: ``(n, 3)`` cubochoric vectors, equal-volume, inside a cube.
    CUBOCHORIC = "cubochoric"


# --------------------------------------------------------------------------- #
# Homochoric coordinates
# --------------------------------------------------------------------------- #


def homochoric_from_quaternions(quaternions: ArrayLike) -> np.ndarray:
    """Homochoric vectors of a batch of quaternions.

    What it does
        Maps each rotation to the point :math:`\\hat{n}\\,[\\tfrac{3}{4}(\\omega -
        \\sin\\omega)]^{1/3}` of the homochoric ball, the coordinate in which the
        volume element of SO(3) is the ordinary Euclidean one.

    When to use it
        Whenever a *density* in orientation space is at stake: kernel-density
        ODF estimation, nearest-neighbour searches that must not be biased by
        rotation angle, and any Monte-Carlo sampling of orientations. Do not use
        it to compose rotations — it is a chart, not an algebra.

    Parameters
    ----------
    quaternions:
        ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order. Normalized on
        entry, so an approximately-unit input is accepted.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` homochoric vectors, each of norm at most
        :data:`HOMOCHORIC_BALL_RADIUS`.

    Notes
    -----
    The identity maps to the origin, and the map is continuous there even though
    the rotation axis is undefined, because the radius vanishes as
    :math:`\\omega \\to 0` like :math:`(\\omega^{3}/8)^{1/3} = \\omega/2`.

    See Also
    --------
    quaternions_from_homochoric : the inverse.
    cubochoric_from_homochoric : on to the equal-volume cube.
    """

    axes, angles = quaternions_to_axes_angles(quaternions)
    radii = np.cbrt(0.75 * (angles - np.sin(angles)))
    return np.ascontiguousarray(axes * radii[:, None])


def quaternions_from_homochoric(homochoric: ArrayLike) -> np.ndarray:
    """Quaternions of a batch of homochoric vectors; inverse of
    :func:`homochoric_from_quaternions`.

    What it does
        Recovers the rotation angle by solving :math:`\\omega - \\sin\\omega =
        \\tfrac{4}{3}\\|\\mathbf{h}\\|^{3}` on :math:`[0, \\pi]`, then rebuilds the
        quaternion from that angle and the direction of ``h``.

    When to use it
        To return to the algebraic representation after any work done in the
        equal-volume chart — after sampling, after interpolating a density, or
        after reading a cubochoric dictionary grid.

    Parameters
    ----------
    homochoric:
        ``(n, 3)`` homochoric vectors. A norm exceeding
        :data:`HOMOCHORIC_BALL_RADIUS` by more than a rounding tolerance is an
        error rather than being clipped silently: it means the caller mixed up
        Rodrigues and homochoric vectors, which this shape cannot otherwise
        distinguish.

    Returns
    -------
    np.ndarray
        ``(n, 4)`` unit quaternions with non-negative scalar part.

    Raises
    ------
    ValueError
        If any input lies outside the homochoric ball.

    Notes
    -----
    The angle equation has no closed-form solution, so it is solved by a
    vectorized bisection on :math:`[0, \\pi]`, where the left-hand side is
    strictly increasing. Sixty iterations bring the bracket below machine
    precision, and the cost is a fixed sixty array operations regardless of the
    batch size — unlike a per-element root finder, which would be a Python loop.
    """

    vectors = as_float_array(homochoric, shape=(None, 3))
    radii = np.linalg.norm(vectors, axis=1)
    if np.any(radii > HOMOCHORIC_BALL_RADIUS * (1.0 + 1e-9)):
        worst = float(np.max(radii))
        raise ValueError(
            f"A homochoric vector of norm {worst:.6f} lies outside the homochoric ball "
            f"of radius {HOMOCHORIC_BALL_RADIUS:.6f}. Rodrigues vectors have the same "
            "shape and are unbounded; check which representation was passed."
        )
    angles = _homochoric_angles(np.clip(radii, 0.0, HOMOCHORIC_BALL_RADIUS))
    finite = radii > _ANGLE_EPS
    axes = np.zeros_like(vectors)
    axes[:, 2] = 1.0
    axes[finite] = vectors[finite] / radii[finite, None]
    return quaternions_from_axes_angles(axes, angles)


def _homochoric_angles(radii: np.ndarray) -> np.ndarray:
    """Solve ``omega - sin(omega) = (4 / 3) r ** 3`` on ``[0, pi]``, vectorized."""

    targets = (4.0 / 3.0) * np.power(radii, 3.0)
    lower = np.zeros_like(targets)
    upper = np.full_like(targets, np.pi)
    for _ in range(60):
        middle = 0.5 * (lower + upper)
        too_small = (middle - np.sin(middle)) < targets
        lower = np.where(too_small, middle, lower)
        upper = np.where(too_small, upper, middle)
    return np.asarray(0.5 * (lower + upper), dtype=np.float64)


# --------------------------------------------------------------------------- #
# Cubochoric coordinates
# --------------------------------------------------------------------------- #


def _sorted_pyramid(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Permute each row so its largest-magnitude component is last, and made positive.

    Returns the permuted rows, the permutation applied (as index arrays), and
    the sign that was stripped from the now-last component. The cube-to-ball map
    is derived on the ``+z`` pyramid and extended to the other five by the
    octahedral equivariance it was built with — a cyclic permutation of the axes
    is a rotation of the cube and of the ball alike, and the mirror ``z -> -z``
    negates the third component of the image and nothing else — so this is the
    whole of the case analysis.
    """

    last = np.argmax(np.abs(vectors), axis=1)
    order = np.empty((vectors.shape[0], 3), dtype=np.intp)
    for column in range(3):
        order[:, column] = (last + 1 + column) % 3
    rows = np.arange(vectors.shape[0])[:, None]
    permuted = np.array(vectors[rows, order], dtype=np.float64)
    signs = np.where(permuted[:, 2] < 0.0, -1.0, 1.0)
    permuted[:, 2] = np.abs(permuted[:, 2])
    return permuted, order, signs


def _unsort_pyramid(values: np.ndarray, order: np.ndarray) -> np.ndarray:
    """Undo the axis permutation of :func:`_sorted_pyramid`."""

    restored = np.empty_like(values)
    rows = np.arange(values.shape[0])[:, None]
    restored[rows, order] = values
    return restored


def homochoric_from_cubochoric(cubochoric: ArrayLike) -> np.ndarray:
    """Homochoric vectors of cubochoric coordinates: the equal-volume cube-to-ball map.

    What it does
        Applies the Rosca-Morawiec-De Graef map carrying the cube of edge
        :data:`CUBOCHORIC_CUBE_EDGE` onto the ball of radius
        :data:`HOMOCHORIC_BALL_RADIUS`, preserving volume at every point.

    When to use it
        When a uniform grid or a uniform random sample of orientations is
        required. Sampling the cube uniformly and mapping through here is the
        only elementary construction that gives a genuinely uniform sample of
        SO(3); sampling Euler angles uniformly does not, and sampling axis-angle
        uniformly does not either.

    Parameters
    ----------
    cubochoric:
        ``(n, 3)`` cubochoric coordinates, each component within
        :data:`CUBOCHORIC_CUBE_HALF_EDGE`.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` homochoric vectors.

    Raises
    ------
    ValueError
        If any coordinate lies outside the cube.

    Notes
    -----
    **Algorithm.** For a point of the ``+z`` pyramid, after rescaling by
    :math:`(\\pi/6)^{1/6}` so that the nested-surface law reads
    :math:`r = \\sqrt{6/\\pi}\\,z`:

    1. The face coordinates are carried into a planar wedge of half-angle
       :math:`\\pi/4` by, for :math:`|y| \\le |x|`,

       .. math::

          \\alpha = \\frac{\\pi y}{12 x}, \\qquad
          k = \\frac{2^{1/4}\\sqrt{6/\\pi}\\; x}{\\sqrt{\\sqrt{2} - \\cos\\alpha}},

       and :math:`(T_{1}, T_{2}) = k\\,(\\sqrt{2}\\cos\\alpha - 1,\\;
       \\sqrt{2}\\sin\\alpha)`; the roles of ``x`` and ``y`` swap on the other
       side of the face diagonal.
    2. The wedge is lifted onto the sphere of radius :math:`R = \\sqrt{6/\\pi}\\,z`
       by the Lambert azimuthal equal-area map,
       :math:`(T_{1}f, T_{2}f, R - \\rho^{2}/2R)` with
       :math:`f = \\sqrt{1 - \\rho^{2}/4R^{2}}` and :math:`\\rho^{2} = T_{1}^{2} +
       T_{2}^{2}`.

    Both steps preserve area, and step 1's prefactor is fixed — not fitted — by
    requiring the face edge to land on the sector boundary
    :math:`\\theta = \\arctan(1/\\cos\\varphi)`.

    See Also
    --------
    cubochoric_from_homochoric : the inverse.
    """

    coordinates = as_float_array(cubochoric, shape=(None, 3))
    if coordinates.size and np.any(
        np.abs(coordinates) > CUBOCHORIC_CUBE_HALF_EDGE * (1.0 + 1e-9)
    ):
        worst = float(np.max(np.abs(coordinates)))
        raise ValueError(
            f"A cubochoric coordinate of {worst:.6f} lies outside the cubochoric cube, "
            f"whose components are bounded by {CUBOCHORIC_CUBE_HALF_EDGE:.6f}."
        )
    if coordinates.size == 0:
        return np.zeros((0, 3), dtype=np.float64)

    scaled = np.clip(
        coordinates, -CUBOCHORIC_CUBE_HALF_EDGE, CUBOCHORIC_CUBE_HALF_EDGE
    ) * _CUBE_SCALE
    permuted, order, signs = _sorted_pyramid(scaled)
    x_component = permuted[:, 0]
    y_component = permuted[:, 1]
    z_component = permuted[:, 2]

    # Face diagonal: which of x, y plays the role of the "radial" coordinate.
    x_dominant = np.abs(y_component) <= np.abs(x_component)
    major = np.where(x_dominant, x_component, y_component)
    minor = np.where(x_dominant, y_component, x_component)
    safe_major = np.where(np.abs(major) < _ANGLE_EPS, 1.0, major)
    alpha = (np.pi / 12.0) * (minor / safe_major)
    scale = _WEDGE_PREFACTOR * major / np.sqrt(_SQRT2 - np.cos(alpha))
    principal = (_SQRT2 * np.cos(alpha) - 1.0) * scale
    secondary = _SQRT2 * np.sin(alpha) * scale
    first = np.where(x_dominant, principal, secondary)
    second = np.where(x_dominant, secondary, principal)

    sphere_radius = np.sqrt(6.0 / np.pi) * z_component
    safe_radius = np.where(sphere_radius < _ANGLE_EPS, 1.0, sphere_radius)
    planar_sq = first * first + second * second
    lift = np.sqrt(np.clip(1.0 - planar_sq / (4.0 * safe_radius * safe_radius), 0.0, None))
    ball = np.column_stack(
        [
            first * lift,
            second * lift,
            sphere_radius - planar_sq / (2.0 * safe_radius),
        ]
    )
    ball[sphere_radius < _ANGLE_EPS] = 0.0
    ball[:, 2] *= signs
    return np.ascontiguousarray(_unsort_pyramid(ball, order))


def cubochoric_from_homochoric(homochoric: ArrayLike) -> np.ndarray:
    """Cubochoric coordinates of homochoric vectors; inverse of
    :func:`homochoric_from_cubochoric`.

    What it does
        Inverts the equal-volume map analytically, step for step: undo the
        Lambert lift, then undo the square-to-wedge map, whose angular part
        inverts in closed form through
        :math:`\\alpha = \\varphi - \\arcsin(\\sin\\varphi/\\sqrt{2})`.

    When to use it
        To place measured orientations on a cubochoric grid — the indexing step
        of dictionary-based EBSD or TEM pattern matching — and to check that a
        sampling scheme really did fill orientation space.

    Parameters
    ----------
    homochoric:
        ``(n, 3)`` homochoric vectors inside the homochoric ball.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` cubochoric coordinates.

    Raises
    ------
    ValueError
        If any input lies outside the homochoric ball.
    """

    vectors = as_float_array(homochoric, shape=(None, 3))
    if vectors.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    radii = np.linalg.norm(vectors, axis=1)
    if np.any(radii > HOMOCHORIC_BALL_RADIUS * (1.0 + 1e-9)):
        worst = float(np.max(radii))
        raise ValueError(
            f"A homochoric vector of norm {worst:.6f} lies outside the homochoric ball "
            f"of radius {HOMOCHORIC_BALL_RADIUS:.6f}."
        )

    permuted, order, signs = _sorted_pyramid(vectors)
    sphere_radius = np.linalg.norm(permuted, axis=1)
    safe_radius = np.where(sphere_radius < _ANGLE_EPS, 1.0, sphere_radius)
    z_component = sphere_radius * np.sqrt(np.pi / 6.0)

    planar_sq = np.clip(2.0 * safe_radius * (safe_radius - permuted[:, 2]), 0.0, None)
    planar = np.sqrt(planar_sq)
    in_plane = np.linalg.norm(permuted[:, :2], axis=1)
    safe_in_plane = np.where(in_plane < _ANGLE_EPS, 1.0, in_plane)
    first = permuted[:, 0] * planar / safe_in_plane
    second = permuted[:, 1] * planar / safe_in_plane

    x_dominant = np.abs(second) <= np.abs(first)
    principal = np.where(x_dominant, first, second)
    secondary = np.where(x_dominant, second, first)
    # Fold the wedge onto its positive-principal half before inverting the
    # angular relation, which is only valid for an azimuth within +/- pi/4.
    major_sign = np.where(principal < 0.0, -1.0, 1.0)
    azimuth = np.arctan2(secondary * major_sign, principal * major_sign)
    alpha = azimuth - np.arcsin(np.clip(np.sin(azimuth) / _SQRT2, -1.0, 1.0))
    wedge_radius = np.sqrt(
        np.clip((3.0 - 2.0 * _SQRT2 * np.cos(alpha)) / (_SQRT2 - np.cos(alpha)), 0.0, None)
    )
    safe_wedge = np.where(wedge_radius < _ANGLE_EPS, 1.0, wedge_radius)
    major = major_sign * np.hypot(principal, secondary) / (_WEDGE_PREFACTOR * safe_wedge)
    minor = major * (12.0 / np.pi) * alpha

    x_component = np.where(x_dominant, major, minor)
    y_component = np.where(x_dominant, minor, major)
    at_face_centre = in_plane < _ANGLE_EPS
    x_component = np.where(at_face_centre, 0.0, x_component)
    y_component = np.where(at_face_centre, 0.0, y_component)

    cube = np.column_stack([x_component, y_component, z_component])
    cube[sphere_radius < _ANGLE_EPS] = 0.0
    cube[:, 2] *= signs
    restored = _unsort_pyramid(cube, order) / _CUBE_SCALE
    return np.ascontiguousarray(
        np.clip(restored, -CUBOCHORIC_CUBE_HALF_EDGE, CUBOCHORIC_CUBE_HALF_EDGE)
    )


def cubochoric_from_quaternions(quaternions: ArrayLike) -> np.ndarray:
    """Cubochoric coordinates of a batch of quaternions.

    The composition of :func:`homochoric_from_quaternions` and
    :func:`cubochoric_from_homochoric`, provided because it is the pairing that
    actually gets used: quaternions are what the library stores, and cubochoric
    is what uniform sampling and dictionary indexing need.
    """

    return cubochoric_from_homochoric(homochoric_from_quaternions(quaternions))


def quaternions_from_cubochoric(cubochoric: ArrayLike) -> np.ndarray:
    """Quaternions of a batch of cubochoric coordinates; inverse of
    :func:`cubochoric_from_quaternions`."""

    return quaternions_from_homochoric(homochoric_from_cubochoric(cubochoric))


# --------------------------------------------------------------------------- #
# Vectorized Euler conversions
# --------------------------------------------------------------------------- #


def _axis_quaternions(axis: str, angles: np.ndarray) -> np.ndarray:
    """``(n, 4)`` quaternions of rotations by ``angles`` about a named axis."""

    half = 0.5 * angles
    quaternions = np.zeros((angles.shape[0], 4), dtype=np.float64)
    quaternions[:, 0] = np.cos(half)
    quaternions[:, {"x": 1, "y": 2, "z": 3}[axis]] = np.sin(half)
    return quaternions


def quaternions_from_euler_angles(
    angles: ArrayLike,
    *,
    convention: str = "bunge",
    degrees: bool = True,
) -> np.ndarray:
    """Quaternions from a batch of Euler-angle triples, vectorized.

    What it does
        Composes the three elementary rotations of the named convention as
        quaternions, in one array operation over the whole batch.

    When to use it
        For bulk conversion of an EBSD scan or a simulated orientation cloud.
        The object-level path, ``EulerSet.to_rotation_set``, builds one
        `Rotation` per row in a Python loop; this does the same arithmetic with
        four array multiplications and matches it exactly.

    Parameters
    ----------
    angles:
        ``(n, 3)`` triples in the order the convention names them.
    convention:
        ``"bunge"`` (ZXZ), ``"matthies"`` or ``"abg"`` (both ZYZ). Aliases such
        as ``"zxz"`` and ``"zyz"`` are accepted.
    degrees:
        Interpret the angles as degrees (default) rather than radians.

    Returns
    -------
    np.ndarray
        ``(n, 4)`` unit quaternions.

    See Also
    --------
    quaternions_to_euler_angles : the inverse.
    `pytex.core.batches.EulerSet` : the typed carrier that keeps the convention
        attached to the numbers.
    """

    triples = as_float_array(angles, shape=(None, 3))
    if degrees:
        triples = np.deg2rad(triples)
    first_axis, second_axis, third_axis = _euler_axes_for_convention(convention)
    return quaternions_multiply(
        quaternions_multiply(
            _axis_quaternions(first_axis, triples[:, 0]),
            _axis_quaternions(second_axis, triples[:, 1]),
        ),
        _axis_quaternions(third_axis, triples[:, 2]),
    )


def quaternions_to_euler_angles(
    quaternions: ArrayLike,
    *,
    convention: str = "bunge",
    degrees: bool = True,
) -> np.ndarray:
    """Euler-angle triples of a batch of quaternions, vectorized.

    What it does
        Extracts the angle triple of the named convention, wrapping each angle
        into ``[0, 2 pi)`` and resolving the gimbal-degenerate cases (Bunge
        ``Phi = 0`` or ``pi``) by setting the third angle to zero — the same
        choice `Rotation.to_euler` makes, so the two agree row for row.

    Parameters
    ----------
    quaternions:
        ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order.
    convention:
        ``"bunge"``, ``"matthies"``, or ``"abg"``.
    degrees:
        Return degrees (default) rather than radians.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` angle triples.
    """

    matrices = quaternions_to_matrices(quaternions)
    radians = np.mod(
        _matrices_to_repeated_axis_euler(matrices, convention=convention), 2.0 * np.pi
    )
    return np.asarray(np.rad2deg(radians) if degrees else radians, dtype=np.float64)


# --------------------------------------------------------------------------- #
# The generic converter
# --------------------------------------------------------------------------- #


def _to_quaternions(values: ArrayLike, kind: RepresentationKind) -> np.ndarray:
    if kind is RepresentationKind.MATRIX:
        return matrices_to_quaternions(values)
    if kind is RepresentationKind.QUATERNION:
        return normalize_quaternions(values)
    if kind is RepresentationKind.AXIS_ANGLE:
        rows = as_float_array(values, shape=(None, 4))
        return quaternions_from_axes_angles(rows[:, :3], rows[:, 3])
    if kind is RepresentationKind.RODRIGUES:
        return quaternions_from_rodrigues(values, frank=False)
    if kind is RepresentationKind.RODRIGUES_FRANK:
        return quaternions_from_rodrigues(values, frank=True)
    if kind is RepresentationKind.EULER_BUNGE:
        return quaternions_from_euler_angles(values, convention="bunge", degrees=True)
    if kind is RepresentationKind.EULER_MATTHIES:
        return quaternions_from_euler_angles(values, convention="matthies", degrees=True)
    if kind is RepresentationKind.HOMOCHORIC:
        return quaternions_from_homochoric(values)
    return quaternions_from_cubochoric(values)


def _from_quaternions(quaternions: np.ndarray, kind: RepresentationKind) -> np.ndarray:
    if kind is RepresentationKind.MATRIX:
        return quaternions_to_matrices(quaternions)
    if kind is RepresentationKind.QUATERNION:
        return quaternions
    if kind is RepresentationKind.AXIS_ANGLE:
        axes, angles = quaternions_to_axes_angles(quaternions)
        return np.column_stack([axes, angles])
    if kind is RepresentationKind.RODRIGUES:
        return quaternions_to_rodrigues(quaternions, frank=False)
    if kind is RepresentationKind.RODRIGUES_FRANK:
        return quaternions_to_rodrigues(quaternions, frank=True)
    if kind is RepresentationKind.EULER_BUNGE:
        return quaternions_to_euler_angles(quaternions, convention="bunge", degrees=True)
    if kind is RepresentationKind.EULER_MATTHIES:
        return quaternions_to_euler_angles(
            quaternions, convention="matthies", degrees=True
        )
    if kind is RepresentationKind.HOMOCHORIC:
        return homochoric_from_quaternions(quaternions)
    return cubochoric_from_quaternions(quaternions)


def convert_orientations(
    values: ArrayLike,
    *,
    source: RepresentationKind | str,
    target: RepresentationKind | str,
) -> np.ndarray:
    """Convert a batch of rotations from one representation to another.

    What it does
        Routes every conversion through the quaternion form, which is the one
        representation with no singularities and no convention branches, so
        there are ten conversions to maintain rather than ninety.

    When to use it
        For bulk numerical conversion where the caller already knows what it
        has and what it wants — file import and export, feeding a sampler,
        preparing a dictionary grid. When the goal is to *report* an
        orientation rather than to compute with it, use
        :func:`orientation_representations`, which returns every form at once
        with the notation and the conventions attached.

    Parameters
    ----------
    values:
        The batch, shaped as the source kind requires; see
        :class:`RepresentationKind`. Angles are in degrees for the Euler kinds
        and in radians for the axis-angle kind, matching the rest of the
        library.
    source, target:
        :class:`RepresentationKind` members or their string values.

    Returns
    -------
    np.ndarray
        The batch in the target representation.

    Raises
    ------
    ValueError
        If a kind is not recognized, or the input violates the domain of the
        source representation.

    Examples
    --------
    Rolling-texture Bunge angles to the cubochoric coordinates a uniform
    dictionary grid is indexed on::

        cube = convert_orientations(
            [[35.0, 45.0, 0.0]], source="euler_bunge", target="cubochoric"
        )

    Notes
    -----
    Round-tripping is exact to floating-point tolerance for every pair, with one
    unavoidable exception: at the gimbal-degenerate Euler angles only the sum or
    difference of the outer angles is determined, so a triple with
    ``Phi = 0`` and ``phi2 != 0`` returns with the same rotation but different
    numbers.
    """

    source_kind = RepresentationKind(str(source))
    target_kind = RepresentationKind(str(target))
    quaternions = _to_quaternions(values, source_kind)
    return _from_quaternions(quaternions, target_kind)


# --------------------------------------------------------------------------- #
# Ideal-orientation (hkl)[uvw] recovery
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class IdealOrientationIndices:
    """The ``(hkl)[uvw]`` a texture component would be named by.

    Purpose
    -------
    The inverse of `pytex.core.orientation.Orientation.from_miller`. An
    orientation is a continuum of numbers; a *name* like ``(110)[1 -1 2]`` is a
    pair of small integers, and only a measure-zero set of orientations has an
    exact one. This object therefore reports the best integer pair **together
    with how far it is from the truth**, because an ideal-orientation label
    quoted without its deviation is a claim the data does not support.

    Attributes
    ----------
    hkl : np.ndarray
        The three plane indices, reduced by their common divisor, of the crystal
        plane closest to the specimen plane (ND by default).
    uvw : np.ndarray
        The three direction indices closest to the specimen direction (RD by
        default).
    plane_deviation_deg : float
        Angle between the true crystal-plane normal along the specimen normal
        and the reported ``(hkl)``.
    direction_deviation_deg : float
        The same for ``[uvw]``.
    hkil, uvtw : np.ndarray or None
        The four-index Miller-Bravais forms, present only for hexagonal and
        trigonal phases, where the international literature writes four indices.
    max_index : int
        The search bound; a larger bound can only reduce the deviations.
    specimen_plane_normal, specimen_direction : np.ndarray
        The specimen axes the indices were referred to, so the label cannot be
        read against the wrong geometry.
    phase_name : str
    """

    hkl: np.ndarray
    uvw: np.ndarray
    plane_deviation_deg: float
    direction_deviation_deg: float
    max_index: int
    specimen_plane_normal: np.ndarray
    specimen_direction: np.ndarray
    phase_name: str
    hkil: np.ndarray | None = None
    uvtw: np.ndarray | None = None

    def __post_init__(self) -> None:
        for name in ("hkl", "uvw"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=np.int64).reshape(3)
            )
        for name in ("specimen_plane_normal", "specimen_direction"):
            object.__setattr__(self, name, as_float_array(getattr(self, name), shape=(3,)))

    @property
    def label(self) -> str:
        """The component name in plain text, ``(hkl)[uvw]``."""

        return (
            f"{format_plane_indices(tuple(int(v) for v in self.hkl), style='plain')}"
            f"{format_direction_indices(tuple(int(v) for v in self.uvw), style='plain')}"
        )

    @property
    def mathtext_label(self) -> str:
        """The component name for a figure, with overbarred negative indices."""

        plane = format_plane_indices(
            tuple(int(v) for v in self.hkl), style="mathtext"
        ).strip("$")
        direction = format_direction_indices(
            tuple(int(v) for v in self.uvw), style="mathtext"
        ).strip("$")
        return f"${plane}{direction}$"

    @property
    def is_exact(self) -> bool:
        """Whether both deviations are below a thousandth of a degree.

        The threshold is deliberately tight: an exact ideal orientation is a
        statement about the crystallography, not about the measurement, and a
        component that misses by a hundredth of a degree is a different
        orientation that happens to be close.
        """

        return (
            self.plane_deviation_deg < 1e-3 and self.direction_deviation_deg < 1e-3
        )

    def describe(self) -> str:
        """Convention-explicit prose naming the component and its deviation."""

        four_index = ""
        if self.hkil is not None and self.uvtw is not None:
            four_index = (
                " In the four-index Miller-Bravais notation the international "
                "literature uses for this lattice, that is "
                f"{format_plane_indices(tuple(int(v) for v in self.hkil), style='plain')}"
                f"{format_direction_indices(tuple(int(v) for v in self.uvtw), style='plain')}."
            )
        verdict = (
            "The orientation is exactly this ideal component to within a "
            "thousandth of a degree."
            if self.is_exact
            else (
                "The component is the nearest label, not an identity: the plane misses "
                f"by {self.plane_deviation_deg:.3f} deg and the direction by "
                f"{self.direction_deviation_deg:.3f} deg, searching indices up to "
                f"{self.max_index}."
            )
        )
        return (
            f"Ideal orientation of {self.phase_name}: {self.label}, meaning the plane "
            "lies in the plane normal to the specimen axis "
            f"{self.specimen_plane_normal.tolist()} and the direction points along "
            f"{self.specimen_direction.tolist()}.{four_index} {verdict}"
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": IDEAL_ORIENTATION_SCHEMA,
            "phase": self.phase_name,
            "hkl": [int(value) for value in self.hkl],
            "uvw": [int(value) for value in self.uvw],
            "hkil": None if self.hkil is None else [int(v) for v in self.hkil],
            "uvtw": None if self.uvtw is None else [int(v) for v in self.uvtw],
            "label": self.label,
            "plane_deviation_deg": self.plane_deviation_deg,
            "direction_deviation_deg": self.direction_deviation_deg,
            "is_exact": self.is_exact,
            "max_index": self.max_index,
            "specimen_plane_normal": self.specimen_plane_normal.tolist(),
            "specimen_direction": self.specimen_direction.tolist(),
        }


def _integer_index_grid(max_index: int) -> np.ndarray:
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1)
    rows = grid.reshape(-1, 3)
    return rows[np.any(rows != 0, axis=1)]


def _nearest_integer_indices(
    target_cartesian: np.ndarray, basis: np.ndarray, max_index: int
) -> tuple[np.ndarray, float]:
    """The integer index triple whose Cartesian image is closest in angle.

    Exhaustive over the index cube, which is a few thousand rows and one
    matrix product — cheaper and far more robust than continued-fraction
    rationalization of three coupled components, which can return a large
    denominator for a direction that has a small, obvious label.
    """

    candidates = _integer_index_grid(max_index)
    cartesian = candidates.astype(np.float64) @ basis.T
    norms = np.linalg.norm(cartesian, axis=1)
    cosines = (cartesian @ target_cartesian) / (
        norms * float(np.linalg.norm(target_cartesian))
    )
    best = int(np.argmax(cosines))
    # From the vectors, not from their cosine: the interesting label is the one
    # that is exact, and arccos cannot report zero for it. See pytex.core._angles.
    deviation = float(
        np.degrees(
            angle_between_unit_vectors_rad(
                cartesian[best] / norms[best],
                target_cartesian / float(np.linalg.norm(target_cartesian)),
            )
        )
    )
    reduced = reduce_indices(candidates[best][None, :])[0]
    return np.asarray(reduced, dtype=np.int64), deviation


def ideal_orientation_indices(
    orientation: Orientation,
    *,
    specimen_plane_normal: str | ArrayLike = "ND",
    specimen_direction: str | ArrayLike = "RD",
    max_index: int = 6,
) -> IdealOrientationIndices:
    """Name an orientation as a ``(hkl)[uvw]`` texture component.

    What it does
        Maps the specimen normal and the specimen reference direction back into
        the crystal, expresses them in the reciprocal and direct bases
        respectively, and finds the integer index triples closest in angle,
        reporting both the labels and the residual angles.

    When to use it
        To report a measured or computed orientation in the language rolling and
        recrystallization texture is written in — ``{110}<112>`` brass,
        ``{112}<111>`` copper, ``{001}<100>`` cube — and to check how far a
        real orientation sits from the ideal component it is being called.

    Parameters
    ----------
    orientation:
        Must carry a phase; the indices are meaningless without a lattice.
    specimen_plane_normal, specimen_direction:
        The specimen axes the label refers to, by name (``"ND"``, ``"RD"``,
        ``"TD"``, ``"x"``, ``"y"``, ``"z"``) or as vectors. The defaults give
        the usual sheet convention, in which ``(hkl)`` lies in the sheet plane
        and ``[uvw]`` points along the rolling direction.
    max_index:
        Largest absolute index searched. Six covers every component in the
        standard rolling-texture tables; raising it can only reduce the reported
        deviations, and will eventually name any orientation with large,
        meaningless indices — which is why the deviation is reported alongside.

    Returns
    -------
    IdealOrientationIndices
        The labels, their deviations, and the four-index forms for hexagonal
        and trigonal phases.

    Raises
    ------
    ValueError
        If the orientation carries no phase, or ``max_index`` is not positive.

    Examples
    --------
    Round-tripping the copper component: build an orientation with
    ``Orientation.from_miller((1, 1, 2), (1, 1, -1), ...)`` and pass it here,
    and the returned ``hkl`` and ``uvw`` reproduce those indices with both
    deviations below a thousandth of a degree.

    See Also
    --------
    `pytex.core.orientation.Orientation.from_miller` : the inverse construction.
    """

    if orientation.phase is None:
        raise ValueError(
            "ideal_orientation_indices needs the phase: Miller indices are components "
            "in a lattice basis, and without the lattice there is nothing to index "
            "against. Build the Orientation with phase=... ."
        )
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")

    phase = orientation.phase
    matrix = np.asarray(orientation.as_matrix(), dtype=np.float64)
    normal_axis = specimen_direction_vector(specimen_plane_normal)
    direction_axis = specimen_direction_vector(specimen_direction)

    plane_normal_crystal = matrix.T @ normal_axis
    direction_crystal = matrix.T @ direction_axis

    reciprocal_basis = as_float_array(
        phase.lattice.reciprocal_basis().matrix, shape=(3, 3)
    )
    direct_basis = as_float_array(phase.lattice.direct_basis().matrix, shape=(3, 3))
    hkl, plane_deviation = _nearest_integer_indices(
        plane_normal_crystal, reciprocal_basis, max_index
    )
    uvw, direction_deviation = _nearest_integer_indices(
        direction_crystal, direct_basis, max_index
    )

    hkil: np.ndarray | None = None
    uvtw: np.ndarray | None = None
    if is_hexagonal_phase(phase):
        hkil = np.asarray(plane_hkl_to_hkil_array(hkl[None, :])[0], dtype=np.int64)
        uvtw = np.asarray(direction_uvw_to_uvtw_array(uvw[None, :])[0], dtype=np.int64)

    return IdealOrientationIndices(
        hkl=hkl,
        uvw=uvw,
        plane_deviation_deg=plane_deviation,
        direction_deviation_deg=direction_deviation,
        max_index=int(max_index),
        specimen_plane_normal=normal_axis,
        specimen_direction=direction_axis,
        phase_name=phase.name,
        hkil=hkil,
        uvtw=uvtw,
    )


# --------------------------------------------------------------------------- #
# The all-representations report
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class OrientationRepresentations:
    """One rotation, written every way at once.

    Purpose
    -------
    The object behind "I have Euler angles; show me everything else". It holds
    the same rotation in all ten numerical representations, plus — when a phase
    is available — the ``(hkl)[uvw]`` component name, and it can print itself as
    convention-explicit prose or as a table.

    When to use
    -----------
    In teaching material, in reports, and whenever a value has to cross a tool
    boundary and the receiving convention is not yet known. For computation
    inside a loop, use the free functions: this object computes all ten forms
    eagerly, which is the right trade for reporting and the wrong one for a hot
    path.

    Attributes
    ----------
    matrix : np.ndarray
        ``(3, 3)``, active convention. For an `Orientation` this is ``g``,
        crystal-to-specimen.
    quaternion : np.ndarray
        ``(4,)`` in ``(w, x, y, z)`` order, scalar part non-negative.
    axis : np.ndarray
        ``(3,)`` unit rotation axis. For the identity it is reported as
        ``[0, 0, 1]`` by convention, and ``angle_deg`` is zero.
    angle_deg : float
        Rotation angle in ``[0, 180]`` degrees.
    rodrigues : np.ndarray
        ``(3,)``; infinite at 180 degrees, which is why the Frank form exists.
    rodrigues_frank : np.ndarray
        ``(4,)`` homogeneous form. The magnitude is a projective coordinate, so
        the 180-degree rotation is exactly representable and exactly invertible
        here, where the 3-vector overflows and loses its axis.
    euler_bunge_deg, euler_matthies_deg : np.ndarray
        ``(3,)`` triples, degrees.
    homochoric, cubochoric : np.ndarray
        ``(3,)`` equal-volume coordinates.
    ideal_indices : IdealOrientationIndices or None
        Present when the source was an `Orientation` carrying a phase.
    phase_name, specimen_frame_name, crystal_frame_name : str or None
        Recorded so the numbers cannot be read against the wrong frames.
    provenance : ProvenanceRecord or None
    """

    matrix: np.ndarray
    quaternion: np.ndarray
    axis: np.ndarray
    angle_deg: float
    rodrigues: np.ndarray
    rodrigues_frank: np.ndarray
    euler_bunge_deg: np.ndarray
    euler_matthies_deg: np.ndarray
    homochoric: np.ndarray
    cubochoric: np.ndarray
    ideal_indices: IdealOrientationIndices | None = None
    phase_name: str | None = None
    specimen_frame_name: str | None = None
    crystal_frame_name: str | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", as_float_array(self.matrix, shape=(3, 3)))
        object.__setattr__(self, "quaternion", as_float_array(self.quaternion, shape=(4,)))
        for name in ("axis", "rodrigues", "homochoric", "cubochoric"):
            object.__setattr__(self, name, as_float_array(getattr(self, name), shape=(3,)))
        object.__setattr__(
            self, "rodrigues_frank", as_float_array(self.rodrigues_frank, shape=(4,))
        )
        for name in ("euler_bunge_deg", "euler_matthies_deg"):
            object.__setattr__(self, name, as_float_array(getattr(self, name), shape=(3,)))

    @property
    def as_rotation(self) -> Rotation:
        """The rotation these numbers all denote."""

        return Rotation(quaternion=self.quaternion)

    def to_table(self) -> str:
        """A fixed-width table of every representation, for a terminal or a notebook.

        The units and the convention are written into the row labels rather than
        assumed, so a table pasted into a report stays unambiguous.
        """

        rows: list[tuple[str, str]] = [
            ("matrix row 1", _format_row(self.matrix[0])),
            ("matrix row 2", _format_row(self.matrix[1])),
            ("matrix row 3", _format_row(self.matrix[2])),
            ("quaternion (w, x, y, z)", _format_row(self.quaternion)),
            ("axis (unit)", _format_row(self.axis)),
            ("angle (deg)", f"{self.angle_deg:.6f}"),
            ("Rodrigues", _format_row(self.rodrigues)),
            ("Rodrigues-Frank", _format_row(self.rodrigues_frank)),
            ("Euler Bunge ZXZ (deg)", _format_row(self.euler_bunge_deg)),
            ("Euler ZYZ (deg)", _format_row(self.euler_matthies_deg)),
            ("homochoric", _format_row(self.homochoric)),
            ("cubochoric", _format_row(self.cubochoric)),
        ]
        if self.ideal_indices is not None:
            rows.append(("ideal (hkl)[uvw]", self.ideal_indices.label))
        width = max(len(label) for label, _ in rows)
        return "\n".join(f"{label:<{width}}  {value}" for label, value in rows)

    def describe(self) -> str:
        """Convention-explicit prose: what the rotation is, in every form."""

        frames = ""
        if self.specimen_frame_name and self.crystal_frame_name:
            frames = (
                " The convention is crystal-to-specimen: the matrix carries a vector "
                f"expressed in the crystal frame '{self.crystal_frame_name}' into the "
                f"specimen frame '{self.specimen_frame_name}'."
            )
        elif self.phase_name is None:
            frames = (
                " No frames are attached: this is a rotation, not an orientation, so "
                "the matrix acts within a single frame under the active convention "
                "v' = R v."
            )
        phi1, capital_phi, phi2 = (float(value) for value in self.euler_bunge_deg)
        alpha, beta, gamma = (float(value) for value in self.euler_matthies_deg)
        head = (
            f"A rotation of {self.angle_deg:.4f} deg about the axis "
            f"[{self.axis[0]:.4f}, {self.axis[1]:.4f}, {self.axis[2]:.4f}]."
            f"{frames}"
        )
        euler = (
            f" Bunge ZXZ (phi1, Phi, phi2) = ({phi1:.4f}, {capital_phi:.4f}, "
            f"{phi2:.4f}) deg; the same rotation in ZYZ (alpha, beta, gamma) = "
            f"({alpha:.4f}, {beta:.4f}, {gamma:.4f}) deg. The two triples are "
            "different numbers for one rotation, not different rotations."
        )
        quaternion = (
            f" Quaternion (w, x, y, z) = ({self.quaternion[0]:.6f}, "
            f"{self.quaternion[1]:.6f}, {self.quaternion[2]:.6f}, "
            f"{self.quaternion[3]:.6f}), reported with a non-negative scalar part "
            "because q and -q are the same rotation."
        )
        homochoric_norm = float(np.linalg.norm(self.homochoric))
        equal_volume = (
            f" Equal-volume coordinates: homochoric norm {homochoric_norm:.6f} "
            f"of a maximum {HOMOCHORIC_BALL_RADIUS:.6f}, and the cubochoric coordinate "
            f"[{self.cubochoric[0]:.6f}, {self.cubochoric[1]:.6f}, "
            f"{self.cubochoric[2]:.6f}] inside a cube of half-edge "
            f"{CUBOCHORIC_CUBE_HALF_EDGE:.6f}. Uniform sampling of that cube is uniform "
            "sampling of orientation space; uniform sampling of Euler angles is not."
        )
        ideal = f" {self.ideal_indices.describe()}" if self.ideal_indices else ""
        return head + euler + quaternion + equal_volume + ideal

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": ORIENTATION_REPRESENTATIONS_SCHEMA,
            "phase": self.phase_name,
            "specimen_frame": self.specimen_frame_name,
            "crystal_frame": self.crystal_frame_name,
            "matrix": self.matrix.tolist(),
            "quaternion_wxyz": self.quaternion.tolist(),
            "axis": self.axis.tolist(),
            "angle_deg": self.angle_deg,
            "rodrigues": self.rodrigues.tolist(),
            "rodrigues_frank": self.rodrigues_frank.tolist(),
            "euler_bunge_deg": self.euler_bunge_deg.tolist(),
            "euler_zyz_deg": self.euler_matthies_deg.tolist(),
            "homochoric": self.homochoric.tolist(),
            "cubochoric": self.cubochoric.tolist(),
            "ideal_orientation": (
                None if self.ideal_indices is None else self.ideal_indices.to_json_dict()
            ),
        }


def _format_row(values: np.ndarray) -> str:
    return "  ".join(f"{float(value):+.6f}" for value in np.ravel(values))


def canonical_quaternions(quaternions: ArrayLike) -> np.ndarray:
    """Pick the ``w >= 0`` representative of each quaternion.

    What it does
        Negates any row whose scalar part is negative, and — when the scalar
        part is exactly zero, the 180-degree rotations, where both signs have
        ``w = 0`` — resolves the remaining tie on the first non-zero vector
        component.

    When to use it
        Before *reporting* or *comparing* quaternions componentwise. ``q`` and
        ``-q`` are the same rotation, so a raw componentwise comparison of two
        equal rotations can report a difference of 2. Composition and rotation
        algebra do not need this and are unaffected by it.

    Parameters
    ----------
    quaternions:
        ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order.

    Returns
    -------
    np.ndarray
        ``(n, 4)`` with the canonical sign chosen row by row.
    """

    values = normalize_quaternions(quaternions)
    if values.size == 0:
        return values
    negative = values[:, 0] < 0.0
    tied = np.isclose(values[:, 0], 0.0, atol=1e-15)
    if np.any(tied):
        leading = np.zeros(values.shape[0], dtype=np.float64)
        for column in (1, 2, 3):
            unresolved = leading == 0.0
            leading = np.where(unresolved, values[:, column], leading)
        negative = np.where(tied, leading < 0.0, negative)
    return np.ascontiguousarray(np.where(negative[:, None], -values, values))


def _representation_fields(quaternion: np.ndarray) -> dict[str, Any]:
    batch = canonical_quaternions(quaternion[None, :])
    quaternion = batch[0]
    axes, angles = quaternions_to_axes_angles(batch)
    return {
        "matrix": quaternions_to_matrices(batch)[0],
        "quaternion": quaternion,
        "axis": axes[0],
        "angle_deg": float(np.degrees(angles[0])),
        "rodrigues": quaternions_to_rodrigues(batch, frank=False)[0],
        "rodrigues_frank": quaternions_to_rodrigues(batch, frank=True)[0],
        "euler_bunge_deg": quaternions_to_euler_angles(batch, convention="bunge")[0],
        "euler_matthies_deg": quaternions_to_euler_angles(batch, convention="matthies")[0],
        "homochoric": homochoric_from_quaternions(batch)[0],
        "cubochoric": cubochoric_from_quaternions(batch)[0],
    }


def rotation_representations(rotation: Rotation) -> OrientationRepresentations:
    """Every numerical form of a bare `Rotation`, in one call.

    Use this for a rotation that is not attached to a crystal — a stage tilt, a
    frame relabelling, a misorientation axis-angle pair — where there is no
    phase and therefore no ``(hkl)[uvw]`` to report.

    See Also
    --------
    orientation_representations : the crystal-attached form.
    """

    quaternion = as_float_array(rotation.quaternion, shape=(4,))
    return OrientationRepresentations(
        provenance=rotation.provenance,
        **_representation_fields(quaternion),
    )


def orientation_representations(
    orientation: Orientation,
    *,
    specimen_plane_normal: str | ArrayLike = "ND",
    specimen_direction: str | ArrayLike = "RD",
    max_index: int = 6,
    include_ideal_indices: bool = True,
) -> OrientationRepresentations:
    """Every numerical form of a crystal orientation, in one call.

    What it does
        Converts the orientation once into each of the ten representations, and
        additionally names it as a ``(hkl)[uvw]`` texture component when the
        orientation carries a phase.

    When to use it
        Whenever an orientation has to be *communicated* rather than computed
        with: reporting a measurement, comparing a PyTex result against a value
        quoted in another convention, or teaching what the representations mean.
        The returned object prints itself with :meth:`~OrientationRepresentations.describe`
        and serializes with :meth:`~OrientationRepresentations.to_json_dict`.

    Parameters
    ----------
    orientation:
        The orientation to describe. Its matrix is ``g``, crystal-to-specimen.
    specimen_plane_normal, specimen_direction:
        The specimen axes the ``(hkl)[uvw]`` label refers to; see
        :func:`ideal_orientation_indices`.
    max_index:
        Index-search bound for the component name.
    include_ideal_indices:
        Set ``False`` to skip the index search, which is the only part of the
        call that costs more than microseconds.

    Returns
    -------
    OrientationRepresentations

    Examples
    --------
    The cube component of a cubic phase is the identity orientation, so its
    Bunge angles are all zero, its rotation angle is zero, its homochoric and
    cubochoric coordinates are the origin, and its ideal indices come back as
    ``(001)[100]``.
    """

    quaternion = as_float_array(orientation.rotation.quaternion, shape=(4,))
    ideal: IdealOrientationIndices | None = None
    if include_ideal_indices and orientation.phase is not None:
        ideal = ideal_orientation_indices(
            orientation,
            specimen_plane_normal=specimen_plane_normal,
            specimen_direction=specimen_direction,
            max_index=max_index,
        )
    return OrientationRepresentations(
        ideal_indices=ideal,
        phase_name=None if orientation.phase is None else orientation.phase.name,
        specimen_frame_name=orientation.specimen_frame.name,
        crystal_frame_name=orientation.crystal_frame.name,
        provenance=orientation.provenance,
        **_representation_fields(quaternion),
    )


@dataclass(frozen=True, slots=True)
class OrientationRepresentationSet:
    """A batch of rotations in every representation, computed vectorized.

    Purpose
    -------
    The batch counterpart of :class:`OrientationRepresentations`. Every array is
    computed once for the whole batch through the vectorized free functions, so
    describing ten thousand orientations costs ten array passes rather than ten
    thousand object constructions.

    When to use
    -----------
    Exporting an orientation cloud in several conventions at once; building the
    comparison tables a validation note needs; feeding a plotting routine that
    wants Rodrigues vectors and a reporting routine that wants Euler angles from
    the same data without converting twice.

    Attributes
    ----------
    quaternions : np.ndarray
        ``(n, 4)``, the storage form everything else is derived from.
    matrices : np.ndarray
        ``(n, 3, 3)``.
    axes : np.ndarray
        ``(n, 3)`` unit axes.
    angles_deg : np.ndarray
        ``(n,)`` in ``[0, 180]``.
    rodrigues, homochoric, cubochoric : np.ndarray
        ``(n, 3)``.
    rodrigues_frank : np.ndarray
        ``(n, 4)``.
    euler_bunge_deg, euler_matthies_deg : np.ndarray
        ``(n, 3)``.
    provenance : ProvenanceRecord or None
    """

    quaternions: np.ndarray
    matrices: np.ndarray
    axes: np.ndarray
    angles_deg: np.ndarray
    rodrigues: np.ndarray
    rodrigues_frank: np.ndarray
    euler_bunge_deg: np.ndarray
    euler_matthies_deg: np.ndarray
    homochoric: np.ndarray
    cubochoric: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    @classmethod
    def from_quaternions(
        cls,
        quaternions: ArrayLike,
        *,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRepresentationSet:
        """Build the full representation batch from ``(n, 4)`` quaternions.

        The quaternions are put in the canonical ``w >= 0`` sign, matching the
        single-orientation report, so that a row of the batch and the report for
        the same rotation agree component for component and not merely up to a
        sign.
        """

        normalized = canonical_quaternions(quaternions)
        axes, angles = quaternions_to_axes_angles(normalized)
        return cls(
            quaternions=normalized,
            matrices=quaternions_to_matrices(normalized),
            axes=axes,
            angles_deg=np.degrees(angles),
            rodrigues=quaternions_to_rodrigues(normalized, frank=False),
            rodrigues_frank=quaternions_to_rodrigues(normalized, frank=True),
            euler_bunge_deg=quaternions_to_euler_angles(normalized, convention="bunge"),
            euler_matthies_deg=quaternions_to_euler_angles(
                normalized, convention="matthies"
            ),
            homochoric=homochoric_from_quaternions(normalized),
            cubochoric=cubochoric_from_quaternions(normalized),
            provenance=provenance,
        )

    @classmethod
    def from_values(
        cls,
        values: ArrayLike,
        *,
        source: RepresentationKind | str,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRepresentationSet:
        """Build the batch from any single representation.

        The bulk analogue of "give me the other nine forms": pass what you have,
        name it, and every other representation comes back computed.
        """

        quaternions = _to_quaternions(values, RepresentationKind(str(source)))
        return cls.from_quaternions(quaternions, provenance=provenance)

    def row(self, index: int) -> OrientationRepresentations:
        """The single-orientation report for one row of the batch."""

        return OrientationRepresentations(
            matrix=self.matrices[index],
            quaternion=self.quaternions[index],
            axis=self.axes[index],
            angle_deg=float(self.angles_deg[index]),
            rodrigues=self.rodrigues[index],
            rodrigues_frank=self.rodrigues_frank[index],
            euler_bunge_deg=self.euler_bunge_deg[index],
            euler_matthies_deg=self.euler_matthies_deg[index],
            homochoric=self.homochoric[index],
            cubochoric=self.cubochoric[index],
            provenance=self.provenance,
        )

    def describe(self) -> str:
        """Convention-explicit prose summarising the batch."""

        if len(self) == 0:
            return "An empty orientation batch: no representations to report."
        angles = np.asarray(self.angles_deg, dtype=np.float64)
        return (
            f"{len(self)} rotations in all ten representations. Rotation angles span "
            f"{float(angles.min()):.3f} to {float(angles.max()):.3f} deg with a mean of "
            f"{float(angles.mean()):.3f} deg. Homochoric norms reach "
            f"{float(np.max(np.linalg.norm(self.homochoric, axis=1))):.4f} of the "
            f"ball radius {HOMOCHORIC_BALL_RADIUS:.4f}, and the cubochoric coordinates "
            f"reach {float(np.max(np.abs(self.cubochoric))):.4f} of the cube half-edge "
            f"{CUBOCHORIC_CUBE_HALF_EDGE:.4f}. Every array is derived from the same "
            "quaternion batch, so the representations cannot disagree with each other."
        )
