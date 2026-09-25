"""The normative geometry of a vertically milled lamella: section 6 of the foundation.

Purpose
-------
Answer, for every orientation at once, the one crystallographic question FIB
planning asks: *how far is the nearest symmetry-equivalent of the target zone
axis from the surface plane, and which way must the lamella face to put it on
the beam?*

The physical constraint
-----------------------
A strictly vertical mill leaves the two faces of the slab perpendicular to the
surface, so the lamella normal ``n_L`` is forced into the surface plane. In the
TEM the beam enters along ``n_L``. The achievable beam directions are therefore
the great circle of in-plane directions, and a grain is usable for a zone axis
exactly when some equivalent of it lies close to that circle. The angular
distance from the circle,

    eps(u) = arcsin |d_S(u) . Z_s|,

**is** the residual tilt the TEM holder must supply.

Orientation convention
----------------------
PyTex orientations map crystal to specimen: ``v_specimen = g v_crystal``
(:class:`pytex.core.orientation.Orientation`). The sample-frame image of a
crystal direction ``u`` is therefore

    d_S = M g u,

with ``M`` the specimen-to-sample rotation of
:class:`pytex.fib.frames.SurfaceGeometry`. The foundation document wrote the
step as ``g^-1 u`` for the opposite convention and asked for it to be checked;
it is checked in ``tests/unit/test_fib_geometry.py`` by a case that fails under
the transposed matrix.

Steps
-----
1. orbit: ``U = {s_i u} union {-U}`` over the point group (decision D8);
2. to sample: ``d_S = M g u`` for every member;
3. out-of-plane angle ``eps = arcsin |d_S . Z_s|``;
4. best member ``u* = argmin eps`` (ties broken deterministically, below);
5. lamella normal ``n_L = normalize(d* - (d* . Z_s) Z_s)``;
6. long axis ``t_L = Z_s x n_L``;
7. azimuth ``theta_S = atan2(n_L . Y_s, n_L . X_s)``.

Steps 1-4 are vectorized over orientations and orbit members with ``einsum``;
nothing here loops over grains.

Canonical choices
-----------------
A plane normal and its reverse describe the same lamella, so ``n_L`` is reported
with ``theta_S`` in ``[0, 180)``. Among members whose ``eps`` ties within
:data:`TIE_TOLERANCE_DEG`, the one with the smallest canonical ``theta_S`` is
taken, and among those the one pointing along ``+n_L`` and then ``+Z_s``. These
rules depend only on the set of options, not on the order of the operators, so
the answer is invariant under a symmetry operator applied to ``g``.

See ``docs/site/theory/fib_lamella_zone_axis_geometry.md`` for the derivation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.lattice import CrystalDirection, Phase, ZoneAxis
from pytex.core.notation import format_direction_family_indices, format_direction_indices
from pytex.fib.frames import SurfaceGeometry, wrap_azimuth_deg

__all__ = [
    "TIE_TOLERANCE_DEG",
    "LamellaGeometry",
    "LamellaGeometryBatch",
    "LamellaOption",
    "TargetOrbit",
    "angle_to_plane_deg",
    "lamella_geometry",
    "lamella_geometry_batch",
    "lamella_options",
    "target_orbit",
]

#: Members whose out-of-plane angles differ by less than this are ties.
TIE_TOLERANCE_DEG = 1e-7

#: Below this in-plane magnitude the lamella normal is undetermined: the target
#: lies along the surface normal and no vertical lamella can view it.
_DEGENERATE_IN_PLANE = 1e-12


def _coerce_indices(target: ZoneAxis | CrystalDirection | ArrayLike) -> np.ndarray:
    if isinstance(target, ZoneAxis):
        values = np.asarray(target.indices)
    elif isinstance(target, CrystalDirection):
        values = np.asarray(target.coordinates)
    else:
        values = np.asarray(target)
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("A target zone axis needs exactly three finite indices [uvw].")
    rounded = np.rint(values)
    if not np.allclose(values, rounded, atol=1e-9):
        raise ValueError("Zone-axis indices [uvw] must be whole numbers.")
    if not np.any(rounded):
        raise ValueError("The zero triplet [000] is not a direction.")
    return rounded.astype(np.int64)


@dataclass(frozen=True, slots=True)
class TargetOrbit:
    """Every crystallographically acceptable version of a target zone axis.

    Purpose
    -------
    Decision D8: any symmetry-equivalent of ``<uvw>`` is acceptable, in either
    sense. The orbit is built once per phase and target and shared by every
    orientation of the map.

    Attributes
    ----------
    phase : Phase
    target_indices : (3,) int array
        The requested ``[uvw]``.
    indices : (m, 3) int array
        Direct-lattice indices of every distinct member.
    cartesian : (m, 3) float array
        Unit crystal-frame Cartesian vectors of the same members.
    both_senses : bool
        Whether ``-U`` was added. True by default (D8); set false only for work
        in which the beam sense is known to matter and has been determined.
    """

    phase: Phase
    target_indices: np.ndarray
    indices: np.ndarray
    cartesian: np.ndarray
    both_senses: bool = True

    def __post_init__(self) -> None:
        for name in ("target_indices", "indices", "cartesian"):
            array = np.ascontiguousarray(getattr(self, name))
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if self.indices.ndim != 2 or self.indices.shape[1] != 3:
            raise ValueError("TargetOrbit.indices must be (m, 3).")
        if self.cartesian.shape != self.indices.shape:
            raise ValueError("TargetOrbit.cartesian must match indices.")

    def __len__(self) -> int:
        return int(self.indices.shape[0])

    @property
    def family_text(self) -> str:
        """The target written as a family ``<uvw>``, through :mod:`pytex.core.notation`."""

        return format_direction_family_indices(
            [int(value) for value in self.target_indices], style="plain"
        )

    def member_text(self, index: int) -> str:
        """Orbit member ``index`` written as a specific direction ``[uvw]``."""

        return format_direction_indices(
            [int(value) for value in self.indices[index]], style="plain"
        )


def target_orbit(
    phase: Phase,
    target: ZoneAxis | CrystalDirection | ArrayLike,
    *,
    both_senses: bool = True,
) -> TargetOrbit:
    """Build the symmetry orbit of a target zone axis.

    Parameters
    ----------
    phase : Phase
        Supplies the lattice (indices to Cartesian) and the point group.
    target : ZoneAxis, CrystalDirection or array_like
        The wanted ``[uvw]``, in direct-lattice indices.
    both_senses : bool, default True
        Add the reversed members (decision D8).

    Returns
    -------
    TargetOrbit

    Notes
    -----
    Member indices are computed exactly as ``D^-1 s_i D u`` with ``D`` the
    direct-basis matrix, which is an integer matrix for every operator of a
    point group expressed in its own lattice, so no rounding tolerance is
    involved beyond removing floating-point noise.

    Examples
    --------
    ``<100>`` of a cubic crystal has six members counting both senses.
    """

    indices = _coerce_indices(target)
    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=np.float64)
    inverse = np.linalg.inv(direct)
    operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
    cartesian = direct @ indices.astype(np.float64)
    images = np.einsum("nij,j->ni", operators, cartesian)
    if both_senses:
        images = np.vstack([images, -images])
    lattice_images = images @ inverse.T
    rounded = np.rint(lattice_images)
    if not np.allclose(lattice_images, rounded, atol=1e-6):
        raise ValueError(
            "The point-group operators do not map the lattice onto itself; the phase's symmetry "
            "and lattice settings are inconsistent."
        )
    members = np.unique(rounded.astype(np.int64), axis=0)
    # A stable, readable order: positive leading entries first, then by indices.
    order = np.lexsort(tuple(-members[:, column] for column in range(2, -1, -1)))
    members = members[order]
    vectors = members.astype(np.float64) @ direct.T
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    return TargetOrbit(
        phase=phase,
        target_indices=indices,
        indices=members,
        cartesian=vectors,
        both_senses=both_senses,
    )


@dataclass(frozen=True, slots=True)
class LamellaGeometry:
    """The solved geometry for one orientation.

    Attributes
    ----------
    eps_deg : float
        ``eps*``: out-of-plane angle of the best member, the residual tilt the
        holder must supply.
    member_index : int
        Row of :attr:`TargetOrbit.indices` chosen.
    member_indices : (3,) int array
        ``[uvw]`` of that member.
    direction_sample : (3,) array
        ``d*``, its unit image in the sample frame.
    normal_sample : (3,) array
        ``n_L`` in sample components, in the surface, azimuth in ``[0, 180)``.
    long_axis_sample : (3,) array
        ``t_L = Z_s x n_L``.
    theta_sample_deg : float
        ``theta_S``: azimuth of ``n_L`` from ``X_s`` towards ``Y_s``.
    tied_members : int
        How many members share ``eps*`` within :data:`TIE_TOLERANCE_DEG`
        (at least one, the chosen member; a member and its reverse count as
        one lamella option each).
    degenerate : bool
        True when every member lies along the surface normal (``eps* = 90``),
        so no vertical lamella can view the target and ``n_L`` is arbitrary.
    """

    eps_deg: float
    member_index: int
    member_indices: np.ndarray
    direction_sample: np.ndarray
    normal_sample: np.ndarray
    long_axis_sample: np.ndarray
    theta_sample_deg: float
    tied_members: int = 1
    degenerate: bool = False

    def __post_init__(self) -> None:
        for name in ("member_indices", "direction_sample", "normal_sample", "long_axis_sample"):
            array = np.ascontiguousarray(getattr(self, name))
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @property
    def long_axis_azimuth_deg(self) -> float:
        """Azimuth of ``t_L``, which is ``theta_S + 90`` wrapped to ``[0, 180)``."""

        return round(float(wrap_azimuth_deg(self.theta_sample_deg + 90.0, period=180.0)), 12)

    @property
    def out_of_plane_sign(self) -> int:
        """``+1`` when ``d*`` rises out of the surface (towards ``+Z_s``), else ``-1``."""

        return 1 if float(self.direction_sample[2]) >= 0.0 else -1

    def sample_to_lamella_matrix(self) -> np.ndarray:
        """Rows ``n_L, t_L, Z_s``: the rotation taking sample to lamella components."""

        return np.vstack(
            [self.normal_sample, self.long_axis_sample, np.array([0.0, 0.0, 1.0])]
        ).astype(np.float64)

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "eps_deg": float(self.eps_deg),
            "member_indices": [int(value) for value in self.member_indices],
            "direction_sample": [float(value) for value in self.direction_sample],
            "normal_sample": [float(value) for value in self.normal_sample],
            "long_axis_sample": [float(value) for value in self.long_axis_sample],
            "theta_sample_deg": float(self.theta_sample_deg),
            "long_axis_azimuth_deg": self.long_axis_azimuth_deg,
            "tied_members": int(self.tied_members),
            "degenerate": bool(self.degenerate),
        }


@dataclass(frozen=True, slots=True)
class LamellaGeometryBatch:
    """The geometry for many orientations, as arrays.

    Attributes
    ----------
    eps_deg : (n,) array
    member_index : (n,) int array
    direction_sample, normal_sample, long_axis_sample : (n, 3) arrays
    theta_sample_deg : (n,) array
    tied_members : (n,) int array
    degenerate : (n,) bool array
    member_eps_deg : (n, m) array
        ``eps`` of every orbit member, for the stereogram and for ranking
        secondary options.
    """

    eps_deg: np.ndarray
    member_index: np.ndarray
    direction_sample: np.ndarray
    normal_sample: np.ndarray
    long_axis_sample: np.ndarray
    theta_sample_deg: np.ndarray
    tied_members: np.ndarray
    degenerate: np.ndarray
    member_eps_deg: np.ndarray

    def __len__(self) -> int:
        return int(self.eps_deg.shape[0])

    def item(self, index: int, orbit: TargetOrbit) -> LamellaGeometry:
        """The single-orientation view of row ``index``."""

        member = int(self.member_index[index])
        return LamellaGeometry(
            eps_deg=float(self.eps_deg[index]),
            member_index=member,
            member_indices=np.asarray(orbit.indices[member]),
            direction_sample=np.asarray(self.direction_sample[index]),
            normal_sample=np.asarray(self.normal_sample[index]),
            long_axis_sample=np.asarray(self.long_axis_sample[index]),
            theta_sample_deg=float(self.theta_sample_deg[index]),
            tied_members=int(self.tied_members[index]),
            degenerate=bool(self.degenerate[index]),
        )


def _as_matrices(orientations: ArrayLike) -> np.ndarray:
    matrices = np.asarray(orientations, dtype=np.float64)
    if matrices.shape == (3, 3):
        matrices = matrices[None]
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3):
        raise ValueError("Orientation matrices must have shape (3, 3) or (n, 3, 3).")
    return matrices


def lamella_geometry_batch(
    orientation_matrices: ArrayLike,
    orbit: TargetOrbit,
    surface: SurfaceGeometry | None = None,
) -> LamellaGeometryBatch:
    """Solve section 6 for a stack of crystal-to-specimen matrices.

    Parameters
    ----------
    orientation_matrices : (n, 3, 3) or (3, 3) array_like
        ``g`` for each orientation, crystal to EBSD specimen, as
        :meth:`pytex.core.orientation.OrientationSet.as_matrices` returns them.
    orbit : TargetOrbit
    surface : SurfaceGeometry, optional
        Defaults to :class:`SurfaceGeometry`.

    Returns
    -------
    LamellaGeometryBatch

    Notes
    -----
    Cost is ``O(n m)`` in memory and time for ``n`` orientations and ``m``
    orbit members (``m <= 48``), with no Python loop over either.
    """

    surface = surface or SurfaceGeometry()
    matrices = _as_matrices(orientation_matrices)
    sample = surface.specimen_to_sample_matrix()
    to_sample = np.einsum("ij,njk->nik", sample, matrices)
    directions = np.einsum("nik,mk->nmi", to_sample, orbit.cartesian)  # (n, m, 3)
    out_of_plane = np.clip(np.abs(directions[..., 2]), 0.0, 1.0)
    member_eps = np.degrees(np.arcsin(out_of_plane))

    in_plane = directions[..., :2]
    in_plane_norm = np.linalg.norm(in_plane, axis=-1)
    safe = np.where(in_plane_norm > _DEGENERATE_IN_PLANE, in_plane_norm, 1.0)
    unit = in_plane / safe[..., None]
    raw_azimuth = np.degrees(np.arctan2(unit[..., 1], unit[..., 0]))
    canonical_azimuth = wrap_azimuth_deg(raw_azimuth, period=180.0)
    canonical_azimuth = np.where(in_plane_norm > _DEGENERATE_IN_PLANE, canonical_azimuth, 0.0)
    # Sense of the member relative to its canonical normal: +1 along +n_L.
    canonical_rad = np.radians(canonical_azimuth)
    along_normal = unit[..., 0] * np.cos(canonical_rad) + unit[..., 1] * np.sin(canonical_rad)

    minimum = member_eps.min(axis=1, keepdims=True)
    tied = member_eps <= minimum + TIE_TOLERANCE_DEG
    # Lexicographic key over tied members: canonical azimuth, then +n_L sense,
    # then +Z_s sense, then orbit row. Encoded as one float per member so the
    # choice is a single argmin, with the azimuth rounded to absorb noise.
    azimuth_key = np.round(canonical_azimuth, 9)
    sense_key = np.where(along_normal >= 0.0, 0.0, 1.0)
    rise_key = np.where(directions[..., 2] >= 0.0, 0.0, 1.0)
    member_rank = np.arange(orbit.cartesian.shape[0], dtype=np.float64)[None, :]
    key = azimuth_key * 1e4 + sense_key * 1e3 + rise_key * 1e2 + member_rank * 1e-3
    key = np.where(tied, key, np.inf)
    chosen = np.argmin(key, axis=1)

    rows = np.arange(matrices.shape[0])
    best_direction = directions[rows, chosen]
    eps = member_eps[rows, chosen]
    degenerate = in_plane_norm[rows, chosen] <= _DEGENERATE_IN_PLANE
    theta = canonical_azimuth[rows, chosen]
    theta_rad = np.radians(theta)
    normal = np.stack([np.cos(theta_rad), np.sin(theta_rad), np.zeros_like(theta_rad)], axis=-1)
    long_axis = np.stack([-np.sin(theta_rad), np.cos(theta_rad), np.zeros_like(theta_rad)], axis=-1)

    # Distinct lamella options among the ties: members differing only in sense
    # share an azimuth and a lamella, so count distinct azimuths.
    tied_azimuths = np.where(tied, azimuth_key, np.nan)
    sorted_azimuths = np.sort(tied_azimuths, axis=1)
    distinct = np.sum(
        np.isfinite(sorted_azimuths)
        & ~np.isclose(
            sorted_azimuths,
            np.concatenate(
                [np.full((matrices.shape[0], 1), -1.0), sorted_azimuths[:, :-1]], axis=1
            ),
            atol=1e-7,
        ),
        axis=1,
    )
    return LamellaGeometryBatch(
        eps_deg=eps,
        member_index=chosen.astype(np.int64),
        direction_sample=best_direction,
        normal_sample=normal,
        long_axis_sample=long_axis,
        theta_sample_deg=theta,
        tied_members=np.maximum(distinct, 1).astype(np.int64),
        degenerate=degenerate,
        member_eps_deg=member_eps,
    )


def lamella_geometry(
    orientation_matrix: ArrayLike,
    orbit: TargetOrbit,
    surface: SurfaceGeometry | None = None,
) -> LamellaGeometry:
    """Solve section 6 for one orientation.

    Parameters
    ----------
    orientation_matrix : (3, 3) array_like
        ``g``, crystal to EBSD specimen.
    orbit : TargetOrbit
    surface : SurfaceGeometry, optional

    Returns
    -------
    LamellaGeometry

    Examples
    --------
    With ``g`` the identity and the cubic target ``<100>``, ``[100]`` lies in
    the surface, so ``eps* = 0`` and the lamella normal is ``X_s``.
    """

    batch = lamella_geometry_batch(np.asarray(orientation_matrix)[None], orbit, surface)
    return batch.item(0, orbit)


def angle_to_plane_deg(direction: ArrayLike, plane_normal: ArrayLike) -> float:
    """Angle between a direction and a plane, ``arcsin |d . n|``, in degrees.

    Used by the optional secondary score (section 7.6): a plane is edge-on to
    the beam when this angle between the beam and the plane is zero.
    """

    d = np.asarray(direction, dtype=np.float64)
    n = np.asarray(plane_normal, dtype=np.float64)
    cosine = abs(float(d @ n)) / (float(np.linalg.norm(d)) * float(np.linalg.norm(n)))
    return math.degrees(math.asin(min(1.0, cosine)))


@dataclass(frozen=True, slots=True)
class LamellaOption:
    """One distinct way to cut a lamella for the target: a member and its azimuth.

    Two members that differ only in sense give the same lamella and are one
    option. The best option is the plan; the next ones are what an operator
    can fall back on when the best azimuth is blocked by a neighbouring grain
    or a scratch.
    """

    member_indices: tuple[int, int, int]
    eps_deg: float
    theta_sample_deg: float

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "member_indices": list(self.member_indices),
            "eps_deg": float(self.eps_deg),
            "theta_sample_deg": float(self.theta_sample_deg),
        }


def lamella_options(
    orientation_matrix: ArrayLike,
    orbit: TargetOrbit,
    surface: SurfaceGeometry | None = None,
    *,
    limit: int = 4,
) -> tuple[LamellaOption, ...]:
    """The best ``limit`` distinct lamella azimuths for one orientation, best first.

    Parameters
    ----------
    orientation_matrix : (3, 3) array_like
        ``g``, crystal to EBSD specimen.
    orbit : TargetOrbit
    surface : SurfaceGeometry, optional
    limit : int, default 4

    Returns
    -------
    tuple of LamellaOption
        Sorted by ``eps``; members sharing a canonical azimuth are merged,
        keeping the one :func:`lamella_geometry_batch` would choose.
    """

    surface = surface or SurfaceGeometry()
    g = _as_matrices(orientation_matrix)[0]
    directions = orbit.cartesian @ (surface.specimen_to_sample_matrix() @ g).T
    eps = np.degrees(np.arcsin(np.clip(np.abs(directions[:, 2]), 0.0, 1.0)))
    azimuth = wrap_azimuth_deg(
        np.degrees(np.arctan2(directions[:, 1], directions[:, 0])), period=180.0
    )
    order = np.lexsort((np.round(azimuth, 9), np.round(eps, 9)))
    options: list[LamellaOption] = []
    for row in order:
        theta = float(azimuth[row])
        if any(
            abs(theta - option.theta_sample_deg) < 1e-6
            or abs(abs(theta - option.theta_sample_deg) - 180.0) < 1e-6
            for option in options
        ):
            continue
        member = tuple(int(value) for value in orbit.indices[row])
        options.append(
            LamellaOption(
                member_indices=(member[0], member[1], member[2]),
                eps_deg=float(eps[row]),
                theta_sample_deg=theta,
            )
        )
        if len(options) >= limit:
            break
    return tuple(options)
