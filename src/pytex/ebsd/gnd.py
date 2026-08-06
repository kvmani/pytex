"""Lattice curvature and geometrically necessary dislocation density.

Plastic deformation that leaves a *net* Burgers vector across a region bends the
lattice, and that bending is measurable: an EBSD map records the orientation
field, and its spatial gradient is the lattice curvature. Nye's dislocation
density tensor converts curvature into the density of dislocations that must be
present to accommodate it — the geometrically necessary dislocations (GNDs).

This is the natural completion of the local-misorientation family already in
:mod:`pytex.ebsd.models` (KAM, GROD, GOS, GAM): those report *how much* the
orientation varies locally, while GND density reports *what dislocation content
that variation implies*, in physical units of line length per unit volume.

Theory
------
For a rotation field with rotation vector :math:`\\boldsymbol{\\omega}(\\mathbf{x})`,
the lattice curvature tensor is

.. math:: \\kappa_{ij} = \\frac{\\partial \\omega_i}{\\partial x_j},

and Nye's dislocation density tensor follows as

.. math:: \\alpha_{ij} = \\kappa_{ji} - \\delta_{ij}\\kappa_{kk}.

For small elastic strains the curvature is dominated by the lattice rotation,
which is what an orientation map measures.

What a surface map cannot see
-----------------------------
A 2-D EBSD map supplies gradients along two directions only, so just **six of
the nine** components of :math:`\\kappa` are measurable, and hence only **five of
the nine** components of :math:`\\alpha`. The two diagonal components
:math:`\\alpha_{11}` and :math:`\\alpha_{22}` additionally depend on the
unmeasurable :math:`\\kappa_{33}` through the trace term.

Every GND density obtained this way is therefore a **lower bound**: dislocation
content that produces no in-plane curvature is invisible, as is all statistically
stored dislocation content, whose Burgers vectors cancel by construction. PyTex
reports the unmeasurable components as ``NaN`` rather than as zero, so the
distinction between "measured to be zero" and "not measured" survives.

Resolution dependence
---------------------
GND density from an orientation map is not a property of the material alone: it
depends on the step size. A finer step resolves sharper gradients and reports a
higher density, because sub-step-scale curvature is averaged away. The step size
must therefore be reported with any GND value, and two maps are only comparable
at the same step size. PyTex keeps the step size in the calculation explicit for
this reason.

See Also
--------
pytex.ebsd.CrystalMap.kernel_average_misorientation_deg : The local-misorientation
    measure the KAM-based estimate is built on.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

from pytex.core._arrays import FloatArray, freeze_array

if TYPE_CHECKING:  # pragma: no cover - typing only
    from pytex.ebsd.models import CrystalMap

#: Metres per map coordinate unit under the default assumption that map
#: coordinates are micrometres, which is the EBSD convention.
_DEFAULT_STEP_SCALE_M = 1e-6

#: Burgers vector magnitudes are quoted in nanometres on this surface; this is
#: the conversion to metres.
_NANOMETRE_M = 1e-9

GNDMethod = Literal["curvature", "kam"]

#: The four Nye-tensor components a 2-D surface map fully determines:
#: ``alpha_01``, ``alpha_02``, ``alpha_10``, ``alpha_12``. Because
#: ``alpha_ij = kappa_ji`` off the diagonal and the third *column* of the
#: curvature is unmeasurable, these are the off-diagonal entries with a first
#: index of 0 or 1. The fifth determinable quantity is the difference
#: ``alpha_00 - alpha_11``, which is not a component and so is not summed here.
_MEASURABLE_NYE_COMPONENTS = np.array(
    [
        [False, True, True],
        [True, False, True],
        [False, False, False],
    ],
    dtype=bool,
)


def _rotation_vectors(left_matrices: np.ndarray, right_matrices: np.ndarray) -> FloatArray:
    """Specimen-frame rotation vectors taking each left orientation to its right.

    The relative rotation is formed in the *specimen* frame as
    ``g_right g_left^T``, because the curvature is a gradient of the lattice
    rotation with respect to specimen coordinates. Its rotation vector is the
    axis scaled by the angle, extracted from the antisymmetric part, which is
    both vectorized and numerically well behaved at the small angles that
    neighbouring measurement points actually exhibit.
    """

    relative = np.einsum("nij,nkj->nik", right_matrices, left_matrices, optimize=True)
    traces = np.trace(relative, axis1=1, axis2=2)
    cosines = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    angles = np.arccos(cosines)
    skew = np.stack(
        [
            relative[:, 2, 1] - relative[:, 1, 2],
            relative[:, 0, 2] - relative[:, 2, 0],
            relative[:, 1, 0] - relative[:, 0, 1],
        ],
        axis=1,
    )
    sines = np.sin(angles)
    # 2 sin(theta) n = skew, so the rotation vector is theta * skew / (2 sin theta).
    # The ratio tends to 1/2 as theta -> 0, which is exactly the regime of
    # neighbouring EBSD points, so the limit is taken explicitly rather than
    # divided through.
    scale = np.where(
        np.abs(sines) > 1e-12,
        np.divide(angles, 2.0 * sines, out=np.full_like(angles, 0.5), where=np.abs(sines) > 1e-12),
        0.5,
    )
    return freeze_array(np.ascontiguousarray(skew * scale[:, None]))


def _require_regular_grid_with_steps(
    crystal_map: CrystalMap,
) -> tuple[int, int, float, float]:
    rows, cols = crystal_map._require_regular_2d_grid()
    if crystal_map.step_sizes is None or len(crystal_map.step_sizes) != 2:
        raise ValueError(
            "Lattice curvature requires a CrystalMap with 2-D step_sizes: the gradient "
            "is per unit length, so the physical step must be known."
        )
    step_x, step_y = (float(value) for value in crystal_map.step_sizes)
    if step_x <= 0.0 or step_y <= 0.0:
        raise ValueError("CrystalMap.step_sizes must be strictly positive.")
    return rows, cols, step_x, step_y


def lattice_curvature_tensor(
    crystal_map: CrystalMap,
    *,
    step_scale_m: float = _DEFAULT_STEP_SCALE_M,
) -> FloatArray:
    """The measurable part of the lattice curvature tensor, per point.

    Purpose
    -------
    The primitive behind GND density: the spatial gradient of the lattice
    rotation field, :math:`\\kappa_{ij} = \\partial\\omega_i/\\partial x_j`, in
    the specimen frame.

    Method
    ------
    Central differences of the specimen-frame rotation vector along the two map
    axes, falling back to one-sided differences at the map edges so that every
    point receives a value. Neighbouring orientations must belong to the same
    phase; pairs crossing a phase boundary yield ``NaN``, because a
    misorientation between different phases is not a lattice rotation.

    Parameters
    ----------
    crystal_map : CrystalMap
        Must be on a regular 2-D grid and carry ``step_sizes``.
    step_scale_m : float
        Metres per map coordinate unit. The default treats map coordinates as
        micrometres, the EBSD convention. Getting this wrong rescales every
        curvature and every derived density, so it is explicit.

    Returns
    -------
    FloatArray
        ``(rows, cols, 3, 3)`` curvature in radians per metre, read-only. The
        third column — the out-of-plane gradient — is ``NaN`` throughout,
        because a surface map cannot measure it. That is a deliberate marker:
        filling it with zeros would silently assert that the lattice is
        unbent in depth.

    Notes
    -----
    Symmetry is *not* reduced here. Between neighbouring points of a
    continuously bent lattice the true relative rotation is small, and applying a
    disorientation reduction would replace it with a symmetry-equivalent
    representative whose axis is unrelated to the physical curvature. Points
    genuinely straddling a grain boundary therefore produce large, meaningless
    curvature; exclude them with a grain segmentation before interpreting the
    result.
    """

    if not np.isfinite(step_scale_m) or step_scale_m <= 0.0:
        raise ValueError("step_scale_m must be finite and strictly positive.")
    rows, cols, step_x, step_y = _require_regular_grid_with_steps(crystal_map)
    matrices = np.asarray(crystal_map.orientations.as_matrices(), dtype=np.float64)
    grid = matrices.reshape(rows, cols, 3, 3)

    phase_ids = crystal_map.phase_id_array
    phase_grid = (
        None if phase_ids is None else np.asarray(phase_ids, dtype=np.int64).reshape(rows, cols)
    )

    curvature = np.full((rows, cols, 3, 3), np.nan, dtype=np.float64)
    # Column 0 is the gradient along the map x axis (varying column index);
    # column 1 along the map y axis (varying row index). Column 2 stays NaN.
    for axis, (extent, spacing) in enumerate(((cols, step_x), (rows, step_y))):
        if extent < 2:
            continue
        gradient = _directional_rotation_gradient(
            grid,
            phase_grid,
            axis=axis,
            spacing=spacing * step_scale_m,
        )
        curvature[:, :, :, axis] = gradient
    return freeze_array(curvature)


def _directional_rotation_gradient(
    grid: np.ndarray,
    phase_grid: np.ndarray | None,
    *,
    axis: int,
    spacing: float,
) -> FloatArray:
    """Rotation-vector gradient along one map axis, by central differences.

    ``axis=0`` differentiates along columns (the map x direction), ``axis=1``
    along rows (map y). Edges fall back to one-sided differences.
    """

    rows, cols = grid.shape[:2]
    gradient = np.full((rows, cols, 3), np.nan, dtype=np.float64)
    # Work on a view where the differentiated axis is last but one, so a single
    # implementation serves both directions.
    if axis == 0:
        lower = grid[:, :-1], grid[:, 1:]
        phases = None if phase_grid is None else (phase_grid[:, :-1], phase_grid[:, 1:])
        step_shape = (rows, cols - 1)
    else:
        lower = grid[:-1, :], grid[1:, :]
        phases = None if phase_grid is None else (phase_grid[:-1, :], phase_grid[1:, :])
        step_shape = (rows - 1, cols)

    left = lower[0].reshape(-1, 3, 3)
    right = lower[1].reshape(-1, 3, 3)
    deltas = np.asarray(_rotation_vectors(left, right), dtype=np.float64).reshape(*step_shape, 3)
    if phases is not None:
        same_phase = phases[0] == phases[1]
        deltas = np.where(same_phase[..., None], deltas, np.nan)

    # Central difference: the average of the two adjacent one-sided steps,
    # divided by the spacing. Edge points keep their single available step.
    if axis == 0:
        gradient[:, 0, :] = deltas[:, 0, :] / spacing
        gradient[:, -1, :] = deltas[:, -1, :] / spacing
        if cols > 2:
            gradient[:, 1:-1, :] = (deltas[:, :-1, :] + deltas[:, 1:, :]) / (2.0 * spacing)
    else:
        gradient[0, :, :] = deltas[0, :, :] / spacing
        gradient[-1, :, :] = deltas[-1, :, :] / spacing
        if rows > 2:
            gradient[1:-1, :, :] = (deltas[:-1, :, :] + deltas[1:, :, :]) / (2.0 * spacing)
    return freeze_array(gradient)


def nye_dislocation_density_tensor(curvature: np.ndarray) -> FloatArray:
    """Nye's dislocation density tensor from a lattice curvature tensor.

    Purpose
    -------
    Applies :math:`\\alpha_{ij} = \\kappa_{ji} - \\delta_{ij}\\kappa_{kk}`, the
    relation that turns measured lattice bending into the dislocation content
    required to produce it.

    Parameters
    ----------
    curvature : np.ndarray
        ``(..., 3, 3)`` curvature in radians per metre, as returned by
        :func:`lattice_curvature_tensor`.

    Returns
    -------
    FloatArray
        ``(..., 3, 3)`` Nye tensor in reciprocal metres, read-only.

    Notes
    -----
    Because the trace :math:`\\kappa_{kk}` contains the unmeasurable
    :math:`\\kappa_{33}`, the diagonal components of the result are ``NaN`` for
    surface-map curvature. Only the five off-diagonal components carrying
    measurable information are finite, which is the honest outcome: five of nine
    components of the Nye tensor are determinable from a 2-D map.
    """

    array = np.asarray(curvature, dtype=np.float64)
    if array.shape[-2:] != (3, 3):
        raise ValueError("curvature must have trailing shape (3, 3).")
    result = np.ascontiguousarray(np.swapaxes(array, -1, -2)).copy()
    trace = np.trace(array, axis1=-2, axis2=-1)
    # The trace term is applied to the diagonal directly rather than as
    # ``trace * identity``. When the curvature carries NaN for its unmeasurable
    # column the trace is NaN, and ``NaN * 0`` is NaN — so the identity form
    # would poison every off-diagonal component, destroying exactly the
    # measurable information this function exists to produce.
    diagonal = np.arange(3)
    result[..., diagonal, diagonal] -= trace[..., None]
    return freeze_array(result)


def geometrically_necessary_dislocation_density(
    crystal_map: CrystalMap,
    *,
    burgers_vector_nm: float,
    method: GNDMethod = "curvature",
    step_scale_m: float = _DEFAULT_STEP_SCALE_M,
    kam_threshold_deg: float | None = 5.0,
    kam_order: int = 1,
) -> FloatArray:
    """Lower-bound geometrically necessary dislocation density, per point.

    Purpose
    -------
    Convert measured lattice curvature into a dislocation density in
    :math:`\\mathrm{m}^{-2}` — the quantity that connects an orientation map to
    work hardening, stored energy, and recrystallization driving force.

    Methods
    -------
    ``"curvature"`` (default) takes the Nye route: it forms the curvature
    tensor, converts it to the dislocation density tensor, and sums the absolute
    values of the measurable components, dividing by the Burgers vector
    magnitude. This uses the full directional information in the orientation
    gradient.

    ``"kam"`` uses the widely quoted scalar estimate
    :math:`\\rho = 2\\theta / (b u)`, with :math:`\\theta` the kernel average
    misorientation in radians and :math:`u` the step size. It is cruder — it
    discards the direction of the gradient — but it is what much of the EBSD
    literature reports, so it is provided for comparability.

    Parameters
    ----------
    crystal_map : CrystalMap
        Regular 2-D grid with ``step_sizes``.
    burgers_vector_nm : float
        Burgers vector magnitude in nanometres: ``0.2556`` for copper,
        ``0.2483`` for alpha-iron, ``0.2863`` for aluminium. Strictly positive.
    method : str
        ``"curvature"`` or ``"kam"``.
    step_scale_m : float
        Metres per map coordinate unit; the default treats them as micrometres.
    kam_threshold_deg : float, optional
        For the ``"kam"`` method, the misorientation above which a neighbour
        pair is excluded as a grain boundary rather than a lattice gradient.
        Without it, boundary pixels report the boundary misorientation and the
        density there is meaningless. Ignored by the curvature method.
    kam_order : int
        Neighbour shell for the ``"kam"`` method.

    Returns
    -------
    FloatArray
        ``(rows, cols)`` density in :math:`\\mathrm{m}^{-2}`, read-only.
        ``NaN`` where the underlying gradient could not be measured — across a
        phase boundary, for example.

    Notes
    -----
    **This is a lower bound, and it is resolution dependent.** Dislocation
    content that produces no in-plane curvature is invisible to a surface map,
    and statistically stored dislocations are invisible by definition, since
    their Burgers vectors cancel. The value also rises as the step size falls,
    because a finer step resolves sharper gradients; report the step size with
    any density, and compare maps only at equal step size.

    Grain boundaries are *not* excluded automatically by the curvature method.
    Across a boundary the relative rotation is not a lattice gradient, and the
    reported density there is an artefact. Mask by a grain segmentation, or use
    the ``"kam"`` method with a threshold, before interpreting a map.
    """

    if not np.isfinite(burgers_vector_nm) or burgers_vector_nm <= 0.0:
        raise ValueError("burgers_vector_nm must be finite and strictly positive.")
    if method not in {"curvature", "kam"}:
        raise ValueError("method must be either 'curvature' or 'kam'.")
    burgers_m = float(burgers_vector_nm) * _NANOMETRE_M

    if method == "kam":
        rows, cols, step_x, step_y = _require_regular_grid_with_steps(crystal_map)
        kam_deg = np.asarray(
            crystal_map.kernel_average_misorientation_deg(
                order=kam_order,
                threshold_deg=kam_threshold_deg,
            ),
            dtype=np.float64,
        ).reshape(rows, cols)
        # The conventional estimate uses one step length; a rectangular grid has
        # two, so the mean is used and the choice is stated rather than hidden.
        step_m = 0.5 * (step_x + step_y) * step_scale_m
        density = 2.0 * np.deg2rad(kam_deg) / (burgers_m * step_m)
        return freeze_array(np.ascontiguousarray(density))

    curvature = np.asarray(lattice_curvature_tensor(crystal_map, step_scale_m=step_scale_m))
    nye = np.asarray(nye_dislocation_density_tensor(curvature))
    # Only components a surface map determines contribute. Since alpha_ij =
    # kappa_ji off the diagonal, and kappa's third column is unmeasured, the
    # determinable off-diagonal components are exactly those with i in {0, 1}:
    # alpha_01, alpha_02, alpha_10, alpha_12. The diagonal is excluded because
    # the trace carries the unmeasurable kappa_22; only the *difference*
    # alpha_00 - alpha_11 = kappa_00 - kappa_11 is determinable, and that fifth
    # piece of information has no place in an L1 sum of components.
    contributions = np.abs(nye[..., _MEASURABLE_NYE_COMPONENTS])
    density = np.sum(contributions, axis=-1) / burgers_m
    return freeze_array(np.ascontiguousarray(density))


__all__ = [
    "GNDMethod",
    "geometrically_necessary_dislocation_density",
    "lattice_curvature_tensor",
    "nye_dislocation_density_tensor",
]
