"""Full-constraint (Taylor) polycrystal plasticity factors.

The Taylor factor ``M`` relates the sum of crystallographic slip needed to
accommodate an imposed macroscopic strain increment to that strain's von Mises
equivalent. Under the Taylor full-constraint hypothesis every grain undergoes
the same strain as the aggregate, and the active slip combination is the one
minimising the total slip:

``minimise  sum_s gamma_s   subject to   sum_s gamma_s M_s = eps,   gamma_s >= 0``

where ``M_s`` is the symmetric Schmid tensor of signed slip system ``s`` in the
sample frame and ``eps`` is the (deviatoric) imposed strain. This is a small
linear program solved per orientation with SciPy's HiGHS solver. All problem
setup -- Schmid tensors, sample-frame rotation, constraint matrices -- is
vectorised over the whole orientation population; only the LP solve itself is
per orientation, as required by the model.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import linprog

from pytex.core.orientation import Orientation, OrientationSet
from pytex.properties.slip import SlipSystemFamily

# Independent components of a symmetric traceless tensor. Matching these five
# also matches the (3, 3) component, since both operands are traceless.
_COMPONENT_INDICES: tuple[tuple[int, int], ...] = ((0, 0), (1, 1), (0, 1), (0, 2), (1, 2))


def uniaxial_strain_tensor(axis: ArrayLike) -> np.ndarray:
    """Deviatoric strain tensor for uniaxial tension along ``axis``.

    Normalised to unit von Mises equivalent strain, so the Taylor factor equals
    the minimal total slip directly.
    """

    unit = np.asarray(axis, dtype=np.float64)
    norm = float(np.linalg.norm(unit))
    if norm == 0.0:
        raise ValueError("axis must be a non-zero vector.")
    unit = unit / norm
    return 1.5 * np.outer(unit, unit) - 0.5 * np.eye(3, dtype=np.float64)


def _equivalent_strain(strain: np.ndarray) -> float:
    return float(np.sqrt(2.0 / 3.0 * np.tensordot(strain, strain)))


def _sample_frame_schmid_tensors(
    family: SlipSystemFamily,
    matrices: np.ndarray,
) -> np.ndarray:
    """Symmetric Schmid tensors of every slip system in the sample frame.

    Returns an array of shape ``(n_orientations, n_systems, 3, 3)``.
    """

    normals = family.plane_normals
    directions = family.slip_directions
    crystal_schmid = 0.5 * (
        directions[:, :, None] * normals[:, None, :]
        + normals[:, :, None] * directions[:, None, :]
    )
    # M_sample = R @ m_crystal @ R^T for every (orientation, system) pair.
    return np.asarray(
        np.einsum(
            "oip,spq,ojq->osij",
            matrices,
            crystal_schmid,
            matrices,
            optimize=True,
        ),
        dtype=np.float64,
    )


def taylor_factors(
    family: SlipSystemFamily,
    orientations: Orientation | OrientationSet,
    *,
    strain_tensor: ArrayLike | None = None,
    tension_axis: ArrayLike | None = None,
) -> np.ndarray | float:
    """Full-constraint Taylor factor for each orientation.

    Provide either an explicit deviatoric ``strain_tensor`` (3x3 symmetric,
    trace ~0) or a ``tension_axis`` for uniaxial tension (default ``z``). Returns
    an array of Taylor factors (or a float for a single `Orientation`); a value
    is ``inf`` when the imposed strain cannot be accommodated by the family.
    """

    if strain_tensor is not None and tension_axis is not None:
        raise ValueError("Provide only one of strain_tensor or tension_axis.")
    if strain_tensor is not None:
        strain = np.asarray(strain_tensor, dtype=np.float64)
        if strain.shape != (3, 3):
            raise ValueError("strain_tensor must have shape (3, 3).")
        if not np.allclose(strain, strain.T, atol=1e-9):
            raise ValueError("strain_tensor must be symmetric.")
        if abs(float(np.trace(strain))) > 1e-6:
            raise ValueError("strain_tensor must be deviatoric (trace ~ 0).")
    else:
        strain = uniaxial_strain_tensor(tension_axis if tension_axis is not None else (0, 0, 1))

    equivalent = _equivalent_strain(strain)
    if equivalent == 0.0:
        raise ValueError("Imposed strain must be non-zero.")

    if isinstance(orientations, Orientation):
        scalar = True
        orientation_set = OrientationSet.from_orientations([orientations])
    else:
        scalar = False
        orientation_set = orientations
    matrices = orientation_set.as_matrices()
    schmid = _sample_frame_schmid_tensors(family, matrices)
    rows = np.array([schmid[:, :, i, j] for (i, j) in _COMPONENT_INDICES])  # (5, n_orient, n_sys)
    target = np.array([strain[i, j] for (i, j) in _COMPONENT_INDICES], dtype=np.float64)
    system_count = schmid.shape[1]
    objective = np.ones(2 * system_count, dtype=np.float64)

    factors = np.empty(matrices.shape[0], dtype=np.float64)
    for index in range(matrices.shape[0]):
        # Signed slip: append the negated Schmid tensors so gamma >= 0 covers
        # both slip senses. Columns are [+systems | -systems].
        a_eq = np.concatenate([rows[:, index, :], -rows[:, index, :]], axis=1)
        result = linprog(objective, A_eq=a_eq, b_eq=target, bounds=(0.0, None), method="highs")
        factors[index] = result.fun / equivalent if result.success else np.inf

    factors = np.ascontiguousarray(factors)
    factors.setflags(write=False)
    return float(factors[0]) if scalar else factors


__all__ = [
    "taylor_factors",
    "uniaxial_strain_tensor",
]
