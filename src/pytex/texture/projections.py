from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import normalize_vectors


def fold_upper_hemisphere(directions: ArrayLike, *, antipodal: bool = True) -> np.ndarray:
    folded = np.array(normalize_vectors(directions), copy=True)
    if antipodal:
        mask = folded[:, 2] < 0.0
        folded[mask] *= -1.0
    folded = np.ascontiguousarray(folded)
    folded.setflags(write=False)
    return folded


def project_directions(
    directions: ArrayLike,
    *,
    method: str = "equal_area",
    antipodal: bool = True,
) -> np.ndarray:
    vectors = fold_upper_hemisphere(directions, antipodal=antipodal)
    denominator = np.clip(1.0 + vectors[:, 2], 1e-12, None)
    if method == "equal_area":
        scale = np.sqrt(2.0 / denominator)
    elif method == "stereographic":
        scale = 1.0 / denominator
    else:
        raise ValueError("Projection method must be either 'equal_area' or 'stereographic'.")
    projected = np.column_stack([vectors[:, 0] * scale, vectors[:, 1] * scale])
    projected = np.ascontiguousarray(projected)
    projected.setflags(write=False)
    return projected
