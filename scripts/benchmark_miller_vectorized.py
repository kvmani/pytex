from __future__ import annotations

import time

import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlaneSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)


def _make_phase() -> Phase:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    return Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def _random_planes(count: int, *, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    indices = rng.integers(-3, 4, size=(count, 3), dtype=np.int64)
    zero_rows = ~np.any(indices != 0, axis=1)
    indices[zero_rows] = np.array([1, 0, 0], dtype=np.int64)
    return indices


def main() -> None:
    phase = _make_phase()
    planes = MillerPlaneSet.from_hkl(_random_planes(100_000), phase=phase)

    start = time.perf_counter()
    normals = planes.normals_cartesian()
    d_spacings = planes.d_spacings_angstrom()
    normals_elapsed = time.perf_counter() - start

    start = time.perf_counter()
    equivalent_indices, mask = planes.symmetry_equivalent_indices()
    symmetry_elapsed = time.perf_counter() - start

    print("PyTex Miller benchmark")
    print(f"planes: {planes.indices.shape[0]}")
    print(
        "normals+d-spacings: "
        f"{normals_elapsed:.3f} s "
        f"(normals shape={normals.shape}, d-spacings shape={d_spacings.shape})"
    )
    print(
        "symmetry expansion: "
        f"{symmetry_elapsed:.3f} s "
        f"(equivalents shape={equivalent_indices.shape}, mask shape={mask.shape})"
    )


if __name__ == "__main__":
    main()
