"""Runnable performance benchmark for the transformation / OR-analysis hot paths.

Development-guide finding 21: the repository keeps performance evidence
runnable, with pinned case definitions (sizes and seeds) so successive runs
are comparable. Results are printed as a table and written as JSON to a local,
git-ignored file (``performance_results.json``); ``benchmarks/`` stays
reserved for schema-validated manifests.

Usage::

    python scripts/benchmark_transformation_performance.py            # full sizes
    python scripts/benchmark_transformation_performance.py --quick    # CI smoke sizes

Timing assertions are deliberately NOT enforced (CI runners are too noisy);
the evidence is the recorded JSON trend, not a pass/fail gate.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    OrientationRelationship,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    fit_orientation_relationship,
    intervariant_misorientations,
    or_deviation,
)
from pytex.experimental import reconstruct_parent_grains

_SEED = 20260717


def _phases() -> tuple[Phase, Phase, ReferenceFrame]:
    parent_frame = ReferenceFrame(
        name="bench_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    child_frame = ReferenceFrame(
        name="bench_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="bench_specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    parent = Phase(
        "austenite",
        lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
        crystal_frame=parent_frame,
    )
    child = Phase(
        "ferrite",
        lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
        crystal_frame=child_frame,
    )
    return parent, child, specimen


def _paired_sets(
    count: int,
) -> tuple[OrientationSet, OrientationSet, OrientationRelationship]:
    parent_phase, child_phase, specimen = _phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variants = ks.generate_variants()
    rng = np.random.default_rng(_SEED)
    eulers = rng.uniform(0.0, 80.0, size=(count, 3))
    parents = OrientationSet.from_orientations(
        [
            Orientation.from_euler(
                *euler, specimen_frame=specimen, symmetry=parent_phase.symmetry,
                phase=parent_phase,
            )
            for euler in eulers
        ]
    )
    picks = rng.integers(0, len(variants), size=count)
    quaternions = np.stack(
        [
            parents[index]
            .rotation.compose(variants[int(picks[index])].parent_to_child_rotation.inverse())
            .quaternion
            for index in range(count)
        ],
        axis=0,
    )
    children = OrientationSet(
        quaternions=quaternions,
        crystal_frame=child_phase.crystal_frame,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    return parents, children, ks


def _time(function: Callable[[], object], *, repeats: int) -> float:
    best = float("inf")
    for _ in range(repeats):
        start = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - start)
    return best


def run(quick: bool) -> dict[str, object]:
    pair_count = 200 if quick else 5_000
    grain_count = 30 if quick else 400
    repeats = 1 if quick else 3
    parents, children, ks = _paired_sets(pair_count)
    grain_parents, grain_children, _ = _paired_sets(grain_count)
    adjacency = np.column_stack(
        [np.arange(grain_count - 1), np.arange(1, grain_count)]
    )
    cases: dict[str, Callable[[], object]] = {
        "intervariant_misorientations_ks_276_pairs": lambda: intervariant_misorientations(ks),
        f"or_deviation_{pair_count}_pairs": lambda: or_deviation(parents, children, ks),
        f"fit_orientation_relationship_{pair_count}_pairs": lambda: (
            fit_orientation_relationship(parents, children, ks)
        ),
        f"reconstruct_parent_grains_{grain_count}_grains": lambda: (
            reconstruct_parent_grains(grain_children, adjacency, ks, tolerance_deg=2.0)
        ),
    }
    results = {
        name: {"best_seconds": round(_time(function, repeats=repeats), 4)}
        for name, function in cases.items()
    }
    return {
        "schema_id": "pytex.benchmarks.transformation_performance",
        "schema_version": "1",
        "seed": _SEED,
        "quick": quick,
        "python": platform.python_version(),
        "machine": platform.machine(),
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="CI smoke sizes")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("performance_results.json"),
        help="JSON output path (full runs only)",
    )
    arguments = parser.parse_args()
    payload = run(arguments.quick)
    for name, entry in payload["results"].items():  # type: ignore[union-attr]
        print(f"{name:55s} {entry['best_seconds']:>10.4f} s")
    if not arguments.quick:
        arguments.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
