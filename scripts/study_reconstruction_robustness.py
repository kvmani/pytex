"""Robustness study for map-scale parent-grain reconstruction.

Characterizes where `pytex.experimental.reconstruct_parent_grains` succeeds and
where it degrades, as the evidence required before the surface can leave
`experimental`. Three axes are swept against planted ground truth:

* **orientation noise** on the child grains (EBSD indexing scatter),
* **edge tolerance**, the one free parameter of the edge test,
* **grain count per parent**, which controls how much averaging the parent
  estimate gets and how often a singleton stays ambiguous.

Reported per cell, over several microstructure seeds:

``partition_exact``
    Fraction of trials whose grouping equals the planted partition exactly.
``parent_error_deg``
    Symmetry-reduced angle between each recovered parent and its planted
    parent, aggregated over trials.
``false_link_rate`` / ``missed_link_rate``
    Cross-parent edges wrongly linked, and same-parent edges wrongly rejected.

Trials whose planted ground truth is genuinely ambiguous at the tolerance under
test — a cross-parent boundary that really does land on the same-parent
fingerprint — are counted separately as ``ambiguous`` and excluded from the
accuracy statistics, because no algorithm can be expected to resolve them from
orientations alone.

Usage::

    python scripts/study_reconstruction_robustness.py            # full sweep
    python scripts/study_reconstruction_robustness.py --quick    # CI smoke

Results print as a table and are written to a local, git-ignored JSON file;
``benchmarks/`` stays reserved for schema-validated manifests.
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    OrientationRelationship,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    boundary_fingerprint_distances_deg,
    intervariant_boundary_fingerprint,
)
from pytex.core.transformation import _symmetry_reduced_angle_between_deg
from pytex.experimental import reconstruct_parent_grains

_IDENTITY = np.eye(3, dtype=np.float64)[None, :, :]


def _phases() -> tuple[Phase, Phase, ReferenceFrame]:
    parent_frame = ReferenceFrame(
        name="study_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    child_frame = ReferenceFrame(
        name="study_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="study_specimen",
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


def _plant(
    relationship: OrientationRelationship,
    specimen: ReferenceFrame,
    *,
    parent_count: int,
    per_parent: int,
    noise_deg: float,
    seed: int,
) -> tuple[OrientationSet, np.ndarray, np.ndarray, list[np.ndarray]]:
    """Random parents, distinct variants each, chain edges plus parent contacts."""

    child_phase = relationship.child_phase
    variants = relationship.generate_variants()
    rng = np.random.default_rng(seed)
    quaternions: list[np.ndarray] = []
    planted: list[int] = []
    parent_matrices: list[np.ndarray] = []
    for _ in range(parent_count):
        quaternion = rng.normal(size=4)
        quaternion /= np.linalg.norm(quaternion)
        parent_rotation = Rotation(quaternion=quaternion)
        parent_matrices.append(parent_rotation.as_matrix())
        picks = rng.choice(len(variants), size=per_parent, replace=False)
        for pick in picks:
            rotation = parent_rotation.compose(
                variants[int(pick)].parent_to_child_rotation.inverse()
            )
            if noise_deg > 0.0:
                axis = rng.normal(size=3)
                axis /= np.linalg.norm(axis)
                rotation = Rotation.from_axis_angle(
                    axis, np.deg2rad(rng.normal(0.0, noise_deg))
                ).compose(rotation)
            quaternions.append(rotation.quaternion)
            planted.append(len(parent_matrices) - 1)
    children = OrientationSet(
        quaternions=np.stack(quaternions, axis=0),
        crystal_frame=child_phase.crystal_frame,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    labels = np.asarray(planted, dtype=np.int64)
    edges: list[tuple[int, int]] = []
    for index in range(parent_count):
        members = np.flatnonzero(labels == index)
        edges.extend(
            (int(members[k]), int(members[k + 1])) for k in range(members.size - 1)
        )
    for index in range(parent_count - 1):
        edges.append(
            (
                int(np.flatnonzero(labels == index)[-1]),
                int(np.flatnonzero(labels == index + 1)[0]),
            )
        )
    return children, np.asarray(edges, dtype=np.int64), labels, parent_matrices


def _partitions_equal(left: np.ndarray, right: np.ndarray) -> bool:
    return bool(
        np.array_equal(
            left[:, None] == left[None, :], right[:, None] == right[None, :]
        )
    )


def _cell(
    *,
    noise_deg: float,
    tolerance_deg: float,
    per_parent: int,
    parent_count: int,
    seeds: range,
) -> dict[str, object]:
    parent_phase, child_phase, specimen = _phases()
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    fingerprint = intervariant_boundary_fingerprint(relationship)
    parent_operators = parent_phase.symmetry.operators

    exact = 0
    judged = 0
    ambiguous = 0
    parent_errors: list[float] = []
    false_links = 0
    cross_total = 0
    missed_links = 0
    same_total = 0

    for seed in seeds:
        children, edges, planted, parent_matrices = _plant(
            relationship,
            specimen,
            parent_count=parent_count,
            per_parent=per_parent,
            noise_deg=noise_deg,
            seed=seed,
        )
        matrices = children.as_matrices()
        relative = np.einsum(
            "eji,ejk->eik", matrices[edges[:, 0]], matrices[edges[:, 1]], optimize=True
        )
        cross = planted[edges[:, 0]] != planted[edges[:, 1]]
        distances = boundary_fingerprint_distances_deg(relative, fingerprint)
        # A cross-parent boundary that genuinely lands on the fingerprint is a
        # physical ambiguity, not an algorithmic failure: exclude the trial.
        if float(distances[cross].min()) <= tolerance_deg:
            ambiguous += 1
            continue

        result = reconstruct_parent_grains(
            children, edges, relationship, tolerance_deg=tolerance_deg
        )
        judged += 1
        if _partitions_equal(result.parent_labels, planted):
            exact += 1
        linked = distances <= tolerance_deg
        false_links += int(np.count_nonzero(linked & cross))
        cross_total += int(np.count_nonzero(cross))
        missed_links += int(np.count_nonzero(~linked & ~cross))
        same_total += int(np.count_nonzero(~cross))

        recovered = result.parent_orientations.as_matrices()
        for index, truth in enumerate(parent_matrices):
            members = np.flatnonzero(planted == index)
            if members.size == 0:
                continue
            cluster = int(result.parent_labels[members[0]])
            parent_errors.append(
                _symmetry_reduced_angle_between_deg(
                    recovered[cluster],
                    truth,
                    child_operators=_IDENTITY,
                    parent_operators=parent_operators,
                )
            )

    errors = np.asarray(parent_errors, dtype=np.float64)
    return {
        "noise_deg": noise_deg,
        "tolerance_deg": tolerance_deg,
        "children_per_parent": per_parent,
        "parent_count": parent_count,
        "trials_judged": judged,
        "trials_ambiguous": ambiguous,
        "partition_exact": (exact / judged) if judged else float("nan"),
        "parent_error_mean_deg": float(errors.mean()) if errors.size else float("nan"),
        "parent_error_max_deg": float(errors.max()) if errors.size else float("nan"),
        "false_link_rate": (false_links / cross_total) if cross_total else 0.0,
        "missed_link_rate": (missed_links / same_total) if same_total else 0.0,
    }


def run(quick: bool) -> dict[str, object]:
    seeds = range(6) if quick else range(25)
    noises = (0.0, 0.5, 2.0) if quick else (0.0, 0.25, 0.5, 1.0, 2.0)
    tolerances = (3.0,) if quick else (1.0, 2.0, 3.0, 5.0)
    sizes = (5,) if quick else (2, 5, 12)

    cells: list[dict[str, object]] = []
    for noise in noises:
        for tolerance in tolerances:
            for per_parent in sizes:
                # A tolerance below the noise cannot admit the noisy same-parent
                # boundaries; skip that degenerate corner rather than report it
                # as an algorithmic failure.
                if tolerance < 2.0 * noise:
                    continue
                cells.append(
                    _cell(
                        noise_deg=noise,
                        tolerance_deg=tolerance,
                        per_parent=per_parent,
                        parent_count=6,
                        seeds=seeds,
                    )
                )
    return {
        "study": "parent_grain_reconstruction_robustness",
        "relationship": "kurdjumov_sachs",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "seeds": len(seeds),
        "cells": cells,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="CI smoke sweep")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reconstruction_robustness.json"),
        help="JSON output path (full runs only)",
    )
    arguments = parser.parse_args()
    payload = run(arguments.quick)

    header = (
        f"{'noise':>6} {'tol':>5} {'n/par':>6} {'exact':>7} {'parent err':>18} "
        f"{'false':>7} {'missed':>7} {'ambig':>6}"
    )
    print(header)
    print("-" * len(header))
    for cell in payload["cells"]:  # type: ignore[union-attr]
        print(
            f"{cell['noise_deg']:6.2f} {cell['tolerance_deg']:5.1f} "
            f"{cell['children_per_parent']:6d} {cell['partition_exact']:6.0%} "
            f"{cell['parent_error_mean_deg']:8.3f} /{cell['parent_error_max_deg']:8.3f} "
            f"{cell['false_link_rate']:6.1%} {cell['missed_link_rate']:6.1%} "
            f"{cell['trials_ambiguous']:6d}"
        )
    if not arguments.quick:
        arguments.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
