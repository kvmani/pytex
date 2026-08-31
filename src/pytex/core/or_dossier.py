"""The orientation-relationship dossier: one declaration, one artifact.

An orientation relationship is reported in the literature as a scatter of
quantities across a paper — a cell table here, an axis/angle there, a variant
count in a caption, a parallelism in the abstract. Reproducing such a report
means recomputing all of it from a sentence. This module turns the declaration
into a single typed object that carries every number, states the convention each
was computed under, and writes itself to disk as a bundle another person can
check.

The blocks follow the flagship program's F17:

1. **Lattice** — both phases' cells, metric tensors, structure matrices,
   volumes and point groups.
2. **Transformation** — the parent-to-child rotation, the direction and plane
   correspondence matrices, and the lattice-correspondence deformation with its
   principal strains and volume change.
3. **Misorientation** — the symmetry-reduced axis/angle representative, the
   variant count and packet grouping, and the intervariant spectrum.
4. **Parallelism** — the defining parallelisms of the chosen variant, plus any
   near-parallelisms discovered among nominated families, in publication
   notation.
5. **Figures** — written by :meth:`ORDossier.export`, not held in memory.

**The rule this module is built on:** it *calls* the existing functions and
never reimplements them. A dossier number that disagreed with the function it
came from would be the exact class of defect this repository exists to prevent,
so every value below is a call into `pytex.core.transformation`,
`pytex.core.lattice` or `pytex.plotting`, and the tests assert the agreement
rather than assuming it.

**What is not here.** The interface block (F16) is not implemented, and the
dossier says so in prose and carries ``null`` for it rather than omitting the
key — a reader must be able to tell "not analysed" from "analysed and empty".
The figure bundle covers the OR stereogram and the variant contact sheet, which
need only the relationship; the variant pole figure and the per-variant SAED
patterns need a measured parent orientation and a diffraction setup, neither of
which a relationship alone supplies.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex, Phase
from pytex.core.miller import canonicalize_sign
from pytex.core.notation import (
    format_direction_indices,
    format_plane_family_indices,
    format_plane_indices,
)
from pytex.core.transformation import (
    OrientationRelationship,
    TransformationVariant,
    find_parallel_directions,
    find_parallel_planes,
    intervariant_misorientation_angles_deg,
    variant_close_packed_groups,
)

__all__ = [
    "OR_DOSSIER_SCHEMA_ID",
    "OR_DOSSIER_SCHEMA_VERSION",
    "ORDossier",
    "ORDossierLatticeBlock",
    "ORDossierMisorientationBlock",
    "ORDossierParallelism",
    "ORDossierParallelismBlock",
    "ORDossierTransformationBlock",
    "or_dossier",
    "or_dossier_schema_path",
]

#: Stable identifier of the dossier JSON contract.
OR_DOSSIER_SCHEMA_ID = "pytex.or_dossier"

#: Version of that contract. Bump on any change to the key set.
OR_DOSSIER_SCHEMA_VERSION = "1.0.0"

#: Default angular tolerance for the *discovered* parallelisms.
#:
#: Read it as the rationalization tolerance it is: the exact child image of a
#: parent plane is parallel to it by construction, so what
#: :func:`~pytex.core.transformation.find_parallel_planes` measures is how far
#: the nearest low-index child index sits from that exact image. Half a degree
#: keeps the pairs a low-index child object really does realize.
DEFAULT_DISCOVERY_TOLERANCE_DEG = 0.5

_CITATION_ITA = (
    "International Tables for Crystallography, Volume A (cell conventions and "
    "point groups)."
)
_CITATION_MORITO = (
    "Morito, Tanaka, Konishi, Furuhara & Maki, Acta Materialia 51 (2003) 1789 "
    "(packet structure and the intervariant table)."
)


def or_dossier_schema_path() -> Path:
    """Path to the JSON schema the dossier's ``to_json()`` is written against."""

    return Path(__file__).resolve().parents[3] / "schemas" / "or_dossier.schema.json"


def _matrix_rows(values: Any) -> list[list[float]]:
    return [[float(entry) for entry in row] for row in np.asarray(values, dtype=np.float64)]


def _index_triple(values: Any) -> tuple[int, int, int]:
    rounded = [round(float(value)) for value in np.asarray(values).reshape(-1)]
    return (rounded[0], rounded[1], rounded[2])


def _plane_triple(values: Any) -> tuple[int, int, int]:
    """A plane's indices under the repository's one sign rule.

    A plane has no sign, and which of ``(1 1 -1)`` and ``(-1 -1 1)`` a symmetry
    image comes back as is an artefact of the arithmetic. Directions keep their
    sign, because there the two spellings are opposite directions.
    """

    canonical = canonicalize_sign(_index_triple(values))[0]
    return (int(canonical[0]), int(canonical[1]), int(canonical[2]))


# --------------------------------------------------------------------------- #
# Block 1 — the lattices
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORDossierLatticeBlock:
    """One phase's cell, metric and symmetry, as the dossier reports them.

    Every quantity is read from the phase rather than recomputed: the cell
    parameters from its ``Lattice``, the structure matrix from
    ``direct_basis()``, the metric tensors from the lattice's own methods, and
    the volume from ``Lattice.volume_angstrom3``. The Cartesian setting of the
    structure matrix is the single one PyTex defines — ``a`` along ``x``, ``b``
    in the ``x-y`` plane — and the dossier states it rather than leaving a
    reader to infer it from the numbers.
    """

    phase_name: str
    point_group: str
    cell_lengths_angstrom: tuple[float, float, float]
    cell_angles_deg: tuple[float, float, float]
    volume_angstrom3: float
    direct_basis: np.ndarray
    reciprocal_basis: np.ndarray
    metric_tensor: np.ndarray
    reciprocal_metric_tensor: np.ndarray

    def to_json(self) -> dict[str, Any]:
        """The block as JSON, matrices as row lists."""

        return {
            "phase_name": self.phase_name,
            "point_group": self.point_group,
            "cell_lengths_angstrom": [float(value) for value in self.cell_lengths_angstrom],
            "cell_angles_deg": [float(value) for value in self.cell_angles_deg],
            "volume_angstrom3": float(self.volume_angstrom3),
            "direct_basis": _matrix_rows(self.direct_basis),
            "reciprocal_basis": _matrix_rows(self.reciprocal_basis),
            "metric_tensor": _matrix_rows(self.metric_tensor),
            "reciprocal_metric_tensor": _matrix_rows(self.reciprocal_metric_tensor),
        }

    def describe(self) -> str:
        """Prose summary: cell, volume and point group, with the frame stated."""

        lengths = ", ".join(f"{value:.4f}" for value in self.cell_lengths_angstrom)
        angles = ", ".join(f"{value:.3f}" for value in self.cell_angles_deg)
        return (
            f"{self.phase_name}: point group {self.point_group}, cell "
            f"a, b, c = {lengths} angstrom and alpha, beta, gamma = {angles} degrees, "
            f"volume {self.volume_angstrom3:.4f} cubic angstrom. The structure matrix "
            "has the crystal axes as its columns in the PyTex Cartesian setting "
            "(a along x, b in the x-y plane); the reciprocal basis is normalized "
            "so that a*_i . a_j = delta_ij."
        )


def _lattice_block(phase: Phase) -> ORDossierLatticeBlock:
    lattice = phase.lattice
    point_group = phase.symmetry.point_group if phase.symmetry is not None else "1"
    return ORDossierLatticeBlock(
        phase_name=phase.name,
        point_group=str(point_group),
        cell_lengths_angstrom=(float(lattice.a), float(lattice.b), float(lattice.c)),
        cell_angles_deg=(
            float(lattice.alpha_deg),
            float(lattice.beta_deg),
            float(lattice.gamma_deg),
        ),
        volume_angstrom3=lattice.volume_angstrom3(),
        direct_basis=np.asarray(lattice.direct_basis().matrix, dtype=np.float64),
        reciprocal_basis=np.asarray(lattice.reciprocal_basis().matrix, dtype=np.float64),
        metric_tensor=np.asarray(lattice.metric_tensor(), dtype=np.float64),
        reciprocal_metric_tensor=np.asarray(
            lattice.reciprocal_metric_tensor(), dtype=np.float64
        ),
    )


# --------------------------------------------------------------------------- #
# Block 2 — the transformation
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORDossierTransformationBlock:
    """The rotation, the two correspondence matrices, and the strain.

    ``correspondence_direct`` carries parent ``[uvw]`` to child ``[uvw]`` and
    ``correspondence_reciprocal`` carries parent ``(hkl)`` to child ``(hkl)``;
    they are inverse transposes of one another, which is what preserves the zone
    law across the mapping. The strain quantities come from
    ``OrientationRelationship.deformation_gradient``, whose report is carried
    whole rather than unpacked, so the dossier and that report cannot disagree.
    """

    rotation_matrix: np.ndarray
    rotation_axis: np.ndarray
    rotation_angle_deg: float
    correspondence_direct: np.ndarray
    correspondence_reciprocal: np.ndarray
    principal_strains_percent: tuple[float, float, float]
    volume_change_percent: float
    polar_rotation_deg: float
    correspondence_denominator: int
    deformation_description: str

    def to_json(self) -> dict[str, Any]:
        """The block as JSON."""

        return {
            "rotation_matrix": _matrix_rows(self.rotation_matrix),
            "rotation_axis": [float(value) for value in self.rotation_axis],
            "rotation_angle_deg": float(self.rotation_angle_deg),
            "correspondence_direct": _matrix_rows(self.correspondence_direct),
            "correspondence_reciprocal": _matrix_rows(self.correspondence_reciprocal),
            "principal_strains_percent": [
                float(value) for value in self.principal_strains_percent
            ],
            "volume_change_percent": float(self.volume_change_percent),
            "polar_rotation_deg": float(self.polar_rotation_deg),
            "correspondence_denominator": int(self.correspondence_denominator),
        }

    def describe(self) -> str:
        """Prose summary: the rotation, then the strain report verbatim."""

        axis = ", ".join(f"{value:.4f}" for value in self.rotation_axis)
        return (
            f"Parent-to-child rotation: {self.rotation_angle_deg:.3f} degrees about "
            f"({axis}) in the parent crystal frame — the rotation as declared, not the "
            "symmetry-reduced representative, which the misorientation block reports. "
            "The direction correspondence carries parent [uvw] to child [uvw] and the "
            "plane correspondence carries parent (hkl) to child (hkl); they are inverse "
            "transposes, which is what preserves the zone law. "
            + self.deformation_description
        )


def _transformation_block(
    relationship: OrientationRelationship, variant: TransformationVariant | None
) -> ORDossierTransformationBlock:
    rotation = (
        relationship.parent_to_child_rotation
        if variant is None
        else variant.parent_to_child_rotation
    )
    deformation = relationship.deformation_gradient(variant=variant)
    strains = (np.asarray(deformation.principal_stretches, dtype=np.float64) - 1.0) * 100.0
    return ORDossierTransformationBlock(
        rotation_matrix=np.asarray(rotation.as_matrix(), dtype=np.float64),
        rotation_axis=np.asarray(rotation.axis, dtype=np.float64),
        rotation_angle_deg=float(rotation.angle_deg),
        correspondence_direct=relationship.correspondence_direct(variant=variant),
        correspondence_reciprocal=relationship.correspondence_reciprocal(variant=variant),
        principal_strains_percent=(float(strains[0]), float(strains[1]), float(strains[2])),
        volume_change_percent=float((deformation.volume_ratio - 1.0) * 100.0),
        polar_rotation_deg=float(deformation.polar_rotation_deg),
        correspondence_denominator=int(deformation.correspondence_denominator),
        deformation_description=deformation.describe(),
    )


# --------------------------------------------------------------------------- #
# Block 3 — the misorientation and the variant family
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORDossierMisorientationBlock:
    """The relationship as a measured quantity, and the family it generates.

    ``angle_deg`` and ``axis`` are the symmetry-reduced *disorientation*
    representative — the form an EBSD boundary measurement is reported in, and
    the form the literature tabulates. ``packet_labels`` groups the variants by
    the parent plane each carries into exact parallelism, which is the packet of
    lath martensite; ``intervariant_angles_deg`` is the distinct spectrum those
    variants can make with each other.
    """

    angle_deg: float
    axis: np.ndarray
    variant_count: int
    packet_plane: tuple[int, int, int] | None
    packet_labels: np.ndarray
    packet_count: int
    intervariant_angles_deg: np.ndarray

    def to_json(self) -> dict[str, Any]:
        """The block as JSON."""

        return {
            "angle_deg": float(self.angle_deg),
            "axis": [float(value) for value in self.axis],
            "variant_count": int(self.variant_count),
            "packet_plane": (
                None if self.packet_plane is None else [int(v) for v in self.packet_plane]
            ),
            "packet_labels": [int(value) for value in self.packet_labels],
            "packet_count": int(self.packet_count),
            "intervariant_angles_deg": [float(value) for value in self.intervariant_angles_deg],
        }

    def describe(self) -> str:
        """Prose summary: disorientation, variant count, packets, spectrum."""

        axis = ", ".join(f"{value:.4f}" for value in self.axis)
        spectrum = ", ".join(f"{value:.2f}" for value in self.intervariant_angles_deg)
        packet_text = (
            "no packet grouping was requested"
            if self.packet_plane is None
            else (
                f"grouped on the parent "
                f"{format_plane_family_indices(self.packet_plane, style='plain')} family, the "
                f"{self.variant_count} variants fall into {self.packet_count} packets"
            )
        )
        return (
            f"Symmetry-reduced misorientation (disorientation): {self.angle_deg:.3f} degrees "
            f"about ({axis}). This is the minimal representative over both point groups, which "
            "is what an EBSD boundary measurement reports and what the literature tabulates. "
            f"The relationship generates {self.variant_count} crystallographically distinct "
            f"variants; {packet_text}. Their distinct intervariant disorientation angles are "
            f"{spectrum} degrees — a measured misorientation histogram of children of one "
            "parent should show peaks at these values and nowhere else."
        )


def _misorientation_block(
    relationship: OrientationRelationship,
    variants: tuple[TransformationVariant, ...],
    packet_plane: CrystalPlane | None,
) -> ORDossierMisorientationBlock:
    misorientation = relationship.misorientation()
    if packet_plane is None:
        labels = np.zeros(0, dtype=np.int64)
        packet_count = 0
        packet_indices: tuple[int, int, int] | None = None
    else:
        labels = np.asarray(
            variant_close_packed_groups(relationship, packet_plane, variants=variants),
            dtype=np.int64,
        )
        packet_count = int(np.unique(labels).size)
        packet_indices = _plane_triple(packet_plane.miller.indices)
    angles = np.asarray(
        intervariant_misorientation_angles_deg(relationship, variants=variants),
        dtype=np.float64,
    )
    upper = angles[np.triu_indices(angles.shape[0], k=1)] if angles.size else np.zeros(0)
    distinct = np.unique(np.round(upper, 2)) if upper.size else np.zeros(0)
    return ORDossierMisorientationBlock(
        angle_deg=float(misorientation.angle_deg),
        axis=np.asarray(misorientation.rotation.axis, dtype=np.float64),
        variant_count=len(variants),
        packet_plane=packet_indices,
        packet_labels=labels + 1 if labels.size else labels,
        packet_count=packet_count,
        intervariant_angles_deg=distinct,
    )


# --------------------------------------------------------------------------- #
# Block 4 — the parallelisms
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORDossierParallelism:
    """One parallel pair, with the deviation that qualifies it.

    ``origin`` distinguishes the two kinds of row, and they must not be read the
    same way. A ``defining`` pair is the statement the relationship was
    constructed from, and its deviation is zero by construction. A
    ``discovered`` pair comes from the parallelism search, whose deviation is
    the **rationalization residual**: the exact child image of a parent plane is
    parallel to it by construction, so what is measured is the angle by which
    the nearest low-index child index misses that exact image.
    """

    kind: str
    origin: str
    parent_indices: tuple[int, int, int]
    child_indices: tuple[int, int, int]
    deviation_deg: float

    def __post_init__(self) -> None:
        if self.kind not in {"plane", "direction"}:
            raise ValueError("ORDossierParallelism.kind must be 'plane' or 'direction'.")
        if self.origin not in {"defining", "discovered"}:
            raise ValueError(
                "ORDossierParallelism.origin must be 'defining' or 'discovered'."
            )
        if not np.isfinite(self.deviation_deg) or self.deviation_deg < 0.0:
            raise ValueError("deviation_deg must be finite and non-negative.")

    @property
    def parent_label(self) -> str:
        """The parent object in publication notation."""

        formatter = format_plane_indices if self.kind == "plane" else format_direction_indices
        return formatter(self.parent_indices, style="plain")

    @property
    def child_label(self) -> str:
        """The child object in publication notation."""

        formatter = format_plane_indices if self.kind == "plane" else format_direction_indices
        return formatter(self.child_indices, style="plain")

    def to_json(self) -> dict[str, Any]:
        """The pair as JSON, indices as lists."""

        return {
            "kind": self.kind,
            "origin": self.origin,
            "parent_indices": [int(value) for value in self.parent_indices],
            "child_indices": [int(value) for value in self.child_indices],
            "parent_label": self.parent_label,
            "child_label": self.child_label,
            "deviation_deg": float(self.deviation_deg),
        }


@dataclass(frozen=True, slots=True)
class ORDossierParallelismBlock:
    """The defining parallelisms, and any discovered near-parallelisms."""

    pairs: tuple[ORDossierParallelism, ...]
    discovery_tolerance_deg: float

    @property
    def defining(self) -> tuple[ORDossierParallelism, ...]:
        """Only the pairs the relationship was declared with."""

        return tuple(pair for pair in self.pairs if pair.origin == "defining")

    @property
    def discovered(self) -> tuple[ORDossierParallelism, ...]:
        """Only the pairs found by the parallelism search."""

        return tuple(pair for pair in self.pairs if pair.origin == "discovered")

    def to_json(self) -> dict[str, Any]:
        """The block as JSON."""

        return {
            "discovery_tolerance_deg": float(self.discovery_tolerance_deg),
            "pairs": [pair.to_json() for pair in self.pairs],
        }

    def describe(self) -> str:
        """Prose summary, stating what each deviation column measures."""

        lines = [
            "Defining parallelisms of this variant (its own symmetry images of the "
            "declared pair, not the relationship's nominal pair):"
        ]
        for pair in self.defining:
            lines.append(
                f"  {pair.parent_label} || {pair.child_label} "
                f"({pair.deviation_deg:.4f} deg, zero by construction)"
            )
        discovered = self.discovered
        if not discovered:
            lines.append("No further families were nominated for the parallelism search.")
            return "\n".join(lines)
        lines.append(
            f"Discovered within {self.discovery_tolerance_deg:.3f} deg. Read this deviation "
            "precisely: the exact child image of a parent object is parallel to it by "
            "construction, so the angle reported is the rationalization residual — how far "
            "the nearest low-index child index sits from that exact image."
        )
        for pair in discovered:
            lines.append(
                f"  {pair.parent_label} || {pair.child_label} "
                f"({pair.deviation_deg:.4f} deg)"
            )
        return "\n".join(lines)


def _defining_pairs(variant: TransformationVariant) -> list[ORDossierParallelism]:
    pairs: list[ORDossierParallelism] = []
    for parent_plane, child_plane in variant.parallel_planes:
        pairs.append(
            ORDossierParallelism(
                kind="plane",
                origin="defining",
                parent_indices=_plane_triple(parent_plane.miller.indices),
                child_indices=_plane_triple(child_plane.miller.indices),
                deviation_deg=0.0,
            )
        )
    for parent_direction, child_direction in variant.parallel_directions:
        pairs.append(
            ORDossierParallelism(
                kind="direction",
                origin="defining",
                parent_indices=_index_triple(parent_direction.coordinates),
                child_indices=_index_triple(child_direction.coordinates),
                deviation_deg=0.0,
            )
        )
    return pairs


def _discovered_pairs(
    relationship: OrientationRelationship,
    variant: TransformationVariant,
    *,
    planes: tuple[CrystalPlane, ...],
    directions: tuple[CrystalDirection, ...],
    tolerance_deg: float,
) -> list[ORDossierParallelism]:
    pairs: list[ORDossierParallelism] = []
    for plane in planes:
        report = find_parallel_planes(
            relationship, plane, tolerance_deg=tolerance_deg, variants=(variant,)
        )
        for match in report.matches:
            pairs.append(
                ORDossierParallelism(
                    kind="plane",
                    origin="discovered",
                    parent_indices=_plane_triple(match.parent_indices),
                    child_indices=_plane_triple(match.child_indices),
                    deviation_deg=float(match.angular_deviation_deg),
                )
            )
    for direction in directions:
        report = find_parallel_directions(
            relationship, direction, tolerance_deg=tolerance_deg, variants=(variant,)
        )
        for match in report.matches:
            pairs.append(
                ORDossierParallelism(
                    kind="direction",
                    origin="discovered",
                    parent_indices=_index_triple(match.parent_indices),
                    child_indices=_index_triple(match.child_indices),
                    deviation_deg=float(match.angular_deviation_deg),
                )
            )
    return pairs


# --------------------------------------------------------------------------- #
# The dossier
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ORDossier:
    """Every number an orientation-relationship declaration implies, in one object.

    Purpose
    -------
    The reproducibility artifact of an OR analysis. Built by :func:`or_dossier`,
    it carries the four numeric blocks, explains itself through
    :meth:`describe`, serializes through :meth:`to_json` against
    ``schemas/or_dossier.schema.json``, and writes a whole bundle — numbers,
    tables and figures — through :meth:`export`.

    Attributes
    ----------
    relationship_name : str
    variant_index : int or None
        Which variant the transformation and parallelism blocks describe.
        ``None`` means the relationship as declared, which is variant 1.
    parent, child : ORDossierLatticeBlock
    transformation : ORDossierTransformationBlock
    misorientation : ORDossierMisorientationBlock
    parallelism : ORDossierParallelismBlock
    interface : None
        The F16 interface block, which is not implemented. It is carried as an
        explicit ``None`` rather than omitted, so a reader can tell "not
        analysed" from "analysed and empty".

    Notes
    -----
    Every value is a call into an existing function. The dossier adds
    aggregation, prose and serialization; it computes no crystallography of its
    own, which is what keeps it from ever disagreeing with the functions a
    reader would check it against.
    """

    relationship: OrientationRelationship
    variant: TransformationVariant
    variant_index: int | None
    parent: ORDossierLatticeBlock
    child: ORDossierLatticeBlock
    transformation: ORDossierTransformationBlock
    misorientation: ORDossierMisorientationBlock
    parallelism: ORDossierParallelismBlock
    interface: None = None

    @property
    def relationship_name(self) -> str:
        """The relationship's own name.

        A property rather than a stored copy: a dossier whose reported name
        could drift from the relationship it was built from would be a small
        version of exactly the defect this module exists to prevent.
        """

        return str(self.relationship.name)

    def describe(self) -> str:
        """Convention-explicit prose over every block, with its citations.

        The order is the order a reader needs: what the two crystals are, what
        the transformation does to the lattice, how the relationship reads as a
        measured misorientation, and what it holds parallel.
        """

        variant_text = (
            "the relationship as declared (variant 1)"
            if self.variant_index is None
            else f"variant {self.variant_index}"
        )
        return "\n\n".join(
            [
                f"Orientation-relationship dossier for '{self.relationship_name}', "
                f"{variant_text}. Angles are in degrees; indices are in the respective "
                "crystal bases, three-index form; matrices act on column vectors.",
                "Lattices.\n  " + self.parent.describe() + "\n  " + self.child.describe(),
                "Transformation.\n  " + self.transformation.describe(),
                "Misorientation and variants.\n  " + self.misorientation.describe(),
                "Parallelisms.\n" + self.parallelism.describe(),
                "Interface. Not analysed: interface crystallography is not implemented, "
                "so no habit plane, misfit or terrace decomposition is reported here. The "
                "absence is stated rather than left to be inferred from a missing section.",
                "Sources.\n  " + _CITATION_ITA + "\n  " + _CITATION_MORITO,
            ]
        )

    def to_json(self) -> dict[str, Any]:
        """The dossier as a JSON-ready dict, against the published schema.

        Keys and their meanings are fixed by ``schemas/or_dossier.schema.json``
        and are kept in lockstep with :meth:`describe`: a value that appears in
        one appears in the other.
        """

        from pytex import __version__

        return {
            "schema_id": OR_DOSSIER_SCHEMA_ID,
            "schema_version": OR_DOSSIER_SCHEMA_VERSION,
            "pytex_version": str(__version__),
            "relationship_name": self.relationship_name,
            "variant_index": self.variant_index,
            "parent": self.parent.to_json(),
            "child": self.child.to_json(),
            "transformation": self.transformation.to_json(),
            "misorientation": self.misorientation.to_json(),
            "parallelism": self.parallelism.to_json(),
            "interface": None,
        }

    def export(self, directory: Any, *, figures: bool = True) -> tuple[Path, ...]:
        """Write the whole bundle: numbers, tables and figures.

        Purpose
        -------
        Produces the directory another person needs in order to check the
        analysis: the JSON against the schema, the parallelism and spectrum
        tables as CSV *and* Markdown, the prose, and the figures as SVG.

        Parameters
        ----------
        directory : path-like
            Created if it does not exist.
        figures : bool
            Write the SVG figures. Set it false to write the numbers alone,
            which is the faster path when the bundle is being consumed by a
            program rather than read.

        Returns
        -------
        tuple of Path
            Every file written, in the order written.

        Notes
        -----
        The figures are the OR stereogram and the variant contact sheet, both of
        which need only the relationship. The variant pole figure and the
        per-variant SAED patterns are deliberately absent: they need a measured
        parent orientation and a diffraction geometry respectively, and neither
        is implied by a relationship.
        """

        target = Path(directory)
        target.mkdir(parents=True, exist_ok=True)
        written: list[Path] = []

        json_path = target / "or_dossier.json"
        json_path.write_text(
            json.dumps(self.to_json(), indent=2, sort_keys=False) + "\n", encoding="utf-8"
        )
        written.append(json_path)

        prose_path = target / "describe.md"
        prose_path.write_text(
            f"# {self.relationship_name}: orientation-relationship dossier\n\n"
            + self.describe()
            + "\n",
            encoding="utf-8",
        )
        written.append(prose_path)

        written.extend(self._write_parallelism_tables(target))
        written.extend(self._write_spectrum_table(target))
        if figures:
            written.extend(self._write_figures(target))
        return tuple(written)

    def _write_parallelism_tables(self, target: Path) -> list[Path]:
        rows = [pair.to_json() for pair in self.parallelism.pairs]
        fields = ["kind", "origin", "parent_label", "child_label", "deviation_deg"]
        csv_path = target / "parallelisms.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        markdown = [
            "| Kind | Origin | Parent | Child | Deviation / deg |",
            "| --- | --- | --- | --- | --- |",
        ]
        markdown.extend(
            f"| {row['kind']} | {row['origin']} | {row['parent_label']} | "
            f"{row['child_label']} | {row['deviation_deg']:.4f} |"
            for row in rows
        )
        markdown_path = target / "parallelisms.md"
        markdown_path.write_text("\n".join(markdown) + "\n", encoding="utf-8")
        return [csv_path, markdown_path]

    def _write_spectrum_table(self, target: Path) -> list[Path]:
        path = target / "intervariant_angles.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["angle_deg"])
            for angle in self.misorientation.intervariant_angles_deg:
                writer.writerow([f"{float(angle):.4f}"])
        return [path]

    def _write_figures(self, target: Path) -> list[Path]:
        """The SVG figures of the dossier.

        Both are drawn by the published plotting functions rather than
        assembled here, so a figure in the bundle is the same figure the
        interactive surfaces draw.
        """

        import matplotlib

        # Written from a script or a server as often as from a session, so the
        # backend is pinned to the file-writing one; an already-chosen backend
        # is left alone.
        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        from pytex.plotting.scene3d import WorldScene3D, render_variant_contact_sheet
        from pytex.plotting.spherical import plot_or_stereogram

        written: list[Path] = []
        stereogram_path = target / "or_stereogram.svg"
        figure = plot_or_stereogram(self.relationship, variant=self.variant)
        try:
            figure.savefig(stereogram_path, format="svg", bbox_inches="tight")
        finally:
            # A leaked figure is a defect in the test suite and a memory leak in
            # a long-running process; the same rule every renderer here follows.
            plt.close(figure)
        written.append(stereogram_path)

        variants = self.relationship.generate_variants()
        scenes = WorldScene3D.variant_scenes(self.relationship, variants=variants)
        sheet_path = target / "variant_contact_sheet.svg"
        sheet = render_variant_contact_sheet(
            scenes,
            variants=variants,
            columns=min(6, len(variants)),
            suptitle=f"{self.relationship_name}: {len(variants)} variants",
        )
        try:
            sheet.savefig(sheet_path, format="svg", bbox_inches="tight")
        finally:
            plt.close(sheet)
        written.append(sheet_path)
        return written


def or_dossier(
    relationship: OrientationRelationship,
    *,
    variant: int | TransformationVariant | None = None,
    planes: Any = (),
    directions: Any = (),
    packet_plane: CrystalPlane | None = None,
    tolerance_deg: float = DEFAULT_DISCOVERY_TOLERANCE_DEG,
) -> ORDossier:
    """Assemble the dossier of an orientation-relationship declaration.

    Purpose
    -------
    Turns "these two phases in the Kurdjumov-Sachs relationship" into every
    number that declaration implies, in one explainable, serializable object.

    Parameters
    ----------
    relationship : OrientationRelationship
    variant : int or TransformationVariant, optional
        Which variant the transformation and parallelism blocks describe.
        ``None`` means the relationship as declared; a one-based ``int`` indexes
        ``generate_variants()``. Variant 1 *is* the declaration, because its
        parent symmetry operator is the identity.
    planes, directions : sequence of CrystalPlane / CrystalDirection
        Parent families to run the parallelism search over, in addition to the
        defining pairs. Their deviations are rationalization residuals; see
        :class:`ORDossierParallelism`.
    packet_plane : CrystalPlane, optional
        The parent family the variants are grouped by. Defaults to the parent
        side of the relationship's first defining plane parallelism, which is
        the family the relationship is built on and therefore the grouping that
        reproduces the published packet structure. Pass a plane explicitly to
        group on another family, or nothing at all if the relationship declares
        no plane parallelism.
    tolerance_deg : float
        Tolerance of the parallelism search.

    Returns
    -------
    ORDossier

    Raises
    ------
    ValueError
        If ``variant`` is an out-of-range index.

    See Also
    --------
    ORDossier.describe : the prose form.
    ORDossier.export : the whole bundle on disk.
    """

    variants = relationship.generate_variants()
    if variant is None:
        resolved = variants[0]
        variant_index: int | None = None
    elif isinstance(variant, int):
        if not 1 <= variant <= len(variants):
            raise ValueError(
                f"variant must be a one-based index in 1..{len(variants)} for "
                f"'{relationship.name}'; got {variant}."
            )
        resolved = variants[variant - 1]
        variant_index = variant
    else:
        resolved = variant
        variant_index = int(variant.variant_index)

    resolved_packet_plane = packet_plane
    if resolved_packet_plane is None and relationship.parallel_planes:
        parent_plane = relationship.parallel_planes[0][0]
        resolved_packet_plane = CrystalPlane(
            MillerIndex(
                np.asarray(parent_plane.miller.indices, dtype=np.int64),
                phase=relationship.parent_phase,
            ),
            phase=relationship.parent_phase,
        )

    pairs = _defining_pairs(resolved)
    # A nominated family contains the defining member, so the search returns the
    # declared pair a second time. Listed twice it reads as two findings, and the
    # second would carry "discovered" where nothing was discovered.
    declared = {(pair.kind, pair.parent_indices, pair.child_indices) for pair in pairs}
    for pair in _discovered_pairs(
        relationship,
        resolved,
        planes=tuple(planes),
        directions=tuple(directions),
        tolerance_deg=tolerance_deg,
    ):
        if (pair.kind, pair.parent_indices, pair.child_indices) in declared:
            continue
        declared.add((pair.kind, pair.parent_indices, pair.child_indices))
        pairs.append(pair)

    return ORDossier(
        relationship=relationship,
        variant=resolved,
        variant_index=variant_index,
        parent=_lattice_block(relationship.parent_phase),
        child=_lattice_block(relationship.child_phase),
        transformation=_transformation_block(relationship, resolved),
        misorientation=_misorientation_block(relationship, variants, resolved_packet_plane),
        parallelism=ORDossierParallelismBlock(
            pairs=tuple(pairs), discovery_tolerance_deg=float(tolerance_deg)
        ),
    )
