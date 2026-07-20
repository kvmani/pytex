"""Composite orientation-relationship SAED pattern assembly.

Purpose: simulate the composite selected-area electron diffraction pattern
formed by a parent phase and any subset of the transformation variants of an
orientation relationship (OR), for an arbitrary parent zone axis. This is the
scientific core of OR verification in TEM practice: the composite pattern
shows which parent and child reflections coincide and how the variants
distribute in-plane.

Geometry conventions (pinned; see
``docs/roadmap/working_notes_composite_saed_program.md``):

- Variant rotations ``V_i`` re-express parent-crystal-frame Cartesian vectors
  in the child crystal frame (``TransformationVariant.map_parent_vector_to_child``).
- The shared detector basis ``(u, v, z)`` is anchored in the **parent**
  crystal frame; each child pattern is simulated with the rotated basis
  ``(V_i u, V_i v, V_i z)`` expressed in the child crystal frame, which is
  algebraically identical to pulling child reciprocal vectors back into the
  parent frame before projection. All phases therefore share one physically
  consistent detector.
- The child zone axis ``z_c = V_i z_p`` is in general irrational; the exact
  Cartesian direction drives the simulation, and the nearest rational
  ``[uvw]`` (with its angular deviation) is reported for labeling.
- Parent and child spot intensities are each normalized to a maximum of 1
  per sub-pattern; kinematic cross-phase intensity ratios are not defined at
  this level of theory and are a rendering choice.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array
from pytex.core.hexagonal import direction_uvw_to_uvtw, plane_hkl_to_hkil
from pytex.core.lattice import (
    CrystalDirection,
    MillerIndex,
    Phase,
    ReciprocalLatticeVector,
    ZoneAxis,
    phases_semantically_match,
)
from pytex.core.provenance import ProvenanceRecord
from pytex.core.transformation import OrientationRelationship, TransformationVariant
from pytex.diffraction.kinematic import (
    KinematicSimulationConfig,
    SpotTable,
    simulate_zone_axis_spots,
    zone_basis_from_axis,
)

# Zone-axis labels beyond index 6 are rarely meaningful in TEM practice; the
# nearest low-index zone with an honest deviation beats a closer high-index one.
_RATIONALIZE_DEFAULT_MAX_INDEX = 6


def _primitive_integer_triples(max_index: int) -> np.ndarray:
    axis = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    nonzero = grid[np.any(grid != 0, axis=1)]
    return np.asarray(nonzero[np.gcd.reduce(np.abs(nonzero), axis=1) == 1], dtype=np.int64)


def _reflection_label(hkl: np.ndarray, *, bravais: bool) -> str:
    indices = plane_hkl_to_hkil(hkl) if bravais else np.asarray(hkl, dtype=np.int64)
    return "(" + " ".join(str(int(value)) for value in indices) + ")"


def is_hexagonal_phase(phase: Phase) -> bool:
    """Whether ``phase`` belongs to the hexagonal crystal system.

    Purpose: decide when to present directions and reflections in four-index
    Miller-Bravais notation (the convention for hexagonal crystals such as
    the alpha-hcp product of the Burgers relationship). Trigonal phases are
    intentionally excluded: PyTex stores them in the hexagonal setting but
    three-index labels remain conventional there.
    """

    return phase.symmetry.to_point_group().crystal_system == "hexagonal"


@dataclass(frozen=True, slots=True)
class RationalizedZoneAxis:
    """Nearest rational ``[uvw]`` to an exact (possibly irrational) zone axis.

    ``deviation_deg`` is the true angular distance between the exact zone
    direction and the Cartesian image of ``indices``; zero means the zone
    axis is exactly rational within the searched index bound. For hexagonal
    phases ``indices_bravais`` carries the equivalent four-index Miller-Bravais
    ``[u v t w]`` form (``t = -(u + v)``), and :meth:`label` prefers it.
    """

    indices: np.ndarray
    deviation_deg: float
    max_index: int
    indices_bravais: np.ndarray | None = None

    def __post_init__(self) -> None:
        indices = as_int_array(self.indices, shape=(3,))
        if not np.any(indices):
            raise ValueError("RationalizedZoneAxis.indices must not be the zero triplet.")
        if not np.isfinite(self.deviation_deg) or self.deviation_deg < 0.0:
            raise ValueError("RationalizedZoneAxis.deviation_deg must be finite and >= 0.")
        if int(self.max_index) != self.max_index or self.max_index < 1:
            raise ValueError("RationalizedZoneAxis.max_index must be a positive integer.")
        indices = np.ascontiguousarray(indices)
        indices.setflags(write=False)
        object.__setattr__(self, "indices", indices)
        if self.indices_bravais is not None:
            bravais = as_int_array(self.indices_bravais, shape=(4,))
            if int(bravais[0] + bravais[1] + bravais[2]) != 0:
                raise ValueError(
                    "RationalizedZoneAxis.indices_bravais must satisfy U + V + T = 0."
                )
            bravais = np.ascontiguousarray(bravais)
            bravais.setflags(write=False)
            object.__setattr__(self, "indices_bravais", bravais)

    def label(self) -> str:
        indices = self.indices_bravais if self.indices_bravais is not None else self.indices
        return "[" + " ".join(str(int(value)) for value in indices) + "]"


def rationalize_zone_axis(
    direction: CrystalDirection, *, max_index: int = _RATIONALIZE_DEFAULT_MAX_INDEX
) -> RationalizedZoneAxis:
    """Nearest primitive integer zone axis to an exact crystal direction.

    Purpose: label the (generally irrational) child zone axes produced by
    mapping a parent zone axis through an OR variant. Candidates are compared
    by the true angle between Cartesian images under the direct basis, so
    the reported deviation is meaningful for any lattice metric.

    Inputs: ``direction`` — the exact zone direction (direct-basis
    coordinates may be non-integer); ``max_index`` — inclusive bound on the
    searched integer components (primitive triples only).

    Output: :class:`RationalizedZoneAxis` with the best ``[uvw]`` and its
    angular deviation in degrees. Sign-sensitive: the triple pointing along
    the direction (not its antipode) wins. For hexagonal phases the result
    also carries the four-index Miller-Bravais form.
    """

    if int(max_index) != max_index or max_index < 1:
        raise ValueError("max_index must be a positive integer.")
    basis_matrix = direction.phase.lattice.direct_basis().matrix
    target = basis_matrix @ direction.coordinates
    magnitude = float(np.linalg.norm(target))
    if np.isclose(magnitude, 0.0):
        raise ValueError("Cannot rationalize a direction with vanishing Cartesian image.")
    target_unit = target / magnitude
    candidates = _primitive_integer_triples(int(max_index))
    images = candidates.astype(np.float64) @ basis_matrix.T
    units = images / np.linalg.norm(images, axis=1)[:, None]
    cosines = units @ target_unit
    best = int(np.argmax(cosines))
    sine = float(np.linalg.norm(np.cross(units[best], target_unit)))
    deviation_deg = float(np.degrees(np.arctan2(sine, cosines[best])))
    best_indices = candidates[best].copy()
    bravais = (
        direction_uvw_to_uvtw(best_indices) if is_hexagonal_phase(direction.phase) else None
    )
    return RationalizedZoneAxis(
        indices=best_indices,
        deviation_deg=deviation_deg,
        max_index=int(max_index),
        indices_bravais=bravais,
    )


@dataclass(frozen=True, slots=True)
class VariantZonePattern:
    """One OR variant's zone-axis pattern inside a composite simulation.

    ``zone_axis_child`` is the exact beam direction re-expressed in the
    child crystal frame (Cartesian-derived, generally irrational);
    ``nearest_zone_axis`` is its rational label. ``spots`` are simulated in
    the shared parent-anchored detector basis, so their ``detector_mm``
    coordinates overlay the parent pattern directly.
    """

    variant: TransformationVariant
    zone_axis_child: CrystalDirection
    nearest_zone_axis: RationalizedZoneAxis
    spots: SpotTable

    def __post_init__(self) -> None:
        child_phase = self.variant.orientation_relationship.child_phase
        if not phases_semantically_match(self.zone_axis_child.phase, child_phase):
            raise ValueError("VariantZonePattern.zone_axis_child must use the child phase.")
        if not phases_semantically_match(self.spots.phase, child_phase):
            raise ValueError("VariantZonePattern.spots must belong to the child phase.")

    @property
    def variant_index(self) -> int:
        return self.variant.variant_index

    def label(self) -> str:
        deviation = self.nearest_zone_axis.deviation_deg
        suffix = "" if deviation < 0.05 else f" ({deviation:.2f} deg off)"
        return (
            f"V{self.variant_index}: {self.nearest_zone_axis.label()}{suffix}"
        )


@dataclass(frozen=True, slots=True)
class CompositeSAEDPattern:
    """Composite parent + OR-variant zone-axis diffraction pattern.

    All sub-patterns share the parent-anchored detector basis
    ``zone_basis_parent`` (columns ``u``, ``v``, zone axis; parent crystal
    Cartesian frame), so every ``detector_mm`` array overlays one detector.
    """

    relationship: OrientationRelationship
    parent_zone_axis: ZoneAxis
    parent_spots: SpotTable | None
    variant_patterns: tuple[VariantZonePattern, ...]
    zone_basis_parent: np.ndarray
    config: KinematicSimulationConfig
    child_config: KinematicSimulationConfig
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if not phases_semantically_match(
            self.parent_zone_axis.phase, self.relationship.parent_phase
        ):
            raise ValueError(
                "CompositeSAEDPattern.parent_zone_axis must use the parent phase."
            )
        basis = as_float_array(self.zone_basis_parent, shape=(3, 3))
        basis = np.ascontiguousarray(basis)
        basis.setflags(write=False)
        object.__setattr__(self, "zone_basis_parent", basis)
        object.__setattr__(self, "variant_patterns", tuple(self.variant_patterns))
        seen: set[int] = set()
        for pattern in self.variant_patterns:
            if pattern.variant_index in seen:
                raise ValueError(
                    "CompositeSAEDPattern variant indices must be unique; got repeated "
                    f"index {pattern.variant_index}."
                )
            seen.add(pattern.variant_index)

    @property
    def variant_indices(self) -> tuple[int, ...]:
        return tuple(pattern.variant_index for pattern in self.variant_patterns)

    def variant_pattern(self, variant_index: int) -> VariantZonePattern:
        for pattern in self.variant_patterns:
            if pattern.variant_index == variant_index:
                return pattern
        raise KeyError(
            f"Variant index {variant_index} is not part of this composite; available: "
            f"{list(self.variant_indices)}."
        )

    def select_variants(self, variant_indices: tuple[int, ...] | list[int]) -> CompositeSAEDPattern:
        """New composite restricted to the requested variant indices (order kept)."""

        selected = tuple(self.variant_pattern(int(index)) for index in variant_indices)
        return replace(self, variant_patterns=selected)

    def iter_spot_tables(self) -> Iterator[tuple[str, SpotTable]]:
        """Yield ``(label, spots)`` pairs: parent first, then each variant."""

        if self.parent_spots is not None:
            yield (f"{self.relationship.parent_phase.name}", self.parent_spots)
        for pattern in self.variant_patterns:
            yield (
                f"{self.relationship.child_phase.name} V{pattern.variant_index}",
                pattern.spots,
            )

    def all_detector_coordinates(self) -> np.ndarray:
        """Stacked ``(N, 2)`` detector coordinates (mm) over every sub-pattern."""

        blocks = [table.detector_mm for _, table in self.iter_spot_tables() if len(table)]
        if not blocks:
            return np.zeros((0, 2), dtype=np.float64)
        return np.vstack(blocks)

    def spot_count(self) -> int:
        return sum(len(table) for _, table in self.iter_spot_tables())

    def describe(self) -> str:
        """Convention-explicit prose summary of the composite pattern."""

        relationship = self.relationship
        parent_label = (
            "[" + " ".join(str(int(v)) for v in self.parent_zone_axis.indices) + "]"
        )
        parent_part = (
            f"{len(self.parent_spots)} parent reflection(s) included"
            if self.parent_spots is not None
            else "parent reflections not included"
        )
        lines = [
            (
                f"Composite kinematic SAED pattern for orientation relationship "
                f"'{relationship.name}' (parent '{relationship.parent_phase.name}' -> "
                f"child '{relationship.child_phase.name}') viewed along parent zone axis "
                f"{parent_label}; {parent_part}; {len(self.variant_patterns)} variant "
                f"pattern(s)."
            ),
            (
                "Conventions: the beam is antiparallel to the zone axis; variant rotations "
                "re-express parent-crystal Cartesian vectors in each child crystal frame; "
                "all sub-patterns share the parent-anchored detector basis, so detector "
                "coordinates overlay one detector. Reflections satisfy the excitation-error "
                f"criterion |s_g| <= {self.config.max_excitation_error_inv_angstrom:g} "
                "1/angstrom at "
                f"{self.config.beam_energy_kev:g} kV (camera constant "
                f"{self.config.camera_constant_mm_angstrom:g} mm*angstrom). Child zone axes "
                "are exact mapped directions; rational labels report their angular deviation. "
                "Each sub-pattern's intensities are max-normalized (kinematic cross-phase "
                "ratios are undefined at this level of theory)."
            ),
        ]
        for pattern in self.variant_patterns:
            lines.append(
                f"Variant {pattern.variant_index}: child zone axis nearest "
                f"{pattern.nearest_zone_axis.label()} (deviation "
                f"{pattern.nearest_zone_axis.deviation_deg:.3f} deg), "
                f"{len(pattern.spots)} reflection(s)."
            )
        return " ".join(lines)


def simulate_composite_saed(
    relationship: OrientationRelationship,
    parent_zone_axis: ZoneAxis,
    *,
    variant_indices: tuple[int, ...] | list[int] | None = None,
    include_parent: bool = True,
    config: KinematicSimulationConfig | None = None,
    child_config: KinematicSimulationConfig | None = None,
    align_parent_g: MillerIndex | None = None,
    in_plane_rotation_deg: float = 0.0,
    rationalize_max_index: int = _RATIONALIZE_DEFAULT_MAX_INDEX,
    provenance: ProvenanceRecord | None = None,
) -> CompositeSAEDPattern:
    """Simulate a composite parent + OR-variant SAED pattern.

    Purpose: the main entry point for OR composite diffraction. For each
    requested transformation variant the parent zone axis is mapped into the
    child crystal frame, the child pattern is simulated on the shared
    parent-anchored detector basis, and the child zone axis is rationalized
    for labeling.

    When to use: predicting or interpreting experimental SAED patterns of
    two-phase microstructures that obey an OR (martensite/austenite,
    precipitates, Burgers-related phases, ...), for any parent orientation
    or zone axis, with any subset of variants.

    Inputs: ``variant_indices`` — 1-based variant indices to include
    (``None`` = all variants from ``relationship.generate_variants()``);
    ``config`` — kinematic settings for the parent (defaults to
    :class:`KinematicSimulationConfig`); ``child_config`` — optional separate
    settings for children (defaults to ``config``); ``align_parent_g`` — a
    parent reflection whose in-plane component is placed along ``+u``;
    ``in_plane_rotation_deg`` — rotates the whole composite counterclockwise
    on the detector; ``rationalize_max_index`` — index bound for child
    zone-axis labels.

    Output: :class:`CompositeSAEDPattern` with parent and per-variant
    :class:`SpotTable` data in one shared detector frame.
    """

    if not phases_semantically_match(parent_zone_axis.phase, relationship.parent_phase):
        raise ValueError("parent_zone_axis.phase must match relationship.parent_phase.")
    parent_config = KinematicSimulationConfig() if config is None else config
    effective_child_config = parent_config if child_config is None else child_config

    align_cartesian: np.ndarray | None = None
    if align_parent_g is not None:
        if not phases_semantically_match(align_parent_g.phase, relationship.parent_phase):
            raise ValueError("align_parent_g.phase must match relationship.parent_phase.")
        align_cartesian = ReciprocalLatticeVector.from_miller_index(
            align_parent_g
        ).cartesian_vector

    zone_unit_parent = parent_zone_axis.unit_vector
    basis_parent = zone_basis_from_axis(
        zone_unit_parent,
        align_g_cartesian=align_cartesian,
        in_plane_rotation_deg=in_plane_rotation_deg,
    )

    parent_spots: SpotTable | None = None
    if include_parent:
        parent_spots = simulate_zone_axis_spots(
            relationship.parent_phase,
            parent_zone_axis,
            config=parent_config,
            basis=basis_parent,
            provenance=provenance,
        )

    variants = relationship.generate_variants()
    if variant_indices is not None:
        requested = [int(index) for index in variant_indices]
        available = {variant.variant_index: variant for variant in variants}
        missing = [index for index in requested if index not in available]
        if missing:
            raise ValueError(
                "Unknown variant indices "
                + ", ".join(str(index) for index in missing)
                + f"; available: 1..{len(variants)}."
            )
        variants = tuple(available[index] for index in requested)

    variant_patterns: list[VariantZonePattern] = []
    for variant in variants:
        rotation = variant.parent_to_child_rotation.as_matrix()
        basis_child = np.ascontiguousarray(rotation @ basis_parent)
        zone_child = CrystalDirection.from_cartesian(
            basis_child[:, 2], phase=relationship.child_phase
        )
        nearest = rationalize_zone_axis(zone_child, max_index=rationalize_max_index)
        spots = simulate_zone_axis_spots(
            relationship.child_phase,
            zone_child,
            config=effective_child_config,
            basis=basis_child,
            provenance=provenance,
        )
        variant_patterns.append(
            VariantZonePattern(
                variant=variant,
                zone_axis_child=zone_child,
                nearest_zone_axis=nearest,
                spots=spots,
            )
        )

    return CompositeSAEDPattern(
        relationship=relationship,
        parent_zone_axis=parent_zone_axis,
        parent_spots=parent_spots,
        variant_patterns=tuple(variant_patterns),
        zone_basis_parent=basis_parent,
        config=parent_config,
        child_config=effective_child_config,
        provenance=provenance,
    )


@dataclass(frozen=True, slots=True)
class SpotCoincidence:
    """One parent/child reflection pair landing at (nearly) the same detector spot.

    ``separation_mm`` is the Euclidean distance between the two spots on the
    shared detector; small separations are what make a composite pattern
    diagnostic of an OR in practice (superimposed reflections).
    ``parent_bravais``/``child_bravais`` request four-index Miller-Bravais
    rendering of the corresponding reflection label, as is conventional for
    hexagonal phases (the alpha-hcp product of the Burgers relationship).
    """

    variant_index: int
    parent_hkl: np.ndarray
    child_hkl: np.ndarray
    parent_detector_mm: np.ndarray
    child_detector_mm: np.ndarray
    separation_mm: float
    parent_bravais: bool = False
    child_bravais: bool = False

    def __post_init__(self) -> None:
        if self.variant_index <= 0:
            raise ValueError("SpotCoincidence.variant_index must be strictly positive.")
        parent_hkl = as_int_array(self.parent_hkl, shape=(3,))
        child_hkl = as_int_array(self.child_hkl, shape=(3,))
        parent_xy = as_float_array(self.parent_detector_mm, shape=(2,))
        child_xy = as_float_array(self.child_detector_mm, shape=(2,))
        if not np.isfinite(self.separation_mm) or self.separation_mm < 0.0:
            raise ValueError("SpotCoincidence.separation_mm must be finite and >= 0.")
        for array in (parent_hkl, child_hkl, parent_xy, child_xy):
            array.setflags(write=False)
        object.__setattr__(self, "parent_hkl", parent_hkl)
        object.__setattr__(self, "child_hkl", child_hkl)
        object.__setattr__(self, "parent_detector_mm", parent_xy)
        object.__setattr__(self, "child_detector_mm", child_xy)

    def label(self) -> str:
        parent = _reflection_label(self.parent_hkl, bravais=self.parent_bravais)
        child = _reflection_label(self.child_hkl, bravais=self.child_bravais)
        return (
            f"{parent}_p || {child}_V{self.variant_index} "
            f"({self.separation_mm:.2f} mm)"
        )


@dataclass(frozen=True, slots=True)
class SpotCoincidenceReport:
    """Quantitative statement of parent/child spot superposition for an OR.

    In TEM practice an OR is verified by exactly this: which child
    reflections fall on (or within a small distance of) parent reflections
    in the composite zone-axis pattern. ``coincidences`` lists every pair
    within ``tolerance_mm`` on the shared detector, per variant.
    """

    relationship_name: str
    parent_zone_label: str
    tolerance_mm: float
    parent_spot_count: int
    variant_spot_counts: tuple[tuple[int, int], ...]
    coincidences: tuple[SpotCoincidence, ...]
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.tolerance_mm) or self.tolerance_mm <= 0.0:
            raise ValueError("tolerance_mm must be finite and strictly positive.")
        if self.parent_spot_count < 0:
            raise ValueError("parent_spot_count must be non-negative.")
        object.__setattr__(self, "coincidences", tuple(self.coincidences))
        object.__setattr__(
            self,
            "variant_spot_counts",
            tuple((int(index), int(count)) for index, count in self.variant_spot_counts),
        )
        known = {index for index, _ in self.variant_spot_counts}
        for coincidence in self.coincidences:
            if coincidence.variant_index not in known:
                raise ValueError(
                    "SpotCoincidenceReport.coincidences reference variant "
                    f"{coincidence.variant_index} missing from variant_spot_counts."
                )
            if coincidence.separation_mm > self.tolerance_mm:
                raise ValueError(
                    "SpotCoincidenceReport.coincidences must respect tolerance_mm."
                )

    @property
    def variant_indices(self) -> tuple[int, ...]:
        return tuple(index for index, _ in self.variant_spot_counts)

    def count_for_variant(self, variant_index: int) -> int:
        return sum(
            1 for item in self.coincidences if item.variant_index == variant_index
        )

    def coincidences_for_variant(self, variant_index: int) -> tuple[SpotCoincidence, ...]:
        return tuple(
            item for item in self.coincidences if item.variant_index == variant_index
        )

    def describe(self) -> str:
        """Prose summary: per-variant superposition counts and their meaning."""

        per_variant = ", ".join(
            f"V{index}: {self.count_for_variant(index)}/{count}"
            for index, count in self.variant_spot_counts
        )
        strongest = "; ".join(item.label() for item in self.coincidences[:4])
        strongest_part = f" Closest pairs: {strongest}." if strongest else ""
        return (
            f"Spot-coincidence report for '{self.relationship_name}' along parent zone "
            f"{self.parent_zone_label}: {len(self.coincidences)} parent/child reflection "
            f"pair(s) separated by <= {self.tolerance_mm:g} mm on the shared detector "
            f"(parent pattern: {self.parent_spot_count} reflection(s); child reflections "
            f"coinciding per variant: {per_variant}).{strongest_part} Superimposed "
            "reflections are the practical TEM signature of the orientation "
            "relationship: variants with high coincidence fractions produce composite "
            "patterns that look like a single decorated zone-axis pattern."
        )


def find_spot_coincidences(
    pattern: CompositeSAEDPattern,
    *,
    tolerance_mm: float = 1.0,
    provenance: ProvenanceRecord | None = None,
) -> SpotCoincidenceReport:
    """Find parent/child reflections that superimpose on the shared detector.

    Purpose: quantify which child reflections of each variant land within
    ``tolerance_mm`` of a parent reflection — the measurable content of an
    OR in a composite SAED pattern.

    Inputs: ``pattern`` — a composite that includes parent spots;
    ``tolerance_mm`` — superposition radius on the detector (1 mm at the
    default 180 mm*angstrom camera constant corresponds to
    |Delta g| = 1/180 = 0.0056 1/angstrom).

    Output: :class:`SpotCoincidenceReport` with every qualifying pair,
    per-variant totals and a ``describe()`` narrative. Pair search uses a
    KD-tree per variant (vectorized; no O(N^2) Python loops).
    """

    from scipy.spatial import cKDTree

    if not np.isfinite(tolerance_mm) or tolerance_mm <= 0.0:
        raise ValueError("tolerance_mm must be finite and strictly positive.")
    parent_table = pattern.parent_spots
    if parent_table is None:
        raise ValueError(
            "find_spot_coincidences requires a composite simulated with "
            "include_parent=True."
        )
    coincidences: list[SpotCoincidence] = []
    variant_counts: list[tuple[int, int]] = []
    parent_is_hexagonal = is_hexagonal_phase(pattern.relationship.parent_phase)
    child_is_hexagonal = is_hexagonal_phase(pattern.relationship.child_phase)
    parent_tree = cKDTree(parent_table.detector_mm) if len(parent_table) else None
    for variant_pattern in pattern.variant_patterns:
        child_table = variant_pattern.spots
        variant_counts.append((variant_pattern.variant_index, len(child_table)))
        if parent_tree is None or not len(child_table):
            continue
        neighbor_lists = parent_tree.query_ball_point(
            child_table.detector_mm, r=tolerance_mm
        )
        for child_row, parent_rows in enumerate(neighbor_lists):
            child_xy = child_table.detector_mm[child_row]
            for parent_row in parent_rows:
                parent_xy = parent_table.detector_mm[parent_row]
                separation = float(np.linalg.norm(child_xy - parent_xy))
                coincidences.append(
                    SpotCoincidence(
                        variant_index=variant_pattern.variant_index,
                        parent_hkl=parent_table.hkl[parent_row],
                        child_hkl=child_table.hkl[child_row],
                        parent_detector_mm=parent_xy,
                        child_detector_mm=child_xy,
                        separation_mm=separation,
                        parent_bravais=parent_is_hexagonal,
                        child_bravais=child_is_hexagonal,
                    )
                )
    coincidences.sort(
        key=lambda item: (
            item.separation_mm,
            item.variant_index,
            tuple(int(v) for v in item.parent_hkl),
            tuple(int(v) for v in item.child_hkl),
        )
    )
    parent_zone_label = "[" + " ".join(
        str(int(v)) for v in pattern.parent_zone_axis.indices
    ) + "]"
    return SpotCoincidenceReport(
        relationship_name=pattern.relationship.name,
        parent_zone_label=parent_zone_label,
        tolerance_mm=float(tolerance_mm),
        parent_spot_count=len(parent_table),
        variant_spot_counts=tuple(variant_counts),
        coincidences=tuple(coincidences),
        provenance=provenance,
    )


def sweep_parent_zone_axes(
    relationship: OrientationRelationship,
    parent_zone_axes: Iterable[ZoneAxis],
    **kwargs: Any,
) -> Iterator[CompositeSAEDPattern]:
    """Lazily simulate composites for a sequence of parent zone axes.

    Purpose: zone-axis surveys (tilt-series planning, variant-visibility
    maps). Keyword arguments are forwarded unchanged to
    :func:`simulate_composite_saed`.
    """

    for zone_axis in parent_zone_axes:
        yield simulate_composite_saed(relationship, zone_axis, **kwargs)


__all__ = [
    "CompositeSAEDPattern",
    "RationalizedZoneAxis",
    "SpotCoincidence",
    "SpotCoincidenceReport",
    "VariantZonePattern",
    "find_spot_coincidences",
    "is_hexagonal_phase",
    "rationalize_zone_axis",
    "simulate_composite_saed",
    "sweep_parent_zone_axes",
]
