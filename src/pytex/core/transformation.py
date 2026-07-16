from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core.batches import VectorSet
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import (
    CrystalDirection,
    CrystalPlane,
    MillerIndex,
    Phase,
    phases_semantically_match,
)
from pytex.core.orientation import (
    Orientation,
    OrientationSet,
    Rotation,
    _plane_direction_rotation_matrices,
    _reduced_pair_disorientation_angles,
)
from pytex.core.provenance import ProvenanceRecord


def _miller_index(values: tuple[int, int, int], *, phase: Phase) -> MillerIndex:
    return MillerIndex(np.asarray(values, dtype=np.int64), phase=phase)


def _crystal_direction(values: tuple[float, float, float], *, phase: Phase) -> CrystalDirection:
    return CrystalDirection(np.asarray(values, dtype=np.float64), phase=phase)


def _coerce_parallel_direction(
    entry: CrystalDirection | ArrayLike, *, phase: Phase, role: str
) -> CrystalDirection:
    """Normalize one parallel-direction entry to a phase-checked CrystalDirection.

    Raw 3-vectors are accepted for backward compatibility and are interpreted
    as Cartesian directions in the crystal frame of ``phase``; they are
    converted to direct-basis coordinates so the stored object keeps index
    meaning.
    """

    if isinstance(entry, CrystalDirection):
        if not phases_semantically_match(entry.phase, phase):
            raise ValueError(
                f"OrientationRelationship parallel {role} directions must belong to {role}_phase."
            )
        return entry
    vector = np.asarray(entry, dtype=np.float64)
    if vector.shape != (3,):
        raise ValueError(
            "OrientationRelationship.parallel_directions entries must each have shape (3,)."
        )
    if np.allclose(vector, 0.0):
        raise ValueError(
            "OrientationRelationship.parallel_directions must not include zero vectors."
        )
    return CrystalDirection.from_cartesian(vector, phase=phase)


def _require_proper_point_group(
    phase: Phase, expected: str, *, role: str, relationship: str
) -> None:
    if phase.symmetry.proper_point_group != expected:
        family = {"432": "cubic", "622": "hexagonal"}.get(expected, expected)
        raise ValueError(
            f"{relationship} correspondence requires a {family} {role} phase "
            f"with proper point group {expected}."
        )


def _require_cubic_phase_for_bain(phase: Phase, *, role: str) -> None:
    _require_proper_point_group(phase, "432", role=role, relationship="Bain")


@dataclass(frozen=True, slots=True)
class OrientationRelationship:
    name: str
    parent_phase: Phase
    child_phase: Phase
    parent_to_child_rotation: Rotation
    parallel_directions: tuple[tuple[CrystalDirection, CrystalDirection], ...] = ()
    parallel_planes: tuple[tuple[CrystalPlane, CrystalPlane], ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        normalized_name = self.name.strip()
        if not normalized_name:
            raise ValueError("OrientationRelationship.name must be non-empty.")
        if phases_semantically_match(self.parent_phase, self.child_phase):
            raise ValueError("OrientationRelationship requires distinct parent and child phases.")
        direction_pairs: list[tuple[CrystalDirection, CrystalDirection]] = []
        for parent_direction, child_direction in self.parallel_directions:
            direction_pairs.append(
                (
                    _coerce_parallel_direction(
                        parent_direction, phase=self.parent_phase, role="parent"
                    ),
                    _coerce_parallel_direction(
                        child_direction, phase=self.child_phase, role="child"
                    ),
                )
            )
        plane_pairs: list[tuple[CrystalPlane, CrystalPlane]] = []
        for parent_plane, child_plane in self.parallel_planes:
            if not phases_semantically_match(parent_plane.phase, self.parent_phase):
                raise ValueError(
                    "OrientationRelationship parallel parent planes must belong to parent_phase."
                )
            if not phases_semantically_match(child_plane.phase, self.child_phase):
                raise ValueError(
                    "OrientationRelationship parallel child planes must belong to child_phase."
                )
            plane_pairs.append((parent_plane, child_plane))
        object.__setattr__(self, "name", normalized_name)
        object.__setattr__(self, "parallel_directions", tuple(direction_pairs))
        object.__setattr__(self, "parallel_planes", tuple(plane_pairs))

    @property
    def parent_crystal_frame(self) -> ReferenceFrame:
        return self.parent_phase.crystal_frame

    @property
    def child_crystal_frame(self) -> ReferenceFrame:
        return self.child_phase.crystal_frame

    @classmethod
    def from_parallel_plane_direction(
        cls,
        *,
        name: str,
        parent_plane: CrystalPlane,
        child_plane: CrystalPlane,
        parent_direction: CrystalDirection,
        child_direction: CrystalDirection,
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        if not phases_semantically_match(parent_plane.phase, parent_direction.phase):
            raise ValueError("parent_plane.phase must match parent_direction.phase.")
        if not phases_semantically_match(child_plane.phase, child_direction.phase):
            raise ValueError("child_plane.phase must match child_direction.phase.")
        matrices = _plane_direction_rotation_matrices(
            crystal_normals=parent_plane.normal[None, :],
            crystal_directions=parent_direction.unit_vector[None, :],
            specimen_normals=child_plane.normal[None, :],
            specimen_directions=child_direction.unit_vector[None, :],
        )
        return cls(
            name=name,
            parent_phase=parent_plane.phase,
            child_phase=child_plane.phase,
            parent_to_child_rotation=Rotation.from_matrix(matrices[0]),
            parallel_directions=((parent_direction, child_direction),),
            parallel_planes=((parent_plane, child_plane),),
            provenance=provenance,
        )

    @classmethod
    def from_bain_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "bain",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        _require_cubic_phase_for_bain(parent_phase, role="parent")
        _require_cubic_phase_for_bain(child_phase, role="child")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(
                _miller_index((0, 0, 1), phase=parent_phase),
                phase=parent_phase,
            ),
            child_plane=CrystalPlane(
                _miller_index((0, 0, 1), phase=child_phase),
                phase=child_phase,
            ),
            parent_direction=_crystal_direction((1.0, 1.0, 0.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 0.0, 0.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_nishiyama_wassermann_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "nishiyama_wassermann",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        _require_cubic_phase_for_bain(parent_phase, role="parent")
        _require_cubic_phase_for_bain(child_phase, role="child")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(
                _miller_index((1, 1, 1), phase=parent_phase),
                phase=parent_phase,
            ),
            child_plane=CrystalPlane(
                _miller_index((0, 1, 1), phase=child_phase),
                phase=child_phase,
            ),
            parent_direction=_crystal_direction((1.0, -1.0, 0.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 0.0, 0.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_kurdjumov_sachs_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "kurdjumov_sachs",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Kurdjumov-Sachs OR: {111}_p || {011}_c, <-101>_p || <-1-11>_c.

        The classic fcc->bcc martensite relationship (24 variants). The
        representative pairing matches the Morito et al. V1 convention.
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Kurdjumov-Sachs"
        )
        _require_proper_point_group(
            child_phase, "432", role="child", relationship="Kurdjumov-Sachs"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((1, 1, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((0, 1, 1), phase=child_phase),
                                     phase=child_phase),
            parent_direction=_crystal_direction((-1.0, 0.0, 1.0), phase=parent_phase),
            child_direction=_crystal_direction((-1.0, -1.0, 1.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_greninger_troiano_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "greninger_troiano",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Greninger-Troiano OR: {111}_p || {011}_c, <5 12 17>_p || <7 17 17>_c.

        Sits between Kurdjumov-Sachs and Nishiyama-Wassermann: the
        representative below is 2.40 deg from KS and 2.86 deg from NW (the
        sign assignment was selected numerically so the pairing lands on the
        published GT orbit rather than a symmetry-inequivalent sibling).
        24 variants.
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Greninger-Troiano"
        )
        _require_proper_point_group(
            child_phase, "432", role="child", relationship="Greninger-Troiano"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((1, 1, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((0, 1, 1), phase=child_phase),
                                     phase=child_phase),
            parent_direction=_crystal_direction((-17.0, 5.0, 12.0), phase=parent_phase),
            child_direction=_crystal_direction((-7.0, 17.0, -17.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_pitsch_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "pitsch",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Pitsch OR: {001}_p || {-101}_c, <110>_p || <111>_c.

        The thin-film fcc->bcc relationship (12 variants); its representative
        lies 5.26 deg from Kurdjumov-Sachs, mirroring the KS-NW separation.
        """

        _require_proper_point_group(parent_phase, "432", role="parent", relationship="Pitsch")
        _require_proper_point_group(child_phase, "432", role="child", relationship="Pitsch")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((0, 0, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((-1, 0, 1), phase=child_phase),
                                     phase=child_phase),
            parent_direction=_crystal_direction((1.0, 1.0, 0.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 1.0, 1.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_burgers_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "burgers",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Burgers OR: {110}_bcc || (0001)_hcp, <-111>_bcc || <11-20>_hcp.

        The bcc->hcp transformation relationship (beta->alpha titanium,
        zirconium; 12 variants). The parent must be cubic (proper group 432)
        and the child hexagonal (proper group 622).
        """

        _require_proper_point_group(parent_phase, "432", role="parent", relationship="Burgers")
        _require_proper_point_group(child_phase, "622", role="child", relationship="Burgers")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((1, 1, 0), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=child_phase),
            parent_direction=_crystal_direction((-1.0, 1.0, 1.0), phase=parent_phase),
            child_direction=CrystalDirection.from_miller_bravais(
                (1, 1, -2, 0), phase=child_phase
            ),
            provenance=provenance,
        )

    def map_parent_vector_to_child(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        if isinstance(vector, VectorSet):
            if vector.reference_frame != self.parent_crystal_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match OrientationRelationship.parent_phase."
                )
            matrix = self.parent_to_child_rotation.as_matrix()
            return VectorSet(
                values=vector.values @ matrix.T,
                reference_frame=self.child_crystal_frame,
                provenance=vector.provenance,
            )
        return self.parent_to_child_rotation.apply(vector)

    def map_child_vector_to_parent(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        inverse = self.parent_to_child_rotation.inverse()
        if isinstance(vector, VectorSet):
            if vector.reference_frame != self.child_crystal_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match OrientationRelationship.child_phase."
                )
            matrix = inverse.as_matrix()
            return VectorSet(
                values=vector.values @ matrix.T,
                reference_frame=self.parent_crystal_frame,
                provenance=vector.provenance,
            )
        return inverse.apply(vector)

    def inverse(self, *, name: str | None = None) -> OrientationRelationship:
        return OrientationRelationship(
            name=name or f"{self.child_phase.name}_to_{self.parent_phase.name}",
            parent_phase=self.child_phase,
            child_phase=self.parent_phase,
            parent_to_child_rotation=self.parent_to_child_rotation.inverse(),
            parallel_directions=tuple(
                (child_direction, parent_direction)
                for parent_direction, child_direction in self.parallel_directions
            ),
            parallel_planes=tuple(
                (child_plane, parent_plane) for parent_plane, child_plane in self.parallel_planes
            ),
            provenance=self.provenance,
        )

    def generate_variants(
        self, *, reduce_by_child_symmetry: bool = True
    ) -> tuple[TransformationVariant, ...]:
        """Enumerate the transformation variants of this relationship.

        With ``reduce_by_child_symmetry`` (the default), variants are the
        crystallographically distinct child orientations produced from one
        fixed parent orientation: parent symmetry operators generate the
        candidate rotations ``R S_p^T``, and candidates that differ only by a
        child symmetry operator collapse into a single variant. This
        reproduces the literature variant counts (Bain 3; NW, Pitsch, Burgers
        12; KS, GT 24). Set it to ``False`` for the historical raw
        enumeration of distinct operator products ``S_c R S_p^T``, which
        counts every symmetry-equivalent description separately.
        """

        parent_symmetry = self.parent_phase.symmetry
        child_symmetry = self.child_phase.symmetry
        parent_operators = (
            parent_symmetry.operators
            if parent_symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        child_operators = (
            child_symmetry.operators
            if child_symmetry is not None
            else np.eye(3, dtype=np.float64)[None, :, :]
        )
        base_matrix = self.parent_to_child_rotation.as_matrix()

        def quaternion_key(matrix: np.ndarray) -> tuple[float, float, float, float]:
            # q and -q describe the same rotation: canonicalize by taking the
            # lexicographically larger of the two ROUNDED sign choices. This
            # stays stable for 180-degree rotations (w ~ 0) and for
            # equal-magnitude components, where pivot- or w-based sign
            # conventions flip on floating-point noise.
            quaternion = Rotation.from_matrix(matrix).quaternion
            rounded_pos = np.round(quaternion, 10) + 0.0
            rounded_neg = np.round(-quaternion, 10) + 0.0
            positive = (
                float(rounded_pos[0]),
                float(rounded_pos[1]),
                float(rounded_pos[2]),
                float(rounded_pos[3]),
            )
            negative = (
                float(rounded_neg[0]),
                float(rounded_neg[1]),
                float(rounded_neg[2]),
                float(rounded_neg[3]),
            )
            return max(positive, negative)

        variants: list[TransformationVariant] = []
        seen: set[tuple[float, float, float, float]] = set()
        if reduce_by_child_symmetry:
            for parent_index, parent_operator in enumerate(parent_operators):
                candidate = base_matrix @ parent_operator.T
                orbit_key = min(
                    quaternion_key(child_operator @ candidate)
                    for child_operator in child_operators
                )
                if orbit_key in seen:
                    continue
                seen.add(orbit_key)
                variants.append(
                    TransformationVariant(
                        orientation_relationship=self,
                        variant_index=len(variants) + 1,
                        parent_operator_index=parent_index,
                        child_operator_index=0,
                        parent_to_child_rotation=Rotation.from_matrix(candidate).canonicalized(),
                        provenance=self.provenance,
                    )
                )
            return tuple(variants)
        for parent_index, parent_operator in enumerate(parent_operators):
            for child_index, child_operator in enumerate(child_operators):
                rotation = Rotation.from_matrix(
                    child_operator @ base_matrix @ parent_operator.T
                ).canonicalized()
                key = quaternion_key(rotation.as_matrix())
                if key in seen:
                    continue
                seen.add(key)
                variants.append(
                    TransformationVariant(
                        orientation_relationship=self,
                        variant_index=len(variants) + 1,
                        parent_operator_index=parent_index,
                        child_operator_index=child_index,
                        parent_to_child_rotation=rotation,
                        provenance=self.provenance,
                    )
                )
        return tuple(variants)


@dataclass(frozen=True, slots=True)
class TransformationVariant:
    orientation_relationship: OrientationRelationship
    variant_index: int
    parent_operator_index: int
    child_operator_index: int
    parent_to_child_rotation: Rotation
    habit_plane_pairs: tuple[tuple[CrystalPlane, CrystalPlane], ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.variant_index <= 0:
            raise ValueError("TransformationVariant.variant_index must be strictly positive.")
        if self.parent_operator_index < 0 or self.child_operator_index < 0:
            raise ValueError("TransformationVariant operator indices must be non-negative.")
        plane_pairs: list[tuple[CrystalPlane, CrystalPlane]] = []
        for parent_plane, child_plane in self.habit_plane_pairs:
            if not phases_semantically_match(
                parent_plane.phase, self.orientation_relationship.parent_phase
            ):
                raise ValueError(
                    "TransformationVariant parent habit planes must belong to parent_phase."
                )
            if not phases_semantically_match(
                child_plane.phase, self.orientation_relationship.child_phase
            ):
                raise ValueError(
                    "TransformationVariant child habit planes must belong to child_phase."
                )
            plane_pairs.append((parent_plane, child_plane))
        object.__setattr__(self, "habit_plane_pairs", tuple(plane_pairs))

    def map_parent_vector_to_child(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        if isinstance(vector, VectorSet):
            if vector.reference_frame != self.orientation_relationship.parent_crystal_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match TransformationVariant.parent_phase."
                )
            matrix = self.parent_to_child_rotation.as_matrix()
            return VectorSet(
                values=vector.values @ matrix.T,
                reference_frame=self.orientation_relationship.child_crystal_frame,
                provenance=vector.provenance,
            )
        return self.parent_to_child_rotation.apply(vector)


@dataclass(frozen=True, slots=True)
class IntervariantMisorientation:
    """The symmetry-reduced misorientation between two transformation variants.

    ``axis_child_frame`` is the unit rotation axis of the minimal
    (disorientation) representative, expressed in the child crystal frame.
    """

    variant_a: int
    variant_b: int
    angle_deg: float
    axis_child_frame: np.ndarray

    def __post_init__(self) -> None:
        if self.variant_a <= 0 or self.variant_b <= 0:
            raise ValueError("IntervariantMisorientation variant indices must be positive.")
        axis = np.asarray(self.axis_child_frame, dtype=np.float64)
        if axis.shape != (3,):
            raise ValueError("axis_child_frame must have shape (3,).")
        norm = float(np.linalg.norm(axis))
        if not np.isclose(norm, 1.0, atol=1e-8):
            raise ValueError("axis_child_frame must be a unit vector.")
        axis = np.ascontiguousarray(axis)
        axis.setflags(write=False)
        object.__setattr__(self, "axis_child_frame", axis)


def _child_operators(relationship: OrientationRelationship) -> np.ndarray:
    symmetry = relationship.child_phase.symmetry
    if symmetry is None:
        return np.eye(3, dtype=np.float64)[None, :, :]
    return np.asarray(symmetry.operators, dtype=np.float64)


def intervariant_misorientation_angles_deg(
    relationship: OrientationRelationship,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
) -> np.ndarray:
    """Pairwise disorientation angles (deg) between all transformation variants.

    Entry ``[i, j]`` is the child-symmetry-reduced misorientation angle
    between variants ``i`` and ``j`` (zero diagonal, symmetric); the row/column
    order follows ``generate_variants()``. For Kurdjumov-Sachs this reproduces
    the published intervariant table (Morito et al.): angles from 10.53 deg up
    to 60.00 deg.
    """

    resolved = relationship.generate_variants() if variants is None else variants
    matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in resolved], axis=0
    )
    count = matrices.shape[0]
    relative = np.einsum("ipq,jrq->ijpr", matrices, matrices, optimize=True).reshape(
        count * count, 3, 3
    )
    operators = _child_operators(relationship)
    angles = _reduced_pair_disorientation_angles(relative, operators, operators)
    result = np.degrees(angles.reshape(count, count))
    result = np.ascontiguousarray(result)
    result.setflags(write=False)
    return result


def intervariant_misorientations(
    relationship: OrientationRelationship,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
) -> tuple[IntervariantMisorientation, ...]:
    """Disorientation angle and axis for every unordered variant pair.

    Returns one `IntervariantMisorientation` per pair ``a < b`` in
    ``generate_variants()`` order, with the axis of the minimal
    symmetry-reduced representative in the child crystal frame.
    """

    resolved = relationship.generate_variants() if variants is None else variants
    operators = _child_operators(relationship)
    results: list[IntervariantMisorientation] = []
    for a in range(len(resolved)):
        matrix_a = resolved[a].parent_to_child_rotation.as_matrix()
        for b in range(a + 1, len(resolved)):
            matrix_b = resolved[b].parent_to_child_rotation.as_matrix()
            relative = matrix_a @ matrix_b.T
            products = np.einsum(
                "aij,jk,blk->abil", operators, relative, operators, optimize=True
            )
            traces = np.trace(products, axis1=2, axis2=3)
            cosines = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
            flat = np.arccos(cosines).reshape(-1)
            best = int(np.argmin(flat))
            representative = Rotation.from_matrix(products.reshape(-1, 3, 3)[best])
            results.append(
                IntervariantMisorientation(
                    variant_a=resolved[a].variant_index,
                    variant_b=resolved[b].variant_index,
                    angle_deg=representative.angle_deg,
                    axis_child_frame=representative.axis,
                )
            )
    return tuple(results)


@dataclass(frozen=True, slots=True)
class PhaseTransformationRecord:
    name: str
    orientation_relationship: OrientationRelationship
    parent_orientation: Orientation
    child_orientations: OrientationSet
    variant_indices: np.ndarray | None = None
    provenance: ProvenanceRecord | None = None
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        normalized_name = self.name.strip()
        if not normalized_name:
            raise ValueError("PhaseTransformationRecord.name must be non-empty.")
        if not phases_semantically_match(
            self.parent_orientation.phase,
            self.orientation_relationship.parent_phase,
        ):
            raise ValueError(
                "PhaseTransformationRecord.parent_orientation.phase must match "
                "the relationship parent phase."
            )
        if not phases_semantically_match(
            self.child_orientations.phase,
            self.orientation_relationship.child_phase,
        ):
            raise ValueError(
                "PhaseTransformationRecord.child_orientations.phase must match "
                "the relationship child phase."
            )
        if self.parent_orientation.specimen_frame != self.child_orientations.specimen_frame:
            raise ValueError("PhaseTransformationRecord orientations must share a specimen frame.")
        if self.variant_indices is not None:
            indices = as_int_array(self.variant_indices, shape=(None,))
            if indices.shape != (len(self.child_orientations),):
                raise ValueError(
                    "PhaseTransformationRecord.variant_indices must have one "
                    "entry per child orientation."
                )
            if np.any(indices <= 0):
                raise ValueError(
                    "PhaseTransformationRecord.variant_indices must be strictly positive."
                )
            object.__setattr__(self, "variant_indices", indices)
        object.__setattr__(self, "name", normalized_name)
        object.__setattr__(self, "notes", tuple(str(note) for note in self.notes))

    @property
    def variant_count(self) -> int:
        if self.variant_indices is None:
            return 0
        return int(np.unique(self.variant_indices).size)

    def predicted_child_orientations(self) -> OrientationSet:
        if self.variant_indices is None:
            predicted_rotations = [self.orientation_relationship.parent_to_child_rotation] * len(
                self.child_orientations
            )
        else:
            variant_lookup = {
                variant.variant_index: variant.parent_to_child_rotation
                for variant in self.orientation_relationship.generate_variants()
            }
            missing = sorted(
                {
                    int(index)
                    for index in np.unique(self.variant_indices)
                    if int(index) not in variant_lookup
                }
            )
            if missing:
                raise ValueError(
                    "PhaseTransformationRecord.variant_indices contain values not produced by "
                    "OrientationRelationship.generate_variants(): "
                    + ", ".join(str(value) for value in missing)
                )
            variant_indices = np.asarray(self.variant_indices, dtype=np.int64)
            predicted_rotations = [
                variant_lookup[int(index)] for index in variant_indices
            ]
        quaternions = np.stack(
            [
                predicted_rotation.compose(self.parent_orientation.rotation).quaternion
                for predicted_rotation in predicted_rotations
            ],
            axis=0,
        )
        return OrientationSet(
            quaternions=quaternions,
            crystal_frame=self.child_orientations.crystal_frame,
            specimen_frame=self.child_orientations.specimen_frame,
            symmetry=self.child_orientations.symmetry,
            phase=self.child_orientations.phase,
            provenance=self.provenance or self.parent_orientation.provenance,
        )
