from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core.batches import VectorSet
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, Phase
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord


def _phase_semantically_matches(left: Phase | None, right: Phase) -> bool:
    if left is None:
        return False
    return (
        left.name == right.name
        and left.crystal_frame == right.crystal_frame
        and left.lattice == right.lattice
        and left.symmetry.point_group == right.symmetry.point_group
    )


@dataclass(frozen=True, slots=True)
class OrientationRelationship:
    name: str
    parent_phase: Phase
    child_phase: Phase
    parent_to_child_rotation: Rotation
    parallel_directions: tuple[tuple[np.ndarray, np.ndarray], ...] = ()
    parallel_planes: tuple[tuple[CrystalPlane, CrystalPlane], ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        normalized_name = self.name.strip()
        if not normalized_name:
            raise ValueError("OrientationRelationship.name must be non-empty.")
        if _phase_semantically_matches(self.parent_phase, self.child_phase):
            raise ValueError("OrientationRelationship requires distinct parent and child phases.")
        direction_pairs: list[tuple[np.ndarray, np.ndarray]] = []
        for parent_direction, child_direction in self.parallel_directions:
            parent_array = np.asarray(parent_direction, dtype=np.float64)
            child_array = np.asarray(child_direction, dtype=np.float64)
            if parent_array.shape != (3,) or child_array.shape != (3,):
                raise ValueError(
                    "OrientationRelationship.parallel_directions entries must each have shape (3,)."
                )
            if np.allclose(parent_array, 0.0) or np.allclose(child_array, 0.0):
                raise ValueError(
                    "OrientationRelationship.parallel_directions must not include zero vectors."
                )
            direction_pairs.append((parent_array, child_array))
        plane_pairs: list[tuple[CrystalPlane, CrystalPlane]] = []
        for parent_plane, child_plane in self.parallel_planes:
            if parent_plane.phase != self.parent_phase:
                raise ValueError(
                    "OrientationRelationship parallel parent planes must belong to parent_phase."
                )
            if child_plane.phase != self.child_phase:
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

    def generate_variants(self) -> tuple[TransformationVariant, ...]:
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
        variants: list[TransformationVariant] = []
        seen: set[tuple[float, float, float, float]] = set()
        for parent_index, parent_operator in enumerate(parent_operators):
            for child_index, child_operator in enumerate(child_operators):
                rotation = Rotation.from_matrix(
                    child_operator @ self.parent_to_child_rotation.as_matrix() @ parent_operator.T
                ).canonicalized()
                rounded = np.round(rotation.quaternion, decimals=12)
                key = (
                    float(rounded[0]),
                    float(rounded[1]),
                    float(rounded[2]),
                    float(rounded[3]),
                )
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
            if parent_plane.phase != self.orientation_relationship.parent_phase:
                raise ValueError(
                    "TransformationVariant parent habit planes must belong to parent_phase."
                )
            if child_plane.phase != self.orientation_relationship.child_phase:
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
        if not _phase_semantically_matches(
            self.parent_orientation.phase,
            self.orientation_relationship.parent_phase,
        ):
            raise ValueError(
                "PhaseTransformationRecord.parent_orientation.phase must match "
                "the relationship parent phase."
            )
        if not _phase_semantically_matches(
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
