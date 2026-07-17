from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

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
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.orientation import (
    Misorientation,
    Orientation,
    OrientationSet,
    Rotation,
    _plane_direction_rotation_matrices,
    _reduced_pair_disorientation_angles,
    matrices_to_quaternions,
)
from pytex.core.provenance import ProvenanceRecord


def _miller_index(values: tuple[int, int, int], *, phase: Phase) -> MillerIndex:
    return MillerIndex(np.asarray(values, dtype=np.int64), phase=phase)


def _crystal_direction(values: tuple[float, float, float], *, phase: Phase) -> CrystalDirection:
    return CrystalDirection(np.asarray(values, dtype=np.float64), phase=phase)


def _index_tuple(values: np.ndarray) -> tuple[int, int, int]:
    array = np.asarray(values, dtype=np.int64).reshape(3)
    return (int(array[0]), int(array[1]), int(array[2]))


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


#: Bound of the integer-triple search used when rationalizing mapped indices.
#: 17 covers the full standard OR catalog, including the Greninger-Troiano
#: <5 12 17> direction family.
DEFAULT_RATIONALIZATION_MAX_INDEX = 17


@lru_cache(maxsize=8)
def _primitive_integer_triples(max_index: int) -> np.ndarray:
    """All primitive (gcd = 1) signed integer triples with entries in [-max_index, max_index]."""

    if max_index < 1:
        raise ValueError("max_index must be at least 1.")
    axis = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    nonzero = grid[np.any(grid != 0, axis=1)]
    primitive = nonzero[np.gcd.reduce(np.abs(nonzero), axis=1) == 1]
    primitive = np.ascontiguousarray(primitive)
    primitive.setflags(write=False)
    return primitive


def _rationalize_components(
    components: np.ndarray, *, basis_matrix: np.ndarray, max_index: int
) -> tuple[np.ndarray, float]:
    """Nearest primitive integer triple to real basis components, by angle.

    Candidate triples are compared with the exact components through their
    Cartesian images under ``basis_matrix`` (direct basis for directions,
    reciprocal basis for plane indices), so the returned residual is the true
    angular deviation in the relevant space. The match is sign-sensitive: the
    triple whose image points closest to the exact image wins.
    """

    target = basis_matrix @ np.asarray(components, dtype=np.float64)
    magnitude = float(np.linalg.norm(target))
    if np.isclose(magnitude, 0.0):
        raise ValueError("Cannot rationalize components with a vanishing Cartesian image.")
    target_unit = target / magnitude
    candidates = _primitive_integer_triples(max_index)
    images = candidates.astype(np.float64) @ basis_matrix.T
    units = images / np.linalg.norm(images, axis=1)[:, None]
    cosines = units @ target_unit
    best = int(np.argmax(cosines))
    # atan2 keeps full precision for near-zero angles where arccos floors out.
    sine = float(np.linalg.norm(np.cross(units[best], target_unit)))
    residual_deg = float(np.degrees(np.arctan2(sine, cosines[best])))
    return candidates[best].copy(), residual_deg


def _index_correspondence(
    components: np.ndarray,
    *,
    rotation: np.ndarray,
    source_basis: np.ndarray,
    target_basis: np.ndarray,
    max_index: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    exact = np.linalg.solve(
        target_basis, rotation @ (source_basis @ np.asarray(components, dtype=np.float64))
    )
    rational, residual_deg = _rationalize_components(
        exact, basis_matrix=target_basis, max_index=max_index
    )
    exact = np.ascontiguousarray(exact)
    exact.setflags(write=False)
    return exact, rational, residual_deg


@dataclass(frozen=True, slots=True)
class DirectionCorrespondence:
    """The image of a crystal direction under an orientation relationship.

    ``target_exact_coordinates`` are the generally irrational direct-basis
    components of the mapped direction in the target phase;
    ``rational_indices`` is the nearest primitive integer ``[uvw]`` (the
    coordinates of ``target``), and ``angular_residual_deg`` is the real-space
    angle between the exact image and its rationalization. ``variant_index``
    is set when the mapping went through a specific transformation variant.
    """

    source: CrystalDirection
    target: CrystalDirection
    target_exact_coordinates: np.ndarray
    rational_indices: np.ndarray
    angular_residual_deg: float
    variant_index: int | None = None

    def __post_init__(self) -> None:
        exact = np.asarray(self.target_exact_coordinates, dtype=np.float64)
        rational = np.asarray(self.rational_indices, dtype=np.int64)
        if exact.shape != (3,) or rational.shape != (3,):
            raise ValueError("DirectionCorrespondence components must have shape (3,).")
        if not np.isfinite(self.angular_residual_deg) or self.angular_residual_deg < 0.0:
            raise ValueError("angular_residual_deg must be finite and non-negative.")
        for array in (exact, rational):
            array.setflags(write=False)
        object.__setattr__(self, "target_exact_coordinates", exact)
        object.__setattr__(self, "rational_indices", rational)

    def describe(self) -> str:
        """Prose summary: source direction, exact image, rationalization, residual."""

        source_idx = np.rint(self.source.coordinates).astype(np.int64)
        source_text = (
            format_direction_indices(_index_tuple(source_idx), style="plain")
            if np.allclose(self.source.coordinates, source_idx, atol=1e-8)
            else np.array2string(self.source.coordinates, precision=4)
        )
        variant_text = (
            f" through variant {self.variant_index}" if self.variant_index is not None else ""
        )
        exact = self.target_exact_coordinates
        rational_text = format_direction_indices(
            _index_tuple(self.rational_indices), style="plain"
        )
        return (
            f"Direction correspondence{variant_text}: parent {source_text} maps to exact "
            f"child components [{exact[0]:.4f} {exact[1]:.4f} {exact[2]:.4f}], rationalized "
            f"to {rational_text} with an "
            f"angular residual of {self.angular_residual_deg:.4f} deg (real-space angle; "
            "indices in the target crystal direct basis)."
        )


@dataclass(frozen=True, slots=True)
class PlaneCorrespondence:
    """The image of a crystal plane under an orientation relationship.

    ``target_exact_indices`` are the generally irrational reciprocal-basis
    components ``(hkl)`` of the mapped plane in the target phase;
    ``rational_indices`` is the nearest primitive integer ``(hkl)`` (the
    Miller indices of ``target``), and ``angular_residual_deg`` is the angle
    between the exact and rationalized plane normals. ``variant_index`` is set
    when the mapping went through a specific transformation variant.
    """

    source: CrystalPlane
    target: CrystalPlane
    target_exact_indices: np.ndarray
    rational_indices: np.ndarray
    angular_residual_deg: float
    variant_index: int | None = None

    def __post_init__(self) -> None:
        exact = np.asarray(self.target_exact_indices, dtype=np.float64)
        rational = np.asarray(self.rational_indices, dtype=np.int64)
        if exact.shape != (3,) or rational.shape != (3,):
            raise ValueError("PlaneCorrespondence components must have shape (3,).")
        if not np.isfinite(self.angular_residual_deg) or self.angular_residual_deg < 0.0:
            raise ValueError("angular_residual_deg must be finite and non-negative.")
        for array in (exact, rational):
            array.setflags(write=False)
        object.__setattr__(self, "target_exact_indices", exact)
        object.__setattr__(self, "rational_indices", rational)

    def describe(self) -> str:
        """Prose summary: source plane, exact image, rationalization, residual."""

        variant_text = (
            f" through variant {self.variant_index}" if self.variant_index is not None else ""
        )
        exact = self.target_exact_indices
        source_text = format_plane_indices(
            _index_tuple(self.source.miller.indices), style="plain"
        )
        rational_text = format_plane_indices(_index_tuple(self.rational_indices), style="plain")
        return (
            f"Plane correspondence{variant_text}: parent "
            f"{source_text} maps to exact "
            f"child components ({exact[0]:.4f} {exact[1]:.4f} {exact[2]:.4f}), rationalized "
            f"to {rational_text} with an angular "
            f"residual of {self.angular_residual_deg:.4f} deg (angle between plane normals; "
            "indices in the target crystal reciprocal basis)."
        )


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
    def from_shoji_nishiyama_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "shoji_nishiyama",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Shoji-Nishiyama OR: {111}_fcc || (0001)_hcp, <-110>_fcc || <11-20>_hcp.

        The fcc->hcp epsilon-martensite relationship (austenite to
        epsilon in high-Mn steels and Co; 4 variants, one per {111}
        close-packed parent plane). The parent must be cubic (proper group
        432) and the child hexagonal (proper group 622).
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Shoji-Nishiyama"
        )
        _require_proper_point_group(
            child_phase, "622", role="child", relationship="Shoji-Nishiyama"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((1, 1, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=child_phase),
            parent_direction=_crystal_direction((-1.0, 1.0, 0.0), phase=parent_phase),
            child_direction=CrystalDirection.from_miller_bravais(
                (1, 1, -2, 0), phase=child_phase
            ),
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

    @classmethod
    def from_pitsch_schrader_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "pitsch_schrader",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Pitsch-Schrader OR: (0001)_hcp || {110}_bcc, <11-20>_hcp || <001>_bcc.

        The hcp->bcc relationship of the ferrite/alpha systems (3 distinct
        variants from one hexagonal parent); its representative sits 5.26 deg
        from the inverse Burgers relationship — the hexagonal analogue of the
        KS-Pitsch separation. The parent must be hexagonal (proper group 622)
        and the child cubic (proper group 432).
        """

        _require_proper_point_group(
            parent_phase, "622", role="parent", relationship="Pitsch-Schrader"
        )
        _require_proper_point_group(
            child_phase, "432", role="child", relationship="Pitsch-Schrader"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((1, 1, 0), phase=child_phase),
                                     phase=child_phase),
            parent_direction=CrystalDirection.from_miller_bravais(
                (1, 1, -2, 0), phase=parent_phase
            ),
            child_direction=_crystal_direction((0.0, 0.0, 1.0), phase=child_phase),
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

    def _resolve_variant_rotation(
        self, variant: TransformationVariant | None
    ) -> tuple[np.ndarray, int | None]:
        if variant is None:
            return self.parent_to_child_rotation.as_matrix(), None
        relationship = variant.orientation_relationship
        if relationship is not self and not (
            relationship.name == self.name
            and phases_semantically_match(relationship.parent_phase, self.parent_phase)
            and phases_semantically_match(relationship.child_phase, self.child_phase)
        ):
            raise ValueError(
                "variant must belong to this OrientationRelationship "
                "(same name, parent phase, and child phase)."
            )
        return variant.parent_to_child_rotation.as_matrix(), variant.variant_index

    def correspondence_direct(
        self, *, variant: TransformationVariant | None = None
    ) -> np.ndarray:
        """Direction-index correspondence matrix (parent ``[uvw]`` to child ``[uvw]``).

        Purpose: the linear map that carries direct-basis direction components
        across the relationship: ``u_child = M @ u_parent`` with
        ``M = A_child^-1 R A_parent`` built from the phases' direct structure
        matrices and the (variant) parent-to-child rotation. It is generally
        not a rotation matrix and generally irrational; use
        ``map_direction_to_child`` for rationalized indices with residuals.

        Inputs: optionally a ``TransformationVariant`` of this relationship.

        Output: a read-only ``(3, 3)`` float matrix.
        """

        rotation, _ = self._resolve_variant_rotation(variant)
        parent_basis = self.parent_phase.lattice.direct_basis().matrix
        child_basis = self.child_phase.lattice.direct_basis().matrix
        matrix = np.linalg.solve(child_basis, rotation @ parent_basis)
        matrix = np.ascontiguousarray(matrix)
        matrix.setflags(write=False)
        return matrix

    def correspondence_reciprocal(
        self, *, variant: TransformationVariant | None = None
    ) -> np.ndarray:
        """Plane-index correspondence matrix (parent ``(hkl)`` to child ``(hkl)``).

        Purpose: the linear map that carries reciprocal-basis plane components
        across the relationship: ``h_child = M* @ h_parent`` with
        ``M* = B_child^-1 R B_parent`` built from the reciprocal structure
        matrices. It equals the inverse-transpose of
        ``correspondence_direct``, which preserves the zone law
        ``h . u`` across the mapping.

        Inputs: optionally a ``TransformationVariant`` of this relationship.

        Output: a read-only ``(3, 3)`` float matrix.
        """

        rotation, _ = self._resolve_variant_rotation(variant)
        parent_basis = self.parent_phase.lattice.reciprocal_basis().matrix
        child_basis = self.child_phase.lattice.reciprocal_basis().matrix
        matrix = np.linalg.solve(child_basis, rotation @ parent_basis)
        matrix = np.ascontiguousarray(matrix)
        matrix.setflags(write=False)
        return matrix

    def map_direction_to_child(
        self,
        direction: CrystalDirection,
        *,
        variant: TransformationVariant | None = None,
        max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    ) -> DirectionCorrespondence:
        """Map a parent crystal direction to its child-phase counterpart.

        Purpose: answers "which child ``[uvw]`` corresponds to this parent
        direction" for the relationship (or one of its variants), returning
        the exact (irrational) child components, the nearest primitive integer
        indices within ``max_index``, and the angular residual between them.

        Inputs: a ``CrystalDirection`` belonging to the parent phase;
        optionally a variant and the rationalization index bound.

        Output: a ``DirectionCorrespondence``.
        """

        if not phases_semantically_match(direction.phase, self.parent_phase):
            raise ValueError("direction.phase must match the relationship parent phase.")
        rotation, variant_index = self._resolve_variant_rotation(variant)
        return self._map_direction(
            direction,
            rotation=rotation,
            source_phase=self.parent_phase,
            target_phase=self.child_phase,
            max_index=max_index,
            variant_index=variant_index,
        )

    def map_direction_to_parent(
        self,
        direction: CrystalDirection,
        *,
        variant: TransformationVariant | None = None,
        max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    ) -> DirectionCorrespondence:
        """Map a child crystal direction back to its parent-phase counterpart.

        The inverse of ``map_direction_to_child``: use it to interpret
        product-phase measurements against parent-frame stereography.
        """

        if not phases_semantically_match(direction.phase, self.child_phase):
            raise ValueError("direction.phase must match the relationship child phase.")
        rotation, variant_index = self._resolve_variant_rotation(variant)
        return self._map_direction(
            direction,
            rotation=rotation.T,
            source_phase=self.child_phase,
            target_phase=self.parent_phase,
            max_index=max_index,
            variant_index=variant_index,
        )

    def _map_direction(
        self,
        direction: CrystalDirection,
        *,
        rotation: np.ndarray,
        source_phase: Phase,
        target_phase: Phase,
        max_index: int,
        variant_index: int | None,
    ) -> DirectionCorrespondence:
        exact, rational, residual_deg = _index_correspondence(
            direction.coordinates,
            rotation=rotation,
            source_basis=source_phase.lattice.direct_basis().matrix,
            target_basis=target_phase.lattice.direct_basis().matrix,
            max_index=max_index,
        )
        return DirectionCorrespondence(
            source=direction,
            target=CrystalDirection(rational.astype(np.float64), phase=target_phase),
            target_exact_coordinates=exact,
            rational_indices=rational,
            angular_residual_deg=residual_deg,
            variant_index=variant_index,
        )

    def map_plane_to_child(
        self,
        plane: CrystalPlane,
        *,
        variant: TransformationVariant | None = None,
        max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    ) -> PlaneCorrespondence:
        """Map a parent crystal plane to its child-phase counterpart.

        Purpose: answers "which child ``(hkl)`` corresponds to this parent
        plane" for the relationship (or one of its variants), returning the
        exact (irrational) child plane components, the nearest primitive
        integer Miller indices within ``max_index``, and the angular residual
        between the exact and rationalized plane normals.

        Inputs: a ``CrystalPlane`` belonging to the parent phase; optionally a
        variant and the rationalization index bound.

        Output: a ``PlaneCorrespondence``.
        """

        if not phases_semantically_match(plane.phase, self.parent_phase):
            raise ValueError("plane.phase must match the relationship parent phase.")
        rotation, variant_index = self._resolve_variant_rotation(variant)
        return self._map_plane(
            plane,
            rotation=rotation,
            source_phase=self.parent_phase,
            target_phase=self.child_phase,
            max_index=max_index,
            variant_index=variant_index,
        )

    def map_plane_to_parent(
        self,
        plane: CrystalPlane,
        *,
        variant: TransformationVariant | None = None,
        max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    ) -> PlaneCorrespondence:
        """Map a child crystal plane back to its parent-phase counterpart.

        The inverse of ``map_plane_to_child``: use it for habit-plane trace
        analysis and for reading product-phase diffraction against the parent.
        """

        if not phases_semantically_match(plane.phase, self.child_phase):
            raise ValueError("plane.phase must match the relationship child phase.")
        rotation, variant_index = self._resolve_variant_rotation(variant)
        return self._map_plane(
            plane,
            rotation=rotation.T,
            source_phase=self.child_phase,
            target_phase=self.parent_phase,
            max_index=max_index,
            variant_index=variant_index,
        )

    def _map_plane(
        self,
        plane: CrystalPlane,
        *,
        rotation: np.ndarray,
        source_phase: Phase,
        target_phase: Phase,
        max_index: int,
        variant_index: int | None,
    ) -> PlaneCorrespondence:
        exact, rational, residual_deg = _index_correspondence(
            plane.miller.indices.astype(np.float64),
            rotation=rotation,
            source_basis=source_phase.lattice.reciprocal_basis().matrix,
            target_basis=target_phase.lattice.reciprocal_basis().matrix,
            max_index=max_index,
        )
        return PlaneCorrespondence(
            source=plane,
            target=CrystalPlane(
                MillerIndex(rational, phase=target_phase), phase=target_phase
            ),
            target_exact_indices=exact,
            rational_indices=rational,
            angular_residual_deg=residual_deg,
            variant_index=variant_index,
        )

    def describe(self) -> str:
        """Convention-explicit prose summary of the relationship.

        States the parent/child phases and point groups, the defining
        parallelisms, the symmetry-reduced misorientation representative, and
        the distinct-variant count. Angles are in degrees; indices are in the
        respective crystal bases (three-index form).
        """

        parent_group = (
            self.parent_phase.symmetry.point_group
            if self.parent_phase.symmetry is not None
            else "1"
        )
        child_group = (
            self.child_phase.symmetry.point_group
            if self.child_phase.symmetry is not None
            else "1"
        )
        lines = [
            f"Orientation relationship '{self.name}': parent phase "
            f"'{self.parent_phase.name}' ({parent_group}) to child phase "
            f"'{self.child_phase.name}' ({child_group})."
        ]
        for parent_plane, child_plane in self.parallel_planes:
            parent_text = format_plane_indices(
                _index_tuple(parent_plane.miller.indices), style="plain"
            )
            child_text = format_plane_indices(
                _index_tuple(child_plane.miller.indices), style="plain"
            )
            lines.append(
                f"Defining plane parallelism: {parent_text} parent || {child_text} child."
            )
        for parent_direction, child_direction in self.parallel_directions:
            parent_idx = np.rint(parent_direction.coordinates).astype(np.int64)
            child_idx = np.rint(child_direction.coordinates).astype(np.int64)
            if np.allclose(parent_direction.coordinates, parent_idx, atol=1e-8) and np.allclose(
                child_direction.coordinates, child_idx, atol=1e-8
            ):
                parent_text = format_direction_indices(_index_tuple(parent_idx), style="plain")
                child_text = format_direction_indices(_index_tuple(child_idx), style="plain")
                lines.append(
                    f"Defining direction parallelism: {parent_text} parent || {child_text} child."
                )
        misorientation = self.misorientation()
        axis = misorientation.rotation.axis
        lines.append(
            "Misorientation representative (child-symmetry and parent-symmetry "
            f"reduced): {misorientation.angle_deg:.2f} deg about "
            f"<{axis[0]:.3f} {axis[1]:.3f} {axis[2]:.3f}> "
            "(axis components identical in both crystal frames)."
        )
        lines.append(
            f"Variants: {len(self.generate_variants())} crystallographically distinct "
            "child orientations from one parent (child-symmetry reduced)."
        )
        return "\n".join(lines)

    def deformation_gradient(
        self,
        *,
        variant: TransformationVariant | None = None,
    ) -> DeformationGradientReport:
        """Lattice-correspondence deformation gradient (Bain-strain analysis).

        Purpose: the third object of the OR doctrine — the physical lattice
        distortion of the transformation. The exact index correspondence's
        images of the three parent basis vectors are rationalized to the
        integer lattice correspondence; the parent-frame deformation is then
        ``F = R^T A_c M_rat A_p^-1`` (rigid rotation removed), whose
        symmetric right-stretch gives the principal transformation strains.
        For Bain with a_fcc = 3.6 and a_bcc = 2.87 the principal stretches
        are the textbook (1.127, 1.127, 0.797) with a +1.3% volume change;
        every KS/NW/GT variant shares the same principal stretches because
        they differ from Bain only by rigid rotation.

        Inputs: optionally a variant and the rationalization bound.

        Output: a ``DeformationGradientReport`` (see its ``describe()``).
        """

        rotation, variant_index = self._resolve_variant_rotation(variant)
        parent_basis = self.parent_phase.lattice.direct_basis().matrix
        child_basis = self.child_phase.lattice.direct_basis().matrix
        exact = np.linalg.solve(child_basis, rotation @ parent_basis)
        # The lattice correspondence is the nearest INTEGER matrix: magnitudes
        # matter for strain, so ray-based rationalization would be wrong here.
        correspondence = np.rint(exact).astype(np.int64)
        if int(round(float(np.linalg.det(correspondence.astype(np.float64))))) == 0:
            raise ValueError(
                "The nearest-integer lattice correspondence is singular; the "
                "relationship's exact correspondence is too far from an integer "
                "matrix for Bain-strain analysis."
            )
        component_error = float(np.max(np.abs(exact - correspondence)))
        gradient = (
            rotation.T
            @ child_basis
            @ correspondence.astype(np.float64)
            @ np.linalg.inv(parent_basis)
        )
        right_cauchy_green = gradient.T @ gradient
        eigenvalues, eigenvectors = np.linalg.eigh(right_cauchy_green)
        order = np.argsort(eigenvalues)[::-1]
        stretches = np.sqrt(eigenvalues[order])
        directions = eigenvectors[:, order].T
        stretch_tensor = eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
        polar_rotation = gradient @ np.linalg.inv(stretch_tensor)
        polar_angle = float(
            _rotation_angles_deg_from_matrices(polar_rotation[None, :, :])[0]
        )
        return DeformationGradientReport(
            relationship_name=self.name,
            variant_index=variant_index,
            deformation_gradient=gradient,
            stretch_tensor=stretch_tensor,
            principal_stretches=stretches,
            principal_directions=directions,
            volume_ratio=float(np.linalg.det(gradient)),
            correspondence=correspondence,
            polar_rotation_deg=polar_angle,
            correspondence_max_component_error=component_error,
            provenance=self.provenance,
        )

    def misorientation(self) -> Misorientation:
        """The relationship as a symmetry-reduced misorientation (disorientation).

        Purpose: expresses the OR the way it is measured and reported in the
        literature — a minimal-angle axis/angle representative reduced by the
        child symmetry (left) and parent symmetry (right). For Kurdjumov-Sachs
        this is the published ~42.85 deg rotation about a <0.968 0.178 0.178>
        axis. The representative's axis components are the same in both
        crystal frames because the axis is the fixed eigenvector of the map.

        Output: a ``Misorientation`` whose ``rotation`` is the deterministic
        fundamental-zone representative; use ``.angle_deg`` and
        ``.rotation.axis`` for reporting, and compare against boundary
        misorientations from EBSD data.
        """

        return Misorientation(
            rotation=self.parent_to_child_rotation,
            left_symmetry=self.child_phase.symmetry,
            right_symmetry=self.parent_phase.symmetry,
            provenance=self.provenance,
        ).disorientation()

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
    count = len(resolved)
    if count < 2:
        return ()
    matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in resolved], axis=0
    )
    row_indices, column_indices = np.triu_indices(count, k=1)
    relative = np.einsum(
        "nij,nkj->nik", matrices[row_indices], matrices[column_indices], optimize=True
    )
    # One symmetry-product tensor for every pair at once; axis order (a, b)
    # matches the historical per-pair enumeration so representative selection
    # is unchanged.
    products = np.einsum(
        "aij,njk,blk->nabil", operators, relative, operators, optimize=True
    )
    traces = np.trace(products, axis1=3, axis2=4)
    cosines = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    angles = np.arccos(cosines).reshape(len(row_indices), -1)
    best = np.argmin(angles, axis=1)
    representatives = products.reshape(len(row_indices), -1, 3, 3)[
        np.arange(len(row_indices)), best
    ]
    results: list[IntervariantMisorientation] = []
    for pair, (a, b) in enumerate(zip(row_indices, column_indices, strict=True)):
        representative = Rotation.from_matrix(representatives[pair])
        results.append(
            IntervariantMisorientation(
                variant_a=resolved[int(a)].variant_index,
                variant_b=resolved[int(b)].variant_index,
                angle_deg=representative.angle_deg,
                axis_child_frame=representative.axis,
            )
        )
    return tuple(results)


@dataclass(frozen=True, slots=True)
class ORDeviationReport:
    """Angular deviation of measured parent/child pairs from a nominal OR.

    ``deviations_deg[i]`` is the child-symmetry-reduced misorientation angle
    between observed child ``i`` and the closest variant prediction from
    parent ``i``; ``best_variant_indices[i]`` is that variant's 1-based index.
    A perfectly obeyed relationship gives zeros; systematic offsets measure
    how far the operative OR sits from the nominal one (e.g. children built
    with Greninger-Troiano deviate ~2.4 deg from Kurdjumov-Sachs).
    """

    relationship_name: str
    deviations_deg: np.ndarray
    best_variant_indices: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        deviations = np.asarray(self.deviations_deg, dtype=np.float64).reshape(-1)
        indices = np.asarray(self.best_variant_indices, dtype=np.int64).reshape(-1)
        if deviations.shape != indices.shape:
            raise ValueError(
                "deviations_deg and best_variant_indices must have the same length."
            )
        if deviations.size == 0:
            raise ValueError("ORDeviationReport requires at least one pair.")
        if np.any(~np.isfinite(deviations)) or np.any(deviations < 0.0):
            raise ValueError("deviations_deg must be finite and non-negative.")
        if np.any(indices <= 0):
            raise ValueError("best_variant_indices must be strictly positive.")
        for array in (deviations, indices):
            array.setflags(write=False)
        object.__setattr__(self, "deviations_deg", deviations)
        object.__setattr__(self, "best_variant_indices", indices)

    @property
    def mean_deviation_deg(self) -> float:
        return float(np.mean(self.deviations_deg))

    @property
    def median_deviation_deg(self) -> float:
        return float(np.median(self.deviations_deg))

    @property
    def max_deviation_deg(self) -> float:
        return float(np.max(self.deviations_deg))

    def describe(self) -> str:
        """Prose summary: pair count, deviation statistics, variants used."""

        used = np.unique(self.best_variant_indices)
        variant_list = ", ".join(str(int(index)) for index in used[:12])
        if used.size > 12:
            variant_list += ", ..."
        return (
            f"Deviation from orientation relationship '{self.relationship_name}' over "
            f"{self.deviations_deg.size} parent/child pair(s): mean "
            f"{self.mean_deviation_deg:.3f} deg, median {self.median_deviation_deg:.3f} deg, "
            f"max {self.max_deviation_deg:.3f} deg (child-symmetry-reduced angles, minimum "
            f"over variants). Best-matching variants used: {variant_list}. Deviations near "
            "zero mean the relationship is obeyed; systematic offsets measure the distance "
            "of the operative relationship from the nominal one."
        )


def or_deviation(
    parent_orientations: OrientationSet,
    child_orientations: OrientationSet,
    relationship: OrientationRelationship,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> ORDeviationReport:
    """Deviation of measured parent/child orientation pairs from a nominal OR.

    Purpose: the quantitative test of "does this transformation follow OR X".
    For each pair ``i``, every variant prediction ``V_k g_parent_i`` is
    compared with the observed child under the child symmetry, and the
    smallest disorientation angle (with the winning variant index) is
    reported. Zero deviations mean the relationship is obeyed exactly;
    the aggregate statistics quantify systematic departure and feed OR
    fitting.

    Inputs: paired ``OrientationSet`` objects (equal length, shared specimen
    frame, phases matching the relationship's parent and child), and
    optionally a precomputed variant tuple.

    Output: an ``ORDeviationReport``.
    """

    if len(parent_orientations) != len(child_orientations):
        raise ValueError(
            "parent_orientations and child_orientations must be paired (equal length)."
        )
    if len(parent_orientations) == 0:
        raise ValueError("or_deviation requires at least one orientation pair.")
    if not phases_semantically_match(parent_orientations.phase, relationship.parent_phase):
        raise ValueError("parent_orientations.phase must match the relationship parent phase.")
    if not phases_semantically_match(child_orientations.phase, relationship.child_phase):
        raise ValueError("child_orientations.phase must match the relationship child phase.")
    if parent_orientations.specimen_frame != child_orientations.specimen_frame:
        raise ValueError("Parent and child orientations must share a specimen frame.")
    resolved = relationship.generate_variants() if variants is None else variants
    variant_matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in resolved], axis=0
    )
    parent_matrices = parent_orientations.as_matrices()
    child_matrices = child_orientations.as_matrices()
    predicted = np.einsum("nij,vkj->nvik", parent_matrices, variant_matrices, optimize=True)
    relative = np.einsum(
        "nji,nvjk->nvik", child_matrices, predicted, optimize=True
    )
    pair_count, variant_count = relative.shape[0], relative.shape[1]
    child_symmetry = relationship.child_phase.symmetry
    operators = (
        child_symmetry.operators
        if child_symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    angles = _reduced_pair_disorientation_angles(
        relative.reshape(pair_count * variant_count, 3, 3), operators, operators
    ).reshape(pair_count, variant_count)
    best_columns = np.argmin(angles, axis=1)
    variant_indices = np.array(
        [resolved[int(column)].variant_index for column in best_columns], dtype=np.int64
    )
    deviations_deg = np.degrees(angles[np.arange(pair_count), best_columns])
    return ORDeviationReport(
        relationship_name=relationship.name,
        deviations_deg=deviations_deg,
        best_variant_indices=variant_indices,
        provenance=provenance or relationship.provenance,
    )


@dataclass(frozen=True, slots=True)
class OrientationRelationshipFitReport:
    """Result of fitting an orientation relationship to measured pairs.

    ``relationship`` is the fitted OR (same phases as the nominal one, name
    suffixed ``_fitted``, defining parallelisms deliberately not carried over
    because the fit is matrix-level). ``residuals_deg`` are the per-pair
    angles between each symmetry-aligned measurement and the fitted rotation;
    ``deviation_from_nominal_deg`` is the symmetry-reduced distance between
    the fitted and nominal relationships.
    """

    relationship: OrientationRelationship
    nominal_name: str
    residuals_deg: np.ndarray
    iterations: int
    converged: bool
    deviation_from_nominal_deg: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        residuals = np.asarray(self.residuals_deg, dtype=np.float64).reshape(-1)
        if residuals.size == 0:
            raise ValueError("OrientationRelationshipFitReport requires at least one pair.")
        if np.any(~np.isfinite(residuals)) or np.any(residuals < 0.0):
            raise ValueError("residuals_deg must be finite and non-negative.")
        if self.iterations <= 0:
            raise ValueError("iterations must be positive.")
        residuals.setflags(write=False)
        object.__setattr__(self, "residuals_deg", residuals)

    @property
    def mean_residual_deg(self) -> float:
        return float(np.mean(self.residuals_deg))

    @property
    def max_residual_deg(self) -> float:
        return float(np.max(self.residuals_deg))

    def describe(self) -> str:
        """Prose summary: fit quality, convergence, and distance from nominal."""

        convergence = (
            f"converged in {self.iterations} iteration(s)"
            if self.converged
            else f"did NOT converge within {self.iterations} iteration(s)"
        )
        misorientation = self.relationship.misorientation()
        axis = misorientation.rotation.axis
        return (
            f"Fitted orientation relationship '{self.relationship.name}' from "
            f"{self.residuals_deg.size} parent/child pair(s), starting from nominal "
            f"'{self.nominal_name}': {convergence}. Fitted misorientation representative: "
            f"{misorientation.angle_deg:.2f} deg about "
            f"<{axis[0]:.3f} {axis[1]:.3f} {axis[2]:.3f}>. Per-pair residuals: mean "
            f"{self.mean_residual_deg:.3f} deg, max {self.max_residual_deg:.3f} deg "
            f"(symmetry-aligned angles). The fit sits "
            f"{self.deviation_from_nominal_deg:.3f} deg from the nominal relationship; "
            "a large value means the operative relationship differs systematically "
            "from the assumed one."
        )


def _rotation_angles_deg_from_matrices(matrices: np.ndarray) -> np.ndarray:
    traces = np.trace(matrices, axis1=-2, axis2=-1)
    cosines = np.clip((traces - 1.0) * 0.5, -1.0, 1.0)
    return np.asarray(np.degrees(np.arccos(cosines)), dtype=np.float64)


def _symmetry_reduced_angle_between_deg(
    left: np.ndarray,
    right: np.ndarray,
    *,
    child_operators: np.ndarray,
    parent_operators: np.ndarray,
) -> float:
    candidates = np.einsum(
        "aij,jk,bkl->abil", child_operators, left, parent_operators, optimize=True
    )
    relative = np.einsum("abij,kj->abik", candidates, right, optimize=True)
    return float(np.min(_rotation_angles_deg_from_matrices(relative.reshape(-1, 3, 3))))


def fit_orientation_relationship(
    parent_orientations: OrientationSet,
    child_orientations: OrientationSet,
    nominal: OrientationRelationship,
    *,
    max_iterations: int = 20,
    convergence_tol_deg: float = 1e-8,
    provenance: ProvenanceRecord | None = None,
) -> OrientationRelationshipFitReport:
    """Fit the operative orientation relationship to measured pairs.

    Purpose: estimates the parent-to-child rotation that best explains
    measured parent/child orientation pairs — the step beyond
    ``or_deviation``: instead of only quantifying departure from a nominal
    relationship, it refines the relationship itself (e.g. starting from
    Kurdjumov-Sachs and recovering the operative Greninger-Troiano-like OR).

    Algorithm: each pair's measured crystal-to-crystal map ``C P^T`` is
    aligned to the current estimate through the parent and child symmetry
    groups (the equivalent description nearest the estimate), the aligned
    rotations are averaged with the quaternion eigen-mean (Markley), and the
    align/average steps iterate to convergence. This is the standard
    symmetry-aware rotation-averaging route used for parent/child OR
    refinement (MTEX ``calcParent2Child`` is the parity reference).

    Inputs: paired ``OrientationSet`` objects matching the nominal
    relationship's phases and sharing a specimen frame; the nominal
    relationship supplies the starting estimate, phases, and symmetry groups.

    Output: an ``OrientationRelationshipFitReport`` (see its ``describe()``).
    """

    if len(parent_orientations) != len(child_orientations):
        raise ValueError(
            "parent_orientations and child_orientations must be paired (equal length)."
        )
    if len(parent_orientations) == 0:
        raise ValueError("fit_orientation_relationship requires at least one pair.")
    if not phases_semantically_match(parent_orientations.phase, nominal.parent_phase):
        raise ValueError("parent_orientations.phase must match the nominal parent phase.")
    if not phases_semantically_match(child_orientations.phase, nominal.child_phase):
        raise ValueError("child_orientations.phase must match the nominal child phase.")
    if parent_orientations.specimen_frame != child_orientations.specimen_frame:
        raise ValueError("Parent and child orientations must share a specimen frame.")
    parent_symmetry = nominal.parent_phase.symmetry
    child_symmetry = nominal.child_phase.symmetry
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
    parent_matrices = parent_orientations.as_matrices()
    child_matrices = child_orientations.as_matrices()
    # Canonical crystal->specimen convention: C = P V^T, so the measured
    # parent-to-child rotation per pair is V = C^T P.
    measured = np.einsum("nji,njk->nik", child_matrices, parent_matrices, optimize=True)
    # All symmetry-equivalent descriptions of every measurement, computed once:
    # S_c (C^T P) S_p over both groups.
    candidates = np.einsum(
        "aij,njk,bkl->nabil", child_operators, measured, parent_operators, optimize=True
    )
    pair_count = measured.shape[0]
    flat_candidates = candidates.reshape(pair_count, -1, 3, 3)
    estimate = nominal.parent_to_child_rotation.as_matrix()
    iterations = 0
    converged = False
    aligned = flat_candidates[:, 0]
    previous_best: np.ndarray | None = None
    while iterations < max_iterations and not converged:
        iterations += 1
        relative = np.einsum("ncij,kj->ncik", flat_candidates, estimate, optimize=True)
        traces = np.trace(relative, axis1=-2, axis2=-1)
        best = np.argmax(traces, axis=1)
        aligned = flat_candidates[np.arange(pair_count), best]
        quaternions = matrices_to_quaternions(aligned)
        scatter = quaternions.T @ quaternions
        eigenvalues, eigenvectors = np.linalg.eigh(scatter)
        mean_quaternion = eigenvectors[:, int(np.argmax(eigenvalues))]
        updated = Rotation(quaternion=mean_quaternion).as_matrix()
        step_angle = _rotation_angles_deg_from_matrices(
            (updated @ estimate.T)[None, :, :]
        )[0]
        estimate = updated
        # Fixed point: alignment assignments stable (the mean is then a
        # deterministic function of them), or the step under the angular
        # tolerance. The assignment test is robust to the ~1e-6 deg
        # matrix->quaternion->matrix round-trip noise floor.
        converged = (
            previous_best is not None and bool(np.array_equal(best, previous_best))
        ) or step_angle <= convergence_tol_deg
        previous_best = best
    residuals = _rotation_angles_deg_from_matrices(
        np.einsum("nij,kj->nik", aligned, estimate, optimize=True)
    )
    fitted_rotation = Rotation.from_matrix(estimate).canonicalized()
    fitted = OrientationRelationship(
        name=f"{nominal.name}_fitted",
        parent_phase=nominal.parent_phase,
        child_phase=nominal.child_phase,
        parent_to_child_rotation=fitted_rotation,
        provenance=provenance or nominal.provenance,
    )
    deviation = _symmetry_reduced_angle_between_deg(
        estimate,
        nominal.parent_to_child_rotation.as_matrix(),
        child_operators=child_operators,
        parent_operators=parent_operators,
    )
    return OrientationRelationshipFitReport(
        relationship=fitted,
        nominal_name=nominal.name,
        residuals_deg=residuals,
        iterations=iterations,
        converged=converged,
        deviation_from_nominal_deg=deviation,
        provenance=provenance or nominal.provenance,
    )


def _integer_index_orbit(
    indices: np.ndarray, *, phase: Phase, reciprocal: bool
) -> np.ndarray:
    """Symmetry-equivalent integer index triples of one plane or direction.

    Operators act on the Cartesian image (reciprocal basis for planes, direct
    basis for directions); the recovered components must be integers for
    crystallographic operators. Antiparallel members collapse to a canonical
    sign (first nonzero component positive), so each returned row names one
    family member once.
    """

    basis = (
        phase.lattice.reciprocal_basis().matrix
        if reciprocal
        else phase.lattice.direct_basis().matrix
    )
    operators = (
        phase.symmetry.operators
        if phase.symmetry is not None
        else np.eye(3, dtype=np.float64)[None, :, :]
    )
    cartesian = basis @ np.asarray(indices, dtype=np.float64)
    images = np.einsum("oij,j->oi", operators, cartesian, optimize=True)
    recovered = np.linalg.solve(basis, images.T).T
    rounded = np.rint(recovered)
    if not np.allclose(recovered, rounded, atol=1e-8):
        raise ValueError(
            "Symmetry orbit did not recover integer indices; the phase symmetry "
            "operators are inconsistent with the lattice."
        )
    members = rounded.astype(np.int64)
    divisors = np.gcd.reduce(np.abs(members), axis=1)
    members = members // np.maximum(divisors, 1)[:, None]
    canonical: dict[tuple[int, int, int], np.ndarray] = {}
    for row in members:
        signed = row.copy()
        nonzero = np.nonzero(signed)[0]
        if signed[nonzero[0]] < 0:
            signed = -signed
        canonical[(int(signed[0]), int(signed[1]), int(signed[2]))] = signed
    orbit = np.stack(list(canonical.values()), axis=0)
    orbit = np.ascontiguousarray(orbit)
    orbit.setflags(write=False)
    return orbit


@dataclass(frozen=True, slots=True)
class ParallelismMatch:
    """One near-parallelism between a parent family member and its child image."""

    variant_index: int
    parent_indices: np.ndarray
    child_indices: np.ndarray
    angular_deviation_deg: float

    def __post_init__(self) -> None:
        parent = np.asarray(self.parent_indices, dtype=np.int64)
        child = np.asarray(self.child_indices, dtype=np.int64)
        if parent.shape != (3,) or child.shape != (3,):
            raise ValueError("ParallelismMatch indices must have shape (3,).")
        if self.variant_index <= 0:
            raise ValueError("ParallelismMatch.variant_index must be positive.")
        if not np.isfinite(self.angular_deviation_deg) or self.angular_deviation_deg < 0.0:
            raise ValueError("angular_deviation_deg must be finite and non-negative.")
        for array in (parent, child):
            array.setflags(write=False)
        object.__setattr__(self, "parent_indices", parent)
        object.__setattr__(self, "child_indices", child)


@dataclass(frozen=True, slots=True)
class ParallelismReport:
    """Near-parallel plane or direction pairs across transformation variants.

    ``matches`` holds every (variant, parent family member, rationalized child
    image) whose angular deviation is within ``tolerance_deg``, sorted by
    variant then deviation. The general machine behind statements like
    "{111} parent is parallel to {011} child within 0 deg in every variant".
    """

    relationship_name: str
    kind: str
    tolerance_deg: float
    matches: tuple[ParallelismMatch, ...]
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.kind not in {"plane", "direction"}:
            raise ValueError("ParallelismReport.kind must be 'plane' or 'direction'.")
        if not np.isfinite(self.tolerance_deg) or self.tolerance_deg < 0.0:
            raise ValueError("tolerance_deg must be finite and non-negative.")
        object.__setattr__(self, "matches", tuple(self.matches))

    def describe(self) -> str:
        formatter = format_plane_indices if self.kind == "plane" else format_direction_indices
        noun = "plane" if self.kind == "plane" else "direction"
        lines = [
            f"Parallel {noun} search under orientation relationship "
            f"'{self.relationship_name}': {len(self.matches)} match(es) within "
            f"{self.tolerance_deg:.3f} deg (angles in degrees, indices in the "
            "respective crystal bases)."
        ]
        for match in self.matches[:12]:
            parent_text = formatter(_index_tuple(match.parent_indices), style="plain")
            child_text = formatter(_index_tuple(match.child_indices), style="plain")
            lines.append(
                f"  variant {match.variant_index}: parent {parent_text} || child {child_text} "
                f"(deviation {match.angular_deviation_deg:.4f} deg)"
            )
        if len(self.matches) > 12:
            lines.append(f"  ... and {len(self.matches) - 12} more.")
        return "\n".join(lines)


def _find_parallels(
    relationship: OrientationRelationship,
    indices: np.ndarray,
    *,
    kind: str,
    tolerance_deg: float,
    include_family: bool,
    variants: tuple[TransformationVariant, ...] | None,
    max_index: int,
    provenance: ProvenanceRecord | None,
) -> ParallelismReport:
    reciprocal = kind == "plane"
    members = (
        _integer_index_orbit(indices, phase=relationship.parent_phase, reciprocal=reciprocal)
        if include_family
        else np.asarray(indices, dtype=np.int64)[None, :]
    )
    resolved = relationship.generate_variants() if variants is None else variants
    parent_phase = relationship.parent_phase
    matches: list[ParallelismMatch] = []
    for variant in resolved:
        for member in members:
            if reciprocal:
                result: PlaneCorrespondence | DirectionCorrespondence = (
                    relationship.map_plane_to_child(
                        CrystalPlane(MillerIndex(member, phase=parent_phase), phase=parent_phase),
                        variant=variant,
                        max_index=max_index,
                    )
                )
            else:
                result = relationship.map_direction_to_child(
                    CrystalDirection(member.astype(np.float64), phase=parent_phase),
                    variant=variant,
                    max_index=max_index,
                )
            if result.angular_residual_deg <= tolerance_deg:
                matches.append(
                    ParallelismMatch(
                        variant_index=variant.variant_index,
                        parent_indices=member,
                        child_indices=result.rational_indices,
                        angular_deviation_deg=result.angular_residual_deg,
                    )
                )
    matches.sort(key=lambda match: (match.variant_index, match.angular_deviation_deg))
    return ParallelismReport(
        relationship_name=relationship.name,
        kind=kind,
        tolerance_deg=tolerance_deg,
        matches=tuple(matches),
        provenance=provenance or relationship.provenance,
    )


def find_parallel_planes(
    relationship: OrientationRelationship,
    parent_plane: CrystalPlane,
    *,
    tolerance_deg: float = 0.5,
    include_family: bool = True,
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    provenance: ProvenanceRecord | None = None,
) -> ParallelismReport:
    """Find child planes (near-)parallel to a parent plane family, per variant.

    Purpose: answers "which child (hkl) is parallel to which member of this
    parent plane family, in which variant, and how exactly" — e.g. under
    Kurdjumov-Sachs each of the 24 variants pairs exactly one {111} parent
    member with a {011} child plane at zero deviation (its close-packed
    plane).

    Inputs: the parent ``CrystalPlane`` (its symmetry family is enumerated
    unless ``include_family=False``), an angular ``tolerance_deg``, and
    optionally a variant tuple and rationalization bound.

    Output: a ``ParallelismReport`` (see its ``describe()``).
    """

    if not phases_semantically_match(parent_plane.phase, relationship.parent_phase):
        raise ValueError("parent_plane.phase must match the relationship parent phase.")
    return _find_parallels(
        relationship,
        parent_plane.miller.indices,
        kind="plane",
        tolerance_deg=tolerance_deg,
        include_family=include_family,
        variants=variants,
        max_index=max_index,
        provenance=provenance,
    )


def find_parallel_directions(
    relationship: OrientationRelationship,
    parent_direction: CrystalDirection,
    *,
    tolerance_deg: float = 0.5,
    include_family: bool = True,
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    provenance: ProvenanceRecord | None = None,
) -> ParallelismReport:
    """Find child directions (near-)parallel to a parent direction family, per variant.

    The direction-space counterpart of ``find_parallel_planes``; antiparallel
    family members collapse to one canonical-sign representative. The parent
    direction must have integer (index-like) coordinates, since the family is
    enumerated as an integer orbit.
    """

    if not phases_semantically_match(parent_direction.phase, relationship.parent_phase):
        raise ValueError("parent_direction.phase must match the relationship parent phase.")
    coordinates = np.asarray(parent_direction.coordinates, dtype=np.float64)
    rounded = np.rint(coordinates)
    if not np.allclose(coordinates, rounded, atol=1e-8):
        raise ValueError(
            "find_parallel_directions requires integer direction coordinates."
        )
    return _find_parallels(
        relationship,
        rounded.astype(np.int64),
        kind="direction",
        tolerance_deg=tolerance_deg,
        include_family=include_family,
        variants=variants,
        max_index=max_index,
        provenance=provenance,
    )


def variant_close_packed_groups(
    relationship: OrientationRelationship,
    parent_plane: CrystalPlane,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
) -> np.ndarray:
    """Group variants by the parent family member they carry into exact parallelism.

    Purpose: the packet classification of martensite crystallography — each
    transformation variant maps exactly one member of the defining parent
    plane family onto its low-index child plane (its close-packed / habit
    packet plane). Variants sharing that member form one packet: for
    Kurdjumov-Sachs and the {111} family this yields the four packets of six
    variants of lath martensite (Morito et al.); for Burgers and the {110}
    family, six groups of two.

    Inputs: the relationship, the defining parent ``CrystalPlane`` (its
    symmetry family is enumerated), and optionally a variant tuple.

    Output: a read-only ``(n_variants,)`` int array of 0-based group labels in
    ``generate_variants()`` order, relabeled by first occurrence.
    """

    if not phases_semantically_match(parent_plane.phase, relationship.parent_phase):
        raise ValueError("parent_plane.phase must match the relationship parent phase.")
    resolved = relationship.generate_variants() if variants is None else variants
    members = _integer_index_orbit(
        parent_plane.miller.indices, phase=relationship.parent_phase, reciprocal=True
    )
    parent_phase = relationship.parent_phase
    member_planes = [
        CrystalPlane(MillerIndex(member, phase=parent_phase), phase=parent_phase)
        for member in members
    ]
    raw_labels = np.empty(len(resolved), dtype=np.int64)
    for index, variant in enumerate(resolved):
        residuals = [
            relationship.map_plane_to_child(
                plane, variant=variant, max_index=max_index
            ).angular_residual_deg
            for plane in member_planes
        ]
        raw_labels[index] = int(np.argmin(residuals))
    relabel: dict[int, int] = {}
    labels = np.empty_like(raw_labels)
    for index, raw in enumerate(raw_labels):
        labels[index] = relabel.setdefault(int(raw), len(relabel))
    labels = np.ascontiguousarray(labels)
    labels.setflags(write=False)
    return labels


@dataclass(frozen=True, slots=True)
class DeformationGradientReport:
    """Lattice-correspondence deformation of a transformation (Bain-strain family).

    ``deformation_gradient`` is the parent-frame map carrying parent lattice
    vectors onto their corresponding child lattice vectors under the integer
    lattice correspondence (nearest-integer matrix of the exact index
    correspondence); ``stretch_tensor`` is its symmetric right-stretch factor
    with ``principal_stretches`` / ``principal_directions`` (parent crystal
    frame) and ``volume_ratio = det F``. ``polar_rotation_deg`` is the angle
    of the residual rotation in the polar decomposition ``F = R_polar U`` —
    zero when the relationship equals the pure correspondence distortion
    (Bain), and the classic rigid-body rotation relative to Bain for KS-class
    relationships. ``correspondence_max_component_error`` is the largest
    entry-wise distance of the exact correspondence from the integer one.
    """

    relationship_name: str
    variant_index: int | None
    deformation_gradient: np.ndarray
    stretch_tensor: np.ndarray
    principal_stretches: np.ndarray
    principal_directions: np.ndarray
    volume_ratio: float
    correspondence: np.ndarray
    polar_rotation_deg: float
    correspondence_max_component_error: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        gradient = np.asarray(self.deformation_gradient, dtype=np.float64)
        stretch = np.asarray(self.stretch_tensor, dtype=np.float64)
        stretches = np.asarray(self.principal_stretches, dtype=np.float64).reshape(-1)
        directions = np.asarray(self.principal_directions, dtype=np.float64)
        correspondence = np.asarray(self.correspondence, dtype=np.int64)
        if gradient.shape != (3, 3) or stretch.shape != (3, 3) or correspondence.shape != (3, 3):
            raise ValueError("Deformation matrices must have shape (3, 3).")
        if stretches.shape != (3,) or directions.shape != (3, 3):
            raise ValueError("Principal quantities must have three entries.")
        if np.any(stretches <= 0.0) or self.volume_ratio <= 0.0:
            raise ValueError("Principal stretches and volume ratio must be positive.")
        if self.polar_rotation_deg < 0.0 or self.correspondence_max_component_error < 0.0:
            raise ValueError("Polar rotation and component error must be non-negative.")
        for array in (gradient, stretch, stretches, directions, correspondence):
            array.setflags(write=False)
        object.__setattr__(self, "deformation_gradient", gradient)
        object.__setattr__(self, "stretch_tensor", stretch)
        object.__setattr__(self, "principal_stretches", stretches)
        object.__setattr__(self, "principal_directions", directions)
        object.__setattr__(self, "correspondence", correspondence)

    def describe(self) -> str:
        """Prose summary: principal strains, volume change, correspondence quality."""

        strains = (self.principal_stretches - 1.0) * 100.0
        variant_text = (
            f" (variant {self.variant_index})" if self.variant_index is not None else ""
        )
        columns = ", ".join(
            format_direction_indices(_index_tuple(self.correspondence[:, i]), style="plain")
            for i in range(3)
        )
        return (
            f"Transformation deformation for '{self.relationship_name}'{variant_text}: "
            f"principal strains {strains[0]:+.2f}%, {strains[1]:+.2f}%, {strains[2]:+.2f}% "
            f"(principal directions in the parent crystal frame), volume change "
            f"{(self.volume_ratio - 1.0) * 100.0:+.2f}%. Integer lattice correspondence maps "
            f"the parent basis to {columns} (child indices; largest entry-wise deviation "
            f"{self.correspondence_max_component_error:.3f}). Residual polar rotation "
            f"{self.polar_rotation_deg:.2f} deg — zero for the pure correspondence "
            "distortion (Bain), and the rigid-body rotation relative to it for "
            "KS-class relationships."
        )


@dataclass(frozen=True, slots=True)
class VariantPoleFigure:
    """Specimen-frame pole positions of a child plane family for every variant.

    ``poles`` holds one specimen-frame unit vector per (variant, family
    member) pair, in variant-major order; ``variant_indices`` gives the
    1-based variant of each row. This is the prediction layer behind variant
    pole figures: overlay it on a measured child pole figure to identify
    which variants are present and how strongly variant selection acts.
    """

    relationship_name: str
    child_plane: CrystalPlane
    poles: VectorSet
    variant_indices: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        indices = np.asarray(self.variant_indices, dtype=np.int64).reshape(-1)
        if indices.shape != (len(self.poles.values),):
            raise ValueError("variant_indices must have one entry per pole row.")
        if indices.size == 0:
            raise ValueError("VariantPoleFigure requires at least one pole.")
        if np.any(indices <= 0):
            raise ValueError("variant_indices must be strictly positive.")
        indices.setflags(write=False)
        object.__setattr__(self, "variant_indices", indices)

    @property
    def variant_count(self) -> int:
        return int(np.unique(self.variant_indices).size)

    def describe(self) -> str:
        """Prose summary: plane family, variants, and pole counts."""

        per_variant = int(
            np.count_nonzero(self.variant_indices == int(self.variant_indices[0]))
        )
        plane_text = format_plane_indices(
            _index_tuple(self.child_plane.miller.indices), style="plain"
        )
        return (
            f"Variant pole figure for orientation relationship "
            f"'{self.relationship_name}': specimen-frame poles of the child "
            f"{plane_text} family for {self.variant_count} variant(s), "
            f"{per_variant} family member(s) per variant "
            f"({self.variant_indices.size} poles total, unit vectors in the "
            "specimen frame; overlay on a measured child pole figure to "
            "identify operative variants)."
        )


def variant_pole_figure(
    parent_orientation: Orientation,
    relationship: OrientationRelationship,
    child_plane: CrystalPlane,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> VariantPoleFigure:
    """Predict the child-plane pole figure of every transformation variant.

    Purpose: given a parent orientation, computes where each variant's child
    plane family lands on the specimen sphere — the standard overlay for
    reading measured product-phase pole figures (which variants formed, and
    whether selection favors some). Uses the canonical composition
    ``g_child = g_parent o V^T`` and the child plane's full symmetry family
    (antipodal-collapsed integer orbit).

    Inputs: the parent ``Orientation`` (phase must match the relationship
    parent phase), the relationship, the child ``CrystalPlane`` whose family
    is plotted, and optionally a variant tuple.

    Output: a ``VariantPoleFigure`` (see its ``describe()``); plot it with
    ``pytex.plot_variant_pole_figure``.
    """

    if not phases_semantically_match(parent_orientation.phase, relationship.parent_phase):
        raise ValueError("parent_orientation.phase must match the relationship parent phase.")
    if not phases_semantically_match(child_plane.phase, relationship.child_phase):
        raise ValueError("child_plane.phase must match the relationship child phase.")
    resolved = relationship.generate_variants() if variants is None else variants
    child_phase = relationship.child_phase
    members = _integer_index_orbit(
        child_plane.miller.indices, phase=child_phase, reciprocal=True
    )
    reciprocal_basis = child_phase.lattice.reciprocal_basis().matrix
    normals = members.astype(np.float64) @ reciprocal_basis.T
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    parent_matrix = parent_orientation.rotation.as_matrix()
    pole_rows: list[np.ndarray] = []
    index_rows: list[int] = []
    for variant in resolved:
        child_matrix = parent_matrix @ variant.parent_to_child_rotation.as_matrix().T
        pole_rows.append(normals @ child_matrix.T)
        index_rows.extend([variant.variant_index] * normals.shape[0])
    poles = VectorSet(
        values=np.concatenate(pole_rows, axis=0),
        reference_frame=parent_orientation.specimen_frame,
        provenance=provenance or relationship.provenance,
    )
    return VariantPoleFigure(
        relationship_name=relationship.name,
        child_plane=child_plane,
        poles=poles,
        variant_indices=np.asarray(index_rows, dtype=np.int64),
        provenance=provenance or relationship.provenance,
    )


def map_direction_across_variants(
    relationship: OrientationRelationship,
    direction: CrystalDirection,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
) -> tuple[DirectionCorrespondence, ...]:
    """Map one parent direction through every transformation variant.

    Purpose: the variant-resolved answer to "what does this parent ``[uvw]``
    become in the product phase" — one ``DirectionCorrespondence`` per variant
    (in ``generate_variants()`` order unless ``variants`` is given), each
    carrying exact components, rationalized indices, and the residual.
    """

    resolved = relationship.generate_variants() if variants is None else variants
    return tuple(
        relationship.map_direction_to_child(direction, variant=variant, max_index=max_index)
        for variant in resolved
    )


def map_plane_across_variants(
    relationship: OrientationRelationship,
    plane: CrystalPlane,
    *,
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
) -> tuple[PlaneCorrespondence, ...]:
    """Map one parent plane through every transformation variant.

    Purpose: the variant-resolved answer to "what does this parent ``(hkl)``
    become in the product phase" — one ``PlaneCorrespondence`` per variant (in
    ``generate_variants()`` order unless ``variants`` is given), each carrying
    exact components, rationalized Miller indices, and the residual.
    """

    resolved = relationship.generate_variants() if variants is None else variants
    return tuple(
        relationship.map_plane_to_child(plane, variant=variant, max_index=max_index)
        for variant in resolved
    )


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
        child_count = len(self.child_orientations)
        if self.variant_indices is None:
            base = self.orientation_relationship.parent_to_child_rotation.as_matrix()
            variant_matrices = np.repeat(base[None, :, :], child_count, axis=0)
        else:
            variant_lookup = {
                variant.variant_index: variant.parent_to_child_rotation.as_matrix()
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
            unique_indices, inverse = np.unique(variant_indices, return_inverse=True)
            unique_matrices = np.stack(
                [variant_lookup[int(index)] for index in unique_indices], axis=0
            )
            variant_matrices = unique_matrices[inverse]
        parent_matrix = self.parent_orientation.rotation.as_matrix()
        # Canonical crystal->specimen orientation convention: the child
        # orientation is g_child = g_parent o V^T, i.e. C = P @ V^T.
        predicted = np.einsum("ij,nkj->nik", parent_matrix, variant_matrices, optimize=True)
        return OrientationSet.from_matrices(
            predicted,
            crystal_frame=self.child_orientations.crystal_frame,
            specimen_frame=self.child_orientations.specimen_frame,
            symmetry=self.child_orientations.symmetry,
            phase=self.child_orientations.phase,
            provenance=self.provenance or self.parent_orientation.provenance,
        )
