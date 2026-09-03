from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike

if TYPE_CHECKING:  # pragma: no cover - import-cycle guard
    # parent_reconstruction imports this module, so the catalog type is
    # available for annotations only; the runtime import is local to the one
    # function that needs it.
    from pytex.core.parent_reconstruction import OrientationRelationshipCatalog

from pytex.core._angles import acute_angle_between_unit_vectors_rad
from pytex.core._arrays import as_int_array
from pytex.core.batches import VectorSet
from pytex.core.frames import ReferenceFrame
from pytex.core.hexagonal import (
    direction_uvw_to_uvtw,
    is_hexagonal_phase,
    plane_hkl_to_hkil,
)
from pytex.core.lattice import (
    CrystalDirection,
    CrystalPlane,
    MillerIndex,
    Phase,
    phases_semantically_match,
)
from pytex.core.notation import (
    format_direction_indices,
    format_plane_family_indices,
    format_plane_indices,
)
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


#: Denominators searched when rationalizing a lattice correspondence, in order.
#:
#: 1 covers the cubic-to-cubic relationships, whose conventional basis vectors
#: map onto child lattice vectors. 2 is needed for reconstructive bcc-to-hcp
#: (Burgers), where the bcc primitive cell is half the hcp one and an atomic
#: shuffle carries the missing half — the correspondence is genuinely
#: half-integer. The list is short on purpose: allowing any denominator would
#: fit numerical noise instead of crystallography.
_CORRESPONDENCE_DENOMINATORS: tuple[int, ...] = (1, 2)

#: How close a correspondence determinant must be to a whole number to count.
#:
#: A lattice correspondence maps one cell of parent lattice points onto a whole
#: number of child cells, so its determinant is an integer — 2 for the
#: fcc-to-bcc family (one fcc cell of four atoms onto two bcc cells), 1 for
#: Burgers in either direction. A candidate whose determinant is fractional is
#: not a lattice correspondence at all, whatever its entry-wise fit, and this is
#: what separates the right denominator from a plausible-looking rounding.
_CORRESPONDENCE_DETERMINANT_TOLERANCE = 1e-6

#: How much better a finer denominator must fit before it is preferred.
#:
#: An integer correspondence is the physically expected case, so a coarser grid
#: wins ties and near-ties. A denominator-2 grid can differ from a
#: denominator-1 grid by at most 0.25 in an entry, so requiring a 0.1 improvement
#: says "only accept the shuffle if it explains a substantial part of the
#: misfit". Without this, a large lattice contraction can make the finer grid fit
#: marginally better with a physically wrong determinant — for a 3.60 to 2.87
#: Bain pair, by 0.009, at determinant 3 instead of 2.
_CORRESPONDENCE_REFINEMENT_MARGIN = 0.1


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
        family = {"432": "cubic", "622": "hexagonal", "222": "orthorhombic"}.get(
            expected, expected
        )
        raise ValueError(
            f"{relationship} correspondence requires a {family} {role} phase "
            f"with proper point group {expected}."
        )


def _require_cubic_phase_for_bain(phase: Phase, *, role: str) -> None:
    _require_proper_point_group(phase, "432", role=role, relationship="Bain")


@dataclass(frozen=True, slots=True)
class OrientationRelationship:
    """The fixed crystallographic relationship between a parent and a child phase.

    Purpose
    -------
    The flagship object of the transformation subsystem. A phase
    transformation does not produce arbitrary child orientations: the product
    lattice inherits a definite relationship to the parent, stated in the
    literature as a pair of parallelisms — Kurdjumov-Sachs as
    ``(111)_gamma || (011)_alpha`` with ``[-101]_gamma || [-1-11]_alpha``.

    This type holds both the rotation that statement implies *and* the
    statement itself, so reports can restate the relationship in the terms it
    was defined by rather than only as an angle.

    Attributes
    ----------
    name : str
        Non-empty identifier, carried into every derived report.
    parent_phase, child_phase : Phase
        Must be semantically distinct; a relationship between a phase and
        itself is rejected.
    parent_to_child_rotation : Rotation
        The rotation taking parent crystal axes to child crystal axes.
    parallel_directions : tuple of (CrystalDirection, CrystalDirection)
        The direction parallelisms defining the relationship, as typed pairs
        so phase and basis meaning are not lost.
    parallel_planes : tuple of (CrystalPlane, CrystalPlane)
        The plane parallelisms.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Because both phases carry symmetry, one relationship generates a family
    of symmetry-equivalent :class:`TransformationVariant` realizations — 24
    for Kurdjumov-Sachs, 12 for Nishiyama-Wassermann.
    """

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
        """Crystal-domain reference frame of the parent phase.
        """

        return self.parent_phase.crystal_frame

    @property
    def child_crystal_frame(self) -> ReferenceFrame:
        """Crystal-domain reference frame of the child (product) phase.
        """

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
        """Build an orientation relationship from the parallelisms that define it.

        Purpose
        -------
        The literature states an OR as a pair of parallelisms — Kurdjumov-Sachs
        as ``(111)_gamma || (011)_alpha`` with ``[-101]_gamma || [-1-11]_alpha``.
        This constructor turns that statement directly into the rotation it
        implies, so a published OR can be reproduced without composing rotations
        by hand.

        Parameters
        ----------
        name : str
            Identifier carried into reports and figures.
        parent_plane, child_plane : CrystalPlane
            The planes made parallel; the child normal is brought onto the
            parent normal.
        parent_direction, child_direction : CrystalDirection
            The in-plane directions made parallel, fixing the remaining
            rotation about the common normal. They need not be exactly
            perpendicular to their plane normals; the normal component is
            removed, so a slightly inconsistent literature statement still
            yields a proper rotation.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        OrientationRelationship
            Carrying both phases, the rotation, and the defining parallelisms —
            which are kept so that reports can restate the OR in the terms it
            was defined by.
        """

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
        """The Bain orientation relationship between two cubic phases.

        Purpose
        -------
        The Bain correspondence is the minimum-strain lattice correspondence for
        the FCC-to-BCC transformation: ``{100}`` planes and ``<100>`` directions
        of the two lattices are held parallel. It is the reference against which
        Kurdjumov-Sachs and Nishiyama-Wassermann are usually described, and the
        starting point for phenomenological martensite theory.

        Both phases must be cubic; this is checked rather than assumed.
        """

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
    def from_cube_on_cube_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "cube_on_cube",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Cube-on-cube OR: (001)_p || (001)_c, [100]_p || [100]_c.

        Purpose
        -------
        The parallel-axes relationship: the two cubic lattices share all three
        crystal axes, so the parent-to-child rotation is the identity. It is the
        relationship of a coherent cubic precipitate grown on a cubic matrix --
        gamma-prime in a nickel superalloy, TiN or NbC in ferrite, an epitaxial
        cubic film on a cubic substrate -- and the baseline every other cubic
        pairing is a departure from.

        It is also the degenerate case worth having explicitly: one variant, no
        misorientation, and a composite diffraction pattern in which every child
        reflection lies on a parent reciprocal-lattice row. A superlattice
        reflection in such a pattern therefore comes from the child's *motif*,
        not from any rotation, which is exactly the reading the identity
        rotation makes available.

        Both phases must be cubic (proper point group 432); this is checked
        rather than assumed, because "cube on cube" says nothing about two
        lattices that have no cube to share.
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Cube-on-cube"
        )
        _require_proper_point_group(child_phase, "432", role="child", relationship="Cube-on-cube")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(
                _miller_index((0, 0, 1), phase=parent_phase), phase=parent_phase
            ),
            child_plane=CrystalPlane(
                _miller_index((0, 0, 1), phase=child_phase), phase=child_phase
            ),
            parent_direction=_crystal_direction((1.0, 0.0, 0.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 0.0, 0.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_fcc_twin_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "fcc_twin",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Coherent FCC twin: (111)_p || (111)_c, [1-10]_p || [-110]_c.

        Purpose
        -------
        The annealing- and deformation-twin relationship of face-centred cubic
        metals, stated as the parallelism it is defined by rather than as the
        angle it comes out as. The statement is the 180 degree rotation about
        the ``<111>`` twin-plane normal -- equivalently, the mirror in the
        composition plane followed by inversion, which is the proper rotation a
        centrosymmetric lattice cannot tell from the mirror.

        Under cubic symmetry that 180 degree rotation reduces to the familiar
        **60 degrees about <111>**, the Sigma 3 coincidence-site boundary, and
        :meth:`misorientation` reports it that way. Four ``{111}`` planes give
        **four variants**, one twin orientation per close-packed plane of the
        parent.

        Why parent and child are two phases
        -----------------------------------
        A twin is not a second material: matrix and twin have the same lattice
        and the same point group. This type requires distinct phases -- an OR
        between a phase and itself is rejected at construction -- so a twin is
        expressed as two phase *objects* differing in name alone, e.g.
        ``Phase("nickel", ...)`` and ``Phase("nickel (twin)", ...)``. That is
        not a workaround: the twin is a distinct orientation domain, it is what
        an EBSD map segments it as, and naming it separately is what lets every
        variant, boundary and composite-diffraction report say which domain a
        quantity belongs to.

        Both phases must be cubic (proper point group 432). The relationship is
        stated for face-centred cubic because that is the structure whose
        stacking makes ``{111}`` the coherent twin plane; the rotation itself is
        a property of the cubic point group, so it applies unchanged to any
        cubic pair whose twin plane is ``{111}``.
        """

        _require_proper_point_group(parent_phase, "432", role="parent", relationship="FCC twin")
        _require_proper_point_group(child_phase, "432", role="child", relationship="FCC twin")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(
                _miller_index((1, 1, 1), phase=parent_phase), phase=parent_phase
            ),
            child_plane=CrystalPlane(
                _miller_index((1, 1, 1), phase=child_phase), phase=child_phase
            ),
            parent_direction=_crystal_direction((1.0, -1.0, 0.0), phase=parent_phase),
            child_direction=_crystal_direction((-1.0, 1.0, 0.0), phase=child_phase),
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
        """The Nishiyama-Wassermann orientation relationship between cubic phases.

        Purpose
        -------
        The NW relationship holds ``{111}_fcc || {110}_bcc`` with
        ``<-1-12>_fcc || <-110>_bcc``, and generates 12 variants where
        Kurdjumov-Sachs generates 24. NW and KS differ by only about 5.26
        degrees, which is why distinguishing them from measured data needs the
        residual statistics that
        :func:`characterize_orientation_relationship` reports.

        Both phases must be cubic; this is checked rather than assumed.
        """

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

    @classmethod
    def from_potter_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "potter",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Potter OR: {10-11}_hcp || {110}_bcc, <11-20>_hcp || <111>_bcc.

        The hcp->bcc precipitate/matrix relationship of Potter (V2N/V3N
        nitrides in alpha-vanadium; J. Less-Common Metals 31 (1973) 299).
        The close-packed directions coincide exactly as in the Burgers
        relationship, but the exact plane parallelism is carried by the
        pyramidal {10-11} plane against a cubic {110} plane; the basal plane
        is left a small, c/a-dependent rotation (~2 deg for typical metal
        ratios) from its Burgers {110} partner about the shared close-packed
        direction. The parent must be hexagonal (proper group 622) and the
        child cubic (proper group 432).
        """

        _require_proper_point_group(parent_phase, "622", role="parent", relationship="Potter")
        _require_proper_point_group(child_phase, "432", role="child", relationship="Potter")
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane.from_miller_bravais((0, 1, -1, 1), phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((1, 1, 0), phase=child_phase),
                                     phase=child_phase),
            parent_direction=CrystalDirection.from_miller_bravais(
                (2, -1, -1, 0), phase=parent_phase
            ),
            child_direction=_crystal_direction((1.0, -1.0, 1.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_bagaryatsky_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "bagaryatsky",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Bagaryatsky OR: (0-11)_bcc || (001)_cem, [1-1-1]_bcc || [100]_cem.

        The classical ferrite->cementite tempering/pearlite relationship
        (Bagaryatskii 1950), stated in the Pnma cementite setting
        (b > a > c; Lipson-Petch axes): [100]_theta || [1-1-1]_alpha,
        [010]_theta || [211]_alpha, [001]_theta || [0-11]_alpha. The parent
        must be cubic (proper group 432) and the child orthorhombic (proper
        group 222). Precise measurements suggest observed orientations are
        actually Isaichev (see ``from_isaichev_correspondence``), a ~3.8 deg
        rotation about the cementite a-axis.
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Bagaryatsky"
        )
        _require_proper_point_group(
            child_phase, "222", role="child", relationship="Bagaryatsky"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((0, -1, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((0, 0, 1), phase=child_phase),
                                     phase=child_phase),
            parent_direction=_crystal_direction((1.0, -1.0, -1.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 0.0, 0.0), phase=child_phase),
            provenance=provenance,
        )

    @classmethod
    def from_isaichev_correspondence(
        cls,
        *,
        parent_phase: Phase,
        child_phase: Phase,
        name: str = "isaichev",
        provenance: ProvenanceRecord | None = None,
    ) -> OrientationRelationship:
        """Isaichev OR: (101)_bcc || (031)_cem, [1-1-1]_bcc || [100]_cem.

        The ferrite->cementite relationship of Isaichev (1947), stated in the
        Pnma cementite setting (b > a > c): it shares the Bagaryatsky
        close-packed direction pairing [100]_theta || [1-1-1]_alpha but pins
        the (031)_theta || (101)_alpha plane parallelism instead, a rotation
        of Bagaryatsky about the cementite a-axis whose magnitude depends on
        the cementite axial ratios (~3.6-3.8 deg for literature lattice
        parameters). Precise diffraction work identifies Isaichev, not
        Bagaryatsky, on tempered martensite. The parent must be cubic (proper
        group 432) and the child orthorhombic (proper group 222).
        """

        _require_proper_point_group(
            parent_phase, "432", role="parent", relationship="Isaichev"
        )
        _require_proper_point_group(
            child_phase, "222", role="child", relationship="Isaichev"
        )
        return cls.from_parallel_plane_direction(
            name=name,
            parent_plane=CrystalPlane(_miller_index((1, 0, 1), phase=parent_phase),
                                      phase=parent_phase),
            child_plane=CrystalPlane(_miller_index((0, 3, 1), phase=child_phase),
                                     phase=child_phase),
            parent_direction=_crystal_direction((1.0, -1.0, -1.0), phase=parent_phase),
            child_direction=_crystal_direction((1.0, 0.0, 0.0), phase=child_phase),
            provenance=provenance,
        )

    def map_parent_vector_to_child(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Map a parent-crystal vector into the child crystal frame.

        Applies the parent-to-child rotation. A ``VectorSet`` must carry the
        parent crystal frame — checked, not assumed — and is returned re-framed
        to the child crystal frame.

        This maps *Cartesian directions*. To answer the index question — which
        child ``(hkl)`` or ``[uvw]`` corresponds to a given parent one — use
        the correspondence surfaces, which additionally rationalize to
        low-integer indices and report the angular residual.
        """

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
        """Map a child-crystal vector into the parent crystal frame.

        The inverse of :meth:`map_parent_vector_to_child`, with the same frame
        checking and the same Cartesian-direction caveat.
        """

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
        correspondence: ArrayLike | None = None,
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

        Reconstructive transformations and the shuffle
            An integer correspondence exists only when the parent's conventional
            basis vectors map onto child *lattice* vectors, which holds for the
            cubic-to-cubic relationships. It fails for bcc to hcp (Burgers): the
            bcc primitive cell has half the volume of the hcp one, so two bcc
            lattice points map onto one hcp lattice point plus one motif atom.
            That missing half is the Burgers **shuffle**, and the correspondence
            is a **half-integer** matrix rather than an integer one.

            The search therefore widens the denominator until it finds a
            non-singular correspondence, and reports which one it used through
            `DeformationGradientReport.correspondence_denominator`: 1 for the
            cubic cases, 2 when a shuffle carries half a cell. Requiring an
            integer would have refused the Burgers case; allowing any rational
            would fit noise, so the denominator is bounded and reported.

            ``correspondence`` may still be supplied explicitly to pin a
            literature lattice correspondence rather than accept the search.

        Inputs: optionally a variant, and optionally an explicit ``(3, 3)``
        lattice ``correspondence``.

        Output: a ``DeformationGradientReport`` (see its ``describe()``).

        Raises
        ------
        ValueError
            If no correspondence is supplied and no correspondence with a
            denominator within the bound is non-singular.
        """

        rotation, variant_index = self._resolve_variant_rotation(variant)
        parent_basis = self.parent_phase.lattice.direct_basis().matrix
        child_basis = self.child_phase.lattice.direct_basis().matrix
        exact = np.linalg.solve(child_basis, rotation @ parent_basis)
        denominator = 1
        correspondence_matrix: np.ndarray | None = None
        if correspondence is not None:
            supplied = np.asarray(correspondence, dtype=np.float64)
            if supplied.shape != (3, 3):
                raise ValueError("correspondence must have shape (3, 3).")
            if np.isclose(float(np.linalg.det(supplied)), 0.0, atol=1e-12):
                raise ValueError("The supplied lattice correspondence is singular.")
            correspondence_matrix = supplied
            component_error = float(np.max(np.abs(exact - supplied)))
            is_integer = np.allclose(supplied, np.rint(supplied), atol=1e-9)
            doubled = supplied * 2.0
            is_half_integer = np.allclose(doubled, np.rint(doubled), atol=1e-9)
            if is_half_integer and not is_integer:
                denominator = 2
        else:
            # Magnitudes matter for strain, so ray-based rationalization would be
            # wrong here: the correspondence is the nearest matrix over a bounded
            # denominator. Denominator 1 covers the cubic cases; 2 is needed when
            # a shuffle carries half a cell, as in Burgers bcc-to-hcp.
            # Three filters, because no single one is sufficient:
            #
            #   * the determinant must be a non-zero whole number — a lattice
            #     correspondence maps a cell onto a whole number of cells, so a
            #     candidate with det 1.5 is not one, however well its entries fit;
            #   * a finer denominator must fit substantially better, not merely
            #     better, or a large lattice contraction lets it win by a hair
            #     with a wrong determinant;
            #   * subject to that, take the coarsest grid, because an integer
            #     correspondence is the physically expected case.
            #
            # The reverse hcp-to-bcc relationship needs all three: its
            # denominator-1 rounding is invertible with integer determinant 2,
            # but fits so badly that using it reports a doubled cell and a
            # nonsensical +96% volume change.
            candidates: list[tuple[float, int, np.ndarray]] = []
            for candidate_denominator in _CORRESPONDENCE_DENOMINATORS:
                scaled = np.rint(exact * candidate_denominator) / candidate_denominator
                determinant = float(np.linalg.det(scaled))
                whole = round(determinant)
                if whole == 0:
                    continue
                if abs(determinant - whole) > _CORRESPONDENCE_DETERMINANT_TOLERANCE:
                    continue
                candidates.append(
                    (float(np.max(np.abs(exact - scaled))), candidate_denominator, scaled)
                )
            if candidates:
                best_deviation = min(item[0] for item in candidates)
                acceptable = [
                    item
                    for item in candidates
                    if item[0] <= best_deviation + _CORRESPONDENCE_REFINEMENT_MARGIN
                ]
                _, denominator, correspondence_matrix = min(
                    acceptable, key=lambda item: item[1]
                )
            if correspondence_matrix is None:
                raise ValueError(
                    "No non-singular lattice correspondence was found with a denominator "
                    f"in {_CORRESPONDENCE_DENOMINATORS}. Pass an explicit "
                    "`correspondence=` matrix (child-basis coordinates of the images of "
                    "the parent basis vectors) taken from the literature."
                )
            component_error = float(np.max(np.abs(exact - correspondence_matrix)))
        if correspondence_matrix is None:  # pragma: no cover - guarded above
            raise ValueError("Failed to resolve a lattice correspondence.")
        correspondence = correspondence_matrix
        gradient = (
            rotation.T
            @ child_basis
            @ np.asarray(correspondence, dtype=np.float64)
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
            correspondence_denominator=denominator,
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
        """The relationship with parent and child roles exchanged.

        Purpose
        -------
        Turns a parent-to-child OR into the child-to-parent one, with the
        rotation inverted and every parallelism pair reversed. Use it to reason
        from the product phase back to the parent — the direction
        parent-grain reconstruction works in.

        Parameters
        ----------
        name : str, optional
            Name for the inverted relationship; defaults to
            ``"<child>_to_<parent>"``.
        """

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


def _symmetry_image_plane(plane: CrystalPlane, operator: np.ndarray) -> CrystalPlane:
    """The image of ``plane`` under a crystal-frame symmetry operator.

    The operator acts on Cartesian crystal-frame vectors, while a plane is
    stored as reciprocal-basis Miller indices, so the image is obtained by
    conjugating with the reciprocal basis: ``B*^-1 S B* (hkl)``. A point-group
    operator maps the lattice onto itself, so the result is integral; it is
    rounded, and a non-integral result means the operator did not belong to
    the phase's symmetry and is rejected rather than silently rationalized.
    """

    basis = plane.phase.lattice.reciprocal_basis().matrix
    indices = np.linalg.solve(basis, operator @ (basis @ plane.miller.indices.astype(np.float64)))
    rounded = np.rint(indices)
    if not np.allclose(indices, rounded, atol=1e-6):
        raise ValueError(
            "Symmetry operator does not map the lattice onto itself: "
            f"({plane.miller.indices.tolist()}) maps to non-integral {indices.tolist()}."
        )
    return CrystalPlane(MillerIndex(rounded, phase=plane.phase), phase=plane.phase)


def _symmetry_image_direction(
    direction: CrystalDirection, operator: np.ndarray
) -> CrystalDirection:
    """The image of ``direction`` under a crystal-frame symmetry operator.

    The direct-basis counterpart of :func:`_symmetry_image_plane`:
    ``B^-1 S B [uvw]``. Directions may legitimately carry non-integral
    coordinates, so an integral result is snapped to exact integers and a
    non-integral one is returned as it stands.
    """

    basis = direction.phase.lattice.direct_basis().matrix
    coordinates = np.asarray(direction.coordinates, dtype=np.float64)
    image = np.linalg.solve(basis, operator @ (basis @ coordinates))
    rounded = np.rint(image)
    if np.allclose(image, rounded, atol=1e-8):
        image = rounded
    return CrystalDirection(image, phase=direction.phase)


@dataclass(frozen=True, slots=True)
class TransformationVariant:
    """One symmetry-equivalent realization of an orientation relationship.

    Purpose
    -------
    A parent grain transforming under a given relationship can produce
    several crystallographically equivalent child orientations, one per
    combination of parent and child symmetry operators. These variants are
    what a martensitic or bainitic microstructure is built from, and which
    variants actually appear — variant selection — is the observable
    signature of the transformation mechanism.

    Attributes
    ----------
    orientation_relationship : OrientationRelationship
        The relationship this variant realizes.
    variant_index : int
        One-based index; strictly positive.
    parent_operator_index, child_operator_index : int
        Which symmetry operators generated this variant. Note that the
        numbering is an enumeration-order artefact and does not by itself
        match a published variant table.
    parent_to_child_rotation : Rotation
        This variant's specific rotation.
    habit_plane_pairs : tuple of (CrystalPlane, CrystalPlane)
        Descriptive only. Habit planes are *not* computed here; populating
        this slot requires invariant-line or phenomenological martensite
        theory, which the library does not yet implement.
    provenance : ProvenanceRecord, optional
    """

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

    @property
    def parent_symmetry_operator(self) -> np.ndarray:
        """The parent point-group operator ``S_p`` that generated this variant."""

        operators = _parent_operators(self.orientation_relationship)
        return np.asarray(operators[self.parent_operator_index], dtype=np.float64)

    @property
    def child_symmetry_operator(self) -> np.ndarray:
        """The child point-group operator ``S_c`` that generated this variant."""

        operators = _child_operators(self.orientation_relationship)
        return np.asarray(operators[self.child_operator_index], dtype=np.float64)

    @property
    def parallel_planes(self) -> tuple[tuple[CrystalPlane, CrystalPlane], ...]:
        """This variant's own plane parallelisms, not the relationship's.

        Purpose
        -------
        A variant is generated as ``V = S_c R S_p^T``, so the parallelism it
        actually realizes is not the nominal pair the relationship was defined
        by: it is that pair carried by the generating symmetry operators,
        ``(S_p n_parent) || (S_c n_child)``. Substituting the nominal pair
        produces a figure or a report that looks right and is wrong for every
        variant but the first, which is why this is a property of the variant
        rather than something a caller re-derives.

        Returns
        -------
        tuple of (CrystalPlane, CrystalPlane)
            Parent- and child-phase planes, one pair per defining parallelism,
            in the same order as ``OrientationRelationship.parallel_planes``.
            ``V`` maps each parent normal exactly onto its child normal.

        See Also
        --------
        parallel_directions : the direction-space counterpart.
        """

        parent_operator = self.parent_symmetry_operator
        child_operator = self.child_symmetry_operator
        return tuple(
            (
                _symmetry_image_plane(parent_plane, parent_operator),
                _symmetry_image_plane(child_plane, child_operator),
            )
            for parent_plane, child_plane in self.orientation_relationship.parallel_planes
        )

    @property
    def parallel_directions(self) -> tuple[tuple[CrystalDirection, CrystalDirection], ...]:
        """This variant's own direction parallelisms; see :attr:`parallel_planes`."""

        parent_operator = self.parent_symmetry_operator
        child_operator = self.child_symmetry_operator
        return tuple(
            (
                _symmetry_image_direction(parent_direction, parent_operator),
                _symmetry_image_direction(child_direction, child_operator),
            )
            for parent_direction, child_direction in (
                self.orientation_relationship.parallel_directions
            )
        )

    def map_parent_vector_to_child(self, vector: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Map a parent-crystal vector into this variant's child crystal frame.

        Each transformation variant is a distinct symmetry-related realization
        of the orientation relationship, so the same parent direction maps to a
        different child direction in each. A ``VectorSet`` must carry the parent
        crystal frame and is returned re-framed to the child crystal frame.
        """

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


def _parent_operators(relationship: OrientationRelationship) -> np.ndarray:
    symmetry = relationship.parent_phase.symmetry
    if symmetry is None:
        return np.eye(3, dtype=np.float64)[None, :, :]
    return np.asarray(symmetry.operators, dtype=np.float64)


#: Edge block size for :func:`boundary_fingerprint_distances_deg`. The distance
#: kernel is a ``(block, 9) @ (9, K)`` GEMM, so the transient is
#: ``block * K * 8`` bytes; 512 keeps that near 40 MB for a cubic-cubic
#: fingerprint while staying large enough to saturate BLAS.
_FINGERPRINT_BLOCK_SIZE = 512


def intervariant_boundary_fingerprint(
    relationship: OrientationRelationship,
) -> np.ndarray:
    r"""Every misorientation two child grains of one parent can exhibit.

    Purpose: the exact admissible set for same-parent child-child boundaries —
    the quantity an edge test must compare against when deciding whether two
    neighbouring product grains descend from a common parent.

    Theory: a child of parent :math:`P` through variant :math:`V_i` is
    :math:`C_i = P V_i^{\mathsf{T}}`, so the crystal-frame boundary
    misorientation of two same-parent children is
    :math:`C_i^{\mathsf{T}} C_j = V_i V_j^{\mathsf{T}}`. Writing
    :math:`V_i = R S_{p,i}` with :math:`R` the parent-to-child rotation and
    :math:`S_p` the parent point group, the set collapses to
    :math:`R G_p R^{\mathsf{T}}` — the parent group *conjugated* by the OR
    rotation. Each child orientation is itself only defined up to its own
    crystal symmetry, so the observable set is the double coset
    :math:`G_c \left(R G_p R^{\mathsf{T}}\right) G_c`, returned here
    deduplicated (``q`` and ``-q`` identified).

    Inputs: the orientation relationship; its parent and child phase symmetry
    supply :math:`G_p` and :math:`G_c` (a phase without symmetry contributes
    the identity only).

    Output: a read-only ``(k, 3, 3)`` array of rotation matrices. Pair it with
    :func:`boundary_fingerprint_distances_deg` to score measured boundaries.

    Note that this is strictly stronger than comparing misorientation *angles*
    against :func:`intervariant_misorientation_angles_deg`: the angle spectrum
    discards the axis, and for a cubic-cubic relationship an angle-only test
    admits a large fraction of entirely unrelated boundaries.
    """

    rotation = relationship.parent_to_child_rotation.as_matrix()
    parent_ops = _parent_operators(relationship)
    child_ops = _child_operators(relationship)
    conjugated = np.einsum("ij,pjk,lk->pil", rotation, parent_ops, rotation, optimize=True)
    products = np.einsum(
        "aij,pjk,bkl->apbil", child_ops, conjugated, child_ops, optimize=True
    ).reshape(-1, 3, 3)
    fingerprint = np.ascontiguousarray(products[_deduplicate_rotations(products)])
    fingerprint.setflags(write=False)
    return fingerprint


#: How close two rotation matrices must be to count as the same group element.
#:
#: The duplicates being merged are the *same* element reached by different
#: operator products, so they agree to the floating-point floor of three chained
#: matrix multiplications (~1e-15). Distinct elements of these groups are
#: separated by ~1e-1 in the same measure, so the tolerance sits in a gap six
#: orders of magnitude wide and its exact value is not delicate.
_ROTATION_DEDUPLICATION_TOLERANCE = 1e-9


def _deduplicate_rotations(matrices: np.ndarray) -> np.ndarray:
    """Indices of one representative per distinct rotation in a set of matrices.

    Deduplicating on **matrices** rather than on quaternions is the point. A
    quaternion needs a sign convention, usually "make the largest-magnitude
    component positive" — and when two components tie in magnitude, which is
    common for the 90 and 180 degree elements of a crystal point group,
    ``argmax`` breaks the tie arbitrarily. Two numerically identical rotations
    then canonicalize to ``q`` and ``-q``, land far apart, and are counted
    twice. Measured on the Kurdjumov-Sachs double coset: 10665 elements by the
    quaternion route against 10584 genuinely distinct rotations, so 81 elements
    were duplicates.

    Rounding to a fixed number of decimals has a second, independent problem: a
    value near a rounding boundary rounds either way depending on floating-point
    noise, so the count of a mathematically fixed set moved with the lattice
    parameters that entered the rotation (10664 for one cubic pair, 10665 for
    another) even though group theory cannot depend on them.

    Sorting the flattened matrices lexicographically makes duplicates adjacent —
    they differ only in the last bits, so nothing distinct can sort between them
    — and one linear pass merges runs whose successive steps stay inside the
    tolerance. Exact, stable, and ``O(n log n)``.

    Duplicates were never a *correctness* problem downstream, because
    :func:`boundary_fingerprint_distances_deg` takes a maximum over the set. They
    mattered because the set's size is quoted as a scientific quantity.
    """

    if matrices.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    from scipy.spatial import cKDTree

    flat = np.ascontiguousarray(matrices.reshape(matrices.shape[0], -1))
    # A spatial query rather than a sort: duplicates are not reliably adjacent
    # in lexicographic order, because a distinct element can agree with them in
    # the leading entries and sort between them. Neighbourhood membership is the
    # actual predicate, so ask for it directly.
    tree = cKDTree(flat)
    keep = np.ones(flat.shape[0], dtype=bool)
    for left, right in tree.query_pairs(_ROTATION_DEDUPLICATION_TOLERANCE):
        # Keep the earliest index of each duplicate run, so the retained order
        # is the enumeration order and the result is reproducible.
        keep[max(left, right)] = False
    return np.asarray(np.flatnonzero(keep), dtype=np.int64)


def boundary_fingerprint_distances_deg(
    relative_matrices: ArrayLike,
    fingerprint: ArrayLike,
) -> np.ndarray:
    r"""Angular distance from each boundary misorientation to a fingerprint set.

    Purpose: scores measured child-child boundary misorientations against the
    admissible same-parent set from
    :func:`intervariant_boundary_fingerprint`. A distance near zero means the
    boundary is consistent with the two grains sharing a parent.

    Algorithm: the distance to the set is
    :math:`\min_F \angle\!\left(M F^{\mathsf{T}}\right)`, and since
    :math:`\operatorname{tr}\!\left(M F^{\mathsf{T}}\right)` is just the
    elementwise product summed, the whole comparison is one
    ``(e, 9) @ (9, k)`` GEMM followed by a row maximum. It is evaluated in
    blocks of ``512`` edges so the transient stays bounded at map scale — the
    unblocked form allocates ``e * k`` floats, which is several gigabytes for
    50 000 edges against a cubic-cubic fingerprint.

    Inputs: ``(e, 3, 3)`` crystal-frame boundary misorientations
    (:math:`M = C_i^{\mathsf{T}} C_j` under the canonical crystal-to-specimen
    orientation convention) and the ``(k, 3, 3)`` fingerprint. The fingerprint
    is already closed under child symmetry on both sides, so ``M`` must **not**
    be symmetry-reduced first.

    Output: an ``(e,)`` array of degrees.
    """

    matrices = np.asarray(relative_matrices, dtype=np.float64)
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3):
        raise ValueError("relative_matrices must have shape (e, 3, 3).")
    reference = np.asarray(fingerprint, dtype=np.float64)
    if reference.ndim != 3 or reference.shape[1:] != (3, 3):
        raise ValueError("fingerprint must have shape (k, 3, 3).")
    if reference.shape[0] == 0:
        raise ValueError("fingerprint must contain at least one rotation.")
    count = matrices.shape[0]
    flat_edges = np.ascontiguousarray(matrices.reshape(count, 9))
    flat_reference = np.ascontiguousarray(reference.reshape(reference.shape[0], 9).T)
    best_traces = np.empty(count, dtype=np.float64)
    for start in range(0, count, _FINGERPRINT_BLOCK_SIZE):
        stop = min(start + _FINGERPRINT_BLOCK_SIZE, count)
        best_traces[start:stop] = (flat_edges[start:stop] @ flat_reference).max(axis=1)
    cosines = np.clip((best_traces - 1.0) * 0.5, -1.0, 1.0)
    distances = np.degrees(np.arccos(cosines))
    distances = np.ascontiguousarray(distances)
    distances.setflags(write=False)
    return distances


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
        """Mean angular departure of the measured pairs from the nominal OR.

        The headline goodness-of-fit number: how well the nominal relationship
        describes the data, in degrees.
        """

        return float(np.mean(self.deviations_deg))

    @property
    def median_deviation_deg(self) -> float:
        """Median angular departure, in degrees.

        More robust than the mean to a few badly indexed pairs, and the more
        honest statistic when the pair set contains outliers.
        """

        return float(np.median(self.deviations_deg))

    @property
    def max_deviation_deg(self) -> float:
        """Largest angular departure over the measured pairs, in degrees.

        The worst-case number. A large maximum with a small median usually means
        misindexed points or a second, unaccounted-for relationship, rather than
        a poor overall fit.
        """

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
        """Mean angular residual of the fitted relationship, in degrees.

        The fit-quality number for a *fitted* OR. Compare it against the
        deviation from the nearest catalogued relationship to judge whether the
        data support a named OR or a genuinely distinct one.
        """

        return float(np.mean(self.residuals_deg))

    @property
    def max_residual_deg(self) -> float:
        """Largest angular residual of the fitted relationship, in degrees.
        """

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


def _symmetry_operator_pair(
    relationship: OrientationRelationship,
) -> tuple[np.ndarray, np.ndarray]:
    """The (parent, child) proper-rotation operator arrays of a relationship.

    Phases without a symmetry specification contribute the identity, so every
    caller can use one code path.
    """

    identity = np.eye(3, dtype=np.float64)[None, :, :]
    parent_symmetry = relationship.parent_phase.symmetry
    child_symmetry = relationship.child_phase.symmetry
    return (
        parent_symmetry.operators if parent_symmetry is not None else identity,
        child_symmetry.operators if child_symmetry is not None else identity,
    )


def _measured_parent_to_child(
    parent_orientations: OrientationSet, child_orientations: OrientationSet
) -> np.ndarray:
    """Per-pair measured parent-to-child rotations ``V_i = C_i^T P_i``.

    The single definition of "the operative rotation of a measured pair" in
    the canonical crystal->specimen convention (``C = P V^T``). Every TX
    surface routes through this helper rather than re-deriving the transpose
    placement.
    """

    return np.asarray(
        np.einsum(
            "nji,njk->nik",
            child_orientations.as_matrices(),
            parent_orientations.as_matrices(),
            optimize=True,
        ),
        dtype=np.float64,
    )


def _fit_from_seed(
    measured: np.ndarray,
    seed: np.ndarray,
    *,
    parent_operators: np.ndarray,
    child_operators: np.ndarray,
    max_iterations: int,
    convergence_tol_deg: float,
) -> tuple[np.ndarray, np.ndarray, int, bool]:
    """Symmetry-aware rotation averaging of measured parent-to-child rotations.

    Alternates (i) aligning every measurement to the symmetry-equivalent
    description ``S_c V S_p`` nearest the current estimate and (ii) replacing
    the estimate with the quaternion eigen-mean (Markley) of the aligned set.
    Shared by ``fit_orientation_relationship`` (seeded by a nominal OR) and
    ``characterize_orientation_relationship`` (seeded by the double-coset
    reduction of the data itself), so there is exactly one implementation of
    the averaging semantics.

    Returns ``(estimate, residuals_deg, iterations, converged)``.
    """

    candidates = np.einsum(
        "aij,njk,bkl->nabil", child_operators, measured, parent_operators, optimize=True
    )
    pair_count = measured.shape[0]
    flat_candidates = candidates.reshape(pair_count, -1, 3, 3)
    estimate = np.asarray(seed, dtype=np.float64)
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
        step_angle = _rotation_angles_deg_from_matrices((updated @ estimate.T)[None, :, :])[0]
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
    return estimate, residuals, iterations, converged


def _double_coset_seed(
    measured: np.ndarray,
    *,
    parent_operators: np.ndarray,
    child_operators: np.ndarray,
) -> np.ndarray:
    """A starting estimate derived from the data alone, with no nominal OR.

    The first measurement is reduced to its minimum-angle (maximum-trace)
    representative within the double coset ``G_c V_0 G_p`` — the disorientation
    description of the relationship that pair shows. Measurements belonging to
    different variants of one relationship differ exactly by a parent symmetry
    operation, which the coset absorbs, so every other pair has an equivalent
    description close to this one and `_fit_from_seed` will find it.

    Only one measurement is reduced, deliberately. Averaging the reduced
    representatives of *all* measurements looks more robust and is not: the
    maximum-trace element is not unique when the relationship's rotation is
    itself symmetric, and different pairs then reduce to different tied
    representatives whose mean is a rotation none of them shows. Bain is the
    concrete failure — 45 deg about <100> with three variants averages to a
    meaningless 27 deg. Seeding from one pair and letting the alignment step
    resolve every other pair against it breaks the ties consistently.

    Any symmetry-equivalent description states the same relationship, so
    downstream comparisons must remain symmetry-reduced (they are).
    """

    candidates = np.einsum(
        "aij,jk,bkl->abil",
        child_operators,
        measured[0],
        parent_operators,
        optimize=True,
    ).reshape(-1, 3, 3)
    traces = np.trace(candidates, axis1=-2, axis2=-1)
    return np.ascontiguousarray(candidates[int(np.argmax(traces))])


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
    parent_operators, child_operators = _symmetry_operator_pair(nominal)
    measured = _measured_parent_to_child(parent_orientations, child_orientations)
    estimate, residuals, iterations, converged = _fit_from_seed(
        measured,
        nominal.parent_to_child_rotation.as_matrix(),
        parent_operators=parent_operators,
        child_operators=child_operators,
        max_iterations=max_iterations,
        convergence_tol_deg=convergence_tol_deg,
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
        """Prose summary of the parallel plane or direction matches found.

        States the relationship, the tolerance, and each match as a
        variant-indexed parallelism with its angular deviation, using
        crystallographic index notation.
        """

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


#: Default angular tolerance for calling two crystallographic objects parallel.
#:
#: Three degrees is the working figure for EBSD-derived orientation relationships:
#: comfortably above the ~0.5 deg orientation-noise floor of a well-calibrated
#: indexed map, and comfortably below the 5.26 deg that separates
#: Kurdjumov-Sachs from Nishiyama-Wassermann, so the two remain distinguishable.
DEFAULT_OR_TOLERANCE_DEG = 3.0


def _canonical_sign_triples(max_index: int) -> np.ndarray:
    """Primitive integer triples with the first nonzero entry positive.

    A plane and its negative describe the same plane, and a direction and its
    reverse the same axis of parallelism, so only one of each antiparallel pair
    is worth testing.
    """

    triples = _primitive_integer_triples(max_index)
    first_nonzero = np.argmax(triples != 0, axis=1)
    leading = triples[np.arange(triples.shape[0]), first_nonzero]
    return np.ascontiguousarray(triples[leading > 0])


def _family_key(
    indices: np.ndarray, *, phase: Phase, reciprocal: bool
) -> tuple[int, int, int]:
    """A deterministic identifier shared by every member of one index family.

    Used to recognize that two parallelism clauses are the same crystallographic
    statement written with different family members, so only one is reported.
    """

    orbit = _integer_index_orbit(
        np.asarray(indices, dtype=np.int64), phase=phase, reciprocal=reciprocal
    )
    return min(_index_tuple(row) for row in orbit)


def _crystallographic_label(indices: np.ndarray, *, phase: Phase, reciprocal: bool) -> str:
    """Format a plane or direction the way the literature writes it for this phase.

    Hexagonal phases take four-index Miller-Bravais labels — three-index
    hexagonal indices hide the symmetry of the family and are not how the
    hcp literature states an orientation relationship — while every other
    system keeps its three-index form.
    """

    triple = _index_tuple(indices)
    if is_hexagonal_phase(phase):
        four = plane_hkl_to_hkil(triple) if reciprocal else direction_uvw_to_uvtw(triple)
        values = tuple(int(value) for value in four)
        return (
            format_plane_indices(values, style="plain")
            if reciprocal
            else format_direction_indices(values, style="plain")
        )
    return (
        format_plane_indices(triple, style="plain")
        if reciprocal
        else format_direction_indices(triple, style="plain")
    )


#: The angle below which two directions are parallel rather than nearly so.
#:
#: A clause that is exactly parallel by construction computes to a deviation of
#: order ``1e-14`` degrees, not to zero, so a bit-exact comparison against the
#: caller's tolerance would refuse it -- and ``tolerance_deg=0.0``, which reads
#: as "only exact parallelism", would accept nothing at all. This is the width
#: of that rounding, five orders below any tolerance a crystallographer would
#: type and far above the residue itself.
_PARALLEL_ROUNDING_DEG = 1e-9


@dataclass(frozen=True, slots=True)
class ORParallelismStatement:
    """One ``(hkl)_p || (hkl)_c`` or ``[uvw]_p || [uvw]_c`` clause of an OR statement.

    This is how the literature states an orientation relationship: a rotation
    matrix is unreadable, but "(111) austenite parallel to (011) ferrite" is the
    working crystallographic fact. ``deviation_deg`` is the angle between the
    exact image of the parent object and the reported child indices; an exact
    defining parallelism reads zero.
    """

    kind: str
    parent_indices: np.ndarray
    child_indices: np.ndarray
    deviation_deg: float
    parent_label: str
    child_label: str

    def __post_init__(self) -> None:
        if self.kind not in {"plane", "direction"}:
            raise ValueError("ORParallelismStatement.kind must be 'plane' or 'direction'.")
        parent = np.asarray(self.parent_indices, dtype=np.int64)
        child = np.asarray(self.child_indices, dtype=np.int64)
        if parent.shape != (3,) or child.shape != (3,):
            raise ValueError("ORParallelismStatement indices must have shape (3,).")
        if not np.isfinite(self.deviation_deg) or self.deviation_deg < 0.0:
            raise ValueError("deviation_deg must be finite and non-negative.")
        for array in (parent, child):
            array.setflags(write=False)
        object.__setattr__(self, "parent_indices", parent)
        object.__setattr__(self, "child_indices", child)

    def as_text(self) -> str:
        """The clause as it would be written in a paper, e.g. ``(111) || (011)``."""

        return f"{self.parent_label} || {self.child_label}"


def _parallelism_statements(
    relationship: OrientationRelationship,
    *,
    kind: str,
    rotation: np.ndarray,
    tolerance_deg: float,
    max_index: int,
    max_statements: int,
    preferred_parent_indices: tuple[np.ndarray, ...] = (),
) -> tuple[ORParallelismStatement, ...]:
    """Extract the best near-parallel index pairs of one kind, vectorized.

    Every canonical-sign primitive parent triple up to ``max_index`` is carried
    through ``rotation`` into the child basis at once, and every candidate child
    triple is compared with every image in a single cosine matrix, rather than
    looping in Python over ~10^2 index pairs and re-running the rationalizer for
    each. Clauses are then deduplicated by index *family*, because the same
    statement written with a different family member is not new information.

    ``preferred_parent_indices`` nominates parent families to report first. This
    matters because a rotation generally satisfies several equally exact
    low-index parallelisms, and which one is *the* statement is a fact about the
    two structures (their close-packed planes and directions), not about the
    rotation — index magnitude alone cannot recover it. When the caller knows
    the relationship is Kurdjumov-Sachs, {111} and <110> are the families worth
    reporting, and the deviations then verify that the fitted rotation really
    does satisfy them.
    """

    reciprocal = kind == "plane"
    parent_phase = relationship.parent_phase
    child_phase = relationship.child_phase
    source_basis = (
        parent_phase.lattice.reciprocal_basis().matrix
        if reciprocal
        else parent_phase.lattice.direct_basis().matrix
    )
    target_basis = (
        child_phase.lattice.reciprocal_basis().matrix
        if reciprocal
        else child_phase.lattice.direct_basis().matrix
    )
    parent_triples = _canonical_sign_triples(max_index)
    child_triples = _canonical_sign_triples(max_index)

    # Cartesian images of every parent triple, rotated into the child frame.
    parent_cartesian = parent_triples.astype(np.float64) @ source_basis.T
    images = parent_cartesian @ rotation.T
    image_units = images / np.linalg.norm(images, axis=1)[:, None]
    # Cartesian images of every candidate child triple, in the same frame.
    child_cartesian = child_triples.astype(np.float64) @ target_basis.T
    child_units = child_cartesian / np.linalg.norm(child_cartesian, axis=1)[:, None]

    # Signed comparison would reject a correct antiparallel description, and the
    # canonical-sign filter has already collapsed each antiparallel pair to one
    # representative, so parallelism is judged on |cos|.
    cosines = np.abs(image_units @ child_units.T)
    best_columns = np.argmax(cosines, axis=1)
    # The winner is found from the cosine matrix, which is what an all-pairs
    # comparison can afford; the deviation itself is then taken from the two
    # vectors. An exactly parallel clause -- which is the whole point of a
    # rational orientation relationship -- must read as zero, and arccos of a
    # cosine cannot report better than about 1e-06 deg. See pytex.core._angles.
    deviations = np.degrees(
        acute_angle_between_unit_vectors_rad(image_units, child_units[best_columns])
    )
    accepted = np.flatnonzero(deviations <= tolerance_deg + _PARALLEL_ROUNDING_DEG)

    parent_keys = {
        int(row): _family_key(parent_triples[row], phase=parent_phase, reciprocal=reciprocal)
        for row in accepted
    }
    preferred_keys = {
        _family_key(np.asarray(indices, dtype=np.int64), phase=parent_phase, reciprocal=reciprocal)
        for indices in preferred_parent_indices
    }
    # Fit quality outranks preference: a nominated family must not promote a
    # visibly worse clause above an exact one. Bucketing the deviation to
    # milli-degrees first lets preference decide among clauses that are equally
    # exact in any meaningful sense, while a genuinely poorer fit still sorts
    # below. The bucket used to have a second job -- absorbing the ~1e-6 deg
    # floor of an arccos-recovered angle -- which it no longer needs.
    ranked = sorted(
        accepted,
        key=lambda row: (
            round(float(deviations[row]), 3),
            0 if parent_keys[int(row)] in preferred_keys else 1,
            round(float(deviations[row]), 9),
            int(np.abs(parent_triples[row]).sum() + np.abs(child_triples[best_columns[row]]).sum()),
            _index_tuple(parent_triples[row]),
        ),
    )
    statements: list[ORParallelismStatement] = []
    seen: set[tuple[tuple[int, int, int], tuple[int, int, int]]] = set()
    for row in ranked:
        parent_indices = parent_triples[row]
        child_indices = child_triples[best_columns[row]]
        key = (
            parent_keys[int(row)],
            _family_key(child_indices, phase=child_phase, reciprocal=reciprocal),
        )
        if key in seen:
            continue
        seen.add(key)
        statements.append(
            ORParallelismStatement(
                kind=kind,
                parent_indices=parent_indices,
                child_indices=child_indices,
                deviation_deg=float(deviations[row]),
                parent_label=_crystallographic_label(
                    parent_indices, phase=parent_phase, reciprocal=reciprocal
                ),
                child_label=_crystallographic_label(
                    child_indices, phase=child_phase, reciprocal=reciprocal
                ),
            )
        )
        if len(statements) >= max_statements:
            break
    return tuple(statements)


def describe_orientation_relationship(
    relationship: OrientationRelationship,
    *,
    tolerance_deg: float = DEFAULT_OR_TOLERANCE_DEG,
    max_index: int = 3,
    max_statements: int = 4,
    variant: TransformationVariant | None = None,
    preferred_parent_planes: Sequence[ArrayLike] | None = None,
    preferred_parent_directions: Sequence[ArrayLike] | None = None,
) -> tuple[tuple[ORParallelismStatement, ...], tuple[ORParallelismStatement, ...]]:
    """Recover the parallel-plane / parallel-direction statement of an OR.

    Purpose: turns a rotation back into crystallography. Given any orientation
    relationship — named, fitted from measurements, or hand-built — this finds
    which low-index parent planes are parallel to which low-index child planes,
    and likewise for directions, so the relationship can be reported the way the
    literature reports it rather than as a matrix.

    When to use: after fitting an OR to measured data (the statement is the
    interpretable answer), when checking that a constructed relationship
    encodes the parallelisms intended, and when teaching what a named OR means.

    A rotation has three degrees of freedom, so one plane parallelism (two
    constraints) plus one in-plane direction parallelism (the third) determines
    the relationship completely. That pair is exactly the classical form, e.g.
    Kurdjumov-Sachs as ``(111) || (011)`` with ``[10-1] || [11-1]``.

    Inputs: the relationship; ``tolerance_deg``, the angle below which two
    objects count as parallel (see :data:`DEFAULT_OR_TOLERANCE_DEG`);
    ``max_index``, the index bound searched on both sides; ``max_statements``,
    how many clauses of each kind to return; an optional ``variant`` to
    describe one specific variant instead of the base relationship; and
    ``preferred_parent_planes`` / ``preferred_parent_directions``, parent index
    families to report first — which default to the relationship's own recorded
    defining parallelisms when it was built from an explicit correspondence.

    Output: ``(plane_statements, direction_statements)``, each ordered by
    preference, then deviation, then total index magnitude, and deduplicated by
    index family so one crystallographic statement is reported once.

    A rotation typically satisfies *several* exact low-index parallelisms at
    once, all of them true. Which one is quoted in the literature is decided by
    the structures — the close-packed plane and direction of the two phases —
    and index magnitude alone cannot recover that, so nominate the families of
    interest through the ``preferred_*`` arguments when they are known.
    `characterize_orientation_relationship` does this automatically from the
    matching catalog relationship's own defining parallelisms.

    Note also that a relationship has several symmetry-equivalent descriptions,
    and the statement recovered is the one belonging to the stored rotation. A
    different but equivalent description (e.g. ``(1-11)`` in place of ``(111)``)
    states the same relationship.

    See also
    --------
    `find_parallel_planes`, `find_parallel_directions` : the per-variant search
        over one nominated parent family.
    `characterize_orientation_relationship` : fits an OR to measurements and
        reports its statement in one call.
    """

    if max_statements < 1:
        raise ValueError("max_statements must be at least 1.")
    rotation = (
        relationship.parent_to_child_rotation.as_matrix()
        if variant is None
        else variant.parent_to_child_rotation.as_matrix()
    )
    # A relationship built from an explicit correspondence already records the
    # families that define it, and those are the ones worth reporting. Falling
    # back to them makes the default answer the defining statement rather than
    # an arbitrary equally-exact alternative.
    if preferred_parent_planes is None:
        preferred_parent_planes = [
            pair[0].miller.indices for pair in relationship.parallel_planes
        ]
    if preferred_parent_directions is None:
        preferred_parent_directions = [
            pair[0].coordinates for pair in relationship.parallel_directions
        ]
    return (
        _parallelism_statements(
            relationship,
            kind="plane",
            rotation=rotation,
            tolerance_deg=tolerance_deg,
            max_index=max_index,
            max_statements=max_statements,
            preferred_parent_indices=tuple(
                np.asarray(indices, dtype=np.int64).reshape(3)
                for indices in preferred_parent_planes
            ),
        ),
        _parallelism_statements(
            relationship,
            kind="direction",
            rotation=rotation,
            tolerance_deg=tolerance_deg,
            max_index=max_index,
            max_statements=max_statements,
            preferred_parent_indices=tuple(
                np.rint(np.asarray(indices, dtype=np.float64).reshape(3)).astype(np.int64)
                for indices in preferred_parent_directions
            ),
        ),
    )


#: How many candidate clauses of each kind the rationalizer considers.
#:
#: The best plane clause is taken outright, but the direction must also lie
#: in it, so the search needs a bench of direction candidates to pick a
#: zone-consistent one from. Sixteen is far past the point where a further
#: candidate could win: they are ranked by deviation first.
_RATIONALIZATION_CANDIDATES = 16


@dataclass(frozen=True, slots=True)
class RationalizedORResult:
    """A measured orientation relationship restated in integers, **with its cost**.

    Purpose
    -------
    A fitted rotation is a measurement; ``(111)_gamma || (011)_alpha`` with
    ``[-101]_gamma || [-1-11]_alpha`` is a crystallographic statement. Turning
    the first into the second is an *idealization*: the integer statement is a
    nearby exact relationship, not the one measured, and the difference is a
    real angle. This type carries both, because a rational OR handed back
    without its cost is exactly the kind of silent boundary error the library
    exists to prevent.

    Attributes
    ----------
    relationship : OrientationRelationship
        The idealized relationship, built from the integer plane and direction
        pair by ``from_parallel_plane_direction``. It is exact by construction —
        the deviations below say how far it sits from the measurement, not how
        far it is from being a relationship.
    plane_statement, direction_statement : ORParallelismStatement
        The clauses it was built from. Their ``deviation_deg`` is the angle
        between the exact image of the parent object under the *measured*
        rotation and the integer child indices reported — the per-clause cost.
    residual_rotation_deg : float
        The symmetry-reduced angle between the measured rotation and the
        idealized one. This is the number to quote: it is the whole cost of the
        idealization, in the same units as the fit residual it should be
        compared against.
    zone_law_deviation_deg : float
        How far the chosen parent direction departs from lying in the chosen
        parent plane. ``from_parallel_plane_direction`` removes the normal
        component of the direction, so a pair that failed the zone law would
        build a relationship the two labels do not describe. Zero for an
        exact-zone pair, and clauses are selected to keep it so.
    max_index, tolerance_deg : int, float
        The search bounds the statement was found under, carried so a reader can
        tell "no better statement exists" from "none was looked for".

    Notes
    -----
    Compare ``residual_rotation_deg`` against the measurement's own scatter. An
    idealization costing less than the scatter is free — the data cannot tell
    the two relationships apart. One costing several times the scatter is a
    claim the data contradicts, however tidy the integers look.
    """

    relationship: OrientationRelationship
    source_relationship_name: str
    plane_statement: ORParallelismStatement
    direction_statement: ORParallelismStatement
    residual_rotation_deg: float
    zone_law_deviation_deg: float
    max_index: int
    tolerance_deg: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        for name in ("residual_rotation_deg", "zone_law_deviation_deg"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if self.plane_statement.kind != "plane":
            raise ValueError("plane_statement must be a plane clause.")
        if self.direction_statement.kind != "direction":
            raise ValueError("direction_statement must be a direction clause.")

    @property
    def statement(self) -> str:
        """The relationship as a paper would write it."""

        return (
            f"{self.plane_statement.as_text()} with {self.direction_statement.as_text()}"
        )

    def describe(self) -> str:
        """Prose: the integer statement, and what idealizing to it cost.

        The cost is stated in the same sentence as the statement, and never in a
        later one, because the two are only meaningful together.
        """

        zone_text = (
            ""
            if self.zone_law_deviation_deg <= 1e-9
            else (
                f" The chosen direction departs from the chosen plane by "
                f"{self.zone_law_deviation_deg:.4f} deg, so the constructed relationship uses "
                "its in-plane component."
            )
        )
        return (
            f"Rationalized orientation relationship from '{self.source_relationship_name}': "
            f"{self.statement}. This is an **idealization**, not the measurement: the integer "
            f"statement is a nearby exact relationship sitting "
            f"{self.residual_rotation_deg:.3f} deg (symmetry-reduced) from the rotation it was "
            f"derived from. Per clause, the plane pair deviates by "
            f"{self.plane_statement.deviation_deg:.4f} deg and the direction pair by "
            f"{self.direction_statement.deviation_deg:.4f} deg, both measured as the angle "
            "between the exact image of the parent object and the integer child indices "
            f"reported. Indices were searched to |index| <= {self.max_index} within "
            f"{self.tolerance_deg:.3f} deg.{zone_text} Compare the residual against the "
            "scatter of the measurement it came from: an idealization costing less than the "
            "scatter is one the data cannot distinguish, and one costing several times the "
            "scatter is a claim the data contradicts."
        )


@dataclass(frozen=True, slots=True)
class ORCharacterizationReport:
    """What orientation relationship a set of measured parent/child pairs shows.

    The one-call answer to "I measured these orientations by EBSD — what is the
    OR?". It carries four separable things, because a trustworthy answer needs
    all of them:

    - ``relationship``: the rotation fitted to the measurements;
    - ``residuals_deg``: how tightly the pairs agree with it (the scatter);
    - ``catalog_names`` / ``catalog_deviations_deg`` / ``best_catalog_name``:
      which named relationship it is, and by how much it beats the runner-up;
    - ``plane_statements`` / ``direction_statements``: the parallelisms, i.e.
      the crystallographic reading of the fitted rotation.

    ``is_conclusive`` is deliberately conservative: a named match is only
    claimed when the winner both fits within ``catalog_tolerance_deg`` and
    leads the runner-up by more than the data scatter and its own misfit.
    """

    relationship: OrientationRelationship
    pair_count: int
    residuals_deg: np.ndarray
    iterations: int
    converged: bool
    catalog_names: tuple[str, ...]
    catalog_deviations_deg: np.ndarray
    best_catalog_name: str | None
    best_catalog_deviation_deg: float
    margin_deg: float
    catalog_tolerance_deg: float
    plane_statements: tuple[ORParallelismStatement, ...]
    direction_statements: tuple[ORParallelismStatement, ...]
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        residuals = np.asarray(self.residuals_deg, dtype=np.float64).reshape(-1)
        deviations = np.asarray(self.catalog_deviations_deg, dtype=np.float64).reshape(-1)
        names = tuple(str(name) for name in self.catalog_names)
        if residuals.size == 0:
            raise ValueError("ORCharacterizationReport requires at least one pair.")
        if self.pair_count != residuals.size:
            raise ValueError("pair_count must equal the number of residuals.")
        if deviations.shape != (len(names),):
            raise ValueError("catalog_deviations_deg must have one entry per catalog name.")
        if self.best_catalog_name is not None and self.best_catalog_name not in names:
            raise ValueError("best_catalog_name must be one of catalog_names.")
        if names and self.best_catalog_name is None:
            raise ValueError("best_catalog_name is required when a catalog was compared.")
        if np.any(~np.isfinite(residuals)) or np.any(residuals < 0.0):
            raise ValueError("residuals_deg must be finite and non-negative.")
        if deviations.size and (np.any(~np.isfinite(deviations)) or np.any(deviations < 0.0)):
            raise ValueError("catalog_deviations_deg must be finite and non-negative.")
        if self.iterations <= 0:
            raise ValueError("iterations must be positive.")
        for array in (residuals, deviations):
            array.setflags(write=False)
        object.__setattr__(self, "residuals_deg", residuals)
        object.__setattr__(self, "catalog_deviations_deg", deviations)
        object.__setattr__(self, "catalog_names", names)
        object.__setattr__(self, "plane_statements", tuple(self.plane_statements))
        object.__setattr__(self, "direction_statements", tuple(self.direction_statements))

    @property
    def mean_residual_deg(self) -> float:
        """Mean angular residual of the characterization, in degrees.
        """

        return float(np.mean(self.residuals_deg))

    @property
    def max_residual_deg(self) -> float:
        """Largest angular residual of the characterization, in degrees.
        """

        return float(np.max(self.residuals_deg))

    def as_rational_relationship(
        self,
        *,
        max_index: int = 4,
        tolerance_deg: float = 3.0,
        name: str | None = None,
    ) -> RationalizedORResult:
        """Restate the fitted relationship in integers, and price the idealization.

        Purpose
        -------
        Turns the measurement into the object a paper states: a genuine
        `OrientationRelationship` built from an integer plane pair and an
        integer direction pair. The result carries the cost of doing so, which
        is the point — an idealization reported without its cost cannot be
        judged against the scatter of the data it came from.

        Parameters
        ----------
        max_index : int
            Largest absolute index considered on either side. Four covers every
            statement the martensite literature is written with. The search is
            re-run at these bounds rather than reusing the report's own clauses,
            so widening it here really does widen it.
        tolerance_deg : float
            How far a clause may deviate and still be considered.
        name : str, optional
            Name for the idealized relationship; defaults to the fitted
            relationship's name with a ``_rationalized`` suffix, so an
            idealization is never later mistaken for the measurement.

        Returns
        -------
        RationalizedORResult

        Raises
        ------
        ValueError
            If no plane clause, or no direction clause lying in the chosen
            plane, was found within the bounds. The two messages differ, because
            "no statement" and "no *consistent* statement" call for different
            responses.

        Notes
        -----
        Two selection rules do the work.

        The **families the characterization already chose** are preferred, so
        widening the index bound cannot quietly swap a canonical statement for
        an equally exact alternative with larger indices.

        The direction is then chosen from the clauses satisfying the **zone
        law** against the chosen plane. This is not a refinement:
        ``from_parallel_plane_direction`` removes the normal component of the
        direction, so a plane and a direction that is not in it would build a
        relationship the two printed labels do not describe. Selecting a
        zone-consistent pair keeps the object and its statement the same thing.
        """

        if max_index < 1:
            raise ValueError("max_index must be at least 1.")
        rotation = self.relationship.parent_to_child_rotation.as_matrix()
        preferred_planes = tuple(
            np.asarray(statement.parent_indices, dtype=np.int64)
            for statement in self.plane_statements[:1]
        )
        preferred_directions = tuple(
            np.asarray(statement.parent_indices, dtype=np.int64)
            for statement in self.direction_statements[:1]
        )
        planes = _parallelism_statements(
            self.relationship,
            kind="plane",
            rotation=rotation,
            tolerance_deg=tolerance_deg,
            max_index=max_index,
            max_statements=_RATIONALIZATION_CANDIDATES,
            preferred_parent_indices=preferred_planes,
        )
        directions = _parallelism_statements(
            self.relationship,
            kind="direction",
            rotation=rotation,
            tolerance_deg=tolerance_deg,
            max_index=max_index,
            max_statements=_RATIONALIZATION_CANDIDATES,
            preferred_parent_indices=preferred_directions,
        )
        if not planes:
            raise ValueError(
                f"No plane parallelism within {tolerance_deg:.3f} deg and |index| <= "
                f"{max_index}, so there is nothing to rationalize. Widen the bounds, or "
                "report the rotation instead of a statement."
            )
        plane_statement = planes[0]
        plane_indices = np.asarray(plane_statement.parent_indices, dtype=np.int64)
        zone_consistent = [
            statement
            for statement in directions
            if int(np.dot(plane_indices, np.asarray(statement.parent_indices, dtype=np.int64)))
            == 0
        ]
        if not zone_consistent:
            raise ValueError(
                f"No direction clause lies in {plane_statement.parent_label} within "
                f"{tolerance_deg:.3f} deg and |index| <= {max_index}, so no integer "
                "plane-and-direction statement describes this relationship. Widen "
                "`max_index`, or state the plane parallelism alone."
            )
        direction_statement = zone_consistent[0]

        parent_phase = self.relationship.parent_phase
        child_phase = self.relationship.child_phase
        parent_plane = CrystalPlane(
            MillerIndex(plane_indices, phase=parent_phase), phase=parent_phase
        )
        child_plane = CrystalPlane(
            MillerIndex(
                np.asarray(plane_statement.child_indices, dtype=np.int64), phase=child_phase
            ),
            phase=child_phase,
        )
        parent_direction = CrystalDirection(
            np.asarray(direction_statement.parent_indices, dtype=np.float64), phase=parent_phase
        )
        child_direction = CrystalDirection(
            np.asarray(direction_statement.child_indices, dtype=np.float64), phase=child_phase
        )
        idealized = OrientationRelationship.from_parallel_plane_direction(
            name=name or f"{self.relationship.name}_rationalized",
            parent_plane=parent_plane,
            child_plane=child_plane,
            parent_direction=parent_direction,
            child_direction=child_direction,
            provenance=self.provenance,
        )

        parent_operators, child_operators = _symmetry_operator_pair(self.relationship)
        residual = _symmetry_reduced_angle_between_deg(
            idealized.parent_to_child_rotation.as_matrix(),
            rotation,
            child_operators=child_operators,
            parent_operators=parent_operators,
        )
        # Zero for a zone-consistent pair, but measured rather than asserted:
        # the selection above is an integer test, and this is the geometric
        # consequence a reader can check.
        cosine = abs(float(np.dot(parent_plane.normal, parent_direction.unit_vector)))
        zone_deviation = abs(90.0 - float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))))
        return RationalizedORResult(
            relationship=idealized,
            source_relationship_name=self.relationship.name,
            plane_statement=plane_statement,
            direction_statement=direction_statement,
            residual_rotation_deg=float(residual),
            zone_law_deviation_deg=float(zone_deviation),
            max_index=int(max_index),
            tolerance_deg=float(tolerance_deg),
            provenance=self.provenance,
        )

    @property
    def matches_catalog(self) -> bool:
        """Whether the fitted OR sits within tolerance of a named relationship."""

        return (
            self.best_catalog_name is not None
            and self.best_catalog_deviation_deg <= self.catalog_tolerance_deg
        )

    @property
    def is_conclusive(self) -> bool:
        """Whether the named identification can be trusted.

        Requires the winner to fit within tolerance *and* to lead the runner-up
        by more than both the measurement scatter and its own misfit — the two
        quantities that could otherwise explain the lead away.
        """

        if not self.matches_catalog:
            return False
        if len(self.catalog_names) == 1:
            return True
        return self.margin_deg > max(self.mean_residual_deg, self.best_catalog_deviation_deg)

    def statement_text(self) -> str:
        """The parallelism statement as one line, e.g. ``(111) || (011), [10-1] || [11-1]``."""

        clauses = [statement.as_text() for statement in self.plane_statements[:1]]
        clauses += [statement.as_text() for statement in self.direction_statements[:1]]
        return ", ".join(clauses) if clauses else "no low-index parallelism within tolerance"

    def describe(self) -> str:
        """Prose summary: fitted rotation, scatter, named match, and parallelisms."""

        misorientation = self.relationship.misorientation()
        axis = misorientation.rotation.axis
        lines = [
            f"Orientation relationship characterized from {self.pair_count} measured "
            f"parent/child orientation pair(s). Fitted parent-to-child rotation: "
            f"{misorientation.angle_deg:.3f} deg about "
            f"<{axis[0]:.4f} {axis[1]:.4f} {axis[2]:.4f}> (disorientation representative, "
            f"parent {self.relationship.parent_phase.name} -> child "
            f"{self.relationship.child_phase.name}). Pair scatter about the fit: mean "
            f"{self.mean_residual_deg:.3f} deg, max {self.max_residual_deg:.3f} deg."
        ]
        if self.pair_count == 1:
            lines.append(
                "  Only one pair was supplied, so the scatter is zero by construction and "
                "says nothing about measurement quality; supply several pairs to estimate it."
            )
        if self.catalog_names:
            order = np.argsort(self.catalog_deviations_deg)
            ranking = "; ".join(
                f"{self.catalog_names[int(index)]}: "
                f"{float(self.catalog_deviations_deg[int(index)]):.3f} deg"
                for index in order[:5]
            )
            verdict = (
                f"identified as '{self.best_catalog_name}'"
                if self.is_conclusive
                else (
                    f"closest to '{self.best_catalog_name}' but NOT conclusively identified"
                    if self.matches_catalog
                    else "matches no catalog relationship within "
                    f"{self.catalog_tolerance_deg:.3f} deg"
                )
            )
            lines.append(
                f"  Catalog comparison ({len(self.catalog_names)} candidate(s)): {verdict} "
                f"at {self.best_catalog_deviation_deg:.3f} deg, leading the runner-up by "
                f"{self.margin_deg:.3f} deg. Ranking: {ranking}."
            )
            if not self.is_conclusive:
                lines.append(
                    "  A lead comparable to the scatter or to the winner's own misfit cannot "
                    "separate the candidates; treat the name as provisional."
                )
        else:
            lines.append(
                "  No relationship catalog was available for this phase pair, so the fitted "
                "rotation is reported without a name."
            )
        if self.plane_statements or self.direction_statements:
            lines.append("  Crystallographic statement of the fitted relationship:")
            for statement in self.plane_statements:
                lines.append(
                    f"    plane     {statement.as_text()} "
                    f"(deviation {statement.deviation_deg:.4f} deg)"
                )
            for statement in self.direction_statements:
                lines.append(
                    f"    direction {statement.as_text()} "
                    f"(deviation {statement.deviation_deg:.4f} deg)"
                )
            lines.append(
                "    One plane clause fixes two of the rotation's three degrees of freedom "
                "and one in-plane direction clause fixes the third; the remaining clauses are "
                "consequences. Indices are those of the stored description; a "
                "symmetry-equivalent description states the same relationship."
            )
        else:
            lines.append(
                "  No low-index plane or direction parallelism was found within tolerance, "
                "which is itself informative: the relationship is not of the classical "
                "parallel-planes type at this index bound."
            )
        return "\n".join(lines)

    def to_json_dict(self) -> dict[str, object]:
        """A plain JSON-ready summary, carrying the same facts as ``describe()``.

        This is a one-way report payload for manifests and downstream tooling,
        not a round-trip contract; use `pytex.contracts` for objects that must
        be reconstructed.
        """

        def _statements(items: tuple[ORParallelismStatement, ...]) -> list[dict[str, object]]:
            return [
                {
                    "kind": statement.kind,
                    "parent_indices": _index_tuple(statement.parent_indices),
                    "child_indices": _index_tuple(statement.child_indices),
                    "parent_label": statement.parent_label,
                    "child_label": statement.child_label,
                    "deviation_deg": statement.deviation_deg,
                }
                for statement in items
            ]

        misorientation = self.relationship.misorientation()
        return {
            "schema": "pytex.or_characterization_report/1",
            "parent_phase": self.relationship.parent_phase.name,
            "child_phase": self.relationship.child_phase.name,
            "pair_count": int(self.pair_count),
            "rotation_angle_deg": float(misorientation.angle_deg),
            "rotation_axis": [float(value) for value in misorientation.rotation.axis],
            "mean_residual_deg": self.mean_residual_deg,
            "max_residual_deg": self.max_residual_deg,
            "converged": bool(self.converged),
            "iterations": int(self.iterations),
            "catalog_names": list(self.catalog_names),
            "catalog_deviations_deg": [float(value) for value in self.catalog_deviations_deg],
            "best_catalog_name": self.best_catalog_name,
            "best_catalog_deviation_deg": float(self.best_catalog_deviation_deg),
            "margin_deg": float(self.margin_deg),
            "catalog_tolerance_deg": float(self.catalog_tolerance_deg),
            "matches_catalog": bool(self.matches_catalog),
            "is_conclusive": bool(self.is_conclusive),
            "statement_text": self.statement_text(),
            "plane_statements": _statements(self.plane_statements),
            "direction_statements": _statements(self.direction_statements),
        }


def characterize_orientation_relationship(
    parent_orientations: OrientationSet,
    child_orientations: OrientationSet,
    *,
    catalog: OrientationRelationshipCatalog | Sequence[OrientationRelationship] | None = None,
    nominal: OrientationRelationship | None = None,
    catalog_tolerance_deg: float = DEFAULT_OR_TOLERANCE_DEG,
    parallelism_tolerance_deg: float = DEFAULT_OR_TOLERANCE_DEG,
    max_index: int = 3,
    max_statements: int = 4,
    max_iterations: int = 20,
    convergence_tol_deg: float = 1e-8,
    provenance: ProvenanceRecord | None = None,
) -> ORCharacterizationReport:
    """Determine the orientation relationship shown by measured orientation pairs.

    Purpose: answers "I measured a parent grain and one or more child grains by
    EBSD — what is the orientation relationship?" in one call, and answers it
    the way the question is meant: not only as a rotation, but as a named
    relationship with an honest confidence verdict and the parallel planes and
    directions that define it.

    When to use: on paired parent/child grain-mean orientations from an EBSD map
    of a partially transformed microstructure (retained austenite plus
    martensite, retained beta plus alpha, a precipitate and its matrix). When no
    parent phase survives, use
    `pytex.experimental.identify_orientation_relationship` instead, which works
    from child-child boundary misorientations alone.

    Algorithm: each pair contributes the measured rotation ``V_i = C_i^T P_i``.
    Without a nominal relationship the starting estimate comes from the data
    itself — every ``V_i`` is reduced to its minimum-angle representative in the
    double coset ``G_c V_i G_p``, which absorbs the parent symmetry operation
    that distinguishes one variant from another, so pairs belonging to different
    variants reduce to the same matrix and can be averaged. The estimate is then
    refined by symmetry-aware rotation averaging (align each measurement to its
    nearest equivalent description, take the quaternion eigen-mean, iterate) —
    the same routine `fit_orientation_relationship` uses. Finally the fit is
    compared with each catalog relationship under both symmetry groups, and its
    parallelisms are extracted by `describe_orientation_relationship`.

    Inputs: paired ``OrientationSet`` objects of equal length sharing a specimen
    frame, one per phase. ``catalog`` accepts an ``OrientationRelationshipCatalog``
    or a tuple of relationships; when omitted, the standard catalog for the two
    crystal systems is used (see `pytex.core.parent_reconstruction.default_relationship_catalog`).
    ``nominal`` overrides the data-derived starting estimate. The two tolerances
    govern the named match and the parallelism search respectively.

    Output: an `ORCharacterizationReport` — read its ``describe()``.

    Examples
    --------
    A parent and one child built through Kurdjumov-Sachs variant 1 are
    identified as Kurdjumov-Sachs at zero deviation, and the recovered statement
    is the defining one, ``{111} || {011}`` with ``<10-1> || <11-1>``.

    See also
    --------
    `fit_orientation_relationship` : the refinement step alone, given a nominal.
    `or_deviation` : how far measurements sit from a relationship already known.
    `describe_orientation_relationship` : the parallelism extraction alone.
    """

    if len(parent_orientations) != len(child_orientations):
        raise ValueError(
            "parent_orientations and child_orientations must be paired (equal length)."
        )
    if len(parent_orientations) == 0:
        raise ValueError("characterize_orientation_relationship requires at least one pair.")
    if parent_orientations.specimen_frame != child_orientations.specimen_frame:
        raise ValueError("Parent and child orientations must share a specimen frame.")
    parent_phase = parent_orientations.phase
    child_phase = child_orientations.phase
    if parent_phase is None or child_phase is None:
        raise ValueError(
            "characterize_orientation_relationship requires both orientation sets to carry "
            "a phase; the relationship is defined between phases."
        )
    if max_statements < 1:
        raise ValueError("max_statements must be at least 1.")

    from pytex.core.parent_reconstruction import (
        OrientationRelationshipCatalog as _Catalog,
    )
    from pytex.core.parent_reconstruction import (
        default_relationship_catalog,
    )

    candidates: tuple[OrientationRelationship, ...]
    if catalog is None:
        resolved_catalog = default_relationship_catalog(
            parent_phase=parent_phase, child_phase=child_phase, provenance=provenance
        )
        candidates = () if resolved_catalog is None else resolved_catalog.relationships
    elif isinstance(catalog, _Catalog):
        candidates = catalog.relationships
    else:
        candidates = tuple(catalog)

    seed_relationship = nominal or OrientationRelationship(
        name="seed",
        parent_phase=parent_phase,
        child_phase=child_phase,
        parent_to_child_rotation=Rotation.identity(),
        provenance=provenance,
    )
    parent_operators, child_operators = _symmetry_operator_pair(seed_relationship)
    measured = _measured_parent_to_child(parent_orientations, child_orientations)
    seed = (
        nominal.parent_to_child_rotation.as_matrix()
        if nominal is not None
        else _double_coset_seed(
            measured, parent_operators=parent_operators, child_operators=child_operators
        )
    )
    estimate, residuals, iterations, converged = _fit_from_seed(
        measured,
        seed,
        parent_operators=parent_operators,
        child_operators=child_operators,
        max_iterations=max_iterations,
        convergence_tol_deg=convergence_tol_deg,
    )
    fitted = OrientationRelationship(
        name=f"{parent_phase.name}_to_{child_phase.name}_fitted",
        parent_phase=parent_phase,
        child_phase=child_phase,
        parent_to_child_rotation=Rotation.from_matrix(estimate).canonicalized(),
        provenance=provenance,
    )
    deviations = np.asarray(
        [
            _symmetry_reduced_angle_between_deg(
                estimate,
                candidate.parent_to_child_rotation.as_matrix(),
                child_operators=child_operators,
                parent_operators=parent_operators,
            )
            for candidate in candidates
        ],
        dtype=np.float64,
    )
    preferred_planes: tuple[np.ndarray, ...] = ()
    preferred_directions: tuple[np.ndarray, ...] = ()
    if deviations.size:
        order = np.argsort(deviations)
        winner = candidates[int(order[0])]
        best_name: str | None = winner.name
        best_deviation = float(deviations[int(order[0])])
        margin = (
            float(deviations[int(order[1])] - deviations[int(order[0])])
            if deviations.size > 1
            else float("inf")
        )
        # The winning relationship already knows which families define it, so
        # its own statement is reported rather than an arbitrary equally-exact
        # low-index alternative. The deviations then verify the statement
        # against the *fitted* rotation instead of asserting it.
        if best_deviation <= catalog_tolerance_deg:
            preferred_planes = tuple(
                np.asarray(pair[0].miller.indices, dtype=np.int64)
                for pair in winner.parallel_planes
            )
            preferred_directions = tuple(
                np.rint(np.asarray(pair[0].coordinates, dtype=np.float64)).astype(np.int64)
                for pair in winner.parallel_directions
            )
    else:
        best_name = None
        best_deviation = float("inf")
        margin = float("inf")
    planes, directions = describe_orientation_relationship(
        fitted,
        tolerance_deg=parallelism_tolerance_deg,
        max_index=max_index,
        max_statements=max_statements,
        preferred_parent_planes=preferred_planes,
        preferred_parent_directions=preferred_directions,
    )
    return ORCharacterizationReport(
        relationship=fitted,
        pair_count=len(parent_orientations),
        residuals_deg=residuals,
        iterations=iterations,
        converged=converged,
        catalog_names=tuple(candidate.name for candidate in candidates),
        catalog_deviations_deg=deviations,
        best_catalog_name=best_name,
        best_catalog_deviation_deg=best_deviation,
        margin_deg=margin,
        catalog_tolerance_deg=float(catalog_tolerance_deg),
        plane_statements=planes,
        direction_statements=directions,
        provenance=provenance,
    )


def orientation_relationship_from_euler(
    parent_euler_deg: ArrayLike,
    child_euler_deg: ArrayLike,
    *,
    parent_phase: Phase,
    child_phase: Phase,
    specimen_frame: ReferenceFrame | None = None,
    convention: str = "bunge",
    catalog: OrientationRelationshipCatalog | Sequence[OrientationRelationship] | None = None,
    nominal: OrientationRelationship | None = None,
    catalog_tolerance_deg: float = DEFAULT_OR_TOLERANCE_DEG,
    parallelism_tolerance_deg: float = DEFAULT_OR_TOLERANCE_DEG,
    max_index: int = 3,
    max_statements: int = 4,
    provenance: ProvenanceRecord | None = None,
) -> ORCharacterizationReport:
    """Determine the OR directly from measured Euler angle triples.

    Purpose: the ergonomic entry point for the common case where the
    measurements are two columns of Euler angles exported from an EBSD system,
    not `OrientationSet` objects. It builds the orientation sets in the given
    convention and defers entirely to `characterize_orientation_relationship`,
    whose options it mirrors.

    Inputs: ``(n, 3)`` parent and child Euler angle arrays **in degrees**, one
    row per parent/child pair and in matching row order; the two phases; and
    optionally the specimen frame the angles refer to (the shared standard
    specimen frame by default) and the Euler ``convention`` (Bunge ZXZ by
    default). The remaining keywords are those of
    `characterize_orientation_relationship`.

    Output: an `ORCharacterizationReport` — read its ``describe()``.
    """

    from pytex.core.frame_catalog import specimen_frame as standard_specimen_frame

    frame = standard_specimen_frame() if specimen_frame is None else specimen_frame
    parent_angles = np.asarray(parent_euler_deg, dtype=np.float64).reshape(-1, 3)
    child_angles = np.asarray(child_euler_deg, dtype=np.float64).reshape(-1, 3)
    parents = OrientationSet.from_euler_angles(
        parent_angles, specimen_frame=frame, phase=parent_phase, convention=convention
    )
    children = OrientationSet.from_euler_angles(
        child_angles, specimen_frame=frame, phase=child_phase, convention=convention
    )
    return characterize_orientation_relationship(
        parents,
        children,
        catalog=catalog,
        nominal=nominal,
        catalog_tolerance_deg=catalog_tolerance_deg,
        parallelism_tolerance_deg=parallelism_tolerance_deg,
        max_index=max_index,
        max_statements=max_statements,
        provenance=provenance,
    )


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
    correspondence_denominator: int = 1
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        gradient = np.asarray(self.deformation_gradient, dtype=np.float64)
        stretch = np.asarray(self.stretch_tensor, dtype=np.float64)
        stretches = np.asarray(self.principal_stretches, dtype=np.float64).reshape(-1)
        directions = np.asarray(self.principal_directions, dtype=np.float64)
        # Float, not int: a reconstructive correspondence is genuinely
        # half-integer (the Burgers shuffle), and coercing to int would silence
        # that by rounding the halves away into a singular matrix.
        correspondence = np.asarray(self.correspondence, dtype=np.float64)
        if gradient.shape != (3, 3) or stretch.shape != (3, 3) or correspondence.shape != (3, 3):
            raise ValueError("Deformation matrices must have shape (3, 3).")
        if stretches.shape != (3,) or directions.shape != (3, 3):
            raise ValueError("Principal quantities must have three entries.")
        if np.any(stretches <= 0.0) or self.volume_ratio <= 0.0:
            raise ValueError("Principal stretches and volume ratio must be positive.")
        if self.polar_rotation_deg < 0.0 or self.correspondence_max_component_error < 0.0:
            raise ValueError("Polar rotation and component error must be non-negative.")
        if self.correspondence_denominator < 1:
            raise ValueError("correspondence_denominator must be a positive integer.")
        if np.isclose(float(np.linalg.det(correspondence)), 0.0, atol=1e-9):
            raise ValueError("The lattice correspondence must be non-singular.")
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
        denominator = self.correspondence_denominator
        scaled = np.rint(np.asarray(self.correspondence) * denominator).astype(np.int64)
        suffix = "" if denominator == 1 else f"/{denominator}"
        columns = ", ".join(
            format_direction_indices(_index_tuple(scaled[:, i]), style="plain") + suffix
            for i in range(3)
        )
        if denominator == 1:
            correspondence_text = (
                f"Integer lattice correspondence maps the parent basis to {columns} "
                f"(child indices; largest entry-wise deviation "
                f"{self.correspondence_max_component_error:.3f})."
            )
        else:
            correspondence_text = (
                f"The lattice correspondence is not integer but has denominator "
                f"{denominator}: it maps the parent basis to {columns} (child indices; "
                f"largest entry-wise deviation "
                f"{self.correspondence_max_component_error:.3f}). A non-unit denominator "
                "means the parent lattice does not map onto the child lattice by strain "
                "alone — an atomic shuffle carries the remainder, as in the "
                "reconstructive Burgers bcc-to-hcp transformation."
            )
        return (
            f"Transformation deformation for '{self.relationship_name}'{variant_text}: "
            f"principal strains {strains[0]:+.2f}%, {strains[1]:+.2f}%, {strains[2]:+.2f}% "
            f"(principal directions in the parent crystal frame), volume change "
            f"{(self.volume_ratio - 1.0) * 100.0:+.2f}%. {correspondence_text} "
            f"Residual polar rotation {self.polar_rotation_deg:.2f} deg — zero for the "
            "pure correspondence distortion (Bain), and the rigid-body rotation relative "
            "to it for KS-class relationships."
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
        """Number of distinct variants contributing poles to this figure.
        """

        return int(np.unique(self.variant_indices).size)

    def describe(self) -> str:
        """Prose summary: plane family, variants, and pole counts."""

        per_variant = int(
            np.count_nonzero(self.variant_indices == int(self.variant_indices[0]))
        )
        # A pole figure plots the whole symmetry-related orbit, so the family
        # brackets {hkl} are the correct notation here, not (hkl).
        plane_text = format_plane_family_indices(
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
class VariantCorrespondenceRow:
    """One (variant, source object, image) entry of a variant correspondence table.

    ``exact_components`` are the generally irrational basis components of the
    image; ``indices`` is the nearest primitive integer triple and
    ``residual_deg`` the angle between the two. A zero residual means the image
    really is that low-index object; a large one means the source has no
    low-index image in this variant, which is itself the answer.

    ``equivalence_group`` labels variants that produce crystallographically
    equivalent images of the *same* source object — members of one index family.
    """

    variant_index: int
    source_indices: np.ndarray
    exact_components: np.ndarray
    indices: np.ndarray
    residual_deg: float
    source_label: str
    image_label: str
    equivalence_group: int

    def __post_init__(self) -> None:
        source = np.asarray(self.source_indices, dtype=np.int64)
        indices = np.asarray(self.indices, dtype=np.int64)
        exact = np.asarray(self.exact_components, dtype=np.float64)
        if source.shape != (3,) or indices.shape != (3,) or exact.shape != (3,):
            raise ValueError("VariantCorrespondenceRow arrays must have shape (3,).")
        if self.variant_index <= 0:
            raise ValueError("variant_index must be positive.")
        if not np.isfinite(self.residual_deg) or self.residual_deg < 0.0:
            raise ValueError("residual_deg must be finite and non-negative.")
        if self.equivalence_group < 0:
            raise ValueError("equivalence_group must be non-negative.")
        for array in (source, indices, exact):
            array.setflags(write=False)
        object.__setattr__(self, "source_indices", source)
        object.__setattr__(self, "indices", indices)
        object.__setattr__(self, "exact_components", exact)


@dataclass(frozen=True, slots=True)
class VariantCorrespondenceTable:
    """What arbitrary parent planes or directions become in every product variant.

    The tabular answer to "I have this ``(hkl)`` (or ``[uvw]``) in the parent —
    what is the parallel plane (direction) in each of the product variants?".
    One row per (source object, variant), carrying the exact image, its nearest
    integer indices, the angular residual between them, and an equivalence-group
    label that collapses variants giving crystallographically equivalent images.

    The grouping is what makes a 24-row table readable. Under Kurdjumov-Sachs
    the austenite ``(111)`` produces exactly four distinct answers across the 24
    variants, six variants each: the ``{011}`` close-packed image at zero
    residual (the packet the variant belongs to) and three higher-index images.
    """

    relationship_name: str
    kind: str
    sense: str
    source_phase_name: str
    image_phase_name: str
    rows: tuple[VariantCorrespondenceRow, ...]
    max_index: int
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.kind not in {"plane", "direction"}:
            raise ValueError("VariantCorrespondenceTable.kind must be 'plane' or 'direction'.")
        if self.sense not in {"parent_to_child", "child_to_parent"}:
            raise ValueError(
                "VariantCorrespondenceTable.sense must be 'parent_to_child' or "
                "'child_to_parent'."
            )
        rows = tuple(self.rows)
        if not rows:
            raise ValueError("VariantCorrespondenceTable requires at least one row.")
        if self.max_index < 1:
            raise ValueError("max_index must be at least 1.")
        object.__setattr__(self, "rows", rows)

    @property
    def source_indices(self) -> tuple[tuple[int, int, int], ...]:
        """The distinct source objects tabulated, in the order supplied."""

        seen: list[tuple[int, int, int]] = []
        for row in self.rows:
            key = _index_tuple(row.source_indices)
            if key not in seen:
                seen.append(key)
        return tuple(seen)

    @property
    def variant_indices(self) -> tuple[int, ...]:
        """The variant indices present in the table, in first-appearance order.
        """

        seen: list[int] = []
        for row in self.rows:
            if row.variant_index not in seen:
                seen.append(row.variant_index)
        return tuple(seen)

    def rows_for(self, source_indices: ArrayLike) -> tuple[VariantCorrespondenceRow, ...]:
        """Every variant's image of one source object, in variant order."""

        key = _index_tuple(np.rint(np.asarray(source_indices, dtype=np.float64)))
        return tuple(row for row in self.rows if _index_tuple(row.source_indices) == key)

    def rows_for_variant(self, variant_index: int) -> tuple[VariantCorrespondenceRow, ...]:
        """Every source object's image in one variant, in source order."""

        return tuple(row for row in self.rows if row.variant_index == int(variant_index))

    def distinct_image_count(self, source_indices: ArrayLike) -> int:
        """How many crystallographically distinct images one source object has.

        The number of index families the variants map it onto — four for the
        Kurdjumov-Sachs ``(111)`` across all 24 variants.
        """

        rows = self.rows_for(source_indices)
        if not rows:
            raise KeyError("No rows for the requested source indices.")
        return len({row.equivalence_group for row in rows})

    def exact_rows(
        self, *, tolerance_deg: float = 1e-6
    ) -> tuple[VariantCorrespondenceRow, ...]:
        """Rows whose image really is the low-index object reported.

        The physically interesting subset: under Kurdjumov-Sachs these are the
        six variants that carry ``(111)`` austenite onto a ``{011}`` ferrite
        plane, i.e. the close-packed (packet) correspondence.
        """

        return tuple(row for row in self.rows if row.residual_deg <= tolerance_deg)

    def to_records(self) -> list[dict[str, object]]:
        """One flat dictionary per row, suitable for a DataFrame or CSV writer."""

        return [
            {
                "relationship": self.relationship_name,
                "kind": self.kind,
                "sense": self.sense,
                "variant": row.variant_index,
                "source_phase": self.source_phase_name,
                "source_indices": " ".join(
                    str(value) for value in _index_tuple(row.source_indices)
                ),
                "source_label": row.source_label,
                "image_phase": self.image_phase_name,
                "image_indices": " ".join(str(value) for value in _index_tuple(row.indices)),
                "image_label": row.image_label,
                "exact_components": " ".join(
                    f"{value:.12g}" for value in row.exact_components
                ),
                "residual_deg": row.residual_deg,
                "equivalence_group": row.equivalence_group,
            }
            for row in self.rows
        ]

    def to_csv(self, path: str | Path) -> Path:
        """Write the table as UTF-8 CSV with a header row and return the path."""

        records = self.to_records()
        output = Path(path)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        return output

    def to_markdown(self) -> str:
        """Render the table as GitHub-flavoured Markdown for reports and notebooks."""

        header = (
            "| variant | source | image | residual (deg) | group |\n"
            "| --- | --- | --- | --- | --- |\n"
        )
        body = "".join(
            f"| {row.variant_index} | {row.source_label} | {row.image_label} | "
            f"{row.residual_deg:.4f} | {row.equivalence_group} |\n"
            for row in self.rows
        )
        return header + body

    def to_json_dict(self) -> dict[str, object]:
        """A JSON-ready payload carrying the same facts as ``describe()``.

        A one-way report payload for manifests and downstream tooling, not a
        round-trip contract; use `pytex.contracts` for objects that must be
        reconstructed.
        """

        return {
            "schema": "pytex.variant_correspondence_table/1",
            "relationship": self.relationship_name,
            "kind": self.kind,
            "sense": self.sense,
            "source_phase": self.source_phase_name,
            "image_phase": self.image_phase_name,
            "max_index": int(self.max_index),
            "rows": self.to_records(),
        }

    def describe(self) -> str:
        """Prose summary: what was mapped, how many distinct answers, and which are exact."""

        noun = "plane" if self.kind == "plane" else "direction"
        origin, destination = (
            ("parent", "child") if self.sense == "parent_to_child" else ("child", "parent")
        )
        lines = [
            f"Variant correspondence table for orientation relationship "
            f"'{self.relationship_name}': {len(self.source_indices)} {origin} "
            f"{noun}(s) mapped through {len(self.variant_indices)} variant(s), "
            f"{len(self.rows)} row(s). Images are given in the {destination} phase "
            f"('{self.image_phase_name}') as exact basis components plus the nearest "
            f"primitive integer indices within max_index {self.max_index}; the residual "
            "is the angle between the two, so a zero residual means the image really is "
            "that low-index object."
        ]
        for source in self.source_indices:
            rows = self.rows_for(source)
            exact = [row for row in rows if row.residual_deg <= 1e-6]
            label = rows[0].source_label
            summary = (
                f"  {label}: {len({row.equivalence_group for row in rows})} "
                f"crystallographically distinct image(s) across {len(rows)} variant(s)"
            )
            if exact:
                variants = ", ".join(str(row.variant_index) for row in exact[:12])
                if len(exact) > 12:
                    variants += ", ..."
                summary += (
                    f"; exactly parallel in {len(exact)} of them "
                    f"(variants {variants}, image {exact[0].image_label})"
                )
            else:
                summary += (
                    "; no variant carries it onto a low-index object exactly, so every "
                    "image is a rationalization with a stated residual"
                )
            lines.append(summary + ".")
        lines.append(
            "  Exactly-parallel images are a property of the relationship and do not depend "
            "on max_index; how the remaining irrational images are grouped does, because a "
            "larger bound splits them across more index families."
        )
        return "\n".join(lines)


def variant_correspondence_table(
    relationship: OrientationRelationship,
    objects: CrystalPlane | CrystalDirection | Sequence[CrystalPlane | CrystalDirection],
    *,
    sense: str = "parent_to_child",
    variants: tuple[TransformationVariant, ...] | None = None,
    max_index: int = DEFAULT_RATIONALIZATION_MAX_INDEX,
    provenance: ProvenanceRecord | None = None,
) -> VariantCorrespondenceTable:
    """Map arbitrary planes or directions through every transformation variant.

    Purpose: the everyday variant-resolved question — given an orientation
    relationship and any parent plane ``(hkl)`` or direction ``[uvw]``, what is
    the parallel plane or direction in each product variant? Used for trace
    analysis, for predicting which variants share a habit or diffraction
    feature, and for teaching what "24 variants" actually means.

    When to use: whenever the answer wanted is a *table* rather than a single
    mapping. `OrientationRelationship.map_plane_to_child` answers for one
    variant; `map_plane_across_variants` returns the raw correspondences; this
    surface adds the grouping, labels, residuals, `describe()` and the CSV,
    Markdown and JSON exports that make the table usable in a report.

    Inputs: the relationship; one `CrystalPlane` or `CrystalDirection`, or a
    sequence of them (all of the same kind and the same phase); ``sense``,
    either ``"parent_to_child"`` (the default; objects belong to the parent
    phase) or ``"child_to_parent"``; an optional explicit ``variants`` tuple
    (default: all, in `generate_variants()` order); and ``max_index``, the bound
    on the rationalized image indices.

    Output: a `VariantCorrespondenceTable` — read its ``describe()``.

    Examples
    --------
    Under Kurdjumov-Sachs the austenite ``(111)`` has four crystallographically
    distinct images across the 24 variants, six variants each. In six of them
    the image is a ``{011}`` ferrite plane at zero residual — the close-packed
    correspondence that defines the packet — and the other eighteen land on
    higher-index rationalizations with residuals of 0.36, 2.48 and 3.69 deg.

    See also
    --------
    `map_plane_across_variants`, `map_direction_across_variants` : the raw
        per-variant correspondences without tabulation.
    `variant_close_packed_groups` : packet labels from the defining plane family.
    `find_parallel_planes` : the search over a whole symmetry family rather than
        one nominated object.
    """

    if sense not in {"parent_to_child", "child_to_parent"}:
        raise ValueError("sense must be 'parent_to_child' or 'child_to_parent'.")
    items: list[CrystalPlane | CrystalDirection] = (
        [objects]
        if isinstance(objects, CrystalPlane | CrystalDirection)
        else list(objects)
    )
    if not items:
        raise ValueError("variant_correspondence_table requires at least one object.")
    kinds = {"plane" if isinstance(item, CrystalPlane) else "direction" for item in items}
    if len(kinds) != 1:
        raise ValueError(
            "variant_correspondence_table requires all objects to be of one kind; "
            "call it once for planes and once for directions."
        )
    kind = kinds.pop()
    reciprocal = kind == "plane"

    forward = sense == "parent_to_child"
    source_phase = relationship.parent_phase if forward else relationship.child_phase
    image_phase = relationship.child_phase if forward else relationship.parent_phase
    for item in items:
        if not phases_semantically_match(item.phase, source_phase):
            raise ValueError(
                f"Every object must belong to the {'parent' if forward else 'child'} "
                f"phase '{source_phase.name}' for sense '{sense}'."
            )

    resolved = relationship.generate_variants() if variants is None else tuple(variants)
    if not resolved:
        raise ValueError("variant_correspondence_table requires at least one variant.")

    rows: list[VariantCorrespondenceRow] = []
    for item in items:
        source_indices = (
            np.asarray(item.miller.indices, dtype=np.int64)
            if isinstance(item, CrystalPlane)
            else np.rint(np.asarray(item.coordinates, dtype=np.float64)).astype(np.int64)
        )
        source_label = _crystallographic_label(
            source_indices, phase=source_phase, reciprocal=reciprocal
        )
        # One equivalence-group numbering per source object: the question
        # "how many distinct answers does *this* object have" is per object.
        group_labels: dict[tuple[int, int, int], int] = {}
        for variant in resolved:
            result: PlaneCorrespondence | DirectionCorrespondence
            if isinstance(item, CrystalPlane):
                plane_result = (
                    relationship.map_plane_to_child(item, variant=variant, max_index=max_index)
                    if forward
                    else relationship.map_plane_to_parent(
                        item, variant=variant, max_index=max_index
                    )
                )
                result = plane_result
                exact = plane_result.target_exact_indices
            else:
                direction_result = (
                    relationship.map_direction_to_child(
                        item, variant=variant, max_index=max_index
                    )
                    if forward
                    else relationship.map_direction_to_parent(
                        item, variant=variant, max_index=max_index
                    )
                )
                result = direction_result
                exact = direction_result.target_exact_coordinates
            indices = np.asarray(result.rational_indices, dtype=np.int64)
            family = _family_key(indices, phase=image_phase, reciprocal=reciprocal)
            group = group_labels.setdefault(family, len(group_labels))
            rows.append(
                VariantCorrespondenceRow(
                    variant_index=variant.variant_index,
                    source_indices=source_indices,
                    exact_components=np.asarray(exact, dtype=np.float64),
                    indices=indices,
                    residual_deg=float(result.angular_residual_deg),
                    source_label=source_label,
                    image_label=_crystallographic_label(
                        indices, phase=image_phase, reciprocal=reciprocal
                    ),
                    equivalence_group=group,
                )
            )
    return VariantCorrespondenceTable(
        relationship_name=relationship.name,
        kind=kind,
        sense=sense,
        source_phase_name=source_phase.name,
        image_phase_name=image_phase.name,
        rows=tuple(rows),
        max_index=max_index,
        provenance=provenance or relationship.provenance,
    )


@dataclass(frozen=True, slots=True)
class PhaseTransformationRecord:
    """A measured parent orientation with its child orientations and their OR.

    Purpose
    -------
    The unit of transformation-crystallography analysis: what transformed
    into what, under which relationship, through which variants. It is the
    input to variant-selection statistics, to deviation analysis, and to
    parent-grain reconstruction.

    Construction cross-checks that the parent and child orientations belong
    to the phases the relationship declares, so a record cannot silently pair
    data with the wrong relationship.

    Attributes
    ----------
    name : str
        Non-empty identifier.
    orientation_relationship : OrientationRelationship
    parent_orientation : Orientation
        Its phase must match the relationship's parent phase.
    child_orientations : OrientationSet
        Their phase must match the relationship's child phase.
    variant_indices : np.ndarray, optional
        One variant index per child, when the assignment has been made.
        ``None`` means unassigned, not "all variant 1".
    notes : tuple of str
        Free-text remarks, including honest limitations of the record.
    provenance : ProvenanceRecord, optional
    """

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
        """Number of distinct transformation variants referenced by this record.

        Zero when no variant assignment was made.
        """

        if self.variant_indices is None:
            return 0
        return int(np.unique(self.variant_indices).size)

    def predicted_child_orientations(self) -> OrientationSet:
        """The child orientations this record's OR and variants predict.

        Purpose
        -------
        The forward model of the transformation: apply each assigned variant to
        its parent orientation and report the child orientation it implies.
        Comparing these against the measured child orientations is what makes
        a variant assignment falsifiable.

        Convention
        ----------
        Composed as ``g_child = g_parent . V^T`` in the crystal-to-specimen
        convention. The transpose is not cosmetic: composing as ``V . g_parent``
        is a different, wrong prediction that self-consistent synthetic tests
        cannot detect but real measured orientations expose immediately.

        Returns
        -------
        OrientationSet
            One predicted orientation per child, in record order. When no
            variant indices are present the base relationship is applied to
            every pair.
        """

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
