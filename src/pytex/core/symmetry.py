from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import (
    as_float_array,
    is_rotation_matrix,
    normalize_vector,
    normalize_vectors,
)
from pytex.core.batches import VectorSet
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.point_groups import (
    PointGroup,
    laue_class_symbol_for,
    normalize_point_group_symbol,
    proper_subgroup_symbol_for,
)
from pytex.core.provenance import ProvenanceRecord

_SPECIMEN_SYMMETRY_POINT_GROUPS = {
    "triclinic": "1",
    "monoclinic": "2",
    "orthorhombic": "222",
    "orthotropic": "222",
}


def _rotation_matrix_from_axis_angle(axis: ArrayLike, angle_deg: float) -> np.ndarray:
    unit_axis = normalize_vector(axis)
    angle_rad = np.deg2rad(angle_deg)
    x, y, z = unit_axis
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    one_minus_c = 1.0 - c
    matrix = np.array(
        [
            [c + x * x * one_minus_c, x * y * one_minus_c - z * s, x * z * one_minus_c + y * s],
            [y * x * one_minus_c + z * s, c + y * y * one_minus_c, y * z * one_minus_c - x * s],
            [z * x * one_minus_c - y * s, z * y * one_minus_c + x * s, c + z * z * one_minus_c],
        ],
        dtype=np.float64,
    )
    return as_float_array(matrix, shape=(3, 3))


def _matrix_key(matrix: np.ndarray) -> tuple[float, ...]:
    rounded = np.round(matrix, decimals=8)
    return tuple(float(value) for value in rounded.ravel())


def _unique_rotation_matrices(matrices: list[np.ndarray]) -> np.ndarray:
    unique: dict[tuple[float, ...], np.ndarray] = {}
    for matrix in matrices:
        if not is_rotation_matrix(matrix):
            raise ValueError("Generated symmetry operator is not a proper rotation matrix.")
        unique[_matrix_key(matrix)] = as_float_array(matrix, shape=(3, 3))
    operators = np.stack(list(unique.values()), axis=0)
    operators = np.ascontiguousarray(operators)
    operators.setflags(write=False)
    return operators


def _group_from_generators(generators: list[np.ndarray]) -> np.ndarray:
    identity = np.eye(3, dtype=np.float64)
    known = {_matrix_key(identity): identity}
    frontier = [identity]
    all_generators = [as_float_array(generator, shape=(3, 3)) for generator in generators]
    while frontier:
        current = frontier.pop()
        for generator in all_generators:
            for candidate in (current @ generator, generator @ current):
                key = _matrix_key(candidate)
                if key not in known:
                    known[key] = candidate
                    frontier.append(candidate)
    return _unique_rotation_matrices(list(known.values()))


def _point_group_generators() -> dict[str, list[np.ndarray]]:
    return {
        "1": [],
        "2": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 180.0)],
        "222": [
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
            _rotation_matrix_from_axis_angle([0.0, 1.0, 0.0], 180.0),
        ],
        "4": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0)],
        "422": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "3": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 120.0)],
        "32": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 120.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "6": [_rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 60.0)],
        "622": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 60.0),
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
        ],
        "23": [
            _rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 180.0),
            _rotation_matrix_from_axis_angle([1.0, 1.0, 1.0], 120.0),
        ],
        "432": [
            _rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 90.0),
            _rotation_matrix_from_axis_angle([1.0, 1.0, 1.0], 120.0),
        ],
    }


_SECTOR_TOLERANCE = 1e-8


def _normalized_proper_point_group(point_group: str) -> str:
    return proper_subgroup_symbol_for(normalize_point_group_symbol(point_group))


def _sector_array(rows: list[list[float]]) -> np.ndarray:
    if not rows:
        return as_float_array(np.zeros((0, 3), dtype=np.float64), shape=(None, 3))
    return normalize_vectors(np.array(rows, dtype=np.float64))


@cache
def _sector_geometry_table() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    # Exact Laue-reduced (antipodal, upper hemisphere) sector polygons per
    # proper rotation group. Vertices are ordered loops; edge i is the great
    # circle arc from vertex i to vertex i+1 with the given inward normal.
    # Sector area is 4*pi / (2 * group order).
    half = 0.5
    root3_half = np.sqrt(3.0) / 2.0
    z_axis = [0.0, 0.0, 1.0]
    x_axis = [1.0, 0.0, 0.0]
    y_axis = [0.0, 1.0, 0.0]
    y_normal = [0.0, 1.0, 0.0]
    octant_vertices = _sector_array([z_axis, x_axis, [0.0, 1.0, 0.0]])
    octant_normals = _sector_array([y_normal, z_axis, x_axis])
    return {
        "1": (_sector_array([]), _sector_array([z_axis])),
        "2": (
            _sector_array([x_axis, [-1.0, 0.0, 0.0]]),
            _sector_array([z_axis, y_normal]),
        ),
        "222": (octant_vertices, octant_normals),
        "4": (octant_vertices, octant_normals),
        "422": (
            _sector_array([z_axis, x_axis, [1.0, 1.0, 0.0]]),
            _sector_array([y_normal, z_axis, [1.0, -1.0, 0.0]]),
        ),
        "3": (
            _sector_array([z_axis, x_axis, [-half, root3_half, 0.0]]),
            _sector_array([y_normal, z_axis, [root3_half, half, 0.0]]),
        ),
        # For D3 the antipodal-composed two-fold axes at azimuths 0/120/240
        # induce mirror lines at azimuths 30/90/150, so the fundamental wedge
        # is [30, 90] degrees rather than a wedge anchored at azimuth 0.
        "32": (
            _sector_array([z_axis, [root3_half, half, 0.0], y_axis]),
            _sector_array([[-half, root3_half, 0.0], z_axis, x_axis]),
        ),
        "6": (
            _sector_array([z_axis, x_axis, [half, root3_half, 0.0]]),
            _sector_array([y_normal, z_axis, [root3_half, -half, 0.0]]),
        ),
        "622": (
            _sector_array([z_axis, x_axis, [root3_half, half, 0.0]]),
            _sector_array([y_normal, z_axis, [half, -root3_half, 0.0]]),
        ),
        "23": (
            _sector_array([z_axis, x_axis, [1.0, 1.0, 1.0]]),
            _sector_array([y_normal, [0.0, -1.0, 1.0], [1.0, -1.0, 0.0]]),
        ),
        "432": (
            _sector_array([z_axis, [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]),
            _sector_array([y_normal, [-1.0, 0.0, 1.0], [1.0, -1.0, 0.0]]),
        ),
    }


def _sector_vertices_for_group(proper_group: str) -> np.ndarray:
    return _sector_geometry_table()[proper_group][0]


def _sector_edge_normals_for_group(proper_group: str) -> np.ndarray:
    return _sector_geometry_table()[proper_group][1]


def _vector_in_fundamental_sector(vector: np.ndarray, proper_group: str) -> bool:
    normals = _sector_edge_normals_for_group(proper_group)
    return bool(np.all(normals @ np.asarray(vector, dtype=np.float64) >= -_SECTOR_TOLERANCE))


def _sector_sort_key(vector: np.ndarray) -> tuple[float, float, float]:
    rounded = np.round(vector, decimals=12)
    return (
        float(rounded[2]),
        float(rounded[0]),
        float(rounded[1]),
    )


@cache
def _operators_for_proper_point_group(proper_group: str) -> np.ndarray:
    operators = _group_from_generators(_point_group_generators()[proper_group])
    operators.setflags(write=False)
    return operators


def _canonical_vector_index(vectors: np.ndarray) -> int:
    rounded = np.round(vectors, decimals=12)
    return int(np.lexsort((rounded[:, 1], rounded[:, 0], rounded[:, 2]))[-1])


@dataclass(frozen=True, slots=True, eq=False)
class SymmetrySpec:
    """A crystal or specimen symmetry group, as explicit rotation operators.

    Purpose
    -------
    The symmetry authority of the library. Every symmetry-aware operation —
    orientation reduction, disorientation, pole-figure families, inverse
    pole figures, variant enumeration — consumes this one type, so no
    subsystem defines its own symmetry model.

    The stored operators are the *proper rotations* of the group: an
    orientation is a rotation, so improper operators cannot map one
    orientation onto another. Use :meth:`to_point_group` when mirrors,
    inversion, and rotoinversions are needed as well.

    Attributes
    ----------
    name : str
    point_group : str
        Hermann-Mauguin symbol. May name a centrosymmetric or Laue group;
        the operators stored are those of the corresponding proper group.
    operators : np.ndarray
        ``(order, 3, 3)`` rotation matrices.
    specimen_symmetry : str, optional
        Name of an associated statistical specimen symmetry.
    reference_frame : ReferenceFrame, optional
        The frame the operators act in. Strongly recommended: without it,
        vectors handed to this specification cannot be frame-checked.
    provenance : ProvenanceRecord, optional
    """

    name: str
    point_group: str
    operators: np.ndarray = field(default_factory=lambda: np.eye(3)[None, :, :])
    specimen_symmetry: str | None = None
    reference_frame: ReferenceFrame | None = None
    provenance: ProvenanceRecord | None = None

    def __eq__(self, other: object) -> bool:
        # The generated dataclass __eq__ would compare the operator ndarray with
        # `==` and raise on distinct-but-equal instances; equality here is the
        # symmetry identity (provenance excluded).
        if not isinstance(other, SymmetrySpec):
            return NotImplemented
        return (
            self.name == other.name
            and self.point_group == other.point_group
            and self.specimen_symmetry == other.specimen_symmetry
            and self.reference_frame == other.reference_frame
            and bool(np.array_equal(self.operators, other.operators))
        )

    def __hash__(self) -> int:
        return hash((self.name, self.point_group, self.specimen_symmetry))

    def __post_init__(self) -> None:
        operators = as_float_array(self.operators, shape=(None, 3, 3))
        for operator in operators:
            if not is_rotation_matrix(operator):
                raise ValueError("All symmetry operators must be proper rotation matrices.")
        if self.reference_frame is not None and self.reference_frame.domain not in {
            FrameDomain.CRYSTAL,
            FrameDomain.SPECIMEN,
        }:
            raise ValueError(
                "SymmetrySpec.reference_frame must belong to the crystal or specimen domain."
            )
        object.__setattr__(self, "operators", operators)

    @property
    def order(self) -> int:
        """Number of symmetry operators in the group.

        24 for cubic ``432``, 12 for hexagonal ``622``, 1 for triclinic. This is
        the factor by which the orientation fundamental region shrinks, and the
        multiplier on the cost of every symmetry-aware reduction.
        """

        return int(self.operators.shape[0])

    @property
    def proper_point_group(self) -> str:
        """The rotation (proper) point group corresponding to this symmetry.

        Diffraction and orientation work is governed by the rotation subgroup:
        an orientation is a rotation, so improper operators cannot map one
        orientation onto another. This normalizes any Hermann-Mauguin symbol —
        including centrosymmetric and Laue symbols — to the proper group whose
        operators are actually stored.
        """

        return _normalized_proper_point_group(self.point_group)

    @property
    def laue_group_symbol(self) -> str:
        """Hermann-Mauguin symbol of the Laue class of this point group.

        The Laue class is the point group with inversion added. It is what a
        diffraction experiment can determine, because Friedel's law makes
        ``+g`` and ``-g`` equal in intensity under kinematic scattering.
        """

        return laue_class_symbol_for(self.point_group)

    @property
    def is_laue(self) -> bool:
        """Whether this symmetry is already a Laue (centrosymmetric) group.
        """

        return normalize_point_group_symbol(self.point_group) == self.laue_group_symbol

    def to_point_group(self) -> PointGroup:
        """The full :class:`~pytex.core.point_groups.PointGroup` for this symmetry.

        ``SymmetrySpec`` carries the proper operators used for orientation
        algebra; the ``PointGroup`` additionally carries improper operators,
        mirror normals, Schoenflies naming, and crystal-system membership.
        """

        return PointGroup.from_symbol(self.point_group)

    def laue_symmetry(self) -> SymmetrySpec:
        """The Laue-class symmetry corresponding to this one.

        Use it when a calculation must respect only what diffraction can
        distinguish — pole figures and inverse pole figures being the standard
        cases.
        """

        return SymmetrySpec.from_point_group(
            self.laue_group_symbol,
            reference_frame=self.reference_frame,
            specimen_symmetry=self.specimen_symmetry,
            provenance=self.provenance,
        )

    @classmethod
    def identity(
        cls,
        *,
        name: str = "identity",
        point_group: str = "1",
        reference_frame: ReferenceFrame | None = None,
    ) -> SymmetrySpec:
        """The trivial symmetry: one operator, point group ``1``.

        The correct explicit choice for a triclinic phase and for "no specimen
        symmetry assumed", which is preferable to leaving symmetry unset.
        """

        return cls(
            name=name,
            point_group=point_group,
            operators=np.eye(3)[None, :, :],
            reference_frame=reference_frame,
        )

    @classmethod
    def from_point_group(
        cls,
        point_group: str,
        *,
        reference_frame: ReferenceFrame | None = None,
        specimen_symmetry: str | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> SymmetrySpec:
        """Build a symmetry specification from a Hermann-Mauguin symbol.

        Purpose
        -------
        The standard constructor: it expands the symbol into the explicit
        rotation-operator array that every symmetry-aware calculation consumes.

        Parameters
        ----------
        point_group : str
            Hermann-Mauguin symbol, for example ``"m-3m"``, ``"432"``,
            ``"6/mmm"``. The stored operators are those of the corresponding
            proper group.
        reference_frame : ReferenceFrame, optional
            The frame the operators act in. Strongly recommended: without it,
            the specification cannot verify that vectors handed to it live in
            the right frame.
        specimen_symmetry : str, optional
            Name of an associated statistical specimen symmetry.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        SymmetrySpec
        """

        proper_group = _normalized_proper_point_group(point_group)
        operators = _operators_for_proper_point_group(proper_group)
        return cls(
            name=proper_group,
            point_group=point_group,
            operators=operators,
            specimen_symmetry=specimen_symmetry,
            reference_frame=reference_frame,
            provenance=provenance,
        )

    @classmethod
    def specimen(
        cls,
        name: str = "triclinic",
        *,
        reference_frame: ReferenceFrame | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> SymmetrySpec:
        """Build a statistical specimen (sample) symmetry by name.

        Purpose
        -------
        Specimen symmetry expresses an assumption about the *process*, not about
        the crystal: rolling is usually taken to impose orthorhombic sample
        symmetry, axisymmetric drawing to impose a fibre symmetry. Imposing it
        averages the ODF over those operations, which sharpens statistics when
        the assumption holds and fabricates symmetry when it does not.

        Parameters
        ----------
        name : str
            Specimen symmetry name; ``"triclinic"`` (no assumption) by default.
            An unsupported name raises and lists the supported ones.
        reference_frame : ReferenceFrame, optional
            Must be a specimen-domain frame.
        provenance : ProvenanceRecord, optional
        """

        normalized = name.strip().lower()
        point_group = _SPECIMEN_SYMMETRY_POINT_GROUPS.get(normalized)
        if point_group is None:
            supported = ", ".join(sorted(_SPECIMEN_SYMMETRY_POINT_GROUPS))
            raise ValueError(
                f"Unsupported specimen symmetry '{name}'. Supported names: {supported}."
            )
        if reference_frame is not None and reference_frame.domain != FrameDomain.SPECIMEN:
            raise ValueError(
                "SymmetrySpec.specimen reference_frame must belong to the specimen domain."
            )
        spec = cls.from_point_group(
            point_group,
            reference_frame=reference_frame,
            specimen_symmetry=normalized,
            provenance=provenance,
        )
        return cls(
            name=f"specimen-{normalized}",
            point_group=spec.point_group,
            operators=spec.operators,
            specimen_symmetry=normalized,
            reference_frame=reference_frame,
            provenance=provenance,
        )

    def apply_to_vectors(self, vectors: ArrayLike | VectorSet) -> np.ndarray:
        """Apply every symmetry operator to the given vectors.

        Parameters
        ----------
        vectors : ArrayLike or VectorSet
            Any array ending in dimension 3. A ``VectorSet`` must carry this
            specification's reference frame, which is checked.

        Returns
        -------
        np.ndarray
            Shape ``(order, ...)``: the operator axis is prepended, so the
            result is the full orbit rather than a reduced representative.
            Read-only.
        """

        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            array = vectors.values
        else:
            array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = np.einsum("oij,...j->o...i", self.operators, array, optimize=True)
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def apply_to_rotation_matrices(self, matrices: ArrayLike, *, side: str = "right") -> np.ndarray:
        """Apply every symmetry operator to a stack of rotation matrices.

        Purpose
        -------
        Orientations and misorientations transform by symmetry on a *side*:
        crystal symmetry acts on the right of a crystal-to-specimen orientation,
        specimen symmetry on the left. Getting the side wrong silently produces
        a different, wrong orbit, so the side is an explicit argument here.

        Parameters
        ----------
        matrices : ArrayLike
            Any array with trailing shape ``(3, 3)``.
        side : str
            ``"right"`` (default) for crystal symmetry, ``"left"`` for specimen
            symmetry.

        Returns
        -------
        np.ndarray
            Shape ``(order, ..., 3, 3)``, read-only.
        """

        rotations = np.asarray(matrices, dtype=np.float64)
        if rotations.shape[-2:] != (3, 3):
            raise ValueError("Input rotation matrices must have trailing shape (3, 3).")
        if side == "right":
            transformed = np.einsum("...ij,ojk->o...ik", rotations, self.operators, optimize=True)
        elif side == "left":
            transformed = np.einsum("oij,...jk->o...ik", self.operators, rotations, optimize=True)
        else:
            raise ValueError("side must be either 'left' or 'right'.")
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def equivalent_vectors(self, vector: ArrayLike, *, antipodal: bool = False) -> np.ndarray:
        """The distinct symmetry-equivalent unit vectors of one direction.

        Duplicates produced by operators that fix the direction are removed, so
        the returned count is the true multiplicity of the direction — 6 for
        cubic ``<100>``, 8 for ``<111>``.

        Parameters
        ----------
        vector : ArrayLike
            A single direction; normalized internally.
        antipodal : bool
            Also include the negated vectors, giving the family without a sense.

        Returns
        -------
        np.ndarray
            ``(m, 3)`` unit vectors, read-only.
        """

        candidates = self.apply_to_vectors(vector)
        candidates = normalize_vectors(candidates)
        if antipodal:
            combined = np.concatenate([candidates, -candidates], axis=0)
            candidates = normalize_vectors(combined)
        unique_vectors: dict[tuple[float, ...], np.ndarray] = {}
        for candidate in candidates:
            key = tuple(np.round(candidate, decimals=8))
            unique_vectors[key] = candidate
        array = np.stack(list(unique_vectors.values()), axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        return array

    def canonicalize_vector(self, vector: ArrayLike, *, antipodal: bool = False) -> np.ndarray:
        """The canonical representative of one direction's symmetry orbit.

        Deterministic: symmetry-equivalent directions always yield the same
        vector, so directions can be compared and grouped. This is a
        *representative*, chosen by a fixed ordering rule, and is not
        necessarily the one inside the fundamental sector — for that use
        :meth:`reduce_vector_to_fundamental_sector`.
        """

        candidates = self.equivalent_vectors(vector, antipodal=antipodal)
        return as_float_array(candidates[_canonical_vector_index(candidates)], shape=(3,))

    def canonicalize_vectors(
        self,
        vectors: ArrayLike | VectorSet,
        *,
        antipodal: bool = False,
    ) -> np.ndarray | VectorSet:
        """Canonical representatives for many directions at once.

        The batch form of :meth:`canonicalize_vector`. A ``VectorSet`` in gives
        a ``VectorSet`` out with its frame and provenance preserved.
        """

        reference_frame = None
        provenance = None
        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            normalized = normalize_vectors(vectors.values)
            reference_frame = vectors.reference_frame
            provenance = vectors.provenance
        else:
            normalized = normalize_vectors(vectors)
        canonicalized = [
            self.canonicalize_vector(vector, antipodal=antipodal) for vector in normalized
        ]
        array = np.stack(canonicalized, axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        if reference_frame is not None:
            return VectorSet(
                values=array,
                reference_frame=reference_frame,
                provenance=provenance,
            )
        return array

    def fundamental_sector(self, *, antipodal: bool = True) -> FundamentalSector:
        """The fundamental sector of this symmetry on the unit sphere.

        Purpose
        -------
        The spherical region containing exactly one representative of every
        direction orbit — the standard stereographic triangle for cubic
        symmetry, and the domain an inverse pole figure is plotted in.

        Parameters
        ----------
        antipodal : bool
            Treat a direction and its reverse as equivalent (default), which
            halves the sector and matches how IPFs are conventionally drawn.

        Returns
        -------
        FundamentalSector
            Carrying the sector vertices and bounding-plane normals.
        """

        return FundamentalSector(
            point_group=self.point_group,
            proper_point_group=self.proper_point_group,
            antipodal=antipodal,
            vertices=_sector_vertices_for_group(self.proper_point_group),
            edge_normals=_sector_edge_normals_for_group(self.proper_point_group),
        )

    def vector_in_fundamental_sector(self, vector: ArrayLike, *, antipodal: bool = True) -> bool:
        """Whether a direction already lies in the fundamental sector.

        With ``antipodal=True`` (default) the direction is first folded onto the
        upper hemisphere, matching the convention of
        :meth:`fundamental_sector`.
        """

        candidate = normalize_vector(vector)
        if antipodal and candidate[2] < 0.0:
            candidate = -candidate
        return _vector_in_fundamental_sector(candidate, self.proper_point_group)

    def reduce_vector_to_fundamental_sector(
        self,
        vector: ArrayLike,
        *,
        antipodal: bool = True,
    ) -> np.ndarray:
        """The representative of a direction inside the fundamental sector.

        Purpose
        -------
        The reduction behind inverse pole figures and IPF colouring: fold a
        direction into the standard triangle so that symmetry-equivalent
        directions land on the same point.

        Parameters
        ----------
        vector : ArrayLike
            A single direction; normalized internally.
        antipodal : bool
            Treat a direction and its reverse as equivalent (default).

        Returns
        -------
        np.ndarray
            The unit representative. When no orbit member satisfies the sector
            test exactly — possible for a direction lying on a sector boundary
            within numerical tolerance — the canonical representative is
            returned instead, so the function always yields a usable direction.
        """

        candidates = self.equivalent_vectors(vector, antipodal=antipodal)
        matching = [
            candidate
            for candidate in candidates
            if _vector_in_fundamental_sector(candidate, self.proper_point_group)
        ]
        if matching:
            selected = max(matching, key=_sector_sort_key)
            return as_float_array(selected, shape=(3,))
        return self.canonicalize_vector(vector, antipodal=antipodal)

    def reduce_vectors_to_fundamental_sector(
        self,
        vectors: ArrayLike | VectorSet,
        *,
        antipodal: bool = True,
    ) -> np.ndarray | VectorSet:
        """Fundamental-sector representatives for many directions at once.

        The batch form of :meth:`reduce_vector_to_fundamental_sector`, and the
        call behind every inverse pole figure computed from a full EBSD map. A
        ``VectorSet`` in gives a ``VectorSet`` out with frame and provenance
        preserved.
        """

        reference_frame = None
        provenance = None
        if isinstance(vectors, VectorSet):
            if self.reference_frame is not None and vectors.reference_frame != self.reference_frame:
                raise ValueError(
                    "VectorSet.reference_frame must match SymmetrySpec.reference_frame."
                )
            normalized = normalize_vectors(vectors.values)
            reference_frame = vectors.reference_frame
            provenance = vectors.provenance
        else:
            normalized = normalize_vectors(vectors)
        reduced = [
            self.reduce_vector_to_fundamental_sector(vector, antipodal=antipodal)
            for vector in normalized
        ]
        array = np.stack(reduced, axis=0)
        array = np.ascontiguousarray(array)
        array.setflags(write=False)
        if reference_frame is not None:
            return VectorSet(
                values=array,
                reference_frame=reference_frame,
                provenance=provenance,
            )
        return array


@dataclass(frozen=True, slots=True)
class FundamentalSector:
    """The spherical region holding one representative of each direction orbit.

    Purpose
    -------
    The standard stereographic triangle, generalized to any point group. It
    is the domain an inverse pole figure is drawn in and the region IPF
    colouring maps over.

    Attributes
    ----------
    point_group, proper_point_group : str
        The group and its rotation subgroup.
    antipodal : bool
        Whether ``+v`` and ``-v`` are treated as equivalent, which halves the
        sector.
    vertices : np.ndarray
        Corner directions of the sector.
    edge_normals : np.ndarray
        Inward normals of the bounding planes, used by :meth:`contains`. An
        empty set means triclinic symmetry, where the whole sphere is the
        sector.
    """

    point_group: str
    proper_point_group: str
    antipodal: bool
    vertices: np.ndarray
    edge_normals: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    def __post_init__(self) -> None:
        vertices = as_float_array(self.vertices, shape=(None, 3))
        if vertices.shape[0] > 0:
            vertices = normalize_vectors(vertices)
        edge_normals = as_float_array(self.edge_normals, shape=(None, 3))
        if edge_normals.shape[0] > 0:
            edge_normals = normalize_vectors(edge_normals)
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "edge_normals", edge_normals)

    def contains(self, vectors: ArrayLike, *, atol: float = _SECTOR_TOLERANCE) -> np.ndarray:
        """Whether each vector lies inside the sector.

        Tests every vector against the sector's bounding-plane normals with a
        tolerance, so directions on a boundary are accepted rather than
        rejected by floating-point noise. Accepts a single ``(3,)`` vector or an
        ``(n, 3)`` array and returns a scalar or ``(n,)`` boolean accordingly. A
        sector with no bounding planes — triclinic symmetry — contains
        everything.
        """

        array = np.asarray(vectors, dtype=np.float64)
        squeeze = array.ndim == 1
        rows = np.atleast_2d(array)
        if rows.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        if self.edge_normals.shape[0] == 0:
            mask = np.ones(rows.shape[0], dtype=bool)
        else:
            mask = np.all(rows @ self.edge_normals.T >= -atol, axis=1)
        if squeeze:
            return np.asarray(mask[0])
        mask = np.ascontiguousarray(mask)
        mask.setflags(write=False)
        return mask

    def center(self) -> np.ndarray:
        """A representative interior direction of the sector.

        The normalized vertex centroid for a genuine polygonal sector, falling
        back to the centroid of the boundary trace for degenerate sectors and to
        ``[0, 0, 1]`` when there are no vertices at all. Used to place labels
        and to seed searches inside the sector.
        """

        if self.vertices.shape[0] == 0:
            return as_float_array([0.0, 0.0, 1.0], shape=(3,))
        if self.vertices.shape[0] >= 3:
            return normalize_vector(self.vertices.sum(axis=0))
        trace = self.boundary_trace()
        return normalize_vector(trace.sum(axis=0))

    def boundary_trace(self, *, samples_per_edge: int = 64) -> np.ndarray:
        """A polyline tracing the sector boundary on the unit sphere.

        Purpose
        -------
        The outline drawn around a standard stereographic triangle. Sampling
        along great-circle edges rather than joining vertices with straight
        lines is what makes the drawn boundary follow the actual spherical
        geometry.

        Parameters
        ----------
        samples_per_edge : int
            Points per edge; at least two. Higher values give a smoother curve.

        Returns
        -------
        np.ndarray
            ``(m, 3)`` unit vectors along the boundary. For triclinic symmetry,
            where there is no bounded sector, the equator is returned.
        """

        if samples_per_edge < 2:
            raise ValueError("boundary_trace requires at least two samples per edge.")
        if self.vertices.shape[0] == 0:
            angles = np.linspace(0.0, 2.0 * np.pi, 4 * samples_per_edge)
            trace = np.column_stack(
                [np.cos(angles), np.sin(angles), np.zeros_like(angles)]
            )
            return as_float_array(np.ascontiguousarray(trace), shape=(None, 3))
        if self.edge_normals.shape[0] != self.vertices.shape[0]:
            raise ValueError(
                "boundary_trace requires one edge normal per vertex of the sector loop."
            )
        segments: list[np.ndarray] = []
        vertex_count = self.vertices.shape[0]
        for index in range(vertex_count):
            start = self.vertices[index]
            end = self.vertices[(index + 1) % vertex_count]
            normal = self.edge_normals[index]
            tangent = np.cross(normal, start)
            tangent_norm = float(np.linalg.norm(tangent))
            if np.isclose(tangent_norm, 0.0):
                raise ValueError(
                    "Sector edge normal must be perpendicular to its starting vertex."
                )
            tangent = tangent / tangent_norm
            sweep = float(
                np.mod(np.arctan2(float(tangent @ end), float(start @ end)), 2.0 * np.pi)
            )
            angles = np.linspace(0.0, sweep, samples_per_edge, endpoint=False)
            segments.append(
                np.outer(np.cos(angles), start) + np.outer(np.sin(angles), tangent)
            )
        segments.append(self.vertices[:1])
        trace = np.concatenate(segments, axis=0)
        return as_float_array(np.ascontiguousarray(trace), shape=(None, 3))
