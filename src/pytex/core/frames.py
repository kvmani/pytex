"""Reference frames and frame-to-frame transforms: the shared geometric foundation.

Every scientific quantity in PyTex — a direction, a plane normal, a pole, a
detector coordinate, an orientation — is only meaningful relative to a frame.
This module fixes *one* model for that, so no subsystem has to invent its own
(`AGENTS.md`, "No subsystem may define its own private frame model").

Three types make up the model:

- `ReferenceFrame` — the identity **and** the geometry of a frame: a name, a
  `FrameDomain` from the fixed repository vocabulary, three axis labels, the
  Cartesian components of those three axes, handedness, convention set, and
  provenance.
- `FrameTransform` — a typed, invertible, composable rigid map between exactly
  two named frames. Source and target are always explicit; a transform can
  never be applied to data living in the wrong frame.
- `FrameGraph` — a registry of frames and the transforms declared between them,
  which resolves the transform between *any* two connected frames by shortest
  path, so a workflow declares only the relationships it actually measured.

Axis-vector convention
----------------------

`ReferenceFrame.axis_vectors` holds the components of the frame's three
labelled axes **in the canonical right-handed Cartesian reference** ``X, Y, Z``
(`pytex.core.frame_catalog.CARTESIAN_FRAME`). The default is the identity
triad, meaning "this frame's axes coincide with the canonical Cartesian axes",
which is the standing convention for the specimen frame (``x, y, z``), the
rolling frame (``RD, TD, ND``), and the default crystal frame (``a, b, c``).

That single convention is what makes `FrameTransform.between_frames` well
defined: two frames whose axis vectors are given in a common reference have a
computable relative rotation without any further declaration.

Axis vectors are *dimensionless orientations*. Physical axis lengths belong to
`pytex.core.lattice.Basis`, which carries a `BasisKind` and a unit; the two are
complementary and must not be conflated.

Transform direction
-------------------

`FrameTransform.rotation_matrix` maps **components in the source frame to
components in the target frame**:

``v_target = R @ v_source + t``

so a transform from the crystal frame to the specimen frame takes crystal-frame
components and returns specimen-frame components.

See also
--------
`pytex.core.frame_catalog` : named standard frames built on this model.
`pytex.plotting.frames` : triads, embeddable gizmos, and documentation SVG.
`docs/standards/notation_and_conventions.md` : the canonical frame chain.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.batches import VectorSet
from pytex.core.conventions import (
    PYTEX_CANONICAL_CONVENTIONS,
    ConventionSet,
    FrameDomain,
    Handedness,
)
from pytex.core.provenance import ProvenanceRecord

if TYPE_CHECKING:  # pragma: no cover - import cycle avoidance only
    from pytex.core.orientation import Rotation

__all__ = [
    "IDENTITY_AXIS_VECTORS",
    "AxisVectors",
    "FrameGraph",
    "FrameTransform",
    "ReferenceFrame",
    "as_axis_vectors",
]

AxisVectors = tuple[tuple[float, float, float], ...]

#: The canonical Cartesian triad, used as the default frame geometry.
IDENTITY_AXIS_VECTORS: AxisVectors = (
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
)

# Degenerate-basis guard. A triad whose determinant falls below this magnitude
# cannot define a frame: components in it would not be uniquely recoverable.
_MIN_ABS_DETERMINANT = 1e-9
_ORTHONORMAL_ATOL = 1e-8


def as_axis_vectors(values: ArrayLike) -> AxisVectors:
    """Normalize any axis-vector input into the hashable ``3 x 3`` float tuple form.

    `ReferenceFrame` stores axis geometry as a tuple of tuples rather than an
    ndarray so that frames stay comparable — frame equality is load-bearing
    throughout PyTex. Use this helper wherever an axis triad arrives as an
    array, a list of lists, or a NumPy matrix and has to be handed to a frame.

    Parameters
    ----------
    values:
        Anything array-like with shape ``(3, 3)``; row ``i`` is the ``i``-th
        axis vector in canonical Cartesian components.

    Returns
    -------
    tuple[tuple[float, float, float], ...]
        The same triad as plain Python floats.

    Raises
    ------
    ValueError
        If the input is not ``(3, 3)`` or contains non-finite entries.
    """

    array = np.asarray(values, dtype=np.float64)
    if array.shape != (3, 3):
        raise ValueError(
            "ReferenceFrame.axis_vectors must contain exactly three 3-component axis vectors."
        )
    if not np.all(np.isfinite(array)):
        raise ValueError("ReferenceFrame.axis_vectors must be finite.")
    return tuple(
        (float(row[0]), float(row[1]), float(row[2])) for row in array
    )


@dataclass(frozen=True, slots=True)
class ReferenceFrame:
    """A named, domain-typed coordinate frame with explicit axis geometry.

    What it does
        Binds together everything needed to interpret a triple of numbers as a
        physical vector: which frame the numbers belong to (`name`, `domain`),
        what the three components are called (`axes`), where those axes actually
        point (`axis_vectors`, in canonical Cartesian components), the chirality
        (`handedness`), the convention set in force, and where the frame came
        from (`provenance`).

    When to use it
        Whenever a stable surface accepts or returns 3-component data. Public
        APIs must not take naked arrays where the frame would be ambiguous. Reach
        for `pytex.core.frame_catalog` first — it already provides the standard
        frames (Cartesian ``X/Y/Z``, sample ``RD/TD/ND``, crystal ``a/b/c``, EBSD
        map, detector, laboratory) — and construct a `ReferenceFrame` directly
        only for a frame the catalog does not cover.

    Parameters
    ----------
    name:
        Stable identifier, unique within a workflow or `FrameGraph`.
    domain:
        A member of the fixed `FrameDomain` vocabulary. New domains may not be
        invented (`docs/standards/notation_and_conventions.md`).
    axes:
        The three axis labels, in order, e.g. ``("RD", "TD", "ND")``.
    handedness:
        `Handedness.RIGHT` (the canonical default) or `Handedness.LEFT`. Must
        agree with the sign of the axis-vector determinant.
    convention:
        The `ConventionSet` under which the frame is interpreted.
    description:
        Free prose describing the frame's physical attachment.
    provenance:
        Optional import/source record.
    metadata:
        Optional string metadata; copied into a read-only mapping.
    axis_vectors:
        Components of the three labelled axes in the canonical Cartesian
        reference, as a hashable ``3 x 3`` tuple of floats. Defaults to the
        identity triad. Need not be orthonormal (an oblique crystal frame is
        legitimate) but must be linearly independent.
    axis_descriptions:
        Optional long names for the three axes, e.g. ``("rolling direction",
        "transverse direction", "normal direction")``. Either empty or length 3.

    Raises
    ------
    ValueError
        If the axis labels are not exactly three, the axis vectors are
        degenerate, the determinant sign contradicts `handedness`, or
        `axis_descriptions` has the wrong length.

    See also
    --------
    `FrameTransform` : the typed map between two frames.
    `pytex.core.frame_catalog` : ready-made standard frames.
    """

    name: str
    domain: FrameDomain
    axes: tuple[str, str, str]
    handedness: Handedness = Handedness.RIGHT
    convention: ConventionSet = PYTEX_CANONICAL_CONVENTIONS
    description: str = ""
    provenance: ProvenanceRecord | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    axis_vectors: AxisVectors = IDENTITY_AXIS_VECTORS
    axis_descriptions: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if len(self.axes) != 3:
            raise ValueError("ReferenceFrame.axes must contain exactly three axis labels.")
        object.__setattr__(self, "axes", tuple(str(label) for label in self.axes))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

        vectors = as_axis_vectors(self.axis_vectors)
        object.__setattr__(self, "axis_vectors", vectors)
        determinant = float(np.linalg.det(np.asarray(vectors, dtype=np.float64).T))
        if abs(determinant) < _MIN_ABS_DETERMINANT:
            raise ValueError(
                f"ReferenceFrame '{self.name}' has linearly dependent axis vectors "
                f"(determinant {determinant:.3e}); a frame must span three dimensions."
            )
        expected_sign = 1.0 if self.handedness is Handedness.RIGHT else -1.0
        if np.sign(determinant) != expected_sign:
            raise ValueError(
                f"ReferenceFrame '{self.name}' declares {self.handedness.value}-handedness but "
                f"its axis vectors have determinant {determinant:.6f}; the sign must match."
            )

        descriptions = tuple(str(item) for item in self.axis_descriptions)
        if descriptions and len(descriptions) != 3:
            raise ValueError(
                "ReferenceFrame.axis_descriptions must be empty or contain exactly three entries."
            )
        object.__setattr__(self, "axis_descriptions", descriptions)

    # ------------------------------------------------------------------ #
    # Geometry
    # ------------------------------------------------------------------ #

    @property
    def basis_matrix(self) -> np.ndarray:
        """The ``(3, 3)`` matrix whose **columns** are the frame's axis vectors.

        Multiplying this matrix by a component triple expressed in this frame
        returns the same vector's canonical Cartesian components:
        ``x_cartesian = frame.basis_matrix @ v_frame``.
        """

        return as_float_array(np.asarray(self.axis_vectors, dtype=np.float64).T, shape=(3, 3))

    @property
    def determinant(self) -> float:
        """Determinant of `basis_matrix`; its sign encodes the handedness."""

        return float(np.linalg.det(self.basis_matrix))

    @property
    def is_right_handed(self) -> bool:
        """Whether the axis triad is right-handed (positive determinant)."""

        return self.determinant > 0.0

    @property
    def is_orthonormal(self) -> bool:
        """Whether the three axis vectors are mutually orthogonal unit vectors.

        Crystal frames for non-cubic phases are legitimately non-orthonormal, so
        this is reported rather than enforced. Visualization normalizes
        non-orthonormal triads for legibility.
        """

        matrix = self.basis_matrix
        return bool(np.allclose(matrix.T @ matrix, np.eye(3), atol=_ORTHONORMAL_ATOL))

    def axis_index(self, axis: str | int) -> int:
        """Resolve an axis label (case-insensitive) or integer index to ``0..2``.

        Parameters
        ----------
        axis:
            An axis label such as ``"RD"`` or ``"nd"``, or an integer index.

        Returns
        -------
        int
            The position of the axis in `axes`.

        Raises
        ------
        KeyError
            If the label does not name an axis of this frame.
        IndexError
            If an integer index is outside ``0..2``.
        """

        if isinstance(axis, int | np.integer):
            index = int(axis)
            if not 0 <= index < 3:
                raise IndexError(f"Axis index {index} is out of range for a three-axis frame.")
            return index
        wanted = str(axis).strip().lower()
        for index, label in enumerate(self.axes):
            if label.lower() == wanted:
                return index
        raise KeyError(
            f"Frame '{self.name}' has no axis '{axis}'; available axes: {', '.join(self.axes)}."
        )

    def axis_vector(self, axis: str | int, *, normalize: bool = True) -> np.ndarray:
        """Return one axis of the frame as canonical Cartesian components.

        Parameters
        ----------
        axis:
            Axis label or index, resolved by `axis_index`.
        normalize:
            When ``True`` (default) the returned vector has unit length, which
            is what direction-facing consumers (arrows, poles, projections)
            want. Set ``False`` to get the stored vector verbatim.

        Returns
        -------
        numpy.ndarray
            A read-only ``(3,)`` array.
        """

        vector = np.asarray(self.axis_vectors[self.axis_index(axis)], dtype=np.float64)
        if normalize:
            norm = float(np.linalg.norm(vector))
            if norm == 0.0:  # pragma: no cover - blocked by the degeneracy check
                raise ValueError(f"Axis '{axis}' of frame '{self.name}' is a zero vector.")
            vector = vector / norm
        return as_float_array(vector, shape=(3,))

    def unit_axis_matrix(self) -> np.ndarray:
        """`basis_matrix` with every column scaled to unit length."""

        matrix = np.asarray(self.basis_matrix, dtype=np.float64)
        norms = np.linalg.norm(matrix, axis=0, keepdims=True)
        return as_float_array(matrix / norms, shape=(3, 3))

    def axis_description(self, axis: str | int) -> str:
        """Long name of one axis, falling back to its label when none is set."""

        index = self.axis_index(axis)
        if self.axis_descriptions:
            return self.axis_descriptions[index]
        return self.axes[index]

    # ------------------------------------------------------------------ #
    # Derivation
    # ------------------------------------------------------------------ #

    def with_axis_vectors(self, axis_vectors: ArrayLike) -> ReferenceFrame:
        """Return a copy of this frame carrying different axis geometry.

        Use it to record that a nominally standard frame is physically rotated —
        for example a specimen mounted with its rolling direction 30 degrees off
        the stage ``X`` axis.
        """

        return self._replace(axis_vectors=as_axis_vectors(axis_vectors))

    def renamed(self, name: str, *, description: str | None = None) -> ReferenceFrame:
        """Return a copy of this frame under a new name (and optional description)."""

        return self._replace(
            name=str(name),
            description=self.description if description is None else str(description),
        )

    def rotated(self, rotation: Rotation, *, name: str | None = None) -> ReferenceFrame:
        """Return a copy of this frame with its axes rotated by ``rotation``.

        The rotation acts on the axis vectors in canonical Cartesian components
        (``a_new = R @ a_old``), so the returned frame keeps the same labels and
        domain while pointing somewhere else in space.

        Parameters
        ----------
        rotation:
            A `pytex.core.orientation.Rotation`.
        name:
            Optional new name; defaults to ``"<name>_rotated"``.
        """

        matrix = np.asarray(rotation.as_matrix(), dtype=np.float64)
        rotated_vectors = np.asarray(self.axis_vectors, dtype=np.float64) @ matrix.T
        return self._replace(
            name=f"{self.name}_rotated" if name is None else str(name),
            axis_vectors=as_axis_vectors(rotated_vectors),
        )

    def _replace(self, **changes: Any) -> ReferenceFrame:
        """Internal copy-with-changes that preserves every unspecified field."""

        fields: dict[str, Any] = {
            "name": self.name,
            "domain": self.domain,
            "axes": self.axes,
            "handedness": self.handedness,
            "convention": self.convention,
            "description": self.description,
            "provenance": self.provenance,
            "metadata": dict(self.metadata),
            "axis_vectors": self.axis_vectors,
            "axis_descriptions": self.axis_descriptions,
        }
        fields.update(changes)
        return ReferenceFrame(**fields)

    # ------------------------------------------------------------------ #
    # Explanation
    # ------------------------------------------------------------------ #

    def describe(self) -> str:
        """Convention-explicit prose summary of the frame.

        States the domain, the axis labels and where they point, the handedness,
        the orthonormality status, and the governing convention set — everything
        a reader needs to interpret components quoted in this frame.
        """

        axis_parts: list[str] = []
        for index, label in enumerate(self.axes):
            vector = self.axis_vector(index, normalize=False)
            long_name = self.axis_description(index)
            suffix = f" ({long_name})" if long_name != label else ""
            axis_parts.append(
                f"{label}{suffix} -> [{vector[0]:.4f} {vector[1]:.4f} {vector[2]:.4f}]"
            )
        geometry = (
            "orthonormal" if self.is_orthonormal else "non-orthonormal (oblique or scaled)"
        )
        description = f" {self.description.strip()}" if self.description.strip() else ""
        provenance = (
            f" Imported from {self.provenance.source_system}." if self.provenance else ""
        )
        return (
            f"Reference frame '{self.name}' in the {self.domain.value} domain, "
            f"{self.handedness.value}-handed and {geometry}, with axes "
            f"{', '.join(axis_parts)} given as components in the canonical Cartesian "
            f"reference (X, Y, Z). Convention set: {self.convention.name} "
            f"({self.convention.angle_convention.value} Euler labeling, quaternions stored "
            f"{''.join(self.convention.quaternion_order)}).{description}{provenance}"
        )


@dataclass(frozen=True, slots=True)
class FrameTransform:
    """A typed rigid map from one reference frame's components to another's.

    What it does
        Stores an orthonormal ``rotation_matrix`` with determinant ``+1`` and a
        ``translation_vector``, together with the `ReferenceFrame` on each end,
        and applies them as ``v_target = R @ v_source + t``. Because both frames
        are named, applying a transform to data from the wrong frame is a
        construction-time error rather than a silent numerical mistake.

    When to use it
        Whenever data crosses a frame boundary: crystal to specimen, specimen to
        map, specimen to detector. For chains that span more than one declared
        relationship, register the transforms in a `FrameGraph` and let it
        compose them.

    Parameters
    ----------
    source, target:
        The frames the transform maps between.
    rotation_matrix:
        ``(3, 3)`` proper rotation, validated at construction.
    translation_vector:
        ``(3,)`` offset applied after the rotation. Zero for pure
        orientation relationships; non-zero when frame origins differ (e.g. a
        map origin offset from the specimen origin).
    provenance:
        Optional record of where the transform was measured or declared.

    See also
    --------
    `FrameGraph.transform_between` : automatic multi-step resolution.
    `ReferenceFrame` : the endpoints.
    """

    source: ReferenceFrame
    target: ReferenceFrame
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray = field(default_factory=lambda: np.zeros(3))
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        rotation = as_float_array(self.rotation_matrix, shape=(3, 3))
        translation = as_float_array(self.translation_vector, shape=(3,))
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8):
            raise ValueError("rotation_matrix must be orthonormal.")
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8):
            raise ValueError("rotation_matrix must have determinant +1.")
        object.__setattr__(self, "rotation_matrix", rotation)
        object.__setattr__(self, "translation_vector", translation)

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def identity(cls, frame: ReferenceFrame) -> FrameTransform:
        """The do-nothing transform from a frame to itself."""

        return cls(source=frame, target=frame, rotation_matrix=np.eye(3))

    @classmethod
    def from_rotation(
        cls,
        rotation: Rotation,
        *,
        source: ReferenceFrame,
        target: ReferenceFrame,
        translation: ArrayLike = (0.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> FrameTransform:
        """Build a transform from a `pytex.core.orientation.Rotation`.

        The rotation matrix is taken to act on source-frame components, i.e.
        ``v_target = R @ v_source``. This is the natural entry point when a
        relationship is already known as a rotation, a quaternion, or an
        orientation matrix.
        """

        return cls(
            source=source,
            target=target,
            rotation_matrix=rotation.as_matrix(),
            translation_vector=np.asarray(translation, dtype=np.float64),
            provenance=provenance,
        )

    @classmethod
    def from_bunge_euler(
        cls,
        phi1: float,
        Phi: float,  # noqa: N803 - Bunge's canonical symbol is capital Phi
        phi2: float,
        *,
        source: ReferenceFrame,
        target: ReferenceFrame,
        degrees: bool = True,
        translation: ArrayLike = (0.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> FrameTransform:
        """Build a transform from Bunge ``(phi1, Phi, phi2)`` Euler angles.

        Uses the repository-canonical Bunge ZXZ labeling
        (`docs/standards/notation_and_conventions.md`). Angles are in degrees
        unless ``degrees=False``.
        """

        from pytex.core.orientation import Rotation

        rotation = Rotation.from_bunge_euler(phi1, Phi, phi2, degrees=degrees)
        return cls.from_rotation(
            rotation,
            source=source,
            target=target,
            translation=translation,
            provenance=provenance,
        )

    @classmethod
    def from_axis_angle(
        cls,
        axis: ArrayLike,
        angle: float,
        *,
        source: ReferenceFrame,
        target: ReferenceFrame,
        degrees: bool = True,
        translation: ArrayLike = (0.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> FrameTransform:
        """Build a transform from a rotation axis and angle.

        Typical use: a specimen re-mounted with a known tilt, or a detector
        rotated by a known angle about the beam.
        """

        from pytex.core.orientation import Rotation

        angle_rad = float(np.deg2rad(angle)) if degrees else float(angle)
        rotation = Rotation.from_axis_angle(axis, angle_rad)
        return cls.from_rotation(
            rotation,
            source=source,
            target=target,
            translation=translation,
            provenance=provenance,
        )

    @classmethod
    def from_axis_correspondence(
        cls,
        source: ReferenceFrame,
        target: ReferenceFrame,
        correspondence: Mapping[str, str],
        *,
        translation: ArrayLike = (0.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> FrameTransform:
        """Build a transform by declaring which target axis each source axis lands on.

        This is the most readable way to express the axis-relabeling conventions
        that differ between EBSD vendors and between tools: rather than writing a
        permutation matrix by hand, state the correspondence in words.

        Parameters
        ----------
        source, target:
            The frames to relate.
        correspondence:
            Maps each source axis label to a target axis label, optionally with a
            leading ``"-"`` for a reversed sense, e.g.
            ``{"x": "RD", "y": "TD", "z": "ND"}`` or
            ``{"x": "TD", "y": "-RD", "z": "ND"}``. All three source axes must
            appear, and the three target axes must be distinct.
        translation:
            Optional origin offset applied after the rotation.

        Returns
        -------
        FrameTransform
            The transform whose rotation carries each source axis direction onto
            its declared target axis direction.

        Raises
        ------
        ValueError
            If the correspondence is incomplete, repeats a target axis, or does
            not describe a proper rotation (for instance an odd permutation with
            no compensating sign flip, which would mirror rather than rotate).

        Notes
        -----
        The declaration is about *component* semantics and is therefore
        independent of where either frame's axes happen to point in the canonical
        Cartesian reference: in its own coordinates a frame's axis ``i`` is the
        standard basis vector ``e_i``, so "source axis ``i`` is target axis ``j``"
        fixes ``R e_i = +/- e_j``, giving a signed permutation matrix.
        """

        if len(correspondence) != 3:
            raise ValueError(
                "from_axis_correspondence requires exactly three source-axis entries; "
                f"received {len(correspondence)}."
            )

        rotation = np.zeros((3, 3), dtype=np.float64)
        seen_sources: set[int] = set()
        seen_targets: set[int] = set()
        for source_axis, target_axis in correspondence.items():
            source_index = source.axis_index(source_axis)
            if source_index in seen_sources:
                raise ValueError(
                    f"Source axis '{source_axis}' appears more than once in the "
                    "correspondence; each source axis may be used exactly once."
                )
            seen_sources.add(source_index)

            label = str(target_axis).strip()
            sign = 1.0
            if label.startswith("-"):
                sign = -1.0
                label = label[1:].strip()
            elif label.startswith("+"):
                label = label[1:].strip()
            target_index = target.axis_index(label)
            if target_index in seen_targets:
                raise ValueError(
                    f"Target axis '{label}' appears more than once in the correspondence; "
                    "each target axis may be used exactly once."
                )
            seen_targets.add(target_index)
            # R e_source = sign * e_target.
            rotation[target_index, source_index] = sign

        if not np.isclose(float(np.linalg.det(rotation)), 1.0, atol=1e-8):
            raise ValueError(
                "The declared axis correspondence does not describe a proper rotation "
                "(determinant is not +1). An odd axis permutation needs a sign flip on one "
                "axis to stay right-handed."
            )
        return cls(
            source=source,
            target=target,
            rotation_matrix=rotation,
            translation_vector=np.asarray(translation, dtype=np.float64),
            provenance=provenance,
        )

    @classmethod
    def between_frames(
        cls,
        source: ReferenceFrame,
        target: ReferenceFrame,
        *,
        translation: ArrayLike = (0.0, 0.0, 0.0),
        provenance: ProvenanceRecord | None = None,
    ) -> FrameTransform:
        """Derive the transform implied by the two frames' own axis geometry.

        Both frames' `ReferenceFrame.axis_vectors` are components in the same
        canonical Cartesian reference, so their relative rotation is already
        determined: ``R = B_target^-1 @ B_source``. Use this when both frames
        were defined against the canonical reference (which every catalog frame
        is) and no separately measured relationship exists.

        Raises
        ------
        ValueError
            If either frame is non-orthonormal, in which case the relationship is
            not a rigid rotation and must be declared explicitly instead.
        """

        for frame in (source, target):
            if not frame.is_orthonormal:
                raise ValueError(
                    f"Frame '{frame.name}' is non-orthonormal, so the frame-to-frame map is not "
                    "a rigid rotation. Declare the relationship explicitly with "
                    "FrameTransform.from_rotation instead."
                )
        rotation = np.linalg.inv(target.basis_matrix) @ source.basis_matrix
        return cls(
            source=source,
            target=target,
            rotation_matrix=rotation,
            translation_vector=np.asarray(translation, dtype=np.float64),
            provenance=provenance,
        )

    # ------------------------------------------------------------------ #
    # Properties
    # ------------------------------------------------------------------ #

    @property
    def is_identity(self) -> bool:
        """Whether the transform leaves components unchanged."""

        return bool(
            np.allclose(self.rotation_matrix, np.eye(3), atol=1e-10)
            and np.allclose(self.translation_vector, 0.0, atol=1e-10)
        )

    @property
    def rotation_angle_deg(self) -> float:
        """The rotation angle in degrees, in ``[0, 180]``."""

        return float(self.as_rotation().angle_deg)

    @property
    def rotation_axis(self) -> np.ndarray:
        """The unit rotation axis as canonical ``(3,)`` components."""

        return np.asarray(self.as_rotation().axis, dtype=np.float64)

    def as_rotation(self) -> Rotation:
        """Return the rotational part as a `pytex.core.orientation.Rotation`."""

        from pytex.core.orientation import Rotation

        return Rotation.from_matrix(self.rotation_matrix)

    # ------------------------------------------------------------------ #
    # Application
    # ------------------------------------------------------------------ #

    def apply_to_vectors(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Map position-like vectors, applying rotation **and** translation.

        Accepts a `VectorSet` (whose frame must equal `source`; the result
        carries `target`) or any array whose last axis has length 3.
        """

        if isinstance(vectors, VectorSet):
            if vectors.reference_frame != self.source:
                raise ValueError("VectorSet.reference_frame must match FrameTransform.source.")
            transformed = vectors.values @ self.rotation_matrix.T + self.translation_vector
            return VectorSet(
                values=transformed,
                reference_frame=self.target,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = array @ self.rotation_matrix.T + self.translation_vector
        transformed = np.ascontiguousarray(transformed)
        transformed.setflags(write=False)
        return transformed

    def apply_to_directions(self, directions: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Map direction-like vectors, applying the rotation **only**.

        Directions, plane normals, and poles are translation-invariant: an origin
        offset must not move them. Use this instead of `apply_to_vectors`
        whenever the quantity is a direction rather than a position.
        """

        if isinstance(directions, VectorSet):
            if directions.reference_frame != self.source:
                raise ValueError("VectorSet.reference_frame must match FrameTransform.source.")
            return VectorSet(
                values=directions.values @ self.rotation_matrix.T,
                reference_frame=self.target,
                provenance=directions.provenance,
            )
        array = np.asarray(directions, dtype=np.float64)
        if array.shape[-1] != 3:
            raise ValueError("Input vectors must end with dimension 3.")
        transformed = np.ascontiguousarray(array @ self.rotation_matrix.T)
        transformed.setflags(write=False)
        return transformed

    def source_axes_in_target(self) -> np.ndarray:
        """The source frame's three axes, as components in the target frame.

        What it does
            Returns the ``(3, 3)`` matrix whose **column** ``i`` holds the target-frame
            components of the source frame's ``i``-th axis. In its own coordinates a
            frame's axis ``i`` is the standard basis vector ``e_i``, so this is simply
            the rotation matrix — but naming it makes the intent explicit at the call
            site and keeps the geometric meaning documented.

        When to use it
            To draw both frames of a relationship in one picture: the target frame
            is the identity triad in its own coordinates, and this matrix places
            the source triad beside it. `pytex.plotting.frames.plot_frame_relationship`
            is built on exactly this.

        Returns
        -------
        numpy.ndarray
            A read-only ``(3, 3)`` array of target-frame components.

        Notes
        -----
        These components are expressed in the **target frame**, not in the
        canonical Cartesian reference that `ReferenceFrame.axis_vectors` uses.
        That is why this returns a matrix rather than a `ReferenceFrame`: wrapping
        it in a frame would quietly break the module's axis-vector convention.
        """

        return as_float_array(self.rotation_matrix, shape=(3, 3))

    # ------------------------------------------------------------------ #
    # Algebra
    # ------------------------------------------------------------------ #

    def inverse(self) -> FrameTransform:
        """The transform back from `target` to `source`."""

        inverse_rotation = self.rotation_matrix.T
        inverse_translation = -(self.translation_vector @ self.rotation_matrix)
        return FrameTransform(
            source=self.target,
            target=self.source,
            rotation_matrix=inverse_rotation,
            translation_vector=inverse_translation,
            provenance=self.provenance,
        )

    def compose(self, previous: FrameTransform) -> FrameTransform:
        """Chain ``previous`` then ``self``; requires ``previous.target == self.source``."""

        if previous.target != self.source:
            raise ValueError(
                "Transform domains do not chain: previous.target must equal self.source."
            )
        rotation = self.rotation_matrix @ previous.rotation_matrix
        translation = previous.translation_vector @ self.rotation_matrix.T + self.translation_vector
        return FrameTransform(
            source=previous.source,
            target=self.target,
            rotation_matrix=rotation,
            translation_vector=translation,
            provenance=self.provenance,
        )

    # ------------------------------------------------------------------ #
    # Explanation
    # ------------------------------------------------------------------ #

    def describe(self) -> str:
        """Convention-explicit prose summary of the transform.

        States the endpoints and their domains, the rotation angle and axis, the
        translation, and the direction in which components are mapped.
        """

        if self.is_identity:
            core = "the identity map (components are unchanged)"
        else:
            axis = self.rotation_axis
            core = (
                f"a rotation of {self.rotation_angle_deg:.4f} deg about the axis "
                f"[{axis[0]:.4f} {axis[1]:.4f} {axis[2]:.4f}]"
            )
        translation = self.translation_vector
        if np.allclose(translation, 0.0, atol=1e-12):
            offset = "no origin offset"
        else:
            offset = (
                f"an origin offset of [{translation[0]:.4f} {translation[1]:.4f} "
                f"{translation[2]:.4f}] in the target frame"
            )
        provenance = (
            f" Declared from {self.provenance.source_system}." if self.provenance else ""
        )
        return (
            f"Frame transform from '{self.source.name}' ({self.source.domain.value}) to "
            f"'{self.target.name}' ({self.target.domain.value}): {core}, with {offset}. "
            f"Applied as v_{self.target.name} = R v_{self.source.name} + t, so it converts "
            f"components expressed in '{self.source.name}' into components expressed in "
            f"'{self.target.name}'.{provenance}"
        )


class FrameGraph:
    """A registry of frames and the transforms declared between them.

    What it does
        Holds `ReferenceFrame` objects by name and `FrameTransform` objects as
        undirected edges (each usable in either direction via
        `FrameTransform.inverse`), then resolves the transform between any two
        connected frames by composing the shortest declared chain.

    When to use it
        When a workflow spans more than one frame relationship — an EBSD dataset
        that knows ``crystal -> specimen`` and ``specimen -> map``, or a
        diffraction setup that knows ``specimen -> detector`` and
        ``detector -> laboratory``. Declare only the relationships you actually
        measured; ask the graph for the rest.

    Why shortest path
        Each composition multiplies rotation matrices, so the fewest hops means
        the least accumulated floating-point error and the clearest provenance.

    Parameters
    ----------
    frames:
        Frames to register up front.
    transforms:
        Transforms to register up front; their endpoint frames are registered
        automatically.
    name:
        A label for the graph, used in `describe`.

    Raises
    ------
    ValueError
        If two different frames are registered under the same name.

    See also
    --------
    `pytex.core.frame_catalog.rolling_frame_graph` : a ready-made example graph.
    """

    __slots__ = ("_adjacency", "_frames", "_transforms", "name")

    def __init__(
        self,
        frames: Iterable[ReferenceFrame] = (),
        transforms: Iterable[FrameTransform] = (),
        *,
        name: str = "frame_graph",
    ) -> None:
        self.name = str(name)
        self._frames: dict[str, ReferenceFrame] = {}
        self._transforms: list[FrameTransform] = []
        self._adjacency: dict[str, list[tuple[str, int, bool]]] = {}
        for frame in frames:
            self.add_frame(frame)
        for transform in transforms:
            self.add_transform(transform)

    # ------------------------------------------------------------------ #
    # Registration
    # ------------------------------------------------------------------ #

    def add_frame(self, frame: ReferenceFrame) -> ReferenceFrame:
        """Register a frame; re-registering an identical frame is a no-op."""

        existing = self._frames.get(frame.name)
        if existing is not None and existing != frame:
            raise ValueError(
                f"Frame name '{frame.name}' is already registered with different definition; "
                "frame names must be unique within a graph."
            )
        self._frames.setdefault(frame.name, frame)
        self._adjacency.setdefault(frame.name, [])
        return frame

    def add_transform(self, transform: FrameTransform) -> FrameTransform:
        """Register a transform and both of its endpoint frames.

        The edge is usable in both directions; the reverse direction is realized
        with `FrameTransform.inverse` when a path traverses it backwards.
        """

        self.add_frame(transform.source)
        self.add_frame(transform.target)
        index = len(self._transforms)
        self._transforms.append(transform)
        self._adjacency[transform.source.name].append((transform.target.name, index, True))
        self._adjacency[transform.target.name].append((transform.source.name, index, False))
        return transform

    # ------------------------------------------------------------------ #
    # Access
    # ------------------------------------------------------------------ #

    def __len__(self) -> int:
        return len(self._frames)

    def __contains__(self, frame: str | ReferenceFrame) -> bool:
        return self.has_frame(frame)

    def has_frame(self, frame: str | ReferenceFrame) -> bool:
        """Whether a frame (by name or object) is registered."""

        name = frame if isinstance(frame, str) else frame.name
        return name in self._frames

    def frame(self, name: str) -> ReferenceFrame:
        """Look up a registered frame by name."""

        try:
            return self._frames[name]
        except KeyError:
            available = ", ".join(sorted(self._frames)) or "<none>"
            raise KeyError(
                f"Frame '{name}' is not registered in graph '{self.name}'. "
                f"Registered frames: {available}."
            ) from None

    def frames(self) -> tuple[ReferenceFrame, ...]:
        """Every registered frame, ordered by name."""

        return tuple(self._frames[key] for key in sorted(self._frames))

    def transforms(self) -> tuple[FrameTransform, ...]:
        """Every registered transform, in registration order."""

        return tuple(self._transforms)

    # ------------------------------------------------------------------ #
    # Resolution
    # ------------------------------------------------------------------ #

    def _resolve_name(self, frame: str | ReferenceFrame) -> str:
        name = frame if isinstance(frame, str) else frame.name
        if name not in self._frames:
            available = ", ".join(sorted(self._frames)) or "<none>"
            raise KeyError(
                f"Frame '{name}' is not registered in graph '{self.name}'. "
                f"Registered frames: {available}."
            )
        if not isinstance(frame, str) and self._frames[name] != frame:
            raise ValueError(
                f"Frame '{name}' is registered with a different definition than the one supplied."
            )
        return name

    def _edge_path(
        self, source: str, target: str
    ) -> tuple[tuple[int, bool], ...]:
        """Breadth-first shortest edge chain from ``source`` to ``target``."""

        if source == target:
            return ()
        previous: dict[str, tuple[str, int, bool]] = {}
        visited = {source}
        queue = [source]
        while queue:
            next_queue: list[str] = []
            for current in queue:
                for neighbour, index, forward in self._adjacency[current]:
                    if neighbour in visited:
                        continue
                    visited.add(neighbour)
                    previous[neighbour] = (current, index, forward)
                    if neighbour == target:
                        chain: list[tuple[int, bool]] = []
                        node = target
                        while node != source:
                            parent, edge_index, edge_forward = previous[node]
                            chain.append((edge_index, edge_forward))
                            node = parent
                        chain.reverse()
                        return tuple(chain)
                    next_queue.append(neighbour)
            queue = next_queue
        raise KeyError(
            f"No declared transform path connects '{source}' to '{target}' in graph "
            f"'{self.name}'. Add the missing relationship with add_transform."
        )

    def path(
        self, source: str | ReferenceFrame, target: str | ReferenceFrame
    ) -> tuple[str, ...]:
        """The shortest chain of frame names from ``source`` to ``target``.

        Returns
        -------
        tuple[str, ...]
            Frame names including both endpoints, e.g.
            ``("crystal", "specimen", "map")``. A single-element tuple is
            returned when source and target are the same frame.

        Raises
        ------
        KeyError
            If either frame is unregistered, or no declared path connects them.
        """

        source_name = self._resolve_name(source)
        target_name = self._resolve_name(target)
        names = [source_name]
        for index, forward in self._edge_path(source_name, target_name):
            transform = self._transforms[index]
            names.append(transform.target.name if forward else transform.source.name)
        return tuple(names)

    def transform_between(
        self, source: str | ReferenceFrame, target: str | ReferenceFrame
    ) -> FrameTransform:
        """Compose the shortest declared chain into one `FrameTransform`.

        Parameters
        ----------
        source, target:
            Registered frames, given as names or objects.

        Returns
        -------
        FrameTransform
            A transform whose `source` and `target` are the requested frames.
            Returns the identity transform when both are the same frame.

        Raises
        ------
        KeyError
            If a frame is unregistered or the frames are not connected.
        """

        source_name = self._resolve_name(source)
        target_name = self._resolve_name(target)
        if source_name == target_name:
            return FrameTransform.identity(self._frames[source_name])
        composed: FrameTransform | None = None
        for index, forward in self._edge_path(source_name, target_name):
            edge = self._transforms[index]
            step = edge if forward else edge.inverse()
            composed = step if composed is None else step.compose(composed)
        if composed is None:  # pragma: no cover - a non-empty edge path always composes
            raise KeyError(
                f"No declared transform path connects '{source_name}' to '{target_name}'."
            )
        return composed

    def convert(
        self,
        vectors: ArrayLike | VectorSet,
        *,
        source: str | ReferenceFrame,
        target: str | ReferenceFrame,
        directions: bool = False,
    ) -> np.ndarray | VectorSet:
        """Convert vector components from one registered frame to another.

        Parameters
        ----------
        vectors:
            A `VectorSet` (whose frame must be the source frame) or an array
            whose last axis has length 3.
        source, target:
            Registered frames, given as names or objects.
        directions:
            When ``True`` the origin offset is ignored, which is correct for
            directions, plane normals, and poles. Defaults to ``False``
            (position-like mapping).
        """

        transform = self.transform_between(source, target)
        if directions:
            return transform.apply_to_directions(vectors)
        return transform.apply_to_vectors(vectors)

    # ------------------------------------------------------------------ #
    # Explanation
    # ------------------------------------------------------------------ #

    def describe(self) -> str:
        """Prose summary of the registered frames and declared relationships."""

        if not self._frames:
            return f"Frame graph '{self.name}' is empty: no frames are registered."
        frame_lines = "; ".join(
            f"'{frame.name}' ({frame.domain.value}: {'/'.join(frame.axes)})"
            for frame in self.frames()
        )
        if self._transforms:
            edge_lines = "; ".join(
                f"'{transform.source.name}' <-> '{transform.target.name}' "
                f"({transform.rotation_angle_deg:.4f} deg)"
                for transform in self._transforms
            )
            edge_text = (
                f" Declared relationships (each usable in both directions): {edge_lines}."
            )
        else:
            edge_text = " No relationships are declared yet, so no frames are connected."
        return (
            f"Frame graph '{self.name}' registers {len(self._frames)} frame(s): "
            f"{frame_lines}.{edge_text} Transforms between connected frames are resolved by "
            "composing the shortest declared chain."
        )

    def __repr__(self) -> str:  # pragma: no cover - debugging convenience
        return (
            f"FrameGraph(name={self.name!r}, frames={len(self._frames)}, "
            f"transforms={len(self._transforms)})"
        )
