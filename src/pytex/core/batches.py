from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_quaternions, normalize_vectors
from pytex.core.provenance import ProvenanceRecord

if TYPE_CHECKING:
    from pytex.core.frames import ReferenceFrame
    from pytex.core.orientation import Rotation


def normalize_euler_convention_name(convention: str) -> str:
    normalized = convention.strip().lower()
    aliases = {
        "bunge": "bunge",
        "bunge_zxz": "bunge",
        "zxz": "bunge",
        "matthies": "matthies",
        "matthies_zyz": "matthies",
        "abg": "abg",
        "abg_zyz": "abg",
        "zyz": "abg",
    }
    resolved = aliases.get(normalized)
    if resolved is None:
        supported = ", ".join(sorted(aliases))
        raise ValueError(
            f"Unsupported Euler convention '{convention}'. Supported conventions: {supported}"
        )
    return resolved


@dataclass(frozen=True, slots=True)
class VectorSet:
    """A batch of 3-vectors that knows which reference frame it lives in.

    Purpose
    -------
    The frame-carrying array primitive. A bare ``(n, 3)`` array cannot say
    whether its rows are crystal directions, specimen directions, or
    laboratory vectors, so operations that mix frames cannot be caught. This
    type makes the frame explicit, and the operations that consume it check
    it rather than assume it.

    Attributes
    ----------
    values : np.ndarray
        ``(n, 3)`` vectors; not normalized, so magnitudes are preserved.
    reference_frame : ReferenceFrame
        The domain-typed frame the vectors live in.
    provenance : ProvenanceRecord, optional
    """

    values: np.ndarray
    reference_frame: ReferenceFrame
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", as_float_array(self.values, shape=(None, 3)))

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | VectorSet:
        selected = self.values[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return VectorSet(
            values=selected,
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        """The underlying ``(n, 3)`` array, stripped of frame meaning.

        Use it only where the frame is already accounted for; passing the bare
        array across an API boundary discards the check that would have caught
        a frame mismatch.
        """

        return self.values

    def normalized(self) -> VectorSet:
        """The same vectors scaled to unit length, in the same frame.
        """

        return VectorSet(
            values=normalize_vectors(self.values),
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )

    def subset(self, indices: ArrayLike) -> VectorSet:
        """The vectors at the given indices, as a new set.

        Accepts integer indices or a boolean mask; frame and provenance are
        preserved.
        """

        return VectorSet(
            values=self.values[np.asarray(indices)],
            reference_frame=self.reference_frame,
            provenance=self.provenance,
        )


@dataclass(frozen=True, slots=True)
class EulerSet:
    """A batch of Euler-angle triples that carries its own convention.

    Purpose
    -------
    An angle triple is meaningless without knowing the convention: the same
    numbers denote different rotations under Bunge ZXZ and Matthies ZYZ. This
    type keeps the convention and the degrees/radians flag with the numbers,
    so a conversion cannot silently apply the wrong axis sequence.

    Attributes
    ----------
    angles : np.ndarray
        ``(n, 3)`` triples in the order the convention names them.
    convention : str
        ``"bunge"`` (default), ``"matthies"``, or ``"abg"``.
    degrees : bool
        Whether the angles are in degrees.
    provenance : ProvenanceRecord, optional
    """

    angles: np.ndarray
    convention: str = "bunge"
    degrees: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "angles", as_float_array(self.angles, shape=(None, 3)))
        object.__setattr__(self, "convention", normalize_euler_convention_name(self.convention))

    def __len__(self) -> int:
        return int(self.angles.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | EulerSet:
        selected = self.angles[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(3,))
        return EulerSet(
            angles=selected,
            convention=self.convention,
            degrees=self.degrees,
            provenance=self.provenance,
        )

    def as_array(self) -> np.ndarray:
        """The underlying ``(n, 3)`` angle array, without the convention metadata.

        Angles alone are ambiguous — the same triple means different rotations
        under Bunge and Matthies conventions — so prefer passing the
        ``EulerSet`` itself wherever the convention still matters.
        """

        return self.angles

    def subset(self, indices: ArrayLike) -> EulerSet:
        """The angle triples at the given indices, as a new set.

        Convention, degrees flag, and provenance are preserved.
        """

        return EulerSet(
            angles=self.angles[np.asarray(indices)],
            convention=self.convention,
            degrees=self.degrees,
            provenance=self.provenance,
        )

    def to_rotation_set(self) -> RotationSet:
        """Convert the angles to rotations under this set's declared convention.

        Because the convention travels with the angles, this conversion cannot
        silently apply the wrong axis sequence.
        """

        return RotationSet.from_euler_set(self)


@dataclass(frozen=True, slots=True)
class QuaternionSet:
    """A batch of unit quaternions in ``(w, x, y, z)`` order.

    Purpose
    -------
    The storage form for rotations, normalized on construction. It carries no
    rotation algebra of its own; convert with :meth:`to_rotation_set` to
    compose, apply, or convert representations.

    Attributes
    ----------
    quaternions : np.ndarray
        ``(n, 4)`` unit quaternions.
    provenance : ProvenanceRecord, optional
    """

    quaternions: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternions", normalize_quaternions(self.quaternions))

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    def __getitem__(self, index: Any) -> np.ndarray | QuaternionSet:
        selected = self.quaternions[index]
        if np.asarray(selected).ndim == 1:
            return as_float_array(selected, shape=(4,))
        return QuaternionSet(quaternions=selected, provenance=self.provenance)

    def as_array(self) -> np.ndarray:
        """The underlying ``(n, 4)`` quaternion array in ``(w, x, y, z)`` order.
        """

        return self.quaternions

    def subset(self, indices: ArrayLike) -> QuaternionSet:
        """The quaternions at the given indices, as a new set.
        """

        return QuaternionSet(
            quaternions=self.quaternions[np.asarray(indices)],
            provenance=self.provenance,
        )

    def to_rotation_set(self) -> RotationSet:
        """Reinterpret these quaternions as rotations.

        A ``QuaternionSet`` is storage; a ``RotationSet`` carries the rotation
        algebra (composition, application to vectors, Euler and Rodrigues
        conversions).
        """

        return RotationSet(quaternions=self.quaternions, provenance=self.provenance)


@dataclass(frozen=True, slots=True)
class RotationSet:
    """A batch of rotations with the full rotation algebra, frame-agnostic.

    Purpose
    -------
    The vectorized counterpart of :class:`~pytex.core.orientation.Rotation`:
    composition, application to vectors, and conversion between axis-angle,
    Euler, matrix, and Rodrigues forms, all as array operations. It carries
    no frame or symmetry meaning; for crystal orientations use
    :class:`~pytex.core.orientation.OrientationSet`, which adds them.

    Attributes
    ----------
    quaternions : np.ndarray
        ``(n, 4)`` unit quaternions in ``(w, x, y, z)`` order.
    provenance : ProvenanceRecord, optional
    """

    quaternions: np.ndarray
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "quaternions", normalize_quaternions(self.quaternions))

    def __len__(self) -> int:
        return int(self.quaternions.shape[0])

    def __getitem__(self, index: Any) -> Rotation | RotationSet:
        selected = self.quaternions[index]
        if np.asarray(selected).ndim == 1:
            from pytex.core.orientation import Rotation

            return Rotation(quaternion=selected, provenance=self.provenance)
        return RotationSet(quaternions=selected, provenance=self.provenance)

    @classmethod
    def from_rotations(cls, rotations: list[Rotation]) -> RotationSet:
        """Build a batch from individual :class:`~pytex.core.orientation.Rotation`
        objects.

        Use this to move from per-object work to the vectorized path before any
        calculation that scales with the number of rotations.
        """

        quaternions = np.stack([rotation.quaternion for rotation in rotations], axis=0)
        provenance = rotations[0].provenance if rotations else None
        return cls(quaternions=quaternions, provenance=provenance)

    @classmethod
    def from_euler_set(cls, euler_set: EulerSet) -> RotationSet:
        """Build a batch from a typed :class:`EulerSet`.

        The convention and degrees flag are taken from the ``EulerSet``, so the
        angles cannot be interpreted under the wrong convention.
        """

        from pytex.core.orientation import Rotation

        quaternions = np.stack(
            [
                Rotation.from_euler(
                    angle1,
                    angle2,
                    angle3,
                    convention=euler_set.convention,
                    degrees=euler_set.degrees,
                ).quaternion
                for angle1, angle2, angle3 in euler_set.angles
            ],
            axis=0,
        )
        return cls(quaternions=quaternions, provenance=euler_set.provenance)

    @classmethod
    def from_axes_angles(
        cls,
        axes: ArrayLike,
        angles_rad: ArrayLike,
        *,
        provenance: ProvenanceRecord | None = None,
    ) -> RotationSet:
        """Build a batch from ``(n, 3)`` rotation axes and ``(n,)`` angles.

        Angles are in radians; axes need not be normalized.
        """

        from pytex.core.orientation import quaternions_from_axes_angles

        return cls(
            quaternions=quaternions_from_axes_angles(axes, angles_rad),
            provenance=provenance,
        )

    @classmethod
    def from_rodrigues(
        cls,
        rodrigues: ArrayLike,
        *,
        frank: bool = False,
        provenance: ProvenanceRecord | None = None,
    ) -> RotationSet:
        """Build a batch from Rodrigues ``(n, 3)`` or Rodrigues-Frank ``(n, 4)`` rows.

        Pass ``frank=True`` for the homogeneous form, which stays finite at a
        rotation angle of ``pi``.
        """

        from pytex.core.orientation import quaternions_from_rodrigues

        return cls(
            quaternions=quaternions_from_rodrigues(rodrigues, frank=frank),
            provenance=provenance,
        )

    @classmethod
    def from_matrices(
        cls,
        matrices: ArrayLike,
        *,
        provenance: ProvenanceRecord | None = None,
    ) -> RotationSet:
        """Build a batch from an ``(n, 3, 3)`` stack of rotation matrices.

        Each matrix must be proper orthogonal.
        """

        from pytex.core.orientation import matrices_to_quaternions

        return cls(quaternions=matrices_to_quaternions(matrices), provenance=provenance)

    def as_quaternion_set(self) -> QuaternionSet:
        """The quaternion storage view of this batch, in ``(w, x, y, z)`` order.
        """

        return QuaternionSet(quaternions=self.quaternions, provenance=self.provenance)

    def as_matrices(self) -> np.ndarray:
        """``(n, 3, 3)`` rotation matrices, in the active convention ``v' = R v``.
        """

        from pytex.core.orientation import quaternions_to_matrices

        return quaternions_to_matrices(self.quaternions)

    def as_euler_set(
        self,
        *,
        convention: str = "bunge",
        degrees: bool = True,
    ) -> EulerSet:
        """Euler angles of the batch as a typed :class:`EulerSet`.

        Parameters
        ----------
        convention : str
            ``"bunge"`` (default), ``"matthies"``, or ``"abg"``.
        degrees : bool
            Emit degrees (default) rather than radians. Both are recorded on the
            returned set.
        """

        from pytex.core.orientation import Rotation

        angles = np.stack(
            [
                Rotation(quaternion=quaternion).to_euler(convention=convention, degrees=degrees)
                for quaternion in self.quaternions
            ],
            axis=0,
        )
        return EulerSet(
            angles=angles,
            convention=convention,
            degrees=degrees,
            provenance=self.provenance,
        )

    def to_axes_angles(self) -> tuple[np.ndarray, np.ndarray]:
        """Rotation axes and angles of the batch.

        Returns
        -------
        tuple of np.ndarray
            ``(n, 3)`` unit axes and ``(n,)`` angles in radians, in ``[0, pi]``.
        """

        from pytex.core.orientation import quaternions_to_axes_angles

        return quaternions_to_axes_angles(self.quaternions)

    def to_rodrigues(self, frank: bool = False) -> np.ndarray:
        """Rodrigues (``(n, 3)``) or Rodrigues-Frank (``(n, 4)``) parameters.

        Rodrigues space is where symmetry fundamental zones are convex
        polyhedra; the Frank form stays finite at a rotation angle of ``pi``.
        """

        from pytex.core.orientation import quaternions_to_rodrigues

        return quaternions_to_rodrigues(self.quaternions, frank=frank)

    def apply(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Rotate vectors by the batch, ``v -> R v``.

        Accepts one shared ``(3,)`` vector applied by every rotation, or an
        ``(n, 3)`` array applied row-wise. A ``VectorSet`` in gives a
        ``VectorSet`` out with its reference frame preserved — a bare rotation
        does not change which frame the vectors live in.
        """

        matrices = self.as_matrices()
        if isinstance(vectors, VectorSet):
            values = vectors.values
            if values.shape[0] != len(self):
                raise ValueError("VectorSet must have the same number of rows as the RotationSet.")
            mapped = np.einsum("nij,nj->ni", matrices, values, optimize=True)
            return VectorSet(
                values=mapped,
                reference_frame=vectors.reference_frame,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape == (3,):
            mapped = np.einsum("nij,j->ni", matrices, array, optimize=True)
        elif array.ndim == 2 and array.shape[1] == 3:
            if array.shape[0] != len(self):
                raise ValueError(
                    "Input vectors must have the same number of rows as the RotationSet."
                )
            mapped = np.einsum("nij,nj->ni", matrices, array, optimize=True)
        else:
            raise ValueError("Input vectors must have shape (3,) or (n, 3).")
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def apply_inverse(self, vectors: ArrayLike | VectorSet) -> np.ndarray | VectorSet:
        """Rotate vectors by the inverse rotations, ``v -> R^T v``.

        Same shape rules as :meth:`apply`.
        """

        inverse_matrices = np.swapaxes(self.as_matrices(), -1, -2)
        if isinstance(vectors, VectorSet):
            values = vectors.values
            if values.shape[0] != len(self):
                raise ValueError("VectorSet must have the same number of rows as the RotationSet.")
            mapped = np.einsum("nij,nj->ni", inverse_matrices, values, optimize=True)
            return VectorSet(
                values=mapped,
                reference_frame=vectors.reference_frame,
                provenance=vectors.provenance,
            )
        array = np.asarray(vectors, dtype=np.float64)
        if array.shape == (3,):
            mapped = np.einsum("nij,j->ni", inverse_matrices, array, optimize=True)
        elif array.ndim == 2 and array.shape[1] == 3:
            if array.shape[0] != len(self):
                raise ValueError(
                    "Input vectors must have the same number of rows as the RotationSet."
                )
            mapped = np.einsum("nij,nj->ni", inverse_matrices, array, optimize=True)
        else:
            raise ValueError("Input vectors must have shape (3,) or (n, 3).")
        mapped = np.ascontiguousarray(mapped)
        mapped.setflags(write=False)
        return mapped

    def subset(self, indices: ArrayLike) -> RotationSet:
        """The rotations at the given indices, as a new batch.
        """

        return RotationSet(
            quaternions=self.quaternions[np.asarray(indices)],
            provenance=self.provenance,
        )
