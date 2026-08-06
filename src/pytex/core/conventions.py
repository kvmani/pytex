from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class BasisKind(StrEnum):
    """Whether components are expressed in the direct or reciprocal basis.

    Outside the cubic system ``[hkl]`` and ``(hkl)`` denote different
    geometric objects, so every stable crystallographic vector representation
    in PyTex declares which basis its components belong to rather than leaving
    it to context.
    """

    DIRECT = "direct"
    RECIPROCAL = "reciprocal"


class FrameDomain(StrEnum):
    """The six reference-frame domains PyTex recognizes.

    The canonical frame chain is
    ``crystal -> specimen -> map -> detector -> laboratory -> reciprocal``.
    Every :class:`~pytex.core.frames.ReferenceFrame` is typed by one of these
    domains, and no subsystem may invent another or silently collapse two of
    them. See ``docs/standards/notation_and_conventions.md``.
    """

    CRYSTAL = "crystal"
    SPECIMEN = "specimen"
    MAP = "map"
    DETECTOR = "detector"
    LABORATORY = "laboratory"
    RECIPROCAL = "reciprocal"


class Handedness(StrEnum):
    """Chirality of a coordinate frame.

    PyTex frames are right-handed throughout; the left-handed member exists so
    that an imported dataset can record the handedness its source system used
    rather than having that information dropped at the boundary.
    """

    RIGHT = "right"
    LEFT = "left"


class AngleConvention(StrEnum):
    """Supported Euler-angle conventions.

    ``BUNGE_ZXZ`` is the ZXZ convention of Bunge, labelled
    ``(phi1, Phi, phi2)``, and is the PyTex default and the texture-community
    standard. ``MATTHIES_ZYZ`` and ``ABG_ZYZ`` are ZYZ conventions labelled
    ``(alpha, beta, gamma)``. The same angle triple denotes different rotations
    under different conventions, which is why the convention travels with the
    angles in :class:`~pytex.core.batches.EulerSet`.
    """

    BUNGE_ZXZ = "bunge_zxz"
    MATTHIES_ZYZ = "matthies_zyz"
    ABG_ZYZ = "abg_zyz"


class ReciprocalConvention(StrEnum):
    KRONECKER_DELTA = "kronecker_delta"


@dataclass(frozen=True, slots=True)
class ConventionSet:
    """The canonical conventions a body of data or results was produced under.

    Purpose
    -------
    Records the Euler convention, quaternion storage order, and reciprocal
    normalization as data rather than as documentation, so that manifests and
    interchange payloads state their conventions explicitly instead of
    relying on a reader knowing them.

    The repository-wide instance is ``PYTEX_CANONICAL_CONVENTIONS``.

    Attributes
    ----------
    name : str
    angle_convention : AngleConvention
        Bunge ZXZ by default.
    quaternion_order : tuple of str
        ``("w", "x", "y", "z")``.
    reciprocal_convention : ReciprocalConvention
        Kronecker-delta normalization, i.e. no ``2*pi`` factor.
    notes : tuple of str
        Human-readable statements of the same conventions.
    """

    name: str
    angle_convention: AngleConvention = AngleConvention.BUNGE_ZXZ
    quaternion_order: tuple[str, str, str, str] = ("w", "x", "y", "z")
    reciprocal_convention: ReciprocalConvention = ReciprocalConvention.KRONECKER_DELTA
    notes: tuple[str, ...] = field(
        default_factory=lambda: (
            "Right-handed Cartesian frames.",
            "Bunge phi1, Phi, phi2 labeling.",
            "Unit quaternions stored in w, x, y, z order.",
            "Reciprocal basis normalized such that a*_i dot a_j = delta_ij.",
        )
    )


PYTEX_CANONICAL_CONVENTIONS = ConventionSet(name="pytex_canonical")
