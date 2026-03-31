from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class BasisKind(StrEnum):
    DIRECT = "direct"
    RECIPROCAL = "reciprocal"


class FrameDomain(StrEnum):
    CRYSTAL = "crystal"
    SPECIMEN = "specimen"
    MAP = "map"
    DETECTOR = "detector"
    LABORATORY = "laboratory"
    RECIPROCAL = "reciprocal"


class Handedness(StrEnum):
    RIGHT = "right"
    LEFT = "left"


class AngleConvention(StrEnum):
    BUNGE_ZXZ = "bunge_zxz"
    MATTHIES_ZYZ = "matthies_zyz"
    ABG_ZYZ = "abg_zyz"


class ReciprocalConvention(StrEnum):
    KRONECKER_DELTA = "kronecker_delta"


@dataclass(frozen=True, slots=True)
class ConventionSet:
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
