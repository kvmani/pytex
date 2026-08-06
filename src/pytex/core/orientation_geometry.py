from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import FundamentalSector, SymmetrySpec

EulerConventionName = Literal["bunge", "roe", "kocks", "abg", "matthies"]


def _wrap_degrees(values: np.ndarray) -> np.ndarray:
    wrapped = np.mod(values, 360.0)
    wrapped[np.isclose(wrapped, 360.0)] = 0.0
    return np.ascontiguousarray(wrapped, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class EulerConventionTransform:
    """A conversion between two named Euler-angle conventions.

    Purpose
    -------
    Data arrives in whichever convention its source system used — Bunge, Roe,
    Kocks, Matthies, ABG. This converts between them through the rotation
    itself rather than by an angle-level formula, so the result stays correct
    at the gimbal-degenerate configurations where such formulae break down.

    Attributes
    ----------
    source, target : str
        Convention names: ``"bunge"``, ``"roe"``, ``"kocks"``, ``"abg"``, or
        ``"matthies"``.
    degrees : bool
        Whether the angles are in degrees.
    provenance : ProvenanceRecord, optional
    """

    source: EulerConventionName
    target: EulerConventionName
    degrees: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        allowed = {"bunge", "roe", "kocks", "abg", "matthies"}
        if self.source not in allowed or self.target not in allowed:
            raise ValueError("Unsupported Euler convention transform.")

    def apply(self, angles: ArrayLike) -> np.ndarray:
        """Convert Euler angles between the source and target conventions.

        Purpose
        -------
        The same rotation is described by different angle triples in Bunge ZXZ
        and Matthies/ABG ZYZ conventions. This converts through the rotation
        itself rather than by an angle-level formula, so it stays correct at the
        gimbal-degenerate configurations where such formulae break down.

        Parameters
        ----------
        angles : ArrayLike
            A ``(3,)`` triple or an ``(n, 3)`` array in the source convention.

        Returns
        -------
        np.ndarray
            The same rotations expressed in the target convention.
        """

        values = np.asarray(angles, dtype=np.float64)
        scalar = values.shape == (3,)
        rows = values.reshape(1, 3) if scalar else as_float_array(values, shape=(None, 3))
        working = np.rad2deg(rows) if not self.degrees else rows.copy()
        bunge = self._to_bunge(working)
        transformed = self._from_bunge(bunge)
        if not self.degrees:
            transformed = np.deg2rad(transformed)
        transformed = np.ascontiguousarray(transformed, dtype=np.float64)
        transformed.setflags(write=False)
        return transformed[0] if scalar else transformed

    def _to_bunge(self, angles_deg: np.ndarray) -> np.ndarray:
        if self.source == "bunge":
            return _wrap_degrees(angles_deg)
        if self.source in {"roe", "abg", "matthies"}:
            converted = np.column_stack(
                [angles_deg[:, 0] + 90.0, angles_deg[:, 1], angles_deg[:, 2] - 90.0]
            )
            return _wrap_degrees(converted)
        converted = np.column_stack(
            [angles_deg[:, 0] + 90.0, angles_deg[:, 1], 90.0 - angles_deg[:, 2]]
        )
        return _wrap_degrees(converted)

    def _from_bunge(self, bunge_deg: np.ndarray) -> np.ndarray:
        if self.target == "bunge":
            return _wrap_degrees(bunge_deg)
        if self.target in {"roe", "abg", "matthies"}:
            converted = np.column_stack(
                [bunge_deg[:, 0] - 90.0, bunge_deg[:, 1], bunge_deg[:, 2] + 90.0]
            )
            return _wrap_degrees(converted)
        converted = np.column_stack(
            [bunge_deg[:, 0] - 90.0, bunge_deg[:, 1], 90.0 - bunge_deg[:, 2]]
        )
        return _wrap_degrees(converted)


@dataclass(frozen=True, slots=True)
class IPFSectorBoundary:
    """The fundamental-sector boundary of an inverse pole figure.

    Purpose
    -------
    The standard stereographic triangle for a given crystal symmetry: the
    region containing one representative of each direction orbit, its corner
    vertices, and the reduction that folds arbitrary directions into it. This
    is what makes IPF colours encode crystallography rather than the
    arbitrary symmetry branch a measurement reported.
    """

    symmetry: SymmetrySpec
    sector: FundamentalSector
    boundary_equations: tuple[str, ...]
    provenance: ProvenanceRecord | None = None

    @classmethod
    def from_symmetry(
        cls,
        symmetry: SymmetrySpec,
        *,
        antipodal: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> IPFSectorBoundary:
        """The inverse-pole-figure sector boundary implied by a crystal symmetry.

        The boundary of the standard stereographic triangle — the region
        containing exactly one representative of each direction orbit, which is
        the domain an inverse pole figure is drawn in.
        """

        sector = symmetry.fundamental_sector(antipodal=antipodal)
        equations = _boundary_equations_for_group(sector.proper_point_group, antipodal=antipodal)
        return cls(
            symmetry=symmetry,
            sector=sector,
            boundary_equations=equations,
            provenance=provenance or symmetry.provenance,
        )

    @property
    def vertices(self) -> np.ndarray:
        """Corner directions of the sector, as unit vectors.

        For cubic symmetry these are the familiar ``[001]``, ``[101]``, ``[111]``
        corners of the standard triangle.
        """

        return self.sector.vertices

    def contains(self, direction: ArrayLike) -> bool:
        """Whether a direction lies inside the sector.

        Boundary directions are accepted within tolerance, so a pole lying
        exactly on an edge is not rejected by floating-point noise.
        """

        return self.symmetry.vector_in_fundamental_sector(
            normalize_vector(direction),
            antipodal=self.sector.antipodal,
        )

    def reduce(self, direction: ArrayLike) -> np.ndarray:
        """Fold a direction into the sector by crystal symmetry.

        The operation behind IPF colouring: symmetry-equivalent directions map to
        the same point of the standard triangle, so the colour encodes the
        crystal direction rather than the arbitrary symmetry branch a
        measurement happened to report.
        """

        return self.symmetry.reduce_vector_to_fundamental_sector(
            direction,
            antipodal=self.sector.antipodal,
        )


def _boundary_equations_for_group(proper_point_group: str, *, antipodal: bool) -> tuple[str, ...]:
    hemisphere = ("z >= 0",) if antipodal else ()
    if proper_point_group in {"23", "432"}:
        return (*hemisphere, "y >= 0", "x >= y", "z >= x")
    if proper_point_group in {"4", "422"}:
        return (*hemisphere, "y >= 0", "x >= y")
    if proper_point_group in {"3", "32"}:
        return (*hemisphere, "x >= 0", "0 <= atan2(y, x) <= 60 deg")
    if proper_point_group in {"6", "622"}:
        return (*hemisphere, "x >= 0", "0 <= atan2(y, x) <= 30 deg")
    if proper_point_group in {"2", "222"}:
        return ("x >= 0", "y >= 0", "z >= 0")
    return (*hemisphere,)


@dataclass(frozen=True, slots=True)
class OrientationFundamentalRegion:
    """The fundamental region of orientation space for a crystal symmetry.

    Purpose
    -------
    The orientation-space counterpart of :class:`IPFSectorBoundary`: the
    region containing one representative of each symmetry-equivalent
    orientation. Reducing into it is the precondition for comparing,
    averaging, or plotting orientations without the symmetry branch
    contaminating the result.
    """

    crystal_symmetry: SymmetrySpec
    specimen_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def reduce(self, orientation: Orientation) -> Orientation:
        """The representative of an orientation inside the fundamental region.

        The orientation-space counterpart of :meth:`IPFSectorBoundary.reduce`;
        see :meth:`~pytex.core.orientation.Orientation.project_to_fundamental_region`.
        """

        if orientation.symmetry is not None and orientation.symmetry != self.crystal_symmetry:
            raise ValueError("orientation.symmetry must match crystal_symmetry.")
        reduced = orientation.project_to_exact_fundamental_region(
            specimen_symmetry=self.specimen_symmetry
        )
        if self.specimen_symmetry is None:
            return reduced
        candidates = self.specimen_symmetry.apply_to_rotation_matrices(
            reduced.rotation.as_matrix(),
            side="left",
        ).reshape(-1, 3, 3)
        candidate_set = OrientationSet.from_matrices(
            candidates,
            crystal_frame=reduced.crystal_frame,
            specimen_frame=reduced.specimen_frame,
            symmetry=self.crystal_symmetry,
            phase=reduced.phase,
            provenance=self.provenance or reduced.provenance,
        )
        identity = Orientation(
            rotation=Rotation.identity(),
            crystal_frame=reduced.crystal_frame,
            specimen_frame=reduced.specimen_frame,
            symmetry=self.crystal_symmetry,
            phase=reduced.phase,
        )
        identity_set = OrientationSet.from_orientations([identity])
        index = int(np.argmin(candidate_set.misorientation_angles_to(identity_set)[:, 0]))
        return candidate_set[index]

    def contains(self, orientation: Orientation, *, tolerance_rad: float = 1e-8) -> bool:
        """Whether an orientation already lies in the fundamental region.

        Parameters
        ----------
        orientation : Orientation
        tolerance_rad : float
            Angular tolerance, so an orientation on a region boundary is not
            rejected by numerical noise.
        """

        return orientation.is_in_fundamental_region(
            specimen_symmetry=self.specimen_symmetry,
            atol=tolerance_rad,
        )


__all__ = [
    "EulerConventionTransform",
    "IPFSectorBoundary",
    "OrientationFundamentalRegion",
]
