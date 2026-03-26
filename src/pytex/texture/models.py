from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector, normalize_vectors
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.projections import project_directions


@dataclass(frozen=True, slots=True)
class KernelSpec:
    name: str = "de_la_vallee_poussin"
    halfwidth_deg: float = 10.0

    def __post_init__(self) -> None:
        if self.name not in {"de_la_vallee_poussin", "von_mises_fisher"}:
            raise ValueError(
                "Kernel name must be either 'de_la_vallee_poussin' or 'von_mises_fisher'."
            )
        if self.halfwidth_deg <= 0.0:
            raise ValueError("Kernel halfwidth must be strictly positive.")

    def evaluate(self, angles_rad: ArrayLike) -> np.ndarray:
        angle_array = np.asarray(angles_rad, dtype=np.float64)
        halfwidth_rad = np.deg2rad(self.halfwidth_deg)
        if self.name == "de_la_vallee_poussin":
            denominator = np.cos(halfwidth_rad / 2.0)
            exponent = (
                1.0
                if np.isclose(denominator, 1.0)
                else float(np.log(0.5) / np.log(denominator))
            )
            exponent = max(1.0, exponent)
            values = np.clip(np.cos(angle_array / 2.0), 0.0, 1.0) ** exponent
        elif self.name == "von_mises_fisher":
            kappa = float(np.log(2.0) / (1.0 - np.cos(halfwidth_rad)))
            values = np.exp(kappa * (np.cos(angle_array) - 1.0))
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values


@dataclass(frozen=True, slots=True)
class PoleFigure:
    pole: CrystalPlane
    sample_directions: np.ndarray
    intensities: np.ndarray
    specimen_frame: ReferenceFrame
    antipodal: bool = True
    sample_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        sample_directions = normalize_vectors(self.sample_directions)
        intensities = as_float_array(self.intensities, shape=(sample_directions.shape[0],))
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError("PoleFigure.specimen_frame must belong to the specimen domain.")
        if np.any(~np.isfinite(intensities)) or np.any(intensities < 0.0):
            raise ValueError("PoleFigure intensities must be finite and non-negative.")
        if (
            self.sample_symmetry is not None
            and self.sample_symmetry.reference_frame != self.specimen_frame
        ):
            raise ValueError(
                "PoleFigure.sample_symmetry.reference_frame must match specimen_frame."
            )
        object.__setattr__(self, "sample_directions", sample_directions)
        object.__setattr__(self, "intensities", intensities)

    @classmethod
    def from_orientations(
        cls,
        orientations: OrientationSet,
        pole: CrystalPlane,
        *,
        weights: ArrayLike | None = None,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
        sample_symmetry: SymmetrySpec | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigure:
        if orientations.crystal_frame != pole.phase.crystal_frame:
            raise ValueError("PoleFigure orientations must use the pole phase crystal frame.")
        if orientations.phase is not None and orientations.phase != pole.phase:
            raise ValueError("PoleFigure orientations and pole must reference the same phase.")
        if orientations.symmetry is not None and orientations.symmetry != pole.phase.symmetry:
            raise ValueError("PoleFigure orientations symmetry must match the pole phase symmetry.")
        intensities = (
            as_float_array(np.ones(len(orientations)), shape=(len(orientations),))
            if weights is None
            else as_float_array(weights, shape=(len(orientations),))
        )
        plane_normals = (
            pole.phase.symmetry.equivalent_vectors(pole.normal)
            if include_symmetry_family
            else pole.normal[None, :]
        )
        specimen_directions = [
            orientations.map_crystal_directions(plane_normal)
            for plane_normal in plane_normals
        ]
        sample_directions = np.vstack(specimen_directions)
        repeated_intensities = np.tile(intensities, len(plane_normals))
        return cls(
            pole=pole,
            sample_directions=sample_directions,
            intensities=repeated_intensities,
            specimen_frame=orientations.specimen_frame,
            antipodal=antipodal,
            sample_symmetry=sample_symmetry,
            provenance=orientations.provenance if provenance is None else provenance,
        )

    def project(self, *, method: str = "equal_area") -> np.ndarray:
        return project_directions(
            self.sample_directions,
            method=method,
            antipodal=self.antipodal,
        )

    def histogram(
        self,
        *,
        bins: int = 72,
        method: str = "equal_area",
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        projected = self.project(method=method)
        radius = np.sqrt(2.0) if method == "equal_area" else 1.0
        histogram, xedges, yedges = np.histogram2d(
            projected[:, 0],
            projected[:, 1],
            bins=bins,
            range=[[-radius, radius], [-radius, radius]],
            weights=self.intensities,
        )
        histogram.setflags(write=False)
        xedges.setflags(write=False)
        yedges.setflags(write=False)
        return histogram, xedges, yedges


@dataclass(frozen=True, slots=True)
class InversePoleFigure:
    sample_direction: np.ndarray
    crystal_directions: np.ndarray
    intensities: np.ndarray
    crystal_frame: ReferenceFrame
    specimen_frame: ReferenceFrame
    antipodal: bool = True
    crystal_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        sample_direction = normalize_vector(self.sample_direction)
        crystal_directions = normalize_vectors(self.crystal_directions)
        intensities = as_float_array(self.intensities, shape=(crystal_directions.shape[0],))
        if self.crystal_frame.domain is not FrameDomain.CRYSTAL:
            raise ValueError("InversePoleFigure.crystal_frame must belong to the crystal domain.")
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError("InversePoleFigure.specimen_frame must belong to the specimen domain.")
        if np.any(~np.isfinite(intensities)) or np.any(intensities < 0.0):
            raise ValueError("InversePoleFigure intensities must be finite and non-negative.")
        if (
            self.crystal_symmetry is not None
            and self.crystal_symmetry.reference_frame != self.crystal_frame
        ):
            raise ValueError(
                "InversePoleFigure.crystal_symmetry.reference_frame must match crystal_frame."
            )
        object.__setattr__(self, "sample_direction", sample_direction)
        object.__setattr__(self, "crystal_directions", crystal_directions)
        object.__setattr__(self, "intensities", intensities)

    @classmethod
    def from_orientations(
        cls,
        orientations: OrientationSet,
        sample_direction: ArrayLike,
        *,
        weights: ArrayLike | None = None,
        reduce_by_symmetry: bool = True,
        antipodal: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> InversePoleFigure:
        normalized_sample_direction = normalize_vector(sample_direction)
        intensities = (
            as_float_array(np.ones(len(orientations)), shape=(len(orientations),))
            if weights is None
            else as_float_array(weights, shape=(len(orientations),))
        )
        crystal_directions = orientations.map_sample_directions_to_crystal(
            normalized_sample_direction
        )
        if reduce_by_symmetry and orientations.symmetry is not None:
            crystal_directions = orientations.symmetry.reduce_vectors_to_fundamental_sector(
                crystal_directions,
                antipodal=antipodal,
            )
        return cls(
            sample_direction=normalized_sample_direction,
            crystal_directions=crystal_directions,
            intensities=intensities,
            crystal_frame=orientations.crystal_frame,
            specimen_frame=orientations.specimen_frame,
            antipodal=antipodal,
            crystal_symmetry=orientations.symmetry,
            provenance=orientations.provenance if provenance is None else provenance,
        )

    def project(self, *, method: str = "equal_area") -> np.ndarray:
        return project_directions(
            self.crystal_directions,
            method=method,
            antipodal=self.antipodal,
        )

    @property
    def sector_vertices(self) -> np.ndarray | None:
        if self.crystal_symmetry is None:
            return None
        return self.crystal_symmetry.fundamental_sector(antipodal=self.antipodal).vertices

    def project_sector_vertices(self, *, method: str = "equal_area") -> np.ndarray | None:
        if self.sector_vertices is None:
            return None
        return project_directions(
            self.sector_vertices,
            method=method,
            antipodal=self.antipodal,
        )


@dataclass(frozen=True, slots=True)
class ODF:
    orientations: OrientationSet
    weights: np.ndarray
    kernel: KernelSpec = field(default_factory=KernelSpec)
    specimen_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        weights = as_float_array(self.weights, shape=(len(self.orientations),))
        if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("ODF weights must be finite and non-negative.")
        if np.isclose(float(np.sum(weights)), 0.0):
            raise ValueError("ODF weights must sum to a positive value.")
        if (
            self.specimen_symmetry is not None
            and self.specimen_symmetry.reference_frame != self.orientations.specimen_frame
        ):
            raise ValueError(
                "ODF.specimen_symmetry.reference_frame must match the "
                "OrientationSet specimen frame."
            )
        object.__setattr__(self, "weights", weights)

    @classmethod
    def from_orientations(
        cls,
        orientations: OrientationSet,
        *,
        weights: ArrayLike | None = None,
        kernel: KernelSpec | None = None,
        specimen_symmetry: SymmetrySpec | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> ODF:
        default_weights = np.ones(len(orientations), dtype=np.float64)
        weights_value = (
            default_weights
            if weights is None
            else np.asarray(weights, dtype=np.float64)
        )
        return cls(
            orientations=orientations,
            weights=weights_value,
            kernel=KernelSpec() if kernel is None else kernel,
            specimen_symmetry=specimen_symmetry,
            provenance=orientations.provenance if provenance is None else provenance,
        )

    @property
    def normalized_weights(self) -> np.ndarray:
        normalized = self.weights / np.sum(self.weights)
        normalized = np.ascontiguousarray(normalized)
        normalized.setflags(write=False)
        return normalized

    def evaluate(
        self,
        orientations: Orientation | OrientationSet,
        *,
        symmetry_aware: bool = True,
    ) -> np.ndarray | float:
        query_set: OrientationSet
        scalar_output = False
        if isinstance(orientations, Orientation):
            query_set = OrientationSet.from_orientations([orientations])
            scalar_output = True
        else:
            query_set = orientations
        if query_set.crystal_frame != self.orientations.crystal_frame:
            raise ValueError("ODF queries must use the same crystal frame as the ODF support.")
        if query_set.specimen_frame != self.orientations.specimen_frame:
            raise ValueError("ODF queries must use the same specimen frame as the ODF support.")
        if query_set.phase is not None and self.orientations.phase is not None:
            if query_set.phase != self.orientations.phase:
                raise ValueError("ODF queries must use the same phase as the ODF support.")
        if (
            symmetry_aware
            and query_set.symmetry is not None
            and self.orientations.symmetry is not None
            and query_set.symmetry != self.orientations.symmetry
        ):
            raise ValueError(
                "Symmetry-aware ODF queries must use the same symmetry as the support."
            )
        angles = query_set.misorientation_angles_to(
            self.orientations,
            symmetry_aware=symmetry_aware,
        )
        kernel_values = self.kernel.evaluate(angles)
        density = kernel_values @ self.normalized_weights
        density = np.ascontiguousarray(density)
        density.setflags(write=False)
        if scalar_output:
            return float(density[0])
        return density

    def reconstruct_pole_figure(
        self,
        pole: CrystalPlane,
        *,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
    ) -> PoleFigure:
        return PoleFigure.from_orientations(
            self.orientations,
            pole,
            weights=self.normalized_weights,
            include_symmetry_family=include_symmetry_family,
            antipodal=antipodal,
            sample_symmetry=self.specimen_symmetry,
            provenance=self.provenance,
        )

    def volume_fraction(
        self,
        center: Orientation,
        *,
        max_angle_deg: float,
        symmetry_aware: bool = True,
    ) -> float:
        query_set = OrientationSet.from_orientations([center])
        angles = query_set.misorientation_angles_to(
            self.orientations,
            symmetry_aware=symmetry_aware,
        )[0]
        mask = angles <= np.deg2rad(max_angle_deg)
        return float(np.sum(self.normalized_weights[mask]))
