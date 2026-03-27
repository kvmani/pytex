from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector, normalize_vectors
from pytex.core.batches import VectorSet
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.projections import project_directions


def _as_direction_array(vectors: np.ndarray | VectorSet) -> np.ndarray:
    if isinstance(vectors, VectorSet):
        return vectors.values
    return vectors


def _orientation_dictionary_response(
    dictionary: OrientationSet,
    pole_figure: PoleFigure,
    kernel: KernelSpec,
    *,
    include_symmetry_family: bool,
) -> np.ndarray:
    if dictionary.crystal_frame != pole_figure.pole.phase.crystal_frame:
        raise ValueError("PoleFigure inversion dictionary must use the pole phase crystal frame.")
    if dictionary.phase is not None and dictionary.phase != pole_figure.pole.phase:
        raise ValueError("PoleFigure inversion dictionary phase must match PoleFigure.pole.phase.")
    if dictionary.symmetry is not None and dictionary.symmetry != pole_figure.pole.phase.symmetry:
        raise ValueError(
            "PoleFigure inversion dictionary symmetry must match PoleFigure.pole.phase.symmetry."
        )
    if dictionary.specimen_frame != pole_figure.specimen_frame:
        raise ValueError(
            "PoleFigure inversion dictionary specimen_frame must match PoleFigure.specimen_frame."
        )
    pole_family = (
        pole_figure.pole.phase.symmetry.equivalent_vectors(pole_figure.pole.normal)
        if include_symmetry_family
        else pole_figure.pole.normal[None, :]
    )
    mapped_families = np.stack(
        [
            _as_direction_array(dictionary.map_crystal_directions(direction))
            for direction in pole_family
        ],
        axis=1,
    )
    cos_angles = np.einsum(
        "mk,nfk->mnf",
        pole_figure.sample_directions,
        mapped_families,
        optimize=True,
    )
    angles = np.arccos(np.clip(cos_angles, -1.0, 1.0))
    response = kernel.evaluate(angles)
    block = np.mean(response, axis=2)
    block = np.ascontiguousarray(block)
    block.setflags(write=False)
    return block


def _projected_gradient_nonnegative_weights(
    system_matrix: np.ndarray,
    observations: np.ndarray,
    *,
    regularization: float,
    max_iterations: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, bool]:
    if regularization < 0.0:
        raise ValueError("regularization must be non-negative.")
    if max_iterations <= 0:
        raise ValueError("max_iterations must be strictly positive.")
    if tolerance <= 0.0:
        raise ValueError("tolerance must be strictly positive.")
    weights = np.full(system_matrix.shape[1], 1.0 / system_matrix.shape[1], dtype=np.float64)
    gram = system_matrix.T @ system_matrix
    rhs = system_matrix.T @ observations
    lipschitz = float(np.linalg.norm(gram + regularization * np.eye(gram.shape[0]), ord=2))
    lipschitz = max(lipschitz, 1e-12)
    history = np.empty(max_iterations, dtype=np.float64)
    converged = False
    for iteration in range(max_iterations):
        gradient = gram @ weights - rhs + regularization * weights
        candidate = np.maximum(weights - gradient / lipschitz, 0.0)
        total = float(np.sum(candidate))
        if np.isclose(total, 0.0):
            candidate = np.full_like(candidate, 1.0 / candidate.size)
        else:
            candidate /= total
        residual = system_matrix @ candidate - observations
        history[iteration] = 0.5 * float(residual @ residual) + 0.5 * regularization * float(
            candidate @ candidate
        )
        delta = np.linalg.norm(candidate - weights)
        scale = max(1.0, float(np.linalg.norm(weights)))
        weights = candidate
        if delta <= tolerance * scale:
            history = history[: iteration + 1]
            converged = True
            break
    else:
        history = history[:max_iterations]
    weights = np.ascontiguousarray(weights)
    weights.setflags(write=False)
    history = np.ascontiguousarray(history)
    history.setflags(write=False)
    return weights, history, converged


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
                1.0 if np.isclose(denominator, 1.0) else float(np.log(0.5) / np.log(denominator))
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
            orientations.map_crystal_directions(plane_normal) for plane_normal in plane_normals
        ]
        sample_directions = np.vstack(
            [
                direction.values if isinstance(direction, VectorSet) else direction
                for direction in specimen_directions
            ]
        )
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
        crystal_direction_array = (
            crystal_directions.values
            if isinstance(crystal_directions, VectorSet)
            else crystal_directions
        )
        return cls(
            sample_direction=normalized_sample_direction,
            crystal_directions=crystal_direction_array,
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
            default_weights if weights is None else np.asarray(weights, dtype=np.float64)
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

    def reconstruct_pole_figures(
        self,
        poles: Sequence[CrystalPlane],
        *,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
    ) -> tuple[PoleFigure, ...]:
        return tuple(
            self.reconstruct_pole_figure(
                pole,
                include_symmetry_family=include_symmetry_family,
                antipodal=antipodal,
            )
            for pole in poles
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

    @classmethod
    def invert_pole_figures(
        cls,
        pole_figures: Sequence[PoleFigure],
        *,
        orientation_dictionary: OrientationSet,
        kernel: KernelSpec | None = None,
        regularization: float = 1e-6,
        include_symmetry_family: bool = True,
        max_iterations: int = 500,
        tolerance: float = 1e-8,
        provenance: ProvenanceRecord | None = None,
    ) -> ODFInversionReport:
        if not pole_figures:
            raise ValueError("ODF inversion requires at least one PoleFigure.")
        specimen_frame = orientation_dictionary.specimen_frame
        for pole_figure in pole_figures:
            if pole_figure.specimen_frame != specimen_frame:
                raise ValueError(
                    "All pole figures and the inversion dictionary must share a specimen frame."
                )
        inversion_kernel = KernelSpec() if kernel is None else kernel
        blocks = [
            _orientation_dictionary_response(
                orientation_dictionary,
                pole_figure,
                inversion_kernel,
                include_symmetry_family=include_symmetry_family,
            )
            for pole_figure in pole_figures
        ]
        system_matrix = np.vstack(blocks)
        observations = np.concatenate([pole_figure.intensities for pole_figure in pole_figures])
        if system_matrix.shape[0] != observations.shape[0]:
            raise ValueError("ODF inversion system matrix and observation vector are inconsistent.")
        weights, objective_history, converged = _projected_gradient_nonnegative_weights(
            system_matrix,
            observations,
            regularization=regularization,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )
        common_sample_symmetry = pole_figures[0].sample_symmetry
        if any(
            pole_figure.sample_symmetry != common_sample_symmetry
            for pole_figure in pole_figures
        ):
            common_sample_symmetry = None
        odf = cls(
            orientations=orientation_dictionary,
            weights=weights,
            kernel=inversion_kernel,
            specimen_symmetry=common_sample_symmetry,
            provenance=orientation_dictionary.provenance if provenance is None else provenance,
        )
        return ODFInversionReport(
            odf=odf,
            residual_norm=float(np.linalg.norm(system_matrix @ weights - observations)),
            objective_history=objective_history,
            iterations=int(objective_history.shape[0]),
            converged=converged,
            regularization=regularization,
            observation_count=int(observations.size),
            dictionary_size=int(len(orientation_dictionary)),
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class ODFInversionReport:
    odf: ODF
    residual_norm: float
    objective_history: np.ndarray
    iterations: int
    converged: bool
    regularization: float
    observation_count: int
    dictionary_size: int
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        history = as_float_array(self.objective_history, shape=(None,))
        if history.size == 0:
            raise ValueError("ODFInversionReport.objective_history must not be empty.")
        if self.residual_norm < 0.0:
            raise ValueError("ODFInversionReport.residual_norm must be non-negative.")
        if self.iterations <= 0:
            raise ValueError("ODFInversionReport.iterations must be strictly positive.")
        if self.regularization < 0.0:
            raise ValueError("ODFInversionReport.regularization must be non-negative.")
        if self.observation_count <= 0:
            raise ValueError("ODFInversionReport.observation_count must be strictly positive.")
        if self.dictionary_size <= 0:
            raise ValueError("ODFInversionReport.dictionary_size must be strictly positive.")
        object.__setattr__(self, "objective_history", history)
