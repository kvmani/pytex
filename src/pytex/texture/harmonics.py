from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import factorial

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, Phase
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.models import KernelSpec, PoleFigure, _pole_density_response_matrix


def _identity_operators() -> np.ndarray:
    operators = np.eye(3, dtype=np.float64)[None, :, :]
    operators.setflags(write=False)
    return operators


def _normalized_weights(values: ArrayLike) -> np.ndarray:
    weights = np.asarray(values, dtype=np.float64)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("Quadrature weights must sum to a positive finite value.")
    normalized = np.ascontiguousarray(weights / total, dtype=np.float64)
    normalized.setflags(write=False)
    return normalized


def _midpoint_axis_values(start_deg: float, stop_deg: float, step_deg: float) -> np.ndarray:
    if step_deg <= 0.0:
        raise ValueError("Quadrature steps must be strictly positive.")
    span = stop_deg - start_deg
    count = int(round(span / step_deg))
    if count <= 0 or not np.isclose(count * step_deg, span, atol=1e-8):
        raise ValueError("Quadrature step must partition the requested angular range exactly.")
    centers = start_deg + (np.arange(count, dtype=np.float64) + 0.5) * step_deg
    centers = np.ascontiguousarray(centers, dtype=np.float64)
    centers.setflags(write=False)
    return centers


def _bunge_quadrature(
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    crystal_symmetry: SymmetrySpec | None,
    phase: Phase | None,
    phi1_step_deg: float,
    big_phi_step_deg: float,
    phi2_step_deg: float,
    provenance: ProvenanceRecord | None,
) -> tuple[OrientationSet, np.ndarray]:
    phi1_values = _midpoint_axis_values(0.0, 360.0, phi1_step_deg)
    big_phi_values = _midpoint_axis_values(0.0, 180.0, big_phi_step_deg)
    phi2_values = _midpoint_axis_values(0.0, 360.0, phi2_step_deg)
    phi1_mesh, big_phi_mesh, phi2_mesh = np.meshgrid(
        phi1_values,
        big_phi_values,
        phi2_values,
        indexing="ij",
    )
    angles_deg = np.column_stack(
        [
            phi1_mesh.reshape(-1),
            big_phi_mesh.reshape(-1),
            phi2_mesh.reshape(-1),
        ]
    )
    orientations = OrientationSet.from_euler_angles(
        angles_deg,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        symmetry=crystal_symmetry,
        phase=phase,
        convention="bunge",
        degrees=True,
        provenance=provenance,
    )
    big_phi_rad = np.deg2rad(angles_deg[:, 1])
    raw_weights = np.sin(big_phi_rad)
    return orientations, _normalized_weights(raw_weights)


@dataclass(frozen=True, slots=True)
class HarmonicBasisTerm:
    degree: int
    sample_order: int
    crystal_order: int
    component: str

    def __post_init__(self) -> None:
        if self.degree < 0:
            raise ValueError("HarmonicBasisTerm.degree must be non-negative.")
        if abs(self.sample_order) > self.degree:
            raise ValueError("sample_order magnitude must not exceed degree.")
        if abs(self.crystal_order) > self.degree:
            raise ValueError("crystal_order magnitude must not exceed degree.")
        if self.component not in {"real", "imag"}:
            raise ValueError("component must be either 'real' or 'imag'.")


def _enumerate_terms(*, degree_bandlimit: int, even_degrees_only: bool) -> tuple[HarmonicBasisTerm, ...]:
    if degree_bandlimit < 0:
        raise ValueError("degree_bandlimit must be non-negative.")
    terms: list[HarmonicBasisTerm] = []
    for degree in range(degree_bandlimit + 1):
        if even_degrees_only and degree % 2:
            continue
        for sample_order in range(-degree, degree + 1):
            for crystal_order in range(-degree, degree + 1):
                terms.append(
                    HarmonicBasisTerm(
                        degree=degree,
                        sample_order=sample_order,
                        crystal_order=crystal_order,
                        component="real",
                    )
                )
                if sample_order != 0 or crystal_order != 0:
                    terms.append(
                        HarmonicBasisTerm(
                            degree=degree,
                            sample_order=sample_order,
                            crystal_order=crystal_order,
                            component="imag",
                        )
                    )
    return tuple(terms)


def _wigner_small_d(degree: int, sample_order: int, crystal_order: int, beta_rad: np.ndarray) -> np.ndarray:
    prefactor = np.sqrt(
        factorial(degree + sample_order)
        * factorial(degree - sample_order)
        * factorial(degree + crystal_order)
        * factorial(degree - crystal_order)
    )
    k_min = max(0, crystal_order - sample_order)
    k_max = min(degree - sample_order, degree + crystal_order)
    cos_half = np.cos(beta_rad / 2.0)
    sin_half = np.sin(beta_rad / 2.0)
    values = np.zeros_like(beta_rad, dtype=np.float64)
    for k in range(k_min, k_max + 1):
        denominator = (
            factorial(degree + crystal_order - k)
            * factorial(k)
            * factorial(sample_order - crystal_order + k)
            * factorial(degree - sample_order - k)
        )
        exponent_cos = 2 * degree + crystal_order - sample_order - 2 * k
        exponent_sin = sample_order - crystal_order + 2 * k
        coefficient = ((-1) ** (k - sample_order + crystal_order)) * prefactor / denominator
        values += coefficient * (cos_half**exponent_cos) * (sin_half**exponent_sin)
    values = np.ascontiguousarray(values, dtype=np.float64)
    values.setflags(write=False)
    return values


def _evaluate_raw_terms(angles_rad: np.ndarray, terms: Sequence[HarmonicBasisTerm]) -> np.ndarray:
    phi1 = angles_rad[:, 0]
    big_phi = angles_rad[:, 1]
    phi2 = angles_rad[:, 2]
    columns: list[np.ndarray] = []
    for term in terms:
        d_values = _wigner_small_d(term.degree, term.sample_order, term.crystal_order, big_phi)
        phase = term.sample_order * phi1 + term.crystal_order * phi2
        if term.component == "real":
            column = d_values * np.cos(phase)
        else:
            column = d_values * np.sin(phase)
        columns.append(np.asarray(column, dtype=np.float64))
    basis = np.column_stack(columns) if columns else np.zeros((angles_rad.shape[0], 0), dtype=np.float64)
    basis = np.ascontiguousarray(basis, dtype=np.float64)
    basis.setflags(write=False)
    return basis


def _symmetry_projected_raw_basis(
    orientations: OrientationSet,
    *,
    terms: Sequence[HarmonicBasisTerm],
    crystal_symmetry: SymmetrySpec | None,
    specimen_symmetry: SymmetrySpec | None,
) -> np.ndarray:
    matrices = orientations.as_matrices()
    specimen_operators = (
        _identity_operators() if specimen_symmetry is None else specimen_symmetry.operators
    )
    crystal_operators = (
        _identity_operators() if crystal_symmetry is None else crystal_symmetry.operators
    )
    transformed = np.einsum(
        "sij,njk,ckl->sncil",
        specimen_operators,
        matrices,
        crystal_operators,
        optimize=True,
    )
    transformed_orientations = OrientationSet.from_matrices(
        transformed.reshape(-1, 3, 3),
        crystal_frame=orientations.crystal_frame,
        specimen_frame=orientations.specimen_frame,
        symmetry=crystal_symmetry,
        phase=orientations.phase,
        provenance=orientations.provenance,
    )
    angles_rad = transformed_orientations.as_euler_set(convention="bunge", degrees=False).angles
    raw_basis = _evaluate_raw_terms(angles_rad, terms)
    projected = raw_basis.reshape(
        specimen_operators.shape[0],
        len(orientations),
        crystal_operators.shape[0],
        raw_basis.shape[1],
    ).mean(axis=(0, 2))
    projected = np.ascontiguousarray(projected, dtype=np.float64)
    projected.setflags(write=False)
    return projected


def _orthonormalize_weighted_basis(
    raw_basis: np.ndarray,
    quadrature_weights: np.ndarray,
    *,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray]:
    gram = raw_basis.T @ (quadrature_weights[:, None] * raw_basis)
    gram = 0.5 * (gram + gram.T)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    keep = eigenvalues > tolerance
    if not np.any(keep):
        raise ValueError(
            "The symmetry-projected harmonic basis is numerically rank-deficient at the requested bandlimit."
        )
    kept_values = eigenvalues[keep]
    kept_vectors = eigenvectors[:, keep]
    transform = kept_vectors / np.sqrt(kept_values)[None, :]
    orthonormal_basis = raw_basis @ transform
    orthonormal_basis = np.ascontiguousarray(orthonormal_basis, dtype=np.float64)
    orthonormal_basis.setflags(write=False)
    transform = np.ascontiguousarray(transform, dtype=np.float64)
    transform.setflags(write=False)
    return orthonormal_basis, transform


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(values * weights))


def _coerce_query_orientations(
    orientations: Orientation | OrientationSet,
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    crystal_symmetry: SymmetrySpec | None,
    phase: Phase | None,
    provenance: ProvenanceRecord | None,
) -> tuple[OrientationSet, bool]:
    if isinstance(orientations, Orientation):
        query_set = OrientationSet.from_orientations([orientations])
        scalar_output = True
    else:
        query_set = orientations
        scalar_output = False
    if query_set.crystal_frame != crystal_frame:
        raise ValueError("HarmonicODF queries must use the same crystal frame as the ODF.")
    if query_set.specimen_frame != specimen_frame:
        raise ValueError("HarmonicODF queries must use the same specimen frame as the ODF.")
    if phase is not None and query_set.phase is not None and query_set.phase != phase:
        raise ValueError("HarmonicODF queries must use the same phase as the ODF.")
    if query_set.phase is None and phase is not None:
        query_set = OrientationSet(
            quaternions=query_set.quaternions,
            crystal_frame=query_set.crystal_frame,
            specimen_frame=query_set.specimen_frame,
            symmetry=crystal_symmetry,
            phase=phase,
            provenance=provenance if query_set.provenance is None else query_set.provenance,
        )
    return query_set, scalar_output


@dataclass(frozen=True, slots=True)
class HarmonicODF:
    coefficients: np.ndarray
    basis_terms: tuple[HarmonicBasisTerm, ...]
    basis_transform: np.ndarray
    quadrature_orientations: OrientationSet
    quadrature_weights: np.ndarray
    quadrature_basis_values: np.ndarray
    degree_bandlimit: int
    crystal_symmetry: SymmetrySpec | None = None
    specimen_symmetry: SymmetrySpec | None = None
    phase: Phase | None = None
    pole_kernel: KernelSpec = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5)
    even_degrees_only: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        coefficients = as_float_array(self.coefficients, shape=(None,))
        basis_transform = as_float_array(self.basis_transform, shape=(len(self.basis_terms), None))
        quadrature_weights = as_float_array(
            self.quadrature_weights,
            shape=(len(self.quadrature_orientations),),
        )
        quadrature_basis_values = as_float_array(
            self.quadrature_basis_values,
            shape=(len(self.quadrature_orientations), coefficients.shape[0]),
        )
        if coefficients.shape[0] != basis_transform.shape[1]:
            raise ValueError("coefficients length must match the retained harmonic basis size.")
        if self.degree_bandlimit < 0:
            raise ValueError("degree_bandlimit must be non-negative.")
        if self.crystal_symmetry is not None:
            if self.crystal_symmetry.reference_frame != self.quadrature_orientations.crystal_frame:
                raise ValueError(
                    "crystal_symmetry.reference_frame must match quadrature_orientations.crystal_frame."
                )
        if self.specimen_symmetry is not None:
            if self.specimen_symmetry.reference_frame != self.quadrature_orientations.specimen_frame:
                raise ValueError(
                    "specimen_symmetry.reference_frame must match quadrature_orientations.specimen_frame."
                )
        if not np.isclose(float(np.sum(quadrature_weights)), 1.0, atol=1e-8):
            quadrature_weights = _normalized_weights(quadrature_weights)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "basis_transform", basis_transform)
        object.__setattr__(self, "quadrature_weights", quadrature_weights)
        object.__setattr__(self, "quadrature_basis_values", quadrature_basis_values)

    @property
    def crystal_frame(self) -> ReferenceFrame:
        return self.quadrature_orientations.crystal_frame

    @property
    def specimen_frame(self) -> ReferenceFrame:
        return self.quadrature_orientations.specimen_frame

    @property
    def basis_size(self) -> int:
        return int(self.coefficients.shape[0])

    @property
    def raw_basis_size(self) -> int:
        return int(len(self.basis_terms))

    @property
    def quadrature_size(self) -> int:
        return len(self.quadrature_orientations)

    @property
    def quadrature_densities(self) -> np.ndarray:
        densities = self.quadrature_basis_values @ self.coefficients
        densities = np.ascontiguousarray(densities, dtype=np.float64)
        densities.setflags(write=False)
        return densities

    @property
    def mean_density(self) -> float:
        return _weighted_mean(self.quadrature_densities, self.quadrature_weights)

    def evaluate(self, orientations: Orientation | OrientationSet) -> np.ndarray | float:
        query_set, scalar_output = _coerce_query_orientations(
            orientations,
            crystal_frame=self.crystal_frame,
            specimen_frame=self.specimen_frame,
            crystal_symmetry=self.crystal_symmetry,
            phase=self.phase,
            provenance=self.provenance,
        )
        raw_basis = _symmetry_projected_raw_basis(
            query_set,
            terms=self.basis_terms,
            crystal_symmetry=self.crystal_symmetry,
            specimen_symmetry=self.specimen_symmetry,
        )
        orthonormal_basis = raw_basis @ self.basis_transform
        values = orthonormal_basis @ self.coefficients
        values = np.ascontiguousarray(values, dtype=np.float64)
        values.setflags(write=False)
        if scalar_output:
            return float(values[0])
        return values

    def evaluate_pole_density(
        self,
        pole: CrystalPlane,
        sample_directions: ArrayLike,
        *,
        include_symmetry_family: bool = True,
    ) -> np.ndarray:
        if self.phase is not None and pole.phase != self.phase:
            raise ValueError("HarmonicODF pole evaluation requires the same phase as the ODF.")
        response = _pole_density_response_matrix(
            self.quadrature_orientations,
            pole=pole,
            sample_directions=sample_directions,
            kernel=self.pole_kernel,
            include_symmetry_family=include_symmetry_family,
        )
        weighted_density = self.quadrature_weights * self.quadrature_densities
        density = response @ weighted_density
        density = np.ascontiguousarray(density, dtype=np.float64)
        density.setflags(write=False)
        return density

    def reconstruct_pole_figure(
        self,
        pole: CrystalPlane,
        *,
        sample_directions: ArrayLike,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigure:
        return PoleFigure(
            pole=pole,
            sample_directions=np.asarray(sample_directions, dtype=np.float64),
            intensities=self.evaluate_pole_density(
                pole,
                sample_directions,
                include_symmetry_family=include_symmetry_family,
            ),
            specimen_frame=self.specimen_frame,
            antipodal=antipodal,
            sample_symmetry=self.specimen_symmetry,
            provenance=self.provenance if provenance is None else provenance,
        )

    @classmethod
    def invert_pole_figures(
        cls,
        pole_figures: Sequence[PoleFigure],
        *,
        degree_bandlimit: int,
        regularization: float = 1e-6,
        include_symmetry_family: bool = True,
        even_degrees_only: bool | None = None,
        specimen_symmetry: SymmetrySpec | None = None,
        pole_kernel: KernelSpec | None = None,
        phi1_step_deg: float = 30.0,
        big_phi_step_deg: float = 30.0,
        phi2_step_deg: float = 30.0,
        basis_tolerance: float = 1e-10,
        provenance: ProvenanceRecord | None = None,
    ) -> HarmonicODFReconstructionReport:
        if not pole_figures:
            raise ValueError("Harmonic ODF inversion requires at least one PoleFigure.")
        if regularization < 0.0:
            raise ValueError("regularization must be non-negative.")
        first = pole_figures[0]
        specimen_frame = first.specimen_frame
        phase = first.pole.phase
        crystal_frame = phase.crystal_frame
        crystal_symmetry = phase.symmetry
        for pole_figure in pole_figures[1:]:
            if pole_figure.specimen_frame != specimen_frame:
                raise ValueError("All pole figures must share the same specimen frame.")
            if pole_figure.pole.phase != phase:
                raise ValueError("All pole figures must reference the same phase.")
        if specimen_symmetry is None:
            common_sample_symmetry = first.sample_symmetry
            if any(pole_figure.sample_symmetry != common_sample_symmetry for pole_figure in pole_figures):
                common_sample_symmetry = None
        else:
            if specimen_symmetry.reference_frame != specimen_frame:
                raise ValueError(
                    "specimen_symmetry.reference_frame must match the pole figure specimen frame."
                )
            common_sample_symmetry = specimen_symmetry
        even_only = (
            all(pole_figure.antipodal for pole_figure in pole_figures)
            if even_degrees_only is None
            else even_degrees_only
        )
        inversion_kernel = (
            KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5)
            if pole_kernel is None
            else pole_kernel
        )
        quadrature_orientations, quadrature_weights = _bunge_quadrature(
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            crystal_symmetry=crystal_symmetry,
            phase=phase,
            phi1_step_deg=phi1_step_deg,
            big_phi_step_deg=big_phi_step_deg,
            phi2_step_deg=phi2_step_deg,
            provenance=provenance,
        )
        basis_terms = _enumerate_terms(
            degree_bandlimit=degree_bandlimit,
            even_degrees_only=even_only,
        )
        raw_basis = _symmetry_projected_raw_basis(
            quadrature_orientations,
            terms=basis_terms,
            crystal_symmetry=crystal_symmetry,
            specimen_symmetry=common_sample_symmetry,
        )
        quadrature_basis_values, basis_transform = _orthonormalize_weighted_basis(
            raw_basis,
            quadrature_weights,
            tolerance=basis_tolerance,
        )
        blocks = []
        for pole_figure in pole_figures:
            response = _pole_density_response_matrix(
                quadrature_orientations,
                pole=pole_figure.pole,
                sample_directions=pole_figure.sample_directions,
                kernel=inversion_kernel,
                include_symmetry_family=include_symmetry_family,
            )
            blocks.append(response @ (quadrature_weights[:, None] * quadrature_basis_values))
        system_matrix = np.vstack(blocks)
        observations = np.concatenate([pole_figure.intensities for pole_figure in pole_figures])
        if regularization > 0.0:
            augmented_matrix = np.vstack(
                [
                    system_matrix,
                    np.sqrt(regularization) * np.eye(system_matrix.shape[1], dtype=np.float64),
                ]
            )
            augmented_observations = np.concatenate(
                [observations, np.zeros(system_matrix.shape[1], dtype=np.float64)]
            )
        else:
            augmented_matrix = system_matrix
            augmented_observations = observations
        coefficients, _, rank, singular_values = np.linalg.lstsq(
            augmented_matrix,
            augmented_observations,
            rcond=None,
        )
        coefficients = np.ascontiguousarray(coefficients, dtype=np.float64)
        coefficients.setflags(write=False)
        harmonic_odf = cls(
            coefficients=coefficients,
            basis_terms=basis_terms,
            basis_transform=basis_transform,
            quadrature_orientations=quadrature_orientations,
            quadrature_weights=quadrature_weights,
            quadrature_basis_values=quadrature_basis_values,
            degree_bandlimit=degree_bandlimit,
            crystal_symmetry=crystal_symmetry,
            specimen_symmetry=common_sample_symmetry,
            phase=phase,
            pole_kernel=inversion_kernel,
            even_degrees_only=even_only,
            provenance=provenance,
        )
        predicted = np.ascontiguousarray(system_matrix @ coefficients, dtype=np.float64)
        predicted.setflags(write=False)
        residual = predicted - observations
        residual_norm = float(np.linalg.norm(residual))
        observation_norm = max(float(np.linalg.norm(observations)), 1e-12)
        condition_number = (
            float(singular_values[0] / singular_values[-1])
            if singular_values.size > 0 and singular_values[-1] > 0.0
            else float("inf")
        )
        return HarmonicODFReconstructionReport(
            odf=harmonic_odf,
            residual_norm=residual_norm,
            relative_residual_norm=float(residual_norm / observation_norm),
            mean_absolute_error=float(np.mean(np.abs(residual))),
            max_absolute_error=float(np.max(np.abs(residual))),
            regularization=regularization,
            observation_count=int(observations.size),
            basis_size=harmonic_odf.basis_size,
            raw_basis_size=harmonic_odf.raw_basis_size,
            quadrature_size=harmonic_odf.quadrature_size,
            degree_bandlimit=degree_bandlimit,
            even_degrees_only=even_only,
            matrix_rank=int(rank),
            condition_number=condition_number,
            predicted_intensities=predicted,
            mean_density=harmonic_odf.mean_density,
            minimum_density=float(np.min(harmonic_odf.quadrature_densities)),
            maximum_density=float(np.max(harmonic_odf.quadrature_densities)),
            crystal_symmetry_order=1 if crystal_symmetry is None else crystal_symmetry.order,
            specimen_symmetry_order=(
                1 if common_sample_symmetry is None else common_sample_symmetry.order
            ),
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class HarmonicODFReconstructionReport:
    odf: HarmonicODF
    residual_norm: float
    relative_residual_norm: float
    mean_absolute_error: float
    max_absolute_error: float
    regularization: float
    observation_count: int
    basis_size: int
    raw_basis_size: int
    quadrature_size: int
    degree_bandlimit: int
    even_degrees_only: bool
    matrix_rank: int
    condition_number: float
    predicted_intensities: np.ndarray
    mean_density: float
    minimum_density: float
    maximum_density: float
    crystal_symmetry_order: int
    specimen_symmetry_order: int
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        predicted = as_float_array(self.predicted_intensities, shape=(self.observation_count,))
        if self.residual_norm < 0.0:
            raise ValueError("residual_norm must be non-negative.")
        if self.relative_residual_norm < 0.0:
            raise ValueError("relative_residual_norm must be non-negative.")
        if self.mean_absolute_error < 0.0:
            raise ValueError("mean_absolute_error must be non-negative.")
        if self.max_absolute_error < 0.0:
            raise ValueError("max_absolute_error must be non-negative.")
        if self.regularization < 0.0:
            raise ValueError("regularization must be non-negative.")
        if self.observation_count <= 0:
            raise ValueError("observation_count must be strictly positive.")
        if self.basis_size <= 0:
            raise ValueError("basis_size must be strictly positive.")
        if self.raw_basis_size < self.basis_size:
            raise ValueError("raw_basis_size must be at least as large as basis_size.")
        if self.quadrature_size <= 0:
            raise ValueError("quadrature_size must be strictly positive.")
        if self.degree_bandlimit < 0:
            raise ValueError("degree_bandlimit must be non-negative.")
        if self.matrix_rank <= 0:
            raise ValueError("matrix_rank must be strictly positive.")
        if self.crystal_symmetry_order <= 0:
            raise ValueError("crystal_symmetry_order must be strictly positive.")
        if self.specimen_symmetry_order <= 0:
            raise ValueError("specimen_symmetry_order must be strictly positive.")
        object.__setattr__(self, "predicted_intensities", predicted)


__all__ = [
    "HarmonicBasisTerm",
    "HarmonicODF",
    "HarmonicODFReconstructionReport",
]
