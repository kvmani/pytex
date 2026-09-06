from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from math import lgamma
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, Phase
from pytex.core.orientation import Orientation, OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.models import (
    KernelSpec,
    PoleFigure,
    _pole_density_response_matrix,
    random_pole_density,
)

if TYPE_CHECKING:
    from pytex.texture.ghosts import GhostCorrectionReport, GhostCorrectionSpec


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
    count = round(span / step_deg)
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
    """One term of the generalized spherical harmonic basis on SO(3).

    Attributes
    ----------
    degree : int
        Harmonic degree ``l``; non-negative.
    sample_order : int
        Specimen-side order ``m``; magnitude must not exceed ``degree``.
    crystal_order : int
        Crystal-side order ``n``; magnitude must not exceed ``degree``.
    component : str
        ``"real"`` or ``"imag"``, since the real-valued basis splits each
        complex term in two.
    """

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


def _enumerate_terms(
    *,
    degree_bandlimit: int,
    even_degrees_only: bool,
) -> tuple[HarmonicBasisTerm, ...]:
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


def _log_factorial(value: int) -> float:
    """``log(value!)`` for a non-negative integer, without forming ``value!``."""

    return float(lgamma(value + 1))


def _wigner_small_d(
    degree: int,
    sample_order: int,
    crystal_order: int,
    beta_rad: np.ndarray,
) -> np.ndarray:
    # The Wigner coefficient is a ratio of factorials, and both halves overflow
    # fast: at degree 7 the numerator already exceeds int64, so evaluating it
    # directly produced a Python big integer that NumPy could only hold as an
    # object and refused to take the square root of. Working in log-gamma keeps
    # every intermediate a float and stays exact to rounding at any degree a
    # texture analysis would use.
    log_prefactor = 0.5 * (
        _log_factorial(degree + sample_order)
        + _log_factorial(degree - sample_order)
        + _log_factorial(degree + crystal_order)
        + _log_factorial(degree - crystal_order)
    )
    k_min = max(0, crystal_order - sample_order)
    k_max = min(degree - sample_order, degree + crystal_order)
    cos_half = np.cos(beta_rad / 2.0)
    sin_half = np.sin(beta_rad / 2.0)
    values = np.zeros_like(beta_rad, dtype=np.float64)
    for k in range(k_min, k_max + 1):
        log_denominator = (
            _log_factorial(degree + crystal_order - k)
            + _log_factorial(k)
            + _log_factorial(sample_order - crystal_order + k)
            + _log_factorial(degree - sample_order - k)
        )
        exponent_cos = 2 * degree + crystal_order - sample_order - 2 * k
        exponent_sin = sample_order - crystal_order + 2 * k
        coefficient = ((-1) ** (k - sample_order + crystal_order)) * np.exp(
            log_prefactor - log_denominator
        )
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
    basis = (
        np.column_stack(columns)
        if columns
        else np.zeros((angles_rad.shape[0], 0), dtype=np.float64)
    )
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


def _column_sign_convention(columns: np.ndarray) -> np.ndarray:
    """Signs that make each column's largest-magnitude entry positive.

    An orthonormal basis is only defined up to the sign of each column, and
    ``numpy.linalg.eigh`` does not pin it: the same Gram matrix yields ``+v`` on
    one LAPACK build and ``-v`` on another. That is not cosmetic here. The
    columns are what the harmonic coefficients multiply, so a flipped constant
    column turns a uniform ODF of density ``+1`` into one of density ``-1``, and
    it flips the sign of every coefficient a stored `basis_transform` reproduces.
    A uniform ODF evaluating to ``-1`` m.r.d. on Linux and ``+1`` on Windows is
    exactly the failure this prevents.

    The convention is the usual one (as in an SVD sign fix): make the entry of
    largest magnitude positive, lowest index winning a tie. It is arbitrary, but
    it is *definite*, and every platform computes the same one.
    """

    dominant = np.argmax(np.abs(columns), axis=0)
    signs = np.sign(columns[dominant, np.arange(columns.shape[1])])
    signs[signs == 0.0] = 1.0
    return np.asarray(signs, dtype=np.float64)


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
            "The symmetry-projected harmonic basis is numerically rank-deficient "
            "at the requested bandlimit."
        )
    kept_values = eigenvalues[keep]
    kept_vectors = eigenvectors[:, keep]
    transform = kept_vectors / np.sqrt(kept_values)[None, :]
    orthonormal_basis = raw_basis @ transform
    signs = _column_sign_convention(orthonormal_basis)
    orthonormal_basis = orthonormal_basis * signs[None, :]
    transform = transform * signs[None, :]
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
    """An orientation distribution represented as a harmonic series on SO(3).

    Purpose
    -------
    The series-expansion representation of the classical Bunge method. Where
    :class:`~pytex.texture.ODF` places weighted kernels on discrete support
    orientations, this stores coefficients of symmetry-projected generalized
    spherical harmonics — compact for symmetric materials, and the natural
    form for convolution and for texture-index and entropy integrals.

    Limits
    ------
    A truncated series is not sign-constrained, so sharply textured or
    under-resolved expansions can evaluate slightly negative. Pole-figure
    data determine only the even-order coefficients under Friedel's law,
    which is the classical ghost problem.

    Attributes
    ----------
    coefficients : np.ndarray
        The symmetry-projected expansion coefficients.
    basis_terms : tuple of HarmonicBasisTerm
        The raw harmonic terms before symmetry projection.
    basis_transform : np.ndarray
        The projection from raw terms onto the symmetry-allowed basis, which
        is what makes the representation compact.
    quadrature_orientations : OrientationSet
    quadrature_weights : np.ndarray
    quadrature_basis_values : np.ndarray
        The integration support and precomputed basis values, so integral
        quantities need no re-evaluation.
    degree_bandlimit : int
        Truncation degree; bounds the angular detail representable.
    crystal_symmetry, specimen_symmetry : SymmetrySpec, optional
    phase : Phase, optional
    pole_kernel : KernelSpec
        Smoothing kernel used for pole-density evaluation.
    even_degrees_only : bool
        Whether odd degrees were excluded — the honest default for
        pole-figure-derived ODFs, which cannot determine them.
    provenance : ProvenanceRecord, optional
    """

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
    pole_kernel: KernelSpec = field(
        default_factory=lambda: KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5)
    )
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
            if (
                self.crystal_symmetry.reference_frame
                != self.quadrature_orientations.crystal_frame
            ):
                raise ValueError(
                    "crystal_symmetry.reference_frame must match "
                    "quadrature_orientations.crystal_frame."
                )
        if self.specimen_symmetry is not None:
            if (
                self.specimen_symmetry.reference_frame
                != self.quadrature_orientations.specimen_frame
            ):
                raise ValueError(
                    "specimen_symmetry.reference_frame must match "
                    "quadrature_orientations.specimen_frame."
                )
        if not np.isclose(float(np.sum(quadrature_weights)), 1.0, atol=1e-8):
            quadrature_weights = _normalized_weights(quadrature_weights)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "basis_transform", basis_transform)
        object.__setattr__(self, "quadrature_weights", quadrature_weights)
        object.__setattr__(self, "quadrature_basis_values", quadrature_basis_values)

    @property
    def crystal_frame(self) -> ReferenceFrame:
        """Crystal-domain frame of the quadrature support.
        """

        return self.quadrature_orientations.crystal_frame

    @property
    def specimen_frame(self) -> ReferenceFrame:
        """Specimen-domain frame of the quadrature support.
        """

        return self.quadrature_orientations.specimen_frame

    @property
    def basis_size(self) -> int:
        """Number of symmetry-projected basis functions actually carried.

        Smaller than :attr:`raw_basis_size`: crystal and specimen symmetry
        eliminate most raw harmonic terms, which is exactly why the harmonic
        representation is compact for symmetric materials.
        """

        return int(self.coefficients.shape[0])

    @property
    def raw_basis_size(self) -> int:
        """Number of raw harmonic terms before symmetry projection.
        """

        return len(self.basis_terms)

    @property
    def quadrature_size(self) -> int:
        """Number of quadrature orientations used for integration.
        """

        return len(self.quadrature_orientations)

    @property
    def quadrature_densities(self) -> np.ndarray:
        """ODF density evaluated at every quadrature orientation, read-only.

        The array all integral quantities — mean density, texture index,
        entropy — are computed from, using the stored quadrature weights.
        """

        densities = self.quadrature_basis_values @ self.coefficients
        densities = np.ascontiguousarray(densities, dtype=np.float64)
        densities.setflags(write=False)
        return densities

    @property
    def mean_density(self) -> float:
        """Quadrature-weighted mean density over SO(3).

        Should be close to 1 for a correctly normalized ODF, since a uniform
        distribution has density 1 in multiples of random. A departure indicates
        a normalization or quadrature problem rather than a texture feature.
        """

        return _weighted_mean(self.quadrature_densities, self.quadrature_weights)

    @property
    def texture_index(self) -> float:
        """The texture index ``J = integral f(g)^2 dg`` (1 for a uniform ODF)."""

        densities = self.quadrature_densities
        return _weighted_mean(densities * densities, self.quadrature_weights)

    def entropy(self, *, floor: float = 1e-12) -> float:
        """Texture entropy ``integral f ln f dg`` (0 for a uniform ODF, positive otherwise).

        Non-positive quadrature densities (harmonic ghost lobes) are floored to
        ``floor`` before taking the logarithm so the integral stays finite.
        """

        densities = np.clip(self.quadrature_densities, floor, None)
        integrand = densities * np.log(densities)
        return _weighted_mean(integrand, self.quadrature_weights)

    def evaluate(self, orientations: Orientation | OrientationSet) -> np.ndarray | float:
        """Evaluate the ODF density at given orientations.

        Parameters
        ----------
        orientations : Orientation or OrientationSet
            Query orientations; must share the ODF's crystal and specimen
            frames. A single orientation returns a scalar.

        Returns
        -------
        float or np.ndarray
            Density in multiples of a random distribution. Because a truncated
            harmonic series is not sign-constrained, sharply textured or
            under-resolved expansions can produce small negative values; treat
            them as truncation artefacts, not as densities.
        """

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
        antipodal: bool = True,
    ) -> np.ndarray:
        """Pole density at chosen specimen directions.

        Parameters
        ----------
        pole : CrystalPlane
            Must match the ODF's phase when one is declared.
        sample_directions : ArrayLike
            ``(n, 3)`` specimen directions.
        include_symmetry_family : bool
            Include the whole ``{hkl}`` family (default), as a measurement does.
        antipodal : bool
            Identify opposite normals (default), which is Friedel's law and what
            a diffraction pole figure measures. Switching it off models a
            hypothetical experiment that can tell ``h`` from ``-h``, in which
            the odd-degree part of the ODF becomes visible and the ghost problem
            does not arise.

        Returns
        -------
        np.ndarray
            ``(n,)`` pole densities in multiples of a random distribution, on the
            same scale as the orientation density itself: a uniform ODF
            (:attr:`mean_density` 1) returns 1.0 at every direction.

        Notes
        -----
        The kernel used to smooth the quadrature sum peaks at 1 rather than
        integrating to 1, so the raw response of a uniform distribution is the
        kernel's spherical mean — about 0.006 at the default 7.5 degree
        halfwidth, not 1. Dividing by
        :func:`~pytex.texture.models.random_pole_density` removes that factor,
        which is the only reason this returns m.r.d. rather than a response.
        """

        if self.phase is not None and pole.phase != self.phase:
            raise ValueError("HarmonicODF pole evaluation requires the same phase as the ODF.")
        response = _pole_density_response_matrix(
            self.quadrature_orientations,
            pole=pole,
            sample_directions=sample_directions,
            kernel=self.pole_kernel,
            include_symmetry_family=include_symmetry_family,
            antipodal=antipodal,
        )
        weighted_density = self.quadrature_weights * self.quadrature_densities
        density = (response @ weighted_density) / random_pole_density(
            self.pole_kernel, antipodal=antipodal
        )
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
        """The pole figure this harmonic ODF predicts for a crystal plane.

        Unlike the discrete ODF's reconstruction, the specimen directions must be
        supplied: a harmonic ODF has no scattered support of its own to reuse, so
        the evaluation grid is the caller's choice.

        Parameters
        ----------
        pole : CrystalPlane
        sample_directions : ArrayLike
            ``(n, 3)`` specimen directions to evaluate on.
        include_symmetry_family : bool
        antipodal : bool
        provenance : ProvenanceRecord, optional
        """

        return PoleFigure(
            pole=pole,
            sample_directions=np.asarray(sample_directions, dtype=np.float64),
            intensities=self.evaluate_pole_density(
                pole,
                sample_directions,
                include_symmetry_family=include_symmetry_family,
                antipodal=antipodal,
            ),
            specimen_frame=self.specimen_frame,
            antipodal=antipodal,
            includes_symmetry_family=include_symmetry_family,
            sample_symmetry=self.specimen_symmetry,
            # The intensities are densities evaluated at the supplied
            # directions, not weights of individual poles.
            sampling="sampled_density",
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
        ghost_correction: GhostCorrectionSpec | bool | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> HarmonicODFReconstructionReport:
        """Estimate a harmonic ODF from measured pole figures.

        Purpose
        -------
        The series-expansion route to PF-to-ODF inversion — the classical Bunge
        method — in which the unknown is a truncated set of symmetry-projected
        harmonic coefficients rather than weights on a discrete support.

        Method and limits
        -----------------
        Builds the pole-density response of each basis function at every measured
        pole-figure point and solves the regularized least-squares system for the
        coefficients. Two limits are intrinsic and are not worked around here:
        pole figures determine only the even-order coefficients under Friedel's
        law, so the odd part is unconstrained — the classical ghost problem — and
        truncation at ``degree_bandlimit`` bounds the angular detail recoverable.
        The truncation is intrinsic. The ghost problem is not left standing: pass
        ``ghost_correction`` to recover an odd part from positivity, and read
        :attr:`HarmonicODFReconstructionReport.ghost_correction` for what that
        inference cost. Without it the result is an even-only, band-limited,
        regularized estimate whose odd part was silently set to zero.

        Parameters
        ----------
        pole_figures : sequence of PoleFigure
            All must share a specimen frame. More independent poles constrain the
            problem better.
        degree_bandlimit : int
            Highest harmonic degree retained. Cost and achievable detail both
            rise with it; so does sensitivity to noise.
        regularization : float
            Tikhonov weight; larger is smoother and more stable.
        include_symmetry_family : bool
            Model the whole ``{hkl}`` family per pole figure (default).
        even_degrees_only : bool, optional
            Restrict to even degrees. ``None`` selects the honest default, since
            odd degrees are not determined by pole-figure data.
        specimen_symmetry : SymmetrySpec, optional
            Statistical specimen symmetry to project onto.
        pole_kernel : KernelSpec, optional
            Smoothing kernel for the response model.
        phi1_step_deg, big_phi_step_deg, phi2_step_deg : float
            Quadrature grid spacing in Bunge Euler space.
        ghost_correction : GhostCorrectionSpec or bool, optional
            Recover the odd part from positivity after the inversion. ``True``
            applies the default :class:`~pytex.texture.GhostCorrectionSpec`;
            a spec applies that spec; ``None`` or ``False`` leaves the odd part
            at zero. The corrected distribution is reached through
            :attr:`~HarmonicODFReconstructionReport.final_odf`, never by
            replacing :attr:`~HarmonicODFReconstructionReport.odf`, so the
            residuals and density diagnostics in the report always describe the
            solution the data alone produced.

        Returns
        -------
        HarmonicODFReconstructionReport
            The estimated ODF with the residual and conditioning information
            needed to judge the fit.

        See Also
        --------
        pytex.texture.correct_ghosts
        """

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
            if any(
                pole_figure.sample_symmetry != common_sample_symmetry
                for pole_figure in pole_figures
            ):
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
        # The forward operator must predict the observations on their own scale.
        # Measured and reconstructed pole figures are in multiples of a random
        # distribution, while the raw kernel response of a uniform ODF is the
        # kernel's spherical mean (0.006 at the default halfwidth), so an
        # unnormalized operator would be fitted by coefficients inflated by its
        # reciprocal — about 163 — and every density the returned ODF reports,
        # including mean_density, texture_index and entropy, would carry that
        # factor. The pole-figure round trip would still close, which is why this
        # scale error is invisible unless the ODF itself is inspected.
        blocks = []
        for pole_figure in pole_figures:
            # Friedel folding is a property of the measurement, so it is read
            # off each figure rather than chosen here; the normalization must
            # follow it, or the fitted densities carry a factor of two.
            kernel_mean = random_pole_density(
                inversion_kernel, antipodal=pole_figure.antipodal
            )
            response = _pole_density_response_matrix(
                quadrature_orientations,
                pole=pole_figure.pole,
                sample_directions=pole_figure.sample_directions,
                kernel=inversion_kernel,
                include_symmetry_family=include_symmetry_family,
                antipodal=pole_figure.antipodal,
            )
            blocks.append(
                (response @ (quadrature_weights[:, None] * quadrature_basis_values)) / kernel_mean
            )
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
        ghost_report: GhostCorrectionReport | None = None
        if ghost_correction:
            # Imported here rather than at module scope: the ghost module builds
            # on this one, and a top-level import would close the cycle.
            from pytex.texture.ghosts import correct_ghosts

            ghost_report = correct_ghosts(
                harmonic_odf,
                spec=None if ghost_correction is True else ghost_correction,
                pole_figures=pole_figures,
                include_symmetry_family=include_symmetry_family,
                provenance=provenance,
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
            coefficient_l2_norm=float(np.linalg.norm(coefficients)),
            coefficient_max_abs=float(np.max(np.abs(coefficients))),
            negative_density_fraction=float(np.mean(harmonic_odf.quadrature_densities < 0.0)),
            ghost_correction=ghost_report,
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class HarmonicODFReconstructionReport:
    """A harmonic PF-to-ODF reconstruction with the evidence to judge it.

    Purpose
    -------
    Pole-figure inversion is ill-posed, so the ODF alone is not a result.
    This adds the residuals against the measured data, the problem
    dimensions, and — importantly — the conditioning of the system, so a
    reader can tell a well-constrained solution from one held up entirely by
    regularization.

    Attributes
    ----------
    odf : HarmonicODF
    residual_norm, relative_residual_norm : float
    mean_absolute_error, max_absolute_error : float
    regularization : float
    observation_count : int
    basis_size, raw_basis_size : int
        Symmetry-projected and raw basis sizes; their ratio shows how much
        the symmetry bought.
    quadrature_size : int
    degree_bandlimit : int
    even_degrees_only : bool
    matrix_rank, condition_number : float
        Conditioning of the least-squares system. A rank below the basis size
        means the data do not determine every coefficient.
    predicted_intensities : np.ndarray
        Pole densities the solution implies, for direct comparison.
    mean_density : float
        Should be near 1 for a correctly normalized ODF.
    ghost_correction : GhostCorrectionReport, optional
        Present when a ghost correction was requested. Its ``odf`` is the
        corrected distribution and its ``describe()`` states what the correction
        assumed and what it cost; every other field of this report continues to
        describe the uncorrected inversion.
    provenance : ProvenanceRecord, optional
    """

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
    coefficient_l2_norm: float = 0.0
    coefficient_max_abs: float = 0.0
    negative_density_fraction: float = 0.0
    ghost_correction: GhostCorrectionReport | None = None
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
        if self.coefficient_l2_norm < 0.0:
            raise ValueError("coefficient_l2_norm must be non-negative.")
        if self.coefficient_max_abs < 0.0:
            raise ValueError("coefficient_max_abs must be non-negative.")
        if not 0.0 <= self.negative_density_fraction <= 1.0:
            raise ValueError("negative_density_fraction must lie in [0, 1].")
        object.__setattr__(self, "predicted_intensities", predicted)

    @property
    def final_odf(self) -> HarmonicODF:
        """The distribution this reconstruction recommends quoting numbers from.

        The ghost-corrected ODF when a correction was applied, and the direct
        inversion output otherwise. Use this — rather than :attr:`odf` — for any
        density, texture index, entropy, volume fraction or Kearns parameter, so
        that requesting a correction actually changes the numbers reported
        downstream instead of only adding a report nobody reads.
        """

        if self.ghost_correction is None:
            return self.odf
        return self.ghost_correction.odf


__all__ = [
    "HarmonicBasisTerm",
    "HarmonicODF",
    "HarmonicODFReconstructionReport",
]
