"""Ghost correction for pole-figure-derived harmonic ODFs.

A pole figure measured under Friedel's law determines only the even-degree part
of the orientation distribution. The odd part is not merely poorly determined,
it is not determined *at all* by any amount of pole-figure data, and the
even-only solution that results carries the classical *ghost* artefacts: false
maxima where the true distribution is empty, and depressed true maxima.

The correction implemented here is the positivity / zero-range family
(Dahms and Bunge). Classically it is stated as an alternating projection between
two convex sets in the quadrature-weighted function space:

* :math:`C_+` — the densities that are physically admissible: non-negative
  everywhere, and identically zero on a declared zero range;
* :math:`C_d` — the densities whose even part equals the measured even part,
  which is an affine set because the even part is exactly what the data fixed.

That characterisation is what makes the answer meaningful — the corrected ODF is
an admissible distribution near the even-only solution, and it reproduces the
measured pole figures exactly as well, because nothing touches the even part.
The *computation* is the equivalent convex minimization rather than the
iteration: the violation of admissibility is a smooth convex function of the
odd coefficients, so a quasi-Newton minimizer reaches the same point in tens of
steps where alternating projection takes thousands. See
``docs/site/algorithms/ghost_correction.md``.

References
----------
Dahms, M. and Bunge, H. J. (1989) "The iterative series-expansion method for
quantitative texture analysis. I. General outline", *Journal of Applied
Crystallography* 22, 439-447, :doi:`10.1107/S0021889889005261`.

Bunge, H. J. (1982) *Texture Analysis in Materials Science*, Butterworths,
chapter 13 (the ghost problem and the zero-range method).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from pytex.core.provenance import ProvenanceRecord
from pytex.texture.harmonics import (
    HarmonicBasisTerm,
    HarmonicODF,
    _enumerate_terms,
    _orthonormalize_weighted_basis,
    _symmetry_projected_raw_basis,
    _weighted_mean,
)
from pytex.texture.models import (
    PoleFigure,
    _pole_density_response_matrix,
    random_pole_density,
)

GHOST_CORRECTION_METHODS = ("positivity", "zero_range")


@dataclass(frozen=True, slots=True)
class GhostCorrectionSpec:
    """Settings for a ghost correction, and nothing that the data supply.

    Purpose
    -------
    Every field here is a *choice*, not a measurement. Keeping them in one
    declared object is what lets :meth:`GhostCorrectionReport.describe` state
    which choices produced the corrected distribution, so a reader can tell the
    part of the answer that came from the specimen from the part that came from
    the analyst.

    Attributes
    ----------
    method : str
        ``"positivity"`` enforces only :math:`f \\ge 0`. ``"zero_range"``
        additionally forces the density to zero wherever the even-only solution
        falls below :attr:`zero_range_threshold`, which is the classical
        zero-range method: an orientation the measurement says is empty stays
        empty, and the odd part is spent making it so.
    zero_range_threshold : float
        The density in multiples of a random distribution (m.r.d.) below which
        the even-only solution is read as declaring an empty range. Used by
        ``"zero_range"`` only. Must be non-negative; a threshold of exactly
        zero makes ``"zero_range"`` weaker than ``"positivity"`` is, not
        stronger, so a strictly positive value is expected.
    max_iterations : int
        Cap on the minimizer's iterations. Reaching it is reported rather than
        raised: a correction stopped early is still a valid lower bound on the
        correction the data admit, provided the reader is told.
    tolerance : float
        Gradient tolerance handed to the minimizer. The gradient has units of
        m.r.d. squared per unit coefficient, so it is a statement about how flat
        the infeasibility surface has become, not about the density directly;
        :attr:`GhostCorrectionReport.infeasibility_after` is the number to read
        for what the residual violation actually is.
    odd_regularization : float
        Weight on the squared norm of the odd coefficients. Positivity alone
        does not determine the odd part: once the density is admissible
        anywhere in the feasible set, every remaining direction is free, and an
        unregularized minimizer stops at whichever feasible point it reaches
        first — a larger odd part than the data force, presented as though the
        data forced it. This term selects the *smallest* odd part that achieves
        admissibility, which is the only defensible choice when the alternative
        is an arbitrary one. It must be small enough not to bias the
        admissibility itself; the default is six orders below a typical
        violation.
    degree_bandlimit : int, optional
        Bandlimit for the odd expansion. ``None`` uses the ODF's own bandlimit,
        which is the honest default: an odd part resolved more finely than the
        even part it corrects would put detail into the answer that no data
        constrain.
    basis_tolerance : float
        Eigenvalue floor when orthonormalizing the odd basis.
    """

    method: str = "positivity"
    zero_range_threshold: float = 0.05
    max_iterations: int = 500
    tolerance: float = 1e-12
    odd_regularization: float = 1e-6
    degree_bandlimit: int | None = None
    basis_tolerance: float = 1e-10

    def __post_init__(self) -> None:
        if self.method not in GHOST_CORRECTION_METHODS:
            raise ValueError(
                "GhostCorrectionSpec.method must be one of "
                f"{GHOST_CORRECTION_METHODS!r}; got {self.method!r}."
            )
        if self.zero_range_threshold < 0.0:
            raise ValueError("zero_range_threshold must be non-negative.")
        if self.max_iterations <= 0:
            raise ValueError("max_iterations must be strictly positive.")
        if self.tolerance <= 0.0:
            raise ValueError("tolerance must be strictly positive.")
        if self.odd_regularization < 0.0:
            raise ValueError("odd_regularization must be non-negative.")
        if self.degree_bandlimit is not None and self.degree_bandlimit < 0:
            raise ValueError("degree_bandlimit must be non-negative.")
        if self.basis_tolerance <= 0.0:
            raise ValueError("basis_tolerance must be strictly positive.")


@dataclass(frozen=True, slots=True)
class GhostCorrectionReport:
    """A ghost-corrected ODF together with what the correction cost.

    Purpose
    -------
    The corrected ODF alone would be indistinguishable from a measured one, and
    it is not: its odd part is an inference from positivity, not an observation.
    This report carries the size of that inference (:attr:`odd_coefficient_norm`
    against :attr:`even_coefficient_norm`), the artefact it removed
    (:attr:`negative_density_fraction_before` against ``_after``), the change it
    made to every derived texture quantity, and the evidence that it did not
    damage the fit to the measured pole figures
    (:attr:`pole_figure_max_change`).

    Attributes
    ----------
    odf : HarmonicODF
        The corrected distribution, carrying both parities.
    method : str
        The method from the spec that produced it.
    iterations, converged : int, bool
        Minimizer iteration count and whether it reported success.
    infeasibility_before, infeasibility_after : float
        Quadrature-weighted L2 norm, in m.r.d., of the part of the density that
        violates the constraints — the negative part, plus the whole density
        inside a declared zero range. This is the quantity the correction
        minimizes, and the gap it cannot close is the honest statement that a
        band-limited expansion of this data cannot be made strictly admissible.
    odd_basis_size, even_basis_size : int
        Sizes of the two symmetry-projected bases. An ``odd_basis_size`` of zero
        means the symmetry admits no odd-degree terms at this bandlimit at all,
        so no ghost correction is possible — or needed.
    zero_range_fraction : float
        Fraction of the quadrature declared empty, for ``"zero_range"``.
    negative_density_fraction_before, negative_density_fraction_after : float
    minimum_density_before, minimum_density_after : float
    maximum_density_before, maximum_density_after : float
    mean_density_before, mean_density_after : float
        The mean must not move: odd-degree harmonics integrate to zero over
        SO(3), so a moved mean is a defect, not a correction.
    texture_index_before, texture_index_after : float
    entropy_before, entropy_after : float
    even_coefficient_norm, odd_coefficient_norm : float
    pole_figure_max_change : float, optional
        Largest change, over every measured direction, in the pole density the
        ODF predicts. Under Friedel's law this must be at the level of the
        quadrature error; a large value means the correction bought positivity
        by moving the fit, which it is not entitled to do.
    provenance : ProvenanceRecord, optional
    """

    odf: HarmonicODF
    method: str
    iterations: int
    converged: bool
    infeasibility_before: float
    infeasibility_after: float
    odd_basis_size: int
    even_basis_size: int
    zero_range_fraction: float
    negative_density_fraction_before: float
    negative_density_fraction_after: float
    minimum_density_before: float
    minimum_density_after: float
    maximum_density_before: float
    maximum_density_after: float
    mean_density_before: float
    mean_density_after: float
    texture_index_before: float
    texture_index_after: float
    entropy_before: float
    entropy_after: float
    even_coefficient_norm: float
    odd_coefficient_norm: float
    pole_figure_max_change: float | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.method not in GHOST_CORRECTION_METHODS:
            raise ValueError(f"Unknown ghost-correction method {self.method!r}.")
        if self.iterations < 0:
            raise ValueError("iterations must be non-negative.")
        if self.infeasibility_before < 0.0 or self.infeasibility_after < 0.0:
            raise ValueError("Infeasibility norms must be non-negative.")
        if self.infeasibility_after > self.infeasibility_before + 1e-9:
            raise ValueError(
                "A ghost correction that increases the infeasibility it minimizes is a defect: "
                "the zero odd part is always available, so the minimum can never be worse."
            )
        if self.odd_basis_size < 0 or self.even_basis_size <= 0:
            raise ValueError("Basis sizes must be non-negative, with a non-empty even basis.")
        for name in (
            "zero_range_fraction",
            "negative_density_fraction_before",
            "negative_density_fraction_after",
        ):
            value = float(getattr(self, name))
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1].")
        if self.even_coefficient_norm < 0.0 or self.odd_coefficient_norm < 0.0:
            raise ValueError("Coefficient norms must be non-negative.")
        if self.pole_figure_max_change is not None and self.pole_figure_max_change < 0.0:
            raise ValueError("pole_figure_max_change must be non-negative.")

    @property
    def ghost_amplitude_ratio(self) -> float:
        """Odd coefficient norm as a fraction of the even one.

        The one number that says how much of the reported distribution was
        inferred rather than measured. Small values mean the even-only solution
        was already nearly admissible; values approaching one mean the reported
        ODF is half an inference.
        """

        if self.even_coefficient_norm <= 0.0:
            return 0.0
        return float(self.odd_coefficient_norm / self.even_coefficient_norm)

    def describe(self) -> str:
        """Prose stating what was corrected, by how much, and at what cost."""

        if self.odd_basis_size == 0:
            return (
                "No ghost correction was applied: the crystal and specimen symmetries admit no "
                f"odd-degree harmonic terms at bandlimit {self.odf.degree_bandlimit}, so the "
                "even-only solution already spans every function the symmetry allows. The "
                "distribution is unchanged, and it has no ghost part to remove."
            )
        method_prose = (
            "positivity alone"
            if self.method == "positivity"
            else "the zero-range method, which additionally holds the density at zero over the "
            f"{self.zero_range_fraction:.1%} of orientation space where the even-only solution "
            "fell below the declared empty-range threshold"
        )
        convergence = (
            f"converged in {self.iterations} iterations"
            if self.converged
            else f"stopped after {self.iterations} iterations without converging, so the "
            "correction reported here is a lower bound on the one the data admit"
        )
        fit = (
            ""
            if self.pole_figure_max_change is None
            else (
                "The pole densities the ODF predicts at the measured directions changed by at "
                f"most {self.pole_figure_max_change:.3e} m.r.d., which is the check that matters: "
                "odd-degree terms are invisible to a Friedel-symmetric pole figure, so a "
                "correction that changed the fit would have been buying positivity with data "
                "agreement it is not entitled to spend. "
            )
        )
        return (
            f"Ghost correction by {method_prose}, {convergence}. The even part determined by the "
            f"pole figures was held fixed and an odd part was added from a "
            f"{self.odd_basis_size}-function symmetry-projected basis, against "
            f"{self.even_basis_size} even functions; its coefficient norm is "
            f"{self.odd_coefficient_norm:.4f} against {self.even_coefficient_norm:.4f} even, a "
            f"ghost amplitude ratio of {self.ghost_amplitude_ratio:.3f}. Negative density fell "
            f"from {self.negative_density_fraction_before:.1%} of orientation space to "
            f"{self.negative_density_fraction_after:.1%}, and the minimum density from "
            f"{self.minimum_density_before:.4f} to {self.minimum_density_after:.4f} m.r.d. The "
            f"maximum rose from {self.maximum_density_before:.4f} to "
            f"{self.maximum_density_after:.4f} m.r.d. — the classical signature of ghost removal, "
            "since the even-only solution depresses true maxima to pay for the false ones. The "
            f"texture index moved from {self.texture_index_before:.4f} to "
            f"{self.texture_index_after:.4f} and the entropy from {self.entropy_before:.4f} to "
            f"{self.entropy_after:.4f}. The mean density is {self.mean_density_after:.6f} against "
            f"{self.mean_density_before:.6f} before, and must not have moved: odd-degree "
            f"harmonics integrate to zero over SO(3). The constraint violation it minimizes — "
            f"the quadrature-weighted norm of the inadmissible part of the density — fell from "
            f"{self.infeasibility_before:.4e} to {self.infeasibility_after:.4e} m.r.d. {fit}"
            "The odd part is an inference from the assumption that a density cannot be negative, "
            "not a measurement; no pole-figure experiment can confirm or refute it."
        )


def _odd_terms(degree_bandlimit: int) -> tuple[HarmonicBasisTerm, ...]:
    """The odd-degree harmonic terms up to the bandlimit."""

    return tuple(
        term
        for term in _enumerate_terms(
            degree_bandlimit=degree_bandlimit,
            even_degrees_only=False,
        )
        if term.degree % 2 == 1
    )


def _odd_basis(
    odf: HarmonicODF,
    *,
    degree_bandlimit: int,
    basis_tolerance: float,
) -> tuple[tuple[HarmonicBasisTerm, ...], np.ndarray, np.ndarray]:
    """Orthonormal odd-degree basis on the ODF's own quadrature.

    Returns the raw terms, the quadrature values, and the transform. A symmetry
    that admits no odd terms yields empty arrays rather than an error: for cubic
    crystal symmetry the first odd-degree invariant appears only at degree 9, so
    an empty odd basis is the ordinary case at modest bandlimits, not a failure.
    """

    terms = _odd_terms(degree_bandlimit)
    quadrature_size = odf.quadrature_size
    if not terms:
        return (), np.zeros((quadrature_size, 0)), np.zeros((0, 0))
    raw_basis = _symmetry_projected_raw_basis(
        odf.quadrature_orientations,
        terms=terms,
        crystal_symmetry=odf.crystal_symmetry,
        specimen_symmetry=odf.specimen_symmetry,
    )
    try:
        values, transform = _orthonormalize_weighted_basis(
            raw_basis,
            odf.quadrature_weights,
            tolerance=basis_tolerance,
        )
    except ValueError:
        # Every odd term was annihilated by the symmetry projection, which is a
        # statement about the symmetry rather than an error in the request.
        return terms, np.zeros((quadrature_size, 0)), np.zeros((len(terms), 0))
    return terms, values, transform


def _combined_odf(
    odf: HarmonicODF,
    *,
    odd_terms: tuple[HarmonicBasisTerm, ...],
    odd_values: np.ndarray,
    odd_transform: np.ndarray,
    odd_coefficients: np.ndarray,
) -> HarmonicODF:
    """The ODF carrying both parities, so ``evaluate`` works off the quadrature."""

    even_terms = odf.basis_terms
    even_transform = odf.basis_transform
    n_even_terms, n_even = even_transform.shape
    n_odd_terms, n_odd = odd_transform.shape
    transform = np.zeros((n_even_terms + n_odd_terms, n_even + n_odd), dtype=np.float64)
    transform[:n_even_terms, :n_even] = even_transform
    transform[n_even_terms:, n_even:] = odd_transform
    return HarmonicODF(
        coefficients=np.concatenate([odf.coefficients, odd_coefficients]),
        basis_terms=even_terms + odd_terms,
        basis_transform=transform,
        quadrature_orientations=odf.quadrature_orientations,
        quadrature_weights=odf.quadrature_weights,
        quadrature_basis_values=np.hstack([odf.quadrature_basis_values, odd_values]),
        degree_bandlimit=max(
            odf.degree_bandlimit,
            max((term.degree for term in odd_terms), default=0),
        ),
        crystal_symmetry=odf.crystal_symmetry,
        specimen_symmetry=odf.specimen_symmetry,
        phase=odf.phase,
        pole_kernel=odf.pole_kernel,
        even_degrees_only=False,
        provenance=odf.provenance,
    )


def _predicted_pole_densities(
    odf: HarmonicODF,
    pole_figures: Sequence[PoleFigure],
    *,
    include_symmetry_family: bool,
) -> np.ndarray:
    """Pole densities the ODF predicts at every measured direction."""

    weighted = odf.quadrature_weights * odf.quadrature_densities
    predictions = []
    for pole_figure in pole_figures:
        kernel_mean = random_pole_density(odf.pole_kernel, antipodal=pole_figure.antipodal)
        response = _pole_density_response_matrix(
            odf.quadrature_orientations,
            pole=pole_figure.pole,
            sample_directions=pole_figure.sample_directions,
            kernel=odf.pole_kernel,
            include_symmetry_family=include_symmetry_family,
            antipodal=pole_figure.antipodal,
        )
        predictions.append((response @ weighted) / kernel_mean)
    return np.concatenate(predictions)


def correct_ghosts(
    odf: HarmonicODF,
    *,
    spec: GhostCorrectionSpec | None = None,
    pole_figures: Sequence[PoleFigure] | None = None,
    include_symmetry_family: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> GhostCorrectionReport:
    """Add the odd part a pole figure cannot measure, from positivity alone.

    Purpose
    -------
    Use this on any ODF obtained from pole figures, before quoting a maximum
    density, a texture index, an entropy, a volume fraction, or a Kearns
    parameter from it. All of those are integrals of a distribution whose odd
    part the even-only solution silently set to zero, which is not a neutral
    choice: it is the choice that produces ghosts.

    When and where
    --------------
    After :meth:`~pytex.texture.HarmonicODF.invert_pole_figures`, or through
    that method's own ``ghost_correction`` argument. It is not applicable to an
    EBSD-derived ODF, which has no ghost problem — individual orientations
    determine both parities — and not to an ODF that already carries odd
    degrees, which this function refuses rather than silently overwriting.

    Parameters
    ----------
    odf : HarmonicODF
        An even-degrees-only harmonic ODF.
    spec : GhostCorrectionSpec, optional
        The correction's settings; the default applies positivity alone.
    pole_figures : sequence of PoleFigure, optional
        The measured figures. Supplying them costs one more forward evaluation
        and buys the check that the correction left the fit alone.
    include_symmetry_family : bool
        Family expansion used when checking against ``pole_figures``; must match
        what the inversion used.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    GhostCorrectionReport
        The corrected ODF and the full cost of the correction.

    Raises
    ------
    ValueError
        If the ODF already carries odd degrees.

    See Also
    --------
    pytex.texture.HarmonicODF.invert_pole_figures
    """

    settings = GhostCorrectionSpec() if spec is None else spec
    if not odf.even_degrees_only:
        raise ValueError(
            "Ghost correction applies to an even-degrees-only ODF: the odd part is what it "
            "supplies. This ODF already carries odd degrees, and overwriting them would "
            "discard whatever determined them."
        )
    bandlimit = (
        odf.degree_bandlimit if settings.degree_bandlimit is None else settings.degree_bandlimit
    )
    odd_terms, odd_values, odd_transform = _odd_basis(
        odf,
        degree_bandlimit=bandlimit,
        basis_tolerance=settings.basis_tolerance,
    )
    weights = odf.quadrature_weights
    even_density = np.asarray(odf.quadrature_densities, dtype=np.float64)
    odd_basis_size = int(odd_values.shape[1])

    zero_mask = np.zeros(even_density.shape, dtype=bool)
    if settings.method == "zero_range":
        zero_mask = even_density < settings.zero_range_threshold

    def violation(density: np.ndarray) -> np.ndarray:
        """The inadmissible part of a density: what the correction must remove.

        Outside a zero range only a negative density is inadmissible. Inside one
        the whole density is, because the measurement's own even part declared
        that range empty and the odd part is being asked to keep it so.
        """

        negative = np.asarray(np.minimum(density, 0.0), dtype=np.float64)
        if settings.method == "zero_range":
            return np.asarray(np.where(zero_mask, density, negative), dtype=np.float64)
        return negative

    def infeasibility(coefficients: np.ndarray) -> float:
        residual = violation(even_density + odd_values @ coefficients)
        return float(np.sqrt(np.sum(weights * residual * residual)))

    def objective(coefficients: np.ndarray) -> tuple[float, np.ndarray]:
        residual = violation(even_density + odd_values @ coefficients)
        value = 0.5 * float(np.sum(weights * residual * residual))
        gradient = odd_values.T @ (weights * residual)
        penalty = settings.odd_regularization
        if penalty > 0.0:
            value += 0.5 * penalty * float(coefficients @ coefficients)
            gradient = gradient + penalty * coefficients
        return value, gradient

    start = np.zeros(odd_basis_size, dtype=np.float64)
    infeasibility_before = infeasibility(start)
    if odd_basis_size == 0:
        coefficients = start
        iterations = 0
        converged = True
        infeasibility_after = infeasibility_before
    else:
        # Smooth convex least-squares in the odd coefficients: the violation is a
        # continuously differentiable function of the density, so a quasi-Newton
        # method reaches the same point that alternating projection between the
        # positivity set and the fixed-even-part set converges to, in tens of
        # iterations instead of thousands.
        solution = minimize(
            objective,
            start,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": settings.max_iterations, "gtol": settings.tolerance},
        )
        coefficients = np.ascontiguousarray(solution.x, dtype=np.float64)
        iterations = int(solution.nit)
        converged = bool(solution.success)
        infeasibility_after = infeasibility(coefficients)

    corrected = _combined_odf(
        odf,
        odd_terms=odd_terms,
        odd_values=odd_values,
        odd_transform=odd_transform,
        odd_coefficients=coefficients,
    )
    corrected_density = corrected.quadrature_densities

    pole_figure_max_change: float | None = None
    if pole_figures:
        before = _predicted_pole_densities(
            odf, pole_figures, include_symmetry_family=include_symmetry_family
        )
        after = _predicted_pole_densities(
            corrected, pole_figures, include_symmetry_family=include_symmetry_family
        )
        pole_figure_max_change = float(np.max(np.abs(after - before)))

    return GhostCorrectionReport(
        odf=corrected,
        method=settings.method,
        iterations=int(iterations),
        converged=bool(converged),
        infeasibility_before=infeasibility_before,
        infeasibility_after=infeasibility_after,
        odd_basis_size=odd_basis_size,
        even_basis_size=odf.basis_size,
        zero_range_fraction=float(np.mean(zero_mask)),
        negative_density_fraction_before=float(np.mean(even_density < 0.0)),
        negative_density_fraction_after=float(np.mean(corrected_density < 0.0)),
        minimum_density_before=float(np.min(even_density)),
        minimum_density_after=float(np.min(corrected_density)),
        maximum_density_before=float(np.max(even_density)),
        maximum_density_after=float(np.max(corrected_density)),
        mean_density_before=_weighted_mean(even_density, weights),
        mean_density_after=_weighted_mean(corrected_density, weights),
        texture_index_before=odf.texture_index,
        texture_index_after=corrected.texture_index,
        entropy_before=odf.entropy(),
        entropy_after=corrected.entropy(),
        even_coefficient_norm=float(np.linalg.norm(odf.coefficients)),
        odd_coefficient_norm=float(np.linalg.norm(coefficients)),
        pole_figure_max_change=pole_figure_max_change,
        provenance=odf.provenance if provenance is None else provenance,
    )


__all__ = [
    "GHOST_CORRECTION_METHODS",
    "GhostCorrectionReport",
    "GhostCorrectionSpec",
    "correct_ghosts",
]
