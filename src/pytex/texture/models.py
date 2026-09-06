from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector, normalize_vectors
from pytex.core.batches import VectorSet
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import CrystalPlane, phases_semantically_match
from pytex.core.notation import format_plane_family_indices, format_plane_indices
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.core.sphere import S2Grid
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.kernels import (
    AbelPoissonKernel,
    DeLaValleePoussinKernel,
    GaussianSO3Kernel,
)
from pytex.texture.projections import project_directions

#: How a pole figure's ``intensities`` are to be read. ``"scattered_poles"`` means
#: they are per-pole weights of a cloud of individual poles, as produced by
#: mapping orientations through a plane normal; the density field they represent
#: is recovered by kernel *density estimation* (a weighted sum). ``"sampled_density"``
#: means they are already pole densities evaluated at the given directions, as a
#: diffractometer raster or an ODF evaluation produces; that field is resampled by
#: kernel *interpolation* (a weighted mean). Applying the wrong estimator to a
#: latitude-longitude raster biases the result towards the poles, where such a
#: raster oversamples, so the reading is recorded rather than assumed.
PoleFigureSampling = Literal["scattered_poles", "sampled_density"]

_POLE_FIGURE_SAMPLINGS = ("scattered_poles", "sampled_density")

#: Resampling estimators for :meth:`PoleFigure.on_grid`.
ResamplingEstimator = Literal["density", "interpolate"]

_RESAMPLING_ESTIMATORS = ("density", "interpolate")

#: Default kernel halfwidth for spherical resampling, in degrees. Five degrees is
#: the conventional texture-measurement step, so a figure resampled at the default
#: is smoothed on the scale it was measured at rather than beyond it.
DEFAULT_RESAMPLING_HALFWIDTH_DEG = 5.0


#: Gauss-Legendre nodes used to integrate a kernel over the sphere. The
#: integrand is smooth in ``cos(omega)``, so this is far more accuracy than any
#: kernel halfwidth in use requires.
_KERNEL_MEAN_NODES = 512


def random_pole_density(kernel: KernelSpec, *, antipodal: bool = False) -> float:
    r"""The pole density a *random* texture produces under this kernel.

    Purpose
    -------
    :meth:`ODF.evaluate_pole_density` returns a kernel-weighted **response**,
    not a value in multiples of a random distribution: the kernel there peaks at
    1 rather than integrating to 1, so a random texture yields the kernel's
    spherical mean instead of unity. That mean is this function, and dividing by
    it is what converts a response into m.r.d.

    When to use
    -----------
    Whenever a discrete ODF's pole densities must be compared with anything on
    the physical scale — a measured pole figure, another ODF, a powder-intensity
    correction. Comparing the raw response directly is a scale error of one to
    two orders of magnitude, not a small one: at a 12 degree halfwidth the
    response of a random texture is about 0.016.

    :class:`~pytex.texture.HarmonicODF` needs no such correction; its pole
    densities are already in m.r.d.

    Method
    ------
    For a random orientation distribution the mapped poles are uniform on the
    sphere, so the expected response is

    .. math:: c = \tfrac{1}{2}\int_{-1}^{1} k(\arccos u)\, \mathrm{d}u ,

    evaluated by Gauss-Legendre quadrature in :math:`u = \cos\omega`.

    Parameters
    ----------
    kernel : KernelSpec
        The kernel the ODF was estimated with.
    antipodal : bool
        Match a response that identifies opposite poles. Friedel's law makes a
        diffraction pole figure blind to the sign of a normal, and a response
        evaluated at ``arccos|cos|`` therefore has a *different* random level:
        the kernel's near lobe is counted for both members of every antipodal
        pair, which asymptotically doubles it. Normalizing a folded response by
        the unfolded mean is a factor-of-two scale error, so this flag must
        match the response being normalized.

    Returns
    -------
    float
        The response of a random texture; strictly positive.
    """

    nodes, weights = np.polynomial.legendre.leggauss(_KERNEL_MEAN_NODES)
    cosines = np.abs(nodes) if antipodal else nodes
    angles = np.arccos(np.clip(cosines, -1.0, 1.0))
    values = np.asarray(kernel.evaluate(angles), dtype=np.float64)
    mean = 0.5 * float(np.sum(values * weights))
    if not np.isfinite(mean) or mean <= 0.0:
        raise ValueError(
            "The ODF smoothing kernel has a non-positive spherical mean, so pole "
            "densities cannot be normalized to multiples of random."
        )
    return mean


#: Largest number of cosine entries held in memory at once while resampling.
#: Resampling forms a (grid x source) cosine matrix; on a fine grid against a
#: large pole cloud that product is far larger than either input, so it is
#: evaluated in blocks. Four million entries is ~32 MB in float64.
_RESAMPLING_BLOCK_ENTRIES = 4_000_000


def _as_direction_array(vectors: np.ndarray | VectorSet) -> np.ndarray:
    if isinstance(vectors, VectorSet):
        return vectors.values
    return vectors


def _spherical_kappa_from_halfwidth(halfwidth_deg: float) -> float:
    """Concentration of the S2 smoothing kernel with the given halfwidth.

    The kernel is ``K(v, u) = exp(kappa * (v.u - 1))``, the von Mises-Fisher
    shape written so that ``K = 1`` at coincidence. Requiring ``K = 1/2`` at an
    angular separation of ``halfwidth_deg`` gives
    ``kappa = ln 2 / (1 - cos(halfwidth))``.
    """

    halfwidth = float(halfwidth_deg)
    if not np.isfinite(halfwidth) or not 0.0 < halfwidth < 180.0:
        raise ValueError("halfwidth_deg must be finite and lie in (0, 180).")
    denominator = 1.0 - float(np.cos(np.deg2rad(halfwidth)))
    if denominator <= 0.0:  # pragma: no cover - guarded by the range check above
        raise ValueError("halfwidth_deg is too small to define a kernel.")
    return float(np.log(2.0) / denominator)


def _spherical_kernel_mean(kappa: float, *, antipodal: bool) -> float:
    """Mean of the S2 smoothing kernel over the sphere.

    Dividing a kernel sum by this value converts it from an unnormalized
    density into multiples of a random distribution, because a uniform
    distribution of poles produces exactly this mean everywhere.

    For ``K(t) = exp(kappa * (t - 1))`` with ``t = cos(angle)``, the spherical
    mean is ``(1/2) * integral of K over t in [-1, 1]``. When opposite poles are
    identified the kernel is evaluated at ``|t|``, which halves the domain and
    doubles the density, giving the second expression.
    """

    if antipodal:
        return float((1.0 - np.exp(-kappa)) / kappa)
    return float((1.0 - np.exp(-2.0 * kappa)) / (2.0 * kappa))


def _resample_directions(
    *,
    target: np.ndarray,
    source: np.ndarray,
    values: np.ndarray,
    kappa: float,
    antipodal: bool,
    estimator: str,
) -> np.ndarray:
    """Evaluate the kernel estimator of ``values`` at ``target`` directions.

    Blocked over target rows so the cosine matrix never has to exist in full.
    """

    kernel_mean = _spherical_kernel_mean(kappa, antipodal=antipodal)
    total_weight = float(np.sum(values))
    block_rows = max(1, _RESAMPLING_BLOCK_ENTRIES // max(1, source.shape[0]))
    estimated = np.empty(target.shape[0], dtype=np.float64)
    for start in range(0, target.shape[0], block_rows):
        stop = min(start + block_rows, target.shape[0])
        cosines = target[start:stop] @ source.T
        if antipodal:
            cosines = np.abs(cosines)
        np.clip(cosines, -1.0, 1.0, out=cosines)
        kernel = np.exp(kappa * (cosines - 1.0))
        if estimator == "density":
            # Kernel density estimate: every pole contributes its own weight,
            # and dividing by the total weight times the kernel's spherical
            # mean puts the result in multiples of random.
            estimated[start:stop] = (kernel @ values) / (total_weight * kernel_mean)
        else:
            # Nadaraya-Watson: a weighted *mean* of nearby samples, which
            # reproduces a constant field exactly and does not inherit the
            # source raster's sampling density.
            denominator = kernel.sum(axis=1)
            supported = denominator > 0.0
            safe = np.where(supported, denominator, 1.0)
            estimated[start:stop] = np.where(supported, (kernel @ values) / safe, 0.0)
    return estimated


def _pole_density_response_matrix(
    dictionary: OrientationSet,
    *,
    pole: CrystalPlane,
    sample_directions: ArrayLike,
    kernel: KernelSpec,
    include_symmetry_family: bool,
    antipodal: bool = False,
) -> np.ndarray:
    """Kernel response of every support orientation at every specimen direction.

    ``antipodal`` folds opposite normals together, which is Friedel's law: a
    diffraction pole figure cannot tell ``h`` from ``-h``, so the response must
    not either. Folding is what makes odd-degree harmonics invisible to the
    forward model, and leaving it out lets an inversion appear to determine an
    odd part that no diffraction measurement can carry. It changes the scale of
    the response, so :func:`random_pole_density` must be called with the same
    flag when converting to multiples of a random distribution.
    """

    direction_array = normalize_vectors(sample_directions)
    if dictionary.crystal_frame != pole.phase.crystal_frame:
        raise ValueError("PoleFigure inversion dictionary must use the pole phase crystal frame.")
    if dictionary.phase is not None and not phases_semantically_match(dictionary.phase, pole.phase):
        raise ValueError("PoleFigure inversion dictionary phase must match PoleFigure.pole.phase.")
    if dictionary.symmetry is not None and dictionary.symmetry != pole.phase.symmetry:
        raise ValueError(
            "PoleFigure inversion dictionary symmetry must match PoleFigure.pole.phase.symmetry."
        )
    pole_family = (
        pole.phase.symmetry.equivalent_vectors(pole.normal)
        if include_symmetry_family
        else pole.normal[None, :]
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
        direction_array,
        mapped_families,
        optimize=True,
    )
    if antipodal:
        cos_angles = np.abs(cos_angles)
    angles = np.arccos(np.clip(cos_angles, -1.0, 1.0))
    response = kernel.evaluate(angles)
    block = np.mean(response, axis=2)
    block = np.ascontiguousarray(block)
    block.setflags(write=False)
    return block


def _orientation_dictionary_response(
    dictionary: OrientationSet,
    pole_figure: PoleFigure,
    kernel: KernelSpec,
    *,
    include_symmetry_family: bool,
) -> np.ndarray:
    if dictionary.specimen_frame != pole_figure.specimen_frame:
        raise ValueError(
            "PoleFigure inversion dictionary specimen_frame must match PoleFigure.specimen_frame."
        )
    return _pole_density_response_matrix(
        dictionary,
        pole=pole_figure.pole,
        sample_directions=pole_figure.sample_directions,
        kernel=kernel,
        include_symmetry_family=include_symmetry_family,
    )


def _project_onto_simplex(vector: np.ndarray) -> np.ndarray:
    r"""Euclidean projection of a vector onto :math:`\{w \ge 0,\ \sum_j w_j = 1\}`.

    Purpose
    -------
    The feasible set of an ODF's dictionary weights. Projected gradient descent
    converges only if the projection is the *nearest* feasible point, and
    clipping negatives then dividing by the sum is not that map: it is a
    radial rescaling, and it makes a Tikhonov term of the form
    :math:`\lambda\|w\|^2` exactly inert, because the term's gradient
    :math:`\lambda w` is parallel to :math:`w` and the rescaling removes any
    change in magnitude. The regularization parameter then has no effect on the
    answer at any value, which is worse than having no parameter.

    Method
    ------
    The standard sort-based algorithm: the projection is
    :math:`\max(w_j - \theta, 0)` for the unique threshold :math:`\theta` making
    the result sum to one, found by scanning the descending sort. Cost is one
    sort per call, negligible beside the gradient evaluation.

    References
    ----------
    Duchi, Shalev-Shwartz, Singer and Chandra, *Efficient projections onto the
    l1-ball for learning in high dimensions*, ICML 2008, Fig. 1.
    """

    descending = np.sort(vector)[::-1]
    cumulative = np.cumsum(descending)
    counts = np.arange(1, vector.size + 1, dtype=np.float64)
    admissible = descending - (cumulative - 1.0) / counts > 0.0
    support = int(np.nonzero(admissible)[0][-1]) + 1 if np.any(admissible) else 1
    threshold = (cumulative[support - 1] - 1.0) / support
    projected: np.ndarray = np.maximum(vector - threshold, 0.0)
    return projected


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
    # Stationarity is measured scale-free, and this is not a refinement. The step
    # length is 1/L, so the raw step size ``||w_next - w||`` is proportional to
    # 1/L: on a system whose operator has large entries the very first step is
    # tiny for that reason alone, and a test of the raw step size against a fixed
    # tolerance then declares convergence immediately and returns the uniform
    # starting guess as the answer. Multiplying the step by L recovers the
    # projected-gradient magnitude, and dividing by ``||A^T b||`` makes the
    # comparison dimensionless, so the same tolerance means the same thing
    # whatever units the pole densities are in.
    gradient_scale = max(float(np.linalg.norm(rhs)), 1e-12)
    history = np.empty(max_iterations, dtype=np.float64)
    converged = False
    for iteration in range(max_iterations):
        gradient = gram @ weights - rhs + regularization * weights
        candidate = _project_onto_simplex(weights - gradient / lipschitz)
        residual = system_matrix @ candidate - observations
        history[iteration] = 0.5 * float(residual @ residual) + 0.5 * regularization * float(
            candidate @ candidate
        )
        stationarity = float(np.linalg.norm(candidate - weights)) * lipschitz / gradient_scale
        weights = candidate
        if stationarity <= tolerance:
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
    """The smoothing kernel and halfwidth used to estimate an ODF.

    Purpose
    -------
    Kernel density estimation on SO(3) needs two decisions, and this makes
    both explicit. The halfwidth is the consequential one: too narrow and the
    ODF reproduces measurement noise, too wide and distinct texture
    components merge. It should reflect the angular resolution of the
    measurement, not be tuned until the figure looks right.

    Attributes
    ----------
    name : str
        ``"de_la_vallee_poussin"`` (default), ``"von_mises_fisher"``,
        ``"gaussian"``, or ``"abel_poisson"``. At equal halfwidth they differ
        mainly in tail weight and in how fast their harmonic coefficients
        decay.
    halfwidth_deg : float
        Angle at which the kernel falls to half its peak value.
    """

    name: str = "de_la_vallee_poussin"
    halfwidth_deg: float = 10.0

    def __post_init__(self) -> None:
        if self.name not in {
            "de_la_vallee_poussin",
            "von_mises_fisher",
            "gaussian",
            "abel_poisson",
        }:
            raise ValueError(
                "Kernel name must be one of 'de_la_vallee_poussin', 'von_mises_fisher', "
                "'gaussian', or 'abel_poisson'."
            )
        if self.halfwidth_deg <= 0.0:
            raise ValueError("Kernel halfwidth must be strictly positive.")

    def evaluate(self, angles_rad: ArrayLike, *, normalized: bool = False) -> np.ndarray:
        """Kernel value at misorientation angles in radians.

        Parameters
        ----------
        angles_rad : ArrayLike
            Misorientation angles.
        normalized : bool
            Scale by the kernel's SO(3) normalization, so that the kernel
            integrates to one over orientation space. Leave it off for relative
            weighting, where only ratios matter.

        See Also
        --------
        as_so3_kernel : The full spectral kernel object, with Chebyshev
            coefficients and bandwidth estimation.
        """

        angle_array = np.asarray(angles_rad, dtype=np.float64)
        halfwidth_rad = np.deg2rad(self.halfwidth_deg)
        if self.name == "de_la_vallee_poussin":
            denominator = np.cos(halfwidth_rad / 2.0)
            exponent = (
                1.0 if np.isclose(denominator, 1.0) else float(np.log(0.5) / np.log(denominator))
            )
            exponent = max(1.0, exponent)
            values = np.clip(np.cos(angle_array / 2.0), 0.0, 1.0) ** exponent
            if normalized:
                dlvp_kernel = self.as_so3_kernel()
                assert isinstance(dlvp_kernel, DeLaValleePoussinKernel)
                values = values * dlvp_kernel.normalization
        elif self.name == "von_mises_fisher":
            kappa = float(np.log(2.0) / (1.0 - np.cos(halfwidth_rad)))
            values = np.exp(kappa * (np.cos(angle_array) - 1.0))
        else:
            kernel = self.as_so3_kernel()
            values = np.asarray(kernel.evaluate(angle_array), dtype=np.float64)
            if not normalized:
                values = values / float(kernel.evaluate(np.zeros(1))[0])
        values = np.ascontiguousarray(values)
        values.setflags(write=False)
        return values

    def as_so3_kernel(
        self,
    ) -> DeLaValleePoussinKernel | GaussianSO3Kernel | AbelPoissonKernel:
        """Return the normalized SO(3) kernel object for this spec.

        Provides the halfwidth <-> bandwidth duality and Chebyshev coefficient
        expansion for every kernel with an SO(3) coefficient representation:
        de la Vallee Poussin (quadrature coefficients), Gaussian
        (Gauss-Weierstrass spectrum), and Abel-Poisson (geometric spectrum).
        The von Mises-Fisher spec remains an S2 evaluation-only kernel.
        """

        if self.name == "gaussian":
            return GaussianSO3Kernel(self.halfwidth_deg)
        if self.name == "abel_poisson":
            return AbelPoissonKernel(self.halfwidth_deg)
        if self.name != "de_la_vallee_poussin":
            raise ValueError(
                "as_so3_kernel() is only defined for the de la Vallee Poussin kernel; "
                f"got '{self.name}'."
            )
        return DeLaValleePoussinKernel(halfwidth_deg=self.halfwidth_deg)

    def bandwidth(self, *, threshold: float = 1e-3, max_bandwidth: int = 512) -> int:
        """Harmonic bandwidth implied by the kernel halfwidth (dVP only)."""

        return self.as_so3_kernel().bandwidth(threshold=threshold, max_bandwidth=max_bandwidth)


@dataclass(frozen=True, slots=True)
class PoleFigure:
    """The distribution of a crystal plane normal over specimen directions.

    Purpose
    -------
    The classical texture representation: where a given ``{hkl}`` points in
    the specimen, as measured by diffraction or computed from orientations.

    Attributes
    ----------
    pole : CrystalPlane
        The plane whose normal is represented.
    sample_directions : np.ndarray
        ``(n, 3)`` specimen-frame directions.
    intensities : np.ndarray
        ``(n,)`` pole densities at those directions.
    specimen_frame : ReferenceFrame
    antipodal : bool
        Whether opposite normals are treated as the same pole — the
        convention that permits a one-hemisphere plot.
    sample_symmetry : SymmetrySpec, optional
        Statistical specimen symmetry imposed, if any. Imposing it is an
        assumption about the process, so it is recorded rather than implied.
    include_symmetry_family : bool
        Whether the whole ``{hkl}`` orbit is represented or only the single
        plane. It decides the correct notation — a family is ``{hkl}``, a
        single plane ``(hkl)`` — which titles and prose read rather than
        assume.
    sampling : str
        How ``intensities`` is to be read: ``"scattered_poles"`` (per-pole
        weights of a pole cloud, the default and what
        :meth:`from_orientations` produces) or ``"sampled_density"`` (pole
        densities already evaluated at ``sample_directions``, what a
        diffractometer raster or an ODF evaluation produces). The two demand
        different resampling estimators, so the reading is recorded rather
        than guessed. See :data:`PoleFigureSampling`.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Two pole figures can be combined arithmetically only once they share a
    support. Use :meth:`on_grid` to resample both onto one
    :class:`~pytex.core.sphere.S2Grid`, and :meth:`normalize_to_mrd` to put
    their magnitudes on the physical multiples-of-random scale, before using
    the arithmetic operators.
    """

    pole: CrystalPlane
    sample_directions: np.ndarray
    intensities: np.ndarray
    specimen_frame: ReferenceFrame
    antipodal: bool = True
    sample_symmetry: SymmetrySpec | None = None
    provenance: ProvenanceRecord | None = None
    #: Whether the plotted poles are the whole symmetry-related orbit of
    #: ``pole`` (the usual case, and what a measured pole figure contains) or
    #: only that single plane. It changes the correct notation for the
    #: quantity: a family is written ``{hkl}``, a single plane ``(hkl)``. Titles
    #: and prose read this rather than assuming.
    includes_symmetry_family: bool = True
    #: How ``intensities`` is to be read; see :data:`PoleFigureSampling`. The
    #: default is the pole-cloud reading, which is what the dominant
    #: constructor :meth:`from_orientations` produces.
    sampling: PoleFigureSampling = "scattered_poles"

    def __post_init__(self) -> None:
        sample_directions = normalize_vectors(self.sample_directions)
        intensities = as_float_array(self.intensities, shape=(sample_directions.shape[0],))
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError("PoleFigure.specimen_frame must belong to the specimen domain.")
        if np.any(~np.isfinite(intensities)) or np.any(intensities < 0.0):
            raise ValueError("PoleFigure intensities must be finite and non-negative.")
        if self.sampling not in _POLE_FIGURE_SAMPLINGS:
            raise ValueError("PoleFigure.sampling must be 'scattered_poles' or 'sampled_density'.")
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
        """Build a pole figure from measured or modelled orientations.

        Purpose
        -------
        Map a crystal plane normal through every orientation into the specimen
        frame — the computation a diffraction pole-figure measurement performs
        physically.

        Parameters
        ----------
        orientations : OrientationSet
            Must share the pole's crystal frame, and its phase when both declare
            one; a mismatch raises.
        pole : CrystalPlane
            The plane whose normal is plotted.
        weights : ArrayLike, optional
            One weight per orientation; uniform when omitted.
        include_symmetry_family : bool
            Plot the whole ``{hkl}`` family (default). A measured pole figure
            always contains the full family, so this is the faithful choice.
        antipodal : bool
            Treat opposite normals as the same pole (default), which is what
            permits a one-hemisphere plot.
        sample_symmetry : SymmetrySpec, optional
            Statistical specimen symmetry to impose. This is an assumption about
            the process, so it is off by default.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        PoleFigure
        """

        if orientations.crystal_frame != pole.phase.crystal_frame:
            raise ValueError("PoleFigure orientations must use the pole phase crystal frame.")
        if orientations.phase is not None and not phases_semantically_match(
            orientations.phase, pole.phase
        ):
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
            includes_symmetry_family=include_symmetry_family,
            # Each row is one pole carrying its orientation's weight, not a
            # density evaluated at that direction.
            sampling="scattered_poles",
        )

    def on_grid(
        self,
        grid: S2Grid,
        *,
        halfwidth_deg: float = DEFAULT_RESAMPLING_HALFWIDTH_DEG,
        estimator: ResamplingEstimator | None = None,
        normalize: bool = True,
        provenance: ProvenanceRecord | None = None,
    ) -> PoleFigure:
        """Resample this pole figure onto a common spherical grid.

        Purpose
        -------
        The prerequisite for comparing or combining two pole figures. Measured
        and computed figures carry whatever support their instrument or their
        orientation set gave them, and two such supports generally share no
        direction at all; arithmetic between them is undefined until both are
        evaluated on one grid. This performs that evaluation.

        When to use
        -----------
        Before :meth:`difference`, the arithmetic operators, or any pointwise
        comparison of two figures — resample both onto the *same*
        :class:`~pytex.core.sphere.S2Grid`. Also useful on its own to convert a
        scattered pole cloud into the smooth density field a contoured plot
        wants, without going through the coarse 2-D binning of
        :meth:`histogram`.

        Method
        ------
        Kernel smoothing on the sphere with
        ``K(v, u) = exp(kappa * (v.u - 1))``, the von Mises-Fisher shape, whose
        concentration follows from ``halfwidth_deg``. Which estimator is correct
        depends on what the intensities mean, so it is taken from
        :attr:`sampling` unless overridden:

        - ``"density"`` (for ``sampling="scattered_poles"``) — kernel density
          estimation, a weighted **sum** over poles, divided by the total weight
          and the kernel's spherical mean so the result is already in multiples
          of random.
        - ``"interpolate"`` (for ``sampling="sampled_density"``) — a
          Nadaraya-Watson weighted **mean** of nearby samples, which reproduces
          a constant field exactly and does not inherit the source raster's own
          sampling density.

        Using the sum where the mean is meant biases the result towards the
        poles of a latitude-longitude raster, which oversamples there; that is
        why the reading is recorded on the object rather than assumed here.

        Parameters
        ----------
        grid : S2Grid
            The target support. Its reference frame must be this figure's
            specimen frame. Prefer :meth:`S2Grid.equispaced` — its equal-area
            weights are what make the m.r.d. normalization and any later
            integration valid.
        halfwidth_deg : float
            Kernel halfwidth, the angle at which the kernel falls to half its
            peak. This is the consequential choice: too small and the result
            reproduces sampling noise, too large and distinct texture
            components merge. Set it from the angular resolution of the
            measurement, not until the figure looks right.
        estimator : str, optional
            Override the estimator implied by :attr:`sampling`. Pass this only
            when the recorded reading is known to be wrong.
        normalize : bool
            Rescale the result so its grid-weighted mean is exactly 1, i.e.
            multiples of a random distribution (default). This is what makes
            two figures numerically comparable. Disable to keep the raw
            estimator output.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        PoleFigure
            A figure on ``grid``'s directions, with ``sampling="sampled_density"``
            because the intensities are now densities evaluated at those
            directions.

        Raises
        ------
        ValueError
            If the grid frame does not match the specimen frame, if the
            estimator name is unknown, or if normalization is requested and the
            resampled field is everywhere zero (nothing to scale).

        Examples
        --------
        A random texture resampled onto an equal-area grid is 1 m.r.d.
        everywhere, to within the kernel's smoothing of the finite pole set;
        that identity is the calibration of this method.
        """

        chosen = "density" if self.sampling == "scattered_poles" else "interpolate"
        if estimator is not None:
            if estimator not in _RESAMPLING_ESTIMATORS:
                raise ValueError("estimator must be 'density' or 'interpolate'.")
            chosen = estimator
        if grid.vectors.reference_frame != self.specimen_frame:
            raise ValueError("PoleFigure.on_grid requires a grid in the figure's specimen frame.")
        if self.intensities.size == 0:
            raise ValueError("PoleFigure.on_grid requires at least one sampled direction.")
        if chosen == "density" and float(np.sum(self.intensities)) <= 0.0:
            raise ValueError("Kernel density estimation requires a positive total pole weight.")

        target = np.asarray(grid.vectors.values, dtype=np.float64)
        estimated = _resample_directions(
            target=target,
            source=self.sample_directions,
            values=self.intensities,
            kappa=_spherical_kappa_from_halfwidth(halfwidth_deg),
            antipodal=self.antipodal,
            estimator=chosen,
        )
        if normalize:
            mean_density = float(np.sum(grid.weights * estimated))
            if not mean_density > 0.0:
                raise ValueError(
                    "Cannot normalize to m.r.d.: the resampled figure is zero everywhere. "
                    "Widen halfwidth_deg, or check that the source directions overlap the grid."
                )
            estimated = estimated / mean_density
        return PoleFigure(
            pole=self.pole,
            sample_directions=target,
            intensities=np.maximum(estimated, 0.0),
            specimen_frame=self.specimen_frame,
            antipodal=self.antipodal,
            sample_symmetry=self.sample_symmetry,
            includes_symmetry_family=self.includes_symmetry_family,
            sampling="sampled_density",
            provenance=self.provenance if provenance is None else provenance,
        )

    def integration_grid(self, *, resolution_deg: float = 5.0) -> S2Grid:
        """An equal-area grid suitable for integrating this figure.

        Purpose: the natural common support for this figure — the hemisphere
        when opposite poles are identified, the whole sphere otherwise — with
        the equal-area weights that make integration over it valid. Use it as
        the ``grid`` argument of :meth:`on_grid` when no particular grid is
        required, so that two figures compared against each other are compared
        on the same one.
        """

        return S2Grid.equispaced(
            resolution_deg,
            reference_frame=self.specimen_frame,
            hemisphere="upper" if self.antipodal else "sphere",
            antipodal=self.antipodal,
        )

    def spherical_mean(
        self,
        *,
        integration_weights: ArrayLike | None = None,
        resolution_deg: float = 5.0,
        halfwidth_deg: float = DEFAULT_RESAMPLING_HALFWIDTH_DEG,
    ) -> float:
        """The mean pole density over the sphere.

        Purpose
        -------
        The quantity that defines the multiples-of-random scale: a pole figure
        is in m.r.d. exactly when this equals 1. It is reported separately from
        :meth:`normalize_to_mrd` so that the scale factor being applied is
        visible rather than implicit.

        Method
        ------
        A mean over the sphere is an integral, and an integral needs solid-angle
        weights, which scattered directions do not carry. Two routes, in order
        of preference:

        - ``integration_weights`` supplied — used directly. This is exact, and
          is the route to take whenever the weights are known: an
          :class:`~pytex.core.sphere.S2Grid` carries them, and a diffractometer
          raster's are proportional to ``sin(tilt)``.
        - Otherwise the figure is resampled onto an equal-area grid and
          integrated there. This is an **estimate**: it inherits the smoothing
          of the resampling kernel, so a very sharp texture sampled coarsely
          will be biased. Supply the weights when accuracy matters.

        Parameters
        ----------
        integration_weights : ArrayLike, optional
            One positive weight per sampled direction. Need not be normalized.
        resolution_deg : float
            Spacing of the fallback integration grid.
        halfwidth_deg : float
            Kernel halfwidth used by the fallback route.

        Returns
        -------
        float
            The mean pole density, in whatever units the intensities carry.
        """

        if integration_weights is not None:
            weights = as_float_array(integration_weights, shape=(self.intensities.shape[0],))
            if np.any(weights <= 0.0):
                raise ValueError("integration_weights must be strictly positive.")
            return float(np.sum(weights * self.intensities) / np.sum(weights))
        grid = self.integration_grid(resolution_deg=resolution_deg)
        resampled = self.on_grid(grid, halfwidth_deg=halfwidth_deg, normalize=False)
        return float(np.sum(grid.weights * resampled.intensities))

    def normalize_to_mrd(
        self,
        *,
        integration_weights: ArrayLike | None = None,
        resolution_deg: float = 5.0,
        halfwidth_deg: float = DEFAULT_RESAMPLING_HALFWIDTH_DEG,
    ) -> PoleFigure:
        """Rescale to multiples of a random distribution.

        Purpose
        -------
        Makes a figure's magnitudes physically meaningful and comparable with
        any other figure. A measured pole figure arrives in detector counts, or
        divided by its own maximum or sum — none of which can be compared
        between two measurements, or against a computed figure, or interpreted
        as "how much stronger than random is this direction". After this, 1
        means random, 2 means twice random, and arithmetic between two figures
        has a defined meaning.

        When to use
        -----------
        On every measured figure before comparing or combining it with another.
        :meth:`on_grid` already normalizes by default, so a resampled figure
        does not need this again.

        Method
        ------
        Divides by :meth:`spherical_mean`; see there for how that mean is
        obtained and when it is exact rather than estimated.

        Parameters
        ----------
        integration_weights : ArrayLike, optional
            Passed to :meth:`spherical_mean`. Supply it when known — the
            normalization is then exact.
        resolution_deg, halfwidth_deg : float
            Fallback integration settings; see :meth:`spherical_mean`.

        Returns
        -------
        PoleFigure
            The same figure with intensities scaled so the mean density is 1.

        Raises
        ------
        ValueError
            If the mean density is not positive, which means there is no scale
            to divide by — an empty or everywhere-zero figure.
        """

        mean_density = self.spherical_mean(
            integration_weights=integration_weights,
            resolution_deg=resolution_deg,
            halfwidth_deg=halfwidth_deg,
        )
        if not mean_density > 0.0:
            raise ValueError("Cannot normalize to m.r.d.: the mean pole density is not positive.")
        return replace(self, intensities=self.intensities / mean_density)

    # ------------------------------------------------------------------
    # Combining figures
    # ------------------------------------------------------------------

    def _require_combinable(self, other: PoleFigure, *, operation: str) -> None:
        """Refuse to combine two figures that do not describe the same quantity.

        Every check here has a way of silently producing a plausible but wrong
        answer if skipped, which is why they raise rather than warn or coerce.
        """

        if self.pole != other.pole:
            raise ValueError(
                f"Cannot {operation} pole figures of different poles: densities of two "
                "different {hkl} are different physical quantities. Compare each pole "
                "against its own counterpart."
            )
        if self.specimen_frame != other.specimen_frame:
            raise ValueError(
                f"Cannot {operation} pole figures in different specimen frames; the same "
                "direction vector means different things in each."
            )
        if self.antipodal != other.antipodal:
            raise ValueError(
                f"Cannot {operation} pole figures with different antipodal conventions: one "
                "identifies opposite normals and the other does not, so their supports are "
                "not the same set of directions."
            )
        if self.includes_symmetry_family != other.includes_symmetry_family:
            raise ValueError(
                f"Cannot {operation} a whole-family {{hkl}} figure with a single-plane (hkl) "
                "figure; they differ by the multiplicity of the family."
            )
        if self.sample_directions.shape != other.sample_directions.shape or not np.allclose(
            self.sample_directions, other.sample_directions, atol=1e-9
        ):
            raise ValueError(
                f"Cannot {operation} pole figures on different supports. Resample both onto "
                "one grid first, e.g. `grid = a.integration_grid(); a.on_grid(grid); "
                "b.on_grid(grid)`."
            )

    def difference(
        self,
        other: PoleFigure,
        *,
        left_label: str = "left",
        right_label: str = "right",
    ) -> PoleFigureDifference:
        """The signed residual figure ``self - other``.

        Purpose
        -------
        The standard comparison of two pole figures: a measurement against the
        figure an ODF recalculates, two measurements of the same specimen, or a
        specimen before and after a treatment. The sign says which way the
        disagreement runs, which is the part a single error number cannot tell
        you.

        Both figures must share a support and should share a scale; put each in
        multiples of random with :meth:`normalize_to_mrd` (or resample with
        :meth:`on_grid`, which normalizes by default) so the residual reads in
        m.r.d. rather than in arbitrary units.

        Parameters
        ----------
        other : PoleFigure
            The figure being subtracted. Must match this one in pole, specimen
            frame, antipodal convention, family flag and support; each mismatch
            raises with the reason.
        left_label, right_label : str
            Names carried into the result so that a plotted or saved residual
            says what was compared.

        Returns
        -------
        PoleFigureDifference
            Read its ``describe()``.
        """

        self._require_combinable(other, operation="subtract")
        return PoleFigureDifference(
            pole=self.pole,
            sample_directions=self.sample_directions,
            values=self.intensities - other.intensities,
            specimen_frame=self.specimen_frame,
            antipodal=self.antipodal,
            left_label=left_label,
            right_label=right_label,
            includes_symmetry_family=self.includes_symmetry_family,
            provenance=self.provenance,
        )

    def __add__(self, other: PoleFigure | float) -> PoleFigure:
        """Sum two pole figures on a shared support, or add a constant.

        Adding figures is how a multi-phase or multi-component pole figure is
        assembled from its parts. The result is a pole density and is therefore
        still a :class:`PoleFigure`.
        """

        if isinstance(other, PoleFigure):
            self._require_combinable(other, operation="add")
            return replace(self, intensities=self.intensities + other.intensities)
        shifted = self.intensities + float(other)
        if np.any(shifted < 0.0):
            raise ValueError(
                "Adding this constant would make pole densities negative, which is not a "
                "pole figure. Use subtraction, which returns a PoleFigureDifference."
            )
        return replace(self, intensities=shifted)

    __radd__ = __add__

    def __sub__(self, other: PoleFigure | float) -> PoleFigureDifference:
        """Signed difference; see :meth:`difference`.

        Subtracting the scalar 1 is the idiomatic way to get the deviation from
        a random distribution, once the figure is in m.r.d.

        Note that this returns a :class:`PoleFigureDifference`, not a
        ``PoleFigure``: the result is signed, and a pole density is not.
        """

        if isinstance(other, PoleFigure):
            return self.difference(other)
        return PoleFigureDifference(
            pole=self.pole,
            sample_directions=self.sample_directions,
            values=self.intensities - float(other),
            specimen_frame=self.specimen_frame,
            antipodal=self.antipodal,
            left_label="figure",
            right_label=f"{float(other):g}",
            includes_symmetry_family=self.includes_symmetry_family,
            provenance=self.provenance,
        )

    def __mul__(self, factor: float) -> PoleFigure:
        """Scale by a non-negative constant.

        Only scalars: the pointwise product of two pole densities is not a pole
        density and has no accepted meaning, so it is not defined. For "how many
        times stronger is this figure than that one" use division.
        """

        value = float(factor)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("PoleFigure can only be scaled by a finite, non-negative factor.")
        return replace(self, intensities=self.intensities * value)

    __rmul__ = __mul__

    def __truediv__(self, other: PoleFigure | float) -> PoleFigure:
        """Ratio to another figure on the same support, or scaling by 1/x.

        The ratio of two pole figures answers "how many times stronger is this
        direction here than there" — the natural comparison when both are in
        m.r.d. It is non-negative, so the result is a :class:`PoleFigure`.

        Raises
        ------
        ValueError
            If any denominator density is zero. A ratio there is undefined, and
            substituting any finite value would invent a comparison the data
            does not support; mask those directions instead.
        """

        if isinstance(other, PoleFigure):
            self._require_combinable(other, operation="divide")
            if np.any(other.intensities <= 0.0):
                raise ValueError(
                    "Cannot divide by a pole figure with zero density: the ratio is undefined "
                    "wherever the denominator vanishes. Restrict both figures to the measured "
                    "region first."
                )
            return replace(self, intensities=self.intensities / other.intensities)
        value = float(other)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("PoleFigure can only be divided by a finite, positive scalar.")
        return replace(self, intensities=self.intensities / value)

    # ------------------------------------------------------------------
    # Transforming a figure
    # ------------------------------------------------------------------

    def rotate(self, rotation: Rotation) -> PoleFigure:
        """Rotate the figure within the specimen frame.

        Purpose: bring two figures measured in different specimen settings into
        a common one — a rolling direction defined differently by two
        laboratories, or a specimen remounted between measurements. Only the
        directions move; the densities travel with them unchanged.

        The rotation is applied to the sampled directions, so the result is no
        longer on the original grid. Resample with :meth:`on_grid` before
        combining it with anything.
        """

        rotated = self.sample_directions @ np.asarray(rotation.as_matrix(), dtype=np.float64).T
        return replace(self, sample_directions=rotated)

    def symmetrize(self, sample_symmetry: SymmetrySpec) -> PoleFigure:
        """Impose a statistical specimen symmetry.

        Purpose
        -------
        Rolled sheet is conventionally assumed orthorhombic about RD/TD/ND, so
        its pole figures are averaged over that symmetry to suppress noise. The
        assumption is about the process, not the data, which is why it is never
        applied implicitly and is recorded on the result.

        Method
        ------
        The support is replicated under every operator of the group and the
        intensities are tiled with it. Under :meth:`on_grid` this yields the
        orbit average for a sampled field and the symmetrized density for a
        pole cloud, so one implementation is correct for both readings.

        Parameters
        ----------
        sample_symmetry : SymmetrySpec
            Must be declared in this figure's specimen frame; a symmetry whose
            operators act in another frame would rotate the poles somewhere
            meaningless.

        Returns
        -------
        PoleFigure
            The replicated figure, carrying the imposed symmetry. Resample it
            with :meth:`on_grid` to obtain the averaged field.
        """

        if sample_symmetry.reference_frame != self.specimen_frame:
            raise ValueError(
                "symmetrize requires a sample_symmetry declared in the figure's specimen frame."
            )
        orbit = sample_symmetry.apply_to_vectors(self.sample_directions)
        directions = np.ascontiguousarray(orbit.reshape(-1, 3), dtype=np.float64)
        return replace(
            self,
            sample_directions=directions,
            intensities=np.tile(self.intensities, sample_symmetry.operators.shape[0]),
            sample_symmetry=sample_symmetry,
        )

    def restrict_polar_range(
        self,
        *,
        max_polar_deg: float,
        min_polar_deg: float = 0.0,
    ) -> PoleFigure:
        """Keep only the directions within a polar-angle band.

        Purpose: an X-ray pole figure is only trustworthy out to the tilt where
        defocusing takes over, typically 70-80 degrees. Everything beyond it is
        instrument, not texture. Discarding it explicitly keeps that
        unmeasured-versus-zero distinction visible instead of letting the
        unreliable rim drive a normalization or a residual.

        The polar angle is measured from the specimen-frame ``+Z`` axis, and
        from the nearer of ``+Z``/``-Z`` when opposite poles are identified.

        Returns
        -------
        PoleFigure
            The retained subset, in the original order.

        Raises
        ------
        ValueError
            If the band is empty or leaves no sampled direction — an empty pole
            figure is not a useful object to hand onwards.
        """

        low, high = float(min_polar_deg), float(max_polar_deg)
        if not 0.0 <= low < high <= 180.0:
            raise ValueError("Require 0 <= min_polar_deg < max_polar_deg <= 180.")
        components = self.sample_directions[:, 2]
        if self.antipodal:
            components = np.abs(components)
        polar_deg = np.degrees(np.arccos(np.clip(components, -1.0, 1.0)))
        keep = (polar_deg >= low) & (polar_deg <= high)
        if not np.any(keep):
            raise ValueError(
                f"No sampled direction lies between {low} and {high} degrees of polar angle."
            )
        return replace(
            self,
            sample_directions=self.sample_directions[keep],
            intensities=self.intensities[keep],
        )

    def project(self, *, method: str = "equal_area") -> np.ndarray:
        """Project the specimen directions onto the plotting plane.

        Parameters
        ----------
        method : str
            ``"equal_area"`` (Schmidt/Lambert, default) preserves area, so
            densities read directly off the plot and pole-density comparisons are
            fair. ``"stereographic"`` (Wulff) preserves angles instead and is the
            right choice for angle-measuring constructions.

        Returns
        -------
        np.ndarray
            ``(n, 2)`` plane coordinates. The projection disc has radius
            ``sqrt(2)`` for equal-area and ``1`` for stereographic.
        """

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
        """Binned pole density on the projection plane.

        Purpose
        -------
        Convert scattered poles into the density field that a contoured pole
        figure displays.

        Parameters
        ----------
        bins : int
            Bins per axis of the square binning grid.
        method : str
            Projection method; see :meth:`project`. Equal-area projection is
            strongly preferred here, since only then do equal bins represent
            equal solid angles and the densities become comparable.

        Returns
        -------
        tuple of np.ndarray
            Intensity-weighted 2-D histogram and the two bin-edge arrays, as
            from ``numpy.histogram2d``.
        """

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
class PoleFigureDifference:
    """The signed difference between two pole figures on a shared support.

    Purpose
    -------
    The product that answers "where, and by how much, do these two figures
    disagree?" — the residual figure that turns an ODF inversion from a number
    into a diagnosis, and the standard way to compare a measurement against a
    model or against another measurement.

    Why this is not a :class:`PoleFigure`
    -------------------------------------
    A pole density is non-negative, and ``PoleFigure`` enforces that. A
    difference is signed, and its sign is the whole content: which regions the
    model over-predicts and which it under-predicts. Coercing it into a
    non-negative type would destroy exactly the information it exists to carry,
    so subtraction returns this instead. Addition, scaling and ratios stay
    non-negative and do return ``PoleFigure``.

    Attributes
    ----------
    pole : CrystalPlane
        The plane both figures represent.
    sample_directions : np.ndarray
        ``(n, 3)`` shared specimen-frame support.
    values : np.ndarray
        ``(n,)`` signed differences, ``left - right``. Positive means the left
        figure is the larger.
    specimen_frame : ReferenceFrame
    antipodal : bool
    left_label, right_label : str
        What was subtracted from what, so a saved or plotted residual is
        self-describing.
    includes_symmetry_family : bool
    provenance : ProvenanceRecord, optional
    """

    pole: CrystalPlane
    sample_directions: np.ndarray
    values: np.ndarray
    specimen_frame: ReferenceFrame
    antipodal: bool = True
    left_label: str = "left"
    right_label: str = "right"
    includes_symmetry_family: bool = True
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        sample_directions = normalize_vectors(self.sample_directions)
        values = as_float_array(self.values, shape=(sample_directions.shape[0],))
        if self.specimen_frame.domain is not FrameDomain.SPECIMEN:
            raise ValueError(
                "PoleFigureDifference.specimen_frame must belong to the specimen domain."
            )
        if np.any(~np.isfinite(values)):
            raise ValueError("PoleFigureDifference values must be finite.")
        object.__setattr__(self, "sample_directions", sample_directions)
        object.__setattr__(self, "values", values)

    def __len__(self) -> int:
        return int(self.values.shape[0])

    @property
    def max_absolute_deviation(self) -> float:
        """Largest disagreement anywhere on the figure, in m.r.d."""

        return float(np.max(np.abs(self.values)))

    @property
    def mean_deviation(self) -> float:
        """Unweighted mean signed difference.

        Compare it against :attr:`rms_deviation`, not against zero. A shape
        disagreement averages out, so its mean is a small fraction of its RMS; a
        normalization error is a constant offset, so its mean accounts for most
        of the RMS. Half is the dividing line :meth:`describe` uses.
        """

        return float(np.mean(self.values))

    @property
    def rms_deviation(self) -> float:
        """Root-mean-square disagreement, the usual single-number summary."""

        return float(np.sqrt(np.mean(self.values**2)))

    def weighted_rms_deviation(self, weights: ArrayLike) -> float:
        """Solid-angle-weighted RMS disagreement.

        Prefer this over :attr:`rms_deviation` whenever the support is a raster
        or any other non-equal-area sampling: an unweighted RMS over such a
        support is dominated by wherever the sampling happens to be densest.
        """

        weight_array = as_float_array(weights, shape=(self.values.shape[0],))
        if np.any(weight_array <= 0.0):
            raise ValueError("weights must be strictly positive.")
        total = float(np.sum(weight_array))
        return float(np.sqrt(np.sum(weight_array * self.values**2) / total))

    def project(self, *, method: str = "equal_area") -> np.ndarray:
        """Project the support onto the plotting plane; see :meth:`PoleFigure.project`."""

        return project_directions(
            self.sample_directions,
            method=method,
            antipodal=self.antipodal,
        )

    def describe(self) -> str:
        """Prose summary: what was compared, how badly it disagrees, and where."""

        indices = [int(value) for value in self.pole.miller.indices]
        notation = (
            format_plane_family_indices(indices, style="plain")
            if self.includes_symmetry_family
            else format_plane_indices(indices, style="plain")
        )
        worst = int(np.argmax(np.abs(self.values)))
        direction = self.sample_directions[worst]
        signed = float(self.values[worst])
        bias = (
            "The mean signed difference accounts for most of the deviation, which indicates the "
            "two figures are not on the same normalization; put both in multiples of random with "
            "normalize_to_mrd before reading the shape of the residual."
            # A constant offset puts all of the deviation into the mean, so the
            # ratio approaches one; a shape disagreement averages out and leaves
            # the mean a small fraction of the RMS. The bar was 5 percent of the
            # RMS, which a well-fitted figure clears routinely -- a residual of
            # 0.002 m.r.d. on a field of order 1 was being reported as a scale
            # error. Half puts the test where the distinction actually is.
            if abs(self.mean_deviation) > 0.5 * max(self.rms_deviation, 1e-12)
            else "The mean signed difference is small next to the RMS deviation, so the two "
            "figures share a scale and the residual describes shape disagreement rather than a "
            "normalization error."
        )
        return (
            f"Difference of {notation} pole figures, {self.left_label} minus {self.right_label}, "
            f"over {len(self)} shared specimen directions. RMS deviation "
            f"{self.rms_deviation:.4f} m.r.d., maximum |deviation| "
            f"{self.max_absolute_deviation:.4f} m.r.d. at specimen direction "
            f"[{direction[0]:.3f} {direction[1]:.3f} {direction[2]:.3f}], where "
            f"{self.left_label} is "
            f"{'higher' if signed > 0.0 else 'lower'} by {abs(signed):.4f} m.r.d. "
            f"Mean signed difference {self.mean_deviation:+.4f} m.r.d. {bias}"
        )


@dataclass(frozen=True, slots=True)
class InversePoleFigure:
    """The distribution of a specimen direction over crystal directions.

    Purpose
    -------
    The complement of a pole figure: which crystal directions align with a
    chosen specimen axis. This is the representation behind IPF colouring and
    behind fibre-texture statements such as "a strong ``<111>`` fibre along
    ND".

    Attributes
    ----------
    sample_direction : np.ndarray
        The specimen axis; normalized on construction.
    crystal_directions : np.ndarray
        ``(n, 3)`` crystal-frame directions, folded into the fundamental
        sector when symmetry reduction was requested.
    intensities : np.ndarray
        ``(n,)`` densities.
    crystal_frame, specimen_frame : ReferenceFrame
    antipodal : bool
    crystal_symmetry : SymmetrySpec, optional
        Needed to draw the standard triangle; without it the figure is an
        unbounded scatter.
    provenance : ProvenanceRecord, optional
    """

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
        """Build an inverse pole figure from measured or modelled orientations.

        Purpose
        -------
        Map a fixed specimen direction through every orientation into the crystal
        frame — the representation behind IPF colouring and behind fibre-texture
        statements such as "a <111> fibre along ND".

        Parameters
        ----------
        orientations : OrientationSet
        sample_direction : ArrayLike
            The specimen axis; normalized internally.
        weights : ArrayLike, optional
            One weight per orientation; uniform when omitted.
        reduce_by_symmetry : bool
            Fold the resulting crystal directions into the symmetry fundamental
            sector (default) — what makes the standard triangle standard.
        antipodal : bool
            Treat a direction and its reverse as equivalent (default).
        provenance : ProvenanceRecord, optional

        Returns
        -------
        InversePoleFigure
        """

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
        """Project the crystal directions onto the plotting plane.

        See :meth:`PoleFigure.project` for the projection methods and their
        respective uses.
        """

        return project_directions(
            self.crystal_directions,
            method=method,
            antipodal=self.antipodal,
        )

    @property
    def sector_vertices(self) -> np.ndarray | None:
        """Corner directions of the symmetry fundamental sector, or ``None``.

        ``None`` when no crystal symmetry is attached, in which case there is no
        standard triangle to draw.
        """

        if self.crystal_symmetry is None:
            return None
        return self.crystal_symmetry.fundamental_sector(antipodal=self.antipodal).vertices

    def project_sector_vertices(self, *, method: str = "equal_area") -> np.ndarray | None:
        """The fundamental-sector corners projected onto the plotting plane.

        Used to draw the standard-triangle outline in the same projection as the
        data. ``None`` when no crystal symmetry is attached.
        """

        if self.sector_vertices is None:
            return None
        return project_directions(
            self.sector_vertices,
            method=method,
            antipodal=self.antipodal,
        )


@dataclass(frozen=True, slots=True)
class ODFSectionData:
    """ODF density sampled on constant-coordinate Bunge-Euler sections.

    ``section_kind`` names the sectioning coordinate: ``"phi2"`` (the
    default, constant-phi2 sections) or ``"sigma"`` (constant
    sigma = phi1 + phi2 sections). Each section is a density grid over
    (phi1, Phi); ``phi2_deg`` holds the per-section value of the sectioning
    coordinate (aliased as `section_values_deg`).
    """

    phi2_deg: np.ndarray
    phi1_deg: np.ndarray
    big_phi_deg: np.ndarray
    densities: np.ndarray
    normalized: bool = False
    section_kind: str = "phi2"

    def __post_init__(self) -> None:
        for name in ("phi2_deg", "phi1_deg", "big_phi_deg", "densities"):
            array = np.ascontiguousarray(np.asarray(getattr(self, name), dtype=np.float64))
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        expected = (self.phi2_deg.shape[0], self.big_phi_deg.shape[0], self.phi1_deg.shape[0])
        if self.densities.shape != expected:
            raise ValueError(
                "ODFSectionData.densities must have shape (n_sections, n_big_phi, n_phi1)."
            )
        if self.section_kind not in {"phi2", "sigma"}:
            raise ValueError("ODFSectionData.section_kind must be 'phi2' or 'sigma'.")

    @property
    def section_values_deg(self) -> np.ndarray:
        """Per-section value of the sectioning coordinate (phi2 or sigma)."""

        return self.phi2_deg

    @property
    def section_count(self) -> int:
        """Number of constant-coordinate sections in this data set."""

        return int(self.phi2_deg.shape[0])

    @property
    def max_density(self) -> float:
        """Largest ODF density over all sections, in multiples of random.

        The number quoted as texture strength, and the natural upper limit for a
        shared contour scale across sections.
        """

        return float(np.max(self.densities))


@dataclass(frozen=True, slots=True)
class ODF:
    """A discrete kernel-density orientation distribution function.

    Purpose
    -------
    The orientation distribution represented as weighted support
    orientations convolved with a smoothing kernel — the natural
    representation for EBSD data, where every indexed point is one support
    orientation. Densities are in multiples of a random distribution, so a
    value of 1 means untextured.

    Attributes
    ----------
    orientations : OrientationSet
        The support. Its resolution bounds the angular detail representable.
    weights : np.ndarray
        One non-negative weight per support orientation, summing to a
        positive value.
    kernel : KernelSpec
        The smoothing kernel and halfwidth.
    specimen_symmetry : SymmetrySpec, optional
        Statistical specimen symmetry imposed, if any.
    provenance : ProvenanceRecord, optional

    See Also
    --------
    pytex.texture.HarmonicODF : The series-expansion representation, compact
        for symmetric materials and the natural form for convolution.
    """

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
        """Build a discrete kernel-density ODF from an orientation set.

        Purpose
        -------
        The standard route from measured orientations to a continuous
        orientation distribution: each orientation contributes a smoothing
        kernel, and the sum is the estimated density.

        Parameters
        ----------
        orientations : OrientationSet
            The support. For EBSD data this is one orientation per indexed point.
        weights : ArrayLike, optional
            One weight per orientation — a confidence index or a grain area, for
            example. Uniform when omitted.
        kernel : KernelSpec, optional
            Smoothing kernel and halfwidth. The halfwidth is the central choice:
            too narrow and the ODF reproduces measurement noise, too wide and
            real texture components are smeared together. It should reflect the
            angular resolution of the measurement, not be tuned for appearance.
        specimen_symmetry : SymmetrySpec, optional
            Statistical specimen symmetry to impose; an assumption about the
            process, so off by default.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        ODF
        """

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
        """Support weights normalized to sum to one, read-only.

        Makes the ODF a probability density over orientation space, so that
        volume fractions read directly as fractions.
        """

        normalized = self.weights / np.sum(self.weights)
        normalized = np.ascontiguousarray(normalized)
        normalized.setflags(write=False)
        return normalized

    def evaluate(
        self,
        orientations: Orientation | OrientationSet,
        *,
        symmetry_aware: bool = True,
        normalized: bool = False,
    ) -> np.ndarray | float:
        """Evaluate the ODF density at given orientations.

        Purpose
        -------
        The density ``f(g)`` in multiples of a random distribution — the number a
        texture-strength statement quotes.

        Parameters
        ----------
        orientations : Orientation or OrientationSet
            Query orientations. Must share the ODF support's crystal and specimen
            frames. A single orientation returns a scalar.
        symmetry_aware : bool
            Reduce query-to-support misorientations by crystal symmetry
            (default). Turning it off makes the density depend on which symmetry
            branch the query happens to be written in, so leave it on unless both
            sides are already reduced.
        normalized : bool
            Use the SO(3)-normalized kernel, so densities are in multiples of
            random rather than in arbitrary units.

        Returns
        -------
        float or np.ndarray
        """

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
            if not phases_semantically_match(query_set.phase, self.orientations.phase):
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
        kernel_values = self.kernel.evaluate(angles, normalized=normalized)
        density = kernel_values @ self.normalized_weights
        density = np.ascontiguousarray(density)
        density.setflags(write=False)
        if scalar_output:
            return float(density[0])
        return density

    def phi2_sections(
        self,
        *,
        phi2_deg: ArrayLike = (0.0, 45.0, 65.0),
        phi1_max_deg: float = 90.0,
        big_phi_max_deg: float = 90.0,
        resolution_deg: float = 5.0,
        normalized: bool = False,
    ) -> ODFSectionData:
        """Sample the ODF density on constant-phi2 Bunge-Euler sections.

        Returns an `ODFSectionData` with the density on a
        ``(n_sections, n_big_phi, n_phi1)`` grid — the standard texture
        visualisation. Ranges default to the cubic ``[0, 90] deg`` domain; widen
        them for lower-symmetry crystals.
        """

        if resolution_deg <= 0.0:
            raise ValueError("resolution_deg must be strictly positive.")
        phi2_values = np.atleast_1d(np.asarray(phi2_deg, dtype=np.float64))
        phi1 = np.arange(0.0, phi1_max_deg + 1e-9, resolution_deg)
        big_phi = np.arange(0.0, big_phi_max_deg + 1e-9, resolution_deg)
        grid_phi1, grid_big_phi = np.meshgrid(phi1, big_phi, indexing="xy")
        sections = []
        for phi2 in phi2_values:
            euler = np.column_stack(
                [
                    grid_phi1.ravel(),
                    grid_big_phi.ravel(),
                    np.full(grid_phi1.size, float(phi2)),
                ]
            )
            query = OrientationSet.from_euler_angles(
                euler,
                crystal_frame=self.orientations.crystal_frame,
                specimen_frame=self.orientations.specimen_frame,
                symmetry=self.orientations.symmetry,
                phase=self.orientations.phase,
                convention="bunge",
                degrees=True,
            )
            density = np.asarray(
                self.evaluate(query, normalized=normalized), dtype=np.float64
            ).reshape(grid_big_phi.shape)
            sections.append(density)
        return ODFSectionData(
            phi2_deg=phi2_values,
            phi1_deg=phi1,
            big_phi_deg=big_phi,
            densities=np.stack(sections, axis=0),
            normalized=normalized,
        )

    def sigma_sections(
        self,
        *,
        sigma_deg: ArrayLike = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
        phi1_max_deg: float = 90.0,
        big_phi_max_deg: float = 90.0,
        resolution_deg: float = 5.0,
        normalized: bool = False,
    ) -> ODFSectionData:
        """Sample the ODF density on constant-sigma Bunge-Euler sections.

        Sigma sections fix ``sigma = phi1 + phi2`` and show the density over
        the (phi1, Phi) plane with ``phi2 = sigma - phi1`` (mod 360 deg) —
        the classic view for gamma-fibre-dominated (e.g. bcc rolling)
        textures. Returns an `ODFSectionData` with
        ``section_kind = "sigma"``.
        """

        if resolution_deg <= 0.0:
            raise ValueError("resolution_deg must be strictly positive.")
        sigma_values = np.atleast_1d(np.asarray(sigma_deg, dtype=np.float64))
        phi1 = np.arange(0.0, phi1_max_deg + 1e-9, resolution_deg)
        big_phi = np.arange(0.0, big_phi_max_deg + 1e-9, resolution_deg)
        grid_phi1, grid_big_phi = np.meshgrid(phi1, big_phi, indexing="xy")
        sections = []
        for sigma in sigma_values:
            phi2_flat = np.mod(float(sigma) - grid_phi1.ravel(), 360.0)
            euler = np.column_stack(
                [
                    grid_phi1.ravel(),
                    grid_big_phi.ravel(),
                    phi2_flat,
                ]
            )
            query = OrientationSet.from_euler_angles(
                euler,
                crystal_frame=self.orientations.crystal_frame,
                specimen_frame=self.orientations.specimen_frame,
                symmetry=self.orientations.symmetry,
                phase=self.orientations.phase,
                convention="bunge",
                degrees=True,
            )
            density = np.asarray(
                self.evaluate(query, normalized=normalized), dtype=np.float64
            ).reshape(grid_big_phi.shape)
            sections.append(density)
        return ODFSectionData(
            phi2_deg=sigma_values,
            phi1_deg=phi1,
            big_phi_deg=big_phi,
            densities=np.stack(sections, axis=0),
            normalized=normalized,
            section_kind="sigma",
        )

    def reconstruct_pole_figure(
        self,
        pole: CrystalPlane,
        *,
        include_symmetry_family: bool = True,
        antipodal: bool = True,
    ) -> PoleFigure:
        """The pole figure this ODF predicts for a given crystal plane.

        Purpose
        -------
        The forward projection from orientation space to pole space. Comparing a
        reconstructed pole figure against the measured one is the standard
        check on an ODF inversion — a good ODF must reproduce the data it came
        from.

        Parameters
        ----------
        pole : CrystalPlane
        include_symmetry_family : bool
            Include the whole ``{hkl}`` family (default), as a measurement does.
        antipodal : bool
            Treat opposite normals as the same pole (default).
        """

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
        """Reconstructed pole figures for several planes.

        The batch form of :meth:`reconstruct_pole_figure`, applying the same
        conventions to every pole so the set is internally consistent.
        """

        return tuple(
            self.reconstruct_pole_figure(
                pole,
                include_symmetry_family=include_symmetry_family,
                antipodal=antipodal,
            )
            for pole in poles
        )

    def evaluate_pole_density(
        self,
        pole: CrystalPlane,
        sample_directions: ArrayLike,
        *,
        include_symmetry_family: bool = True,
        antipodal: bool = False,
    ) -> np.ndarray:
        """Pole density at chosen specimen directions, without binning.

        Purpose
        -------
        Evaluate the pole figure exactly where it is wanted — along a fibre, at a
        measured pole position, or on a supplied grid — instead of reconstructing
        a whole figure and interpolating.

        Parameters
        ----------
        pole : CrystalPlane
            Must share the ODF support's crystal frame and phase.
        sample_directions : ArrayLike
            ``(n, 3)`` specimen directions.
        include_symmetry_family : bool
            Include the whole ``{hkl}`` family (default).
        antipodal : bool
            Identify opposite normals, as Friedel's law does. Off by default,
            because this method returns the raw kernel response and the caller
            owns the normalization: switching it on without also passing
            ``antipodal=True`` to :func:`random_pole_density` is a
            factor-of-two scale error.
            :meth:`~pytex.texture.HarmonicODF.evaluate_pole_density`, which
            normalizes internally and therefore owns the convention itself,
            folds by default.

        Returns
        -------
        np.ndarray
            ``(n,)`` pole densities.
        """

        if self.orientations.crystal_frame != pole.phase.crystal_frame:
            raise ValueError("Pole density evaluation requires the same crystal frame as the pole.")
        if self.orientations.phase is not None and not phases_semantically_match(
            self.orientations.phase, pole.phase
        ):
            raise ValueError("Pole density evaluation requires the same phase as the ODF support.")
        response = _pole_density_response_matrix(
            self.orientations,
            pole=pole,
            sample_directions=sample_directions,
            kernel=self.kernel,
            include_symmetry_family=include_symmetry_family,
            antipodal=antipodal,
        )
        density = response @ self.normalized_weights
        density = np.ascontiguousarray(density, dtype=np.float64)
        density.setflags(write=False)
        return density

    def volume_fraction(
        self,
        center: Orientation,
        *,
        max_angle_deg: float,
        symmetry_aware: bool = True,
    ) -> float:
        """Fraction of the texture lying within a given angle of an orientation.

        Purpose
        -------
        The quantitative statement behind "the cube component makes up 18 percent
        of the texture": the summed normalized weight of all support orientations
        within ``max_angle_deg`` of the component centre.

        Parameters
        ----------
        center : Orientation
            The component's ideal orientation.
        max_angle_deg : float
            Angular radius defining the component. The result depends strongly on
            this radius, so it must be reported alongside the fraction for the
            number to mean anything.
        symmetry_aware : bool
            Use symmetry-reduced disorientation angles (default).

        Returns
        -------
        float
            A fraction in ``[0, 1]``.

        Notes
        -----
        Computed by hard cut-off on the discrete support, not by integrating the
        smoothed density, so it is a support-weight fraction rather than a
        kernel-integrated volume. The two agree closely when the support is dense
        relative to the kernel halfwidth.
        """

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
        """Estimate an ODF from measured pole figures (PF-to-ODF inversion).

        Purpose
        -------
        The classical inverse problem of quantitative texture analysis: X-ray and
        neutron diffraction measure pole figures, but the orientation
        distribution is what physical models need.

        Method and limits
        -----------------
        Builds the pole-density response of every dictionary orientation to every
        measured pole-figure point and solves the resulting non-negative,
        regularized least-squares system iteratively. The problem is
        ill-posed — pole figures are projections and lose the odd-order harmonic
        information — so the solution depends on the dictionary, the kernel, and
        the regularization, and several pole figures from different planes are
        needed to constrain it. Ghost correction and zero-range methods are not
        applied here; treat the result as a regularized estimate, not a unique
        inversion.

        Parameters
        ----------
        pole_figures : sequence of PoleFigure
            At least one; all must share a specimen frame with the dictionary.
            More independent poles give a better-conditioned problem. The
            intensities are taken to be pole densities in **multiples of a random
            distribution**, which is what a measured figure carries and what
            :meth:`PoleFigure.on_grid` produces. A scattered pole cloud whose
            intensities are per-pole *weights* is not on that scale; resample it
            onto a grid first.
        orientation_dictionary : OrientationSet
            The support the ODF is expressed on; its resolution bounds the
            achievable angular detail.
        kernel : KernelSpec, optional
            Smoothing kernel for the response model.
        regularization : float
            Tikhonov weight. Larger values give smoother, more stable, less
            detailed solutions.
        include_symmetry_family : bool
            Model the whole ``{hkl}`` family per pole figure (default), matching
            what the measurement contains.
        max_iterations : int
            Iteration cap for the solver.
        tolerance : float
            Convergence tolerance.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        ODFInversionReport
            The estimated ODF together with the residuals and convergence
            information needed to judge whether the fit is usable.
        """

        if not pole_figures:
            raise ValueError("ODF inversion requires at least one PoleFigure.")
        specimen_frame = orientation_dictionary.specimen_frame
        for pole_figure in pole_figures:
            if pole_figure.specimen_frame != specimen_frame:
                raise ValueError(
                    "All pole figures and the inversion dictionary must share a specimen frame."
                )
        inversion_kernel = KernelSpec() if kernel is None else kernel
        # The observations are pole densities in multiples of a random
        # distribution -- that is what a measured figure carries and what
        # PoleFigure.on_grid produces -- while the dictionary response is the raw
        # kernel sum, whose value for a random texture is the kernel's spherical
        # mean rather than one. The operator is put on the observations' scale
        # here. Without this the system is unfittable rather than merely
        # mis-scaled: the weights are constrained to sum to one, so the model
        # cannot absorb a factor of 64 (at a 12 degree halfwidth) into its
        # amplitude, and the solver stalls at a relative residual near 1 while
        # reporting convergence.
        kernel_mean = random_pole_density(inversion_kernel)
        blocks = [
            _orientation_dictionary_response(
                orientation_dictionary,
                pole_figure,
                inversion_kernel,
                include_symmetry_family=include_symmetry_family,
            )
            / kernel_mean
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
        predicted = np.ascontiguousarray(system_matrix @ weights, dtype=np.float64)
        residual = predicted - observations
        residual_norm = float(np.linalg.norm(residual))
        common_sample_symmetry = pole_figures[0].sample_symmetry
        if any(
            pole_figure.sample_symmetry != common_sample_symmetry for pole_figure in pole_figures
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
            residual_norm=residual_norm,
            objective_history=objective_history,
            iterations=int(objective_history.shape[0]),
            converged=converged,
            regularization=regularization,
            observation_count=int(observations.size),
            dictionary_size=len(orientation_dictionary),
            relative_residual_norm=float(
                residual_norm / max(float(np.linalg.norm(observations)), 1e-12)
            ),
            mean_absolute_error=float(np.mean(np.abs(residual))),
            max_absolute_error=float(np.max(np.abs(residual))),
            predicted_intensities=predicted,
            dictionary_coverage_ratio=float(
                observations.size / max(1, len(orientation_dictionary))
            ),
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class ODFInversionReport:
    """A PF-to-ODF inversion result together with the evidence to judge it.

    Purpose
    -------
    Pole-figure inversion is ill-posed, so an ODF alone is not a result. This
    carries the residuals, the convergence history, and the problem
    dimensions, so a reader can see whether the solution actually reproduces
    the measured data and how strongly it was regularized.

    Attributes
    ----------
    odf : ODF
        The estimated distribution.
    residual_norm, relative_residual_norm : float
        Absolute and normalized misfit against the measured pole densities.
    objective_history : np.ndarray
        Objective value per iteration; a history still descending at the end
        means the solver stopped early rather than converged.
    iterations : int
    converged : bool
    regularization : float
        The Tikhonov weight used; larger means smoother and less detailed.
    observation_count, dictionary_size : int
        The problem dimensions. A dictionary much larger than the
        observation count is underdetermined and leans on regularization.
    mean_absolute_error, max_absolute_error : float
    predicted_intensities : np.ndarray
        Pole densities the solution implies, for direct comparison against
        the measurement.
    provenance : ProvenanceRecord, optional
    """

    odf: ODF
    residual_norm: float
    objective_history: np.ndarray
    iterations: int
    converged: bool
    regularization: float
    observation_count: int
    dictionary_size: int
    relative_residual_norm: float
    mean_absolute_error: float
    max_absolute_error: float
    predicted_intensities: np.ndarray
    dictionary_coverage_ratio: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        history = as_float_array(self.objective_history, shape=(None,))
        predicted = as_float_array(self.predicted_intensities, shape=(self.observation_count,))
        if history.size == 0:
            raise ValueError("ODFInversionReport.objective_history must not be empty.")
        if self.residual_norm < 0.0:
            raise ValueError("ODFInversionReport.residual_norm must be non-negative.")
        if self.relative_residual_norm < 0.0:
            raise ValueError("ODFInversionReport.relative_residual_norm must be non-negative.")
        if self.mean_absolute_error < 0.0:
            raise ValueError("ODFInversionReport.mean_absolute_error must be non-negative.")
        if self.max_absolute_error < 0.0:
            raise ValueError("ODFInversionReport.max_absolute_error must be non-negative.")
        if self.iterations <= 0:
            raise ValueError("ODFInversionReport.iterations must be strictly positive.")
        if self.regularization < 0.0:
            raise ValueError("ODFInversionReport.regularization must be non-negative.")
        if self.observation_count <= 0:
            raise ValueError("ODFInversionReport.observation_count must be strictly positive.")
        if self.dictionary_size <= 0:
            raise ValueError("ODFInversionReport.dictionary_size must be strictly positive.")
        if self.dictionary_coverage_ratio <= 0.0:
            raise ValueError("ODFInversionReport.dictionary_coverage_ratio must be positive.")
        object.__setattr__(self, "objective_history", history)
        object.__setattr__(self, "predicted_intensities", predicted)
