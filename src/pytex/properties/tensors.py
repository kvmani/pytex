"""Elastic stiffness and compliance tensors with directional properties.

`StiffnessTensor` and `ComplianceTensor` store the full rank-4 elastic tensor in
crystal-frame Cartesian coordinates and round-trip to the 6x6 Voigt matrix.
Stiffness uses no Voigt factors; compliance carries the standard 1/2 and 1/4
factors on shear terms so that ``S = C^-1`` in Voigt form is consistent with the
tensor contraction ``1/E(n) = n_i n_j n_k n_l S_ijkl``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike

if TYPE_CHECKING:
    from pytex.core.orientation import OrientationSet

# Voigt index pairs: Voigt 0..5 -> (i, j) tensor indices.
_VOIGT_PAIRS: tuple[tuple[int, int], ...] = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


def _voigt_index(i: int, j: int) -> int:
    return _VOIGT_PAIRS.index((i, j) if (i, j) in _VOIGT_PAIRS else (j, i))


def _compliance_factor(voigt_index: int) -> float:
    return 1.0 if voigt_index < 3 else 2.0


def _voigt_to_tensor(matrix: np.ndarray, *, compliance: bool) -> np.ndarray:
    tensor = np.zeros((3, 3, 3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            capital = _voigt_index(i, j)
            for k in range(3):
                for m in range(3):
                    lower = _voigt_index(k, m)
                    value = matrix[capital, lower]
                    if compliance:
                        value /= _compliance_factor(capital) * _compliance_factor(lower)
                    tensor[i, j, k, m] = value
    return tensor


def _tensor_to_voigt(tensor: np.ndarray, *, compliance: bool) -> np.ndarray:
    matrix = np.zeros((6, 6), dtype=np.float64)
    for capital, (i, j) in enumerate(_VOIGT_PAIRS):
        for lower, (k, m) in enumerate(_VOIGT_PAIRS):
            value = tensor[i, j, k, m]
            if compliance:
                value *= _compliance_factor(capital) * _compliance_factor(lower)
            matrix[capital, lower] = value
    return matrix


@dataclass(frozen=True, slots=True)
class ElasticTensor:
    """Base rank-4 elastic tensor (crystal-frame Cartesian, units of GPa or 1/GPa)."""

    tensor: np.ndarray
    _compliance: bool = False

    def __post_init__(self) -> None:
        array = np.ascontiguousarray(np.asarray(self.tensor, dtype=np.float64))
        if array.shape != (3, 3, 3, 3):
            raise ValueError("ElasticTensor.tensor must have shape (3, 3, 3, 3).")
        array.setflags(write=False)
        object.__setattr__(self, "tensor", array)

    def voigt_matrix(self) -> np.ndarray:
        """The ``6x6`` Voigt matrix of this fourth-rank tensor, read-only.

        The compact engineering form in which elastic constants are tabulated.
        Note that the Voigt convention differs between stiffness and compliance —
        compliance carries factors of 2 and 4 on the shear terms — and the
        conversion here follows whichever the concrete subclass is, so the two
        never get mixed up.
        """

        matrix = _tensor_to_voigt(self.tensor, compliance=self._compliance)
        matrix = np.ascontiguousarray(matrix)
        matrix.setflags(write=False)
        return matrix

    def rotate(self, rotation_matrix: ArrayLike) -> ElasticTensor:
        """Rotate the tensor by a proper rotation ``R`` (crystal -> new frame)."""

        rotation = np.asarray(rotation_matrix, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError("rotation_matrix must have shape (3, 3).")
        rotated = np.einsum(
            "ip,jq,kr,ls,pqrs->ijkl",
            rotation,
            rotation,
            rotation,
            rotation,
            self.tensor,
            optimize=True,
        )
        return type(self)(tensor=rotated, _compliance=self._compliance)


@dataclass(frozen=True, slots=True)
class StiffnessTensor(ElasticTensor):
    """Elastic stiffness tensor ``C`` (stress = C : strain), Voigt units of GPa."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "_compliance", False)
        ElasticTensor.__post_init__(self)

    @classmethod
    def from_voigt(cls, matrix: ArrayLike) -> StiffnessTensor:
        """Stiffness tensor from a ``6x6`` Voigt matrix.

        The matrix must be symmetric; this is checked, since an asymmetric
        stiffness matrix violates the existence of a strain-energy density and
        usually indicates a transcription error.
        """

        voigt = np.asarray(matrix, dtype=np.float64)
        if voigt.shape != (6, 6):
            raise ValueError("Voigt stiffness matrix must have shape (6, 6).")
        if not np.allclose(voigt, voigt.T, atol=1e-9):
            raise ValueError("Voigt stiffness matrix must be symmetric.")
        return cls(tensor=_voigt_to_tensor(voigt, compliance=False))

    @classmethod
    def cubic(cls, c11: float, c12: float, c44: float) -> StiffnessTensor:
        """Cubic stiffness tensor from the three independent constants.

        Cubic symmetry leaves only ``C11``, ``C12``, and ``C44`` independent. The
        Zener anisotropy ratio ``2*C44 / (C11 - C12)`` measures how far the
        crystal is from isotropy: it equals 1 for tungsten and about 3.8 for
        copper, which is why directional-modulus surfaces are near-spherical for
        one and strongly lobed for the other.

        Units are the caller's; they carry through to every derived modulus.
        """

        matrix = np.zeros((6, 6), dtype=np.float64)
        for index in range(3):
            matrix[index, index] = c11
            matrix[index + 3, index + 3] = c44
        for i in range(3):
            for j in range(3):
                if i != j:
                    matrix[i, j] = c12
        return cls.from_voigt(matrix)

    @classmethod
    def hexagonal(
        cls,
        c11: float,
        c12: float,
        c13: float,
        c33: float,
        c44: float,
    ) -> StiffnessTensor:
        """Hexagonal stiffness tensor from the five independent constants.

        Hexagonal (transversely isotropic) symmetry leaves ``C11``, ``C12``,
        ``C13``, ``C33``, and ``C44`` independent; ``C66 = (C11 - C12) / 2``
        follows from the symmetry and is applied here rather than being asked
        for.
        """

        c66 = 0.5 * (c11 - c12)
        matrix = np.array(
            [
                [c11, c12, c13, 0.0, 0.0, 0.0],
                [c12, c11, c13, 0.0, 0.0, 0.0],
                [c13, c13, c33, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, c44, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, c66],
            ],
            dtype=np.float64,
        )
        return cls.from_voigt(matrix)

    @classmethod
    def isotropic(cls, *, youngs_modulus: float, poisson_ratio: float) -> StiffnessTensor:
        """Isotropic stiffness tensor from Young's modulus and Poisson's ratio.

        Purpose
        -------
        The reference case with no directional dependence — useful as a baseline
        against which a real crystal's anisotropy is measured, and as a stand-in
        when single-crystal constants are unavailable.

        Parameters
        ----------
        youngs_modulus : float
            Young's modulus; sets the unit of every derived quantity.
        poisson_ratio : float
            Poisson's ratio. Thermodynamic stability requires it to lie in
            ``(-1, 0.5)``; values approaching 0.5 make the tensor
            near-incompressible and numerically stiff.
        """

        e = youngs_modulus
        nu = poisson_ratio
        factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu))
        c11 = factor * (1.0 - nu)
        c12 = factor * nu
        c44 = 0.5 * e / (1.0 + nu)
        return cls.cubic(c11, c12, c44)

    def compliance(self) -> ComplianceTensor:
        """The compliance tensor, obtained by inverting the Voigt stiffness matrix.

        Inverse of :meth:`ComplianceTensor.stiffness`. Directional moduli are
        naturally expressed through compliance, which is why the directional
        queries below delegate to it.
        """

        voigt_compliance = np.linalg.inv(self.voigt_matrix())
        return ComplianceTensor(tensor=_voigt_to_tensor(voigt_compliance, compliance=True))

    def youngs_modulus(self, direction: ArrayLike) -> np.ndarray | float:
        """Directional Young's modulus; see
        :meth:`ComplianceTensor.youngs_modulus`.
        """

        return self.compliance().youngs_modulus(direction)

    def linear_compressibility(self, direction: ArrayLike) -> np.ndarray | float:
        """Directional linear compressibility; see
        :meth:`ComplianceTensor.linear_compressibility`.
        """

        return self.compliance().linear_compressibility(direction)

    def shear_modulus(
        self, plane_normal: ArrayLike, shear_direction: ArrayLike
    ) -> np.ndarray | float:
        """Directional shear modulus; see
        :meth:`ComplianceTensor.shear_modulus`.
        """

        return self.compliance().shear_modulus(plane_normal, shear_direction)

    def poisson_ratio(
        self, direction: ArrayLike, transverse_direction: ArrayLike
    ) -> np.ndarray | float:
        """Directional Poisson's ratio; see
        :meth:`ComplianceTensor.poisson_ratio`.
        """

        return self.compliance().poisson_ratio(direction, transverse_direction)


def _unit_rows(direction: ArrayLike, name: str) -> tuple[np.ndarray, bool]:
    """Normalize a single ``(3,)`` vector or an ``(n, 3)`` batch to unit rows."""

    unit = np.asarray(direction, dtype=np.float64)
    scalar = unit.ndim == 1
    vectors = np.atleast_2d(unit)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    if np.any(norms == 0.0):
        raise ValueError(f"{name} must be a non-zero vector.")
    return vectors / norms, scalar


def _orthonormal_inplane_basis(normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return unit vectors (u, v) spanning the plane orthogonal to each unit normal."""

    seed = np.zeros_like(normals)
    seed[np.arange(normals.shape[0]), np.argmin(np.abs(normals), axis=1)] = 1.0
    u = np.cross(normals, seed)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v = np.cross(normals, u)
    return u, v


def _planar_quadratic_extrema(
    matrices: np.ndarray, normals: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Exact (min, max) of ``m^T Q m`` over unit ``m`` orthogonal to each normal.

    ``matrices`` has shape ``(n, 3, 3)``. Restricted to the plane, the quadratic
    form becomes a symmetric 2x2 eigenproblem in an orthonormal in-plane basis,
    so the extrema are its eigenvalues (closed form, no angular sweep).
    """

    u, v = _orthonormal_inplane_basis(normals)
    q_uu = np.einsum("nj,njl,nl->n", u, matrices, u, optimize=True)
    q_vv = np.einsum("nj,njl,nl->n", v, matrices, v, optimize=True)
    q_uv = 0.5 * (
        np.einsum("nj,njl,nl->n", u, matrices, v, optimize=True)
        + np.einsum("nj,njl,nl->n", v, matrices, u, optimize=True)
    )
    center = 0.5 * (q_uu + q_vv)
    radius = np.sqrt((0.5 * (q_uu - q_vv)) ** 2 + q_uv**2)
    return center - radius, center + radius


@dataclass(frozen=True, slots=True)
class ComplianceTensor(ElasticTensor):
    """Elastic compliance tensor ``S`` (strain = S : stress), Voigt units of 1/GPa."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "_compliance", True)
        ElasticTensor.__post_init__(self)

    @classmethod
    def from_voigt(cls, matrix: ArrayLike) -> ComplianceTensor:
        """Compliance tensor from a ``6x6`` Voigt compliance matrix.

        The compliance Voigt convention carries factors of 2 and 4 on the shear
        terms; the conversion applies them, so a compliance matrix must not be
        passed to :meth:`StiffnessTensor.from_voigt` or the shear components will
        be wrong by those factors.
        """

        voigt = np.asarray(matrix, dtype=np.float64)
        if voigt.shape != (6, 6):
            raise ValueError("Voigt compliance matrix must have shape (6, 6).")
        return cls(tensor=_voigt_to_tensor(voigt, compliance=True))

    def stiffness(self) -> StiffnessTensor:
        """The stiffness tensor, obtained by inverting the Voigt compliance matrix.
        """

        voigt_stiffness = np.linalg.inv(self.voigt_matrix())
        return StiffnessTensor(tensor=_voigt_to_tensor(voigt_stiffness, compliance=False))

    def youngs_modulus(self, direction: ArrayLike) -> np.ndarray | float:
        """Directional Young's modulus ``E(n) = 1 / (n_i n_j n_k n_l S_ijkl)``.

        Accepts a single ``(3,)`` direction (returns a float) or a batch of
        ``(n, 3)`` directions (returns an ``(n,)`` array), fully vectorised.
        """

        vectors, scalar = _unit_rows(direction, "direction")
        inverse_modulus = np.einsum(
            "ni,nj,nk,nl,ijkl->n", vectors, vectors, vectors, vectors, self.tensor, optimize=True
        )
        if np.any(inverse_modulus <= 0.0):
            raise ValueError("Non-physical compliance produced a non-positive modulus.")
        modulus = 1.0 / inverse_modulus
        return float(modulus[0]) if scalar else np.ascontiguousarray(modulus)

    def linear_compressibility(self, direction: ArrayLike) -> np.ndarray | float:
        """Linear compressibility ``beta(n) = n_i n_j S_ijkk`` under hydrostatic load.

        Accepts a single ``(3,)`` direction (returns a float) or a batch of
        ``(n, 3)`` directions (returns an ``(n,)`` array).
        """

        vectors, scalar = _unit_rows(direction, "direction")
        beta = np.einsum("ni,nj,ijkk->n", vectors, vectors, self.tensor, optimize=True)
        return float(beta[0]) if scalar else np.ascontiguousarray(beta)

    def shear_modulus(
        self, plane_normal: ArrayLike, shear_direction: ArrayLike
    ) -> np.ndarray | float:
        """Directional shear modulus ``G(n, m) = 1 / (4 n_i m_j n_k m_l S_ijkl)``.

        ``plane_normal`` is the shear-plane normal ``n`` and ``shear_direction``
        the in-plane shear direction ``m``; the two must be orthogonal. Accepts
        single ``(3,)`` vectors (returns a float) or matching ``(n, 3)`` batches
        (returns an ``(n,)`` array), fully vectorised.
        """

        normals, scalar_n = _unit_rows(plane_normal, "plane_normal")
        shears, scalar_m = _unit_rows(shear_direction, "shear_direction")
        if normals.shape != shears.shape:
            raise ValueError("plane_normal and shear_direction must have matching shapes.")
        if np.any(np.abs(np.einsum("ni,ni->n", normals, shears)) > 1e-8):
            raise ValueError(
                "shear_direction must be orthogonal to plane_normal (lie in the shear plane)."
            )
        inverse_modulus = 4.0 * np.einsum(
            "ni,nj,nk,nl,ijkl->n", normals, shears, normals, shears, self.tensor, optimize=True
        )
        if np.any(inverse_modulus <= 0.0):
            raise ValueError("Non-physical compliance produced a non-positive shear modulus.")
        modulus = 1.0 / inverse_modulus
        return float(modulus[0]) if scalar_n and scalar_m else np.ascontiguousarray(modulus)

    def poisson_ratio(
        self, direction: ArrayLike, transverse_direction: ArrayLike
    ) -> np.ndarray | float:
        """Directional Poisson's ratio for uniaxial stress along ``n``.

        ``nu(n, m) = - (m_i m_j n_k n_l S_ijkl) / (n_i n_j n_k n_l S_ijkl)`` is
        the negative ratio of the transverse strain along ``m`` to the axial
        strain along ``n``; the two directions must be orthogonal. Accepts
        single ``(3,)`` vectors (returns a float) or matching ``(n, 3)`` batches
        (returns an ``(n,)`` array), fully vectorised.
        """

        normals, scalar_n = _unit_rows(direction, "direction")
        transverse, scalar_m = _unit_rows(transverse_direction, "transverse_direction")
        if normals.shape != transverse.shape:
            raise ValueError("direction and transverse_direction must have matching shapes.")
        if np.any(np.abs(np.einsum("ni,ni->n", normals, transverse)) > 1e-8):
            raise ValueError("transverse_direction must be orthogonal to direction.")
        axial = np.einsum(
            "ni,nj,nk,nl,ijkl->n", normals, normals, normals, normals, self.tensor, optimize=True
        )
        if np.any(axial <= 0.0):
            raise ValueError("Non-physical compliance produced a non-positive axial strain.")
        coupling = np.einsum(
            "ni,nj,nk,nl,ijkl->n",
            transverse,
            transverse,
            normals,
            normals,
            self.tensor,
            optimize=True,
        )
        ratio = -coupling / axial
        return float(ratio[0]) if scalar_n and scalar_m else np.ascontiguousarray(ratio)


def _normalized_weights(count: int, weights: ArrayLike | None) -> np.ndarray:
    if weights is None:
        return np.full(count, 1.0 / count, dtype=np.float64)
    values = np.asarray(weights, dtype=np.float64)
    if values.shape != (count,):
        raise ValueError("weights must provide one value per orientation.")
    if np.any(values < 0.0):
        raise ValueError("weights must be non-negative.")
    total = float(values.sum())
    if np.isclose(total, 0.0):
        raise ValueError("weights must not sum to zero.")
    return values / total


def homogenize_elastic(
    stiffness: StiffnessTensor,
    orientations: OrientationSet,
    *,
    weights: ArrayLike | None = None,
    scheme: str = "hill",
) -> StiffnessTensor:
    """Orientation-weighted polycrystal elastic average of a single crystal.

    Rotates the single-crystal stiffness into the sample frame for every
    orientation and averages under the requested scheme:

    - ``"voigt"``: arithmetic mean of the stiffness tensors (uniform strain),
    - ``"reuss"``: inverse of the mean compliance (uniform stress),
    - ``"hill"``: the Voigt-Reuss-Hill average ``(C_voigt + C_reuss) / 2``.

    Returns the aggregate `StiffnessTensor` in the sample frame.
    """

    if scheme not in {"voigt", "reuss", "hill"}:
        raise ValueError("scheme must be one of 'voigt', 'reuss', or 'hill'.")
    matrices = orientations.as_matrices()
    weight_values = _normalized_weights(matrices.shape[0], weights)
    crystal_stiffness = stiffness.tensor
    crystal_compliance = stiffness.compliance().tensor
    rotated_c = np.einsum(
        "nip,njq,nkr,nls,pqrs->nijkl",
        matrices,
        matrices,
        matrices,
        matrices,
        crystal_stiffness,
        optimize=True,
    )
    voigt_tensor = np.einsum("n,nijkl->ijkl", weight_values, rotated_c, optimize=True)
    if scheme == "voigt":
        return StiffnessTensor(tensor=voigt_tensor)
    rotated_s = np.einsum(
        "nip,njq,nkr,nls,pqrs->nijkl",
        matrices,
        matrices,
        matrices,
        matrices,
        crystal_compliance,
        optimize=True,
    )
    mean_compliance = np.einsum("n,nijkl->ijkl", weight_values, rotated_s, optimize=True)
    reuss_stiffness = ComplianceTensor(tensor=mean_compliance).stiffness()
    if scheme == "reuss":
        return reuss_stiffness
    hill_voigt = 0.5 * (voigt_tensor + reuss_stiffness.tensor)
    return StiffnessTensor(tensor=hill_voigt)


@dataclass(frozen=True, slots=True)
class DirectionalModulusSurface:
    """A directional elastic property sampled on a spherical (theta, phi) grid."""

    theta: np.ndarray
    phi: np.ndarray
    values: np.ndarray
    property_name: str = "youngs_modulus"

    def __post_init__(self) -> None:
        for name in ("theta", "phi", "values"):
            array = np.ascontiguousarray(np.asarray(getattr(self, name), dtype=np.float64))
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if self.values.shape != (self.theta.shape[0], self.phi.shape[0]):
            raise ValueError("values must have shape (len(theta), len(phi)).")

    @property
    def minimum(self) -> float:
        """Smallest directional property value on the sampled surface.

        For Young's modulus of a cubic metal this is the soft direction —
        ``<100>`` when the Zener ratio exceeds 1, ``<111>`` when it is below.
        """

        return float(np.min(self.values))

    @property
    def maximum(self) -> float:
        """Largest directional property value on the sampled surface.
        """

        return float(np.max(self.values))

    @property
    def anisotropy_ratio(self) -> float:
        """Ratio of maximum to minimum directional value (1 for isotropy)."""

        return self.maximum / self.minimum

    def cartesian_surface(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (X, Y, Z) grids with radius scaled by the property value."""

        polar, azimuth = np.meshgrid(self.theta, self.phi, indexing="ij")
        radius = self.values
        x = radius * np.sin(polar) * np.cos(azimuth)
        y = radius * np.sin(polar) * np.sin(azimuth)
        z = radius * np.cos(polar)
        return x, y, z


def _directional_property_surface(
    property_fn: Callable[[np.ndarray], np.ndarray | float],
    property_name: str,
    *,
    n_theta: int,
    n_phi: int,
) -> DirectionalModulusSurface:
    if n_theta < 2 or n_phi < 2:
        raise ValueError("n_theta and n_phi must each be at least 2.")
    theta = np.linspace(0.0, np.pi, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)
    polar, azimuth = np.meshgrid(theta, phi, indexing="ij")
    directions = np.stack(
        [
            np.sin(polar) * np.cos(azimuth),
            np.sin(polar) * np.sin(azimuth),
            np.cos(polar),
        ],
        axis=-1,
    ).reshape(-1, 3)
    values = np.asarray(property_fn(directions), dtype=np.float64)
    return DirectionalModulusSurface(
        theta, phi, values.reshape(n_theta, n_phi), property_name=property_name
    )


def youngs_modulus_surface(
    tensor: StiffnessTensor | ComplianceTensor,
    *,
    n_theta: int = 90,
    n_phi: int = 180,
) -> DirectionalModulusSurface:
    """Sample the directional Young's modulus over the unit sphere.

    Returns a `DirectionalModulusSurface` with the modulus on a structured
    ``(n_theta, n_phi)`` grid, ready for 3D plotting or anisotropy analysis.
    """

    return _directional_property_surface(
        tensor.youngs_modulus, "youngs_modulus", n_theta=n_theta, n_phi=n_phi
    )


def linear_compressibility_surface(
    tensor: StiffnessTensor | ComplianceTensor,
    *,
    n_theta: int = 90,
    n_phi: int = 180,
) -> DirectionalModulusSurface:
    """Sample the directional linear compressibility over the unit sphere."""

    return _directional_property_surface(
        tensor.linear_compressibility, "linear_compressibility", n_theta=n_theta, n_phi=n_phi
    )


def shear_modulus_surface(
    tensor: StiffnessTensor | ComplianceTensor,
    *,
    mode: str = "min",
    n_theta: int = 90,
    n_phi: int = 180,
) -> DirectionalModulusSurface:
    """Sample the extremal directional shear modulus over the unit sphere.

    For each shear-plane normal ``n``, the shear modulus still depends on the
    in-plane shear direction ``m``; this returns, per normal, the exact minimum
    (``mode="min"``) or maximum (``mode="max"``) of ``G(n, m)`` over all
    in-plane directions. Because ``1/G`` is a quadratic form in ``m``, the
    extremes are the eigenvalues of the form projected onto the plane -- no
    angular sweep or sampling error.
    """

    if mode not in {"min", "max"}:
        raise ValueError("mode must be 'min' or 'max'.")
    compliance = tensor.compliance() if isinstance(tensor, StiffnessTensor) else tensor

    def extremal_shear(directions: np.ndarray) -> np.ndarray:
        normals = np.asarray(directions, dtype=np.float64)
        forms = 4.0 * np.einsum(
            "ni,nk,ijkl->njl", normals, normals, compliance.tensor, optimize=True
        )
        low, high = _planar_quadratic_extrema(forms, normals)
        inverse_modulus = high if mode == "min" else low
        if np.any(inverse_modulus <= 0.0):
            raise ValueError("Non-physical compliance produced a non-positive shear modulus.")
        return 1.0 / inverse_modulus

    return _directional_property_surface(
        extremal_shear, f"shear_modulus_{mode}", n_theta=n_theta, n_phi=n_phi
    )


def poisson_ratio_surface(
    tensor: StiffnessTensor | ComplianceTensor,
    *,
    mode: str = "min",
    n_theta: int = 90,
    n_phi: int = 180,
) -> DirectionalModulusSurface:
    """Sample the extremal directional Poisson's ratio over the unit sphere.

    For each loading direction ``n``, Poisson's ratio still depends on the
    transverse direction ``m``; this returns, per direction, the exact minimum
    (``mode="min"``) or maximum (``mode="max"``) of ``nu(n, m)`` over all
    transverse directions. Because ``nu`` is a quadratic form in ``m``, the
    extremes are the eigenvalues of the form projected onto the transverse
    plane -- no angular sweep. Negative minima flag auxetic response.
    """

    if mode not in {"min", "max"}:
        raise ValueError("mode must be 'min' or 'max'.")
    compliance = tensor.compliance() if isinstance(tensor, StiffnessTensor) else tensor

    def extremal_poisson(directions: np.ndarray) -> np.ndarray:
        normals = np.asarray(directions, dtype=np.float64)
        axial = np.einsum(
            "ni,nj,nk,nl,ijkl->n",
            normals,
            normals,
            normals,
            normals,
            compliance.tensor,
            optimize=True,
        )
        if np.any(axial <= 0.0):
            raise ValueError("Non-physical compliance produced a non-positive axial strain.")
        forms = -np.einsum(
            "nk,nl,ijkl->nij", normals, normals, compliance.tensor, optimize=True
        ) / axial[:, None, None]
        low, high = _planar_quadratic_extrema(forms, normals)
        return low if mode == "min" else high

    return _directional_property_surface(
        extremal_poisson, f"poisson_ratio_{mode}", n_theta=n_theta, n_phi=n_phi
    )


__all__ = [
    "ComplianceTensor",
    "DirectionalModulusSurface",
    "ElasticTensor",
    "StiffnessTensor",
    "homogenize_elastic",
    "linear_compressibility_surface",
    "poisson_ratio_surface",
    "shear_modulus_surface",
    "youngs_modulus_surface",
]
