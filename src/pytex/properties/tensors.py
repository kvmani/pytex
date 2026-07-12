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
        voigt = np.asarray(matrix, dtype=np.float64)
        if voigt.shape != (6, 6):
            raise ValueError("Voigt stiffness matrix must have shape (6, 6).")
        if not np.allclose(voigt, voigt.T, atol=1e-9):
            raise ValueError("Voigt stiffness matrix must be symmetric.")
        return cls(tensor=_voigt_to_tensor(voigt, compliance=False))

    @classmethod
    def cubic(cls, c11: float, c12: float, c44: float) -> StiffnessTensor:
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
        e = youngs_modulus
        nu = poisson_ratio
        factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu))
        c11 = factor * (1.0 - nu)
        c12 = factor * nu
        c44 = 0.5 * e / (1.0 + nu)
        return cls.cubic(c11, c12, c44)

    def compliance(self) -> ComplianceTensor:
        voigt_compliance = np.linalg.inv(self.voigt_matrix())
        return ComplianceTensor(tensor=_voigt_to_tensor(voigt_compliance, compliance=True))

    def youngs_modulus(self, direction: ArrayLike) -> np.ndarray | float:
        return self.compliance().youngs_modulus(direction)

    def linear_compressibility(self, direction: ArrayLike) -> np.ndarray | float:
        return self.compliance().linear_compressibility(direction)


@dataclass(frozen=True, slots=True)
class ComplianceTensor(ElasticTensor):
    """Elastic compliance tensor ``S`` (strain = S : stress), Voigt units of 1/GPa."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "_compliance", True)
        ElasticTensor.__post_init__(self)

    @classmethod
    def from_voigt(cls, matrix: ArrayLike) -> ComplianceTensor:
        voigt = np.asarray(matrix, dtype=np.float64)
        if voigt.shape != (6, 6):
            raise ValueError("Voigt compliance matrix must have shape (6, 6).")
        return cls(tensor=_voigt_to_tensor(voigt, compliance=True))

    def stiffness(self) -> StiffnessTensor:
        voigt_stiffness = np.linalg.inv(self.voigt_matrix())
        return StiffnessTensor(tensor=_voigt_to_tensor(voigt_stiffness, compliance=False))

    def youngs_modulus(self, direction: ArrayLike) -> np.ndarray | float:
        """Directional Young's modulus ``E(n) = 1 / (n_i n_j n_k n_l S_ijkl)``.

        Accepts a single ``(3,)`` direction (returns a float) or a batch of
        ``(n, 3)`` directions (returns an ``(n,)`` array), fully vectorised.
        """

        unit = np.asarray(direction, dtype=np.float64)
        scalar = unit.ndim == 1
        vectors = np.atleast_2d(unit)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        if np.any(norms == 0.0):
            raise ValueError("direction must be a non-zero vector.")
        vectors = vectors / norms
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

        unit = np.asarray(direction, dtype=np.float64)
        scalar = unit.ndim == 1
        vectors = np.atleast_2d(unit)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        if np.any(norms == 0.0):
            raise ValueError("direction must be a non-zero vector.")
        vectors = vectors / norms
        beta = np.einsum("ni,nj,ijkk->n", vectors, vectors, self.tensor, optimize=True)
        return float(beta[0]) if scalar else np.ascontiguousarray(beta)


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
        return float(np.min(self.values))

    @property
    def maximum(self) -> float:
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


__all__ = [
    "ComplianceTensor",
    "DirectionalModulusSurface",
    "ElasticTensor",
    "StiffnessTensor",
    "homogenize_elastic",
    "linear_compressibility_surface",
    "youngs_modulus_surface",
]
