"""Elastic stiffness and compliance tensors with directional properties.

`StiffnessTensor` and `ComplianceTensor` store the full rank-4 elastic tensor in
crystal-frame Cartesian coordinates and round-trip to the 6x6 Voigt matrix.
Stiffness uses no Voigt factors; compliance carries the standard 1/2 and 1/4
factors on shear terms so that ``S = C^-1`` in Voigt form is consistent with the
tensor contraction ``1/E(n) = n_i n_j n_k n_l S_ijkl``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

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

    def youngs_modulus(self, direction: ArrayLike) -> float:
        return self.compliance().youngs_modulus(direction)


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

    def youngs_modulus(self, direction: ArrayLike) -> float:
        """Directional Young's modulus ``E(n) = 1 / (n_i n_j n_k n_l S_ijkl)``."""

        unit = np.asarray(direction, dtype=np.float64)
        norm = float(np.linalg.norm(unit))
        if norm == 0.0:
            raise ValueError("direction must be a non-zero vector.")
        unit = unit / norm
        inverse_modulus = float(
            np.einsum("i,j,k,l,ijkl->", unit, unit, unit, unit, self.tensor, optimize=True)
        )
        if inverse_modulus <= 0.0:
            raise ValueError("Non-physical compliance produced a non-positive modulus.")
        return 1.0 / inverse_modulus


__all__ = [
    "ComplianceTensor",
    "ElasticTensor",
    "StiffnessTensor",
]
