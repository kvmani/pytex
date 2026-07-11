"""Full crystallographic point-group model for all 32 point groups.

Point groups are expressed in the crystal Cartesian frame with the principal
axis along +Z and the secondary axis along +X (hexagonal setting for trigonal
and hexagonal groups). Operators include improper operations (mirrors,
inversion, rotoinversions) as orthogonal matrices with determinant -1, unlike
`SymmetrySpec`, which carries only the proper rotation subgroup used for
orientation reduction.

Axis-placement rule: every group is oriented so that its proper rotations
coincide operator-for-operator with the canonical proper point group of its
rotation subgroup (secondary 2-fold axes along +X). For -6m2 this corresponds
to the -62m axis setting; both symbols are accepted.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, freeze_array, normalize_vector

_OPERATOR_TOLERANCE = 1e-8

_CRYSTAL_SYSTEMS = (
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
)


def _rotation(axis: ArrayLike, angle_deg: float) -> np.ndarray:
    unit_axis = normalize_vector(axis)
    angle_rad = np.deg2rad(angle_deg)
    x, y, z = unit_axis
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    one_minus_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_minus_c, x * y * one_minus_c - z * s, x * z * one_minus_c + y * s],
            [y * x * one_minus_c + z * s, c + y * y * one_minus_c, y * z * one_minus_c - x * s],
            [z * x * one_minus_c - y * s, z * y * one_minus_c + x * s, c + z * z * one_minus_c],
        ],
        dtype=np.float64,
    )


def _mirror(normal: ArrayLike) -> np.ndarray:
    unit_normal = normalize_vector(normal)
    return np.eye(3, dtype=np.float64) - 2.0 * np.outer(unit_normal, unit_normal)


def _rotoinversion(axis: ArrayLike, angle_deg: float) -> np.ndarray:
    return -_rotation(axis, angle_deg)


_INVERSION = -np.eye(3, dtype=np.float64)

_X = (1.0, 0.0, 0.0)
_Y = (0.0, 1.0, 0.0)
_Z = (0.0, 0.0, 1.0)
_DIAGONAL_111 = (1.0, 1.0, 1.0)


def _is_orthogonal_operator(matrix: np.ndarray) -> bool:
    if matrix.shape != (3, 3):
        return False
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=_OPERATOR_TOLERANCE):
        return False
    determinant = float(np.linalg.det(matrix))
    return bool(np.isclose(abs(determinant), 1.0, atol=_OPERATOR_TOLERANCE))


def _matrix_key(matrix: np.ndarray) -> tuple[float, ...]:
    rounded = np.round(matrix, decimals=8) + 0.0
    return tuple(float(value) for value in rounded.ravel())


def _group_from_orthogonal_generators(generators: tuple[np.ndarray, ...]) -> np.ndarray:
    identity = np.eye(3, dtype=np.float64)
    known: dict[tuple[float, ...], np.ndarray] = {_matrix_key(identity): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            for candidate in (current @ generator, generator @ current):
                if not _is_orthogonal_operator(candidate):
                    raise ValueError(
                        "Point-group generators produced a non-orthogonal operator."
                    )
                key = _matrix_key(candidate)
                if key not in known:
                    known[key] = candidate
                    frontier.append(candidate)
    stacked = np.stack(list(known.values()), axis=0)
    return freeze_array(np.ascontiguousarray(stacked))


@dataclass(frozen=True)
class _PointGroupDefinition:
    hermann_mauguin: str
    schoenflies: str
    crystal_system: str
    laue_class: str
    proper_subgroup: str
    order: int
    generators: tuple[np.ndarray, ...]


def _definitions() -> dict[str, _PointGroupDefinition]:
    def define(
        hermann_mauguin: str,
        schoenflies: str,
        crystal_system: str,
        laue_class: str,
        proper_subgroup: str,
        order: int,
        generators: tuple[np.ndarray, ...],
    ) -> _PointGroupDefinition:
        return _PointGroupDefinition(
            hermann_mauguin=hermann_mauguin,
            schoenflies=schoenflies,
            crystal_system=crystal_system,
            laue_class=laue_class,
            proper_subgroup=proper_subgroup,
            order=order,
            generators=generators,
        )

    return {
        "1": define("1", "C1", "triclinic", "-1", "1", 1, ()),
        "-1": define("-1", "Ci", "triclinic", "-1", "1", 2, (_INVERSION,)),
        "2": define("2", "C2", "monoclinic", "2/m", "2", 2, (_rotation(_Z, 180.0),)),
        "m": define("m", "Cs", "monoclinic", "2/m", "1", 2, (_mirror(_Z),)),
        "2/m": define(
            "2/m", "C2h", "monoclinic", "2/m", "2", 4, (_rotation(_Z, 180.0), _INVERSION)
        ),
        "222": define(
            "222",
            "D2",
            "orthorhombic",
            "mmm",
            "222",
            4,
            (_rotation(_X, 180.0), _rotation(_Y, 180.0)),
        ),
        "mm2": define(
            "mm2", "C2v", "orthorhombic", "mmm", "2", 4, (_mirror(_X), _mirror(_Y))
        ),
        "mmm": define(
            "mmm",
            "D2h",
            "orthorhombic",
            "mmm",
            "222",
            8,
            (_rotation(_X, 180.0), _rotation(_Y, 180.0), _INVERSION),
        ),
        "4": define("4", "C4", "tetragonal", "4/m", "4", 4, (_rotation(_Z, 90.0),)),
        "-4": define("-4", "S4", "tetragonal", "4/m", "2", 4, (_rotoinversion(_Z, 90.0),)),
        "4/m": define(
            "4/m", "C4h", "tetragonal", "4/m", "4", 8, (_rotation(_Z, 90.0), _INVERSION)
        ),
        "422": define(
            "422",
            "D4",
            "tetragonal",
            "4/mmm",
            "422",
            8,
            (_rotation(_Z, 90.0), _rotation(_X, 180.0)),
        ),
        "4mm": define(
            "4mm", "C4v", "tetragonal", "4/mmm", "4", 8, (_rotation(_Z, 90.0), _mirror(_X))
        ),
        "-42m": define(
            "-42m",
            "D2d",
            "tetragonal",
            "4/mmm",
            "222",
            8,
            (_rotoinversion(_Z, 90.0), _rotation(_X, 180.0)),
        ),
        "4/mmm": define(
            "4/mmm",
            "D4h",
            "tetragonal",
            "4/mmm",
            "422",
            16,
            (_rotation(_Z, 90.0), _rotation(_X, 180.0), _INVERSION),
        ),
        "3": define("3", "C3", "trigonal", "-3", "3", 3, (_rotation(_Z, 120.0),)),
        "-3": define("-3", "C3i", "trigonal", "-3", "3", 6, (_rotation(_Z, 120.0), _INVERSION)),
        "32": define(
            "32",
            "D3",
            "trigonal",
            "-3m",
            "32",
            6,
            (_rotation(_Z, 120.0), _rotation(_X, 180.0)),
        ),
        "3m": define(
            "3m", "C3v", "trigonal", "-3m", "3", 6, (_rotation(_Z, 120.0), _mirror(_X))
        ),
        "-3m": define(
            "-3m",
            "D3d",
            "trigonal",
            "-3m",
            "32",
            12,
            (_rotation(_Z, 120.0), _rotation(_X, 180.0), _INVERSION),
        ),
        "6": define("6", "C6", "hexagonal", "6/m", "6", 6, (_rotation(_Z, 60.0),)),
        "-6": define("-6", "C3h", "hexagonal", "6/m", "3", 6, (_rotoinversion(_Z, 60.0),)),
        "6/m": define(
            "6/m", "C6h", "hexagonal", "6/m", "6", 12, (_rotation(_Z, 60.0), _INVERSION)
        ),
        "622": define(
            "622",
            "D6",
            "hexagonal",
            "6/mmm",
            "622",
            12,
            (_rotation(_Z, 60.0), _rotation(_X, 180.0)),
        ),
        "6mm": define(
            "6mm", "C6v", "hexagonal", "6/mmm", "6", 12, (_rotation(_Z, 60.0), _mirror(_X))
        ),
        "-6m2": define(
            "-6m2",
            "D3h",
            "hexagonal",
            "6/mmm",
            "32",
            12,
            (_rotoinversion(_Z, 60.0), _rotation(_X, 180.0)),
        ),
        "6/mmm": define(
            "6/mmm",
            "D6h",
            "hexagonal",
            "6/mmm",
            "622",
            24,
            (_rotation(_Z, 60.0), _rotation(_X, 180.0), _INVERSION),
        ),
        "23": define(
            "23",
            "T",
            "cubic",
            "m-3",
            "23",
            12,
            (_rotation(_X, 180.0), _rotation(_DIAGONAL_111, 120.0)),
        ),
        "m-3": define(
            "m-3",
            "Th",
            "cubic",
            "m-3",
            "23",
            24,
            (_rotation(_X, 180.0), _rotation(_DIAGONAL_111, 120.0), _INVERSION),
        ),
        "432": define(
            "432",
            "O",
            "cubic",
            "m-3m",
            "432",
            24,
            (_rotation(_Z, 90.0), _rotation(_DIAGONAL_111, 120.0)),
        ),
        "-43m": define(
            "-43m",
            "Td",
            "cubic",
            "m-3m",
            "23",
            24,
            (_rotoinversion(_Z, 90.0), _rotation(_DIAGONAL_111, 120.0)),
        ),
        "m-3m": define(
            "m-3m",
            "Oh",
            "cubic",
            "m-3m",
            "432",
            48,
            (_rotation(_Z, 90.0), _rotation(_DIAGONAL_111, 120.0), _INVERSION),
        ),
    }


@cache
def _definition_table() -> dict[str, _PointGroupDefinition]:
    return _definitions()


@cache
def _symbol_aliases() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for symbol, definition in _definition_table().items():
        aliases[symbol.lower()] = symbol
        aliases[definition.schoenflies.lower()] = symbol
    aliases["m3"] = "m-3"
    aliases["m3m"] = "m-3m"
    aliases["43m"] = "-43m"
    aliases["s6"] = "-3"
    aliases["3/m"] = "-6"
    aliases["-62m"] = "-6m2"
    aliases["62m"] = "-6m2"
    return aliases


def normalize_point_group_symbol(symbol: str) -> str:
    normalized = symbol.replace(" ", "").lower()
    resolved = _symbol_aliases().get(normalized)
    if resolved is None:
        supported = ", ".join(sorted(_definition_table()))
        raise ValueError(
            f"Unsupported point-group symbol '{symbol}'. Supported Hermann-Mauguin "
            f"symbols: {supported} (Schoenflies names are also accepted)."
        )
    return resolved


def all_point_group_symbols() -> tuple[str, ...]:
    return tuple(_definition_table())


def proper_subgroup_symbol_for(symbol: str) -> str:
    return _definition_table()[normalize_point_group_symbol(symbol)].proper_subgroup


def laue_class_symbol_for(symbol: str) -> str:
    return _definition_table()[normalize_point_group_symbol(symbol)].laue_class


def laue_class_symbols() -> tuple[str, ...]:
    return tuple(dict.fromkeys(entry.laue_class for entry in _definition_table().values()))


@cache
def _operators_for_symbol(symbol: str) -> np.ndarray:
    definition = _definition_table()[symbol]
    operators = _group_from_orthogonal_generators(definition.generators)
    if operators.shape[0] != definition.order:
        raise ValueError(
            f"Point group '{symbol}' generated {operators.shape[0]} operators, "
            f"expected {definition.order}."
        )
    return operators


def _canonical_direction_sign(vector: np.ndarray) -> np.ndarray:
    for component in vector:
        if not np.isclose(component, 0.0, atol=_OPERATOR_TOLERANCE):
            return vector if component > 0.0 else -vector
    return vector


@dataclass(frozen=True, slots=True)
class PointGroup:
    hermann_mauguin: str
    schoenflies: str
    crystal_system: str
    operators: np.ndarray

    def __post_init__(self) -> None:
        operators = as_float_array(self.operators, shape=(None, 3, 3))
        for operator in operators:
            if not _is_orthogonal_operator(np.asarray(operator)):
                raise ValueError(
                    "All point-group operators must be orthogonal with determinant +1 or -1."
                )
        if self.crystal_system not in _CRYSTAL_SYSTEMS:
            supported = ", ".join(_CRYSTAL_SYSTEMS)
            raise ValueError(
                f"Unsupported crystal system '{self.crystal_system}'. "
                f"Supported systems: {supported}."
            )
        object.__setattr__(self, "operators", operators)

    @classmethod
    def from_symbol(cls, symbol: str) -> PointGroup:
        canonical = normalize_point_group_symbol(symbol)
        definition = _definition_table()[canonical]
        return cls(
            hermann_mauguin=definition.hermann_mauguin,
            schoenflies=definition.schoenflies,
            crystal_system=definition.crystal_system,
            operators=_operators_for_symbol(canonical),
        )

    @property
    def order(self) -> int:
        return int(self.operators.shape[0])

    @property
    def determinants(self) -> np.ndarray:
        return freeze_array(np.ascontiguousarray(np.linalg.det(self.operators)))

    @property
    def rotations(self) -> np.ndarray:
        proper = self.operators[self.determinants > 0.0]
        return freeze_array(np.ascontiguousarray(proper))

    @property
    def improper_operators(self) -> np.ndarray:
        improper = self.operators[self.determinants < 0.0]
        return freeze_array(np.ascontiguousarray(improper))

    @property
    def is_proper(self) -> bool:
        return bool(np.all(self.determinants > 0.0))

    @property
    def is_centrosymmetric(self) -> bool:
        return bool(
            any(
                np.allclose(operator, _INVERSION, atol=_OPERATOR_TOLERANCE)
                for operator in self.operators
            )
        )

    @property
    def is_laue(self) -> bool:
        return self.hermann_mauguin == self.laue_class_symbol

    @property
    def laue_class_symbol(self) -> str:
        return _definition_table()[self.hermann_mauguin].laue_class

    @property
    def proper_subgroup_symbol(self) -> str:
        return _definition_table()[self.hermann_mauguin].proper_subgroup

    def laue_class(self) -> PointGroup:
        return PointGroup.from_symbol(self.laue_class_symbol)

    def proper_subgroup(self) -> PointGroup:
        return PointGroup.from_symbol(self.proper_subgroup_symbol)

    def mirror_normals(self) -> np.ndarray:
        normals: dict[tuple[float, ...], np.ndarray] = {}
        for operator in self.improper_operators:
            if not np.isclose(float(np.trace(operator)), 1.0, atol=_OPERATOR_TOLERANCE):
                continue
            eigenvalues, eigenvectors = np.linalg.eigh(operator)
            normal = _canonical_direction_sign(
                normalize_vector(eigenvectors[:, int(np.argmin(eigenvalues))])
            )
            normals[tuple(np.round(normal, decimals=8))] = normal
        if not normals:
            return freeze_array(np.zeros((0, 3), dtype=np.float64))
        stacked = np.stack(list(normals.values()), axis=0)
        return freeze_array(np.ascontiguousarray(stacked))

    def equivalent_directions(self, vector: ArrayLike, *, antipodal: bool = False) -> np.ndarray:
        unit = normalize_vector(vector)
        candidates = np.einsum("oij,j->oi", self.operators, unit, optimize=True)
        if antipodal:
            candidates = np.concatenate([candidates, -candidates], axis=0)
        unique: dict[tuple[float, ...], np.ndarray] = {}
        for candidate in candidates:
            unique[tuple(np.round(candidate, decimals=8))] = candidate
        stacked = np.stack(list(unique.values()), axis=0)
        return freeze_array(np.ascontiguousarray(stacked))

    def to_symmetry_spec(self, **kwargs: object) -> object:
        from pytex.core.symmetry import SymmetrySpec

        return SymmetrySpec.from_point_group(self.hermann_mauguin, **kwargs)  # type: ignore[arg-type]
