from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np

from pytex.core._arrays import FloatArray, IntArray, as_float_array, as_int_array
from pytex.core.lattice import CrystalPlane, Phase, ZoneAxis

_ROUND_DECIMALS = 12


def _direct_basis_matrix(phase: Phase) -> FloatArray:
    return as_float_array(phase.lattice.direct_basis().matrix, shape=(3, 3))


def _reciprocal_basis_matrix(phase: Phase) -> FloatArray:
    return as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))


def _require_matching_phases(
    left: Phase,
    right: Phase,
    *,
    left_name: str,
    right_name: str,
) -> None:
    if left != right:
        raise ValueError(f"{left_name}.phase must match {right_name}.phase.")


def _as_index_rows(values: Any, *, name: str, columns: int) -> IntArray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (columns,):
            raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
        array = array[None, :]
    elif array.ndim == 2:
        if array.shape[1] != columns:
            raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
    else:
        raise ValueError(f"{name} must have shape ({columns},) or (n, {columns}).")
    return as_int_array(array, shape=(None, columns))


def _validate_nonzero_rows(indices: IntArray, *, name: str) -> None:
    if np.any(~np.any(indices != 0, axis=1)):
        raise ValueError(f"{name} must not contain the zero triplet.")


def _rowwise_gcd(values: IntArray) -> IntArray:
    if values.shape[0] == 0:
        return as_int_array(np.ones(0, dtype=np.int64), shape=(0,))
    abs_values = np.abs(values)
    divisors = np.gcd.reduce(abs_values, axis=1)
    divisors[divisors == 0] = 1
    return as_int_array(divisors, shape=(values.shape[0],))


def _first_nonzero_signs(values: IntArray) -> IntArray:
    non_zero = values != 0
    first_index = np.argmax(non_zero, axis=1)
    rows = np.arange(values.shape[0], dtype=np.int64)
    selected = values[rows, first_index]
    signs = np.sign(selected).astype(np.int64, copy=False)
    signs[signs == 0] = 1
    return as_int_array(signs, shape=(values.shape[0],))


def reduce_indices(values: Any, *, name: str = "indices") -> IntArray:
    """Divide index rows by their greatest common divisor.

    ``(2, 2, 0)`` becomes ``(1, 1, 0)``: the same plane or direction written
    in lowest terms. Accepts a ``(3,)`` triple or an ``(n, 3)`` array and
    always returns the ``(n, 3)`` form. The zero triplet is rejected, since
    it denotes neither a plane nor a direction.
    """

    rows = _as_index_rows(values, name=name, columns=3)
    _validate_nonzero_rows(rows, name=name)
    divisors = _rowwise_gcd(rows)
    reduced = rows // divisors[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 3))


def canonicalize_sign(values: Any, *, name: str = "indices") -> IntArray:
    """Fix the overall sign of index rows by their first nonzero component.

    Maps ``(-1, 1, 0)`` and ``(1, -1, 0)`` to the same representative, which
    is what makes antipodal families comparable. Use it when two index sets
    must be tested for equality up to inversion; use
    :func:`canonicalize_family_indices` when the rows should also be reduced
    to lowest terms first.
    """

    rows = _as_index_rows(values, name=name, columns=3)
    _validate_nonzero_rows(rows, name=name)
    signs = _first_nonzero_signs(rows)
    canonical = rows * signs[:, None]
    return as_int_array(canonical, shape=(rows.shape[0], 3))


def canonicalize_family_indices(
    values: Any,
    *,
    antipodal: bool,
    name: str = "indices",
) -> IntArray:
    """Reduce index rows to lowest terms and, optionally, canonicalize sign.

    The single entry point for putting Miller indices into a comparable
    form. ``antipodal=True`` is right for planes and for direction families
    quoted without a sense; ``antipodal=False`` keeps ``[uvw]`` distinct
    from ``[-u-v-w]``, which matters for slip directions and for
    Burgers-vector bookkeeping.
    """

    reduced = reduce_indices(values, name=name)
    if antipodal:
        return canonicalize_sign(reduced, name=name)
    return reduced


def antipodal_keys(values: Any, *, name: str = "indices") -> IntArray:
    """Comparable keys for index rows under inversion.

    Shorthand for ``canonicalize_family_indices(values, antipodal=True)``.
    Two rows describing the same plane, or the same direction up to sense,
    produce identical keys, so the result can be used for grouping,
    deduplication, and set membership.
    """

    return canonicalize_family_indices(values, antipodal=True, name=name)


def plane_hkl_to_hkil_array(hkl: Any) -> IntArray:
    """Three-index planes ``(hkl)`` to hexagonal four-index ``(hkil)``.

    Inserts the redundant third index ``i = -(h + k)`` required by the
    Miller-Bravais convention, which exists so that the three symmetry-
    equivalent prismatic planes of a hexagonal crystal receive indices that
    are permutations of one another. Accepts ``(3,)`` or ``(n, 3)`` and
    returns ``(n, 4)``.
    """

    rows = _as_index_rows(hkl, name="hkl", columns=3)
    h = rows[:, 0]
    k = rows[:, 1]
    ell = rows[:, 2]
    converted = np.column_stack([h, k, -(h + k), ell])
    return as_int_array(converted, shape=(rows.shape[0], 4))


def plane_hkil_to_hkl_array(hkil: Any) -> IntArray:
    """Hexagonal four-index planes ``(hkil)`` to three-index ``(hkl)``.

    Drops the redundant ``i``. The constraint ``i = -(h + k)`` is *checked*,
    not assumed, so a mistyped four-index plane raises instead of silently
    producing wrong three-index values. Accepts ``(4,)`` or ``(n, 4)`` and
    returns ``(n, 3)``.
    """

    array = np.ascontiguousarray(np.asarray(hkil, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (4,):
            raise ValueError("hkil must have shape (4,) or (n, 4).")
        array = array[None, :]
    elif array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("hkil must have shape (4,) or (n, 4).")
    rows = as_int_array(array, shape=(None, 4))
    if np.any(rows[:, 2] != -(rows[:, 0] + rows[:, 1])):
        raise ValueError("Hexagonal four-index planes must satisfy i = -(h + k).")
    converted = np.column_stack([rows[:, 0], rows[:, 1], rows[:, 3]])
    return as_int_array(converted, shape=(rows.shape[0], 3))


def direction_uvw_to_uvtw_array(uvw: Any) -> IntArray:
    """Three-index directions ``[uvw]`` to hexagonal four-index ``[UVTW]``.

    Applies the standard transformation ``U = (2u - v)/3``,
    ``V = (2v - u)/3``, ``T = -(U + V)``, ``W = w``. Where the division by
    three would not give integers, the row is scaled by three instead of
    rounded, and the result is then reduced by its greatest common divisor,
    so the returned indices are always exact integers describing the same
    direction. Accepts ``(3,)`` or ``(n, 3)`` and returns ``(n, 4)``.
    """

    rows = _as_index_rows(uvw, name="uvw", columns=3)
    two_u_minus_v = 2 * rows[:, 0] - rows[:, 1]
    two_v_minus_u = 2 * rows[:, 1] - rows[:, 0]
    minus_u_minus_v = -(rows[:, 0] + rows[:, 1])
    divisible_by_three = (
        (two_u_minus_v % 3 == 0)
        & (two_v_minus_u % 3 == 0)
        & (minus_u_minus_v % 3 == 0)
    )
    reduced_form = np.column_stack(
        [two_u_minus_v // 3, two_v_minus_u // 3, minus_u_minus_v // 3, rows[:, 2]]
    )
    expanded_form = np.column_stack(
        [two_u_minus_v, two_v_minus_u, minus_u_minus_v, 3 * rows[:, 2]]
    )
    scaled = np.where(divisible_by_three[:, None], reduced_form, expanded_form)
    reduced = scaled // _rowwise_gcd(scaled)[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 4))


def direction_uvtw_to_uvw_array(uvtw: Any) -> IntArray:
    """Hexagonal four-index directions ``[UVTW]`` to three-index ``[uvw]``.

    Applies ``u = 2U + V``, ``v = 2V + U``, ``w = W`` and reduces to lowest
    terms. The redundancy constraint ``U + V + T = 0`` is checked, not
    assumed. Accepts ``(4,)`` or ``(n, 4)`` and returns ``(n, 3)``.
    """

    array = np.ascontiguousarray(np.asarray(uvtw, dtype=np.int64))
    if array.ndim == 1:
        if array.shape != (4,):
            raise ValueError("UVTW must have shape (4,) or (n, 4).")
        array = array[None, :]
    elif array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("UVTW must have shape (4,) or (n, 4).")
    rows = as_int_array(array, shape=(None, 4))
    if np.any(np.sum(rows[:, :3], axis=1) != 0):
        raise ValueError("Hexagonal four-index directions must satisfy U + V + T = 0.")
    converted = np.column_stack(
        [2 * rows[:, 0] + rows[:, 1], 2 * rows[:, 1] + rows[:, 0], rows[:, 3]]
    )
    reduced = converted // _rowwise_gcd(converted)[:, None]
    return as_int_array(reduced, shape=(rows.shape[0], 3))


def zone_law_value_hkil_uvtw(hkil: Any, uvtw: Any) -> IntArray:
    """The zone-law value ``hu + kv + lw`` for four-index hexagonal inputs.

    Purpose
    -------
    A direction lies in a plane exactly when the zone-law value is zero.
    This four-index form converts both arguments to three-index form first,
    because the zone law is *not* the naive dot product of the four-index
    sets. Working in four-index form directly is a classic error; this
    function exists so callers need not repeat the conversion themselves.

    Parameters
    ----------
    hkil, uvtw : ArrayLike
        ``(4,)`` or ``(n, 4)`` index rows. One row broadcasts against many.

    Returns
    -------
    IntArray
        ``(n,)`` integer zone-law values; zero means the direction lies in
        the plane.
    """

    planes = plane_hkil_to_hkl_array(hkil)
    directions = direction_uvtw_to_uvw_array(uvtw)
    plane_rows, direction_rows = _broadcast_rows(
        planes,
        directions,
        left_name="hkil",
        right_name="UVTW",
    )
    values = np.einsum("ni,ni->n", plane_rows, direction_rows, optimize=True)
    return as_int_array(values.astype(np.int64), shape=(values.shape[0],))


def _family_unique_rows(values: IntArray) -> tuple[IntArray, IntArray]:
    if values.shape[0] == 0:
        empty_rows = np.empty((0, values.shape[1]), dtype=np.int64)
        empty_inverse = np.empty(0, dtype=np.int64)
        return (
            as_int_array(empty_rows, shape=(0, values.shape[1])),
            as_int_array(empty_inverse, shape=(0,)),
        )
    unique_rows, inverse = np.unique(values, axis=0, return_inverse=True)
    return (
        as_int_array(unique_rows, shape=(unique_rows.shape[0], unique_rows.shape[1])),
        as_int_array(inverse, shape=(inverse.shape[0],)),
    )


def _freeze_bool_array(values: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype=bool)
    array.setflags(write=False)
    return array


def _broadcast_rows(
    left: np.ndarray,
    right: np.ndarray,
    *,
    left_name: str,
    right_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    if left.shape[0] == right.shape[0]:
        return left, right
    if left.shape[0] == 1:
        return np.broadcast_to(left, right.shape), right
    if right.shape[0] == 1:
        return left, np.broadcast_to(right, left.shape)
    raise ValueError(
        f"{left_name} and {right_name} must have matching lengths or one row for broadcasting."
    )


def _cartesian_from_direct_indices(indices: IntArray, phase: Phase) -> FloatArray:
    basis = _direct_basis_matrix(phase)
    cartesian = np.asarray(indices, dtype=np.float64) @ basis.T
    return as_float_array(cartesian, shape=(indices.shape[0], 3))


def _cartesian_from_reciprocal_indices(indices: IntArray, phase: Phase) -> FloatArray:
    basis = _reciprocal_basis_matrix(phase)
    cartesian = np.asarray(indices, dtype=np.float64) @ basis.T
    return as_float_array(cartesian, shape=(indices.shape[0], 3))


def _unit_vectors(vectors: FloatArray, *, name: str) -> FloatArray:
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(np.isclose(norms, 0.0)):
        raise ValueError(f"{name} contains a zero cartesian vector.")
    unit = vectors / norms[:, None]
    return as_float_array(unit, shape=(vectors.shape[0], 3))


def _structured_sort_rows(values: IntArray) -> IntArray:
    if values.size == 0:
        return values
    structured = values.view(dtype=[("c0", np.int64), ("c1", np.int64), ("c2", np.int64)]).reshape(
        values.shape[0], values.shape[1]
    )
    sorted_structured = np.sort(structured, axis=1)
    sorted_values = sorted_structured.view(np.int64).reshape(values.shape[0], values.shape[1], 3)
    return as_int_array(sorted_values, shape=(values.shape[0], values.shape[1], 3))


def _pack_unique_rows(values: IntArray) -> tuple[IntArray, np.ndarray]:
    if values.shape[0] == 0:
        empty_values = np.empty((0, 0, 3), dtype=np.int64)
        empty_mask = np.empty((0, 0), dtype=bool)
        return as_int_array(empty_values, shape=(0, 0, 3)), _freeze_bool_array(empty_mask)
    sorted_values = _structured_sort_rows(values)
    unique_mask = np.ones(sorted_values.shape[:2], dtype=bool)
    if sorted_values.shape[1] > 1:
        unique_mask[:, 1:] = np.any(
            sorted_values[:, 1:, :] != sorted_values[:, :-1, :],
            axis=2,
        )
    counts = np.sum(unique_mask, axis=1)
    max_count = int(np.max(counts)) if counts.size else 0
    packed = np.zeros((sorted_values.shape[0], max_count, 3), dtype=np.int64)
    packed_mask = np.zeros((sorted_values.shape[0], max_count), dtype=bool)
    if max_count:
        row_ids = np.broadcast_to(
            np.arange(sorted_values.shape[0], dtype=np.int64)[:, None],
            unique_mask.shape,
        )
        destination = np.cumsum(unique_mask, axis=1) - 1
        packed[row_ids[unique_mask], destination[unique_mask]] = sorted_values[unique_mask]
        packed_mask[row_ids[unique_mask], destination[unique_mask]] = True
    return (
        as_int_array(packed, shape=(packed.shape[0], packed.shape[1], 3)),
        _freeze_bool_array(packed_mask),
    )


def _recover_integer_indices(
    cartesian_vectors: np.ndarray,
    *,
    phase: Phase,
    reciprocal: bool,
    tol: float,
) -> IntArray:
    basis = _reciprocal_basis_matrix(phase) if reciprocal else _direct_basis_matrix(phase)
    inverse = np.linalg.inv(np.asarray(basis, dtype=np.float64))
    coordinates = np.einsum("...j,ij->...i", cartesian_vectors, inverse, optimize=True)
    rounded = np.rint(coordinates)
    if not np.allclose(coordinates, rounded, atol=tol, rtol=0.0):
        raise ValueError(
            "Symmetry operators produced coordinates that could not be recovered as "
            "integer Miller indices."
        )
    return as_int_array(rounded.astype(np.int64), shape=coordinates.shape)


def _symmetry_equivalent_indices(
    indices: IntArray,
    *,
    phase: Phase,
    reciprocal: bool,
    unique: bool,
    antipodal: bool,
    tol: float,
) -> tuple[IntArray, np.ndarray]:
    if indices.shape[0] == 0:
        empty_values = np.empty((0, 0, 3), dtype=np.int64)
        empty_mask = np.empty((0, 0), dtype=bool)
        return as_int_array(empty_values, shape=(0, 0, 3)), _freeze_bool_array(empty_mask)
    cartesian = (
        _cartesian_from_reciprocal_indices(indices, phase)
        if reciprocal
        else _cartesian_from_direct_indices(indices, phase)
    )
    transformed = np.einsum(
        "oij,nj->noi",
        phase.symmetry.operators,
        np.asarray(cartesian, dtype=np.float64),
        optimize=True,
    )
    recovered = _recover_integer_indices(
        transformed,
        phase=phase,
        reciprocal=reciprocal,
        tol=tol,
    )
    flattened = recovered.reshape(-1, 3)
    canonical = canonicalize_family_indices(
        flattened,
        antipodal=antipodal,
        name="symmetry_equivalent_indices",
    ).reshape(recovered.shape)
    if unique:
        return _pack_unique_rows(as_int_array(canonical, shape=recovered.shape))
    mask = np.ones(canonical.shape[:2], dtype=bool)
    return (
        as_int_array(canonical, shape=canonical.shape),
        _freeze_bool_array(mask),
    )


def _pairwise_angles_from_unit_vectors(
    left: FloatArray,
    right: FloatArray,
    *,
    left_name: str,
    right_name: str,
    antipodal: bool,
) -> FloatArray:
    broadcast_left, broadcast_right = _broadcast_rows(
        left,
        right,
        left_name=left_name,
        right_name=right_name,
    )
    dots = np.einsum("ni,ni->n", broadcast_left, broadcast_right, optimize=True)
    if antipodal:
        dots = np.abs(dots)
    return as_float_array(np.arccos(np.clip(dots, -1.0, 1.0)), shape=(dots.shape[0],))


def _pairwise_matrix_angles_from_unit_vectors(
    left: FloatArray,
    right: FloatArray,
    *,
    antipodal: bool,
) -> FloatArray:
    dots = np.asarray(left, dtype=np.float64) @ np.asarray(right, dtype=np.float64).T
    if antipodal:
        dots = np.abs(dots)
    return as_float_array(
        np.arccos(np.clip(dots, -1.0, 1.0)),
        shape=(left.shape[0], right.shape[0]),
    )


def _projection_vectors(
    directions: FloatArray,
    plane_normals: FloatArray,
) -> tuple[FloatArray, np.ndarray]:
    broadcast_directions, broadcast_normals = _broadcast_rows(
        directions,
        plane_normals,
        left_name="directions",
        right_name="planes",
    )
    normal_components = np.einsum(
        "ni,ni->n",
        broadcast_directions,
        broadcast_normals,
        optimize=True,
    )
    projected = broadcast_directions - normal_components[:, None] * broadcast_normals
    degenerate_mask = np.isclose(np.linalg.norm(projected, axis=1), 0.0)
    projected[degenerate_mask] = 0.0
    return (
        as_float_array(projected, shape=(projected.shape[0], 3)),
        _freeze_bool_array(degenerate_mask),
    )


@dataclass(frozen=True, slots=True)
class MillerPlane:
    """A crystal plane ``(hkl)`` with the full index algebra attached.

    Purpose
    -------
    Where :class:`~pytex.core.lattice.CrystalPlane` holds the geometry, this
    adds the *index* operations: symmetry families ``{hkl}``, reduction to
    lowest terms, antipodal keys for comparison, and four-index conversion.
    It is the type to reach for when reasoning about which planes are
    equivalent, rather than about where one plane points.

    Attributes
    ----------
    indices : np.ndarray
        Integer ``(h, k, l)``; the zero triplet is rejected.
    phase : Phase
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("MillerPlane indices must not be the zero triplet.")

    @classmethod
    def from_hkl(cls, indices: Any, *, phase: Phase) -> MillerPlane:
        """Plane from three-index ``(hkl)`` Miller indices on a phase.
        """

        return cls(indices=indices, phase=phase)

    @classmethod
    def from_hkil(cls, indices: Any, *, phase: Phase) -> MillerPlane:
        """Plane from hexagonal four-index ``(hkil)`` Miller-Bravais indices.

        The redundancy constraint ``i = -(h + k)`` is checked on conversion.
        """

        return cls(indices=plane_hkil_to_hkl_array(indices)[0], phase=phase)

    @classmethod
    def from_crystal_plane(cls, plane: CrystalPlane) -> MillerPlane:
        """Adopt the indices and phase of a :class:`~pytex.core.lattice.CrystalPlane`.

        ``MillerPlane`` is the index-algebra view (families, symmetry orbits,
        reductions); ``CrystalPlane`` is the lattice-geometry view. This
        converts between them without losing phase meaning.
        """

        return cls(indices=plane.miller.indices, phase=plane.phase)

    @property
    def hkil(self) -> IntArray:
        """The four-index ``(hkil)`` form of this plane.

        Meaningful for hexagonal and trigonal phases, where it makes the
        symmetry of the prismatic and pyramidal families visible; computed
        formally for any phase.
        """

        return as_int_array(plane_hkl_to_hkil_array(self.indices)[0], shape=(4,))

    @property
    def reduced_indices(self) -> IntArray:
        """The indices divided by their greatest common divisor.

        ``(2 2 0)`` reports ``(1 1 0)``. Note that the reduced form denotes the
        same *plane* but a different *reflection*: ``(220)`` and ``(110)`` have
        different d-spacings and are distinct diffraction events.
        """

        return as_int_array(reduce_indices(self.indices)[0], shape=(3,))

    @property
    def antipodal_key(self) -> IntArray:
        """A sign-canonical key identifying this plane's family membership.

        Two planes that differ only by inversion give the same key, so keys can
        be compared, grouped, and deduplicated.
        """

        return as_int_array(antipodal_keys(self.indices)[0], shape=(3,))

    @property
    def reciprocal_vector_cartesian(self) -> FloatArray:
        """The reciprocal-lattice vector ``g = h a* + k b* + l c*``, in Cartesian
        crystal-frame coordinates and in inverse angstroms.

        Not normalized: its magnitude is ``1 / d``, which is what the
        diffraction layer needs. For the direction alone use
        :attr:`normal_cartesian`.
        """

        return as_float_array(
            _cartesian_from_reciprocal_indices(self.indices[None, :], self.phase)[0],
            shape=(3,),
        )

    @property
    def normal_cartesian(self) -> FloatArray:
        """Unit normal of the plane in the Cartesian crystal frame.

        The normal is obtained through the *reciprocal* basis, which is the only
        correct route in a non-cubic lattice: outside the cubic system the
        direction ``[hkl]`` is not parallel to the normal of ``(hkl)``.
        """

        return as_float_array(
            _unit_vectors(self.reciprocal_vector_cartesian[None, :], name="MillerPlane")[0],
            shape=(3,),
        )

    @property
    def d_spacing_angstrom(self) -> float:
        """Interplanar spacing ``d`` in angstroms.

        Computed as the reciprocal of the reciprocal-vector magnitude, so it is
        correct for every crystal system without a per-system formula. This is
        the quantity Bragg's law consumes.
        """

        magnitude = float(np.linalg.norm(self.reciprocal_vector_cartesian))
        if np.isclose(magnitude, 0.0):
            raise ValueError("MillerPlane reciprocal vector magnitude must be non-zero.")
        return 1.0 / magnitude

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        """Symmetry-equivalent indices of every plane, as a padded array.

        Purpose
        -------
        Enumerate the symmetry family {hkl} of each plane: the indices that
        the crystal's point group makes crystallographically indistinguishable.

        Parameters
        ----------
        unique : bool
            Collapse repeated members of the orbit (default ``True``). Repeats
            occur whenever an operator maps a plane onto itself.
        antipodal : bool
            Treat ``h`` and ``-h`` as the same family member. Only ``True`` is accepted for
            planes: a plane and its opposite normal are the same plane, so
            ``antipodal=False`` would be a meaningless request and raises.
        tol : float
            Tolerance on recovering integer indices after applying the Cartesian
            symmetry operators. Exceeding it raises rather than silently
            rounding to wrong indices.

        Returns
        -------
        tuple of (IntArray, np.ndarray)
            ``(n, m, 3)`` indices and an ``(n, m)`` boolean validity mask. Rows
            are padded to the longest family in the batch, so the mask — not the
            array shape — says how many members each row really has.
        """

        if not antipodal:
            raise ValueError("Miller planes are always treated with antipodal equivalence.")
        return _symmetry_equivalent_indices(
            self.indices[None, :],
            phase=self.phase,
            reciprocal=True,
            unique=unique,
            antipodal=True,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> MillerPlaneSet:
        """The symmetry family ``{hkl}`` of this plane, as a plane set.

        The typed form of :meth:`symmetry_equivalent_indices`: the padding mask
        is already applied, so every row of the returned set is a real family
        member. This is the set a pole figure of ``{hkl}`` plots.
        """

        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return MillerPlaneSet(indices=equivalent_indices[0, mask[0]], phase=self.phase)

    def to_crystal_plane(self) -> CrystalPlane:
        """The lattice-geometry view of this plane.

        See :meth:`from_crystal_plane` for why the two views exist.
        """

        return CrystalPlane.from_miller_bravais(self.hkil, phase=self.phase)

    def to_miller_bravais(self) -> MillerBravaisPlane:
        """The four-index :class:`MillerBravaisPlane` view of this plane.
        """

        return MillerBravaisPlane(indices=self.hkil, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerBravaisPlane:
    """A hexagonal plane in four-index ``(hkil)`` form.

    Purpose
    -------
    The four-index notation used for hexagonal and trigonal crystals, where
    the redundant ``i = -(h + k)`` makes the three symmetry-equivalent
    prismatic planes receive indices that are permutations of one another —
    something the three-index form cannot express.

    The redundancy constraint is checked on construction. Note that the zone
    law is *not* the naive dot product of four-index sets, which is why
    :meth:`zone_law_value` converts to three-index form first.

    Attributes
    ----------
    indices : np.ndarray
        Integer ``(h, k, i, l)`` satisfying ``i = -(h + k)``.
    phase : Phase
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        rows = plane_hkil_to_hkl_array(self.indices)
        hkil = plane_hkl_to_hkil_array(rows)[0]
        if not np.any(rows[0]):
            raise ValueError("MillerBravaisPlane indices must not be the zero quartet.")
        object.__setattr__(self, "indices", as_int_array(hkil, shape=(4,)))

    @classmethod
    def from_hkil(cls, indices: Any, *, phase: Phase) -> MillerBravaisPlane:
        """Four-index hexagonal plane from ``(hkil)``, checking ``i = -(h + k)``.
        """

        return cls(indices=indices, phase=phase)

    @classmethod
    def from_hkl(cls, indices: Any, *, phase: Phase) -> MillerBravaisPlane:
        """Four-index hexagonal plane from three-index ``(hkl)``.
        """

        return cls(indices=plane_hkl_to_hkil_array(indices)[0], phase=phase)

    @property
    def hkl(self) -> IntArray:
        """The three-index ``(hkl)`` form, with the redundant ``i`` dropped.
        """

        return as_int_array(plane_hkil_to_hkl_array(self.indices)[0], shape=(3,))

    @property
    def reduced_indices(self) -> IntArray:
        """The four indices divided by their greatest common divisor.
        """

        return as_int_array(plane_hkl_to_hkil_array(reduce_indices(self.hkl))[0], shape=(4,))

    def zone_law_value(self, direction: MillerBravaisDirection) -> int:
        """``hu + kv + lw`` for this plane and a four-index direction.

        Zero means the direction lies in the plane. Both operands are converted
        to three-index form first, because the zone law is not the naive dot
        product of four-index sets.
        """

        if direction.phase != self.phase:
            raise ValueError("direction.phase must match MillerBravaisPlane.phase.")
        return int(zone_law_value_hkil_uvtw(self.indices, direction.indices)[0])

    def contains_direction(self, direction: MillerBravaisDirection) -> bool:
        """Whether a four-index direction lies in this plane (zone-law value zero).
        """

        return self.zone_law_value(direction) == 0

    def to_miller_plane(self) -> MillerPlane:
        """The three-index :class:`MillerPlane` view of this plane.
        """

        return MillerPlane.from_hkl(self.hkl, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerDirection:
    """A crystal direction ``[uvw]`` with the full index algebra attached.

    Purpose
    -------
    The direction counterpart of :class:`MillerPlane`: symmetry families
    ``<uvw>``, reduction, antipodal keys, four-index conversion, and
    conversion to a zone axis.

    Unlike planes, directions have a genuine sense, so antipodal treatment is
    optional here rather than mandatory — ``[uvw]`` and ``[-u-v-w]`` are the
    same *family member* but not the same Burgers vector.

    Attributes
    ----------
    indices : np.ndarray
        Integer ``(u, v, w)``; the zero triplet is rejected.
    phase : Phase
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(3,)))
        if not np.any(self.indices):
            raise ValueError("MillerDirection indices must not be the zero triplet.")

    @classmethod
    def from_uvw(cls, indices: Any, *, phase: Phase) -> MillerDirection:
        """Direction from three-index ``[uvw]`` indices on a phase.
        """

        return cls(indices=indices, phase=phase)

    @classmethod
    def from_UVTW(cls, indices: Any, *, phase: Phase) -> MillerDirection:  # noqa: N802
        """Direction from hexagonal four-index ``[UVTW]`` indices.

        The redundancy constraint ``U + V + T = 0`` is checked on conversion.
        """

        return cls(indices=direction_uvtw_to_uvw_array(indices)[0], phase=phase)

    @classmethod
    def from_zone_axis(cls, zone_axis: ZoneAxis) -> MillerDirection:
        """Adopt the indices and phase of a :class:`~pytex.core.lattice.ZoneAxis`.

        A zone axis *is* a lattice direction; this gives the same vector the
        index-algebra surface (families, symmetry orbits, reductions).
        """

        return cls(indices=zone_axis.indices, phase=zone_axis.phase)

    @property
    def UVTW(self) -> IntArray:  # noqa: N802
        """The hexagonal four-index ``[UVTW]`` form of this direction.
        """

        return as_int_array(direction_uvw_to_uvtw_array(self.indices)[0], shape=(4,))

    @property
    def reduced_indices(self) -> IntArray:
        """The indices divided by their greatest common divisor.
        """

        return as_int_array(reduce_indices(self.indices)[0], shape=(3,))

    @property
    def antipodal_key(self) -> IntArray:
        """A sign-canonical key identifying this direction up to sense.

        Use it when ``[uvw]`` and ``[-u-v-w]`` should compare equal — for a
        direction family quoted without a sense. Do not use it where the sense
        carries physics, as with a Burgers vector or a shear direction.
        """

        return as_int_array(antipodal_keys(self.indices)[0], shape=(3,))

    @property
    def direct_vector_cartesian(self) -> FloatArray:
        """The lattice vector ``u a + v b + w c`` in Cartesian crystal-frame
        coordinates and in angstroms.

        Not normalized: its length is the repeat distance along the direction,
        which is what Burgers-vector magnitudes and structure-factor phases
        need.
        """

        return as_float_array(
            _cartesian_from_direct_indices(self.indices[None, :], self.phase)[0],
            shape=(3,),
        )

    @property
    def unit_vector_cartesian(self) -> FloatArray:
        """Unit vector along the direction, in the Cartesian crystal frame.

        Obtained through the *direct* basis — the correct route for directions,
        just as plane normals must go through the reciprocal basis.
        """

        return as_float_array(
            _unit_vectors(self.direct_vector_cartesian[None, :], name="MillerDirection")[0],
            shape=(3,),
        )

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        """Symmetry-equivalent indices of every direction, as a padded array.

        Purpose
        -------
        Enumerate the symmetry family <uvw> of each direction: the indices that
        the crystal's point group makes crystallographically indistinguishable.

        Parameters
        ----------
        unique : bool
            Collapse repeated members of the orbit (default ``True``). Repeats
            occur whenever an operator maps a direction onto itself.
        antipodal : bool
            Treat ``h`` and ``-h`` as the same family member. Directions have a genuine sense, so
            ``antipodal=False`` is meaningful and keeps ``[uvw]`` and
            ``[-u-v-w]`` distinct.
        tol : float
            Tolerance on recovering integer indices after applying the Cartesian
            symmetry operators. Exceeding it raises rather than silently
            rounding to wrong indices.

        Returns
        -------
        tuple of (IntArray, np.ndarray)
            ``(n, m, 3)`` indices and an ``(n, m)`` boolean validity mask. Rows
            are padded to the longest family in the batch, so the mask — not the
            array shape — says how many members each row really has.
        """

        return _symmetry_equivalent_indices(
            self.indices[None, :],
            phase=self.phase,
            reciprocal=False,
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> MillerDirectionSet:
        """The symmetry family ``<uvw>`` of this direction, as a direction set.

        The typed form of :meth:`symmetry_equivalent_indices`, with the padding
        mask already applied. This is the set an inverse pole figure of
        ``<uvw>`` plots, and the set slip-system enumeration expands.
        """

        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return MillerDirectionSet(indices=equivalent_indices[0, mask[0]], phase=self.phase)

    def to_zone_axis(self) -> ZoneAxis:
        """The zone-axis view of this direction, for diffraction work.
        """

        return ZoneAxis(indices=self.indices, phase=self.phase)

    def to_miller_bravais(self) -> MillerBravaisDirection:
        """The four-index :class:`MillerBravaisDirection` view of this direction.
        """

        return MillerBravaisDirection(indices=self.UVTW, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerBravaisDirection:
    """A hexagonal direction in four-index ``[UVTW]`` form.

    Purpose
    -------
    The four-index direction notation for hexagonal and trigonal crystals,
    which makes the ``<11-20>`` family's symmetry visible in the indices
    themselves. The constraint ``U + V + T = 0`` is checked on construction.

    Attributes
    ----------
    indices : np.ndarray
        Integer ``(U, V, T, W)`` satisfying ``U + V + T = 0``.
    phase : Phase
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        rows = direction_uvtw_to_uvw_array(self.indices)
        uvtw = direction_uvw_to_uvtw_array(rows)[0]
        if not np.any(rows[0]):
            raise ValueError("MillerBravaisDirection indices must not be the zero quartet.")
        object.__setattr__(self, "indices", as_int_array(uvtw, shape=(4,)))

    @classmethod
    def from_UVTW(cls, indices: Any, *, phase: Phase) -> MillerBravaisDirection:  # noqa: N802
        """Four-index hexagonal direction from ``[UVTW]``, checking ``U + V + T = 0``.
        """

        return cls(indices=indices, phase=phase)

    @classmethod
    def from_uvw(cls, indices: Any, *, phase: Phase) -> MillerBravaisDirection:
        """Four-index hexagonal direction from three-index ``[uvw]``.
        """

        return cls(indices=direction_uvw_to_uvtw_array(indices)[0], phase=phase)

    @property
    def uvw(self) -> IntArray:
        """The three-index ``[uvw]`` form of this direction.
        """

        return as_int_array(direction_uvtw_to_uvw_array(self.indices)[0], shape=(3,))

    @property
    def reduced_indices(self) -> IntArray:
        """The four indices divided by their greatest common divisor.
        """

        return as_int_array(direction_uvw_to_uvtw_array(reduce_indices(self.uvw))[0], shape=(4,))

    def zone_law_value(self, plane: MillerBravaisPlane) -> int:
        """``hu + kv + lw`` for this direction and a four-index plane.

        Zero means the direction lies in the plane. Both operands are converted
        to three-index form first.
        """

        if plane.phase != self.phase:
            raise ValueError("plane.phase must match MillerBravaisDirection.phase.")
        return int(zone_law_value_hkil_uvtw(plane.indices, self.indices)[0])

    def lies_in_plane(self, plane: MillerBravaisPlane) -> bool:
        """Whether this direction lies in a four-index plane (zone-law value zero).
        """

        return self.zone_law_value(plane) == 0

    def to_miller_direction(self) -> MillerDirection:
        """The three-index :class:`MillerDirection` view of this direction.
        """

        return MillerDirection.from_uvw(self.uvw, phase=self.phase)


@dataclass(frozen=True, slots=True)
class MillerPlaneSet:
    """A batch of crystal planes on one phase, with vectorized index algebra.

    Purpose
    -------
    The array form of :class:`MillerPlane`. Family expansion, d-spacings,
    normals, deduplication, and the full pairwise interplanar-angle table are
    all computed as array operations — the form a diffraction-pattern
    indexing or a pole-figure family expansion needs.

    Attributes
    ----------
    indices : np.ndarray
        ``(n, 3)`` integer indices; no row may be the zero triplet.
    phase : Phase
        Shared by every plane in the set.
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(None, 3)))
        _validate_nonzero_rows(cast(IntArray, self.indices), name="MillerPlaneSet.indices")

    @classmethod
    def from_hkl(cls, indices: Any, *, phase: Phase) -> MillerPlaneSet:
        """Plane set from an ``(n, 3)`` array of ``(hkl)`` indices on one phase.
        """

        return cls(indices=_as_index_rows(indices, name="hkl", columns=3), phase=phase)

    @classmethod
    def from_hkil(cls, indices: Any, *, phase: Phase) -> MillerPlaneSet:
        """Plane set from an ``(n, 4)`` array of ``(hkil)`` indices on one phase.
        """

        return cls(indices=plane_hkil_to_hkl_array(indices), phase=phase)

    def to_hkil(self) -> IntArray:
        """``(n, 4)`` hexagonal four-index form of the whole set.
        """

        return plane_hkl_to_hkil_array(self.indices)

    def reduce_indices(self) -> IntArray:
        """``(n, 3)`` indices divided row-wise by their greatest common divisor.
        """

        return reduce_indices(self.indices)

    def canonical_indices(self) -> IntArray:
        """``(n, 3)`` reduced, sign-canonical indices, comparable up to inversion.

        Two rows describing the same plane produce identical output rows, which
        is what :meth:`unique` groups on.
        """

        return antipodal_keys(self.indices)

    def unique(self) -> tuple[MillerPlaneSet, IntArray]:
        """The distinct planes in the set, up to reduction and inversion.

        Returns
        -------
        tuple of (MillerPlaneSet, IntArray)
            The deduplicated set, and an ``(n,)`` inverse-index array mapping
            every original row to its position in the deduplicated set — the
            same contract as ``numpy.unique(..., return_inverse=True)``, so
            per-row quantities can be scattered back onto the original ordering.
        """

        unique_indices, inverse = _family_unique_rows(self.canonical_indices())
        return MillerPlaneSet(indices=unique_indices, phase=self.phase), inverse

    def reciprocal_vectors_cartesian(self) -> FloatArray:
        """``(n, 3)`` reciprocal-lattice vectors in the Cartesian crystal frame,
        in inverse angstroms. Magnitudes are ``1 / d``.
        """

        return _cartesian_from_reciprocal_indices(cast(IntArray, self.indices), self.phase)

    def normals_cartesian(self) -> FloatArray:
        """``(n, 3)`` unit plane normals in the Cartesian crystal frame.

        Computed through the reciprocal basis, the only correct route outside
        the cubic system.
        """

        return _unit_vectors(self.reciprocal_vectors_cartesian(), name="MillerPlaneSet")

    def d_spacings_angstrom(self) -> FloatArray:
        """``(n,)`` interplanar spacings in angstroms.

        The vectorized form of :attr:`MillerPlane.d_spacing_angstrom`, and the
        input to a powder-pattern peak-position calculation.
        """

        magnitudes = np.linalg.norm(self.reciprocal_vectors_cartesian(), axis=1)
        if np.any(np.isclose(magnitudes, 0.0)):
            raise ValueError("MillerPlaneSet reciprocal vectors must be non-zero.")
        return as_float_array(1.0 / magnitudes, shape=(self.indices.shape[0],))

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        """Symmetry-equivalent indices of every plane, as a padded array.

        Purpose
        -------
        Enumerate the symmetry family {hkl} of each plane: the indices that
        the crystal's point group makes crystallographically indistinguishable.

        Parameters
        ----------
        unique : bool
            Collapse repeated members of the orbit (default ``True``). Repeats
            occur whenever an operator maps a plane onto itself.
        antipodal : bool
            Treat ``h`` and ``-h`` as the same family member. Only ``True`` is accepted for
            planes: a plane and its opposite normal are the same plane, so
            ``antipodal=False`` would be a meaningless request and raises.
        tol : float
            Tolerance on recovering integer indices after applying the Cartesian
            symmetry operators. Exceeding it raises rather than silently
            rounding to wrong indices.

        Returns
        -------
        tuple of (IntArray, np.ndarray)
            ``(n, m, 3)`` indices and an ``(n, m)`` boolean validity mask. Rows
            are padded to the longest family in the batch, so the mask — not the
            array shape — says how many members each row really has.
        """

        if not antipodal:
            raise ValueError("Miller planes are always treated with antipodal equivalence.")
        return _symmetry_equivalent_indices(
            cast(IntArray, self.indices),
            phase=self.phase,
            reciprocal=True,
            unique=unique,
            antipodal=True,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[MillerPlaneSet, ...]:
        """One plane set per input plane, holding that plane's ``{hkl}`` family.

        The typed form of :meth:`symmetry_equivalent_indices`, with the padding
        mask applied per row, so families of different sizes are returned as
        separate correctly sized sets rather than one padded array.
        """

        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return tuple(
            MillerPlaneSet(indices=equivalent_indices[row, mask[row]], phase=self.phase)
            for row in range(equivalent_indices.shape[0])
        )

    def to_crystal_planes(self) -> tuple[CrystalPlane, ...]:
        """The lattice-geometry view of every plane in the set.
        """

        return tuple(
            MillerPlane(indices=row, phase=self.phase).to_crystal_plane()
            for row in self.indices
        )

    def angle_matrix_rad(self, other: MillerPlaneSet | None = None) -> FloatArray:
        """Pairwise interplanar angles, in radians.

        Purpose
        -------
        The full angle table between two plane sets — the quantity used to index
        a diffraction pattern by matching measured inter-spot angles against
        computed ones, and to check that indexed poles are mutually consistent.

        Parameters
        ----------
        other : MillerPlaneSet, optional
            The second set; defaults to this set, giving the self-angle matrix
            with a zero diagonal. Must be on the same phase.

        Returns
        -------
        FloatArray
            ``(len(self), len(other))`` angles in radians. Antipodal
            equivalence is applied — a normal and its opposite describe the same
            plane — so every entry lies in ``[0, pi/2]``.
        """

        target = self if other is None else other
        _require_matching_phases(
            self.phase,
            target.phase,
            left_name="MillerPlaneSet",
            right_name="other",
        )
        return _pairwise_matrix_angles_from_unit_vectors(
            self.normals_cartesian(),
            target.normals_cartesian(),
            antipodal=True,
        )


@dataclass(frozen=True, slots=True)
class MillerDirectionSet:
    """A batch of crystal directions on one phase, with vectorized index algebra.

    Purpose
    -------
    The array form of :class:`MillerDirection`, and the form slip-system
    enumeration and inverse-pole-figure family expansion consume.

    Attributes
    ----------
    indices : np.ndarray
        ``(n, 3)`` integer indices; no row may be the zero triplet.
    phase : Phase
        Shared by every direction in the set.
    """

    indices: np.ndarray
    phase: Phase

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", as_int_array(self.indices, shape=(None, 3)))
        _validate_nonzero_rows(cast(IntArray, self.indices), name="MillerDirectionSet.indices")

    @classmethod
    def from_uvw(cls, indices: Any, *, phase: Phase) -> MillerDirectionSet:
        """Direction set from an ``(n, 3)`` array of ``[uvw]`` indices on one phase.
        """

        return cls(indices=_as_index_rows(indices, name="uvw", columns=3), phase=phase)

    @classmethod
    def from_UVTW(cls, indices: Any, *, phase: Phase) -> MillerDirectionSet:  # noqa: N802
        """Direction set from an ``(n, 4)`` array of ``[UVTW]`` indices on one phase.
        """

        return cls(indices=direction_uvtw_to_uvw_array(indices), phase=phase)

    def to_UVTW(self) -> IntArray:  # noqa: N802
        """``(n, 4)`` hexagonal four-index form of the whole set.
        """

        return direction_uvw_to_uvtw_array(self.indices)

    def reduce_indices(self) -> IntArray:
        """``(n, 3)`` indices divided row-wise by their greatest common divisor.
        """

        return reduce_indices(self.indices)

    def canonical_indices(self, *, antipodal: bool = True) -> IntArray:
        """``(n, 3)`` reduced indices, optionally made sign-canonical.

        Pass ``antipodal=False`` to keep ``[uvw]`` distinct from ``[-u-v-w]``,
        which matters wherever the sense of the direction carries physics.
        """

        return canonicalize_family_indices(self.indices, antipodal=antipodal)

    def unique(self, *, antipodal: bool = True) -> tuple[MillerDirectionSet, IntArray]:
        """The distinct directions in the set, up to reduction and (by default)
        inversion.

        Returns
        -------
        tuple of (MillerDirectionSet, IntArray)
            The deduplicated set and an ``(n,)`` inverse-index array mapping
            every original row to its position in it.
        """

        unique_indices, inverse = _family_unique_rows(self.canonical_indices(antipodal=antipodal))
        return MillerDirectionSet(indices=unique_indices, phase=self.phase), inverse

    def direct_vectors_cartesian(self) -> FloatArray:
        """``(n, 3)`` lattice vectors in the Cartesian crystal frame, in angstroms.

        Not normalized; lengths are the repeat distances along each direction.
        """

        return _cartesian_from_direct_indices(cast(IntArray, self.indices), self.phase)

    def unit_vectors_cartesian(self) -> FloatArray:
        """``(n, 3)`` unit direction vectors in the Cartesian crystal frame.
        """

        return _unit_vectors(self.direct_vectors_cartesian(), name="MillerDirectionSet")

    def symmetry_equivalent_indices(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[IntArray, np.ndarray]:
        """Symmetry-equivalent indices of every direction, as a padded array.

        Purpose
        -------
        Enumerate the symmetry family <uvw> of each direction: the indices that
        the crystal's point group makes crystallographically indistinguishable.

        Parameters
        ----------
        unique : bool
            Collapse repeated members of the orbit (default ``True``). Repeats
            occur whenever an operator maps a direction onto itself.
        antipodal : bool
            Treat ``h`` and ``-h`` as the same family member. Directions have a genuine sense, so
            ``antipodal=False`` is meaningful and keeps ``[uvw]`` and
            ``[-u-v-w]`` distinct.
        tol : float
            Tolerance on recovering integer indices after applying the Cartesian
            symmetry operators. Exceeding it raises rather than silently
            rounding to wrong indices.

        Returns
        -------
        tuple of (IntArray, np.ndarray)
            ``(n, m, 3)`` indices and an ``(n, m)`` boolean validity mask. Rows
            are padded to the longest family in the batch, so the mask — not the
            array shape — says how many members each row really has.
        """

        return _symmetry_equivalent_indices(
            cast(IntArray, self.indices),
            phase=self.phase,
            reciprocal=False,
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )

    def symmetry_equivalents(
        self,
        *,
        unique: bool = True,
        antipodal: bool = True,
        tol: float = 1e-10,
    ) -> tuple[MillerDirectionSet, ...]:
        """One direction set per input direction, holding that direction's
        ``<uvw>`` family.

        The typed form of :meth:`symmetry_equivalent_indices`, with the padding
        mask applied per row.
        """

        equivalent_indices, mask = self.symmetry_equivalent_indices(
            unique=unique,
            antipodal=antipodal,
            tol=tol,
        )
        return tuple(
            MillerDirectionSet(indices=equivalent_indices[row, mask[row]], phase=self.phase)
            for row in range(equivalent_indices.shape[0])
        )

    def to_zone_axes(self) -> tuple[ZoneAxis, ...]:
        """The zone-axis view of every direction in the set.
        """

        return tuple(
            MillerDirection(indices=row, phase=self.phase).to_zone_axis()
            for row in self.indices
        )

    def angle_matrix_rad(
        self,
        other: MillerDirectionSet | None = None,
        *,
        antipodal: bool = True,
    ) -> FloatArray:
        """Pairwise angles between crystal directions, in radians.

        Parameters
        ----------
        other : MillerDirectionSet, optional
            The second set; defaults to this set, giving the self-angle matrix
            with a zero diagonal. Must be on the same phase.
        antipodal : bool
            When ``True`` (default) a direction and its reverse are treated as
            equivalent and angles lie in ``[0, pi/2]``. Set ``False`` to keep the
            sense, giving angles in ``[0, pi]`` — the right choice for slip
            directions and Burgers vectors.

        Returns
        -------
        FloatArray
            ``(len(self), len(other))`` angles in radians.
        """

        target = self if other is None else other
        _require_matching_phases(
            self.phase,
            target.phase,
            left_name="MillerDirectionSet",
            right_name="other",
        )
        return _pairwise_matrix_angles_from_unit_vectors(
            self.unit_vectors_cartesian(),
            target.unit_vectors_cartesian(),
            antipodal=antipodal,
        )


def _plane_normals(planes: MillerPlane | MillerPlaneSet) -> tuple[FloatArray, Phase]:
    if isinstance(planes, MillerPlane):
        return planes.normal_cartesian[None, :], planes.phase
    return planes.normals_cartesian(), planes.phase


def _direction_units(directions: MillerDirection | MillerDirectionSet) -> tuple[FloatArray, Phase]:
    if isinstance(directions, MillerDirection):
        return directions.unit_vector_cartesian[None, :], directions.phase
    return directions.unit_vectors_cartesian(), directions.phase


def angle_plane_plane_rad(
    left: MillerPlane | MillerPlaneSet,
    right: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    r"""Angle between crystallographic planes, computed from their normals.

    Purpose
    -------
    Return the angle (radians) between the normals of ``(hkl)`` planes, evaluated
    in the shared crystal frame through the reciprocal metric. Plane families are
    treated with antipodal equivalence, so a normal and its opposite describe the
    same plane and the returned angle lies in ``[0, pi/2]``.

    When to use
    -----------
    Use this when relating indexed poles in a pole figure, checking Kikuchi band
    intersections, or confirming that a phase's frame and symmetry are wired
    correctly. For a cubic phase the ``(100)`` vs ``(110)`` angle is exactly
    ``45`` degrees for any lattice parameter, a convenient first sanity check.

    Parameters
    ----------
    left, right : MillerPlane or MillerPlaneSet
        Planes on the same phase. Two sets give row-aligned pairwise angles.

    Returns
    -------
    float or FloatArray
        A ``float`` when both arguments are single planes; otherwise a 1-D array
        of angles in radians.

    See Also
    --------
    angle_dir_dir_rad : Angle between crystallographic directions.
    angle_dir_plane_normal_rad : Angle between a direction and a plane normal.

    Notes
    -----
    Verified live by the ``cubic-angle-100-110`` and ``hex-angle-basal-prism``
    worked examples; see the executable-examples gallery under ``docs/site/examples/``.
    """

    left_normals, left_phase = _plane_normals(left)
    right_normals, right_phase = _plane_normals(right)
    _require_matching_phases(left_phase, right_phase, left_name="left", right_name="right")
    angles = _pairwise_angles_from_unit_vectors(
        left_normals,
        right_normals,
        left_name="left",
        right_name="right",
        antipodal=True,
    )
    if isinstance(left, MillerPlane) and isinstance(right, MillerPlane):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_dir_rad(
    left: MillerDirection | MillerDirectionSet,
    right: MillerDirection | MillerDirectionSet,
    *,
    antipodal: bool = True,
) -> float | FloatArray:
    r"""Angle between crystallographic directions, computed in the crystal frame.

    Purpose
    -------
    Return the angle (radians) between ``[uvw]`` directions, evaluated with the
    direct-space metric. With ``antipodal=True`` (the default) a direction and its
    reverse are equivalent and the angle lies in ``[0, pi/2]``; set
    ``antipodal=False`` to keep signed directions and allow angles up to ``pi``.

    When to use
    -----------
    Use this in slip-system and Schmid-factor work (angle between a slip direction
    and a loading axis), when measuring the opening of a zone, or when validating
    hexagonal four-index handling. In cubic metrics the ``[110]`` vs ``[111]``
    angle is ``arccos(sqrt(2/3)) = 35.2644`` degrees.

    Parameters
    ----------
    left, right : MillerDirection or MillerDirectionSet
        Directions on the same phase. Two sets give row-aligned pairwise angles.
    antipodal : bool, optional
        If ``True`` (default) treat a direction and its reverse as equivalent.

    Returns
    -------
    float or FloatArray
        A ``float`` when both arguments are single directions; otherwise a 1-D
        array of angles in radians.

    See Also
    --------
    angle_plane_plane_rad : Angle between plane normals.

    Notes
    -----
    Verified live by the ``cubic-angle-dir-110-111`` worked example; see the
    executable-examples gallery under ``docs/site/examples/``.
    """

    left_units, left_phase = _direction_units(left)
    right_units, right_phase = _direction_units(right)
    _require_matching_phases(left_phase, right_phase, left_name="left", right_name="right")
    angles = _pairwise_angles_from_unit_vectors(
        left_units,
        right_units,
        left_name="left",
        right_name="right",
        antipodal=antipodal,
    )
    if isinstance(left, MillerDirection) and isinstance(right, MillerDirection):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_plane_normal_rad(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    r"""Angle between a direction and a plane normal.

    Purpose
    -------
    Return the angle (radians) between ``[uvw]`` direction(s) and the normal(s) of
    ``(hkl)`` plane(s), in ``[0, pi/2]`` under antipodal equivalence. This is the
    complement of the inclination of the direction to the plane itself.

    When to use
    -----------
    Use this to test whether a direction lies in a plane (angle ``= pi/2``) or is
    parallel to its normal (angle ``= 0``), for example when checking zone-axis
    conditions or resolving a Burgers vector against a slip-plane normal. For the
    inclination to the plane surface, use :func:`angle_dir_plane_inclination_rad`.

    Parameters
    ----------
    directions : MillerDirection or MillerDirectionSet
        Directions on the same phase as ``planes``.
    planes : MillerPlane or MillerPlaneSet
        Planes on the same phase as ``directions``.

    Returns
    -------
    float or FloatArray
        A ``float`` for a single direction and plane; otherwise a 1-D array.

    See Also
    --------
    angle_dir_plane_inclination_rad : Inclination of a direction to a plane.
    """

    direction_units, direction_phase = _direction_units(directions)
    plane_normals, plane_phase = _plane_normals(planes)
    _require_matching_phases(
        direction_phase,
        plane_phase,
        left_name="directions",
        right_name="planes",
    )
    angles = _pairwise_angles_from_unit_vectors(
        direction_units,
        plane_normals,
        left_name="directions",
        right_name="planes",
        antipodal=True,
    )
    if isinstance(directions, MillerDirection) and isinstance(planes, MillerPlane):
        return float(angles[0])
    return as_float_array(angles, shape=(angles.shape[0],))


def angle_dir_plane_inclination_rad(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> float | FloatArray:
    r"""Inclination of a direction to a plane surface (radians).

    Purpose
    -------
    Return ``pi/2`` minus the direction-to-normal angle, i.e. the angle between a
    direction and the plane itself. A direction lying in the plane gives ``0``; a
    direction along the plane normal gives ``pi/2``.

    When to use
    -----------
    Use this when the natural quantity is how steeply a direction rises out of a
    plane, such as the inclination of a slip or trace direction to a surface.

    Parameters
    ----------
    directions : MillerDirection or MillerDirectionSet
        Directions on the same phase as ``planes``.
    planes : MillerPlane or MillerPlaneSet
        Planes on the same phase as ``directions``.

    Returns
    -------
    float or FloatArray
        A ``float`` for a single direction and plane; otherwise a 1-D array.

    See Also
    --------
    angle_dir_plane_normal_rad : Angle between a direction and a plane normal.
    """

    normal_angles = angle_dir_plane_normal_rad(directions, planes)
    if isinstance(normal_angles, float):
        return float((np.pi / 2.0) - normal_angles)
    return as_float_array((np.pi / 2.0) - normal_angles, shape=(normal_angles.shape[0],))


def project_directions_onto_planes(
    directions: MillerDirection | MillerDirectionSet,
    planes: MillerPlane | MillerPlaneSet,
) -> tuple[FloatArray, np.ndarray]:
    """Components of directions lying within given planes.

    Purpose
    -------
    Remove the plane-normal component from each direction, leaving the
    in-plane part. This is the geometric step behind slip-trace analysis,
    behind resolving a direction into a habit or boundary plane, and behind
    checking how far an intended in-plane direction actually departs from
    the plane.

    Parameters
    ----------
    directions : MillerDirection or MillerDirectionSet
    planes : MillerPlane or MillerPlaneSet
        Must be on the same phase as ``directions``. One row broadcasts
        against many.

    Returns
    -------
    tuple of (FloatArray, np.ndarray)
        ``(n, 3)`` projected Cartesian vectors in the crystal frame — *not*
        normalized, so their length reports how much of the direction
        survived the projection — and an ``(n,)`` boolean mask flagging
        degenerate rows where the direction was parallel to the plane normal
        and nothing survived. Degenerate rows are returned as exact zeros.
    """

    direction_units, direction_phase = _direction_units(directions)
    plane_normals, plane_phase = _plane_normals(planes)
    _require_matching_phases(
        direction_phase,
        plane_phase,
        left_name="directions",
        right_name="planes",
    )
    return _projection_vectors(direction_units, plane_normals)


__all__ = [
    "MillerBravaisDirection",
    "MillerBravaisPlane",
    "MillerDirection",
    "MillerDirectionSet",
    "MillerPlane",
    "MillerPlaneSet",
    "angle_dir_dir_rad",
    "angle_dir_plane_inclination_rad",
    "angle_dir_plane_normal_rad",
    "angle_plane_plane_rad",
    "antipodal_keys",
    "canonicalize_family_indices",
    "canonicalize_sign",
    "direction_uvtw_to_uvw_array",
    "direction_uvw_to_uvtw_array",
    "plane_hkil_to_hkl_array",
    "plane_hkl_to_hkil_array",
    "project_directions_onto_planes",
    "reduce_indices",
    "zone_law_value_hkil_uvtw",
]
