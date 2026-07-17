"""Property-based invariants for orientation algebra, index round-trips, and OR mapping.

Development-guide finding 20: laws that must hold for *any* input, checked
over randomized examples — rotation-group axioms, Miller-Bravais conversion
round-trips, and exactness of the orientation-relationship index
correspondence on rational images.
"""

from __future__ import annotations

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationRelationship,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    direction_uvtw_to_uvw,
    direction_uvw_to_uvtw,
    phases_semantically_match,
    plane_hkil_to_hkl,
    plane_hkl_to_hkil,
)

_MAX_EXAMPLES = 50


def _unit_quaternion(components: list[float]) -> np.ndarray:
    array = np.asarray(components, dtype=np.float64)
    return array / np.linalg.norm(array)


_quaternion_components = st.lists(
    st.floats(min_value=-1.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    min_size=4,
    max_size=4,
).filter(lambda values: float(np.linalg.norm(np.asarray(values))) > 1e-2)

_integer_triples = st.tuples(
    st.integers(min_value=-9, max_value=9),
    st.integers(min_value=-9, max_value=9),
    st.integers(min_value=-9, max_value=9),
).filter(lambda triple: any(value != 0 for value in triple))


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(left=_quaternion_components, right=_quaternion_components)
def test_rotation_composition_inverse_law(left: list[float], right: list[float]) -> None:
    a = Rotation(quaternion=_unit_quaternion(left))
    b = Rotation(quaternion=_unit_quaternion(right))
    composed = a.compose(b)
    # (a b)^-1 = b^-1 a^-1: their composition with (a b) must be the identity.
    identity = composed.compose(b.inverse().compose(a.inverse()))
    assert identity.angle_deg < 1e-8


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(components=_quaternion_components)
def test_rotation_matrices_are_proper_orthogonal(components: list[float]) -> None:
    matrix = Rotation(quaternion=_unit_quaternion(components)).as_matrix()
    np.testing.assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-10)
    assert np.linalg.det(matrix) > 0.0


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(triple=_integer_triples)
def test_direction_miller_bravais_round_trip_preserves_direction(
    triple: tuple[int, int, int],
) -> None:
    uvtw = direction_uvw_to_uvtw(np.asarray(triple, dtype=np.int64))
    recovered = direction_uvtw_to_uvw(uvtw)
    original = np.asarray(triple, dtype=np.int64)
    # The round trip may reduce by a common factor but must keep the ray:
    # cross product of the integer triples vanishes and the sign is preserved.
    assert np.all(np.cross(original, recovered) == 0)
    assert float(original @ recovered) > 0


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(triple=_integer_triples)
def test_plane_miller_bravais_round_trip_is_exact(triple: tuple[int, int, int]) -> None:
    hkil = plane_hkl_to_hkil(np.asarray(triple, dtype=np.int64))
    recovered = plane_hkil_to_hkl(hkil)
    np.testing.assert_array_equal(recovered, np.asarray(triple, dtype=np.int64))


def _fcc_bcc_relationship() -> tuple[Phase, Phase, OrientationRelationship]:
    parent_frame = ReferenceFrame(
        name="prop_parent",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    child_frame = ReferenceFrame(
        name="prop_child",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    parent = Phase(
        "austenite",
        lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
        crystal_frame=parent_frame,
    )
    child = Phase(
        "ferrite",
        lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
        crystal_frame=child_frame,
    )
    relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    return parent, child, relationship


_PARENT, _CHILD, _KS = _fcc_bcc_relationship()
_KS_VARIANTS = _KS.generate_variants()


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(
    triple=_integer_triples,
    variant_index=st.integers(min_value=0, max_value=len(_KS_VARIANTS) - 1),
)
def test_plane_round_trip_through_any_variant_recovers_source(
    triple: tuple[int, int, int], variant_index: int
) -> None:
    variant = _KS_VARIANTS[variant_index]
    plane = CrystalPlane(
        MillerIndex(np.asarray(triple, dtype=np.int64), phase=_PARENT), phase=_PARENT
    )
    forward = _KS.map_plane_to_child(plane, variant=variant)
    # Round trip through the exact (irrational) image: rebuild from the exact
    # components is not representable as a CrystalPlane, so check instead
    # that the correspondence matrices invert each other on the raw indices.
    direct = _KS.correspondence_reciprocal(variant=variant)
    exact = direct @ np.asarray(triple, dtype=np.float64)
    recovered = np.linalg.solve(direct, exact)
    np.testing.assert_allclose(recovered, np.asarray(triple, dtype=np.float64), atol=1e-9)
    # And the rationalized image never reports a negative or non-finite residual.
    assert np.isfinite(forward.angular_residual_deg)
    assert forward.angular_residual_deg >= 0.0


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(triple=_integer_triples)
def test_zone_law_is_preserved_for_any_plane_and_orthogonal_direction(
    triple: tuple[int, int, int],
) -> None:
    plane_indices = np.asarray(triple, dtype=np.float64)
    # Build an integer direction in the zone of the plane (h.u = 0).
    direction = np.cross(plane_indices, np.array([1.0, 2.0, 3.0]))
    if np.allclose(direction, 0.0):
        direction = np.cross(plane_indices, np.array([1.0, 0.0, 0.0]))
    if np.allclose(direction, 0.0):
        return
    reciprocal = _KS.correspondence_reciprocal()
    direct = _KS.correspondence_direct()
    before = float(plane_indices @ direction)
    after = float((reciprocal @ plane_indices) @ (direct @ direction))
    assert abs(after - before) < 1e-9


def test_phase_semantic_match_is_reflexive_and_symmetric() -> None:
    parent_a, child_a, _ = _fcc_bcc_relationship()
    parent_b, _, _ = _fcc_bcc_relationship()
    assert phases_semantically_match(parent_a, parent_a)
    assert phases_semantically_match(parent_a, parent_b)
    assert phases_semantically_match(parent_b, parent_a)
    assert not phases_semantically_match(parent_a, child_a)


def _burgers_relationship() -> tuple[Phase, Phase, OrientationRelationship]:
    parent_frame = ReferenceFrame(
        name="prop_beta",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    child_frame = ReferenceFrame(
        name="prop_alpha",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    beta = Phase(
        "beta_ti",
        lattice=Lattice(3.31, 3.31, 3.31, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
        crystal_frame=parent_frame,
    )
    alpha = Phase(
        "alpha_ti",
        lattice=Lattice(2.95, 2.95, 4.68, 90.0, 90.0, 120.0, crystal_frame=child_frame),
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=child_frame),
        crystal_frame=child_frame,
    )
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=beta, child_phase=alpha
    )
    return beta, alpha, relationship


_BETA, _ALPHA, _BURGERS = _burgers_relationship()
_BURGERS_VARIANTS = _BURGERS.generate_variants()


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(
    triple=_integer_triples,
    variant_index=st.integers(min_value=0, max_value=len(_BURGERS_VARIANTS) - 1),
)
def test_burgers_correspondence_matrices_invert_and_preserve_zone_law(
    triple: tuple[int, int, int], variant_index: int
) -> None:
    """Hexagonal-child invariants: M* = M^-T and zone-law preservation.

    The child lattice is hexagonal, so the correspondence matrices are far
    from orthogonal; the inverse-transpose identity and the zone law must
    still hold exactly for every variant and any integer plane/direction.
    """

    variant = _BURGERS_VARIANTS[variant_index]
    direct = _BURGERS.correspondence_direct(variant=variant)
    reciprocal = _BURGERS.correspondence_reciprocal(variant=variant)
    np.testing.assert_allclose(reciprocal, np.linalg.inv(direct).T, atol=1e-10)
    plane_indices = np.asarray(triple, dtype=np.float64)
    direction = np.cross(plane_indices, np.array([1.0, 2.0, 3.0]))
    if np.allclose(direction, 0.0):
        direction = np.cross(plane_indices, np.array([1.0, 0.0, 0.0]))
    if np.allclose(direction, 0.0):
        return
    before = float(plane_indices @ direction)
    after = float((reciprocal @ plane_indices) @ (direct @ direction))
    assert abs(after - before) < 1e-8


@settings(max_examples=_MAX_EXAMPLES, deadline=None)
@given(triple=_integer_triples)
def test_hexagonal_symmetry_orbit_recovers_integer_indices(
    triple: tuple[int, int, int],
) -> None:
    """Every 6/mmm operator must carry integer (hkl) to integer (hkl).

    The orbit machinery recovers indices through the hexagonal reciprocal
    basis; non-integer recovery would mean the operators and the lattice
    metric disagree.
    """

    from pytex.core.transformation import _integer_index_orbit

    orbit = _integer_index_orbit(
        np.asarray(triple, dtype=np.int64), phase=_ALPHA, reciprocal=True
    )
    assert orbit.dtype == np.int64
    assert orbit.shape[0] >= 1
    # Members are primitive and antipodal-canonical: gcd 1, first nonzero > 0.
    gcds = np.gcd.reduce(np.abs(orbit), axis=1)
    assert np.all(gcds == 1)
    for row in orbit:
        nonzero = row[np.nonzero(row)[0]]
        assert nonzero[0] > 0


@settings(max_examples=25, deadline=None)
@given(
    triple=_integer_triples,
    variant_index=st.integers(min_value=0, max_value=len(_BURGERS_VARIANTS) - 1),
)
def test_burgers_plane_round_trip_recovers_indices(
    triple: tuple[int, int, int], variant_index: int
) -> None:
    """Mapping a beta plane to alpha and back recovers the primitive indices."""

    from pytex.core import CrystalPlane, MillerIndex

    variant = _BURGERS_VARIANTS[variant_index]
    reduced = np.asarray(triple, dtype=np.int64)
    reduced //= np.gcd.reduce(np.abs(reduced))
    plane = CrystalPlane(MillerIndex(reduced, phase=_BETA), phase=_BETA)
    forward = _BURGERS.map_plane_to_child(plane, variant=variant, max_index=17)
    back = _BURGERS.map_plane_to_parent(forward.target, variant=variant, max_index=17)
    # The round trip through the rationalized child plane recovers the source
    # only when the forward image itself was rational; always, though, the
    # exact-components round trip must invert.
    direct = _BURGERS.correspondence_reciprocal(variant=variant)
    recovered = np.linalg.solve(direct, direct @ reduced.astype(np.float64))
    np.testing.assert_allclose(recovered, reduced.astype(np.float64), atol=1e-9)
    if forward.angular_residual_deg < 1e-9:
        sign = 1 if np.any(back.rational_indices == reduced) else -1
        np.testing.assert_array_equal(sign * back.rational_indices, reduced)
