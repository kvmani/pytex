"""Tests for the core reference-frame model: geometry, invariants, and transforms."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    FrameDomain,
    FrameTransform,
    Handedness,
    ReferenceFrame,
    Rotation,
    VectorSet,
)
from pytex.core.frame_catalog import (
    SAMPLE_RD_TD_ND_FRAME,
    SPECIMEN_FRAME,
    crystal_frame,
    map_frame,
    sample_frame,
    specimen_frame,
)
from pytex.core.frames import IDENTITY_AXIS_VECTORS, as_axis_vectors


def make_frame(name: str, domain: FrameDomain) -> ReferenceFrame:
    return ReferenceFrame(
        name=name,
        domain=domain,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )


# ---------------------------------------------------------------------------
# ReferenceFrame: geometry and invariants
# ---------------------------------------------------------------------------


def test_default_frame_is_the_identity_triad() -> None:
    frame = make_frame("specimen", FrameDomain.SPECIMEN)
    assert frame.axis_vectors == IDENTITY_AXIS_VECTORS
    assert_allclose(frame.basis_matrix, np.eye(3), atol=1e-12)
    assert frame.is_orthonormal
    assert frame.is_right_handed
    assert frame.determinant == pytest.approx(1.0)


def test_basis_matrix_columns_are_the_axis_vectors() -> None:
    frame = specimen_frame(
        axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert_allclose(frame.basis_matrix[:, 0], [0.0, 1.0, 0.0], atol=1e-12)
    assert_allclose(frame.basis_matrix[:, 1], [-1.0, 0.0, 0.0], atol=1e-12)


def test_basis_matrix_maps_frame_components_to_cartesian() -> None:
    # A frame whose x axis points along canonical Y: a vector with components
    # (1, 0, 0) in that frame is the canonical Y direction.
    frame = specimen_frame(
        axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    )
    cartesian = frame.basis_matrix @ np.array([1.0, 0.0, 0.0])
    assert_allclose(cartesian, [0.0, 1.0, 0.0], atol=1e-12)


def test_frame_rejects_linearly_dependent_axes() -> None:
    with pytest.raises(ValueError, match="linearly dependent"):
        specimen_frame(axis_vectors=[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


def test_frame_rejects_handedness_contradicting_the_axis_determinant() -> None:
    with pytest.raises(ValueError, match="the sign must match"):
        ReferenceFrame(
            name="mirrored",
            domain=FrameDomain.SPECIMEN,
            axes=("x", "y", "z"),
            handedness=Handedness.RIGHT,
            axis_vectors=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
        )


def test_left_handed_frame_is_accepted_when_declared() -> None:
    frame = ReferenceFrame(
        name="mirrored",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.LEFT,
        axis_vectors=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
    )
    assert not frame.is_right_handed
    assert frame.determinant < 0.0


def test_frame_rejects_non_finite_axis_vectors() -> None:
    with pytest.raises(ValueError, match="must be finite"):
        specimen_frame(axis_vectors=[[np.nan, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


def test_frame_rejects_wrong_axis_label_count() -> None:
    with pytest.raises(ValueError, match="exactly three axis labels"):
        ReferenceFrame(name="bad", domain=FrameDomain.SPECIMEN, axes=("x", "y"))


def test_frame_rejects_partial_axis_descriptions() -> None:
    with pytest.raises(ValueError, match="empty or contain exactly three"):
        ReferenceFrame(
            name="partial",
            domain=FrameDomain.SPECIMEN,
            axes=("x", "y", "z"),
            axis_descriptions=("only one",),
        )


def test_oblique_frame_is_allowed_but_reported_non_orthonormal() -> None:
    frame = crystal_frame(
        axis_vectors=[[1.0, 0.0, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.0]]
    )
    assert not frame.is_orthonormal
    assert frame.is_right_handed


def test_frames_remain_comparable_and_new_fields_do_not_break_equality() -> None:
    assert specimen_frame() == specimen_frame()
    assert specimen_frame() != map_frame()
    assert specimen_frame() != specimen_frame(axis_vectors=np.diag([1.0, -1.0, -1.0]))


# ---------------------------------------------------------------------------
# ReferenceFrame: axis access
# ---------------------------------------------------------------------------


def test_axis_index_resolves_labels_case_insensitively_and_integers() -> None:
    frame = SAMPLE_RD_TD_ND_FRAME
    assert frame.axis_index("RD") == 0
    assert frame.axis_index("td") == 1
    assert frame.axis_index(2) == 2


def test_axis_index_rejects_unknown_labels_and_out_of_range_indices() -> None:
    frame = SAMPLE_RD_TD_ND_FRAME
    with pytest.raises(KeyError, match="has no axis"):
        frame.axis_index("ED")
    with pytest.raises(IndexError, match="out of range"):
        frame.axis_index(3)


def test_axis_vector_normalizes_by_default() -> None:
    frame = crystal_frame(axis_vectors=[[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]])
    assert_allclose(frame.axis_vector("a"), [1.0, 0.0, 0.0], atol=1e-12)
    assert_allclose(frame.axis_vector("a", normalize=False), [3.0, 0.0, 0.0], atol=1e-12)


def test_unit_axis_matrix_scales_every_column_to_unit_length() -> None:
    frame = crystal_frame(axis_vectors=[[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]])
    assert_allclose(np.linalg.norm(frame.unit_axis_matrix(), axis=0), np.ones(3), atol=1e-12)


def test_axis_description_falls_back_to_the_label() -> None:
    assert SAMPLE_RD_TD_ND_FRAME.axis_description("RD") == "rolling direction"
    assert specimen_frame().axis_description("x") == "x"


def test_returned_arrays_are_read_only() -> None:
    frame = specimen_frame()
    assert not frame.basis_matrix.flags.writeable
    assert not frame.axis_vector("x").flags.writeable


# ---------------------------------------------------------------------------
# ReferenceFrame: derivation
# ---------------------------------------------------------------------------


def test_with_axis_vectors_preserves_every_other_field() -> None:
    original = sample_frame(name="sheet", description="rolled sheet")
    rotated = original.with_axis_vectors(np.diag([1.0, -1.0, -1.0]))
    assert rotated.name == original.name
    assert rotated.axes == original.axes
    assert rotated.description == original.description
    assert rotated.axis_descriptions == original.axis_descriptions
    assert rotated.axis_vectors != original.axis_vectors


def test_renamed_changes_only_the_name_and_optional_description() -> None:
    renamed = specimen_frame().renamed("stage", description="stage-fixed")
    assert renamed.name == "stage"
    assert renamed.description == "stage-fixed"
    assert renamed.domain is FrameDomain.SPECIMEN


def test_rotated_turns_the_axes_by_the_given_rotation() -> None:
    frame = sample_frame()
    turned = frame.rotated(Rotation.from_axis_angle((0.0, 0.0, 1.0), np.pi / 2.0))
    assert_allclose(turned.axis_vector("RD"), [0.0, 1.0, 0.0], atol=1e-12)
    assert_allclose(turned.axis_vector("ND"), [0.0, 0.0, 1.0], atol=1e-12)
    assert turned.name == "sample_rd_td_nd_rotated"
    assert turned.is_right_handed


def test_rotated_accepts_an_explicit_name() -> None:
    turned = sample_frame().rotated(Rotation.identity(), name="mounted")
    assert turned.name == "mounted"


def test_as_axis_vectors_rejects_wrong_shapes() -> None:
    with pytest.raises(ValueError, match="exactly three 3-component"):
        as_axis_vectors([[1.0, 0.0], [0.0, 1.0]])


# ---------------------------------------------------------------------------
# ReferenceFrame: describe()
# ---------------------------------------------------------------------------


def test_describe_states_domain_handedness_axes_and_convention() -> None:
    text = SAMPLE_RD_TD_ND_FRAME.describe()
    assert "sample_rd_td_nd" in text
    assert "specimen domain" in text
    assert "right-handed" in text
    assert "orthonormal" in text
    assert "rolling direction" in text
    assert "canonical Cartesian reference" in text
    assert "bunge_zxz" in text


def test_describe_flags_a_non_orthonormal_frame() -> None:
    frame = crystal_frame(
        axis_vectors=[[1.0, 0.0, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.0]]
    )
    assert "non-orthonormal" in frame.describe()


# ---------------------------------------------------------------------------
# FrameTransform: legacy behaviour that must not regress
# ---------------------------------------------------------------------------


def test_identity_transform_round_trip() -> None:
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    transform = FrameTransform.identity(specimen)
    vectors = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, 6.0]])
    assert_allclose(transform.apply_to_vectors(vectors), vectors)
    assert_allclose(transform.inverse().apply_to_vectors(vectors), vectors)


def test_transform_composition_matches_stepwise_application() -> None:
    crystal = make_frame("crystal", FrameDomain.CRYSTAL)
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    lab = make_frame("lab", FrameDomain.LABORATORY)

    first = FrameTransform(
        source=crystal,
        target=specimen,
        rotation_matrix=np.eye(3),
        translation_vector=np.array([1.0, 0.0, 0.0]),
    )
    second = FrameTransform(
        source=specimen,
        target=lab,
        rotation_matrix=np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        translation_vector=np.array([0.0, 2.0, 0.0]),
    )
    chained = second.compose(first)
    vector = np.array([[1.0, 2.0, 3.0]])
    stepwise = second.apply_to_vectors(first.apply_to_vectors(vector))
    assert_allclose(chained.apply_to_vectors(vector), stepwise)


def test_transform_updates_vector_set_frame() -> None:
    crystal = make_frame("crystal", FrameDomain.CRYSTAL)
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    transform = FrameTransform(
        source=crystal,
        target=specimen,
        rotation_matrix=np.eye(3),
    )
    vectors = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=crystal)
    transformed = transform.apply_to_vectors(vectors)
    assert isinstance(transformed, VectorSet)
    assert transformed.reference_frame == specimen


def test_transform_rejects_a_non_rotation_matrix() -> None:
    frame = make_frame("specimen", FrameDomain.SPECIMEN)
    with pytest.raises(ValueError, match="orthonormal"):
        FrameTransform(
            source=frame, target=frame, rotation_matrix=np.diag([1.0, 1.0, 2.0])
        )
    with pytest.raises(ValueError, match="determinant"):
        FrameTransform(
            source=frame, target=frame, rotation_matrix=np.diag([1.0, 1.0, -1.0])
        )


def test_compose_rejects_a_broken_chain() -> None:
    crystal = make_frame("crystal", FrameDomain.CRYSTAL)
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    lab = make_frame("lab", FrameDomain.LABORATORY)
    first = FrameTransform.identity(crystal)
    second = FrameTransform(source=specimen, target=lab, rotation_matrix=np.eye(3))
    with pytest.raises(ValueError, match="do not chain"):
        second.compose(first)


def test_apply_to_vectors_rejects_a_vector_set_from_another_frame() -> None:
    crystal = make_frame("crystal", FrameDomain.CRYSTAL)
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    transform = FrameTransform(source=crystal, target=specimen, rotation_matrix=np.eye(3))
    stray = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=specimen)
    with pytest.raises(ValueError, match="must match FrameTransform.source"):
        transform.apply_to_vectors(stray)


# ---------------------------------------------------------------------------
# FrameTransform: constructors
# ---------------------------------------------------------------------------


def test_from_rotation_uses_the_rotation_matrix_directly() -> None:
    rotation = Rotation.from_axis_angle((0.0, 0.0, 1.0), np.pi / 2.0)
    transform = FrameTransform.from_rotation(
        rotation, source=specimen_frame(), target=sample_frame()
    )
    assert_allclose(transform.rotation_matrix, rotation.as_matrix(), atol=1e-12)
    assert transform.rotation_angle_deg == pytest.approx(90.0)


def test_from_bunge_euler_matches_the_rotation_helper() -> None:
    transform = FrameTransform.from_bunge_euler(
        30.0, 20.0, 10.0, source=crystal_frame(), target=specimen_frame()
    )
    assert_allclose(
        transform.rotation_matrix,
        Rotation.from_bunge_euler(30.0, 20.0, 10.0).as_matrix(),
        atol=1e-12,
    )


def test_from_axis_angle_accepts_degrees_and_radians() -> None:
    degrees = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0), 45.0, source=specimen_frame(), target=sample_frame()
    )
    radians = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0),
        np.pi / 4.0,
        source=specimen_frame(),
        target=sample_frame(),
        degrees=False,
    )
    assert_allclose(degrees.rotation_matrix, radians.rotation_matrix, atol=1e-12)
    assert degrees.rotation_angle_deg == pytest.approx(45.0)


def test_axis_correspondence_maps_each_source_axis_onto_its_declared_partner() -> None:
    transform = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "TD", "y": "-RD", "z": "ND"}
    )
    # specimen x is the sample TD axis, so its sample components are (0, 1, 0).
    assert_allclose(
        np.asarray(transform.apply_to_directions(np.array([1.0, 0.0, 0.0]))),
        [0.0, 1.0, 0.0],
        atol=1e-12,
    )
    # specimen y is the reversed RD axis, so its sample components are (-1, 0, 0).
    assert_allclose(
        np.asarray(transform.apply_to_directions(np.array([0.0, 1.0, 0.0]))),
        [-1.0, 0.0, 0.0],
        atol=1e-12,
    )


def test_axis_correspondence_is_about_components_not_where_frames_point() -> None:
    # The declaration fixes component semantics, so it must not depend on either
    # frame's canonical-Cartesian axis geometry.
    plain = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "RD", "y": "TD", "z": "ND"}
    )
    tilted_target = sample_frame(
        axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    )
    tilted = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME, tilted_target, {"x": "RD", "y": "TD", "z": "ND"}
    )
    assert_allclose(plain.rotation_matrix, tilted.rotation_matrix, atol=1e-12)
    assert plain.is_identity


def test_axis_correspondence_accepts_an_explicit_plus_sign() -> None:
    transform = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "+RD", "y": "+TD", "z": "+ND"}
    )
    assert transform.is_identity


def test_axis_correspondence_rejects_incomplete_or_repeated_declarations() -> None:
    with pytest.raises(ValueError, match="exactly three source-axis entries"):
        FrameTransform.from_axis_correspondence(
            SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "RD", "y": "TD"}
        )
    with pytest.raises(ValueError, match="more than once"):
        FrameTransform.from_axis_correspondence(
            SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "RD", "y": "RD", "z": "ND"}
        )


def test_axis_correspondence_rejects_an_improper_permutation() -> None:
    # Swapping two axes without a sign flip is a mirror, not a rotation.
    with pytest.raises(ValueError, match="not a proper rotation|determinant is not"):
        FrameTransform.from_axis_correspondence(
            SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "TD", "y": "RD", "z": "ND"}
        )


def test_between_frames_derives_the_relationship_from_axis_geometry() -> None:
    tilted = sample_frame(
        axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    )
    transform = FrameTransform.between_frames(SPECIMEN_FRAME, tilted)
    # Specimen x is canonical X, which is the tilted frame's -TD direction.
    assert_allclose(
        np.asarray(transform.apply_to_directions(np.array([1.0, 0.0, 0.0]))),
        [0.0, -1.0, 0.0],
        atol=1e-12,
    )


def test_between_frames_refuses_non_orthonormal_frames() -> None:
    oblique = crystal_frame(
        axis_vectors=[[1.0, 0.0, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.0]]
    )
    with pytest.raises(ValueError, match="non-orthonormal"):
        FrameTransform.between_frames(oblique, specimen_frame())


# ---------------------------------------------------------------------------
# FrameTransform: properties and application
# ---------------------------------------------------------------------------


def test_identity_flag_and_rotation_properties() -> None:
    identity = FrameTransform.identity(specimen_frame())
    assert identity.is_identity
    assert identity.rotation_angle_deg == pytest.approx(0.0)

    turned = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0), 60.0, source=specimen_frame(), target=sample_frame()
    )
    assert not turned.is_identity
    assert turned.rotation_angle_deg == pytest.approx(60.0)
    assert_allclose(np.abs(turned.rotation_axis), [0.0, 0.0, 1.0], atol=1e-12)


def test_a_translation_alone_is_not_the_identity() -> None:
    transform = FrameTransform(
        source=specimen_frame(),
        target=map_frame(),
        rotation_matrix=np.eye(3),
        translation_vector=np.array([0.5, 0.0, 0.0]),
    )
    assert not transform.is_identity


def test_directions_ignore_the_origin_offset_but_positions_do_not() -> None:
    transform = FrameTransform(
        source=specimen_frame(),
        target=map_frame(),
        rotation_matrix=np.eye(3),
        translation_vector=np.array([10.0, 0.0, 0.0]),
    )
    direction = np.array([[1.0, 0.0, 0.0]])
    assert_allclose(np.asarray(transform.apply_to_directions(direction)), direction, atol=1e-12)
    assert_allclose(
        np.asarray(transform.apply_to_vectors(direction)), [[11.0, 0.0, 0.0]], atol=1e-12
    )


def test_apply_to_directions_preserves_vector_set_typing() -> None:
    transform = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0), 90.0, source=specimen_frame(), target=map_frame()
    )
    vectors = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=specimen_frame())
    mapped = transform.apply_to_directions(vectors)
    assert isinstance(mapped, VectorSet)
    assert mapped.reference_frame == map_frame()
    assert_allclose(mapped.values[0], [0.0, 1.0, 0.0], atol=1e-12)


def test_apply_to_directions_rejects_a_vector_set_from_another_frame() -> None:
    transform = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0), 90.0, source=specimen_frame(), target=map_frame()
    )
    stray = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=map_frame())
    with pytest.raises(ValueError, match="must match FrameTransform.source"):
        transform.apply_to_directions(stray)


def test_apply_helpers_reject_arrays_that_do_not_end_in_three() -> None:
    transform = FrameTransform.identity(specimen_frame())
    with pytest.raises(ValueError, match="end with dimension 3"):
        transform.apply_to_vectors(np.zeros((2, 4)))
    with pytest.raises(ValueError, match="end with dimension 3"):
        transform.apply_to_directions(np.zeros((2, 4)))


def test_source_axes_in_target_returns_the_rotation_columns() -> None:
    transform = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME, SAMPLE_RD_TD_ND_FRAME, {"x": "TD", "y": "-RD", "z": "ND"}
    )
    axes = transform.source_axes_in_target()
    # Column 0 is specimen x expressed in sample components: the TD axis.
    assert_allclose(axes[:, 0], [0.0, 1.0, 0.0], atol=1e-12)
    assert_allclose(axes, transform.rotation_matrix, atol=1e-12)


def test_as_rotation_round_trips_through_the_rotation_type() -> None:
    transform = FrameTransform.from_bunge_euler(
        15.0, 25.0, 35.0, source=crystal_frame(), target=specimen_frame()
    )
    assert_allclose(transform.as_rotation().as_matrix(), transform.rotation_matrix, atol=1e-12)


def test_inverse_undoes_the_rotation_and_the_translation() -> None:
    transform = FrameTransform(
        source=specimen_frame(),
        target=map_frame(),
        rotation_matrix=Rotation.from_axis_angle((0.0, 0.0, 1.0), 0.7).as_matrix(),
        translation_vector=np.array([1.0, -2.0, 3.0]),
    )
    points = np.array([[1.0, 2.0, 3.0], [-1.0, 0.5, 2.0]])
    round_tripped = transform.inverse().apply_to_vectors(transform.apply_to_vectors(points))
    assert_allclose(np.asarray(round_tripped), points, atol=1e-12)
    assert transform.inverse().source == map_frame()


# ---------------------------------------------------------------------------
# FrameTransform: describe()
# ---------------------------------------------------------------------------


def test_transform_describe_states_endpoints_angle_axis_and_direction() -> None:
    transform = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0), 42.0, source=specimen_frame(), target=sample_frame()
    )
    text = transform.describe()
    assert "specimen" in text
    assert "sample_rd_td_nd" in text
    assert "42.0000 deg" in text
    assert "no origin offset" in text
    assert "v_sample_rd_td_nd = R v_specimen + t" in text


def test_transform_describe_reports_identity_and_translation() -> None:
    identity_text = FrameTransform.identity(specimen_frame()).describe()
    assert "identity map" in identity_text

    offset = FrameTransform(
        source=specimen_frame(),
        target=map_frame(),
        rotation_matrix=np.eye(3),
        translation_vector=np.array([1.0, 2.0, 3.0]),
    )
    assert "origin offset of" in offset.describe()
