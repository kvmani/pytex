"""Tests for the standard reference-frame catalog and the frame graph.

The identity-preservation tests here are load-bearing: frame equality gates
`VectorSet`, `FrameTransform`, `Orientation`, and `SymmetrySpec` consistency
checks, so a catalog frame must compare equal to the hand-built frame the
repository used before the catalog existed. If one of those tests fails, a
catalog default has drifted and cross-module frame identity is broken.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import FrameDomain, Handedness, ReferenceFrame
from pytex.core.frame_catalog import (
    CARTESIAN_FRAME,
    CRYSTAL_FRAME,
    DETECTOR_FRAME,
    LABORATORY_FRAME,
    MAP_FRAME,
    SAMPLE_RD_TD_ND_FRAME,
    SPECIMEN_FRAME,
    STANDARD_FRAMES,
    cartesian_frame,
    crystal_frame,
    detector_frame,
    get_standard_frame,
    laboratory_frame,
    list_standard_frames,
    map_frame,
    reciprocal_frame_for,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)
from pytex.core.frames import FrameGraph, FrameTransform

# ---------------------------------------------------------------------------
# Identity preservation
# ---------------------------------------------------------------------------


def test_catalog_crystal_frame_equals_the_historical_hand_built_frame() -> None:
    historical = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
    )
    assert crystal_frame() == historical
    assert CRYSTAL_FRAME == historical


def test_catalog_specimen_and_map_frames_equal_the_historical_hand_built_frames() -> None:
    assert specimen_frame() == ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
    )
    assert map_frame() == ReferenceFrame(
        name="map",
        domain=FrameDomain.MAP,
        axes=("x", "y", "z"),
    )


def test_frame_built_with_explicit_right_handedness_still_compares_equal() -> None:
    # Handedness.RIGHT is the default, so passing it explicitly (as several call
    # sites historically did) must not change frame identity.
    assert crystal_frame() == ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def test_default_ebsd_frames_use_the_catalog() -> None:
    from pytex.adapters.scan_files import default_ebsd_frames

    crystal, specimen, scan_map = default_ebsd_frames()
    assert (crystal, specimen, scan_map) == (crystal_frame(), specimen_frame(), map_frame())


# ---------------------------------------------------------------------------
# Catalog contents
# ---------------------------------------------------------------------------


def test_every_standard_frame_uses_the_fixed_domain_vocabulary() -> None:
    for frame in STANDARD_FRAMES.values():
        assert isinstance(frame.domain, FrameDomain)
        assert frame.is_right_handed
        assert frame.is_orthonormal


def test_standard_frames_have_distinct_names() -> None:
    names = [frame.name for frame in STANDARD_FRAMES.values()]
    assert len(names) == len(set(names))


def test_sample_frame_carries_rolling_geometry_semantics() -> None:
    frame = SAMPLE_RD_TD_ND_FRAME
    assert frame.axes == ("RD", "TD", "ND")
    assert frame.domain is FrameDomain.SPECIMEN
    assert frame.axis_description("RD") == "rolling direction"
    assert frame.axis_description("ND") == "normal direction"
    assert_allclose(frame.axis_vector("TD"), [0.0, 1.0, 0.0], atol=1e-12)


def test_cartesian_frame_is_the_identity_reference() -> None:
    assert CARTESIAN_FRAME.axes == ("X", "Y", "Z")
    assert_allclose(CARTESIAN_FRAME.basis_matrix, np.eye(3), atol=1e-12)


def test_detector_and_laboratory_frames_are_domain_typed() -> None:
    assert DETECTOR_FRAME.domain is FrameDomain.DETECTOR
    assert DETECTOR_FRAME.axes == ("u", "v", "n")
    assert DETECTOR_FRAME.axis_description("n") == "detector-plane normal"
    assert LABORATORY_FRAME.domain is FrameDomain.LABORATORY


def test_get_standard_frame_is_case_insensitive_and_reports_unknown_slugs() -> None:
    assert get_standard_frame("SAMPLE") is SAMPLE_RD_TD_ND_FRAME
    assert get_standard_frame("map") is MAP_FRAME
    with pytest.raises(KeyError, match="not a standard PyTex frame"):
        get_standard_frame("beamline")


def test_list_standard_frames_matches_the_mapping_keys() -> None:
    assert set(list_standard_frames()) == set(STANDARD_FRAMES)


def test_builders_accept_custom_names_and_provenance() -> None:
    parent = crystal_frame("parent")
    child = crystal_frame("child")
    assert parent != child
    assert parent.domain is child.domain is FrameDomain.CRYSTAL


def test_builders_reject_wrong_axis_label_counts() -> None:
    with pytest.raises(ValueError, match="exactly three axis labels"):
        specimen_frame(axes=("x", "y"))


def test_frame_builders_accept_explicit_axis_geometry() -> None:
    tilted = sample_frame(axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert_allclose(tilted.axis_vector("RD"), [0.0, 1.0, 0.0], atol=1e-12)
    assert tilted.is_right_handed


# ---------------------------------------------------------------------------
# Reciprocal frames
# ---------------------------------------------------------------------------


def test_reciprocal_frame_stars_the_axis_labels_and_switches_domain() -> None:
    reciprocal = reciprocal_frame_for(CRYSTAL_FRAME)
    assert reciprocal.domain is FrameDomain.RECIPROCAL
    assert reciprocal.axes == ("a*", "b*", "c*")
    assert reciprocal.name == "crystal_reciprocal"
    assert "dual to a" in reciprocal.axis_description("a*")


def test_lattice_reciprocal_basis_uses_the_catalog_frame() -> None:
    from pytex.core.lattice import Lattice

    lattice = Lattice(
        a=3.6,
        b=3.6,
        c=3.6,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=crystal_frame(),
    )
    basis = lattice.reciprocal_basis()
    assert basis.frame == reciprocal_frame_for(lattice.crystal_frame)
    assert basis.frame.domain is FrameDomain.RECIPROCAL


# ---------------------------------------------------------------------------
# Frame graph
# ---------------------------------------------------------------------------


def test_rolling_frame_graph_resolves_the_full_chain() -> None:
    graph = rolling_frame_graph(rd_offset_deg=30.0)
    assert graph.path("cartesian", "sample_rd_td_nd") == (
        "cartesian",
        "specimen",
        "sample_rd_td_nd",
    )
    transform = graph.transform_between("cartesian", "sample_rd_td_nd")
    assert transform.rotation_angle_deg == pytest.approx(30.0, abs=1e-9)


def test_rolling_frame_graph_with_zero_offset_is_the_identity_relationship() -> None:
    graph = rolling_frame_graph()
    transform = graph.transform_between("cartesian", "sample_rd_td_nd")
    assert transform.is_identity


def test_rolling_offset_maps_the_stage_axis_onto_the_expected_sample_components() -> None:
    # With RD rotated +90 degrees about ND from specimen x, a vector along
    # specimen x is the sample -TD direction, i.e. components (0, -1, 0).
    graph = rolling_frame_graph(rd_offset_deg=90.0)
    mapped = graph.convert(
        np.array([1.0, 0.0, 0.0]),
        source="specimen",
        target="sample_rd_td_nd",
        directions=True,
    )
    assert_allclose(np.asarray(mapped), [0.0, -1.0, 0.0], atol=1e-12)


def test_graph_resolves_transform_in_both_directions() -> None:
    graph = rolling_frame_graph(rd_offset_deg=25.0)
    forward = graph.transform_between("specimen", "sample_rd_td_nd")
    backward = graph.transform_between("sample_rd_td_nd", "specimen")
    assert_allclose(
        forward.rotation_matrix @ backward.rotation_matrix, np.eye(3), atol=1e-12
    )


def test_graph_identity_for_same_frame() -> None:
    graph = rolling_frame_graph()
    transform = graph.transform_between("specimen", "specimen")
    assert transform.is_identity
    assert transform.source == transform.target == SPECIMEN_FRAME


def test_graph_reports_unknown_and_disconnected_frames() -> None:
    graph = rolling_frame_graph()
    graph.add_frame(detector_frame())
    with pytest.raises(KeyError, match="is not registered"):
        graph.transform_between("specimen", "nowhere")
    with pytest.raises(KeyError, match="No declared transform path"):
        graph.transform_between("specimen", "detector")


def test_graph_rejects_conflicting_definitions_of_one_frame_name() -> None:
    graph = FrameGraph(name="conflict")
    graph.add_frame(specimen_frame())
    with pytest.raises(ValueError, match="already registered with different definition"):
        graph.add_frame(specimen_frame(axes=("a", "b", "c")))


def test_graph_registering_an_identical_frame_twice_is_a_no_op() -> None:
    graph = FrameGraph(frames=[specimen_frame(), specimen_frame()])
    assert len(graph) == 1


def test_graph_chooses_the_shortest_path_when_a_shortcut_exists() -> None:
    graph = rolling_frame_graph(rd_offset_deg=40.0)
    direct = FrameTransform.from_axis_angle(
        (0.0, 0.0, 1.0),
        40.0,
        source=cartesian_frame(),
        target=sample_frame(
            axis_vectors=rolling_frame_graph(rd_offset_deg=40.0)
            .frame("sample_rd_td_nd")
            .axis_vectors
        ),
    )
    graph.add_transform(direct)
    assert graph.path("cartesian", "sample_rd_td_nd") == ("cartesian", "sample_rd_td_nd")


def test_graph_membership_and_lookup() -> None:
    graph = rolling_frame_graph()
    assert "specimen" in graph
    assert SPECIMEN_FRAME in graph
    assert graph.frame("specimen") == SPECIMEN_FRAME
    assert not graph.has_frame("detector")
    with pytest.raises(KeyError, match="Registered frames"):
        graph.frame("detector")


def test_graph_frames_and_transforms_are_exposed_as_tuples() -> None:
    graph = rolling_frame_graph()
    assert isinstance(graph.frames(), tuple)
    assert isinstance(graph.transforms(), tuple)
    assert [frame.name for frame in graph.frames()] == sorted(
        frame.name for frame in graph.frames()
    )
    assert len(graph.transforms()) == 2


def test_graph_convert_accepts_frame_objects_and_vector_sets() -> None:
    from pytex.core.batches import VectorSet

    graph = rolling_frame_graph(rd_offset_deg=90.0)
    vectors = VectorSet(values=[[1.0, 0.0, 0.0]], reference_frame=SPECIMEN_FRAME)
    converted = graph.convert(
        vectors,
        source=SPECIMEN_FRAME,
        target=graph.frame("sample_rd_td_nd"),
        directions=True,
    )
    assert isinstance(converted, VectorSet)
    assert converted.reference_frame == graph.frame("sample_rd_td_nd")
    assert_allclose(converted.values[0], [0.0, -1.0, 0.0], atol=1e-12)


def test_graph_describe_names_frames_and_relationships() -> None:
    text = rolling_frame_graph(rd_offset_deg=15.0).describe()
    assert "rolling_geometry" in text
    assert "sample_rd_td_nd" in text
    assert "specimen" in text
    assert "shortest declared chain" in text


def test_empty_graph_describe_says_so() -> None:
    assert "is empty" in FrameGraph(name="blank").describe()


def test_graph_without_relationships_says_frames_are_unconnected() -> None:
    graph = FrameGraph(frames=[specimen_frame(), detector_frame()], name="loose")
    assert "No relationships are declared" in graph.describe()


def test_graph_accepts_transforms_at_construction() -> None:
    transform = FrameTransform.from_axis_correspondence(
        SPECIMEN_FRAME,
        SAMPLE_RD_TD_ND_FRAME,
        {"x": "RD", "y": "TD", "z": "ND"},
    )
    graph = FrameGraph(transforms=[transform], name="declared")
    assert len(graph) == 2
    assert graph.transform_between("specimen", "sample_rd_td_nd").is_identity


def test_laboratory_frame_builder_is_usable_in_a_graph() -> None:
    graph = FrameGraph(name="diffraction")
    graph.add_transform(
        FrameTransform.from_axis_angle(
            (0.0, 0.0, 1.0),
            12.0,
            source=laboratory_frame(),
            target=detector_frame(),
        )
    )
    assert graph.transform_between("laboratory", "detector").rotation_angle_deg == pytest.approx(
        12.0
    )
