"""The orientation-relationship stereogram (F18): pairs, tie-lines, circles."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    plot_or_stereogram,
)
from pytex.core.transformation import OrientationRelationship
from pytex.plotting.spherical import (
    ORStereogramPair,
    build_or_stereogram_figure_spec,
    or_stereogram_pairs,
)


def _cubic(name: str, a: float) -> Phase:
    frame = ReferenceFrame(name, FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(a, a, a, 90.0, 90.0, 90.0, crystal_frame=frame)
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
    )


def _ks() -> OrientationRelationship:
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=_cubic("austenite", 3.6), child_phase=_cubic("ferrite", 2.87)
    )


def _nw() -> OrientationRelationship:
    return OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=_cubic("austenite", 3.6), child_phase=_cubic("ferrite", 2.87)
    )


# --------------------------------------------------------------------------- #
# The pairs
# --------------------------------------------------------------------------- #


def test_defining_pairs_coincide_on_the_net_for_every_variant() -> None:
    relationship = _ks()
    for variant in range(1, 25):
        pairs = or_stereogram_pairs(relationship, variant=variant)
        assert len(pairs) == 2  # one plane pair, one direction pair
        for pair in pairs:
            assert pair.deviation_deg == pytest.approx(0.0, abs=1e-9)
            np.testing.assert_allclose(pair.parent_vector, pair.child_vector, atol=1e-9)


def test_variant_none_is_variant_one_which_is_the_relationship_as_stated() -> None:
    relationship = _ks()
    assert or_stereogram_pairs(relationship, variant=None) == or_stereogram_pairs(
        relationship, variant=1
    )
    nominal_plane = relationship.parallel_planes[0][0]
    pair = or_stereogram_pairs(relationship, variant=1)[0]
    np.testing.assert_allclose(pair.parent_vector, nominal_plane.normal, atol=1e-12)


def test_labels_follow_the_variant_rather_than_repeating_variant_ones() -> None:
    relationship = _ks()
    labels = {or_stereogram_pairs(relationship, variant=k)[0].label for k in range(1, 25)}
    # the 24 variants name the four {111} parent members, not one
    assert len(labels) == 4
    assert "(111) ∥ (011)" in labels


def test_nominated_family_adds_non_zero_deviation_tie_lines() -> None:
    relationship = _nw()
    parent = relationship.parent_phase
    pairs = or_stereogram_pairs(
        relationship,
        parent_planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        tolerance_deg=6.0,
    )
    deviations = sorted(pair.deviation_deg for pair in pairs)
    assert deviations[0] == pytest.approx(0.0, abs=1e-9)
    assert deviations[-1] > 0.1  # the non-parallel {111} members carry real angles
    assert all(deviation <= 6.0 for deviation in deviations)


def test_nominated_family_does_not_duplicate_the_defining_pair() -> None:
    relationship = _nw()
    parent = relationship.parent_phase
    pairs = or_stereogram_pairs(
        relationship,
        parent_planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        tolerance_deg=6.0,
    )
    exact = [pair for pair in pairs if pair.deviation_deg < 1e-9 and pair.kind == "plane"]
    assert len(exact) == 1


def test_nominated_directions_are_accepted_too() -> None:
    relationship = _ks()
    parent = relationship.parent_phase
    pairs = or_stereogram_pairs(
        relationship,
        parent_directions=[CrystalDirection([1.0, 0.0, -1.0], phase=parent)],
        tolerance_deg=3.0,
    )
    assert any(pair.kind == "direction" for pair in pairs)
    assert all(pair.deviation_deg <= 3.0 for pair in pairs)


def test_pair_rejects_bad_kind_and_degenerate_vectors() -> None:
    with pytest.raises(ValueError, match="kind must be"):
        ORStereogramPair(
            kind="bogus",
            label="x",
            parent_vector=np.array([0.0, 0.0, 1.0]),
            child_vector=np.array([0.0, 0.0, 1.0]),
            deviation_deg=0.0,
        )
    with pytest.raises(ValueError, match="non-zero finite"):
        ORStereogramPair(
            kind="plane",
            label="x",
            parent_vector=np.zeros(3),
            child_vector=np.array([0.0, 0.0, 1.0]),
            deviation_deg=0.0,
        )
    with pytest.raises(ValueError, match="deviation_deg"):
        ORStereogramPair(
            kind="plane",
            label="x",
            parent_vector=np.array([0.0, 0.0, 1.0]),
            child_vector=np.array([0.0, 0.0, 1.0]),
            deviation_deg=-1.0,
        )


# --------------------------------------------------------------------------- #
# The figure
# --------------------------------------------------------------------------- #


def test_figure_spec_carries_open_parent_and_filled_child_markers() -> None:
    spec = build_or_stereogram_figure_spec(_ks())
    assert len(spec.marker_layers) == 2
    parent_layer, child_layer = spec.marker_layers
    assert parent_layer.label == "parent"
    assert list(parent_layer.facecolors) == ["none"] * parent_layer.points.shape[0]
    assert child_layer.label == "child"
    assert parent_layer.points.shape == child_layer.points.shape
    # a parallelism plots as two symbols on top of each other
    np.testing.assert_allclose(parent_layer.points, child_layer.points, atol=1e-9)


def test_plane_pairs_contribute_two_great_circles_each() -> None:
    relationship = _ks()
    with_circles = build_or_stereogram_figure_spec(relationship)
    without = build_or_stereogram_figure_spec(relationship, show_great_circles=False)
    plane_pairs = sum(
        1 for pair in or_stereogram_pairs(relationship) if pair.kind == "plane"
    )
    # at least two -- a trace that crosses the rim is split into segments
    assert len(with_circles.line_layers) - len(without.line_layers) >= 2 * plane_pairs


def test_tie_lines_are_drawn_and_can_be_switched_off() -> None:
    relationship = _nw()
    parent = relationship.parent_phase
    planes = [CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)]
    with_ties = build_or_stereogram_figure_spec(
        relationship, parent_planes=planes, tolerance_deg=6.0
    )
    without = build_or_stereogram_figure_spec(
        relationship, parent_planes=planes, tolerance_deg=6.0, show_tie_lines=False
    )
    assert len(with_ties.line_layers) > len(without.line_layers)


def test_labels_state_the_indices_and_the_deviation() -> None:
    spec = build_or_stereogram_figure_spec(_ks())
    texts = [layer.text for layer in spec.text_layers]
    assert any(text.startswith("(111) ∥ (011)") and text.endswith("0.00 deg") for text in texts)
    bare = build_or_stereogram_figure_spec(_ks(), label_pairs=False)
    assert bare.text_layers == ()


def test_every_drawn_point_lies_inside_the_projection_boundary() -> None:
    relationship = _nw()
    parent = relationship.parent_phase
    spec = build_or_stereogram_figure_spec(
        relationship,
        parent_planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        tolerance_deg=6.0,
    )
    limit = spec.boundary_circle_radius + 1e-6
    for layer in spec.marker_layers:
        assert np.all(np.linalg.norm(layer.points, axis=1) <= limit)
    for layer in spec.line_layers:
        assert np.all(np.linalg.norm(layer.points, axis=1) <= limit)


def test_equal_area_net_is_available_and_rescales_the_boundary() -> None:
    stereographic = build_or_stereogram_figure_spec(_ks(), method="stereographic")
    equal_area = build_or_stereogram_figure_spec(_ks(), method="equal_area")
    assert equal_area.boundary_circle_radius > stereographic.boundary_circle_radius


def test_relationship_without_parallelisms_refuses_to_draw_nothing() -> None:
    parent, child = _cubic("austenite", 3.6), _cubic("ferrite", 2.87)
    bare = OrientationRelationship(
        name="bare",
        parent_phase=parent,
        child_phase=child,
        parent_to_child_rotation=_ks().parent_to_child_rotation,
    )
    with pytest.raises(ValueError, match="nothing to draw"):
        build_or_stereogram_figure_spec(bare)


def test_plot_or_stereogram_renders() -> None:
    figure = plot_or_stereogram(_ks(), variant=17, title="KS variant 17")
    assert figure.axes[0].get_title() == "KS variant 17"
    plt.close(figure)


def test_every_variants_parallelism_plots_as_one_point_and_one_circle() -> None:
    # The defect this pins: variants 7 and 9 of Kurdjumov-Sachs have a defining
    # direction lying in the equatorial plane, and the variant rotation returns
    # the child copy at z = -8e-16. Folded independently, the two ends of a
    # zero-deviation tie-line land on opposite rims -- a diameter across a
    # figure whose entire claim is that the two poles coincide.
    relationship = _ks()
    worst_pole = 0.0
    worst_circle = 0.0
    for variant in range(1, 25):
        spec = build_or_stereogram_figure_spec(
            relationship, variant=variant, include_wulff_net=False, show_tie_lines=False
        )
        parent_layer, child_layer = spec.marker_layers
        worst_pole = max(
            worst_pole,
            float(np.max(np.linalg.norm(parent_layer.points - child_layer.points, axis=1))),
        )
        parent_circle, child_circle = (
            spec.line_layers[0].points,
            spec.line_layers[1].points,
        )
        worst_circle = max(
            worst_circle,
            float(np.max(np.linalg.norm(parent_circle - child_circle, axis=1))),
        )
    assert worst_pole < 1e-9
    assert worst_circle < 1e-9


def test_a_real_straddle_splits_the_tie_line_rather_than_crossing_the_net() -> None:
    relationship = _nw()
    parent = relationship.parent_phase
    spec = build_or_stereogram_figure_spec(
        relationship,
        parent_planes=[CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=parent), phase=parent)],
        tolerance_deg=6.0,
        include_wulff_net=False,
    )
    # no segment this figure draws may span more than the disc radius in one
    # step; the Wulff net's own meridians are excluded because their rim jumps
    # predate this figure and belong to `generate_stereonet_grid`
    for layer in spec.line_layers:
        steps = np.linalg.norm(np.diff(layer.points, axis=0), axis=1)
        assert np.all(steps <= spec.boundary_circle_radius + 1e-9)
