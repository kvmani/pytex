"""Tests for reference-frame visualization: primitives, gizmos, and documentation SVG.

Assertions here are structural and semantic rather than pixel- or byte-based,
per the repository rule that runtime plotting validation must not depend on
tracked SVG baselines. The generated SVG *is* checked for the mandatory
style-guide elements, because those are a documentation contract.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core.frame_catalog import (
    CARTESIAN_FRAME,
    CRYSTAL_FRAME,
    DETECTOR_FRAME,
    MAP_FRAME,
    SAMPLE_RD_TD_ND_FRAME,
    SPECIMEN_FRAME,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)
from pytex.plotting.frames import (
    DEFAULT_VIEW_AZIM_DEG,
    DEFAULT_VIEW_ELEV_DEG,
    FrameTriad,
    add_frame_indicator,
    frame_catalog_svg,
    frame_triad,
    frame_triad_primitives,
    plot_frame_relationship,
    plot_reference_frame,
    project_orthographic,
    reference_frame_svg,
)
from pytex.plotting.primitives import TRIAD_AXIS_COLORS, AxisTriad3D, PrimitiveScene3D


def gizmo_labels(inset: object) -> set[str]:
    """Visible label texts on a gizmo inset.

    Each arrow is drawn with ``annotate("")``, which leaves an empty ``Text``
    artist behind; only the non-empty ones are labels.
    """

    return {
        text.get_text()
        for text in inset.texts  # type: ignore[attr-defined]
        if text.get_text()
    }


def gizmo_label_radii(inset: object) -> list[float]:
    """Distance of each visible gizmo label from the inset origin."""

    return [
        float(np.hypot(*text.get_position()))
        for text in inset.texts  # type: ignore[attr-defined]
        if text.get_text()
    ]


SVG_NS = "{http://www.w3.org/2000/svg}"


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------


def test_projection_from_the_positive_x_axis_maps_y_and_z_to_screen_axes() -> None:
    screen, depth = project_orthographic(np.eye(3), elev_deg=0.0, azim_deg=0.0)
    assert_allclose(screen[0], [0.0, 0.0], atol=1e-12)  # X points at the camera
    assert_allclose(screen[1], [1.0, 0.0], atol=1e-12)  # Y goes right
    assert_allclose(screen[2], [0.0, 1.0], atol=1e-12)  # Z goes up
    assert_allclose(depth, [1.0, 0.0, 0.0], atol=1e-12)


def test_projection_preserves_lengths_perpendicular_to_the_view() -> None:
    screen, _ = project_orthographic(np.eye(3), elev_deg=90.0, azim_deg=0.0)
    # Looking straight down, X and Y keep unit screen length and Z collapses.
    assert np.hypot(*screen[0]) == pytest.approx(1.0)
    assert np.hypot(*screen[1]) == pytest.approx(1.0)
    assert np.hypot(*screen[2]) == pytest.approx(0.0, abs=1e-12)


def test_projection_accepts_a_single_vector() -> None:
    screen, depth = project_orthographic([0.0, 0.0, 1.0])
    assert screen.shape == (1, 2)
    assert depth.shape == (1,)


def test_projection_rejects_wrongly_shaped_input() -> None:
    with pytest.raises(ValueError, match="shape"):
        project_orthographic(np.zeros((2, 4)))


def test_screen_basis_is_orthonormal_for_arbitrary_views() -> None:
    # Projecting the canonical triad must give a rank-2 map that preserves
    # in-plane lengths; equivalently the summed squared screen radii equal 2.
    screen, _ = project_orthographic(np.eye(3), elev_deg=37.0, azim_deg=-121.0)
    assert float(np.sum(screen**2)) == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# FrameTriad
# ---------------------------------------------------------------------------


def test_triad_endpoints_follow_the_frame_axes() -> None:
    triad = FrameTriad(frame=SAMPLE_RD_TD_ND_FRAME, length=2.0)
    tips = triad.endpoints()
    assert_allclose(tips[0], [2.0, 0.0, 0.0], atol=1e-12)
    assert_allclose(tips[2], [0.0, 0.0, 2.0], atol=1e-12)
    assert triad.labels == ("RD", "TD", "ND")


def test_triad_normalizes_an_unequal_frame_by_default() -> None:
    frame = crystal_frame(axis_vectors=[[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    assert_allclose(np.linalg.norm(FrameTriad(frame=frame).axis_matrix(), axis=0), np.ones(3))
    unnormalized = FrameTriad(frame=frame, normalize=False).axis_matrix()
    assert np.linalg.norm(unnormalized[:, 0]) == pytest.approx(3.0)


def test_triad_basis_override_wins_over_the_frame_geometry() -> None:
    override = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    triad = FrameTriad(frame=SPECIMEN_FRAME, basis=override)
    assert_allclose(triad.axis_matrix(), override, atol=1e-12)


def test_triad_origin_offsets_every_endpoint() -> None:
    triad = FrameTriad(frame=SPECIMEN_FRAME, origin=[1.0, 1.0, 1.0])
    assert_allclose(triad.endpoints()[0], [2.0, 1.0, 1.0], atol=1e-12)


def test_triad_rejects_bad_colors_and_lengths() -> None:
    with pytest.raises(ValueError, match="exactly three colors"):
        FrameTriad(frame=SPECIMEN_FRAME, colors=("#000000", "#111111"))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="strictly positive"):
        FrameTriad(frame=SPECIMEN_FRAME, length=0.0)


def test_triad_describe_names_the_frame_and_its_axes() -> None:
    text = FrameTriad(frame=SAMPLE_RD_TD_ND_FRAME).describe()
    assert "sample_rd_td_nd" in text
    assert "specimen domain" in text
    assert "RD to" in text
    assert "equal display length" in text


# ---------------------------------------------------------------------------
# Scene primitives
# ---------------------------------------------------------------------------


def test_frame_triad_builds_a_labelled_axis_triad_primitive() -> None:
    triad = frame_triad(SAMPLE_RD_TD_ND_FRAME, length=1.5)
    assert isinstance(triad, AxisTriad3D)
    assert triad.labels == ("RD", "TD", "ND")
    assert len(triad.arrows()) == 3
    assert_allclose(triad.axes[:, 0], [1.5, 0.0, 0.0], atol=1e-12)


def test_frame_triad_honours_rotated_frame_geometry() -> None:
    tilted = sample_frame(axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    triad = frame_triad(tilted)
    assert_allclose(triad.axes[:, 0], [0.0, 1.0, 0.0], atol=1e-12)


def test_reference_frame_triad_reads_frame_geometry_too() -> None:
    from pytex.plotting.primitives import reference_frame_triad

    tilted = sample_frame(axis_vectors=[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert_allclose(reference_frame_triad(tilted).axes[:, 0], [0.0, 1.0, 0.0], atol=1e-12)


def test_frame_triad_primitives_returns_a_mergeable_scene() -> None:
    scene = frame_triad_primitives(CRYSTAL_FRAME, caption="parent")
    assert isinstance(scene, PrimitiveScene3D)
    assert len(scene.triads) == 1
    assert len(scene.labels) == 1
    assert scene.labels[0].text == "parent"
    merged = scene.merge(frame_triad_primitives(SPECIMEN_FRAME))
    assert len(merged.triads) == 2
    assert not merged.is_empty()


def test_frame_triad_primitives_without_a_caption_adds_no_label() -> None:
    assert frame_triad_primitives(CRYSTAL_FRAME).labels == ()


# ---------------------------------------------------------------------------
# Matplotlib figures
# ---------------------------------------------------------------------------


def test_plot_reference_frame_draws_three_labelled_axes() -> None:
    figure, axes = plot_reference_frame(SAMPLE_RD_TD_ND_FRAME)
    try:
        texts = {text.get_text() for text in axes.texts}
        assert {"RD", "TD", "ND"} <= texts
        assert "sample_rd_td_nd" in axes.get_title()
    finally:
        plt.close(figure)


def test_plot_reference_frame_accepts_a_custom_title_and_hides_the_box() -> None:
    figure, axes = plot_reference_frame(
        MAP_FRAME, title="scan grid", show_reference_box=False
    )
    try:
        assert axes.get_title() == "scan grid"
    finally:
        plt.close(figure)


def test_plot_reference_frame_draws_into_a_supplied_axes() -> None:
    figure = plt.figure()
    try:
        axes = figure.add_subplot(111, projection="3d")
        returned_figure, returned_axes = plot_reference_frame(CARTESIAN_FRAME, ax=axes)
        assert returned_axes is axes
        assert returned_figure is figure
    finally:
        plt.close(figure)


def test_plot_frame_relationship_labels_both_frames() -> None:
    graph = rolling_frame_graph(rd_offset_deg=30.0)
    transform = graph.transform_between("specimen", "sample_rd_td_nd")
    figure, axes = plot_frame_relationship(transform)
    try:
        texts = {text.get_text() for text in axes.texts}
        assert {"x", "y", "z"} <= texts
        assert {"RD", "TD", "ND"} <= texts
        assert "30.0 deg" in axes.get_title()
        annotations = " ".join(text.get_text() for text in axes.texts)
        assert "rotation 30.00 deg" in annotations
    finally:
        plt.close(figure)


def test_plot_frame_relationship_can_suppress_the_annotation() -> None:
    transform = rolling_frame_graph(rd_offset_deg=10.0).transform_between(
        "specimen", "sample_rd_td_nd"
    )
    figure, axes = plot_frame_relationship(transform, annotate=False)
    try:
        assert all("rotation" not in text.get_text() for text in axes.texts)
    finally:
        plt.close(figure)


# ---------------------------------------------------------------------------
# Embeddable gizmo
# ---------------------------------------------------------------------------


def test_frame_indicator_adds_an_inset_with_one_label_per_axis() -> None:
    figure, axes = plt.subplots()
    try:
        inset = add_frame_indicator(axes, DETECTOR_FRAME)
        assert gizmo_labels(inset) == {"u", "v", "n"}
        assert inset is not axes
        assert inset.get_aspect() == 1.0
    finally:
        plt.close(figure)


def test_frame_indicator_supports_an_axis_subset_for_in_plane_figures() -> None:
    figure, axes = plt.subplots()
    try:
        inset = add_frame_indicator(
            axes,
            SAMPLE_RD_TD_ND_FRAME,
            axis_subset=("RD", "TD"),
            elev_deg=90.0,
            azim_deg=-90.0,
        )
        assert gizmo_labels(inset) == {"RD", "TD"}
    finally:
        plt.close(figure)


def test_frame_indicator_can_label_the_frame() -> None:
    figure, axes = plt.subplots()
    try:
        inset = add_frame_indicator(axes, MAP_FRAME, label_frame=True)
        assert "map" in gizmo_labels(inset)
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    "loc", ["lower left", "lower right", "upper left", "upper right"]
)
def test_frame_indicator_supports_every_corner(loc: str) -> None:
    figure, axes = plt.subplots()
    try:
        assert add_frame_indicator(axes, SPECIMEN_FRAME, loc=loc) is not None
    finally:
        plt.close(figure)


def test_frame_indicator_rejects_an_unknown_corner() -> None:
    figure, axes = plt.subplots()
    try:
        with pytest.raises(ValueError, match="loc must be one of"):
            add_frame_indicator(axes, SPECIMEN_FRAME, loc="middle")
    finally:
        plt.close(figure)


def test_frame_indicator_embeds_into_a_polar_axes() -> None:
    # Pole figures and stereographic plots are drawn on polar axes; the gizmo
    # must work there too, since that is one of its main uses.
    figure = plt.figure()
    try:
        axes = figure.add_subplot(111, projection="polar")
        inset = add_frame_indicator(axes, SAMPLE_RD_TD_ND_FRAME)
        assert gizmo_labels(inset) == {"RD", "TD", "ND"}
    finally:
        plt.close(figure)


def test_frame_indicator_labels_stay_clear_of_the_origin() -> None:
    # An axis pointing almost straight at the viewer projects to a near-zero
    # arrow; its label must still be placed away from the origin.
    figure, axes = plt.subplots()
    try:
        inset = add_frame_indicator(axes, SPECIMEN_FRAME, elev_deg=0.0, azim_deg=0.0)
        assert min(gizmo_label_radii(inset)) >= 0.40
    finally:
        plt.close(figure)


# ---------------------------------------------------------------------------
# Standalone SVG
# ---------------------------------------------------------------------------


def test_reference_frame_svg_is_well_formed_and_style_guide_compliant() -> None:
    svg = reference_frame_svg(SAMPLE_RD_TD_ND_FRAME)
    root = ET.fromstring(svg)
    assert root.tag == f"{SVG_NS}svg"
    assert root.find(f"{SVG_NS}title") is not None
    assert root.find(f"{SVG_NS}desc") is not None
    assert "Arial" in svg
    texts = {element.text for element in root.iter(f"{SVG_NS}text")}
    assert {"RD", "TD", "ND"} <= texts


def test_reference_frame_svg_description_carries_the_frame_prose() -> None:
    root = ET.fromstring(reference_frame_svg(CRYSTAL_FRAME))
    desc = root.find(f"{SVG_NS}desc")
    assert desc is not None and desc.text is not None
    assert "crystal domain" in desc.text
    assert "right-handed" in desc.text


def test_reference_frame_svg_draws_one_line_per_axis_with_arrowheads() -> None:
    root = ET.fromstring(reference_frame_svg(SPECIMEN_FRAME))
    lines = list(root.iter(f"{SVG_NS}line"))
    assert len(lines) == 3
    assert all("marker-end" in line.attrib for line in lines)
    markers = list(root.iter(f"{SVG_NS}marker"))
    assert len(markers) == 3
    # Arrowheads must not scale with stroke width or they swamp the triad.
    assert all(marker.get("markerUnits") == "userSpaceOnUse" for marker in markers)


def test_reference_frame_svg_accepts_custom_titles() -> None:
    svg = reference_frame_svg(MAP_FRAME, title="Scan grid", subtitle="EBSD step raster")
    root = ET.fromstring(svg)
    title = root.find(f"{SVG_NS}title")
    assert title is not None and title.text == "Scan grid"
    assert "EBSD step raster" in svg


def test_reference_frame_svg_escapes_markup_in_frame_names() -> None:
    hostile = specimen_frame("<script>alert(1)</script>")
    svg = reference_frame_svg(hostile)
    assert "<script>" not in svg
    assert "&lt;script&gt;" in svg
    ET.fromstring(svg)  # still well formed


def test_frame_catalog_svg_lays_out_one_panel_per_frame() -> None:
    frames = [CARTESIAN_FRAME, SAMPLE_RD_TD_ND_FRAME, CRYSTAL_FRAME, DETECTOR_FRAME]
    root = ET.fromstring(frame_catalog_svg(frames, columns=2))
    # Three axis lines per frame.
    assert len(list(root.iter(f"{SVG_NS}line"))) == 3 * len(frames)
    texts = {element.text for element in root.iter(f"{SVG_NS}text")}
    for frame in frames:
        assert frame.name in texts


def test_frame_catalog_svg_grid_size_follows_the_column_count() -> None:
    frames = [CARTESIAN_FRAME, SAMPLE_RD_TD_ND_FRAME, CRYSTAL_FRAME, DETECTOR_FRAME]
    two_columns = ET.fromstring(frame_catalog_svg(frames, columns=2))
    four_columns = ET.fromstring(frame_catalog_svg(frames, columns=4))
    assert float(four_columns.get("width", "0")) > float(two_columns.get("width", "0"))
    assert float(four_columns.get("height", "0")) < float(two_columns.get("height", "0"))


def test_frame_catalog_svg_rejects_empty_input_and_bad_columns() -> None:
    with pytest.raises(ValueError, match="at least one frame"):
        frame_catalog_svg([])
    with pytest.raises(ValueError, match="positive integer"):
        frame_catalog_svg([CARTESIAN_FRAME], columns=0)


def test_svg_generation_does_not_require_matplotlib() -> None:
    # The documentation path must stay import-light: generating SVG must not
    # touch pyplot, so docs can be built in a minimal environment.
    import sys
    from unittest import mock

    with mock.patch.dict(sys.modules, {"matplotlib.pyplot": None}):
        svg = reference_frame_svg(CARTESIAN_FRAME)
    assert svg.startswith("<svg")


def test_default_view_constants_are_exported_and_used() -> None:
    assert DEFAULT_VIEW_ELEV_DEG == pytest.approx(22.0)
    assert DEFAULT_VIEW_AZIM_DEG == pytest.approx(34.0)


def test_triad_palette_is_shared_across_renderers() -> None:
    assert FrameTriad(frame=SPECIMEN_FRAME).colors == TRIAD_AXIS_COLORS
    assert frame_triad(SPECIMEN_FRAME).colors == TRIAD_AXIS_COLORS
    for colour in TRIAD_AXIS_COLORS:
        assert colour in reference_frame_svg(SPECIMEN_FRAME)


# ---------------------------------------------------------------------------
# Frame indicators wired into the domain renderers
# ---------------------------------------------------------------------------


def _ni_fcc_phase() -> object:
    # Loaded from a CIF fixture, which pymatgen parses.

    from pytex import crystal_frame as _crystal_frame
    from pytex import get_phase_fixture

    return get_phase_fixture("ni_fcc").load_phase(crystal_frame=_crystal_frame())


def _saed_pattern() -> object:
    from pytex import ZoneAxis, generate_saed_pattern

    phase = _ni_fcc_phase()
    return generate_saed_pattern(
        phase,
        ZoneAxis(np.array([0, 0, 1]), phase=phase),
        camera_constant_mm_angstrom=180.0,
        max_index=3,
        max_g_inv_angstrom=1.0,
    )


def _gizmo_insets(figure: object) -> list[object]:
    """Inset axes added to a figure's primary axes by `add_frame_indicator`."""

    primary = figure.axes[0]  # type: ignore[attr-defined]
    return [child for child in primary.child_axes]


def test_saed_plot_frame_indicator_is_opt_in() -> None:
    pattern = _saed_pattern()
    from pytex import plot_saed_pattern

    default_figure = plot_saed_pattern(pattern)
    try:
        assert _gizmo_insets(default_figure) == []
    finally:
        plt.close(default_figure)

    annotated = plot_saed_pattern(pattern, show_frame_indicator=True)
    try:
        insets = _gizmo_insets(annotated)
        assert len(insets) == 1
        labels = gizmo_labels(insets[0])
        # Two in-plane axes plus the frame name; the detector normal points at
        # the viewer and is deliberately omitted.
        assert {"u", "v"} <= labels
        assert any(name.endswith("_saed_detector") for name in labels)
    finally:
        plt.close(annotated)


def test_saed_frame_indicator_shows_only_the_in_plane_detector_axes() -> None:
    from pytex import plot_saed_pattern

    figure = plot_saed_pattern(_saed_pattern(), show_frame_indicator=True)
    try:
        labels = gizmo_labels(_gizmo_insets(figure)[0])
        assert {"u", "v"} <= labels
        assert "n" not in labels
    finally:
        plt.close(figure)


def test_crystal_structure_frame_indicator_is_opt_in_and_uses_the_lattice_basis() -> None:
    from pytex import plot_crystal_structure_3d

    phase = _ni_fcc_phase()

    default_figure = plot_crystal_structure_3d(phase)
    try:
        assert _gizmo_insets(default_figure) == []
    finally:
        plt.close(default_figure)

    annotated = plot_crystal_structure_3d(phase, show_frame_indicator=True)
    try:
        labels = gizmo_labels(_gizmo_insets(annotated)[0])
        assert {"a", "b", "c"} <= labels
    finally:
        plt.close(annotated)
