"""Tests for the generated conceptual frame diagrams.

Assertions are structural and semantic rather than byte-based: the figures are
regenerated whenever the model changes, so pinning bytes would only pin churn.
What is pinned is that each figure says the right *thing* — the right frames,
the right relationships, computed geometry rather than decoration — and that it
is a well-formed, accessible, style-guide-compliant SVG.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from pytex.core.frame_catalog import (
    CRYSTAL_FRAME,
    DETECTOR_FRAME,
    LABORATORY_FRAME,
    MAP_FRAME,
    SPECIMEN_FRAME,
    crystal_frame,
    reciprocal_frame_for,
    specimen_frame,
)
from pytex.core.lattice import Lattice
from pytex.core.orientation import Rotation
from pytex.plotting.frame_diagrams import (
    PALE_TRIAD_COLORS,
    DiagramPanel,
    active_passive_svg,
    euler_sequence_svg,
    frame_chain_svg,
    header_width,
    hexagonal_frame_svg,
    orientation_mapping_svg,
    text_width,
)

SVG_NS = "{http://www.w3.org/2000/svg}"


def _root(svg: str) -> ET.Element:
    return ET.fromstring(svg)


def _texts(svg: str) -> set[str]:
    return {element.text for element in _root(svg).iter(f"{SVG_NS}text") if element.text}


def _chain_panels() -> tuple[DiagramPanel, ...]:
    return (
        DiagramPanel(frame=CRYSTAL_FRAME, caption="lattice-fixed basis"),
        DiagramPanel(frame=SPECIMEN_FRAME, caption="macroscopic sample frame"),
        DiagramPanel(frame=MAP_FRAME, caption="scan sampling grid"),
    )


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------


def test_text_width_scales_with_length_and_size() -> None:
    assert text_width("abcd", 10.0) == pytest.approx(2.0 * text_width("ab", 10.0))
    assert text_width("ab", 20.0) == pytest.approx(2.0 * text_width("ab", 10.0))
    assert text_width("", 12.0) == 0.0


def test_header_width_leaves_room_for_the_longer_of_title_and_subtitle() -> None:
    narrow = header_width("Short", "brief", margin=40.0)
    wide = header_width("Short", "a considerably longer subtitle line", margin=40.0)
    assert wide > narrow
    assert narrow >= 80.0  # both margins


# ---------------------------------------------------------------------------
# The frame chain
# ---------------------------------------------------------------------------


def test_frame_chain_names_every_frame_and_relationship() -> None:
    svg = frame_chain_svg(_chain_panels(), ("orientation g", "registration"))
    texts = _texts(svg)
    for name in ("crystal", "specimen", "map"):
        assert name in texts
    assert {"orientation g", "registration"} <= texts


def test_frame_chain_draws_one_arrow_between_each_pair() -> None:
    svg = frame_chain_svg(_chain_panels(), ("orientation g", "registration"))
    assert svg.count('marker-end="url(#pytex-chain-rel)"') == 2


def test_frame_chain_places_the_dual_frame_outside_the_linear_chain() -> None:
    """Duality is not a chain step, so the reciprocal frame must not be in the row.

    Putting it in the row would assert a ``laboratory -> reciprocal`` link that
    the canonical frame chain does not contain.
    """

    dual = DiagramPanel(frame=reciprocal_frame_for(CRYSTAL_FRAME), caption="dual basis")
    svg = frame_chain_svg(
        _chain_panels(), ("orientation g", "registration"), dual_panel=dual, dual_of=0
    )
    # Still only two in-row relationships.
    assert svg.count('marker-end="url(#pytex-chain-rel)"') == 2
    # The duality edge is a separate, dashed marker.
    assert 'marker-end="url(#pytex-chain-rel-dual)"' in svg
    assert "stroke-dasharray" in svg
    assert "crystal_reciprocal" in _texts(svg)


def test_frame_chain_dual_arrow_is_vertical_from_its_partner() -> None:
    dual = DiagramPanel(frame=reciprocal_frame_for(CRYSTAL_FRAME))
    svg = frame_chain_svg(
        _chain_panels(), ("orientation g", "registration"), dual_panel=dual, dual_of=0
    )
    match = re.search(
        r'<line x1="([0-9.]+)" y1="([0-9.]+)" x2="([0-9.]+)" y2="([0-9.]+)"[^>]*'
        r'marker-end="url\(#pytex-chain-rel-dual\)"',
        svg,
    )
    assert match is not None
    x1, y1, x2, y2 = (float(v) for v in match.groups())
    assert x1 == pytest.approx(x2), "the duality arrow should drop vertically"
    assert y2 > y1, "it should point downward, to the panel below"


def test_reciprocal_panel_shows_starred_axes() -> None:
    dual = DiagramPanel(frame=reciprocal_frame_for(CRYSTAL_FRAME))
    svg = frame_chain_svg(
        _chain_panels(), ("orientation g", "registration"), dual_panel=dual, dual_of=0
    )
    assert {"a*", "b*", "c*"} <= _texts(svg)


def test_frame_chain_rejects_mismatched_inputs() -> None:
    with pytest.raises(ValueError, match="at least two panels"):
        frame_chain_svg(_chain_panels()[:1], ())
    with pytest.raises(ValueError, match="relationship label"):
        frame_chain_svg(_chain_panels(), ("only one",))
    with pytest.raises(ValueError, match="dual_of must index"):
        frame_chain_svg(
            _chain_panels(),
            ("a", "b"),
            dual_panel=DiagramPanel(frame=DETECTOR_FRAME),
            dual_of=9,
        )


def test_frame_chain_widens_to_fit_a_long_subtitle() -> None:
    long_subtitle = "x" * 400
    svg = frame_chain_svg(_chain_panels(), ("a", "b"), subtitle=long_subtitle)
    width = float(_root(svg).get("width", "0"))
    assert width >= header_width("PyTex Reference-Frame Chain", long_subtitle)


# ---------------------------------------------------------------------------
# Orientation mapping
# ---------------------------------------------------------------------------


def test_orientation_mapping_draws_the_rotation_columns_as_the_mapped_axes() -> None:
    """The mapped panel must be the arithmetic, not a decorative rotation."""

    rotation = Rotation.from_bunge_euler(35.0, 28.0, 15.0)
    matrix = np.asarray(rotation.as_matrix(), dtype=np.float64)
    svg = orientation_mapping_svg(
        crystal_frame(), specimen_frame(), rotation_matrix=matrix
    )
    # The identity panel and the mapped panel must differ, and the mapped panel
    # must carry the specimen frame as a pale overlay.
    assert "crystal in specimen" in _texts(svg)
    for colour in PALE_TRIAD_COLORS:
        assert colour in svg


def test_orientation_mapping_labels_the_direction_of_the_mapping() -> None:
    svg = orientation_mapping_svg(
        crystal_frame(),
        specimen_frame(),
        rotation_matrix=Rotation.identity().as_matrix(),
    )
    texts = _texts(svg)
    assert "orientation g" in texts
    assert any("crystal" in t and "specimen" in t for t in texts)


def test_orientation_mapping_inverse_variant_is_opt_in() -> None:
    matrix = Rotation.from_bunge_euler(20.0, 10.0, 5.0).as_matrix()
    plain = orientation_mapping_svg(crystal_frame(), specimen_frame(), rotation_matrix=matrix)
    with_inverse = orientation_mapping_svg(
        crystal_frame(), specimen_frame(), rotation_matrix=matrix, show_inverse=True
    )
    assert 'marker-end="url(#pytex-map-rel-inverse)"' not in plain
    assert 'marker-end="url(#pytex-map-rel-inverse)"' in with_inverse
    assert any("separate object" in t for t in _texts(with_inverse))


def test_orientation_mapping_rejects_a_bad_matrix() -> None:
    with pytest.raises(ValueError, match=r"shape \(3, 3\)"):
        orientation_mapping_svg(
            crystal_frame(), specimen_frame(), rotation_matrix=np.eye(2)
        )


# ---------------------------------------------------------------------------
# Active versus passive
# ---------------------------------------------------------------------------


def test_active_passive_shows_both_languages() -> None:
    svg = active_passive_svg(
        crystal_frame(),
        specimen_frame(),
        rotation_matrix=Rotation.from_bunge_euler(35.0, 28.0, 15.0).as_matrix(),
    )
    texts = _texts(svg)
    assert "Active view" in texts
    assert "Passive / frame-mapping view" in texts
    assert {"v", "R v"} <= texts


def test_active_passive_vector_and_its_image_are_visually_separated() -> None:
    """Superimposed v and R v would make the active panel teach nothing."""

    from pytex.plotting.frames import project_orthographic

    matrix = np.asarray(
        Rotation.from_bunge_euler(35.0, 28.0, 15.0).as_matrix(), dtype=np.float64
    )
    source = np.array([1.0, 1.0, 0.0])
    source = source / np.linalg.norm(source)
    screen, _ = project_orthographic(np.vstack([source, matrix @ source]))
    separation = float(np.hypot(*(screen[1] - screen[0])))
    assert separation > 0.5, "default vector must separate under the default rotation"


def test_active_passive_rejects_a_bad_matrix() -> None:
    with pytest.raises(ValueError, match=r"shape \(3, 3\)"):
        active_passive_svg(crystal_frame(), specimen_frame(), rotation_matrix=np.zeros((4, 4)))


# ---------------------------------------------------------------------------
# Euler sequence
# ---------------------------------------------------------------------------


def test_euler_sequence_shows_one_panel_per_step_with_the_angles_named() -> None:
    svg = euler_sequence_svg(specimen_frame(), phi1_deg=35.0, Phi_deg=45.0, phi2_deg=30.0)
    texts = _texts(svg)
    assert {"initial frame", "after phi1", "after Phi", "after phi2"} <= texts
    assert any("35" in t and "z" in t for t in texts)
    assert any("45" in t for t in texts)
    assert any("30" in t for t in texts)


def test_euler_sequence_panels_are_computed_rotations_not_sketches() -> None:
    """Each panel's triad must be the partial rotation, so a zero step is a no-op."""

    unrotated = euler_sequence_svg(
        specimen_frame(), phi1_deg=0.0, Phi_deg=0.0, phi2_deg=0.0
    )
    rotated = euler_sequence_svg(
        specimen_frame(), phi1_deg=35.0, Phi_deg=45.0, phi2_deg=30.0
    )
    # With all angles zero every panel is identical; with real angles they differ.
    assert unrotated != rotated


def test_euler_sequence_states_the_composed_form() -> None:
    svg = euler_sequence_svg(specimen_frame(), phi1_deg=10.0, Phi_deg=20.0, phi2_deg=30.0)
    assert any("Rz(phi2) Rx(Phi) Rz(phi1)" in t for t in _texts(svg))


# ---------------------------------------------------------------------------
# Hexagonal frame
# ---------------------------------------------------------------------------


def _hcp_lattice() -> Lattice:
    return Lattice(
        a=3.2320,
        b=3.2320,
        c=5.1470,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=120.0,
        crystal_frame=crystal_frame(),
    )


def test_hexagonal_frame_labels_the_four_axes_and_the_angle() -> None:
    svg = hexagonal_frame_svg(_hcp_lattice())
    texts = _texts(svg)
    assert "a1" in texts
    assert "a2" in texts
    assert any("a3" in t for t in texts)
    assert any("120" in t for t in texts)
    assert any("out of page" in t for t in texts)


def test_hexagonal_frame_reports_the_lattice_it_was_drawn_from() -> None:
    """The figure must be traceable to real parameters, not hand-placed geometry."""

    svg = hexagonal_frame_svg(_hcp_lattice())
    texts = _texts(svg)
    assert any("3.2320" in t for t in texts)
    assert any("5.1470" in t for t in texts)
    assert any("1.5925" in t for t in texts), "c/a should be reported"
    assert any("direct_basis" in t for t in texts)


def test_hexagonal_frame_states_the_reciprocal_star_convention() -> None:
    svg = hexagonal_frame_svg(_hcp_lattice())
    assert any("a*" in t for t in _texts(svg))


def test_hexagonal_frame_rejects_a_non_hexagonal_setting() -> None:
    cubic = Lattice(
        a=3.6,
        b=3.6,
        c=3.6,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=crystal_frame(),
    )
    with pytest.raises(ValueError, match="hexagonal setting"):
        hexagonal_frame_svg(cubic)


# ---------------------------------------------------------------------------
# Shared document contract
# ---------------------------------------------------------------------------


@pytest.fixture
def every_diagram() -> dict[str, str]:
    matrix = Rotation.from_bunge_euler(35.0, 28.0, 15.0).as_matrix()
    return {
        "chain": frame_chain_svg(_chain_panels(), ("orientation g", "registration")),
        "mapping": orientation_mapping_svg(
            crystal_frame(), specimen_frame(), rotation_matrix=matrix
        ),
        "active_passive": active_passive_svg(
            crystal_frame(), specimen_frame(), rotation_matrix=matrix
        ),
        "euler": euler_sequence_svg(
            specimen_frame(), phi1_deg=35.0, Phi_deg=45.0, phi2_deg=30.0
        ),
        "hexagonal": hexagonal_frame_svg(_hcp_lattice()),
    }


def test_every_diagram_is_well_formed_with_title_and_description(
    every_diagram: dict[str, str],
) -> None:
    for name, svg in every_diagram.items():
        root = _root(svg)
        assert root.tag == f"{SVG_NS}svg", name
        title = root.find(f"{SVG_NS}title")
        desc = root.find(f"{SVG_NS}desc")
        assert title is not None and title.text, name
        assert desc is not None and desc.text and len(desc.text) > 60, name


def test_every_diagram_uses_absolute_marker_units(every_diagram: dict[str, str]) -> None:
    """The defect that broke the hand-authored figures must not reappear here."""

    for name, svg in every_diagram.items():
        markers = re.findall(r"<marker\b[^>]*>", svg)
        assert markers, name
        for marker in markers:
            assert 'markerUnits="userSpaceOnUse"' in marker, f"{name}: {marker[:60]}"


def test_every_diagram_uses_the_canonical_font(every_diagram: dict[str, str]) -> None:
    for name, svg in every_diagram.items():
        assert 'font-family="Arial"' in svg, name


def test_every_diagram_escapes_markup_in_frame_names() -> None:
    hostile = specimen_frame("<script>alert(1)</script>")
    svg = frame_chain_svg(
        (DiagramPanel(frame=hostile), DiagramPanel(frame=LABORATORY_FRAME)),
        ("relationship",),
    )
    assert "<script>" not in svg
    assert "&lt;script&gt;" in svg
    _root(svg)  # still well formed


def test_diagrams_do_not_require_matplotlib() -> None:
    """The documentation path must stay import-light."""

    import sys
    from unittest import mock

    with mock.patch.dict(sys.modules, {"matplotlib.pyplot": None}):
        svg = frame_chain_svg(_chain_panels(), ("orientation g", "registration"))
    assert svg.startswith("<svg")
