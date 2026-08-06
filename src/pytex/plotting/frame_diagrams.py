"""Generated conceptual SVG diagrams built from real reference-frame geometry.

`pytex.plotting.frames` draws *one* frame. This module composes several frames,
their declared relationships, and the prose that explains them into the
multi-panel teaching diagrams the documentation needs — the frame chain, the
orientation-mapping semantics, the hexagonal crystal frame.

Why generated rather than hand-authored
---------------------------------------

These figures were previously drawn by hand, and hand-authored SVG drifts from
the model it illustrates. It also gets details wrong in ways that are invisible
in source and glaring on screen: every one of the replaced figures declared its
arrowhead markers without ``markerUnits="userSpaceOnUse"``, so SVG scaled the
head by the line's stroke width and a 12-unit arrowhead rendered at 48 units,
swamping the triads it was supposed to annotate.

Generating them from `pytex.core.frame_catalog` frames and `pytex.core.lattice`
bases means the axis directions in a figure are the axis directions in the
model, and the drawing rules are fixed in one place.

Everything here is pure Python — **no matplotlib** — so documentation figures
can be produced in a minimal environment.

Layout
------

Text is placed using an Arial advance-width estimate (`text_width`), so panels
size themselves to their content and labels do not collide. All geometry is in
SVG user units with ``y`` growing downward; the triad projection comes from
`pytex.plotting.frames.project_orthographic`, so a triad here and a triad in a
catalog figure are drawn identically.

See also
--------
`pytex.plotting.frames` : the single-frame renderers this builds on.
`scripts/generate_reference_frame_figures.py` : writes these into `docs/figures/`.
`docs/standards/visualization_style_guide.md` : the canonical tokens.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice
from pytex.plotting.frames import (
    DEFAULT_VIEW_AZIM_DEG,
    DEFAULT_VIEW_ELEV_DEG,
    FrameTriad,
    _svg_arrow_markers,
    _svg_triad_group,
    project_orthographic,
)
from pytex.plotting.primitives import TRIAD_AXIS_COLORS
from pytex.plotting.svg_primitives import (
    ACCENT as _ACCENT,
)
from pytex.plotting.svg_primitives import (
    MUTED as _MUTED,
)
from pytex.plotting.svg_primitives import (
    PANEL as _PANEL,
)
from pytex.plotting.svg_primitives import (
    PANEL_STROKE as _PANEL_STROKE,
)
from pytex.plotting.svg_primitives import (
    VIOLET as _VIOLET,
)
from pytex.plotting.svg_primitives import (
    arrow_marker as _plain_marker,
)
from pytex.plotting.svg_primitives import (
    callout as _callout,
)
from pytex.plotting.svg_primitives import (
    card as _card,
)
from pytex.plotting.svg_primitives import (
    document as _document,
)
from pytex.plotting.svg_primitives import (
    header_width,
    text_width,
)
from pytex.plotting.svg_primitives import (
    relationship_arrow as _relationship_arrow,
)
from pytex.plotting.svg_primitives import (
    text as _text,
)

#: Desaturated triad colours for a reference frame drawn behind a mapped one.
PALE_TRIAD_COLORS: tuple[str, str, str] = ("#93b4f5", "#8fd6bd", "#f2a0a0")

__all__ = [
    "PALE_TRIAD_COLORS",
    "DiagramPanel",
    "active_passive_svg",
    "euler_sequence_svg",
    "frame_chain_svg",
    "header_width",
    "hexagonal_frame_svg",
    "orientation_mapping_svg",
    "text_width",
]

# Frame-domain accent colours, so a domain keeps one identity across figures.
_DOMAIN_TINTS: dict[str, tuple[str, str]] = {
    "crystal": ("#eef4ff", "#b9cdf2"),
    "specimen": ("#fdf6e8", "#e6d4a8"),
    "map": ("#fdeef1", "#eebfc8"),
    "detector": ("#eaf4fb", "#b3d5ea"),
    "laboratory": ("#eef7f0", "#b7d9c0"),
    "reciprocal": ("#f3eefc", "#cdbdee"),
}


@dataclass(frozen=True, slots=True)
class DiagramPanel:
    """One frame panel in a composed concept diagram.

    Parameters
    ----------
    frame:
        The frame whose triad the panel shows.
    caption:
        A short line under the triad saying what the frame is attached to.
    basis:
        Optional ``(3, 3)`` override whose columns are the axis vectors to draw,
        for panels showing a frame in another frame's coordinates.
    title:
        Optional panel heading, defaulting to the frame name. Set it when two
        panels show the same frame in different coordinates, so the headings
        stay distinguishable.
    subtitle:
        Optional second heading line, defaulting to the domain and axis labels.
    overlay:
        Optional reference triad drawn *underneath* the main one in pale
        colours — used to show a mapped frame against the frame it was mapped
        into, so the rotation is visible in one panel instead of two.
    """

    frame: ReferenceFrame
    caption: str = ""
    basis: np.ndarray | None = field(default=None)
    title: str | None = None
    subtitle: str | None = None
    overlay: DiagramPanel | None = None


def _panel(
    panel: DiagramPanel,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    marker_id: str,
    scale: float,
    colors: tuple[str, str, str] = TRIAD_AXIS_COLORS,
    elev_deg: float = DEFAULT_VIEW_ELEV_DEG,
    azim_deg: float = DEFAULT_VIEW_AZIM_DEG,
    title_size: float = 15.0,
    axis_font: float = 12.0,
) -> str:
    """Render one titled frame panel."""

    frame = panel.frame
    fill, stroke = _DOMAIN_TINTS.get(frame.domain.value, (_PANEL, _PANEL_STROKE))
    centre_x = x + width / 2.0
    # Offset below the panel mid-line so a tall axis label clears the subtitle.
    centre_y = y + height / 2.0 + 30.0
    triad = FrameTriad(frame=frame, length=1.0, colors=colors, basis=panel.basis)
    parts = [
        _card(x, y, width, height, fill=fill, stroke=stroke),
        _text(
            centre_x,
            y + 26.0,
            panel.title if panel.title is not None else frame.name,
            size=title_size,
            anchor="middle",
        ),
        _text(
            centre_x,
            y + 44.0,
            panel.subtitle
            if panel.subtitle is not None
            else f"{frame.domain.value} · {'/'.join(frame.axes)}",
            size=11.0,
            fill=_MUTED,
            anchor="middle",
        ),
    ]
    if panel.overlay is not None:
        # Reference triad first, so the mapped axes draw on top of it.
        reference = panel.overlay
        parts.append(
            _svg_triad_group(
                FrameTriad(
                    frame=reference.frame,
                    length=1.0,
                    colors=PALE_TRIAD_COLORS,
                    basis=reference.basis,
                ),
                centre=(centre_x, centre_y),
                scale=scale,
                elev_deg=elev_deg,
                azim_deg=azim_deg,
                fontsize=axis_font - 1.0,
                marker_id=f"{marker_id}-pale",
                label_scale=1.22,
            )
        )
    parts.append(
        _svg_triad_group(
            triad,
            centre=(centre_x, centre_y),
            scale=scale,
            elev_deg=elev_deg,
            azim_deg=azim_deg,
            fontsize=axis_font,
            marker_id=marker_id,
        )
    )
    if panel.caption:
        parts.append(
            _text(centre_x, y + height - 14.0, panel.caption, size=11.0,
                  fill=_MUTED, anchor="middle")
        )
    return "\n".join(parts)


# --------------------------------------------------------------------------- #
# Figure: the canonical frame chain
# --------------------------------------------------------------------------- #


def frame_chain_svg(
    panels: Sequence[DiagramPanel],
    relationships: Sequence[str],
    *,
    title: str = "PyTex Reference-Frame Chain",
    subtitle: str = (
        "Crystal, specimen, map, detector, laboratory and reciprocal are distinct domains. "
        "The arrows name scientific relationships, not permission to collapse the frames."
    ),
    dual_panel: DiagramPanel | None = None,
    dual_label: str = "duality",
    dual_note: str | None = None,
    dual_of: int = 0,
    panel_width: float = 176.0,
    panel_height: float = 208.0,
    gap: float = 118.0,
) -> str:
    """Render the canonical frame chain as a row of related frame panels.

    What it does
        Places one panel per frame domain, each showing that frame's own triad,
        and connects consecutive panels with a labelled relationship arrow. An
        optional dashed note beneath records the crystal-to-reciprocal duality,
        which is not part of the linear chain.

    When to use it
        As the canonical `docs/figures/reference_frames.svg` asset, and anywhere
        the whole frame vocabulary has to be seen at once.

    Parameters
    ----------
    panels:
        The frames of the **linear** chain, in order.
    relationships:
        Labels for the arrows between consecutive panels; must be one shorter
        than ``panels``.
    dual_panel:
        A frame that is *not* part of the linear chain but is dual to one of its
        members — the reciprocal frame. It is drawn below the chain and joined to
        panel ``dual_of`` by a dashed arrow, because duality is a separate
        relationship rather than another link in the chain. Putting it in the row
        would assert a chain step that does not exist.
    dual_label:
        Label for that dashed arrow.
    dual_note:
        Optional explanatory line under the dual panel.
    dual_of:
        Index of the chain panel the dual panel belongs to.

    Raises
    ------
    ValueError
        If fewer than two panels are supplied, the relationship count does not
        match, or ``dual_of`` is out of range.
    """

    if len(panels) < 2:
        raise ValueError("frame_chain_svg requires at least two panels.")
    if len(relationships) != len(panels) - 1:
        raise ValueError(
            f"frame_chain_svg needs {len(panels) - 1} relationship label(s) for "
            f"{len(panels)} panels; received {len(relationships)}."
        )

    if dual_panel is not None and not 0 <= dual_of < len(panels):
        raise ValueError(
            f"dual_of must index one of the {len(panels)} chain panels; received {dual_of}."
        )

    margin = 40.0
    header = 96.0
    width = max(
        2.0 * margin + len(panels) * panel_width + (len(panels) - 1) * gap,
        header_width(title, subtitle, margin=margin),
    )
    if dual_panel is not None:
        dual_gap = 64.0
        footer = dual_gap + panel_height + (34.0 if dual_note else 18.0)
    else:
        footer = 24.0
    height = header + panel_height + footer

    marker_id = "pytex-chain-axis"
    relation_marker = "pytex-chain-rel"
    scale = 0.20 * min(panel_width, panel_height)

    body: list[str] = []
    for index, panel in enumerate(panels):
        left = margin + index * (panel_width + gap)
        body.append(
            _panel(
                panel,
                x=left,
                y=header,
                width=panel_width,
                height=panel_height,
                marker_id=marker_id,
                scale=scale,
            )
        )
        if index < len(relationships):
            arrow_y = header + panel_height / 2.0
            body.append(
                _relationship_arrow(
                    left + panel_width + 14.0,
                    arrow_y,
                    left + panel_width + gap - 14.0,
                    arrow_y,
                    marker_id=relation_marker,
                    label=relationships[index],
                )
            )

    if dual_panel is not None:
        dual_left = margin + dual_of * (panel_width + gap)
        dual_top = header + panel_height + 64.0
        # Duality is a separate relationship, not another chain link, so the
        # dual frame hangs below its partner rather than extending the row.
        body.append(
            _relationship_arrow(
                dual_left + panel_width / 2.0,
                header + panel_height + 8.0,
                dual_left + panel_width / 2.0,
                dual_top - 8.0,
                marker_id=f"{relation_marker}-dual",
                colour=_VIOLET,
                dashed=True,
            )
        )
        body.append(
            _text(
                dual_left + panel_width / 2.0 + 12.0,
                header + panel_height + 40.0,
                dual_label,
                size=12.0,
                fill=_VIOLET,
            )
        )
        body.append(
            _panel(
                dual_panel,
                x=dual_left,
                y=dual_top,
                width=panel_width,
                height=panel_height,
                marker_id=marker_id,
                scale=scale,
            )
        )
        if dual_note:
            body.append(
                _text(
                    dual_left + panel_width + 24.0,
                    dual_top + panel_height / 2.0,
                    dual_note,
                    size=12.5,
                    fill=_VIOLET,
                )
            )

    description = (
        f"{subtitle} Frames shown in chain order: "
        + "; ".join(
            f"{p.frame.name} ({p.frame.domain.value}, axes {'/'.join(p.frame.axes)})"
            for p in panels
        )
        + ". Relationships: "
        + "; ".join(relationships)
        + "."
    )
    if dual_panel is not None:
        description += (
            f" Shown separately below the chain: {dual_panel.frame.name} "
            f"({dual_panel.frame.domain.value}, axes {'/'.join(dual_panel.frame.axes)}), "
            f"dual to {panels[dual_of].frame.name}."
        )
    markers = "\n".join(
        [
            _svg_arrow_markers(marker_id, TRIAD_AXIS_COLORS),
            _plain_marker(relation_marker, _MUTED),
            _plain_marker(f"{relation_marker}-dual", _VIOLET),
        ]
    )
    return _document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=markers,
    )


# --------------------------------------------------------------------------- #
# Figure: orientation mapping semantics
# --------------------------------------------------------------------------- #


def orientation_mapping_svg(
    crystal: ReferenceFrame,
    specimen: ReferenceFrame,
    *,
    rotation_matrix: np.ndarray,
    show_inverse: bool = False,
    title: str = "Orientation Maps The Crystal Frame Into The Specimen Frame",
    subtitle: str = (
        "An orientation is not an anonymous rotation: it is a typed mapping with a source "
        "frame and a target frame. The inverse relationship is separate and must be asked for."
    ),
) -> str:
    """Render the orientation-mapping semantics figure.

    What it does
        Shows the crystal triad, the specimen triad, and — between them — the
        crystal axes **as they land in the specimen frame**, which is what an
        orientation actually computes. The middle panel is drawn from the
        supplied rotation matrix, so the picture is the arithmetic.

    When to use it
        As the canonical `docs/figures/reference_frames_vectors.svg` asset and on
        the reference-frame concept pages.

    Parameters
    ----------
    crystal, specimen:
        The two frames of the mapping.
    rotation_matrix:
        The ``(3, 3)`` orientation matrix mapping crystal-frame components to
        specimen-frame components. Its columns are the crystal axes in specimen
        coordinates, which is exactly what the right-hand panel draws.
    show_inverse:
        Add a dashed return arrow for the inverse mapping, labelled to make the
        point that it is a **separate object** rather than an implicit default.
        Use it for the figure whose subject is the direction of the mapping;
        leave it off when the subject is the vocabulary.
    """

    matrix = np.asarray(rotation_matrix, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError("rotation_matrix must have shape (3, 3).")

    margin = 40.0
    header = 96.0
    panel_width = 300.0
    panel_height = 286.0
    # The gap grows if the heading needs more room than the panels do, so the
    # subtitle can never be clipped by the document edge.
    minimum_width = header_width(title, subtitle, margin=margin)
    gap = max(208.0, minimum_width - 2.0 * margin - 2.0 * panel_width)
    width = 2.0 * margin + 2.0 * panel_width + gap
    scale = 0.24 * min(panel_width, panel_height)

    marker_id = "pytex-map-axis"
    relation_marker = "pytex-map-rel"

    source_panel = DiagramPanel(
        frame=crystal,
        caption="the frame the orientation maps from",
    )
    # One panel, two triads: the specimen axes as the pale reference, and the
    # crystal axes where the orientation puts them. That is the whole statement.
    result_panel = DiagramPanel(
        frame=crystal,
        title=f"{crystal.name} in {specimen.name}",
        subtitle=f"components in {specimen.name} · {'/'.join(specimen.axes)}",
        basis=matrix,
        caption="crystal axes solid; specimen axes pale",
        overlay=DiagramPanel(frame=specimen),
    )

    body: list[str] = [
        _panel(
            source_panel,
            x=margin,
            y=header,
            width=panel_width,
            height=panel_height,
            marker_id=marker_id,
            scale=scale,
            title_size=17.0,
            axis_font=14.0,
        ),
        _panel(
            result_panel,
            x=margin + panel_width + gap,
            y=header,
            width=panel_width,
            height=panel_height,
            marker_id=marker_id,
            scale=scale,
            title_size=17.0,
            axis_font=14.0,
        ),
    ]

    arrow_y = header + panel_height / 2.0
    body.append(
        _relationship_arrow(
            margin + panel_width + 24.0,
            arrow_y,
            margin + panel_width + gap - 24.0,
            arrow_y,
            marker_id=relation_marker,
            label="orientation g",
            sublabel="crystal \u2192 specimen",
            colour=_ACCENT,
            label_size=14.0,
        )
    )

    note_top = header + panel_height + 26.0
    if show_inverse:
        # The inverse arc sits below the panels, visually returning from the
        # target to the source, so its separateness is the thing you see.
        arc_y = header + panel_height + 34.0
        left_x = margin + panel_width / 2.0
        right_x = margin + panel_width + gap + panel_width / 2.0
        body.append(
            f'  <path d="M{right_x:.1f},{header + panel_height - 4.0:.1f} '
            f'C{right_x:.1f},{arc_y + 34.0:.1f} {left_x:.1f},{arc_y + 34.0:.1f} '
            f'{left_x:.1f},{header + panel_height + 4.0:.1f}" fill="none" '
            f'stroke="#9f1239" stroke-width="1.8" stroke-dasharray="7 5" '
            f'marker-end="url(#{relation_marker}-inverse)"/>'
        )
        body.append(
            _text(
                0.5 * (left_x + right_x),
                arc_y + 52.0,
                "the inverse mapping is a separate object: ask for it explicitly",
                size=13.0,
                fill="#9f1239",
                anchor="middle",
            )
        )
        note_top = arc_y + 74.0

    note, _note_width, note_height = _callout(
        margin,
        note_top,
        "Reading rule",
        [
            "The right-hand solid triad is the left-hand triad re-expressed: its axes are",
            "the columns of the orientation matrix, quoted in specimen components.",
            "The pale triad is the specimen frame those components are quoted in.",
            "Reversing the mapping is a different object, never an implicit default.",
        ],
        width=width - 2.0 * margin,
    )
    body.append(note)
    height = note_top + note_height + 28.0

    description = (
        f"{subtitle} Crystal frame '{crystal.name}' with axes {'/'.join(crystal.axes)} is shown "
        f"on the left. On the right the same axes are drawn in specimen frame '{specimen.name}' "
        f"(axes {'/'.join(specimen.axes)}, pale) as the columns of the orientation matrix."
    )
    markers = "\n".join(
        [
            _svg_arrow_markers(marker_id, TRIAD_AXIS_COLORS),
            _svg_arrow_markers(f"{marker_id}-pale", PALE_TRIAD_COLORS),
            _plain_marker(relation_marker, _ACCENT),
            _plain_marker(f"{relation_marker}-inverse", "#9f1239"),
        ]
    )
    return _document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=markers,
    )


# --------------------------------------------------------------------------- #
# Figure: the hexagonal crystal frame
# --------------------------------------------------------------------------- #


def hexagonal_frame_svg(
    lattice: Lattice,
    *,
    title: str = "Canonical Hexagonal Crystal Frame",
    subtitle: str = (
        "PyTex fixes a1 along the crystal x axis, a2 at 120 degrees in the basal plane, "
        "and c along the crystal z axis. The basal axes are drawn from the lattice itself."
    ),
    width: float | None = None,
    height: float = 520.0,
) -> str:
    """Render the canonical hexagonal frame, viewed down the ``c`` axis.

    What it does
        Takes a real hexagonal `pytex.core.lattice.Lattice`, reads ``a1`` and
        ``a2`` from its direct basis, derives the redundant fourth axis
        ``a3 = -(a1 + a2)``, and draws the basal plane with the 120-degree angle
        marked, alongside the four-index notation rules.

    When to use it
        As the canonical `docs/figures/hcp_reference_frame.svg` asset and on the
        hexagonal-conventions pages.

    Why it is generated
        The basal directions come from `Lattice.direct_basis`, so the figure
        cannot disagree with the conversion helpers it illustrates.

    Raises
    ------
    ValueError
        If the lattice is not hexagonal in setting (``gamma`` must be 120
        degrees with ``alpha = beta = 90``).
    """

    if not (
        np.isclose(lattice.gamma_deg, 120.0, atol=1e-6)
        and np.isclose(lattice.alpha_deg, 90.0, atol=1e-6)
        and np.isclose(lattice.beta_deg, 90.0, atol=1e-6)
    ):
        raise ValueError(
            "hexagonal_frame_svg requires a hexagonal setting with alpha = beta = 90 deg "
            f"and gamma = 120 deg; received alpha={lattice.alpha_deg}, "
            f"beta={lattice.beta_deg}, gamma={lattice.gamma_deg}."
        )

    if width is None:
        width = max(940.0, header_width(title, subtitle))

    basis = np.asarray(lattice.direct_basis().matrix, dtype=np.float64)
    a1 = basis[:, 0] / np.linalg.norm(basis[:, 0])
    a2 = basis[:, 1] / np.linalg.norm(basis[:, 1])
    a3 = -(a1 + a2)
    a3 = a3 / np.linalg.norm(a3)
    c_axis = basis[:, 2] / np.linalg.norm(basis[:, 2])

    # Viewed straight down c, so the basal plane lies in the page and the 120
    # degree relationships are true angles rather than foreshortened ones.
    axes = np.vstack([a1, a2, a3, c_axis])
    screen, _ = project_orthographic(axes, elev_deg=90.0, azim_deg=-90.0)

    centre = (272.0, 300.0)
    radius = 128.0
    marker_id = "pytex-hex"
    colours = (TRIAD_AXIS_COLORS[0], TRIAD_AXIS_COLORS[1], "#b45309", TRIAD_AXIS_COLORS[2])
    labels = ("a1", "a2", "a3 = -(a1 + a2)", "c")

    body: list[str] = [
        _card(40.0, 96.0, 464.0, height - 136.0, fill=_PANEL, stroke=_PANEL_STROKE),
        f'  <circle cx="{centre[0]:.1f}" cy="{centre[1]:.1f}" r="{radius:.1f}" '
        f'fill="none" stroke="#c9d6ea" stroke-width="1.4"/>',
    ]

    for index in range(3):
        end_x = centre[0] + radius * float(screen[index, 0])
        end_y = centre[1] - radius * float(screen[index, 1])
        body.append(
            f'  <line x1="{centre[0]:.1f}" y1="{centre[1]:.1f}" x2="{end_x:.1f}" '
            f'y2="{end_y:.1f}" stroke="{colours[index]}" stroke-width="3.2" '
            f'stroke-linecap="round" marker-end="url(#{marker_id}-{index})"/>'
        )
        label_x = centre[0] + (radius + 34.0) * float(screen[index, 0])
        label_y = centre[1] - (radius + 34.0) * float(screen[index, 1])
        body.append(
            _text(label_x, label_y + 5.0, labels[index], size=15.0,
                  fill=colours[index], anchor="middle")
        )

    # The c axis points at the viewer in this projection, so it is drawn as a
    # ringed dot rather than a degenerate zero-length arrow.
    body.append(
        f'  <circle cx="{centre[0]:.1f}" cy="{centre[1]:.1f}" r="7" fill="none" '
        f'stroke="{colours[3]}" stroke-width="2.4"/>'
    )
    body.append(
        f'  <circle cx="{centre[0]:.1f}" cy="{centre[1]:.1f}" r="2.6" fill="{colours[3]}"/>'
    )
    # a1 runs right and a3 runs down-left, and the 120 degree arc sweeps the
    # upper region, so the lower-right is the only clear place for this label.
    body.append(
        _text(centre[0] + 16.0, centre[1] + 26.0, "c (out of page)", size=13.0, fill=colours[3])
    )

    # 120 degree arc between a1 and a2.
    arc_r = 58.0
    start = (centre[0] + arc_r * float(screen[0, 0]), centre[1] - arc_r * float(screen[0, 1]))
    end = (centre[0] + arc_r * float(screen[1, 0]), centre[1] - arc_r * float(screen[1, 1]))
    body.append(
        f'  <path d="M{start[0]:.1f},{start[1]:.1f} A{arc_r:.1f},{arc_r:.1f} 0 0 0 '
        f'{end[0]:.1f},{end[1]:.1f}" fill="none" stroke="#b45309" stroke-width="1.8"/>'
    )
    body.append(_text(centre[0] + 4.0, centre[1] - arc_r - 12.0, "120°", size=13.0,
                      fill="#b45309", anchor="middle"))

    note_top = 132.0
    note, _note_width, note_height = _callout(
        552.0,
        note_top,
        "Notation fixed in docs and code",
        [
            "3-index direct-space direction: [u v w]",
            "4-index Weber / Miller-Bravais direction: [U V T W]",
            "Constraint: U + V + T = 0",
            "Plane notation: (h k i l), with i = -(h + k)",
            "Reciprocal basis axes carry the star: a*, b*, c*",
            "Internal storage stays canonical 3-component",
        ],
        width=348.0,
    )
    body.append(note)
    provenance_top = note_top + note_height + 28.0
    body.append(
        _text(
            552.0,
            provenance_top,
            f"Drawn from a = {lattice.a:.4f} angstrom, c = {lattice.c:.4f} angstrom",
            size=12.0,
            fill=_MUTED,
        )
    )
    body.append(
        _text(
            552.0,
            provenance_top + 20.0,
            f"c/a = {lattice.c / lattice.a:.4f}, gamma = {lattice.gamma_deg:.1f} deg",
            size=12.0,
            fill=_MUTED,
        )
    )
    body.append(
        _text(
            552.0,
            provenance_top + 40.0,
            "Basal directions read from Lattice.direct_basis(), not hand-placed.",
            size=12.0,
            fill=_MUTED,
        )
    )

    description = (
        f"{subtitle} Basal axes a1, a2 and the derived a3 = -(a1 + a2) are read from the "
        f"lattice direct basis with a = {lattice.a:.4f} angstrom and c = {lattice.c:.4f} "
        "angstrom; the c axis points out of the page in this basal-plane view."
    )
    markers = "\n".join(_plain_marker(f"{marker_id}-{i}", colours[i]) for i in range(3))
    return _document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=markers,
    )


# --------------------------------------------------------------------------- #
# Figure: active versus passive rotation language
# --------------------------------------------------------------------------- #


def active_passive_svg(
    crystal: ReferenceFrame,
    specimen: ReferenceFrame,
    *,
    rotation_matrix: np.ndarray,
    vector: Sequence[float] = (1.0, 1.0, 0.0),
    title: str = "Active Versus Passive Rotation Views",
    subtitle: str = (
        "The same mathematics, two languages. PyTex keeps the frame-mapping language primary "
        "because it is the safer one at scientific workflow boundaries."
    ),
) -> str:
    """Render the active-versus-passive rotation figure.

    What it does
        Left panel, *active*: the frame is held fixed and a vector moves, drawn
        as the original vector and its rotated image. Right panel, *passive*: the
        vector is held fixed and the basis changes, drawn as the crystal triad
        against the pale specimen triad. Both panels use the same rotation, so
        the equivalence is visible rather than asserted.

    When to use it
        As the canonical `docs/figures/active_passive_rotation.svg` asset, on the
        pages that fix orientation language.

    Parameters
    ----------
    crystal, specimen:
        The two frames of the passive view.
    rotation_matrix:
        The ``(3, 3)`` rotation shared by both views.
    vector:
        The vector the active view rotates, in frame components. Choose one whose
        rotated image separates clearly in the projected view, or the two arrows
        superimpose and the panel stops teaching anything; the default is checked
        against the default rotation for exactly that.
    """

    matrix = np.asarray(rotation_matrix, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError("rotation_matrix must have shape (3, 3).")
    source = np.asarray(vector, dtype=np.float64).reshape(3)
    source = source / np.linalg.norm(source)
    rotated = matrix @ source

    margin = 40.0
    header = 96.0
    panel_width = 384.0
    panel_height = 300.0
    gap = 56.0
    width = max(
        2.0 * margin + 2.0 * panel_width + gap,
        header_width(title, subtitle, margin=margin),
    )
    scale = 0.20 * min(panel_width, panel_height)
    marker_id = "pytex-ap-axis"
    vector_marker = "pytex-ap-vec"

    body: list[str] = []

    # --- active view: fixed basis, moving vector -------------------------- #
    left = margin
    body.append(
        _panel(
            DiagramPanel(
                frame=specimen,
                title="Active view",
                subtitle="the basis is held fixed and the vector moves",
                caption="v and its rotated image R v, in one unchanged frame",
            ),
            x=left,
            y=header,
            width=panel_width,
            height=panel_height,
            marker_id=f"{marker_id}-pale",
            scale=scale,
            colors=PALE_TRIAD_COLORS,
            title_size=17.0,
            axis_font=13.0,
        )
    )
    centre = (left + panel_width / 2.0, header + panel_height / 2.0 + 30.0)
    screen, _ = project_orthographic(np.vstack([source, rotated]))
    for index, (label, colour) in enumerate((("v", "#b45309"), ("R v", "#9f1239"))):
        end = (
            centre[0] + 1.25 * scale * float(screen[index, 0]),
            centre[1] - 1.25 * scale * float(screen[index, 1]),
        )
        body.append(
            f'  <line x1="{centre[0]:.1f}" y1="{centre[1]:.1f}" x2="{end[0]:.1f}" '
            f'y2="{end[1]:.1f}" stroke="{colour}" stroke-width="3.4" stroke-linecap="round" '
            f'marker-end="url(#{vector_marker}-{index})"/>'
        )
        body.append(
            _text(
                centre[0] + 1.52 * scale * float(screen[index, 0]),
                centre[1] - 1.52 * scale * float(screen[index, 1]) + 5.0,
                label,
                size=15.0,
                fill=colour,
                anchor="middle",
            )
        )

    # --- passive view: fixed vector, changing basis ------------------------ #
    right = margin + panel_width + gap
    body.append(
        _panel(
            DiagramPanel(
                frame=crystal,
                title="Passive / frame-mapping view",
                subtitle="the object is fixed and the basis describing it changes",
                basis=matrix,
                caption="crystal axes solid; specimen axes pale",
                overlay=DiagramPanel(frame=specimen),
            ),
            x=right,
            y=header,
            width=panel_width,
            height=panel_height,
            marker_id=marker_id,
            scale=scale,
            title_size=17.0,
            axis_font=13.0,
        )
    )

    note, _note_width, note_height = _callout(
        margin,
        header + panel_height + 26.0,
        "Which language PyTex uses",
        [
            "PyTex documents an Orientation in the passive, frame-mapping language, because it",
            "keeps the source and target frames expressed rather than implied.",
            "The active view stays mathematically valid; it is simply easier to misread once",
            "crystal and specimen frames are both in play.",
        ],
        width=width - 2.0 * margin,
    )
    body.append(note)
    height = header + panel_height + 26.0 + note_height + 28.0

    description = (
        f"{subtitle} Left: the specimen frame '{specimen.name}' held fixed while a vector is "
        "rotated. Right: the same rotation read as a change of basis, with crystal frame "
        f"'{crystal.name}' drawn against the pale specimen frame."
    )
    markers = "\n".join(
        [
            _svg_arrow_markers(marker_id, TRIAD_AXIS_COLORS),
            _svg_arrow_markers(f"{marker_id}-pale", PALE_TRIAD_COLORS),
            _plain_marker(f"{vector_marker}-0", "#b45309"),
            _plain_marker(f"{vector_marker}-1", "#9f1239"),
        ]
    )
    return _document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=markers,
    )


# --------------------------------------------------------------------------- #
# Figure: the Bunge Euler sequence
# --------------------------------------------------------------------------- #


def euler_sequence_svg(
    frame: ReferenceFrame,
    *,
    phi1_deg: float,
    Phi_deg: float,  # noqa: N803 - Bunge's canonical symbol is capital Phi
    phi2_deg: float,
    title: str = "Bunge Euler Sequence",
    subtitle: str = (
        "Three ordered rotations, ZXZ: phi1 about z, then Phi about the rotated x, then phi2 "
        "about the final z. Each panel is the frame after that step, computed not sketched."
    ),
) -> str:
    """Render the Bunge ZXZ Euler sequence as a row of computed frame states.

    What it does
        Builds the three partial rotations with
        `pytex.core.orientation.Rotation.from_bunge_euler` and draws the frame
        after each one, so the panels *are* the arithmetic. The angle labels on
        the arrows name the rotation and the axis it turns about.

    When to use it
        As the canonical `docs/figures/bunge_euler_geometry.svg` asset.

    Parameters
    ----------
    frame:
        The frame being rotated, used for the axis labels.
    phi1_deg, Phi_deg, phi2_deg:
        The Bunge angles, in degrees.
    """

    from pytex.core.orientation import Rotation

    steps = (
        ("initial frame", Rotation.identity()),
        ("after phi1", Rotation.from_bunge_euler(phi1_deg, 0.0, 0.0)),
        ("after Phi", Rotation.from_bunge_euler(phi1_deg, Phi_deg, 0.0)),
        ("after phi2", Rotation.from_bunge_euler(phi1_deg, Phi_deg, phi2_deg)),
    )
    relationships = (
        f"phi1 = {phi1_deg:g}° about z",
        f"Phi = {Phi_deg:g}° about x'",
        f"phi2 = {phi2_deg:g}° about z''",
    )

    margin = 40.0
    header = 96.0
    panel_width = 196.0
    panel_height = 216.0
    gap = 132.0
    width = max(
        2.0 * margin + 4.0 * panel_width + 3.0 * gap,
        header_width(title, subtitle, margin=margin),
    )
    scale = 0.21 * min(panel_width, panel_height)
    marker_id = "pytex-euler-axis"
    relation_marker = "pytex-euler-rel"

    body: list[str] = []
    for index, (caption, rotation) in enumerate(steps):
        left = margin + index * (panel_width + gap)
        body.append(
            _panel(
                DiagramPanel(
                    frame=frame,
                    title=caption,
                    subtitle="" if index == 0 else "axes after this step",
                    basis=np.asarray(rotation.as_matrix(), dtype=np.float64),
                ),
                x=left,
                y=header,
                width=panel_width,
                height=panel_height,
                marker_id=marker_id,
                scale=scale,
            )
        )
        if index < len(relationships):
            arrow_y = header + panel_height / 2.0
            body.append(
                _relationship_arrow(
                    left + panel_width + 12.0,
                    arrow_y,
                    left + panel_width + gap - 12.0,
                    arrow_y,
                    marker_id=relation_marker,
                    label=relationships[index],
                    colour=_ACCENT,
                    label_size=12.0,
                )
            )

    note, _note_width, note_height = _callout(
        margin,
        header + panel_height + 26.0,
        "Reading rule",
        [
            "The angle names are not positional placeholders: they belong to this ordered ZXZ",
            "construction, which is why PyTex keeps a named Bunge entry point rather than",
            "treating every three-angle input as interchangeable.",
            "Composed form: R = Rz(phi2) Rx(Phi) Rz(phi1).",
        ],
        width=width - 2.0 * margin,
    )
    body.append(note)
    height = header + panel_height + 26.0 + note_height + 28.0

    description = (
        f"{subtitle} Panels show frame '{frame.name}' (axes {'/'.join(frame.axes)}) initially and "
        f"after each of phi1 = {phi1_deg:g}, Phi = {Phi_deg:g}, phi2 = {phi2_deg:g} degrees."
    )
    markers = "\n".join(
        [_svg_arrow_markers(marker_id, TRIAD_AXIS_COLORS), _plain_marker(relation_marker, _ACCENT)]
    )
    return _document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=markers,
    )
