"""Shared building blocks for generated documentation SVG.

Every canonical figure in `docs/figures/` that is produced by code — the
reference-frame diagrams, the algorithm flow sheets — is assembled from the
elements here: the palette tokens from
`docs/standards/visualization_style_guide.md`, an Arial advance-width estimate
so labels can be placed without colliding, rounded cards, labelled arrows,
callout boxes, and a document wrapper that emits the ``<title>`` and ``<desc>``
the style guide requires.

Two rules are enforced by construction rather than left to the author.

**Arrowheads never scale with stroke width.** SVG markers default to
``markerUnits="strokeWidth"``, multiplying the head by the line's stroke width,
so a 12-unit head on a 4-unit stroke renders at 48 units. `arrow_marker` always
sets ``userSpaceOnUse``. Across the previously hand-authored figure set this one
default produced arrowheads occupying up to 125% of the lines they annotated.

**Text is measured before it is placed.** SVG carries no offline text metrics,
so `text_width` estimates the rendered width from Arial's average advance and
callers size panels to their content. It is an estimate: leave a margin rather
than butting elements against the result.

Everything here is pure Python — no matplotlib — so documentation figures can be
produced in a minimal environment.

See also
--------
`pytex.plotting.frame_diagrams` : the reference-frame figures.
`pytex.plotting.algorithm_diagrams` : the algorithm flow sheets.
`docs/standards/visualization_style_guide.md` : the canonical tokens.
"""

from __future__ import annotations

from collections.abc import Sequence
from xml.sax.saxutils import escape

#: Canonical documentation tokens (docs/standards/visualization_style_guide.md).
INK = "#07122f"
MUTED = "#40506f"
PAPER = "#fbfdff"
PANEL = "#ffffff"
PANEL_STROKE = "#d7e0ef"
ACCENT = "#2563eb"
TEAL = "#0f9f9f"
VIOLET = "#7c3aed"
AMBER = "#f59e0b"
ROSE = "#e11d48"
GREEN = "#16a34a"

#: The one font used by every canonical SVG asset.
#:
#: Arial alone, with no fallback stack. The stack was dropped deliberately: a
#: figure that renders in Arial on one machine and Helvetica or a generic
#: sans-serif on another is not a canonical asset, because the text metrics in
#: :mod:`pytex.plotting._svg_text` — which is what keeps labels inside their
#: cards without a renderer — are Arial's. Naming the fallbacks invited exactly
#: the substitution the metrics cannot account for.
FONT = "Arial"

#: Fill/stroke pairs derived from the tokens, for cards that carry a role.
#:
#: Colour groups concepts; it never carries information a label does not also
#: carry, so a reader who cannot distinguish the hues loses nothing.
ROLE_TINTS: dict[str, tuple[str, str]] = {
    "input": ("#eef7f0", "#a7ccb4"),
    "compute": ("#eef4ff", "#b9cdf2"),
    "decision": ("#fdf6e8", "#e6d4a8"),
    "reject": ("#fdeef1", "#eebfc8"),
    "output": ("#f3eefc", "#cdbdee"),
    "note": ("#f7faff", "#c8d8f2"),
}

#: Text colour to pair with each role tint.
ROLE_INK: dict[str, str] = {
    "input": "#124a2a",
    "compute": "#0f2f6b",
    "decision": "#6b4a06",
    "reject": "#8a1030",
    "output": "#42227e",
    "note": INK,
}

# Arial advance widths are ~0.52 em on average for mixed-case text; the digits
# and capitals that dominate these labels run slightly wider.
_AVERAGE_ADVANCE_EM = 0.55

__all__ = [
    "ACCENT",
    "AMBER",
    "FONT",
    "GREEN",
    "INK",
    "MUTED",
    "PANEL",
    "PANEL_STROKE",
    "PAPER",
    "ROLE_INK",
    "ROLE_TINTS",
    "ROSE",
    "TEAL",
    "VIOLET",
    "arrow_marker",
    "callout",
    "card",
    "document",
    "header_width",
    "relationship_arrow",
    "text",
    "text_width",
    "wrap_text",
]


def text_width(content: str, font_size: float) -> float:
    """Estimate the rendered width of an Arial string, in SVG user units.

    Used to size panels and place labels so nothing collides. It is an estimate,
    not a measurement — SVG has no text metrics available offline — so callers
    should leave a margin rather than butt elements against the result.
    """

    return _AVERAGE_ADVANCE_EM * float(font_size) * len(str(content))


def wrap_text(content: str, font_size: float, max_width: float) -> list[str]:
    """Greedy word wrap to an estimated pixel width.

    Returns at least one line even when a single word exceeds ``max_width``,
    because dropping content is worse than overflowing it — and the figure
    layout guards will catch the overflow.
    """

    words = str(content).split()
    if not words:
        return [""]
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = f"{current} {word}"
        if text_width(candidate, font_size) <= max_width:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def text(
    x: float,
    y: float,
    content: str,
    *,
    size: float = 13.0,
    fill: str = INK,
    anchor: str = "start",
    weight: str | None = None,
) -> str:
    """One SVG text element with the canonical font."""

    bold = f' font-weight="{weight}"' if weight else ""
    return (
        f'  <text x="{x:.1f}" y="{y:.1f}" font-size="{size:.0f}" font-family="{FONT}" '
        f'fill="{fill}" text-anchor="{anchor}"{bold}>{escape(str(content))}</text>'
    )


def card(
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    fill: str = PANEL,
    stroke: str = PANEL_STROKE,
    radius: float = 8.0,
    stroke_width: float = 1.0,
    dashed: bool = False,
) -> str:
    """A rounded panel card; 8px radius per the style guide."""

    dash = ' stroke-dasharray="6 5"' if dashed else ""
    return (
        f'  <rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" '
        f'rx="{radius:.0f}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"{dash}/>'
    )


def relationship_arrow(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    marker_id: str,
    label: str | None = None,
    sublabel: str | None = None,
    colour: str = MUTED,
    dashed: bool = False,
    label_size: float = 12.0,
) -> str:
    """A labelled relationship arrow between two elements.

    The label sits above the line and the optional sublabel below it, both
    centred on the midpoint, so a relationship reads without crossing the arrow
    it describes.
    """

    dash = ' stroke-dasharray="7 5"' if dashed else ""
    parts = [
        f'  <line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{colour}" stroke-width="2"{dash} marker-end="url(#{marker_id})"/>'
    ]
    mid_x = 0.5 * (x1 + x2)
    mid_y = 0.5 * (y1 + y2)
    if label is not None:
        parts.append(text(mid_x, mid_y - 9.0, label, size=label_size, fill=INK, anchor="middle"))
    if sublabel is not None:
        offset = 20.0 if label is not None else 14.0
        parts.append(
            text(
                mid_x,
                mid_y + offset - 9.0,
                sublabel,
                size=label_size - 1.0,
                fill=MUTED,
                anchor="middle",
            )
        )
    return "\n".join(parts)


def callout(
    x: float,
    y: float,
    title: str,
    lines: Sequence[str],
    *,
    width: float | None = None,
    title_size: float = 14.0,
    line_size: float = 12.0,
    line_height: float = 19.0,
    fill: str = "#f7faff",
    stroke: str = "#c8d8f2",
) -> tuple[str, float, float]:
    """A titled note box. Returns the fragment and its rendered size."""

    if width is None:
        widest = max([text_width(title, title_size)] + [text_width(t, line_size) for t in lines])
        width = widest + 36.0
    height = 34.0 + line_height * len(lines) + 12.0
    parts = [card(x, y, width, height, fill=fill, stroke=stroke)]
    parts.append(text(x + 16.0, y + 25.0, title, size=title_size, fill=INK))
    for index, line in enumerate(lines):
        parts.append(
            text(x + 16.0, y + 48.0 + line_height * index, line, size=line_size, fill=MUTED)
        )
    return "\n".join(parts), width, height


def header_width(title: str, subtitle: str, *, margin: float = 40.0) -> float:
    """Minimum document width that fits the heading and subtitle without clipping."""

    return margin * 2.0 + max(text_width(title, 25.0), text_width(subtitle, 13.0))


def arrow_marker(marker_id: str, colour: str) -> str:
    """An arrowhead that does **not** scale with stroke width.

    ``markerUnits="userSpaceOnUse"`` is the whole point: without it SVG
    multiplies the head by the line's stroke width, which is the defect that
    motivated generating these figures in the first place.
    """

    return (
        f'    <marker id="{marker_id}" markerUnits="userSpaceOnUse" markerWidth="12" '
        f'markerHeight="9" refX="11" refY="4.5" orient="auto">\n'
        f'      <path d="M0,0 L12,4.5 L0,9 z" fill="{colour}"/>\n'
        f"    </marker>"
    )


def document(
    *,
    width: float,
    height: float,
    title: str,
    description: str,
    subtitle: str,
    body: str,
    marker_defs: str,
) -> str:
    """Wrap fragments into a complete, accessible SVG document.

    The ``<title>`` and ``<desc>`` are mandatory: they are what a screen reader
    announces and what the figure guard tests check for.
    """

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img">\n'
        f"  <title>{escape(title)}</title>\n"
        f"  <desc>{escape(description)}</desc>\n"
        f"  <defs>\n{marker_defs}\n  </defs>\n"
        f'  <rect width="{width:.0f}" height="{height:.0f}" fill="{PAPER}"/>\n'
        f"{text(40.0, 46.0, title, size=25.0, fill=INK)}\n"
        f"{text(40.0, 72.0, subtitle, size=13.0, fill=MUTED)}\n"
        f"{body}\n"
        f"</svg>\n"
    )
