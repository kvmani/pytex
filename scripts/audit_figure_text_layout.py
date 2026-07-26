"""Find text that overflows its box, overflows the canvas, or collides, in canonical SVGs.

Why this exists
---------------

A figure can be perfectly well-formed SVG and still be unreadable: a caption
wider than the panel it sits in spills over the panel border, two labels sharing
a baseline overprint each other, and a title wider than the canvas is simply cut
off. None of that is visible in the source, and none of it is caught by an XML
parser.

What it measures
----------------

For every ``<text>`` element the auditor resolves the ``translate(...)``
transforms of its ancestor ``<g>`` elements — much of the text in hand-authored
figures lives inside such groups, and an auditor that ignores them checks almost
nothing — then estimates the rendered width from the font size and string length
and reports three defects:

``canvas``
    the text extends past the document's own ``viewBox``.
``box``
    the text extends past the ``<rect>`` that encloses its anchor point. This is
    the "text running outside its box" case.
``collision``
    two texts on the same or nearly the same baseline overlap horizontally.
    "Nearly" matters: two runs half a line apart still overprint each other, and
    an exact-baseline test misses them.
``over-card``
    a text run is painted across a filled panel it does not belong to — a
    centred header sitting over the top card of a radial layout, say. This is
    invisible to a containment test, because the text is not *in* that box.

Widths come from `pytex.plotting._svg_text`, the Helvetica advance-width table,
so generated and hand-authored figures are judged by one ruler and capital-heavy
strings are not under-measured. Font sizes are resolved from the element's
``font-size`` attribute, from the ``class`` it carries against the document's
``<style>`` block, or from an ancestor — these figures declare most of their
typography in CSS classes, and an auditor that ignores that measures almost
nothing correctly.

A tolerance is applied: the auditor finds defects for a human to confirm, it is
not a pixel-accurate renderer.

Usage::

    python scripts/audit_figure_text_layout.py             # report every defect
    python scripts/audit_figure_text_layout.py --check      # exit 1 if any remain
    python scripts/audit_figure_text_layout.py --figure X   # one figure, verbose
"""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from pytex.plotting._svg_text import text_width  # noqa: E402

FIGURES_DIR = REPO_ROOT / "docs" / "figures"
SVG_NS = "{http://www.w3.org/2000/svg}"

#: Fallback font size when neither an attribute, a class, nor an ancestor sets one.
DEFAULT_FONT_SIZE = 12.0

#: Slack in user units before a near-miss is called a defect. Generous, because
#: the width model is an estimate and a one-pixel graze is not a defect.
TOLERANCE = 6.0

_TRANSLATE_RE = re.compile(r"translate\(\s*(-?[0-9.]+)[ ,]+(-?[0-9.]+)\s*\)")
# `.name{font:700 20px Arial...}` or `.name{font-size:15px}`
_CLASS_RULE_RE = re.compile(r"\.([A-Za-z0-9_-]+)\s*\{([^}]*)\}")
_FONT_SHORTHAND_RE = re.compile(r"font\s*:\s*([^;}]*)")
_FONT_SIZE_RE = re.compile(r"font-size\s*:\s*([0-9.]+)px")
_SIZE_IN_SHORTHAND_RE = re.compile(r"(?:^|\s)([0-9.]+)px")


def _stylesheet(root: ET.Element) -> dict[str, tuple[float, bool]]:
    """Font size and boldness per CSS class, from the document's <style> block."""

    styles: dict[str, tuple[float, bool]] = {}
    for style in root.iter(f"{SVG_NS}style"):
        css = "".join(style.itertext())
        for name, body in _CLASS_RULE_RE.findall(css):
            size: float | None = None
            bold = False
            shorthand = _FONT_SHORTHAND_RE.search(body)
            if shorthand:
                declaration = shorthand.group(1)
                bold = bool(re.search(r"\b(bold|[5-9]00)\b", declaration))
                size_match = _SIZE_IN_SHORTHAND_RE.search(declaration)
                if size_match:
                    size = float(size_match.group(1))
            explicit = _FONT_SIZE_RE.search(body)
            if explicit:
                size = float(explicit.group(1))
            if re.search(r"font-weight\s*:\s*(bold|[5-9]00)", body):
                bold = True
            if size is not None:
                styles[name] = (size, bold)
    return styles


@dataclass(frozen=True)
class TextSpan:
    """One rendered text run, resolved into absolute figure coordinates."""

    left: float
    right: float
    baseline: float
    size: float
    content: str

    @property
    def width(self) -> float:
        return self.right - self.left


@dataclass(frozen=True)
class Defect:
    kind: str
    figure: str
    detail: str


def _translate(element: ET.Element) -> tuple[float, float]:
    transform = element.get("transform", "")
    match = _TRANSLATE_RE.search(transform)
    if match is None:
        return 0.0, 0.0
    return float(match.group(1)), float(match.group(2))


def _resolve_font(
    element: ET.Element,
    styles: dict[str, tuple[float, bool]],
    inherited: tuple[float, bool],
) -> tuple[float, bool]:
    """Font size and boldness for one element: attribute, then class, then parent."""

    size, bold = inherited
    class_name = element.get("class")
    if class_name:
        for token in class_name.split():
            if token in styles:
                size, bold = styles[token]
    attribute = element.get("font-size")
    if attribute:
        try:
            size = float(re.sub(r"[^0-9.]", "", attribute) or size)
        except ValueError:  # pragma: no cover - malformed attribute
            pass
    weight = element.get("font-weight", "")
    if weight and re.match(r"bold|[5-9]00", weight):
        bold = True
    return size, bold


def _collect(
    element: ET.Element,
    offset: tuple[float, float],
    texts: list[TextSpan],
    boxes: list[tuple[float, float, float, float]],
    styles: dict[str, tuple[float, bool]],
    font: tuple[float, bool],
) -> None:
    """Walk the tree accumulating translate offsets, fonts, texts and rects."""

    dx, dy = offset
    for child in element:
        cdx, cdy = _translate(child)
        child_offset = (dx + cdx, dy + cdy)
        child_font = _resolve_font(child, styles, font)
        tag = child.tag
        if tag == f"{SVG_NS}text":
            size, bold = child_font
            anchor = child.get("text-anchor", "start")
            base_x = float(child.get("x", "0") or 0) + child_offset[0]
            base_y = float(child.get("y", "0") or 0) + child_offset[1]
            spans = list(child.iter(f"{SVG_NS}tspan"))
            if any(span.get("x") is not None for span in spans):
                # Wrapped text: a tspan carrying its own `x` starts a new line.
                # A tspan without one (a subscript inside a formula, say) is an
                # inline continuation and belongs to the line it sits in —
                # treating those as separate lines reports phantom collisions
                # between the parts of a single equation.
                cursor = base_y
                current: list[str] = []
                current_x = base_x

                def flush(
                    line_parts: list[str],
                    line_x: float,
                    line_y: float,
                    _size: float = size,
                    _bold: bool = bold,
                    _anchor: str = anchor,
                ) -> None:
                    content = "".join(line_parts).strip()
                    if not content:
                        return
                    estimated = text_width(content, _size, bold=_bold)
                    left = line_x
                    if _anchor == "middle":
                        left -= estimated / 2.0
                    elif _anchor == "end":
                        left -= estimated
                    texts.append(TextSpan(left, left + estimated, line_y, _size, content))

                if (child.text or "").strip():
                    current.append(child.text or "")
                for span in spans:
                    own_x = span.get("x")
                    if own_x is not None:
                        flush(current, current_x, cursor)
                        current = []
                        current_x = float(own_x) + child_offset[0]
                        cursor += float(span.get("dy", "0") or 0)
                    current.append(span.text or "")
                    if span.tail:
                        current.append(span.tail)
                flush(current, current_x, cursor)
            elif spans:
                # Inline tspans only: measure the whole run as one line.
                content = "".join(child.itertext()).strip()
                if content:
                    estimated = text_width(content, size, bold=bold)
                    left = base_x
                    if anchor == "middle":
                        left -= estimated / 2.0
                    elif anchor == "end":
                        left -= estimated
                    texts.append(TextSpan(left, left + estimated, base_y, size, content))
            else:
                content = (child.text or "").strip()
                if content:
                    estimated = text_width(content, size, bold=bold)
                    left = base_x
                    if anchor == "middle":
                        left -= estimated / 2.0
                    elif anchor == "end":
                        left -= estimated
                    texts.append(TextSpan(left, left + estimated, base_y, size, content))
        elif tag == f"{SVG_NS}rect":
            try:
                boxes.append(
                    (
                        float(child.get("x", "0") or 0) + child_offset[0],
                        float(child.get("y", "0") or 0) + child_offset[1],
                        float(child.get("width", "0") or 0),
                        float(child.get("height", "0") or 0),
                    )
                )
            except ValueError:  # pragma: no cover - malformed attribute
                pass
        _collect(child, child_offset, texts, boxes, styles, child_font)


def _enclosing_box(
    span: TextSpan, boxes: list[tuple[float, float, float, float]], canvas: tuple[float, float]
) -> tuple[float, float, float, float] | None:
    """The smallest rect that plausibly *contains* the text, if any.

    Containment needs more than the anchor point falling inside the rect. A
    centred figure title whose baseline happens to land a few units below the
    top edge of some card is not inside that card — its ascender is above the
    card entirely. Requiring the ascender to clear the top edge rejects those
    coincidences while still catching every real caption-in-a-card.

    The full-canvas background rect is skipped: overflowing it is the ``canvas``
    defect, reported separately.
    """

    anchor_x = 0.5 * (span.left + span.right)
    anchor_y = span.baseline
    ascender = 0.75 * span.size
    candidates = [
        box
        for box in boxes
        if box[2] * box[3] > 0
        and not (box[2] >= canvas[0] - 1.0 and box[3] >= canvas[1] - 1.0)
        and box[0] - 1.0 <= anchor_x <= box[0] + box[2] + 1.0
        and box[1] + ascender <= anchor_y <= box[1] + box[3] + 1.0
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda b: b[2] * b[3])


def audit(path: Path) -> list[Defect]:
    """Every layout defect found in one figure."""

    root = ET.fromstring(path.read_text(encoding="utf-8"))
    viewbox = root.get("viewBox")
    if viewbox:
        parts = viewbox.split()
        canvas = (float(parts[2]), float(parts[3]))
    else:
        canvas = (float(root.get("width", "0") or 0), float(root.get("height", "0") or 0))

    texts: list[TextSpan] = []
    boxes: list[tuple[float, float, float, float]] = []
    styles = _stylesheet(root)
    _collect(root, (0.0, 0.0), texts, boxes, styles, (DEFAULT_FONT_SIZE, False))

    defects: list[Defect] = []
    name = path.name

    for span in texts:
        if span.right > canvas[0] + TOLERANCE:
            defects.append(
                Defect(
                    "canvas",
                    name,
                    f"{span.content[:46]!r} extends {span.right - canvas[0]:.0f}u past the canvas",
                )
            )
            continue
        box = _enclosing_box(span, boxes, canvas)
        if box is None:
            continue
        box_right = box[0] + box[2]
        if span.right > box_right + TOLERANCE:
            defects.append(
                Defect(
                    "box",
                    name,
                    f"{span.content[:46]!r} extends {span.right - box_right:.0f}u past its box",
                )
            )

    ordered = sorted(texts, key=lambda s: (s.baseline, s.left))
    for index, first in enumerate(ordered):
        for second in ordered[index + 1 :]:
            # Runs closer than about half a line height overprint each other even
            # when their baselines are not identical.
            vertical_gap = second.baseline - first.baseline
            if vertical_gap > 0.55 * max(first.size, second.size):
                break
            if second.left < first.right - TOLERANCE and first.left < second.right - TOLERANCE:
                defects.append(
                    Defect(
                        "collision",
                        name,
                        f"{first.content[:28]!r} overlaps {second.content[:28]!r}",
                    )
                )

    # Text painted across a panel it is not inside. A run may legitimately sit
    # inside several nested panels (a card within a section), so every box that
    # contains it is excluded, not merely the tightest one.
    for span in texts:
        anchor_x = 0.5 * (span.left + span.right)
        containers = {
            id(box)
            for box in boxes
            if box[0] - 1.0 <= anchor_x <= box[0] + box[2] + 1.0
            and box[1] - 1.0 <= span.baseline <= box[1] + box[3] + 1.0
        }
        for box in boxes:
            if id(box) in containers or box[2] * box[3] <= 0:
                continue
            if box[2] >= canvas[0] - 1.0 and box[3] >= canvas[1] - 1.0:
                continue
            if box[2] * box[3] < 2000.0:
                continue  # decorative rules and accent bars, not panels
            horizontal = span.left < box[0] + box[2] - TOLERANCE and span.right > box[0] + TOLERANCE
            top, bottom = span.baseline - 0.75 * span.size, span.baseline + 0.25 * span.size
            vertical = top < box[1] + box[3] - 1.0 and bottom > box[1] + 1.0
            if horizontal and vertical:
                defects.append(
                    Defect(
                        "over-card",
                        name,
                        f"{span.content[:34]!r} is painted across a panel it is not inside",
                    )
                )
                break
    return defects


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="exit 1 if any defect remains")
    parser.add_argument("--figure", help="audit a single figure by name")
    args = parser.parse_args()

    paths = sorted(FIGURES_DIR.glob("*.svg"))
    if args.figure:
        paths = [p for p in paths if args.figure in p.name]
        if not paths:
            print(f"no figure matching {args.figure!r}")
            return 1

    total = 0
    for path in paths:
        defects = audit(path)
        if not defects:
            continue
        total += len(defects)
        counts: dict[str, int] = {}
        for defect in defects:
            counts[defect.kind] = counts.get(defect.kind, 0) + 1
        summary = ", ".join(f"{count} {kind}" for kind, count in sorted(counts.items()))
        print(f"\n{path.name}  ({summary})")
        for defect in defects:
            print(f"    [{defect.kind}] {defect.detail}")

    print(f"\n{total} defect(s) across {len(paths)} figure(s).")
    if args.check and total:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
