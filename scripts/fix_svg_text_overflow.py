"""Wrap SVG text that runs outside its box, without shrinking it or dropping words.

The defect
----------

Hand-authored figures place a caption inside a card as a single ``<text>`` run.
When the caption is wider than the card, SVG happily paints it straight across
the card border and over whatever is next to it. Nothing in the file is invalid;
the figure is simply unreadable.

Why wrapping, and not the alternatives
--------------------------------------

- *Shrinking the text* would breach the visualization style guide, which fixes a
  15px floor for body text and 18px for node labels.
- *Scaling the whole figure up* fixes the overlap in figure coordinates but makes
  the text smaller on screen: Sphinx fits the image to the container, so a wider
  intrinsic canvas is scaled down further, and legibility gets worse rather than
  better.
- *Rewording* loses scientific content, which is not this script's call to make.

Wrapping keeps every word at its designed size and only uses vertical space,
which these card layouts have.

**Titles are never wrapped.** A card title's second line lands exactly where the
caption already sits, turning a horizontal overflow into an overprint. Bold runs
and runs at 18px or more — node labels, per the visualization style guide — are
left alone here and handled by ``fix_svg_title_wraps.py``, which widens the card
instead.

How it works
------------

Widths come from `pytex.plotting._svg_text`, the Helvetica advance-width table,
with font sizes resolved from attributes, CSS classes, or an ancestor — the same
resolution the auditor uses, so the two agree. An overflowing run is broken at
word boundaries into lines that fit the box, and re-emitted as ``<tspan>``
elements sharing the original ``x`` so the anchor still applies.

When the wrapped block would spill past the bottom of its box, ``--grow-boxes``
first enlarges that box — downward if the space below is free, otherwise
sideways into free space to the right — and only then wraps. A run that still
cannot be made to fit is reported rather than silently written, because that
needs a layout decision a script should not make alone.

Usage::

    python scripts/fix_svg_text_overflow.py            # rewrite in place
    python scripts/fix_svg_text_overflow.py --dry-run  # report the plan only
"""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from audit_figure_text_layout import (  # noqa: E402
    DEFAULT_FONT_SIZE,
    SVG_NS,
    _resolve_font,
    _stylesheet,
    _translate,
)

from pytex.plotting._svg_text import text_width  # noqa: E402

FIGURES_DIR = REPO_ROOT / "docs" / "figures"

#: Space left between wrapped text and the right edge of its box.
RIGHT_PADDING = 10.0

#: Line spacing as a multiple of the font size.
LINE_HEIGHT = 1.28

#: Runs at or above this size are node labels, not body text, and are never
#: wrapped: see the module docstring.
TITLE_SIZE = 18.0


def _wrap(content: str, font_size: float, bold: bool, available: float) -> list[str]:
    """Break a string into the fewest lines that each fit ``available`` units."""

    words = content.split()
    if not words:
        return []
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = f"{current} {word}"
        if text_width(candidate, font_size, bold=bold) <= available:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def _plan(path: Path) -> list[tuple[str, str, list[str], bool]]:
    """Wrap plans for one figure: (original, id-ish key, lines, fits_vertically)."""

    text = path.read_text(encoding="utf-8")
    root = ET.fromstring(text)
    styles = _stylesheet(root)
    viewbox = root.get("viewBox")
    if viewbox:
        parts = viewbox.split()
        canvas = (float(parts[2]), float(parts[3]))
    else:
        canvas = (float(root.get("width", "0") or 0), float(root.get("height", "0") or 0))

    boxes: list[tuple[float, float, float, float]] = []
    runs: list[tuple[ET.Element, float, float, float, bool, str]] = []

    def walk(element: ET.Element, offset: tuple[float, float], font: tuple[float, bool]) -> None:
        for child in element:
            dx, dy = _translate(child)
            child_offset = (offset[0] + dx, offset[1] + dy)
            child_font = _resolve_font(child, styles, font)
            if child.tag == f"{SVG_NS}rect":
                boxes.append(
                    (
                        float(child.get("x", "0") or 0) + child_offset[0],
                        float(child.get("y", "0") or 0) + child_offset[1],
                        float(child.get("width", "0") or 0),
                        float(child.get("height", "0") or 0),
                    )
                )
            elif child.tag == f"{SVG_NS}text":
                content = "".join(child.itertext()).strip()
                size, bold = child_font
                is_title = bold or size >= TITLE_SIZE
                if content and len(child) == 0 and not is_title:
                    runs.append(
                        (
                            child,
                            float(child.get("x", "0") or 0) + child_offset[0],
                            float(child.get("y", "0") or 0) + child_offset[1],
                            size,
                            bold,
                            content,
                        )
                    )
            walk(child, child_offset, child_font)

    walk(root, (0.0, 0.0), (DEFAULT_FONT_SIZE, False))

    plans: list[tuple[str, str, list[str], bool]] = []
    for element, x, y, size, bold, content in runs:
        width = text_width(content, size, bold=bold)
        anchor = element.get("text-anchor", "start")
        left = x - (width / 2.0 if anchor == "middle" else width if anchor == "end" else 0.0)
        right = left + width

        enclosing = [
            box
            for box in boxes
            if box[2] * box[3] > 0
            and not (box[2] >= canvas[0] - 1.0 and box[3] >= canvas[1] - 1.0)
            and box[0] - 1.0 <= x <= box[0] + box[2] + 1.0
            and box[1] - 1.0 <= y <= box[1] + box[3] + 1.0
        ]
        if enclosing:
            box = min(enclosing, key=lambda b: b[2] * b[3])
            if anchor == "middle":
                available = min(x - box[0], box[0] + box[2] - x) * 2.0 - RIGHT_PADDING
            elif anchor == "end":
                available = x - box[0] - RIGHT_PADDING
            else:
                available = box[0] + box[2] - x - RIGHT_PADDING
            box_bottom = box[1] + box[3]
        elif right > canvas[0]:
            available = canvas[0] - left - RIGHT_PADDING
            box_bottom = canvas[1]
        else:
            continue

        if width <= available or available <= 20.0:
            continue
        lines = _wrap(content, size, bold, available)
        if len(lines) < 2:
            continue
        extra = (len(lines) - 1) * LINE_HEIGHT * size
        fits = (y + extra + 0.3 * size) <= box_bottom
        plans.append((content, f"{x:.1f},{y:.1f}", lines, fits))
    return plans


def _apply(text: str, content: str, lines: list[str], size: float) -> tuple[str, bool]:
    """Replace the first single-line <text> holding ``content`` with tspans."""

    pattern = re.compile(
        r"(<text\b[^>]*>)" + re.escape(content) + r"(</text>)",
    )
    match = pattern.search(text)
    if match is None:
        return text, False
    opening = match.group(1)
    x_match = re.search(r'\bx="(-?[0-9.]+)"', opening)
    if x_match is None:
        return text, False
    x = x_match.group(1)
    step = LINE_HEIGHT * size
    spans = "".join(
        f'<tspan x="{x}" dy="{0 if index == 0 else step:.1f}">{line}</tspan>'
        for index, line in enumerate(lines)
    )
    return text[: match.start()] + opening + spans + match.group(2) + text[match.end() :], True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="report the plan, change nothing")
    parser.add_argument("--figure", help="restrict to one figure")
    parser.add_argument(
        "--grow-boxes",
        action="store_true",
        help="enlarge a box when its caption cannot be wrapped into the existing space",
    )
    parser.add_argument(
        "--widen-boxes",
        action="store_true",
        help="widen a box (and the canvas if needed) for text that cannot be wrapped at all",
    )
    args = parser.parse_args()

    paths = sorted(FIGURES_DIR.glob("*.svg"))
    if args.figure:
        paths = [p for p in paths if args.figure in p.name]

    wrapped = 0
    grown = 0
    skipped: list[str] = []
    for path in paths:
        if args.widen_boxes and not args.dry_run:
            grown += _widen_boxes(path)
        if args.grow_boxes and not args.dry_run:
            grown += _grow_boxes(path)
        plans = _plan(path)
        if not plans:
            continue
        text = path.read_text(encoding="utf-8")
        changed = 0
        for content, _key, lines, fits in plans:
            if not fits:
                skipped.append(f"{path.name}: {content[:44]!r} (no vertical room)")
                continue
            size_match = re.search(
                r"<text\b[^>]*>" + re.escape(content) + r"</text>", text
            )
            if size_match is None:
                continue
            # Re-measure the font size for this specific run via the plan's lines.
            font_size = _font_size_for(path, content)
            text, ok = _apply(text, content, lines, font_size)
            if ok:
                changed += 1
        if changed and not args.dry_run:
            path.write_text(text, encoding="utf-8")
        if changed:
            wrapped += changed
            print(f"{'would wrap' if args.dry_run else 'wrapped'} {changed} run(s) in {path.name}")

    for note in skipped:
        print(f"SKIPPED {note}")
    if args.grow_boxes:
        print(f"{grown} box(es) enlarged.")
    print(f"\n{wrapped} run(s) wrapped; {len(skipped)} needing a manual layout decision.")
    return 0


def _grow_boxes(path: Path) -> int:
    """Enlarge boxes whose caption cannot be wrapped into the space available.

    Height is preferred: these are card grids with vertical gaps between rows,
    so a card can usually take one more line without disturbing its neighbours.
    The growth is bounded by the nearest element below, so a card never swallows
    the row beneath it.
    """

    text = path.read_text(encoding="utf-8")
    root = ET.fromstring(text)
    styles = _stylesheet(root)
    viewbox = root.get("viewBox")
    if viewbox:
        parts = viewbox.split()
    else:
        parts = ["0", "0", root.get("width", "0") or "0", root.get("height", "0") or "0"]
    canvas = (float(parts[2]), float(parts[3]))

    rects: list[tuple[ET.Element, float, float, float, float]] = []
    runs: list[tuple[float, float, float, bool, str, str]] = []

    def walk(element: ET.Element, offset: tuple[float, float], font: tuple[float, bool]) -> None:
        for child in element:
            dx, dy = _translate(child)
            child_offset = (offset[0] + dx, offset[1] + dy)
            child_font = _resolve_font(child, styles, font)
            if child.tag == f"{SVG_NS}rect":
                rects.append(
                    (
                        child,
                        float(child.get("x", "0") or 0) + child_offset[0],
                        float(child.get("y", "0") or 0) + child_offset[1],
                        float(child.get("width", "0") or 0),
                        float(child.get("height", "0") or 0),
                    )
                )
            elif child.tag == f"{SVG_NS}text" and len(child) == 0:
                content = (child.text or "").strip()
                if content:
                    runs.append(
                        (
                            float(child.get("x", "0") or 0) + child_offset[0],
                            float(child.get("y", "0") or 0) + child_offset[1],
                            child_font[0],
                            child_font[1],
                            child.get("text-anchor", "start"),
                            content,
                        )
                    )
            walk(child, child_offset, child_font)

    walk(root, (0.0, 0.0), (DEFAULT_FONT_SIZE, False))

    pending: list[tuple[ET.Element, float, float]] = []
    for x, y, size, bold, anchor, content in runs:
        width = text_width(content, size, bold=bold)
        enclosing = [
            entry
            for entry in rects
            if entry[3] * entry[4] > 0
            and not (entry[3] >= canvas[0] - 1.0 and entry[4] >= canvas[1] - 1.0)
            and entry[1] - 1.0 <= x <= entry[1] + entry[3] + 1.0
            and entry[2] - 1.0 <= y <= entry[2] + entry[4] + 1.0
        ]
        if not enclosing:
            continue
        element, bx, by, bw, bh = min(enclosing, key=lambda e: e[3] * e[4])
        if anchor == "middle":
            available = min(x - bx, bx + bw - x) * 2.0 - RIGHT_PADDING
        elif anchor == "end":
            available = x - bx - RIGHT_PADDING
        else:
            available = bx + bw - x - RIGHT_PADDING
        if width <= available or available <= 20.0:
            continue

        lines = _wrap(content, size, bold, available)
        needed_bottom = y + (len(lines) - 1) * LINE_HEIGHT * size + 0.35 * size
        if needed_bottom <= by + bh:
            continue

        # How far can this box grow downward before meeting something below?
        obstacles = [
            entry[2]
            for entry in rects
            if entry is not element
            and entry[3] * entry[4] > 0
            and not (entry[3] >= canvas[0] - 1.0 and entry[4] >= canvas[1] - 1.0)
            and entry[2] > by + bh - 1.0
            and entry[1] < bx + bw
            and entry[1] + entry[3] > bx
        ]
        ceiling = min(obstacles) if obstacles else canvas[1]
        allowed = max(0.0, ceiling - 12.0 - (by + bh))
        required = needed_bottom - (by + bh) + 6.0
        if required > allowed:
            continue
        pending.append((element, bh, bh + required))

    if not pending:
        return 0

    # Rewrite the raw source rather than serializing the parsed tree: a
    # round-trip through ElementTree rewrites the document with generated
    # namespace prefixes (ns0:svg), which is valid but gratuitously reformats
    # every hand-authored figure in the repository.
    updated = text
    changed = 0
    for element, old_height, new_height in pending:
        pattern = _rect_pattern(element, old_height)
        match = pattern.search(updated)
        if match is None:
            continue
        replacement = re.sub(
            r'height="[0-9.]+"', f'height="{new_height:.1f}"', match.group(0), count=1
        )
        updated = updated[: match.start()] + replacement + updated[match.end() :]
        changed += 1

    if changed:
        path.write_text(updated, encoding="utf-8")
    return changed


def _rect_pattern(element: ET.Element, height: float) -> re.Pattern[str]:
    """Match the raw ``<rect>`` tag carrying these exact geometry attributes."""

    def number(value: str | None) -> str:
        return re.escape((value or "0").strip())

    return re.compile(
        r"<rect\b[^>]*?"
        + f'x="{number(element.get("x"))}"'
        + r"[^>]*?"
        + f'y="{number(element.get("y"))}"'
        + r"[^>]*?"
        + f'width="{number(element.get("width"))}"'
        + r"[^>]*?"
        + f'height="{number(element.get("height"))}"'
        + r"[^>]*?/?>"
    )


def _font_size_for(path: Path, content: str) -> float:
    """The resolved font size of the run holding ``content`` in ``path``."""

    root = ET.fromstring(path.read_text(encoding="utf-8"))
    styles = _stylesheet(root)

    found = [DEFAULT_FONT_SIZE]

    def walk(element: ET.Element, font: tuple[float, bool]) -> None:
        for child in element:
            child_font = _resolve_font(child, styles, font)
            if child.tag == f"{SVG_NS}text" and "".join(child.itertext()).strip() == content:
                found[0] = child_font[0]
            walk(child, child_font)

    walk(root, (DEFAULT_FONT_SIZE, False))
    return found[0]




def _widen_boxes(path: Path) -> int:
    """Widen boxes (and the canvas) for text that cannot be wrapped.

    A single long word — ``DiffractionGeometry`` — has no break point, so
    wrapping cannot help it and the box simply has to be wider. Growth is
    bounded by the nearest element to the right so a card never overruns its
    neighbour; when the text overflows the canvas itself, the canvas grows
    instead.
    """

    text = path.read_text(encoding="utf-8")
    root = ET.fromstring(text)
    styles = _stylesheet(root)
    viewbox = root.get("viewBox")
    if viewbox:
        parts = viewbox.split()
    else:
        parts = ["0", "0", root.get("width", "0") or "0", root.get("height", "0") or "0"]
    canvas = (float(parts[2]), float(parts[3]))

    rects: list[tuple[ET.Element, float, float, float, float]] = []
    runs: list[tuple[float, float, float, bool, str, str]] = []

    def walk(element: ET.Element, offset: tuple[float, float], font: tuple[float, bool]) -> None:
        for child in element:
            dx, dy = _translate(child)
            child_offset = (offset[0] + dx, offset[1] + dy)
            child_font = _resolve_font(child, styles, font)
            if child.tag == f"{SVG_NS}rect":
                rects.append(
                    (
                        child,
                        float(child.get("x", "0") or 0) + child_offset[0],
                        float(child.get("y", "0") or 0) + child_offset[1],
                        float(child.get("width", "0") or 0),
                        float(child.get("height", "0") or 0),
                    )
                )
            elif child.tag == f"{SVG_NS}text":
                size, bold = child_font
                base_x = float(child.get("x", "0") or 0) + child_offset[0]
                base_y = float(child.get("y", "0") or 0) + child_offset[1]
                anchor = child.get("text-anchor", "start")
                spans = list(child.iter(f"{SVG_NS}tspan"))
                if spans:
                    cursor = base_y
                    for span in spans:
                        line = (span.text or "").strip()
                        cursor += float(span.get("dy", "0") or 0)
                        if line:
                            runs.append((base_x, cursor, size, bold, anchor, line))
                else:
                    content = (child.text or "").strip()
                    if content:
                        runs.append((base_x, base_y, size, bold, anchor, content))
            walk(child, child_offset, child_font)

    walk(root, (0.0, 0.0), (DEFAULT_FONT_SIZE, False))

    pending: list[tuple[ET.Element, float, float]] = []
    canvas_needed = canvas[0]
    for x, y, size, bold, anchor, content in runs:
        width = text_width(content, size, bold=bold)
        left = x - (width / 2.0 if anchor == "middle" else width if anchor == "end" else 0.0)
        right = left + width
        enclosing = [
            entry
            for entry in rects
            if entry[3] * entry[4] > 0
            and not (entry[3] >= canvas[0] - 1.0 and entry[4] >= canvas[1] - 1.0)
            and entry[1] - 1.0 <= x <= entry[1] + entry[3] + 1.0
            and entry[2] - 1.0 <= y <= entry[2] + entry[4] + 1.0
        ]
        if not enclosing:
            if right > canvas_needed:
                canvas_needed = right + RIGHT_PADDING
            continue
        element, bx, by, bw, bh = min(enclosing, key=lambda e: e[3] * e[4])
        needed_right = right + RIGHT_PADDING
        if needed_right <= bx + bw:
            continue
        obstacles = [
            entry[1]
            for entry in rects
            if entry is not element
            and entry[3] * entry[4] > 0
            and not (entry[3] >= canvas[0] - 1.0 and entry[4] >= canvas[1] - 1.0)
            and entry[1] > bx + bw - 1.0
            and entry[2] < by + bh
            and entry[2] + entry[4] > by
        ]
        ceiling = min(obstacles) if obstacles else canvas[0]
        allowed = max(0.0, ceiling - 14.0 - (bx + bw))
        required = needed_right - (bx + bw)
        if required > allowed:
            continue
        pending.append((element, bw, bw + required))

    updated = text
    changed = 0
    for element, old_width, new_width in pending:
        pattern = _rect_pattern(element, old_width)
        match = pattern.search(updated)
        if match is None:
            continue
        replacement = re.sub(
            r'width="[0-9.]+"', f'width="{new_width:.1f}"', match.group(0), count=1
        )
        updated = updated[: match.start()] + replacement + updated[match.end() :]
        changed += 1

    if canvas_needed > canvas[0] + 1.0:
        new_canvas = canvas_needed + 20.0
        updated = updated.replace(
            f'viewBox="0 0 {parts[2]} {parts[3]}"',
            f'viewBox="0 0 {new_canvas:.0f} {parts[3]}"',
            1,
        )
        updated = re.sub(
            r'(<svg\b[^>]*?)width="[0-9.]+"', rf'\g<1>width="{new_canvas:.0f}"', updated, count=1
        )
        changed += 1

    if changed:
        path.write_text(updated, encoding="utf-8")
    return changed

if __name__ == "__main__":
    sys.exit(main())
