"""Un-wrap card titles and widen their cards instead.

Wrapping is the right repair for a caption, which has a blank line beneath it.
It is the wrong repair for a card *title*: the second line lands exactly where
the caption already sits, turning a horizontal overflow into an overprint. A
title is short, so the honest fix is to make the card wide enough for it.

This pass finds every wrapped title — a bold run, or one at 18px or more, which
is what the visualization style guide calls a node label — restores it to a
single line, and widens the enclosing card to fit, bounded by the nearest card
to the right so nothing overruns its neighbour.

Run it after ``fix_svg_text_overflow.py``; ``audit_figure_text_layout.py``
should then report clean.

Usage::

    python scripts/fix_svg_title_wraps.py
    python scripts/fix_svg_title_wraps.py --dry-run
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

#: A run at or above this size counts as a node label rather than body text.
TITLE_SIZE = 18.0

#: Space kept between a title and the right edge of its card.
CARD_PADDING = 16.0


def _wrapped_titles(path: Path) -> list[tuple[str, str, float, float, float]]:
    """(raw_block, joined_title, x, y, size) for each wrapped title run."""

    root = ET.fromstring(path.read_text(encoding="utf-8"))
    styles = _stylesheet(root)
    found: list[tuple[str, str, float, float, float]] = []

    def walk(element: ET.Element, offset: tuple[float, float], font: tuple[float, bool]) -> None:
        for child in element:
            dx, dy = _translate(child)
            child_offset = (offset[0] + dx, offset[1] + dy)
            child_font = _resolve_font(child, styles, font)
            if child.tag == f"{SVG_NS}text":
                spans = [s for s in child.iter(f"{SVG_NS}tspan") if s.get("x") is not None]
                size, bold = child_font
                if len(spans) > 1 and (bold or size >= TITLE_SIZE):
                    title = " ".join(
                        part for span in spans if (part := (span.text or "").strip())
                    )
                    found.append(
                        (
                            ET.tostring(child, encoding="unicode"),
                            title,
                            float(child.get("x", "0") or 0) + child_offset[0],
                            float(child.get("y", "0") or 0) + child_offset[1],
                            size,
                        )
                    )
            walk(child, child_offset, child_font)

    walk(root, (0.0, 0.0), (DEFAULT_FONT_SIZE, False))
    return found


def _cards(path: Path) -> list[tuple[float, float, float, float, str]]:
    """Absolute (x, y, w, h, raw_tag) for every rect in the figure."""

    text = path.read_text(encoding="utf-8")
    out: list[tuple[float, float, float, float, str]] = []
    for match in re.finditer(
        r'<rect\b[^>]*?x="(-?[0-9.]+)"[^>]*?y="(-?[0-9.]+)"[^>]*?'
        r'width="([0-9.]+)"[^>]*?height="([0-9.]+)"[^>]*?/?>',
        text,
    ):
        x, y, w, h = (float(v) for v in match.groups())
        out.append((x, y, w, h, match.group(0)))
    return out


def fix_figure(path: Path, *, dry_run: bool) -> tuple[int, int]:
    """Un-wrap titles and widen cards. Returns (titles restored, cards widened)."""

    titles = _wrapped_titles(path)
    if not titles:
        return 0, 0
    text = path.read_text(encoding="utf-8")
    cards = _cards(path)
    canvas_match = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', text)
    canvas_w = float(canvas_match.group(1)) if canvas_match else 0.0

    restored = 0
    widened = 0
    for _raw, title, x, y, size in titles:
        # Restore the single line by rewriting the tspan block in the raw source.
        pattern = re.compile(
            r'<text\b([^>]*\bx="' + re.escape(f"{x:g}") + r'"[^>]*)>'
            r"(?:<tspan[^>]*>[^<]*</tspan>)+</text>"
        )
        match = pattern.search(text)
        if match is None:
            # x may be written with a decimal point; retry on the y coordinate.
            pattern = re.compile(
                r'<text\b([^>]*\by="' + re.escape(f"{y:g}") + r'"[^>]*)>'
                r"(?:<tspan[^>]*>[^<]*</tspan>)+</text>"
            )
            match = pattern.search(text)
        if match is None:
            continue
        text = text[: match.start()] + f"<text{match.group(1)}>{title}</text>" + text[match.end() :]
        restored += 1

        # Widen the card that holds it, bounded by the next card to the right.
        containing = [
            card
            for card in cards
            if card[2] * card[3] > 2000.0
            and not (card[2] >= canvas_w - 1.0)
            and card[0] - 1.0 <= x <= card[0] + card[2] + 1.0
            and card[1] - 1.0 <= y <= card[1] + card[3] + 1.0
        ]
        if not containing:
            continue
        cx, cy, cw, ch, raw_tag = min(containing, key=lambda c: c[2] * c[3])
        needed = x - cx + text_width(title, size, bold=True) + CARD_PADDING
        if needed <= cw:
            continue
        obstacles = [
            card[0]
            for card in cards
            if card is not None
            and card[4] != raw_tag
            and card[2] * card[3] > 2000.0
            and card[0] > cx + cw - 1.0
            and card[1] < cy + ch
            and card[1] + card[3] > cy
        ]
        ceiling = min(obstacles) if obstacles else canvas_w
        allowed = cw + max(0.0, ceiling - 14.0 - (cx + cw))
        new_width = min(needed, allowed)
        if new_width <= cw + 0.5:
            continue
        replacement = re.sub(
            r'width="[0-9.]+"', f'width="{new_width:.1f}"', raw_tag, count=1
        )
        text = text.replace(raw_tag, replacement, 1)
        widened += 1

    if (restored or widened) and not dry_run:
        path.write_text(text, encoding="utf-8")
    return restored, widened


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--figure")
    args = parser.parse_args()

    paths = sorted(FIGURES_DIR.glob("*.svg"))
    if args.figure:
        paths = [p for p in paths if args.figure in p.name]

    total_titles = total_cards = 0
    for path in paths:
        restored, widened = fix_figure(path, dry_run=args.dry_run)
        if restored or widened:
            print(f"{path.name}: {restored} title(s) un-wrapped, {widened} card(s) widened")
            total_titles += restored
            total_cards += widened
    print(f"\n{total_titles} title(s) restored; {total_cards} card(s) widened.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
