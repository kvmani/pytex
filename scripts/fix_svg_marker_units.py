"""Stop hand-authored SVG arrowheads from scaling with their line's stroke width.

The problem
-----------

SVG markers default to ``markerUnits="strokeWidth"``, which multiplies the
marker's geometry by the stroke width of the line it terminates. A figure that
declares a 12-unit arrowhead and draws a ``stroke-width="4"`` line therefore
renders a **48-unit** arrowhead. Across ``docs/figures/`` this produced
arrowheads occupying 11% to 125% of the lines they terminate — in the worst case
the head was longer than the whole arrow, and the triads it annotated were
unreadable.

The fix
-------

For every marker, switch to ``markerUnits="userSpaceOnUse"`` (absolute units) and
pre-multiply the marker's own geometry so the arrowhead keeps a sensible size:

``scale = min(median stroke width using this marker,
              cap so the head is at most 25% of the median line length)``

The first term preserves each figure's existing visual weight — a diagram whose
arrows already looked right keeps looking the same. The second term rescues the
figures where the head had run away, by bounding it against the geometry it
actually annotates.

Usage::

    python scripts/fix_svg_marker_units.py            # rewrite in place
    python scripts/fix_svg_marker_units.py --check    # report only, exit 1 if any

``--check`` is what CI and `tests/unit/test_figure_markers.py` use: a figure that
reintroduces a stroke-scaled marker fails rather than silently shipping.

Generated figures (`scripts/generate_reference_frame_figures.py`) already emit
``markerUnits="userSpaceOnUse"`` and are left untouched.
"""

from __future__ import annotations

import argparse
import math
import re
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"

#: An arrowhead may not exceed this fraction of the line it terminates.
MAX_HEAD_FRACTION = 0.25

#: Fallback stroke width when a marker's users declare none.
DEFAULT_STROKE = 2.0

_MARKER_RE = re.compile(r"<marker\b[^>]*>", re.DOTALL)
_ATTR_RE = re.compile(r'(\b(?:markerWidth|markerHeight|refX|refY)=")([0-9.]+)(")')
_NUMBER_RE = re.compile(r"-?\d+(?:\.\d+)?")


def _marker_id(tag: str) -> str | None:
    match = re.search(r'id="([^"]+)"', tag)
    return match.group(1) if match else None


def _marker_width(tag: str) -> float:
    match = re.search(r'markerWidth="([0-9.]+)"', tag)
    return float(match.group(1)) if match else 10.0


def _users(text: str, marker_id: str) -> tuple[list[float], list[float]]:
    """Stroke widths and line lengths of the elements using ``marker_id``."""

    strokes: list[float] = []
    lengths: list[float] = []
    for element in re.finditer(r"<(?:line|path|polyline)\b[^>]*>", text, re.DOTALL):
        tag = element.group(0)
        if f"url(#{marker_id})" not in tag:
            continue
        stroke = re.search(r'stroke-width="([0-9.]+)"', tag)
        strokes.append(float(stroke.group(1)) if stroke else DEFAULT_STROKE)
        coords = re.search(
            r'x1="(-?[0-9.]+)"[^>]*y1="(-?[0-9.]+)"[^>]*x2="(-?[0-9.]+)"[^>]*y2="(-?[0-9.]+)"',
            tag,
        )
        if coords:
            x1, y1, x2, y2 = (float(v) for v in coords.groups())
            lengths.append(math.hypot(x2 - x1, y2 - y1))
    return strokes, lengths


def _scale_for(text: str, tag: str) -> float:
    """The factor to pre-multiply a marker's geometry by."""

    marker_id = _marker_id(tag)
    if marker_id is None:
        return DEFAULT_STROKE
    strokes, lengths = _users(text, marker_id)
    scale = statistics.median(strokes) if strokes else DEFAULT_STROKE
    if lengths:
        cap = MAX_HEAD_FRACTION * statistics.median(lengths) / max(_marker_width(tag), 1e-6)
        scale = min(scale, cap)
    return max(scale, 0.5)


def _scale_marker_block(block: str, scale: float) -> str:
    """Multiply a marker's declared size, reference point and path data."""

    def scale_attr(match: re.Match[str]) -> str:
        return f"{match.group(1)}{float(match.group(2)) * scale:.2f}{match.group(3)}"

    scaled = _ATTR_RE.sub(scale_attr, block)

    def scale_path(match: re.Match[str]) -> str:
        data = _NUMBER_RE.sub(lambda n: f"{float(n.group(0)) * scale:.2f}", match.group(2))
        return f"{match.group(1)}{data}{match.group(3)}"

    return re.sub(r'(\bd=")([^"]*)(")', scale_path, scaled)


def fix_text(text: str) -> tuple[str, int]:
    """Return the corrected SVG text and the number of markers changed."""

    changed = 0
    out = text
    for tag in _MARKER_RE.findall(text):
        if "markerUnits" in tag:
            continue
        marker_id = _marker_id(tag)
        if marker_id is None:
            continue
        start = out.index(tag)
        end = out.index("</marker>", start) + len("</marker>")
        block = out[start:end]
        scale = _scale_for(text, tag)
        scaled = _scale_marker_block(block, scale)
        scaled = scaled.replace("<marker ", '<marker markerUnits="userSpaceOnUse" ', 1)
        out = out[:start] + scaled + out[end:]
        changed += 1
    return out, changed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="report offending figures without rewriting; exit 1 if any are found",
    )
    args = parser.parse_args()

    offenders: list[str] = []
    rewritten: list[str] = []
    for path in sorted(FIGURES_DIR.glob("*.svg")):
        text = path.read_text(encoding="utf-8")
        if "<marker" not in text:
            continue
        fixed, changed = fix_text(text)
        if not changed:
            continue
        relative = str(path.relative_to(REPO_ROOT))
        if args.check:
            offenders.append(relative)
            continue
        path.write_text(fixed, encoding="utf-8")
        rewritten.append(f"{relative} ({changed} marker(s))")

    if args.check:
        for name in offenders:
            print(f"stroke-scaled marker in {name}")
        if offenders:
            print(
                f"\n{len(offenders)} figure(s) declare markers without "
                'markerUnits="userSpaceOnUse". Run: python scripts/fix_svg_marker_units.py'
            )
            return 1
        print("All figure markers use absolute units.")
        return 0

    for name in rewritten:
        print(f"fixed {name}")
    print(f"{len(rewritten)} figure(s) updated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
