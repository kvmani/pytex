"""Structural guards for the canonical SVG figures in ``docs/figures/``.

Two classes of defect are pinned here because both are invisible in the SVG
source and glaring on screen:

1. **Stroke-scaled arrowheads.** SVG markers default to
   ``markerUnits="strokeWidth"``, multiplying the arrowhead by the stroke width
   of the line it terminates. A 12-unit head on a ``stroke-width="4"`` line
   renders at 48 units. Across the figure set this produced arrowheads occupying
   up to 125% of the line they annotated.
2. **Text that overflows or collides.** A caption wider than its figure, or two
   labels sharing a baseline, is a broken figure even though the file parses.

The generated figures (`scripts/generate_reference_frame_figures.py`) must pass
both checks strictly. Hand-authored figures must pass the marker check; their
text layout is checked only for the frame-related subset that has been
regenerated, with the remainder recorded as a known follow-on.
"""

from __future__ import annotations

import itertools
import math
import re
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"

#: Figures produced by scripts/generate_algorithm_figures.py.
ALGORITHM_FIGURES = (
    "or_determination_algorithm.svg",
    "variant_correspondence_algorithm.svg",
    "composite_saed_algorithm.svg",
    "saed_indexing_algorithm.svg",
)

#: Figures produced by scripts/generate_reference_frame_figures.py.
FRAME_FIGURES = (
    "reference_frame_catalog.svg",
    "sample_frame_rd_td_nd.svg",
    "reference_frames.svg",
    "reference_frames_vectors.svg",
    "orientation_mapping_semantics.svg",
    "active_passive_rotation.svg",
    "bunge_euler_geometry.svg",
    "hcp_reference_frame.svg",
)

#: Figures produced by scripts/generate_class_model_figures.py.
CLASS_MODEL_FIGURES = (
    "class_model_architecture.svg",
    "class_hierarchy.svg",
    "class_model_core.svg",
    "class_model_transformation.svg",
    "class_model_texture.svg",
    "class_model_ebsd.svg",
    "class_model_diffraction.svg",
    "class_model_tem.svg",
)

#: Every figure written by a generator, and therefore held to the strict layout
#: contract. Hand-authored figures pass the marker check only.
GENERATED_FIGURES = FRAME_FIGURES + ALGORITHM_FIGURES + CLASS_MODEL_FIGURES

# Text is measured with the library's own Helvetica advance-width table rather
# than an average character width. The average model is wrong by 20-30% for
# strings dominated by narrow lowercase or by capitals, which on a figure of
# side-by-side cards reports overlaps between labels that do not touch. One
# ruler is used everywhere: this guard, `scripts/audit_figure_text_layout.py`,
# and the figure generators themselves.
from pytex.plotting._svg_text import text_width  # noqa: E402

_TEXT_RE = re.compile(
    r'<text[^>]*\bx="(-?[0-9.]+)"[^>]*\by="(-?[0-9.]+)"[^>]*font-size="([0-9.]+)"[^>]*>'
    r"([^<]*)</text>"
)


def _figures() -> list[Path]:
    return sorted(FIGURES_DIR.glob("*.svg"))


def _markers(text: str) -> list[str]:
    return re.findall(r"<marker\b[^>]*>", text, re.DOTALL)


# ---------------------------------------------------------------------------
# Arrowhead scaling
# ---------------------------------------------------------------------------


def test_every_figure_marker_uses_absolute_units() -> None:
    """No figure may let its arrowheads scale with stroke width."""

    offenders = [
        f"{path.name}: {tag[:70]}"
        for path in _figures()
        for tag in _markers(path.read_text(encoding="utf-8"))
        if "markerUnits" not in tag
    ]
    assert not offenders, (
        "These markers scale with stroke width and will render oversized arrowheads. "
        "Run `python scripts/fix_svg_marker_units.py`. Offenders: " + "; ".join(offenders)
    )


def test_marker_fix_script_reports_clean() -> None:
    """The maintenance script's --check mode agrees with the test above."""

    import subprocess
    import sys

    completed = subprocess.run(
        [sys.executable, "scripts/fix_svg_marker_units.py", "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def _head_to_line_ratios(text: str) -> list[float]:
    """Arrowhead length divided by arrow length, for every marker-terminated line."""

    widths = {
        (re.search(r'id="([^"]+)"', tag) or [None, ""])[1]: float(
            (re.search(r'markerWidth="([0-9.]+)"', tag) or [None, "10"])[1]
        )
        for tag in _markers(text)
    }
    ratios: list[float] = []
    for match in re.finditer(
        r'<line[^>]*x1="(-?[0-9.]+)"[^>]*y1="(-?[0-9.]+)"[^>]*x2="(-?[0-9.]+)"'
        r'[^>]*y2="(-?[0-9.]+)"[^>]*marker-end="url\(#([^)]+)\)"',
        text,
    ):
        x1, y1, x2, y2, marker_id = match.groups()
        length = math.hypot(float(x2) - float(x1), float(y2) - float(y1))
        if length > 1.0:
            ratios.append(widths.get(marker_id, 10.0) / length)
    return ratios


def test_arrowheads_stay_small_relative_to_the_lines_they_terminate() -> None:
    """The typical arrowhead must not dominate the arrow it terminates.

    This is a **runaway detector**, not a style rule. The ratio alone cannot
    separate handsome from ugly, because it depends on the drawing scale: a
    small triad legitimately carries a head worth a third of its short axis,
    and those figures have been checked by eye. What it does catch is the
    failure mode that was actually present — markers scaling with stroke width,
    which drove medians to 0.62, 0.68 and 1.11, the last being an arrowhead
    longer than its own arrow.

    The bound is on each figure's **median** so one short connector stub cannot
    trip it.
    """

    failures: list[str] = []
    for path in _figures():
        ratios = _head_to_line_ratios(path.read_text(encoding="utf-8"))
        if not ratios:
            continue
        median = statistics.median(ratios)
        if median > 0.5:
            failures.append(f"{path.name}: median head/line = {median:.2f}")
    assert not failures, (
        "Arrowheads dominate their arrows in: "
        + "; ".join(failures)
        + ". Run `python scripts/fix_svg_marker_units.py`."
    )


# ---------------------------------------------------------------------------
# Document validity and accessibility
# ---------------------------------------------------------------------------


def test_every_figure_is_well_formed_xml() -> None:
    for path in _figures():
        try:
            ET.fromstring(path.read_text(encoding="utf-8"))
        except ET.ParseError as error:  # pragma: no cover - fails loudly when hit
            pytest.fail(f"{path.name} is not well-formed SVG: {error}")


@pytest.mark.parametrize("name", GENERATED_FIGURES)
def test_generated_figures_carry_title_and_description(name: str) -> None:
    """The style guide requires both; <desc> is also the accessible text."""

    namespace = "{http://www.w3.org/2000/svg}"
    root = ET.fromstring((FIGURES_DIR / name).read_text(encoding="utf-8"))
    title = root.find(f"{namespace}title")
    desc = root.find(f"{namespace}desc")
    assert title is not None and title.text, f"{name} has no <title>"
    assert desc is not None and desc.text, f"{name} has no <desc>"
    assert len(desc.text) > 60, f"{name} has an uninformative <desc>"


@pytest.mark.parametrize("name", GENERATED_FIGURES)
def test_generated_figures_use_the_canonical_font(name: str) -> None:
    text = (FIGURES_DIR / name).read_text(encoding="utf-8")
    assert 'font-family="Arial"' in text
    # No fallback stack: the layout metrics are Arial's, so a substituted face
    # would render text the figure was not laid out for. See the style guide.
    assert "Helvetica" not in text
    assert "serif" not in text


# ---------------------------------------------------------------------------
# Text layout
# ---------------------------------------------------------------------------


def _layout_defects(path: Path) -> tuple[list[str], list[str]]:
    """Estimated text overflows and same-baseline collisions in one figure."""

    text = path.read_text(encoding="utf-8")
    viewbox = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', text)
    if viewbox is None:
        return [], []
    width = float(viewbox.group(1))

    spans: list[tuple[float, float, float, str]] = []
    for match in _TEXT_RE.finditer(text):
        x, y, size, body = match.groups()
        content = body.strip()
        if not content:
            continue
        anchor_match = re.search(r'text-anchor="([a-z]+)"', match.group(0))
        anchor = anchor_match.group(1) if anchor_match else "start"
        estimated = text_width(content, float(size))
        left = float(x)
        if anchor == "middle":
            left -= estimated / 2.0
        elif anchor == "end":
            left -= estimated
        spans.append((left, left + estimated, float(y), content))

    overflows = [f"{s[3][:40]!r} by {s[1] - width:.0f}px" for s in spans if s[1] > width + 2.0]

    collisions: list[str] = []
    baselines: dict[int, list[tuple[float, float, float, str]]] = {}
    for span in spans:
        baselines.setdefault(round(span[2]), []).append(span)
    for group in baselines.values():
        group.sort()
        for first, second in itertools.pairwise(group):
            if second[0] < first[1] - 2.0:
                collisions.append(f"{first[3][:24]!r}/{second[3][:24]!r}")
    return overflows, collisions


@pytest.mark.parametrize("name", GENERATED_FIGURES)
def test_generated_figures_have_no_text_overflow_or_collisions(name: str) -> None:
    """Generated layout must fit: panels size themselves from estimated text width."""

    overflows, collisions = _layout_defects(FIGURES_DIR / name)
    assert not overflows, f"{name} text overflows the canvas: {'; '.join(overflows)}"
    assert not collisions, f"{name} has colliding labels: {'; '.join(collisions)}"


def test_generated_figure_list_matches_the_generator() -> None:
    """`FRAME_FIGURES` must not drift from what the frame script actually writes."""

    source = (REPO_ROOT / "scripts" / "generate_reference_frame_figures.py").read_text(
        encoding="utf-8"
    )
    written = set(re.findall(r'FIGURES_DIR / "([^"]+\.svg)"', source))
    assert written == set(FRAME_FIGURES), (
        "The generator writes a different figure set than this test tracks: "
        f"generator={sorted(written)}, test={sorted(FRAME_FIGURES)}"
    )


def test_algorithm_figure_list_matches_the_generator() -> None:
    """`ALGORITHM_FIGURES` must not drift from what the algorithm script writes."""

    source = (REPO_ROOT / "scripts" / "generate_algorithm_figures.py").read_text(encoding="utf-8")
    written = set(re.findall(r'"([a-z0-9_]+\.svg)": ', source))
    assert written == set(ALGORITHM_FIGURES), (
        "The generator writes a different figure set than this test tracks: "
        f"generator={sorted(written)}, test={sorted(ALGORITHM_FIGURES)}"
    )


def test_class_model_figure_list_matches_the_generator() -> None:
    """`CLASS_MODEL_FIGURES` must not drift from what the atlas generator writes."""

    import sys

    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from generate_class_model_figures import build_figures

    assert set(build_figures()) == set(CLASS_MODEL_FIGURES)


def test_algorithm_figures_are_deterministic() -> None:
    """Regenerating must reproduce the committed bytes exactly.

    A generated asset whose bytes move between runs cannot be checked for drift,
    which defeats the reason these figures are generated rather than drawn. This
    caught `builtins.hash` being used for marker ids, which Python randomizes per
    process.
    """

    import sys

    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from generate_algorithm_figures import (
        composite_saed_figure,
        or_determination_figure,
        saed_indexing_figure,
        variant_correspondence_figure,
    )

    produced = {
        "or_determination_algorithm.svg": or_determination_figure(),
        "variant_correspondence_algorithm.svg": variant_correspondence_figure(),
        "composite_saed_algorithm.svg": composite_saed_figure(),
        "saed_indexing_algorithm.svg": saed_indexing_figure(),
    }
    for name, svg in produced.items():
        committed = (FIGURES_DIR / name).read_text(encoding="utf-8")
        assert svg == committed, (
            f"{name} differs from the committed asset. Run "
            "`python scripts/generate_algorithm_figures.py`."
        )


# ---------------------------------------------------------------------------
# Text layout across every figure, not only the generated ones
# ---------------------------------------------------------------------------


def _layout_defects_by_kind() -> dict[str, list[str]]:
    """Every layout defect the auditor finds, grouped by kind."""

    import sys

    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from audit_figure_text_layout import audit

    grouped: dict[str, list[str]] = {}
    for path in _figures():
        for defect in audit(path):
            grouped.setdefault(defect.kind, []).append(f"{defect.figure}: {defect.detail}")
    return grouped


def test_no_figure_text_runs_outside_its_box_or_canvas() -> None:
    """The reported defect: captions painted across their panel border.

    ``box`` and ``canvas`` overflows are held at zero. They are the ones a reader
    sees as text spilling out of a card or being cut off at the edge, and every
    one of them has been repaired — 5 remained after the automated passes and
    were fixed individually.
    """

    grouped = _layout_defects_by_kind()
    offenders = grouped.get("box", []) + grouped.get("canvas", [])
    assert not offenders, (
        "Text overflows its box or the canvas in: "
        + "; ".join(offenders)
        + ". Run `python scripts/audit_figure_text_layout.py` to inspect, then "
        "`fix_svg_text_overflow.py` / `fix_svg_title_wraps.py` to repair."
    )


#: Overlaps that remain in dense hand-authored teaching figures. These are
#: tracked rather than zeroed because each needs an individual layout decision
#: in a figure whose geometry is meaningful (an axis label beside a formula, a
#: callout over a sphere). The count is a ceiling: it must not grow.
KNOWN_OVERLAP_CEILING = 20


def test_figure_text_overlaps_do_not_grow() -> None:
    """Collisions and text-over-panel are capped so repairs cannot regress.

    The starting point of this program was 86 layout defects measured with
    accurate Helvetica metrics. Captions overflowing cards and titles overflowing
    cards are now zero; what remains are overlaps inside dense diagrams, held
    under a ceiling so no change can quietly add more.
    """

    grouped = _layout_defects_by_kind()
    remaining = grouped.get("collision", []) + grouped.get("over-card", [])
    assert len(remaining) <= KNOWN_OVERLAP_CEILING, (
        f"Figure text overlaps rose to {len(remaining)}, above the tracked ceiling of "
        f"{KNOWN_OVERLAP_CEILING}: " + "; ".join(remaining)
    )
