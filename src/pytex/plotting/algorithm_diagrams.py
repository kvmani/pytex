"""Generated flow-sheet SVG for the library's scientific algorithms.

An algorithm page needs a figure that shows the *shape* of the computation:
what enters, which stages transform it, where a decision can reject a candidate,
and what leaves. This module builds those flow sheets from a declarative
description, so the documentation figures share one visual language with the
reference-frame diagrams and cannot drift into hand-drawn approximations.

Why a lane layout
-----------------

The style guide forbids long single-column flowsheets: content that has phases
should show its phases. `algorithm_flow_svg` therefore lays stages out in
**lanes** — a row per phase, cards within a row — so a reader sees the algorithm's
structure before reading a single label. Side notes attach to a stage without
entering the flow, which is where constraints, tolerances and failure modes go:
they belong beside the step they govern, not in a separate legend.

Colour carries role (`input`, `compute`, `decision`, `reject`, `output`) using
the canonical tokens, but never carries information the label does not also
carry, so the figure survives grayscale printing and colour-vision deficiency.

See also
--------
`pytex.plotting.svg_primitives` : the shared cards, arrows, and document wrapper.
`scripts/generate_algorithm_figures.py` : writes these into `docs/figures/`.
`docs/standards/visualization_style_guide.md` : the canonical tokens.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from pytex.plotting.svg_primitives import (
    ACCENT,
    MUTED,
    ROLE_INK,
    ROLE_TINTS,
    VIOLET,
    arrow_marker,
    callout,
    card,
    document,
    header_width,
    text,
    wrap_text,
)

__all__ = [
    "AlgorithmStage",
    "SideNote",
    "algorithm_flow_svg",
]

_ROLES = tuple(ROLE_TINTS)


@dataclass(frozen=True, slots=True)
class AlgorithmStage:
    """One step of an algorithm, as it appears in a flow sheet.

    Parameters
    ----------
    label:
        The step's name, in repository terminology.
    detail:
        Lines under the label saying what the step computes. Keep each short;
        the figure is a map, not the prose.
    role:
        One of ``input``, ``compute``, ``decision``, ``reject``, ``output``.
        Chooses the card tint. ``decision`` marks a step that can reject a
        candidate; ``reject`` marks where rejected candidates go.
    formula:
        Optional single-line mathematical statement shown in the card, for the
        step whose mathematics is the point of the figure.
    """

    label: str
    detail: Sequence[str] = ()
    role: str = "compute"
    formula: str | None = None

    def __post_init__(self) -> None:
        if self.role not in _ROLES:
            raise ValueError(
                f"AlgorithmStage.role must be one of {', '.join(_ROLES)}; got '{self.role}'."
            )
        if not str(self.label).strip():
            raise ValueError("AlgorithmStage.label must be non-empty.")
        object.__setattr__(self, "detail", tuple(self.detail))


@dataclass(frozen=True, slots=True)
class SideNote:
    """A constraint, tolerance or failure mode attached to one stage.

    Placed beside the stage it governs rather than in a legend, because a
    constraint read three inches from the step it constrains is a constraint the
    reader has to re-associate.

    ``stage_index`` is the index of the stage within its lane's ``stages`` list,
    counted across the whole figure in declaration order.
    """

    stage_index: int
    title: str
    lines: Sequence[str] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "lines", tuple(self.lines))
        if self.stage_index < 0:
            raise ValueError("SideNote.stage_index must be non-negative.")


@dataclass(frozen=True, slots=True)
class _Lane:
    """One phase of the algorithm: a caption plus the stages it contains."""

    caption: str
    stages: tuple[AlgorithmStage, ...] = field(default_factory=tuple)


def _stage_card(
    stage: AlgorithmStage,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
) -> str:
    """Render one stage card with its label, formula and detail lines."""

    fill, stroke = ROLE_TINTS[stage.role]
    ink = ROLE_INK[stage.role]
    parts = [card(x, y, width, height, fill=fill, stroke=stroke, stroke_width=1.4)]
    cursor = y + 27.0
    for line in wrap_text(stage.label, 15.0, width - 28.0):
        parts.append(text(x + 14.0, cursor, line, size=15.0, fill=ink, weight="bold"))
        cursor += 19.0
    if stage.formula is not None:
        cursor += 3.0
        # Wrapped like every other string: a formula that overflows its card is
        # the defect the figure layout guard exists to catch, and measuring here
        # is cheaper than repairing figures one at a time.
        for line in wrap_text(stage.formula, 13.0, width - 28.0):
            parts.append(text(x + 14.0, cursor, line, size=13.0, fill=ACCENT))
            cursor += 19.0
    for line in stage.detail:
        for wrapped in wrap_text(line, 12.0, width - 28.0):
            parts.append(text(x + 14.0, cursor, wrapped, size=12.0, fill=MUTED))
            cursor += 16.0
    return "\n".join(parts)


def _stage_height(stage: AlgorithmStage, width: float) -> float:
    """The height a stage card needs for its content, before lane alignment."""

    lines = len(wrap_text(stage.label, 15.0, width - 28.0))
    height = 27.0 + 19.0 * lines
    if stage.formula is not None:
        height += 3.0 + 19.0 * len(wrap_text(stage.formula, 13.0, width - 28.0))
    for line in stage.detail:
        height += 16.0 * len(wrap_text(line, 12.0, width - 28.0))
    return height + 14.0


def algorithm_flow_svg(
    lanes: Sequence[tuple[str, Sequence[AlgorithmStage]]],
    *,
    title: str,
    subtitle: str,
    description: str,
    notes: Sequence[SideNote] = (),
    footer: Sequence[str] = (),
    stage_width: float = 236.0,
    lane_gap: float = 46.0,
    stage_gap: float = 34.0,
    note_width: float = 300.0,
) -> str:
    """Render an algorithm as a lane-based flow sheet.

    What it does
        Lays each phase out as a horizontal lane of stage cards, connects
        consecutive stages within a lane and consecutive lanes with arrows, and
        attaches constraint notes to the right of the lane that owns them.
        Card heights are computed from their content and equalized within a
        lane, so the figure reads as rows rather than as ragged boxes.

    When to use it
        For any documented algorithm whose structure is worth seeing before it
        is read — which is most of them. Geometry figures (what a vector or an
        angle *is*) belong in a dedicated renderer instead.

    Parameters
    ----------
    lanes:
        ``(caption, stages)`` pairs, one per phase, in execution order.
    title, subtitle, description:
        Figure heading, the line under it, and the ``<desc>`` text a screen
        reader announces. All three are mandatory: a figure that does not say
        what it shows is not a documentation asset.
    notes:
        Constraints and failure modes, each attached to a stage by its index in
        declaration order across the whole figure.
    footer:
        Closing lines, for the honest statement of what the figure omits.
    stage_width:
        Width of every stage card. Fixed rather than per-card so lanes align.

    Returns
    -------
    str
        A complete SVG document with ``<title>`` and ``<desc>``.
    """

    resolved = tuple(_Lane(caption=caption, stages=tuple(stages)) for caption, stages in lanes)
    if not resolved:
        raise ValueError("algorithm_flow_svg requires at least one lane.")
    if any(not lane.stages for lane in resolved):
        raise ValueError("Every lane must contain at least one stage.")

    total_stages = sum(len(lane.stages) for lane in resolved)
    for note in notes:
        if note.stage_index >= total_stages:
            raise ValueError(
                f"SideNote.stage_index {note.stage_index} exceeds the {total_stages} stages."
            )

    widest_lane = max(len(lane.stages) for lane in resolved)
    margin = 40.0
    flow_width = widest_lane * stage_width + (widest_lane - 1) * stage_gap
    has_notes = bool(notes)
    body_width = flow_width + (note_width + 40.0 if has_notes else 0.0)
    width = max(margin * 2.0 + body_width, header_width(title, subtitle))

    # Group notes by the lane that owns their stage, so a lane's notes stack
    # beside it rather than floating at the figure's edge.
    lane_of_stage: list[int] = []
    for lane_index, lane in enumerate(resolved):
        lane_of_stage.extend([lane_index] * len(lane.stages))
    notes_by_lane: dict[int, list[SideNote]] = {}
    for note in notes:
        notes_by_lane.setdefault(lane_of_stage[note.stage_index], []).append(note)

    parts: list[str] = []
    marker = "pytex-algo-arrow"
    cursor_y = 108.0

    for lane_index, lane in enumerate(resolved):
        lane_stage_height = max(_stage_height(stage, stage_width) for stage in lane.stages)
        note_block_height = sum(
            34.0 + 19.0 * len(note.lines) + 12.0 + 14.0
            for note in notes_by_lane.get(lane_index, [])
        )
        lane_height = max(lane_stage_height, note_block_height)

        # Lane captions take the documentation/teaching token, so a caption is
        # never mistaken for a stage label.
        parts.append(text(margin, cursor_y - 12.0, lane.caption, size=13.0, fill=VIOLET))

        for position, stage in enumerate(lane.stages):
            x = margin + position * (stage_width + stage_gap)
            parts.append(
                _stage_card(stage, x=x, y=cursor_y, width=stage_width, height=lane_stage_height)
            )
            if position > 0:
                previous_x = margin + (position - 1) * (stage_width + stage_gap)
                parts.append(
                    f'  <line x1="{previous_x + stage_width:.1f}" '
                    f'y1="{cursor_y + lane_stage_height / 2.0:.1f}" '
                    f'x2="{x - 4.0:.1f}" y2="{cursor_y + lane_stage_height / 2.0:.1f}" '
                    f'stroke="{MUTED}" stroke-width="2" marker-end="url(#{marker})"/>'
                )

        note_y = cursor_y
        for note in notes_by_lane.get(lane_index, []):
            fragment, _, note_height = callout(
                margin + flow_width + 40.0,
                note_y,
                note.title,
                note.lines,
                width=note_width,
            )
            parts.append(fragment)
            note_y += note_height + 14.0

        if lane_index + 1 < len(resolved):
            arrow_x = margin + stage_width / 2.0
            parts.append(
                f'  <line x1="{arrow_x:.1f}" y1="{cursor_y + lane_height:.1f}" '
                f'x2="{arrow_x:.1f}" y2="{cursor_y + lane_height + lane_gap - 6.0:.1f}" '
                f'stroke="{ACCENT}" stroke-width="2.4" marker-end="url(#{marker}-flow)"/>'
            )
        cursor_y += lane_height + lane_gap

    footer_y = cursor_y + 6.0
    for line in footer:
        parts.append(text(margin, footer_y, line, size=12.0, fill=MUTED))
        footer_y += 18.0

    height = footer_y + 22.0
    markers = "\n".join([arrow_marker(marker, MUTED), arrow_marker(f"{marker}-flow", ACCENT)])
    return document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(parts),
        marker_defs=markers,
    )
