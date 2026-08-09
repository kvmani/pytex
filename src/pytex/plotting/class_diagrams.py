"""UML-style class and object-model diagrams, rendered as canonical SVG.

These are the figures behind the Class & Object Model Atlas: an architecture
overview, the class hierarchy, and one object-model diagram per scientific
domain. They are drawn on the same stack as the reference-frame and algorithm
figures (`pytex.plotting.svg_primitives`, the style-guide tokens, and the
Helvetica advance-width metrics in `pytex.plotting._svg_text`), so the atlas
looks like the rest of the documentation rather than like a tool's default
output.

Why a renderer rather than Graphviz
-----------------------------------
The usual route — ``pyreverse`` or ``sphinx.ext.inheritance_diagram`` — renders
through the Graphviz ``dot`` binary. That is a system dependency the docs build
does not otherwise carry, and its default styling ignores the visualization
style guide entirely. The layout done here is a compact Sugiyama pipeline:
assign layers by longest path over the acyclic part of the graph, order each
layer by repeated barycentre sweeps to reduce crossings, then place nodes with a
median-and-separate pass. It is deterministic — these figures are committed
assets compared byte-for-byte by the test suite.

What the notation means
-----------------------
A **filled diamond** marks composition: a required field, so the owner cannot
exist without the part. A **dashed line with a plain arrowhead** marks
association: an optional or defaulted field. A **hollow triangle** at the upper
end marks inheritance, in the usual UML direction (the triangle touches the
base class). Multiplicity ``*`` on an edge means the field holds a sequence.

See also
--------
`scripts/class_model.py` : reads the model these diagrams draw.
`scripts/generate_class_model_figures.py` : writes the committed assets.
`docs/standards/visualization_style_guide.md` : the canonical tokens.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from pytex.plotting._svg_text import text_width
from pytex.plotting.svg_primitives import (
    ACCENT,
    AMBER,
    GREEN,
    INK,
    MUTED,
    PANEL,
    ROSE,
    TEAL,
    VIOLET,
    arrow_marker,
    card,
    document,
    text,
)

__all__ = [
    "ROLE_COLOURS",
    "ClassBox",
    "ClassEdge",
    "class_diagram_svg",
]

#: Accent colour per semantic role. Roles are the PyTex subpackages plus the
#: kinds a reader must be able to tell apart at a glance. Colour groups; the
#: card's module line carries the same information in text, so nothing is lost
#: to a reader who cannot distinguish the hues.
ROLE_COLOURS: dict[str, str] = {
    "core": ACCENT,
    "texture": VIOLET,
    "ebsd": TEAL,
    "diffraction": AMBER,
    "tem": "#0ea5e9",
    "plotting": "#db2777",
    "properties": GREEN,
    "experimental": ROSE,
    "adapters": "#b45309",
    "external": MUTED,
    "governance": "#475569",
}

#: Very light fill derived from each accent, so a card reads as tinted paper
#: rather than as a saturated block.
_ROLE_FILLS: dict[str, str] = {
    "core": "#eef4ff",
    "texture": "#f5efff",
    "ebsd": "#e9f7f7",
    "diffraction": "#fdf6e8",
    "tem": "#eaf6fe",
    "plotting": "#fdeef6",
    "properties": "#eef7f0",
    "experimental": "#fdeef1",
    "adapters": "#fbf2e6",
    "external": "#f4f6fa",
    "governance": "#f1f5f9",
}

# Type geometry. Every dimension here is measured, never guessed: card widths
# come from `text_width`, so a long field name widens its card instead of
# spilling out of it.
_TITLE_SIZE = 18.0
_MODULE_SIZE = 13.0
_FIELD_SIZE = 15.0
_STEREOTYPE_SIZE = 13.0
_PAD_X = 15.0
_LINE_HEIGHT = 20.0
_HEADER_HEIGHT = 52.0
_MIN_CARD_WIDTH = 160.0
_MAX_CARD_WIDTH = 340.0
_NODE_GAP = 34.0
_LAYER_GAP = 100.0
_MARGIN = 40.0
_BODY_TOP = 110.0
_EDGE_LABEL_SIZE = 13.0


@dataclass(frozen=True)
class ClassBox:
    """One class card in a diagram."""

    key: str
    """Identity used by edges; usually the fully qualified class name."""

    name: str
    module: str = ""
    """Shown under the name, e.g. ``core.lattice``. Empty hides the line."""

    stereotype: str = ""
    """UML stereotype such as ``frozen dataclass`` or ``StrEnum``."""

    fields: tuple[str, ...] = ()
    """Attribute lines, already formatted as ``name: Type``."""

    role: str = "core"
    footnote: str = ""
    """One dimmed line under the fields, e.g. ``+ 6 more``."""

    emphasis: bool = False
    """Draw with a heavier border: the anchor object of the diagram."""


@dataclass(frozen=True)
class ClassEdge:
    """A relation drawn between two cards."""

    source: str
    target: str
    kind: str = "composition"
    """``composition``, ``association``, ``inheritance`` or ``flow``."""

    label: str = ""
    multiplicity: str = ""


@dataclass
class _Placed:
    box: ClassBox
    width: float
    height: float
    layer: int
    order: int
    x: float = 0.0
    y: float = 0.0
    lines: tuple[str, ...] = ()

    @property
    def centre(self) -> float:
        return self.x + self.width / 2.0

    @property
    def bottom(self) -> float:
        return self.y + self.height


# ---------------------------------------------------------------------------
# Card measurement
# ---------------------------------------------------------------------------


def _truncate(content: str, size: float, limit: float) -> str:
    """Shorten a line until it measures under ``limit``, marking the cut."""

    if text_width(content, size) <= limit:
        return content
    trimmed = content
    while trimmed and text_width(f"{trimmed}…", size) > limit:
        trimmed = trimmed[:-1]
    return f"{trimmed}…"


def _measure(box: ClassBox) -> _Placed:
    """Size a card to its own text, then clip any line that still overruns."""

    # The header is a floor, not a candidate: a class name is not something a
    # reader can reconstruct from a truncation, so a long name widens its card
    # past the field cap rather than being cut. Field lines, whose type suffix is
    # recoverable from the API docs, are clipped to keep cards comparable.
    header_need = text_width(box.name, _TITLE_SIZE, bold=True) + 2.0 * _PAD_X
    inner_limit = _MAX_CARD_WIDTH - 2.0 * _PAD_X
    lines = tuple(_truncate(line, _FIELD_SIZE, inner_limit) for line in box.fields)

    candidates = [text_width(line, _FIELD_SIZE) for line in lines]
    if box.module:
        candidates.append(text_width(box.module, _MODULE_SIZE))
    if box.stereotype:
        candidates.append(text_width(f"«{box.stereotype}»", _STEREOTYPE_SIZE))
    if box.footnote:
        candidates.append(text_width(box.footnote, _FIELD_SIZE))

    width = (max(candidates) if candidates else 0.0) + 2.0 * _PAD_X
    width = min(max(width, _MIN_CARD_WIDTH), _MAX_CARD_WIDTH)
    width = max(width, header_need)

    rows = len(lines) + (1 if box.footnote else 0)
    height = _HEADER_HEIGHT + (14.0 + _LINE_HEIGHT * rows if rows else 0.0)
    if box.stereotype:
        height += 17.0
    return _Placed(box=box, width=width, height=height, layer=0, order=0, lines=lines)


# ---------------------------------------------------------------------------
# Layered layout
# ---------------------------------------------------------------------------


def _acyclic(keys: Sequence[str], edges: Sequence[ClassEdge]) -> list[ClassEdge]:
    """Drop the edges that close a cycle, keeping the rest for layering.

    Mutual references are real in a scientific model (a report holding its
    configuration, a configuration naming its report type). They are still drawn
    — only their contribution to the *layering* is dropped, since a layered
    drawing has no consistent level for a cycle.
    """

    outgoing: dict[str, list[ClassEdge]] = {key: [] for key in keys}
    for edge in edges:
        if edge.source in outgoing and edge.target in outgoing:
            outgoing[edge.source].append(edge)

    state: dict[str, int] = {key: 0 for key in keys}
    kept: list[ClassEdge] = []

    def visit(key: str) -> None:
        state[key] = 1
        for edge in outgoing[key]:
            if state[edge.target] == 1:
                continue  # back edge
            kept.append(edge)
            if state[edge.target] == 0:
                visit(edge.target)
        state[key] = 2

    for key in keys:
        if state[key] == 0:
            visit(key)
    return kept


def _assign_layers(keys: Sequence[str], edges: Sequence[ClassEdge]) -> dict[str, int]:
    """Longest-path layering: a node sits one level below its deepest owner."""

    forward = _acyclic(keys, edges)
    incoming: dict[str, list[str]] = {key: [] for key in keys}
    for edge in forward:
        incoming[edge.target].append(edge.source)

    layer: dict[str, int] = {}

    def depth(key: str, guard: frozenset[str] = frozenset()) -> int:
        if key in layer:
            return layer[key]
        if key in guard:  # pragma: no cover - cycles already removed
            return 0
        parents = incoming[key]
        value = 0 if not parents else 1 + max(depth(p, guard | {key}) for p in parents)
        layer[key] = value
        return value

    for key in keys:
        depth(key)
    return layer


def _order_layers(
    layers: Mapping[int, list[_Placed]],
    edges: Sequence[ClassEdge],
    nodes: Mapping[str, _Placed],
    *,
    sweeps: int = 6,
) -> None:
    """Barycentre sweeps: repeatedly sort each layer by its neighbours' positions.

    Two passes in each direction are enough for graphs of this size, and the
    result is deterministic because ties break on the existing order.
    """

    neighbours_up: dict[str, list[str]] = {key: [] for key in nodes}
    neighbours_down: dict[str, list[str]] = {key: [] for key in nodes}
    for edge in edges:
        if edge.source not in nodes or edge.target not in nodes:
            continue
        neighbours_up[edge.target].append(edge.source)
        neighbours_down[edge.source].append(edge.target)

    indices = sorted(layers)
    for sweep in range(sweeps):
        downward = sweep % 2 == 0
        sequence = indices if downward else list(reversed(indices))
        for index in sequence:
            related = neighbours_up if downward else neighbours_down
            row = layers[index]
            positions = {node.box.key: position for position, node in enumerate(row)}

            def barycentre(
                node: _Placed,
                _related: Mapping[str, list[str]] = related,
                _positions: Mapping[str, int] = positions,
            ) -> float:
                ranks = [
                    nodes[other].order
                    for other in _related[node.box.key]
                    if other in nodes and nodes[other].layer != node.layer
                ]
                if not ranks:
                    return float(_positions[node.box.key])
                return sum(ranks) / len(ranks)

            row.sort(key=lambda node: (barycentre(node), node.box.key))
            for position, node in enumerate(row):
                node.order = position


def _place(
    layers: Mapping[int, list[_Placed]],
    edges: Sequence[ClassEdge],
    nodes: Mapping[str, _Placed],
    *,
    passes: int = 4,
) -> None:
    """Assign coordinates: pack each layer, then pull nodes toward their kin.

    The refinement is the classic median-and-separate heuristic — move each node
    toward the average centre of the nodes it connects to, then restore the
    minimum gap left-to-right and right-to-left so the ordering survives.
    """

    indices = sorted(layers)
    y = _BODY_TOP
    for index in indices:
        row = layers[index]
        cursor = 0.0
        for node in row:
            node.x = cursor
            node.y = y
            cursor += node.width + _NODE_GAP
        y += max(node.height for node in row) + _LAYER_GAP

    widest = max(
        sum(node.width for node in row) + _NODE_GAP * (len(row) - 1) for row in layers.values()
    )
    for row in layers.values():
        span = sum(node.width for node in row) + _NODE_GAP * (len(row) - 1)
        shift = (widest - span) / 2.0
        for node in row:
            node.x += shift

    linked: dict[str, list[str]] = {key: [] for key in nodes}
    for edge in edges:
        if edge.source in nodes and edge.target in nodes:
            linked[edge.source].append(edge.target)
            linked[edge.target].append(edge.source)

    for step in range(passes):
        sequence = indices if step % 2 == 0 else list(reversed(indices))
        for index in sequence:
            row = layers[index]
            for node in row:
                partners = [
                    nodes[other]
                    for other in linked[node.box.key]
                    if nodes[other].layer != node.layer
                ]
                if partners:
                    target = sum(p.centre for p in partners) / len(partners)
                    node.x = target - node.width / 2.0
            row.sort(key=lambda node: (node.x, node.box.key))
            for position, node in enumerate(row):
                node.order = position
            for position in range(1, len(row)):
                minimum = row[position - 1].x + row[position - 1].width + _NODE_GAP
                row[position].x = max(row[position].x, minimum)
            for position in range(len(row) - 2, -1, -1):
                maximum = row[position + 1].x - _NODE_GAP - row[position].width
                row[position].x = min(row[position].x, maximum)

    _pack_components(layers, edges, nodes)

    left = min(node.x for row in layers.values() for node in row)
    for row in layers.values():
        for node in row:
            node.x += _MARGIN - left


def _pack_components(
    layers: Mapping[int, list[_Placed]],
    edges: Sequence[ClassEdge],
    nodes: Mapping[str, _Placed],
) -> None:
    """Slide disconnected groups together so they do not leave a canyon between them.

    Each group is pulled under its own parents, which on a figure with several
    unrelated families (the class hierarchy, where the tensors and the tilt
    envelopes share nothing) leaves a wide empty strip in the middle. Packing
    moves whole groups, never their internals, so no relation is redrawn.
    """

    parent: dict[str, str] = {key: key for key in nodes}

    def find(key: str) -> str:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    for edge in edges:
        if edge.source in nodes and edge.target in nodes:
            parent[find(edge.source)] = find(edge.target)

    groups: dict[str, list[_Placed]] = {}
    for key, node in nodes.items():
        groups.setdefault(find(key), []).append(node)
    if len(groups) < 2:
        return

    ordered = sorted(groups.values(), key=lambda group: min(node.x for node in group))
    cursor = min(node.x for node in ordered[0])
    for group in ordered:
        start = min(node.x for node in group)
        shift = cursor - start
        for node in group:
            node.x += shift
        cursor = max(node.x + node.width for node in group) + _NODE_GAP * 2.0


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def _card_svg(node: _Placed) -> str:
    """One class card: name, module, stereotype, then its fields."""

    box = node.box
    accent = ROLE_COLOURS.get(box.role, ACCENT)
    fill = _ROLE_FILLS.get(box.role, PANEL)
    parts = [
        card(
            node.x,
            node.y,
            node.width,
            node.height,
            fill=fill,
            stroke=accent,
            stroke_width=2.0 if box.emphasis else 1.2,
        ),
        # Accent bar: a 4-unit rule under the header, tying the card to its role.
        f'  <rect x="{node.x + 1.0:.1f}" y="{node.y + 1.0:.1f}" width="{node.width - 2.0:.1f}" '
        f'height="4" rx="2" fill="{accent}" opacity="0.85"/>',
        text(node.x + _PAD_X, node.y + 29.0, box.name, size=_TITLE_SIZE, fill=INK, weight="bold"),
    ]
    cursor = node.y + 46.0
    if box.module:
        parts.append(text(node.x + _PAD_X, cursor, box.module, size=_MODULE_SIZE, fill=MUTED))
        cursor += 18.0
    if box.stereotype:
        parts.append(
            text(
                node.x + _PAD_X,
                cursor,
                f"«{box.stereotype}»",
                size=_STEREOTYPE_SIZE,
                fill=accent,
            )
        )
        cursor += 17.0
    if node.lines or box.footnote:
        divider = cursor - 13.0
        parts.append(
            f'  <line x1="{node.x + 1.0:.1f}" y1="{divider:.1f}" '
            f'x2="{node.x + node.width - 1.0:.1f}" y2="{divider:.1f}" '
            f'stroke="{accent}" stroke-width="1" opacity="0.35"/>'
        )
        cursor += 9.0
        for line in node.lines:
            parts.append(text(node.x + _PAD_X, cursor, line, size=_FIELD_SIZE, fill=INK))
            cursor += _LINE_HEIGHT
        if box.footnote:
            parts.append(text(node.x + _PAD_X, cursor, box.footnote, size=_FIELD_SIZE, fill=MUTED))
    return "\n".join(parts)


def _ports(
    nodes: Mapping[str, _Placed], edges: Sequence[ClassEdge]
) -> dict[int, tuple[float, float]]:
    """Spread each node's connectors across its edge instead of stacking them.

    Every edge leaving one point makes a fan that hides which relation is which.
    Connectors are distributed across the middle 70% of the card edge, ordered by
    where the other end sits, so lines do not cross needlessly near the card.
    """

    out_ports: dict[int, tuple[float, float]] = {}
    grouped_out: dict[str, list[int]] = {}
    grouped_in: dict[str, list[int]] = {}
    for index, edge in enumerate(edges):
        grouped_out.setdefault(edge.source, []).append(index)
        grouped_in.setdefault(edge.target, []).append(index)

    def spread(node: _Placed, count: int, position: int) -> float:
        usable = node.width * 0.7
        if count == 1:
            return node.centre
        step = usable / (count - 1)
        return node.centre - usable / 2.0 + step * position

    for key, indices in grouped_out.items():
        node = nodes[key]
        ordered = sorted(indices, key=lambda i: (nodes[edges[i].target].centre, edges[i].label))
        for position, index in enumerate(ordered):
            out_ports[index] = (spread(node, len(ordered), position), 0.0)
    for key, indices in grouped_in.items():
        node = nodes[key]
        ordered = sorted(indices, key=lambda i: (nodes[edges[i].source].centre, edges[i].label))
        for position, index in enumerate(ordered):
            x, _ = out_ports[index]
            out_ports[index] = (x, spread(node, len(ordered), position))
    return out_ports


_EDGE_STYLE: dict[str, tuple[str, bool, str]] = {
    # kind: (colour, dashed, marker-start decoration)
    "composition": (INK, False, "diamond"),
    "association": (MUTED, True, ""),
    "inheritance": ("#334155", False, "triangle"),
    "flow": (ACCENT, False, ""),
}


def _edge_svg(
    edge: ClassEdge,
    source: _Placed,
    target: _Placed,
    ports: tuple[float, float],
) -> str:
    """One relation, routed as a gentle curve between two card edges."""

    colour, dashed, decoration = _EDGE_STYLE.get(edge.kind, _EDGE_STYLE["composition"])
    dash = ' stroke-dasharray="6 5"' if dashed else ""
    start_x, end_x = ports
    marker_start = f' marker-start="url(#cd-{decoration})"' if decoration else ""

    # A relation spanning several layers has to travel past the cards between its
    # ends. Drawing those long runs lighter than the adjacent-layer ones splits
    # the picture into a readable foreground of local structure and a background
    # of long-range references, instead of one uniform thicket.
    span = abs(target.layer - source.layer)
    opacity = 1.0 if span <= 1 else (0.55 if span == 2 else 0.4)
    weight = 1.6 if span <= 1 else 1.3
    fade = f' opacity="{opacity:.2f}"' if opacity < 1.0 else ""

    if target.layer > source.layer:
        y1, y2 = source.bottom, target.y
        lift = max(24.0, (y2 - y1) * 0.42)
        path = (
            f"M{start_x:.1f},{y1:.1f} C{start_x:.1f},{y1 + lift:.1f} "
            f"{end_x:.1f},{y2 - lift:.1f} {end_x:.1f},{y2:.1f}"
        )
    else:
        # Same-layer or upward relation: route around the right-hand side so it
        # never runs underneath the cards it connects.
        y1 = source.y + source.height / 2.0
        y2 = target.y + target.height / 2.0
        x1 = source.x + source.width
        x2 = target.x + target.width
        bulge = 46.0
        path = (
            f"M{x1:.1f},{y1:.1f} C{x1 + bulge:.1f},{y1:.1f} "
            f"{x2 + bulge:.1f},{y2:.1f} {x2:.1f},{y2:.1f}"
        )

    return (
        f'  <path d="{path}" fill="none" stroke="{colour}" stroke-width="{weight}"{dash}'
        f'{fade} marker-end="url(#cd-arrow)"{marker_start}/>'
    )


def _edge_label(edge: ClassEdge) -> str:
    """The text drawn beside a relation, or an empty string for none."""

    if edge.kind == "inheritance":
        return ""
    if edge.multiplicity == "*":
        return f"{edge.label} *" if edge.label else "*"
    return edge.label


def _label_svg(
    label: str,
    x: float,
    y: float,
    occupied: list[tuple[float, float, float]],
) -> str:
    """Place a relation label in the gap between layers, or skip it if it collides.

    Skipping is deliberate. Every relation's field name is already printed inside
    its owning card, so a label that cannot be placed cleanly costs the reader
    nothing, whereas two labels overprinting cost them the diagram.
    """

    half = text_width(label, _EDGE_LABEL_SIZE) / 2.0 + 3.0
    for row in range(3):
        candidate_y = y + row * 18.0
        clash = any(
            abs(candidate_y - oy) < 15.0 and abs(x - ox) < half + oh for ox, oy, oh in occupied
        )
        if clash:
            continue
        occupied.append((x, candidate_y, half))
        return (
            f'  <rect x="{x - half:.1f}" y="{candidate_y - 11.5:.1f}" width="{2 * half:.1f}" '
            f'height="13" rx="3" fill="#ffffff" opacity="0.86"/>\n'
            + text(x, candidate_y, label, size=_EDGE_LABEL_SIZE, fill=MUTED, anchor="middle")
        )
    return ""


def _markers() -> str:
    """Arrowhead, hollow inheritance triangle, and filled composition diamond."""

    return "\n".join(
        [
            arrow_marker("cd-arrow", MUTED),
            (
                '    <marker id="cd-triangle" markerUnits="userSpaceOnUse" markerWidth="14" '
                'markerHeight="12" refX="1" refY="6" orient="auto">\n'
                '      <path d="M13,0 L1,6 L13,12 z" fill="#ffffff" stroke="#334155" '
                'stroke-width="1.4"/>\n'
                "    </marker>"
            ),
            (
                '    <marker id="cd-diamond" markerUnits="userSpaceOnUse" markerWidth="15" '
                'markerHeight="10" refX="1" refY="5" orient="auto">\n'
                f'      <path d="M1,5 L8,1 L15,5 L8,9 z" fill="{INK}"/>\n'
                "    </marker>"
            ),
        ]
    )


_LEGEND_ENTRIES: tuple[tuple[str, str], ...] = (
    ("composition", "composition - a required field: the owner cannot exist without it"),
    ("association", "association - an optional or defaulted field"),
    ("inheritance", "inheritance - the triangle touches the base class"),
    ("flow", "typed references crossing a package boundary"),
)


def _legend_svg(
    x: float,
    y: float,
    kinds: Sequence[str],
    notes: Sequence[str],
    *,
    starred: bool,
) -> tuple[str, float]:
    """A legend explaining only the notation the figure actually uses."""

    rows = [(kind, label) for kind, label in _LEGEND_ENTRIES if kind in kinds]
    tail = ["*  after a field name: the field holds a sequence"] if starred else []
    width = (
        max(
            [text_width("Notation", 15.0, bold=True)]
            + [text_width(label, 13.5) + 66.0 for _, label in rows]
            + [text_width(row, 13.5) + 28.0 for row in tail]
            + [text_width(note, 13.0) + 20.0 for note in notes]
        )
        + 26.0
    )
    line_count = len(rows) + len(tail)
    height = 40.0 + 21.0 * line_count + (10.0 + 18.0 * len(notes) if notes else 0.0) + 12.0
    parts = [
        card(x, y, width, height, fill="#f7faff", stroke="#c8d8f2"),
        text(x + 16.0, y + 26.0, "Notation", size=15.0, fill=INK, weight="bold"),
    ]
    cursor = y + 52.0
    for kind, label in rows:
        colour, dashed, decoration = _EDGE_STYLE[kind]
        dash = ' stroke-dasharray="6 5"' if dashed else ""
        marker_start = f' marker-start="url(#cd-{decoration})"' if decoration else ""
        parts.append(
            f'  <line x1="{x + 18.0:.1f}" y1="{cursor - 4.0:.1f}" x2="{x + 58.0:.1f}" '
            f'y2="{cursor - 4.0:.1f}" stroke="{colour}" stroke-width="1.6"{dash} '
            f'marker-end="url(#cd-arrow)"{marker_start}/>'
        )
        parts.append(text(x + 76.0, cursor, label, size=13.5, fill=MUTED))
        cursor += 21.0
    for row in tail:
        parts.append(text(x + 18.0, cursor, row, size=13.5, fill=MUTED))
        cursor += 21.0
    if notes:
        cursor += 6.0
        for note in notes:
            parts.append(text(x + 16.0, cursor, note, size=13.0, fill=MUTED))
            cursor += 18.0
    return "\n".join(parts), height


def class_diagram_svg(
    *,
    title: str,
    subtitle: str,
    description: str,
    boxes: Sequence[ClassBox],
    edges: Sequence[ClassEdge],
    notes: Sequence[str] = (),
    legend: bool = True,
) -> str:
    """Lay out and render one class or object-model diagram.

    What it does
        Places ``boxes`` in layers derived from ``edges`` — owners above the
        objects they own, base classes above their subclasses — orders each layer
        to reduce crossings, and draws the relations in UML notation.

    When to use it
        For the generated atlas figures. The layout is deterministic, so the
        result can be committed and compared byte-for-byte.

    Parameters
    ----------
    title, subtitle, description:
        Figure heading, one-line context, and the ``<desc>`` text a screen reader
        announces.
    boxes:
        The class cards. Order affects only tie-breaking.
    edges:
        Relations between card keys. Edges naming an absent card are ignored, so
        a caller may pass a whole domain's relations and a filtered card list.
    notes:
        Extra lines for the legend, e.g. what the view deliberately omits.
    legend:
        Draw the notation key. Turn it off for a figure whose edges are uniform.

    Returns
    -------
    str
        A complete SVG document, with the ``<title>`` and ``<desc>`` the
        visualization style guide requires.
    """

    if not boxes:
        raise ValueError("a class diagram needs at least one box")

    nodes = {box.key: _measure(box) for box in boxes}
    keys = list(nodes)
    drawn = [edge for edge in edges if edge.source in nodes and edge.target in nodes]

    layer_of = _assign_layers(keys, drawn)
    layers: dict[int, list[_Placed]] = {}
    for key, node in nodes.items():
        node.layer = layer_of[key]
        layers.setdefault(node.layer, []).append(node)
    for index in layers:
        layers[index].sort(key=lambda node: node.box.key)
        for position, node in enumerate(layers[index]):
            node.order = position

    _order_layers(layers, drawn, nodes)
    _place(layers, drawn, nodes)

    ports = _ports(nodes, drawn)
    body: list[str] = []
    labels: list[str] = []
    occupied: list[tuple[float, float, float]] = []

    # A label may only sit in the empty band *between* two layers. The bottom of
    # its own source card is not that band: cards in one layer have different
    # heights, so the space under a short card is occupied by its taller
    # neighbours, and a label placed there lands on top of one of them.
    band_top = {index: min(node.y for node in row) for index, row in layers.items()}
    band_bottom = {index: max(node.bottom for node in row) for index, row in layers.items()}

    for index, edge in enumerate(drawn):
        source, target = nodes[edge.source], nodes[edge.target]
        body.append(_edge_svg(edge, source, target, ports[index]))
        label = _edge_label(edge)
        if label and target.layer > source.layer:
            start_x = ports[index][0]
            gap_top = band_bottom[source.layer]
            gap_bottom = band_top.get(source.layer + 1, gap_top)
            # The label goes at the relation's own end, in the empty band just
            # below its owner. Labelling the midpoint instead loses every
            # long-range relation — and on the architecture view the numbers on
            # those long edges are the entire point of the figure.
            if gap_bottom - gap_top >= 58.0:
                labels.append(_label_svg(label, start_x, gap_top + 24.0, occupied))

    for node in nodes.values():
        body.append(_card_svg(node))
    body.extend(part for part in labels if part)

    content_right = max(node.x + node.width for node in nodes.values())
    content_bottom = max(node.bottom for node in nodes.values())

    if legend:
        kinds = [edge.kind for edge in drawn]
        legend_svg, legend_height = _legend_svg(
            _MARGIN,
            content_bottom + 34.0,
            kinds,
            notes,
            starred=any(edge.multiplicity == "*" for edge in drawn),
        )
        body.append(legend_svg)
        content_bottom += 34.0 + legend_height

    width = max(content_right + _MARGIN, _MARGIN * 2 + text_width(title, 25.0), 720.0)
    height = content_bottom + _MARGIN

    return document(
        width=width,
        height=height,
        title=title,
        description=description,
        subtitle=subtitle,
        body="\n".join(body),
        marker_defs=_markers(),
    )
