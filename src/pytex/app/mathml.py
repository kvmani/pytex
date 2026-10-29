# ruff: noqa: RUF001
"""Typeset equations for results: a small MathML tree, and its two printed forms.

Purpose
-------
Some results are only believable when the reader can follow the calculation:
the equation, the numbers substituted into it, the number that came out. Prose
cannot carry a 3x3 matrix, and a picture of an equation cannot be copied or
read by a screen reader. This module builds equations as a **MathML element
tree** in plain JSON, which

- the browser renders natively and offline (MathML Core; no script, no CDN),
- :func:`to_markup` writes into the printable HTML report, and
- :func:`to_latex` writes into the Markdown report as ``$$ ... $$``.

One tree, three surfaces, so the equation a reader judged on screen is the one
in the file.

Why a tree and not a TeX string
-------------------------------
A result crosses the service boundary as JSON and is read by a browser that
does not typeset TeX (``tests/unit/test_app_manifest.py`` rejects mathtext in a
result for that reason). A tree also makes the browser side safe by
construction: it builds elements from a whitelist of tags and attributes and
never parses markup, so no string in a result can become active content.

Where symbols come from
-----------------------
Registered quantities are written by name, :func:`symbol` ``("phi_1")``, and
converted from the mathtext form :mod:`pytex.core.symbols` holds; Miller indices
are converted from what :mod:`pytex.core.notation` writes. Neither is typed
here, which is the rule ``AGENTS.md`` sets for every surface.

Wire form
---------
A node is ``{"tag": ..., "text": ..., "font": ..., "attrs": {...}, "children": [...]}``
with only the keys it needs. ``font`` is one of :data:`FONTS` and is a property
of the glyphs, not of the tag: MathML Core renders only ``mathvariant="normal"``
reliably, so bold and sans-serif are applied by the renderer as CSS, and by
:func:`to_latex` as ``\\mathbf`` and ``\\mathsf``.
"""

from __future__ import annotations

import html
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

__all__ = [
    "ALLOWED_ATTRIBUTES",
    "ALLOWED_TAGS",
    "FONTS",
    "Node",
    "equation",
    "fenced",
    "frac",
    "from_mathtext",
    "matrix",
    "mi",
    "miller",
    "mn",
    "mo",
    "mtext",
    "number",
    "row",
    "space",
    "sub",
    "subsup",
    "sup",
    "symbol",
    "to_latex",
    "to_markup",
    "under",
    "validate",
    "vector",
]

Node = dict[str, Any]

#: The MathML Core elements a result may use. Anything else is refused when the
#: tree is built, and again by the browser renderer.
ALLOWED_TAGS = frozenset(
    {
        "math",
        "mrow",
        "mi",
        "mn",
        "mo",
        "mtext",
        "mspace",
        "msub",
        "msup",
        "msubsup",
        "mfrac",
        "msqrt",
        "mover",
        "munder",
        "mtable",
        "mtr",
        "mtd",
    }
)

#: Attributes a node may carry, all presentational.
ALLOWED_ATTRIBUTES = frozenset(
    {"display", "stretchy", "accent", "width", "lspace", "rspace", "form", "largeop"}
)

#: Glyph styles. ``None`` is MathML's default: a single-letter identifier in
#: italic, anything else upright.
FONTS = frozenset({"bold", "upright", "sans", "bold-italic"})

_LEAVES = frozenset({"mi", "mn", "mo", "mtext"})

#: The minus sign. A hyphen is not a minus: it is shorter, sits lower, and a
#: column of numbers using it does not align.
MINUS = "−"


# -- leaves ---------------------------------------------------------------------


def _leaf(tag: str, text: str, *, font: str | None = None, **attrs: str) -> Node:
    node: Node = {"tag": tag, "text": str(text)}
    if font is not None:
        if font not in FONTS:
            raise ValueError(f"Unknown font {font!r}; expected one of {sorted(FONTS)}.")
        node["font"] = font
    if attrs:
        node["attrs"] = dict(attrs)
    return node


def mi(text: str, *, font: str | None = None) -> Node:
    """An identifier: a variable, a function name, a label."""

    return _leaf("mi", text, font=font)


def mn(text: str) -> Node:
    """A number, as already formatted text."""

    return _leaf("mn", text)


def mo(text: str, **attrs: str) -> Node:
    """An operator, fence or separator."""

    return _leaf("mo", text, **attrs)


def mtext(text: str) -> Node:
    """Words inside an equation: upright, spaced as text."""

    return _leaf("mtext", text)


def space(width_em: float = 0.5) -> Node:
    """Horizontal space, in em."""

    return {"tag": "mspace", "attrs": {"width": f"{float(width_em):g}em"}}


def number(value: float, digits: int = 4) -> Node:
    """A real number to ``digits`` decimals, with a true minus sign.

    Negative zero is printed as zero: ``-0.0000`` in a matrix reads as a
    measurement of something negative.
    """

    if not math.isfinite(value):
        raise ValueError("Only finite numbers can be typeset.")
    text = f"{value:.{digits}f}"
    if text.startswith("-") and float(text) == 0.0:
        text = text[1:]
    return mn(text.replace("-", MINUS))


# -- structure ------------------------------------------------------------------


def _children(items: Iterable[Node | str]) -> list[Node]:
    nodes: list[Node] = []
    for item in items:
        nodes.append(mo(item) if isinstance(item, str) else item)
    return nodes


def row(*items: Node | str) -> Node:
    """A horizontal group. Bare strings are operators: ``row(a, "=", b)``."""

    return {"tag": "mrow", "children": _children(items)}


def sub(base: Node | str, script: Node | str) -> Node:
    """``base`` with a subscript."""

    return {"tag": "msub", "children": _children((base, script))}


def sup(base: Node | str, script: Node | str) -> Node:
    """``base`` with a superscript."""

    return {"tag": "msup", "children": _children((base, script))}


def subsup(base: Node | str, lower: Node | str, upper: Node | str) -> Node:
    """``base`` with both scripts."""

    return {"tag": "msubsup", "children": _children((base, lower, upper))}


def frac(numerator: Node | str, denominator: Node | str) -> Node:
    """A built-up fraction."""

    return {"tag": "mfrac", "children": _children((numerator, denominator))}


def under(base: Node | str, script: Node | str) -> Node:
    """A limit written under its operator: ``min`` over a set."""

    return {"tag": "munder", "children": _children((base, script))}


def fenced(content: Node, open_: str = "(", close: str = ")") -> Node:
    """``content`` between stretchy fences."""

    return row(
        mo(open_, stretchy="true", form="prefix"),
        content,
        mo(close, stretchy="true", form="postfix"),
    )


def matrix(values: Sequence[Sequence[float]] | Any, *, digits: int = 4) -> Node:
    """A numeric matrix in square brackets, every entry to ``digits`` decimals."""

    rows = [list(entry) for entry in values]
    if not rows or any(len(entry) != len(rows[0]) for entry in rows):
        raise ValueError("A matrix needs rows of one length.")
    table = {
        "tag": "mtable",
        "children": [
            {
                "tag": "mtr",
                "children": [
                    {"tag": "mtd", "children": [number(float(value), digits)]} for value in entry
                ],
            }
            for entry in rows
        ],
    }
    return fenced(table, "[", "]")


def vector(values: Sequence[float] | Any, *, digits: int = 4) -> Node:
    """A column vector."""

    return matrix([[float(value)] for value in values], digits=digits)


def equation(*items: Node | str) -> Node:
    """The root of one displayed equation."""

    node: Node = {"tag": "math", "attrs": {"display": "block"}, "children": [row(*items)]}
    validate(node)
    return node


# -- registered symbols and crystallographic notation ---------------------------

_GREEK: Mapping[str, str] = {
    "alpha": "α",
    "beta": "β",
    "gamma": "γ",
    "delta": "δ",
    "epsilon": "ϵ",
    "varepsilon": "ε",
    "eta": "η",
    "theta": "θ",
    "lambda": "λ",
    "mu": "μ",
    "nu": "ν",
    "xi": "ξ",
    "pi": "π",
    "rho": "ρ",
    "sigma": "σ",
    "tau": "τ",
    "phi": "ϕ",
    "varphi": "φ",
    "chi": "χ",
    "psi": "ψ",
    "omega": "ω",
    "Gamma": "Γ",
    "Delta": "Δ",
    "Theta": "Θ",
    "Lambda": "Λ",
    "Sigma": "Σ",
    "Phi": "Φ",
    "Psi": "Ψ",
    "Omega": "Ω",
}
_GREEK_TO_LATEX: Mapping[str, str] = {glyph: name for name, glyph in _GREEK.items()}

_FONT_COMMANDS: Mapping[str, str] = {
    "mathbf": "bold",
    "mathrm": "upright",
    "mathsf": "sans",
    "boldsymbol": "bold-italic",
}

_TOKEN = re.compile(r"\\[A-Za-z]+|\\[;,! {}]|[A-Za-z]|[0-9]+(?:\.[0-9]+)?|\S")


class _Parser:
    """Recursive descent over the small mathtext subset the registries emit."""

    def __init__(self, source: str) -> None:
        self.tokens = _TOKEN.findall(source)
        self.position = 0

    def peek(self) -> str | None:
        return self.tokens[self.position] if self.position < len(self.tokens) else None

    def take(self) -> str:
        token = self.peek()
        if token is None:
            raise ValueError("Unexpected end of mathtext.")
        self.position += 1
        return token

    def group(self, font: str | None) -> list[Node]:
        if self.peek() == "{":
            self.take()
            nodes = self.sequence(font, closing="}")
            self.take()
            return nodes
        return [self.atom(font)]

    def sequence(self, font: str | None, closing: str | None = None) -> list[Node]:
        nodes: list[Node] = []
        while (token := self.peek()) is not None and token != closing:
            if token in {"_", "^"}:
                if not nodes:
                    raise ValueError("A script needs a base.")
                self.take()
                script = _single(self.group(None))
                base = nodes.pop()
                if token == "^" and base.get("tag") == "msub":
                    nodes.append(subsup(base["children"][0], base["children"][1], script))
                elif token == "_" and base.get("tag") == "msup":
                    nodes.append(subsup(base["children"][0], script, base["children"][1]))
                else:
                    nodes.append((sub if token == "_" else sup)(base, script))
                continue
            nodes.append(self.atom(font))
        return nodes

    def atom(self, font: str | None) -> Node:
        token = self.take()
        if token == "{":
            nodes = self.sequence(font, closing="}")
            self.take()
            return _single(nodes)
        if token.startswith("\\"):
            name = token[1:]
            if name in _FONT_COMMANDS:
                return _single(self.group(_FONT_COMMANDS[name]))
            if name == "bar":
                return {
                    "tag": "mover",
                    "attrs": {"accent": "true"},
                    "children": [_single(self.group(font)), mo("¯", stretchy="true")],
                }
            if name in {";", ","}:
                return space(0.28 if name == ";" else 0.17)
            if name in {"{", "}"}:
                return mo(name, stretchy="false")
            if name == "langle":
                return mo("⟨", stretchy="false")
            if name == "rangle":
                return mo("⟩", stretchy="false")
            if name in _GREEK:
                glyph = _GREEK[name]
                upright = glyph.isupper() and font is None
                return mi(glyph, font="upright" if upright else font)
            raise ValueError(f"Unsupported mathtext command {token!r}.")
        if token[0].isdigit():
            return mn(token)
        if token.isalpha():
            return mi(token, font=font)
        if token == "-":
            return mo(MINUS)
        if token == "*":
            return mo("*")
        return mo(token, stretchy="false") if token in "()[]<>" else mo(token)


def _single(nodes: list[Node]) -> Node:
    return nodes[0] if len(nodes) == 1 else row(*nodes)


def from_mathtext(source: str) -> Node:
    """Convert a registry mathtext body (no ``$``) into a MathML node.

    Supports what :mod:`pytex.core.symbols` and :mod:`pytex.core.notation` emit:
    letters, numbers, Greek commands, ``\\mathbf``/``\\mathrm``/``\\mathsf``/
    ``\\boldsymbol``, ``\\bar``, thin spaces, braces and scripts. Anything else
    is refused rather than shown wrongly.
    """

    parser = _Parser(source.strip().strip("$"))
    nodes = parser.sequence(None)
    if parser.peek() is not None:
        raise ValueError(f"Unbalanced mathtext: {source!r}.")
    return _single(nodes)


def symbol(name: str) -> Node:
    """A registered symbol, by its registry name."""

    from pytex.core.symbols import symbol_latex

    return from_mathtext(symbol_latex(name))


def miller(indices: Sequence[int], *, family: str = "plane", scope: str = "specific") -> Node:
    """Miller or Miller-Bravais indices, overbarred, as :mod:`pytex.core.notation` writes them."""

    from pytex.core.notation import format_miller_indices

    text = format_miller_indices(
        tuple(int(value) for value in indices), family=family, style="mathtext", scope=scope
    )
    return from_mathtext(text)


# -- validation and the two printed forms ----------------------------------------


def validate(node: Mapping[str, Any]) -> None:
    """Refuse a tree with an unknown tag, attribute or font, or a malformed node."""

    tag = node.get("tag")
    if tag not in ALLOWED_TAGS:
        raise ValueError(f"MathML tag {tag!r} is not allowed in a result.")
    for key in node:
        if key not in {"tag", "text", "font", "attrs", "children"}:
            raise ValueError(f"Unknown MathML node key {key!r}.")
    for attribute in node.get("attrs", {}):
        if attribute not in ALLOWED_ATTRIBUTES:
            raise ValueError(f"MathML attribute {attribute!r} is not allowed in a result.")
    font = node.get("font")
    if font is not None and font not in FONTS:
        raise ValueError(f"Unknown font {font!r}.")
    if tag in _LEAVES:
        if not isinstance(node.get("text"), str) or "children" in node:
            raise ValueError(f"<{tag}> holds text and no children.")
        return
    if "text" in node:
        raise ValueError(f"<{tag}> holds children, not text.")
    for child in node.get("children", ()):
        validate(child)


_FONT_STYLE: Mapping[str, str] = {
    "bold": "font-weight:bold",
    "upright": "",
    "sans": "font-family:sans-serif",
    "bold-italic": "font-weight:bold;font-style:italic",
}


def to_markup(node: Mapping[str, Any]) -> str:
    """MathML markup for an HTML document, every text escaped."""

    tag = str(node["tag"])
    attrs = dict(node.get("attrs", {}))
    if tag == "math":
        attrs.setdefault("xmlns", "http://www.w3.org/1998/Math/MathML")
    font = node.get("font")
    if font is not None:
        attrs["mathvariant"] = "normal"
        style = _FONT_STYLE[str(font)]
        if style:
            attrs["style"] = style
    opening = tag + "".join(
        f' {key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items()
    )
    if tag in _LEAVES:
        return f"<{opening}>{html.escape(str(node['text']))}</{tag}>"
    inner = "".join(to_markup(child) for child in node.get("children", ()))
    return f"<{opening}>{inner}</{tag}>"


_OPERATOR_LATEX: Mapping[str, str] = {
    MINUS: "-",
    "×": r"\times",
    "⋅": r"\cdot",
    "∈": r"\in",
    "≤": r"\le",
    "≥": r"\ge",
    "≈": r"\approx",
    "∥": r"\parallel",
    "→": r"\to",
    "⟨": r"\langle",
    "⟩": r"\rangle",
    "{": r"\{",
    "}": r"\}",
    "‖": r"\|",
    "⁡": "",
    "⁢": "",
}

_FONT_LATEX: Mapping[str, str] = {
    "bold": r"\mathbf",
    "upright": r"\mathrm",
    "sans": r"\mathsf",
    "bold-italic": r"\boldsymbol",
}


def _latex_text(text: str) -> str:
    return "".join(
        "\\" + _GREEK_TO_LATEX[char] + " " if char in _GREEK_TO_LATEX else char for char in text
    ).strip()


def to_latex(node: Mapping[str, Any]) -> str:
    """LaTeX for the Markdown report (KaTeX/MathJax compatible)."""

    tag = node["tag"]
    children = list(node.get("children", ()))
    if tag == "mi":
        text = str(node["text"])
        latex = _latex_text(text)
        font = node.get("font")
        if font is not None:
            return f"{_FONT_LATEX[str(font)]}{{{latex}}}"
        return latex if len(text) == 1 else rf"\mathrm{{{latex}}}"
    if tag == "mn":
        return str(node["text"]).replace(MINUS, "-")
    if tag == "mo":
        text = str(node["text"])
        if text in _OPERATOR_LATEX:
            return _OPERATOR_LATEX[text]
        if text.isalpha():
            return rf"\operatorname{{{text}}}"
        return text
    if tag == "mtext":
        return rf"\text{{{node['text']}}}"
    if tag == "mspace":
        return r"\;"
    if tag in {"math", "mrow"}:
        return _latex_row(children)
    if tag == "msub":
        return f"{{{to_latex(children[0])}}}_{{{to_latex(children[1])}}}"
    if tag == "msup":
        return f"{{{to_latex(children[0])}}}^{{{to_latex(children[1])}}}"
    if tag == "msubsup":
        return (
            f"{{{to_latex(children[0])}}}_{{{to_latex(children[1])}}}^{{{to_latex(children[2])}}}"
        )
    if tag == "mfrac":
        return rf"\frac{{{to_latex(children[0])}}}{{{to_latex(children[1])}}}"
    if tag == "msqrt":
        return rf"\sqrt{{{_latex_row(children)}}}"
    if tag == "mover":
        return rf"\bar{{{to_latex(children[0])}}}"
    if tag == "munder":
        base = children[0]
        if base.get("tag") == "mo" and str(base.get("text")).isalpha():
            return rf"\{base['text']}_{{{to_latex(children[1])}}}"
        return rf"\underset{{{to_latex(children[1])}}}{{{to_latex(base)}}}"
    if tag == "mtable":
        rows = [
            " & ".join(_latex_row(cell.get("children", ())) for cell in entry.get("children", ()))
            for entry in children
        ]
        return r"\begin{matrix}" + r" \\ ".join(rows) + r"\end{matrix}"
    raise ValueError(f"No LaTeX form for <{tag}>.")


def _latex_row(children: Iterable[Mapping[str, Any]]) -> str:
    parts: list[str] = []
    for child in children:
        if child.get("tag") == "mo" and (child.get("attrs") or {}).get("stretchy") == "true":
            fence = str(child["text"])
            form = (child.get("attrs") or {}).get("form")
            latex_fence = {"{": r"\{", "}": r"\}"}.get(fence, fence)
            parts.append((r"\left" if form == "prefix" else r"\right") + latex_fence)
            continue
        parts.append(to_latex(child))
    return " ".join(part for part in parts if part)
