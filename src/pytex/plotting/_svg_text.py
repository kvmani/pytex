"""Arial/Helvetica text metrics for laying out SVG without a renderer.

PyTex generates and audits SVG figures offline, where no text-measurement API is
available. Guessing an average character width is not good enough: it is wrong
by 30% or more for capital- or narrow-letter-heavy strings, and every layout
decision built on it inherits that error.

This module carries the standard Helvetica advance-width table (units per 1000
em), which Arial matches closely by design — the two are metrically compatible,
which is exactly why the canonical font stack lists them together. Widths are
therefore accurate to within a percent or so for the strings these figures use,
which is enough to size a panel or to decide that a caption overflows its card.

Both the figure generator (`pytex.plotting.frame_diagrams`) and the layout
auditor (``scripts/audit_figure_text_layout.py``) measure through this module, so
generated and hand-authored figures are judged by one ruler.

References
----------
Adobe Systems, *Helvetica* AFM metrics, as distributed with the Core 14 PostScript
font set; Arial is metrically compatible with Helvetica by design.
"""

from __future__ import annotations

__all__ = ["BOLD_WIDTH_FACTOR", "DEFAULT_ADVANCE", "advance_widths", "text_width"]

#: Advance width used for any character absent from the table, in units/1000 em.
DEFAULT_ADVANCE = 556

#: Helvetica-Bold runs wider than the regular face; applying one factor is far
#: closer than ignoring weight, without carrying a second full table.
BOLD_WIDTH_FACTOR = 1.07

_WIDTHS: dict[str, int] = {
    " ": 278, "!": 278, '"': 355, "#": 556, "$": 556, "%": 889, "&": 667, "'": 191,
    "(": 333, ")": 333, "*": 389, "+": 584, ",": 278, "-": 333, ".": 278, "/": 278,
    "0": 556, "1": 556, "2": 556, "3": 556, "4": 556, "5": 556, "6": 556, "7": 556,
    "8": 556, "9": 556, ":": 278, ";": 278, "<": 584, "=": 584, ">": 584, "?": 556,
    "@": 1015,
    "A": 667, "B": 667, "C": 722, "D": 722, "E": 667, "F": 611, "G": 778, "H": 722,
    "I": 278, "J": 500, "K": 667, "L": 556, "M": 833, "N": 722, "O": 778, "P": 667,
    "Q": 778, "R": 722, "S": 667, "T": 611, "U": 722, "V": 667, "W": 944, "X": 667,
    "Y": 667, "Z": 611,
    "[": 278, "\\": 278, "]": 278, "^": 469, "_": 556, "`": 333,
    "a": 556, "b": 556, "c": 500, "d": 556, "e": 556, "f": 278, "g": 556, "h": 556,
    "i": 222, "j": 222, "k": 500, "l": 222, "m": 833, "n": 556, "o": 556, "p": 556,
    "q": 556, "r": 333, "s": 500, "t": 278, "u": 556, "v": 500, "w": 722, "x": 500,
    "y": 500, "z": 500,
    "{": 334, "|": 260, "}": 334, "~": 584,
    # Punctuation and symbols that appear in these figures.
    "\u00b0": 400,   # degree sign
    "\u00b7": 278,   # middle dot
    "\u2013": 556,   # en dash
    "\u2014": 1000,  # em dash
    "\u2018": 222, "\u2019": 222, "\u201c": 333, "\u201d": 333,  # curly quotes
    "\u2192": 800,   # rightwards arrow
    "\u2194": 1000,  # left-right arrow
    "\u2260": 584,   # not equal
    "\u2264": 584, "\u2265": 584,  # less/greater or equal
    "\u2208": 584,   # element of
    "\u22c5": 278,   # dot operator
    "\u00d7": 584,   # multiplication sign
    "\u221a": 549,   # square root
    "\u03b1": 556, "\u03b2": 556, "\u03b3": 500,  # alpha, beta, gamma
    "\u03c6": 556, "\u03a6": 778,                  # phi, Phi
    "\u03c9": 722, "\u03b8": 556, "\u0394": 612, "\u03c3": 556,  # omega, theta, Delta, sigma
}


def advance_widths() -> dict[str, int]:
    """A copy of the advance-width table, in units per 1000 em."""

    return dict(_WIDTHS)


def text_width(text: str, font_size: float, *, bold: bool = False) -> float:
    """Rendered width of a string in Arial/Helvetica, in SVG user units.

    What it does
        Sums per-character advance widths from the Helvetica metric table and
        scales by the font size, so a string of capitals measures wider than the
        same number of narrow lowercase letters — the distinction an average
        character width throws away.

    When to use it
        Whenever SVG layout has to be decided offline: sizing a panel to its
        caption, deciding whether a label fits inside its box, or placing two
        labels without overlap.

    Parameters
    ----------
    text:
        The string to measure.
    font_size:
        Font size in SVG user units.
    bold:
        Apply the bold-face widening factor.

    Returns
    -------
    float
        The estimated advance width. Accurate to about a percent for Latin text;
        leave a margin rather than butting elements against the result.

    Examples
    --------
    Capitals are much wider than narrow lowercase, which an average-width model
    cannot express: ``text_width("MMMM", 10)`` is far larger than
    ``text_width("llll", 10)``.
    """

    total = sum(_WIDTHS.get(character, DEFAULT_ADVANCE) for character in str(text))
    width = total * float(font_size) / 1000.0
    return width * BOLD_WIDTH_FACTOR if bold else width
