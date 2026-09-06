"""The registered display symbols, in the one form a screen can render.

Purpose
-------
:mod:`pytex.core.notation` writes Miller indices; this module writes everything
else the symbol registry names --- the Bunge angles, the lattice parameters, the
diffraction angles, the tilt axes --- as the international literature writes
them, so a control in the workbench is labelled ``φ₁`` rather than ``phi1`` and
a figure axis and a form label cannot disagree about which quantity they mean.

Why a module rather than a literal
----------------------------------
``AGENTS.md`` fixes symbol meaning centrally
(``docs/standards/terminology_and_symbol_registry.md``) and forbids formatting
notation inline. A Greek letter typed into a service declaration is inline
formatting: it is unreviewable against the registry, it is spelled three
different ways within a week (``phi1``, ``Phi 1``, ``φ1``), and nothing fails
when it drifts. Naming the *registered symbol* instead makes the spelling a
lookup, and an unregistered name a construction-time error rather than a typo
that ships.

Why the label survives
----------------------
A symbol is a shorthand for people who already know the quantity, and the
workbench is also used by people meeting it for the first time. So a symbol
never replaces the words: it is what the control *shows*, and the words are what
the control is *called* --- carried into the accessible name, the tooltip and the
help text. Compactness on screen is not bought with a form that cannot be read
by a newcomer or by a screen reader.

Two rendered forms
------------------
``text`` is Unicode, for HTML labels, plain-text output and log messages.
``latex`` is a mathtext body without delimiters, for matplotlib axis labels and
MyST prose. Both are given explicitly rather than one derived from the other,
because the mapping is not mechanical: ``φ₁`` is ``\\varphi_1`` and not
``\\phi_1``, the registry having fixed the variant form for the first Bunge
angle.

Examples
--------
>>> from pytex.core.symbols import symbol, symbol_text
>>> symbol_text("phi_1")
'φ₁'
>>> symbol("two_theta").latex
'2\\\\theta'
>>> symbol("alpha_tilt").meaning
'Primary (α) holder tilt of a double-tilt TEM stage.'
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

__all__ = [
    "Symbol",
    "has_symbol",
    "registered_symbols",
    "symbol",
    "symbol_latex",
    "symbol_text",
]


@dataclass(frozen=True)
class Symbol:
    """One registered symbol in its renderable forms.

    Attributes
    ----------
    name : str
        Registry key. Lower-case ASCII with underscores, so it can be typed in
        a service declaration and matched against the registry mechanically.
    text : str
        Unicode display form, e.g. ``"φ₁"``. Safe in HTML text content, in a
        terminal, and in a log line.
    latex : str
        Mathtext body without ``$`` delimiters, e.g. ``r"\\varphi_1"``.
    meaning : str
        One sentence, taken from the terminology and symbol registry, saying
        what the quantity is. Read by tests that check the two documents agree.
    """

    name: str
    text: str
    latex: str
    meaning: str


def _entry(name: str, text: str, latex: str, meaning: str) -> tuple[str, Symbol]:
    return name, Symbol(name=name, text=text, latex=latex, meaning=meaning)


#: Every symbol the application may name. Adding one here is the act of
#: registering it; ``docs/standards/terminology_and_symbol_registry.md`` remains
#: the prose definition and the two are checked against each other in tests.
_SYMBOLS: Mapping[str, Symbol] = MappingProxyType(
    dict(
        (
            # -- orientation -------------------------------------------------
            _entry(
                "phi_1", "φ₁", r"\varphi_1", "First Bunge Euler angle, about the specimen z axis."
            ),
            _entry("Phi", "Φ", r"\Phi", "Second Bunge Euler angle, about the new x axis."),
            _entry("phi_2", "φ₂", r"\varphi_2", "Third Bunge Euler angle, about the new z axis."),
            _entry("omega", "ω", r"\omega", "Rotation angle of an axis-angle pair."),
            # -- lattice -----------------------------------------------------
            _entry("a", "a", "a", "Direct-lattice edge length along the a axis."),
            _entry("b", "b", "b", "Direct-lattice edge length along the b axis."),
            _entry("c", "c", "c", "Direct-lattice edge length along the c axis."),
            _entry("alpha", "α", r"\alpha", "Lattice angle between the b and c axes."),
            _entry("beta", "β", r"\beta", "Lattice angle between the c and a axes."),
            _entry("gamma", "γ", r"\gamma", "Lattice angle between the a and b axes."),
            _entry("d_spacing", "d", "d", "Interplanar spacing of a reflecting plane family."),
            # -- diffraction -------------------------------------------------
            _entry("wavelength", "λ", r"\lambda", "Radiation wavelength."),
            _entry(
                "camera_constant", "Lλ", r"L\lambda", "Camera constant of a diffraction pattern."
            ),
            _entry("camera_length", "L", "L", "Effective camera length of a diffraction geometry."),
            _entry(
                "accelerating_voltage", "V", "V", "Accelerating voltage of the electron source."
            ),
            _entry("s_1", "s₁", "s_{1}", "Deviation parameter of the first CBED fringe minimum."),
            _entry("s_2", "s₂", "s_{2}", "Deviation parameter of the second CBED fringe minimum."),
            _entry("s_3", "s₃", "s_{3}", "Deviation parameter of the third CBED fringe minimum."),
            _entry("s_4", "s₄", "s_{4}", "Deviation parameter of the fourth CBED fringe minimum."),
            _entry("s_5", "s₅", "s_{5}", "Deviation parameter of the fifth CBED fringe minimum."),
            # -- stage and detector -------------------------------------------
            _entry(
                "alpha_tilt", "α", r"\alpha", "Primary (α) holder tilt of a double-tilt TEM stage."
            ),
            _entry(
                "beta_tilt", "β", r"\beta", "Secondary (β) holder tilt of a double-tilt TEM stage."
            ),
            _entry(
                "pattern_centre_x",
                "x*",
                "x^{*}",
                "Pattern-centre abscissa in detector-width fractions.",
            ),
            _entry(
                "pattern_centre_y",
                "y*",
                "y^{*}",
                "Pattern-centre ordinate in detector-height fractions.",
            ),
            _entry(
                "detector_distance", "z*", "z^{*}", "Detector distance in detector-width fractions."
            ),
            # -- texture -------------------------------------------------------
            _entry("halfwidth", "ψ", r"\psi", "Halfwidth of an ODF kernel."),
            # -- misc ----------------------------------------------------------
        )
    )
)


def registered_symbols() -> Mapping[str, Symbol]:
    """Return the whole registry, keyed by name.

    When to use
    -----------
    In tests and in tooling that must enumerate what is registered --- the
    manifest test that checks every symbol a parameter names exists, and the
    documentation test that checks the registry and this module agree. Callers
    that want one symbol should use :func:`symbol`, which fails usefully.

    Returns
    -------
    Mapping[str, Symbol]
        A read-only view. The registry is immutable at runtime: a symbol is
        added by editing this module, which is the reviewable act.
    """

    return _SYMBOLS


def has_symbol(name: str) -> bool:
    """Whether ``name`` is a registered symbol."""

    return name in _SYMBOLS


def symbol(name: str) -> Symbol:
    """Look one symbol up by its registered name.

    Parameters
    ----------
    name : str
        Registry key, e.g. ``"phi_1"``.

    Returns
    -------
    Symbol
        Its Unicode form, its mathtext form, and its one-line meaning.

    Raises
    ------
    KeyError
        If the name is not registered. The message lists the near matches,
        because the usual cause is a spelling (``"phi1"`` for ``"phi_1"``) and
        the usual fix is one character.

    Examples
    --------
    >>> symbol("Phi").text
    'Φ'
    """

    try:
        return _SYMBOLS[name]
    except KeyError:
        stem = name.lower().replace("_", "")
        near = sorted(key for key in _SYMBOLS if key.lower().replace("_", "").startswith(stem[:3]))
        hint = f" Did you mean: {', '.join(near)}?" if near else ""
        raise KeyError(
            f"{name!r} is not a registered symbol. Register it in "
            f"pytex.core.symbols and in docs/standards/"
            f"terminology_and_symbol_registry.md before using it.{hint}"
        ) from None


def symbol_text(name: str) -> str:
    """The Unicode display form of a registered symbol."""

    return symbol(name).text


def symbol_latex(name: str) -> str:
    """The mathtext body of a registered symbol, without ``$`` delimiters."""

    return symbol(name).latex
