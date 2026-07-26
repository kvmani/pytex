"""Human-facing crystallographic notation: brackets, overbars, and reciprocal stars.

This module is the single place PyTex turns crystallographic quantities into
text. Every figure label, `describe()` string, docstring example, and generated
documentation line goes through it, so the repository cannot drift into writing
`(1 -1 0)` in one plot and `(1bar 1 0)` in another.

The conventions implemented here are fixed centrally in
`docs/standards/notation_and_conventions.md` and anchored to the IUCr
*International Tables for Crystallography*.

Bracket families
----------------

============================  ==================  ===================
quantity                      specific            symmetry family
============================  ==================  ===================
lattice plane                 ``(hkl)``           ``{hkl}``
lattice direction             ``[uvw]``           ``<uvw>``
Miller-Bravais plane          ``(hkil)``          ``{hkil}``
Miller-Bravais direction      ``[uvtw]``          ``<uvtw>``
============================  ==================  ===================

Negative indices and separators
-------------------------------

A negative index is written with an **overbar** over the digit, not with a
leading minus sign, in publication-facing output. The ``"mathtext"`` style
emits ``\\bar{1}``; the ``"plain"`` style falls back to ``-1`` because a plain
terminal has no overbar.

Components are concatenated only when that is unambiguous. Bare concatenation is
the classical form — ``(110)`` — but it fails as soon as a component needs more
than one character: ``[1-10]`` could be read as ``1, -1, 0`` or as ``1, -10``,
and ``(1210)`` as ``1, 2, 1, 0`` or ``12, 1, 0``. A separator is therefore
inserted whenever any component is negative in plain style, or any component has
more than one digit. Mathtext uses a thin space (``\\;``) so the result stays
typographically correct.

The reciprocal star
-------------------

The star marks the **basis**, never the indices:

- direct basis vectors are ``a, b, c``;
- **reciprocal basis vectors are** ``a*, b*, c*``, and so are the axis labels of
  a reciprocal-domain `pytex.core.frames.ReferenceFrame` and the reciprocal
  lattice parameters;
- Miller indices ``(hkl)`` are *already* reciprocal-basis components by
  definition, so they are **not** starred — starring them would name a
  different quantity;
- a reciprocal-lattice vector is ``g_hkl = h a* + k b* + l c*``.

Use `format_reciprocal_axis_label` (or `format_reciprocal_axis_labels`) rather
than appending ``"*"`` by hand: the helpers are idempotent, so re-labelling an
already-starred axis cannot produce ``a**``.

See also
--------
`pytex.core.frame_catalog.reciprocal_frame_for` : builds starred reciprocal frames.
`docs/standards/hexagonal_and_trigonal_conventions.md` : the four-index rules.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

__all__ = [
    "RECIPROCAL_STAR",
    "format_direction_family_indices",
    "format_direction_indices",
    "format_miller_indices",
    "format_plane_family_indices",
    "format_plane_indices",
    "format_reciprocal_axis_label",
    "format_reciprocal_axis_labels",
    "format_reciprocal_lattice_vector",
    "is_reciprocal_axis_label",
    "strip_reciprocal_star",
]

#: The character marking a reciprocal-space basis vector or axis (IUCr convention).
RECIPROCAL_STAR = "*"

_FAMILIES = ("direction", "plane")
_SCOPES = ("specific", "family")
_STYLES = ("mathtext", "plain")

# Opening/closing delimiters keyed by (family, scope). The family forms are the
# international standard: braces for a plane family, angle brackets for a
# direction family.
_DELIMITERS: dict[tuple[str, str], tuple[str, str]] = {
    ("direction", "specific"): ("[", "]"),
    ("direction", "family"): ("<", ">"),
    ("plane", "specific"): ("(", ")"),
    ("plane", "family"): ("{", "}"),
}

# Mathtext needs escaped braces and \langle/\rangle rather than bare < >.
_MATHTEXT_DELIMITERS: dict[tuple[str, str], tuple[str, str]] = {
    ("direction", "specific"): ("[", "]"),
    ("direction", "family"): (r"\langle ", r" \rangle"),
    ("plane", "specific"): ("(", ")"),
    ("plane", "family"): (r"\{", r"\}"),
}


def _normalize_indices(indices: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(int(value) for value in indices)
    if len(normalized) not in {3, 4}:
        raise ValueError("Miller notation helpers require a 3-index or 4-index tuple.")
    return normalized


def _mathtext_component(value: int) -> str:
    if value < 0:
        return rf"\bar{{{abs(value)}}}"
    return str(value)


def _plain_component(value: int) -> str:
    return str(value)


def _needs_separator(values: tuple[int, ...], *, style: str) -> bool:
    """Whether concatenating these components would be ambiguous.

    Multi-digit components are always ambiguous when run together. In plain
    style a negative component is ambiguous too, because the minus sign reads
    as a separator: ``[1-10]`` could be ``1, -1, 0`` or ``1, -10``. Mathtext
    renders negatives as overbars, so only the digit count matters there.
    """

    if any(abs(value) > 9 for value in values):
        return True
    return style == "plain" and any(value < 0 for value in values)


def _validate_style(style: str) -> str:
    if style not in _STYLES:
        raise ValueError(f"style must be one of {', '.join(_STYLES)}.")
    return style


def format_miller_indices(
    indices: Sequence[int],
    *,
    family: str,
    style: str = "mathtext",
    scope: str = "specific",
) -> str:
    """Render Miller (or Miller-Bravais) indices in canonical bracket notation.

    What it does
        Wraps an index tuple in the internationally standard delimiters for its
        kind and renders negative components with an overbar in ``"mathtext"``.

    When to use it
        For every human-facing rendering of indices — plot labels, `describe()`
        prose, generated documentation. Never format indices inline; that is how
        a repository ends up with three spellings of the same plane.

    Parameters
    ----------
    indices:
        Three (Miller) or four (Miller-Bravais) integer indices.
    family:
        ``"plane"`` or ``"direction"``.
    style:
        ``"mathtext"`` (matplotlib/LaTeX math, overbars for negatives, wrapped
        in ``$``) or ``"plain"`` (bare ASCII with leading minus signs).
    scope:
        ``"specific"`` for one plane or direction — ``(hkl)`` / ``[uvw]`` — or
        ``"family"`` for the symmetry-related set — ``{hkl}`` / ``<uvw>``.

    Returns
    -------
    str
        The formatted indices.

    Raises
    ------
    ValueError
        If ``indices`` is not length 3 or 4, or ``family``, ``style``, or
        ``scope`` is not a recognized value.

    Examples
    --------
    ``format_miller_indices((1, 1, 0), family="plane", style="plain")`` gives
    ``"(110)"``. A negative component forces a separator, so
    ``format_miller_indices((1, -1, 0), family="plane", style="plain")`` gives
    ``"(1 -1 0)"``, and with ``scope="family"`` it gives ``"{1 -1 0}"``.
    """

    normalized = _normalize_indices(indices)
    if family not in _FAMILIES:
        raise ValueError(f"family must be one of {', '.join(_FAMILIES)}.")
    if scope not in _SCOPES:
        raise ValueError(f"scope must be one of {', '.join(_SCOPES)}.")
    _validate_style(style)

    if style == "mathtext":
        opener, closer = _MATHTEXT_DELIMITERS[(family, scope)]
        separator = r"\;" if _needs_separator(normalized, style=style) else ""
        payload = separator.join(_mathtext_component(value) for value in normalized)
        return f"${opener}{payload}{closer}$"
    opener, closer = _DELIMITERS[(family, scope)]
    separator = " " if _needs_separator(normalized, style=style) else ""
    payload = separator.join(_plain_component(value) for value in normalized)
    return f"{opener}{payload}{closer}"


def format_direction_indices(indices: Sequence[int], *, style: str = "mathtext") -> str:
    """Render one specific lattice direction as ``[uvw]`` (or ``[uvtw]``)."""

    return format_miller_indices(indices, family="direction", style=style)


def format_plane_indices(indices: Sequence[int], *, style: str = "mathtext") -> str:
    """Render one specific lattice plane as ``(hkl)`` (or ``(hkil)``)."""

    return format_miller_indices(indices, family="plane", style=style)


def format_direction_family_indices(
    indices: Sequence[int], *, style: str = "mathtext"
) -> str:
    """Render a symmetry-related direction family as ``<uvw>`` (or ``<uvtw>``).

    Use this whenever the quantity is the whole orbit rather than one member —
    a slip-direction family, a multiplicity count, a fibre axis family.
    """

    return format_miller_indices(indices, family="direction", style=style, scope="family")


def format_plane_family_indices(indices: Sequence[int], *, style: str = "mathtext") -> str:
    """Render a symmetry-related plane family as ``{hkl}`` (or ``{hkil}``).

    Use this for slip-plane families, powder-reflection families, and any
    multiplicity-bearing quantity, where naming a single member would misstate
    the science.
    """

    return format_miller_indices(indices, family="plane", style=style, scope="family")


# --------------------------------------------------------------------------- #
# Reciprocal-space notation
# --------------------------------------------------------------------------- #


def strip_reciprocal_star(label: str) -> str:
    """Remove a trailing reciprocal star from an axis label, if present.

    Makes the starring helpers idempotent, so a label that has already been
    starred cannot acquire a second star when it passes through another layer.
    """

    text = str(label).strip()
    while text.endswith(RECIPROCAL_STAR):
        text = text[: -len(RECIPROCAL_STAR)].rstrip()
    return text


def is_reciprocal_axis_label(label: str) -> bool:
    """Whether an axis label already carries the reciprocal star."""

    return str(label).strip().endswith(RECIPROCAL_STAR)


def format_reciprocal_axis_label(label: str, *, style: str = "plain") -> str:
    """Render a reciprocal-space axis label with the IUCr star.

    What it does
        Turns ``"a"`` into ``"a*"`` (plain) or ``"$a^{*}$"`` (mathtext),
        idempotently — an already-starred label is returned in the requested
        style rather than starred twice.

    When to use it
        For the axis labels of any reciprocal-domain frame or basis, and for any
        prose or figure text naming a reciprocal basis vector. The star marks the
        **basis**: do not apply it to Miller indices, which are already
        reciprocal-basis components.

    Parameters
    ----------
    label:
        The direct-space axis label (``"a"``) or an already-starred one.
    style:
        ``"plain"`` for ``a*``; ``"mathtext"`` for ``$a^{*}$``.

    Returns
    -------
    str
        The starred label.
    """

    _validate_style(style)
    base = strip_reciprocal_star(label)
    if style == "mathtext":
        return f"${base}^{{{RECIPROCAL_STAR}}}$"
    return f"{base}{RECIPROCAL_STAR}"


def format_reciprocal_axis_labels(
    labels: Iterable[str], *, style: str = "plain"
) -> tuple[str, ...]:
    """Star every label in an axis triad; see `format_reciprocal_axis_label`.

    Returns
    -------
    tuple[str, ...]
        The starred labels, in input order.
    """

    return tuple(format_reciprocal_axis_label(label, style=style) for label in labels)


def format_reciprocal_lattice_vector(
    indices: Sequence[int], *, style: str = "mathtext"
) -> str:
    """Render a reciprocal-lattice vector as ``g_hkl``.

    What it does
        Produces the standard subscripted symbol for the reciprocal-lattice
        vector belonging to a set of Miller indices — ``g_{110}`` — with
        overbars on negative indices in ``"mathtext"``.

    When to use it
        For diffraction-facing labels where the quantity is the scattering
        vector itself rather than the plane: SAED spot annotations, Ewald-sphere
        constructions, structure-factor tables.

    Notes
    -----
    ``g_hkl = h a* + k b* + l c*``: the vector is a combination of the
    **starred** reciprocal basis vectors, which is why the indices themselves
    carry no star.
    """

    normalized = _normalize_indices(indices)
    _validate_style(style)
    if style == "mathtext":
        separator = r"\;" if _needs_separator(normalized, style=style) else ""
        payload = separator.join(_mathtext_component(value) for value in normalized)
        return f"$g_{{{payload}}}$"
    separator = " " if _needs_separator(normalized, style=style) else ""
    payload = separator.join(_plain_component(value) for value in normalized)
    return f"g_{payload}"
