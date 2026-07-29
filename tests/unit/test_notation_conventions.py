"""Tests for the central crystallographic notation surface.

These pin the internationally standard conventions fixed in
`docs/standards/notation_and_conventions.md`: bracket families, overbars for
negative indices, unambiguous separators, and the reciprocal star. They are
convention tests, not implementation tests — a failure here means the repository
would print something a crystallographer would read as a different quantity.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from pytex.core.notation import (
    RECIPROCAL_STAR,
    format_direction_family_indices,
    format_direction_indices,
    format_miller_indices,
    format_plane_family_indices,
    format_plane_indices,
    format_reciprocal_axis_label,
    format_reciprocal_axis_labels,
    format_reciprocal_lattice_vector,
    is_reciprocal_axis_label,
    strip_reciprocal_star,
)

# ---------------------------------------------------------------------------
# Bracket families
# ---------------------------------------------------------------------------


def test_specific_and_family_brackets_follow_the_international_convention() -> None:
    assert format_plane_indices((1, 1, 0), style="plain") == "(110)"
    assert format_plane_family_indices((1, 1, 0), style="plain") == "{110}"
    assert format_direction_indices((1, 1, 0), style="plain") == "[110]"
    assert format_direction_family_indices((1, 1, 0), style="plain") == "<110>"


def test_four_index_miller_bravais_forms_use_the_same_brackets() -> None:
    assert format_plane_indices((1, 1, -2, 1), style="plain") == "(1 1 -2 1)"
    assert format_plane_family_indices((0, 0, 0, 1), style="plain") == "{0001}"
    assert format_direction_family_indices((1, 1, -2, 0), style="plain") == "<1 1 -2 0>"


def test_mathtext_family_delimiters_are_escaped_for_the_renderer() -> None:
    assert format_plane_family_indices((1, 1, 1)) == r"$\{111\}$"
    assert format_direction_family_indices((1, 1, 0)) == r"$\langle 110 \rangle$"


# ---------------------------------------------------------------------------
# Overbars and separators
# ---------------------------------------------------------------------------


def test_negative_indices_use_an_overbar_in_mathtext() -> None:
    assert format_plane_indices((1, -1, 0)) == r"$(1\bar{1}0)$"
    assert format_direction_indices((1, 1, -2, 0)) == r"$[11\bar{2}0]$"


def test_plain_style_separates_components_when_concatenation_is_ambiguous() -> None:
    # "[1-10]" could be read as [1, -1, 0] or as [1, -10].
    assert format_direction_indices((1, -1, 0), style="plain") == "[1 -1 0]"
    # No negatives and single digits: the classical concatenated form is safe.
    assert format_direction_indices((1, 1, 0), style="plain") == "[110]"


def test_multi_digit_components_are_separated_in_both_styles() -> None:
    # "(1210)" could be read as (1, 2, 1, 0) or (12, 1, 0).
    assert format_plane_indices((12, 1, 0), style="plain") == "(12 1 0)"
    assert format_plane_indices((12, 1, 0)) == r"$(12\;1\;0)$"


def test_mathtext_does_not_separate_single_digit_negatives() -> None:
    # The overbar already delimits the component, so no separator is needed.
    assert r"\;" not in format_plane_indices((1, -1, 0))


@pytest.mark.parametrize(
    "indices",
    [(1, 1, 0), (1, -1, 0), (12, 1, 0), (1, 1, -2, 1), (2, -1, -1, 0), (0, 0, 0, -1)],
)
@pytest.mark.parametrize(
    "formatter",
    [
        format_plane_indices,
        format_plane_family_indices,
        format_direction_indices,
        format_direction_family_indices,
    ],
)
def test_every_mathtext_form_parses_in_matplotlib(indices, formatter) -> None:  # type: ignore[no-untyped-def]
    """A label the renderer cannot parse is a broken figure, not a broken string."""

    figure, axes = plt.subplots()
    try:
        axes.set_title(formatter(indices))
        figure.canvas.draw()
    finally:
        plt.close(figure)


def test_formatter_rejects_bad_arguments() -> None:
    with pytest.raises(ValueError, match="3-index or 4-index"):
        format_plane_indices((1, 1))
    with pytest.raises(ValueError, match="family must be one of"):
        format_miller_indices((1, 1, 1), family="pole")
    with pytest.raises(ValueError, match="scope must be one of"):
        format_miller_indices((1, 1, 1), family="plane", scope="orbit")
    with pytest.raises(ValueError, match="style must be one of"):
        format_plane_indices((1, 1, 1), style="html")


# ---------------------------------------------------------------------------
# The reciprocal star
# ---------------------------------------------------------------------------


def test_reciprocal_axis_labels_carry_the_star() -> None:
    assert format_reciprocal_axis_label("a") == "a*"
    assert format_reciprocal_axis_labels(("a", "b", "c")) == ("a*", "b*", "c*")
    assert format_reciprocal_axis_label("a", style="mathtext") == "$a^{*}$"


def test_starring_is_idempotent() -> None:
    # Passing an already-starred label through another layer must not yield "a**".
    assert format_reciprocal_axis_label("a*") == "a*"
    assert format_reciprocal_axis_label(format_reciprocal_axis_label("a")) == "a*"
    assert format_reciprocal_axis_labels(("a*", "b*", "c*")) == ("a*", "b*", "c*")


def test_star_helpers_report_and_strip() -> None:
    assert is_reciprocal_axis_label("a*")
    assert not is_reciprocal_axis_label("a")
    assert strip_reciprocal_star("a*") == "a"
    assert strip_reciprocal_star("a") == "a"
    assert RECIPROCAL_STAR == "*"


def test_miller_indices_are_never_starred() -> None:
    """The star marks the basis, not the indices.

    ``(hkl)`` are already reciprocal-basis components by definition; starring
    them would name a different quantity. This guards against a well-meaning
    future change that "makes reciprocal quantities starred" too broadly.
    """

    for text in (
        format_plane_indices((1, 1, 1), style="plain"),
        format_plane_family_indices((1, 1, 1), style="plain"),
        format_plane_indices((1, 1, 1)),
    ):
        assert RECIPROCAL_STAR not in text


def test_reciprocal_lattice_vector_uses_the_g_symbol() -> None:
    assert format_reciprocal_lattice_vector((1, 1, 0)) == "$g_{110}$"
    assert format_reciprocal_lattice_vector((1, 1, 0), style="plain") == "g_110"
    assert format_reciprocal_lattice_vector((1, -1, 0)) == r"$g_{1\bar{1}0}$"
    assert format_reciprocal_lattice_vector((1, -1, 0), style="plain") == "g_1 -1 0"


# ---------------------------------------------------------------------------
# Repository-wide consistency
# ---------------------------------------------------------------------------


def test_reciprocal_frames_built_anywhere_carry_starred_axes() -> None:
    from pytex.core.frame_catalog import CRYSTAL_FRAME, reciprocal_frame_for
    from pytex.core.lattice import Lattice

    assert reciprocal_frame_for(CRYSTAL_FRAME).axes == ("a*", "b*", "c*")

    lattice = Lattice(
        a=3.6,
        b=3.6,
        c=3.6,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=CRYSTAL_FRAME,
    )
    assert lattice.reciprocal_basis().frame.axes == ("a*", "b*", "c*")


def test_direct_basis_axes_are_not_starred() -> None:
    from pytex.core.frame_catalog import CRYSTAL_FRAME
    from pytex.core.lattice import Lattice

    lattice = Lattice(
        a=3.6,
        b=3.6,
        c=3.6,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=CRYSTAL_FRAME,
    )
    assert lattice.direct_basis().frame.axes == ("a", "b", "c")


def test_pole_figure_titles_match_what_is_actually_plotted() -> None:
    """A family pole figure is ``{hkl}``; a single-pole one is ``(hkl)``."""

    from pytex.plotting.builders import _pole_figure_label

    class _Stub:
        def __init__(self, family: bool) -> None:
            self.includes_symmetry_family = family
            self.pole = type(
                "P", (), {"miller": type("M", (), {"indices": (1, 1, 1)})()}
            )()

    assert _pole_figure_label(_Stub(True)) == "{111}"  # type: ignore[arg-type]
    assert _pole_figure_label(_Stub(False)) == "(111)"  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Policy: the conventions are fixed centrally and not re-implemented inline
# ---------------------------------------------------------------------------


def test_notation_conventions_are_documented_in_the_standards() -> None:
    from pathlib import Path

    repo_root = Path(__file__).resolve().parents[2]
    standard = (repo_root / "docs/standards/notation_and_conventions.md").read_text(
        encoding="utf-8"
    )
    # The standard writes notation as rendered mathematics, so each concept is
    # pinned by any of its accepted spellings (plain text or LaTeX).
    concepts = (
        ("a*", "a^{*}"),
        ("{hkl}", r"\{hkl\}"),
        ("<uvw>", r"\langle uvw \rangle"),
        ("(hkl)",),
        ("[uvw]",),
        ("overbar",),
        ("zone law",),
    )
    for spellings in concepts:
        assert any(
            token in standard.lower() or token in standard for token in spellings
        ), spellings

    agents = (repo_root / "AGENTS.md").read_text(encoding="utf-8")
    assert "pytex.core.notation" in agents
    assert "reciprocal" in agents.lower()

    registry = (repo_root / "docs/standards/terminology_and_symbol_registry.md").read_text(
        encoding="utf-8"
    )
    assert "reciprocal star" in registry
    assert "symmetry family" in registry


def test_no_module_formats_miller_indices_inline() -> None:
    """Index brackets must come from `pytex.core.notation`, not string literals.

    An inline ``"(" + " ".join(str(int(v)) ...) + ")"`` is how a repository
    drifts into several spellings of the same plane. The check targets that
    integer-index pattern specifically: bracketing *decimal* components (an
    irrational zone axis, say) is not index notation and stays allowed. The
    composite-SAED module keeps its own documented spot formatter, which has
    domain-specific spacing rules, and the notation module is obviously exempt.
    """

    import re
    from pathlib import Path

    repo_root = Path(__file__).resolve().parents[2]
    exempt = {"notation.py", "composite_saed.py"}
    # An opening bracket literal, joined integer components, a closing bracket.
    pattern = re.compile(r"""["'][\(\[]["']\s*\+.*str\(int\(""", re.DOTALL)
    offenders: list[str] = []
    for path in sorted((repo_root / "src" / "pytex").rglob("*.py")):
        if path.name in exempt:
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            if pattern.search(line):
                offenders.append(f"{path.relative_to(repo_root)}: {line.strip()}")
    assert not offenders, "Inline index formatting found: " + "; ".join(offenders)


def test_no_module_appends_a_reciprocal_star_by_hand() -> None:
    """Starring must go through the idempotent helper, or ``a**`` becomes possible."""

    from pathlib import Path

    repo_root = Path(__file__).resolve().parents[2]
    offenders: list[str] = []
    for path in sorted((repo_root / "src" / "pytex").rglob("*.py")):
        if path.name == "notation.py":
            continue
        text = path.read_text(encoding="utf-8")
        if 'f"{label}*"' in text or "f'{label}*'" in text:
            offenders.append(str(path.relative_to(repo_root)))
    assert not offenders, "Hand-rolled reciprocal starring found in: " + ", ".join(offenders)
