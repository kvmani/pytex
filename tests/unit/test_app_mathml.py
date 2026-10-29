# ruff: noqa: RUF001
"""The equation layer: a whitelisted MathML tree and its HTML and LaTeX forms.

The expectations are typographic identities rather than recorded output: a
negative index is overbarred, a minus is a minus sign, a registered symbol is
the registry's, and the tree a result carries can only hold tags the browser
renderer accepts.
"""

from __future__ import annotations

import pytest

from pytex.app import mathml as mm
from pytex.app.export import result_to_html, result_to_markdown, result_to_xlsx
from pytex.app.results import AppResult, Column, ResultEquation, ResultStage, ResultTable


class TestLeaves:
    def test_a_negative_number_uses_the_minus_sign(self) -> None:
        assert mm.number(-0.5, 2)["text"] == "−0.50"

    def test_negative_zero_is_printed_as_zero(self) -> None:
        assert mm.number(-0.00001, 4)["text"] == "0.0000"

    def test_non_finite_numbers_are_refused(self) -> None:
        with pytest.raises(ValueError, match="finite"):
            mm.number(float("nan"))

    def test_an_unknown_font_is_refused(self) -> None:
        with pytest.raises(ValueError, match="font"):
            mm.mi("x", font="fraktur")


class TestRegistry:
    def test_the_first_bunge_angle_is_varphi_sub_one(self) -> None:
        node = mm.symbol("phi_1")
        assert node["tag"] == "msub"
        assert node["children"][0]["text"] == "φ"
        assert node["children"][1]["text"] == "1"

    def test_an_upper_case_greek_symbol_is_upright(self) -> None:
        assert mm.symbol("Phi")["font"] == "upright"

    def test_bold_and_scripts_parse(self) -> None:
        node = mm.from_mathtext(r"\mathbf{S}_{p,k}^{\mathsf{T}}")
        assert node["tag"] == "msubsup"
        assert node["children"][0] == {"tag": "mi", "text": "S", "font": "bold"}
        assert node["children"][2]["font"] == "sans"

    def test_unsupported_commands_are_refused(self) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            mm.from_mathtext(r"\int x")


class TestMiller:
    def test_a_negative_index_is_overbarred(self) -> None:
        node = mm.miller((1, -1, 0))
        overbars = [child for child in node["children"] if child["tag"] == "mover"]
        assert len(overbars) == 1
        assert overbars[0]["children"][0]["text"] == "1"
        assert r"\bar{1}" in mm.to_latex(node)

    def test_a_direction_family_is_in_angle_brackets(self) -> None:
        markup = mm.to_markup(mm.miller((1, 1, 0), family="direction", scope="family"))
        assert "⟨" in markup and "⟩" in markup


class TestValidation:
    def test_a_script_tag_cannot_enter_a_tree(self) -> None:
        with pytest.raises(ValueError, match="not allowed"):
            mm.validate({"tag": "script", "text": "alert(1)"})

    def test_an_event_attribute_cannot_enter_a_tree(self) -> None:
        with pytest.raises(ValueError, match="not allowed"):
            mm.validate({"tag": "mi", "text": "x", "attrs": {"onclick": "x"}})

    def test_markup_escapes_text(self) -> None:
        markup = mm.to_markup(mm.equation(mm.mtext("<b>&</b>")))
        assert "<b>" not in markup
        assert "&lt;b&gt;&amp;&lt;/b&gt;" in markup

    def test_an_equation_must_be_rooted_at_math(self) -> None:
        with pytest.raises(ValueError, match="math"):
            ResultEquation(math=mm.mi("x"))


class TestLatex:
    def test_a_matrix_is_a_bracketed_matrix_environment(self) -> None:
        latex = mm.to_latex(
            mm.equation(mm.symbol("phi_1"), "=", mm.matrix([[1, -2], [3, 4]], digits=1))
        )
        assert latex == (
            r"{\varphi}_{1} = \left[ \begin{matrix}1.0 & -2.0 \\ 3.0 & 4.0\end{matrix} \right]"
        )

    def test_a_minimum_over_a_set_is_a_subscripted_min(self) -> None:
        node = mm.under(mm.mo("min"), mm.row(mm.mi("x"), "∈", mm.mi("X")))
        assert mm.to_latex(node) == r"\min_{x \in X}"


def _result() -> AppResult:
    equation = ResultEquation(
        lead="The measured rotation is",
        math=mm.equation(mm.from_mathtext(r"\mathbf{V}"), "=", mm.matrix([[1, 0], [0, -1]])),
        tag="3",
        note="A proper rotation.",
    )
    stage = ResultStage(
        key="step",
        title="Step 1",
        summary="Compute it.",
        section="derivation",
        equations=(equation,),
    )
    table = ResultTable(columns=(Column("k", "K"),), rows=({"k": 1},))
    return AppResult(title="T", summary="S", stages=(stage,), table=table)


class TestExports:
    def test_a_stage_serializes_its_equations(self) -> None:
        payload = _result().to_json()
        equation = payload["stages"][0]["equations"][0]
        assert equation["tag"] == "3"
        assert equation["math"]["tag"] == "math"

    def test_the_html_report_embeds_mathml(self) -> None:
        page = result_to_html(_result().to_json()).decode()
        assert '<math display="block" xmlns="http://www.w3.org/1998/Math/MathML">' in page
        assert "(3)" in page
        assert "Step-by-step derivation" in page

    def test_the_markdown_report_writes_numbered_latex(self) -> None:
        text = result_to_markdown(_result().to_json()).decode()
        assert "$$" in text
        assert r"\tag{3}" in text
        assert "## Step-by-step derivation" in text

    def test_the_workbook_records_the_equation(self) -> None:
        assert result_to_xlsx(_result().to_json())[:2] == b"PK"
