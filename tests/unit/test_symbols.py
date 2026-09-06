"""The registered display symbols.

What is being protected here is the rule that a symbol is *looked up*, never
typed. A Greek letter written inline into a service declaration or a figure
label is unreviewable against
``docs/standards/terminology_and_symbol_registry.md``, and drifts silently. So
the tests below check the two things that make the lookup worth having: that a
name which is not registered fails loudly, and that every registered entry is
complete enough to render in both of the places PyTex renders symbols.
"""

from __future__ import annotations

import re

import pytest

from pytex.core.symbols import (
    Symbol,
    has_symbol,
    registered_symbols,
    symbol,
    symbol_latex,
    symbol_text,
)


class TestLookup:
    def test_a_registered_symbol_returns_all_of_its_forms(self) -> None:
        entry = symbol("phi_1")
        assert isinstance(entry, Symbol)
        assert entry.name == "phi_1"
        assert entry.text == "φ₁"
        assert entry.latex == r"\varphi_1"

    def test_the_first_bunge_angle_uses_the_variant_phi(self) -> None:
        # The registry fixes \varphi rather than \phi for the Bunge angles. It
        # is not a stylistic preference: \phi_1 and \varphi_1 render as
        # different glyphs, and a figure that mixes them reads as two
        # quantities.
        assert symbol_latex("phi_1") == r"\varphi_1"
        assert symbol_latex("phi_2") == r"\varphi_2"

    def test_an_unregistered_name_is_refused_with_the_near_matches(self) -> None:
        with pytest.raises(KeyError) as excinfo:
            symbol("phi1")
        message = excinfo.value.args[0]
        assert "not a registered symbol" in message
        # The usual cause is one character, so the message must offer the
        # candidates rather than only reporting the miss.
        assert "phi_1" in message

    def test_has_symbol_agrees_with_lookup(self) -> None:
        assert has_symbol("Phi")
        assert not has_symbol("definitely_not_registered")

    def test_the_registry_cannot_be_mutated_by_a_caller(self) -> None:
        with pytest.raises(TypeError):
            registered_symbols()["injected"] = symbol("Phi")  # type: ignore[index]


class TestEveryEntry:
    """Completeness rules, applied to the whole registry rather than samples."""

    @pytest.mark.parametrize("name", sorted(registered_symbols()))
    def test_every_form_is_present_and_non_empty(self, name: str) -> None:
        entry = symbol(name)
        assert entry.name == name
        assert entry.text.strip()
        assert entry.latex.strip()
        assert entry.meaning.strip()

    @pytest.mark.parametrize("name", sorted(registered_symbols()))
    def test_the_meaning_reads_as_a_sentence(self, name: str) -> None:
        meaning = symbol(name).meaning
        assert meaning[0].isupper() or meaning[0].isdigit()
        assert meaning.endswith(".")

    @pytest.mark.parametrize("name", sorted(registered_symbols()))
    def test_the_mathtext_form_carries_no_delimiters(self, name: str) -> None:
        # Callers wrap it themselves — a matplotlib label needs one pair of
        # dollars around the whole expression, and an entry that brought its
        # own would produce ``$$\varphi_1$$``, which renders as literal text.
        assert "$" not in symbol_latex(name)

    @pytest.mark.parametrize("name", sorted(registered_symbols()))
    def test_the_unicode_form_is_short_enough_to_label_a_control(self, name: str) -> None:
        # The whole point of a symbol on a control is that it is narrower than
        # the words. Four characters is `2θ` with room to spare; anything
        # longer is an abbreviation, not a symbol, and belongs in the label.
        assert len(symbol_text(name)) <= 4

    @pytest.mark.parametrize("name", sorted(registered_symbols()))
    def test_the_name_is_typable_in_a_declaration(self, name: str) -> None:
        # Names are written by hand in service declarations, so they stay ASCII
        # and free of anything a keyboard makes awkward.
        assert re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name), name


class TestAgreementWithTheProseRegistry:
    """The module and the standards document must describe the same symbols."""

    def test_the_bunge_triple_is_written_the_way_the_document_writes_it(self) -> None:
        # ``docs/standards/terminology_and_symbol_registry.md`` states the
        # triple as $(\phi_1, \Phi, \phi_2)$. The three entries must therefore
        # exist and be distinct; a registry with two of them, or with the same
        # glyph twice, would let a form label two angles identically.
        triple = [symbol_text(name) for name in ("phi_1", "Phi", "phi_2")]
        assert len(set(triple)) == 3

    def test_the_holder_tilts_and_the_lattice_angles_are_separate_entries(self) -> None:
        # Both are written alpha and beta, and they are unrelated quantities.
        # Registering them separately is what lets a help string and a test say
        # which one a control means, even though the glyph is the same.
        assert symbol("alpha_tilt").meaning != symbol("alpha").meaning
        assert symbol("alpha_tilt").text == symbol("alpha").text
