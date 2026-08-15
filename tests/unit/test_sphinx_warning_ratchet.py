from __future__ import annotations

from scripts.check_sphinx_warnings import warning_count_from_output


def test_sphinx_warning_count_uses_the_authoritative_summary() -> None:
    output = (
        "source.py:1: WARNING: first\nsource.py:2: WARNING: second\nbuild succeeded, 3 warnings.\n"
    )
    assert warning_count_from_output(output) == 3


def test_sphinx_warning_count_falls_back_to_warning_lines() -> None:
    output = "source.py:1: WARNING: first\nplain output\nsource.py:2: WARNING: second\n"
    assert warning_count_from_output(output) == 2


def test_sphinx_warning_count_accepts_a_warning_free_build() -> None:
    assert warning_count_from_output("build succeeded.\n") == 0
