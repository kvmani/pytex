from __future__ import annotations

import warnings

import pytest

from pytex._deprecation import deprecated, deprecation_message, warn_deprecated


def test_deprecation_message_shape() -> None:
    message = deprecation_message(
        "pytex.old_name", since="0.2.0", removal="0.4.0", replacement="pytex.new_name"
    )
    assert "pytex.old_name is deprecated since PyTex 0.2.0" in message
    assert "may be removed in PyTex 0.4.0" in message
    assert "Use pytex.new_name instead." in message
    bare = deprecation_message("pytex.old_name", since="0.2.0", removal="0.4.0")
    assert "Use" not in bare


def test_decorator_warns_and_preserves_behavior() -> None:
    @deprecated(since="0.2.0", removal="0.4.0", replacement="new_add")
    def old_add(a: int, b: int) -> int:
        """Add two numbers."""

        return a + b

    with pytest.warns(DeprecationWarning, match="old_add is deprecated since PyTex 0.2.0"):
        assert old_add(2, 3) == 5
    # wraps metadata preserved; policy marker attached
    assert old_add.__doc__ == "Add two numbers."
    assert "0.4.0" in old_add.__deprecated__  # type: ignore[attr-defined]


def test_warn_deprecated_implicates_caller() -> None:
    def deprecated_surface() -> None:
        warn_deprecated(
            "deprecated_surface", since="0.2.0", removal="0.4.0", stacklevel=3
        )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        deprecated_surface()
    assert len(caught) == 1
    # stacklevel=3 attributes the warning to this test file, not _deprecation.py
    assert caught[0].filename == __file__
