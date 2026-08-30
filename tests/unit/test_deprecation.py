from __future__ import annotations

import warnings

import pytest

from pytex._deprecation import deprecated, deprecation_message, warn_deprecated

#: The removal version these examples pretend to be scheduled for.
#:
#: Deliberately unreachable. It was a plausible next version, and when the real
#: version reached it the release-metadata guard -- which requires the version
#: literal to appear only in `_version.py` -- failed on this file. The guard was
#: right; a version-shaped literal in a test is exactly what it is looking for.
_EXAMPLE_REMOVAL = "99.0.0"


def test_deprecation_message_shape() -> None:
    message = deprecation_message(
        "pytex.old_name", since="0.2.3", removal=_EXAMPLE_REMOVAL, replacement="pytex.new_name"
    )
    assert "pytex.old_name is deprecated since PyTex 0.2.3" in message
    assert f"may be removed in PyTex {_EXAMPLE_REMOVAL}" in message
    assert "Use pytex.new_name instead." in message
    bare = deprecation_message("pytex.old_name", since="0.2.3", removal=_EXAMPLE_REMOVAL)
    assert "Use" not in bare


def test_decorator_warns_and_preserves_behavior() -> None:
    @deprecated(since="0.2.3", removal=_EXAMPLE_REMOVAL, replacement="new_add")
    def old_add(a: int, b: int) -> int:
        """Add two numbers."""

        return a + b

    with pytest.warns(DeprecationWarning, match="old_add is deprecated since PyTex 0.2.3"):
        assert old_add(2, 3) == 5
    # wraps metadata preserved; policy marker attached
    assert old_add.__doc__ == "Add two numbers."
    assert _EXAMPLE_REMOVAL in old_add.__deprecated__  # type: ignore[attr-defined]


def test_warn_deprecated_implicates_caller() -> None:
    def deprecated_surface() -> None:
        warn_deprecated("deprecated_surface", since="0.2.3", removal=_EXAMPLE_REMOVAL, stacklevel=3)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        deprecated_surface()
    assert len(caught) == 1
    # stacklevel=3 attributes the warning to this test file, not _deprecation.py
    assert caught[0].filename == __file__
