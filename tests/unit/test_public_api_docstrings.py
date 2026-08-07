"""The docstring ratchet on the public API surface.

`AGENTS.md` and `docs/standards/documentation_architecture.md` require every
documented public surface to state its purpose, when to use it, its inputs, and
its outputs. These tests enforce the precondition for that contract: that the
surface is documented at all. They are deliberately structural — they assert
presence and minimum substance, not wording — so they cannot be satisfied by a
placeholder and cannot go stale when prose is improved.

An export added without a docstring fails here rather than reaching users
undocumented.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import pytex

# Section headers of the NumPy docstring style used across the repository, plus
# the repository's own prose headings. A docstring that opens one of these must
# be reachable as a member docstring, so they are excluded from the
# "first line is a sentence" check below.
_MINIMUM_SUMMARY_WORDS = 3


def _docstring_of(member: Any) -> str:
    if isinstance(member, property):
        return (member.fget.__doc__ or "") if member.fget is not None else ""
    if isinstance(member, staticmethod | classmethod):
        return member.__func__.__doc__ or ""
    if callable(member):
        return member.__doc__ or ""
    return ""


def _public_class_members(obj: type) -> list[tuple[str, Any]]:
    """Public methods, properties, and class methods defined on ``obj`` itself.

    Inherited members are excluded: they are documented where they are defined,
    and re-checking them here would demand redundant docstrings on subclasses.
    """

    members: list[tuple[str, Any]] = []
    for name, member in vars(obj).items():
        if name.startswith("_"):
            continue
        if isinstance(member, property | staticmethod | classmethod) or callable(member):
            members.append((name, member))
    return members


def _exported_surfaces() -> list[tuple[str, Any]]:
    """Every public name and class member reachable through ``pytex.__all__``."""

    surfaces: list[tuple[str, Any]] = []
    for name in pytex.__all__:
        obj = getattr(pytex, name)
        if not (inspect.isclass(obj) or inspect.isfunction(obj)):
            continue
        surfaces.append((name, obj))
        if inspect.isclass(obj):
            surfaces.extend(
                (f"{name}.{member_name}", member)
                for member_name, member in _public_class_members(obj)
            )
    return surfaces


_SURFACES = _exported_surfaces()


def test_public_api_surface_is_non_trivial() -> None:
    """Guard the guard: a broken collector must not silently pass everything."""

    assert len(_SURFACES) > 500


@pytest.mark.parametrize("label,surface", _SURFACES, ids=[label for label, _ in _SURFACES])
def test_every_public_surface_has_a_docstring(label: str, surface: Any) -> None:
    docstring = (
        surface.__doc__ or "" if inspect.isclass(surface) or inspect.isfunction(surface)
        else _docstring_of(surface)
    )
    assert docstring.strip(), (
        f"Public surface 'pytex.{label}' has no docstring. Every name reachable "
        "through pytex.__all__ is part of the documented product surface; see "
        "the docstring contract in docs/standards/documentation_architecture.md."
    )


@pytest.mark.parametrize("label,surface", _SURFACES, ids=[label for label, _ in _SURFACES])
def test_every_public_docstring_opens_with_a_summary(label: str, surface: Any) -> None:
    """The first line must be a real summary, not a section header or a stub."""

    docstring = (
        surface.__doc__ or "" if inspect.isclass(surface) or inspect.isfunction(surface)
        else _docstring_of(surface)
    )
    summary = docstring.strip().splitlines()[0].strip()
    assert not summary.startswith(("Parameters", "Returns", "Purpose", "----")), (
        f"'pytex.{label}' opens its docstring with a section header rather than a "
        "one-line summary of what the surface produces."
    )
    assert len(summary.split()) >= _MINIMUM_SUMMARY_WORDS, (
        f"'pytex.{label}' has a {len(summary.split())}-word docstring summary; a "
        "placeholder is not documentation."
    )
