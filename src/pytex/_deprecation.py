"""Uniform deprecation warnings implementing the API stability standard.

Every PyTex deprecation goes through this module so warnings carry the same
grep-able shape (name, since, removal, replacement), point at the caller via
``stacklevel``, and can be asserted in tests. See
``docs/standards/api_stability_and_deprecation.md`` for the policy: stable
surfaces warn for at least two minor releases before removal.
"""

from __future__ import annotations

import functools
import warnings
from collections.abc import Callable
from typing import Any, TypeVar

_F = TypeVar("_F", bound=Callable[..., Any])


def deprecation_message(
    name: str,
    *,
    since: str,
    removal: str,
    replacement: str | None = None,
) -> str:
    """Build the canonical deprecation message for ``name``."""

    message = (
        f"{name} is deprecated since PyTex {since} and may be removed in "
        f"PyTex {removal}."
    )
    if replacement:
        message += f" Use {replacement} instead."
    return message


def warn_deprecated(
    name: str,
    *,
    since: str,
    removal: str,
    replacement: str | None = None,
    stacklevel: int = 3,
) -> None:
    """Emit the canonical `DeprecationWarning` for a deprecated surface.

    ``stacklevel=3`` implicates the caller of the deprecated public function
    (this helper -> deprecated function -> user code); adjust when the call
    chain is deeper.
    """

    warnings.warn(
        deprecation_message(name, since=since, removal=removal, replacement=replacement),
        DeprecationWarning,
        stacklevel=stacklevel,
    )


def deprecated(
    *,
    since: str,
    removal: str,
    replacement: str | None = None,
) -> Callable[[_F], _F]:
    """Decorator marking a callable as deprecated per the stability policy.

    The wrapped callable keeps its behavior and signature; every call emits
    the canonical `DeprecationWarning` naming the replacement (when given)
    and the first release allowed to remove the surface.
    """

    def decorate(func: _F) -> _F:
        qualified = f"{func.__module__}.{func.__qualname__}"

        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            warn_deprecated(
                qualified,
                since=since,
                removal=removal,
                replacement=replacement,
                stacklevel=3,
            )
            return func(*args, **kwargs)

        wrapper.__deprecated__ = deprecation_message(  # type: ignore[attr-defined]
            qualified, since=since, removal=removal, replacement=replacement
        )
        return wrapper  # type: ignore[return-value]

    return decorate


__all__ = ["deprecated", "deprecation_message", "warn_deprecated"]
