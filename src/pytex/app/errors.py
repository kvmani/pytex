"""Errors that reach a user's screen.

An application error is not a stack trace. It is a sentence the person at the
keyboard can act on, and the whole point of this module is to make that sentence
mandatory rather than optional.

Every failure a service raises deliberately carries three things: a stable
machine-readable ``code`` the frontend can branch on, a ``message`` written for a
researcher rather than for a maintainer, and a ``hint`` naming the next thing to
try. Anything raised that is *not* a :class:`ServiceError` is a bug in PyTex, and
the transport layer reports it as such rather than dressing it up as user error.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "DependencyMissingError",
    "InvalidInputError",
    "ServiceError",
    "UnknownOperationError",
    "UnsupportedRequestError",
]


class ServiceError(Exception):
    """A failure the user caused and can correct.

    Purpose
    -------
    Distinguishes "you asked for something impossible" from "PyTex is broken".
    The transport layer turns a :class:`ServiceError` into a ``400``-class JSON
    envelope carrying the message verbatim, and turns anything else into a
    ``500`` with a generic message plus a logged traceback.

    Parameters
    ----------
    message : str
        What went wrong, in one sentence, addressed to the user. State the
        offending value when there is one.
    code : str
        Stable dotted identifier, e.g. ``"input.invalid"``. The frontend may
        branch on it; it must not change once released.
    hint : str, optional
        What to do about it. Almost always worth writing.
    details : dict, optional
        Structured context — the offending field name, allowed values, bounds.
        Rendered by the frontend next to the control that caused the failure.

    Examples
    --------
    >>> raise ServiceError(
    ...     "Lattice parameter a must be positive; got -1.0 Angstrom.",
    ...     code="input.invalid",
    ...     hint="Enter the edge length of the unit cell in Angstrom.",
    ...     details={"field": "a"},
    ... )  # doctest: +SKIP
    """

    #: HTTP status the transport layer should use for this class of error.
    status: int = 400

    def __init__(
        self,
        message: str,
        *,
        code: str = "service.error",
        hint: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.code = code
        self.hint = hint
        self.details: dict[str, Any] = dict(details or {})

    def to_json(self) -> dict[str, Any]:
        """Return the JSON body the frontend renders."""

        payload: dict[str, Any] = {"code": self.code, "message": self.message}
        if self.hint is not None:
            payload["hint"] = self.hint
        if self.details:
            payload["details"] = self.details
        return payload


class InvalidInputError(ServiceError):
    """A parameter is missing, malformed, or out of range."""

    status = 400

    def __init__(
        self,
        message: str,
        *,
        field: str | None = None,
        hint: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        merged = dict(details or {})
        if field is not None:
            merged["field"] = field
        super().__init__(message, code="input.invalid", hint=hint, details=merged)


class UnknownOperationError(ServiceError):
    """The requested operation is not registered."""

    status = 404

    def __init__(self, operation: str, *, known: tuple[str, ...] = ()) -> None:
        super().__init__(
            f"No operation named {operation!r} is registered.",
            code="operation.unknown",
            hint="Fetch /api/manifest for the list of available operations.",
            details={"operation": operation, "known": list(known)},
        )


class UnsupportedRequestError(ServiceError):
    """The request is well formed but asks for something this build cannot do."""

    status = 409

    def __init__(self, message: str, *, hint: str | None = None) -> None:
        super().__init__(message, code="request.unsupported", hint=hint)


class DependencyMissingError(ServiceError):
    """An optional dependency is needed for this request and is not installed.

    Optional extras are optional in PyTex, so an operation that needs one fails
    with a named install command rather than an :class:`ImportError` traceback.
    """

    status = 501

    def __init__(self, package: str, *, purpose: str, extra: str | None = None) -> None:
        install = f"pip install pytex[{extra}]" if extra else f"pip install {package}"
        super().__init__(
            f"{purpose} requires the optional dependency {package!r}, which is not installed.",
            code="dependency.missing",
            hint=f"Install it with: {install}",
            details={"package": package, "extra": extra, "install": install},
        )
