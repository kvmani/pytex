"""The wire format shared by every shell.

One envelope, used by the HTTP server, the desktop bridge, and the tests, so
there is exactly one definition of what a successful call and a failed call look
like. Keeping it here rather than in the server means a shell that is not HTTP at
all still speaks the same language.

The envelope
------------
Success::

    {"ok": true, "operation": "calc.plane_angles", "result": {...}}

Failure::

    {"ok": false, "operation": "calc.plane_angles",
     "error": {"code": "input.invalid", "message": "...", "hint": "..."}}

``ok`` is always present and is the only field a client must check. The result
object is whatever the operation documents in its manifest entry.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

import numpy as np

from pytex.app.errors import ServiceError
from pytex.app.registry import REGISTRY, ServiceRegistry

__all__ = [
    "APP_ENVELOPE_SCHEMA",
    "dumps",
    "error_envelope",
    "execute",
    "success_envelope",
    "to_jsonable",
]

#: Schema identifier of the call envelope.
APP_ENVELOPE_SCHEMA = "pytex.app_envelope/1"


def to_jsonable(value: Any) -> Any:
    """Convert a result object into something :mod:`json` can serialise.

    Purpose
    -------
    Services compute with NumPy, and NumPy scalars and arrays are not JSON. This
    is the one place that conversion happens, so no service has to remember to
    call ``.tolist()`` and no ``float32`` leaks into a response.

    Non-finite floats are converted to ``None`` rather than to JavaScript's
    ``NaN`` literal, which is not valid JSON and which every strict parser
    rejects. A missing value is honest; an unparseable document is not.

    Parameters
    ----------
    value : object
        Any nesting of mappings, sequences, NumPy arrays and scalars, and
        objects exposing ``to_json()``.

    Returns
    -------
    object
        A structure of ``dict``, ``list``, ``str``, ``int``, ``float``, ``bool``
        and ``None`` only.
    """

    if value is None or isinstance(value, str | bool | int):
        return value
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, np.generic):
        return to_jsonable(value.item())
    if isinstance(value, np.ndarray):
        return [to_jsonable(item) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    to_json = getattr(value, "to_json", None)
    if callable(to_json):
        return to_jsonable(to_json())
    if isinstance(value, tuple | list | set | frozenset):
        return [to_jsonable(item) for item in value]
    return str(value)


def dumps(value: Any, *, indent: int | None = None) -> str:
    """Serialise a result to JSON text, converting NumPy types on the way."""

    return json.dumps(to_jsonable(value), indent=indent, allow_nan=False)


def success_envelope(operation: str, result: Any) -> dict[str, Any]:
    """Wrap a result in the success envelope."""

    return {"ok": True, "operation": operation, "result": to_jsonable(result)}


def error_envelope(operation: str, error: ServiceError) -> dict[str, Any]:
    """Wrap a deliberate failure in the error envelope."""

    return {"ok": False, "operation": operation, "error": error.to_json()}


def execute(
    operation: str,
    request: Mapping[str, Any] | None = None,
    *,
    registry: ServiceRegistry | None = None,
) -> tuple[dict[str, Any], int]:
    """Run one operation and return its envelope with an HTTP-style status.

    Purpose
    -------
    The single dispatch path every shell uses. Deliberate failures become an
    error envelope with the operation's own status; unexpected exceptions are
    re-raised, because they are PyTex defects and the transport layer must log
    them rather than present them as user error.

    Parameters
    ----------
    operation : str
        Registered operation identifier.
    request : mapping, optional
        Raw, unvalidated parameters as they arrived from the client.
    registry : ServiceRegistry, optional
        Defaults to the application-wide :data:`~pytex.app.registry.REGISTRY`.

    Returns
    -------
    tuple of (dict, int)
        The envelope and the status code to send with it.
    """

    active = registry if registry is not None else REGISTRY
    try:
        result = active.call(operation, request)
    except ServiceError as error:
        return error_envelope(operation, error), error.status
    return success_envelope(operation, result), 200
