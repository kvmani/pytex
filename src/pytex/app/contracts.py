"""The wire format shared by every shell.

One envelope, used by the HTTP server, the desktop bridge, and the tests, so
there is exactly one definition of what a successful call and a failed call look
like. Keeping it here rather than in the server means a shell that is not HTTP at
all still speaks the same language.

The envelope
------------
Success::

    {"ok": true, "operation": "calc.plane_angles", "result": {...}, "log": [...]}

Failure::

    {"ok": false, "operation": "calc.plane_angles", "log": [...],
     "error": {"code": "input.invalid", "message": "...", "hint": "..."}}

``ok`` is always present and is the only field a client must check. The result
object is whatever the operation documents in its manifest entry.

``log`` is the narration this call produced, in
:mod:`pytex.app.logbook` wire form — always present, often empty. It travels
with the envelope rather than being polled separately so that a message about a
calculation can never arrive before the calculation's own result, or after the
user has moved on to the next one.
"""

from __future__ import annotations

import json
import time
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.errors import ServiceError
from pytex.app.logbook import APP_LOG, LogRecord, collecting
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


def _log_payload(records: Sequence[LogRecord] | None) -> list[dict[str, Any]]:
    """Render captured records for the wire, oldest first."""

    return [record.to_json() for record in (records or ())]


def success_envelope(
    operation: str, result: Any, *, log: Sequence[LogRecord] | None = None
) -> dict[str, Any]:
    """Wrap a result in the success envelope."""

    return {
        "ok": True,
        "operation": operation,
        "result": to_jsonable(result),
        "log": _log_payload(log),
    }


def error_envelope(
    operation: str, error: ServiceError, *, log: Sequence[LogRecord] | None = None
) -> dict[str, Any]:
    """Wrap a deliberate failure in the error envelope."""

    return {
        "ok": False,
        "operation": operation,
        "error": error.to_json(),
        "log": _log_payload(log),
    }


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

    Notes
    -----
    Every call narrates itself into :data:`~pytex.app.logbook.APP_LOG`: one
    record when it starts, one when it finishes or fails. That happens here
    rather than in each shell so the console tells the same story whether the
    call arrived over HTTP, through the desktop bridge, or from a test — and
    rather than in each service, so no operation can be added without appearing
    in the log.
    """

    active = registry if registry is not None else REGISTRY
    started = time.monotonic()
    # The registered title, not the dotted id: the console is read by the person
    # who pressed the button, and they pressed "Powder XRD pattern" rather than
    # "xrd.powder_pattern". An unknown id has no title, and saying so plainly is
    # better than inventing one.
    try:
        title = active.get(operation).title
    except ServiceError:
        title = operation
    with collecting() as records:
        APP_LOG.info(f"{title} started.", source=operation)
        try:
            result = active.call(operation, request)
        except ServiceError as error:
            # A ServiceError is a message already written for the user, so the
            # log repeats it verbatim rather than paraphrasing: the console and
            # the toast beside the control must not describe the same rejection
            # in two different sentences.
            detail: dict[str, Any] = {"code": error.code}
            if error.hint is not None:
                detail["hint"] = error.hint
            if "field" in error.details:
                detail["field"] = error.details["field"]
            APP_LOG.error(error.message, source=operation, detail=detail)
            return error_envelope(operation, error, log=records), error.status
        elapsed_ms = (time.monotonic() - started) * 1000.0
        APP_LOG.success(
            f"{title} completed in {elapsed_ms:.0f} ms.",
            source=operation,
            detail={"duration_ms": round(elapsed_ms, 1)},
        )
        return success_envelope(operation, result, log=records), 200
