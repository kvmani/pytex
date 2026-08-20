"""Files a user opens in the workbench, on their way to a library reader.

Purpose
-------
The library's importers read *paths* — `read_ang`, `read_ctf`,
`read_xrdml_pole_figure` all open a file — and the workbench has a browser at one
end and, often, a server on another machine at the other. This module is the one
place that gap is crossed: the browser sends the file's contents as an ordinary
parameter, and :func:`uploaded_file` materialises them as a temporary path for
the reader to open.

Why the text and not a multipart upload
---------------------------------------
An operation is a JSON request. Giving one operation a second, binary transport
would mean two request paths, two validation stories and two places for an error
to come back from, for files that are text and rarely more than a few tens of
megabytes. The text rides in the request the operation already has, and
`MAX_REQUEST_BYTES` in the server bounds it.

Binary files ride the same request
----------------------------------
Not every scan format is text: an EDAX OIM scan is HDF5, and `read_oh5` opens
it. Rather than add the second transport this module exists to avoid, a binary
file arrives base64-encoded in the same JSON field, as ``{"name": ...,
"data_base64": ...}``, and :func:`uploaded_file` writes the decoded bytes.
One request path, one validation story; the cost is the 4/3 inflation base64
carries, which `MAX_REQUEST_BYTES` bounds along with everything else. A scan
saved with its diffraction patterns is far larger than that limit, which is a
property of the file rather than of this transport: strip the patterns on
export, or read the file with `pytex.adapters.read_scan` in a script.

Why a temporary file rather than a parser that takes text
----------------------------------------------------------
Because the readers are the library's, and they are the specification. Adding a
text-taking twin of each would be a second implementation of the same format to
keep in step, which is exactly the duplication the adapters exist to prevent.
The file is written into the system temporary directory and removed as soon as
the reader has finished with it; nothing is retained between requests.

When and where to use it
------------------------
In any service handler that accepts a user data file. Declare the parameter as
an :class:`~pytex.app.registry.ObjectParameter`, then::

    with uploaded_file(request["scan_file"], field="scan_file",
                       suffixes=SCAN_FILE_SUFFIXES) as (path, name):
        result = read_scan(path)
"""

from __future__ import annotations

import base64
import binascii
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from pytex.app.errors import InvalidInputError

__all__ = [
    "describe_upload",
    "uploaded_bytes",
    "uploaded_file",
    "uploaded_name_and_text",
]


def uploaded_name_and_text(
    payload: Any,
    *,
    field: str,
    suffixes: Sequence[str],
) -> tuple[str, str]:
    """Validate an upload payload and return its file name and text.

    Parameters
    ----------
    payload : Any
        What the browser sent: ``{"name": "map.ang", "text": "..."}``.
    field : str
        The parameter name, so an error can be shown beside the right control.
    suffixes : Sequence[str]
        Accepted file extensions, lower case and including the dot.

    Returns
    -------
    tuple of (str, str)
        The file name as the user knows it, and its contents.

    Raises
    ------
    InvalidInputError
        If the payload is not an object, names a file kind this operation does
        not read, or carries no text.
    """

    name = _validated_name(payload, field=field, suffixes=suffixes)
    text = payload.get("text")
    if not isinstance(text, str) or not text.strip():
        raise InvalidInputError(
            f"The file {name or 'that was opened'} is empty.",
            field=field,
            hint="Check that the file downloaded completely before opening it.",
        )
    return name, text


def uploaded_bytes(
    payload: Any,
    *,
    field: str,
    suffixes: Sequence[str],
) -> tuple[str, bytes]:
    """Validate a base64 upload payload and return its file name and bytes.

    The binary counterpart of :func:`uploaded_name_and_text`, for formats a
    reader opens as bytes rather than as text — an EDAX OIM HDF5 scan above all.

    Parameters
    ----------
    payload : Any
        What the browser sent: ``{"name": "map.oh5", "data_base64": "..."}``.
    field : str
        The parameter name, so an error can be shown beside the right control.
    suffixes : Sequence[str]
        Accepted file extensions, lower case and including the dot.

    Returns
    -------
    tuple of (str, bytes)
        The file name as the user knows it, and its decoded contents.

    Raises
    ------
    InvalidInputError
        If the payload is not an object, names a file kind this operation does
        not read, or carries a field that is empty or not valid base64.
    """

    name = _validated_name(payload, field=field, suffixes=suffixes)
    encoded = payload.get("data_base64")
    if not isinstance(encoded, str) or not encoded.strip():
        raise InvalidInputError(
            f"The file {name or 'that was opened'} is empty.",
            field=field,
            hint="Check that the file downloaded completely before opening it.",
        )
    try:
        data = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise InvalidInputError(
            f"The file {name or 'that was opened'} did not arrive intact.",
            field=field,
            hint="Open it again; the browser encodes it for the request each time.",
        ) from error
    if not data:
        raise InvalidInputError(
            f"The file {name or 'that was opened'} is empty.",
            field=field,
            hint="Check that the file downloaded completely before opening it.",
        )
    return name, data


@contextmanager
def uploaded_file(
    payload: Any,
    *,
    field: str,
    suffixes: Sequence[str],
) -> Iterator[tuple[Path, str]]:
    """Materialise an upload as a temporary path, and remove it afterwards.

    Takes a text payload (``text``) or a base64 binary one (``data_base64``),
    whichever the browser sent for this file's format, and writes the same kind
    of temporary file either way — so a caller passes the path to its reader
    without knowing which of the two arrived.

    Yields
    ------
    tuple of (Path, str)
        The path a library reader can open, and the file name the user knows it
        by — which is the one to put in a title or an error, because the
        temporary name means nothing to them.
    """

    if isinstance(payload, Mapping) and "data_base64" in payload:
        name, data = uploaded_bytes(payload, field=field, suffixes=suffixes)
    else:
        name, text = uploaded_name_and_text(payload, field=field, suffixes=suffixes)
        # Written as bytes rather than as text, which is what the text branch
        # meant by `newline=""`: the reader must see the file the browser sent,
        # not a copy with this platform's line endings substituted into it.
        data = text.encode("utf-8")
    handle = tempfile.NamedTemporaryFile(
        mode="wb", suffix=Path(name).suffix or suffixes[0], delete=False
    )
    try:
        handle.write(data)
        handle.close()
        yield Path(handle.name), name
    finally:
        Path(handle.name).unlink(missing_ok=True)


def _validated_name(payload: Any, *, field: str, suffixes: Sequence[str]) -> str:
    """The file's name, once it is an object naming a kind this operation reads.

    Refusing an unreadable extension here rather than letting the reader fail on
    the content is what turns "list index out of range" into "this operation
    reads .ang, .ctf, .oh5 or .h5 files".
    """

    if not isinstance(payload, Mapping):
        raise InvalidInputError(
            "No file has been opened.",
            field=field,
            hint=f"Choose a {_readable_list(suffixes)} file.",
        )
    name = str(payload.get("name") or "").strip()
    suffix = Path(name).suffix.lower()
    if suffix not in {value.lower() for value in suffixes}:
        raise InvalidInputError(
            f"{name or 'That file'} is not a file this reads: its extension is "
            f"{suffix or 'missing'}.",
            field=field,
            hint=f"This operation reads {_readable_list(suffixes)} files.",
        )
    return name


def describe_upload(name: str, text: str) -> dict[str, Any]:
    """A provenance block naming what was read and how much of it.

    Attached to a result so that a figure made from a user's own file says which
    file, rather than looking like the built-in example beside it.
    """

    return {
        "name": name,
        "characters": len(text),
        "lines": text.count("\n") + 1,
    }


def _readable_list(suffixes: Sequence[str]) -> str:
    values = list(suffixes)
    if len(values) == 1:
        return values[0]
    return ", ".join(values[:-1]) + f" or {values[-1]}"
