"""Files a user opens in the workbench, on their way to a library reader.

Purpose
-------
The library's importers read *paths* — `read_ang`, `read_ctf`,
`read_xrdml_pole_figure` all open a file — and the workbench has a browser at one
end and, often, a server on another machine at the other. This module is the one
place that gap is crossed: the browser sends the file's text as an ordinary
parameter, and :func:`uploaded_file` materialises it as a temporary path for the
reader to open.

Why the text and not a multipart upload
---------------------------------------
An operation is a JSON request. Giving one operation a second, binary transport
would mean two request paths, two validation stories and two places for an error
to come back from, for files that are text and rarely more than a few tens of
megabytes. The text rides in the request the operation already has, and
`MAX_REQUEST_BYTES` in the server bounds it.

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
                       suffixes=(".ang", ".ctf")) as (path, name):
        result = read_ang(path)
"""

from __future__ import annotations

import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from pytex.app.errors import InvalidInputError

__all__ = ["describe_upload", "uploaded_file", "uploaded_name_and_text"]


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
        If the payload is not an object, has no name or text, or names a file
        kind this operation does not read. Refusing an unreadable extension here
        rather than letting the reader fail on the content is what turns "list
        index out of range" into "this operation reads .ang and .ctf files".
    """

    if not isinstance(payload, Mapping):
        raise InvalidInputError(
            "No file has been opened.",
            field=field,
            hint=f"Choose a {_readable_list(suffixes)} file.",
        )
    name = str(payload.get("name") or "").strip()
    text = payload.get("text")
    if not isinstance(text, str) or not text.strip():
        raise InvalidInputError(
            f"The file {name or 'that was opened'} is empty.",
            field=field,
            hint="Check that the file downloaded completely before opening it.",
        )
    suffix = Path(name).suffix.lower()
    if suffix not in {value.lower() for value in suffixes}:
        raise InvalidInputError(
            f"{name or 'That file'} is not a file this reads: its extension is "
            f"{suffix or 'missing'}.",
            field=field,
            hint=f"This operation reads {_readable_list(suffixes)} files.",
        )
    return name, text


@contextmanager
def uploaded_file(
    payload: Any,
    *,
    field: str,
    suffixes: Sequence[str],
) -> Iterator[tuple[Path, str]]:
    """Materialise an upload as a temporary path, and remove it afterwards.

    Yields
    ------
    tuple of (Path, str)
        The path a library reader can open, and the file name the user knows it
        by — which is the one to put in a title or an error, because the
        temporary name means nothing to them.
    """

    name, text = uploaded_name_and_text(payload, field=field, suffixes=suffixes)
    suffix = Path(name).suffix or suffixes[0]
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, encoding="utf-8", newline="", delete=False
    )
    try:
        handle.write(text)
        handle.close()
        yield Path(handle.name), name
    finally:
        Path(handle.name).unlink(missing_ok=True)


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
