"""Getting results out: CSV, XLSX, and JSON, from any result at all.

Because every operation returns the same shape — a title, prose, an optional
table, the inputs it came from — export is written once here and works for every
operation that exists or ever will. No operation contains export code, and none
has to be remembered when a new format is added.

Three formats, three jobs:

- **CSV** is the lowest common denominator: one row per reported entity, full
  precision, openable by anything.
- **XLSX** is the same table plus a second sheet recording the inputs and the
  provenance, because a spreadsheet that has lost the conditions it was computed
  under is a spreadsheet nobody can defend in review.
- **JSON** is the whole result object, schema-tagged and round-trippable back
  into the application, so a figure in a paper can be regenerated from the file
  that produced it.

Why there is an ``.xlsx`` writer here
-------------------------------------
An ``.xlsx`` file is a zip of XML documents, and the subset needed to write one
sheet of strings and numbers is small enough to implement in a page and a half.
That is worth doing rather than adding ``openpyxl`` to the runtime, because the
application's deployment target is a host that may never reach PyPI (see
Decision 3 in ``docs/architecture/application_platform.md``). The writer is
deliberately minimal: no formatting, no formulas, no charts. It writes numbers as
numbers and text as text, which is the whole requirement for re-plottable data.
"""

from __future__ import annotations

import csv
import io
import json
import zipfile
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from typing import Any
from xml.sax.saxutils import escape

from pytex.app.contracts import to_jsonable
from pytex.app.errors import InvalidInputError

__all__ = [
    "EXPORT_FORMATS",
    "export_result",
    "result_to_csv",
    "result_to_json",
    "result_to_xlsx",
    "write_xlsx",
]

#: The formats any table-bearing result can be exported as.
EXPORT_FORMATS: dict[str, dict[str, str]] = {
    "csv": {
        "label": "CSV",
        "mime": "text/csv; charset=utf-8",
        "extension": "csv",
        "description": "One row per reported entity, at full precision.",
    },
    "xlsx": {
        "label": "Excel workbook",
        "mime": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        "extension": "xlsx",
        "description": "The table, plus a sheet recording the inputs it was computed from.",
    },
    "json": {
        "label": "JSON",
        "mime": "application/json; charset=utf-8",
        "extension": "json",
        "description": "The complete result, round-trippable back into the application.",
    },
}


def _table_of(result: Mapping[str, Any]) -> Mapping[str, Any]:
    table = result.get("table")
    if not table or not table.get("columns"):
        raise InvalidInputError(
            "This result has no table to export.",
            field="format",
            hint="Export it as JSON, which carries the whole result including its inputs.",
        )
    return dict(table)


def result_to_csv(result: Mapping[str, Any]) -> bytes:
    """Write the result table as CSV.

    Full precision, not display precision: the on-screen table rounds so it can
    be read, and a file that inherited that rounding would be useless for
    re-plotting.
    """

    table = _table_of(result)
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\r\n")
    writer.writerow(
        [
            f"{column['label']} ({column['units']})" if column.get("units") else column["label"]
            for column in table["columns"]
        ]
    )
    for row in table["rows"]:
        writer.writerow([_cell(row.get(column["key"])) for column in table["columns"]])
    return buffer.getvalue().encode("utf-8-sig")


def result_to_json(result: Mapping[str, Any]) -> bytes:
    """Write the whole result, including the inputs that produced it."""

    payload = dict(to_jsonable(result))
    payload["exported_utc"] = datetime.now(UTC).isoformat(timespec="seconds")
    return json.dumps(payload, indent=2).encode("utf-8")


def result_to_xlsx(result: Mapping[str, Any]) -> bytes:
    """Write the result table plus a provenance sheet as an ``.xlsx`` workbook."""

    table = _table_of(result)
    headers = [
        f"{column['label']} ({column['units']})" if column.get("units") else column["label"]
        for column in table["columns"]
    ]
    rows: list[list[Any]] = [
        [row.get(column["key"]) for column in table["columns"]] for row in table["rows"]
    ]

    provenance: list[list[Any]] = [
        ["Result", result.get("title", "")],
        ["Summary", result.get("summary", "")],
        ["Exported (UTC)", datetime.now(UTC).isoformat(timespec="seconds")],
    ]
    for note in result.get("notes", ()):
        provenance.append(["Note", note])
    for citation in result.get("citations", ()):
        provenance.append(["Source", citation])
    for key, value in sorted((result.get("inputs") or {}).items()):
        provenance.append(
            [f"Input: {key}", value if isinstance(value, str | int | float) else json.dumps(value)]
        )

    return write_xlsx(
        {
            "Data": {"headers": headers, "rows": rows},
            "Provenance": {"headers": ["Field", "Value"], "rows": provenance},
        }
    )


def export_result(result: Mapping[str, Any], *, fmt: str) -> tuple[bytes, str, str]:
    """Export a result in one of :data:`EXPORT_FORMATS`.

    Returns
    -------
    tuple of (bytes, str, str)
        The payload, its MIME type, and a suggested filename.
    """

    if fmt not in EXPORT_FORMATS:
        raise InvalidInputError(
            f"{fmt!r} is not an export format.",
            field="format",
            hint="Available: " + ", ".join(sorted(EXPORT_FORMATS)) + ".",
        )
    writers = {"csv": result_to_csv, "xlsx": result_to_xlsx, "json": result_to_json}
    payload = writers[fmt](result)
    spec = EXPORT_FORMATS[fmt]
    return (
        payload,
        spec["mime"],
        f"{_slug(str(result.get('title', 'pytex-result')))}.{spec['extension']}",
    )


def _slug(text: str) -> str:
    cleaned = "".join(character if character.isalnum() else "-" for character in text.lower())
    while "--" in cleaned:
        cleaned = cleaned.replace("--", "-")
    return cleaned.strip("-")[:60] or "pytex-result"


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


# --------------------------------------------------------------------------
# A minimal .xlsx writer
# --------------------------------------------------------------------------

_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
{sheets}
</Types>"""

_ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""


def write_xlsx(sheets: Mapping[str, Mapping[str, Any]]) -> bytes:
    """Write a minimal ``.xlsx`` workbook from plain rows.

    Purpose
    -------
    Produce a spreadsheet with no third-party dependency. An ``.xlsx`` file is a
    zip of XML parts; this writes the four that a reader requires plus one sheet
    part per sheet.

    Parameters
    ----------
    sheets : mapping
        Sheet name to ``{"headers": [...], "rows": [[...], ...]}``. Sheet names
        are truncated to Excel's 31-character limit and stripped of the
        characters Excel forbids, because a workbook that will not open is worse
        than a truncated tab label.

    Returns
    -------
    bytes
        The workbook.

    Notes
    -----
    Values are written inline rather than through a shared-strings table. That
    is a larger file for text-heavy sheets and a much simpler one to verify;
    for result tables, which are mostly numbers, the difference is small.
    """

    if not sheets:
        raise ValueError("A workbook needs at least one sheet.")

    names = [_sheet_name(name) for name in sheets]
    parts: dict[str, bytes] = {}

    overrides = "\n".join(
        f'<Override PartName="/xl/worksheets/sheet{index + 1}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for index in range(len(names))
    )
    parts["[Content_Types].xml"] = _CONTENT_TYPES.format(sheets=overrides).encode("utf-8")
    parts["_rels/.rels"] = _ROOT_RELS.encode("utf-8")

    sheet_entries = "".join(
        f'<sheet name="{escape(name)}" sheetId="{index + 1}" r:id="rId{index + 1}"/>'
        for index, name in enumerate(names)
    )
    parts["xl/workbook.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{sheet_entries}</sheets></workbook>"
    ).encode()

    relationships = "".join(
        f'<Relationship Id="rId{index + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{index + 1}.xml"/>'
        for index in range(len(names))
    )
    parts["xl/_rels/workbook.xml.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{relationships}</Relationships>"
    ).encode()

    for index, (_, content) in enumerate(sheets.items()):
        rows = [list(content.get("headers", []))] if content.get("headers") else []
        rows.extend([list(row) for row in content.get("rows", [])])
        parts[f"xl/worksheets/sheet{index + 1}.xml"] = _sheet_xml(rows)

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, payload in parts.items():
            archive.writestr(name, payload)
    return buffer.getvalue()


def _sheet_name(name: str) -> str:
    forbidden = set(r"[]:*?/\\")
    cleaned = "".join(" " if character in forbidden else character for character in str(name))
    return cleaned[:31] or "Sheet"


def _sheet_xml(rows: Sequence[Sequence[Any]]) -> bytes:
    body: list[str] = []
    for row_index, row in enumerate(rows, start=1):
        cells: list[str] = []
        for column_index, value in enumerate(row):
            reference = f"{_column_letter(column_index)}{row_index}"
            cells.append(_cell_xml(reference, value))
        body.append(f'<row r="{row_index}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f"<sheetData>{''.join(body)}</sheetData></worksheet>"
    ).encode()


def _cell_xml(reference: str, value: Any) -> str:
    if value is None or value == "":
        return f'<c r="{reference}"/>'
    if isinstance(value, bool):
        # Excel's boolean type exists, but a spreadsheet of "yes"/"no" is what a
        # reader expects from a column headed "Allowed".
        return _inline_string(reference, "yes" if value else "no")
    if isinstance(value, int | float):
        if not _finite(value):
            return f'<c r="{reference}"/>'
        return f'<c r="{reference}"><v>{value!r}</v></c>'
    return _inline_string(reference, str(value))


def _finite(value: float) -> bool:
    return value == value and value not in (float("inf"), float("-inf"))


def _inline_string(reference: str, text: str) -> str:
    return (
        f'<c r="{reference}" t="inlineStr"><is><t xml:space="preserve">{escape(text)}</t></is></c>'
    )


def _column_letter(index: int) -> str:
    letters = ""
    current = index
    while True:
        letters = chr(ord("A") + current % 26) + letters
        current = current // 26 - 1
        if current < 0:
            break
    return letters
