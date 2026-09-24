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

import base64
import csv
import html
import io
import json
import zipfile
from collections.abc import Callable, Mapping, Sequence
from datetime import UTC, datetime
from typing import Any
from xml.sax.saxutils import escape

from pytex.app.contracts import to_jsonable
from pytex.app.errors import InvalidInputError
from pytex.app.results import REPORT_SECTIONS

__all__ = [
    "EXPORT_FORMATS",
    "export_result",
    "result_to_bundle",
    "result_to_csv",
    "result_to_html",
    "result_to_json",
    "result_to_markdown",
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
    "md": {
        "label": "Report",
        "mime": "text/markdown; charset=utf-8",
        "extension": "md",
        "description": "A readable report: what was computed, from what, with the sources.",
    },
    "html": {
        "label": "Printable report",
        "mime": "text/html; charset=utf-8",
        "extension": "html",
        "description": (
            "The report as one self-contained web page with every figure: open it in any "
            "browser, or print it to PDF."
        ),
    },
    "zip": {
        "label": "Report + figures",
        "mime": "application/zip",
        "extension": "zip",
        "description": (
            "The report (Markdown and printable HTML) with every figure as a separate SVG "
            "file, and the complete result as JSON."
        ),
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

    sheets: dict[str, Mapping[str, Any]] = {
        "Data": {"headers": headers, "rows": rows},
        "Provenance": {"headers": ["Field", "Value"], "rows": provenance},
    }
    # Every intermediate stage is written too: one sheet of stage metrics, and
    # one sheet per stage table. A workbook holding only the final table would
    # drop exactly the evidence the stages exist to make reviewable.
    stages = list(result.get("stages") or ())
    if stages:
        metric_rows: list[list[Any]] = []
        for stage in stages:
            metric_rows.append([stage.get("title", ""), "Summary", stage.get("summary", ""), ""])
            for metric in stage.get("metrics") or ():
                metric_rows.append(
                    [
                        stage.get("title", ""),
                        metric.get("label", ""),
                        _cell(metric.get("value")),
                        metric.get("units", ""),
                    ]
                )
        sheets["Stages"] = {"headers": ["Stage", "Quantity", "Value", "Units"], "rows": metric_rows}
        for number, stage in enumerate(stages, start=1):
            stage_table = stage.get("table") or {}
            stage_columns = stage_table.get("columns") or []
            if not stage_columns:
                continue
            sheets[f"Stage {number} {stage.get('key', '')}"] = {
                "headers": [
                    f"{column['label']} ({column['units']})"
                    if column.get("units")
                    else column["label"]
                    for column in stage_columns
                ],
                "rows": [
                    [_cell(row.get(column["key"])) for column in stage_columns]
                    for row in stage_table.get("rows") or ()
                ],
            }
    return write_xlsx(sheets)


def result_to_markdown(
    result: Mapping[str, Any],
    *,
    figure_link: Callable[[Mapping[str, Any]], str] | None = None,
) -> bytes:
    """Write the result as a report a person can read without a spreadsheet.

    Purpose
    -------
    The other formats are for machines or for grids. CSV is a table with no
    account of where it came from; JSON is complete and unreadable; a workbook is
    a table plus a sheet of key-value pairs. None of them is the thing to paste
    into a notebook entry or attach to an email, which is a page saying *what was
    computed, from what, how far to trust it, and on whose authority*.

    That is what this writes, in the order a reader needs it: the answer and its
    reliability, the warnings, the figures, then the stages grouped as evidence,
    diagnostics, method and audit details (see
    :data:`~pytex.app.results.REPORT_SECTIONS`), then the exact inputs and the
    citations. Markdown because it is readable as plain text, renders
    everywhere, and survives being pasted into anything.

    Parameters
    ----------
    result : mapping
        The result payload.
    figure_link : callable, optional
        Maps a figure payload to the image target written into the report. By
        default each figure is embedded as a base64 SVG data URL, so the single
        ``.md`` file is self-contained; the zip bundle passes a relative path.

    A result with no table still exports - the prose and the provenance are the
    point - which is why this does not go through ``_table_of``.
    """

    link = figure_link or _embedded_figure
    lines: list[str] = [f"# {result.get('title', 'PyTex result')}", ""]
    summary = str(result.get("summary") or "").strip()
    if summary:
        lines += [summary, ""]

    highlights = list(result.get("highlights") or ())
    if highlights:
        lines += ["## Result and reliability", ""]
        lines += _metric_table(highlights)
    warnings = [str(item) for item in (result.get("warnings") or ())]
    if warnings:
        lines += ["## Warnings", ""]
        lines += [f"- **Warning.** {warning}" for warning in warnings]
        lines.append("")

    for figure in result.get("figures") or ():
        lines += _figure_block(figure, link)

    stages = list(result.get("stages") or ())
    sectioned = any(stage.get("section") for stage in stages)
    table_lines = _data_table(result)

    if sectioned:
        for section, heading in REPORT_SECTIONS.items():
            members = [stage for stage in stages if stage.get("section") == section]
            if section == "audit" and table_lines:
                members_lines = table_lines
                table_lines = []
            else:
                members_lines = []
            if not members and not members_lines:
                continue
            if section == "result" and highlights:
                heading = "Result in detail"
            lines += [f"## {heading}", ""]
            for stage in members:
                lines += _stage_block(stage, link)
            lines += members_lines
        loose = [stage for stage in stages if not stage.get("section")]
        if loose:
            lines += ["## Further stages", ""]
            for stage in loose:
                lines += _stage_block(stage, link)
        lines += table_lines
    else:
        lines += table_lines
        if stages:
            lines += ["## How the result was reached", ""]
            for stage in stages:
                lines += _stage_block(stage, link)

    notes = [str(note) for note in (result.get("notes") or ())]
    if notes:
        lines += ["## Notes", ""]
        lines += [f"- {note}" for note in notes]
        lines.append("")

    inputs = result.get("inputs") or {}
    if inputs:
        lines += ["## Inputs", "", "| Field | Value |", "| --- | --- |"]
        for key, value in sorted(inputs.items()):
            rendered = value if isinstance(value, str | int | float | bool) else json.dumps(value)
            lines.append(f"| {key} | {_markdown_cell(rendered)} |")
        lines.append("")

    citations = [str(item) for item in (result.get("citations") or ())]
    if citations:
        lines += ["## Sources", ""]
        lines += [f"- {citation}" for citation in citations]
        lines.append("")

    lines += [
        "---",
        "",
        f"Produced by PyTex, exported {datetime.now(UTC).isoformat(timespec='seconds')}.",
        "",
    ]
    return "\n".join(lines).encode("utf-8")


_HTML_STYLE = """
:root { color-scheme: light; }
body { font-family: "Source Serif 4", Georgia, "Times New Roman", serif; color: #111827;
  background: #ffffff; max-width: 60rem; margin: 2.5rem auto; padding: 0 1.25rem;
  line-height: 1.5; font-size: 11pt; }
h1 { font-size: 1.7rem; line-height: 1.25; margin: 0 0 0.4rem; }
h2 { font-size: 1.3rem; border-bottom: 2px solid #111827; padding-bottom: 0.2rem;
  margin-top: 2.2rem; }
h3 { font-size: 1.1rem; margin-top: 1.6rem; }
h1, h2, h3, figcaption, table, .meta { font-family: "Inter", "Segoe UI", Arial, sans-serif; }
.meta { color: #4b5563; font-size: 0.85rem; margin-bottom: 1.2rem; }
.summary { font-size: 1.05rem; }
.warning { border-left: 4px solid #dc2626; background: #fef2f2; padding: 0.5rem 0.8rem;
  margin: 0.5rem 0; }
.check { color: #b45309; font-size: 0.85rem; font-weight: 600; }
table { border-collapse: collapse; margin: 0.8rem 0; font-size: 0.8rem; width: 100%; }
caption { text-align: left; font-style: italic; color: #374151; padding-bottom: 0.3rem; }
th, td { border-bottom: 1px solid #d1d5db; padding: 0.25rem 0.5rem; text-align: left;
  vertical-align: top; }
th { border-bottom: 2px solid #6b7280; }
td.number { text-align: right; font-variant-numeric: tabular-nums; }
figure { margin: 1.2rem 0; page-break-inside: avoid; break-inside: avoid; }
figure img { max-width: 100%; height: auto; display: block; margin: 0 auto; }
figcaption { font-size: 0.85rem; color: #374151; margin-top: 0.4rem; }
figcaption strong { color: #111827; }
.explanation { white-space: pre-line; background: #f9fafb; border: 1px solid #e5e7eb;
  padding: 0.6rem 0.9rem; font-size: 0.92rem; }
.explanation::before { content: "How to read this. "; font-weight: 600; }
.sources li { font-size: 0.9rem; }
@media print {
  body { margin: 0; max-width: none; font-size: 10pt; }
  h2 { page-break-after: avoid; break-after: avoid; }
  h3 { page-break-after: avoid; break-after: avoid; }
  a { color: inherit; text-decoration: none; }
}
"""


def result_to_html(result: Mapping[str, Any]) -> bytes:
    """Write the result as one self-contained, printable web page.

    Purpose
    -------
    The Markdown report is readable as text and renders in most viewers, but
    it is not a document a reader can simply open, read with its figures in
    place, and print. This is the same report, in the same order (see
    :func:`result_to_markdown`), as a single HTML file: every figure embedded as
    an SVG data URL, the tables as tables, a stylesheet that paginates sensibly
    when printed to PDF, and no external resource of any kind, so it reads the
    same offline, attached to an e-mail, or archived with a data set.

    Nothing is recomputed: every number and every word comes from the result
    payload, exactly as the other formats write it.
    """

    def text(value: Any) -> str:
        return html.escape(str(value if value is not None else ""))

    def figure_block(figure: Mapping[str, Any]) -> list[str]:
        title = text(figure.get("title") or figure.get("key") or "Figure")
        parts = [
            "<figure>",
            f'<img src="{_embedded_figure(figure)}" alt="{title}">',
            f"<figcaption><strong>{title}.</strong> {text(figure.get('caption') or '')}",
        ]
        interpretation = str(figure.get("interpretation") or "").strip()
        if interpretation:
            parts.append(f"<br><em>What it shows.</em> {text(interpretation)}")
        parts.append("</figcaption></figure>")
        return parts

    def metric_table(metrics: Sequence[Mapping[str, Any]]) -> list[str]:
        rows = [
            "<tr>"
            f"<td>{text(metric.get('label'))}</td>"
            f'<td class="number">{text(_format_metric(metric.get("value")))}</td>'
            f"<td>{text(metric.get('units') or '')}</td>"
            "</tr>"
            for metric in metrics
        ]
        return [
            "<table><thead><tr><th>Quantity</th><th>Value</th><th>Units</th></tr></thead>",
            "<tbody>",
            *rows,
            "</tbody></table>",
        ]

    def rows_table(
        columns: Sequence[Mapping[str, Any]],
        rows: Sequence[Mapping[str, Any]],
        caption: str = "",
    ) -> list[str]:
        head = "".join(
            f"<th>{text(column['label'])}"
            + (f" ({text(column['units'])})" if column.get("units") else "")
            + "</th>"
            for column in columns
        )
        body = []
        for row in rows:
            cells = []
            for column in columns:
                value = row.get(column["key"])
                numeric = bool(column.get("numeric")) and isinstance(value, int | float)
                if numeric and isinstance(value, float) and column.get("digits") is not None:
                    shown = f"{value:.{int(column['digits'])}f}"
                else:
                    shown = _markdown_cell(value).replace("\\|", "|")
                cells.append(
                    f'<td class="number">{text(shown)}</td>'
                    if numeric
                    else f"<td>{text(shown)}</td>"
                )
            body.append("<tr>" + "".join(cells) + "</tr>")
        return [
            "<table>",
            f"<caption>{text(caption)}</caption>" if caption else "",
            f"<thead><tr>{head}</tr></thead>",
            "<tbody>",
            *body,
            "</tbody></table>",
        ]

    def stage_block(stage: Mapping[str, Any]) -> list[str]:
        marker = (
            ' <span class="check">(check this)</span>'
            if str(stage.get("status") or "ok") == "warning"
            else ""
        )
        parts = [f"<h3>{text(stage.get('title', ''))}{marker}</h3>"]
        summary = str(stage.get("summary") or "").strip()
        if summary:
            parts.append(f"<p>{text(summary)}</p>")
        metrics = list(stage.get("metrics") or ())
        if metrics:
            parts += metric_table(metrics)
        for figure in stage.get("figures") or ():
            parts += figure_block(figure)
        table = stage.get("table") or {}
        if table.get("columns") and table.get("rows"):
            parts += rows_table(table["columns"], table["rows"], str(table.get("caption") or ""))
        explanation = str(stage.get("explanation") or "").strip()
        if explanation:
            parts.append(f'<div class="explanation">{text(explanation)}</div>')
        return parts

    title = text(result.get("title", "PyTex result"))
    exported = datetime.now(UTC).isoformat(timespec="seconds")
    body: list[str] = [
        f"<h1>{title}</h1>",
        f'<div class="meta">Produced by PyTex · exported {exported}</div>',
    ]
    summary = str(result.get("summary") or "").strip()
    if summary:
        body.append(f'<p class="summary">{text(summary)}</p>')
    highlights = list(result.get("highlights") or ())
    if highlights:
        body += ["<h2>Result and reliability</h2>", *metric_table(highlights)]
    warnings = [str(item) for item in (result.get("warnings") or ())]
    if warnings:
        body.append("<h2>Warnings</h2>")
        body += [f'<p class="warning"><strong>Warning.</strong> {text(w)}</p>' for w in warnings]
    for figure in result.get("figures") or ():
        body += figure_block(figure)

    stages = list(result.get("stages") or ())
    table = result.get("table") or {}
    data_table = (
        [
            "<h3>Data</h3>",
            *rows_table(table["columns"], table["rows"], str(table.get("caption") or "")),
        ]
        if table.get("columns") and table.get("rows")
        else []
    )
    if any(stage.get("section") for stage in stages):
        for section, heading in REPORT_SECTIONS.items():
            members = [stage for stage in stages if stage.get("section") == section]
            extra = data_table if section == "audit" else []
            if section == "audit":
                data_table = []
            if not members and not extra:
                continue
            if section == "result" and highlights:
                heading = "Result in detail"
            body.append(f"<h2>{text(heading)}</h2>")
            for stage in members:
                body += stage_block(stage)
            body += extra
        loose = [stage for stage in stages if not stage.get("section")]
        if loose:
            body.append("<h2>Further stages</h2>")
            for stage in loose:
                body += stage_block(stage)
        body += data_table
    else:
        body += data_table
        if stages:
            body.append("<h2>How the result was reached</h2>")
            for stage in stages:
                body += stage_block(stage)

    notes = [str(note) for note in (result.get("notes") or ())]
    if notes:
        body += ["<h2>Notes</h2>", "<ul>", *(f"<li>{text(note)}</li>" for note in notes), "</ul>"]
    inputs = result.get("inputs") or {}
    if inputs:
        rows = []
        for key, value in sorted(inputs.items()):
            rendered = value if isinstance(value, str | int | float | bool) else json.dumps(value)
            rows.append(f"<tr><td>{text(key)}</td><td>{text(rendered)}</td></tr>")
        body += [
            "<h2>Inputs</h2>",
            "<table><thead><tr><th>Field</th><th>Value</th></tr></thead><tbody>",
            *rows,
            "</tbody></table>",
        ]
    citations = [str(item) for item in (result.get("citations") or ())]
    if citations:
        body += [
            "<h2>Sources</h2>",
            '<ol class="sources">',
            *(f"<li>{text(citation)}</li>" for citation in citations),
            "</ol>",
        ]
    document = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{title}</title>",
        f"<style>{_HTML_STYLE}</style>",
        "</head>",
        "<body>",
        *(line for line in body if line),
        "</body>",
        "</html>",
        "",
    ]
    return "\n".join(document).encode("utf-8")


def result_to_bundle(result: Mapping[str, Any]) -> bytes:
    """Write a zip holding the report, every figure as its own file, and the JSON.

    Purpose
    -------
    The single Markdown file embeds its figures, which keeps it self-contained
    but makes a figure awkward to lift into a manuscript. The bundle is the same
    report with each figure written beside it as ``figures/<key>.svg`` and
    linked by relative path, plus the complete result as JSON, so every number
    and every picture in the report can be traced to the file that holds it.
    """

    figures = _all_figures(result)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            "report.md",
            result_to_markdown(
                result, figure_link=lambda figure: f"figures/{_slug(str(figure['key']))}.svg"
            ),
        )
        for figure in figures:
            archive.writestr(
                f"figures/{_slug(str(figure['key']))}.svg", str(figure["svg"]).encode("utf-8")
            )
        archive.writestr("report.html", result_to_html(result))
        archive.writestr("result.json", result_to_json(result))
    return buffer.getvalue()


def _all_figures(result: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    figures = list(result.get("figures") or ())
    for stage in result.get("stages") or ():
        figures.extend(stage.get("figures") or ())
    return figures


def _embedded_figure(figure: Mapping[str, Any]) -> str:
    encoded = base64.b64encode(str(figure.get("svg", "")).encode("utf-8")).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def _figure_block(figure: Mapping[str, Any], link: Callable[[Mapping[str, Any]], str]) -> list[str]:
    title = str(figure.get("title") or figure.get("key") or "Figure")
    lines = [f"#### Figure: {title}", "", f"![{title}]({link(figure)})", ""]
    caption = str(figure.get("caption") or "").strip()
    if caption:
        lines += [f"*{caption}*", ""]
    interpretation = str(figure.get("interpretation") or "").strip()
    if interpretation:
        lines += [f"**What it shows.** {interpretation}", ""]
    return lines


def _metric_table(metrics: Sequence[Mapping[str, Any]]) -> list[str]:
    lines = ["| Quantity | Value | Units |", "| --- | --- | --- |"]
    for metric in metrics:
        lines.append(
            f"| {_markdown_cell(metric.get('label'))} | "
            f"{_markdown_cell(_format_metric(metric.get('value')))} | "
            f"{_markdown_cell(metric.get('units'))} |"
        )
    lines.append("")
    return lines


def _format_metric(value: Any) -> Any:
    """Round a float for reading; full precision stays in the JSON and XLSX exports."""

    if isinstance(value, float) and not isinstance(value, bool):
        if value != value:  # NaN
            return "n/a"
        return f"{value:.6g}"
    return value


def _rows_table(
    columns: Sequence[Mapping[str, Any]], rows: Sequence[Mapping[str, Any]]
) -> list[str]:
    headers = [
        f"{column['label']} / {column['units']}" if column.get("units") else column["label"]
        for column in columns
    ]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append(
            "| " + " | ".join(_markdown_cell(row.get(column["key"])) for column in columns) + " |"
        )
    lines.append("")
    return lines


def _data_table(result: Mapping[str, Any]) -> list[str]:
    table = result.get("table") or {}
    columns = table.get("columns") or []
    rows = table.get("rows") or []
    if not (columns and rows):
        return []
    lines = ["### Data", ""]
    caption = str(table.get("caption") or "").strip()
    if caption:
        lines += [caption, ""]
    return lines + _rows_table(columns, rows)


def _stage_block(stage: Mapping[str, Any], link: Callable[[Mapping[str, Any]], str]) -> list[str]:
    status = str(stage.get("status") or "ok")
    marker = " (check this)" if status == "warning" else ""
    lines = [f"### {stage.get('title', '')}{marker}", ""]
    stage_summary = str(stage.get("summary") or "").strip()
    if stage_summary:
        lines += [stage_summary, ""]
    metrics = list(stage.get("metrics") or ())
    if metrics:
        lines += _metric_table(metrics)
    for figure in stage.get("figures") or ():
        lines += _figure_block(figure, link)
    stage_table = stage.get("table") or {}
    stage_columns = stage_table.get("columns") or []
    stage_rows = stage_table.get("rows") or []
    if stage_columns and stage_rows:
        stage_caption = str(stage_table.get("caption") or "").strip()
        if stage_caption:
            lines += [stage_caption, ""]
        lines += _rows_table(stage_columns, stage_rows)
    explanation = str(stage.get("explanation") or "").strip()
    if explanation:
        lines += [f"*How to read this.* {explanation}", ""]
    return lines


def _markdown_cell(value: Any) -> str:
    """Render one cell, escaping the pipe that would otherwise split the row."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value).replace("|", r"\|").replace("\n", " ")


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
    writers = {
        "csv": result_to_csv,
        "xlsx": result_to_xlsx,
        "json": result_to_json,
        "md": result_to_markdown,
        "html": result_to_html,
        "zip": result_to_bundle,
    }
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
