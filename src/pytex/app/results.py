"""The shape every result has, so export never has to be written twice.

A result object in this application is not free-form JSON. It carries a title, a
prose summary, an optional **table**, structured data, and the inputs it came
from. The uniformity is what makes the export surface possible: any operation
that fills in a table is exportable to CSV and XLSX by generic code, and any
result at all is exportable to JSON with its provenance attached, without the
operation knowing that export exists.

The prose summary is the application's face of the explainable-results doctrine
in ``AGENTS.md``: a result that cannot say in words what it means, with its
conventions stated, is not finished.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

__all__ = ["APP_RESULT_SCHEMA", "AppResult", "Column", "ResultTable"]

#: Schema identifier of the result payload.
APP_RESULT_SCHEMA = "pytex.app_result/1"


@dataclass(frozen=True)
class Column:
    """One column of a result table.

    Attributes
    ----------
    key : str
        Key in each row mapping.
    label : str
        Header shown in the UI and written to CSV/XLSX.
    units : str, optional
        Rendered in the header, e.g. ``"°"`` or ``"Å"``.
    help_text : str, optional
        Explains the column when its meaning is not obvious from the header.
    numeric : bool
        Right-aligns the column and formats it with ``digits``.
    digits : int, optional
        Decimal places for display. Export keeps full precision regardless.
    """

    key: str
    label: str
    units: str | None = None
    help_text: str | None = None
    numeric: bool = False
    digits: int | None = None

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this column."""

        payload: dict[str, Any] = {"key": self.key, "label": self.label, "numeric": self.numeric}
        if self.units is not None:
            payload["units"] = self.units
        if self.help_text is not None:
            payload["help"] = self.help_text
        if self.digits is not None:
            payload["digits"] = self.digits
        return payload


@dataclass(frozen=True)
class ResultTable:
    """Tabular numbers behind a result, one row per reported entity."""

    columns: tuple[Column, ...]
    rows: tuple[Mapping[str, Any], ...] = ()
    caption: str | None = None

    def __post_init__(self) -> None:
        keys = {column.key for column in self.columns}
        for index, row in enumerate(self.rows):
            missing = keys - set(row)
            if missing:
                raise ValueError(
                    f"Row {index} of the result table is missing column(s): "
                    f"{', '.join(sorted(missing))}."
                )

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this table."""

        payload: dict[str, Any] = {
            "columns": [column.to_json() for column in self.columns],
            "rows": [dict(row) for row in self.rows],
        }
        if self.caption is not None:
            payload["caption"] = self.caption
        return payload


@dataclass(frozen=True)
class AppResult:
    """A complete answer: what it is, what it says, and the numbers behind it.

    Attributes
    ----------
    title : str
        Heading for the result panel.
    summary : str
        Prose the user can read and, ideally, paste into a notebook. State the
        conventions the numbers depend on.
    table : ResultTable, optional
        The re-plottable numbers. Fill this in whenever the result *is* a set of
        rows; the export surface needs nothing else to write CSV and XLSX.
    data : mapping
        Structured extras that are not rows: matrices, scene descriptions,
        nested reports.
    inputs : mapping
        Echo of the resolved inputs, so a saved JSON export regenerates the run.
    notes : sequence of str
        Caveats worth showing beside the result — an absent space group, a
        degenerate family, a tolerance that was hit.
    citations : sequence of str
        Sources for the science behind this particular answer.
    """

    title: str
    summary: str
    table: ResultTable | None = None
    data: Mapping[str, Any] = field(default_factory=dict)
    inputs: Mapping[str, Any] = field(default_factory=dict)
    notes: Sequence[str] = ()
    citations: Sequence[str] = ()

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this result."""

        payload: dict[str, Any] = {
            "schema": APP_RESULT_SCHEMA,
            "title": self.title,
            "summary": self.summary,
            "data": dict(self.data),
            "inputs": dict(self.inputs),
            "notes": list(self.notes),
            "citations": list(self.citations),
        }
        if self.table is not None:
            payload["table"] = self.table.to_json()
        return payload

    def describe(self) -> str:
        """Return the result as prose, including its caveats.

        The explainable-results surface: what a user copies into a lab notebook
        when they want the answer *and* the conditions it holds under.
        """

        lines = [self.title, "", self.summary]
        if self.notes:
            lines.extend(["", *(f"Note: {note}" for note in self.notes)])
        if self.citations:
            lines.extend(["", "Sources: " + "; ".join(self.citations)])
        return "\n".join(lines)
