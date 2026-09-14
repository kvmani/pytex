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

__all__ = [
    "APP_RESULT_SCHEMA",
    "STAGE_STATUSES",
    "AppResult",
    "Column",
    "ResultMetric",
    "ResultStage",
    "ResultTable",
]

#: How a stage judged its own outcome. ``"ok"`` means the stage did what it is
#: for; ``"warning"`` means it completed but left something a reader must check
#: before believing what follows; ``"info"`` is a stage that records a setting
#: or an input rather than a computation.
STAGE_STATUSES = ("ok", "warning", "info")

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
class ResultMetric:
    """One named number a stage produced, with its units and what it means.

    Attributes
    ----------
    label : str
        What the number is, in words a reader already has.
    value : Any
        The value, kept at full precision; the renderer rounds for display.
    units : str, optional
        Units of ``value``.
    help_text : str, optional
        How to judge the number: what a good value looks like, and what a bad
        one says about the analysis.
    """

    label: str
    value: Any
    units: str | None = None
    help_text: str | None = None

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this metric."""

        payload: dict[str, Any] = {"label": self.label, "value": self.value}
        if self.units is not None:
            payload["units"] = self.units
        if self.help_text is not None:
            payload["help"] = self.help_text
        return payload


@dataclass(frozen=True)
class ResultStage:
    """One intermediate step on the way to a result, reported in its own right.

    Purpose
    -------
    A final number is only as believable as the steps behind it. An analysis
    that detects peaks, indexes them, and then fits a cell can go wrong at any
    of the three, and a reader who sees only the cell cannot tell which. A
    stage carries what that step produced — its own numbers, its own table —
    and a sentence saying how to read them, so the chain of evidence is on the
    page rather than inside the program.

    Attributes
    ----------
    key : str
        Stable identifier, unique within a result, so tests and exports can
        refer to a stage without depending on its title.
    title : str
        Heading, normally numbered in the order the stages ran.
    summary : str
        What the stage found, as prose.
    metrics : tuple of ResultMetric
        The stage's headline numbers.
    table : ResultTable, optional
        The stage's rows, for instance every detected peak.
    explanation : str, optional
        How to read the stage: what each number means and what to look for.
    status : str
        One of :data:`STAGE_STATUSES`.
    """

    key: str
    title: str
    summary: str
    metrics: tuple[ResultMetric, ...] = ()
    table: ResultTable | None = None
    explanation: str | None = None
    status: str = "ok"

    def __post_init__(self) -> None:
        if not self.key.strip() or not self.title.strip():
            raise ValueError("A result stage needs a non-empty key and title.")
        if self.status not in STAGE_STATUSES:
            raise ValueError(f"A result stage status must be one of {STAGE_STATUSES}.")
        object.__setattr__(self, "metrics", tuple(self.metrics))

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this stage."""

        payload: dict[str, Any] = {
            "key": self.key,
            "title": self.title,
            "summary": self.summary,
            "status": self.status,
            "metrics": [metric.to_json() for metric in self.metrics],
        }
        if self.table is not None:
            payload["table"] = self.table.to_json()
        if self.explanation is not None:
            payload["explanation"] = self.explanation
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
    stages : sequence of ResultStage
        The intermediate steps behind the answer, in the order they ran. Empty
        for a result that is a single computation. Serialized only when
        present, so a result without stages keeps its existing wire form.
    """

    title: str
    summary: str
    table: ResultTable | None = None
    data: Mapping[str, Any] = field(default_factory=dict)
    inputs: Mapping[str, Any] = field(default_factory=dict)
    notes: Sequence[str] = ()
    citations: Sequence[str] = ()
    stages: Sequence[ResultStage] = ()

    def __post_init__(self) -> None:
        keys = [stage.key for stage in self.stages]
        if len(keys) != len(set(keys)):
            raise ValueError("Result stage keys must be unique within one result.")

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
        if self.stages:
            payload["stages"] = [stage.to_json() for stage in self.stages]
        return payload

    def describe(self) -> str:
        """Return the result as prose, including its caveats.

        The explainable-results surface: what a user copies into a lab notebook
        when they want the answer *and* the conditions it holds under.
        """

        lines = [self.title, "", self.summary]
        if self.stages:
            lines.extend(["", "How the result was reached:"])
            lines.extend(f"- {stage.title}: {stage.summary}" for stage in self.stages)
        if self.notes:
            lines.extend(["", *(f"Note: {note}" for note in self.notes)])
        if self.citations:
            lines.extend(["", "Sources: " + "; ".join(self.citations)])
        return "\n".join(lines)
