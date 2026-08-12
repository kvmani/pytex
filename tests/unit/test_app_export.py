"""Export, which is written once and must therefore work for every result.

The `.xlsx` writer gets the most attention here because it is the one format
this repository implements itself. The tests read the produced workbook back
through :mod:`zipfile` and :mod:`xml.etree`, checking the parts a reader
requires and that numbers arrive as numbers rather than as text — the failure
that makes a spreadsheet look right and sort wrong.
"""

from __future__ import annotations

import csv
import io
import json
import zipfile
from typing import Any
from xml.etree import ElementTree

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.export import (
    EXPORT_FORMATS,
    export_result,
    result_to_csv,
    result_to_json,
    result_to_xlsx,
    write_xlsx,
)

_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


@pytest.fixture(scope="module")
def result() -> dict[str, Any]:
    return REGISTRY.call("calc.d_spacings", {"phase": {"builtin": "ni_fcc"}, "max_index": 3})


def read_sheet(payload: bytes, index: int = 1) -> list[list[str]]:
    """Read a sheet back as rows of strings, whatever the cell type."""

    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        tree = ElementTree.fromstring(archive.read(f"xl/worksheets/sheet{index}.xml"))
    rows: list[list[str]] = []
    for row in tree.iter(f"{_MAIN}row"):
        values: list[str] = []
        for cell in row.iter(f"{_MAIN}c"):
            inline = cell.find(f"{_MAIN}is/{_MAIN}t")
            numeric = cell.find(f"{_MAIN}v")
            if inline is not None:
                values.append(inline.text or "")
            elif numeric is not None:
                values.append(numeric.text or "")
            else:
                values.append("")
        rows.append(values)
    return rows


class TestCsv:
    def test_header_carries_the_units(self, result: dict[str, Any]) -> None:
        text = result_to_csv(result).decode("utf-8-sig")
        header = next(csv.reader(io.StringIO(text)))
        assert "d (Å)" in header

    def test_every_row_is_written(self, result: dict[str, Any]) -> None:
        text = result_to_csv(result).decode("utf-8-sig")
        rows = list(csv.reader(io.StringIO(text)))
        assert len(rows) == len(result["table"]["rows"]) + 1

    def test_full_precision_is_kept(self, result: dict[str, Any]) -> None:
        text = result_to_csv(result).decode("utf-8-sig")
        rows = list(csv.reader(io.StringIO(text)))
        header = rows[0]
        column = header.index("d (Å)")
        first = float(rows[1][column])
        assert first == pytest.approx(float(result["table"]["rows"][0]["d_angstrom"]), abs=0.0)

    def test_a_field_containing_a_comma_is_quoted(self) -> None:
        payload = {
            "title": "t",
            "summary": "s",
            "table": {
                "columns": [{"key": "a", "label": "A"}],
                "rows": [{"a": "one, two"}],
            },
        }
        text = result_to_csv(payload).decode("utf-8-sig")
        assert '"one, two"' in text


class TestJson:
    def test_the_inputs_travel_with_the_result(self, result: dict[str, Any]) -> None:
        payload = json.loads(result_to_json(result))
        assert payload["inputs"]["phase"]["name"] == "Nickel (fcc)"
        assert payload["schema"] == "pytex.app_result/1"

    def test_an_export_timestamp_is_added(self, result: dict[str, Any]) -> None:
        payload = json.loads(result_to_json(result))
        assert payload["exported_utc"].endswith("+00:00")

    def test_a_result_without_a_table_still_exports_as_json(self) -> None:
        payload = {"title": "t", "summary": "s", "data": {"x": 1}}
        assert json.loads(result_to_json(payload))["data"] == {"x": 1}


class TestXlsx:
    def test_the_workbook_holds_the_parts_a_reader_requires(self, result: dict[str, Any]) -> None:
        with zipfile.ZipFile(io.BytesIO(result_to_xlsx(result))) as archive:
            names = set(archive.namelist())
            assert {
                "[Content_Types].xml",
                "_rels/.rels",
                "xl/workbook.xml",
                "xl/_rels/workbook.xml.rels",
                "xl/worksheets/sheet1.xml",
            } <= names
            assert archive.testzip() is None

    def test_the_data_sheet_matches_the_table(self, result: dict[str, Any]) -> None:
        rows = read_sheet(result_to_xlsx(result))
        assert len(rows) == len(result["table"]["rows"]) + 1
        assert rows[0][0] == "Family"

    def test_numbers_are_written_as_numbers(self, result: dict[str, Any]) -> None:
        with zipfile.ZipFile(io.BytesIO(result_to_xlsx(result))) as archive:
            sheet = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
        # A numeric cell has a bare <v>; a text cell is t="inlineStr". Getting
        # this wrong produces a spreadsheet that looks right and sorts wrong.
        assert "<v>" in sheet
        first_row = read_sheet(result_to_xlsx(result))[1]
        assert float(first_row[4]) == pytest.approx(float(result["table"]["rows"][0]["d_angstrom"]))

    def test_the_provenance_sheet_records_the_inputs(self, result: dict[str, Any]) -> None:
        rows = read_sheet(result_to_xlsx(result), index=2)
        fields = {row[0] for row in rows if row}
        assert "Result" in fields
        assert any(field.startswith("Input: ") for field in fields)
        assert any(field == "Source" for field in fields)

    def test_a_second_sheet_is_declared_in_the_workbook(self, result: dict[str, Any]) -> None:
        with zipfile.ZipFile(io.BytesIO(result_to_xlsx(result))) as archive:
            workbook = archive.read("xl/workbook.xml").decode("utf-8")
        assert 'name="Data"' in workbook
        assert 'name="Provenance"' in workbook

    def test_special_characters_are_escaped(self) -> None:
        payload = write_xlsx({"S": {"headers": ["a<b & c>"], "rows": [['"quoted"']]}})
        rows = read_sheet(payload)
        assert rows[0][0] == "a<b & c>"
        assert rows[1][0] == '"quoted"'

    def test_a_forbidden_sheet_name_is_repaired(self) -> None:
        payload = write_xlsx({"a/b:c*d?e[f]g" * 4: {"headers": ["x"], "rows": [[1]]}})
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            workbook = archive.read("xl/workbook.xml").decode("utf-8")
        # Excel refuses to open a workbook whose sheet name is too long or holds
        # a forbidden character, so both are repaired rather than passed through.
        name = workbook.split('name="')[1].split('"')[0]
        assert len(name) <= 31
        assert not set(name) & set(r"[]:*?/\\")

    def test_booleans_are_readable_words(self) -> None:
        rows = read_sheet(write_xlsx({"S": {"headers": ["ok"], "rows": [[True], [False]]}}))
        assert rows[1][0] == "yes"
        assert rows[2][0] == "no"

    def test_an_empty_workbook_is_refused(self) -> None:
        with pytest.raises(ValueError, match="at least one sheet"):
            write_xlsx({})

    def test_columns_past_z_are_addressed_correctly(self) -> None:
        headers = [f"c{index}" for index in range(30)]
        payload = write_xlsx({"S": {"headers": headers, "rows": [list(range(30))]}})
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            sheet = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
        assert 'r="AA1"' in sheet
        assert 'r="AD1"' in sheet


class TestExportResult:
    @pytest.mark.parametrize("fmt", sorted(EXPORT_FORMATS))
    def test_every_declared_format_produces_bytes(self, result: dict[str, Any], fmt: str) -> None:
        payload, mime, filename = export_result(result, fmt=fmt)
        assert payload
        assert mime
        assert filename.endswith(EXPORT_FORMATS[fmt]["extension"])

    def test_the_filename_comes_from_the_title(self, result: dict[str, Any]) -> None:
        _, _, filename = export_result(result, fmt="csv")
        assert filename == "d-spacings-of-nickel-fcc.csv"

    def test_an_unknown_format_lists_the_known_ones(self, result: dict[str, Any]) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            export_result(result, fmt="pdf")
        assert "csv" in (excinfo.value.hint or "")

    def test_a_tableless_result_says_to_use_json(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            export_result({"title": "t", "summary": "s"}, fmt="csv")
        assert "JSON" in (excinfo.value.hint or "")


class TestExportRoute:
    def test_the_server_serves_a_download(self, result: dict[str, Any]) -> None:
        import json as json_module
        import urllib.request

        from pytex.app.server import create_server

        server = create_server("127.0.0.1", 0)
        server.serve_in_background()
        try:
            request = urllib.request.Request(
                f"{server.url}/api/export",
                data=json_module.dumps({"result": result, "format": "xlsx"}).encode(),
                headers={"Content-Type": "application/json"},
            )
            with urllib.request.urlopen(request) as response:
                payload = response.read()
                assert response.headers["Content-Disposition"].startswith("attachment;")
                assert response.headers["Content-Type"].endswith("spreadsheetml.sheet")
            assert zipfile.ZipFile(io.BytesIO(payload)).testzip() is None
        finally:
            server.shutdown()
            server.server_close()
