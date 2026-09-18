"""The shared figure contract: drawn once, shown, exported and reported identically.

A figure in the workbench is part of a result, not decoration: it travels in
the payload, is embedded in the Markdown report and written as its own file in
the report bundle. These tests hold the pieces that make that true — the
figure object, the drawing helpers, the report writer, the bundle — to what
they promise.
"""

from __future__ import annotations

import base64
import io
import json
import re
import zipfile
from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.app.export import EXPORT_FORMATS, export_result, result_to_bundle, result_to_markdown
from pytex.app.figures import (
    draw_correlation,
    draw_image,
    draw_normalized_residuals,
    draw_observed_model,
    draw_residuals,
    draw_ticks,
    render_figure,
)
from pytex.app.results import (
    REPORT_SECTIONS,
    AppResult,
    ResultFigure,
    ResultMetric,
    ResultStage,
)

_SVG = "{http://www.w3.org/2000/svg}"


def _scatter(figure: object) -> None:
    axes = figure.subplots()  # type: ignore[attr-defined]
    axes.plot([0, 1, 2], [0, 1, 4], "o")


def test_render_figure_returns_self_contained_parseable_svg() -> None:
    figure = render_figure(_scatter, key="demo", title="Demo", caption="Three points.")
    assert figure.svg.startswith("<svg")
    root = ElementTree.fromstring(figure.svg)
    assert root.tag == f"{_SVG}svg"
    # Text is drawn as paths and nothing is fetched: renders the same offline.
    assert "<text" not in figure.svg
    assert not re.search(r'(?:href|src)="(?!#|data:)', figure.svg)
    assert figure.width_in == pytest.approx(6.4)


def test_render_figure_is_deterministic() -> None:
    first = render_figure(_scatter, key="demo", title="Demo", caption="c")
    second = render_figure(_scatter, key="demo", title="Demo", caption="c")
    assert first.svg == second.svg


def test_render_figure_leaves_no_pyplot_figure_open() -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    render_figure(_scatter, key="demo", title="Demo", caption="c")
    assert plt.get_fignums() == before


def test_every_drawing_helper_renders() -> None:
    x = np.linspace(20.0, 80.0, 200)
    observed = 100.0 + 50.0 * np.exp(-((x - 44.0) ** 2) / 0.5)

    def draw(figure: object) -> None:
        axes = figure.subplots(2, 3)  # type: ignore[attr-defined]
        draw_observed_model(axes[0, 0], x, observed, observed * 0.99, difference_axes=axes[1, 0])
        draw_ticks(axes[0, 0], [44.0, 51.0], labels=["(111)", "(200)"])
        draw_residuals(axes[0, 1], [40, 50], [0.1, -0.2], [0.05, 0.05], labels=["a", "b"])
        draw_normalized_residuals(
            axes[1, 1], [40, 50, 60], [0.5, -2.5, 4.0], labels=["a", "b", "c"]
        )
        draw_correlation(figure, axes[0, 2], np.array([[1.0, -0.9], [-0.9, 1.0]]), ["A", "D"])
        draw_image(figure, axes[1, 2], np.eye(4), colorbar_label="counts")

    figure = render_figure(draw, key="all", title="All", caption="Every helper.", height_in=5.0)
    ElementTree.fromstring(figure.svg)


def test_result_figure_rejects_non_svg_and_blank_key() -> None:
    with pytest.raises(ValueError, match="SVG"):
        ResultFigure(key="k", title="t", svg="<png/>", caption="c")
    with pytest.raises(ValueError, match="key"):
        ResultFigure(key=" ", title="t", svg="<svg/>", caption="c")


def test_stage_section_is_validated() -> None:
    assert list(REPORT_SECTIONS) == ["result", "evidence", "diagnostics", "method", "audit"]
    with pytest.raises(ValueError, match="section"):
        ResultStage(key="k", title="t", summary="s", section="appendix")


def test_figure_keys_must_be_unique_across_stages() -> None:
    figure = ResultFigure(key="same", title="t", svg="<svg/>", caption="c")
    stage = ResultStage(key="s", title="t", summary="s", figures=(figure,))
    with pytest.raises(ValueError, match="figure keys"):
        AppResult(title="t", summary="s", stages=(stage,), figures=(figure,))


def _sectioned_result() -> dict[str, object]:
    figure = render_figure(_scatter, key="scan", title="The scan", caption="Every point.")
    stages = (
        ResultStage(key="audit_one", title="Raw numbers", summary="Audit.", section="audit"),
        ResultStage(key="method_one", title="How", summary="Method.", section="method"),
        ResultStage(
            key="evidence_one",
            title="Peaks",
            summary="Evidence.",
            section="evidence",
            figures=(figure,),
        ),
        ResultStage(
            key="diag_one",
            title="Residuals",
            summary="Diagnostics.",
            section="diagnostics",
            status="warning",
        ),
    )
    return AppResult(
        title="Demo result",
        summary="The answer is 42.",
        highlights=(ResultMetric("Answer", 42.0, "units"),),
        warnings=("The fit is poor.",),
        stages=stages,
    ).to_json()


def test_wire_form_carries_new_fields_only_when_present() -> None:
    bare = AppResult(title="t", summary="s").to_json()
    for key in ("highlights", "warnings", "figures", "stages"):
        assert key not in bare
    payload = _sectioned_result()
    assert payload["highlights"] == [{"label": "Answer", "value": 42.0, "units": "units"}]
    assert payload["warnings"] == ["The fit is poor."]
    stages = payload["stages"]
    assert isinstance(stages, list)
    evidence = next(stage for stage in stages if stage["key"] == "evidence_one")
    assert evidence["section"] == "evidence"
    assert evidence["figures"][0]["svg"].startswith("<svg")


def test_markdown_reads_result_then_evidence_diagnostics_method_audit() -> None:
    text = result_to_markdown(_sectioned_result()).decode("utf-8")
    order = [
        text.index("## Result and reliability"),
        text.index("## Warnings"),
        text.index("## Evidence"),
        text.index("## Diagnostics"),
        text.index("## Method"),
        text.index("## Audit details"),
    ]
    assert order == sorted(order)
    assert "**Warning.** The fit is poor." in text
    assert "### Residuals (check this)" in text
    # The figure is embedded, so the single file is self-contained.
    match = re.search(r"!\[The scan\]\(data:image/svg\+xml;base64,([A-Za-z0-9+/=]+)\)", text)
    assert match is not None
    assert base64.b64decode(match.group(1)).decode("utf-8").startswith("<svg")
    assert "*Every point.*" in text


def test_bundle_holds_report_figures_and_json() -> None:
    payload = _sectioned_result()
    archive = zipfile.ZipFile(io.BytesIO(result_to_bundle(payload)))
    names = set(archive.namelist())
    assert names == {"report.md", "figures/scan.svg", "result.json"}
    report = archive.read("report.md").decode("utf-8")
    assert "![The scan](figures/scan.svg)" in report
    ElementTree.fromstring(archive.read("figures/scan.svg"))
    assert json.loads(archive.read("result.json"))["title"] == "Demo result"


def test_bundle_is_an_offered_export_format() -> None:
    assert "zip" in EXPORT_FORMATS
    data, mime, filename = export_result(_sectioned_result(), fmt="zip")
    assert mime == "application/zip"
    assert filename == "demo-result.zip"
    assert zipfile.is_zipfile(io.BytesIO(data))


def test_unsectioned_results_keep_their_report_layout() -> None:
    stage = ResultStage(key="one", title="1. First", summary="First step.")
    text = result_to_markdown(AppResult(title="t", summary="s", stages=(stage,)).to_json())
    assert "## How the result was reached" in text.decode("utf-8")
