"""The measured texture analysis, checked against answers fixed before it runs.

The demonstration figures are measured from a known model texture with 3%
counting noise, which makes each view checkable:

* the recalculated figures must reproduce the measured ones (a small RP factor),
  and the difference figure must be exactly recalculated minus measured;
* the volume fractions must name the components the figures came from, and a
  component the model does not contain must sit near or below random;
* the Euler range of the sections follows from the crystal and sample symmetry;
* axial symmetry makes every tilt ring of the symmetrized figure constant;
* a change of view reuses the inversion rather than solving it again.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

FCC: dict[str, Any] = {
    "phase": {"builtin": "ni_fcc"},
    "poles": [[1, 1, 1], [2, 0, 0], [2, 2, 0]],
    "sample_symmetry": "orthorhombic",
    "dictionary_count": 800,
}
ZR: dict[str, Any] = {
    "phase": {"builtin": "zr_hcp"},
    "poles": [[0, 0, 2], [1, 0, 0], [1, 0, 1]],
    "sample_symmetry": "orthorhombic",
    "dictionary_count": 800,
}
FIXTURE = Path("fixtures/xrdml/synthetic_random_standard.xrdml")


def analyse(base: dict[str, Any], **overrides: Any) -> dict[str, Any]:
    return REGISTRY.call("texture.analysis", {**base, **overrides})


@pytest.fixture(scope="module")
def fcc() -> dict[str, Any]:
    return analyse(FCC)


class TestOneRequestEveryView:
    def test_every_view_is_in_the_one_result(self, fcc: dict[str, Any]) -> None:
        data = fcc["data"]
        assert data["demonstration"] is True
        assert [figure["label"] for figure in data["figures"]] == ["{111}", "{200}", "{220}"]
        assert data["odf"]["sections"]
        assert data["volume_fractions"]["rows"]
        assert [stage["key"] for stage in fcc["stages"]] == [
            "measured",
            "symmetry",
            "inversion",
            "recalculated",
            "sections",
            "fractions",
        ]

    def test_each_point_carries_measured_recalculated_and_their_difference(
        self, fcc: dict[str, Any]
    ) -> None:
        for figure in fcc["data"]["figures"]:
            points = figure["points"]
            assert len(points) == 16 * 72
            measured = np.array([point["measured"] for point in points])
            recalculated = np.array([point["recalculated"] for point in points])
            difference = np.array([point["difference"] for point in points])
            assert np.allclose(difference, recalculated - measured)

    def test_the_odf_reproduces_the_figures_it_came_from(self, fcc: dict[str, Any]) -> None:
        """The inverted figures are fitted closely, and because the demonstration
        specimen really is orthorhombic, the measured ones nearly as closely."""

        # Measured at 6.1% and 6.7% with the default 800-orientation dictionary;
        # 10% is where the service itself stops calling a fit good.
        assert fcc["data"]["mean_fit_rp_percent"] < 10.0
        assert fcc["data"]["mean_rp_percent"] < 12.0
        for row in fcc["table"]["rows"]:
            assert row["rp_fit_percent"] <= row["rp_percent"] + 1e-9
            assert row["recalculated_max"] == pytest.approx(row["measured_max"], rel=0.35)

    def test_the_volume_fractions_name_the_components_the_figures_came_from(
        self, fcc: dict[str, Any]
    ) -> None:
        rows = {row["component"]: row for row in fcc["data"]["volume_fractions"]["rows"]}
        for name in ("Cube", "Goss", "Brass", "Copper", "S"):
            assert rows[name]["times_random"] > 1.5, name
        # Rotated cube is not in the model: it must not look like a component.
        assert rows["Rotated cube"]["times_random"] < rows["Cube"]["times_random"]
        for row in rows.values():
            assert row["random_percent"] == pytest.approx(
                100.0 * 24 * (np.radians(15.0) - np.sin(np.radians(15.0))) / np.pi
            )

    def test_the_standard_cubic_sections_over_the_cubic_box(self, fcc: dict[str, Any]) -> None:
        odf = fcc["data"]["odf"]
        assert [section["value_deg"] for section in odf["sections"]] == [0.0, 45.0, 65.0]
        assert odf["ranges"] == {
            "phi1_max_deg": 90.0,
            "big_phi_max_deg": 90.0,
            "phi2_max_deg": 90.0,
        }
        assert odf["horizontal_coordinate"] == "phi1"

    def test_the_prose_names_the_symmetry_and_the_fit(self, fcc: dict[str, Any]) -> None:
        assert "orthorhombic (rolling) sample symmetry" in fcc["summary"]
        assert "RP factor" in fcc["summary"]


class TestViewsReuseTheInversion:
    def test_a_different_section_reuses_the_inversion(self, fcc: dict[str, Any]) -> None:
        plate = analyse(FCC, section_preset="labotex", section_resolution_deg=15.0)
        assert plate["data"]["cached_inversion"] is True
        assert (
            plate["data"]["odf"]["inversion_residual"] == fcc["data"]["odf"]["inversion_residual"]
        )
        values = [section["value_deg"] for section in plate["data"]["odf"]["sections"]]
        assert values == [float(value) for value in range(0, 91, 5)]

    def test_phi1_and_custom_sections(self) -> None:
        result = analyse(
            FCC,
            section_kind="phi1",
            section_preset="custom",
            section_values="0, 30",
            section_resolution_deg=15.0,
        )
        odf = result["data"]["odf"]
        assert odf["horizontal_coordinate"] == "phi2"
        assert [section["value_deg"] for section in odf["sections"]] == [0.0, 30.0]

    def test_a_section_outside_the_symmetry_range_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            analyse(FCC, section_preset="custom", section_values="0, 120")
        assert excinfo.value.details["field"] == "section_values"

    def test_chosen_values_with_nothing_typed_are_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            analyse(FCC, section_preset="custom", section_values="")
        assert excinfo.value.details["field"] == "section_values"


class TestSampleSymmetryIsExposed:
    def test_axial_symmetry_makes_every_tilt_ring_constant(self) -> None:
        result = analyse(ZR, sample_symmetry="axial")
        assert result["data"]["sample_symmetry"] == "axial"
        for figure in result["data"]["figures"]:
            rings: dict[float, list[float]] = {}
            for point in figure["points"]:
                rings.setdefault(round(point["polar_deg"], 3), []).append(point["symmetrized"])
            for values in rings.values():
                assert max(values) - min(values) < 1e-9

    def test_hexagonal_sections_span_sixty_degrees_of_phi2(self) -> None:
        result = analyse(ZR, sample_symmetry="axial")
        assert result["data"]["odf"]["ranges"]["phi2_max_deg"] == 60.0
        assert [s["value_deg"] for s in result["data"]["odf"]["sections"]] == [0.0, 30.0]

    def test_an_assumption_the_specimen_breaks_shows_as_a_large_symmetry_change(self) -> None:
        """The split-basal demonstration is orthorhombic, not axial."""

        rolled = analyse(ZR, sample_symmetry="orthorhombic")
        axial = analyse(ZR, sample_symmetry="axial")
        rolled_change = max(row["symmetry_rms"] for row in rolled["table"]["rows"])
        axial_change = max(row["symmetry_rms"] for row in axial["table"]["rows"])
        assert axial_change > 3.0 * rolled_change

    def test_hexagonal_ideal_orientations_are_reported(self) -> None:
        rows = {row["component"]: row for row in analyse(ZR)["data"]["volume_fractions"]["rows"]}
        assert set(rows) == {
            "Basal",
            "Basal, 30° to TD",
            "Basal, 30° to RD",
            "c along TD",
            "c along RD",
        }
        assert rows["Basal, 30° to TD"]["times_random"] > rows["c along RD"]["times_random"]


class TestMeasuredFiles:
    def files(self, count: int) -> dict[str, Any]:
        text = FIXTURE.read_text(encoding="utf-8")
        return {"items": [{"name": f"ni-{index}.xrdml", "text": text} for index in range(count)]}

    def test_opened_files_are_analysed_instead_of_the_demonstration(self) -> None:
        result = analyse(FCC, files=self.files(3), dictionary_count=200)
        data = result["data"]
        assert data["demonstration"] is False
        assert [figure["file"] for figure in data["figures"]] == [
            "ni-0.xrdml",
            "ni-1.xrdml",
            "ni-2.xrdml",
        ]
        assert all(figure["count"] == 12 for figure in data["figures"])
        assert result["inputs"]["files"] == ["ni-0.xrdml", "ni-1.xrdml", "ni-2.xrdml"]

    def test_an_unreadable_file_is_named(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            analyse(FCC, files={"items": [{"name": "broken.xrdml", "text": "<nope/>"}]})
        assert "broken.xrdml" in str(excinfo.value)


def test_the_harmonic_route_reports_a_texture_index() -> None:
    result = analyse(FCC, odf_method="harmonic", odf_bandlimit=4, section_resolution_deg=15.0)
    odf = result["data"]["odf"]
    assert odf["method"] == "harmonic"
    assert odf["texture_index"] >= 1.0
    assert result["data"]["figures"][0]["points"][0]["recalculated"] >= -1.0
