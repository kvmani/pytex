"""Kearns parameters from three measured scans, against a texture whose f is known.

The demonstration scans are generated from a pilgered-tube texture whose exact
Kearns parameters are the mean of cos^2 over its basal poles. The route sees
only three noisy diffractograms, so recovering those values end to end - peak
detection, matching, random intensities, tilt profile, quadrature - is the claim
under test. It carries the diffractogram route's own interpolation error, so the
tolerance is the one Kearns' method earns, not the exact route's.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.registry import REGISTRY

OPERATION = "kearns.from_three_sections"


def run(**overrides: Any) -> dict[str, Any]:
    spec = next(item for item in REGISTRY.operations() if item.id == OPERATION)
    request = {p.name: p.default for p in spec.parameters if p.default is not None}
    request.update(overrides)
    return REGISTRY.call(OPERATION, request)


@pytest.fixture(scope="module")
def demo() -> dict[str, Any]:
    return run()


def by_label(result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {section["label"]: section for section in result["data"]["sections"]}


class TestTheTriadIsRecovered:
    def test_each_value_is_close_to_the_exact_value_of_the_model(
        self, demo: dict[str, Any]
    ) -> None:
        truth = demo["data"]["truth"]
        for section in demo["data"]["sections"]:
            # Measured within 0.006 of the exact values; 0.02 is the route's budget
            # for interpolating over unevenly spaced reflections.
            assert section["f"] == pytest.approx(truth[section["key"]], abs=0.02), section["key"]

    def test_the_texture_puts_basal_poles_near_radial_and_transverse(
        self, demo: dict[str, Any]
    ) -> None:
        sections = by_label(demo)
        assert sections["r"]["f"] > sections["a"]["f"]
        assert sections["t"]["f"] > sections["a"]["f"]

    def test_the_closure_check_is_genuine(self, demo: dict[str, Any]) -> None:
        data = demo["data"]
        total = sum(section["f"] for section in data["sections"])
        assert data["triad_sum"] == pytest.approx(total)
        assert data["orientation_tensor"] is None
        assert abs(total - 1.0) < 0.1
        assert sum(data["normalized"]) == pytest.approx(1.0)

    def test_the_labels_are_fa_fr_ft_for_a_tube(self, demo: dict[str, Any]) -> None:
        assert [d["label"] for d in demo["data"]["directions"]] == ["a", "r", "t"]
        assert "f_a = " in demo["summary"] and "f_t = " in demo["summary"]

    def test_a_plate_names_the_same_triad_by_its_axes(self) -> None:
        result = run(geometry="plate")
        assert [d["label"] for d in result["data"]["directions"]] == ["RD", "ND", "TD"]


class TestTheEvidenceIsShown:
    def test_each_section_carries_its_diffractogram_peaks_and_reflections(
        self, demo: dict[str, Any]
    ) -> None:
        for section in demo["data"]["sections"]:
            pattern = section["pattern"]
            assert len(pattern["two_theta_deg"]) == len(pattern["intensity"]) > 1000
            assert len(pattern["background"]) == len(pattern["intensity"])
            assert section["peaks"]
            used = [row for row in section["reflections"] if row["used"]]
            assert len(used) >= 6
            planes = {row["plane"] for row in section["reflections"]}
            assert "(0 0 0 2)" in planes or "(0002)" in planes

    def test_matched_peaks_sit_where_the_structure_predicts(self, demo: dict[str, Any]) -> None:
        for section in demo["data"]["sections"]:
            for row in section["reflections"]:
                if row["status"] == "used":
                    assert abs(row["offset_deg"]) < 0.05

    def test_the_quadrature_contributions_sum_to_f(self, demo: dict[str, Any]) -> None:
        for section in demo["data"]["sections"]:
            contributions = sum(node["contribution"] for node in section["quadrature"])
            assert contributions == pytest.approx(section["f"], abs=1e-9)

    def test_every_step_is_a_stage(self, demo: dict[str, Any]) -> None:
        keys = [stage["key"] for stage in demo["stages"]]
        assert keys == [
            "axial_reflections",
            "axial_quadrature",
            "radial_reflections",
            "radial_quadrature",
            "transverse_reflections",
            "transverse_quadrature",
            "triad",
        ]
        assert all(stage.get("table") for stage in demo["stages"][:6])

    def test_the_density_is_measured_over_random(self, demo: dict[str, Any]) -> None:
        for section in demo["data"]["sections"]:
            for row in section["reflections"]:
                if row["status"] == "used":
                    assert row["density"] == pytest.approx(row["measured"] / row["random"])


class TestInputs:
    def xy(self, result: dict[str, Any], key: str) -> dict[str, str]:
        section = next(s for s in result["data"]["sections"] if s["key"] == key)
        lines = [
            f"{t:.3f} {i:.1f}"
            for t, i in zip(
                section["pattern"]["two_theta_deg"], section["pattern"]["intensity"], strict=True
            )
        ]
        return {"name": f"{key}.xy", "text": "\n".join(lines)}

    def test_uploaded_scans_reproduce_the_demonstration(self, demo: dict[str, Any]) -> None:
        files = {key: self.xy(demo, key) for key in ("axial", "radial", "transverse")}
        result = run(scan_files=files)
        assert result["data"]["demonstration"] is False
        assert result["data"]["truth"] is None
        for uploaded, generated in zip(
            result["data"]["sections"], demo["data"]["sections"], strict=True
        ):
            assert uploaded["file"] == f"{uploaded['key']}.xy"
            assert uploaded["f"] == pytest.approx(generated["f"], abs=0.01)

    def test_a_missing_section_is_named(self, demo: dict[str, Any]) -> None:
        files = {key: self.xy(demo, key) for key in ("axial", "radial")}
        with pytest.raises(InvalidInputError) as excinfo:
            run(scan_files=files)
        assert "transverse" in str(excinfo.value)

    def test_a_measured_random_standard_is_required_when_chosen(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            run(random_source="measured")
        assert excinfo.value.details["field"] == "random_file"

    def test_a_cubic_phase_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            run(phase={"builtin": "ni_fcc"})

    def test_peak_height_is_available_and_changes_little(self, demo: dict[str, Any]) -> None:
        heights = run(intensity_measure="height")
        for a, b in zip(heights["data"]["sections"], demo["data"]["sections"], strict=True):
            assert a["f"] == pytest.approx(b["f"], abs=0.05)
            assert np.isfinite(a["f"])
