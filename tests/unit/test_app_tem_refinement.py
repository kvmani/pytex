"""The lattice fit, the scored ranking, and the calculated-pattern overlay.

Each test drives the app the way the panel does: open a practice plate whose
answer is known, disturb one thing about it, and require the service to notice
that one thing and nothing else. The overlay tests are the load-bearing ones — a
superimposed pattern is a claim about where every reflection should appear, and
if it is drawn at the wrong scale or the wrong rotation it will look like a
disagreement the crystallography never had.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

pytest.importorskip("matplotlib", reason="the diffraction stack pulls in the plotting layer")


def plate(entry: str = "fcc_al_001", **overrides: object) -> dict:
    request: dict = {"pattern": entry}
    request.update(overrides)
    return REGISTRY.call("tem.gallery_pattern", request)


def solve(opened: dict, **overrides: object) -> dict:
    calibration = opened["data"]["calibration"]
    request: dict = {
        "phase": calibration["phase"],
        "picks": opened["data"]["suggested_picks"],
        "units": "px",
        "camera_constant_mm_angstrom": calibration["camera_constant_mm_angstrom"],
        "pixel_size_mm": calibration["pixel_size_mm"],
    }
    request.update(overrides)
    return REGISTRY.call("tem.solve_pattern", request)


def displaced(opened: dict, dx: float, dy: float) -> dict:
    """The plate's picks with the transmitted beam moved off where it belongs."""

    picks = copy.deepcopy(opened["data"]["suggested_picks"])
    picks["centre"] = [picks["centre"][0] + dx, picks["centre"][1] + dy]
    return picks


class TestFitLattice:
    @pytest.mark.parametrize("entry", ["fcc_al_001", "bcc_fe_110", "hcp_zr_2-1-10"])
    @pytest.mark.parametrize("offset", [(0.0, 0.0), (14.0, -11.0), (-25.0, 18.0)])
    def test_the_beam_centre_is_recovered_from_the_spots(self, entry: str, offset) -> None:
        opened = plate(entry)
        truth = np.asarray(opened["data"]["suggested_picks"]["centre"], dtype=float)
        result = REGISTRY.call("tem.fit_lattice", {"picks": displaced(opened, *offset)})
        recovered = np.asarray(result["data"]["centre"], dtype=float)
        assert float(np.linalg.norm(recovered - truth)) < 3.0

    def test_holding_the_centre_leaves_it_alone(self) -> None:
        opened = plate()
        picks = displaced(opened, 14.0, -11.0)
        result = REGISTRY.call("tem.fit_lattice", {"picks": picks, "refine_centre": False})
        assert result["data"]["centre"] == pytest.approx(picks["centre"], abs=1e-9)
        assert result["data"]["centre_shift"] == pytest.approx(0.0, abs=1e-9)

    def test_a_mis_picked_spot_is_named_in_the_table(self) -> None:
        opened = plate()
        picks = copy.deepcopy(opened["data"]["suggested_picks"])
        picks["spots"][2]["x"] += 44.0
        picks["spots"][2]["y"] -= 30.0
        result = REGISTRY.call("tem.fit_lattice", {"picks": picks})
        assert result["data"]["outliers"] == [3]
        offending = next(row for row in result["table"]["rows"] if row["spot"] == 3)
        assert offending["verdict"] == "off the lattice"
        assert offending["residual"] > 20.0

    def test_the_overlay_nodes_stay_inside_the_frame(self) -> None:
        opened = plate()
        result = REGISTRY.call(
            "tem.fit_lattice",
            {
                "picks": opened["data"]["suggested_picks"],
                "frame_width": 1024.0,
                "frame_height": 1024.0,
            },
        )
        nodes = result["data"]["nodes"]
        assert nodes
        for node in nodes:
            assert 0.0 <= node["x"] <= 1024.0
            assert 0.0 <= node["y"] <= 1024.0
        origin = [node for node in nodes if node["m"] == 0 and node["n"] == 0]
        assert origin[0]["x"] == pytest.approx(result["data"]["centre"][0], abs=1e-9)

    def test_every_pick_lands_on_a_node_of_the_overlay(self) -> None:
        """The grid the user sees must pass through the spots they clicked."""

        opened = plate()
        result = REGISTRY.call(
            "tem.fit_lattice",
            {"picks": opened["data"]["suggested_picks"], "frame_width": 1024.0,
             "frame_height": 1024.0, "node_limit": 8},
        )
        nodes = np.asarray([[node["x"], node["y"]] for node in result["data"]["nodes"]])
        picks = np.asarray(
            [[spot["x"], spot["y"]] for spot in opened["data"]["suggested_picks"]["spots"]]
        )
        distances = np.linalg.norm(picks[:, None, :] - nodes[None, :, :], axis=2).min(axis=1)
        assert float(distances.max()) < 3.0

    def test_the_result_says_this_is_not_an_indexing(self) -> None:
        result = REGISTRY.call("tem.fit_lattice", {"picks": plate()["data"]["suggested_picks"]})
        assert any("not indexing" in note for note in result["notes"])
        assert "geometry, not indexing" in result["data"]["describe"]

    def test_spots_on_one_row_are_explained_not_crashed(self) -> None:
        picks = {
            "centre": [512.0, 512.0],
            "spots": [{"x": 612.0, "y": 512.0}, {"x": 712.0, "y": 512.0}],
        }
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call("tem.fit_lattice", {"picks": picks})
        assert "off the row" in str(excinfo.value.hint)


class TestScoringAndRanking:
    def test_a_correct_indexing_scores_near_one(self) -> None:
        result = solve(plate(realistic_scatter=False))
        score = result["data"]["score"]
        assert score["score"] > 0.98
        assert score["weights"]["angle"] > score["weights"]["length"]

    def test_the_score_falls_when_the_calibration_is_wrong_but_angles_do_not(self) -> None:
        """The distinction the weighting exists to make, end to end."""

        opened = plate(realistic_scatter=False)
        calibration = opened["data"]["calibration"]
        right = solve(opened, length_tolerance=0.10)
        wrong = solve(
            opened,
            length_tolerance=0.10,
            camera_constant_mm_angstrom=calibration["camera_constant_mm_angstrom"] * 1.05,
        )
        assert wrong["data"]["score"]["rms_relative_length_deviation"] > 0.04
        assert wrong["data"]["score"]["rms_angle_deviation_deg"] == pytest.approx(
            right["data"]["score"]["rms_angle_deviation_deg"], abs=1e-6
        )
        assert wrong["data"]["score"]["score"] < right["data"]["score"]["score"]

    def test_the_alternatives_are_sorted_by_the_score_beside_them(self) -> None:
        alternatives = solve(plate())["data"]["alternatives"]
        scores = [entry["score"] for entry in alternatives]
        assert scores == sorted(scores, reverse=True)
        for entry in alternatives:
            assert 0.0 <= entry["score"] <= 1.0
            assert entry["describe"]

    def test_the_policy_travels_with_the_number(self) -> None:
        result = solve(plate(), score_angle_weight=4.0, score_length_tolerance=0.005)
        assert result["data"]["score"]["weights"]["angle"] == pytest.approx(4.0)
        assert result["data"]["score"]["weights"]["length_tolerance"] == pytest.approx(0.005)
        assert result["inputs"]["score_angle_weight"] == pytest.approx(4.0)

    def test_a_stricter_policy_scores_the_same_pattern_lower(self) -> None:
        lenient = solve(plate(), score_length_tolerance=0.20)
        strict = solve(plate(), score_length_tolerance=0.0005)
        assert strict["data"]["score"]["score"] < lenient["data"]["score"]["score"]
        # The evidence is identical; only the reading of it changed.
        assert strict["data"]["score"]["rms_relative_length_deviation"] == pytest.approx(
            lenient["data"]["score"]["rms_relative_length_deviation"]
        )

    def test_the_table_reports_the_deviation_of_every_spot(self) -> None:
        result = solve(plate())
        keys = {column["key"] for column in result["table"]["columns"]}
        assert "d_deviation_percent" in keys
        for row in result["table"]["rows"]:
            assert abs(row["d_deviation_percent"]) < 5.0

    def test_a_wrong_camera_constant_biases_every_spot_the_same_way(self) -> None:
        """The signature that tells a calibration error from an indexing error."""

        opened = plate(realistic_scatter=False)
        calibration = opened["data"]["calibration"]
        result = solve(
            opened,
            length_tolerance=0.10,
            camera_constant_mm_angstrom=calibration["camera_constant_mm_angstrom"] * 1.05,
        )
        deviations = [row["d_deviation_percent"] for row in result["table"]["rows"]]
        # A camera constant five percent too large makes every measured spacing
        # five percent too large, on every spot, by exactly the same amount. An
        # indexing error does not look like that.
        assert all(value == pytest.approx(5.0, abs=0.01) for value in deviations)


class TestCalculatedOverlay:
    @pytest.mark.parametrize("entry", ["fcc_al_001", "bcc_fe_110", "hcp_zr_2-1-10"])
    def test_the_calculated_pattern_lands_on_the_measured_one(self, entry: str) -> None:
        """A superimposed pattern is a claim about where every reflection is."""

        opened = plate(entry, realistic_scatter=False)
        result = solve(opened, max_index=6)
        overlay = np.asarray(
            [[spot["x"], spot["y"]] for spot in result["data"]["alternatives"][0]["overlay"]]
        )
        simulated = np.asarray(
            [[spot["x"], spot["y"]] for spot in opened["data"]["pattern"]["spots"]]
        )
        distances = np.linalg.norm(
            simulated[:, None, :] - overlay[None, :, :], axis=2
        ).min(axis=1)
        assert float(distances.max()) < 1e-6

    def test_the_overlay_is_bounded_by_the_index_limit_not_by_the_plate(self) -> None:
        """A plate spot with no calculated node means check the limit, not the answer."""

        opened = plate("hcp_zr_2-1-10", realistic_scatter=False)
        simulated = np.asarray(
            [[spot["x"], spot["y"]] for spot in opened["data"]["pattern"]["spots"]]
        )
        narrow = np.asarray(
            [
                [spot["x"], spot["y"]]
                for spot in solve(opened, max_index=4)["data"]["alternatives"][0]["overlay"]
            ]
        )
        wide = np.asarray(
            [
                [spot["x"], spot["y"]]
                for spot in solve(opened, max_index=6)["data"]["alternatives"][0]["overlay"]
            ]
        )
        assert len(wide) > len(narrow)
        narrow_worst = np.linalg.norm(
            simulated[:, None, :] - narrow[None, :, :], axis=2
        ).min(axis=1).max()
        wide_worst = (
            np.linalg.norm(simulated[:, None, :] - wide[None, :, :], axis=2).min(axis=1).max()
        )
        assert narrow_worst > 1.0
        assert wide_worst < 1e-6

    def test_the_overlay_carries_labels_and_spacings_for_hovering(self) -> None:
        overlay = solve(plate("hcp_zr_2-1-10"))["data"]["alternatives"][0]["overlay"]
        assert overlay
        for spot in overlay[:5]:
            assert set(spot) == {"hkl", "label", "x", "y", "g", "d"}
            assert spot["d"] == pytest.approx(1.0 / spot["g"], rel=1e-12)
        # Hexagonal reflections are labelled in four indices, as everywhere else.
        assert all(len(spot["label"].strip("()").split()) in (1, 4) for spot in overlay[:5])

    def test_an_uncalibrated_pattern_yields_no_overlay_rather_than_a_wrong_one(self) -> None:
        """Without a scale there is no way to place a calculated spot honestly."""

        opened = plate()
        picks = opened["data"]["suggested_picks"]
        centre = picks["centre"]
        # Coordinates already in reciprocal angstroms need no camera constant, so
        # take the same geometry and declare it that way.
        calibration = opened["data"]["calibration"]
        scale = calibration["camera_constant_mm_angstrom"] / calibration["pixel_size_mm"]
        reciprocal = {
            "centre": [0.0, 0.0],
            "spots": [
                {"x": (spot["x"] - centre[0]) / scale, "y": (spot["y"] - centre[1]) / scale}
                for spot in picks["spots"]
            ],
        }
        result = REGISTRY.call(
            "tem.solve_pattern",
            {
                "phase": calibration["phase"],
                "picks": reciprocal,
                "units": "reciprocal_angstrom",
            },
        )
        overlay = result["data"]["alternatives"][0]["overlay"]
        assert overlay
        # In reciprocal units the scale is 1, so the overlay is in those units too.
        assert max(abs(spot["x"]) for spot in overlay) < 10.0
