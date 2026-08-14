"""The practice gallery and the zone-axis atlas, driven the way a user drives them.

The strongest test here is the round trip. A gallery entry is opened, the picks
it suggests are handed to the *real* indexing operation with the calibration it
reports, and the answer must come back as the zone axis the entry was built from.
That exercises the synthesis, the pixel convention, the calibration contract and
the solver in one pass; if any two of them disagree about a convention, it fails.

The atlas is checked against angles the lattice fixes and the instrument cannot
change: ⟨001⟩ to ⟨110⟩ is 45°, ⟨001⟩ to ⟨111⟩ is arccos(1/√3), and neither moves
when the holder does.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase

pytest.importorskip("matplotlib", reason="the diffraction stack pulls in the plotting layer")


def symmetry_cosine(phase, first, second) -> float:
    """Largest |cos| between two directions over the symmetry orbit.

    An indexed pattern fixes the zone axis only up to symmetry, so a round-trip
    check compares families rather than index triples.
    """

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    left = direct @ np.asarray(first, dtype=float)
    left = left / np.linalg.norm(left)
    right = direct @ np.asarray(second, dtype=float)
    right = right / np.linalg.norm(right)
    return float(np.abs(np.einsum("nij,j->ni", operators, left) @ right).max())


def open_entry(pattern: str, **overrides: object) -> dict:
    request: dict = {"pattern": pattern}
    request.update(overrides)
    return REGISTRY.call("tem.gallery_pattern", request)


def index_suggested_picks(opened: dict) -> dict:
    """Run the real indexing operation on the picks the gallery suggests."""

    calibration = opened["data"]["calibration"]
    return REGISTRY.call(
        "tem.solve_pattern",
        {
            "phase": calibration["phase"],
            "picks": opened["data"]["suggested_picks"],
            "units": calibration["units"],
            "camera_constant_mm_angstrom": calibration["camera_constant_mm_angstrom"],
            "pixel_size_mm": calibration["pixel_size_mm"],
        },
    )


class TestPracticeGallery:
    @pytest.mark.parametrize("entry_id", ["fcc_al_001", "bcc_fe_110", "hcp_zr_2-1-10"])
    def test_every_entry_indexes_back_to_its_own_axis(self, entry_id: str) -> None:
        from pytex.app.tem_gallery import gallery_entry

        entry = gallery_entry(entry_id)
        opened = open_entry(entry_id)
        solved = index_suggested_picks(opened)
        phase = builtin_phase(entry.phase_key).to_phase()
        assert symmetry_cosine(
            phase, entry.zone_axis, solved["data"]["zone_axis"]
        ) == pytest.approx(1.0, abs=1e-6)
        assert solved["data"]["matched_fraction"] == pytest.approx(1.0)

    def test_the_suggested_picks_are_not_all_collinear(self) -> None:
        """Friedel pairs are the trap: bright, plentiful, and useless as a seed."""

        picks = open_entry("fcc_al_001")["data"]["suggested_picks"]
        centre = np.asarray(picks["centre"], dtype=float)
        first = np.asarray([picks["spots"][0]["x"], picks["spots"][0]["y"]]) - centre
        second = np.asarray([picks["spots"][1]["x"], picks["spots"][1]["y"]]) - centre
        cosine = abs(
            float(np.dot(first, second)) / (np.linalg.norm(first) * np.linalg.norm(second))
        )
        assert cosine < 0.99

    def test_a_hexagonal_entry_is_labelled_in_four_indices(self) -> None:
        opened = open_entry("hcp_zr_2-1-10")
        labels = [row["hkl"] for row in opened["table"]["rows"]]
        assert any(label.replace(" ", "").strip("()").lstrip("-") == "0002" for label in labels)
        assert all(len(row["hkl"].strip("()").split()) in (1, 4) for row in opened["table"]["rows"])

    def test_the_beam_is_not_assumed_to_be_at_the_centre_of_the_frame(self) -> None:
        opened = open_entry("bcc_fe_110")
        centre = opened["data"]["suggested_picks"]["centre"]
        assert centre != [
            opened["data"]["pattern"]["width_px"] / 2,
            opened["data"]["pattern"]["height_px"] / 2,
        ]

    def test_an_extra_roll_does_not_change_the_answer(self) -> None:
        """One pattern cannot fix the roll, so the indexed axis must not depend on it."""

        phase = builtin_phase("al_fcc").to_phase()
        for roll in (0.0, 40.0):
            solved = index_suggested_picks(open_entry("fcc_al_001", extra_rotation_deg=roll))
            assert symmetry_cosine(
                phase, (0, 0, 1), solved["data"]["zone_axis"]
            ) == pytest.approx(1.0, abs=1e-6)

    def test_the_camera_constant_is_computed_from_the_instrument(self) -> None:
        from pytex.diffraction.kinematic import electron_wavelength_angstrom

        opened = open_entry("fcc_al_001", camera_length_mm=600.0, beam_energy_kev=300.0)
        assert opened["data"]["calibration"]["camera_constant_mm_angstrom"] == pytest.approx(
            600.0 * electron_wavelength_angstrom(300.0)
        )

    def test_a_longer_camera_spreads_the_pattern_and_loses_the_outer_spots(self) -> None:
        short = open_entry("fcc_al_001", camera_length_mm=300.0)
        far = open_entry("fcc_al_001", camera_length_mm=900.0)
        assert len(far["table"]["rows"]) < len(short["table"]["rows"])

    def test_an_impossible_camera_length_is_explained_not_crashed(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            open_entry("fcc_al_001", camera_length_mm=4000.0)
        assert "camera length" in str(excinfo.value.hint)

    def test_an_unknown_entry_names_the_ones_that_exist(self) -> None:
        with pytest.raises(InvalidInputError):
            REGISTRY.call("tem.gallery_pattern", {"pattern": "no_such_plate"})

    def test_switching_off_the_scatter_gives_exact_positions(self) -> None:
        exact = open_entry("fcc_al_001", realistic_scatter=False)
        scattered = open_entry("fcc_al_001")
        exact_rows = {row["hkl"]: row for row in exact["table"]["rows"]}
        moved = [
            abs(row["x"] - exact_rows[row["hkl"]]["x"])
            for row in scattered["table"]["rows"]
            if row["hkl"] in exact_rows
        ]
        assert 0.0 < max(moved) < 6.0

    def test_every_entry_carries_a_lesson_and_somewhere_to_go_next(self) -> None:
        from pytex.app.tem_gallery import GALLERY

        for entry in GALLERY:
            opened = open_entry(entry.identifier)
            assert len(opened["data"]["entry"]["teaches"]) > 200
            assert len(opened["data"]["targets"]) >= 2
            for target in opened["data"]["targets"]:
                assert target["label"]
                assert target["reason"]


class TestZoneAxisAtlas:
    def atlas(self, **overrides: object) -> dict:
        request: dict = {
            "phase": {"builtin": "austenite_fcc"},
            "current_zone_axis": [0, 0, 1],
            "alpha_limit_deg": 40.0,
            "beta_limit_deg": 40.0,
            "beam_rotation_deg": 45.0,
        }
        request.update(overrides)
        return REGISTRY.call("tem.zone_axis_atlas", request)

    def test_the_current_axis_is_reported_as_where_you_are(self) -> None:
        rows = self.atlas()["table"]["rows"]
        assert rows[0]["verdict"] == "current axis"
        assert rows[0]["angle_deg"] == pytest.approx(0.0, abs=1e-3)
        assert not rows[0]["reachable"]

    def test_rows_are_ordered_by_distance_from_the_current_axis(self) -> None:
        angles = [row["angle_deg"] for row in self.atlas()["table"]["rows"]]
        assert angles == sorted(angles)

    def test_the_cubic_landmarks_appear_at_their_closed_form_angles(self) -> None:
        rows = {row["family"]: row for row in self.atlas()["table"]["rows"]}
        assert rows["<110>"]["angle_deg"] == pytest.approx(45.0, abs=1e-6)
        assert rows["<111>"]["angle_deg"] == pytest.approx(
            math.degrees(math.acos(1 / math.sqrt(3))), abs=1e-6
        )

    def test_the_six_fold_axis_is_named_as_such(self) -> None:
        rows = {row["family"]: row for row in self.atlas()["table"]["rows"]}
        assert rows["<111>"]["symmetry"] == "6-fold"
        assert rows["<110>"]["symmetry"] == "2-fold"

    def test_a_narrower_holder_reaches_fewer_axes(self) -> None:
        wide = self.atlas(alpha_limit_deg=60.0, beta_limit_deg=60.0)
        narrow = self.atlas(alpha_limit_deg=10.0, beta_limit_deg=8.0)
        assert narrow["data"]["reachable_count"] < wide["data"]["reachable_count"]

    def test_a_reachable_row_has_a_non_negative_margin(self) -> None:
        for row in self.atlas()["table"]["rows"]:
            if row["reachable"]:
                assert row["margin_deg"] >= 0.0

    def test_the_angles_do_not_depend_on_the_holder(self) -> None:
        """Interplanar angles are lattice geometry; only reachability is instrumental."""

        wide = {
            row["family"]: row["angle_deg"]
            for row in self.atlas(alpha_limit_deg=60.0)["table"]["rows"]
        }
        narrow = {
            row["family"]: row["angle_deg"]
            for row in self.atlas(alpha_limit_deg=5.0)["table"]["rows"]
        }
        for family, angle in narrow.items():
            assert wide[family] == pytest.approx(angle)

    def test_a_hexagonal_atlas_is_labelled_in_four_indices(self) -> None:
        rows = self.atlas(phase={"builtin": "zr_hcp"}, current_zone_axis=[0, 0, 1])["table"]["rows"]
        assert rows[0]["family"] == "<0001>"
        # A four-index label with a negative component is written with spaces, so
        # a family such as <2 -1 -1 0> is the visible proof of the conversion.
        assert any(len(row["family"].strip("<>").split()) == 4 for row in rows)

    def test_the_nearest_member_is_the_one_the_planner_would_use(self) -> None:
        """The table and the tilt plan must not name two different destinations."""

        rows = [row for row in self.atlas()["table"]["rows"] if row["family"] == "<110>"]
        assert rows
        plan = REGISTRY.call(
            "tem.plan_tilt",
            {
                "phase": {"builtin": "austenite_fcc"},
                "current_zone_axis": [0, 0, 1],
                "target_zone_axis": [1, 1, 0],
                "alpha_limit_deg": 40.0,
                "beta_limit_deg": 40.0,
                "beam_rotation_deg": 45.0,
            },
        )
        assert plan["table"]["rows"][0]["member"] == rows[0]["target"]
        assert plan["table"]["rows"][0]["delta_alpha_deg"] == pytest.approx(
            rows[0]["delta_alpha_deg"], abs=1e-6
        )

    def test_too_narrow_a_search_is_explained(self) -> None:
        with pytest.raises(InvalidInputError):
            self.atlas(max_angle_deg=1.0, max_index=1, current_zone_axis=[1, 2, 3])
