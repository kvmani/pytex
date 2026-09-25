"""Map-wide ranking, the preparability raster, the footprint fit and the JSON contract.

A synthetic two-grain map carries the known answers: the left grain has a cube
orientation, so a ``<001>`` zone axis lies in the surface there (``eps* = 0``);
the right grain has ``[111]`` along the surface normal, so every ``<001>`` member
rises ``asin(1/sqrt 3) = 35.26`` degrees out of it --- beyond a +/-30 degree
holder for every mounting rotation that matters.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.ebsd.models import CrystalMap
from pytex.fib import (
    FeasibilityClass,
    ImageRegistration,
    LamellaRankingWeights,
    LamellaSpec,
    fit_footprint,
    guaranteed_radius_deg,
    preparability_raster,
    rank_grains,
    work_order_html,
)
from pytex.tem.stage import EllipticalEnvelope, RectangularEnvelope

SCHEMA = json.loads(
    (Path(__file__).resolve().parents[2] / "schemas" / "fib_lamella_plan.schema.json").read_text(
        encoding="utf-8"
    )
)
EPS_111 = math.degrees(math.asin(1.0 / math.sqrt(3.0)))


def _phase() -> Phase:
    frame = frame_catalog.crystal_frame()
    lattice = Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=frame)
    return Phase(
        "nickel",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=frame),
        crystal_frame=frame,
    )


def _111_along_z() -> np.ndarray:
    """A rotation taking the crystal [111] onto the specimen z axis."""

    target = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    axis = np.cross(target, [0.0, 0.0, 1.0])
    angle = math.acos(float(target[2]))
    return Rotation.from_axis_angle(axis / np.linalg.norm(axis), angle).as_matrix()


def _two_grain_map(rows: int = 40, cols: int = 60, step: float = 1.0) -> CrystalMap:
    phase = _phase()
    matrices = np.empty((rows, cols, 3, 3))
    matrices[:, : cols // 2] = np.eye(3)
    matrices[:, cols // 2 :] = _111_along_z()
    orientations = OrientationSet.from_matrices(
        matrices.reshape(-1, 3, 3),
        crystal_frame=phase.crystal_frame,
        specimen_frame=frame_catalog.specimen_frame(),
        phase=phase,
        symmetry=phase.symmetry,
    )
    yy, xx = np.mgrid[0:rows, 0:cols]
    coordinates = np.column_stack([xx.ravel() * step, yy.ravel() * step])
    return CrystalMap(
        coordinates=coordinates,
        orientations=orientations,
        map_frame=frame_catalog.specimen_frame(),
        grid_shape=(rows, cols),
        grid_kind="square",
        step_sizes=(step, step),
    )


@pytest.fixture(scope="module")
def segmentation():
    return _two_grain_map().segment_grains(max_misorientation_deg=5.0)


@pytest.fixture(scope="module")
def report(segmentation):
    return rank_grains(segmentation, (0, 0, 1))


class TestRaster:
    def test_the_raster_carries_the_known_residuals(self, segmentation) -> None:
        raster = preparability_raster(segmentation.crystal_map, (0, 0, 1))
        assert raster.eps_deg.shape == (40, 60)
        np.testing.assert_allclose(raster.eps_deg[:, :30], 0.0, atol=1e-9)
        np.testing.assert_allclose(raster.eps_deg[:, 30:], EPS_111, atol=1e-9)
        assert raster.fraction_below(10.0) == pytest.approx(0.5)
        assert raster.extent_um == (-0.5, 59.5, -0.5, 39.5)

    def test_a_full_size_map_is_computed_without_a_per_point_loop(self) -> None:
        from pytex.app.ebsd_gallery import build_map

        big = build_map("equiaxed_polycrystal", grid=200, step_um=1.0)
        raster = preparability_raster(big, (1, 1, 1))
        assert raster.eps_deg.shape == (200, 200)
        assert np.isfinite(raster.eps_deg).all()
        assert float(np.nanmax(raster.eps_deg)) <= EPS_111 + 1e-9  # the <111> bound for cubic

    def test_a_hexagonal_scan_is_refused_with_a_reason(self) -> None:
        from types import SimpleNamespace

        # The topology is checked before anything else is read from the map.
        with pytest.raises(ValueError, match="square-grid"):
            preparability_raster(SimpleNamespace(grid_kind="hexagonal"), (0, 0, 1))


class TestRanking:
    def test_the_in_plane_grain_ranks_first_and_is_guaranteed(self, report, segmentation) -> None:
        best = report.best()
        left = int(segmentation.label_grid[0, 0])
        assert best.grain_id == left
        assert best.eps_deg == pytest.approx(0.0, abs=1e-9)
        assert best.feasibility is FeasibilityClass.GUARANTEED
        assert best.sweep.solver == "navigation"  # the best plan is re-solved exactly

    def test_the_out_of_plane_grain_is_not_guaranteed(self, report) -> None:
        other = report.plans[1]
        assert other.eps_deg == pytest.approx(EPS_111, abs=1e-9)
        assert other.feasibility is not FeasibilityClass.GUARANTEED
        assert report.counts()["guaranteed"] == 1

    def test_the_footprint_fits_inside_the_grain(self, report, segmentation) -> None:
        placement = report.best().placement
        assert placement is not None and placement.fits
        # At theta_S = 0 the normal is along x: the 15 um length runs along y
        # in a 40 um tall, 30 um wide grain, so the clearance is set by the
        # length: (40 - 15) / 2 = 12.5 um, resolved to the 1 um step.
        assert report.best().geometry.theta_sample_deg == pytest.approx(0.0, abs=1e-9)
        assert placement.margin_um == pytest.approx(12.0, abs=1.0)
        label_grid = segmentation.label_grid
        for x, y in placement.corners_scan_um:
            assert label_grid[round(y), round(x)] == report.best().grain_id

    def test_identity_registration_places_the_site_at_its_scan_position(self, report) -> None:
        placement = report.best().placement
        np.testing.assert_allclose(placement.center_image, placement.center_scan_um)

    def test_the_weights_are_stated_in_the_prose(self, report) -> None:
        text = report.describe()
        assert "Ranking weights: eps 0.5, fit 0.2, reliability 0.1" in text
        assert "What is undetermined at planning time" in text
        assert set(report.best().score_terms) == {"eps", "fit", "reliability", "edge", "neighbour"}

    def test_changing_the_weights_changes_the_score_not_the_geometry(self, segmentation) -> None:
        heavy = rank_grains(
            segmentation, (0, 0, 1), weights=LamellaRankingWeights(eps=0.1, fit=2.0), shortlist=2
        )
        light = rank_grains(segmentation, (0, 0, 1), shortlist=2)
        assert heavy.best().score != light.best().score
        assert heavy.best().geometry.theta_sample_deg == light.best().geometry.theta_sample_deg

    def test_ranking_is_stable_across_repeated_runs(self, segmentation, report) -> None:
        again = rank_grains(segmentation, (0, 0, 1))
        assert [plan.grain_id for plan in again.plans] == [plan.grain_id for plan in report.plans]
        assert [plan.score for plan in again.plans] == [plan.score for plan in report.plans]

    def test_a_long_lamella_that_does_not_fit_is_flagged(self, segmentation) -> None:
        report = rank_grains(
            segmentation, (0, 0, 1), spec=LamellaSpec(length_um=45.0, width_um=2.0), shortlist=1
        )
        placement = report.best().placement
        assert not placement.fits
        assert placement.largest_length_um == pytest.approx(40.0, abs=1.0)
        assert "footprint_does_not_fit" in {flag.code for flag in report.best().risk_flags}

    def test_a_registered_image_moves_the_site_and_the_azimuth(self, segmentation) -> None:
        registration = ImageRegistration.similarity(
            rotation_deg=90.0, scale=2.0, translation=(100, 0)
        )
        plan = rank_grains(segmentation, (0, 0, 1), registration=registration, shortlist=1).best()
        cx, cy = plan.placement.center_scan_um
        np.testing.assert_allclose(plan.placement.center_image, [100.0 - 2.0 * cy, 2.0 * cx])
        assert plan.theta_image_deg == pytest.approx(90.0, abs=1e-9)


class TestEnvelopeRadius:
    def test_a_symmetric_rectangle_gives_its_smaller_half_range(self) -> None:
        assert guaranteed_radius_deg(RectangularEnvelope(-30, 30, -25, 25)) == 25.0

    def test_a_circle_gives_slightly_less_than_its_radius(self) -> None:
        # The exact tilt curve alpha = asin(cos phi sin e), beta = atan(sin phi tan e)
        # bulges outside the circle alpha^2 + beta^2 = e^2 between the axes, so a
        # 20 degree circular envelope guarantees a little less than 20 degrees --
        # which the small-angle picture would have missed.
        radius = guaranteed_radius_deg(EllipticalEnvelope(20.0, 20.0))
        assert 19.5 < radius < 20.0


class TestFootprint:
    def test_a_hand_computed_rectangle_is_reproduced(self) -> None:
        # A 40 x 20 um grain (step 0.5), lamella 15 x 2 um with its normal along
        # y (theta 90): clearance limited by the short side, (20 - 2) / 2 = 9 um.
        mask = np.zeros((60, 100), dtype=bool)
        mask[10:50, 10:90] = True
        result = fit_footprint(
            mask, step_um=(0.5, 0.5), theta_sample_deg=90.0, length_um=15.0, width_um=2.0
        )
        assert result.fits and result.margin_um == pytest.approx(9.0)
        np.testing.assert_allclose(result.center_scan_um, (24.75, 14.75), atol=1e-9)
        rotated = fit_footprint(
            mask, step_um=(0.5, 0.5), theta_sample_deg=0.0, length_um=15.0, width_um=2.0
        )
        assert rotated.margin_um == pytest.approx(2.5)

    def test_the_largest_length_is_reported_when_it_does_not_fit(self) -> None:
        mask = np.zeros((60, 100), dtype=bool)
        mask[10:50, 10:90] = True
        result = fit_footprint(
            mask, step_um=(0.5, 0.5), theta_sample_deg=0.0, length_um=25.0, width_um=2.0
        )
        assert not result.fits and result.largest_length_um == pytest.approx(20.0)
        assert "does NOT fit" in result.describe()

    def test_an_empty_mask_is_refused(self) -> None:
        with pytest.raises(ValueError):
            fit_footprint(
                np.zeros((5, 5), bool),
                step_um=(1, 1),
                theta_sample_deg=0.0,
                length_um=1.0,
                width_um=0.5,
            )


class TestContract:
    def test_the_report_validates_against_its_schema(self, report) -> None:
        jsonschema.validate(json.loads(json.dumps(report.to_json_dict())), SCHEMA)

    def test_a_single_plan_validates_against_the_plan_definition(self, report) -> None:
        plan_schema = {"$defs": SCHEMA["$defs"], "$ref": "#/$defs/plan"}
        jsonschema.validate(json.loads(json.dumps(report.best().to_json_dict())), plan_schema)

    def test_the_schema_is_listed_in_the_schema_index(self) -> None:
        index = (Path(__file__).resolve().parents[2] / "schemas" / "README.md").read_text(
            encoding="utf-8"
        )
        assert "fib_lamella_plan.schema.json" in index

    def test_the_work_order_is_self_contained_html(self, report) -> None:
        page = work_order_html(report.best())
        assert page.startswith("<!DOCTYPE html>")
        assert "http" not in page.split("<footer>")[0]
        assert "UNCALIBRATED" in page and "+/-" in page
