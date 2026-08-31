"""The EBSD scan summary, distributions and discrete figures.

The practice datasets are constructions with known answers, which is what makes
these assertions checks rather than snapshots:

* The **bicrystal** is two orientations with a deformation gradient across one of
  them. It has exactly two grains at any sensible threshold, one boundary between
  them, and a 56 x 56 square grid at a 0.5 micron step — so the point count, the
  scanned area and the grain count are all known before the code runs.
* The **sigma-3 twin** dataset is a coherent twin, whose boundary misorientation
  is 60 degrees about <111> by definition of the relationship. Its
  misorientation-angle distribution must therefore put essentially all of its
  boundary length at 60 degrees.
* A **histogram** has arithmetic that is true of every histogram: the bin edges
  span the data, the counts sum to the population, and the cumulative column
  ends at one.
* A **stereographic projection** of unit directions onto one hemisphere lands
  inside the unit disc, whatever the orientations are.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.registry import REGISTRY


def summary(**overrides):
    request = {"dataset": "bicrystal_gradient"}
    request.update(overrides)
    return REGISTRY.call("ebsd.scan_summary", request)


def distribution(**overrides):
    request = {"dataset": "equiaxed_polycrystal"}
    request.update(overrides)
    return REGISTRY.call("ebsd.distribution", request)


def figure(**overrides):
    request = {"dataset": "equiaxed_polycrystal"}
    request.update(overrides)
    return REGISTRY.call("ebsd.discrete_figure", request)


# ------------------------------------------------------------- the summary


def test_the_summary_reports_the_geometry_the_dataset_was_built_with() -> None:
    result = summary(grid_points=56)
    acquisition = result["data"]["acquisition"]
    assert acquisition["point_count"] == 56 * 56
    assert acquisition["grid_kind"] == "square"
    assert acquisition["grid_shape"] == [56, 56]
    assert acquisition["step_sizes_um"] == [0.5, 0.5]
    assert acquisition["is_multiphase"] is False


def test_the_summary_counts_the_two_grains_the_bicrystal_has() -> None:
    result = summary()
    microstructure = result["data"]["microstructure"]
    assert microstructure["grain_count"] == 2
    # And says at what threshold, because the count means nothing without it.
    assert microstructure["grain_threshold_deg"] == 5.0
    assert "5 deg grain threshold" in result["table"]["caption"]


def test_the_summary_sections_every_row_it_reports() -> None:
    result = summary()
    groups = {row["group"] for row in result["table"]["rows"]}
    assert groups == {"Acquisition", "Indexing quality", "Phases", "Microstructure"}
    for row in result["table"]["rows"]:
        assert row["metric"] and row["value"]
        # Every quantity says what it means; a table of bare numbers is a table
        # nobody can act on.
        assert len(row["note"]) > 10


def test_the_indexed_fraction_moves_with_the_threshold_it_is_quoted_at() -> None:
    channels = summary()["data"]["channels"]
    assert "confidence_index" in channels
    lenient = summary(confidence_threshold=0.0)["data"]["indexed_fraction"]
    strict = summary(confidence_threshold=0.99)["data"]["indexed_fraction"]
    assert lenient == pytest.approx(1.0)
    assert strict <= lenient
    # And the threshold travels with the answer.
    assert summary(confidence_threshold=0.5)["data"]["confidence_threshold"] == 0.5


def test_a_quality_channel_reports_its_shape_and_not_only_its_mean() -> None:
    statistics = summary()["data"]["channels"]["confidence_index"]
    for key in ("mean", "median", "std", "minimum", "maximum", "p05", "p95"):
        assert key in statistics
    assert statistics["minimum"] <= statistics["p05"] <= statistics["median"]
    assert statistics["median"] <= statistics["p95"] <= statistics["maximum"]


# -------------------------------------------------------- the distributions


@pytest.mark.parametrize(
    "quantity",
    ["grain_diameter", "grain_area", "misorientation_angle", "kam", "grod", "confidence_index"],
)
def test_every_histogram_is_arithmetically_a_histogram(quantity: str) -> None:
    result = distribution(quantity=quantity, bins=12)
    rows = result["table"]["rows"]
    assert len(rows) == 12
    edges = result["data"]["edges"]
    assert len(edges) == 13
    assert edges == sorted(edges)

    statistics = result["data"]["statistics"]
    assert edges[0] == pytest.approx(statistics["minimum"])
    assert edges[-1] == pytest.approx(statistics["maximum"])
    assert statistics["minimum"] <= statistics["mean"] <= statistics["maximum"]

    assert sum(row["fraction"] for row in rows) == pytest.approx(1.0)
    assert rows[-1]["cumulative"] == pytest.approx(1.0)
    for previous, row in pairwise(rows):
        assert row["cumulative"] >= previous["cumulative"]
        assert row["lower"] == pytest.approx(previous["upper"])


def test_a_coherent_twin_puts_its_boundary_at_sixty_degrees() -> None:
    """The defining property of a sigma-3 twin, read off the distribution."""

    result = REGISTRY.call(
        "ebsd.distribution",
        {"dataset": "sigma3_twin", "quantity": "misorientation_angle", "bins": 30},
    )
    assert result["data"]["weighted_by_length"] is True
    assert result["data"]["statistics"]["mean"] == pytest.approx(60.0, abs=0.5)
    # Every boundary is at the same angle, so the quantity has no spread at all.
    # That is a distribution and must still be drawable: it comes back as one
    # bin holding the whole boundary length, rather than as thirty bins of which
    # two share it because the value landed on their common edge.
    rows = result["table"]["rows"]
    assert len(rows) == 1
    assert rows[0]["lower"] <= 60.0 <= rows[0]["upper"]
    assert rows[0]["fraction"] == pytest.approx(1.0)


def test_the_misorientation_distribution_carries_its_random_reference() -> None:
    result = distribution(quantity="misorientation_angle", bins=18)
    reference = result["data"]["reference"]
    assert reference, "a measured misorientation distribution is read against a random one"
    assert len(reference["fractions"]) == 18
    assert sum(reference["fractions"]) == pytest.approx(1.0)
    assert reference["pair_count"] > 0
    # Random pairs of a cubic material have a mean disorientation near 45
    # degrees — a property of the cubic group, not of this dataset.
    assert 30.0 < reference["mean"] < 55.0
    # The quantities that are not misorientations have no such reference.
    assert distribution(quantity="kam")["data"]["reference"] is None


def test_a_grain_distribution_counts_grains_and_a_point_one_counts_points() -> None:
    grains = distribution(quantity="grain_diameter")
    points = distribution(quantity="kam")
    assert grains["data"]["population"] == "grains"
    assert points["data"]["population"] == "points"
    assert points["data"]["statistics"]["count"] > grains["data"]["statistics"]["count"]


def test_the_grain_threshold_changes_what_is_being_counted() -> None:
    coarse = distribution(quantity="grain_diameter", grain_threshold_deg=15.0)
    fine = distribution(quantity="grain_diameter", grain_threshold_deg=1.0)
    # A finer threshold resolves subgrains, so there are more grains and they are
    # smaller. This is the dependence the notes exist to state.
    assert fine["data"]["statistics"]["count"] >= coarse["data"]["statistics"]["count"]
    assert fine["data"]["grain_threshold_deg"] == 1.0


def test_an_unknown_dataset_is_refused_rather_than_substituted() -> None:
    with pytest.raises(InvalidInputError) as excinfo:
        REGISTRY.call("ebsd.distribution", {"dataset": "no_such_scan"})
    assert excinfo.value.details["field"] == "dataset"


# ----------------------------------------------------- the discrete figures


def test_every_projected_point_lands_inside_the_unit_disc() -> None:
    for kind in ("pole", "inverse"):
        result = figure(kind=kind, max_points=200)
        points = result["data"]["points"]
        assert points
        for point in points:
            assert math.hypot(point["x"], point["y"]) <= 1.0 + 1e-4


def test_a_pole_figure_draws_the_whole_family_from_each_measurement() -> None:
    result = figure(kind="pole", pole="1 1 1", max_points=100)
    # Cubic {111} has four antipodal pairs, so one measurement contributes more
    # than one point: the figure is of a family, as a measured one inevitably is.
    assert result["data"]["drawn_points"] > result["data"]["measurement_points"]


def test_subsampling_is_reported_rather_than_silent() -> None:
    small = figure(kind="inverse", max_points=20000)
    assert small["data"]["subsampled"] is False
    assert small["data"]["measurement_points"] == small["data"]["scan_points"]

    large = figure(kind="inverse", max_points=100)
    assert large["data"]["subsampled"] is True
    assert large["data"]["measurement_points"] == 100
    assert "subsampled from" in large["summary"]


def test_the_same_scan_gives_the_same_scatter_twice() -> None:
    """A figure that reshuffles on redraw cannot be compared with the one beside it."""

    first = figure(kind="inverse", max_points=150)["data"]["points"]
    second = figure(kind="inverse", max_points=150)["data"]["points"]
    assert first == second
