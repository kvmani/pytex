"""Fitting a 2D lattice to picked spots, checked against constructions.

Every case here is built from a lattice whose answer is known before the fit
runs: nodes are generated from stated basis vectors about a stated centre, then
the fit is asked to recover them from the points alone. The interesting cases are
the ones where a plausible wrong answer exists — a sub-lattice that explains the
same points, a centre half a spacing out, a spot clicked somewhere it should not
be — because those are what the overlay is for.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.diffraction.lattice_fit import (
    DEFAULT_INLIER_FRACTION,
    PlanarLatticeFit,
    fit_planar_lattice,
)

CENTRE = np.array([512.0, 384.0])


def nodes_from(basis, indices, *, centre=CENTRE) -> np.ndarray:
    """Lattice points at stated integer indices, about a stated centre."""

    return centre + np.asarray(indices, dtype=float) @ np.asarray(basis, dtype=float)


SQUARE = np.array([[100.0, 0.0], [0.0, 100.0]])
OBLIQUE = np.array([[120.0, 0.0], [40.0, 90.0]])
RING = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, -1), (2, 0), (0, 2)]


# ------------------------------------------------------------- the geometry


def test_an_exact_lattice_is_recovered_exactly() -> None:
    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE)
    assert fit.centre == pytest.approx(CENTRE, abs=1e-9)
    assert fit.rms_residual == pytest.approx(0.0, abs=1e-9)
    assert fit.basis_angle_deg == pytest.approx(90.0, abs=1e-9)
    assert fit.basis_lengths == pytest.approx((100.0, 100.0), abs=1e-9)
    assert fit.inlier_count == len(RING)


def test_an_oblique_lattice_keeps_its_shape() -> None:
    """Cell area is basis-independent; the reduced angle is the readable invariant."""

    fit = fit_planar_lattice(nodes_from(OBLIQUE, RING), CENTRE)
    assert abs(float(np.linalg.det(fit.basis))) == pytest.approx(
        abs(float(np.linalg.det(OBLIQUE))), rel=1e-9
    )
    assert fit.rms_residual == pytest.approx(0.0, abs=1e-9)
    # Gauss reduction puts the included angle of any 2D lattice in [60, 120].
    assert 60.0 - 1e-9 <= fit.basis_angle_deg <= 120.0 + 1e-9


def test_the_reported_basis_is_reduced_so_its_angle_means_something() -> None:
    """A square lattice must not be reported as 135 degrees, which is also true."""

    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE)
    assert fit.basis_angle_deg == pytest.approx(90.0, abs=1e-9)
    assert fit.basis_lengths == pytest.approx((100.0, 100.0), abs=1e-9)


@pytest.mark.parametrize("offset", [(0.0, 0.0), (7.0, -4.0), (30.0, 25.0), (-44.0, 12.0)])
def test_a_badly_picked_centre_is_recovered_from_the_spots(offset) -> None:
    """The point of the exercise: the centre is over-determined by the picks."""

    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE + np.asarray(offset))
    assert fit.centre == pytest.approx(CENTRE, abs=1e-6)
    assert fit.supplied_centre == pytest.approx(CENTRE + np.asarray(offset), abs=1e-12)
    assert fit.centre_shift == pytest.approx(float(np.linalg.norm(offset)), abs=1e-6)


def test_holding_the_centre_leaves_it_exactly_where_it_was_put() -> None:
    supplied = CENTRE + np.array([9.0, -6.0])
    fit = fit_planar_lattice(nodes_from(SQUARE, RING), supplied, refine_centre=False)
    assert fit.centre == pytest.approx(supplied, abs=1e-12)
    assert fit.centre_shift == pytest.approx(0.0, abs=1e-12)


def test_a_centre_wrong_by_an_exact_lattice_vector_cannot_be_detected() -> None:
    """An honest limit, worth stating: geometry alone cannot pick out the beam.

    Move the assumed centre by three lattice vectors and every spot is still an
    exact node of a lattice about the new origin. The fit is perfect, the
    residuals are zero, and nothing in the geometry says anything is wrong —
    because nothing is, geometrically. What identifies the transmitted beam is
    that it is the brightest thing on the plate, which is a judgement about
    intensity, not about position.
    """

    displaced = CENTRE + 3.0 * SQUARE[0]
    fit = fit_planar_lattice(nodes_from(SQUARE, RING), displaced)
    assert fit.centre == pytest.approx(displaced, abs=1e-6)
    assert fit.rms_residual == pytest.approx(0.0, abs=1e-9)


def test_a_centre_beyond_half_a_spacing_is_not_silently_accepted() -> None:
    """The one thing the fit must never do here is look confident.

    Off by more than half a spacing and not by a lattice vector, no candidate can
    both keep the origin near the pick and explain the spots. The fit is then
    held by its leash, or it reports that no lattice explains most of the picks;
    either way it must say so, and both notes point at the transmitted beam,
    which is what the user has to fix.
    """

    fit = fit_planar_lattice(
        nodes_from(SQUARE, RING), CENTRE + np.array([61.0, 58.0]), centre_leash_fraction=0.5
    )
    assert fit.notes
    assert any(
        "half a lattice spacing" in note or "transmitted beam" in note for note in fit.notes
    )


def test_the_leash_holds_the_centre_when_a_candidate_wants_to_run() -> None:
    fit = fit_planar_lattice(
        nodes_from(SQUARE, RING), CENTRE + np.array([61.0, 58.0]), centre_leash_fraction=0.5
    )
    # Whatever it settles on, it stays inside half of its own shortest spacing.
    assert fit.centre_shift <= 0.5 * min(fit.basis_lengths) + 1e-6


def test_the_coarsest_adequate_lattice_wins() -> None:
    """A half-cell explains these points too, and asserts nodes nothing requires."""

    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE)
    area = abs(float(np.linalg.det(fit.basis)))
    assert area == pytest.approx(100.0 * 100.0, rel=1e-9)


def test_a_friedel_pair_cannot_seed_the_basis_alone() -> None:
    """g and -g are collinear through the beam; the seed must look past them."""

    picks = nodes_from(SQUARE, [(1, 0), (-1, 0), (2, 0), (-2, 0), (0, 1), (0, -1)])
    fit = fit_planar_lattice(picks, CENTRE)
    assert fit.rms_residual == pytest.approx(0.0, abs=1e-9)
    assert fit.centre == pytest.approx(CENTRE, abs=1e-6)


def test_spots_on_one_row_are_refused_with_a_usable_message() -> None:
    picks = nodes_from(SQUARE, [(1, 0), (2, 0), (3, 0), (-1, 0)])
    with pytest.raises(ValueError, match="row rather than a lattice"):
        fit_planar_lattice(picks, CENTRE)


# -------------------------------------------------------------- bad picks


def test_a_mis_picked_spot_is_named() -> None:
    """The overlay exists to make this visible; the fit must agree with the eye."""

    picks = nodes_from(SQUARE, RING)
    picks[3] += np.array([38.0, 24.0])
    fit = fit_planar_lattice(picks, CENTRE)
    assert [spot.index for spot in fit.outliers] == [3]
    assert fit.centre == pytest.approx(CENTRE, abs=1.0)
    # One bad pick must not drag the lattice: the rest still fit tightly.
    assert fit.rms_residual < 1.0


def test_one_bad_pick_does_not_blame_its_neighbour() -> None:
    """Plain least squares spreads the error and accuses the wrong spot."""

    for moved in range(len(RING)):
        picks = nodes_from(SQUARE, RING)
        picks[moved] += np.array([0.0, 41.0])
        fit = fit_planar_lattice(picks, CENTRE)
        assert [spot.index for spot in fit.outliers] == [moved]


def test_the_same_spot_picked_twice_is_reported_not_averaged() -> None:
    picks = nodes_from(SQUARE, RING)
    picks = np.vstack([picks, picks[2] + np.array([1.0, 1.0])])
    fit = fit_planar_lattice(picks, CENTRE)
    assert any("same lattice node" in note for note in fit.notes)
    flagged = {spot.index for spot in fit.outliers}
    assert {2, len(picks) - 1} <= flagged


def test_an_outlier_is_labelled_not_removed() -> None:
    picks = nodes_from(SQUARE, RING)
    picks[1] += np.array([45.0, 0.0])
    fit = fit_planar_lattice(picks, CENTRE)
    assert len(fit.spots) == len(RING)
    assert fit.spots[1].inlier is False
    assert fit.spots[1].residual > 20.0


# ------------------------------------------------------------ the overlay


def test_node_positions_cover_the_grid_and_respect_the_frame() -> None:
    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE)
    unbounded = fit.node_positions(max_index=2)
    assert unbounded.shape == (25, 4)
    bounded = fit.node_positions(max_index=6, bounds=(1024.0, 768.0))
    assert bounded.shape[1] == 4
    assert np.all(bounded[:, 0] >= 0.0) and np.all(bounded[:, 0] <= 1024.0)
    assert np.all(bounded[:, 1] >= 0.0) and np.all(bounded[:, 1] <= 768.0)
    origin = bounded[(bounded[:, 2] == 0) & (bounded[:, 3] == 0)]
    assert origin[0, :2] == pytest.approx(CENTRE, abs=1e-9)


def test_node_positions_rejects_an_empty_grid() -> None:
    fit = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE)
    with pytest.raises(ValueError, match="strictly positive"):
        fit.node_positions(max_index=0)


# ---------------------------------------------------------- the reporting


def test_describe_states_what_the_fit_does_not_prove() -> None:
    text = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE).describe()
    assert "geometry, not indexing" in text
    assert "necessary" in text


def test_describe_names_the_spots_that_do_not_fit() -> None:
    picks = nodes_from(SQUARE, RING)
    picks[4] += np.array([50.0, 10.0])
    text = fit_planar_lattice(picks, CENTRE).describe()
    assert "do not lie on the lattice" in text
    assert "spot(s) 5" in text


def test_json_carries_the_fit_and_every_spot() -> None:
    payload = fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE).to_json()
    assert set(payload) == {
        "centre",
        "supplied_centre",
        "centre_shift",
        "basis",
        "basis_lengths",
        "basis_angle_deg",
        "rms_residual",
        "inlier_count",
        "iterations",
        "notes",
        "spots",
        "describe",
    }
    assert len(payload["spots"]) == len(RING)
    assert set(payload["spots"][0]) == {
        "index",
        "x",
        "y",
        "m",
        "n",
        "predicted_x",
        "predicted_y",
        "residual",
        "inlier",
    }


# ---------------------------------------------------------- the guardrails


@pytest.mark.parametrize(
    "kwargs",
    [
        {"inlier_fraction": 0.0},
        {"inlier_fraction": 1.5},
        {"centre_leash_fraction": -0.1},
        {"centre_leash_fraction": 3.0},
        {"max_iterations": 0},
    ],
)
def test_impossible_settings_are_refused(kwargs) -> None:
    with pytest.raises(ValueError):
        fit_planar_lattice(nodes_from(SQUARE, RING), CENTRE, **kwargs)


def test_one_spot_cannot_define_a_lattice() -> None:
    with pytest.raises(ValueError, match="At least two spots"):
        fit_planar_lattice(np.array([[600.0, 400.0]]), CENTRE)


def test_the_default_tolerance_is_a_picking_precision_not_a_lattice_spacing() -> None:
    """Four percent of the shortest spot separation is a few pixels on a plate."""

    assert 0.0 < DEFAULT_INLIER_FRACTION <= 0.1


# ------------------------------------------------- against simulated plates


@pytest.mark.parametrize("entry", ["fcc_al_001", "bcc_fe_110", "hcp_zr_2-1-10"])
@pytest.mark.parametrize("offset", [(0.0, 0.0), (12.0, -9.0), (-28.0, 16.0)])
def test_a_simulated_plate_refines_its_own_beam_centre(entry: str, offset) -> None:
    """End to end on the practice plates, whose true centre is known."""

    pytest.importorskip("matplotlib", reason="the diffraction stack pulls in the plotting layer")
    from pytex.app import REGISTRY

    opened = REGISTRY.call("tem.gallery_pattern", {"pattern": entry})
    picks = opened["data"]["suggested_picks"]
    truth = np.asarray(picks["centre"], dtype=float)
    positions = np.asarray([[spot["x"], spot["y"]] for spot in picks["spots"]], dtype=float)
    fit = fit_planar_lattice(positions, truth + np.asarray(offset))
    # Within the sub-pixel scatter the plate itself carries.
    assert float(np.linalg.norm(fit.centre - truth)) < 3.0
    assert isinstance(fit, PlanarLatticeFit)
