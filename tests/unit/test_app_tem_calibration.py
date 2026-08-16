"""The pixel-scale calibration: an image that carries a scale bar and nothing else.

The camera equation uses one number, the ratio ``pixel size / camera constant``.
A user with a printed plate or a shared micrograph often knows that ratio — from
a scale bar, or from a standard — and neither of the two numbers it is made of.
The *pixels with a measured scale* mode takes the ratio directly, and the claim
these tests make is that taking it directly changes nothing: the same pattern
indexes to the same answer as it does through a camera constant chosen to have
that quotient. Not bit-for-bit — ``(r * pixel) / camera`` and
``r * (pixel / camera)`` are the same number evaluated in a different order, and
they agree to about 1e-8 relative rather than exactly — but far inside any
tolerance the indexing uses, and far inside the precision the answer is read at.
"""

from __future__ import annotations

import math

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.services.tem import _picking_scale, measured_pattern_from_picks


def practice_plate() -> dict:
    return REGISTRY.call("tem.gallery_pattern", {"pattern": "fcc_al_001"})["data"]["pattern"]


def picks_from(pattern: dict, count: int = 6) -> dict:
    centre = (pattern["width_px"] / 2.0, pattern["height_px"] / 2.0)
    ordered = sorted(
        pattern["spots"],
        key=lambda spot: (spot["x"] - centre[0]) ** 2 + (spot["y"] - centre[1]) ** 2,
    )
    return {
        "centre": [centre[0], centre[1]],
        "spots": [{"x": spot["x"], "y": spot["y"]} for spot in ordered[1 : count + 1]],
    }


def test_a_measured_scale_reproduces_the_camera_constant_answer_exactly() -> None:
    pattern = practice_plate()
    picks = picks_from(pattern)
    camera_constant = pattern["camera_constant_mm_angstrom"]
    pixel_size = pattern["pixel_size_mm"]
    base = {"phase": {"builtin": "al_fcc"}, "picks": picks}

    through_camera = REGISTRY.call(
        "tem.solve_pattern",
        {
            **base,
            "units": "px",
            "camera_constant_mm_angstrom": camera_constant,
            "pixel_size_mm": pixel_size,
        },
    )
    through_scale = REGISTRY.call(
        "tem.solve_pattern",
        {
            **base,
            "units": "px_scale",
            "reciprocal_per_px_angstrom": pixel_size / camera_constant,
        },
    )

    assert through_scale["data"]["zone_axis"] == through_camera["data"]["zone_axis"] == [0, 0, 1]
    assert through_scale["data"]["score"]["score"] == pytest.approx(
        through_camera["data"]["score"]["score"], rel=1e-7
    )
    # Compared per spot on the *measured* spacing, which is what the calibration
    # decides. Which symmetry-equivalent triple a spot is labelled with is not:
    # the solution is determined up to a symmetry operation, so (020) and
    # (2̅0 0) can name the same spot in two equally correct runs.
    scaled_rows = {row["spot"]: row for row in through_scale["table"]["rows"]}
    camera_rows = {row["spot"]: row for row in through_camera["table"]["rows"]}
    assert scaled_rows.keys() == camera_rows.keys()
    for spot, row in scaled_rows.items():
        assert row["d_observed"] == pytest.approx(camera_rows[spot]["d_observed"], rel=1e-7)
        assert row["d_calculated"] == pytest.approx(camera_rows[spot]["d_calculated"], rel=1e-7)


def test_the_scale_means_what_it_says() -> None:
    """One pixel spans exactly the stated amount of reciprocal space.

    Checked against the geometry rather than against the other calibration
    path: a spot ``r`` pixels from the beam is at ``|g| = r s``, so the spacing
    it reports must be ``1 / (r s)``. The pattern object is built directly,
    because an arbitrary scale indexes to nothing — which is the point: the
    conversion has to be right before the indexing can be.
    """

    pattern = practice_plate()
    picks = picks_from(pattern, count=4)
    scale = 0.004
    measured, centre, _ = measured_pattern_from_picks(
        {
            "picks": picks,
            "units": "px_scale",
            "reciprocal_per_px_angstrom": scale,
        }
    )
    assert centre == (picks["centre"][0], picks["centre"][1])
    for spot, original in zip(measured.spots, picks["spots"], strict=True):
        radius = math.hypot(original["x"] - centre[0], original["y"] - centre[1])
        assert math.hypot(*spot.position) == pytest.approx(radius * scale, rel=1e-12)


def test_the_calculated_overlay_is_drawn_at_the_same_scale() -> None:
    """Picking units per inverse angstrom is the inverse of the scale.

    A predicted pattern drawn at a different scale from the measured one it is
    superimposed on would look like a disagreement the crystallography never
    had.
    """

    assert _picking_scale({"units": "px_scale", "reciprocal_per_px_angstrom": 0.05}) == 20.0
    # An unusable calibration returns zero rather than dividing by it, which is
    # how the caller knows not to draw.
    assert _picking_scale({"units": "px_scale", "reciprocal_per_px_angstrom": 0.0}) == 0.0
    assert _picking_scale({"units": "px_scale"}) == 0.0


def test_a_missing_scale_is_refused_and_names_its_field() -> None:
    pattern = practice_plate()
    with pytest.raises(InvalidInputError) as raised:
        REGISTRY.call(
            "tem.solve_pattern",
            {
                "phase": {"builtin": "al_fcc"},
                "picks": picks_from(pattern),
                "units": "px_scale",
                "reciprocal_per_px_angstrom": 0.0,
            },
        )
    assert raised.value.details["field"] == "reciprocal_per_px_angstrom"
    assert "Calibrate" in (raised.value.hint or "")


def test_the_mode_is_offered_in_the_manifest() -> None:
    """A calibration a user cannot select is not one."""

    operation = next(
        entry for entry in REGISTRY.manifest()["operations"] if entry["id"] == "tem.solve_pattern"
    )
    units = next(entry for entry in operation["parameters"] if entry["name"] == "units")
    keys = [
        option[0] if isinstance(option, list) else option["value"] for option in units["options"]
    ]
    assert "px_scale" in keys
    scale = next(
        entry for entry in operation["parameters"] if entry["name"] == "reciprocal_per_px_angstrom"
    )
    assert scale["units"] == "Å⁻¹/px"
    assert "Calibrate" in scale["help"]
