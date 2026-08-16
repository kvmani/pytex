"""The Kikuchi overlay is tested as geometry, not as a drawing.

Every assertion here is against something the overlay must satisfy by
construction of the diffraction geometry, and each of those is checkable
without running the service: a band is as wide as the ``000 -> g`` spot distance
for its own plane, its centre line is perpendicular to ``g``, at an exact zone
axis every band of the zone runs through the transmitted beam, and the band
joining two zone axes is the plane both of them lie in.

The first of those is the free self-check the feature gives a user — the band
and the spot are two measurements of the same ``|g|`` — so it is pinned rather
than left implicit.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

#: Å⁻¹ per pixel of the aluminium practice plate, and its beam position.
SCALE = 0.0023924
CENTRE = (512.0, 512.0)


def rotation_z(degrees: float) -> np.ndarray:
    angle = math.radians(degrees)
    return np.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def overlay(matrix: np.ndarray | None = None, **overrides: object) -> dict:
    """Run the overlay for aluminium down [001] unless told otherwise."""

    crystal_to_pattern = np.eye(3) if matrix is None else np.asarray(matrix, dtype=float)
    request: dict[str, object] = {
        "phase": {"builtin": "al_fcc"},
        "orientation": {"crystal_to_pattern": [float(v) for v in crystal_to_pattern.reshape(-1)]},
        "units": "px_scale",
        "reciprocal_per_px_angstrom": SCALE,
        "centre_x": CENTRE[0],
        "centre_y": CENTRE[1],
        "frame_width": 1024.0,
        "frame_height": 1024.0,
    }
    request.update(overrides)
    return REGISTRY.call("tem.kikuchi_overlay", request)


def band_named(result: dict, label: str) -> dict:
    matches = [band for band in result["data"]["bands"] if band["label"] == label]
    assert matches, f"no band labelled {label} in {[b['label'] for b in result['data']['bands']]}"
    return matches[0]


def test_band_width_equals_the_spot_radius_for_the_same_plane() -> None:
    """The band for ``(hkl)`` is as wide as the ``000 -> g`` distance.

    ``width = L * 2 theta_B`` and ``r_g = L * 2 theta_B`` are the same length,
    which is why a user can check the overlay against the plate they just
    looked at. The exact Kossel-cone construction makes the ratio
    ``1 / cos(theta_B)`` rather than exactly one; at 200 kV that is 2e-5, and
    the tolerance below is set to see a factor error, not that.
    """

    result = overlay()
    zone_bands = [band for band in result["data"]["bands"] if band["in_zone"]]
    assert zone_bands, "a zone-axis pattern must show the bands of its own zone"
    for band in zone_bands:
        assert band["radius_px"] == pytest.approx(band["g_inv_angstrom"] / SCALE, rel=1e-12)
        assert band["width_px"] == pytest.approx(band["radius_px"], rel=1e-3)
        # And the direction of the departure is the one the geometry predicts.
        wavelength = result["data"]["wavelength_angstrom"]
        theta = math.asin(min(1.0, band["g_inv_angstrom"] * wavelength / 2.0))
        ratio = band["width_px"] / band["radius_px"]
        assert ratio == pytest.approx(1.0 / math.cos(theta), rel=1e-6)


def test_the_drawn_edges_sit_at_half_the_band_width_from_the_centre() -> None:
    """The polylines drawn are the same geometry as the width that is reported.

    ``width_px`` is computed from the Kossel cones analytically; the edges are
    drawn from the sampled cones. Nothing forces the two to agree unless they
    are the same construction, so this checks that they are — and it also
    catches the edges being drawn as a pair of straight lines, because a real
    edge is a hyperbola and widens away from its closest approach.
    """

    result = overlay()
    for band in result["data"]["bands"]:
        a, b, c = band["line_px"]
        for edge in band["edges"]:
            points = np.asarray([point for run in edge for point in run], dtype=float)
            assert points.size, band["label"]
            offsets = a * points[:, 0] + b * points[:, 1] + c
            assert np.all(np.sign(offsets) == np.sign(offsets[0]))
            closest = float(np.min(np.abs(offsets)))
            assert closest == pytest.approx(band["width_px"] / 2.0, rel=1e-9)
            # Away from the foot the band is wider: the edge is a conic, not a
            # line parallel to the centre.
            assert float(np.max(np.abs(offsets))) > closest


def test_the_band_centre_line_is_perpendicular_to_g() -> None:
    """A plane and its normal are one object: the trace is perpendicular to g."""

    result = overlay()
    for band in result["data"]["bands"]:
        a, b, _ = band["line_px"]
        gx, gy = band["g_direction_px"]
        # `line_px` has (a, b) along the line's normal, which must be g's own
        # in-plane direction; the cross product of two parallel vectors is zero.
        assert abs(a * gy - b * gx) == pytest.approx(0.0, abs=1e-9)


def test_every_band_of_the_zone_passes_through_the_transmitted_beam() -> None:
    """At an exact zone axis the plane contains the beam, so its trace hits 000.

    This is also why the overlay must not label at the pole: at 000 every band
    of the zone crosses, which is the least informative point of the figure.
    """

    result = overlay()
    for band in result["data"]["bands"]:
        a, b, c = band["line_px"]
        distance = abs(a * CENTRE[0] + b * CENTRE[1] + c)
        # The beam is [001] here, so a plane belongs to the zone exactly when
        # g . [001] = l = 0. The condition comes from the lattice, not from the
        # service's own flag, which is checked against it.
        in_zone = band["hkl"][2] == 0
        assert band["in_zone"] is in_zone, band["label"]
        if in_zone:
            assert distance == pytest.approx(0.0, abs=1e-6), band["label"]
        # Labels are placed away from the pole, where the bands separate.
        label_x, label_y = band["label_at"]
        assert math.hypot(label_x - CENTRE[0], label_y - CENTRE[1]) > 60.0


def test_the_connecting_band_from_001_to_011_is_100() -> None:
    """Derived from ``g . [001] = 0`` and ``g . [011] = 0``, not from output.

    Both zone axes lie in the plane whose normal is their cross product,
    ``[001] x [011] = [-100]``, so the plane to follow is ``(100)``. Aluminium
    is face-centred, so the band the crystal actually produces is its first
    allowed order ``(200)`` — the same physical band, and the one a user should
    be sent looking for.
    """

    result = overlay(target_zone_axis="0 1 1")
    connecting = result["data"]["connecting"]
    assert connecting is not None
    assert connecting["plane_hkl"] == [1, 0, 0]
    assert connecting["plane_label"] == "(100)"
    assert connecting["hkl"] == [2, 0, 0]
    assert "follow (200) toward" in connecting["text"]
    assert "[011]" in connecting["text"]
    # The named band is one of the drawn bands, flagged rather than duplicated.
    assert any(band["connecting"] for band in result["data"]["bands"])


def test_a_target_that_no_single_band_reaches_is_said_so() -> None:
    """A zone axis parallel to the current one spans no plane."""

    result = overlay(target_zone_axis="0 0 2")
    assert result["data"]["connecting"] is None
    assert "no single band connects these zones" in result["data"]["connecting_note"].lower()


def test_rotating_the_orientation_rotates_the_bands_with_the_spots() -> None:
    """The overlay is rigid with the crystal, which is the navigation claim."""

    angle = 37.0
    plain = overlay()
    turned = overlay(rotation_z(angle))
    rotation = rotation_z(angle)[:2, :2]
    for band in plain["data"]["bands"]:
        moved = band_named(turned, band["label"])
        expected = rotation @ np.asarray(band["g_direction_px"], dtype=float)
        assert_allclose(moved["g_direction_px"], expected, atol=1e-9)
        assert moved["radius_px"] == pytest.approx(band["radius_px"], rel=1e-12)
        assert moved["width_px"] == pytest.approx(band["width_px"], rel=1e-12)


def test_only_bands_crossing_the_visible_field_are_drawn() -> None:
    """A plane perpendicular to the beam has its trace at infinity."""

    result = overlay()
    half_diagonal = math.hypot(1024.0, 1024.0) / 2.0
    for band in result["data"]["bands"]:
        a, b, c = band["line_px"]
        assert abs(a * CENTRE[0] + b * CENTRE[1] + c) <= half_diagonal + 1e-9
    # (002) is perpendicular to a [001] beam: it is a spot, never a band here.
    assert all(band["hkl"] != [0, 0, 2] for band in result["data"]["bands"])


def test_the_result_states_what_is_and_is_not_modelled() -> None:
    """Positions are geometry; contrast is not, and the overlay is a prediction."""

    result = overlay()
    prose = " ".join(result["notes"]) + result["data"]["describe"]
    lowered = prose.lower()
    assert "excess" in lowered and "deficient" in lowered
    assert "holz" in lowered
    assert "thin" in lowered
    assert "prediction" in lowered or "predicted" in lowered


def test_an_orientation_that_is_not_a_rotation_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        overlay(np.diag([1.0, 1.0, 2.0]))


def test_an_uncalibrated_pattern_is_refused_rather_than_drawn_at_a_guess() -> None:
    with pytest.raises(InvalidInputError):
        overlay(reciprocal_per_px_angstrom=0.0)
