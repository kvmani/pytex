"""Synthetic SAED patterns and the zone-axis atlas, checked against known answers.

Every assertion here has a source outside the code under test. Spot radii are
checked against ``r = (camera constant) / d``, which is the calibration identity
the whole indexing chain rests on. Interplanar angles are checked against the
closed forms of the cubic system, where ⟨001⟩ to ⟨011⟩ is exactly 45° and ⟨001⟩
to ⟨111⟩ is exactly ``arccos(1/sqrt(3))``. Family sizes are the orbit sizes the
cubic point group dictates. And the patterns are round-tripped through the real
indexer: a pattern built from a zone axis must index back to that zone axis, or
the construction and the solver disagree about the convention they share.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.core.lattice import ZoneAxis
from pytex.tem.atlas import pattern_rotational_order, zone_axis_atlas
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

CAMERA_CONSTANT = 10.0317  # 200 kV, 400 mm camera length.
PIXEL_SIZE_MM = 0.024


def phase_of(identifier: str):
    return builtin_phase(identifier).to_phase()


def image_of(identifier: str, indices, **overrides):
    phase = phase_of(identifier)
    axis = ZoneAxis(indices=indices, phase=phase)
    settings = {
        "camera_constant_mm_angstrom": CAMERA_CONSTANT,
        "raster": DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=PIXEL_SIZE_MM),
    }
    settings.update(overrides)
    return synthesize_saed_image(phase, axis, **settings)


# --------------------------------------------------------------- the raster


def test_detector_raster_defaults_to_its_geometric_centre() -> None:
    raster = DetectorRaster(width_px=512, height_px=256, pixel_size_mm=0.05)
    assert raster.centre == (256.0, 128.0)


def test_detector_raster_half_extent_reaches_the_far_corner() -> None:
    raster = DetectorRaster(width_px=100, height_px=100, pixel_size_mm=0.1, centre_px=(0.0, 0.0))
    assert raster.half_extent_mm() == pytest.approx(math.hypot(10.0, 10.0))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"width_px": 0},
        {"height_px": -4},
        {"pixel_size_mm": 0.0},
        {"centre_px": (1.0, float("nan"))},
    ],
)
def test_detector_raster_rejects_impossible_geometry(kwargs) -> None:
    settings = {"width_px": 64, "height_px": 64, "pixel_size_mm": 0.05}
    settings.update(kwargs)
    with pytest.raises(ValueError):
        DetectorRaster(**settings)


# ------------------------------------------------------------- the geometry


def test_every_spot_sits_at_the_camera_constant_over_its_spacing() -> None:
    """The one identity a diffraction pattern is calibrated by."""

    image = image_of("al_fcc", (0, 0, 1))
    assert image.spots
    centre = np.asarray(image.centre_px)
    for spot in image.spots:
        radius_px = float(np.linalg.norm(np.asarray(spot.position_px) - centre))
        radius_mm = radius_px * PIXEL_SIZE_MM
        assert radius_mm == pytest.approx(CAMERA_CONSTANT / spot.d_spacing_angstrom, rel=1e-12)
        assert spot.d_spacing_angstrom == pytest.approx(1.0 / spot.g_inv_angstrom, rel=1e-12)


def test_every_spot_obeys_the_zone_law() -> None:
    """``hu + kv + lw = 0`` is what makes a reflection part of the pattern."""

    for identifier, axis in (("al_fcc", (0, 0, 1)), ("fe_bcc", (1, 1, 0)), ("zr_hcp", (1, 0, 0))):
        image = image_of(identifier, axis)
        assert image.spots
        for spot in image.spots:
            assert int(np.dot(spot.miller_indices, np.asarray(axis))) == 0


def test_fcc_shows_only_unmixed_indices() -> None:
    """Face centring extinguishes every reflection of mixed parity."""

    image = image_of("al_fcc", (0, 0, 1))
    for spot in image.spots:
        parities = {int(value) % 2 for value in spot.miller_indices}
        assert len(parities) == 1


def test_bcc_shows_only_even_index_sums() -> None:
    """Body centring extinguishes ``h + k + l`` odd."""

    image = image_of("fe_bcc", (1, 1, 0))
    for spot in image.spots:
        assert int(np.sum(spot.miller_indices)) % 2 == 0


def test_spots_outside_the_frame_are_not_simulated() -> None:
    image = image_of("al_fcc", (0, 0, 1))
    for spot in image.spots:
        assert 0.0 <= float(spot.position_px[0]) <= 1024.0
        assert 0.0 <= float(spot.position_px[1]) <= 1024.0


def test_an_off_centre_beam_moves_every_spot_with_it() -> None:
    """The workflow must not assume the beam is at the middle of the frame."""

    centred = image_of("al_fcc", (0, 0, 1))
    offset = image_of(
        "al_fcc",
        (0, 0, 1),
        raster=DetectorRaster(
            width_px=1024, height_px=1024, pixel_size_mm=PIXEL_SIZE_MM, centre_px=(400.0, 620.0)
        ),
    )
    assert offset.centre_px == (400.0, 620.0)
    common = {tuple(int(v) for v in spot.miller_indices) for spot in centred.spots} & {
        tuple(int(v) for v in spot.miller_indices) for spot in offset.spots
    }
    assert common
    centred_by_index = {tuple(int(v) for v in s.miller_indices): s for s in centred.spots}
    offset_by_index = {tuple(int(v) for v in s.miller_indices): s for s in offset.spots}
    for indices in common:
        delta = np.asarray(offset_by_index[indices].position_px) - np.asarray(
            centred_by_index[indices].position_px
        )
        assert delta == pytest.approx(np.array([400.0 - 512.0, 620.0 - 512.0]), abs=1e-9)


def test_the_in_plane_roll_rotates_the_pattern_rigidly() -> None:
    """A roll about the beam moves spots along circles about the beam.

    It does not preserve the *set* of visible reflections, and should not: the
    frame is a square, so rolling the pattern carries reflections through the
    corners and out through the edges exactly as a real plate does. What must
    hold is that every reflection appearing in both keeps its radius.
    """

    upright = image_of("al_fcc", (0, 0, 1))
    rolled = image_of("al_fcc", (0, 0, 1), in_plane_rotation_deg=31.0)
    centre = np.asarray(upright.centre_px)
    upright_by_index = {tuple(int(v) for v in s.miller_indices): s for s in upright.spots}
    shared = 0
    for spot in rolled.spots:
        source = upright_by_index.get(tuple(int(v) for v in spot.miller_indices))
        if source is None:
            continue
        shared += 1
        assert np.linalg.norm(np.asarray(spot.position_px) - centre) == pytest.approx(
            np.linalg.norm(np.asarray(source.position_px) - centre), rel=1e-12
        )
    assert shared >= 4


def test_jitter_is_deterministic_and_bounded() -> None:
    """A gallery entry must look the same every time it is opened."""

    first = image_of("al_fcc", (0, 0, 1), position_jitter_px=1.0, rng_seed=7)
    second = image_of("al_fcc", (0, 0, 1), position_jitter_px=1.0, rng_seed=7)
    exact = image_of("al_fcc", (0, 0, 1))
    assert [tuple(s.position_px) for s in first.spots] == [
        tuple(s.position_px) for s in second.spots
    ]
    exact_by_index = {tuple(int(v) for v in s.miller_indices): s for s in exact.spots}
    displacements = [
        float(np.linalg.norm(np.asarray(s.position_px) - exact_by_index[
            tuple(int(v) for v in s.miller_indices)
        ].position_px))
        for s in first.spots
    ]
    assert max(displacements) < 6.0
    assert max(displacements) > 0.0


def test_intensities_are_normalized_and_ordered() -> None:
    image = image_of("al_fcc", (0, 0, 1))
    intensities = [spot.relative_intensity for spot in image.spots]
    assert max(intensities) == pytest.approx(1.0)
    assert intensities == sorted(intensities, reverse=True)


def test_brighter_spots_are_drawn_larger_but_bounded() -> None:
    image = image_of("zr_hcp", (1, 0, 0))
    radii = [spot.apparent_radius_px for spot in image.spots]
    assert min(radii) > 0.0
    assert max(radii) / min(radii) < 4.0


def test_independent_seed_spots_span_the_zone() -> None:
    """Friedel pairs are collinear through the beam and cannot seed indexing."""

    image = image_of("al_fcc", (0, 0, 1))
    seeds = image.independent_seed_spots(4)
    assert len(seeds) >= 2
    centre = np.asarray(image.centre_px)
    first = np.asarray(seeds[0].position_px) - centre
    second = np.asarray(seeds[1].position_px) - centre
    cosine = abs(float(np.dot(first, second)) / (np.linalg.norm(first) * np.linalg.norm(second)))
    assert cosine < 0.999


def test_max_spots_caps_the_returned_pattern() -> None:
    image = image_of("zr_hcp", (1, 0, 0), max_spots=5)
    assert len(image.spots) == 5


def test_describe_states_the_phase_axis_and_kinematic_limit() -> None:
    text = image_of("al_fcc", (0, 0, 1)).describe()
    assert "Aluminium" in text
    assert "[001]" in text
    assert "kinematic" in text


def test_json_round_trips_the_answer_with_the_picture() -> None:
    payload = image_of("fe_bcc", (1, 1, 0)).to_json()
    assert payload["zone_axis"] == [1, 1, 0]
    assert payload["spots"]
    assert set(payload["spots"][0]) == {
        "hkl",
        "label",
        "x",
        "y",
        "g_inv_angstrom",
        "d_angstrom",
        "intensity",
        "radius_px",
        "double_diffraction",
        "double_diffraction_origin",
        "double_diffraction_parents",
    }
    # A purely kinematic pattern says so on every spot rather than by omitting
    # the field: a renderer that has to distinguish "not double diffraction"
    # from "this service is too old to say" cannot be written safely.
    assert all(spot["double_diffraction"] is False for spot in payload["spots"])
    assert all(spot["double_diffraction_origin"] == "" for spot in payload["spots"])
    assert all(spot["double_diffraction_parents"] is None for spot in payload["spots"])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"intensity_floor": 1.0},
        {"max_spots": 0},
        {"base_radius_px": 0.0},
        {"bloom": -1.0},
        {"position_jitter_px": -0.5},
    ],
)
def test_synthesis_rejects_impossible_settings(kwargs) -> None:
    with pytest.raises(ValueError):
        image_of("al_fcc", (0, 0, 1), **kwargs)


def test_synthesis_rejects_a_zone_axis_of_another_phase() -> None:
    other = phase_of("fe_bcc")
    with pytest.raises(ValueError, match="must match phase"):
        synthesize_saed_image(phase_of("al_fcc"), ZoneAxis(indices=(0, 0, 1), phase=other))


# ------------------------------------------------------- the indexing round trip


def solve_image(image, *, length_tolerance=0.03):
    from pytex.diffraction.solving import (
        MeasuredSAEDPattern,
        MeasuredSpot,
        PatternCalibration,
        solve_saed_pattern,
    )

    calibration = PatternCalibration(
        units="px",
        centre=image.centre_px,
        camera_constant_mm_angstrom=image.camera_constant_mm_angstrom,
        pixel_size_mm=image.raster.pixel_size_mm,
    )
    pattern = MeasuredSAEDPattern(
        name="synthetic",
        spots=tuple(
            MeasuredSpot(position=tuple(float(v) for v in spot.position_px))
            for spot in image.spots[:10]
        ),
        calibration=calibration,
    )
    return solve_saed_pattern(
        pattern,
        [image.phase],
        max_index=4,
        length_tolerance_relative=length_tolerance,
        angle_tolerance_deg=2.0,
    )


def symmetry_angle_deg(phase, first, second) -> float:
    """Smallest angle between two lattice directions over the symmetry orbit.

    An indexed pattern fixes the zone axis only up to symmetry — a bcc [110]
    pattern is indistinguishable from a [101] one — so a round-trip check must
    compare families, not index triples.
    """

    direct = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    operators = np.asarray(phase.symmetry.operators, dtype=float)
    left = direct @ np.asarray(first, dtype=float)
    left = left / np.linalg.norm(left)
    right = direct @ np.asarray(second, dtype=float)
    right = right / np.linalg.norm(right)
    cosines = np.abs(np.einsum("nij,j->ni", operators, left) @ right)
    return float(math.degrees(math.acos(float(np.clip(cosines.max(), -1.0, 1.0)))))


@pytest.mark.parametrize(
    ("identifier", "axis"),
    [("al_fcc", (0, 0, 1)), ("fe_bcc", (1, 1, 0)), ("zr_hcp", (1, 0, 0))],
)
def test_a_synthesized_pattern_indexes_back_to_its_own_zone_axis(identifier, axis) -> None:
    """The construction and the solver must agree about their shared convention."""

    image = image_of(identifier, axis, in_plane_rotation_deg=17.0)
    report = solve_image(image)
    assert report.solutions
    best = report.best()
    assert best.phase_name == image.phase.name
    assert symmetry_angle_deg(image.phase, best.zone_axis.indices, axis) == pytest.approx(
        0.0, abs=1e-6
    )
    assert best.matched_fraction == pytest.approx(1.0)


def test_a_jittered_pattern_still_indexes_with_a_realistic_residual() -> None:
    image = image_of("al_fcc", (0, 0, 1), position_jitter_px=1.0, rng_seed=3)
    report = solve_image(image)
    best = report.best()
    assert best.matched_fraction == pytest.approx(1.0)
    assert 0.0 < best.mean_residual_inv_angstrom < 0.01


def test_the_roll_about_the_beam_does_not_change_the_indexed_axis() -> None:
    """One indexed pattern cannot fix the roll — so the answer must not depend on it."""

    phase = phase_of("al_fcc")
    for roll in (0.0, 23.0, 61.0):
        report = solve_image(image_of("al_fcc", (0, 0, 1), in_plane_rotation_deg=roll))
        assert symmetry_angle_deg(
            phase, report.best().zone_axis.indices, (0, 0, 1)
        ) == pytest.approx(0.0, abs=1e-6)


# ------------------------------------------------------------------- the atlas


def test_pattern_rotational_order_reads_a_square_and_a_hexagon() -> None:
    square = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    assert pattern_rotational_order(square, np.ones(4)) == 4
    hexagon = np.array(
        [[math.cos(k * math.pi / 3), math.sin(k * math.pi / 3)] for k in range(6)]
    )
    assert pattern_rotational_order(hexagon, np.ones(6)) == 6
    assert pattern_rotational_order(np.array([[1.0, 0.0], [-1.0, 0.0]]), np.ones(2)) == 2
    assert pattern_rotational_order(np.array([[1.0, 0.0]]), np.ones(1)) == 1


def test_pattern_rotational_order_respects_intensity() -> None:
    """Four spots in a square are not four-fold if one of them is brighter."""

    square = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    assert pattern_rotational_order(square, np.array([1.0, 0.5, 1.0, 0.5])) == 2


def test_cubic_pattern_symmetry_matches_the_zone_axis() -> None:
    phase = phase_of("al_fcc")
    atlas = zone_axis_atlas(phase, max_index=1)
    by_label = {entry.label: entry for entry in atlas.entries}
    assert by_label["[100]"].rotational_order == 4
    assert by_label["[110]"].rotational_order == 2
    assert by_label["[111]"].rotational_order == 6


def test_cubic_family_sizes_are_the_point_group_orbits() -> None:
    atlas = zone_axis_atlas(phase_of("al_fcc"), max_index=1)
    sizes = {entry.label: entry.family_size for entry in atlas.entries}
    assert sizes["[100]"] == 3
    assert sizes["[110]"] == 6
    assert sizes["[111]"] == 4
    families = {entry.label: entry.family_label for entry in atlas.entries}
    assert families["[110]"] == "<110>"


def test_cubic_interplanar_angles_match_the_closed_form() -> None:
    phase = phase_of("al_fcc")
    axis = ZoneAxis(indices=(0, 0, 1), phase=phase)
    atlas = zone_axis_atlas(phase, current_zone_axis=axis, max_index=1)
    angles = {entry.label: entry.angle_from_current_deg for entry in atlas.entries}
    assert angles["[001]"] == pytest.approx(0.0, abs=1e-9)
    assert angles["[110]"] == pytest.approx(45.0, abs=1e-9)
    assert angles["[111]"] == pytest.approx(math.degrees(math.acos(1 / math.sqrt(3))), abs=1e-9)


def test_the_current_family_is_labelled_with_the_indices_the_user_has() -> None:
    """[001] and [100] are one cubic family; the row must say where you are."""

    phase = phase_of("al_fcc")
    atlas = zone_axis_atlas(
        phase, current_zone_axis=ZoneAxis(indices=(0, 0, 1), phase=phase), max_index=1
    )
    assert atlas.entries[0].label == "[001]"
    assert atlas.entries[0].angle_from_current_deg == pytest.approx(0.0, abs=1e-9)


def test_the_nearest_member_is_the_one_worth_tilting_to() -> None:
    phase = phase_of("al_fcc")
    axis = ZoneAxis(indices=(0, 0, 1), phase=phase)
    atlas = zone_axis_atlas(phase, current_zone_axis=axis, max_index=1)
    entry = next(item for item in atlas.entries if item.label == "[110]")
    member = np.asarray(entry.nearest_member, dtype=float)
    cosine = abs(float(np.dot(member, [0.0, 0.0, 1.0])) / np.linalg.norm(member))
    assert math.degrees(math.acos(cosine)) == pytest.approx(45.0, abs=1e-9)
    assert int(np.dot(entry.nearest_member, [0, 0, 1])) != 0


def test_hexagonal_basal_and_prism_axes_are_ninety_degrees_apart() -> None:
    """True for every c/a, which is why the pair is the standard hcp check."""

    phase = phase_of("zr_hcp")
    atlas = zone_axis_atlas(
        phase, current_zone_axis=ZoneAxis(indices=(0, 0, 1), phase=phase), max_index=1
    )
    prism = next(entry for entry in atlas.entries if entry.label == "[100]")
    assert prism.angle_from_current_deg == pytest.approx(90.0, abs=1e-9)


def test_the_angle_filter_drops_distant_families() -> None:
    phase = phase_of("al_fcc")
    axis = ZoneAxis(indices=(0, 0, 1), phase=phase)
    atlas = zone_axis_atlas(phase, current_zone_axis=axis, max_index=2, max_angle_deg=30.0)
    assert atlas.entries
    assert all(entry.angle_from_current_deg <= 30.0 for entry in atlas.entries)


def test_entries_are_ranked_by_distance_then_richness() -> None:
    phase = phase_of("fe_bcc")
    axis = ZoneAxis(indices=(1, 1, 0), phase=phase)
    atlas = zone_axis_atlas(phase, current_zone_axis=axis, max_index=3, limit=12)
    angles = [entry.angle_from_current_deg for entry in atlas.entries]
    assert angles == sorted(angles)
    assert len(atlas.entries) == 12


def test_without_a_current_axis_the_richest_pattern_comes_first() -> None:
    atlas = zone_axis_atlas(phase_of("al_fcc"), max_index=2, limit=5)
    counts = [entry.reflection_count for entry in atlas.entries]
    assert counts == sorted(counts, reverse=True)
    assert all(math.isnan(entry.angle_from_current_deg) for entry in atlas.entries)


def test_atlas_describe_names_its_bounds_and_its_ranking() -> None:
    phase = phase_of("al_fcc")
    text = zone_axis_atlas(
        phase, current_zone_axis=ZoneAxis(indices=(0, 0, 1), phase=phase), max_index=1
    ).describe()
    assert "Aluminium" in text
    assert "kinematic" in text
    assert "[001]" in text


def test_atlas_rejects_inconsistent_requests() -> None:
    phase = phase_of("al_fcc")
    other = ZoneAxis(indices=(0, 0, 1), phase=phase_of("fe_bcc"))
    with pytest.raises(ValueError, match="must match phase"):
        zone_axis_atlas(phase, current_zone_axis=other)
    with pytest.raises(ValueError, match="requires a current_zone_axis"):
        zone_axis_atlas(phase, max_angle_deg=20.0)
    with pytest.raises(ValueError, match="strictly positive"):
        zone_axis_atlas(phase, max_index=0)
    with pytest.raises(ValueError, match="strictly positive"):
        zone_axis_atlas(phase, limit=0)
