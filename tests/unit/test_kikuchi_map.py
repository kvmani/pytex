"""Tests for the stereographic Kikuchi map and its routing.

The identities pinned here are integer or closed-form, not tolerances against a
prior run: the Weiss zone law is exact arithmetic, the band width is
``2 arcsin(lambda / 2d)``, and the angles between low-index cubic and hexagonal
zone axes have exact values.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import FrameDomain, ReferenceFrame, get_phase_fixture
from pytex.diffraction.kikuchi_map import (
    KikuchiMapBand,
    KikuchiMapZoneAxis,
    StereographicKikuchiMap,
    compute_kikuchi_map,
    plan_kikuchi_route,
)

# Every test here starts from a phase loaded out of its CIF fixture, and
# CIF-backed phase creation is pymatgen's job. Declared at module scope because
# it is the whole module: without the optional `adapters` extra these are not
# failures but tests that cannot be run, and reporting them as failures hid a
# real defect underneath them for as long as it lasted.
pytest.importorskip(
    "pymatgen",
    reason="loading a phase from its CIF fixture needs the 'adapters' extra",
)

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


@pytest.fixture(scope="module")
def nickel():
    return get_phase_fixture("ni_fcc").load_phase(crystal_frame=CRYSTAL)


@pytest.fixture(scope="module")
def zirconium():
    return get_phase_fixture("zr_hcp").load_phase(crystal_frame=CRYSTAL)


@pytest.fixture(scope="module")
def cubic_map(nickel):
    return compute_kikuchi_map(
        nickel, beam_energy_kev=200.0, max_index=4, zone_axis_max_index=3
    )


@pytest.fixture(scope="module")
def hexagonal_map(zirconium):
    return compute_kikuchi_map(
        zirconium, beam_energy_kev=200.0, max_index=3, zone_axis_max_index=3
    )


# ---------------------------------------------------------------------------
# Geometry: the two exact identities
# ---------------------------------------------------------------------------


def test_every_recorded_band_satisfies_the_weiss_zone_law(cubic_map, hexagonal_map) -> None:
    """A band belongs to a zone axis only if ``hu + kv + lw`` is exactly zero.

    This is integer arithmetic, so the check admits no tolerance. It is the one
    statement that would break if the map frame, the cross products, or the
    rationalization of a zone axis disagreed about which frame they are in.
    """

    for kikuchi_map in (cubic_map, hexagonal_map):
        for axis in kikuchi_map.zone_axes:
            direction = np.asarray(axis.indices, dtype=np.int64)
            for position in axis.band_indices:
                band = kikuchi_map.bands[position]
                plane = np.asarray(band.indices, dtype=np.int64)
                assert int(plane @ direction) == 0, (axis.indices, band.indices)


def test_band_width_is_twice_the_bragg_angle_and_grows_as_the_spacing_falls(cubic_map) -> None:
    r"""``2 theta_B = 2 arcsin(lambda / 2d)``, and it is a *decreasing* function of ``d``.

    The second half is the part that is easy to state backwards. Width goes as
    ``lambda / d`` for small angles, so the widest bands on a map come from the
    smallest spacings — the high-index planes — while the strongest come from the
    largest. Two different orderings.
    """

    wavelength = cubic_map.wavelength_angstrom
    for band in cubic_map.bands:
        expected = 2.0 * math.degrees(math.asin(wavelength / (2.0 * band.d_spacing_angstrom)))
        assert band.angular_width_deg == pytest.approx(expected, rel=1e-12)

    spacings = np.array([band.d_spacing_angstrom for band in cubic_map.bands])
    widths = np.array([band.angular_width_deg for band in cubic_map.bands])
    order = np.argsort(spacings)
    assert np.all(np.diff(widths[order]) <= 1e-12), "width must fall as spacing rises"
    widest = max(cubic_map.bands, key=lambda band: band.angular_width_deg)
    finest = min(cubic_map.bands, key=lambda band: band.d_spacing_angstrom)
    assert widest.indices == finest.indices


def test_wavelength_scales_every_band_width_and_nothing_else(nickel) -> None:
    """The accelerating voltage enters only through the wavelength.

    Halving the wavelength halves every Bragg angle to first order and leaves the
    band network — which planes, which traces, which zone axes — untouched. That
    separation is worth pinning, because it is why one map serves any voltage once
    the widths are rescaled.
    """

    low = compute_kikuchi_map(nickel, beam_energy_kev=100.0, max_index=3, zone_axis_max_index=2)
    high = compute_kikuchi_map(nickel, beam_energy_kev=300.0, max_index=3, zone_axis_max_index=2)
    assert [band.indices for band in low.bands] == [band.indices for band in high.bands]
    assert [axis.indices for axis in low.zone_axes] == [axis.indices for axis in high.zone_axes]
    ratio = np.array(
        [
            high_band.angular_width_deg / low_band.angular_width_deg
            for low_band, high_band in zip(low.bands, high.bands, strict=True)
        ]
    )
    expected = high.wavelength_angstrom / low.wavelength_angstrom
    assert_allclose(ratio, expected, rtol=1e-4)


# ---------------------------------------------------------------------------
# The cubic map against the standard projection
# ---------------------------------------------------------------------------


def test_the_map_centre_is_a_zone_axis_whose_bands_are_the_vertical_ones(cubic_map) -> None:
    """Centred on [001], the axis at the centre of the projection is [001].

    Its bands are exactly the planes with ``l = 0``: the zone law for [001] reads
    ``l = 0``, so the traces through the centre of the map are the ones whose
    normals lie in the plane of the projection.
    """

    centre = cubic_map.zone_axis_for_direction([0, 0, 1])
    assert centre is not None
    assert centre.polar_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert_allclose(centre.projected(), [0.0, 0.0], atol=1e-12)
    for band in cubic_map.bands_through([0, 0, 1]):
        assert band.indices[2] == 0


def test_the_strongest_fcc_bands_are_the_close_packed_planes(cubic_map) -> None:
    """{111} then {200} then {220}, which is the fcc reflection sequence.

    The ordering comes from the Mott-Bethe structure factor rather than from the
    atomic-number proxy. That matters: the proxy is s-independent, so on a
    monatomic fixture with no Debye-Waller factor it gives every allowed
    reflection the same intensity and the band ordering carries no information.
    """

    strongest = cubic_map.bands[0]
    assert sorted(abs(value) for value in strongest.indices) == [1, 1, 1]
    assert strongest.relative_intensity == pytest.approx(1.0)
    by_family: dict[tuple[int, ...], float] = {}
    for band in cubic_map.bands:
        key = tuple(sorted(abs(value) for value in band.indices))
        by_family.setdefault(key, band.relative_intensity)
    assert by_family[(1, 1, 1)] > by_family[(0, 0, 2)] > by_family[(0, 2, 2)]


def test_higher_orders_are_one_band_unless_asked_for(nickel) -> None:
    """All orders of a reflection share one centre line, so they are one band.

    Keeping (222) beside (111) would draw coincident lines and would double the
    band count at every zone axis they both pass through — and that count is the
    number an operator reads to judge whether an intersection is unmistakable.
    """

    folded = compute_kikuchi_map(nickel, max_index=4, zone_axis_max_index=2)
    unfolded = compute_kikuchi_map(
        nickel, max_index=4, zone_axis_max_index=2, include_higher_orders=True
    )
    assert folded.band_for_plane([1, 1, 1]) is not None
    assert folded.band_for_plane([2, 2, 2]) is None
    assert unfolded.band_for_plane([2, 2, 2]) is not None
    assert unfolded.band_count > folded.band_count

    first = unfolded.band_for_plane([1, 1, 1])
    second = unfolded.band_for_plane([2, 2, 2])
    assert first is not None and second is not None
    # Same trace, different width: the centre lines coincide and the Bragg angles
    # do not.
    assert_allclose(np.abs(first.normal_map), np.abs(second.normal_map), atol=1e-12)
    assert second.angular_width_deg > first.angular_width_deg


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------


def test_a_single_band_joins_the_cubic_zone_axes_at_the_exact_angle(cubic_map) -> None:
    r"""[001] to [111] is one band and exactly :math:`\arccos(1/\sqrt{3})`.

    Both zone axes are perpendicular to (1-10), so that band runs between them,
    and the great-circle arc along it is the geodesic. The angle is closed-form:
    :math:`54.735610^\circ`.
    """

    shared = cubic_map.shared_bands([0, 0, 1], [1, 1, 1])
    assert shared, "a (1-10)-type band must join [001] and [111]"
    assert all(sum(band.indices) == 0 for band in shared)

    route = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=60.0)
    assert route.reachable
    assert route.hop_count == 1
    exact = math.degrees(math.acos(1.0 / math.sqrt(3.0)))
    assert route.total_tilt_deg == pytest.approx(exact, abs=1e-9)
    assert route.direct_tilt_deg == pytest.approx(exact, abs=1e-9)
    assert "follow the" in route.legs[0].describe()


def test_route_waypoints_lie_on_the_followed_band_between_the_endpoints(cubic_map) -> None:
    """The landmarks an operator expects must really be on the road.

    Each waypoint has to satisfy the zone law for the followed band and lie inside
    the arc, which is what makes it a confirmation that the tilt is tracking.
    """

    route = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=60.0)
    leg = route.legs[0]
    assert leg.waypoint_indices, "[112] and its neighbours lie on this arc"
    band = np.asarray(leg.band_indices, dtype=np.int64)
    start = cubic_map.zone_axis_for_direction(leg.start_indices)
    end = cubic_map.zone_axis_for_direction(leg.end_indices)
    assert start is not None and end is not None
    span = math.degrees(
        math.acos(abs(float(np.dot(start.direction_map, end.direction_map))))
    )
    for indices in leg.waypoint_indices:
        assert int(band @ np.asarray(indices, dtype=np.int64)) == 0
        axis = cubic_map.zone_axis_for_direction(indices)
        assert axis is not None
        from_start = math.degrees(
            math.acos(min(1.0, abs(float(np.dot(start.direction_map, axis.direction_map)))))
        )
        assert 0.0 < from_start < span + 1e-9


def test_a_long_hop_is_split_and_the_split_costs_no_extra_travel(cubic_map) -> None:
    """Capping the leg length must not lengthen a route that lies on one band.

    Splitting an arc at points on the same great circle leaves the total arc
    length unchanged, so the multi-hop route costs nothing in travel and buys the
    chance to re-index part way.
    """

    single = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=60.0)
    split = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=30.0)
    assert split.reachable
    assert split.hop_count > single.hop_count
    assert split.total_tilt_deg == pytest.approx(single.total_tilt_deg, abs=1e-6)
    assert split.total_tilt_deg == pytest.approx(split.direct_tilt_deg, abs=1e-6)


def test_an_impossible_leg_budget_is_reported_rather_than_raised(cubic_map) -> None:
    """An unreachable target returns a describable answer, not an exception.

    The caller asked a reasonable question and the honest reply is "not with that
    budget", together with what to change.
    """

    route = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=1.0)
    assert not route.reachable
    assert route.legs[0].band_indices is None
    text = route.describe()
    assert "No route" in text
    assert "zone-axis index bound" in text
    assert route.to_json_dict()["reachable"] is False


def test_routing_between_identical_axes_is_a_no_op(cubic_map) -> None:
    route = cubic_map.route_to([0, 0, 1], [0, 0, 1])
    assert route.reachable
    assert route.legs == ()
    assert route.total_tilt_deg == 0.0
    assert "no tilt required" in route.describe()


# ---------------------------------------------------------------------------
# Hexagonal notation and geometry
# ---------------------------------------------------------------------------


def test_the_hexagonal_map_names_itself_in_four_index_notation(hexagonal_map) -> None:
    """A hexagonal direction is written ``[uvtw]`` and a plane ``(hkil)``.

    The three-index form hides the symmetry — the members of one family do not
    look like permutations of each other — which is why the literature does not
    use it for hexagonal crystals.
    """

    text = hexagonal_map.describe()
    assert "[0001]" in text
    basal = hexagonal_map.zone_axis_for_direction([0, 0, 1])
    assert basal is not None
    assert "[0001]" in basal.describe()
    # [100] in three indices is [2-1-10] in four: u = 1, v = 0, t = -(u+v) = -1,
    # and the four-index form is scaled so the first three sum to zero.
    prismatic = hexagonal_map.zone_axis_for_direction([1, 0, 0])
    assert prismatic is not None
    assert prismatic.indices == (1, 0, 0)
    assert "[2 -1 -1 0]" in prismatic.describe()
    # And a plane: (100) in three indices is (10-10) in four.
    band = hexagonal_map.band_for_plane([1, 0, 0])
    if band is not None:
        assert "(1 0 -1 0)" in band.describe()


def test_the_basal_and_prismatic_hexagonal_axes_are_exactly_ninety_degrees_apart(
    hexagonal_map,
) -> None:
    """[0001] and [2-1-10] are perpendicular in any hexagonal lattice.

    The c axis is perpendicular to every direction in the basal plane whatever the
    ``c/a`` ratio, so this angle is a lattice-independent right angle — the
    cleanest available check that the direct basis and the map frame agree.
    """

    route = hexagonal_map.route_to([0, 0, 1], [1, 0, 0], max_leg_deg=95.0)
    assert route.reachable
    assert route.direct_tilt_deg == pytest.approx(90.0, abs=1e-9)
    assert route.total_tilt_deg == pytest.approx(90.0, abs=1e-6)
    assert "[0001]" in route.describe()


def test_the_strongest_hexagonal_band_is_the_pyramidal_one(hexagonal_map) -> None:
    """{10-11} outranks {11-20} and {0002} in hcp, as the reflection tables have it."""

    strongest = hexagonal_map.bands[0]
    assert sorted(abs(value) for value in strongest.indices) == [0, 1, 1]
    assert strongest.d_spacing_angstrom == pytest.approx(2.4589, abs=1e-3)


# ---------------------------------------------------------------------------
# Reporting surfaces and guards
# ---------------------------------------------------------------------------


def test_describe_states_the_conventions_and_the_limits(cubic_map) -> None:
    text = cubic_map.describe()
    assert "Kossel cones" in text
    assert "kinematic" in text
    assert "smallest spacing, not the strongest" in text
    assert f"{cubic_map.band_count} bands" in text


def test_the_json_payload_matches_the_map(cubic_map) -> None:
    payload = cubic_map.to_json_dict()
    assert payload["schema"] == "pytex.stereographic_kikuchi_map.v1"
    assert len(payload["bands"]) == cubic_map.band_count
    assert len(payload["zone_axes"]) == cubic_map.zone_axis_count
    assert payload["centre"] == [0, 0, 1]
    first = payload["bands"][0]
    assert first["angular_width_deg"] == pytest.approx(2.0 * first["bragg_angle_deg"])
    for entry in payload["zone_axes"]:
        assert entry["order"] == len(entry["bands"])


def test_a_band_contains_the_zone_axes_recorded_on_it(cubic_map) -> None:
    """Membership and the zone-axis list must be the same statement.

    A zone axis on a band's centre line is at 90 degrees to the plane normal, so
    it is inside the band by any positive Bragg angle; a direction along the
    normal is as far outside as it is possible to be.
    """

    # Bands are ordered by intensity, so position zero is the strongest band.
    position = 0
    band = cubic_map.bands[position]
    on_it = [axis for axis in cubic_map.zone_axes if position in axis.band_indices]
    assert on_it, "the strongest band must carry some zone axes"
    for axis in on_it:
        assert band.contains_direction(axis.direction_map)
    assert not band.contains_direction(band.normal_map)


def test_traces_are_bounded_and_the_centre_line_passes_through_its_zone_axes(
    cubic_map,
) -> None:
    """Every trace lies inside the unit disc, and a centre line hits its own axes.

    Boundedness is the property that makes the stereographic projection usable as
    an atlas: a hemisphere maps into the unit circle, whereas the gnomonic
    projection of the same hemisphere is unbounded. And a band's centre line must
    pass through every zone axis recorded on it, since both statements are the
    zone law seen once as geometry and once as arithmetic.
    """

    position = 0
    band = cubic_map.bands[position]
    centre = band.centre_trace(samples=721)
    lower, upper = band.edge_traces(samples=721)
    for trace in (centre, lower, upper):
        assert np.all(np.linalg.norm(trace, axis=1) <= 1.0 + 1e-9)
    assert lower.shape == upper.shape == centre.shape

    for axis in cubic_map.zone_axes:
        if position not in axis.band_indices:
            continue
        point = axis.projected()
        gap = float(np.min(np.linalg.norm(centre - point[None, :], axis=1)))
        # 721 samples over a full turn is a half-degree step; the trace is smooth,
        # so the nearest sample is within that of any point on it.
        assert gap < 5e-3, (axis.indices, gap)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_index": 0}, "max_index must be strictly positive"),
        ({"zone_axis_max_index": 0}, "zone_axis_max_index must be strictly positive"),
        ({"min_relative_intensity": 1.5}, "min_relative_intensity must lie in"),
        ({"max_bands": 0}, "max_bands must be strictly positive"),
        ({"min_d_spacing_angstrom": 0.0}, "min_d_spacing_angstrom must be strictly positive"),
        ({"min_zone_axis_order": 1}, "a crossing needs two bands"),
        ({"max_polar_angle_deg": 0.0}, "max_polar_angle_deg must lie in"),
        ({"horizontal_direction": (0, 0, 1)}, "parallel to the map centre"),
    ],
)
def test_invalid_map_requests_are_rejected_with_an_explanation(nickel, kwargs, message) -> None:
    with pytest.raises(ValueError, match=message):
        compute_kikuchi_map(nickel, **kwargs)


def test_routing_to_a_direction_that_is_not_on_the_map_says_what_to_change(cubic_map) -> None:
    with pytest.raises(ValueError, match="zone_axis_max_index"):
        plan_kikuchi_route(cubic_map, [0, 0, 1], [7, 9, 11])


def test_a_zero_leg_budget_is_rejected(cubic_map) -> None:
    with pytest.raises(ValueError, match="max_leg_deg must be strictly positive"):
        cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=0.0)


def test_construction_guards_on_the_records(cubic_map) -> None:
    band = cubic_map.bands[0]
    with pytest.raises(ValueError, match=r"must lie in \(0, 90\)"):
        KikuchiMapBand(
            plane=band.plane,
            normal_map=band.normal_map,
            bragg_angle_deg=0.0,
            d_spacing_angstrom=band.d_spacing_angstrom,
            relative_intensity=band.relative_intensity,
            family_multiplicity=band.family_multiplicity,
        )
    axis = cubic_map.zone_axes[0]
    with pytest.raises(ValueError, match="at least two bands"):
        KikuchiMapZoneAxis(
            phase=cubic_map.phase,
            indices=axis.indices,
            direction_map=axis.direction_map,
            band_indices=(0,),
            polar_angle_deg=axis.polar_angle_deg,
        )
    with pytest.raises(ValueError, match="must be orthonormal"):
        StereographicKikuchiMap(
            phase=cubic_map.phase,
            beam_energy_kev=cubic_map.beam_energy_kev,
            wavelength_angstrom=cubic_map.wavelength_angstrom,
            centre_indices=cubic_map.centre_indices,
            horizontal_indices=cubic_map.horizontal_indices,
            view_matrix=np.full((3, 3), 0.5),
            bands=cubic_map.bands,
            zone_axes=cubic_map.zone_axes,
        )
