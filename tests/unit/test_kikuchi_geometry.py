"""Kikuchi band geometry and the gnomonic projection.

Expected values here are derived analytically, not copied from a prior run of
the code under test. Three independent anchors are used:

1. **The defining property of the gnomonic projection** — a great circle maps to
   a straight line — which the band centre trace must satisfy exactly.
2. **Bragg's law**, ``sin(theta_B) = lambda / 2d``, evaluated from cell
   parameters rather than from the simulator.
3. **Pure geometry** — a direction 45 degrees from the detector normal must land
   at gnomonic radius ``tan(45 deg) = 1`` — which pins the projection scale
   without reference to any crystallography.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import (
    CRYSTAL_FRAME,
    DETECTOR_FRAME,
    LABORATORY_FRAME,
    SPECIMEN_FRAME,
)
from pytex.core.lattice import Lattice, Phase
from pytex.core.miller import MillerPlane
from pytex.core.orientation import Orientation
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.kikuchi import (
    GnomonicProjection,
    KikuchiBand,
    simulate_kikuchi_pattern,
)
from pytex.diffraction.models import DiffractionGeometry

#: Nickel, a = 3.52387 A, FCC. Chosen because its {111} and {002} spacings and
#: hence its Kikuchi band widths at 20 kV are standard textbook values.
_NICKEL_A = 3.52387


def _cubic_phase(*, lattice_parameter: float = _NICKEL_A, symbol: str = "Fm-3m") -> Phase:
    lattice = Lattice(
        a=lattice_parameter,
        b=lattice_parameter,
        c=lattice_parameter,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=CRYSTAL_FRAME,
    )
    return Phase(
        name="test-cubic",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
        space_group_symbol=symbol,
        space_group_number=225,
    )


def _geometry(*, beam_energy_kev: float = 20.0, tilt: tuple[float, float, float] = (0.0, 0.0, 0.0)):
    return DiffractionGeometry(
        detector_frame=DETECTOR_FRAME,
        specimen_frame=SPECIMEN_FRAME,
        laboratory_frame=LABORATORY_FRAME,
        beam_energy_kev=beam_energy_kev,
        camera_length_mm=15.0,
        pattern_center=np.array([0.5, 0.5, 0.6]),
        detector_pixel_size_um=(50.0, 50.0),
        detector_shape=(480, 640),
        tilt_degrees=tilt,
    )


def _identity_orientation(phase: Phase) -> Orientation:
    return Orientation.from_euler(0.0, 0.0, 0.0, specimen_frame=SPECIMEN_FRAME, phase=phase)


# --------------------------------------------------------------------------- #
# GnomonicProjection
# --------------------------------------------------------------------------- #


def test_direction_at_45_degrees_projects_to_unit_gnomonic_radius() -> None:
    """Pure geometry: gnomonic radius is the tangent of the angle from the normal.

    Independent of any crystallography, so it pins the projection scale on its
    own.
    """

    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    basis = geometry.detector_basis_lab
    for angle_deg in (0.0, 15.0, 30.0, 45.0, 60.0):
        angle = np.deg2rad(angle_deg)
        direction = np.cos(angle) * basis[:, 2] + np.sin(angle) * basis[:, 0]
        coordinates, valid = projection.project_directions(direction[None, :])
        assert bool(valid[0])
        assert coordinates[0, 0] == pytest.approx(np.tan(angle), abs=1e-12)
        assert coordinates[0, 1] == pytest.approx(0.0, abs=1e-12)


def test_gnomonic_projection_round_trips_through_directions_and_pixels() -> None:
    projection = GnomonicProjection(geometry=_geometry())
    rng = np.random.default_rng(20260805)
    coordinates = rng.uniform(-0.8, 0.8, size=(64, 2))

    recovered, valid = projection.project_directions(projection.unproject(coordinates))
    assert bool(np.all(valid))
    assert np.allclose(recovered, coordinates, atol=1e-12)

    pixels = projection.to_detector_px(coordinates)
    assert np.allclose(projection.from_detector_px(pixels), coordinates, atol=1e-12)


def test_directions_away_from_the_detector_are_flagged_invalid_not_projected() -> None:
    """A backward direction must be reported as invalid rather than mirrored."""

    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    basis = geometry.detector_basis_lab
    forward = basis[:, 2]
    coordinates, valid = projection.project_directions(np.stack([forward, -forward]))
    assert bool(valid[0])
    assert not bool(valid[1])
    assert np.all(np.isnan(coordinates[1]))


def test_pattern_center_projects_to_the_gnomonic_origin_and_its_pixel() -> None:
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    coordinates, valid = projection.project_directions(geometry.detector_basis_lab[:, 2][None, :])
    assert bool(valid[0])
    assert np.allclose(coordinates[0], (0.0, 0.0), atol=1e-12)
    assert np.allclose(
        projection.to_detector_px(coordinates[:1])[0],
        geometry.pattern_center_px,
        atol=1e-9,
    )


def test_detector_containment_matches_the_physical_sensor() -> None:
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    corners = projection.detector_corner_coordinates()
    assert bool(np.all(projection.contains(corners)))
    # Ten gnomonic units is far beyond any detector at this distance.
    assert not bool(projection.contains(np.array([[10.0, 10.0]]))[0])


# --------------------------------------------------------------------------- #
# Band geometry
# --------------------------------------------------------------------------- #


def test_band_angular_width_is_twice_the_bragg_angle_from_braggs_law() -> None:
    """Independent anchor: d and theta_B computed from the cell, not the simulator."""

    phase = _cubic_phase()
    geometry = _geometry()
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)
    wavelength = geometry.electron_wavelength_angstrom

    band = pattern.band_for_plane((1, 1, 1))
    assert band is not None
    expected_d = _NICKEL_A / np.sqrt(3.0)
    assert band.d_spacing_angstrom == pytest.approx(expected_d, rel=1e-12)
    expected_bragg = np.arcsin(wavelength / (2.0 * expected_d))
    assert band.bragg_angle_rad == pytest.approx(expected_bragg, rel=1e-12)
    assert band.angular_width_rad == pytest.approx(2.0 * expected_bragg, rel=1e-12)

    band = pattern.band_for_plane((0, 0, 2))
    assert band is not None
    assert band.d_spacing_angstrom == pytest.approx(_NICKEL_A / 2.0, rel=1e-12)


def test_nickel_111_band_width_matches_the_published_20kv_value() -> None:
    """A {111} Ni band is about 2.42 degrees wide at 20 kV.

    Checked against the value implied by the standard 20 kV electron wavelength
    of 0.0859 A and the tabulated Ni lattice parameter, computed here by hand.
    """

    phase = _cubic_phase()
    geometry = _geometry(beam_energy_kev=20.0)
    assert geometry.electron_wavelength_angstrom == pytest.approx(0.0859, abs=5e-5)

    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)
    band = pattern.band_for_plane((1, 1, 1))
    assert band is not None
    assert np.rad2deg(band.angular_width_rad) == pytest.approx(2.42, abs=0.01)


def test_wider_bands_come_from_larger_d_spacings() -> None:
    """Band width is a direct measurement of the lattice, so the order is fixed."""

    phase = _cubic_phase()
    pattern = simulate_kikuchi_pattern(
        _geometry(), phase, _identity_orientation(phase), max_index=2
    )
    ordered = pattern.widest_bands(len(pattern.bands))
    spacings = [band.d_spacing_angstrom for band in ordered]
    widths = [band.angular_width_rad for band in ordered]
    assert spacings == sorted(spacings, reverse=True)
    assert widths == sorted(widths)


@pytest.mark.parametrize("tilt", [(0.0, 0.0, 0.0), (10.0, -7.0, 3.0), (25.0, 0.0, 0.0)])
def test_band_centre_traces_are_exactly_straight_in_gnomonic_coordinates(
    tilt: tuple[float, float, float],
) -> None:
    """The defining property of the gnomonic projection, at arbitrary detector tilt.

    A lattice-plane trace is a great circle, and great circles map to straight
    lines. This is what makes gnomonic coordinates the right frame for Kikuchi
    band detection, and it must hold exactly regardless of how the detector is
    tilted.
    """

    phase = _cubic_phase()
    geometry = _geometry(tilt=tilt)
    projection = GnomonicProjection(geometry=geometry)
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)

    checked = 0
    for band in pattern.bands:
        trace = band.center_trace(projection, samples=721)
        if trace.shape[0] == 0:
            continue
        a, b, c = band.center_line_coefficients(projection)
        residual = np.abs(a * trace[:, 0] + b * trace[:, 1] + c)
        # The gnomonic coordinate diverges for near-grazing rays, so the line
        # residual must be judged relative to the coordinate magnitude rather
        # than against a fixed absolute floor.
        scale = 1.0 + np.abs(trace).max(axis=1)
        assert (residual / scale).max() < 1e-14
        checked += 1
    assert checked > 0


def test_band_edges_lie_exactly_on_their_kossel_cones() -> None:
    """Edges are conics, not lines: they are checked on the cone, not on a chord."""

    phase = _cubic_phase()
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)

    checked = 0
    for band in pattern.bands:
        for sign, edge in zip((1.0, -1.0), band.edge_traces(projection, samples=361), strict=True):
            if edge.shape[0] == 0:
                continue
            directions = projection.unproject(edge)
            dots = directions @ band.plane_normal_lab
            assert np.abs(dots - sign * np.sin(band.bragg_angle_rad)).max() < 1e-12
            checked += 1
    assert checked > 0


def test_band_edges_are_not_straight_lines() -> None:
    """Guard against a small-angle approximation creeping in.

    The Kossel cones are not great circles, so their gnomonic traces are conics.
    A test suite that only checked "edges are parallel to the centre" would pass
    for a wrong, linearized implementation; this one would not.
    """

    phase = _cubic_phase()
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    # A deliberately large Bragg angle makes the curvature unmistakable, and is
    # constructed directly rather than through a simulation so the test states
    # its own premise.
    band = KikuchiBand(
        plane=MillerPlane(indices=np.array([1, 1, 1]), phase=phase),
        plane_normal_lab=np.array([0.0, 0.6, 0.8]),
        bragg_angle_rad=np.deg2rad(20.0),
        d_spacing_angstrom=2.0,
        intensity=1.0,
    )
    edge, _ = band.edge_traces(projection, samples=1001)
    assert edge.shape[0] > 100
    # Fit a straight line to the edge and require a real misfit.
    fit = np.polyfit(edge[:, 0], edge[:, 1], 1)
    residual = np.abs(np.polyval(fit, edge[:, 0]) - edge[:, 1]).max()
    assert residual > 1e-3


def test_band_width_is_symmetric_about_the_centre_line_in_angle() -> None:
    """Both edges sit exactly theta_B from the plane, on opposite sides."""

    phase = _cubic_phase()
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)
    band = pattern.widest_bands(1)[0]

    lower, upper = band.edge_traces(projection, samples=181)
    for sign, edge in ((1.0, lower), (-1.0, upper)):
        if edge.shape[0] == 0:
            continue
        directions = projection.unproject(edge)
        angles = np.arcsin(np.clip(directions @ band.plane_normal_lab, -1.0, 1.0))
        assert np.allclose(angles, sign * band.bragg_angle_rad, atol=1e-12)


def test_band_width_at_pattern_center_matches_the_edge_traces() -> None:
    """The analytic width agrees with the edge geometry it summarizes."""

    phase = _cubic_phase()
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)

    for band in pattern.bands[:6]:
        a, b, _ = band.center_line_coefficients(projection)
        in_plane_norm = float(np.hypot(a, b))
        if in_plane_norm < 1e-9:
            continue
        width = band.width_at_pattern_center(projection)
        if not np.isfinite(width):
            continue
        # Walk the perpendicular through the pattern centre and locate the two
        # points at exactly the Bragg angle to the plane.
        unit = np.array([a, b]) / in_plane_norm
        offsets = np.linspace(-4.0, 4.0, 400_001)
        directions = projection.unproject(offsets[:, None] * unit[None, :])
        dots = directions @ band.plane_normal_lab
        target = float(np.sin(band.bragg_angle_rad))
        crossings = [
            offsets[int(np.argmin(np.abs(dots - value)))] for value in (target, -target)
        ]
        assert abs(abs(crossings[0] - crossings[1]) - width) / width < 1e-3


# --------------------------------------------------------------------------- #
# Zone axes
# --------------------------------------------------------------------------- #


def test_zone_axis_lies_on_the_centre_line_of_every_band_containing_it() -> None:
    """The defining property of a zone axis, checked against the zone law."""

    phase = _cubic_phase()
    geometry = _geometry()
    projection = GnomonicProjection(geometry=geometry)
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)
    assert pattern.zone_axes

    for axis in pattern.zone_axes:
        matched = 0
        for band in pattern.bands:
            if int(np.asarray(band.plane.indices) @ axis.indices) != 0:
                continue
            a, b, c = band.center_line_coefficients(projection)
            assert abs(a * axis.coordinates[0] + b * axis.coordinates[1] + c) < 1e-9
            matched += 1
        assert matched == axis.band_count
        assert matched >= 2


def test_cube_axis_projects_to_45_degrees_for_the_identity_orientation() -> None:
    """A closed-form check that the crystal-to-laboratory chain is wired correctly.

    With the identity orientation, an untilted detector, and a cubic phase, the
    ``[011]`` zone axis lies 45 degrees from the detector normal and must project
    to gnomonic radius ``tan(45 deg) = 1``. A transposed or misordered rotation
    anywhere in the crystal to specimen to laboratory chain would move it.
    """

    phase = _cubic_phase()
    geometry = _geometry()
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)
    axes = {tuple(int(v) for v in axis.indices): axis for axis in pattern.zone_axes}
    assert (0, 1, 1) in axes
    radius = float(np.hypot(*axes[(0, 1, 1)].coordinates))
    assert radius == pytest.approx(1.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# Simulation contract
# --------------------------------------------------------------------------- #


def test_fcc_centring_absences_are_applied_to_the_band_list() -> None:
    """Only all-even or all-odd ``(hkl)`` may produce a band in an FCC lattice."""

    phase = _cubic_phase(symbol="Fm-3m")
    pattern = simulate_kikuchi_pattern(
        _geometry(), phase, _identity_orientation(phase), max_index=3
    )
    for band in pattern.bands:
        h, k, ell = (int(value) for value in band.plane.indices)
        assert h % 2 == k % 2 == ell % 2, (h, k, ell)


def test_antipodal_planes_produce_a_single_band() -> None:
    """``(hkl)`` and ``(-h-k-l)`` are the same plane, so they are one band."""

    phase = _cubic_phase()
    pattern = simulate_kikuchi_pattern(
        _geometry(), phase, _identity_orientation(phase), max_index=2
    )
    keys = {tuple(int(value) for value in band.plane.indices) for band in pattern.bands}
    for key in keys:
        assert tuple(-value for value in key) not in keys


def test_simulation_rejects_inconsistent_frames_and_phases() -> None:
    phase = _cubic_phase()
    other = _cubic_phase(lattice_parameter=4.0)
    geometry = _geometry()
    with pytest.raises(ValueError, match=r"orientation\.phase must match phase"):
        simulate_kikuchi_pattern(geometry, phase, _identity_orientation(other))
    with pytest.raises(ValueError, match="max_index must be strictly positive"):
        simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=0)


def test_higher_beam_energy_narrows_every_band() -> None:
    """Shorter wavelength means a smaller Bragg angle at fixed d-spacing."""

    phase = _cubic_phase()
    orientation = _identity_orientation(phase)
    low = simulate_kikuchi_pattern(_geometry(beam_energy_kev=10.0), phase, orientation, max_index=2)
    high = simulate_kikuchi_pattern(
        _geometry(beam_energy_kev=30.0), phase, orientation, max_index=2
    )
    for band in low.bands:
        counterpart = high.band_for_plane(band.plane.indices)
        assert counterpart is not None
        assert counterpart.angular_width_rad < band.angular_width_rad


def test_describe_states_the_conventions_and_the_kinematic_limitation() -> None:
    """`describe()` is a tested surface, per the explainable-results doctrine."""

    phase = _cubic_phase()
    geometry = _geometry()
    pattern = simulate_kikuchi_pattern(geometry, phase, _identity_orientation(phase), max_index=2)

    text = pattern.describe()
    assert "Bunge Euler" in text
    assert "gnomonic" in text
    assert "kinematic" in text.lower()
    assert f"{geometry.beam_energy_kev:.1f} keV" in text
    widest = pattern.widest_bands(1)[0]
    assert f"{widest.d_spacing_angstrom:.4f}" in text

    band_text = widest.describe()
    assert f"{np.rad2deg(widest.angular_width_rad):.4f}" in band_text
    assert "Kossel" in band_text

    axis_text = pattern.zone_axes[0].describe()
    assert "Zone axis" in axis_text
    assert f"{pattern.zone_axes[0].band_count}" in axis_text
