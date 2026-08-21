"""The EBSD camera geometry, and the pattern the workbench simulates on it.

Expected values are derived from the stated convention and from Bragg's law, not
copied from a run of the code under test. Three anchors carry the file:

1. **The convention itself.** With the beam along the laboratory z axis, the
   stage tilting the specimen normal towards the camera, and the camera axis at
   elevation ``eps`` above the plane perpendicular to the beam, the angle
   between the specimen normal and the camera axis is ``90 - (sigma - eps)``.
   The specimen normal therefore projects at gnomonic radius
   ``tan(90 - sigma + eps)`` — pure geometry, independent of any crystallography.
2. **Bragg's law**, ``sin(theta_B) = lambda / 2d``, evaluated from the lattice
   parameter rather than from the simulator, which fixes every band width.
3. **The definition of a zone axis as a direction**, which is why the axes a
   pattern reports must be coprime.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app.errors import InvalidInputError
from pytex.app.registry import REGISTRY
from pytex.diffraction.kikuchi import GnomonicProjection
from pytex.diffraction.models import DiffractionGeometry

_OPERATION = "ebsd.simulate_kikuchi_pattern"

#: Nickel, a = 3.52387 A, fcc. Its {111} and {200} spacings are standard values.
_NICKEL_A = 3.52387


def _defaults(**overrides: object) -> dict[str, object]:
    """The request the panel sends when nothing has been touched."""

    spec = next(entry for entry in REGISTRY.operations() if entry.id == _OPERATION)
    request: dict[str, object] = {
        parameter.name: parameter.default
        for parameter in spec.parameters
        if parameter.default is not None
    }
    request["phase"] = {"builtin": "ni_fcc"}
    request.update(overrides)
    return request


# --------------------------------------------------------------------------- #
# DiffractionGeometry.for_ebsd
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("tilt", "elevation"),
    [(70.0, 0.0), (70.0, 10.0), (60.0, 0.0), (30.0, 0.0), (75.0, -5.0)],
)
def test_the_specimen_normal_projects_where_the_convention_says(
    tilt: float, elevation: float
) -> None:
    """``tan(90 - sigma + eps)``, from the stated frame alone.

    This is the one number that pins the whole geometry: it fixes the sense of
    the stage tilt, the sense of the elevation, and the direction the camera
    faces, all at once. A sign error in any of them changes it.
    """

    geometry = DiffractionGeometry.for_ebsd(
        sample_tilt_deg=tilt, detector_elevation_deg=elevation
    )
    projection = GnomonicProjection(geometry=geometry)
    normal_lab = np.asarray(geometry.specimen_vectors_to_lab(np.array([[0.0, 0.0, 1.0]])))
    coordinates, valid = projection.project_directions(normal_lab)

    expected = math.tan(math.radians(90.0 - tilt + elevation))
    assert bool(valid[0])
    assert float(np.hypot(*coordinates[0])) == pytest.approx(expected, rel=1e-12)


def test_an_untilted_specimen_shows_its_normal_edge_on_to_the_camera() -> None:
    """The limit of the same formula: at zero tilt and zero elevation it diverges.

    The specimen normal then lies *in* the screen plane, which is the geometric
    reason a scan is run tilted at all — untilted, the surface the camera has to
    look at is edge-on to it.
    """

    geometry = DiffractionGeometry.for_ebsd(sample_tilt_deg=0.0, detector_elevation_deg=0.0)
    _, valid = GnomonicProjection(geometry=geometry).project_directions(
        np.asarray(geometry.specimen_vectors_to_lab(np.array([[0.0, 0.0, 1.0]])))
    )
    assert not bool(valid[0])


def test_the_specimen_normal_faces_the_beam_rather_than_following_it() -> None:
    """A surface the beam strikes has its normal pointing back up the column.

    The opposite sign would put the specimen normal into the specimen, and every
    pattern would be simulated from the wrong hemisphere while still looking
    like a plausible band network.
    """

    geometry = DiffractionGeometry.for_ebsd(sample_tilt_deg=70.0)
    normal_lab = np.asarray(geometry.specimen_vectors_to_lab(np.array([[0.0, 0.0, 1.0]])))[0]
    assert float(normal_lab @ geometry.beam_direction_lab) < 0.0
    # And it tips towards the camera, which sits on -y.
    assert float(normal_lab[1]) < 0.0
    assert float(normal_lab @ geometry.beam_direction_lab) == pytest.approx(
        -math.cos(math.radians(70.0)), rel=1e-12
    )


def test_the_beam_has_no_image_on_an_unelevated_screen() -> None:
    """At zero elevation the beam is parallel to the screen — as it must be.

    An EBSD camera never sees the incident beam. A geometry that projected it
    onto the screen would be one where the camera faces the wrong way.
    """

    geometry = DiffractionGeometry.for_ebsd(detector_elevation_deg=0.0)
    projection = GnomonicProjection(geometry=geometry)
    _, valid = projection.project_directions(geometry.beam_direction_lab.reshape(1, 3))
    assert not bool(valid[0])

    raised = DiffractionGeometry.for_ebsd(detector_elevation_deg=12.0)
    coordinates, valid = GnomonicProjection(geometry=raised).project_directions(
        raised.beam_direction_lab.reshape(1, 3)
    )
    assert bool(valid[0])
    assert float(np.hypot(*coordinates[0])) == pytest.approx(
        math.tan(math.radians(90.0 - 12.0)), rel=1e-12
    )


def test_the_camera_distance_is_the_pattern_centre_fraction_of_the_screen_width() -> None:
    """``z*`` is a fraction of the screen *width*, which is what calibration reports."""

    geometry = DiffractionGeometry.for_ebsd(
        pattern_center=(0.5, 0.5, 0.65),
        detector_shape=(480, 640),
        detector_pixel_size_um=(50.0, 50.0),
    )
    assert geometry.camera_length_mm == pytest.approx(0.65 * 640 * 0.05)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"sample_tilt_deg": 90.0}, "sample_tilt_deg"),
        ({"sample_tilt_deg": -1.0}, "sample_tilt_deg"),
        ({"detector_elevation_deg": 90.0}, "detector_elevation_deg"),
        ({"pattern_center": (0.5, 0.5, 0.0)}, "z component"),
    ],
)
def test_an_impossible_camera_is_refused_at_construction(
    kwargs: dict[str, object], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        DiffractionGeometry.for_ebsd(**kwargs)  # type: ignore[arg-type]


def test_a_zone_axis_is_a_direction_so_the_reported_axes_are_coprime() -> None:
    """``[002]`` is the same axis as ``[001]``: listing both counts one hub twice."""

    from pytex.app.phases import phase_from_request
    from pytex.core.frame_catalog import SPECIMEN_FRAME
    from pytex.core.orientation import Orientation
    from pytex.diffraction.kikuchi import simulate_kikuchi_pattern

    _, phase = phase_from_request({"builtin": "ni_fcc"})
    geometry = DiffractionGeometry.for_ebsd()
    orientation = Orientation.from_euler(
        0.0, 0.0, 0.0, degrees=True, specimen_frame=SPECIMEN_FRAME, phase=phase
    )
    pattern = simulate_kikuchi_pattern(geometry, phase, orientation, max_index=3, max_bands=20)

    assert pattern.zone_axes
    seen = set()
    for axis in pattern.zone_axes:
        indices = tuple(int(value) for value in axis.indices)
        assert math.gcd(*(abs(value) for value in indices)) == 1, indices
        assert indices not in seen
        seen.add(indices)


# --------------------------------------------------------------------------- #
# The workbench operation
# --------------------------------------------------------------------------- #


def test_every_band_width_is_twice_its_bragg_angle_from_the_lattice() -> None:
    """Bragg's law evaluated from ``a = 3.52387 A``, not from the simulator."""

    result = REGISTRY.call(_OPERATION, _defaults())
    wavelength = float(result["data"]["wavelength_angstrom"])
    assert result["data"]["bands"]

    for band in result["data"]["bands"]:
        h, k, ell = band["hkl"]
        spacing = _NICKEL_A / math.sqrt(h * h + k * k + ell * ell)
        assert band["d_angstrom"] == pytest.approx(spacing, rel=1e-9)
        expected = 2.0 * math.degrees(math.asin(wavelength / (2.0 * spacing)))
        assert band["width_deg"] == pytest.approx(expected, rel=1e-9)


def test_lowering_the_voltage_widens_every_band_and_moves_none_of_them() -> None:
    """The voltage is a band-width control and nothing else geometric.

    The widths must scale exactly as ``arcsin(lambda/2d)`` does, and the band
    *set* must be unchanged, because neither the lattice nor the orientation nor
    the camera moved.
    """

    at_20 = REGISTRY.call(_OPERATION, _defaults(beam_energy_kev=20.0))["data"]
    at_10 = REGISTRY.call(_OPERATION, _defaults(beam_energy_kev=10.0))["data"]

    wide = {tuple(band["hkl"]): band for band in at_10["bands"]}
    assert at_10["wavelength_angstrom"] > at_20["wavelength_angstrom"]
    for band in at_20["bands"]:
        other = wide[tuple(band["hkl"])]
        expected = 2.0 * math.degrees(
            math.asin(at_10["wavelength_angstrom"] / (2.0 * band["d_angstrom"]))
        )
        assert other["width_deg"] == pytest.approx(expected, rel=1e-9)
        assert other["width_deg"] > band["width_deg"]


def test_the_summary_states_the_geometry_that_was_simulated() -> None:
    """The specimen-normal radius is the check a reader can make on the geometry."""

    result = REGISTRY.call(_OPERATION, _defaults(sample_tilt_deg=70.0))
    assert result["data"]["specimen_normal_gnomonic_radius"] == pytest.approx(
        math.tan(math.radians(20.0))
    )
    assert "70" in result["summary"]
    assert "kinematic" in result["summary"]


def test_nothing_is_drawn_further_than_one_frame_outside_the_screen() -> None:
    """Traces are clipped, so a Kossel conic cannot draw a chord across the picture."""

    data = REGISTRY.call(_OPERATION, _defaults())["data"]
    width = float(data["width_px"])
    height = float(data["height_px"])
    margin = max(width, height)

    drawn = 0
    for band in data["bands"]:
        runs = [*band["centre"], *band["edges"][0], *band["edges"][1]]
        for run in runs:
            assert len(run) >= 2
            drawn += 1
            for x, y in run:
                assert -margin <= x <= width + margin
                assert -margin <= y <= height + margin
    assert drawn > 0


def test_bands_arrive_strongest_first_so_the_drawing_can_stop_anywhere() -> None:
    data = REGISTRY.call(_OPERATION, _defaults())["data"]
    intensities = [band["intensity"] for band in data["bands"]]
    assert intensities == sorted(intensities, reverse=True)


def test_a_field_of_view_that_sees_no_band_is_refused_with_a_hint() -> None:
    """A camera pushed far enough away sees a patch of sky with nothing in it.

    The failure must be a stated one. Returning an empty pattern would look like
    a phase with no bands rather than a geometry with no field of view.
    """

    with pytest.raises(InvalidInputError) as error:
        REGISTRY.call(
            _OPERATION,
            _defaults(
                detector_distance=3.0,
                detector_width_px=64,
                detector_height_px=64,
                pixel_size_um=1.0,
                max_bands=1,
                max_index=1,
            ),
        )
    assert "screen" in str(error.value)
