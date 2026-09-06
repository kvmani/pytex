"""The Kikuchi pattern an EBSD camera would record, simulated from the lattice.

Why this sits in the EBSD workspace and not beside the TEM simulator
--------------------------------------------------------------------
The physics is shared — a band is a pair of Kossel cones about a plane normal,
and :mod:`pytex.diffraction.kikuchi` computes it once for both techniques. The
*geometry* is not shared at all. A TEM Kikuchi overlay is drawn on a plate whose
detector faces the beam, at 200 kV, on a specimen the beam goes through. An EBSD
pattern comes off a specimen tilted seventy degrees, onto a screen the beam
never reaches, at twenty kilovolts, and it is configured by a microscopist in
terms of stage tilt, camera elevation and a pattern centre — not camera length.
An operation that asked for those in TEM terms would be answering a question
nobody at an EBSD stage is asking.

So the panel is stated in EBSD terms, and every one of them is turned into the
one canonical :class:`~pytex.diffraction.models.DiffractionGeometry` by
:meth:`~pytex.diffraction.models.DiffractionGeometry.for_ebsd`, which owns the
frame convention. Nothing here defines a frame of its own.

What the simulation is, and what it is not
------------------------------------------
Geometric and kinematic. Band positions and widths are exact for the stated
geometry — a band's angular width is exactly twice its Bragg angle — and those
positions are what orientation determination uses. Intensities are a kinematic
``|F|^2`` proxy: there is no excess/deficiency asymmetry across a band, no
dynamical contrast, no inelastic background, and no detector response. The
result says so, because a simulated pattern that looked photographic would
invite comparisons it cannot support.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ExampleScenario,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.app.services.traces import clipped_runs

__all__: tuple[str, ...] = ()

_CITATION_SCHWARTZ = (
    "Schwartz, Kumar, Adams & Field (eds.), Electron Backscatter Diffraction in Materials "
    "Science, 2nd ed., Springer 2009, Chs. 1-5."
)
_CITATION_BRITTON = (
    "Britton, Jiang, Guo, Vilalta-Clemente, Wallis, Hansen, Winkelmann & Wilkinson, "
    "'Tutorial: Crystal orientations and EBSD - Or which way is up?', Materials "
    "Characterization 117 (2016) 113, doi:10.1016/j.matchar.2016.04.008."
)
_CITATION_RANDLE = (
    "Randle & Engler, Introduction to Texture Analysis, 2nd ed., CRC Press 2010, Ch. 6."
)

#: Points sampled along each trace before clipping. A band edge is a conic and
#: the centre line curves in detector pixels on a tilted screen, so the polyline
#: must be fine enough that neither shows its vertices.
_TRACE_SAMPLES = 721

_BAND_COLUMNS = (
    Column("plane", "Plane"),
    Column("d_angstrom", "d", units="Å", numeric=True, digits=4),
    Column("bragg_angle_deg", "θ_B", units="deg", numeric=True, digits=4),
    Column("width_deg", "Band width 2θ_B", units="deg", numeric=True, digits=4),
    Column(
        "width_px",
        "Apparent width",
        units="px",
        numeric=True,
        digits=1,
        help_text=(
            "Measured across the pattern centre. A band away from the pattern centre appears "
            "wider than this, because the gnomonic projection stretches with distance from the "
            "origin."
        ),
    ),
    Column("intensity", "Relative intensity", numeric=True, digits=4),
)

_AXIS_COLUMNS = (
    Column("axis", "Zone axis"),
    Column("bands", "Bands meeting", numeric=True),
    Column("x_px", "u", units="px", numeric=True, digits=1),
    Column("y_px", "v", units="px", numeric=True, digits=1),
    Column("on_detector", "On the screen"),
)


@REGISTRY.operation(
    "ebsd.simulate_kikuchi_pattern",
    title="Simulate an EBSD Kikuchi pattern",
    summary="The pattern a tilted specimen of this phase would throw onto the camera.",
    help_text=(
        "The forward problem behind every indexed point in a scan: given a phase, an "
        "orientation and the camera geometry, which bands land where.\n\n"
        "**Configured the way the microscope is.** Stage tilt (70° on almost every EBSD "
        "stage), the camera's elevation and azimuth, the pattern centre and the camera "
        "distance as the fraction of the screen width that calibration reports. The "
        "accelerating voltage is 20 kV rather than the TEM's 200, and that single number is "
        "why an EBSD pattern looks as it does: the wavelength is about three times longer, so "
        "every band is about three times wider and a screen tens of degrees across shows a "
        "network rather than a few fine lines.\n\n"
        "**The frame is stated, not implied.** The beam is the laboratory z axis, the stage "
        "tilt axis is x, and the camera sits on -y. An untilted specimen faces the beam; "
        "tilting the stage tips the specimen normal towards the camera. That normal then "
        "projects at a gnomonic radius of tan(90° - tilt + elevation) — about 0.36 at 70° — "
        "which is the arithmetic to check a suspicious geometry against.\n\n"
        "**What it is for.** To see what an orientation *looks like* before or instead of "
        "measuring it; to check by eye that an indexing solution's bands fall where the "
        "measured ones do; and to teach the geometry, because every quantity in it can be "
        "changed here and the consequence seen immediately.\n\n"
        "**What it is not.** Kinematic and geometric: band positions and widths are exact, "
        "intensities are a |F|² proxy. There is no excess/deficiency asymmetry across a band, "
        "no dynamical contrast, no inelastic background and no detector response, so this is "
        "a map of where the bands are, not a photograph of one."
    ),
    parameters=(
        phase_parameter(help_text="The phase whose bands are simulated."),
        NumberParameter(
            name="phi1_deg",
            label="First Bunge angle",
            symbol="phi_1",
            help_text="First Bunge angle, about the specimen z axis. The orientation is "
            "crystal-to-specimen, as every EBSD file reports it.",
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="Phi_deg",
            label="Second Bunge angle",
            symbol="Phi",
            help_text="Second Bunge angle, about the new x axis.",
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="phi2_deg",
            label="Third Bunge angle",
            symbol="phi_2",
            help_text="Third Bunge angle, about the new z axis.",
            units="deg",
            default=0.0,
            group="Orientation (Bunge)",
            row="Orientation (Bunge)",
            field_width="short",
        ),
        NumberParameter(
            name="sample_tilt_deg",
            label="Sample tilt",
            help_text=(
                "Stage tilt about the tilt axis. Seventy degrees is the standard working "
                "point: it is a compromise between backscatter yield, which rises with tilt, "
                "and the foreshortening and depth-of-field loss that also rise with it."
            ),
            units="deg",
            default=70.0,
            minimum=0.0,
            maximum=89.0,
            group="Camera",
        ),
        NumberParameter(
            name="detector_elevation_deg",
            label="Camera elevation",
            help_text=(
                "How far the camera axis is raised above the plane perpendicular to the beam. "
                "Zero puts the screen normal horizontal; a few degrees is typical of a camera "
                "inserted above the nominal port."
            ),
            units="deg",
            default=0.0,
            minimum=-40.0,
            maximum=40.0,
            group="Camera",
            row="Camera",
        ),
        NumberParameter(
            name="detector_azimuth_deg",
            label="Camera azimuth",
            help_text="Rotation of the camera about the beam, from the nominal port.",
            units="deg",
            default=0.0,
            minimum=-180.0,
            maximum=180.0,
            group="Camera",
            row="Camera",
        ),
        NumberParameter(
            name="pattern_centre_x",
            label="Pattern centre, across the screen",
            help_text=(
                "Where the camera axis meets the screen, as a fraction of the screen width "
                "from its left edge."
            ),
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            group="Camera",
            symbol="pattern_centre_x",
            row="Pattern centre",
        ),
        NumberParameter(
            name="pattern_centre_y",
            label="Pattern centre, down the screen",
            help_text=(
                "The same, as a fraction of the screen height from its top edge. Vendors "
                "differ over this origin, so a value copied from a vendor file may need "
                "converting to the frame stated in the help above."
            ),
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            group="Camera",
            symbol="pattern_centre_y",
            row="Pattern centre",
        ),
        NumberParameter(
            name="detector_distance",
            label="Camera distance",
            help_text=(
                "Specimen-to-screen distance, as a fraction of the screen width. This is the "
                "quantity a pattern-centre calibration reports; 0.5 to 0.8 covers the usual "
                "working range, and a larger value narrows the field of view."
            ),
            default=0.65,
            minimum=0.1,
            maximum=3.0,
            group="Camera",
            symbol="detector_distance",
            row="Pattern centre",
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text=(
                "Sets the electron wavelength, and with it every band width: the width is "
                "2·arcsin(λ/2d), so halving the voltage widens every band by about 40 per cent."
            ),
            units="kV",
            default=20.0,
            minimum=5.0,
            maximum=40.0,
            group="Camera",
        ),
        IntegerParameter(
            name="detector_width_px",
            label="Screen width",
            help_text="Camera width in pixels, before any binning.",
            units="px",
            default=640,
            minimum=64,
            maximum=2048,
            group="Camera",
            row="Screen size",
        ),
        IntegerParameter(
            name="detector_height_px",
            label="Screen height",
            help_text="Camera height in pixels, before any binning.",
            units="px",
            default=480,
            minimum=64,
            maximum=2048,
            group="Camera",
            row="Screen size",
        ),
        NumberParameter(
            name="pixel_size_um",
            label="Pixel size",
            help_text=(
                "Screen pixel pitch. With the width in pixels it fixes the physical screen "
                "size, and so the field of view in degrees."
            ),
            units="µm",
            default=50.0,
            minimum=1.0,
            maximum=500.0,
            group="Camera",
            advanced=True,
        ),
        IntegerParameter(
            name="max_bands",
            label="Bands drawn",
            help_text="Keep the strongest this many bands. A real pattern shows tens.",
            default=24,
            minimum=1,
            maximum=120,
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest |h|, |k| or |l| enumerated. Raising it admits narrower, weaker "
            "high-order bands at cubic cost.",
            default=4,
            minimum=1,
            maximum=8,
            advanced=True,
        ),
        IntegerParameter(
            name="zone_axis_max_index",
            label="Zone-axis limit",
            help_text="Largest |u|, |v| or |w| enumerated when locating the zone axes the "
            "bands meet at.",
            default=3,
            minimum=1,
            maximum=6,
            advanced=True,
        ),
    ),
    returns=(
        "One row per band, with the zone axes as a second table; the band and axis geometry "
        "in screen pixels under `data`."
    ),
    panel="ebsd_kikuchi",
    citations=(_CITATION_SCHWARTZ, _CITATION_BRITTON, _CITATION_RANDLE),
    tags=("EBSD", "Kikuchi", "simulation", "pattern centre", "detector", "teaching"),
)
def _simulate_kikuchi_pattern(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.frame_catalog import SPECIMEN_FRAME
    from pytex.core.notation import format_direction_indices
    from pytex.core.orientation import Orientation
    from pytex.diffraction.kikuchi import simulate_kikuchi_pattern
    from pytex.diffraction.models import DiffractionGeometry

    spec, phase = phase_from_request(request["phase"])
    width_px = int(request["detector_width_px"])
    height_px = int(request["detector_height_px"])
    pixel_size = float(request["pixel_size_um"])
    tilt = float(request["sample_tilt_deg"])
    elevation = float(request["detector_elevation_deg"])
    distance = float(request["detector_distance"])
    energy = float(request["beam_energy_kev"])

    geometry = DiffractionGeometry.for_ebsd(
        beam_energy_kev=energy,
        sample_tilt_deg=tilt,
        detector_elevation_deg=elevation,
        detector_azimuth_deg=float(request["detector_azimuth_deg"]),
        pattern_center=(
            float(request["pattern_centre_x"]),
            float(request["pattern_centre_y"]),
            distance,
        ),
        detector_shape=(height_px, width_px),
        detector_pixel_size_um=(pixel_size, pixel_size),
    )
    orientation = Orientation.from_euler(
        float(request["phi1_deg"]),
        float(request["Phi_deg"]),
        float(request["phi2_deg"]),
        degrees=True,
        specimen_frame=SPECIMEN_FRAME,
        phase=phase,
    )
    try:
        pattern = simulate_kikuchi_pattern(
            geometry,
            phase,
            orientation,
            max_index=int(request["max_index"]),
            max_bands=int(request["max_bands"]),
            zone_axis_max_index=int(request["zone_axis_max_index"]),
        )
    except ValueError as error:
        raise InvalidInputError(
            f"No pattern could be simulated: {error}",
            field="max_index",
            hint="Raise the index limit, or check that the phase carries an atomic basis.",
        ) from error

    projection = pattern.projection
    to_px = projection.to_detector_px

    bands: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for band in pattern.bands:
        indices = tuple(int(value) for value in band.plane.indices)
        centre = clipped_runs(
            np.asarray(to_px(band.center_trace(projection, samples=_TRACE_SAMPLES)), dtype=float),
            width_px,
            height_px,
        )
        edges = [
            clipped_runs(np.asarray(to_px(edge), dtype=float), width_px, height_px)
            for edge in band.edge_traces(projection, samples=_TRACE_SAMPLES)
        ]
        if not centre and not any(edges):
            # The band's trace misses the screen entirely. Reporting it as a row
            # with nothing drawn would be a band the reader cannot find.
            continue
        gnomonic_width = band.width_at_pattern_center(projection)
        width_in_px = (
            float(gnomonic_width * geometry.camera_length_mm / (pixel_size / 1000.0))
            if math.isfinite(gnomonic_width)
            else float("inf")
        )
        label = plane_label(indices, spec=spec)
        bands.append(
            {
                "hkl": list(indices),
                "label": label,
                "d_angstrom": float(band.d_spacing_angstrom),
                "bragg_angle_deg": float(np.degrees(band.bragg_angle_rad)),
                "width_deg": float(np.degrees(band.angular_width_rad)),
                "width_px": width_in_px,
                "intensity": float(band.intensity),
                "centre": centre,
                "edges": edges,
            }
        )
        rows.append(
            {
                "plane": label,
                "d_angstrom": float(band.d_spacing_angstrom),
                "bragg_angle_deg": float(np.degrees(band.bragg_angle_rad)),
                "width_deg": float(np.degrees(band.angular_width_rad)),
                "width_px": width_in_px,
                "intensity": float(band.intensity),
            }
        )

    axes: list[dict[str, Any]] = []
    axis_rows: list[dict[str, Any]] = []
    for axis in pattern.zone_axes:
        point = np.asarray(to_px(np.asarray(axis.coordinates, dtype=float).reshape(1, 2)))[0]
        indices = tuple(int(value) for value in axis.indices)
        label = format_direction_indices(indices, style="plain")
        axes.append(
            {
                "uvw": list(indices),
                "label": label,
                "x": float(point[0]),
                "y": float(point[1]),
                "order": int(axis.band_count),
                "on_detector": bool(axis.on_detector),
            }
        )
        axis_rows.append(
            {
                "axis": label,
                "bands": int(axis.band_count),
                "x_px": float(point[0]),
                "y_px": float(point[1]),
                "on_detector": "yes" if axis.on_detector else "no",
            }
        )

    if not bands:
        raise InvalidInputError(
            "Every simulated band misses the screen at this geometry.",
            field="detector_distance",
            hint=(
                "A very large camera distance narrows the field of view. Reduce z*, or move "
                "the pattern centre back towards the middle of the screen."
            ),
        )

    on_screen = sum(1 for axis in axes if axis["on_detector"])
    normal_radius = math.tan(math.radians(90.0 - tilt + elevation))
    field_deg = 2.0 * math.degrees(
        math.atan(0.5 * width_px * (pixel_size / 1000.0) / geometry.camera_length_mm)
    )
    result = AppResult(
        title=f"Simulated EBSD pattern of {spec.name}",
        summary=(
            f"{len(bands)} band(s) on a {width_px}x{height_px} screen at {energy:g} kV, with the "
            f"stage tilted {tilt:g}° and the camera {distance:g} screen-widths away — a field of "
            f"view about {field_deg:.0f}° across, in which {on_screen} zone axis/axes fall. "
            f"The specimen normal projects at gnomonic radius {normal_radius:.3f} from the "
            "pattern centre, which is tan(90° - tilt + elevation) and fixes the geometry. Band "
            "positions and widths are exact; the intensities are kinematic."
        ),
        table=ResultTable(
            columns=_BAND_COLUMNS,
            rows=tuple(rows),
            caption=(
                f"Bands drawn for {spec.name}, strongest first. Width is 2θ_B, so the widest "
                "band has the smallest spacing — the opposite ordering to the intensities."
            ),
        ),
        data={
            "width_px": width_px,
            "height_px": height_px,
            "pattern_centre_px": [
                float(geometry.pattern_center_px[0]),
                float(geometry.pattern_center_px[1]),
            ],
            "camera_length_mm": float(geometry.camera_length_mm),
            "pixel_size_mm": pixel_size / 1000.0,
            "wavelength_angstrom": float(pattern.wavelength_angstrom),
            "field_of_view_deg": field_deg,
            "specimen_normal_gnomonic_radius": normal_radius,
            "bands": bands,
            "zone_axes": axes,
            # The axes travel as their own rows and columns rather than as a
            # second result table, because a result carries one table: the bands
            # are what the operation is about, and the panel draws the axes on
            # the pattern where they belong.
            "zone_axis_rows": axis_rows,
            "columns": [column.to_json() for column in _BAND_COLUMNS],
            "zone_axis_columns": [column.to_json() for column in _AXIS_COLUMNS],
        },
        inputs={
            "phase": spec.name,
            "euler_deg": [
                float(request["phi1_deg"]),
                float(request["Phi_deg"]),
                float(request["phi2_deg"]),
            ],
            "sample_tilt_deg": tilt,
            "detector_elevation_deg": elevation,
            "detector_distance": distance,
            "beam_energy_kev": energy,
        },
        notes=(
            "Kinematic and geometric. Band positions and widths are exact for this geometry; "
            "the intensities are a |F|² proxy, with no excess/deficiency asymmetry, no "
            "dynamical contrast and no background.",
            "Band centre lines are exactly straight in gnomonic coordinates and curve in screen "
            "pixels once the camera is tilted, which is why they are drawn as sampled polylines "
            "rather than as straight segments.",
            "The apparent width in the table is measured across the pattern centre. A band far "
            "from the pattern centre appears wider, because the gnomonic projection stretches "
            "with distance from the origin.",
        ),
        citations=(_CITATION_SCHWARTZ, _CITATION_BRITTON, _CITATION_RANDLE),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="ebsd_kikuchi.example.nickel_cube",
            title="Nickel in the cube orientation, 70° tilt",
            panel="ebsd_kikuchi",
            summary="A cube-oriented fcc grain on a standard stage, seen by a standard camera.",
            teaches=(
                "What the standard geometry actually shows: with the specimen normal 70° from "
                "the beam, the [001] pole is *not* at the pattern centre — the crystal "
                "direction there is about twenty degrees away from it, which is why an EBSD "
                "pattern of a cube-oriented grain does not look like a four-fold zone-axis "
                "pattern."
            ),
            operation="ebsd.simulate_kikuchi_pattern",
            request={
                "phase": {"builtin": "ni_fcc"},
                "phi1_deg": 0.0,
                "Phi_deg": 0.0,
                "phi2_deg": 0.0,
                "sample_tilt_deg": 70.0,
                "detector_distance": 0.65,
                "beam_energy_kev": 20.0,
                "max_bands": 24,
            },
        ),
        ExampleScenario(
            id="ebsd_kikuchi.example.voltage",
            title="The same grain at 10 kV",
            panel="ebsd_kikuchi",
            summary="Half the accelerating voltage, on the same geometry.",
            teaches=(
                "That the accelerating voltage is a band-width control and nothing else "
                "geometric: the wavelength rises by about 40 per cent, every band widens by "
                "the same factor, and the bands stay exactly where they were. It is the one "
                "parameter that changes how a pattern *looks* without changing what it means."
            ),
            operation="ebsd.simulate_kikuchi_pattern",
            request={
                "phase": {"builtin": "ni_fcc"},
                "sample_tilt_deg": 70.0,
                "detector_distance": 0.65,
                "beam_energy_kev": 10.0,
                "max_bands": 24,
            },
        ),
        ExampleScenario(
            id="ebsd_kikuchi.example.hexagonal",
            title="Zirconium with the basal pole near the beam",
            panel="ebsd_kikuchi",
            summary="An hcp phase at Bunge (0°, 30°, 0°) on the same camera.",
            teaches=(
                "How a hexagonal pattern differs from a cubic one at a glance: the wide, "
                "strong prismatic and pyramidal bands make a six-fold hub where the cubic "
                "pattern makes four-fold and three-fold ones, and the c/a ratio shows up "
                "directly as the difference in width between the basal and prismatic bands."
            ),
            operation="ebsd.simulate_kikuchi_pattern",
            request={
                "phase": {"builtin": "zr_hcp"},
                "phi1_deg": 0.0,
                "Phi_deg": 30.0,
                "phi2_deg": 0.0,
                "sample_tilt_deg": 70.0,
                "detector_distance": 0.6,
                "beam_energy_kev": 20.0,
                "max_bands": 28,
            },
        ),
    )
)
