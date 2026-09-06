# ruff: noqa: RUF001, RUF002
"""Convergent-beam electron diffraction for the shared web and desktop workbench.

The application layer implements no physics. It validates a human-scale request,
calls :mod:`pytex.diffraction.cbed`, and turns the returned pattern into the
common ``AppResult`` contract plus one thing the parallel-beam panels never
needed: a **raster image**.

Why a raster, and not spots
---------------------------
A SAED pattern is a set of points, and a point is a circle in an SVG. A CBED
pattern is not: each reflection is a *disc* whose interior carries the rocking
curve, and the fringes inside that disc are the measurement. Drawing the discs
as outlines would throw away the entire content of the technique. So the discs
are rasterised here, into one intensity image with a stated extent in
millimetres, and the panel draws the image with the disc outlines and labels
over it — vector where the geometry lives, raster where the intensity does.

The image travels as base64 8-bit greyscale, the same way
:mod:`pytex.app.services.crystal` sends a PNG, because JSON has no bytes.

What the numbers mean, and do not
---------------------------------
Every limitation is the simulator's own and is reported with the result rather
than hidden: the two-beam method computes each disc independently and therefore
displays a symmetry that belongs to the method and not to the crystal, which is
why symmetry determination refuses to run on it. See
:class:`pytex.diffraction.cbed.CBEDPattern` for the full statement.
"""

from __future__ import annotations

import base64
from typing import Any, cast

import numpy as np

from pytex.app.errors import InvalidInputError, UnsupportedRequestError
from pytex.app.logbook import APP_LOG, ProgressReporter
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label

__all__: tuple[str, ...] = ()

_CITATION_WILLIAMS_CARTER = (
    "Williams & Carter, Transmission Electron Microscopy, 2nd ed., Part 2 (CBED), Springer 2009."
)
_CITATION_BUXTON = (
    "Buxton, Eades, Steeds & Rackham, Phil. Trans. R. Soc. A 281 (1976) 171, "
    "doi:10.1098/rsta.1976.0024."
)
_CITATION_KELLY = (
    "Kelly, Jostsons, Blake & Napier, Phys. Status Solidi A 31 (1975) 771, "
    "doi:10.1002/pssa.2210310251."
)
_CITATION_SPENCE_ZUO = "Spence & Zuo, Electron Microdiffraction, Plenum 1992, Chs. 3-4."

#: Side of the rasterised pattern, in pixels. Large enough that the fringes
#: inside a disc survive — at 512 a disc of typical size spans roughly 80
#: pixels, more than the sampling of the rocking curve behind it — and small
#: enough that the payload stays a quarter-megabyte before base64.
_IMAGE_SIZE = 512

_DISC_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column("d_angstrom", "d", units="Å", numeric=True, digits=4),
    Column("g_inv_angstrom", "|g|", units="Å⁻¹", numeric=True, digits=4),
    Column("x_mm", "Disc centre x", units="mm", numeric=True, digits=3),
    Column("y_mm", "Disc centre y", units="mm", numeric=True, digits=3),
    Column(
        "extinction_distance_angstrom",
        "ξ_g",
        units="Å",
        numeric=True,
        digits=1,
        help_text=(
            "Extinction distance. The thickness over which the beam is exchanged once with the "
            "transmitted beam; it sets the fringe spacing inside this disc."
        ),
    ),
    Column("structure_factor_amplitude", "|F|", units="Å", numeric=True, digits=4),
    Column(
        "mean_intensity",
        "Mean disc intensity",
        numeric=True,
        digits=4,
        help_text="Averaged over the illuminated area of the disc, relative to the incident beam.",
    ),
)

_HOLZ_COLUMNS = (
    Column("order", "Laue zone", numeric=True),
    Column("radius_mm", "Ring radius", units="mm", numeric=True, digits=3),
    Column("radius_inv_angstrom", "Ring radius", units="Å⁻¹", numeric=True, digits=4),
)

_FRINGE_COLUMNS = (
    Column("order", "Minimum n", numeric=True),
    Column(
        "excitation_error_inv_angstrom",
        "s_n",
        units="Å⁻¹",
        numeric=True,
        digits=5,
        help_text="Deviation parameter at the nth fringe minimum, as read off the disc.",
    ),
    Column("inverse_order_squared", "1/n²", numeric=True, digits=5),
    Column("s_over_n_squared", "(s_n/n)²", units="Å⁻²", numeric=True, digits=8),
)


def _config_from_request(request: dict[str, Any]) -> Any:
    """Build the instrument settings, translating a refusal into a user message.

    The scientific constructor validates hard — a zero convergence angle is a
    parallel beam, which is a different technique rather than a degenerate CBED
    — and its message already says why. Re-raising it as an
    :class:`InvalidInputError` puts it beside the control that caused it instead
    of turning a stated limit into a 500.
    """

    from pytex.diffraction.cbed import ConvergentBeamConfig

    method = str(request["method"])
    laue = (0,) if method == "two-beam" or not request.get("include_holz_beams") else (0, 1)
    try:
        return ConvergentBeamConfig(
            beam_energy_kev=float(request["beam_energy_kev"]),
            convergence_semi_angle_mrad=float(request["convergence_semi_angle_mrad"]),
            thickness_angstrom=float(request["thickness_nm"]) * 10.0,
            camera_constant_mm_angstrom=float(request["camera_constant_mm_angstrom"]),
            max_index=int(request["max_index"]),
            g_max_inv_angstrom=float(request["g_max_inv_angstrom"]),
            disc_samples=int(request["disc_samples"]),
            method=cast(Any, method),
            laue_zones=laue,
        )
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            hint="Adjust the instrument settings in the Beam and Detector groups.",
        ) from error


def _rasterise(pattern: Any, *, size: int = _IMAGE_SIZE) -> dict[str, Any]:
    """Draw the discs into one greyscale image with a stated extent.

    Purpose
    -------
    Turns the per-disc intensity maps into the single picture the panel shows.

    Convention
    ----------
    A point inside a disc corresponds to one incident direction in the
    convergence cone, and it is drawn at the offset from the disc centre that
    that direction subtends: a tilt of ``(θ_u, θ_v)`` lands at
    ``(θ_u, θ_v) / α`` of the disc radius, with ``u`` along detector *x* and
    ``v`` along detector *y* — the same zone basis
    :func:`pytex.diffraction.cbed.simulate_cbed_pattern` builds the pattern in.
    This is the standard CBED disc coordinate: position within a disc *is*
    incident-beam direction.

    Overlapping discs add. In the Kossel regime the overlap is real and the
    total is what a detector integrates; what a two-beam simulation cannot
    supply there is the *interference* between the overlapping beams, which is
    why the result reports the regime rather than leaving the reader to assume
    the overlap is modelled.

    Returns
    -------
    dict
        ``width``, ``height``, ``extent_mm``, ``encoding``, ``data`` (base64
        8-bit greyscale, row-major from the top), and ``peak_intensity``, which
        is the physical value the byte 255 stands for.
    """

    discs = pattern.discs
    radius = float(pattern.config.disc_radius_mm)
    centres = np.asarray([disc.centre_mm for disc in discs], dtype=np.float64)
    extent = float(np.max(np.abs(centres)) + radius) * 1.06 if len(centres) else radius * 2.0

    axis = np.linspace(-extent, extent, size)
    # Row 0 is the top of the image, so y descends down the rows: without this
    # the pattern is drawn mirrored, and a mirrored diffraction pattern is a
    # different crystal.
    grid_x, grid_y = np.meshgrid(axis, -axis, indexing="xy")
    canvas = np.zeros((size, size), dtype=np.float64)

    reporter = ProgressReporter(
        "cbed.raster",
        total=len(discs),
        source="cbed.pattern",
        label="Drawing the convergent-beam discs",
    )
    for disc in discs:
        local_x = (grid_x - float(disc.centre_mm[0])) / radius
        local_y = (grid_y - float(disc.centre_mm[1])) / radius
        inside = (local_x * local_x + local_y * local_y) <= 1.0
        if not np.any(inside):
            reporter.advance()
            continue
        samples = disc.intensity.shape[0]
        # `intensity[i, j]` is sampled at tilt (axis[i], axis[j]) with the u
        # index first, so u indexes rows of the array and v indexes columns.
        index_u = np.clip(
            np.rint((local_x[inside] + 1.0) * 0.5 * (samples - 1)).astype(np.int64), 0, samples - 1
        )
        index_v = np.clip(
            np.rint((local_y[inside] + 1.0) * 0.5 * (samples - 1)).astype(np.int64), 0, samples - 1
        )
        values = np.asarray(disc.intensity, dtype=np.float64)[index_u, index_v]
        canvas[inside] += np.nan_to_num(values, nan=0.0)
        reporter.advance()
    reporter.finish()

    peak = float(canvas.max())
    # A pattern of exactly zero intensity is possible — every reflection
    # forbidden — and dividing by its peak would fill the frame with NaN.
    scaled = np.zeros_like(canvas) if peak <= 0.0 else canvas / peak
    bytes_out = np.clip(np.rint(scaled * 255.0), 0, 255).astype(np.uint8)
    return {
        "width": size,
        "height": size,
        "extent_mm": extent,
        "encoding": "base64-gray8",
        "data": base64.b64encode(bytes_out.tobytes()).decode("ascii"),
        "peak_intensity": peak,
    }


def _disc_rows(pattern: Any, spec: Any) -> list[dict[str, Any]]:
    """One row per drawn disc, in the order the pattern holds them."""

    rows: list[dict[str, Any]] = []
    for disc in pattern.discs:
        indices = tuple(int(value) for value in np.asarray(disc.miller_indices, dtype=int))
        magnitude = float(np.linalg.norm(np.asarray(disc.g_detector_inv_angstrom, dtype=float)))
        intensity = np.asarray(disc.intensity, dtype=float)
        finite = intensity[np.isfinite(intensity)]
        rows.append(
            {
                "hkl_label": "(000) direct beam"
                if disc.is_transmitted
                else plane_label(indices, spec=spec),
                "indices": list(indices),
                "d_angstrom": float(1.0 / magnitude) if magnitude > 0.0 else None,
                "g_inv_angstrom": magnitude,
                "x_mm": float(disc.centre_mm[0]),
                "y_mm": float(disc.centre_mm[1]),
                "extinction_distance_angstrom": (
                    None
                    if not np.isfinite(disc.extinction_distance_angstrom)
                    else float(disc.extinction_distance_angstrom)
                ),
                "structure_factor_amplitude": float(abs(disc.structure_factor_angstrom)),
                "mean_intensity": float(finite.mean()) if finite.size else 0.0,
            }
        )
    return rows


@REGISTRY.operation(
    "cbed.pattern",
    title="CBED pattern",
    summary=(
        "Simulate a zone-axis convergent-beam pattern, with the rocking curve inside each disc."
    ),
    help_text=(
        "Simulates what a convergent-beam exposure of the chosen phase and zone axis would show. "
        "Each zeroth-Laue-zone reflection becomes a disc of angular radius equal to the "
        "convergence semi-angle, and the disc is filled with the diffracted intensity at the "
        "excitation error of every incident direction in the illumination cone. Position inside a "
        "disc is therefore incident-beam direction, which is what makes the fringes readable.\n\n"
        "**Convergence angle is the parameter that matters.** Small enough and the discs stay "
        "apart — the Kossel-Moellenstedt regime, where each disc is an independent rocking curve "
        "and the fringes give the foil thickness. Larger and the discs overlap into the Kossel "
        "regime, where the overlaps carry interference between beams.\n\n"
        "**Two-beam** computes each disc on its own: cheap, and exactly the model the thickness "
        "measurement inverts. It is also symmetric in the excitation error by construction, so a "
        "two-beam pattern displays a symmetry belonging to the method rather than to the crystal, "
        "and symmetry determination is refused on it. **Bloch wave** solves the coupled many-beam "
        "problem, which is the only method whose relative intensities and symmetry mean anything "
        "— and the only one that can determine a point group, including whether the crystal is "
        "centrosymmetric, which kinematic diffraction cannot decide at all. It costs cubically in "
        "the number of beams, so reduce the disc sampling before widening the cut-off.\n\n"
        "Not modelled: inelastic background, probe aberrations, specimen bend or wedge, and — in "
        "the two-beam path — absorption, so fringes do not decay with thickness as they do on a "
        "real plate."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "The phase whose cell, symmetry and atomic sites generate the pattern. A unit "
                "cell with atom positions is required: the extinction distance that sets the "
                "fringe spacing is computed from the structure factor."
            ),
            builtin="si_diamond",
        ),
        IndicesParameter(
            name="zone_axis",
            label="Zone axis [uvw]",
            help_text=(
                "The beam direction in crystal indices. Low-index zones give the many-beam "
                "conditions CBED is used at; [001], [011] and [111] are the conventional choices."
            ),
            default=[0, 0, 1],
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text="Sets the electron wavelength and hence the Ewald-sphere curvature.",
            units="kV",
            default=200.0,
            minimum=20.0,
            maximum=1500.0,
            group="Beam",
        ),
        NumberParameter(
            name="convergence_semi_angle_mrad",
            label="Convergence semi-angle α",
            help_text=(
                "Half-angle of the illumination cone. It is the disc radius, and therefore what "
                "decides whether the discs stay separated or overlap. The discs separate when the "
                "disc diameter 2α/λ falls below the closest reciprocal-lattice spacing of the "
                "zone, so the usable angle depends on the material and the zone, not on the "
                "instrument alone."
            ),
            units="mrad",
            default=3.0,
            minimum=0.1,
            maximum=60.0,
            group="Beam",
        ),
        NumberParameter(
            name="thickness_nm",
            label="Foil thickness",
            help_text=(
                "Specimen thickness along the beam. It sets how many fringes cross each disc: "
                "roughly one more minimum per extinction distance of thickness."
            ),
            units="nm",
            default=100.0,
            minimum=1.0,
            maximum=2000.0,
            group="Beam",
        ),
        ChoiceParameter(
            name="method",
            label="Method",
            help_text=(
                "How disc intensities are computed. Start with two-beam; move to Bloch wave when "
                "relative intensities or symmetry matter."
            ),
            options=(
                (
                    "two-beam",
                    "Two-beam (fast)",
                    "Each disc computed independently by the closed-form rocking curve.",
                ),
                (
                    "bloch",
                    "Bloch wave (many-beam)",
                    "Coupled solution; the only method whose intensities and symmetry are "
                    "meaningful.",
                ),
            ),
            default="two-beam",
        ),
        NumberParameter(
            name="camera_constant_mm_angstrom",
            label="Camera constant Lλ",
            help_text="Detector scale, the same quantity the SAED indexing path is calibrated in.",
            units="mm·Å",
            default=180.0,
            minimum=1.0,
            maximum=5000.0,
            group="Detector",
            symbol="camera_constant",
        ),
        NumberParameter(
            name="g_max_inv_angstrom",
            label="Radial cut-off",
            help_text="Largest |g| drawn; the equivalent of the recorded detector extent.",
            units="Å⁻¹",
            default=1.2,
            minimum=0.1,
            maximum=6.0,
            group="Detector",
        ),
        IntegerParameter(
            name="max_index",
            label="Largest Miller index",
            help_text="Bound on the reflections enumerated for the zeroth Laue zone.",
            default=4,
            minimum=1,
            maximum=12,
            group="Detector",
            advanced=True,
        ),
        IntegerParameter(
            name="disc_samples",
            label="Samples across a disc",
            help_text=(
                "Resolution of the rocking curve inside each disc. It also bounds how precisely a "
                "fringe minimum — and therefore a fitted thickness — can be located. The Bloch "
                "path costs cubically in the beam count and linearly in this, so reduce it first."
            ),
            default=81,
            minimum=9,
            maximum=241,
            group="Detector",
            advanced=True,
        ),
        BooleanParameter(
            name="include_holz_beams",
            label="Admit HOLZ beams (Bloch only)",
            help_text=(
                "Include first-order Laue zone reflections in the coupled beam set. This is what "
                "produces HOLZ deficiency lines inside the discs, and what breaks the projection "
                "symmetry that would otherwise make every pattern look centrosymmetric — so it is "
                "required for an honest symmetry determination. It multiplies the beam count, and "
                "the cost with its cube."
            ),
            default=False,
            required=False,
            group="Detector",
            advanced=True,
        ),
        BooleanParameter(
            name="determine_point_group",
            label="Determine the point group",
            help_text=(
                "Read the pattern's symmetry and report which crystal point groups are consistent "
                "with it. Requires the Bloch method: a two-beam pattern is symmetric in the "
                "excitation error by construction, so its symmetry is the method's and not the "
                "crystal's."
            ),
            default=False,
            required=False,
            group="Detector",
            advanced=True,
        ),
    ),
    returns=(
        "A rasterised pattern under `data.image`, the disc list as a table, the overlap regime, "
        "the HOLZ ring radii, and — when asked and available — the point-group determination."
    ),
    panel="cbed",
    tags=("cbed", "convergent beam", "tem", "thickness", "symmetry", "holz", "point group"),
    citations=(_CITATION_WILLIAMS_CARTER, _CITATION_SPENCE_ZUO, _CITATION_BUXTON),
)
def _pattern(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.cbed import simulate_cbed_pattern

    spec, phase = phase_from_request(request["phase"])
    indices = tuple(int(value) for value in request["zone_axis"])
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise InvalidInputError(
            f"CBED needs the atom positions, because the extinction distance does: phase "
            f"'{phase.name}' carries no unit cell.",
            details={"field": "phase"},
            hint="Choose a built-in phase, or give the custom phase an atomic basis.",
        )
    config = _config_from_request(request)
    zone_axis = ZoneAxis(indices=np.asarray(indices, dtype=int), phase=phase)
    zone_text = f"[{' '.join(str(value) for value in indices)}]"

    method = str(request["method"])
    APP_LOG.info(
        f"Simulating {spec.name} down {zone_text} at {config.convergence_semi_angle_mrad:g} mrad "
        f"by the {'Bloch-wave' if method == 'bloch' else 'two-beam'} method.",
        source="cbed.pattern",
        detail={"method": method, "disc_samples": config.disc_samples},
    )
    if method == "bloch":
        APP_LOG.notice(
            "The Bloch-wave solution costs cubically in the beam count; this may take a while.",
            source="cbed.pattern",
        )

    pattern = simulate_cbed_pattern(phase, zone_axis, config=config)
    image = _rasterise(pattern)
    rows = _disc_rows(pattern, spec)

    separated = pattern.is_kossel_moellenstedt
    APP_LOG.notice(
        f"{len(pattern.discs)} discs, {'separated' if separated else 'overlapping'} "
        f"({pattern.regime}); nearest centres {pattern.nearest_disc_separation_mm:.3f} mm apart "
        f"and the disc diameter is {2 * config.disc_radius_mm:.3f} mm.",
        source="cbed.pattern",
        detail={"regime": pattern.regime, "discs": len(pattern.discs)},
    )

    determination = _determination(request, pattern)
    holz_rows = [
        {
            "order": int(order),
            "radius_mm": float(radius),
            "radius_inv_angstrom": float(radius) / config.camera_constant_mm_angstrom,
        }
        for order, radius in zip(pattern.holz_orders, pattern.holz_radii_mm, strict=True)
    ]

    notes = [
        (
            "Position inside a disc is incident-beam direction: a point at a fraction of the disc "
            "radius from its centre is the direction tilted by that fraction of the convergence "
            "semi-angle."
        ),
    ]
    if method == "two-beam":
        notes.append(
            "Two-beam: each disc is computed independently, no absorption enters, and every disc "
            "is symmetric in the excitation error by construction. The pattern therefore displays "
            "a symmetry that belongs to the method rather than to the crystal, which is why "
            "symmetry determination is refused on it."
        )
    else:
        notes.append(
            "Bloch wave: the discs are mutually consistent and the transmitted disc is a genuine "
            "bright-field intensity. Still unmodelled are inelastic background, probe "
            "aberrations, and specimen bend or wedge."
        )
    if not separated:
        notes.append(
            "The discs overlap, so this is the Kossel regime. Overlapping areas are drawn as the "
            "sum of the contributing discs; the interference between the overlapping beams, which "
            "is what an experiment actually records there, is not modelled. Reduce the "
            "convergence semi-angle to separate them."
        )

    result = AppResult(
        title=f"CBED of {spec.name} down {zone_text}",
        summary=(
            f"{len(pattern.discs)} discs of radius {config.disc_radius_mm:.3f} mm at "
            f"{config.beam_energy_kev:g} kV, {config.convergence_semi_angle_mrad:g} mrad "
            f"convergence and {config.thickness_angstrom / 10:g} nm thickness. The discs are "
            f"{'separated' if separated else 'overlapping'}, which is the "
            f"{pattern.regime.replace('-', '–')} regime"
            + (
                ", so each disc is an independent rocking curve and its fringes measure the foil "
                "thickness."
                if separated
                else ", so a disc is no longer an independent rocking curve and the fringe "
                "thickness method does not apply."
            )
            + (
                f" {len(holz_rows)} higher-order Laue zone ring(s) fall within the cut-off."
                if holz_rows
                else " No higher-order Laue zone ring falls within the cut-off."
            )
        ),
        table=ResultTable(
            columns=_DISC_COLUMNS,
            rows=tuple(rows),
            caption=f"Discs drawn for {spec.name} down {zone_text}, strongest structure factor "
            "first.",
        ),
        data={
            "image": image,
            "discs": rows,
            "columns": [column.to_json() for column in _DISC_COLUMNS],
            "disc_radius_mm": float(config.disc_radius_mm),
            "regime": pattern.regime,
            "separated": bool(separated),
            "nearest_disc_separation_mm": float(pattern.nearest_disc_separation_mm),
            "camera_constant_mm_angstrom": float(config.camera_constant_mm_angstrom),
            "zone_axis_label": zone_text,
            "phase_name": spec.name,
            "method": method,
            "holz_rings": holz_rows,
            "holz_columns": [column.to_json() for column in _HOLZ_COLUMNS],
            "point_group": determination,
        },
        inputs={
            "phase": spec.name,
            "zone_axis": list(indices),
            "beam_energy_kev": float(config.beam_energy_kev),
            "convergence_semi_angle_mrad": float(config.convergence_semi_angle_mrad),
            "thickness_nm": float(config.thickness_angstrom / 10.0),
            "camera_constant_mm_angstrom": float(config.camera_constant_mm_angstrom),
            "method": method,
            "disc_samples": int(config.disc_samples),
            "g_max_inv_angstrom": float(config.g_max_inv_angstrom),
            "max_index": int(config.max_index),
        },
        notes=tuple(notes),
        citations=(_CITATION_WILLIAMS_CARTER, _CITATION_SPENCE_ZUO, _CITATION_BUXTON),
    )
    return result.to_json()


def _determination(request: dict[str, Any], pattern: Any) -> dict[str, Any] | None:
    """Read the pattern's symmetry, or explain why that would be meaningless.

    A two-beam pattern is symmetric in the excitation error by construction, so
    determining a point group from one would report the method's symmetry as the
    crystal's. The scientific layer refuses; this turns that refusal into a
    message beside the control that asked for it.
    """

    if not request.get("determine_point_group"):
        return None
    if str(request["method"]) != "bloch":
        raise UnsupportedRequestError(
            "Symmetry determination needs the Bloch-wave method: a two-beam pattern is symmetric "
            "in the excitation error by construction, so its symmetry is the simulation's rather "
            "than the crystal's.",
            details={"field": "determine_point_group"},
            hint="Set Method to Bloch wave, and admit HOLZ beams for a determination that can "
            "distinguish a centre of symmetry.",
        )
    APP_LOG.info("Reading the pattern's symmetry.", source="cbed.pattern")
    try:
        determination = pattern.determine_point_group(
            require_holz=bool(request.get("include_holz_beams"))
        )
    except (ValueError, RuntimeError) as error:
        raise UnsupportedRequestError(
            str(error),
            details={"field": "determine_point_group"},
            hint="Admit HOLZ beams, or widen the radial cut-off so more discs enter the pattern.",
        ) from error
    payload = cast(dict[str, Any], determination.to_json_dict())
    APP_LOG.notice(
        f"Point groups consistent with the pattern: "
        f"{', '.join(payload.get('point_groups', [])) or 'none identified'}.",
        source="cbed.pattern",
    )
    return payload


@REGISTRY.operation(
    "cbed.thickness_from_fringes",
    title="Thickness from CBED fringes",
    summary="Foil thickness and extinction distance from the fringe minima of one disc.",
    help_text=(
        "The measurement CBED is used for most often. Inside a separated diffracted disc the "
        "intensity oscillates with the deviation parameter, and the nth minimum satisfies "
        "`t·s_eff,n = n` with `s_eff² = s² + ξ_g⁻²`. Rearranged, `(s_n/n)²` is linear in `1/n²` "
        "with intercept `t⁻²` and slope `−ξ_g⁻²`, so a single least-squares fit returns the "
        "thickness *and* the extinction distance — the extinction distance being the check, since "
        "it can be compared with the value the structure predicts.\n\n"
        "Enter the deviation parameters of the minima in order, outermost first or innermost "
        "first as you read them; what matters is that they are consecutive. The fit assigns the "
        "orders, searching over which order the first minimum carries, because misassigning n by "
        "one is the classic way this measurement goes wrong.\n\n"
        "This inverts the two-beam model, so it applies to a disc that is genuinely an "
        "independent rocking curve: separated discs, in the Kossel-Moellenstedt regime."
    ),
    parameters=(
        NumberParameter(
            name="s1",
            label="First minimum",
            help_text="Deviation parameter of the first fringe minimum used.",
            units="Å⁻¹",
            default=0.0071,
            group="Minima",
            symbol="s_1",
            row="Minima",
            field_width="short",
        ),
        NumberParameter(
            name="s2",
            label="Second minimum",
            help_text="Deviation parameter of the next minimum.",
            units="Å⁻¹",
            default=0.0128,
            group="Minima",
            symbol="s_2",
            row="Minima",
            field_width="short",
        ),
        NumberParameter(
            name="s3",
            label="Third minimum",
            help_text=(
                "Deviation parameter of the third minimum. Leave it blank to fit from two minima "
                "only — though a third is what makes the straight line worth drawing, since two "
                "points always lie on one."
            ),
            units="Å⁻¹",
            required=False,
            group="Minima",
            symbol="s_3",
            row="Minima",
            field_width="short",
        ),
        NumberParameter(
            name="s4",
            label="Fourth minimum",
            help_text="Deviation parameter of the fourth minimum, if resolved.",
            units="Å⁻¹",
            required=False,
            group="Minima",
            symbol="s_4",
            row="Minima",
            field_width="short",
        ),
        NumberParameter(
            name="s5",
            label="Fifth minimum",
            help_text="Deviation parameter of the fifth minimum, if resolved.",
            units="Å⁻¹",
            required=False,
            group="Minima",
            symbol="s_5",
            row="Minima",
            field_width="short",
        ),
        IntegerParameter(
            name="first_order",
            label="Order of the first minimum",
            help_text=(
                "Leave at 0 to let the fit search for it, which is the safe choice. Set it only "
                "when the order is known independently."
            ),
            default=0,
            minimum=0,
            maximum=6,
            required=False,
            advanced=True,
        ),
    ),
    returns="The fitted thickness and extinction distance, with the linearised fit as a table.",
    panel="cbed",
    tags=("cbed", "thickness", "extinction distance", "two-beam", "kelly"),
    citations=(_CITATION_KELLY, _CITATION_WILLIAMS_CARTER),
)
def _thickness_from_fringes(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.diffraction.cbed import thickness_from_fringe_minima

    minima = [
        float(request[key])
        for key in ("s1", "s2", "s3", "s4", "s5")
        if request.get(key) is not None
    ]
    if len(minima) < 2:
        raise InvalidInputError(
            f"A thickness fit needs at least two fringe minima; {len(minima)} were given.",
            details={"field": "s2"},
            hint="Read a second minimum off the disc and enter its deviation parameter.",
        )
    first_order = int(request.get("first_order") or 0) or None
    try:
        report = thickness_from_fringe_minima(minima, first_order=first_order)
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            hint="Check that the minima are consecutive and in order, and that none is zero.",
        ) from error

    payload = report.to_json_dict()
    orders = [int(report.first_order + offset) for offset in range(len(minima))]
    rows = [
        {
            "order": order,
            "excitation_error_inv_angstrom": value,
            "inverse_order_squared": 1.0 / (order * order),
            "s_over_n_squared": (value / order) ** 2,
        }
        for order, value in zip(orders, minima, strict=True)
    ]
    APP_LOG.notice(
        f"Fitted thickness {report.thickness_angstrom / 10.0:.1f} nm with extinction distance "
        f"{report.extinction_distance_angstrom:.1f} Å from {len(minima)} minima.",
        source="cbed.thickness_from_fringes",
        detail={"minima": len(minima), "first_order": report.first_order},
    )

    result = AppResult(
        title="Foil thickness from CBED fringe minima",
        summary=report.describe(),
        table=ResultTable(
            columns=_FRINGE_COLUMNS,
            rows=tuple(rows),
            caption="The linearised fit: (s_n/n)² against 1/n², whose intercept is t⁻² and whose "
            "slope is −ξ_g⁻².",
        ),
        data=dict(
            payload,
            fit_points=rows,
            columns=[column.to_json() for column in _FRINGE_COLUMNS],
            thickness_nm=float(report.thickness_angstrom) / 10.0,
        ),
        inputs={"minima_inv_angstrom": minima, "first_order": first_order},
        notes=(
            "This inverts the two-beam model, so it holds for a disc that is genuinely an "
            "independent rocking curve — a separated disc in the Kossel-Moellenstedt regime.",
            "The fitted extinction distance is the check on the result. Compare it with the value "
            "the structure predicts for that reflection; a large disagreement means the fringe "
            "orders were misassigned, which is this measurement's usual failure.",
        ),
        citations=(_CITATION_KELLY, _CITATION_WILLIAMS_CARTER),
    )
    return result.to_json()


@REGISTRY.operation(
    "cbed.holz_rings",
    title="HOLZ ring radii",
    summary="Higher-order Laue zone ring radii and the reciprocal-lattice repeat along the beam.",
    help_text=(
        "A zone-axis pattern is blind to one dimension: every zeroth-zone reflection is "
        "perpendicular to the beam, so nothing in it measures the lattice repeat *along* the "
        "beam. The higher-order Laue zone rings do. Ring `n` has projected radius "
        "`G_n ≈ √(2nH/λ)`, where `H = 1/|r_uvw|` is the reciprocal-lattice layer spacing along "
        "the zone axis, so a measured ring radius converts directly into `H`.\n\n"
        "This is how CBED measures local strain and composition: a change in `H` is a change in "
        "that lattice parameter, read from a ring whose radius can be measured to a fraction of a "
        "percent."
    ),
    parameters=(
        phase_parameter(
            help_text="The phase whose lattice sets the layer spacing along the zone axis.",
            builtin="si_diamond",
        ),
        IndicesParameter(
            name="zone_axis",
            label="Zone axis [uvw]",
            help_text="Beam direction. The layer spacing is measured along it.",
            default=[0, 0, 1],
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text="Sets the wavelength, and so the ring radius for a given layer spacing.",
            units="kV",
            default=200.0,
            minimum=20.0,
            maximum=1500.0,
        ),
        IntegerParameter(
            name="orders",
            label="Rings to report",
            help_text="How many Laue zones to enumerate outward from the zeroth.",
            default=2,
            minimum=1,
            maximum=6,
        ),
        NumberParameter(
            name="camera_constant_mm_angstrom",
            label="Camera constant Lλ",
            help_text="Used only to report the radii in millimetres beside the reciprocal units.",
            units="mm·Å",
            default=180.0,
            minimum=1.0,
            maximum=5000.0,
            advanced=True,
            symbol="camera_constant",
        ),
    ),
    returns="Ring radii in reciprocal angstrom and millimetres, with the layer spacing H.",
    panel="cbed",
    tags=("cbed", "holz", "laue zone", "lattice parameter", "strain"),
    citations=(_CITATION_WILLIAMS_CARTER, _CITATION_SPENCE_ZUO),
)
def _holz_rings(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.cbed import holz_ring_radii_inv_angstrom

    spec, phase = phase_from_request(request["phase"])
    indices = tuple(int(value) for value in request["zone_axis"])
    zone_axis = ZoneAxis(indices=np.asarray(indices, dtype=int), phase=phase)
    zone_text = f"[{' '.join(str(value) for value in indices)}]"
    camera = float(request["camera_constant_mm_angstrom"])

    # The scientific function returns the *orders* first and the radii second.
    orders_array, radii_array = holz_ring_radii_inv_angstrom(
        phase,
        zone_axis,
        beam_energy_kev=float(request["beam_energy_kev"]),
        orders=int(request["orders"]),
    )
    orders = np.asarray(orders_array, dtype=int)
    radii = np.asarray(radii_array, dtype=float)
    if radii.size == 0:
        raise InvalidInputError(
            f"No higher-order Laue zone of {spec.name} down {zone_text} carries an allowed "
            "reflection within the enumerated index range.",
            details={"field": "orders"},
            hint="Try a lower-index zone axis, or ask for more rings.",
        )

    # The layer spacing is read back out of the radii rather than re-derived
    # from the lattice: G_n = sqrt(2nH/lambda) inverts exactly, and taking it
    # from the same numbers the table shows makes the two impossible to
    # disagree — which a second, independent derivation of H would allow.
    from pytex.diffraction.cbed import ConvergentBeamConfig

    wavelength = ConvergentBeamConfig(
        beam_energy_kev=float(request["beam_energy_kev"])
    ).wavelength_angstrom
    layer_spacing = float(radii[0] ** 2 * wavelength / (2.0 * float(orders[0])))
    repeat = 1.0 / layer_spacing

    rows = [
        {
            "order": int(order),
            "radius_inv_angstrom": float(radius),
            "radius_mm": float(radius) * camera,
        }
        for order, radius in zip(orders, radii, strict=True)
    ]
    APP_LOG.notice(
        f"Layer spacing along {zone_text} is {layer_spacing:.5f} Å⁻¹ "
        f"(repeat {repeat:.4f} Å); first ring at {rows[0]['radius_inv_angstrom']:.4f} Å⁻¹.",
        source="cbed.holz_rings",
    )

    result = AppResult(
        title=f"HOLZ rings of {spec.name} down {zone_text}",
        summary=(
            f"The reciprocal-lattice layer spacing along {zone_text} is {layer_spacing:.5f} Å⁻¹, "
            f"which is the reciprocal of the {repeat:.4f} Å real-space repeat along the beam. At "
            f"{float(request['beam_energy_kev']):g} kV the first enumerated ring falls at "
            f"{rows[0]['radius_inv_angstrom']:.4f} Å⁻¹, which is {rows[0]['radius_mm']:.2f} mm at "
            f"a camera constant of {camera:g} mm·Å."
        ),
        table=ResultTable(
            columns=_HOLZ_COLUMNS,
            rows=tuple(rows),
            caption=f"Projected HOLZ ring radii for {spec.name} down {zone_text}.",
        ),
        data={
            "rings": rows,
            "columns": [column.to_json() for column in _HOLZ_COLUMNS],
            "layer_spacing_inv_angstrom": layer_spacing,
            "real_space_repeat_angstrom": repeat,
            "zone_axis_label": zone_text,
            "phase_name": spec.name,
        },
        inputs={
            "phase": spec.name,
            "zone_axis": list(indices),
            "beam_energy_kev": float(request["beam_energy_kev"]),
            "orders": int(request["orders"]),
            "camera_constant_mm_angstrom": camera,
        },
        notes=(
            "Radii are the kinematic projection √(2nH/λ). Dynamical shifts of the HOLZ lines "
            "themselves are not included, and it is those line positions, rather than the ring "
            "radius, that a high-accuracy lattice-parameter measurement uses.",
            "Only Laue zones carrying an allowed reflection within the enumerated index range are "
            "listed; a zone absent from the table is not necessarily absent from the crystal.",
        ),
        citations=(_CITATION_WILLIAMS_CARTER, _CITATION_SPENCE_ZUO),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="cbed.example.silicon_001_separated",
            title="Silicon [001] with separated discs",
            panel="cbed",
            summary="Diamond-cubic silicon down [001] at a 3 mrad probe and 100 nm thickness.",
            teaches=(
                "The Kossel-Moellenstedt regime: the discs stay apart, so each one is an "
                "independent rocking curve and the fringes across it are what measure the foil "
                "thickness."
            ),
            operation="cbed.pattern",
            request={
                "phase": {"builtin": "si_diamond"},
                "zone_axis": [0, 0, 1],
                "beam_energy_kev": 200.0,
                "convergence_semi_angle_mrad": 3.0,
                "thickness_nm": 100.0,
                "method": "two-beam",
                "camera_constant_mm_angstrom": 180.0,
                "g_max_inv_angstrom": 1.2,
                "max_index": 4,
                "disc_samples": 81,
            },
        ),
        ExampleScenario(
            id="cbed.example.silicon_001_overlapping",
            title="Silicon [001] with overlapping discs",
            panel="cbed",
            summary="The same zone and thickness at a four-times wider convergence angle.",
            teaches=(
                "Why the convergence semi-angle is the CBED parameter: at 12 mrad the discs merge "
                "into the Kossel regime, a disc is no longer an independent rocking curve, and "
                "the fringe thickness method stops applying."
            ),
            operation="cbed.pattern",
            request={
                "phase": {"builtin": "si_diamond"},
                "zone_axis": [0, 0, 1],
                "beam_energy_kev": 200.0,
                "convergence_semi_angle_mrad": 12.0,
                "thickness_nm": 100.0,
                "method": "two-beam",
                "camera_constant_mm_angstrom": 180.0,
                "g_max_inv_angstrom": 1.2,
                "max_index": 4,
                "disc_samples": 81,
            },
        ),
        ExampleScenario(
            id="cbed.example.thickness_from_minima",
            title="Foil thickness from three fringe minima",
            panel="cbed",
            summary="A two-beam thickness fit from three consecutive minima of one disc.",
            teaches=(
                "That one straight line gives two unknowns: the intercept of (s_n/n)² against "
                "1/n² is the thickness and its slope is the extinction distance, so the "
                "measurement carries its own consistency check."
            ),
            operation="cbed.thickness_from_fringes",
            request={"s1": 0.0071, "s2": 0.0128, "s3": 0.0184, "first_order": 0},
        ),
        ExampleScenario(
            id="cbed.example.holz_silicon_001",
            title="Silicon [001] Laue-zone rings",
            panel="cbed",
            summary="Where the first higher-order Laue zone rings fall at 200 kV.",
            teaches=(
                "That the ring radius measures the lattice repeat along the beam — the one "
                "dimension a zone-axis spot pattern is blind to, because every zeroth-zone "
                "reflection is perpendicular to the zone axis."
            ),
            operation="cbed.holz_rings",
            request={
                "phase": {"builtin": "si_diamond"},
                "zone_axis": [0, 0, 1],
                "beam_energy_kev": 200.0,
                "orders": 2,
                "camera_constant_mm_angstrom": 180.0,
            },
        ),
    )
)
