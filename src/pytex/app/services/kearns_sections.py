"""Kearns parameters from three measured theta-2theta scans, in one go.

What it does
    Takes three symmetric theta-2theta scans, one from each principal section of a
    tube or a plate, and returns the three Kearns parameters - ``f_a``, ``f_r``
    and ``f_t`` for a tube (the effective fractions of basal poles along the
    axial, radial and transverse directions) - together with every step between
    the scans and the numbers:

    1. the diffractogram of each section, with the reflections the phase
       predicts, the peaks found in it, and which were used;
    2. each reflection's measured and random-powder intensity, their ratio (the
       basal-pole density at that reflection's tilt to ``[0001]``) and its tilt;
    3. the tilt profile interpolated from those densities, and Kearns' quadrature
       over it, node by node, whose contributions sum to ``f``;
    4. the three values, their sum and the values normalised to sum to one.

Which scan is which
    Kearns' ``f`` along a direction ``d`` is measured on a section whose
    *surface normal* is ``d``: the symmetric scan diffracts only from planes
    parallel to the surface. So the scan labelled *axial* is taken on a
    cross-section of the tube (its surface normal is the tube axis), *radial* on
    the tube's outer surface, flattened, and *transverse* on a longitudinal
    section tangent to the wall. For a plate the same three slots are RD, ND and
    TD. Getting this assignment wrong swaps two values and is the commonest
    mistake the method allows.

Why the closure check means something here
    The three values are measured on three different surfaces, independently, so
    their sum is not 1 by construction. It departs from 1 because the available
    reflections sample the tilt range unevenly and because each section carries
    its own background and absorption errors; Kearns (1965) found sums of 0.94 to
    1.06. The departure is shown, and the normalised values ``f_i / sum`` beside
    the raw ones, as Mani Krishna *et al.* (2011) report them.

See ``docs/site/theory/kearns_parameter_and_basal_pole_texture.md``.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.logbook import APP_LOG
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    ExampleScenario,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
)
from pytex.app.results import AppResult, Column, ResultMetric, ResultStage, ResultTable
from pytex.app.services.calculator import phase_parameter
from pytex.app.services.kearns import (
    _CITATION_KEARNS,
    _CITATION_MANI,
    _ISOTROPIC,
    _basal_pole,
    _reflection_label,
    _report_payload,
    _specimen_frame,
)
from pytex.app.services.kearns_figures import (
    section_scan_figure,
    tilt_profile_figure,
    triad_figure,
)

__all__: tuple[str, ...] = ()

#: The three section slots, in the order the result reports them: the key the
#: request uses, the subscript of f, the direction's name for a tube and for a
#: plate, and the specimen axis the direction is.
_SECTIONS: tuple[tuple[str, str, str, str, tuple[float, float, float]], ...] = (
    ("axial", "a", "axial", "RD", (1.0, 0.0, 0.0)),
    ("radial", "r", "radial", "ND", (0.0, 0.0, 1.0)),
    ("transverse", "t", "transverse", "TD", (0.0, 1.0, 0.0)),
)

_RADIATION_OPTIONS = (
    # Greek alpha is the published name of the characteristic line, not a confusable a.
    ("cu_ka_doublet", "Cu Kα1/Kα2", "Common laboratory copper doublet."),  # noqa: RUF001
    ("cu_ka", "Cu Kα (single averaged line)", "One copper line without splitting."),  # noqa: RUF001
    ("mo_ka_doublet", "Mo Kα1/Kα2", "Short-wavelength molybdenum doublet."),  # noqa: RUF001
    ("co_ka_doublet", "Co Kα1/Kα2", "Useful for reducing Fe fluorescence."),  # noqa: RUF001
)

_PEAK_GROUP = "Peaks and intensities"
_DEMO_GROUP = "Demonstration scans"

#: Points kept per diffractogram in the payload: enough to draw every peak.
_MAX_PLOT_POINTS = 4000

_STATUS_TEXT = {
    "used": "used",
    "not_detected": "not detected: intensity taken as 0",
    "overlapped": "overlaps another reflection: excluded",
    "weak_random": "too weak in a random powder: excluded",
    "no_random_peak": "not found in the random standard: excluded",
}


def _radiation(name: str) -> Any:
    from pytex.diffraction.xrd import RadiationSpec

    return {
        "cu_ka_doublet": RadiationSpec.cu_ka_doublet,
        "cu_ka": RadiationSpec.cu_ka,
        "mo_ka_doublet": RadiationSpec.mo_ka_doublet,
        "co_ka_doublet": RadiationSpec.co_ka,
    }[name]()


def _expected_reflections(phase: Any, radiation: Any, low: float, high: float) -> list[Any]:
    """Reflection families the phase diffracts into the window, strongest intensities attached."""

    from pytex.diffraction.xrd import generate_powder_reflections

    reflections = generate_powder_reflections(
        phase, radiation=radiation, two_theta_range_deg=(low, high), max_index=6
    )
    return sorted(reflections, key=lambda reflection: reflection.two_theta_deg)


def _read_pattern(item: Any, *, field: str, radiation: Any) -> Any:
    from pytex.app.services.xrd import PATTERN_FILE_SUFFIXES
    from pytex.app.uploads import uploaded_file
    from pytex.diffraction.xrd_measurement import read_powder_pattern

    with uploaded_file(item, field=field, suffixes=PATTERN_FILE_SUFFIXES) as (path, name):
        try:
            return read_powder_pattern(path, name=name, radiation=radiation)
        except Exception as error:  # every reader failure means the same thing to a user
            raise InvalidInputError(
                f"{name} could not be read as a theta-2theta scan: {error}",
                field=field,
                hint="Open a two-column .xy/.csv/.dat file or a PANalytical .xrdml line scan.",
            ) from error


def _demonstration_patterns(
    phase: Any, radiation: Any, expected: list[Any], request: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, float], str]:
    """Three scans a diffractometer would record from a known tube texture.

    The texture is the one pilgered zirconium tubing is known for: basal poles
    tilted about 30 degrees from the radial direction towards the transverse
    (hoop) direction, with a spread. Its Kearns parameters are computed exactly,
    as the mean of cos^2 over the basal poles of the orientations, so the route
    can be read against the truth.
    """

    from pytex.core.lattice import CrystalPlane, MillerIndex
    from pytex.core.orientation import OrientationSet, Rotation
    from pytex.diffraction.xrd import _pseudo_voigt_profile
    from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
    from pytex.texture.kearns import basal_tilt_angle_deg

    spread = math.radians(float(request["demonstration_spread_deg"]))
    generator = np.random.default_rng(int(request["demonstration_seed"]))
    centres = ((180.0, 30.0, 0.0), (0.0, 30.0, 0.0))
    # Scattered in rotation space - a random axis, a Gaussian angle - rather than
    # by adding noise to each Euler angle. Euler noise piles poles up at Phi = 0,
    # where the coordinates are singular, and that spike of density in a region of
    # vanishing solid angle is not a texture any specimen has.
    count = 4000
    blocks = []
    for centre in centres:
        base = Rotation.from_euler(*centre, convention="bunge", degrees=True).as_matrix()
        axes = generator.normal(size=(count, 3))
        axes /= np.linalg.norm(axes, axis=1, keepdims=True)
        angles = np.abs(generator.normal(0.0, spread, size=count))
        skew = np.zeros((count, 3, 3))
        skew[:, 0, 1], skew[:, 0, 2] = -axes[:, 2], axes[:, 1]
        skew[:, 1, 0], skew[:, 1, 2] = axes[:, 2], -axes[:, 0]
        skew[:, 2, 0], skew[:, 2, 1] = -axes[:, 1], axes[:, 0]
        rotation = (
            np.eye(3)[None, :, :]
            + np.sin(angles)[:, None, None] * skew
            + (1.0 - np.cos(angles))[:, None, None] * np.einsum("nij,njk->nik", skew, skew)
        )
        blocks.append(np.einsum("ij,njk->nik", base, rotation))
    orientations = OrientationSet.from_matrices(
        np.concatenate(blocks, axis=0),
        specimen_frame=_specimen_frame(),
        crystal_frame=phase.crystal_frame,
        symmetry=phase.symmetry,
        phase=phase,
    )
    mapped = orientations.map_crystal_directions(np.array([0.0, 0.0, 1.0]))
    poles = np.array(getattr(mapped, "values", mapped), dtype=float)
    poles = poles / np.linalg.norm(poles, axis=1, keepdims=True)

    low = float(request["two_theta_min_deg"])
    high = float(request["two_theta_max_deg"])
    grid = np.arange(low, high + 1e-9, 0.02)
    strongest = max(reflection.intensity for reflection in expected)
    truth: dict[str, float] = {}
    patterns: dict[str, Any] = {}
    for key, _sub, name, _plate, vector in _SECTIONS:
        axis = np.asarray(vector)
        cosine = np.abs(poles @ axis)
        truth[key] = float(np.mean(cosine * cosine))
        # Basal-pole density at a tilt from this section's normal, in multiples
        # of random: the share of poles in a band of +/-4 degrees about it over
        # the band's share of the hemisphere. A band rather than 1-degree bins,
        # so the density near zero tilt is not a ratio of two tiny numbers.
        pole_tilts = np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0)))
        profile = 150.0 + 60.0 * np.exp(-(grid - low) / 25.0)
        for reflection in expected:
            plane = CrystalPlane(
                miller=MillerIndex(np.asarray(reflection.miller_indices, dtype=int), phase=phase),
                phase=phase,
            )
            tilt = float(basal_tilt_angle_deg(plane))
            lower, upper = max(tilt - 4.0, 0.0), min(tilt + 4.0, 90.0)
            share = float(np.mean((pole_tilts >= lower) & (pole_tilts <= upper)))
            area = math.cos(math.radians(lower)) - math.cos(math.radians(upper))
            pole_density = share / area
            area = 60000.0 * (reflection.intensity / strongest) * pole_density
            fwhm = 0.08 + 0.001 * reflection.two_theta_deg
            profile = profile + area * _pseudo_voigt_profile(
                grid, reflection.two_theta_deg, fwhm, 0.5
            )
            if radiation.kalpha2_wavelength_angstrom is not None:
                theta = math.radians(reflection.two_theta_deg / 2.0)
                ratio = radiation.kalpha2_wavelength_angstrom / radiation.wavelength_angstrom
                partner = 2.0 * math.degrees(math.asin(min(1.0, ratio * math.sin(theta))))
                profile = profile + area * radiation.kalpha2_relative_intensity * (
                    _pseudo_voigt_profile(grid, partner, fwhm, 0.5)
                )
        counts_noisy = generator.poisson(np.clip(profile, 0.0, None)).astype(float)
        patterns[key] = {
            "file": f"demonstration {name} section",
            "pattern": MeasuredPowderPattern(
                name=f"demonstration {name} section",
                two_theta_deg=grid,
                intensity=counts_noisy,
                radiation=radiation,
                synthetic=True,
            ),
        }
    return patterns, truth, "basal poles 30° either side of radial towards transverse"


def _fit_peaks(pattern: Any, *, radiation: Any, request: dict[str, Any], field: str) -> Any:
    from pytex.diffraction.xrd_peaks import detect_and_fit_peaks

    low = max(float(request["two_theta_min_deg"]), float(pattern.two_theta_deg[0]))
    high = min(float(request["two_theta_max_deg"]), float(pattern.two_theta_deg[-1]))
    if not low < high:
        raise InvalidInputError(
            f"{pattern.name} does not cover the angular range {low:g}-{high:g}°.",
            field=field,
            hint="Widen the angular range, or check the scan is the one you meant to open.",
        )
    try:
        return detect_and_fit_peaks(
            pattern,
            radiation=radiation,
            prominence_sigma=float(request["prominence_sigma"]),
            two_theta_range_deg=(low, high),
        )
    except ValueError as error:
        raise InvalidInputError(
            f"No peaks were found in {pattern.name}: {error}",
            field=field,
            hint="Lower the detection threshold, or check the angular range covers the peaks.",
        ) from error


def _match(peaks: Any, expected: list[Any], tolerance: float) -> dict[int, int | None]:
    """Nearest fitted peak within tolerance for each expected reflection; overlaps marked -1."""

    positions = np.array([peak.two_theta_deg for peak in peaks], dtype=float)
    matched: dict[int, int | None] = {}
    claims: dict[int, list[int]] = {}
    for index, reflection in enumerate(expected):
        if positions.size == 0:
            matched[index] = None
            continue
        nearest = int(np.argmin(np.abs(positions - reflection.two_theta_deg)))
        if abs(positions[nearest] - reflection.two_theta_deg) <= tolerance:
            matched[index] = nearest
            claims.setdefault(nearest, []).append(index)
        else:
            matched[index] = None
    for claimants in claims.values():
        if len(claimants) > 1:
            for index in claimants:
                matched[index] = -1
    return matched


def _peak_value(peak: Any, measure: str) -> float:
    return float(peak.integrated_intensity if measure == "integrated" else peak.height)


def _section(
    key: str,
    label: str,
    direction: tuple[float, float, float],
    source: dict[str, Any],
    *,
    phase: Any,
    radiation: Any,
    expected: list[Any],
    random_values: dict[int, float] | None,
    request: dict[str, Any],
) -> dict[str, Any]:
    """Everything one section contributes: peaks, reflection table, profile, f."""

    from pytex.core.lattice import CrystalPlane, MillerIndex
    from pytex.diffraction.xrd_background import estimate_background
    from pytex.texture.kearns import (
        DiffractogramReflection,
        basal_tilt_angle_deg,
        basal_tilt_profile,
        kearns_from_diffractogram,
    )

    pattern = source["pattern"]
    peaks = _fit_peaks(pattern, radiation=radiation, request=request, field="scan_files")
    matched = _match(peaks, expected, float(request["match_tolerance_deg"]))
    measure = str(request["intensity_measure"])
    strongest = max(reflection.intensity for reflection in expected)
    weak_limit = float(request["min_random_percent"]) / 100.0 * strongest

    rows: list[dict[str, Any]] = []
    reflections: list[Any] = []
    for index, reflection in enumerate(expected):
        plane = CrystalPlane(
            miller=MillerIndex(np.asarray(reflection.miller_indices, dtype=int), phase=phase),
            phase=phase,
        )
        tilt = float(basal_tilt_angle_deg(plane))
        match = matched[index]
        random_intensity = (
            float(reflection.intensity) if random_values is None else random_values.get(index)
        )
        fitted = peaks[match] if match is not None and match >= 0 else None
        if reflection.intensity < weak_limit:
            status = "weak_random"
        elif match == -1:
            status = "overlapped"
        elif random_intensity is None or random_intensity <= 0.0:
            status = "no_random_peak"
        elif fitted is None:
            status = "not_detected"
        else:
            status = "used"
        measured = _peak_value(fitted, measure) if fitted is not None else 0.0
        used = status in {"used", "not_detected"}
        if used:
            reflections.append(
                DiffractogramReflection(
                    plane=plane,
                    intensity=max(measured, 0.0),
                    random_intensity=float(random_intensity),  # type: ignore[arg-type]
                )
            )
        rows.append(
            {
                "plane": _reflection_label(DiffractogramReflection(plane=plane, intensity=0.0)),
                "basal_tilt_deg": tilt,
                "two_theta_expected_deg": float(reflection.two_theta_deg),
                "two_theta_fitted_deg": None if fitted is None else float(fitted.two_theta_deg),
                "offset_deg": (
                    None
                    if fitted is None
                    else float(fitted.two_theta_deg - reflection.two_theta_deg)
                ),
                "fwhm_deg": None if fitted is None else float(fitted.fwhm_deg),
                "measured": float(measured),
                "random": None if random_intensity is None else float(random_intensity),
                "density": (
                    float(measured / random_intensity) if used and random_intensity else None
                ),
                "status": status,
                "status_text": _STATUS_TEXT[status],
                "used": used,
            }
        )

    distinct_tilts = {round(float(r.basal_tilt_deg), 3) for r in reflections}
    if len(reflections) < 2 or len(distinct_tilts) < 2:
        raise InvalidInputError(
            f"The {label} section has fewer than two usable reflections at different basal "
            "tilts, so no tilt profile can be built from it.",
            field="scan_files",
            hint=(
                "Check the angular range covers (0002) and (10-10) at least, loosen the match "
                "tolerance, or lower the detection threshold."
            ),
        )

    normalization = str(request["normalization"])
    bin_width = float(request["bin_width_deg"])
    report = kearns_from_diffractogram(
        reflections,
        specimen_frame=_specimen_frame(),
        normalization=normalization,
        bin_width_deg=bin_width,
        direction=direction,
        direction_label=label,
    )
    polar, profile, profile_notes = basal_tilt_profile(
        reflections, normalization=normalization, bin_width_deg=bin_width
    )
    polar = np.asarray(polar, dtype=float)
    profile = np.asarray(profile, dtype=float)
    weights = profile * np.sin(np.deg2rad(polar))
    total = float(weights.sum())
    fractions = weights / total if total > 0.0 else np.zeros_like(weights)
    cos2 = np.cos(np.deg2rad(polar)) ** 2
    quadrature = [
        {
            "polar_deg": float(polar[i]),
            "density": float(profile[i]),
            "sin": float(math.sin(math.radians(polar[i]))),
            "volume_fraction": float(fractions[i]),
            "cos2": float(cos2[i]),
            "contribution": float(fractions[i] * cos2[i]),
        }
        for i in range(polar.size)
    ]

    background = np.asarray(
        estimate_background(pattern, method="snip", half_window_deg=2.0).background, dtype=float
    )
    two_theta = np.asarray(pattern.two_theta_deg, dtype=float)
    intensity = np.asarray(pattern.intensity, dtype=float)
    stride = max(1, math.ceil(two_theta.size / _MAX_PLOT_POINTS))
    value = float(report.value(label))
    return {
        "key": key,
        "label": label,
        "file": source["file"],
        "direction": list(direction),
        "f": value,
        "report": report,
        "pattern": {
            "two_theta_deg": [float(v) for v in two_theta[::stride]],
            "intensity": [float(v) for v in intensity[::stride]],
            "background": [float(v) for v in background[::stride]],
        },
        "peaks": [
            {
                "two_theta_deg": float(peak.two_theta_deg),
                "height": float(peak.height),
                "integrated_intensity": float(peak.integrated_intensity),
                "fwhm_deg": float(peak.fwhm_deg),
            }
            for peak in peaks
        ],
        "reflections": rows,
        "profile": {
            "polar_deg": [float(v) for v in polar],
            "intensity": [float(v) for v in profile],
            "notes": list(profile_notes),
        },
        "quadrature": quadrature,
        "used_count": sum(1 for row in rows if row["used"]),
        "notes": list(report.notes),
    }


_REFLECTION_COLUMNS = (
    Column("plane", "Reflection", help_text="The reflecting plane, in Miller-Bravais form."),
    Column(
        "basal_tilt_deg",
        "Tilt to [0001]",
        units="°",
        numeric=True,
        digits=2,
        help_text="Angle between this plane's normal and the c axis, from the phase metric.",
    ),
    Column("two_theta_expected_deg", "2θ expected", units="°", numeric=True, digits=3),
    Column("two_theta_fitted_deg", "2θ fitted", units="°", numeric=True, digits=3),
    Column("measured", "Measured I", numeric=True, digits=1),
    Column(
        "random",
        "Random I",
        numeric=True,
        digits=2,
        help_text="The same reflection in a random powder: calculated or measured.",
    ),
    Column(
        "density",
        "I / I random",
        numeric=True,
        digits=4,
        help_text="The basal-pole density at this reflection's tilt, before normalisation.",
    ),
    Column("status_text", "Status"),
)

_QUADRATURE_COLUMNS = (
    Column("polar_deg", "Tilt φ", units="°", numeric=True, digits=1),
    Column("density", "Pole density I(φ)", numeric=True, digits=4),
    Column("sin", "sin φ", numeric=True, digits=4),
    Column(
        "volume_fraction",
        "Volume fraction",
        numeric=True,
        digits=5,
        help_text="I sin φ, normalised: the fraction of crystals at this tilt.",
    ),
    Column("cos2", "cos²φ", numeric=True, digits=5),
    Column(
        "contribution",
        "Contribution to f",
        numeric=True,
        digits=5,
        help_text="Volume fraction times cos²φ. These sum to f.",
    ),
)


@REGISTRY.operation(
    "kearns.from_three_sections",
    title="Kearns f from three measured scans",
    summary=(
        "Open a theta-2theta scan of each principal section and get f along all three "
        "directions in one go, with every peak, intensity and quadrature step shown."
    ),
    help_text=(
        "Kearns' 1965 route applied to three measured scans at once, one per principal section, "
        "so the whole triad comes out of one calculation that can be checked step by step.\n\n"
        "**Which scan is which.** f along a direction is measured on a section whose *surface "
        "normal* is that direction, because a symmetric scan only diffracts from planes parallel "
        "to the surface. For a tube: *axial* is a cross-section (its normal is the tube axis), "
        "*radial* is the outer surface, flattened, and *transverse* is a longitudinal section "
        "tangent to the wall. For a plate the three are RD, ND and TD.\n\n"
        "**What is done to each scan.** Peaks are detected against the instrument noise and "
        "fitted with pseudo-Voigt profiles (K-alpha2 modelled). Each reflection the phase predicts "
        "in the angular range is matched to the nearest fitted peak within the tolerance; its "
        "intensity - the fitted area by default - is divided by the same reflection's intensity "
        "in a random powder, which is the basal-pole density at that reflection's tilt to "
        "[0001]. Reflections that overlap another, or that are too weak in a random powder to "
        "measure, are excluded and say so; a reflection with a random intensity worth measuring "
        "but no peak is counted as zero, because a missing peak is texture too.\n\n"
        "**The random powder.** By default its intensities are calculated from the structure, "
        "with multiplicity and the Lorentz-polarisation factor. A measured random standard is "
        "better where absorption or the instrument matter; open one and choose it.\n\n"
        "**The integration.** The densities are interpolated onto tilt bins and integrated by "
        "Kearns' Eq. (5): f is the cos²-weighted mean of the volume fraction I(φ) sin φ. The "
        "table of every node is shown, and its contributions sum to f.\n\n"
        "**The closure check is real here.** The three values come from three independent "
        "measurements, so their sum is not 1 by construction; its departure measures the "
        "method's systematic error. The values normalised to sum to one are shown beside the "
        "raw ones.\n\n"
        "**With no scans open**, three demonstration scans of a known pilgered-tube texture are "
        "analysed, and the exact f of that texture is shown for comparison."
    ),
    parameters=(
        ObjectParameter(
            name="scan_files",
            label="Section scans",
            help_text=(
                'The three opened scans, as `{"axial": {"name": ..., "text": ...}, "radial": ..., '
                '"transverse": ...}`. Supplied by the three **Open section scans** controls.'
            ),
            required=False,
        ),
        ObjectParameter(
            name="random_file",
            label="Random standard scan",
            help_text=(
                "A scan of a random powder of the same phase, measured under the same conditions. "
                "Used when the random intensities come from a measured standard."
            ),
            required=False,
        ),
        phase_parameter(
            label="Phase",
            help_text="The hexagonal phase the scans were measured on.",
            builtin="zr_hcp",
        ),
        ChoiceParameter(
            name="geometry",
            label="Specimen",
            help_text=(
                "Names the three directions. A tube reports f along axial, radial and transverse; "
                "a plate along RD, ND and TD. The calculation is the same."
            ),
            options=(
                ("tube", "Tube (axial, radial, transverse)", "Pressure tubes and cladding."),
                ("plate", "Plate or sheet (RD, ND, TD)", "Rolled product."),
            ),
            default="tube",
        ),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text=(
                "The radiation the scans were measured with. It places the expected reflections "
                "and the K-alpha2 partners, and weights the calculated random intensities."
            ),
            options=_RADIATION_OPTIONS,
            default="cu_ka_doublet",
        ),
        ChoiceParameter(
            name="random_source",
            label="Random intensities",
            help_text=(
                "Where each reflection's random-powder intensity comes from. Calculated uses "
                "the structure; measured reads the random standard scan."
            ),
            options=(
                ("calculated", "Calculated from the structure", "Multiplicity, |F|² and LP."),
                ("measured", "Measured random standard", "The opened random standard scan."),
            ),
            default="calculated",
        ),
        NumberParameter(
            name="two_theta_min_deg",
            label="Start angle",
            help_text="Lower edge of the angular range searched for reflections.",
            units="° 2θ",
            default=28.0,
            minimum=5.0,
            maximum=170.0,
            group=_PEAK_GROUP,
            row="Angular range",
        ),
        NumberParameter(
            name="two_theta_max_deg",
            label="End angle",
            help_text="Upper edge of the angular range searched for reflections.",
            units="° 2θ",
            default=100.0,
            minimum=10.0,
            maximum=175.0,
            group=_PEAK_GROUP,
            row="Angular range",
        ),
        ChoiceParameter(
            name="intensity_measure",
            label="Peak intensity",
            help_text=(
                "Integrated intensity (the fitted area) is what Kearns used and what is "
                "insensitive to peak broadening; peak height is offered for comparison."
            ),
            options=(
                ("integrated", "Integrated (fitted area)", "Kearns' choice."),
                ("height", "Peak height", "Sensitive to broadening differences."),
            ),
            default="integrated",
            group=_PEAK_GROUP,
        ),
        NumberParameter(
            name="match_tolerance_deg",
            label="Match tolerance",
            help_text=(
                "How far a fitted peak may sit from a predicted reflection and still be taken as "
                "it. Wide enough for specimen displacement, narrow enough to keep neighbours "
                "apart."
            ),
            units="° 2θ",
            default=0.3,
            minimum=0.02,
            maximum=2.0,
            group=_PEAK_GROUP,
        ),
        NumberParameter(
            name="min_random_percent",
            label="Weakest random line",
            help_text=(
                "Reflections weaker than this percentage of the strongest one in a random powder "
                "are excluded: they cannot be measured reliably."
            ),
            units="%",
            default=1.0,
            minimum=0.0,
            maximum=50.0,
            group=_PEAK_GROUP,
        ),
        NumberParameter(
            name="prominence_sigma",
            label="Detection threshold",
            help_text="Peak detection threshold, in robust noise standard deviations.",
            default=4.0,
            minimum=1.0,
            maximum=50.0,
            group=_PEAK_GROUP,
        ),
        ChoiceParameter(
            name="normalization",
            label="Normalisation",
            help_text="How intensity ratios become pole densities. Both give the same f.",
            options=(
                ("random_standard", "Against the random powder", "Kearns' own choice."),
                ("harris", "Harris texture coefficients", "Mean coefficient of 1."),
            ),
            default="random_standard",
            group="Integration",
        ),
        NumberParameter(
            name="bin_width_deg",
            label="Tilt bin width",
            help_text=(
                "Width of the tilt bins the densities are interpolated onto. Kearns used 10°."
            ),
            units="°",
            default=10.0,
            minimum=2.0,
            maximum=30.0,
            group="Integration",
        ),
        NumberParameter(
            name="demonstration_spread_deg",
            label="Demonstration spread",
            help_text="Spread of the demonstration texture about its ideal basal tilts.",
            units="°",
            default=15.0,
            minimum=3.0,
            maximum=40.0,
            group=_DEMO_GROUP,
            group_collapsed=True,
        ),
        IntegerParameter(
            name="demonstration_seed",
            label="Demonstration seed",
            help_text="Seeds the demonstration texture and its counting noise.",
            default=5,
            minimum=0,
            maximum=1_000_000,
            group=_DEMO_GROUP,
            group_collapsed=True,
            field_width="short",
        ),
    ),
    returns=(
        "One row per direction with f and its normalised value; each section's diffractogram, "
        "peaks, reflection table, tilt profile and quadrature under `data.sections`."
    ),
    panel="kearns",
    citations=(_CITATION_KEARNS, _CITATION_MANI),
    tags=("kearns", "XRD", "theta-2theta", "three sections", "tube", "Fa", "Fr", "Ft", "triad"),
)
def _from_three_sections(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.texture.kearns import KearnsReport

    spec, phase = phase_from_request(request["phase"])
    pole = _basal_pole(phase)
    radiation = _radiation(str(request["radiation"]))
    low = float(request["two_theta_min_deg"])
    high = float(request["two_theta_max_deg"])
    if not low < high:
        raise InvalidInputError(
            "The start angle must be below the end angle.", field="two_theta_min_deg"
        )
    expected = _expected_reflections(phase, radiation, low, high)
    if len(expected) < 2:
        raise InvalidInputError(
            f"{spec.name} has fewer than two reflections between {low:g} and {high:g}°.",
            field="two_theta_max_deg",
        )
    geometry = str(request["geometry"])

    files = request.get("scan_files") or {}
    if not isinstance(files, dict):
        raise InvalidInputError(
            "The section scans arrived in an unexpected shape.", field="scan_files"
        )
    opened: dict[str, dict[str, Any]] = {
        key: dict(files[key]) for key, *_ in _SECTIONS if isinstance(files.get(key), dict)
    }
    truth: dict[str, float] | None = None
    if not opened:
        sources, truth, texture_name = _demonstration_patterns(phase, radiation, expected, request)
        demonstration = True
    else:
        missing = [key for key, *_ in _SECTIONS if key not in opened]
        if missing:
            raise InvalidInputError(
                "A triad needs all three sections; no scan was opened for "
                + ", ".join(missing)
                + ".",
                field="scan_files",
                hint="Open one scan per section. The single-scan route computes one direction.",
            )
        sources = {
            key: {
                "file": str(opened[key].get("name", key)),
                "pattern": _read_pattern(opened[key], field="scan_files", radiation=radiation),
            }
            for key in opened
        }
        texture_name = ""
        demonstration = False

    random_values: dict[int, float] | None = None
    random_note = "calculated from the structure (multiplicity, |F|² and Lorentz-polarisation)"
    if str(request["random_source"]) == "measured":
        random_item = request.get("random_file")
        if not random_item:
            raise InvalidInputError(
                "Measured random intensities were chosen, but no random standard scan is open.",
                field="random_file",
                hint="Open the random standard, or use intensities calculated from the structure.",
            )
        random_pattern = _read_pattern(random_item, field="random_file", radiation=radiation)
        random_peaks = _fit_peaks(
            random_pattern, radiation=radiation, request=request, field="random_file"
        )
        random_match = _match(random_peaks, expected, float(request["match_tolerance_deg"]))
        random_values = {
            index: _peak_value(random_peaks[m], str(request["intensity_measure"]))
            for index, m in random_match.items()
            if m is not None and m >= 0
        }
        random_note = f"measured on the random standard {random_item.get('name', '')}"

    sections = []
    for key, sub, tube_name, plate_name, vector in _SECTIONS:
        label = sub if geometry == "tube" else plate_name
        sections.append(
            _section(
                key,
                label,
                vector,
                sources[key],
                phase=phase,
                radiation=radiation,
                expected=expected,
                random_values=random_values,
                request=request,
            )
        )
        sections[-1]["name"] = tube_name if geometry == "tube" else plate_name

    values = np.array([section["f"] for section in sections], dtype=float)
    total = float(values.sum())
    triad = KearnsReport(
        values=values,
        directions=np.array([section["direction"] for section in sections], dtype=float),
        direction_labels=tuple(section["label"] for section in sections),
        method="diffractogram",
        pole=pole,
        specimen_frame=_specimen_frame(),
        orientation_tensor=None,
        diagnostics={
            "triad_departure": abs(total - 1.0),
            **{f"reflections_used_{s['key']}": float(s["used_count"]) for s in sections},
        },
        notes=(
            "Each value is Kearns' diffractogram route on its own section; the three sections "
            "were measured independently, so their sum is a genuine check.",
        ),
    )
    APP_LOG.info(
        "Kearns triad from three sections: "
        + ", ".join(f"f_{s['label']}={s['f']:.4f}" for s in sections),
        source="kearns.from_three_sections",
    )

    rows = [
        {
            "direction": f"f_{section['label']} ({section['name']})",
            "f": section["f"],
            "normalized": section["f"] / total if total > 0.0 else float("nan"),
            "vs_random": section["f"] / _ISOTROPIC,
            "reflections": section["used_count"],
            "file": section["file"],
            **({"truth": truth[section["key"]]} if truth is not None else {}),
        }
        for section in sections
    ]
    columns: tuple[Column, ...] = (
        Column("direction", "Direction"),
        Column("f", "f (measured)", numeric=True, digits=4),
        Column(
            "normalized",
            "f / sum",
            numeric=True,
            digits=4,
            help_text=(
                "The three values scaled to sum to one, as Mani Krishna et al. (2011) report."
            ),
        ),
        Column("vs_random", "f / (1/3)", numeric=True, digits=3),
        Column("reflections", "Reflections used", numeric=True),
        *(
            (
                Column(
                    "truth",
                    "Exact f of the model",
                    numeric=True,
                    digits=4,
                    help_text="Mean of cos² over the basal poles of the demonstration texture.",
                ),
            )
            if truth is not None
            else ()
        ),
        Column("file", "Scan"),
    )

    triad_text = ", ".join(f"f_{s['label']} = {s['f']:.4f}" for s in sections)
    departure = abs(total - 1.0)
    summary = (
        f"{triad_text} from three theta-2theta scans of {spec.name}"
        + (
            f" (demonstration scans of a tube with {texture_name})"
            if demonstration
            else f" ({', '.join(s['file'] for s in sections)})"
        )
        + f". They sum to {total:.4f}, {departure:.4f} from 1; because the three sections were "
        "measured independently that departure is the measurement's systematic error, and "
        "normalised to sum to one the values are "
        + ", ".join(f"{s['f'] / total:.4f}" for s in sections)
        + f". Random intensities were {random_note}; "
        + ", ".join(f"{s['used_count']} reflections" for s in sections)
        + " were used in the three sections."
        + (
            " The exact values of the model texture are "
            + ", ".join(f"f_{s['label']} = {truth[s['key']]:.4f}" for s in sections)
            + "."
            if truth is not None
            else ""
        )
    )

    stages: list[ResultStage] = []
    for number, section in enumerate(sections, start=1):
        excluded = [row for row in section["reflections"] if not row["used"]]
        stages.append(
            ResultStage(
                key=f"{section['key']}_reflections",
                title=(
                    f"{2 * number - 1}. {section['name'].capitalize()} section: "
                    "peaks and intensities"
                ),
                summary=(
                    f"{len(section['peaks'])} peaks fitted in {section['file']}; "
                    f"{section['used_count']} of {len(section['reflections'])} predicted "
                    f"reflections used, {len(excluded)} excluded."
                ),
                metrics=(
                    ResultMetric("Peaks fitted", len(section["peaks"])),
                    ResultMetric("Reflections used", section["used_count"]),
                ),
                table=ResultTable(columns=_REFLECTION_COLUMNS, rows=tuple(section["reflections"])),
                figures=(
                    section_scan_figure(
                        key=f"{section['key']}_scan",
                        name=section["name"],
                        pattern=section["pattern"],
                        reflections=section["reflections"],
                    ),
                ),
                explanation=(
                    "Measured over random intensity is the basal-pole density at the reflection's "
                    "tilt to [0001]. A reflection with a large 2θ offset has been matched to the "
                    "wrong peak or the specimen is displaced; check it on the diffractogram."
                ),
            )
        )
        stages.append(
            ResultStage(
                key=f"{section['key']}_quadrature",
                title=f"{2 * number}. {section['name'].capitalize()} section: Kearns' quadrature",
                summary=(
                    f"f_{section['label']} = {section['f']:.4f}: the contributions below sum to it."
                ),
                metrics=(ResultMetric(f"f_{section['label']}", section["f"]),),
                table=ResultTable(columns=_QUADRATURE_COLUMNS, rows=tuple(section["quadrature"])),
                figures=(
                    tilt_profile_figure(
                        key=f"{section['key']}_profile",
                        polar_deg=section["profile"]["polar_deg"],
                        density=section["profile"]["intensity"],
                        direction=section["label"],
                        f_value=float(section["f"]),
                        points=[
                            {
                                "tilt": row["basal_tilt_deg"],
                                "density": row["density"],
                                "label": row["plane"],
                            }
                            for row in section["reflections"]
                            if row["used"] and row["density"] is not None
                        ],
                        title=f"{section['name'].capitalize()} section: how f_{section['label']} "
                        "is built",
                    ),
                ),
                explanation=(
                    "f = sum of I(φ) sin φ cos²φ over sum of I(φ) sin φ. "
                    "The sin φ factor turns "
                    "pole density into volume fraction, so low tilts weigh little however intense."
                ),
            )
        )
    stages.append(
        ResultStage(
            key="triad",
            title="7. The triad and its closure",
            summary=f"Sum {total:.4f}; departure from 1 is {departure:.4f}.",
            metrics=(
                ResultMetric("Sum", total),
                ResultMetric("Departure from 1", departure),
            ),
            explanation=(
                "Kearns (1965) found sums between 0.94 and 1.06 from the unevenly spaced "
                "reflections available. A larger departure points at one section: compare their "
                "reflection tables and backgrounds."
            ),
            status="ok" if departure <= 0.1 else "warning",
            figures=(
                triad_figure(
                    key="kearns_triad",
                    labels=[section["label"] for section in sections],
                    values=[float(section["f"]) for section in sections],
                    closure_by_construction=False,
                    truth=(
                        None
                        if truth is None
                        else {section["label"]: truth[section["key"]] for section in sections}
                    ),
                ),
            ),
        )
    )

    payload = _report_payload(triad, spec=spec)
    payload["geometry"] = geometry
    payload["demonstration"] = demonstration
    payload["truth"] = truth
    payload["normalized"] = [s["f"] / total for s in sections] if total > 0 else None
    payload["sections"] = [
        {k: v for k, v in section.items() if k != "report"} for section in sections
    ]

    notes = [
        "Every section is Kearns' diffractogram route: its reflections sample basal tilt "
        "unevenly, which is why the triad does not sum exactly to one.",
        "A reflection with a random intensity worth measuring but no detected peak is counted "
        "as zero intensity, because the absence of a peak is information about the texture; "
        "check such rows against the diffractogram.",
        "Absorption and defocusing differ between sections of different geometry; a measured "
        "random standard of each geometry corrects for them, calculated intensities do not.",
    ]
    for section in sections:
        notes.extend(f"{section['name'].capitalize()}: {note}" for note in section["notes"][:1])

    result = AppResult(
        title=f"Kearns parameters of {spec.name} from three sections",
        summary=summary,
        table=ResultTable(
            columns=columns, rows=tuple(rows), caption="The triad, raw and normalised."
        ),
        data=payload,
        inputs={
            "phase": spec.to_json(),
            "geometry": geometry,
            "radiation": str(request["radiation"]),
            "random_source": str(request["random_source"]),
            "intensity_measure": str(request["intensity_measure"]),
            "two_theta_range_deg": [low, high],
            "match_tolerance_deg": float(request["match_tolerance_deg"]),
            "min_random_percent": float(request["min_random_percent"]),
            "prominence_sigma": float(request["prominence_sigma"]),
            "normalization": str(request["normalization"]),
            "bin_width_deg": float(request["bin_width_deg"]),
            "files": {s["key"]: s["file"] for s in sections},
        },
        notes=tuple(notes),
        citations=(_CITATION_KEARNS, _CITATION_MANI),
        stages=tuple(stages),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="kearns.example.three_sections",
            title="A zirconium tube, three sections at once",
            panel="kearns",
            summary="f along axial, radial and transverse from three scans of a known texture.",
            teaches=(
                "Read the three diffractograms first: the radial section's (0002) is strong and "
                "its (11-20) weak, the transverse section the reverse. Then open a reflection "
                "table and follow one row from the fitted area, through its random intensity, "
                "to the pole density at its tilt, and a quadrature table whose contributions sum "
                "to f. The triad sums close to, but not exactly, one - three independent "
                "measurements - and the exact values of the model texture are in the table to "
                "read the method against."
            ),
            operation="kearns.from_three_sections",
            request={"phase": {"builtin": "zr_hcp"}, "geometry": "tube"},
        ),
    )
)
