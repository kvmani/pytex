# ruff: noqa: RUF001
"""Powder X-ray diffraction for the shared web and desktop workbench.

The application layer does not implement diffraction physics.  It validates a
human-scale request, calls :func:`pytex.diffraction.xrd.generate_xrd_pattern`,
and turns the resulting reflection objects and sampled profile into the common
``AppResult`` contract used by tables, hover cards, and exports.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Literal, cast

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.logbook import APP_LOG
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.app.uploads import uploaded_file
from pytex.diffraction.rietveld import _scaled_phase, refine_rietveld
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_background import estimate_background
from pytex.diffraction.xrd_corrections import specimen_displacement_shift_deg
from pytex.diffraction.xrd_instrument import (
    InstrumentBroadening,
    calibrate_instrument_broadening,
    deconvolve_instrument_width,
    scherrer_size_nm,
    williamson_hall,
)
from pytex.diffraction.xrd_lattice_parameter import (
    crystal_system_of,
    determine_lattice_parameters_from_pattern,
    extrapolation_values,
)
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern, read_powder_pattern
from pytex.diffraction.xrd_phase_identification import identify_phase_from_pattern

__all__: tuple[str, ...] = ()

_CITATION_CULLITY = "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Chs. 3–4."
_CITATION_BEARDEN = "Bearden, Rev. Mod. Phys. 39 (1967) 78, doi:10.1103/RevModPhys.39.78."

_RADIATION = {
    "cu_ka_doublet": RadiationSpec.cu_ka_doublet,
    "cu_ka": RadiationSpec.cu_ka,
    "mo_ka_doublet": RadiationSpec.mo_ka_doublet,
    "co_ka_doublet": RadiationSpec.co_ka,
}

_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column("two_theta_deg", "2θ", units="°", numeric=True, digits=4),
    Column("d_angstrom", "d", units="Å", numeric=True, digits=5),
    Column(
        "relative_intensity",
        "Relative intensity",
        numeric=True,
        digits=4,
        help_text="Kinematic integrated intensity, normalized to the strongest Kα1 reflection.",
    ),
    Column("multiplicity", "Multiplicity", numeric=True),
    Column("structure_factor_amplitude", "|F|", numeric=True, digits=4),
    Column("lorentz_polarization", "L·P", numeric=True, digits=4),
)


def _powder_label(indices: tuple[int, int, int], *, spec: Any) -> str:
    """Use the conventional descending positive representative for cubic peaks."""

    display_indices = indices
    if spec.crystal_system == "cubic":
        display_indices = cast(
            tuple[int, int, int],
            tuple(sorted((abs(value) for value in indices), reverse=True)),
        )
    return plane_label(display_indices, spec=spec)


@REGISTRY.operation(
    "xrd.powder_pattern",
    title="Powder XRD pattern",
    summary="Structure-aware powder peaks and a broadened diffractogram with attributable indices.",
    help_text=(
        "Simulates a powder X-ray diffractogram from the selected phase's canonical lattice, "
        "symmetry and atomic basis. Peak positions follow Bragg's law; systematic absences and "
        "structure-factor contrast come from the structure; multiplicity and the powder "
        "Lorentz–polarization factor contribute to integrated intensity.\n\n"
        "Choose a single Kα line for a clean teaching pattern or a Kα1/Kα2 doublet to reproduce "
        "the high-angle splitting of laboratory data. Gaussian and pseudo-Voigt profiles model "
        "displayed peak breadth without changing the underlying reflection list.\n\n"
        "This is a kinematic, background-free simulation for indexing, phase-identification "
        "teaching and method development. It is not a Rietveld refinement, quantitative phase "
        "analysis, or a calibrated instrument response."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "The crystalline phase whose cell, symmetry and atomic sites generate the pattern. "
                "Built-in phases are ready-to-run; a custom phase must include an appropriate "
                "basis "
                "for structure-sensitive intensities."
            ),
            builtin="ni_fcc",
        ),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text=(
                "Laboratory characteristic radiation. Doublet choices include the weaker Kα2 line; "
                "Mo radiation has a shorter wavelength and moves peaks to lower angles."
            ),
            options=(
                ("cu_ka_doublet", "Cu Kα1/Kα2", "Common laboratory copper doublet."),
                ("cu_ka", "Cu Kα (single averaged line)", "One copper line without splitting."),
                ("mo_ka_doublet", "Mo Kα1/Kα2", "Short-wavelength molybdenum doublet."),
                ("co_ka_doublet", "Co Kα1/Kα2", "Useful for reducing Fe fluorescence."),
            ),
            default="cu_ka_doublet",
        ),
        NumberParameter(
            name="two_theta_min_deg",
            label="Start angle",
            help_text="Lower edge of the displayed and enumerated 2θ interval.",
            units="° 2θ",
            default=20.0,
            minimum=0.0,
            maximum=175.0,
            group="Scan",
            row="Angular range",
        ),
        NumberParameter(
            name="two_theta_max_deg",
            label="End angle",
            help_text="Upper edge of the displayed and enumerated 2θ interval.",
            units="° 2θ",
            default=120.0,
            minimum=1.0,
            maximum=180.0,
            group="Scan",
            row="Angular range",
        ),
        ChoiceParameter(
            name="profile",
            label="Peak profile",
            help_text=(
                "Gaussian is compact; pseudo-Voigt mixes Gaussian and Lorentzian tails and more "
                "closely resembles many laboratory peak shapes."
            ),
            options=(
                ("gaussian", "Gaussian", "Compact symmetric profile."),
                ("pseudo_voigt", "Pseudo-Voigt", "Gaussian–Lorentzian mixture with tails."),
            ),
            default="gaussian",
            group="Profile",
        ),
        NumberParameter(
            name="fwhm_deg",
            label="Peak FWHM",
            help_text=(
                "Constant full width at half maximum applied to every reflection. It is a display "
                "profile here, not a crystallite-size or microstrain refinement."
            ),
            units="° 2θ",
            default=0.15,
            minimum=0.01,
            maximum=5.0,
            group="Profile",
        ),
        NumberParameter(
            name="pseudo_voigt_eta",
            label="Lorentzian fraction η",
            help_text="Pseudo-Voigt mixing: 0 is Gaussian and 1 is Lorentzian.",
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            group="Profile",
        ),
        ChoiceParameter(
            name="intensity_model",
            label="Scattering model",
            help_text=(
                "Tabulated angle-dependent X-ray form factors are preferred. Constant atomic "
                "number is the simpler legacy proxy; unit amplitude isolates geometry and "
                "multiplicity."
            ),
            options=(
                (
                    "xray_tabulated",
                    "Tabulated X-ray form factors",
                    "Angle-dependent scattering factors; preferred for realistic contrast.",
                ),
                (
                    "xray_atomic_number",
                    "Constant atomic-number proxy",
                    "Legacy approximation with no angular falloff.",
                ),
                ("unit", "Unit structure factor", "Geometry and multiplicity only."),
            ),
            default="xray_tabulated",
            advanced=True,
        ),
        NumberParameter(
            name="resolution_deg",
            label="Angular sampling",
            help_text=(
                "Spacing of the continuous profile grid; smaller values produce more samples."
            ),
            units="° 2θ",
            default=0.02,
            minimum=0.005,
            maximum=0.25,
            advanced=True,
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest absolute h, k or l enumerated; raise it for wide high-angle scans.",
            default=6,
            minimum=1,
            maximum=12,
            advanced=True,
        ),
    ),
    returns="One row per Kα1 reflection family plus the sampled normalized intensity profile.",
    panel="xrd",
    citations=(_CITATION_CULLITY, _CITATION_BEARDEN),
    tags=("XRD", "powder diffraction", "Bragg", "phase identification", "K alpha"),
)
def _powder_pattern(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    minimum = float(request["two_theta_min_deg"])
    maximum = float(request["two_theta_max_deg"])
    if minimum >= maximum:
        raise InvalidInputError(
            "The start angle must be smaller than the end angle.",
            field="two_theta_max_deg",
            hint="Increase the end angle or reduce the start angle.",
        )

    radiation_key = str(request["radiation"])
    radiation = _RADIATION[radiation_key]()
    try:
        pattern = generate_xrd_pattern(
            phase,
            radiation=radiation,
            two_theta_range_deg=(minimum, maximum),
            resolution_deg=float(request["resolution_deg"]),
            max_index=int(request["max_index"]),
            intensity_model=cast(
                Literal["xray_atomic_number", "xray_tabulated", "unit"],
                request["intensity_model"],
            ),
            broadening_fwhm_deg=float(request["fwhm_deg"]),
            profile=cast(Literal["gaussian", "pseudo_voigt"], request["profile"]),
            pseudo_voigt_eta=float(request["pseudo_voigt_eta"]),
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The powder pattern could not be generated: {error}",
            hint="Check the angular range, profile width, sampling and phase definition.",
        ) from error
    if not pattern.reflections:
        raise InvalidInputError(
            "No reflections fall inside this angular window.",
            field="two_theta_max_deg",
            hint="Widen the scan range, change the radiation, or raise the index limit.",
        )

    strongest = max(reflection.intensity for reflection in pattern.reflections) or 1.0
    rows = tuple(
        {
            "hkl_label": _powder_label(tuple(reflection.miller_indices), spec=spec),
            "h": int(reflection.miller_indices[0]),
            "k": int(reflection.miller_indices[1]),
            "l": int(reflection.miller_indices[2]),
            "two_theta_deg": float(reflection.two_theta_deg),
            "d_angstrom": float(reflection.d_spacing_angstrom),
            "relative_intensity": float(reflection.intensity / strongest),
            "integrated_intensity": float(reflection.intensity),
            "multiplicity": int(reflection.multiplicity),
            "structure_factor_amplitude": float(reflection.structure_factor_amplitude),
            "lorentz_polarization": float(reflection.lorentz_polarization_factor or 0.0),
        }
        for reflection in pattern.reflections
    )

    profile_name = "pseudo-Voigt" if request["profile"] == "pseudo_voigt" else "Gaussian"
    doublet = radiation.kalpha2_wavelength_angstrom is not None
    result = AppResult(
        title=f"Powder XRD of {spec.name}",
        summary=(
            f"{len(rows)} Kα1 reflection families from {minimum:g}° to {maximum:g}° 2θ using "
            f"{radiation.name} radiation and a {float(request['fwhm_deg']):g}° {profile_name} "
            f"profile. Peak positions follow Bragg's law; displayed intensity is normalized to "
            "the profile maximum."
        ),
        table=ResultTable(
            columns=_COLUMNS,
            rows=rows,
            caption=f"Indexed Kα1 powder reflections of {spec.name}.",
        ),
        data={
            "two_theta_deg": pattern.two_theta_grid_deg.tolist(),
            "intensity": pattern.intensity_grid.tolist(),
            "reflections": list(rows),
            "columns": [column.to_json() for column in _COLUMNS],
            "radiation_name": radiation.name,
            "wavelength_angstrom": radiation.wavelength_angstrom,
            "kalpha2_wavelength_angstrom": radiation.kalpha2_wavelength_angstrom,
            "doublet": doublet,
            "phase_name": spec.name,
        },
        inputs={
            "phase": spec.to_json(),
            "radiation": radiation_key,
            "two_theta_min_deg": minimum,
            "two_theta_max_deg": maximum,
            "profile": request["profile"],
            "fwhm_deg": float(request["fwhm_deg"]),
            "pseudo_voigt_eta": float(request["pseudo_voigt_eta"]),
            "intensity_model": request["intensity_model"],
            "resolution_deg": float(request["resolution_deg"]),
            "max_index": int(request["max_index"]),
        },
        notes=(
            "The simulation is kinematic and background-free. It omits absorption, fluorescence, "
            "specimen displacement, axial divergence and a calibrated instrument response.",
            "The reflection table lists the primary Kα1 families. When a doublet is selected, the "
            "weaker Kα2 contribution is present in the continuous profile but is not duplicated "
            "in the table.",
            "Relative heights are not suitable for quantitative phase analysis or Rietveld "
            "refinement.",
        ),
        citations=(_CITATION_CULLITY, _CITATION_BEARDEN),
    )
    return result.to_json()


# ---------------------------------------------------------------------------
# Quantitative analysis: background, refinement, and size/strain.
#
# These three operations work on a *measured* profile, which the simulation
# operation above does not. Where the scan comes from is therefore the first
# question each of them has to answer, and it is answered by one shared control
# rather than three, so a scan pasted into one view means the same thing in the
# next.
# ---------------------------------------------------------------------------

_CITATION_RIETVELD_METHOD = (
    "Rietveld, J. Appl. Crystallogr. 2 (1969) 65, doi:10.1107/S0021889869006558."
)
_CITATION_TOBY_R = "Toby, Powder Diffr. 21 (2006) 67, doi:10.1154/1.2179804."
_CITATION_SNIP_METHOD = (
    "Ryan et al., Nucl. Instrum. Methods B 34 (1988) 396, doi:10.1016/0168-583X(88)90063-8."
)
_CITATION_CAGLIOTI_UVW = (
    "Caglioti, Paoletti & Ricci, Nucl. Instrum. 3 (1958) 223, doi:10.1016/0369-643X(58)90029-X."
)
_CITATION_WH = "Williamson & Hall, Acta Metall. 1 (1953) 22, doi:10.1016/0001-6160(53)90006-6."

#: How a demonstration scan departs from the ideal. These are stated here, used
#: by the synthesizer, and quoted back in every result that uses it, so a
#: simulated measurement can never be mistaken for a real one or its answer
#: mistaken for a discovery.
_DEMO_LATTICE_SCALE = 1.003
_DEMO_ZERO_SHIFT_DEG = 0.05
_DEMO_FWHM_DEG = 0.14
_DEMO_PEAK_COUNTS = 20000.0
_DEMO_BACKGROUND_COUNTS = 150.0

#: The most rows any profile table carries. The full curve always travels in
#: `data` for the plot; the table is what a reader scrolls and what the CSV
#: export writes, and thirty thousand rows serves neither.
_MAX_PROFILE_ROWS = 1200


PATTERN_FILE_SUFFIXES: tuple[str, ...] = (".xy", ".xrdml", ".csv", ".dat", ".txt")


def _scan_parameters(*, group: str = "Measurement") -> tuple[Any, ...]:
    """Return the shared "where does the scan come from" controls.

    One declaration, used by every analysis operation, so the meaning of an
    experimental scan cannot drift between views.
    """

    return (
        ChoiceParameter(
            name="data_source",
            label="Scan source",
            help_text=(
                "Where the measured profile comes from. A demonstration scan is generated from "
                "the selected phase with a deliberate cell dilation, zero-point error, peak "
                "width and curved background, then given Poisson counting noise. An experimental "
                "pattern file (.xy, .xrdml, .csv, .dat) can be opened and loaded directly."
            ),
            options=(
                (
                    "demonstration",
                    "Generate a demonstration scan",
                    "Synthetic, with known answers, for learning and method development.",
                ),
                (
                    "file",
                    "Use an experimental pattern file",
                    "An experimental powder diffractogram (.xy, .xrdml, .csv, .dat).",
                ),
                (
                    "paste",
                    "Use the pasted scan",
                    "Two columns of your own data: 2θ and intensity (legacy fallback).",
                ),
            ),
            default="demonstration",
            group=group,
        ),
        ObjectParameter(
            name="scan_file",
            label="Pattern file",
            help_text=(
                "An experimental powder diffraction pattern (.xy, .xrdml, .csv, .dat). "
                "Opened through **Open a pattern file** in the workbench rail."
            ),
            required=False,
            group=group,
        ),
        TextParameter(
            name="scan",
            label="Pasted scan",
            help_text=(
                "Two numbers per line: 2θ in degrees, then intensity. Blank lines and "
                "lines beginning with `#` are ignored, and commas count as separators so a "
                "pasted `.xy` or CSV export works unchanged.\n\n"
                "Paste the **raw** scan. Do not subtract a background first: the refinement "
                "fits the background jointly with everything else, and subtracting beforehand "
                "discards the correlation between background and scale that the reported "
                "uncertainties depend on."
            ),
            multiline=True,
            required=False,
            default="",
            placeholder="30.00  152\n30.02  148\n30.04  157",
            advanced=True,
            group=group,
        ),
        IntegerParameter(
            name="demonstration_seed",
            label="Demonstration noise seed",
            help_text=(
                "Seed for the counting noise of the demonstration scan, so a result is "
                "reproducible. Ignored when an experimental pattern is used."
            ),
            default=20260905,
            minimum=0,
            maximum=2**31 - 1,
            advanced=True,
            group=group,
            field_width="short",
        ),
    )


def _parse_scan(text: str) -> tuple[np.ndarray, np.ndarray]:
    """Read a pasted two-column scan into angles and intensities."""

    angles: list[float] = []
    intensities: list[float] = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.replace(",", " ").split()
        if len(fields) < 2:
            raise InvalidInputError(
                f"Line {number} of the pasted scan has only one number.",
                field="scan",
                hint="Each line needs a 2θ value and an intensity, separated by "
                "whitespace or a comma.",
            )
        try:
            angle, intensity = float(fields[0]), float(fields[1])
        except ValueError as error:
            raise InvalidInputError(
                f"Line {number} of the pasted scan is not two numbers: {stripped!r}.",
                field="scan",
                hint="Remove header text, or comment it out with a leading `#`.",
            ) from error
        angles.append(angle)
        intensities.append(intensity)
    if len(angles) < 20:
        raise InvalidInputError(
            f"The pasted scan holds {len(angles)} points; at least 20 are needed to analyse a "
            "profile.",
            field="scan",
            hint="Paste the whole scan, not a summary of its peaks.",
        )
    return np.asarray(angles, dtype=float), np.asarray(intensities, dtype=float)


def _measured_from_request(
    request: dict[str, Any], phase: Any, radiation: RadiationSpec
) -> tuple[MeasuredPowderPattern, bool]:
    """Return the profile to analyse, and whether it was generated rather than measured."""

    data_source = str(request.get("data_source", "demonstration"))
    scan_file = request.get("scan_file")

    if data_source == "file" or scan_file:
        if not scan_file:
            raise InvalidInputError(
                "No experimental pattern file was provided.",
                field="scan_file",
                hint=(
                    "Choose an experimental pattern file (.xy, .xrdml, .csv, .dat) using the "
                    "pattern loader, or switch the scan source back to the demonstration scan."
                ),
            )
        with uploaded_file(scan_file, field="scan_file", suffixes=PATTERN_FILE_SUFFIXES) as (
            path,
            name,
        ):
            try:
                pattern = read_powder_pattern(path, name=name, radiation=radiation)
            except Exception as error:
                raise InvalidInputError(
                    f"The pattern file {name} could not be read: {error}",
                    field="scan_file",
                    hint=(
                        "Verify that the file is a valid 2-column .xy/.csv or a PANalytical .xrdml "
                        "powder scan with increasing 2θ angles."
                    ),
                ) from error
        APP_LOG.info(
            f"Loaded experimental pattern '{pattern.name}' ({len(pattern)} points).",
            source="xrd",
            detail={
                "points": len(pattern),
                "two_theta_min": float(pattern.two_theta_deg[0]),
                "two_theta_max": float(pattern.two_theta_deg[-1]),
            },
        )
        return pattern, False

    if data_source == "paste":
        text = str(request.get("scan") or "")
        if not text.strip():
            raise InvalidInputError(
                "No scan was pasted.",
                field="scan",
                hint="Paste two columns of data, or switch the scan source back to the "
                "demonstration scan.",
            )
        angles, intensities = _parse_scan(text)
        try:
            return (
                MeasuredPowderPattern(
                    name="pasted scan",
                    two_theta_deg=angles,
                    intensity=intensities,
                    radiation=radiation,
                ),
                False,
            )
        except ValueError as error:
            raise InvalidInputError(
                f"The pasted scan could not be read as a profile: {error}",
                field="scan",
                hint="Angles must increase strictly down the column and intensities must be "
                "finite and non-negative.",
            ) from error

    dilated = _scaled_phase(phase, _DEMO_LATTICE_SCALE)
    try:
        ideal = generate_xrd_pattern(
            dilated,
            radiation=radiation,
            two_theta_range_deg=(30.0, 130.0),
            resolution_deg=0.02,
            broadening_fwhm_deg=_DEMO_FWHM_DEG,
            intensity_model="xray_tabulated",
        )
    except ValueError as error:
        raise InvalidInputError(
            f"A demonstration scan could not be generated for this phase: {error}",
            field="phase",
            hint="Choose a phase with an atomic basis, or open an experimental pattern instead.",
        ) from error
    angles = ideal.two_theta_grid_deg + _DEMO_ZERO_SHIFT_DEG
    noiseless = (
        _DEMO_PEAK_COUNTS * ideal.intensity_grid
        + _DEMO_BACKGROUND_COUNTS
        + 400.0 * np.exp(-0.5 * ((angles - 34.0) / 6.0) ** 2)
    )
    generator = np.random.default_rng(int(request["demonstration_seed"]))
    return (
        MeasuredPowderPattern(
            name="demonstration scan (synthetic)",
            two_theta_deg=angles,
            intensity=generator.poisson(noiseless).astype(float),
            radiation=radiation,
            synthetic=True,
        ),
        True,
    )


def _demonstration_notes(phase_a: float) -> tuple[str, ...]:
    """Return the notes every result built on a generated scan must carry."""

    return (
        "This scan was generated, not measured. It is a synthetic profile of the selected phase "
        f"with the cell dilated by {_DEMO_LATTICE_SCALE:g} (a = "
        f"{phase_a * _DEMO_LATTICE_SCALE:.5f} \u00c5 against a tabulated "
        f"{phase_a:.5f} \u00c5), a {_DEMO_ZERO_SHIFT_DEG:g}\u00b0 detector zero error, a "
        f"{_DEMO_FWHM_DEG:g}\u00b0 peak width and a curved background, then Poisson counting "
        "noise.",
        "Because those departures are known, the numbers below can be checked rather than "
        "believed. That is what the demonstration scan is for; switch the scan source to a "
        "pasted scan to analyse real data.",
    )


def _decimate(count: int) -> np.ndarray:
    """Return row indices for a table that stays readable on a long scan."""

    if count <= _MAX_PROFILE_ROWS:
        return np.arange(count)
    return np.unique(np.linspace(0, count - 1, _MAX_PROFILE_ROWS).round().astype(int))


_BACKGROUND_COLUMNS = (
    Column("two_theta_deg", "2\u03b8", units="\u00b0", numeric=True, digits=4),
    Column("observed", "Observed", numeric=True, digits=2),
    Column("background", "Background", numeric=True, digits=2),
    Column(
        "subtracted",
        "Observed \u2212 background",
        numeric=True,
        digits=2,
        help_text="Clipped at zero: where the estimate crosses above the data, the model is "
        "what went wrong, not the measurement.",
    ),
)


@REGISTRY.operation(
    "xrd.background",
    title="Background estimation",
    summary="Separate the slowly varying background of a raw scan from its Bragg intensity.",
    help_text=(
        "A raw diffractogram is Bragg intensity plus a background built from air scatter, "
        "sample fluorescence, incoherent scattering, the holder and detector noise. Every "
        "quantitative use of the scan \u2014 integrated intensities, peak widths, refinement "
        "\u2014 depends on telling the two apart, and the separation is a modelling choice "
        "rather than a measurement.\n\n"
        "SNIP clips each point against the mean of its neighbours at a growing separation, in a "
        "domain that makes the clip insensitive to count level. It assumes only that the "
        "background varies more slowly with angle than the peaks do, so it follows curved and "
        "structured backgrounds that no low-order polynomial can. Use it to see what the "
        "background is doing.\n\n"
        "The Chebyshev fit discards points lying above the current curve and refits, "
        "repeatedly. The asymmetry is deliberate: peaks are one-sided excursions, so a "
        "symmetric rejection would drag the curve up into the peak feet. It yields the "
        "coefficients a refinement carries.\n\n"
        "Neither estimator is told where the peaks are, which is what allows a background to be "
        "used to find them. Estimate a background to look at it \u2014 but do not subtract it "
        "before a refinement, which fits it jointly."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "Only used to generate the demonstration scan. A pasted scan is analysed "
                "without reference to any phase, which is the point: the background estimate "
                "must not depend on knowing what the material is."
            ),
            builtin="ni_fcc",
        ),
        *_scan_parameters(),
        ChoiceParameter(
            name="method",
            label="Estimator",
            help_text=(
                "SNIP is non-parametric and follows curved backgrounds. The Chebyshev fit is "
                "parametric and produces coefficients a refinement can carry."
            ),
            options=(
                ("snip", "SNIP peak clipping", "Non-parametric; follows structured backgrounds."),
                (
                    "chebyshev",
                    "Chebyshev polynomial",
                    "Parametric; the family used inside whole-profile refinement.",
                ),
            ),
            default="snip",
            group="Estimator",
        ),
        NumberParameter(
            name="half_window_deg",
            label="SNIP clipping window",
            help_text=(
                "Half-width of the clipping window. Set it comfortably wider than the broadest "
                "peak and comfortably narrower than the curvature of the background: too small "
                "leaves peak feet behind as background, too large flattens real background "
                "structure."
            ),
            units="\u00b0 2\u03b8",
            default=2.0,
            minimum=0.05,
            maximum=20.0,
            group="Estimator",
        ),
        IntegerParameter(
            name="degree",
            label="Chebyshev degree",
            help_text=(
                "Polynomial order. Four to eight covers most laboratory scans; a high degree "
                "will absorb genuine broad features, including an amorphous halo you may want "
                "to see."
            ),
            default=6,
            minimum=0,
            maximum=20,
            group="Estimator",
        ),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text="Recorded with the profile, and used to generate a demonstration scan.",
            options=(
                ("cu_ka", "Cu K\u03b1 (single averaged line)", "One copper line."),
                ("cu_ka_doublet", "Cu K\u03b11/K\u03b12", "Common laboratory copper doublet."),
                ("co_ka_doublet", "Co K\u03b11/K\u03b12", "Reduces Fe fluorescence."),
                ("mo_ka_doublet", "Mo K\u03b11/K\u03b12", "Short-wavelength molybdenum."),
            ),
            default="cu_ka",
            advanced=True,
        ),
    ),
    returns="The estimated background and the background-subtracted profile on the measured grid.",
    panel="xrd",
    citations=(_CITATION_SNIP_METHOD, _CITATION_CULLITY),
    tags=("XRD", "background", "SNIP", "Chebyshev", "experimental data"),
)
def _background(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    radiation = _RADIATION[str(request["radiation"])]()
    measured, generated = _measured_from_request(request, phase, radiation)
    method = cast(Literal["snip", "chebyshev"], request["method"])
    try:
        estimate = estimate_background(
            measured,
            method=method,
            half_window_deg=float(request["half_window_deg"]),
            degree=int(request["degree"]),
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The background could not be estimated: {error}",
            field="half_window_deg" if method == "snip" else "degree",
            hint="Widen or narrow the clipping window, or lower the polynomial degree.",
        ) from error

    subtracted = estimate.subtracted_intensity()
    rows_index = _decimate(estimate.point_count)
    rows = tuple(
        {
            "two_theta_deg": float(estimate.two_theta_deg[index]),
            "observed": float(estimate.observed_intensity[index]),
            "background": float(estimate.background[index]),
            "subtracted": float(subtracted[index]),
        }
        for index in rows_index
    )
    caption = f"Background of {measured.name} by the {method} estimator."
    if rows_index.size < estimate.point_count:
        caption += (
            f" Sampled to {rows_index.size} of {estimate.point_count} points for display; the "
            "plot and the JSON export carry every point."
        )
    notes = [
        "The background is a modelling choice, not a measurement. An estimate that clips into "
        "weak reflections removes intensity a refinement will then be unable to account for.",
        "Do not subtract a background before a Rietveld refinement. The refinement fits the "
        "background jointly, and subtracting first discards the correlation between background "
        "and scale that the reported uncertainties depend on.",
    ]
    if generated:
        notes = [*_demonstration_notes(phase.lattice.a), *notes]
    result = AppResult(
        title=f"Background of {measured.name}",
        summary=(
            f"The {method} estimator assigns "
            f"{100.0 * estimate.background_fraction:.1f}% of the total measured signal to "
            f"background over {estimate.point_count} points from "
            f"{estimate.two_theta_deg[0]:.2f}\u00b0 to {estimate.two_theta_deg[-1]:.2f}\u00b0 "
            f"2\u03b8, with a level between {float(np.min(estimate.background)):.4g} and "
            f"{float(np.max(estimate.background)):.4g}."
        ),
        table=ResultTable(columns=_BACKGROUND_COLUMNS, rows=rows, caption=caption),
        data={
            "two_theta_deg": estimate.two_theta_deg.tolist(),
            "observed": estimate.observed_intensity.tolist(),
            "background": estimate.background.tolist(),
            "subtracted": subtracted.tolist(),
            "columns": [column.to_json() for column in _BACKGROUND_COLUMNS],
            "method": method,
            "background_fraction": estimate.background_fraction,
            "synthetic": generated,
            "describe": estimate.describe(),
            "phase_name": spec.name,
        },
        inputs={
            "phase": spec.to_json(),
            "data_source": request["data_source"],
            "method": method,
            "half_window_deg": float(request["half_window_deg"]),
            "degree": int(request["degree"]),
            "radiation": request["radiation"],
            "demonstration_seed": int(request["demonstration_seed"]),
        },
        notes=tuple(notes),
        citations=(_CITATION_SNIP_METHOD,),
    )
    return result.to_json()


_REFINEMENT_COLUMNS = (
    Column("parameter", "Parameter"),
    Column("value", "Refined value", help_text="Uncertainty on the last quoted digits."),
    Column("initial", "Started at", numeric=True, digits=6),
    Column("shift", "Moved by", numeric=True, digits=6),
    Column("units", "Units"),
    Column("meaning", "What it means"),
)


@REGISTRY.operation(
    "xrd.rietveld",
    title="Rietveld refinement",
    summary="Fit the whole measured profile and report how well the model accounts for it.",
    help_text=(
        "Rietveld's method fits the measured profile point by point against a pattern "
        "calculated from a structural model, rather than reducing the scan to integrated "
        "intensities first. Overlapping reflections cannot be separated reliably, but they can "
        "be calculated, so the overlap never has to be resolved.\n\n"
        "Refine the detector zero alongside the cell. A zero-point error and a cell dilation "
        "both move every peak, so holding the zero at nought lets the cell silently absorb it "
        "and return a confidently wrong lattice parameter with a small uncertainty attached.\n\n"
        "Read the fit from R_wp against R_exp, not from R_wp alone: R_exp is the value R_wp "
        "would take if the only remaining misfit were counting noise, so their ratio \u2014 the "
        "goodness of fit \u2014 is what removes the flattering effect of a large background. A "
        "goodness of fit near 1 means the model explains the data down to the noise; below 1 "
        "means it is following noise. Then read the difference curve, which is the most "
        "informative single output a refinement has.\n\n"
        "This refines the profile, the cell dilation, the zero point, the peak width, the "
        "texture strength and the background against a known structure. It does not refine "
        "atomic coordinates, occupancies or anisotropic displacement parameters: those need "
        "constraints and restraints this surface does not offer, and refining them without "
        "that apparatus buys a lower R factor and a structure nobody should publish."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "The structural model to test and adjust. It needs an atomic basis for "
                "structure-sensitive intensities; without one only geometry and multiplicity "
                "contribute, and the intensity ratios will not match any real scan."
            ),
            builtin="ni_fcc",
        ),
        *_scan_parameters(),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text="The radiation the scan was collected with. A wrong choice moves every "
            "calculated peak and the refinement will try to absorb it into the cell.",
            options=(
                ("cu_ka", "Cu K\u03b1 (single averaged line)", "One copper line."),
                ("cu_ka_doublet", "Cu K\u03b11/K\u03b12", "Common laboratory copper doublet."),
                ("co_ka_doublet", "Co K\u03b11/K\u03b12", "Reduces Fe fluorescence."),
                ("mo_ka_doublet", "Mo K\u03b11/K\u03b12", "Short-wavelength molybdenum."),
            ),
            default="cu_ka",
        ),
        BooleanParameter(
            name="refine_cell",
            label="Refine the cell",
            help_text=(
                "Vary an isotropic dilation of the unit cell. Symmetry-preserving in every "
                "crystal system, which is why it is the safe general cell parameter."
            ),
            default=True,
            group="What to refine",
        ),
        BooleanParameter(
            name="refine_zero_shift",
            label="Refine the detector zero",
            help_text=(
                "Vary the zero-point error. Refine it whenever you refine the cell: they are "
                "correlated, and the one held fixed is absorbed by the one that is not."
            ),
            default=True,
            group="What to refine",
        ),
        BooleanParameter(
            name="refine_width",
            label="Refine the peak width",
            help_text=(
                "Vary the angle-independent Caglioti W. A width held at the wrong value leaves "
                "a residual with the derivative shape of the peak at every reflection, which "
                "shows up in the Durbin\u2013Watson statistic long before it shows up in R_wp."
            ),
            default=True,
            group="What to refine",
        ),
        BooleanParameter(
            name="refine_width_trend",
            label="Refine the width's angular trend",
            help_text=(
                "Also vary Caglioti U and V, so the width may follow "
                "U tan\u00b2\u03b8 + V tan\u03b8 + W across the scan. Worth it on a wide "
                "scan; on a narrow one U, V and W cannot be told apart and the uncertainties "
                "will say so."
            ),
            default=False,
            group="What to refine",
        ),
        BooleanParameter(
            name="refine_texture",
            label="Refine texture strength",
            help_text=(
                "Vary a March\u2013Dollase coefficient about the stated axis. Below 1 is "
                "plate-like and enhances the preferred reflections; above 1 is needle-like and "
                "suppresses them."
            ),
            default=False,
            group="What to refine",
        ),
        IndicesParameter(
            name="texture_axis",
            label="Preferred orientation (hkl)",
            help_text=(
                "The plane whose normals cluster along the specimen axis \u2014 (001) for a "
                "basal-textured sheet, (100) for many pressed powders. A texture strength "
                "about an unstated axis has no meaning, so this is required before texture may "
                "be refined."
            ),
            default=(1, 1, 1),
            group="What to refine",
        ),
        IntegerParameter(
            name="background_degree",
            label="Background degree",
            help_text=(
                "Chebyshev order, refined jointly with everything else. Four to eight covers "
                "most laboratory scans. If the difference curve wanders slowly under the whole "
                "pattern and the Durbin\u2013Watson statistic is low while the peaks look "
                "well fitted, the background is too stiff \u2014 raise it."
            ),
            default=8,
            minimum=0,
            maximum=20,
            group="Background",
        ),
        NumberParameter(
            name="starting_fwhm_deg",
            label="Starting peak width",
            help_text="Where the width refinement begins, as a full width at half maximum.",
            units="\u00b0 2\u03b8",
            default=0.10,
            minimum=0.01,
            maximum=3.0,
            advanced=True,
            group="Background",
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest absolute h, k or l enumerated; raise it for wide high-angle scans.",
            default=6,
            minimum=1,
            maximum=12,
            advanced=True,
        ),
    ),
    returns=(
        "One row per refinement parameter, plus the observed, calculated, background and "
        "difference profiles and the agreement indices."
    ),
    panel="xrd",
    citations=(_CITATION_RIETVELD_METHOD, _CITATION_TOBY_R, _CITATION_CULLITY),
    tags=("XRD", "Rietveld", "refinement", "lattice parameter", "R factor", "experimental data"),
)
def _rietveld(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    radiation = _RADIATION[str(request["radiation"])]()
    measured, generated = _measured_from_request(request, phase, radiation)

    refine: list[str] = ["scale"]
    if bool(request["refine_cell"]):
        refine.append("lattice_scale")
    if bool(request["refine_zero_shift"]):
        refine.append("zero_shift_deg")
    if bool(request["refine_width"]):
        refine.append("caglioti_w")
    if bool(request["refine_width_trend"]):
        refine.extend(("caglioti_u", "caglioti_v"))
    texture_axis = tuple(int(value) for value in request["texture_axis"])
    if bool(request["refine_texture"]):
        refine.append("march_coefficient")

    try:
        result_object = refine_rietveld(
            measured,
            phase,
            radiation=radiation,
            instrument=InstrumentBroadening.ideal(float(request["starting_fwhm_deg"])),
            refine=tuple(refine),
            background_degree=int(request["background_degree"]),
            preferred_orientation_plane=cast(tuple[int, int, int], texture_axis),
            max_index=int(request["max_index"]),
            intensity_model="xray_tabulated",
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The refinement could not be run: {error}",
            hint="Check that the phase, the radiation and the scan belong together, and that "
            "the scan covers enough reflections to constrain what is being refined.",
        ) from error

    rows = tuple(
        {
            "parameter": parameter.name,
            "value": parameter.format() if parameter.refined else f"{parameter.value:.6g}",
            "initial": float(parameter.initial_value),
            "shift": float(parameter.shift),
            "units": parameter.units or "",
            "meaning": parameter.description
            + ("" if parameter.refined else " Held fixed in this refinement."),
        }
        for parameter in result_object.parameters
    )
    profile_index = _decimate(result_object.point_count)
    notes = [
        "This refinement varies scale, cell dilation, zero point, profile width, background and "
        "texture strength only. It does not refine atomic coordinates, occupancies or "
        "anisotropic displacement parameters, so it tests and adjusts a structural model rather "
        "than determining a structure.",
        "Standard uncertainties are precision estimates conditional on the model being right. "
        "They say nothing about whether it is right, and are routinely optimistic by a factor of "
        "two or three because neighbouring profile points are correlated.",
        "A parameter reported with an enormous uncertainty is one the data cannot determine "
        "separately from another; that is the intended signal, not a numerical failure.",
    ]
    if not result_object.converged:
        notes.insert(
            0,
            "The optimizer did not converge. The values below are wherever it halted and must "
            "not be quoted.",
        )
    if generated:
        notes = [*_demonstration_notes(phase.lattice.a), *notes]
    result = AppResult(
        title=f"Rietveld refinement of {spec.name}",
        summary=(
            f"R_wp = {100.0 * result_object.weighted_profile_r_factor:.3f}% against an R_exp of "
            f"{100.0 * result_object.expected_r_factor:.3f}%, a goodness of fit of "
            f"{result_object.goodness_of_fit:.3f} and R_Bragg = "
            f"{100.0 * result_object.bragg_r_factor:.3f}%, varying "
            f"{result_object.refined_parameter_count} parameters against "
            f"{result_object.point_count} observations. Refined cell: a = "
            f"{result_object.phase.lattice.a:.5f} \u00c5."
        ),
        table=ResultTable(
            columns=_REFINEMENT_COLUMNS,
            rows=rows,
            caption=f"Refinement parameters for {spec.name}.",
        ),
        data={
            "two_theta_deg": result_object.two_theta_deg[profile_index].tolist(),
            "observed": result_object.observed_intensity[profile_index].tolist(),
            "calculated": result_object.calculated_intensity[profile_index].tolist(),
            "background": result_object.background_intensity[profile_index].tolist(),
            "residual": result_object.residual_intensity[profile_index].tolist(),
            "reflections": [
                {
                    "hkl_label": _powder_label(
                        cast(
                            tuple[int, int, int],
                            tuple(int(value) for value in reflection.miller_indices),
                        ),
                        spec=spec,
                    ),
                    "two_theta_deg": float(reflection.two_theta_deg),
                    "d_angstrom": float(reflection.d_spacing_angstrom),
                }
                for reflection in result_object.reflections
            ],
            "columns": [column.to_json() for column in _REFINEMENT_COLUMNS],
            "profile_r_factor": result_object.profile_r_factor,
            "weighted_profile_r_factor": result_object.weighted_profile_r_factor,
            "expected_r_factor": result_object.expected_r_factor,
            "bragg_r_factor": result_object.bragg_r_factor,
            "goodness_of_fit": result_object.goodness_of_fit,
            "durbin_watson": result_object.durbin_watson,
            "weight_model": result_object.weight_model,
            "converged": result_object.converged,
            "refined_lattice_a": result_object.phase.lattice.a,
            "synthetic": generated,
            "describe": result_object.describe(),
            "phase_name": spec.name,
        },
        inputs={
            "phase": spec.to_json(),
            "data_source": request["data_source"],
            "radiation": request["radiation"],
            "refine": list(refine),
            "background_degree": int(request["background_degree"]),
            "starting_fwhm_deg": float(request["starting_fwhm_deg"]),
            "texture_axis": list(texture_axis),
            "max_index": int(request["max_index"]),
            "demonstration_seed": int(request["demonstration_seed"]),
        },
        notes=tuple(notes),
        citations=(_CITATION_RIETVELD_METHOD, _CITATION_TOBY_R),
    )
    return result.to_json()


_SIZE_STRAIN_COLUMNS = (
    Column("two_theta_deg", "2\u03b8", units="\u00b0", numeric=True, digits=4),
    Column("observed_fwhm_deg", "Measured FWHM", units="\u00b0", numeric=True, digits=5),
    Column("instrument_fwhm_deg", "Instrument FWHM", units="\u00b0", numeric=True, digits=5),
    Column(
        "sample_fwhm_deg",
        "Sample FWHM",
        units="\u00b0",
        numeric=True,
        digits=5,
        help_text="What remains after the instrumental width is deconvolved; the only width a "
        "size or strain may be quoted from.",
    ),
    Column(
        "scherrer_size_nm",
        "Scherrer size",
        units="nm",
        numeric=True,
        digits=3,
        help_text="Attributes all remaining broadening to size, so it is a lower bound whenever "
        "strain is present.",
    ),
    Column("abscissa", "4 sin\u03b8", numeric=True, digits=5),
    Column("ordinate", "\u03b2 cos\u03b8", units="rad", numeric=True, digits=6),
)


@REGISTRY.operation(
    "xrd.size_strain",
    title="Crystallite size and microstrain",
    summary="Calibrate the instrumental width from a standard, then separate size from strain.",
    help_text=(
        "A measured peak is wider than the sample makes it. Until the instrumental contribution "
        "is removed, a crystallite size read off a width is not an underestimate of the truth "
        "\u2014 it is a measurement of the diffractometer.\n\n"
        "So this takes two lists of widths. The first is from a line-profile standard (NIST SRM "
        "660 LaB\u2086, SRM 640 silicon, or any specimen whose own broadening is negligible), "
        "and fits the Caglioti resolution function U tan\u00b2\u03b8 + V tan\u03b8 + W to it. "
        "That fit is exactly linear in U, V and W, so it needs no starting values and cannot "
        "fail to converge. The second is from your specimen, and is deconvolved against that "
        "calibration.\n\n"
        "What is left contains size and strain together, and one peak cannot separate them "
        "because both widen it. Their angular dependences differ \u2014 size broadening as "
        "1/cos\u03b8, strain as tan\u03b8 \u2014 so several peaks can: fitting "
        "\u03b2 cos\u03b8 = K\u03bb/D + 4\u03b5 sin\u03b8 puts the size in the intercept "
        "and the strain in the slope.\n\n"
        "Watch the R\u00b2. A poor straight line means the uniform-deformation assumption does "
        "not hold for this specimen \u2014 most often anisotropic broadening \u2014 and the "
        "size and strain should not then be quoted."
    ),
    parameters=(
        TextParameter(
            name="standard_peaks",
            label="Standard peak widths",
            help_text=(
                "Two numbers per line: 2\u03b8 in degrees, then the measured FWHM in degrees, "
                "for a line-profile standard. At least three peaks, spanning as wide an "
                "angular range as possible \u2014 three peaks within ten degrees of each other "
                "determine the parabola no better than one does.\n\n"
                "These are the LaB\u2086 lines at Cu K\u03b1, at the widths a representative "
                "Bragg\u2013Brentano instrument gives them, paired below with a specimen "
                "whose answer is known in advance: the untouched form returns exactly "
                "25 nm and 0.2% microstrain. Replace both columns with your own."
            ),
            multiline=True,
            default=(
                "21.36  0.0757\n30.38  0.0754\n37.44  0.0754\n43.51  0.0756\n48.96  0.0760\n5"
                "3.99  0.0764\n63.22  0.0777\n67.55  0.0785\n71.75  0.0794\n75.84  0.0806\n79"
                ".87  0.0818\n83.85  0.0833\n87.79  0.0850"
            ),
            placeholder="21.36  0.0784",
            group="Standard",
        ),
        TextParameter(
            name="sample_peaks",
            label="Sample peak widths",
            help_text=(
                "Two numbers per line: 2\u03b8 in degrees, then the measured FWHM in degrees, "
                "for the specimen under study, measured the same way as the standard's. At "
                "least two peaks, and in practice at least four across a wide range: the "
                "intercept of a line fitted over a narrow range is almost unconstrained."
            ),
            multiline=True,
            default=(
                "21.36  0.4167\n30.38  0.4600\n37.44  0.4966\n43.51  0.5305\n48.96  0.5630\n5"
                "3.99  0.5950\n63.22  0.6598\n67.55  0.6933\n71.75  0.7280\n75.84  0.7642\n79"
                ".87  0.8023\n83.85  0.8428\n87.79  0.8861"
            ),
            placeholder="21.36  0.4167",
            group="Specimen",
        ),
        NumberParameter(
            name="wavelength_angstrom",
            label="Wavelength",
            help_text="The radiation the widths were measured at.",
            units="\u00c5",
            default=1.5406,
            minimum=0.1,
            maximum=10.0,
            group="Specimen",
            symbol="wavelength",
        ),
        ChoiceParameter(
            name="mode",
            label="Deconvolution model",
            help_text=(
                "How the instrumental width is removed, which depends on peak shape. Gaussian "
                "subtracts in quadrature and suits a well-behaved laboratory instrument; "
                "Lorentzian subtracts directly and suits size-dominated broadening, which is "
                "Lorentzian-like; pseudo-Voigt separates the two components and treats each in "
                "its own algebra."
            ),
            options=(
                ("gaussian", "Gaussian", "Widths subtract in quadrature."),
                ("lorentzian", "Lorentzian", "Widths subtract directly."),
                ("pseudo_voigt", "Pseudo-Voigt", "Components separated, then recombined."),
            ),
            default="gaussian",
            group="Specimen",
        ),
        NumberParameter(
            name="shape_factor",
            label="Scherrer shape factor K",
            help_text=(
                "0.9 by convention for roughly spherical crystallites measured by FWHM. Every "
                "reported Scherrer size assumes a value of K, and most do not say which."
            ),
            default=0.9,
            minimum=0.5,
            maximum=1.5,
            advanced=True,
            group="Specimen",
        ),
    ),
    returns=(
        "One row per specimen reflection with its measured, instrumental and sample-only "
        "widths, plus the fitted crystallite size and microstrain."
    ),
    panel="xrd",
    citations=(_CITATION_WH, _CITATION_CAGLIOTI_UVW),
    tags=("XRD", "crystallite size", "microstrain", "Williamson-Hall", "Scherrer", "broadening"),
)
def _size_strain(request: dict[str, Any]) -> dict[str, Any]:
    def _peaks(text: str, field: str) -> tuple[np.ndarray, np.ndarray]:
        angles: list[float] = []
        widths: list[float] = []
        for number, line in enumerate(str(text).splitlines(), start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.replace(",", " ").split()
            if len(fields) < 2:
                raise InvalidInputError(
                    f"Line {number} needs an angle and a width.",
                    field=field,
                    hint="Two numbers per line: 2\u03b8 in degrees, then FWHM in degrees.",
                )
            try:
                angles.append(float(fields[0]))
                widths.append(float(fields[1]))
            except ValueError as error:
                raise InvalidInputError(
                    f"Line {number} is not two numbers: {stripped!r}.",
                    field=field,
                    hint="Remove any header row, or comment it out with a leading `#`.",
                ) from error
        return np.asarray(angles, dtype=float), np.asarray(widths, dtype=float)

    standard_angles, standard_widths = _peaks(request["standard_peaks"], "standard_peaks")
    sample_angles, sample_widths = _peaks(request["sample_peaks"], "sample_peaks")
    wavelength = float(request["wavelength_angstrom"])
    shape_factor = float(request["shape_factor"])
    mode = cast(Literal["gaussian", "lorentzian", "pseudo_voigt"], request["mode"])

    try:
        instrument = calibrate_instrument_broadening(
            standard_angles, standard_widths, name="calibrated from the pasted standard"
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The instrumental resolution function could not be fitted: {error}",
            field="standard_peaks",
            hint="Give at least three standard peaks spanning a wide angular range.",
        ) from error
    try:
        sample_only = deconvolve_instrument_width(
            sample_widths, instrument, sample_angles, mode=mode
        )
        analysis = williamson_hall(
            sample_angles,
            sample_only,
            wavelength_angstrom=wavelength,
            shape_factor=shape_factor,
        )
    except ValueError as error:
        raise InvalidInputError(
            f"Size and strain could not be separated: {error}",
            field="sample_peaks",
            hint="Specimen peaks must be wider than the instrument's own width at the same "
            "angle, and must span enough range to constrain the intercept.",
        ) from error

    instrument_widths = np.asarray(instrument.fwhm_deg(analysis.two_theta_deg), dtype=float)
    scherrer = scherrer_size_nm(
        analysis.sample_fwhm_deg,
        analysis.two_theta_deg,
        wavelength_angstrom=wavelength,
        shape_factor=shape_factor,
    )
    observed_sorted = np.asarray(
        [
            float(sample_widths[int(np.argmin(np.abs(sample_angles - angle)))])
            for angle in analysis.two_theta_deg
        ]
    )
    rows = tuple(
        {
            "two_theta_deg": float(analysis.two_theta_deg[index]),
            "observed_fwhm_deg": float(observed_sorted[index]),
            "instrument_fwhm_deg": float(instrument_widths[index]),
            "sample_fwhm_deg": float(analysis.sample_fwhm_deg[index]),
            "scherrer_size_nm": float(scherrer[index]),
            "abscissa": float(analysis.abscissa[index]),
            "ordinate": float(analysis.ordinate[index]),
        }
        for index in range(analysis.reflection_count)
    )
    notes = [
        "The instrumental width is a property of the diffractometer and its slits. A "
        "calibration only applies to widths measured in the same configuration and read the "
        "same way as the standard's.",
        "Scherrer sizes are quoted per reflection because their disagreement is informative: "
        "when they fall systematically with angle, the extra broadening is strain and the "
        "Williamson-Hall line is the answer to read.",
    ]
    if analysis.r_squared < 0.8:
        notes.insert(
            0,
            "The straight-line fit is poor, so the uniform-deformation assumption is not "
            "supported by these reflections and the size and strain should not be quoted.",
        )
    if analysis.microstrain < 0.0:
        notes.insert(
            0,
            "The fitted strain is negative, which the model cannot mean physically. Usually the "
            "instrumental width has been over-subtracted, or size broadening dominates and the "
            "slope is fitting noise.",
        )
    result = AppResult(
        title="Crystallite size and microstrain",
        summary=(
            f"D = {analysis.crystallite_size_nm:.4g} nm and \u03b5 = "
            f"{analysis.microstrain:.4g} ({100.0 * analysis.microstrain:.3g}%) from "
            f"{analysis.reflection_count} reflections, R\u00b2 = {analysis.r_squared:.5f}. "
            f"The instrument contributes FWHM\u00b2 = {instrument.caglioti_u:.6g} "
            f"tan\u00b2\u03b8 + {instrument.caglioti_v:.6g} tan\u03b8 + "
            f"{instrument.caglioti_w:.6g} deg\u00b2."
        ),
        table=ResultTable(
            columns=_SIZE_STRAIN_COLUMNS,
            rows=rows,
            caption="Specimen reflections, their widths at each stage, and the Williamson-Hall "
            "coordinates.",
        ),
        data={
            "abscissa": analysis.abscissa.tolist(),
            "ordinate": analysis.ordinate.tolist(),
            "slope": analysis.slope,
            "intercept": analysis.intercept,
            "r_squared": analysis.r_squared,
            "crystallite_size_nm": analysis.crystallite_size_nm,
            "microstrain": analysis.microstrain,
            "shape_factor": shape_factor,
            "wavelength_angstrom": wavelength,
            "caglioti_u": instrument.caglioti_u,
            "caglioti_v": instrument.caglioti_v,
            "caglioti_w": instrument.caglioti_w,
            "columns": [column.to_json() for column in _SIZE_STRAIN_COLUMNS],
            "describe": analysis.describe(),
            "instrument_describe": instrument.describe(),
        },
        inputs={
            "standard_peaks": request["standard_peaks"],
            "sample_peaks": request["sample_peaks"],
            "wavelength_angstrom": wavelength,
            "mode": mode,
            "shape_factor": shape_factor,
        },
        notes=tuple(notes),
        citations=(_CITATION_WH, _CITATION_CAGLIOTI_UVW),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="xrd.example.nickel_doublet",
            title="Nickel standard with a resolved Cu doublet",
            panel="xrd",
            summary="fcc nickel from 40–145° 2θ with narrow pseudo-Voigt peaks.",
            teaches=(
                "The fcc 111/200/220/311 sequence is the indexing anchor; at high angle the "
                "weaker Kα2 partner separates visibly to the right of Kα1."
            ),
            operation="xrd.powder_pattern",
            request={
                "phase": {"builtin": "ni_fcc"},
                "radiation": "cu_ka_doublet",
                "two_theta_min_deg": 40.0,
                "two_theta_max_deg": 145.0,
                "profile": "pseudo_voigt",
                "fwhm_deg": 0.08,
                "pseudo_voigt_eta": 0.35,
                "resolution_deg": 0.01,
            },
        ),
        ExampleScenario(
            id="xrd.example.silicon",
            title="Silicon phase-identification fingerprint",
            panel="xrd",
            summary="Diamond-cubic silicon with Cu Kα and tabulated form factors.",
            teaches=(
                "Diamond glide extinctions remove many geometrically possible lines. The remaining "
                "peak positions and intensity pattern form the familiar silicon fingerprint."
            ),
            operation="xrd.powder_pattern",
            request={
                "phase": {"builtin": "si_diamond"},
                "radiation": "cu_ka",
                "two_theta_min_deg": 20.0,
                "two_theta_max_deg": 120.0,
                "profile": "gaussian",
                "fwhm_deg": 0.12,
            },
        ),
        ExampleScenario(
            id="xrd.example.molybdenum_nickel",
            title="What a shorter wavelength changes",
            panel="xrd",
            summary="The nickel standard under Mo Kα1/Kα2 radiation.",
            teaches=(
                "The crystal is unchanged, but the shorter Mo wavelength moves every family to "
                "lower 2θ and places more reciprocal-lattice points inside a fixed scan window."
            ),
            operation="xrd.powder_pattern",
            request={
                "phase": {"builtin": "ni_fcc"},
                "radiation": "mo_ka_doublet",
                "two_theta_min_deg": 15.0,
                "two_theta_max_deg": 80.0,
                "fwhm_deg": 0.1,
            },
        ),
        ExampleScenario(
            id="xrd.example.background_snip",
            title="What the background is doing",
            panel="xrd",
            summary="A generated nickel scan with a curved low-angle background, clipped by SNIP.",
            teaches=(
                "The background here is not flat — it carries a broad low-angle hump of the "
                "kind an amorphous holder or a fluorescing specimen produces. SNIP follows it "
                "because it assumes only that the background varies more slowly than the peaks; "
                "drop the clipping window towards a peak width and watch it start eating the "
                "peaks instead."
            ),
            operation="xrd.background",
            request={
                "phase": {"builtin": "ni_fcc"},
                "data_source": "demonstration",
                "method": "snip",
                "half_window_deg": 2.0,
                "radiation": "cu_ka",
            },
        ),
        ExampleScenario(
            id="xrd.example.lattice_average_fails",
            title="Averaging a lattice parameter, and why it fails",
            panel="xrd",
            summary=(
                "The naive per-reflection mean on a scan with a known detector zero error."
            ),
            teaches=(
                "The demonstration scan carries a cell dilated by 1.003 and a 0.05\u00b0 "
                "detector zero error, so the answer is known before the calculation starts. "
                "Averaging a lattice parameter over the reflections lands about 3000 "
                "microangstrom low \u2014 8 parts in 10\u2074, when elastic strains of "
                "engineering interest are 1 part in 10\u2074. Nothing about the arithmetic is "
                "wrong: \u0394d/d = \u2212cot\u03b8\u00b7\u0394\u03b8 makes a fixed "
                "angular error a \u03b8-dependent spacing error, and averaging a bias does not "
                "remove it. Now run the companion example."
            ),
            operation="xrd.lattice_parameters",
            request={
                "phase": {"builtin": "ni_fcc"},
                "radiation": "cu_ka_doublet",
                "method": "average",
                "extrapolation": "none",
                "specimen_displacement_mm": 0.0,
            },
        ),
        ExampleScenario(
            id="xrd.example.lattice_cohen_extrapolation",
            title="The same scan, extrapolated to \u03b8 = 90\u00b0",
            panel="xrd",
            summary=(
                "Cohen least squares with the extrapolation function that matches the "
                "aberration."
            ),
            teaches=(
                "Identical data, one changed assumption. A detector zero error is a constant "
                "\u0394(2\u03b8), and since \u0394(sin\u00b2\u03b8) = "
                "sin\u03b8\u00b7cos\u03b8\u00b7\u0394(2\u03b8) = "
                "sin\u00b2\u03b8\u00b7cot\u03b8\u00b7\u0394(2\u03b8), the matching "
                "extrapolation function is cot\u03b8 \u2014 which vanishes at \u03b8 = "
                "90\u00b0, so the fitted cell is the extrapolated one. The determined a now "
                "lands within about 10 microangstrom of the truth, 3 parts in 10\u2076: better "
                "than the average by more than two orders of magnitude.\n\n"
                "In the plot, the scatter of the points about the line is the random error and "
                "the *slope* is the systematic one. Averaging the points would land on their "
                "mean; the answer is the intercept. Try Nelson\u2013Riley and "
                "cos\u00b2\u03b8/sin\u03b8 as well: they are the wrong shape for a zero "
                "error and leave most of it behind, which is the point of choosing the "
                "function to match the aberration."
            ),
            operation="xrd.lattice_parameters",
            request={
                "phase": {"builtin": "ni_fcc"},
                "radiation": "cu_ka_doublet",
                "method": "cohen",
                "extrapolation": "cot_theta",
                "specimen_displacement_mm": 0.0,
            },
        ),
        ExampleScenario(
            id="xrd.example.lattice_hexagonal_le_bail",
            title="A hexagonal cell from the whole pattern",
            panel="xrd",
            summary="Titanium a and c together, by Le Bail decomposition.",
            teaches=(
                "Outside the cubic system a lattice parameter *per reflection* does not exist: "
                "one reflection cannot determine both a and c, so the average method refuses "
                "the phase rather than returning a number. The joint solution is the only kind "
                "available, and hexagonal patterns overlap badly enough that fitting individual "
                "peaks runs out of resolvable lines before the reflection list runs out of "
                "reflections.\n\n"
                "Le Bail decomposition uses every measured point and extracts the reflection "
                "intensities rather than modelling them, so neither texture nor a wrong atomic "
                "basis can bias the cell. Watch c/a, which is the quantity hexagonal work turns "
                "on, and switch the systematic term to 'Neither' to see the goodness of fit "
                "fail loudly rather than quietly."
            ),
            operation="xrd.lattice_parameters",
            request={
                "phase": {"builtin": "ti_hcp"},
                "radiation": "cu_ka_doublet",
                "method": "le_bail",
                "systematic": "zero",
                "specimen_displacement_mm": 0.0,
            },
        ),
        ExampleScenario(
            id="xrd.example.rietveld_zero_and_cell",
            title="A cell and a zero error, told apart",
            panel="xrd",
            summary="Refining a generated scan whose cell and detector zero are both wrong.",
            teaches=(
                "The scan was built with the cell dilated by 1.003 and the detector zero out "
                "by 0.05 degrees, and the refinement recovers both. Turn off 'Refine the "
                "detector zero' and run again: the cell takes a visibly wrong value, because it "
                "is now the only parameter left that can move the peaks.\n\n"
                "Notice that the fit degrades too, rather than the error being absorbed "
                "silently. A zero error displaces every peak by the same angle while a cell "
                "dilation displaces them in proportion to tan(theta), so over a wide scan the "
                "two cannot fully trade places. That is exactly why refining both works — "
                "and why the same experiment over a narrow angular range would be far less "
                "forgiving, with the two parameters nearly indistinguishable and the "
                "uncertainties saying so."
            ),
            operation="xrd.rietveld",
            request={
                "phase": {"builtin": "ni_fcc"},
                "data_source": "demonstration",
                "radiation": "cu_ka",
                "refine_cell": True,
                "refine_zero_shift": True,
                "refine_width": True,
                "background_degree": 8,
            },
        ),
        ExampleScenario(
            id="xrd.example.rietveld_width_wrong",
            title="A good R factor hiding a bad model",
            panel="xrd",
            summary="The same scan with the peak width held at the wrong value.",
            teaches=(
                "Every index degrades, but not equally. What moves furthest are the two that "
                "are sensitive to the *shape* of the misfit rather than its size — R_Bragg "
                "and the Durbin–Watson statistic, which collapses. The difference curve "
                "acquires the derivative shape of the peak, at every reflection, which is why a "
                "refinement is read from that curve and not from one number."
            ),
            operation="xrd.rietveld",
            request={
                "phase": {"builtin": "ni_fcc"},
                "data_source": "demonstration",
                "radiation": "cu_ka",
                "refine_cell": True,
                "refine_zero_shift": True,
                "refine_width": False,
                "starting_fwhm_deg": 0.1,
                "background_degree": 8,
            },
        ),
        ExampleScenario(
            id="xrd.example.size_strain",
            title="Twenty-five nanometres and two parts per thousand",
            panel="xrd",
            summary="LaB6 standard widths calibrate the instrument; the specimen gives up its "
            "size and strain.",
            teaches=(
                "The specimen widths were built from a 25 nm crystallite size and a 0.2% "
                "microstrain, and the analysis returns exactly that — so the method can be "
                "checked rather than trusted. Note the Scherrer column falling steadily with "
                "angle: that fall *is* the strain, expressed as the disagreement between "
                "estimates that assume it away."
            ),
            operation="xrd.size_strain",
            request={"mode": "gaussian", "wavelength_angstrom": 1.5406, "shape_factor": 0.9},
        ),
        ExampleScenario(
            id="xrd.example.zirconium",
            title="Alpha-zirconium hexagonal powder pattern",
            panel="xrd",
            summary="hcp zirconium with Cu Kα over a conventional laboratory range.",
            teaches=(
                "The basal and prismatic families separate because a and c are independent in "
                "the hexagonal metric; labels use canonical four-index Miller–Bravais notation."
            ),
            operation="xrd.powder_pattern",
            request={
                "phase": {"builtin": "zr_hcp"},
                "radiation": "cu_ka_doublet",
                "two_theta_min_deg": 25.0,
                "two_theta_max_deg": 120.0,
                "fwhm_deg": 0.16,
            },
        ),
    )
)


_CITATION_COHEN_LS = "Cohen, Rev. Sci. Instrum. 6 (1935) 68, doi:10.1063/1.1751937."
_CITATION_NELSON_RILEY_FN = (
    "Nelson & Riley, Proc. Phys. Soc. 57 (1945) 160, doi:10.1088/0959-5309/57/3/302."
)
_CITATION_LE_BAIL_METHOD = (
    "Le Bail, Duroy & Fourquet, Mater. Res. Bull. 23 (1988) 447, "
    "doi:10.1016/0025-5408(88)90019-0."
)
_CITATION_CULLITY_CH11 = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Ch. 11 "
    "(Precise Parameter Measurements)."
)

_LATTICE_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column(
        "two_theta_observed_deg",
        "2\u03b8 observed",
        units="\u00b0",
        numeric=True,
        digits=5,
        help_text="Fitted K\u03b11 position after profile fitting.",
    ),
    Column(
        "standard_uncertainty_mdeg",
        "\u03c3(2\u03b8)",
        units="m\u00b0",
        numeric=True,
        digits=2,
        help_text="Position uncertainty from the profile fit, in millidegrees.",
    ),
    Column(
        "two_theta_calculated_deg",
        "2\u03b8 from the cell",
        units="\u00b0",
        numeric=True,
        digits=5,
        help_text="Where the determined cell puts this reflection.",
    ),
    Column(
        "residual_mdeg",
        "Residual",
        units="m\u00b0",
        numeric=True,
        digits=2,
        help_text="Observed minus calculated. Structure here is structure the model missed.",
    ),
    Column(
        "systematic_shift_mdeg",
        "Removed by the drift term",
        units="m\u00b0",
        numeric=True,
        digits=2,
        help_text=(
            "How far the refined systematic-error term moved this reflection. Compare it with "
            "\u03c3(2\u03b8): if it is much larger, the correction did real work."
        ),
    ),
    Column(
        "d_observed_angstrom",
        "d observed",
        units="\u00c5",
        numeric=True,
        digits=6,
    ),
)


@REGISTRY.operation(
    "xrd.lattice_parameters",
    title="Determine lattice parameters",
    summary=(
        "Precise cell determination with the systematic error refined away, not averaged over."
    ),
    help_text=(
        "This determines a unit cell. It does not refine a structure: the atomic basis is held "
        "fixed and only the cell and the errors of the instrument are varied, which is what "
        "stops texture and an imperfect structural model from leaking into the answer.\n\n"
        "Averaging a lattice parameter over several reflections is the intuitive method and it "
        "does not work. Differentiating Bragg's law gives \u0394d/d = \u2212cot\u03b8 "
        "\u0394\u03b8, so a fixed angular error produces a *\u03b8-dependent* error in the "
        "spacing. The errors that dominate a laboratory scan \u2014 a detector zero offset, a "
        "specimen a few tens of micrometres off the diffractometer axis \u2014 are systematic, "
        "so averaging divides the random scatter by \u221aN and leaves the bias untouched. Run "
        "the average method against the extrapolation methods on the same demonstration scan "
        "and read the difference: it is roughly three orders of magnitude.\n\n"
        "Cohen's method is the default. Because sin\u00b2\u03b8 = (\u03bb\u00b2/4)\u00b7"
        "h\u1d40G*h is *linear* in the reciprocal metric tensor, one solution covers cubic "
        "through triclinic with no starting guess and an analytic covariance \u2014 so the "
        "uncertainties quoted here are real, not decorative. A systematic-error coefficient is "
        "refined alongside the cell against the chosen extrapolation function, every one of "
        "which vanishes at \u03b8 = 90\u00b0. That is why extrapolation works, and why the "
        "highest-angle reflections carry nearly all the weight.\n\n"
        "Le Bail decomposition uses every measured point instead of a handful of fitted "
        "positions, and is the method for a hexagonal or lower-symmetry pattern where the "
        "reflections overlap. Its intensities are extracted rather than modelled, so neither "
        "texture nor a wrong basis can bias the cell.\n\n"
        "For residual stress this supplies the spacings, not the stress. A symmetric "
        "\u03b8\u20132\u03b8 scan measures planes parallel to the surface only; a stress "
        "needs several specimen tilts and the X-ray elastic constants of the reflection used."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "The phase whose cell is determined. Its symmetry decides how many cell "
                "parameters may vary \u2014 one for cubic, two for hexagonal \u2014 and its "
                "reflections are what the measured peaks are indexed against."
            ),
            builtin="ni_fcc",
        ),
        *_scan_parameters(),
        NumberParameter(
            name="specimen_displacement_mm",
            label="Injected specimen displacement",
            help_text=(
                "Added to the demonstration scan as a *known* aberration, so the methods can be "
                "judged against a value that is known in advance. A real specimen is routinely "
                "50 \u00b5m off the axis. Ignored for a pasted scan."
            ),
            units="mm",
            default=0.0,
            minimum=-1.0,
            maximum=1.0,
            group="Measurement",
        ),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text=(
                "The wavelength every spacing is referred to. A wrong choice scales every "
                "lattice parameter by the wavelength ratio."
            ),
            options=(
                ("cu_ka", "Cu K\u03b1 (single averaged line)", "One copper line."),
                ("cu_ka_doublet", "Cu K\u03b11/K\u03b12", "Common laboratory copper doublet."),
                ("co_ka_doublet", "Co K\u03b11/K\u03b12", "Reduces Fe fluorescence."),
                ("mo_ka_doublet", "Mo K\u03b11/K\u03b12", "Short-wavelength molybdenum."),
            ),
            default="cu_ka_doublet",
        ),
        ChoiceParameter(
            name="method",
            label="Method",
            help_text=(
                "The naive average is offered so its failure can be seen rather than described. "
                "Cohen is the default. Le Bail is for overlapped patterns."
            ),
            options=(
                (
                    "cohen",
                    "Cohen least squares (recommended)",
                    "Linear in the metric tensor; every crystal system; real uncertainties.",
                ),
                (
                    "average",
                    "Average over reflections (teaching comparison)",
                    "Cubic only. Cannot remove a systematic error; run it to see by how much.",
                ),
                (
                    "le_bail",
                    "Le Bail whole-pattern decomposition",
                    "Every measured point; handles overlapped reflections.",
                ),
            ),
            default="cohen",
            group="Method",
        ),
        ChoiceParameter(
            name="extrapolation",
            label="Extrapolation function",
            help_text=(
                "The angular form the systematic error is assumed to take. cos\u00b2\u03b8/"
                "sin\u03b8 is the exact form for specimen displacement on a Bragg\u2013Brentano "
                "instrument; Nelson\u2013Riley approximates displacement and absorption "
                "together and is the usual choice. All of them vanish at \u03b8 = 90\u00b0, "
                "which is what makes the fitted cell the extrapolated one. Ignored by the "
                "average and Le Bail methods."
            ),
            options=(
                (
                    "nelson_riley",
                    "Nelson\u2013Riley",
                    "\u00bd(cos\u00b2\u03b8/sin\u03b8 + cos\u00b2\u03b8/\u03b8); "
                    "the standard choice.",
                ),
                (
                    "cos_squared_over_sin",
                    "cos\u00b2\u03b8/sin\u03b8",
                    "Exact for specimen displacement.",
                ),
                (
                    "cot_theta",
                    "cot\u03b8",
                    "Exact for a detector zero error; the fitted D is the zero itself.",
                ),
                (
                    "bradley_jay",
                    "Bradley\u2013Jay (cos\u00b2\u03b8)",
                    "Gives Cohen's classical sin\u00b2(2\u03b8) drift column.",
                ),
                (
                    "none",
                    "None",
                    "No correction. Run it to see what the correction was worth.",
                ),
            ),
            default="nelson_riley",
            group="Method",
        ),
        ChoiceParameter(
            name="systematic",
            label="Le Bail systematic term",
            help_text=(
                "Which aberration the whole-pattern fit refines. Exactly one: a zero and a "
                "displacement differ only as constant against cos\u03b8, and over one scan's "
                "angular range that difference is comparable to the noise, so refining both is "
                "ill-conditioned. Zero belongs to a calibrated instrument; displacement belongs "
                "to the specimen."
            ),
            options=(
                ("displacement", "Specimen displacement", "Refined in millimetres."),
                ("zero", "Detector zero", "Refined in degrees 2\u03b8."),
                ("none", "Neither", "Run it to see the fit fail."),
            ),
            default="displacement",
            group="Method",
        ),
        NumberParameter(
            name="minimum_two_theta_deg",
            label="Discard reflections below",
            help_text=(
                "The crudest defence against systematic error: cot\u03b8 shrinks towards "
                "back-reflection, so high-angle reflections are intrinsically more precise. It "
                "helps, and it is no substitute for refining the drift term."
            ),
            units="\u00b0 2\u03b8",
            default=0.0,
            minimum=0.0,
            maximum=175.0,
            group="Method",
        ),
        NumberParameter(
            name="expected_fwhm_deg",
            label="Expected peak width",
            help_text=(
                "Sets the scale of the matched filter that finds the peaks and the size of each "
                "fit window. It only needs to be right to within about a factor of two."
            ),
            units="\u00b0 2\u03b8",
            default=0.14,
            minimum=0.01,
            maximum=2.0,
            advanced=True,
            group="Peak finding",
        ),
        NumberParameter(
            name="prominence_sigma",
            label="Detection threshold",
            help_text=(
                "How far above the noise a feature must rise to be treated as a peak, in robust "
                "standard deviations of the matched-filter response."
            ),
            default=5.0,
            minimum=1.0,
            maximum=50.0,
            advanced=True,
            group="Peak finding",
        ),
        NumberParameter(
            name="tolerance_deg",
            label="Indexing tolerance",
            help_text=(
                "How far a measured peak may sit from a calculated reflection and still be "
                "matched to it. Wider than any uncorrected zero or displacement error, and "
                "narrower than the spacing between neighbouring lines."
            ),
            units="\u00b0 2\u03b8",
            default=0.3,
            minimum=0.01,
            maximum=3.0,
            advanced=True,
            group="Peak finding",
        ),
        IntegerParameter(
            name="max_index",
            label="Maximum Miller index",
            help_text="Largest |h|, |k|, |l| enumerated when predicting reflections.",
            default=6,
            minimum=1,
            maximum=12,
            advanced=True,
            group="Peak finding",
        ),
    ),
    returns=(
        "The determined cell with standard uncertainties, the refined systematic-error term, "
        "and the per-reflection residuals behind both."
    ),
    panel="xrd",
    citations=(
        _CITATION_CULLITY_CH11,
        _CITATION_COHEN_LS,
        _CITATION_NELSON_RILEY_FN,
        _CITATION_LE_BAIL_METHOD,
    ),
    tags=(
        "XRD",
        "lattice parameter",
        "precise",
        "Cohen",
        "Nelson-Riley",
        "Le Bail",
        "strain",
        "stress",
    ),
)
def _lattice_parameters(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    radiation = _RADIATION[str(request["radiation"])]()
    measured, generated = _measured_from_request(request, phase, radiation)

    displacement = float(request["specimen_displacement_mm"])
    if generated and displacement != 0.0:
        axis = np.asarray(measured.two_theta_deg, dtype=float)
        measured = MeasuredPowderPattern(
            name=measured.name,
            two_theta_deg=axis
            + specimen_displacement_shift_deg(
                axis, displacement_mm=displacement, goniometer_radius_mm=240.0
            ),
            intensity=measured.intensity,
            radiation=radiation,
            synthetic=True,
        )

    method = cast(Literal["cohen", "average", "le_bail"], request["method"])
    extrapolation = cast(Any, request["extrapolation"])
    instrument = InstrumentBroadening.ideal(float(request["expected_fwhm_deg"]))

    try:
        system = crystal_system_of(phase)
    except ValueError as error:
        raise InvalidInputError(
            f"This phase has no cell parameterization here: {error}",
            field="phase",
            hint="Choose a phase in one of the seven crystal systems.",
        ) from error

    floor = float(request["minimum_two_theta_deg"])
    try:
        result, indexing = determine_lattice_parameters_from_pattern(
            measured,
            phase,
            method=method,
            extrapolation=extrapolation,
            radiation=radiation,
            instrument=instrument,
            systematic=cast(Any, request["systematic"]),
            tolerance_deg=float(request["tolerance_deg"]),
            max_index=int(request["max_index"]),
            prominence_sigma=float(request["prominence_sigma"]),
            minimum_two_theta_deg=floor if floor > 0.0 else None,
            phase_name=spec.name,
        )
    except ValueError as error:
        message = str(error)
        if "cubic cell" in message:
            field, hint = (
                "method",
                "The average method needs a cubic phase; every other cell needs a joint "
                "solution, so choose Cohen least squares or Le Bail.",
            )
        elif "detected" in message or "fitted" in message:
            field, hint = (
                "prominence_sigma",
                "Lower the detection threshold, or check the expected peak width.",
            )
        elif "angular restriction" in message:
            field, hint = (
                "minimum_two_theta_deg",
                "Lower the angular floor so more reflections survive.",
            )
        else:
            field, hint = (
                "tolerance_deg",
                "Widen the indexing tolerance, or check that the phase and radiation are right.",
            )
        raise InvalidInputError(
            f"The lattice parameters could not be determined: {message}",
            field=field,
            hint=hint,
        ) from error

    shifts = result.systematic_shift_deg
    calculated = result.two_theta_deg - result.residual_two_theta_deg
    # Keyed on the angle, not on position. An angular floor filters the
    # determination's reflections without filtering the indexing's, so a
    # positional lookup would attach each uncertainty to the wrong reflection
    # exactly when the operator restricts the range -- which is the case they
    # would restrict it for.
    sigma_by_angle = (
        {}
        if indexing is None
        else {
            round(item.peak.two_theta_deg, 9): (
                1000.0 * item.peak.two_theta_standard_uncertainty_deg
            )
            for item in indexing.reflections
        }
    )
    wavelength = float(radiation.wavelength_angstrom)
    rows = tuple(
        {
            "hkl_label": _powder_label(indices, spec=spec),
            "two_theta_observed_deg": float(angle),
            "standard_uncertainty_mdeg": float(
                sigma_by_angle.get(round(float(angle), 9), 0.0)
            ),
            "two_theta_calculated_deg": float(calculated[index]),
            "residual_mdeg": 1000.0 * float(result.residual_two_theta_deg[index]),
            "systematic_shift_mdeg": 1000.0 * float(shifts[index]),
            "d_observed_angstrom": float(
                wavelength / (2.0 * np.sin(np.deg2rad(0.5 * float(angle))))
            ),
        }
        for index, (indices, angle) in enumerate(
            zip(result.miller_indices, result.two_theta_deg, strict=True)
        )
    )

    # The plot each method deserves. The whole-pattern branch is tested first
    # and not last: a Le Bail fit measures no individual peak position, so a
    # per-reflection lattice parameter computed from its output would be
    # computed from the *calculated* angles -- a plot of the model against
    # itself, which would look convincing and mean nothing. Its diagnostic is
    # the difference curve. A cubic cell determined from fitted positions gets
    # the classical extrapolation, where a lattice parameter per reflection
    # does exist and the intercept is the answer. Everything else gets the
    # residual against angle, which works in every crystal system.
    if result.profile_two_theta_deg is not None:
        assert result.profile_observed is not None
        assert result.profile_calculated is not None
        keep = _decimate(int(result.profile_two_theta_deg.size))
        plot = {
            "plot_kind": "profile",
            "abscissa": result.profile_two_theta_deg[keep].tolist(),
            "ordinate": result.profile_observed[keep].tolist(),
            "calculated": result.profile_calculated[keep].tolist(),
            "difference": (
                result.profile_observed[keep] - result.profile_calculated[keep]
            ).tolist(),
            "abscissa_label": "2\u03b8 (\u00b0)",
            "ordinate_label": "intensity, background removed",
            "line_slope": 0.0,
            "line_intercept": 0.0,
            "determined": float(result.a),
        }
    elif system == "cubic" and result.two_theta_deg.size >= 2:
        sums = np.array(
            [sum(value**2 for value in indices) for indices in result.miller_indices],
            dtype=float,
        )
        spacings = wavelength / (2.0 * np.sin(np.deg2rad(0.5 * result.two_theta_deg)))
        per_reflection = spacings * np.sqrt(sums)
        function = (
            "nelson_riley" if result.extrapolation == "none" else result.extrapolation
        )
        abscissa = extrapolation_values(result.two_theta_deg, function=cast(Any, function))
        slope, intercept = np.polyfit(abscissa, per_reflection, 1)
        plot = {
            "plot_kind": "extrapolation",
            "abscissa": abscissa.tolist(),
            "ordinate": per_reflection.tolist(),
            "abscissa_label": _EXTRAPOLATION_LABELS[function],
            "ordinate_label": "a from each reflection (\u00c5)",
            "line_slope": float(slope),
            "line_intercept": float(intercept),
            "determined": float(result.a),
        }
    else:
        plot = {
            "plot_kind": "residual",
            "abscissa": result.two_theta_deg.tolist(),
            "ordinate": (1000.0 * result.residual_two_theta_deg).tolist(),
            "abscissa_label": "2\u03b8 (\u00b0)",
            "ordinate_label": "observed \u2212 calculated (m\u00b0)",
            "line_slope": 0.0,
            "line_intercept": 0.0,
            "determined": float(result.a),
        }

    notes: list[str] = []
    if generated:
        notes.extend(_demonstration_notes(phase.lattice.a))
        if displacement != 0.0:
            lowest = float(result.two_theta_deg[0])
            worst = 1000.0 * abs(
                float(
                    specimen_displacement_shift_deg(
                        [lowest], displacement_mm=displacement, goniometer_radius_mm=240.0
                    )[0]
                )
            )
            notes.append(
                f"A {1000.0 * displacement:.0f} \u00b5m specimen displacement was injected on "
                "top of that, at a 240 mm goniometer radius. It moves the lowest-angle "
                f"reflection by {worst:.0f} m\u00b0 and the highest by less, which is exactly "
                "the \u03b8-dependence that averaging cannot remove. Note that the scan now "
                "carries two aberrations of different angular form \u2014 a constant zero and "
                "a cos\u03b8 displacement \u2014 and no single extrapolation function removes "
                "both. That is why precise work calibrates the zero against a standard first."
            )
        notes.append(
            "The cell to recover is a = "
            f"{phase.lattice.a * _DEMO_LATTICE_SCALE:.5f} \u00c5. Compare it with the value "
            "above, then switch the method and watch it move."
        )
    notes.append(
        "Read the drift column against \u03c3(2\u03b8). A correction much larger than the "
        "position uncertainties did real work; one much smaller means the specimen was already "
        "well aligned and every method should agree."
    )
    notes.append(
        "This is a lattice parameter, not a stress. A symmetric scan measures the spacing of "
        "planes parallel to the specimen surface only; a stress needs several specimen tilts "
        "and the X-ray elastic constants of the reflection used."
    )

    if method == "le_bail":
        # A whole-pattern decomposition never fits an individual peak position,
        # so it has no observed angle, no position uncertainty and no residual
        # to report per reflection. Printing zeros in those columns would be a
        # claim; the narrower column set is the honest one.
        columns: tuple[Column, ...] = _LE_BAIL_COLUMNS
        rows = tuple(
            {
                key: row[key]
                for key in ("hkl_label", "two_theta_calculated_deg", "d_observed_angstrom")
            }
            for row in rows
        )
        caption = (
            f"{result.reflection_count} reflections modelled by the whole-pattern fit, at the "
            "angles the determined cell puts them. A Le Bail fit measures no individual peak "
            "position, so it reports none."
        )
    else:
        columns = _LATTICE_COLUMNS
        caption = (
            f"{result.reflection_count} reflections behind the determined cell, with the "
            "residual and the systematic shift removed from each."
        )
    result_payload = AppResult(
        title=f"Lattice parameters of {spec.name}",
        summary=(
            f"a = {result.a:.6f} \u00b1 {result.a_standard_uncertainty:.6f} \u00c5"
            + (
                ""
                if system == "cubic"
                else f", c = {result.c:.6f} \u00b1 {result.c_standard_uncertainty:.6f} "
                f"\u00c5, c/a = {result.axial_ratio:.6f}"
            )
            + f" from {result.reflection_count} reflections of a {system} cell, a relative "
            f"uncertainty of {result.relative_uncertainty:.1e}."
        ),
        table=ResultTable(columns=columns, rows=rows, caption=caption),
        data={
            **plot,
            "columns": [column.to_json() for column in columns],
            "method": method,
            "crystal_system": system,
            "a": float(result.a),
            "b": float(result.b),
            "c": float(result.c),
            "a_standard_uncertainty": float(result.a_standard_uncertainty),
            "c_standard_uncertainty": float(result.c_standard_uncertainty),
            "axial_ratio": float(result.axial_ratio),
            "relative_uncertainty": float(result.relative_uncertainty),
            "drift_coefficient": float(result.drift_coefficient),
            "drift_standard_uncertainty": float(result.drift_standard_uncertainty),
            "extrapolation": result.extrapolation,
            "reduced_chi_squared": float(result.reduced_chi_squared),
            "weighted_profile_r": result.weighted_profile_r,
            "reflection_count": int(result.reflection_count),
            "strain_relative_to_reference": result.strain_relative_to_reference,
            "figure_of_merit_m": (
                None if indexing is None else float(indexing.figure_of_merit_m()[0])
            ),
            "synthetic": generated,
            "describe": result.describe(),
            "phase_name": spec.name,
        },
        inputs={
            "phase": spec.to_json(),
            "data_source": request["data_source"],
            "specimen_displacement_mm": displacement,
            "radiation": request["radiation"],
            "method": method,
            "extrapolation": request["extrapolation"],
            "systematic": request["systematic"],
            "minimum_two_theta_deg": float(request["minimum_two_theta_deg"]),
            "expected_fwhm_deg": float(request["expected_fwhm_deg"]),
            "prominence_sigma": float(request["prominence_sigma"]),
            "tolerance_deg": float(request["tolerance_deg"]),
            "max_index": int(request["max_index"]),
            "demonstration_seed": int(request["demonstration_seed"]),
        },
        notes=tuple(notes),
        citations=(_CITATION_CULLITY_CH11, _CITATION_COHEN_LS, _CITATION_NELSON_RILEY_FN),
    )
    return result_payload.to_json()


_LE_BAIL_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column(
        "two_theta_calculated_deg",
        "2\u03b8 from the cell",
        units="\u00b0",
        numeric=True,
        digits=5,
        help_text="Where the determined cell puts this reflection.",
    ),
    Column("d_observed_angstrom", "d", units="\u00c5", numeric=True, digits=6),
)

_EXTRAPOLATION_LABELS = {
    "nelson_riley": "\u00bd(cos\u00b2\u03b8/sin\u03b8 + cos\u00b2\u03b8/\u03b8)",
    "cos_squared_over_sin": "cos\u00b2\u03b8 / sin\u03b8",
    "cot_theta": "cot\u03b8",
    "bradley_jay": "cos\u00b2\u03b8",
    "none": "cos\u00b2\u03b8 / sin\u03b8",
}


# ---------------------------------------------------------------------------
# Phase identification. The one operation on this panel that does not require
# the answer as an input: every other analysis is told the phase and measures
# something about it, while this one is told several and decides between them.
# ---------------------------------------------------------------------------

_CITATION_HANAWALT = (
    "Hanawalt, Rinn & Frevel, Ind. Eng. Chem. Anal. Ed. 10 (1938) 457, doi:10.1021/ac50125a001."
)
_CITATION_SMITH_SNYDER_FN = (
    "Smith & Snyder, J. Appl. Crystallogr. 12 (1979) 60, doi:10.1107/S002188987901178X."
)
_CITATION_DOLLASE_MARCH = (
    "Dollase, J. Appl. Crystallogr. 19 (1986) 267, doi:10.1107/S0021889886089458."
)
_CITATION_CULLITY_CH14 = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 14."
)

#: The three weightings the operation offers, and what each is for. Exposing
#: four raw numbers would put a scoring model on the control rail; exposing the
#: three specimen situations a laboratory actually meets puts the *decision*
#: there instead, which is the thing the operator knows and the software does
#: not.
_WEIGHTING_PRESETS: dict[str, dict[str, float]] = {
    "standard": {
        "explained_intensity_fraction": 0.40,
        "completeness": 0.25,
        "position_score": 0.20,
        "intensity_agreement": 0.15,
    },
    "textured": {
        "explained_intensity_fraction": 0.40,
        "completeness": 0.30,
        "position_score": 0.30,
        "intensity_agreement": 0.00,
    },
    "positions_only": {
        "explained_intensity_fraction": 0.50,
        "completeness": 0.00,
        "position_score": 0.50,
        "intensity_agreement": 0.00,
    },
}

_CANDIDATE_COLUMNS = (
    Column("rank", "Rank", numeric=True),
    Column("phase_name", "Candidate"),
    Column(
        "score",
        "Score",
        numeric=True,
        digits=3,
        help_text=(
            "Weighted mean of the four criteria, in [0, 1]. Read the criteria beside it: which "
            "one a candidate fails says more than its total does."
        ),
    ),
    Column(
        "explained",
        "Intensity explained",
        units="%",
        numeric=True,
        digits=1,
        help_text=(
            "Share of the measured integrated intensity carried by peaks this candidate "
            "indexed. A strong unindexed peak is the signature of a second phase."
        ),
    ),
    Column(
        "completeness",
        "Lines seen",
        units="%",
        numeric=True,
        digits=1,
        help_text=(
            "Share of the candidate's own strong reflections, inside the measured range, that "
            "were actually observed. This is what separates two cells differing by a centring: "
            "a centring is a statement about which lines are absent."
        ),
    ),
    Column(
        "position",
        "Position",
        numeric=True,
        digits=3,
        help_text=(
            "1 - mean|Δ2θ| / tolerance. How far inside the matching window the lines "
            "landed, not merely whether they landed inside it."
        ),
    ),
    Column(
        "intensity",
        "Intensity",
        numeric=True,
        digits=3,
        help_text=(
            "Bounded similarity of observed and calculated relative intensities. Weighted least "
            "of the four: preferred orientation moves intensities without moving positions."
        ),
    ),
    Column("indexed", "Peaks indexed"),
    Column(
        "cell_dilation_percent",
        "Cell",
        units="%",
        numeric=True,
        digits=3,
        help_text=(
            "How far the candidate's cell had to be dilated to place its lines, as a "
            "percentage. A few hundredths is the ordinary difference between a tabulated cell "
            "and a real solid solution; a value at the edge of the search range means the "
            "candidate had to be stretched to fit, and should be read with suspicion."
        ),
    ),
    Column(
        "figure_of_merit_m",
        "M",
        numeric=True,
        digits=1,
        help_text="de Wolff's figure of merit for this candidate's assignment.",
    ),
    Column("source", "From"),
)


def _candidate_phases(request: dict[str, Any]) -> tuple[list[tuple[str, Any]], dict[str, str]]:
    """Resolve the candidate list into named phases and their provenance.

    One malformed entry names itself in the error. A user who has opened five
    CIFs needs to be told *which* one the server could not read, and a message
    that only says "a phase could not be parsed" makes them close all five.
    """

    payload = request.get("candidates") or {}
    entries = payload.get("phases") if isinstance(payload, Mapping) else None
    if not isinstance(entries, list) or not entries:
        raise InvalidInputError(
            "No candidate phases were offered.",
            field="candidates",
            hint=(
                "Add at least two candidates — built-in phases or .cif files you open "
                "— so there is something to choose between. One candidate is a check on "
                "that phase rather than an identification."
            ),
        )

    named: list[tuple[str, Any]] = []
    sources: dict[str, str] = {}
    seen: dict[str, int] = {}
    for position, entry in enumerate(entries, start=1):
        if not isinstance(entry, Mapping):
            raise InvalidInputError(
                f"Candidate {position} is not a phase.",
                field="candidates",
                hint="Remove it and add the phase again.",
            )
        try:
            spec, phase = phase_from_request(entry.get("phase"))
        except InvalidInputError as error:
            raise InvalidInputError(
                f"Candidate {position} could not be read: {error.message}",
                field="candidates",
                hint=error.hint,
            ) from error
        label = str(entry.get("label") or spec.name or f"candidate {position}")
        count = seen.get(label, 0)
        seen[label] = count + 1
        if count:
            label = f"{label} ({count + 1})"
        named.append((label, phase))
        sources[label] = str(spec.source or "built-in catalogue")
    return named, sources


@REGISTRY.operation(
    "xrd.phase_identification",
    title="Identify the phase",
    summary="Rank several candidate structures against a measured scan, and say which fits.",
    help_text=(
        "Every other analysis on this panel is *told* the phase and measures something about "
        "it. This one is told several and decides between them, which is the step that comes "
        "first on a specimen whose identity is suspected rather than established.\n\n"
        "The peaks are detected and profile-fitted, then each candidate is indexed against them "
        "by a global one-to-one assignment and scored on four criteria. Read the criteria, not "
        "only the total — *which* one a candidate fails is the diagnosis:\n\n"
        "- **Intensity explained.** How much of the measured intensity the candidate accounts "
        "for. A strong peak it cannot explain means something else is in the specimen.\n"
        "- **Lines seen.** How many of the candidate's own strong reflections actually appeared. "
        "This is the criterion that separates a face-centred cell from a body-centred one, "
        "because a centring is a claim about which lines are *absent* rather than present.\n"
        "- **Position.** How far *inside* the matching tolerance the lines landed. A candidate "
        "sitting at the edge of the window throughout has the wrong cell dimensions even though "
        "every peak was formally indexed.\n"
        "- **Intensity.** Whether the relative intensities track the calculated ones. It carries "
        "the least weight on purpose: preferred orientation, microabsorption and a coarse powder "
        "move measured intensities by factors without moving a single peak position, so no "
        "candidate is rejected on intensity alone.\n\n"
        "A ranking always has a winner, which is not the same as having an answer. The verdict "
        "therefore states two things separately: whether the best candidate explains the pattern "
        "in absolute terms, and whether it beats the runner-up by enough to be distinguished "
        "from it. When it does not, the honest readings are printed — that none of the "
        "candidates offered accounts for this scan, or that this scan does not tell the top two "
        "apart and a longer count at high angle, a different wavelength or chemistry is needed.\n\n"
        "This is not a database search. The candidates must be supplied; the operation ranks "
        "what it is given and cannot propose a phase nobody thought of. Nor does it quantify a "
        "mixture: when several candidates each explain part of the pattern, the next step is a "
        "multi-phase Rietveld refinement, not a higher score."
    ),
    parameters=(
        ObjectParameter(
            name="candidates",
            label="Candidate phases",
            help_text=(
                "The structures to choose between. Add built-in phases or open .cif files; two "
                "or more make it an identification rather than a check. A candidate that cannot "
                "be indexed at all is scored zero with the reason stated rather than aborting "
                "the comparison, so one unreadable structure among five costs you only that one."
            ),
            editor="phase_candidates",
            default={
                "phases": [
                    {"phase": {"builtin": "ni_fcc"}},
                    {"phase": {"builtin": "fe_bcc"}},
                ]
            },
        ),
        phase_parameter(
            label="Demonstration specimen",
            help_text=(
                "Only used to generate the demonstration scan — it is the phase the "
                "synthetic specimen is *made of*, and the answer the ranking should recover. An "
                "experimental scan is analysed without reference to it, which is the point: an "
                "identification that consulted a declared phase would not be one."
            ),
            builtin="ni_fcc",
        ),
        *_scan_parameters(),
        ChoiceParameter(
            name="weighting",
            label="Evidence weighting",
            help_text=(
                "Which evidence to trust for *this* specimen. Balanced suits a well-prepared "
                "random powder. Choose the textured weighting for a rolled sheet, a coating or "
                "anything with a rolling or fibre texture, where measured intensities are moved "
                "by orientation rather than by structure. Positions only is the strictest "
                "setting: it asks solely where the lines are and how much intensity is left "
                "unexplained."
            ),
            options=(
                (
                    "standard",
                    "Balanced",
                    "All four criteria, with intensity weighted least.",
                ),
                (
                    "textured",
                    "Textured specimen",
                    "Ignore intensities: orientation moves them, not structure.",
                ),
                (
                    "positions_only",
                    "Positions only",
                    "Line positions and unexplained intensity alone.",
                ),
            ),
            default="standard",
            group="Scoring",
        ),
        NumberParameter(
            name="tolerance_deg",
            label="Matching tolerance",
            help_text=(
                "How far a calculated line may sit from a measured peak and still be matched. "
                "Set it wider than the instrument's uncorrected zero-point and displacement "
                "errors and narrower than the spacing between neighbouring calculated lines. It "
                "also sets the scale of the position criterion, so widening it to rescue a "
                "candidate judges every match against the laxer standard it was admitted under."
            ),
            units="° 2θ",
            default=0.3,
            minimum=0.01,
            maximum=3.0,
            group="Scoring",
            field_width="short",
        ),
        NumberParameter(
            name="prominence_sigma",
            label="Detection threshold",
            help_text=(
                "How far above the local noise a feature must rise to be fitted as a peak, in "
                "units of that noise. Lower it to admit weak lines, at the cost of admitting "
                "background structure with them. Peak detection is where an identification most "
                "often goes wrong, so the fitted peaks are listed with the result."
            ),
            units="σ",
            default=5.0,
            minimum=1.0,
            maximum=30.0,
            group="Scoring",
            field_width="short",
        ),
        NumberParameter(
            name="minimum_two_theta_deg",
            label="Ignore below",
            help_text=(
                "Discard everything below this angle. The low-angle end of a laboratory scan "
                "often carries a beam-stop shadow or an air-scatter rise that is not "
                "diffraction, and fitting it as peaks penalizes every candidate equally and "
                "wrongly."
            ),
            units="° 2θ",
            default=0.0,
            minimum=0.0,
            maximum=90.0,
            advanced=True,
            group="Scoring",
            field_width="short",
        ),
        NumberParameter(
            name="strong_line_threshold",
            label="Strong-line threshold",
            help_text=(
                "Relative intensity above which a predicted line is one the operator would "
                "expect to see, and so counts towards the lines-seen criterion. Raise it for a "
                "noisy scan in which weak calculated lines genuinely could not have been "
                "detected."
            ),
            default=0.05,
            minimum=0.005,
            maximum=0.5,
            advanced=True,
            group="Scoring",
            field_width="short",
        ),
        NumberParameter(
            name="cell_scale_range",
            label="Cell dilation searched",
            help_text=(
                "How far each candidate's cell may be uniformly dilated before matching, as a "
                "fraction. A CIF records the cell of somebody else's specimen; yours is a "
                "different composition, at a different temperature, possibly stressed, and by "
                "Δ2θ = 2·e·tanθ a difference of three parts in a thousand moves a "
                "back-reflection line by more than half a degree. Without this the true phase "
                "loses exactly the high-angle lines that would have confirmed it. A uniform "
                "dilation preserves every ratio of d spacings — which is what indexing "
                "tests — so it cannot make a wrong structure fit, and the factor each "
                "candidate needed is reported. Set it to zero to match the CIF cells exactly."
            ),
            default=0.02,
            minimum=0.0,
            maximum=0.1,
            advanced=True,
            group="Scoring",
            field_width="short",
        ),
        IntegerParameter(
            name="max_index",
            label="Largest index enumerated",
            help_text="Largest |h|, |k|, |l| generated for every candidate.",
            default=6,
            minimum=2,
            maximum=12,
            advanced=True,
            group="Scoring",
            field_width="short",
        ),
        NumberParameter(
            name="minimum_score",
            label="Acceptance threshold",
            help_text=(
                "The score the best candidate must reach before the identification is called "
                "conclusive. Below it, the verdict says that none of the candidates offered "
                "accounts for the pattern."
            ),
            default=0.55,
            minimum=0.1,
            maximum=0.95,
            advanced=True,
            group="Verdict",
            field_width="short",
        ),
        NumberParameter(
            name="decisive_margin",
            label="Decisive margin",
            help_text=(
                "The lead over the runner-up below which the top two are reported as not "
                "distinguished by this scan."
            ),
            default=0.05,
            minimum=0.005,
            maximum=0.5,
            advanced=True,
            group="Verdict",
            field_width="short",
        ),
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text=(
                "The wavelength the scan was measured with. It converts every angle into a "
                "spacing, so a wrong choice moves every candidate's calculated lines together "
                "and can make the true phase look wrong."
            ),
            options=(
                ("cu_ka", "Cu Kα (single averaged line)", "One copper line."),
                ("cu_ka_doublet", "Cu Kα1/Kα2", "Common laboratory copper doublet."),
                ("co_ka_doublet", "Co Kα1/Kα2", "Reduces Fe fluorescence."),
                ("mo_ka_doublet", "Mo Kα1/Kα2", "Short-wavelength molybdenum."),
            ),
            default="cu_ka",
            advanced=True,
        ),
    ),
    returns=(
        "The ranked candidates with their four criteria, the verdict on the best match, the "
        "fitted peaks, and the calculated line positions of every candidate for overlay."
    ),
    panel="xrd",
    citations=(
        _CITATION_HANAWALT,
        _CITATION_SMITH_SNYDER_FN,
        _CITATION_DOLLASE_MARCH,
        _CITATION_CULLITY_CH14,
    ),
    tags=(
        "XRD",
        "phase identification",
        "search match",
        "CIF",
        "indexing",
        "experimental data",
    ),
)
def _phase_identification(request: dict[str, Any]) -> dict[str, Any]:
    spec, demonstration_phase = phase_from_request(request["phase"])
    radiation = _RADIATION[str(request["radiation"])]()
    measured, generated = _measured_from_request(request, demonstration_phase, radiation)
    named, sources = _candidate_phases(request)

    floor = float(request["minimum_two_theta_deg"])
    try:
        identification, table = identify_phase_from_pattern(
            measured,
            named,
            radiation=radiation,
            sources=sources,
            prominence_sigma=float(request["prominence_sigma"]),
            minimum_two_theta_deg=floor if floor > 0.0 else None,
            tolerance_deg=float(request["tolerance_deg"]),
            max_index=int(request["max_index"]),
            strong_line_threshold=float(request["strong_line_threshold"]),
            cell_scale_range=float(request["cell_scale_range"]),
            weights=_WEIGHTING_PRESETS[str(request["weighting"])],
            minimum_score=float(request["minimum_score"]),
            decisive_margin=float(request["decisive_margin"]),
            name=f"{measured.name} phase identification",
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The pattern could not be identified: {error}",
            field="prominence_sigma",
            hint=(
                "If no peak was detected, lower the detection threshold. If every candidate was "
                "refused, check the radiation and the angular range of the scan."
            ),
        ) from error

    best = identification.best
    rows = []
    for position, candidate in enumerate(identification.candidates, start=1):
        merit = (
            None if candidate.indexing is None else float(candidate.indexing.figure_of_merit_m()[0])
        )
        unindexed = 0 if candidate.indexing is None else len(candidate.indexing.unindexed_peaks)
        rows.append(
            {
                "rank": position,
                "phase_name": candidate.phase_name,
                "score": float(candidate.score),
                "explained": 100.0 * float(candidate.explained_intensity_fraction),
                "completeness": 100.0 * float(candidate.completeness),
                "position": float(candidate.position_score),
                "intensity": float(candidate.intensity_agreement),
                "indexed": (
                    candidate.rejection
                    if candidate.indexing is None
                    else f"{candidate.indexed_count} of {candidate.indexed_count + unindexed}"
                ),
                "cell_dilation_percent": 1.0e2 * (float(candidate.cell_scale) - 1.0),
                "figure_of_merit_m": merit,
                "source": candidate.source,
            }
        )

    # Every candidate's calculated line positions travel with the result, so the
    # plot can overlay the runner-up on the scan as well as the winner. Reading
    # *where* the loser's lines fall is how a user checks the ranking rather
    # than trusting it.
    overlays = []
    for candidate in identification.candidates:
        if candidate.indexing is None:
            continue
        overlays.append(
            {
                "phase_name": candidate.phase_name,
                "score": float(candidate.score),
                "two_theta_deg": [
                    float(item.two_theta_calculated_deg) for item in candidate.indexing
                ],
                "labels": [
                    _powder_label(item.miller_indices, spec=spec) for item in candidate.indexing
                ],
                "relative_intensity": [
                    float(item.relative_intensity_calculated) for item in candidate.indexing
                ],
            }
        )

    verdict = (
        "conclusive and decisive"
        if identification.is_conclusive and identification.is_decisive
        else (
            "not decisive"
            if identification.is_conclusive
            else "not conclusive"
        )
    )
    if identification.is_conclusive and identification.is_decisive:
        summary = (
            f"{best.phase_name} at a score of {best.score:.3f}"
            + (
                ""
                if identification.runner_up is None
                else f", clear of {identification.runner_up.phase_name} by "
                f"{identification.margin:.3f}"
            )
            + f", from {identification.peak_count} fitted peaks."
        )
    elif identification.is_conclusive:
        summary = (
            f"{best.phase_name} scores best at {best.score:.3f}, but "
            f"{identification.runner_up.phase_name if identification.runner_up else 'the next'} "
            "is too close to be told apart on this scan."
        )
    else:
        summary = (
            f"No candidate accounts for this pattern. The best, {best.phase_name}, reaches only "
            f"{best.score:.3f} against a {identification.minimum_score:.2f} threshold."
        )

    notes: list[str] = []
    if generated:
        notes.append(
            "This scan was generated, not measured: a synthetic profile of the demonstration "
            f"specimen ({spec.name}) with its cell dilated by {_DEMO_LATTICE_SCALE:g}, a "
            f"{_DEMO_ZERO_SHIFT_DEG:g}° detector zero error, a curved background and "
            "Poisson counting noise. The answer is therefore known in advance, which is what "
            "makes it useful for learning: the dilation is larger than many real alloying "
            "effects, so watch how the matching tolerance has to accommodate it."
        )
    if best.indexing is not None and best.indexing.unindexed_peaks:
        notes.append(
            f"{len(best.indexing.unindexed_peaks)} measured peaks are not explained by the best "
            "match. Check whether any is strong: a strong unexplained peak means a second phase, "
            "and the quantitative step is then a multi-phase Rietveld refinement rather than a "
            "further search."
        )
    notes.append(
        "Read the losing candidates' criteria, not only their totals. A candidate that fails on "
        "lines-seen while scoring well on position has the right cell metric and the wrong "
        "centring or basis — a different fault from one that fails on position, which has "
        "the wrong cell dimensions."
    )
    notes.append(
        "The candidates must be supplied. This ranks what it is given and cannot propose a phase "
        "nobody offered, so a low best score is as likely to mean the right structure is missing "
        "from the list as that the scan is poor."
    )

    return AppResult(
        title=f"Phase identification of {measured.name}",
        summary=summary,
        table=ResultTable(
            columns=_CANDIDATE_COLUMNS,
            rows=tuple(rows),
            caption=(
                f"{len(identification.candidates)} candidates ranked against "
                f"{identification.peak_count} fitted peaks within "
                f"{float(request['tolerance_deg']):.2f}° 2θ."
            ),
        ),
        data={
            "two_theta_deg": [float(value) for value in measured.two_theta_deg],
            "observed": [float(value) for value in measured.intensity],
            "peaks": [
                {
                    "two_theta_deg": float(peak.two_theta_deg),
                    "height": float(peak.height),
                    "integrated_intensity": float(peak.integrated_intensity),
                    "fwhm_deg": float(peak.fwhm_deg),
                }
                for peak in table
            ],
            "candidates": [item.to_json() for item in identification.candidates],
            "overlays": overlays,
            "best_phase_name": best.phase_name,
            "best_score": float(best.score),
            "margin": float(identification.margin),
            "is_conclusive": identification.is_conclusive,
            "is_decisive": identification.is_decisive,
            "verdict": verdict,
            "peak_count": int(identification.peak_count),
            "synthetic": generated,
            "describe": identification.describe(),
            "best_describe": best.describe(),
            "columns": [column.to_json() for column in _CANDIDATE_COLUMNS],
        },
        inputs={
            "candidates": request["candidates"],
            "phase": spec.to_json(),
            "data_source": request["data_source"],
            "radiation": request["radiation"],
            "weighting": request["weighting"],
            "tolerance_deg": float(request["tolerance_deg"]),
            "prominence_sigma": float(request["prominence_sigma"]),
            "minimum_two_theta_deg": floor,
            "strong_line_threshold": float(request["strong_line_threshold"]),
            "cell_scale_range": float(request["cell_scale_range"]),
            "max_index": int(request["max_index"]),
            "minimum_score": float(request["minimum_score"]),
            "decisive_margin": float(request["decisive_margin"]),
            "demonstration_seed": int(request["demonstration_seed"]),
        },
        notes=tuple(notes),
        citations=(
            _CITATION_HANAWALT,
            _CITATION_SMITH_SNYDER_FN,
            _CITATION_DOLLASE_MARCH,
            _CITATION_CULLITY_CH14,
        ),
    ).to_json()
