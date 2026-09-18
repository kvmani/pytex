# ruff: noqa: RUF001
"""Powder X-ray diffraction for the shared web and desktop workbench.

The application layer does not implement diffraction physics.  It validates a
human-scale request, calls :func:`pytex.diffraction.xrd.generate_xrd_pattern`,
and turns the resulting reflection objects and sampled profile into the common
``AppResult`` contract used by tables, hover cards, and exports.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
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
from pytex.app.results import (
    AppResult,
    Column,
    ResultFigure,
    ResultMetric,
    ResultStage,
    ResultTable,
)
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.app.services.xrd_figures import (
    background_figure,
    identification_figures,
    pattern_figures,
    rietveld_figures,
    size_strain_figures,
    size_strain_highlights,
    williamson_hall_uncertainties,
)
from pytex.app.services.xrd_lattice_report import (
    EXTRAPOLATION_NAMES,
    correction_figure,
    correlation_figure,
    cross_check_figure,
    extrapolation_figure,
    indexing_figure,
    lattice_highlights,
    lattice_warnings,
    le_bail_figure,
    normalized_residual_figure,
    normalized_residuals,
    peak_overview_figure,
    peak_quality_figure,
    peak_windows_figure,
    residual_figure,
    scan_figure,
)
from pytex.app.uploads import uploaded_file
from pytex.diffraction.rietveld import _scaled_phase, refine_rietveld
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_background import estimate_background
from pytex.diffraction.xrd_corrections import specimen_displacement_shift_deg
from pytex.diffraction.xrd_indexing import PeakIndexing
from pytex.diffraction.xrd_instrument import (
    InstrumentBroadening,
    calibrate_instrument_broadening,
    deconvolve_instrument_width,
    scherrer_size_nm,
    williamson_hall,
)
from pytex.diffraction.xrd_lattice_parameter import (
    LatticeParameterResult,
    crystal_system_of,
    determine_lattice_parameters,
    extrapolation_values,
    lattice_parameter_pipeline,
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


def _powder_label(indices: tuple[int, int, int], *, spec: Any, style: str = "plain") -> str:
    """Use the conventional descending positive representative for cubic peaks.

    ``style="mathtext"`` is for figure labels, where a negative index is drawn
    overbarred as publication output requires; tables and prose keep plain text.
    """

    display_indices = indices
    if spec.crystal_system == "cubic":
        display_indices = cast(
            tuple[int, int, int],
            tuple(sorted((abs(value) for value in indices), reverse=True)),
        )
    return plane_label(display_indices, spec=spec, style=style)


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
    result = replace(
        result,
        figures=pattern_figures(
            pattern.two_theta_grid_deg,
            pattern.intensity_grid,
            rows,
            labels=[
                _powder_label(tuple(reflection.miller_indices), spec=spec, style="mathtext")
                for reflection in pattern.reflections
            ],
            radiation=str(radiation.name),
        ),
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
    result = replace(
        result,
        figures=(
            background_figure(
                estimate.two_theta_deg,
                estimate.observed_intensity,
                estimate.background,
                method=method,
                fraction=float(estimate.background_fraction),
            ),
        ),
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
    result = replace(
        result,
        figures=rietveld_figures(
            result_object,
            labels=[
                _powder_label(
                    cast(
                        tuple[int, int, int],
                        tuple(int(value) for value in reflection.miller_indices),
                    ),
                    spec=spec,
                    style="mathtext",
                )
                for reflection in result_object.reflections
            ],
        ),
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
    statistics = williamson_hall_uncertainties(
        analysis.abscissa,
        analysis.ordinate,
        wavelength_angstrom=wavelength,
        shape_factor=shape_factor,
    )
    result = replace(
        result,
        highlights=size_strain_highlights(statistics, float(analysis.r_squared)),
        figures=size_strain_figures(
            standard_angles=standard_angles,
            standard_widths=standard_widths,
            instrument=instrument,
            sample_angles=np.asarray(analysis.two_theta_deg, dtype=float),
            observed_widths=observed_sorted,
            sample_widths=np.asarray(analysis.sample_fwhm_deg, dtype=float),
            abscissa=np.asarray(analysis.abscissa, dtype=float),
            ordinate=np.asarray(analysis.ordinate, dtype=float),
            scherrer_nm=np.asarray(scherrer, dtype=float),
            statistics=statistics,
        ),
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
        "normalized_residual",
        "Residual / \u03c3",
        numeric=True,
        digits=2,
        help_text=(
            "(2\u03b8obs \u2212 2\u03b8calc) / \u03c3(2\u03b8). Within \u00b12 for about 95 % of "
            "reflections when the model and the uncertainties are right; beyond \u00b13 is an "
            "outlier worth checking."
        ),
    ),
    Column(
        "systematic_shift_mdeg",
        "Systematic correction",
        units="m\u00b0",
        numeric=True,
        digits=2,
        help_text=(
            "The angle-dependent correction fitted with the cell, at this reflection. Compare "
            "it with \u03c3(2\u03b8): if it is much larger, the correction did real work."
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
        result, indexing, peak_table, passes = lattice_parameter_pipeline(
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
    labels = [_powder_label(indices, spec=spec) for indices in result.miller_indices]
    figure_labels = [
        _powder_label(indices, spec=spec, style="mathtext") for indices in result.miller_indices
    ]
    sigma_mdeg = np.array(
        [sigma_by_angle.get(round(float(angle), 9), 0.0) for angle in result.two_theta_deg],
        dtype=float,
    )
    normalized = (
        None
        if method == "le_bail"
        else normalized_residuals(result, sigma_mdeg / 1000.0)
    )
    rows = tuple(
        {
            "hkl_label": labels[index],
            "two_theta_observed_deg": float(angle),
            "standard_uncertainty_mdeg": float(sigma_mdeg[index]),
            "two_theta_calculated_deg": float(calculated[index]),
            "residual_mdeg": 1000.0 * float(result.residual_two_theta_deg[index]),
            "normalized_residual": (
                None
                if normalized is None or not np.isfinite(normalized[index])
                else float(normalized[index])
            ),
            "systematic_shift_mdeg": 1000.0 * float(shifts[index]),
            "d_observed_angstrom": float(
                wavelength / (2.0 * np.sin(np.deg2rad(0.5 * float(angle))))
            ),
        }
        for index, angle in enumerate(result.two_theta_deg)
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
    extrapolation_plot = None
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
        # sigma(a_i) / a_i = cot(theta) sigma(theta), with sigma(theta) half of
        # sigma(2 theta): the first-order propagation of Bragg's law.
        half_angle = np.deg2rad(0.5 * result.two_theta_deg)
        per_sigma = per_reflection * np.abs(np.cos(half_angle) / np.sin(half_angle)) * (
            np.deg2rad(0.5 * sigma_mdeg / 1000.0)
        )
        extrapolation_plot = extrapolation_figure(
            abscissa=abscissa,
            per_reflection=per_reflection,
            per_reflection_sigma=per_sigma,
            labels=figure_labels,
            slope=float(slope),
            intercept=float(intercept),
            reported=float(result.a),
            reported_sigma=float(result.a_standard_uncertainty),
            function_label=_EXTRAPOLATION_LABELS[function],
        )
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
                f"A specimen displacement of {1000.0 * displacement:.0f} \u00b5m (240 mm "
                f"goniometer radius) was added to the demonstration scan. It moves the "
                f"lowest-angle reflection by {worst:.0f} m\u00b0 and higher-angle ones by less. "
                "The scan also carries a constant zero offset. The two errors have different "
                "angular forms, so no single correction function removes both exactly; careful "
                "work calibrates the zero against a standard first."
            )
        notes.append(
            "The true cell of the demonstration scan is a = "
            f"{phase.lattice.a * _DEMO_LATTICE_SCALE:.5f} \u00c5. Compare the reported value "
            "with it, then change the method and see how it moves."
        )
    notes.append(
        "The \u00b1 values are precision: the standard uncertainty of this fit on this scan. "
        "Accuracy also needs the instrument calibrated against a certified standard such as "
        "NIST SRM 640 (silicon) or SRM 660 (LaB\u2086)."
    )
    notes.append(
        "This is a lattice parameter, not a stress. A symmetric \u03b8\u20132\u03b8 scan "
        "measures only planes parallel to the specimen surface; a stress needs measurements at "
        "several specimen tilts and the X-ray elastic constants of the reflection used."
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
            f"The {result.reflection_count} reflections behind the determined cell: observed "
            "and calculated position, the residual in m\u00b0 and in units of \u03c3, and the "
            "systematic correction at each."
        )

    used_angles = {round(float(angle), 9) for angle in result.two_theta_deg}
    used_peaks = (
        []
        if indexing is None
        else [
            item.peak
            for item in indexing.reflections
            if round(float(item.peak.two_theta_deg), 9) in used_angles
        ]
    )
    warnings = lattice_warnings(
        result,
        normalized=normalized,
        labels=labels,
        peaks_used=used_peaks,
        figure_of_merit=None if indexing is None else float(indexing.figure_of_merit_m()[0]),
        unindexed_count=0 if indexing is None else len(indexing.unindexed_peaks),
        generated=generated,
    )
    profile_points = (
        None if result.profile_two_theta_deg is None else int(result.profile_two_theta_deg.size)
    )
    highlights = lattice_highlights(result, system=system, profile_points=profile_points)
    cell_text = f"a = {result.a:.6f} \u00b1 {result.a_standard_uncertainty:.6f} \u00c5" + (
        ""
        if system == "cubic"
        else f", c = {result.c:.6f} \u00b1 {result.c_standard_uncertainty:.6f} \u00c5 "
        f"(c/a = {result.axial_ratio:.6f})"
    )
    if method == "le_bail":
        how = (
            f"from a Le Bail whole-pattern fit of {profile_points} profile points modelling "
            f"{result.reflection_count} reflections (profile reduced \u03c7\u00b2 = "
            f"{result.reduced_chi_squared:.2f})"
        )
    else:
        if method == "average":
            technique = "by averaging the value from each reflection"
        elif result.extrapolation == "none":
            technique = "by Cohen least squares without a systematic correction"
        else:
            technique = (
                "by Cohen least squares with a "
                f"{EXTRAPOLATION_NAMES[result.extrapolation]} systematic correction"
            )
        free = len(result.free_parameter_names) + (
            0 if method == "average" or result.extrapolation == "none" else 1
        )
        how = (
            f"from {result.reflection_count} reflections {technique}, leaving "
            f"{result.reflection_count - free} degrees of freedom and a lattice-fit reduced "
            f"\u03c7\u00b2 of {result.reduced_chi_squared:.2f}"
        )
    summary = (
        f"{cell_text}, {how}. The \u00b1 values are one standard uncertainty: the precision "
        f"of this fit (relative precision {result.relative_uncertainty:.1e}), not its accuracy."
        + (
            f" {len(warnings)} check{'s' if len(warnings) != 1 else ''} below "
            f"{'need' if len(warnings) != 1 else 'needs'} attention before the value is used."
            if warnings
            else " No reliability check failed."
        )
    )

    result_payload = AppResult(
        title=f"Lattice parameters of {spec.name}",
        summary=summary,
        highlights=highlights,
        warnings=warnings,
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
            # The same number under the name that says what it is: a change
            # against the tabulated cell, which is an elastic strain only when
            # that cell is the stress-free one of this material.
            "relative_change_from_reference": result.strain_relative_to_reference,
            "warnings": list(warnings),
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
        stages=_lattice_stages(
            request=request,
            measured=measured,
            generated=generated,
            radiation=radiation,
            displacement=displacement,
            spec=spec,
            phase=phase,
            system=system,
            result=result,
            indexing=indexing,
            peak_table=peak_table,
            passes=passes,
            labels=figure_labels,
            sigma_mdeg=sigma_mdeg,
            normalized=normalized,
            extrapolation_plot=extrapolation_plot,
        ),
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

_METHOD_LABELS = {
    "cohen": "Cohen least squares",
    "average": "Average over reflections",
    "le_bail": "Le Bail whole-pattern decomposition",
}

_STAGE_PEAK_COLUMNS = (
    Column("two_theta_deg", "2θ", units="°", numeric=True, digits=4),
    Column(
        "sigma_mdeg",
        "σ(2θ)",
        units="m°",
        numeric=True,
        digits=2,
        help_text="Standard uncertainty of the fitted centre, from the fit covariance.",
    ),
    Column("height", "Height", numeric=True, digits=1),
    Column("integrated_intensity", "Integrated intensity", numeric=True, digits=1),
    Column("fwhm_deg", "FWHM", units="°", numeric=True, digits=4),
    Column(
        "eta",
        "η",
        numeric=True,
        digits=3,
        help_text="Lorentzian fraction of the pseudo-Voigt: 0 is Gaussian, 1 is Lorentzian.",
    ),
    Column(
        "reduced_chi_squared",
        "χ²ν",
        numeric=True,
        digits=2,
        help_text="Reduced chi-squared of this peak's profile fit; about 1 is ideal.",
    ),
    Column("converged", "Converged"),
)

_STAGE_PASS_COLUMNS = (
    Column("pass", "Pass", numeric=True),
    Column(
        "a_angstrom",
        "a",
        units="Å",
        numeric=True,
        digits=6,
        help_text=(
            "Cell edge determined from this pass's assignment. For a pass that was not taken, "
            "the cell it was indexed against."
        ),
    ),
    Column("indexed_count", "Indexed", numeric=True),
    Column("unindexed_count", "Unindexed", numeric=True),
    Column(
        "figure_of_merit_m",
        "M",
        numeric=True,
        digits=1,
        help_text=(
            "de Wolff's figure of merit for the pass: above 10 plausible, above 20 convincing."
        ),
    ),
    Column("mean_delta_mdeg", "Mean |Δ2θ|", units="m°", numeric=True, digits=1),
    Column("outcome", "Outcome"),
)

_STAGE_ASSIGNMENT_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column("two_theta_observed_deg", "2θ observed", units="°", numeric=True, digits=4),
    Column("sigma_mdeg", "σ(2θ)", units="m°", numeric=True, digits=2),
    Column(
        "two_theta_calculated_deg",
        "2θ calculated",
        units="°",
        numeric=True,
        digits=4,
        help_text="From the cell of the final indexing pass, before any systematic correction.",
    ),
    Column("delta_mdeg", "Δ2θ", units="m°", numeric=True, digits=1),
    Column("d_observed_angstrom", "d observed", units="Å", numeric=True, digits=5),
    Column("d_calculated_angstrom", "d calculated", units="Å", numeric=True, digits=5),
    Column("multiplicity", "Multiplicity", numeric=True),
    Column(
        "relative_intensity",
        "I calculated",
        numeric=True,
        digits=3,
        help_text="Calculated relative intensity; used to rank lines, never to fit the cell.",
    ),
)

_STAGE_CROSS_CHECK_COLUMNS = (
    Column("method", "Method"),
    Column("a_angstrom", "a", units="Å", numeric=True, digits=6),
    Column("sigma_angstrom", "σ(a)", units="Å", numeric=True, digits=6),
    Column(
        "difference_ppm",
        "Difference from reported a",
        units="ppm",
        numeric=True,
        digits=0,
        help_text="(a_method − a_reported) / a_reported, in parts per million.",
    ),
    Column("reduced_chi_squared", "χ²ν", numeric=True, digits=3),
)


def _scan_stage(
    *,
    request: Mapping[str, Any],
    measured: MeasuredPowderPattern,
    generated: bool,
    radiation: RadiationSpec,
    displacement: float | None,
) -> ResultStage:
    """The first stage of every scan analysis: what was read, and at what wavelength."""

    wavelength = float(radiation.wavelength_angstrom)
    axis = np.asarray(measured.two_theta_deg, dtype=float)
    step = float(np.median(np.diff(axis))) if axis.size > 1 else float("nan")
    if generated:
        source = "Generated demonstration scan"
    else:
        kind = {"file": "Pattern file", "paste": "Pasted scan"}.get(
            str(request.get("data_source")), "Measured scan"
        )
        source = f"{kind} '{measured.name}'"
    metrics = [
        ResultMetric("Points", int(axis.size)),
        ResultMetric("First 2θ", float(axis[0]), "°"),
        ResultMetric("Last 2θ", float(axis[-1]), "°"),
        ResultMetric(
            "Median step",
            step,
            "°",
            "Needs to be several times smaller than the peak width for a centre to be located "
            "to a small fraction of a step.",
        ),
        ResultMetric("Radiation", str(radiation.name)),
        ResultMetric("Wavelength λ", wavelength, "Å"),
    ]
    if displacement is not None:
        metrics.append(ResultMetric("Injected specimen displacement", displacement, "mm"))
    return ResultStage(
        key="scan",
        title="1. The scan as read",
        summary=(
            f"{source}: {axis.size} points from {axis[0]:.3f}° to {axis[-1]:.3f}° 2θ in steps "
            f"of {step:.4f}°, referred to {radiation.name} (λ = {wavelength:.6f} Å)."
        ),
        metrics=tuple(metrics),
        explanation=(
            "Every spacing below is computed from this wavelength, so the wrong radiation "
            "scales every cell edge by the ratio of the two wavelengths. The step matters "
            "because a profile fit can only place a centre to a small fraction of a step "
            "when there are at least about five points across a peak's full width at half "
            "maximum."
        ),
        status="info",
    )


def _peak_stage_rows(peaks: Any) -> tuple[dict[str, Any], ...]:
    """One row per fitted peak, in the columns of :data:`_STAGE_PEAK_COLUMNS`."""

    return tuple(
        {
            "two_theta_deg": float(peak.two_theta_deg),
            "sigma_mdeg": 1000.0 * float(peak.two_theta_standard_uncertainty_deg),
            "height": float(peak.height),
            "integrated_intensity": float(peak.integrated_intensity),
            "fwhm_deg": float(peak.fwhm_deg),
            "eta": float(peak.eta),
            "reduced_chi_squared": float(peak.reduced_chi_squared),
            "converged": bool(peak.converged),
        }
        for peak in peaks
    )


def _lattice_stages(
    *,
    request: Mapping[str, Any],
    measured: MeasuredPowderPattern,
    generated: bool,
    radiation: RadiationSpec,
    displacement: float,
    spec: Any,
    phase: Any,
    system: str,
    result: LatticeParameterResult,
    indexing: PeakIndexing | None,
    peak_table: Any,
    passes: tuple[Mapping[str, float | bool], ...],
    labels: list[str],
    sigma_mdeg: np.ndarray,
    normalized: np.ndarray | None,
    extrapolation_plot: ResultFigure | None,
) -> tuple[ResultStage, ...]:
    """Report a lattice-parameter determination as a readable, checkable report.

    The stages are returned in the order the report is read, not the order the
    computation ran: the determined cell first, then the evidence (the scan,
    the fitted peaks, the indexing), then the diagnostics (the lattice-fit
    residuals and systematic correction, the peak-fit quality, the other
    methods), then the method, then the audit trail (the indexing passes). Each
    stage names its section, carries its own numbers and figures, and says how
    to read them, because the step that went wrong is not visible from the cell
    alone.
    """

    axis = np.asarray(measured.two_theta_deg, dtype=float)
    counts = np.asarray(measured.intensity, dtype=float)
    unit = str(measured.intensity_unit)
    source = (
        "Generated demonstration scan"
        if generated
        else {"file": "Pattern file", "paste": "Pasted scan"}.get(
            str(request.get("data_source")), "Measured scan"
        )
        + f" '{measured.name}'"
    )
    unindexed_angles = (
        [] if indexing is None else [float(peak.two_theta_deg) for peak in indexing.unindexed_peaks]
    )
    detected_angles = [] if peak_table is None else [float(p.two_theta_deg) for p in peak_table]
    scan = replace(
        _scan_stage(
            request=request,
            measured=measured,
            generated=generated,
            radiation=radiation,
            displacement=displacement if generated else None,
        ),
        title="Measured scan",
        section="evidence",
        figures=(
            scan_figure(
                axis,
                counts,
                unit=unit,
                peak_angles=detected_angles,
                unindexed_angles=unindexed_angles,
                source=source,
            ),
        ),
    )
    evidence: list[ResultStage] = [scan]
    diagnostics: list[ResultStage] = []
    method_stages: list[ResultStage] = []
    audit: list[ResultStage] = []

    if result.method == "le_bail":
        systematic = str(request.get("systematic", "none"))
        units = {"displacement": "mm", "zero": "° 2θ"}.get(systematic)
        whole_metrics = [
            ResultMetric("Reflections modelled", int(result.reflection_count)),
            ResultMetric(
                "Profile-fit reduced χ²",
                float(result.reduced_chi_squared),
                None,
                "Calculated profile against every measured point; about 1 is ideal.",
            ),
            ResultMetric(
                "R_wp (background removed)",
                result.weighted_profile_r,
                None,
                "Computed on the background-subtracted profile, so it is systematically higher "
                "than a Rietveld program's R_wp on the raw scan and must not be compared with it.",
            ),
            ResultMetric("Systematic term refined", systematic),
        ]
        if units is not None:
            whole_metrics.append(
                ResultMetric("Refined systematic value", float(result.drift_coefficient), units)
            )
        evidence.append(
            ResultStage(
                key="whole_pattern",
                title="Whole-pattern (Le Bail) fit",
                summary=(
                    f"The Le Bail fit modelled {result.reflection_count} reflections against "
                    f"every measured point and reached a profile reduced χ² of "
                    f"{result.reduced_chi_squared:.3f}"
                    + (
                        ""
                        if result.weighted_profile_r is None
                        else f" with R_wp = {100.0 * result.weighted_profile_r:.2f} %"
                    )
                    + "."
                ),
                metrics=tuple(whole_metrics),
                explanation=(
                    "No individual peak is located. The intensity of each reflection is "
                    "extracted from the measured profile while the cell, the peak shapes and one "
                    "systematic term are refined, so neither texture nor the atomic basis can "
                    "bias the cell. The diagnostic is therefore the difference curve, not a "
                    "per-reflection residual: structure left in it is structure the model did "
                    "not describe."
                ),
                status="warning" if result.reduced_chi_squared > 3.0 else "ok",
                section="evidence",
                figures=(
                    le_bail_figure(
                        result,
                        calculated_angles=[float(value) for value in result.two_theta_deg],
                        calculated_labels=labels,
                    ),
                ),
            )
        )
    else:
        assert indexing is not None
        assert peak_table is not None
        peaks = list(peak_table)
        converged = int(peak_table.converged_count)
        parameter_count = len(result.free_parameter_names) + (
            0 if result.extrapolation == "none" else 1
        )
        shape = str(peak_table.settings.get("shape", "pseudo_voigt")).replace("_", "-")
        doublet = bool(peak_table.settings.get("model_doublet", False))
        sigmas = np.array([1000.0 * peak.two_theta_standard_uncertainty_deg for peak in peaks])
        widths = np.array([peak.fwhm_deg for peak in peaks])
        label_by_angle = {
            round(float(item.peak.two_theta_deg), 9): _powder_label(
                item.miller_indices, spec=spec, style="mathtext"
            )
            for item in indexing.reflections
        }
        peak_labels = [label_by_angle.get(round(float(p.two_theta_deg), 9), "") for p in peaks]
        calculated_angles = [float(item.two_theta_calculated_deg) for item in indexing.reflections]
        calculated_labels = [
            _powder_label(item.miller_indices, spec=spec, style="mathtext")
            for item in indexing.reflections
        ]
        evidence.append(
            ResultStage(
                key="peaks",
                title="Peak positions from profile fitting",
                summary=(
                    f"{len(peaks)} peaks were found and each was fitted with a {shape} profile"
                    + (" including its Kα₂ partner" if doublet else "")
                    + f"; {converged} of {len(peaks)} fits converged. The median uncertainty of "
                    f"a peak position is {float(np.median(sigmas)):.2f} m° and the median width "
                    f"{float(np.median(widths)):.4f}°."
                ),
                metrics=(
                    ResultMetric("Peaks detected", len(peaks)),
                    ResultMetric("Fits converged", converged),
                    ResultMetric(
                        "Detection threshold", float(request["prominence_sigma"]), "σ"
                    ),
                    ResultMetric("Expected width", float(request["expected_fwhm_deg"]), "°"),
                    ResultMetric("Median σ(2θ)", float(np.median(sigmas)), "m°"),
                    ResultMetric("Median FWHM", float(np.median(widths)), "°"),
                ),
                table=ResultTable(
                    columns=_STAGE_PEAK_COLUMNS,
                    rows=_peak_stage_rows(peaks),
                    caption=(
                        "Every peak kept, in ascending 2θ. χ²ν here is the peak-profile fit, "
                        "not the lattice fit."
                    ),
                ),
                explanation=(
                    "Peaks are found by filtering the scan with a kernel matched to the expected "
                    "width and keeping maxima above the threshold. Each peak is then fitted in "
                    "its own window with a pseudo-Voigt profile on a straight local background. "
                    "σ(2θ) is the standard uncertainty of the fitted centre; it becomes the "
                    "weight of that reflection in the lattice fit, so an imprecise peak counts "
                    "for little. The peak-profile χ²ν compares one profile with the counts in "
                    "its window: near 1 is a good description; far above 1 flags an overlapped, "
                    "asymmetric or badly backgrounded peak whose position deserves suspicion."
                ),
                status=(
                    "warning"
                    if converged < len(peaks) or len(peaks) <= parameter_count
                    else "ok"
                ),
                section="evidence",
                figures=(
                    peak_overview_figure(
                        axis,
                        counts,
                        unit=unit,
                        peaks=peaks,
                        radiation=radiation,
                        calculated_angles=calculated_angles,
                        calculated_labels=calculated_labels,
                    ),
                    peak_windows_figure(
                        axis, counts, peaks=peaks, labels=peak_labels, radiation=radiation
                    ),
                ),
            )
        )

        merit_m, count_m = indexing.figure_of_merit_m()
        merit_f, count_f = indexing.figure_of_merit_f()
        unindexed = [float(peak.two_theta_deg) for peak in indexing.unindexed_peaks]
        residuals = indexing.delta_two_theta_deg
        same_sign = residuals.size >= 3 and bool(
            np.all(residuals > 0.0) or np.all(residuals < 0.0)
        )
        assignment_summary = (
            f"{indexing.indexed_count} of {indexing.indexed_count + len(unindexed)} peaks were "
            f"matched one-to-one to reflections of {spec.name} within "
            f"±{indexing.tolerance_deg:.3f}°. The figures of merit are de Wolff "
            f"M_{count_m} = {merit_m:.1f} and Smith–Snyder F_{count_f} = {merit_f:.1f}; the "
            f"mean |Δ2θ| is {1000.0 * indexing.mean_absolute_delta_two_theta_deg:.1f} m°."
        )
        if same_sign:
            assignment_summary += (
                " Every Δ2θ has the same sign: the signature of an uncorrected zero or "
                "displacement error, which the lattice fit corrects, rather than of a wrong cell."
            )
        if unindexed:
            assignment_summary += (
                f" {len(unindexed)} peak{'s' if len(unindexed) != 1 else ''} at "
                + ", ".join(f"{value:.3f}°" for value in unindexed)
                + " stayed unindexed; a strong one belongs to something the phase does not "
                "describe — a second phase, a Kβ or tungsten line, or the sample holder."
            )
        if indexing.unobserved_indices:
            assignment_summary += (
                f" {len(indexing.unobserved_indices)} calculated reflections above the intensity "
                "threshold were not observed."
            )
        evidence.append(
            ResultStage(
                key="assignment",
                title="Peak indexing",
                summary=assignment_summary,
                metrics=(
                    ResultMetric("Indexed fraction", 100.0 * indexing.indexed_fraction, "%"),
                    ResultMetric(
                        f"de Wolff M_{count_m}",
                        float(merit_m),
                        None,
                        "Above 10 is a plausible cell, above 20 a convincing one.",
                    ),
                    ResultMetric(f"Smith–Snyder F_{count_f}", float(merit_f)),
                    ResultMetric("Unindexed peaks", len(unindexed)),
                    ResultMetric("Unobserved strong lines", len(indexing.unobserved_indices)),
                ),
                table=ResultTable(
                    columns=_STAGE_ASSIGNMENT_COLUMNS,
                    rows=tuple(
                        {
                            "hkl_label": _powder_label(item.miller_indices, spec=spec),
                            "two_theta_observed_deg": float(item.peak.two_theta_deg),
                            "sigma_mdeg": 1000.0
                            * float(item.peak.two_theta_standard_uncertainty_deg),
                            "two_theta_calculated_deg": float(item.two_theta_calculated_deg),
                            "delta_mdeg": 1000.0 * float(item.delta_two_theta_deg),
                            "d_observed_angstrom": float(item.d_observed_angstrom),
                            "d_calculated_angstrom": float(item.d_calculated_angstrom),
                            "multiplicity": int(item.multiplicity),
                            "relative_intensity": float(item.relative_intensity_calculated),
                        }
                        for item in indexing.reflections
                    ),
                    caption="The final assignment, one row per indexed peak.",
                ),
                explanation=(
                    "Δ2θ is observed minus calculated against the cell of the final indexing "
                    "pass, before any systematic correction, so a smooth trend with angle is "
                    "expected when a zero or displacement error is present — the lattice fit "
                    "removes it. What must not appear is one reflection far off the trend: that "
                    "is a misassignment. M and F both penalize a cell that fits only because it "
                    "predicts many lines. Calculated intensities rank lines and are never used "
                    "to fit the cell."
                ),
                status="warning" if (merit_m < 10.0 or unindexed) else "ok",
                section="evidence",
                figures=(
                    indexing_figure(
                        observed_angles=[
                            float(item.peak.two_theta_deg) for item in indexing.reflections
                        ],
                        observed_heights=[float(item.peak.height) for item in indexing.reflections],
                        observed_sigma_mdeg=[
                            1000.0 * float(item.peak.two_theta_standard_uncertainty_deg)
                            for item in indexing.reflections
                        ],
                        calculated_angles=calculated_angles,
                        calculated_intensities=[
                            float(item.relative_intensity_calculated)
                            for item in indexing.reflections
                        ],
                        delta_mdeg=[
                            1000.0 * float(item.delta_two_theta_deg)
                            for item in indexing.reflections
                        ],
                        labels=calculated_labels,
                        unindexed_angles=unindexed,
                        unindexed_heights=[float(p.height) for p in indexing.unindexed_peaks],
                        tolerance_deg=float(indexing.tolerance_deg),
                    ),
                ),
            )
        )

        floor_value = float(request["minimum_two_theta_deg"])
        dropped = indexing.indexed_count - result.reflection_count
        names = result.correlation_parameter_names
        degrees = result.reflection_count - (
            1 if result.method == "average" else len(names)
        )
        chi = float(result.reduced_chi_squared)
        finite = (
            np.array([], dtype=float)
            if normalized is None
            else normalized[np.isfinite(normalized)]
        )
        within_two = int(np.count_nonzero(np.abs(finite) <= 2.0))
        beyond_three = int(np.count_nonzero(np.abs(finite) > 3.0))
        shifts = 1000.0 * np.abs(result.systematic_shift_deg)
        largest_shift = float(np.max(shifts)) if shifts.size else 0.0
        fit_figures: list[ResultFigure] = [
            residual_figure(result, sigma_mdeg, labels),
        ]
        if normalized is not None:
            fit_figures.append(normalized_residual_figure(result, normalized, labels))
        if result.method == "cohen" and result.extrapolation != "none":
            fit_figures.append(
                correction_figure(
                    result,
                    two_theta_range=(float(axis[0]), float(axis[-1])),
                    sigma_mdeg=sigma_mdeg,
                    labels=labels,
                )
            )
        correlation_plot = correlation_figure(result)
        if correlation_plot is not None:
            fit_figures.append(correlation_plot)
        quality = (
            "the residuals match the peak uncertainties"
            if 0.3 <= chi <= 3.0
            else (
                "the residuals are larger than the peak uncertainties allow"
                if chi > 3.0
                else "the residuals are smaller than the peak uncertainties suggest"
            )
        )
        diagnostics.append(
            ResultStage(
                key="lattice_fit",
                title="Lattice-fit residuals and systematic correction",
                summary=(
                    f"After the fit, {within_two} of {finite.size} reflections lie within ±2σ "
                    f"of the calculated positions and {beyond_three} beyond ±3σ. The "
                    f"lattice-fit reduced χ² is {chi:.2f} on {degrees} degrees of freedom, so "
                    f"{quality}."
                    + (
                        f" The systematic correction moves reflections by up to "
                        f"{largest_shift:.1f} m°, against a median peak uncertainty of "
                        f"{float(np.median(sigma_mdeg)):.2f} m°."
                        if result.method == "cohen" and result.extrapolation != "none"
                        else ""
                    )
                ),
                metrics=(
                    ResultMetric("Reflections used", int(result.reflection_count)),
                    ResultMetric("Degrees of freedom", int(degrees)),
                    ResultMetric(
                        "Lattice-fit reduced χ²",
                        chi,
                        None,
                        "Σ[(2θobs − 2θcalc)/σ(2θ)]² / degrees of freedom. About 1: residuals "
                        "match the position uncertainties. Much larger: an unmodelled error or a "
                        "misassignment. Much smaller: overstated uncertainties.",
                    ),
                    ResultMetric("Within ±2σ", within_two),
                    ResultMetric("Beyond ±3σ", beyond_three),
                    ResultMetric(
                        "RMS residual",
                        float(np.sqrt(np.mean(np.square(1000.0 * result.residual_two_theta_deg))))
                        if result.two_theta_deg.size
                        else 0.0,
                        "m°",
                    ),
                    ResultMetric("Largest systematic correction", largest_shift, "m°"),
                ),
                explanation=(
                    "The residual is the observed position minus the position the fitted cell "
                    "and systematic correction predict. Divided by its own σ(2θ) it should "
                    "behave like a standard normal variable: about 95 % inside ±2, nearly all "
                    "inside ±3. The systematic correction is an angle-dependent term fitted "
                    "together with the cell; its form matches specimen displacement and "
                    "absorption, but the fit cannot say which aberration caused it, so it is a "
                    "correction, not a measured displacement. The correlation map shows how "
                    "strongly the scan couples each pair of refined parameters."
                ),
                status=(
                    "warning"
                    if chi > 3.0 or beyond_three > 0 or degrees < 1
                    else "ok"
                ),
                section="diagnostics",
                figures=tuple(fit_figures),
            )
        )

        used = {round(float(angle), 9) for angle in result.two_theta_deg}
        chi_peaks = np.array([float(peak.reduced_chi_squared) for peak in peaks])
        poor = int(np.count_nonzero(chi_peaks > 10.0))
        diagnostics.append(
            ResultStage(
                key="peak_quality",
                title="Peak-fit quality",
                summary=(
                    f"The peak-profile χ²ν ranges from {float(np.min(chi_peaks)):.2f} to "
                    f"{float(np.max(chi_peaks)):.2f} (median {float(np.median(chi_peaks)):.2f})"
                    + (
                        f"; {poor} peak{'s are' if poor != 1 else ' is'} fitted poorly (χ²ν > 10)"
                        if poor
                        else ""
                    )
                    + f". Widths run from {float(np.min(widths)):.4f}° to "
                    f"{float(np.max(widths)):.4f}°."
                ),
                metrics=(
                    ResultMetric("Median peak-fit χ²ν", float(np.median(chi_peaks))),
                    ResultMetric("Peaks with χ²ν > 10", poor),
                    ResultMetric("Fits not converged", len(peaks) - converged),
                ),
                explanation=(
                    "These χ² values judge how well a pseudo-Voigt describes each peak's counts, "
                    "and are unrelated to the lattice-fit χ², which judges positions against the "
                    "cell. A poorly fitted peak can still give a good position, but its σ(2θ) is "
                    "then only a lower bound. Widths growing smoothly with angle are normal "
                    "instrument and specimen broadening."
                ),
                status="warning" if poor or converged < len(peaks) else "ok",
                section="diagnostics",
                figures=(
                    peak_quality_figure(
                        angles=[float(peak.two_theta_deg) for peak in peaks],
                        sigma_mdeg=sigmas,
                        fwhm_deg=widths,
                        reduced_chi_squared=chi_peaks,
                        converged=[bool(peak.converged) for peak in peaks],
                        used=[round(float(peak.two_theta_deg), 9) in used for peak in peaks],
                    ),
                ),
            )
        )

        if result.method == "average":
            ls_summary = (
                f"The lattice parameter was computed separately from each of "
                f"{result.reflection_count} reflections and averaged. No systematic term can be "
                "refined this way, so any angle-dependent error stays inside the mean."
            )
            ls_metrics: tuple[ResultMetric, ...] = (
                ResultMetric("Method", _METHOD_LABELS[result.method]),
                ResultMetric("Reflections", int(result.reflection_count)),
                ResultMetric("Reduced χ²", float(result.reduced_chi_squared)),
            )
            ls_table = None
            strongest = 0.0
        else:
            significance = (
                abs(result.drift_coefficient) / result.drift_standard_uncertainty
                if result.drift_standard_uncertainty > 0.0
                else float("nan")
            )
            matrix = result.parameter_correlation
            strongest = 0.0
            if matrix is not None and matrix.shape[0] > 1:
                off = np.abs(matrix - np.eye(matrix.shape[0]))
                strongest = float(np.max(off))
            drift_text = (
                "No systematic correction was refined, so any zero, displacement or "
                "transparency error is inside the cell."
                if result.extrapolation == "none"
                else (
                    f"The angle-dependent systematic correction, D·sin²θ·f(θ) with "
                    f"f = {_EXTRAPOLATION_LABELS[result.extrapolation]}, refined to "
                    f"D = {result.drift_coefficient:.3e} ± "
                    f"{result.drift_standard_uncertainty:.1e} ({significance:.1f} σ)."
                )
            )
            ls_summary = (
                f"Weighted least squares in sin²θ solved for {len(names)} parameter"
                f"{'s' if len(names) != 1 else ''} ({', '.join(names)}) from "
                f"{result.reflection_count} reflections. {drift_text}"
            )
            ls_metrics = (
                ResultMetric("Method", _METHOD_LABELS[result.method]),
                ResultMetric("Extrapolation function", result.extrapolation),
                ResultMetric("Refined parameters", ", ".join(names)),
                ResultMetric("Observations", int(result.reflection_count)),
                ResultMetric("Degrees of freedom", int(degrees)),
                ResultMetric("Drift coefficient D", float(result.drift_coefficient)),
                ResultMetric("σ(D)", float(result.drift_standard_uncertainty)),
                ResultMetric(
                    "|D| / σ(D)",
                    significance,
                    None,
                    "Above about 2, the systematic term is significantly different from zero.",
                ),
                ResultMetric("Largest systematic correction", largest_shift, "m°"),
                ResultMetric("Lattice-fit reduced χ²", chi),
                ResultMetric("Angular floor", floor_value, "°"),
                ResultMetric("Reflections discarded by the floor", int(dropped)),
            )
            ls_table = None
            if matrix is not None:
                keys = [f"p{index}" for index in range(len(names))]
                ls_table = ResultTable(
                    columns=(
                        Column("parameter", "Parameter"),
                        *(
                            Column(key, name, numeric=True, digits=3)
                            for key, name in zip(keys, names, strict=True)
                        ),
                    ),
                    rows=tuple(
                        {
                            "parameter": name,
                            **{
                                key: float(matrix[row, column])
                                for column, key in enumerate(keys)
                            },
                        }
                        for row, name in enumerate(names)
                    ),
                    caption="Correlation matrix of the refined parameters (−1 to +1).",
                )
        method_stages.append(
            ResultStage(
                key="least_squares",
                title="How the cell was calculated",
                summary=ls_summary,
                metrics=ls_metrics,
                table=ls_table,
                explanation=(
                    "Bragg's law is linear in the reciprocal metric tensor G*: "
                    "sin²θ = (λ²/4)·hᵀG*h. Each reflection gives one linear equation in the "
                    "components the crystal system leaves free (A = a*², …), and one more "
                    "column, D·sin²θ·f(θ), carries the angle-dependent systematic error through "
                    "the extrapolation function f, which vanishes at θ = 90°. Each equation is "
                    "weighted by 1/σ(sin²θ), propagated from the peak's σ(2θ), so the precise "
                    "high-angle reflections dominate. The covariance is the inverse normal "
                    "matrix multiplied by the lattice-fit reduced χ², so the quoted "
                    "uncertainties grow when the residuals are larger than the peak "
                    "uncertainties predict (and shrink when they are smaller)."
                ),
                status=(
                    "warning"
                    if result.reduced_chi_squared > 3.0 or strongest > 0.98
                    else "ok"
                ),
                section="method",
            )
        )

        first = int(passes[0]["indexed_count"]) if passes else indexing.indexed_count
        final = indexing.indexed_count
        recovered = final - first
        audit.append(
            ResultStage(
                key="passes",
                title="Indexing passes",
                summary=(
                    f"{len(passes)} pass{'es' if len(passes) != 1 else ''} ran. The first "
                    f"indexed {first} reflections against the tabulated cell of {spec.name}; the "
                    f"final assignment holds {final}. "
                    + (
                        f"Re-indexing against the determined cell recovered {recovered} "
                        "reflections the starting cell had placed outside the tolerance — "
                        "typically the high-angle ones that carry most of the precision."
                        if recovered > 0
                        else "Re-indexing found nothing the first pass had missed, so the "
                        "starting cell was already close enough."
                    )
                ),
                metrics=(
                    ResultMetric("Indexing tolerance", float(request["tolerance_deg"]), "°"),
                    ResultMetric("Largest Miller index", int(request["max_index"])),
                    ResultMetric("Reflections recovered by re-indexing", recovered),
                ),
                table=ResultTable(
                    columns=_STAGE_PASS_COLUMNS,
                    rows=tuple(
                        {
                            "pass": int(item["pass"]),
                            "a_angstrom": float(item["a_angstrom"]),
                            "indexed_count": int(item["indexed_count"]),
                            "unindexed_count": int(item["unindexed_count"]),
                            "figure_of_merit_m": float(item["figure_of_merit_m"]),
                            "mean_delta_mdeg": 1000.0
                            * float(item["mean_absolute_delta_two_theta_deg"]),
                            "outcome": (
                                "taken" if item["accepted"] else "not taken: indexed no more"
                            ),
                        }
                        for item in passes
                    ),
                ),
                explanation=(
                    "A starting cell wrong by a fraction e misplaces a reflection by "
                    "Δ(2θ) = 2e·tanθ, which grows towards back-reflection, so the first pass can "
                    "silently drop exactly the high-angle lines that matter most. Each later "
                    "pass indexes against the cell the previous one determined, and the loop "
                    "stops at the first pass that indexes no more reflections than the one "
                    "before it."
                ),
                section="audit",
            )
        )

    strain = result.strain_relative_to_reference
    cell_metrics = [
        ResultMetric("a", float(result.a), "Å"),
        ResultMetric("σ(a)", float(result.a_standard_uncertainty), "Å"),
    ]
    if system != "cubic":
        cell_metrics += [
            ResultMetric("c", float(result.c), "Å"),
            ResultMetric("σ(c)", float(result.c_standard_uncertainty), "Å"),
            ResultMetric("c/a", float(result.axial_ratio)),
        ]
    if system in {"orthorhombic", "monoclinic", "triclinic"}:
        cell_metrics += [
            ResultMetric("b", float(result.b), "Å"),
            ResultMetric("σ(b)", float(result.b_standard_uncertainty), "Å"),
        ]
    if system in {"monoclinic", "triclinic"}:
        cell_metrics += [
            ResultMetric("α", float(result.alpha_deg), "°"),
            ResultMetric("β", float(result.beta_deg), "°"),
            ResultMetric("γ", float(result.gamma_deg), "°"),
        ]
    cell_metrics.append(
        ResultMetric(
            "Relative precision σ(a)/a",
            float(result.relative_uncertainty),
            None,
            "About 1e-5 is strain-grade; 1e-4 composition-grade; 1e-3 identification-grade.",
        )
    )
    if result.reference_lattice is not None:
        cell_metrics.append(
            ResultMetric(
                "Reference a of the selected phase", float(result.reference_lattice.a), "Å"
            )
        )
    if strain is not None:
        cell_metrics.append(
            ResultMetric(
                "Change from the reference cell (a)",
                float(strain),
                None,
                "(a − a_ref)/a_ref. Not an elastic strain unless a_ref is the stress-free cell "
                "of this material measured on this instrument.",
            )
        )
    grade = (
        "strain-grade (about 1e-5)"
        if result.relative_uncertainty < 5.0e-5
        else (
            "composition-grade (about 1e-4)"
            if result.relative_uncertainty < 5.0e-4
            else "identification-grade (about 1e-3)"
        )
    )
    cell_stage = ResultStage(
        key="cell",
        title="Determined cell parameters",
        summary=(
            f"a = {result.a:.6f} ± {result.a_standard_uncertainty:.6f} Å for the {system} "
            f"cell: a relative precision of {result.relative_uncertainty:.1e}, which is "
            f"{grade}."
            + (
                ""
                if strain is None
                else f" This is {strain:+.3e} relative to the tabulated cell of the selected "
                "phase — a comparison with the database value, not a measured elastic strain."
            )
        ),
        metrics=tuple(cell_metrics),
        explanation=(
            "The ± value is one standard uncertainty propagated from the fit covariance through "
            "the reciprocal-to-direct cell conversion. It is precision — how reproducible the "
            "number is on this scan — not accuracy. Accuracy needs the instrument calibrated "
            "against a certified standard such as NIST SRM 640 silicon or SRM 660 lanthanum "
            "hexaboride. The change from the reference cell includes composition, temperature "
            "and calibration differences; it becomes an elastic strain only when the reference "
            "is the stress-free cell of the same material on the same instrument, and even then "
            "a symmetric scan gives the strain normal to the surface only, not a stress."
        ),
        section="result",
        figures=() if extrapolation_plot is None else (extrapolation_plot,),
    )

    if indexing is not None and result.method != "le_bail":
        floor_arg = float(request["minimum_two_theta_deg"]) or None
        alternatives: list[tuple[str, LatticeParameterResult]] = []
        attempts: list[tuple[str, str, str]] = [
            ("Cohen, no systematic term", "cohen", "none"),
            ("Cohen, Nelson–Riley term", "cohen", "nelson_riley"),
        ]
        if system == "cubic":
            attempts.append(("Average over reflections", "average", "none"))
        for label, alt_method, alt_extrapolation in attempts:
            if alt_method == result.method and alt_extrapolation == result.extrapolation:
                continue
            try:
                alternatives.append(
                    (
                        label,
                        determine_lattice_parameters(
                            indexing,
                            phase,
                            method=cast(Any, alt_method),
                            extrapolation=cast(Any, alt_extrapolation),
                            minimum_two_theta_deg=floor_arg,
                        ),
                    )
                )
            except ValueError:
                continue
        if alternatives:
            reported = float(result.a)
            rows_cross: list[dict[str, Any]] = [
                {
                    "method": (
                        f"Reported: {_METHOD_LABELS[result.method]}, "
                        f"{EXTRAPOLATION_NAMES[result.extrapolation]}"
                    ),
                    "a_angstrom": reported,
                    "sigma_angstrom": float(result.a_standard_uncertainty),
                    "difference_ppm": 0.0,
                    "reduced_chi_squared": float(result.reduced_chi_squared),
                }
            ]
            rows_cross += [
                {
                    "method": label,
                    "a_angstrom": float(item.a),
                    "sigma_angstrom": float(item.a_standard_uncertainty),
                    "difference_ppm": 1.0e6 * (float(item.a) - reported) / reported,
                    "reduced_chi_squared": float(item.reduced_chi_squared),
                }
                for label, item in alternatives
            ]
            widest = max(abs(float(row["difference_ppm"])) for row in rows_cross)
            diagnostics.append(
                ResultStage(
                    key="cross_check",
                    title="Comparison with other methods",
                    summary=(
                        f"The same {indexing.indexed_count} assigned reflections through "
                        f"{len(alternatives)} other method"
                        f"{'s' if len(alternatives) != 1 else ''} give values of a up to "
                        f"{widest:.0f} ppm from the reported one."
                    ),
                    metrics=(ResultMetric("Largest difference", widest, "ppm"),),
                    table=ResultTable(columns=_STAGE_CROSS_CHECK_COLUMNS, rows=tuple(rows_cross)),
                    explanation=(
                        "These alternatives reuse this run's peaks and assignment, so every "
                        "difference is due to the method alone. A large gap between the fits "
                        "with and without a systematic term says the term did real work; a gap "
                        "within a few σ(a) says the specimen was well aligned. The average is "
                        "shown for a cubic cell as a teaching comparison: it cannot remove an "
                        "angle-dependent error, so its disagreement grows with the aberration. "
                        "Prefer the method whose lattice-fit χ²ν is closest to 1."
                    ),
                    status="info",
                    section="diagnostics",
                    figures=(cross_check_figure(rows_cross),),
                )
            )
    return (cell_stage, *evidence, *diagnostics, *method_stages, *audit)


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

    identified = AppResult(
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
        stages=_identification_stages(
            request=request,
            measured=measured,
            generated=generated,
            radiation=radiation,
            spec=spec,
            identification=identification,
            table=table,
        ),
    )
    return replace(identified, figures=identification_figures(identified.data)).to_json()


_CRITERION_LABELS = {
    "explained_intensity_fraction": "Intensity explained",
    "completeness": "Lines seen",
    "position_score": "Position",
    "intensity_agreement": "Intensity agreement",
}

_STAGE_CELL_SEARCH_COLUMNS = (
    Column("phase_name", "Candidate"),
    Column(
        "cell_dilation_percent",
        "Cell dilation",
        units="%",
        numeric=True,
        digits=3,
        help_text="Uniform scale applied to the candidate's cell before matching, minus one.",
    ),
    Column(
        "at_edge",
        "At the search limit",
        help_text="A dilation at the limit means the candidate had to be stretched to fit.",
    ),
    Column("indexed", "Peaks indexed"),
)

_STAGE_MATCH_COLUMNS = (
    Column("hkl_label", "Reflection"),
    Column("two_theta_observed_deg", "2θ observed", units="°", numeric=True, digits=4),
    Column("two_theta_calculated_deg", "2θ calculated", units="°", numeric=True, digits=4),
    Column("delta_mdeg", "Δ2θ", units="m°", numeric=True, digits=1),
    Column(
        "intensity_observed",
        "I observed",
        numeric=True,
        digits=3,
        help_text="Integrated intensity relative to the strongest indexed peak.",
    ),
    Column(
        "intensity_calculated",
        "I calculated",
        numeric=True,
        digits=3,
        help_text="Calculated relative intensity of the reflection.",
    ),
)

_STAGE_SCORE_COLUMNS = (
    Column("rank", "Rank", numeric=True),
    Column("phase_name", "Candidate"),
    *(
        Column(
            f"from_{name}",
            f"From {label.lower()}",
            numeric=True,
            digits=3,
            help_text=(
                "This criterion's weighted contribution to the score: weight × value, divided by "
                "the total weight of the criteria defined for this candidate."
            ),
        )
        for name, label in _CRITERION_LABELS.items()
    ),
    Column("score", "Score", numeric=True, digits=3),
)


def _identification_stages(
    *,
    request: Mapping[str, Any],
    measured: MeasuredPowderPattern,
    generated: bool,
    radiation: RadiationSpec,
    spec: Any,
    identification: Any,
    table: Any,
) -> tuple[ResultStage, ...]:
    """Report every step of a phase identification as its own stage.

    Detection, the per-candidate cell search, the assignment of the leading
    candidates, the criterion-by-criterion scores, and the decision against its
    thresholds. The ranking alone cannot say whether a candidate lost on its
    cell, its centring or its intensities; the stages can.
    """

    stages: list[ResultStage] = [
        _scan_stage(
            request=request,
            measured=measured,
            generated=generated,
            radiation=radiation,
            displacement=None,
        )
    ]
    peaks = list(table)
    prominence = float(request["prominence_sigma"])
    floor = float(request["minimum_two_theta_deg"])
    stages.append(
        ResultStage(
            key="peaks",
            title="2. Peak detection and profile fitting",
            summary=(
                f"{len(peaks)} converged peak fits were kept above {prominence:g} robust noise "
                "standard deviations"
                + (f" and above {floor:.2f}° 2θ" if floor > 0.0 else "")
                + ". Every candidate is judged against this one list, so a missed or spurious "
                "peak affects all of them alike."
            ),
            metrics=(
                ResultMetric("Peaks kept", len(peaks)),
                ResultMetric("Detection threshold", prominence, "σ"),
                ResultMetric("Angular floor", floor, "°"),
                ResultMetric(
                    "Total integrated intensity",
                    float(np.sum(table.integrated_intensity)) if peaks else 0.0,
                ),
            ),
            table=ResultTable(
                columns=_STAGE_PEAK_COLUMNS,
                rows=_peak_stage_rows(peaks),
                caption="The measured lines every candidate is scored against.",
            ),
            explanation=(
                "Peaks are found by a width-matched Ricker filter on the variance-stabilized "
                "scan and each is fitted with a pseudo-Voigt on a linear local background; only "
                "converged fits are kept. The integrated intensity of each peak is what the "
                "'intensity explained' criterion sums, so a background ripple promoted to a peak "
                "lowers every candidate's score, and a weak line lost below the threshold can "
                "hide a second phase."
            ),
            status="warning" if len(peaks) < 3 else "ok",
        )
    )

    candidates = list(identification.candidates)
    search = float(request["cell_scale_range"])
    cell_rows = []
    at_edge_count = 0
    for candidate in candidates:
        dilation = float(candidate.cell_scale) - 1.0
        at_edge = search > 0.0 and abs(dilation) >= 0.98 * search
        if at_edge and candidate.indexing is not None:
            at_edge_count += 1
        unindexed = 0 if candidate.indexing is None else len(candidate.indexing.unindexed_peaks)
        cell_rows.append(
            {
                "phase_name": candidate.phase_name,
                "cell_dilation_percent": 100.0 * dilation,
                "at_edge": "yes" if at_edge else "no",
                "indexed": (
                    candidate.rejection
                    if candidate.indexing is None
                    else f"{candidate.indexed_count} of {candidate.indexed_count + unindexed}"
                ),
            }
        )
    stages.append(
        ResultStage(
            key="cell_search",
            title="3. Cell dilation search",
            summary=(
                (
                    f"Each candidate's cell was scaled uniformly within ±{100.0 * search:.2f} % "
                    "to place its lines best before matching."
                    if search > 0.0
                    else "No cell dilation was searched: every candidate was matched at its "
                    "tabulated cell."
                )
                + (
                    f" {at_edge_count} indexed candidate"
                    f"{'s' if at_edge_count != 1 else ''} reached the limit of the search."
                    if at_edge_count
                    else ""
                )
            ),
            metrics=(ResultMetric("Search half-range", 100.0 * search, "%"),),
            table=ResultTable(columns=_STAGE_CELL_SEARCH_COLUMNS, rows=tuple(cell_rows)),
            explanation=(
                "A tabulated cell belongs to someone else's specimen; a real solid solution, "
                "temperature or residual stress shifts it by a fraction of a per cent, which "
                "displaces back-reflection lines by more than any sensible tolerance. The search "
                "tries a grid of uniform scales and keeps the one minimizing the total "
                "tolerance-clipped distance from each measured peak to its nearest calculated "
                "line. A uniform scale preserves every ratio of d spacings, so it cannot make a "
                "wrong structure fit; a dilation at the limit of the range, however, means the "
                "candidate was stretched as far as allowed and should be read with suspicion."
            ),
            status="warning" if at_edge_count else "ok",
        )
    )

    def match_stage(key: str, title: str, candidate: Any) -> ResultStage:
        indexing = candidate.indexing
        reflections = list(indexing.reflections)
        strongest = max(
            (float(item.peak.integrated_intensity) for item in reflections), default=0.0
        )
        unindexed = [float(peak.two_theta_deg) for peak in indexing.unindexed_peaks]
        merit_m, count_m = indexing.figure_of_merit_m()
        summary = (
            f"{candidate.phase_name} indexed {indexing.indexed_count} of "
            f"{indexing.indexed_count + len(unindexed)} peaks within "
            f"±{indexing.tolerance_deg:.3f}°, with a mean |Δ2θ| of "
            f"{1000.0 * indexing.mean_absolute_delta_two_theta_deg:.1f} m° and de Wolff "
            f"M_{count_m} = {merit_m:.1f}."
        )
        if unindexed:
            summary += (
                " Unindexed: " + ", ".join(f"{value:.3f}°" for value in unindexed) + "."
            )
        if candidate.strongest_unobserved_relative_intensity > 0.0:
            summary += (
                f" The strongest predicted line that was not observed has "
                f"{100.0 * candidate.strongest_unobserved_relative_intensity:.0f} % relative "
                "intensity."
            )
        return ResultStage(
            key=key,
            title=title,
            summary=summary,
            metrics=(
                ResultMetric("Peaks indexed", int(indexing.indexed_count)),
                ResultMetric("Peaks unindexed", len(unindexed)),
                ResultMetric(
                    "Strongest unexplained peak",
                    100.0 * float(candidate.strongest_unexplained_fraction),
                    "% of strongest",
                    "Near 100 % means the most prominent feature of the scan is unexplained.",
                ),
                ResultMetric(
                    "Strongest unobserved line",
                    100.0 * float(candidate.strongest_unobserved_relative_intensity),
                    "% relative",
                    "A strong predicted line that is absent argues against the candidate.",
                ),
                ResultMetric(f"de Wolff M_{count_m}", float(merit_m)),
            ),
            table=ResultTable(
                columns=_STAGE_MATCH_COLUMNS,
                rows=tuple(
                    {
                        "hkl_label": _powder_label(item.miller_indices, spec=spec),
                        "two_theta_observed_deg": float(item.peak.two_theta_deg),
                        "two_theta_calculated_deg": float(item.two_theta_calculated_deg),
                        "delta_mdeg": 1000.0 * float(item.delta_two_theta_deg),
                        "intensity_observed": (
                            float(item.peak.integrated_intensity) / strongest
                            if strongest > 0.0
                            else 0.0
                        ),
                        "intensity_calculated": float(item.relative_intensity_calculated),
                    }
                    for item in reflections
                ),
                caption=f"How {candidate.phase_name}'s dilated cell accounts for the peaks.",
            ),
            explanation=(
                "Peaks are assigned one-to-one to calculated reflections within the tolerance, "
                "after the cell dilation of the previous stage. Δ2θ that scatters about zero is "
                "a good match; a trend with angle means the cell or a zero error is still off. "
                "Compare the two intensity columns row by row: preferred orientation moves "
                "intensities without moving positions, so disagreement there is weaker evidence "
                "than a missing strong line."
            ),
            status="warning" if unindexed else "ok",
        )

    indexed = [candidate for candidate in candidates if candidate.indexing is not None]
    if indexed:
        stages.append(
            match_stage("best_match", "4. Assignment of the leading candidate", indexed[0])
        )
    if len(indexed) > 1:
        stages.append(match_stage("runner_up_match", "5. Assignment of the runner-up", indexed[1]))

    weights = dict(candidates[0].weights) if candidates else {}
    score_rows = []
    for position, candidate in enumerate(candidates, start=1):
        criteria = dict(candidate.criteria)
        defined = {
            name: value for name, value in criteria.items() if np.isfinite(value)
        }
        total_weight = sum(float(candidate.weights[name]) for name in defined)
        row: dict[str, Any] = {"rank": position, "phase_name": candidate.phase_name}
        for name in _CRITERION_LABELS:
            if candidate.indexing is None or name not in defined or total_weight <= 0.0:
                row[f"from_{name}"] = None
            else:
                row[f"from_{name}"] = (
                    float(candidate.weights[name]) * float(defined[name]) / total_weight
                )
        row["score"] = float(candidate.score)
        score_rows.append(row)
    number = len(stages) + 1
    weighting = str(request.get("weighting", "standard"))
    best = identification.best
    stages.append(
        ResultStage(
            key="scoring",
            title=f"{number}. Criterion scores",
            summary=(
                f"Scores combine four criteria with the '{weighting}' weighting ("
                + ", ".join(
                    f"{_CRITERION_LABELS[name].lower()} {100.0 * float(value):.0f} %"
                    for name, value in weights.items()
                )
                + f"). {best.phase_name} leads at {best.score:.3f}."
            ),
            metrics=tuple(
                ResultMetric(f"Weight: {_CRITERION_LABELS[name]}", float(value))
                for name, value in weights.items()
            ),
            table=ResultTable(
                columns=_STAGE_SCORE_COLUMNS,
                rows=tuple(score_rows),
                caption=(
                    "Each criterion's weighted contribution; the contributions sum to the score."
                ),
            ),
            explanation=(
                "Intensity explained is the share of the measured integrated intensity carried "
                "by the peaks the candidate indexed. Lines seen is the share of the candidate's "
                "own strong reflections inside the scan that were observed — the criterion that "
                "separates two cells differing only by a centring. Position is "
                "1 − mean|Δ2θ| / tolerance. Intensity agreement is one minus the Bray–Curtis "
                "dissimilarity of observed and calculated relative intensities. A criterion that "
                "cannot be measured for a candidate — intensity agreement with fewer than two "
                "indexed lines — is left out and the remaining weights renormalized, rather than "
                "counted as a failure. Read which column a losing candidate falls short in: "
                "position points to the wrong cell size, lines seen to the wrong centring or "
                "basis, intensity explained to a phase that is missing from the list."
            ),
        )
    )

    runner = identification.runner_up
    stages.append(
        ResultStage(
            key="decision",
            title=f"{number + 1}. Decision",
            summary=(
                f"The best score {best.score:.3f} "
                + (
                    "reaches"
                    if identification.is_conclusive
                    else "falls short of"
                )
                + f" the {identification.minimum_score:.2f} threshold"
                + (
                    ""
                    if runner is None
                    else f", and leads {runner.phase_name} by {identification.margin:.3f} "
                    f"against a required margin of {identification.decisive_margin:.3f}"
                )
                + ". The identification is "
                + (
                    "conclusive and decisive."
                    if identification.is_conclusive and identification.is_decisive
                    else "conclusive but not decisive."
                    if identification.is_conclusive
                    else "not conclusive."
                )
            ),
            metrics=(
                ResultMetric("Best score", float(best.score)),
                ResultMetric("Minimum score", float(identification.minimum_score)),
                ResultMetric(
                    "Runner-up score", None if runner is None else float(runner.score)
                ),
                ResultMetric("Margin", float(identification.margin)),
                ResultMetric("Required margin", float(identification.decisive_margin)),
                ResultMetric("Conclusive", bool(identification.is_conclusive)),
                ResultMetric("Decisive", bool(identification.is_decisive)),
            ),
            explanation=(
                "Conclusive means the best candidate explains the pattern well enough to be "
                "believed at all; decisive means, in addition, that no other candidate comes "
                "within the margin. A conclusive but indecisive result says the scan cannot "
                "separate the two leaders — extend the angular range, count longer, or add a "
                "technique that can. With a single candidate this is a check of that phase, not "
                "an identification, and a low best score is as likely to mean the right "
                "structure is missing from the list as that the scan is poor."
            ),
            status=(
                "ok"
                if identification.is_conclusive and identification.is_decisive
                else "warning"
            ),
        )
    )
    return tuple(stages)
