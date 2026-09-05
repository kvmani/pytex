# ruff: noqa: RUF001
"""Powder X-ray diffraction for the shared web and desktop workbench.

The application layer does not implement diffraction physics.  It validates a
human-scale request, calls :func:`pytex.diffraction.xrd.generate_xrd_pattern`,
and turns the resulting reflection objects and sampled profile into the common
``AppResult`` contract used by tables, hover cards, and exports.
"""

from __future__ import annotations

from typing import Any, Literal, cast

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.diffraction.rietveld import _scaled_phase, refine_rietveld
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_background import estimate_background
from pytex.diffraction.xrd_instrument import (
    InstrumentBroadening,
    calibrate_instrument_broadening,
    deconvolve_instrument_width,
    scherrer_size_nm,
    williamson_hall,
)
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

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


def _scan_parameters(*, group: str = "Measurement") -> tuple[Any, ...]:
    """Return the shared "where does the scan come from" controls.

    One declaration, used by every analysis operation, so the meaning of a
    pasted scan cannot drift between views.
    """

    return (
        ChoiceParameter(
            name="data_source",
            label="Scan source",
            help_text=(
                "Where the measured profile comes from. A demonstration scan is generated from "
                "the selected phase with a deliberate cell dilation, zero-point error, peak "
                "width and curved background, then given Poisson counting noise \u2014 so every "
                "answer below can be checked against a value that is known in advance. It is "
                "marked synthetic in the result and is not a measurement of anything."
            ),
            options=(
                (
                    "demonstration",
                    "Generate a demonstration scan",
                    "Synthetic, with known answers, for learning and method development.",
                ),
                (
                    "paste",
                    "Use the pasted scan",
                    "Two columns of your own data: 2\u03b8 and intensity.",
                ),
            ),
            default="demonstration",
            group=group,
        ),
        TextParameter(
            name="scan",
            label="Pasted scan",
            help_text=(
                "Two numbers per line: 2\u03b8 in degrees, then intensity. Blank lines and "
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
            group=group,
        ),
        IntegerParameter(
            name="demonstration_seed",
            label="Demonstration noise seed",
            help_text=(
                "Seed for the counting noise of the demonstration scan, so a result is "
                "reproducible. Ignored when a pasted scan is used."
            ),
            default=20260905,
            minimum=0,
            maximum=2**31 - 1,
            advanced=True,
            group=group,
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
                hint="Each line needs a 2\u03b8 value and an intensity, separated by "
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

    if str(request["data_source"]) == "paste":
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
            hint="Choose a phase with an atomic basis, or paste a scan instead.",
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
