# ruff: noqa: RUF001
"""Powder X-ray diffraction for the shared web and desktop workbench.

The application layer does not implement diffraction physics.  It validates a
human-scale request, calls :func:`pytex.diffraction.xrd.generate_xrd_pattern`,
and turns the resulting reflection objects and sampled profile into the common
``AppResult`` contract used by tables, hover cards, and exports.
"""

from __future__ import annotations

from typing import Any, Literal, cast

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    ExampleScenario,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern

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
