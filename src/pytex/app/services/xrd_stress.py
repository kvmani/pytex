# ruff: noqa: RUF001
"""Residual stress by the sin²ψ method, for the shared web and desktop workbench.

The application layer does not implement the science. It validates a
human-scale request, resolves the stress-free spacing and the diffraction
elastic constants, calls
:func:`pytex.diffraction.xrd_residual_stress.residual_stress_pipeline`, and
turns the result into the common ``AppResult`` contract: headline numbers,
warnings, a staged report with figures, and a per-measurement table. The
"Printable report" and "Report + figures" exports of that result are the downloadable
report.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any, cast

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.radiation import radiation_from_request, radiation_parameters
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    DocumentationLink,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
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
from pytex.app.services.xrd_stress_report import (
    COMPONENT_LABELS,
    budget_figure,
    dec_figure,
    geometry_figure,
    mohr_figure,
    normalized_strain_residuals,
    peak_fits_figure,
    peak_quality_figure,
    peak_shift_figure,
    sigma_phi_figure,
    sin2psi_figure,
    splitting_figure,
    strain_figure,
    strain_residual_figure,
    stress_correlation_figure,
    stress_highlights,
    stress_warnings,
)
from pytex.core.symbols import symbol_text
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_residual_stress import (
    SINGLE_CRYSTAL_STIFFNESS_GPA,
    DiffractionElasticConstants,
    ResidualStressResult,
    Sin2PsiMeasurement,
    determine_residual_stress,
    indices_of_orientations,
    parse_stress_peak_positions,
    parse_stress_scans,
    residual_stress_pipeline,
    simulate_sin2psi_measurement,
    single_crystal_stiffness,
    strain_design_matrix,
)
from pytex.properties.tensors import StiffnessTensor

__all__: tuple[str, ...] = ()

_CITATION_MACHERAUCH = "Macherauch & Müller, Z. angew. Phys. 13 (1961) 305."
_CITATION_NOYAN_COHEN = (
    "Noyan & Cohen, Residual Stress: Measurement by Diffraction and Interpretation, "
    "Springer (1987), doi:10.1007/978-1-4613-9570-6."
)
_CITATION_HAUK = (
    "Hauk (ed.), Structural and Residual Stress Analysis by Nondestructive Methods, "
    "Elsevier (1997), doi:10.1016/B978-0-444-82476-9.X5000-2."
)
_CITATION_DOELLE = "Dölle, J. Appl. Cryst. 12 (1979) 489, doi:10.1107/S0021889879013169."
_CITATION_WELZEL = (
    "Welzel, Ligot, Lamparter, Vermeulen & Mittemeijer, J. Appl. Cryst. 38 (2005) 1, "
    "doi:10.1107/S0021889804029516."
)
_CITATION_KROENER = "Kröner, Z. Physik 151 (1958) 504, doi:10.1007/BF01337948."
_CITATION_NPL = (
    "Fitzpatrick et al., Determination of Residual Stresses by X-ray Diffraction, "
    "NPL Measurement Good Practice Guide No. 52, issue 2 (2005)."
)
_CITATION_GUM = (
    "JCGM 100:2008, Evaluation of measurement data — Guide to the expression of "
    "uncertainty in measurement."
)
_CITATION_SIMMONS = (
    "Simmons & Wang, Single Crystal Elastic Constants and Calculated Aggregate Properties, "
    "MIT Press (1971)."
)


#: The built-in phases whose single-crystal stiffness is tabulated.
_STIFFNESS_OF_PHASE = {key: key for key in SINGLE_CRYSTAL_STIFFNESS_GPA}

_MATERIAL_LABELS = {
    "fe_bcc": "Ferrite (bcc Fe)",
    "ni_fcc": "Nickel",
    "al_fcc": "Aluminium",
    "cu_fcc": "Copper",
    "w_bcc": "Tungsten",
    "ti_hcp": "α-Titanium",
}

_MODEL_LABELS = {
    "kroener": "Kröner (self-consistent)",
    "reuss": "Reuss (uniform stress)",
    "voigt": "Voigt (uniform strain)",
    "hill": "Neerfeld–Hill (Reuss–Voigt mean)",
    "isotropic": "Isotropic (E, ν)",
    "user": "Entered directly (S₁, ½S₂)",
}

_STATE_LABELS = {
    "biaxial": "Biaxial (σ11, σ22, σ12)",
    "biaxial_shear": "Biaxial with out-of-plane shear (+ σ13, σ23)",
    "triaxial": "Triaxial (all six components)",
}

_METHOD_LABELS = {
    "pseudo_voigt": "Pseudo-Voigt profile fit (Kα1 + Kα2)",
    "parabola": "Parabola through the peak top",
    "centroid": "Centroid above half maximum",
    "given": "Supplied peak positions",
}

_DEMO_PHIS = "0, 45, 90"
_DEMO_PSIS = "-45, -37.8, -30, -21.1, 0, 21.1, 30, 37.8, 45"


def _psi() -> str:
    return symbol_text("stress_tilt")


def _phi() -> str:
    return symbol_text("stress_azimuth")


@REGISTRY.operation(
    "xrd.residual_stress",
    title="Residual stress (sin²ψ)",
    summary=(
        "The in-plane stress tensor from peak shifts at several tilts and azimuths, with every "
        "intermediate result, its uncertainty budget, and a downloadable report."
    ),
    help_text=(
        "A symmetric θ–2θ scan sees only the planes parallel to the surface, and a stress in "
        "the surface reaches them only through the Poisson contraction — small, and "
        "indistinguishable from a change of composition. Tilting the specimen by ψ brings "
        "inclined planes into reflection, and it is the *change* of their spacing with "
        "inclination that measures the stress.\n\n"
        "The fundamental equation: the strain along the scattering vector at azimuth φ and "
        "tilt ψ is ε = ½S₂·(m·σ·m) + S₁·tr σ, which under plane stress is "
        "ε = ½S₂·σφ·sin²ψ + S₁·(σ11 + σ22). So **d is linear in sin²ψ and the slope is "
        "d₀·½S₂·σφ**. Because the stress comes from a slope, an error in d₀ only rescales it; "
        "d₀ matters for the absolute strain, and so for σ11 + σ22 through S₁.\n\n"
        "Three azimuths not 180° apart give the three in-plane components. Tilts of both "
        "signs expose out-of-plane shear, which splits the ψ > 0 and ψ < 0 branches "
        "(ψ-splitting). The **diffraction elastic constants** S₁ and ½S₂ belong to the "
        "reflection: computed here from single-crystal stiffness under the Reuss, Voigt, "
        "Neerfeld–Hill or Kröner grain-interaction model, or from E and ν.\n\n"
        "Each peak is located by a Kα1/Kα2 pseudo-Voigt fit (or, as a cross-check, by a "
        "parabola or a centroid after Kα2 stripping) with a standard uncertainty. All positions "
        "then enter one weighted least-squares fit for the tensor, beside the classical line "
        "per azimuth. The uncertainty budget separates the statistical part from d₀ and the "
        "elastic constants, and a Monte Carlo propagation checks it.\n\n"
        "**Data.** The demonstration generates a measurement of a *known* stress, so every "
        "number can be checked. Measured data are pasted or opened as a table: either whole "
        "scans (four columns φ ψ 2θ intensity, one row per point) or peak positions already "
        "located elsewhere (φ ψ 2θ and optionally u(2θ))."
    ),
    parameters=(
        phase_parameter(
            help_text=(
                "The phase the stress is measured in. Its cell, with the stress-free lattice "
                "parameter below, gives d₀ of the reflection; its symmetry decides how the "
                "elastic constants are averaged."
            ),
            builtin="fe_bcc",
        ),
        IndicesParameter(
            name="reflection",
            label="Stress reflection",
            help_text=(
                "The (hkl) measured at every tilt. Choose a reflection at high 2θ — strain "
                "sensitivity grows as tanθ — such as ferrite (211) with Cr Kα at 156°."
            ),
            default=[2, 1, 1],
        ),
        *radiation_parameters(
            help_text=(
                "Cr Kα is conventional for ferritic steel, Co or Mn Kα for austenite, Cu Kα for "
                "nickel alloys at (420). The wavelength converts every angle to a spacing."
            ),
            default="cr_ka_doublet",
        ),
        NumberParameter(
            name="a0_angstrom",
            label="Stress-free lattice parameter a",
            help_text=(
                "a₀ of the unstressed material, from which d₀ of the reflection is computed. "
                "Leave empty to use the phase's tabulated cell. Measure it on stress-free "
                "material of the same composition: a d₀ error moves σ11 + σ22 through S₁."
            ),
            units="Å",
            required=False,
            minimum=0.5,
            maximum=50.0,
            symbol="a_zero",
            row="stress_free_cell",
            group="Stress-free reference",
        ),
        NumberParameter(
            name="c0_angstrom",
            label="Stress-free lattice parameter c",
            help_text=(
                "c₀ for a hexagonal, tetragonal or orthorhombic cell. Leave empty to scale the "
                "tabulated c with a₀. Ignored for cubic phases."
            ),
            units="Å",
            required=False,
            minimum=0.5,
            maximum=50.0,
            symbol="c_zero",
            row="stress_free_cell",
            group="Stress-free reference",
        ),
        NumberParameter(
            name="a0_uncertainty_angstrom",
            label="Standard uncertainty of the stress-free lattice parameter",
            help_text=(
                "u(a₀), carried into the budget as u(d₀)/d₀ = u(a₀)/a₀. A part in 10⁴ is a "
                "good laboratory d₀; composition differences are often worse."
            ),
            units="Å",
            default=0.0003,
            minimum=0.0,
            maximum=0.1,
            group="Stress-free reference",
        ),
        BooleanParameter(
            name="refine_d0",
            label="Determine d₀ from the data (plane stress)",
            help_text=(
                "Refine d₀ with the stress under the assumption σ33 = 0, instead of using the "
                "value above. Useful when no stress-free reference exists; not possible with "
                "the triaxial evaluation, where σ33 and d₀ are inseparable."
            ),
            default=False,
            group="Stress-free reference",
        ),
        ChoiceParameter(
            name="dec_model",
            label="Elastic-constant model",
            help_text=(
                "How S₁ and ½S₂ of the reflection are obtained. Kröner's self-consistent model "
                "is usually closest to measurement for untextured cubic metals; Reuss and Voigt "
                "bound it; the isotropic constants ignore the reflection entirely."
            ),
            options=tuple(
                (key, _MODEL_LABELS[key], description)
                for key, description in (
                    ("kroener", "Eshelby sphere in the self-consistent effective medium."),
                    ("reuss", "Every grain carries the macroscopic stress."),
                    ("voigt", "Every grain carries the macroscopic strain."),
                    ("hill", "The mean of the Reuss and Voigt constants."),
                    ("isotropic", "From Young's modulus and Poisson's ratio."),
                    ("user", "Enter S₁ and ½S₂ directly."),
                )
            ),
            default="kroener",
            group="Elastic constants",
        ),
        ChoiceParameter(
            name="stiffness_source",
            label="Single-crystal stiffness",
            help_text=(
                "The single-crystal constants the grain-interaction models average. "
                "'From the phase' uses the tabulated constants of the built-in phase chosen."
            ),
            options=(
                ("auto", "From the phase", "Tabulated for Fe, Ni, Al, Cu, W and α-Ti."),
                *(
                    (key, _MATERIAL_LABELS[key], "Tabulated at room temperature.")
                    for key in SINGLE_CRYSTAL_STIFFNESS_GPA
                ),
                ("custom", "Cubic, entered below", "C₁₁, C₁₂ and C₄₄ from the fields below."),
            ),
            default="auto",
            group="Elastic constants",
        ),
        NumberParameter(
            name="dec_uncertainty_percent",
            label="Relative uncertainty of the elastic constants",
            help_text=(
                "u(S)/S for both constants, carried into the budget. The spread between the "
                "Reuss and Voigt constants of the reflection is a fair guide; 5 % is typical."
            ),
            units="%",
            default=5.0,
            minimum=0.0,
            maximum=50.0,
            group="Elastic constants",
        ),
        NumberParameter(
            name="c11_gpa",
            label="Single-crystal stiffness C11",
            help_text="For the 'Cubic, entered below' stiffness.",
            units="GPa",
            default=231.4,
            minimum=1.0,
            maximum=2000.0,
            symbol="stiffness_c11",
            row="cubic_stiffness",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="c12_gpa",
            label="Single-crystal stiffness C12",
            help_text="For the 'Cubic, entered below' stiffness.",
            units="GPa",
            default=134.7,
            minimum=0.0,
            maximum=2000.0,
            symbol="stiffness_c12",
            row="cubic_stiffness",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="c44_gpa",
            label="Single-crystal stiffness C44",
            help_text="For the 'Cubic, entered below' stiffness.",
            units="GPa",
            default=116.4,
            minimum=1.0,
            maximum=2000.0,
            symbol="stiffness_c44",
            row="cubic_stiffness",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="youngs_modulus_gpa",
            label="Young's modulus",
            help_text="For the isotropic model.",
            units="GPa",
            default=210.0,
            minimum=1.0,
            maximum=2000.0,
            symbol="youngs_modulus",
            row="isotropic_constants",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="poisson_ratio",
            label="Poisson's ratio",
            help_text="For the isotropic model.",
            default=0.28,
            minimum=-0.99,
            maximum=0.49,
            symbol="poisson_ratio",
            row="isotropic_constants",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="s1_per_tpa",
            label="Diffraction elastic constant S1",
            help_text="For entering the constants directly, in 1/TPa (10⁻⁶/MPa).",
            units="TPa⁻¹",
            default=-1.25,
            minimum=-50.0,
            maximum=0.0,
            symbol="dec_s1",
            row="direct_constants",
            advanced=True,
            group="Elastic constants",
        ),
        NumberParameter(
            name="half_s2_per_tpa",
            label="Diffraction elastic constant half S2",
            help_text="For entering the constants directly, in 1/TPa (10⁻⁶/MPa).",
            units="TPa⁻¹",
            default=5.8,
            minimum=0.01,
            maximum=100.0,
            symbol="dec_half_s2",
            row="direct_constants",
            advanced=True,
            group="Elastic constants",
        ),
        ChoiceParameter(
            name="data_source",
            label="Data",
            help_text=(
                "A demonstration of a known stress, or your own measurement pasted (or opened "
                "with **Open a measurement file**) in the box below."
            ),
            options=(
                (
                    "demonstration",
                    "Generate a demonstration measurement",
                    "Synthetic scans of a known stress, with counting noise.",
                ),
                (
                    "scans",
                    "Measured scans (φ ψ 2θ intensity)",
                    "One row per measured point; the peaks are located here.",
                ),
                (
                    "positions",
                    "Measured peak positions (φ ψ 2θ [u])",
                    "Peaks already located elsewhere; optional fourth column u(2θ).",
                ),
            ),
            default="demonstration",
            group="Measurement",
        ),
        ChoiceParameter(
            name="geometry",
            label="Tilt geometry",
            help_text=(
                "ω-tilting tilts in the diffraction plane, so absorption changes across the "
                "profile and the peak broadens with tilt; χ-tilting tilts about the line of "
                "the diffraction plane and the surface, and does neither."
            ),
            options=(
                ("omega", "ω (iso-inclination)", "Tilt in the diffraction plane."),
                ("chi", "χ (side-inclination)", "Tilt perpendicular to it."),
            ),
            default="omega",
            group="Measurement",
        ),
        TextParameter(
            name="measurement",
            label="Measured data",
            help_text=(
                "Scans: four numbers per line, φ ψ 2θ intensity, one line per measured point "
                "(a scan is the lines sharing φ and ψ). Peak positions: φ ψ 2θ and optionally "
                "u(2θ), one line per tilt. Angles in degrees; whitespace, commas or "
                "semicolons separate; `#` starts a comment; a header line is skipped. Ignored "
                "for the demonstration."
            ),
            multiline=True,
            required=False,
            default="",
            placeholder="0  -45  155.20  412\n0  -45  155.25  436",
            advanced=True,
            group="Measurement",
        ),
        NumberParameter(
            name="true_sigma_11_mpa",
            label="Demonstration stress sigma 11",
            help_text="The stress the demonstration is generated with, along S1.",
            units="MPa",
            default=-350.0,
            minimum=-3000.0,
            maximum=3000.0,
            symbol="sigma_11",
            row="true_in_plane",
            group="Demonstration",
        ),
        NumberParameter(
            name="true_sigma_22_mpa",
            label="Demonstration stress sigma 22",
            help_text="The stress the demonstration is generated with, along S2.",
            units="MPa",
            default=-150.0,
            minimum=-3000.0,
            maximum=3000.0,
            symbol="sigma_22",
            row="true_in_plane",
            group="Demonstration",
        ),
        NumberParameter(
            name="true_sigma_12_mpa",
            label="Demonstration stress sigma 12",
            help_text="The in-plane shear the demonstration is generated with.",
            units="MPa",
            default=60.0,
            minimum=-3000.0,
            maximum=3000.0,
            symbol="sigma_12",
            row="true_in_plane",
            group="Demonstration",
        ),
        NumberParameter(
            name="true_sigma_13_mpa",
            label="Demonstration stress sigma 13",
            help_text=(
                "Out-of-plane shear in the demonstration. Non-zero values split the ψ > 0 and "
                "ψ < 0 branches."
            ),
            units="MPa",
            default=0.0,
            minimum=-1000.0,
            maximum=1000.0,
            symbol="sigma_13",
            row="true_shear",
            advanced=True,
            group="Demonstration",
        ),
        NumberParameter(
            name="true_sigma_23_mpa",
            label="Demonstration stress sigma 23",
            help_text="Out-of-plane shear in the demonstration.",
            units="MPa",
            default=0.0,
            minimum=-1000.0,
            maximum=1000.0,
            symbol="sigma_23",
            row="true_shear",
            advanced=True,
            group="Demonstration",
        ),
        TextParameter(
            name="demo_azimuths",
            label="Demonstration azimuths",
            help_text="Azimuths φ in degrees, comma separated.",
            default=_DEMO_PHIS,
            max_length=120,
            advanced=True,
            group="Demonstration",
        ),
        TextParameter(
            name="demo_tilts",
            label="Demonstration tilts",
            help_text=(
                "Tilts ψ in degrees, comma separated. The default steps equally in sin²ψ from "
                "0 to 0.5 on both sides."
            ),
            default=_DEMO_PSIS,
            max_length=200,
            advanced=True,
            group="Demonstration",
        ),
        NumberParameter(
            name="demo_fwhm_deg",
            label="Demonstration peak width",
            help_text="FWHM at ψ = 0; it grows as 1/cosψ under ω-tilting.",
            units="° 2θ",
            default=1.2,
            minimum=0.05,
            maximum=6.0,
            advanced=True,
            group="Demonstration",
        ),
        TextParameter(
            name="demo_bad_points",
            label="Demonstration bad measurements",
            help_text=(
                "φ ψ pairs (degrees, separated by semicolons) at which the demonstration's peak is "
                "displaced by 0.25° 2θ, as a specimen-height error at that one tilt would do: a "
                "known bad point to find in the plot and exclude."
            ),
            required=False,
            default="",
            max_length=400,
            placeholder="45 30",
            advanced=True,
            group="Demonstration",
        ),
        IntegerParameter(
            name="demo_seed",
            label="Demonstration noise seed",
            help_text="Seed of the counting noise, so a demonstration is reproducible.",
            default=20260924,
            minimum=0,
            maximum=2**31 - 1,
            advanced=True,
            group="Demonstration",
            field_width="short",
        ),
        ChoiceParameter(
            name="stress_state",
            label="Stress state",
            help_text=(
                "Biaxial assumes no stress on the free surface's normal (σi3 = 0), which holds "
                "within the few micrometres X-rays reach. Add the out-of-plane shear when the "
                "branches split; the triaxial evaluation also frees σ33 and needs an exact d₀."
            ),
            options=tuple(
                (key, _STATE_LABELS[key], description)
                for key, description in (
                    ("biaxial", "Needs three azimuths."),
                    ("biaxial_shear", "Needs tilts of both signs."),
                    ("triaxial", "Needs an exact d₀."),
                )
            ),
            default="biaxial",
            group="Evaluation",
        ),
        ChoiceParameter(
            name="peak_method",
            label="Peak location",
            help_text=(
                "The profile fit models the Kα doublet and is the default. The parabola and the "
                "centroid are classical cross-checks, applied after Kα2 stripping."
            ),
            options=(
                ("pseudo_voigt", _METHOD_LABELS["pseudo_voigt"], "Recommended."),
                ("parabola", _METHOD_LABELS["parabola"], "The top 20 % of the stripped peak."),
                ("centroid", _METHOD_LABELS["centroid"], "Continuous window at half maximum."),
            ),
            default="pseudo_voigt",
            group="Evaluation",
        ),
        TextParameter(
            name="excluded_points",
            label="Excluded measurements",
            help_text=(
                "Measurements left out of every fit, as φ ψ pairs in degrees separated by "
                "semicolons, e.g. `0 -45; 90 30`. Click a point in the d against sin²ψ plot to "
                "add or remove it, then press Refit. An excluded point stays in the plot and the "
                "report, with the strain the tensor fitted without it predicts, so the decision "
                "can be checked and reversed. Exclude a bad peak location; do not exclude a "
                "point only because it disagrees with a straight line — curvature is information."
            ),
            required=False,
            default="",
            max_length=4000,
            placeholder="0 -45; 90 30",
            group="Evaluation",
        ),
        BooleanParameter(
            name="lpa_correction",
            label="Correct for Lorentz-polarization and absorption",
            help_text=(
                "Divide each profile by the LPA factor before locating the peak. Broad peaks "
                "are skewed by it, differently at each tilt under ω-tilting."
            ),
            default=True,
            group="Evaluation",
        ),
        NumberParameter(
            name="window_deg",
            label="Fit window",
            help_text=(
                "Width of scan, centred on the expected peak, that the peak is located in. "
                "Wide enough for the peak and some background on both sides."
            ),
            units="° 2θ",
            default=8.0,
            minimum=0.5,
            maximum=40.0,
            advanced=True,
            group="Evaluation",
        ),
        NumberParameter(
            name="expected_fwhm_deg",
            label="Expected peak width",
            help_text="Starting width of the profile fit; right to within a factor of two.",
            units="° 2θ",
            default=1.2,
            minimum=0.02,
            maximum=8.0,
            advanced=True,
            group="Evaluation",
        ),
        IntegerParameter(
            name="monte_carlo_draws",
            label="Monte Carlo draws",
            help_text=(
                "Draws of every input from its uncertainty, refitted each time, as a check on "
                "the linear propagation. Zero skips it."
            ),
            default=2000,
            minimum=0,
            maximum=100000,
            advanced=True,
            group="Evaluation",
            field_width="short",
        ),
    ),
    returns=(
        "The stress tensor with its uncertainty budget, the per-azimuth d against sin²ψ lines, "
        "every located peak, and a staged report with publication figures."
    ),
    panel="xrd",
    citations=(
        _CITATION_MACHERAUCH,
        _CITATION_NOYAN_COHEN,
        _CITATION_HAUK,
        _CITATION_DOELLE,
        _CITATION_WELZEL,
        _CITATION_NPL,
        _CITATION_GUM,
    ),
    tags=(
        "XRD",
        "residual stress",
        "stress",
        "strain",
        "sin2psi",
        "sin²ψ",
        "X-ray elastic constants",
        "psi splitting",
        "Kröner",
    ),
    documentation=DocumentationLink(
        "Residual stress by the sin²ψ method", "algorithms/residual_stress_sin2psi"
    ),
)
def _residual_stress(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    radiation = radiation_from_request(request)
    hkl = tuple(int(value) for value in request["reflection"])
    if len(hkl) != 3:
        raise InvalidInputError(
            "The stress reflection needs three Miller indices.", field="reflection"
        )
    reflection = (hkl[0], hkl[1], hkl[2])
    label = plane_label(reflection, spec=spec)

    # -- the stress-free spacing -------------------------------------------
    lattice = phase.lattice
    a0 = request.get("a0_angstrom")
    c0 = request.get("c0_angstrom")
    scale_a = 1.0 if a0 is None else float(a0) / lattice.a
    scale_c = scale_a if c0 is None else float(c0) / lattice.c
    cubic = spec.crystal_system == "cubic"
    if cubic:
        scale_c = scale_a
    stress_free = replace(
        lattice, a=lattice.a * scale_a, b=lattice.b * scale_a, c=lattice.c * scale_c
    )
    reciprocal = stress_free.reciprocal_basis().matrix @ np.asarray(reflection, dtype=float)
    magnitude = float(np.linalg.norm(reciprocal))
    if magnitude == 0.0:
        raise InvalidInputError("The reflection (0 0 0) has no spacing.", field="reflection")
    d0 = 1.0 / magnitude
    tabulated_reciprocal = lattice.reciprocal_basis().matrix @ np.asarray(reflection, dtype=float)
    d0_tabulated = 1.0 / float(np.linalg.norm(tabulated_reciprocal))
    u_d0 = d0 * float(request["a0_uncertainty_angstrom"]) / stress_free.a
    if radiation.wavelength_angstrom >= 2.0 * d0:
        raise InvalidInputError(
            f"The {label} reflection (d₀ = {d0:.4f} Å) cannot diffract "
            f"{radiation.name} (λ = {radiation.wavelength_angstrom:.4f} Å): λ > 2d.",
            field="reflection",
            hint="Choose a lower-index reflection or a shorter wavelength.",
        )
    bragg = float(np.rad2deg(2.0 * np.arcsin(radiation.wavelength_angstrom / (2.0 * d0))))

    # -- the elastic constants ---------------------------------------------
    dec, stiffness, stiffness_note = _elastic_constants(request, spec, phase, reflection, label)

    # -- the measurement ---------------------------------------------------
    source = str(request["data_source"])
    geometry = cast(Any, request["geometry"])
    measurement: Sin2PsiMeasurement | None = None
    true_stress: np.ndarray | None = None
    stress_state = cast(Any, request["stress_state"])
    common = {
        "dec": dec,
        "stress_state": stress_state,
        "d0_uncertainty_angstrom": u_d0,
        "refine_d0": bool(request["refine_d0"]),
        "monte_carlo_draws": int(request["monte_carlo_draws"]),
        "seed": int(request["demo_seed"]),
        "reflection_label": label,
        "phase_name": spec.name,
    }
    pairs = _orientation_pairs(request.get("excluded_points"))
    try:
        if source == "positions":
            peaks = parse_stress_peak_positions(str(request.get("measurement") or ""))
            result = determine_residual_stress(
                peaks,
                wavelength_angstrom=radiation.wavelength_angstrom,
                d0_angstrom=d0,
                excluded=_excluded_indices(peaks, pairs),
                **common,
            )
        else:
            if source == "demonstration":
                true_stress = _demonstration_stress(request)
                measurement = simulate_sin2psi_measurement(
                    d0_angstrom=d0_tabulated,
                    stress_mpa=true_stress,
                    dec=dec,
                    phi_deg=_angles(request["demo_azimuths"], field="demo_azimuths"),
                    psi_deg=_angles(request["demo_tilts"], field="demo_tilts"),
                    radiation=radiation,
                    fwhm_deg=float(request["demo_fwhm_deg"]),
                    geometry=geometry,
                    apply_lpa=True,
                    seed=int(request["demo_seed"]),
                    name=f"demonstration: {spec.name} {label}",
                    corrupted=_orientation_pairs(
                        request.get("demo_bad_points"), field="demo_bad_points"
                    ),
                )
            else:
                measurement = parse_stress_scans(
                    str(request.get("measurement") or ""),
                    radiation=radiation,
                    geometry=geometry,
                    name="pasted measurement",
                )
            result = residual_stress_pipeline(
                measurement,
                d0_angstrom=d0,
                excluded=_excluded_indices(measurement.scans, pairs),
                expected_two_theta_deg=bragg,
                window_deg=float(request["window_deg"]),
                peak_method=cast(Any, request["peak_method"]),
                expected_fwhm_deg=float(request["expected_fwhm_deg"]),
                lpa_correction=bool(request["lpa_correction"]),
                **common,
            )
    except ValueError as error:
        message = str(error)
        field = "measurement" if source != "demonstration" else "demo_tilts"
        hint = (
            "Check the table: angles in degrees, one row per point (scans) or per tilt "
            "(positions), at least three tilts."
        )
        if "below the" in message:
            field, hint = "geometry", "Use χ-tilting or smaller tilts for this reflection."
        elif "triaxial" in message:
            field, hint = "refine_d0", "Switch off the d₀ refinement for a triaxial evaluation."
        elif "remain after the exclusions" in message:
            field, hint = "excluded_points", "Include some of the excluded measurements again."
        elif "cover the reflection" in message:
            field, hint = (
                "window_deg",
                f"The scans must include the reflection near 2θ = {bragg:.2f}°; check the "
                "reflection, the radiation and the stress-free lattice parameter.",
            )
        raise InvalidInputError(
            f"The residual stress could not be determined: {message}", field=field, hint=hint
        ) from error

    return _build_result(
        request=request,
        spec=spec,
        label=label,
        radiation=radiation,
        bragg=bragg,
        d0_tabulated=d0_tabulated,
        result=result,
        measurement=measurement,
        true_stress=true_stress,
        stiffness=stiffness,
        stiffness_note=stiffness_note,
        reflection=reflection,
        cubic=cubic,
        phase=phase,
    ).to_json()


def _orientation_pairs(text: Any, *, field: str = "excluded_points") -> list[tuple[float, float]]:
    """Read ``φ ψ; φ ψ`` pairs; commas may separate the two angles of a pair."""

    pairs: list[tuple[float, float]] = []
    for chunk in str(text or "").replace("\n", ";").split(";"):
        fields = chunk.replace(",", " ").split()
        if not fields:
            continue
        if len(fields) != 2:
            raise InvalidInputError(
                f"{chunk.strip()!r} is not a φ ψ pair.",
                field=field,
                hint="Two angles in degrees per measurement, pairs separated by semicolons.",
            )
        try:
            pairs.append((float(fields[0]), float(fields[1])))
        except ValueError as error:
            raise InvalidInputError(f"{chunk.strip()!r} is not a φ ψ pair.", field=field) from error
    return pairs


def _excluded_indices(items: Any, pairs: list[tuple[float, float]]) -> tuple[int, ...]:
    if not pairs:
        return ()
    try:
        return indices_of_orientations(items, pairs, tolerance_deg=0.01)
    except ValueError as error:
        raise InvalidInputError(
            str(error),
            field="excluded_points",
            hint="Each pair must name a measured azimuth and tilt; remove the pair or fix it.",
        ) from error


def _angles(text: Any, *, field: str) -> list[float]:
    values: list[float] = []
    for item in str(text).replace(";", ",").replace(" ", ",").split(","):
        if not item.strip():
            continue
        try:
            values.append(float(item))
        except ValueError as error:
            raise InvalidInputError(
                f"{item!r} is not an angle.", field=field, hint="Degrees, comma separated."
            ) from error
    if not values:
        raise InvalidInputError("No angles were given.", field=field)
    return values


def _demonstration_stress(request: dict[str, Any]) -> np.ndarray:
    s11 = float(request["true_sigma_11_mpa"])
    s22 = float(request["true_sigma_22_mpa"])
    s12 = float(request["true_sigma_12_mpa"])
    s13 = float(request["true_sigma_13_mpa"])
    s23 = float(request["true_sigma_23_mpa"])
    return np.array([[s11, s12, s13], [s12, s22, s23], [s13, s23, 0.0]])


def _elastic_constants(
    request: dict[str, Any],
    spec: Any,
    phase: Any,
    reflection: tuple[int, int, int],
    label: str,
) -> tuple[DiffractionElasticConstants, StiffnessTensor | None, str]:
    model = str(request["dec_model"])
    relative = float(request["dec_uncertainty_percent"]) / 100.0
    if model == "isotropic":
        dec = DiffractionElasticConstants.isotropic(
            float(request["youngs_modulus_gpa"]),
            float(request["poisson_ratio"]),
            hkl=reflection,
            relative_standard_uncertainty=relative,
        )
        return dec, None, "isotropic constants; no single-crystal stiffness used"
    if model == "user":
        try:
            dec = DiffractionElasticConstants(
                s1_per_tpa=float(request["s1_per_tpa"]),
                half_s2_per_tpa=float(request["half_s2_per_tpa"]),
                model="user",
                hkl=reflection,
                source="entered directly",
                relative_standard_uncertainty=relative,
            )
        except ValueError as error:
            raise InvalidInputError(str(error), field="s1_per_tpa") from error
        return dec, None, "constants entered directly"
    source = str(request["stiffness_source"])
    if source == "custom":
        if spec.crystal_system != "cubic":
            raise InvalidInputError(
                "Cubic single-crystal constants cannot describe this non-cubic phase.",
                field="stiffness_source",
                hint="Choose a tabulated material of the same symmetry, or the isotropic model.",
            )
        stiffness = StiffnessTensor.cubic(
            float(request["c11_gpa"]), float(request["c12_gpa"]), float(request["c44_gpa"])
        )
        note = (
            f"cubic constants entered: C11 = {float(request['c11_gpa']):g}, "
            f"C12 = {float(request['c12_gpa']):g}, C44 = {float(request['c44_gpa']):g} GPa"
        )
    else:
        key = source
        if source == "auto":
            builtin = (request.get("phase") or {}).get("builtin")
            key = _STIFFNESS_OF_PHASE.get(str(builtin), "")
            if not key:
                raise InvalidInputError(
                    f"No single-crystal stiffness is tabulated for {spec.name}.",
                    field="stiffness_source",
                    hint=(
                        "Choose a tabulated material, enter cubic constants, or use the "
                        "isotropic model with E and ν."
                    ),
                )
        stiffness = single_crystal_stiffness(key)
        system, constants = SINGLE_CRYSTAL_STIFFNESS_GPA[key]
        names = ("C11", "C12", "C44") if system == "cubic" else ("C11", "C12", "C13", "C33", "C44")
        note = (
            f"{_MATERIAL_LABELS[key]}: "
            + ", ".join(f"{name} = {value:g}" for name, value in zip(names, constants, strict=True))
            + " GPa"
        )
    dec = DiffractionElasticConstants.for_reflection(
        phase,
        reflection,
        stiffness,
        model=cast(Any, model),
        source=note,
        relative_standard_uncertainty=relative,
    )
    return dec, stiffness, note


# ---------------------------------------------------------------------------
# The result
# ---------------------------------------------------------------------------

_POINT_COLUMNS = (
    Column("phi_deg", "φ", units="°", numeric=True, digits=1, help_text="Azimuth from S1."),
    Column("psi_deg", "ψ", units="°", numeric=True, digits=2, help_text="Tilt from S3."),
    Column("sin2psi", "sin²ψ", numeric=True, digits=4),
    Column(
        "two_theta_deg", "2θ", units="°", numeric=True, digits=4, help_text="Located Kα1 position."
    ),
    Column("two_theta_uncertainty_mdeg", "u(2θ)", units="m°", numeric=True, digits=2),
    Column("d_angstrom", "d", units="Å", numeric=True, digits=6),
    Column("d_uncertainty_angstrom", "u(d)", units="Å", numeric=True, digits=6),
    Column(
        "strain_micro",
        "ε",
        units="10⁻⁶",
        numeric=True,
        digits=1,
        help_text="(d − d₀)/d₀ in microstrain.",
    ),
    Column("strain_uncertainty_micro", "u(ε)", units="10⁻⁶", numeric=True, digits=1),
    Column("fitted_strain_micro", "ε from the tensor", units="10⁻⁶", numeric=True, digits=1),
    Column("used", "Used", help_text="Whether the measurement entered the fits."),
    Column(
        "normalized_residual",
        "Residual / u",
        numeric=True,
        digits=2,
        help_text="(ε − ε from the tensor)/u(ε); within ±2 for about 95 % of rows.",
    ),
)

_SCAN_COLUMNS = (
    Column("phi_deg", "φ", units="°", numeric=True, digits=1),
    Column("psi_deg", "ψ", units="°", numeric=True, digits=2),
    Column("sin2psi", "sin²ψ", numeric=True, digits=4),
    Column("points", "Points", numeric=True, digits=0),
    Column("two_theta_min_deg", "2θ from", units="°", numeric=True, digits=2),
    Column("two_theta_max_deg", "2θ to", units="°", numeric=True, digits=2),
    Column("maximum_counts", "Maximum counts", numeric=True, digits=0),
)

_PEAK_COLUMNS = (
    Column("phi_deg", "φ", units="°", numeric=True, digits=1),
    Column("psi_deg", "ψ", units="°", numeric=True, digits=2),
    Column("two_theta_deg", "2θ", units="°", numeric=True, digits=4),
    Column("two_theta_uncertainty_mdeg", "u(2θ)", units="m°", numeric=True, digits=2),
    Column("fwhm_deg", "FWHM", units="°", numeric=True, digits=3),
    Column("height", "Height", numeric=True, digits=0),
    Column("reduced_chi_squared", "Peak-fit χ²ν", numeric=True, digits=2),
    Column("status", "Status"),
)

_LINE_COLUMNS = (
    Column("phi_deg", "φ", units="°", numeric=True, digits=1),
    Column("tilts", "Tilts", numeric=True, digits=0),
    Column("intercept_angstrom", "d at sin²ψ = 0", units="Å", numeric=True, digits=6),
    Column(
        "slope_mangstrom",
        "Slope",
        units="mÅ",
        numeric=True,
        digits=4,
        help_text="d(d)/d(sin²ψ) = d₀·½S₂·σφ.",
    ),
    Column("slope_uncertainty_mangstrom", "u(slope)", units="mÅ", numeric=True, digits=4),
    Column("sigma_phi_mpa", "σφ", units="MPa", numeric=True, digits=1),
    Column("sigma_phi_uncertainty_mpa", "u(σφ)", units="MPa", numeric=True, digits=1),
    Column(
        "tau_phi_mpa",
        "τφ",
        units="MPa",
        numeric=True,
        digits=1,
        help_text="From the ψ-splitting term; empty without tilts of both signs.",
    ),
    Column("tau_phi_uncertainty_mpa", "u(τφ)", units="MPa", numeric=True, digits=1),
    Column(
        "curvature_t",
        "Curvature / u",
        numeric=True,
        digits=2,
        help_text="The sin⁴ψ coefficient in units of its uncertainty; beyond ±3 is curved.",
    ),
    Column("reduced_chi_squared", "Line χ²ν", numeric=True, digits=2),
    Column("r_squared", "R²", numeric=True, digits=5),
)

_EXCLUDED_COLUMNS = (
    Column("phi_deg", "φ", units="°", numeric=True, digits=1),
    Column("psi_deg", "ψ", units="°", numeric=True, digits=2),
    Column("two_theta_deg", "2θ", units="°", numeric=True, digits=4),
    Column("normalized_residual", "Deleted residual / u", numeric=True, digits=2),
    Column("state", "State"),
)

_TENSOR_COLUMNS = (
    Column("component", "Component"),
    Column("value_mpa", "Value", units="MPa", numeric=True, digits=1),
    Column("combined_mpa", "u combined", units="MPa", numeric=True, digits=1),
    Column("statistical_mpa", "u statistical", units="MPa", numeric=True, digits=1),
    Column("d0_mpa", "u from d₀", units="MPa", numeric=True, digits=1),
    Column("dec_mpa", "u from S₁, ½S₂", units="MPa", numeric=True, digits=1),
    Column("monte_carlo_mpa", "u Monte Carlo", units="MPa", numeric=True, digits=1),
    Column(
        "true_mpa",
        "Generated with",
        units="MPa",
        numeric=True,
        digits=1,
        help_text="The demonstration's known stress; empty for measured data.",
    ),
)


def _number(value: float) -> float | None:
    return float(value) if math.isfinite(value) else None


def _build_result(
    *,
    request: dict[str, Any],
    spec: Any,
    label: str,
    radiation: RadiationSpec,
    bragg: float,
    d0_tabulated: float,
    result: ResidualStressResult,
    measurement: Sin2PsiMeasurement | None,
    true_stress: np.ndarray | None,
    stiffness: StiffnessTensor | None,
    stiffness_note: str,
    reflection: tuple[int, int, int],
    cubic: bool,
    phase: Any,
) -> AppResult:
    fit = result.tensor
    normalized = normalized_strain_residuals(result)
    rows = []
    for index, peak in enumerate(result.peaks):
        rows.append(
            {
                "phi_deg": peak.phi_deg,
                "psi_deg": peak.psi_deg,
                "sin2psi": peak.sin2psi,
                "two_theta_deg": peak.two_theta_deg,
                "two_theta_uncertainty_mdeg": 1000.0 * peak.two_theta_uncertainty_deg,
                "d_angstrom": float(result.d_angstrom[index]),
                "d_uncertainty_angstrom": float(result.d_uncertainty_angstrom[index]),
                "strain_micro": 1e6 * float(result.strain[index]),
                "strain_uncertainty_micro": 1e6 * float(result.strain_uncertainty[index]),
                "fitted_strain_micro": (
                    None if fit is None else 1e6 * float(fit.fitted_strain[index])
                ),
                "normalized_residual": (None if fit is None else _number(float(normalized[index]))),
                "used": "yes" if bool(result.included_mask[index]) else "excluded",
            }
        )

    warnings = stress_warnings(
        result, measurement=measurement, stress_state=str(request["stress_state"])
    )
    highlights = stress_highlights(result)

    if fit is not None:
        principal = fit.in_plane_principal()
        parts = ", ".join(
            f"{COMPONENT_LABELS[name]} = {value:.0f} ± {u:.0f} MPa"
            for name, value, u in zip(
                fit.component_names, fit.values_mpa, fit.combined_uncertainty_mpa, strict=True
            )
        )
        answer = (
            f"{parts}; in-plane principal stresses σI = {principal['sigma_I_mpa']:.0f} and "
            f"σII = {principal['sigma_II_mpa']:.0f} MPa, σI at {principal['angle_deg']:.0f}° "
            f"from S1"
        )
    else:
        answer = "; ".join(
            f"σ{_phi()} = {line.sigma_phi_mpa:.0f} ± {line.sigma_phi_uncertainty_mpa:.0f} MPa at "
            f"{_phi()} = {line.phi_deg:g}°"
            for line in result.regressions
        )
    summary = (
        f"{answer}. From the {label} reflection of {spec.name} with {radiation.name} at "
        f"2θ ≈ {bragg:.1f}°: {len(result.peaks)} peak positions at "
        f"{len(result.regressions)} azimuth(s), d₀ = {result.d0_angstrom:.6f} Å"
        + (" (refined under plane stress)" if result.d0_refined else "")
        + (
            f"; {len(result.excluded_indices)} of {len(result.peaks)} measurements excluded "
            "by the analyst"
            if result.excluded_indices
            else ""
        )
        + f", elastic constants from {_MODEL_LABELS[result.dec.model].split(' (')[0]} "
        f"(½S₂ = {result.dec.half_s2_per_tpa:.3f} TPa⁻¹). The ± values are combined standard "
        "uncertainties; tensile stress is positive."
        + (
            f" {len(warnings)} check{'s' if len(warnings) != 1 else ''} below "
            f"{'need' if len(warnings) != 1 else 'needs'} attention."
            if warnings
            else " No reliability check failed."
        )
    )

    stages = _stages(
        request=request,
        spec=spec,
        label=label,
        radiation=radiation,
        bragg=bragg,
        d0_tabulated=d0_tabulated,
        result=result,
        measurement=measurement,
        true_stress=true_stress,
        stiffness=stiffness,
        stiffness_note=stiffness_note,
        reflection=reflection,
        cubic=cubic,
        phase=phase,
    )

    notes: list[str] = []
    if measurement is not None and measurement.synthetic and true_stress is not None:
        notes.append(
            "The scans were generated, not measured: the reflection of the tabulated cell "
            f"(d₀ = {d0_tabulated:.6f} Å) strained by σ11 = {true_stress[0, 0]:g}, σ22 = "
            f"{true_stress[1, 1]:g}, σ12 = {true_stress[0, 1]:g}, σ13 = {true_stress[0, 2]:g}, "
            f"σ23 = {true_stress[1, 2]:g} MPa through the same elastic constants, as a Kα1/Kα2 "
            "pseudo-Voigt doublet broadened as 1/cosψ, shaped by the LPA factor, with Poisson "
            "counting noise. Enter a different stress-free lattice parameter to see what a "
            "wrong d₀ does: the slopes barely move, σ11 + σ22 does."
        )
    notes.append(
        "This is the macroscopic stress averaged over the depth the X-rays reach (a few "
        "micrometres in steel with Cr Kα), for an untextured material. Curved or oscillating "
        "d against sin²ψ would indicate a stress gradient or texture, which this evaluation "
        "does not model."
    )
    notes.append(
        "Download the full report as Printable report (one web page with every figure, to "
        "read or print to PDF) or as Report + figures (the same report in Markdown and HTML, "
        "each figure as a separate SVG, and the complete result as JSON)."
    )

    data = _plot_data(result)
    data.update(
        {
            "describe": result.describe(),
            "library": result.to_json(),
            "reflection": label,
            "phase_name": spec.name,
            "bragg_two_theta_deg": bragg,
            "synthetic": bool(measurement is not None and measurement.synthetic),
            "warnings": list(warnings),
        }
    )
    return AppResult(
        title=f"Residual stress in {spec.name} from {label}",
        summary=summary,
        highlights=highlights,
        warnings=warnings,
        table=ResultTable(
            columns=_POINT_COLUMNS,
            rows=tuple(rows),
            caption=(
                "Every measurement: the located peak, its spacing and strain with standard "
                "uncertainties, the strain the fitted tensor predicts, and the residual in units "
                "of u."
            ),
        ),
        data=data,
        inputs={key: value for key, value in request.items() if key != "measurement"}
        | {"measurement_lines": len(str(request.get("measurement") or "").splitlines())},
        notes=tuple(notes),
        citations=(
            _CITATION_MACHERAUCH,
            _CITATION_NOYAN_COHEN,
            _CITATION_HAUK,
            _CITATION_DOELLE,
            _CITATION_WELZEL,
            _CITATION_KROENER,
            _CITATION_NPL,
            _CITATION_GUM,
            _CITATION_SIMMONS,
        ),
        stages=stages,
    )


def _plot_data(result: ResidualStressResult) -> dict[str, Any]:
    """What the workbench view draws: every measurement and the lines, per azimuth.

    Every point is sent, the excluded ones too, each with its ``included`` flag
    and its strain residual in units of u, so the view can draw what was left
    out, flag what looks like an outlier, and let the analyst change the
    selection by clicking.
    """

    normalized = normalized_strain_residuals(result)
    lines = {round(line.phi_deg % 360.0, 6): line for line in result.regressions}
    azimuths = np.round(result.phi_deg % 360.0, 6)
    series = []
    for azimuth in sorted(set(azimuths.tolist())):
        members = np.flatnonzero(azimuths == azimuth)
        order = members[np.argsort(result.psi_deg[members])]
        entry: dict[str, Any] = {
            "phi_deg": float(azimuth),
            "psi_deg": result.psi_deg[order].tolist(),
            "measured_phi_deg": result.phi_deg[order].tolist(),
            "sin2psi": result.sin2psi[order].tolist(),
            "d_angstrom": result.d_angstrom[order].tolist(),
            "d_uncertainty_angstrom": result.d_uncertainty_angstrom[order].tolist(),
            "included": result.included_mask[order].tolist(),
            "deleted_residual": (
                [None] * int(order.size)
                if result.deleted_residuals is None
                else [_number(float(value)) for value in result.deleted_residuals[order]]
            ),
            "normalized_residual": (
                [_number(float(value)) for value in normalized[order]]
                if normalized.size
                else [None] * int(order.size)
            ),
        }
        line = lines.get(float(azimuth))
        top = float(np.max(result.sin2psi[order]))
        psi_grid = np.rad2deg(np.arcsin(np.sqrt(np.linspace(0.0, top, 40))))
        entry["grid_sin2psi"] = (np.sin(np.deg2rad(psi_grid)) ** 2).tolist()
        if line is not None:
            entry.update(
                {
                    "intercept_angstrom": line.intercept_angstrom,
                    "slope_angstrom": line.slope_angstrom,
                    "sigma_phi_mpa": line.sigma_phi_mpa,
                    "sigma_phi_uncertainty_mpa": line.sigma_phi_uncertainty_mpa,
                }
            )
        if result.tensor is not None:
            design = strain_design_matrix(
                np.full_like(psi_grid, azimuth),
                psi_grid,
                s1_per_tpa=result.dec.s1_per_tpa,
                half_s2_per_tpa=result.dec.half_s2_per_tpa,
                components=result.tensor.component_names,
            )
            entry["tensor_d_angstrom"] = (
                result.d0_angstrom * (1.0 + design @ result.tensor.values_mpa)
            ).tolist()
        series.append(entry)
    payload: dict[str, Any] = {
        "plot_kind": "sin2psi",
        "series": series,
        "d0_angstrom": result.d0_angstrom,
    }
    if result.tensor is not None:
        payload["tensor"] = {
            name: {"value_mpa": float(value), "uncertainty_mpa": float(u)}
            for name, value, u in zip(
                result.tensor.component_names,
                result.tensor.values_mpa,
                result.tensor.combined_uncertainty_mpa,
                strict=True,
            )
        }
        payload["reduced_chi_squared"] = result.tensor.reduced_chi_squared
    payload["excluded"] = [
        {"phi_deg": result.peaks[index].phi_deg, "psi_deg": result.peaks[index].psi_deg}
        for index in result.excluded_indices
    ]
    payload["suggested_outliers"] = [
        {"phi_deg": result.peaks[index].phi_deg, "psi_deg": result.peaks[index].psi_deg}
        for index in result.suggested_outliers()
    ]
    return payload


def _stages(
    *,
    request: dict[str, Any],
    spec: Any,
    label: str,
    radiation: RadiationSpec,
    bragg: float,
    d0_tabulated: float,
    result: ResidualStressResult,
    measurement: Sin2PsiMeasurement | None,
    true_stress: np.ndarray | None,
    stiffness: StiffnessTensor | None,
    stiffness_note: str,
    reflection: tuple[int, int, int],
    cubic: bool,
    phase: Any,
) -> tuple[ResultStage, ...]:
    fit = result.tensor
    stages: list[ResultStage] = []

    # -- Result: the tensor ------------------------------------------------
    if fit is not None:
        index = {
            "sigma_11": (0, 0),
            "sigma_22": (1, 1),
            "sigma_12": (0, 1),
            "sigma_13": (0, 2),
            "sigma_23": (1, 2),
            "sigma_33": (2, 2),
        }
        tensor_rows = []
        for position, name in enumerate(fit.component_names):
            i, j = index[name]
            tensor_rows.append(
                {
                    "component": COMPONENT_LABELS[name],
                    "value_mpa": float(fit.values_mpa[position]),
                    "combined_mpa": float(fit.combined_uncertainty_mpa[position]),
                    "statistical_mpa": float(fit.budget_mpa["statistical"][position]),
                    "d0_mpa": (
                        float(fit.budget_mpa["d0"][position]) if "d0" in fit.budget_mpa else None
                    ),
                    "dec_mpa": (
                        float(fit.budget_mpa["elastic constants"][position])
                        if "elastic constants" in fit.budget_mpa
                        else None
                    ),
                    "monte_carlo_mpa": (
                        None
                        if fit.monte_carlo_uncertainty_mpa is None
                        else float(fit.monte_carlo_uncertainty_mpa[position])
                    ),
                    "true_mpa": None if true_stress is None else float(true_stress[i, j]),
                }
            )
        principal = fit.in_plane_principal()
        equivalent, equivalent_u = fit.von_mises()
        figures = [figure for figure in (sigma_phi_figure(result), mohr_figure(result)) if figure]
        stages.append(
            ResultStage(
                key="stress_tensor",
                title="The stress tensor",
                section="result",
                summary=(
                    f"{_STATE_LABELS[fit.stress_state]} evaluation of all {len(result.peaks)} "
                    f"positions in one weighted least-squares fit: strain-fit χ²ν = "
                    f"{fit.reduced_chi_squared:.2f} on {fit.degrees_of_freedom} degrees of "
                    "freedom."
                    + (
                        " The demonstration's generating stress is listed beside each "
                        "component, so the answer can be judged against the truth."
                        if true_stress is not None
                        else ""
                    )
                ),
                metrics=(
                    ResultMetric("σI", round(principal["sigma_I_mpa"], 2), units="MPa"),
                    ResultMetric("σII", round(principal["sigma_II_mpa"], 2), units="MPa"),
                    ResultMetric(
                        "Direction of σI from S1", round(principal["angle_deg"], 2), units="°"
                    ),
                    ResultMetric("von Mises equivalent", round(equivalent, 2), units="MPa"),
                    ResultMetric("u(von Mises)", round(equivalent_u, 2), units="MPa"),
                ),
                table=ResultTable(
                    columns=_TENSOR_COLUMNS,
                    rows=tuple(tensor_rows),
                    caption="Each free component with its standard uncertainty, source by source.",
                ),
                explanation=(
                    "The components are in the specimen frame S1 S2 S3 (S3 the surface normal), "
                    "tensile positive. The principal stresses are the eigenvalues of the in-plane "
                    "2 × 2 block, their direction ½·atan2(2σ12, σ11 − σ22) from S1; their "
                    "uncertainties are the combined covariance propagated through those "
                    "functions. The von Mises equivalent is √(3/2 s:s) of the deviator."
                ),
                figures=tuple(figures),
            )
        )
        budget = budget_figure(result)
        stages.append(
            ResultStage(
                key="uncertainty_budget",
                title="Uncertainty budget",
                section="result",
                summary=_budget_summary(result),
                metrics=tuple(
                    ResultMetric(
                        f"Largest source for {COMPONENT_LABELS[name]}",
                        _largest_source(fit.budget_mpa, k),
                    )
                    for k, name in enumerate(fit.component_names)
                ),
                explanation=(
                    "Statistical: (AᵀWA)⁻¹ from the peak-position uncertainties, multiplied by "
                    "χ²ν when that exceeds one (the Birge ratio), so unexplained scatter widens "
                    "it. d₀: the exact sensitivity of the linear solution to d₀ times u(d₀). "
                    "Elastic constants: the sensitivities to S₁ and ½S₂ times their standard "
                    "uncertainties. The sources are independent and add in quadrature (GUM). The "
                    "Monte Carlo column redraws every input from its uncertainty and refits; "
                    "agreement with the combined column confirms that the linearization holds."
                ),
                figures=(budget,) if budget is not None else (),
            )
        )
    else:
        stages.append(
            ResultStage(
                key="stress_tensor",
                title="The stress tensor",
                section="result",
                status="warning",
                summary=(
                    "No tensor was determined: "
                    + (result.tensor_unavailable_reason or "the data do not constrain it.")
                    + " The stress along each measured azimuth is still read from its slope."
                ),
            )
        )

    # -- Evidence: the measurement ------------------------------------------
    if measurement is not None:
        scan_rows = tuple(
            {
                "phi_deg": scan.phi_deg,
                "psi_deg": scan.psi_deg,
                "sin2psi": float(np.sin(np.deg2rad(scan.psi_deg)) ** 2),
                "points": len(scan.pattern),
                "two_theta_min_deg": float(scan.pattern.two_theta_deg[0]),
                "two_theta_max_deg": float(scan.pattern.two_theta_deg[-1]),
                "maximum_counts": float(np.max(scan.pattern.intensity)),
            }
            for scan in measurement.scans
        )
        stages.append(
            ResultStage(
                key="measurement",
                title="The measurement",
                section="evidence",
                status="info",
                summary=(
                    f"{len(measurement)} scans of {label} with {radiation.name} "
                    f"({measurement.geometry}-tilting): "
                    f"{len({round(s.phi_deg % 360.0, 6) for s in measurement.scans})} azimuth(s), "
                    f"{len({round(s.psi_deg, 6) for s in measurement.scans})} tilt(s). The "
                    f"stress-free reflection lies at 2θ = {bragg:.3f}°."
                ),
                table=ResultTable(
                    columns=_SCAN_COLUMNS, rows=scan_rows, caption="Each scan as measured."
                ),
                explanation=(
                    "Tilts equally spaced in sin²ψ spread the points evenly along the line whose "
                    "slope is wanted; tilts of both signs test for shear; azimuths not 180° apart "
                    "fix the in-plane tensor."
                ),
                figures=(geometry_figure(result), peak_shift_figure(measurement, result)),
            )
        )
    else:
        stages.append(
            ResultStage(
                key="measurement",
                title="The measurement",
                section="evidence",
                status="info",
                summary=(
                    f"{len(result.peaks)} supplied peak positions of {label} with "
                    f"{radiation.name}; no scans, so no profile was located here."
                    + (
                        " No position carried an uncertainty, so the stress uncertainty is "
                        "taken from the scatter alone."
                        if all(peak.uncertainty_is_nominal for peak in result.peaks)
                        else ""
                    )
                ),
                figures=(geometry_figure(result),),
            )
        )

    # -- Evidence: the peaks ------------------------------------------------
    peak_rows = tuple(
        {
            "phi_deg": peak.phi_deg,
            "psi_deg": peak.psi_deg,
            "two_theta_deg": peak.two_theta_deg,
            "two_theta_uncertainty_mdeg": 1000.0 * peak.two_theta_uncertainty_deg,
            "fwhm_deg": _number(peak.fwhm_deg),
            "height": _number(peak.height),
            "reduced_chi_squared": _number(peak.reduced_chi_squared),
            "status": "converged" if peak.converged else "check",
        }
        for peak in result.peaks
    )
    methods = {peak.method for peak in result.peaks}
    method_label = ", ".join(_METHOD_LABELS[method] for method in sorted(methods))
    u_values = np.array([peak.two_theta_uncertainty_deg for peak in result.peaks])
    peak_figures: list[ResultFigure] = []
    if measurement is not None:
        peak_figures.append(peak_fits_figure(measurement, result))
        peak_figures.append(peak_quality_figure(result))
    stages.append(
        ResultStage(
            key="peak_positions",
            title="Peak positions",
            section="evidence",
            status="ok" if all(peak.converged for peak in result.peaks) else "warning",
            summary=(
                f"{len(result.peaks)} positions by {method_label.lower()}"
                + (
                    ", after the LPA correction"
                    if any(p.lpa_corrected for p in result.peaks)
                    else ""
                )
                + (", after Kα2 stripping" if any(p.kalpha2_stripped for p in result.peaks) else "")
                + f". Median u(2θ) = {1000.0 * float(np.median(u_values)):.2f} m°."
            ),
            metrics=(
                ResultMetric(
                    "Median u(2θ)", round(1000.0 * float(np.median(u_values)), 3), units="m°"
                ),
                ResultMetric(
                    "Largest u(2θ)", round(1000.0 * float(np.max(u_values)), 3), units="m°"
                ),
                ResultMetric("Located cleanly", sum(1 for p in result.peaks if p.converged)),
            ),
            table=ResultTable(columns=_PEAK_COLUMNS, rows=peak_rows, caption="Every located peak."),
            explanation=(
                "u(2θ) is the standard uncertainty of the position: from the profile-fit "
                "covariance scaled by that window's χ²ν, or propagated from the counting "
                "statistics for the parabola and the centroid. It becomes u(d) = d·cotθ·u(θ), "
                "the weight of the point in every fit that follows. The peak-fit χ²ν judges one "
                "profile against its counts; it is a different number from the strain-fit χ²ν."
            ),
            figures=tuple(peak_figures),
        )
    )

    # -- Evidence: the sin²ψ lines ------------------------------------------
    line_rows = tuple(
        {
            "phi_deg": line.phi_deg,
            "tilts": int(line.psi_deg.size),
            "intercept_angstrom": line.intercept_angstrom,
            "slope_mangstrom": 1000.0 * line.slope_angstrom,
            "slope_uncertainty_mangstrom": 1000.0 * line.slope_uncertainty_angstrom,
            "sigma_phi_mpa": line.sigma_phi_mpa,
            "sigma_phi_uncertainty_mpa": line.sigma_phi_uncertainty_mpa,
            "tau_phi_mpa": _number(line.tau_phi_mpa),
            "tau_phi_uncertainty_mpa": _number(line.tau_phi_uncertainty_mpa),
            "curvature_t": _number(line.curvature_t),
            "reduced_chi_squared": _number(line.reduced_chi_squared),
            "r_squared": _number(line.r_squared),
        }
        for line in result.regressions
    )
    stages.append(
        ResultStage(
            key="sin2psi_lines",
            title="d against sin²ψ at each azimuth",
            section="evidence",
            summary=(
                "At every azimuth d falls on a straight line in sin²ψ whose slope is "
                "d₀·½S₂·σφ: "
                + "; ".join(
                    f"{_phi()} = {line.phi_deg:g}°: σ{_phi()} = {line.sigma_phi_mpa:.0f} ± "
                    f"{line.sigma_phi_uncertainty_mpa:.0f} MPa"
                    for line in result.regressions
                )
                + "."
            ),
            table=ResultTable(
                columns=_LINE_COLUMNS, rows=line_rows, caption="The weighted line at each azimuth."
            ),
            explanation=(
                "Weighted least squares of d = c₀ + c₁·sin²ψ (+ c₂·sin2ψ when both signs of ψ "
                "were measured), with weights 1/u(d)²; uncertainties scaled by √χ²ν when it "
                "exceeds one. σφ − σ33 = c₁/(d₀·½S₂) and τφ = c₂/(d₀·½S₂). This is the "
                "classical evaluation, independent of the joint tensor fit, and the two should "
                "agree azimuth by azimuth."
            ),
            figures=(sin2psi_figure(result), strain_figure(result)),
        )
    )

    # -- Diagnostics --------------------------------------------------------
    diagnostic_figures = [
        figure
        for figure in (strain_residual_figure(result), stress_correlation_figure(result))
        if figure is not None
    ]
    if fit is not None:
        normalized = normalized_strain_residuals(result)
        outside = int(np.sum(np.abs(normalized) > 3.0))
        stages.append(
            ResultStage(
                key="fit_residuals",
                title="Residuals and correlations of the stress fit",
                section="diagnostics",
                status=("warning" if fit.reduced_chi_squared > 3.0 or outside > 0 else "ok"),
                summary=(
                    f"Strain-fit χ²ν = {fit.reduced_chi_squared:.2f}; {outside} of "
                    f"{normalized.size} residuals beyond ±3u; root-mean-square residual "
                    f"{1e6 * float(np.sqrt(np.mean(fit.residual_strain**2))):.1f} × 10⁻⁶."
                ),
                metrics=(
                    ResultMetric("Strain-fit χ²ν", round(fit.reduced_chi_squared, 4)),
                    ResultMetric("Degrees of freedom", fit.degrees_of_freedom),
                    ResultMetric("Residuals beyond ±3u", outside),
                ),
                explanation=(
                    "Under a correct model with correct uncertainties the normalized residuals "
                    "scatter as a unit normal. The correlation matrix says which components the "
                    "chosen directions separate well; σ11 and σ22 share the d₀ and S₁ terms."
                ),
                figures=tuple(diagnostic_figures),
            )
        )
    suggested = result.suggested_outliers()
    if result.excluded_indices or suggested:
        excluded_rows = tuple(
            {
                "phi_deg": result.peaks[index].phi_deg,
                "psi_deg": result.peaks[index].psi_deg,
                "two_theta_deg": result.peaks[index].two_theta_deg,
                "normalized_residual": (
                    None
                    if result.deleted_residuals is None
                    else _number(float(result.deleted_residuals[index]))
                ),
                "state": "excluded" if not result.included_mask[index] else "suggested",
            }
            for index in (*result.excluded_indices, *suggested)
        )
        stages.append(
            ResultStage(
                key="excluded_measurements",
                title="Excluded and suspect measurements",
                section="diagnostics",
                status="warning" if suggested else "info",
                summary=(
                    (
                        f"{len(result.excluded_indices)} measurement(s) were excluded by the "
                        "analyst and entered no fit. "
                        if result.excluded_indices
                        else ""
                    )
                    + (
                        f"{len(suggested)} included measurement(s) lie more than 3.5 of their own "
                        "standard uncertainties from the tensor fitted without them: candidates "
                        "to inspect in the peak-fit figure, and to exclude "
                        "by clicking them in the d against sin²ψ plot if the peak location is "
                        "bad."
                        if suggested
                        else "No included measurement lies more than 3.5u from the tensor "
                        "fitted without it."
                    )
                ),
                table=ResultTable(
                    columns=_EXCLUDED_COLUMNS,
                    rows=excluded_rows,
                    caption=(
                        "Deleted residual: each point's strain residual against the tensor "
                        "fitted without it, in units of its standard uncertainty (Birge-scaled). "
                        "A large value confirms an exclusion and a small one questions it."
                    ),
                ),
                explanation=(
                    "Exclude a measurement when its peak was located badly — a fit that missed "
                    "the peak, a spurious reflection, a detector artefact — and say why in the "
                    "record. Do not exclude points only because they spoil a straight line: "
                    "systematic curvature or oscillation of d against sin²ψ is a stress "
                    "gradient or texture, and removing it hides the physics. The deleted "
                    "residual judges each point against a fit that never saw it, so one bad "
                    "point cannot hide itself by pulling the fit, nor make its neighbours look "
                    "bad; a good point exceeds 3.5 about once in 2000."
                ),
            )
        )
    split = splitting_figure(result)
    curved = [
        line
        for line in result.regressions
        if math.isfinite(line.curvature_t) and abs(line.curvature_t) > 3.0
    ]
    stages.append(
        ResultStage(
            key="linearity",
            title="Linearity and ψ-splitting",
            section="diagnostics",
            status="warning" if curved else "ok",
            summary=(
                (
                    "d is straight in sin²ψ at every azimuth (no sin⁴ψ term beyond 3u)."
                    if not curved
                    else "d is curved in sin²ψ at "
                    + ", ".join(f"{_phi()} = {line.phi_deg:g}°" for line in curved)
                    + ": a gradient or texture the model does not describe."
                )
                + (
                    " Tilts of both signs were measured; the branch means and half-differences "
                    "are plotted."
                    if split is not None
                    else " No tilt was measured at both signs, so shear is untested."
                )
            ),
            metrics=tuple(
                ResultMetric(
                    f"Curvature / u at {_phi()} = {line.phi_deg:g}°", _number(line.curvature_t)
                )
                for line in result.regressions
            ),
            explanation=(
                "Curvature test: the line is refitted with an added sin⁴ψ term and its "
                "coefficient compared with its own uncertainty. ψ-splitting: at each |ψ| "
                "measured at both signs, a₁ = (d₊ + d₋)/2 is shear-free and a₂ = (d₊ − d₋)/2 = "
                "d₀·½S₂·τφ·sin|2ψ| (Dölle–Hauk)."
            ),
            figures=(split,) if split is not None else (),
        )
    )

    # -- Method --------------------------------------------------------------
    stages.append(
        ResultStage(
            key="theory",
            title="Theory: from peak shift to stress",
            section="method",
            status="info",
            summary=(
                "The spacing measured along a direction is a strain gauge along that direction; "
                "the tilt orients the gauge. Hooke's law through the reflection's diffraction "
                "elastic constants turns the strains into a stress."
            ),
            explanation=(
                "1. Bragg's law gives the spacing, d = λ/(2 sinθ), and the lattice strain along "
                "the scattering vector is ε(φ,ψ) = (d − d₀)/d₀.\n"
                "2. The scattering vector in the specimen frame is m = (cosφ sinψ, sinφ sinψ, "
                "cosψ), so ε(φ,ψ) = m·ε·m = ε11 cos²φ sin²ψ + ε12 sin2φ sin²ψ + ε22 sin²φ sin²ψ "
                "+ ε33 cos²ψ + ε13 cosφ sin2ψ + ε23 sinφ sin2ψ.\n"
                "3. For a quasi-isotropic polycrystal the grains that diffract respond through "
                "two constants: ε(φ,ψ) = ½S₂·m·σ·m + S₁·tr σ (isotropic: ½S₂ = (1 + ν)/E, "
                "S₁ = −ν/E).\n"
                "4. Expanded: ε(φ,ψ) = ½S₂[(σφ − σ33) sin²ψ + τφ sin2ψ + σ33] + S₁(σ11 + σ22 + "
                "σ33), with σφ = σ11 cos²φ + σ12 sin2φ + σ22 sin²φ and τφ = σ13 cosφ + σ23 sinφ.\n"
                "5. Plane stress (σi3 = 0): d(φ,ψ) = d₀[1 + ½S₂ σφ sin²ψ + S₁(σ11 + σ22)] — "
                "linear in sin²ψ with slope d₀·½S₂·σφ.\n"
                "6. The slope is insensitive to d₀: replacing d₀ by d(ψ = 0) changes σφ by a "
                "part in 10³. The intercept is not: it carries S₁(σ11 + σ22), which is how d₀ "
                "enters the tensor."
            ),
        )
    )
    dec = result.dec
    dec_figures: list[ResultFigure] = []
    if stiffness is not None:
        candidates: list[tuple[int, int, int]] = (
            [
                (1, 1, 0),
                (2, 0, 0),
                (2, 1, 1),
                (2, 2, 0),
                (3, 1, 0),
                (2, 2, 2),
                (3, 2, 1),
                (4, 0, 0),
                (4, 2, 0),
                (3, 3, 1),
            ]
            if cubic
            else [
                (1, 0, 0),
                (0, 0, 2),
                (1, 0, 1),
                (1, 0, 2),
                (1, 1, 0),
                (1, 0, 3),
                (1, 1, 2),
                (2, 0, 1),
                (0, 0, 4),
                (2, 0, 3),
            ]
        )
        reciprocal = phase.lattice.reciprocal_basis().matrix
        dec_figures.append(
            dec_figure(
                stiffness,
                normals=[reciprocal @ np.asarray(item, dtype=float) for item in candidates],
                labels=[plane_label(item, spec=spec) for item in candidates],
                chosen_normal=reciprocal @ np.asarray(reflection, dtype=float),
                chosen_label=label,
                used=dec,
                gammas=[_gamma(item) for item in candidates] if cubic else None,
                chosen_gamma=_gamma(reflection) if cubic else None,
            )
        )
    stages.append(
        ResultStage(
            key="elastic_constants",
            title="Diffraction elastic constants",
            section="method",
            status="info",
            summary=(
                f"S₁ = {dec.s1_per_tpa:.4f} TPa⁻¹ and ½S₂ = {dec.half_s2_per_tpa:.4f} TPa⁻¹ for "
                f"{label}, from {_MODEL_LABELS[dec.model]}; {stiffness_note}. Equivalent to "
                f"E = {dec.youngs_modulus_gpa:.1f} GPa and ν = {dec.poisson_ratio:.3f} for this "
                f"reflection. Relative uncertainty {100 * dec.relative_standard_uncertainty:.1f} %."
            ),
            metrics=(
                ResultMetric("S₁", round(dec.s1_per_tpa, 5), units="TPa⁻¹"),
                ResultMetric("½S₂", round(dec.half_s2_per_tpa, 5), units="TPa⁻¹"),
                ResultMetric("E of the reflection", round(dec.youngs_modulus_gpa, 2), units="GPa"),
                ResultMetric("ν of the reflection", round(dec.poisson_ratio, 4)),
            ),
            explanation=(
                "Each grain-interaction model gives the strain of a grain per unit macroscopic "
                "stress. Reuss: the single-crystal compliance (all grains carry the same stress). "
                "Voigt: the compliance of the Voigt-averaged aggregate (all grains share the "
                "strain). Kröner: each grain an Eshelby sphere in the self-consistent effective "
                "medium, A = (C_g + C*)⁻¹(C_eff + C*), with C_eff iterated to the orientation "
                "average of C_g A. The grains that diffract have the plane normal along the "
                "scattering vector and every rotation about it; averaging over that rotation "
                "leaves S₁ (strain along the normal per unit transverse stress) and S₁ + ½S₂ "
                "(per unit stress along the normal). All the tensor algebra is done in Mandel "
                "form, where the matrix inverse is the tensor inverse."
            ),
            figures=tuple(dec_figures),
        )
    )
    stages.append(
        ResultStage(
            key="algorithm",
            title="Algorithm",
            section="method",
            status="info",
            summary="Six steps, each reported above as its own stage.",
            explanation=(
                "1. Locate the peak in every scan within the fit window about the stress-free "
                "Bragg angle: optionally divide by the LPA factor; fit a Kα1/Kα2 pseudo-Voigt "
                "with a straight background (or strip Kα2 by Rachinger's recursion and fit a "
                "parabola or take a continuous centroid); keep the position and u(2θ).\n"
                "2. Convert to d and u(d) = d·cotθ·u(θ), then to ε = (d − d₀)/d₀.\n"
                "3. At each azimuth, fit d against sin²ψ (with a sin2ψ term when both signs of ψ "
                "were measured) and test for curvature.\n"
                "4. Fit all strains at once, ε = Aσ, by weighted linear least squares for the "
                "free components of the chosen stress state (A from step 3 of the theory); "
                "with d₀ refined, fit d = d₀ + A(d₀σ), which is linear in (d₀, d₀σ).\n"
                "5. Propagate: statistical (Birge-scaled), d₀ and elastic-constant "
                "sensitivities, combined in quadrature; cross-check by Monte Carlo.\n"
                "6. Derive the principal stresses, their direction and the von Mises equivalent, "
                "with the combined covariance propagated through each."
            ),
        )
    )
    return tuple(stages)


def _largest_source(budget: Any, component: int) -> str:
    """The budget source contributing most to one component."""

    return str(max(budget, key=lambda source: float(budget[source][component])))


def _gamma(indices: tuple[int, int, int]) -> float:
    """The cubic orientation parameter (h²k² + k²l² + l²h²)/(h² + k² + l²)²."""

    first, second, third = (float(value) for value in indices)
    square = first * first + second * second + third * third
    return (
        first * first * second * second
        + second * second * third * third
        + third * third * first * first
    ) / (square * square)


def _budget_summary(result: ResidualStressResult) -> str:
    fit = result.tensor
    assert fit is not None
    total = fit.combined_uncertainty_mpa
    parts = []
    for source, values in fit.budget_mpa.items():
        share = float(np.mean(values**2 / np.maximum(total**2, 1e-30)))
        parts.append(f"{source} {100.0 * share:.0f} %")
    text = "Share of the combined variance, averaged over the components: " + ", ".join(parts)
    if fit.monte_carlo_uncertainty_mpa is not None:
        ratio = fit.monte_carlo_uncertainty_mpa / np.maximum(total, 1e-30)
        text += (
            f". The Monte Carlo standard deviations are {float(np.min(ratio)):.2f}–"
            f"{float(np.max(ratio)):.2f} times the combined uncertainties over "
            f"{fit.monte_carlo_draws} draws."
        )
    return text + "."


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="xrd.residual_stress.ferrite_shot_peened",
            title="Residual stress: ferrite (211), Cr Kα",
            panel="xrd",
            summary=(
                "A compressive biaxial stress in ferrite, measured at three azimuths and nine "
                "tilts of both signs — the textbook sin²ψ measurement, with a known answer."
            ),
            teaches=(
                "Every d against sin²ψ line falls with tilt — compressive stress — and the three "
                "slopes fix the in-plane tensor. The recovered components match the generating "
                "stress within their uncertainties, and the budget shows u(d₀) dominating "
                "σ11 and σ22 while barely touching σ12."
            ),
            operation="xrd.residual_stress",
            request={"phase": {"builtin": "fe_bcc"}, "reflection": [2, 1, 1]},
        ),
        ExampleScenario(
            id="xrd.residual_stress.psi_splitting",
            title="Residual stress: ψ-splitting from out-of-plane shear",
            panel="xrd",
            summary=(
                "The same ferrite measurement with σ13 = 60 MPa: the ψ > 0 and ψ < 0 branches "
                "separate, and the evaluation with shear components recovers it."
            ),
            teaches=(
                "Out-of-plane shear opens the d against sin²ψ line into two branches. The a₂ "
                "half-difference is a straight line in sin|2ψ| whose slope is τφ, largest at "
                "φ = 0° where τφ = σ13."
            ),
            operation="xrd.residual_stress",
            request={
                "phase": {"builtin": "fe_bcc"},
                "reflection": [2, 1, 1],
                "true_sigma_13_mpa": 60.0,
                "stress_state": "biaxial_shear",
            },
        ),
        ExampleScenario(
            id="xrd.residual_stress.outlier",
            title="Residual stress: find and exclude a bad measurement",
            panel="xrd",
            summary=(
                "The ferrite measurement with one bad point — its peak displaced by 0.25° at "
                "φ = 45°, ψ = 30° — to find in the d against sin²ψ plot and exclude by clicking."
            ),
            teaches=(
                "The bad point sits far off its line, is ringed as a suggested outlier, and "
                "drags σ12 and the strain-fit χ²ν up. Click it in the plot and press Refit: it "
                "turns into a red cross, χ²ν falls back towards one and the tensor returns to the "
                "generating stress, while the point stays on the page for review."
            ),
            operation="xrd.residual_stress",
            request={
                "phase": {"builtin": "fe_bcc"},
                "reflection": [2, 1, 1],
                "demo_bad_points": "45 30",
            },
        ),
        ExampleScenario(
            id="xrd.residual_stress.synchrotron",
            title="Residual stress: synchrotron, 0.5 Å, χ-tilting",
            panel="xrd",
            summary=(
                "The ferrite (211) measurement at a synchrotron: a monochromatic 0.5 Å (24.8 keV) "
                "beam, polarized in the orbit plane, with χ-tilting and sharp peaks."
            ),
            teaches=(
                "At 0.5 Å the reflection sits near 2θ = 25°, so the peak shifts are small, and "
                "ω-tilting to 45° would take the beam below the surface — which is why "
                "synchrotron stress work uses χ-tilting. No Kα2 line is fitted, and the stress "
                "matches the laboratory measurement of the same specimen."
            ),
            operation="xrd.residual_stress",
            request={
                "phase": {"builtin": "fe_bcc"},
                "reflection": [2, 1, 1],
                "radiation": "monochromatic",
                "wavelength_angstrom": 0.5,
                "polarization_fraction": 0.95,
                "geometry": "chi",
                "demo_fwhm_deg": 0.3,
                "expected_fwhm_deg": 0.3,
                "window_deg": 3.0,
            },
        ),
        ExampleScenario(
            id="xrd.residual_stress.nickel_cu",
            title="Residual stress: nickel (420), Cu Kα, Reuss constants",
            panel="xrd",
            summary=(
                "A tensile stress in nickel on the (420) reflection at 2θ ≈ 155° with Cu Kα, "
                "using the Reuss constants and the parabola peak location."
            ),
            teaches=(
                "Rising lines mean tension. Change the elastic-constant model to Voigt and "
                "watch every stress scale by the ratio of the ½S₂ values — the elastic "
                "constants are a systematic, not a statistical, uncertainty."
            ),
            operation="xrd.residual_stress",
            request={
                "phase": {"builtin": "ni_fcc"},
                "reflection": [4, 2, 0],
                "radiation": "cu_ka_doublet",
                "dec_model": "reuss",
                "peak_method": "parabola",
                "true_sigma_11_mpa": 250.0,
                "true_sigma_22_mpa": 120.0,
                "true_sigma_12_mpa": -40.0,
                "demo_fwhm_deg": 0.8,
                "expected_fwhm_deg": 0.8,
            },
        ),
    )
)
