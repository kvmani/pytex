"""HRTEM simulation service: multislice and CTF evaluation for double-corrected TEM."""

from __future__ import annotations

import math
from typing import Any

from pytex.adapters.abtem import is_abtem_available, simulate_hrem
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    DocumentationLink,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
)


def _aberrations_from_request(
    request: dict[str, Any],
    voltage: float,
    defocus: float,
    cs_um: float,
    mode: DoubleCorrectionMode,
) -> MicroscopeAberrations:
    """Build the lens state from a workbench request.

    Both HRTEM operations declare the same aberration controls, so both read them
    the same way: a residual term the user leaves at zero simply does not enter the
    wave aberration. The coherence defaults preserve the behaviour the simulation
    operation had before those terms were exposed — a tighter focal spread under
    double correction — so an existing request without them is unchanged.
    """
    default_focal_spread = 5.0 if mode == DoubleCorrectionMode.DOUBLE_CORRECTED else 25.0
    aperture = float(request.get("aperture_cutoff_mrad", 0.0) or 0.0)
    return MicroscopeAberrations(
        energy_kev=voltage,
        defocus_angstrom=defocus,
        cs_mm=cs_um * 1e-3,
        c5_mm=float(request.get("c5_mm", 0.0)),
        astigmatism_angstrom=float(request.get("astigmatism_angstrom", 0.0)),
        astigmatism_angle_deg=float(request.get("astigmatism_angle_deg", 0.0)),
        coma_angstrom=float(request.get("coma_angstrom", 0.0)),
        coma_angle_deg=float(request.get("coma_angle_deg", 0.0)),
        trefoil_angstrom=float(request.get("trefoil_angstrom", 0.0)),
        trefoil_angle_deg=float(request.get("trefoil_angle_deg", 0.0)),
        focal_spread_angstrom=float(
            request.get("focal_spread_angstrom", default_focal_spread) or default_focal_spread
        ),
        convergence_semiangle_mrad=float(request.get("convergence_semiangle_mrad", 0.2)),
        aperture_cutoff_mrad=aperture if aperture > 0.0 else None,
        mode=mode,
    )


def _residual_rows(aberrations: MicroscopeAberrations) -> tuple[dict[str, str], ...]:
    """Result-table rows for the aberrations that actually shaped this result.

    A term the user left at zero is omitted rather than printed as a zero: the table
    should say what limited the image, not recite the whole coefficient list.
    """
    rows: list[dict[str, str]] = []
    if aberrations.c5_mm != 0.0:
        rows.append(
            {
                "metric": "Fifth-order spherical aberration C₅",
                "value": f"{aberrations.c5_mm:.4f}",
                "units": "mm",
            }
        )
    labels = {
        "astigmatism_2fold": "Two-fold astigmatism C₁₂",
        "axial_coma": "Axial coma C₂₁",
        "trefoil": "Trefoil C₂₃",
    }
    for key, _symbol, amplitude, angle in aberrations.residual_aberration_terms():
        rows.append({"metric": labels[key], "value": f"{amplitude:.1f}", "units": "Å"})
        rows.append({"metric": f"{labels[key]} azimuth", "value": f"{angle:.1f}", "units": "°"})
    return tuple(rows)


@REGISTRY.operation(
    "tem.simulate_hrem",
    title="HRTEM simulation",
    summary="Simulate high-resolution TEM phase-contrast micrographs and Thon-ring power spectra.",
    help_text=(
        "Simulate high-resolution transmission electron microscopy micrographs using multislice "
        "potential slicing or pure-Python phase-object transmission. Supports double-corrected "
        "electron optics, spherical and chromatic aberration damping, and specimen models "
        "including oriented crystal slabs, atomic vacancies, dislocation cores, "
        "and amorphous films."
    ),
    panel="tem_hrem",
    documentation=DocumentationLink("HRTEM multislice and CTF", "theory/hrem_multislice_and_ctf"),
    returns=(
        "Micrograph intensity array, base64 PNG preview, Thon rings FFT power spectrum, "
        "contrast metrics, and explainable optics report."
    ),
    parameters=(
        phase_parameter(
            builtin="si_diamond",
            help_text="Crystal structure for oriented crystalline slabs or defect lattices.",
        ),
        ChoiceParameter(
            name="sample_type",
            label="Sample morphology",
            options=(
                (
                    "crystalline",
                    "Perfect Crystal",
                    "Ideal periodic crystal slab oriented along zone axis",
                ),
                ("vacancy", "Point Defect Vacancy", "Crystal slab containing atomic vacancies"),
                (
                    "dislocation",
                    "Edge Dislocation Core",
                    "Crystal slab with Volterra displacement field",
                ),
                (
                    "amorphous",
                    "Amorphous Carbon Foil",
                    "Dense random packing for Thon rings and envelope calibration",
                ),
            ),
            default="crystalline",
            row="specimen",
            field_width="medium",
            help_text="Specimen microstructure morphology.",
        ),
        IndicesParameter(
            name="zone_axis",
            label="Zone axis [uvw]",
            default=[0, 0, 1],
            row="specimen",
            help_text="Incident electron beam direction in direct crystal basis.",
        ),
        IntegerParameter(
            name="supercell_xy",
            label="Supercell size",
            default=2,
            minimum=1,
            maximum=10,
            row="specimen",
            field_width="tiny",
            help_text="Number of unit-cell repeats along transverse x and y directions.",
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            default=200.0,
            minimum=20.0,
            maximum=1000.0,
            units="kV",
            row="optics_primary",
            field_width="short",
            help_text="Accelerating potential defining relativistic electron wavelength.",
        ),
        ChoiceParameter(
            name="mode",
            label="Optics mode",
            options=(
                ("uncorrected", "Uncorrected TEM", "Conventional uncorrected objective lens"),
                ("cs_corrected", "Cs-corrected", "Spherical aberration corrected"),
                ("double_corrected", "Double-corrected", "Cs and Cc aberration corrected"),
                ("ncsi", "Negative Cs Imaging", "Negative Cs imaging with overfocus"),
            ),
            default="double_corrected",
            row="optics_primary",
            field_width="medium",
            help_text="Aberration correction regime.",
        ),
        NumberParameter(
            name="defocus_angstrom",
            label="Objective lens defocus",
            default=-30.0,
            units="Å",
            symbol="defocus",
            row="aberrations",
            field_width="short",
            help_text="Objective lens defocus in Angstrom (underfocus is negative).",
        ),
        NumberParameter(
            name="cs_um",
            label="Spherical aberration coefficient",
            default=0.0,
            units="µm",
            symbol="spherical_aberration",
            row="aberrations",
            field_width="short",
            help_text="Spherical aberration coefficient Cs in micrometers.",
        ),
        NumberParameter(
            name="c5_mm",
            label="Fifth-order spherical aberration",
            default=0.0,
            units="mm",
            symbol="spherical_aberration_5th",
            row="aberrations",
            field_width="short",
            help_text=(
                "Fifth-order spherical aberration C5 in millimetres. Round, so it shifts phase "
                "transfer equally at every azimuth; it is what limits a Cs-corrected lens once "
                "the third-order term is tuned out."
            ),
        ),
        NumberParameter(
            name="astigmatism_angstrom",
            label="Two-fold astigmatism",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="astigmatism_2fold",
            row="astigmatism",
            field_width="short",
            help_text=(
                "Two-fold astigmatism amplitude C12 in Angstrom. Along its own azimuth it acts "
                "as a defocus offset of +C12 and across it as -C12, so it splits the point "
                "resolution between two orthogonal directions."
            ),
        ),
        NumberParameter(
            name="astigmatism_angle_deg",
            label="Astigmatism azimuth",
            default=0.0,
            minimum=0.0,
            maximum=180.0,
            units="°",
            symbol="astigmatism_2fold_azimuth",
            row="astigmatism",
            field_width="short",
            help_text=(
                "Azimuth phi12 of the two-fold astigmatism axis in the back focal plane. "
                "The term has period 180°, so azimuths beyond that repeat."
            ),
        ),
        NumberParameter(
            name="coma_angstrom",
            label="Axial coma",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="axial_coma",
            row="coma",
            field_width="short",
            help_text=(
                "Axial coma amplitude C21 in Angstrom. Single-fold, so it transfers differently "
                "in opposite directions and shifts image detail asymmetrically."
            ),
        ),
        NumberParameter(
            name="coma_angle_deg",
            label="Coma azimuth",
            default=0.0,
            minimum=0.0,
            maximum=360.0,
            units="°",
            symbol="axial_coma_azimuth",
            row="coma",
            field_width="short",
            help_text="Azimuth phi21 of the axial coma axis in the back focal plane.",
        ),
        NumberParameter(
            name="trefoil_angstrom",
            label="Trefoil",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="trefoil",
            row="trefoil",
            field_width="short",
            help_text=(
                "Three-fold astigmatism (trefoil) amplitude C23 in Angstrom. The dominant "
                "residual of many hexapole correctors."
            ),
        ),
        NumberParameter(
            name="trefoil_angle_deg",
            label="Trefoil azimuth",
            default=0.0,
            minimum=0.0,
            maximum=120.0,
            units="°",
            symbol="trefoil_azimuth",
            row="trefoil",
            field_width="short",
            help_text=(
                "Azimuth phi23 of the trefoil axis in the back focal plane. "
                "The term has period 120°, so azimuths beyond that repeat."
            ),
        ),
        NumberParameter(
            name="focal_spread_angstrom",
            label="Chromatic focal spread",
            default=10.0,
            minimum=0.0,
            maximum=100.0,
            units="Å",
            symbol="focal_spread",
            row="coherence",
            field_width="short",
            help_text="Temporal coherence 1/e focal spread in Angstrom.",
        ),
        NumberParameter(
            name="convergence_semiangle_mrad",
            label="Convergence semi-angle",
            default=0.2,
            minimum=0.0,
            maximum=5.0,
            units="mrad",
            row="coherence",
            field_width="short",
            help_text="Illumination convergence semi-angle in mrad.",
        ),
        NumberParameter(
            name="aperture_cutoff_mrad",
            label="Objective aperture cutoff",
            default=0.0,
            minimum=0.0,
            maximum=100.0,
            units="mrad",
            symbol="semiangle_cutoff",
            row="coherence",
            field_width="short",
            help_text=(
                "Objective aperture cutoff semi-angle in mrad. Zero means no aperture is "
                "inserted and transfer is limited by coherence alone."
            ),
        ),
        NumberParameter(
            name="sampling_angstrom",
            label="Pixel sampling pitch",
            default=0.2,
            minimum=0.05,
            maximum=0.5,
            units="Å/px",
            field_width="short",
            help_text="Real-space pixel sampling pitch in Angstrom per pixel.",
        ),
    ),
)
def _simulate_hrem(request: dict[str, Any]) -> dict[str, Any]:
    _phase_spec, phase = phase_from_request(request["phase"])
    sample_type = str(request.get("sample_type", "crystalline"))
    raw_zone = [int(x) for x in request.get("zone_axis", [0, 0, 1])]
    zone_axis = (raw_zone[0], raw_zone[1], raw_zone[2]) if len(raw_zone) >= 3 else (0, 0, 1)
    supercell_xy = int(request.get("supercell_xy", 2))
    voltage = float(request.get("beam_energy_kev", 200.0))
    mode_str = str(request.get("mode", "double_corrected"))
    defocus = float(request.get("defocus_angstrom", -30.0))
    cs_um = float(request.get("cs_um", 0.0))
    sampling = float(request.get("sampling_angstrom", 0.2))

    if sample_type == "amorphous":
        snap = AtomicSnapshot.amorphous_sample(
            species="C",
            density_g_cm3=1.8,
            dimensions_angstrom=(16.0, 16.0, 8.0),
            seed=42,
        )
    elif sample_type == "vacancy":
        snap = AtomicSnapshot.crystalline_with_vacancy(
            phase,
            supercell=(supercell_xy, supercell_xy, 2),
            vacancy_count=1,
            zone_axis=zone_axis,
        )
    elif sample_type == "dislocation":
        snap = AtomicSnapshot.crystalline_with_dislocation(
            phase,
            supercell=(supercell_xy, supercell_xy, 2),
            dislocation_type="edge",
        )
    else:
        snap = AtomicSnapshot.from_phase(
            phase,
            supercell=(supercell_xy, supercell_xy, 2),
            zone_axis=zone_axis,
        )

    mode = DoubleCorrectionMode(mode_str)
    aberr = _aberrations_from_request(request, voltage, defocus, cs_um, mode)

    result = simulate_hrem(snap, aberr, sampling_angstrom=sampling)

    table = ResultTable(
        columns=(
            Column("metric", "Metric", help_text="Physical or contrast diagnostic parameter"),
            Column("value", "Value", help_text="Measured numerical value"),
            Column("units", "Units", help_text="Dimensional units of the quantity"),
        ),
        rows=(
            {
                "metric": "Relativistic wavelength",
                "value": f"{aberr.wavelength_angstrom:.5f}",
                "units": "Å",
            },
            {"metric": "Accelerating voltage", "value": f"{aberr.energy_kev:.1f}", "units": "kV"},
            {
                "metric": "Objective lens defocus",
                "value": f"{aberr.defocus_angstrom:.1f}",
                "units": "Å",
            },
            {"metric": "Spherical aberration Cs", "value": f"{aberr.cs_um:.2f}", "units": "µm"},
            *_residual_rows(aberr),
            {
                "metric": "First zero point resolution",
                "value": f"{result.point_resolution_angstrom:.3f}",
                "units": "Å",
            },
            {
                "metric": "Information limit (1/e²)",
                "value": f"{result.information_limit_angstrom:.3f}",
                "units": "Å",
            },
            {
                "metric": "Michelson contrast",
                "value": f"{result.michelson_contrast * 100.0:.1f}",
                "units": "%",
            },
            {
                "metric": "Weber contrast",
                "value": f"{result.weber_contrast * 100.0:.1f}",
                "units": "%",
            },
            {
                "metric": "Micrograph field of view",
                "value": f"{result.extent_angstrom[0]:.1f} × {result.extent_angstrom[1]:.1f}",
                "units": "Å",
            },
            {
                "metric": "Image grid dimensions",
                "value": f"{result.image.shape[1]} × {result.image.shape[0]}",
                "units": "px",
            },
        ),
        caption="High-resolution electron microscopy imaging and optical performance metrics",
    )

    engine_note = (
        "Simulated using abTEM multislice wave propagation with Kirkland potential slicing."
        if is_abtem_available()
        else "Simulated using pure-Python phase-object transmission."
    )

    mode_display_labels = {
        "uncorrected": "Uncorrected TEM",
        "cs_corrected": "Cs-corrected",
        "double_corrected": "Double-corrected",
        "ncsi": "Negative Cs Imaging",
    }
    mode_label = mode_display_labels.get(mode.value, mode.value)

    summary_text = (
        f"Simulated HRTEM micrograph for {snap.label} under {mode_label} electron optics at "
        f"{aberr.energy_kev:.1f} kV (wavelength λ = {aberr.wavelength_angstrom:.5f} Å). "
        f"With defocus Δf = {aberr.defocus_angstrom:.1f} Å and Cs = {aberr.cs_um:.2f} µm, "
        f"the point resolution is {result.point_resolution_angstrom:.2f} Å and the information "
        f"limit is {result.information_limit_angstrom:.2f} Å. The simulated image exhibits a "
        f"Michelson contrast of {result.michelson_contrast * 100.0:.1f}% across a field of view of "
        f"{result.extent_angstrom[0]:.1f} × {result.extent_angstrom[1]:.1f} Å."
    )
    if aberr.has_azimuthal_aberrations:
        residual_written = ", ".join(
            f"{symbol} = {amplitude:.1f} Å at {angle:.1f}°"
            for _key, symbol, amplitude, angle in aberr.residual_aberration_terms()
        )
        summary_text += (
            f" Residual non-round aberrations {residual_written} entered the wave aberration, so "
            "the image is not resolved equally in every direction; the CTF view reports the "
            "azimuthal spread."
        )

    app_res = AppResult(
        title=f"HRTEM Simulation: {snap.label}",
        summary=summary_text,
        table=table,
        data={
            "image_png": result.to_png_base64(),
            "power_spectrum_png": result.to_power_spectrum_base64(),
            "michelson_contrast": result.michelson_contrast,
            "weber_contrast": result.weber_contrast,
            "point_resolution_angstrom": result.point_resolution_angstrom,
            "information_limit_angstrom": result.information_limit_angstrom,
            "pixel_size_angstrom": result.pixel_size_angstrom,
            "extent_angstrom": list(result.extent_angstrom),
            "natoms": snap.natoms,
            "sample_label": snap.label,
            "has_azimuthal_aberrations": aberr.has_azimuthal_aberrations,
            "residual_aberrations": [
                {
                    "key": key,
                    "symbol": symbol,
                    "amplitude_angstrom": amplitude,
                    "azimuth_deg": angle,
                }
                for key, symbol, amplitude, angle in aberr.residual_aberration_terms()
            ],
            "description": result.describe(),
        },
        inputs=dict(request),
        notes=(engine_note,),
        citations=(
            "Kirkland (2010), Advanced Computing in Electron Microscopy, 2nd ed.",
            "Cowley & Moodie (1957), Acta Crystallogr. 10, 609-619.",
            "Madsen et al. (2021), abTEM: An open-source framework, ChemPhysChem 22, 1-13.",
        ),
    )
    return app_res.to_json()


@REGISTRY.operation(
    "tem.ctf_calculator",
    title="CTF calculator",
    summary=(
        "Calculate objective lens Contrast Transfer Function (CTF), damping envelopes, "
        "and resolution limits."
    ),
    help_text=(
        "Calculate the 1D Contrast Transfer Function and Frank partial coherence damping envelopes "
        "for high-resolution electron microscopy. Computes point resolution at the Scherzer zero "
        "crossing, temporal and spatial coherence damping, and the information limit."
    ),
    panel="tem_hrem",
    documentation=DocumentationLink("HRTEM multislice and CTF", "theory/hrem_multislice_and_ctf"),
    returns=(
        "1D CTF curves, phase shift, damping envelopes, point resolution, and information limit."
    ),
    parameters=(
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            default=200.0,
            minimum=20.0,
            maximum=1000.0,
            units="kV",
            row="optics_primary",
            field_width="short",
            help_text="Accelerating potential in kV.",
        ),
        ChoiceParameter(
            name="mode",
            label="Optics mode",
            options=(
                ("uncorrected", "Uncorrected TEM", "Conventional uncorrected objective lens"),
                ("cs_corrected", "Cs-corrected", "Spherical aberration corrected"),
                ("double_corrected", "Double-corrected", "Cs and Cc aberration corrected"),
                ("ncsi", "Negative Cs Imaging", "Negative Cs imaging with overfocus"),
            ),
            default="double_corrected",
            row="optics_primary",
            field_width="medium",
            help_text="Aberration correction mode.",
        ),
        NumberParameter(
            name="defocus_angstrom",
            label="Objective lens defocus",
            default=-30.0,
            units="Å",
            symbol="defocus",
            row="aberrations",
            field_width="short",
            help_text="Objective lens defocus in Angstrom.",
        ),
        NumberParameter(
            name="cs_um",
            label="Spherical aberration coefficient",
            default=0.0,
            units="µm",
            symbol="spherical_aberration",
            row="aberrations",
            field_width="short",
            help_text="Spherical aberration coefficient Cs in micrometers.",
        ),
        NumberParameter(
            name="c5_mm",
            label="Fifth-order spherical aberration",
            default=0.0,
            units="mm",
            symbol="spherical_aberration_5th",
            row="aberrations",
            field_width="short",
            help_text=(
                "Fifth-order spherical aberration C5 in millimetres. Round, so it shifts phase "
                "transfer equally at every azimuth; it is what limits a Cs-corrected lens once "
                "the third-order term is tuned out."
            ),
        ),
        NumberParameter(
            name="astigmatism_angstrom",
            label="Two-fold astigmatism",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="astigmatism_2fold",
            row="astigmatism",
            field_width="short",
            help_text=(
                "Two-fold astigmatism amplitude C12 in Angstrom. Along its own azimuth it acts "
                "as a defocus offset of +C12 and across it as -C12, so it splits the point "
                "resolution between two orthogonal directions."
            ),
        ),
        NumberParameter(
            name="astigmatism_angle_deg",
            label="Astigmatism azimuth",
            default=0.0,
            minimum=0.0,
            maximum=180.0,
            units="°",
            symbol="astigmatism_2fold_azimuth",
            row="astigmatism",
            field_width="short",
            help_text=(
                "Azimuth phi12 of the two-fold astigmatism axis in the back focal plane. "
                "The term has period 180°, so azimuths beyond that repeat."
            ),
        ),
        NumberParameter(
            name="coma_angstrom",
            label="Axial coma",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="axial_coma",
            row="coma",
            field_width="short",
            help_text=(
                "Axial coma amplitude C21 in Angstrom. Single-fold, so it transfers differently "
                "in opposite directions and shifts image detail asymmetrically."
            ),
        ),
        NumberParameter(
            name="coma_angle_deg",
            label="Coma azimuth",
            default=0.0,
            minimum=0.0,
            maximum=360.0,
            units="°",
            symbol="axial_coma_azimuth",
            row="coma",
            field_width="short",
            help_text="Azimuth phi21 of the axial coma axis in the back focal plane.",
        ),
        NumberParameter(
            name="trefoil_angstrom",
            label="Trefoil",
            default=0.0,
            minimum=0.0,
            units="Å",
            symbol="trefoil",
            row="trefoil",
            field_width="short",
            help_text=(
                "Three-fold astigmatism (trefoil) amplitude C23 in Angstrom. The dominant "
                "residual of many hexapole correctors."
            ),
        ),
        NumberParameter(
            name="trefoil_angle_deg",
            label="Trefoil azimuth",
            default=0.0,
            minimum=0.0,
            maximum=120.0,
            units="°",
            symbol="trefoil_azimuth",
            row="trefoil",
            field_width="short",
            help_text=(
                "Azimuth phi23 of the trefoil axis in the back focal plane. "
                "The term has period 120°, so azimuths beyond that repeat."
            ),
        ),
        NumberParameter(
            name="focal_spread_angstrom",
            label="Chromatic focal spread",
            default=10.0,
            minimum=0.0,
            maximum=100.0,
            units="Å",
            symbol="focal_spread",
            row="coherence",
            field_width="short",
            help_text="Temporal coherence 1/e focal spread in Angstrom.",
        ),
        NumberParameter(
            name="convergence_semiangle_mrad",
            label="Convergence semi-angle",
            default=0.2,
            minimum=0.0,
            maximum=5.0,
            units="mrad",
            row="coherence",
            field_width="short",
            help_text="Illumination convergence semi-angle in mrad.",
        ),
        NumberParameter(
            name="aperture_cutoff_mrad",
            label="Objective aperture cutoff",
            default=25.0,
            minimum=0.0,
            maximum=100.0,
            units="mrad",
            symbol="semiangle_cutoff",
            row="coherence",
            field_width="short",
            help_text="Objective aperture cutoff semi-angle in mrad.",
        ),
        NumberParameter(
            name="max_q_inv_angstrom",
            label="Frequency range",
            default=2.5,
            minimum=0.5,
            maximum=10.0,
            units="Å⁻¹",
            row="cut",
            field_width="short",
            help_text=(
                "Upper spatial frequency of the profile. A corrected instrument transfers well "
                "beyond the 2.5 Å⁻¹ that suits an uncorrected one."
            ),
        ),
        NumberParameter(
            name="azimuth_deg",
            label="Cut azimuth",
            default=0.0,
            minimum=0.0,
            maximum=360.0,
            units="°",
            symbol="ctf_azimuth",
            row="cut",
            field_width="short",
            help_text=(
                "Azimuth of the radial cut through the back focal plane. A round lens transfers "
                "identically at every azimuth; once a residual term is nonzero this selects one "
                "direction, and the azimuthal band shows the spread across all of them."
            ),
        ),
    ),
)
def _calculate_ctf(request: dict[str, Any]) -> dict[str, Any]:
    voltage = float(request.get("beam_energy_kev", 200.0))
    mode_str = str(request.get("mode", "double_corrected"))
    defocus = float(request.get("defocus_angstrom", -30.0))
    cs_um = float(request.get("cs_um", 0.0))
    azimuth = float(request.get("azimuth_deg", 0.0))

    mode = DoubleCorrectionMode(mode_str)
    aberr = _aberrations_from_request(
        {"focal_spread_angstrom": 10.0, "aperture_cutoff_mrad": 25.0, **request},
        voltage,
        defocus,
        cs_um,
        mode,
    )

    max_q = float(request.get("max_q_inv_angstrom", 2.5))
    ctf = aberr.evaluate_ctf_1d(max_q_inv_angstrom=max_q, num_points=500, azimuth_deg=azimuth)
    band = aberr.evaluate_ctf_azimuthal(
        max_q_inv_angstrom=max_q, num_radial=500, num_azimuthal=72
    )

    first_zero_str = (
        f"{ctf.point_resolution_angstrom:.3f}"
        if not math.isnan(ctf.point_resolution_angstrom)
        else "N/A"
    )
    info_limit_str = (
        f"{ctf.information_limit_angstrom:.3f}"
        if not math.isnan(ctf.information_limit_angstrom)
        else "N/A"
    )
    best_over_azimuth, worst_over_azimuth = band.point_resolution_range_angstrom

    table = ResultTable(
        columns=(
            Column("parameter", "Parameter", help_text="Optical characteristic"),
            Column("value", "Value", help_text="Computed value"),
            Column("units", "Units", help_text="Physical units"),
        ),
        rows=(
            {
                "parameter": "Accelerating voltage",
                "value": f"{aberr.energy_kev:.1f}",
                "units": "kV",
            },
            {
                "parameter": "Relativistic wavelength",
                "value": f"{aberr.wavelength_angstrom:.5f}",
                "units": "Å",
            },
            {
                "parameter": "Interaction parameter σ",
                "value": f"{aberr.interaction_parameter_inv_v_angstrom:.6f}",
                "units": "V⁻¹Å⁻¹",
            },
            {
                "parameter": "Objective lens defocus",
                "value": f"{aberr.defocus_angstrom:.1f}",
                "units": "Å",
            },
            {"parameter": "Spherical aberration Cs", "value": f"{aberr.cs_um:.2f}", "units": "µm"},
            {
                "parameter": "First zero crossing frequency",
                "value": f"{ctf.first_zero_q_inv_angstrom:.3f}",
                "units": "Å⁻¹",
            },
            {"parameter": "First zero point resolution", "value": first_zero_str, "units": "Å"},
            {
                "parameter": "Information limit frequency",
                "value": f"{ctf.information_limit_q_inv_angstrom:.3f}",
                "units": "Å⁻¹",
            },
            {"parameter": "Information limit d-spacing", "value": info_limit_str, "units": "Å"},
            {
                "parameter": "Scherzer defocus Δf_Sch",
                "value": f"{aberr.scherzer_defocus_angstrom:.1f}",
                "units": "Å",
            },
            {
                "parameter": "Scherzer point resolution",
                "value": f"{aberr.scherzer_resolution_angstrom:.2f}",
                "units": "Å",
            },
            *(
                {"parameter": row["metric"], "value": row["value"], "units": row["units"]}
                for row in _residual_rows(aberr)
            ),
            *(
                (
                    {
                        "parameter": "Cut azimuth θ",
                        "value": f"{ctf.azimuth_deg:.1f}",
                        "units": "°",
                    },
                    {
                        "parameter": "Point resolution over azimuth",
                        "value": f"{best_over_azimuth:.3f} – {worst_over_azimuth:.3f}",
                        "units": "Å",
                    },
                    {
                        "parameter": "Resolution anisotropy",
                        "value": f"{band.resolution_anisotropy_angstrom:.3f}",
                        "units": "Å",
                    },
                    {
                        "parameter": "Coarsest transfer azimuth",
                        "value": f"{band.worst_azimuth_deg:.1f}",
                        "units": "°",
                    },
                )
                if aberr.has_azimuthal_aberrations and not math.isnan(best_over_azimuth)
                else ()
            ),
        ),
        caption="Objective lens Contrast Transfer Function optical properties",
    )

    mode_display_labels = {
        "uncorrected": "Uncorrected TEM",
        "cs_corrected": "Cs-corrected",
        "double_corrected": "Double-corrected",
        "ncsi": "Negative Cs Imaging",
    }
    mode_label = mode_display_labels.get(mode.value, mode.value)

    summary_text = (
        f"Objective lens Contrast Transfer Function for {mode_label} optics at "
        f"{aberr.energy_kev:.1f} kV (wavelength λ = {aberr.wavelength_angstrom:.5f} Å). "
        f"Defocus Δf = {aberr.defocus_angstrom:.1f} Å and Cs = {aberr.cs_um:.2f} µm "
        f"produce a point resolution of {first_zero_str} Å and an information limit "
        f"of {info_limit_str} Å."
    )
    if aberr.has_azimuthal_aberrations:
        residual_written = ", ".join(
            f"{symbol} = {amplitude:.1f} Å at {angle:.1f}°"
            for _key, symbol, amplitude, angle in aberr.residual_aberration_terms()
        )
        summary_text += (
            f" That figure is the cut at azimuth θ = {ctf.azimuth_deg:.1f}°, not the whole lens: "
            f"the residual non-round aberrations {residual_written} make phase transfer depend on "
            "direction."
        )
        if not math.isnan(best_over_azimuth):
            summary_text += (
                f" Across azimuth the point resolution runs from {best_over_azimuth:.2f} Å to "
                f"{worst_over_azimuth:.2f} Å, an anisotropy of "
                f"{band.resolution_anisotropy_angstrom:.2f} Å, coarsest at "
                f"{band.worst_azimuth_deg:.1f}°."
            )

    app_res = AppResult(
        title=f"CTF: {aberr.energy_kev:.0f} kV ({mode_label})",
        summary=summary_text,
        table=table,
        data={
            "spatial_frequencies": ctf.spatial_frequencies_inv_angstrom.tolist(),
            "transfer_function": ctf.transfer_function.tolist(),
            "undamped_ctf": ctf.ctf_undamped.tolist(),
            "total_envelope": ctf.total_envelope.tolist(),
            "temporal_envelope": ctf.temporal_envelope.tolist(),
            "spatial_envelope": ctf.spatial_envelope.tolist(),
            "point_resolution_angstrom": ctf.point_resolution_angstrom,
            "information_limit_angstrom": ctf.information_limit_angstrom,
            "azimuth_deg": ctf.azimuth_deg,
            "has_azimuthal_aberrations": aberr.has_azimuthal_aberrations,
            "azimuthal_transfer_min": band.transfer_min.tolist(),
            "azimuthal_transfer_max": band.transfer_max.tolist(),
            "azimuthal_point_resolution_best_angstrom": best_over_azimuth,
            "azimuthal_point_resolution_worst_angstrom": worst_over_azimuth,
            "resolution_anisotropy_angstrom": band.resolution_anisotropy_angstrom,
            "worst_azimuth_deg": band.worst_azimuth_deg,
            "azimuths_deg": band.azimuths_deg.tolist(),
            "azimuthal_point_resolution_angstrom": [
                None if math.isnan(value) else float(value)
                for value in band.point_resolution_angstrom
            ],
            "residual_aberrations": [
                {
                    "key": key,
                    "symbol": symbol,
                    "amplitude_angstrom": amplitude,
                    "azimuth_deg": angle,
                }
                for key, symbol, amplitude, angle in aberr.residual_aberration_terms()
            ],
            "description": ctf.describe(),
            "azimuthal_description": band.describe(),
        },
        inputs=dict(request),
        notes=(
            "Temporal coherence envelope modelled via Frank's Gaussian focal spread approximation.",
            "Spatial coherence envelope modelled via convergence angle illumination integration.",
            (
                "Frank's envelopes are isotropic, so the azimuthal spread reported here is that "
                "of the transfer oscillation, not of the coherence damping."
            ),
        ),
        citations=(
            "Frank (1973), An envelope for the transfer function, Optik 38, 519-536.",
            "Scherzer (1949), The theoretical resolution limit, J. Appl. Phys. 20, 20-29.",
            "Kirkland (2010), Advanced Computing in Electron Microscopy, 2nd ed.",
        ),
    )
    return app_res.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="hrem.example.si_dumbbell",
            title="Silicon [110] atomic dumbbells",
            panel="tem_hrem",
            summary=(
                "Diamond-cubic silicon oriented along the [110] zone axis under "
                "double-corrected electron optics."
            ),
            teaches=(
                "Notice the resolved silicon atomic column dumbbells separated by 1.36 Å "
                "along the [110] projection. Spherical and chromatic aberration correction "
                "extends the passband well beyond the Scherzer boundary, transferring high "
                "spatial frequencies without contrast reversals."
            ),
            operation="tem.simulate_hrem",
            request={
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 2,
                "beam_energy_kev": 200.0,
                "mode": "double_corrected",
                "defocus_angstrom": -30.0,
                "cs_um": 0.0,
                "sampling_angstrom": 0.2,
            },
        ),
        ExampleScenario(
            id="hrem.example.graphene_vacancy",
            title="Carbon vacancy in negative Cs imaging mode",
            panel="tem_hrem",
            summary=(
                "Diamond cubic carbon lattice containing an atomic vacancy imaged under "
                "Negative Cs Imaging (NCSI) optics."
            ),
            teaches=(
                "Notice that negative spherical aberration combined with a small overfocus "
                "produces bright atom contrast on a dark background. The missing atom column "
                "produces a sharp intensity drop with minimal Fresnel fringing around the "
                "defect core."
            ),
            operation="tem.simulate_hrem",
            request={
                "phase": {"builtin": "diamond"},
                "sample_type": "vacancy",
                "zone_axis": [0, 0, 1],
                "supercell_xy": 2,
                "beam_energy_kev": 200.0,
                "mode": "ncsi",
                "defocus_angstrom": 50.0,
                "cs_um": -15.0,
                "sampling_angstrom": 0.2,
            },
        ),
        ExampleScenario(
            id="hrem.example.amorphous_thon_rings",
            title="Amorphous carbon foil and Thon rings",
            panel="tem_hrem",
            summary=(
                "Amorphous carbon thin foil demonstrating Thon-ring formation in the "
                "fast Fourier transform power spectrum."
            ),
            teaches=(
                "Notice the concentric dark rings of zero contrast in the 2D power spectrum. "
                "These Thon rings correspond directly to the zero crossings of the objective lens "
                "contrast transfer function, providing a direct experimental method to calibrate "
                "defocus, astigmatism, and coherence damping."
            ),
            operation="tem.simulate_hrem",
            request={
                "phase": {"builtin": "diamond"},
                "sample_type": "amorphous",
                "zone_axis": [0, 0, 1],
                "supercell_xy": 2,
                "beam_energy_kev": 300.0,
                "mode": "uncorrected",
                "defocus_angstrom": -600.0,
                "cs_um": 1000.0,
                "sampling_angstrom": 0.25,
            },
        ),
        ExampleScenario(
            id="hrem.example.ctf_double_corrected",
            title="Double-corrected 300 kV objective lens transfer",
            panel="tem_hrem",
            summary=(
                "Contrast Transfer Function and partial coherence damping envelopes for a "
                "double-corrected 300 kV microscope."
            ),
            teaches=(
                "Notice that reducing both spherical aberration Cs and focal spread extends "
                "the transfer function envelope past 1.5 reciprocal angstroms. The information "
                "limit surpasses the classical Scherzer point resolution, enabling true "
                "sub-angstrom phase-contrast imaging."
            ),
            operation="tem.ctf_calculator",
            request={
                "beam_energy_kev": 300.0,
                "mode": "double_corrected",
                "defocus_angstrom": -20.0,
                "cs_um": 0.0,
                "focal_spread_angstrom": 5.0,
                "convergence_semiangle_mrad": 0.1,
                "aperture_cutoff_mrad": 35.0,
            },
        ),
        ExampleScenario(
            id="hrem.example.ctf_residual_astigmatism",
            title="Residual two-fold astigmatism after correction",
            panel="tem_hrem",
            summary=(
                "A corrected 300 kV lens carrying 20 Å of residual two-fold astigmatism, "
                "showing how phase transfer then depends on direction."
            ),
            teaches=(
                "Notice that the transfer curve is now one cut through a lens that transfers "
                "differently in different directions, and that the shaded band spans every "
                "azimuth. Two-fold astigmatism enters the wave aberration as a cosine of twice "
                "the azimuth, so along its own axis it acts as a defocus offset of +C₁₂ and "
                "across it as −C₁₂: the point resolution splits between two orthogonal "
                "directions 90° apart. This is why correcting Cs alone does not settle the "
                "resolution of a corrected instrument — once the round terms are tuned out, the "
                "residual non-round terms are what remains, and the anisotropy is the quantity "
                "worth minimising."
            ),
            operation="tem.ctf_calculator",
            request={
                "beam_energy_kev": 300.0,
                "mode": "double_corrected",
                "defocus_angstrom": -50.0,
                "cs_um": 1.0,
                "astigmatism_angstrom": 20.0,
                "astigmatism_angle_deg": 30.0,
                "focal_spread_angstrom": 8.0,
                "convergence_semiangle_mrad": 0.1,
                "aperture_cutoff_mrad": 0.0,
                "max_q_inv_angstrom": 2.5,
                "azimuth_deg": 30.0,
            },
        ),
    )
)
