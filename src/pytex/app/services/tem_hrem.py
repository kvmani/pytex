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
    aberr = MicroscopeAberrations(
        energy_kev=voltage,
        defocus_angstrom=defocus,
        cs_mm=cs_um * 1e-3,
        focal_spread_angstrom=5.0 if mode == DoubleCorrectionMode.DOUBLE_CORRECTED else 25.0,
        mode=mode,
    )

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
            name="focal_spread_angstrom",
            label="Chromatic focal spread",
            default=10.0,
            minimum=0.0,
            maximum=100.0,
            units="Å",
            symbol="chromatic_aberration",
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
            field_width="short",
            help_text="Objective aperture cutoff semi-angle in mrad.",
        ),
    ),
)
def _calculate_ctf(request: dict[str, Any]) -> dict[str, Any]:
    voltage = float(request.get("beam_energy_kev", 200.0))
    mode_str = str(request.get("mode", "double_corrected"))
    defocus = float(request.get("defocus_angstrom", -30.0))
    cs_um = float(request.get("cs_um", 0.0))
    focal_spread = float(request.get("focal_spread_angstrom", 10.0))
    alpha_s = float(request.get("convergence_semiangle_mrad", 0.2))
    aperture = float(request.get("aperture_cutoff_mrad", 25.0))

    mode = DoubleCorrectionMode(mode_str)
    aberr = MicroscopeAberrations(
        energy_kev=voltage,
        defocus_angstrom=defocus,
        cs_mm=cs_um * 1e-3,
        focal_spread_angstrom=focal_spread,
        convergence_semiangle_mrad=alpha_s,
        aperture_cutoff_mrad=aperture,
        mode=mode,
    )

    ctf = aberr.evaluate_ctf_1d(max_q_inv_angstrom=2.5, num_points=500)

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
            "description": ctf.describe(),
        },
        inputs=dict(request),
        notes=(
            "Temporal coherence envelope modelled via Frank's Gaussian focal spread approximation.",
            "Spatial coherence envelope modelled via convergence angle illumination integration.",
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
    )
)
