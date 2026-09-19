"""HRTEM simulation service: multislice and CTF evaluation for double-corrected TEM."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import replace
from typing import Any

import numpy as np

from pytex.adapters.abtem import is_abtem_available, simulate_hrem
from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    ChoiceParameter,
    DocumentationLink,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ObjectParameter,
    Parameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import phase_parameter, plane_label
from pytex.app.services.tem_figures import (
    beam_thickness_figure,
    focal_contrast_figure,
    hrem_spectrum_figure,
    hrtem_tableau_figure,
)
from pytex.app.uploads import describe_upload, uploaded_name_and_text
from pytex.core._chemistry import atomic_number
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    _array_png_data_url,
)
from pytex.diffraction.multislice import (
    HREMEngine,
    MultisliceExitWave,
    MultisliceGrid,
    PotentialParametrization,
    TemporalCoherence,
    multislice,
    multislice_summary,
    periodic_slab,
)

#: File kinds the imported-structure specimen reads.
STRUCTURE_FILE_SUFFIXES = (".xyz", ".extxyz")

#: The largest specimen the workbench simulates. The phase-object fallback
#: places atoms one at a time, and a request that would run for many minutes
#: is refused with a hint rather than left to time out.
MAX_SPECIMEN_ATOMS = 20000


def _specimen_from_request(
    request: Mapping[str, Any], phase: Any
) -> tuple[AtomicSnapshot, dict[str, Any]]:
    """Build the specimen a simulation request describes, and say how it was built.

    Returns the snapshot and a record of the specimen for the report: where it
    came from, the thickness asked for, the thickness delivered and the repeats
    used. A crystal slab is a whole number of unit cells, so the delivered
    thickness is generally a little more than the request, and the report says
    so rather than quoting the request as though it were exact.
    """

    sample_type = str(request.get("sample_type", "crystalline"))
    raw_zone = [int(x) for x in request.get("zone_axis", [0, 0, 1])]
    zone_axis = (raw_zone[0], raw_zone[1], raw_zone[2]) if len(raw_zone) >= 3 else (0, 0, 1)
    lateral = int(request.get("supercell_xy", 2))
    thickness = float(request.get("thickness_angstrom", 10.0))
    info: dict[str, Any] = {
        "sample_type": sample_type,
        "requested_thickness_angstrom": thickness,
        "unit_cell_repeats": None,
        "structure_file": None,
        "zone_axis_cell": None,
        "_zone_axis_cell": None,
    }

    if sample_type == "imported":
        name, text = uploaded_name_and_text(
            request.get("structure_file"), field="structure_file", suffixes=STRUCTURE_FILE_SUFFIXES
        )
        lines = text.splitlines()
        if len(lines) < 3:
            raise InvalidInputError(
                f"{name} is too short to be an XYZ structure.",
                field="structure_file",
                hint="An XYZ file holds an atom count, a comment line, then one line per atom.",
            )
        try:
            # Rejoined with newlines so the reader always receives content: a
            # single-line upload must never be mistaken for a path on the server.
            snapshot = AtomicSnapshot.from_xyz("\n".join(lines) + "\n").prepared_for_imaging(
                periodic_xy="lattice=" in lines[1].lower()
            )
        except ValueError as error:
            raise InvalidInputError(
                f"{name} could not be read as an XYZ structure: {error}",
                field="structure_file",
                hint=(
                    "Check the atom count on the first line, one 'element x y z' line per atom, "
                    "and an orthogonal Lattice entry if the comment line carries one."
                ),
            ) from error
        unknown = []
        for species in sorted(set(snapshot.species)):
            try:
                atomic_number(species)
            except ValueError:
                unknown.append(species)
        if unknown:
            raise InvalidInputError(
                f"{name} names species that are not elements: {', '.join(unknown)}.",
                field="structure_file",
                hint=(
                    "Write element symbols in the first column. A molecular-dynamics dump that "
                    "numbers atom types needs them mapped to elements before export."
                ),
            )
        label = snapshot.label if snapshot.label != "Imported XYZ snapshot" else ""
        snapshot = replace(snapshot, label=f"{name}: {label}" if label else name)
        info["structure_file"] = describe_upload(name, text)
        info["slab_thickness_angstrom"] = float(np.ptp(snapshot.positions[:, 2]))
    elif sample_type == "amorphous":
        snapshot = AtomicSnapshot.amorphous_sample(
            species="C",
            density_g_cm3=1.8,
            dimensions_angstrom=(16.0, 16.0, thickness),
            seed=42,
        )
        info["slab_thickness_angstrom"] = thickness
    else:
        # A crystal is built in an exact lattice-periodic box: the multislice
        # grid is periodic, and a box that is not a lattice period puts a seam at
        # every edge and a forbidden reflection in every pattern.
        slab_axis = (0, 0, 1) if sample_type == "dislocation" else zone_axis
        try:
            if sample_type == "dislocation":
                repeats_abc, slab = AtomicSnapshot.repeats_for_thickness(
                    phase, slab_axis, thickness, lateral
                )
                snapshot = AtomicSnapshot.crystalline_with_dislocation(
                    phase, supercell=repeats_abc, dislocation_type="edge"
                )
                repeats: tuple[int, ...] = repeats_abc
            else:
                snapshot, cell, repeats = periodic_slab(
                    phase, zone_axis, (lateral, lateral), thickness_angstrom=thickness
                )
                slab = float(snapshot.cell[2, 2])
                info["zone_axis_cell"] = cell.describe()
                info["_zone_axis_cell"] = cell
                if sample_type == "vacancy":
                    centre = 0.5 * np.diag(snapshot.cell)
                    nearest = int(np.argmin(np.linalg.norm(snapshot.positions - centre, axis=1)))
                    keep = np.arange(snapshot.natoms) != nearest
                    snapshot = replace(
                        snapshot,
                        species=tuple(s for s, k in zip(snapshot.species, keep, strict=True) if k),
                        positions=snapshot.positions[keep],
                        label=f"{snapshot.label} with one vacancy",
                    )
        except ValueError as error:
            raise InvalidInputError(
                f"The specimen could not be built: {error}",
                field="zone_axis",
                hint=(
                    "Give a non-zero, low-index zone axis and a positive thickness. A zone axis "
                    "with no perpendicular lattice vectors cannot be put in a periodic box."
                ),
            ) from error
        info["unit_cell_repeats"] = list(repeats)
        info["slab_thickness_angstrom"] = slab

    if snapshot.natoms > MAX_SPECIMEN_ATOMS:
        imported = sample_type == "imported"
        raise InvalidInputError(
            f"The specimen holds {snapshot.natoms} atoms, more than the {MAX_SPECIMEN_ATOMS} the "
            "workbench simulates.",
            field="structure_file" if imported else "thickness_angstrom",
            hint=(
                "Crop the structure to a smaller region before exporting it."
                if imported
                else "Reduce the thickness or the supercell size."
            ),
        )
    info["atom_count"] = snapshot.natoms
    return snapshot, info


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


_SIMULATE_PARAMETERS: tuple[Parameter, ...] = (
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
            (
                "imported",
                "Imported structure (.xyz)",
                "Atomic coordinates from a simulation, opened as an .xyz file in the rail",
            ),
        ),
        default="crystalline",
        row="specimen",
        field_width="medium",
        help_text=(
            "Specimen microstructure morphology. Choose the imported structure to simulate "
            "atomic coordinates from a molecular-dynamics, DFT or structure-building code, "
            "opened as an .xyz file; the crystal phase then serves only the other specimen "
            "types."
        ),
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
        row="specimen_size",
        field_width="tiny",
        help_text="Number of unit-cell repeats along transverse x and y directions.",
    ),
    NumberParameter(
        name="thickness_angstrom",
        label="Specimen thickness",
        default=10.0,
        minimum=2.0,
        maximum=500.0,
        units="Å",
        symbol="foil_thickness",
        row="specimen_size",
        field_width="short",
        help_text=(
            "Thickness of the specimen along the beam. A crystal slab is built from whole "
            "unit-cell repeats along the lattice vector closest to the beam, so the "
            "thickness simulated is the nearest whole number of repeats at or above this "
            "value, and the result reports it. For the amorphous foil it is the foil depth. "
            "An imported structure keeps the thickness of its own coordinates. In the "
            "multislice engine a thicker specimen propagates the wave through more slices; "
            "the pure-Python phase-object fallback projects every atom into one plane, so "
            "there thickness only strengthens the projected potential."
        ),
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
        row="sampling",
        field_width="short",
        help_text=(
            "Real-space pixel pitch. The multislice band limit is a third of its inverse, so "
            "0.2 Å represents scattering to about 40 mrad at 200 kV and 0.05 Å to about 170 mrad; "
            "a strongly scattering specimen needs the finer pitch."
        ),
    ),
)

_STRUCTURE_FILE_PARAMETER = ObjectParameter(
    name="structure_file",
    label="Structure file",
    help_text=(
        "An atomic structure written by a simulation, as .xyz or extended .xyz. Each "
        "atom line is an element symbol and Cartesian x, y, z in ångströms, with the "
        "beam along +z. An extended-XYZ comment line may carry "
        "`Lattice=\"ax ay az bx by bz cx cy cz\"`, which sets an orthogonal periodic box, "
        "and `Properties=species:S:1:pos:R:3:...`, which locates the columns. Without a "
        "lattice a box is fitted around the atoms with a 1 Å margin. Opened through "
        "**Open a structure file** in the workbench rail; used when the sample morphology "
        "is the imported structure."
    ),
    required=False,
)

#: A workbench multislice is refused beyond this many lateral grid points ...
MAX_WORKBENCH_GRID_POINTS = 1024 * 1024
#: ... and beyond this much work, counted as grid points x slices x frozen-phonon
#: configurations: about a minute of FFTs, after which a request is better run
#: from Python than from a browser tab waiting on it.
MAX_WORKBENCH_WORK = 400_000_000

_ENGINE_LABELS = {
    "multislice": "PyTex multislice",
    "phase_object": "single phase object",
    "abtem": "abTEM multislice",
}

#: The multislice controls shared by the micrograph and the series operations.
#: `slice_thickness_angstrom` follows `sampling_angstrom` so the two share a row.
_MULTISLICE_PARAMETERS: tuple[Parameter, ...] = (
    NumberParameter(
        name="slice_thickness_angstrom",
        label="Slice thickness",
        default=1.5,
        minimum=0.2,
        maximum=10.0,
        units="Å",
        symbol="slice_thickness",
        row="sampling",
        field_width="short",
        help_text=(
            "Thickness Δz of each multislice slice along the beam. The splitting error grows "
            "with it, so 1-2 Å is usual; a value that divides the crystal period along the beam "
            "makes every period's slices identical, and they are then built only once."
        ),
    ),
    ChoiceParameter(
        name="engine",
        label="Simulation engine",
        options=(
            (
                "multislice",
                "PyTex multislice",
                "abTEM's algorithm in PyTex: dynamical scattering through every slice",
            ),
            (
                "phase_object",
                "Single phase object",
                "Every atom projected into one plane; valid only for a very thin, weak specimen",
            ),
            ("abtem", "abTEM", "abTEM itself, when it is installed on this server"),
        ),
        default="multislice",
        row="engine",
        field_width="medium",
        help_text=(
            "How the exit wave is computed. The PyTex multislice propagates the wave through the "
            "specimen slice by slice, so thickness, channelling and dynamical scattering are all "
            "present; the phase object projects the whole specimen into one plane."
        ),
    ),
    ChoiceParameter(
        name="parametrization",
        label="Scattering factors",
        options=(
            ("lobato", "Lobato & Van Dyck", "Lobato & Van Dyck (2014), abTEM's default"),
            ("kirkland", "Kirkland", "Kirkland (2010), Appendix C"),
            ("mott_bethe", "Mott-Bethe", "PyTex's X-ray table through the Mott-Bethe relation"),
        ),
        default="lobato",
        row="engine",
        field_width="medium",
        help_text=(
            "The independent-atom electron scattering factor the projected potential is built "
            "from. The three agree within a few percent; Mott-Bethe is the one the Bloch-wave "
            "(CBED) solver uses."
        ),
    ),
    ChoiceParameter(
        name="temporal_coherence",
        label="Focal spread treatment",
        options=(
            (
                "quasi_coherent",
                "Frank envelope",
                "Damp the wave transfer by Frank's envelope: exact for linear image terms only",
            ),
            (
                "focal_integration",
                "Focal integration",
                "Average images over the defocus spread: exact, slower",
            ),
        ),
        default="quasi_coherent",
        advanced=True,
        field_width="medium",
        help_text=(
            "The focal spread can damp the transferred wave (Frank's envelope, fast, exact for a "
            "weak object) or be integrated exactly as an average of images over defocus, which "
            "keeps the non-linear interference of equivalent beams that the envelope suppresses."
        ),
    ),
    NumberParameter(
        name="tilt_x_mrad",
        label="Beam tilt along x",
        default=0.0,
        minimum=-50.0,
        maximum=50.0,
        units="mrad",
        row="tilt",
        advanced=True,
        field_width="short",
        help_text="Tilt of the incident beam from the zone axis towards +x, in mrad.",
    ),
    NumberParameter(
        name="tilt_y_mrad",
        label="Beam tilt along y",
        default=0.0,
        minimum=-50.0,
        maximum=50.0,
        units="mrad",
        row="tilt",
        advanced=True,
        field_width="short",
        help_text="Tilt of the incident beam from the zone axis towards +y, in mrad.",
    ),
    IntegerParameter(
        name="frozen_phonons",
        label="Frozen-phonon configurations",
        default=1,
        minimum=1,
        maximum=32,
        row="phonons",
        advanced=True,
        field_width="tiny",
        help_text=(
            "Number of thermally displaced configurations whose images are averaged. One means "
            "a static lattice; eight or more give thermal diffuse scattering."
        ),
    ),
    NumberParameter(
        name="thermal_rms_angstrom",
        label="Thermal RMS displacement",
        default=0.08,
        minimum=0.0,
        maximum=0.3,
        units="Å",
        symbol="thermal_displacement",
        row="phonons",
        advanced=True,
        field_width="short",
        help_text=(
            "One-axis RMS displacement of every atom in the frozen-phonon ensemble, "
            "sqrt(B / 8π²); about 0.076 Å for silicon at room temperature."
        ),
    ),
)


def _multislice_options(request: Mapping[str, Any]) -> dict[str, Any]:
    configurations = int(request.get("frozen_phonons", 1))
    rms = float(request.get("thermal_rms_angstrom", 0.08))
    thermal = configurations > 1 and rms > 0.0
    return {
        "parametrization": PotentialParametrization(str(request.get("parametrization", "lobato"))),
        "tilt_mrad": (
            float(request.get("tilt_x_mrad", 0.0)),
            float(request.get("tilt_y_mrad", 0.0)),
        ),
        "frozen_phonon_sigma_angstrom": rms if thermal else None,
        "frozen_phonon_configurations": configurations if thermal else 1,
        "seed": 0,
    }


def _run_multislice(
    snapshot: AtomicSnapshot,
    energy_kev: float,
    request: Mapping[str, Any],
    exit_depths_angstrom: Any = None,
) -> MultisliceExitWave:
    """Run the PyTex multislice for a workbench request, refusing what would not finish."""

    sampling = float(request.get("sampling_angstrom", 0.2))
    slice_thickness = float(request.get("slice_thickness_angstrom", 1.5))
    options = _multislice_options(request)
    grid = MultisliceGrid.from_sampling(
        (float(snapshot.cell[0, 0]), float(snapshot.cell[1, 1])), sampling
    )
    ny, nx = grid.shape
    if ny * nx > MAX_WORKBENCH_GRID_POINTS:
        raise InvalidInputError(
            f"A {nx} × {ny} grid is more than the workbench computes.",
            field="sampling_angstrom",
            hint="Coarsen the sampling or reduce the supercell size.",
        )
    slices = math.ceil(float(snapshot.cell[2, 2]) / slice_thickness)
    work = ny * nx * slices * int(options["frozen_phonon_configurations"])
    if work > MAX_WORKBENCH_WORK:
        raise InvalidInputError(
            f"{slices} slices on a {nx} × {ny} grid"
            + (
                f" for {options['frozen_phonon_configurations']} configurations"
                if options["frozen_phonon_configurations"] > 1
                else ""
            )
            + " is more than the workbench computes.",
            field="slice_thickness_angstrom",
            hint=(
                "Use thicker slices, fewer frozen-phonon configurations, a coarser sampling, "
                "a smaller supercell or a thinner specimen."
            ),
        )
    try:
        return multislice(
            snapshot,
            energy_kev,
            grid=grid,
            slice_thickness_angstrom=slice_thickness,
            exit_depths_angstrom=exit_depths_angstrom,
            **options,
        )
    except ValueError as error:
        raise InvalidInputError(
            f"The multislice could not run: {error}",
            field="parametrization",
            hint="Every species needs tabulated scattering factors in the chosen table.",
        ) from error


def _engine_note(engine: HREMEngine, exit_wave: MultisliceExitWave | None) -> str:
    if engine is HREMEngine.MULTISLICE and exit_wave is not None:
        return (
            "Simulated by PyTex multislice (the algorithm of abTEM, validated against it and "
            f"against Bloch waves): {exit_wave.parametrization.value} potentials projected into "
            f"{exit_wave.num_slices} slices of {exit_wave.slice_thickness_angstrom:.2f} Å, "
            "band-limited transmission and Fresnel propagation."
        )
    if engine is HREMEngine.ABTEM:
        return "Simulated by abTEM multislice through the PyTex adapter."
    return (
        "Simulated as a single phase object: every atom is projected into one plane, so no "
        "propagation of the wave within the specimen is modelled. Valid only for a very thin, "
        "weakly scattering specimen; choose the PyTex multislice for anything thicker."
    )


def _multislice_rows(exit_wave: MultisliceExitWave | None) -> tuple[dict[str, str], ...]:
    if exit_wave is None:
        return ()
    ny, nx = exit_wave.grid.shape
    rows = [
        {
            "metric": "Slices (distinct potentials)",
            "value": f"{exit_wave.num_slices} ({exit_wave.unique_slices})",
            "units": "",
        },
        {
            "metric": "Slice thickness Δz",
            "value": f"{exit_wave.slice_thickness_angstrom:.3f}",
            "units": "Å",
        },
        {
            "metric": "Band limit (largest scattering angle)",
            "value": (
                f"{exit_wave.grid.antialias_cutoff_inv_angstrom:.2f} Å⁻¹ "
                f"({exit_wave.grid.max_scattering_angle_mrad(exit_wave.energy_kev):.0f} mrad)"
            ),
            "units": "",
        },
        {
            "metric": "Retained intensity at the exit surface",
            "value": f"{float(exit_wave.retained_intensity[-1]):.4f}",
            "units": "",
        },
        {"metric": "Scattering factors", "value": exit_wave.parametrization.value, "units": ""},
    ]
    if exit_wave.configurations > 1:
        rows.append(
            {
                "metric": "Frozen-phonon configurations",
                "value": str(exit_wave.configurations),
                "units": "",
            }
        )
    if any(exit_wave.tilt_mrad):
        rows.append(
            {
                "metric": "Beam tilt (x, y)",
                "value": f"{exit_wave.tilt_mrad[0]:.2f}, {exit_wave.tilt_mrad[1]:.2f}",
                "units": "mrad",
            }
        )
    del ny, nx
    return tuple(rows)


def _sampling_note(exit_wave: MultisliceExitWave | None) -> tuple[str, ...]:
    if exit_wave is None or float(exit_wave.retained_intensity[-1]) >= 0.95:
        return ()
    return (
        f"Only {float(exit_wave.retained_intensity[-1]) * 100:.1f} % of the beam intensity "
        "remains inside the band limit at the exit surface; the rest was scattered beyond "
        f"{exit_wave.grid.max_scattering_angle_mrad(exit_wave.energy_kev):.0f} mrad. Refine the "
        "sampling before trusting fine detail.",
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
        *_SIMULATE_PARAMETERS,
        *_MULTISLICE_PARAMETERS,
        _STRUCTURE_FILE_PARAMETER,
    ),
)
def _simulate_hrem(request: dict[str, Any]) -> dict[str, Any]:
    _phase_spec, phase = phase_from_request(request["phase"])
    voltage = float(request.get("beam_energy_kev", 200.0))
    mode_str = str(request.get("mode", "double_corrected"))
    defocus = float(request.get("defocus_angstrom", -30.0))
    cs_um = float(request.get("cs_um", 0.0))
    sampling = float(request.get("sampling_angstrom", 0.2))

    snap, specimen = _specimen_from_request(request, phase)
    imported = specimen["sample_type"] == "imported"
    slab_thickness = float(specimen["slab_thickness_angstrom"])
    repeats = specimen["unit_cell_repeats"]

    mode = DoubleCorrectionMode(mode_str)
    aberr = _aberrations_from_request(request, voltage, defocus, cs_um, mode)

    engine = HREMEngine(str(request.get("engine", "multislice")))
    coherence = TemporalCoherence(str(request.get("temporal_coherence", "quasi_coherent")))
    exit_wave: MultisliceExitWave | None = None
    if engine is HREMEngine.MULTISLICE:
        exit_wave = _run_multislice(snap, voltage, request)
        result = exit_wave.hrem_result(aberr, temporal_coherence=coherence)
    elif engine is HREMEngine.ABTEM:
        if not is_abtem_available():
            raise InvalidInputError(
                "abTEM is not installed on this server.",
                field="engine",
                hint="Choose the PyTex multislice, which implements the same algorithm.",
            )
        result = simulate_hrem(
            snap,
            aberr,
            sampling_angstrom=sampling,
            slice_thickness_angstrom=float(request.get("slice_thickness_angstrom", 1.5)),
            engine=HREMEngine.ABTEM,
        )
    else:
        result = simulate_hrem(
            snap, aberr, sampling_angstrom=sampling, engine=HREMEngine.PHASE_OBJECT
        )

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
                "metric": "Specimen source",
                "value": (
                    specimen["structure_file"]["name"] if imported else str(specimen["sample_type"])
                ),
                "units": "",
            },
            {
                "metric": (
                    "Specimen thickness (atom-centre span along the beam)"
                    if imported
                    else "Specimen thickness along the beam"
                ),
                "value": f"{slab_thickness:.2f}",
                "units": "Å",
            },
            *(
                ()
                if repeats is None
                else (
                    {
                        "metric": (
                            "Unit-cell repeats along a × b × c"
                            if specimen["zone_axis_cell"] is None
                            else "Periodic-cell repeats across × across × along the beam"
                        ),
                        "value": " × ".join(str(value) for value in repeats),
                        "units": "",
                    },
                )
            ),
            {"metric": "Atoms in the specimen", "value": str(snap.natoms), "units": ""},
            {"metric": "Simulation engine", "value": _ENGINE_LABELS[engine.value], "units": ""},
            *_multislice_rows(exit_wave),
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

    engine_note = _engine_note(engine, exit_wave)

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
    if imported:
        summary_text += (
            f" The specimen is the {snap.natoms}-atom structure read from "
            f"{specimen['structure_file']['name']}, whose atom centres span "
            f"{slab_thickness:.1f} Å along the beam; the thickness control does not apply to it."
        )
    elif repeats is not None:
        summary_text += (
            f" The slab is {slab_thickness:.1f} Å thick along the beam "
            f"({' × '.join(str(value) for value in repeats)} "
            f"{'unit' if specimen['zone_axis_cell'] is None else 'periodic'} cells, "
            f"{snap.natoms} atoms), "
            f"the nearest whole number of repeats at or above the "
            f"{specimen['requested_thickness_angstrom']:.1f} Å requested."
        )
    else:
        summary_text += (
            f" The amorphous foil is {slab_thickness:.1f} Å thick ({snap.natoms} atoms)."
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

    image_rows, image_columns = result.image.shape
    notes: list[str] = [engine_note, *_sampling_note(exit_wave)]
    if specimen["zone_axis_cell"]:
        notes.append(str(specimen["zone_axis_cell"]))
    inputs = dict(request)
    if imported:
        # The provenance of the file, not its text: a saved result names what was
        # read without carrying a second copy of a possibly large structure.
        inputs["structure_file"] = specimen["structure_file"]

    app_res = AppResult(
        title=f"HRTEM Simulation: {snap.label}",
        summary=summary_text,
        table=table,
        data={
            # One PNG pixel per simulation pixel, row 0 at the bottom, so zooming
            # in the viewer shows the real sampling rather than a resampling.
            "image_png": result.to_png_base64(native_resolution=True),
            "power_spectrum_png": result.to_power_spectrum_base64(native_resolution=True),
            "image_shape_px": [int(image_rows), int(image_columns)],
            "nyquist_inv_angstrom": [
                0.5 * image_columns / float(result.extent_angstrom[0]),
                0.5 * image_rows / float(result.extent_angstrom[1]),
            ],
            "specimen_thickness_angstrom": slab_thickness,
            "requested_thickness_angstrom": float(specimen["requested_thickness_angstrom"]),
            "unit_cell_repeats": repeats,
            "structure_file": specimen["structure_file"],
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
            "engine": engine.value,
            "multislice": None if exit_wave is None else multislice_summary(exit_wave),
            "description": result.describe()
            + ("" if exit_wave is None else " " + exit_wave.describe()),
        },
        inputs=inputs,
        notes=tuple(notes),
        citations=(
            "Kirkland (2010), Advanced Computing in Electron Microscopy, 2nd ed.",
            "Cowley & Moodie (1957), Acta Crystallogr. 10, 609-619.",
            "Madsen et al. (2021), abTEM: An open-source framework, ChemPhysChem 22, 1-13.",
        ),
    )
    app_res = replace(
        app_res,
        figures=(
            hrem_spectrum_figure(
                power_spectrum=np.asarray(result.power_spectrum, dtype=float),
                pixel_size_angstrom=float(result.pixel_size_angstrom),
                frequencies=np.asarray(result.ctf.spatial_frequencies_inv_angstrom, dtype=float),
                transfer=np.asarray(result.ctf.transfer_function, dtype=float),
                envelope=np.asarray(result.ctf.total_envelope, dtype=float),
                point_resolution=result.point_resolution_angstrom,
                information_limit=result.information_limit_angstrom,
            ),
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


#: A series is refused beyond this many defoci or this many images in all.
MAX_SERIES_DEFOCI = 25
MAX_SERIES_IMAGES = 96

#: Depths at which the beam intensities are sampled for the thickness plot.
_BEAM_CURVE_DEPTHS = 48

_SERIES_PARAMETERS: tuple[Parameter, ...] = (
    NumberParameter(
        name="defocus_start_angstrom",
        label="Defocus from",
        default=-150.0,
        minimum=-5000.0,
        maximum=5000.0,
        units="Å",
        row="focal_range",
        field_width="short",
        help_text="First defocus of the series; negative is underfocus.",
    ),
    NumberParameter(
        name="defocus_stop_angstrom",
        label="Defocus to",
        default=150.0,
        minimum=-5000.0,
        maximum=5000.0,
        units="Å",
        row="focal_range",
        field_width="short",
        help_text="Last defocus of the series, at or above the first.",
    ),
    NumberParameter(
        name="defocus_step_angstrom",
        label="Defocus step",
        default=50.0,
        minimum=1.0,
        maximum=2000.0,
        units="Å",
        row="focal_range",
        field_width="short",
        help_text=f"Spacing of the defoci; at most {MAX_SERIES_DEFOCI} images per thickness.",
    ),
    IntegerParameter(
        name="thickness_steps",
        label="Thickness rows",
        default=3,
        minimum=1,
        maximum=8,
        field_width="tiny",
        help_text=(
            "Number of thicknesses in the tableau, evenly spaced up to the specimen thickness. "
            "One gives a single focal series at the full thickness."
        ),
    ),
)


def _beam_curves(
    exit_wave: MultisliceExitWave,
    specimen: Mapping[str, Any],
    spec: Any,
) -> tuple[list[dict[str, Any]], list[tuple[str, list[float]]]]:
    """The transmitted beam and the four strongest diffracted beams against thickness.

    Only a crystal in a periodic zone-axis box has Bragg beams on the grid, at
    Fourier indices that are multiples of the lateral repeats. Each is labelled
    by its Miller indices, recovered through the box's rotation.
    """

    cell = specimen.get("_zone_axis_cell")
    repeats = specimen.get("unit_cell_repeats")
    if cell is None or repeats is None:
        return [], []
    m1, m2 = int(repeats[0]), int(repeats[1])
    ny, nx = exit_wave.grid.shape
    spectra = np.fft.fft2(exit_wave.waves, norm="forward")
    intensity = np.mean(np.abs(spectra) ** 2, axis=0)  # (depths, ny, nx)
    rows = np.arange(0, ny, m2)
    cols = np.arange(0, nx, m1)
    lattice = intensity[:, rows][:, :, cols]
    strongest = np.max(lattice, axis=0)
    strongest[0, 0] = -1.0
    picks = [(0, 0)] + [
        tuple(int(v) for v in np.unravel_index(int(i), strongest.shape))
        for i in np.argsort(strongest, axis=None)[::-1][:4]
    ]
    t1, t2, _ = cell.lengths_angstrom
    basis = spec.to_phase().lattice.direct_basis()
    direct = np.vstack([basis.vector(0), basis.vector(1), basis.vector(2)])
    beams: list[dict[str, Any]] = []
    curves: list[tuple[str, list[float]]] = []
    for r, c in picks:
        i = c if c <= len(cols) // 2 else c - len(cols)
        j = r if r <= len(rows) // 2 else r - len(rows)
        g_box = np.array([i / t1, j / t2, 0.0])
        hkl = np.rint(direct @ (cell.rotation.T @ g_box)).astype(int)
        plain = plane_label(hkl, spec=spec, style="plain")
        values = [float(v) for v in lattice[:, r, c]]
        beams.append({"label": plain, "intensity": values})
        curves.append((plane_label(hkl, spec=spec, style="mathtext"), values))
    return beams, curves


@REGISTRY.operation(
    "tem.hrtem_series",
    title="HRTEM focal and thickness series",
    summary=(
        "A focal series, or a defocus-thickness tableau, of HRTEM images from one multislice "
        "run, with the beam intensities against thickness."
    ),
    help_text=(
        "Propagate the wave once through the specimen with the PyTex multislice, keep the exit "
        "wave at several thicknesses, and image each through the objective lens at every "
        "defocus of the series. The exit wave does not depend on the lens, so a focal series "
        "costs one multislice run. One thickness row gives a focal series; several give the "
        "defocus-thickness tableau used to match experimental images."
    ),
    panel="tem_hrem",
    documentation=DocumentationLink("Multislice HRTEM", "algorithms/multislice_hrtem"),
    returns=(
        "Images per thickness and defocus as PNGs, their RMS contrast, the beam intensities "
        "against thickness for a crystal, and the multislice diagnostics."
    ),
    parameters=(
        *(
            parameter
            for parameter in _SIMULATE_PARAMETERS
            if parameter.name not in {"defocus_angstrom"}
        ),
        # A series always runs the PyTex multislice, so the engine choice is
        # absent and the scattering factors no longer share its row.
        *(
            replace(parameter, row=None) if parameter.name == "parametrization" else parameter
            for parameter in _MULTISLICE_PARAMETERS
            if parameter.name != "engine"
        ),
        *_SERIES_PARAMETERS,
        _STRUCTURE_FILE_PARAMETER,
    ),
)
def _hrtem_series(request: dict[str, Any]) -> dict[str, Any]:
    spec, phase = phase_from_request(request["phase"])
    voltage = float(request.get("beam_energy_kev", 200.0))
    mode = DoubleCorrectionMode(str(request.get("mode", "double_corrected")))
    cs_um = float(request.get("cs_um", 0.0))
    start = float(request.get("defocus_start_angstrom", -150.0))
    stop = float(request.get("defocus_stop_angstrom", 150.0))
    step = float(request.get("defocus_step_angstrom", 50.0))
    if stop < start:
        raise InvalidInputError(
            "The last defocus is below the first.",
            field="defocus_stop_angstrom",
            hint="Give the series from its most negative (underfocus) value upwards.",
        )
    defoci = np.arange(start, stop + 0.5 * step, step)
    if defoci.size > MAX_SERIES_DEFOCI:
        raise InvalidInputError(
            f"The series asks for {defoci.size} defoci, more than {MAX_SERIES_DEFOCI}.",
            field="defocus_step_angstrom",
            hint="Increase the step or narrow the range.",
        )
    rows = int(request.get("thickness_steps", 3))
    if rows * defoci.size > MAX_SERIES_IMAGES:
        raise InvalidInputError(
            f"{rows} thicknesses by {defoci.size} defoci is more than {MAX_SERIES_IMAGES} images.",
            field="thickness_steps",
            hint="Use fewer thickness rows or fewer defoci.",
        )

    snap, specimen = _specimen_from_request(request, phase)
    aberr = _aberrations_from_request(request, voltage, 0.0, cs_um, mode)
    coherence = TemporalCoherence(str(request.get("temporal_coherence", "quasi_coherent")))
    depth = float(snap.cell[2, 2])
    tableau_depths = depth * np.arange(1, rows + 1) / rows
    curve_depths = depth * np.arange(1, _BEAM_CURVE_DEPTHS + 1) / _BEAM_CURVE_DEPTHS
    exit_wave = _run_multislice(
        snap, voltage, request, exit_depths_angstrom=np.concatenate([tableau_depths, curve_depths])
    )
    indices = sorted({exit_wave.depth_index(float(t)) for t in tableau_depths})
    thicknesses = [float(exit_wave.depths_angstrom[i]) for i in indices]
    series = [
        exit_wave.focal_series(aberr, defoci, depth_index=i, temporal_coherence=coherence)
        for i in indices
    ]
    images = np.stack([s.images for s in series])
    contrasts = np.array([s.contrasts for s in series])
    beams, curves = _beam_curves(exit_wave, specimen, spec)

    tiles = [
        {
            "thickness_angstrom": thicknesses[i],
            "defocus_angstrom": float(defoci[j]),
            "contrast": float(contrasts[i, j]),
            "png": _array_png_data_url(images[i, j], "gray"),
        }
        for i in range(len(thicknesses))
        for j in range(defoci.size)
    ]
    best = np.unravel_index(int(np.argmax(contrasts)), contrasts.shape)
    kind = "Focal series" if len(thicknesses) == 1 else "Defocus-thickness tableau"
    summary = (
        f"{kind} of {snap.label}: {len(thicknesses)} thickness"
        f"{'es' if len(thicknesses) != 1 else ''} "
        f"({', '.join(f'{t:.1f}' for t in thicknesses)} Å) by {defoci.size} defoci from "
        f"{defoci[0]:.0f} to {defoci[-1]:.0f} Å at {aberr.energy_kev:.0f} kV "
        f"(Cs = {aberr.cs_um:.1f} µm), all from one multislice run of {exit_wave.num_slices} "
        f"slices. The strongest RMS contrast, {contrasts[best] * 100:.1f} %, is at "
        f"t = {thicknesses[best[0]]:.1f} Å and Δf = {defoci[best[1]]:.0f} Å."
    )
    table = ResultTable(
        columns=(
            Column("thickness", "Thickness", units="Å", help_text="Exit-surface depth"),
            Column("defocus", "Defocus", units="Å", help_text="Objective defocus"),
            Column("contrast", "RMS contrast", units="%", help_text="Standard deviation / mean"),
            Column("mean", "Mean intensity", help_text="Mean of the image"),
        ),
        rows=tuple(
            {
                "thickness": f"{thicknesses[i]:.1f}",
                "defocus": f"{defoci[j]:.0f}",
                "contrast": f"{contrasts[i, j] * 100:.2f}",
                "mean": f"{float(np.mean(images[i, j])):.4f}",
            }
            for i in range(len(thicknesses))
            for j in range(defoci.size)
        ),
        caption="Every image of the series and its contrast",
    )
    figures = [
        hrtem_tableau_figure(
            images=images,
            defoci_angstrom=[float(d) for d in defoci],
            thicknesses_angstrom=thicknesses,
            extent_angstrom=exit_wave.grid.extent_angstrom,
        ),
        focal_contrast_figure(
            defoci_angstrom=[float(d) for d in defoci],
            thicknesses_angstrom=thicknesses,
            contrasts=contrasts,
        ),
    ]
    if curves:
        figures.append(
            beam_thickness_figure(
                thickness_angstrom=[float(t) for t in exit_wave.depths_angstrom],
                beams=curves,
            )
        )
    notes = [
        _engine_note(HREMEngine.MULTISLICE, exit_wave),
        (
            "Every image shares one exit wave, so the series differs only through the objective "
            "lens; each tile is scaled to its own grey range."
        ),
        *_sampling_note(exit_wave),
    ]
    if specimen.get("zone_axis_cell"):
        notes.append(str(specimen["zone_axis_cell"]))
    ny, nx = exit_wave.grid.shape
    app_res = AppResult(
        title=f"HRTEM {kind.lower()}: {snap.label}",
        summary=summary,
        table=table,
        data={
            "tiles": tiles,
            "defoci_angstrom": [float(d) for d in defoci],
            "thicknesses_angstrom": thicknesses,
            "contrasts": contrasts.tolist(),
            "extent_angstrom": list(exit_wave.grid.extent_angstrom),
            "image_shape_px": [int(ny), int(nx)],
            "pixel_size_angstrom": float(exit_wave.grid.sampling_angstrom[0]),
            "beam_thickness_angstrom": [float(t) for t in exit_wave.depths_angstrom],
            "beams": beams,
            "multislice": multislice_summary(exit_wave),
            "unit_cell_repeats": specimen["unit_cell_repeats"],
            "natoms": snap.natoms,
            "sample_label": snap.label,
            "description": exit_wave.describe() + " " + series[-1].describe(),
        },
        inputs=_inputs_for_report(request, specimen),
        notes=tuple(notes),
        citations=(
            "Cowley & Moodie (1957), Acta Crystallogr. 10, 609-619.",
            "Kirkland (2010), Advanced Computing in Electron Microscopy, 2nd ed.",
            "Madsen & Susi (2021), The abTEM code, Open Research Europe 1, 24.",
            "Coene et al. (1992), Phys. Rev. Lett. 69, 3743-3746.",
            "O'Keefe & Kilaas (1988), Scanning Microscopy Suppl. 2, 225-244.",
        ),
    )
    return replace(app_res, figures=tuple(figures)).to_json()


def _inputs_for_report(request: Mapping[str, Any], specimen: Mapping[str, Any]) -> dict[str, Any]:
    inputs = dict(request)
    if specimen["sample_type"] == "imported":
        # The provenance of the file, not its text: a saved result names what was
        # read without carrying a second copy of a possibly large structure.
        inputs["structure_file"] = specimen["structure_file"]
    return inputs


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
            id="hrem.example.si_focal_series",
            title="Through-focus series of silicon [110]",
            panel="tem_hrem",
            summary=(
                "A Cs-corrected 300 kV focal series of a 40 Å silicon [110] slab, nine images "
                "from one multislice run."
            ),
            teaches=(
                "Notice that the dumbbells change from dark to bright and back as the defocus "
                "passes through zero, and that the contrast curve has its minimum close to "
                "Gaussian focus. Every image comes from the same exit wave: defocusing the lens "
                "is the same as letting the exit wave propagate through vacuum, which is why a "
                "focal series costs one multislice run and why an experimental focal series can "
                "be inverted to recover the exit wave."
            ),
            operation="tem.hrtem_series",
            request={
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 2,
                "thickness_angstrom": 40.0,
                "beam_energy_kev": 300.0,
                "mode": "cs_corrected",
                "cs_um": 5.0,
                "sampling_angstrom": 0.1,
                "defocus_start_angstrom": -200.0,
                "defocus_stop_angstrom": 200.0,
                "defocus_step_angstrom": 50.0,
                "thickness_steps": 1,
            },
        ),
        ExampleScenario(
            id="hrem.example.si_defocus_thickness_map",
            title="Defocus-thickness tableau of silicon [110]",
            panel="tem_hrem",
            summary=(
                "Four thicknesses up to 160 Å by five defoci, with the beam intensities against "
                "thickness, for a Cs-corrected 200 kV lens."
            ),
            teaches=(
                "Notice that image contrast depends on thickness as strongly as on focus: down a "
                "column of the tableau the dumbbells fade and reverse as the transmitted beam "
                "hands its intensity to the diffracted beams and takes it back (the beam plot). "
                "This dynamical exchange is why a single HRTEM image cannot be read as a "
                "projection of the atoms, and why experimental images are matched against a "
                "tableau like this one."
            ),
            operation="tem.hrtem_series",
            request={
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 1,
                "thickness_angstrom": 160.0,
                "beam_energy_kev": 200.0,
                "mode": "cs_corrected",
                "cs_um": 10.0,
                "sampling_angstrom": 0.1,
                "slice_thickness_angstrom": 1.92,
                "defocus_start_angstrom": -100.0,
                "defocus_stop_angstrom": 100.0,
                "defocus_step_angstrom": 50.0,
                "thickness_steps": 4,
            },
        ),
        ExampleScenario(
            id="hrem.example.si_mistilt",
            title="A crystal tilted 10 mrad off its zone axis",
            panel="tem_hrem",
            summary=(
                "Silicon [110], 80 Å thick, with the beam tilted 10 mrad from the zone axis "
                "towards the dumbbell axis."
            ),
            teaches=(
                "Notice that the dumbbells smear and lose their symmetry although the lens is "
                "unchanged. The tilted beam runs along the columns only near the entrance "
                "surface and drifts off them with depth, so the exit wave records each column "
                "sheared by a distance proportional to the thickness. A few mrad of mistilt is "
                "the most common reason a thick-specimen image refuses to match a simulation "
                "made exactly on the zone axis."
            ),
            operation="tem.simulate_hrem",
            request={
                "phase": {"builtin": "si_diamond"},
                "sample_type": "crystalline",
                "zone_axis": [1, 1, 0],
                "supercell_xy": 2,
                "thickness_angstrom": 80.0,
                "beam_energy_kev": 300.0,
                "mode": "cs_corrected",
                "defocus_angstrom": -40.0,
                "cs_um": 5.0,
                "sampling_angstrom": 0.1,
                "tilt_x_mrad": 10.0,
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
