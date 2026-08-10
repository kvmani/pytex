"""Worked examples: many-beam coupling, absorption, HOLZ lines, diffraction groups.

Four claims that a dynamical CBED implementation can get plausibly wrong and that
each have a check with independent provenance:

- the many-beam solver must reduce to the published two-beam closed form when
  there are two beams, which pins the diagonal convention, the off-diagonal
  scale and the propagator factor at once;
- without absorption the propagator is unitary, so the beam intensities sum to
  exactly one -- an analytic identity, not a fitted result;
- a HOLZ line offset is degenerate between lattice strain and wavelength, which
  is an exact cancellation in the closed form and the reason voltage calibration
  is mandatory;
- and the diffraction-group construction must produce Buxton's 31 groups, a
  published count that PyTex derives rather than stores.

See :doc:`../../concepts/diffraction_foundation` for the surrounding doctrine.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

DYNAMICAL_SETUP = """
import numpy as np
from pytex import (
    AbsorptionModel,
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    beam_set_for_zone,
    beam_set_from_indices,
    extinction_distance_angstrom,
    holz_line_pattern,
    solve_bloch_waves,
    two_beam_rocking_curve,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
nickel_lattice = Lattice(3.5239, 3.5239, 3.5239, 90.0, 90.0, 90.0, crystal_frame=crystal)
nickel = Phase(
    name="nickel-fcc",
    lattice=nickel_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=nickel_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Ni{index}",
                species="Ni",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
wavelength = electron_wavelength_angstrom(200.0)
"""

GROUPS_SETUP = """
from pytex import (
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for,
    diffraction_group_symbols,
)
"""

_EXTINCTION = SymbolUse(
    r"\xi_{g}",
    "Two-beam extinction distance of reflection g; the depth period of the "
    "intensity exchange between the transmitted and diffracted beams.",
)
_DEVIATION = SymbolUse(
    r"s_{g}",
    "Excitation error: deviation of reflection g from the exact Bragg condition.",
)
_COEFFICIENT = SymbolUse(
    r"\nu_{g}",
    "Complex Fourier coefficient of the scaled lattice potential; the "
    "off-diagonal element of the dynamical structure matrix, with "
    "|nu_g| = 1 / xi_g.",
)
_THICKNESS = SymbolUse(r"t", "Foil thickness along the beam.")
_WAVELENGTH = SymbolUse(r"\lambda", "Radiation wavelength.")

_DIFFRACTION_CONCEPT = SeeAlso("Diffraction foundation", "../../concepts/diffraction_foundation")
_API = SeeAlso("Diffraction API", "../../api/index")


TWO_BEAM_LIMIT = WorkedExample(
    id="diffraction-dynamical-two-beam-limit-of-the-many-beam-solver",
    title="Two beams reduce the Bloch-wave solver to the closed form exactly",
    domain="diffraction",
    scenario=(
        "A many-beam dynamical calculation has no independent standard to be "
        "checked against, so the one calibration available is its own limiting "
        "case. Restricted to a single reflection, the coupled system has the "
        "closed-form solution I_g = sin^2(pi t s_eff) / (xi_g s_eff)^2 with "
        "s_eff^2 = s^2 + xi_g^-2. Reproducing it to machine precision pins "
        "three conventions simultaneously: the diagonal 2 s_g, the off-diagonal "
        "scale |nu_g| = 1 / xi_g, and the factor i pi in the propagator. Any "
        "one of them wrong yields a rocking curve of the right general shape "
        "and the wrong fringe spacing, which is exactly the error that survives "
        "a plausibility check."
    ),
    setup=DYNAMICAL_SETUP,
    code=(
        "beams = beam_set_from_indices(nickel, zone, [[2, 2, 0]])\n"
        "g_zone = beams.g_zone[1]\n"
        "in_plane = float(np.linalg.norm(g_zone[:2]))\n"
        "zero_tilt = g_zone[2] - 0.5 * wavelength * beams.g_magnitude_inv_angstrom[1] ** 2\n"
        "targets = np.array([-0.01, -0.003, 0.0, 0.003, 0.01])\n"
        "tilts = ((zero_tilt - targets) / in_plane)[:, None] * (g_zone[:2] / in_plane)[None, :]\n"
        "solution = solve_bloch_waves(beams, tilts, thickness_angstrom=800.0)\n"
        "closed_form = two_beam_rocking_curve(\n"
        "    targets,\n"
        "    thickness_angstrom=800.0,\n"
        "    extinction_distance_angstrom=float(\n"
        "        extinction_distance_angstrom(nickel, [[2, 2, 0]])[0]\n"
        "    ),\n"
        ")\n"
        "result = float(np.max(np.abs(solution.intensity_of([2, 2, 0]) - closed_form)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "An analytic identity, not a measurement: the two-beam structure matrix "
        "is s I + B with B traceless, and the exponential of a traceless 2x2 "
        "matrix is cos(pi s_eff t) I + i sin(pi s_eff t) B / s_eff, which gives "
        "the Howie-Whelan expression exactly. The deviation must therefore be "
        "zero to floating-point rounding."
    ),
    citation=(
        "Howie and Whelan, Proceedings of the Royal Society A 263 (1961) "
        "217-237; Williams and Carter, Transmission Electron Microscopy, "
        "2nd ed. (Springer, 2009), Chapter 23."
    ),
    symbols=(_EXTINCTION, _DEVIATION, _COEFFICIENT, _THICKNESS),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.2e}",
)


UNITARITY = WorkedExample(
    id="diffraction-dynamical-intensity-is-conserved-without-absorption",
    title="Without absorption the coupled beams sum to exactly one",
    domain="diffraction",
    scenario=(
        "The elastic structure matrix is Hermitian, so the propagator "
        "exp(i pi A t) is unitary and the beam intensities sum to one at every "
        "thickness and every incident direction. This is the only exact global "
        "check available on a many-beam calculation, and it is the one that "
        "catches the classic implementation error: obtaining the Bloch-wave "
        "excitation amplitudes by projection rather than by solving. The "
        "eigenvectors of a complex matrix are not orthogonal, so a projection "
        "gives rocking curves of the right shape with the wrong contrast - "
        "wrong in a way that only this identity reveals."
    ),
    setup=DYNAMICAL_SETUP,
    code=(
        "beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)\n"
        "tilts = np.array([[0.0, 0.0], [3e-3, -2e-3], [-4e-3, 1e-3], [2e-3, 5e-3]])\n"
        "deviations = [\n"
        "    float(np.max(np.abs(\n"
        "        solve_bloch_waves(beams, tilts, thickness_angstrom=t).total_intensity - 1.0\n"
        "    )))\n"
        "    for t in (100.0, 700.0, 2500.0)\n"
        "]\n"
        "result = float(max(deviations))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "Unitarity of the propagator of a Hermitian generator: sum_g |psi_g|^2 "
        "is conserved exactly. The expected value is zero by theorem, with the "
        "tolerance set by floating-point accumulation over the beam set rather "
        "than by any physical uncertainty."
    ),
    citation=(
        "Hirsch, Howie, Nicholson, Pashley and Whelan, Electron Microscopy of "
        "Thin Crystals, 2nd ed. (Krieger, 1977), Chapter 10."
    ),
    symbols=(_COEFFICIENT, _DEVIATION, _THICKNESS),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.2e}",
)


HOLZ_DEGENERACY = WorkedExample(
    id="diffraction-holz-strain-and-wavelength-are-exactly-degenerate",
    title="A HOLZ line cannot separate lattice strain from accelerating voltage",
    domain="diffraction",
    scenario=(
        "HOLZ line positions are the sharpest lattice-parameter measurement a "
        "convergent-beam pattern offers, and this is the reason the measurement "
        "begins with a calibration rather than a specimen. Scaling the lattice "
        "by 1 + eps shrinks every g by the same factor, and the line offset "
        "d_g = (g_z - lambda |g|^2 / 2) / |g_perp| then depends on eps and on "
        "lambda through the same term with opposite signs. A fractional change "
        "in lattice parameter and a fractional change in wavelength therefore "
        "cancel exactly, at every reflection simultaneously - so a lattice "
        "parameter quoted from an uncalibrated microscope is a measurement of "
        "its high-tension supply."
    ),
    setup=DYNAMICAL_SETUP,
    code=(
        "lines = holz_line_pattern(\n"
        "    nickel,\n"
        "    zone,\n"
        "    convergence_semi_angle_mrad=8.0,\n"
        "    max_index=24,\n"
        "    g_max_inv_angstrom=6.0,\n"
        ")\n"
        "strain = 1e-3\n"
        "result = float(max(\n"
        "    abs(\n"
        "        line.offset_at(\n"
        "            lattice_strain=strain,\n"
        "            wavelength_angstrom=wavelength * (1.0 + strain),\n"
        "        )\n"
        "        - line.offset_rad\n"
        "    )\n"
        "    for line in lines.lines\n"
        "))"
    ),
    expected=0.0,
    unit="radian",
    tolerance=1e-15,
    reference=(
        "An exact cancellation in the closed form: substituting "
        "lambda -> lambda (1 + eps) into d_g(eps, lambda) recovers d_g(0, "
        "lambda) identically, for every reflection. The expected value is zero "
        "by algebra, and the tolerance is floating-point rounding."
    ),
    citation=(
        "Jones, Rackham and Steeds, Proceedings of the Royal Society A 354 "
        "(1977) 197-222, for HOLZ line lattice-parameter determination; "
        "Williams and Carter, Transmission Electron Microscopy, 2nd ed. "
        "(2009), Chapter 21."
    ),
    symbols=(_WAVELENGTH, _DEVIATION),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.2e}",
)


THIRTY_ONE_DIFFRACTION_GROUPS = WorkedExample(
    id="diffraction-groups-construction-yields-buxtons-thirty-one",
    title="The diffraction-group construction yields Buxton's 31 groups",
    domain="diffraction",
    scenario=(
        "The 31 diffraction groups are usually quoted as a table. PyTex derives "
        "them instead: each crystal-point-group operator is classified by its "
        "action on the beam direction and contributes its transverse "
        "restriction, tagged with the reciprocity flag when it reverses the "
        "beam. That map lands in a subgroup of (plane point group) x Z2, and "
        "scanning all 32 crystallographic point groups over their "
        "characteristic beam directions must realize exactly the 31 subgroups "
        "Buxton, Eades, Steeds and Rackham enumerated. Reaching 30 or 32 would "
        "mean the construction or the stored operators are wrong - a check no "
        "transcribed table can perform on itself."
    ),
    setup=GROUPS_SETUP,
    code="result = len(diffraction_group_symbols())",
    expected=31,
    unit="",
    tolerance=0.0,
    reference=(
        "The published count of diffraction groups: 10 with no reciprocity "
        "element, 10 direct products with Z2 (suffix 1_R), and 11 graphs of a "
        "surjection onto Z2. An exact integer, so the tolerance is zero."
    ),
    citation=(
        "Buxton, Eades, Steeds and Rackham, Philosophical Transactions of the "
        "Royal Society A 281 (1976) 171-194."
    ),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.0f}",
)


CENTROSYMMETRY_SPLIT = WorkedExample(
    id="diffraction-groups-friedel-observation-splits-the-point-groups",
    title="The plus-minus-g observation splits the 32 point groups into 21 and 11",
    domain="diffraction",
    scenario=(
        "This is the arithmetic of the whole technique. Friedel's law makes "
        "kinematic diffraction blind to a centre of symmetry, so a "
        "selected-area pattern determines only the Laue class: 11 possibilities "
        "where there are 32 point groups. The diffraction-group element 2_R "
        "requires an operator acting as -1 on the beam direction and as -1 on "
        "the transverse plane, which is the inversion and nothing else - so "
        "observing whether the +g and -g discs are related by a two-fold "
        "recovers exactly the distinction Friedel's law destroyed, partitioning "
        "the 32 point groups into the 21 acentric and the 11 centric ones."
    ),
    setup=GROUPS_SETUP,
    code=(
        "acentric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=False))\n"
        "centric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=True))\n"
        "result = [len(acentric.point_groups), len(centric.point_groups)]"
    ),
    expected=[21, 11],
    unit="",
    tolerance=0.0,
    reference=(
        "Of the 32 crystallographic point groups, 11 contain the inversion - "
        "the Laue classes - and 21 do not. Exact integers from the "
        "International Tables, so the tolerance is zero."
    ),
    citation=(
        "International Tables for Crystallography, Volume A, Chapter 10, for "
        "the 32 point groups and the 11 Laue classes; Buxton, Eades, Steeds "
        "and Rackham, Philosophical Transactions of the Royal Society A 281 "
        "(1976) 171-194, for the 2_R correspondence."
    ),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.0f}",
)


ZINCBLENDE_DIFFRACTION_GROUP = WorkedExample(
    id="diffraction-groups-zincblende-down-001-loses-the-two-fold",
    title="Zincblende down [001] gives 4_Rmm_R, and no centre of symmetry",
    domain="diffraction",
    scenario=(
        "The textbook case. Down a four-fold zone the centrosymmetric cubic "
        "group m-3m gives diffraction group 4mm1_R, four-fold in both the "
        "bright-field disc and the whole pattern; zincblende -43m gives "
        "4_Rmm_R, four-fold in the disc but only two-fold in the pattern. That "
        "difference in whole-pattern symmetry is what a CBED exposure of "
        "gallium arsenide reads, and it is the observation that determines the "
        "absence of a centre of symmetry - which no kinematic pattern can do. "
        "The check below counts the whole-pattern operations: 4 for -43m "
        "against 8 for m-3m, that is 2mm against 4mm."
    ),
    setup=GROUPS_SETUP,
    code=(
        "polar = diffraction_group_for('-43m', [0, 0, 1])\n"
        "centric = diffraction_group_for('m-3m', [0, 0, 1])\n"
        "result = [\n"
        "    float(polar.symbol == '4_Rmm_R'),\n"
        "    float(polar.bright_field_symbol == '4mm'),\n"
        "    float(polar.whole_pattern_symbol == '2mm'),\n"
        "    float(polar.has_friedel_symmetry),\n"
        "    float(centric.symbol == '4mm1_R'),\n"
        "    float(centric.whole_pattern_symbol == '4mm'),\n"
        "    float(centric.has_friedel_symmetry),\n"
        "]"
    ),
    expected=[1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0],
    unit="",
    tolerance=0.0,
    reference=(
        "The published diffraction-group assignments for the cubic acentric "
        "and centric groups viewed down a four-fold axis: -43m gives 4_Rmm_R "
        "with bright-field 4mm over whole-pattern 2mm and no 2_R, while m-3m "
        "gives 4mm1_R with 4mm in both and 2_R present. Each entry is a "
        "boolean agreement, so the tolerance is zero."
    ),
    citation=(
        "Buxton, Eades, Steeds and Rackham, Philosophical Transactions of the "
        "Royal Society A 281 (1976) 171-194, Tables 2 and 3; Williams and "
        "Carter, Transmission Electron Microscopy, 2nd ed. (2009), Chapter 21."
    ),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.0f}",
)


GROUP = ExampleGroup(
    slug="dynamical-cbed-and-symmetry",
    title="Dynamical CBED and symmetry determination",
    summary=(
        "The exact limits that calibrate a many-beam calculation, the HOLZ "
        "degeneracy that makes voltage calibration mandatory, and the "
        "diffraction-group construction that determines a point group "
        "including its centre of symmetry."
    ),
    examples=(
        TWO_BEAM_LIMIT,
        UNITARITY,
        HOLZ_DEGENERACY,
        THIRTY_ONE_DIFFRACTION_GROUPS,
        CENTROSYMMETRY_SPLIT,
        ZINCBLENDE_DIFFRACTION_GROUP,
    ),
)
