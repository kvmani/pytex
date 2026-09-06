"""Worked examples: what a diffraction pole figure cannot determine, and the correction.

A lattice plane has no sense and Friedel's law makes +h and -h scatter
identically, so a measured pole figure is centrosymmetric whatever the ODF
is. That kills every odd-degree harmonic of the ODF exactly - the ghost
problem - and it removes close to half the basis.

Both facts are checked here rather than asserted: a deliberately one-sided
orientation population still yields a pole set closed under negation, and
the harmonic term count shows the fraction discarded approaching 1/2.

See :doc:`../../theory/ghost_problem_and_odd_harmonics`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

GHOST_SETUP = """
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.texture.models import PoleFigure

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
"""

_PF = SymbolUse(
    r"P_{\mathbf{h}}(\mathbf{y})",
    "Pole density of plane family h along specimen direction y.",
)
_ELL = SymbolUse(r"\ell", "Degree of a generalized spherical harmonic term.")

_THEORY = SeeAlso(
    "The ghost problem", "../../theory/ghost_problem_and_odd_harmonics"
)
_HARMONIC = SeeAlso(
    "Harmonic ODF reconstruction", "../../theory/harmonic_odf_reconstruction"
)


POLE_FIGURE_IS_CENTROSYMMETRIC = WorkedExample(
    id="ghost-pole-figure-is-centrosymmetric",
    title="An asymmetric texture still gives a centrosymmetric pole figure",
    domain="texture",
    scenario=(
        "Build a deliberately one-sided orientation population - 200 "
        "orientations whose first Euler angle is confined to 0-40 degrees, "
        "with nothing symmetric about it - and generate its {111} pole figure. "
        "Every pole in the result has its antipode also present. The "
        "centrosymmetry is a property of the measurement, not of the specimen: "
        "a lattice plane has no sense and Friedel's law makes +h and -h "
        "scatter identically. This is the root cause of the ghost problem, and "
        "the reason no amount of pole-figure data can recover the odd part of "
        "an ODF. The example returns the fraction of the first 300 poles whose "
        "antipode is in the set, which must be 1."
    ),
    setup=GHOST_SETUP,
    code=(
        "angles = np.linspace(0.0, 40.0, 200)\n"
        "orientations = OrientationSet.from_euler_angles(\n"
        "    np.column_stack(\n"
        "        [angles, 25.0 * np.ones_like(angles), np.zeros_like(angles)]\n"
        "    ),\n"
        "    crystal_frame=crystal,\n"
        "    specimen_frame=specimen,\n"
        "    symmetry=symmetry,\n"
        ")\n"
        "figure = PoleFigure.from_orientations(orientations, pole)\n"
        "directions = figure.sample_directions\n"
        "poles = np.asarray(\n"
        "    getattr(directions, 'values', directions)\n"
        ")\n"
        "paired = sum(\n"
        "    1\n"
        "    for row in poles[:300]\n"
        "    if np.min(np.linalg.norm(poles + row, axis=1)) < 1e-9\n"
        ")\n"
        "result = paired / 300.0"
    ),
    expected=1.0,
    unit="",
    tolerance=0.0,
    reference=(
        "Analytic: a plane normal enters a pole figure as an axis and Friedel's "
        "law equates +h with -h, so the pole set is closed under negation for "
        "any orientation distribution whatsoever."
    ),
    citation=(
        "Matthies, On the reproducibility of the orientation distribution "
        "function of texture samples from pole figures (ghost phenomena), "
        "Phys. Status Solidi B 92 (1979) K135-K138."
    ),
    symbols=(_PF,),
    see_also=(_THEORY, _HARMONIC),
    result_format="{:.6f}",
)


ODD_DEGREES_ARE_HALF_THE_BASIS = WorkedExample(
    id="ghost-odd-degrees-are-half-the-harmonic-basis",
    title="Excluding odd degrees discards nearly half the harmonic basis",
    domain="texture",
    scenario=(
        "Count the generalized spherical harmonic terms retained and discarded "
        "when odd degrees are excluded, at bandlimit 22. Because pole figures "
        "annihilate every odd degree, those coefficients are not poorly "
        "determined but wholly undetermined, and the count says how much of "
        "the ODF a diffraction measurement is silent about: 15147 of 32407 "
        "terms, or 46.7 percent, tending to one half as the bandlimit grows. "
        "The example returns the discarded fraction."
    ),
    setup="""
from pytex.texture.harmonics import _enumerate_terms
""",
    code=(
        "even = len(\n"
        "    _enumerate_terms(degree_bandlimit=22, even_degrees_only=True)\n"
        ")\n"
        "every = len(\n"
        "    _enumerate_terms(degree_bandlimit=22, even_degrees_only=False)\n"
        ")\n"
        "result = (every - even) / every"
    ),
    expected=0.46739,
    unit="",
    tolerance=1e-4,
    reference=(
        "Exact term count: sum over degrees 0..22 of (2l+1)^2 is 32407, of "
        "which the even degrees contribute 17260, leaving 15147 discarded. The "
        "fraction tends to 1/2 as the bandlimit grows."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science: Mathematical Methods "
        "(Butterworths 1969) - the generalized spherical harmonic expansion."
    ),
    symbols=(_ELL,),
    see_also=(_THEORY, _HARMONIC),
    result_format="{:.5f}",
)


INVARIANT_SETUP = """
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
cubic = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)


def invariant_count(symmetry, degree):
    \"\"\"Dimension of the degree-l invariant subspace of a rotation group.

    Character theory: the number of invariants is the group average of the
    character of the degree-l representation of SO(3),
    chi_l(theta) = sin((l + 1/2) theta) / sin(theta / 2).
    \"\"\"

    operators = np.asarray(symmetry.operators)
    cosines = np.clip((np.trace(operators, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.arccos(cosines)
    half = np.sin(theta / 2.0)
    identity = np.abs(half) < 1e-12
    chi = np.where(
        identity,
        2.0 * degree + 1.0,
        np.sin((degree + 0.5) * theta) / np.where(identity, 1.0, half),
    )
    return float(np.mean(chi))
"""

CUBIC_FIRST_ODD_INVARIANT = WorkedExample(
    id="ghost-cubic-first-odd-invariant-is-degree-nine",
    title="A cubic material has no odd harmonic to correct below degree 9",
    domain="texture",
    scenario=(
        "Ghost correction supplies the odd part of an ODF that a pole figure "
        "cannot measure. Whether there is an odd part to supply at all is a "
        "question about the symmetry, not about the data: the crystal rotation "
        "group admits odd-degree terms only where it has an odd-degree "
        "invariant. Counting those invariants by character theory - the group "
        "average of the SO(3) character - gives the classical answer for the "
        "cubic rotation group 432: nothing at degrees 1, 3, 5 or 7, and the "
        "first odd invariant at degree 9. A cubic ODF expanded to degree 6 or 8 "
        "therefore has no ghost part, and PyTex's correction reports that "
        "rather than a correction of size zero. Lower symmetries admit odd "
        "terms much earlier - degree 7 for hexagonal 622, degree 3 for "
        "orthorhombic 222."
    ),
    setup=INVARIANT_SETUP,
    code=(
        "odd_degrees = range(1, 16, 2)\n"
        "result = min(\n"
        "    degree for degree in odd_degrees if invariant_count(cubic, degree) > 0.5\n"
        ")"
    ),
    expected=9.0,
    unit="",
    tolerance=0.0,
    reference=(
        "Standard result of cubic harmonic analysis: the cubic rotation group "
        "has invariants at degrees 0, 4, 6, 8, 9, 10, ... and the lowest "
        "odd-degree cubic harmonic is degree 9. Independently derivable from "
        "the Molien series of the octahedral group."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science: Mathematical Methods "
        "(Butterworths 1969), tables of cubic symmetric generalized spherical "
        "harmonics."
    ),
    symbols=(_ELL,),
    see_also=(_THEORY, _HARMONIC),
    result_format="{:.0f}",
)


CORRECTION_SETUP = """
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    HarmonicODF,
    KernelSpec,
    Lattice,
    MillerIndex,
    ODF,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    SymmetrySpec,
    random_pole_density,
)
from pytex.diffraction.stereonets import spherical_angles_to_directions

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
# Orthorhombic, because 222 admits odd-degree terms from degree 3; a cubic
# material has no ghost part below degree 9 and nothing to demonstrate here.
symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal)
phase = Phase(
    name="orthorhombic-demo",
    lattice=Lattice(3.0, 4.0, 5.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
# A single broad component: broad enough that a degree-4 expansion represents
# it, so the demonstration is of the ghost problem and not of truncation.
truth = ODF(
    orientations=OrientationSet.from_euler_angles(
        np.array([[35.0, 45.0, 20.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    ),
    weights=np.array([1.0]),
    kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=50.0),
)
polar, azimuth = np.meshgrid(
    np.arange(0.0, 91.0, 15.0), np.arange(0.0, 360.0, 15.0), indexing="ij"
)
directions = spherical_angles_to_directions(polar, azimuth).reshape(-1, 3)
scale = random_pole_density(truth.kernel, antipodal=True)


def measured(indices):
    pole = CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase)
    return PoleFigure(
        pole=pole,
        sample_directions=directions,
        intensities=truth.evaluate_pole_density(pole, directions, antipodal=True) / scale,
        specimen_frame=specimen,
        antipodal=True,
        includes_symmetry_family=True,
        sampling="sampled_density",
    )


pole_figures = [measured(indices) for indices in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0])]
report = HarmonicODF.invert_pole_figures(
    pole_figures,
    degree_bandlimit=4,
    regularization=1e-6,
    pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=10.0),
    phi1_step_deg=15.0,
    big_phi_step_deg=15.0,
    phi2_step_deg=15.0,
    ghost_correction=True,
)
correction = report.ghost_correction
"""

CORRECTION_RESTORES_POSITIVITY = WorkedExample(
    id="ghost-correction-restores-a-non-negative-density",
    title="Ghost correction removes the negative density the even part leaves",
    domain="texture",
    scenario=(
        "The even-degree half of an ODF is all a pole figure determines, and "
        "on its own it is not a density: here it falls to about -0.46 multiples "
        "of random over part of orientation space, which is the classical ghost "
        "artefact. Correction holds that even part fixed - it is what the data "
        "say - and adds the smallest odd part that makes the whole "
        "non-negative. Because the added part is the smallest one that works, "
        "the constraint ends up exactly tight: the minimum density of the "
        "corrected distribution sits on zero rather than comfortably above it. "
        "The maximum rises at the same time, from 3.75 to 4.22 m.r.d., which is "
        "the other half of the ghost signature: the even-only solution pays for "
        "its false lobes by depressing the true peak."
    ),
    setup=CORRECTION_SETUP,
    code=(
        "result = float(np.min(correction.odf.quadrature_densities))"
    ),
    expected=0.0,
    unit="m.r.d.",
    tolerance=1e-3,
    reference=(
        "Analytic: an orientation distribution is a probability density and "
        "cannot be negative, and the correction minimizes the norm of the odd "
        "part subject to that constraint. A minimum-norm feasible point lies on "
        "the boundary of the feasible set whenever the unconstrained point is "
        "infeasible, so the minimum density is zero to within the quadrature "
        "resolution."
    ),
    citation=(
        "Dahms and Bunge, The iterative series-expansion method for "
        "quantitative texture analysis. I. General outline, J. Appl. Cryst. 22 "
        "(1989) 439-447."
    ),
    symbols=(_ELL,),
    see_also=(_THEORY, _HARMONIC),
    result_format="{:.4f}",
)


CORRECTION_KEEPS_THE_FIT = WorkedExample(
    id="ghost-correction-leaves-the-measured-fit-untouched",
    title="The odd part a correction adds is invisible to the pole figures",
    domain="texture",
    scenario=(
        "A correction that improved the density by moving the fit would be "
        "spending data agreement it has no right to spend, so this is the check "
        "that matters. Under Friedel's law the forward operator annihilates "
        "every odd-degree harmonic exactly, and adding one therefore leaves "
        "every predicted pole density where it was. The example returns the "
        "largest change over all four figures and every measured direction, "
        "against intensities that reach 2.4 m.r.d.; it is nonzero only because "
        "the odd basis is orthonormalized on a discrete quadrature, so its "
        "orthogonality to the even part is exact only in the continuum. "
        "Refining the quadrature drives it to zero."
    ),
    setup=CORRECTION_SETUP,
    code=(
        "result = float(correction.pole_figure_max_change)"
    ),
    expected=0.0,
    unit="m.r.d.",
    tolerance=5e-3,
    reference=(
        "Analytic: a Friedel-symmetric pole figure is the integral of the ODF "
        "over a kernel even under h -> -h, and an odd-degree generalized "
        "spherical harmonic integrates to zero against an even kernel. The "
        "change is therefore exactly zero in the continuum; the tolerance is "
        "the quadrature discretization error at a 15 degree Bunge step."
    ),
    citation=(
        "Matthies, On the reproducibility of the orientation distribution "
        "function of texture samples from pole figures (ghost phenomena), "
        "Phys. Status Solidi B 92 (1979) K135-K138."
    ),
    symbols=(_PF, _ELL),
    see_also=(_THEORY, _HARMONIC),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="ghost-problem",
    title="The ghost problem, and its correction",
    summary=(
        "What diffraction pole figures cannot determine, and what positivity "
        "can recover of it: an asymmetric texture still gives a pole set closed "
        "under negation, excluding the odd harmonic degrees discards nearly "
        "half the basis, a cubic material has no odd term below degree 9, and "
        "the correction removes the negative density without moving the fit."
    ),
    examples=(
        POLE_FIGURE_IS_CENTROSYMMETRIC,
        ODD_DEGREES_ARE_HALF_THE_BASIS,
        CUBIC_FIRST_ODD_INVARIANT,
        CORRECTION_RESTORES_POSITIVITY,
        CORRECTION_KEEPS_THE_FIT,
    ),
)
