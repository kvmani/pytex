"""Worked examples: what a diffraction pole figure cannot determine.

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


GROUP = ExampleGroup(
    slug="ghost-problem",
    title="The ghost problem",
    summary=(
        "What diffraction pole figures cannot determine: an asymmetric texture "
        "still gives a pole set closed under negation, and excluding the odd "
        "harmonic degrees that centrosymmetry annihilates discards nearly half "
        "the basis."
    ),
    examples=(
        POLE_FIGURE_IS_CENTROSYMMETRIC,
        ODD_DEGREES_ARE_HALF_THE_BASIS,
    ),
)
