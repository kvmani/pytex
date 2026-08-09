"""Worked examples: the CBED absolute scale, and thickness from fringes.

Convergent-beam diffraction is where the *absolute* scale of an electron
structure factor first matters. Disc geometry and rocking curves look right for
any extinction distance; only the fringe spacing knows the difference. So the
first example checks the extinction distance against a published table, and the
second checks that the standard two-beam thickness analysis inverts its own
defining relation exactly.

See :doc:`../../concepts/diffraction_foundation` for the surrounding
diffraction doctrine.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

CBED_SETUP = """
import numpy as np
from pytex import (
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    extinction_distance_angstrom,
    thickness_from_fringe_minima,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
aluminium_lattice = Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal)
aluminium = Phase(
    name="aluminium-fcc",
    lattice=aluminium_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=aluminium_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Al{index}",
                species="Al",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
"""

_EXTINCTION = SymbolUse(
    r"\xi_{g}",
    "Two-beam extinction distance of reflection g; the depth period of the "
    "intensity exchange between the transmitted and diffracted beams.",
)
_THICKNESS = SymbolUse(
    r"t",
    "Foil thickness along the beam.",
)
_DEVIATION = SymbolUse(
    r"s_{g}",
    "Excitation error: deviation of reflection g from the exact Bragg condition.",
)

_DIFFRACTION_CONCEPT = SeeAlso("Diffraction foundation", "../../concepts/diffraction_foundation")
_API = SeeAlso("Diffraction API", "../../api/index")


EXTINCTION_DISTANCES = WorkedExample(
    id="diffraction-cbed-aluminium-extinction-distances-at-100kv",
    title="Extinction distances of aluminium at 100 kV match the published table",
    domain="diffraction",
    scenario=(
        "Every dynamical quantity in a CBED pattern is measured in units of "
        "the extinction distance, so its absolute scale has to be right - and "
        "a wrong scale is invisible in the geometry, showing only as fringes "
        "at the wrong spacing. Two things set that scale: the Mott-Bethe "
        "conversion of the X-ray form factor into an electron scattering "
        "factor in angstrom, and the relativistic factor gamma = 1 + E/m0c^2, "
        "which is 1.20 at 100 kV. Aluminium is the calibration case because "
        "the fitted scattering-factor parametrization is most accurate for "
        "light elements; for heavy elements the same calculation is only "
        "indicative, which is why CBED practice measures the extinction "
        "distance rather than tabulating it."
    ),
    setup=CBED_SETUP,
    code=(
        "result = extinction_distance_angstrom(\n"
        "    aluminium, [(1, 1, 1), (2, 0, 0), (2, 2, 0)], beam_energy_kev=100.0\n"
        ")"
    ),
    expected=[556.0, 673.0, 1057.0],
    unit="angstrom",
    tolerance=10.0,
    reference=(
        "Published two-beam extinction distances for aluminium at 100 kV, "
        "556, 673 and 1057 angstrom for {111}, {200} and {220}. The tolerance "
        "of 10 angstrom is about 1.5 percent, the accuracy the fitted "
        "scattering-factor parametrization supports for a light element."
    ),
    citation=(
        "Williams and Carter, Transmission Electron Microscopy, 2nd ed. "
        "(Springer, 2009), Table 23.1."
    ),
    symbols=(_EXTINCTION,),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.1f}",
)


THICKNESS_FROM_FRINGES = WorkedExample(
    id="diffraction-cbed-two-beam-thickness-inverts-the-fringe-relation",
    title="The Kelly plot recovers both the thickness and the extinction distance from fringe positions",
    domain="diffraction",
    scenario=(
        "Reading the dark fringes off a single CBED disc gives the local foil "
        "thickness - and, from the same straight-line fit, the extinction "
        "distance, so the thickness does not inherit the error of a tabulated "
        "constant. Here the fringe positions are generated from the two-beam "
        "relation for a chosen thickness of 2000 angstrom and extinction "
        "distance of 500 angstrom, and the fit is asked to recover both. "
        "Because the input is the closed-form relation rather than a "
        "simulation, this tests the inversion itself and nothing else."
    ),
    setup=CBED_SETUP,
    code=(
        "thickness, extinction = 2000.0, 500.0\n"
        "orders = np.arange(5, 11, dtype=float)\n"
        "minima = np.sqrt((orders / thickness) ** 2 - extinction**-2)\n"
        "report = thickness_from_fringe_minima(minima, first_order=5)\n"
        "result = np.array([\n"
        "    report.thickness_angstrom,\n"
        "    report.extinction_distance_angstrom,\n"
        "])"
    ),
    expected=[2000.0, 500.0],
    unit="angstrom",
    tolerance=1e-6,
    reference=(
        "An exact inversion of the two-beam minimum condition "
        "t s_eff,n = n with s_eff^2 = s^2 + xi^-2, which rearranges to "
        "(s_n/n)^2 = 1/t^2 - (1/xi^2)(1/n^2). The generated minima lie on that "
        "line by construction, so a least-squares fit returns the intercept "
        "1/t^2 and the slope -1/xi^2 to machine precision."
    ),
    citation=(
        "Kelly, Jostsons, Blake and Napier, Physica Status Solidi (a) 31 "
        "(1975) 771-780, for the linearization; Williams and Carter, "
        "Transmission Electron Microscopy, 2nd ed. (2009), Chapter 23."
    ),
    symbols=(_THICKNESS, _EXTINCTION, _DEVIATION),
    see_also=(_DIFFRACTION_CONCEPT, _API),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="convergent-beam-diffraction",
    title="Convergent-beam diffraction",
    summary=(
        "The absolute scale of the two-beam extinction distance, checked "
        "against a published table, and the fringe analysis that measures a "
        "foil thickness without needing one."
    ),
    examples=(EXTINCTION_DISTANCES, THICKNESS_FROM_FRINGES),
)
