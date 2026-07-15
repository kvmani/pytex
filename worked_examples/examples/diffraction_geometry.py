"""Worked examples: diffraction geometry from crystallographic spacings.

These examples connect the PyTex crystallographic core to a diffraction
observable. The interplanar spacing d_hkl is computed from the metric tensor by
PyTex, and Bragg's law lambda = 2 d sin(theta) then fixes the powder scattering
angle 2*theta. The reference is the well-known Ni(111) reflection position for
Cu K-alpha1 radiation, near 44.5 degrees.

See :doc:`../../workflows/xrd_generation` and the diffraction theory notes.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

NICKEL_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom
"""

_DSPACING = SymbolUse(r"d_{hkl}", "Interplanar spacing of the (hkl) family.")
_THETA = SymbolUse(r"\theta", "Bragg half-angle.")
_TWO_THETA = SymbolUse(r"2\theta", "Powder-diffraction scattering angle.")
_LAMBDA = SymbolUse(r"\lambda", "Radiation wavelength.")

_XRD_WORKFLOW = SeeAlso("Powder XRD generation", "../../workflows/xrd_generation")
_DIFF_CONCEPT = SeeAlso("Diffraction foundation", "../../concepts/diffraction_foundation")


NI_111_TWO_THETA = WorkedExample(
    id="diffraction-ni-111-two-theta",
    title="Ni(111) powder reflection angle for Cu K-alpha1",
    domain="diffraction",
    scenario=(
        "You are calibrating or interpreting a powder pattern and need to predict where the Ni(111) "
        "peak should appear with a copper source. PyTex supplies the interplanar spacing from the "
        "lattice metric; Bragg's law then gives the scattering angle. The result should land on the "
        "textbook Ni(111) position near 44.5 degrees for Cu K-alpha1."
    ),
    setup=NICKEL_SETUP,
    code=(
        "d_111 = MillerPlane.from_hkl([1, 1, 1], phase=nickel).d_spacing_angstrom\n"
        "theta = np.arcsin(cu_ka1 / (2.0 * d_111))\n"
        "result = float(np.degrees(2.0 * theta))"
    ),
    expected=44.496,
    unit="deg",
    tolerance=5e-3,
    reference=(
        "d_111 = 3.52387 / sqrt(3) = 2.03451 angstrom; with lambda = 1.5406 angstrom, "
        "2*theta = 2*arcsin(lambda / (2 d)) = 44.50 degrees, matching standard Ni powder data."
    ),
    citation="ICDD PDF 04-0850 (nickel); Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed.",
    symbols=(_DSPACING, _LAMBDA, _THETA, _TWO_THETA),
    see_also=(_XRD_WORKFLOW, _DIFF_CONCEPT),
    result_format="{:.3f}",
)


GROUP = ExampleGroup(
    slug="diffraction",
    title="Diffraction geometry",
    summary=(
        "Powder scattering angles derived from PyTex interplanar spacings via Bragg's law, checked "
        "against a standard reference reflection position."
    ),
    examples=(NI_111_TWO_THETA,),
)

__all__ = ["GROUP"]
