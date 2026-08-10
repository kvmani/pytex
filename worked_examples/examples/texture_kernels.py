"""Worked examples: SO(3) texture kernels and the m.r.d. scale.

The kernel surface is verifiable through two analytic identities: a
normalized SO(3) kernel has zeroth Chebyshev (character) coefficient exactly
one, and its halfwidth is defined by ``psi(halfwidth) = psi(0) / 2``. Both
are computed live for the spectral Gaussian (Gauss-Weierstrass) kernel.

A third identity fixes the scale on which pole densities are reported: a
uniform orientation distribution sends poles uniformly over the sphere, so
its pole figure is flat at exactly one multiple of a random distribution, in
every direction and for every plane. Preserving that identity is what
``random_pole_density`` is for, and checking it is what catches a kernel
response being mistaken for a density.

See :doc:`../../concepts/orientation_texture` for the surrounding ODF
doctrine.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

KERNEL_SETUP = """
import numpy as np
from pytex import GaussianSO3Kernel
"""

_PSI = SymbolUse(
    r"\psi(\omega)",
    "SO(3) radial kernel density as a function of the rotation angle.",
)

_TEXTURE_CONCEPT = SeeAlso("Orientations and texture", "../../concepts/orientation_texture")
_API = SeeAlso("Texture API", "../../api/index")

MRD_SETUP = """
import numpy as np
from pytex import (
    ODF,
    CrystalPlane,
    FrameDomain,
    KernelSpec,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    random_pole_density,
)

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

_POLE_DENSITY = SymbolUse(
    r"P_{\mathbf{h}}(\mathbf{y})",
    "Pole density of plane normal h along specimen direction y, in multiples "
    "of a random distribution.",
)


GAUSSIAN_KERNEL_IDENTITIES = WorkedExample(
    id="texture-gaussian-kernel-normalization-and-halfwidth",
    title="The Gaussian SO(3) kernel is normalized and honors its halfwidth",
    domain="texture",
    scenario=(
        "Construct a Gaussian (Gauss-Weierstrass) kernel with a 10 degree "
        "halfwidth and verify the two defining identities: the zeroth "
        "Chebyshev coefficient equals one (the kernel integrates to one over "
        "SO(3) with the normalized Haar measure), and the density at the "
        "halfwidth equals half the peak density."
    ),
    setup=KERNEL_SETUP,
    code=(
        "kernel = GaussianSO3Kernel(10.0)\n"
        "a0 = float(kernel.chebyshev_coefficients(0)[0])\n"
        "ratio = float(\n"
        "    kernel.evaluate(np.array([np.deg2rad(10.0)]))[0]\n"
        "    / kernel.evaluate(np.array([0.0]))[0]\n"
        ")\n"
        "result = np.array([a0, ratio])"
    ),
    expected=[1.0, 0.5],
    unit="",
    tolerance=1e-6,
    reference=(
        "A_0 = 1 is the SO(3) normalization identity for character "
        "expansions, and psi(halfwidth) = psi(0)/2 is the definition of the "
        "kernel halfwidth (both analytic identities)."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982), harmonic "
        "expansion of ODF kernels."
    ),
    symbols=(_PSI,),
    see_also=(_TEXTURE_CONCEPT, _API),
    result_format="{:.6f}",
)


UNIFORM_ODF_POLE_DENSITY = WorkedExample(
    id="texture-uniform-odf-pole-density-is-one-mrd",
    title="A uniform ODF gives a pole figure flat at one m.r.d.",
    domain="texture",
    scenario=(
        "Build an ODF from an equispaced grid over the cubic fundamental "
        "region with equal weights - a texture-free aggregate - and evaluate "
        "its {111} pole density along three unrelated specimen directions. A "
        "discrete ODF's evaluate_pole_density returns a kernel-weighted "
        "response, whose peak is one rather than whose integral is one, so "
        "converting it to multiples of a random distribution means dividing by "
        "random_pole_density: the response a random texture produces. Skipping "
        "that division is a scale error of about two orders of magnitude, not "
        "a small one."
    ),
    setup=MRD_SETUP,
    code=(
        "dictionary = OrientationSet.from_equispaced_so3_grid(\n"
        "    10.0,\n"
        "    crystal_frame=crystal,\n"
        "    specimen_frame=specimen,\n"
        "    symmetry=symmetry,\n"
        "    phase=phase,\n"
        ")\n"
        "kernel = KernelSpec(name='de_la_vallee_poussin', halfwidth_deg=15.0)\n"
        "odf = ODF.from_orientations(dictionary, kernel=kernel)\n"
        "directions = np.array(\n"
        "    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0] / np.sqrt(3.0)]\n"
        ")\n"
        "response = np.asarray(odf.evaluate_pole_density(pole, directions))\n"
        "result = response / random_pole_density(kernel)"
    ),
    expected=[1.0, 1.0, 1.0],
    unit="m.r.d.",
    tolerance=1e-3,
    reference=(
        "Analytic identity: a uniform orientation distribution maps poles "
        "uniformly onto the sphere, so every pole density equals one multiple "
        "of a random distribution by the definition of the m.r.d. scale. The "
        "residual deviation is the finite orientation grid, not the scale."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982), Sec. 4 - "
        "normalization of the ODF and of pole figures to multiples of a random "
        "distribution."
    ),
    symbols=(_POLE_DENSITY,),
    see_also=(_TEXTURE_CONCEPT, _API),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="texture",
    title="Texture kernels",
    summary=(
        "Analytic identities of the SO(3) kernel surface - normalization "
        "(A_0 = 1) and the halfwidth definition - together with the m.r.d. "
        "scale on which pole densities are reported, all computed live."
    ),
    examples=(GAUSSIAN_KERNEL_IDENTITIES, UNIFORM_ODF_POLE_DENSITY),
)
