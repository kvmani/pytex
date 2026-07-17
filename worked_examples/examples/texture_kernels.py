"""Worked examples: SO(3) texture kernels.

The kernel surface is verifiable through two analytic identities: a
normalized SO(3) kernel has zeroth Chebyshev (character) coefficient exactly
one, and its halfwidth is defined by ``psi(halfwidth) = psi(0) / 2``. Both
are computed live for the spectral Gaussian (Gauss-Weierstrass) kernel.

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


GROUP = ExampleGroup(
    slug="texture",
    title="Texture kernels",
    summary=(
        "Analytic identities of the SO(3) kernel surface: normalization "
        "(A_0 = 1) and the halfwidth definition, computed live."
    ),
    examples=(GAUSSIAN_KERNEL_IDENTITIES,),
)
