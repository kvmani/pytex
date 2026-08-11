"""Worked examples: Schmid factors and Taylor factors against known answers.

Every expected value here is either an exact closed form or the classic
published result:

* the Schmid factor of fcc octahedral slip under [001] tension is exactly
  1/sqrt(6), with eight systems equally stressed and four carrying nothing;
* the full-constraint Taylor factors of the [001] and [111] orientations are
  exactly sqrt(6) and 3 sqrt(6)/2;
* the average Taylor factor of a randomly textured fcc aggregate in tension
  is Taylor's 1938 value of about 3.06.

The single-orientation values test the linear-programming formulation, the
five-component constraint, and the signed-slip doubling together.

See :doc:`../../theory/schmid_and_taylor_plasticity`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

PLASTICITY_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.properties import fcc_octahedral_slip, taylor_factors

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
family = fcc_octahedral_slip(phase)

cube = Orientation(
    rotation=Rotation.identity(),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
"""

_M_SCHMID = SymbolUse(r"m", "Schmid factor, cos(phi) cos(lambda); bounded by 1/2.")
_M_TAYLOR = SymbolUse(
    r"M", "Full-constraint Taylor factor: minimum total slip per unit equivalent strain."
)

_THEORY = SeeAlso(
    "Schmid factors and the Taylor factor",
    "../../theory/schmid_and_taylor_plasticity",
)


CUBE_SCHMID_FACTOR = WorkedExample(
    id="plasticity-fcc-cube-schmid-factor",
    title="Eight fcc systems share a Schmid factor of exactly 1/sqrt(6)",
    domain="properties",
    scenario=(
        "Resolve [001] tension onto the twelve {111}<110> systems of a "
        "cube-oriented fcc grain. The magnitudes take only two values: eight "
        "systems at exactly 1/sqrt(6) = 0.408248 and four at zero. The "
        "eightfold degeneracy is why a cube-oriented grain has no single "
        "preferred slip system, and it is the origin of the Taylor ambiguity - "
        "many different five-system combinations accommodate the same strain "
        "at the same cost. The example returns the largest magnitude."
    ),
    setup=PLASTICITY_SETUP,
    code=(
        "factors = np.asarray(\n"
        "    family.schmid_factors(cube, (0.0, 0.0, 1.0))\n"
        ").ravel()\n"
        "result = float(np.abs(factors).max())"
    ),
    expected=0.4082482904638631,
    unit="",
    tolerance=1e-9,
    reference=(
        "Analytic: for [001] tension the {111}<110> systems give "
        "m = (t.n)(t.d) = (1/sqrt(3))(1/sqrt(2)) = 1/sqrt(6) on the eight "
        "systems that are stressed at all."
    ),
    citation=(
        "Schmid and Boas, Kristallplastizitaet (1935); Kocks, Tome and Wenk, "
        "Texture and Anisotropy (CUP 1998)."
    ),
    symbols=(_M_SCHMID,),
    see_also=(_THEORY,),
    result_format="{:.6f}",
)


CUBE_TAYLOR_FACTOR = WorkedExample(
    id="plasticity-fcc-cube-taylor-factor",
    title="The [001] fcc Taylor factor is exactly sqrt(6)",
    domain="properties",
    scenario=(
        "Solve the full-constraint Taylor problem for a cube-oriented fcc "
        "grain in uniaxial tension along [001]. Minimising total slip subject "
        "to matching the five independent components of the imposed deviatoric "
        "strain gives exactly sqrt(6) = 2.449490. Reproducing a closed form to "
        "six figures exercises the linear-programming formulation, the "
        "five-constraint reduction (the sixth component is implied because "
        "both sides are traceless), and the signed-slip doubling that lets "
        "non-negative variables represent either shear sense."
    ),
    setup=PLASTICITY_SETUP,
    code=(
        "result = float(\n"
        "    taylor_factors(family, cube, tension_axis=(0.0, 0.0, 1.0))\n"
        ")"
    ),
    expected=2.449489742783178,
    unit="",
    tolerance=1e-6,
    reference=(
        "Analytic: the full-constraint Taylor factor of the cube orientation "
        "for {111}<110> slip in uniaxial tension is sqrt(6)."
    ),
    citation=(
        "Taylor, Plastic strain in metals, J. Inst. Metals 62 (1938) 307-324; "
        "Bishop and Hill, Phil. Mag. 42 (1951) 414-427."
    ),
    symbols=(_M_TAYLOR,),
    see_also=(_THEORY,),
    result_format="{:.6f}",
)


RANDOM_TEXTURE_TAYLOR_FACTOR = WorkedExample(
    id="plasticity-random-fcc-taylor-factor",
    title="A random fcc texture has an average Taylor factor near 3.06",
    domain="properties",
    scenario=(
        "Average the full-constraint Taylor factor over 2000 Haar-random "
        "orientations of an fcc aggregate in uniaxial tension, recovering "
        "Taylor's 1938 result of about 3.06. This is the number that converts "
        "a single-crystal critical resolved shear stress into a polycrystal "
        "flow stress. The spread matters as much as the mean: over a random "
        "texture M runs from roughly 2.29 to 3.67, so the hardest orientation "
        "is about 60 percent harder than the softest, which is why a textured "
        "sheet can differ substantially in flow stress from a random one at "
        "the same composition and grain size."
    ),
    setup=PLASTICITY_SETUP
    + """
rng = np.random.default_rng(5)
u1, u2, u3 = rng.random(2000), rng.random(2000), rng.random(2000)
quaternions = np.stack(
    [
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ],
    axis=-1,
)
orientations = OrientationSet.from_quaternions(
    quaternions,
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
""",
    code=(
        "factors = np.asarray(\n"
        "    taylor_factors(\n"
        "        family, orientations, tension_axis=(0.0, 0.0, 1.0)\n"
        "    )\n"
        ")\n"
        "result = float(factors.mean())"
    ),
    expected=3.06,
    unit="",
    tolerance=0.03,
    reference=(
        "Taylor's 1938 value for a randomly oriented fcc aggregate deforming "
        "by {111}<110> slip in uniaxial tension, M ~ 3.06. Tolerance covers "
        "the Monte-Carlo standard error at n = 2000 (about 0.009)."
    ),
    citation=(
        "Taylor, Plastic strain in metals, J. Inst. Metals 62 (1938) 307-324."
    ),
    symbols=(_M_TAYLOR,),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


GROUP = ExampleGroup(
    slug="schmid-and-taylor",
    title="Schmid and Taylor plasticity factors",
    summary=(
        "Slip geometry against exact answers: eight fcc systems share a Schmid "
        "factor of 1/sqrt(6) under [001] tension, the cube orientation's "
        "full-constraint Taylor factor is exactly sqrt(6), and a random fcc "
        "texture averages Taylor's 1938 value of 3.06."
    ),
    examples=(
        CUBE_SCHMID_FACTOR,
        CUBE_TAYLOR_FACTOR,
        RANDOM_TEXTURE_TAYLOR_FACTOR,
    ),
)
