"""Worked examples: the random-disorientation baseline against closed forms.

An MDF only means something relative to what a texture-free aggregate would
have produced, so the baseline is the part that has to be right. Two of the
three checks here are against analytic results rather than against another
simulation:

* with no symmetry the disorientation-angle density is exactly
  (1 - cos w)/pi, whose mean is pi/2 + 2/pi = 126.4756 degrees;
* the largest possible cubic disorientation is a property of the Rodrigues
  fundamental zone, 2*arctan(sqrt(23 - 16*sqrt(2))) = 62.7994 degrees.

The third records the cubic mean, which has no elementary closed form and is
the number most often quoted loosely - the mean is 40.7 degrees while the
median is 42.3, and they are not interchangeable.

See :doc:`../../theory/random_disorientation_baseline` for the derivations.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

MDF_SETUP = """
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec
from pytex.core.misorientation_distribution import (
    random_disorientation_angles_deg,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


def sampled_angles(point_group, total, chunk=20000, seed0=100):
    # The reduction materialises n*|G|^2 matrices, so a large baseline is
    # generated in chunks and concatenated; the samples are independent.
    symmetry = SymmetrySpec.from_point_group(
        point_group, reference_frame=crystal
    )
    parts = [
        random_disorientation_angles_deg(symmetry, chunk, seed=seed0 + i)
        for i in range(total // chunk)
    ]
    return np.concatenate(parts)
"""

_OMEGA = SymbolUse(
    r"\omega",
    "Rotation (disorientation) angle of a misorientation.",
)
_RHO = SymbolUse(
    r"\boldsymbol{\rho}",
    "Rodrigues vector n tan(omega/2); the chart in which the cubic "
    "fundamental zone is a cube intersected with an octahedron.",
)

_THEORY = SeeAlso(
    "Random disorientation baseline", "../../theory/random_disorientation_baseline"
)
_ZONE = SeeAlso(
    "Orientation space and disorientation",
    "../../theory/orientation_space_and_disorientation",
)


TRICLINIC_MEAN_IS_ANALYTIC = WorkedExample(
    id="mdf-triclinic-mean-disorientation-angle",
    title="With no symmetry the mean disorientation angle is pi/2 + 2/pi",
    domain="orientation",
    scenario=(
        "Sample disorientation angles between independent uniformly (Haar) "
        "distributed orientations with triclinic symmetry, so no reduction "
        "occurs, and compare the mean against the exact value. The Haar "
        "measure in axis-angle coordinates carries the factor (1 - cos w), "
        "which is the volume of the shell of rotations at angle w, so the "
        "density is (1 - cos w)/pi and its first moment integrates to "
        "pi/2 + 2/pi = 126.4756 degrees. This is the strongest available check "
        "on the sampler because the target is analytic rather than another "
        "simulation."
    ),
    setup=MDF_SETUP,
    code=(
        "angles = sampled_angles('1', 200000, seed0=200)\n"
        "result = float(angles.mean())"
    ),
    expected=126.4756,
    unit="deg",
    tolerance=0.25,
    reference=(
        "Analytic: integral of w*(1 - cos w)/pi over [0, pi] equals "
        "pi/2 + 2/pi radians = 126.4756 degrees. Tolerance covers the "
        "Monte-Carlo standard error at n = 2e5."
    ),
    citation=(
        "Morawiec, Orientations and Rotations (Springer 2004) - the invariant "
        "measure on SO(3) in axis-angle coordinates."
    ),
    symbols=(_OMEGA,),
    see_also=(_THEORY, _ZONE),
    result_format="{:.4f}",
)


CUBIC_MAXIMUM_IS_A_ZONE_VERTEX = WorkedExample(
    id="mdf-cubic-maximum-disorientation-angle",
    title="The largest cubic disorientation is 2*arctan(sqrt(23-16*sqrt2))",
    domain="orientation",
    scenario=(
        "The maximum cubic disorientation is not a sampling outcome but a "
        "property of the Rodrigues fundamental zone, which for cubic symmetry "
        "is the cube |rho_i| <= sqrt(2)-1 intersected with the octahedron "
        "sum|rho_i| <= 1. The angle increases with |rho|, so the maximum sits "
        "at the vertex farthest from the origin, (sqrt2-1, sqrt2-1, 3-2sqrt2), "
        "which meets both constraints with equality. Its magnitude is "
        "sqrt(23 - 16 sqrt2), giving 62.7994 degrees about <1, 1, sqrt2-1>. "
        "This example evaluates that closed form; a sampled maximum converges "
        "to it only slowly from below, which is why the exact value is the one "
        "to quote."
    ),
    setup=MDF_SETUP,
    code=(
        "rho_max = np.array(\n"
        "    [np.sqrt(2) - 1, np.sqrt(2) - 1, 3 - 2 * np.sqrt(2)]\n"
        ")\n"
        "result = float(\n"
        "    np.degrees(2 * np.arctan(np.linalg.norm(rho_max)))\n"
        ")"
    ),
    expected=62.7994,
    unit="deg",
    tolerance=1e-3,
    reference=(
        "Analytic: |rho| = sqrt(23 - 16*sqrt(2)) at the cubic fundamental-zone "
        "vertex, so omega = 2*arctan(sqrt(23 - 16*sqrt(2))) = 62.79943 degrees."
    ),
    citation=(
        "Mackenzie, Second paper on statistics associated with the random "
        "disorientation of cubes, Biometrika 45 (1958) 229-240, "
        "DOI 10.1093/biomet/45.1-2.229."
    ),
    symbols=(_RHO, _OMEGA),
    see_also=(_THEORY, _ZONE),
    result_format="{:.4f}",
)


CUBIC_LOW_ANGLE_FALSE_POSITIVE_RATE = WorkedExample(
    id="mdf-cubic-random-low-angle-fraction",
    title="2.2 percent of random cubic boundaries are low-angle by chance",
    domain="orientation",
    scenario=(
        "Compute the fraction of random cubic disorientations falling below "
        "the conventional 15 degree low-angle threshold. This is the null "
        "hypothesis a low-angle-boundary fraction has to beat: a texture-free "
        "aggregate already delivers about 2.2 percent low-angle boundaries "
        "from geometry alone, so a map reporting 3 percent has demonstrated "
        "essentially nothing while one reporting 30 percent has. The value is "
        "recomputed rather than tabulated because it depends on the point "
        "group."
    ),
    setup=MDF_SETUP,
    code=(
        "angles = sampled_angles('m-3m', 300000, seed0=300)\n"
        "result = float((angles < 15.0).mean())"
    ),
    expected=0.0223,
    unit="",
    tolerance=0.002,
    reference=(
        "Monte-Carlo estimate over 3e5 independent Haar-random misorientations "
        "reduced by the 24 proper cubic rotations; cross-checked against an "
        "independently written quaternion implementation. Tolerance covers the "
        "binomial standard error and the seed."
    ),
    citation=(
        "Randle and Engler, Introduction to Texture Analysis - the "
        "conventional 15 degree low-angle boundary threshold."
    ),
    symbols=(_OMEGA,),
    see_also=(_THEORY, _ZONE),
    result_format="{:.4f}",
)


GROUP = ExampleGroup(
    slug="random-disorientation",
    title="Random disorientation baseline",
    summary=(
        "The null hypothesis every MDF claim is measured against: the exact "
        "(1 - cos w)/pi density with no symmetry, the exact cubic maximum "
        "62.7994 degrees from the Rodrigues zone vertex, and the 2.2 percent "
        "of random cubic boundaries that are low-angle by chance."
    ),
    examples=(
        TRICLINIC_MEAN_IS_ANALYTIC,
        CUBIC_MAXIMUM_IS_A_ZONE_VERTEX,
        CUBIC_LOW_ANGLE_FALSE_POSITIVE_RATE,
    ),
)
