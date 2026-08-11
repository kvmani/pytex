"""Worked examples: why a pole-figure raster needs solid-angle weights.

A tilt/rotation raster is a regular grid in polar angle and azimuth, which is
*not* a uniform sampling of the sphere: the same azimuthal point count crowds
into a vanishing ring at the pole and spreads over a full circle at the
equator. Averaging such a raster without weights is a mean with respect to
d(psi) d(phi) rather than sin(psi) d(psi) d(phi), and over-counts the pole.

Both averages of cos^2(psi) over a hemisphere are elementary - 1/2 unweighted
and 1/3 correct - so the bias is exactly +50 percent, and it is a bias rather
than a discretisation error: refining the raster does not reduce it. These
examples compute both, at two raster steps, against those closed forms.

See :doc:`../../theory/pole_figure_arithmetic_and_mrd`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

RASTER_SETUP = """
import numpy as np
from pytex.core.sphere import raster_solid_angle_weights


def raster_polar_angles(step_deg):
    polar = np.arange(0.0, 90.0 + step_deg, step_deg)
    azimuth = np.arange(0.0, 360.0, step_deg)
    grid, _ = np.meshgrid(polar, azimuth, indexing="ij")
    return grid.ravel()
"""

_WEIGHT = SymbolUse(
    r"w_i",
    "Solid-angle integration weight of sampled direction i; weights sum to 1 "
    "over the sampled region.",
)
_PSI = SymbolUse(
    r"\psi", "Polar (tilt) angle of a specimen direction from the specimen +Z axis."
)

_THEORY = SeeAlso(
    "Pole-figure arithmetic and the m.r.d. scale",
    "../../theory/pole_figure_arithmetic_and_mrd",
)


UNWEIGHTED_RASTER_MEAN_IS_BIASED = WorkedExample(
    id="pole-figure-raster-unweighted-mean-is-biased",
    title="An unweighted raster mean of cos^2 is 1/2, not 1/3",
    domain="texture",
    scenario=(
        "Average cos^2(psi) over a 5 degree tilt/rotation raster without "
        "weights. Because the raster is uniform in polar angle rather than in "
        "solid angle, the result is the elementary integral of cos^2 over "
        "[0, pi/2] divided by pi/2, which is exactly 1/2 - while the true "
        "spherical mean is 1/3. The unweighted answer is 3/2 times the correct "
        "one, a 50 percent error, and it is a bias rather than a "
        "discretisation error: the companion example shows it is unchanged by "
        "refining the raster."
    ),
    setup=RASTER_SETUP,
    code=(
        "polar = raster_polar_angles(5.0)\n"
        "field = np.cos(np.radians(polar)) ** 2\n"
        "result = float(field.mean())"
    ),
    expected=0.5,
    unit="",
    tolerance=1e-9,
    reference=(
        "Analytic: an unweighted mean over a raster uniform in psi evaluates "
        "the integral of cos^2(psi) d(psi) over [0, pi/2] divided by pi/2, "
        "which is exactly 1/2. The correct spherical mean is 1/3."
    ),
    citation=(
        "Randle and Engler, Introduction to Texture Analysis - measured "
        "pole-figure rasters and their sampling geometry."
    ),
    symbols=(_PSI,),
    see_also=(_THEORY,),
    result_format="{:.6f}",
)


REFINING_THE_RASTER_DOES_NOT_HELP = WorkedExample(
    id="pole-figure-raster-bias-survives-refinement",
    title="Halving the raster step leaves the unweighted error at exactly 50 percent",
    domain="texture",
    scenario=(
        "Repeat the unweighted average at a 2.5 degree step, half the previous "
        "one. The answer is again exactly 1/2. This is the signature of a "
        "biased estimator: a discretisation error would halve, while a bias "
        "from integrating against the wrong measure is independent of "
        "resolution. Buying a finer scan cannot fix it; only weighting can. "
        "The example returns the ratio of the fine-step mean to the "
        "coarse-step mean, which must be exactly 1."
    ),
    setup=RASTER_SETUP,
    code=(
        "coarse = np.cos(np.radians(raster_polar_angles(5.0))) ** 2\n"
        "fine = np.cos(np.radians(raster_polar_angles(2.5))) ** 2\n"
        "result = float(fine.mean() / coarse.mean())"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-9,
    reference=(
        "Analytic: both unweighted means equal 1/2 exactly, independent of the "
        "raster step, so their ratio is 1. A discretisation error would not "
        "behave this way."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982) - normalization "
        "of pole figures as spherical integrals."
    ),
    symbols=(_PSI,),
    see_also=(_THEORY,),
    result_format="{:.9f}",
)


WEIGHTED_RASTER_MEAN_CONVERGES = WorkedExample(
    id="pole-figure-raster-weighted-mean-converges",
    title="Solid-angle weights recover the spherical mean, and refine correctly",
    domain="texture",
    scenario=(
        "Average the same field using raster_solid_angle_weights, which give "
        "each ring the solid angle of the band midway to its neighbours, "
        "cos(lower) - cos(upper), shared among its points. The weighted mean "
        "approaches the exact 1/3, and unlike the unweighted one it improves "
        "with resolution: the error falls from 4.1 percent at 5 degrees to 2.1 "
        "percent at 2.5 degrees. The example returns both weighted means so "
        "the convergence is visible."
    ),
    setup=RASTER_SETUP,
    code=(
        "means = []\n"
        "for step in (5.0, 2.5):\n"
        "    polar = raster_polar_angles(step)\n"
        "    weights = raster_solid_angle_weights(polar)\n"
        "    field = np.cos(np.radians(polar)) ** 2\n"
        "    means.append(float((weights * field).sum()))\n"
        "result = np.array(means)"
    ),
    expected=[0.31960, 0.32627],
    unit="",
    tolerance=1e-4,
    reference=(
        "Converging on the exact spherical mean 1/3 = 0.33333 from below, with "
        "the ring-band quadrature error halving as the step halves."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982) - spherical "
        "integration of pole densities."
    ),
    symbols=(_WEIGHT, _PSI),
    see_also=(_THEORY,),
    result_format="{:.5f}",
)


GROUP = ExampleGroup(
    slug="pole-figure-sampling",
    title="Pole-figure raster sampling",
    summary=(
        "Why a measured pole figure needs solid-angle weights: the unweighted "
        "mean of cos^2 over a tilt raster is exactly 1/2 against a true "
        "spherical mean of 1/3, a 50 percent bias that survives halving the "
        "raster step, while the weighted mean converges."
    ),
    examples=(
        UNWEIGHTED_RASTER_MEAN_IS_BIASED,
        REFINING_THE_RASTER_DOES_NOT_HELP,
        WEIGHTED_RASTER_MEAN_CONVERGES,
    ),
)
