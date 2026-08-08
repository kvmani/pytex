"""Worked examples: pole-figure resampling, m.r.d. scale, and arithmetic.

Two pole figures cannot be compared until they share a support and a scale.
Both steps are pinned here by exact identities rather than by tolerances: a
weighted-mean interpolator reproduces a constant field exactly, an m.r.d.
figure has unit mean density by definition, and densities add.

See :doc:`../../concepts/orientation_texture` for the surrounding ODF doctrine.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

POLE_FIGURE_SETUP = """
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    S2Grid,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("x", "y", "z"), handedness=Handedness.RIGHT
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
pole = CrystalPlane(MillerIndex((1, 1, 1), phase=phase), phase=phase)
grid = S2Grid.equispaced(10.0, reference_frame=specimen, hemisphere="upper", antipodal=True)
"""

_POLE_DENSITY = SymbolUse(
    r"P_{hkl}(\mathbf{y})",
    "Pole density of the family in multiples of a random distribution.",
)
_WEIGHTS = SymbolUse(
    r"w_i",
    "Solid-angle integration weight of sampled direction i; the weights sum to one.",
)
_DIFFERENCE = SymbolUse(
    r"\Delta P",
    "Signed difference of two pole densities on a shared support.",
)

_TEXTURE_CONCEPT = SeeAlso("Orientations and texture", "../../concepts/orientation_texture")
_API = SeeAlso("Texture API", "../../api/index")


MRD_NORMALIZATION = WorkedExample(
    id="texture-pole-figure-mrd-unit-mean-density",
    title="A pole figure in m.r.d. has unit mean density, and its deviation from random sums to zero",
    domain="texture",
    scenario=(
        "Resample a two-component texture onto an equal-area grid with the "
        "default m.r.d. normalization, then verify the two identities that "
        "define the scale: the solid-angle-weighted mean density is exactly "
        "one, and the weighted integral of the deviation from random is "
        "exactly zero. Neither holds for a figure normalized by its maximum "
        "or by its sum, which is why those scales cannot be compared between "
        "measurements."
    ),
    setup=POLE_FIGURE_SETUP,
    code=(
        "orientations = OrientationSet.from_euler_angles(\n"
        "    np.array([[0.0, 0.0, 0.0], [35.0, 20.0, 10.0]]),\n"
        "    specimen_frame=specimen,\n"
        "    phase=phase,\n"
        ")\n"
        "figure = PoleFigure.from_orientations(orientations, pole).on_grid(\n"
        "    grid, halfwidth_deg=15.0\n"
        ")\n"
        "mean_density = float(np.sum(grid.weights * figure.intensities))\n"
        "deviation = figure - 1.0\n"
        "deviation_integral = float(np.sum(grid.weights * deviation.values))\n"
        "result = np.array([mean_density, deviation_integral])"
    ),
    expected=[1.0, 0.0],
    unit="m.r.d.",
    tolerance=1e-12,
    reference=(
        "Definition of the multiples-of-random scale: sum_i w_i P_i = 1 with "
        "solid-angle weights summing to one. The deviation identity follows "
        "immediately, since sum_i w_i (P_i - 1) = 1 - 1 = 0."
    ),
    citation=(
        "Bunge, Texture Analysis in Materials Science (1982), normalization of "
        "pole figures to multiples of a random distribution."
    ),
    symbols=(_POLE_DENSITY, _WEIGHTS, _DIFFERENCE),
    see_also=(_TEXTURE_CONCEPT, _API),
    result_format="{:.12f}",
)


RESAMPLING_AND_ADDITION = WorkedExample(
    id="texture-pole-figure-resampling-and-addition-identities",
    title="Resampling preserves a constant field, and pole densities add",
    domain="texture",
    scenario=(
        "Resample a pole figure that is 2.5 m.r.d. everywhere onto a coarser "
        "grid, and add two normalized figures. The interpolating estimator is "
        "a weighted mean, so a constant field passes through it unchanged for "
        "any kernel halfwidth — the partition-of-unity property that "
        "distinguishes it from the summing estimator used for pole clouds. "
        "Densities then add pointwise, so the mean of a sum of two m.r.d. "
        "figures is exactly two."
    ),
    setup=POLE_FIGURE_SETUP,
    code=(
        "flat = PoleFigure(\n"
        "    pole=pole,\n"
        "    sample_directions=grid.vectors.values,\n"
        "    intensities=np.full(len(grid), 2.5),\n"
        "    specimen_frame=specimen,\n"
        "    antipodal=True,\n"
        "    sampling='sampled_density',\n"
        ")\n"
        "coarse = S2Grid.equispaced(\n"
        "    18.0, reference_frame=specimen, hemisphere='upper', antipodal=True\n"
        ")\n"
        "resampled = flat.on_grid(coarse, halfwidth_deg=7.0, normalize=False)\n"
        "constant_field = float(np.max(np.abs(resampled.intensities - 2.5)))\n"
        "\n"
        "single = PoleFigure.from_orientations(\n"
        "    OrientationSet.from_euler_angles(\n"
        "        np.zeros((1, 3)), specimen_frame=specimen, phase=phase\n"
        "    ),\n"
        "    pole,\n"
        ").on_grid(grid, halfwidth_deg=15.0)\n"
        "total = single + single\n"
        "sum_mean = float(np.sum(grid.weights * total.intensities))\n"
        "result = np.array([constant_field, sum_mean])"
    ),
    expected=[0.0, 2.0],
    unit="m.r.d.",
    tolerance=1e-12,
    reference=(
        "Partition of unity of the Nadaraya-Watson estimator: sum_i K_i f_i / "
        "sum_i K_i = c whenever every f_i = c, for any kernel and any query "
        "direction. Linearity of the solid-angle mean then gives "
        "mean(P + P) = 2 for a figure of unit mean."
    ),
    citation=(
        "Nadaraya (1964), Theory of Probability and its Applications 9:141, "
        "and Watson (1964), Sankhya A 26:359, for the weighted-mean estimator."
    ),
    symbols=(_POLE_DENSITY, _WEIGHTS),
    see_also=(_TEXTURE_CONCEPT, _API),
    result_format="{:.12f}",
)


GROUP = ExampleGroup(
    slug="pole-figure-arithmetic",
    title="Pole-figure arithmetic",
    summary=(
        "Exact identities behind comparing two pole figures: the "
        "multiples-of-random scale, resampling onto a shared support, and "
        "the additivity of pole densities."
    ),
    examples=(MRD_NORMALIZATION, RESAMPLING_AND_ADDITION),
)
