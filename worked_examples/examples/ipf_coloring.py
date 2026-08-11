"""Worked examples: how an inverse-pole-figure colour is actually computed.

IPF colouring is the most reproduced image in texture analysis and the least
specified: textbooks show the coloured triangle without stating the function
that carries a direction to a colour. These examples pin that function down
against values derived by hand rather than recorded from a previous run.

Three identities are checked. The sector corners must take the primaries
exactly, because that is what makes the key readable. A general direction
must take the colour that the closed form

    beta = (dz - dx, sqrt(2)(dx - dy), sqrt(3) dy)
    c_i  = (beta_i / max_j beta_j) ** (1 / gamma_s)

predicts, which for [113] is exactly (1, 0, 3/4). And every symmetric
equivalent of a direction must take the *same* colour, since otherwise the map
would show contrast where there is no crystallography.

See :doc:`../../theory/ipf_color_keys` for the derivation.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

IPF_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    IPFColorKey,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction="z")
"""

_BETA = SymbolUse(
    r"\boldsymbol{\beta}",
    "Barycentric weights of a direction in the fundamental-sector corner basis.",
)
_GAMMA_S = SymbolUse(
    r"\gamma_{s}",
    "IPF saturation parameter; channels are raised to the power 1/gamma_s.",
)

_THEORY = SeeAlso("IPF colour keys", "../../theory/ipf_color_keys")
_WORKFLOW = SeeAlso("IPF colour workflow", "../../workflows/ipf_colors")


SECTOR_CORNERS_ARE_PRIMARIES = WorkedExample(
    id="ipf-cubic-sector-corners-are-primaries",
    title="The cubic sector corners colour to exactly red, green and blue",
    domain="visualization",
    scenario=(
        "Colour the three corners of the cubic standard triangle - [001], "
        "[101] and [111]. Each corner has barycentric weights equal to a "
        "standard basis vector, so after the saturation power and the "
        "max-channel renormalization it must return exactly one primary. This "
        "is the identity that makes an IPF legend readable without consulting "
        "a lookup table, and it is what fails first if the sector corners or "
        "the colour basis are mis-ordered."
    ),
    setup=IPF_SETUP,
    code=(
        "corners = np.array(\n"
        "    [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]\n"
        ")\n"
        "result = key.colors_from_crystal_directions(corners)"
    ),
    expected=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    unit="",
    tolerance=1e-12,
    reference=(
        "Analytic identity: a sector corner has barycentric weights equal to a "
        "standard basis vector, which the colour map carries to the "
        "corresponding primary for any saturation exponent."
    ),
    citation=(
        "Nolze and Hielscher, Orientations - perfectly colored, J. Appl. "
        "Cryst. 49 (2016) 1786-1802, DOI 10.1107/S1600576716012942."
    ),
    symbols=(_BETA,),
    see_also=(_THEORY, _WORKFLOW),
    result_format="{:.6f}",
)


CLOSED_FORM_COLOUR_113 = WorkedExample(
    id="ipf-cubic-closed-form-colour-113",
    title="The [113] direction colours to exactly (1, 0, 3/4)",
    domain="visualization",
    scenario=(
        "Colour a direction lying on the [001]-[111] edge of the cubic "
        "triangle, where the whole chain can be followed by hand. For "
        "[113]/sqrt(11) the closed-form weights are beta = (2, 0, sqrt(3))/"
        "sqrt(11); the largest is beta_1, so the colour is "
        "(1, 0, sqrt(3)/2) raised to the power 1/gamma_s = 2, giving exactly "
        "(1, 0, 3/4). Agreement here exercises the symmetry reduction, the "
        "barycentric solve, the saturation power and the renormalization "
        "together, against a number derived rather than recorded."
    ),
    setup=IPF_SETUP,
    code=(
        "direction = np.array([[1.0, 1.0, 3.0]])\n"
        "direction = direction / np.linalg.norm(direction)\n"
        "result = key.colors_from_crystal_directions(direction)[0]"
    ),
    expected=[1.0, 0.0, 0.75],
    unit="",
    tolerance=1e-12,
    reference=(
        "Closed form: beta = (dz - dx, sqrt(2)(dx - dy), sqrt(3) dy) evaluated "
        "at (1,1,3)/sqrt(11) gives (2, 0, sqrt(3))/sqrt(11), and "
        "(beta / max beta) ** 2 = (1, 0, 3/4) exactly."
    ),
    citation=(
        "International Tables for Crystallography, Vol. A - the m-3m "
        "asymmetric unit that fixes the sector corners."
    ),
    symbols=(_BETA, _GAMMA_S),
    see_also=(_THEORY, _WORKFLOW),
    result_format="{:.6f}",
)


SYMMETRIC_EQUIVALENTS_SHARE_A_COLOUR = WorkedExample(
    id="ipf-symmetric-equivalents-share-one-colour",
    title="All 24 cubic equivalents of a direction take one colour",
    domain="visualization",
    scenario=(
        "Generate every symmetric equivalent of a general direction under "
        "m-3m and colour them all. The spread across the orbit must be zero: "
        "symmetry-equivalent directions are the same physical direction, so a "
        "colouring that separated them would paint contrast where there is no "
        "crystallography, and the picture would depend on which equivalent "
        "index a file happened to store. The example reports the maximum "
        "channel spread over the orbit, which is bounded by the rotation "
        "arithmetic rather than by a tolerance in the colouring."
    ),
    setup=IPF_SETUP,
    code=(
        "direction = np.array([0.3, 0.1, 0.9])\n"
        "direction = direction / np.linalg.norm(direction)\n"
        "orbit = symmetry.equivalent_vectors(direction)\n"
        "orbit = np.asarray(\n"
        "    orbit.values if hasattr(orbit, 'values') else orbit\n"
        ").reshape(-1, 3)\n"
        "colors = key.colors_from_crystal_directions(orbit)\n"
        "result = float(np.abs(colors - colors[0]).max())"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "Analytic identity: the colour is a function of the symmetry-reduced "
        "direction alone, so it is constant on a symmetry orbit by "
        "construction. The expected spread is exactly zero."
    ),
    citation=(
        "Nolze and Hielscher, Orientations - perfectly colored, J. Appl. "
        "Cryst. 49 (2016) 1786-1802, DOI 10.1107/S1600576716012942."
    ),
    symbols=(_BETA,),
    see_also=(_THEORY, _WORKFLOW),
    result_format="{:.2e}",
)


GROUP = ExampleGroup(
    slug="ipf-coloring",
    title="Inverse-pole-figure colouring",
    summary=(
        "What an IPF colour actually is, checked against hand-derived values: "
        "the sector corners colour to exact primaries, a direction on the "
        "[001]-[111] edge colours to exactly (1, 0, 3/4) by the closed form, "
        "and every symmetric equivalent shares one colour."
    ),
    examples=(
        SECTOR_CORNERS_ARE_PRIMARIES,
        CLOSED_FORM_COLOUR_113,
        SYMMETRIC_EQUIVALENTS_SHARE_A_COLOUR,
    ),
)
