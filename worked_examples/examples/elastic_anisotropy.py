"""Worked examples: cubic elastic anisotropy against closed forms.

Three checks, all against results derived rather than recorded:

* the directional Young's modulus of a cubic crystal depends on direction
  only through J = n1^2 n2^2 + n2^2 n3^2 + n3^2 n1^2, so [110] and [112] -
  which share J = 1/4 - are exactly equally stiff, for every cubic material;
* for a randomly textured cubic aggregate the Voigt and Reuss *bulk* moduli
  are identical, so the whole bound gap lives in the shear modulus;
* the numerically homogenized aggregate reproduces the closed-form Voigt and
  Reuss shear moduli.

Copper single-crystal constants (Simmons and Wang) are used throughout.

See :doc:`../../theory/elastic_anisotropy_and_homogenization`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

ELASTIC_SETUP = """
import numpy as np
from pytex.properties import StiffnessTensor

# Copper single-crystal stiffness, GPa.
C11, C12, C44 = 168.4, 121.4, 75.4
stiffness = StiffnessTensor.cubic(C11, C12, C44)
compliance = stiffness.compliance()
voigt = np.asarray(compliance.voigt_matrix())
S11, S12, S44 = voigt[0, 0], voigt[0, 1], voigt[3, 3]


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)
"""

AGGREGATE_SETUP = ELASTIC_SETUP + """
from pytex.core import (
    FrameDomain,
    OrientationSet,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.properties import homogenize_elastic

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)

# Haar-uniform orientations (Shoemake), so the aggregate is texture-free.
rng = np.random.default_rng(17)
u1, u2, u3 = rng.random(40000), rng.random(40000), rng.random(40000)
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
"""

_E = SymbolUse(
    r"E(\hat{\mathbf{n}})",
    "Young's modulus along the crystal direction n.",
)
_K = SymbolUse(r"K_{V}, K_{R}", "Voigt and Reuss bulk moduli of an aggregate.")
_MU = SymbolUse(r"\mu_{V}, \mu_{R}", "Voigt and Reuss aggregate shear moduli.")

_THEORY = SeeAlso(
    "Elastic anisotropy and homogenization",
    "../../theory/elastic_anisotropy_and_homogenization",
)


YOUNGS_MODULUS_EXTREMES = WorkedExample(
    id="elastic-cubic-youngs-modulus-110-equals-112",
    title="[110] and [112] are exactly equally stiff in any cubic crystal",
    domain="properties",
    scenario=(
        "Evaluate Young's modulus of copper along [100], [110], [112] and "
        "[111]. For cubic symmetry the direction enters only through "
        "J = n1^2 n2^2 + n2^2 n3^2 + n3^2 n1^2, which is 0 along <100>, 1/3 "
        "along <111>, and exactly 1/4 along both [110] and [112]. Those two "
        "directions therefore have identical stiffness - not approximately, "
        "and not only for copper. The example returns all four moduli so the "
        "extremes and the coincidence are checked together."
    ),
    setup=ELASTIC_SETUP,
    code=(
        "directions = [\n"
        "    [1.0, 0.0, 0.0],\n"
        "    [1.0, 1.0, 0.0],\n"
        "    [1.0, 1.0, 2.0],\n"
        "    [1.0, 1.0, 1.0],\n"
        "]\n"
        "result = np.array(\n"
        "    [float(compliance.youngs_modulus(unit(d))) for d in directions]\n"
        ")"
    ),
    expected=[66.6888, 130.3376, 130.3376, 191.1497],
    unit="GPa",
    tolerance=1e-3,
    reference=(
        "Closed form 1/E = S11 - 2(S11 - S12 - S44/2) J with J = 0, 1/4, 1/4, "
        "1/3 respectively; S11, S12, S44 from the analytic cubic inverse of "
        "(C11, C12, C44)."
    ),
    citation=(
        "Nye, Physical Properties of Crystals (OUP); Simmons and Wang, Single "
        "Crystal Elastic Constants (MIT Press 1971) for the copper constants."
    ),
    symbols=(_E,),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


CUBIC_BULK_MODULUS_BOUNDS_COINCIDE = WorkedExample(
    id="elastic-cubic-voigt-reuss-bulk-moduli-coincide",
    title="The Voigt and Reuss bulk moduli of a cubic aggregate are equal",
    domain="properties",
    scenario=(
        "Compute the Voigt and Reuss bulk moduli of a randomly textured cubic "
        "aggregate from the closed forms and take their difference. It is "
        "exactly zero: a cubic crystal responds isotropically to hydrostatic "
        "pressure, so the uniform-stress and uniform-strain assumptions cannot "
        "disagree about dilatation. The whole Voigt-Reuss gap therefore lives "
        "in the shear modulus, and reporting a 'Hill bulk modulus' for a cubic "
        "aggregate implies an uncertainty that does not exist. The example "
        "returns the difference K_V - K_R, which must vanish."
    ),
    setup=ELASTIC_SETUP,
    code=(
        "k_voigt = (C11 + 2.0 * C12) / 3.0\n"
        "k_reuss = 1.0 / (3.0 * (S11 + 2.0 * S12))\n"
        "result = float(k_voigt - k_reuss)"
    ),
    expected=0.0,
    unit="GPa",
    tolerance=1e-9,
    reference=(
        "Analytic identity: K_V = (C11 + 2 C12)/3 and K_R = 1/(3(S11 + 2 S12)) "
        "are equal for cubic symmetry, since S11 + 2 S12 = 3/(C11 + 2 C12)."
    ),
    citation=(
        "Hill, The elastic behaviour of a crystalline aggregate, Proc. Phys. "
        "Soc. A 65 (1952) 349-354, DOI 10.1088/0370-1298/65/5/307."
    ),
    symbols=(_K,),
    see_also=(_THEORY,),
    result_format="{:.3e}",
)


AGGREGATE_MATCHES_CLOSED_FORM = WorkedExample(
    id="elastic-random-aggregate-matches-voigt-reuss-closed-form",
    title="A homogenized random aggregate reproduces the Voigt and Reuss shear moduli",
    domain="properties",
    scenario=(
        "Homogenize copper over 40000 Haar-random orientations under the Voigt "
        "and Reuss schemes and read the aggregate shear modulus back as C44 of "
        "the averaged tensor, comparing with the closed forms "
        "mu_V = (C11 - C12 + 3 C44)/5 and "
        "mu_R = 5/(4(S11 - S12) + 3 S44). This exercises the rank-four "
        "rotation, the weighted average, and the compliance inversion "
        "together. The residual is finite-sample texture in the random "
        "orientation set, not an error in the averaging."
    ),
    setup=AGGREGATE_SETUP,
    code=(
        "moduli = []\n"
        "for scheme in ('voigt', 'reuss'):\n"
        "    aggregate = homogenize_elastic(\n"
        "        stiffness, orientations, scheme=scheme\n"
        "    )\n"
        "    moduli.append(float(np.asarray(aggregate.voigt_matrix())[3, 3]))\n"
        "result = np.array(moduli)"
    ),
    expected=[54.6400, 40.0339],
    unit="GPa",
    tolerance=0.05,
    reference=(
        "Closed forms for a randomly textured cubic aggregate: "
        "mu_V = (C11 - C12 + 3 C44)/5 = 54.6400 GPa and "
        "mu_R = 5/(4(S11 - S12) + 3 S44) = 40.0339 GPa. Tolerance covers the "
        "finite-sample texture of 40000 random orientations."
    ),
    citation=(
        "Hill, Proc. Phys. Soc. A 65 (1952) 349-354; Simmons and Wang, Single "
        "Crystal Elastic Constants (MIT Press 1971)."
    ),
    symbols=(_MU,),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


GROUP = ExampleGroup(
    slug="elastic-anisotropy",
    title="Elastic anisotropy and homogenization",
    summary=(
        "Cubic elasticity against closed forms: [110] and [112] are exactly "
        "equally stiff, the Voigt and Reuss bulk moduli of a cubic aggregate "
        "are identical so the entire bound gap is in the shear modulus, and a "
        "numerically homogenized random aggregate reproduces both shear bounds."
    ),
    examples=(
        YOUNGS_MODULUS_EXTREMES,
        CUBIC_BULK_MODULUS_BOUNDS_COINCIDE,
        AGGREGATE_MATCHES_CLOSED_FORM,
    ),
)
