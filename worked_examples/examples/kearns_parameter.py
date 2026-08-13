"""Worked examples: the Kearns parameter and the identities that calibrate it.

The Kearns orientation parameter ``f`` is the volume-weighted mean of ``cos^2``
of the angle between each crystal's basal pole and a specimen direction, and
the effective fraction of basal poles aligned with it. Written as the quadratic
form of the pole orientation tensor ``A = <c c^T>``, two facts the zirconium
literature states empirically become identities: the values along an
orthonormal triad sum to ``tr(A) = 1`` for every texture, and a random texture
gives ``1/3`` in every direction.

Those two, plus Kearns' own tabulated calculation, are what these examples
check. None of the expected values is a prior output of this code: three are
elementary consequences of the definition and one is the number printed in
Kearns' 1965 report.

See :doc:`../../theory/kearns_parameter_and_basal_pole_texture`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

ZIRCONIUM_SETUP = """
import numpy as np
from scipy.spatial.transform import Rotation

from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.core.lattice import CrystalPlane
from pytex.texture.kearns import kearns_from_orientations, kearns_from_tilt_profile

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
SPECIMEN = ReferenceFrame(
    "sample_rd_td_nd", FrameDomain.SPECIMEN, ("RD", "TD", "ND"), Handedness.RIGHT
)
SYMMETRY = SymmetrySpec.from_point_group("6/mmm", reference_frame=CRYSTAL)
LATTICE = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=CRYSTAL)
ZIRCONIUM = Phase("alpha_zr", lattice=LATTICE, symmetry=SYMMETRY, crystal_frame=CRYSTAL)
BASAL = CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=ZIRCONIUM)


def zirconium_orientations(matrices):
    return OrientationSet.from_matrices(
        np.asarray(matrices, dtype=np.float64).reshape(-1, 3, 3),
        crystal_frame=CRYSTAL,
        specimen_frame=SPECIMEN,
        phase=ZIRCONIUM,
        symmetry=SYMMETRY,
    )
"""

_F = SymbolUse(
    r"f",
    "Kearns orientation parameter along a specimen direction: the volume-weighted "
    "mean of cos^2 of the angle between each crystal's basal pole and that "
    "direction. Exactly 1/3 for a random texture.",
)
_TENSOR = SymbolUse(
    r"\mathbf{A}",
    "Pole orientation tensor of the basal-pole distribution in the specimen frame, "
    "with f(d) = d^T A d. Unit trace, which is why an orthonormal triad's Kearns "
    "parameters sum identically to 1.",
)
_PHI = SymbolUse(
    r"\phi",
    "Tilt of a crystal's basal pole [0001] from the specimen reference direction.",
)
_TILT_PROFILE = SymbolUse(
    r"I(\phi)",
    "Basal-pole density averaged over the full 360 degrees of azimuth about the "
    "reference direction; Kearns' Eq. (5) integrates it against sin(phi) cos^2(phi).",
)

_THEORY = SeeAlso(
    "The Kearns parameter and basal-pole texture",
    "../../theory/kearns_parameter_and_basal_pole_texture",
)
_TUTORIAL = SeeAlso(
    "Tutorial 31: Kearns parameter estimation",
    "../../tutorials/notebooks/31_kearns_parameter",
)


RANDOM_TEXTURE_GIVES_ONE_THIRD = WorkedExample(
    id="kearns-random-texture-is-one-third",
    title="A random texture gives f = 1/3 along every direction",
    domain="texture",
    scenario=(
        "Resolve the basal poles of 20000 uniformly random alpha-zirconium "
        "orientations onto the specimen axes. A texture-free aggregate sends "
        "basal poles uniformly over the sphere, so the pole orientation tensor "
        "is the isotropic I/3 and f is 1/3 in every direction. This is the null "
        "hypothesis every measured f is quoted against, and the value that fixes "
        "the scale: 0 means no basal poles along the direction, 1/3 means "
        "random, 1 means all of them. The example returns the largest departure "
        "from 1/3 over the three axes, which is orientation-sampling error alone."
    ),
    setup=ZIRCONIUM_SETUP,
    code=(
        "orientations = zirconium_orientations(Rotation.random(20000, random_state=7).as_matrix())\n"
        "report = kearns_from_orientations(orientations, pole=BASAL)\n"
        "result = float(np.max(np.abs(np.asarray(report.values) - 1.0 / 3.0)))"
    ),
    expected=0.0,
    unit="",
    tolerance=6e-3,
    reference=(
        "Analytic: a uniform distribution of unit vectors has second-moment "
        "tensor I/3, so f = d.(I/3).d = 1/3 for every unit d. The residual is "
        "Monte-Carlo error in the 20000-orientation sample, falling as 1/sqrt(n)."
    ),
    citation=(
        "J. J. Kearns, Thermal Expansion and Preferred Orientation in Zircaloy, "
        "WAPD-TM-472 (1965), Section II: 'a value of 1/3 in each direction "
        "defines the isotropic case'."
    ),
    symbols=(_F, _TENSOR),
    see_also=(_THEORY, _TUTORIAL),
    result_format="{:.6f}",
)


TRIAD_SUM_IS_EXACTLY_ONE = WorkedExample(
    id="kearns-triad-sum-is-exactly-one",
    title="f_RD + f_TD + f_ND = 1 exactly, for any texture and any triad",
    domain="texture",
    scenario=(
        "Take a deliberately awkward texture -- a lopsided mixture of three "
        "unrelated orientation clusters, nothing like a fibre -- and sum its "
        "Kearns parameters over the specimen axes. The answer is 1, not "
        "approximately 1: the sum is the trace of the pole orientation tensor, "
        "and every basal pole is a unit vector, so the trace is 1 whatever the "
        "distribution. This is the single most useful check available on a "
        "Kearns measurement, because it needs no reference specimen: a measured "
        "triad that misses 1 is reporting the systematic error of the "
        "measurement -- unmeasured tilt range, a wrong random standard, an "
        "unbalanced background -- and nothing about the material."
    ),
    setup=ZIRCONIUM_SETUP,
    code=(
        "clusters = Rotation.from_euler(\n"
        '    "zxz",\n'
        "    [[13.0, 71.0, 5.0], [200.0, 12.0, 87.0], [95.0, 140.0, 31.0]],\n"
        "    degrees=True,\n"
        ").as_matrix()\n"
        "orientations = zirconium_orientations(clusters)\n"
        "report = kearns_from_orientations(orientations, pole=BASAL, weights=[5.0, 2.0, 1.0])\n"
        "result = float(report.triad_sum)"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "Analytic: sum over an orthonormal triad of d.A.d equals tr(A) = "
        "<c.c> = 1, since every basal pole c is a unit vector. Exact for every "
        "texture and every orthonormal triad, so the tolerance is floating-point "
        "rather than statistical."
    ),
    citation=(
        "J. J. Kearns, WAPD-TM-472 (1965), Section II: 'the sum of f in the "
        "three principal directions must be unity'."
    ),
    symbols=(_F, _TENSOR),
    see_also=(_THEORY, _TUTORIAL),
    result_format="{:.12f}",
)


IDEAL_BASAL_GIRDLE = WorkedExample(
    id="kearns-ideal-basal-girdle",
    title="A basal girdle in the RD-TD plane gives f = (1/2, 1/2, 0)",
    domain="texture",
    scenario=(
        "Spread the basal poles uniformly around the RD-TD great circle -- the "
        "limiting case of the transverse texture that pilgering and extrusion "
        "produce in zirconium tubing, with every c axis in the plane and none "
        "along the normal. The Kearns parameters follow from the mean of cos^2 "
        "around a circle: 1/2 in the two in-plane directions and exactly 0 along "
        "the normal, summing to 1 as they must. It is worth contrasting with the "
        "single crystal, which gives (0, 0, 1): both are extreme textures, but a "
        "girdle can never exceed 1/2 in any direction. The example returns the "
        "three values."
    ),
    setup=ZIRCONIUM_SETUP,
    code=(
        "azimuth = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)\n"
        "c_axes = np.stack(\n"
        "    [np.cos(azimuth), np.sin(azimuth), np.zeros_like(azimuth)], axis=1\n"
        ")\n"
        "first = np.cross(np.array([0.0, 0.0, 1.0]), c_axes)\n"
        "second = np.cross(c_axes, first)\n"
        "matrices = np.stack([first, second, c_axes], axis=2)\n"
        "report = kearns_from_orientations(zirconium_orientations(matrices), pole=BASAL)\n"
        "result = np.asarray(report.values)"
    ),
    expected=[0.5, 0.5, 0.0],
    unit="",
    tolerance=1e-10,
    reference=(
        "Analytic: the mean of cos^2 of the angle to an in-plane axis, over a "
        "uniform great circle, is 1/2; the angle to the normal is 90 degrees "
        "everywhere, so that component is 0."
    ),
    citation=(
        "J. L. Baron et al., Textures and Microstructures 12 (1990) 125-140, "
        "doi:10.1155/TSM.12.125 -- the T-type (transverse) texture of zircaloy "
        "tubing, whose basal poles concentrate in the plane perpendicular to the "
        "tube axis."
    ),
    symbols=(_F, _TENSOR),
    see_also=(_THEORY, _TUTORIAL),
    result_format="{:.6f}",
)


KEARNS_TABLE_3_LONGITUDINAL = WorkedExample(
    id="kearns-1965-table-3-longitudinal-section",
    title="Kearns' own tabulated calculation reproduces at f = 0.488",
    domain="texture",
    scenario=(
        "Integrate Kearns' Equation (5) over the basal-pole tilt profile he "
        "tabulates for the longitudinal section of a swaged Zircaloy-2 rod: the "
        "azimuthally averaged (0001) pole density at the midpoints of ten-degree "
        "tilt bins. The sin(phi) factor in the integrand converts pole density "
        "to volume fraction -- the band of orientations at tilt phi has "
        "circumference proportional to sin(phi) -- which is why the profile at "
        "high tilt matters even where the density is small, and why the volume "
        "fraction vanishes at the figure centre however intense the pole is "
        "there. The result must be the 0.488 printed in his Table 3."
    ),
    setup=ZIRCONIUM_SETUP,
    code=(
        "tilt_deg = np.arange(5.0, 90.0, 10.0)\n"
        "density = np.array([3.27, 2.71, 1.69, 1.35, 1.17, 0.97, 0.73, 0.62, 0.55])\n"
        "report = kearns_from_tilt_profile(\n"
        "    tilt_deg, density, pole=BASAL, specimen_frame=SPECIMEN\n"
        ")\n"
        'result = report.value("ND")'
    ),
    expected=0.488,
    unit="",
    tolerance=1e-3,
    reference=(
        "Kearns (1965) Table 3, longitudinal-section block: his I_phi column is "
        "the input above and his tabulated total is f = 0.488. The tolerance is "
        "the rounding of his three-significant-figure entries, not a fitted "
        "margin. His transverse-section block does not reproduce -- its 70-80 "
        "degree row lists 0.0214 where 0.353 cos^2(75 deg) = 0.0237, so the "
        "quoted 0.0508 should be 0.0532 -- which is why the longitudinal block "
        "is the one pinned here."
    ),
    citation=(
        "J. J. Kearns, Thermal Expansion and Preferred Orientation in Zircaloy, "
        "WAPD-TM-472, Bettis Atomic Power Laboratory (1965), Eq. (5) and Table 3."
    ),
    symbols=(_F, _PHI, _TILT_PROFILE),
    see_also=(_THEORY, _TUTORIAL),
    result_format="{:.4f}",
)


GROUP = ExampleGroup(
    slug="kearns-parameter",
    title="The Kearns parameter",
    summary=(
        "The scalar texture index the zirconium industry specifies components "
        "against, checked against the identities that calibrate it -- 1/3 for a "
        "random texture, an exact sum of 1 over any orthonormal triad, "
        "(1/2, 1/2, 0) for an ideal basal girdle -- and against the tabulated "
        "calculation in Kearns' own 1965 report."
    ),
    examples=(
        RANDOM_TEXTURE_GIVES_ONE_THIRD,
        TRIAD_SUM_IS_EXACTLY_ONE,
        IDEAL_BASAL_GIRDLE,
        KEARNS_TABLE_3_LONGITUDINAL,
    ),
)
