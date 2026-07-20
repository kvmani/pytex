"""Worked examples: composite OR SAED simulation.

These examples validate the two numerical cornerstones of the composite
diffraction surface: the relativistic electron wavelength that fixes the
Ewald-sphere radius, and the exactness of the Kurdjumov-Sachs child-zone
mapping (the defining direction parallelism reproduced by the variant
machinery to machine precision).

See ``docs/roadmap/working_notes_composite_saed_program.md`` and the OR
concept documentation.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

WAVELENGTH_SETUP = """
from pytex.diffraction.kinematic import electron_wavelength_angstrom
"""

COMPOSITE_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

parent_frame = ReferenceFrame(
    "parent_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
child_frame = ReferenceFrame(
    "child_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
martensite = Phase(
    "martensite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=martensite
)
"""

_LAMBDA = SymbolUse(r"\lambda", "Radiation wavelength.")

_DIFF_CONCEPT = SeeAlso("Diffraction foundation", "../../concepts/diffraction_foundation")
_OR_CONCEPT = SeeAlso(
    "Orientation relationships", "../../concepts/orientation_relationships"
)


ELECTRON_WAVELENGTH_200KV = WorkedExample(
    id="composite-electron-wavelength-200kv",
    title="Relativistic electron wavelength at 200 kV",
    domain="diffraction",
    scenario=(
        "Every kinematic TEM computation starts from the electron wavelength, which fixes the "
        "Ewald-sphere radius k = 1/lambda and hence every excitation error. The relativistic "
        "formula lambda = h / sqrt(2 m0 e V (1 + e V / (2 m0 c^2))) must reproduce the standard "
        "tabulated value at a 200 kV accelerating voltage."
    ),
    setup=WAVELENGTH_SETUP,
    code="result = electron_wavelength_angstrom(200.0)",
    expected=0.02508,
    unit="angstrom",
    tolerance=5e-6,
    reference=(
        "The standard relativistic electron wavelength at 200 kV is 2.508 pm = 0.02508 angstrom."
    ),
    citation=(
        "De Graef, Introduction to Conventional Transmission Electron Microscopy, "
        "Cambridge University Press, 2003, Table 2.2."
    ),
    symbols=(_LAMBDA,),
    see_also=(_DIFF_CONCEPT,),
    result_format="{:.5f}",
)


KS_EXACT_CHILD_ZONE = WorkedExample(
    id="composite-ks-exact-child-zone",
    title="KS maps the parent [0 1 -1] zone exactly onto a <1 1 1> child zone",
    domain="diffraction",
    scenario=(
        "The Kurdjumov-Sachs relationship is defined by the parallelism <-1 0 1>_fcc || "
        "<-1 -1 1>_bcc. When the composite SAED simulator maps a parent [0 1 -1] zone axis "
        "(a member of the <-1 0 1> family) through all 24 variants, at least one variant's "
        "child zone axis must land exactly on a <1 1 1>-type direction: the minimal angular "
        "deviation between mapped and rational child zones over the variants is zero."
    ),
    setup=COMPOSITE_SETUP,
    code=(
        "zone = ZoneAxis(np.array([0, 1, -1]), phase=austenite)\n"
        "composite = simulate_composite_saed(ks, zone, include_parent=False)\n"
        "result = min(\n"
        "    pattern.nearest_zone_axis.deviation_deg\n"
        "    for pattern in composite.variant_patterns\n"
        ")"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "The defining KS direction parallelism makes the mapped child zone rational, so the "
        "deviation of the best variant is exactly 0 degrees."
    ),
    citation=(
        "Kurdjumov and Sachs, Z. Physik 64 (1930) 325; Morito et al., Acta Materialia 51 "
        "(2003) 1789 (variant conventions)."
    ),
    symbols=(),
    see_also=(_OR_CONCEPT, _DIFF_CONCEPT),
)


BURGERS_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

beta_frame = ReferenceFrame(
    "beta_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
alpha_frame = ReferenceFrame(
    "alpha_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
# Beta-titanium (bcc) and alpha-titanium (hcp), room-temperature parameters.
beta_ti = Phase(
    "beta-titanium",
    lattice=Lattice(3.3065, 3.3065, 3.3065, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
    crystal_frame=beta_frame,
)
alpha_ti = Phase(
    "alpha-titanium",
    lattice=Lattice(2.9508, 2.9508, 4.6855, 90.0, 90.0, 120.0, crystal_frame=alpha_frame),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=alpha_frame),
    crystal_frame=alpha_frame,
)
burgers = OrientationRelationship.from_burgers_correspondence(
    parent_phase=beta_ti, child_phase=alpha_ti
)
"""

_BURGERS_CONCEPT = SeeAlso(
    "Composite OR diffraction workflow", "../../workflows/composite_or_diffraction"
)


BURGERS_EXACT_BASAL_ZONE = WorkedExample(
    id="composite-burgers-exact-basal-zone",
    title="Burgers maps the parent <110> zone exactly onto the hcp [0001] basal zone",
    domain="diffraction",
    scenario=(
        "The Burgers relationship governing the beta->alpha transformation of titanium, "
        "zirconium and hafnium is defined by the plane parallelism {110}_bcc || (0001)_hcp. "
        "Viewing a beta crystal down a <110> zone axis must therefore look straight down the "
        "hcp c-axis for the variants whose basal plane is that particular {110}: the minimal "
        "angular deviation between the mapped child zone and a rational [0001] zone must be "
        "exactly zero."
    ),
    setup=BURGERS_SETUP,
    code=(
        "zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)\n"
        "composite = simulate_composite_saed(burgers, zone, include_parent=False)\n"
        "result = min(\n"
        "    pattern.nearest_zone_axis.deviation_deg\n"
        "    for pattern in composite.variant_patterns\n"
        ")"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "The defining Burgers plane parallelism {110}_bcc || (0001)_hcp makes the mapped child "
        "zone exactly rational, so the deviation of the best variant is 0 degrees."
    ),
    citation="Burgers, Physica 1 (1934) 561.",
    symbols=(),
    see_also=(_OR_CONCEPT, _BURGERS_CONCEPT),
)


BURGERS_BASAL_COINCIDENCE = WorkedExample(
    id="composite-burgers-110-0002-coincidence",
    title="Burgers {110}_bcc and (0002)_hcp reflections nearly superimpose",
    domain="diffraction",
    scenario=(
        "The practical TEM signature of the Burgers relationship is that the beta {110} "
        "reflection lands almost exactly on the alpha (0002) reflection, because the plane "
        "parallelism pairs two nearly equal interplanar spacings: d(110)_bcc = a/sqrt(2) = "
        "2.3381 angstrom against d(0002)_hcp = c/2 = 2.3428 angstrom. At a 180 mm*angstrom "
        "camera constant the residual detector separation is well under a spot diameter, so "
        "the composite pattern reads as a single decorated pattern. This computes that "
        "separation from the simulated composite."
    ),
    setup=BURGERS_SETUP,
    code=(
        "from pytex.diffraction.composite import find_spot_coincidences\n"
        "\n"
        "zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)\n"
        "composite = simulate_composite_saed(burgers, zone)\n"
        "report = find_spot_coincidences(composite, tolerance_mm=1.0)\n"
        "result = report.coincidences[0].separation_mm"
    ),
    expected=0.15450,
    unit="mm",
    tolerance=1e-4,
    reference=(
        "Analytically the separation is (sqrt(2)/a_bcc - 2/c_hcp) * camera_constant = "
        "(1.414214/3.3065 - 2/4.6855) * 180 = 0.15450 mm."
    ),
    citation="Burgers, Physica 1 (1934) 561; lattice parameters from standard Ti data.",
    symbols=(),
    see_also=(_BURGERS_CONCEPT, _DIFF_CONCEPT),
    result_format="{:.5f}",
)


GROUP = ExampleGroup(
    slug="composite-diffraction",
    title="Composite OR diffraction",
    summary=(
        "Numerical cornerstones of composite orientation-relationship SAED simulation: the "
        "relativistic electron wavelength against the standard 200 kV value, the exactness of "
        "the Kurdjumov-Sachs child-zone mapping, and the two defining Burgers beta->alpha "
        "signatures (exact basal zone and the {110}/(0002) near-coincidence)."
    ),
    examples=(
        ELECTRON_WAVELENGTH_200KV,
        KS_EXACT_CHILD_ZONE,
        BURGERS_EXACT_BASAL_ZONE,
        BURGERS_BASAL_COINCIDENCE,
    ),
)

__all__ = ["GROUP"]
