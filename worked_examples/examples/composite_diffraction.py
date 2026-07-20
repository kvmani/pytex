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


GROUP = ExampleGroup(
    slug="composite-diffraction",
    title="Composite OR diffraction",
    summary=(
        "Numerical cornerstones of composite orientation-relationship SAED simulation: the "
        "relativistic electron wavelength against the standard 200 kV value, and the exactness "
        "of the Kurdjumov-Sachs child-zone mapping."
    ),
    examples=(ELECTRON_WAVELENGTH_200KV, KS_EXACT_CHILD_ZONE),
)

__all__ = ["GROUP"]
