"""Worked examples: simulated SAED plates and the zone-axis atlas.

Every value checked here is fixed by the lattice and the camera constant alone,
which is the point: a practice pattern is only useful if the geometry on it is
the geometry a microscope would record. The calibration identity ``r = L*lambda
/ d`` is what turns a picture into a measurement; the hcp prism-zone aspect ratio
is a calibration-free measurement of ``c/a``; and the basal-to-prism angle is
exactly 90 degrees in any hexagonal lattice whatever the axial ratio.

See :doc:`../../workflows/tem_pattern_indexing` for the workflow these support,
and :doc:`../../theory/saed_ratio_angle_indexing` for the indexing method.
"""

from __future__ import annotations

import math

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

PATTERN_SETUP = """
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
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Aluminium, a = 4.0495 A (Wyckoff, Crystal Structures Vol. 1).
aluminium = Phase(
    "aluminium-fcc",
    lattice=Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# Alpha zirconium, a = 3.232 A, c = 5.147 A, so c/a = 1.5925.
zirconium = Phase(
    "zirconium-hcp",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
# 200 kV and a 400 mm camera length: L*lambda = 10.0317 mm.A.
CAMERA_CONSTANT = 10.031748206
RASTER = DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=0.024)
"""

_D = SymbolUse("d", "Interplanar spacing of a reflecting plane.")
_G = SymbolUse("g", "Reciprocal-lattice vector of a reflection; |g| = 1/d.")
_LAMBDA = SymbolUse(r"\lambda", "Electron wavelength at the accelerating voltage.")
_ZONE = SymbolUse(r"[uvw]", "Lattice direction brought parallel to the beam.")

_WORKFLOW = SeeAlso("TEM pattern indexing workflow", "../../workflows/tem_pattern_indexing")
_INDEXING = SeeAlso("Ratio and angle indexing", "../../theory/saed_ratio_angle_indexing")
_NAVIGATION = SeeAlso(
    "TEM specimen tilt navigation", "../../theory/tem_specimen_tilt_navigation"
)


CALIBRATION_IDENTITY = WorkedExample(
    id="saed-practice-camera-constant-identity",
    title="Where the 200 reflection of aluminium lands on the detector",
    domain="tem",
    scenario=(
        "The single identity a diffraction pattern is calibrated by: a reflection sits at a "
        "distance from the transmitted beam equal to the camera constant divided by its "
        "d-spacing. Everything downstream — the indexed answer, the lattice parameter, the phase "
        "identification — inherits this one relation, which is why a camera constant taken from "
        "the wrong camera length produces a self-consistent pattern of the wrong material. Here "
        "the simulated plate is asked where its strongest reflection is, and the answer must be "
        "the one the definition gives."
    ),
    setup=PATTERN_SETUP,
    code=(
        "image = synthesize_saed_image(\n"
        "    aluminium,\n"
        "    ZoneAxis([0, 0, 1], phase=aluminium),\n"
        "    camera_constant_mm_angstrom=CAMERA_CONSTANT,\n"
        "    raster=RASTER,\n"
        ")\n"
        "spot = next(\n"
        "    entry\n"
        "    for entry in image.spots\n"
        "    if sorted(abs(int(value)) for value in entry.miller_indices) == [0, 0, 2]\n"
        ")\n"
        "centre = np.asarray(image.centre_px)\n"
        "radius_px = float(np.linalg.norm(np.asarray(spot.position_px) - centre))\n"
        "result = radius_px * image.raster.pixel_size_mm"
    ),
    expected=10.031748206 / (4.0495 / 2.0),
    unit="mm",
    tolerance=1e-9,
    reference=(
        "r = (L*lambda) / d with d_200 = a/2 = 2.02475 A for a = 4.0495 A, giving "
        "r = 10.031748206 / 2.02475 = 4.95456 mm. Analytic from the lattice parameter "
        "and the definition of the camera constant; no program output enters it."
    ),
    citation=(
        "Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, 2nd ed., "
        "Springer, DOI: 10.1007/978-0-387-76501-3, chapter 18 (the camera equation "
        "R d = L lambda); Wyckoff, R. W. G., Crystal Structures Vol. 1 (1963) for a."
    ),
    symbols=(_D, _G, _LAMBDA),
    see_also=(_WORKFLOW, _INDEXING),
    result_format="{:.5f}",
)


PRISM_ZONE_AXIAL_RATIO = WorkedExample(
    id="saed-practice-hcp-prism-axial-ratio",
    title="Reading c/a off one hcp prism-zone pattern",
    domain="tem",
    scenario=(
        "The hcp [2-1-10] pattern is a rectangle whose two shortest vectors are 0002 along c* and "
        "01-10 perpendicular to it. Their lengths are 2/c and 2/(sqrt(3) a), so their ratio is "
        "sqrt(3) a / c and depends on nothing else — not on the camera constant, not on the "
        "accelerating voltage, not on the exposure. That makes this one pattern a "
        "calibration-free measurement of the axial ratio, and the standard way to tell zirconium "
        "(1.0876) from titanium (1.0908) or magnesium (1.0668) on the microscope. Here the ratio "
        "is measured off the simulated plate exactly as it would be measured off a real one."
    ),
    setup=PATTERN_SETUP,
    code=(
        "image = synthesize_saed_image(\n"
        "    zirconium,\n"
        "    ZoneAxis([1, 0, 0], phase=zirconium),\n"
        "    camera_constant_mm_angstrom=CAMERA_CONSTANT,\n"
        "    raster=RASTER,\n"
        ")\n"
        "by_indices = {\n"
        "    tuple(int(value) for value in spot.miller_indices): spot for spot in image.spots\n"
        "}\n"
        "result = by_indices[(0, 0, 2)].g_inv_angstrom / by_indices[(0, 1, 0)].g_inv_angstrom"
    ),
    expected=math.sqrt(3.0) * 3.232 / 5.147,
    unit="",
    tolerance=1e-9,
    reference=(
        "|g_0002| / |g_01-10| = (2/c) / (2 / (sqrt(3) a)) = sqrt(3) a / c = "
        "sqrt(3) * 3.232 / 5.147 = 1.08757, from the hexagonal reciprocal metric alone. "
        "Independent of the camera constant, which cancels in the ratio."
    ),
    citation=(
        "Edington, J. W., Practical Electron Microscopy in Materials Science, Macmillan "
        "(1975), on interpreting hexagonal zone-axis patterns; lattice parameters from "
        "Wyckoff, R. W. G., Crystal Structures Vol. 1 (1963)."
    ),
    symbols=(_D, _G, _ZONE),
    see_also=(_WORKFLOW, _INDEXING),
    result_format="{:.5f}",
)


BASAL_TO_PRISM_ANGLE = WorkedExample(
    id="saed-practice-atlas-basal-to-prism",
    title="Basal to prism is 90 degrees for every hexagonal metal",
    domain="tem",
    scenario=(
        "The zone-axis atlas exists to answer 'where should I go next', and the first thing it has "
        "to get right is how far away each candidate is. The hexagonal basal-to-prism pair is the "
        "cleanest possible check: the c axis is perpendicular to every a axis by the definition of "
        "the hexagonal cell, so the angle is exactly 90 degrees whatever the axial ratio, whatever "
        "the metal. It is also the practical lesson the pair teaches, because 90 degrees is beyond "
        "any conventional double-tilt holder in one move."
    ),
    setup=PATTERN_SETUP,
    code=(
        "from pytex.tem.atlas import zone_axis_atlas\n"
        "\n"
        "atlas = zone_axis_atlas(\n"
        "    zirconium,\n"
        "    current_zone_axis=ZoneAxis([0, 0, 1], phase=zirconium),\n"
        "    max_index=1,\n"
        ")\n"
        "prism = next(entry for entry in atlas.entries if entry.label == '[100]')\n"
        "result = prism.angle_from_current_deg"
    ),
    expected=90.0,
    unit="deg",
    tolerance=1e-9,
    reference=(
        "In a hexagonal lattice alpha = beta = 90 degrees, so c is orthogonal to both a "
        "axes by construction of the cell. The angle between [0001] and [2-1-10] is "
        "therefore exactly 90 degrees for any c/a."
    ),
    citation=(
        "International Tables for Crystallography, Volume A, on the hexagonal cell "
        "setting; Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, "
        "2nd ed., DOI: 10.1007/978-0-387-76501-3, chapter 18."
    ),
    symbols=(_ZONE,),
    see_also=(_WORKFLOW, _NAVIGATION),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="saed_practice_patterns",
    title="Simulated SAED plates and the zone-axis atlas",
    summary=(
        "The geometry a practice diffraction pattern must reproduce if indexing it is to teach "
        "anything: the camera-constant identity that places every reflection, the hcp prism-zone "
        "aspect ratio that measures c/a without any calibration at all, and the basal-to-prism "
        "angle the zone-axis atlas has to report as exactly 90 degrees."
    ),
    examples=(
        CALIBRATION_IDENTITY,
        PRISM_ZONE_AXIAL_RATIO,
        BASAL_TO_PRISM_ANGLE,
    ),
)

__all__ = ["GROUP"]
