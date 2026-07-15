"""Worked examples: orientation representations and disorientation angles.

These examples exercise the ``Orientation`` object: the round-trip consistency
of its representation conversions and the symmetry-reduced disorientation angle
between two orientations. The reference values are exact identities (a round
trip must return the input) and a textbook coincidence-site value (the Sigma 3
twin).

See :doc:`../../concepts/orientation_texture` and the theory note on
orientation space and disorientation.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

ORIENTATION_SETUP = """
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cube = Orientation.from_euler(
    0.0, 0.0, 0.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic
)
"""

_OMEGA = SymbolUse(r"\omega", "Disorientation angle: minimal misorientation angle over symmetry.")
_EULER = SymbolUse(r"(\phi_1, \Phi, \phi_2)", "Bunge Euler angles.")
_ROT = SymbolUse(r"\mathbf{R}", "Active rotation matrix.")

_TEXTURE_CONCEPT = SeeAlso("Orientations and texture", "../../concepts/orientation_texture")
_ORI_API = SeeAlso("Orientation constructors and helpers", "../../api/index")


EULER_MATRIX_ROUNDTRIP = WorkedExample(
    id="orientation-euler-matrix-roundtrip",
    title="Euler -> matrix -> Euler round trip returns the input orientation",
    domain="orientation",
    scenario=(
        "Before trusting any orientation pipeline you must confirm that converting Bunge Euler "
        "angles to a rotation matrix and back is lossless. Here we take a non-degenerate "
        "orientation, convert to a matrix, rebuild an Orientation from that matrix, and measure the "
        "disorientation to the original. A correct implementation returns exactly zero."
    ),
    setup=ORIENTATION_SETUP,
    code=(
        "g = Orientation.from_euler(\n"
        "    30.0, 40.0, 50.0, specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic\n"
        ")\n"
        "rebuilt = Orientation.from_matrix(\n"
        "    g.as_matrix(), specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic\n"
        ")\n"
        "result = float(np.degrees(g.distance_to(rebuilt)))"
    ),
    expected=0.0,
    unit="deg",
    tolerance=1e-9,
    reference="A representation round trip is an identity map; the disorientation must be 0 degrees.",
    citation="Bunge, Texture Analysis in Materials Science, 1982, Chapter 2.",
    symbols=(_EULER, _ROT, _OMEGA),
    see_also=(_TEXTURE_CONCEPT, _ORI_API),
)

SIGMA3_DISORIENTATION = WorkedExample(
    id="orientation-sigma3-disorientation",
    title="Sigma 3 twin disorientation is 60 degrees about <111>",
    domain="orientation",
    scenario=(
        "The coherent annealing twin in FCC metals is the Sigma 3 boundary, a 60-degree rotation "
        "about a <111> axis. Computing the symmetry-reduced disorientation between the cube "
        "orientation and its 60-degree/<111> partner is the standard validation that cubic "
        "crystal symmetry is applied correctly during misorientation reduction."
    ),
    setup=ORIENTATION_SETUP,
    code=(
        "twin = Orientation.from_axis_angle(\n"
        "    (1, 1, 1), np.radians(60.0),\n"
        "    specimen_frame=specimen, symmetry=cubic.symmetry, phase=cubic,\n"
        ")\n"
        "result = float(np.degrees(cube.distance_to(twin)))"
    ),
    expected=60.0,
    unit="deg",
    tolerance=1e-6,
    reference="The Sigma 3 coincidence-site boundary is a 60-degree rotation about <111>.",
    citation="Grimmer, Bollmann and Warrington, Acta Cryst. A30 (1974) 197; Randle, The Role of the CSL.",
    symbols=(_OMEGA,),
    see_also=(
        _TEXTURE_CONCEPT,
        SeeAlso("Symmetry and fundamental regions", "../../concepts/symmetry_and_fundamental_regions"),
    ),
)


GROUP = ExampleGroup(
    slug="orientation",
    title="Orientations and disorientation angles",
    summary=(
        "Round-trip consistency of orientation representations and symmetry-reduced disorientation "
        "angles, checked against exact identities and the Sigma 3 twin reference."
    ),
    examples=(
        EULER_MATRIX_ROUNDTRIP,
        SIGMA3_DISORIENTATION,
    ),
)

__all__ = ["GROUP"]
