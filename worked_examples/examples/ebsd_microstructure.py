"""Worked examples: microstructural quantities from an orientation map.

These examples connect the EBSD layer to a physical, dimensioned quantity. A
lattice bent by a known amount contains a calculable density of geometrically
necessary dislocations, and Nye's relation fixes that density in closed form —
so the expected value comes from the planted curvature and Burgers vector, not
from a prior run of the code.

See :doc:`../../workflows/ebsd_kam` and the GND theory note.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

TILTED_MAP_SETUP = """
import numpy as np
from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    geometrically_necessary_dislocation_density,
    lattice_curvature_tensor,
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
    axes=("RD", "TD", "ND"),
    handedness=Handedness.RIGHT,
)
map_frame = ReferenceFrame(
    name="map",
    domain=FrameDomain.MAP,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
copper = Phase(
    "copper",
    lattice=Lattice(3.615, 3.615, 3.615, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)

# A 9 x 11 map, 0.5 um step, whose lattice rotates about the specimen normal by
# a known 0.8 degrees per micrometre along x. The gradient is planted through
# Rotation.from_axis_angle, independently of the curvature code.
ROWS, COLS, STEP_UM, GRADIENT_DEG_PER_UM = 9, 11, 0.5, 0.8
BURGERS_NM = 0.2556  # copper, a / sqrt(2)

coordinates = []
quaternions = []
for row in range(ROWS):
    for col in range(COLS):
        x_um, y_um = col * STEP_UM, row * STEP_UM
        coordinates.append((x_um, y_um))
        angle = np.deg2rad(GRADIENT_DEG_PER_UM * x_um)
        quaternions.append(Rotation.from_axis_angle([0.0, 0.0, 1.0], angle).quaternion)

bent_map = CrystalMap(
    coordinates=np.asarray(coordinates, dtype=np.float64),
    orientations=OrientationSet.from_quaternions(
        np.asarray(quaternions), specimen_frame=specimen, phase=copper
    ),
    map_frame=map_frame,
    grid_shape=(ROWS, COLS),
    step_sizes=(STEP_UM, STEP_UM),
)
"""

_KAPPA = SymbolUse(r"\kappa_{ij}", "Lattice curvature tensor, in radians per metre.")
_RHO_GND = SymbolUse(
    r"\rho_{\mathrm{GND}}",
    "Geometrically necessary dislocation density, in reciprocal square metres.",
)
_BURGERS = SymbolUse("b", "Burgers vector magnitude.")

_KAM_WORKFLOW = SeeAlso("EBSD local misorientation", "../../workflows/ebsd_kam")
_EBSD_CONCEPT = SeeAlso("EBSD foundation", "../../concepts/ebsd_foundation")
_EBSD_GRAPH_THEORY = SeeAlso(
    "Multiphase and hex-grid EBSD graphs",
    "../../theory/multiphase_ebsd_graph_workflows",
)


HEX_GRID_KAM = WorkedExample(
    id="ebsd-hex-grid-six-neighbor-kam",
    title="Six-neighbour KAM on a staggered hexagonal scan",
    domain="ebsd",
    scenario=(
        "EDAX/TSL .ang scans can alternate long and short rows. Treating those points as a "
        "rectangle invents pixels and changes every local average, so the reader preserves the "
        "3/2/3 row topology and KAM uses the six-neighbour graph directly."
    ),
    setup="""
from pathlib import Path
import numpy as np
from pytex.adapters import read_ang

hex_map = read_ang(Path("fixtures/ebsd/synthetic_hex_grid.ang")).crystal_map
""",
    code=(
        "graph = hex_map.neighbor_graph()\n"
        "kam = hex_map.kernel_average_misorientation_deg()\n"
        "result = np.concatenate([[len(graph.pairs)], kam])"
    ),
    expected=[13.0, 3.0, 6.0, 3.0, 1.2, 1.2, 0.0, 0.0, 0.0],
    unit="degree (KAM entries)",
    tolerance=2e-7,
    reference=(
        "The synthetic scan has five horizontal edges and four edges across each of its two "
        "staggered row boundaries: 5 + 4 + 4 = 13. Only point 1 is rotated by 6 degrees. Its "
        "four incident edges therefore give KAM 6 degrees; points 0 and 2 average one 6-degree "
        "edge over degree 2 (3 degrees), points 3 and 4 over degree 5 (1.2 degrees), and the "
        "remaining points see only zero-angle edges."
    ),
    citation=(
        "MTEX, 'Gridded EBSD Data', Hexagonal Grids section, documenting EBSDhex as the distinct "
        "hexagonal-grid topology and warning that square resampling distorts the unit cells: "
        "https://mtex-toolbox.github.io/EBSDGrid.html. Fixture metadata states explicitly that "
        "the values are analytically constructed and not experimental."
    ),
    see_also=(_KAM_WORKFLOW, _EBSD_GRAPH_THEORY, _EBSD_CONCEPT),
    result_format="{:.6f}",
)


PLANTED_LATTICE_CURVATURE = WorkedExample(
    id="ebsd-planted-lattice-curvature",
    title="Lattice curvature recovered from a known orientation gradient",
    domain="ebsd",
    scenario=(
        "Before trusting a curvature map on real data it is worth confirming that the "
        "measurement recovers a gradient that was put there deliberately. A lattice rotated about "
        "the specimen normal by a constant 0.8 degrees per micrometre along x has exactly one "
        "nonzero curvature component, kappa_20 = d(omega_z)/dx, and its value follows from unit "
        "conversion alone."
    ),
    setup=TILTED_MAP_SETUP,
    code=(
        "curvature = lattice_curvature_tensor(bent_map)\n"
        "result = float(curvature[4, 5, 2, 0])"
    ),
    expected=13962.634015954637,
    unit="rad/m",
    tolerance=1e-6,
    reference=(
        "0.8 degrees per micrometre is 0.8 * pi / 180 radians per 1e-6 metre, that is "
        "0.013962634015954637 / 1e-6 = 13962.634015954637 rad/m. The identity is exact; the "
        "tolerance is numerical only. Note that the third column of the curvature tensor is NaN "
        "throughout, because a surface map cannot measure the out-of-plane gradient."
    ),
    citation=(
        "Nye, J. F., Acta Metall. 1, 153-162 (1953), DOI: 10.1016/0001-6160(53)90054-6; "
        "Pantleon, W., Scripta Mater. 58, 994-997 (2008), "
        "DOI: 10.1016/j.scriptamat.2008.01.050."
    ),
    symbols=(_KAPPA,),
    see_also=(_KAM_WORKFLOW, _EBSD_CONCEPT),
    result_format="{:.6f}",
)


GND_DENSITY_FROM_CURVATURE = WorkedExample(
    id="ebsd-gnd-density-from-curvature",
    title="GND density of a bent copper lattice",
    domain="ebsd",
    scenario=(
        "A bent lattice must contain dislocations of one sign to accommodate the bending, and "
        "their density follows from the curvature and the Burgers vector. This is the quantity "
        "that links an orientation map to stored energy and work hardening. For a pure tilt the "
        "Nye relation collapses to rho = (d(theta)/dx) / b, so the answer is available in closed "
        "form and the map's own number can be checked against it."
    ),
    setup=TILTED_MAP_SETUP,
    code=(
        "density = geometrically_necessary_dislocation_density(\n"
        "    bent_map, burgers_vector_nm=BURGERS_NM\n"
        ")\n"
        "result = float(density[4, 5])"
    ),
    expected=5.4626893255690666e13,
    unit="1/m^2",
    tolerance=1e6,
    reference=(
        "For a single-axis tilt the only nonzero Nye component is alpha_02 = kappa_20, so "
        "rho = kappa_20 / b = 13962.634015954637 / (0.2556e-9) = 5.46269e13 m^-2. That lands in "
        "the 1e13 to 1e14 range expected of lightly deformed copper. The value is a lower bound "
        "and is resolution dependent: content producing no in-plane curvature is invisible to a "
        "surface map, and a finer step would resolve sharper gradients and report more."
    ),
    citation=(
        "Nye, J. F., Acta Metall. 1, 153-162 (1953), DOI: 10.1016/0001-6160(53)90054-6; "
        "Kysar, J. W. et al., Int. J. Plasticity 26, 1097-1123 (2010), "
        "DOI: 10.1016/j.ijplas.2010.03.009."
    ),
    symbols=(_KAPPA, _RHO_GND, _BURGERS),
    see_also=(_KAM_WORKFLOW, _EBSD_CONCEPT),
    result_format="{:.6e}",
)


GROUP = ExampleGroup(
    slug="ebsd",
    title="EBSD microstructure",
    summary=(
        "Hex-grid KAM, lattice curvature, and geometrically necessary dislocation density "
        "recovered from analytically planted topology or orientation gradients."
    ),
    examples=(
        HEX_GRID_KAM,
        PLANTED_LATTICE_CURVATURE,
        GND_DENSITY_FROM_CURVATURE,
    ),
)

__all__ = ["GROUP"]
