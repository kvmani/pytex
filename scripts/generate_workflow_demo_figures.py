"""Regenerate the demo figures embedded in the workflow and concept pages.

``docs/figures/crystal_visualization_demo.png``, ``powder_xrd_demo.svg``, and
``saed_demo.svg`` are the first picture a reader sees on the crystal-visualization,
XRD, and SAED pages. They were committed as one-off exports and drifted away from
what the library actually renders: the NaCl figure had lost its per-element
colours, bonds, and unit cell and showed only overlapping grey blobs.

This script regenerates all three from the public API, so the pages show the
current renderer and the figures can be refreshed whenever it improves::

    python scripts/generate_workflow_demo_figures.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from pytex import (
    AtomicSite,
    CrystalDirection,
    CrystalDirectionOverlay,
    CrystalPlane,
    CrystalPlaneOverlay,
    FrameDomain,
    Lattice,
    MillerIndex,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    build_crystal_scene,
    generate_saed_pattern,
    generate_xrd_pattern,
    get_phase_fixture,
    plot_crystal_structure_3d,
    plot_saed_pattern,
    plot_xrd_pattern,
)

warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF")
warnings.filterwarnings("ignore", message="No _symmetry_equiv_pos_as_xyz")

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_ROOT = REPO_ROOT / "docs" / "figures"

CRYSTAL_FRAME = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))

# Shared publication style. `show_title` is off so each figure carries one
# descriptive title rather than the renderer's default plus a caption.
CRYSTAL_STYLE = {
    "common": {"figure": {"figsize": [7.0, 6.4], "dpi": 150}},
    "crystal": {
        "show_title": False,
        "atom_radius_scale": 0.40,
        "bond_radius_scale": 0.15,
        "bond_alpha": 0.95,
        "bond_color": "#9aa7b6",
        "lattice_color": "#1e293b",
        "lattice_linewidth": 1.6,
        "depth_cue_strength": 0.18,
    },
}


def rocksalt_sodium_chloride() -> Phase:
    """The NaCl (halite) structure: the textbook two-species rocksalt cell."""

    lattice = Lattice(
        a=5.6402,
        b=5.6402,
        c=5.6402,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=90.0,
        crystal_frame=CRYSTAL_FRAME,
    )
    cation = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    anion = cation + np.array([0.5, 0.0, 0.0])
    sites = tuple(
        AtomicSite(label=f"Na{index + 1}", species="Na", fractional_coordinates=position)
        for index, position in enumerate(cation)
    ) + tuple(
        AtomicSite(label=f"Cl{index + 1}", species="Cl", fractional_coordinates=position % 1.0)
        for index, position in enumerate(anion)
    )
    return Phase(
        name="NaCl",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
        crystal_frame=CRYSTAL_FRAME,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group_symbol="Fm-3m",
        space_group_number=225,
        chemical_formula="NaCl",
    )


def write_crystal_demo() -> Path:
    phase = rocksalt_sodium_chloride()
    scene = build_crystal_scene(
        phase,
        repeats=(1, 1, 1),
        show_unit_cells=True,
        # Na-Cl is 2.82 A across the cube edge; the default covalent-radius sum
        # stops just short of it, so the octahedral coordination needs a nudge.
        bond_tolerance_angstrom=0.65,
        plane_overlays=(
            CrystalPlaneOverlay(
                plane=CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
                label_indices=(1, 1, 1),
                color="#f97316",
                alpha=0.18,
            ),
        ),
        direction_overlays=(
            CrystalDirectionOverlay(
                direction=CrystalDirection(np.array([1.0, 1.0, 0.0]), phase=phase),
                anchor_fractional=np.array([0.0, 0.0, 0.0]),
                label_indices=(1, 1, 0),
                color="#2563eb",
                # A full-length face diagonal ends on the far corner atom and
                # parks its label on top of that sphere.
                arrow_length_scale=0.62,
            ),
        ),
        style_overrides=CRYSTAL_STYLE,
    )
    figure = plot_crystal_structure_3d(
        scene,
        projection="persp",
        style_overrides=CRYSTAL_STYLE,
        elev_deg=20.0,
        azim_deg=32.0,
        show_legend=True,
    )
    figure.axes[0].set_title(r"NaCl (rocksalt) — $(111)$ plane and $[110]$ direction", pad=14)
    # PNG, not SVG. A lit ball-and-stick render is a shaded raster image: every
    # sphere is thousands of individually coloured mesh quads, so the SVG of this
    # one figure came to 11 MB of vector polygons that browsers choke on. Only
    # the line-art figures under docs/figures benefit from staying vector.
    path = FIGURE_ROOT / "crystal_visualization_demo.png"
    figure.savefig(path, format="png", dpi=200, bbox_inches="tight")
    plt.close(figure)
    return path


def write_xrd_demo() -> Path:
    phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=CRYSTAL_FRAME)
    pattern = generate_xrd_pattern(
        phase,
        radiation=RadiationSpec.cu_ka(),
        two_theta_range_deg=(20.0, 120.0),
        max_index=6,
        broadening_fwhm_deg=0.18,
    )
    figure = plot_xrd_pattern(pattern, theme="journal")
    path = FIGURE_ROOT / "powder_xrd_demo.svg"
    figure.savefig(path, format="svg", bbox_inches="tight")
    plt.close(figure)
    return path


def write_saed_demo() -> Path:
    phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=CRYSTAL_FRAME)
    pattern = generate_saed_pattern(
        phase,
        ZoneAxis(indices=np.array([0, 0, 1]), phase=phase),
        camera_constant_mm_angstrom=180.0,
        max_index=4,
        max_g_inv_angstrom=2.6,
    )
    # The journal theme keeps the title readable; the old export used a dark
    # background with a dark title, which made the heading invisible.
    figure = plot_saed_pattern(
        pattern,
        theme="journal",
        show_frame_indicator=True,
        frame_indicator_loc="lower left",
    )
    path = FIGURE_ROOT / "saed_demo.svg"
    figure.savefig(path, format="svg", bbox_inches="tight")
    plt.close(figure)
    return path


def main() -> None:
    for writer in (write_crystal_demo, write_xrd_demo, write_saed_demo):
        print("wrote", writer().relative_to(REPO_ROOT))


if __name__ == "__main__":
    main()
