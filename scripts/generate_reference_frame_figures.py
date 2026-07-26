"""Generate the canonical reference-frame documentation figures.

The figures written here are *generated assets*, not hand-authored ones: they
are produced by the same `pytex.plotting.frames` code path the library exposes
to users, so a documentation figure can never drift from the frame model it
illustrates. Re-run this script whenever the standard frame catalog changes.

Usage::

    python scripts/generate_reference_frame_figures.py

Outputs (tracked as canonical documentation assets):

- ``docs/figures/reference_frame_catalog.svg`` — every standard frame, one panel
  each, sharing one viewing direction so the panels are directly comparable.
- ``docs/figures/sample_frame_rd_td_nd.svg`` — the rolling-geometry sample frame
  on its own, for pages that only need RD/TD/ND.
- ``docs/figures/reference_frames.svg`` — the canonical frame chain.
- ``docs/figures/reference_frames_vectors.svg`` — the frame vocabulary and the
  crystal-to-specimen mapping.
- ``docs/figures/orientation_mapping_semantics.svg`` — the same mapping with the
  inverse drawn as a separate, explicitly-requested relationship.
- ``docs/figures/active_passive_rotation.svg`` — the two rotation languages.
- ``docs/figures/bunge_euler_geometry.svg`` — the Bunge ZXZ sequence, one panel
  per computed step.
- ``docs/figures/hcp_reference_frame.svg`` — the hexagonal crystal frame, drawn
  from a real hexagonal lattice.

The last six replaced hand-authored assets. Every one of them declared its
arrowhead markers without ``markerUnits="userSpaceOnUse"``, so SVG scaled the
head by the line's stroke width and 12-unit arrowheads rendered at 48 units,
swamping the triads they annotated. Generating them removes that whole class of
defect and keeps the drawn axes identical to the modelled axes.
"""

from __future__ import annotations

from pathlib import Path

from pytex.core.frame_catalog import (
    CARTESIAN_FRAME,
    CRYSTAL_FRAME,
    DETECTOR_FRAME,
    LABORATORY_FRAME,
    MAP_FRAME,
    SAMPLE_RD_TD_ND_FRAME,
    SPECIMEN_FRAME,
    crystal_frame,
    reciprocal_frame_for,
    specimen_frame,
)
from pytex.core.lattice import Lattice
from pytex.core.orientation import Rotation
from pytex.plotting.frame_diagrams import (
    DiagramPanel,
    active_passive_svg,
    euler_sequence_svg,
    frame_chain_svg,
    hexagonal_frame_svg,
    orientation_mapping_svg,
)
from pytex.plotting.frames import frame_catalog_svg, reference_frame_svg

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"

CATALOG_ORDER = (
    CARTESIAN_FRAME,
    SPECIMEN_FRAME,
    SAMPLE_RD_TD_ND_FRAME,
    CRYSTAL_FRAME,
    MAP_FRAME,
    DETECTOR_FRAME,
    LABORATORY_FRAME,
)


def write_catalog_figure() -> Path:
    """Write the multi-panel standard-frame catalog figure."""

    svg = frame_catalog_svg(
        CATALOG_ORDER,
        columns=4,
        title="PyTex Standard Reference Frames",
        subtitle=(
            "One panel per catalog frame. Axis vectors are components in the canonical "
            "Cartesian reference; all panels share one viewing direction."
        ),
    )
    destination = FIGURES_DIR / "reference_frame_catalog.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_sample_frame_figure() -> Path:
    """Write the standalone rolling-geometry sample-frame figure."""

    # Width is set by the subtitle, which is the widest element in the figure.
    svg = reference_frame_svg(
        SAMPLE_RD_TD_ND_FRAME,
        width=560.0,
        height=400.0,
        title="Sample Frame (RD, TD, ND)",
        subtitle=(
            "RD rolling direction; TD transverse direction; ND sheet normal - a right-handed triad"
        ),
    )
    destination = FIGURES_DIR / "sample_frame_rd_td_nd.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_frame_chain_figure() -> Path:
    """Write the canonical frame-chain figure."""

    # The linear chain is exactly the one fixed in
    # docs/standards/notation_and_conventions.md. The reciprocal frame is NOT a
    # link in it — it is dual to the crystal frame — so it hangs below the
    # crystal panel rather than extending the row.
    panels = (
        DiagramPanel(frame=CRYSTAL_FRAME, caption="lattice-fixed basis"),
        DiagramPanel(frame=SPECIMEN_FRAME, caption="macroscopic sample frame"),
        DiagramPanel(frame=MAP_FRAME, caption="scan sampling grid"),
        DiagramPanel(frame=DETECTOR_FRAME, caption="pattern-plane coordinates"),
        DiagramPanel(frame=LABORATORY_FRAME, caption="instrument geometry"),
    )
    relationships = (
        "orientation g",
        "registration",
        "detector geometry",
        "instrument pose",
    )
    svg = frame_chain_svg(
        panels,
        relationships,
        dual_panel=DiagramPanel(
            frame=reciprocal_frame_for(CRYSTAL_FRAME),
            caption="dual basis, starred axes",
        ),
        dual_label="duality",
        dual_of=0,
        dual_note=(
            "Crystal geometry is also expressible in the reciprocal basis a*, b*, c*. "
            "The star marks the basis, never the Miller indices: (hkl) are already "
            "reciprocal-basis components."
        ),
    )
    destination = FIGURES_DIR / "reference_frames.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_orientation_mapping_figure() -> Path:
    """Write the orientation-mapping semantics figure."""

    # A concrete, reproducible orientation so the middle panel shows real
    # arithmetic rather than a decorative rotation.
    rotation = Rotation.from_bunge_euler(35.0, 28.0, 15.0)
    svg = orientation_mapping_svg(
        crystal_frame(),
        specimen_frame(),
        rotation_matrix=rotation.as_matrix(),
    )
    destination = FIGURES_DIR / "reference_frames_vectors.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_orientation_semantics_figure() -> Path:
    """Write the orientation-direction figure, with the inverse shown as separate."""

    rotation = Rotation.from_bunge_euler(35.0, 28.0, 15.0)
    svg = orientation_mapping_svg(
        crystal_frame(),
        specimen_frame(),
        rotation_matrix=rotation.as_matrix(),
        show_inverse=True,
        title="Orientation Mapping Semantics",
        subtitle=(
            "PyTex fixes orientation meaning as crystal frame to specimen frame. Reversing it is "
            "a different object and must be requested explicitly, never assumed."
        ),
    )
    destination = FIGURES_DIR / "orientation_mapping_semantics.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_active_passive_figure() -> Path:
    """Write the active-versus-passive rotation-language figure."""

    rotation = Rotation.from_bunge_euler(35.0, 28.0, 15.0)
    svg = active_passive_svg(
        crystal_frame(),
        specimen_frame(),
        rotation_matrix=rotation.as_matrix(),
    )
    destination = FIGURES_DIR / "active_passive_rotation.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_euler_sequence_figure() -> Path:
    """Write the Bunge ZXZ sequence, computed step by step."""

    svg = euler_sequence_svg(specimen_frame(), phi1_deg=35.0, Phi_deg=45.0, phi2_deg=30.0)
    destination = FIGURES_DIR / "bunge_euler_geometry.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def write_hexagonal_frame_figure() -> Path:
    """Write the canonical hexagonal-frame figure from a real hcp lattice."""

    # Zirconium alpha-hcp lattice parameters, matching the zr_hcp phase fixture.
    lattice = Lattice(
        a=3.2320,
        b=3.2320,
        c=5.1470,
        alpha_deg=90.0,
        beta_deg=90.0,
        gamma_deg=120.0,
        crystal_frame=crystal_frame(),
    )
    svg = hexagonal_frame_svg(lattice)
    destination = FIGURES_DIR / "hcp_reference_frame.svg"
    destination.write_text(svg, encoding="utf-8")
    return destination


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    written = (
        write_catalog_figure(),
        write_sample_frame_figure(),
        write_frame_chain_figure(),
        write_orientation_mapping_figure(),
        write_orientation_semantics_figure(),
        write_active_passive_figure(),
        write_euler_sequence_figure(),
        write_hexagonal_frame_figure(),
    )
    for destination in written:
        print(f"wrote {destination.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
