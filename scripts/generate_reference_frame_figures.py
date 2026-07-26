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


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for destination in (write_catalog_figure(), write_sample_frame_figure()):
        print(f"wrote {destination.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
