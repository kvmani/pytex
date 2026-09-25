"""The figures of FIB lamella planning (section 9 of the foundation document).

Two kinds of check. The canonical frame-convention SVG is a generated asset: it
must be byte-identical to its generator's output, accessible, and drawn from the
library's own geometry. The six runtime figures must render in both the light
and the dark theme without a warning and without leaving a figure open.
"""

from __future__ import annotations

import importlib.util
import io
import math
import sys
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest

from pytex.fib import default_tem_stage, rank_grains
from pytex.plotting.fib_figures import FIB_FIGURE_KINDS, FIB_THEMES, fib_figure, reach_curve_deg
from pytex.tem.stage import RectangularEnvelope

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "generate_fib_lamella_figures.py"


def _generator():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("generate_fib_lamella_figures", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestCanonicalFigure:
    def test_the_committed_figure_is_the_generators_output(self) -> None:
        generator = _generator()
        committed = generator.FRAMES_FIGURE.read_text(encoding="utf-8").replace("\r\n", "\n")
        assert committed == generator.frames_svg()

    def test_the_figure_is_well_formed_accessible_and_in_arial(self) -> None:
        root = ElementTree.fromstring(_generator().frames_svg())
        namespace = "{http://www.w3.org/2000/svg}"
        assert root.find(f"{namespace}title") is not None
        assert root.find(f"{namespace}desc") is not None
        for marker in root.iter(f"{namespace}marker"):
            assert marker.get("markerUnits") == "userSpaceOnUse"
        for text in root.iter(f"{namespace}text"):
            assert text.get("font-family") == "Arial"

    def test_the_drawn_lamella_is_the_library_solution(self) -> None:
        generator = _generator()
        geometry = generator.example_geometry()
        assert geometry.eps_deg == pytest.approx(generator.EPS_DEG, abs=1e-9)
        # The target's in-plane part is at AZIMUTH_DEG, so the lamella normal is
        # along it (its canonical azimuth in [0, 180)).
        assert geometry.theta_sample_deg == pytest.approx(generator.AZIMUTH_DEG % 180.0, abs=1e-9)
        assert f"θ = {generator.AZIMUTH_DEG:.0f}°" in generator.frames_svg()


@pytest.fixture(scope="module")
def planned():
    from pytex.app.ebsd_gallery import build_map

    crystal_map = build_map("equiaxed_polycrystal", grid=48, step_um=2.0)
    segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
    report = rank_grains(segmentation, (1, 1, 1), shortlist=3)
    rows, cols = crystal_map.grid_shape
    image_quality = np.asarray(crystal_map.get_property("image_quality")).reshape(rows, cols)
    return report, segmentation, image_quality


def _arguments(kind: str, planned) -> dict:  # type: ignore[no-untyped-def]
    report, segmentation, image_quality = planned
    plan = report.best()
    return {
        "plan_view": {
            "plan": plan,
            "background": image_quality,
            "extent": report.raster.extent_um,
            "grain_mask": segmentation.label_grid == plan.grain_id,
        },
        "section": {"plan": plan},
        "stereogram": {"plan": plan, "guaranteed_radius_deg": report.guaranteed_radius_deg},
        "phi_feasibility": {"plan": plan, "envelope": default_tem_stage().envelope},
        "preparability": {
            "raster": report.raster,
            "plans": report.plans,
            "guaranteed_radius_deg": report.guaranteed_radius_deg,
        },
        "saed": {"plan": plan},
    }[kind]


class TestRuntimeFigures:
    @pytest.mark.parametrize("theme", sorted(FIB_THEMES))
    @pytest.mark.parametrize("kind", FIB_FIGURE_KINDS)
    def test_every_figure_renders_in_both_themes(self, planned, kind, theme) -> None:
        import matplotlib.pyplot as plt

        before = set(plt.get_fignums())
        figure = fib_figure(kind, _arguments(kind, planned), theme=theme)
        buffer = io.StringIO()
        figure.savefig(buffer, format="svg")
        svg = buffer.getvalue()
        assert "<svg" in svg and len(svg) > 2000
        # The theme's paper colour is the figure background.
        assert FIB_THEMES[theme]["paper"].lower() in svg.lower()
        # Nothing was registered with pyplot, so nothing can leak.
        assert set(plt.get_fignums()) == before

    def test_an_unknown_theme_or_kind_is_refused(self, planned) -> None:
        with pytest.raises(ValueError, match="theme"):
            fib_figure("section", {"plan": planned[0].best()}, theme="sepia")
        with pytest.raises(ValueError, match="kind"):
            fib_figure("hologram", {})

    def test_the_plan_view_needs_a_placement(self, planned) -> None:
        from dataclasses import replace

        bare = replace(planned[0].best(), placement=None)
        with pytest.raises(ValueError, match="placement"):
            fib_figure("plan_view", {"plan": bare})


class TestReachCurve:
    def test_a_symmetric_rectangle_reaches_its_half_range_on_the_axes(self) -> None:
        envelope = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)
        reach = reach_curve_deg(envelope, np.array([0.0, 90.0, 180.0, 270.0]))
        np.testing.assert_allclose(reach, 30.0, atol=0.1)

    def test_the_diagonal_reaches_further_than_the_axes(self) -> None:
        # Both tilts share the work along a diagonal of the rectangle.
        envelope = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)
        reach = reach_curve_deg(envelope, np.array([0.0, 45.0]))
        assert reach[1] > reach[0] + 5.0
        # The exact corner: alpha = beta = 30 is a tilt of acos(cos^2 30) from zero.
        corner = math.degrees(math.acos(math.cos(math.radians(30.0)) ** 2))
        assert reach[1] <= corner + 0.2
