from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pytex.plotting._visual_regression import (
    build_visual_regression_figures,
    save_visual_regression_svg,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
VISUAL_REGRESSION_ROOT = REPO_ROOT / "fixtures" / "visual_regression"


def test_pinned_visual_regression_svgs_match_current_runtime_output(tmp_path: Path) -> None:
    for name, figure in build_visual_regression_figures().items():
        generated_path = save_visual_regression_svg(figure, tmp_path / f"{name}.svg")
        pinned_path = VISUAL_REGRESSION_ROOT / f"{name}.svg"
        assert pinned_path.exists()
        assert generated_path.read_bytes() == pinned_path.read_bytes()
        plt.close(figure)
