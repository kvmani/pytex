from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pytex.plotting._visual_regression import (
    build_visual_regression_figures,
    save_visual_regression_svg,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> int:
    output_root = _repo_root() / "fixtures" / "visual_regression"
    for name, figure in build_visual_regression_figures().items():
        save_visual_regression_svg(figure, output_root / f"{name}.svg")
        plt.close(figure)
    print(f"Wrote visual regression artifacts to {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
