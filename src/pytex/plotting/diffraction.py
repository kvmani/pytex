from __future__ import annotations

from typing import Any

import numpy as np

from pytex.diffraction.saed import SAEDPattern
from pytex.diffraction.xrd import PowderPattern
from pytex.plotting.styles import resolve_style


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, object


def plot_xrd_pattern(
    pattern: PowderPattern,
    *,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    xrd_style = style["xrd"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=common["figure"]["facecolor"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(common["figure"]["axes_facecolor"])
    axes.plot(
        pattern.two_theta_grid_deg,
        pattern.intensity_grid,
        color=xrd_style["line_color"],
        linewidth=float(common["line"]["width"]),
    )
    axes.fill_between(
        pattern.two_theta_grid_deg,
        pattern.intensity_grid,
        color=xrd_style["fill_color"],
        alpha=0.55,
    )
    axes.set_xlabel(r"2$\theta$ (deg)")
    axes.set_ylabel("normalized intensity")
    axes.set_title(f"{pattern.phase.name} Powder XRD ({pattern.radiation.name})")
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    if xrd_style.get("annotate_peaks", True):
        ranked = sorted(
            pattern.reflections,
            key=lambda reflection: reflection.intensity,
            reverse=True,
        )
        for reflection in ranked[: int(xrd_style.get("max_labels", 12))]:
            label = "(" + " ".join(str(int(value)) for value in reflection.miller_indices) + ")"
            axes.axvline(
                reflection.two_theta_deg,
                color=xrd_style["peak_color"],
                alpha=0.25,
                linewidth=0.8,
            )
            axes.text(
                reflection.two_theta_deg,
                1.02,
                label,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=float(common["font"]["size"]) - 1.0,
                color=xrd_style["peak_color"],
                transform=axes.get_xaxis_transform(),
            )
    fig.tight_layout()
    return fig


def plot_saed_pattern(
    pattern: SAEDPattern,
    *,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    ax: Any | None = None,
) -> Any:
    plt, _ = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    saed_style = style["saed"]
    if ax is None:
        fig, axes = plt.subplots(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=saed_style["background"],
        )
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(saed_style["background"])
    if pattern.spots:
        coordinates = np.vstack([spot.detector_coordinates for spot in pattern.spots])
        intensities = np.array([spot.intensity for spot in pattern.spots], dtype=np.float64)
        if np.max(intensities) > 0.0:
            sizes = float(saed_style["spot_scale"]) * intensities / np.max(intensities)
        else:
            sizes = np.full_like(intensities, float(saed_style["spot_scale"]))
        axes.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=sizes,
            color=saed_style["spot_color"],
            edgecolors="white",
            linewidths=0.35,
            alpha=0.95,
        )
        if saed_style.get("annotate_spots", True):
            for spot in pattern.spots[: int(saed_style.get("max_labels", 16))]:
                if spot.label is None:
                    continue
                axes.text(
                    float(spot.detector_coordinates[0]),
                    float(spot.detector_coordinates[1]),
                    spot.label,
                    fontsize=float(common["font"]["size"]) - 1.0,
                    color=saed_style["label_color"],
                    ha="left",
                    va="bottom",
                )
    extent = pattern.detector_extent_mm()
    axes.axhline(0.0, color=saed_style["ring_color"], linewidth=0.8, alpha=0.45)
    axes.axvline(0.0, color=saed_style["ring_color"], linewidth=0.8, alpha=0.45)
    axes.set_xlim(-extent, extent)
    axes.set_ylim(-extent, extent)
    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("detector u (mm)")
    axes.set_ylabel("detector v (mm)")
    axes.set_title(
        "SAED Pattern "
        + "("
        + " ".join(str(int(value)) for value in pattern.zone_axis.indices)
        + " zone axis)"
    )
    axes.grid(alpha=float(common["figure"]["grid_alpha"]))
    fig.tight_layout()
    return fig
