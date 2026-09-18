"""Figures of how a Kearns parameter is reached, for every Kearns route.

Purpose
-------
A Kearns parameter is a single number, ``f = ⟨cos²φ⟩`` of the basal poles about
a specimen direction, and it hides everything a reader needs to judge it: which
tilts were measured, how the pole density varies with tilt, and which tilts
actually carry the result once the ``sin φ`` volume weighting is applied. These
figures show those steps. The running Kearns sum is an exact rearrangement of
the quadrature the routes use — ``f = Σ I sinφ cos²φ / Σ I sinφ`` — so its last
value *is* the reported ``f``; ``tests/unit/test_app_kearns_figures.py`` holds it
to that.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.figures import COLORS, draw_ticks, render_figure
from pytex.app.results import ResultFigure

__all__ = [
    "orientation_statistics",
    "orientation_tilt_figure",
    "running_kearns_sum",
    "section_scan_figure",
    "tilt_profile_figure",
    "triad_figure",
]

_RANDOM = 1.0 / 3.0


def running_kearns_sum(
    polar_deg: Sequence[float] | np.ndarray, density: Sequence[float] | np.ndarray
) -> np.ndarray:
    """``Σ_{j≤k} I_j sinφ_j cos²φ_j / Σ_all I sinφ`` at each tilt node, in ascending tilt.

    The last value is Kearns' ``f`` for the profile, exactly as
    :func:`pytex.texture.kearns.kearns_from_tilt_profile` computes it.
    """

    tilts = np.asarray(polar_deg, dtype=float)
    values = np.asarray(density, dtype=float)
    order = np.argsort(tilts)
    radians = np.deg2rad(tilts[order])
    volume = values[order] * np.sin(radians)
    total = float(np.sum(volume))
    if total <= 0.0:
        return np.full(tilts.size, np.nan)
    cumulative: np.ndarray = np.cumsum(volume * np.cos(radians) ** 2) / total
    return cumulative


def tilt_profile_figure(
    *,
    key: str,
    polar_deg: Sequence[float] | np.ndarray,
    density: Sequence[float] | np.ndarray,
    direction: str,
    f_value: float,
    points: Sequence[Mapping[str, Any]] = (),
    title: str | None = None,
) -> ResultFigure:
    """The basal-pole tilt profile, and how the Kearns sum builds up across it."""

    tilts = np.asarray(polar_deg, dtype=float)
    values = np.asarray(density, dtype=float)
    order = np.argsort(tilts)
    tilts = tilts[order]
    values = values[order]
    running = running_kearns_sum(tilts, values)
    radians = np.deg2rad(tilts)
    weight = values * np.sin(radians)
    weight_share = weight / float(np.sum(weight)) if np.sum(weight) > 0 else weight
    step = float(np.median(np.diff(tilts))) if tilts.size > 1 else 5.0

    def draw(figure: Any) -> None:
        top, middle, bottom = figure.subplots(3, 1, sharex=True, height_ratios=[1.6, 1.0, 1.2])
        top.bar(
            tilts,
            values,
            width=0.85 * step,
            color=COLORS["band_2"],
            edgecolor=COLORS["difference"],
            lw=0.6,
            label="Tilt profile I(φ)",
        )
        if points:
            x = [float(point["tilt"]) for point in points]
            y = [float(point["density"]) for point in points]
            top.plot(x, y, "o", ms=4.5, color=COLORS["data"], label="Measured reflection")
            # Only reflections that carry weight are named; the rest would be a
            # pile of overlapping labels along the zero line at high tilt.
            strongest = max(y) if y else 0.0
            for point in points:
                if point.get("label") and float(point["density"]) >= 0.1 * strongest:
                    top.annotate(
                        str(point["label"]),
                        (float(point["tilt"]), float(point["density"])),
                        xytext=(4, 3),
                        textcoords="offset points",
                        fontsize=6.5,
                        color=COLORS["background"],
                    )
        top.axhline(1.0, color=COLORS["guide"], lw=0.7, ls="--")
        top.set_ylabel("Pole density")
        top.legend(loc="best", frameon=False, fontsize=7)
        middle.bar(
            tilts,
            weight_share,
            width=0.85 * step,
            color=COLORS["band_3"],
            edgecolor=COLORS["model"],
            lw=0.6,
        )
        middle.set_ylabel("Volume share\nI sinφ / ΣI sinφ")
        bottom.plot(
            tilts,
            running,
            "o-",
            ms=3.5,
            lw=1.2,
            color=COLORS["model"],
            label=f"Running sum → f_{direction} = {f_value:.4f}",
        )
        bottom.axhline(_RANDOM, color=COLORS["guide"], lw=0.8, ls="--", label="Random, 1/3")
        bottom.axhline(f_value, color=COLORS["accent"], lw=0.8, ls=":")
        bottom.set_ylabel("Kearns sum")
        bottom.set_xlabel("Basal tilt φ from the direction (°)")
        bottom.set_xlim(0.0, 90.0)
        bottom.legend(loc="best", frameon=False, fontsize=7)

    dominant = float(tilts[int(np.argmax(weight_share))]) if tilts.size else float("nan")
    return render_figure(
        draw,
        key=key,
        title=title or f"How f_{direction} is built from the tilt profile",
        caption=(
            "Top: the azimuthally averaged basal-pole density against tilt from the direction "
            "(bars; 1 is random), with the measured reflections it was interpolated from. "
            "Middle: each tilt's share of the volume, I·sinφ normalized. Bottom: the running "
            "Kearns sum Σ I sinφ cos²φ / Σ I sinφ, whose last value is the reported f."
        ),
        interpretation=(
            f"The tilts near {dominant:.0f}° carry the most volume. Because of the sinφ factor, "
            "an intense pole near φ = 0 contributes little volume, and the high-tilt bins decide "
            "f even where the density is low — so a profile that stops short of 90°, or "
            "interpolates across a gap, biases f most there."
        ),
        height_in=4.6,
    )


def triad_figure(
    *,
    key: str,
    labels: Sequence[str],
    values: Sequence[float],
    closure_by_construction: bool,
    uncertainties: Sequence[float] | None = None,
    truth: Mapping[str, float] | None = None,
) -> ResultFigure:
    """The Kearns parameter along each direction against the random value 1/3."""

    heights = np.asarray(values, dtype=float)
    total = float(np.sum(heights))

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(labels))
        axes.bar(
            positions,
            heights,
            width=0.55,
            color=COLORS["band_2"],
            edgecolor=COLORS["difference"],
            lw=0.8,
            yerr=None if uncertainties is None else np.asarray(uncertainties, dtype=float),
            capsize=4,
            label="f",
        )
        for x, value in zip(positions, heights, strict=True):
            axes.text(x, value + 0.02, f"{value:.4f}", ha="center", fontsize=8)
        if truth:
            known = [truth.get(label) for label in labels]
            axes.plot(
                [x for x, value in zip(positions, known, strict=True) if value is not None],
                [value for value in known if value is not None],
                "D",
                ms=6,
                color=COLORS["accent"],
                label="Known answer",
            )
        axes.axhline(_RANDOM, color=COLORS["guide"], lw=0.9, ls="--", label="Random, 1/3")
        axes.set_xticks(positions, labels=[f"f_{label}" for label in labels])
        axes.set_ylim(0.0, max(1.0, float(np.max(heights)) + 0.12))
        axes.set_ylabel("Kearns parameter")
        axes.legend(loc="upper right", frameon=False, fontsize=7)

    closure = (
        f"The three values sum to {total:.4f}. "
        + (
            "That is closure by construction — they are the diagonal of one orientation tensor "
            "whose trace is 1 — so the sum checks the arithmetic, not the measurement."
            if closure_by_construction
            else "They were measured independently, so the departure from 1 measures the "
            "systematic error of the measurement (Kearns found 0.94 to 1.06)."
        )
        if len(labels) == 3
        else ""
    )
    return render_figure(
        draw,
        key=key,
        title="Kearns parameters against the random value",
        caption=(
            "The Kearns parameter along each specimen direction (bars"
            + (", ±1 standard error" if uncertainties is not None else "")
            + "), with the untextured value 1/3 dashed"
            + (" and the known answer of the demonstration (diamonds)" if truth else "")
            + "."
        ),
        interpretation=(
            "A bar above 1/3 means basal poles concentrate along that direction; below, they "
            "avoid it. " + closure
        ).strip(),
        width_in=4.8,
        height_in=3.0,
    )


def orientation_statistics(
    orientations: Any, pole: Any, directions: np.ndarray, *, include_symmetry_family: bool
) -> dict[str, np.ndarray]:
    """Per-grain ``cos²φ`` along each direction, and the mean and its sampling standard error.

    Each orientation's family members are averaged first, because the members of
    one grain are not independent samples; the standard error is then
    ``std(cos²)/√N`` over grains. The means are the route's ``f`` values.
    """

    normals = (
        pole.phase.symmetry.equivalent_vectors(pole.normal)
        if include_symmetry_family
        else pole.normal[None, :]
    )
    axes = np.asarray(directions, dtype=float)
    per_grain = np.zeros((len(orientations), axes.shape[0]))
    for normal in normals:
        mapped = orientations.map_crystal_directions(normal)
        vectors = np.asarray(getattr(mapped, "values", mapped), dtype=float)
        vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
        per_grain += np.square(vectors @ axes.T)
    per_grain /= len(normals)
    count = per_grain.shape[0]
    return {
        "cos2": per_grain,
        "mean": per_grain.mean(axis=0),
        "standard_error": (
            per_grain.std(axis=0, ddof=1) / np.sqrt(count)
            if count > 1
            else np.full(axes.shape[0], np.nan)
        ),
    }


def orientation_tilt_figure(
    *, cos2: np.ndarray, labels: Sequence[str], values: Sequence[float]
) -> ResultFigure:
    """How the basal poles of the orientation set are tilted from each specimen axis."""

    tilt = np.rad2deg(np.arccos(np.sqrt(np.clip(cos2, 0.0, 1.0))))
    edges = np.linspace(0.0, 90.0, 19)
    centres = 0.5 * (edges[1:] + edges[:-1])
    random = np.cos(np.deg2rad(edges[:-1])) - np.cos(np.deg2rad(edges[1:]))
    palette = (COLORS["difference"], COLORS["model"], COLORS["accent"])

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        for column, label in enumerate(labels):
            counts, _ = np.histogram(tilt[:, column], bins=edges)
            share = counts / max(counts.sum(), 1)
            axes.step(
                centres,
                share,
                where="mid",
                lw=1.4,
                color=palette[column % 3],
                label=f"from {label} (f = {values[column]:.3f})",
            )
        axes.plot(centres, random, "--", lw=0.9, color=COLORS["guide"], label="Random (∝ sin φ)")
        axes.set_xlabel("Basal-pole tilt φ from the axis (°)")
        axes.set_ylabel("Fraction of grains")
        axes.set_xlim(0.0, 90.0)
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="basal_tilt_distribution",
        title="Where the basal poles point",
        caption=(
            "For each specimen axis, the fraction of orientations whose basal pole is tilted by "
            "φ from it (5° bins), against the random distribution, which is proportional to "
            "sin φ. f is the mean of cos²φ over these grains."
        ),
        interpretation=(
            "An axis whose distribution crowds towards φ = 0 collects the basal poles and has "
            "a large f; one crowding towards 90° has a small f. The difference from the dashed "
            "curve is the texture."
        ),
        height_in=3.0,
    )


def section_scan_figure(
    *, key: str, name: str, pattern: Mapping[str, Any], reflections: Sequence[Mapping[str, Any]]
) -> ResultFigure:
    """One section's measured scan, its background, and every predicted reflection's fate."""

    axis = np.asarray(pattern["two_theta_deg"], dtype=float)
    intensity = np.asarray(pattern["intensity"], dtype=float)
    background = np.asarray(pattern["background"], dtype=float)
    used = [row for row in reflections if row.get("used")]
    unused = [row for row in reflections if not row.get("used")]

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        axes.plot(
            axis,
            np.sqrt(np.clip(intensity, 0.0, None)),
            lw=0.6,
            color=COLORS["data"],
            label="Measured (√counts)",
        )
        axes.plot(
            axis,
            np.sqrt(np.clip(background, 0.0, None)),
            lw=1.0,
            ls="--",
            color=COLORS["background"],
            label="Background",
        )
        draw_ticks(
            axes,
            [float(row["two_theta_expected_deg"]) for row in used],
            labels=[str(row["plane"]) for row in used],
            color=COLORS["accent"],
            label="Reflection used",
        )
        if unused:
            draw_ticks(
                axes,
                [float(row["two_theta_expected_deg"]) for row in unused],
                color=COLORS["warning"],
                label="Reflection not used",
            )
        axes.set_xlabel("2θ (°)")
        axes.set_ylabel("√intensity")
        axes.legend(loc="upper right", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key=key,
        title=f"{name.capitalize()} section: the scan and its reflections",
        caption=(
            "The measured scan and its background on a square-root scale (so weak reflections "
            "are visible), with every predicted reflection: green were used, red were excluded "
            "(too weak in the random standard, or no peak within the matching tolerance)."
        ),
        interpretation=(
            "Every green tick should sit under a peak. A peak with no tick is a reflection of "
            "something else; a red tick under a clear peak was excluded by the matching rules "
            "and deserves a look at the tolerance."
        ),
        height_in=2.8,
    )
