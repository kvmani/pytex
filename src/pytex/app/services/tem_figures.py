# ruff: noqa: RUF001, RUF002
"""Figures of how a TEM pattern was measured, fitted and indexed.

Purpose
-------
The TEM panels draw the pattern itself, with the picks and the fitted lattice
over it. What the pattern cannot show at its own scale is how well each pick
agrees with the model: a residual of half a pixel is invisible on a 2000-pixel
pattern. These figures put the misfit at a readable scale — residual vectors
magnified, observed spacings against calculated ones, the candidates' scores
broken down — and, for a simulated micrograph, the rotational average of its
power spectrum against the transfer function of the lens that made it. Every
number is read from the operation's own result.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.figures import COLORS, render_figure
from pytex.app.fitstats import fit_line
from pytex.app.results import ResultFigure

__all__ = [
    "hrem_spectrum_figure",
    "lattice_fit_figure",
    "radial_average",
    "solve_figures",
    "thickness_figure",
]


def lattice_fit_figure(
    rows: Sequence[Mapping[str, Any]], *, centre: Sequence[float], rms: float
) -> ResultFigure:
    """Every pick with its lattice node, residuals magnified, and residual per spot."""

    x = np.array([float(row["x"]) for row in rows])
    y = np.array([float(row["y"]) for row in rows])
    px = np.array([float(row["predicted_x"]) for row in rows])
    py = np.array([float(row["predicted_y"]) for row in rows])
    residual = np.array([float(row["residual"]) for row in rows])
    inlier = np.array([row["verdict"] == "on the lattice" for row in rows])
    span = float(max(np.ptp(np.r_[x, px, centre[0]]), np.ptp(np.r_[y, py, centre[1]]), 1.0))
    largest = float(np.max(residual)) if residual.size else 0.0
    magnification = 0.08 * span / largest if largest > 0 else 1.0

    def draw(figure: Any) -> None:
        left, right = figure.subplots(1, 2, width_ratios=[1.3, 1.0])
        left.plot(px, py, "+", ms=9, mew=1.2, color=COLORS["model"], label="Fitted lattice node")
        left.plot(
            x[inlier], y[inlier], "o", ms=4, color=COLORS["data"], label="Pick on the lattice"
        )
        if np.any(~inlier):
            left.plot(
                x[~inlier],
                y[~inlier],
                "x",
                ms=6,
                color=COLORS["warning"],
                label="Pick off the lattice",
            )
        left.quiver(
            px,
            py,
            (x - px) * magnification,
            (y - py) * magnification,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            width=0.004,
            color=COLORS["difference"],
        )
        left.plot(
            [centre[0]], [centre[1]], "s", ms=6, color=COLORS["accent"], label="Transmitted beam"
        )
        for index, row in enumerate(rows):
            left.annotate(
                str(row["spot"]),
                (x[index], y[index]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=6.5,
                color=COLORS["background"],
            )
        left.set_aspect("equal")
        left.invert_yaxis()
        left.set_xlabel("x (picking units)")
        left.set_ylabel("y (picking units)")
        left.legend(loc="best", frameon=False, fontsize=6.5)
        positions = np.arange(len(rows))
        right.bar(
            positions,
            residual,
            color=[COLORS["data"] if flag else COLORS["warning"] for flag in inlier],
        )
        right.axhline(rms, color=COLORS["model"], lw=1.0, ls="--", label=f"r.m.s. {rms:.2f}")
        right.set_xticks(positions, labels=[str(row["spot"]) for row in rows], fontsize=7)
        right.set_xlabel("Spot")
        right.set_ylabel("Distance to node (picking units)")
        right.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="lattice_fit_residuals",
        title="Each pick against its lattice node",
        caption=(
            "Left: picks (dots; red crosses are off the lattice) and the nodes they were "
            f"assigned to (+), with the pick-minus-node residual drawn ×{magnification:.0f} "
            "as blue arrows so it is visible at the scale of the pattern. Right: the length of "
            "each residual, with the r.m.s. of the inliers dashed."
        ),
        interpretation=(
            "Random arrows of similar length are picking noise. Arrows that all point the same "
            "way are a mis-set transmitted beam; arrows that grow with distance from it are a "
            "wrong lattice spacing or a distorted (elliptical) pattern — the camera calibration "
            "or projector distortion, not the crystal."
        ),
        height_in=3.4,
    )


def solve_figures(
    rows: Sequence[Mapping[str, Any]],
    alternatives: Sequence[Mapping[str, Any]],
    *,
    best_label: str,
) -> tuple[ResultFigure, ...]:
    """Observed against calculated spacings for the best solution; every candidate's score."""

    indexed = [row for row in rows if row.get("d_calculated")]
    figures: list[ResultFigure] = []
    if indexed:
        observed = np.array([float(row["d_observed"]) for row in indexed])
        calculated = np.array([float(row["d_calculated"]) for row in indexed])

        def draw_parity(figure: Any) -> None:
            left, right = figure.subplots(1, 2)
            top = float(max(observed.max(), calculated.max())) * 1.08
            left.plot([0.0, top], [0.0, top], "--", lw=0.8, color=COLORS["guide"])
            left.plot(calculated, observed, "o", ms=5, color=COLORS["data"])
            for row, cx, ox in zip(indexed, calculated, observed, strict=True):
                left.annotate(
                    str(row["hkl"]),
                    (cx, ox),
                    xytext=(4, -8),
                    textcoords="offset points",
                    fontsize=6.5,
                    color=COLORS["background"],
                )
            left.set_xlim(0.0, top)
            left.set_ylim(0.0, top)
            left.set_aspect("equal")
            left.set_xlabel("d calculated (Å)")
            left.set_ylabel("d measured (Å)")
            deviation = 100.0 * (observed - calculated) / calculated
            right.axhline(0.0, color=COLORS["guide"], lw=0.8)
            right.bar(np.arange(len(indexed)), deviation, color=COLORS["difference"])
            right.set_xticks(
                np.arange(len(indexed)), labels=[str(row["spot"]) for row in indexed], fontsize=7
            )
            right.set_xlabel("Spot")
            right.set_ylabel("(d_meas − d_calc)/d_calc (%)")

        figures.append(
            render_figure(
                draw_parity,
                key="solve_spacings",
                title=f"Measured against calculated spacings ({best_label})",
                caption=(
                    "Left: each indexed spot's measured interplanar spacing against the "
                    "spacing of the reflection it was indexed as, with the 1:1 line. Right: "
                    "the relative deviation of each spot."
                ),
                interpretation=(
                    "Deviations of the same sign on every spot are a camera-constant error "
                    "and scale all spacings together; deviations that change with direction "
                    "are an elliptical pattern. Scatter of a percent or two is ordinary for "
                    "a SAED pattern calibrated against its own camera constant."
                ),
                height_in=3.0,
            )
        )
    if alternatives:

        def draw_scores(figure: Any) -> None:
            axes = figure.subplots()
            names = [f"{item['phase']} {item.get('zone_axis', '')}" for item in alternatives]
            criteria = (
                ("length_agreement", "Spacings", COLORS["difference"]),
                ("angle_agreement", "Angles", COLORS["model"]),
                ("coverage_agreement", "Coverage", COLORS["accent"]),
            )
            positions = np.arange(len(names))
            width = 0.8 / len(criteria)
            for index, (key, label, color) in enumerate(criteria):
                axes.bar(
                    positions + (index - 1) * width,
                    [float(item.get(key) or 0.0) for item in alternatives],
                    width,
                    color=color,
                    label=label,
                )
            axes.plot(
                positions,
                [float(item.get("score") or 0.0) for item in alternatives],
                "D",
                ms=6,
                color=COLORS["data"],
                label="Score",
            )
            axes.set_xticks(positions, labels=names, rotation=25, ha="right", fontsize=7)
            axes.set_ylim(0.0, 1.05)
            axes.set_ylabel("Agreement (0 to 1)")
            axes.legend(loc="upper right", frameon=False, fontsize=7, ncols=2)

        figures.append(
            render_figure(
                draw_scores,
                key="solve_candidates",
                title="Every candidate solution, scored",
                caption=(
                    "The agreement of each candidate's spacings, angles and spot coverage "
                    "with the picks (bars), and the weighted score (diamonds), in rank order."
                ),
                interpretation=(
                    "A clear winner leads on all three criteria. Two candidates with similar "
                    "scores are two indexings the pattern cannot tell apart — tilt to a second "
                    "zone axis before believing either."
                ),
                height_in=3.2,
            )
        )
    return tuple(figures)


def radial_average(
    image: np.ndarray, *, pixel_size_angstrom: float, bins: int = 120
) -> tuple[np.ndarray, np.ndarray]:
    """Rotational average of a centred 2-D spectrum against spatial frequency (Å⁻¹).

    The spectrum is assumed centred with zero frequency at the middle pixel, as
    the HRTEM result stores it; frequency is ``|k|/(N·dx)`` per axis.
    """

    values = np.asarray(image, dtype=float)
    rows, columns = values.shape
    qy = (np.arange(rows) - rows // 2) / (rows * pixel_size_angstrom)
    qx = (np.arange(columns) - columns // 2) / (columns * pixel_size_angstrom)
    radius = np.hypot(*np.meshgrid(qx, qy))
    limit = float(min(np.max(np.abs(qx)), np.max(np.abs(qy))))
    edges = np.linspace(0.0, limit, bins + 1)
    index = np.digitize(radius.ravel(), edges) - 1
    valid = (index >= 0) & (index < bins)
    totals = np.bincount(index[valid], weights=values.ravel()[valid], minlength=bins)
    counts = np.bincount(index[valid], minlength=bins)
    centres = 0.5 * (edges[1:] + edges[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        average = totals / counts
    keep = counts > 0
    return centres[keep], average[keep]


def hrem_spectrum_figure(
    *,
    power_spectrum: np.ndarray,
    pixel_size_angstrom: float,
    frequencies: np.ndarray,
    transfer: np.ndarray,
    envelope: np.ndarray,
    point_resolution: float | None,
    information_limit: float | None,
) -> ResultFigure:
    """The simulated micrograph's spectrum, rotationally averaged, against the lens CTF."""

    q, average = radial_average(power_spectrum, pixel_size_angstrom=pixel_size_angstrom)
    frequency = np.asarray(frequencies, dtype=float)
    ctf = np.asarray(transfer, dtype=float)
    damping = np.asarray(envelope, dtype=float)

    def draw(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True)
        top.plot(q, average, lw=1.0, color=COLORS["data"])
        top.set_ylabel("Power spectrum,\nrotational average")
        bottom.plot(frequency, ctf, lw=1.1, color=COLORS["difference"], label="CTF sin χ · E")
        bottom.plot(frequency, damping, lw=1.0, ls="--", color=COLORS["model"], label="Envelope E")
        bottom.plot(frequency, -damping, lw=1.0, ls="--", color=COLORS["model"])
        bottom.axhline(0.0, color=COLORS["guide"], lw=0.7)
        for axes in (top, bottom):
            if point_resolution and point_resolution > 0:
                axes.axvline(1.0 / point_resolution, color=COLORS["accent"], lw=0.8, ls=":")
            if information_limit and information_limit > 0:
                axes.axvline(1.0 / information_limit, color=COLORS["warning"], lw=0.8, ls=":")
        bottom.set_ylabel("Transfer")
        bottom.set_xlabel("Spatial frequency q (Å⁻¹)")
        bottom.set_xlim(0.0, float(min(q.max(initial=1.0), frequency.max(initial=1.0))))
        bottom.legend(loc="upper right", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="hrem_spectrum_vs_ctf",
        title="What the lens passed: spectrum against transfer function",
        caption=(
            "Top: the rotational average of the simulated micrograph's power spectrum. "
            "Bottom: the objective-lens transfer function with its coherence envelope, for the "
            "same settings. Dotted lines: point resolution (green) and information limit (red)."
        ),
        interpretation=(
            "Spatial frequencies where the transfer crosses zero are missing from the image, "
            "and frequencies beyond the information limit are damped away; a Bragg spacing that "
            "falls in a zero or past the limit is not imaged at all, whatever the crystal. For a "
            "crystal the spectrum is dominated by Bragg peaks; an amorphous edge shows the "
            "rings of the transfer function directly."
        ),
        height_in=3.6,
    )


def thickness_figure(
    *,
    minima: Sequence[float],
    first_order: int,
    thickness_angstrom: float,
    extinction_angstrom: float,
) -> tuple[ResultFigure, dict[str, float]]:
    """The linearised two-beam fit, with every rejected order assignment for comparison.

    Returns the figure and the ordinary least-squares uncertainties of the chosen
    line propagated to the thickness and the extinction distance
    (``t = I^-1/2`` so ``σ(t)/t = σ(I)/2I``; likewise ``ξ = (−S)^-1/2``).
    """

    values = np.sort(np.abs(np.asarray(minima, dtype=float)))

    def line_for(order: int) -> tuple[np.ndarray, np.ndarray]:
        orders = np.arange(order, order + values.size, dtype=float)
        return 1.0 / orders**2, (values / orders) ** 2

    x, y = line_for(first_order)
    line = fit_line(x, y)
    statistics = {
        "sigma_thickness_angstrom": thickness_angstrom
        * line.sigma_intercept
        / (2 * line.intercept),
        "sigma_extinction_angstrom": extinction_angstrom * line.sigma_slope / (2 * abs(line.slope)),
    }

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        # The neighbouring assignments only: they are the ones a reader could
        # confuse with the chosen one, and further ones would set the scale.
        neighbours = [order for order in (first_order - 1, first_order + 1) if order >= 1]
        for index, order in enumerate(neighbours):
            cx, cy = line_for(order)
            axes.plot(
                cx,
                1e6 * cy,
                "o--",
                ms=3,
                lw=0.7,
                color=COLORS["guide"],
                label="Neighbouring order assignments" if index == 0 else None,
            )
        grid = np.linspace(0.0, float(x.max()) * 1.05, 50)
        if line.has_uncertainty:
            band = line.band(grid)
            axes.fill_between(
                grid,
                1e6 * (line.intercept + line.slope * grid - band),
                1e6 * (line.intercept + line.slope * grid + band),
                color=COLORS["band_3"],
                lw=0,
                label="±1σ of the line",
            )
        axes.plot(
            grid,
            1e6 * (line.intercept + line.slope * grid),
            lw=1.3,
            color=COLORS["model"],
            label=f"Fit, first minimum n = {first_order}",
        )
        axes.plot(x, 1e6 * y, "o", ms=6, color=COLORS["data"], label="Measured minima")
        axes.plot(
            [0.0], [1e6 * line.intercept], "s", ms=6, color=COLORS["accent"], label="Intercept 1/t²"
        )
        axes.set_xlabel("1/n²")
        axes.set_ylabel("(s_n/n)²  (10⁻⁶ Å⁻²)")
        axes.set_xlim(0.0, float(x.max()) * 1.6)
        top = max(float(np.max(y)), float(line.intercept)) * 1e6
        axes.set_ylim(0.0, top * 1.4)
        axes.legend(loc="best", frameon=False, fontsize=7)

    uncertainty = (
        f" t = {thickness_angstrom / 10:.1f} ± {statistics['sigma_thickness_angstrom'] / 10:.1f} "
        f"nm, ξg = {extinction_angstrom / 10:.1f} ± "
        f"{statistics['sigma_extinction_angstrom'] / 10:.1f} nm from the scatter about the line."
        if line.has_uncertainty
        else " With two minima the line passes through both points exactly, so the data cannot "
        "estimate an uncertainty; measure a third minimum for one."
    )
    return (
        render_figure(
            draw,
            key="thickness_fit",
            title="Two-beam thickness fit",
            caption=(
                "(s_n/n)² against 1/n² for the measured fringe minima (points), the straight "
                f"line for the chosen assignment, first minimum n = {first_order} (orange), "
                "and the same minima under the neighbouring assignments n ± 1 (grey). The "
                "intercept is 1/t² and the slope −1/ξg²." + uncertainty
            ),
            interpretation=(
                "Only the correct order assignment makes the points collinear; the others "
                "curve. With three or more minima the choice is made by the data; with two it "
                "rests on the first order you supply. Compare the fitted ξg with the value the "
                "structure predicts — a large difference means the orders are misassigned."
            ),
            height_in=3.2,
        ),
        statistics,
    )
