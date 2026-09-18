# ruff: noqa: RUF001, RUF002
"""Figures of the intermediate results of the XRD operations other than the cell.

Purpose
-------
Every XRD analysis turns a measured pattern into a number through a model, and
the reader can only trust the number after seeing the model laid over the data.
These builders draw that picture for each operation — the background under the
scan, the Rietveld profile and its weighted residuals, the instrument
calibration and the Williamson–Hall line, the candidates' lines against the
peaks, the reflections behind a simulated pattern — from the operation's own
result objects. Nothing here is computed that the operation did not compute,
except where a figure states a standard derived statistic (an ordinary
least-squares standard error, a residual divided by its weight), and those are
pinned by ``tests/unit/test_app_xrd_figures.py``.

The lattice-parameter report has its own module,
:mod:`pytex.app.services.xrd_lattice_report`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from pytex.app.figures import (
    COLORS,
    draw_normalized_residuals,
    draw_observed_model,
    draw_ticks,
    render_figure,
)
from pytex.app.results import ResultFigure, ResultMetric

__all__ = [
    "background_figure",
    "identification_figures",
    "pattern_figures",
    "rietveld_figures",
    "rietveld_weighted_residuals",
    "size_strain_figures",
    "size_strain_highlights",
    "williamson_hall_uncertainties",
]

_TWO_THETA = "2θ (°)"


# ---------------------------------------------------------------------------
# Background
# ---------------------------------------------------------------------------


def background_figure(
    two_theta_deg: np.ndarray,
    observed: np.ndarray,
    background: np.ndarray,
    *,
    method: str,
    fraction: float,
) -> ResultFigure:
    """The estimated background under the measured scan, and what is left above it."""

    axis = np.asarray(two_theta_deg, dtype=float)
    measured = np.asarray(observed, dtype=float)
    level = np.asarray(background, dtype=float)
    remainder = measured - level

    def draw(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[2.2, 1.0])
        top.plot(axis, measured, lw=0.6, color=COLORS["data"], label="Measured")
        top.plot(axis, level, lw=1.4, color=COLORS["model"], label=f"Background ({method})")
        top.set_yscale("log")
        top.set_ylabel("Intensity (log scale)")
        top.legend(loc="upper right", frameon=False)
        bottom.plot(axis, remainder, lw=0.6, color=COLORS["difference"])
        bottom.axhline(0.0, color=COLORS["guide"], lw=0.6)
        negative = remainder < 0.0
        if np.any(negative):
            bottom.plot(axis[negative], remainder[negative], ".", ms=1.5, color=COLORS["warning"])
        bottom.set_ylabel("Measured − background")
        bottom.set_xlabel(_TWO_THETA)
        bottom.set_xlim(float(axis[0]), float(axis[-1]))

    below = float(np.mean(remainder < 0.0)) if remainder.size else 0.0
    return render_figure(
        draw,
        key="background",
        title="Background under the scan",
        caption=(
            f"Measured scan (black, logarithmic scale so the background is visible) and the "
            f"{method} background estimate (orange); below, the measured intensity minus the "
            f"background, with negative points in red. The estimate assigns {100 * fraction:.1f} "
            "% of the total signal to background."
        ),
        interpretation=(
            f"{100.0 * below:.1f} % of points fall below the estimate. Between peaks the "
            "remainder should scatter about zero with the counting noise; a systematic dip under "
            "a peak means the estimate is clipping into the reflection and removing real "
            "intensity, and a hump between peaks means it is too low."
        ),
        height_in=3.8,
    )


# ---------------------------------------------------------------------------
# Rietveld
# ---------------------------------------------------------------------------


def rietveld_weighted_residuals(result: Any) -> np.ndarray | None:
    """``(y_obs − y_calc)·√w`` at every point, or ``None`` for unit weights.

    With the weights the refinement minimized, the sum of squares of these
    values divided by the degrees of freedom is the square of the goodness of
    fit. For unit weights the residual has no statistical scale, so no
    normalized residual is drawn rather than one with an invented scale.
    """

    model = str(result.weight_model)
    observed = np.asarray(result.observed_intensity, dtype=float)
    residual = np.asarray(result.residual_intensity, dtype=float)
    if model == "poisson":
        weights = 1.0 / np.clip(observed, 1.0, None)
    elif model == "inverse_variance":
        measured = result.measured
        if measured.standard_uncertainty is None:
            return None
        full_axis = np.asarray(measured.two_theta_deg, dtype=float)
        index = np.searchsorted(full_axis, np.asarray(result.two_theta_deg, dtype=float))
        index = np.clip(index, 0, full_axis.size - 1)
        sigma = np.asarray(measured.standard_uncertainty, dtype=float)[index]
        weights = 1.0 / np.square(np.clip(sigma, 1e-12, None))
    else:
        return None
    weighted: np.ndarray = residual * np.sqrt(weights)
    return weighted


def rietveld_figures(result: Any, *, labels: Sequence[str]) -> tuple[ResultFigure, ...]:
    """Observed/calculated/background/difference, weighted residuals, parameter shifts."""

    axis = np.asarray(result.two_theta_deg, dtype=float)
    observed = np.asarray(result.observed_intensity, dtype=float)
    calculated = np.asarray(result.calculated_intensity, dtype=float)
    background = np.asarray(result.background_intensity, dtype=float)
    positions = [float(item.two_theta_deg) for item in result.reflections]

    def draw_profile(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[3.0, 1.0])
        draw_observed_model(
            top,
            axis,
            observed,
            calculated,
            difference_axes=bottom,
            observed_label="Observed",
            model_label="Calculated",
            observed_as_points=True,
        )
        top.plot(axis, background, lw=0.9, ls="--", color=COLORS["background"], label="Background")
        ceiling = float(np.max(observed)) if observed.size else 1.0
        top.set_ylim(bottom=-0.18 * ceiling, top=1.08 * ceiling)
        draw_ticks(top, positions, labels=labels if len(labels) <= 40 else None, height=0.06)
        top.set_yticks([tick for tick in top.get_yticks() if 0.0 <= tick <= 1.08 * ceiling])
        top.set_ylabel("Intensity")
        top.legend(loc="upper right", frameon=False)
        bottom.set_xlabel(_TWO_THETA)
        bottom.set_xlim(float(axis[0]), float(axis[-1]))

    figures = [
        render_figure(
            draw_profile,
            key="rietveld_profile",
            title="Rietveld fit: observed, calculated, difference",
            caption=(
                "Measured points (black), calculated profile (orange), refined background "
                "(dashed), reflection positions of the refined cell (ticks) and observed minus "
                f"calculated below. R_wp = {100 * result.weighted_profile_r_factor:.2f} %, "
                f"R_exp = {100 * result.expected_r_factor:.2f} %, goodness of fit "
                f"{result.goodness_of_fit:.3g}."
            ),
            interpretation=(
                "A flat difference curve at the level of the noise is a model that describes "
                "the pattern. A derivative-shaped wiggle at every peak is a cell or zero error; "
                "a symmetric dip or peak is a width or shape error; a difference that follows "
                "the intensities of some reflections only is texture or a wrong structure."
            ),
            height_in=4.2,
        )
    ]

    weighted = rietveld_weighted_residuals(result)
    if weighted is not None:

        def draw_weighted(figure: Any) -> None:
            axes = figure.subplots()
            draw_normalized_residuals(
                axes,
                axis,
                weighted,
                xlabel=_TWO_THETA,
                ylabel="(y_obs − y_calc) / σ",
                marker_size=1.6,
            )

        inside = float(np.mean(np.abs(weighted) <= 2.0)) if weighted.size else 0.0
        figures.append(
            render_figure(
                draw_weighted,
                key="rietveld_weighted_residuals",
                title="Weighted residuals",
                caption=(
                    "Each point's residual divided by its counting standard uncertainty "
                    f"({result.weight_model} weights, the ones the refinement minimized), with "
                    "±2σ and ±3σ bands. The mean square of these values per degree of freedom "
                    f"is the square of the goodness of fit, {result.goodness_of_fit:.3g}² = "
                    f"{result.goodness_of_fit**2:.3g}."
                ),
                interpretation=(
                    f"{100.0 * inside:.1f} % of points lie within ±2σ (about 95 % is expected "
                    "for a model that describes the pattern within the noise). Runs of "
                    "same-sign residuals, which the Durbin–Watson statistic "
                    f"({result.durbin_watson:.3g}; 2 means uncorrelated) measures, show where "
                    "the misfit is systematic."
                ),
                height_in=3.0,
            )
        )

    refined = [item for item in result.parameters if item.refined]
    shifts = [
        (item.name, (item.value - item.initial_value) / item.standard_uncertainty)
        for item in refined
        if item.standard_uncertainty
    ]
    if shifts:

        def draw_shifts(figure: Any) -> None:
            axes = figure.subplots()
            names = [name for name, _ in shifts]
            values = np.array([value for _, value in shifts])
            colors = [COLORS["warning"] if abs(value) > 3.0 else COLORS["data"] for value in values]
            positions_y = np.arange(len(names))
            axes.barh(positions_y, values, color=colors, height=0.6)
            axes.axvline(0.0, color=COLORS["guide"], lw=0.8)
            for level in (-3.0, 3.0):
                axes.axvline(level, color=COLORS["guide"], lw=0.7, ls=":")
            axes.set_yticks(positions_y, labels=names)
            axes.invert_yaxis()
            axes.set_xscale("symlog", linthresh=3.0)
            axes.set_xlabel("(refined − starting value) / σ")

        figures.append(
            render_figure(
                draw_shifts,
                key="rietveld_parameter_shifts",
                title="How far each parameter moved",
                caption=(
                    "Each refined parameter's shift from its starting value in units of its "
                    "own standard uncertainty (symmetric-log scale; dotted lines at ±3σ). Red "
                    "bars moved by more than 3σ."
                ),
                interpretation=(
                    "A parameter that moved many σ was genuinely determined by the data; one "
                    "that barely moved either started at the answer or is not constrained by "
                    "this pattern. The uncertainties are the fit covariance scaled by the "
                    "goodness of fit."
                ),
                height_in=0.7 + 0.32 * len(shifts),
            )
        )
    return tuple(figures)


# ---------------------------------------------------------------------------
# Size and strain
# ---------------------------------------------------------------------------


def williamson_hall_uncertainties(
    abscissa: np.ndarray,
    ordinate: np.ndarray,
    *,
    wavelength_angstrom: float,
    shape_factor: float,
) -> dict[str, float]:
    """Ordinary least-squares standard errors of the Williamson–Hall line, propagated.

    The line is fitted unweighted, so the standard errors follow from the
    scatter about it: ``s² = Σr²/(n − 2)``, ``σ(slope)² = s²/Sxx`` and
    ``σ(intercept)² = s²(1/n + x̄²/Sxx)``. The crystallite size is
    ``D = Kλ/intercept``, so ``σ(D)/D = σ(intercept)/intercept``; the microstrain
    is the slope itself. With two reflections the line has no scatter to
    estimate, and every uncertainty is ``nan``.
    """

    x = np.asarray(abscissa, dtype=float)
    y = np.asarray(ordinate, dtype=float)
    n = x.size
    slope, intercept = np.polyfit(x, y, 1)
    if n <= 2:
        nan = float("nan")
        return {
            "slope": float(slope),
            "intercept": float(intercept),
            "sigma_slope": nan,
            "sigma_intercept": nan,
            "size_nm": float(shape_factor * wavelength_angstrom / intercept / 10.0),
            "sigma_size_nm": nan,
            "covariance": nan,
        }
    residual = y - (slope * x + intercept)
    variance = float(np.sum(residual**2) / (n - 2))
    mean = float(np.mean(x))
    sxx = float(np.sum((x - mean) ** 2))
    sigma_slope = float(np.sqrt(variance / sxx))
    sigma_intercept = float(np.sqrt(variance * (1.0 / n + mean**2 / sxx)))
    size = float(shape_factor * wavelength_angstrom / intercept / 10.0)
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "sigma_slope": sigma_slope,
        "sigma_intercept": sigma_intercept,
        "size_nm": size,
        "sigma_size_nm": abs(size) * sigma_intercept / abs(float(intercept)),
        "covariance": -mean * variance / sxx,
    }


def size_strain_figures(
    *,
    standard_angles: np.ndarray,
    standard_widths: np.ndarray,
    instrument: Any,
    sample_angles: np.ndarray,
    observed_widths: np.ndarray,
    sample_widths: np.ndarray,
    abscissa: np.ndarray,
    ordinate: np.ndarray,
    scherrer_nm: np.ndarray,
    statistics: Mapping[str, float],
) -> tuple[ResultFigure, ...]:
    """Instrument calibration, width decomposition, Williamson–Hall line, Scherrer sizes."""

    std_x = np.asarray(standard_angles, dtype=float)
    std_y = np.asarray(standard_widths, dtype=float)
    low = float(min(np.min(std_x), np.min(sample_angles))) - 2.0
    high = float(max(np.max(std_x), np.max(sample_angles))) + 2.0
    grid = np.linspace(max(low, 1.0), min(high, 179.0), 300)
    curve = np.asarray(instrument.fwhm_deg(grid), dtype=float)

    def draw_calibration(figure: Any) -> None:
        axes = figure.subplots()
        axes.plot(std_x, std_y, "o", ms=5, color=COLORS["data"], label="Standard peaks")
        axes.plot(grid, curve, lw=1.4, color=COLORS["model"], label="Fitted Caglioti function")
        axes.set_xlabel(_TWO_THETA)
        axes.set_ylabel("FWHM (°)")
        axes.legend(loc="best", frameon=False)

    fitted_at_standard = np.asarray(instrument.fwhm_deg(std_x), dtype=float)
    rms = float(np.sqrt(np.mean((std_y - fitted_at_standard) ** 2))) if std_y.size else 0.0

    def draw_widths(figure: Any) -> None:
        axes = figure.subplots()
        axes.plot(
            sample_angles,
            observed_widths,
            "o",
            ms=5,
            color=COLORS["data"],
            label="Measured specimen width",
        )
        axes.plot(
            grid, curve, lw=1.0, ls="--", color=COLORS["background"], label="Instrument width"
        )
        axes.plot(
            sample_angles,
            sample_widths,
            "s",
            ms=5,
            color=COLORS["model"],
            label="Specimen-only width (deconvolved)",
        )
        axes.set_xlabel(_TWO_THETA)
        axes.set_ylabel("FWHM (°)")
        axes.legend(loc="best", frameon=False)

    x = np.asarray(abscissa, dtype=float)
    y = np.asarray(ordinate, dtype=float)
    slope = statistics["slope"]
    intercept = statistics["intercept"]
    line_grid = np.linspace(0.0, float(np.max(x)) * 1.05, 100)
    variance_line = (
        statistics["sigma_intercept"] ** 2
        + (line_grid**2) * statistics["sigma_slope"] ** 2
        + 2.0 * line_grid * statistics["covariance"]
    )
    band = np.sqrt(np.maximum(variance_line, 0.0)) if np.isfinite(variance_line).all() else None

    def draw_wh(figure: Any) -> None:
        axes = figure.subplots()
        if band is not None:
            axes.fill_between(
                line_grid,
                1000.0 * (intercept + slope * line_grid - band),
                1000.0 * (intercept + slope * line_grid + band),
                color=COLORS["band_3"],
                lw=0,
                label="±1σ of the fitted line",
            )
        axes.plot(
            line_grid,
            1000.0 * (intercept + slope * line_grid),
            lw=1.4,
            color=COLORS["model"],
            label="Least-squares line",
        )
        axes.plot(x, 1000.0 * y, "o", ms=5, color=COLORS["data"], label="Reflection")
        axes.errorbar(
            [0.0],
            [1000.0 * intercept],
            yerr=[1000.0 * statistics["sigma_intercept"]]
            if np.isfinite(statistics["sigma_intercept"])
            else None,
            fmt="s",
            ms=5,
            capsize=3,
            color=COLORS["accent"],
            label="Intercept Kλ/D",
        )
        axes.set_xlabel("4 sinθ")
        axes.set_ylabel("β cosθ (mrad)")
        axes.set_xlim(left=-0.05 * float(np.max(x)))
        axes.legend(loc="best", frameon=False, fontsize=7)

    def draw_scherrer(figure: Any) -> None:
        axes = figure.subplots()
        axes.plot(
            sample_angles,
            scherrer_nm,
            "o",
            ms=5,
            color=COLORS["data"],
            label="Scherrer size of each reflection",
        )
        size = statistics["size_nm"]
        sigma = statistics["sigma_size_nm"]
        axes.axhline(size, color=COLORS["model"], lw=1.2, label="Williamson–Hall size")
        if np.isfinite(sigma):
            axes.axhspan(size - sigma, size + sigma, color=COLORS["band_3"], lw=0)
        axes.set_xlabel(_TWO_THETA)
        axes.set_ylabel("Crystallite size (nm)")
        axes.legend(loc="best", frameon=False)

    trend = "fall" if np.polyfit(sample_angles, scherrer_nm, 1)[0] < 0 else "rise"
    return (
        render_figure(
            draw_calibration,
            key="instrument_calibration",
            title="Instrument resolution from the standard",
            caption=(
                "FWHM of each standard peak (points) and the fitted Caglioti function "
                f"FWHM² = U tan²θ + V tanθ + W (line). RMS misfit {rms:.4f}°."
            ),
            interpretation=(
                "The curve is the width the diffractometer alone gives each angle. It must pass "
                "through the standard's points and span the specimen's angular range; outside "
                "that range it is an extrapolation."
            ),
            height_in=3.0,
        ),
        render_figure(
            draw_widths,
            key="width_decomposition",
            title="Measured, instrumental and specimen-only widths",
            caption=(
                "Measured specimen FWHM (circles), instrument FWHM from the calibration "
                "(dashed) and the specimen-only width after deconvolution (squares)."
            ),
            interpretation=(
                "Where the measured width is close to the instrument width, the specimen-only "
                "width is a small difference of two large numbers and is poorly determined; "
                "those reflections carry little information about size or strain."
            ),
            height_in=3.0,
        ),
        render_figure(
            draw_wh,
            key="williamson_hall",
            title="Williamson–Hall plot",
            caption=(
                "β cosθ against 4 sinθ for each reflection (β the specimen-only FWHM in "
                "radians), the unweighted least-squares line and its ±1σ band from the scatter "
                f"about it. Slope ε = {slope:.3e} ± {statistics['sigma_slope']:.1e}; intercept "
                f"Kλ/D = {1000 * intercept:.3f} ± {1000 * statistics['sigma_intercept']:.3f} "
                "mrad."
            ),
            interpretation=(
                "The intercept is size broadening and the slope is strain broadening. Points "
                "scattered far from a line mean the uniform-size, uniform-strain model does "
                "not hold (anisotropic broadening, stacking faults); a negative slope is "
                "unphysical and usually means the instrument width was over-subtracted."
            ),
            height_in=3.2,
        ),
        render_figure(
            draw_scherrer,
            key="scherrer_sizes",
            title="Scherrer size of each reflection",
            caption=(
                "Crystallite size from each reflection alone by the Scherrer equation, "
                "against angle, with the Williamson–Hall size and its ±1σ band."
            ),
            interpretation=(
                f"The single-reflection sizes {trend} with angle. Sizes that fall with angle "
                "mean part of the broadening is strain, which the Scherrer equation wrongly "
                "reads as small crystallites; that is why the Williamson–Hall intercept is the "
                "size to quote."
            ),
            height_in=3.0,
        ),
    )


# ---------------------------------------------------------------------------
# Phase identification
# ---------------------------------------------------------------------------


def identification_figures(data: Mapping[str, Any]) -> tuple[ResultFigure, ...]:
    """The scan with every candidate's lines under it, and the scores broken down."""

    axis = np.asarray(data["two_theta_deg"], dtype=float)
    observed = np.asarray(data["observed"], dtype=float)
    overlays = list(data.get("overlays") or ())
    candidates = list(data.get("candidates") or ())
    peaks = list(data.get("peaks") or ())
    palette = ["#d97706", "#2563eb", "#0f766e", "#7c3aed", "#db2777", "#65a30d", "#0891b2"]

    def draw_lines(figure: Any) -> None:
        rows = max(1, len(overlays))
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[2.4, 0.35 * rows + 0.4])
        top.plot(axis, observed, lw=0.6, color=COLORS["data"], label="Measured")
        if peaks:
            angles = [float(peak["two_theta_deg"]) for peak in peaks]
            heights = np.interp(angles, axis, observed)
            top.plot(
                angles, heights * 1.04, "v", ms=4, color=COLORS["accent"], label="Detected peak"
            )
        top.set_ylabel("Intensity")
        top.legend(loc="upper right", frameon=False)
        for row, overlay in enumerate(overlays):
            color = palette[row % len(palette)]
            positions = np.asarray(overlay.get("two_theta_deg") or [], dtype=float)
            strengths = np.asarray(overlay.get("relative_intensity") or [], dtype=float)
            if positions.size:
                scale = strengths / (float(np.max(strengths)) or 1.0) if strengths.size else 1.0
                bottom.vlines(positions, row, row + 0.8 * scale, color=color, lw=1.3)
        bottom.set_yticks(
            [row + 0.4 for row in range(len(overlays))],
            labels=[
                f"{overlay['phase_name']} ({float(overlay['score']):.2f})" for overlay in overlays
            ],
        )
        bottom.set_ylim(-0.2, max(1, len(overlays)))
        bottom.invert_yaxis()
        bottom.set_xlabel(_TWO_THETA)
        bottom.set_xlim(float(axis[0]), float(axis[-1]))

    figures = [
        render_figure(
            draw_lines,
            key="identification_lines",
            title="Measured peaks against every candidate's lines",
            caption=(
                "Measured scan with detected peaks (top), and below it one row per candidate, "
                "in rank order, with its calculated lines at the cell it was fitted to (line "
                "height = calculated relative intensity; the score is in brackets)."
            ),
            interpretation=(
                "The right phase puts a line under every strong peak and few lines where there "
                "is none. A strong peak with no line under any candidate is a phase not on the "
                "list; a candidate whose strong lines fall on empty scan is ruled out whatever "
                "its score."
            ),
            height_in=3.4 + 0.3 * len(overlays),
        )
    ]
    criteria = (
        ("explained_intensity_fraction", "Intensity explained"),
        ("completeness", "Lines seen"),
        ("position_score", "Position"),
        ("intensity_agreement", "Intensity agreement"),
    )
    if candidates:

        def draw_scores(figure: Any) -> None:
            axes = figure.subplots()
            names = [str(item.get("phase_name")) for item in candidates]
            width = 0.8 / len(criteria)
            positions = np.arange(len(names))
            for index, (key, label) in enumerate(criteria):
                values = [float(item.get(key) or 0.0) for item in candidates]
                axes.bar(
                    positions + (index - 1.5) * width,
                    values,
                    width,
                    color=palette[index],
                    label=label,
                )
            scores = [float(item.get("score") or 0.0) for item in candidates]
            axes.plot(positions, scores, "D", ms=6, color=COLORS["data"], label="Total score")
            axes.set_xticks(positions, labels=names, rotation=20, ha="right")
            axes.set_ylim(0.0, 1.05)
            axes.set_ylabel("Criterion value (0 to 1)")
            axes.legend(loc="upper right", frameon=False, fontsize=7, ncols=2)

        figures.append(
            render_figure(
                draw_scores,
                key="identification_scores",
                title="What each score is made of",
                caption=(
                    "The four criteria behind each candidate's score (bars) and the weighted "
                    "total (diamonds), candidates in rank order."
                ),
                interpretation=(
                    "Read the criterion a candidate fails, not only its total: a low intensity "
                    "agreement with good positions is texture, not a wrong phase; low lines-seen "
                    "with good positions is a different centring of the same cell."
                ),
                height_in=3.2,
            )
        )
    return tuple(figures)


# ---------------------------------------------------------------------------
# Simulated pattern
# ---------------------------------------------------------------------------


def pattern_figures(
    two_theta_deg: np.ndarray,
    intensity: np.ndarray,
    rows: Sequence[Mapping[str, Any]],
    *,
    labels: Sequence[str],
    radiation: str,
) -> tuple[ResultFigure, ...]:
    """The simulated pattern with its reflections, and the factors behind each intensity."""

    axis = np.asarray(two_theta_deg, dtype=float)
    profile = np.asarray(intensity, dtype=float)
    positions = np.array([float(row["two_theta_deg"]) for row in rows])
    relative = np.array([float(row["relative_intensity"]) for row in rows])
    structure = np.array([float(row["structure_factor_amplitude"]) for row in rows])
    multiplicity = np.array([float(row["multiplicity"]) for row in rows])
    lorentz = np.array([float(row["lorentz_polarization"]) for row in rows])

    def draw_pattern(figure: Any) -> None:
        axes = figure.subplots()
        ceiling = float(np.max(profile)) if profile.size else 1.0
        axes.plot(axis, profile / ceiling, lw=0.9, color=COLORS["data"], label="Profile")
        axes.vlines(
            positions,
            0.0,
            relative,
            color=COLORS["model"],
            lw=1.2,
            label="Integrated intensity (relative)",
        )
        for x, y, text in zip(positions, relative, labels, strict=False):
            if y >= 0.03:
                axes.annotate(
                    text,
                    (x, y),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    fontsize=6.5,
                    rotation=90,
                    color=COLORS["background"],
                )
        axes.set_ylim(0.0, 1.25)
        axes.set_xlabel(_TWO_THETA)
        axes.set_ylabel("Relative intensity")
        axes.legend(loc="upper right", frameon=False)

    def draw_factors(figure: Any) -> None:
        axes = figure.subplots()
        with np.errstate(divide="ignore"):
            series = (
                ("|F|²", np.square(structure), "o"),
                ("Multiplicity", multiplicity, "s"),
                ("Lorentz–polarization", lorentz, "^"),
                (
                    "Integrated intensity",
                    relative * float(np.max(np.square(structure)) or 1.0),
                    "D",
                ),
            )
        for (label, values, marker), color in zip(
            series,
            (COLORS["model"], COLORS["accent"], COLORS["difference"], COLORS["data"]),
            strict=True,
        ):
            finite = values > 0
            axes.plot(positions[finite], values[finite], marker, ms=4, color=color, label=label)
        axes.set_yscale("log")
        axes.set_xlabel(_TWO_THETA)
        axes.set_ylabel("Factor (log scale)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return (
        render_figure(
            draw_pattern,
            key="powder_pattern",
            title="Simulated powder pattern",
            caption=(
                f"Calculated profile for {radiation} (line, normalized to its maximum) and the "
                "integrated intensity of each reflection (sticks), labelled where visible."
            ),
            interpretation=(
                "Overlapping reflections merge into one profile peak, so the tallest peak need "
                "not belong to the strongest single reflection. Reflections allowed by the "
                "lattice but absent here are extinguished by the structure factor."
            ),
            height_in=3.2,
        ),
        render_figure(
            draw_factors,
            key="intensity_factors",
            title="What sets each reflection's intensity",
            caption=(
                "Per reflection: the squared structure-factor amplitude |F|², the multiplicity, "
                "the Lorentz–polarization factor, and the resulting integrated intensity (scaled "
                "onto the |F|² axis), on a logarithmic scale."
            ),
            interpretation=(
                "The integrated intensity is the product of these factors (and a temperature "
                "factor where one is modelled). A strong line with a small |F|² owes its "
                "strength to multiplicity or to the Lorentz factor at low angle."
            ),
            height_in=3.0,
        ),
    )


def size_strain_highlights(
    statistics: Mapping[str, float], r_squared: float
) -> tuple[ResultMetric, ...]:
    """Size and strain with the ordinary least-squares uncertainties of the line."""

    return (
        ResultMetric(
            "Crystallite size D",
            f"{statistics['size_nm']:.4g} ± {statistics['sigma_size_nm']:.2g}",
            "nm",
            "From the Williamson–Hall intercept; ± from the scatter of the points about the line.",
        ),
        ResultMetric(
            "Microstrain ε",
            f"{statistics['slope']:.3e} ± {statistics['sigma_slope']:.1e}",
            None,
            "The Williamson–Hall slope.",
        ),
        ResultMetric("R² of the line", float(r_squared)),
    )
