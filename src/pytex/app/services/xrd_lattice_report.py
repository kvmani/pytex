# ruff: noqa: RUF001, RUF002
"""The report of a lattice-parameter determination: headline, warnings and figures.

Purpose
-------
:func:`pytex.diffraction.xrd_lattice_parameter.lattice_parameter_pipeline` does
the science. This module says what it found, in the order a reader needs it —
the cell and how far to trust it, the evidence, the diagnostics, the method —
and draws the figures that let the reader check each step. It computes no new
science: every number here is either read from the pipeline's result objects
or is an algebraic rearrangement of them (a residual divided by its own
uncertainty; the fitted correction evaluated on a grid), and the tests hold
those rearrangements to the quantities they rearrange.

Vocabulary kept distinct on purpose
-----------------------------------
- **Peak-fit chi-squared** judges one profile against the counts in one window.
  **Lattice-fit chi-squared** judges the fitted peak *positions* against the
  cell. They answer different questions and are never reported under one name.
- **Precision** is the standard uncertainty from the fit: how reproducible the
  number is on this scan. **Accuracy** needs a calibration against a certified
  standard; no fit can supply it.
- The **relative change from the reference cell** compares with the tabulated
  cell of the selected phase. It includes composition, temperature and
  instrument calibration differences, and is an elastic strain only when the
  reference is the stress-free cell of the same material on the same instrument.
- The Cohen / Nelson–Riley term is an **angle-dependent systematic
  correction**. Its angular form matches specimen displacement and absorption,
  but the fit cannot tell which aberration produced it, so it is not a measured
  displacement.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from matplotlib.ticker import MaxNLocator

from pytex.app.figures import (
    COLORS,
    draw_correlation,
    draw_normalized_residuals,
    draw_observed_model,
    draw_residuals,
    draw_ticks,
    render_figure,
)
from pytex.app.results import ResultFigure, ResultMetric
from pytex.diffraction.xrd_lattice_parameter import (
    LatticeParameterResult,
    extrapolation_values,
)
from pytex.diffraction.xrd_peaks import kalpha_doublet_parameters

__all__ = [
    "EXTRAPOLATION_NAMES",
    "correction_figure",
    "correlation_figure",
    "cross_check_figure",
    "extrapolation_figure",
    "indexing_figure",
    "lattice_highlights",
    "lattice_warnings",
    "le_bail_figure",
    "normalized_residual_figure",
    "normalized_residuals",
    "peak_quality_figure",
    "peak_windows_figure",
    "residual_figure",
    "scan_figure",
    "strongest_correlation",
    "systematic_correction_curve",
]

#: Plain names of the extrapolation functions, for prose.
EXTRAPOLATION_NAMES = {
    "nelson_riley": "Nelson–Riley",
    "cos_squared_over_sin": "cos²θ/sinθ",
    "cot_theta": "cotθ",
    "bradley_jay": "Bradley–Jay cos²θ",
    "none": "no correction",
}

#: A lattice-fit reduced chi-squared outside this band gets a warning.
CHI_SQUARED_HIGH = 3.0
CHI_SQUARED_LOW = 0.3
#: A correlation this close to +/-1 means two parameters are barely separable.
CORRELATION_LIMIT = 0.95
#: Fewer degrees of freedom than this and chi-squared is itself very uncertain.
MINIMUM_DEGREES = 3


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------


def normalized_residuals(
    result: LatticeParameterResult, sigma_two_theta_deg: np.ndarray
) -> np.ndarray:
    """Return ``(2θobs − 2θcalc) / σ(2θ)`` for each reflection of a position fit.

    ``2θcalc`` includes the fitted systematic correction, so this is the
    residual the lattice fit left. For Cohen's method the sum of squares of
    these values divided by the degrees of freedom *is* the reported reduced
    chi-squared, because the residual and its uncertainty are converted from
    ``sin²θ`` to ``2θ`` by the same derivative.
    """

    sigma = np.asarray(sigma_two_theta_deg, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        values = np.asarray(result.residual_two_theta_deg, dtype=float) / sigma
    return np.where(sigma > 0.0, values, np.nan)


def systematic_correction_curve(
    result: LatticeParameterResult, two_theta_deg: Sequence[float] | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The fitted systematic correction and its standard uncertainty at any angle.

    The same expression as
    :attr:`~pytex.diffraction.xrd_lattice_parameter.LatticeParameterResult.systematic_shift_deg`,
    evaluated on an arbitrary grid so it can be drawn as a curve: the drift
    column ``D·sin²θ·f(θ)`` converted from ``sin²θ`` to degrees ``2θ``. The
    correction is linear in ``D``, so its uncertainty is ``|shift|·σ(D)/|D|``.

    Returns
    -------
    tuple of numpy.ndarray
        ``(shift_deg, sigma_deg)``; both zero when no term was refined.
    """

    grid = np.asarray(two_theta_deg, dtype=float)
    if result.extrapolation == "none" or result.drift_coefficient == 0.0:
        zeros = np.zeros_like(grid)
        return zeros, zeros
    function = extrapolation_values(grid, function=result.extrapolation)
    theta = np.deg2rad(0.5 * grid)
    delta = result.drift_coefficient * np.square(np.sin(theta)) * function
    shift = np.rad2deg(2.0 * delta / np.sin(2.0 * theta))
    sigma = np.abs(shift) * result.drift_standard_uncertainty / abs(result.drift_coefficient)
    return shift, sigma


def strongest_correlation(result: LatticeParameterResult) -> tuple[float, str, str] | None:
    """The largest off-diagonal correlation, signed, and the two parameters it links."""

    matrix = result.parameter_correlation
    if matrix is None or matrix.shape[0] < 2:
        return None
    off = np.abs(matrix - np.eye(matrix.shape[0]))
    row, column = np.unravel_index(int(np.argmax(off)), off.shape)
    names = result.correlation_parameter_names
    return float(matrix[row, column]), names[int(row)], names[int(column)]


def _degrees_of_freedom(result: LatticeParameterResult) -> int:
    if result.method == "average":
        return int(result.reflection_count) - 1
    if result.method == "le_bail":
        return -1
    free = len(result.free_parameter_names) + (0 if result.extrapolation == "none" else 1)
    return int(result.reflection_count) - free


def _cell_text(value: float, sigma: float, digits: int = 6) -> str:
    return f"{value:.{digits}f} ± {sigma:.{digits}f}"


# ---------------------------------------------------------------------------
# Headline and warnings
# ---------------------------------------------------------------------------


def lattice_highlights(
    result: LatticeParameterResult,
    *,
    system: str,
    profile_points: int | None = None,
) -> tuple[ResultMetric, ...]:
    """The answer and the numbers that decide how far to trust it, in that order."""

    metrics: list[ResultMetric] = [
        ResultMetric(
            "a",
            _cell_text(result.a, result.a_standard_uncertainty),
            "Å",
            "Cell edge ± one standard uncertainty (precision on this scan, not accuracy).",
        )
    ]
    if system in {"orthorhombic", "monoclinic", "triclinic"}:
        metrics.append(ResultMetric("b", _cell_text(result.b, result.b_standard_uncertainty), "Å"))
    if system != "cubic":
        metrics.append(ResultMetric("c", _cell_text(result.c, result.c_standard_uncertainty), "Å"))
        metrics.append(ResultMetric("c/a", f"{result.axial_ratio:.6f}"))
    if system in {"monoclinic", "triclinic"}:
        metrics.append(
            ResultMetric(
                "α, β, γ",
                f"{result.alpha_deg:.4f}, {result.beta_deg:.4f}, {result.gamma_deg:.4f}",
                "°",
            )
        )
    metrics.append(
        ResultMetric(
            "Relative precision σ(a)/a",
            float(result.relative_uncertainty),
            None,
            "About 1e-5 is strain-grade; 1e-4 composition-grade; 1e-3 identification-grade.",
        )
    )
    if result.method == "le_bail":
        metrics.append(ResultMetric("Reflections modelled", int(result.reflection_count)))
        if profile_points is not None:
            metrics.append(ResultMetric("Profile points fitted", int(profile_points)))
        metrics.append(
            ResultMetric(
                "Profile-fit reduced χ²",
                float(result.reduced_chi_squared),
                None,
                "Whole-pattern fit of the calculated profile to every measured point.",
            )
        )
        if result.weighted_profile_r is not None:
            metrics.append(
                ResultMetric(
                    "R_wp (background removed)",
                    float(result.weighted_profile_r),
                    None,
                    "Computed on the background-subtracted profile; higher than an R_wp on raw "
                    "counts and not comparable with one.",
                )
            )
    else:
        degrees = _degrees_of_freedom(result)
        metrics += [
            ResultMetric("Reflections used", int(result.reflection_count)),
            ResultMetric(
                "Degrees of freedom",
                int(degrees),
                None,
                "Reflections minus refined parameters. Below about 3, χ² and the "
                "uncertainties are themselves poorly determined.",
            ),
            ResultMetric(
                "Lattice-fit reduced χ²",
                float(result.reduced_chi_squared),
                None,
                "Peak positions against the cell, in units of their own uncertainties. "
                "About 1 is ideal. Not the same as a peak-profile χ².",
            ),
        ]
        strongest = strongest_correlation(result)
        if strongest is not None:
            value, first, second = strongest
            metrics.append(
                ResultMetric(
                    "Strongest parameter correlation",
                    f"{value:+.3f} ({first}–{second})",
                    None,
                    "Beyond ±0.95 the scan barely separates the two parameters.",
                )
            )
        if result.extrapolation != "none" and result.method == "cohen":
            significance = (
                abs(result.drift_coefficient) / result.drift_standard_uncertainty
                if result.drift_standard_uncertainty > 0.0
                else float("nan")
            )
            largest = (
                float(np.max(np.abs(result.systematic_shift_deg))) * 1000.0
                if result.two_theta_deg.size
                else 0.0
            )
            metrics.append(
                ResultMetric(
                    "Systematic correction",
                    f"{EXTRAPOLATION_NAMES[result.extrapolation]}: D = "
                    f"{result.drift_coefficient:.3e} ± {result.drift_standard_uncertainty:.1e} "
                    f"({significance:.1f}σ), up to {largest:.1f} m° in 2θ",
                    None,
                    "An angle-dependent correction fitted with the cell. It is not a measured "
                    "specimen displacement.",
                )
            )
        else:
            metrics.append(ResultMetric("Systematic correction", "none refined"))
    change = result.strain_relative_to_reference
    if change is not None and result.reference_lattice is not None:
        metrics.append(
            ResultMetric(
                "Change from the reference cell (a)",
                float(change),
                None,
                f"(a − a_ref)/a_ref against the tabulated a_ref = "
                f"{result.reference_lattice.a:.6f} Å. Not an elastic strain unless the "
                "reference is the stress-free cell of this material on this instrument.",
            )
        )
    return tuple(metrics)


def lattice_warnings(
    result: LatticeParameterResult,
    *,
    normalized: np.ndarray | None = None,
    labels: Sequence[str] = (),
    peaks_used: Sequence[Any] = (),
    figure_of_merit: float | None = None,
    unindexed_count: int = 0,
    generated: bool = False,
) -> tuple[str, ...]:
    """Every reason, found in this run, to distrust the reported cell.

    Each warning says what was found, what it usually means, and what to do.
    ``generated`` is accepted for symmetry with the notes and raises nothing:
    a demonstration scan is announced in the notes, and is not a defect of the
    fit.
    """

    warnings: list[str] = []
    chi = float(result.reduced_chi_squared)
    if result.method == "le_bail":
        if chi > CHI_SQUARED_HIGH:
            warnings.append(
                f"The whole-pattern fit has a reduced χ² of {chi:.2f}: the calculated profile "
                "does not describe the measured one within the noise. Look for structure in the "
                "difference curve — an extra phase, a wrong peak shape, or a systematic term "
                "that is not the right one."
            )
        return tuple(warnings)

    degrees = _degrees_of_freedom(result)
    if degrees <= 0:
        warnings.append(
            f"The fit is exactly determined ({result.reflection_count} reflections for as many "
            "parameters). There is no redundancy, so χ² and every quoted uncertainty are "
            "meaningless. Add reflections (lower the detection threshold or the angular floor) "
            "or refine fewer parameters."
        )
    elif degrees < MINIMUM_DEGREES:
        warnings.append(
            f"Only {degrees} degree{'s' if degrees != 1 else ''} of freedom. With so few "
            "surplus reflections the reduced χ² and the uncertainties are themselves very "
            "uncertain; treat σ as indicative."
        )
    if result.method == "average":
        warnings.append(
            "The average method cannot remove an angle-dependent systematic error; the result "
            "carries whatever zero or displacement error the scan has. Use Cohen's method for a "
            "reportable value."
        )
    elif result.extrapolation == "none":
        warnings.append(
            "No systematic correction was refined. Any zero-offset, displacement or "
            "transparency error is absorbed into the cell."
        )
    if degrees > 0 and chi > CHI_SQUARED_HIGH:
        warnings.append(
            f"Poor lattice fit: reduced χ² = {chi:.2f}. The peak positions scatter about the "
            f"cell {np.sqrt(chi):.1f} times more than their own uncertainties allow. The quoted "
            "uncertainties have been enlarged by that factor, but the cause — an unmodelled "
            "aberration, a misassigned reflection, or underestimated peak uncertainties — "
            "should be found."
        )
    elif degrees >= MINIMUM_DEGREES and chi < CHI_SQUARED_LOW:
        warnings.append(
            f"Reduced χ² = {chi:.2f} is well below 1: the peak-position uncertainties are "
            "probably overstated. The quoted cell uncertainties have been scaled down with it."
        )
    if normalized is not None:
        outliers = [
            (label, float(value))
            for label, value in zip(labels, normalized, strict=False)
            if np.isfinite(value) and abs(value) > 3.0
        ]
        if outliers:
            listed = ", ".join(f"{label} ({value:+.1f}σ)" for label, value in outliers)
            warnings.append(
                f"Reflection{'s' if len(outliers) != 1 else ''} beyond ±3σ after the fit: "
                f"{listed}. Check each for misindexing or an overlapped peak, and consider "
                "excluding it."
            )
    strongest = strongest_correlation(result)
    if strongest is not None and abs(strongest[0]) > CORRELATION_LIMIT:
        value, first, second = strongest
        warnings.append(
            f"{first} and {second} are correlated at {value:+.3f}: over this angular range the "
            "scan can hardly tell them apart, which is why the cell uncertainty is larger than "
            "the peak scatter alone suggests. Reflections at higher angle break the "
            "correlation."
        )
    if result.method == "cohen" and result.extrapolation != "none":
        if result.drift_standard_uncertainty > 0.0:
            significance = abs(result.drift_coefficient) / result.drift_standard_uncertainty
            if significance < 2.0:
                warnings.append(
                    f"The systematic correction is not significant ({significance:.1f}σ). It "
                    "still costs precision; refitting without it may give a smaller "
                    "uncertainty if the instrument is known to be well aligned."
                )
    unconverged = [peak for peak in peaks_used if not bool(getattr(peak, "converged", True))]
    if unconverged:
        warnings.append(
            f"{len(unconverged)} of the peaks used did not converge in profile fitting; their "
            "positions should not be trusted."
        )
    poor = [peak for peak in peaks_used if float(getattr(peak, "reduced_chi_squared", 1.0)) > 10.0]
    if poor:
        angles = ", ".join(f"{float(peak.two_theta_deg):.2f}°" for peak in poor)
        warnings.append(
            f"The profile model fits {len(poor)} peak{'s' if len(poor) != 1 else ''} poorly "
            f"(peak-fit χ²ν > 10 at {angles}). Their positions may be biased by overlap, "
            "asymmetry or a curved background."
        )
    if figure_of_merit is not None and figure_of_merit < 10.0:
        warnings.append(
            f"The indexing figure of merit is low (M = {figure_of_merit:.1f}; above 10 is "
            "plausible, above 20 convincing). Check the phase and the radiation."
        )
    if unindexed_count:
        warnings.append(
            f"{unindexed_count} detected peak{'s were' if unindexed_count != 1 else ' was'} not "
            "indexed by this phase — a second phase, a Kβ or tungsten line, or the holder."
        )
    change = result.strain_relative_to_reference
    if change is not None and abs(change) > 0.01:
        warnings.append(
            f"The cell differs from the reference cell by {100.0 * change:+.2f} %. That is far "
            "more than any elastic strain; check the phase, the radiation and the wavelength."
        )
    return tuple(warnings)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def _two_theta_label() -> str:
    return "2θ (°)"


def scan_figure(
    two_theta_deg: np.ndarray,
    intensity: np.ndarray,
    *,
    unit: str,
    peak_angles: Sequence[float] = (),
    unindexed_angles: Sequence[float] = (),
    source: str,
) -> ResultFigure:
    """The measured scan exactly as read, with the detected peaks marked."""

    axis = np.asarray(two_theta_deg, dtype=float)
    counts = np.asarray(intensity, dtype=float)

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        axes.plot(axis, counts, lw=0.6, color=COLORS["data"], label="Measured")
        top = float(np.max(counts)) if counts.size else 1.0
        indexed = [angle for angle in peak_angles if angle not in set(unindexed_angles)]
        if indexed:
            heights = np.interp(indexed, axis, counts)
            axes.plot(
                indexed,
                heights + 0.04 * top,
                "v",
                ms=4,
                color=COLORS["accent"],
                label="Detected peak",
            )
        if unindexed_angles:
            heights = np.interp(list(unindexed_angles), axis, counts)
            axes.plot(
                list(unindexed_angles),
                heights + 0.04 * top,
                "v",
                ms=5,
                color=COLORS["unassigned"],
                label="Detected, not indexed",
            )
        axes.set_xlabel(_two_theta_label())
        axes.set_ylabel(f"Intensity ({unit})")
        axes.set_xlim(float(axis[0]), float(axis[-1]))
        axes.legend(loc="upper right", frameon=False)

    step = float(np.median(np.diff(axis))) if axis.size > 1 else float("nan")
    return render_figure(
        draw,
        key="scan",
        title="Measured XRD scan",
        caption=(
            f"{source}: {axis.size} points, {axis[0]:.2f}° to {axis[-1]:.2f}° 2θ, step "
            f"{step:.4f}°. Triangles mark the peaks detection found; red ones were not indexed."
        ),
        interpretation=(
            "This is the input exactly as read. Every peak used later appears here; a strong "
            "unmarked feature is a peak detection missed, and a red triangle is intensity the "
            "selected phase does not explain."
        ),
        height_in=3.2,
    )


def _peak_curves(
    axis: np.ndarray, peaks: Sequence[Any], doublet: tuple[float, float] | None
) -> tuple[np.ndarray, list[tuple[np.ndarray, np.ndarray]]]:
    """The fitted model on a fine grid per window, and on the scan axis inside windows only."""

    model = np.full_like(axis, np.nan)
    curves: list[tuple[np.ndarray, np.ndarray]] = []
    for peak in peaks:
        low, high = peak.window_deg
        inside = (axis >= low) & (axis <= high)
        if not np.any(inside):
            continue
        grid = np.linspace(low, high, 240)
        curves.append((grid, peak.evaluate(grid, doublet=doublet)))
        model[inside] = peak.evaluate(axis[inside], doublet=doublet)
    return model, curves


def peak_overview_figure(
    two_theta_deg: np.ndarray,
    intensity: np.ndarray,
    *,
    unit: str,
    peaks: Sequence[Any],
    radiation: Any,
    calculated_angles: Sequence[float],
    calculated_labels: Sequence[str],
) -> ResultFigure:
    """Every fitted profile over the scan, the difference inside each window, and the indexing."""

    axis = np.asarray(two_theta_deg, dtype=float)
    counts = np.asarray(intensity, dtype=float)
    doublet = kalpha_doublet_parameters(radiation)
    model_on_axis, curves = _peak_curves(axis, peaks, doublet)
    difference = counts - model_on_axis

    def draw(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[3.0, 1.0])
        top.plot(axis, counts, lw=0.6, color=COLORS["data"], label="Measured")
        for index, (grid, values) in enumerate(curves):
            top.plot(
                grid,
                values,
                lw=1.2,
                color=COLORS["model"],
                label="Fitted profile + local background" if index == 0 else None,
            )
        draw_ticks(
            top,
            calculated_angles,
            labels=calculated_labels,
            label="Calculated reflection (final cell)",
        )
        top.set_ylabel(f"Intensity ({unit})")
        # The reflection ticks get their own band below zero intensity, so they
        # and their labels never sit on top of the scan's baseline.
        ceiling = float(np.nanmax(counts)) if counts.size else 1.0
        top.set_ylim(bottom=-0.3 * ceiling, top=1.08 * ceiling)
        top.set_yticks([tick for tick in top.get_yticks() if 0.0 <= tick <= 1.08 * ceiling])
        top.legend(loc="upper right", frameon=False)
        bottom.plot(axis, difference, lw=0.6, color=COLORS["difference"])
        bottom.axhline(0.0, color=COLORS["guide"], lw=0.6)
        bottom.set_ylabel("Obs − fit")
        bottom.set_xlabel(_two_theta_label())
        bottom.set_xlim(float(axis[0]), float(axis[-1]))

    return render_figure(
        draw,
        key="fitted_peaks",
        title="Fitted peak profiles over the scan",
        caption=(
            "Top: measured scan (black), each fitted pseudo-Voigt profile with its K-α₂ partner "
            "and straight local background (orange), and the reflections of the final cell "
            "(ticks, labelled). Bottom: measured minus fitted, inside the fit windows only."
        ),
        interpretation=(
            "Every orange curve should sit on the black one. A difference with a sharp "
            "positive–negative swing at a peak means the fitted centre is off; a broad hump "
            "means an overlapped or asymmetric peak. Ticks without a peak are reflections "
            "predicted but not observed."
        ),
        height_in=4.0,
    )


def peak_windows_figure(
    two_theta_deg: np.ndarray,
    intensity: np.ndarray,
    *,
    peaks: Sequence[Any],
    labels: Sequence[str],
    radiation: Any,
    maximum_panels: int = 12,
) -> ResultFigure:
    """One small panel per fitted peak: the points, the fit, the background, the numbers."""

    axis = np.asarray(two_theta_deg, dtype=float)
    counts = np.asarray(intensity, dtype=float)
    doublet = kalpha_doublet_parameters(radiation)
    chosen = list(range(len(peaks)))
    if len(chosen) > maximum_panels:
        chosen = sorted({round(value) for value in np.linspace(0, len(peaks) - 1, maximum_panels)})
    columns = min(4, max(1, len(chosen)))
    rows = int(np.ceil(len(chosen) / columns))

    def draw(figure: Any) -> None:
        grid = figure.subplots(rows, columns, squeeze=False)
        for slot, axes in enumerate(grid.flat):
            if slot >= len(chosen):
                axes.set_visible(False)
                continue
            peak = peaks[chosen[slot]]
            low, high = peak.window_deg
            inside = (axis >= low) & (axis <= high)
            fine = np.linspace(low, high, 240)
            axes.plot(axis[inside], counts[inside], "o", ms=2.2, color=COLORS["data"])
            axes.plot(fine, peak.evaluate(fine, doublet=doublet), lw=1.1, color=COLORS["model"])
            background = peak.background_intercept + peak.background_slope * (
                fine - peak.two_theta_deg
            )
            axes.plot(fine, background, lw=0.7, ls="--", color=COLORS["background"])
            axes.axvline(peak.two_theta_deg, lw=0.6, color=COLORS["guide"])
            status = "" if peak.converged else " (not converged)"
            axes.set_title(labels[chosen[slot]] or f"{peak.two_theta_deg:.2f}°", fontsize=8)
            axes.text(
                0.02,
                0.97,
                f"2θ = {peak.two_theta_deg:.4f}°\n"
                f"σ = {1000.0 * peak.two_theta_standard_uncertainty_deg:.2f} m°\n"
                f"χ²ν = {peak.reduced_chi_squared:.2f}{status}",
                transform=axes.transAxes,
                va="top",
                fontsize=6.5,
                color=COLORS["warning"] if not peak.converged else COLORS["data"],
            )
            axes.tick_params(labelsize=6.5)
            axes.set_yticks([])
        figure.supxlabel(_two_theta_label(), fontsize=8)

    shown = len(chosen)
    return render_figure(
        draw,
        key="peak_windows",
        title="Each peak fit, close up",
        caption=(
            f"{shown} of {len(peaks)} fitted peaks"
            + (" (evenly spaced in angle)" if shown < len(peaks) else "")
            + ": measured points (black), fitted profile (orange), local background (dashed), "
            "fitted centre (grey line). σ is the standard uncertainty of the centre; χ²ν is the "
            "peak-profile reduced chi-squared of that window."
        ),
        interpretation=(
            "A good fit follows the points through the top and both tails. The peak-profile "
            "χ²ν judges this picture, not the cell: a value far above 1 flags a shape or "
            "background problem that can bias the centre, and its σ is then a lower bound."
        ),
        height_in=1.7 * rows + 0.5,
    )


def indexing_figure(
    *,
    observed_angles: Sequence[float],
    observed_heights: Sequence[float],
    observed_sigma_mdeg: Sequence[float],
    calculated_angles: Sequence[float],
    calculated_intensities: Sequence[float],
    delta_mdeg: Sequence[float],
    labels: Sequence[str],
    unindexed_angles: Sequence[float],
    unindexed_heights: Sequence[float],
    tolerance_deg: float,
) -> ResultFigure:
    """Observed against calculated lines, and the misfit of each assignment before correction."""

    heights = np.asarray(observed_heights, dtype=float)
    scale = float(np.max(heights)) if heights.size else 1.0
    calc = np.asarray(calculated_intensities, dtype=float)
    calc_scale = float(np.max(calc)) if calc.size else 1.0

    def draw(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[1.3, 1.0])
        top.vlines(
            observed_angles,
            0.0,
            heights / scale,
            color=COLORS["data"],
            lw=1.4,
            label="Observed peak (height)",
        )
        top.vlines(
            calculated_angles,
            0.0,
            -calc / calc_scale,
            color=COLORS["model"],
            lw=1.4,
            label="Calculated line (intensity)",
        )
        if len(unindexed_angles):
            top.vlines(
                unindexed_angles,
                0.0,
                np.asarray(unindexed_heights, dtype=float) / scale,
                color=COLORS["unassigned"],
                lw=1.6,
                label="Not indexed",
            )
        for angle, text in zip(calculated_angles, labels, strict=False):
            top.text(
                angle,
                0.03,
                text,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=COLORS["background"],
            )
        top.axhline(0.0, color=COLORS["guide"], lw=0.6)
        top.set_ylim(-1.15, 1.35)
        top.set_yticks([-1, 0, 1], labels=["1", "0", "1"])
        top.set_ylabel("calc.  |  obs.")
        top.legend(loc="upper right", frameon=False, ncols=3, fontsize=7)
        draw_residuals(
            bottom,
            observed_angles,
            delta_mdeg,
            observed_sigma_mdeg,
            xlabel=_two_theta_label(),
            ylabel="Δ2θ obs − calc (m°)",
        )
        limit = 1000.0 * tolerance_deg
        values = np.abs(np.asarray(delta_mdeg, dtype=float))
        if values.size and float(np.max(values)) > 0.3 * limit:
            bottom.axhline(limit, color=COLORS["warning"], lw=0.7, ls=":")
            bottom.axhline(-limit, color=COLORS["warning"], lw=0.7, ls=":")

    return render_figure(
        draw,
        key="indexing",
        title="Peak indexing",
        caption=(
            "Top: observed peaks (black, up) against the calculated lines of the phase (orange, "
            "down), labelled by reflection; red lines were not indexed. Bottom: observed minus "
            "calculated position for each assignment, against the cell of the final indexing "
            "pass and before any systematic correction, with ±1σ(2θ) error bars. Dotted red "
            f"lines, when shown, are the ±{tolerance_deg:.2f}° indexing tolerance."
        ),
        interpretation=(
            "A smooth trend of Δ2θ with angle is the signature of a zero or displacement error, "
            "which the lattice fit then removes. A single reflection far off the trend is a "
            "likely misassignment."
        ),
        height_in=3.8,
    )


def peak_quality_figure(
    *,
    angles: Sequence[float],
    sigma_mdeg: Sequence[float] | np.ndarray,
    fwhm_deg: Sequence[float] | np.ndarray,
    reduced_chi_squared: Sequence[float] | np.ndarray,
    converged: Sequence[bool],
    used: Sequence[bool],
) -> ResultFigure:
    """Per-peak fit diagnostics against angle: position uncertainty, width, profile χ²."""

    x = np.asarray(angles, dtype=float)
    mask_used = np.asarray(used, dtype=bool)
    mask_bad = ~np.asarray(converged, dtype=bool)

    def draw(figure: Any) -> None:
        first, second, third = figure.subplots(3, 1, sharex=True)
        for axes, values, label in (
            (first, np.asarray(sigma_mdeg, dtype=float), "σ(2θ) (m°)"),
            (second, np.asarray(fwhm_deg, dtype=float), "FWHM (°)"),
            (third, np.asarray(reduced_chi_squared, dtype=float), "Peak-fit χ²ν"),
        ):
            axes.plot(
                x[mask_used],
                values[mask_used],
                "o",
                ms=4,
                color=COLORS["data"],
                label="Used in the lattice fit",
            )
            if np.any(~mask_used):
                axes.plot(
                    x[~mask_used],
                    values[~mask_used],
                    "o",
                    ms=4,
                    mfc="none",
                    color=COLORS["background"],
                    label="Not used",
                )
            if np.any(mask_bad):
                axes.plot(
                    x[mask_bad],
                    values[mask_bad],
                    "x",
                    ms=6,
                    color=COLORS["warning"],
                    label="Did not converge",
                )
            axes.set_ylabel(label)
        third.axhline(1.0, color=COLORS["guide"], lw=0.8, ls="--")
        chi = np.asarray(reduced_chi_squared, dtype=float)
        if chi.size and float(np.nanmax(chi)) / max(float(np.nanmin(chi)), 1e-9) > 30.0:
            third.set_yscale("log")
        first.legend(loc="best", frameon=False, fontsize=7)
        third.set_xlabel(_two_theta_label())

    return render_figure(
        draw,
        key="peak_quality",
        title="Peak-fit diagnostics",
        caption=(
            "For every fitted peak: the standard uncertainty of its centre, its full width at "
            "half maximum, and the reduced χ² of its profile fit (dashed line at 1). Open "
            "circles were detected but not used in the lattice fit; crosses did not converge."
        ),
        interpretation=(
            "σ(2θ) sets each reflection's weight in the lattice fit, so a peak with a small σ "
            "dominates it. Widths should vary smoothly with angle; an outlier is an overlap or "
            "a second phase. These χ² values judge profile shapes, not the cell."
        ),
        height_in=4.6,
    )


def residual_figure(
    result: LatticeParameterResult,
    sigma_mdeg: np.ndarray,
    labels: Sequence[str],
) -> ResultFigure:
    """Observed minus calculated positions after the lattice fit, with ±1σ bars."""

    residual = 1000.0 * np.asarray(result.residual_two_theta_deg, dtype=float)

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        draw_residuals(
            axes,
            result.two_theta_deg,
            residual,
            sigma_mdeg,
            labels=labels,
            xlabel=_two_theta_label(),
            ylabel="2θobs − 2θcalc (m°)",
        )

    rms = float(np.sqrt(np.mean(np.square(residual)))) if residual.size else 0.0
    return render_figure(
        draw,
        key="residuals",
        title="Final residuals",
        caption=(
            "Observed position minus the position calculated from the fitted cell and "
            f"systematic correction, for the {result.reflection_count} reflections used, with "
            f"±1σ(2θ) error bars from the peak fits. RMS residual {rms:.2f} m°."
        ),
        interpretation=(
            "Points should scatter about zero with no trend, and most bars should cross zero. "
            "A trend with angle means the systematic correction has the wrong angular form."
        ),
        height_in=3.0,
    )


def normalized_residual_figure(
    result: LatticeParameterResult, normalized: np.ndarray, labels: Sequence[str]
) -> ResultFigure:
    """Residuals in units of their own uncertainty, against the ±2σ and ±3σ guides."""

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        draw_normalized_residuals(
            axes,
            result.two_theta_deg,
            normalized,
            labels=labels,
            xlabel=_two_theta_label(),
            ylabel="(2θobs − 2θcalc) / σ(2θ)",
        )

    finite = normalized[np.isfinite(normalized)]
    inside = int(np.count_nonzero(np.abs(finite) <= 2.0))
    beyond = int(np.count_nonzero(np.abs(finite) > 3.0))
    return render_figure(
        draw,
        key="normalized_residuals",
        title="Normalized residuals",
        caption=(
            "Each residual divided by its own standard uncertainty. Blue band ±2σ, yellow band "
            f"±3σ. {inside} of {finite.size} reflections lie within ±2σ; {beyond} lie beyond ±3σ "
            "(red diamonds). The sum of their squares divided by the degrees of freedom is the "
            f"lattice-fit reduced χ², {result.reduced_chi_squared:.2f}."
        ),
        interpretation=(
            "With correct uncertainties about 95 % of points fall inside ±2σ. Many points "
            "outside it means the peak uncertainties are too small or the model is incomplete; "
            "a single point beyond ±3σ is a candidate misassignment."
        ),
        height_in=3.0,
    )


def correction_figure(
    result: LatticeParameterResult,
    *,
    two_theta_range: tuple[float, float],
    sigma_mdeg: np.ndarray,
    labels: Sequence[str],
) -> ResultFigure:
    """The fitted angle-dependent systematic correction, as a curve with its ±1σ band."""

    low = max(float(two_theta_range[0]), 5.0)
    high = min(float(two_theta_range[1]), 170.0)
    grid = np.linspace(low, high, 400)
    curve, band = systematic_correction_curve(result, grid)
    at_reflections = 1000.0 * np.asarray(result.systematic_shift_deg, dtype=float)
    name = EXTRAPOLATION_NAMES[result.extrapolation]

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        axes.fill_between(
            grid,
            1000.0 * (curve - band),
            1000.0 * (curve + band),
            color=COLORS["band_3"],
            lw=0,
            label="±1σ of the correction",
        )
        axes.plot(
            grid, 1000.0 * curve, lw=1.4, color=COLORS["model"], label=f"Fitted correction ({name})"
        )
        axes.errorbar(
            result.two_theta_deg,
            at_reflections,
            yerr=sigma_mdeg,
            fmt="o",
            ms=4,
            capsize=2.5,
            color=COLORS["data"],
            ecolor=COLORS["background"],
            label="At each reflection, with that peak's σ(2θ)",
        )
        for x, y, text in zip(result.two_theta_deg, at_reflections, labels, strict=False):
            axes.annotate(
                text,
                (x, y),
                xytext=(5, 3),
                textcoords="offset points",
                ha="left",
                fontsize=6.5,
                color=COLORS["background"],
            )
        axes.axhline(0.0, color=COLORS["guide"], lw=0.8)
        axes.set_xlabel(_two_theta_label())
        axes.set_ylabel("Correction to 2θ (m°)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    significance = (
        abs(result.drift_coefficient) / result.drift_standard_uncertainty
        if result.drift_standard_uncertainty > 0.0
        else float("nan")
    )
    largest = float(np.max(np.abs(at_reflections))) if at_reflections.size else 0.0
    typical = float(np.median(sigma_mdeg)) if len(sigma_mdeg) else float("nan")
    verdict = (
        "much larger than the peak uncertainties, so the correction did real work and an "
        "uncorrected cell would be biased"
        if largest > 3.0 * typical
        else "comparable with the peak uncertainties, so the correction changes the cell little"
    )
    return render_figure(
        draw,
        key="systematic_correction",
        title="Systematic correction against 2θ",
        caption=(
            f"The angle-dependent term refined with the cell, D·sin²θ·f(θ) with f = {name}, "
            f"converted to degrees 2θ. D = {result.drift_coefficient:.3e} ± "
            f"{result.drift_standard_uncertainty:.1e} ({significance:.1f}σ). The band is ±1σ "
            "from σ(D); error bars are each reflection's own σ(2θ) for scale."
        ),
        interpretation=(
            f"The correction reaches {largest:.1f} m°, {verdict}. It falls to zero at "
            "2θ = 180°, which is why the fit extrapolates there. Its angular form is that of "
            "specimen displacement and absorption, but the fit cannot say which aberration "
            "caused it: read it as a correction, not as a measured displacement."
        ),
        height_in=3.1,
    )


def correlation_figure(result: LatticeParameterResult) -> ResultFigure | None:
    """Heat map of the refined-parameter correlation matrix, when there are two or more."""

    matrix = result.parameter_correlation
    if matrix is None or matrix.shape[0] < 2:
        return None
    names = result.correlation_parameter_names
    strongest = strongest_correlation(result)
    assert strongest is not None
    value, first, second = strongest

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        draw_correlation(figure, axes, matrix, names)

    size = 1.2 + 0.9 * len(names)
    return render_figure(
        draw,
        key="correlation",
        title="Parameter correlations",
        caption=(
            "Correlation coefficients of the refined least-squares parameters: the free "
            "components of the reciprocal metric tensor (A = a*², …) and, when refined, the "
            "systematic-correction coefficient D."
        ),
        interpretation=(
            f"The strongest correlation is {value:+.3f}, between {first} and {second}. "
            + (
                "Beyond ±0.95 the scan barely separates them, and each one's uncertainty is "
                "inflated accordingly; more high-angle reflections reduce it."
                if abs(value) > CORRELATION_LIMIT
                else "Below ±0.95 the parameters are adequately separated by this scan."
            )
        ),
        width_in=min(6.4, size + 1.6),
        height_in=size,
    )


def extrapolation_figure(
    *,
    abscissa: np.ndarray,
    per_reflection: np.ndarray,
    per_reflection_sigma: np.ndarray,
    labels: Sequence[str],
    slope: float,
    intercept: float,
    reported: float,
    reported_sigma: float,
    function_label: str,
) -> ResultFigure:
    """The classical extrapolation plot for a cubic cell, with error bars."""

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        axes.errorbar(
            abscissa,
            per_reflection,
            yerr=per_reflection_sigma,
            fmt="o",
            ms=4,
            capsize=2.5,
            color=COLORS["data"],
            ecolor=COLORS["background"],
            label="a from one reflection ±1σ",
        )
        top = float(np.max(abscissa)) if abscissa.size else 1.0
        grid = np.linspace(0.0, top * 1.05, 50)
        axes.plot(
            grid,
            intercept + slope * grid,
            lw=1.0,
            color=COLORS["model"],
            label="Unweighted straight-line guide",
        )
        axes.errorbar(
            [0.0],
            [reported],
            yerr=[reported_sigma],
            fmt="s",
            ms=5,
            capsize=3,
            color=COLORS["accent"],
            label="Reported a (weighted fit)",
        )
        for x, y, text in zip(abscissa, per_reflection, labels, strict=False):
            axes.annotate(
                text,
                (x, y),
                xytext=(5, 3),
                textcoords="offset points",
                ha="left",
                fontsize=6.5,
                color=COLORS["background"],
            )
        axes.set_xlabel(f"Extrapolation function f(θ) = {function_label}")
        axes.set_ylabel("a (Å)")
        axes.set_xlim(left=-0.03 * top)
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="extrapolation",
        title="Extrapolation to θ = 90°",
        caption=(
            "The cell edge computed from each reflection alone, against the extrapolation "
            "function, with ±1σ from that peak's σ(2θ) (σ(a)/a = cotθ·σ(θ)). The line is an "
            "unweighted guide; the square at f = 0 is the reported value from the weighted "
            "least-squares fit."
        ),
        interpretation=(
            "The slope is the systematic error; the scatter is the random error. Averaging the "
            "points would give their mean, not the intercept — the difference between the two "
            "is what the correction is worth."
        ),
        height_in=3.2,
    )


def cross_check_figure(rows: Sequence[dict[str, Any]]) -> ResultFigure:
    """The cell edge from each method on the same assignment, with ±1σ."""

    names = [str(row["method"]) for row in rows]
    values = np.array([float(row["a_angstrom"]) for row in rows])
    sigmas = np.array([float(row["sigma_angstrom"]) for row in rows])

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        positions = np.arange(len(rows))
        axes.errorbar(
            values,
            positions,
            xerr=sigmas,
            fmt="o",
            ms=5,
            capsize=3,
            color=COLORS["data"],
            ecolor=COLORS["background"],
        )
        axes.plot(values[:1], positions[:1], "s", ms=7, color=COLORS["accent"])
        axes.axvline(values[0], color=COLORS["accent"], lw=0.7, ls="--")
        axes.set_yticks(positions, labels=names)
        axes.set_ylim(len(rows) - 0.5, -0.5)
        axes.set_xlabel("a (Å)")
        axes.ticklabel_format(axis="x", useOffset=False)
        axes.xaxis.set_major_locator(MaxNLocator(nbins=4))

    return render_figure(
        draw,
        key="cross_check",
        title="Same peaks, different methods",
        caption=(
            "The cell edge ±1σ from each method applied to this run's peaks and assignment. "
            "The square and dashed line are the reported value."
        ),
        interpretation=(
            "Because the peaks and assignment are shared, every difference is due to the "
            "method alone. A large gap between fits with and without a systematic term means "
            "the correction mattered; error bars that do not overlap mean the methods disagree "
            "by more than their precision."
        ),
        height_in=0.6 + 0.45 * len(rows),
    )


def le_bail_figure(
    result: LatticeParameterResult,
    *,
    calculated_angles: Sequence[float],
    calculated_labels: Sequence[str],
) -> ResultFigure:
    """Observed, calculated and difference profiles of a whole-pattern fit."""

    assert result.profile_two_theta_deg is not None
    assert result.profile_observed is not None
    assert result.profile_calculated is not None
    axis = np.asarray(result.profile_two_theta_deg, dtype=float)
    observed = np.asarray(result.profile_observed, dtype=float)
    calculated = np.asarray(result.profile_calculated, dtype=float)

    def draw(figure: Any) -> None:
        top, bottom = figure.subplots(2, 1, sharex=True, height_ratios=[3.0, 1.0])
        draw_observed_model(
            top,
            axis,
            observed,
            calculated,
            difference_axes=bottom,
            observed_label="Observed (background removed)",
            model_label="Calculated (Le Bail)",
        )
        draw_ticks(top, calculated_angles, labels=calculated_labels, label="Reflection")
        top.set_ylabel("Intensity")
        bottom.set_xlabel(_two_theta_label())
        bottom.set_xlim(float(axis[0]), float(axis[-1]))

    r_wp = result.weighted_profile_r
    return render_figure(
        draw,
        key="le_bail_profile",
        title="Whole-pattern fit: observed, calculated, difference",
        caption=(
            "Background-subtracted measured profile (black), the Le Bail calculated profile "
            "(orange), reflection positions of the determined cell (ticks), and observed minus "
            f"calculated below. Reduced χ² {result.reduced_chi_squared:.2f}"
            + ("" if r_wp is None else f", R_wp {100.0 * r_wp:.2f} %")
            + "."
        ),
        interpretation=(
            "A Le Bail fit locates no individual peak, so there are no per-reflection "
            "residuals; the difference curve is the diagnostic. A flat difference means the "
            "cell and profile describe the scan; a derivative-shaped wiggle at every peak "
            "means the positions are off; an isolated peak in the difference is unmodelled "
            "intensity."
        ),
        height_in=4.0,
    )
