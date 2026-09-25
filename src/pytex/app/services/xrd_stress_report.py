# ruff: noqa: RUF001, RUF002
"""The report of a sin²ψ residual-stress evaluation: headline, warnings and figures.

Purpose
-------
:func:`pytex.diffraction.xrd_residual_stress.determine_residual_stress` does the
science. This module says what it found, in the order a reader needs it — the
stress and how far to trust it, the evidence (the raw peaks and the d against
sin²ψ lines), the diagnostics, the method — and draws the figures that let the
reader check each step. It computes no new science: every number is read from
the result objects or is an algebraic rearrangement of them (a residual divided
by its own uncertainty; the fitted tensor evaluated along an azimuth), and the
tests hold those rearrangements to the quantities they rearrange.

Vocabulary kept distinct on purpose
-----------------------------------
- **Peak-fit χ²ν** judges one profile against the counts of one scan.
  **Strain-fit χ²ν** judges the whole set of spacings against the stress
  tensor. They answer different questions and are never reported under one name.
- **Statistical** uncertainty comes from the peak positions. **d₀** and the
  **elastic constants** are systematic: repeating the measurement does not
  reduce them, and they are listed separately in the budget for that reason.
- The stress is the **macroscopic (type I) stress** averaged over the depth
  the X-rays reach, for an untextured material. It is not the stress at the
  surface, and not a grain-scale stress.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, Literal

import numpy as np

from pytex.app.figures import (
    COLORS,
    draw_correlation,
    draw_normalized_residuals,
    render_figure,
)
from pytex.app.results import ResultFigure, ResultMetric
from pytex.core.symbols import symbol_latex, symbol_text
from pytex.diffraction.xrd_peaks import kalpha_doublet_parameters
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants,
    ResidualStressResult,
    Sin2PsiMeasurement,
    StressTensorFit,
    lpa_factor,
    strain_design_matrix,
)
from pytex.properties.tensors import StiffnessTensor

__all__ = [
    "COMPONENT_LABELS",
    "budget_figure",
    "dec_figure",
    "geometry_figure",
    "mohr_figure",
    "normalized_strain_residuals",
    "peak_fits_figure",
    "peak_quality_figure",
    "peak_shift_figure",
    "sigma_phi_figure",
    "sin2psi_figure",
    "splitting_figure",
    "strain_figure",
    "strain_residual_figure",
    "stress_correlation_figure",
    "stress_highlights",
    "stress_warnings",
    "tensor_sigma_phi",
]

#: Registered symbol of each stress component, keyed as the library names it.
_COMPONENT_SYMBOLS = {
    "sigma_11": "sigma_11",
    "sigma_22": "sigma_22",
    "sigma_12": "sigma_12",
    "sigma_13": "sigma_13",
    "sigma_23": "sigma_23",
    "sigma_33": "sigma_33",
}

#: Display text of each component, for tables and prose.
COMPONENT_LABELS = {name: symbol_text(key) for name, key in _COMPONENT_SYMBOLS.items()}

#: One colour per azimuth, in the order the azimuths are drawn.
_AZIMUTH_COLORS = ("#1d4ed8", "#b45309", "#0f766e", "#7c3aed", "#be123c", "#4d7c0f")

#: A strain-fit reduced chi-squared outside this band gets a warning.
CHI_SQUARED_HIGH = 3.0
CHI_SQUARED_LOW = 0.3
#: |t| of the sin⁴ψ term above which d is called non-linear in sin²ψ.
CURVATURE_T_LIMIT = 3.0
#: A correlation this close to ±1 means two components are barely separable.
CORRELATION_LIMIT = 0.95


def _psi() -> str:
    return symbol_text("stress_tilt")


def _phi() -> str:
    return symbol_text("stress_azimuth")


def _math(name: str) -> str:
    return f"${symbol_latex(name)}$"


def _sin2psi_label() -> str:
    return f"$\\sin^{{2}}{symbol_latex('stress_tilt')}$"


def _azimuth_color(index: int) -> str:
    return _AZIMUTH_COLORS[index % len(_AZIMUTH_COLORS)]


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------


def normalized_strain_residuals(result: ResidualStressResult) -> np.ndarray:
    """Return ``(ε_obs − ε_fit) / u(ε)`` for every measurement of a tensor fit.

    The sum of squares of these values over the degrees of freedom *is* the
    reported strain-fit χ²ν, because the fit was weighted by ``1/u(ε)²``.
    Empty when no tensor was determined.
    """

    if result.tensor is None:
        return np.zeros(0)
    values: np.ndarray = result.tensor.residual_strain / result.strain_uncertainty
    return values


def tensor_sigma_phi(
    fit: StressTensorFit, phi_deg: Sequence[float] | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``σ_φ`` of the fitted tensor along each azimuth and its uncertainty.

    ``σ_φ = σ11 cos²φ + σ12 sin2φ + σ22 sin²φ`` is linear in the components, so
    its variance is ``g Σ gᵀ`` with ``g`` the gradient and ``Σ`` the combined
    covariance: the band a reader should see around the curve.
    """

    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))
    gradient = np.zeros((phi.size, len(fit.component_names)))
    names = fit.component_names
    gradient[:, names.index("sigma_11")] = np.cos(phi) ** 2
    gradient[:, names.index("sigma_22")] = np.sin(phi) ** 2
    gradient[:, names.index("sigma_12")] = np.sin(2.0 * phi)
    values = gradient @ fit.values_mpa
    variance = np.einsum("ni,ij,nj->n", gradient, fit.combined_covariance_mpa2, gradient)
    return values, np.sqrt(np.maximum(variance, 0.0))


# ---------------------------------------------------------------------------
# Headline and warnings
# ---------------------------------------------------------------------------


def stress_highlights(result: ResidualStressResult) -> tuple[ResultMetric, ...]:
    """The numbers a reader needs first: the stress, its reliability, its basis."""

    metrics: list[ResultMetric] = []
    fit = result.tensor
    if fit is not None:
        for name, value, u in zip(
            fit.component_names, fit.values_mpa, fit.combined_uncertainty_mpa, strict=True
        ):
            metrics.append(
                ResultMetric(
                    COMPONENT_LABELS[name],
                    f"{value:.1f} ± {u:.1f}",
                    units="MPa",
                    help_text=(
                        "Value ± combined standard uncertainty (statistical, d₀ and elastic "
                        "constants in quadrature). Tensile positive."
                    ),
                )
            )
        principal = fit.in_plane_principal()
        metrics.append(
            ResultMetric(
                "In-plane principal stresses σI, σII",
                f"{principal['sigma_I_mpa']:.1f} ± {principal['sigma_I_uncertainty_mpa']:.1f}, "
                f"{principal['sigma_II_mpa']:.1f} ± {principal['sigma_II_uncertainty_mpa']:.1f}",
                units="MPa",
            )
        )
        metrics.append(
            ResultMetric(
                f"Direction of σI from S1 (azimuth {_phi()})",
                f"{principal['angle_deg']:.1f} ± {principal['angle_uncertainty_deg']:.1f}",
                units="°",
                help_text="Undefined, and its uncertainty large, when σI ≈ σII.",
            )
        )
        equivalent, equivalent_u = fit.von_mises()
        metrics.append(
            ResultMetric(
                "von Mises equivalent stress", f"{equivalent:.1f} ± {equivalent_u:.1f}", units="MPa"
            )
        )
        metrics.append(
            ResultMetric(
                "Strain-fit χ²ν",
                round(fit.reduced_chi_squared, 3),
                help_text=(
                    "Near 1 when the stress model and the peak-position uncertainties agree. "
                    "Above 3: scatter the model does not explain (texture, gradients, a poor "
                    "peak fit). The statistical uncertainty is already scaled by √χ²ν when it "
                    "exceeds 1."
                ),
            )
        )
        metrics.append(ResultMetric("Degrees of freedom", fit.degrees_of_freedom))
    else:
        for line in result.regressions:
            metrics.append(
                ResultMetric(
                    f"σ{_phi()} at {_phi()} = {line.phi_deg:g}°",
                    f"{line.sigma_phi_mpa:.1f} ± {line.sigma_phi_uncertainty_mpa:.1f}",
                    units="MPa",
                )
            )
    metrics.append(
        ResultMetric(
            "Stress-free spacing d₀" + (" (refined)" if result.d0_refined else ""),
            f"{result.d0_angstrom:.6f} ± {result.d0_uncertainty_angstrom:.6f}",
            units="Å",
        )
    )
    metrics.append(
        ResultMetric(
            f"Elastic constants {symbol_text('dec_s1')}, {symbol_text('dec_half_s2')}",
            f"{result.dec.s1_per_tpa:.3f}, {result.dec.half_s2_per_tpa:.3f}",
            units="TPa⁻¹",
        )
    )
    used = int(np.count_nonzero(result.included_mask))
    metrics.append(
        ResultMetric(
            "Measurements used",
            f"{used} of {len(result.peaks)}",
            help_text=(
                "Measurements left out by the analyst (Excluded measurements) enter no fit but "
                "stay in the plots, as red crosses, and in the tables."
            ),
        )
    )
    return tuple(metrics)


def stress_warnings(
    result: ResidualStressResult,
    *,
    measurement: Sin2PsiMeasurement | None,
    stress_state: str,
) -> tuple[str, ...]:
    """Reasons, found by the evaluation itself, to distrust this result."""

    warnings: list[str] = []
    fit = result.tensor
    if fit is None:
        warnings.append(
            "No stress tensor: " + (result.tensor_unavailable_reason or "it is not determined.")
        )
    else:
        if fit.reduced_chi_squared > CHI_SQUARED_HIGH:
            warnings.append(
                f"The strain-fit χ²ν is {fit.reduced_chi_squared:.1f}: the spacings scatter "
                "about the stress model by more than the peak-position uncertainties allow. "
                "The statistical uncertainty has been widened by √χ²ν, but look at the "
                "residuals for structure — curvature or oscillation of d against sin²ψ means "
                "texture or a stress gradient, which the model does not describe."
            )
        elif fit.reduced_chi_squared < CHI_SQUARED_LOW and fit.degrees_of_freedom >= 3:
            warnings.append(
                f"The strain-fit χ²ν is {fit.reduced_chi_squared:.2f}: the peak-position "
                "uncertainties look overstated, so the statistical uncertainty is conservative."
            )
        # Whether the measured directions separate two components is a property
        # of the design, so it is judged on the statistical covariance. The
        # combined covariance also carries d0 and the elastic constants, which
        # shift components together by construction and would make every pair
        # look inseparable whenever u(d0) dominates.
        internal = fit.covariance_internal_mpa2
        scale = np.sqrt(np.diag(internal))
        correlation = internal / np.outer(scale, scale)
        names = fit.component_names
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                if abs(correlation[i, j]) > CORRELATION_LIMIT:
                    warnings.append(
                        f"{COMPONENT_LABELS[names[i]]} and {COMPONENT_LABELS[names[j]]} are "
                        f"correlated at {correlation[i, j]:+.3f} by the measurement itself: the "
                        "measured directions barely separate them, and only their combination is "
                        "well determined. More azimuths, or tilts further out, separate them."
                    )
        budget = fit.budget_mpa
        if "d0" in budget and np.any(budget["d0"] > 2.0 * budget["statistical"]):
            warnings.append(
                "The uncertainty of d₀ dominates the budget. Through S₁ it moves σ11 + σ22 "
                "together; the slopes — the per-azimuth stresses — are nearly immune to it. "
                "Measure d₀ on stress-free material (powder, or a cut-free coupon) of the same "
                "composition, or refine it under the plane-stress assumption."
            )
    for line in result.regressions:
        if math.isfinite(line.curvature_t) and abs(line.curvature_t) > CURVATURE_T_LIMIT:
            warnings.append(
                f"At {_phi()} = {line.phi_deg:g}° d is curved in sin²ψ (the sin⁴ψ term is "
                f"{line.curvature_t:.1f} times its uncertainty). The sin²ψ method assumes a "
                "straight line; curvature points to a stress gradient within the penetration "
                "depth or, if it oscillates, to texture."
            )
        if (
            stress_state == "biaxial"
            and line.has_splitting_term
            and abs(line.tau_phi_mpa) > 3.0 * line.tau_phi_uncertainty_mpa
            and abs(line.tau_phi_mpa) > 10.0
        ):
            warnings.append(
                f"At {_phi()} = {line.phi_deg:g}° the ψ > 0 and ψ < 0 branches split: "
                f"τ{_phi()} = {line.tau_phi_mpa:.0f} ± {line.tau_phi_uncertainty_mpa:.0f} MPa. "
                "A biaxial evaluation assumes no out-of-plane shear; choose the evaluation with "
                "shear components."
            )
    tilts = np.abs(result.psi_deg)
    top = float(np.max(np.sin(np.deg2rad(tilts)) ** 2)) if tilts.size else 0.0
    if top < 0.3:
        warnings.append(
            f"The largest tilt reaches sin²ψ = {top:.2f}. The slope is determined over a short "
            "lever arm; tilts to sin²ψ ≈ 0.5 (ψ ≈ 45°) or beyond reduce its uncertainty."
        )
    two_theta = result.two_theta_deg
    if two_theta.size and float(np.median(two_theta)) < 120.0:
        warnings.append(
            f"The reflection lies at 2θ ≈ {float(np.median(two_theta)):.0f}°. Strain "
            "sensitivity goes as tanθ, so a reflection above about 140° (Cr Kα on ferrite "
            "(211), Cu Kα on nickel (420)) gives several times better precision. At a "
            "synchrotron's short wavelength every reflection is at low angle; there the remedy "
            "is counting statistics, a higher-index reflection, and χ-tilting, since ω-tilting "
            "soon takes the beam below the surface."
        )
    failed = [peak for peak in result.peaks if not peak.converged]
    if failed:
        warnings.append(
            f"{len(failed)} peak location(s) did not converge cleanly "
            f"({', '.join(f'{_phi()} {p.phi_deg:g}°, {_psi()} {p.psi_deg:g}°' for p in failed[:4])}"
            f"{'…' if len(failed) > 4 else ''}); check them in the peak-fit figure."
        )
    if measurement is not None and measurement.synthetic:
        warnings.append(
            "The scans were generated, not measured: a demonstration with a known answer."
        )
    return tuple(warnings)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def geometry_figure(result: ResidualStressResult) -> ResultFigure:
    """The measured directions on the specimen hemisphere."""

    phi = np.deg2rad(result.phi_deg)
    psi = result.psi_deg

    def draw(figure: Any) -> None:
        axes = figure.add_subplot(projection="polar")
        # A negative tilt at azimuth φ is the direction at φ + 180° and |ψ|.
        azimuth = np.where(psi < 0.0, phi + np.pi, phi)
        radius = np.abs(psi)
        for index, value in enumerate(sorted(set(np.round(result.phi_deg % 360.0, 6)))):
            members = np.isclose(np.round(result.phi_deg % 360.0, 6), value)
            color = _azimuth_color(index)
            positive = members & (psi >= 0.0)
            negative = members & (psi < 0.0)
            axes.plot(
                azimuth[positive],
                radius[positive],
                "o",
                ms=5,
                color=color,
                label=f"{_phi()} = {value:g}°",
            )
            axes.plot(azimuth[negative], radius[negative], "o", ms=5, mfc="white", color=color)
        axes.set_theta_zero_location("E")
        axes.set_theta_direction(1)
        axes.set_rmax(max(60.0, float(np.max(np.abs(psi))) + 5.0))
        axes.set_rlabel_position(22.5)
        axes.set_xticks(np.deg2rad([0, 90, 180, 270]), labels=["S1", "S2", "−S1", "−S2"])
        axes.grid(True, lw=0.4, color=COLORS["guide"])
        axes.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="measurement_directions",
        title="Measured directions",
        caption=(
            f"Each scattering vector as a point on the specimen hemisphere seen from above: "
            f"azimuth {_phi()} from S1 towards S2, radius the tilt |{_psi()}| in degrees from the "
            f"surface normal. Filled markers are {_psi()} > 0; open markers {_psi()} < 0, drawn "
            f"at {_phi()} + 180°, which is where that direction actually points."
        ),
        interpretation=(
            "Three azimuths not 180° apart determine the in-plane tensor. Points on both sides "
            "of the centre along one azimuth are what exposes out-of-plane shear."
        ),
        width_in=5.2,
        height_in=3.6,
    )


def peak_fits_figure(
    measurement: Sin2PsiMeasurement, result: ResidualStressResult, *, maximum_rows: int = 11
) -> ResultFigure:
    """Every raw scan with its located peak: one column per azimuth, one row per tilt."""

    doublet = kalpha_doublet_parameters(measurement.radiation)
    azimuths = sorted({round(scan.phi_deg % 360.0, 6) for scan in measurement.scans})
    tilts = sorted({round(scan.psi_deg, 6) for scan in measurement.scans})
    if len(tilts) > maximum_rows:
        tilts = [tilts[round(i)] for i in np.linspace(0, len(tilts) - 1, maximum_rows)]
    peaks = {(round(p.phi_deg % 360.0, 6), round(p.psi_deg, 6)): p for p in result.peaks}
    scans = {(round(s.phi_deg % 360.0, 6), round(s.psi_deg, 6)): s for s in measurement.scans}
    excluded = {
        (round(result.peaks[i].phi_deg % 360.0, 6), round(result.peaks[i].psi_deg, 6))
        for i in result.excluded_indices
    }
    rows, columns = len(tilts), len(azimuths)

    def draw(figure: Any) -> None:
        grid = figure.subplots(rows, columns, squeeze=False, sharex="col")
        for row, psi in enumerate(tilts):
            for column, phi in enumerate(azimuths):
                axes = grid[row, column]
                scan = scans.get((phi, psi))
                peak = peaks.get((phi, psi))
                if scan is None or peak is None:
                    axes.set_visible(False)
                    continue
                axis = np.asarray(scan.pattern.two_theta_deg)
                counts = np.asarray(scan.pattern.intensity, dtype=float)
                if peak.lpa_corrected:
                    # Shown on the scale the peak was located on: divided by the
                    # LPA factor normalized at the same reference angle.
                    fraction = measurement.radiation.polarization_perpendicular_fraction
                    factor = lpa_factor(
                        axis,
                        psi_deg=peak.psi_deg,
                        geometry=measurement.geometry,
                        perpendicular_fraction=fraction,
                    )
                    reference = lpa_factor(
                        [peak.lpa_reference_deg],
                        psi_deg=peak.psi_deg,
                        geometry=measurement.geometry,
                        perpendicular_fraction=fraction,
                    )[0]
                    counts = counts / (factor / reference)
                axes.plot(axis, counts, ".", ms=1.6, color=COLORS["data"])
                fit = peak.peak_fit
                if fit is not None:
                    fine = np.linspace(fit.window_deg[0], fit.window_deg[1], 400)
                    axes.plot(
                        fine, fit.evaluate(fine, doublet=doublet), lw=1.0, color=COLORS["model"]
                    )
                axes.axvline(peak.two_theta_deg, lw=0.9, ls="--", color=COLORS["difference"])
                axes.set_yticks([])
                axes.tick_params(labelsize=6)
                axes.text(
                    0.02,
                    0.95,
                    f"{_psi()} = {peak.psi_deg:g}°\n2θ = {peak.two_theta_deg:.3f}°\n"
                    f"u = {1000.0 * peak.two_theta_uncertainty_deg:.1f} m°"
                    + ("\nexcluded" if (phi, psi) in excluded else ""),
                    transform=axes.transAxes,
                    va="top",
                    fontsize=5.5,
                    color=(
                        COLORS["data"]
                        if peak.converged and (phi, psi) not in excluded
                        else COLORS["warning"]
                    ),
                )
                if row == 0:
                    axes.set_title(f"{_phi()} = {phi:g}°", fontsize=8)
        figure.supxlabel("2θ (°)", fontsize=8)

    return render_figure(
        draw,
        key="peak_fits",
        title="Every scan and its located peak",
        caption=(
            f"The measured counts of each scan (points), one column per azimuth {_phi()} and "
            f"one row per tilt {_psi()}, with the fitted profile (orange; for the profile-fit "
            "method, Kα1 and Kα2 with a straight background) and the located Kα1 position "
            "(blue dashed line). u is the standard uncertainty of that position. When the LPA "
            "correction was applied the counts are shown divided by it, on the scale the peak "
            "was located on."
        ),
        interpretation=(
            "The located position should sit on the peak in every panel, and the peak should "
            "move steadily from row to row. A panel where the curve misses the points, or the "
            "line misses the peak, is a bad location that will show up as an outlier below."
        ),
        width_in=6.4,
        height_in=0.95 * rows + 0.6,
    )


def peak_shift_figure(
    measurement: Sin2PsiMeasurement, result: ResidualStressResult
) -> ResultFigure:
    """The normalized profiles of each azimuth overlaid, coloured by tilt."""

    azimuths = sorted({round(scan.phi_deg % 360.0, 6) for scan in measurement.scans})
    peaks = {(round(p.phi_deg % 360.0, 6), round(p.psi_deg, 6)): p for p in result.peaks}

    def draw(figure: Any) -> None:
        import matplotlib
        from matplotlib import cm
        from matplotlib.colors import Normalize

        grid = figure.subplots(1, len(azimuths), squeeze=False, sharey=True)
        tilts = np.array([scan.psi_deg for scan in measurement.scans])
        scale = Normalize(vmin=float(tilts.min()), vmax=float(tilts.max()))
        colormap = matplotlib.colormaps["coolwarm"]
        for column, phi in enumerate(azimuths):
            axes = grid[0, column]
            for scan in measurement.scans:
                if round(scan.phi_deg % 360.0, 6) != phi:
                    continue
                axis = np.asarray(scan.pattern.two_theta_deg)
                counts = np.asarray(scan.pattern.intensity, dtype=float)
                base = float(np.percentile(counts, 5))
                top = float(np.max(counts)) - base
                color = colormap(scale(scan.psi_deg))
                axes.plot(axis, (counts - base) / (top if top > 0 else 1.0), lw=0.7, color=color)
                peak = peaks.get((phi, round(scan.psi_deg, 6)))
                if peak is not None:
                    axes.axvline(peak.two_theta_deg, lw=0.5, color=color, alpha=0.8)
            axes.set_title(f"{_phi()} = {phi:g}°", fontsize=8)
            axes.set_xlabel("2θ (°)")
        grid[0, 0].set_ylabel("Normalized counts")
        mappable = cm.ScalarMappable(norm=scale, cmap=colormap)
        colorbar = figure.colorbar(mappable, ax=list(grid.flat), shrink=0.9)
        colorbar.set_label(f"{_psi()} (°)")

    return render_figure(
        draw,
        key="peak_shift",
        title="The peak moves with tilt",
        caption=(
            f"The scans of each azimuth, background-offset and scaled to unit height, coloured "
            f"by tilt {_psi()}, with the located positions as thin lines of the same colour."
        ),
        interpretation=(
            "This is the raw signal of the method. A compressive in-plane stress contracts the "
            "planes inclined to the surface, so the peak moves to higher 2θ as the tilt grows; "
            "a tensile stress moves it the other way. The width growing with tilt is "
            "defocusing, not stress."
        ),
        height_in=2.8,
    )


def sin2psi_figure(result: ResidualStressResult) -> ResultFigure:
    """The d against sin²ψ lines: the picture the method is named for."""

    lines = result.regressions
    columns = min(3, max(1, len(lines)))
    rows = math.ceil(len(lines) / columns) if lines else 1
    fit = result.tensor

    def draw(figure: Any) -> None:
        grid = figure.subplots(rows, columns, squeeze=False)
        handles: list[Any] = []
        names: list[str] = []
        for slot, axes in enumerate(grid.flat):
            if slot >= len(lines):
                axes.set_visible(False)
                continue
            line = lines[slot]
            color = _azimuth_color(slot)
            positive = line.psi_deg >= 0.0
            axes.errorbar(
                line.sin2psi[positive],
                line.d_angstrom[positive],
                yerr=line.d_uncertainty_angstrom[positive],
                fmt="o",
                ms=4,
                color=color,
                capsize=2,
                label=f"{_psi()} ≥ 0",
            )
            if np.any(~positive):
                axes.errorbar(
                    line.sin2psi[~positive],
                    line.d_angstrom[~positive],
                    yerr=line.d_uncertainty_angstrom[~positive],
                    fmt="o",
                    ms=4,
                    mfc="white",
                    color=color,
                    capsize=2,
                    label=f"{_psi()} < 0",
                )
            left_out = _excluded_at(result, line.phi_deg)
            if left_out.size:
                axes.errorbar(
                    result.sin2psi[left_out],
                    result.d_angstrom[left_out],
                    yerr=result.d_uncertainty_angstrom[left_out],
                    fmt="x",
                    ms=6,
                    mew=1.4,
                    color=COLORS["warning"],
                    capsize=2,
                    label="Excluded by the analyst",
                )
            top = float(np.max(line.sin2psi))
            psi_grid = np.rad2deg(np.arcsin(np.sqrt(np.linspace(0.0, top, 60))))
            grid_x = np.sin(np.deg2rad(psi_grid)) ** 2
            axes.plot(
                grid_x,
                line.intercept_angstrom + line.slope_angstrom * grid_x,
                lw=1.1,
                color=COLORS["model"],
                label="Line fit at this azimuth",
            )
            if line.has_splitting_term:
                for sign in (1.0, -1.0):
                    axes.plot(
                        grid_x,
                        line.fitted_d(sign * psi_grid),
                        lw=0.7,
                        ls=":",
                        color=COLORS["model"],
                    )
            if fit is not None:
                # With out-of-plane shear the tensor predicts two branches; they
                # coincide without it, so drawing both costs nothing then.
                signs = (1.0, -1.0) if np.any(line.psi_deg < 0.0) else (1.0,)
                for sign in signs:
                    model = _tensor_d_curve(result, line.phi_deg, sign * psi_grid)
                    axes.plot(
                        grid_x,
                        model,
                        lw=0.9,
                        ls="--" if sign > 0 else "-.",
                        color=COLORS["accent"],
                        label=(
                            f"Stress tensor, {_psi()} > 0"
                            if sign > 0
                            else f"Stress tensor, {_psi()} < 0"
                        )
                        if len(signs) > 1
                        else "Stress tensor",
                    )
            axes.set_title(f"{_phi()} = {line.phi_deg:g}°", fontsize=8)
            axes.text(
                0.03,
                0.04,
                f"slope = {line.slope_angstrom * 1e3:.4f} mÅ\n"
                f"σ{_phi()} = {line.sigma_phi_mpa:.0f} ± {line.sigma_phi_uncertainty_mpa:.0f} MPa",
                transform=axes.transAxes,
                fontsize=6.5,
                va="bottom",
            )
            axes.set_xlabel(_sin2psi_label())
            axes.ticklabel_format(axis="y", useOffset=False)
            axes.tick_params(labelsize=7)
            if slot % columns == 0:
                axes.set_ylabel("d (Å)")
            for handle, name in zip(*axes.get_legend_handles_labels(), strict=True):
                if name not in names:
                    handles.append(handle)
                    names.append(name)
        figure.legend(
            handles,
            names,
            loc="outside lower center",
            ncol=min(5, len(names)),
            frameon=False,
            fontsize=6.5,
        )

    return render_figure(
        draw,
        key="d_vs_sin2psi",
        title="d against sin²ψ — the slope is the stress",
        caption=(
            f"The spacing d = λ/(2 sinθ) at every tilt, ±u(d), one panel per azimuth; filled "
            f"markers {_psi()} ≥ 0, open {_psi()} < 0. Orange: the weighted straight line at "
            "that azimuth (dotted: its two branches when a ψ-splitting term was fitted). Green "
            "dashed: what the jointly fitted stress tensor predicts."
        ),
        interpretation=(
            "The slope is d₀·½S₂·σφ, so a falling line is compressive and a rising one tensile. "
            "The points should scatter about a straight line within their error bars; the two "
            "branches coincide unless there is out-of-plane shear; and the green curve should "
            "follow the orange one at every azimuth if one stress tensor explains all of them."
        ),
        height_in=2.5 * rows + 0.8,
    )


def _excluded_at(result: ResidualStressResult, phi_deg: float) -> np.ndarray:
    """Indices of the excluded measurements taken at one azimuth."""

    azimuth = round(phi_deg % 360.0, 6)
    return np.array(
        [
            index
            for index in result.excluded_indices
            if round(result.peaks[index].phi_deg % 360.0, 6) == azimuth
        ],
        dtype=int,
    )


def _tensor_d_curve(
    result: ResidualStressResult, phi_deg: float, psi_deg: np.ndarray
) -> np.ndarray:
    fit = result.tensor
    assert fit is not None
    design = strain_design_matrix(
        np.full_like(psi_deg, phi_deg),
        psi_deg,
        s1_per_tpa=result.dec.s1_per_tpa,
        half_s2_per_tpa=result.dec.half_s2_per_tpa,
        components=fit.component_names,
    )
    curve: np.ndarray = result.d0_angstrom * (1.0 + design @ fit.values_mpa)
    return curve


def strain_figure(result: ResidualStressResult) -> ResultFigure:
    """Strain against sin²ψ for every azimuth on one axis, with the tensor model."""

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        azimuths = [line.phi_deg for line in result.regressions]
        phi = np.round(result.phi_deg % 360.0, 6)
        for index, value in enumerate(azimuths):
            members = np.isclose(phi, value)
            color = _azimuth_color(index)
            left_out = members & ~result.included_mask
            members = members & result.included_mask
            if np.any(left_out):
                axes.plot(
                    result.sin2psi[left_out],
                    1e6 * result.strain[left_out],
                    "x",
                    ms=6,
                    mew=1.4,
                    color=COLORS["warning"],
                )
            axes.errorbar(
                result.sin2psi[members],
                1e6 * result.strain[members],
                yerr=1e6 * result.strain_uncertainty[members],
                fmt="o",
                ms=3.5,
                color=color,
                capsize=1.5,
                label=f"{_phi()} = {value:g}°",
            )
            if result.tensor is not None:
                top = float(np.max(result.sin2psi[members]))
                psi_grid = np.rad2deg(np.arcsin(np.sqrt(np.linspace(0.0, top, 60))))
                model = _tensor_d_curve(result, value, psi_grid) / result.d0_angstrom - 1.0
                axes.plot(np.sin(np.deg2rad(psi_grid)) ** 2, 1e6 * model, lw=1.0, color=color)
        axes.axhline(0.0, lw=0.6, color=COLORS["guide"])
        axes.set_xlabel(_sin2psi_label())
        axes.set_ylabel("Lattice strain (d − d₀)/d₀ (10⁻⁶)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="strain_vs_sin2psi",
        title="Lattice strain against sin²ψ",
        caption=(
            "Strain ε = (d − d₀)/d₀ in microstrain, ±u(ε), all azimuths together; lines are "
            "the stress tensor's prediction for each azimuth."
        ),
        interpretation=(
            "Where the lines cross zero strain lies the 'strain-free direction'. Unlike the "
            "slopes, the absolute level of every point depends on d₀: an error in d₀ shifts all "
            "of them together, and through S₁ it moves σ11 + σ22 — the reason d₀ appears in "
            "the uncertainty budget."
        ),
        height_in=3.2,
    )


def splitting_figure(result: ResidualStressResult) -> ResultFigure | None:
    """The ψ-splitting evaluation: branch mean and half-difference."""

    data = [(line, *line.branch_averages()) for line in result.regressions]
    data = [row for row in data if row[1].size]
    if not data:
        return None

    def draw(figure: Any) -> None:
        left, right = figure.subplots(1, 2)
        for index, (line, s2, a1, u1, a2, u2) in enumerate(data):
            color = _azimuth_color(index)
            left.errorbar(
                s2,
                a1,
                yerr=u1,
                fmt="o",
                ms=3.5,
                color=color,
                capsize=1.5,
                label=f"{_phi()} = {line.phi_deg:g}°",
            )
            abscissa = np.sin(2.0 * np.arcsin(np.sqrt(s2)))
            right.errorbar(
                abscissa, 1e3 * a2, yerr=1e3 * u2, fmt="o", ms=3.5, color=color, capsize=1.5
            )
            if line.has_splitting_term:
                grid = np.linspace(0.0, 1.0, 20)
                right.plot(grid, 1e3 * line.splitting_angstrom * grid, lw=0.9, color=color)
        left.set_xlabel(_sin2psi_label())
        left.set_ylabel("a₁ = (d₊ + d₋)/2 (Å)")
        left.ticklabel_format(axis="y", useOffset=False)
        right.axhline(0.0, lw=0.6, color=COLORS["guide"])
        right.set_xlabel(f"$\\sin|2{symbol_latex('stress_tilt')}|$")
        right.set_ylabel("a₂ = (d₊ − d₋)/2 (mÅ)")
        left.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="psi_splitting",
        title="ψ-splitting: the two branches separated",
        caption=(
            "Left: the mean a₁ of the spacings at +ψ and −ψ, free of out-of-plane shear, against "
            "sin²ψ. Right: their half-difference a₂ against sin|2ψ|, with the fitted "
            "splitting term as a line through the origin."
        ),
        interpretation=(
            "a₂ = d₀·½S₂·τφ·sin|2ψ|: a₂ flat at zero means no out-of-plane shear, a straight "
            "line through the origin measures it. This is the Dölle–Hauk evaluation, and the "
            "reason tilts of both signs are worth the extra time."
        ),
        height_in=2.8,
    )


def sigma_phi_figure(result: ResidualStressResult) -> ResultFigure | None:
    """σφ from each azimuth's slope, against the fitted tensor as a curve in φ."""

    fit = result.tensor
    if fit is None:
        return None

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        grid = np.linspace(0.0, 180.0, 181)
        values, uncertainty = tensor_sigma_phi(fit, grid)
        axes.fill_between(
            grid,
            values - uncertainty,
            values + uncertainty,
            color=COLORS["band_2"],
            lw=0,
            label="Tensor ± u",
        )
        axes.plot(grid, values, lw=1.2, color=COLORS["model"], label="Stress tensor")
        for index, line in enumerate(result.regressions):
            axes.errorbar(
                [line.phi_deg % 180.0],
                [line.sigma_phi_mpa],
                yerr=[line.sigma_phi_uncertainty_mpa],
                fmt="o",
                ms=5,
                color=_azimuth_color(index),
                capsize=2.5,
            )
        axes.axhline(0.0, lw=0.6, color=COLORS["guide"])
        principal = fit.in_plane_principal()
        axes.axvline(
            principal["angle_deg"] % 180.0,
            lw=0.7,
            ls=":",
            color=COLORS["accent"],
            label="Direction of σI",
        )
        axes.set_xlim(0.0, 180.0)
        axes.set_xticks(range(0, 181, 30))
        axes.set_xlabel(f"Azimuth {_math('stress_azimuth')} (°)")
        axes.set_ylabel(f"$\\sigma_{{{symbol_latex('stress_azimuth')}}}$ (MPa)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="sigma_phi",
        title="Normal stress along each azimuth",
        caption=(
            "Points: σφ from the slope of d against sin²ψ at each measured azimuth, ± its "
            "statistical uncertainty. Curve: σ11 cos²φ + σ12 sin2φ + σ22 sin²φ of the jointly "
            "fitted tensor, with its combined uncertainty as a band. The dotted line marks the "
            "direction of the larger principal stress."
        ),
        interpretation=(
            "An in-plane stress is a quadratic form in the azimuth, so σφ must follow this "
            "curve; three azimuths fix it and further ones test it. The curve's extremes are "
            "the principal stresses."
        ),
        height_in=3.0,
    )


def mohr_figure(result: ResidualStressResult) -> ResultFigure | None:
    """Mohr's circle of the in-plane stress."""

    fit = result.tensor
    if fit is None:
        return None
    s11, _ = fit.component("sigma_11")
    s22, _ = fit.component("sigma_22")
    s12, _ = fit.component("sigma_12")
    principal = fit.in_plane_principal()

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        centre = 0.5 * (s11 + s22)
        radius = math.hypot(0.5 * (s11 - s22), s12)
        angle = np.linspace(0.0, 2.0 * np.pi, 361)
        axes.plot(
            centre + radius * np.cos(angle), radius * np.sin(angle), lw=1.2, color=COLORS["model"]
        )
        axes.plot([s11, s22], [s12, -s12], "-", lw=0.8, color=COLORS["guide"])
        axes.plot([s11], [s12], "o", ms=5, color=_azimuth_color(0), label="(σ11, σ12) — face S1")
        axes.plot([s22], [-s12], "s", ms=5, color=_azimuth_color(2), label="(σ22, −σ12) — face S2")
        axes.plot(
            [principal["sigma_I_mpa"], principal["sigma_II_mpa"]],
            [0.0, 0.0],
            "D",
            ms=5,
            color=COLORS["accent"],
            label="σI, σII",
        )
        axes.axhline(0.0, lw=0.6, color=COLORS["guide"])
        axes.axvline(0.0, lw=0.6, color=COLORS["guide"])
        axes.set_aspect("equal", adjustable="datalim")
        axes.set_xlabel("Normal stress (MPa)")
        axes.set_ylabel("Shear stress (MPa)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="mohr_circle",
        title="Mohr's circle of the in-plane stress",
        caption=(
            f"The in-plane stress state: centre (σ11 + σ22)/2 = {0.5 * (s11 + s22):.0f} MPa, "
            f"radius {math.hypot(0.5 * (s11 - s22), s12):.0f} MPa. The circle meets the normal-"
            "stress axis at the principal stresses."
        ),
        interpretation=(
            "A circle entirely left of zero is compressive in every in-plane direction, the "
            "usual result of shot peening or grinding; one straddling zero is tensile along "
            "some directions and compressive along others."
        ),
        width_in=4.6,
        height_in=3.4,
    )


def budget_figure(result: ResidualStressResult) -> ResultFigure | None:
    """The uncertainty budget per component, source by source."""

    fit = result.tensor
    if fit is None:
        return None
    sources = list(fit.budget_mpa)
    labels = {
        "statistical": "Statistical (peak positions)",
        "d0": "d₀",
        "elastic constants": "Elastic constants",
    }

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        names = [COMPONENT_LABELS[name] for name in fit.component_names]
        count = len(sources) + 1 + (1 if fit.monte_carlo_uncertainty_mpa is not None else 0)
        width = 0.8 / count
        positions = np.arange(len(names))
        palette = [COLORS["difference"], COLORS["model"], COLORS["accent"]]
        for index, source in enumerate(sources):
            axes.bar(
                positions + (index - count / 2 + 0.5) * width,
                fit.budget_mpa[source],
                width,
                color=palette[index % len(palette)],
                label=labels.get(source, source),
            )
        slot = len(sources)
        axes.bar(
            positions + (slot - count / 2 + 0.5) * width,
            fit.combined_uncertainty_mpa,
            width,
            color=COLORS["data"],
            label="Combined",
        )
        if fit.monte_carlo_uncertainty_mpa is not None:
            slot += 1
            axes.bar(
                positions + (slot - count / 2 + 0.5) * width,
                fit.monte_carlo_uncertainty_mpa,
                width,
                color="white",
                edgecolor=COLORS["data"],
                hatch="///",
                label=f"Monte Carlo ({fit.monte_carlo_draws} draws)",
            )
        axes.set_xticks(positions, labels=names)
        axes.set_ylabel("Standard uncertainty (MPa)")
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="uncertainty_budget",
        title="Uncertainty budget",
        caption=(
            "Standard uncertainty of each stress component from each source, their combination "
            "in quadrature, and — hatched — the standard deviation over Monte Carlo draws of "
            "every input from its own uncertainty."
        ),
        interpretation=(
            "The tallest coloured bar is the one worth reducing: counting time or more tilts "
            "for the statistical part, a stress-free reference for d₀, measured constants for "
            "the elastic part. The hatched bar matching the black one confirms that the linear "
            "propagation is adequate."
        ),
        height_in=3.0,
    )


def stress_correlation_figure(result: ResidualStressResult) -> ResultFigure | None:
    """The correlation matrix of the stress components."""

    fit = result.tensor
    if fit is None or len(fit.component_names) < 2:
        return None

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        draw_correlation(
            figure, axes, fit.correlation, [COMPONENT_LABELS[name] for name in fit.component_names]
        )

    size = len(fit.component_names)
    return render_figure(
        draw,
        key="stress_correlation",
        title="Correlation of the stress components",
        caption="Correlation coefficients from the combined covariance of the tensor fit.",
        interpretation=(
            "This is the combined covariance, so it includes the systematic sources: an error "
            "in d₀ or S₁ moves σ11 and σ22 together, which correlates them strongly whenever "
            "those sources dominate the budget. That is a statement about the reference, not "
            "about the measured directions; the warning above, if any, is judged on the "
            "statistical part alone."
        ),
        width_in=2.2 + 0.55 * size,
        height_in=1.6 + 0.5 * size,
    )


def strain_residual_figure(result: ResidualStressResult) -> ResultFigure | None:
    """Normalized strain residuals of the tensor fit against sin²ψ."""

    if result.tensor is None:
        return None
    normalized = normalized_strain_residuals(result)

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        kept = result.included_mask
        draw_normalized_residuals(
            axes,
            result.sin2psi[kept],
            normalized[kept],
            xlabel=_sin2psi_label(),
            ylabel="(ε_obs − ε_fit) / u(ε)",
        )
        if not np.all(kept):
            axes.plot(
                result.sin2psi[~kept],
                normalized[~kept],
                "x",
                ms=7,
                mew=1.5,
                color=COLORS["warning"],
                label="Excluded (not fitted)",
            )
            low, high = axes.get_ylim()
            span = max(abs(low), abs(high), float(np.max(np.abs(normalized[~kept]))) * 1.1)
            axes.set_ylim(-span, span)
            axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="strain_residuals",
        title="Residuals of the stress fit",
        caption=(
            "The strain each measurement leaves after the tensor fit, in units of its own "
            "standard uncertainty, against sin²ψ; bands at ±2 and ±3. Excluded measurements, "
            "which were not fitted, are red crosses."
        ),
        interpretation=(
            "About 95 % inside ±2 when the model and the peak uncertainties are right. A "
            "systematic arch against sin²ψ is curvature the straight-line model cannot follow."
        ),
        height_in=2.8,
    )


def peak_quality_figure(result: ResidualStressResult) -> ResultFigure:
    """Peak width, height and position uncertainty against tilt."""

    def draw(figure: Any) -> None:
        widths, heights, precision = figure.subplots(1, 3)
        phi = np.round(result.phi_deg % 360.0, 6)
        for index, value in enumerate(sorted(set(phi))):
            members = np.isclose(phi, value)
            color = _azimuth_color(index)
            fwhm = np.array([peak.fwhm_deg for peak in result.peaks])[members]
            height = np.array([peak.height for peak in result.peaks])[members]
            u = np.array([peak.two_theta_uncertainty_deg for peak in result.peaks])[members]
            psi = result.psi_deg[members]
            order = np.argsort(psi)
            widths.plot(
                psi[order],
                fwhm[order],
                "o-",
                ms=3,
                lw=0.8,
                color=color,
                label=f"{_phi()} = {value:g}°",
            )
            heights.plot(psi[order], height[order], "o-", ms=3, lw=0.8, color=color)
            precision.plot(psi[order], 1000.0 * u[order], "o-", ms=3, lw=0.8, color=color)
        widths.set_ylabel("FWHM (°)")
        heights.set_ylabel("Peak height (counts)")
        precision.set_ylabel("u(2θ) (m°)")
        for axes in (widths, heights, precision):
            axes.set_xlabel(f"{_math('stress_tilt')} (°)")
            axes.tick_params(labelsize=7)
        widths.legend(loc="best", frameon=False, fontsize=6)

    return render_figure(
        draw,
        key="peak_quality",
        title="Peak width, height and precision against tilt",
        caption=(
            "For every located peak: the fitted width, the height above background, and the "
            "standard uncertainty of the position, against the tilt ψ."
        ),
        interpretation=(
            "Under ω-tilting the width grows and the height falls with |ψ| (defocusing and "
            "absorption), so the high-tilt positions are the least precise, and the weighted "
            "fit trusts them least. A width that changes with the sign of ψ, or differs between "
            "azimuths, points to a gradient or to texture rather than to geometry."
        ),
        height_in=2.5,
    )


def dec_figure(
    stiffness: StiffnessTensor,
    *,
    normals: Sequence[np.ndarray],
    labels: Sequence[str],
    chosen_normal: np.ndarray,
    chosen_label: str,
    used: DiffractionElasticConstants,
    gammas: Sequence[float] | None = None,
    chosen_gamma: float | None = None,
) -> ResultFigure:
    """½S₂ of several reflections under each grain-interaction model.

    ``normals`` are plane normals in the crystal Cartesian frame. For a cubic
    crystal pass the orientation parameters ``gammas`` too: the Reuss constant
    is linear in Γ, and the plot is then a set of lines rather than of points.
    """

    models: tuple[Literal["reuss", "voigt", "hill", "kroener"], ...] = (
        "reuss",
        "hill",
        "kroener",
        "voigt",
    )
    names = {"reuss": "Reuss", "hill": "Neerfeld–Hill", "kroener": "Kröner", "voigt": "Voigt"}
    table = {
        model: np.array(
            [
                DiffractionElasticConstants.from_single_crystal(
                    stiffness, normal, model=model
                ).half_s2_per_tpa
                for normal in normals
            ]
        )
        for model in models
    }
    cubic = gammas is not None and chosen_gamma is not None
    if cubic:
        abscissa = np.asarray(gammas, dtype=float)
        chosen_x = float(chosen_gamma)  # type: ignore[arg-type]
    else:
        abscissa = np.arange(len(normals), dtype=float)
        chosen_x = float(len(normals))
    del chosen_normal

    def draw(figure: Any) -> None:
        axes = figure.subplots()
        order = np.argsort(abscissa)
        palette = {
            "reuss": COLORS["difference"],
            "hill": COLORS["background"],
            "kroener": COLORS["model"],
            "voigt": COLORS["accent"],
        }
        for model in models:
            axes.plot(
                abscissa[order],
                table[model][order],
                "o-",
                ms=3.5,
                lw=1.0,
                color=palette[model],
                label=names[model],
            )
        # Reflections sharing Γ share every cubic constant, so they share a
        # label rather than print on top of one another.
        grouped: dict[float, list[int]] = {}
        for index, x in enumerate(abscissa):
            grouped.setdefault(round(float(x), 9), []).append(index)
        for members in grouped.values():
            first = members[0]
            axes.annotate(
                ", ".join(labels[index] for index in members),
                (abscissa[first], max(table[model][first] for model in models)),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                fontsize=6.5,
                color=COLORS["background"],
            )
        axes.plot(
            [chosen_x],
            [used.half_s2_per_tpa],
            "*",
            ms=11,
            color=COLORS["warning"],
            label=f"Used: {chosen_label}",
        )
        if cubic:
            axes.set_xlabel("Cubic orientation parameter Γ = (h²k² + k²l² + l²h²)/(h² + k² + l²)²")
        else:
            axes.set_xticks([*abscissa.tolist(), chosen_x], labels=[*labels, chosen_label])
            axes.set_xlabel("Reflection")
        axes.set_ylabel(f"{_math('dec_half_s2')} (TPa⁻¹)")
        axes.margins(y=0.12)
        axes.legend(loc="best", frameon=False, fontsize=7)

    return render_figure(
        draw,
        key="elastic_constants",
        title="Diffraction elastic constants by reflection and model",
        caption=(
            "½S₂ of low-index reflections from the single-crystal stiffness under four grain-"
            "interaction models; the star marks the constant this evaluation used."
            + (" For a cubic crystal the Reuss constant is linear in Γ." if cubic else "")
        ),
        interpretation=(
            "The Reuss and Voigt models bound the aggregate; Kröner's self-consistent model "
            "sits between them and is usually closest to measurement. The spread at the "
            "chosen reflection is a fair measure of how uncertain the constant is, and the "
            "stress scales as its inverse."
        ),
        height_in=3.1,
    )
