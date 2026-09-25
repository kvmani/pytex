"""Which grain to cut: map-wide ranking and the preparability raster (section 7.2).

The primary workflow. Given a segmented EBSD map and a target zone axis, rank
every grain by how good a lamella it would make, and draw, for every point of
the map, how far the target lies from the surface plane there --- the
*preparability raster* an operator reads over the SEM image to see at a glance
where the good sites are.

The score
---------
Transparent and weighted, with every weight a user-visible input stated back in
:meth:`pytex.fib.report.LamellaPlanReport.describe`:

``eps``
    ``1 - (eps* + U) / eps_limit``, clipped to ``[0, 1]``: the upper bound of
    the residual tilt against the largest residual the holder reaches for every
    mounting rotation. The primary term.
``fit``
    ``margin / margin_target`` clipped to ``[0, 1]`` when the rectangle fits
    inside the grain, else ``0``.
``reliability``
    ``exp(-GOS / gos_scale)``: grain orientation spread as a proxy for how well
    one orientation describes the grain.
``edge``
    ``distance to the map edge / edge_target``, clipped: sites near the edge
    of the scan are where the least is known about the surroundings.
``neighbour``
    ``1`` when the grain lies wholly inside the scan, ``0`` when the scan edge
    truncates it --- its extent and subsurface are then partly unseen.
``secondary``
    Optional: ``1 - angle/90`` for the secondary edge-on plane (section 7.6).

``score = sum(w_i term_i) / sum(w_i)``. Candidates are ordered first by
feasibility class (guaranteed, probabilistic, unreachable) and then by score, so
no score can lift an unreachable grain above a reachable one.

Cost
----
The raster, the per-grain residuals and a first extent check are vectorized
over every point and every grain. The full plan --- the ``phi`` sweep and the
footprint fit --- runs only for a shortlist, whose length is an input.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord
from pytex.fib.frames import ChamberGeometry, ImageRegistration, SurfaceGeometry
from pytex.fib.geometry import TargetOrbit, lamella_geometry_batch, target_orbit
from pytex.fib.placement import fit_footprint
from pytex.fib.planning import (
    FeasibilityClass,
    LamellaPlan,
    LamellaSpec,
    MountModel,
    RiskFlag,
    SecondaryPlane,
    UncertaintyInputs,
    _curve_inside,
    default_tem_stage,
    plan_lamella,
)
from pytex.tem.stage import RectangularEnvelope, StageModel, TiltEnvelope

__all__ = [
    "LamellaRankingWeights",
    "PreparabilityRaster",
    "guaranteed_radius_deg",
    "preparability_raster",
    "rank_grains",
]

_FEASIBILITY_ORDER = {
    FeasibilityClass.GUARANTEED: 0,
    FeasibilityClass.PROBABILISTIC: 1,
    FeasibilityClass.UNREACHABLE: 2,
}


@dataclass(frozen=True, slots=True)
class LamellaRankingWeights:
    """Weights and scales of the ranking score; all user-visible.

    Attributes
    ----------
    eps, fit, reliability, edge, neighbour, secondary : float
        Non-negative weights; normalized by their sum (``secondary`` counts only
        when a secondary plane was given).
    margin_target_um : float, default 2
        Clearance at which the fit term saturates.
    gos_scale_deg : float, default 1
        GOS at which the reliability term falls to ``1/e``.
    edge_target_um : float, default 10
        Distance from the map edge at which the edge term saturates.
    """

    eps: float = 0.5
    fit: float = 0.2
    reliability: float = 0.1
    edge: float = 0.1
    neighbour: float = 0.1
    secondary: float = 0.0
    margin_target_um: float = 2.0
    gos_scale_deg: float = 1.0
    edge_target_um: float = 10.0

    def __post_init__(self) -> None:
        for name in ("eps", "fit", "reliability", "edge", "neighbour", "secondary"):
            if float(getattr(self, name)) < 0.0:
                raise ValueError(f"LamellaRankingWeights.{name} must be non-negative.")
        if self.eps + self.fit + self.reliability + self.edge + self.neighbour <= 0.0:
            raise ValueError("At least one ranking weight must be positive.")
        for name in ("margin_target_um", "gos_scale_deg", "edge_target_um"):
            if float(getattr(self, name)) <= 0.0:
                raise ValueError(f"LamellaRankingWeights.{name} must be positive.")

    def terms(self, values: dict[str, float], *, with_secondary: bool) -> float:
        """The normalized weighted sum of the named terms."""

        weights = {
            "eps": self.eps,
            "fit": self.fit,
            "reliability": self.reliability,
            "edge": self.edge,
            "neighbour": self.neighbour,
        }
        if with_secondary:
            weights["secondary"] = self.secondary
        total = sum(weights.values())
        return float(sum(weights[key] * values.get(key, 0.0) for key in weights) / total)

    def describe(self) -> str:
        """The weights and scales in one sentence."""

        return (
            f"Ranking weights: eps {self.eps:g}, fit {self.fit:g}, "
            f"reliability {self.reliability:g}, "
            f"edge {self.edge:g}, neighbour {self.neighbour:g}, secondary {self.secondary:g} "
            f"(normalized by their sum); the fit term saturates at {self.margin_target_um:g} um "
            f"clearance, reliability is exp(-GOS / {self.gos_scale_deg:g} deg), the edge term "
            f"saturates {self.edge_target_um:g} um from the map edge."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "eps": float(self.eps),
            "fit": float(self.fit),
            "reliability": float(self.reliability),
            "edge": float(self.edge),
            "neighbour": float(self.neighbour),
            "secondary": float(self.secondary),
            "margin_target_um": float(self.margin_target_um),
            "gos_scale_deg": float(self.gos_scale_deg),
            "edge_target_um": float(self.edge_target_um),
        }


@dataclass(frozen=True, slots=True)
class PreparabilityRaster:
    """``eps*`` at every point of the map, for display over the SEM image.

    Attributes
    ----------
    eps_deg : (rows, cols) array
        NaN where the point belongs to another phase.
    theta_sample_deg : (rows, cols) array
        The lamella-normal azimuth each point would need.
    step_um : (dx, dy)
    origin_um : (x0, y0)
    target_text : str
    """

    eps_deg: np.ndarray
    theta_sample_deg: np.ndarray
    step_um: tuple[float, float]
    origin_um: tuple[float, float]
    target_text: str

    def __post_init__(self) -> None:
        for name in ("eps_deg", "theta_sample_deg"):
            array = np.ascontiguousarray(getattr(self, name), dtype=np.float64)
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @property
    def extent_um(self) -> tuple[float, float, float, float]:
        """``(x_min, x_max, y_min, y_max)`` of the pixel edges."""

        rows, cols = self.eps_deg.shape
        dx, dy = self.step_um
        x0, y0 = self.origin_um
        return (x0 - dx / 2, x0 + (cols - 0.5) * dx, y0 - dy / 2, y0 + (rows - 0.5) * dy)

    def fraction_below(self, eps_deg: float) -> float:
        """Fraction of the phase's points with ``eps*`` at or below a threshold."""

        finite = np.isfinite(self.eps_deg)
        if not finite.any():
            return 0.0
        return float(np.mean(self.eps_deg[finite] <= eps_deg))

    def to_json_dict(self) -> dict[str, Any]:
        """Summary payload (the raster itself travels as an image, not as JSON)."""

        finite = self.eps_deg[np.isfinite(self.eps_deg)]
        return {
            "shape": list(self.eps_deg.shape),
            "step_um": list(self.step_um),
            "origin_um": list(self.origin_um),
            "target_text": self.target_text,
            "eps_min_deg": float(finite.min()) if finite.size else None,
            "eps_median_deg": float(np.median(finite)) if finite.size else None,
            "eps_max_deg": float(finite.max()) if finite.size else None,
            "fraction_below_10_deg": self.fraction_below(10.0),
            "fraction_below_20_deg": self.fraction_below(20.0),
        }


def _grid(crystal_map: Any) -> tuple[int, int, tuple[float, float], tuple[float, float]]:
    if crystal_map.grid_kind == "hexagonal":
        raise ValueError(
            "FIB lamella planning needs a square-grid scan; resample the hexagonal scan to a "
            "square grid first."
        )
    rows, cols = crystal_map._require_regular_2d_grid()
    steps = crystal_map.step_sizes
    if steps is None or len(steps) < 2:
        raise ValueError("The scan declares no step sizes; placement needs micrometres.")
    coordinates = np.asarray(crystal_map.coordinates, dtype=np.float64)
    origin = (float(coordinates[0, 0]), float(coordinates[0, 1]))
    return int(rows), int(cols), (float(steps[0]), float(steps[1])), origin


def _phase_points(crystal_map: Any, phase: Phase | None) -> tuple[Phase, np.ndarray]:
    """The phase to plan for and the per-point mask of its points."""

    count = len(crystal_map.orientations)
    if phase is None:
        phase = crystal_map.primary_phase
        if phase is None:
            raise ValueError("The map declares no phase; pass phase=... explicitly.")
    if crystal_map.is_multiphase:
        mask = np.asarray(crystal_map.phase_mask(phase), dtype=bool)
    else:
        mask = np.ones(count, dtype=bool)
    if not mask.any():
        raise ValueError(f"No point of the map belongs to phase {phase.name!r}.")
    return phase, mask


def preparability_raster(
    crystal_map: Any,
    target: Any,
    *,
    phase: Phase | None = None,
    surface: SurfaceGeometry | None = None,
    orbit: TargetOrbit | None = None,
) -> PreparabilityRaster:
    """Compute ``eps*`` at every point of a square-grid map.

    Parameters
    ----------
    crystal_map : CrystalMap
    target : array_like
        ``[uvw]``.
    phase : Phase, optional
        Defaults to the map's primary phase; points of other phases are NaN.
    surface : SurfaceGeometry, optional
    orbit : TargetOrbit, optional

    Returns
    -------
    PreparabilityRaster
    """

    rows, cols, step, origin = _grid(crystal_map)
    phase, mask = _phase_points(crystal_map, phase)
    orbit = orbit or target_orbit(phase, target)
    matrices = crystal_map.orientations.as_matrices()
    batch = lamella_geometry_batch(matrices, orbit, surface)
    eps = np.where(mask, batch.eps_deg, np.nan).reshape(rows, cols)
    theta = np.where(mask, batch.theta_sample_deg, np.nan).reshape(rows, cols)
    return PreparabilityRaster(
        eps_deg=eps,
        theta_sample_deg=theta,
        step_um=step,
        origin_um=origin,
        target_text=orbit.family_text,
    )


def guaranteed_radius_deg(envelope: TiltEnvelope, *, samples: int = 360) -> float:
    """Largest residual reachable for every mounting rotation on this envelope.

    For a rectangular envelope symmetric about zero this is ``min(alpha_max,
    beta_max)`` exactly; for anything else it is found by bisection on the exact
    tilt curve used by the guaranteed class.
    """

    if isinstance(envelope, RectangularEnvelope):
        return float(
            min(
                -envelope.alpha_min_deg,
                envelope.alpha_max_deg,
                -envelope.beta_min_deg,
                envelope.beta_max_deg,
            )
        )
    phi = np.arange(samples) * (360.0 / samples)

    def inside(eps: float) -> bool:
        return bool(all(np.all(_curve_inside(envelope, eps, 1, phi, flip)) for flip in (1, -1)))

    low, high = 0.0, 89.9
    if not inside(low):
        return 0.0
    for _ in range(40):
        middle = 0.5 * (low + high)
        if inside(middle):
            low = middle
        else:
            high = middle
    return float(low)


def _grain_extents(
    labels_flat: np.ndarray,
    points_xy: np.ndarray,
    theta_by_label: np.ndarray,
    surface: SurfaceGeometry,
    count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extent of every grain along its own ``t_L`` and ``n_L``, vectorized."""

    sample = surface.scan_to_sample_xy(points_xy)
    valid = labels_flat >= 0
    theta = np.radians(theta_by_label[np.clip(labels_flat, 0, None)])
    along = -sample[:, 0] * np.sin(theta) + sample[:, 1] * np.cos(theta)
    across = sample[:, 0] * np.cos(theta) + sample[:, 1] * np.sin(theta)
    minimum_t = np.full(count, np.inf)
    maximum_t = np.full(count, -np.inf)
    minimum_n = np.full(count, np.inf)
    maximum_n = np.full(count, -np.inf)
    ids = labels_flat[valid]
    np.minimum.at(minimum_t, ids, along[valid])
    np.maximum.at(maximum_t, ids, along[valid])
    np.minimum.at(minimum_n, ids, across[valid])
    np.maximum.at(maximum_n, ids, across[valid])
    return maximum_t - minimum_t, maximum_n - minimum_n


def rank_grains(
    segmentation: Any,
    target: Any,
    *,
    phase: Phase | None = None,
    surface: SurfaceGeometry | None = None,
    registration: ImageRegistration | None = None,
    chamber: ChamberGeometry | None = None,
    spec: LamellaSpec | None = None,
    stage: StageModel | None = None,
    mount: MountModel | None = None,
    uncertainty: UncertaintyInputs | None = None,
    weights: LamellaRankingWeights | None = None,
    secondary: SecondaryPlane | None = None,
    both_senses: bool = True,
    shortlist: int = 8,
    min_grain_points: int = 10,
    navigation_for_best: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> Any:
    """Rank the grains of a segmented map for a target zone axis.

    Purpose
    -------
    The primary workflow of section 7.2: which grain to cut, where in it, at
    what azimuth, and how sure that is --- for every candidate at once, ordered.

    Parameters
    ----------
    segmentation : GrainSegmentation
        Of a square-grid map with step sizes.
    target : array_like
        ``[uvw]``.
    phase : Phase, optional
        Defaults to the map's primary phase.
    surface, registration, chamber, spec, stage, mount, uncertainty : optional
        As for :func:`pytex.fib.planning.plan_lamella`.
    weights : LamellaRankingWeights, optional
    secondary : SecondaryPlane, optional
    both_senses : bool, default True
    shortlist : int, default 8
        How many grains receive the full plan and footprint fit.
    min_grain_points : int, default 10
        Grains with fewer points are not considered (too small to be sure of,
        and too small to hold a lamella at any practical step).
    navigation_for_best : bool, default True
        Re-solve the best plan's ``phi`` sweep with the TEM navigation solver,
        whatever ``mount.solver`` is.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    LamellaPlanReport
    """

    from pytex.fib.report import LamellaPlanReport

    crystal_map = segmentation.crystal_map
    rows, cols, step, origin = _grid(crystal_map)
    phase, phase_mask = _phase_points(crystal_map, phase)
    surface = surface or SurfaceGeometry()
    registration = registration or ImageRegistration.identity()
    chamber = chamber or ChamberGeometry()
    spec = spec or LamellaSpec()
    stage = stage or default_tem_stage()
    mount = mount or MountModel(solver="closed_form")
    uncertainty = uncertainty or UncertaintyInputs()
    weights = weights or LamellaRankingWeights()
    if shortlist < 1:
        raise ValueError("shortlist must be at least 1.")

    orbit = target_orbit(phase, target, both_senses=both_senses)
    raster = preparability_raster(crystal_map, target, phase=phase, surface=surface, orbit=orbit)

    labels = np.asarray(segmentation.labels, dtype=np.int64)
    grains = {int(grain.grain_id): grain for grain in segmentation.grains}
    grain_count = int(labels.max()) + 1 if labels.size else 0
    sizes = np.bincount(labels[labels >= 0], minlength=grain_count)
    phase_points = np.bincount(labels[phase_mask & (labels >= 0)], minlength=grain_count)
    means = segmentation.grain_mean_orientations()
    spreads = segmentation.grain_orientation_spread_deg()

    candidate_ids = np.array(
        [
            gid
            for gid in sorted(means)
            if sizes[gid] >= min_grain_points and phase_points[gid] * 2 >= sizes[gid]
        ],
        dtype=np.int64,
    )
    notes: list[str] = []
    if candidate_ids.size == 0:
        return LamellaPlanReport(
            plans=(),
            orbit=orbit,
            phase_name=str(phase.name),
            raster=raster,
            weights=weights,
            grain_count=len(grains),
            candidate_count=0,
            guaranteed_radius_deg=guaranteed_radius_deg(stage.envelope),
            notes=(f"No grain of {phase.name} has at least {min_grain_points} points.",),
            provenance=provenance,
        )

    matrices = np.stack([means[int(gid)].as_matrix() for gid in candidate_ids])
    batch = lamella_geometry_batch(matrices, orbit, surface)
    theta_by_label = np.zeros(grain_count)
    theta_by_label[candidate_ids] = batch.theta_sample_deg
    extent_t, extent_n = _grain_extents(
        labels,
        np.asarray(crystal_map.coordinates, dtype=np.float64)[:, :2],
        theta_by_label,
        surface,
        grain_count,
    )
    spread = np.array([spreads.get(int(gid), 0.0) for gid in candidate_ids])
    combined = np.sqrt(
        uncertainty.ebsd_accuracy_deg**2
        + spread**2
        + registration.residual_deg**2
        + uncertainty.mount_repeatability_deg**2
    )
    radius = guaranteed_radius_deg(stage.envelope)
    eps_upper = batch.eps_deg + uncertainty.coverage_factor * combined
    eps_term = np.clip(1.0 - eps_upper / max(radius, 1e-9), 0.0, 1.0)
    extent_ok = (extent_t[candidate_ids] + step[0] >= spec.length_um) & (
        extent_n[candidate_ids] + step[0] >= spec.width_um
    )
    reliability = np.exp(-spread / weights.gos_scale_deg)
    pre_score = weights.eps * eps_term + weights.fit * extent_ok + weights.reliability * reliability
    # Reachable-looking grains first (eps within the envelope), then pre-score.
    order = np.lexsort((-pre_score, batch.eps_deg > radius))
    shortlisted = candidate_ids[order[:shortlist]]

    label_grid = labels.reshape(rows, cols)
    border = np.unique(
        np.concatenate([label_grid[0], label_grid[-1], label_grid[:, 0], label_grid[:, -1]])
    )
    truncated = set(int(value) for value in border if value >= 0)

    plans: list[LamellaPlan] = []
    for position, gid in enumerate(shortlisted):
        gid = int(gid)
        grain = grains[gid]
        orientation = means[gid]
        plan = plan_lamella(
            orientation.as_matrix(),
            target,
            phase=phase,
            surface=surface,
            registration=registration,
            chamber=chamber,
            spec=spec,
            stage=stage,
            mount=mount,
            uncertainty=uncertainty,
            grain_spread_deg=float(spreads.get(gid, 0.0)),
            grain_id=gid,
            location_um=(float(grain.mean_coordinate[0]), float(grain.mean_coordinate[1])),
            secondary=secondary,
            orbit=orbit,
            crosscheck=position == 0,
            provenance=provenance,
        )
        mask = label_grid == gid
        rr, cc = np.nonzero(mask)
        r0, r1 = max(int(rr.min()) - 1, 0), min(int(rr.max()) + 2, rows)
        c0, c1 = max(int(cc.min()) - 1, 0), min(int(cc.max()) + 2, cols)
        placement = fit_footprint(
            mask[r0:r1, c0:c1],
            step_um=step,
            origin_um=(origin[0] + c0 * step[0], origin[1] + r0 * step[1]),
            theta_sample_deg=plan.geometry.theta_sample_deg,
            length_um=spec.length_um,
            width_um=spec.width_um,
            surface=surface,
            registration=registration,
        )
        # The edge distance must be to the whole map, not to the cropped window.
        centre = placement.center_scan_um
        x_max = origin[0] + (cols - 1) * step[0]
        y_max = origin[1] + (rows - 1) * step[1]
        edge = max(
            0.0,
            min(centre[0] - origin[0], x_max - centre[0], centre[1] - origin[1], y_max - centre[1]),
        )
        placement = replace(placement, edge_distance_um=float(edge))
        extra: list[RiskFlag] = []
        if not placement.fits:
            extra.append(
                RiskFlag(
                    "footprint_does_not_fit",
                    "warning",
                    f"The {spec.length_um:g} um lamella does not fit inside grain {gid} at this "
                    f"azimuth (longest that fits: {placement.largest_length_um:.1f} um).",
                )
            )
        elif placement.margin_um < 1.0:
            extra.append(
                RiskFlag(
                    "small_margin",
                    "warning",
                    f"Only {placement.margin_um:.2f} um clearance to the grain boundary.",
                )
            )
        if gid in truncated:
            extra.append(
                RiskFlag(
                    "grain_truncated",
                    "warning",
                    f"Grain {gid} touches the scan edge, so part of it is unseen.",
                )
            )
        plan = plan.with_placement(placement, extra)
        terms = {
            "eps": float(np.clip(1.0 - plan.budget.eps_upper_deg / max(radius, 1e-9), 0.0, 1.0)),
            "fit": (
                float(np.clip(placement.margin_um / weights.margin_target_um, 0.0, 1.0))
                if placement.fits
                else 0.0
            ),
            "reliability": float(math.exp(-(plan.grain_spread_deg or 0.0) / weights.gos_scale_deg)),
            "edge": float(np.clip(edge / weights.edge_target_um, 0.0, 1.0)),
            "neighbour": 0.0 if gid in truncated else 1.0,
        }
        if plan.secondary is not None:
            terms["secondary"] = float(max(0.0, 1.0 - plan.secondary.angle_deg / 90.0))
        score = weights.terms(terms, with_secondary=plan.secondary is not None)
        plans.append(plan.with_ranking(score, terms))

    plans.sort(key=lambda item: (_FEASIBILITY_ORDER[item.feasibility], -(item.score or 0.0)))
    if navigation_for_best and plans and mount.solver != "navigation" and orbit.both_senses:
        best = plans[0]
        resolved = plan_lamella(
            means[int(best.grain_id)].as_matrix(),  # type: ignore[arg-type]
            target,
            phase=phase,
            surface=surface,
            registration=registration,
            chamber=chamber,
            spec=spec,
            stage=stage,
            mount=replace(mount, solver="navigation"),
            uncertainty=uncertainty,
            grain_spread_deg=best.grain_spread_deg,
            grain_id=best.grain_id,
            location_um=best.location_um,
            secondary=secondary,
            orbit=orbit,
            crosscheck=True,
            provenance=provenance,
        )
        extra_flags = [flag for flag in best.risk_flags if flag.code in _PLACEMENT_FLAGS]
        resolved = resolved.with_placement(best.placement, extra_flags)
        assert best.score is not None and best.score_terms is not None
        plans[0] = resolved.with_ranking(best.score, best.score_terms)
    if len(candidate_ids) > shortlist:
        notes.append(
            f"{len(candidate_ids)} grains were screened with the vectorized residual and extent "
            f"check; the best {shortlist} received the full plan and footprint fit."
        )
    return LamellaPlanReport(
        plans=tuple(plans),
        orbit=orbit,
        phase_name=str(phase.name),
        raster=raster,
        weights=weights,
        grain_count=len(grains),
        candidate_count=int(candidate_ids.size),
        guaranteed_radius_deg=radius,
        notes=tuple(notes),
        provenance=provenance,
    )


_PLACEMENT_FLAGS = frozenset({"footprint_does_not_fit", "small_margin", "grain_truncated"})
