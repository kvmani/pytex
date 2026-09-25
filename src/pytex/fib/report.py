"""The ranked answer, and the page an operator takes to the FIB.

Purpose
-------
:class:`LamellaPlanReport` is what map-wide ranking returns and what the
workbench shows: the candidates in order, the preparability raster behind them,
the weights that ordered them, and --- per the explainable-results doctrine ---
a prose account (:meth:`LamellaPlanReport.describe`) and a JSON payload
(:meth:`LamellaPlanReport.to_json_dict`, schema
``schemas/fib_lamella_plan.schema.json``) that state the same numbers.

:func:`work_order_html` renders one plan as a self-contained, single-page HTML
work order for printing, carrying the calibration caveat on ``theta_ion`` and
the expected TEM residual with its uncertainty next to the numbers they qualify.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from typing import Any

from pytex import __version__
from pytex.core.provenance import ProvenanceRecord
from pytex.fib.geometry import TargetOrbit
from pytex.fib.planning import (
    UNVALIDATED_CHAIN_NOTE,
    FeasibilityClass,
    LamellaPlan,
)
from pytex.fib.selection import LamellaRankingWeights, PreparabilityRaster

__all__ = [
    "LAMELLA_PLAN_REPORT_SCHEMA",
    "LamellaPlanReport",
    "work_order_html",
]

#: Schema identifier of the report payload; also the ``$id`` of
#: ``schemas/fib_lamella_plan.schema.json``.
LAMELLA_PLAN_REPORT_SCHEMA = "pytex.fib_lamella_plan_report/1"

_WHAT_IS_UNDETERMINED = (
    "What is undetermined at planning time, and what to do about it: (1) the in-plane mounting "
    "rotation phi and the front/back flip - read them off the first TEM pattern, then solve the "
    "tilt with the TEM navigation workspace; (2) the FIB rotation sign s_R and offset R0 - measure "
    "them once per instrument with the fiducial protocol; (3) the subsurface - "
    "the grain is assumed "
    "columnar over the milling depth, so check the lamella's first image for a buried boundary."
)


@dataclass(frozen=True, slots=True)
class LamellaPlanReport:
    """Ranked lamella plans for one map and one target.

    Attributes
    ----------
    plans : tuple of LamellaPlan
        Best first: feasibility class, then ranking score.
    orbit : TargetOrbit
    phase_name : str
    raster : PreparabilityRaster or None
    weights : LamellaRankingWeights
    grain_count : int
        Grains in the segmentation.
    candidate_count : int
        Grains large enough and of the right phase to be screened.
    guaranteed_radius_deg : float
        Largest residual the holder reaches for every mounting rotation.
    notes : tuple of str
    provenance : ProvenanceRecord or None
    """

    plans: tuple[LamellaPlan, ...]
    orbit: TargetOrbit
    phase_name: str
    raster: PreparabilityRaster | None
    weights: LamellaRankingWeights
    grain_count: int
    candidate_count: int
    guaranteed_radius_deg: float
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "plans", tuple(self.plans))
        object.__setattr__(self, "notes", tuple(self.notes))

    def __len__(self) -> int:
        return len(self.plans)

    def best(self) -> LamellaPlan:
        """The recommended plan, or raise with a diagnosis."""

        if not self.plans:
            raise ValueError(
                f"No candidate grain for {self.orbit.family_text} in {self.phase_name}: "
                + (" ".join(self.notes) or "the map has no usable grain.")
            )
        return self.plans[0]

    def counts(self) -> dict[str, int]:
        """Plans per feasibility class."""

        return {
            item.value: sum(1 for plan in self.plans if plan.feasibility is item)
            for item in FeasibilityClass
        }

    def describe(self) -> str:
        """The whole answer as prose: recommendation, alternatives, weights, caveats."""

        head = (
            f"FIB lamella planning for zone axis {self.orbit.family_text} in {self.phase_name}: "
            f"{self.candidate_count} of {self.grain_count} grains screened, {len(self.plans)} "
            f"planned in full. The holder reaches a residual of up to "
            f"{self.guaranteed_radius_deg:.1f} deg for every mounting rotation. "
        )
        counts = self.counts()
        head += (
            f"{counts['guaranteed']} guaranteed, {counts['probabilistic']} probabilistic, "
            f"{counts['unreachable']} unreachable. "
        )
        if self.raster is not None:
            summary = self.raster.to_json_dict()
            if summary["eps_min_deg"] is not None:
                head += (
                    f"Across the map eps* runs from {summary['eps_min_deg']:.1f} to "
                    f"{summary['eps_max_deg']:.1f} deg (median {summary['eps_median_deg']:.1f}); "
                    f"{100 * summary['fraction_below_10_deg']:.0f}% of points lie within 10 deg. "
                )
        if not self.plans:
            body = "No plan could be made. " + " ".join(self.notes)
        else:
            best = self.plans[0]
            body = "RECOMMENDED: " + best.describe() + " "
            if len(self.plans) > 1:
                alternatives = "; ".join(
                    f"grain {plan.grain_id}: eps* {plan.eps_deg:.1f} deg, "
                    f"theta_S {plan.geometry.theta_sample_deg:.1f} deg, {plan.feasibility.value}, "
                    f"score {plan.score if plan.score is not None else float('nan'):.3f}"
                    for plan in self.plans[1:6]
                )
                body += f"Alternatives: {alternatives}. "
        notes = (" ".join(self.notes) + " ") if self.notes else ""
        return (
            head
            + body
            + self.weights.describe()
            + " "
            + notes
            + _WHAT_IS_UNDETERMINED
            + " "
            + UNVALIDATED_CHAIN_NOTE
        ).strip()

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "schema": LAMELLA_PLAN_REPORT_SCHEMA,
            "pytex_version": __version__,
            "phase": self.phase_name,
            "target": [int(value) for value in self.orbit.target_indices],
            "target_text": self.orbit.family_text,
            "orbit_size": len(self.orbit),
            "grain_count": int(self.grain_count),
            "candidate_count": int(self.candidate_count),
            "guaranteed_radius_deg": float(self.guaranteed_radius_deg),
            "counts": self.counts(),
            "raster": None if self.raster is None else self.raster.to_json_dict(),
            "weights": self.weights.to_json_dict(),
            "plans": [plan.to_json_dict() for plan in self.plans],
            "undetermined": _WHAT_IS_UNDETERMINED,
            "validation_status": UNVALIDATED_CHAIN_NOTE,
            "notes": list(self.notes),
        }


_WORK_ORDER_STYLE = """
body { font-family: "DejaVu Sans", Arial, sans-serif; color: #111827; margin: 16mm;
  font-size: 10pt; }
h1 { font-size: 15pt; margin: 0 0 2mm 0; }
p.sub { margin: 0 0 4mm 0; color: #374151; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #9ca3af; padding: 1.6mm 2.2mm; vertical-align: top; text-align: left; }
th { width: 34%; background: #f3f4f6; font-weight: 600; }
tr.critical td, tr.critical th { background: #fee2e2; }
tr.warning td, tr.warning th { background: #fef3c7; }
.figures { display: flex; gap: 4mm; margin-top: 4mm; }
.figures div { flex: 1; }
.figures svg { width: 100%; height: auto; }
footer { margin-top: 4mm; font-size: 8pt; color: #4b5563; }
@page { size: A4; margin: 12mm; }
@media print { body { margin: 0; } }
"""


def work_order_html(
    plan: LamellaPlan,
    *,
    title: str | None = None,
    figures: tuple[tuple[str, str], ...] = (),
) -> str:
    """One plan as a self-contained, printable, single-page HTML work order.

    Parameters
    ----------
    plan : LamellaPlan
    title : str, optional
        Heading; defaults to the target and grain.
    figures : tuple of (caption, svg markup)
        Up to two small figures (the plan-view overlay and the section
        schematic are the useful ones) placed under the table.

    Returns
    -------
    str
        A complete HTML document with no external resource.
    """

    heading = title or (
        f"FIB lamella work order - {plan.orbit.family_text} in {plan.phase_name}"
        + (f", grain {plan.grain_id}" if plan.grain_id is not None else "")
    )
    rows = []
    for item, value in plan.work_order_lines():
        css = ""
        if item.startswith("Risk (critical)") or (
            item == "Calibration" and not plan.chamber.calibrated
        ):
            css = ' class="critical"'
        elif item.startswith("Risk (warning)"):
            css = ' class="warning"'
        rows.append(f"<tr{css}><th>{html.escape(item)}</th><td>{html.escape(value)}</td></tr>")
    figure_blocks = "".join(
        f"<div>{svg}<p>{html.escape(caption)}</p></div>" for caption, svg in figures[:2]
    )
    return (
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
        f"<title>{html.escape(heading)}</title><style>{_WORK_ORDER_STYLE}</style></head><body>"
        f"<h1>{html.escape(heading)}</h1>"
        f"<p class='sub'>Feasibility: <strong>{html.escape(plan.feasibility.value)}</strong>. "
        f"Residual tilt eps* = {plan.eps_deg:.2f} +/- {plan.budget.expanded_deg:.2f} deg.</p>"
        f"<table>{''.join(rows)}</table>"
        + (f"<div class='figures'>{figure_blocks}</div>" if figure_blocks else "")
        + f"<footer>Generated by PyTex {html.escape(__version__)} (pytex.fib). "
        "Conventions: sample frame X_s, Y_s in the surface, Z_s outward; azimuths from X_s "
        "towards Y_s. The mounting rotation phi is unknown until the lamella is on the "
        "holder; no single (alpha, beta) is implied.</footer></body></html>"
    )
