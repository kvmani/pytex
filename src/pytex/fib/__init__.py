"""FIB lamella planning: which grain to cut, where, at what azimuth, for a TEM zone axis.

The question
------------
Given an EBSD scan of a bulk surface and a target TEM zone axis ``<uvw>``: which
grain should be cut, where exactly in the SEM field of view should the lamella be
placed, at what in-plane azimuth, and what residual tilt will the TEM holder then
have to supply? This subpackage answers it reproducibly, with every convention
stated.

The governing constraint
------------------------
A lamella milled strictly vertically has its faces perpendicular to the surface,
so its normal --- the TEM beam direction --- lies in the surface plane. A grain is
usable for ``<uvw>`` exactly when some symmetry equivalent lies close to that
plane, and the angle by which it misses, ``eps*``, is the tilt the holder must
supply. The crystallography is small; the bookkeeping of five frames (crystal,
sample, SEM image, FIB ion view, lamella, TEM holder) is the product.

Where to start
--------------
- :func:`rank_grains` --- the primary workflow: a segmented map in, ranked
  :class:`LamellaPlan` objects in a :class:`LamellaPlanReport` out, with the
  preparability raster.
- :func:`plan_lamella` --- one orientation (a grain chosen by eye, or a point).
- :func:`calibrate_chamber_from_fiducials` --- the one-time FIB azimuth
  calibration without which ``theta_ion`` must not be used to mill.

Every plan reports three feasibility answers kept apart (guaranteed for every
mounting rotation, the fraction of mounting rotations that work, unreachable),
because the in-plane rotation of a lamella on the TEM grid is not controlled.

Scope of version 1 (non-goals)
------------------------------
Stated here so scope creep is visible:

- tilted, wedge or pre-tilted milling --- milling is strictly vertical;
- bicrystal or interface lamellae needing a zone axis in both crystals;
- "boundary plane edge-on" as a co-equal objective --- an optional secondary
  score ranks candidates but never changes the primary solve;
- 3D or serial-section EBSD --- the subsurface is assumed columnar, and every
  plan says so;
- instrument-native pattern files for any FIB vendor --- the output is a
  printable work order and a JSON manifest (``schemas/fib_lamella_plan.schema.json``).

The end-to-end chain is validated analytically and against the TEM tilt solver
of :mod:`pytex.tem`, but not yet against a lamella cut on a real instrument; every
plan carries that statement.

See ``docs/architecture/fib_lamella_planning_foundation.md`` (specification),
``docs/site/theory/fib_lamella_zone_axis_geometry.md`` (derivation),
``docs/site/workflows/fib_lamella_planning.md`` (workbench walk-through) and
``docs/site/workflows/fib_azimuth_calibration.md`` (instrument protocol).
"""

from __future__ import annotations

from pytex.fib.frames import (
    EBSD_SPECIMEN_FRAME,
    FIB_SAMPLE_FRAME,
    ION_VIEW_FRAME,
    LAMELLA_FRAME,
    SEM_IMAGE_FRAME,
    SUPPORTED_COLUMN_ANGLES_DEG,
    ChamberGeometry,
    FiducialObservation,
    ImageRegistration,
    SurfaceGeometry,
    calibrate_chamber_from_fiducials,
    lamella_frame_graph,
    wrap_azimuth_deg,
)
from pytex.fib.geometry import (
    TIE_TOLERANCE_DEG,
    LamellaGeometry,
    LamellaGeometryBatch,
    LamellaOption,
    TargetOrbit,
    angle_to_plane_deg,
    lamella_geometry,
    lamella_geometry_batch,
    lamella_options,
    target_orbit,
)
from pytex.fib.placement import PlacementResult, fit_footprint, rectangle_corners
from pytex.fib.planning import (
    LAMELLA_PLAN_SCHEMA,
    UNVALIDATED_CHAIN_NOTE,
    FeasibilityClass,
    LamellaPlan,
    LamellaSpec,
    MountModel,
    PhiSweep,
    RiskFlag,
    SecondaryPlane,
    SecondaryScore,
    UncertaintyBudget,
    UncertaintyComponent,
    UncertaintyInputs,
    default_tem_stage,
    plan_lamella,
    sweep_mounting_rotation,
)
from pytex.fib.report import LAMELLA_PLAN_REPORT_SCHEMA, LamellaPlanReport, work_order_html
from pytex.fib.selection import (
    LamellaRankingWeights,
    PreparabilityRaster,
    guaranteed_radius_deg,
    preparability_raster,
    rank_grains,
)

__all__ = [
    "EBSD_SPECIMEN_FRAME",
    "FIB_SAMPLE_FRAME",
    "ION_VIEW_FRAME",
    "LAMELLA_FRAME",
    "LAMELLA_PLAN_REPORT_SCHEMA",
    "LAMELLA_PLAN_SCHEMA",
    "SEM_IMAGE_FRAME",
    "SUPPORTED_COLUMN_ANGLES_DEG",
    "TIE_TOLERANCE_DEG",
    "UNVALIDATED_CHAIN_NOTE",
    "ChamberGeometry",
    "FeasibilityClass",
    "FiducialObservation",
    "ImageRegistration",
    "LamellaGeometry",
    "LamellaGeometryBatch",
    "LamellaOption",
    "LamellaPlan",
    "LamellaPlanReport",
    "LamellaRankingWeights",
    "LamellaSpec",
    "MountModel",
    "PhiSweep",
    "PlacementResult",
    "PreparabilityRaster",
    "RiskFlag",
    "SecondaryPlane",
    "SecondaryScore",
    "SurfaceGeometry",
    "TargetOrbit",
    "UncertaintyBudget",
    "UncertaintyComponent",
    "UncertaintyInputs",
    "angle_to_plane_deg",
    "calibrate_chamber_from_fiducials",
    "default_tem_stage",
    "fit_footprint",
    "guaranteed_radius_deg",
    "lamella_frame_graph",
    "lamella_geometry",
    "lamella_geometry_batch",
    "lamella_options",
    "plan_lamella",
    "preparability_raster",
    "rank_grains",
    "rectangle_corners",
    "sweep_mounting_rotation",
    "target_orbit",
    "work_order_html",
    "wrap_azimuth_deg",
]
