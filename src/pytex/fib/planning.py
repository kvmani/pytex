"""One lamella, planned end to end: geometry, mounting, uncertainty, verdict.

Purpose
-------
Turn a measured orientation and a target zone axis into a :class:`LamellaPlan`:
the azimuth to cut at in every frame an operator will type it into, the
rectangle to mill, the residual tilt the TEM holder will have to supply, how
sure that number is, and whether the holder can supply it.

The unknown mounting rotation
-----------------------------
There is no deterministic lift-out convention (decision D6). Once the lamella is
welded to the grid, its in-plane rotation ``phi`` about the lamella normal, and
whether it was turned front-to-back, are unknown until it is looked at in the
TEM. The residual ``eps*`` is therefore a known *angle* whose decomposition into
holder tilts ``(alpha, beta)`` is not known at planning time. A plan reports
three different things and never conflates them:

1. **guaranteed** --- ``eps*`` plus its expanded uncertainty is inside the holder
   envelope by the safety margin for *every* ``phi``;
2. **probabilistic** --- the fraction of ``phi`` for which some orbit member is
   inside the envelope, from an exact solve at each sampled ``phi``;
3. **unreachable** --- no ``phi`` works.

A plan never reports one ``(alpha, beta)`` pair as if ``phi`` were known. The
small-angle picture ``alpha ~ eps cos phi, beta ~ eps sin phi`` explains the
shape of the feasible arcs, but it is not how they are computed.

Solvers
-------
``MountModel.solver = "navigation"`` (the default) solves every sampled ``phi``
with :func:`pytex.tem.navigation.plan_tilt_to_zone_axis`, the validated TEM
solver, through the full forward-validated stage model. ``"closed_form"``
evaluates the same exact trigonometry --- the four branches of
:func:`pytex.tem.navigation.solve_tilts_for_direction` --- vectorized over
``phi`` and orbit members, which is what map-wide ranking needs; the two agree
to rounding (``tests/unit/test_fib_planning.py``).

See ``docs/architecture/fib_lamella_planning_foundation.md`` section 7.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.lattice import Phase
from pytex.core.notation import format_direction_indices, format_plane_family_indices
from pytex.core.orientation import Orientation
from pytex.core.provenance import ProvenanceRecord
from pytex.fib.frames import (
    ChamberGeometry,
    ImageRegistration,
    SurfaceGeometry,
    lamella_frame_graph,
    wrap_azimuth_deg,
)
from pytex.fib.geometry import (
    LamellaGeometry,
    LamellaOption,
    TargetOrbit,
    angle_to_plane_deg,
    lamella_geometry,
    lamella_options,
    target_orbit,
)
from pytex.tem.stage import DoubleTiltStage, RectangularEnvelope, StageModel, TiltEnvelope

__all__ = [
    "LAMELLA_PLAN_SCHEMA",
    "UNVALIDATED_CHAIN_NOTE",
    "FeasibilityClass",
    "LamellaPlan",
    "LamellaSpec",
    "MountModel",
    "PhiSweep",
    "RiskFlag",
    "SecondaryPlane",
    "SecondaryScore",
    "UncertaintyBudget",
    "UncertaintyComponent",
    "UncertaintyInputs",
    "default_tem_stage",
    "plan_lamella",
    "sweep_mounting_rotation",
]

#: Schema identifier of one plan's payload.
LAMELLA_PLAN_SCHEMA = "pytex.fib_lamella_plan/1"

#: Carried by every plan until a real lamella has closed the loop (section 12).
UNVALIDATED_CHAIN_NOTE = (
    "The EBSD-to-FIB-to-TEM chain has been validated analytically and against PyTex's TEM tilt "
    "solver, but not yet against a lamella cut on a real instrument with its achieved TEM tilt "
    "recorded; treat the azimuths as a plan to check, not as a guarantee."
)

_SUBSURFACE_NOTE = (
    "The grain is assumed columnar: its orientation at the surface is taken to hold over the "
    "milling depth. EBSD sees only the top tens of nanometres, so a grain boundary below the "
    "surface would put a different crystal in the lamella."
)


def default_tem_stage() -> DoubleTiltStage:
    """Decision D7: an ideal double-tilt holder with a +/-30 x +/-30 degree envelope."""

    return DoubleTiltStage(
        envelope=RectangularEnvelope(-30.0, 30.0, -30.0, 30.0),
        name="double tilt, +/-30 deg x +/-30 deg (default)",
    )


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class LamellaSpec:
    """The slab to mill, in micrometres (thickness in nanometres).

    Attributes
    ----------
    length_um : float, default 15
        ``L``, along the long axis ``t_L``.
    width_um : float, default 2
        ``W``, along the normal ``n_L``: trench to trench, including the
        protective strap.
    depth_um : float, default 8
        ``D``, along ``-Z_s``.
    final_thickness_nm : float, default 80
        ``t``, the electron-transparent thickness after thinning.
    """

    length_um: float = 15.0
    width_um: float = 2.0
    depth_um: float = 8.0
    final_thickness_nm: float = 80.0

    def __post_init__(self) -> None:
        for name in ("length_um", "width_um", "depth_um", "final_thickness_nm"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"LamellaSpec.{name} must be positive and finite.")
        if self.width_um > self.length_um:
            raise ValueError("LamellaSpec.width_um must not exceed length_um.")
        if self.final_thickness_nm / 1000.0 >= self.width_um:
            raise ValueError("LamellaSpec.final_thickness_nm must be thinner than the width.")

    def describe(self) -> str:
        """The dimensions in one sentence."""

        return (
            f"Lamella L = {self.length_um:g} um long (along t_L), W = {self.width_um:g} um wide "
            f"(along n_L, trench to trench with the strap), D = {self.depth_um:g} um deep, thinned "
            f"to t = {self.final_thickness_nm:g} nm."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "length_um": float(self.length_um),
            "width_um": float(self.width_um),
            "depth_um": float(self.depth_um),
            "final_thickness_nm": float(self.final_thickness_nm),
        }


@dataclass(frozen=True, slots=True)
class MountModel:
    """The lamella's unknown placement on the TEM grid (decision D6).

    Attributes
    ----------
    phi_samples : int, default 72
        How many mounting rotations ``phi`` in ``[0, 360)`` are solved. Five
        degrees is finer than the arcs a holder envelope cuts.
    margin_deg : float, default 5
        Safety margin inside the envelope for the *guaranteed* class.
    solver : {"navigation", "closed_form"}
        See the module docstring.

    Notes
    -----
    The lamella-to-holder rotation for a mounting hypothesis ``(phi, s)``
    sends ``n_L`` to ``s z_H`` (the beam axis at zero tilt), ``t_L`` to
    ``(cos phi, sin phi, 0)``, and ``Z_s`` to whatever completes a proper
    rotation. ``s = +1`` is the *front* branch, ``s = -1`` the lamella turned
    over about its long axis (*back*), which reverses the beam sense.
    """

    phi_samples: int = 72
    margin_deg: float = 5.0
    solver: str = "navigation"

    def __post_init__(self) -> None:
        if int(self.phi_samples) < 4:
            raise ValueError("MountModel.phi_samples must be at least 4.")
        if not 0.0 <= float(self.margin_deg) < 45.0:
            raise ValueError("MountModel.margin_deg must lie in [0, 45).")
        if self.solver not in ("navigation", "closed_form"):
            raise ValueError("MountModel.solver must be 'navigation' or 'closed_form'.")

    @property
    def phi_deg(self) -> np.ndarray:
        """The sampled mounting rotations, in degrees."""

        return np.arange(int(self.phi_samples), dtype=np.float64) * (360.0 / self.phi_samples)

    @staticmethod
    def lamella_to_holder(phi_deg: ArrayLike, flip: int = 1) -> np.ndarray:
        """Lamella-to-holder rotation(s) for mounting rotation ``phi`` and branch ``flip``.

        Returns ``(3, 3)`` for scalar ``phi`` and ``(n, 3, 3)`` otherwise.
        """

        if int(flip) not in (-1, 1):
            raise ValueError("flip must be +1 (front) or -1 (back).")
        phi = np.radians(np.asarray(phi_deg, dtype=np.float64))
        c, s = np.cos(phi), np.sin(phi)
        zero, one = np.zeros_like(c), np.ones_like(c)
        if flip > 0:
            rows = [[zero, c, -s], [zero, s, c], [one, zero, zero]]
        else:
            rows = [[zero, c, s], [zero, s, -c], [-one, zero, zero]]
        matrix = np.moveaxis(np.array(rows, dtype=np.float64), (0, 1), (-2, -1))
        return np.asarray(matrix, dtype=np.float64)

    def describe(self) -> str:
        """The mounting model in one sentence."""

        solver_text = (
            "TEM navigation solver" if self.solver == "navigation" else "closed-form tilt solution"
        )
        return (
            f"Mounting: the in-plane rotation phi about the lamella normal and the front/back flip "
            f"are unknown until the lamella is on the holder, so {self.phi_samples} values of phi "
            f"in [0, 360) deg are solved for both branches with the "
            f"{solver_text}; "
            f"the guaranteed class keeps a {self.margin_deg:g} deg margin inside the envelope."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "phi_samples": int(self.phi_samples),
            "margin_deg": float(self.margin_deg),
            "solver": self.solver,
        }


@dataclass(frozen=True, slots=True)
class UncertaintyInputs:
    """User-supplied terms of the uncertainty budget (section 7.5).

    Attributes
    ----------
    ebsd_accuracy_deg : float, default 0.5
        Standard uncertainty of an indexed EBSD orientation.
    mount_repeatability_deg : float, default 2
        Standard uncertainty added by lift-out, welding and holder seating.
    coverage_factor : float, default 2
        ``k`` of the expanded uncertainty ``U = k u``; 2 is roughly 95 percent.
    """

    ebsd_accuracy_deg: float = 0.5
    mount_repeatability_deg: float = 2.0
    coverage_factor: float = 2.0

    def __post_init__(self) -> None:
        if self.ebsd_accuracy_deg < 0.0 or self.mount_repeatability_deg < 0.0:
            raise ValueError("Uncertainty inputs must be non-negative.")
        if not 0.0 < self.coverage_factor <= 5.0:
            raise ValueError("UncertaintyInputs.coverage_factor must lie in (0, 5].")

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "ebsd_accuracy_deg": float(self.ebsd_accuracy_deg),
            "mount_repeatability_deg": float(self.mount_repeatability_deg),
            "coverage_factor": float(self.coverage_factor),
        }


@dataclass(frozen=True, slots=True)
class SecondaryPlane:
    """An optional plane to bring towards edge-on (section 7.6), for ranking only.

    Attributes
    ----------
    label : str
        What the plane is, e.g. ``"habit plane"`` or ``"grain boundary"``.
    crystal_indices : (h, k, l) or None
        A crystallographic plane family ``{hkl}`` of the lamella's own crystal.
    sample_normal : (3,) or None
        A plane normal in sample components, e.g. a boundary normal inferred
        from its surface trace under a vertical-boundary assumption.
    """

    label: str
    crystal_indices: tuple[int, int, int] | None = None
    sample_normal: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        if (self.crystal_indices is None) == (self.sample_normal is None):
            raise ValueError("SecondaryPlane needs exactly one of crystal_indices, sample_normal.")
        if self.crystal_indices is not None and not any(self.crystal_indices):
            raise ValueError("SecondaryPlane.crystal_indices must not be (000).")
        if self.sample_normal is not None and not np.any(self.sample_normal):
            raise ValueError("SecondaryPlane.sample_normal must be non-zero.")


@dataclass(frozen=True, slots=True)
class SecondaryScore:
    """How close the secondary plane comes to edge-on at the target zone axis.

    ``angle_deg`` is the angle between the beam (along the chosen member) and
    the plane; zero is exactly edge-on. It informs ranking and never changes
    the primary solve.
    """

    label: str
    angle_deg: float
    description: str

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {"label": self.label, "angle_deg": float(self.angle_deg), "plane": self.description}


# --------------------------------------------------------------------------- #
# Results
# --------------------------------------------------------------------------- #


class FeasibilityClass(StrEnum):
    """The three verdicts of section 7.4, kept distinct."""

    #: Reachable for every mounting rotation, with margin, at the upper bound.
    GUARANTEED = "guaranteed"
    #: Reachable for some mounting rotations only; the fraction says how many.
    PROBABILISTIC = "probabilistic"
    #: No mounting rotation works.
    UNREACHABLE = "unreachable"


@dataclass(frozen=True, slots=True)
class RiskFlag:
    """One reason to be careful with a plan.

    Attributes
    ----------
    code : str
        Stable identifier, e.g. ``"uncalibrated_azimuth"``.
    severity : {"info", "warning", "critical"}
    message : str
    """

    code: str
    severity: str
    message: str

    def __post_init__(self) -> None:
        if self.severity not in ("info", "warning", "critical"):
            raise ValueError("RiskFlag.severity must be info, warning or critical.")

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {"code": self.code, "severity": self.severity, "message": self.message}


@dataclass(frozen=True, slots=True)
class UncertaintyComponent:
    """One term of the budget: a standard uncertainty in degrees and its source."""

    name: str
    value_deg: float
    source: str

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {"name": self.name, "value_deg": float(self.value_deg), "source": self.source}


@dataclass(frozen=True, slots=True)
class UncertaintyBudget:
    """Combined standard uncertainty of ``eps*`` and the upper bound the verdict uses.

    The components are combined in quadrature with unit sensitivity, which is
    an upper bound: an orientation error of ``delta`` moves the target by at
    most ``delta``, and an azimuth error moves the residual by at most as much.
    """

    components: tuple[UncertaintyComponent, ...]
    combined_deg: float
    coverage_factor: float
    eps_deg: float

    @property
    def expanded_deg(self) -> float:
        """``U = k u``."""

        return float(self.coverage_factor * self.combined_deg)

    @property
    def eps_upper_deg(self) -> float:
        """``eps* + U``, the value the feasibility verdict is judged on."""

        return float(self.eps_deg + self.expanded_deg)

    @property
    def eps_lower_deg(self) -> float:
        """``max(0, eps* - U)``."""

        return float(max(0.0, self.eps_deg - self.expanded_deg))

    def describe(self) -> str:
        """The budget as one sentence per term, then the combination."""

        terms = "; ".join(
            f"{item.name} {item.value_deg:.2f} deg ({item.source})" for item in self.components
        )
        return (
            f"Uncertainty budget (standard uncertainties, combined in quadrature with unit "
            f"sensitivity): {terms}. Combined u(eps*) = {self.combined_deg:.2f} deg; expanded "
            f"U = {self.coverage_factor:g} u = {self.expanded_deg:.2f} deg, so "
            f"eps* = {self.eps_deg:.2f} +/- {self.expanded_deg:.2f} deg and the verdict is judged "
            f"on the upper bound {self.eps_upper_deg:.2f} deg."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload."""

        return {
            "components": [item.to_json_dict() for item in self.components],
            "combined_deg": float(self.combined_deg),
            "coverage_factor": float(self.coverage_factor),
            "expanded_deg": self.expanded_deg,
            "eps_deg": float(self.eps_deg),
            "eps_upper_deg": self.eps_upper_deg,
            "eps_lower_deg": self.eps_lower_deg,
        }


def _budget(
    eps_deg: float,
    inputs: UncertaintyInputs,
    *,
    grain_spread_deg: float | None,
    registration: ImageRegistration,
) -> UncertaintyBudget:
    components = [
        UncertaintyComponent("EBSD orientation accuracy", inputs.ebsd_accuracy_deg, "user input"),
        UncertaintyComponent(
            "intragranular spread",
            float(grain_spread_deg) if grain_spread_deg is not None else 0.0,
            "grain orientation spread (GOS) from the map"
            if grain_spread_deg is not None
            else "no grain given; a single point carries none",
        ),
        UncertaintyComponent(
            "registration residual",
            float(registration.residual_deg),
            "affine control-point fit" if registration.method == "affine" else "declared, exact",
        ),
        UncertaintyComponent(
            "stage and mount repeatability", inputs.mount_repeatability_deg, "user input"
        ),
    ]
    combined = float(math.sqrt(sum(item.value_deg**2 for item in components)))
    return UncertaintyBudget(
        components=tuple(components),
        combined_deg=combined,
        coverage_factor=float(inputs.coverage_factor),
        eps_deg=float(eps_deg),
    )


@dataclass(frozen=True, slots=True)
class PhiSweep:
    """Reachability of the target over the unknown mounting rotation.

    Attributes
    ----------
    phi_deg : (n,) array
    reachable_front, reachable_back : (n,) bool arrays
        Whether some orbit member is inside the envelope at each ``phi``, for
        the front and back branches.
    alpha_front_deg, beta_front_deg, alpha_back_deg, beta_back_deg : (n,) arrays
        The least-tilt reachable solution at each ``phi`` (NaN where none).
    solver : str
    conservative_fraction : float
        Fraction of ``phi`` for which the *chosen* member, displaced to the
        upper bound ``eps* + U``, is still inside the envelope.
    guaranteed_at_point : bool
        Whether ``eps* + margin`` is inside for every ``phi``.
    guaranteed_at_upper : bool
        Whether ``eps* + U + margin`` is inside for every ``phi``.
    """

    phi_deg: np.ndarray
    reachable_front: np.ndarray
    reachable_back: np.ndarray
    alpha_front_deg: np.ndarray
    beta_front_deg: np.ndarray
    alpha_back_deg: np.ndarray
    beta_back_deg: np.ndarray
    solver: str
    conservative_fraction: float
    guaranteed_at_point: bool
    guaranteed_at_upper: bool

    def __post_init__(self) -> None:
        for name in (
            "phi_deg",
            "reachable_front",
            "reachable_back",
            "alpha_front_deg",
            "beta_front_deg",
            "alpha_back_deg",
            "beta_back_deg",
        ):
            array = np.ascontiguousarray(getattr(self, name))
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @property
    def fraction_front(self) -> float:
        """Fraction of mounting rotations reachable on the front branch."""

        return float(np.mean(self.reachable_front))

    @property
    def fraction_back(self) -> float:
        """Fraction of mounting rotations reachable on the back branch."""

        return float(np.mean(self.reachable_back))

    @property
    def fraction(self) -> float:
        """Fraction over both branches, each taken as equally likely."""

        return 0.5 * (self.fraction_front + self.fraction_back)

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload (the per-``phi`` arrays included, NaN as null)."""

        def clean(values: np.ndarray) -> list[float | None]:
            return [None if not np.isfinite(value) else float(value) for value in values]

        return {
            "solver": self.solver,
            "phi_deg": [float(value) for value in self.phi_deg],
            "reachable_front": [bool(value) for value in self.reachable_front],
            "reachable_back": [bool(value) for value in self.reachable_back],
            "alpha_front_deg": clean(self.alpha_front_deg),
            "beta_front_deg": clean(self.beta_front_deg),
            "alpha_back_deg": clean(self.alpha_back_deg),
            "beta_back_deg": clean(self.beta_back_deg),
            "fraction_front": self.fraction_front,
            "fraction_back": self.fraction_back,
            "fraction": self.fraction,
            "conservative_fraction": float(self.conservative_fraction),
            "guaranteed_at_point": bool(self.guaranteed_at_point),
            "guaranteed_at_upper": bool(self.guaranteed_at_upper),
        }


# --------------------------------------------------------------------------- #
# The phi sweep
# --------------------------------------------------------------------------- #


def _inside(envelope: TiltEnvelope, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Vectorized envelope membership, equal to ``margin_deg >= 0`` element-wise."""

    finite = np.isfinite(alpha) & np.isfinite(beta)
    if isinstance(envelope, RectangularEnvelope):
        inside = (
            (alpha >= envelope.alpha_min_deg)
            & (alpha <= envelope.alpha_max_deg)
            & (beta >= envelope.beta_min_deg)
            & (beta <= envelope.beta_max_deg)
        )
        return np.asarray(inside & finite)
    margin = np.vectorize(
        lambda a, b: envelope.margin_deg(float(a), float(b)) if np.isfinite(a) else -1.0,
        otypes=[np.float64],
    )(alpha, beta)
    return np.asarray((margin >= 0.0) & finite)


def _closed_form_branches(w: np.ndarray, *, allow_reverse: bool) -> tuple[np.ndarray, np.ndarray]:
    """The branches of ``solve_tilts_for_direction``, vectorized over ``w[..., 3]``.

    Returns alpha and beta arrays with a trailing branch axis of length two or
    four. Degenerate directions (along the beta axis) give NaN, which no
    envelope contains, matching the ``alpha = +/-90`` the scalar solver returns.
    """

    w = w / np.linalg.norm(w, axis=-1, keepdims=True)
    w1, w2, w3 = w[..., 0], w[..., 1], w[..., 2]
    rho = np.hypot(w1, w3)
    beta = np.degrees(np.arctan2(-w1, w3))
    beta_opposite = np.where(beta <= 0.0, beta + 180.0, beta - 180.0)
    alpha = np.degrees(np.arctan2(w2, rho))
    alpha_opposite = np.degrees(np.arctan2(w2, -rho))
    alphas = [alpha, alpha_opposite]
    betas = [beta, beta_opposite]
    if allow_reverse:
        alphas += [np.degrees(np.arctan2(-w2, -rho)), np.degrees(np.arctan2(-w2, rho))]
        betas += [beta, beta_opposite]
    alpha_all = np.stack(alphas, axis=-1)
    beta_all = np.stack(betas, axis=-1)
    degenerate = (rho < 1e-9)[..., None]
    return np.where(degenerate, np.nan, alpha_all), np.where(degenerate, np.nan, beta_all)


def _tilt_magnitude_deg(alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Angle between the beam at ``(alpha, beta)`` and at zero tilt."""

    cosine = np.cos(np.radians(alpha)) * np.cos(np.radians(beta))
    return np.asarray(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _crystal_to_lamella(
    orientation_matrix: np.ndarray, surface: SurfaceGeometry, geometry: LamellaGeometry
) -> np.ndarray:
    graph = lamella_frame_graph(surface, sample_to_lamella=geometry.sample_to_lamella_matrix())
    specimen_to_lamella = graph.transform_between("specimen", "lamella").rotation_matrix
    return np.asarray(specimen_to_lamella @ orientation_matrix, dtype=np.float64)


def _curve_inside(
    envelope: TiltEnvelope, eps_deg: float, rise_sign: int, phi: np.ndarray, flip: int
) -> np.ndarray:
    """Whether the tilt that removes a residual ``eps_deg`` is inside, at each ``phi``.

    The chosen member sits at ``(cos eps, 0, s sin eps)`` in the lamella frame;
    mounted at ``phi`` it needs exactly ``alpha = asin(sigma cos phi sin eps)``,
    ``beta = atan2(sigma sin phi sin eps, cos eps)`` (front branch), which is the
    closed form evaluated for that direction, not a small-angle approximation.
    """

    if eps_deg >= 90.0:
        return np.zeros(phi.shape, dtype=bool)
    eps = math.radians(eps_deg)
    member = np.array([math.cos(eps), 0.0, rise_sign * math.sin(eps)])
    w = np.einsum("nij,j->ni", MountModel.lamella_to_holder(phi, flip), member)
    alpha, beta = _closed_form_branches(w, allow_reverse=True)
    return np.asarray(np.any(_inside(envelope, alpha, beta), axis=-1))


def sweep_mounting_rotation(
    orientation_matrix: ArrayLike,
    orbit: TargetOrbit,
    geometry: LamellaGeometry,
    *,
    surface: SurfaceGeometry | None = None,
    stage: StageModel | None = None,
    mount: MountModel | None = None,
    eps_upper_deg: float | None = None,
) -> PhiSweep:
    """Solve the holder tilts for every sampled mounting rotation and both branches.

    Parameters
    ----------
    orientation_matrix : (3, 3) array_like
        ``g``, crystal to EBSD specimen.
    orbit : TargetOrbit
    geometry : LamellaGeometry
        The solved lamella for this orientation (fixes ``n_L`` and ``t_L``).
    surface : SurfaceGeometry, optional
    stage : StageModel, optional
        Defaults to :func:`default_tem_stage`.
    mount : MountModel, optional
    eps_upper_deg : float, optional
        Upper bound of ``eps*`` for the conservative quantities; defaults to
        ``eps*`` itself.

    Returns
    -------
    PhiSweep
    """

    surface = surface or SurfaceGeometry()
    stage = stage or default_tem_stage()
    mount = mount or MountModel()
    g = np.asarray(orientation_matrix, dtype=np.float64)
    phi = mount.phi_deg
    crystal_to_lamella = _crystal_to_lamella(g, surface, geometry)
    allow_reverse = bool(orbit.both_senses)

    results: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for flip in (1, -1):
        if mount.solver == "navigation" and allow_reverse:
            results[flip] = _navigation_sweep(
                crystal_to_lamella, orbit, stage, phi, flip, allow_reverse=allow_reverse
            )
        else:
            results[flip] = _closed_form_sweep(
                crystal_to_lamella, orbit, stage.envelope, phi, flip, allow_reverse=allow_reverse
            )

    rise = geometry.out_of_plane_sign
    upper = geometry.eps_deg if eps_upper_deg is None else float(eps_upper_deg)
    envelope = stage.envelope
    guaranteed_point = bool(
        np.all(_curve_inside(envelope, geometry.eps_deg + mount.margin_deg, rise, phi, 1))
        and np.all(_curve_inside(envelope, geometry.eps_deg + mount.margin_deg, rise, phi, -1))
    )
    guaranteed_upper = bool(
        np.all(_curve_inside(envelope, upper + mount.margin_deg, rise, phi, 1))
        and np.all(_curve_inside(envelope, upper + mount.margin_deg, rise, phi, -1))
    )
    conservative = 0.5 * (
        float(np.mean(_curve_inside(envelope, upper, rise, phi, 1)))
        + float(np.mean(_curve_inside(envelope, upper, rise, phi, -1)))
    )
    solver = "navigation" if (mount.solver == "navigation" and allow_reverse) else "closed_form"
    return PhiSweep(
        phi_deg=phi,
        reachable_front=results[1][0],
        reachable_back=results[-1][0],
        alpha_front_deg=results[1][1],
        beta_front_deg=results[1][2],
        alpha_back_deg=results[-1][1],
        beta_back_deg=results[-1][2],
        solver=solver,
        conservative_fraction=conservative,
        guaranteed_at_point=guaranteed_point,
        guaranteed_at_upper=guaranteed_upper,
    )


def _closed_form_sweep(
    crystal_to_lamella: np.ndarray,
    orbit: TargetOrbit,
    envelope: TiltEnvelope,
    phi: np.ndarray,
    flip: int,
    *,
    allow_reverse: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    members_lamella = orbit.cartesian @ crystal_to_lamella.T  # (m, 3)
    mounts = MountModel.lamella_to_holder(phi, flip)  # (n, 3, 3)
    w = np.einsum("nij,mj->nmi", mounts, members_lamella)  # (n, m, 3)
    alpha, beta = _closed_form_branches(w, allow_reverse=allow_reverse)  # (n, m, b)
    inside = _inside(envelope, alpha, beta)
    magnitude = np.where(inside, _tilt_magnitude_deg(alpha, beta), np.inf)
    flat = magnitude.reshape(magnitude.shape[0], -1)
    best = np.argmin(flat, axis=1)
    rows = np.arange(flat.shape[0])
    reachable = np.isfinite(flat[rows, best])
    alpha_best = np.where(reachable, alpha.reshape(flat.shape)[rows, best], np.nan)
    beta_best = np.where(reachable, beta.reshape(flat.shape)[rows, best], np.nan)
    return reachable, alpha_best, beta_best


def _navigation_sweep(
    crystal_to_lamella: np.ndarray,
    orbit: TargetOrbit,
    stage: StageModel,
    phi: np.ndarray,
    flip: int,
    *,
    allow_reverse: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from pytex.tem.navigation import plan_tilt_to_zone_axis
    from pytex.tem.reconstruction import HOLDER_FRAME, CurrentState
    from pytex.tem.stage import StagePosition

    reachable = np.zeros(phi.shape, dtype=bool)
    alpha_best = np.full(phi.shape, np.nan)
    beta_best = np.full(phi.shape, np.nan)
    mounts = MountModel.lamella_to_holder(phi, flip)
    target = [int(value) for value in orbit.target_indices]
    for index in range(phi.shape[0]):
        crystal_to_holder = mounts[index] @ crystal_to_lamella
        orientation = Orientation.from_matrix(
            crystal_to_holder, specimen_frame=HOLDER_FRAME, phase=orbit.phase
        )
        current = CurrentState.from_orientation(orientation, StagePosition(0.0, 0.0))
        report = plan_tilt_to_zone_axis(
            current,
            target,
            stage,
            allow_reverse=allow_reverse,
            include_paths=False,
            max_solutions=8,
        )
        if report.is_reachable:
            best = min(
                report.solutions,
                key=lambda item: float(
                    _tilt_magnitude_deg(
                        np.array(item.position.alpha_deg), np.array(item.position.beta_deg)
                    )
                ),
            )
            reachable[index] = True
            alpha_best[index] = best.position.alpha_deg
            beta_best[index] = best.position.beta_deg
    return reachable, alpha_best, beta_best


def _crosscheck_residual_deg(
    crystal_to_lamella: np.ndarray, orbit: TargetOrbit, *, allow_reverse: bool
) -> float:
    """The TEM solver's own answer for the residual, at ``phi = 0``, front branch.

    A holder wide enough to reach anything is used, so the least tilt the
    solver finds is the angle between the beam at zero tilt (``n_L``) and the
    nearest orbit member --- which must equal ``eps*``.
    """

    from pytex.tem.navigation import plan_tilt_to_zone_axis
    from pytex.tem.reconstruction import HOLDER_FRAME, CurrentState
    from pytex.tem.stage import StagePosition

    wide = DoubleTiltStage(envelope=RectangularEnvelope(-89.0, 89.0, -89.0, 89.0), name="wide")
    crystal_to_holder = MountModel.lamella_to_holder(0.0, 1) @ crystal_to_lamella
    orientation = Orientation.from_matrix(
        crystal_to_holder, specimen_frame=HOLDER_FRAME, phase=orbit.phase
    )
    current = CurrentState.from_orientation(orientation, StagePosition(0.0, 0.0))
    report = plan_tilt_to_zone_axis(
        current,
        [int(value) for value in orbit.target_indices],
        wide,
        allow_reverse=allow_reverse,
        include_paths=False,
        max_solutions=len(orbit) * 2,
    )
    if not report.is_reachable:
        return float("nan")
    return float(
        min(
            float(
                _tilt_magnitude_deg(
                    np.array(item.position.alpha_deg), np.array(item.position.beta_deg)
                )
            )
            for item in report.solutions
        )
    )


# --------------------------------------------------------------------------- #
# The plan
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class LamellaPlan:
    """One candidate lamella, with everything needed to cut and to judge it.

    Purpose
    -------
    The single-grain answer of section 7.1. It carries the geometry in every
    frame an operator types it into, the placement (when a grain mask was
    given), the three feasibility answers of section 7.4 kept apart, the
    uncertainty budget, and the risk flags --- and it can say all of that in
    prose (:meth:`describe`), as JSON (:meth:`to_json_dict`) and as a
    printable work order (:meth:`work_order_lines`).

    Attributes
    ----------
    phase_name : str
    orbit : TargetOrbit
    geometry : LamellaGeometry
    theta_image_deg : float
        Azimuth of ``n_L`` in the SEM image frame.
    theta_ion_deg : float
        Pattern rotation for the FIB, from the (possibly uncalibrated) chamber.
    spec, surface, registration, chamber, mount : settings used
    stage_description : str
    sweep : PhiSweep
    budget : UncertaintyBudget
    feasibility : FeasibilityClass
    crosscheck_residual_deg : float
        ``eps*`` recomputed by the TEM navigation solver; NaN if not run.
    options : tuple of LamellaOption
        The best distinct lamella azimuths, the plan's own first.
    risk_flags : tuple of RiskFlag
    grain_id : int or None
    location_um : (x, y) or None
        Where the plan applies, in scan micrometres.
    grain_spread_deg : float or None
    placement : PlacementResult or None
        From :mod:`pytex.fib.placement`.
    secondary : SecondaryScore or None
    score : float or None
        Ranking score in ``[0, 1]`` when produced by map-wide ranking.
    score_terms : mapping or None
        The weighted terms behind ``score``.
    member_directions_sample : (m, 3) array or None
        Every orbit member's unit direction in the sample frame, row for row
        with :attr:`TargetOrbit.indices`; what the stereogram draws.
    provenance : ProvenanceRecord or None
    """

    phase_name: str
    orbit: TargetOrbit
    geometry: LamellaGeometry
    theta_image_deg: float
    theta_ion_deg: float
    spec: LamellaSpec
    surface: SurfaceGeometry
    registration: ImageRegistration
    chamber: ChamberGeometry
    mount: MountModel
    stage_description: str
    sweep: PhiSweep
    budget: UncertaintyBudget
    feasibility: FeasibilityClass
    crosscheck_residual_deg: float
    options: tuple[LamellaOption, ...] = ()
    risk_flags: tuple[RiskFlag, ...] = ()
    grain_id: int | None = None
    location_um: tuple[float, float] | None = None
    grain_spread_deg: float | None = None
    placement: Any = None
    secondary: SecondaryScore | None = None
    score: float | None = None
    score_terms: dict[str, float] | None = None
    member_directions_sample: np.ndarray | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", tuple(self.options))
        object.__setattr__(self, "risk_flags", tuple(self.risk_flags))
        if self.member_directions_sample is not None:
            directions = np.ascontiguousarray(self.member_directions_sample, dtype=np.float64)
            directions.setflags(write=False)
            object.__setattr__(self, "member_directions_sample", directions)

    # -- derived quantities --------------------------------------------------

    @property
    def eps_deg(self) -> float:
        """``eps*``."""

        return float(self.geometry.eps_deg)

    @property
    def member_text(self) -> str:
        """The chosen orbit member as ``[uvw]``."""

        return format_direction_indices(
            [int(value) for value in self.geometry.member_indices], style="plain"
        )

    @property
    def long_axis_ion_deg(self) -> float:
        """Pattern rotation of the lamella *long axis*: ``theta_ion + 90`` mod 180."""

        return round(float(wrap_azimuth_deg(self.theta_ion_deg + 90.0, period=180.0)), 12)

    @property
    def is_feasible(self) -> bool:
        """Guaranteed or probabilistic."""

        return self.feasibility is not FeasibilityClass.UNREACHABLE

    def with_ranking(self, score: float, terms: dict[str, float]) -> LamellaPlan:
        """A copy carrying the ranking score and its terms."""

        from dataclasses import replace

        return replace(self, score=float(score), score_terms=dict(terms))

    def with_placement(self, placement: Any, extra_flags: Sequence[RiskFlag] = ()) -> LamellaPlan:
        """A copy carrying a footprint placement and any flags it raised."""

        from dataclasses import replace

        return replace(
            self, placement=placement, risk_flags=tuple(self.risk_flags) + tuple(extra_flags)
        )

    # -- explanation ---------------------------------------------------------

    def describe(self) -> str:
        """The plan as convention-explicit prose.

        Every number stated here is also in :meth:`to_json_dict`.
        """

        where = ""
        if self.grain_id is not None:
            where = f"grain {self.grain_id}"
            if self.grain_spread_deg is not None:
                where += f" (GOS {self.grain_spread_deg:.2f} deg)"
        if self.location_um is not None:
            where += (" at " if where else "the point at ") + (
                f"scan ({self.location_um[0]:.2f}, {self.location_um[1]:.2f}) um"
            )
        where = where or "the given orientation"
        geometry = self.geometry
        head = (
            f"Lamella for zone axis {self.orbit.family_text} in {self.phase_name}, {where}. "
            f"Of the {len(self.orbit)} members of the orbit (every symmetry equivalent"
            f"{', both senses' if self.orbit.both_senses else ', one sense only'}), "
            f"{self.member_text} lies closest to the surface plane, eps* = {self.eps_deg:.2f} deg "
            "out of it; that angle is the residual tilt the TEM holder must supply, because a "
            "vertically milled lamella can only put in-plane directions on the beam. "
        )
        if geometry.degenerate:
            head += (
                "Every member lies along the surface normal, so no vertical lamella can view this "
                "zone axis; the azimuth below is arbitrary. "
            )
        if geometry.tied_members > 1:
            head += (
                f"{geometry.tied_members} distinct azimuths reach the same eps*; the one with the "
                "smallest azimuth is given. "
            )
        azimuths = (
            f"Cut with the lamella normal n_L at theta_S = "
            f"{geometry.theta_sample_deg:.2f} deg in the "
            f"sample frame (from X_s towards Y_s; long axis t_L at "
            f"{geometry.long_axis_azimuth_deg:.2f} deg), theta_I = {self.theta_image_deg:.2f} deg "
            f"in the SEM image, and pattern rotation theta_ion = {self.theta_ion_deg:.2f} deg for "
            f"the normal ({self.long_axis_ion_deg:.2f} deg for the long axis). "
            f"{self.chamber.calibration_caveat()} "
        )
        placement = ""
        if self.placement is not None:
            placement = self.placement.describe() + " "
        sweep = self.sweep
        verdict = {
            FeasibilityClass.GUARANTEED: (
                "GUARANTEED: at the upper bound eps* + U the target stays inside the holder "
                f"envelope with a {self.mount.margin_deg:g} deg margin for every mounting rotation."
            ),
            FeasibilityClass.PROBABILISTIC: (
                f"PROBABILISTIC: reachable for {100 * sweep.fraction:.0f}% of mounting rotations "
                f"(front {100 * sweep.fraction_front:.0f}%, "
                f"back {100 * sweep.fraction_back:.0f}%), "
                f"{100 * sweep.conservative_fraction:.0f}% at the upper bound; not guaranteed."
            ),
            FeasibilityClass.UNREACHABLE: (
                "UNREACHABLE: no mounting rotation brings any orbit member inside the envelope."
            ),
        }[self.feasibility]
        tem = (
            f"TEM: {self.stage_description}. The in-plane mounting rotation phi and the "
            "front/back flip are unknown until the lamella is on the holder, so the residual is a "
            "known angle whose split into alpha and beta is not; roughly alpha ~ eps cos phi and "
            "beta ~ eps sin phi (a small-angle picture only - every number here comes from the "
            f"exact solver). {verdict} "
        )
        if math.isfinite(self.crosscheck_residual_deg):
            tem += (
                f"Cross-check: PyTex's TEM navigation solver finds a least tilt of "
                f"{self.crosscheck_residual_deg:.4f} deg for the mounted lamella, against "
                f"eps* = {self.eps_deg:.4f} deg. "
            )
        tem += (
            "The flip reverses the beam sense (the back branch looks down the reverse of the "
            "front's axis); for kinematic and centrosymmetric conditions that is immaterial, for "
            "CBED of a polar or non-centrosymmetric phase it is not. "
        )
        secondary = ""
        if self.secondary is not None:
            secondary = (
                f"Secondary (ranking only): the {self.secondary.label} "
                f"{self.secondary.description} "
                f"is {self.secondary.angle_deg:.1f} deg from edge-on. "
            )
        flags = ""
        if self.risk_flags:
            flags = (
                "Risks: "
                + " ".join(f"[{flag.severity}] {flag.message}" for flag in self.risk_flags)
                + " "
            )
        score = ""
        if self.score is not None and self.score_terms is not None:
            score = (
                f"Ranking score {self.score:.3f} from "
                + ", ".join(f"{key} {value:.3f}" for key, value in self.score_terms.items())
                + ". "
            )
        return (
            head
            + azimuths
            + placement
            + tem
            + self.budget.describe()
            + " "
            + secondary
            + score
            + flags
            + self.surface.describe()
            + " "
            + self.registration.describe()
            + " "
            + self.spec.describe()
            + " Sources: symmetry orbits after International Tables for Crystallography, "
            "Vol. A (doi:10.1107/97809553602060000100); the double-tilt closed form after De "
            "Graef, Introduction to Conventional Transmission Electron Microscopy "
            "(doi:10.1017/CBO9780511615092) and Williams & Carter, Transmission Electron "
            "Microscopy (doi:10.1007/978-0-387-76501-3)."
        ).strip()

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, in lockstep with :meth:`describe`."""

        return {
            "schema": LAMELLA_PLAN_SCHEMA,
            "phase": self.phase_name,
            "target": [int(value) for value in self.orbit.target_indices],
            "target_text": self.orbit.family_text,
            "orbit_size": len(self.orbit),
            "both_senses": bool(self.orbit.both_senses),
            "grain_id": self.grain_id,
            "grain_spread_deg": self.grain_spread_deg,
            "location_um": None if self.location_um is None else list(self.location_um),
            "member_text": self.member_text,
            "geometry": self.geometry.to_json_dict(),
            "theta_image_deg": float(self.theta_image_deg),
            "theta_ion_deg": float(self.theta_ion_deg),
            "long_axis_ion_deg": self.long_axis_ion_deg,
            "feasibility": self.feasibility.value,
            "crosscheck_residual_deg": (
                None
                if not math.isfinite(self.crosscheck_residual_deg)
                else float(self.crosscheck_residual_deg)
            ),
            "stage": self.stage_description,
            "sweep": self.sweep.to_json_dict(),
            "uncertainty": self.budget.to_json_dict(),
            "options": [option.to_json_dict() for option in self.options],
            "placement": None if self.placement is None else self.placement.to_json_dict(),
            "secondary": None if self.secondary is None else self.secondary.to_json_dict(),
            "score": self.score,
            "score_terms": self.score_terms,
            "risk_flags": [flag.to_json_dict() for flag in self.risk_flags],
            "spec": self.spec.to_json_dict(),
            "surface": self.surface.to_json_dict(),
            "registration": self.registration.to_json_dict(),
            "chamber": self.chamber.to_json_dict(),
            "mount": self.mount.to_json_dict(),
            "validation_status": UNVALIDATED_CHAIN_NOTE,
        }

    def work_order_lines(self) -> list[tuple[str, str]]:
        """The printable work order as ``(item, value)`` rows.

        Everything an operator at the FIB needs, in the order they need it, with
        the calibration caveat on ``theta_ion`` and the expected TEM residual
        with its uncertainty.
        """

        rows: list[tuple[str, str]] = [
            ("Target zone axis", f"{self.orbit.family_text} (member {self.member_text})"),
            ("Phase", self.phase_name),
        ]
        if self.grain_id is not None:
            rows.append(("Grain", str(self.grain_id)))
        placement = self.placement
        if placement is not None:
            cx, cy = placement.center_scan_um
            rows.append(("Site centre (EBSD scan)", f"x = {cx:.2f} um, y = {cy:.2f} um"))
            ix, iy = placement.center_image
            rows.append(
                (
                    "Site centre (SEM image)",
                    f"x = {ix:.2f}, y = {iy:.2f} {self.registration.image_units}",
                )
            )
        elif self.location_um is not None:
            rows.append(
                (
                    "Site (EBSD scan)",
                    f"x = {self.location_um[0]:.2f} um, y = {self.location_um[1]:.2f} um",
                )
            )
        rows += [
            (
                "Stage",
                f"tilt to T = {self.chamber.column_angle_deg:g} deg, rotation "
                f"R_stage = {self.chamber.stage_rotation_deg:g} deg",
            ),
            (
                "Pattern rotation theta_ion",
                f"{self.theta_ion_deg:.2f} deg (normal), "
                f"{self.long_axis_ion_deg:.2f} deg (long axis)",
            ),
            ("Calibration", self.chamber.calibration_caveat()),
            (
                "Lamella normal azimuth",
                f"theta_S = {self.geometry.theta_sample_deg:.2f} deg (sample), "
                f"theta_I = {self.theta_image_deg:.2f} deg (SEM image)",
            ),
            (
                "Rectangle",
                f"L = {self.spec.length_um:g} um x W = {self.spec.width_um:g} um, depth "
                f"D = {self.spec.depth_um:g} um, "
                f"final thickness {self.spec.final_thickness_nm:g} nm",
            ),
        ]
        if placement is not None:
            rows.append(
                (
                    "Clearance",
                    f"{placement.margin_um:.2f} um inside the grain"
                    if placement.fits
                    else "DOES NOT FIT; largest length that fits "
                    f"{placement.largest_length_um:.2f} um",
                )
            )
        rows += [
            (
                "Expected TEM residual tilt",
                f"eps* = {self.eps_deg:.2f} +/- {self.budget.expanded_deg:.2f} deg "
                f"(k = {self.budget.coverage_factor:g})",
            ),
            (
                "Feasibility",
                f"{self.feasibility.value}; reachable for {100 * self.sweep.fraction:.0f}% of "
                "mounting rotations",
            ),
            ("Holder", self.stage_description),
            ("Validation status", UNVALIDATED_CHAIN_NOTE),
        ]
        for flag in self.risk_flags:
            if flag.severity != "info":
                rows.append((f"Risk ({flag.severity})", flag.message))
        return rows


def _risk_flags(
    *,
    geometry: LamellaGeometry,
    chamber: ChamberGeometry,
    sweep: PhiSweep,
    budget: UncertaintyBudget,
    orbit: TargetOrbit,
    grain_spread_deg: float | None,
    crosscheck: float,
    feasibility: FeasibilityClass,
) -> list[RiskFlag]:
    flags: list[RiskFlag] = []
    if not chamber.calibrated:
        flags.append(
            RiskFlag(
                "uncalibrated_azimuth",
                "critical",
                "The FIB pattern rotation uses an uncalibrated sign and offset; calibrate with a "
                "fiducial before milling.",
            )
        )
    if not chamber.is_standard_column_angle:
        flags.append(
            RiskFlag(
                "nonstandard_column_angle",
                "warning",
                f"The column angle {chamber.column_angle_deg:g} deg is neither 52 nor 54 deg.",
            )
        )
    if sweep.guaranteed_at_point and not sweep.guaranteed_at_upper:
        flags.append(
            RiskFlag(
                "point_estimate_only",
                "warning",
                "Guaranteed only at the point estimate of eps*; its upper bound "
                f"{budget.eps_upper_deg:.1f} deg leaves the margin.",
            )
        )
    if feasibility is FeasibilityClass.PROBABILISTIC and sweep.fraction < 0.5:
        flags.append(
            RiskFlag(
                "low_mount_fraction",
                "warning",
                f"Only {100 * sweep.fraction:.0f}% of mounting rotations reach the target; "
                "expect to "
                "need luck, or a rotation holder.",
            )
        )
    if geometry.degenerate:
        flags.append(
            RiskFlag(
                "degenerate_normal",
                "critical",
                "Every orbit member is along the surface normal; "
                "a vertical lamella cannot view it.",
            )
        )
    if geometry.tied_members > 1:
        flags.append(
            RiskFlag(
                "tied_azimuths",
                "info",
                f"{geometry.tied_members} azimuths are equally good; pick by the surroundings.",
            )
        )
    if grain_spread_deg is not None and grain_spread_deg > 2.0:
        flags.append(
            RiskFlag(
                "high_spread",
                "warning",
                f"The grain's orientation spread is {grain_spread_deg:.1f} deg; the lamella may "
                "not share the mean orientation.",
            )
        )
    if math.isfinite(crosscheck) and abs(crosscheck - geometry.eps_deg) > 1e-6:
        flags.append(
            RiskFlag(
                "crosscheck_mismatch",
                "critical",
                f"The TEM solver finds {crosscheck:.4f} deg where the geometry gives "
                f"{geometry.eps_deg:.4f} deg; do not use this plan.",
            )
        )
    if not orbit.both_senses:
        flags.append(
            RiskFlag(
                "beam_sense_matters",
                "warning",
                "Only one sense of the zone axis is accepted, so the front/back flip decides "
                "success; mark the lamella's top edge before lift-out.",
            )
        )
    flags.append(RiskFlag("subsurface_assumed", "info", _SUBSURFACE_NOTE))
    flags.append(RiskFlag("unvalidated_chain", "info", UNVALIDATED_CHAIN_NOTE))
    return flags


def _secondary_score(
    secondary: SecondaryPlane,
    orbit: TargetOrbit,
    geometry: LamellaGeometry,
) -> SecondaryScore:
    member = orbit.cartesian[geometry.member_index]
    if secondary.crystal_indices is not None:
        phase = orbit.phase
        hkl = np.asarray(secondary.crystal_indices, dtype=np.float64)
        reciprocal = np.asarray(phase.lattice.reciprocal_basis().matrix, dtype=np.float64)
        normal = reciprocal @ hkl
        operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
        normals = np.einsum("nij,j->ni", operators, normal)
        cosines = np.abs(normals @ member) / np.linalg.norm(normals, axis=1)
        angle = float(np.degrees(np.arcsin(np.clip(cosines.min(), 0.0, 1.0))))
        description = "family " + format_plane_family_indices(
            [int(value) for value in secondary.crystal_indices], style="plain"
        )
    else:
        assert secondary.sample_normal is not None
        angle = angle_to_plane_deg(geometry.direction_sample, secondary.sample_normal)
        n = np.asarray(secondary.sample_normal, dtype=np.float64)
        n = n / np.linalg.norm(n)
        description = f"with sample-frame normal ({n[0]:.3f}, {n[1]:.3f}, {n[2]:.3f})"
    return SecondaryScore(label=secondary.label, angle_deg=angle, description=description)


def _coerce_orientation(
    orientation: Orientation | ArrayLike, phase: Phase | None
) -> tuple[np.ndarray, Phase]:
    if isinstance(orientation, Orientation):
        matrix = np.asarray(orientation.as_matrix(), dtype=np.float64)
        resolved = phase or orientation.phase
    else:
        matrix = np.asarray(orientation, dtype=np.float64)
        resolved = phase
    if resolved is None:
        raise ValueError(
            "plan_lamella needs a phase: pass an Orientation carrying one, or phase=..."
        )
    if matrix.shape != (3, 3):
        raise ValueError("An orientation matrix must be 3x3.")
    return matrix, resolved


def plan_lamella(
    orientation: Orientation | ArrayLike,
    target: ArrayLike,
    *,
    phase: Phase | None = None,
    surface: SurfaceGeometry | None = None,
    registration: ImageRegistration | None = None,
    chamber: ChamberGeometry | None = None,
    spec: LamellaSpec | None = None,
    stage: StageModel | None = None,
    mount: MountModel | None = None,
    uncertainty: UncertaintyInputs | None = None,
    grain_spread_deg: float | None = None,
    grain_id: int | None = None,
    location_um: tuple[float, float] | None = None,
    secondary: SecondaryPlane | None = None,
    both_senses: bool = True,
    orbit: TargetOrbit | None = None,
    crosscheck: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> LamellaPlan:
    """Plan one lamella for one orientation and one target zone axis.

    Purpose
    -------
    The single-grain solve of section 7.1: where to cut, at what azimuth in
    every frame, and whether the TEM holder will then reach the target --- for
    a grain someone has already chosen by eye, or for a point.

    When to use
    -----------
    For map-wide selection use :func:`pytex.fib.selection.rank_grains`, which
    calls this for its shortlist and adds the footprint placement. Call this
    directly when the grain is already chosen, or from a script.

    Parameters
    ----------
    orientation : Orientation or (3, 3) array_like
        Crystal to EBSD specimen (``v_specimen = g v_crystal``).
    target : array_like
        ``[uvw]`` in direct-lattice indices.
    phase : Phase, optional
        Required when ``orientation`` is a bare matrix.
    surface, registration, chamber, spec, stage, mount, uncertainty : optional
        The settings of :mod:`pytex.fib.frames` and this module, each defaulting
        to the documented default.
    grain_spread_deg : float, optional
        GOS of the grain, entering the uncertainty budget.
    grain_id : int, optional
    location_um : (x, y), optional
        Scan coordinates of the site, for the work order.
    secondary : SecondaryPlane, optional
        Informational edge-on score.
    both_senses : bool, default True
        Decision D8. False only when the beam sense has been determined and
        matters.
    orbit : TargetOrbit, optional
        Reuse a prebuilt orbit (map-wide ranking builds it once).
    crosscheck : bool, default True
        Recompute ``eps*`` through the TEM navigation solver.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    LamellaPlan

    Examples
    --------
    A cubic crystal at the identity orientation already has ``[100]`` in the
    surface, so the residual is zero and the lamella normal is ``X_s``.
    """

    matrix, resolved_phase = _coerce_orientation(orientation, phase)
    surface = surface or SurfaceGeometry()
    registration = registration or ImageRegistration.identity()
    chamber = chamber or ChamberGeometry()
    spec = spec or LamellaSpec()
    stage = stage or default_tem_stage()
    mount = mount or MountModel()
    uncertainty = uncertainty or UncertaintyInputs()
    orbit = orbit or target_orbit(resolved_phase, target, both_senses=both_senses)

    geometry = lamella_geometry(matrix, orbit, surface)
    budget = _budget(
        geometry.eps_deg,
        uncertainty,
        grain_spread_deg=grain_spread_deg,
        registration=registration,
    )
    sweep = sweep_mounting_rotation(
        matrix,
        orbit,
        geometry,
        surface=surface,
        stage=stage,
        mount=mount,
        eps_upper_deg=budget.eps_upper_deg,
    )
    if sweep.guaranteed_at_upper:
        feasibility = FeasibilityClass.GUARANTEED
    elif sweep.fraction > 0.0:
        feasibility = FeasibilityClass.PROBABILISTIC
    else:
        feasibility = FeasibilityClass.UNREACHABLE

    crosscheck_value = float("nan")
    if crosscheck and orbit.both_senses and not geometry.degenerate:
        crosscheck_value = _crosscheck_residual_deg(
            _crystal_to_lamella(matrix, surface, geometry), orbit, allow_reverse=True
        )

    flags = _risk_flags(
        geometry=geometry,
        chamber=chamber,
        sweep=sweep,
        budget=budget,
        orbit=orbit,
        grain_spread_deg=grain_spread_deg,
        crosscheck=crosscheck_value,
        feasibility=feasibility,
    )
    envelope_description = stage.envelope.describe()
    stage_name = getattr(stage, "name", "stage")
    return LamellaPlan(
        phase_name=str(resolved_phase.name),
        orbit=orbit,
        geometry=geometry,
        theta_image_deg=registration.image_azimuth_deg(geometry.theta_sample_deg, surface),
        theta_ion_deg=chamber.ion_azimuth_deg(geometry.theta_sample_deg),
        spec=spec,
        surface=surface,
        registration=registration,
        chamber=chamber,
        mount=mount,
        stage_description=f"{stage_name}; {envelope_description}",
        sweep=sweep,
        budget=budget,
        feasibility=feasibility,
        crosscheck_residual_deg=crosscheck_value,
        options=lamella_options(matrix, orbit, surface),
        risk_flags=tuple(flags),
        grain_id=grain_id,
        location_um=location_um,
        grain_spread_deg=grain_spread_deg,
        secondary=None if secondary is None else _secondary_score(secondary, orbit, geometry),
        member_directions_sample=orbit.cartesian @ (surface.specimen_to_sample_matrix() @ matrix).T,
        provenance=provenance,
    )
