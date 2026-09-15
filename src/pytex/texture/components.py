"""Named ideal texture components and component volume fractions.

Component Euler angles are Bunge (phi1, Phi, phi2) in degrees and give one
symmetry-representative of the component; volume-fraction assignment is
symmetry-aware, so the choice of representative does not matter. Angle values
follow the standard fcc rolling-texture tables (Randle & Engler).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.models import ODF, KernelSpec

ODF_COMPONENT_FIT_SCHEMA = "pytex.texture.odf_component_fit"


@dataclass(frozen=True, slots=True)
class TextureComponent:
    """A named ideal texture orientation, as Euler angles.

    Purpose
    -------
    The catalogue entries the literature names — cube, Goss, brass, copper,
    S — so that a component can be referred to by name and turned into a
    concrete orientation on a specific phase and specimen frame, rather than
    having its angles retyped at each use.

    Attributes
    ----------
    name : str
        Conventional component name.
    bunge_euler_deg : tuple of float
        The ideal orientation as Bunge ``(phi1, Phi, phi2)`` in degrees.
    Remaining attributes record the component's Miller description and any
    notes.
    """

    name: str
    bunge_euler_deg: tuple[float, float, float]
    miller_label: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("TextureComponent.name must be a non-empty string.")
        angles = tuple(float(value) for value in self.bunge_euler_deg)
        if len(angles) != 3 or not all(np.isfinite(angles)):
            raise ValueError(
                "TextureComponent.bunge_euler_deg must be three finite angles in degrees."
            )
        object.__setattr__(self, "bunge_euler_deg", angles)

    def orientation(
        self,
        *,
        specimen_frame: ReferenceFrame,
        crystal_frame: ReferenceFrame | None = None,
        symmetry: SymmetrySpec | None = None,
        phase: Phase | None = None,
    ) -> Orientation:
        """This named component as a concrete :class:`~pytex.core.orientation.Orientation`.

        Purpose
        -------
        Turn a catalogued component — cube, Goss, brass, copper, S — into an
        orientation on a specific phase and specimen frame, so it can be used as
        a volume-fraction centre, plotted, or compared against measured data.

        Parameters
        ----------
        specimen_frame : ReferenceFrame
            The specimen-domain frame. Required, because the component's Euler
            angles are defined relative to specimen axes.
        crystal_frame : ReferenceFrame, optional
            Required unless ``phase`` supplies it.
        symmetry : SymmetrySpec, optional
            Inferred from ``phase`` when omitted.
        phase : Phase, optional
            Supplies crystal frame and symmetry.
        """

        if crystal_frame is None:
            if phase is None:
                raise ValueError("crystal_frame is required when phase is not provided.")
            crystal_frame = phase.crystal_frame
        rotation = Rotation.from_euler(
            *self.bunge_euler_deg,
            convention="bunge",
            degrees=True,
        )
        return Orientation(
            rotation=rotation,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
        )


@dataclass(frozen=True, slots=True)
class ODFComponentFit:
    """Explainable non-negative fit of named components to an ODF density.

    Fractions are constrained to be non-negative and, together with the
    optional random fraction, sum to one. Observed and predicted normalized
    densities remain attached so fit quality is inspectable rather than
    inferred from the fractions alone.
    """

    components: tuple[TextureComponent, ...]
    fractions: np.ndarray
    random_fraction: float
    kernel: KernelSpec
    observed_density: np.ndarray
    predicted_density: np.ndarray
    solver_message: str = ""

    def __post_init__(self) -> None:
        components = tuple(self.components)
        if not components:
            raise ValueError("ODFComponentFit.components must not be empty.")
        if len({component.name for component in components}) != len(components):
            raise ValueError("ODFComponentFit component names must be unique.")
        fractions = np.ascontiguousarray(np.asarray(self.fractions, dtype=np.float64))
        observed = np.ascontiguousarray(np.asarray(self.observed_density, dtype=np.float64))
        predicted = np.ascontiguousarray(np.asarray(self.predicted_density, dtype=np.float64))
        if fractions.shape != (len(components),):
            raise ValueError("ODFComponentFit.fractions must provide one value per component.")
        if observed.ndim != 1 or predicted.shape != observed.shape or observed.size == 0:
            raise ValueError(
                "ODFComponentFit observed and predicted density must be non-empty 1-D arrays "
                "with identical shape."
            )
        random_fraction = float(self.random_fraction)
        if (
            np.any(~np.isfinite(fractions))
            or np.any(fractions < 0.0)
            or not np.isfinite(random_fraction)
            or random_fraction < 0.0
        ):
            raise ValueError("ODFComponentFit fractions must be finite and non-negative.")
        if not np.isclose(float(np.sum(fractions)) + random_fraction, 1.0, atol=1e-8):
            raise ValueError("ODFComponentFit component and random fractions must sum to one.")
        if np.any(~np.isfinite(observed)) or np.any(~np.isfinite(predicted)):
            raise ValueError("ODFComponentFit densities must be finite.")
        for array in (fractions, observed, predicted):
            array.setflags(write=False)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "fractions", fractions)
        object.__setattr__(self, "random_fraction", random_fraction)
        object.__setattr__(self, "observed_density", observed)
        object.__setattr__(self, "predicted_density", predicted)

    @property
    def residual(self) -> np.ndarray:
        """Observed minus predicted normalized density on the fit support."""

        values = np.ascontiguousarray(self.observed_density - self.predicted_density)
        values.setflags(write=False)
        return values

    @property
    def rms_residual(self) -> float:
        """Unweighted root-mean-square density residual."""

        return float(np.sqrt(np.mean(self.residual**2)))

    @property
    def max_absolute_residual(self) -> float:
        """Largest absolute density residual on the evaluation support."""

        return float(np.max(np.abs(self.residual)))

    @property
    def r_squared(self) -> float:
        """Unweighted coefficient of determination on the evaluation support."""

        centered = self.observed_density - float(np.mean(self.observed_density))
        denominator = float(centered @ centered)
        if np.isclose(denominator, 0.0):
            return 1.0 if np.allclose(self.residual, 0.0) else 0.0
        return 1.0 - float(self.residual @ self.residual) / denominator

    def fraction_for(self, component_name: str) -> float:
        """Fitted volume fraction for one named component."""

        for component, fraction in zip(self.components, self.fractions, strict=True):
            if component.name == component_name:
                return float(fraction)
        raise KeyError(f"Unknown fitted component {component_name!r}.")

    def describe(self) -> str:
        """Summarize fractions, residual evidence, conventions, and limits."""

        terms = ", ".join(
            f"{component.name}={float(fraction):.4f}"
            for component, fraction in zip(self.components, self.fractions, strict=True)
        )
        return (
            f"Named-component ODF fit on {self.observed_density.size} orientation(s): {terms}; "
            f"random={self.random_fraction:.4f}. Non-negative fractions sum to one and use "
            f"the {self.kernel.name} kernel at {self.kernel.halfwidth_deg:g} deg halfwidth. "
            f"Unweighted normalized-density RMS residual={self.rms_residual:.6g}, maximum "
            f"absolute residual={self.max_absolute_residual:.6g}, R^2={self.r_squared:.6g}. "
            "Fractions describe only this declared component basis and random term; omitted "
            "components, kernel choice, and non-uniform evaluation sampling can change them."
        )


CUBE = TextureComponent("cube", (0.0, 0.0, 0.0), miller_label="{001}<100>")
ROTATED_CUBE = TextureComponent("rotated_cube", (45.0, 0.0, 0.0), miller_label="{001}<110>")
GOSS = TextureComponent("goss", (0.0, 45.0, 0.0), miller_label="{011}<100>")
BRASS = TextureComponent("brass", (35.264389682754654, 45.0, 0.0), miller_label="{011}<211>")
COPPER = TextureComponent(
    "copper",
    (90.0, 35.264389682754654, 45.0),
    miller_label="{112}<111>",
)
S_COMPONENT = TextureComponent("s", (58.98, 36.7, 63.43), miller_label="{123}<634>")
ROTATED_GOSS = TextureComponent("rotated_goss", (90.0, 45.0, 0.0), miller_label="{011}<011>")

STANDARD_FCC_ROLLING_COMPONENTS: tuple[TextureComponent, ...] = (
    CUBE,
    GOSS,
    BRASS,
    COPPER,
    S_COMPONENT,
)

STANDARD_BCC_ROLLING_COMPONENTS: tuple[TextureComponent, ...] = (
    ROTATED_CUBE,
    GOSS,
    ROTATED_GOSS,
)


def component_volume_fractions(
    orientations: OrientationSet,
    components: tuple[TextureComponent, ...] | list[TextureComponent] | None = None,
    *,
    tolerance_deg: float = 15.0,
    weights: ArrayLike | None = None,
) -> dict[str, float]:
    """Return the weight fraction of orientations within tolerance of each component.

    Distances are symmetry-aware disorientation angles, so components defined
    by any symmetry-representative Euler triple behave identically. Components
    are evaluated independently: overlapping components can both claim the
    same orientation, and the fractions need not sum to one.
    """

    resolved_components = tuple(
        components if components is not None else STANDARD_FCC_ROLLING_COMPONENTS
    )
    if not resolved_components:
        raise ValueError("component_volume_fractions requires at least one component.")
    if not 0.0 < float(tolerance_deg) <= 62.8:
        raise ValueError("tolerance_deg must lie in (0, 62.8] degrees.")
    count = len(orientations)
    if count == 0:
        raise ValueError("component_volume_fractions requires at least one orientation.")
    if weights is None:
        weight_values = np.full(count, 1.0 / count, dtype=np.float64)
    else:
        weight_values = np.asarray(weights, dtype=np.float64)
        if weight_values.shape != (count,):
            raise ValueError("weights must provide one value per orientation.")
        if np.any(weight_values < 0.0):
            raise ValueError("weights must be non-negative.")
        total = float(weight_values.sum())
        if np.isclose(total, 0.0):
            raise ValueError("weights must not sum to zero.")
        weight_values = weight_values / total

    component_orientations = [
        component.orientation(
            specimen_frame=orientations.specimen_frame,
            crystal_frame=orientations.crystal_frame,
            symmetry=orientations.symmetry,
            phase=orientations.phase,
        )
        for component in resolved_components
    ]
    component_set = OrientationSet.from_orientations(component_orientations)
    angles_deg = np.rad2deg(orientations.misorientation_angles_to(component_set))
    tolerance = float(tolerance_deg)
    return {
        component.name: float(weight_values[angles_deg[:, index] <= tolerance].sum())
        for index, component in enumerate(resolved_components)
    }


def fit_odf_components(
    odf: ODF,
    components: tuple[TextureComponent, ...] | list[TextureComponent] | None = None,
    *,
    evaluation_orientations: OrientationSet | None = None,
    kernel: KernelSpec | None = None,
    include_random: bool = True,
) -> ODFComponentFit:
    """Fit a non-negative named-component mixture to normalized ODF density.

    Purpose
    -------
    Turn a measured or reconstructed ODF into interpretable named fractions
    while retaining the density residual that determines whether the chosen
    component basis is adequate.

    Method
    ------
    Each named ideal orientation contributes one symmetry-aware normalized
    kernel-density column. An optional constant column represents random
    texture. SciPy SLSQP minimizes the unweighted density residual subject to
    non-negative coefficients summing exactly to one. A rank-deficient design
    raises because its fractions are not identifiable.

    Parameters
    ----------
    odf : ODF
        ODF whose normalized density is fitted.
    components : sequence of TextureComponent, optional
        Defaults to the standard FCC rolling catalogue.
    evaluation_orientations : OrientationSet, optional
        Density sampling support. Defaults to ``odf.orientations``; supply a
        richer, approximately uniform support when the ODF support is sparse or
        strongly non-uniform.
    kernel : KernelSpec, optional
        Component peak shape; defaults to ``odf.kernel``.
    include_random : bool
        Add a constant density-one random-texture term.

    Returns
    -------
    ODFComponentFit
        Fractions, predicted/observed density, residual metrics, and explanation.
    """

    from scipy.optimize import minimize

    resolved = tuple(components if components is not None else STANDARD_FCC_ROLLING_COMPONENTS)
    if not resolved:
        raise ValueError("fit_odf_components requires at least one component.")
    if len({component.name for component in resolved}) != len(resolved):
        raise ValueError("fit_odf_components component names must be unique.")
    query = odf.orientations if evaluation_orientations is None else evaluation_orientations
    if query.crystal_frame != odf.orientations.crystal_frame:
        raise ValueError("evaluation orientations must use the ODF crystal frame.")
    if query.specimen_frame != odf.orientations.specimen_frame:
        raise ValueError("evaluation orientations must use the ODF specimen frame.")
    component_orientations = [
        component.orientation(
            specimen_frame=odf.orientations.specimen_frame,
            crystal_frame=odf.orientations.crystal_frame,
            symmetry=odf.orientations.symmetry,
            phase=odf.orientations.phase,
        )
        for component in resolved
    ]
    component_set = OrientationSet.from_orientations(component_orientations)
    angles = query.misorientation_angles_to(component_set, symmetry_aware=True)
    component_kernel = odf.kernel if kernel is None else kernel
    design = np.asarray(component_kernel.evaluate(angles, normalized=True), dtype=np.float64)
    if include_random:
        design = np.column_stack([design, np.ones(len(query), dtype=np.float64)])
    parameter_count = design.shape[1]
    if len(query) < parameter_count or np.linalg.matrix_rank(design) < parameter_count:
        raise ValueError(
            "Component design is rank deficient on the evaluation support; provide more or "
            "better-distributed evaluation orientations or fewer components."
        )
    observed = np.asarray(odf.evaluate(query, normalized=True), dtype=np.float64)

    def objective(coefficients: np.ndarray) -> float:
        residual = design @ coefficients - observed
        return 0.5 * float(residual @ residual)

    def gradient(coefficients: np.ndarray) -> np.ndarray:
        return np.asarray(design.T @ (design @ coefficients - observed), dtype=np.float64)

    initial = np.full(parameter_count, 1.0 / parameter_count, dtype=np.float64)
    solution = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * parameter_count,
        constraints={"type": "eq", "fun": lambda values: float(np.sum(values) - 1.0)},
        options={"ftol": 1e-12, "maxiter": 1000},
    )
    if not solution.success:
        raise RuntimeError(f"ODF component fit did not converge: {solution.message}")
    coefficients = np.maximum(np.asarray(solution.x, dtype=np.float64), 0.0)
    coefficients /= float(np.sum(coefficients))
    predicted = design @ coefficients
    component_fractions = coefficients[: len(resolved)]
    random_fraction = float(coefficients[-1]) if include_random else 0.0
    return ODFComponentFit(
        components=resolved,
        fractions=component_fractions,
        random_fraction=random_fraction,
        kernel=component_kernel,
        observed_density=observed,
        predicted_density=predicted,
        solver_message=str(solution.message),
    )


#: Ideal orientations of hexagonal close-packed metals, as the zirconium and
#: titanium literature names them: by where the basal pole ``[0001]`` lies.
#: With Bunge angles the crystal ``[0001]`` appears in specimen axes at
#: ``(sin phi1 sin Phi, -cos phi1 sin Phi, cos Phi)``, so ``Phi`` is the basal
#: tilt from ND and ``phi1`` the direction it tilts towards. Pinned by
#: ``tests/unit/test_texture_sections_and_sample_symmetry.py``.
HCP_BASAL = TextureComponent(
    "basal",
    (0.0, 0.0, 0.0),
    miller_label="(0001) parallel to ND",
    notes="Basal poles along the normal direction.",
)
HCP_BASAL_TD_SPLIT = TextureComponent(
    "basal_30_td",
    (180.0, 30.0, 0.0),
    miller_label="(0001) 30 deg from ND towards TD",
    notes="The split-basal texture of cold-rolled zirconium and titanium sheet.",
)
HCP_BASAL_RD_SPLIT = TextureComponent(
    "basal_30_rd",
    (90.0, 30.0, 0.0),
    miller_label="(0001) 30 deg from ND towards RD",
    notes="Basal poles tilted towards the rolling direction.",
)
HCP_C_ALONG_TD = TextureComponent(
    "c_along_td",
    (0.0, 90.0, 0.0),
    miller_label="(0001) parallel to TD",
    notes="Basal poles in the transverse (hoop) direction, as in pilgered tube.",
)
HCP_C_ALONG_RD = TextureComponent(
    "c_along_rd",
    (90.0, 90.0, 0.0),
    miller_label="(0001) parallel to RD",
    notes="Basal poles along the rolling (axial) direction.",
)

STANDARD_HCP_COMPONENTS: tuple[TextureComponent, ...] = (
    HCP_BASAL,
    HCP_BASAL_TD_SPLIT,
    HCP_BASAL_RD_SPLIT,
    HCP_C_ALONG_TD,
    HCP_C_ALONG_RD,
)


def random_component_fraction(tolerance_deg: float, symmetry_order: int = 1) -> float:
    """Volume fraction of a random texture within a misorientation of an ideal orientation.

    Purpose
    -------
    The null value a component fraction is read against. "12 percent cube" means
    nothing until it is set beside what a texture-free specimen gives for the
    same tolerance.

    Method
    ------
    Under the invariant (Haar) measure the rotation angle ``omega`` of a uniform
    rotation has density ``(1 - cos omega) / pi`` on ``[0, pi]``, so a ball of
    radius ``w`` holds ``(w - sin w) / pi`` of SO(3). An ideal orientation has
    ``|G|`` crystal-symmetry equivalents, each with its own ball; while the balls
    do not overlap - ``w`` below half the smallest symmetry rotation, 45 degrees
    for cubic and 30 for hexagonal - the fraction is ``|G| (w - sin w) / pi``.

    Parameters
    ----------
    tolerance_deg : float
        Ball radius ``w`` in degrees.
    symmetry_order : int
        Number of proper crystal symmetry operators ``|G|``.

    Returns
    -------
    float
        The random fraction, capped at 1.
    """

    omega = np.deg2rad(float(tolerance_deg))
    return float(min(1.0, int(symmetry_order) * (omega - np.sin(omega)) / np.pi))


def _ball_perturbations(tolerance_deg: float, count: int, seed: int) -> np.ndarray:
    """Rotation matrices distributed uniformly (Haar) inside a ball about the identity."""

    generator = np.random.default_rng(seed)
    omega_max = np.deg2rad(tolerance_deg)
    table = np.linspace(0.0, omega_max, 2049)
    cdf = table - np.sin(table)
    cdf /= cdf[-1]
    angles = np.interp(generator.uniform(0.0, 1.0, size=count), cdf, table)
    axes = generator.normal(size=(count, 3))
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    skew = np.zeros((count, 3, 3), dtype=np.float64)
    skew[:, 0, 1] = -axes[:, 2]
    skew[:, 0, 2] = axes[:, 1]
    skew[:, 1, 0] = axes[:, 2]
    skew[:, 1, 2] = -axes[:, 0]
    skew[:, 2, 0] = -axes[:, 1]
    skew[:, 2, 1] = axes[:, 0]
    sine = np.sin(angles)[:, None, None]
    versine = (1.0 - np.cos(angles))[:, None, None]
    squared = np.einsum("nij,njk->nik", skew, skew)
    return np.asarray(np.eye(3)[None, :, :] + sine * skew + versine * squared, dtype=np.float64)


def odf_component_volume_fractions(
    odf: Any,
    components: tuple[TextureComponent, ...] | list[TextureComponent],
    *,
    tolerance_deg: float = 15.0,
    sample_count: int = 1500,
    seed: int = 0,
) -> list[dict[str, Any]]:
    """Volume fraction of a continuous ODF within a tolerance of each ideal orientation.

    Purpose
    -------
    The quantitative reading of an ODF that a pole figure supports only by eye:
    "the basal component is 31 percent of this texture, 4.2 times random". It
    works on a reconstructed ODF - discrete or harmonic - rather than on counted
    grains, which is what a measured pole figure provides.

    Method
    ------
    The fraction is the integral of the density over the ``|G|`` symmetry-
    equivalent balls of radius ``w`` about the ideal orientation. Because the
    density is symmetry-invariant this is ``|G|`` times the integral over one
    ball, which equals the ball's volume ``(w - sin w) / pi`` times the mean
    density inside it. The mean is estimated from ``sample_count`` orientations
    drawn uniformly (Haar) inside the ball, so the result integrates the smoothed
    density rather than cutting the discrete support at a hard edge. A random
    texture therefore returns :func:`random_component_fraction` exactly, and
    ``times_random`` is the mean density in the ball in m.r.d.

    Parameters
    ----------
    odf : ODF or HarmonicODF
    components : sequence of TextureComponent
    tolerance_deg : float
        Ball radius. Keep it below half the smallest symmetry rotation (45 degrees
        cubic, 30 hexagonal), or neighbouring balls overlap and are counted twice.
    sample_count : int
        Orientations sampled per ball.
    seed : int
        Seeds the sampling, so a result is reproducible.

    Returns
    -------
    list of dict
        One entry per component: ``component``, ``miller``, ``fraction``,
        ``random_fraction`` and ``times_random``. Components are independent, so
        overlapping balls may both claim the same volume.
    """

    if not components:
        raise ValueError("odf_component_volume_fractions requires at least one component.")
    if not 0.0 < float(tolerance_deg) <= 62.8:
        raise ValueError("tolerance_deg must lie in (0, 62.8] degrees.")
    harmonic = hasattr(odf, "quadrature_orientations")
    if harmonic:
        crystal_frame, specimen_frame = odf.crystal_frame, odf.specimen_frame
        symmetry, phase = odf.crystal_symmetry, odf.phase
    else:
        support = odf.orientations
        crystal_frame, specimen_frame = support.crystal_frame, support.specimen_frame
        symmetry, phase = support.symmetry, support.phase
    order = 1 if symmetry is None else int(np.asarray(symmetry.operators).shape[0])
    perturbations = _ball_perturbations(float(tolerance_deg), int(sample_count), int(seed))
    random_fraction = random_component_fraction(float(tolerance_deg), order)
    block = 1 if harmonic else max(1, 400_000 // max(len(odf.orientations) * order, 1))
    results: list[dict[str, Any]] = []
    for component in components:
        centre = Rotation.from_euler(
            *component.bunge_euler_deg, convention="bunge", degrees=True
        ).as_matrix()
        matrices = np.einsum("ij,njk->nik", centre, perturbations)
        densities = np.empty(matrices.shape[0], dtype=np.float64)
        step = matrices.shape[0] if harmonic else block
        for start in range(0, matrices.shape[0], step):
            stop = min(start + step, matrices.shape[0])
            query = OrientationSet.from_matrices(
                matrices[start:stop],
                crystal_frame=crystal_frame,
                specimen_frame=specimen_frame,
                symmetry=symmetry,
                phase=phase,
            )
            if harmonic:
                densities[start:stop] = np.asarray(odf.evaluate(query), dtype=np.float64)
            else:
                densities[start:stop] = (
                    np.asarray(odf.evaluate(query, normalized=True), dtype=np.float64) / order
                )
        mean_density = float(np.mean(densities))
        results.append(
            {
                "component": component.name,
                "miller": component.miller_label,
                "fraction": float(min(1.0, max(0.0, random_fraction * mean_density))),
                "random_fraction": random_fraction,
                "times_random": mean_density,
            }
        )
    return results


__all__ = [
    "BRASS",
    "COPPER",
    "CUBE",
    "GOSS",
    "HCP_BASAL",
    "HCP_BASAL_RD_SPLIT",
    "HCP_BASAL_TD_SPLIT",
    "HCP_C_ALONG_RD",
    "HCP_C_ALONG_TD",
    "ODF_COMPONENT_FIT_SCHEMA",
    "ROTATED_CUBE",
    "ROTATED_GOSS",
    "STANDARD_BCC_ROLLING_COMPONENTS",
    "STANDARD_FCC_ROLLING_COMPONENTS",
    "STANDARD_HCP_COMPONENTS",
    "S_COMPONENT",
    "ODFComponentFit",
    "TextureComponent",
    "component_volume_fractions",
    "fit_odf_components",
    "odf_component_volume_fractions",
    "random_component_fraction",
]
