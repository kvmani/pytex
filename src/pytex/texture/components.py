"""Named ideal texture components and component volume fractions.

Component Euler angles are Bunge (phi1, Phi, phi2) in degrees and give one
symmetry-representative of the component; volume-fraction assignment is
symmetry-aware, so the choice of representative does not matter. Angle values
follow the standard fcc rolling-texture tables (Randle & Engler).
"""

from __future__ import annotations

from dataclasses import dataclass

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


__all__ = [
    "BRASS",
    "COPPER",
    "CUBE",
    "GOSS",
    "ODF_COMPONENT_FIT_SCHEMA",
    "ROTATED_CUBE",
    "ROTATED_GOSS",
    "STANDARD_BCC_ROLLING_COMPONENTS",
    "STANDARD_FCC_ROLLING_COMPONENTS",
    "S_COMPONENT",
    "ODFComponentFit",
    "TextureComponent",
    "component_volume_fractions",
    "fit_odf_components",
]
