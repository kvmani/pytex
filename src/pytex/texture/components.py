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


__all__ = [
    "BRASS",
    "COPPER",
    "CUBE",
    "GOSS",
    "ROTATED_CUBE",
    "ROTATED_GOSS",
    "STANDARD_BCC_ROLLING_COMPONENTS",
    "STANDARD_FCC_ROLLING_COMPONENTS",
    "S_COMPONENT",
    "TextureComponent",
    "component_volume_fractions",
]
