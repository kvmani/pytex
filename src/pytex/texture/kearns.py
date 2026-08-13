"""The Kearns orientation parameter ``f`` and the pole orientation tensor.

The Kearns parameter is the single most used scalar index of texture in the
zirconium industry: ``f_d`` is the *effective fraction of basal poles aligned
with specimen direction* ``d``, defined by Kearns (WAPD-TM-472, 1965) as
``f_d = sum_i V_i cos^2(phi_i)``, the volume-weighted mean of ``cos^2`` of the
angle between each crystal's ``[0001]`` axis and ``d``. Its usefulness comes
from Kearns' Equation (1): any property that varies as
``P = P_par cos^2(phi) + P_perp sin^2(phi)`` over a hexagonal single crystal —
thermal expansion, irradiation growth, and the second-rank part of elastic and
creep response — averages over a polycrystal to ``P_d = f_d P_par + (1 - f_d) P_perp``.
One number then replaces the whole orientation distribution for that class of
property.

This module implements the four estimation routes the literature uses, all of
them as estimators of one object.

The unifying object
-------------------
Every route estimates the second-moment tensor of the basal-pole direction
distribution,

.. math::

    \\mathbf{A} = \\left\\langle \\mathbf{c}\\,\\mathbf{c}^{\\mathsf{T}} \\right\\rangle ,
    \\qquad
    f(\\mathbf{d}) = \\mathbf{d}^{\\mathsf{T}} \\mathbf{A}\\, \\mathbf{d} ,

with ``c`` the unit basal-pole direction in the specimen frame. Two properties
that the literature states as separate empirical facts follow identically:

- ``f_RD + f_TD + f_ND = tr(A) = 1`` for any orthonormal specimen triad,
  because every ``c`` is a unit vector. This is exact, not approximate, and so
  a triad sum that misses 1 is a defect in the *data or the correction*, never
  a property of the texture.
- A random texture has ``A = I/3`` and hence ``f = 1/3`` in every direction.

Working with ``A`` rather than with three separately integrated numbers is what
lets this module report a triad sum as a *diagnostic of the measurement* and
give ``f`` along any direction, not only the three principal ones.

The four routes
---------------
============================  ===========================================================
:func:`kearns_from_orientations`  Discrete orientations (EBSD, simulation). Exact: no
                                  binning, no interpolation, no truncation.
:func:`kearns_from_odf`           A fitted :class:`~pytex.texture.ODF`, with the closed-form
                                  correction for the kernel's shrinkage of ``A`` toward
                                  isotropy (:func:`kernel_axis_shrinkage`).
:func:`kearns_from_pole_figure`   A measured or computed basal pole figure, integrating
                                  Baron *et al.* Eq. (5) with solid-angle weights.
:func:`kearns_from_diffractogram` The original Kearns route: relative peak intensities of a
                                  theta-2theta scan, assigned to their known basal tilt
                                  angles, integrated by his Eq. (5).
============================  ===========================================================

References
----------
J. J. Kearns, *Thermal expansion and preferred orientation in Zircaloy*,
WAPD-TM-472, Bettis Atomic Power Laboratory (1965) — the defining report;
Eqs. (1)-(7) and Table 3.

J. L. Baron *et al.*, *Interlaboratories tests of textures of Zircaloy-4 tubes.
Part 1: pole figure measurements and calculation of Kearns coefficients*,
Textures and Microstructures **12** (1990) 125-140, doi:10.1155/TSM.12.125 —
the pole-figure route, Eqs. (4)-(5), and the incomplete-figure pseudo-norm.

R. A. Holt and S. A. Aldridge, *J. Nucl. Mater.* **135** (1985) 246-259 —
``F_d = sum V(theta) cos^2(theta)``, the resolved-basal-pole form used
throughout the CANDU pressure-tube literature.

K. V. Mani Krishna *et al.*, *J. Nucl. Mater.* **414** (2011) 492-497,
doi:10.1016/j.jnucmat.2011.04.065 — comparison of the routes, their dependence
on the measured cross-section, and the normalization the diffractogram route
needs.

See Also
--------
:doc:`../../theory/kearns_parameter_and_basal_pole_texture` : the derivation.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, freeze_array, normalize_vectors
from pytex.core.conventions import FrameDomain
from pytex.core.frames import ReferenceFrame
from pytex.core.hexagonal import is_hexagonal_phase, plane_hkl_to_hkil
from pytex.core.lattice import CrystalPlane, Phase
from pytex.core.notation import format_miller_indices
from pytex.core.orientation import OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.core.sphere import directions_to_spherical_angles, raster_solid_angle_weights
from pytex.texture.models import ODF, KernelSpec, PoleFigure

__all__ = [
    "KEARNS_ISOTROPIC_VALUE",
    "DiffractogramReflection",
    "KearnsMethod",
    "KearnsReport",
    "basal_tilt_angle_deg",
    "basal_tilt_profile",
    "harris_texture_coefficients",
    "kearns_from_diffractogram",
    "kearns_from_odf",
    "kearns_from_orientations",
    "kearns_from_pole_figure",
    "kearns_from_tilt_profile",
    "kernel_axis_shrinkage",
    "pole_orientation_tensor",
]

#: The Kearns parameter of a random (untextured) polycrystal, in every
#: direction. It is ``1/3`` because the three parameters of an orthonormal triad
#: sum identically to 1 and are equal by isotropy. Quote it whenever a measured
#: ``f`` needs a null hypothesis.
KEARNS_ISOTROPIC_VALUE: float = 1.0 / 3.0

#: Which experimental route produced a :class:`KearnsReport`. Recorded on the
#: report because the routes have genuinely different systematic errors and a
#: number quoted without its route is not comparable with another laboratory's.
KearnsMethod = Literal["orientations", "odf", "pole_figure", "tilt_profile", "diffractogram"]

_DIFFRACTOGRAM_NORMALIZATIONS = ("random_standard", "harris")

# The mean of cos^2(psi) over a sphere, used as the isotropic reference of the
# per-direction Kearns value; identical to KEARNS_ISOTROPIC_VALUE but named
# separately where the geometric meaning is what matters.
_ISOTROPIC_TENSOR = np.eye(3, dtype=np.float64) / 3.0


def _plane_label(plane: CrystalPlane, *, style: str = "plain") -> str:
    """One plane in the notation its phase's literature uses.

    Hexagonal phases get the four-index Miller-Bravais form, because that is
    what every hexagonal-symmetry argument and the whole zirconium literature
    write; other systems get three indices. Produced by
    :mod:`pytex.core.notation` rather than formatted inline.
    """

    indices = tuple(int(value) for value in plane.miller.indices)
    if is_hexagonal_phase(plane.phase):
        indices = tuple(int(value) for value in plane_hkl_to_hkil(indices))
    return format_miller_indices(indices, family="plane", style=style, scope="specific")


def _require_specimen_frame(frame: ReferenceFrame, *, argument: str) -> ReferenceFrame:
    if frame.domain is not FrameDomain.SPECIMEN:
        raise ValueError(f"{argument} must belong to the specimen domain.")
    return frame


def _default_triad(frame: ReferenceFrame) -> tuple[np.ndarray, tuple[str, ...]]:
    """The frame's own three axes, which is what ``f`` is almost always quoted along."""

    return np.eye(3, dtype=np.float64), tuple(str(label) for label in frame.axes)


def _resolve_directions(
    directions: ArrayLike | None,
    direction_labels: Sequence[str] | None,
    frame: ReferenceFrame,
) -> tuple[np.ndarray, tuple[str, ...]]:
    if directions is None:
        if direction_labels is not None:
            raise ValueError("direction_labels cannot be given without directions.")
        return _default_triad(frame)
    resolved = normalize_vectors(np.atleast_2d(np.asarray(directions, dtype=np.float64)))
    if resolved.ndim != 2 or resolved.shape[1] != 3:
        raise ValueError("directions must be an (m, 3) array of specimen-frame vectors.")
    if direction_labels is None:
        labels = tuple(f"d{index}" for index in range(resolved.shape[0]))
    else:
        labels = tuple(str(label) for label in direction_labels)
        if len(labels) != resolved.shape[0]:
            raise ValueError("direction_labels must have one label per direction.")
    return resolved, labels


def _is_orthonormal_triad(directions: np.ndarray) -> bool:
    if directions.shape != (3, 3):
        return False
    return bool(np.allclose(directions @ directions.T, np.eye(3), atol=1e-8))


def pole_orientation_tensor(
    directions: ArrayLike,
    weights: ArrayLike | None = None,
) -> np.ndarray:
    r"""The weighted second-moment tensor of a direction distribution.

    Computes :math:`\mathbf{A} = \langle \mathbf{v}\,\mathbf{v}^{\mathsf{T}}\rangle`.

    Purpose
    -------
    The object every Kearns route estimates. For basal poles it is the *Kearns
    tensor*: its quadratic form along any specimen direction is that direction's
    Kearns parameter, its trace is identically 1, and its eigenvalues are the
    Kearns parameters along the texture's own principal axes — the largest ``f``
    the specimen can exhibit in any direction, and the smallest.

    When to use
    -----------
    Directly, when ``f`` is wanted along several directions or when the
    principal axes of the basal-pole distribution matter (a pilgered tube's
    basal maxima are not generally on the RD/TD/ND axes). The
    ``kearns_from_*`` functions call it and wrap the result in a
    :class:`KearnsReport`.

    Method
    ------
    ``A = sum_i w_i v_i v_i^T / sum_i w_i`` with ``v_i`` normalized. The tensor
    is invariant under ``v -> -v``, so antipodal data need no folding and a pole
    figure measured on one hemisphere is as good as one measured on two.

    Parameters
    ----------
    directions : ArrayLike
        ``(n, 3)`` direction vectors; normalized internally, so magnitudes are
        ignored and only orientation matters.
    weights : ArrayLike, optional
        ``(n,)`` non-negative weights — volume fractions, pole densities times
        solid angles, or ODF weights. Uniform when omitted.

    Returns
    -------
    np.ndarray
        ``(3, 3)`` symmetric positive-semidefinite tensor of unit trace,
        read-only.

    Raises
    ------
    ValueError
        If the shapes disagree, if any weight is negative or non-finite, or if
        the weights sum to zero (nothing to average).

    See Also
    --------
    pytex.core.sphere.SphericalVectorSet.orientation_tensor : the unweighted
        form, for directional statistics of an unweighted direction set.

    Examples
    --------
    Three mutually perpendicular directions of equal weight give ``I/3``,
    whose quadratic form is ``1/3`` along every direction — the random-texture
    value.
    """

    values = normalize_vectors(np.atleast_2d(np.asarray(directions, dtype=np.float64)))
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("pole_orientation_tensor requires an (n, 3) array of directions.")
    if values.shape[0] == 0:
        raise ValueError("pole_orientation_tensor requires at least one direction.")
    if weights is None:
        weight_array = np.ones(values.shape[0], dtype=np.float64)
    else:
        weight_array = as_float_array(weights, shape=(values.shape[0],))
        if np.any(~np.isfinite(weight_array)) or np.any(weight_array < 0.0):
            raise ValueError("pole_orientation_tensor weights must be finite and non-negative.")
    total = float(np.sum(weight_array))
    if not total > 0.0:
        raise ValueError("pole_orientation_tensor weights must sum to a positive value.")
    tensor = np.einsum("n,ni,nj->ij", weight_array, values, values) / total
    tensor = 0.5 * (tensor + tensor.T)
    return freeze_array(np.ascontiguousarray(tensor, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class KearnsReport:
    """Kearns parameters along chosen specimen directions, with their evidence.

    Purpose
    -------
    A bare ``f`` is not comparable between laboratories: the route, the
    normalization, the measured tilt range, and the pole all change it by more
    than the differences that texture work cares about. This carries the number
    together with everything needed to judge it, and its :meth:`describe`
    renders that as prose a report can quote.

    Attributes
    ----------
    values : np.ndarray
        ``(m,)`` Kearns parameters, one per direction, each in ``[0, 1]``.
    directions : np.ndarray
        ``(m, 3)`` unit specimen directions the parameters refer to.
    direction_labels : tuple of str
        Names for those directions — the specimen frame's axis labels by
        default, so a report reads ``f_RD`` rather than ``f_0``.
    method : KearnsMethod
        Which route produced the values.
    pole : CrystalPlane
        The pole resolved along the directions — ``(0001)`` for the classical
        Kearns parameter. Recorded because the same machinery answers the
        question for any pole, and a non-basal answer must not be mistaken for
        the Kearns parameter.
    specimen_frame : ReferenceFrame
    orientation_tensor : np.ndarray, optional
        The ``(3, 3)`` tensor the values came from, when the route determines
        it. The diffractogram and tilt-profile routes determine only one
        diagonal element per measured section and leave this ``None``.
    diagnostics : Mapping[str, float]
        Route-specific numbers needed to judge the result — measured solid-angle
        coverage, the maximum tilt reached, the kernel shrinkage factor removed,
        the spread of diffraction-vector directions, and so on.
    notes : tuple of str
        Warnings and convention statements that belong with the number.
    provenance : ProvenanceRecord, optional

    See Also
    --------
    :doc:`../../theory/kearns_parameter_and_basal_pole_texture`
    """

    values: np.ndarray
    directions: np.ndarray
    direction_labels: tuple[str, ...]
    method: KearnsMethod
    pole: CrystalPlane
    specimen_frame: ReferenceFrame
    orientation_tensor: np.ndarray | None = None
    diagnostics: Mapping[str, float] = field(default_factory=dict)
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        directions = normalize_vectors(np.atleast_2d(np.asarray(self.directions, dtype=np.float64)))
        values = as_float_array(self.values, shape=(directions.shape[0],))
        if np.any(~np.isfinite(values)):
            raise ValueError("KearnsReport values must be finite.")
        if np.any(values < -1e-9) or np.any(values > 1.0 + 1e-9):
            raise ValueError(
                "KearnsReport values must lie in [0, 1]; f is a mean of cos^2 and cannot leave it."
            )
        if len(self.direction_labels) != directions.shape[0]:
            raise ValueError("KearnsReport needs one direction label per direction.")
        _require_specimen_frame(self.specimen_frame, argument="KearnsReport.specimen_frame")
        if self.orientation_tensor is not None:
            tensor = as_float_array(self.orientation_tensor, shape=(3, 3))
            if not np.isclose(float(np.trace(tensor)), 1.0, atol=1e-8):
                raise ValueError(
                    "KearnsReport.orientation_tensor must have unit trace; "
                    f"got {float(np.trace(tensor)):.6f}."
                )
            object.__setattr__(self, "orientation_tensor", freeze_array(tensor))
        object.__setattr__(self, "directions", freeze_array(directions))
        object.__setattr__(self, "values", freeze_array(values))
        object.__setattr__(self, "direction_labels", tuple(str(x) for x in self.direction_labels))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))
        object.__setattr__(self, "notes", tuple(str(note) for note in self.notes))

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def value(self, label: str) -> float:
        """The Kearns parameter along the direction named ``label``.

        Lets a caller ask for ``f_ND`` by name instead of by position, which is
        how the number is quoted in specifications.
        """

        try:
            index = self.direction_labels.index(label)
        except ValueError:
            available = ", ".join(self.direction_labels)
            raise KeyError(f"No direction labelled '{label}'. Available: {available}.") from None
        return float(self.values[index])

    @property
    def is_orthonormal_triad(self) -> bool:
        """Whether the directions form an orthonormal triad, so the sum rule applies."""

        return _is_orthonormal_triad(np.asarray(self.directions, dtype=np.float64))

    @property
    def triad_sum(self) -> float | None:
        """Sum of the three values when the directions are an orthonormal triad.

        The single most informative check on a Kearns measurement: the sum is
        identically 1 for any texture, so a departure measures the systematic
        error of the *measurement* — unmeasured tilt range, a wrong random
        standard, an unbalanced background — and nothing about the material.
        ``None`` when the directions are not a triad, where no such rule exists.
        """

        if not self.is_orthonormal_triad:
            return None
        return float(np.sum(self.values))

    @property
    def principal_values(self) -> np.ndarray | None:
        """Eigenvalues of the orientation tensor, ascending; ``None`` without one.

        The Kearns parameters along the texture's own principal axes: the
        largest is the greatest ``f`` obtainable in any direction and the
        smallest the least, so together they bound every direction's value.
        """

        if self.orientation_tensor is None:
            return None
        eigenvalues = np.linalg.eigvalsh(np.asarray(self.orientation_tensor, dtype=np.float64))
        return freeze_array(np.ascontiguousarray(eigenvalues, dtype=np.float64))

    def to_json(self) -> dict[str, Any]:
        """A JSON-serializable dictionary of the report.

        Kept in lockstep with :meth:`describe`: every number the prose states
        appears here under the same name, so a machine-read result and a
        human-read one cannot disagree.
        """

        payload: dict[str, Any] = {
            "method": self.method,
            "pole": _plane_label(self.pole),
            "phase": self.pole.phase.name,
            "specimen_frame": self.specimen_frame.name,
            "directions": [
                {
                    "label": label,
                    "vector": [float(component) for component in vector],
                    "f": float(value),
                }
                for label, vector, value in zip(
                    self.direction_labels, self.directions, self.values, strict=True
                )
            ],
            "isotropic_value": KEARNS_ISOTROPIC_VALUE,
            "triad_sum": self.triad_sum,
            "diagnostics": {key: float(value) for key, value in self.diagnostics.items()},
            "notes": list(self.notes),
        }
        if self.orientation_tensor is not None:
            payload["orientation_tensor"] = [
                [float(component) for component in row] for row in self.orientation_tensor
            ]
            principal = self.principal_values
            assert principal is not None
            payload["principal_values"] = [float(value) for value in principal]
        return payload

    def to_json_string(self, *, indent: int = 2) -> str:
        """The JSON contract of :meth:`to_json` rendered as text."""

        return json.dumps(self.to_json(), indent=indent, sort_keys=False)

    def describe(self) -> str:
        """Convention-explicit prose describing the result and its reliability.

        Purpose
        -------
        A Kearns parameter is quoted in specifications and compared between
        laboratories, so the number alone is not a result. This states the
        route, the pole, the values against the ``1/3`` random reference, the
        triad sum that tests the measurement, and every diagnostic and caveat
        the route carries.
        """

        method_prose = {
            "orientations": (
                "resolved directly from discrete orientations, with no binning, "
                "interpolation, or tilt truncation"
            ),
            "odf": "resolved from a fitted orientation distribution function",
            "pole_figure": (
                "integrated over a pole figure after Baron et al. (1990) Eq. (5), "
                "with solid-angle weights"
            ),
            "tilt_profile": (
                "integrated over a basal-pole tilt profile after Kearns (1965) Eq. (5)"
            ),
            "diffractogram": (
                "estimated from relative diffraction-peak intensities after "
                "Kearns (1965), the inverse-pole-figure route"
            ),
        }[self.method]
        lines = [
            f"Kearns parameter f for the {_plane_label(self.pole)} pole of "
            f"{self.pole.phase.name}, {method_prose}.",
            "",
            "f is the volume-weighted mean of cos^2(angle between the pole and the "
            "direction): the effective fraction of poles aligned with that direction. "
            f"A random texture gives {KEARNS_ISOTROPIC_VALUE:.4f} in every direction.",
            "",
        ]
        for label, value in zip(self.direction_labels, self.values, strict=True):
            comparison = (
                "at the random value"
                if abs(value - KEARNS_ISOTROPIC_VALUE) < 5e-3
                else (
                    f"{value / KEARNS_ISOTROPIC_VALUE:.2f} times random"
                    if value > KEARNS_ISOTROPIC_VALUE
                    else f"{value / KEARNS_ISOTROPIC_VALUE:.2f} times random (depleted)"
                )
            )
            lines.append(f"  f_{label} = {value:.4f}  ({comparison})")
        triad_sum = self.triad_sum
        if triad_sum is not None:
            deviation = abs(triad_sum - 1.0)
            verdict = (
                "consistent with the exact sum rule"
                if deviation <= 0.02
                else "a systematic error of the measurement, not a property of the material"
            )
            lines += [
                "",
                f"Triad sum = {triad_sum:.4f}. The sum over any orthonormal specimen "
                f"triad is identically 1 for every texture, so the departure of "
                f"{deviation:.4f} is {verdict}.",
            ]
        principal = self.principal_values
        if principal is not None:
            lines += [
                "",
                "Principal Kearns values (eigenvalues of the pole orientation tensor, "
                f"ascending): {principal[0]:.4f}, {principal[1]:.4f}, {principal[2]:.4f}. "
                "No direction can give an f outside this range.",
            ]
        if self.diagnostics:
            lines += ["", "Diagnostics:"]
            lines += [f"  {key} = {value:.6g}" for key, value in sorted(self.diagnostics.items())]
        if self.notes:
            lines += ["", "Notes:"]
            lines += [f"  - {note}" for note in self.notes]
        return "\n".join(lines)


def _pole_directions_from_orientations(
    orientations: OrientationSet,
    pole: CrystalPlane,
    *,
    weights: ArrayLike | None,
    include_symmetry_family: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Specimen-frame pole directions and their weights, one row per (grain, family member)."""

    if orientations.crystal_frame != pole.phase.crystal_frame:
        raise ValueError("The orientations must use the pole phase's crystal frame.")
    base_weights = (
        np.ones(len(orientations), dtype=np.float64)
        if weights is None
        else as_float_array(weights, shape=(len(orientations),))
    )
    if np.any(~np.isfinite(base_weights)) or np.any(base_weights < 0.0):
        raise ValueError("Orientation weights must be finite and non-negative.")
    normals = (
        pole.phase.symmetry.equivalent_vectors(pole.normal)
        if include_symmetry_family
        else pole.normal[None, :]
    )
    mapped = []
    for normal in normals:
        directions = orientations.map_crystal_directions(normal)
        mapped.append(np.asarray(getattr(directions, "values", directions), dtype=np.float64))
    stacked = np.vstack(mapped)
    tiled = np.tile(base_weights, len(normals))
    return stacked, tiled


def kearns_from_orientations(
    orientations: OrientationSet,
    *,
    pole: CrystalPlane,
    weights: ArrayLike | None = None,
    directions: ArrayLike | None = None,
    direction_labels: Sequence[str] | None = None,
    include_symmetry_family: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> KearnsReport:
    """Kearns parameters from discrete orientations — the exact route.

    Purpose
    -------
    Evaluate ``f`` directly from an orientation set: EBSD data, a simulated
    texture, or the support of a model. This is the reference implementation
    against which the diffraction routes are judged, because it involves no
    binning, no interpolation, no tilt truncation, and no normalization
    assumption — the definition ``f = <cos^2 phi>`` is evaluated as written.

    When to use
    -----------
    For EBSD-derived Kearns parameters, which Mani Krishna *et al.* (2011) found
    the most consistent route for recrystallized microstructures; and as the
    ground truth whenever a pole-figure or diffractogram implementation is being
    validated on a simulated texture.

    Method
    ------
    Every orientation maps the pole normal into the specimen frame; the weighted
    second-moment tensor of those directions is
    :func:`pole_orientation_tensor`, and ``f_d = d^T A d``. With
    ``include_symmetry_family``, the whole ``{hkl}`` orbit is mapped and
    averaged, which is what a diffraction measurement sees. For ``(0001)`` in
    ``6/mmm`` the orbit is ``+/-[0001]`` and the flag makes no difference, since
    ``c c^T`` is invariant under sign.

    Parameters
    ----------
    orientations : OrientationSet
        Crystal-to-specimen orientations, in the pole phase's crystal frame.
    pole : CrystalPlane
        The pole to resolve — ``(0001)`` for the Kearns parameter proper.
    weights : ArrayLike, optional
        One non-negative weight per orientation: grain area or volume for EBSD
        data, where an unweighted mean would count a large grain and a small one
        equally. Uniform when omitted.
    directions : ArrayLike, optional
        ``(m, 3)`` specimen directions to resolve along. Defaults to the
        specimen frame's own three axes.
    direction_labels : Sequence[str], optional
        Names for those directions; required to be one per direction when
        ``directions`` is given.
    include_symmetry_family : bool
        Average over the pole's symmetry orbit (default), which is what a
        diffraction measurement of ``{hkl}`` records.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KearnsReport
        Carrying the full orientation tensor, so the triad sum is exactly 1 up
        to floating point and any other direction can be queried.

    Raises
    ------
    ValueError
        If the orientation set does not use the pole phase's crystal frame, if
        weights are negative or misshapen, or if the direction labels do not
        match the directions.

    Examples
    --------
    A uniformly random orientation set gives ``f = 1/3`` along every direction
    to within sampling error, and its triad sum is exactly 1.

    See Also
    --------
    kearns_from_odf : the same quantity from a fitted distribution, with the
        kernel-shrinkage correction that a fitted ODF needs.
    """

    frame = _require_specimen_frame(
        orientations.specimen_frame, argument="orientations.specimen_frame"
    )
    resolved, labels = _resolve_directions(directions, direction_labels, frame)
    pole_directions, pole_weights = _pole_directions_from_orientations(
        orientations,
        pole,
        weights=weights,
        include_symmetry_family=include_symmetry_family,
    )
    tensor = pole_orientation_tensor(pole_directions, pole_weights)
    values = np.einsum("mi,ij,mj->m", resolved, tensor, resolved)
    return KearnsReport(
        values=values,
        directions=resolved,
        direction_labels=labels,
        method="orientations",
        pole=pole,
        specimen_frame=frame,
        orientation_tensor=tensor,
        diagnostics={
            "orientation_count": float(len(orientations)),
            "pole_count": float(pole_directions.shape[0]),
        },
        notes=(
            "Exact evaluation of the definition: no binning, interpolation, or tilt "
            "truncation enters, so the only error is the orientation statistics.",
        ),
        provenance=orientations.provenance if provenance is None else provenance,
    )


def kernel_axis_shrinkage(kernel: KernelSpec, *, quadrature_points: int = 4096) -> float:
    r"""How much an ODF kernel shrinks a pole orientation tensor toward isotropy.

    Purpose
    -------
    A fitted ODF is a *smoothed* estimate: convolving the true distribution with
    a kernel of halfwidth ``h`` spreads every crystal's pole over a cone, which
    pulls the pole orientation tensor — and therefore every Kearns parameter —
    toward the isotropic ``1/3``. The bias is not small: at a 15 degree
    halfwidth it removes several percent of the departure from random, and it is
    a *bias*, so more data does not reduce it. This returns the factor that
    quantifies it, so it can be divided out exactly rather than tolerated.

    Method
    ------
    Write the kernel-smeared pole direction as ``v = R c`` with ``R`` a rotation
    of angle ``omega`` about a uniformly distributed axis. Rodrigues' formula
    gives ``c . R c = cos(omega) + t (1 - cos(omega))`` with ``t = (a . c)^2``,
    and for a uniform axis ``E[t] = 1/3``, ``E[t^2] = 1/5``. Averaging
    ``(c . Rc)^2`` over the axis and then over the kernel's angular density on
    SO(3) — the kernel value times the Haar factor ``(1 - cos omega)/pi`` —
    gives

    .. math::

        \rho = \left\langle \cos^{2}\beta \right\rangle, \qquad
        \mathbf{A}_{\text{smoothed}} = \rho\,\mathbf{A}
        + \tfrac{1-\rho}{2}\left(\mathbf{I} - \mathbf{A}\right),

    so that ``f_smoothed = 1/3 + beta (f - 1/3)`` with the *shrinkage factor*
    ``beta = (3 rho - 1) / 2`` returned by this function's caller through
    :func:`kearns_from_odf`. The result is independent of crystal symmetry:
    every symmetry operation conjugates the kernel into an equally isotropic
    one, so the symmetrized kernel shrinks the tensor by exactly the same
    factor.

    Parameters
    ----------
    kernel : KernelSpec
        The kernel and halfwidth the ODF was estimated with.
    quadrature_points : int
        Nodes of the Simpson quadrature over ``omega`` in ``[0, pi]``. The
        default is far more than the smooth integrand needs.

    Returns
    -------
    float
        ``rho``: 1 for a zero-halfwidth kernel, falling toward ``1/3`` as the
        kernel approaches uniform, where all orientation information is gone.
        At usable halfwidths it stays well above ``1/3`` — 0.958 at 10 degrees,
        0.718 at 30 for the de la Vallee Poussin kernel. Beyond about 90 degrees
        it can fall to or below ``1/3``, which is not a regime to correct in but
        a sign the halfwidth is meaningless; :func:`kearns_from_odf` raises
        rather than dividing by the resulting non-positive shrinkage.

    Examples
    --------
    A very narrow kernel returns essentially 1 (no shrinkage); a kernel wide
    enough to be effectively uniform returns ``1/3``, at which point the
    correction is undefined because no texture survives the smoothing.

    See Also
    --------
    kearns_from_odf : applies the correction.
    """

    if quadrature_points < 8:
        raise ValueError("kernel_axis_shrinkage needs at least 8 quadrature points.")
    # A narrow kernel concentrates the integrand into the first few degrees, and
    # a grid uniform over [0, pi] cannot resolve it: at a 0.5 degree halfwidth a
    # 4096-point uniform grid misses the peak entirely. Refine a panel scaled to
    # the halfwidth and merge it with the coarse one, so accuracy does not
    # depend on the halfwidth being large.
    inner_bound = float(min(np.pi, 12.0 * np.deg2rad(kernel.halfwidth_deg)))
    omega = np.union1d(
        np.linspace(0.0, np.pi, quadrature_points),
        np.linspace(0.0, inner_bound, quadrature_points),
    )
    cos_omega = np.cos(omega)
    complement = 1.0 - cos_omega
    # <cos^2(beta)> over a uniformly distributed rotation axis, at fixed omega.
    axis_average = cos_omega**2 + (2.0 / 3.0) * cos_omega * complement + complement**2 / 5.0
    # The SO(3) Haar density in the rotation angle, (1 - cos omega)/pi.
    measure = np.asarray(kernel.evaluate(omega), dtype=np.float64) * complement / np.pi
    numerator = float(np.trapezoid(measure * axis_average, omega))
    denominator = float(np.trapezoid(measure, omega))
    if not denominator > 0.0:
        raise ValueError("The kernel integrates to zero over SO(3); its halfwidth is degenerate.")
    return float(numerator / denominator)


def kearns_from_odf(
    odf: ODF,
    *,
    pole: CrystalPlane,
    directions: ArrayLike | None = None,
    direction_labels: Sequence[str] | None = None,
    include_symmetry_family: bool = True,
    deconvolve_kernel: bool = False,
    provenance: ProvenanceRecord | None = None,
) -> KearnsReport:
    """Kearns parameters from a fitted orientation distribution function.

    Purpose
    -------
    The route Mani Krishna *et al.* (2011) call the ODF method: reconstruct an
    ODF — typically by inverting several incomplete pole figures — and resolve
    the basal poles of its orientation weights onto the specimen directions.
    Its advantage over the pole-figure route is that it does not require a
    strong ``(0002)`` peak in the measured section: alternative reflections plus
    the inversion supply the basal information, which is why it stays usable in
    an ND-TD section of strongly basal-textured zirconium where ``(0002)``
    intensity is negligible.

    Method
    ------
    A discrete ODF holds two distinguishable objects, and they have different
    Kearns parameters:

    - the **support tensor** ``A_support``, the weighted second moment of the
      support orientations' poles — what :func:`kearns_from_orientations`
      returns for the same weighted set;
    - the **density tensor** ``A_density``, the second moment under the
      continuous density ``sum_j w_j K(g; g_j)`` that the ODF object *is*.

    Convolving with the kernel shrinks every departure from isotropy by the
    closed-form factor ``beta = (3 rho - 1)/2`` from
    :func:`kernel_axis_shrinkage`, so the two are related exactly by

    ``A_density = I/3 + beta (A_support - I/3)``,

    with no numerical deconvolution needed in either direction. Which one to
    report depends on where the ODF came from, so it is a parameter rather than
    a hidden choice: see ``deconvolve_kernel``.

    Parameters
    ----------
    odf : ODF
        A discrete kernel-density ODF.
    pole : CrystalPlane
        The pole to resolve — ``(0001)`` for the Kearns parameter.
    directions, direction_labels
        As in :func:`kearns_from_orientations`.
    include_symmetry_family : bool
        Average over the pole's symmetry orbit (default).
    deconvolve_kernel : bool
        ``False`` (default) reports the density tensor: the ``f`` of the
        distribution this ODF object represents, which is what binning the ODF
        and summing ``V cos^2 phi`` gives, and the right choice for an ODF
        fitted by pole-figure inversion — there the weights were chosen so that
        the *smoothed* density matches the measurement, so the smoothed density
        is the model of the material.

        ``True`` reports the support tensor instead, undoing the smoothing. Use
        it when the kernel is estimation blur rather than part of the model —
        an ODF built by :meth:`~pytex.texture.ODF.from_orientations` from
        measured EBSD orientations, where the support *is* the data and the
        kernel only makes it continuous.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KearnsReport
        With ``kernel_halfwidth_deg`` and ``kernel_shrinkage_factor``
        diagnostics, so the size of the difference between the two readings is
        visible whichever was chosen.

    Raises
    ------
    ValueError
        If the ODF support does not use the pole phase's crystal frame.

    Examples
    --------
    For an ODF built from a single orientation with a 10 degree kernel, the
    density reading gives ``f = 0.937`` along the basal pole and the support
    reading gives exactly 1. The whole of that gap is the kernel, not the
    material — which is why an ODF's halfwidth must be reported alongside any
    ``f`` taken from it.

    See Also
    --------
    kernel_axis_shrinkage : the shrinkage factor and its derivation.
    """

    frame = _require_specimen_frame(
        odf.orientations.specimen_frame, argument="odf.orientations.specimen_frame"
    )
    resolved, labels = _resolve_directions(directions, direction_labels, frame)
    pole_directions, pole_weights = _pole_directions_from_orientations(
        odf.orientations,
        pole,
        weights=odf.normalized_weights,
        include_symmetry_family=include_symmetry_family,
    )
    support_tensor = np.asarray(
        pole_orientation_tensor(pole_directions, pole_weights), dtype=np.float64
    )
    rho = kernel_axis_shrinkage(odf.kernel)
    shrinkage = (3.0 * rho - 1.0) / 2.0
    notes: list[str] = []
    if deconvolve_kernel:
        tensor = support_tensor
        notes.append(
            "Reported from the ODF's support orientations, treating the "
            f"{odf.kernel.halfwidth_deg:g} degree kernel as estimation blur rather than "
            "part of the model. The density this ODF represents would give values "
            f"shrunk toward 1/3 by a factor {shrinkage:.4f}."
        )
    else:
        tensor = _ISOTROPIC_TENSOR + shrinkage * (support_tensor - _ISOTROPIC_TENSOR)
        notes.append(
            "Reported from the continuous density this ODF represents, whose "
            f"{odf.kernel.halfwidth_deg:g} degree kernel shrinks every departure from 1/3 "
            f"by a factor {shrinkage:.4f}. Pass deconvolve_kernel=True for the support's "
            "own value, which is the right reading when the ODF was estimated from "
            "measured orientations."
        )
    tensor = 0.5 * (tensor + tensor.T)
    values = np.clip(np.einsum("mi,ij,mj->m", resolved, tensor, resolved), 0.0, 1.0)
    return KearnsReport(
        values=values,
        directions=resolved,
        direction_labels=labels,
        method="odf",
        pole=pole,
        specimen_frame=frame,
        orientation_tensor=tensor,
        diagnostics={
            "support_size": float(len(odf.orientations)),
            "kernel_halfwidth_deg": float(odf.kernel.halfwidth_deg),
            "kernel_shrinkage_factor": shrinkage,
        },
        notes=tuple(notes),
        provenance=odf.provenance if provenance is None else provenance,
    )


def kearns_from_pole_figure(
    pole_figure: PoleFigure,
    *,
    directions: ArrayLike | None = None,
    direction_labels: Sequence[str] | None = None,
    integration_weights: ArrayLike | None = None,
    provenance: ProvenanceRecord | None = None,
) -> KearnsReport:
    r"""Kearns parameters by integrating a pole figure — Baron *et al.* Eq. (5).

    Purpose
    -------
    The standard industrial route: measure a ``(0002)`` pole figure and
    integrate it. Baron *et al.* (1990) define

    .. math::

        f_{i} = \frac{\int_{0}^{\alpha_{\max}}\!\!\int_{0}^{2\pi}
                 P(\alpha,\beta)\cos^{2}\alpha_{i}\,\sin\alpha\,
                 \mathrm{d}\alpha\,\mathrm{d}\beta}
                {\int_{0}^{\alpha_{\max}}\!\!\int_{0}^{2\pi}
                 P(\alpha,\beta)\,\sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}

    with ``alpha_i`` the angle between the pole direction and specimen direction
    ``i``. Dividing by the measured integral rather than by ``2 pi`` is the
    *pseudo-norm* of Kern and Bergmann: it makes an incomplete figure usable,
    at the price of assuming the unmeasured cap has the same mean as the
    measured one.

    When to use
    -----------
    Whenever a basal pole figure exists and the section carries real ``(0002)``
    intensity. It is *not* usable in a section where the basal peak is
    negligible — the ND-TD section of strongly basal-textured zirconium — where
    the normalization divides by noise; Mani Krishna *et al.* (2011) traced the
    inconsistent ``f_RD`` values from that section to exactly this. Use
    :func:`kearns_from_odf` there instead.

    Method
    ------
    ``f_d = d^T A d`` with ``A`` the solid-angle-weighted second-moment tensor
    of the measured pole directions, which is the ratio above written as a
    tensor. The weights follow the figure's ``sampling`` attribute:

    - ``"scattered_poles"`` — each row is one pole carrying its own weight, so
      the intensities *are* the weights and no quadrature is involved.
    - ``"sampled_density"`` — the rows are densities sampled on a raster, so
      each is multiplied by its solid angle from
      :func:`~pytex.core.sphere.raster_solid_angle_weights` before the sum.

    Using the wrong one over-counts the pole of a tilt raster by up to 50
    percent, which is why the reading is taken from the figure rather than
    assumed.

    Parameters
    ----------
    pole_figure : PoleFigure
        The measured or computed figure. Background subtraction and the
        defocusing correction belong *before* this call; see
        :class:`~pytex.texture.PoleFigureCorrectionSpec`.
    directions, direction_labels
        As in :func:`kearns_from_orientations`.
    integration_weights : ArrayLike, optional
        Explicit per-point solid-angle weights, one per sampled direction.
        Supply these when the support is neither a pole cloud nor a regular
        tilt raster — for instance an :class:`~pytex.core.sphere.S2Grid`, whose
        own ``weights`` are the correct choice.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KearnsReport
        With ``max_polar_deg`` and ``measured_solid_angle_fraction``
        diagnostics. The latter is the fraction of the hemisphere actually
        integrated; a figure truncated at 75 degrees covers only 74 percent of
        it, and the note records that the missing cap is assumed to resemble
        the measured one.

    Raises
    ------
    ValueError
        If explicit weights do not match the sampled directions, or if the
        intensities are everywhere zero, leaving nothing to normalize.

    Examples
    --------
    A pole figure computed from a random texture integrates to ``1/3`` in every
    direction and to a triad sum of exactly 1 — the calibration of this
    function.

    See Also
    --------
    pytex.core.sphere.raster_solid_angle_weights : the quadrature weights.
    """

    frame = _require_specimen_frame(
        pole_figure.specimen_frame, argument="pole_figure.specimen_frame"
    )
    resolved, labels = _resolve_directions(directions, direction_labels, frame)
    sample_directions = np.asarray(pole_figure.sample_directions, dtype=np.float64)
    intensities = np.asarray(pole_figure.intensities, dtype=np.float64)
    polar_deg, _ = directions_to_spherical_angles(sample_directions)
    polar_deg = np.asarray(polar_deg, dtype=np.float64).reshape(-1)
    if pole_figure.antipodal:
        # Antipodal data may be recorded on either hemisphere; v v^T does not
        # care, but the ring structure the quadrature groups by does.
        polar_deg = np.minimum(polar_deg, 180.0 - polar_deg)

    notes: list[str] = []
    if integration_weights is not None:
        quadrature = as_float_array(integration_weights, shape=(intensities.shape[0],))
        if np.any(~np.isfinite(quadrature)) or np.any(quadrature < 0.0):
            raise ValueError("integration_weights must be finite and non-negative.")
        weights = intensities * quadrature
        notes.append("Integration used the explicitly supplied per-point solid-angle weights.")
    elif pole_figure.sampling == "scattered_poles":
        weights = intensities
        notes.append(
            "The figure records scattered poles, so its intensities are per-pole weights "
            "and no solid-angle quadrature is applied."
        )
    else:
        # An antipodal figure lives on the hemisphere, so its outermost ring
        # owns no solid angle beyond the equator; saying so keeps the quadrature
        # error at the equator from dominating the integral.
        quadrature = np.asarray(
            raster_solid_angle_weights(
                polar_deg,
                polar_max_deg=90.0 if pole_figure.antipodal else None,
            ),
            dtype=np.float64,
        )
        weights = intensities * quadrature
        notes.append(
            "Sampled densities were weighted by the solid angle of each tilt ring; an "
            "unweighted mean over a tilt raster over-counts the figure centre."
        )
    if not float(np.sum(weights)) > 0.0:
        raise ValueError(
            "The pole figure carries no positive weight over its measured region, so f is "
            "undefined. In a section where the measured pole has negligible intensity, use "
            "kearns_from_odf instead."
        )

    tensor = pole_orientation_tensor(sample_directions, weights)
    values = np.einsum("mi,ij,mj->m", resolved, tensor, resolved)

    max_polar = float(np.max(polar_deg))
    covered = float(1.0 - np.cos(np.deg2rad(min(max_polar, 90.0))))
    if max_polar < 89.0:
        notes.append(
            f"The figure reaches only {max_polar:.1f} degrees of tilt, covering "
            f"{covered:.3f} of the hemisphere. Following Kern and Bergmann's pseudo-norm, "
            "the result normalizes over the measured cap alone and therefore assumes the "
            "unmeasured cap resembles it. Averaging two perpendicular sections, as Mani "
            "Krishna et al. (2011) recommend, is the usual mitigation."
        )
    return KearnsReport(
        values=values,
        directions=resolved,
        direction_labels=labels,
        method="pole_figure",
        pole=pole_figure.pole,
        specimen_frame=frame,
        orientation_tensor=tensor,
        diagnostics={
            "sampled_point_count": float(sample_directions.shape[0]),
            "max_polar_deg": max_polar,
            "measured_solid_angle_fraction": covered,
        },
        notes=tuple(notes),
        provenance=pole_figure.provenance if provenance is None else provenance,
    )


def kearns_from_tilt_profile(
    polar_deg: ArrayLike,
    intensity: ArrayLike,
    *,
    pole: CrystalPlane,
    specimen_frame: ReferenceFrame,
    direction: ArrayLike | None = None,
    direction_label: str | None = None,
    provenance: ProvenanceRecord | None = None,
) -> KearnsReport:
    r"""Kearns parameter from a basal-pole tilt profile — Kearns (1965) Eq. (5).

    Purpose
    -------
    The innermost calculation of the original method, exposed on its own.
    Kearns showed that ``f`` needs only the pole density averaged over the full
    360 degrees of azimuth, ``I(phi)``, as a function of tilt ``phi`` from the
    reference direction — the azimuthal detail of a pole figure is irrelevant to
    a second-rank property. Then

    .. math::

        f = \int_{0}^{\pi/2} I(\phi)\,\sin\phi\,\cos^{2}\phi\,\mathrm{d}\phi
        \Big/ \int_{0}^{\pi/2} I(\phi)\,\sin\phi\,\mathrm{d}\phi .

    When to use
    -----------
    When the tilt profile is what is available: a pole-figure ring average, a
    published ``I(phi)`` curve, or the output of :func:`basal_tilt_profile`. Use
    it also to reproduce a literature calculation, since it is exactly the
    tabulated procedure of Kearns' Table 3.

    Method
    ------
    The ``sin phi`` factor converts pole density to volume fraction — the ring
    of orientations at tilt ``phi`` has circumference proportional to
    ``sin phi``, which is why the profile at high tilt matters even where the
    density is small, and why ``V`` vanishes at ``phi = 0`` however intense the
    pole is there. The quotient makes the result independent of the profile's
    scale, so ``I`` may be in times-random, counts, or anything proportional.
    Samples are treated as bin representatives, as in Kearns' Table 3, so the
    common bin width cancels; the nodes must therefore be equally spaced.

    Parameters
    ----------
    polar_deg : ArrayLike
        ``(n,)`` tilt angles in degrees, within ``[0, 90]``, equally spaced.
        Bin midpoints (5, 15, ... 85 for 10 degree bins) are the classical
        choice and avoid the endpoints, where the integrand vanishes or the
        measurement is worst.
    intensity : ArrayLike
        ``(n,)`` azimuthally averaged pole densities at those tilts, in any
        units proportional to density.
    pole : CrystalPlane
        The pole the profile describes.
    specimen_frame : ReferenceFrame
    direction : ArrayLike, optional
        The specimen direction the tilts are measured from. Defaults to the
        frame's third axis, which is the section normal of a reflection
        measurement.
    direction_label : str, optional
        Its name; defaults to that axis's label.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KearnsReport
        With no orientation tensor: one profile determines one direction's
        value, not the whole tensor, so no triad sum is available.

    Raises
    ------
    ValueError
        If the tilts leave ``[0, 90]``, are not equally spaced, if any intensity
        is negative, or if the profile is everywhere zero.

    Examples
    --------
    The longitudinal-section profile of Kearns' Table 3 — intensities
    ``3.27, 2.71, 1.69, 1.35, 1.17, 0.97, 0.73, 0.62, 0.55`` at the midpoints of
    10 degree bins — returns 0.4879 against his tabulated 0.488.

    See Also
    --------
    kearns_from_diffractogram : builds the profile from peak intensities and
        calls this.
    """

    tilts = as_float_array(np.asarray(polar_deg, dtype=np.float64).reshape(-1), shape=(None,))
    values = as_float_array(
        np.asarray(intensity, dtype=np.float64).reshape(-1), shape=(tilts.shape[0],)
    )
    if tilts.size < 2:
        raise ValueError("A tilt profile needs at least two samples.")
    if np.any(tilts < -1e-9) or np.any(tilts > 90.0 + 1e-9):
        raise ValueError("Tilt angles must lie in [0, 90] degrees.")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Tilt-profile intensities must be finite and non-negative.")
    order = np.argsort(tilts)
    tilts = tilts[order]
    values = values[order]
    steps = np.diff(tilts)
    if not np.allclose(steps, steps[0], atol=1e-6):
        raise ValueError(
            "kearns_from_tilt_profile expects equally spaced tilt nodes, because Kearns' "
            "procedure treats each sample as a bin representative and lets the common bin "
            "width cancel. Resample the profile onto a uniform grid first."
        )
    radians = np.deg2rad(tilts)
    volume = values * np.sin(radians)
    total = float(np.sum(volume))
    if not total > 0.0:
        raise ValueError("The tilt profile is everywhere zero, so f is undefined.")
    f_value = float(np.sum(volume * np.cos(radians) ** 2) / total)

    frame = _require_specimen_frame(specimen_frame, argument="specimen_frame")
    if direction is None:
        axis = np.array([0.0, 0.0, 1.0])
        label = frame.axes[2] if direction_label is None else direction_label
    else:
        axis = normalize_vectors(np.atleast_2d(np.asarray(direction, dtype=np.float64)))[0]
        label = "d0" if direction_label is None else direction_label
    return KearnsReport(
        values=np.array([f_value]),
        directions=axis[None, :],
        direction_labels=(str(label),),
        method="tilt_profile",
        pole=pole,
        specimen_frame=frame,
        orientation_tensor=None,
        diagnostics={
            "profile_node_count": float(tilts.size),
            "tilt_step_deg": float(steps[0]),
            "max_tilt_deg": float(tilts[-1]),
        },
        notes=(
            "One tilt profile determines f along one direction only; the orientation "
            "tensor, and hence the triad sum, needs a profile per section.",
        ),
        provenance=provenance,
    )


def basal_tilt_angle_deg(plane: CrystalPlane) -> float:
    """Angle in degrees between a plane's normal and the crystal's ``c`` axis.

    Purpose
    -------
    The quantity that makes the diffractogram route work: a reflection observed
    in a symmetric scan comes from grains whose ``(hkil)`` normal is along the
    section normal, and those grains therefore have their basal pole at this
    fixed angle to it. Kearns tabulated these angles by hand for alpha-Zr; here
    they are computed from the phase's own reciprocal metric, so they follow the
    lattice parameters used rather than a transcribed constant.

    Parameters
    ----------
    plane : CrystalPlane
        Any plane of a hexagonal (or other uniaxial) phase.

    Returns
    -------
    float
        The angle in ``[0, 90]`` degrees, with antipodal equivalence, so
        ``(0002)`` gives 0 and ``(10-10)`` gives 90.

    Examples
    --------
    For alpha-Zr with ``c/a = 1.593``, ``(10-11)`` gives 61.5 degrees, matching
    the 61.4 of Kearns' Table 2.
    """

    phase: Phase = plane.phase
    c_axis = phase.lattice.direct_basis().matrix @ np.array([0.0, 0.0, 1.0])
    c_axis = c_axis / float(np.linalg.norm(c_axis))
    normal = np.asarray(plane.normal, dtype=np.float64)
    cosine = float(np.clip(abs(float(np.dot(normal, c_axis))), 0.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def harris_texture_coefficients(intensities: ArrayLike) -> np.ndarray:
    """Normalize relative pole densities to a mean of one — Kearns' Eqs. (17)-(19).

    Purpose
    -------
    The diffractogram route needs pole densities in times-random units, which
    requires a random-powder standard measured under identical conditions. When
    none exists, the Harris texture coefficient substitutes an assumption for
    the standard: that the *mean* of ``I/I0`` over the measured reflections is
    one. Dividing by that mean then removes the unknown proportionality constant
    that relates the observed and random intensity scales.

    When to use
    -----------
    Only when no random standard is available. Kearns measured one and showed
    the assumption fails: over his samples the mean of ``I/I0`` ran from 1.02 to
    1.56, averaging 1.23, so texture coefficients underestimated the true pole
    densities by about 23 percent. Mani Krishna *et al.* (2011) found the same
    at the level of ``f``, which did not sum to 1 over the three sections until
    normalized. Prefer a measured standard; use this, and its sum rule, when
    there is none.

    Parameters
    ----------
    intensities : ArrayLike
        ``(n,)`` relative intensities, ``I/I0`` if a standard exists or raw
        integrated intensities if not.

    Returns
    -------
    np.ndarray
        The same values divided by their arithmetic mean, so they average to 1,
        read-only.

    Raises
    ------
    ValueError
        If any value is negative or non-finite, or if the mean is not positive.
    """

    values = as_float_array(np.asarray(intensities, dtype=np.float64).reshape(-1), shape=(None,))
    if values.size == 0:
        raise ValueError("harris_texture_coefficients requires at least one intensity.")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Intensities must be finite and non-negative.")
    mean = float(np.mean(values))
    if not mean > 0.0:
        raise ValueError("Intensities must have a positive mean to normalize against.")
    return freeze_array(np.ascontiguousarray(values / mean, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class DiffractogramReflection:
    """One measured peak of a theta-2theta scan, ready for the Kearns route.

    Purpose
    -------
    Bundles what the diffractogram route needs per reflection and nothing else:
    which plane it is, how strong it was, and how strong the same peak is in a
    random powder. Keeping the random intensity beside the measured one makes
    the normalization visible instead of hidden in a scale factor.

    Attributes
    ----------
    plane : CrystalPlane
        The reflecting plane; its tilt to ``[0001]`` is computed, not supplied.
    intensity : float
        Integrated intensity of the peak in the textured specimen, in any
        consistent units.
    random_intensity : float, optional
        Integrated intensity of the same peak from a random powder measured
        under identical conditions. When every reflection carries one, the
        densities are genuinely in times-random units; when none does, the
        Harris normalization is used instead and its assumption is recorded.
    specimen_tilt_deg : float
        Angle between this reflection's diffraction vector and the section
        normal. Zero in a symmetric Bragg-Brentano scan, which is what Kearns'
        derivation assumes. In a **fixed-omega detector scan** it is
        ``theta - omega`` and grows through the pattern, so the reflections no
        longer all probe the same specimen direction; record it and the report
        will say so rather than silently mixing directions.
    """

    plane: CrystalPlane
    intensity: float
    random_intensity: float | None = None
    specimen_tilt_deg: float = 0.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("DiffractogramReflection.intensity must be finite and non-negative.")
        if self.random_intensity is not None and not self.random_intensity > 0.0:
            raise ValueError(
                "DiffractogramReflection.random_intensity must be strictly positive when given."
            )
        if not np.isfinite(self.specimen_tilt_deg):
            raise ValueError("DiffractogramReflection.specimen_tilt_deg must be finite.")

    @property
    def basal_tilt_deg(self) -> float:
        """Angle between this plane's normal and ``[0001]``, from the phase metric."""

        return basal_tilt_angle_deg(self.plane)


def _reflection_densities(
    reflections: Sequence[DiffractogramReflection],
    normalization: str,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    if normalization not in _DIFFRACTOGRAM_NORMALIZATIONS:
        raise ValueError(
            f"normalization must be 'random_standard' or 'harris'; got '{normalization}'."
        )
    if len(reflections) < 2:
        raise ValueError(
            "The diffractogram route needs at least two reflections at different basal tilts."
        )
    tilts = np.array([reflection.basal_tilt_deg for reflection in reflections], dtype=np.float64)
    measured = np.array([reflection.intensity for reflection in reflections], dtype=np.float64)
    notes: list[str] = []
    if normalization == "random_standard":
        missing = [
            _plane_label(reflection.plane)
            for reflection in reflections
            if reflection.random_intensity is None
        ]
        if missing:
            raise ValueError(
                "normalization='random_standard' needs a random_intensity on every "
                f"reflection; missing for {', '.join(missing)}. Measure a random powder "
                "standard, or pass normalization='harris' and accept its assumption."
            )
        standard = np.array(
            [float(reflection.random_intensity or 0.0) for reflection in reflections],
            dtype=np.float64,
        )
        densities = measured / standard
        notes.append(
            "Pole densities are in times-random units from a measured random-powder "
            "standard, which is the normalization Kearns (1965) used."
        )
    else:
        base = measured
        if all(reflection.random_intensity is not None for reflection in reflections):
            base = measured / np.array(
                [float(reflection.random_intensity or 0.0) for reflection in reflections],
                dtype=np.float64,
            )
        densities = np.asarray(harris_texture_coefficients(base), dtype=np.float64)
        notes.append(
            "Harris texture coefficients were used: the pole densities are scaled so that "
            "their mean over the measured reflections is 1. Kearns (1965) measured that "
            "mean at 1.02 to 1.56 (average 1.23) rather than 1, so densities normalized "
            "this way are systematically low by that much."
        )
    return tilts, densities, notes


def basal_tilt_profile(
    reflections: Sequence[DiffractogramReflection],
    *,
    normalization: str = "random_standard",
    bin_width_deg: float = 10.0,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Build the ``I(phi)`` tilt profile from measured diffraction peaks.

    Purpose
    -------
    The step that turns a diffractogram into something Kearns' Eq. (5) can
    integrate. Each reflection's normalized intensity is the basal-pole density
    at that reflection's own tilt to ``[0001]``; those scattered points are
    interpolated onto the uniform tilt grid the quadrature needs.

    Method
    ------
    Reflections sharing a tilt are averaged, then the density is linearly
    interpolated in ``phi`` onto the midpoints of ``bin_width_deg`` bins
    spanning ``[0, 90]``. This is the "linear variation of diffraction
    intensities across different rotation angles" that Mani Krishna *et al.*
    (2011) identify as the route's central approximation, and the reason its
    three values need normalizing before they sum to 1: the available
    reflections are not uniformly distributed over tilt, so the interpolated
    curve is better resolved where they cluster. Outside the measured tilt
    range the nearest measured density is held constant, which is the graphical
    extension Kearns drew by hand.

    Parameters
    ----------
    reflections : Sequence[DiffractogramReflection]
        At least two, at different basal tilts.
    normalization : str
        ``"random_standard"`` (default) divides by a measured random powder;
        ``"harris"`` scales the densities to a mean of one instead. See
        :func:`harris_texture_coefficients`.
    bin_width_deg : float
        Width of the tilt bins. Kearns used 10 degrees.

    Returns
    -------
    tuple
        ``(tilt_midpoints_deg, density, notes)`` — the profile on bin midpoints
        and the caveats the normalization carries.

    Raises
    ------
    ValueError
        If fewer than two reflections are given, if all share one tilt, if the
        bin width does not divide 90, or if a required random intensity is
        missing.
    """

    if bin_width_deg <= 0.0 or not np.isclose(90.0 / bin_width_deg, round(90.0 / bin_width_deg)):
        raise ValueError("bin_width_deg must be positive and divide 90 degrees exactly.")
    tilts, densities, notes = _reflection_densities(reflections, normalization)
    unique_tilts, inverse = np.unique(np.round(tilts, 6), return_inverse=True)
    if unique_tilts.size < 2:
        raise ValueError(
            "All reflections sit at the same basal tilt, so no profile can be interpolated. "
            "Include reflections spanning a range of tilts, e.g. (0002), (10-11), (10-12) "
            "and (11-20) for alpha-Zr."
        )
    averaged = np.zeros(unique_tilts.size, dtype=np.float64)
    np.add.at(averaged, inverse, densities)
    counts = np.bincount(inverse, minlength=unique_tilts.size).astype(np.float64)
    averaged /= counts
    bin_count = round(90.0 / bin_width_deg)
    midpoints = (np.arange(bin_count, dtype=np.float64) + 0.5) * bin_width_deg
    profile = np.interp(midpoints, unique_tilts, averaged)
    extra: list[str] = list(notes)
    if unique_tilts[0] > midpoints[0] or unique_tilts[-1] < midpoints[-1]:
        extra.append(
            f"The measured reflections span {unique_tilts[0]:.1f} to "
            f"{unique_tilts[-1]:.1f} degrees of basal tilt; outside that range the nearest "
            "measured density is held constant, as Kearns extended his curves graphically."
        )
    return (
        freeze_array(np.ascontiguousarray(midpoints)),
        freeze_array(np.ascontiguousarray(profile)),
        tuple(extra),
    )


def kearns_from_diffractogram(
    reflections: Sequence[DiffractogramReflection],
    *,
    specimen_frame: ReferenceFrame,
    basal_pole: CrystalPlane | None = None,
    normalization: str = "random_standard",
    bin_width_deg: float = 10.0,
    direction: ArrayLike | None = None,
    direction_label: str | None = None,
    provenance: ProvenanceRecord | None = None,
) -> KearnsReport:
    """Kearns parameter from a theta-2theta scan — the original 1965 route.

    Purpose
    -------
    Kearns' own method, and still the only one that needs no texture goniometer:
    measure the integrated intensities of every ``(hkil)`` peak on a flat
    section, compare each with a random powder, and read the result as the
    basal-pole density at that reflection's known tilt to ``[0001]``. The
    intensities are exactly an inverse pole figure of the section normal, which
    is why Mani Krishna *et al.* (2011) call this the IPF route.

    When to use
    -----------
    When no goniometer is available, when the material is too deformed for EBSD
    to index reliably, or as an independent check on a pole-figure result.
    Its cost is that ``f`` must be measured on all three principal sections to
    obtain the triad, which for thin product forms such as clad tubing means
    stacking slices, and that the three values do not sum to 1 without
    normalization.

    Method
    ------
    :func:`basal_tilt_profile` builds ``I(phi)`` from the reflections, then
    :func:`kearns_from_tilt_profile` integrates Kearns' Eq. (5) over it. Both
    steps are exposed separately so the intermediate profile — the most
    informative plot in the whole method — can be inspected.

    Parameters
    ----------
    reflections : Sequence[DiffractogramReflection]
        The measured peaks. Include reflections spanning the full tilt range:
        the ``sin phi`` factor makes the high-tilt reflections
        (``(11-20)`` at 90 degrees, ``(10-13)`` at 31.5) matter far more than
        their intensity suggests, and omitting them biases ``f`` upward.
    specimen_frame : ReferenceFrame
    basal_pole : CrystalPlane, optional
        The pole the profile describes, ``(0002)`` of the reflections' phase by
        default.
    normalization : str
        ``"random_standard"`` (default) or ``"harris"``.
    bin_width_deg : float
        Tilt bin width for the profile; Kearns used 10 degrees.
    direction, direction_label
        The specimen direction the section normal represents. Defaults to the
        frame's third axis.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    KearnsReport
        With ``diffraction_vector_tilt_spread_deg`` among the diagnostics.

    Raises
    ------
    ValueError
        If the reflections do not share a phase, if fewer than two distinct
        basal tilts are present, or if a required random intensity is missing.

    Notes
    -----
    **Fixed-omega scans.** Kearns' derivation assumes a symmetric scan, in which
    only planes parallel to the specimen surface diffract. A fixed-omega
    detector scan — common on modern four-circle instruments, where the file
    records ``scanAxis="2Theta"`` with a single ``Omega`` position — violates
    that: the diffraction vector sits ``theta - omega`` from the surface normal
    and moves through the pattern, so different reflections probe different
    specimen directions. Set ``specimen_tilt_deg`` on each reflection and this
    function reports the spread; a spread of more than a few degrees means the
    profile mixes directions and the result is not the section-normal ``f``.

    Examples
    --------
    Applied to Kearns' own Table 2 intensities for the longitudinal section of
    a swaged Zircaloy-2 rod, the route reproduces his tabulated ``f = 0.488``.
    """

    if not reflections:
        raise ValueError("kearns_from_diffractogram requires at least one reflection.")
    phase = reflections[0].plane.phase
    for reflection in reflections[1:]:
        if reflection.plane.phase != phase:
            raise ValueError("Every reflection must belong to the same phase.")
    pole = (
        CrystalPlane.from_miller_bravais((0, 0, 0, 2), phase=phase)
        if basal_pole is None
        else basal_pole
    )
    midpoints, profile, notes = basal_tilt_profile(
        reflections,
        normalization=normalization,
        bin_width_deg=bin_width_deg,
    )
    report = kearns_from_tilt_profile(
        midpoints,
        profile,
        pole=pole,
        specimen_frame=specimen_frame,
        direction=direction,
        direction_label=direction_label,
        provenance=provenance,
    )
    tilt_spread = float(
        np.ptp(np.array([r.specimen_tilt_deg for r in reflections], dtype=np.float64))
    )
    extra = list(notes)
    if tilt_spread > 2.0:
        extra.append(
            f"The reflections' diffraction vectors span {tilt_spread:.1f} degrees of "
            "specimen tilt, so they do not all probe the section normal. This is the "
            "signature of a fixed-omega detector scan; Kearns' derivation assumes a "
            "symmetric scan, in which the spread is zero."
        )
    basal_tilts = sorted({round(r.basal_tilt_deg, 1) for r in reflections})
    return KearnsReport(
        values=report.values,
        directions=report.directions,
        direction_labels=report.direction_labels,
        method="diffractogram",
        pole=pole,
        specimen_frame=report.specimen_frame,
        orientation_tensor=None,
        diagnostics={
            **dict(report.diagnostics),
            "reflection_count": float(len(reflections)),
            "distinct_basal_tilt_count": float(len(basal_tilts)),
            "min_basal_tilt_deg": float(basal_tilts[0]),
            "max_basal_tilt_deg": float(basal_tilts[-1]),
            "diffraction_vector_tilt_spread_deg": tilt_spread,
        },
        notes=(
            *extra,
            "One section determines f along one direction; the three principal sections "
            "are needed for a triad, and Kearns (1965) noted that their sum departs from 1 "
            "because the available reflections are not uniformly spread over tilt.",
        ),
        provenance=provenance,
    )
