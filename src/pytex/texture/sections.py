"""ODF sections in any Euler coordinate, over the ranges symmetry actually needs.

What it does
    Slices an orientation distribution - discrete (:class:`~pytex.texture.ODF`)
    or harmonic (:class:`~pytex.texture.HarmonicODF`) - into two-dimensional
    sections at constant ``phi2`` (the default and the convention of most texture
    papers), constant ``phi1``, or constant ``sigma = phi1 + phi2``, and returns
    the densities in multiples of a random distribution for both representations.

Why the ranges are computed rather than fixed
    The box of Euler space a section has to cover is set by the symmetries, not
    by habit. A cubic crystal repeats every 90 degrees in ``phi2``, a hexagonal
    one every 60; an orthorhombic specimen folds ``phi1`` into ``[0, 90]`` where
    a specimen with no symmetry needs ``[0, 360]``. A section drawn over the
    cubic box for a hexagonal phase repeats a third of itself, and one drawn over
    ``[0, 90]`` of ``phi1`` for a triclinic specimen hides three quarters of the
    texture. :func:`euler_section_ranges` derives each limit from the operators,
    so the plotted box is the asymmetric one the conventions of Bunge (1982,
    section 4.2) and LaboTex describe.

References
    Bunge, *Texture Analysis in Materials Science* (1982), sections 2.3 and 4.2;
    Randle and Engler, *Introduction to Texture Analysis*, 2nd ed., section 5.2.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core.orientation import OrientationSet
from pytex.core.symmetry import SymmetrySpec
from pytex.texture.models import ODFSectionData

__all__ = [
    "ODF_SECTION_KINDS",
    "euler_section_ranges",
    "odf_sections",
]

#: The section coordinates offered, in the order texture software offers them.
ODF_SECTION_KINDS: tuple[str, ...] = ("phi2", "phi1", "sigma")

#: Query orientations evaluated against the support per block. Bounds the
#: misorientation work array of a discrete ODF, which scales with the product
#: of query count, support size and crystal symmetry order.
_BLOCK_PAIRS = 400_000


def _rotation_order_about_z(operators: np.ndarray) -> int:
    """How many operators fix the third axis: the order of the principal rotation."""

    return int(np.sum(np.isclose(operators[:, 2, 2], 1.0, atol=1e-8)))


def _has_two_fold_perpendicular_to_z(operators: np.ndarray) -> bool:
    return bool(np.any(np.isclose(operators[:, 2, 2], -1.0, atol=1e-8)))


def euler_section_ranges(
    crystal_symmetry: SymmetrySpec | None,
    specimen_symmetry: SymmetrySpec | str | None = None,
) -> dict[str, float]:
    """The Bunge-Euler box an ODF section must cover under given symmetries.

    Purpose
    -------
    Decide how much of each Euler angle a set of sections has to show so that
    every distinct orientation appears once - no more, which would repeat part
    of the texture, and no less, which would hide part of it.

    Method
    ------
    * ``phi2`` spans ``360 / n``, where ``n`` is the order of the crystal's
      rotation axis along its third axis: 90 degrees for cubic and tetragonal,
      60 for hexagonal, 120 for trigonal, 180 for a two-fold, 360 for none.
    * ``Phi`` spans 90 degrees when a two-fold axis lies perpendicular to the
      third axis of the crystal or of the specimen, and 180 otherwise.
    * ``phi1`` spans 360 degrees for a triclinic specimen, 180 for monoclinic
      and 90 for orthorhombic. For an **axial** specimen the density does not
      depend on ``phi1`` at all; 90 degrees is returned so a section is still
      drawn, and every column of it is the same.

    Parameters
    ----------
    crystal_symmetry : SymmetrySpec, optional
        The crystal's proper rotations. ``None`` is treated as triclinic.
    specimen_symmetry : SymmetrySpec or str, optional
        A specimen symmetry or its name. ``None`` is treated as triclinic.

    Returns
    -------
    dict
        ``phi1_max_deg``, ``big_phi_max_deg`` and ``phi2_max_deg``.
    """

    crystal_operators = (
        np.eye(3)[None, :, :]
        if crystal_symmetry is None
        else np.asarray(crystal_symmetry.operators, dtype=np.float64)
    )
    if isinstance(specimen_symmetry, SymmetrySpec):
        specimen_name = (specimen_symmetry.specimen_symmetry or "triclinic").lower()
    else:
        specimen_name = (specimen_symmetry or "triclinic").strip().lower()
    if specimen_name in {"fibre", "fiber", "cylindrical"}:
        specimen_name = "axial"
    phi1_max = {"triclinic": 360.0, "monoclinic": 180.0, "orthorhombic": 90.0,
                "orthotropic": 90.0, "axial": 90.0}.get(specimen_name, 360.0)
    order = max(_rotation_order_about_z(crystal_operators), 1)
    specimen_two_fold = specimen_name in {"orthorhombic", "orthotropic", "axial"}
    big_phi_max = (
        90.0
        if _has_two_fold_perpendicular_to_z(crystal_operators) or specimen_two_fold
        else 180.0
    )
    return {
        "phi1_max_deg": phi1_max,
        "big_phi_max_deg": big_phi_max,
        "phi2_max_deg": 360.0 / order,
    }


def _odf_frames(odf: Any) -> tuple[Any, Any, Any, Any]:
    if hasattr(odf, "quadrature_orientations"):
        return odf.crystal_frame, odf.specimen_frame, odf.crystal_symmetry, odf.phase
    support = odf.orientations
    return support.crystal_frame, support.specimen_frame, support.symmetry, support.phase


def _evaluate_mrd(odf: Any, euler_deg: np.ndarray) -> np.ndarray:
    """Density in multiples of random at Bunge angles, for either ODF representation.

    A harmonic ODF evaluates in m.r.d. already. A discrete ODF's
    ``evaluate(normalized=True)`` folds each query into the fundamental zone,
    which is ``1/|G|`` of SO(3), so a random texture reads ``|G|``; dividing by
    the crystal symmetry order puts it on the same scale.
    """

    crystal_frame, specimen_frame, symmetry, phase = _odf_frames(odf)
    harmonic = hasattr(odf, "quadrature_orientations")
    support_size = 1 if harmonic else len(odf.orientations)
    order = 1 if symmetry is None else int(np.asarray(symmetry.operators).shape[0])
    block = max(1, _BLOCK_PAIRS // max(support_size * (1 if harmonic else order), 1))
    values = np.empty(euler_deg.shape[0], dtype=np.float64)
    for start in range(0, euler_deg.shape[0], block):
        stop = min(start + block, euler_deg.shape[0])
        query = OrientationSet.from_euler_angles(
            euler_deg[start:stop],
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            symmetry=symmetry,
            phase=phase,
            convention="bunge",
            degrees=True,
        )
        if harmonic:
            values[start:stop] = np.asarray(odf.evaluate(query), dtype=np.float64)
        else:
            values[start:stop] = (
                np.asarray(odf.evaluate(query, normalized=True), dtype=np.float64) / order
            )
    return values


def odf_sections(
    odf: Any,
    *,
    kind: str = "phi2",
    values_deg: ArrayLike | None = None,
    specimen_symmetry: SymmetrySpec | str | None = None,
    step_deg: float = 5.0,
    resolution_deg: float = 5.0,
    ranges: dict[str, float] | None = None,
) -> ODFSectionData:
    """Sample an ODF on sections of constant ``phi2``, ``phi1`` or ``sigma``.

    Purpose
    -------
    The standard way an orientation distribution is read: a stack of
    two-dimensional maps over Euler space, drawn over the box the symmetries
    require, in multiples of a random distribution whichever representation the
    ODF uses.

    When to use
    -----------
    ``phi2`` sections are the default of the texture literature and of LaboTex;
    for cubic metals the ``phi2 = 0, 45, 65`` trio carries the standard rolling
    components. ``phi1`` sections suit textures organised along ``phi2``, and
    ``sigma`` sections the gamma fibre of rolled bcc steels. Leave ``values_deg``
    unset for the LaboTex-style plate: every section at ``step_deg`` over the
    whole range.

    Parameters
    ----------
    odf : ODF or HarmonicODF
    kind : str
        ``"phi2"`` (default), ``"phi1"`` or ``"sigma"``.
    values_deg : ArrayLike, optional
        The constant value of each section. Default: ``0, step, 2 step, ...``
        through the range of the sectioning coordinate.
    specimen_symmetry : SymmetrySpec or str, optional
        Used only to size the box when ``ranges`` is not given; the ODF's own
        declared specimen symmetry is used when this is omitted.
    step_deg : float
        Spacing of the default section values.
    resolution_deg : float
        Grid spacing inside each section.
    ranges : dict, optional
        Explicit ``phi1_max_deg``, ``big_phi_max_deg`` and ``phi2_max_deg``,
        overriding :func:`euler_section_ranges`.

    Returns
    -------
    ODFSectionData
        ``densities[i, Phi, across]`` in m.r.d.; the across coordinate is named by
        :attr:`~pytex.texture.ODFSectionData.horizontal_coordinate`.

    Raises
    ------
    ValueError
        For an unknown kind or non-positive spacings.
    """

    if kind not in ODF_SECTION_KINDS:
        raise ValueError(f"kind must be one of {', '.join(ODF_SECTION_KINDS)}; got {kind!r}.")
    if step_deg <= 0.0 or resolution_deg <= 0.0:
        raise ValueError("step_deg and resolution_deg must be strictly positive.")
    _crystal_frame, _specimen_frame, crystal_symmetry, _phase = _odf_frames(odf)
    if ranges is None:
        declared = specimen_symmetry
        if declared is None:
            declared = getattr(odf, "specimen_symmetry", None)
        ranges = euler_section_ranges(crystal_symmetry, declared)
    phi1_max = float(ranges["phi1_max_deg"])
    big_phi_max = float(ranges["big_phi_max_deg"])
    phi2_max = float(ranges["phi2_max_deg"])

    section_max = {"phi2": phi2_max, "phi1": phi1_max, "sigma": phi1_max + phi2_max}[kind]
    if values_deg is None:
        # Both edges of the box are drawn, as LaboTex and the printed atlases do:
        # the last section repeats the first under symmetry, and showing it lets
        # a component on the boundary be read without wrapping round.
        count = int(np.floor(section_max / step_deg + 1e-9))
        section_values = np.arange(count + 1) * step_deg
        section_values = section_values[section_values <= section_max + 1e-9]
    else:
        section_values = np.atleast_1d(np.asarray(values_deg, dtype=np.float64))

    across_max = phi2_max if kind == "phi1" else phi1_max
    across = np.arange(0.0, across_max + 1e-9, resolution_deg)
    big_phi = np.arange(0.0, big_phi_max + 1e-9, resolution_deg)
    grid_across, grid_big_phi = np.meshgrid(across, big_phi, indexing="xy")
    flat_across = grid_across.ravel()
    flat_big_phi = grid_big_phi.ravel()

    blocks = []
    for value in section_values:
        constant = np.full(flat_across.size, float(value))
        if kind == "phi2":
            euler = np.column_stack([flat_across, flat_big_phi, constant])
        elif kind == "phi1":
            euler = np.column_stack([constant, flat_big_phi, flat_across])
        else:
            euler = np.column_stack(
                [flat_across, flat_big_phi, np.mod(constant - flat_across, 360.0)]
            )
        blocks.append(euler)
    densities = _evaluate_mrd(odf, np.concatenate(blocks, axis=0)).reshape(
        section_values.size, big_phi.size, across.size
    )
    return ODFSectionData(
        phi2_deg=section_values,
        phi1_deg=across,
        big_phi_deg=big_phi,
        densities=densities,
        normalized=True,
        section_kind=kind,
    )
