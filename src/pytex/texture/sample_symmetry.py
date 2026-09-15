"""Statistical sample symmetry imposed on measured pole figures.

What it does
    Averages a measured pole figure over the operations of an assumed sample
    (specimen) symmetry, on the figure's own measured directions. The result is
    the figure a specimen with exactly that symmetry would give, and it is the
    form in which a sample-symmetry assumption enters a pole-figure inversion.

When to use it
    When the process that made the specimen imposes a symmetry and the
    measurement's statistics are worth improving by using it: rolled sheet is
    conventionally orthorhombic about RD, TD and ND; drawn wire, extruded rod and
    a tube read along its axis are **axial** (fibre) about their axis. Imposing a
    symmetry the specimen does not have fabricates it, so the choice is recorded
    on the returned figure rather than implied.

The axial case is treated exactly
    Axial symmetry is the continuous group of every rotation about the fibre axis
    together with the two-fold axes perpendicular to it (the Curie group
    ``infinity/mm``, whose proper part is ``infinity 2 2``). No finite operator
    list represents it, so it is not approximated by one here: an axially
    symmetric figure depends on the polar angle alone, and the average over the
    continuous group is the azimuthal mean on each ring of constant polar angle,
    taken with circular trapezoid weights so an unevenly spaced ring is still
    integrated correctly. The two-fold axes add nothing further for an antipodal
    figure, because they map a pole at polar angle ``theta`` to ``180 - theta``,
    which the antipodal fold returns to ``theta``.

Finite groups are averaged on the measured support
    For triclinic (nothing), monoclinic and orthorhombic symmetry the orbit of
    every measured direction is formed and each image is read from the nearest
    measured direction. Every operation of those groups maps a polar angle onto
    the same polar angle once the antipodal fold is applied, so the images stay
    inside a measured cap truncated at any tilt: no value is extrapolated from
    outside the measured range.

References
    Bunge, *Texture Analysis in Materials Science* (1982), section 4.2 (sample
    symmetry); Randle and Engler, *Introduction to Texture Analysis*, 2nd ed.,
    section 5.3; Kocks, Tome and Wenk, *Texture and Anisotropy* (1998), chapter 1.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from pytex.core.symmetry import SymmetrySpec
from pytex.texture.models import PoleFigure

__all__ = [
    "SAMPLE_SYMMETRY_NAMES",
    "impose_sample_symmetry",
    "sample_symmetry_spec",
]

#: The sample symmetries the texture workflows offer, in increasing order of
#: what they assume. ``"axial"`` is the fibre symmetry about the specimen's third
#: axis (ND for sheet, the axis of a rod).
SAMPLE_SYMMETRY_NAMES: tuple[str, ...] = ("triclinic", "monoclinic", "orthorhombic", "axial")


def sample_symmetry_spec(
    symmetry: str | SymmetrySpec, *, pole_figure: PoleFigure | None = None
) -> SymmetrySpec:
    """Resolve a sample-symmetry name to a specimen-frame :class:`SymmetrySpec`.

    Parameters
    ----------
    symmetry : str or SymmetrySpec
        One of :data:`SAMPLE_SYMMETRY_NAMES`, or an already-built specification.
    pole_figure : PoleFigure, optional
        Supplies the specimen frame the specification is declared in.

    Returns
    -------
    SymmetrySpec
    """

    if isinstance(symmetry, SymmetrySpec):
        return symmetry
    frame = None if pole_figure is None else pole_figure.specimen_frame
    return SymmetrySpec.specimen(str(symmetry), reference_frame=frame)


def _ring_labels(polar_deg: np.ndarray, tolerance_deg: float) -> np.ndarray:
    """Group polar angles into rings: consecutive sorted angles closer than the tolerance."""

    order = np.argsort(polar_deg, kind="stable")
    sorted_polar = polar_deg[order]
    breaks = np.concatenate([[0], (np.diff(sorted_polar) > tolerance_deg).astype(np.int64)])
    labels_sorted = np.cumsum(breaks)
    labels = np.empty_like(labels_sorted)
    labels[order] = labels_sorted
    return labels


def _axial_average(
    directions: np.ndarray, values: np.ndarray, *, ring_tolerance_deg: float | None
) -> tuple[np.ndarray, int]:
    polar = np.degrees(np.arccos(np.clip(directions[:, 2], -1.0, 1.0)))
    if ring_tolerance_deg is None:
        unique = np.unique(np.round(polar, 6))
        gaps = np.diff(unique)
        gaps = gaps[gaps > 1e-6]
        # Half the smallest spacing between distinct polar angles: a goniometer
        # raster's tilt rings are separated by its tilt step, and a scattered
        # support has no rings to preserve, so its points fall into bands.
        ring_tolerance_deg = 0.5 * float(gaps.min()) if gaps.size else 1e-3
    labels = _ring_labels(polar, float(ring_tolerance_deg))
    azimuth = np.arctan2(directions[:, 1], directions[:, 0]) % (2.0 * np.pi)
    averaged = np.empty_like(values)
    ring_count = int(labels.max()) + 1
    for label in range(ring_count):
        members = np.flatnonzero(labels == label)
        if members.size == 1 or np.max(np.sin(np.radians(polar[members]))) < 1e-9:
            # A single point, or the pole itself, is its own ring average.
            averaged[members] = float(np.mean(values[members]))
            continue
        order = members[np.argsort(azimuth[members], kind="stable")]
        phi = azimuth[order]
        # Circular trapezoid weights: each point owns half of the arc to each
        # neighbour, so an unevenly sampled ring is integrated, not merely averaged.
        forward = np.diff(np.concatenate([phi, [phi[0] + 2.0 * np.pi]]))
        backward = np.roll(forward, 1)
        weights = 0.5 * (forward + backward)
        total = float(weights.sum())
        mean = (
            float(np.dot(weights, values[order]) / total)
            if total > 0.0
            else float(np.mean(values[order]))
        )
        averaged[members] = mean
    return averaged, ring_count


def impose_sample_symmetry(
    pole_figure: PoleFigure,
    symmetry: str | SymmetrySpec,
    *,
    ring_tolerance_deg: float | None = None,
) -> PoleFigure:
    """Average a measured pole figure over a sample symmetry, on its own directions.

    Purpose
    -------
    Make a sample-symmetry assumption explicit and apply it once, before a
    figure is inverted or compared: the returned figure is what a specimen with
    that symmetry would give, on exactly the directions that were measured.

    When to use
    -----------
    Before an ODF inversion of figures from a specimen whose process fixes a
    symmetry - orthorhombic for rolled sheet, axial for drawn or extruded
    product. Do not use it to make a noisy figure look tidy: a symmetry the
    specimen lacks is fabricated, not revealed.

    Parameters
    ----------
    pole_figure : PoleFigure
        A figure of sampled densities, such as a diffractometer raster. It keeps
        its directions; only the intensities change.
    symmetry : str or SymmetrySpec
        ``"triclinic"`` (no assumption), ``"monoclinic"``, ``"orthorhombic"`` or
        ``"axial"``, or a specimen-frame :class:`SymmetrySpec`.
    ring_tolerance_deg : float, optional
        Axial symmetry only: polar angles closer than this belong to one ring.
        Inferred as half the smallest spacing between distinct polar angles.

    Returns
    -------
    PoleFigure
        The symmetrized figure, with ``sample_symmetry`` recorded.

    Notes
    -----
    Axial averaging is exact (the continuous group); finite groups read every
    orbit image from the nearest measured direction. The area-weighted mean of
    the figure is preserved by the axial average exactly and by the finite-group
    average to the accuracy of the nearest-neighbour lookup, because every
    operation maps the measured support onto itself when the raster is closed
    under it.
    """

    spec = sample_symmetry_spec(symmetry, pole_figure=pole_figure)
    directions = np.asarray(pole_figure.sample_directions, dtype=np.float64)
    values = np.asarray(pole_figure.intensities, dtype=np.float64)
    if pole_figure.antipodal:
        directions = np.where(directions[:, 2:3] < 0.0, -directions, directions)

    name = (spec.specimen_symmetry or "").lower()
    if name == "axial":
        averaged, _rings = _axial_average(
            directions, values, ring_tolerance_deg=ring_tolerance_deg
        )
    else:
        operators = np.asarray(spec.operators, dtype=np.float64)
        images = np.einsum("oij,nj->oni", operators, directions, optimize=True)
        if pole_figure.antipodal:
            images = np.where(images[..., 2:3] < 0.0, -images, images)
        tree = cKDTree(directions)
        _distances, indices = tree.query(images.reshape(-1, 3))
        averaged = values[indices.reshape(operators.shape[0], -1)].mean(axis=0)

    return PoleFigure(
        pole=pole_figure.pole,
        sample_directions=pole_figure.sample_directions,
        intensities=np.clip(averaged, 0.0, None),
        specimen_frame=pole_figure.specimen_frame,
        antipodal=pole_figure.antipodal,
        sample_symmetry=SymmetrySpec(
            name=spec.name,
            point_group=spec.point_group,
            operators=spec.operators,
            specimen_symmetry=spec.specimen_symmetry,
            reference_frame=pole_figure.specimen_frame,
            provenance=spec.provenance,
        ),
        provenance=pole_figure.provenance,
        includes_symmetry_family=pole_figure.includes_symmetry_family,
        sampling=pole_figure.sampling,
    )
