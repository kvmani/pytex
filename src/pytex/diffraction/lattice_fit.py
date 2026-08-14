"""Fitting a two-dimensional lattice to picked diffraction spots.

Why this exists
---------------
A zone-axis diffraction pattern is a *plane* of the reciprocal lattice, so every
spot on it sits at ``c + m a + n b`` for integers ``m``, ``n`` and two basis
vectors the crystal fixes. That constraint is available before any
crystallography: it needs no phase, no camera constant and no zone axis. Only the
picks.

Two things fall out of imposing it, and both matter more than they sound.

**The beam centre stops being a guess.** Picking the transmitted beam by eye is
the largest avoidable error in the whole workflow — it biases every d-spacing at
once, and in a way that still produces a self-consistent, wrong answer. But if
the spots lie on a lattice the centre is *over-determined*: four or more spots
give more equations than unknowns, and least squares returns the centre that best
explains all of them at once. The user's pick becomes a starting point rather
than a measurement.

**A mis-picked spot becomes visible.** Draw the fitted lattice over the pattern
and a spot clicked one node off, or clicked on a dust particle, stops matching
the grid. That is a judgement a person makes instantly from a picture and
struggles to make from a residual column.

Method
------
Seed, then iterate assign-and-refine, then trim:

1. **Seed** from the shortest *differences between picked spots* whose
   directions differ. Differences, not offsets from the picked centre, and
   several candidates rather than one — see :func:`_candidate_bases` for the two
   failures those choices avoid.
2. **Assign** every spot to a node by rounding ``inv(B) (p - c)`` to integers.
3. **Refine** ``a``, ``b`` and ``c`` together by linear least squares with the
   integers held fixed, then clamp ``c`` to stay within half a lattice spacing of
   where it was picked — see :func:`fit_planar_lattice` on why a refinement may
   adjust the centre but must not relabel which node it is.
4. **Trim.** Refit without the worst outlier, up to twice, and score every spot
   against the trimmed fit. Ordinary least squares spreads one bad pick across
   every residual and then points at the wrong spot; trimming is what makes the
   answer name the pick that is actually wrong.

The refinement is linear because the model is linear in exactly the parameters
being fitted, the integers being known by then. There is no optimizer, no
starting-value sensitivity, and no convergence criterion beyond "the integers
stopped moving".

Limits
------
Purely geometric. This does not know what phase the pattern is, and a good fit is
*not* an indexing: it says the picks are consistent with some 2D lattice, which
is necessary for a correct indexing and nowhere near sufficient. Use it to place
the centre and to check the picks; use
:func:`pytex.diffraction.solving.solve_saed_pattern` to find out what the pattern
is.

The lattice fitted is the one the *picked* spots generate, which can be finer
than the lattice of the obvious strong reflections — a hexagonal zone routinely
has weaker reflections at half-integer positions relative to its evident
rectangle. That is correct behaviour, and it is one reason the overlay is worth
looking at rather than merely trusting.

Double diffraction and a picked spot that is not a reflection both show up the
same way: as a point the lattice does not explain. The fit reports them rather
than absorbing them.

See :doc:`/workflows/tem_pattern_indexing` for the workflow, and
:doc:`/theory/saed_ratio_angle_indexing` for the indexing that follows.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array

__all__ = [
    "DEFAULT_CENTRE_LEASH_FRACTION",
    "DEFAULT_INLIER_FRACTION",
    "LatticeFitSpot",
    "PlanarLatticeFit",
    "fit_planar_lattice",
]

#: Residual tolerance, as a fraction of the shortest separation between two
#: picked spots, for a spot to count as explained by the lattice.
#:
#: This is a stand-in for *picking precision*, not for lattice spacing: it is how
#: close a click has to land to a node to be the same thing. Four percent of the
#: shortest spot separation is a few pixels on a typical plate, which is about
#: how well a spot centre can be clicked.
#:
#: It is deliberately **not** a fraction of the fitted basis. Tying the tolerance
#: to each candidate lattice would make the node density seen at that tolerance
#: identical for a cell and for the same cell halved, which is exactly the
#: distinction the ranking exists to make. See :meth:`_Attempt.evidence`.
DEFAULT_INLIER_FRACTION = 0.04

#: How far the refinement may move the centre, as a fraction of the shortest
#: basis vector. See :func:`fit_planar_lattice`.
DEFAULT_CENTRE_LEASH_FRACTION = 0.5

#: Cosine above which two picked directions are treated as the same direction.
#:
#: Friedel pairs make this necessary rather than merely tidy: ``g`` and ``-g``
#: are both picked routinely and are collinear through the beam, so a seed basis
#: chosen by length alone would be degenerate.
_COLLINEAR_COSINE = 0.999

#: Most spots the trimming stage will set aside.
_MAX_TRIMMED = 2

#: Below this, a trial basis is treated as degenerate and its fit abandoned.
#:
#: Compared against the product of the two vector lengths, so it is a sine of the
#: angle between them rather than an absolute area — a scale-free test that means
#: the same thing for a pattern picked in pixels and one picked in millimetres.
_DEGENERATE_SINE = 1e-6

#: How many trial bases to fit before choosing between them.
#:
#: Each costs a handful of 2x2 solves on a few dozen points, so the search is
#: free in any practical sense; six is enough to survive two bad picks creating
#: spuriously short differences between them.
_SEED_CANDIDATES = 6


@dataclass(frozen=True, slots=True)
class LatticeFitSpot:
    """One picked spot, as the fitted lattice explains it.

    Attributes
    ----------
    index : int
        Position in the picked list, so a row traces back to a click.
    position : np.ndarray
        Where it was picked, in the picking coordinates.
    lattice_indices : tuple of int
        The ``(m, n)`` node it was assigned to, relative to the fitted centre.
    predicted : np.ndarray
        Where that node lies under the fitted basis.
    residual : float
        Distance from the pick to the node, in the picking coordinates' units.
    inlier : bool
        Whether the residual is inside the tolerance. An outlier is reported,
        never silently dropped: a spot the lattice cannot explain is a finding.
    """

    index: int
    position: np.ndarray
    lattice_indices: tuple[int, int]
    predicted: np.ndarray
    residual: float
    inlier: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", as_float_array(self.position, shape=(2,)))
        object.__setattr__(self, "predicted", as_float_array(self.predicted, shape=(2,)))
        if self.index < 0:
            raise ValueError("LatticeFitSpot.index must be non-negative.")

    def to_json(self) -> dict[str, Any]:
        """The row as JSON-ready data."""

        return {
            "index": int(self.index),
            "x": float(self.position[0]),
            "y": float(self.position[1]),
            "m": int(self.lattice_indices[0]),
            "n": int(self.lattice_indices[1]),
            "predicted_x": float(self.predicted[0]),
            "predicted_y": float(self.predicted[1]),
            "residual": float(self.residual),
            "inlier": bool(self.inlier),
        }


@dataclass(frozen=True, slots=True)
class PlanarLatticeFit:
    """A two-dimensional lattice fitted to picked spots, with a refined centre.

    Attributes
    ----------
    centre : np.ndarray
        The transmitted-beam position the fit prefers.
    supplied_centre : np.ndarray
        What the user picked, kept so the correction can be reported.
    basis : np.ndarray
        ``(2, 2)`` with the two basis vectors as **rows**, in the picking
        coordinates.
    spots : tuple of LatticeFitSpot
    rms_residual : float
        Root-mean-square distance from an inlier pick to its node.
    iterations : int
        Assign-and-refine rounds before the integers stopped moving.
    notes : tuple of str
        Findings a user must see: two picks on one node, a centre that reached
        its leash, and anything else the fit detects but cannot resolve.
    """

    centre: np.ndarray
    supplied_centre: np.ndarray
    basis: np.ndarray
    spots: tuple[LatticeFitSpot, ...]
    rms_residual: float
    iterations: int
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "notes", tuple(self.notes))
        object.__setattr__(self, "centre", as_float_array(self.centre, shape=(2,)))
        object.__setattr__(
            self, "supplied_centre", as_float_array(self.supplied_centre, shape=(2,))
        )
        object.__setattr__(self, "basis", as_float_array(self.basis, shape=(2, 2)))
        object.__setattr__(self, "spots", tuple(self.spots))

    @property
    def inlier_count(self) -> int:
        """How many picks the lattice explains."""

        return sum(1 for spot in self.spots if spot.inlier)

    @property
    def outliers(self) -> tuple[LatticeFitSpot, ...]:
        """The picks the lattice does not explain — the ones worth looking at."""

        return tuple(spot for spot in self.spots if not spot.inlier)

    @property
    def centre_shift(self) -> float:
        """How far the refinement moved the centre from where it was picked."""

        return float(np.linalg.norm(self.centre - self.supplied_centre))

    @property
    def basis_lengths(self) -> tuple[float, float]:
        """Lengths of the two basis vectors, shortest first."""

        lengths = sorted(float(np.linalg.norm(row)) for row in self.basis)
        return (lengths[0], lengths[1])

    @property
    def basis_angle_deg(self) -> float:
        """Angle between the two basis vectors, in degrees.

        Reported in ``[0, 180)``: the two vectors span a cell, and which of the
        supplementary pair is quoted is a convention, not a measurement.
        """

        first, second = self.basis[0], self.basis[1]
        cosine = float(np.dot(first, second)) / (
            float(np.linalg.norm(first)) * float(np.linalg.norm(second))
        )
        return float(math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0)))))

    def node_positions(
        self, *, max_index: int = 6, bounds: tuple[float, float] | None = None
    ) -> np.ndarray:
        """Positions of the lattice nodes, for drawing the overlay.

        Parameters
        ----------
        max_index : int
            Largest ``|m|`` and ``|n|`` generated.
        bounds : tuple of float, optional
            ``(width, height)`` of the image. Nodes outside it are dropped, so an
            overlay carries no points nobody can see.

        Returns
        -------
        np.ndarray
            ``(k, 4)`` of ``x, y, m, n``.
        """

        if max_index <= 0:
            raise ValueError("max_index must be strictly positive.")
        span = np.arange(-max_index, max_index + 1)
        grid = np.stack(np.meshgrid(span, span, indexing="ij"), axis=-1).reshape(-1, 2)
        positions = self.centre + grid.astype(float) @ self.basis
        stacked = np.column_stack([positions, grid.astype(float)])
        if bounds is not None:
            width, height = float(bounds[0]), float(bounds[1])
            inside = (
                (stacked[:, 0] >= 0.0)
                & (stacked[:, 0] <= width)
                & (stacked[:, 1] >= 0.0)
                & (stacked[:, 1] <= height)
            )
            stacked = stacked[inside]
        return np.ascontiguousarray(stacked)

    def describe(self) -> str:
        """A prose account of the fit, stating what it does and does not prove."""

        short, long_ = self.basis_lengths
        head = (
            f"A two-dimensional lattice fitted to {len(self.spots)} picked spot(s): basis vectors "
            f"of {short:.2f} and {long_:.2f} picking units at {self.basis_angle_deg:.2f} degrees "
            f"to each other, explaining {self.inlier_count} of them with an r.m.s. residual of "
            f"{self.rms_residual:.3f} units after {self.iterations} refinement round(s). "
        )
        centre = (
            f"The transmitted beam refines to ({self.centre[0]:.2f}, {self.centre[1]:.2f}), "
            f"{self.centre_shift:.2f} units from where it was picked. "
            if self.centre_shift > 1e-9
            else f"The transmitted beam was held at ({self.centre[0]:.2f}, {self.centre[1]:.2f}). "
        )
        outliers = (
            ""
            if not self.outliers
            else (
                f"{len(self.outliers)} pick(s) do not lie on the lattice: spot(s) "
                + ", ".join(str(spot.index + 1) for spot in self.outliers)
                + ". A double-diffraction reflection, a second phase, or a mis-click are the "
                "usual causes, in that order. "
            )
        )
        remarks = (" ".join(self.notes) + " ") if self.notes else ""
        return (
            head
            + centre
            + outliers
            + remarks
            + (
                "This is geometry, not indexing: a good fit says the picks are consistent with "
                "some two-dimensional lattice, which is necessary for a correct indexing and far "
                "from sufficient."
            )
        )

    def to_json(self) -> dict[str, Any]:
        """The fit as JSON-ready data."""

        return {
            "centre": [float(value) for value in self.centre],
            "supplied_centre": [float(value) for value in self.supplied_centre],
            "centre_shift": float(self.centre_shift),
            "basis": [[float(value) for value in row] for row in self.basis],
            "basis_lengths": list(self.basis_lengths),
            "basis_angle_deg": float(self.basis_angle_deg),
            "rms_residual": float(self.rms_residual),
            "inlier_count": int(self.inlier_count),
            "iterations": int(self.iterations),
            "notes": list(self.notes),
            "spots": [spot.to_json() for spot in self.spots],
            "describe": self.describe(),
        }


def _candidate_bases(picked: np.ndarray, centre: np.ndarray, limit: int) -> list[np.ndarray]:
    """Trial bases to try, shortest first, as (2, 2) arrays of row vectors.

    **Differences between spots, not offsets from the picked centre.** A spot
    sits at ``c + v`` with ``v`` a lattice vector, so a *difference* of two spots
    is a lattice vector however badly ``c`` was picked, while an *offset from the
    picked centre* is one only if the pick was already right. Seeding from
    offsets fails in exactly the case this module exists to repair: a centre half
    a spacing out makes some offsets spuriously short, the shortest pair
    generates a sub-lattice, and the fit then explains every spot perfectly by
    halving the cell around the wrong centre — machine-precision residuals, an
    uncorrected centre, and nothing anywhere to say something is wrong.

    **Several candidates, not one.** The shortest difference is a genuine lattice
    vector only if every pick is genuine. One spot clicked forty pixels off
    creates short differences that are not lattice vectors at all, and a seed
    taken from one of those produces a confident fit to a lattice that does not
    exist. So the shortest few independent pairs are all tried, and the caller
    keeps whichever *fit* comes out best. A spurious seed loses on inlier count,
    which is the observable that a wrong lattice cannot fake.
    """

    if picked.shape[0] >= 3:
        differences = (picked[:, None, :] - picked[None, :, :]).reshape(-1, 2)
    else:
        differences = picked - centre
    lengths = np.linalg.norm(differences, axis=1)
    order = [index for index in np.argsort(lengths) if lengths[index] > 1e-12]

    # Deduplicate by direction and length, so `a` and `-a` do not both start a
    # candidate and produce the same fit twice.
    unique: list[np.ndarray] = []
    for index in order:
        vector = differences[index]
        if any(
            abs(abs(float(np.dot(vector, kept))) / (lengths[index] * np.linalg.norm(kept)) - 1.0)
            < 1e-9
            and abs(lengths[index] - float(np.linalg.norm(kept))) < 1e-9
            for kept in unique
        ):
            continue
        unique.append(vector)
        if len(unique) >= 2 * limit + 2:
            break

    bases: list[np.ndarray] = []
    for first_index, first in enumerate(unique):
        for second in unique[first_index + 1 :]:
            cosine = abs(
                float(np.dot(first, second))
                / (float(np.linalg.norm(first)) * float(np.linalg.norm(second)))
            )
            if cosine >= _COLLINEAR_COSINE:
                continue
            bases.append(np.vstack([first, second]))
            break
        if len(bases) >= limit:
            break
    if not bases:
        raise ValueError(
            "The picked spots lie on a single line through the transmitted beam, so they span a "
            "row rather than a lattice. Pick a spot off that row."
        )
    return bases


def _gauss_reduce(basis: np.ndarray) -> np.ndarray:
    """Return the shortest pair of vectors generating the same 2D lattice.

    A lattice has infinitely many bases — ``(a, b)`` and ``(a, b - a)`` generate
    the same nodes — so any quantity read off an arbitrary fitted basis is a
    property of the arithmetic rather than of the crystal. Fit a square lattice
    and the raw answer may well come back as two vectors 135 degrees apart, which
    is correct and useless: the pattern is square.

    Gauss reduction removes that freedom. Repeatedly subtract from the longer
    vector the integer multiple of the shorter that shortens it most, swapping
    when the order changes, until neither can be shortened. The result is the
    reduced basis, unique up to sign and exchange, and its lengths and included
    angle are then genuine invariants a user can read: the angle lands in
    [60, 120] degrees, and 90 means the lattice really is rectangular.

    The lattice itself is untouched — same nodes, same overlay, different labels.
    """

    first, second = np.array(basis[0], dtype=float), np.array(basis[1], dtype=float)
    for _ in range(64):
        if np.dot(first, first) > np.dot(second, second):
            first, second = second, first
        denominator = float(np.dot(first, first))
        if not denominator > 0.0 or not np.all(np.isfinite(second)):
            break
        step = round(float(np.dot(second, first)) / denominator)
        if step == 0:
            break
        second = second - step * first
    return np.ascontiguousarray(np.vstack([first, second]))


def _solve_basis_and_centre(
    positions: np.ndarray, indices: np.ndarray, *, refine_centre: bool
) -> tuple[np.ndarray, np.ndarray]:
    """Least squares for ``a``, ``b`` and optionally ``c`` in ``p = c + m a + n b``.

    Linear in every fitted parameter, the integers being known by this point.
    ``x`` and ``y`` share one design matrix and differ only in the right-hand
    side, which is why this is a single ``lstsq`` on a two-column right-hand side
    rather than a loop over coordinates.
    """

    if refine_centre:
        design = np.column_stack(
            [np.ones(len(indices)), indices[:, 0].astype(float), indices[:, 1].astype(float)]
        )
        solution, *_ = np.linalg.lstsq(design, positions, rcond=None)
        return np.ascontiguousarray(solution[1:3]), np.ascontiguousarray(solution[0])
    solution, *_ = np.linalg.lstsq(indices.astype(float), positions, rcond=None)
    return np.ascontiguousarray(solution), np.zeros(2)


def _converge(
    picked: np.ndarray,
    supplied_centre: np.ndarray,
    seed: np.ndarray,
    *,
    refine_centre: bool,
    leash_fraction: float,
    max_iterations: int,
) -> tuple[np.ndarray, np.ndarray, int, bool]:
    """Run assign-and-refine to convergence on one subset of the picks.

    Returns the basis, the centre, the number of rounds, and whether the centre
    reached its leash.
    """

    basis = seed
    centre = np.array(supplied_centre, dtype=float)
    assignment: np.ndarray | None = None
    iterations = 0
    leashed = False

    for round_index in range(int(max_iterations)):
        iterations = round_index + 1
        # A refinement can collapse a marginal basis onto a line. That is a dead
        # trial, not an error: the caller has other seeds to try, and one of them
        # is usually the right lattice.
        # `not (sine >= tol)` rather than `sine < tol`: a degenerate refinement
        # can produce a non-finite basis, and every comparison against NaN is
        # false, so the natural spelling lets exactly the worst case through and
        # the failure surfaces later as a NaN cast warning in unrelated code.
        sine = abs(float(np.linalg.det(basis))) / (
            float(np.linalg.norm(basis[0])) * float(np.linalg.norm(basis[1])) or 1.0
        )
        if not np.all(np.isfinite(basis)) or not sine >= _DEGENERATE_SINE:
            raise _DegenerateBasisError
        # `solve` rather than `inv`: these are the coordinates of the offsets in
        # the basis, and saying so keeps the intent visible. A near-singular
        # basis can return non-finite coordinates without raising, so the result
        # is checked rather than trusted — casting a NaN to int is undefined and
        # surfaces far away from here.
        try:
            coordinates = np.linalg.solve(basis.T, (picked - centre).T).T
        except np.linalg.LinAlgError as error:
            raise _DegenerateBasisError from error
        # Finite is not enough: a near-singular basis returns coordinates of
        # order 1e300, which are finite, cast to int as garbage, and warn about
        # it a long way from the cause. A lattice index of a million is already
        # absurd for a diffraction pattern.
        if not np.all(np.isfinite(coordinates)) or float(np.abs(coordinates).max()) > 1e6:
            raise _DegenerateBasisError
        next_assignment = np.rint(coordinates).astype(int)
        if assignment is not None and np.array_equal(next_assignment, assignment):
            break
        assignment = next_assignment
        if np.all(assignment == 0):
            # This trial basis is longer than the whole picked set, so every
            # spot rounds to the origin. Another candidate is shorter; only if
            # all of them fail is there something to tell the user.
            raise _DegenerateBasisError
        observations = picked if refine_centre else picked - centre
        basis, refined = _solve_basis_and_centre(
            observations, assignment, refine_centre=refine_centre
        )
        # Reduce here rather than at the end, so that "half a lattice spacing"
        # below means half of the *shortest* spacing rather than half of whatever
        # pair the arithmetic happened to produce. Same lattice either way.
        basis = _gauss_reduce(basis)
        if refine_centre:
            leash = leash_fraction * min(
                float(np.linalg.norm(basis[0])), float(np.linalg.norm(basis[1]))
            )
            offset = refined - supplied_centre
            distance = float(np.linalg.norm(offset))
            if distance > leash > 0.0:
                refined = supplied_centre + offset * (leash / distance)
                leashed = True
            centre = refined
        if not np.all(np.isfinite(centre)):
            raise _DegenerateBasisError

    return basis, centre, iterations, leashed


class _DegenerateBasisError(Exception):
    """A trial basis collapsed onto a line, so this seed has no fit to offer."""


@dataclass(frozen=True, slots=True)
class _Attempt:
    """One completed fit from one trial basis, ready to be ranked against others."""

    basis: np.ndarray
    centre: np.ndarray
    assignment: np.ndarray
    predicted: np.ndarray
    residuals: np.ndarray
    inliers: np.ndarray
    rms: float
    iterations: int
    leashed: bool

    @property
    def cell_area(self) -> float:
        """Area of the fitted unit cell."""

        return abs(float(np.linalg.det(self.basis)))

    def evidence(self, tolerance: float) -> float:
        """How surprising this fit is, if the picks were placed at random.

        Counting inliers alone cannot choose between lattices, because a finer
        lattice always explains at least as many points as a coarser one — halve
        the cell and every spot still lands on a node, along with a mis-picked
        spot and anything else. Inlier count therefore rewards exactly the models
        that assert the most and predict the least.

        The fix is to ask how *unlikely* each hit is. A cell of area ``A``
        scattered with nodes puts about ``pi t^2 / A`` of the plane within
        tolerance ``t`` of some node, so one inlier carries
        ``log(A / (pi t^2))`` nats of evidence and a whole fit carries that times
        its inlier count. A lattice half as coarse doubles its node density and
        pays ``log 4`` per inlier for the privilege, which one extra explained
        spot does not usually cover.

        The tolerance must be the *same* for every candidate or the comparison is
        empty — a tolerance proportional to each candidate's own basis makes the
        density ratio scale-invariant and the criterion blind to the very thing
        it is meant to detect. See :func:`fit_planar_lattice`, which derives one
        tolerance from the data and passes it to all of them.
        """

        density_ratio = self.cell_area / max(math.pi * tolerance * tolerance, 1e-12)
        return float(self.inliers.sum()) * math.log(max(density_ratio, 1.0 + 1e-12))


def _fit_from_seed(
    picked: np.ndarray,
    supplied_centre: np.ndarray,
    seed: np.ndarray,
    *,
    refine_centre: bool,
    tolerance: float,
    centre_leash_fraction: float,
    max_iterations: int,
) -> _Attempt:
    """Converge from one trial basis, trimming the worst spots as it goes.

    Least squares gives every point equal say, so one spot forty pixels out of
    place tilts the whole lattice and the largest residual can then land on an
    innocent neighbour rather than on the culprit. Refitting without the worst
    offender puts the blame where it belongs. Two removals is the limit: past
    that, the thing being called an outlier is the model.
    """

    keep = np.ones(picked.shape[0], dtype=bool)
    for _ in range(_MAX_TRIMMED + 1):
        basis, centre, iterations, leashed = _converge(
            picked[keep],
            supplied_centre,
            seed,
            refine_centre=refine_centre,
            leash_fraction=centre_leash_fraction,
            max_iterations=max_iterations,
        )
        coordinates = np.linalg.solve(basis.T, (picked - centre).T).T
        assignment = np.rint(coordinates).astype(int)
        predicted = centre + assignment.astype(float) @ basis
        residuals = np.linalg.norm(picked - predicted, axis=1)
        worst = int(np.argmax(np.where(keep, residuals, -1.0)))
        if residuals[worst] <= tolerance or int(keep.sum()) <= 4:
            break
        keep[worst] = False

    inliers = residuals <= tolerance
    rms = (
        float(np.sqrt(np.mean(residuals[inliers] ** 2))) if bool(np.any(inliers)) else float("inf")
    )
    return _Attempt(
        basis=basis,
        centre=centre,
        assignment=assignment,
        predicted=predicted,
        residuals=residuals,
        inliers=inliers,
        rms=rms,
        iterations=iterations,
        leashed=leashed,
    )


def fit_planar_lattice(
    positions: Any,
    centre: Any,
    *,
    refine_centre: bool = True,
    inlier_fraction: float = DEFAULT_INLIER_FRACTION,
    centre_leash_fraction: float = DEFAULT_CENTRE_LEASH_FRACTION,
    max_iterations: int = 8,
) -> PlanarLatticeFit:
    """Fit a 2D lattice to picked spots, optionally refining the beam centre.

    Purpose
    -------
    Impose the one constraint a zone-axis pattern satisfies before any
    crystallography is done — that its spots lie on a plane lattice — and use it
    to place the transmitted beam and to expose picks that do not belong.

    When and where to use it
    ------------------------
    Between picking and indexing, every time. The overlay it supports turns "the
    centre looks about right" into a judgement anyone can make, and the refined
    centre is over-determined by the picks rather than guessed from one of them.
    It is also the cheapest check that a set of picks is worth indexing at all: a
    set that will not fit a lattice will not index either, and finding that out
    here costs nothing.

    Parameters
    ----------
    positions : array_like
        ``(n, 2)`` picked spot coordinates, **not** including the transmitted
        beam. At least two, and not all collinear with the centre.
    centre : array_like
        The picked transmitted-beam position, in the same coordinates.
    refine_centre : bool
        Solve for the centre together with the basis. With four or more
        non-degenerate spots the system is over-determined and the answer beats
        any single pick. Turn it off to hold a centre established another way — a
        measured beam-stop position, for instance.
    inlier_fraction : float
        Residual tolerance as a fraction of the shorter basis vector. See
        :data:`DEFAULT_INLIER_FRACTION`.
    centre_leash_fraction : float
        How far the refinement may move the centre, as a fraction of the shorter
        basis vector. **A refinement may adjust the centre; it must not relabel
        which node the centre is.** Past half a lattice spacing a fit is no
        longer correcting the pick, it is choosing a different origin — and with
        few spots it will do exactly that, settling on a self-consistent
        assignment around the wrong node while the residuals stay small enough to
        look convincing. The leash stops that, and reaching it is reported in
        ``notes`` rather than hidden.
    max_iterations : int
        Cap on assign-and-refine rounds. Reaching it is not an error; the
        assignment simply had not settled, which ``iterations`` reports.

    Returns
    -------
    PlanarLatticeFit

    Raises
    ------
    ValueError
        If fewer than two spots are given, if they are all collinear with the
        centre, or if a bound is outside its admissible range.

    Notes
    -----
    Outliers are *labelled*, never removed from the reported output. A mis-picked
    spot the lattice rejects is the single most useful thing this can tell a user,
    and a fit that quietly dropped it would hide its own best finding. Trimming is
    different and narrower: the fit is *recomputed* without the worst one or two
    spots, so that ordinary least squares does not spread one bad pick across
    every residual and point at an innocent neighbour. Every spot, trimmed or
    not, is then scored against that cleaner fit and reported.

    Examples
    --------
    A square lattice with the centre picked two units off is recovered exactly::

        >>> import numpy as np
        >>> nodes = np.array([[10.0, 0.0], [0.0, 10.0], [-10.0, 0.0], [0.0, -10.0],
        ...                   [10.0, 10.0], [-10.0, -10.0]])
        >>> fit = fit_planar_lattice(nodes + np.array([100.0, 100.0]), (102.0, 99.0))
        >>> bool(np.allclose(fit.centre, [100.0, 100.0], atol=1e-9))
        True
        >>> round(fit.basis_angle_deg, 6)
        90.0
    """

    picked = as_float_array(np.asarray(positions, dtype=float).reshape(-1, 2), shape=(None, 2))
    supplied_centre = as_float_array(centre, shape=(2,))
    if picked.shape[0] < 2:
        raise ValueError("At least two spots are needed to define a two-dimensional lattice.")
    if not 0.0 < inlier_fraction <= 1.0:
        raise ValueError("inlier_fraction must lie in (0, 1].")
    if not 0.0 <= centre_leash_fraction <= 2.0:
        raise ValueError("centre_leash_fraction must lie in [0, 2].")
    if max_iterations <= 0:
        raise ValueError("max_iterations must be strictly positive.")

    # One tolerance, derived from the data and shared by every candidate. The
    # shortest difference between two picks is an upper bound on the lattice
    # spacing, so a fraction of it is a length in the same units as the picking
    # error — and, crucially, it does not move when a candidate lattice does.
    differences = (picked[:, None, :] - picked[None, :, :]).reshape(-1, 2)
    lengths = np.linalg.norm(differences, axis=1)
    shortest = float(lengths[lengths > 1e-12].min())
    tolerance = inlier_fraction * shortest

    attempts: list[_Attempt] = []
    for seed in _candidate_bases(picked, supplied_centre, _SEED_CANDIDATES):
        try:
            attempts.append(
                _fit_from_seed(
                    picked,
                    supplied_centre,
                    seed,
                    refine_centre=refine_centre,
                    tolerance=tolerance,
                    centre_leash_fraction=centre_leash_fraction,
                    max_iterations=max_iterations,
                )
            )
        except _DegenerateBasisError:
            continue
    if not attempts:
        raise ValueError(
            "No trial basis survived refinement: every candidate either collapsed onto a line or "
            "was longer than the whole picked set. The spots are too nearly collinear, or too few "
            "of them are distinct, to define a lattice. Pick a spot off the row you have."
        )
    # Rank by evidence, not by inlier count: see `_Attempt.evidence` for why
    # counting hits alone always rewards the finer lattice, and what replaces it.
    # Residual breaks ties between models that are equally surprising.
    attempts.sort(key=lambda item: (-item.evidence(tolerance), round(item.rms, 6)))
    best = attempts[0]
    # Safety net. The evidence criterion assumes some candidate is roughly right;
    # when none is, it can crown a coarse lattice that explains almost nothing.
    # Falling back to the most-explained fit keeps the overlay useful, and the
    # note says plainly that the lattice is not established.
    poorly_determined = float(best.inliers.sum()) < 0.6 * picked.shape[0]
    if poorly_determined:
        best = max(attempts, key=lambda item: (int(item.inliers.sum()), -round(item.rms, 6)))
    basis = best.basis
    working_centre = best.centre
    assignment = best.assignment
    predicted = best.predicted
    residuals = best.residuals
    inliers = best.inliers
    iterations = best.iterations
    leashed = best.leashed

    notes: list[str] = []
    if poorly_determined:
        notes.append(
            "No candidate lattice explains most of these picks, so the overlay is a best effort "
            "rather than a fit. Check the transmitted beam first, then look for a pick that is "
            "not a reflection; two spots on one row and nothing off it will also do this."
        )
    if leashed:
        notes.append(
            "The centre refinement reached its limit of half a lattice spacing from the picked "
            "position. Past that a fit stops correcting the pick and starts choosing a different "
            "node as the origin, so it is held here: nudge the beam closer and refine again."
        )

    # Two picks on one node cannot both be right, and the fit must not average
    # them into a plausible-looking answer. It happens when one reflection is
    # clicked twice, and when the centre is more than half a spacing out — the
    # second being exactly what a user is trying to fix, so it is reported rather
    # than quietly resolved.
    seen: dict[tuple[int, int], int] = {}
    duplicated: set[int] = set()
    for index in range(assignment.shape[0]):
        key = (int(assignment[index, 0]), int(assignment[index, 1]))
        if key in seen:
            duplicated.update({seen[key], index})
        else:
            seen[key] = index
    if duplicated:
        inliers = inliers & ~np.isin(np.arange(picked.shape[0]), sorted(duplicated))
        notes.append(
            "Spot(s) "
            + ", ".join(str(index + 1) for index in sorted(duplicated))
            + " were assigned to the same lattice node, so at most one of them can be right. "
            "Either one reflection was picked twice, or the transmitted beam is more than half a "
            "lattice spacing from where it was marked."
        )

    rms = (
        float(np.sqrt(np.mean(residuals[inliers] ** 2))) if bool(np.any(inliers)) else float("inf")
    )
    spots = tuple(
        LatticeFitSpot(
            index=int(index),
            position=picked[index],
            lattice_indices=(int(assignment[index, 0]), int(assignment[index, 1])),
            predicted=predicted[index],
            residual=float(residuals[index]),
            inlier=bool(inliers[index]),
        )
        for index in range(picked.shape[0])
    )
    return PlanarLatticeFit(
        centre=working_centre,
        supplied_centre=supplied_centre,
        basis=basis,
        spots=spots,
        rms_residual=rms,
        iterations=iterations,
        notes=tuple(notes),
    )
