"""Assigning Miller indices to measured reflections, and scoring the assignment.

Indexing is the step that turns a list of angles into crystallography. Until a
peak carries an ``(hkl)`` it contributes nothing to a lattice parameter, because
``sin^2(theta) = (lambda^2 / 4) h^T G* h`` needs the ``h`` as much as it needs
the angle.

This module does the *known-phase* case: a candidate phase is proposed, its
reflections are enumerated from its own cell and symmetry, and the measured
peaks are matched to them. That covers the overwhelming majority of laboratory
work -- a specimen whose phase is known or strongly suspected, measured to
determine its lattice parameter, its texture, or its stress. Determining an
unknown cell from scratch (ITO, TREOR, DICVOL and their successors) is a
different and much harder problem and is deliberately not attempted here.

Two decisions are worth stating.

**The assignment is global, not greedy.** Walking the peak list and taking the
nearest calculated line each time can assign two peaks to the same reflection,
or pair a peak with a near neighbour and leave the true partner for a worse
match further on. The Hungarian algorithm
(:func:`scipy.optimize.linear_sum_assignment`) instead minimizes the *total*
discrepancy over all one-to-one pairings at once, which is both the correct
formulation and, at these sizes, faster than the sorting a careful greedy pass
would need.

**The figures of merit are reported even when the phase was assumed.** de
Wolff's ``M_N`` and Smith and Snyder's ``F_N`` were designed to rank candidate
cells produced by an autoindexer, but their content is more general: they ask
how well a cell accounts for the observed lines *and* how much of its
explanatory power comes from simply predicting many lines. A cell that
"explains" every peak because it predicts a reflection everywhere has a poor
figure of merit, and that warning is as useful for a known phase as for an
unknown one.

References
----------
de Wolff, P. M., *J. Appl. Crystallogr.* **1** (1968) 108-113,
doi:10.1107/S002188986800508X -- the figure of merit ``M_20``.

Smith, G. S. & Snyder, R. L., *J. Appl. Crystallogr.* **12** (1979) 60-65,
doi:10.1107/S002188987901178X -- the figure of merit ``F_N``.

Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Ch. 10 -- indexing powder patterns of cubic and
non-cubic materials.

Kuhn, H. W., *Naval Res. Logist. Quart.* **2** (1955) 83-97,
doi:10.1002/nav.3800020109 -- the assignment algorithm.

Werner, P.-E., Eriksson, L. & Westdahl, M., *J. Appl. Crystallogr.* **18**
(1985) 367-370, doi:10.1107/S0021889885010512; Boultif, A. & Louer, D.,
*J. Appl. Crystallogr.* **37** (2004) 724-731,
doi:10.1107/S0021889804014876 -- the unknown-cell autoindexers this module
deliberately does not implement.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from pytex.core._arrays import as_float_array
from pytex.core.lattice import Phase
from pytex.diffraction.xrd import RadiationSpec, generate_powder_reflections
from pytex.diffraction.xrd_peaks import PeakFit, PeakTable

INDEXED_REFLECTION_SCHEMA = "pytex.diffraction.indexed_reflection"
PEAK_INDEXING_SCHEMA = "pytex.diffraction.peak_indexing"

_CITATION_DE_WOLFF = (
    "de Wolff, J. Appl. Crystallogr. 1 (1968) 108, doi:10.1107/S002188986800508X."
)
_CITATION_SMITH_SNYDER = (
    "Smith & Snyder, J. Appl. Crystallogr. 12 (1979) 60, doi:10.1107/S002188987901178X."
)
_CITATION_CULLITY = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 10."
)


def _triple(indices: Any) -> tuple[int, int, int]:
    """Return three Miller indices as a fixed-length tuple of Python ints.

    ``PowderReflection`` stores its representative as a NumPy integer array;
    the report type declares a three-tuple, and stating the length here keeps
    that declaration true rather than merely intended.
    """

    values = tuple(int(value) for value in indices)
    if len(values) != 3:
        raise ValueError(f"Expected three Miller indices, received {len(values)}.")
    return (values[0], values[1], values[2])


@dataclass(frozen=True, slots=True)
class IndexedReflection:
    """One measured peak paired with the reflection it was assigned to.

    Purpose
    -------
    Hold the pairing *and* its discrepancy, so that an indexing result can be
    audited line by line rather than trusted wholesale. The residual
    ``delta_two_theta_deg`` is the quantity a reader should scan first: a
    systematic sign pattern in it is a displaced specimen or a wrong zero, and
    a single large value is a misassignment or a second phase.

    Attributes
    ----------
    peak : PeakFit
        The measured reflection, carrying its own position uncertainty.
    miller_indices : tuple[int, int, int]
        The family representative ``(hkl)`` assigned to it.
    multiplicity : int
        Number of symmetry-equivalent planes in the family.
    two_theta_calculated_deg : float
        Where the candidate cell puts this reflection.
    d_observed_angstrom, d_calculated_angstrom : float
        Observed and calculated interplanar spacings.
    relative_intensity_calculated : float
        Kinematic relative intensity of the assigned family, for judging
        whether a strong observed peak was matched to a line that should have
        been invisible.
    """

    peak: PeakFit
    miller_indices: tuple[int, int, int]
    multiplicity: int
    two_theta_calculated_deg: float
    d_observed_angstrom: float
    d_calculated_angstrom: float
    relative_intensity_calculated: float

    def __post_init__(self) -> None:
        if len(self.miller_indices) != 3:
            raise ValueError("IndexedReflection.miller_indices must hold three integers.")
        if self.multiplicity <= 0:
            raise ValueError("IndexedReflection.multiplicity must be strictly positive.")
        if self.d_observed_angstrom <= 0.0 or self.d_calculated_angstrom <= 0.0:
            raise ValueError("IndexedReflection spacings must be strictly positive.")
        object.__setattr__(self, "miller_indices", _triple(self.miller_indices))

    @property
    def delta_two_theta_deg(self) -> float:
        """Return observed minus calculated position, in degrees."""

        return float(self.peak.two_theta_deg - self.two_theta_calculated_deg)

    @property
    def delta_q(self) -> float:
        """Return observed minus calculated ``Q = 1 / d^2``, in inverse square angstrom.

        ``Q`` rather than ``2 theta`` is the discrepancy measure of de Wolff's
        figure of merit, because ``Q`` is linear in the cell parameters and so
        weights a given cell error equally at every angle.
        """

        return float(
            1.0 / self.d_observed_angstrom**2 - 1.0 / self.d_calculated_angstrom**2
        )

    @property
    def normalized_residual(self) -> float:
        """Return the position residual in units of its own standard uncertainty.

        A value far above about three means the discrepancy is not explained by
        counting statistics, and therefore points at the cell, an uncorrected
        aberration, or a wrong assignment.
        """

        return float(
            self.delta_two_theta_deg / self.peak.two_theta_standard_uncertainty_deg
        )

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this pairing."""

        return {
            "schema": INDEXED_REFLECTION_SCHEMA,
            "miller_indices": list(self.miller_indices),
            "multiplicity": int(self.multiplicity),
            "two_theta_observed_deg": float(self.peak.two_theta_deg),
            "two_theta_calculated_deg": float(self.two_theta_calculated_deg),
            "two_theta_standard_uncertainty_deg": float(
                self.peak.two_theta_standard_uncertainty_deg
            ),
            "delta_two_theta_deg": self.delta_two_theta_deg,
            "normalized_residual": self.normalized_residual,
            "d_observed_angstrom": float(self.d_observed_angstrom),
            "d_calculated_angstrom": float(self.d_calculated_angstrom),
            "delta_q_inv_angstrom_squared": self.delta_q,
            "relative_intensity_calculated": float(self.relative_intensity_calculated),
        }


@dataclass(frozen=True, slots=True)
class PeakIndexing:
    """An assignment of measured peaks to the reflections of a candidate phase.

    Purpose
    -------
    Be the bridge between :class:`~pytex.diffraction.xrd_peaks.PeakTable` and
    lattice-parameter determination, and be explicit about what it did *not*
    explain: unindexed peaks and unobserved calculated lines are both carried,
    because both are evidence.

    Attributes
    ----------
    name : str
        A human name for the result.
    phase_name : str
        The candidate phase's name.
    reflections : tuple[IndexedReflection, ...]
        The accepted pairings, in ascending observed angle.
    unindexed_peaks : tuple[PeakFit, ...]
        Measured peaks with no calculated line inside the tolerance. A strong
        one is the signature of a second phase.
    unobserved_indices : tuple[tuple[int, int, int], ...]
        Calculated reflections above the intensity threshold that no peak was
        matched to. Systematic absence of a family can indicate texture, or a
        centring the candidate phase does not declare.
    tolerance_deg : float
        The matching tolerance the assignment was made under.
    radiation : RadiationSpec | None
    settings : Mapping[str, float | str]
        The settings that produced the result.
    """

    name: str
    phase_name: str
    reflections: tuple[IndexedReflection, ...]
    unindexed_peaks: tuple[PeakFit, ...] = ()
    unobserved_indices: tuple[tuple[int, int, int], ...] = ()
    tolerance_deg: float = 0.2
    radiation: RadiationSpec | None = None
    settings: Mapping[str, float | str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("PeakIndexing.name must be non-empty.")
        ordered = tuple(
            sorted(self.reflections, key=lambda item: item.peak.two_theta_deg)
        )
        object.__setattr__(self, "reflections", ordered)
        object.__setattr__(self, "unindexed_peaks", tuple(self.unindexed_peaks))
        object.__setattr__(
            self,
            "unobserved_indices",
            tuple(_triple(item) for item in self.unobserved_indices),
        )
        object.__setattr__(self, "settings", MappingProxyType(dict(self.settings)))

    def __len__(self) -> int:
        return len(self.reflections)

    def __iter__(self) -> Any:
        return iter(self.reflections)

    @property
    def indexed_count(self) -> int:
        """Return the number of peaks that received indices."""

        return len(self.reflections)

    @property
    def indexed_fraction(self) -> float:
        """Return the fraction of measured peaks that received indices."""

        total = len(self.reflections) + len(self.unindexed_peaks)
        return 0.0 if total == 0 else len(self.reflections) / total

    @property
    def mean_absolute_delta_two_theta_deg(self) -> float:
        """Return the mean absolute position discrepancy, in degrees."""

        if not self.reflections:
            return float("nan")
        return float(
            np.mean([abs(item.delta_two_theta_deg) for item in self.reflections])
        )

    @property
    def mean_absolute_delta_q(self) -> float:
        """Return the mean absolute ``Q`` discrepancy, in inverse square angstrom."""

        if not self.reflections:
            return float("nan")
        return float(np.mean([abs(item.delta_q) for item in self.reflections]))

    @property
    def delta_two_theta_deg(self) -> np.ndarray:
        """Return the signed position residuals as a read-only array."""

        return as_float_array(
            [item.delta_two_theta_deg for item in self.reflections], shape=(None,)
        )

    def figure_of_merit_m(self, *, count: int = 20) -> tuple[float, int]:
        """Return de Wolff's ``M_N`` and the ``N`` actually used.

        Purpose
        -------
        Score how convincingly the candidate cell explains the observed lines,
        penalizing a cell that succeeds only by predicting a great many of them.

        Method
        ------
        de Wolff (1968) defines

        ``M_N = Q_N / (2 <|Delta Q|> N_poss)``

        where ``Q = 1 / d^2``, ``Q_N`` is the ``Q`` of the ``N``-th observed
        line, ``<|Delta Q|>`` is the mean absolute discrepancy in ``Q`` over
        the first ``N`` lines, and ``N_poss`` is the number of *distinct*
        reflections the cell predicts below ``Q_N``. The last factor is the
        part that does the work: halving the discrepancies doubles ``M``, but
        so does a cell that predicts half as many lines while fitting equally
        well.

        Conventionally ``N = 20``, and ``M_20 > 10`` is taken as a plausible
        cell, ``M_20 > 20`` as a good one. When fewer than ``N`` lines are
        indexed the figure is computed over what there is and the ``N`` used is
        returned alongside it, because ``M_7`` and ``M_20`` are not comparable
        and reporting a bare number would invite treating them as if they were.

        Parameters
        ----------
        count
            The nominal ``N``. Clamped to the number of indexed lines.

        Returns
        -------
        tuple[float, int]
            The figure of merit and the ``N`` it was computed over. The figure
            is ``nan`` when nothing is indexed.
        """

        used = min(int(count), len(self.reflections))
        if used < 1:
            return (float("nan"), 0)
        subset = self.reflections[:used]
        q_limit = 1.0 / subset[-1].d_observed_angstrom ** 2
        mean_delta = float(np.mean([abs(item.delta_q) for item in subset]))
        possible = int(self.settings.get("possible_lines", used))
        if mean_delta <= 0.0 or possible < 1:
            return (float("inf"), used)
        return (float(q_limit / (2.0 * mean_delta * possible)), used)

    def figure_of_merit_f(self, *, count: int = 30) -> tuple[float, int]:
        """Return Smith and Snyder's ``F_N`` and the ``N`` actually used.

        Purpose
        -------
        Score the same thing as ``M_N`` but in the angular units an operator
        reads off the diffractometer, which makes it easier to sanity-check
        against the instrument's own resolution.

        Method
        ------
        Smith & Snyder (1979) define

        ``F_N = N / (<|Delta 2 theta|> N_poss)``

        with ``<|Delta 2 theta|>`` in degrees and ``N_poss`` again the number
        of calculated lines up to the ``N``-th observed one. It is normally
        quoted as ``F_N = value (mean discrepancy, N_poss)``, and this method's
        return value plus
        :attr:`mean_absolute_delta_two_theta_deg` supply both parts.

        Parameters
        ----------
        count
            The nominal ``N``. Clamped to the number of indexed lines.

        Returns
        -------
        tuple[float, int]
            The figure of merit and the ``N`` it was computed over.
        """

        used = min(int(count), len(self.reflections))
        if used < 1:
            return (float("nan"), 0)
        subset = self.reflections[:used]
        mean_delta = float(np.mean([abs(item.delta_two_theta_deg) for item in subset]))
        possible = int(self.settings.get("possible_lines", used))
        if mean_delta <= 0.0 or possible < 1:
            return (float("inf"), used)
        return (float(used / (mean_delta * possible)), used)

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this indexing."""

        merit_m, count_m = self.figure_of_merit_m()
        merit_f, count_f = self.figure_of_merit_f()
        return {
            "schema": PEAK_INDEXING_SCHEMA,
            "name": self.name,
            "phase_name": self.phase_name,
            "tolerance_deg": float(self.tolerance_deg),
            "indexed_count": self.indexed_count,
            "indexed_fraction": self.indexed_fraction,
            "mean_absolute_delta_two_theta_deg": self.mean_absolute_delta_two_theta_deg,
            "mean_absolute_delta_q_inv_angstrom_squared": self.mean_absolute_delta_q,
            "figure_of_merit_m": {"value": merit_m, "count": count_m},
            "figure_of_merit_f": {"value": merit_f, "count": count_f},
            "reflections": [item.to_json() for item in self.reflections],
            "unindexed_two_theta_deg": [
                float(peak.two_theta_deg) for peak in self.unindexed_peaks
            ],
            "unobserved_indices": [list(item) for item in self.unobserved_indices],
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this indexing."""

        if not self.reflections:
            return (
                f"Indexing '{self.name}' assigned no reflection of '{self.phase_name}' to any of "
                f"the {len(self.unindexed_peaks)} measured peaks within "
                f"{self.tolerance_deg:.3f} degrees. Either the phase is wrong, the tolerance is "
                f"too tight, or an uncorrected zero or displacement error exceeds it."
            )
        merit_m, count_m = self.figure_of_merit_m()
        merit_f, count_f = self.figure_of_merit_f()
        residuals = self.delta_two_theta_deg
        same_sign = bool(np.all(residuals > 0.0) or np.all(residuals < 0.0))
        systematic = (
            " Every residual carries the same sign, which is the signature of an uncorrected "
            "zero-point or specimen-displacement error rather than of a wrong cell; determining "
            "the lattice parameter with a refined systematic-error term will absorb it."
            if same_sign and len(residuals) >= 3
            else ""
        )
        leftover = (
            f" {len(self.unindexed_peaks)} measured peaks were not indexed; a strong one is the "
            "signature of a second phase."
            if self.unindexed_peaks
            else " Every measured peak was indexed."
        )
        missing = (
            f" {len(self.unobserved_indices)} calculated reflections above the intensity "
            "threshold were not observed, which can indicate texture or a centring the candidate "
            "phase does not declare."
            if self.unobserved_indices
            else ""
        )
        quality = (
            "which is a convincing account of the data"
            if merit_m > 20.0
            else (
                "which is plausible but not conclusive"
                if merit_m > 10.0
                else "which is weak: the cell explains the lines no better than a cell that simply "
                "predicts many of them would"
            )
        )
        return (
            f"Indexing '{self.name}' assigned {self.indexed_count} of "
            f"{self.indexed_count + len(self.unindexed_peaks)} measured peaks to reflections of "
            f"'{self.phase_name}' by global one-to-one assignment within "
            f"{self.tolerance_deg:.3f} degrees 2*theta. The mean absolute position discrepancy is "
            f"{self.mean_absolute_delta_two_theta_deg:.5f} degrees. de Wolff's figure of merit is "
            f"M_{count_m} = {merit_m:.1f}, {quality}; Smith and Snyder's is "
            f"F_{count_f} = {merit_f:.1f}.{systematic}{leftover}{missing} "
            f"{_CITATION_DE_WOLFF} {_CITATION_SMITH_SNYDER} {_CITATION_CULLITY}"
        )


def index_peaks(
    table: PeakTable,
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    tolerance_deg: float = 0.3,
    max_index: int = 6,
    minimum_relative_intensity: float = 0.001,
    phase_name: str | None = None,
    name: str | None = None,
) -> PeakIndexing:
    """Assign Miller indices to measured peaks from a candidate phase.

    Purpose
    -------
    Turn fitted angles into indexed reflections, which is the input every
    lattice-parameter method in this library requires, and report how well the
    candidate cell actually accounted for the pattern.

    Method
    ------
    1. The phase's reflection families are enumerated over the measured angular
       range by
       :func:`~pytex.diffraction.xrd.generate_powder_reflections`, which applies
       the symmetry and the systematic absences of the phase rather than a
       generic ``(hkl)`` list. Families weaker than
       ``minimum_relative_intensity`` are dropped, because matching a strong
       observed peak to a line that should be invisible is not an explanation.
    2. A cost matrix of ``|2 theta_obs - 2 theta_calc|`` is formed, with pairs
       beyond ``tolerance_deg`` made prohibitively expensive.
    3. The Hungarian algorithm minimizes the total cost over all one-to-one
       pairings. This is global: a greedy nearest-line pass can assign two
       peaks to the same reflection or pair a peak with a near neighbour and
       strand the true partner, and neither failure is visible in the result it
       produces.
    4. Pairs whose discrepancy still exceeds the tolerance are rejected, and
       the peaks and reflections involved are reported as unindexed and
       unobserved.

    Parameters
    ----------
    table
        Fitted peaks. Use
        :meth:`~pytex.diffraction.xrd_peaks.PeakTable.filter_converged` first if
        any fit failed.
    phase
        The candidate phase, supplying cell, symmetry and -- for the intensity
        filter -- the atomic basis.
    radiation
        Radiation to calculate positions with. Falls back to the table's.
    tolerance_deg
        Maximum accepted ``|2 theta_obs - 2 theta_calc|``. It must be wider
        than any uncorrected zero or displacement error, and narrower than the
        spacing between neighbouring calculated lines; the default of 0.3
        degrees suits a laboratory scan of a well-aligned instrument.
    max_index
        Largest ``|h|``, ``|k|``, ``|l|`` enumerated.
    minimum_relative_intensity
        Drop calculated families weaker than this fraction of the strongest.
    phase_name, name
        Names for the report.

    Returns
    -------
    PeakIndexing
        The pairing, its residuals, and its figures of merit.

    Raises
    ------
    ValueError
        If the table is empty, carries no radiation and none was supplied, the
        tolerance is not positive, or the phase predicts no reflection in the
        measured range.

    See Also
    --------
    pytex.diffraction.xrd_lattice_parameter.determine_lattice_parameters :
        consumes this result.

    Notes
    -----
    An indexing that succeeds is not by itself evidence that the phase is
    right. Read the residual column: a run of same-signed residuals is an
    instrument error, not a cell error, and
    :meth:`PeakIndexing.describe` says so when it sees one.
    """

    if len(table) == 0:
        raise ValueError("index_peaks was given an empty peak table.")
    if not np.isfinite(tolerance_deg) or tolerance_deg <= 0.0:
        raise ValueError("index_peaks requires a finite, positive tolerance_deg.")
    spec = radiation if radiation is not None else table.radiation
    if spec is None:
        raise ValueError(
            "index_peaks needs a radiation: the peak table declares none and none was passed."
        )

    observed = table.two_theta_deg
    margin = max(2.0 * float(tolerance_deg), 1.0)
    low = max(float(observed[0]) - margin, 1.0e-3)
    high = min(float(observed[-1]) + margin, 179.999)
    calculated = generate_powder_reflections(
        phase,
        radiation=spec,
        two_theta_range_deg=(low, high),
        max_index=max_index,
    )
    if not calculated:
        raise ValueError(
            f"The candidate phase predicts no reflection between {low:.2f} and {high:.2f} "
            "degrees 2*theta. Check the phase, the radiation and max_index."
        )
    strongest = max(item.intensity for item in calculated)
    if strongest > 0.0:
        calculated = tuple(
            item
            for item in calculated
            if item.intensity >= minimum_relative_intensity * strongest
        )
    if not calculated:
        raise ValueError(
            "Every predicted reflection fell below minimum_relative_intensity; lower it or "
            "check that the phase carries an atomic basis."
        )

    calculated_angles = np.array([item.two_theta_deg for item in calculated])
    cost = np.abs(observed[:, None] - calculated_angles[None, :])
    # Forbidden pairs must be expensive rather than infinite: linear_sum_assignment
    # requires a feasible complete assignment on the rectangular problem, and an
    # infinity anywhere in a row makes that impossible to express.
    penalty = float(cost.max()) + 1.0e3
    cost = np.where(cost <= float(tolerance_deg), cost, penalty)

    peak_rows, reflection_columns = linear_sum_assignment(cost)

    wavelength = float(spec.wavelength_angstrom)
    accepted: list[IndexedReflection] = []
    matched_peaks: set[int] = set()
    matched_reflections: set[int] = set()
    for row, column in zip(peak_rows, reflection_columns, strict=True):
        if cost[row, column] >= penalty:
            continue
        peak = table.peaks[int(row)]
        reflection = calculated[int(column)]
        accepted.append(
            IndexedReflection(
                peak=peak,
                miller_indices=_triple(reflection.miller_indices),
                multiplicity=int(reflection.multiplicity),
                two_theta_calculated_deg=float(reflection.two_theta_deg),
                d_observed_angstrom=peak.d_spacing_angstrom(wavelength),
                d_calculated_angstrom=float(reflection.d_spacing_angstrom),
                relative_intensity_calculated=float(reflection.intensity / strongest)
                if strongest > 0.0
                else 0.0,
            )
        )
        matched_peaks.add(int(row))
        matched_reflections.add(int(column))

    unindexed = tuple(
        peak for index, peak in enumerate(table.peaks) if index not in matched_peaks
    )
    unobserved = tuple(
        _triple(item.miller_indices)
        for index, item in enumerate(calculated)
        if index not in matched_reflections
    )

    # N_poss for the figures of merit: distinct calculated lines up to the
    # highest indexed observed line, which is what both definitions specify.
    if accepted:
        highest = max(item.peak.two_theta_deg for item in accepted)
        possible = int(np.count_nonzero(calculated_angles <= highest))
    else:
        possible = len(calculated)

    return PeakIndexing(
        name=name or f"{table.name} indexed",
        phase_name=phase_name or str(getattr(phase, "name", None) or "candidate phase"),
        reflections=tuple(accepted),
        unindexed_peaks=unindexed,
        unobserved_indices=unobserved,
        tolerance_deg=float(tolerance_deg),
        radiation=spec,
        settings={
            "possible_lines": float(max(possible, 1)),
            "max_index": float(max_index),
            "minimum_relative_intensity": float(minimum_relative_intensity),
            "source_table": table.name,
        },
    )


__all__ = [
    "INDEXED_REFLECTION_SCHEMA",
    "PEAK_INDEXING_SCHEMA",
    "IndexedReflection",
    "PeakIndexing",
    "index_peaks",
]
