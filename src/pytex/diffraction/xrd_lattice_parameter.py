"""Precise lattice-parameter *determination* from a measured powder pattern.

This module answers one question: given a diffractogram of a phase whose
structure is known, what is its unit cell, and how well is it known? It is not
a structure refinement. Nothing here varies an atomic coordinate, a thermal
parameter, or a site occupancy; the structure is held fixed and only the cell
and the instrument's own errors are determined. That restriction is the point,
because it is what stops texture and an imperfect structural model from leaking
into the answer.

Why this needs its own module
-----------------------------
The obvious approach -- compute ``a`` from each reflection and average -- is
worth understanding precisely because it fails, and it fails for a reason no
amount of extra data repairs. Differentiating Bragg's law gives

``Delta d / d = -cot(theta) Delta theta``

so a fixed angular error produces a *theta-dependent* error in the spacing. The
errors that dominate a laboratory scan are systematic, not random: a detector
zero offset, a specimen a few tens of micrometres off the diffractometer axis,
a beam that penetrates before it diffracts. Averaging over reflections divides
the random scatter by ``sqrt(N)`` and leaves the systematic part untouched. The
result is a lattice parameter good to about one part in ``10^3`` -- roughly two
orders of magnitude short of what stress and thermal-expansion work needs.

The methods offered here
------------------------
``"average"``
    The naive per-reflection mean, offered deliberately so the comparison
    against the others can be made on the reader's own data in one click. It
    exists only for cubic cells, and *that* is instructive too: outside the
    cubic system a "lattice parameter per reflection" is not defined at all,
    because one reflection cannot determine two cell parameters.

``"cohen"`` (default)
    Weighted linear least squares in ``sin^2(theta)``, with a systematic-error
    coefficient refined alongside the cell. This is Cullity's Chapter 11
    treatment generalized: because

    ``sin^2(theta) = (lambda^2 / 4) h^T G* h``

    is *linear* in the components of the reciprocal metric tensor ``G*``, one
    code path covers cubic through triclinic with no starting guess, no local
    minima, and an analytic covariance matrix that yields real standard
    uncertainties on every cell parameter.

``"le_bail"``
    Whole-pattern decomposition. Every measured point contributes, reflection
    intensities are *extracted* rather than modelled, and overlapped
    reflections -- the normal situation in a hexagonal pattern -- are handled
    properly. Because the intensities are free, neither texture nor a wrong
    atomic basis can bias the cell.

The systematic-error term
-------------------------
The key identity that unifies the classical graphical extrapolations with the
least-squares treatment: if a systematic aberration produces a fractional
spacing error ``Delta d / d = -K f(theta)`` for some extrapolation function
``f``, then since ``Delta(sin^2 theta) / sin^2(theta) = -2 Delta d / d``,

``sin^2(theta)_observed = (lambda^2 / 4) h^T G* h + D sin^2(theta) f(theta)``

with ``D = 2 K``. So *whatever function you would plot ``a`` against in the
classical method is the same function that appears, multiplied by
``sin^2(theta)``, as a design column here.* Taking ``f = cos^2(theta)``
reproduces Cohen's classical ``sin^2(2 theta)`` drift column exactly, since
``sin^2(theta) cos^2(theta) = sin^2(2 theta) / 4``. Because every offered
``f`` vanishes at ``theta = 90`` degrees, the fitted cell is the extrapolated
one.

What is deliberately not here
-----------------------------
Refining the detector zero and the specimen displacement *together* from a
single specimen scan. They are separable in principle -- one is constant and
the other goes as ``cos(theta)`` -- and badly correlated in practice over the
angular range of one pattern. Zero belongs to a calibrated instrument
(:class:`~pytex.diffraction.xrd_instrument.InstrumentBroadening` carries one);
displacement belongs to the specimen.

Also not here: the sin-squared-psi analysis that residual-stress measurement
actually requires. A single symmetric scan determines the spacing along one
direction only. This module supplies the precise spacings that analysis
consumes; it does not compute a stress, and
:meth:`LatticeParameterResult.describe` says so.

References
----------
Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Ch. 11 "Precise Parameter Measurements" -- the
extrapolation functions, Cohen's method, and the precision each reaches.

Nelson, J. B. & Riley, D. P., *Proc. Phys. Soc.* **57** (1945) 160-177,
doi:10.1088/0959-5309/57/3/302 -- the extrapolation function.

Bradley, A. J. & Jay, A. H., *Proc. Phys. Soc.* **44** (1932) 563-579,
doi:10.1088/0959-5309/44/5/305 -- the ``cos^2(theta)`` extrapolation.

Cohen, M. U., *Rev. Sci. Instrum.* **6** (1935) 68-74,
doi:10.1063/1.1751937 -- least squares with a drift term in place of a
graphical extrapolation.

Le Bail, A., Duroy, H. & Fourquet, J. L., *Mater. Res. Bull.* **23** (1988)
447-452, doi:10.1016/0025-5408(88)90019-0 -- whole-pattern intensity
extraction.

Pawley, G. S., *J. Appl. Crystallogr.* **14** (1981) 357-361,
doi:10.1107/S0021889881009618 -- the alternative whole-pattern method.

Wilson, A. J. C., *Mathematical Theory of X-ray Powder Diffractometry*, Philips
Technical Library (1963) -- the aberration forms the drift term absorbs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
from scipy.optimize import least_squares

from pytex.core._arrays import as_float_array
from pytex.core.lattice import Lattice, Phase
from pytex.diffraction.xrd import RadiationSpec, generate_powder_reflections
from pytex.diffraction.xrd_background import estimate_background
from pytex.diffraction.xrd_indexing import PeakIndexing, index_peaks
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import (
    PeakTable,
    detect_and_fit_peaks,
    pseudo_voigt_profile,
)

LATTICE_PARAMETER_SCHEMA = "pytex.diffraction.lattice_parameter_result"

LatticeMethod = Literal["cohen", "average", "le_bail"]
ExtrapolationFunction = Literal["nelson_riley", "bradley_jay", "cos_squared_over_sin", "none"]
SystematicTerm = Literal["zero", "displacement", "none"]

LATTICE_METHODS: tuple[LatticeMethod, ...] = ("cohen", "average", "le_bail")
EXTRAPOLATION_FUNCTIONS: tuple[ExtrapolationFunction, ...] = (
    "nelson_riley",
    "bradley_jay",
    "cos_squared_over_sin",
    "none",
)
SYSTEMATIC_TERMS: tuple[SystematicTerm, ...] = ("zero", "displacement", "none")

_CITATION_CULLITY = (
    "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11."
)
_CITATION_NELSON_RILEY = (
    "Nelson & Riley, Proc. Phys. Soc. 57 (1945) 160, doi:10.1088/0959-5309/57/3/302."
)
_CITATION_COHEN = "Cohen, Rev. Sci. Instrum. 6 (1935) 68, doi:10.1063/1.1751937."
_CITATION_LE_BAIL = (
    "Le Bail, Duroy & Fourquet, Mater. Res. Bull. 23 (1988) 447, "
    "doi:10.1016/0025-5408(88)90019-0."
)

#: Free reciprocal-metric-tensor parameters per crystal system, as the names
#: that appear in a report and the 6-by-n matrix that expands them into the
#: six independent components ``(G*11, G*22, G*33, G*12, G*13, G*23)``.
#:
#: The hexagonal row is the one worth reading: ``a* = b*`` and ``gamma* = 60``
#: degrees, so ``G*12 = a* b* cos(gamma*) = G*11 / 2``, and the quadratic form
#: collapses to ``A (h^2 + h k + k^2) + C l^2`` -- the familiar hexagonal
#: expression, obtained here as a constraint rather than as a special case.
_SYSTEM_CONSTRAINTS: dict[str, tuple[tuple[str, ...], tuple[tuple[float, ...], ...]]] = {
    "cubic": (
        ("a*^2",),
        ((1.0,), (1.0,), (1.0,), (0.0,), (0.0,), (0.0,)),
    ),
    "tetragonal": (
        ("a*^2", "c*^2"),
        ((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)),
    ),
    "hexagonal": (
        ("a*^2", "c*^2"),
        ((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.0), (0.0, 0.0), (0.0, 0.0)),
    ),
    "trigonal": (
        ("a*^2", "c*^2"),
        ((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.0), (0.0, 0.0), (0.0, 0.0)),
    ),
    "orthorhombic": (
        ("a*^2", "b*^2", "c*^2"),
        (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ),
    ),
    "monoclinic": (
        ("a*^2", "b*^2", "c*^2", "a*c*cos(beta*)"),
        (
            (1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
            (0.0, 0.0, 0.0, 0.0),
        ),
    ),
    "triclinic": (
        (
            "a*^2",
            "b*^2",
            "c*^2",
            "a*b*cos(gamma*)",
            "a*c*cos(beta*)",
            "b*c*cos(alpha*)",
        ),
        tuple(tuple(1.0 if row == column else 0.0 for column in range(6)) for row in range(6)),
    ),
}


def crystal_system_of(phase: Phase) -> str:
    """Return the crystal system that constrains a phase's cell.

    Purpose
    -------
    Decide how many cell parameters a determination may vary, from the phase's
    own symmetry rather than from a guess based on the numerical values of its
    cell edges.

    Method
    ------
    The point group derived from the phase's symmetry names the system
    directly. A cell whose edges happen to be equal is *not* cubic unless its
    symmetry says so, and determining it as though it were would impose a
    constraint the crystal does not obey.

    Parameters
    ----------
    phase
        The phase whose system is wanted.

    Returns
    -------
    str
        One of the seven crystal systems, lower case.

    Raises
    ------
    ValueError
        If the phase's symmetry does not resolve to a supported system.
    """

    system = str(phase.symmetry.to_point_group().crystal_system).strip().lower()
    if system not in _SYSTEM_CONSTRAINTS:
        raise ValueError(
            f"Crystal system '{system}' has no cell parameterization here. Supported systems: "
            + ", ".join(sorted(_SYSTEM_CONSTRAINTS))
            + "."
        )
    return system


def _constraint_matrix(system: str) -> tuple[tuple[str, ...], np.ndarray]:
    names, rows = _SYSTEM_CONSTRAINTS[system]
    return names, np.asarray(rows, dtype=np.float64)


def _quadratic_rows(indices: np.ndarray) -> np.ndarray:
    """Return the ``[h^2, k^2, l^2, 2hk, 2hl, 2kl]`` rows of a set of indices.

    This is ``h^T G* h`` written as a dot product with the six independent
    components of ``G*``, which is what makes the whole determination linear.
    """

    h = indices[:, 0]
    k = indices[:, 1]
    l_index = indices[:, 2]
    return np.column_stack(
        [
            h * h,
            k * k,
            l_index * l_index,
            2.0 * h * k,
            2.0 * h * l_index,
            2.0 * k * l_index,
        ]
    ).astype(np.float64)


def _reciprocal_tensor(components: np.ndarray) -> np.ndarray:
    """Assemble a symmetric ``G*`` from its six independent components."""

    g11, g22, g33, g12, g13, g23 = (float(value) for value in components)
    return np.array(
        [[g11, g12, g13], [g12, g22, g23], [g13, g23, g33]], dtype=np.float64
    )


def _cell_from_reciprocal_tensor(components: np.ndarray) -> tuple[float, ...]:
    """Return ``(a, b, c, alpha, beta, gamma)`` from the six ``G*`` components.

    The direct metric tensor is the inverse of the reciprocal one, and the
    direct cell reads straight off it: ``a = sqrt(G11)`` and
    ``cos(alpha) = G23 / (b c)``. Going through the tensor rather than through
    per-system algebra means one inversion serves every crystal system.
    """

    reciprocal = _reciprocal_tensor(components)
    determinant = float(np.linalg.det(reciprocal))
    if not np.isfinite(determinant) or determinant <= 0.0:
        raise ValueError(
            "The determined reciprocal metric tensor is not positive definite, so it does not "
            "describe a real cell. The indexing, the tolerance, or the angular range is at fault."
        )
    direct = np.linalg.inv(reciprocal)
    lengths = np.sqrt(np.diag(direct))
    if np.any(~np.isfinite(lengths)) or np.any(lengths <= 0.0):
        raise ValueError("The determined cell has a non-positive edge length.")
    cos_alpha = direct[1, 2] / (lengths[1] * lengths[2])
    cos_beta = direct[0, 2] / (lengths[0] * lengths[2])
    cos_gamma = direct[0, 1] / (lengths[0] * lengths[1])
    angles = np.rad2deg(np.arccos(np.clip([cos_alpha, cos_beta, cos_gamma], -1.0, 1.0)))
    return (
        float(lengths[0]),
        float(lengths[1]),
        float(lengths[2]),
        float(angles[0]),
        float(angles[1]),
        float(angles[2]),
    )


def extrapolation_values(
    two_theta_deg: Any, *, function: ExtrapolationFunction = "nelson_riley"
) -> np.ndarray:
    """Return the extrapolation function ``f(theta)`` at each angle.

    Purpose
    -------
    Supply the abscissa of the classical ``a`` versus ``f(theta)`` plot, and --
    multiplied by ``sin^2(theta)`` -- the drift column of the least-squares
    determination. Having one implementation serve both is what makes the
    graphical and algebraic treatments provably the same method.

    Method
    ------
    Every function below tends to zero as ``theta`` tends to 90 degrees, which
    is the whole point: at back-reflection the geometric aberrations vanish, so
    extrapolating there removes them.

    * ``"nelson_riley"``:
      ``f = (cos^2(theta) / sin(theta) + cos^2(theta) / theta) / 2``, with
      ``theta`` in radians. The standard choice, empirically linear over a wide
      range because it approximates the sum of the absorption and displacement
      aberrations.
    * ``"bradley_jay"``: ``f = cos^2(theta)``. Simpler, and adequate when only
      high-angle reflections are used. Multiplied by ``sin^2(theta)`` it gives
      Cohen's classical ``sin^2(2 theta) / 4`` drift column.
    * ``"cos_squared_over_sin"``: ``f = cos^2(theta) / sin(theta)``. This is
      the exact form for specimen displacement on a Bragg-Brentano
      diffractometer, and is the right choice when displacement is known to
      dominate.
    * ``"none"``: ``f = 0``. No systematic correction, which is worth running
      to see how much the correction was worth.

    Parameters
    ----------
    two_theta_deg
        Angles in degrees.
    function
        One of :data:`EXTRAPOLATION_FUNCTIONS`.

    Returns
    -------
    np.ndarray
        ``f(theta)``, dimensionless, same shape as the input.

    Raises
    ------
    ValueError
        If the function is unknown, or an angle is at 0 or 180 degrees where
        the ``1 / sin`` forms diverge.
    """

    if function not in EXTRAPOLATION_FUNCTIONS:
        raise ValueError(f"extrapolation_values requires function in {EXTRAPOLATION_FUNCTIONS}.")
    angles = np.asarray(two_theta_deg, dtype=np.float64)
    theta = np.deg2rad(0.5 * angles)
    if function == "none":
        return np.zeros_like(theta)
    if function == "bradley_jay":
        return np.square(np.cos(theta))
    if np.any(np.sin(theta) < 1.0e-9) or np.any(theta < 1.0e-9):
        raise ValueError(
            "The Nelson-Riley and cos^2/sin extrapolation functions diverge at 2*theta = 0; "
            "restrict the angular range."
        )
    if function == "cos_squared_over_sin":
        return np.square(np.cos(theta)) / np.sin(theta)
    return 0.5 * (
        np.square(np.cos(theta)) / np.sin(theta) + np.square(np.cos(theta)) / theta
    )


@dataclass(frozen=True, slots=True)
class LatticeParameterResult:
    """A determined unit cell, its uncertainties, and how it was obtained.

    Purpose
    -------
    Report a lattice parameter the way a measurement must be reported: with a
    standard uncertainty, with the systematic correction that was applied and
    how much it mattered, and with enough of the residual structure visible
    that a reader can judge whether to believe it.

    Attributes
    ----------
    method : str
        ``"cohen"``, ``"average"`` or ``"le_bail"``.
    phase_name : str
        The phase whose cell was determined.
    crystal_system : str
        The system that constrained the free parameters.
    a, b, c : float
        Determined cell edges in angstrom.
    alpha_deg, beta_deg, gamma_deg : float
        Determined cell angles in degrees.
    a_standard_uncertainty, b_standard_uncertainty, c_standard_uncertainty : float
        Standard uncertainties on the edges, in angstrom, propagated from the
        fit covariance.
    free_parameter_names : tuple[str, ...]
        The reciprocal-metric-tensor parameters that were varied.
    extrapolation : str
        The systematic-error function used, or ``"none"``.
    drift_coefficient, drift_standard_uncertainty : float
        The fitted ``D`` and its uncertainty. ``D`` is dimensionless and its
        physical meaning is direct: the fractional spacing error it removed at
        angle ``theta`` is ``-D f(theta) / 2``.
    reflection_count : int
        Number of indexed reflections used (or, for Le Bail, modelled).
    reduced_chi_squared : float
        Weighted residual sum of squares per degree of freedom.
    residual_two_theta_deg : np.ndarray
        Observed minus recalculated positions after determination, in degrees.
        Structure left here is structure the model did not describe.
    miller_indices : tuple[tuple[int, int, int], ...]
        The reflections behind ``residual_two_theta_deg``, in the same order.
    two_theta_deg : np.ndarray
        The observed positions those residuals belong to.
    reference_lattice : Lattice | None
        The starting cell, retained so a strain can be quoted against it.
    settings : Mapping[str, float | str]
        The settings the determination ran with.
    """

    method: LatticeMethod
    phase_name: str
    crystal_system: str
    a: float
    b: float
    c: float
    alpha_deg: float
    beta_deg: float
    gamma_deg: float
    a_standard_uncertainty: float
    b_standard_uncertainty: float
    c_standard_uncertainty: float
    free_parameter_names: tuple[str, ...]
    extrapolation: ExtrapolationFunction
    drift_coefficient: float
    drift_standard_uncertainty: float
    reflection_count: int
    reduced_chi_squared: float
    residual_two_theta_deg: np.ndarray
    miller_indices: tuple[tuple[int, int, int], ...]
    two_theta_deg: np.ndarray
    reference_lattice: Lattice | None = None
    settings: Mapping[str, float | str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.method not in LATTICE_METHODS:
            raise ValueError(f"LatticeParameterResult.method must be one of {LATTICE_METHODS}.")
        if self.extrapolation not in EXTRAPOLATION_FUNCTIONS:
            raise ValueError(
                f"LatticeParameterResult.extrapolation must be one of {EXTRAPOLATION_FUNCTIONS}."
            )
        for name in ("a", "b", "c"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"Determined cell edge '{name}' must be strictly positive.")
        for name in ("alpha_deg", "beta_deg", "gamma_deg"):
            if not 0.0 < getattr(self, name) < 180.0:
                raise ValueError(f"Determined cell angle '{name}' must lie in (0, 180).")
        if self.reflection_count < 1:
            raise ValueError("A determination needs at least one reflection.")
        residuals = as_float_array(self.residual_two_theta_deg, shape=(None,))
        positions = as_float_array(self.two_theta_deg, shape=(None,))
        if residuals.shape != positions.shape:
            raise ValueError("Residuals must align with the positions they belong to.")
        if len(self.miller_indices) != residuals.size:
            raise ValueError("Residuals must align with their Miller indices.")
        object.__setattr__(self, "residual_two_theta_deg", residuals)
        object.__setattr__(self, "two_theta_deg", positions)
        object.__setattr__(self, "settings", MappingProxyType(dict(self.settings)))

    @property
    def relative_uncertainty(self) -> float:
        """Return ``sigma(a) / a``, the figure that decides whether this is precise.

        A naive per-reflection average reaches about ``1e-3``; a Cohen
        determination or a Le Bail fit with a refined drift term reaches about
        ``1e-5``. Elastic strains of engineering interest are ``1e-4`` to
        ``1e-3``, so the difference is the difference between measuring a
        strain and not.
        """

        return float(self.a_standard_uncertainty / self.a)

    @property
    def axial_ratio(self) -> float:
        """Return ``c / a``, the quantity hexagonal and tetragonal work turns on."""

        return float(self.c / self.a)

    @property
    def strain_relative_to_reference(self) -> float | None:
        """Return ``(a - a_ref) / a_ref`` against the starting cell, if there is one.

        This is the lattice strain along the ``a`` axis relative to whatever
        cell the phase was defined with. It is a *strain*, not a stress: a
        symmetric theta-2theta scan measures the spacing of planes parallel to
        the specimen surface only, and converting that to a stress requires
        measurements at several specimen tilts and the X-ray elastic constants
        of the reflection used.
        """

        if self.reference_lattice is None:
            return None
        return float((self.a - self.reference_lattice.a) / self.reference_lattice.a)

    @property
    def systematic_shift_deg(self) -> np.ndarray:
        """Return the angular shift the drift term removed at each reflection.

        Purpose
        -------
        Make the correction auditable. A drift coefficient is abstract; the
        number of millidegrees it moved each peak is not, and comparing it with
        the position uncertainties says immediately whether refining it was
        worth doing.
        """

        if self.extrapolation == "none" or self.drift_coefficient == 0.0:
            return np.zeros_like(self.two_theta_deg)
        function = extrapolation_values(self.two_theta_deg, function=self.extrapolation)
        theta = np.deg2rad(0.5 * self.two_theta_deg)
        # Delta(sin^2 theta) = D sin^2(theta) f, and d(sin^2 theta)/d(2 theta)
        # = sin(2 theta) / 2, so the angular equivalent follows directly.
        delta_sin_squared = self.drift_coefficient * np.square(np.sin(theta)) * function
        shift: np.ndarray = np.rad2deg(2.0 * delta_sin_squared / np.sin(2.0 * theta))
        return shift

    def to_lattice(self) -> Lattice:
        """Return the determined cell as a :class:`~pytex.core.lattice.Lattice`.

        Raises
        ------
        ValueError
            If no reference lattice is carried, since the crystal frame must
            come from somewhere and inventing one would be a silent
            convention choice.
        """

        if self.reference_lattice is None:
            raise ValueError(
                "to_lattice needs the reference lattice for its crystal frame; this result was "
                "built without one."
            )
        return replace(
            self.reference_lattice,
            a=self.a,
            b=self.b,
            c=self.c,
            alpha_deg=self.alpha_deg,
            beta_deg=self.beta_deg,
            gamma_deg=self.gamma_deg,
        )

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable contract for this determination."""

        return {
            "schema": LATTICE_PARAMETER_SCHEMA,
            "method": self.method,
            "phase_name": self.phase_name,
            "crystal_system": self.crystal_system,
            "cell": {
                "a": float(self.a),
                "b": float(self.b),
                "c": float(self.c),
                "alpha_deg": float(self.alpha_deg),
                "beta_deg": float(self.beta_deg),
                "gamma_deg": float(self.gamma_deg),
            },
            "standard_uncertainty": {
                "a": float(self.a_standard_uncertainty),
                "b": float(self.b_standard_uncertainty),
                "c": float(self.c_standard_uncertainty),
            },
            "relative_uncertainty": self.relative_uncertainty,
            "axial_ratio": self.axial_ratio,
            "free_parameter_names": list(self.free_parameter_names),
            "extrapolation": self.extrapolation,
            "drift_coefficient": float(self.drift_coefficient),
            "drift_standard_uncertainty": float(self.drift_standard_uncertainty),
            "reflection_count": int(self.reflection_count),
            "reduced_chi_squared": float(self.reduced_chi_squared),
            "strain_relative_to_reference": self.strain_relative_to_reference,
            "reflections": [
                {
                    "miller_indices": list(indices),
                    "two_theta_deg": float(angle),
                    "residual_two_theta_deg": float(residual),
                    "systematic_shift_deg": float(shift),
                }
                for indices, angle, residual, shift in zip(
                    self.miller_indices,
                    self.two_theta_deg,
                    self.residual_two_theta_deg,
                    self.systematic_shift_deg,
                    strict=True,
                )
            ],
        }

    def describe(self) -> str:
        """Return convention-explicit scientific prose about this determination."""

        method_prose = {
            "cohen": (
                "weighted linear least squares in sin^2(theta), which is linear in the components "
                "of the reciprocal metric tensor G* and therefore has no starting guess, no local "
                "minimum, and an analytic covariance"
            ),
            "average": (
                "an unweighted mean of the lattice parameter computed separately from each "
                "reflection, which removes random scatter as 1/sqrt(N) but cannot remove a "
                "systematic error at all, because that error depends on theta"
            ),
            "le_bail": (
                "whole-pattern decomposition with iteratively extracted reflection intensities, "
                "so that every measured point contributes and neither texture nor the atomic "
                "basis can bias the cell"
            ),
        }[self.method]
        parameters = ", ".join(self.free_parameter_names)
        cell = f"a = {self.a:.6f} +/- {self.a_standard_uncertainty:.6f} angstrom"
        if self.crystal_system not in {"cubic"}:
            cell += f", c = {self.c:.6f} +/- {self.c_standard_uncertainty:.6f} angstrom"
            cell += f", c/a = {self.axial_ratio:.6f}"
        if self.crystal_system in {"orthorhombic", "monoclinic", "triclinic"}:
            cell += f", b = {self.b:.6f} +/- {self.b_standard_uncertainty:.6f} angstrom"
        if self.crystal_system in {"monoclinic", "triclinic"}:
            cell += f", beta = {self.beta_deg:.5f} degrees"

        shifts = self.systematic_shift_deg
        if self.extrapolation == "none":
            drift = (
                " No systematic-error term was refined, so any uncorrected zero-point, specimen "
                "displacement or transparency error is inside the quoted cell rather than "
                "removed from it."
            )
        else:
            largest = float(np.max(np.abs(shifts))) if shifts.size else 0.0
            significance = (
                "significantly different from zero"
                if abs(self.drift_coefficient) > 2.0 * self.drift_standard_uncertainty
                else "not significantly different from zero, so the specimen and instrument were "
                "already well aligned"
            )
            drift = (
                f" A systematic-error term was refined against the {self.extrapolation} "
                f"extrapolation function, giving D = {self.drift_coefficient:.3e} +/- "
                f"{self.drift_standard_uncertainty:.1e}, which is {significance}. It moved the "
                f"lowest-angle reflection by up to {1000.0 * largest:.1f} millidegrees; the "
                f"function vanishes at theta = 90 degrees, so the reported cell is the "
                f"extrapolated one."
            )

        precision = (
            "This reaches the 1e-5 level that thermal-expansion and strain work needs"
            if self.relative_uncertainty < 5.0e-5
            else (
                "This is about 1e-4, useful for phase and composition work but marginal for strain"
                if self.relative_uncertainty < 5.0e-4
                else "This is about 1e-3, which is phase-identification precision and is not "
                "sufficient for strain measurement"
            )
        )
        fit_quality = (
            "the model describes the positions within their stated uncertainties"
            if 0.3 <= self.reduced_chi_squared <= 3.0
            else (
                "the residuals are larger than the position uncertainties allow, so an "
                "aberration remains uncorrected, the phase is not quite right, or the "
                "uncertainties are underestimated"
                if self.reduced_chi_squared > 3.0
                else "the residuals are smaller than the stated uncertainties, which usually "
                "means the position uncertainties are overestimated"
            )
        )
        strain = self.strain_relative_to_reference
        strain_prose = (
            ""
            if strain is None
            else (
                f" Relative to the reference cell the lattice strain along a is "
                f"{strain:+.3e}. This is a strain along the scattering vector of a symmetric "
                f"scan, that is normal to the specimen surface; converting it to a stress needs "
                f"measurements at several specimen tilts and the X-ray elastic constants of the "
                f"reflections used, which this result does not attempt."
            )
        )
        return (
            f"Lattice parameters of '{self.phase_name}' determined by {method_prose}. The "
            f"{self.crystal_system} system leaves {len(self.free_parameter_names)} free "
            f"reciprocal-cell parameters ({parameters}), determined from "
            f"{self.reflection_count} reflections. {cell}. The relative uncertainty is "
            f"{self.relative_uncertainty:.2e}. {precision}.{drift} The reduced chi-squared is "
            f"{self.reduced_chi_squared:.3f}, so {fit_quality}.{strain_prose} "
            f"{_CITATION_CULLITY} {_CITATION_COHEN}"
        )


def _propagate_cell_uncertainty(
    components: np.ndarray,
    covariance: np.ndarray,
    constraint: np.ndarray,
) -> tuple[float, float, float]:
    """Return ``sigma(a), sigma(b), sigma(c)`` by propagating the fit covariance.

    The map from the free reciprocal parameters to the direct cell edges is
    smooth but not linear (it goes through a matrix inversion), so the
    covariance is propagated with a numerically differentiated Jacobian:
    ``sigma^2 = J C J^T``. Finite differences are used rather than an analytic
    derivative because the analytic form differs per crystal system and the
    numerical one demonstrably does not.
    """

    free = covariance.shape[0]
    base = _cell_from_reciprocal_tensor(components)[:3]
    jacobian = np.zeros((3, free), dtype=np.float64)
    parameters = np.linalg.lstsq(constraint, components, rcond=None)[0]
    for index in range(free):
        step = max(abs(parameters[index]) * 1.0e-6, 1.0e-12)
        shifted = parameters.copy()
        shifted[index] += step
        moved = _cell_from_reciprocal_tensor(constraint @ shifted)[:3]
        jacobian[:, index] = (np.asarray(moved) - np.asarray(base)) / step
    variances = np.diag(jacobian @ covariance @ jacobian.T)
    safe = np.sqrt(np.maximum(variances, 0.0))
    return (float(safe[0]), float(safe[1]), float(safe[2]))


def determine_lattice_parameters(
    indexing: PeakIndexing,
    phase: Phase,
    *,
    method: LatticeMethod = "cohen",
    extrapolation: ExtrapolationFunction = "nelson_riley",
    minimum_two_theta_deg: float | None = None,
    name: str | None = None,
) -> LatticeParameterResult:
    """Determine a unit cell from indexed peak positions.

    Purpose
    -------
    Produce the most precise lattice parameter the measured positions support,
    with a standard uncertainty that can be compared against another
    measurement, and with the systematic error removed rather than averaged
    over.

    Method
    ------
    Bragg's law in the form that makes this a *linear* problem:

    ``sin^2(theta) = (lambda^2 / 4) h^T G* h``

    where ``G*`` is the reciprocal metric tensor. The quadratic form expands as
    a dot product with the six independent components of ``G*`` against the row
    ``[h^2, k^2, l^2, 2hk, 2hl, 2kl]``, and the crystal system supplies a
    constraint matrix reducing those six to the one (cubic) to six (triclinic)
    genuinely free parameters. For the hexagonal system, ``a* = b*`` and
    ``gamma* = 60`` degrees give ``G*12 = G*11 / 2``, and the quadratic form
    collapses to the familiar ``A (h^2 + hk + k^2) + C l^2`` -- obtained as a
    consequence of the symmetry rather than written in as a special case.

    A systematic-error column is appended:

    ``sin^2(theta)_obs = (lambda^2 / 4) h^T G* h + D sin^2(theta) f(theta)``

    Observations are weighted by ``1 / sigma^2(sin^2 theta)``, with
    ``sigma(sin^2 theta) = sin(2 theta) sigma(2 theta) / 2`` propagated from
    each peak's own fitted position uncertainty. High-angle reflections
    therefore dominate, which is correct: ``Delta d / d = -cot(theta)
    Delta theta`` makes them intrinsically more informative, and this weighting
    is the quantitative form of that fact.

    The normal equations are solved once. The covariance is
    ``(X^T W X)^-1`` scaled by the reduced chi-squared, and the cell
    uncertainties are propagated through the reciprocal-to-direct inversion
    with a numerically differentiated Jacobian.

    ``method="average"`` instead computes ``a`` separately from every
    reflection and takes the unweighted mean, reporting the standard error of
    that mean. It is offered for comparison and is available for cubic cells
    only, because outside the cubic system one reflection cannot determine two
    cell parameters and a "lattice parameter per reflection" does not exist.

    Parameters
    ----------
    indexing
        Indexed reflections, from
        :func:`~pytex.diffraction.xrd_indexing.index_peaks`.
    phase
        The phase whose symmetry constrains the cell, and whose lattice is
        retained as the reference for a strain.
    method
        ``"cohen"`` (default) or ``"average"``. For ``"le_bail"`` call
        :func:`determine_lattice_parameters_le_bail`, which needs the whole
        profile rather than a peak list.
    extrapolation
        The systematic-error function, one of
        :data:`EXTRAPOLATION_FUNCTIONS`. Ignored by ``"average"``, which by
        construction cannot correct a systematic error.
    minimum_two_theta_deg
        Discard reflections below this angle. Restricting a determination to
        high angles is the crudest systematic-error defence and is worth
        comparing against a refined drift term.
    name
        Ignored placeholder for symmetry with the other constructors; the
        phase name is used.

    Returns
    -------
    LatticeParameterResult
        The determined cell, its uncertainties, the drift term, and the
        per-reflection residuals.

    Raises
    ------
    ValueError
        If the method or extrapolation is unknown, the indexing carries no
        radiation, ``"average"`` is asked for a non-cubic cell, or there are
        fewer reflections than free parameters.

    See Also
    --------
    determine_lattice_parameters_le_bail : the whole-pattern alternative.
    nelson_riley_extrapolation : the classical graphical form of the same idea.

    Notes
    -----
    Read :attr:`LatticeParameterResult.systematic_shift_deg` before believing
    the result. If the drift term moved the peaks by far more than their
    position uncertainties, it did real work and the uncorrected cell would
    have been wrong; if it moved them by less, the specimen was well aligned
    and the two answers should agree.
    """

    if method not in {"cohen", "average"}:
        raise ValueError(
            "determine_lattice_parameters handles 'cohen' and 'average'; for 'le_bail' call "
            "determine_lattice_parameters_le_bail, which needs the measured profile."
        )
    if extrapolation not in EXTRAPOLATION_FUNCTIONS:
        raise ValueError(
            f"determine_lattice_parameters requires extrapolation in {EXTRAPOLATION_FUNCTIONS}."
        )
    if indexing.radiation is None:
        raise ValueError("determine_lattice_parameters needs the indexing to carry a radiation.")

    reflections = list(indexing.reflections)
    if minimum_two_theta_deg is not None:
        reflections = [
            item
            for item in reflections
            if item.peak.two_theta_deg >= float(minimum_two_theta_deg)
        ]
    if not reflections:
        raise ValueError(
            "No indexed reflection survives the angular restriction, so there is nothing to "
            "determine a cell from."
        )

    system = crystal_system_of(phase)
    names, constraint = _constraint_matrix(system)
    wavelength = float(indexing.radiation.wavelength_angstrom)

    angles = np.array([item.peak.two_theta_deg for item in reflections], dtype=np.float64)
    uncertainties = np.array(
        [item.peak.two_theta_standard_uncertainty_deg for item in reflections],
        dtype=np.float64,
    )
    indices = np.array([item.miller_indices for item in reflections], dtype=np.float64)
    miller = tuple(item.miller_indices for item in reflections)

    theta = np.deg2rad(0.5 * angles)
    observed = np.square(np.sin(theta))
    # sigma(sin^2 theta) = |d(sin^2 theta) / d(2 theta)| sigma(2 theta)
    #                    = sin(2 theta) sigma(2 theta) / 2, in radians.
    sigma = 0.5 * np.abs(np.sin(2.0 * theta)) * np.deg2rad(uncertainties)
    sigma = np.maximum(sigma, 1.0e-12)

    if method == "average":
        if len(names) != 1:
            raise ValueError(
                f"The 'average' method is defined only for a cubic cell. The {system} system "
                f"has {len(names)} free cell parameters, and a single reflection cannot "
                "determine more than one, so a lattice parameter per reflection does not exist. "
                "Use method='cohen'."
            )
        # One reflection, one parameter: solve each independently, then average.
        rows = _quadratic_rows(indices) @ constraint
        per_reflection = observed / (0.25 * wavelength**2 * rows[:, 0])
        parameters = np.array([float(np.mean(per_reflection))])
        count = per_reflection.size
        spread = float(np.std(per_reflection, ddof=1)) if count > 1 else 0.0
        standard_error = spread / np.sqrt(count) if count > 1 else float("inf")
        covariance = np.array([[standard_error**2]])
        drift = 0.0
        drift_uncertainty = 0.0
        used_extrapolation: ExtrapolationFunction = "none"
        calculated = 0.25 * wavelength**2 * rows[:, 0] * parameters[0]
        residual_sin_squared = observed - calculated
        reduced_chi_squared = (
            float(np.sum(np.square(residual_sin_squared / sigma)) / max(count - 1, 1))
            if count > 1
            else float("nan")
        )
    else:
        design = 0.25 * wavelength**2 * (_quadratic_rows(indices) @ constraint)
        used_extrapolation = extrapolation
        if extrapolation != "none":
            drift_column = observed * extrapolation_values(angles, function=extrapolation)
            design = np.column_stack([design, drift_column])
        free = design.shape[1]
        if angles.size < free:
            raise ValueError(
                f"A {system} determination with the '{extrapolation}' term needs at least "
                f"{free} indexed reflections and was given {angles.size}."
            )
        weights = 1.0 / sigma
        weighted_design = design * weights[:, None]
        weighted_observed = observed * weights
        solution, *_ = np.linalg.lstsq(weighted_design, weighted_observed, rcond=None)
        residual_sin_squared = observed - design @ solution
        degrees_of_freedom = max(angles.size - free, 1)
        reduced_chi_squared = float(
            np.sum(np.square(residual_sin_squared / sigma)) / degrees_of_freedom
        )
        normal = weighted_design.T @ weighted_design
        try:
            full_covariance = np.linalg.inv(normal) * reduced_chi_squared
        except np.linalg.LinAlgError as error:  # pragma: no cover - singular designs
            raise ValueError(
                "The least-squares design matrix is singular: the indexed reflections do not "
                "determine every free cell parameter. This happens when, for example, no "
                "reflection with a non-zero l index was indexed in a hexagonal pattern."
            ) from error
        if extrapolation != "none":
            parameters = solution[:-1]
            covariance = full_covariance[:-1, :-1]
            drift = float(solution[-1])
            drift_uncertainty = float(np.sqrt(max(full_covariance[-1, -1], 0.0)))
        else:
            parameters = solution
            covariance = full_covariance
            drift = 0.0
            drift_uncertainty = 0.0
        count = int(angles.size)

    components = constraint @ parameters
    cell = _cell_from_reciprocal_tensor(components)
    sigma_a, sigma_b, sigma_c = _propagate_cell_uncertainty(components, covariance, constraint)

    # Residuals converted back to degrees for a reader who thinks in angles.
    residual_deg = np.rad2deg(
        2.0 * residual_sin_squared / np.maximum(np.abs(np.sin(2.0 * theta)), 1.0e-12)
    )

    return LatticeParameterResult(
        method=method,
        phase_name=indexing.phase_name,
        crystal_system=system,
        a=cell[0],
        b=cell[1],
        c=cell[2],
        alpha_deg=cell[3],
        beta_deg=cell[4],
        gamma_deg=cell[5],
        a_standard_uncertainty=sigma_a,
        b_standard_uncertainty=sigma_b,
        c_standard_uncertainty=sigma_c,
        free_parameter_names=names,
        extrapolation=used_extrapolation,
        drift_coefficient=drift,
        drift_standard_uncertainty=drift_uncertainty,
        reflection_count=count,
        reduced_chi_squared=reduced_chi_squared,
        residual_two_theta_deg=residual_deg,
        miller_indices=miller,
        two_theta_deg=angles,
        reference_lattice=phase.lattice,
        settings={
            "wavelength_angstrom": wavelength,
            "source_indexing": indexing.name,
            "minimum_two_theta_deg": float(minimum_two_theta_deg or 0.0),
        },
    )


def nelson_riley_extrapolation(
    indexing: PeakIndexing,
    phase: Phase,
    *,
    function: ExtrapolationFunction = "nelson_riley",
) -> dict[str, np.ndarray | float]:
    """Return the classical ``a`` versus ``f(theta)`` extrapolation, for teaching.

    Purpose
    -------
    Reproduce the plot that made precise parameter measurement possible before
    least squares was routine, and that still explains *why* it works better
    than any amount of prose: the per-reflection values fall on a straight
    line, and where that line meets ``f(theta) = 0`` is the answer.

    Method
    ------
    A lattice parameter is computed from each reflection independently and
    plotted against the extrapolation function. The straight-line fit is
    unweighted, exactly as the graphical construction is, and its intercept is
    the extrapolated parameter. Comparing that intercept with the value from
    :func:`determine_lattice_parameters` shows what weighting and a joint
    solution buy.

    Parameters
    ----------
    indexing
        Indexed reflections.
    phase
        The phase; must be cubic, since a per-reflection lattice parameter
        exists only when the cell has one free parameter.
    function
        The extrapolation function for the abscissa.

    Returns
    -------
    dict
        ``"extrapolation_function"`` and ``"lattice_parameter"`` arrays for the
        plot, plus the fitted ``"intercept"`` and ``"slope"``.

    Raises
    ------
    ValueError
        If the phase is not cubic, or fewer than two reflections are indexed.
    """

    system = crystal_system_of(phase)
    if system != "cubic":
        raise ValueError(
            f"A per-reflection lattice parameter exists only for a cubic cell; the {system} "
            "system needs a joint solution. Use determine_lattice_parameters."
        )
    if indexing.radiation is None:
        raise ValueError("nelson_riley_extrapolation needs the indexing to carry a radiation.")
    if len(indexing.reflections) < 2:
        raise ValueError("An extrapolation needs at least two indexed reflections.")

    wavelength = float(indexing.radiation.wavelength_angstrom)
    angles = np.array(
        [item.peak.two_theta_deg for item in indexing.reflections], dtype=np.float64
    )
    sums = np.array(
        [sum(value**2 for value in item.miller_indices) for item in indexing.reflections],
        dtype=np.float64,
    )
    spacings = wavelength / (2.0 * np.sin(np.deg2rad(0.5 * angles)))
    per_reflection = spacings * np.sqrt(sums)
    abscissa = extrapolation_values(angles, function=function)
    slope, intercept = np.polyfit(abscissa, per_reflection, 1)
    return {
        "extrapolation_function": abscissa,
        "lattice_parameter": per_reflection,
        "two_theta_deg": angles,
        "intercept": float(intercept),
        "slope": float(slope),
    }


def determine_lattice_parameters_le_bail(
    measured: MeasuredPowderPattern,
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    instrument: InstrumentBroadening | None = None,
    systematic: SystematicTerm = "zero",
    goniometer_radius_mm: float = 240.0,
    cycles: int = 12,
    max_index: int = 6,
    background_half_window_deg: float = 2.0,
    two_theta_range_deg: tuple[float, float] | None = None,
) -> LatticeParameterResult:
    """Determine a cell by whole-pattern decomposition with extracted intensities.

    Purpose
    -------
    Use every measured point rather than a handful of fitted positions, and
    handle overlapped reflections properly. This is the method for a hexagonal
    or lower-symmetry pattern, where single-peak fitting runs out of resolvable
    peaks long before the reflection list runs out of reflections.

    Method
    ------
    Le Bail, Duroy & Fourquet (1988). The calculated profile is a sum of
    pseudo-Voigt peaks at the positions the current cell puts them, each scaled
    by an intensity that is *not* a refined parameter:

    1. Every reflection starts with equal intensity.
    2. **Extraction.** The observed intensity at each point is partitioned
       between the reflections that overlap it, in proportion to their current
       contributions:

       ``I_k <- sum_i y_obs(i) I_k P_k(i) / sum_j I_j P_j(i)``

       This is the whole trick. A Pawley fit instead treats the intensities as
       free least-squares parameters, whose normal matrix becomes singular
       exactly when two reflections overlap completely -- which is the case the
       method exists to handle. Le Bail's partition is stable there, because
       two exactly coincident reflections simply split the intensity in their
       current ratio and neither the cell nor the fit notices.
    3. **Refinement.** With the intensities held at their extracted values, the
       cell parameters, the systematic term, the Caglioti width coefficients
       and the mixing parameter are refined by bounded Levenberg-Marquardt.
    4. Repeat.

    Each reflection's profile is the whole K-alpha multiplet, not a single
    line: when the radiation declares a K-alpha2 wavelength, a partner peak is
    added at the position Bragg's law puts it at, sharing the width and mixing
    and scaled by the tabulated intensity ratio. It costs no parameter.
    Omitting it is not a cosmetic error -- the unmodelled alpha2 peak is a
    residual as large as itself, and the refinement absorbs it into the cell.

    The background's *shape* is removed first with SNIP rather than refined
    simultaneously, so that a flexible background cannot absorb peak intensity
    and shift the cell with it. What SNIP leaves is a small level offset, and
    that is refined -- as a straight line in the reduced angular coordinate,
    two parameters, deliberately far too stiff to follow a peak. Refusing to
    model it is not the conservative choice: on a clean synthetic pattern the
    unmodelled offset accounted for about ninety per cent of the total misfit,
    all of it between the peaks.

    Exactly one systematic term is refined -- a detector zero, *or* a specimen
    displacement, *or* neither. Refining both from a single specimen scan is
    ill-conditioned: they differ only as constant against ``cos(theta)``, and
    over one pattern's angular range that difference is comparable to the
    noise. Zero belongs to a calibrated instrument; displacement belongs to the
    specimen.

    Parameters
    ----------
    measured
        The measured profile, background included.
    phase
        The phase whose symmetry constrains the cell and whose lattice is the
        starting point.
    radiation
        Falls back to ``measured.radiation``.
    instrument
        Starting Caglioti coefficients. Defaults to a laboratory
        Bragg-Brentano resolution function.
    systematic
        ``"zero"``, ``"displacement"`` or ``"none"``.
    goniometer_radius_mm
        Needed to interpret a refined displacement in millimetres.
    cycles
        Number of extract-then-refine cycles. Convergence is normally reached
        in five to ten.
    max_index
        Largest index enumerated.
    background_half_window_deg
        SNIP half-window.
    two_theta_range_deg
        Restrict the fit to this window.

    Returns
    -------
    LatticeParameterResult
        The determined cell with uncertainties from the final Jacobian, and
        residuals evaluated at the modelled reflection positions.

    Raises
    ------
    ValueError
        If the systematic term is unknown, no radiation is available, the
        phase predicts no reflection in range, or the number of cycles is not
        positive.

    See Also
    --------
    determine_lattice_parameters : the peak-position method, which should agree
        with this one to within their combined uncertainties.

    Notes
    -----
    Le Bail intensities are *extracted*, not measured: for two reflections that
    overlap completely, the partition between them is whatever ratio the
    iteration started with. They are therefore fine for describing the profile
    and unfit for structural work, which is precisely why this method cannot be
    biased by texture -- it never claims to know what the intensities should
    have been.
    """

    if systematic not in SYSTEMATIC_TERMS:
        raise ValueError(
            f"determine_lattice_parameters_le_bail requires systematic in {SYSTEMATIC_TERMS}."
        )
    if cycles < 1:
        raise ValueError("determine_lattice_parameters_le_bail requires at least one cycle.")
    spec = radiation if radiation is not None else measured.radiation
    if spec is None:
        raise ValueError(
            "determine_lattice_parameters_le_bail needs a radiation: the pattern declares none."
        )

    system = crystal_system_of(phase)
    names, constraint = _constraint_matrix(system)
    wavelength = float(spec.wavelength_angstrom)

    background = estimate_background(
        measured, method="snip", half_window_deg=background_half_window_deg
    )
    axis = np.asarray(measured.two_theta_deg, dtype=np.float64)
    raw = np.asarray(measured.intensity, dtype=np.float64)
    stated = (
        None
        if measured.standard_uncertainty is None
        else np.asarray(measured.standard_uncertainty, dtype=np.float64)
    )
    # Not clipped at zero: clipping a symmetric noise excursion is a one-sided
    # operation and biases the background level upward, which the peak model
    # then cannot represent.
    observed = raw - background.background
    if two_theta_range_deg is not None:
        window = (axis >= float(two_theta_range_deg[0])) & (
            axis <= float(two_theta_range_deg[1])
        )
        axis = axis[window]
        raw = raw[window]
        observed = observed[window]
        stated = None if stated is None else stated[window]
    if axis.size < 32:
        raise ValueError("A whole-pattern fit needs at least thirty-two measured points.")

    # Enumerate once, over a padded range so a shifting cell cannot move a
    # reflection into or out of the model mid-refinement and change the
    # function being minimized.
    reflections = generate_powder_reflections(
        phase,
        radiation=spec,
        two_theta_range_deg=(
            max(float(axis[0]) - 2.0, 1.0e-3),
            min(float(axis[-1]) + 2.0, 179.999),
        ),
        max_index=max_index,
    )
    if not reflections:
        raise ValueError(
            "The phase predicts no reflection inside the measured range, so there is no "
            "whole-pattern model to fit."
        )
    indices = np.array(
        [[int(value) for value in item.miller_indices] for item in reflections],
        dtype=np.float64,
    )
    miller = tuple(
        (int(item.miller_indices[0]), int(item.miller_indices[1]), int(item.miller_indices[2]))
        for item in reflections
    )
    rows = _quadratic_rows(indices) @ constraint

    resolution = (
        instrument if instrument is not None else InstrumentBroadening.laboratory_bragg_brentano()
    )
    start_parameters = np.linalg.lstsq(
        constraint, phase.lattice.reciprocal_metric_tensor()[
            np.array([0, 1, 2, 0, 0, 1]), np.array([0, 1, 2, 1, 2, 2])
        ],
        rcond=None,
    )[0]

    free_cell = start_parameters.size
    # [cell..., systematic, U, V, W, eta, residual background constant and slope]
    scale = float(np.max(np.abs(observed))) or 1.0
    start = np.concatenate(
        [
            start_parameters,
            np.array([0.0]),
            np.array(
                [resolution.caglioti_u, resolution.caglioti_v, resolution.caglioti_w, 0.5]
            ),
            np.array([0.0, 0.0]),
        ]
    )
    lower = np.concatenate(
        [
            start_parameters * 0.9,
            np.array([-1.0]),
            np.array([-0.5, -0.5, 1.0e-5, 0.0]),
            np.array([-0.05 * scale, -0.05 * scale]),
        ]
    )
    upper = np.concatenate(
        [
            start_parameters * 1.1,
            np.array([1.0]),
            np.array([0.5, 0.5, 1.0, 1.0]),
            np.array([0.05 * scale, 0.05 * scale]),
        ]
    )
    if systematic == "none":
        lower[free_cell] = -1.0e-12
        upper[free_cell] = 1.0e-12

    # A straight line in the reduced angular coordinate, deliberately too stiff
    # to follow a peak. SNIP removes the background's shape; what is left is a
    # small level offset, and refusing to model it puts that offset into the
    # residual of every point between the peaks -- which on a clean synthetic
    # pattern is where nearly all of the misfit was found to live.
    reduced_axis = (axis - axis[0]) / max(axis[-1] - axis[0], 1.0e-9) * 2.0 - 1.0

    def residual_background(parameters: np.ndarray) -> np.ndarray:
        level: np.ndarray = parameters[-2] + parameters[-1] * reduced_axis
        return level

    # Weights come from the *measured* counts, never from the
    # background-subtracted profile. Subtraction removes signal, not variance:
    # a background point with 150 counts still has a standard deviation near
    # 12, and weighting it as though its uncertainty were one -- which is what
    # using the subtracted value does -- over-weights every point between the
    # peaks by an order of magnitude and hands the fit to the background.
    weights = (
        1.0 / stated if stated is not None else 1.0 / np.sqrt(np.maximum(raw, 1.0))
    )

    def positions_and_shapes(
        parameters: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, float]:
        cell_parameters = parameters[:free_cell]
        sin_squared = 0.25 * wavelength**2 * (rows @ cell_parameters)
        visible = (sin_squared > 1.0e-9) & (sin_squared < 1.0)
        angles = np.full(sin_squared.shape, np.nan)
        angles[visible] = np.rad2deg(2.0 * np.arcsin(np.sqrt(sin_squared[visible])))
        offset = float(parameters[free_cell])
        if systematic == "displacement":
            theta = np.deg2rad(0.5 * angles)
            angles = angles + np.rad2deg(-2.0 * offset * np.cos(theta) / goniometer_radius_mm)
        elif systematic == "zero":
            angles = angles + offset
        widths_squared = (
            parameters[free_cell + 1] * np.square(np.tan(np.deg2rad(0.5 * angles)))
            + parameters[free_cell + 2] * np.tan(np.deg2rad(0.5 * angles))
            + parameters[free_cell + 3]
        )
        widths = np.sqrt(np.maximum(widths_squared, 1.0e-6))
        return angles, widths, float(parameters[free_cell + 4])

    # The measured profile of one reflection is the whole K-alpha multiplet, not
    # a single line. Modelling only alpha1 against doublet data leaves a
    # residual as large as the alpha2 peak itself, which the refinement then
    # tries to absorb into the cell -- a systematic lattice-parameter error, not
    # merely a poor fit. The partner is not a free parameter: Bragg's law at
    # fixed d puts it at sin(theta_2) = (lambda_2 / lambda_1) sin(theta_1).
    doublet: tuple[float, float] | None = None
    if spec.kalpha2_wavelength_angstrom is not None:
        doublet = (
            float(spec.kalpha2_wavelength_angstrom / spec.wavelength_angstrom),
            float(spec.kalpha2_relative_intensity),
        )

    def _accumulate(
        stack: np.ndarray, index: int, centre: float, width: float, eta: float, scale: float
    ) -> None:
        near = np.abs(axis - centre) <= 8.0 * width
        if not np.any(near):
            return
        stack[index, near] += scale * pseudo_voigt_profile(
            axis[near], centre_deg=centre, fwhm_deg=width, eta=eta
        )

    def profiles(parameters: np.ndarray) -> np.ndarray:
        angles, widths, eta = positions_and_shapes(parameters)
        stack = np.zeros((angles.size, axis.size), dtype=np.float64)
        for index, (centre, width) in enumerate(zip(angles, widths, strict=True)):
            if not np.isfinite(centre):
                continue
            _accumulate(stack, index, float(centre), float(width), eta, 1.0)
            if doublet is not None:
                ratio, relative = doublet
                argument = ratio * np.sin(np.deg2rad(0.5 * float(centre)))
                if abs(argument) <= 1.0:
                    partner = float(np.rad2deg(2.0 * np.arcsin(argument)))
                    _accumulate(stack, index, partner, float(width), eta, relative)
        # Normalize each reflection's profile to unit sum over the measured
        # points. This makes the extracted quantity the *integrated* intensity,
        # which is what the Le Bail partition actually computes: the partition
        # sums observed counts, so using its output as the amplitude of a
        # unit-height profile overstates every peak by the number of points
        # under it -- about a factor of twenty at this step size.
        totals = stack.sum(axis=1)
        np.divide(
            stack,
            np.where(totals > 0.0, totals, 1.0)[:, None],
            out=stack,
            where=(totals > 0.0)[:, None],
        )
        return stack

    intensities = np.full(
        len(reflections), float(np.sum(np.maximum(observed, 0.0))) / len(reflections)
    )
    parameters = start.copy()
    solution = None
    for _ in range(int(cycles)):
        # --- Extraction: partition the observed profile between reflections.
        stack = profiles(parameters)
        peak_only = np.maximum(observed - residual_background(parameters), 0.0)
        contributions = intensities[:, None] * stack
        total = contributions.sum(axis=0)
        share = np.divide(
            contributions,
            np.where(total > 0.0, total, 1.0),
            out=np.zeros_like(contributions),
            where=total > 0.0,
        )
        extracted = share @ peak_only
        # A reflection that the current cell places outside the data receives
        # nothing; keep it alive at a small value so it can return if the cell
        # moves, rather than dropping permanently out of the model.
        floor = 1.0e-6 * float(np.max(extracted)) if np.any(extracted > 0.0) else 1.0
        intensities = np.maximum(extracted, floor)

        # --- Refinement: cell and profile, with the intensities held fixed.
        def residual(values: np.ndarray, fixed: np.ndarray = intensities) -> np.ndarray:
            calculated = fixed @ profiles(values) + residual_background(values)
            weighted: np.ndarray = (calculated - observed) * weights
            return weighted

        solution = least_squares(
            residual,
            x0=parameters,
            bounds=(lower, upper),
            method="trf",
            # The cell parameters are of order 0.08 inverse square angstrom and
            # the profile coefficients of order 0.01 degrees squared; without
            # Jacobian scaling the trust region is set by whichever happens to
            # be largest and the cell barely moves.
            x_scale="jac",
            max_nfev=600,
        )
        parameters = solution.x

    if solution is None:  # pragma: no cover - cycles >= 1 is validated above
        raise ValueError("The Le Bail refinement performed no cycle.")

    degrees_of_freedom = max(axis.size - parameters.size, 1)
    reduced_chi_squared = float(np.sum(np.square(solution.fun)) / degrees_of_freedom)
    jacobian = np.asarray(solution.jac, dtype=np.float64)
    try:
        full_covariance = np.linalg.inv(jacobian.T @ jacobian) * reduced_chi_squared
    except np.linalg.LinAlgError:  # pragma: no cover - degenerate patterns only
        full_covariance = np.full((parameters.size, parameters.size), np.nan)
    covariance = np.nan_to_num(full_covariance[:free_cell, :free_cell], nan=0.0)
    systematic_uncertainty = float(
        np.sqrt(max(float(np.nan_to_num(full_covariance[free_cell, free_cell])), 0.0))
    )

    cell_parameters = parameters[:free_cell]
    components = constraint @ cell_parameters
    cell = _cell_from_reciprocal_tensor(components)
    sigma_a, sigma_b, sigma_c = _propagate_cell_uncertainty(
        components, covariance, constraint
    )

    angles, _, _ = positions_and_shapes(parameters)
    modelled = np.isfinite(angles) & (angles >= axis[0]) & (angles <= axis[-1])
    strong = intensities >= 0.01 * float(np.max(intensities))
    keep = modelled & strong
    return LatticeParameterResult(
        method="le_bail",
        phase_name=str(getattr(phase, "name", None) or "candidate phase"),
        crystal_system=system,
        a=cell[0],
        b=cell[1],
        c=cell[2],
        alpha_deg=cell[3],
        beta_deg=cell[4],
        gamma_deg=cell[5],
        a_standard_uncertainty=sigma_a,
        b_standard_uncertainty=sigma_b,
        c_standard_uncertainty=sigma_c,
        free_parameter_names=names,
        extrapolation="none",
        drift_coefficient=float(parameters[free_cell]),
        drift_standard_uncertainty=systematic_uncertainty,
        reflection_count=int(np.count_nonzero(keep)),
        reduced_chi_squared=reduced_chi_squared,
        residual_two_theta_deg=np.zeros(int(np.count_nonzero(keep))),
        miller_indices=tuple(item for item, flag in zip(miller, keep, strict=True) if flag),
        two_theta_deg=angles[keep],
        reference_lattice=phase.lattice,
        settings={
            "wavelength_angstrom": wavelength,
            "systematic": systematic,
            "cycles": float(cycles),
            "point_count": float(axis.size),
            "goniometer_radius_mm": float(goniometer_radius_mm),
        },
    )


def determine_lattice_parameters_from_pattern(
    measured: MeasuredPowderPattern,
    phase: Phase,
    *,
    method: LatticeMethod = "cohen",
    extrapolation: ExtrapolationFunction = "nelson_riley",
    instrument: InstrumentBroadening | None = None,
    tolerance_deg: float = 0.3,
    max_index: int = 6,
    prominence_sigma: float = 5.0,
    minimum_two_theta_deg: float | None = None,
    phase_name: str | None = None,
) -> tuple[LatticeParameterResult, PeakIndexing | None]:
    """Run the whole pipeline: detect, fit, index, determine.

    Purpose
    -------
    Provide the single call an operator makes, without hiding the stages: the
    indexing is returned alongside the result so the reflection assignment and
    its figures of merit can be inspected, which is where a wrong answer is
    visible.

    Parameters
    ----------
    measured
        The measured pattern.
    phase
        The candidate phase.
    method
        ``"cohen"``, ``"average"`` or ``"le_bail"``.
    extrapolation
        Systematic-error function for the position methods.
    instrument
        Calibrated resolution function, used for detection and fitting.
    tolerance_deg
        Indexing tolerance.
    max_index
        Largest index enumerated.
    prominence_sigma
        Detection threshold in noise standard deviations.
    minimum_two_theta_deg
        Discard reflections below this angle.
    phase_name
        Name for the report.

    Returns
    -------
    tuple[LatticeParameterResult, PeakIndexing | None]
        The determination, and the indexing behind it. The indexing is
        ``None`` for ``"le_bail"``, which never forms a peak list.

    Raises
    ------
    ValueError
        Propagated from any stage; the message names the stage that failed.
    """

    if method == "le_bail":
        return (
            determine_lattice_parameters_le_bail(
                measured, phase, instrument=instrument, max_index=max_index
            ),
            None,
        )
    table: PeakTable = detect_and_fit_peaks(
        measured, instrument=instrument, prominence_sigma=prominence_sigma
    )
    indexing = index_peaks(
        table,
        phase,
        tolerance_deg=tolerance_deg,
        max_index=max_index,
        phase_name=phase_name,
    )
    result = determine_lattice_parameters(
        indexing,
        phase,
        method=method,
        extrapolation=extrapolation,
        minimum_two_theta_deg=minimum_two_theta_deg,
    )
    return (result, indexing)


__all__ = [
    "EXTRAPOLATION_FUNCTIONS",
    "LATTICE_METHODS",
    "LATTICE_PARAMETER_SCHEMA",
    "SYSTEMATIC_TERMS",
    "ExtrapolationFunction",
    "LatticeMethod",
    "LatticeParameterResult",
    "SystematicTerm",
    "crystal_system_of",
    "determine_lattice_parameters",
    "determine_lattice_parameters_from_pattern",
    "determine_lattice_parameters_le_bail",
    "extrapolation_values",
    "nelson_riley_extrapolation",
]

