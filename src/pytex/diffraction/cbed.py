"""Convergent-beam electron diffraction: discs, rocking curves, and thickness.

Selected-area diffraction illuminates the specimen with a *parallel* beam, so a
reflection is a spot and the pattern carries the projected symmetry of the zone
and nothing else. Converge the beam instead — focus it to a probe with a
semi-angle of a few milliradians — and every spot spreads into a **disc**, with
each point of the disc corresponding to a different incident direction. The disc
is therefore a *map of the rocking curve*: one CBED exposure contains the
intensity as a function of deviation from the Bragg condition, which a parallel
beam could only obtain by tilting the specimen through a series of exposures.

Three things follow, and this module implements all three.

**Thickness.** The two-beam rocking curve has zeros at known deviations, and
their spacing depends on the foil thickness. Measuring the fringe positions in a
single disc yields the local thickness — and, in the same fit, the extinction
distance, so the answer does not depend on a tabulated constant. This is the
standard Kelly *et al.* analysis, implemented as
:func:`thickness_from_fringe_minima`.

**Lattice parameter along the beam.** The convergence also makes higher-order
Laue zones visible as rings, whose radii depend on the reciprocal-lattice layer
spacing along the zone axis — the one lattice dimension a zone-axis SAED pattern
cannot see at all (:func:`holz_ring_radii_inv_angstrom`) — and as sharp *lines*
inside the bright-field disc, which measure the lattice parameters themselves to
a part in ten thousand (`pytex.diffraction.holz`).

**Point and space group.** The symmetry *within* the discs and of the whole
pattern determines the diffraction group, and hence the point group including
the presence or absence of a centre of symmetry, which Friedel's law hides from
kinematic SAED. :meth:`CBEDPattern.determine_point_group` performs that
determination from a simulated pattern, and
`pytex.diffraction.diffraction_groups` carries the group theory behind it.

Two simulation methods, and why the choice matters
--------------------------------------------------

:class:`ConvergentBeamConfig` selects between them.

``"two-beam"`` (the default) computes each disc as an independent two-beam
rocking curve. It is cheap, it is exactly the model the thickness analysis
inverts, and it is the right choice for teaching disc geometry, the
Kossel-Moellenstedt regime and fringe counting. It is the wrong choice for
anything about *relative* intensities or symmetry, because the discs of one
pattern never exchange intensity: each is symmetric in :math:`s` by
construction, so the pattern reports symmetry the crystal may not have.

``"bloch"`` solves the coupled many-beam problem through
`pytex.diffraction.dynamical`, optionally with absorption and with higher-order
Laue zone reflections in the beam set. It is the only method whose symmetry
means anything, and :meth:`CBEDPattern.symmetry_observations` refuses to run on
a two-beam pattern for that reason.

What this module models, and what it does not
---------------------------------------------

Modelled: disc geometry from the convergence semi-angle, the excitation error
across each disc, either the two-beam rocking curve or the full many-beam
dynamical solution with an absorptive potential, the Kossel-Moellenstedt versus
Kossel overlap regime, HOLZ ring radii and HOLZ line geometry, thickness
extraction from fringe minima, and the diffraction-group symmetry determination.

Not modelled: inelastic background, probe aberrations, specimen bending, wedge
or strain gradients, and Bethe perturbation of weak beams (so a beam set
containing a whole HOLZ ring is expensive). The absorptive potential is
phenomenological rather than computed from absorptive form factors. Those limits
are stated on :class:`CBEDPattern` and repeated by its ``describe()``.

The geometry
------------

Let the zone axis point toward the gun, so the beam propagates along
:math:`-\\hat{z}`, and let the incident direction be tilted by
:math:`(\\theta_x, \\theta_y)` within the convergence cone
:math:`|\\boldsymbol{\\theta}| \\le \\alpha`. For a reflection
:math:`\\mathbf{g} = (g_u, g_v, g_z)` in the zone basis, the excitation error is

.. math::

   s_g(\\boldsymbol{\\theta}) = g_z
       - \\theta_x g_u - \\theta_y g_v - \\tfrac{1}{2}\\lambda |\\mathbf{g}|^{2},

which reduces at zero tilt to the parallel-beam expression
:math:`s_g = g_z - \\lambda g^{2}/2` used by
`pytex.diffraction.kinematic`. Two consequences are worth reading off directly:

- :math:`s_g` varies **linearly across the disc, along the direction of**
  :math:`\\mathbf{g}_\\perp`. The fringes are therefore straight and
  perpendicular to :math:`\\mathbf{g}` — which is what a Kossel-Moellenstedt
  pattern looks like.
- A tilt of :math:`\\theta` displaces the diffracted beam by
  :math:`\\theta/\\lambda` in reciprocal space, so the disc radius is
  :math:`\\alpha/\\lambda`, and on the detector
  :math:`(L\\lambda)\\,\\alpha/\\lambda = L\\alpha`. Discs of neighbouring
  reflections touch when :math:`2\\alpha/\\lambda = |\\mathbf{g}|_{\\min}`, which
  separates the two classical regimes.

See Also
--------
`pytex.diffraction.kinematic` : the parallel-beam zone-axis engine this shares
    its conventions with.
`docs/tex/algorithms/convergent_beam_electron_diffraction.tex` : the derivations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from numpy.typing import ArrayLike
from scipy import ndimage

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.diffraction_groups import (
    DiffractionGroup,
    PointGroupDetermination,
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for_zone_axis,
    plane_point_group_symbol,
)
from pytex.diffraction.dynamical import (
    AbsorptionModel,
    BeamSet,
    beam_set_for_zone,
    beam_set_from_indices,
    solve_bloch_waves,
)
from pytex.diffraction.holz import HOLZLinePattern, holz_line_pattern
from pytex.diffraction.kinematic import (
    centering_allowed_mask,
    electron_wavelength_angstrom,
    zone_basis_from_axis,
)
from pytex.diffraction.physics import ReflectionCondition
from pytex.diffraction.scattering import electron_structure_factor_angstrom

#: How a CBED pattern's disc intensities are computed.
SimulationMethod = Literal["two-beam", "bloch"]

__all__ = [
    "CBED_PATTERN_SCHEMA",
    "CBED_THICKNESS_SCHEMA",
    "CBEDDisc",
    "CBEDPattern",
    "ConvergentBeamConfig",
    "SimulationMethod",
    "TwoBeamThicknessReport",
    "electron_structure_factor_angstrom",
    "extinction_distance_angstrom",
    "fringe_minimum_excitation_errors",
    "holz_ring_radii_inv_angstrom",
    "simulate_cbed_pattern",
    "thickness_from_fringe_minima",
    "two_beam_rocking_curve",
]

#: Schema identifier of the CBED-pattern payload.
CBED_PATTERN_SCHEMA = "pytex.cbed_pattern/1"

#: Schema identifier of the two-beam thickness payload.
CBED_THICKNESS_SCHEMA = "pytex.cbed_two_beam_thickness/1"

_DEBYE_WALLER_DENOMINATOR = 16.0 * np.pi * np.pi


# --------------------------------------------------------------------------- #
# Extinction distances
# --------------------------------------------------------------------------- #


def extinction_distance_angstrom(
    phase: Phase,
    hkl: ArrayLike,
    *,
    beam_energy_kev: float = 200.0,
) -> np.ndarray:
    """Two-beam extinction distances in angstrom.

    What it does
        Returns

        .. math::

           \\xi_{g} = \\frac{\\pi V_{c}\\cos\\theta_{B}}{\\lambda\\,|F_{g}|},

        the length over which the intensity oscillates completely between the
        transmitted and the diffracted beam under exact Bragg conditions. It is
        the single parameter that sets the scale of every dynamical effect in a
        two-beam calculation.

    When to use it
        To predict thickness-fringe spacing, to judge whether a foil is thin
        enough for the kinematic approximation (it is not, once
        :math:`t \\gtrsim \\xi_{g}/3`), and as the starting value for a CBED
        thickness fit — though the fit itself returns a *measured*
        :math:`\\xi_{g}`, which is the number to prefer.

    Parameters
    ----------
    phase:
        Must carry a unit cell.
    hkl:
        ``(n, 3)`` Miller indices, or one triple. A forbidden reflection has
        :math:`F_{g} = 0` and an infinite extinction distance, which is returned
        as ``inf`` rather than raising: it is the correct answer.
    beam_energy_kev:
        Accelerating voltage.

    Returns
    -------
    np.ndarray
        ``(n,)`` extinction distances in angstrom.

    Notes
    -----
    Accuracy is limited by the fitted scattering-factor parametrization, which
    is good to a few percent for light elements and degrades for heavy ones. As
    a calibration point, aluminium ``{111}`` at 100 kV returns 555 angstrom
    against the tabulated 556 (Williams and Carter, *Transmission Electron
    Microscopy*, 2nd ed., Table 23.1). Do not quote these to three figures for a
    heavy element; measure them instead, with
    :func:`thickness_from_fringe_minima`.
    """

    indices = np.atleast_2d(np.asarray(hkl, dtype=np.int64))
    structure_factors = electron_structure_factor_angstrom(
        phase, indices, beam_energy_kev=beam_energy_kev
    )
    wavelength = electron_wavelength_angstrom(beam_energy_kev)
    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_magnitude = np.linalg.norm(indices.astype(np.float64) @ reciprocal.T, axis=1)
    sin_theta = np.clip(wavelength * g_magnitude / 2.0, -1.0, 1.0)
    cos_theta = np.sqrt(np.clip(1.0 - sin_theta * sin_theta, 0.0, 1.0))

    volume = abs(float(np.linalg.det(phase.lattice.direct_basis().matrix)))
    magnitude = np.abs(structure_factors)
    with np.errstate(divide="ignore"):
        distances = np.where(
            magnitude > 1e-12,
            np.pi * volume * cos_theta / (wavelength * np.maximum(magnitude, 1e-300)),
            np.inf,
        )
    return np.asarray(distances, dtype=np.float64)


# --------------------------------------------------------------------------- #
# The two-beam rocking curve
# --------------------------------------------------------------------------- #


def two_beam_rocking_curve(
    excitation_error_inv_angstrom: ArrayLike,
    *,
    thickness_angstrom: float,
    extinction_distance_angstrom: float,
) -> np.ndarray:
    """Diffracted intensity of the two-beam dynamical solution.

    What it does
        Evaluates

        .. math::

           I_{g}(s) = \\frac{\\sin^{2}(\\pi t s_{\\mathrm{eff}})}
                            {(\\xi_{g}s_{\\mathrm{eff}})^{2}},
           \\qquad
           s_{\\mathrm{eff}} = \\sqrt{s^{2} + \\xi_{g}^{-2}},

        the Howie-Whelan two-beam result without absorption. The transmitted
        intensity is :math:`1 - I_{g}`, exactly: the two beams exchange
        intensity and nothing is lost.

    When to use it
        To model a CBED disc, to predict thickness-fringe spacing in a wedge,
        and to see where the kinematic approximation stops being true — it is
        the :math:`\\xi_{g} \\to \\infty` limit of this expression, in which
        :math:`s_{\\mathrm{eff}} \\to s` and the prefactor becomes
        :math:`(\\pi t/\\xi_{g})^{2}`.

    Parameters
    ----------
    excitation_error_inv_angstrom:
        Deviation parameter ``s``, any shape, in 1/angstrom.
    thickness_angstrom:
        Foil thickness ``t``, strictly positive.
    extinction_distance_angstrom:
        ``xi_g``, strictly positive. ``inf`` is accepted and yields zero
        intensity, the correct answer for a forbidden reflection.

    Returns
    -------
    np.ndarray
        Intensity in ``[0, 1]``, the shape of the input.

    Notes
    -----
    At exact Bragg incidence, :math:`s = 0` and
    :math:`I_{g} = \\sin^{2}(\\pi t/\\xi_{g})`: complete oscillation between the
    beams with thickness, which is what makes :math:`\\xi_{g}` measurable. The
    minima away from :math:`s = 0` occur where
    :math:`t\\,s_{\\mathrm{eff}} = n` for integer ``n``, and those are the
    positions :func:`thickness_from_fringe_minima` reads.
    """

    if not np.isfinite(thickness_angstrom) or thickness_angstrom <= 0.0:
        raise ValueError("thickness_angstrom must be finite and strictly positive.")
    if extinction_distance_angstrom <= 0.0 or math.isnan(extinction_distance_angstrom):
        raise ValueError("extinction_distance_angstrom must be strictly positive.")

    s_values = np.asarray(excitation_error_inv_angstrom, dtype=np.float64)
    if math.isinf(extinction_distance_angstrom):
        return np.zeros_like(s_values)
    inverse_xi = 1.0 / extinction_distance_angstrom
    s_effective = np.sqrt(s_values * s_values + inverse_xi * inverse_xi)
    argument = np.pi * thickness_angstrom * s_effective
    return np.asarray(
        np.square(np.sin(argument)) / np.square(extinction_distance_angstrom * s_effective),
        dtype=np.float64,
    )


def fringe_minimum_excitation_errors(
    excitation_error_inv_angstrom: ArrayLike,
    intensity: ArrayLike,
) -> np.ndarray:
    """Excitation errors of the interior local minima of a rocking curve.

    What it does
        Returns the ``s`` values at which a sampled intensity profile is
        strictly lower than both of its neighbours — the dark fringes an
        experimenter measures off a CBED disc.

    When to use it
        As the bridge between a disc (simulated or measured) and
        :func:`thickness_from_fringe_minima`, which needs the minimum positions
        and nothing else.

    Parameters
    ----------
    excitation_error_inv_angstrom, intensity:
        Matching one-dimensional arrays, with ``s`` monotonically increasing.

    Returns
    -------
    np.ndarray
        The ``s`` values of the interior minima, in input order.

    Notes
    -----
    The result is only as accurate as the sampling: a minimum lying between two
    samples is reported at the nearer sample. That quantization propagates
    directly into the fitted thickness, which is why the sampling density is a
    parameter of :class:`ConvergentBeamConfig` rather than a hidden constant.
    """

    s_values = as_float_array(excitation_error_inv_angstrom, shape=(None,))
    values = as_float_array(intensity, shape=(None,))
    if s_values.shape != values.shape:
        raise ValueError("excitation errors and intensity must have the same shape.")
    if values.size < 3:
        return np.zeros(0, dtype=np.float64)
    interior = values[1:-1]
    is_minimum = (interior < values[:-2]) & (interior < values[2:])
    return np.asarray(s_values[1:-1][is_minimum], dtype=np.float64)


@dataclass(frozen=True, slots=True)
class TwoBeamThicknessReport:
    """Foil thickness and extinction distance measured from CBED fringes.

    Purpose
    -------
    The output of the standard two-beam CBED thickness determination. Its value
    is that it returns **both** unknowns from the same data: the fringe
    positions over-determine the pair :math:`(t, \\xi_{g})`, so the thickness
    does not inherit the error of a tabulated extinction distance.

    Attributes
    ----------
    thickness_angstrom : float
        The fitted foil thickness ``t``.
    extinction_distance_angstrom : float
        The fitted ``xi_g``, to be compared with the tabulated value as a check
        on the whole analysis.
    first_order : int
        The integer ``n`` assigned to the innermost measured minimum. Getting
        this wrong is the classical failure of the method, and is why the
        assignment is reported rather than assumed.
    orders : np.ndarray
        The full integer sequence used.
    excitation_errors_inv_angstrom : np.ndarray
        The measured minimum positions, in input order.
    r_squared : float
        Coefficient of determination of the straight-line fit. A wrong order
        assignment shows as curvature and therefore as a fall in this number.
    residuals : np.ndarray
        Fit residuals in the transformed coordinates.
    provenance : ProvenanceRecord or None
    """

    thickness_angstrom: float
    extinction_distance_angstrom: float
    first_order: int
    orders: np.ndarray
    excitation_errors_inv_angstrom: np.ndarray
    r_squared: float
    residuals: np.ndarray
    provenance: ProvenanceRecord | None = None

    def describe(self) -> str:
        """Convention-explicit prose: the measurement and what could invalidate it."""

        return (
            f"Two-beam CBED analysis of {len(self.orders)} fringe minima, assigned orders "
            f"n = {self.first_order}..{int(self.orders[-1])}: foil thickness "
            f"{self.thickness_angstrom:.1f} A and extinction distance "
            f"{self.extinction_distance_angstrom:.1f} A, from a straight-line fit of "
            f"(s_n/n)^2 against 1/n^2 with R^2 = {self.r_squared:.6f}. Both quantities "
            "come from the same fringe positions, so the thickness does not depend on a "
            "tabulated extinction distance. The assignment of the innermost minimum is "
            "the assumption to check: choosing n too small curves the plot, which shows "
            "as a fall in R^2, and scales the reported thickness."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": CBED_THICKNESS_SCHEMA,
            "thickness_angstrom": self.thickness_angstrom,
            "extinction_distance_angstrom": self.extinction_distance_angstrom,
            "first_order": self.first_order,
            "orders": [int(value) for value in self.orders],
            "excitation_errors_inv_angstrom": [
                float(value) for value in self.excitation_errors_inv_angstrom
            ],
            "r_squared": self.r_squared,
            "residuals": [float(value) for value in self.residuals],
        }


def thickness_from_fringe_minima(
    excitation_errors_inv_angstrom: ArrayLike,
    *,
    first_order: int | None = None,
    max_first_order: int = 6,
    provenance: ProvenanceRecord | None = None,
) -> TwoBeamThicknessReport:
    """Foil thickness and extinction distance from CBED fringe minima.

    What it does
        Implements the linearization of Kelly *et al.* (1975). The two-beam
        minima satisfy :math:`t\\,s_{\\mathrm{eff},n} = n` with
        :math:`s_{\\mathrm{eff}}^{2} = s^{2} + \\xi_{g}^{-2}`, so

        .. math::

           \\left(\\frac{s_{n}}{n}\\right)^{2}
             = \\frac{1}{t^{2}} - \\frac{1}{\\xi_{g}^{2}}\\,\\frac{1}{n^{2}} .

        Plotting :math:`(s_n/n)^2` against :math:`1/n^{2}` gives a straight line
        whose intercept is :math:`t^{-2}` and whose slope is
        :math:`-\\xi_{g}^{-2}`: one least-squares fit returns both unknowns.

    When to use it
        On the dark fringes of one CBED disc, simulated or measured, once their
        deviation parameters are known. This is the routine measurement behind
        every quantitative TEM technique that needs a local thickness — GND
        density, precipitate volume fraction, EELS ``t/lambda`` calibration.

    Parameters
    ----------
    excitation_errors_inv_angstrom:
        The measured minimum positions, at least two of them. Signs are
        ignored: only the magnitudes enter, and they are sorted ascending so
        the innermost minimum takes the lowest order.
    first_order:
        The integer to assign to the innermost minimum. Leave as ``None`` to
        search ``1..max_first_order`` and take the assignment giving a physical
        fit (positive intercept, negative slope) with the best :math:`R^{2}`.
    max_first_order:
        Upper bound of that search.
    provenance:
        Optional record.

    Returns
    -------
    TwoBeamThicknessReport

    Raises
    ------
    ValueError
        If fewer than two minima are supplied, or if no order assignment yields
        a physically valid fit — which happens when the data are not two-beam
        fringes at all, and is a more useful outcome than a plausible number.

    Notes
    -----
    The order assignment is the method's known weakness. When the innermost
    visible minimum is not in fact ``n = 1`` — common when the disc is small or
    the foil thin — assuming so inflates the thickness. The search here exploits
    the fact that a wrong assignment *curves* the plot, so the correct one is
    identifiable from the data whenever three or more minima are visible; with
    exactly two minima the fit is exact for every assignment and the caller must
    supply ``first_order`` from the physics.
    """

    values = np.sort(np.abs(as_float_array(excitation_errors_inv_angstrom, shape=(None,))))
    if values.size < 2:
        raise ValueError(
            "At least two fringe minima are needed: the fit determines two unknowns, "
            f"the thickness and the extinction distance, and {values.size} minimum/minima "
            "cannot."
        )
    if max_first_order < 1:
        raise ValueError("max_first_order must be at least 1.")

    candidates = [first_order] if first_order is not None else list(range(1, max_first_order + 1))
    best: tuple[float, int, float, float, np.ndarray, np.ndarray] | None = None
    for candidate in candidates:
        if candidate is None or candidate < 1:
            raise ValueError("first_order must be a positive integer.")
        orders = np.arange(candidate, candidate + values.size, dtype=np.float64)
        abscissa = 1.0 / (orders * orders)
        ordinate = np.square(values / orders)
        slope, intercept = np.polyfit(abscissa, ordinate, 1)
        if intercept <= 0.0 or slope >= 0.0:
            continue
        predicted = slope * abscissa + intercept
        residuals = ordinate - predicted
        spread = float(np.sum(np.square(ordinate - ordinate.mean())))
        r_squared = 1.0 if spread <= 0.0 else 1.0 - float(np.sum(np.square(residuals))) / spread
        if best is None or r_squared > best[0]:
            best = (r_squared, candidate, float(slope), float(intercept), orders, residuals)

    if best is None:
        raise ValueError(
            "No order assignment gives a physical two-beam fit: every candidate produced a "
            "non-positive intercept (1/t^2) or a non-negative slope (-1/xi^2). The minima "
            "supplied are not two-beam thickness fringes, or their order is beyond "
            f"max_first_order = {max_first_order}."
        )

    r_squared, candidate, slope, intercept, orders, residuals = best
    return TwoBeamThicknessReport(
        thickness_angstrom=float(1.0 / math.sqrt(intercept)),
        extinction_distance_angstrom=float(1.0 / math.sqrt(-slope)),
        first_order=int(candidate),
        orders=np.asarray(orders, dtype=np.int64),
        excitation_errors_inv_angstrom=values,
        r_squared=float(r_squared),
        residuals=np.asarray(residuals, dtype=np.float64),
        provenance=provenance,
    )


# --------------------------------------------------------------------------- #
# HOLZ rings
# --------------------------------------------------------------------------- #


def holz_ring_radii_inv_angstrom(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    beam_energy_kev: float = 200.0,
    orders: int = 2,
    max_index: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Radii of the higher-order Laue zone rings, and the layer spacing.

    What it does
        Returns the projected radius of each HOLZ ring,

        .. math::

           G_{n} \\simeq \\sqrt{\\frac{2 n H}{\\lambda}},

        where :math:`H = 1/|\\mathbf{r}_{uvw}|` is the reciprocal-lattice layer
        spacing along the zone axis, together with the layer orders that
        actually carry allowed reflections.

    When to use it
        To measure the lattice repeat **along the beam** — the one dimension a
        zone-axis pattern is blind to, because every ZOLZ reflection is
        perpendicular to the zone axis. A HOLZ ring radius converts directly
        into :math:`H`, and a change in :math:`H` is a change in that lattice
        parameter, which is how CBED measures local strain and composition.

    Parameters
    ----------
    phase, zone_axis:
        The zone axis must belong to the phase.
    beam_energy_kev:
        Accelerating voltage; the radius scales as :math:`\\lambda^{-1/2}`.
    orders:
        How many Laue zones above the zeroth to report.
    max_index:
        Bound of the reflection search used to decide which layers are
        systematically absent. A centred lattice can extinguish an entire layer,
        and reporting a ring that cannot appear would be worse than reporting
        none.

    Returns
    -------
    tuple of np.ndarray
        The layer orders that carry allowed reflections, and their radii in
        1/angstrom.

    Notes
    -----
    The square-root expression is the small-angle approximation to the
    intersection of the Ewald sphere with the ``n``-th layer, and is accurate to
    better than a percent for the first few zones at ordinary TEM voltages. The
    layer *spacing* it returns is exact.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if orders < 1:
        raise ValueError("orders must be at least 1.")

    direct = as_float_array(phase.lattice.direct_basis().matrix, shape=(3, 3))
    indices = np.asarray(zone_axis.indices, dtype=np.int64)
    divisor = int(np.gcd.reduce(np.abs(indices[indices != 0]))) if np.any(indices) else 1
    repeat = direct @ (indices.astype(np.float64) / max(divisor, 1))
    repeat_length = float(np.linalg.norm(repeat))
    if repeat_length <= 0.0:
        raise ValueError("The zone axis has zero length in the direct lattice.")
    layer_spacing = 1.0 / repeat_length

    condition = ReflectionCondition.from_phase(phase)
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid != 0, axis=1)]
    allowed = grid[centering_allowed_mask(grid, condition)]
    layer_index = allowed @ (indices // max(divisor, 1))

    wavelength = electron_wavelength_angstrom(beam_energy_kev)
    present: list[int] = []
    radii: list[float] = []
    for order in range(1, orders + 1):
        if not np.any(np.abs(layer_index) == order):
            continue
        present.append(order)
        radii.append(math.sqrt(2.0 * order * layer_spacing / wavelength))
    return np.asarray(present, dtype=np.int64), np.asarray(radii, dtype=np.float64)


# --------------------------------------------------------------------------- #
# Pattern simulation
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ConvergentBeamConfig:
    """Everything the microscope contributes to a CBED pattern.

    Purpose
    -------
    Separates the instrument settings from the crystallography, so that the same
    phase and zone can be re-simulated at another convergence angle or thickness
    without rebuilding anything else.

    Attributes
    ----------
    beam_energy_kev : float
        Accelerating voltage.
    convergence_semi_angle_mrad : float
        The half-angle of the illumination cone, ``alpha``. This is *the* CBED
        parameter: it sets the disc radius, and therefore whether the pattern is
        Kossel-Moellenstedt or Kossel.
    thickness_angstrom : float
        Foil thickness, which sets the number of fringes across a disc.
    camera_constant_mm_angstrom : float
        ``L * lambda``, the same quantity the SAED path is calibrated in.
    max_index : int
        Largest absolute Miller index enumerated.
    g_max_inv_angstrom : float
        Radial cut-off, equivalent to the recorded detector extent.
    max_excitation_error_inv_angstrom : float
        Zero-tilt selection half-width for reflections, as in
        `pytex.diffraction.kinematic`. Keeps the pattern to the zeroth Laue zone.
    disc_samples : int
        Number of samples across a disc diameter. The rocking curve is read at
        this resolution, so it also bounds how precisely fringe minima — and
        therefore a fitted thickness — can be located.
    apply_centering_absences : bool
        Remove reflections forbidden by the lattice centering.
    method : {"two-beam", "bloch"}
        How disc intensities are computed. ``"two-beam"`` gives each disc its own
        independent rocking curve — cheap, and exactly the model
        :func:`thickness_from_fringe_minima` inverts. ``"bloch"`` solves the
        coupled many-beam problem through `pytex.diffraction.dynamical`, which
        costs ``O(m n^3)`` but is the only method whose relative intensities and
        symmetry mean anything.
    absorption : AbsorptionModel or None
        The imaginary optical potential. Only accepted with ``method="bloch"``,
        because the two-beam closed form has no absorptive term; asking for it in
        the two-beam path raises rather than silently ignoring it.
    laue_zones : tuple of int
        Which Laue zones enter the **beam set**. Must contain ``0``, which
        carries the discs the pattern is drawn from. Adding ``1`` (or ``1`` and
        ``-1``) admits higher-order reflections, which is what produces HOLZ
        deficiency lines inside the discs and, crucially, what breaks the
        projection symmetry that would otherwise make every pattern look
        centrosymmetric. Only meaningful with ``method="bloch"``.
    holz_max_index : int
        Largest absolute Miller index enumerated when searching for higher-order
        reflections. These have large indices — the first HOLZ ring of a cubic
        metal sits near ``|g| ~ 5`` inverse angstrom — so the bound is separate
        from ``max_index``, which governs only the drawn zeroth-zone discs.
    holz_g_max_inv_angstrom : float
        Radial cut-off for the same search.
    """

    beam_energy_kev: float = 200.0
    convergence_semi_angle_mrad: float = 6.0
    thickness_angstrom: float = 1000.0
    camera_constant_mm_angstrom: float = 180.0
    max_index: int = 4
    g_max_inv_angstrom: float = 1.2
    max_excitation_error_inv_angstrom: float = 0.02
    disc_samples: int = 81
    apply_centering_absences: bool = True
    method: SimulationMethod = "two-beam"
    absorption: AbsorptionModel | None = None
    laue_zones: tuple[int, ...] = (0,)
    holz_max_index: int = 20
    holz_g_max_inv_angstrom: float = 6.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.beam_energy_kev) or self.beam_energy_kev <= 0.0:
            raise ValueError("beam_energy_kev must be finite and strictly positive.")
        if (
            not np.isfinite(self.convergence_semi_angle_mrad)
            or self.convergence_semi_angle_mrad <= 0.0
        ):
            raise ValueError(
                "convergence_semi_angle_mrad must be finite and strictly positive; a "
                "zero convergence angle is a parallel beam, which is SAED and is "
                "simulated by pytex.diffraction.kinematic."
            )
        if not np.isfinite(self.thickness_angstrom) or self.thickness_angstrom <= 0.0:
            raise ValueError("thickness_angstrom must be finite and strictly positive.")
        if self.camera_constant_mm_angstrom <= 0.0:
            raise ValueError("camera_constant_mm_angstrom must be strictly positive.")
        if self.max_index <= 0:
            raise ValueError("max_index must be strictly positive.")
        if self.g_max_inv_angstrom <= 0.0:
            raise ValueError("g_max_inv_angstrom must be strictly positive.")
        if self.max_excitation_error_inv_angstrom < 0.0:
            raise ValueError("max_excitation_error_inv_angstrom must be non-negative.")
        if self.disc_samples < 3:
            raise ValueError("disc_samples must be at least 3 to resolve a rocking curve.")
        if self.method not in {"two-beam", "bloch"}:
            raise ValueError(
                f"method must be 'two-beam' or 'bloch'; got '{self.method}'. Symmetry "
                "analysis requires 'bloch', because a two-beam disc is symmetric in s by "
                "construction and would report symmetry the crystal may not have."
            )
        if self.absorption is not None and self.method != "bloch":
            raise ValueError(
                "Absorption is only meaningful in the coupled calculation: the two-beam disc "
                "path evaluates a closed form that has no absorptive term. Set "
                "method='bloch' or leave absorption as None."
            )
        if not self.laue_zones:
            raise ValueError("laue_zones must name at least one Laue zone, normally 0.")
        if 0 not in self.laue_zones:
            raise ValueError(
                "laue_zones must include 0: the zeroth Laue zone carries the discs the "
                "pattern is drawn from, and a calculation without it would have nothing to "
                "plot."
            )
        if self.method != "bloch" and set(self.laue_zones) != {0}:
            raise ValueError(
                "Higher-order Laue zone beams only affect a coupled calculation. In the "
                "two-beam path each disc is independent, so a HOLZ reflection could not "
                "change it. Set method='bloch'."
            )
        if self.holz_max_index <= 0:
            raise ValueError("holz_max_index must be strictly positive.")
        if self.holz_g_max_inv_angstrom <= 0.0:
            raise ValueError("holz_g_max_inv_angstrom must be strictly positive.")

    @property
    def includes_holz_beams(self) -> bool:
        """Whether any higher-order Laue zone is in the beam set.

        The single most important thing to know before reading symmetry off a
        simulated pattern: without it, the calculation samples the *projected*
        potential, whose symmetry is at least as high as the crystal's and is
        often strictly higher.
        """

        return any(int(zone) != 0 for zone in self.laue_zones)

    @property
    def wavelength_angstrom(self) -> float:
        """Relativistic electron wavelength at this accelerating voltage."""

        return electron_wavelength_angstrom(self.beam_energy_kev)

    @property
    def convergence_semi_angle_rad(self) -> float:
        """The convergence semi-angle in radians."""

        return float(self.convergence_semi_angle_mrad) * 1e-3

    @property
    def disc_radius_inv_angstrom(self) -> float:
        """Disc radius in reciprocal space, ``alpha / lambda``."""

        return self.convergence_semi_angle_rad / self.wavelength_angstrom

    @property
    def disc_radius_mm(self) -> float:
        """Disc radius on the detector, ``(L lambda) alpha / lambda = L alpha``."""

        return self.camera_constant_mm_angstrom * self.disc_radius_inv_angstrom


@dataclass(frozen=True, slots=True)
class CBEDDisc:
    """One reflection's disc: its geometry and its rocking-curve map.

    Attributes
    ----------
    miller_indices : np.ndarray
        The reflection ``(hkl)``; ``(0, 0, 0)`` for the transmitted disc.
    centre_mm : np.ndarray
        Disc centre on the detector, in millimetres.
    radius_mm : float
        Disc radius, the same for every disc of one pattern.
    g_detector_inv_angstrom : np.ndarray
        The in-plane part of ``g``, before camera scaling.
    excitation_error_inv_angstrom : np.ndarray
        ``(m, m)`` deviation parameter across the disc, ``nan`` outside it.
    intensity : np.ndarray
        ``(m, m)`` two-beam intensity, ``nan`` outside the disc.
    extinction_distance_angstrom : float
        The ``xi_g`` used. ``inf`` for the transmitted disc, which has none.
    structure_factor_angstrom : complex
    tilt_axis_mrad : np.ndarray
        ``(m,)`` sample positions along one axis of the disc, in milliradians,
        so a profile can be plotted against the physical tilt.
    label : str
        Rendered index label, per the repository notation standard.
    """

    miller_indices: np.ndarray
    centre_mm: np.ndarray
    radius_mm: float
    g_detector_inv_angstrom: np.ndarray
    excitation_error_inv_angstrom: np.ndarray
    intensity: np.ndarray
    extinction_distance_angstrom: float
    structure_factor_angstrom: complex
    tilt_axis_mrad: np.ndarray
    label: str

    @property
    def is_transmitted(self) -> bool:
        """Whether this is the direct (000) disc."""

        return bool(np.all(self.miller_indices == 0))

    def radial_profile(self) -> tuple[np.ndarray, np.ndarray]:
        """The rocking curve along the disc's own ``g`` direction.

        Purpose
        -------
        The excitation error varies linearly across the disc along ``g``, so a
        cut along that direction *is* the rocking curve — this is the profile an
        experimenter measures to count fringes.

        Returns
        -------
        tuple of np.ndarray
            The excitation errors and intensities along the diameter parallel to
            ``g``, with the samples outside the disc removed. For the
            transmitted disc, whose own ``g`` is zero, the cut is taken along
            the first axis instead.
        """

        samples = self.intensity.shape[0]
        middle = samples // 2
        direction = np.asarray(self.g_detector_inv_angstrom, dtype=np.float64)
        along_u = self.is_transmitted or abs(direction[0]) >= abs(direction[1])
        if along_u:
            s_line = self.excitation_error_inv_angstrom[:, middle]
            intensity_line = self.intensity[:, middle]
        else:
            s_line = self.excitation_error_inv_angstrom[middle, :]
            intensity_line = self.intensity[middle, :]
        inside = np.isfinite(s_line) & np.isfinite(intensity_line)
        order = np.argsort(s_line[inside])
        return s_line[inside][order], intensity_line[inside][order]


@dataclass(frozen=True, slots=True)
class CBEDPattern:
    """A simulated convergent-beam pattern, with its limits attached.

    Purpose
    -------
    What a CBED exposure of one zone would show: the disc lattice, the rocking
    curve inside each disc, the overlap regime, and where the HOLZ rings fall.

    Limits
    ------
    They depend on ``config.method``, and the object states which was used.

    With ``"two-beam"``, **each disc is an independent calculation**. Real CBED
    is many-beam: the discs of one pattern share intensity, and at a zone axis
    the two-beam picture is at its worst. No absorption enters, so the fringes do
    not decay with thickness as they do in practice, and — the trap worth naming
    — every disc is symmetric in ``s`` by construction, so the pattern displays
    symmetry that belongs to the *method* and not to the crystal.
    :meth:`symmetry_observations` therefore refuses to run on such a pattern.

    With ``"bloch"``, the discs are mutually consistent, absorption is available,
    and higher-order Laue zone beams can be admitted. What remains unmodelled is
    inelastic background, probe aberrations, specimen bending or wedge, and the
    absorptive potential's *magnitude*, which is phenomenological.

    Attributes
    ----------
    phase : Phase
    zone_axis : ZoneAxis
    config : ConvergentBeamConfig
    discs : tuple of CBEDDisc
        The transmitted disc first, then the diffracted discs by decreasing
        structure-factor magnitude.
    zone_basis_crystal : np.ndarray
        ``(3, 3)`` with columns ``u``, ``v`` and the zone-axis unit vector.
    holz_orders : np.ndarray
    holz_radii_mm : np.ndarray
    holz_lines : HOLZLinePattern or None
        Geometry of the first-order Laue zone lines, if any fall inside the
        illumination cone. Positions are exact and independent of the simulation
        method; whether the lines are *visible* in the disc intensities depends
        on the beam set. See `pytex.diffraction.holz`.
    beam_set : BeamSet or None
        The coupled beam set, for ``method="bloch"``; ``None`` otherwise.
    nearest_disc_separation_mm : float
        Centre-to-centre distance of the closest pair of discs; ``inf`` when
        only the transmitted disc is present.
    provenance : ProvenanceRecord or None
    """

    phase: Phase
    zone_axis: ZoneAxis
    config: ConvergentBeamConfig
    discs: tuple[CBEDDisc, ...]
    zone_basis_crystal: np.ndarray
    holz_orders: np.ndarray
    holz_radii_mm: np.ndarray
    nearest_disc_separation_mm: float
    holz_lines: HOLZLinePattern | None = None
    beam_set: BeamSet | None = None
    provenance: ProvenanceRecord | None = None

    @property
    def transmitted_disc(self) -> CBEDDisc:
        """The direct (000) disc, which is always present and always first."""

        return self.discs[0]

    @property
    def is_kossel_moellenstedt(self) -> bool:
        """Whether the discs are separated, the Kossel-Moellenstedt regime.

        True when the convergence angle is small enough that neighbouring discs
        do not touch. This is the regime in which each disc can be read as an
        independent rocking curve, and therefore the regime thickness
        determination requires. Above it the discs overlap (the Kossel regime)
        and the overlapping regions show interference between beams rather than
        a single rocking curve.
        """

        return 2.0 * self.config.disc_radius_mm < self.nearest_disc_separation_mm

    @property
    def regime(self) -> str:
        """``"kossel-moellenstedt"`` or ``"kossel"``; see :attr:`is_kossel_moellenstedt`."""

        return "kossel-moellenstedt" if self.is_kossel_moellenstedt else "kossel"

    def disc_for(self, hkl: ArrayLike) -> CBEDDisc:
        """The disc of a named reflection.

        Raises
        ------
        KeyError
            If the reflection is not in the pattern — because it is forbidden,
            outside the ``g`` cut-off, or outside the excitation-error window.
        """

        wanted = np.asarray(hkl, dtype=np.int64).reshape(3)
        for disc in self.discs:
            if np.array_equal(disc.miller_indices, wanted):
                return disc
        raise KeyError(
            f"Reflection {tuple(int(v) for v in wanted)} is not in this pattern. It may be "
            "forbidden by the lattice centering, beyond g_max_inv_angstrom, or outside "
            "max_excitation_error_inv_angstrom."
        )

    # ----------------------------------------------------------------- #
    # Symmetry: the point-group determination
    # ----------------------------------------------------------------- #

    def predicted_diffraction_group(self) -> DiffractionGroup:
        """The diffraction group this zone *should* show, from the declared point group.

        Purpose
        -------
        The forward direction: what the crystallography says the pattern must
        look like. Comparing it against :meth:`symmetry_observations` is how a
        simulation is validated, and comparing it against a *recorded* pattern is
        how a candidate structure is tested.

        Returns
        -------
        DiffractionGroup

        See Also
        --------
        `pytex.diffraction.diffraction_groups.diffraction_group_for_zone_axis`
        """

        return diffraction_group_for_zone_axis(self.phase, self.zone_axis)

    def symmetry_observations(
        self,
        *,
        tolerance: float = 0.05,
        require_holz: bool = True,
        friedel_pair_two_fold: bool | None = None,
    ) -> SymmetryObservations:
        """Measure the pattern's symmetry the way an experimenter reads it.

        What it does
            Tests candidate plane operations against the computed intensities and
            keeps those that survive:

            - **Bright field.** Resamples the transmitted disc under each
              candidate operation and compares.
            - **Whole pattern.** Requires the operation to permute the disc
              centres *and* to map each disc's intensity map onto its partner's.

            The survivors are closed into a group and named. Together these two
            observations determine the diffraction group in most cases,
            centrosymmetry included: down a four-fold zone,
            ``BF = 4mm, WP = 4mm`` is ``4mm1_R`` and centric, while
            ``BF = 4mm, WP = 2mm`` is ``4_Rmm_R`` and acentric — the classic
            silicon-versus-gallium-arsenide comparison, which this measurement
            separates with residuals of ``0.00`` against ``0.32``.

        When to use it
            To close the loop — simulate a pattern, read its symmetry back, and
            hand the result to
            `pytex.diffraction.diffraction_groups.determine_point_group`. The
            same observations are what an operator records at the microscope, so
            the report format is the same either way.

        Parameters
        ----------
        tolerance:
            Largest accepted mismatch, as a fraction of each map's mean absolute
            deviation. Operations that are not grid-aligned — three- and six-fold
            rotations and their mirrors — carry resampling error that falls with
            ``config.disc_samples``; use at least 101 samples on a trigonal or
            hexagonal zone.
        require_holz:
            Refuse to report when the beam set contains no higher-order Laue zone
            reflection. Pass ``False`` only when the *projection* symmetry is
            what is wanted, and know that it is generally higher than the
            crystal's own. Without HOLZ beams, gallium arsenide down ``[001]``
            reports the ``4mm`` whole-pattern symmetry of silicon.
        friedel_pair_two_fold:
            Passed straight through to the report. **This library does not
            measure it**; see the note below.

        Returns
        -------
        SymmetryObservations

        Raises
        ------
        ValueError
            If the pattern was computed with ``method="two-beam"``, whose discs
            are symmetric in ``s`` by construction and would report symmetry that
            belongs to the method; or if ``require_holz`` is set and the beam set
            is confined to the zeroth Laue zone.

        Notes
        -----
        **Candidate operations are derived from the pattern, not assumed.**
        Rotations of order 2, 3, 4 and 6 are tested, and mirror lines are taken
        from the azimuths of the disc centres and their perpendiculars — the only
        orientations at which a mirror could permute the discs. Testing a dense
        sweep of angles instead would find spurious mirrors in a smooth map. The
        survivors are then closed under multiplication before being named,
        because ``{1, R_2, R_3, R_6}`` is four matrices but a six-fold group.

        **Why the ``+-g`` relation is not measured here.** Buxton's ``2_R``
        observation compares the ``+g`` and ``-g`` *dark-field* discs, each
        recorded with its own reflection at the Bragg condition — two exposures
        at different specimen tilts, related by the reciprocity theorem. It is
        *not* a two-fold rotation of a single zone-axis pattern. Taking it to be
        one gives a test that fails: the excitation errors satisfy
        ``s_{-g}(-theta) - s_g(theta) = -2 g_z``, which vanishes only in the
        zeroth Laue zone, so once higher-order beams are admitted the two-fold is
        broken for a centrosymmetric crystal too. That is confirmed numerically
        — the residual *grows* with the beam set for centric and acentric
        structures alike, so it is physics and not truncation — and it is why
        this method leaves the field to the caller rather than reporting a number
        that would sometimes be wrong. The determination does not need it: the
        bright-field and whole-pattern symmetries settle centrosymmetry at any
        zone whose diffraction groups differ in them.
        """

        if self.config.method != "bloch":
            raise ValueError(
                "Symmetry cannot be read from a two-beam pattern. Each disc there is an "
                "independent rocking curve, symmetric in s by construction, so every +-g "
                "pair matches and the pattern reports a centre of symmetry whatever the "
                "crystal is. Re-simulate with ConvergentBeamConfig(method='bloch')."
            )
        if require_holz and not self.config.includes_holz_beams:
            raise ValueError(
                "The beam set is confined to the zeroth Laue zone, so this calculation "
                "samples the potential projected along the beam. That projection is "
                "frequently centrosymmetric when the crystal is not - zincblende down [111] "
                "is the standard example - and reading symmetry from it would report the "
                "projection's symmetry as the crystal's. Add higher-order Laue zones with "
                "ConvergentBeamConfig(laue_zones=(0, 1, -1)), or pass require_holz=False if "
                "the projection symmetry is genuinely what is wanted."
            )
        if not 0.0 < tolerance < 1.0:
            raise ValueError("tolerance must lie strictly between 0 and 1.")
        if len(self.discs) < 2:
            raise ValueError(
                "The pattern has no diffracted discs, so there is no symmetry to read: a "
                "single transmitted disc is invariant under everything. Widen "
                "g_max_inv_angstrom or max_excitation_error_inv_angstrom until reflections "
                "are selected."
            )

        floor = self._contrast_floor()
        candidates = self._candidate_plane_operations()
        bright = [
            operation
            for operation in candidates
            if self._disc_maps_agree(
                self.transmitted_disc, self.transmitted_disc, operation, tolerance, floor
            )
        ]
        whole = [
            operation
            for operation in candidates
            if self._whole_pattern_agrees(operation, tolerance, floor)
        ]
        return SymmetryObservations(
            bright_field=plane_point_group_symbol(_close_plane_group(bright)),
            whole_pattern=plane_point_group_symbol(_close_plane_group(whole)),
            friedel_pair_two_fold=friedel_pair_two_fold,
        )

    def determine_point_group(
        self,
        *,
        tolerance: float = 0.05,
        require_holz: bool = True,
        friedel_pair_two_fold: bool | None = None,
        **kwargs: Any,
    ) -> PointGroupDetermination:
        """Measure the pattern's symmetry and turn it into a point group.

        Purpose
        -------
        The end-to-end capability in one call: simulate, read the symmetry, and
        report which crystal point groups are consistent with it — including
        whether the crystal has a centre of symmetry, which kinematic
        diffraction cannot determine at all.

        Parameters
        ----------
        tolerance, require_holz, friedel_pair_two_fold:
            As in :meth:`symmetry_observations`.
        **kwargs:
            Forwarded to
            `pytex.diffraction.diffraction_groups.determine_point_group`, chiefly
            ``candidate_point_groups`` for narrowing with prior knowledge.

        Returns
        -------
        PointGroupDetermination

        Notes
        -----
        On a *simulated* pattern this is a self-consistency check: the answer
        should contain the point group the phase declared, and
        :meth:`predicted_diffraction_group` says which diffraction group it
        should have gone through. On a *measured* pattern it is the
        determination itself.

        **Choose the zone before trusting the answer.** Two point groups that
        share a diffraction group at one beam direction generally differ at
        another; `pytex.diffraction.diffraction_groups.diffraction_group_table`
        names the directions. Down ``[001]`` the cubic pair ``m-3m`` and
        ``-43m`` differ in whole-pattern symmetry and are separated outright;
        down ``[111]`` they share both bright-field and whole-pattern symmetry
        and are not.
        """

        observations = self.symmetry_observations(
            tolerance=tolerance,
            require_holz=require_holz,
            friedel_pair_two_fold=friedel_pair_two_fold,
        )
        return determine_point_group(observations, **kwargs)

    def _candidate_plane_operations(self) -> list[np.ndarray]:
        """Plane operations that could conceivably permute the discs."""

        operations = [np.eye(2)]
        for order in (2, 3, 4, 6):
            angle = 2.0 * math.pi / order
            operations.append(
                np.array(
                    [
                        [math.cos(angle), -math.sin(angle)],
                        [math.sin(angle), math.cos(angle)],
                    ]
                )
            )
        angles: set[int] = set()
        for disc in self.discs:
            if disc.is_transmitted:
                continue
            azimuth = math.atan2(float(disc.centre_mm[1]), float(disc.centre_mm[0]))
            for offset in (0.0, math.pi / 2.0):
                angles.add(round(math.degrees((azimuth + offset) % math.pi) * 10))
        for tenths in sorted(angles):
            doubled = 2.0 * math.radians(tenths / 10.0)
            operations.append(
                np.array(
                    [
                        [math.cos(doubled), math.sin(doubled)],
                        [math.sin(doubled), -math.cos(doubled)],
                    ]
                )
            )
        return operations

    def _interior_mask(self) -> np.ndarray:
        """Samples safely inside the disc, away from the edge where resampling fails."""

        alpha = self.config.convergence_semi_angle_rad
        axis = self.transmitted_disc.tilt_axis_mrad * 1e-3
        grid_u, grid_v = np.meshgrid(axis, axis, indexing="ij")
        return np.asarray(grid_u * grid_u + grid_v * grid_v <= (0.9 * alpha) ** 2, dtype=bool)

    def _resample(self, values: np.ndarray, operation: np.ndarray) -> np.ndarray:
        """``values`` evaluated at ``operation @ theta`` on the same grid."""

        axis = self.transmitted_disc.tilt_axis_mrad * 1e-3
        samples = axis.size
        alpha = self.config.convergence_semi_angle_rad
        grid_u, grid_v = np.meshgrid(axis, axis, indexing="ij")
        stacked = np.stack([grid_u.reshape(-1), grid_v.reshape(-1)])
        mapped = operation @ stacked
        indices = (mapped + alpha) / (2.0 * alpha) * (samples - 1)
        filled = np.nan_to_num(values, nan=0.0)
        return np.asarray(
            ndimage.map_coordinates(filled, indices, order=3, mode="nearest").reshape(
                samples, samples
            ),
            dtype=np.float64,
        )

    def _contrast_floor(self) -> float:
        """The smallest per-disc contrast a symmetry test will take seriously.

        Each comparison is normalized by the disc's own mean absolute deviation,
        because the physics that breaks a symmetry frequently lives in the
        *weak* discs — for gallium arsenide down ``[001]`` the four-fold is
        broken in the near-forbidden ``{200}`` reflections, whose contrast is
        half a percent of the strongest disc's, and normalizing everything by
        the brightest disc would hide it completely.

        That per-disc normalization needs a floor, and this is it. A
        **systematically absent** reflection has an identically zero disc whose
        residual floating-point noise has a mean absolute deviation of order
        ``1e-30``; dividing by that turns rounding error into a catastrophic
        symmetry violation. Silicon down ``[001]`` has four such discs, the
        absent ``{200}``, and they were enough to destroy the four-fold before
        the floor was added. The floor is ``1e-6`` of the strongest disc's
        contrast: a disc below that carries no recordable information, so it
        cannot testify either way.
        """

        interior = self._interior_mask()
        scales = [
            float(np.mean(np.abs(values - values.mean())))
            for values in (
                np.nan_to_num(disc.intensity, nan=0.0)[interior] for disc in self.discs
            )
        ]
        return 1e-6 * max(scales) if scales else 0.0

    def _disc_maps_agree(
        self,
        source: CBEDDisc,
        target: CBEDDisc,
        operation: np.ndarray,
        tolerance: float,
        floor: float,
    ) -> bool:
        """Whether ``I_target(T theta) == I_source(theta)`` inside the disc.

        The residual is the mean absolute deviation over the disc interior,
        relative to the source disc's own contrast (floored — see
        :meth:`_contrast_floor`). A worst-case (max) criterion is not usable
        here: HOLZ lines are far narrower than the tilt sampling can resolve, so
        resampling them at a rotation that is not grid-aligned produces large
        errors along a handful of thin loci while the map as a whole is
        symmetric. An L1 measure weights those loci by their area, which is what
        an eye reading the pattern does too.
        """

        interior = self._interior_mask()
        reference = np.nan_to_num(source.intensity, nan=0.0)[interior]
        scale = max(float(np.mean(np.abs(reference - reference.mean()))), floor)
        if scale <= 0.0:
            return True
        resampled = self._resample(target.intensity, operation)[interior]
        return float(np.mean(np.abs(resampled - reference))) <= tolerance * scale

    def _whole_pattern_agrees(
        self, operation: np.ndarray, tolerance: float, floor: float
    ) -> bool:
        """Whether the operation permutes the discs and maps each map onto its partner."""

        centres = np.stack([disc.centre_mm for disc in self.discs])
        radius = max(self.config.disc_radius_mm, 1e-9)
        for index, disc in enumerate(self.discs):
            image = operation @ centres[index]
            distances = np.linalg.norm(centres - image[None, :], axis=1)
            partner = int(np.argmin(distances))
            if float(distances[partner]) > 0.05 * radius:
                return False
            if not self._disc_maps_agree(
                disc, self.discs[partner], operation, tolerance, floor
            ):
                return False
        return True

    def describe(self) -> str:
        """Convention-explicit prose: the pattern, the regime, and the limits."""

        zone_label = format_direction_indices(
            tuple(int(value) for value in self.zone_axis.indices), style="plain"
        )
        diffracted = len(self.discs) - 1
        overlap = (
            "The discs are separated (Kossel-Moellenstedt), so each one can be read as an "
            "independent rocking curve and a thickness can be fitted from its fringes."
            if self.is_kossel_moellenstedt
            else (
                "The discs overlap (Kossel regime): the convergence semi-angle exceeds half "
                "the smallest reflection spacing, so overlapping regions carry interference "
                "between beams rather than a single rocking curve, and fringe counting is "
                "no longer valid there."
            )
        )
        holz = (
            "No higher-order Laue zone carries an allowed reflection within the search "
            "bound."
            if self.holz_orders.size == 0
            else (
                "HOLZ rings at "
                + ", ".join(
                    f"n = {int(order)}: {radius:.1f} mm"
                    for order, radius in zip(self.holz_orders, self.holz_radii_mm, strict=True)
                )
                + ". Their radii measure the reciprocal-lattice layer spacing along the "
                "beam, which no zone-axis SAED pattern can see."
            )
        )
        lines = (
            ""
            if self.holz_lines is None or not self.holz_lines.bright_field_lines
            else (
                f" {len(self.holz_lines.bright_field_lines)} first-order HOLZ lines cross "
                "the bright-field disc; their positions are exact and measure the lattice "
                "parameters, subject to the wavelength degeneracy described by "
                "pytex.diffraction.holz."
            )
        )
        if self.config.method == "bloch":
            beams = self.beam_set.size if self.beam_set is not None else 0
            holz_beams = (
                int(np.count_nonzero(self.beam_set.holz_mask)) if self.beam_set is not None else 0
            )
            absorption = (
                "without absorption"
                if self.config.absorption is None or not self.config.absorption.is_absorbing
                else (
                    "with an absorptive potential of ratio "
                    f"{self.config.absorption.reflection_ratio:.3f}"
                )
            )
            method = (
                f"Intensities come from a coupled {beams}-beam Bloch-wave calculation "
                f"({holz_beams} higher-order Laue zone beams) {absorption}, so the discs of "
                "this pattern are mutually consistent. "
                + (
                    "Higher-order Laue zone beams are present, so the projection symmetry is "
                    "broken and the pattern's symmetry can be read: see "
                    "symmetry_observations() and determine_point_group()."
                    if holz_beams
                    else (
                        "The beam set is confined to the zeroth Laue zone, so the calculation "
                        "carries the symmetry of the *projected* potential, which is often "
                        "higher than the crystal's; symmetry_observations() refuses to run "
                        "on it unless asked explicitly."
                    )
                )
            )
        else:
            method = (
                "Every disc here is an independent two-beam calculation without absorption: "
                "the geometry, the regime and the fringe positions are meaningful, but the "
                "relative intensities of different discs are not, and the pattern's symmetry "
                "belongs to the method rather than to the crystal - which is why "
                "symmetry_observations() refuses to run on it. Re-simulate with "
                "ConvergentBeamConfig(method='bloch', laue_zones=(0, 1, -1)) for that."
            )
        return (
            f"Convergent-beam pattern of {self.phase.name} down {zone_label} at "
            f"{self.config.beam_energy_kev:.0f} kV with a convergence semi-angle of "
            f"{self.config.convergence_semi_angle_mrad:.2f} mrad and a foil thickness of "
            f"{self.config.thickness_angstrom:.0f} A. {diffracted} diffracted discs plus "
            f"the transmitted disc, each of radius {self.config.disc_radius_mm:.2f} mm, "
            f"with the closest disc centres {self.nearest_disc_separation_mm:.2f} mm "
            f"apart. {overlap} {holz}{lines} {method}"
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": CBED_PATTERN_SCHEMA,
            "phase": self.phase.name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "beam_energy_kev": self.config.beam_energy_kev,
            "convergence_semi_angle_mrad": self.config.convergence_semi_angle_mrad,
            "method": self.config.method,
            "laue_zones": [int(zone) for zone in self.config.laue_zones],
            "beam_count": self.beam_set.size if self.beam_set is not None else None,
            "holz_beam_count": (
                int(np.count_nonzero(self.beam_set.holz_mask))
                if self.beam_set is not None
                else None
            ),
            "absorption": (
                None
                if self.config.absorption is None
                else {
                    "mean_ratio": self.config.absorption.mean_ratio,
                    "reflection_ratio": self.config.absorption.reflection_ratio,
                }
            ),
            "bright_field_holz_line_count": (
                0 if self.holz_lines is None else len(self.holz_lines.bright_field_lines)
            ),
            "predicted_diffraction_group": self.predicted_diffraction_group().symbol,
            "thickness_angstrom": self.config.thickness_angstrom,
            "camera_constant_mm_angstrom": self.config.camera_constant_mm_angstrom,
            "disc_radius_mm": self.config.disc_radius_mm,
            "nearest_disc_separation_mm": self.nearest_disc_separation_mm,
            "regime": self.regime,
            "holz_orders": [int(value) for value in self.holz_orders],
            "holz_radii_mm": [float(value) for value in self.holz_radii_mm],
            "discs": [
                {
                    "hkl": [int(value) for value in disc.miller_indices],
                    "centre_mm": disc.centre_mm.tolist(),
                    "extinction_distance_angstrom": disc.extinction_distance_angstrom,
                    "structure_factor_modulus_angstrom": float(
                        abs(disc.structure_factor_angstrom)
                    ),
                }
                for disc in self.discs
            ],
        }


def simulate_cbed_pattern(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    config: ConvergentBeamConfig | None = None,
    in_plane_rotation_deg: float = 0.0,
    provenance: ProvenanceRecord | None = None,
) -> CBEDPattern:
    """Simulate a zone-axis convergent-beam pattern.

    What it does
        Selects the zeroth-Laue-zone reflections as the parallel-beam engine
        does, turns each into a disc of radius ``alpha / lambda``, and fills the
        disc with the intensity at the excitation error of each incident
        direction in the convergence cone — by the two-beam closed form, or, with
        ``config.method = "bloch"``, by the coupled many-beam solution of
        `pytex.diffraction.dynamical`. It also computes the first-order HOLZ line
        geometry, whose positions do not depend on the method.

    When to use it
        To see what a CBED exposure of a given phase, zone, voltage,
        convergence angle and thickness would look like; to check whether an
        intended convergence angle keeps the discs separated; to generate the
        fringes that :func:`thickness_from_fringe_minima` reads back — the round
        trip that demonstrates the thickness method end to end; and, with the
        Bloch method and higher-order Laue zones admitted, to determine the
        crystal point group from the pattern's symmetry.

    Parameters
    ----------
    phase:
        Must carry a unit cell, since the extinction distances need one.
    zone_axis:
        The beam direction in crystal indices; must belong to ``phase``.
    config:
        Instrument and specimen settings; the defaults are a 200 kV instrument,
        a 6 mrad probe and a 100 nm foil.
    in_plane_rotation_deg:
        Rotates the rendered pattern about the zone axis, as in
        `pytex.diffraction.kinematic.zone_basis_from_axis`.
    provenance:
        Optional record.

    Returns
    -------
    CBEDPattern
        With the transmitted disc first. Read
        :attr:`CBEDPattern.is_kossel_moellenstedt` before treating any disc as a
        rocking curve.

    Raises
    ------
    ValueError
        If the zone axis belongs to a different phase, or the phase has no unit
        cell.

    Notes
    -----
    **Algorithm.**

    1. Build the right-handed zone basis ``(u, v, z)`` with ``z`` along the zone
       axis, pointing toward the gun.
    2. Enumerate ``hkl`` within ``max_index``, drop those forbidden by the
       lattice centering, and keep those with ``|g| <= g_max`` and a zero-tilt
       excitation error within the window — the zeroth Laue zone.
    3. For each surviving reflection, compute ``F_g`` on an absolute scale and
       hence ``xi_g``.
    4. Sample the convergence cone on a square grid of tilts, mask to the disc,
       and evaluate ``s_g(theta) = g_z - theta . g_perp - lambda g^2 / 2``.
    5. Fill the discs.

       - ``"two-beam"``: each disc gets the closed-form rocking curve at its own
         ``s`` values, and the transmitted disc gets ``1 - I_g`` for the
         strongest diffracted reflection — the two-beam complement.
       - ``"bloch"``: one coupled calculation supplies every disc at once,
         including the transmitted one, which is then a genuine bright-field
         intensity rather than a complement. The beam set is the drawn
         zeroth-zone reflections plus, when ``laue_zones`` asks for them, the
         higher-order reflections whose Bragg loci cross the cone. Only tilts
         inside the cone are solved.

    6. Compute the first-order HOLZ line geometry (`pytex.diffraction.holz`).

    **Cost.** The Bloch path is ``O(m n^3)`` in the number of in-cone tilts ``m``
    and beams ``n``. Admitting a HOLZ ring can multiply ``n`` several-fold, so
    reduce ``disc_samples`` or tighten ``max_excitation_error_inv_angstrom``
    rather than expecting the default sampling to stay cheap.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    settings = config or ConvergentBeamConfig()
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(
            "CBED simulation needs the atom positions, because the extinction distance "
            f"does: phase '{phase.name}' carries no unit cell."
        )

    wavelength = settings.wavelength_angstrom
    zone_unit = normalize_vector(np.asarray(zone_axis.unit_vector, dtype=np.float64))
    basis = zone_basis_from_axis(zone_unit, in_plane_rotation_deg=in_plane_rotation_deg)

    limit = settings.max_index
    values = np.arange(-limit, limit + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid != 0, axis=1)]
    if settings.apply_centering_absences:
        grid = grid[centering_allowed_mask(grid, ReflectionCondition.from_phase(phase))]

    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_cartesian = grid.astype(np.float64) @ reciprocal.T
    g_zone = g_cartesian @ basis
    g_magnitude = np.linalg.norm(g_cartesian, axis=1)
    zero_tilt_s = g_zone[:, 2] - 0.5 * wavelength * g_magnitude * g_magnitude

    keep = (g_magnitude <= settings.g_max_inv_angstrom) & (
        np.abs(zero_tilt_s) <= settings.max_excitation_error_inv_angstrom
    )
    grid = grid[keep]
    g_zone = g_zone[keep]
    g_magnitude = g_magnitude[keep]

    structure_factors = (
        electron_structure_factor_angstrom(phase, grid, beam_energy_kev=settings.beam_energy_kev)
        if grid.size
        else np.zeros(0, dtype=np.complex128)
    )
    extinction = (
        extinction_distance_angstrom(phase, grid, beam_energy_kev=settings.beam_energy_kev)
        if grid.size
        else np.zeros(0, dtype=np.float64)
    )
    order = np.argsort(-np.abs(structure_factors), kind="stable") if grid.size else np.zeros(0, int)
    grid = grid[order]
    g_zone = g_zone[order]
    g_magnitude = g_magnitude[order]
    structure_factors = structure_factors[order]
    extinction = extinction[order]

    # The tilt grid: a square of incident directions, masked to the cone.
    alpha = settings.convergence_semi_angle_rad
    tilt_axis = np.linspace(-alpha, alpha, settings.disc_samples)
    tilt_u, tilt_v = np.meshgrid(tilt_axis, tilt_axis, indexing="ij")
    outside = (tilt_u * tilt_u + tilt_v * tilt_v) > alpha * alpha

    def build_disc(
        indices: np.ndarray,
        g_components: np.ndarray,
        magnitude: float,
        factor: complex,
        xi_g: float,
        intensity_map: np.ndarray,
        s_map: np.ndarray,
    ) -> CBEDDisc:
        centre = settings.camera_constant_mm_angstrom * g_components[:2]
        return CBEDDisc(
            miller_indices=np.asarray(indices, dtype=np.int64),
            centre_mm=np.asarray(centre, dtype=np.float64),
            radius_mm=settings.disc_radius_mm,
            g_detector_inv_angstrom=np.asarray(g_components[:2], dtype=np.float64),
            excitation_error_inv_angstrom=np.where(outside, np.nan, s_map),
            intensity=np.where(outside, np.nan, intensity_map),
            extinction_distance_angstrom=float(xi_g),
            structure_factor_angstrom=complex(factor),
            tilt_axis_mrad=np.asarray(tilt_axis * 1e3, dtype=np.float64),
            label=format_plane_indices(
                tuple(int(value) for value in indices), style="mathtext"
            ),
        )

    excitation_maps = [
        g_zone[row, 2]
        - tilt_u * g_zone[row, 0]
        - tilt_v * g_zone[row, 1]
        - 0.5 * wavelength * g_magnitude[row] ** 2
        for row in range(grid.shape[0])
    ]

    beam_set: BeamSet | None = None
    if settings.method == "bloch":
        intensity_maps, transmitted_intensity, beam_set = _bloch_disc_intensities(
            phase,
            zone_axis,
            settings,
            grid,
            in_plane_rotation_deg=in_plane_rotation_deg,
            tilt_u=tilt_u,
            tilt_v=tilt_v,
            outside=outside,
        )
    else:
        intensity_maps = [
            two_beam_rocking_curve(
                s_map,
                thickness_angstrom=settings.thickness_angstrom,
                extinction_distance_angstrom=float(extinction[row]),
            )
            if np.isfinite(extinction[row])
            else np.zeros_like(s_map)
            for row, s_map in enumerate(excitation_maps)
        ]
        transmitted_intensity = (
            1.0 - intensity_maps[0] if intensity_maps else np.ones_like(tilt_u)
        )

    discs = [
        build_disc(
            grid[row],
            g_zone[row],
            float(g_magnitude[row]),
            complex(structure_factors[row]),
            float(extinction[row]),
            intensity_maps[row],
            excitation_maps[row],
        )
        for row in range(grid.shape[0])
    ]
    discs.insert(
        0,
        build_disc(
            np.zeros(3, dtype=np.int64),
            np.zeros(3, dtype=np.float64),
            0.0,
            complex(0.0, 0.0),
            float("inf"),
            transmitted_intensity,
            np.zeros_like(tilt_u),
        ),
    )

    centres = np.stack([disc.centre_mm for disc in discs])
    if len(centres) > 1:
        differences = centres[:, None, :] - centres[None, :, :]
        distances = np.linalg.norm(differences, axis=2)
        np.fill_diagonal(distances, np.inf)
        nearest = float(np.min(distances))
    else:
        nearest = float("inf")

    holz_orders, holz_radii = holz_ring_radii_inv_angstrom(
        phase, zone_axis, beam_energy_kev=settings.beam_energy_kev
    )
    lines = holz_line_pattern(
        phase,
        zone_axis,
        beam_energy_kev=settings.beam_energy_kev,
        convergence_semi_angle_mrad=settings.convergence_semi_angle_mrad,
        camera_constant_mm_angstrom=settings.camera_constant_mm_angstrom,
        max_index=settings.holz_max_index,
        g_max_inv_angstrom=settings.holz_g_max_inv_angstrom,
        in_plane_rotation_deg=in_plane_rotation_deg,
        apply_centering_absences=settings.apply_centering_absences,
    )
    return CBEDPattern(
        phase=phase,
        zone_axis=zone_axis,
        config=settings,
        discs=tuple(discs),
        zone_basis_crystal=basis,
        holz_orders=holz_orders,
        holz_radii_mm=settings.camera_constant_mm_angstrom * holz_radii,
        nearest_disc_separation_mm=nearest,
        holz_lines=lines,
        beam_set=beam_set,
        provenance=provenance,
    )


def _close_plane_group(operations: list[np.ndarray]) -> np.ndarray:
    """Close a set of surviving plane operations under multiplication.

    Purpose
    -------
    The operations tested against a pattern are *generators* — one rotation per
    order, one matrix per mirror azimuth — not a group. Naming the survivors
    directly would count ``{1, R_2, R_3, R_6}`` as four rotations and report a
    four-fold axis where the crystal has a six-fold. Closing first is not a
    tidiness step: it is what makes the name correct.

    No re-verification is needed. A product of two symmetries of a map is a
    symmetry of that map, so every element of the closure is one already.

    Raises
    ------
    ValueError
        If the closure exceeds twelve elements, the largest crystallographic
        plane point group. That means a spurious operation survived the
        intensity test, which is a sampling or tolerance problem and should be
        reported rather than named.
    """

    closure: dict[tuple[float, ...], np.ndarray] = {}
    for operation in [np.eye(2), *operations]:
        closure.setdefault(tuple(np.round(np.asarray(operation).reshape(-1), 6) + 0.0), operation)
    changed = True
    while changed:
        changed = False
        for first in list(closure.values()):
            for second in list(closure.values()):
                product = first @ second
                key = tuple(np.round(product.reshape(-1), 6) + 0.0)
                if key not in closure:
                    closure[key] = product
                    changed = True
        if len(closure) > 12:
            raise ValueError(
                f"The measured symmetry operations generate {len(closure)} elements, more "
                "than the twelve of the largest crystallographic plane point group. A "
                "spurious operation passed the intensity test: increase disc_samples so the "
                "sharpest features are resolved, or tighten the tolerance."
            )
    return np.stack(list(closure.values()))


def _bloch_disc_intensities(
    phase: Phase,
    zone_axis: ZoneAxis,
    settings: ConvergentBeamConfig,
    drawn_indices: np.ndarray,
    *,
    in_plane_rotation_deg: float,
    tilt_u: np.ndarray,
    tilt_v: np.ndarray,
    outside: np.ndarray,
) -> tuple[list[np.ndarray], np.ndarray, BeamSet]:
    """Disc intensities from one coupled many-beam calculation.

    The beam set is the drawn zeroth-zone reflections plus, when
    ``settings.laue_zones`` asks for them, the higher-order reflections whose
    Bragg loci cross the illumination cone. Only tilts inside the cone are
    solved: the square sampling grid wastes 21 percent of its points on the
    corners, and at ``O(n^3)`` per tilt that is worth not paying.
    """

    holz_zones = tuple(int(zone) for zone in settings.laue_zones if int(zone) != 0)
    indices = drawn_indices
    if holz_zones:
        higher = beam_set_for_zone(
            phase,
            zone_axis,
            beam_energy_kev=settings.beam_energy_kev,
            max_index=settings.holz_max_index,
            g_max_inv_angstrom=settings.holz_g_max_inv_angstrom,
            max_excitation_error_inv_angstrom=settings.max_excitation_error_inv_angstrom,
            convergence_semi_angle_mrad=settings.convergence_semi_angle_mrad,
            laue_zones=holz_zones,
            in_plane_rotation_deg=in_plane_rotation_deg,
            apply_centering_absences=settings.apply_centering_absences,
        )
        indices = np.concatenate([drawn_indices, higher.miller_indices], axis=0)

    beams = beam_set_from_indices(
        phase,
        zone_axis,
        indices,
        beam_energy_kev=settings.beam_energy_kev,
        in_plane_rotation_deg=in_plane_rotation_deg,
    )

    inside = ~outside
    tilts = np.stack([tilt_u[inside], tilt_v[inside]], axis=1)
    solution = solve_bloch_waves(
        beams,
        tilts,
        thickness_angstrom=settings.thickness_angstrom,
        absorption=settings.absorption,
    )

    def scatter(column: np.ndarray) -> np.ndarray:
        filled = np.zeros_like(tilt_u)
        filled[inside] = column
        return filled

    maps = [
        scatter(solution.intensities[:, beams.index_of(row)]) for row in drawn_indices
    ]
    return maps, scatter(solution.transmitted_intensity), beams
