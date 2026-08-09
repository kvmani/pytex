"""Convergent-beam electron diffraction: discs, rocking curves, and thickness.

Selected-area diffraction illuminates the specimen with a *parallel* beam, so a
reflection is a spot and the pattern carries the projected symmetry of the zone
and nothing else. Converge the beam instead — focus it to a probe with a
semi-angle of a few milliradians — and every spot spreads into a **disc**, with
each point of the disc corresponding to a different incident direction. The disc
is therefore a *map of the rocking curve*: one CBED exposure contains the
intensity as a function of deviation from the Bragg condition, which a parallel
beam could only obtain by tilting the specimen through a series of exposures.

Three things follow, and this module implements the first two of them.

**Thickness.** The two-beam rocking curve has zeros at known deviations, and
their spacing depends on the foil thickness. Measuring the fringe positions in a
single disc yields the local thickness — and, in the same fit, the extinction
distance, so the answer does not depend on a tabulated constant. This is the
standard Kelly *et al.* analysis, implemented as
:func:`thickness_from_fringe_minima`.

**Lattice parameter along the beam.** The convergence also makes higher-order
Laue zones visible as rings, whose radii depend on the reciprocal-lattice layer
spacing along the zone axis — the one lattice dimension a zone-axis SAED pattern
cannot see at all. See :func:`holz_ring_radii_inv_angstrom`.

**Point and space group.** The symmetry *within* the discs and of the whole
pattern determines the diffraction group, and hence the point group including
the presence or absence of a centre of symmetry, which Friedel's law hides from
kinematic SAED. That determination needs full dynamical intensities and is
**not** implemented here; the pattern object says so rather than implying it.

What this module models, and what it does not
---------------------------------------------

Modelled: disc geometry from the convergence semi-angle, the excitation error
across each disc, the two-beam dynamical rocking curve with its extinction
distance, the Kossel-Moellenstedt versus Kossel overlap regime, HOLZ ring radii,
and thickness extraction from fringe minima.

Not modelled: many-beam dynamical coupling (each disc is computed in its own
two-beam approximation, so the discs of one pattern are not mutually
consistent), absorption, HOLZ *line* positions within the bright-field disc,
inelastic background, and probe aberrations. Those limits are stated on
:class:`CBEDPattern` and repeated by its ``describe()``, so a caller cannot
mistake a simulated disc for a dynamical calculation.

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
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.notation import format_direction_indices, format_plane_indices
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.kinematic import (
    centering_allowed_mask,
    electron_wavelength_angstrom,
    zone_basis_from_axis,
)
from pytex.diffraction.physics import ReflectionCondition
from pytex.diffraction.scattering import electron_scattering_factors

__all__ = [
    "CBED_PATTERN_SCHEMA",
    "CBED_THICKNESS_SCHEMA",
    "CBEDDisc",
    "CBEDPattern",
    "ConvergentBeamConfig",
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

#: Electron rest energy in keV, for the relativistic correction ``gamma``.
_ELECTRON_REST_ENERGY_KEV = 510.99895

_DEBYE_WALLER_DENOMINATOR = 16.0 * np.pi * np.pi


# --------------------------------------------------------------------------- #
# Structure factors on an absolute scale, and extinction distances
# --------------------------------------------------------------------------- #


def electron_structure_factor_angstrom(
    phase: Phase,
    hkl: ArrayLike,
    *,
    beam_energy_kev: float = 200.0,
) -> np.ndarray:
    """Electron structure factors ``F_g`` in angstrom, relativistically corrected.

    What it does
        Sums the Mott-Bethe electron scattering factors over the unit cell,

        .. math::

           F_{g} = \\gamma \\sum_{j} o_{j}\\, f_{e}^{j}(s)\\,
                   e^{-B_{j}s^{2}}\\, e^{2\\pi i\\,\\mathbf{g}\\cdot\\mathbf{r}_{j}},
           \\qquad s = |\\mathbf{g}|/2,

        with the relativistic factor
        :math:`\\gamma = 1 + E/m_{0}c^{2}` applied because the incident electron
        is not slow: at 200 kV it is 1.39, and omitting it would make every
        extinction distance 39 percent too long.

    When to use it
        Whenever the *absolute* scale matters — extinction distances, dynamical
        rocking curves, thickness determination. For relative spot intensities
        within one zone, `pytex.diffraction.kinematic.electron_structure_factors`
        is the cheaper atomic-number proxy and is sufficient.

    Parameters
    ----------
    phase:
        Must carry a unit cell with sites; a bare lattice has no structure
        factor and raises rather than returning a fabricated one.
    hkl:
        ``(n, 3)`` integer Miller indices, or one triple.
    beam_energy_kev:
        Accelerating voltage, for the wavelength and for ``gamma``.

    Returns
    -------
    np.ndarray
        ``(n,)`` complex structure factors in angstrom.

    Raises
    ------
    ValueError
        If the phase has no unit-cell sites.
    """

    indices = np.atleast_2d(np.asarray(hkl, dtype=np.int64))
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (3,) or (n, 3).")
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(
            "Electron structure factors on an absolute scale need the atom positions: "
            f"phase '{phase.name}' carries no unit cell. Load the phase from a CIF, or "
            "use the relative-intensity path in pytex.diffraction.kinematic."
        )
    if not np.isfinite(beam_energy_kev) or beam_energy_kev <= 0.0:
        raise ValueError("beam_energy_kev must be finite and strictly positive.")

    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_cartesian = indices.astype(np.float64) @ reciprocal.T
    s_values = np.linalg.norm(g_cartesian, axis=1) / 2.0

    gamma = 1.0 + float(beam_energy_kev) / _ELECTRON_REST_ENERGY_KEV
    total = np.zeros(indices.shape[0], dtype=np.complex128)
    for site in phase.unit_cell.sites:
        factors = electron_scattering_factors(site.species, s_values)
        b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
        damping = np.exp(-b_iso * s_values * s_values)
        phase_argument = indices.astype(np.float64) @ site.fractional_coordinates
        total += (
            float(site.occupancy)
            * factors
            * damping
            * np.exp(2.0j * np.pi * phase_argument)
        )
    return np.asarray(gamma * total, dtype=np.complex128)


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
    **Each disc is an independent two-beam calculation.** Real CBED is
    many-beam: the discs of one pattern share intensity, and at a zone axis
    the two-beam picture is at its worst. Absorption is not modelled, so the
    fringes do not decay with thickness as they do in practice. HOLZ *lines*
    inside the bright-field disc — the sharp features used for lattice-parameter
    metrology — are not simulated; only the ring radii are given. And the
    symmetry analysis that determines the diffraction group, and with it the
    point group, is not implemented at all.

    What the object is therefore good for: geometry, regime, fringe counting,
    thickness methodology, and teaching. What it is not good for: quantitative
    intensity comparison against an experimental pattern.

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
        return (
            f"Convergent-beam pattern of {self.phase.name} down {zone_label} at "
            f"{self.config.beam_energy_kev:.0f} kV with a convergence semi-angle of "
            f"{self.config.convergence_semi_angle_mrad:.2f} mrad and a foil thickness of "
            f"{self.config.thickness_angstrom:.0f} A. {diffracted} diffracted discs plus "
            f"the transmitted disc, each of radius {self.config.disc_radius_mm:.2f} mm, "
            f"with the closest disc centres {self.nearest_disc_separation_mm:.2f} mm "
            f"apart. {overlap} {holz} "
            "Every disc here is an independent two-beam calculation without absorption: "
            "the geometry, the regime and the fringe positions are meaningful, but the "
            "relative intensities of different discs are not, and the diffraction-group "
            "symmetry determination is not implemented."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": CBED_PATTERN_SCHEMA,
            "phase": self.phase.name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "beam_energy_kev": self.config.beam_energy_kev,
            "convergence_semi_angle_mrad": self.config.convergence_semi_angle_mrad,
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
        disc with the two-beam rocking curve evaluated at the excitation error
        of each incident direction in the convergence cone.

    When to use it
        To see what a CBED exposure of a given phase, zone, voltage,
        convergence angle and thickness would look like; to check whether an
        intended convergence angle keeps the discs separated; and to generate
        the fringes that :func:`thickness_from_fringe_minima` reads back — the
        round trip that demonstrates the thickness method end to end.

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
    5. Fill the disc with the two-beam intensity at those ``s`` values. The
       transmitted disc is filled with ``1 - I_g`` for the strongest diffracted
       reflection, the two-beam complement.
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

    discs: list[CBEDDisc] = []
    strongest_intensity: np.ndarray | None = None
    for row in range(grid.shape[0]):
        s_map = (
            g_zone[row, 2]
            - tilt_u * g_zone[row, 0]
            - tilt_v * g_zone[row, 1]
            - 0.5 * wavelength * g_magnitude[row] ** 2
        )
        xi_g = float(extinction[row])
        intensity_map = (
            two_beam_rocking_curve(
                s_map,
                thickness_angstrom=settings.thickness_angstrom,
                extinction_distance_angstrom=xi_g,
            )
            if np.isfinite(xi_g)
            else np.zeros_like(s_map)
        )
        if strongest_intensity is None:
            strongest_intensity = intensity_map
        discs.append(
            build_disc(
                grid[row], g_zone[row], float(g_magnitude[row]),
                complex(structure_factors[row]), xi_g, intensity_map, s_map,
            )
        )

    transmitted_intensity = (
        1.0 - strongest_intensity if strongest_intensity is not None
        else np.ones_like(tilt_u)
    )
    transmitted = build_disc(
        np.zeros(3, dtype=np.int64),
        np.zeros(3, dtype=np.float64),
        0.0,
        complex(0.0, 0.0),
        float("inf"),
        transmitted_intensity,
        np.zeros_like(tilt_u),
    )
    discs.insert(0, transmitted)

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
    return CBEDPattern(
        phase=phase,
        zone_axis=zone_axis,
        config=settings,
        discs=tuple(discs),
        zone_basis_crystal=basis,
        holz_orders=holz_orders,
        holz_radii_mm=settings.camera_constant_mm_angstrom * holz_radii,
        nearest_disc_separation_mm=nearest,
        provenance=provenance,
    )
