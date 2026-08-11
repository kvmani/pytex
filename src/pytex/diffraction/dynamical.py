"""Many-beam dynamical electron scattering: Bloch waves with absorption.

A two-beam calculation answers "how does reflection ``g`` behave if nothing else
is excited?". At a zone axis that premise is false: dozens of reflections are
simultaneously near the Bragg condition, they exchange intensity with one
another, and the answer for any one of them depends on all the others. This
module solves the coupled problem instead, by the Bloch-wave method, and adds
the imaginary part of the crystal potential that makes the result a physically
realizable one.

The equations
-------------

In the column approximation the beam amplitudes
:math:`\\psi = (\\psi_{0}, \\psi_{g_1}, \\dots)` obey

.. math::

   \\frac{\\mathrm{d}\\psi}{\\mathrm{d}z} = i\\pi\\,\\mathbf{A}\\,\\psi,
   \\qquad
   A_{gg} = 2 s_{g} + \\frac{i}{\\xi'_{0}},
   \\qquad
   A_{gh} = \\nu_{g-h} + \\frac{i}{\\xi'_{g-h}}\\;(g \\neq h),

where :math:`\\nu_{g} = \\lambda F_{g} / (\\pi V_{c}\\cos\\theta_{g})` is the
complex Fourier coefficient of the scaled lattice potential, so that
:math:`|\\nu_{g}| = 1/\\xi_{g}` reproduces the extinction distance of
`pytex.diffraction.cbed.extinction_distance_angstrom` exactly. With ``A``
independent of ``z`` the solution is one matrix exponential,

.. math::

   \\psi(t) = \\exp(i\\pi\\mathbf{A}t)\\,\\psi(0),
   \\qquad \\psi(0) = \\mathbf{e}_{0},

evaluated here by eigen-decomposition, because the eigenvalues
:math:`\\gamma_{j}` *are* the Bloch-wave excitations — the dispersion surface —
and their imaginary parts are the per-branch absorption coefficients that
produce the Borrmann effect.

Three properties are exact and are asserted by the tests rather than assumed:

- **Two beams reproduce the closed form.** For a single reflection the
  eigenvalues are :math:`s \\pm \\sqrt{s^{2} + \\xi_{g}^{-2}}` and the solution
  is :math:`I_{g} = \\sin^{2}(\\pi t s_{\\mathrm{eff}})/(\\xi_{g}s_{\\mathrm{eff}})^{2}`,
  the expression `pytex.diffraction.cbed.two_beam_rocking_curve` evaluates.
- **Without absorption, intensity is conserved.** The elastic ``A`` is
  Hermitian, so :math:`\\exp(i\\pi\\mathbf{A}t)` is unitary and
  :math:`\\sum_{g} I_{g} = 1` at every thickness and every tilt.
- **Normal absorption is a scalar.** The term :math:`i/\\xi'_{0}` sits on every
  diagonal element, so it factors out of the exponential as
  :math:`\\exp(-2\\pi t/\\xi'_{0})` and cannot change a relative intensity, a
  fringe position, or a symmetry. Only the *off-diagonal* absorptive terms do
  that, and what they produce is anomalous absorption.

Friedel's law, and why CBED escapes it
--------------------------------------

:math:`A_{gh} = \\nu_{g-h}` and :math:`A_{hg} = \\nu_{h-g} = \\nu_{g-h}^{*}` for
a real potential, so ``A`` is Hermitian but **symmetric only when every**
:math:`\\nu_{g}` **in the beam set is real** — which happens exactly when the
structure those beams sample has a centre of symmetry, with the origin on it.
Relabelling the beams by :math:`g \\mapsto -g` and the incident tilt by
:math:`\\boldsymbol{\\theta} \\mapsto -\\boldsymbol{\\theta}` turns ``A`` into
:math:`\\mathbf{A}^{\\mathsf{T}}` (exactly, for a zeroth-Laue-zone set, because
:math:`s_{g}(\\boldsymbol{\\theta}) = s_{-g}(-\\boldsymbol{\\theta})` there), so
:math:`M = \\exp(i\\pi\\mathbf{A}t)` becomes :math:`M^{\\mathsf{T}}` and

    :math:`I_{g}(\\boldsymbol{\\theta}) = I_{-g}(-\\boldsymbol{\\theta})
    \\iff |M_{g0}| = |M_{0g}| \\iff M` symmetric :math:`\\iff` the sampled
    structure is centrosymmetric.

This is Friedel's law recovered as a *theorem about the propagator* rather than
as a kinematic accident, and it is what
`pytex.diffraction.diffraction_groups` turns into a point-group determination.

The word "sampled" carries the whole difficulty of the technique. A beam set
confined to the zeroth Laue zone samples the potential *projected* along the
beam, and that projection is frequently centrosymmetric even when the crystal is
not: for zincblende down :math:`[111]` every ZOLZ coefficient is real, so a
projection calculation reports Friedel's law to machine precision and cannot see
the polarity. Admit the first-order Laue zone and the coefficients acquire
phases, the propagator loses its symmetry, and the same pair of discs differs by
tens of percent. **Higher-order Laue zone interaction is not a refinement here;
it is the entire mechanism**, and a symmetry conclusion drawn from a ZOLZ-only
calculation is worthless. :meth:`BeamSet.holz_mask` exists so that this can be
checked rather than assumed.

What this module does not do
----------------------------

Every beam in the set is solved exactly; there is no Bethe perturbation of weak
beams, so a calculation that needs a whole HOLZ ring pays for it in full
(:math:`O(m n^{3})` for ``m`` incident directions and ``n`` beams). The economy
that works is a tighter ``max_excitation_error_inv_angstrom``, which removes
beams that were barely coupled, rather than coarser tilt sampling.

The specimen is a perfect parallel-sided slab in the column approximation:
no wedge, no bending, no strain gradient, no surface relaxation, and no probe
aberration. The absorptive potential is phenomenological — see
:class:`AbsorptionModel`.

See Also
--------
`pytex.diffraction.cbed` : disc geometry and the two-beam surface this extends.
`pytex.diffraction.holz` : the higher-order reflections whose lines this feeds.
`pytex.diffraction.diffraction_groups` : the symmetry analysis built on it.
`docs/site/theory/dynamical_cbed_and_symmetry_determination.md` : derivations.
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
from pytex.diffraction.scattering import electron_structure_factor_angstrom

__all__ = [
    "BLOCH_WAVE_SOLUTION_SCHEMA",
    "AbsorptionModel",
    "BeamSet",
    "BlochWaveSolution",
    "beam_set_for_zone",
    "beam_set_from_indices",
    "potential_coefficients_inv_angstrom",
    "solve_bloch_waves",
    "structure_matrix",
]

#: Schema identifier of the Bloch-wave solution payload.
BLOCH_WAVE_SOLUTION_SCHEMA = "pytex.bloch_wave_solution/1"

_MINIMUM_STRUCTURE_FACTOR = 1e-12


# --------------------------------------------------------------------------- #
# Absorption
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class AbsorptionModel:
    """The imaginary part of the crystal potential, as absorption ratios.

    Purpose
    -------
    Electrons are removed from the coherent elastic wavefield by processes the
    Bloch-wave calculation does not follow — chiefly thermal diffuse scattering,
    and then plasmon and core losses. Their effect is represented, exactly as in
    the standard treatment, by adding an imaginary term to every Fourier
    coefficient of the potential: :math:`\\xi_{g}^{-1} \\to \\xi_{g}^{-1} +
    i\\,\\xi_{g}'^{-1}`. This object carries the two numbers that fix those
    imaginary parts.

    What is and is not claimed
    --------------------------
    The *structure* of the model is not an approximation: an imaginary optical
    potential is the correct way to represent loss from a coherent wavefield,
    and everything interesting — the Borrmann effect, the asymmetry of the
    bright-field rocking curve, the decay of thickness fringes — follows from it
    as a **derived** consequence of the eigenvector structure, not as an applied
    correction.

    The *magnitudes* are phenomenological. They are given here as the ratios
    :math:`\\xi_{0}/\\xi'_{0}` and :math:`\\xi_{g}/\\xi'_{g}`, whose customary
    working value for a metal near 100-200 kV is about ``0.1`` (Hirsch, Howie,
    Nicholson, Pashley and Whelan, *Electron Microscopy of Thin Crystals*,
    Ch. 12). A first-principles absorptive form factor — the Einstein-model
    thermal-diffuse integral of Hall and Hirsch, as parametrized by Bird and
    King — is **not** implemented, so a ratio taken from this default should be
    treated as an order-of-magnitude physical model and not as a measured
    constant. Where the answer must not depend on it, say so: see
    :attr:`BlochWaveSolution.normal_absorption_factor`.

    Attributes
    ----------
    mean_ratio : float
        :math:`\\xi_{0}/\\xi'_{0}`, the *normal* absorption. It multiplies the
        whole wavefield by :math:`\\exp(-2\\pi t/\\xi'_{0})` and can therefore
        change no relative intensity, no fringe position and no symmetry.
    reflection_ratio : float
        :math:`\\xi_{g}/\\xi'_{g}`, the *anomalous* absorption. This one is
        responsible for every qualitative effect of absorption, because it
        absorbs the Bloch-wave branches at different rates.

    Raises
    ------
    ValueError
        If either ratio is negative, or if ``reflection_ratio`` exceeds
        ``mean_ratio``. The latter is not a stylistic restriction: the
        absorptive matrix has eigenvalues :math:`\\xi_{0}'^{-1} \\pm
        \\xi_{g}'^{-1}`, and since :math:`\\xi_{0} < \\xi_{g}` always, a
        reflection ratio larger than the mean ratio can make one of them
        negative — a Bloch wave that *gains* intensity, which no absorption
        process can do.

    Examples
    --------
    >>> AbsorptionModel.none().is_absorbing
    False
    >>> AbsorptionModel().reflection_ratio
    0.1
    """

    mean_ratio: float = 0.1
    reflection_ratio: float = 0.1

    def __post_init__(self) -> None:
        if not np.isfinite(self.mean_ratio) or self.mean_ratio < 0.0:
            raise ValueError("mean_ratio must be finite and non-negative.")
        if not np.isfinite(self.reflection_ratio) or self.reflection_ratio < 0.0:
            raise ValueError("reflection_ratio must be finite and non-negative.")
        if self.reflection_ratio > self.mean_ratio:
            raise ValueError(
                "reflection_ratio must not exceed mean_ratio: the absorptive matrix would "
                "then have a negative eigenvalue, meaning a Bloch wave that gains intensity "
                f"with depth. Got reflection_ratio={self.reflection_ratio} > "
                f"mean_ratio={self.mean_ratio}."
            )

    @classmethod
    def none(cls) -> AbsorptionModel:
        """The absorption-free model, for checking intensity conservation.

        With it, the structure matrix is Hermitian and the propagator unitary,
        so :math:`\\sum_{g} I_{g} = 1` exactly. That identity is the numerical
        check that the many-beam solver is correct, and it is unavailable once
        absorption is switched on.
        """

        return cls(mean_ratio=0.0, reflection_ratio=0.0)

    @property
    def is_absorbing(self) -> bool:
        """Whether any absorption at all is applied."""

        return self.mean_ratio > 0.0 or self.reflection_ratio > 0.0

    def describe(self) -> str:
        """Convention-explicit prose: what the model does and what it cannot claim."""

        if not self.is_absorbing:
            return (
                "Absorption is switched off. The structure matrix is Hermitian, the "
                "propagator is unitary, and the beam intensities sum to exactly one at every "
                "thickness. Thickness fringes therefore never decay, which is the one "
                "respect in which the calculation is visibly not an experiment."
            )
        return (
            f"Absorption as an imaginary optical potential with xi_0/xi'_0 = "
            f"{self.mean_ratio:.3f} (normal) and xi_g/xi'_g = {self.reflection_ratio:.3f} "
            "(anomalous). Normal absorption multiplies the entire wavefield by "
            "exp(-2 pi t / xi'_0) and changes no relative intensity, fringe position or "
            "symmetry; the anomalous term absorbs the Bloch-wave branches at different "
            "rates and is what makes the bright-field rocking curve asymmetric while the "
            "dark-field one stays symmetric. The ratios are the customary phenomenological "
            "values of Hirsch et al. rather than computed absorptive form factors, so "
            "quantitative absorption strengths from this model are indicative only."
        )


# --------------------------------------------------------------------------- #
# Potential coefficients
# --------------------------------------------------------------------------- #


def potential_coefficients_inv_angstrom(
    phase: Phase,
    hkl: ArrayLike,
    *,
    beam_energy_kev: float = 200.0,
) -> np.ndarray:
    """Complex Fourier coefficients ``nu_g`` of the scaled lattice potential.

    What it does
        Returns

        .. math::

           \\nu_{g} = \\frac{\\lambda F_{g}}{\\pi V_{c}\\cos\\theta_{g}},

        the off-diagonal element of the dynamical structure matrix. Its modulus
        is the reciprocal extinction distance, :math:`|\\nu_{g}| = 1/\\xi_{g}`,
        so this is the same physical quantity that
        `pytex.diffraction.cbed.extinction_distance_angstrom` returns — but
        **with its phase kept**. The phase is not decoration: it is what
        distinguishes a centrosymmetric structure (all :math:`\\nu_{g}` real)
        from a non-centrosymmetric one, and therefore what makes a point-group
        determination possible.

    When to use it
        Whenever the coupling between two beams is needed rather than the
        strength of one reflection: building a structure matrix, judging which
        reflections couple strongly enough to keep in a beam set, or inspecting
        the structure-factor phases of a candidate space group.

    Parameters
    ----------
    phase:
        Must carry a unit cell with sites.
    hkl:
        ``(n, 3)`` integer Miller indices, or one triple. ``(0, 0, 0)`` is
        allowed and returns the mean inner potential coefficient
        :math:`\\nu_{0} = 1/\\xi_{0}`.
    beam_energy_kev:
        Accelerating voltage; enters through the wavelength and the
        relativistic correction inside :math:`F_{g}`.

    Returns
    -------
    np.ndarray
        ``(n,)`` complex coefficients in 1/angstrom.

    Notes
    -----
    The :math:`\\cos\\theta_{g}` factor is carried so that the two-beam limit of
    the many-beam solver reproduces
    `pytex.diffraction.cbed.two_beam_rocking_curve` to machine precision rather
    than to within a Bragg-angle correction. Because it depends only on
    :math:`|g|`, keeping it does not spoil the index symmetry
    :math:`\\nu_{-g} = \\nu_{g}^{*}` that the Hermiticity of the structure
    matrix rests on.

    See Also
    --------
    `pytex.diffraction.scattering.electron_structure_factor_angstrom` : supplies
        :math:`F_{g}`.
    """

    indices = np.atleast_2d(np.asarray(hkl, dtype=np.int64))
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (3,) or (n, 3).")

    structure_factors = electron_structure_factor_angstrom(
        phase, indices, beam_energy_kev=beam_energy_kev
    )
    wavelength = electron_wavelength_angstrom(beam_energy_kev)
    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_magnitude = np.linalg.norm(indices.astype(np.float64) @ reciprocal.T, axis=1)
    sin_theta = np.clip(wavelength * g_magnitude / 2.0, -1.0, 1.0)
    cos_theta = np.sqrt(np.clip(1.0 - sin_theta * sin_theta, 1e-12, 1.0))

    volume = abs(float(np.linalg.det(phase.lattice.direct_basis().matrix)))
    coefficients = wavelength * structure_factors / (np.pi * volume * cos_theta)
    return np.asarray(coefficients, dtype=np.complex128)


# --------------------------------------------------------------------------- #
# Beam selection
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class BeamSet:
    """The reflections a dynamical calculation follows, in a fixed order.

    Purpose
    -------
    A many-beam calculation is only defined once the beams are named, and the
    choice is a scientific one: too few and the coupling is wrong, too many and
    nothing is gained but cost. This object records exactly which reflections
    were kept, which Laue zone each belongs to, and the geometry needed to
    evaluate excitation errors at any incident tilt — so a result can be
    reproduced and a convergence test can be stated.

    The transmitted beam ``(0 0 0)`` is always present and always first, because
    it is the incident boundary condition.

    Attributes
    ----------
    phase : Phase
    zone_axis : ZoneAxis
        Beam direction in crystal indices; the unit vector points toward the
        gun, as everywhere else in `pytex.diffraction`.
    beam_energy_kev : float
    miller_indices : np.ndarray
        ``(n, 3)`` integers, ``(0, 0, 0)`` first.
    g_zone : np.ndarray
        ``(n, 3)`` components of each ``g`` in the zone basis: the first two are
        the in-plane part seen on the detector, the third is the component along
        the beam.
    g_magnitude_inv_angstrom : np.ndarray
        ``(n,)`` lengths of ``g``.
    laue_zone : np.ndarray
        ``(n,)`` integer Laue-zone order ``n = g . [uvw]`` with ``[uvw]``
        reduced by its greatest common divisor. ``0`` is the ZOLZ.
    zone_basis_crystal : np.ndarray
        ``(3, 3)`` with columns ``u``, ``v`` and the zone-axis unit vector.
    """

    phase: Phase
    zone_axis: ZoneAxis
    beam_energy_kev: float
    miller_indices: np.ndarray
    g_zone: np.ndarray
    g_magnitude_inv_angstrom: np.ndarray
    laue_zone: np.ndarray
    zone_basis_crystal: np.ndarray

    @property
    def size(self) -> int:
        """Number of beams, including the transmitted one."""

        return int(self.miller_indices.shape[0])

    @property
    def wavelength_angstrom(self) -> float:
        """Relativistic electron wavelength at this accelerating voltage."""

        return electron_wavelength_angstrom(self.beam_energy_kev)

    @property
    def zolz_mask(self) -> np.ndarray:
        """``(n,)`` boolean: beams in the zeroth Laue zone."""

        return np.asarray(self.laue_zone == 0, dtype=bool)

    @property
    def holz_mask(self) -> np.ndarray:
        """``(n,)`` boolean: beams in a higher-order Laue zone.

        A calculation with no ``True`` here is a *projection* calculation: it
        carries an extra symmetry that the real crystal does not have, which is
        the single most important caveat in CBED symmetry work. See
        `pytex.diffraction.diffraction_groups`.
        """

        return np.asarray(self.laue_zone != 0, dtype=bool)

    def index_of(self, hkl: ArrayLike) -> int:
        """Position of a reflection in the beam ordering.

        Raises
        ------
        KeyError
            If the reflection is not in the set.
        """

        wanted = np.asarray(hkl, dtype=np.int64).reshape(3)
        matches = np.flatnonzero(np.all(self.miller_indices == wanted, axis=1))
        if matches.size == 0:
            raise KeyError(
                f"Reflection {tuple(int(value) for value in wanted)} is not in this beam set. "
                "It was excluded by the centering absences, the |g| cut-off, or the "
                "excitation-error window."
            )
        return int(matches[0])

    def excitation_errors(self, tilts_rad: ArrayLike) -> np.ndarray:
        """Excitation errors of every beam at every incident tilt.

        What it does
            Evaluates

            .. math::

               s_{g}(\\boldsymbol{\\theta}) = g_{z}
                   - \\theta_{x}g_{u} - \\theta_{y}g_{v}
                   - \\tfrac{1}{2}\\lambda|\\mathbf{g}|^{2},

            the same expression `pytex.diffraction.cbed` uses for a disc, so a
            dynamical calculation and a two-beam disc are evaluated at
            identical geometry and any difference between them is physics.

        Parameters
        ----------
        tilts_rad:
            ``(m, 2)`` incident-beam tilts in radians, or one pair.

        Returns
        -------
        np.ndarray
            ``(m, n)`` excitation errors in 1/angstrom.
        """

        tilts = np.atleast_2d(as_float_array(tilts_rad, shape=(None, 2)))
        wavelength = self.wavelength_angstrom
        zero_tilt = self.g_zone[:, 2] - 0.5 * wavelength * np.square(
            self.g_magnitude_inv_angstrom
        )
        return np.asarray(
            zero_tilt[None, :] - tilts @ self.g_zone[:, :2].T, dtype=np.float64
        )

    def describe(self) -> str:
        """Convention-explicit prose: what was kept, and whether it is a projection."""

        zone_label = format_direction_indices(
            tuple(int(value) for value in self.zone_axis.indices), style="plain"
        )
        zolz = int(np.count_nonzero(self.zolz_mask))
        holz = int(np.count_nonzero(self.holz_mask))
        orders = sorted({int(value) for value in self.laue_zone if value != 0})
        holz_note = (
            "The set is confined to the zeroth Laue zone, so the calculation is a "
            "projection along the beam: it carries the symmetry of the projected "
            "potential, which is at least as high as the crystal's own and can be strictly "
            "higher. Add higher-order Laue zones before drawing any symmetry conclusion."
            if holz == 0
            else (
                f"Higher-order Laue zones {orders} contribute {holz} beams, so the "
                "projection symmetry is broken and the calculation can distinguish "
                "symmetries that a ZOLZ-only pattern cannot."
            )
        )
        return (
            f"Beam set of {self.size} reflections for {self.phase.name} down {zone_label} "
            f"at {self.beam_energy_kev:.0f} kV: the transmitted beam, {zolz - 1} further "
            f"zeroth-Laue-zone reflections, and {holz} higher-order ones. {holz_note}"
        )


def beam_set_for_zone(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    beam_energy_kev: float = 200.0,
    max_index: int = 4,
    g_max_inv_angstrom: float = 1.6,
    max_excitation_error_inv_angstrom: float = 0.05,
    convergence_semi_angle_mrad: float = 0.0,
    laue_zones: tuple[int, ...] = (0,),
    in_plane_rotation_deg: float = 0.0,
    apply_centering_absences: bool = True,
) -> BeamSet:
    """Choose the reflections a dynamical calculation will follow.

    What it does
        Enumerates ``hkl`` within ``max_index``, discards those forbidden by the
        lattice centering, keeps those in the requested Laue zones and inside
        the ``|g|`` cut-off, and then applies the selection that matters: a
        reflection is kept when it comes within
        ``max_excitation_error_inv_angstrom`` of the Bragg condition *for some
        incident direction in the illumination cone*,

        .. math::

           \\min_{|\\boldsymbol{\\theta}| \\le \\alpha} |s_{g}(\\boldsymbol{\\theta})|
             = \\max\\bigl(0,\\; |s_{g}(\\mathbf{0})| - \\alpha|\\mathbf{g}_{\\perp}|\\bigr).

    When to use it
        Before every call to :func:`solve_bloch_waves`. The two knobs that
        change the answer rather than the cost are
        ``max_excitation_error_inv_angstrom`` (how far off Bragg a beam may be
        and still couple) and ``laue_zones``.

    Why the cone enters the selection
        A higher-order Laue zone reflection is far from Bragg at the centre of
        the pattern and *exactly* at Bragg somewhere inside the bright-field
        disc — that locus is a HOLZ line. Selecting on the zero-tilt excitation
        error alone would discard every HOLZ reflection and, with it, every HOLZ
        line and the whole symmetry-breaking mechanism. Pass the convergence
        semi-angle actually used.

    Parameters
    ----------
    phase, zone_axis:
        The zone axis must belong to the phase, and the phase must carry a unit
        cell.
    beam_energy_kev:
        Accelerating voltage.
    max_index:
        Largest absolute Miller index enumerated.
    g_max_inv_angstrom:
        Radial cut-off on ``|g|``.
    max_excitation_error_inv_angstrom:
        Coupling window, applied to the minimum of ``|s_g|`` over the cone.
    convergence_semi_angle_mrad:
        Illumination semi-angle ``alpha``. Zero gives the parallel-beam
        selection used by `pytex.diffraction.kinematic`.
    laue_zones:
        Which Laue zones to admit; ``(0,)`` is a projection calculation,
        ``(0, 1)`` adds the first HOLZ.
    in_plane_rotation_deg:
        Rotation of the zone basis about the beam, as in
        `pytex.diffraction.kinematic.zone_basis_from_axis`.
    apply_centering_absences:
        Remove reflections forbidden by the lattice centering.

    Returns
    -------
    BeamSet
        With ``(0, 0, 0)`` first, then reflections sorted by increasing
        ``|g|`` so the ordering is deterministic and readable.

    Raises
    ------
    ValueError
        If the zone axis belongs to a different phase, if the phase has no unit
        cell, or if a Laue zone order is negative in a way that cannot occur.

    Notes
    -----
    **Algorithm.**

    1. Build the zone basis ``(u, v, z)`` with ``z`` along the zone axis.
    2. Enumerate the ``hkl`` grid, drop centering-forbidden reflections.
    3. Compute the Laue-zone order ``n = g . [uvw]_reduced`` and keep the
       requested zones.
    4. Keep ``|g| <= g_max``.
    5. Keep reflections whose minimum ``|s_g|`` over the cone is inside the
       window.
    6. Sort by ``|g|``, prepend the transmitted beam.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(
            "A dynamical calculation needs the atom positions, because the potential "
            f"coefficients do: phase '{phase.name}' carries no unit cell."
        )
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if g_max_inv_angstrom <= 0.0:
        raise ValueError("g_max_inv_angstrom must be strictly positive.")
    if max_excitation_error_inv_angstrom < 0.0:
        raise ValueError("max_excitation_error_inv_angstrom must be non-negative.")
    if convergence_semi_angle_mrad < 0.0:
        raise ValueError("convergence_semi_angle_mrad must be non-negative.")
    if not laue_zones:
        raise ValueError("laue_zones must name at least one Laue zone, normally 0.")

    wavelength = electron_wavelength_angstrom(beam_energy_kev)
    zone_unit = normalize_vector(np.asarray(zone_axis.unit_vector, dtype=np.float64))
    basis = zone_basis_from_axis(zone_unit, in_plane_rotation_deg=in_plane_rotation_deg)

    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid != 0, axis=1)]
    if apply_centering_absences:
        grid = grid[centering_allowed_mask(grid, ReflectionCondition.from_phase(phase))]

    axis_indices = np.asarray(zone_axis.indices, dtype=np.int64)
    non_zero = np.abs(axis_indices[axis_indices != 0])
    divisor = int(np.gcd.reduce(non_zero)) if non_zero.size else 1
    reduced_axis = axis_indices // max(divisor, 1)
    laue_order = grid @ reduced_axis
    grid = grid[np.isin(laue_order, np.asarray(laue_zones, dtype=np.int64))]
    laue_order = grid @ reduced_axis

    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_cartesian = grid.astype(np.float64) @ reciprocal.T
    g_zone = g_cartesian @ basis
    g_magnitude = np.linalg.norm(g_cartesian, axis=1)

    zero_tilt_s = g_zone[:, 2] - 0.5 * wavelength * g_magnitude * g_magnitude
    alpha = float(convergence_semi_angle_mrad) * 1e-3
    in_plane = np.linalg.norm(g_zone[:, :2], axis=1)
    reachable_s = np.maximum(np.abs(zero_tilt_s) - alpha * in_plane, 0.0)

    keep = (g_magnitude <= g_max_inv_angstrom) & (
        reachable_s <= max_excitation_error_inv_angstrom
    )
    grid = grid[keep]
    g_zone = g_zone[keep]
    g_magnitude = g_magnitude[keep]
    laue_order = laue_order[keep]

    order = np.lexsort((grid[:, 2], grid[:, 1], grid[:, 0], np.round(g_magnitude, 9)))
    return _assemble_beam_set(
        phase, zone_axis, float(beam_energy_kev), grid[order], basis, reduced_axis
    )


def beam_set_from_indices(
    phase: Phase,
    zone_axis: ZoneAxis,
    hkl: ArrayLike,
    *,
    beam_energy_kev: float = 200.0,
    in_plane_rotation_deg: float = 0.0,
) -> BeamSet:
    """Build a beam set from an explicit list of reflections.

    What it does
        Takes exactly the reflections named — no enumeration, no cut-off, no
        excitation-error window — and computes their zone-basis geometry. The
        transmitted beam is prepended if it is not already present, and
        duplicates are removed while preserving the order given.

    When to use it
        For the calculations whose whole point is a *chosen* beam set: a
        two-beam comparison against the closed form, a systematic row
        ``{g, 2g, 3g, ...}``, a convergence study that adds one reflection at a
        time, or reproducing a published calculation whose beam list is stated.
        For ordinary work use :func:`beam_set_for_zone`, which selects on the
        physics.

    Parameters
    ----------
    phase, zone_axis:
        The zone axis must belong to the phase, and the phase must carry a unit
        cell.
    hkl:
        ``(n, 3)`` integer Miller indices, or one triple. Reflections need not
        lie in the zeroth Laue zone and need not obey the zone law: their
        excitation errors will simply be large.
    beam_energy_kev, in_plane_rotation_deg:
        As in :func:`beam_set_for_zone`.

    Returns
    -------
    BeamSet

    Raises
    ------
    ValueError
        If the zone axis belongs to a different phase, or the phase has no unit
        cell.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(
            "A dynamical calculation needs the atom positions, because the potential "
            f"coefficients do: phase '{phase.name}' carries no unit cell."
        )

    indices = np.atleast_2d(np.asarray(hkl, dtype=np.int64))
    if indices.ndim != 2 or indices.shape[1] != 3:
        raise ValueError("hkl must have shape (3,) or (n, 3).")
    indices = indices[np.any(indices != 0, axis=1)]
    _, first_seen = np.unique(indices, axis=0, return_index=True)
    indices = indices[np.sort(first_seen)]

    zone_unit = normalize_vector(np.asarray(zone_axis.unit_vector, dtype=np.float64))
    basis = zone_basis_from_axis(zone_unit, in_plane_rotation_deg=in_plane_rotation_deg)
    axis_indices = np.asarray(zone_axis.indices, dtype=np.int64)
    non_zero = np.abs(axis_indices[axis_indices != 0])
    divisor = int(np.gcd.reduce(non_zero)) if non_zero.size else 1
    reduced_axis = axis_indices // max(divisor, 1)
    return _assemble_beam_set(
        phase, zone_axis, float(beam_energy_kev), indices, basis, reduced_axis
    )


def _assemble_beam_set(
    phase: Phase,
    zone_axis: ZoneAxis,
    beam_energy_kev: float,
    grid: np.ndarray,
    basis: np.ndarray,
    reduced_axis: np.ndarray,
) -> BeamSet:
    """Shared geometry assembly for the two public beam-set constructors."""

    reciprocal = as_float_array(phase.lattice.reciprocal_basis().matrix, shape=(3, 3))
    g_cartesian = grid.astype(np.float64) @ reciprocal.T
    g_zone = g_cartesian @ basis
    g_magnitude = np.linalg.norm(g_cartesian, axis=1)
    laue_order = grid @ reduced_axis
    return BeamSet(
        phase=phase,
        zone_axis=zone_axis,
        beam_energy_kev=beam_energy_kev,
        miller_indices=np.concatenate([np.zeros((1, 3), dtype=np.int64), grid], axis=0),
        g_zone=np.concatenate([np.zeros((1, 3), dtype=np.float64), g_zone], axis=0),
        g_magnitude_inv_angstrom=np.concatenate(
            [np.zeros(1, dtype=np.float64), g_magnitude], axis=0
        ),
        laue_zone=np.concatenate([np.zeros(1, dtype=np.int64), laue_order], axis=0),
        zone_basis_crystal=basis,
    )


# --------------------------------------------------------------------------- #
# The structure matrix
# --------------------------------------------------------------------------- #


def structure_matrix(
    beams: BeamSet,
    tilts_rad: ArrayLike,
    *,
    absorption: AbsorptionModel | None = None,
) -> np.ndarray:
    """Build the dynamical structure matrix ``A`` at each incident tilt.

    What it does
        Assembles

        .. math::

           A_{gg} = 2 s_{g}(\\boldsymbol{\\theta}) + i\\,r_{0}\\nu_{0},
           \\qquad
           A_{gh} = \\nu_{g-h} + i\\,r_{g}\\nu_{g-h}\\;(g \\neq h),

        where :math:`\\nu` are the potential coefficients of
        :func:`potential_coefficients_inv_angstrom` and :math:`r_{0}, r_{g}` the
        ratios of :class:`AbsorptionModel`. Every distinct difference vector
        ``g - h`` is evaluated once and scattered back, so the cost is set by
        the number of distinct differences rather than by ``n^2``.

    When to use it
        Directly, to inspect the coupling — the off-diagonal magnitudes are the
        reciprocal extinction distances and say at a glance which beams matter.
        Otherwise :func:`solve_bloch_waves` calls it for you.

    Parameters
    ----------
    beams:
        The beam set; its ordering fixes the matrix ordering.
    tilts_rad:
        ``(m, 2)`` incident tilts in radians, or one pair.
    absorption:
        Defaults to :meth:`AbsorptionModel.none`, giving a Hermitian matrix.

    Returns
    -------
    np.ndarray
        ``(m, n, n)`` complex matrices in 1/angstrom.

    Notes
    -----
    The mean inner potential contributes a *real* constant :math:`\\nu_{0}` to
    every diagonal element. It is deliberately omitted: a constant added to the
    diagonal multiplies the whole solution by one phase factor and changes no
    intensity, and carrying it would only invite the reader to attach meaning to
    an unobservable. Its *imaginary* partner :math:`r_{0}\\nu_{0}` is kept,
    because that one is observable — it is the normal absorption.
    """

    model = absorption or AbsorptionModel.none()
    indices = beams.miller_indices
    size = beams.size

    differences = indices[:, None, :] - indices[None, :, :]
    flat = differences.reshape(-1, 3)
    unique, inverse = np.unique(flat, axis=0, return_inverse=True)
    unique_coefficients = potential_coefficients_inv_angstrom(
        beams.phase, unique, beam_energy_kev=beams.beam_energy_kev
    )
    elastic = unique_coefficients[np.reshape(inverse, -1)].reshape(size, size)

    zero_coefficient = potential_coefficients_inv_angstrom(
        beams.phase, [[0, 0, 0]], beam_energy_kev=beams.beam_energy_kev
    )[0]

    matrix = elastic * (1.0 + 1.0j * model.reflection_ratio)
    diagonal_absorption = 1.0j * model.mean_ratio * float(np.real(zero_coefficient))
    np.einsum("ii->i", matrix)[...] = diagonal_absorption

    excitation = beams.excitation_errors(tilts_rad)
    stacked = np.broadcast_to(matrix, (excitation.shape[0], size, size)).copy()
    diagonal_view = np.einsum("...ii->...i", stacked)
    diagonal_view += 2.0 * excitation
    return np.asarray(stacked, dtype=np.complex128)


# --------------------------------------------------------------------------- #
# The Bloch-wave solution
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class BlochWaveSolution:
    """Beam intensities from a many-beam calculation, with the wavefield kept.

    Purpose
    -------
    The output of :func:`solve_bloch_waves`. It reports the intensity of every
    beam at every requested incident tilt, and — when asked — retains the
    eigenbasis, so a thickness series costs one extra matrix product rather than
    a new eigen-decomposition, and so the dispersion surface and the per-branch
    absorption can be inspected directly.

    Attributes
    ----------
    beams : BeamSet
    tilts_rad : np.ndarray
        ``(m, 2)`` incident tilts.
    thickness_angstrom : float
    absorption : AbsorptionModel
    intensities : np.ndarray
        ``(m, n)`` beam intensities at ``thickness_angstrom``, in the beam
        ordering of ``beams``. Column ``0`` is the transmitted beam.
    eigenvalues : np.ndarray
        ``(m, n)`` complex Bloch-wave excitations :math:`\\gamma_{j}`. The real
        parts are the dispersion surface; the imaginary parts are the
        absorption coefficient of each branch, and their spread *is* the
        anomalous absorption.
    eigenvectors : np.ndarray or None
        ``(m, n, n)`` Bloch-wave coefficients, columns indexed by branch, or
        ``None`` when ``keep_eigenbasis`` was false.
    excitation_amplitudes : np.ndarray or None
        ``(m, n)`` amplitudes with which the incident plane wave excites each
        branch, or ``None``.
    provenance : ProvenanceRecord or None
    """

    beams: BeamSet
    tilts_rad: np.ndarray
    thickness_angstrom: float
    absorption: AbsorptionModel
    intensities: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray | None = None
    excitation_amplitudes: np.ndarray | None = None
    provenance: ProvenanceRecord | None = None

    @property
    def transmitted_intensity(self) -> np.ndarray:
        """``(m,)`` intensity of the direct beam: the bright-field signal."""

        return np.asarray(self.intensities[:, 0], dtype=np.float64)

    @property
    def total_intensity(self) -> np.ndarray:
        """``(m,)`` sum over beams.

        Exactly ``1`` at every tilt and thickness when absorption is off — the
        unitarity check on the whole calculation. Strictly below ``1`` when it
        is on, and the shortfall is the absorbed fraction.
        """

        return np.asarray(np.sum(self.intensities, axis=1), dtype=np.float64)

    @property
    def normal_absorption_factor(self) -> float:
        """The scalar :math:`\\exp(-2\\pi t/\\xi'_{0})` contained in every intensity.

        Purpose
        -------
        Makes the separation explicit and checkable. Dividing
        :attr:`intensities` by this number removes normal absorption entirely
        and leaves only the anomalous part, which is the part that carries
        physics. Any conclusion that survives this division does not depend on
        the phenomenological ``mean_ratio`` of :class:`AbsorptionModel`.
        """

        zero_coefficient = potential_coefficients_inv_angstrom(
            self.beams.phase, [[0, 0, 0]], beam_energy_kev=self.beams.beam_energy_kev
        )[0]
        exponent = (
            2.0
            * math.pi
            * self.thickness_angstrom
            * self.absorption.mean_ratio
            * float(np.real(zero_coefficient))
        )
        return float(math.exp(-exponent))

    def intensity_of(self, hkl: ArrayLike) -> np.ndarray:
        """``(m,)`` intensity of one named reflection across the tilts."""

        return np.asarray(
            self.intensities[:, self.beams.index_of(hkl)], dtype=np.float64
        )

    def intensities_at(self, thickness_angstrom: ArrayLike) -> np.ndarray:
        """Re-evaluate intensities at other thicknesses from the stored eigenbasis.

        Purpose
        -------
        The expensive step is the eigen-decomposition, and it does not depend on
        thickness. A thickness series — a wedge profile, a fringe count, a
        convergence check — therefore costs one matrix product per thickness.

        Parameters
        ----------
        thickness_angstrom:
            Scalar or ``(k,)`` thicknesses in angstrom, strictly positive.

        Returns
        -------
        np.ndarray
            ``(k, m, n)`` intensities, or ``(m, n)`` for a scalar input.

        Raises
        ------
        ValueError
            If the eigenbasis was not kept; pass ``keep_eigenbasis=True`` to
            :func:`solve_bloch_waves`.
        """

        if self.eigenvectors is None or self.excitation_amplitudes is None:
            raise ValueError(
                "The eigenbasis was not retained, so this solution can only report the "
                "thickness it was solved at. Call solve_bloch_waves(..., "
                "keep_eigenbasis=True) to enable thickness series."
            )
        thicknesses = np.atleast_1d(np.asarray(thickness_angstrom, dtype=np.float64))
        if np.any(~np.isfinite(thicknesses)) or np.any(thicknesses <= 0.0):
            raise ValueError("thickness_angstrom must be finite and strictly positive.")
        phases = np.exp(1.0j * np.pi * thicknesses[:, None, None] * self.eigenvalues[None])
        amplitudes = np.einsum(
            "mgj,tmj->tmg", self.eigenvectors, phases * self.excitation_amplitudes[None]
        )
        intensities = np.abs(amplitudes) ** 2
        if np.isscalar(thickness_angstrom) or np.asarray(thickness_angstrom).ndim == 0:
            return np.asarray(intensities[0], dtype=np.float64)
        return np.asarray(intensities, dtype=np.float64)

    def describe(self) -> str:
        """Convention-explicit prose: the calculation, its checks, and its limits."""

        zone_label = format_direction_indices(
            tuple(int(value) for value in self.beams.zone_axis.indices), style="plain"
        )
        totals = self.total_intensity
        strongest = int(np.argmax(np.mean(self.intensities[:, 1:], axis=0))) + 1
        strongest_label = format_plane_indices(
            tuple(int(value) for value in self.beams.miller_indices[strongest]),
            style="plain",
        )
        holz = int(np.count_nonzero(self.beams.holz_mask))
        conservation = (
            "Absorption is off, so the sum over beams is unity to "
            f"{float(np.max(np.abs(totals - 1.0))):.2e} at every tilt — the unitarity check "
            "that the coupled solution is right."
            if not self.absorption.is_absorbing
            else (
                f"With absorption the beams carry {float(np.min(totals)):.4f} to "
                f"{float(np.max(totals)):.4f} of the incident intensity; the remainder is "
                f"absorbed, of which the scalar normal part is a factor "
                f"{self.normal_absorption_factor:.4f}."
            )
        )
        return (
            f"Many-beam Bloch-wave solution for {self.beams.phase.name} down {zone_label} "
            f"at {self.beams.beam_energy_kev:.0f} kV, {self.beams.size} coupled beams "
            f"({holz} of them higher-order Laue zone), foil thickness "
            f"{self.thickness_angstrom:.0f} A, evaluated at {self.tilts_rad.shape[0]} "
            f"incident directions. The strongest diffracted beam is {strongest_label}. "
            f"{conservation} Each beam here is coupled to every other through the full "
            "structure matrix, so unlike a two-beam disc the intensities of different "
            "reflections in one pattern are mutually consistent."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": BLOCH_WAVE_SOLUTION_SCHEMA,
            "phase": self.beams.phase.name,
            "zone_axis": [int(value) for value in self.beams.zone_axis.indices],
            "beam_energy_kev": self.beams.beam_energy_kev,
            "thickness_angstrom": self.thickness_angstrom,
            "beam_count": self.beams.size,
            "holz_beam_count": int(np.count_nonzero(self.beams.holz_mask)),
            "absorption": {
                "mean_ratio": self.absorption.mean_ratio,
                "reflection_ratio": self.absorption.reflection_ratio,
            },
            "normal_absorption_factor": self.normal_absorption_factor,
            "tilt_count": int(self.tilts_rad.shape[0]),
            "beams": [
                {
                    "hkl": [int(value) for value in row],
                    "laue_zone": int(zone),
                    "mean_intensity": float(np.mean(column)),
                }
                for row, zone, column in zip(
                    self.beams.miller_indices,
                    self.beams.laue_zone,
                    self.intensities.T,
                    strict=True,
                )
            ],
        }


def solve_bloch_waves(
    beams: BeamSet,
    tilts_rad: ArrayLike,
    *,
    thickness_angstrom: float,
    absorption: AbsorptionModel | None = None,
    keep_eigenbasis: bool = False,
    provenance: ProvenanceRecord | None = None,
) -> BlochWaveSolution:
    """Solve the coupled many-beam problem at a set of incident directions.

    What it does
        Diagonalizes the structure matrix and propagates the incident plane wave
        through the foil:

        .. math::

           \\mathbf{A} = \\mathbf{C}\\,\\mathrm{diag}(\\gamma_{j})\\,\\mathbf{C}^{-1},
           \\qquad
           \\psi_{g}(t) = \\sum_{j} C_{gj}\\,\\alpha_{j}\\,e^{i\\pi\\gamma_{j}t},
           \\qquad
           \\boldsymbol{\\alpha} = \\mathbf{C}^{-1}\\mathbf{e}_{0},

        and reports :math:`I_{g} = |\\psi_{g}|^{2}`. The Bloch waves are the
        eigenvectors: each is a wavefield that propagates through the crystal
        unchanged in shape, attenuated at its own rate
        :math:`\\mathrm{Im}\\,\\gamma_{j}`.

    When to use it
        Whenever the two-beam approximation is not defensible, which at a zone
        axis is always: to compute CBED disc intensities that are consistent
        across the pattern, to produce HOLZ deficiency lines, to model
        thickness fringes with realistic decay, and as the forward model behind
        `pytex.diffraction.diffraction_groups`.

    Parameters
    ----------
    beams:
        From :func:`beam_set_for_zone`. Its ordering is the ordering of every
        returned array.
    tilts_rad:
        ``(m, 2)`` incident tilts in radians, or one pair. For a CBED disc these
        sample the convergence cone; for a rocking curve they run along a line.
    thickness_angstrom:
        Foil thickness, strictly positive.
    absorption:
        Defaults to :meth:`AbsorptionModel.none`.
    keep_eigenbasis:
        Retain the eigenvectors and excitation amplitudes so
        :meth:`BlochWaveSolution.intensities_at` can serve a thickness series.
        Costs ``m * n^2`` complex numbers, which is why it is off by default.
    provenance:
        Optional record.

    Returns
    -------
    BlochWaveSolution

    Raises
    ------
    ValueError
        If the thickness is not positive, or the tilt array is malformed.

    Notes
    -----
    **Algorithm.**

    1. Build ``A`` at every tilt (:func:`structure_matrix`).
    2. Eigen-decompose. ``A`` is complex symmetric for a centrosymmetric
       structure and merely complex otherwise, so the general non-Hermitian
       routine is used; the eigenvectors are not orthogonal in general and the
       excitation amplitudes must be obtained by solving, not by projection.
       Using an orthogonality that does not hold is the classic silent error in
       a Bloch-wave implementation, and it produces plausible-looking rocking
       curves with the wrong contrast.
    3. Solve ``C alpha = e_0`` for the excitation amplitudes.
    4. Recombine at the requested thickness.

    **Cost.** One eigen-decomposition per tilt, so ``O(m n^3)``. Halving the
    excitation-error window is usually a better economy than reducing the tilt
    sampling, because it removes beams that were barely coupled.
    """

    if not np.isfinite(thickness_angstrom) or thickness_angstrom <= 0.0:
        raise ValueError("thickness_angstrom must be finite and strictly positive.")
    model = absorption or AbsorptionModel.none()
    tilts = np.atleast_2d(as_float_array(tilts_rad, shape=(None, 2)))

    matrices = structure_matrix(beams, tilts, absorption=model)
    eigenvalues, eigenvectors = np.linalg.eig(matrices)

    incident = np.zeros((tilts.shape[0], beams.size), dtype=np.complex128)
    incident[:, 0] = 1.0
    amplitudes = np.linalg.solve(eigenvectors, incident[..., None])[..., 0]

    propagated = amplitudes * np.exp(1.0j * np.pi * thickness_angstrom * eigenvalues)
    wavefield = np.einsum("mgj,mj->mg", eigenvectors, propagated)
    intensities = np.asarray(np.abs(wavefield) ** 2, dtype=np.float64)

    return BlochWaveSolution(
        beams=beams,
        tilts_rad=tilts,
        thickness_angstrom=float(thickness_angstrom),
        absorption=model,
        intensities=intensities,
        eigenvalues=np.asarray(eigenvalues, dtype=np.complex128),
        eigenvectors=(
            np.asarray(eigenvectors, dtype=np.complex128) if keep_eigenbasis else None
        ),
        excitation_amplitudes=(
            np.asarray(amplitudes, dtype=np.complex128) if keep_eigenbasis else None
        ),
        provenance=provenance,
    )
