"""HOLZ lines: the sharp features that measure a lattice parameter.

A zone-axis pattern is blind to the one lattice dimension along the beam, and
nearly blind to small changes in the other two: a spot moves by the same
fractional amount the lattice does, which for a strain of :math:`10^{-4}` is
invisible. Higher-order Laue zone lines are the way out. They are the loci, in
the illumination cone, where a HOLZ reflection is *exactly* at the Bragg
condition — dark **deficiency** lines in the bright-field disc, matching bright
**excess** lines in the corresponding diffracted disc — and because a HOLZ
reflection has a large :math:`|\\mathbf{g}|`, its Bragg condition is met over a
very narrow range of incident direction. The line is therefore sharp, and a
small change in the lattice moves it by much more than its own width.

The geometry is exact
---------------------

For a reflection :math:`\\mathbf{g}` the excitation error is affine in the
incident tilt,

.. math::

   s_{g}(\\boldsymbol{\\theta}) = g_{z}
       - \\boldsymbol{\\theta}\\cdot\\mathbf{g}_{\\perp}
       - \\tfrac{1}{2}\\lambda|\\mathbf{g}|^{2},

so :math:`s_{g} = 0` is a **straight line** in the plane of incident directions,

.. math::

   \\boldsymbol{\\theta}\\cdot\\hat{\\mathbf{g}}_{\\perp} = d_{g},
   \\qquad
   d_{g} = \\frac{g_{z} - \\tfrac{1}{2}\\lambda|\\mathbf{g}|^{2}}
                 {|\\mathbf{g}_{\\perp}|},

with unit normal :math:`\\hat{\\mathbf{g}}_{\\perp}` and signed distance
:math:`d_{g}` from the pattern centre. No approximation beyond the one already
in :math:`s_{g}` itself enters, so the line positions this module reports are as
accurate as the excitation error is. The line crosses the bright-field disc when
:math:`|d_{g}| \\le \\alpha`.

The metrology, and the trap in it
---------------------------------

Scale the lattice by :math:`1 + \\varepsilon` and every :math:`\\mathbf{g}`
shrinks by the same factor. Substituting gives the line offset in closed form,

.. math::

   d_{g}(\\varepsilon, \\lambda) =
     \\frac{g_{z}}{|\\mathbf{g}_{\\perp}|}
     - \\frac{\\lambda|\\mathbf{g}|^{2}}{2(1+\\varepsilon)|\\mathbf{g}_{\\perp}|},
   \\qquad
   \\frac{\\partial d_{g}}{\\partial\\varepsilon}\\Bigg|_{0}
     = \\frac{\\lambda|\\mathbf{g}|^{2}}{2|\\mathbf{g}_{\\perp}|},

and the wavelength enters the *same* term with the opposite sign. Setting
:math:`\\lambda \\to \\lambda(1+\\varepsilon)` therefore cancels a lattice strain
:math:`\\varepsilon` **exactly**, at every reflection simultaneously:

    A fractional change in lattice parameter and a fractional change in
    wavelength are indistinguishable from HOLZ line positions.

This is not a limitation of the model; it is why quantitative HOLZ metrology
begins by calibrating the accelerating voltage against a standard of known
lattice parameter, and why a measured "lattice parameter" from an uncalibrated
microscope is a measurement of the high-tension supply. The degeneracy is
asserted in the tests rather than described.

What is here and what is not
----------------------------

This module is **geometry**: line positions, their chords inside the discs,
their angular width for a given foil thickness, their intersections, and their
sensitivity to strain and to voltage. It says nothing about how dark a line is.
Line *contrast* requires the coupled dynamical solution with those HOLZ
reflections in the beam set, which is
`pytex.diffraction.dynamical.solve_bloch_waves`; the two are joined in
`pytex.diffraction.cbed.simulate_cbed_pattern`.

See Also
--------
`pytex.diffraction.cbed.holz_ring_radii_inv_angstrom` : the rings, which measure
    the layer spacing; these are the lines, which measure the lattice.
`docs/tex/algorithms/dynamical_cbed_and_symmetry_determination.tex` : derivations.
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

__all__ = [
    "HOLZ_LINE_PATTERN_SCHEMA",
    "HOLZLine",
    "HOLZLineIntersection",
    "HOLZLinePattern",
    "holz_line_pattern",
]

#: Schema identifier of the HOLZ-line-pattern payload.
HOLZ_LINE_PATTERN_SCHEMA = "pytex.holz_line_pattern/1"

_PARALLEL_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class HOLZLine:
    """One higher-order Laue zone line, as an exact locus in the tilt plane.

    Purpose
    -------
    The Bragg locus of a single HOLZ reflection, carried in the coordinates an
    experiment actually uses: incident tilt in radians, and detector position in
    millimetres. Because the locus is a straight line in closed form, everything
    about it — where it crosses the bright-field disc, how wide it is, how far it
    moves for a given strain — is available analytically rather than by sampling
    a simulated pattern.

    Attributes
    ----------
    miller_indices : np.ndarray
        The HOLZ reflection ``(hkl)``.
    laue_zone : int
        Its Laue-zone order ``n = g . [uvw]``, never zero.
    normal_tilt : np.ndarray
        ``(2,)`` unit normal of the line in the tilt plane, which is
        ``g_perp / |g_perp|`` in the zone basis. The line runs perpendicular to
        ``g``, as every Bragg locus does.
    offset_rad : float
        Signed distance ``d_g`` of the line from the pattern centre, in radians
        of incident tilt. This is the number a metrology measurement reads.
    g_perp_inv_angstrom : float
        ``|g_perp|``, the in-plane length of ``g``. It sets both the sharpness of
        the line and the lever arm of the measurement.
    g_magnitude_inv_angstrom : float
        ``|g|``.
    g_zone_inv_angstrom : float
        ``g_z``, the component along the beam.
    wavelength_angstrom : float
    convergence_semi_angle_rad : float
        The illumination cone half-angle the line is reported against.
    camera_constant_mm_angstrom : float
    disc_centre_mm : np.ndarray
        ``(2,)`` detector position of this reflection's disc, where the *excess*
        line appears. The *deficiency* line appears at the same tilt coordinates
        in the bright-field disc, that is, translated to the origin.
    """

    miller_indices: np.ndarray
    laue_zone: int
    normal_tilt: np.ndarray
    offset_rad: float
    g_perp_inv_angstrom: float
    g_magnitude_inv_angstrom: float
    g_zone_inv_angstrom: float
    wavelength_angstrom: float
    convergence_semi_angle_rad: float
    camera_constant_mm_angstrom: float
    disc_centre_mm: np.ndarray

    @property
    def label(self) -> str:
        """Rendered index label, per the repository notation standard."""

        return format_plane_indices(
            tuple(int(value) for value in self.miller_indices), style="mathtext"
        )

    @property
    def crosses_bright_field(self) -> bool:
        """Whether the line falls inside the illumination cone.

        A HOLZ reflection whose Bragg locus lies outside the cone contributes no
        visible line, however strongly it would scatter if it were reached.
        """

        return abs(self.offset_rad) <= self.convergence_semi_angle_rad

    @property
    def strain_sensitivity_rad(self) -> float:
        """``d(offset)/d(strain)`` at zero strain: ``lambda |g|^2 / (2 |g_perp|)``.

        Purpose
        -------
        The lever arm of the measurement, in radians of line movement per unit
        of isotropic lattice strain. Comparing it with
        :meth:`angular_width_rad` says directly what strain resolution a given
        foil thickness supports, which is the question an experiment design has
        to answer before the microscope is switched on.
        """

        return (
            self.wavelength_angstrom
            * self.g_magnitude_inv_angstrom**2
            / (2.0 * self.g_perp_inv_angstrom)
        )

    def offset_at(
        self,
        *,
        lattice_strain: float = 0.0,
        wavelength_angstrom: float | None = None,
    ) -> float:
        """The line offset for a strained lattice, a different voltage, or both.

        What it does
            Evaluates the closed form

            .. math::

               d_{g}(\\varepsilon, \\lambda) =
                 \\frac{g_{z}}{|\\mathbf{g}_{\\perp}|}
                 - \\frac{\\lambda|\\mathbf{g}|^{2}}
                        {2(1+\\varepsilon)|\\mathbf{g}_{\\perp}|},

            exactly rather than to first order, so it stays correct for the
            large excursions used to *illustrate* the effect as well as for the
            small ones used to measure it.

        When to use it
            To predict how far a line moves for a candidate strain, to compare a
            measured shift against a model, and — most usefully — to demonstrate
            the wavelength/strain degeneracy that makes voltage calibration
            mandatory.

        Parameters
        ----------
        lattice_strain:
            Isotropic fractional change ``eps`` in the lattice parameters, so
            every spacing becomes ``(1 + eps)`` times its original value. Must be
            greater than ``-1``.
        wavelength_angstrom:
            Electron wavelength; defaults to the one the line was built with.

        Returns
        -------
        float
            The offset in radians of incident tilt.

        Notes
        -----
        Passing ``wavelength_angstrom = lambda * (1 + eps)`` together with
        ``lattice_strain = eps`` returns the unstrained offset to machine
        precision, for every reflection at once. That is the degeneracy stated
        in the module docstring, and it is why this method takes both arguments
        rather than pretending they are independent measurements.
        """

        if not np.isfinite(lattice_strain) or lattice_strain <= -1.0:
            raise ValueError(
                "lattice_strain must be finite and greater than -1; a strain of -1 or "
                "below would collapse the lattice."
            )
        wavelength = (
            self.wavelength_angstrom if wavelength_angstrom is None else float(wavelength_angstrom)
        )
        if not np.isfinite(wavelength) or wavelength <= 0.0:
            raise ValueError("wavelength_angstrom must be finite and strictly positive.")
        return float(
            self.g_zone_inv_angstrom / self.g_perp_inv_angstrom
            - wavelength
            * self.g_magnitude_inv_angstrom**2
            / (2.0 * (1.0 + lattice_strain) * self.g_perp_inv_angstrom)
        )

    def angular_width_rad(self, thickness_angstrom: float) -> float:
        """Half-width of the line, from the thickness broadening of the rocking curve.

        What it does
            The two-beam rocking curve has its first zero at
            :math:`|s| = 1/t`, and :math:`s` varies across the disc at the rate
            :math:`|\\mathbf{g}_{\\perp}|` per radian of tilt, so the line's
            angular half-width is

            .. math::

               \\Delta\\theta = \\frac{1}{t\\,|\\mathbf{g}_{\\perp}|}.

        When to use it
            To decide whether an experiment can resolve the shift it is trying to
            measure. The ratio :meth:`strain_sensitivity_rad` divided by this
            width is the strain that moves the line by its own width, and it
            improves in proportion to the foil thickness — which is why HOLZ
            metrology wants a *thick* specimen, the opposite of most TEM work.

        Parameters
        ----------
        thickness_angstrom:
            Foil thickness, strictly positive.

        Returns
        -------
        float
            Half-width in radians of incident tilt.
        """

        if not np.isfinite(thickness_angstrom) or thickness_angstrom <= 0.0:
            raise ValueError("thickness_angstrom must be finite and strictly positive.")
        return float(1.0 / (thickness_angstrom * self.g_perp_inv_angstrom))

    def resolvable_strain(self, thickness_angstrom: float) -> float:
        """The strain that shifts this line by one half-width.

        A single-number figure of merit for the line, combining
        :meth:`strain_sensitivity_rad` and :meth:`angular_width_rad`. It scales
        as :math:`1/(t\\,\\lambda|\\mathbf{g}|^{2})`, so a thicker foil and a
        higher-order reflection both help, and it is the reason HOLZ lines rather
        than ZOLZ spots are used for this measurement at all.
        """

        return float(self.angular_width_rad(thickness_angstrom) / self.strain_sensitivity_rad)

    def chord_tilt_rad(self) -> np.ndarray | None:
        """The two endpoints of the line inside the illumination cone.

        Returns
        -------
        np.ndarray or None
            ``(2, 2)`` tilt coordinates of the chord endpoints, or ``None`` when
            the line misses the cone entirely.
        """

        alpha = self.convergence_semi_angle_rad
        if abs(self.offset_rad) > alpha:
            return None
        half_chord = math.sqrt(max(alpha * alpha - self.offset_rad * self.offset_rad, 0.0))
        normal = np.asarray(self.normal_tilt, dtype=np.float64)
        tangent = np.array([-normal[1], normal[0]], dtype=np.float64)
        centre = self.offset_rad * normal
        return np.stack([centre - half_chord * tangent, centre + half_chord * tangent])

    def deficiency_chord_mm(self) -> np.ndarray | None:
        """The chord of the dark line inside the **bright-field** disc, in millimetres.

        The bright-field disc is centred on the optic axis, so the detector
        coordinate of an incident tilt is ``(C / lambda) * theta`` with ``C`` the
        camera constant. Returns ``None`` when the line misses the disc.
        """

        chord = self.chord_tilt_rad()
        if chord is None:
            return None
        return np.asarray(
            chord * (self.camera_constant_mm_angstrom / self.wavelength_angstrom),
            dtype=np.float64,
        )

    def excess_chord_mm(self) -> np.ndarray | None:
        """The chord of the bright line inside this reflection's **own** disc.

        The excess line sits at the same incident tilts as the deficiency line —
        it is the same Bragg condition — so on the detector it is the deficiency
        chord translated by the disc centre. Seeing the pair, dark in the direct
        disc and bright in the diffracted one, is how a HOLZ line is recognised
        as such rather than as a bend contour.
        """

        chord = self.deficiency_chord_mm()
        if chord is None:
            return None
        return np.asarray(chord + self.disc_centre_mm[None, :], dtype=np.float64)

    def describe(self) -> str:
        """Convention-explicit prose: what the line is and what it can measure."""

        inside = (
            f"crosses the bright-field disc at {self.offset_rad * 1e3:.3f} mrad from the "
            "centre"
            if self.crosses_bright_field
            else (
                f"lies {self.offset_rad * 1e3:.3f} mrad from the centre, outside the "
                f"{self.convergence_semi_angle_rad * 1e3:.2f} mrad cone, so it produces no "
                "visible line at this convergence angle"
            )
        )
        indices = format_plane_indices(
            tuple(int(value) for value in self.miller_indices), style="plain"
        )
        return (
            f"HOLZ line of {indices} "
            f"in Laue zone n = {self.laue_zone}: normal to g, {inside}. Its offset moves by "
            f"{self.strain_sensitivity_rad * 1e3:.3f} mrad per unit of isotropic lattice "
            "strain, and by an exactly compensating amount for a fractional change in "
            "wavelength, so the line position alone cannot separate a strained lattice from "
            "a mis-set accelerating voltage."
        )


@dataclass(frozen=True, slots=True)
class HOLZLineIntersection:
    """Where two HOLZ lines cross, and how fast that crossing moves.

    Purpose
    -------
    Line intersections, not individual lines, are what a HOLZ measurement
    normally reads: a crossing is a point rather than a locus, and when the two
    lines are nearly parallel it moves much faster than either line does. The
    amplification is :math:`1/|\\sin\\phi|` for an intersection angle
    :math:`\\phi`, and it is reported here so that a near-degenerate pair can be
    recognised as the sensitive one rather than as a nuisance.

    Attributes
    ----------
    first_indices, second_indices : np.ndarray
        The two reflections.
    position_tilt_rad : np.ndarray
        ``(2,)`` intersection point in the tilt plane.
    angle_deg : float
        The angle between the two lines, in ``(0, 90]`` degrees.
    strain_sensitivity_rad : float
        Speed of the intersection point under isotropic lattice strain, in
        radians per unit strain. Compare with the individual
        :attr:`HOLZLine.strain_sensitivity_rad` values to see the amplification.
    inside_bright_field : bool
    """

    first_indices: np.ndarray
    second_indices: np.ndarray
    position_tilt_rad: np.ndarray
    angle_deg: float
    strain_sensitivity_rad: float
    inside_bright_field: bool


@dataclass(frozen=True, slots=True)
class HOLZLinePattern:
    """Every HOLZ line of one zone-axis setting, with its metrology properties.

    Purpose
    -------
    What the sharp features of a convergent-beam pattern would look like, and
    what they can measure. The object is purely geometric — it has no
    intensities — which is deliberate: line *positions* are exact and are what
    metrology uses, while line *contrast* needs the coupled dynamical
    calculation and carries all of that calculation's approximations.

    Attributes
    ----------
    phase : Phase
    zone_axis : ZoneAxis
    beam_energy_kev : float
    convergence_semi_angle_mrad : float
    camera_constant_mm_angstrom : float
    lines : tuple of HOLZLine
        Sorted by increasing ``|offset|``, so the lines nearest the pattern
        centre — the ones an operator sees first — come first.
    zone_basis_crystal : np.ndarray
        ``(3, 3)`` with columns ``u``, ``v`` and the zone-axis unit vector.
    provenance : ProvenanceRecord or None
    """

    phase: Phase
    zone_axis: ZoneAxis
    beam_energy_kev: float
    convergence_semi_angle_mrad: float
    camera_constant_mm_angstrom: float
    lines: tuple[HOLZLine, ...]
    zone_basis_crystal: np.ndarray
    provenance: ProvenanceRecord | None = None

    @property
    def wavelength_angstrom(self) -> float:
        """Relativistic electron wavelength at this accelerating voltage."""

        return electron_wavelength_angstrom(self.beam_energy_kev)

    @property
    def bright_field_lines(self) -> tuple[HOLZLine, ...]:
        """The lines that actually cross the bright-field disc."""

        return tuple(line for line in self.lines if line.crosses_bright_field)

    @property
    def miller_indices(self) -> np.ndarray:
        """``(n, 3)`` reflections producing the lines, in the same order.

        These are exactly the HOLZ reflections a dynamical calculation must
        include for the lines to appear in a simulated disc; pass them to
        `pytex.diffraction.dynamical.beam_set_from_indices`, or use the same
        convergence semi-angle with
        `pytex.diffraction.dynamical.beam_set_for_zone`, which selects on the
        same criterion.
        """

        if not self.lines:
            return np.zeros((0, 3), dtype=np.int64)
        return np.stack([line.miller_indices for line in self.lines])

    def line_for(self, hkl: ArrayLike) -> HOLZLine:
        """The line of a named reflection.

        Raises
        ------
        KeyError
            If that reflection produces no line in this pattern.
        """

        wanted = np.asarray(hkl, dtype=np.int64).reshape(3)
        for line in self.lines:
            if np.array_equal(line.miller_indices, wanted):
                return line
        raise KeyError(
            f"Reflection {tuple(int(value) for value in wanted)} has no HOLZ line here. It "
            "may be in the zeroth Laue zone, forbidden by the lattice centering, or outside "
            "the search bounds."
        )

    def intersections(
        self, *, bright_field_only: bool = True, minimum_angle_deg: float = 1.0
    ) -> tuple[HOLZLineIntersection, ...]:
        """Every crossing of two lines, with its strain sensitivity.

        What it does
            Solves each pair of line equations. Two lines with unit normals
            :math:`\\hat{n}_{1}, \\hat{n}_{2}` and offsets :math:`d_{1}, d_{2}`
            meet where both constraints hold, and differentiating that
            :math:`2\\times 2` system with respect to strain gives the speed of
            the crossing directly.

        When to use it
            When designing or interpreting a HOLZ measurement. A pair meeting at
            a small angle moves as :math:`1/\\sin\\phi` times faster than the
            lines themselves, which is where the technique's remarkable
            sensitivity comes from — and also where its dependence on an
            accurately known voltage becomes acute.

        Parameters
        ----------
        bright_field_only:
            Keep only crossings inside the illumination cone, which are the only
            ones actually recorded.
        minimum_angle_deg:
            Discard pairs closer to parallel than this. Two exactly parallel
            lines have no crossing, and two nearly parallel ones have a crossing
            whose position is numerically meaningless.

        Returns
        -------
        tuple of HOLZLineIntersection
            Sorted by decreasing strain sensitivity: the most useful crossings
            first.
        """

        if not 0.0 < minimum_angle_deg < 90.0:
            raise ValueError("minimum_angle_deg must lie strictly between 0 and 90.")
        alpha = float(self.convergence_semi_angle_mrad) * 1e-3
        threshold = math.sin(math.radians(minimum_angle_deg))
        found: list[HOLZLineIntersection] = []
        for first_index, first in enumerate(self.lines):
            for second in self.lines[first_index + 1 :]:
                matrix = np.stack([first.normal_tilt, second.normal_tilt])
                determinant = float(np.linalg.det(matrix))
                if abs(determinant) < max(threshold, _PARALLEL_TOLERANCE):
                    continue
                offsets = np.array([first.offset_rad, second.offset_rad], dtype=np.float64)
                position = np.linalg.solve(matrix, offsets)
                derivative = np.linalg.solve(
                    matrix,
                    np.array(
                        [first.strain_sensitivity_rad, second.strain_sensitivity_rad],
                        dtype=np.float64,
                    ),
                )
                inside = bool(float(np.linalg.norm(position)) <= alpha)
                if bright_field_only and not inside:
                    continue
                found.append(
                    HOLZLineIntersection(
                        first_indices=first.miller_indices,
                        second_indices=second.miller_indices,
                        position_tilt_rad=np.asarray(position, dtype=np.float64),
                        angle_deg=float(math.degrees(math.asin(min(abs(determinant), 1.0)))),
                        strain_sensitivity_rad=float(np.linalg.norm(derivative)),
                        inside_bright_field=inside,
                    )
                )
        return tuple(sorted(found, key=lambda item: -item.strain_sensitivity_rad))

    def describe(self) -> str:
        """Convention-explicit prose: the lines, what they measure, and the degeneracy."""

        zone_label = format_direction_indices(
            tuple(int(value) for value in self.zone_axis.indices), style="plain"
        )
        visible = self.bright_field_lines
        if not visible:
            coverage = (
                "No line falls inside the illumination cone at this convergence semi-angle, "
                "so the bright-field disc carries no HOLZ detail. Widen the probe or accept "
                "that this zone shows none."
            )
        else:
            best = min(visible, key=lambda line: line.resolvable_strain(1000.0))
            best_indices = format_plane_indices(
                tuple(int(value) for value in best.miller_indices), style="plain"
            )
            coverage = (
                f"{len(visible)} of them cross the bright-field disc. The sharpest is "
                f"{best_indices}, "
                f"moving {best.strain_sensitivity_rad * 1e3:.3f} mrad per unit lattice strain "
                f"against a half-width of {best.angular_width_rad(1000.0) * 1e3:.4f} mrad in a "
                f"1000 A foil: a strain of {best.resolvable_strain(1000.0):.2e} shifts it by "
                "its own width."
            )
        return (
            f"HOLZ line geometry for {self.phase.name} down {zone_label} at "
            f"{self.beam_energy_kev:.0f} kV with a {self.convergence_semi_angle_mrad:.2f} mrad "
            f"probe: {len(self.lines)} lines from higher-order Laue zone reflections, and "
            f"{coverage} Line positions are exact given the excitation-error expression, but "
            "an isotropic lattice strain and a fractional wavelength change of the same size "
            "shift every line by exactly compensating amounts, so nothing here separates a "
            "strained lattice from an uncalibrated accelerating voltage. Contrast is not "
            "modelled: these are loci, not intensities."
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": HOLZ_LINE_PATTERN_SCHEMA,
            "phase": self.phase.name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "beam_energy_kev": self.beam_energy_kev,
            "convergence_semi_angle_mrad": self.convergence_semi_angle_mrad,
            "camera_constant_mm_angstrom": self.camera_constant_mm_angstrom,
            "line_count": len(self.lines),
            "bright_field_line_count": len(self.bright_field_lines),
            "lines": [
                {
                    "hkl": [int(value) for value in line.miller_indices],
                    "laue_zone": line.laue_zone,
                    "offset_mrad": line.offset_rad * 1e3,
                    "normal_tilt": line.normal_tilt.tolist(),
                    "crosses_bright_field": line.crosses_bright_field,
                    "strain_sensitivity_mrad": line.strain_sensitivity_rad * 1e3,
                }
                for line in self.lines
            ],
        }


def holz_line_pattern(
    phase: Phase,
    zone_axis: ZoneAxis,
    *,
    beam_energy_kev: float = 200.0,
    convergence_semi_angle_mrad: float = 8.0,
    camera_constant_mm_angstrom: float = 180.0,
    laue_zones: tuple[int, ...] = (1,),
    max_index: int = 24,
    g_max_inv_angstrom: float = 6.0,
    offset_margin_rad: float = 0.0,
    in_plane_rotation_deg: float = 0.0,
    apply_centering_absences: bool = True,
    provenance: ProvenanceRecord | None = None,
) -> HOLZLinePattern:
    """Compute the HOLZ lines of a zone-axis convergent-beam setting.

    What it does
        Enumerates higher-order Laue zone reflections, discards those forbidden
        by the lattice centering or beyond the ``|g|`` cut-off, and converts each
        into its exact Bragg locus in the plane of incident directions. The
        result carries every line, marked with whether it crosses the
        bright-field disc, so a caller can see what a larger probe would add.

    When to use it
        To predict which HOLZ lines a given zone will show; to choose a zone axis
        and convergence angle for a lattice-parameter or strain measurement; to
        identify a line in a recorded pattern by its position; and to quantify
        how much of a shift a candidate strain would produce before committing
        microscope time to measuring it.

    Parameters
    ----------
    phase, zone_axis:
        The zone axis must belong to the phase. Unlike a dynamical calculation
        this needs no unit cell: line positions are lattice geometry, and the
        atom positions affect only whether a line is visible, not where it is.
    beam_energy_kev:
        Accelerating voltage. HOLZ line positions are notoriously sensitive to
        it; see :meth:`HOLZLine.offset_at`.
    convergence_semi_angle_mrad:
        Illumination cone half-angle, which decides which lines are recorded.
    camera_constant_mm_angstrom:
        ``L * lambda``, for the detector-plane chords.
    laue_zones:
        Which Laue zones to search. ``(1,)`` is the first-order zone, the one
        normally used; ``0`` is rejected, because a zeroth-zone reflection has
        ``g_z = 0`` and its Bragg locus is a *disc boundary* effect rather than a
        line inside the pattern.
    max_index:
        Largest absolute Miller index enumerated. HOLZ reflections have large
        indices — the first ring of a cubic metal sits near ``|g| ~ 5`` inverse
        angstrom — so this bound must be generous or the search finds nothing.
    g_max_inv_angstrom:
        Radial cut-off.
    offset_margin_rad:
        Keep lines up to this far outside the illumination cone as well. Useful
        for showing what a wider probe would reveal; the default keeps only what
        would be recorded.
    in_plane_rotation_deg:
        Rotation of the zone basis about the beam, as in
        `pytex.diffraction.kinematic.zone_basis_from_axis`.
    apply_centering_absences:
        Remove reflections forbidden by the lattice centering. A forbidden
        reflection has no line, and reporting one would be worse than reporting
        none.
    provenance:
        Optional record.

    Returns
    -------
    HOLZLinePattern
        Lines sorted by increasing distance from the pattern centre.

    Raises
    ------
    ValueError
        If the zone axis belongs to a different phase, if ``0`` appears in
        ``laue_zones``, or if a bound is non-positive.

    Notes
    -----
    **Algorithm.**

    1. Build the zone basis ``(u, v, z)`` with ``z`` along the zone axis.
    2. Enumerate ``hkl`` within ``max_index``; drop centering-forbidden ones.
    3. Keep reflections whose Laue-zone order ``n = g . [uvw]`` is requested and
       whose ``|g|`` is inside the cut-off.
    4. For each, form the line normal ``g_perp / |g_perp|`` and the offset
       ``(g_z - lambda |g|^2 / 2) / |g_perp|``.
    5. Keep the lines within ``alpha + offset_margin`` of the centre, and sort.

    **Why no unit cell is required.** The offset depends only on the
    reciprocal-lattice geometry and the wavelength. Structure factors decide
    whether a line is *visible* and how dark it is, and that question belongs to
    `pytex.diffraction.dynamical`.
    """

    if zone_axis.phase != phase:
        raise ValueError("zone_axis.phase must match phase.")
    if any(int(zone) == 0 for zone in laue_zones):
        raise ValueError(
            "Laue zone 0 has no HOLZ lines: a zeroth-zone reflection has g_z = 0, so its "
            "Bragg locus is a fixed offset -lambda|g|^2/2 independent of the tilt direction "
            "and appears as a shift of the whole disc rather than as a line across it. Pass "
            "the higher-order zones, normally (1,)."
        )
    if not laue_zones:
        raise ValueError("laue_zones must name at least one Laue zone, normally 1.")
    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    if g_max_inv_angstrom <= 0.0:
        raise ValueError("g_max_inv_angstrom must be strictly positive.")
    if convergence_semi_angle_mrad <= 0.0:
        raise ValueError(
            "convergence_semi_angle_mrad must be strictly positive; a parallel beam samples "
            "one incident direction and therefore shows no line at all."
        )
    if offset_margin_rad < 0.0:
        raise ValueError("offset_margin_rad must be non-negative.")
    if camera_constant_mm_angstrom <= 0.0:
        raise ValueError("camera_constant_mm_angstrom must be strictly positive.")

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
    in_plane = np.linalg.norm(g_zone[:, :2], axis=1)

    finite = (g_magnitude <= g_max_inv_angstrom) & (in_plane > _PARALLEL_TOLERANCE)
    grid = grid[finite]
    g_zone = g_zone[finite]
    g_magnitude = g_magnitude[finite]
    laue_order = laue_order[finite]
    in_plane = in_plane[finite]

    offsets = (g_zone[:, 2] - 0.5 * wavelength * g_magnitude * g_magnitude) / in_plane
    alpha = float(convergence_semi_angle_mrad) * 1e-3
    keep = np.abs(offsets) <= alpha + offset_margin_rad
    order = np.argsort(np.abs(offsets[keep]), kind="stable")

    lines = tuple(
        HOLZLine(
            miller_indices=np.asarray(grid[keep][index], dtype=np.int64),
            laue_zone=int(laue_order[keep][index]),
            normal_tilt=np.asarray(
                g_zone[keep][index, :2] / in_plane[keep][index], dtype=np.float64
            ),
            offset_rad=float(offsets[keep][index]),
            g_perp_inv_angstrom=float(in_plane[keep][index]),
            g_magnitude_inv_angstrom=float(g_magnitude[keep][index]),
            g_zone_inv_angstrom=float(g_zone[keep][index, 2]),
            wavelength_angstrom=wavelength,
            convergence_semi_angle_rad=alpha,
            camera_constant_mm_angstrom=float(camera_constant_mm_angstrom),
            disc_centre_mm=np.asarray(
                camera_constant_mm_angstrom * g_zone[keep][index, :2], dtype=np.float64
            ),
        )
        for index in order
    )

    return HOLZLinePattern(
        phase=phase,
        zone_axis=zone_axis,
        beam_energy_kev=float(beam_energy_kev),
        convergence_semi_angle_mrad=float(convergence_semi_angle_mrad),
        camera_constant_mm_angstrom=float(camera_constant_mm_angstrom),
        lines=lines,
        zone_basis_crystal=basis,
        provenance=provenance,
    )
