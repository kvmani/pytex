"""Worked examples: multislice HRTEM simulation.

Six numbers a multislice implementation must reproduce, each with an expected
value that comes from a closed form or from a different method, never from a
previous run of the engine:

- the projected potential of one atom integrates to (h^2 / 2 pi m0 e) f_e(0);
- the independent-atom mean inner potential of silicon follows from the same
  constant and the unit-cell volume;
- the band limit fixes the largest scattering angle a grid represents;
- imaging at defocus is Fresnel propagation of the exit wave by that distance;
- multislice and Bloch waves agree on identical potentials;
- integrating the focal spread reproduces Frank's temporal envelope for a
  linear (weak) fringe.

See :doc:`../../theory/multislice_hrtem` for the derivations and
:doc:`../../algorithms/multislice_hrtem` for the algorithm.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

MULTISLICE_SETUP = """
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
"""

_V_SLICE = SymbolUse(r"v_{n}", "Projected potential of slice n, in V Å.")
_FE = SymbolUse(r"f_{e}(s)", "Electron atomic scattering factor in ångström.")
_V0 = SymbolUse(r"V_{0}", "Mean inner potential of a crystal, in volts.")
_GMAX = SymbolUse(r"g_{\max}", "Band limit of the multislice grid, two thirds of Nyquist.")
_LAMBDA = SymbolUse(r"\lambda", "Relativistic electron wavelength.")
_DEFOCUS = SymbolUse(r"\Delta f", "Objective lens defocus; negative is underfocus.")
_DZ = SymbolUse(r"\Delta z", "Slice thickness of the multislice calculation.")
_EC = SymbolUse(r"E_{c}", "Temporal coherence (focal spread) envelope.")
_SPREAD = SymbolUse(r"\Delta", "Focal spread: standard deviation of the defocus distribution.")

_THEORY = SeeAlso("Multislice HRTEM theory", "../../theory/multislice_hrtem")
_ALGORITHM = SeeAlso("Multislice HRTEM algorithm", "../../algorithms/multislice_hrtem")


POTENTIAL_INTEGRAL = WorkedExample(
    id="multislice-projected-potential-integral",
    title="The projected potential of one atom integrates to (h²/2πm₀e) f_e(0)",
    domain="diffraction",
    scenario=(
        "The independent-atom potential is the Fourier transform of the electron scattering "
        "factor scaled by h^2/(2 pi m0 e) = 2 pi a0 e = 47.878 V Å^2, so the zero-frequency "
        "coefficient of a projected potential - its integral over the plane - is that "
        "constant times f_e(0). Place one silicon atom anywhere in a periodic 8 x 7 Å cell, "
        "build its projected potential from the Lobato-Van Dyck parametrization, and "
        "integrate it over the cell in V Å^3."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "grid = MultisliceGrid((140, 160), (8.0, 7.0))\n"
        "v = slice_potential(['Si'], np.array([[2.3, 4.1]]), grid, 'lobato')\n"
        "dx, dy = grid.sampling_angstrom\n"
        "result = float(np.sum(v) * dx * dy)"
    ),
    expected=279.414,
    unit="V Å^3",
    tolerance=1e-3,
    reference=(
        "2 pi a0 e x f_e(0) with a0 = 0.529177 Å, e^2/(4 pi eps0) = 14.39965 eV Å (CODATA) and "
        "f_e(0) = 2 sum(a_i) = 5.8360 Å from the Lobato-Van Dyck coefficients of Si: "
        "47.8776 x 5.8360 = 279.414 V Å^3."
    ),
    citation=(
        "Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, 2nd ed., "
        "eq. 5.9; Lobato, I. & Van Dyck, D. (2014). Acta Cryst. A70, 636-649."
    ),
    symbols=(_V_SLICE, _FE),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.3f}",
)


MEAN_INNER_POTENTIAL = WorkedExample(
    id="multislice-silicon-mean-inner-potential",
    title="Independent-atom mean inner potential of silicon",
    domain="diffraction",
    scenario=(
        "Averaged over the crystal, the potential of eight silicon atoms per cubic cell of "
        "edge a = 5.43102 Å is V0 = (h^2/2 pi m0 e) 8 f_e(0) / a^3 in the independent-atom "
        "model. Build a periodic silicon slab along [001], slice it, and report the mean "
        "inner potential the sliced potential carries, in volts."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "phase = phase_from_request({'builtin': 'si_diamond'})[1]\n"
        "snap, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=3)\n"
        "grid = MultisliceGrid.from_sampling((snap.cell[0, 0], snap.cell[1, 1]), 0.1)\n"
        "sliced = SlicedPotential.from_snapshot(snap, grid, phase.lattice.a / 4)\n"
        "result = sliced.mean_inner_potential_volt()"
    ),
    expected=13.954,
    unit="V",
    tolerance=1e-3,
    reference=(
        "47.8776 V Å^2 x 8 x 5.8360 Å / (5.43102 Å)^3 = 13.954 V. The independent-atom "
        "model ignores bonding, which lowers the measured value to about 12 V: the gap is a "
        "known limitation of the model, not of the multislice."
    ),
    citation=(
        "Kirkland (2010), sec. 5.4; Gajdardziska-Josifovska, M. et al. (1993). "
        "Ultramicroscopy 50, 285-299 (measured mean inner potentials)."
    ),
    symbols=(_V0, _FE),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.3f}",
)


BAND_LIMIT_ANGLE = WorkedExample(
    id="multislice-band-limit-scattering-angle",
    title="The largest scattering angle a 0.05 Å grid represents at 200 kV",
    domain="diffraction",
    scenario=(
        "Multislice band-limits the transmission function and the propagator to two thirds "
        "of the Nyquist frequency 1/(2 dx) so the product of wave and transmission function "
        "cannot alias. The largest scattering semi-angle the calculation represents is then "
        "alpha_max = lambda g_max. Report it in mrad for a 0.05 Å pixel at 200 kV."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "grid = MultisliceGrid.from_sampling((10.0, 10.0), 0.05)\n"
        "result = grid.max_scattering_angle_mrad(200.0)"
    ),
    expected=167.196,
    unit="mrad",
    tolerance=1e-3,
    reference=(
        "alpha_max = lambda (2/3) / (2 dx) = 0.0250793 Å x 6.6667 Å^-1 = 167.196 mrad, with "
        "the relativistic wavelength at 200 kV."
    ),
    citation="Kirkland (2010), sec. 6.8, the two-thirds band-limit rule.",
    symbols=(_GMAX, _LAMBDA),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.3f}",
)


DEFOCUS_IS_PROPAGATION = WorkedExample(
    id="multislice-defocus-is-fresnel-propagation",
    title="Imaging at defocus is Fresnel propagation of the exit wave",
    domain="diffraction",
    scenario=(
        "With no spherical aberration and full coherence, the objective lens multiplies the "
        "exit-wave spectrum by exp(-i pi lambda Delta f g^2), which is exactly the free-space "
        "propagator over a distance Delta f. A focal series therefore needs one multislice "
        "run. Image a silicon [110] slab at Delta f = -150 Å and report the largest "
        "difference from |exit wave propagated by -150 Å|^2."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "phase = phase_from_request({'builtin': 'si_diamond'})[1]\n"
        "snap, _, _ = periodic_slab(phase, (1, 1, 0), (1, 1), thickness_angstrom=30.0)\n"
        "wave = multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=1.0)\n"
        "lens = MicroscopeAberrations(energy_kev=200.0, defocus_angstrom=-150.0, cs_mm=0.0,\n"
        "                             focal_spread_angstrom=0.0, convergence_semiangle_mrad=0.0)\n"
        "g2 = wave.grid.frequency_magnitude() ** 2\n"
        "free = np.fft.ifft2(np.fft.fft2(wave.wave())\n"
        "                    * np.exp(-1j * math.pi * wave.wavelength_angstrom * -150.0 * g2))\n"
        "result = float(np.max(np.abs(wave.image(lens) - np.abs(free) ** 2)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-10,
    reference=(
        "Identity of the defocus term of the wave aberration, pi lambda Delta f g^2, with the "
        "phase of the Fresnel propagator over Delta f: the two images are the same function."
    ),
    citation="Kirkland (2010), sec. 3.3 and eq. 6.92.",
    symbols=(_DEFOCUS, _LAMBDA),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.1e}",
)


BLOCH_AGREEMENT = WorkedExample(
    id="multislice-agrees-with-bloch-waves",
    title="Multislice and Bloch waves agree on the same crystal potential",
    domain="diffraction",
    scenario=(
        "Multislice and the Bloch-wave method are two exact solutions of the same "
        "high-energy Schrodinger equation. Give both the same Mott-Bethe potential for silicon "
        "along [001] at 200 kV, keep only the zero-order Laue zone (the projected potential of "
        "one period, here split into eight equal slices so that the splitting error is "
        "negligible), propagate through 37 periods (201 Å), and report the largest difference "
        "between the two methods' intensities of the transmitted beam, 220 and 400."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "import scipy.fft\n"
        "from pytex.core.lattice import ZoneAxis\n"
        "from pytex.diffraction.dynamical import beam_set_for_zone, solve_bloch_waves\n"
        "from pytex.diffraction.hrem import relativistic_interaction_parameter_inv_v_angstrom\n"
        "from pytex.diffraction.multislice import antialias_aperture, fresnel_propagator\n"
        "phase = phase_from_request({'builtin': 'si_diamond'})[1]\n"
        "a = phase.lattice.a\n"
        "cell, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1))\n"
        "grid = MultisliceGrid((64, 64), (a, a))\n"
        "v = slice_potential(cell.species, cell.positions[:, :2], grid, 'mott_bethe')\n"
        "sigma = relativistic_interaction_parameter_inv_v_angstrom(200.0)\n"
        "tau = scipy.fft.ifft2(scipy.fft.fft2(np.exp(1j * sigma * v / 8))\n"
        "                      * antialias_aperture(grid))\n"
        "propagator = fresnel_propagator(grid, 200.0, a / 8)\n"
        "psi = np.ones(grid.shape, dtype=complex)\n"
        "for _ in range(37 * 8):\n"
        "    psi = scipy.fft.ifft2(scipy.fft.fft2(psi * tau) * propagator)\n"
        "spectrum = scipy.fft.fft2(psi, norm='forward')\n"
        "hkl = [(0, 0, 0), (2, 2, 0), (4, 0, 0)]\n"
        "ms = np.array([abs(spectrum[k, h]) ** 2 for h, k, _ in hkl])\n"
        "beams = beam_set_for_zone(phase, ZoneAxis(indices=(0, 0, 1), phase=phase),\n"
        "                          beam_energy_kev=200.0, max_index=24,\n"
        "                          g_max_inv_angstrom=3.9,\n"
        "                          max_excitation_error_inv_angstrom=0.6)\n"
        "bloch = solve_bloch_waves(beams, [[0.0, 0.0]], thickness_angstrom=37 * a)\n"
        "bw = np.array([float(bloch.intensity_of(h)[0]) for h in hkl])\n"
        "result = float(np.max(np.abs(ms - bw)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1.5e-3,
    reference=(
        "Two exact solutions of one equation on one potential must agree. The tolerance "
        "covers the Bloch solver's finite beam set and its cos(theta) normalisation. With one "
        "slice per period instead of eight, the multislice differs by about 1 % at 108 Å: the "
        "splitting error of a 5.4 Å slice, which is why slices are kept thin."
    ),
    citation=(
        "Kirkland (2010), ch. 6 (multislice) and ch. 7 (Bloch waves); "
        "Self, P. G. et al. (1983). Ultramicroscopy 11, 35-52."
    ),
    symbols=(_DZ, _V_SLICE),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.4f}",
)


FOCAL_INTEGRATION_ENVELOPE = WorkedExample(
    id="multislice-focal-integration-reproduces-frank-envelope",
    title="Integrating the focal spread reproduces Frank's temporal envelope",
    domain="diffraction",
    scenario=(
        "For a linear (weak) fringe the incoherent average over a Gaussian defocus spread of "
        "standard deviation Delta multiplies its contrast by Frank's envelope "
        "E_c = exp(-pi^2 lambda^2 Delta^2 g^4 / 2). Image a weak 1 Å^-1 fringe at 300 kV "
        "and Delta = 30 Å by explicit Gauss-Hermite focal integration and report the "
        "fringe contrast relative to the coherent contrast."
    ),
    setup=MULTISLICE_SETUP,
    code=(
        "grid = MultisliceGrid((64, 200), (10.0, 6.4))\n"
        "x = np.arange(200) * 0.05\n"
        "eps = 1e-4\n"
        "wave = np.broadcast_to(1.0 + eps * np.cos(2 * math.pi * x), (64, 200))\n"
        "lens = MicroscopeAberrations(energy_kev=300.0, defocus_angstrom=0.0, cs_mm=0.0,\n"
        "                             focal_spread_angstrom=30.0, convergence_semiangle_mrad=0.0)\n"
        "image = hrtem_image(wave, grid, lens,\n"
        "                    temporal_coherence=TemporalCoherence.FOCAL_INTEGRATION)\n"
        "fringe = 2 * abs(np.fft.fft2(image)[0, 10]) / image.size\n"
        "result = float(fringe / (2 * eps))"
    ),
    expected=0.17881,
    unit="",
    tolerance=1e-4,
    reference=(
        "exp(-0.5 pi^2 lambda^2 Delta^2 g^4) with lambda = 0.0196875 Å (300 kV), Delta = 30 Å, "
        "g = 1 Å^-1: exp(-1.72146) = 0.17881."
    ),
    citation=(
        "Frank, J. (1973). The envelope of electron microscopic transfer functions for "
        "partially coherent illumination. Optik 38, 519-536."
    ),
    symbols=(_EC, _SPREAD, _LAMBDA),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.5f}",
)


GROUP = ExampleGroup(
    slug="multislice-hrtem",
    title="Multislice HRTEM simulation",
    summary=(
        "The multislice engine against closed forms and a second method: the absolute scale "
        "of the projected potential and the mean inner potential, the band limit, defocus as "
        "Fresnel propagation, agreement with Bloch waves on the same potential, and focal "
        "integration against Frank's envelope."
    ),
    examples=(
        POTENTIAL_INTEGRAL,
        MEAN_INNER_POTENTIAL,
        BAND_LIMIT_ANGLE,
        DEFOCUS_IS_PROPAGATION,
        BLOCH_AGREEMENT,
        FOCAL_INTEGRATION_ENVELOPE,
    ),
)
