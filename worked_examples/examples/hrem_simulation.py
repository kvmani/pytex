"""Worked examples: High-Resolution Transmission Electron Microscopy (HRTEM) simulation.

Four physical and mathematical invariants that calibrate an HRTEM simulation and
contrast transfer function (CTF) implementation:

- Relativistic electron wavelength, Scherzer defocus, and Scherzer point
  resolution: at 300 keV with Cs = 1.0 mm, analytical de Broglie wavelength
  and Scherzer balance formulas are reproduced exactly.
- CTF passband zero-crossing: under Scherzer underfocus, the wave aberration phase
  chi(q) vanishes identically at q_cross = sqrt(2.4 / sqrt(Cs lambda^3)), producing
  sin(chi) = 0 to floating-point precision at the edge of the first passband.
- Chromatic aberration damping envelope: in a double-corrected microscope with Cc
  correction reducing focal spread from 35 Å to 5 Å, the 1/e^2 information limit
  frequency extends by exactly sqrt(35 / 5) = sqrt(7).
- Phase contrast inversion: weak-phase object simulation correctly produces dark
  atomic columns under Scherzer conditions (negative phase contrast) and bright
  atomic columns under Negative Cs Imaging (NCSI) conditions.

See :doc:`../../theory/hrem_multislice_and_ctf` for the underlying physics.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

HREM_SETUP = """
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
"""

_DEFOCUS = SymbolUse(
    r"\Delta f",
    "Defocus of the objective lens; negative values correspond to underfocus.",
)
_CS = SymbolUse(
    r"C_{s}",
    "Third-order spherical aberration coefficient of the objective lens.",
)
_CC = SymbolUse(
    r"C_{c}",
    "Chromatic aberration coefficient of the objective lens.",
)
_ALPHA_OBJ = SymbolUse(
    r"\alpha_{\mathrm{obj}}",
    "Objective aperture semiangle cutoff.",
)
_LAMBDA = SymbolUse(
    r"\lambda",
    "Relativistic electron radiation wavelength.",
)

_THEORY = SeeAlso("HRTEM multislice and CTF theory", "../../theory/hrem_multislice_and_ctf")
_API = SeeAlso("Diffraction API", "../../api/index")


SCHERZER_OPTICS = WorkedExample(
    id="hrem-scherzer-optics-300kv",
    title="Relativistic wavelength, Scherzer defocus, and point resolution at 300 kV",
    domain="diffraction",
    scenario=(
        "In conventional high-resolution transmission electron microscopy (HRTEM), "
        "Scherzer defocus balances the phase shifts induced by the spherical aberration "
        "of the objective lens against the defocus term to maximize the bandwidth of "
        "constant phase contrast. For an accelerating voltage E_0 = 300 keV and third-order "
        "spherical aberration C_s = 1.0 mm, compute the relativistic electron wavelength "
        "lambda (in pm), the Scherzer underfocus Delta f_Sch = -1.2 sqrt(C_s lambda) (in nm), "
        "and the Scherzer point resolution d_Sch = 0.64 (C_s lambda^3)^(1/4) (in Å)."
    ),
    setup=HREM_SETUP,
    code=(
        "optics = MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0)\n"
        "result = [\n"
        "    round(optics.wavelength_angstrom * 1e3, 4),\n"
        "    round(abs(optics.scherzer_defocus_angstrom) * 0.1, 3),\n"
        "    round(optics.scherzer_resolution_angstrom, 3),\n"
        "]"
    ),
    expected=[19.6875, 53.245, 1.892],
    unit="",
    tolerance=1e-3,
    reference=(
        "Analytical relativistic de Broglie wavelength for 300 keV electrons "
        "(lambda = 1.96875 pm), Scherzer defocus Delta f_Sch = -53.245 nm, and "
        "Scherzer point resolution d_Sch = 1.892 Å."
    ),
    citation=(
        "Scherzer, O. (1949). The theoretical resolution limit of the electron microscope. "
        "J. Appl. Phys. 20, 20-29; Williams & Carter, Transmission Electron Microscopy, "
        "2nd ed. (Springer, 2009), Chapter 28."
    ),
    symbols=(_LAMBDA, _DEFOCUS, _CS),
    see_also=(_THEORY, _API),
    result_format="{:.3f}",
)


SCHERZER_ZERO_CROSSING = WorkedExample(
    id="hrem-ctf-scherzer-zero-crossing",
    title="Analytical zero crossing of the CTF at the Scherzer passband edge",
    domain="diffraction",
    scenario=(
        "Under Scherzer underfocus Delta f_Sch = -1.2 sqrt(C_s lambda), the wave aberration "
        "phase chi(q) = pi Delta f lambda q^2 + 0.5 pi C_s lambda^3 q^4 possesses an exact "
        "non-trivial root at spatial frequency q_cross = sqrt(2.4 / sqrt(C_s lambda^3)). "
        "At this frequency, chi(q_cross) = 0 and sin(chi(q_cross)) vanishes identically to "
        "floating-point precision, marking the physical edge of the primary broad phase-contrast "
        "passband before rapid spatial frequency oscillations begin."
    ),
    setup=HREM_SETUP,
    code=(
        "optics = MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0)\n"
        "lam = optics.wavelength_angstrom\n"
        "cs = optics.cs_angstrom\n"
        "q_cross = np.sqrt(2.4 / np.sqrt(cs * (lam**3)))\n"
        "chi = float(optics.wave_aberration(np.array([q_cross]))[0])\n"
        "result = float(abs(np.sin(chi)))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "An exact analytical identity: Delta f + 0.5 Cs lambda^2 q^2 = 0 implies "
        "chi(q_cross) = 0 and sin(chi) = 0 to machine precision."
    ),
    citation=(
        "Reimer, L. & Kohl, H. Transmission Electron Microscopy: Physics of Image Formation, "
        "5th ed. (Springer, 2008), Chapter 6."
    ),
    symbols=(_LAMBDA, _DEFOCUS, _CS),
    see_also=(_THEORY, _API),
    result_format="{:.2e}",
)


DOUBLE_CORRECTED_INFO_LIMIT = WorkedExample(
    id="hrem-double-corrected-information-limit",
    title="Double-corrected information limit extension by chromatic aberration damping reduction",
    domain="diffraction",
    scenario=(
        "In a conventional 300 kV TEM with focal spread Delta = 35 Å (dominated by Cc and "
        "energy spread Delta E), chromatic damping envelope E_c(q) = exp(-0.5 pi^2 lambda^2 Delta^2 q^4) "
        "rapidly suppresses high spatial frequencies. In a double-corrected instrument where "
        "chromatic aberration is corrected to Delta = 5 Å, the 1/e^2 temporal information cutoff "
        "frequency q_info = sqrt(2 / (pi lambda Delta)) scales as 1 / sqrt(Delta). The ratio of "
        "double-corrected to uncorrected cutoff frequency is an exact algebraic ratio sqrt(35 / 5) = sqrt(7)."
    ),
    setup=HREM_SETUP,
    code=(
        "uncorr = MicroscopeAberrations(energy_kev=300.0, cs_mm=1.0, focal_spread_angstrom=35.0)\n"
        "corr = MicroscopeAberrations(energy_kev=300.0, cs_mm=0.0, focal_spread_angstrom=5.0)\n"
        "q_uncorr = np.sqrt(2.0 / (np.pi * uncorr.wavelength_angstrom * 35.0))\n"
        "q_corr = np.sqrt(2.0 / (np.pi * corr.wavelength_angstrom * 5.0))\n"
        "ratio = q_corr / q_uncorr\n"
        "expected_ratio = np.sqrt(35.0 / 5.0)\n"
        "result = float(abs(ratio - expected_ratio))"
    ),
    expected=0.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "An algebraic scaling identity: for temporal coherence damping, the 1/e^2 "
        "cutoff frequency scales inversely with sqrt(Delta). The ratio of cutoff "
        "frequencies between Delta = 35 Å and Delta = 5 Å is sqrt(7) exactly."
    ),
    citation=(
        "Rose, H. (2009). Historical aspects of aberration correction. J. Electron Microsc. 58, 77-85; "
        "Haider, M. et al. (1998). Electron microscopy image enhanced. Nature 392, 768-769."
    ),
    symbols=(_LAMBDA, _CC, _DEFOCUS),
    see_also=(_THEORY, _API),
    result_format="{:.2e}",
)


NCSI_CONTRAST_INVERSION = WorkedExample(
    id="hrem-ncsi-contrast-inversion",
    title="Negative Cs Imaging (NCSI) atomic column contrast inversion",
    domain="diffraction",
    scenario=(
        "Under conventional Scherzer conditions (positive Cs, underfocus Delta f < 0), "
        "sin(chi) < 0 across the passband, causing atomic columns to appear as dark minima "
        "on a bright background in bright-field HRTEM. Under Negative Spherical Aberration "
        "Imaging (NCSI, Cs < 0, overfocus Delta f > 0), sin(chi) > 0 across the passband, "
        "causing atomic columns to appear as bright maxima on a dark background. Verify this "
        "contrast inversion for an isolated gold atom: the center intensity is below the "
        "mean background (< 1.0) for Scherzer and above the mean background (> 1.0) for NCSI."
    ),
    setup=HREM_SETUP,
    code=(
        "snap = AtomicSnapshot(\n"
        "    species=('Au',),\n"
        "    positions=np.array([[2.0, 2.0, 1.0]]),\n"
        "    cell=np.diag([4.0, 4.0, 4.0]),\n"
        ")\n"
        "scherzer_res = pure_python_phase_object_simulation(\n"
        "    snap, MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0), sampling_angstrom=0.1\n"
        ")\n"
        "ncsi_res = pure_python_phase_object_simulation(\n"
        "    snap, MicroscopeAberrations.ncsi(energy_kev=300.0), sampling_angstrom=0.1\n"
        ")\n"
        "ny, nx = scherzer_res.image.shape\n"
        "cx, cy = nx // 2, ny // 2\n"
        "scherzer_dark = float(scherzer_res.image[cy, cx] < 1.0)\n"
        "ncsi_bright = float(ncsi_res.image[cy, cx] > 1.0)\n"
        "result = [scherzer_dark, ncsi_bright]"
    ),
    expected=[1.0, 1.0],
    unit="",
    tolerance=0.0,
    reference=(
        "Physical principle of NCSI: phase contrast sign reversal between "
        "conventional underfocus HRTEM (dark atoms) and NCSI (bright atoms)."
    ),
    citation=(
        "Jia, C. L., Lentzen, M. & Urban, K. (2003). Atomic-resolution imaging of oxygen "
        "in perovskite ceramics. Science 299, 870-873; Urban, K. W. (2008). Studying "
        "microstructure with aberration-corrected transmission electron microscopy. "
        "Science 321, 506-510."
    ),
    symbols=(_CS, _DEFOCUS),
    see_also=(_THEORY, _API),
    result_format="{:.0f}",
)


GROUP = ExampleGroup(
    slug="hrem-simulation-and-ctf",
    title="HRTEM simulation and contrast transfer function optics",
    summary=(
        "Optics and contrast transfer of high-resolution transmission electron microscopy: "
        "relativistic electron wavelength, Scherzer defocus and resolution, analytical "
        "zero-crossing of the CTF, chromatic aberration damping reduction in double-corrected "
        "TEM, and phase contrast sign inversion in Negative Cs Imaging (NCSI)."
    ),
    examples=(
        SCHERZER_OPTICS,
        SCHERZER_ZERO_CROSSING,
        DOUBLE_CORRECTED_INFO_LIMIT,
        NCSI_CONTRAST_INVERSION,
    ),
)
