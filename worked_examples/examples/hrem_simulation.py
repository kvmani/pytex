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

_C12 = SymbolUse(
    r"C_{12}",
    "Two-fold astigmatism amplitude of the objective lens.",
)
_PHI12 = SymbolUse(
    r"\varphi_{12}",
    "Azimuth of the two-fold astigmatism axis in the back focal plane.",
)
_THETA_Q = SymbolUse(
    r"\theta_{q}",
    "Azimuth in the back focal plane at which a transfer profile is cut.",
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


ASTIGMATISM_EQUIVALENT_DEFOCUS = WorkedExample(
    id="hrem-astigmatism-equivalent-defocus",
    title="Two-fold astigmatism acts as a defocus offset along its own azimuth",
    domain="diffraction",
    scenario=(
        "Two-fold astigmatism enters the wave aberration as "
        "pi lambda q^2 C_12 cos(2(theta - phi_12)). Along its own azimuth, theta = phi_12, "
        "the cosine equals +1 and the term is algebraically indistinguishable from adding "
        "C_12 to the defocus; ninety degrees away the cosine equals -1 and it subtracts the "
        "same amount. For a 300 kV lens at Delta f = -50 Å with Cs = 1 um and C_12 = 20 Å at "
        "phi_12 = 30 degrees, the point resolution of the cut along phi_12 must therefore "
        "equal that of a round lens at Delta f = -30 Å, and the cut across it that of a "
        "round lens at Delta f = -70 Å. Compute the difference between each cut and its "
        "equivalent round lens."
    ),
    setup=HREM_SETUP,
    code=(
        "astigmatic = MicroscopeAberrations(\n"
        "    energy_kev=300.0,\n"
        "    defocus_angstrom=-50.0,\n"
        "    cs_mm=0.001,\n"
        "    astigmatism_angstrom=20.0,\n"
        "    astigmatism_angle_deg=30.0,\n"
        ")\n"
        "def round_lens(defocus):\n"
        "    return MicroscopeAberrations(energy_kev=300.0, defocus_angstrom=defocus, cs_mm=0.001)\n"
        "def d0(lens, azimuth=0.0):\n"
        "    return lens.evaluate_ctf_1d(2.5, 4000, azimuth_deg=azimuth).point_resolution_angstrom\n"
        "result = [\n"
        "    abs(d0(astigmatic, 30.0) - d0(round_lens(-30.0))),\n"
        "    abs(d0(astigmatic, 120.0) - d0(round_lens(-70.0))),\n"
        "]"
    ),
    expected=[0.0, 0.0],
    unit="Å",
    tolerance=1e-12,
    reference=(
        "An exact algebraic identity of the wave aberration function: at theta = phi_12 the "
        "C_12 term reduces to pi lambda q^2 C_12, which is the defocus term with "
        "Delta f -> Delta f + C_12, so the two lenses have identical chi and hence identical "
        "zero crossings."
    ),
    citation=(
        "Krivanek, O. L., Dellby, N. & Lupini, A. R. (1999). Towards sub-A electron beams. "
        "Ultramicroscopy 78, 1-11; Kirkland, E. J. Advanced Computing in Electron Microscopy, "
        "2nd ed. (Springer, 2010), Chapter 3."
    ),
    symbols=(_C12, _PHI12, _THETA_Q, _DEFOCUS),
    see_also=(_THEORY, _API),
    result_format="{:.2e}",
)


ASTIGMATIC_RESOLUTION_SPLIT = WorkedExample(
    id="hrem-astigmatic-resolution-split",
    title="Orthogonal split of point resolution under two-fold astigmatism",
    domain="diffraction",
    scenario=(
        "Because the C_12 term varies as cos(2 theta), its period in azimuth is 180 degrees "
        "and its two extremes lie 90 degrees apart. A lens is therefore resolved best along "
        "one direction and worst along the perpendicular one, and the separation between "
        "those two directions is fixed by the multiplicity of the aberration rather than by "
        "its size. Sample the transfer function of the same 300 kV lens over 180 azimuths and "
        "report the angular separation, modulo 180 degrees, between the finest and the "
        "coarsest point resolution."
    ),
    setup=HREM_SETUP,
    code=(
        "lens = MicroscopeAberrations(\n"
        "    energy_kev=300.0,\n"
        "    defocus_angstrom=-50.0,\n"
        "    cs_mm=0.001,\n"
        "    astigmatism_angstrom=20.0,\n"
        "    astigmatism_angle_deg=30.0,\n"
        ")\n"
        "band = lens.evaluate_ctf_azimuthal(2.5, 4000, 180)\n"
        "resolutions = band.point_resolution_angstrom\n"
        "best = float(band.azimuths_deg[int(np.nanargmin(resolutions))])\n"
        "result = float(abs(band.worst_azimuth_deg - best) % 180.0)"
    ),
    expected=90.0,
    unit="deg",
    tolerance=1.0,
    reference=(
        "The azimuthal multiplicity of two-fold astigmatism: cos(2 theta) has period 180 "
        "degrees, so its maximum and minimum are separated by exactly 90 degrees. The "
        "tolerance is one degree, the spacing of the 180-point azimuthal grid."
    ),
    citation=(
        "Krivanek, O. L., Dellby, N. & Lupini, A. R. (1999). Towards sub-A electron beams. "
        "Ultramicroscopy 78, 1-11."
    ),
    symbols=(_C12, _PHI12, _THETA_Q),
    see_also=(_THEORY, _API),
    result_format="{:.1f}",
)


GROUP = ExampleGroup(
    slug="hrem-simulation-and-ctf",
    title="HRTEM simulation and contrast transfer function optics",
    summary=(
        "Optics and contrast transfer of high-resolution transmission electron microscopy: "
        "relativistic electron wavelength, Scherzer defocus and resolution, analytical "
        "zero-crossing of the CTF, chromatic aberration damping reduction in double-corrected "
        "TEM, phase contrast sign inversion in Negative Cs Imaging (NCSI), and the "
        "directional splitting of point resolution that a residual two-fold astigmatism "
        "imposes on a corrected lens."
    ),
    examples=(
        SCHERZER_OPTICS,
        SCHERZER_ZERO_CROSSING,
        DOUBLE_CORRECTED_INFO_LIMIT,
        NCSI_CONTRAST_INVERSION,
        ASTIGMATISM_EQUIVALENT_DEFOCUS,
        ASTIGMATIC_RESOLUTION_SPLIT,
    ),
)
