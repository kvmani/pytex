"""Worked examples: residual stress by the sin^2(psi) method.

Every expected value here has independent provenance. The isotropic constants,
the strain-free direction and the principal stresses are closed-form algebra
that can be done by hand; the Reuss constant of ferrite (211) is the textbook
cubic formula evaluated on tabulated single-crystal constants, worked through
in the scenario; the Kroener modulus is the positive root of Kroener's cubic;
the slope identity is the sin^2(psi) law itself; and the end-to-end example is
checked against the stress its synthetic data were generated with.

None of them is a copied output of this code.

See :doc:`../../theory/residual_stress_sin2psi`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

_THEORY = SeeAlso("Residual stress by the sin^2(psi) method", "../../theory/residual_stress_sin2psi")
_ALGORITHM = SeeAlso(
    "Residual stress: the algorithm and the report", "../../algorithms/residual_stress_sin2psi"
)

_PSI = SymbolUse(r"\psi", "Tilt of the scattering vector from the surface normal.")
_PHI = SymbolUse(r"\varphi", "Azimuth of the tilt plane, from S1 towards S2.")
_S1 = SymbolUse(r"S_{1}", "First diffraction elastic constant of the reflection.")
_HALF_S2 = SymbolUse(r"\tfrac{1}{2}S_{2}", "Second diffraction elastic constant of the reflection.")
_D0 = SymbolUse(r"d_{0}", "Stress-free interplanar spacing.")
_SIGMA_PHI = SymbolUse(r"\sigma_{\varphi}", "Normal stress along azimuth phi in the surface.")
_GAMMA = SymbolUse(r"\Gamma", "Cubic orientation parameter of a reflection.")

_SETUP_EXACT = (
    "import math\n"
    "import numpy as np\n"
    "from pytex.diffraction.xrd_residual_stress import (\n"
    "    DiffractionElasticConstants, StressPeak, determine_residual_stress,\n"
    "    measurement_direction,\n"
    ")\n"
    "WAVELENGTH = 2.2897\n"
    "def exact_peaks(stress, dec, d0, phis, psis):\n"
    "    peaks = []\n"
    "    for phi in phis:\n"
    "        for psi in psis:\n"
    "            m = measurement_direction(phi, psi)\n"
    "            strain = 1e-6 * (dec.half_s2_per_tpa * float(m @ stress @ m)\n"
    "                             + dec.s1_per_tpa * float(np.trace(stress)))\n"
    "            two_theta = math.degrees(2 * math.asin(WAVELENGTH / (2 * d0 * (1 + strain))))\n"
    "            peaks.append(StressPeak(phi_deg=phi, psi_deg=psi, two_theta_deg=two_theta,\n"
    "                                    two_theta_uncertainty_deg=1e-3, method='given'))\n"
    "    return peaks\n"
)


ISOTROPIC_CONSTANTS = WorkedExample(
    id="stress-isotropic-half-s2",
    title="The isotropic diffraction elastic constant of a steel",
    domain="residual-stress",
    scenario=(
        "For an elastically isotropic solid Hooke's law projected on the scattering vector gives "
        "1/2 S2 = (1 + nu)/E. For a steel with E = 210 GPa and nu = 0.28 that is "
        "1.28 / 210 GPa = 6.0952 x 10^-3 / GPa = 6.0952 / TPa: a 100 MPa stress along the "
        "scattering vector strains the planes by 6.1 x 10^-4 through this term."
    ),
    setup="from pytex.diffraction.xrd_residual_stress import DiffractionElasticConstants\n",
    code=(
        "result = DiffractionElasticConstants.isotropic(210.0, 0.28).half_s2_per_tpa"
    ),
    expected=6.0952,
    unit="TPa^-1",
    tolerance=1e-4,
    reference="1.28 / 210 = 0.0060952 per GPa = 6.0952 per TPa, by hand.",
    citation=(
        "Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6."
    ),
    symbols=(_HALF_S2,),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


REUSS_FERRITE_211 = WorkedExample(
    id="stress-reuss-ferrite-211",
    title="The Reuss constant 1/2 S2 of ferrite (211), by the cubic formula",
    domain="residual-stress",
    scenario=(
        "Under the Reuss model a cubic reflection has 1/2 S2 = S11 - S12 - 3 S0 Gamma with "
        "S0 = S11 - S12 - S44/2 and Gamma = (h^2 k^2 + k^2 l^2 + l^2 h^2)/(h^2 + k^2 + l^2)^2. For "
        "ferrite, C11 = 231.4, C12 = 134.7 and C44 = 116.4 GPa (Simmons & Wang). With "
        "(C11 - C12)(C11 + 2 C12) = 96.7 x 500.8 = 48427.4 GPa^2: S11 = 366.1 / 48427.4 = "
        "7.5598 /TPa, S12 = -134.7 / 48427.4 = -2.7815 /TPa and S44 = 1/116.4 = 8.5911 /TPa, so "
        "S0 = 6.0458 /TPa. For (211), Gamma = (4 + 1 + 4)/36 = 1/4, and 1/2 S2 = 10.3413 - "
        "3 x 6.0458 / 4 = 5.8069 /TPa. The library averages the single-crystal compliance over "
        "rotations about the plane normal in Mandel form; it must land on the same number."
    ),
    setup=(
        "from pytex.diffraction.xrd_residual_stress import (\n"
        "    DiffractionElasticConstants, single_crystal_stiffness,\n"
        ")\n"
    ),
    code=(
        "result = DiffractionElasticConstants.from_single_crystal(\n"
        "    single_crystal_stiffness('fe_bcc'), [2, 1, 1], model='reuss'\n"
        ").half_s2_per_tpa"
    ),
    expected=5.8069,
    unit="TPa^-1",
    tolerance=2e-4,
    reference=(
        "By hand from the cubic Reuss formula: S11 - S12 = 10.3413, 3 S0 Gamma = 4.5344, "
        "1/2 S2 = 5.8069 /TPa."
    ),
    citation=(
        "Welzel et al., J. Appl. Cryst. 38 (2005) 1, doi:10.1107/S0021889804029516; "
        "Simmons & Wang, Single Crystal Elastic Constants, MIT Press (1971)."
    ),
    symbols=(_HALF_S2, _GAMMA),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


KROENER_SHEAR_MODULUS = WorkedExample(
    id="stress-kroener-shear-modulus-ferrite",
    title="Kroener's self-consistent shear modulus of a ferrite aggregate",
    domain="residual-stress",
    scenario=(
        "For cubic crystals Kroener's self-consistent shear modulus is the positive root of "
        "G^3 + (5 C11 + 4 C12)/8 G^2 - C44 (7 C11 - 4 C12)/8 G - C44 (C11 - C12)(C11 + 2 C12)/8 "
        "= 0. For ferrite the coefficients are 211.975, -15728.55 and -704618.09 (GPa units), "
        "and the positive root is G = 82.448 GPa, between the Reuss (74.47) and Voigt (89.18) "
        "shear moduli. "
        "The library does not solve this cubic: it iterates the Eshelby-sphere concentration "
        "tensor in the general (any crystal system) scheme, which must converge to the root."
    ),
    setup=(
        "import numpy as np\n"
        "from pytex.diffraction.xrd_residual_stress import (\n"
        "    _kroener_grain_compliance, _to_mandel, single_crystal_stiffness,\n"
        ")\n"
    ),
    code=(
        "stiffness = np.asarray(single_crystal_stiffness('fe_bcc').tensor)\n"
        "_, bulk, shear = _kroener_grain_compliance(_to_mandel(stiffness))\n"
        "result = shear"
    ),
    expected=82.448,
    unit="GPa",
    tolerance=1e-3,
    reference=(
        "Substituting G = 82.448 into the cubic: 560 454.5 + 1 440 936.7 - 1 296 787.5 - "
        "704 618.1 = -14.4, against a derivative 3G^2 + 2 alpha G + beta = 39 618 per GPa, so "
        "the root is 82.448 + 0.0004 GPa."
    ),
    citation="Kroener, Z. Physik 151 (1958) 504, doi:10.1007/BF01337948.",
    symbols=(),
    see_also=(_THEORY,),
    result_format="{:.3f}",
)


SLOPE_IS_THE_STRESS = WorkedExample(
    id="stress-sin2psi-slope-identity",
    title="The slope of d against sin^2(psi) is d0 1/2 S2 sigma_phi",
    domain="residual-stress",
    scenario=(
        "Under plane stress d = d0 [1 + 1/2 S2 sigma_phi sin^2(psi) + S1 (sigma_11 + sigma_22)], "
        "which is exactly linear in sin^2(psi). For a uniaxial sigma_11 = -300 MPa measured at "
        "phi = 0 with 1/2 S2 = 5.8 /TPa and d0 = 1.1702 angstrom, the slope is "
        "1.1702 x 5.8 x 10^-6 x (-300) = -2.03615 x 10^-3 angstrom. Exact peak positions are "
        "generated from the fundamental equation, and the per-azimuth regression must return "
        "that slope."
    ),
    setup=_SETUP_EXACT,
    code=(
        "dec = DiffractionElasticConstants(s1_per_tpa=-1.25, half_s2_per_tpa=5.8)\n"
        "stress = np.diag([-300.0, 0.0, 0.0])\n"
        "peaks = exact_peaks(stress, dec, 1.1702, (0.0,), (0.0, 15.0, 25.0, 35.0, 45.0))\n"
        "fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,\n"
        "                                d0_angstrom=1.1702, dec=dec)\n"
        "result = 1000.0 * fit.regressions[0].slope_angstrom"
    ),
    expected=-2.03615,
    unit="mÅ",
    tolerance=2e-5,
    reference="1.1702 x 5.8e-6 x (-300) = -2.03615e-3 angstrom, by hand.",
    citation="Macherauch & Mueller, Z. angew. Phys. 13 (1961) 305.",
    symbols=(_D0, _HALF_S2, _SIGMA_PHI, _PSI),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.5f}",
)


STRAIN_FREE_DIRECTION = WorkedExample(
    id="stress-strain-free-direction",
    title="The strain-free tilt of an equibiaxial stress",
    domain="residual-stress",
    scenario=(
        "Under an equibiaxial stress the sin^2(psi) line crosses d = d0 at "
        "sin^2(psi*) = -2 S1 / (1/2 S2), independent of the stress. For isotropic constants that "
        "is 2 nu / (1 + nu); with nu = 0.28, 0.56 / 1.28 = 0.4375 (psi* = 41.4 degrees). The "
        "tilt is read off the fitted line as (d0 - intercept) / slope."
    ),
    setup=_SETUP_EXACT,
    code=(
        "dec = DiffractionElasticConstants.isotropic(210.0, 0.28)\n"
        "stress = np.diag([-400.0, -400.0, 0.0])\n"
        "peaks = exact_peaks(stress, dec, 1.1702, (0.0, 45.0, 90.0),\n"
        "                    (0.0, 15.0, 25.0, 35.0, 45.0))\n"
        "fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,\n"
        "                                d0_angstrom=1.1702, dec=dec)\n"
        "line = fit.regressions[0]\n"
        "result = (1.1702 - line.intercept_angstrom) / line.slope_angstrom"
    ),
    expected=0.4375,
    unit="",
    tolerance=1e-4,
    reference="2 x 0.28 / 1.28 = 0.4375, by hand.",
    citation="Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6.",
    symbols=(_S1, _HALF_S2, _PSI),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


PRINCIPAL_STRESS = WorkedExample(
    id="stress-in-plane-principal",
    title="The larger in-plane principal stress from three azimuths",
    domain="residual-stress",
    scenario=(
        "Three azimuths, 0, 45 and 90 degrees, fix sigma_11, sigma_22 and sigma_12, and the "
        "in-plane principal stresses are (sigma_11 + sigma_22)/2 +/- sqrt(((sigma_11 - "
        "sigma_22)/2)^2 + sigma_12^2). For sigma_11 = -350, sigma_22 = -150 and sigma_12 = 60 MPa "
        "the larger is -250 + sqrt(100^2 + 60^2) = -250 + 116.619 = -133.381 MPa."
    ),
    setup=_SETUP_EXACT,
    code=(
        "dec = DiffractionElasticConstants(s1_per_tpa=-1.23, half_s2_per_tpa=5.69)\n"
        "stress = np.array([[-350.0, 60.0, 0.0], [60.0, -150.0, 0.0], [0.0, 0.0, 0.0]])\n"
        "peaks = exact_peaks(stress, dec, 1.1702, (0.0, 45.0, 90.0),\n"
        "                    (-45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0))\n"
        "fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,\n"
        "                                d0_angstrom=1.1702, dec=dec)\n"
        "result = fit.tensor.in_plane_principal()['sigma_I_mpa']"
    ),
    expected=-133.381,
    unit="MPa",
    tolerance=1e-3,
    reference="-250 + sqrt(13600) = -250 + 116.619 = -133.381 MPa, by hand.",
    citation="Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6.",
    symbols=(_PHI,),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.3f}",
)


END_TO_END_FERRITE = WorkedExample(
    id="stress-end-to-end-ferrite",
    title="sigma_11 of a shot-peened ferrite, from simulated Cr K-alpha scans",
    domain="residual-stress",
    scenario=(
        "Twenty-one scans of ferrite (211) with Cr K-alpha - three azimuths, seven tilts of both "
        "signs - generated from sigma_11 = -350, sigma_22 = -150, sigma_12 = 60 MPa through the "
        "Kroener constants, as K-alpha1/K-alpha2 pseudo-Voigt doublets broadened as 1/cos(psi), "
        "shaped by the LPA factor and given Poisson noise. Each peak is located by the doublet "
        "profile fit after the LPA correction, and the tensor fitted to all 21 positions. The "
        "answer must be the stress the data were generated with, to within the statistical "
        "uncertainty of a few MPa."
    ),
    setup=(
        "import math\n"
        "from pytex.diffraction.xrd_residual_stress import (\n"
        "    DiffractionElasticConstants, residual_stress_pipeline,\n"
        "    simulate_sin2psi_measurement, single_crystal_stiffness,\n"
        ")\n"
    ),
    code=(
        "d0 = 2.8665 / math.sqrt(6.0)\n"
        "dec = DiffractionElasticConstants.from_single_crystal(\n"
        "    single_crystal_stiffness('fe_bcc'), [2, 1, 1], model='kroener')\n"
        "scans = simulate_sin2psi_measurement(\n"
        "    d0_angstrom=d0, dec=dec, seed=1,\n"
        "    stress_mpa={'sigma_11': -350.0, 'sigma_22': -150.0, 'sigma_12': 60.0})\n"
        "fit = residual_stress_pipeline(scans, d0_angstrom=d0, dec=dec, window_deg=8.0)\n"
        "result = fit.tensor.component('sigma_11')[0]"
    ),
    expected=-350.0,
    unit="MPa",
    tolerance=8.0,
    reference=(
        "The generating stress. The tolerance is about four statistical standard uncertainties "
        "of the profile-fit route on this data set."
    ),
    citation="Macherauch & Mueller, Z. angew. Phys. 13 (1961) 305.",
    symbols=(_PSI, _PHI),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.1f}",
)


GROUP = ExampleGroup(
    slug="residual-stress",
    title="Residual stress by the sin^2(psi) method",
    summary=(
        "The diffraction elastic constants checked against the isotropic formula, the cubic "
        "Reuss formula worked by hand and the root of Kroener's cubic; the sin^2(psi) law's "
        "slope, strain-free tilt and principal stresses checked against closed-form algebra on "
        "exact data; and an end-to-end evaluation of noisy simulated scans checked against the "
        "stress they were generated with."
    ),
    examples=(
        ISOTROPIC_CONSTANTS,
        REUSS_FERRITE_211,
        KROENER_SHEAR_MODULUS,
        SLOPE_IS_THE_STRESS,
        STRAIN_FREE_DIRECTION,
        PRINCIPAL_STRESS,
        END_TO_END_FERRITE,
    ),
)
