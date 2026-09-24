"""Residual stress from lattice strain: the sin-squared-psi method.

This module answers one question: given the position of one reflection
measured at several specimen tilts ``psi`` and azimuths ``phi``, and the
stress-free spacing ``d0`` of that reflection, what is the stress in the
irradiated surface layer, and how well is it known?

Why this needs its own module
-----------------------------
A symmetric theta-2theta scan measures the spacing of planes parallel to the
surface only. A stress lying *in* the surface is visible there only through the
Poisson contraction, which is small and inseparable from a composition or
temperature change of ``d0``. Tilting the specimen brings planes inclined to the
surface into the diffracting condition, and it is the *change* of their spacing
with inclination that measures the stress. The spacing is the strain gauge; the
tilt is what orients it.

The measured strain
-------------------
The scattering vector at azimuth ``phi`` (measured in the surface from the
specimen axis ``S1``) and tilt ``psi`` (from the surface normal ``S3``) is the
specimen-frame unit vector

``m = (cos(phi) sin(psi), sin(phi) sin(psi), cos(psi))``

and a reflection whose planes are normal to ``m`` measures the normal strain

``epsilon_phi_psi = (d_phi_psi - d0) / d0 = m . epsilon . m``.

For a quasi-isotropic (untextured) aggregate the strain seen by the grains that
diffract is linear in the macroscopic stress through two *diffraction elastic
constants* of the reflection, ``S1(hkl)`` and ``1/2 S2(hkl)``:

``epsilon_phi_psi = 1/2 S2 m . sigma . m + S1 tr(sigma)``

which, written out, is the fundamental equation of X-ray stress analysis:

``epsilon_phi_psi = 1/2 S2 [(sigma_phi - sigma_33) sin^2(psi) + tau_phi sin(2 psi)]
                   + 1/2 S2 sigma_33 + S1 (sigma_11 + sigma_22 + sigma_33)``

with ``sigma_phi = sigma_11 cos^2(phi) + sigma_12 sin(2 phi) + sigma_22 sin^2(phi)``
and ``tau_phi = sigma_13 cos(phi) + sigma_23 sin(phi)``. For an isotropic
material ``1/2 S2 = (1 + nu)/E`` and ``S1 = -nu/E``.

Three consequences are the whole method:

- Under plane stress (``sigma_i3 = 0``) ``d`` is **linear in sin^2(psi)**, and the
  slope is ``d0 1/2 S2 sigma_phi``. The stress is read from a slope, so a
  small error in ``d0`` scales it by ``1 + delta d0/d0`` and nothing more -- the
  reason the method works on real materials whose ``d0`` is not known to
  better than a part in ``10^4``.
- A shear stress ``tau_phi`` adds a term **odd in psi**: ``d`` against
  ``sin^2(psi)`` splits into two branches for ``psi > 0`` and ``psi < 0``
  (psi-splitting). Measuring both signs of ``psi`` exposes it; averaging the
  branches removes it.
- Three azimuths not related by 180 degrees determine the three in-plane
  components, because ``sigma_phi`` is a quadratic form in ``(cos(phi), sin(phi))``.

What the module offers
----------------------
:class:`DiffractionElasticConstants`
    ``S1`` and ``1/2 S2`` from ``(E, nu)``, or from single-crystal stiffness
    under the Reuss, Voigt, Neerfeld-Hill or Kroener grain-interaction model,
    for any crystal system, by averaging the grain response about the plane
    normal.
:func:`locate_stress_peak`
    One peak position per tilt, by a pseudo-Voigt doublet fit (default), a
    parabola through the peak top, or a centroid, with an optional
    Lorentz-polarization-absorption correction, and a standard uncertainty.
:func:`fit_sin2psi_lines`
    The classical picture: ``d`` against ``sin^2(psi)`` per azimuth, with slope,
    psi-splitting term, a curvature test and ``sigma_phi``.
:func:`determine_residual_stress`
    All tilts and azimuths in one weighted linear least-squares for the stress
    tensor (biaxial, biaxial with shear, or triaxial), with a full uncertainty
    budget: profile-fit statistics, scatter, ``d0``, the elastic constants, and
    a Monte Carlo cross-check of the linear propagation.
:func:`simulate_sin2psi_measurement`
    A measurement with a *known* stress, for learning, method development and
    the tests.

Standard uncertainty is written ``u(x)`` throughout, following the GUM, because
``sigma`` is the stress.

What is deliberately not here
-----------------------------
Stress gradients (the ``tau``-method and the scattering-vector method),
textured materials (the crystallite-group and ``f_ij``-stress-factor methods)
and grain-interaction models beyond Kroener (Vook-Witt, Eshelby-Kroener for
non-spherical grains). Each shows up in the diagnostics of this module as
curvature or oscillation of ``d`` against ``sin^2(psi)``, and the report says so;
none is silently fitted around.

References
----------
Macherauch, E. & Mueller, P., *Z. angew. Phys.* **13** (1961) 305-312 -- the
sin^2(psi) method in its standard form.

Noyan, I. C. & Cohen, J. B., *Residual Stress: Measurement by Diffraction and
Interpretation*, Springer (1987), doi:10.1007/978-1-4613-9570-6.

Hauk, V. (ed.), *Structural and Residual Stress Analysis by Nondestructive
Methods*, Elsevier (1997), doi:10.1016/B978-0-444-82476-9.X5000-2.

Doelle, H., *J. Appl. Crystallogr.* **12** (1979) 489-501,
doi:10.1107/S0021889879013169 -- psi-splitting and the triaxial evaluation.

Welzel, U., Ligot, J., Lamparter, P., Vermeulen, A. C. & Mittemeijer, E. J.,
*J. Appl. Crystallogr.* **38** (2005) 1-29, doi:10.1107/S0021889804029516 --
diffraction elastic constants and grain-interaction models.

Kroener, E., *Z. Physik* **151** (1958) 504-518, doi:10.1007/BF01337948 --
the self-consistent grain-interaction model.

Reuss, A., *Z. angew. Math. Mech.* **9** (1929) 49-58,
doi:10.1002/zamm.19290090104.

Hill, R., *Proc. Phys. Soc. A* **65** (1952) 349-354,
doi:10.1088/0370-1298/65/5/307.

Eshelby, J. D., *Proc. R. Soc. Lond. A* **241** (1957) 376-396,
doi:10.1098/rspa.1957.0133 -- the constraint tensor of a spherical inclusion.

Fitzpatrick, M. E. et al., *Determination of Residual Stresses by X-ray
Diffraction*, NPL Measurement Good Practice Guide No. 52, issue 2 (2005).

SAE HS-784, *Residual Stress Measurement by X-Ray Diffraction* (2003) -- the
Lorentz-polarization-absorption correction and peak-location practice.

JCGM 100:2008, *Evaluation of measurement data -- Guide to the expression of
uncertainty in measurement* -- the uncertainty budget.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Literal

import numpy as np

from pytex.core.lattice import Phase
from pytex.diffraction.xrd import RadiationSpec
from pytex.diffraction.xrd_corrections import strip_kalpha2
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import (
    PeakFit,
    fit_peaks,
    kalpha_doublet_parameters,
    pseudo_voigt_profile,
)
from pytex.properties.tensors import StiffnessTensor

__all__ = [
    "DEC_MODELS",
    "PEAK_LOCATION_METHODS",
    "RESIDUAL_STRESS_SCHEMA",
    "SINGLE_CRYSTAL_STIFFNESS_GPA",
    "STRESS_COMPONENTS",
    "STRESS_STATES",
    "TILT_GEOMETRIES",
    "DecModel",
    "DiffractionElasticConstants",
    "PeakLocationMethod",
    "ResidualStressResult",
    "Sin2PsiMeasurement",
    "Sin2PsiRegression",
    "StressPeak",
    "StressScan",
    "StressState",
    "StressTensorFit",
    "TiltGeometry",
    "determine_residual_stress",
    "fit_sin2psi_lines",
    "kroener_shear_modulus_cubic",
    "locate_stress_peak",
    "lpa_factor",
    "measurement_direction",
    "parse_stress_peak_positions",
    "parse_stress_scans",
    "residual_stress_pipeline",
    "simulate_sin2psi_measurement",
    "single_crystal_stiffness",
    "strain_design_matrix",
]

RESIDUAL_STRESS_SCHEMA = "pytex.diffraction.residual_stress_result"

StressState = Literal["biaxial", "biaxial_shear", "triaxial"]
DecModel = Literal["isotropic", "reuss", "voigt", "hill", "kroener", "user"]
PeakLocationMethod = Literal["pseudo_voigt", "parabola", "centroid", "given"]
TiltGeometry = Literal["omega", "chi"]

STRESS_STATES: tuple[StressState, ...] = ("biaxial", "biaxial_shear", "triaxial")
DEC_MODELS: tuple[DecModel, ...] = ("isotropic", "reuss", "voigt", "hill", "kroener", "user")
PEAK_LOCATION_METHODS: tuple[PeakLocationMethod, ...] = (
    "pseudo_voigt",
    "parabola",
    "centroid",
    "given",
)
TILT_GEOMETRIES: tuple[TiltGeometry, ...] = ("omega", "chi")

#: The six independent stress components, in the order every design matrix,
#: covariance and table of this module uses.
STRESS_COMPONENTS: tuple[str, ...] = (
    "sigma_11",
    "sigma_22",
    "sigma_12",
    "sigma_13",
    "sigma_23",
    "sigma_33",
)
_COMPONENT_INDEX: dict[str, tuple[int, int]] = {
    "sigma_11": (0, 0),
    "sigma_22": (1, 1),
    "sigma_12": (0, 1),
    "sigma_13": (0, 2),
    "sigma_23": (1, 2),
    "sigma_33": (2, 2),
}
_STATE_COMPONENTS: dict[str, tuple[str, ...]] = {
    "biaxial": STRESS_COMPONENTS[:3],
    "biaxial_shear": STRESS_COMPONENTS[:5],
    "triaxial": STRESS_COMPONENTS,
}

#: One per 10^6: a strain in microstrain times a stress in MPa times a
#: compliance in 1/TPa all meet at this factor.
_PER_TPA_TO_PER_MPA = 1.0e-6

#: Single-crystal stiffness constants (GPa) for the materials residual stress is
#: most often measured on, at room temperature. Cubic entries give
#: ``(C11, C12, C44)``; the hexagonal entry ``(C11, C12, C13, C33, C44)``.
#: Sources: Simmons & Wang, *Single Crystal Elastic Constants and Calculated
#: Aggregate Properties*, MIT Press (1971) for the cubic metals; Fisher &
#: Renken, Phys. Rev. 135 (1964) A482, doi:10.1103/PhysRev.135.A482 for
#: alpha-titanium.
SINGLE_CRYSTAL_STIFFNESS_GPA: Mapping[str, tuple[str, tuple[float, ...]]] = MappingProxyType(
    {
        "fe_bcc": ("cubic", (231.4, 134.7, 116.4)),
        "ni_fcc": ("cubic", (246.5, 147.3, 124.7)),
        "al_fcc": ("cubic", (107.3, 60.9, 28.3)),
        "cu_fcc": ("cubic", (168.4, 121.4, 75.4)),
        "w_bcc": ("cubic", (522.4, 204.4, 160.6)),
        "ti_hcp": ("hexagonal", (162.4, 92.0, 69.0, 180.7, 46.7)),
    }
)

_CITATION_MACHERAUCH = "Macherauch & Mueller, Z. angew. Phys. 13 (1961) 305."
_CITATION_NOYAN_COHEN = (
    "Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6."
)
_CITATION_DOELLE = "Doelle, J. Appl. Cryst. 12 (1979) 489, doi:10.1107/S0021889879013169."
_CITATION_WELZEL = "Welzel et al., J. Appl. Cryst. 38 (2005) 1, doi:10.1107/S0021889804029516."
_CITATION_KROENER = "Kroener, Z. Physik 151 (1958) 504, doi:10.1007/BF01337948."
_CITATION_GUM = "JCGM 100:2008, Guide to the expression of uncertainty in measurement."


def single_crystal_stiffness(material: str) -> StiffnessTensor:
    """Return the tabulated single-crystal stiffness of a named material.

    Parameters
    ----------
    material
        A key of :data:`SINGLE_CRYSTAL_STIFFNESS_GPA`, e.g. ``"fe_bcc"``.

    Returns
    -------
    StiffnessTensor
        In GPa, in the crystal Cartesian frame of
        :meth:`pytex.core.lattice.Lattice.direct_basis`.

    Raises
    ------
    KeyError
        For a material that is not tabulated; the message lists those that are.
    """

    try:
        system, constants = SINGLE_CRYSTAL_STIFFNESS_GPA[material]
    except KeyError:
        known = ", ".join(sorted(SINGLE_CRYSTAL_STIFFNESS_GPA))
        raise KeyError(f"No tabulated stiffness for {material!r}; known: {known}.") from None
    if system == "cubic":
        return StiffnessTensor.cubic(*constants)
    return StiffnessTensor.hexagonal(*constants)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def measurement_direction(phi_deg: Any, psi_deg: Any) -> np.ndarray:
    """Return the specimen-frame scattering-vector directions for ``(phi, psi)``.

    Purpose
    -------
    Fix the one convention every other function of the module depends on: the
    azimuth ``phi`` is measured in the specimen surface from ``S1`` towards
    ``S2``, and the tilt ``psi`` from the surface normal ``S3``.

    Parameters
    ----------
    phi_deg, psi_deg
        Azimuths and tilts in degrees; broadcast against each other.

    Returns
    -------
    numpy.ndarray
        Unit vectors ``(cos phi sin psi, sin phi sin psi, cos psi)``, shape
        ``(..., 3)``.
    """

    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))
    psi = np.deg2rad(np.asarray(psi_deg, dtype=float))
    phi, psi = np.broadcast_arrays(phi, psi)
    return np.stack((np.cos(phi) * np.sin(psi), np.sin(phi) * np.sin(psi), np.cos(psi)), axis=-1)


def _projection_coefficients(phi_deg: np.ndarray, psi_deg: np.ndarray) -> np.ndarray:
    """Return ``m.sigma.m`` coefficients of the six components, shape ``(n, 6)``."""

    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))
    psi = np.deg2rad(np.asarray(psi_deg, dtype=float))
    sin2 = np.sin(psi) ** 2
    return np.stack(
        (
            np.cos(phi) ** 2 * sin2,
            np.sin(phi) ** 2 * sin2,
            np.sin(2.0 * phi) * sin2,
            np.cos(phi) * np.sin(2.0 * psi),
            np.sin(phi) * np.sin(2.0 * psi),
            np.cos(psi) ** 2,
        ),
        axis=-1,
    )


_TRACE_ROW = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0])


def strain_design_matrix(
    phi_deg: Any,
    psi_deg: Any,
    *,
    s1_per_tpa: float,
    half_s2_per_tpa: float,
    components: Sequence[str] = STRESS_COMPONENTS,
) -> np.ndarray:
    """Return the matrix that maps stress components (MPa) to measured strain.

    Purpose
    -------
    The fundamental equation ``epsilon = 1/2 S2 m.sigma.m + S1 tr(sigma)`` is
    linear in the stress components, so every measurement is one row of a
    design matrix and the stress tensor is a linear least-squares solution:
    no starting guess, no local minima, and an analytic covariance.

    Parameters
    ----------
    phi_deg, psi_deg
        Azimuth and tilt of each measurement, in degrees.
    s1_per_tpa, half_s2_per_tpa
        Diffraction elastic constants in ``1/TPa`` (``10^-6 / MPa``).
    components
        The free components, a subset of :data:`STRESS_COMPONENTS`, in the
        column order wanted.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, len(components))``; the dimensionless strain produced by
        one MPa of each component.
    """

    phi = np.atleast_1d(np.asarray(phi_deg, dtype=float))
    psi = np.atleast_1d(np.asarray(psi_deg, dtype=float))
    coefficients = _projection_coefficients(phi, psi)
    full = (half_s2_per_tpa * coefficients + s1_per_tpa * _TRACE_ROW) * _PER_TPA_TO_PER_MPA
    try:
        columns = [STRESS_COMPONENTS.index(name) for name in components]
    except ValueError as error:
        raise ValueError(f"Unknown stress component in {tuple(components)}.") from error
    return np.ascontiguousarray(full[:, columns])


# ---------------------------------------------------------------------------
# Diffraction elastic constants
# ---------------------------------------------------------------------------

_MANDEL_PAIRS: tuple[tuple[int, int], ...] = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
_MANDEL_WEIGHTS = np.array([1.0, 1.0, 1.0, math.sqrt(2.0), math.sqrt(2.0), math.sqrt(2.0)])
_MANDEL_J = np.zeros((6, 6))
_MANDEL_J[:3, :3] = 1.0 / 3.0
_MANDEL_K = np.eye(6) - _MANDEL_J


def _to_mandel(tensor: np.ndarray) -> np.ndarray:
    """Return the 6x6 Mandel (orthonormal) form of a minor-symmetric rank-4 tensor.

    Mandel form, not Voigt form, because it is the one in which the inverse of
    the matrix is the matrix of the inverse tensor for stiffness and compliance
    alike -- no factors of two to remember, which is where hand-rolled
    micromechanics usually goes wrong.
    """

    matrix = np.empty((6, 6))
    for row, (i, j) in enumerate(_MANDEL_PAIRS):
        for column, (k, m) in enumerate(_MANDEL_PAIRS):
            matrix[row, column] = tensor[i, j, k, m]
    return matrix * np.outer(_MANDEL_WEIGHTS, _MANDEL_WEIGHTS)


def _mandel_dyad(vectors: np.ndarray) -> np.ndarray:
    """Return the Mandel vectors of ``a a`` for unit vectors ``a``, shape ``(..., 6)``."""

    a = np.asarray(vectors, dtype=float)
    return np.stack(
        (
            a[..., 0] ** 2,
            a[..., 1] ** 2,
            a[..., 2] ** 2,
            math.sqrt(2.0) * a[..., 1] * a[..., 2],
            math.sqrt(2.0) * a[..., 0] * a[..., 2],
            math.sqrt(2.0) * a[..., 0] * a[..., 1],
        ),
        axis=-1,
    )


def _isotropic_moduli(mandel: np.ndarray) -> tuple[float, float]:
    """Return ``(3K, 2G)``-type coefficients of the isotropic projection.

    Any rank-4 tensor averaged over all orientations becomes ``a J + b K``,
    with ``J = (1/3) 1 (x) 1`` the spherical and ``K = I - J`` the deviatoric
    projector. ``a = J::T`` and ``b = K::T / 5``. For a stiffness ``a = 3K``
    and ``b = 2G``; for a compliance ``a = 1/(3K)`` and ``b = 1/(2G)``.
    """

    a = float(np.sum(_MANDEL_J * mandel))
    b = float(np.sum(_MANDEL_K * mandel) / 5.0)
    return a, b


def _isotropic_mandel(bulk: float, shear: float) -> np.ndarray:
    return 3.0 * bulk * _MANDEL_J + 2.0 * shear * _MANDEL_K


def kroener_shear_modulus_cubic(c11: float, c12: float, c44: float) -> float:
    """Return the Kroener self-consistent shear modulus of a cubic aggregate.

    Purpose
    -------
    The closed-form check of the general self-consistent iteration used by
    :meth:`DiffractionElasticConstants.from_single_crystal`: for cubic crystals
    the Kroener shear modulus ``G`` is the positive root of

    ``G^3 + alpha G^2 + beta G + gamma = 0``

    with ``alpha = (5 C11 + 4 C12)/8``, ``beta = -C44 (7 C11 - 4 C12)/8`` and
    ``gamma = -C44 (C11 - C12)(C11 + 2 C12)/8``.

    Parameters
    ----------
    c11, c12, c44
        Single-crystal stiffness constants in any consistent unit.

    Returns
    -------
    float
        ``G`` in the unit of the constants.
    """

    alpha = (5.0 * c11 + 4.0 * c12) / 8.0
    beta = -c44 * (7.0 * c11 - 4.0 * c12) / 8.0
    gamma = -c44 * (c11 - c12) * (c11 + 2.0 * c12) / 8.0
    roots = np.roots([1.0, alpha, beta, gamma])
    real = [float(root.real) for root in roots if abs(root.imag) < 1e-9 and root.real > 0.0]
    if len(real) != 1:
        raise ValueError("The Kroener cubic has no unique positive root for these constants.")
    return real[0]


def _kroener_grain_compliance(stiffness_mandel: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Return the Kroener grain strain-per-macro-stress map and the aggregate moduli.

    The self-consistent scheme: a spherical grain of stiffness ``C_g`` in an
    effective isotropic matrix ``C_eff`` strains by ``A_g : epsilon_macro`` with
    ``A_g = (C_g + C*)^-1 : (C_eff + C*)``, where ``C*`` is Hill's constraint
    tensor of a sphere, ``K* = 4G/3`` and ``G* = G (9K + 8G) / (6 (K + 2G))``.
    ``C_eff`` is the orientation average of ``C_g : A_g``, iterated to a fixed
    point from the Hill average. Returns ``A_g : S_eff``.
    """

    compliance = np.linalg.inv(stiffness_mandel)
    voigt_a, voigt_b = _isotropic_moduli(stiffness_mandel)
    reuss_a, reuss_b = _isotropic_moduli(compliance)
    bulk = 0.5 * (voigt_a / 3.0 + 1.0 / (3.0 * reuss_a))
    shear = 0.5 * (voigt_b / 2.0 + 1.0 / (2.0 * reuss_b))
    concentration = np.eye(6)
    for _ in range(500):
        bulk_star = 4.0 * shear / 3.0
        shear_star = shear * (9.0 * bulk + 8.0 * shear) / (6.0 * (bulk + 2.0 * shear))
        constraint = _isotropic_mandel(bulk_star, shear_star)
        effective = _isotropic_mandel(bulk, shear)
        concentration = np.linalg.solve(stiffness_mandel + constraint, effective + constraint)
        new_a, new_b = _isotropic_moduli(stiffness_mandel @ concentration)
        new_bulk, new_shear = new_a / 3.0, new_b / 2.0
        converged = abs(new_bulk - bulk) <= 1e-13 * abs(bulk) and abs(
            new_shear - shear
        ) <= 1e-13 * abs(shear)
        bulk, shear = new_bulk, new_shear
        if converged:
            break
    effective = _isotropic_mandel(bulk, shear)
    bulk_star = 4.0 * shear / 3.0
    shear_star = shear * (9.0 * bulk + 8.0 * shear) / (6.0 * (bulk + 2.0 * shear))
    constraint = _isotropic_mandel(bulk_star, shear_star)
    concentration = np.linalg.solve(stiffness_mandel + constraint, effective + constraint)
    grain: np.ndarray = concentration @ np.linalg.inv(effective)
    return grain, bulk, shear


def _decs_from_grain_map(grain_map: np.ndarray, normal: np.ndarray) -> tuple[float, float]:
    """Return ``(S1, 1/2 S2)`` of a grain strain-per-stress map about a plane normal.

    The grains that diffract have the plane normal ``n`` along the scattering
    vector and every rotation ``omega`` about it. Averaging over ``omega``
    makes the response transversely isotropic about ``n``, so two numbers
    describe it: the strain along ``n`` per unit stress along a perpendicular
    ``t(omega)`` is ``S1``, and per unit stress along ``n`` it is
    ``S1 + 1/2 S2``.
    """

    n = normal / np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    first = np.cross(n, helper)
    first /= np.linalg.norm(first)
    second = np.cross(n, first)
    # 360 equally spaced angles integrate every trigonometric polynomial of
    # degree below 360 exactly; the integrand here is of degree 4 in omega.
    omega = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
    perpendicular = np.cos(omega)[:, None] * first + np.sin(omega)[:, None] * second
    along = _mandel_dyad(n)
    across = _mandel_dyad(perpendicular)
    s1 = float(np.mean(across @ grain_map.T @ along))
    total = float(along @ grain_map @ along)
    return s1, total - s1


@dataclass(frozen=True, slots=True)
class DiffractionElasticConstants:
    """The two constants that turn a lattice strain into a stress.

    Purpose
    -------
    Carry ``S1(hkl)`` and ``1/2 S2(hkl)`` with the model they came from, so a
    reported stress always says which elastic constants it rests on. They are
    reflection-specific for an elastically anisotropic crystal: in ferrite
    ``1/2 S2`` of ``{200}`` is about 40 % larger than that of ``{222}``
    under the Reuss model, and a stress computed with the wrong one is wrong by
    the same factor.

    Attributes
    ----------
    s1_per_tpa : float
        ``S1`` in ``1/TPa`` (``10^-6 / MPa``). Negative for every real
        material: it is the Poisson term ``-nu/E`` of an isotropic solid.
    half_s2_per_tpa : float
        ``1/2 S2`` in ``1/TPa``; ``(1 + nu)/E`` for an isotropic solid. Strictly
        positive.
    model : str
        One of :data:`DEC_MODELS`.
    hkl : tuple of int, optional
        The reflection the constants belong to.
    source : str
        Where the constants came from, in words: the material, the single-crystal
        constants, or "user supplied".
    relative_standard_uncertainty : float
        ``u(S)/S``, applied to both constants independently in the uncertainty
        budget. Zero leaves the elastic constants out of the budget, which is
        a claim that they are exact and should be made knowingly.
    """

    s1_per_tpa: float
    half_s2_per_tpa: float
    model: DecModel = "user"
    hkl: tuple[int, int, int] | None = None
    source: str = "user supplied"
    relative_standard_uncertainty: float = 0.0

    def __post_init__(self) -> None:
        if self.model not in DEC_MODELS:
            raise ValueError(f"DiffractionElasticConstants.model must be one of {DEC_MODELS}.")
        if not (math.isfinite(self.s1_per_tpa) and math.isfinite(self.half_s2_per_tpa)):
            raise ValueError("Diffraction elastic constants must be finite.")
        if self.half_s2_per_tpa <= 0.0:
            raise ValueError("1/2 S2 must be strictly positive (it is 1/(2 G) of the grain).")
        if self.half_s2_per_tpa + 3.0 * self.s1_per_tpa <= 0.0:
            raise ValueError(
                "1/2 S2 + 3 S1 must be positive: it is the compressibility 1/(3K) of the "
                "aggregate, and a negative value describes a material that expands under "
                "pressure."
            )
        if not (
            math.isfinite(self.relative_standard_uncertainty)
            and self.relative_standard_uncertainty >= 0.0
        ):
            raise ValueError("relative_standard_uncertainty must be finite and non-negative.")
        if self.hkl is not None:
            object.__setattr__(self, "hkl", tuple(int(value) for value in self.hkl))

    @classmethod
    def isotropic(
        cls,
        youngs_modulus_gpa: float,
        poisson_ratio: float,
        *,
        hkl: tuple[int, int, int] | None = None,
        relative_standard_uncertainty: float = 0.0,
    ) -> DiffractionElasticConstants:
        """Constants of an elastically isotropic material.

        ``S1 = -nu/E`` and ``1/2 S2 = (1 + nu)/E``: the same for every
        reflection. Adequate for tungsten, which is nearly isotropic, and a
        known approximation for everything else; the Reuss-to-Voigt spread of
        :meth:`from_single_crystal` says by how much.
        """

        if youngs_modulus_gpa <= 0.0:
            raise ValueError("Young's modulus must be positive.")
        if not -1.0 < poisson_ratio < 0.5:
            raise ValueError("Poisson's ratio must lie in (-1, 0.5).")
        scale = 1000.0 / youngs_modulus_gpa
        return cls(
            s1_per_tpa=-poisson_ratio * scale,
            half_s2_per_tpa=(1.0 + poisson_ratio) * scale,
            model="isotropic",
            hkl=hkl,
            source=f"isotropic, E = {youngs_modulus_gpa:g} GPa, nu = {poisson_ratio:g}",
            relative_standard_uncertainty=relative_standard_uncertainty,
        )

    @classmethod
    def from_single_crystal(
        cls,
        stiffness: StiffnessTensor,
        normal: Any,
        *,
        model: Literal["reuss", "voigt", "hill", "kroener"] = "kroener",
        hkl: tuple[int, int, int] | None = None,
        source: str = "single-crystal stiffness",
        relative_standard_uncertainty: float = 0.0,
    ) -> DiffractionElasticConstants:
        """Constants of one reflection of an untextured polycrystal.

        Purpose
        -------
        Compute ``S1(hkl)`` and ``1/2 S2(hkl)`` from single-crystal stiffness
        under a stated grain-interaction model, for any crystal system.

        Method
        ------
        Each model gives the strain of a grain per unit macroscopic stress, a
        rank-4 map ``P_g`` in the crystal frame:

        - **Reuss**: every grain carries the macroscopic stress, ``P_g = S_g``.
          Reflection dependent; the lower bound of stiffness.
        - **Voigt**: every grain carries the macroscopic strain, ``P_g`` is the
          compliance of the Voigt-averaged aggregate. Reflection independent.
        - **Neerfeld-Hill**: the arithmetic mean of the Reuss and Voigt
          constants.
        - **Kroener**: each grain is an Eshelby sphere in the self-consistent
          effective medium; see :func:`_kroener_grain_compliance`. Usually the
          closest of the four to measurement for untextured cubic metals.

        The diffracting grains have the plane normal along the scattering
        vector and every rotation about it; ``P_g`` is averaged over that
        rotation, which leaves the two constants.

        Parameters
        ----------
        stiffness
            Single-crystal stiffness in GPa, in the crystal Cartesian frame.
        normal
            Plane normal of the reflection in that frame, e.g.
            :attr:`pytex.core.lattice.CrystalPlane.normal`. Need not be unit.
        model
            ``"reuss"``, ``"voigt"``, ``"hill"`` or ``"kroener"``.
        hkl, source, relative_standard_uncertainty
            Recorded on the result.
        """

        direction = np.asarray(normal, dtype=float).reshape(3)
        if not np.all(np.isfinite(direction)) or np.linalg.norm(direction) == 0.0:
            raise ValueError("The plane normal must be a finite, non-zero vector.")
        stiffness_mandel = _to_mandel(np.asarray(stiffness.tensor, dtype=float))
        compliance_mandel = np.linalg.inv(stiffness_mandel)
        if model == "reuss":
            s1, half_s2 = _decs_from_grain_map(compliance_mandel, direction)
        elif model == "voigt":
            a, b = _isotropic_moduli(stiffness_mandel)
            s1, half_s2 = _decs_from_grain_map(
                np.linalg.inv(_isotropic_mandel(a / 3.0, b / 2.0)), direction
            )
        elif model == "hill":
            reuss = _decs_from_grain_map(compliance_mandel, direction)
            a, b = _isotropic_moduli(stiffness_mandel)
            voigt = _decs_from_grain_map(
                np.linalg.inv(_isotropic_mandel(a / 3.0, b / 2.0)), direction
            )
            s1, half_s2 = 0.5 * (reuss[0] + voigt[0]), 0.5 * (reuss[1] + voigt[1])
        elif model == "kroener":
            grain, _, _ = _kroener_grain_compliance(stiffness_mandel)
            s1, half_s2 = _decs_from_grain_map(grain, direction)
        else:
            raise ValueError("model must be 'reuss', 'voigt', 'hill' or 'kroener'.")
        return cls(
            s1_per_tpa=1000.0 * s1,
            half_s2_per_tpa=1000.0 * half_s2,
            model=model,
            hkl=hkl,
            source=source,
            relative_standard_uncertainty=relative_standard_uncertainty,
        )

    @classmethod
    def for_reflection(
        cls,
        phase: Phase,
        hkl: tuple[int, int, int],
        stiffness: StiffnessTensor,
        *,
        model: Literal["reuss", "voigt", "hill", "kroener"] = "kroener",
        source: str = "single-crystal stiffness",
        relative_standard_uncertainty: float = 0.0,
    ) -> DiffractionElasticConstants:
        """:meth:`from_single_crystal` with the plane normal resolved from ``(hkl)``.

        The normal goes through the reciprocal basis of ``phase``, which is the
        only correct route outside the cubic system.
        """

        reciprocal = phase.lattice.reciprocal_basis().matrix
        normal = reciprocal @ np.asarray(hkl, dtype=float)
        return cls.from_single_crystal(
            stiffness,
            normal,
            model=model,
            hkl=hkl,
            source=source,
            relative_standard_uncertainty=relative_standard_uncertainty,
        )

    @property
    def youngs_modulus_gpa(self) -> float:
        """The reflection's effective Young's modulus ``1/(S1 + 1/2 S2)``, in GPa."""

        return 1000.0 / (self.s1_per_tpa + self.half_s2_per_tpa)

    @property
    def poisson_ratio(self) -> float:
        """The reflection's effective Poisson ratio ``-S1/(S1 + 1/2 S2)``."""

        return -self.s1_per_tpa / (self.s1_per_tpa + self.half_s2_per_tpa)

    def with_uncertainty(self, relative_standard_uncertainty: float) -> DiffractionElasticConstants:
        """Return a copy carrying a different relative standard uncertainty."""

        return DiffractionElasticConstants(
            s1_per_tpa=self.s1_per_tpa,
            half_s2_per_tpa=self.half_s2_per_tpa,
            model=self.model,
            hkl=self.hkl,
            source=self.source,
            relative_standard_uncertainty=relative_standard_uncertainty,
        )

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable form."""

        return {
            "s1_per_tpa": self.s1_per_tpa,
            "half_s2_per_tpa": self.half_s2_per_tpa,
            "model": self.model,
            "hkl": None if self.hkl is None else list(self.hkl),
            "source": self.source,
            "relative_standard_uncertainty": self.relative_standard_uncertainty,
            "youngs_modulus_gpa": self.youngs_modulus_gpa,
            "poisson_ratio": self.poisson_ratio,
        }

    def describe(self) -> str:
        """Return the constants as prose, with their model and meaning."""

        names = {
            "isotropic": "isotropic elasticity",
            "reuss": "the Reuss (uniform stress) grain-interaction model",
            "voigt": "the Voigt (uniform strain) grain-interaction model",
            "hill": "the Neerfeld-Hill average of the Reuss and Voigt models",
            "kroener": "the Kroener self-consistent grain-interaction model",
            "user": "user-supplied values",
        }
        reflection = (
            "" if self.hkl is None else f" of the ({''.join(map(str, self.hkl))}) reflection"
        )
        uncertainty = (
            f" A relative standard uncertainty of {100 * self.relative_standard_uncertainty:.1f} %"
            " is carried into the stress uncertainty budget."
            if self.relative_standard_uncertainty > 0.0
            else " No uncertainty is attached to them, so the budget treats them as exact."
        )
        return (
            f"Diffraction elastic constants{reflection}: S1 = {self.s1_per_tpa:.4f} /TPa and "
            f"1/2 S2 = {self.half_s2_per_tpa:.4f} /TPa, from {names[self.model]} "
            f"({self.source}). They correspond to an effective Young's modulus of "
            f"{self.youngs_modulus_gpa:.1f} GPa and Poisson ratio {self.poisson_ratio:.3f} for "
            "this reflection." + uncertainty + " " + _CITATION_WELZEL
        )


# ---------------------------------------------------------------------------
# Measurement containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StressScan:
    """One scan of the stress reflection at one specimen orientation.

    Attributes
    ----------
    phi_deg : float
        Azimuth in the specimen surface, from ``S1`` towards ``S2``.
    psi_deg : float
        Tilt of the scattering vector from the surface normal. Its *sign*
        matters when shear stresses are present; see :func:`fit_sin2psi_lines`.
    pattern : MeasuredPowderPattern
        The measured window around the reflection.
    """

    phi_deg: float
    psi_deg: float
    pattern: MeasuredPowderPattern

    def __post_init__(self) -> None:
        if not math.isfinite(self.phi_deg):
            raise ValueError("StressScan.phi_deg must be finite.")
        if not -90.0 < self.psi_deg < 90.0:
            raise ValueError("StressScan.psi_deg must lie strictly inside (-90, 90) degrees.")
        object.__setattr__(self, "phi_deg", float(self.phi_deg))
        object.__setattr__(self, "psi_deg", float(self.psi_deg))


@dataclass(frozen=True, slots=True)
class Sin2PsiMeasurement:
    """A set of scans of one reflection at several ``(phi, psi)``.

    Attributes
    ----------
    name : str
        What was measured.
    scans : tuple of StressScan
        At least three, and at least two distinct tilts.
    radiation : RadiationSpec
        The radiation every scan was taken with.
    geometry : str
        ``"omega"`` (iso-inclination: tilt in the diffraction plane) or
        ``"chi"`` (side-inclination: tilt about the line of the diffraction
        plane and the surface). It changes the absorption factor across the
        profile, and nothing else here.
    synthetic : bool
        Whether the scans were generated rather than measured.
    true_stress_mpa : numpy.ndarray, optional
        The stress a synthetic measurement was generated with, so an analysis
        can be judged against a known answer.
    metadata : mapping
        Free-form provenance.
    """

    name: str
    scans: tuple[StressScan, ...]
    radiation: RadiationSpec
    geometry: TiltGeometry = "omega"
    synthetic: bool = False
    true_stress_mpa: np.ndarray | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Sin2PsiMeasurement.name must be non-empty.")
        scans = tuple(self.scans)
        if len(scans) < 3:
            raise ValueError("A sin^2(psi) measurement needs at least three scans.")
        if len({round(scan.psi_deg, 6) for scan in scans}) < 2:
            raise ValueError("A sin^2(psi) measurement needs at least two distinct tilts.")
        if self.geometry not in TILT_GEOMETRIES:
            raise ValueError(f"geometry must be one of {TILT_GEOMETRIES}.")
        object.__setattr__(self, "scans", scans)
        if self.true_stress_mpa is not None:
            stress = np.asarray(self.true_stress_mpa, dtype=float)
            if stress.shape != (3, 3) or not np.allclose(stress, stress.T):
                raise ValueError("true_stress_mpa must be a symmetric 3x3 tensor.")
            stress = stress.copy()
            stress.setflags(write=False)
            object.__setattr__(self, "true_stress_mpa", stress)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def __len__(self) -> int:
        return len(self.scans)

    @property
    def phi_deg(self) -> np.ndarray:
        """Azimuth of each scan, in degrees."""

        return np.array([scan.phi_deg for scan in self.scans])

    @property
    def psi_deg(self) -> np.ndarray:
        """Tilt of each scan, in degrees."""

        return np.array([scan.psi_deg for scan in self.scans])

    def describe(self) -> str:
        """Return the measurement plan as prose."""

        azimuths = sorted({round(scan.phi_deg % 360.0, 6) for scan in self.scans})
        tilts = sorted({round(scan.psi_deg, 6) for scan in self.scans})
        kind = "Synthetic" if self.synthetic else "Measured"
        return (
            f"{kind} sin^2(psi) data set '{self.name}': {len(self.scans)} scans with "
            f"{self.radiation.name}, {len(azimuths)} azimuth(s) phi = "
            f"{', '.join(f'{value:g}' for value in azimuths)} degrees and {len(tilts)} tilt(s) "
            f"from {tilts[0]:g} to {tilts[-1]:g} degrees ({self.geometry} tilting)."
        )


# ---------------------------------------------------------------------------
# Peak location
# ---------------------------------------------------------------------------


def lpa_factor(
    two_theta_deg: Any, *, psi_deg: float, geometry: TiltGeometry = "omega"
) -> np.ndarray:
    """Return the Lorentz-polarization-absorption factor across a profile.

    Purpose
    -------
    A stress peak is broad -- often several degrees on a hardened steel -- and
    across that width the Lorentz-polarization factor and, under omega tilting,
    the absorption factor change enough to skew the profile and move its
    fitted centre by an amount that depends on ``psi``. Dividing the measured
    intensity by this factor before locating the peak removes the skew.

    Method
    ------
    ``LP = (1 + cos^2(2 theta)) / sin^2(theta)``, the form used for broad
    stress profiles (SAE HS-784). Under omega tilting the absorption factor of
    an infinitely thick specimen is ``A = 1 - tan(psi) cot(theta)``; under chi
    tilting it does not vary across the profile and is omitted.

    Parameters
    ----------
    two_theta_deg
        Angles across the profile.
    psi_deg
        The tilt of the scan.
    geometry
        ``"omega"`` or ``"chi"``.

    Returns
    -------
    numpy.ndarray
        The factor, positive, at each angle.
    """

    theta = np.deg2rad(0.5 * np.asarray(two_theta_deg, dtype=float))
    factor = (1.0 + np.cos(2.0 * theta) ** 2) / np.sin(theta) ** 2
    if geometry == "omega":
        absorption = 1.0 - np.tan(np.deg2rad(psi_deg)) / np.tan(theta)
        if np.any(absorption <= 0.0):
            raise ValueError(
                "Under omega tilting this tilt takes the incident or diffracted beam below the "
                "specimen surface (1 - tan(psi) cot(theta) <= 0); use chi tilting or a smaller "
                "tilt."
            )
        factor = factor * absorption
    elif geometry != "chi":
        raise ValueError(f"geometry must be one of {TILT_GEOMETRIES}.")
    result: np.ndarray = factor
    return result


@dataclass(frozen=True, slots=True)
class StressPeak:
    """The located position of the stress reflection at one ``(phi, psi)``.

    Attributes
    ----------
    phi_deg, psi_deg : float
        The specimen orientation.
    two_theta_deg : float
        Peak position; the K-alpha1 position when the doublet was modelled.
    two_theta_uncertainty_deg : float
        Standard uncertainty ``u(2 theta)`` of the position.
    method : str
        One of :data:`PEAK_LOCATION_METHODS`; ``"given"`` for a position read
        from a table rather than located here.
    fwhm_deg : float
        Width of the fitted profile, or ``nan`` when the method measures none.
    height : float
        Peak height above background, or ``nan``.
    reduced_chi_squared : float
        Of the profile fit, or ``nan``.
    point_count : int
        Measured points used to locate the peak; zero for ``"given"``.
    converged : bool
        Whether the location succeeded cleanly.
    lpa_corrected : bool
        Whether the intensities were divided by :func:`lpa_factor` first.
    kalpha2_stripped : bool
        Whether the K-alpha2 line was removed by Rachinger's recursion first.
    uncertainty_is_nominal : bool
        The position came with no uncertainty and carries a placeholder; the
        stress uncertainty is then taken from the scatter alone.
    peak_fit : PeakFit, optional
        The profile fit, for ``"pseudo_voigt"``.
    """

    phi_deg: float
    psi_deg: float
    two_theta_deg: float
    two_theta_uncertainty_deg: float
    method: PeakLocationMethod
    fwhm_deg: float = float("nan")
    height: float = float("nan")
    reduced_chi_squared: float = float("nan")
    point_count: int = 0
    converged: bool = True
    lpa_corrected: bool = False
    kalpha2_stripped: bool = False
    uncertainty_is_nominal: bool = False
    peak_fit: PeakFit | None = None

    def __post_init__(self) -> None:
        if self.method not in PEAK_LOCATION_METHODS:
            raise ValueError(f"StressPeak.method must be one of {PEAK_LOCATION_METHODS}.")
        if not 0.0 < self.two_theta_deg < 180.0:
            raise ValueError("StressPeak.two_theta_deg must lie strictly inside (0, 180).")
        if not (
            math.isfinite(self.two_theta_uncertainty_deg) and self.two_theta_uncertainty_deg > 0.0
        ):
            raise ValueError("StressPeak.two_theta_uncertainty_deg must be finite and positive.")
        if not -90.0 < self.psi_deg < 90.0:
            raise ValueError("StressPeak.psi_deg must lie strictly inside (-90, 90) degrees.")

    @property
    def sin2psi(self) -> float:
        """``sin^2(psi)``, the abscissa of the stress plot."""

        return float(np.sin(np.deg2rad(self.psi_deg)) ** 2)

    def d_spacing_angstrom(self, wavelength_angstrom: float) -> float:
        """Return ``d = lambda / (2 sin theta)``."""

        return float(wavelength_angstrom / (2.0 * np.sin(np.deg2rad(0.5 * self.two_theta_deg))))

    def d_spacing_uncertainty_angstrom(self, wavelength_angstrom: float) -> float:
        """Return ``u(d) = d cot(theta) u(theta)``, first-order Bragg propagation."""

        theta = np.deg2rad(0.5 * self.two_theta_deg)
        d = self.d_spacing_angstrom(wavelength_angstrom)
        return float(d / np.tan(theta) * np.deg2rad(0.5 * self.two_theta_uncertainty_deg))

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable form."""

        return {
            "phi_deg": self.phi_deg,
            "psi_deg": self.psi_deg,
            "sin2psi": self.sin2psi,
            "two_theta_deg": self.two_theta_deg,
            "two_theta_uncertainty_deg": self.two_theta_uncertainty_deg,
            "method": self.method,
            "fwhm_deg": None if math.isnan(self.fwhm_deg) else self.fwhm_deg,
            "height": None if math.isnan(self.height) else self.height,
            "reduced_chi_squared": (
                None if math.isnan(self.reduced_chi_squared) else self.reduced_chi_squared
            ),
            "point_count": self.point_count,
            "converged": self.converged,
            "lpa_corrected": self.lpa_corrected,
            "kalpha2_stripped": self.kalpha2_stripped,
            "uncertainty_is_nominal": self.uncertainty_is_nominal,
        }

    def describe(self) -> str:
        """Return the located peak as prose."""

        names = {
            "pseudo_voigt": "a pseudo-Voigt profile fit",
            "parabola": "a parabola through the peak top",
            "centroid": "the centroid above half the net maximum",
            "given": "a supplied position",
        }
        return (
            f"At phi = {self.phi_deg:g} and psi = {self.psi_deg:g} degrees (sin^2(psi) = "
            f"{self.sin2psi:.4f}) the peak lies at 2theta = {self.two_theta_deg:.4f} +/- "
            f"{self.two_theta_uncertainty_deg:.4f} degrees, from {names[self.method]}"
            + (" after the LPA correction." if self.lpa_corrected else ".")
        )


def _window(
    pattern: MeasuredPowderPattern, centre_deg: float | None, width_deg: float | None
) -> tuple[np.ndarray, np.ndarray]:
    axis = np.asarray(pattern.two_theta_deg, dtype=float)
    counts = np.asarray(pattern.intensity, dtype=float)
    if centre_deg is None or width_deg is None:
        return axis, counts
    keep = np.abs(axis - centre_deg) <= 0.5 * width_deg
    if np.count_nonzero(keep) < 7:
        raise ValueError(
            f"Only {np.count_nonzero(keep)} measured points lie within {width_deg:g} degrees of "
            f"{centre_deg:.3f} degrees 2theta; the scan does not cover the reflection."
        )
    return axis[keep], counts[keep]


def _linear_background(axis: np.ndarray, values: np.ndarray, fraction: float = 0.1) -> np.ndarray:
    """Return the straight line through the mean of each end of the window."""

    count = max(2, round(fraction * axis.size))
    x0, y0 = float(np.mean(axis[:count])), float(np.mean(values[:count]))
    x1, y1 = float(np.mean(axis[-count:])), float(np.mean(values[-count:]))
    slope = (y1 - y0) / (x1 - x0) if x1 != x0 else 0.0
    return y0 + slope * (axis - x0)


def locate_stress_peak(
    scan: StressScan,
    *,
    radiation: RadiationSpec | None = None,
    method: PeakLocationMethod = "pseudo_voigt",
    expected_two_theta_deg: float | None = None,
    window_deg: float | None = None,
    expected_fwhm_deg: float = 1.0,
    lpa_correction: bool = True,
    geometry: TiltGeometry = "omega",
    top_fraction: float = 0.8,
    strip_doublet: bool = True,
) -> StressPeak:
    """Locate the stress reflection in one scan, with a standard uncertainty.

    Purpose
    -------
    Turn a measured profile into the one number the stress analysis consumes
    -- the peak position -- together with the uncertainty that weights it.

    Method
    ------
    The window is restricted to ``window_deg`` about ``expected_two_theta_deg``
    when both are given. With ``lpa_correction`` the intensities are divided by
    :func:`lpa_factor` (normalized at the window centre, so the counts keep
    their scale) and their Poisson uncertainties are divided with them. Then:

    ``"pseudo_voigt"``
        :func:`pytex.diffraction.xrd_peaks.fit_peaks` with the K-alpha2
        partner fixed by Bragg's law and a straight local background. The
        position uncertainty is the fit covariance scaled by the window's
        reduced chi-squared. The default, and the only method that models the
        doublet.
    ``"parabola"``
        Weighted least squares of a parabola through the points above
        ``top_fraction`` of the net maximum; the vertex ``-b/(2c)`` and its
        uncertainty from the covariance of ``(b, c)``. The classical method for
        broad peaks. It sees whatever profile it is given, so the K-alpha2 line
        is stripped first (``strip_doublet``) by
        :func:`pytex.diffraction.xrd_corrections.strip_kalpha2`.
    ``"centroid"``
        The first moment of the net intensity above half its maximum, with the
        Poisson propagation of its uncertainty. Robust and biased by any
        asymmetry; offered as a cross-check.

    A constant bias in the position -- a detector zero, a wavelength error --
    moves every tilt alike and leaves the *slope* of ``d`` against
    ``sin^2(psi)`` unchanged to first order, which is why the stress survives
    such a bias and the absolute strain does not.

    Parameters
    ----------
    scan
        The scan to analyse.
    radiation
        Source of the K-alpha1 wavelength and the doublet. Falls back to the
        pattern's own radiation.
    method
        ``"pseudo_voigt"``, ``"parabola"`` or ``"centroid"``.
    expected_two_theta_deg, window_deg
        Where the reflection is expected, and the width of scan to use.
    expected_fwhm_deg
        Starting width for the profile fit.
    lpa_correction, geometry
        Apply :func:`lpa_factor` for this tilt geometry.
    top_fraction
        Fraction of the net maximum above which the parabola is fitted.
    strip_doublet
        For the parabola and the centroid, remove the K-alpha2 line with
        :func:`pytex.diffraction.xrd_corrections.strip_kalpha2` first. The
        profile fit models the doublet instead and ignores this.

    Returns
    -------
    StressPeak

    Raises
    ------
    ValueError
        If the window holds too few points or the location fails outright.
    """

    if method == "given":
        raise ValueError("'given' positions are constructed directly, not located.")
    if method not in PEAK_LOCATION_METHODS:
        raise ValueError(f"method must be one of {PEAK_LOCATION_METHODS}.")
    radiation = radiation or scan.pattern.radiation
    axis, counts = _window(scan.pattern, expected_two_theta_deg, window_deg)
    if scan.pattern.standard_uncertainty is not None:
        _, stated = _window(
            MeasuredPowderPattern(
                name="uncertainty",
                two_theta_deg=scan.pattern.two_theta_deg,
                intensity=scan.pattern.standard_uncertainty,
            ),
            expected_two_theta_deg,
            window_deg,
        )
        uncertainty = np.asarray(stated, dtype=float)
    else:
        uncertainty = np.sqrt(np.maximum(counts, 1.0))
    if lpa_correction:
        factor = lpa_factor(axis, psi_deg=scan.psi_deg, geometry=geometry)
        factor = factor / factor[axis.size // 2]
        counts = counts / factor
        uncertainty = uncertainty / factor

    if method == "pseudo_voigt":
        corrected = MeasuredPowderPattern(
            name=f"phi {scan.phi_deg:g}, psi {scan.psi_deg:g}",
            two_theta_deg=axis,
            intensity=np.maximum(counts, 0.0),
            standard_uncertainty=np.maximum(uncertainty, 1e-12),
            intensity_unit=scan.pattern.intensity_unit,
            radiation=radiation,
            synthetic=scan.pattern.synthetic,
        )
        background = _linear_background(axis, counts)
        start = float(axis[int(np.argmax(counts - background))])
        half_span = 0.5 * float(axis[-1] - axis[0])
        table = fit_peaks(
            corrected,
            [start],
            radiation=radiation,
            expected_fwhm_deg=expected_fwhm_deg,
            model_doublet=kalpha_doublet_parameters(radiation) is not None,
            window_fwhm=max(1.0, half_span / expected_fwhm_deg),
        )
        if len(table) == 0:
            raise ValueError(
                f"The profile fit at phi = {scan.phi_deg:g}, psi = {scan.psi_deg:g} returned "
                "no peak."
            )
        fit = min(table, key=lambda item: abs(item.two_theta_deg - start))
        return StressPeak(
            phi_deg=scan.phi_deg,
            psi_deg=scan.psi_deg,
            two_theta_deg=float(fit.two_theta_deg),
            two_theta_uncertainty_deg=float(fit.two_theta_standard_uncertainty_deg),
            method="pseudo_voigt",
            fwhm_deg=float(fit.fwhm_deg),
            height=float(fit.height),
            reduced_chi_squared=float(fit.reduced_chi_squared),
            point_count=int(fit.point_count),
            converged=bool(fit.converged),
            lpa_corrected=lpa_correction,
            peak_fit=fit,
        )

    stripped = bool(strip_doublet and kalpha_doublet_parameters(radiation) is not None)
    if stripped:
        # Rachinger's recursion, with the variance carried through it. At the
        # back-reflection angles stress work uses, the doublet is blended into
        # a profile whose width grows with tilt, so a parabola or a centroid of
        # the blend is biased differently at every tilt and the bias does not
        # cancel in the slope: on a Cr K-alpha ferrite (211) measurement it
        # moves the stress by more than 100 MPa.
        single = strip_kalpha2(
            MeasuredPowderPattern(
                name="stress window",
                two_theta_deg=axis,
                intensity=np.maximum(counts, 0.0),
                standard_uncertainty=np.maximum(uncertainty, 1e-12),
                intensity_unit=scan.pattern.intensity_unit,
                radiation=radiation,
            ),
            radiation=radiation,
        )
        counts = np.asarray(single.intensity, dtype=float)
        if single.standard_uncertainty is not None:
            uncertainty = np.asarray(single.standard_uncertainty, dtype=float)
    net = counts - _linear_background(axis, counts)
    # Where the peak is and where its edges fall are decided on a lightly
    # smoothed copy; the location itself uses the measured points. Decided on
    # the raw points, a single noisy count at the threshold moves an edge by a
    # whole step, and the centroid with it, by far more than its propagated
    # uncertainty admits.
    step = float(np.median(np.diff(axis)))
    span = max(3, round(0.25 * expected_fwhm_deg / step) | 1)
    smooth = np.convolve(net, np.ones(span) / span, mode="same")
    peak_index = int(np.argmax(smooth))
    maximum = float(smooth[peak_index])
    if maximum <= 0.0:
        raise ValueError(
            f"No net intensity above background at phi = {scan.phi_deg:g}, psi = {scan.psi_deg:g}."
        )
    if method == "parabola":
        # The contiguous run of points above the threshold that contains the
        # maximum: a stray noise spike elsewhere must not join the fit.
        above = smooth >= top_fraction * maximum
        lo = peak_index
        while lo > 0 and above[lo - 1]:
            lo -= 1
        hi = peak_index
        while hi < axis.size - 1 and above[hi + 1]:
            hi += 1
        if hi - lo + 1 < 5:
            lo, hi = max(0, peak_index - 2), min(axis.size - 1, peak_index + 2)
        x = axis[lo : hi + 1]
        y = net[lo : hi + 1]
        weight = 1.0 / np.maximum(uncertainty[lo : hi + 1], 1e-12) ** 2
        centre = float(x.mean())
        design = np.stack((np.ones_like(x), x - centre, (x - centre) ** 2), axis=1)
        normal = design.T @ (design * weight[:, None])
        covariance = np.linalg.inv(normal)
        coefficients = covariance @ (design.T @ (weight * y))
        residual = y - design @ coefficients
        dof = max(1, x.size - 3)
        chi2 = float(np.sum(weight * residual**2) / dof)
        covariance = covariance * max(1.0, chi2)
        _, b, c = coefficients
        if c >= 0.0:
            raise ValueError(
                f"The parabola at phi = {scan.phi_deg:g}, psi = {scan.psi_deg:g} opens upward; "
                "the peak top is not resolved. Lower top_fraction or use the profile fit."
            )
        vertex = -b / (2.0 * c)
        gradient = np.array([0.0, -1.0 / (2.0 * c), b / (2.0 * c * c)])
        u_vertex = float(np.sqrt(gradient @ covariance @ gradient))
        return StressPeak(
            phi_deg=scan.phi_deg,
            psi_deg=scan.psi_deg,
            two_theta_deg=float(centre + vertex),
            two_theta_uncertainty_deg=max(u_vertex, 1e-9),
            method="parabola",
            height=float(coefficients[0] - b * b / (4.0 * c)),
            reduced_chi_squared=chi2,
            point_count=int(x.size),
            converged=bool(x[0] <= centre + vertex <= x[-1]),
            lpa_corrected=lpa_correction,
            kalpha2_stripped=stripped,
        )

    keep = smooth >= 0.5 * maximum
    lo = peak_index
    while lo > 0 and keep[lo - 1]:
        lo -= 1
    hi = peak_index
    while hi < axis.size - 1 and keep[hi + 1]:
        hi += 1
    # The half-maximum crossings, interpolated between the bracketing points
    # so the window width is continuous in the data too.
    level = 0.5 * maximum

    def crossing(inner: int, outer: int) -> float:
        if outer < 0 or outer >= axis.size:
            return float(axis[inner])
        y_in, y_out = float(smooth[inner]), float(smooth[outer])
        if y_in == y_out:
            return float(axis[inner])
        fraction = (y_in - level) / (y_in - y_out)
        return float(axis[inner] + fraction * (axis[outer] - axis[inner]))

    half_width = 0.5 * (crossing(hi, hi + 1) - crossing(lo, lo - 1))
    centre = float(axis[peak_index])
    weights = np.zeros(axis.size)
    for _ in range(8):
        weights, point_count = _window_centroid_weights(axis, centre, half_width)
        numerator = weights[0] @ net
        denominator = weights[1] @ net
        if denominator <= 0.0:
            raise ValueError(
                f"No net intensity inside the centroid window at phi = {scan.phi_deg:g}, "
                f"psi = {scan.psi_deg:g}."
            )
        updated = float(numerator / denominator)
        moved = abs(updated - centre)
        centre = updated
        if moved < 1e-9:
            break
    # At a fixed window the centroid is a ratio of two linear functionals of
    # the intensities, F = (a . I)/(b . I), with dF/dI_i = (a_i - F b_i)/(b . I).
    # But the window is centred on the answer, so the centroid is the fixed
    # point c = F(I, c), and dc/dI = (dF/dI) / (1 - dF/dc). With the edges at
    # half maximum dF/dc is about 0.6: ignoring it understates u by a factor
    # of about 2.5, which is exactly what a chi-squared of 6 then reports.
    gradient = (weights[0] - centre * weights[1]) / float(weights[1] @ net)
    shift = 1e-4 * half_width
    ahead, _ = _window_centroid_weights(axis, centre + shift, half_width)
    behind, _ = _window_centroid_weights(axis, centre - shift, half_width)
    feedback = (
        float(ahead[0] @ net) / float(ahead[1] @ net)
        - float(behind[0] @ net) / float(behind[1] @ net)
    ) / (2.0 * shift)
    if feedback < 0.95:
        gradient = gradient / (1.0 - feedback)
    variance = np.maximum(uncertainty, 1e-12) ** 2
    u_centroid = float(np.sqrt(np.sum(gradient**2 * variance)))
    return StressPeak(
        phi_deg=scan.phi_deg,
        psi_deg=scan.psi_deg,
        two_theta_deg=centre,
        two_theta_uncertainty_deg=max(u_centroid, 1e-9),
        method="centroid",
        fwhm_deg=2.0 * half_width,
        height=maximum,
        point_count=int(point_count),
        converged=bool(point_count >= 3 and moved < 1e-6),
        lpa_corrected=lpa_correction,
        kalpha2_stripped=stripped,
    )


def _linear_basis_integrals(
    start: float, stop: float, anchor: float, sign: float, width: float
) -> tuple[float, float]:
    """Integrals of ``sign (x - anchor)/width`` and of ``x`` times it over ``[start, stop]``."""

    zeroth = sign * ((stop - anchor) ** 2 - (start - anchor) ** 2) / (2.0 * width)
    moment = sign * ((stop**3 - start**3) / 3.0 - anchor * (stop**2 - start**2) / 2.0) / width
    return zeroth, moment


def _window_centroid_weights(
    axis: np.ndarray, centre: float, half_width: float
) -> tuple[np.ndarray, int]:
    """Return the linear functionals ``(a, b)`` of a continuous windowed centroid.

    The profile is taken as the piecewise-linear interpolant of the measured
    points and integrated exactly over ``[centre - half_width, centre +
    half_width]``: ``a . I`` is the first moment and ``b . I`` the area. Both
    are continuous in ``centre``, so the centroid does not jump when a point
    crosses the window edge -- which a sum over the points inside the window
    does, by a whole step times the edge intensity.
    """

    lower = max(centre - half_width, float(axis[0]))
    upper = min(centre + half_width, float(axis[-1]))
    first = np.zeros(axis.size)
    area = np.zeros(axis.size)
    inside = 0
    for index in range(axis.size - 1):
        x0, x1 = float(axis[index]), float(axis[index + 1])
        left, right = max(x0, lower), min(x1, upper)
        if right <= left:
            continue
        inside += 1
        width = x1 - x0
        # On [x0, x1] the interpolant is I0 (x1 - x)/w + I1 (x - x0)/w; each
        # basis function, and x times it, is integrated over [left, right].
        area0, first0 = _linear_basis_integrals(left, right, x1, -1.0, width)
        area1, first1 = _linear_basis_integrals(left, right, x0, 1.0, width)
        area[index] += area0
        area[index + 1] += area1
        first[index] += first0
        first[index + 1] += first1
    return np.stack((first, area)), inside + 1


# ---------------------------------------------------------------------------
# Per-azimuth sin^2(psi) regressions
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Sin2PsiRegression:
    """The classical ``d`` against ``sin^2(psi)`` line at one azimuth.

    Purpose
    -------
    The picture the method is named for, with the numbers read from it. The
    slope is ``d0 1/2 S2 (sigma_phi - sigma_33)``; with both signs of ``psi``
    measured, a ``sin(2 psi)`` term measures the shear ``tau_phi`` that splits
    the two branches; and a curvature test says whether the line is a line.

    Attributes
    ----------
    phi_deg : float
        The azimuth.
    psi_deg, sin2psi, d_angstrom, d_uncertainty_angstrom : numpy.ndarray
        The points, in measurement order.
    intercept_angstrom, intercept_uncertainty_angstrom : float
        ``d`` at ``sin^2(psi) = 0``, the spacing along the surface normal.
    slope_angstrom, slope_uncertainty_angstrom : float
        ``d(d)/d(sin^2 psi)``.
    splitting_angstrom, splitting_uncertainty_angstrom : float
        Coefficient of ``sin(2 psi)``; ``nan`` without both signs of ``psi``.
    curvature_angstrom, curvature_t : float
        Coefficient of an added ``sin^4(psi)`` term and its value in units of its
        own uncertainty. ``|t| > 3`` means ``d`` is not linear in
        ``sin^2(psi)`` -- texture, a steep stress gradient, or a grain-size
        effect -- and the slope is an average of something that is not a slope.
        ``nan`` when too few tilts remain to test it.
    sigma_phi_mpa, sigma_phi_uncertainty_mpa : float
        ``slope / (d0 1/2 S2)``: the normal stress along the azimuth, minus
        ``sigma_33`` (zero under plane stress).
    tau_phi_mpa, tau_phi_uncertainty_mpa : float
        ``splitting / (d0 1/2 S2)``, or ``nan``.
    reduced_chi_squared : float
        Of the line (or line plus splitting).
    r_squared : float
        Weighted coefficient of determination of the line.
    degrees_of_freedom : int
    """

    phi_deg: float
    psi_deg: np.ndarray
    sin2psi: np.ndarray
    d_angstrom: np.ndarray
    d_uncertainty_angstrom: np.ndarray
    intercept_angstrom: float
    intercept_uncertainty_angstrom: float
    slope_angstrom: float
    slope_uncertainty_angstrom: float
    splitting_angstrom: float
    splitting_uncertainty_angstrom: float
    curvature_angstrom: float
    curvature_t: float
    sigma_phi_mpa: float
    sigma_phi_uncertainty_mpa: float
    tau_phi_mpa: float
    tau_phi_uncertainty_mpa: float
    reduced_chi_squared: float
    r_squared: float
    degrees_of_freedom: int

    @property
    def has_splitting_term(self) -> bool:
        """Whether both signs of ``psi`` were measured, so shear was fitted."""

        return math.isfinite(self.splitting_angstrom)

    def fitted_d(self, psi_deg: Any) -> np.ndarray:
        """Evaluate the fitted line (with its splitting term) at tilts ``psi_deg``."""

        psi = np.deg2rad(np.asarray(psi_deg, dtype=float))
        value = self.intercept_angstrom + self.slope_angstrom * np.sin(psi) ** 2
        if self.has_splitting_term:
            value = value + self.splitting_angstrom * np.sin(2.0 * psi)
        result: np.ndarray = value
        return result

    def branch_averages(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(sin2psi, a1, u(a1), a2, u(a2))`` over tilts measured at both signs.

        ``a1 = (d+ + d-)/2`` is free of shear and linear in ``sin^2(psi)``;
        ``a2 = (d+ - d-)/2`` is ``d0 1/2 S2 tau_phi sin|2 psi|``. They are the
        two plots of the Doelle-Hauk evaluation. Empty arrays when no tilt was
        measured at both signs.
        """

        magnitudes = np.round(np.abs(self.psi_deg), 6)
        rows: list[tuple[float, float, float, float, float]] = []
        for value in np.unique(magnitudes[magnitudes > 0.0]):
            plus = (magnitudes == value) & (self.psi_deg > 0.0)
            minus = (magnitudes == value) & (self.psi_deg < 0.0)
            if not (np.any(plus) and np.any(minus)):
                continue
            d_plus = float(np.mean(self.d_angstrom[plus]))
            d_minus = float(np.mean(self.d_angstrom[minus]))
            u_plus = float(np.sqrt(np.sum(self.d_uncertainty_angstrom[plus] ** 2))) / int(
                np.count_nonzero(plus)
            )
            u_minus = float(np.sqrt(np.sum(self.d_uncertainty_angstrom[minus] ** 2))) / int(
                np.count_nonzero(minus)
            )
            u = 0.5 * math.hypot(u_plus, u_minus)
            rows.append(
                (
                    float(np.sin(np.deg2rad(value)) ** 2),
                    0.5 * (d_plus + d_minus),
                    u,
                    0.5 * (d_plus - d_minus),
                    u,
                )
            )
        if not rows:
            empty = np.zeros(0)
            return empty, empty, empty, empty, empty
        table = np.asarray(rows)
        return table[:, 0], table[:, 1], table[:, 2], table[:, 3], table[:, 4]

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable form."""

        def number(value: float) -> float | None:
            return None if not math.isfinite(value) else float(value)

        return {
            "phi_deg": self.phi_deg,
            "psi_deg": self.psi_deg.tolist(),
            "sin2psi": self.sin2psi.tolist(),
            "d_angstrom": self.d_angstrom.tolist(),
            "d_uncertainty_angstrom": self.d_uncertainty_angstrom.tolist(),
            "intercept_angstrom": self.intercept_angstrom,
            "intercept_uncertainty_angstrom": self.intercept_uncertainty_angstrom,
            "slope_angstrom": self.slope_angstrom,
            "slope_uncertainty_angstrom": self.slope_uncertainty_angstrom,
            "splitting_angstrom": number(self.splitting_angstrom),
            "splitting_uncertainty_angstrom": number(self.splitting_uncertainty_angstrom),
            "curvature_angstrom": number(self.curvature_angstrom),
            "curvature_t": number(self.curvature_t),
            "sigma_phi_mpa": self.sigma_phi_mpa,
            "sigma_phi_uncertainty_mpa": self.sigma_phi_uncertainty_mpa,
            "tau_phi_mpa": number(self.tau_phi_mpa),
            "tau_phi_uncertainty_mpa": number(self.tau_phi_uncertainty_mpa),
            "reduced_chi_squared": number(self.reduced_chi_squared),
            "r_squared": number(self.r_squared),
            "degrees_of_freedom": self.degrees_of_freedom,
        }

    def describe(self) -> str:
        """Return the regression as prose."""

        text = (
            f"At phi = {self.phi_deg:g} degrees, d = {self.intercept_angstrom:.6f} + "
            f"{self.slope_angstrom:.3e} sin^2(psi) angstrom over {self.psi_deg.size} tilts, "
            f"giving sigma_phi = {self.sigma_phi_mpa:.1f} +/- "
            f"{self.sigma_phi_uncertainty_mpa:.1f} MPa (reduced chi-squared "
            f"{self.reduced_chi_squared:.2f})."
        )
        if self.has_splitting_term:
            text += (
                f" Both signs of psi were measured; the psi-splitting term gives "
                f"tau_phi = {self.tau_phi_mpa:.1f} +/- {self.tau_phi_uncertainty_mpa:.1f} MPa."
            )
        if math.isfinite(self.curvature_t) and abs(self.curvature_t) > 3.0:
            text += (
                f" The line is curved (curvature at {self.curvature_t:.1f} of its standard "
                "uncertainty): texture or a stress gradient, which this evaluation does not "
                "model."
            )
        return text


def _weighted_least_squares(
    design: np.ndarray, values: np.ndarray, sigma: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Return ``(coefficients, covariance, chi2, dof)``; covariance unscaled."""

    weight = 1.0 / np.asarray(sigma, dtype=float) ** 2
    normal = design.T @ (design * weight[:, None])
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (weight * values))
    residual = values - design @ coefficients
    chi2 = float(np.sum(weight * residual**2))
    return coefficients, covariance, chi2, int(values.size - design.shape[1])


def _azimuth_groups(phi_deg: np.ndarray) -> list[tuple[float, np.ndarray]]:
    wrapped = np.round(np.mod(phi_deg, 360.0), 6)
    wrapped[wrapped == 360.0] = 0.0
    return [(float(value), np.flatnonzero(wrapped == value)) for value in np.unique(wrapped)]


def fit_sin2psi_lines(
    phi_deg: Any,
    psi_deg: Any,
    d_angstrom: Any,
    d_uncertainty_angstrom: Any,
    *,
    d0_angstrom: float,
    dec: DiffractionElasticConstants,
) -> tuple[Sin2PsiRegression, ...]:
    """Fit ``d`` against ``sin^2(psi)`` at every azimuth separately.

    Purpose
    -------
    The classical evaluation, kept beside the global tensor fit because it is
    the one a reader can check by eye: at each azimuth the points should lie on
    a straight line, and the two branches ``psi > 0`` and ``psi < 0`` should
    coincide unless there is shear.

    Method
    ------
    Weighted least squares of ``d = c0 + c1 sin^2(psi) [+ c2 sin(2 psi)]`` with
    weights ``1/u(d)^2``; the ``sin(2 psi)`` column only when both signs of
    ``psi`` are present. Uncertainties are scaled by ``sqrt(chi2_nu)`` when
    that exceeds one (the Birge ratio), so scatter beyond counting statistics
    widens them. Then ``sigma_phi - sigma_33 = c1 / (d0 1/2 S2)`` and
    ``tau_phi = c2 / (d0 1/2 S2)``. The curvature test refits with an added
    ``sin^4(psi)`` column.

    Parameters
    ----------
    phi_deg, psi_deg, d_angstrom, d_uncertainty_angstrom
        One entry per measurement.
    d0_angstrom
        Stress-free spacing; it only scales the slope here.
    dec
        Diffraction elastic constants of the reflection.

    Returns
    -------
    tuple of Sin2PsiRegression
        One per distinct azimuth (modulo 360 degrees), in increasing ``phi``.
        Azimuths with fewer than two distinct tilts are skipped.
    """

    phi = np.asarray(phi_deg, dtype=float)
    psi = np.asarray(psi_deg, dtype=float)
    d = np.asarray(d_angstrom, dtype=float)
    u = np.asarray(d_uncertainty_angstrom, dtype=float)
    if not (phi.shape == psi.shape == d.shape == u.shape):
        raise ValueError("phi, psi, d and u(d) must have one entry per measurement.")
    scale = d0_angstrom * dec.half_s2_per_tpa * _PER_TPA_TO_PER_MPA
    regressions: list[Sin2PsiRegression] = []
    for azimuth, members in _azimuth_groups(phi):
        tilts = psi[members]
        if len(np.unique(np.round(np.sin(np.deg2rad(tilts)) ** 2, 9))) < 2:
            continue
        sin2 = np.sin(np.deg2rad(tilts)) ** 2
        columns = [np.ones_like(sin2), sin2]
        split = bool(np.any(tilts > 1e-9) and np.any(tilts < -1e-9) and tilts.size >= 4)
        if split:
            columns.append(np.sin(2.0 * np.deg2rad(tilts)))
        design = np.stack(columns, axis=1)
        coefficients, covariance, chi2, dof = _weighted_least_squares(
            design, d[members], u[members]
        )
        reduced = chi2 / dof if dof > 0 else float("nan")
        birge = max(1.0, reduced) if math.isfinite(reduced) else 1.0
        scaled = covariance * birge
        weight = 1.0 / u[members] ** 2
        mean = float(np.sum(weight * d[members]) / np.sum(weight))
        total = float(np.sum(weight * (d[members] - mean) ** 2))
        r_squared = 1.0 - chi2 / total if total > 0.0 else float("nan")
        curvature, curvature_t = float("nan"), float("nan")
        if tilts.size - design.shape[1] - 1 >= 1 and len(np.unique(np.round(sin2, 9))) >= 3:
            extended = np.column_stack((design, sin2**2))
            ext_coefficients, ext_covariance, ext_chi2, ext_dof = _weighted_least_squares(
                extended, d[members], u[members]
            )
            ext_reduced = ext_chi2 / ext_dof if ext_dof > 0 else 1.0
            ext_u = math.sqrt(ext_covariance[-1, -1] * max(1.0, ext_reduced))
            curvature = float(ext_coefficients[-1])
            curvature_t = curvature / ext_u if ext_u > 0.0 else float("nan")
        splitting = float(coefficients[2]) if split else float("nan")
        u_splitting = math.sqrt(scaled[2, 2]) if split else float("nan")
        regressions.append(
            Sin2PsiRegression(
                phi_deg=azimuth,
                psi_deg=tilts.copy(),
                sin2psi=sin2,
                d_angstrom=d[members].copy(),
                d_uncertainty_angstrom=u[members].copy(),
                intercept_angstrom=float(coefficients[0]),
                intercept_uncertainty_angstrom=math.sqrt(scaled[0, 0]),
                slope_angstrom=float(coefficients[1]),
                slope_uncertainty_angstrom=math.sqrt(scaled[1, 1]),
                splitting_angstrom=splitting,
                splitting_uncertainty_angstrom=u_splitting,
                curvature_angstrom=curvature,
                curvature_t=curvature_t,
                sigma_phi_mpa=float(coefficients[1]) / scale,
                sigma_phi_uncertainty_mpa=math.sqrt(scaled[1, 1]) / scale,
                tau_phi_mpa=splitting / scale if split else float("nan"),
                tau_phi_uncertainty_mpa=u_splitting / scale if split else float("nan"),
                reduced_chi_squared=reduced,
                r_squared=r_squared,
                degrees_of_freedom=dof,
            )
        )
    return tuple(regressions)


# ---------------------------------------------------------------------------
# The global tensor fit
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StressTensorFit:
    """The stress tensor from all tilts and azimuths at once, with its budget.

    Attributes
    ----------
    stress_state : str
        One of :data:`STRESS_STATES`.
    component_names : tuple of str
        The free components, in the order of every array below.
    values_mpa : numpy.ndarray
        The fitted components.
    covariance_internal_mpa2 : numpy.ndarray
        ``(A^T W A)^-1`` from the propagated peak-position uncertainties alone.
    reduced_chi_squared : float
        Of the strain fit. Near one when the model and the peak uncertainties
        are right.
    degrees_of_freedom : int
    budget_mpa : mapping of str to numpy.ndarray
        Standard uncertainty of each component from each source:
        ``"statistical"`` (propagated, Birge-scaled when chi-squared exceeds
        one), ``"d0"`` and ``"elastic constants"``.
    combined_covariance_mpa2 : numpy.ndarray
        The covariance of all sources together; its diagonal is the combined
        standard uncertainty squared.
    fitted_strain, residual_strain : numpy.ndarray
        Per measurement, in the order of the result's points.
    monte_carlo_uncertainty_mpa : numpy.ndarray, optional
        Standard deviation of the components over Monte Carlo draws of every
        input from its uncertainty, a check on the linear propagation.
    monte_carlo_draws : int
    """

    stress_state: StressState
    component_names: tuple[str, ...]
    values_mpa: np.ndarray
    covariance_internal_mpa2: np.ndarray
    reduced_chi_squared: float
    degrees_of_freedom: int
    budget_mpa: Mapping[str, np.ndarray]
    combined_covariance_mpa2: np.ndarray
    fitted_strain: np.ndarray
    residual_strain: np.ndarray
    monte_carlo_uncertainty_mpa: np.ndarray | None = None
    monte_carlo_draws: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "budget_mpa", MappingProxyType(dict(self.budget_mpa)))

    @property
    def combined_uncertainty_mpa(self) -> np.ndarray:
        """Combined standard uncertainty of each component."""

        result: np.ndarray = np.sqrt(np.diag(self.combined_covariance_mpa2))
        return result

    @property
    def correlation(self) -> np.ndarray:
        """Correlation matrix of the components from the combined covariance."""

        u = self.combined_uncertainty_mpa
        with np.errstate(divide="ignore", invalid="ignore"):
            result: np.ndarray = self.combined_covariance_mpa2 / np.outer(u, u)
        return result

    def component(self, name: str) -> tuple[float, float]:
        """Return ``(value, combined u)`` of one component; ``(0, 0)`` if not free."""

        if name not in self.component_names:
            if name not in STRESS_COMPONENTS:
                raise KeyError(f"Unknown stress component {name!r}.")
            return 0.0, 0.0
        index = self.component_names.index(name)
        return float(self.values_mpa[index]), float(self.combined_uncertainty_mpa[index])

    @property
    def tensor_mpa(self) -> np.ndarray:
        """The full symmetric 3x3 stress tensor; components not fitted are zero."""

        tensor = np.zeros((3, 3))
        for name, value in zip(self.component_names, self.values_mpa, strict=True):
            i, j = _COMPONENT_INDEX[name]
            tensor[i, j] = tensor[j, i] = value
        return tensor

    def _derived(self, function: Any) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate a function of the components and propagate the covariance."""

        values = np.asarray(function(self.values_mpa), dtype=float)
        steps = np.maximum(1e-6 * np.abs(self.values_mpa), 1e-6)
        jacobian = np.empty((values.size, self.values_mpa.size))
        for column, step in enumerate(steps):
            shift = np.zeros_like(self.values_mpa)
            shift[column] = step
            upper = np.asarray(function(self.values_mpa + shift), dtype=float)
            lower = np.asarray(function(self.values_mpa - shift), dtype=float)
            jacobian[:, column] = (upper - lower).ravel() / (2.0 * step)
        covariance = jacobian @ self.combined_covariance_mpa2 @ jacobian.T
        return values.ravel(), np.sqrt(np.maximum(np.diag(covariance), 0.0))

    def _in_plane(self, values: np.ndarray) -> np.ndarray:
        names = self.component_names
        s11 = values[names.index("sigma_11")]
        s22 = values[names.index("sigma_22")]
        s12 = values[names.index("sigma_12")]
        mean = 0.5 * (s11 + s22)
        radius = math.hypot(0.5 * (s11 - s22), s12)
        angle = 0.5 * math.degrees(math.atan2(2.0 * s12, s11 - s22))
        return np.array([mean + radius, mean - radius, angle])

    def in_plane_principal(self) -> dict[str, float]:
        """Return the in-plane principal stresses, their direction, and uncertainties.

        ``sigma_I >= sigma_II`` are the eigenvalues of the in-plane 2x2 block and
        ``angle_deg`` is the azimuth of ``sigma_I`` from ``S1``,
        ``1/2 atan2(2 sigma_12, sigma_11 - sigma_22)``. The angle is undefined
        when the two principal stresses are equal, and its uncertainty then
        diverges -- a correct statement, not a numerical accident.
        """

        values, uncertainties = self._derived(self._in_plane)
        return {
            "sigma_I_mpa": float(values[0]),
            "sigma_II_mpa": float(values[1]),
            "angle_deg": float(values[2]),
            "sigma_I_uncertainty_mpa": float(uncertainties[0]),
            "sigma_II_uncertainty_mpa": float(uncertainties[1]),
            "angle_uncertainty_deg": float(uncertainties[2]),
        }

    def _von_mises(self, values: np.ndarray) -> np.ndarray:
        tensor = np.zeros((3, 3))
        for name, value in zip(self.component_names, values, strict=True):
            i, j = _COMPONENT_INDEX[name]
            tensor[i, j] = tensor[j, i] = value
        deviator = tensor - np.trace(tensor) / 3.0 * np.eye(3)
        return np.array([math.sqrt(1.5 * float(np.sum(deviator * deviator)))])

    def von_mises(self) -> tuple[float, float]:
        """Return the von Mises equivalent stress and its standard uncertainty.

        The uncertainty is the linear propagation, which is poor near a stress
        state of zero equivalent stress (the function has a cusp there).
        """

        values, uncertainties = self._derived(self._von_mises)
        return float(values[0]), float(uncertainties[0])

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable form."""

        principal = self.in_plane_principal()
        equivalent, equivalent_u = self.von_mises()
        return {
            "stress_state": self.stress_state,
            "component_names": list(self.component_names),
            "values_mpa": self.values_mpa.tolist(),
            "combined_uncertainty_mpa": self.combined_uncertainty_mpa.tolist(),
            "budget_mpa": {key: value.tolist() for key, value in self.budget_mpa.items()},
            "combined_covariance_mpa2": self.combined_covariance_mpa2.tolist(),
            "covariance_internal_mpa2": self.covariance_internal_mpa2.tolist(),
            "reduced_chi_squared": self.reduced_chi_squared,
            "degrees_of_freedom": self.degrees_of_freedom,
            "tensor_mpa": self.tensor_mpa.tolist(),
            "in_plane_principal": principal,
            "von_mises_mpa": equivalent,
            "von_mises_uncertainty_mpa": equivalent_u,
            "monte_carlo_uncertainty_mpa": (
                None
                if self.monte_carlo_uncertainty_mpa is None
                else self.monte_carlo_uncertainty_mpa.tolist()
            ),
            "monte_carlo_draws": self.monte_carlo_draws,
        }


@dataclass(frozen=True, slots=True)
class ResidualStressResult:
    """A complete sin^2(psi) evaluation: the peaks, the lines, the tensor.

    Attributes
    ----------
    wavelength_angstrom : float
        K-alpha1 wavelength every spacing is referred to.
    d0_angstrom, d0_uncertainty_angstrom : float
        The stress-free spacing used, and its standard uncertainty.
    d0_refined : bool
        Whether ``d0`` was determined from the data under plane stress rather
        than given.
    dec : DiffractionElasticConstants
    peaks : tuple of StressPeak
        One per measurement, in measurement order.
    d_angstrom, d_uncertainty_angstrom, strain, strain_uncertainty : numpy.ndarray
        Per measurement.
    regressions : tuple of Sin2PsiRegression
        One per azimuth.
    tensor : StressTensorFit, optional
        ``None`` when the measurement does not determine the requested tensor;
        ``tensor_unavailable_reason`` then says why.
    tensor_unavailable_reason : str, optional
    reflection_label : str
        E.g. ``"(211)"``.
    phase_name : str, optional
    """

    wavelength_angstrom: float
    d0_angstrom: float
    d0_uncertainty_angstrom: float
    d0_refined: bool
    dec: DiffractionElasticConstants
    peaks: tuple[StressPeak, ...]
    d_angstrom: np.ndarray
    d_uncertainty_angstrom: np.ndarray
    strain: np.ndarray
    strain_uncertainty: np.ndarray
    regressions: tuple[Sin2PsiRegression, ...]
    tensor: StressTensorFit | None
    tensor_unavailable_reason: str | None = None
    reflection_label: str = ""
    phase_name: str | None = None

    @property
    def phi_deg(self) -> np.ndarray:
        """Azimuth of each measurement."""

        return np.array([peak.phi_deg for peak in self.peaks])

    @property
    def psi_deg(self) -> np.ndarray:
        """Tilt of each measurement."""

        return np.array([peak.psi_deg for peak in self.peaks])

    @property
    def sin2psi(self) -> np.ndarray:
        """``sin^2(psi)`` of each measurement."""

        result: np.ndarray = np.sin(np.deg2rad(self.psi_deg)) ** 2
        return result

    @property
    def two_theta_deg(self) -> np.ndarray:
        """Located peak position of each measurement."""

        return np.array([peak.two_theta_deg for peak in self.peaks])

    def to_json(self) -> dict[str, Any]:
        """Return the JSON-serializable form, schema-tagged."""

        return {
            "schema": RESIDUAL_STRESS_SCHEMA,
            "wavelength_angstrom": self.wavelength_angstrom,
            "d0_angstrom": self.d0_angstrom,
            "d0_uncertainty_angstrom": self.d0_uncertainty_angstrom,
            "d0_refined": self.d0_refined,
            "dec": self.dec.to_json(),
            "reflection": self.reflection_label,
            "phase_name": self.phase_name,
            "peaks": [peak.to_json() for peak in self.peaks],
            "d_angstrom": self.d_angstrom.tolist(),
            "d_uncertainty_angstrom": self.d_uncertainty_angstrom.tolist(),
            "strain": self.strain.tolist(),
            "strain_uncertainty": self.strain_uncertainty.tolist(),
            "regressions": [item.to_json() for item in self.regressions],
            "tensor": None if self.tensor is None else self.tensor.to_json(),
            "tensor_unavailable_reason": self.tensor_unavailable_reason,
        }

    def describe(self) -> str:
        """Return the evaluation as convention-explicit, citation-backed prose."""

        reflection = f" {self.reflection_label}" if self.reflection_label else ""
        of = f" of {self.phase_name}" if self.phase_name else ""
        azimuths = len(self.regressions)
        lines = [
            f"Residual stress{of} by the sin^2(psi) method from the{reflection} reflection: "
            f"{len(self.peaks)} peak positions at {azimuths} azimuth(s), wavelength "
            f"{self.wavelength_angstrom:.6f} angstrom.",
            (
                f"The stress-free spacing d0 = {self.d0_angstrom:.6f} +/- "
                f"{self.d0_uncertainty_angstrom:.6f} angstrom was "
                + (
                    "determined from the data under the plane-stress assumption sigma_33 = 0."
                    if self.d0_refined
                    else "given."
                )
            ),
            self.dec.describe(),
            "Convention: phi is measured in the surface from S1 towards S2 and psi from the "
            "surface normal S3; strain is (d - d0)/d0; tensile stress is positive; "
            "uncertainties are one standard uncertainty u.",
        ]
        if self.tensor is None:
            lines.append(
                "The stress tensor was not determined: "
                + (self.tensor_unavailable_reason or "the data do not constrain it.")
            )
        else:
            fit = self.tensor
            parts = [
                f"{name.replace('sigma_', 'sigma')} = {value:.1f} +/- {u:.1f} MPa"
                for name, value, u in zip(
                    fit.component_names,
                    fit.values_mpa,
                    fit.combined_uncertainty_mpa,
                    strict=True,
                )
            ]
            principal = fit.in_plane_principal()
            lines.append(
                f"Stress tensor ({fit.stress_state.replace('_', ' ')}): " + "; ".join(parts) + "."
            )
            lines.append(
                f"In-plane principal stresses sigma_I = {principal['sigma_I_mpa']:.1f} and "
                f"sigma_II = {principal['sigma_II_mpa']:.1f} MPa, sigma_I at "
                f"{principal['angle_deg']:.1f} degrees from S1. Strain-fit reduced chi-squared "
                f"{fit.reduced_chi_squared:.2f} on {fit.degrees_of_freedom} degrees of freedom."
            )
        lines.extend(item.describe() for item in self.regressions)
        lines.append(
            "This is the macroscopic (type I) stress of the irradiated layer, averaged over the "
            "penetration depth, for an untextured material. Curvature of d against "
            "sin^2(psi) would indicate texture or a stress gradient, which this evaluation does "
            "not model. "
            + " ".join(
                (_CITATION_MACHERAUCH, _CITATION_NOYAN_COHEN, _CITATION_DOELLE, _CITATION_GUM)
            )
        )
        return "\n".join(lines)


def _identifiable(design: np.ndarray) -> bool:
    if design.shape[0] < design.shape[1]:
        return False
    singular = np.linalg.svd(design, compute_uv=False)
    return bool(singular[-1] > 1e-9 * singular[0])


def _unidentifiable_reason(
    state: StressState, phi: np.ndarray, psi: np.ndarray, refine_d0: bool
) -> str:
    azimuths = {round(float(value) % 180.0, 6) for value in phi}
    if len(azimuths) < 3:
        return (
            f"the in-plane tensor needs at least three azimuths not related by 180 degrees; "
            f"{len(azimuths)} were measured. The per-azimuth stresses sigma_phi below are "
            "still valid."
        )
    if state in ("biaxial_shear", "triaxial") and not (np.any(psi > 0) and np.any(psi < 0)):
        return "the shear components sigma_13 and sigma_23 need tilts of both signs."
    if state == "triaxial" and refine_d0:
        return "a triaxial evaluation needs an exact d0; it cannot be refined with sigma_33."
    return "the measured orientations do not constrain every requested component."


def determine_residual_stress(
    peaks: Sequence[StressPeak],
    *,
    wavelength_angstrom: float,
    d0_angstrom: float,
    dec: DiffractionElasticConstants,
    stress_state: StressState = "biaxial",
    d0_uncertainty_angstrom: float = 0.0,
    refine_d0: bool = False,
    monte_carlo_draws: int = 0,
    seed: int = 0,
    reflection_label: str = "",
    phase_name: str | None = None,
) -> ResidualStressResult:
    """Determine the residual stress tensor from located peak positions.

    Purpose
    -------
    Reduce every measured ``(phi, psi, 2 theta)`` to one stress tensor with a
    defensible uncertainty, and keep the per-azimuth lines beside it.

    Method
    ------
    1. ``d = lambda / (2 sin theta)`` and ``u(d) = d cot(theta) u(theta)``.
    2. ``epsilon = (d - d0)/d0`` and ``u(epsilon) = u(d)/d0``.
    3. Weighted linear least squares of ``epsilon = A x`` for the free
       components ``x`` of the stress state, with ``A`` from
       :func:`strain_design_matrix`. With ``refine_d0`` (plane stress only)
       the model is ``d = d0 + A (d0 x)``, linear in ``(d0, d0 x)``.
    4. Uncertainty budget. Statistical: ``(A^T W A)^-1``, scaled by the reduced
       chi-squared when it exceeds one. ``d0``: the exact sensitivity
       ``dx/dd0`` of the linear solution times ``u(d0)``. Elastic constants:
       the sensitivities to ``S1`` and ``1/2 S2`` times their standard
       uncertainties. Sources are independent and combine in quadrature
       (GUM). Optionally, a Monte Carlo propagation of all inputs checks the
       linearization.

    Parameters
    ----------
    peaks
        Located peaks, e.g. from :func:`locate_stress_peak`.
    wavelength_angstrom
        K-alpha1 wavelength.
    d0_angstrom, d0_uncertainty_angstrom
        Stress-free spacing and its standard uncertainty.
    dec
        Diffraction elastic constants of the reflection.
    stress_state
        ``"biaxial"`` (sigma_11, sigma_22, sigma_12), ``"biaxial_shear"``
        (+ sigma_13, sigma_23; needs both signs of psi) or ``"triaxial"`` (all
        six; needs an exact ``d0``).
    refine_d0
        Determine ``d0`` from the data under ``sigma_33 = 0``.
    monte_carlo_draws, seed
        Size and seed of the Monte Carlo cross-check; zero skips it.
    reflection_label, phase_name
        Recorded on the result.

    Returns
    -------
    ResidualStressResult

    Raises
    ------
    ValueError
        For fewer than three peaks, a non-positive ``d0`` or wavelength, or
        ``refine_d0`` with a triaxial state.
    """

    if stress_state not in STRESS_STATES:
        raise ValueError(f"stress_state must be one of {STRESS_STATES}.")
    if len(peaks) < 3:
        raise ValueError("At least three peak positions are needed.")
    if d0_angstrom <= 0.0 or wavelength_angstrom <= 0.0:
        raise ValueError("d0 and the wavelength must be positive.")
    if d0_uncertainty_angstrom < 0.0:
        raise ValueError("u(d0) must be non-negative.")
    if refine_d0 and stress_state == "triaxial":
        raise ValueError(
            "d0 cannot be refined in a triaxial evaluation: sigma_33 and d0 enter the strain "
            "identically at every tilt, so the data cannot separate them."
        )
    peaks = tuple(peaks)
    phi = np.array([peak.phi_deg for peak in peaks])
    psi = np.array([peak.psi_deg for peak in peaks])
    two_theta = np.array([peak.two_theta_deg for peak in peaks])
    u_two_theta = np.array([peak.two_theta_uncertainty_deg for peak in peaks])
    theta = np.deg2rad(0.5 * two_theta)
    d = wavelength_angstrom / (2.0 * np.sin(theta))
    u_d = d / np.tan(theta) * np.deg2rad(0.5 * u_two_theta)

    components = _STATE_COMPONENTS[stress_state]
    design = strain_design_matrix(
        phi,
        psi,
        s1_per_tpa=dec.s1_per_tpa,
        half_s2_per_tpa=dec.half_s2_per_tpa,
        components=components,
    )

    nominal = all(peak.uncertainty_is_nominal for peak in peaks)
    tensor: StressTensorFit | None = None
    reason: str | None = None
    d0 = float(d0_angstrom)
    u_d0 = float(d0_uncertainty_angstrom)

    augmented = np.column_stack((np.ones(d.size), design)) if refine_d0 else design
    if not _identifiable(augmented) or augmented.shape[0] <= augmented.shape[1]:
        reason = _unidentifiable_reason(stress_state, phi, psi, refine_d0)
        if augmented.shape[0] <= augmented.shape[1] and _identifiable(augmented):
            reason = (
                f"{augmented.shape[1]} unknowns need more than {augmented.shape[0]} measurements "
                "to leave a degree of freedom for the uncertainty."
            )
    elif refine_d0:
        coefficients, covariance_p, chi2, dof = _weighted_least_squares(augmented, d, u_d)
        d0 = float(coefficients[0])
        values = coefficients[1:] / d0
        jacobian = np.zeros((values.size, coefficients.size))
        jacobian[:, 0] = -coefficients[1:] / d0**2
        jacobian[:, 1:] = np.eye(values.size) / d0
        covariance = jacobian @ covariance_p @ jacobian.T
        reduced = chi2 / dof
        u_d0 = math.sqrt(covariance_p[0, 0] * max(1.0, reduced))
        tensor = _assemble_fit(
            stress_state,
            components,
            values,
            covariance,
            reduced,
            dof,
            phi,
            psi,
            d,
            u_d,
            d0,
            None,
            dec,
            monte_carlo_draws,
            seed,
            refine_d0=True,
            external_only=nominal,
        )
    else:
        strain = d / d0 - 1.0
        u_strain = u_d / d0
        values, covariance, chi2, dof = _weighted_least_squares(design, strain, u_strain)
        tensor = _assemble_fit(
            stress_state,
            components,
            values,
            covariance,
            chi2 / dof,
            dof,
            phi,
            psi,
            d,
            u_d,
            d0,
            u_d0,
            dec,
            monte_carlo_draws,
            seed,
            refine_d0=False,
            external_only=nominal,
        )

    strain = d / d0 - 1.0
    u_strain = u_d / d0
    regressions = fit_sin2psi_lines(phi, psi, d, u_d, d0_angstrom=d0, dec=dec)
    return ResidualStressResult(
        wavelength_angstrom=float(wavelength_angstrom),
        d0_angstrom=d0,
        d0_uncertainty_angstrom=u_d0,
        d0_refined=bool(refine_d0 and tensor is not None),
        dec=dec,
        peaks=peaks,
        d_angstrom=d,
        d_uncertainty_angstrom=u_d,
        strain=strain,
        strain_uncertainty=u_strain,
        regressions=regressions,
        tensor=tensor,
        tensor_unavailable_reason=reason,
        reflection_label=reflection_label,
        phase_name=phase_name,
    )


def _solve_strain(
    phi: np.ndarray,
    psi: np.ndarray,
    d: np.ndarray,
    u_d: np.ndarray,
    d0: float,
    s1: float,
    half_s2: float,
    components: Sequence[str],
) -> np.ndarray:
    design = strain_design_matrix(
        phi, psi, s1_per_tpa=s1, half_s2_per_tpa=half_s2, components=components
    )
    values, _, _, _ = _weighted_least_squares(design, d / d0 - 1.0, u_d / d0)
    return values


def _solve_refined(
    phi: np.ndarray,
    psi: np.ndarray,
    d: np.ndarray,
    u_d: np.ndarray,
    s1: float,
    half_s2: float,
    components: Sequence[str],
) -> np.ndarray:
    design = strain_design_matrix(
        phi, psi, s1_per_tpa=s1, half_s2_per_tpa=half_s2, components=components
    )
    augmented = np.column_stack((np.ones(d.size), design))
    coefficients, _, _, _ = _weighted_least_squares(augmented, d, u_d)
    result: np.ndarray = coefficients[1:] / coefficients[0]
    return result


def _assemble_fit(
    stress_state: StressState,
    components: tuple[str, ...],
    values: np.ndarray,
    covariance: np.ndarray,
    reduced: float,
    dof: int,
    phi: np.ndarray,
    psi: np.ndarray,
    d: np.ndarray,
    u_d: np.ndarray,
    d0: float,
    u_d0: float | None,
    dec: DiffractionElasticConstants,
    monte_carlo_draws: int,
    seed: int,
    *,
    refine_d0: bool,
    external_only: bool = False,
) -> StressTensorFit:
    """Build the budget, the combined covariance and the Monte Carlo check."""

    # Birge scaling: scatter beyond the propagated uncertainties widens them,
    # and scatter below them is not allowed to narrow them. When the positions
    # carried only placeholder uncertainties the scatter is all there is, and
    # the scaling applies in both directions.
    birge = reduced if external_only else max(1.0, reduced)
    statistical = covariance * birge
    budget: dict[str, np.ndarray] = {"statistical": np.sqrt(np.diag(statistical))}
    combined = statistical.copy()

    def solve(d_values: np.ndarray, d0_value: float, s1: float, half_s2: float) -> np.ndarray:
        if refine_d0:
            return _solve_refined(phi, psi, d_values, u_d, s1, half_s2, components)
        return _solve_strain(phi, psi, d_values, u_d, d0_value, s1, half_s2, components)

    if not refine_d0 and u_d0 is not None:
        # The solution is linear in epsilon = d/d0 - 1, so the derivative is
        # exact: dx/dd0 = G depsilon/dd0 with depsilon_k/dd0 = -d_k/d0^2.
        design = strain_design_matrix(
            phi,
            psi,
            s1_per_tpa=dec.s1_per_tpa,
            half_s2_per_tpa=dec.half_s2_per_tpa,
            components=components,
        )
        weight = 1.0 / (u_d / d0) ** 2
        gain = covariance @ (design.T * weight)
        sensitivity = gain @ (-d / d0**2)
        budget["d0"] = np.abs(sensitivity) * u_d0
        combined = combined + np.outer(sensitivity, sensitivity) * u_d0**2
    if dec.relative_standard_uncertainty > 0.0:
        total = np.zeros((values.size, values.size))
        for which in ("s1", "half_s2"):
            base_s1, base_half = dec.s1_per_tpa, dec.half_s2_per_tpa
            u_value = dec.relative_standard_uncertainty * (
                abs(base_s1) if which == "s1" else base_half
            )
            step = 1e-4 * (abs(base_s1) if which == "s1" else base_half)
            if which == "s1":
                upper = solve(d, d0, base_s1 + step, base_half)
                lower = solve(d, d0, base_s1 - step, base_half)
            else:
                upper = solve(d, d0, base_s1, base_half + step)
                lower = solve(d, d0, base_s1, base_half - step)
            sensitivity = (upper - lower) / (2.0 * step)
            total = total + np.outer(sensitivity, sensitivity) * u_value**2
        budget["elastic constants"] = np.sqrt(np.diag(total))
        combined = combined + total

    design_final = strain_design_matrix(
        phi,
        psi,
        s1_per_tpa=dec.s1_per_tpa,
        half_s2_per_tpa=dec.half_s2_per_tpa,
        components=components,
    )
    fitted = design_final @ values
    residual = (d / d0 - 1.0) - fitted

    monte_carlo: np.ndarray | None = None
    draws = int(monte_carlo_draws)
    if draws > 0:
        generator = np.random.default_rng(seed)
        scale = math.sqrt(birge)
        d_draws = d + generator.standard_normal((draws, d.size)) * (u_d * scale)
        d0_draws = (
            d0 + generator.standard_normal(draws) * (u_d0 or 0.0)
            if not refine_d0
            else np.full(draws, d0)
        )
        relative = dec.relative_standard_uncertainty
        s1_draws = dec.s1_per_tpa + generator.standard_normal(draws) * relative * abs(
            dec.s1_per_tpa
        )
        half_draws = dec.half_s2_per_tpa * (1.0 + generator.standard_normal(draws) * relative)
        coefficients = _projection_coefficients(phi, psi)
        columns = [STRESS_COMPONENTS.index(name) for name in components]
        projection = coefficients[:, columns]
        trace = _TRACE_ROW[columns]
        # One batched solve for every draw: A_k = 1/2 S2_k P + S1_k T.
        design_draws = (
            half_draws[:, None, None] * projection[None, :, :]
            + s1_draws[:, None, None] * trace[None, None, :]
        ) * _PER_TPA_TO_PER_MPA
        if refine_d0:
            design_draws = np.concatenate((np.ones((draws, d.size, 1)), design_draws), axis=2)
            target = d_draws
            weight = 1.0 / u_d**2
        else:
            target = d_draws / d0_draws[:, None] - 1.0
            weight = 1.0 / (u_d / d0) ** 2
        normal = np.einsum("kni,n,knj->kij", design_draws, weight, design_draws)
        right = np.einsum("kni,n,kn->ki", design_draws, weight, target)
        solutions = np.linalg.solve(normal, right[..., None])[..., 0]
        if refine_d0:
            solutions = solutions[:, 1:] / solutions[:, :1]
        monte_carlo = np.std(solutions, axis=0, ddof=1)

    return StressTensorFit(
        stress_state=stress_state,
        component_names=components,
        values_mpa=np.asarray(values, dtype=float),
        covariance_internal_mpa2=covariance,
        reduced_chi_squared=float(reduced),
        degrees_of_freedom=int(dof),
        budget_mpa=budget,
        combined_covariance_mpa2=combined,
        fitted_strain=fitted,
        residual_strain=residual,
        monte_carlo_uncertainty_mpa=monte_carlo,
        monte_carlo_draws=draws,
    )


# ---------------------------------------------------------------------------
# Pipeline, synthetic data and table readers
# ---------------------------------------------------------------------------


def residual_stress_pipeline(
    measurement: Sin2PsiMeasurement,
    *,
    d0_angstrom: float,
    dec: DiffractionElasticConstants,
    expected_two_theta_deg: float | None = None,
    window_deg: float | None = None,
    peak_method: PeakLocationMethod = "pseudo_voigt",
    expected_fwhm_deg: float = 1.0,
    lpa_correction: bool = True,
    strip_doublet: bool = True,
    stress_state: StressState = "biaxial",
    d0_uncertainty_angstrom: float = 0.0,
    refine_d0: bool = False,
    monte_carlo_draws: int = 0,
    seed: int = 0,
    reflection_label: str = "",
    phase_name: str | None = None,
) -> ResidualStressResult:
    """Locate every peak of a measurement and determine the stress.

    The one call a routine evaluation makes: :func:`locate_stress_peak` on each
    scan, then :func:`determine_residual_stress`. Both halves stay available
    separately for the evaluation that needs to edit the peak list between
    them. ``expected_two_theta_deg`` defaults to the Bragg angle of ``d0``.
    """

    wavelength = measurement.radiation.wavelength_angstrom
    if expected_two_theta_deg is None:
        ratio = wavelength / (2.0 * d0_angstrom)
        if not 0.0 < ratio < 1.0:
            raise ValueError("The reflection is beyond the Ewald limit at this wavelength.")
        expected_two_theta_deg = float(np.rad2deg(2.0 * np.arcsin(ratio)))
    peaks = tuple(
        locate_stress_peak(
            scan,
            radiation=measurement.radiation,
            method=peak_method,
            expected_two_theta_deg=expected_two_theta_deg,
            window_deg=window_deg,
            expected_fwhm_deg=expected_fwhm_deg,
            lpa_correction=lpa_correction,
            geometry=measurement.geometry,
            strip_doublet=strip_doublet,
        )
        for scan in measurement.scans
    )
    return determine_residual_stress(
        peaks,
        wavelength_angstrom=wavelength,
        d0_angstrom=d0_angstrom,
        dec=dec,
        stress_state=stress_state,
        d0_uncertainty_angstrom=d0_uncertainty_angstrom,
        refine_d0=refine_d0,
        monte_carlo_draws=monte_carlo_draws,
        seed=seed,
        reflection_label=reflection_label,
        phase_name=phase_name,
    )


def _as_stress_tensor(stress: Any) -> np.ndarray:
    if isinstance(stress, Mapping):
        tensor = np.zeros((3, 3))
        for name, value in stress.items():
            if name not in _COMPONENT_INDEX:
                raise ValueError(f"Unknown stress component {name!r}.")
            i, j = _COMPONENT_INDEX[name]
            tensor[i, j] = tensor[j, i] = float(value)
        return tensor
    tensor = np.asarray(stress, dtype=float)
    if tensor.shape != (3, 3) or not np.allclose(tensor, tensor.T):
        raise ValueError("The stress must be a symmetric 3x3 tensor or a component mapping.")
    return tensor


def simulate_sin2psi_measurement(
    *,
    d0_angstrom: float,
    stress_mpa: Any,
    dec: DiffractionElasticConstants,
    phi_deg: Iterable[float] = (0.0, 45.0, 90.0),
    psi_deg: Iterable[float] = (-45.0, -35.3, -24.1, 0.0, 24.1, 35.3, 45.0),
    radiation: RadiationSpec | None = None,
    fwhm_deg: float = 1.2,
    eta: float = 0.5,
    peak_counts: float = 3000.0,
    background_counts: float = 150.0,
    window_half_width_deg: float = 5.0,
    step_deg: float = 0.05,
    geometry: TiltGeometry = "omega",
    apply_lpa: bool = True,
    seed: int = 0,
    name: str = "synthetic sin^2(psi) measurement",
) -> Sin2PsiMeasurement:
    """Generate a measurement of a known stress, with counting noise.

    Purpose
    -------
    Give the method a case whose answer is known before it runs: every
    estimator and every uncertainty can then be checked against the truth,
    which no real specimen allows.

    Method
    ------
    At each ``(phi, psi)`` the strain is ``1/2 S2 m.sigma.m + S1 tr(sigma)``,
    the spacing ``d0 (1 + epsilon)`` and the K-alpha1 angle from Bragg's law.
    The profile is a pseudo-Voigt doublet (the K-alpha2 partner at its
    Bragg-law position and tabulated ratio), broadened as ``1/cos(psi)`` under
    omega tilting (defocusing), multiplied by :func:`lpa_factor` when
    ``apply_lpa`` is set, on a flat background, then given Poisson noise.

    Parameters
    ----------
    d0_angstrom
        Stress-free spacing.
    stress_mpa
        3x3 tensor or mapping such as ``{"sigma_11": -300.0}``.
    dec
        Diffraction elastic constants used to generate the strain.
    phi_deg, psi_deg
        Azimuths and tilts; every combination is measured.
    radiation
        Defaults to Cr K-alpha, the conventional choice for ferritic steel.
    fwhm_deg, eta, peak_counts, background_counts
        Profile width at ``psi = 0``, Lorentzian fraction, peak height and
        background, in counts.
    window_half_width_deg, step_deg
        Scan range about the peak and step.
    geometry, apply_lpa
        Tilt geometry and whether the LPA factor shapes the profile.
    seed
        Seed of the counting noise.
    name
        Recorded on the measurement.

    Returns
    -------
    Sin2PsiMeasurement
        With ``synthetic=True`` and ``true_stress_mpa`` set.
    """

    radiation = radiation or RadiationSpec.cr_ka()
    stress = _as_stress_tensor(stress_mpa)
    generator = np.random.default_rng(seed)
    phis = [float(value) for value in phi_deg]
    psis = [float(value) for value in psi_deg]
    ratio = kalpha_doublet_parameters(radiation)
    scans: list[StressScan] = []
    for phi in phis:
        for psi in psis:
            m = measurement_direction(phi, psi)
            strain = (
                dec.half_s2_per_tpa * float(m @ stress @ m)
                + dec.s1_per_tpa * float(np.trace(stress))
            ) * _PER_TPA_TO_PER_MPA
            d = d0_angstrom * (1.0 + strain)
            centre = float(np.rad2deg(2.0 * np.arcsin(radiation.wavelength_angstrom / (2.0 * d))))
            axis = np.arange(
                centre - window_half_width_deg,
                centre + window_half_width_deg + 0.5 * step_deg,
                step_deg,
            )
            width = fwhm_deg / math.cos(math.radians(psi)) if geometry == "omega" else fwhm_deg
            profile = pseudo_voigt_profile(axis, centre_deg=centre, fwhm_deg=width, eta=eta)
            if ratio is not None:
                partner = float(
                    np.rad2deg(2.0 * np.arcsin(ratio[0] * np.sin(np.deg2rad(0.5 * centre))))
                )
                profile = profile + ratio[1] * pseudo_voigt_profile(
                    axis, centre_deg=partner, fwhm_deg=width, eta=eta
                )
            if apply_lpa:
                factor = lpa_factor(axis, psi_deg=psi, geometry=geometry)
                profile = profile * factor / factor[axis.size // 2]
            expected = peak_counts * profile + background_counts
            counts = generator.poisson(expected).astype(float)
            scans.append(
                StressScan(
                    phi_deg=phi,
                    psi_deg=psi,
                    pattern=MeasuredPowderPattern(
                        name=f"phi {phi:g}, psi {psi:g}",
                        two_theta_deg=axis,
                        intensity=counts,
                        radiation=radiation,
                        synthetic=True,
                    ),
                )
            )
    return Sin2PsiMeasurement(
        name=name,
        scans=tuple(scans),
        radiation=radiation,
        geometry=geometry,
        synthetic=True,
        true_stress_mpa=stress,
        metadata={"d0_angstrom": f"{d0_angstrom:.8f}", "seed": str(seed)},
    )


def _numeric_rows(text: str, *, minimum_columns: int, what: str) -> list[list[float]]:
    rows: list[list[float]] = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        fields = stripped.replace(",", " ").replace(";", " ").split()
        try:
            values = [float(item) for item in fields]
        except ValueError:
            if not rows:
                # A header line before the data is common in exported tables.
                continue
            raise ValueError(f"Line {number} of the {what} is not numeric: {stripped!r}.") from None
        if len(values) < minimum_columns:
            raise ValueError(
                f"Line {number} of the {what} has {len(values)} column(s); at least "
                f"{minimum_columns} are needed."
            )
        rows.append(values)
    if not rows:
        raise ValueError(f"The {what} holds no numeric rows.")
    return rows


def parse_stress_scans(
    text: str,
    *,
    radiation: RadiationSpec,
    geometry: TiltGeometry = "omega",
    name: str = "measured sin^2(psi) scans",
) -> Sin2PsiMeasurement:
    """Read scans from a four-column table ``phi psi 2theta intensity``.

    One row per measured point; the rows of one scan share ``(phi, psi)`` and
    need not be contiguous. Whitespace, commas or semicolons separate columns,
    ``#`` starts a comment, and a header line before the data is skipped. The
    long format is what most diffractometer software can export for a whole
    stress measurement in one file.
    """

    rows = np.asarray(_numeric_rows(text, minimum_columns=4, what="scan table"))
    keys = np.round(rows[:, :2], 6)
    scans: list[StressScan] = []
    unique = np.unique(keys, axis=0)
    for phi, psi in unique:
        members = np.all(keys == (phi, psi), axis=1)
        block = rows[members]
        order = np.argsort(block[:, 2])
        block = block[order]
        if block.shape[0] < 10:
            raise ValueError(
                f"The scan at phi = {phi:g}, psi = {psi:g} has {block.shape[0]} points; at least "
                "10 are needed to locate a peak."
            )
        scans.append(
            StressScan(
                phi_deg=float(phi),
                psi_deg=float(psi),
                pattern=MeasuredPowderPattern(
                    name=f"phi {phi:g}, psi {psi:g}",
                    two_theta_deg=block[:, 2],
                    intensity=np.maximum(block[:, 3], 0.0),
                    radiation=radiation,
                ),
            )
        )
    return Sin2PsiMeasurement(name=name, scans=tuple(scans), radiation=radiation, geometry=geometry)


def parse_stress_peak_positions(text: str) -> tuple[StressPeak, ...]:
    """Read located peaks from a table ``phi psi 2theta [u(2theta)]``.

    For peak positions already located by other software. Without a fourth
    column every position gets the same nominal uncertainty of 0.001 degrees,
    so the fit is unweighted and its uncertainties come entirely from the
    scatter (the Birge scaling): they are then external, not propagated.
    """

    rows = _numeric_rows(text, minimum_columns=3, what="peak-position table")
    peaks = []
    for values in rows:
        uncertainty = values[3] if len(values) >= 4 and values[3] > 0.0 else 1e-3
        peaks.append(
            StressPeak(
                phi_deg=values[0],
                psi_deg=values[1],
                two_theta_deg=values[2],
                two_theta_uncertainty_deg=uncertainty,
                method="given",
                uncertainty_is_nominal=not (len(values) >= 4 and values[3] > 0.0),
            )
        )
    return tuple(peaks)
