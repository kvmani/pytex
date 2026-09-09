"""High-Resolution Electron Microscopy (HREM / HRTEM) simulation physics and models.

This module provides the canonical crystallographic and electron-optics data model
for high-resolution transmission electron microscopy simulations, with support for:

1. Double-corrected microscope optics (spherical and chromatic aberration correction,
   higher-order aberrations up to 5th order, partial coherence damping envelopes,
   and standard imaging regimes such as Scherzer, Lichte, Cs-corrected, and NCSI).
2. Atomic structure snapshots representing crystalline slabs (oriented along any
   zone axis), crystals with defects (vacancies, dislocations), and amorphous foils.
3. Contrast Transfer Function (CTF) evaluation, 2D exit-wave propagation, image
   intensity synthesis, FFT power spectra (Thon rings), and explainable report surfaces.

Scientific References
---------------------
- Scherzer, O. (1949). The theoretical resolution limit of the electron microscope.
  J. Appl. Phys. 20, 20-29.
- Cowley, J. M. & Moodie, A. F. (1957). The scattering of electrons by atoms and
  crystals. I. A new theoretical approach. Acta Crystallogr. 10, 609-619.
- Frank, J. (1973). An envelope for the transfer function of the electron microscope.
  Optik 38, 519-536.
- Haider, M. et al. (1998). Electron microscopy image with a spherical-aberration-corrected
  objective lens. Nature 392, 768-770.
- Kabius, B. et al. (2009). First application of an un-monochromated, Cs and Cc
  corrected transmission electron microscope: High resolution imaging and EELS.
  Microsc. Microanal. 15 (Suppl 2), 1150-1151.
- Urban, K. W. et al. (2009). Negative spherical aberration imaging in transmission
  electron microscopy. Phil. Trans. R. Soc. A 367, 3735-3753.
- Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, 2nd ed.,
  Springer.
- Madsen, J. et al. (2021). abTEM: An open-source framework for simulation of
  transmission electron microscopy. ChemPhysChem 22, 1-13.
"""

from __future__ import annotations

import base64
import io
import math
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from pytex.core._arrays import as_float_array
from pytex.core.lattice import Phase

if TYPE_CHECKING:
    import ase

# Physical constants (SI units and CODATA 2018 values)
_ELECTRON_MASS_KG = 9.1093837015e-31
_ELEMENTARY_CHARGE_C = 1.602176634e-19
_PLANCK_CONSTANT_JS = 6.62607015e-34
_SPEED_OF_LIGHT_M_S = 2.99792458e8


def relativistic_wavelength_angstrom(energy_kev: float) -> float:
    """Compute relativistic electron de Broglie wavelength in Angstrom.

    Parameters
    ----------
    energy_kev : float
        Electron kinetic energy / accelerating voltage in keV.

    Returns
    -------
    float
        Relativistic de Broglie wavelength in Angstrom (Å).

    Raises
    ------
    ValueError
        If ``energy_kev`` is non-positive.

    Examples
    --------
    >>> round(relativistic_wavelength_angstrom(200.0), 6)
    0.025079
    >>> round(relativistic_wavelength_angstrom(300.0), 6)
    0.019687
    >>> round(relativistic_wavelength_angstrom(80.0), 6)
    0.041756
    """
    if energy_kev <= 0.0:
        raise ValueError(f"Accelerating energy must be positive, got {energy_kev} keV.")
    voltage_v = energy_kev * 1000.0
    ev = _ELEMENTARY_CHARGE_C * voltage_v
    m0 = _ELECTRON_MASS_KG
    c = _SPEED_OF_LIGHT_M_S
    # Relativistic momentum p = sqrt(2*m0*E_k*(1 + E_k/(2*m0*c^2)))
    rel_factor = 1.0 + ev / (2.0 * m0 * c * c)
    momentum_kg_m_s = math.sqrt(2.0 * m0 * ev * rel_factor)
    wavelength_m = _PLANCK_CONSTANT_JS / momentum_kg_m_s
    return wavelength_m * 1e10


def relativistic_interaction_parameter_inv_v_angstrom(energy_kev: float) -> float:
    """Compute interaction parameter sigma in 1 / (V * Angstrom).

    Parameters
    ----------
    energy_kev : float
        Accelerating voltage in keV.

    Returns
    -------
    float
        Interaction constant sigma in V^-1 * Å^-1.
    """
    wavelength_angstrom = relativistic_wavelength_angstrom(energy_kev)
    voltage_v = energy_kev * 1000.0
    m0_c2 = _ELECTRON_MASS_KG * (_SPEED_OF_LIGHT_M_S**2)
    ev = _ELEMENTARY_CHARGE_C * voltage_v
    rel_factor = (m0_c2 + ev) / (2.0 * m0_c2 + ev)
    return (2.0 * math.pi / (wavelength_angstrom * voltage_v)) * rel_factor


class DoubleCorrectionMode(StrEnum):
    """Microscope lens correction regime."""

    UNCORRECTED = "uncorrected"
    CS_CORRECTED = "cs_corrected"
    DOUBLE_CORRECTED = "double_corrected"
    NCSI = "ncsi"


@dataclass(frozen=True, slots=True)
class MicroscopeAberrations:
    """Optical parameters and wave aberration coefficients of a TEM objective lens.

    Supports aberrations up to 5th order, partial temporal coherence
    (focal spread / chromatic aberration), partial spatial coherence
    (convergence semi-angle), and objective aperture cutoff.

    Attributes
    ----------
    energy_kev : float
        Accelerating voltage / beam energy in keV (e.g. 200.0, 300.0).
    defocus_angstrom : float
        Defocus Delta f (C10) in Angstrom. Convention: negative is underfocus
        (Scherzer underfocus), positive is overfocus.
    cs_mm : float
        3rd-order spherical aberration coefficient Cs (C30) in millimeters.
    c5_mm : float
        5th-order spherical aberration coefficient C50 in millimeters.
    astigmatism_angstrom : float
        2-fold astigmatism amplitude C12 in Angstrom.
    astigmatism_angle_deg : float
        Azimuthal angle phi12 of 2-fold astigmatism in degrees.
    trefoil_angstrom : float
        3-fold astigmatism (trefoil) amplitude C23 in Angstrom.
    trefoil_angle_deg : float
        Azimuthal angle phi23 of trefoil in degrees.
    coma_angstrom : float
        Axial coma amplitude C21 in Angstrom.
    coma_angle_deg : float
        Azimuthal angle phi21 of axial coma in degrees.
    focal_spread_angstrom : float
        Standard deviation Delta of the Gaussian focal spread (temporal coherence
        envelope) in Angstrom.
    convergence_semiangle_mrad : float
        Illumination convergence semi-angle alpha_s (spatial coherence envelope)
        in milliradians.
    aperture_cutoff_mrad : float | None
        Objective aperture semi-angle cutoff in milliradians. None means no aperture.
    mode : DoubleCorrectionMode
        Nominal correction regime descriptor.
    """

    energy_kev: float = 200.0
    defocus_angstrom: float = 0.0
    cs_mm: float = 0.0
    c5_mm: float = 0.0
    astigmatism_angstrom: float = 0.0
    astigmatism_angle_deg: float = 0.0
    trefoil_angstrom: float = 0.0
    trefoil_angle_deg: float = 0.0
    coma_angstrom: float = 0.0
    coma_angle_deg: float = 0.0
    focal_spread_angstrom: float = 20.0
    convergence_semiangle_mrad: float = 0.2
    aperture_cutoff_mrad: float | None = None
    mode: DoubleCorrectionMode = DoubleCorrectionMode.DOUBLE_CORRECTED

    @property
    def wavelength_angstrom(self) -> float:
        """Relativistic electron wavelength in Angstrom."""
        return relativistic_wavelength_angstrom(self.energy_kev)

    @property
    def interaction_parameter_inv_v_angstrom(self) -> float:
        """Relativistic electron interaction parameter sigma in V^-1 Å^-1."""
        return relativistic_interaction_parameter_inv_v_angstrom(self.energy_kev)

    @property
    def cs_um(self) -> float:
        """Spherical aberration in micrometers."""
        return self.cs_mm * 1e3

    @property
    def cs_angstrom(self) -> float:
        """Spherical aberration in Angstrom."""
        return self.cs_mm * 1e7

    @property
    def c5_angstrom(self) -> float:
        """5th-order spherical aberration in Angstrom."""
        return self.c5_mm * 1e7

    @property
    def has_azimuthal_aberrations(self) -> bool:
        """Whether any non-round (azimuth-dependent) aberration coefficient is nonzero.

        Round aberrations — defocus, Cs, C5 — shift the wave aberration by the same
        amount at every azimuth, so a single radial cut of the transfer function
        describes the lens completely. Two-fold astigmatism, axial coma and trefoil
        do not: chi becomes a function of azimuth, and a radial cut describes one
        direction in the image only. This predicate is what tells a caller whether
        :meth:`evaluate_ctf_1d` is the whole story or one section through it.
        """
        return (
            self.astigmatism_angstrom != 0.0
            or self.coma_angstrom != 0.0
            or self.trefoil_angstrom != 0.0
        )

    def residual_aberration_terms(self) -> tuple[tuple[str, str, float, float], ...]:
        """Nonzero non-round aberrations as ``(key, symbol, amplitude_A, azimuth_deg)``.

        Purpose
        -------
        Report the aberrations that a corrected instrument is actually limited by,
        in the order of the Krivanek coefficient naming used by aberration
        correctors, so that a result object can state which terms shaped it.

        Returns
        -------
        tuple
            One entry per nonzero non-round term: its registry key, its Krivanek
            label, its amplitude in Angstrom and its azimuth in degrees. Empty when
            the lens is round, in which case :meth:`evaluate_ctf_1d` is complete.
        """
        terms: list[tuple[str, str, float, float]] = []
        if self.astigmatism_angstrom != 0.0:
            terms.append(
                ("astigmatism_2fold", "C12", self.astigmatism_angstrom, self.astigmatism_angle_deg)
            )
        if self.coma_angstrom != 0.0:
            terms.append(("axial_coma", "C21", self.coma_angstrom, self.coma_angle_deg))
        if self.trefoil_angstrom != 0.0:
            terms.append(("trefoil", "C23", self.trefoil_angstrom, self.trefoil_angle_deg))
        return tuple(terms)

    @property
    def scherzer_defocus_angstrom(self) -> float:
        """Scherzer defocus Delta f_Sch = -1.2 * sqrt(Cs * lambda) in Angstrom."""
        if self.cs_mm <= 0.0:
            return 0.0
        return -1.2 * math.sqrt(self.cs_angstrom * self.wavelength_angstrom)

    @property
    def scherzer_resolution_angstrom(self) -> float:
        """Scherzer point resolution d_Sch = 0.64 * (Cs * lambda^3)^(1/4) in Angstrom."""
        if self.cs_mm <= 0.0:
            return 0.0
        return 0.64 * math.pow(self.cs_angstrom * (self.wavelength_angstrom**3), 0.25)

    @property
    def lichte_defocus_angstrom(self) -> float:
        """Lichte least-aberration defocus Delta f = -0.75 * sqrt(Cs * lambda) in Angstrom."""
        if self.cs_mm <= 0.0:
            return 0.0
        return -0.75 * math.sqrt(self.cs_angstrom * self.wavelength_angstrom)

    @classmethod
    def conventional_tem(
        cls,
        energy_kev: float = 200.0,
        cs_mm: float = 1.0,
        focal_spread_angstrom: float = 35.0,
        convergence_semiangle_mrad: float = 0.5,
    ) -> MicroscopeAberrations:
        """Preset for uncorrected conventional TEM at Scherzer defocus."""
        wavelength = relativistic_wavelength_angstrom(energy_kev)
        cs_angstrom = cs_mm * 1e7
        scherzer_defocus = -1.2 * math.sqrt(cs_angstrom * wavelength)
        return cls(
            energy_kev=energy_kev,
            defocus_angstrom=scherzer_defocus,
            cs_mm=cs_mm,
            focal_spread_angstrom=focal_spread_angstrom,
            convergence_semiangle_mrad=convergence_semiangle_mrad,
            mode=DoubleCorrectionMode.UNCORRECTED,
        )

    @classmethod
    def cs_corrected(
        cls,
        energy_kev: float = 200.0,
        cs_um: float = 5.0,
        defocus_angstrom: float = -30.0,
        focal_spread_angstrom: float = 25.0,
        convergence_semiangle_mrad: float = 0.2,
        aperture_cutoff_mrad: float | None = 25.0,
    ) -> MicroscopeAberrations:
        """Preset for Cs-corrected TEM with small residual Cs and moderate focal spread."""
        return cls(
            energy_kev=energy_kev,
            defocus_angstrom=defocus_angstrom,
            cs_mm=cs_um * 1e-3,
            focal_spread_angstrom=focal_spread_angstrom,
            convergence_semiangle_mrad=convergence_semiangle_mrad,
            aperture_cutoff_mrad=aperture_cutoff_mrad,
            mode=DoubleCorrectionMode.CS_CORRECTED,
        )

    @classmethod
    def double_corrected(
        cls,
        energy_kev: float = 300.0,
        cs_um: float = 0.0,
        defocus_angstrom: float = -10.0,
        focal_spread_angstrom: float = 5.0,
        convergence_semiangle_mrad: float = 0.1,
        aperture_cutoff_mrad: float | None = 35.0,
    ) -> MicroscopeAberrations:
        """Preset for double-corrected (Cs + Cc) TEM with near-zero Cs and narrow focal spread."""
        return cls(
            energy_kev=energy_kev,
            defocus_angstrom=defocus_angstrom,
            cs_mm=cs_um * 1e-3,
            focal_spread_angstrom=focal_spread_angstrom,
            convergence_semiangle_mrad=convergence_semiangle_mrad,
            aperture_cutoff_mrad=aperture_cutoff_mrad,
            mode=DoubleCorrectionMode.DOUBLE_CORRECTED,
        )

    @classmethod
    def ncsi(
        cls,
        energy_kev: float = 200.0,
        cs_um: float = -15.0,
        defocus_angstrom: float = 50.0,
        focal_spread_angstrom: float = 15.0,
        convergence_semiangle_mrad: float = 0.2,
        aperture_cutoff_mrad: float | None = 30.0,
    ) -> MicroscopeAberrations:
        """Preset for Negative Cs Imaging (NCSI): Cs < 0 and Delta f > 0."""
        return cls(
            energy_kev=energy_kev,
            defocus_angstrom=defocus_angstrom,
            cs_mm=cs_um * 1e-3,
            focal_spread_angstrom=focal_spread_angstrom,
            convergence_semiangle_mrad=convergence_semiangle_mrad,
            aperture_cutoff_mrad=aperture_cutoff_mrad,
            mode=DoubleCorrectionMode.NCSI,
        )

    def wave_aberration(self, q: np.ndarray, theta: np.ndarray | None = None) -> np.ndarray:
        """Evaluate the wave aberration function chi(q, theta) in radians.

        Parameters
        ----------
        q : np.ndarray
            Spatial frequencies in 1/Å.
        theta : np.ndarray, optional
            Azimuthal angles in radians. Defaults to zero for 1D radial evaluation.

        Returns
        -------
        np.ndarray
            Phase aberration chi in radians.
        """
        lam = self.wavelength_angstrom
        # Radially symmetric terms:
        # chi = pi*Delta_f*lambda*q^2 + 0.5*pi*Cs*lambda^3*q^4 + (1/3)*pi*C5*lambda^5*q^6
        chi = (
            math.pi * self.defocus_angstrom * lam * (q**2)
            + 0.5 * math.pi * self.cs_angstrom * (lam**3) * (q**4)
            + (1.0 / 3.0) * math.pi * self.c5_angstrom * (lam**5) * (q**6)
        )
        if theta is not None:
            if self.astigmatism_angstrom != 0.0:
                phi12 = math.radians(self.astigmatism_angle_deg)
                chi = chi + (
                    math.pi
                    * self.astigmatism_angstrom
                    * lam
                    * (q**2)
                    * np.cos(2.0 * (theta - phi12))
                )
            if self.trefoil_angstrom != 0.0:
                phi23 = math.radians(self.trefoil_angle_deg)
                chi = chi + (
                    (2.0 / 3.0)
                    * math.pi
                    * self.trefoil_angstrom
                    * (lam**2)
                    * (q**3)
                    * np.cos(3.0 * (theta - phi23))
                )
            if self.coma_angstrom != 0.0:
                phi21 = math.radians(self.coma_angle_deg)
                chi = chi + (
                    (2.0 / 3.0)
                    * math.pi
                    * self.coma_angstrom
                    * (lam**2)
                    * (q**3)
                    * np.cos(theta - phi21)
                )
        return chi

    def temporal_envelope(self, q: np.ndarray) -> np.ndarray:
        """Evaluate Frank's temporal coherence damping envelope E_c(q).

        E_c(q) = exp(-0.5 * pi^2 * lambda^2 * Delta^2 * q^4)
        """
        lam = self.wavelength_angstrom
        delta = self.focal_spread_angstrom
        arg = -0.5 * (math.pi**2) * (lam**2) * (delta**2) * (q**4)
        return np.exp(np.clip(arg, -100.0, 0.0))

    def spatial_envelope(self, q: np.ndarray) -> np.ndarray:
        """Evaluate Frank's spatial coherence damping envelope E_s(q).

        E_s(q) = exp(-pi^2 * alpha_s^2 * (Delta_f * q + Cs * lambda^2 * q^3)^2)
        """
        alpha_rad = self.convergence_semiangle_mrad * 1e-3
        lam = self.wavelength_angstrom
        df = self.defocus_angstrom
        cs = self.cs_angstrom
        deriv = df * q + cs * (lam**2) * (q**3)
        arg = -(math.pi**2) * (alpha_rad**2) * (deriv**2)
        return np.exp(np.clip(arg, -100.0, 0.0))

    def aperture_mask(self, q: np.ndarray) -> np.ndarray:
        """Evaluate objective aperture transmission mask A(q)."""
        if self.aperture_cutoff_mrad is None or self.aperture_cutoff_mrad <= 0.0:
            return np.ones_like(q, dtype=np.float64)
        q_cutoff = (self.aperture_cutoff_mrad * 1e-3) / self.wavelength_angstrom
        return (q <= q_cutoff).astype(np.float64)

    def evaluate_ctf_1d(
        self,
        max_q_inv_angstrom: float = 2.0,
        num_points: int = 500,
        azimuth_deg: float = 0.0,
    ) -> CTF1D:
        """Compute the Contrast Transfer Function along one azimuth of the back focal plane.

        Purpose
        -------
        Evaluate the damped transfer function T(q) = E(q) sin(chi(q, theta)) along a
        radial cut at a stated azimuth, together with the coherence envelopes, the
        first zero crossing and the information limit.

        When to use
        -----------
        For a round lens — defocus, Cs and C5 only — the cut is azimuth-independent
        and describes the instrument completely. When
        :attr:`has_azimuthal_aberrations` is true this is one section through an
        anisotropic transfer function, and :meth:`evaluate_ctf_azimuthal` gives the
        spread across azimuth that a single cut cannot show.

        Parameters
        ----------
        max_q_inv_angstrom : float
            Upper spatial frequency of the radial grid, in 1/Angstrom.
        num_points : int
            Number of radial samples.
        azimuth_deg : float
            Azimuth theta of the cut, in degrees, measured in the back focal plane
            from the same reference as the aberration azimuths phi12, phi21 and
            phi23. The non-round terms are evaluated at this azimuth.

        Returns
        -------
        CTF1D
            The transfer profile, its envelopes, and the azimuth it was cut at.

        See Also
        --------
        evaluate_ctf_azimuthal : the transfer function over the full azimuth range.
        """
        q = np.linspace(0.0, max_q_inv_angstrom, num_points)
        theta = np.full_like(q, math.radians(azimuth_deg))
        chi = self.wave_aberration(q, theta)
        sin_chi = np.sin(chi)
        ec = self.temporal_envelope(q)
        es = self.spatial_envelope(q)
        ap = self.aperture_mask(q)
        total_env = ec * es * ap
        transfer = total_env * sin_chi

        first_zero = _first_zero_crossing(q, sin_chi)

        # Information limit where total envelope drops to 1/e^2 ~ 0.135
        info_limit_q = float("nan")
        below_thresh = np.where((ec * es) < math.exp(-2.0))[0]
        if len(below_thresh) > 0 and below_thresh[0] > 0:
            idx = below_thresh[0]
            q1, q2 = q[idx - 1], q[idx]
            y1, y2 = (ec * es)[idx - 1], (ec * es)[idx]
            target = math.exp(-2.0)
            info_limit_q = float(q1 + (target - y1) * (q2 - q1) / (y2 - y1))
        elif len(below_thresh) == 0:
            info_limit_q = max_q_inv_angstrom

        return CTF1D(
            spatial_frequencies_inv_angstrom=q,
            phase_shift_rad=chi,
            ctf_undamped=sin_chi,
            temporal_envelope=ec,
            spatial_envelope=es,
            total_envelope=total_env,
            transfer_function=transfer,
            first_zero_q_inv_angstrom=first_zero,
            information_limit_q_inv_angstrom=info_limit_q,
            aberrations=self,
            azimuth_deg=float(azimuth_deg),
        )

    def evaluate_ctf_azimuthal(
        self,
        max_q_inv_angstrom: float = 2.0,
        num_radial: int = 500,
        num_azimuthal: int = 72,
    ) -> AzimuthalCTF:
        """Compute the transfer function over the full azimuth range of the back focal plane.

        Purpose
        -------
        Resolve the anisotropy that two-fold astigmatism, axial coma and trefoil
        impose on phase transfer. Each is evaluated on the same radial grid at
        ``num_azimuthal`` azimuths, and the result carries both the per-azimuth
        profiles and the envelope of best and worst transfer across azimuth.

        When to use
        -----------
        Whenever a residual non-round aberration is nonzero, which for a corrected
        instrument is the normal case: correction drives Cs toward zero and the
        residual terms then set the achievable resolution, differently in different
        directions. For a round lens the band collapses onto the single radial cut
        and this adds nothing over :meth:`evaluate_ctf_1d`.

        Parameters
        ----------
        max_q_inv_angstrom : float
            Upper spatial frequency of the radial grid, in 1/Angstrom.
        num_radial : int
            Number of radial samples per azimuth.
        num_azimuthal : int
            Number of azimuths sampled uniformly over [0, 2*pi).

        Returns
        -------
        AzimuthalCTF
            Per-azimuth transfer, the best/worst transfer band, and the range of
            point resolution over azimuth.

        Notes
        -----
        The coherence envelopes are Frank's isotropic forms, which are derived for a
        round lens: the spatial envelope uses the radial derivative of the round part
        of chi. The anisotropy reported here is therefore the anisotropy of the
        transfer oscillation, not of the damping.
        """
        if num_azimuthal < 1:
            raise ValueError(f"num_azimuthal must be at least 1, got {num_azimuthal}")
        q = np.linspace(0.0, max_q_inv_angstrom, num_radial)
        azimuths = np.linspace(0.0, 2.0 * math.pi, num_azimuthal, endpoint=False)
        # A round lens ignores theta entirely, so broadcast the radial profile back
        # onto the azimuth axis: the grid shape must not depend on which terms are set.
        chi = np.broadcast_to(
            self.wave_aberration(q[None, :], azimuths[:, None]), (azimuths.size, q.size)
        ).copy()
        sin_chi = np.sin(chi)
        envelope = self.temporal_envelope(q) * self.spatial_envelope(q) * self.aperture_mask(q)
        transfer = envelope[None, :] * sin_chi

        first_zeros = np.array(
            [_first_zero_crossing(q, sin_chi[index]) for index in range(azimuths.size)]
        )
        return AzimuthalCTF(
            spatial_frequencies_inv_angstrom=q,
            azimuths_rad=azimuths,
            phase_shift_rad=chi,
            transfer_function=transfer,
            total_envelope=envelope,
            first_zero_q_inv_angstrom=first_zeros,
            aberrations=self,
        )

    def describe(self) -> str:
        """Prose explanation of the objective lens state per explainable-results doctrine."""
        lam = self.wavelength_angstrom
        lines = [
            f"TEM objective lens operating at {self.energy_kev:.1f} kV (relativistic wavelength "
            f"λ = {lam:.5f} Å).",
            f"Regime: {self.mode.value} with defocus Δf = {self.defocus_angstrom:.1f} Å and "
            f"Cs = {self.cs_um:.2f} µm ({self.cs_mm:.4f} mm).",
        ]
        if self.mode == DoubleCorrectionMode.UNCORRECTED and self.cs_mm > 0.0:
            lines.append(
                f"Scherzer defocus is Δf_Sch = {self.scherzer_defocus_angstrom:.1f} Å, "
                f"giving a point resolution of d_Sch = {self.scherzer_resolution_angstrom:.2f} Å."
            )
        elif self.mode == DoubleCorrectionMode.DOUBLE_CORRECTED:
            lines.append(
                "Double-corrected condition (spherical and chromatic correction): Cs is reduced to "
                f"{self.cs_um:.2f} µm and focal spread to {self.focal_spread_angstrom:.1f} Å, "
                "extending the information limit beyond the Scherzer boundary."
            )
        elif self.mode == DoubleCorrectionMode.NCSI:
            lines.append(
                f"Negative Cs Imaging (NCSI) condition: Cs = {self.cs_um:.1f} µm (< 0) "
                f"with overfocus Δf = +{self.defocus_angstrom:.1f} Å, producing bright atom "
                "columns on a dark background with maximum contrast and minimal Fresnel fringing."
            )
        if self.c5_mm != 0.0:
            lines.append(
                f"Fifth-order spherical aberration C5 = {self.c5_mm:.4f} mm is retained; it is "
                "round, so it shifts phase transfer equally at every azimuth."
            )
        residual = self.residual_aberration_terms()
        if residual:
            written = ", ".join(
                f"{symbol} = {amplitude:.1f} Å at azimuth {angle:.1f}°"
                for _, symbol, amplitude, angle in residual
            )
            lines.append(
                f"Residual non-round aberrations are present: {written}. Phase transfer therefore "
                "depends on azimuth, and a single radial cut of the transfer function describes "
                "one direction rather than the lens."
            )
        if self.focal_spread_angstrom > 0.0:
            lines.append(
                f"Coherence damping: focal spread Δ = {self.focal_spread_angstrom:.1f} Å, "
                f"convergence semi-angle α_s = {self.convergence_semiangle_mrad:.2f} mrad."
            )
        if self.aperture_cutoff_mrad is not None:
            lines.append(
                f"Objective aperture: semi-angle cutoff = {self.aperture_cutoff_mrad:.1f} mrad."
            )
        return " ".join(lines)


@dataclass(frozen=True, slots=True)
class CTF1D:
    """1D Contrast Transfer Function and coherence damping profile.

    Attributes
    ----------
    spatial_frequencies_inv_angstrom : np.ndarray
        Spatial frequency grid q in 1/Å.
    phase_shift_rad : np.ndarray
        Wave aberration chi(q) in radians.
    ctf_undamped : np.ndarray
        Undamped transfer sin(chi(q)).
    temporal_envelope : np.ndarray
        Frank's temporal coherence damping envelope Ec(q).
    spatial_envelope : np.ndarray
        Frank's spatial coherence damping envelope Es(q).
    total_envelope : np.ndarray
        Total envelope Ec(q) * Es(q) * A(q).
    transfer_function : np.ndarray
        Damped transfer function T(q) = total_envelope * sin(chi(q)).
    first_zero_q_inv_angstrom : float
        First zero crossing of sin(chi) in 1/Å.
    information_limit_q_inv_angstrom : float
        Spatial frequency where damping reaches exp(-2) ~ 0.135.
    aberrations : MicroscopeAberrations
        Associated lens parameters.
    azimuth_deg : float
        Azimuth of the radial cut, in degrees. Meaningful only when the lens
        carries a non-round aberration; for a round lens every azimuth gives the
        same profile.
    """

    spatial_frequencies_inv_angstrom: np.ndarray
    phase_shift_rad: np.ndarray
    ctf_undamped: np.ndarray
    temporal_envelope: np.ndarray
    spatial_envelope: np.ndarray
    total_envelope: np.ndarray
    transfer_function: np.ndarray
    first_zero_q_inv_angstrom: float
    information_limit_q_inv_angstrom: float
    aberrations: MicroscopeAberrations
    azimuth_deg: float = 0.0

    @property
    def q_inv_angstrom(self) -> np.ndarray:
        """Alias for spatial_frequencies_inv_angstrom."""
        return self.spatial_frequencies_inv_angstrom

    @property
    def transfer(self) -> np.ndarray:
        """Alias for transfer_function."""
        return self.transfer_function

    @property
    def sin_chi(self) -> np.ndarray:
        """Alias for ctf_undamped."""
        return self.ctf_undamped

    @property
    def point_resolution_angstrom(self) -> float:
        """Point resolution d = 1 / q_first_zero in Angstrom."""
        if math.isnan(self.first_zero_q_inv_angstrom) or self.first_zero_q_inv_angstrom <= 0.0:
            return float("nan")
        return 1.0 / self.first_zero_q_inv_angstrom

    @property
    def first_zero_d_spacing_angstrom(self) -> float:
        """Alias for point_resolution_angstrom."""
        return self.point_resolution_angstrom

    @property
    def information_limit_angstrom(self) -> float:
        """Information limit d_info = 1 / q_info in Angstrom."""
        if (
            math.isnan(self.information_limit_q_inv_angstrom)
            or self.information_limit_q_inv_angstrom <= 0.0
        ):
            return float("nan")
        return 1.0 / self.information_limit_q_inv_angstrom

    @property
    def information_limit_d_spacing_angstrom(self) -> float:
        """Alias for information_limit_angstrom."""
        return self.information_limit_angstrom

    def describe(self) -> str:
        """Scientific prose description of the CTF properties."""
        lines = [
            f"1D Contrast Transfer Function computed at {self.aberrations.energy_kev:.1f} kV with "
            f"defocus Δf = {self.aberrations.defocus_angstrom:.1f} Å and "
            f"Cs = {self.aberrations.cs_um:.2f} µm."
        ]
        if not math.isnan(self.point_resolution_angstrom):
            lines.append(
                f"The first zero crossing sits at q = {self.first_zero_q_inv_angstrom:.3f} Å⁻¹, "
                f"corresponding to a point resolution of {self.point_resolution_angstrom:.3f} Å."
            )
        else:
            lines.append(
                "No zero crossing occurs within the passband, indicating a non-oscillating "
                "or aperture-limited transfer regime."
            )
        if not math.isnan(self.information_limit_angstrom):
            lines.append(
                "The 1/e² information limit sits at q = "
                f"{self.information_limit_q_inv_angstrom:.3f} Å⁻¹ "
                f"(d_info = {self.information_limit_angstrom:.3f} Å), governed by the focal "
                f"spread of {self.aberrations.focal_spread_angstrom:.1f} Å."
            )
        residual = self.aberrations.residual_aberration_terms()
        if residual:
            written = ", ".join(
                f"{symbol} = {amplitude:.1f} Å at {angle:.1f}°"
                for _, symbol, amplitude, angle in residual
            )
            lines.append(
                f"This profile is the cut at azimuth θ = {self.azimuth_deg:.1f}°, not the whole "
                f"lens: the non-round residual aberrations {written} make phase transfer depend on "
                "azimuth, so the point resolution above applies to this direction alone."
            )
        return " ".join(lines)


def _first_zero_crossing(q: np.ndarray, sin_chi: np.ndarray) -> float:
    """Locate the first zero of sin(chi) above q = 0.05 1/Angstrom by linear interpolation.

    The trivial zero at the origin is excluded: chi(0) = 0 for every lens, and it
    carries no information about resolution. Returns NaN when the profile does not
    cross zero inside the sampled range, which is the non-oscillating regime of a
    well-corrected lens rather than a failure.
    """
    sign_changes = np.where(np.diff(np.sign(sin_chi)))[0]
    for idx in sign_changes:
        if q[idx] <= 0.05:
            continue
        q1, q2 = float(q[idx]), float(q[idx + 1])
        y1, y2 = float(sin_chi[idx]), float(sin_chi[idx + 1])
        if y2 == y1:
            return q1
        return q1 - y1 * (q2 - q1) / (y2 - y1)
    return float("nan")


@dataclass(frozen=True, slots=True)
class AzimuthalCTF:
    """Contrast transfer resolved over azimuth for a lens with non-round aberrations.

    Purpose
    -------
    A radial cut of the transfer function describes a round lens completely. Once
    two-fold astigmatism, axial coma or trefoil is present, chi depends on azimuth
    and the achievable resolution differs by direction — which is what an operator
    tuning a corrector is looking at when they judge the symmetry of a Thon-ring
    pattern. This object holds the transfer at each sampled azimuth together with
    the band it sweeps out.

    Attributes
    ----------
    spatial_frequencies_inv_angstrom : np.ndarray
        Radial grid q in 1/Angstrom, shape ``(num_radial,)``.
    azimuths_rad : np.ndarray
        Sampled azimuths in radians over [0, 2*pi), shape ``(num_azimuthal,)``.
    phase_shift_rad : np.ndarray
        Wave aberration chi(q, theta) in radians, shape ``(num_azimuthal, num_radial)``.
    transfer_function : np.ndarray
        Damped transfer at each azimuth, shape ``(num_azimuthal, num_radial)``.
    total_envelope : np.ndarray
        Coherence and aperture envelope, shape ``(num_radial,)``. Frank's envelopes
        are isotropic, so one radial profile applies at every azimuth.
    first_zero_q_inv_angstrom : np.ndarray
        First zero crossing per azimuth in 1/Angstrom, shape ``(num_azimuthal,)``.
        NaN where the cut does not cross zero within the sampled range.
    aberrations : MicroscopeAberrations
        The lens state these profiles were computed from.
    """

    spatial_frequencies_inv_angstrom: np.ndarray
    azimuths_rad: np.ndarray
    phase_shift_rad: np.ndarray
    transfer_function: np.ndarray
    total_envelope: np.ndarray
    first_zero_q_inv_angstrom: np.ndarray
    aberrations: MicroscopeAberrations

    @property
    def azimuths_deg(self) -> np.ndarray:
        """Sampled azimuths in degrees."""
        return np.asarray(np.degrees(self.azimuths_rad), dtype=np.float64)

    @property
    def transfer_min(self) -> np.ndarray:
        """Lowest transfer over azimuth at each spatial frequency."""
        return np.asarray(np.min(self.transfer_function, axis=0), dtype=np.float64)

    @property
    def transfer_max(self) -> np.ndarray:
        """Highest transfer over azimuth at each spatial frequency."""
        return np.asarray(np.max(self.transfer_function, axis=0), dtype=np.float64)

    @property
    def point_resolution_angstrom(self) -> np.ndarray:
        """Point resolution per azimuth, 1/q_first_zero, in Angstrom."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(
                self.first_zero_q_inv_angstrom > 0.0,
                1.0 / self.first_zero_q_inv_angstrom,
                np.nan,
            )

    @property
    def point_resolution_range_angstrom(self) -> tuple[float, float]:
        """Best and worst point resolution over azimuth, in Angstrom.

        Both entries are NaN when no azimuth produces a zero crossing inside the
        sampled range.
        """
        values = self.point_resolution_angstrom
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return (float("nan"), float("nan"))
        return (float(np.min(finite)), float(np.max(finite)))

    @property
    def resolution_anisotropy_angstrom(self) -> float:
        """Spread between the best and worst point resolution over azimuth, in Angstrom.

        Zero for a round lens. It is the single number that says how much the
        residual non-round aberrations cost, and in a corrected instrument it is
        usually the quantity worth minimising.
        """
        best, worst = self.point_resolution_range_angstrom
        if math.isnan(best):
            return float("nan")
        return worst - best

    @property
    def worst_azimuth_deg(self) -> float:
        """Azimuth of the coarsest point resolution, in degrees; NaN if none resolves."""
        values = self.point_resolution_angstrom
        if not np.any(np.isfinite(values)):
            return float("nan")
        return float(self.azimuths_deg[int(np.nanargmax(values))])

    def describe(self) -> str:
        """Scientific prose description of the azimuthal transfer behaviour."""
        residual = self.aberrations.residual_aberration_terms()
        best, worst = self.point_resolution_range_angstrom
        if not residual:
            return (
                "The lens carries no non-round aberration, so phase transfer is isotropic: every "
                "azimuth gives the same profile and the point resolution of "
                f"{best:.3f} Å applies in all directions."
                if not math.isnan(best)
                else (
                    "The lens carries no non-round aberration, so phase transfer is isotropic. No "
                    "zero crossing occurs within the sampled passband."
                )
            )
        written = ", ".join(
            f"{symbol} = {amplitude:.1f} Å at {angle:.1f}°"
            for _, symbol, amplitude, angle in residual
        )
        lines = [
            f"Phase transfer is anisotropic: the residual aberrations {written} make the wave "
            f"aberration depend on azimuth, evaluated here at {self.azimuths_rad.size} azimuths."
        ]
        if math.isnan(best):
            lines.append(
                "No azimuth produces a zero crossing inside the sampled passband, so no "
                "directional point resolution can be quoted over this range."
            )
        else:
            lines.append(
                f"Point resolution ranges from {best:.3f} Å at best to {worst:.3f} Å at worst "
                f"(azimuth {self.worst_azimuth_deg:.1f}°), an anisotropy of "
                f"{self.resolution_anisotropy_angstrom:.3f} Å."
            )
        lines.append(
            "The coherence envelopes are Frank's isotropic forms, so this anisotropy is that of "
            "the transfer oscillation, not of the damping."
        )
        return " ".join(lines)


@dataclass(frozen=True, slots=True)
class AtomicSnapshot:
    """Atomic coordinates, chemical species, and simulation bounding cell.

    Represents a discrete atomic configuration suitable for multislice HREM
    simulation: single crystals, crystals containing vacancies or dislocations,
    amorphous foils, or snapshots imported from MD simulations (XYZ, POSCAR, CIF).

    Attributes
    ----------
    species : tuple of str
        Chemical element symbol for each atom (length N).
    positions : np.ndarray
        Cartesian coordinates in Angstrom (N, 3).
    cell : np.ndarray
        Simulation cell vectors as a 3x3 matrix [[ax, ay, az], [bx, by, bz], [cx, cy, cz]] in Å.
    periodicity : tuple of bool
        Periodic boundary conditions along x, y, z (default True, True, False for thin foil).
    label : str
        Descriptive label of the snapshot.
    """

    species: tuple[str, ...]
    positions: np.ndarray
    cell: np.ndarray
    periodicity: tuple[bool, bool, bool] = (True, True, False)
    label: str = "Atomic Snapshot"

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", tuple(self.species))
        pos = as_float_array(self.positions)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"positions must have shape (N, 3), got {pos.shape}.")
        if len(self.species) != pos.shape[0]:
            raise ValueError(
                f"species length ({len(self.species)}) must match positions count ({pos.shape[0]})."
            )
        object.__setattr__(self, "positions", pos)
        cell_mat = as_float_array(self.cell)
        if cell_mat.shape != (3, 3):
            raise ValueError(f"cell must have shape (3, 3), got {cell_mat.shape}.")
        object.__setattr__(self, "cell", cell_mat)

    @property
    def natoms(self) -> int:
        """Total number of atoms in the snapshot."""
        return len(self.species)

    @property
    def dimensions_angstrom(self) -> tuple[float, float, float]:
        """Dimensions of the simulation box (Lx, Ly, Lz) in Angstrom."""
        return (
            float(np.linalg.norm(self.cell[0])),
            float(np.linalg.norm(self.cell[1])),
            float(np.linalg.norm(self.cell[2])),
        )

    @property
    def thickness_angstrom(self) -> float:
        """Specimen thickness along beam direction (Z) in Angstrom."""
        if self.natoms == 0:
            return 0.0
        z = self.positions[:, 2]
        span = float(np.max(z) - np.min(z))
        # Add atomic sphere margin (~ 2 Å) or use cell height if bounded
        cell_z = float(self.cell[2, 2])
        return cell_z if cell_z > span else span + 2.0

    @classmethod
    def from_phase(
        cls,
        phase: Phase,
        supercell: tuple[int, int, int] = (2, 2, 4),
        zone_axis: tuple[int, int, int] = (0, 0, 1),
        label: str | None = None,
    ) -> AtomicSnapshot:
        """Build an oriented supercell snapshot from a PyTex Phase.

        Parameters
        ----------
        phase : Phase
            Crystallographic phase with unit cell and atomic sites.
        supercell : tuple of (int, int, int)
            Number of unit cell repetitions (nx, ny, nz).
        zone_axis : tuple of (int, int, int)
            Crystallographic zone axis [uvw] to align along the electron beam (+Z).
        label : str, optional
            Snapshot descriptor.
        """
        if phase.unit_cell is None or not phase.unit_cell.sites:
            raise ValueError(
                f"Phase '{phase.name}' must have a populated unit_cell with atomic sites."
            )

        basis = phase.lattice.direct_basis()
        a_vec = basis.vector(0)
        b_vec = basis.vector(1)
        c_vec = basis.vector(2)
        basis_matrix = np.vstack([a_vec, b_vec, c_vec])

        # Generate fractional coordinates in supercell
        nx, ny, nz = supercell
        sites = phase.unit_cell.sites
        species_list: list[str] = []
        positions_list: list[np.ndarray] = []

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    shift = np.array([ix, iy, iz], dtype=np.float64)
                    for site in sites:
                        frac = site.fractional_coordinates + shift
                        # Cartesian pos before supercell scaling
                        cart = frac @ basis_matrix
                        species_list.append(site.species)
                        positions_list.append(cart)

        positions = np.array(positions_list, dtype=np.float64)
        super_cell = np.vstack([a_vec * nx, b_vec * ny, c_vec * nz])

        # If zone_axis != (0, 0, 1), rotate so zone axis points along +Z
        za = np.asarray(zone_axis, dtype=np.float64)
        za_cart = za @ basis_matrix
        za_norm = za_cart / np.linalg.norm(za_cart)
        target_z = np.array([0.0, 0.0, 1.0], dtype=np.float64)

        dot = float(np.clip(np.dot(za_norm, target_z), -1.0, 1.0))
        if abs(dot - 1.0) > 1e-6:
            if abs(dot + 1.0) < 1e-6:
                rot_mat = np.diag([1.0, -1.0, -1.0])
            else:
                rot_axis = np.cross(za_norm, target_z)
                rot_axis = rot_axis / np.linalg.norm(rot_axis)
                angle = math.acos(dot)
                k = np.array(
                    [
                        [0.0, -rot_axis[2], rot_axis[1]],
                        [rot_axis[2], 0.0, -rot_axis[0]],
                        [-rot_axis[1], rot_axis[0], 0.0],
                    ]
                )
                rot_mat = np.eye(3) + math.sin(angle) * k + (1.0 - math.cos(angle)) * (k @ k)
            positions = positions @ rot_mat.T
            super_cell = super_cell @ rot_mat.T

        # Shift positions so min x, y, z are non-negative with a small buffer
        min_pos = np.min(positions, axis=0)
        positions -= min_pos - np.array([0.5, 0.5, 0.5])
        max_pos = np.max(positions, axis=0)
        box_cell = np.diag(max_pos + 1.0)

        za_str = f"[{zone_axis[0]} {zone_axis[1]} {zone_axis[2]}]"
        desc = label or f"{phase.name} supercell ({nx}x{ny}x{nz}) along {za_str}"
        return cls(
            species=tuple(species_list),
            positions=positions,
            cell=box_cell,
            periodicity=(True, True, False),
            label=desc,
        )

    @classmethod
    def crystalline_with_vacancy(
        cls,
        phase: Phase,
        supercell: tuple[int, int, int] = (3, 3, 4),
        zone_axis: tuple[int, int, int] = (0, 0, 1),
        vacancy_count: int = 1,
    ) -> AtomicSnapshot:
        """Create a crystalline snapshot with one or more point vacancies."""
        base = cls.from_phase(phase, supercell=supercell, zone_axis=zone_axis)
        n = base.natoms
        if vacancy_count >= n:
            raise ValueError(f"Cannot create {vacancy_count} vacancies in a system with {n} atoms.")

        # Remove atom closest to center
        center = np.mean(base.positions, axis=0)
        dists = np.linalg.norm(base.positions - center, axis=1)
        remove_indices = np.argsort(dists)[:vacancy_count]

        keep_mask = np.ones(n, dtype=bool)
        keep_mask[remove_indices] = False

        new_species = [s for i, s in enumerate(base.species) if keep_mask[i]]
        new_positions = base.positions[keep_mask]
        return cls(
            species=tuple(new_species),
            positions=new_positions,
            cell=base.cell,
            periodicity=base.periodicity,
            label=f"{base.label} with {vacancy_count} vacancy",
        )

    @classmethod
    def crystalline_with_dislocation(
        cls,
        phase: Phase,
        supercell: tuple[int, int, int] = (4, 4, 3),
        burgers_vector_angstrom: float = 2.5,
        dislocation_type: str = "edge",
    ) -> AtomicSnapshot:
        """Create a crystalline slab with a Volterra dislocation line along Z.

        Parameters
        ----------
        phase : Phase
            Parent phase.
        supercell : tuple
            Supercell repetitions.
        burgers_vector_angstrom : float
            Burgers vector magnitude b in Angstrom.
        dislocation_type : {"edge", "screw"}
            Type of dislocation.
        """
        base = cls.from_phase(phase, supercell=supercell)
        pos = base.positions.copy()
        center = np.mean(pos[:, :2], axis=0)
        x = pos[:, 0] - center[0]
        y = pos[:, 1] - center[1]
        theta = np.arctan2(y, x)
        b = burgers_vector_angstrom
        nu = 0.33  # Poisson ratio

        if dislocation_type == "screw":
            # u_z = b * theta / (2 * pi)
            pos[:, 2] += b * theta / (2.0 * math.pi)
        else:
            # Edge dislocation along Z with Burgers vector along X:
            # u_x = (b / (2*pi)) * [theta + sin(2*theta)/(4*(1-nu))]
            # u_y = -(b / (2*pi)) * [(1-2*nu)/(2*(1-nu))*ln(r) + cos(2*theta)/(4*(1-nu))]
            r = np.clip(np.hypot(x, y), 0.5, None)
            factor = b / (2.0 * math.pi)
            pos[:, 0] += factor * (theta + np.sin(2.0 * theta) / (4.0 * (1.0 - nu)))
            pos[:, 1] -= factor * (
                ((1.0 - 2.0 * nu) / (2.0 * (1.0 - nu))) * np.log(r)
                + np.cos(2.0 * theta) / (4.0 * (1.0 - nu))
            )

        return cls(
            species=base.species,
            positions=pos,
            cell=base.cell,
            periodicity=(False, False, False),
            label=f"{base.label} with {dislocation_type} dislocation (b={b:.2f} Å)",
        )

    @classmethod
    def amorphous_sample(
        cls,
        species: str = "C",
        density_g_cm3: float = 2.0,
        dimensions_angstrom: tuple[float, float, float] = (20.0, 20.0, 15.0),
        min_distance_angstrom: float = 1.45,
        seed: int = 42,
    ) -> AtomicSnapshot:
        """Generate a dense random-packed amorphous sample snapshot.

        Simulates an amorphous foil (such as amorphous carbon or silicon)
        widely used for microscope alignment and Thon ring CTF calibration.

        Parameters
        ----------
        species : str
            Element symbol (e.g. 'C', 'Si').
        density_g_cm3 : float
            Target mass density in g/cm³.
        dimensions_angstrom : tuple of (Lx, Ly, Lz)
            Box size in Angstrom.
        min_distance_angstrom : float
            Hard-sphere exclusion distance preventing unphysical atomic overlap.
        seed : int
            RNG seed for reproducibility.
        """
        rng = np.random.default_rng(seed)
        lx, ly, lz = dimensions_angstrom
        volume_angstrom3 = lx * ly * lz
        volume_cm3 = volume_angstrom3 * 1e-24

        # Atomic mass from atomic number approximation
        atomic_weights = {"H": 1.008, "C": 12.011, "Si": 28.085, "Ge": 72.63, "Au": 196.97}
        m_atom_g = atomic_weights.get(species, 12.0) / 6.02214076e23
        target_natoms = round((density_g_cm3 * volume_cm3) / m_atom_g)
        target_natoms = max(10, min(target_natoms, 2000))

        positions: list[np.ndarray] = []
        r2_min = min_distance_angstrom**2
        max_attempts = target_natoms * 100
        attempts = 0

        while len(positions) < target_natoms and attempts < max_attempts:
            attempts += 1
            pt = rng.uniform(
                low=[0.5, 0.5, 0.5],
                high=[lx - 0.5, ly - 0.5, lz - 0.5],
            )
            if positions:
                cur = np.array(positions)
                # Check distance with periodic wrapping in XY
                dx = np.abs(cur[:, 0] - pt[0])
                dx = np.minimum(dx, lx - dx)
                dy = np.abs(cur[:, 1] - pt[1])
                dy = np.minimum(dy, ly - dy)
                dz = np.abs(cur[:, 2] - pt[2])
                dist2 = dx**2 + dy**2 + dz**2
                if np.min(dist2) < r2_min:
                    continue
            positions.append(pt)

        pos_arr = np.array(positions, dtype=np.float64)
        cell_mat = np.diag([lx, ly, lz])
        return cls(
            species=tuple([species] * len(positions)),
            positions=pos_arr,
            cell=cell_mat,
            periodicity=(True, True, False),
            label=f"Amorphous {species} ({len(positions)} atoms, {density_g_cm3:.2f} g/cm³)",
        )

    @classmethod
    def from_xyz(cls, text_or_path: str | Path, cell: np.ndarray | None = None) -> AtomicSnapshot:
        """Load one XYZ snapshot, preserving its comment and Cartesian angstrom coordinates.

        Use this for an atomic configuration exported by another tool. Supply a
        multiline XYZ string or a filesystem path, and optionally a 3x3 cell in
        angstroms. Without a cell, a bounding box is inferred; XYZ does not define
        lattice periodicity. Returns an :class:`AtomicSnapshot`. A blank comment is
        valid, but nonpositive counts, truncated/multiple frames and nonfinite
        coordinates raise ``ValueError``. Missing paths raise ``FileNotFoundError``.
        """
        if isinstance(text_or_path, Path) or (
            "\n" not in text_or_path and "\r" not in text_or_path
        ):
            content = Path(text_or_path).read_text(encoding="utf-8")
        else:
            content = text_or_path
        lines = content.splitlines()
        if len(lines) < 3:
            raise ValueError("XYZ data must contain at least 3 lines (count, comment, atom).")
        natoms = int(lines[0])
        atom_lines = [line for line in lines[2:] if line.strip()]
        if natoms <= 0 or len(atom_lines) != natoms:
            raise ValueError("XYZ atom count must be positive and match exactly one frame.")
        comment = lines[1].strip()
        species_list: list[str] = []
        pos_list: list[list[float]] = []
        for line in atom_lines:
            parts = line.split()
            if len(parts) < 4:
                raise ValueError("Each XYZ atom needs a species and three coordinates.")
            species_list.append(parts[0])
            pos_list.append([float(parts[1]), float(parts[2]), float(parts[3])])
        pos_arr = np.array(pos_list, dtype=np.float64)
        if not np.all(np.isfinite(pos_arr)):
            raise ValueError("XYZ coordinates must be finite.")
        if cell is not None:
            box = as_float_array(cell)
        else:
            min_p = np.min(pos_arr, axis=0)
            max_p = np.max(pos_arr, axis=0)
            box = np.diag(np.maximum(max_p - min_p + 2.0, [5.0, 5.0, 5.0]))
        return cls(
            species=tuple(species_list),
            positions=pos_arr,
            cell=box,
            periodicity=(True, True, False),
            label=comment or "Imported XYZ snapshot",
        )

    def to_xyz(self) -> str:
        """Export snapshot to standard XYZ file format string."""
        lines = [str(self.natoms), self.label]
        for s, p in zip(self.species, self.positions, strict=True):
            lines.append(f"{s:<3} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
        return "\n".join(lines) + "\n"

    def to_ase(self) -> ase.Atoms:
        """Convert snapshot to an ASE Atoms object."""
        from ase import Atoms

        return Atoms(
            symbols=list(self.species),
            positions=self.positions,
            cell=self.cell,
            pbc=self.periodicity,
        )

    @classmethod
    def from_ase(cls, atoms: ase.Atoms, label: str = "ASE Snapshot") -> AtomicSnapshot:
        """Create snapshot from an ASE Atoms object."""
        return cls(
            species=tuple(atoms.get_chemical_symbols()),
            positions=np.asarray(atoms.get_positions(), dtype=np.float64),
            cell=np.asarray(atoms.get_cell(), dtype=np.float64),
            periodicity=tuple(bool(p) for p in atoms.get_pbc()),  # type: ignore[arg-type]
            label=label,
        )

    def describe(self) -> str:
        """Convention-explicit description of the atomic snapshot."""
        dims = self.dimensions_angstrom
        unique_species = sorted(set(self.species))
        counts = {s: self.species.count(s) for s in unique_species}
        comp_str = ", ".join(f"{s}: {counts[s]}" for s in unique_species)
        return (
            f"{self.label}: {self.natoms} atoms ({comp_str}) in a "
            f"{dims[0]:.2f} x {dims[1]:.2f} x {dims[2]:.2f} Å cell. "
            f"Foil thickness: {self.thickness_angstrom:.2f} Å. "
            f"PBC: ({self.periodicity[0]}, {self.periodicity[1]}, {self.periodicity[2]})."
        )


@dataclass(frozen=True, slots=True)
class HREMSimulationResult:
    """Simulated High-Resolution TEM image intensity, exit wave, and CTF analytics.

    Attributes
    ----------
    image : np.ndarray
        2D simulated intensity image array (Ny, Nx).
    exit_wave : np.ndarray | None
        2D complex exit wavefunction at the specimen bottom plane (Ny, Nx).
    pixel_size_angstrom : float
        Sampling interval dx = dy in Angstrom per pixel.
    extent_angstrom : tuple of (float, float)
        Real-space field of view (Lx, Ly) in Angstrom.
    power_spectrum : np.ndarray
        2D log-intensity FFT power spectrum |F(q)|^2 (Thon rings).
    ctf : CTF1D
        1D contrast transfer function profile and coherence envelopes.
    aberrations : MicroscopeAberrations
        Lens optical settings used during simulation.
    snapshot : AtomicSnapshot
        Specimen snapshot that was simulated.
    """

    image: np.ndarray
    exit_wave: np.ndarray | None
    pixel_size_angstrom: float
    extent_angstrom: tuple[float, float]
    power_spectrum: np.ndarray
    ctf: CTF1D
    aberrations: MicroscopeAberrations
    snapshot: AtomicSnapshot

    @property
    def min_intensity(self) -> float:
        """Minimum pixel intensity in the image."""
        return float(np.min(self.image))

    @property
    def max_intensity(self) -> float:
        """Maximum pixel intensity in the image."""
        return float(np.max(self.image))

    @property
    def mean_intensity(self) -> float:
        """Mean pixel intensity in the image."""
        return float(np.mean(self.image))

    @property
    def weber_contrast(self) -> float:
        """Weber contrast (I_max - I_min) / I_mean."""
        mean_val = self.mean_intensity
        if mean_val <= 0.0:
            return 0.0
        return float((self.max_intensity - self.min_intensity) / mean_val)

    @property
    def michelson_contrast(self) -> float:
        """Michelson contrast (I_max - I_min) / (I_max + I_min)."""
        denom = self.max_intensity + self.min_intensity
        if denom <= 0.0:
            return 0.0
        return float((self.max_intensity - self.min_intensity) / denom)

    @property
    def point_resolution_angstrom(self) -> float:
        """First-zero point resolution in Angstrom from CTF."""
        return self.ctf.point_resolution_angstrom

    @property
    def information_limit_angstrom(self) -> float:
        """Information limit in Angstrom from CTF envelope damping."""
        return self.ctf.information_limit_angstrom

    def line_profile(
        self,
        start_px: tuple[int, int] | None = None,
        end_px: tuple[int, int] | None = None,
        num_points: int = 100,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract a 1D real-space intensity line profile across the image.

        Parameters
        ----------
        start_px : tuple of (int, int), optional
            (y0, x0) starting pixel coordinate. Defaults to (Ny//2, 0).
        end_px : tuple of (int, int), optional
            (y1, x1) ending pixel coordinate. Defaults to (Ny//2, Nx-1).
        num_points : int
            Number of interpolated points along the line.

        Returns
        -------
        distance_angstrom : np.ndarray
            Distance from start in Angstrom.
        intensity : np.ndarray
            Sampled intensity values.
        """
        ny, nx = self.image.shape
        y0, x0 = start_px if start_px is not None else (ny // 2, 0)
        y1, x1 = end_px if end_px is not None else (ny // 2, nx - 1)

        y_coords = np.linspace(y0, y1, num_points)
        x_coords = np.linspace(x0, x1, num_points)
        # Bilinear interpolation
        y_int = np.clip(y_coords.astype(int), 0, ny - 2)
        x_int = np.clip(x_coords.astype(int), 0, nx - 2)
        dy = y_coords - y_int
        dx = x_coords - x_int

        intensities = (
            (1.0 - dy) * (1.0 - dx) * self.image[y_int, x_int]
            + dy * (1.0 - dx) * self.image[y_int + 1, x_int]
            + (1.0 - dy) * dx * self.image[y_int, x_int + 1]
            + dy * dx * self.image[y_int + 1, x_int + 1]
        )
        dist_px = np.hypot(x_coords - x0, y_coords - y0)
        return dist_px * self.pixel_size_angstrom, intensities

    def to_png_base64(self, colormap: str = "gray") -> str:
        """Render the simulated image as a base64-encoded PNG data URL."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.axis("off")
        ax.imshow(
            self.image,
            cmap=colormap,
            origin="lower",
            extent=(0, self.extent_angstrom[0], 0, self.extent_angstrom[1]),
        )
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("ascii")
        return f"data:image/png;base64,{encoded}"

    def to_power_spectrum_base64(self) -> str:
        """Render the 2D FFT power spectrum (Thon rings) as a base64 PNG data URL."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax.axis("off")
        max_q = 0.5 / self.pixel_size_angstrom
        ax.imshow(
            self.power_spectrum,
            cmap="inferno",
            origin="lower",
            extent=(-max_q, max_q, -max_q, max_q),
        )
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("ascii")
        return f"data:image/png;base64,{encoded}"

    def describe(self) -> str:
        """Convention-explicit description of the simulated HREM result."""
        lx, ly = self.extent_angstrom
        ny, nx = self.image.shape
        contrast_pct = self.michelson_contrast * 100.0
        lines = [
            f"Simulated High-Resolution TEM image for {self.snapshot.label}.",
            (
                f"Microscope: {self.aberrations.energy_kev:.1f} kV, "
                f"{self.aberrations.mode.value} regime, "
                f"defocus Δf = {self.aberrations.defocus_angstrom:.1f} Å, "
                f"Cs = {self.aberrations.cs_um:.2f} µm."
            ),
            (
                f"Grid: {nx} x {ny} pixels (pixel size: {self.pixel_size_angstrom:.3f} Å/px, "
                f"field of view: {lx:.1f} x {ly:.1f} Å)."
            ),
            (
                f"Image intensity: min = {self.min_intensity:.3f}, "
                f"max = {self.max_intensity:.3f}, "
                f"Michelson contrast = {contrast_pct:.1f}%."
            ),
        ]
        if not math.isnan(self.ctf.point_resolution_angstrom):
            lines.append(f"Point resolution: {self.ctf.point_resolution_angstrom:.3f} Å.")
        if not math.isnan(self.ctf.information_limit_angstrom):
            lines.append(f"Information limit: {self.ctf.information_limit_angstrom:.3f} Å.")
        return " ".join(lines)


def pure_python_phase_object_simulation(
    snapshot: AtomicSnapshot,
    aberrations: MicroscopeAberrations,
    sampling_angstrom: float = 0.1,
) -> HREMSimulationResult:
    """Pure-Python phase-object HREM simulation fallback.

    Evaluates the projected potential using Kirkland parameterized electron
    scattering factors, phase grating transmission t(x,y) = exp(i * sigma * V_p(x,y)),
    Fourier-space objective lens CTF filtering, and intensity formation.

    Parameters
    ----------
    snapshot : AtomicSnapshot
        Specimen atoms.
    aberrations : MicroscopeAberrations
        Objective lens aberrations.
    sampling_angstrom : float
        Real-space pixel size in Angstrom.

    Returns
    -------
    HREMSimulationResult
        Computed HREM result.
    """
    lx, ly = snapshot.dimensions_angstrom[:2]
    nx = max(32, round(lx / sampling_angstrom))
    ny = max(32, round(ly / sampling_angstrom))
    dx = lx / nx
    dy = ly / ny

    # Kirkland Gaussian approximation for projected atomic potential V_p(r):
    # V_p(r) ~ Z^(2/3) * exp(-r^2 / (2*r0^2))
    x_grid = np.linspace(0.0, lx, nx, endpoint=False)
    y_grid = np.linspace(0.0, ly, ny, endpoint=False)
    xx, yy = np.meshgrid(x_grid, y_grid)
    v_proj = np.zeros((ny, nx), dtype=np.float64)

    atomic_z = {
        "H": 1,
        "C": 6,
        "N": 7,
        "O": 8,
        "Al": 13,
        "Si": 14,
        "Fe": 26,
        "Ni": 28,
        "Zr": 40,
        "Au": 79,
    }

    for s, p in zip(snapshot.species, snapshot.positions, strict=True):
        z_num = atomic_z.get(s, 14)
        peak_v = 40.0 * (z_num**0.7)  # in Volt*Angstrom
        width2 = 0.35**2  # radius squared in Angstrom^2

        # Compute minimum distance with periodic wrapping in x and y
        rx = np.abs(xx - p[0])
        rx = np.minimum(rx, lx - rx)
        ry = np.abs(yy - p[1])
        ry = np.minimum(ry, ly - ry)
        r2 = rx**2 + ry**2
        mask = r2 < (3.0 * 0.35) ** 2
        v_proj[mask] += peak_v * np.exp(-r2[mask] / (2.0 * width2))

    # Phase object approximation: psi_exit = exp(i * sigma * V_proj)
    sigma = relativistic_interaction_parameter_inv_v_angstrom(aberrations.energy_kev)
    exit_wave = np.exp(1j * sigma * v_proj)

    # Fourier transform to reciprocal space
    psi_q = np.fft.fft2(exit_wave)
    qx = np.fft.fftfreq(nx, d=dx)
    qy = np.fft.fftfreq(ny, d=dy)
    qxx, qyy = np.meshgrid(qx, qy)
    qq = np.hypot(qxx, qyy)
    theta = np.arctan2(qyy, qxx)

    # Evaluate CTF
    chi = aberrations.wave_aberration(qq, theta)
    ec = aberrations.temporal_envelope(qq)
    es = aberrations.spatial_envelope(qq)
    ap = aberrations.aperture_mask(qq)
    # Objective transfer function H(q) = A(q) * E(q) * exp(-i * chi(q))
    h_q = ap * ec * es * np.exp(-1j * chi)

    # Filtered wave in image plane
    image_wave = np.fft.ifft2(psi_q * h_q)
    intensity = np.real(image_wave * np.conj(image_wave))

    # 2D FFT Power Spectrum (Thon rings)
    ps = np.fft.fftshift(np.log10(np.abs(np.fft.fft2(intensity)) ** 2 + 1e-6))
    # Normalize power spectrum to [0, 1]
    ps = (ps - np.min(ps)) / (np.max(ps) - np.min(ps) + 1e-12)

    ctf_1d = aberrations.evaluate_ctf_1d(max_q_inv_angstrom=float(np.max(qq)))

    return HREMSimulationResult(
        image=intensity,
        exit_wave=exit_wave,
        pixel_size_angstrom=float(dx),
        extent_angstrom=(float(lx), float(ly)),
        power_spectrum=ps,
        ctf=ctf_1d,
        aberrations=aberrations,
        snapshot=snapshot,
    )
