from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_int_array
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import Phase
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.xrd import RadiationSpec

ScatteringModelName = Literal["unit", "atomic_number", "xray_gaussian_proxy"]
IntensityCorrectionName = Literal["none", "lorentz_polarization"]


def lorentz_polarization_factor(two_theta_rad: float) -> float:
    """The combined Lorentz-polarization factor for an unpolarized beam.

    Purpose
    -------
    Powder intensities are not proportional to ``|F|^2`` alone: the time a
    reflection spends in the diffracting condition and the polarization of
    the scattered beam both depend on angle. This is the standard correction
    ``(1 + cos^2(2*theta)) / (sin^2(theta) cos(theta))`` for an unpolarized
    incident beam and a conventional powder diffractometer without a
    monochromator.

    Parameters
    ----------
    two_theta_rad : float
        Scattering angle in radians, strictly between 0 and ``pi``. Values at
        the ends raise, since the factor diverges there.

    Notes
    -----
    A monochromated instrument needs a different polarization term; this
    function does not model one.
    """

    theta = 0.5 * float(two_theta_rad)
    if not np.isfinite(two_theta_rad) or two_theta_rad <= 0.0 or two_theta_rad >= np.pi:
        raise ValueError("two_theta_rad must be finite and satisfy 0 < two_theta_rad < pi.")
    sin_theta = max(float(np.sin(theta)), 1e-8)
    cos_theta = max(float(np.cos(theta)), 1e-8)
    cos_two_theta = float(np.cos(two_theta_rad))
    return float((1.0 + cos_two_theta * cos_two_theta) / (sin_theta * sin_theta * cos_theta))


@dataclass(frozen=True, slots=True)
class ScatteringFactorTable:
    """The atomic scattering model used for structure-factor calculations.

    Purpose
    -------
    Selects how scattering power varies with species and scattering angle.

    Limits
    ------
    The available models are deliberate proxies, not tabulated Cromer-Mann or
    Doyle-Turner factors: ``"unit"`` counts sites, ``"atomic_number"`` uses
    ``Z`` with no angular falloff, and the remaining model applies a smooth
    monotonic decay. They give correct systematic absences and correct
    relative intensities among reflections of similar angle, but they are not
    quantitative across a wide angular range.

    Attributes
    ----------
    model : str
        The scattering model name.
    provenance : ProvenanceRecord, optional
    """

    model: ScatteringModelName = "atomic_number"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.model not in {"unit", "atomic_number", "xray_gaussian_proxy"}:
            raise ValueError(
                "ScatteringFactorTable.model must be 'unit', 'atomic_number', or "
                "'xray_gaussian_proxy'."
            )

    def scattering_factor(self, species: str, g_magnitude_inv_angstrom: float) -> float:
        """Atomic scattering factor for a species at a given ``|g|``.

        Parameters
        ----------
        species : str
            Element symbol.
        g_magnitude_inv_angstrom : float
            Scattering-vector magnitude, finite and non-negative.

        Returns
        -------
        float
            Depends on the table's ``model``: ``"unit"`` returns 1 (structure
            factors then count sites rather than scattering power);
            ``"atomic_number"`` returns ``Z``, the correct forward-scattering
            limit but with no angular falloff; the remaining model applies a
            smooth monotonic decay with ``|g|``.

        Notes
        -----
        These are deliberate proxies, not tabulated Cromer-Mann or Doyle-Turner
        factors. They give correct relative intensities among reflections of
        similar angle and correct systematic absences, but they are not
        quantitative across a wide angular range.
        """

        if not np.isfinite(g_magnitude_inv_angstrom) or g_magnitude_inv_angstrom < 0.0:
            raise ValueError("g_magnitude_inv_angstrom must be finite and non-negative.")
        if self.model == "unit":
            return 1.0
        z = float(atomic_number(species))
        if self.model == "atomic_number":
            return z
        # Smooth monotonic proxy. This is not a tabulated Cromer-Mann replacement.
        return float(z * np.exp(-0.02 * g_magnitude_inv_angstrom * g_magnitude_inv_angstrom))


@dataclass(frozen=True, slots=True)
class StructureFactor:
    """The complex amplitude ``F(hkl)`` scattered by one unit cell.

    Purpose
    -------
    ``F = sum_j f_j exp(2 pi i (h x_j + k y_j + l z_j))`` over the atomic
    basis. Its modulus squared sets reflection intensity, and its vanishing
    produces the systematic absences that identify a structure.

    Attributes
    ----------
    miller_indices : np.ndarray
        The reflection.
    value : complex
        The structure-factor amplitude.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Thermal (Debye-Waller) and anomalous-dispersion terms are not included.
    """

    miller_indices: np.ndarray
    value: complex
    amplitude: float
    phase_rad: float
    phase: Phase
    scattering_table: ScatteringFactorTable = field(default_factory=ScatteringFactorTable)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        if not np.isfinite(self.value.real) or not np.isfinite(self.value.imag):
            raise ValueError("StructureFactor.value must have finite real and imaginary parts.")
        if not np.isfinite(self.amplitude) or self.amplitude < 0.0:
            raise ValueError("StructureFactor.amplitude must be finite and non-negative.")
        if not np.isfinite(self.phase_rad):
            raise ValueError("StructureFactor.phase_rad must be finite.")

    @classmethod
    def from_phase_hkl(
        cls,
        phase: Phase,
        hkl: ArrayLike,
        *,
        scattering_table: ScatteringFactorTable | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> StructureFactor:
        """Structure factor ``F(hkl)`` of a phase from its atomic basis.

        Purpose
        -------
        The amplitude scattered by one unit cell,
        ``F = sum_j f_j exp(2 pi i (h x_j + k y_j + l z_j))``. Its modulus
        squared sets reflection intensity, and its vanishing is what produces
        systematic absences.

        Parameters
        ----------
        phase : Phase
            Must carry a unit cell with atomic sites for a meaningful result.
        hkl : ArrayLike
            The reflection indices.
        scattering_table : ScatteringFactorTable, optional
            Atomic scattering model; the default is used when omitted.
        provenance : ProvenanceRecord, optional

        Returns
        -------
        StructureFactor
            Carrying the complex amplitude.

        Notes
        -----
        A phase with no unit cell has no atomic basis to sum over; the factor is
        then taken as ``1 + 0j``, which means lattice-only reasoning applies and
        absences from the atomic basis cannot be detected. Thermal
        (Debye-Waller) and anomalous-dispersion terms are not included.
        """

        indices = as_int_array(hkl, shape=(3,))
        table = ScatteringFactorTable() if scattering_table is None else scattering_table
        if phase.unit_cell is None or not phase.unit_cell.sites:
            value = 1.0 + 0.0j
        else:
            reciprocal = phase.lattice.reciprocal_basis().matrix
            g_cart = reciprocal @ indices.astype(np.float64)
            g_magnitude = float(np.linalg.norm(g_cart))
            total = 0.0j
            for site in phase.unit_cell.sites:
                factor = table.scattering_factor(site.species, g_magnitude)
                occupancy = float(site.occupancy)
                b_iso = 0.0 if site.b_iso is None else float(site.b_iso)
                debye_waller = np.exp(
                    -(b_iso * g_magnitude * g_magnitude) / max(16.0 * np.pi * np.pi, 1e-12)
                )
                phase_argument = float(
                    np.dot(indices.astype(np.float64), site.fractional_coordinates)
                )
                phase_factor = np.exp(
                    2.0j * np.pi * phase_argument
                )
                total += occupancy * factor * debye_waller * phase_factor
            value = complex(total)
        return cls(
            miller_indices=indices,
            value=value,
            amplitude=float(abs(value)),
            phase_rad=float(np.angle(value)),
            phase=phase,
            scattering_table=table,
            provenance=provenance or table.provenance or phase.provenance,
        )


@dataclass(frozen=True, slots=True)
class ReflectionCondition:
    """The lattice-centring rule deciding which reflections are allowed.

    Purpose
    -------
    Centring absences are the most visible systematic absences in a pattern:
    a body-centred metal shows only reflections with ``h+k+l`` even. Applying
    them is what keeps a simulation from listing forbidden reflections as
    present.

    Limits
    ------
    Centring only. Glide-plane and screw-axis absences depend on the full
    space group rather than on its centring letter, and are not applied.

    Attributes
    ----------
    centering : str
        Lattice centring letter — ``P``, ``I``, ``F``, ``A``, ``B``, ``C``,
        or ``R``.
    provenance : ProvenanceRecord, optional
    """

    centering: str = "P"
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        centering = self.centering.strip().upper()
        if centering not in {"P", "I", "F", "A", "B", "C", "R"}:
            raise ValueError("ReflectionCondition.centering must be P, I, F, A, B, C, or R.")
        object.__setattr__(self, "centering", centering)

    @classmethod
    def from_phase(cls, phase: Phase) -> ReflectionCondition:
        """The lattice-centring reflection condition implied by a phase.

        Reads the leading letter of the space-group symbol — ``P``, ``I``, ``F``,
        ``A``, ``B``, ``C``, ``R`` — which is what determines the centring
        absences. A phase with no space group is treated as primitive, meaning
        no centring absences are applied; that is a real limitation on the
        result, not a safe default, since simulating a body-centred metal
        without it lists forbidden reflections as present.
        """

        symbol = phase.space_group_symbol
        if symbol is None:
            symbol = phase.space_group.symbol if phase.space_group else "P"
        centering = symbol.strip()[:1].upper() if symbol.strip() else "P"
        return cls(centering=centering, provenance=phase.provenance)

    def is_allowed(self, hkl: ArrayLike) -> bool:
        """Whether a reflection survives this lattice's centring condition.

        Applies the standard integral conditions: ``h+k+l`` even for ``I``;
        ``h, k, l`` all even or all odd for ``F``; ``k+l``, ``h+l``, ``h+k`` even
        for ``A``, ``B``, ``C``; ``-h+k+l`` divisible by three for ``R``; and no
        condition for ``P``.

        These are centring (lattice) absences only. Glide-plane and screw-axis
        absences, which depend on the full space group rather than its centring
        letter, are not applied here.
        """

        h, k, ell = (int(value) for value in as_int_array(hkl, shape=(3,)))
        if self.centering == "P":
            return True
        if self.centering == "I":
            return (h + k + ell) % 2 == 0
        if self.centering == "F":
            return h % 2 == k % 2 == ell % 2
        if self.centering == "A":
            return (k + ell) % 2 == 0
        if self.centering == "B":
            return (h + ell) % 2 == 0
        if self.centering == "C":
            return (h + k) % 2 == 0
        return (-h + k + ell) % 3 == 0


@dataclass(frozen=True, slots=True)
class DiffractionIntensityModel:
    """How structure factors are turned into observable powder intensities.

    Purpose
    -------
    Combines the structure factor, the reflection multiplicity, and the
    angle-dependent Lorentz-polarization factor. It also holds the reflection
    condition, so absences are applied consistently with the intensity
    calculation rather than separately.

    Limits
    ------
    Kinematic and relative. Absorption, preferred orientation, extinction,
    and thermal factors are not modelled, so these intensities are not
    quantitative phase-fraction inputs.

    Attributes
    ----------
    scattering_table : ScatteringFactorTable, optional
    reflection_condition : ReflectionCondition, optional
        Derived from the phase when omitted.
    provenance : ProvenanceRecord, optional
    """

    scattering_table: ScatteringFactorTable = field(default_factory=ScatteringFactorTable)
    correction: IntensityCorrectionName = "lorentz_polarization"
    reflection_condition: ReflectionCondition | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if self.correction not in {"none", "lorentz_polarization"}:
            raise ValueError("correction must be either 'none' or 'lorentz_polarization'.")

    def intensity(
        self,
        phase: Phase,
        hkl: ArrayLike,
        *,
        two_theta_rad: float,
        multiplicity: int = 1,
        radiation: RadiationSpec | None = None,
    ) -> float:
        """Powder intensity of one reflection under this model.

        Purpose
        -------
        Combine the structure factor, multiplicity, and the angle-dependent
        Lorentz-polarization factor into the relative intensity a powder pattern
        shows.

        Parameters
        ----------
        phase : Phase
        hkl : ArrayLike
            The reflection.
        two_theta_rad : float
            Scattering angle, strictly between 0 and ``pi``.
        multiplicity : int
            Number of symmetry-equivalent reflections contributing at this
            angle; a positive integer.
        radiation : RadiationSpec, optional
            Accepted for interface symmetry and currently unused, because the
            scattering models here are not wavelength-dependent.

        Returns
        -------
        float
            Relative intensity; zero for a reflection forbidden by the centring
            condition.

        Notes
        -----
        Kinematic and relative. Absorption, preferred orientation, extinction,
        and thermal factors are not modelled, so these intensities are not
        quantitative phase-fraction inputs.
        """

        del radiation
        if not np.isfinite(two_theta_rad) or two_theta_rad <= 0.0 or two_theta_rad >= np.pi:
            raise ValueError("two_theta_rad must be finite and satisfy 0 < two_theta_rad < pi.")
        if int(multiplicity) != multiplicity or multiplicity <= 0:
            raise ValueError("multiplicity must be a positive integer.")
        condition = self.reflection_condition or ReflectionCondition.from_phase(phase)
        if not condition.is_allowed(hkl):
            return 0.0
        structure_factor = StructureFactor.from_phase_hkl(
            phase,
            hkl,
            scattering_table=self.scattering_table,
            provenance=self.provenance,
        )
        correction = (
            1.0
            if self.correction == "none"
            else lorentz_polarization_factor(two_theta_rad)
        )
        return float(max(0.0, multiplicity * structure_factor.amplitude**2 * correction))


def phase_centering_is_declared(phase: Phase) -> bool:
    """Whether a phase actually states its lattice centering, or `P` was assumed.

    Purpose: makes a silent failure mode visible.
    `ReflectionCondition.from_phase` reads the centering from the first letter
    of the space-group symbol and falls back to primitive when the phase carries
    no symbol at all. A body-centred phase supplied without that metadata is
    therefore simulated as primitive, and its pattern shows reflections that the
    real structure forbids — with nothing in the output to say so.

    When to use: before trusting a simulated reflection list, and in any report
    or manifest that records how a pattern was produced. Simulation surfaces in
    this package use it to label the applied centering as *declared* or
    *assumed*.

    Inputs: the phase.

    Output: ``True`` when the phase carries a non-empty space-group symbol
    (directly or through its `SpaceGroupSpec`), ``False`` when the primitive
    default was assumed.

    See also
    --------
    `ReflectionCondition.from_phase` : the centering derivation itself.
    """

    symbol = phase.space_group_symbol
    if symbol is None and phase.space_group is not None:
        symbol = phase.space_group.symbol
    return bool(symbol and symbol.strip())


__all__ = [
    "DiffractionIntensityModel",
    "ReflectionCondition",
    "ScatteringFactorTable",
    "StructureFactor",
    "lorentz_polarization_factor",
    "phase_centering_is_declared",
]
