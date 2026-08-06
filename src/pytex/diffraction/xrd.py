from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from pytex.core._arrays import as_float_array, as_int_array
from pytex.core._chemistry import atomic_number
from pytex.core.lattice import Phase
from pytex.core.miller import MillerPlane
from pytex.core.provenance import ProvenanceRecord
from pytex.diffraction.preferred_orientation import PreferredOrientationModel

_TWO_THETA_TOLERANCE_DEG = 1e-8


def _gaussian_profile(x: np.ndarray, center: float, sigma: float) -> np.ndarray:
    exponent = -0.5 * ((x - center) / sigma) ** 2
    profile = np.exp(exponent) / (sigma * np.sqrt(2.0 * np.pi))
    profile = np.ascontiguousarray(profile, dtype=np.float64)
    profile.setflags(write=False)
    return profile


def _lorentzian_profile(x: np.ndarray, center: float, gamma: float) -> np.ndarray:
    """Area-normalized Lorentzian with half-width-at-half-maximum ``gamma``."""

    profile = gamma / (np.pi * ((x - center) ** 2 + gamma * gamma))
    profile = np.ascontiguousarray(profile, dtype=np.float64)
    profile.setflags(write=False)
    return profile


def _pseudo_voigt_profile(x: np.ndarray, center: float, fwhm: float, eta: float) -> np.ndarray:
    """Area-normalized pseudo-Voigt: ``eta`` Lorentzian + ``1 - eta`` Gaussian.

    ``eta = 0`` recovers the pure Gaussian, ``eta = 1`` the pure Lorentzian,
    both sharing the same full width at half maximum.
    """

    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gaussian = _gaussian_profile(x, center, sigma)
    lorentzian = _lorentzian_profile(x, center, 0.5 * fwhm)
    profile = eta * lorentzian + (1.0 - eta) * gaussian
    profile = np.ascontiguousarray(profile, dtype=np.float64)
    profile.setflags(write=False)
    return profile


def _caglioti_fwhm_deg(
    two_theta_deg: float, caglioti_uvw: tuple[float, float, float]
) -> float:
    """Caglioti instrumental width: FWHM^2 = U tan^2(theta) + V tan(theta) + W (deg^2)."""

    theta = np.deg2rad(0.5 * two_theta_deg)
    tangent = float(np.tan(theta))
    u, v, w = caglioti_uvw
    fwhm_squared = u * tangent * tangent + v * tangent + w
    if fwhm_squared <= 0.0:
        raise ValueError(
            "caglioti_uvw produced a non-positive squared FWHM at "
            f"two-theta = {two_theta_deg:.3f} deg."
        )
    return float(np.sqrt(fwhm_squared))


def _lorentz_polarization(two_theta_rad: float) -> float:
    theta = 0.5 * two_theta_rad
    sin_theta = max(float(np.sin(theta)), 1e-8)
    cos_theta = max(float(np.cos(theta)), 1e-8)
    cos_two_theta = float(np.cos(two_theta_rad))
    return float((1.0 + cos_two_theta * cos_two_theta) / (sin_theta * sin_theta * cos_theta))


def _structure_factors_xray(
    phase: Phase, hkls: np.ndarray, *, tabulated: bool = False
) -> np.ndarray:
    """X-ray structure factors for a batch of ``hkls`` (shape (n, 3)) -> (n,) complex.

    Vectorised over both reflections and unit-cell sites. With ``tabulated``
    the atomic-number proxy ``f = Z`` is replaced by the tabulated
    angle-dependent form factors ``f(s)`` from
    `pytex.diffraction.scattering`.
    """

    hkls_float = np.atleast_2d(np.asarray(hkls, dtype=np.float64))
    if phase.unit_cell is None or not phase.unit_cell.sites:
        return np.ones(hkls_float.shape[0], dtype=np.complex128)
    reciprocal = phase.lattice.reciprocal_basis().matrix
    g_cart = hkls_float @ reciprocal.T  # (n, 3): reciprocal @ hkl per row
    g_sq = np.sum(g_cart * g_cart, axis=1)  # (n,)
    sites = phase.unit_cell.sites
    if tabulated:
        from pytex.diffraction.scattering import xray_form_factor_matrix

        s_values = 0.5 * np.sqrt(g_sq)  # s = sin(theta)/lambda = |g| / 2
        z = xray_form_factor_matrix(
            tuple(site.species for site in sites), s_values
        )  # (n, s)
    else:
        z = np.array([float(atomic_number(site.species)) for site in sites], dtype=np.float64)
    occupancy = np.array([float(site.occupancy) for site in sites], dtype=np.float64)
    b_iso = np.array(
        [0.0 if site.b_iso is None else float(site.b_iso) for site in sites], dtype=np.float64
    )
    fractional = np.array(
        [np.asarray(site.fractional_coordinates, dtype=np.float64) for site in sites]
    )  # (s, 3)
    debye_waller = np.exp(
        -(b_iso[None, :] * g_sq[:, None]) / max(16.0 * np.pi * np.pi, 1e-12)
    )  # (n, s)
    phase_factor = np.exp(2.0j * np.pi * (hkls_float @ fractional.T))  # (n, s)
    form_factors = z if z.ndim == 2 else np.broadcast_to(z[None, :], phase_factor.shape)
    contributions = occupancy[None, :] * form_factors * debye_waller * phase_factor
    return np.asarray(np.sum(contributions, axis=1), dtype=np.complex128)


def _structure_factor_xray(phase: Phase, hkl: np.ndarray) -> complex:
    return complex(_structure_factors_xray(phase, np.asarray(hkl)[None, :])[0])


def _equivalent_hkls(phase: Phase, hkl: np.ndarray) -> tuple[np.ndarray, ...]:
    reciprocal_basis = phase.lattice.reciprocal_basis().matrix
    inverse_basis = np.linalg.inv(reciprocal_basis)
    g_cart = reciprocal_basis @ hkl.astype(np.float64)
    equivalents: set[tuple[int, int, int]] = set()
    for operator in phase.symmetry.operators:
        transformed = operator @ g_cart
        transformed_hkl = inverse_basis @ transformed
        rounded = np.rint(transformed_hkl).astype(np.int64)
        if np.allclose(transformed_hkl, rounded, atol=1e-6):
            equivalents.add((int(rounded[0]), int(rounded[1]), int(rounded[2])))
            negatives = tuple([-int(value) for value in rounded])
            equivalents.add((negatives[0], negatives[1], negatives[2]))
    return tuple(np.array(item, dtype=np.int64) for item in sorted(equivalents))


def _reflection_multiplicity(phase: Phase, hkl: np.ndarray) -> int:
    return len(_equivalent_hkls(phase, hkl))


def _reflection_family_key(phase: Phase, hkl: np.ndarray) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (int(equivalent[0]), int(equivalent[1]), int(equivalent[2]))
        for equivalent in _equivalent_hkls(phase, hkl)
    )


def _reflection_family_representative(
    family_key: tuple[tuple[int, int, int], ...],
) -> np.ndarray:
    representative = min(
        family_key,
        key=lambda indices: (
            sum(value < 0 for value in indices),
            tuple(abs(value) for value in indices),
            indices,
        ),
    )
    return np.asarray(representative, dtype=np.int64)


def _enumerate_hkls(max_index: int) -> np.ndarray:
    values = range(-max_index, max_index + 1)
    hkls = [
        (h, k, ell)
        for h in values
        for k in values
        for ell in values
        if not (h == 0 and k == 0 and ell == 0)
    ]
    array = np.asarray(hkls, dtype=np.int64)
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class RadiationSpec:
    """A diffraction radiation line (or K-alpha doublet).

    ``wavelength_angstrom`` is the primary (K-alpha1) line. When
    ``kalpha2_wavelength_angstrom`` is set, pattern generation superposes the
    K-alpha2 contribution weighted by ``kalpha2_relative_intensity`` (the
    conventional 0.5). ``anode`` records the X-ray tube target;
    ``kind`` distinguishes X-ray from neutron radiation.
    """

    name: str
    wavelength_angstrom: float
    kalpha2_wavelength_angstrom: float | None = None
    kalpha2_relative_intensity: float = 0.5
    anode: str | None = None
    kind: str = "xray"

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("RadiationSpec.name must be non-empty.")
        if self.wavelength_angstrom <= 0.0:
            raise ValueError("RadiationSpec.wavelength_angstrom must be strictly positive.")
        if self.kalpha2_wavelength_angstrom is not None:
            if self.kalpha2_wavelength_angstrom <= 0.0:
                raise ValueError(
                    "RadiationSpec.kalpha2_wavelength_angstrom must be strictly positive."
                )
            if not 0.0 < self.kalpha2_relative_intensity <= 1.0:
                raise ValueError(
                    "RadiationSpec.kalpha2_relative_intensity must lie in (0, 1]."
                )
        if self.kind not in {"xray", "neutron"}:
            raise ValueError("RadiationSpec.kind must be 'xray' or 'neutron'.")

    @classmethod
    def cu_ka(cls) -> RadiationSpec:
        """Cu K-alpha radiation, K-alpha wavelength.

        The most common laboratory X-ray source. Uses a single averaged
        wavelength; for doublet-resolved peak shapes use :meth:`cu_ka_doublet`.
        """

        return cls(name="Cu Ka", wavelength_angstrom=1.5406, anode="Cu")

    @classmethod
    def mo_ka(cls) -> RadiationSpec:
        """Mo K-alpha radiation, K-alpha wavelength.

        Shorter wavelength than Cu, so more reflections fall inside a given
        ``2*theta`` range.
        """

        return cls(name="Mo Ka", wavelength_angstrom=0.71073, anode="Mo")

    # Bearden (1967) K-alpha1 / K-alpha2 wavelengths for the common anodes.
    @classmethod
    def cu_ka_doublet(cls) -> RadiationSpec:
        """Cu radiation, K-alpha1/K-alpha2 doublet wavelengths (Bearden 1967).

        Models the K-alpha1/K-alpha2 pair explicitly, which is what produces the
        characteristic peak splitting at high ``2*theta``.
        """

        return cls(
            name="Cu Ka1/Ka2",
            wavelength_angstrom=1.540562,
            kalpha2_wavelength_angstrom=1.544390,
            anode="Cu",
        )

    @classmethod
    def mo_ka_doublet(cls) -> RadiationSpec:
        """Mo radiation, K-alpha1/K-alpha2 doublet wavelengths (Bearden 1967).
        """

        return cls(
            name="Mo Ka1/Ka2",
            wavelength_angstrom=0.709300,
            kalpha2_wavelength_angstrom=0.713590,
            anode="Mo",
        )

    @classmethod
    def co_ka(cls) -> RadiationSpec:
        """Co radiation, K-alpha1/K-alpha2 doublet wavelengths (Bearden 1967).

        Preferred for iron-bearing samples, where Cu radiation excites strong Fe
        fluorescence and raises the background.
        """

        return cls(
            name="Co Ka1/Ka2",
            wavelength_angstrom=1.788965,
            kalpha2_wavelength_angstrom=1.792850,
            anode="Co",
        )

    @classmethod
    def cr_ka(cls) -> RadiationSpec:
        """Cr radiation, K-alpha1/K-alpha2 doublet wavelengths (Bearden 1967).

        Long wavelength; conventional for residual-stress measurement, where the
        large ``2*theta`` of accessible reflections improves strain sensitivity.
        """

        return cls(
            name="Cr Ka1/Ka2",
            wavelength_angstrom=2.289700,
            kalpha2_wavelength_angstrom=2.293606,
            anode="Cr",
        )

    @classmethod
    def fe_ka(cls) -> RadiationSpec:
        """Fe radiation, K-alpha1/K-alpha2 doublet wavelengths (Bearden 1967).
        """

        return cls(
            name="Fe Ka1/Ka2",
            wavelength_angstrom=1.936042,
            kalpha2_wavelength_angstrom=1.939980,
            anode="Fe",
        )

    @classmethod
    def neutron(cls, wavelength_angstrom: float, *, name: str = "neutron") -> RadiationSpec:
        """A neutron radiation specification at a given wavelength.

        Marked as neutron rather than X-ray radiation. Note that neutron
        scattering lengths differ fundamentally from X-ray form factors — they do
        not scale with atomic number and can be negative — so intensity models
        built on X-ray scattering factors do not transfer to neutron data.
        """

        return cls(name=name, wavelength_angstrom=wavelength_angstrom, kind="neutron")


@dataclass(frozen=True, slots=True)
class PowderReflection:
    """One reflection family in a powder diffraction pattern.

    Attributes
    ----------
    miller_indices : np.ndarray
        The family representative ``(hkl)``.
    two_theta_deg : float
        Bragg angle for the pattern's radiation.
    d_spacing_angstrom : float
    multiplicity : int
        Number of symmetry-equivalent reflections contributing at this angle.
        This is why a high-multiplicity family can outshine one with a larger
        structure factor.
    intensity : float
        Relative intensity, kinematic.
    """

    miller_indices: np.ndarray
    d_spacing_angstrom: float
    two_theta_deg: float
    intensity: float
    structure_factor_amplitude: float
    multiplicity: int
    structure_factor_real: float | None = None
    structure_factor_imag: float | None = None
    lorentz_polarization_factor: float | None = None
    intensity_model: str = "xray_atomic_number"

    def __post_init__(self) -> None:
        object.__setattr__(self, "miller_indices", as_int_array(self.miller_indices, shape=(3,)))
        if self.d_spacing_angstrom <= 0.0:
            raise ValueError("PowderReflection.d_spacing_angstrom must be strictly positive.")
        if self.two_theta_deg < 0.0:
            raise ValueError("PowderReflection.two_theta_deg must be non-negative.")
        if not np.isfinite(self.intensity) or self.intensity < 0.0:
            raise ValueError("PowderReflection.intensity must be finite and non-negative.")
        if (
            not np.isfinite(self.structure_factor_amplitude)
            or self.structure_factor_amplitude < 0.0
        ):
            raise ValueError(
                "PowderReflection.structure_factor_amplitude must be finite and non-negative."
            )
        if self.multiplicity <= 0:
            raise ValueError("PowderReflection.multiplicity must be strictly positive.")
        for name, value in (
            ("structure_factor_real", self.structure_factor_real),
            ("structure_factor_imag", self.structure_factor_imag),
            ("lorentz_polarization_factor", self.lorentz_polarization_factor),
        ):
            if value is not None and not np.isfinite(value):
                raise ValueError(f"PowderReflection.{name} must be finite when provided.")
        if not self.intensity_model.strip():
            raise ValueError("PowderReflection.intensity_model must be non-empty.")


@dataclass(frozen=True, slots=True)
class PowderPattern:
    """A simulated continuous powder diffraction pattern.

    Purpose
    -------
    Reflection positions convolved with a peak profile onto a ``2*theta``
    grid — the form a measured diffractogram takes, and the form a comparison
    against measurement needs.

    Limits
    ------
    Kinematic intensities with no background model, no absorption, no
    preferred-orientation correction, and no size or strain broadening beyond
    the configured profile. Suitable for phase identification and peak
    indexing, not for quantitative phase analysis.

    Attributes
    ----------
    two_theta_deg : np.ndarray
        The angular grid.
    intensity : np.ndarray
        Profile-convolved intensity on that grid.
    reflections : tuple of PowderReflection
        The underlying reflection list, retained so peaks stay attributable
        to their indices.
    radiation : RadiationSpec
    phase : Phase
    provenance : ProvenanceRecord, optional
    """

    phase: Phase
    radiation: RadiationSpec
    reflections: tuple[PowderReflection, ...]
    two_theta_grid_deg: np.ndarray
    intensity_grid: np.ndarray
    broadening_fwhm_deg: float | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "reflections", tuple(self.reflections))
        object.__setattr__(
            self,
            "two_theta_grid_deg",
            as_float_array(self.two_theta_grid_deg, shape=(None,)),
        )
        object.__setattr__(
            self,
            "intensity_grid",
            as_float_array(self.intensity_grid, shape=(None,)),
        )
        if self.two_theta_grid_deg.shape != self.intensity_grid.shape:
            raise ValueError("PowderPattern grid arrays must have the same shape.")
        if np.any(~np.isfinite(self.intensity_grid)) or np.any(self.intensity_grid < 0.0):
            raise ValueError("PowderPattern intensity_grid must be finite and non-negative.")


def generate_powder_reflections(
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    two_theta_range_deg: tuple[float, float] = (5.0, 120.0),
    max_index: int = 6,
    intensity_model: Literal["xray_atomic_number", "xray_tabulated", "unit"] = (
        "xray_atomic_number"
    ),
) -> tuple[PowderReflection, ...]:
    """Enumerate the powder reflections of a phase in a ``2*theta`` window.

    Purpose
    -------
    The reflection list behind a simulated powder pattern: which ``(hkl)``
    diffract, at what angle, with what multiplicity and relative intensity.

    Method
    ------
    All reflections up to ``max_index`` are enumerated and filtered
    vectorially by the Bragg condition and the requested angular window,
    then grouped into symmetry families so each family contributes one
    entry carrying its multiplicity.

    Parameters
    ----------
    phase : Phase
        Lattice, symmetry, and — for structure factors — the atomic basis.
    radiation : RadiationSpec, optional
        Defaults to Cu K-alpha.
    two_theta_range_deg : tuple of float
        Inclusive angular window; must satisfy ``0 <= min < max <= 180``.
    max_index : int
        Largest absolute Miller index enumerated. Too small silently omits
        high-angle reflections, so raise it if the pattern looks truncated.
    intensity_model : str
        ``"xray_atomic_number"`` (default), ``"xray_tabulated"``, or
        ``"unit"``.

    Returns
    -------
    tuple of PowderReflection
        One entry per reflection family, with ``2*theta``, d-spacing,
        multiplicity, and relative intensity.

    See Also
    --------
    generate_xrd_pattern : Convolves these reflections with a peak profile
        to produce a continuous pattern.
    """

    if max_index <= 0:
        raise ValueError("max_index must be strictly positive.")
    radiation_spec = RadiationSpec.cu_ka() if radiation is None else radiation
    min_two_theta, max_two_theta = (float(two_theta_range_deg[0]), float(two_theta_range_deg[1]))
    if not 0.0 <= min_two_theta < max_two_theta <= 180.0:
        raise ValueError("two_theta_range_deg must satisfy 0 <= min < max <= 180.")

    # Vectorised geometric filtering over all candidate reflections at once.
    all_hkls = _enumerate_hkls(max_index)
    g_cartesian = all_hkls.astype(np.float64) @ phase.lattice.reciprocal_basis().matrix.T
    g_magnitudes = np.linalg.norm(g_cartesian, axis=1)
    nonzero = ~np.isclose(g_magnitudes, 0.0)
    safe_magnitudes = np.where(nonzero, g_magnitudes, 1.0)
    d_spacings = 1.0 / safe_magnitudes
    arguments = radiation_spec.wavelength_angstrom / (2.0 * d_spacings)
    two_theta_all = np.rad2deg(2.0 * np.arcsin(np.clip(arguments, -1.0, 1.0)))
    keep = (
        nonzero
        & (arguments > 0.0)
        & (arguments <= 1.0)
        & (two_theta_all >= min_two_theta - _TWO_THETA_TOLERANCE_DEG)
        & (two_theta_all <= max_two_theta + _TWO_THETA_TOLERANCE_DEG)
    )
    candidate_hkls = all_hkls[keep]
    visited_families: set[tuple[tuple[int, int, int], ...]] = set()
    representative_hkls: list[np.ndarray] = []
    for hkl in candidate_hkls:
        family_key = _reflection_family_key(phase, hkl)
        if family_key in visited_families:
            continue
        visited_families.add(family_key)
        representative_hkls.append(_reflection_family_representative(family_key))
    surviving_hkls = np.asarray(representative_hkls, dtype=np.int64).reshape((-1, 3))
    surviving_g = surviving_hkls.astype(np.float64) @ phase.lattice.reciprocal_basis().matrix.T
    surviving_d = 1.0 / np.linalg.norm(surviving_g, axis=1)
    surviving_arguments = radiation_spec.wavelength_angstrom / (2.0 * surviving_d)
    surviving_two_theta = np.rad2deg(2.0 * np.arcsin(surviving_arguments))
    if intensity_model == "unit":
        structure_factors = np.ones(surviving_hkls.shape[0], dtype=np.complex128)
    else:
        structure_factors = _structure_factors_xray(
            phase, surviving_hkls, tabulated=intensity_model == "xray_tabulated"
        )

    reflections: list[PowderReflection] = []
    for index in range(surviving_hkls.shape[0]):
        hkl = surviving_hkls[index]
        structure_factor = complex(structure_factors[index])
        amplitude = float(abs(structure_factor))
        multiplicity = _reflection_multiplicity(phase, hkl)
        two_theta_deg = float(surviving_two_theta[index])
        lorentz_polarization = _lorentz_polarization(np.deg2rad(two_theta_deg))
        intensity = multiplicity * amplitude * amplitude * lorentz_polarization
        if intensity <= 1e-14:
            continue
        reflections.append(
            PowderReflection(
                miller_indices=hkl,
                d_spacing_angstrom=float(surviving_d[index]),
                two_theta_deg=two_theta_deg,
                intensity=intensity,
                structure_factor_amplitude=amplitude,
                multiplicity=multiplicity,
                structure_factor_real=float(structure_factor.real),
                structure_factor_imag=float(structure_factor.imag),
                lorentz_polarization_factor=lorentz_polarization,
                intensity_model=intensity_model,
            )
        )
    # Round the primary key so reflections at the same 2-theta (symmetry
    # equivalents whose 2-theta differs only by floating-point noise) tie and
    # order deterministically by Miller indices, independent of FP rounding.
    reflections.sort(
        key=lambda reflection: (
            round(reflection.two_theta_deg, 9),
            tuple(int(value) for value in reflection.miller_indices),
        )
    )
    return tuple(reflections)


def _accumulate_reflection_profiles(
    grid: np.ndarray,
    reflections: tuple[PowderReflection, ...],
    *,
    weight: float,
    broadening_fwhm_deg: float | None,
    profile: str,
    pseudo_voigt_eta: float,
    caglioti_uvw: tuple[float, float, float] | None,
) -> np.ndarray:
    """Sum weighted peak profiles for one radiation line onto the 2-theta grid."""

    intensity_grid = np.zeros_like(grid)
    if broadening_fwhm_deg is None and caglioti_uvw is None:
        for reflection in reflections:
            index = int(np.argmin(np.abs(grid - reflection.two_theta_deg)))
            intensity_grid[index] += weight * reflection.intensity
        return intensity_grid
    for reflection in reflections:
        fwhm = (
            _caglioti_fwhm_deg(reflection.two_theta_deg, caglioti_uvw)
            if caglioti_uvw is not None
            else float(broadening_fwhm_deg)  # type: ignore[arg-type]
        )
        if profile == "pseudo_voigt":
            peak = _pseudo_voigt_profile(
                grid, reflection.two_theta_deg, fwhm, pseudo_voigt_eta
            )
        else:
            sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            peak = _gaussian_profile(grid, reflection.two_theta_deg, sigma)
        intensity_grid += weight * reflection.intensity * peak
    return intensity_grid


def apply_preferred_orientation(
    reflections: Sequence[PowderReflection],
    model: PreferredOrientationModel,
    *,
    phase: Phase,
) -> tuple[PowderReflection, ...]:
    """Scale powder reflection intensities by a preferred-orientation model.

    Purpose
    -------
    Apply a texture correction to an already-computed reflection list, leaving
    every other field — indices, spacing, angle, multiplicity, structure
    factor — unchanged, so the corrected list stays comparable with the
    uncorrected one reflection by reflection.

    Parameters
    ----------
    reflections : sequence of PowderReflection
        The random-powder reflection list.
    model : PreferredOrientationModel
        Any object satisfying the protocol, such as
        :class:`~pytex.diffraction.preferred_orientation.MarchDollaseModel` or
        :class:`~pytex.diffraction.preferred_orientation.ODFPreferredOrientationModel`.
    phase : Phase
        The phase the reflections belong to; needed to build the typed planes
        the model consumes.

    Returns
    -------
    tuple of PowderReflection
        Corrected reflections in the input order.

    Notes
    -----
    Powder intensities are relative, so the corrected list is deliberately not
    renormalized: the *ratios* between reflections carry the texture
    information, and rescaling to a fixed maximum would discard part of the
    effect being modelled.
    """

    if not reflections:
        return ()
    planes = [
        MillerPlane(indices=reflection.miller_indices, phase=phase)
        for reflection in reflections
    ]
    factors = np.asarray(model.factors(planes), dtype=np.float64)
    if factors.shape != (len(reflections),):
        raise ValueError(
            "A preferred-orientation model must return exactly one factor per reflection."
        )
    if np.any(~np.isfinite(factors)) or np.any(factors < 0.0):
        raise ValueError("Preferred-orientation factors must be finite and non-negative.")
    return tuple(
        PowderReflection(
            miller_indices=reflection.miller_indices,
            d_spacing_angstrom=reflection.d_spacing_angstrom,
            two_theta_deg=reflection.two_theta_deg,
            intensity=float(reflection.intensity * factor),
            structure_factor_amplitude=reflection.structure_factor_amplitude,
            multiplicity=reflection.multiplicity,
            structure_factor_real=reflection.structure_factor_real,
            structure_factor_imag=reflection.structure_factor_imag,
            lorentz_polarization_factor=reflection.lorentz_polarization_factor,
            intensity_model=reflection.intensity_model,
        )
        for reflection, factor in zip(reflections, factors, strict=True)
    )


def generate_xrd_pattern(
    phase: Phase,
    *,
    radiation: RadiationSpec | None = None,
    two_theta_range_deg: tuple[float, float] = (5.0, 120.0),
    resolution_deg: float = 0.02,
    max_index: int = 6,
    intensity_model: Literal["xray_atomic_number", "xray_tabulated", "unit"] = (
        "xray_atomic_number"
    ),
    broadening_fwhm_deg: float | None = 0.15,
    profile: Literal["gaussian", "pseudo_voigt"] = "gaussian",
    pseudo_voigt_eta: float = 0.5,
    caglioti_uvw: tuple[float, float, float] | None = None,
    preferred_orientation: PreferredOrientationModel | None = None,
    provenance: ProvenanceRecord | None = None,
) -> PowderPattern:
    """Simulate a powder XRD pattern for ``phase``.

    Peak shapes: ``profile="gaussian"`` (default) or ``"pseudo_voigt"`` with
    mixing ``pseudo_voigt_eta``. Angular widths come from the constant
    ``broadening_fwhm_deg`` or, when ``caglioti_uvw`` is given, from the
    Caglioti relation ``FWHM^2 = U tan^2(theta) + V tan(theta) + W``. When the
    radiation carries a K-alpha2 line, its pattern is superposed at
    ``kalpha2_relative_intensity``; the returned ``reflections`` remain the
    K-alpha1 list.

    Preferred orientation: pass a model such as
    :class:`~pytex.diffraction.preferred_orientation.MarchDollaseModel` or
    :class:`~pytex.diffraction.preferred_orientation.ODFPreferredOrientationModel`
    as ``preferred_orientation`` to scale the reflection intensities by texture.
    The returned ``reflections`` are then the *corrected* ones, so the pattern
    and its reflection list cannot disagree. Omitting it assumes a randomly
    oriented powder, which is an assumption about the specimen and not a neutral
    default.
    """

    if resolution_deg <= 0.0:
        raise ValueError("resolution_deg must be strictly positive.")
    if broadening_fwhm_deg is not None and broadening_fwhm_deg <= 0.0:
        raise ValueError("broadening_fwhm_deg must be strictly positive when provided.")
    if profile not in {"gaussian", "pseudo_voigt"}:
        raise ValueError("profile must be 'gaussian' or 'pseudo_voigt'.")
    if not 0.0 <= pseudo_voigt_eta <= 1.0:
        raise ValueError("pseudo_voigt_eta must lie in [0, 1].")
    if caglioti_uvw is not None and len(caglioti_uvw) != 3:
        raise ValueError("caglioti_uvw must provide (U, V, W).")
    radiation_spec = RadiationSpec.cu_ka() if radiation is None else radiation
    reflections = generate_powder_reflections(
        phase,
        radiation=radiation_spec,
        two_theta_range_deg=two_theta_range_deg,
        max_index=max_index,
        intensity_model=intensity_model,
    )
    if preferred_orientation is not None:
        reflections = apply_preferred_orientation(
            reflections, preferred_orientation, phase=phase
        )
    min_two_theta, max_two_theta = map(float, two_theta_range_deg)
    grid = np.arange(min_two_theta, max_two_theta + 0.5 * resolution_deg, resolution_deg)
    intensity_grid = _accumulate_reflection_profiles(
        grid,
        reflections,
        weight=1.0,
        broadening_fwhm_deg=broadening_fwhm_deg,
        profile=profile,
        pseudo_voigt_eta=pseudo_voigt_eta,
        caglioti_uvw=caglioti_uvw,
    )
    if radiation_spec.kalpha2_wavelength_angstrom is not None:
        kalpha2_spec = RadiationSpec(
            name=f"{radiation_spec.name} (Ka2)",
            wavelength_angstrom=radiation_spec.kalpha2_wavelength_angstrom,
            anode=radiation_spec.anode,
            kind=radiation_spec.kind,
        )
        kalpha2_reflections = generate_powder_reflections(
            phase,
            radiation=kalpha2_spec,
            two_theta_range_deg=two_theta_range_deg,
            max_index=max_index,
            intensity_model=intensity_model,
        )
        if preferred_orientation is not None:
            # The K-alpha2 line diffracts from the same textured specimen, so it
            # carries the same correction; omitting it here would make the
            # doublet ratio depend on the texture.
            kalpha2_reflections = apply_preferred_orientation(
                kalpha2_reflections, preferred_orientation, phase=phase
            )
        intensity_grid += _accumulate_reflection_profiles(
            grid,
            kalpha2_reflections,
            weight=radiation_spec.kalpha2_relative_intensity,
            broadening_fwhm_deg=broadening_fwhm_deg,
            profile=profile,
            pseudo_voigt_eta=pseudo_voigt_eta,
            caglioti_uvw=caglioti_uvw,
        )
    if (broadening_fwhm_deg is not None or caglioti_uvw is not None) and intensity_grid.size > 0:
        intensity_grid *= resolution_deg
    if np.max(intensity_grid) > 0.0:
        intensity_grid /= float(np.max(intensity_grid))
    intensity_grid = np.ascontiguousarray(intensity_grid, dtype=np.float64)
    intensity_grid.setflags(write=False)
    grid = np.ascontiguousarray(grid, dtype=np.float64)
    grid.setflags(write=False)
    return PowderPattern(
        phase=phase,
        radiation=radiation_spec,
        reflections=reflections,
        two_theta_grid_deg=grid,
        intensity_grid=intensity_grid,
        broadening_fwhm_deg=broadening_fwhm_deg,
        provenance=provenance,
    )
