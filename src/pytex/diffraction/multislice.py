"""Multislice simulation of high-resolution TEM: potentials, exit waves and image series.

This module is PyTex's own implementation of the conventional Fourier-space
multislice algorithm, written to reproduce the algorithm of abTEM (Madsen &
Susi, 2021) step for step so that the two can be compared number for number,
without needing abTEM installed. It turns an `AtomicSnapshot` into the
electron exit wave at any set of depths, and the exit wave into HRTEM images
through the objective lens of `MicroscopeAberrations` - singly, as a focal
series, as a thickness series, or as a defocus-thickness map.

The pipeline
------------

1. **Projected potential.** The snapshot is cut into slices of thickness
   :math:`\\Delta z` along the beam (+z). Each atom is projected entirely into
   the slice holding its centre ("infinite projection"), and the slice's
   projected potential is assembled in Fourier space from the independent-atom
   electron scattering factor :math:`f_e(g)`,

   .. math::

      \\tilde v_n(\\mathbf g) = \\frac{h^2}{2\\pi m_0 e}\\,\\frac{1}{A}
          \\sum_{j \\in n} f_{e,j}(|\\mathbf g|)\\,
          e^{-2\\pi^2 u_j^2 g^2}\\, e^{-2\\pi i\\,\\mathbf g\\cdot\\mathbf r_j},

   in V Å, where :math:`A` is the area of the cell and :math:`u_j` an optional
   static Debye-Waller RMS displacement. :math:`f_e` comes from the Lobato &
   Van Dyck (2014) or Kirkland (2010) fits, or from the Mott-Bethe transform of
   PyTex's X-ray form factors.
2. **Transmission.** :math:`t_n = \\exp(i\\sigma v_n)`, band-limited to two
   thirds of the Nyquist frequency so that the product :math:`\\psi t` cannot
   alias.
3. **Propagation.** :math:`\\psi_{n+1} = \\mathcal F^{-1}\\{\\mathcal F\\{\\psi_n t_n\\}
   P\\}` with the Fresnel propagator
   :math:`P = \\exp(-i\\pi\\lambda g^2\\Delta z - 2\\pi i\\Delta z\\,\\mathbf g\\cdot
   \\tan\\boldsymbol\\theta)` carrying the same band limit and an optional beam
   tilt :math:`\\boldsymbol\\theta`.
4. **Imaging.** The objective lens multiplies the Fourier transform of the exit
   wave by :math:`A(g)E(g)e^{-i\\chi(\\mathbf g)}` (quasi-coherent), or the
   temporal coherence is integrated exactly as an incoherent sum over a
   Gaussian defocus spread.
5. **Thermal diffuse scattering** (optional) by frozen phonons: independent
   Gaussian displacements of every atom, the whole calculation repeated per
   configuration, and *intensities* averaged.

The theory, derivations and validation are in
``docs/site/theory/multislice_hrtem.md``.

References
----------
- Cowley, J. M. & Moodie, A. F. (1957). Acta Cryst. 10, 609-619.
- Ishizuka, K. & Uyeda, N. (1977). Acta Cryst. A33, 740-749.
- Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, 2nd ed.,
  Springer.
- Lobato, I. & Van Dyck, D. (2014). Acta Cryst. A70, 636-649.
- Loane, R. F., Xu, P. & Silcox, J. (1991). Acta Cryst. A47, 267-278.
- Madsen, J. & Susi, T. (2021). The abTEM code: transmission electron
  microscopy from first principles. Open Research Europe 1, 24.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from enum import StrEnum
from functools import lru_cache
from importlib import resources
from typing import Any

import numpy as np
import scipy.fft
import scipy.special

from pytex.core._chemistry import atomic_number
from pytex.core.notation import format_miller_indices
from pytex.core.progress import report
from pytex.diffraction.hrem import (
    _ELECTRON_MASS_KG,
    _ELEMENTARY_CHARGE_C,
    _PLANCK_CONSTANT_JS,
    AtomicSnapshot,
    HREMSimulationResult,
    MicroscopeAberrations,
    relativistic_interaction_parameter_inv_v_angstrom,
    relativistic_wavelength_angstrom,
)

__all__ = [
    "ANTIALIAS_CUTOFF_FRACTION",
    "POTENTIAL_PREFACTOR_V_ANGSTROM2",
    "DefocusThicknessMap",
    "FocalSeries",
    "HREMEngine",
    "MultisliceExitWave",
    "MultisliceGrid",
    "PotentialParametrization",
    "SlicedPotential",
    "TemporalCoherence",
    "ZoneAxisCell",
    "antialias_aperture",
    "focal_integration_nodes",
    "fresnel_propagator",
    "hrtem_image",
    "multislice",
    "parametrized_electron_scattering_factor",
    "periodic_slab",
    "simulate_multislice_hrem",
    "slice_potential",
    "zone_axis_cell",
]

#: :math:`h^2 / (2\pi m_0 e)` in V Å\ :sup:`2`: turns an electron scattering
#: factor (Å) into a potential (V Å\ :sup:`3` per unit reciprocal volume). Equal
#: to :math:`2\pi a_0 e / (4\pi\varepsilon_0)` = 47.878 V Å\ :sup:`2` (Kirkland 2010,
#: eq. 5.9).
POTENTIAL_PREFACTOR_V_ANGSTROM2 = (
    _PLANCK_CONSTANT_JS**2 / (2.0 * math.pi * _ELECTRON_MASS_KG * _ELEMENTARY_CHARGE_C) * 1e20
)

#: The band limit as a fraction of the Nyquist frequency :math:`1/(2\Delta x)`.
#: Two thirds is the largest limit for which the product of two band-limited
#: functions - the wave and the transmission function - cannot alias back into
#: the band (Kirkland 2010, sec. 6.8). abTEM uses the same value.
ANTIALIAS_CUTOFF_FRACTION = 2.0 / 3.0

#: Width of the raised-cosine edge of the band limit, in units of
#: :math:`1/\Delta x`, as abTEM's default ``antialias.taper``.
ANTIALIAS_TAPER = 0.01

#: The most Gauss-Hermite nodes focal integration uses. With the node rule of
#: `focal_integration_nodes` it resolves a defocus phase rate a <= 44.5; see
#: `hrtem_image`.
MAX_FOCAL_INTEGRATION_POINTS = 512

_DATA_PACKAGE = "pytex.diffraction._data"
_TABLE_FILENAME = "electron_potential_parametrizations.json"

#: Above this many grid points a request is refused: a 4096 x 4096 complex wave
#: is 256 MB, and the algorithm holds several.
MAX_GRID_POINTS = 4096 * 4096


class PotentialParametrization(StrEnum):
    """Which independent-atom electron scattering factor builds the potential.

    ``LOBATO``
        Lobato & Van Dyck (2014): five hydrogen-like terms, fitted to
        Hartree-Fock charge densities to high accuracy at all angles and with the
        correct asymptotics. abTEM's default, and PyTex's.
    ``KIRKLAND``
        Kirkland (2010), Appendix C: three Lorentzians and three Gaussians,
        fitted to relativistic Hartree-Fock scattering factors.
    ``MOTT_BETHE``
        The Mott-Bethe transform of PyTex's tabulated X-ray form factors,
        `pytex.diffraction.scattering.electron_scattering_factors` - the same
        :math:`f_e` the Bloch-wave solver of `pytex.diffraction.dynamical`
        uses, which is what lets the two methods be compared on identical input.
    """

    LOBATO = "lobato"
    KIRKLAND = "kirkland"
    MOTT_BETHE = "mott_bethe"


class TemporalCoherence(StrEnum):
    """How partial temporal coherence (focal spread) enters an HRTEM image.

    ``QUASI_COHERENT``
        Frank's envelope :math:`E_c(g)` multiplies the transfer of the *wave*.
        Exact for the linear (weak-object) image terms, and what abTEM applies;
        it misrepresents the non-linear terms of a strong object.
    ``FOCAL_INTEGRATION``
        The image is the Gauss-Hermite weighted incoherent sum of images at
        defoci spread about the nominal one with the focal-spread standard
        deviation. Exact for every term, at the cost of one image per node.
    """

    QUASI_COHERENT = "quasi_coherent"
    FOCAL_INTEGRATION = "focal_integration"


class HREMEngine(StrEnum):
    """Which engine `pytex.adapters.abtem.simulate_hrem` uses.

    ``MULTISLICE``
        PyTex's multislice (this module) - the default.
    ``ABTEM``
        abTEM, through `pytex.adapters.abtem`.
    ``PHASE_OBJECT``
        The single-plane phase-object approximation.
    """

    MULTISLICE = "multislice"
    ABTEM = "abtem"
    PHASE_OBJECT = "phase_object"


@lru_cache(maxsize=1)
def _parametrization_table() -> dict[str, dict[str, np.ndarray]]:
    payload = json.loads(
        resources.files(_DATA_PACKAGE).joinpath(_TABLE_FILENAME).read_text(encoding="utf-8")
    )
    return {
        name: {symbol: np.asarray(rows, dtype=np.float64) for symbol, rows in table.items()}
        for name, table in payload["coefficients"].items()
    }


def parametrized_electron_scattering_factor(
    species: str,
    g_inv_angstrom: np.ndarray | float,
    parametrization: PotentialParametrization | str = PotentialParametrization.LOBATO,
) -> np.ndarray:
    """Independent-atom electron scattering factor :math:`f_e(g)` in ångström.

    What it does
        Evaluates the chosen parametrization at the scattering-vector magnitude
        :math:`g = 2\\sin\\theta/\\lambda`. The zero-angle value fixes the
        integrated projected potential of the atom,
        :math:`\\int v\\,\\mathrm d^2r = (h^2/2\\pi m_0 e)\\,f_e(0)`, and hence the
        mean inner potential of a crystal.

    When to use it
        To inspect or compare the potentials a multislice run uses; the engine
        calls it for every species on its Fourier grid.

    Parameters
    ----------
    species:
        Element symbol.
    g_inv_angstrom:
        :math:`g` in Å\\ :sup:`-1`, any shape.
    parametrization:
        A `PotentialParametrization` or its value.

    Returns
    -------
    np.ndarray
        :math:`f_e(g)` in Å, the shape of ``g_inv_angstrom``.

    Raises
    ------
    ValueError
        If the element has no coefficients in the chosen table.
    """

    method = PotentialParametrization(parametrization)
    g = np.asarray(g_inv_angstrom, dtype=np.float64)
    g2 = g * g
    if method is PotentialParametrization.MOTT_BETHE:
        from pytex.diffraction.scattering import electron_scattering_factors, tabulated_species

        if species not in tabulated_species():
            raise ValueError(f"No Mott-Bethe scattering factor is tabulated for {species!r}.")
        flat_mb = electron_scattering_factors(species, 0.5 * g.reshape(-1))
        return np.asarray(flat_mb).reshape(g.shape)
    table = _parametrization_table()[method.value]
    if species not in table:
        raise ValueError(
            f"No {method.value} electron scattering factor is tabulated for {species!r}."
        )
    p = table[species]
    if method is PotentialParametrization.LOBATO:
        a, b = p[0][:, None], p[1][:, None]
        flat = g2.reshape(1, -1)
        f = np.sum(a * (2.0 + b * flat) / (1.0 + b * flat) ** 2, axis=0)
    else:
        a, b, c, d = (row[:, None] for row in p)
        flat = g2.reshape(1, -1)
        f = np.sum(a / (flat + b) + c * np.exp(-d * flat), axis=0)
    return np.asarray(f.reshape(g.shape))


@dataclass(frozen=True, slots=True)
class MultisliceGrid:
    """The lateral sampling grid shared by the potential, the wave and the image.

    Arrays on this grid are indexed ``[row, column] = [y, x]`` with row 0 at
    ``y = 0``, the convention of `HREMSimulationResult`. The grid is periodic:
    a multislice calculation on it is the calculation for an infinite lateral
    repetition of the cell.

    Attributes
    ----------
    shape:
        ``(ny, nx)`` grid points.
    extent_angstrom:
        ``(Lx, Ly)`` of the periodic cell in Å.
    """

    shape: tuple[int, int]
    extent_angstrom: tuple[float, float]

    def __post_init__(self) -> None:
        ny, nx = (int(v) for v in self.shape)
        lx, ly = (float(v) for v in self.extent_angstrom)
        if ny < 2 or nx < 2:
            raise ValueError(f"A multislice grid needs at least 2 x 2 points, got {self.shape}.")
        if ny * nx > MAX_GRID_POINTS:
            raise ValueError(
                f"A {nx} x {ny} grid exceeds the {MAX_GRID_POINTS}-point limit; "
                "coarsen the sampling or reduce the cell."
            )
        if not (lx > 0.0 and ly > 0.0 and math.isfinite(lx) and math.isfinite(ly)):
            raise ValueError(f"The cell extent must be positive and finite, got {(lx, ly)}.")
        object.__setattr__(self, "shape", (ny, nx))
        object.__setattr__(self, "extent_angstrom", (lx, ly))

    @classmethod
    def from_sampling(
        cls, extent_angstrom: tuple[float, float], sampling_angstrom: float
    ) -> MultisliceGrid:
        """The grid whose pixel is at most ``sampling_angstrom`` on each axis.

        The point counts are ``ceil(L / sampling)`` rounded up to a size the FFT
        factorizes quickly, so the delivered pixel is never coarser than asked.
        """

        if not sampling_angstrom > 0.0:
            raise ValueError(f"The sampling must be positive, got {sampling_angstrom}.")
        lx, ly = extent_angstrom
        nx = scipy.fft.next_fast_len(max(2, math.ceil(lx / sampling_angstrom - 1e-9)))
        ny = scipy.fft.next_fast_len(max(2, math.ceil(ly / sampling_angstrom - 1e-9)))
        return cls(shape=(ny, nx), extent_angstrom=(lx, ly))

    @property
    def sampling_angstrom(self) -> tuple[float, float]:
        """Pixel size ``(dx, dy)`` in Å."""
        return (self.extent_angstrom[0] / self.shape[1], self.extent_angstrom[1] / self.shape[0])

    def spatial_frequencies(self) -> tuple[np.ndarray, np.ndarray]:
        """Fourier-grid frequencies ``(gx, gy)`` in Å\\ :sup:`-1`, broadcastable to ``shape``."""
        dx, dy = self.sampling_angstrom
        gx = np.fft.fftfreq(self.shape[1], d=dx)[None, :]
        gy = np.fft.fftfreq(self.shape[0], d=dy)[:, None]
        return gx, gy

    def frequency_magnitude(self) -> np.ndarray:
        """:math:`|\\mathbf g|` on the Fourier grid in Å\\ :sup:`-1`."""
        gx, gy = self.spatial_frequencies()
        return np.asarray(np.hypot(gx, gy))

    @property
    def antialias_cutoff_inv_angstrom(self) -> float:
        """The band limit :math:`g_{\\max} = \\tfrac23 / (2\\max(\\Delta x, \\Delta y))`."""
        return ANTIALIAS_CUTOFF_FRACTION / (2.0 * max(self.sampling_angstrom))

    def max_scattering_angle_mrad(self, energy_kev: float) -> float:
        """Largest scattering semi-angle the band-limited grid represents, in mrad."""
        return 1e3 * relativistic_wavelength_angstrom(energy_kev) * (
            self.antialias_cutoff_inv_angstrom
        )


def antialias_aperture(grid: MultisliceGrid) -> np.ndarray:
    """The band-limiting aperture on the Fourier grid, 1 inside and 0 outside.

    A raised-cosine edge of width ``ANTIALIAS_TAPER / max(dx, dy)`` ends at the
    cutoff :math:`g_{\\max}` of `MultisliceGrid.antialias_cutoff_inv_angstrom`,
    exactly as abTEM's ``antialias_aperture``.
    """

    cutoff = grid.antialias_cutoff_inv_angstrom
    taper = ANTIALIAS_TAPER / max(grid.sampling_angstrom)
    g = grid.frequency_magnitude()
    edge = 0.5 * (1.0 + np.cos(math.pi * (g - cutoff + taper) / taper))
    aperture = np.where(g > cutoff - taper, edge, 1.0)
    aperture[g > cutoff] = 0.0
    return aperture


def fresnel_propagator(
    grid: MultisliceGrid,
    energy_kev: float,
    thickness_angstrom: float,
    tilt_mrad: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Band-limited Fresnel propagator over ``thickness_angstrom`` of vacuum.

    .. math::

       P(\\mathbf g) = A_{\\mathrm{aa}}(g)\\,
       \\exp\\!\\left(-i\\pi\\lambda g^2\\Delta z
       - 2\\pi i\\,\\Delta z\\,(g_x\\tan\\theta_x + g_y\\tan\\theta_y)\\right).

    Multiplying the Fourier transform of a wave by it advances the wave by
    :math:`\\Delta z` along +z in the paraxial (Fresnel) approximation. Two
    propagators compose: ``P(dz1) * P(dz2) == P(dz1 + dz2)`` inside the band.
    """

    lam = relativistic_wavelength_angstrom(energy_kev)
    gx, gy = grid.spatial_frequencies()
    phase = -math.pi * lam * thickness_angstrom * (gx * gx + gy * gy)
    tx, ty = (math.tan(1e-3 * float(t)) for t in tilt_mrad)
    if tx != 0.0 or ty != 0.0:
        phase = phase - 2.0 * math.pi * thickness_angstrom * (gx * tx + gy * ty)
    return np.asarray(antialias_aperture(grid) * np.exp(1j * phase))


def slice_potential(
    species: Sequence[str],
    positions_xy: np.ndarray,
    grid: MultisliceGrid,
    parametrization: PotentialParametrization | str = PotentialParametrization.LOBATO,
    debye_waller_sigma_angstrom: Mapping[str, float] | float | None = None,
) -> np.ndarray:
    """Projected potential :math:`v(x, y)` of one slice, in V Å.

    What it does
        Sums the infinite projections of the atoms given, in Fourier space, on
        the periodic ``grid``. The result integrates to
        :math:`(h^2/2\\pi m_0 e)\\sum_j f_{e,j}(0)` over the cell, the identity
        the tests assert.

    Parameters
    ----------
    species:
        Element symbol per atom.
    positions_xy:
        ``(N, 2)`` lateral positions in Å.
    grid:
        The sampling grid.
    parametrization:
        Scattering-factor fit.
    debye_waller_sigma_angstrom:
        Static RMS displacement :math:`u` per element (or one value for all):
        each scattering factor is damped by :math:`e^{-2\\pi^2u^2g^2}`, the
        time-averaged potential of a vibrating atom.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` real array in V Å.
    """

    positions = np.asarray(positions_xy, dtype=np.float64).reshape(-1, 2)
    ny, nx = grid.shape
    if positions.shape[0] == 0:
        return np.zeros((ny, nx), dtype=np.float64)
    method = PotentialParametrization(parametrization)
    gx, gy = grid.spatial_frequencies()
    g = np.hypot(gx, gy)
    lx, ly = grid.extent_angstrom
    total = np.zeros((ny, nx), dtype=np.complex128)
    labels = np.asarray(species)
    for symbol in sorted(set(species)):
        chosen = positions[labels == symbol]
        # The structure factor of the species factorizes over x and y, so it is
        # one (ny, N) @ (N, nx) product rather than an (ny, nx, N) broadcast.
        ex = np.exp(-2j * math.pi * chosen[:, 0:1] * gx)  # (N, nx)
        ey = np.exp(-2j * math.pi * chosen[:, 1:2] * gy.reshape(1, -1))  # (N, ny)
        structure = ey.T @ ex
        f = parametrized_electron_scattering_factor(symbol, g, method)
        sigma = _sigma_for(symbol, debye_waller_sigma_angstrom)
        if sigma > 0.0:
            f = f * np.exp(-2.0 * math.pi**2 * sigma**2 * g * g)
        total += f * structure
    total *= POTENTIAL_PREFACTOR_V_ANGSTROM2 / (lx * ly)
    return np.asarray(scipy.fft.ifft2(total, norm="forward").real, dtype=np.float64)


def _sigma_for(symbol: str, sigmas: Mapping[str, float] | float | None) -> float:
    if sigmas is None:
        return 0.0
    if isinstance(sigmas, Mapping):
        return float(sigmas.get(symbol, 0.0))
    return float(sigmas)


@dataclass(frozen=True, slots=True)
class SlicedPotential:
    """A snapshot cut into slices along the beam, ready for multislice.

    The atoms of each slice are kept rather than the potential arrays: a slice's
    potential is built when the wave reaches it, and slices whose projected atom
    sets are identical - every repeat of a crystal whose period the slicing
    divides - are recognised by `slice_keys` and built once.

    Attributes
    ----------
    snapshot:
        The specimen, in an orthogonal box with the beam along +z.
    grid:
        Lateral sampling.
    slice_edges_angstrom:
        ``num_slices + 1`` z edges from 0 to the box height.
    slice_of_atom:
        Slice index of each atom (the one holding its centre).
    parametrization:
        Scattering-factor fit.
    debye_waller_sigma_angstrom:
        Static RMS displacement per element, or one value for all, or None.
    """

    snapshot: AtomicSnapshot
    grid: MultisliceGrid
    slice_edges_angstrom: np.ndarray
    slice_of_atom: np.ndarray
    parametrization: PotentialParametrization = PotentialParametrization.LOBATO
    debye_waller_sigma_angstrom: Mapping[str, float] | float | None = None
    slice_keys: tuple[int, ...] = field(default=(), compare=False)

    @classmethod
    def from_snapshot(
        cls,
        snapshot: AtomicSnapshot,
        grid: MultisliceGrid,
        slice_thickness_angstrom: float = 1.0,
        parametrization: PotentialParametrization | str = PotentialParametrization.LOBATO,
        debye_waller_sigma_angstrom: Mapping[str, float] | float | None = None,
    ) -> SlicedPotential:
        """Slice ``snapshot`` into equal slices no thicker than ``slice_thickness_angstrom``.

        The box height :math:`L_z` is divided into
        :math:`n = \\lceil L_z/\\Delta z \\rceil` equal slices (abTEM's rule), so
        the delivered thickness :math:`L_z/n` never exceeds the request.
        """

        if not slice_thickness_angstrom > 0.0:
            raise ValueError(
                f"The slice thickness must be positive, got {slice_thickness_angstrom}."
            )
        lz = _box_height(snapshot)
        count = max(1, math.ceil(lz / slice_thickness_angstrom - 1e-9))
        edges = np.linspace(0.0, lz, count + 1)
        z = snapshot.positions[:, 2]
        # An atom on a slice boundary belongs to the slice that begins there; the
        # tolerance keeps floating-point rounding of z from deciding otherwise.
        index = np.clip(np.floor(z / (lz / count) + 1e-6).astype(np.int64), 0, count - 1)
        method = PotentialParametrization(parametrization)
        for symbol in set(snapshot.species):
            atomic_number(symbol)
            parametrized_electron_scattering_factor(symbol, 0.0, method)
        return cls(
            snapshot=snapshot,
            grid=grid,
            slice_edges_angstrom=edges,
            slice_of_atom=index,
            parametrization=method,
            debye_waller_sigma_angstrom=debye_waller_sigma_angstrom,
            slice_keys=_slice_keys(snapshot, index, count, grid),
        )

    @property
    def num_slices(self) -> int:
        """Number of slices."""
        return int(self.slice_edges_angstrom.size - 1)

    @property
    def slice_thicknesses_angstrom(self) -> np.ndarray:
        """Thickness of each slice in Å."""
        return np.diff(self.slice_edges_angstrom)

    @property
    def thickness_angstrom(self) -> float:
        """Total thickness traversed, the box height, in Å."""
        return float(self.slice_edges_angstrom[-1])

    @property
    def unique_slice_count(self) -> int:
        """Number of distinct slice potentials the calculation builds."""
        return len(set(self.slice_keys))

    def projected_potential(self, index: int) -> np.ndarray:
        """Projected potential of slice ``index`` in V Å, ``(ny, nx)``."""
        chosen = self.slice_of_atom == index
        return slice_potential(
            [s for s, keep in zip(self.snapshot.species, chosen, strict=True) if keep],
            self.snapshot.positions[chosen, :2],
            self.grid,
            self.parametrization,
            self.debye_waller_sigma_angstrom,
        )

    def total_projected_potential(self) -> np.ndarray:
        """Projected potential of the whole specimen in V Å (the phase-object limit)."""
        return slice_potential(
            self.snapshot.species,
            self.snapshot.positions[:, :2],
            self.grid,
            self.parametrization,
            self.debye_waller_sigma_angstrom,
        )

    def mean_inner_potential_volt(self) -> float:
        """Mean potential of the filled box in V: :math:`\\langle v \\rangle / L_z`.

        For a slab filling the box this is the mean inner potential of the
        crystal in the independent-atom approximation.
        """

        prefactor = POTENTIAL_PREFACTOR_V_ANGSTROM2
        f0 = sum(
            float(parametrized_electron_scattering_factor(s, 0.0, self.parametrization))
            for s in self.snapshot.species
        )
        lx, ly = self.grid.extent_angstrom
        return prefactor * f0 / (lx * ly * self.thickness_angstrom)


def _direction(indices: Any) -> str:
    """A specific direction [uvw] in the plain notation of `pytex.core.notation`."""
    return format_miller_indices(
        tuple(int(v) for v in indices), family="direction", style="plain", scope="specific"
    )


@dataclass(frozen=True, slots=True)
class ZoneAxisCell:
    """An orthogonal, lattice-periodic cell with the zone axis along +z.

    A multislice calculation is periodic in x and y, so the lateral box must be
    an exact repeat of the crystal - otherwise every cell edge is a seam, the
    reflections sit off the Fourier grid, and forbidden reflections appear.
    This cell is built from three mutually perpendicular *lattice* vectors,
    so it tiles the crystal exactly.

    Attributes
    ----------
    zone_axis:
        The requested ``[uvw]``.
    vectors_uvw:
        ``(3, 3)`` integer rows ``t1, t2, t3`` in direct-lattice components;
        ``t3`` is the primitive lattice vector along ``[uvw]``, and
        ``(t1 x t2) . t3 > 0``.
    vectors_cartesian:
        The same rows in Å, in the crystal's Cartesian frame.
    """

    zone_axis: tuple[int, int, int]
    vectors_uvw: np.ndarray
    vectors_cartesian: np.ndarray

    @property
    def lengths_angstrom(self) -> tuple[float, float, float]:
        """Lengths ``|t1|, |t2|, |t3|`` - the box edges of one cell."""
        n = np.linalg.norm(self.vectors_cartesian, axis=1)
        return (float(n[0]), float(n[1]), float(n[2]))

    @property
    def rotation(self) -> np.ndarray:
        """Rotation from crystal Cartesian coordinates to the box frame (rows t1, t2, t3 unit)."""
        return np.asarray(
            self.vectors_cartesian / np.linalg.norm(self.vectors_cartesian, axis=1)[:, None]
        )

    @property
    def cells_per_box(self) -> int:
        """Number of crystallographic unit cells in one box, ``|det(vectors_uvw)|``."""
        return round(abs(float(np.linalg.det(self.vectors_uvw))))

    def describe(self) -> str:
        """The cell as a crystallographer writes it."""
        t1, t2, t3 = (_direction(row) for row in self.vectors_uvw)
        l1, l2, l3 = self.lengths_angstrom
        return (
            f"Orthogonal periodic cell for zone axis "
            f"{_direction(self.zone_axis)}: x ∥ {t1} ({l1:.4f} Å), "
            f"y ∥ {t2} ({l2:.4f} Å), z ∥ {t3} ({l3:.4f} Å), "
            f"{self.cells_per_box} unit cells per box."
        )


def zone_axis_cell(
    phase: Any, zone_axis: tuple[int, int, int], max_index: int = 6
) -> ZoneAxisCell:
    """The smallest orthogonal lattice-periodic cell with ``[uvw]`` along +z.

    What it does
        Takes ``t3`` as the primitive lattice vector along ``[uvw]`` and
        searches integer vectors with components up to ``max_index`` for the
        shortest ``t1`` perpendicular to it and then the shortest ``t2``
        perpendicular to both. For cubic ``[110]`` this gives ``[001]``,
        ``[1-10]``, ``[110]``; for hexagonal ``[001]`` it gives ``[100]``,
        ``[120]``, ``[001]``.

    When to use it
        Before any periodic image or exit-wave simulation of a crystal. The
        box a naive rotation of the conventional cell produces is not periodic
        along x and y, and the simulation then sees a seam at every edge.

    Raises
    ------
    ValueError
        For a zero zone axis, or when no orthogonal lattice pair exists within
        ``max_index`` (a low-symmetry lattice at a general zone axis): an
        approximate, strained box would simulate a different crystal.
    """

    uvw = np.asarray(zone_axis, dtype=np.int64)
    if uvw.shape != (3,) or not np.any(uvw):
        raise ValueError(f"The zone axis must be three integers, not all zero; got {zone_axis}.")
    uvw = uvw // math.gcd(*(int(v) for v in uvw))
    basis = phase.lattice.direct_basis()
    b = np.vstack([basis.vector(0), basis.vector(1), basis.vector(2)])
    t3 = uvw @ b
    span = np.arange(-max_index, max_index + 1)
    candidates = np.stack(np.meshgrid(span, span, span, indexing="ij"), axis=-1).reshape(-1, 3)
    candidates = candidates[np.any(candidates != 0, axis=1)]
    cart = candidates @ b
    lengths = np.linalg.norm(cart, axis=1)

    def perpendicular(to: np.ndarray) -> np.ndarray:
        return np.asarray(np.abs(cart @ to) <= 1e-6 * lengths * float(np.linalg.norm(to)))

    def shortest(mask: np.ndarray) -> int:
        chosen = np.flatnonzero(mask)
        if chosen.size == 0:
            raise ValueError(
                f"No lattice vector with components up to {max_index} is perpendicular to "
                f"the zone axis {tuple(int(v) for v in uvw)}; this lattice has no orthogonal "
                "periodic cell for it within that range."
            )
        # Shortest first; among equals prefer fewer negative components, then
        # the lexicographically largest, so the choice is deterministic.
        order = sorted(
            chosen.tolist(),
            key=lambda i: (
                round(float(lengths[i]), 6),
                int(np.sum(candidates[i] < 0)),
                tuple(int(v) for v in -candidates[i]),
            ),
        )
        return int(order[0])

    i1 = shortest(perpendicular(t3))
    i2 = shortest(perpendicular(t3) & perpendicular(cart[i1]))
    v1, v2 = candidates[i1].copy(), candidates[i2].copy()
    if float(np.dot(np.cross(v1 @ b, v2 @ b), t3)) < 0.0:
        v2 = -v2
    vectors = np.vstack([v1, v2, uvw])
    return ZoneAxisCell(
        zone_axis=(int(zone_axis[0]), int(zone_axis[1]), int(zone_axis[2])),
        vectors_uvw=vectors,
        vectors_cartesian=vectors @ b,
    )


def periodic_slab(
    phase: Any,
    zone_axis: tuple[int, int, int] = (0, 0, 1),
    lateral_repeats: tuple[int, int] = (1, 1),
    *,
    thickness_angstrom: float | None = None,
    beam_repeats: int | None = None,
    max_index: int = 6,
    label: str | None = None,
) -> tuple[AtomicSnapshot, ZoneAxisCell, tuple[int, int, int]]:
    """A crystal slab in an exact orthogonal periodic box, ``[uvw]`` along +z.

    What it does
        Tiles `zone_axis_cell` ``lateral_repeats`` times across the beam and
        ``beam_repeats`` times along it (or the fewest times that reach
        ``thickness_angstrom``), fills the box with every atom of the phase's
        unit cell whose fractional box coordinates lie in ``[0, 1)``, and
        expresses it in the box frame. The atom count is checked against the
        number of unit cells the box holds, so a slab is never silently
        incomplete.

    Returns
    -------
    tuple
        ``(snapshot, cell, repeats)`` with ``repeats = (m1, m2, m3)``.
    """

    if phase.unit_cell is None or not phase.unit_cell.sites:
        raise ValueError(f"Phase '{phase.name}' must carry a unit cell with atomic sites.")
    cell = zone_axis_cell(phase, zone_axis, max_index)
    m1, m2 = (int(v) for v in lateral_repeats)
    if m1 < 1 or m2 < 1:
        raise ValueError(f"Lateral repeats must be at least 1, got {lateral_repeats}.")
    l1, l2, l3 = cell.lengths_angstrom
    if beam_repeats is not None:
        m3 = int(beam_repeats)
        if m3 < 1:
            raise ValueError(f"Beam repeats must be at least 1, got {beam_repeats}.")
    elif thickness_angstrom is not None:
        if not thickness_angstrom > 0.0:
            raise ValueError("The specimen thickness must be strictly positive.")
        m3 = max(1, math.ceil(thickness_angstrom / l3 - 1e-9))
    else:
        m3 = 1
    box_uvw = cell.vectors_uvw * np.array([[m1], [m2], [m3]])
    inverse = np.linalg.inv(box_uvw.astype(np.float64))
    corners = (
        np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=np.float64)
        @ box_uvw
    )
    lo = np.floor(corners.min(axis=0)).astype(np.int64) - 1
    hi = np.ceil(corners.max(axis=0)).astype(np.int64) + 1
    translations = np.stack(
        np.meshgrid(*(np.arange(a, b + 1) for a, b in zip(lo, hi, strict=True)), indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    sites = phase.unit_cell.sites
    eps = 1e-7
    extent = np.array([m1 * l1, m2 * l2, m3 * l3])
    species: list[str] = []
    positions: list[np.ndarray] = []
    for site in sites:
        u = (translations + np.asarray(site.fractional_coordinates, dtype=np.float64)) @ inverse
        inside = np.all((u >= -eps) & (u < 1.0 - eps), axis=1)
        chosen = np.clip(u[inside], 0.0, None)
        positions.append(chosen * extent)
        species.extend([site.species] * int(chosen.shape[0]))
    xyz = np.vstack(positions)
    expected = cell.cells_per_box * m1 * m2 * m3 * len(sites)
    if xyz.shape[0] != expected:
        raise ValueError(
            f"The periodic slab holds {xyz.shape[0]} atoms where {expected} were expected; "
            "the unit cell's sites may not be a complete cell."
        )
    snapshot = AtomicSnapshot(
        species=tuple(species),
        positions=xyz,
        cell=np.diag(extent),
        periodicity=(True, True, False),
        label=label
        or f"{phase.name} {_direction(cell.zone_axis)} slab ({m1}×{m2}×{m3} periodic cells)",
    )
    return snapshot, cell, (m1, m2, m3)


def _box_height(snapshot: AtomicSnapshot) -> float:
    height = float(snapshot.cell[2, 2])
    if height > 0.0:
        return height
    return max(float(np.ptp(snapshot.positions[:, 2])) + 2.0, 1.0)


def _slice_keys(
    snapshot: AtomicSnapshot, index: np.ndarray, count: int, grid: MultisliceGrid
) -> tuple[int, ...]:
    """A hash per slice of its projected atoms, equal for identical projections.

    Under infinite projection only an atom's species and lateral position
    matter, so two slices with the same sorted ``(species, x, y)`` set have the
    same potential. Positions are compared on a grid of 1e-4 Å.
    """

    lx, ly = grid.extent_angstrom
    xy = np.round(
        np.mod(snapshot.positions[:, :2], (lx, ly)) * 1e4
    ).astype(np.int64)
    species = np.asarray(snapshot.species)
    keys: list[int] = []
    for n in range(count):
        chosen = np.flatnonzero(index == n)
        items = sorted(
            zip(species[chosen].tolist(), xy[chosen, 0].tolist(), xy[chosen, 1].tolist(),
                strict=True)
        )
        keys.append(hash(tuple(items)))
    return tuple(keys)


def _transmission(
    potential: np.ndarray, sigma: float, aperture: np.ndarray
) -> np.ndarray:
    t = np.exp(1j * sigma * potential)
    return np.asarray(scipy.fft.ifft2(scipy.fft.fft2(t, workers=-1) * aperture, workers=-1))


def _thermal_snapshot(
    snapshot: AtomicSnapshot,
    sigmas: Mapping[str, float] | float,
    rng: np.random.Generator,
) -> AtomicSnapshot:
    rms = np.array([_sigma_for(s, sigmas) for s in snapshot.species], dtype=np.float64)
    moved = snapshot.positions + rng.normal(size=snapshot.positions.shape) * rms[:, None]
    lz = _box_height(snapshot)
    moved[:, 2] = np.clip(moved[:, 2], 0.0, np.nextafter(lz, 0.0))
    return replace(snapshot, positions=moved)


@dataclass(frozen=True, slots=True)
class MultisliceExitWave:
    """Exit waves of a multislice run, at every depth asked for and every configuration.

    Attributes
    ----------
    waves:
        Complex array ``(configurations, depths, ny, nx)``. The incident plane
        wave has amplitude 1, so ``|psi|^2`` averages to 1 in vacuum.
    depths_angstrom:
        Depth below the entrance surface of each stored plane, in Å.
    grid:
        The lateral grid.
    energy_kev:
        Beam energy.
    tilt_mrad:
        Beam tilt ``(theta_x, theta_y)`` in mrad.
    slice_thickness_angstrom:
        The (uniform) slice thickness used.
    num_slices, unique_slices:
        Slices traversed, and distinct slice potentials built (per configuration).
    parametrization:
        Scattering-factor fit.
    retained_intensity:
        Mean of ``|psi|^2`` at each depth, averaged over configurations. It
        falls below 1 only by what the band limit removes, so a value well
        below 1 says the sampling is too coarse for the scattering.
    frozen_phonon_sigma_angstrom:
        RMS displacement of the frozen-phonon ensemble, or None for a static
        lattice.
    snapshot:
        The specimen.
    """

    waves: np.ndarray
    depths_angstrom: np.ndarray
    grid: MultisliceGrid
    energy_kev: float
    tilt_mrad: tuple[float, float]
    slice_thickness_angstrom: float
    num_slices: int
    unique_slices: int
    parametrization: PotentialParametrization
    retained_intensity: np.ndarray
    frozen_phonon_sigma_angstrom: Mapping[str, float] | float | None
    snapshot: AtomicSnapshot

    @property
    def configurations(self) -> int:
        """Number of frozen-phonon configurations (1 for a static lattice)."""
        return int(self.waves.shape[0])

    @property
    def wavelength_angstrom(self) -> float:
        """Relativistic electron wavelength."""
        return relativistic_wavelength_angstrom(self.energy_kev)

    def depth_index(self, depth_angstrom: float | None) -> int:
        """Index of the stored plane nearest ``depth_angstrom`` (the last if None)."""
        if depth_angstrom is None:
            return int(self.depths_angstrom.size - 1)
        return int(np.argmin(np.abs(self.depths_angstrom - float(depth_angstrom))))

    def wave(self, depth_index: int = -1) -> np.ndarray:
        """The configuration-averaged (elastic, coherent) exit wave at one depth."""
        return np.asarray(np.mean(self.waves[:, depth_index], axis=0))

    def diffraction_pattern(self, depth_index: int = -1) -> np.ndarray:
        """Far-field intensity :math:`\\langle|\\Psi(\\mathbf g)|^2\\rangle`, zero at the centre.

        Normalised so the pattern sums to the retained intensity. With frozen
        phonons the average is of intensities, so it holds the thermal diffuse
        background between the Bragg beams.
        """

        ny, nx = self.grid.shape
        spectra = np.abs(scipy.fft.fft2(self.waves[:, depth_index], workers=-1)) ** 2
        return np.asarray(np.fft.fftshift(np.mean(spectra, axis=0)) / (nx * ny) ** 2)

    def beam_intensities(
        self, g_xy_inv_angstrom: np.ndarray, depth_index: int | None = None
    ) -> np.ndarray:
        """Intensity of the beams at lateral reciprocal vectors ``g`` (Å\\ :sup:`-1`).

        Each ``g`` is taken at the nearest Fourier-grid point; a beam lies on
        the grid exactly when the cell is a whole number of lattice periods.
        With ``depth_index`` None the result has one row per stored depth.
        """

        g = np.asarray(g_xy_inv_angstrom, dtype=np.float64).reshape(-1, 2)
        lx, ly = self.grid.extent_angstrom
        ny, nx = self.grid.shape
        cols = np.mod(np.rint(g[:, 0] * lx).astype(np.int64), nx)
        rows = np.mod(np.rint(g[:, 1] * ly).astype(np.int64), ny)
        spectra = scipy.fft.fft2(self.waves, workers=-1, norm="forward")
        intensity = np.mean(np.abs(spectra[..., rows, cols]) ** 2, axis=0)
        if depth_index is None:
            return np.asarray(intensity)
        return np.asarray(intensity[depth_index])

    def image(
        self,
        aberrations: MicroscopeAberrations,
        depth_index: int = -1,
        temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
        focal_integration_points: int | None = None,
    ) -> np.ndarray:
        """HRTEM image intensity at one depth; see `hrtem_image`."""
        return hrtem_image(
            self.waves[:, depth_index],
            self.grid,
            aberrations,
            temporal_coherence=temporal_coherence,
            focal_integration_points=focal_integration_points,
        )

    def focal_series(
        self,
        aberrations: MicroscopeAberrations,
        defoci_angstrom: Iterable[float],
        depth_index: int = -1,
        temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
    ) -> FocalSeries:
        """Images at each defocus from this exit wave - no new multislice run.

        The exit wave does not depend on the objective lens, so a focal series
        of any length costs one multislice calculation and one pair of FFTs
        per image.
        """

        defoci = np.asarray(list(defoci_angstrom), dtype=np.float64)
        if defoci.size == 0:
            raise ValueError("A focal series needs at least one defocus.")
        images = np.stack(
            [
                self.image(
                    replace(aberrations, defocus_angstrom=float(df)),
                    depth_index,
                    temporal_coherence,
                )
                for df in defoci
            ]
        )
        return FocalSeries(
            defoci_angstrom=defoci,
            images=images,
            thickness_angstrom=float(self.depths_angstrom[depth_index]),
            aberrations=aberrations,
            temporal_coherence=TemporalCoherence(temporal_coherence),
            grid=self.grid,
            label=self.snapshot.label,
        )

    def defocus_thickness_map(
        self,
        aberrations: MicroscopeAberrations,
        defoci_angstrom: Iterable[float],
        temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
    ) -> DefocusThicknessMap:
        """Images at every (stored depth, defocus) pair: the classic HRTEM tableau."""
        defoci = np.asarray(list(defoci_angstrom), dtype=np.float64)
        if defoci.size == 0:
            raise ValueError("A defocus-thickness map needs at least one defocus.")
        images = np.stack(
            [
                self.focal_series(aberrations, defoci, index, temporal_coherence).images
                for index in range(self.depths_angstrom.size)
            ]
        )
        return DefocusThicknessMap(
            defoci_angstrom=defoci,
            thicknesses_angstrom=np.asarray(self.depths_angstrom, dtype=np.float64),
            images=images,
            aberrations=aberrations,
            temporal_coherence=TemporalCoherence(temporal_coherence),
            grid=self.grid,
            label=self.snapshot.label,
        )

    def hrem_result(
        self,
        aberrations: MicroscopeAberrations,
        depth_index: int = -1,
        temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
    ) -> HREMSimulationResult:
        """The image at one depth packaged as an `HREMSimulationResult`.

        Adds the power spectrum of the image and the one-dimensional CTF of the
        lens, so a multislice image reads the same as any other engine's.
        """

        image = self.image(aberrations, depth_index, temporal_coherence)
        dx, _dy = self.grid.sampling_angstrom
        fluctuation = image - np.mean(image)
        spectrum = np.fft.fftshift(np.abs(scipy.fft.fft2(fluctuation, workers=-1)) ** 2)
        return HREMSimulationResult(
            image=image,
            exit_wave=self.wave(depth_index),
            pixel_size_angstrom=float(dx),
            extent_angstrom=self.grid.extent_angstrom,
            power_spectrum=np.log10(1.0 + spectrum),
            ctf=aberrations.evaluate_ctf_1d(max_q_inv_angstrom=1.0 / (2.0 * dx)),
            aberrations=aberrations,
            snapshot=self.snapshot,
        )

    def describe(self) -> str:
        """Convention-explicit account of the calculation and its diagnostics."""
        ny, nx = self.grid.shape
        dx, dy = self.grid.sampling_angstrom
        lx, ly = self.grid.extent_angstrom
        lines = [
            (
                f"Multislice exit wave of {self.snapshot.label} ({self.snapshot.natoms} atoms) "
                f"at {self.energy_kev:.1f} kV (λ = {self.wavelength_angstrom:.5f} Å), "
                f"plane-wave illumination along +z"
                + (
                    f" tilted by ({self.tilt_mrad[0]:.2f}, {self.tilt_mrad[1]:.2f}) mrad."
                    if any(self.tilt_mrad)
                    else "."
                )
            ),
            (
                f"Grid {nx} × {ny} over {lx:.2f} × {ly:.2f} Å (Δx = {dx:.4f}, Δy = {dy:.4f} Å), "
                f"periodic; band limit {self.grid.antialias_cutoff_inv_angstrom:.3f} Å⁻¹ "
                f"(two thirds of Nyquist), i.e. scattering up to "
                f"{self.grid.max_scattering_angle_mrad(self.energy_kev):.1f} mrad."
            ),
            (
                f"{self.num_slices} slices of Δz = {self.slice_thickness_angstrom:.4f} Å "
                f"({self.unique_slices} distinct), infinite projection, "
                f"{self.parametrization.value} scattering factors (Cowley & Moodie 1957; "
                "Kirkland 2010)."
            ),
            (
                f"Stored depths: {', '.join(f'{t:.1f}' for t in self.depths_angstrom)} Å. "
                f"Retained intensity at the exit plane {float(self.retained_intensity[-1]):.4f}"
                " (1 is lossless; the band limit removes the rest)."
            ),
        ]
        if self.frozen_phonon_sigma_angstrom is not None and self.configurations > 1:
            lines.append(
                f"Thermal diffuse scattering by {self.configurations} frozen-phonon "
                f"configurations (Loane, Xu & Silcox 1991), intensities averaged."
            )
        if float(self.retained_intensity[-1]) < 0.95:
            lines.append(
                "More than 5 % of the intensity was scattered beyond the band limit: "
                "refine the sampling before trusting high-angle detail."
            )
        return " ".join(lines)


def multislice(
    snapshot: AtomicSnapshot,
    energy_kev: float,
    *,
    sampling_angstrom: float | None = 0.05,
    grid: MultisliceGrid | None = None,
    slice_thickness_angstrom: float = 1.0,
    parametrization: PotentialParametrization | str = PotentialParametrization.LOBATO,
    tilt_mrad: tuple[float, float] = (0.0, 0.0),
    exit_depths_angstrom: Iterable[float] | None = None,
    debye_waller_sigma_angstrom: Mapping[str, float] | float | None = None,
    frozen_phonon_sigma_angstrom: Mapping[str, float] | float | None = None,
    frozen_phonon_configurations: int = 1,
    seed: int | None = 0,
) -> MultisliceExitWave:
    """Propagate a plane wave through ``snapshot`` by the multislice algorithm.

    What it does
        Slices the specimen, and for each slice multiplies the wave by the
        band-limited transmission function :math:`\\exp(i\\sigma v_n)` and
        advances it by the slice thickness with the Fresnel propagator - the
        order abTEM uses, so the exit plane is the bottom of the box. The wave
        is stored at the slice boundaries nearest each requested depth, so a
        thickness series costs a single pass.

    When to use it
        For any HRTEM, exit-wave or electron-diffraction question about a
        specimen thicker than a weak phase object - a few ångströms of light
        atoms. The images come from `MultisliceExitWave.image`,
        `MultisliceExitWave.focal_series` and
        `MultisliceExitWave.defocus_thickness_map`.

    Parameters
    ----------
    snapshot:
        Specimen in an orthogonal box, beam along +z, periodic in x and y
        (`AtomicSnapshot.prepared_for_imaging` provides one).
    energy_kev:
        Beam energy.
    sampling_angstrom:
        Largest acceptable pixel; ignored when ``grid`` is given.
    grid:
        An explicit grid, e.g. to match another code point for point.
    slice_thickness_angstrom:
        Largest acceptable slice thickness.
    parametrization:
        Scattering-factor fit.
    tilt_mrad:
        Beam tilt ``(theta_x, theta_y)`` in mrad, entered in the propagator.
    exit_depths_angstrom:
        Depths at which to keep the wave; the bottom of the box if None.
    debye_waller_sigma_angstrom:
        Static RMS displacement damping the potential (time-averaged lattice).
    frozen_phonon_sigma_angstrom:
        RMS displacement per element (or one value) of a frozen-phonon ensemble.
    frozen_phonon_configurations:
        Number of configurations; 1 means a static lattice.
    seed:
        Seed of the displacement generator, so an ensemble is reproducible.

    Returns
    -------
    MultisliceExitWave

    Raises
    ------
    ValueError
        For a non-positive energy, sampling or slice thickness, an element
        without scattering factors, or a grid beyond `MAX_GRID_POINTS`.
    """

    if not energy_kev > 0.0:
        raise ValueError(f"The beam energy must be positive, got {energy_kev}.")
    if grid is None:
        if sampling_angstrom is None:
            raise ValueError("Give either a sampling or a grid.")
        extent = (float(snapshot.cell[0, 0]), float(snapshot.cell[1, 1]))
        grid = MultisliceGrid.from_sampling(extent, sampling_angstrom)
    configurations = int(frozen_phonon_configurations)
    if configurations < 1:
        raise ValueError(f"At least one configuration is needed, got {configurations}.")
    thermal = frozen_phonon_sigma_angstrom is not None and configurations > 1
    if not thermal:
        configurations = 1

    base = SlicedPotential.from_snapshot(
        snapshot, grid, slice_thickness_angstrom, parametrization, debye_waller_sigma_angstrom
    )
    edges = base.slice_edges_angstrom
    if exit_depths_angstrom is None:
        stops = np.array([base.num_slices], dtype=np.int64)
    else:
        wanted = np.asarray(list(exit_depths_angstrom), dtype=np.float64)
        if wanted.size == 0 or np.any(wanted <= 0.0):
            raise ValueError("Exit depths must be positive.")
        stops = np.unique(
            np.clip(np.rint(np.interp(wanted, edges, np.arange(edges.size))), 1, base.num_slices)
        ).astype(np.int64)

    sigma = relativistic_interaction_parameter_inv_v_angstrom(energy_kev)
    aperture = antialias_aperture(grid)
    dz = float(base.slice_thicknesses_angstrom[0])
    propagator = fresnel_propagator(grid, energy_kev, dz, tilt_mrad)
    ny, nx = grid.shape
    waves = np.empty((configurations, stops.size, ny, nx), dtype=np.complex128)
    retained = np.zeros(stops.size, dtype=np.float64)
    rng = np.random.default_rng(seed)
    total_steps = configurations * base.num_slices
    step = 0
    unique = base.unique_slice_count

    for configuration in range(configurations):
        sliced = base
        if thermal:
            assert frozen_phonon_sigma_angstrom is not None
            sliced = SlicedPotential.from_snapshot(
                _thermal_snapshot(snapshot, frozen_phonon_sigma_angstrom, rng),
                grid,
                slice_thickness_angstrom,
                parametrization,
                debye_waller_sigma_angstrom,
            )
            unique = max(unique, sliced.unique_slice_count)
        repeats = Counter(sliced.slice_keys)
        cache: dict[int, np.ndarray] = {}
        psi = np.ones((ny, nx), dtype=np.complex128)
        stored = 0
        for n in range(sliced.num_slices):
            key = sliced.slice_keys[n]
            t = cache.get(key)
            if t is None:
                t = _transmission(sliced.projected_potential(n), sigma, aperture)
                if repeats[key] > 1:
                    cache[key] = t
            psi = scipy.fft.ifft2(
                scipy.fft.fft2(psi * t, workers=-1) * propagator, workers=-1
            )
            step += 1
            if step % max(1, total_steps // 50) == 0:
                report(step / total_steps, stage="Propagating the wave through the slices")
            if stored < stops.size and n + 1 == stops[stored]:
                waves[configuration, stored] = psi
                retained[stored] += float(np.mean(np.abs(psi) ** 2)) / configurations
                stored += 1

    return MultisliceExitWave(
        waves=waves,
        depths_angstrom=edges[stops],
        grid=grid,
        energy_kev=float(energy_kev),
        tilt_mrad=(float(tilt_mrad[0]), float(tilt_mrad[1])),
        slice_thickness_angstrom=dz,
        num_slices=base.num_slices,
        unique_slices=unique,
        parametrization=base.parametrization,
        retained_intensity=retained,
        frozen_phonon_sigma_angstrom=frozen_phonon_sigma_angstrom if thermal else None,
        snapshot=snapshot,
    )


def _lens_transfer(grid: MultisliceGrid, aberrations: MicroscopeAberrations) -> np.ndarray:
    gx, gy = grid.spatial_frequencies()
    g = np.hypot(gx, gy)
    theta = np.arctan2(np.broadcast_to(gy, g.shape), np.broadcast_to(gx, g.shape))
    chi = aberrations.wave_aberration(g, theta)
    envelope = (
        aberrations.aperture_mask(g)
        * aberrations.temporal_envelope(g)
        * aberrations.spatial_envelope(g)
    )
    return np.asarray(envelope * np.exp(-1j * chi))


def focal_integration_nodes(
    phase_rate: float, requested: int | None = None
) -> tuple[int, float]:
    """Gauss-Hermite node count for a defocus phase rate, and the rate it resolves.

    An :math:`n`-node Gauss-Hermite rule integrates :math:`e^{-x^2 + iax}`
    to better than :math:`10^{-8}` for :math:`a \\le 2\\sqrt{n - 16}` (checked
    numerically up to 512 nodes). Given the largest rate ``phase_rate`` a
    calculation needs, this returns the smallest such :math:`n` (at least 24,
    at most ``MAX_FOCAL_INTEGRATION_POINTS``), or ``requested`` if given, and
    the rate :math:`2\\sqrt{n - 16}` that count resolves.
    """

    if requested is not None:
        count = int(requested)
        if count < 1:
            raise ValueError(f"Focal integration needs at least one node, got {requested}.")
    else:
        count = max(24, math.ceil(phase_rate * phase_rate / 4.0) + 16)
        count = min(count, MAX_FOCAL_INTEGRATION_POINTS)
    return count, 2.0 * math.sqrt(max(count - 16, 0))


def hrtem_image(
    exit_waves: np.ndarray,
    grid: MultisliceGrid,
    aberrations: MicroscopeAberrations,
    *,
    temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
    focal_integration_points: int | None = None,
) -> np.ndarray:
    """Image intensity of exit wave(s) through the objective lens.

    What it does
        For each exit wave :math:`\\psi` (one per frozen-phonon configuration)
        forms :math:`|\\mathcal F^{-1}\\{\\mathcal F\\{\\psi\\}H\\}|^2` with

        .. math::

           H(\\mathbf g) = A(g)\\,E_s(g)\\,E_c(g)\\,e^{-i\\chi(\\mathbf g)}

        and averages over configurations. With
        ``TemporalCoherence.FOCAL_INTEGRATION`` :math:`E_c` is dropped and the
        focal spread :math:`\\Delta` is instead integrated exactly,

        .. math::

           I = \\sum_k \\frac{w_k}{\\sqrt\\pi}
               \\bigl|\\mathcal F^{-1}\\{\\mathcal F\\{\\psi\\}
               H_{\\Delta f + \\sqrt2\\,\\Delta\\,x_k}\\}\\bigr|^2,

        over Gauss-Hermite nodes :math:`x_k` and weights :math:`w_k`. The
        defocus phase of a pair of beams varies across the spread as
        :math:`e^{iax}` with :math:`a = \\sqrt2\\,\\pi\\lambda\\Delta
        |g^2 - g'^2|`, and an :math:`n`-node rule integrates that exactly only
        up to :math:`a \\approx 2\\sqrt{n - 16}`. The node count is therefore
        chosen from the highest frequency the aperture and the grid pass
        (`focal_integration_nodes`), up to ``MAX_FOCAL_INTEGRATION_POINTS``;
        frequencies beyond what that many nodes resolve - where Frank's
        envelope is below :math:`e^{-495}` - are excluded, which drops only the
        interference of two such beams.

    Parameters
    ----------
    exit_waves:
        ``(ny, nx)`` or ``(configurations, ny, nx)`` complex exit wave(s).
    grid:
        Their grid.
    aberrations:
        Objective lens and coherence.
    temporal_coherence:
        See `TemporalCoherence`.
    focal_integration_points:
        Gauss-Hermite nodes for focal integration; chosen automatically when
        None.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` intensity; a vacuum wave of amplitude 1 images to 1.
    """

    waves = np.asarray(exit_waves, dtype=np.complex128)
    if waves.ndim == 2:
        waves = waves[None]
    spectra = scipy.fft.fft2(waves, workers=-1)
    mode = TemporalCoherence(temporal_coherence)
    if mode is TemporalCoherence.QUASI_COHERENT or aberrations.focal_spread_angstrom <= 0.0:
        lenses = [(1.0, _lens_transfer(grid, aberrations))]
    else:
        spread = aberrations.focal_spread_angstrom
        g = grid.frequency_magnitude()
        passed = g[aberrations.aperture_mask(g) > 0.0]
        rate = math.sqrt(2.0) * math.pi * aberrations.wavelength_angstrom * spread
        count, resolved = focal_integration_nodes(
            rate * float(np.max(passed)) ** 2 if passed.size else 0.0, focal_integration_points
        )
        nodes, weights = scipy.special.roots_hermite(count)
        coherent_band = rate * g * g <= resolved
        lenses = [
            (
                float(w) / math.sqrt(math.pi),
                coherent_band
                * _lens_transfer(
                    grid,
                    replace(
                        aberrations,
                        defocus_angstrom=aberrations.defocus_angstrom
                        + math.sqrt(2.0) * spread * float(x),
                        focal_spread_angstrom=0.0,
                    ),
                ),
            )
            for x, w in zip(nodes, weights, strict=True)
        ]
    image = np.zeros(grid.shape, dtype=np.float64)
    for weight, lens in lenses:
        image_waves = scipy.fft.ifft2(spectra * lens, workers=-1)
        image += weight * np.mean(np.abs(image_waves) ** 2, axis=0)
    return image


def _contrast(image: np.ndarray) -> float:
    mean = float(np.mean(image))
    return float(np.std(image) / mean) if mean > 0.0 else 0.0


@dataclass(frozen=True, slots=True)
class FocalSeries:
    """HRTEM images of one exit wave at a sequence of defoci.

    Attributes
    ----------
    defoci_angstrom:
        Defocus of each image (PyTex convention: negative is underfocus).
    images:
        ``(defoci, ny, nx)`` intensities.
    thickness_angstrom:
        Specimen thickness at which the exit wave was taken.
    aberrations:
        The lens, apart from the defocus that varies.
    temporal_coherence:
        How focal spread entered.
    grid:
        Lateral grid.
    label:
        Specimen label.
    """

    defoci_angstrom: np.ndarray
    images: np.ndarray
    thickness_angstrom: float
    aberrations: MicroscopeAberrations
    temporal_coherence: TemporalCoherence
    grid: MultisliceGrid
    label: str

    @property
    def contrasts(self) -> np.ndarray:
        """RMS contrast :math:`\\sigma_I/\\langle I\\rangle` of each image."""
        return np.array([_contrast(image) for image in self.images])

    @property
    def mean_intensities(self) -> np.ndarray:
        """Mean intensity of each image (1 when nothing is lost to the aperture)."""
        return np.asarray(np.mean(self.images, axis=(1, 2)))

    @property
    def minimum_contrast_defocus_angstrom(self) -> float:
        """The defocus of least contrast - near Gaussian focus for a weak phase object."""
        return float(self.defoci_angstrom[int(np.argmin(self.contrasts))])

    @property
    def maximum_contrast_defocus_angstrom(self) -> float:
        """The defocus of greatest RMS contrast in the series."""
        return float(self.defoci_angstrom[int(np.argmax(self.contrasts))])

    def describe(self) -> str:
        """Convention-explicit account of the series."""
        c = self.contrasts
        step = (
            f" in steps of {float(np.mean(np.diff(self.defoci_angstrom))):.1f} Å"
            if self.defoci_angstrom.size > 1
            else ""
        )
        return (
            f"Focal series of {self.label} at {self.thickness_angstrom:.1f} Å thickness: "
            f"{self.defoci_angstrom.size} images from Δf = {self.defoci_angstrom[0]:.1f} to "
            f"{self.defoci_angstrom[-1]:.1f} Å{step} (negative is underfocus), "
            f"{self.aberrations.energy_kev:.0f} kV, Cs = {self.aberrations.cs_um:.1f} µm, "
            f"{self.temporal_coherence.value.replace('_', ' ')} temporal coherence. "
            f"RMS contrast runs from {float(c.min()) * 100:.1f} % (Δf = "
            f"{self.minimum_contrast_defocus_angstrom:.1f} Å) to {float(c.max()) * 100:.1f} % "
            f"(Δf = {self.maximum_contrast_defocus_angstrom:.1f} Å). All images share one exit "
            "wave, so they differ only by the objective lens (Coene et al. 1992)."
        )


@dataclass(frozen=True, slots=True)
class DefocusThicknessMap:
    """HRTEM images over a grid of thicknesses (rows) and defoci (columns).

    Attributes
    ----------
    defoci_angstrom, thicknesses_angstrom:
        The two axes.
    images:
        ``(thicknesses, defoci, ny, nx)`` intensities.
    aberrations, temporal_coherence, grid, label:
        As in `FocalSeries`.
    """

    defoci_angstrom: np.ndarray
    thicknesses_angstrom: np.ndarray
    images: np.ndarray
    aberrations: MicroscopeAberrations
    temporal_coherence: TemporalCoherence
    grid: MultisliceGrid
    label: str

    @property
    def contrasts(self) -> np.ndarray:
        """``(thicknesses, defoci)`` RMS contrast."""
        return np.array([[_contrast(image) for image in row] for row in self.images])

    def describe(self) -> str:
        """Convention-explicit account of the map."""
        c = self.contrasts
        i, j = np.unravel_index(int(np.argmax(c)), c.shape)
        return (
            f"Defocus-thickness map of {self.label}: {self.thicknesses_angstrom.size} "
            f"thicknesses from {self.thicknesses_angstrom[0]:.1f} to "
            f"{self.thicknesses_angstrom[-1]:.1f} Å by {self.defoci_angstrom.size} defoci from "
            f"{self.defoci_angstrom[0]:.1f} to {self.defoci_angstrom[-1]:.1f} Å, "
            f"{self.aberrations.energy_kev:.0f} kV, Cs = {self.aberrations.cs_um:.1f} µm. "
            f"Greatest RMS contrast {float(c[i, j]) * 100:.1f} % at t = "
            f"{self.thicknesses_angstrom[i]:.1f} Å, Δf = {self.defoci_angstrom[j]:.1f} Å. "
            "Matching an experimental image against such a tableau is how thickness and "
            "defocus are determined in practice (O'Keefe & Kilaas 1988)."
        )


def simulate_multislice_hrem(
    snapshot: AtomicSnapshot,
    aberrations: MicroscopeAberrations,
    sampling_angstrom: float = 0.1,
    slice_thickness_angstrom: float = 1.0,
    *,
    parametrization: PotentialParametrization | str = PotentialParametrization.LOBATO,
    tilt_mrad: tuple[float, float] = (0.0, 0.0),
    temporal_coherence: TemporalCoherence | str = TemporalCoherence.QUASI_COHERENT,
    frozen_phonon_sigma_angstrom: Mapping[str, float] | float | None = None,
    frozen_phonon_configurations: int = 1,
    seed: int | None = 0,
) -> HREMSimulationResult:
    """One HRTEM image by PyTex multislice, as an `HREMSimulationResult`.

    The engine behind ``simulate_hrem(..., engine="multislice")`` and the
    workbench: `multislice` for the exit wave, then `hrtem_image` through
    ``aberrations``, then the power spectrum of the image.
    """

    exit_wave = multislice(
        snapshot,
        aberrations.energy_kev,
        sampling_angstrom=sampling_angstrom,
        slice_thickness_angstrom=slice_thickness_angstrom,
        parametrization=parametrization,
        tilt_mrad=tilt_mrad,
        frozen_phonon_sigma_angstrom=frozen_phonon_sigma_angstrom,
        frozen_phonon_configurations=frozen_phonon_configurations,
        seed=seed,
    )
    report(1.0, stage="Applying the objective lens transfer function")
    return exit_wave.hrem_result(aberrations, temporal_coherence=temporal_coherence)


def multislice_summary(exit_wave: MultisliceExitWave) -> dict[str, Any]:
    """Machine-readable diagnostics of a run, for reports and the workbench."""
    ny, nx = exit_wave.grid.shape
    return {
        "grid_px": [nx, ny],
        "sampling_angstrom": list(exit_wave.grid.sampling_angstrom),
        "slice_thickness_angstrom": exit_wave.slice_thickness_angstrom,
        "num_slices": exit_wave.num_slices,
        "unique_slices": exit_wave.unique_slices,
        "configurations": exit_wave.configurations,
        "parametrization": exit_wave.parametrization.value,
        "band_limit_inv_angstrom": exit_wave.grid.antialias_cutoff_inv_angstrom,
        "max_scattering_angle_mrad": exit_wave.grid.max_scattering_angle_mrad(
            exit_wave.energy_kev
        ),
        "depths_angstrom": [float(t) for t in exit_wave.depths_angstrom],
        "retained_intensity": [float(v) for v in exit_wave.retained_intensity],
    }
