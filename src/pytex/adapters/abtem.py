"""Adapter layer connecting PyTex canonical HREM data models to abTEM and ASE.

This module provides bidirectional conversion between PyTex crystallographic snapshots
and ASE/abTEM representations, enabling state-of-the-art multislice High-Resolution
Transmission Electron Microscopy (HRTEM) simulations for double-corrected microscopes.

Scientific References
---------------------
- Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, 2nd ed., Springer.
- Madsen, J. et al. (2021). abTEM: An open-source framework for simulation of
  transmission electron microscopy. ChemPhysChem 22, 1-13.
"""

from __future__ import annotations

import importlib.util
import math
from typing import TYPE_CHECKING, Any

import numpy as np

from pytex.diffraction.hrem import (
    AtomicSnapshot,
    HREMSimulationResult,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
)

if TYPE_CHECKING:
    import abtem
    import ase


def is_abtem_available() -> bool:
    """Check whether both abTEM and ASE packages are importable.

    Returns
    -------
    bool
        True if both abtem and ase can be imported, False otherwise.
    """
    return (
        importlib.util.find_spec("abtem") is not None
        and importlib.util.find_spec("ase") is not None
    )


def to_ase_atoms(snapshot: AtomicSnapshot) -> ase.Atoms:
    """Convert a canonical PyTex AtomicSnapshot into an ASE Atoms object.

    Parameters
    ----------
    snapshot : AtomicSnapshot
        The PyTex atomic structure snapshot.

    Returns
    -------
    ase.Atoms
        The corresponding ASE Atoms structure with cell, positions, and boundary conditions.
    """
    return snapshot.to_ase()


def from_ase_atoms(atoms: ase.Atoms, label: str = "Imported from ASE") -> AtomicSnapshot:
    """Create a canonical PyTex AtomicSnapshot from an ASE Atoms object.

    Parameters
    ----------
    atoms : ase.Atoms
        The ASE Atoms instance.
    label : str, default="Imported from ASE"
        Descriptive label for the snapshot.

    Returns
    -------
    AtomicSnapshot
        The canonical PyTex snapshot.
    """
    return AtomicSnapshot.from_ase(atoms, label=label)


def to_abtem_ctf(aberrations: MicroscopeAberrations) -> abtem.transfer.CTF:
    """Construct an abTEM Contrast Transfer Function (CTF) from canonical PyTex aberrations.

    Translates electron kinetic energy (eV), defocus (Å), spherical aberration (Å),
    higher-order aberrations, focal spread (Å, 1/e width), angular spread (mrad),
    and objective aperture cutoff (mrad) into abTEM's optical parameter convention.

    Parameters
    ----------
    aberrations : MicroscopeAberrations
        Canonical PyTex microscope aberration parameters.

    Returns
    -------
    abtem.transfer.CTF
        An abTEM CTF filter ready to be applied to multislice exit waves.

    Raises
    ------
    ImportError
        If abTEM is not installed in the active environment.
    """
    if not is_abtem_available():
        raise ImportError(
            "abTEM is not installed. Install it with `pip install abtem ase` to use CTF conversion."
        )

    import abtem.transfer

    energy_ev = float(aberrations.energy_kev * 1000.0)
    defocus_angstrom = float(aberrations.defocus_angstrom)
    cs_angstrom = float(aberrations.cs_angstrom)
    c5_angstrom = float(aberrations.c5_mm * 1e7)
    c12_angstrom = float(aberrations.astigmatism_angstrom)
    phi12_rad = float(math.radians(aberrations.astigmatism_angle_deg))

    # Objective aperture semiangle cutoff in mrad (abTEM uses inf for open aperture)
    if aberrations.aperture_cutoff_mrad is not None and aberrations.aperture_cutoff_mrad > 0.0:
        semiangle_cutoff = float(aberrations.aperture_cutoff_mrad)
    else:
        semiangle_cutoff = float("inf")

    # Temporal coherence: 1/e focal spread in Å
    focal_spread = float(aberrations.focal_spread_angstrom)

    # Spatial coherence: angular spread in mrad
    angular_spread = float(aberrations.convergence_semiangle_mrad)

    kwargs: dict[str, Any] = {
        "energy": energy_ev,
        "semiangle_cutoff": semiangle_cutoff,
        "focal_spread": focal_spread,
        "angular_spread": angular_spread,
        "C10": defocus_angstrom,
        "C30": cs_angstrom,
    }

    if abs(c5_angstrom) > 1e-6:
        kwargs["C50"] = c5_angstrom
    if abs(c12_angstrom) > 1e-6:
        kwargs["C12"] = c12_angstrom
        kwargs["phi12"] = phi12_rad

    return abtem.transfer.CTF(**kwargs)


def simulate_hrem_multislice(
    snapshot: AtomicSnapshot,
    aberrations: MicroscopeAberrations,
    sampling_angstrom: float = 0.1,
    slice_thickness_angstrom: float = 1.0,
) -> HREMSimulationResult:
    """Simulate High-Resolution TEM using abTEM multislice wave propagation.

    Computes the electrostatic projected potential slices from atomic positions,
    propagates a relativistic plane wave through the specimen via the multislice
    algorithm, filters the complex exit wave by the objective lens contrast transfer
    function (CTF), and measures real-space intensity.

    Parameters
    ----------
    snapshot : AtomicSnapshot
        The atomic structure snapshot to simulate.
    aberrations : MicroscopeAberrations
        Microscope optical parameters and aberration coefficients.
    sampling_angstrom : float, default=0.1
        Real-space pixel sampling in Angstrom.
    slice_thickness_angstrom : float, default=1.0
        Multislice z-slice thickness in Angstrom.

    Returns
    -------
    HREMSimulationResult
        The simulated micrograph, complex exit wave, power spectrum (Thon rings),
        and explainable diagnostics.

    Raises
    ------
    ImportError
        If abTEM or ASE is not installed.
    """
    if not is_abtem_available():
        raise ImportError(
            "abTEM and ASE are required for multislice simulation. "
            "Install them via `pip install abtem ase` or use "
            "`simulate_hrem(..., prefer_abtem=False)`."
        )

    import abtem

    atoms = snapshot.to_ase()

    # Ensure cell has finite thickness along z for multislice slicing
    if atoms.cell[2, 2] <= 0.0 or math.isclose(atoms.cell[2, 2], 0.0):
        cell = atoms.cell.copy()
        cell[2, 2] = max(snapshot.thickness_angstrom, 10.0)
        atoms.set_cell(cell)

    potential = abtem.Potential(
        atoms,
        sampling=sampling_angstrom,
        slice_thickness=slice_thickness_angstrom,
    )

    energy_ev = float(aberrations.energy_kev * 1000.0)
    plane_wave = abtem.PlaneWave(energy=energy_ev)

    exit_wave_calc = plane_wave.multislice(potential)
    ctf_filter = to_abtem_ctf(aberrations)
    image_wave = exit_wave_calc.apply_ctf(ctf_filter)

    # Compute arrays
    intensity_raw = image_wave.intensity().compute().array
    exit_wave_raw = exit_wave_calc.compute().array

    # Squeeze potential batch/ensemble dimensions
    intensity_2d = np.squeeze(np.asarray(intensity_raw, dtype=np.float64))
    exit_wave_2d = np.squeeze(np.asarray(exit_wave_raw, dtype=np.complex128))

    if intensity_2d.ndim != 2:
        intensity_2d = intensity_2d.reshape((intensity_2d.shape[-2], intensity_2d.shape[-1]))
    if exit_wave_2d.ndim != 2:
        exit_wave_2d = exit_wave_2d.reshape((exit_wave_2d.shape[-2], exit_wave_2d.shape[-1]))

    ny, nx = intensity_2d.shape
    extent_x = float(atoms.cell[0, 0]) if atoms.cell[0, 0] > 0 else float(nx * sampling_angstrom)
    extent_y = float(atoms.cell[1, 1]) if atoms.cell[1, 1] > 0 else float(ny * sampling_angstrom)
    pixel_size = float(extent_x / nx)

    # 2D power spectrum (log-scaled Thon rings)
    diff = intensity_2d - np.mean(intensity_2d)
    fft_val = np.fft.fftshift(np.fft.fft2(diff))
    power_spectrum = np.log10(1.0 + np.abs(fft_val) ** 2)

    ctf_1d = aberrations.evaluate_ctf_1d(max_q_inv_angstrom=1.0 / (2.0 * pixel_size))

    return HREMSimulationResult(
        image=intensity_2d,
        exit_wave=exit_wave_2d,
        pixel_size_angstrom=pixel_size,
        extent_angstrom=(extent_x, extent_y),
        power_spectrum=power_spectrum,
        ctf=ctf_1d,
        aberrations=aberrations,
        snapshot=snapshot,
    )


def simulate_hrem(
    snapshot: AtomicSnapshot,
    aberrations: MicroscopeAberrations,
    sampling_angstrom: float = 0.1,
    slice_thickness_angstrom: float = 1.0,
    prefer_abtem: bool = True,
) -> HREMSimulationResult:
    """Unified HREM simulation entry point with automatic multislice / pure-Python dispatch.

    If `prefer_abtem` is True and the `abtem` package is installed, the simulation
    employs full multislice propagation through atomic potential slices. Otherwise,
    it falls back cleanly to the pure-Python phase-object transmission simulation.

    Parameters
    ----------
    snapshot : AtomicSnapshot
        The atomic structure snapshot.
    aberrations : MicroscopeAberrations
        Microscope optical parameters.
    sampling_angstrom : float, default=0.1
        Lateral real-space pixel sampling in Angstrom.
    slice_thickness_angstrom : float, default=1.0
        Multislice z-slice thickness in Angstrom (used if multislice is enabled).
    prefer_abtem : bool, default=True
        Whether to prioritize abTEM multislice if abTEM is installed.

    Returns
    -------
    HREMSimulationResult
        The simulated high-resolution micrograph and associated metadata.
    """
    if prefer_abtem and is_abtem_available():
        return simulate_hrem_multislice(
            snapshot=snapshot,
            aberrations=aberrations,
            sampling_angstrom=sampling_angstrom,
            slice_thickness_angstrom=slice_thickness_angstrom,
        )

    return pure_python_phase_object_simulation(
        snapshot=snapshot,
        aberrations=aberrations,
        sampling_angstrom=sampling_angstrom,
    )
