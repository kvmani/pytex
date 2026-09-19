"""PyTex multislice against abTEM, number for number, when abTEM is installed.

abTEM is optional, so these tests skip without it; the analytic and Bloch-wave
validation in ``test_multislice.py`` does not need it. What agreement to expect:

- the parametrization tables are identical (PyTex's is generated from abTEM's);
- with atoms on grid points the exit waves agree to abTEM's single precision;
- with atoms off the grid they agree to about a percent, because abTEM places
  each atom by bilinear delta interpolation followed by an approximate sinc
  correction, while PyTex synthesizes the exact phase factor
  :math:`e^{-2\\pi i\\mathbf g\\cdot\\mathbf r}`;
- the objective lens (including the sqrt(2) focal-spread conversion of the
  adapter) reproduces abTEM's quasi-coherent image.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

abtem = pytest.importorskip("abtem")

from pytex.adapters.abtem import to_abtem_ctf  # noqa: E402
from pytex.app.phases import phase_from_request  # noqa: E402
from pytex.diffraction.hrem import (  # noqa: E402
    MicroscopeAberrations,
    relativistic_interaction_parameter_inv_v_angstrom,
    relativistic_wavelength_angstrom,
)
from pytex.diffraction.multislice import (  # noqa: E402
    MultisliceGrid,
    SlicedPotential,
    _parametrization_table,
    multislice,
    periodic_slab,
)


def _abtem_exit_wave(snapshot, grid, slice_thickness):  # type: ignore[no-untyped-def]
    ny, nx = grid.shape
    potential = abtem.Potential(
        snapshot.to_ase(), gpts=(nx, ny), slice_thickness=slice_thickness, projection="infinite"
    )
    wave = abtem.PlaneWave(energy=200e3).multislice(potential)
    return wave, np.asarray(wave.compute().array).T  # abTEM indexes [x, y]


def test_parametrization_tables_are_abtems() -> None:
    source = Path(abtem.__file__).resolve().parent / "parametrizations" / "data"
    for name in ("lobato", "kirkland"):
        theirs = json.loads((source / f"{name}.json").read_text(encoding="utf-8"))
        ours = _parametrization_table()[name]
        assert set(ours) == set(theirs)
        for symbol, rows in theirs.items():
            np.testing.assert_array_equal(ours[symbol], np.asarray(rows, dtype=np.float64))


def test_wavelength_and_interaction_parameter_match() -> None:
    from abtem.core.energy import energy2sigma, energy2wavelength

    for kev in (80.0, 200.0, 300.0):
        assert relativistic_wavelength_angstrom(kev) == pytest.approx(
            energy2wavelength(kev * 1e3), rel=1e-7
        )
        assert relativistic_interaction_parameter_inv_v_angstrom(kev) == pytest.approx(
            energy2sigma(kev * 1e3), rel=1e-7
        )


def test_projected_potential_mean_matches() -> None:
    phase = phase_from_request({"builtin": "si_diamond"})[1]
    snap, _, _ = periodic_slab(phase, (1, 1, 0), (1, 1), beam_repeats=1)
    grid = MultisliceGrid((96, 64), (snap.cell[0, 0], snap.cell[1, 1]))
    ny, nx = grid.shape
    theirs = abtem.Potential(
        snap.to_ase(), gpts=(nx, ny), slice_thickness=snap.cell[2, 2], projection="infinite"
    ).build().compute().array
    ours = SlicedPotential.from_snapshot(snap, grid, snap.cell[2, 2]).projected_potential(0)
    assert float(np.mean(ours)) == pytest.approx(float(np.mean(theirs)), rel=1e-5)


def test_exit_wave_matches_with_atoms_on_the_grid() -> None:
    phase = phase_from_request({"builtin": "si_diamond"})[1]
    a = phase.lattice.a
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (2, 2), beam_repeats=10)
    grid = MultisliceGrid((128, 128), (2 * a, 2 * a))  # every atom 16 pixels apart
    ours = multislice(snap, 200.0, grid=grid, slice_thickness_angstrom=a / 4)
    _, theirs = _abtem_exit_wave(snap, grid, a / 4)
    error = np.linalg.norm(ours.wave() - theirs) / np.linalg.norm(theirs)
    assert error < 2e-3


def test_exit_wave_matches_with_atoms_off_the_grid() -> None:
    phase = phase_from_request({"builtin": "si_diamond"})[1]
    snap, _, _ = periodic_slab(phase, (1, 1, 0), (2, 1), thickness_angstrom=40.0)
    grid = MultisliceGrid.from_sampling((snap.cell[0, 0], snap.cell[1, 1]), 0.08)
    ours = multislice(snap, 200.0, grid=grid, slice_thickness_angstrom=1.0)
    _, theirs = _abtem_exit_wave(snap, grid, 1.0)
    error = np.linalg.norm(ours.wave() - theirs) / np.linalg.norm(theirs)
    assert error < 2e-2


def test_image_through_the_objective_lens_matches() -> None:
    phase = phase_from_request({"builtin": "si_diamond"})[1]
    a = phase.lattice.a
    snap, _, _ = periodic_slab(phase, (0, 0, 1), (2, 2), beam_repeats=6)
    grid = MultisliceGrid((128, 128), (2 * a, 2 * a))
    lens = MicroscopeAberrations(
        energy_kev=200.0, defocus_angstrom=-60.0, cs_mm=0.02, focal_spread_angstrom=25.0,
        convergence_semiangle_mrad=0.3, aperture_cutoff_mrad=30.0,
    )
    ours = multislice(snap, 200.0, grid=grid, slice_thickness_angstrom=a / 4).image(lens)
    wave, _ = _abtem_exit_wave(snap, grid, a / 4)
    theirs = np.asarray(wave.apply_ctf(to_abtem_ctf(lens)).intensity().compute().array).T
    # abTEM's hard aperture edge and single precision differ slightly at the cutoff.
    np.testing.assert_allclose(ours, theirs, atol=5e-3 * float(np.ptp(theirs)))
