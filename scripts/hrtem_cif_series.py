"""Generate focal and thickness HRTEM series from an MD snapshot stored as CIF.

Run from a PyTex checkout (or an environment where PyTex is installed):

    python scripts/hrtem_cif_series.py snapshot.cif

The CIF cell must be orthogonal.  ``--beam-axis`` chooses which CIF cell vector
is parallel to the electron beam (+z).  Intensities are saved losslessly in an
NPZ file; PNG files use one common display range so contrast remains comparable.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
import numpy as np
from pymatgen.core import Structure

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pytex.diffraction.hrem import AtomicSnapshot, DoubleCorrectionMode, MicroscopeAberrations
from pytex.diffraction.multislice import TemporalCoherence, multislice


@dataclass(frozen=True)
class RunParameters:
    """Inputs recorded beside every generated series."""

    cif: str
    beam_axis: str
    energy_kev: float
    sampling_angstrom: float
    slice_thickness_angstrom: float
    defocus_start_angstrom: float
    defocus_stop_angstrom: float
    defocus_step_angstrom: float
    thickness_count: int
    cs_um: float
    c5_mm: float
    focal_spread_angstrom: float
    convergence_semiangle_mrad: float
    aperture_cutoff_mrad: float
    temporal_coherence: str


def load_cif_snapshot(path: Path, beam_axis: str) -> AtomicSnapshot:
    """Load an ordered orthogonal CIF snapshot and put ``beam_axis`` along +z."""

    structure = Structure.from_file(path)
    if any(not site.is_ordered for site in structure):
        raise ValueError("The CIF contains disordered/partially occupied sites.")
    angles = np.asarray(structure.lattice.angles, dtype=np.float64)
    if not np.allclose(angles, 90.0, atol=1.0e-4):
        raise ValueError(
            "The CIF cell is not orthogonal. Export an orthogonal MD supercell before HRTEM "
            f"simulation (cell angles are {angles.tolist()} degrees)."
        )

    # Fractional coordinates make this independent of pymatgen's Cartesian
    # representation of the orthogonal cell.  The selected cell vector becomes z.
    permutations = {"a": (1, 2, 0), "b": (2, 0, 1), "c": (0, 1, 2)}
    order = permutations[beam_axis]
    lengths = np.asarray(structure.lattice.abc, dtype=np.float64)[list(order)]
    fractional = np.mod(np.asarray(structure.frac_coords, dtype=np.float64), 1.0)
    positions = fractional[:, list(order)] * lengths
    species = tuple(site.specie.symbol for site in structure)
    return AtomicSnapshot(
        species=species,
        positions=positions,
        cell=np.diag(lengths),
        periodicity=(True, True, False),
        label=f"{path.name} (beam along CIF {beam_axis})",
    ).prepared_for_imaging(periodic_xy=True)


def inclusive_range(start: float, stop: float, step: float) -> np.ndarray:
    """Return an inclusive increasing range without silently overshooting."""

    if step <= 0.0 or stop < start:
        raise ValueError("Range step must be positive and stop must be at least start.")
    return np.arange(start, stop + 0.5 * step, step, dtype=np.float64)


def save_image(path: Path, image: np.ndarray, vmin: float, vmax: float) -> None:
    """Save one image without axes, using the common series display limits."""

    plt.imsave(path, image, cmap="gray", vmin=vmin, vmax=vmax, origin="lower")


def save_tableau(
    path: Path,
    images: np.ndarray,
    thicknesses: np.ndarray,
    defoci: np.ndarray,
    extent: tuple[float, float],
    vmin: float,
    vmax: float,
) -> None:
    """Save a labelled thickness-by-defocus overview."""

    rows, cols = images.shape[:2]
    fig, axes = plt.subplots(rows, cols, figsize=(2.3 * cols, 2.3 * rows), squeeze=False)
    for i in range(rows):
        for j in range(cols):
            axes[i, j].imshow(
                images[i, j],
                cmap="gray",
                vmin=vmin,
                vmax=vmax,
                origin="lower",
                extent=(0.0, extent[0], 0.0, extent[1]),
            )
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            if i == 0:
                axes[i, j].set_title(f"Δf = {defoci[j]:.0f} Å", fontsize=9)
            if j == 0:
                axes[i, j].set_ylabel(f"t = {thicknesses[i]:.1f} Å", fontsize=9)
    fig.suptitle("HRTEM defocus–thickness tableau (common intensity scale)")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    """Run one multislice calculation and write both requested series."""

    cif = args.cif.resolve()
    output = (args.output or Path("outputs") / f"hrtem_{cif.stem}").resolve()
    output.mkdir(parents=True, exist_ok=True)
    focal_dir = output / "focal_series"
    thickness_dir = output / "thickness_series"
    focal_dir.mkdir(exist_ok=True)
    thickness_dir.mkdir(exist_ok=True)

    snapshot = load_cif_snapshot(cif, args.beam_axis)
    defoci = inclusive_range(args.defocus_start, args.defocus_stop, args.defocus_step)
    if args.thickness_count < 1:
        raise ValueError("--thickness-count must be at least 1.")
    if defoci.size * args.thickness_count > 256:
        raise ValueError("The series is limited to 256 thickness-defocus image pairs.")
    maximum_thickness = float(snapshot.cell[2, 2])
    first_depth = min(
        maximum_thickness,
        max(args.slice_thickness, maximum_thickness / args.thickness_count),
    )
    requested_depths = np.linspace(first_depth, maximum_thickness, args.thickness_count)

    lens = MicroscopeAberrations(
        energy_kev=args.energy,
        defocus_angstrom=0.0,
        cs_mm=args.cs_um * 1.0e-3,
        c5_mm=args.c5_mm,
        focal_spread_angstrom=args.focal_spread,
        convergence_semiangle_mrad=args.convergence,
        aperture_cutoff_mrad=args.aperture if args.aperture > 0.0 else None,
        mode=DoubleCorrectionMode.DOUBLE_CORRECTED,
    )
    coherence = TemporalCoherence(args.temporal_coherence)
    exit_wave = multislice(
        snapshot,
        args.energy,
        sampling_angstrom=args.sampling,
        slice_thickness_angstrom=args.slice_thickness,
        exit_depths_angstrom=requested_depths,
    )
    tableau = exit_wave.defocus_thickness_map(lens, defoci, coherence)
    images = tableau.images
    thicknesses = tableau.thicknesses_angstrom
    contrasts = tableau.contrasts

    # One robust common display range preserves relative contrast across all PNGs.
    vmin, vmax = (float(x) for x in np.percentile(images, [0.5, 99.5]))
    if vmax <= vmin:
        vmin, vmax = float(images.min()), float(images.max())

    for j, defocus in enumerate(defoci):
        save_image(
            focal_dir / f"defocus_{defocus:+08.1f}A.png",
            images[-1, j],
            vmin,
            vmax,
        )
    focus_index = int(np.argmin(np.abs(defoci - args.thickness_defocus)))
    for i, thickness in enumerate(thicknesses):
        save_image(
            thickness_dir / f"thickness_{thickness:08.2f}A.png",
            images[i, focus_index],
            vmin,
            vmax,
        )

    np.savez_compressed(
        output / "hrtem_series.npz",
        images=images,
        defoci_angstrom=defoci,
        thicknesses_angstrom=thicknesses,
        contrasts=contrasts,
        sampling_angstrom=np.asarray(exit_wave.grid.sampling_angstrom),
        extent_angstrom=np.asarray(exit_wave.grid.extent_angstrom),
    )
    with (output / "series.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("thickness_angstrom", "defocus_angstrom", "rms_contrast"))
        for i, thickness in enumerate(thicknesses):
            for j, defocus in enumerate(defoci):
                writer.writerow((f"{thickness:.6g}", f"{defocus:.6g}", f"{contrasts[i, j]:.8g}"))

    params = RunParameters(
        cif=str(cif),
        beam_axis=args.beam_axis,
        energy_kev=args.energy,
        sampling_angstrom=args.sampling,
        slice_thickness_angstrom=args.slice_thickness,
        defocus_start_angstrom=args.defocus_start,
        defocus_stop_angstrom=args.defocus_stop,
        defocus_step_angstrom=args.defocus_step,
        thickness_count=args.thickness_count,
        cs_um=args.cs_um,
        c5_mm=args.c5_mm,
        focal_spread_angstrom=args.focal_spread,
        convergence_semiangle_mrad=args.convergence,
        aperture_cutoff_mrad=args.aperture,
        temporal_coherence=args.temporal_coherence,
    )
    (output / "run_parameters.json").write_text(
        json.dumps(asdict(params), indent=2) + "\n", encoding="utf-8"
    )
    (output / "simulation_report.txt").write_text(
        snapshot.describe() + "\n\n" + exit_wave.describe() + "\n\n" + tableau.describe() + "\n",
        encoding="utf-8",
    )
    save_tableau(
        output / "defocus_thickness_tableau.png",
        images,
        thicknesses,
        defoci,
        exit_wave.grid.extent_angstrom,
        vmin,
        vmax,
    )
    print(f"Saved HRTEM series to {output}")
    return output


def parser() -> argparse.ArgumentParser:
    """Command-line interface, with conservative double-corrected defaults."""

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cif", type=Path, help="CIF containing one ordered MD atomic snapshot")
    p.add_argument("--output", type=Path, help="Output directory (default: outputs/hrtem_<name>)")
    p.add_argument("--beam-axis", choices=("a", "b", "c"), default="c")
    p.add_argument("--energy", type=float, default=300.0, help="Beam energy in keV")
    p.add_argument("--sampling", type=float, default=0.10, help="Maximum image pixel size in Å")
    p.add_argument(
        "--slice-thickness", type=float, default=1.0, help="Maximum slice thickness in Å"
    )
    p.add_argument("--defocus-start", type=float, default=-100.0, help="First defocus in Å")
    p.add_argument("--defocus-stop", type=float, default=100.0, help="Last defocus in Å")
    p.add_argument("--defocus-step", type=float, default=25.0, help="Defocus step in Å")
    p.add_argument(
        "--thickness-count", type=int, default=5, help="Stored thicknesses up to full cell"
    )
    p.add_argument(
        "--thickness-defocus", type=float, default=0.0, help="Defocus used for thickness PNGs"
    )
    p.add_argument("--cs-um", type=float, default=0.0, help="Residual third-order Cs in µm")
    p.add_argument("--c5-mm", type=float, default=0.0, help="Residual fifth-order C5 in mm")
    p.add_argument("--focal-spread", type=float, default=5.0, help="Gaussian focal-spread σ in Å")
    p.add_argument("--convergence", type=float, default=0.2, help="Convergence semi-angle in mrad")
    p.add_argument(
        "--aperture",
        type=float,
        default=30.0,
        help="Objective aperture semi-angle in mrad; 0 disables",
    )
    p.add_argument(
        "--temporal-coherence",
        choices=("quasi_coherent", "focal_integration"),
        default="quasi_coherent",
    )
    return p


if __name__ == "__main__":
    run(parser().parse_args())
