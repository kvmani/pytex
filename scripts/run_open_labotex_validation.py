from __future__ import annotations

import json
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    KernelSpec,
    Lattice,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    SymmetrySpec,
    invert_labotex_pole_figures,
    load_labotex_pole_figures,
    plot_odf,
    plot_pole_figure,
)

EXERCISES_ZIP_URL = "https://labosoft.com.pl/download/Exercises.zip"
OUTPUT_ROOT = Path("inspection_outputs/open_labotex_validation")
DATASET_ROOT = OUTPUT_ROOT / "datasets"
PLOTS_ROOT = OUTPUT_ROOT / "plots"

DATASETS = {
    "cu_fcc": {
        "filename": "Exercise2-Cu.PPF",
        "phase_name": "copper_demo",
        "lattice": (1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        "point_group": "m-3m",
        "dictionary_steps": (15.0, 10.0, 15.0),
        "kernel_halfwidth_deg": 8.0,
    },
    "fe_bcc": {
        "filename": "Exercise4-Fe.epf",
        "phase_name": "ferrite_demo",
        "lattice": (1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        "point_group": "m-3m",
        "dictionary_steps": (20.0, 15.0, 20.0),
        "kernel_halfwidth_deg": 10.0,
    },
    "al_fcc": {
        "filename": "Exercise5-Al.epf",
        "phase_name": "aluminium_demo",
        "lattice": (1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        "point_group": "m-3m",
        "dictionary_steps": (20.0, 15.0, 20.0),
        "kernel_halfwidth_deg": 10.0,
    },
}


def _phase_from_spec(spec: dict[str, object], *, crystal_frame: ReferenceFrame) -> Phase:
    a, b, c, alpha, beta, gamma = spec["lattice"]
    symmetry = SymmetrySpec.from_point_group(
        str(spec["point_group"]), reference_frame=crystal_frame
    )
    lattice = Lattice(
        float(a),
        float(b),
        float(c),
        float(alpha),
        float(beta),
        float(gamma),
        crystal_frame=crystal_frame,
    )
    return Phase(
        name=str(spec["phase_name"]),
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal_frame,
    )


def _download_exercises_zip(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(EXERCISES_ZIP_URL, destination)
    return destination


def _extract_exercises(zip_path: Path, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as archive:
        archive.extractall(destination)
    return destination


def main() -> int:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    zip_path = _download_exercises_zip(DATASET_ROOT / "Exercises.zip")
    extract_root = _extract_exercises(zip_path, DATASET_ROOT / "Exercises")

    summary: dict[str, object] = {
        "source_url": EXERCISES_ZIP_URL,
        "datasets": {},
        "note": (
            "These inspection outputs use open LaboTex exercise datasets downloaded on demand. "
            "They are suitable for software-path validation and measured-vs-reconstructed "
            "comparison, but they are not bundled in-repo because redistribution terms are not "
            "stated explicitly on the source site."
        ),
    }

    for dataset_id, spec in DATASETS.items():
        phase = _phase_from_spec(spec, crystal_frame=crystal)
        data_path = extract_root / str(spec["filename"])
        measured = load_labotex_pole_figures(
            data_path,
            phase=phase,
            specimen_frame=specimen,
            intensity_normalization="max",
        )
        phi1_step_deg, big_phi_step_deg, phi2_step_deg = spec["dictionary_steps"]
        dictionary = OrientationSet.from_bunge_grid(
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
            phi1_step_deg=phi1_step_deg,
            big_phi_step_deg=big_phi_step_deg,
            phi2_step_deg=phi2_step_deg,
        )
        report = invert_labotex_pole_figures(
            [data_path],
            phase=phase,
            specimen_frame=specimen,
            orientation_dictionary=dictionary,
            kernel=KernelSpec(
                name="von_mises_fisher", halfwidth_deg=float(spec["kernel_halfwidth_deg"])
            ),
            intensity_normalization="max",
            include_symmetry_family=False,
            max_iterations=800,
            tolerance=1e-8,
        )
        predicted_figures: list[PoleFigure] = []
        difference_figures: list[PoleFigure] = []
        for measured_pf, predicted_vector in zip(
            measured,
            np.split(report.predicted_intensities, len(measured)),
            strict=True,
        ):
            predicted_vector = np.ascontiguousarray(predicted_vector, dtype=np.float64)
            predicted_scale = max(float(np.max(predicted_vector)), 1e-12)
            normalized_predicted = predicted_vector / predicted_scale
            predicted = PoleFigure(
                pole=measured_pf.pole,
                sample_directions=measured_pf.sample_directions,
                intensities=np.ascontiguousarray(normalized_predicted, dtype=np.float64),
                specimen_frame=measured_pf.specimen_frame,
                antipodal=measured_pf.antipodal,
                sample_symmetry=measured_pf.sample_symmetry,
                provenance=measured_pf.provenance,
            )
            difference = PoleFigure(
                pole=measured_pf.pole,
                sample_directions=measured_pf.sample_directions,
                intensities=np.abs(measured_pf.intensities - predicted.intensities),
                specimen_frame=measured_pf.specimen_frame,
                antipodal=measured_pf.antipodal,
                sample_symmetry=measured_pf.sample_symmetry,
                provenance=measured_pf.provenance,
            )
            predicted_figures.append(predicted)
            difference_figures.append(difference)

        dataset_plot_root = PLOTS_ROOT / dataset_id
        dataset_plot_root.mkdir(parents=True, exist_ok=True)
        for index, (measured_pf, predicted_pf, difference_pf) in enumerate(
            zip(measured, predicted_figures, difference_figures, strict=True),
            start=1,
        ):
            h, k, ell = measured_pf.pole.miller.indices.tolist()
            measured_plot = plot_pole_figure(
                measured_pf,
                kind="contour",
                method="equal_area",
                bins=81,
                sigma_bins=1.25,
                levels=12,
                title=f"{dataset_id} measured ({h}{k}{ell})",
            )
            measured_plot.savefig(
                dataset_plot_root / f"pf_{index:02d}_{h}{k}{ell}_measured.png",
                dpi=220,
                bbox_inches="tight",
                facecolor="white",
            )
            plt.close(measured_plot)

            predicted_plot = plot_pole_figure(
                predicted_pf,
                kind="contour",
                method="equal_area",
                bins=81,
                sigma_bins=1.25,
                levels=12,
                title=f"{dataset_id} reconstructed ({h}{k}{ell})",
            )
            predicted_plot.savefig(
                dataset_plot_root / f"pf_{index:02d}_{h}{k}{ell}_reconstructed.png",
                dpi=220,
                bbox_inches="tight",
                facecolor="white",
            )
            plt.close(predicted_plot)

            difference_plot = plot_pole_figure(
                difference_pf,
                kind="contour",
                method="equal_area",
                bins=81,
                sigma_bins=1.25,
                levels=12,
                title=f"{dataset_id} abs difference ({h}{k}{ell})",
            )
            difference_plot.savefig(
                dataset_plot_root / f"pf_{index:02d}_{h}{k}{ell}_difference.png",
                dpi=220,
                bbox_inches="tight",
                facecolor="white",
            )
            plt.close(difference_plot)

        odf_plot = plot_odf(
            report.odf,
            kind="sections",
            section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0),
            section_phi1_steps=91,
            section_big_phi_steps=61,
            levels=10,
            title=f"{dataset_id} ODF sections",
        )
        odf_plot.savefig(
            dataset_plot_root / "odf_sections.png",
            dpi=220,
            bbox_inches="tight",
            facecolor="white",
        )
        plt.close(odf_plot)

        summary["datasets"][dataset_id] = {
            "source_file": str(data_path).replace("\\", "/"),
            "pole_count": len(measured),
            "poles": [pf.pole.miller.indices.tolist() for pf in measured],
            "dictionary_size": report.dictionary_size,
            "observation_count": report.observation_count,
            "dictionary_steps": [phi1_step_deg, big_phi_step_deg, phi2_step_deg],
            "kernel_halfwidth_deg": float(spec["kernel_halfwidth_deg"]),
            "iterations": report.iterations,
            "converged": report.converged,
            "relative_residual_norm": report.relative_residual_norm,
            "mean_absolute_error": report.mean_absolute_error,
            "max_absolute_error": report.max_absolute_error,
            "plot_comparison_normalization": "per_pole_figure_max",
        }

    (PLOTS_ROOT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
