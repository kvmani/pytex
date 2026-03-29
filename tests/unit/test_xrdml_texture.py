from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pytex import (
    FrameDomain,
    Handedness,
    KernelSpec,
    Lattice,
    ODF,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    invert_xrdml_pole_figures,
    load_xrdml_pole_figure,
    read_xrdml_pole_figure,
    spherical_angles_to_directions,
)
from pytex.core.lattice import CrystalPlane, MillerIndex

REPO_ROOT = Path(__file__).resolve().parents[2]


def make_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame(
        "specimen",
        FrameDomain.SPECIMEN,
        ("x", "y", "z"),
        Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def write_xrdml_fixture(
    path: Path,
    *,
    phi_deg: np.ndarray,
    psi_deg: np.ndarray,
    intensity_grid: np.ndarray,
    wavelength_angstrom: float = 1.5406,
    two_theta_deg: float = 45.0,
    omega_deg: float = 22.5,
) -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/1.3" status="Completed">',
        "  <sample>",
        "    <name>synthetic_pf</name>",
        "  </sample>",
        '  <xrdMeasurement measurementType="Area measurement" measurementStepAxis="Psi" sampleMode="Reflection">',
        "    <usedWavelength>",
        f"      <kAlpha1 unit=\"Angstrom\">{wavelength_angstrom:.6f}</kAlpha1>",
        "    </usedWavelength>",
    ]
    for phi_row, psi_row, intensity_row in zip(phi_deg, psi_deg, intensity_grid, strict=True):
        lines.extend(
            [
                '    <scan mode="Continuous" scanAxis="Phi" status="Completed">',
                "      <dataPoints>",
                "        <commonCountingTime>1.0</commonCountingTime>",
                "        <intensities>"
                + " ".join(f"{value:.8f}" for value in intensity_row)
                + "</intensities>",
                '        <positions axis="2Theta">',
                f"          <commonPosition>{two_theta_deg:.4f}</commonPosition>",
                "        </positions>",
                '        <positions axis="Omega">',
                f"          <commonPosition>{omega_deg:.4f}</commonPosition>",
                "        </positions>",
                '        <positions axis="Phi">',
                f"          <startPosition>{float(phi_row[0]):.4f}</startPosition>",
                f"          <endPosition>{float(phi_row[-1]):.4f}</endPosition>",
                "        </positions>",
                '        <positions axis="Psi">',
                f"          <commonPosition>{float(psi_row[0]):.4f}</commonPosition>",
                "        </positions>",
                "      </dataPoints>",
                "    </scan>",
            ]
        )
    lines.extend(["  </xrdMeasurement>", "</xrdMeasurements>"])
    path.write_text("\n".join(lines), encoding="utf-8")


def test_read_xrdml_pole_figure_parses_grid_and_metadata(tmp_path: Path) -> None:
    phi = np.tile(np.array([0.0, 90.0, 180.0, 270.0]), (3, 1))
    psi = np.repeat(np.array([[0.0], [25.0], [55.0]]), 4, axis=1)
    intensities = np.array(
        [
            [1.0, 0.8, 0.7, 0.8],
            [0.9, 1.1, 0.95, 1.05],
            [0.4, 0.6, 0.5, 0.55],
        ],
        dtype=np.float64,
    )
    path = tmp_path / "synthetic_pf.xrdml"
    write_xrdml_fixture(path, phi_deg=phi, psi_deg=psi, intensity_grid=intensities)
    measurement = read_xrdml_pole_figure(path)
    assert measurement.shape == (3, 4)
    assert measurement.sample_mode == "Reflection"
    assert measurement.measurement_axis == "Psi"
    assert measurement.scan_axis == "Phi"
    assert measurement.wavelength_angstrom is not None
    assert_allclose(measurement.sample_directions[0], [0.0, 0.0, 1.0], atol=1e-8)
    assert_allclose(measurement.flattened_intensities[:4], intensities[0], atol=1e-8)


def test_load_xrdml_pole_figure_builds_pytex_pole_figure(tmp_path: Path) -> None:
    _, specimen, phase = make_context()
    phi = np.tile(np.array([0.0, 120.0, 240.0]), (2, 1))
    psi = np.repeat(np.array([[0.0], [35.0]]), 3, axis=1)
    intensities = np.ones((2, 3), dtype=np.float64)
    path = tmp_path / "single_pf.xrdml"
    write_xrdml_fixture(path, phi_deg=phi, psi_deg=psi, intensity_grid=intensities)
    pole = CrystalPlane(MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = load_xrdml_pole_figure(path, pole=pole, specimen_frame=specimen, antipodal=False)
    assert pole_figure.sample_directions.shape == (6, 3)
    assert pole_figure.pole == pole
    assert pole_figure.specimen_frame == specimen


def test_invert_xrdml_pole_figures_recovers_synthetic_dictionary_weights(tmp_path: Path) -> None:
    crystal, specimen, phase = make_context()
    dictionary = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(35.0, 25.0, 10.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
        ]
    )
    kernel = KernelSpec(name="von_mises_fisher", halfwidth_deg=9.0)
    true_odf = ODF.from_orientations(dictionary, weights=[3.0, 1.0], kernel=kernel)
    phi_values = np.array([0.0, 60.0, 120.0, 180.0, 240.0, 300.0], dtype=np.float64)
    psi_values = np.array([0.0, 25.0, 45.0, 65.0], dtype=np.float64)
    phi_grid = np.tile(phi_values, (psi_values.size, 1))
    psi_grid = np.repeat(psi_values[:, None], phi_values.size, axis=1)
    directions = spherical_angles_to_directions(psi_grid, phi_grid).reshape(-1, 3)
    poles = (
        CrystalPlane(MillerIndex([1, 0, 0], phase=phase), phase=phase),
        CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
    )
    paths = []
    for index, pole in enumerate(poles):
        intensities = true_odf.evaluate_pole_density(
            pole,
            directions,
            include_symmetry_family=False,
        ).reshape(psi_grid.shape)
        path = tmp_path / f"synthetic_{index}.xrdml"
        write_xrdml_fixture(path, phi_deg=phi_grid, psi_deg=psi_grid, intensity_grid=intensities)
        paths.append(path)
    report = invert_xrdml_pole_figures(
        paths,
        poles=poles,
        specimen_frame=specimen,
        orientation_dictionary=dictionary,
        kernel=kernel,
        regularization=1e-8,
        include_symmetry_family=False,
        max_iterations=1000,
        tolerance=1e-10,
    )
    assert_allclose(report.odf.normalized_weights, true_odf.normalized_weights, atol=5e-2)


def test_real_xrdml_fixture_from_xrayutilities_is_readable() -> None:
    fixture_path = REPO_ROOT / "fixtures" / "xrdml" / "polefig_Ge113_xrayutilities.xrdml.bz2"
    measurement = read_xrdml_pole_figure(fixture_path)
    assert measurement.shape == (91, 1199)
    assert measurement.wavelength_angstrom is not None
    assert measurement.sample_mode == "Reflection"
