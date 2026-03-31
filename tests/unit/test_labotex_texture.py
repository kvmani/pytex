from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pytex import (
    ODF,
    FrameDomain,
    Handedness,
    KernelSpec,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
    invert_labotex_pole_figures,
    load_labotex_pole_figures,
    read_labotex_pole_figures,
)
from pytex.adapters.labotex import LaboTexPoleFigureDescriptor
from pytex.core.lattice import CrystalPlane, MillerIndex


def make_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("cubic_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def write_labotex_fixture(
    path: Path,
    *,
    title: str,
    descriptors: list[LaboTexPoleFigureDescriptor],
    intensity_grids: list[np.ndarray],
    suffix: str,
) -> None:
    lines = [
        title,
        "Synthetic validation case",
        "Structure Code      a         b         c       ALPHA     BETA    GAMMA",
        "            7         1         1         1        90        90        90",
        f" {len(descriptors)}  number of Pole figures. Information depth- "
        ",Mode: continous scan *SP",
        "    2TH   A-first   A-last  A-step   B-first   B-last   B-step   "
        "scan   H  K  L  Background/Figure",
    ]
    for descriptor in descriptors:
        lines.append(
            f"{descriptor.two_theta_deg:7.3f}    "
            f"{descriptor.alpha_start_deg:4.1f}      "
            f"{descriptor.alpha_end_deg:4.1f}     "
            f"{descriptor.alpha_step_deg:4.1f}      "
            f"{descriptor.beta_start_deg:4.1f}     "
            f"{descriptor.beta_end_deg:4.1f}     "
            f"{descriptor.beta_step_deg:4.1f}       "
            f"{descriptor.scan_index:d}      "
            f"{descriptor.hkl[0]} {descriptor.hkl[1]} {descriptor.hkl[2]}             "
            f"{descriptor.background_flag}   {descriptor.figure_flag}"
        )
    for grid in intensity_grids:
        flattened = np.asarray(grid, dtype=np.float64).reshape(-1)
        for start in range(0, flattened.size, 8):
            lines.append(" ".join(f"{value:10.4f}" for value in flattened[start : start + 8]))
    destination = path.with_suffix(suffix)
    destination.write_text("\n".join(lines), encoding="utf-8")


def test_read_labotex_pole_figures_parses_ppf_layout(tmp_path: Path) -> None:
    descriptor_111 = LaboTexPoleFigureDescriptor(
        55.0, 0.0, 10.0, 5.0, 0.0, 15.0, 5.0, 0, (1, 1, 1), 1, 1
    )
    descriptor_200 = LaboTexPoleFigureDescriptor(
        65.0, 0.0, 10.0, 5.0, 0.0, 15.0, 5.0, 0, (2, 0, 0), 1, 1
    )
    grids = [
        np.arange(np.prod(descriptor_111.shape), dtype=np.float64).reshape(descriptor_111.shape),
        (10.0 + np.arange(np.prod(descriptor_200.shape), dtype=np.float64)).reshape(
            descriptor_200.shape
        ),
    ]
    write_labotex_fixture(
        tmp_path / "cu_case",
        title="Cu pole figures",
        descriptors=[descriptor_111, descriptor_200],
        intensity_grids=grids,
        suffix=".PPF",
    )
    measurement = read_labotex_pole_figures(tmp_path / "cu_case.PPF")
    assert measurement.format_kind == "PPF"
    assert len(measurement.descriptors) == 2
    assert measurement.descriptors[0].hkl == (1, 1, 1)
    assert measurement.intensity_grids[1].shape == descriptor_200.shape
    assert measurement.metadata["pole_figure_count"] == "2"


def test_load_labotex_pole_figures_builds_pytex_pole_figures(tmp_path: Path) -> None:
    _, specimen, phase = make_context()
    descriptor = LaboTexPoleFigureDescriptor(
        55.0, 0.0, 10.0, 5.0, 0.0, 15.0, 5.0, 0, (1, 1, 1), 1, 1
    )
    grid = np.arange(np.prod(descriptor.shape), dtype=np.float64).reshape(descriptor.shape) + 1.0
    write_labotex_fixture(
        tmp_path / "al_case",
        title="Al pole figures",
        descriptors=[descriptor],
        intensity_grids=[grid],
        suffix=".epf",
    )
    pole_figures = load_labotex_pole_figures(
        tmp_path / "al_case.epf",
        phase=phase,
        specimen_frame=specimen,
        intensity_normalization="max",
    )
    assert len(pole_figures) == 1
    assert pole_figures[0].pole.miller.indices.tolist() == [1, 1, 1]
    assert_allclose(np.max(pole_figures[0].intensities), 1.0, atol=1e-12)


def test_invert_labotex_pole_figures_recovers_synthetic_dictionary_weights(tmp_path: Path) -> None:
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
    descriptors = [
        LaboTexPoleFigureDescriptor(55.0, 0.0, 20.0, 10.0, 0.0, 30.0, 10.0, 0, (1, 0, 0), 1, 1),
        LaboTexPoleFigureDescriptor(65.0, 0.0, 20.0, 10.0, 0.0, 30.0, 10.0, 0, (1, 1, 1), 1, 1),
    ]
    planes = [
        CrystalPlane(MillerIndex(descriptor.hkl, phase=phase), phase=phase)
        for descriptor in descriptors
    ]
    grids = []
    for descriptor, plane in zip(descriptors, planes, strict=True):
        densities = true_odf.evaluate_pole_density(
            plane,
            descriptor.sample_directions,
            include_symmetry_family=False,
        )
        grids.append(densities.reshape(descriptor.shape))
    write_labotex_fixture(
        tmp_path / "fe_case",
        title="Fe pole figures",
        descriptors=descriptors,
        intensity_grids=grids,
        suffix=".epf",
    )
    report = invert_labotex_pole_figures(
        [tmp_path / "fe_case.epf"],
        phase=phase,
        specimen_frame=specimen,
        orientation_dictionary=dictionary,
        kernel=kernel,
        regularization=1e-8,
        include_symmetry_family=False,
        intensity_normalization="none",
        max_iterations=1000,
        tolerance=1e-10,
    )
    assert_allclose(report.odf.normalized_weights, true_odf.normalized_weights, atol=5e-2)
    assert report.relative_residual_norm <= 6e-2
