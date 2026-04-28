from __future__ import annotations

import numpy as np
import pytest

from pytex import (
    AtomicSite,
    CrystalMap,
    CrystalPlane,
    EBSDTextureWorkflow,
    EBSDTextureWorkflowResult,
    EulerConventionTransform,
    FrameDomain,
    Handedness,
    IPFSectorBoundary,
    Lattice,
    MillerIndex,
    ODFReconstructionConfig,
    Orientation,
    OrientationQualityWeights,
    OrientationRelationship,
    OrientationSet,
    ParentReconstructionConfig,
    ParentReconstructionReport,
    Phase,
    PhaseTransformationRecord,
    PoleFigure,
    PoleFigureCorrectionSpec,
    PoleFigureResidualReport,
    RadiationSpec,
    ReferenceFrame,
    ReflectionCondition,
    Rotation,
    ScatteringFactorTable,
    StructureFactor,
    SymmetrySpec,
    UnitCell,
    VariantSelectionReport,
    from_json_contract,
    lorentz_polarization_factor,
    reconstruct_parent_orientation,
    to_json_contract,
)
from pytex.diffraction import DiffractionIntensityModel


def _phase(name: str = "ni-fcc") -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    lattice = Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    unit_cell = UnitCell(
        lattice=lattice,
        sites=(
            AtomicSite("Ni1", "Ni", np.array([0.0, 0.0, 0.0])),
            AtomicSite("Ni2", "Ni", np.array([0.0, 0.5, 0.5])),
            AtomicSite("Ni3", "Ni", np.array([0.5, 0.0, 0.5])),
            AtomicSite("Ni4", "Ni", np.array([0.5, 0.5, 0.0])),
        ),
    )
    phase = Phase(
        name=name,
        lattice=lattice,
        symmetry=symmetry,
        crystal_frame=crystal,
        unit_cell=unit_cell,
        space_group_symbol="Fm-3m",
    )
    return crystal, specimen, phase


def _orientations(count: int = 3) -> OrientationSet:
    crystal, specimen, phase = _phase()
    angles = np.array(
        [
            [0.0, 0.0, 0.0],
            [20.0, 30.0, 10.0],
            [70.0, 15.0, 25.0],
        ],
        dtype=np.float64,
    )[:count]
    return OrientationSet.from_euler_angles(
        angles,
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )


def test_texture_reconstruction_v2_correction_residual_and_contracts() -> None:
    orientations = _orientations()
    phase = orientations.phase
    assert phase is not None
    pole = CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(orientations, pole, weights=[1.0, 2.0, 3.0])
    correction = PoleFigureCorrectionSpec(scale=2.0, background=0.25)
    corrected = correction.apply(pole_figure)
    assert np.all(corrected.intensities >= 0.0)

    config = ODFReconstructionConfig(algorithm="discrete", correction=correction)
    report = config.reconstruct([pole_figure], orientation_dictionary=orientations)
    assert report.observation_count == len(corrected.intensities)
    residual = PoleFigureResidualReport.from_odf(report.odf, corrected)
    assert residual.observation_count == len(corrected.intensities)

    for obj in (corrected, report.odf, report):
        payload = to_json_contract(obj)
        assert to_json_contract(from_json_contract(payload)) == payload


def test_orientation_geometry_and_euler_convention_surface() -> None:
    crystal, _, phase = _phase()
    sector = IPFSectorBoundary.from_symmetry(phase.symmetry)
    assert sector.contains([0.0, 0.0, 1.0])
    assert sector.boundary_equations

    bunge_to_kocks = EulerConventionTransform(source="bunge", target="kocks")
    assert np.allclose(
        bunge_to_kocks.apply([301.0, 36.7, 26.63]),
        np.array([211.0, 36.7, 63.37]),
    )
    bunge_to_roe = EulerConventionTransform(source="bunge", target="roe")
    assert np.allclose(
        bunge_to_roe.apply([301.0, 36.7, 26.63]),
        np.array([211.0, 36.7, 116.63]),
    )
    assert sector.symmetry.reference_frame == crystal


def test_structure_to_diffraction_physics_layer() -> None:
    _, _, phase = _phase()
    condition = ReflectionCondition.from_phase(phase)
    assert not condition.is_allowed([1, 0, 0])
    assert condition.is_allowed([1, 1, 1])

    table = ScatteringFactorTable(model="atomic_number")
    structure_factor = StructureFactor.from_phase_hkl(phase, [1, 1, 1], scattering_table=table)
    assert structure_factor.amplitude > 0.0

    model = DiffractionIntensityModel(scattering_table=table, reflection_condition=condition)
    intensity = model.intensity(
        phase,
        [1, 1, 1],
        two_theta_rad=np.deg2rad(44.0),
        multiplicity=8,
        radiation=RadiationSpec.cu_ka(),
    )
    assert intensity > 0.0
    assert model.intensity(phase, [1, 0, 0], two_theta_rad=np.deg2rad(20.0)) == 0.0


def test_diffraction_physics_rejects_invalid_inputs() -> None:
    _, _, phase = _phase()
    table = ScatteringFactorTable(model="unit")
    model = DiffractionIntensityModel(scattering_table=table)

    with pytest.raises(ValueError, match="g_magnitude"):
        table.scattering_factor("Ni", -1.0)
    with pytest.raises(ValueError, match="two_theta_rad"):
        lorentz_polarization_factor(0.0)
    with pytest.raises(ValueError, match="two_theta_rad"):
        model.intensity(phase, [1, 1, 1], two_theta_rad=np.pi)
    with pytest.raises(ValueError, match="multiplicity"):
        model.intensity(phase, [1, 1, 1], two_theta_rad=np.deg2rad(44.0), multiplicity=0)
    with pytest.raises(ValueError, match="finite real and imaginary"):
        StructureFactor(
            miller_indices=[1, 1, 1],
            value=complex(np.nan, 0.0),
            amplitude=1.0,
            phase_rad=0.0,
            phase=phase,
        )


def test_workflow_grade_ebsd_to_texture_pipeline_and_contracts() -> None:
    orientations = _orientations()
    coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    crystal_map = CrystalMap(
        coordinates=coordinates,
        orientations=orientations,
        map_frame=orientations.specimen_frame,
    )
    phase = orientations.phase
    assert phase is not None
    pole = CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase)
    workflow = EBSDTextureWorkflow(
        poles=(pole,),
        weights=OrientationQualityWeights([0.2, 0.3, 0.5]),
    )
    result = workflow.run(crystal_map)
    assert result.orientation_count == 3
    assert np.isclose(np.sum(result.weights), 1.0)
    assert len(result.texture_report.pole_figures) == 1

    for obj in (crystal_map, result.texture_report):
        payload = to_json_contract(obj)
        assert to_json_contract(from_json_contract(payload)) == payload


def test_ebsd_texture_workflow_rejects_invalid_weight_and_threshold_states() -> None:
    orientations = _orientations()
    coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    crystal_map = CrystalMap(
        coordinates=coordinates,
        orientations=orientations,
        map_frame=orientations.specimen_frame,
    )
    phase = orientations.phase
    assert phase is not None
    pole = CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase)
    report = crystal_map.texture_report(poles=(pole,))

    with pytest.raises(ValueError, match="segmentation_threshold_deg"):
        EBSDTextureWorkflow(segmentation_threshold_deg=0.0)
    with pytest.raises(ValueError, match="sample_directions"):
        EBSDTextureWorkflow(sample_directions=())
    with pytest.raises(ValueError, match="No valid positive orientation weights"):
        OrientationQualityWeights([1.0, 1.0, 1.0], valid_mask=[False, False, False]).for_count(3)
    with pytest.raises(ValueError, match="finite and non-negative"):
        EBSDTextureWorkflowResult(
            crystal_map=crystal_map,
            texture_report=report,
            odf=report.odf,
            weights=np.array([1.0, np.nan, 0.0]),
        )


def test_stable_parent_reconstruction_track() -> None:
    crystal, specimen, parent_phase = _phase("parent")
    child_phase = Phase(
        name="child",
        lattice=parent_phase.lattice,
        symmetry=parent_phase.symmetry,
        crystal_frame=crystal,
    )
    relationship = OrientationRelationship(
        name="identity_demo",
        parent_phase=parent_phase,
        child_phase=child_phase,
        parent_to_child_rotation=Rotation.identity(),
    )
    parent = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=parent_phase.symmetry,
        phase=parent_phase,
    )
    child = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    candidate_parents = OrientationSet.from_orientations(
        [
            parent,
            Orientation(
                rotation=Rotation.from_bunge_euler(45.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=parent_phase.symmetry,
                phase=parent_phase,
            ),
        ]
    )
    record = PhaseTransformationRecord(
        name="demo",
        orientation_relationship=relationship,
        parent_orientation=parent,
        child_orientations=OrientationSet.from_orientations([child]),
    )
    report = reconstruct_parent_orientation(
        record,
        candidate_parents,
        config=ParentReconstructionConfig(ambiguity_tolerance_deg=0.1),
    )
    assert report.best_index == 0
    assert report.best_score_deg <= 1e-8


def test_parent_reconstruction_reports_reject_invalid_states() -> None:
    crystal, specimen, parent_phase = _phase("parent")
    child_phase = Phase(
        name="child",
        lattice=parent_phase.lattice,
        symmetry=parent_phase.symmetry,
        crystal_frame=crystal,
    )
    parent = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=parent_phase.symmetry,
        phase=parent_phase,
    )
    child = Orientation(
        rotation=Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=child_phase.symmetry,
        phase=child_phase,
    )
    candidate_parents = OrientationSet.from_orientations([parent])
    record = PhaseTransformationRecord(
        name="demo",
        orientation_relationship=OrientationRelationship(
            name="identity_demo",
            parent_phase=parent_phase,
            child_phase=child_phase,
            parent_to_child_rotation=Rotation.identity(),
        ),
        parent_orientation=parent,
        child_orientations=OrientationSet.from_orientations([child]),
    )

    with pytest.raises(ValueError, match="ambiguity_tolerance_deg"):
        ParentReconstructionConfig(ambiguity_tolerance_deg=np.inf)
    with pytest.raises(ValueError, match="variant_indices"):
        VariantSelectionReport(variant_indices=[0], scores_deg=[0.0])
    with pytest.raises(ValueError, match="scores_deg"):
        ParentReconstructionReport(
            record=record,
            candidate_parents=candidate_parents,
            scores_deg=np.array([np.nan]),
            best_index=0,
            best_score_deg=0.0,
            ambiguous_indices=(0,),
            reduction="mean",
            symmetry_aware=True,
        )
