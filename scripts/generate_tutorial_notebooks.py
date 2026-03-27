from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent

REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_ROOT = REPO_ROOT / "docs" / "site" / "tutorials" / "notebooks"


def markdown_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": dedent(source).strip() + "\n",
    }


def code_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": dedent(source).strip() + "\n",
    }


def notebook(title: str, cells: list[dict[str, object]]) -> dict[str, object]:
    normalized_cells: list[dict[str, object]] = []
    for index, cell in enumerate(cells, start=1):
        normalized_cell = dict(cell)
        normalized_cell.setdefault("id", f"cell-{index:03d}")
        normalized_cells.append(normalized_cell)
    return {
        "cells": normalized_cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "version": "3.11",
            },
            "title": title,
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def write_notebook(path: Path, title: str, cells: list[dict[str, object]]) -> None:
    payload = notebook(title, cells)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def build_notebooks() -> dict[str, dict[str, object]]:
    common_setup = """
    import numpy as np

    from pytex import (
        AcquisitionGeometry,
        BenchmarkManifest,
        CalibrationRecord,
        CrystalMap,
        CrystalPlane,
        DiffractionGeometry,
        EulerSet,
        ExperimentManifest,
        FrameDomain,
        FrameTransform,
        Handedness,
        InversePoleFigure,
        KernelSpec,
        KinematicSimulation,
        Lattice,
        MeasurementQuality,
        MillerIndex,
        ODF,
        Orientation,
        OrientationRelationship,
        OrientationSet,
        Phase,
        PhaseTransformationRecord,
        PoleFigure,
        ReferenceFrame,
        Rotation,
        ScatteringSetup,
        SymmetrySpec,
        TransformationVariant,
        ValidationManifest,
        VectorSet,
        WorkflowResultManifest,
        ZoneAxis,
        plot_odf,
        plot_inverse_pole_figure,
        plot_orientations,
        plot_pole_figure,
        plot_symmetry_elements,
        plot_symmetry_orbit,
        plot_vector_set,
    )


    def make_context():
        crystal = ReferenceFrame(
            "crystal",
            FrameDomain.CRYSTAL,
            ("a", "b", "c"),
            Handedness.RIGHT,
        )
        specimen = ReferenceFrame(
            "specimen",
            FrameDomain.SPECIMEN,
            ("x", "y", "z"),
            Handedness.RIGHT,
        )
        map_frame = ReferenceFrame(
            "map",
            FrameDomain.MAP,
            ("i", "j", "k"),
            Handedness.RIGHT,
        )
        detector = ReferenceFrame(
            "detector",
            FrameDomain.DETECTOR,
            ("u", "v", "n"),
            Handedness.RIGHT,
        )
        lab = ReferenceFrame(
            "lab",
            FrameDomain.LABORATORY,
            ("X", "Y", "Z"),
            Handedness.RIGHT,
        )
        symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
        lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
        phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
        return crystal, specimen, map_frame, detector, lab, phase
    """

    notebooks: dict[str, dict[str, object]] = {}

    notebooks["01_reference_frames_and_transforms.ipynb"] = {
        "title": "Reference Frames And Transforms",
        "cells": [
            markdown_cell(
                """
                # Reference Frames And Transforms

                This notebook explains the first semantic layer of PyTex: named frames, explicit
                domains, and reusable rigid transforms. The goal is not only to show the API, but
                to fix the meaning of every vector before later orientation, EBSD, or diffraction
                workflows are attempted.

                ## Learning Goals

                - construct `ReferenceFrame` objects and understand frame domains
                - express rigid mappings with `FrameTransform`
                - work with `VectorSet` so frame meaning stays attached to batched data
                - connect the code to the mathematical map

                ## Mathematical Contract

                PyTex treats a frame transform as

                $$
                \\mathbf{v}^{(t)} = \\mathbf{R}_{s\\rightarrow t}\\,\\mathbf{v}^{(s)} + \\mathbf{t}_{s\\rightarrow t},
                $$

                where the source and target frames are explicit objects, not comments or naming
                conventions.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()

                specimen_to_map = FrameTransform(
                    source=specimen,
                    target=map_frame,
                    rotation_matrix=np.eye(3),
                    translation_vector=np.array([0.0, 0.0, 0.0]),
                )

                vectors = VectorSet(
                    values=np.array(
                        [
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                        ]
                    ),
                    reference_frame=specimen,
                )

                mapped = specimen_to_map.apply_to_vectors(vectors)

                print(vectors.reference_frame.name, "->", mapped.reference_frame.name)
                print(mapped.values)
                """
            ),
            markdown_cell(
                """
                ## Why The `VectorSet` Matters

                A raw `(n, 3)` array is not enough when the frame matters. PyTex therefore treats
                batched vectors as a first-class semantic object:

                - the storage remains NumPy-backed and vectorized
                - the shared frame is explicit
                - later transforms can reject mismatches early

                This is the pattern that the rest of the library follows: preserve vectorization,
                but never let the science collapse back into anonymous arrays.
                """
            ),
            code_cell(
                """
                try:
                    wrong_vectors = VectorSet(values=np.array([[1.0, 0.0, 0.0]]), reference_frame=crystal)
                    specimen_to_map.apply_to_vectors(wrong_vectors)
                except ValueError as exc:
                    print(type(exc).__name__)
                    print(exc)
                """
            ),
        ],
    }

    notebooks["02_rotations_orientations_and_batch_primitives.ipynb"] = {
        "title": "Rotations, Orientations, And Batch Primitives",
        "cells": [
            markdown_cell(
                """
                # Rotations, Orientations, And Batch Primitives

                PyTex distinguishes between a mathematical rotation and an orientation that lives
                between a crystal frame and a specimen frame. This notebook also shows how the batch
                primitives let large vectorized workflows stay semantically typed.

                ## Key Mathematical Ideas

                A unit quaternion is stored in `w, x, y, z` order and defines the active rotation

                $$
                \\mathbf{v}' = \\mathbf{R}(q)\\,\\mathbf{v}.
                $$

                An `Orientation` then wraps that rotation together with the frame pair
                `(crystal_frame, specimen_frame)`.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, *_ = make_context()

                eulers = np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [45.0, 35.0, 15.0],
                        [90.0, 0.0, 0.0],
                    ]
                )

                euler_set = EulerSet(eulers, convention="bunge", degrees=True)
                rotation_set = euler_set.to_rotation_set()
                quaternion_set = rotation_set.as_quaternion_set()
                orientation_set = OrientationSet.from_euler_angles(
                    euler_set,
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
                )

                print("Euler angles")
                print(euler_set.as_array())
                print("Quaternions")
                print(quaternion_set.as_array())
                """
            ),
            code_cell(
                """
                crystal_vectors = VectorSet(
                    values=np.array(
                        [
                            [1.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                        ]
                    ),
                    reference_frame=crystal,
                )

                specimen_vectors = orientation_set.map_crystal_directions(crystal_vectors)
                print(specimen_vectors.reference_frame.name)
                print(specimen_vectors.values)
                """
            ),
            markdown_cell(
                """
                ## Batch Semantics

                The important point is not just that these operations are vectorized. It is that the
                shared convention or frame survives the operation:

                - `EulerSet` retains the Euler convention
                - `QuaternionSet` retains canonical quaternion normalization
                - `RotationSet` rotates many vectors at once
                - `OrientationSet` carries crystal/specimen meaning while doing batched direction maps
                """
            ),
        ],
    }

    notebooks["03_symmetry_and_fundamental_regions.ipynb"] = {
        "title": "Symmetry And Fundamental Regions",
        "cells": [
            markdown_cell(
                """
                # Symmetry And Fundamental Regions

                Symmetry is where many texture workflows become ambiguous in other tools. PyTex
                treats `SymmetrySpec` as a reusable scientific object and exposes explicit
                vector- and orientation-reduction paths.

                ## Mathematical Idea

                A symmetry orbit acts on a direction or orientation by group multiplication. A
                reduction then selects one representative from that orbit according to an explicit
                rule rather than a hidden plotting convention.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, *_ = make_context()
                symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)

                directions = VectorSet(
                    values=np.array(
                        [
                            [1.0, 1.0, 1.0],
                            [-1.0, 1.0, 1.0],
                            [1.0, -1.0, 1.0],
                            [1.0, 1.0, -1.0],
                        ]
                    ),
                    reference_frame=crystal,
                ).normalized()

                reduced = symmetry.reduce_vectors_to_fundamental_sector(directions, antipodal=True)
                print(reduced.values)
                """
            ),
            code_cell(
                """
                orientations = OrientationSet.from_euler_angles(
                    np.array(
                        [
                            [0.0, 0.0, 0.0],
                            [90.0, 0.0, 0.0],
                            [35.0, 20.0, 10.0],
                        ]
                    ),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=symmetry,
                )

                projected = orientations.projected_to_fundamental_region()
                keys = projected.exact_fundamental_region_keys()
                print(keys)
                """
            ),
            markdown_cell(
                """
                The combination of `SymmetrySpec` plus explicit batch reduction is what lets PyTex
                explain exactly how an inverse pole figure sector or orientation fundamental region
                was obtained.
                """
            ),
        ],
    }

    notebooks["04_phases_lattices_space_groups_and_cif.ipynb"] = {
        "title": "Phases, Lattices, Space Groups, And CIF Import",
        "cells": [
            markdown_cell(
                """
                # Phases, Lattices, Space Groups, And CIF Import

                PyTex keeps structure semantics explicit by separating point-group reduction
                semantics from structure-facing space-group identity.

                ## Key Rule

                - `SymmetrySpec` is the orientation-reduction surface
                - `SpaceGroupSpec` is the structure-definition surface
                - `Phase` owns both when the structure source provides them
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()

                print("Phase:", phase.name)
                print("Point group:", phase.symmetry.point_group)
                print("Direct basis:")
                print(phase.lattice.direct_basis().matrix)
                print("Reciprocal basis:")
                print(phase.lattice.reciprocal_basis().matrix)
                """
            ),
            code_cell(
                '''
                nacl_cif = """
                data_NaCl
                _symmetry_space_group_name_H-M 'F m -3 m'
                _symmetry_Int_Tables_number 225
                _cell_length_a 5.6402
                _cell_length_b 5.6402
                _cell_length_c 5.6402
                _cell_angle_alpha 90
                _cell_angle_beta 90
                _cell_angle_gamma 90
                loop_
                  _atom_site_label
                  _atom_site_type_symbol
                  _atom_site_fract_x
                  _atom_site_fract_y
                  _atom_site_fract_z
                  Na1 Na 0.0 0.0 0.0
                  Cl1 Cl 0.5 0.5 0.5
                """

                try:
                    imported_phase = Phase.from_cif_string(nacl_cif, crystal_frame=make_context()[0])
                    print(imported_phase.name)
                    print(imported_phase.space_group.symbol, imported_phase.space_group.number)
                    print(imported_phase.symmetry.point_group)
                except ImportError:
                    print("Install the adapters extra to execute CIF-backed examples.")
                '''
            ),
        ],
    }

    notebooks["05_multimodal_acquisition_and_manifests.ipynb"] = {
        "title": "Multimodal Acquisition, Calibration, And Manifests",
        "cells": [
            markdown_cell(
                """
                # Multimodal Acquisition, Calibration, And Manifests

                PyTex now has a stable shared experiment layer. This notebook shows how acquisition
                semantics, calibration state, quality metadata, and machine-readable manifests fit
                together.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()

                acquisition = AcquisitionGeometry(
                    specimen_frame=specimen,
                    modality="ebsd",
                    map_frame=map_frame,
                    specimen_to_map=FrameTransform(
                        source=specimen,
                        target=map_frame,
                        rotation_matrix=np.eye(3),
                    ),
                    calibration_record=CalibrationRecord(
                        source="stage-fit",
                        status="calibrated",
                        residual_error=0.1,
                    ),
                    measurement_quality=MeasurementQuality(
                        confidence=0.95,
                        valid_fraction=0.99,
                        uncertainty={"tilt_deg": 0.1},
                    ),
                )

                experiment = ExperimentManifest.from_acquisition_geometry(
                    acquisition,
                    source_system="pytex",
                    phase=phase,
                    referenced_files=("scan.ang",),
                )
                print(experiment.to_dict()["schema_id"])
                print(experiment.to_dict()["modality"])
                """
            ),
            code_cell(
                """
                benchmark = BenchmarkManifest(
                    benchmark_id="ebsd_regular_grid_demo",
                    subsystem="ebsd",
                    baseline_kind="internal_plus_mtex",
                    workflows=("kam", "segmentation", "grod"),
                    tolerances={"misorientation_atol_deg": 1e-6},
                )

                validation = ValidationManifest(
                    campaign_name="texture_validation_demo",
                    subsystem="texture",
                    baseline_kind="mtex_plus_internal",
                    status="implemented",
                    reference_ids=("MTEX", "Bunge"),
                    linked_benchmark_ids=(benchmark.benchmark_id,),
                )

                workflow = WorkflowResultManifest(
                    result_id="demo_result",
                    workflow_name="ebsd_pipeline",
                    modality="ebsd",
                    produced_by="pytex",
                    input_manifest_ids=(experiment.schema_id,),
                    artifact_paths=("results/demo_figure.svg",),
                )

                print(benchmark.to_dict()["workflows"])
                print(validation.to_dict()["reference_ids"])
                print(workflow.to_dict()["artifact_paths"])
                """
            ),
        ],
    }

    notebooks["06_texture_odf_and_pole_figure_inversion.ipynb"] = {
        "title": "Texture, ODF Evaluation, And Pole-Figure Inversion",
        "cells": [
            markdown_cell(
                """
                # Texture, ODF Evaluation, And Pole-Figure Inversion

                This notebook demonstrates the current PyTex texture model:

                - pole-figure synthesis from orientation sets
                - inverse-pole-figure synthesis
                - kernel-weighted ODF evaluation
                - discrete pole-figure inversion through an explicit orientation dictionary
                - contour pole figures and Bunge-section ODF plotting through the runtime API

                ## Discrete Forward Model

                For measurement directions $\\mathbf{s}_m$ and dictionary orientations $g_j$, the
                current inversion builds a forward matrix

                $$
                A_{mj} = \\frac{1}{|\\mathcal{H}|}\\sum_{h\\in\\mathcal{H}} K\\bigl(\\angle(\\mathbf{s}_m, g_j h)\\bigr),
                $$

                then solves a regularized nonnegative least-squares problem over the dictionary
                weights.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, *_ , phase = make_context()
                symmetry = phase.symmetry

                support = OrientationSet.from_euler_angles(
                    np.array(
                        [
                            [0.0, 0.0, 0.0],
                            [90.0, 0.0, 0.0],
                            [35.0, 25.0, 10.0],
                        ]
                    ),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=symmetry,
                    phase=phase,
                )

                odf = ODF.from_orientations(
                    support,
                    weights=[4.0, 2.0, 1.0],
                    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
                )

                query = Orientation(
                    rotation=Rotation.identity(),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=symmetry,
                    phase=phase,
                )
                print("ODF(query) =", odf.evaluate(query))
                print("Volume fraction within 15 deg =", odf.volume_fraction(query, max_angle_deg=15.0))
                """
            ),
            code_cell(
                """
                poles = (
                    CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase),
                    CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase),
                )

                pole_figures = odf.reconstruct_pole_figures(
                    poles,
                    include_symmetry_family=False,
                    antipodal=False,
                )

                inversion = ODF.invert_pole_figures(
                    pole_figures,
                    orientation_dictionary=support,
                    kernel=odf.kernel,
                    regularization=1e-8,
                    include_symmetry_family=False,
                )

                print("Converged:", inversion.converged)
                print("Estimated weights:", inversion.odf.normalized_weights)
                print("Residual norm:", inversion.residual_norm)
                """
            ),
            code_cell(
                """
                ipf = InversePoleFigure.from_orientations(
                    support,
                    [0.0, 0.0, 1.0],
                    reduce_by_symmetry=True,
                    antipodal=True,
                )
                print(ipf.crystal_directions)
                """
            ),
            code_cell(
                """
                contour_pf = plot_pole_figure(
                    pole_figures[0],
                    kind="contour",
                    bins=81,
                    sigma_bins=1.5,
                    levels=12,
                    title="Pole Figure Contours",
                )
                odf_sections = plot_odf(
                    inversion.odf,
                    kind="sections",
                    section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
                    section_phi1_steps=121,
                    section_big_phi_steps=61,
                    levels=10,
                    title="Estimated ODF Bunge Sections",
                )
                contour_pf
                """
            ),
            code_cell("odf_sections"),
            markdown_cell(
                """
                The contour pole figure is built from a smoothed projected density grid over the
                reconstructed pole data. The ODF section plot is a kernel-smoothed Bunge-section
                inspection view over the discrete support rather than a harmonic ODF expansion.
                """
            ),
        ],
    }

    notebooks["07_ebsd_regular_grid_workflows.ipynb"] = {
        "title": "EBSD Regular-Grid Workflows",
        "cells": [
            markdown_cell(
                """
                # EBSD Regular-Grid Workflows

                This notebook walks through the current EBSD foundation:

                - regular-grid neighborhood construction
                - KAM
                - grain segmentation
                - GROD
                - grain-boundary and grain-graph summaries
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()
                orientations = OrientationSet.from_euler_angles(
                    np.array(
                        [
                            [0.0, 0.0, 0.0],
                            [2.0, 0.0, 0.0],
                            [25.0, 0.0, 0.0],
                            [27.0, 0.0, 0.0],
                        ]
                    ),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=phase.symmetry,
                    phase=phase,
                )

                crystal_map = CrystalMap(
                    coordinates=np.array(
                        [
                            [0.0, 0.0],
                            [1.0, 0.0],
                            [0.0, 1.0],
                            [1.0, 1.0],
                        ]
                    ),
                    orientations=orientations,
                    map_frame=map_frame,
                    grid_shape=(2, 2),
                    step_sizes=(1.0, 1.0),
                )

                print("Neighbor pairs")
                print(crystal_map.neighbor_pairs())
                print("KAM")
                print(crystal_map.kernel_average_misorientation_deg(symmetry_aware=False))
                """
            ),
            code_cell(
                """
                segmentation = crystal_map.segment_grains(
                    max_misorientation_deg=5.0,
                    symmetry_aware=False,
                )
                print("Labels")
                print(segmentation.label_grid)
                print("GROD")
                print(segmentation.grod_map_deg())

                boundary_network = segmentation.boundary_network(min_misorientation_deg=5.0)
                print("Boundary count:", boundary_network.count)
                print("Graph adjacency")
                print(boundary_network.grain_graph().adjacency_matrix)
                """
            ),
            code_cell(
                """
                experiment = crystal_map.to_experiment_manifest(source_system="pytex")
                print(experiment.to_dict()["modality"])
                print(experiment.to_dict()["metadata"])
                """
            ),
        ],
    }

    notebooks["08_diffraction_geometry_and_kinematic_spots.ipynb"] = {
        "title": "Diffraction Geometry And Kinematic Spots",
        "cells": [
            markdown_cell(
                """
                # Diffraction Geometry And Kinematic Spots

                PyTex currently treats diffraction through explicit geometry, reciprocal-space
                primitives, and kinematic spot simulation. This notebook shows the geometry layer
                and the current simulation surface.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()

                geometry = DiffractionGeometry(
                    detector_frame=detector,
                    specimen_frame=specimen,
                    laboratory_frame=lab,
                    beam_energy_kev=200.0,
                    camera_length_mm=150.0,
                    pattern_center=np.array([0.5, 0.5, 0.7]),
                    detector_pixel_size_um=(50.0, 50.0),
                    detector_shape=(512, 512),
                    scattering_setup=ScatteringSetup(laboratory_frame=lab, beam_energy_kev=200.0),
                )

                plane = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
                print("Ring radius (mm):", geometry.ring_radius_mm_for_plane(plane))
                """
            ),
            code_cell(
                """
                orientation = Orientation(
                    rotation=Rotation.identity(),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=phase.symmetry,
                    phase=phase,
                )

                simulation = KinematicSimulation.simulate_spots(
                    geometry,
                    phase,
                    np.array(
                        [
                            [1, 0, 0],
                            [1, 1, 0],
                            [1, 1, 1],
                            [2, 0, 0],
                        ]
                    ),
                    orientation=orientation,
                )

                print("Accepted spots:", len(simulation.accepted_spots()))
                print("Reflection families:", len(simulation.reflection_families))
                """
            ),
            code_cell(
                """
                experiment = geometry.to_experiment_manifest(source_system="pytex", phase=phase)
                print(experiment.to_dict()["modality"])
                print(experiment.to_dict()["metadata"]["camera_length_mm"])
                """
            ),
        ],
    }

    notebooks["09_phase_transformation_foundations.ipynb"] = {
        "title": "Phase Transformation Foundations",
        "cells": [
            markdown_cell(
                """
                # Phase Transformation Foundations

                PyTex now has the first stable primitive family for phase-transformation work:

                - `OrientationRelationship`
                - `TransformationVariant`
                - `PhaseTransformationRecord`

                This notebook explains the semantics rather than pretending the repo already has a
                full parent-reconstruction platform.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, parent_phase = make_context()
                child_crystal = ReferenceFrame(
                    "child_crystal",
                    FrameDomain.CRYSTAL,
                    ("a", "b", "c"),
                    Handedness.RIGHT,
                )
                child_symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=child_crystal)
                child_lattice = Lattice(3.2, 3.2, 3.2, 90.0, 90.0, 90.0, crystal_frame=child_crystal)
                child_phase = Phase(
                    "child_phase",
                    lattice=child_lattice,
                    symmetry=child_symmetry,
                    crystal_frame=child_crystal,
                )

                relationship = OrientationRelationship(
                    name="demo_or",
                    parent_phase=parent_phase,
                    child_phase=child_phase,
                    parent_to_child_rotation=Rotation.from_axis_angle([0.0, 0.0, 1.0], np.pi / 4.0),
                )

                variants = relationship.generate_variants()
                print("Variant count:", len(variants))
                """
            ),
            code_cell(
                """
                parent_orientation = Orientation(
                    rotation=Rotation.identity(),
                    crystal_frame=parent_phase.crystal_frame,
                    specimen_frame=specimen,
                    symmetry=parent_phase.symmetry,
                    phase=parent_phase,
                )
                child_orientations = OrientationSet.from_orientations(
                    [
                        Orientation(
                            rotation=relationship.parent_to_child_rotation,
                            crystal_frame=child_phase.crystal_frame,
                            specimen_frame=specimen,
                            symmetry=child_phase.symmetry,
                            phase=child_phase,
                        )
                    ]
                )
                record = PhaseTransformationRecord(
                    name="demo_record",
                    orientation_relationship=relationship,
                    parent_orientation=parent_orientation,
                    child_orientations=child_orientations,
                    variant_indices=np.array([1]),
                )
                print(record.variant_count)
                """
            ),
        ],
    }

    notebooks["10_plotting_semantic_primitives.ipynb"] = {
        "title": "Plotting Semantic Primitives",
        "cells": [
            markdown_cell(
                """
                # Plotting Semantic Primitives

                PyTex plotting now treats semantic objects as the plotting inputs rather than
                expecting users to manually convert everything into anonymous arrays.

                The runtime plotting surface returns normal Matplotlib figures. Canonical
                documentation assets may still be exported to SVG when they belong in the repo.
                """
            ),
            code_cell(common_setup),
            code_cell(
                """
                crystal, specimen, map_frame, detector, lab, phase = make_context()

                vectors = VectorSet(
                    values=np.array(
                        [
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                        ]
                    ),
                    reference_frame=crystal,
                )
                figure = plot_vector_set(vectors, normalize=True)
                figure
                """
            ),
            code_cell(
                """
                symmetry = phase.symmetry
                orbit_figure = plot_symmetry_orbit(symmetry, np.array([1.0, 0.0, 0.0]))
                elements_figure = plot_symmetry_elements(symmetry)
                orbit_figure
                """
            ),
            code_cell(
                """
                orientations = OrientationSet.from_euler_angles(
                    np.array(
                        [
                            [0.0, 0.0, 0.0],
                            [45.0, 35.0, 10.0],
                            [90.0, 0.0, 0.0],
                        ]
                    ),
                    crystal_frame=crystal,
                    specimen_frame=specimen,
                    symmetry=symmetry,
                    phase=phase,
                )
                orientation_figure = plot_orientations(orientations, representation="euler")
                orientation_figure
                """
            ),
            code_cell(
                """
                odf = ODF.from_orientations(orientations, weights=[4.0, 2.0, 1.0])
                pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
                pole_figure = odf.reconstruct_pole_figure(
                    pole,
                    include_symmetry_family=False,
                    antipodal=False,
                )
                ipf = InversePoleFigure.from_orientations(
                    orientations,
                    [0.0, 0.0, 1.0],
                    reduce_by_symmetry=True,
                    antipodal=True,
                )
                pf_figure = plot_pole_figure(
                    pole_figure,
                    kind="contour",
                    bins=81,
                    sigma_bins=1.5,
                    levels=12,
                )
                ipf_figure = plot_inverse_pole_figure(ipf)
                odf_figure = plot_odf(
                    odf,
                    kind="sections",
                    section_phi2_deg=(0.0, 30.0, 60.0),
                    section_phi1_steps=91,
                    section_big_phi_steps=46,
                    levels=10,
                )
                pf_figure
                """
            ),
            code_cell("ipf_figure"),
            code_cell("odf_figure"),
            markdown_cell(
                """
                The plotting notebook now shows the three main texture inspection surfaces
                implemented today:

                - contour pole figures
                - inverse pole figures
                - Bunge-section ODF plots

                All of them are produced by the same runtime plotting API that ordinary user code
                calls, rather than by notebook-specific plotting helpers.
                """
            ),
        ],
    }

    return notebooks


def main() -> int:
    notebooks = build_notebooks()
    NOTEBOOK_ROOT.mkdir(parents=True, exist_ok=True)
    for filename, payload in notebooks.items():
        write_notebook(NOTEBOOK_ROOT / filename, payload["title"], payload["cells"])
    print(f"Wrote {len(notebooks)} tutorial notebooks to {NOTEBOOK_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
