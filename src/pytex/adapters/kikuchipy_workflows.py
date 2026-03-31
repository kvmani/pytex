from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.adapters.ebsd import (
    EBSDImportManifest,
    NormalizedEBSDDataset,
    normalize_kikuchipy_dataset,
)
from pytex.adapters.manifests import ExperimentManifest
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase
from pytex.core.orientation import OrientationSet
from pytex.core.provenance import ProvenanceRecord
from pytex.ebsd import CrystalMap


@dataclass(frozen=True, slots=True)
class KikuchiPyWorkflowResult:
    dataset: NormalizedEBSDDataset
    experiment_manifest: ExperimentManifest
    xmap: Any | None = None
    index_data: Any | None = None
    band_data: Any | None = None


def _resolve_frames(
    *,
    frames: tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]
    | dict[str, ReferenceFrame]
    | None = None,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
) -> tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]:
    resolved_crystal = crystal_frame
    resolved_specimen = specimen_frame
    resolved_map = map_frame
    if frames is not None:
        if isinstance(frames, dict):
            resolved_crystal = frames.get("crystal_frame", resolved_crystal)
            resolved_specimen = frames.get("specimen_frame", resolved_specimen)
            resolved_map = frames.get("map_frame", resolved_map)
        else:
            if len(frames) != 3:
                raise ValueError("frames must contain crystal, specimen, and map frames.")
            resolved_crystal, resolved_specimen, resolved_map = frames
    if resolved_crystal is None or resolved_specimen is None or resolved_map is None:
        raise ValueError(
            "crystal_frame, specimen_frame, and map_frame are required unless frames supplies them."
        )
    return resolved_crystal, resolved_specimen, resolved_map


def _first_string_attribute(obj: Any, names: tuple[str, ...]) -> str | None:
    for name in names:
        value = getattr(obj, name, None)
        if isinstance(value, str) and value:
            return value
    return None


def _coordinates_from_xmap(
    xmap: Any,
) -> tuple[np.ndarray, tuple[int, int] | None, tuple[float, float] | None]:
    x = getattr(xmap, "x", None)
    y = getattr(xmap, "y", None)
    if x is None and hasattr(xmap, "col"):
        x = xmap.col
    if y is None and hasattr(xmap, "row"):
        y = xmap.row
    if x is None:
        raise ValueError("Unable to extract x coordinates from the KikuchiPy crystal map.")
    x_array = np.asarray(x, dtype=np.float64).reshape(-1)
    if y is None:
        y_array = np.zeros_like(x_array)
    else:
        y_array = np.asarray(y, dtype=np.float64).reshape(-1)
    if x_array.shape != y_array.shape:
        raise ValueError("KikuchiPy crystal-map x and y coordinates must have matching lengths.")
    grid_shape_raw = tuple(int(value) for value in np.atleast_1d(getattr(xmap, "shape", ())))
    grid_shape = grid_shape_raw if len(grid_shape_raw) == 2 else None
    dx = getattr(xmap, "dx", None)
    dy = getattr(xmap, "dy", None)
    step_sizes = None
    if grid_shape is not None and dx is not None and dy is not None:
        step_sizes = (float(dx), float(dy))
    return np.column_stack([x_array, y_array]), grid_shape, step_sizes


def _phase_metadata_from_inputs(
    xmap: Any,
    phase: Phase | None,
) -> tuple[str, str]:
    if phase is not None:
        return phase.name, phase.symmetry.point_group
    phase_name = _first_string_attribute(xmap, ("phase_name",))
    point_group = _first_string_attribute(xmap, ("point_group",))
    if phase_name is None or point_group is None:
        raise ValueError(
            "normalize_ebsd() requires a PyTex phase or a KikuchiPy/orix crystal map exposing "
            "phase_name and point_group metadata."
        )
    return phase_name, point_group


def _dataset_from_xmap(
    xmap: Any,
    *,
    crystal_frame: ReferenceFrame,
    specimen_frame: ReferenceFrame,
    map_frame: ReferenceFrame,
    phase: Phase | None = None,
    source_file: str | None = None,
    metadata: dict[str, str] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> NormalizedEBSDDataset:
    rotations = getattr(xmap, "rotations", None)
    if rotations is None:
        raise ValueError("KikuchiPy crystal maps must expose a rotations attribute.")
    quaternions = np.asarray(rotations, dtype=np.float64)
    if quaternions.ndim != 2 or quaternions.shape[1] != 4:
        raise ValueError("KikuchiPy crystal-map rotations must have quaternion shape (n, 4).")
    coordinates, grid_shape, step_sizes = _coordinates_from_xmap(xmap)
    phase_name, point_group = _phase_metadata_from_inputs(xmap, phase)
    resolved_provenance = provenance or ProvenanceRecord.minimal(
        "kikuchipy",
        note="Normalized from a KikuchiPy/orix crystal map.",
    )
    orientations = OrientationSet.from_quaternions(
        quaternions,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        phase=phase,
        provenance=resolved_provenance,
    )
    crystal_map = CrystalMap(
        coordinates=coordinates,
        orientations=orientations,
        map_frame=map_frame,
        grid_shape=grid_shape,
        step_sizes=step_sizes,
        provenance=resolved_provenance,
    )
    manifest = EBSDImportManifest(
        source_system="kikuchipy",
        source_file=source_file
        or _first_string_attribute(xmap, ("source_file", "path", "filename"))
        or "kikuchipy_in_memory",
        phase_name=phase_name,
        point_group=point_group,
        crystal_frame={
            "name": crystal_frame.name,
            "axis_1": crystal_frame.axes[0],
            "axis_2": crystal_frame.axes[1],
            "axis_3": crystal_frame.axes[2],
        },
        specimen_frame={
            "name": specimen_frame.name,
            "axis_1": specimen_frame.axes[0],
            "axis_2": specimen_frame.axes[1],
            "axis_3": specimen_frame.axes[2],
        },
        metadata=dict(metadata or {}),
    )
    return NormalizedEBSDDataset(crystal_map=crystal_map, manifest=manifest)


def _experiment_manifest_for_dataset(dataset: NormalizedEBSDDataset) -> ExperimentManifest:
    referenced_files = (
        ()
        if dataset.manifest.source_file.endswith("_in_memory")
        else (dataset.manifest.source_file,)
    )
    return dataset.crystal_map.to_experiment_manifest(
        source_system=dataset.manifest.source_system,
        referenced_files=referenced_files,
        metadata=dict(dataset.manifest.metadata),
    )


def _dataset_with_overrides(
    dataset: NormalizedEBSDDataset,
    *,
    phase: Phase | None = None,
    source_file: str | None = None,
    metadata: dict[str, str] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> NormalizedEBSDDataset:
    resolved_provenance = provenance or dataset.crystal_map.provenance
    resolved_phase = phase or dataset.crystal_map.orientations.phase
    if resolved_phase is not None:
        orientations = OrientationSet.from_quaternions(
            dataset.crystal_map.orientations.quaternions,
            crystal_frame=dataset.crystal_map.orientations.crystal_frame,
            specimen_frame=dataset.crystal_map.orientations.specimen_frame,
            phase=resolved_phase,
            provenance=resolved_provenance,
        )
    else:
        orientations = OrientationSet.from_quaternions(
            dataset.crystal_map.orientations.quaternions,
            crystal_frame=dataset.crystal_map.orientations.crystal_frame,
            specimen_frame=dataset.crystal_map.orientations.specimen_frame,
            symmetry=dataset.crystal_map.orientations.symmetry,
            provenance=resolved_provenance,
        )
    crystal_map = CrystalMap(
        coordinates=dataset.crystal_map.coordinates,
        orientations=orientations,
        map_frame=dataset.crystal_map.map_frame,
        grid_shape=dataset.crystal_map.grid_shape,
        step_sizes=dataset.crystal_map.step_sizes,
        acquisition_geometry=dataset.crystal_map.acquisition_geometry,
        calibration_record=dataset.crystal_map.calibration_record,
        measurement_quality=dataset.crystal_map.measurement_quality,
        provenance=resolved_provenance,
    )
    manifest = EBSDImportManifest(
        source_system=dataset.manifest.source_system,
        source_file=source_file or dataset.manifest.source_file,
        phase_name=(
            resolved_phase.name
            if resolved_phase is not None
            else dataset.manifest.phase_name
        ),
        point_group=(
            resolved_phase.symmetry.point_group
            if resolved_phase is not None
            else dataset.manifest.point_group
        ),
        orientation_convention=dataset.manifest.orientation_convention,
        angle_unit=dataset.manifest.angle_unit,
        crystal_frame=dict(dataset.manifest.crystal_frame),
        specimen_frame=dict(dataset.manifest.specimen_frame),
        metadata={**dataset.manifest.metadata, **(metadata or {})},
    )
    return NormalizedEBSDDataset(crystal_map=crystal_map, manifest=manifest)


def normalize_ebsd(
    ebsd_signal: Any,
    *,
    frames: tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]
    | dict[str, ReferenceFrame]
    | None = None,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    source_file: str | None = None,
    metadata: dict[str, str] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> NormalizedEBSDDataset:
    resolved_frames = _resolve_frames(
        frames=frames,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
    )
    crystal_frame, specimen_frame, map_frame = resolved_frames
    xmap = getattr(ebsd_signal, "xmap", ebsd_signal)
    if hasattr(xmap, "rotations") and (hasattr(xmap, "x") or hasattr(xmap, "col")):
        return _dataset_from_xmap(
            xmap,
            crystal_frame=crystal_frame,
            specimen_frame=specimen_frame,
            map_frame=map_frame,
            phase=phase,
            source_file=source_file,
            metadata=metadata,
            provenance=provenance,
        )
    dataset = normalize_kikuchipy_dataset(
        ebsd_signal,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
    )
    if phase is None and metadata is None and source_file is None and provenance is None:
        return dataset
    return _dataset_with_overrides(
        dataset,
        phase=phase,
        source_file=source_file,
        metadata=metadata,
        provenance=provenance,
    )


def index_hough(
    ebsd_signal: Any,
    *,
    phase_list: Any,
    indexer: Any,
    frames: tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]
    | dict[str, ReferenceFrame]
    | None = None,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    source_file: str | None = None,
    metadata: dict[str, str] | None = None,
    provenance: ProvenanceRecord | None = None,
    **hough_kwargs: Any,
) -> KikuchiPyWorkflowResult:
    result = ebsd_signal.hough_indexing(phase_list=phase_list, indexer=indexer, **hough_kwargs)
    if isinstance(result, tuple):
        xmap = result[0]
        index_data = result[1] if len(result) > 1 else None
        band_data = result[2] if len(result) > 2 else None
    else:
        xmap = result
        index_data = None
        band_data = None
    dataset = normalize_ebsd(
        xmap,
        frames=frames,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
        phase=phase,
        source_file=source_file,
        metadata={**(metadata or {}), "workflow": "hough_indexing"},
        provenance=provenance,
    )
    return KikuchiPyWorkflowResult(
        dataset=dataset,
        experiment_manifest=_experiment_manifest_for_dataset(dataset),
        xmap=xmap,
        index_data=index_data,
        band_data=band_data,
    )


def refine_orientations(
    ebsd_signal: Any,
    xmap: Any | KikuchiPyWorkflowResult,
    *,
    frames: tuple[ReferenceFrame, ReferenceFrame, ReferenceFrame]
    | dict[str, ReferenceFrame]
    | None = None,
    crystal_frame: ReferenceFrame | None = None,
    specimen_frame: ReferenceFrame | None = None,
    map_frame: ReferenceFrame | None = None,
    phase: Phase | None = None,
    source_file: str | None = None,
    metadata: dict[str, str] | None = None,
    provenance: ProvenanceRecord | None = None,
    **refine_kwargs: Any,
) -> KikuchiPyWorkflowResult:
    seed_xmap = xmap.xmap if isinstance(xmap, KikuchiPyWorkflowResult) else xmap
    refined_xmap = ebsd_signal.refine_orientation(seed_xmap, **refine_kwargs)
    dataset = normalize_ebsd(
        refined_xmap,
        frames=frames,
        crystal_frame=crystal_frame,
        specimen_frame=specimen_frame,
        map_frame=map_frame,
        phase=phase,
        source_file=source_file,
        metadata={**(metadata or {}), "workflow": "orientation_refinement"},
        provenance=provenance,
    )
    return KikuchiPyWorkflowResult(
        dataset=dataset,
        experiment_manifest=_experiment_manifest_for_dataset(dataset),
        xmap=refined_xmap,
    )
