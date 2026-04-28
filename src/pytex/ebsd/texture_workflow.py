from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array
from pytex.core.lattice import CrystalPlane, Phase
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec
from pytex.ebsd.models import CrystalMap, CrystalMapPhase, GrainSegmentation, TextureReport
from pytex.texture import ODF, KernelSpec


@dataclass(frozen=True, slots=True)
class OrientationQualityWeights:
    weights: np.ndarray
    valid_mask: np.ndarray | None = None
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        weights = np.asarray(self.weights, dtype=np.float64).reshape(-1)
        if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("OrientationQualityWeights.weights must be finite and non-negative.")
        if np.isclose(float(np.sum(weights)), 0.0):
            raise ValueError("OrientationQualityWeights.weights must sum to a positive value.")
        weights = np.ascontiguousarray(weights, dtype=np.float64)
        weights.setflags(write=False)
        object.__setattr__(self, "weights", weights)
        if self.valid_mask is not None:
            mask = np.asarray(self.valid_mask, dtype=bool).reshape(-1)
            if mask.shape != weights.shape:
                raise ValueError("valid_mask must match weights shape.")
            mask = np.ascontiguousarray(mask, dtype=bool)
            mask.setflags(write=False)
            object.__setattr__(self, "valid_mask", mask)

    @classmethod
    def uniform(cls, count: int) -> OrientationQualityWeights:
        if count <= 0:
            raise ValueError("count must be positive.")
        return cls(np.ones(count, dtype=np.float64))

    def for_count(self, count: int) -> np.ndarray:
        if count != self.weights.shape[0]:
            raise ValueError("OrientationQualityWeights length must match the selected map view.")
        weights = self.weights.copy()
        if self.valid_mask is not None:
            weights = np.where(self.valid_mask, weights, 0.0)
        total = float(np.sum(weights))
        if total <= 0.0:
            raise ValueError("No valid positive orientation weights remain after masking.")
        normalized = np.ascontiguousarray(weights / total, dtype=np.float64)
        normalized.setflags(write=False)
        return normalized


@dataclass(frozen=True, slots=True)
class EBSDTextureWorkflowResult:
    crystal_map: CrystalMap
    texture_report: TextureReport
    odf: ODF
    weights: np.ndarray
    segmentation: GrainSegmentation | None = None
    experiment_manifest: Any | None = None
    metadata: dict[str, str] = field(default_factory=dict)
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        weights = as_float_array(self.weights, shape=(len(self.crystal_map.orientations),))
        if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("EBSDTextureWorkflowResult.weights must be finite and non-negative.")
        if float(np.sum(weights)) <= 0.0:
            raise ValueError("EBSDTextureWorkflowResult.weights must sum to a positive value.")
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "metadata", {str(k): str(v) for k, v in self.metadata.items()})

    @property
    def orientation_count(self) -> int:
        return len(self.crystal_map.orientations)

    @property
    def grain_count(self) -> int:
        return 0 if self.segmentation is None else len(self.segmentation.grains)

    def summary(self) -> dict[str, Any]:
        return {
            "orientation_count": self.orientation_count,
            "grain_count": self.grain_count,
            "phase_summary": self.crystal_map.phase_summary(),
            "weight_sum": float(np.sum(self.weights)),
            **self.metadata,
        }


@dataclass(frozen=True, slots=True)
class EBSDTextureWorkflow:
    phase: int | str | Phase | CrystalMapPhase | None = None
    poles: tuple[CrystalPlane | ArrayLike, ...] = ()
    sample_directions: tuple[str | ArrayLike, ...] = ("x", "y", "z")
    weights: OrientationQualityWeights | None = None
    kernel: KernelSpec = field(default_factory=KernelSpec)
    specimen_symmetry: SymmetrySpec | None = None
    include_symmetry_family: bool = True
    reduce_by_symmetry: bool = True
    antipodal: bool = True
    segment_grains: bool = False
    segmentation_threshold_deg: float = 5.0
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if (
            not np.isfinite(self.segmentation_threshold_deg)
            or self.segmentation_threshold_deg <= 0.0
        ):
            raise ValueError("segmentation_threshold_deg must be positive and finite.")
        if len(self.sample_directions) == 0:
            raise ValueError("sample_directions must contain at least one direction.")

    def run(self, crystal_map: CrystalMap) -> EBSDTextureWorkflowResult:
        phase_view = crystal_map.select_phase(self.phase) if self.phase is not None else crystal_map
        weight_source = self.weights or OrientationQualityWeights.uniform(
            len(phase_view.orientations)
        )
        weights = weight_source.for_count(len(phase_view.orientations))
        segmentation = (
            phase_view.segment_grains(max_misorientation_deg=self.segmentation_threshold_deg)
            if self.segment_grains
            else None
        )
        report = phase_view.texture_report(
            poles=self.poles,
            sample_directions=self.sample_directions,
            weights=weights,
            kernel=self.kernel,
            specimen_symmetry=self.specimen_symmetry,
            include_symmetry_family=self.include_symmetry_family,
            reduce_by_symmetry=self.reduce_by_symmetry,
            antipodal=self.antipodal,
            provenance=self.provenance or phase_view.provenance,
        )
        manifest = phase_view.to_experiment_manifest(
            source_system="pytex.ebsd_texture_workflow",
            metadata={"workflow": "ebsd_texture"},
        )
        return EBSDTextureWorkflowResult(
            crystal_map=phase_view,
            texture_report=report,
            odf=report.odf,
            weights=weights,
            segmentation=segmentation,
            experiment_manifest=manifest,
            metadata={
                "workflow": "ebsd_texture",
                "weighted": str(self.weights is not None).lower(),
            },
            provenance=self.provenance or phase_view.provenance,
        )


__all__ = [
    "EBSDTextureWorkflow",
    "EBSDTextureWorkflowResult",
    "OrientationQualityWeights",
]
