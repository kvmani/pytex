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
    """Per-point weights and validity mask for a texture calculation.

    Purpose
    -------
    Not every indexed point deserves equal weight: confidence index, image
    quality, or band contrast all indicate how much a measurement should
    count. This carries the weights together with an optional validity mask,
    and normalizes them at the point of use — checking the length against the
    map view, so weights cannot be misaligned with the points they describe.

    Attributes
    ----------
    weights : np.ndarray
        One weight per point.
    valid_mask : np.ndarray, optional
        Points to exclude entirely, zeroed rather than dropped so indexing
        stays aligned.
    """

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
        """Equal weights for ``count`` orientations.

        The explicit no-weighting choice, preferable to passing ``None`` where the
        weighting decision should be visible in the record.
        """

        if count <= 0:
            raise ValueError("count must be positive.")
        return cls(np.ones(count, dtype=np.float64))

    def for_count(self, count: int) -> np.ndarray:
        """Normalized weights for a map view of exactly ``count`` orientations.

        Purpose
        -------
        Apply the validity mask and normalize to unit sum, checking that the
        weight vector matches the view it is being applied to. The length check
        is the point: silently broadcasting or truncating weights would attach
        the wrong quality value to the wrong measurement point.

        Raises
        ------
        ValueError
            When the length does not match, or when masking leaves no positive
            weight at all — in which case there is no texture to compute.
        """

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
    """The outcome of an EBSD texture workflow run.

    Attributes
    ----------
    crystal_map : CrystalMap
        The map view that was analysed — the phase-selected sub-map when a
        phase was declared.
    report : TextureReport
        The ODF, pole figures, and inverse pole figures.
    segmentation : GrainSegmentation, optional
        Present only when grain segmentation was requested.
    weights : np.ndarray
        The normalized weights actually applied, recorded so the result is
        auditable.
    metadata : dict
        Additional run metadata.
    """

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
        """Number of orientations in the analysed map view.
        """

        return len(self.crystal_map.orientations)

    @property
    def grain_count(self) -> int:
        """Number of segmented grains; zero when segmentation was not requested.
        """

        return 0 if self.segmentation is None else len(self.segmentation.grains)

    def summary(self) -> dict[str, Any]:
        """Compact machine-readable summary of the workflow run.

        Orientation and grain counts, per-phase point counts, the summed weight,
        and any metadata the workflow recorded — enough to audit what was
        analysed without re-running it.
        """

        return {
            "orientation_count": self.orientation_count,
            "grain_count": self.grain_count,
            "phase_summary": self.crystal_map.phase_summary(),
            "weight_sum": float(np.sum(self.weights)),
            **self.metadata,
        }


@dataclass(frozen=True, slots=True)
class EBSDTextureWorkflow:
    """A declared, reproducible EBSD texture-analysis procedure.

    Purpose
    -------
    Holds the analysis choices — which phase, which poles, which specimen
    directions, which weighting, whether to segment grains — as one
    configuration object, so a study's decisions are stated once and applied
    consistently rather than being re-specified at each call site.

    Attributes
    ----------
    phase : optional
        Which phase to analyse; required in effect for a multiphase map.
    poles : tuple
        Planes to produce pole figures for.
    sample_directions : tuple
        Specimen axes for inverse pole figures; ``("x", "y", "z")`` by
        default.
    weights : OrientationQualityWeights, optional
        Quality weighting; uniform when omitted.
    segment_grains : bool
    segmentation_threshold_deg : float
        The grain-boundary criterion, when segmenting.
    Remaining attributes carry the kernel, symmetry, and plotting options
    passed through to the texture report.
    """

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
        """Execute the configured texture workflow on a crystal map.

        Purpose
        -------
        The one-call reproducible pipeline: select the phase, resolve and
        normalize the quality weights, optionally segment grains, and produce the
        texture report — all under one declared configuration, so a study's
        analysis choices live in the workflow object rather than scattered across
        call sites.

        Parameters
        ----------
        crystal_map : CrystalMap
            The map to analyse. When the workflow declares a phase, that phase's
            sub-map is selected first.

        Returns
        -------
        EBSDTextureWorkflowResult
            The texture report, the segmentation when requested, the weights
            actually used, and the summary metadata.
        """

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
