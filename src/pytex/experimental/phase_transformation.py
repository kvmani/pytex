"""Experimental parent-scoring facade over the shared core implementation.

The import path remains available for existing experiments. Sharing the implementation
with the stable bounded-candidate report does not promote map-scale reconstruction.
"""

from pytex.core._parent_scoring import (
    ParentReconstructionResult,
    score_parent_orientations,
)
from pytex.core._parent_scoring import (
    ReductionMode as ReductionMode,
)

__all__ = ["ParentReconstructionResult", "score_parent_orientations"]
