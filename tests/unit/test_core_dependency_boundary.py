"""Stable imports must remain independent of experimental research modules."""

import subprocess
import sys
from pathlib import Path


def test_importing_core_does_not_load_experimental_modules() -> None:
    source = Path(__file__).resolve().parents[2] / "src"
    script = (
        f"import sys; sys.path.insert(0, {str(source)!r}); import pytex.core; "
        "assert not any(name.startswith('pytex.experimental') for name in sys.modules)"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_experimental_facade_uses_the_shared_implementation() -> None:
    from pytex.core._parent_scoring import ParentReconstructionResult, score_parent_orientations
    from pytex.experimental import phase_transformation

    assert phase_transformation.score_parent_orientations is score_parent_orientations
    assert phase_transformation.ParentReconstructionResult is ParentReconstructionResult
