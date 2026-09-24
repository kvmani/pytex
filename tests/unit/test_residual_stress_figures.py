"""The canonical sin^2(psi) geometry figure is exactly what its generator writes.

The figure is computed from ``measurement_direction``, the function the
analysis uses, so a committed figure that differs from the generator's output
would be a documentation of a convention the code no longer follows.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "generate_residual_stress_figures.py"


def _generator():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("generate_residual_stress_figures", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_committed_figure_is_the_generators_output() -> None:
    generator = _generator()
    committed = generator.GEOMETRY_FIGURE.read_text(encoding="utf-8").replace("\r\n", "\n")
    assert committed == generator.geometry_svg()


def test_the_figure_is_well_formed_and_accessible() -> None:
    generator = _generator()
    root = ElementTree.fromstring(generator.geometry_svg())
    namespace = "{http://www.w3.org/2000/svg}"
    assert root.tag == f"{namespace}svg"
    assert root.find(f"{namespace}title") is not None
    for marker in root.iter(f"{namespace}marker"):
        # Stroke-scaled markers render several times their intended size.
        assert marker.get("markerUnits") == "userSpaceOnUse"


def test_the_drawn_vector_is_the_analysis_convention() -> None:
    generator = _generator()
    from pytex.diffraction.xrd_residual_stress import measurement_direction

    m = measurement_direction(generator.PHI_DEG, generator.PSI_DEG)
    phi, psi = math.radians(generator.PHI_DEG), math.radians(generator.PSI_DEG)
    assert np.allclose(
        m, [math.cos(phi) * math.sin(psi), math.sin(phi) * math.sin(psi), math.cos(psi)]
    )
    # The azimuth is measured from S1 towards S2 and the tilt from S3.
    assert m[0] > 0.0 and m[1] > 0.0 and m[2] == np.cos(psi)
