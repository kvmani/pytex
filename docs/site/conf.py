from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"

sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(REPO_ROOT))

project = "PyTex"
author = "PyTex contributors"
release = "0.1.0.dev0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinxcontrib.mermaid",
]

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = ["_build", "README.md"]

html_theme = "furo"
html_title = "PyTex"
html_static_path = []

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

root_doc = "index"
