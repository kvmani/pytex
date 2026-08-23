from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"

sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(REPO_ROOT))

project = "PyTex"
author = "PyTex contributors"
# Read, not restated. A hard-coded release here is a second version literal
# that nothing checks, so it drifts silently the first time the package is
# bumped and every built page then carries the wrong number.
from pytex._version import __version__  # noqa: E402

release = __version__

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "special-members": "__len__, __getitem__",
    "show-inheritance": True,
    "member-order": "bysource",
    # The exhaustive reference is the only autodoc consumer. It renders the
    # complete surface, while the curated API guide owns searchable index
    # entries; registering both produced hundreds of duplicate-object warnings.
    "no-index": True,
}
autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_typehints_format = "short"
autodoc_preserve_defaults = True
napoleon_numpy_docstring = True
napoleon_google_docstring = True
templates_path = ["_templates"]
exclude_patterns = ["_build", "README.md"]

html_theme = "furo"
html_title = "PyTex"
html_static_path = ["_static"]
html_css_files = ["architecture.css"]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
]

# Headings down to the third level get anchors, so a page can link to a *section*
# of another page. Without this every `other_page.md#some-section` link resolves
# to the page and warns about the fragment, which is how the TEM indexing page
# came to point at a section of the workbench guide that the build could not
# find even though the heading was right there.
myst_heading_anchors = 3

# Notebooks are committed *without* outputs, so the site must produce them at
# build time rather than render stored ones. "cache" executes each notebook once
# and reuses the result until its source changes, which keeps incremental builds
# cheap. Raising on error means a tutorial that no longer runs fails the docs
# build instead of publishing a traceback as if it were a result.
nb_execution_mode = "cache"
nb_execution_raise_on_error = True
nb_execution_timeout = 300

root_doc = "index"
