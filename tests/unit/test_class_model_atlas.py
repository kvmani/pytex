"""Guards for the Class & Object Model Atlas.

Three things have to stay true for the atlas to be worth trusting:

1. The model is read from the source, not written down. If ``get_type_hints``
   silently fails for a module, the relations it declares vanish from every
   diagram with no visible error — so an unresolved annotation is a test failure,
   not a warning.
2. The committed SVGs are what the generator currently produces. A figure that
   has drifted from the code is worse than no figure, because a reader cannot
   tell the difference.
3. The renderer's layout is deterministic. Byte-comparison in (2) is only
   possible if two runs agree.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from class_model import build_model, package_flow, short_type_name  # noqa: E402
from generate_class_model_figures import VIEWS, build_figures  # noqa: E402

from pytex.plotting.class_diagrams import ClassBox, ClassEdge, class_diagram_svg  # noqa: E402


@pytest.fixture(scope="module")
def model():  # type: ignore[no-untyped-def]
    return build_model()


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------


def test_every_annotation_resolves(model) -> None:  # type: ignore[no-untyped-def]
    """An unresolved annotation silently deletes relations from the atlas.

    Several modules import the types they reference under ``TYPE_CHECKING`` to
    break an import cycle, so the names are absent at runtime. The model resolves
    those through the library's own class index; if that ever stops working the
    affected classes lose their edges without any other symptom.
    """

    assert model.unresolved == ()


def test_model_reads_the_declared_core_composition(model) -> None:  # type: ignore[no-untyped-def]
    """`Phase` composes a lattice, a symmetry, and a crystal frame — all required."""

    relations = {
        (rel.target, rel.label): rel
        for rel in model.relations
        if rel.source == "pytex.core.lattice.Phase"
    }
    for target, label in (
        ("pytex.core.lattice.Lattice", "lattice"),
        ("pytex.core.symmetry.SymmetrySpec", "symmetry"),
        ("pytex.core.frames.ReferenceFrame", "crystal_frame"),
    ):
        assert (target, label) in relations, f"Phase.{label} is missing from the model"
        assert relations[(target, label)].kind == "composition"
        assert relations[(target, label)].multiplicity == "1"

    optional = relations[("pytex.core.lattice.UnitCell", "unit_cell")]
    assert optional.kind == "association"
    assert optional.multiplicity == "0..1"


def test_sequence_fields_are_recorded_as_many(model) -> None:
    """A crystal map holds a *sequence* of phases, and the model must say so."""

    entries = [
        rel
        for rel in model.relations
        if rel.source == "pytex.ebsd.models.CrystalMap"
        and rel.target == "pytex.ebsd.models.CrystalMapPhase"
    ]
    assert entries, "CrystalMap.phase_entries is missing from the model"
    assert all(rel.multiplicity == "*" for rel in entries)


def test_inheritance_is_read_from_the_bases(model) -> None:
    """The tilt-envelope family is the library's one multi-subclass hierarchy."""

    subclasses = {
        rel.target
        for rel in model.of_kind("inheritance")
        if rel.source == "pytex.tem.stage.TiltEnvelope"
    }
    assert subclasses == {
        "pytex.tem.stage.EllipticalEnvelope",
        "pytex.tem.stage.MaskedEnvelope",
        "pytex.tem.stage.PolygonEnvelope",
        "pytex.tem.stage.RectangularEnvelope",
    }


def test_the_library_is_composition_first(model) -> None:
    """The claim the atlas makes in prose, pinned as a measurement.

    If PyTex ever grows a deep hierarchy this fails, and the atlas pages that say
    "PyTex composes rather than subclasses" must be rewritten rather than left to
    quietly become false.
    """

    dataclasses = sum(1 for entry in model.entries.values() if entry.kind == "dataclass")
    assert dataclasses > 0.85 * len(model.entries)
    assert len(model.of_kind("inheritance")) < 0.1 * len(model.entries)
    assert len(model.of_kind("composition")) > 100


def test_no_domain_package_is_referenced_by_the_core(model) -> None:
    """The layering claim on the architecture figure, checked rather than asserted."""

    offenders = [
        (source, target, count) for source, target, count in package_flow(model) if source == "core"
    ]
    assert not offenders, f"pytex.core holds typed references to domain packages: {offenders}"


@pytest.mark.parametrize(
    ("annotation", "expected"),
    [
        (int, "int"),
        (tuple[str, ...], "str[]"),
        (str | None, "str?"),
        (dict[str, int], "{str: int}"),
        (tuple[float, float], "tuple[float, float]"),
    ],
)
def test_short_type_names(annotation: object, expected: str) -> None:
    """Card lines must stay readable, so annotations are rendered compactly."""

    assert short_type_name(annotation) == expected


# ---------------------------------------------------------------------------
# The renderer
# ---------------------------------------------------------------------------


def _diagram(edges: tuple[ClassEdge, ...]) -> str:
    boxes = [
        ClassBox(key="A", name="Owner", module="core.a", fields=("part: Part",)),
        ClassBox(key="B", name="Part", module="core.b", fields=("value: float",)),
        ClassBox(key="C", name="Other", module="core.c"),
    ]
    return class_diagram_svg(
        title="Test",
        subtitle="subtitle",
        description="description",
        boxes=boxes,
        edges=edges,
    )


def test_owner_is_drawn_above_the_part() -> None:
    """Layering is the diagram's main statement: what depends on what."""

    svg = _diagram((ClassEdge(source="A", target="B", label="part"),))
    owner = re.search(r'<text x="[0-9.]+" y="([0-9.]+)"[^>]*>Owner</text>', svg)
    part = re.search(r'<text x="[0-9.]+" y="([0-9.]+)"[^>]*>Part</text>', svg)
    assert owner is not None and part is not None
    assert float(owner.group(1)) < float(part.group(1))


def test_diagram_is_deterministic() -> None:
    """Two runs must agree, or the committed figures cannot be compared."""

    edges = (
        ClassEdge(source="A", target="B", label="part"),
        ClassEdge(source="C", target="B", kind="association", label="other"),
    )
    assert _diagram(edges) == _diagram(edges)


def test_diagram_carries_the_accessible_heading() -> None:
    """The style guide requires a title and description on every canonical SVG."""

    svg = _diagram((ClassEdge(source="A", target="B"),))
    assert "<title>Test</title>" in svg
    assert "<desc>description</desc>" in svg
    assert 'markerUnits="userSpaceOnUse"' in svg


def test_edges_naming_absent_boxes_are_ignored() -> None:
    """Callers may pass a whole domain's relations with a filtered card list."""

    svg = _diagram((ClassEdge(source="A", target="ZZZ", label="ghost"),))
    assert ">ghost<" not in svg


def test_empty_diagram_is_rejected() -> None:
    with pytest.raises(ValueError):
        class_diagram_svg(title="t", subtitle="s", description="d", boxes=(), edges=())


# ---------------------------------------------------------------------------
# The committed assets
# ---------------------------------------------------------------------------

#: Every figure the atlas generator writes.
ATLAS_FIGURES = (
    "class_model_architecture.svg",
    "class_hierarchy.svg",
    *(f"class_model_{spec.key}.svg" for spec in VIEWS),
)


@pytest.fixture(scope="module")
def figures() -> dict[str, str]:
    """Every atlas figure, rendered once for the whole module."""

    return build_figures()


def test_generator_writes_exactly_the_tracked_figures(figures: dict[str, str]) -> None:
    assert set(figures) == set(ATLAS_FIGURES)


@pytest.mark.parametrize("name", ATLAS_FIGURES)
def test_committed_figures_match_the_generator(name: str, figures: dict[str, str]) -> None:
    """A generated figure that has drifted from the code is a defect."""

    committed = (FIGURES_DIR / name).read_text(encoding="utf-8")
    assert figures[name] == committed, (
        f"{name} differs from the committed asset. Run "
        "`python scripts/generate_class_model_figures.py`."
    )


ATLAS_PAGE = REPO_ROOT / "docs" / "site" / "architecture" / "class_model_atlas.md"


def test_the_page_does_not_hand_transcribe_the_model_counts(model) -> None:  # type: ignore[no-untyped-def]
    """The three numbers the page states must be the three the model reports.

    Documentation numbers are not allowed to be typed in from memory here any
    more than anywhere else in PyTex. This is the one place the atlas page names
    figures in prose, so it is the one place that can go stale.
    """

    page = ATLAS_PAGE.read_text(encoding="utf-8")
    match = re.search(
        r"Of (\d+) public classes, (\d+) are dataclasses and only (\d+)\s+inheritance relations",
        page,
    )
    assert match is not None, "the atlas page no longer states the counts in the expected form"
    classes, dataclasses_stated, inheritance = (int(value) for value in match.groups())

    assert classes == len(model.entries)
    assert dataclasses_stated == sum(
        1 for entry in model.entries.values() if entry.kind == "dataclass"
    )
    assert inheritance == len(model.of_kind("inheritance"))


def test_the_atlas_covers_every_scientific_domain() -> None:
    """The atlas is per-domain by design; one enormous graph was the thing to avoid."""

    keys = {spec.key for spec in VIEWS}
    assert {"core", "texture", "ebsd", "diffraction"} <= keys
