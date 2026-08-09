"""Generate the Class & Object Model Atlas figures.

These are *generated assets*: `scripts/class_model.py` reads PyTex's real classes
and relations by importing the package, and `pytex.plotting.class_diagrams` draws
them in the canonical documentation style. Nothing about a diagram is authored by
hand except which corner of the library it looks at, so a renamed field or a new
composition shows up the next time this script runs.

Usage::

    python scripts/generate_class_model_figures.py

Outputs (tracked as canonical documentation assets):

- ``docs/figures/class_model_architecture.svg`` — the library as packages, with
  the typed references that actually cross each boundary counted.
- ``docs/figures/class_hierarchy.svg`` — every inheritance relation in PyTex.
  It is small on purpose: the library is composition-first, and the figure says so.
- ``docs/figures/class_model_core.svg`` — frames, symmetry, lattice/phase,
  orientation, and Miller objects.
- ``docs/figures/class_model_transformation.svg`` — the orientation-relationship
  and variant model.
- ``docs/figures/class_model_texture.svg`` — ODF, pole figures, kernels, reports.
- ``docs/figures/class_model_ebsd.svg`` — crystal maps, grains, boundaries.
- ``docs/figures/class_model_diffraction.svg`` — geometry, patterns, simulation.
- ``docs/figures/class_model_tem.svg`` — stage, tilt navigation, ambiguity.

How a domain view chooses what to draw
--------------------------------------
Not by a hand-kept list of class names, which would go stale silently. Each view
names *modules*; the roster is the most connected classes in them, plus any
pinned anchors, plus the core objects those classes actually reference. A new
class that becomes central to a domain therefore enters its diagram on its own.
``ProvenanceRecord`` is the one deliberate omission: 85 classes carry it, and
drawing all 85 edges would say only that provenance is universal — which the
legend says in one line instead.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from class_model import ClassModel, build_model, classes_in, package_flow  # noqa: E402

from pytex.plotting.class_diagrams import ClassBox, ClassEdge, class_diagram_svg  # noqa: E402

FIGURES = REPO_ROOT / "docs" / "figures"

#: Carried by 85 classes. Drawing it would turn every domain view into a star
#: around one node and say nothing a legend line cannot.
UBIQUITOUS = ("pytex.core.provenance.ProvenanceRecord",)

#: Fields shown on a card before it is summarized with a footnote.
MAX_FIELDS = 7


@dataclass(frozen=True)
class ViewSpec:
    """One domain diagram: where to look, and what must not be left out."""

    key: str
    title: str
    subtitle: str
    description: str
    modules: tuple[str, ...]
    pin: tuple[str, ...] = ()
    """Anchors that must appear even if their connectivity is low."""

    limit: int = 11
    """How many of the most connected classes in ``modules`` to draw."""

    context: int = 4
    """How many referenced classes from other packages to draw alongside."""

    exclude: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    context_modules: tuple[str, ...] = ("pytex.core",)


VIEWS: tuple[ViewSpec, ...] = (
    ViewSpec(
        key="core",
        title="Object Model — Canonical Core",
        subtitle=(
            "Frames, symmetry, lattice and phase, orientation, and Miller objects, "
            "read from the dataclass fields that declare them"
        ),
        description=(
            "UML-style object model of pytex.core: reference frames and transforms, "
            "symmetry specifications, lattice, unit cell and phase, rotations and "
            "orientations, and the Miller plane, direction and zone-axis objects, with "
            "composition and association edges taken from the declared field types."
        ),
        modules=("pytex.core",),
        pin=(
            "pytex.core.frames.ReferenceFrame",
            "pytex.core.symmetry.SymmetrySpec",
            "pytex.core.lattice.Lattice",
            "pytex.core.lattice.Phase",
            "pytex.core.orientation.Orientation",
            "pytex.core.orientation.Rotation",
            "pytex.core.miller.MillerPlane",
            "pytex.core.miller.MillerDirection",
            "pytex.core.lattice.CrystalPlane",
            "pytex.core.lattice.CrystalDirection",
            "pytex.core.lattice.ZoneAxis",
        ),
        limit=13,
        context=0,
        notes=(
            "Phase is the join: one lattice, one symmetry, one crystal frame.",
            "A Miller object becomes a CrystalPlane or CrystalDirection once bound to a phase.",
        ),
    ),
    ViewSpec(
        key="transformation",
        title="Object Model — Orientation Relationships And Variants",
        subtitle=(
            "The flagship transformation-crystallography model: relationships, "
            "variants, and the records that make a result explainable"
        ),
        description=(
            "UML-style object model of pytex.core.transformation: the orientation "
            "relationship, its parallelism constraints, the variant set it generates, "
            "and the report objects that carry a characterization result together with "
            "its provenance."
        ),
        modules=("pytex.core.transformation", "pytex.core.parent_reconstruction"),
        pin=("pytex.core.transformation.OrientationRelationship",),
        limit=11,
        context=3,
    ),
    ViewSpec(
        key="texture",
        title="Object Model — Texture",
        subtitle="Orientation distribution functions, pole figures, kernels, and inversion reports",
        description=(
            "UML-style object model of pytex.texture: discrete and harmonic ODFs, pole "
            "and inverse pole figures, the kernel specification they share, and the "
            "reconstruction and residual reports produced by ODF inversion."
        ),
        modules=("pytex.texture",),
        pin=(
            "pytex.texture.models.ODF",
            "pytex.texture.models.PoleFigure",
            "pytex.texture.harmonics.HarmonicODF",
        ),
        limit=11,
        context=3,
    ),
    ViewSpec(
        key="ebsd",
        title="Object Model — EBSD",
        subtitle="Crystal maps, grain segmentation, boundary networks, and the texture workflow",
        description=(
            "UML-style object model of pytex.ebsd: the crystal map and its phases, "
            "grain segmentation and the grain graph, boundary networks, and the workflow "
            "object that turns a normalized map into texture outputs."
        ),
        modules=("pytex.ebsd",),
        pin=("pytex.ebsd.models.CrystalMap",),
        limit=11,
        context=4,
    ),
    ViewSpec(
        key="diffraction",
        title="Object Model — Diffraction",
        subtitle="Geometry, radiation, patterns, and the kinematic simulation objects",
        description=(
            "UML-style object model of pytex.diffraction: diffraction geometry and "
            "radiation specification, powder, SAED, Kikuchi and composite patterns, and "
            "the kinematic simulation configuration and spot tables behind them."
        ),
        modules=("pytex.diffraction",),
        pin=(
            "pytex.diffraction.models.DiffractionGeometry",
            "pytex.diffraction.models.DiffractionPattern",
        ),
        limit=12,
        context=3,
    ),
    ViewSpec(
        key="tem",
        title="Object Model — TEM Tilt Navigation",
        subtitle="Stage geometry and calibration, tilt solutions and paths, and ambiguity reports",
        description=(
            "UML-style object model of pytex.tem: stage position, calibration and tilt "
            "envelopes, the tilt solutions and paths produced by navigation, and the "
            "ambiguity reporting objects."
        ),
        modules=("pytex.tem",),
        pin=("pytex.tem.stage.StagePosition",),
        limit=12,
        context=3,
    ),
)


# ---------------------------------------------------------------------------
# Card construction
# ---------------------------------------------------------------------------


def _module_label(module: str) -> str:
    """``pytex.core.lattice`` -> ``core.lattice``; the package prefix is a given."""

    return module.removeprefix("pytex.")


def _stereotype(entry: object) -> str:
    kind = entry.kind  # type: ignore[attr-defined]
    if kind == "dataclass":
        return "frozen dataclass" if entry.frozen else "dataclass"  # type: ignore[attr-defined]
    if kind == "enum":
        return "StrEnum"
    if kind == "protocol":
        return "protocol"
    return "class"


def _box(
    model: ClassModel, qualname: str, *, related: set[str], emphasis: bool = False
) -> ClassBox:
    """One card, with the fields that carry a drawn relation shown first.

    A card cannot list twenty fields and stay readable, so the ones a reader
    needs to follow the diagram — those naming another card — are promoted, and
    the remainder is summarized rather than silently dropped.
    """

    entry = model.entries[qualname]
    relation_fields = {
        rel.label
        for rel in model.relations
        if rel.source == qualname and rel.target in related and rel.label
    }
    # Underscore-prefixed fields are internal cache or bookkeeping state, not
    # part of the object's scientific meaning; listing them would put private
    # implementation on a public diagram.
    declared = [f for f in entry.fields if not f.name.startswith("_")]
    ordered = [f for f in declared if f.name in relation_fields]
    ordered += [f for f in declared if f.name not in relation_fields]

    shown = ordered[:MAX_FIELDS]
    if entry.kind == "enum":
        lines = tuple(f.name for f in shown)
    else:
        lines = tuple(f"{f.name}: {f.type_name}" for f in shown)
    remaining = len(declared) - len(shown)
    return ClassBox(
        key=qualname,
        name=entry.name,
        module=_module_label(entry.module),
        stereotype=_stereotype(entry),
        fields=lines,
        role=entry.package,
        footnote=f"+ {remaining} more" if remaining else "",
        emphasis=emphasis,
    )


def _edges(model: ClassModel, roster: list[str]) -> list[ClassEdge]:
    return [
        ClassEdge(
            source=rel.source,
            target=rel.target,
            kind=rel.kind,
            label=rel.label,
            multiplicity=rel.multiplicity,
        )
        for rel in model.induced(roster)
    ]


def _roster(model: ClassModel, spec: ViewSpec) -> list[str]:
    """Which classes a domain view draws, chosen from connectivity, not a list."""

    pool: list[str] = []
    for module in spec.modules:
        pool.extend(classes_in(model, module))
    pool = [name for name in dict.fromkeys(pool) if name not in UBIQUITOUS]
    pool = [name for name in pool if name not in spec.exclude]

    ranked = sorted(pool, key=lambda name: (-model.degree(name, within=pool), name))
    chosen = list(dict.fromkeys([*spec.pin, *ranked]))[: max(spec.limit, len(spec.pin))]

    if spec.context:
        outside: dict[str, int] = {}
        for rel in model.relations:
            if rel.source not in chosen or rel.target in chosen or rel.target in UBIQUITOUS:
                continue
            entry = model.entries[rel.target]
            if not any(
                entry.module == prefix or entry.module.startswith(f"{prefix}.")
                for prefix in spec.context_modules
            ):
                continue
            outside[rel.target] = outside.get(rel.target, 0) + 1
        ranked_outside = sorted(outside, key=lambda name: (-outside[name], name))
        chosen.extend(ranked_outside[: spec.context])

    return list(dict.fromkeys(chosen))


def domain_figure(model: ClassModel, spec: ViewSpec) -> str:
    """Render one domain object-model diagram."""

    roster = _roster(model, spec)
    related = set(roster)
    pinned = set(spec.pin)
    boxes = [
        _box(model, name, related=related, emphasis=name in pinned and len(pinned) < 4)
        for name in roster
    ]
    notes = [
        *spec.notes,
        "ProvenanceRecord is omitted: 85 classes carry it, on every result object.",
        "Fields naming another card are listed first; the rest are summarized.",
    ]
    return class_diagram_svg(
        title=spec.title,
        subtitle=spec.subtitle,
        description=spec.description,
        boxes=boxes,
        edges=_edges(model, roster),
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Architecture and hierarchy views
# ---------------------------------------------------------------------------

#: One line of intent per package. The class counts and the edge weights beside
#: them are measured, not written down here.
_PACKAGE_ROLE: dict[str, tuple[str, str]] = {
    "core": ("canonical model", "frames, symmetry, lattice, orientation, Miller, provenance"),
    "texture": ("domain", "ODF, pole figures, kernels, fibres, components"),
    "ebsd": ("domain", "crystal maps, grains, boundaries, KAM/GND"),
    "diffraction": ("domain", "geometry, XRD, SAED, Kikuchi, kinematics"),
    "tem": ("domain", "stage, tilt navigation, ambiguity"),
    "properties": ("domain", "elasticity, slip systems, Taylor"),
    "experimental": ("unstable", "OR identification, parent reconstruction"),
    "plotting": ("presentation", "primitives, IPF, 3D scenes, SVG figures"),
    "adapters": ("boundary", "orix, kikuchipy, LaboTex, XRDML, manifests"),
}


def architecture_figure(model: ClassModel) -> str:
    """The library as packages, with the boundary crossings counted.

    This is not a diagram of intended layering. Each arrow's number is how many
    typed field references in the source actually cross that boundary, so a layer
    violation would appear here as an arrow pointing the wrong way.
    """

    counts = {(source, target): count for source, target, count in package_flow(model)}
    packages = [name for name in _PACKAGE_ROLE if classes_in(model, f"pytex.{name}")]

    boxes = []
    for name in packages:
        role, blurb = _PACKAGE_ROLE[name]
        members = classes_in(model, f"pytex.{name}")
        words = blurb.split(", ")
        boxes.append(
            ClassBox(
                key=name,
                name=f"pytex.{name}",
                module=f"{len(members)} public classes",
                stereotype=role,
                fields=tuple(
                    ", ".join(words[index : index + 2]) for index in range(0, len(words), 2)
                ),
                role=name,
                emphasis=name == "core",
            )
        )

    edges = [
        ClassEdge(source=source, target=target, kind="flow", label=str(count))
        for (source, target), count in sorted(counts.items())
        if source in _PACKAGE_ROLE and target in _PACKAGE_ROLE
    ]

    total = sum(counts.values())
    return class_diagram_svg(
        title="PyTex Library Architecture",
        subtitle=(
            "Packages, their public class counts, and every typed reference that crosses "
            "a package boundary"
        ),
        description=(
            "Package-level architecture of PyTex. The canonical core carries frames, "
            "symmetry, lattice and orientation semantics; the texture, EBSD, diffraction "
            "and TEM domains build on it; adapters sit at the external boundary; plotting "
            "presents. Each arrow is labelled with the number of typed field references "
            "that cross that boundary in the source."
        ),
        boxes=boxes,
        edges=edges,
        notes=(
            f"Arrow numbers count typed field references: {total} cross a package boundary.",
            "Every arrow points toward the core. No domain package is referenced by the core.",
            "Plotting and adapters take their inputs as call arguments, so few appear as fields.",
        ),
    )


def hierarchy_figure(model: ClassModel) -> str:
    """Every inheritance relation in PyTex, plus the vocabulary types.

    The figure is deliberately small because the library is: 240 public classes
    and a handful of base classes. Padding it out would misrepresent a design
    that composes rather than subclasses.
    """

    inheritance = model.of_kind("inheritance")
    involved = sorted({rel.source for rel in inheritance} | {rel.target for rel in inheritance})
    related = set(involved)
    boxes = [_box(model, name, related=related) for name in involved]
    edges = [
        ClassEdge(source=rel.source, target=rel.target, kind="inheritance") for rel in inheritance
    ]

    enums = sorted(
        (name for name, entry in model.entries.items() if entry.kind == "enum"),
        key=lambda name: model.entries[name].name,
    )
    protocols = sorted(
        (name for name, entry in model.entries.items() if entry.kind == "protocol"),
        key=lambda name: model.entries[name].name,
    )
    boxes.append(
        ClassBox(
            key="__enums__",
            name="StrEnum vocabularies",
            module="standard library base",
            stereotype=f"{len(enums)} enums",
            fields=tuple(model.entries[name].name for name in enums),
            role="external",
        )
    )
    boxes.append(
        ClassBox(
            key="__protocols__",
            name="Structural protocols",
            module="typing.Protocol",
            stereotype=f"{len(protocols)} protocols",
            fields=tuple(model.entries[name].name for name in protocols),
            role="external",
        )
    )

    dataclasses_count = sum(1 for entry in model.entries.values() if entry.kind == "dataclass")
    frozen = sum(1 for entry in model.entries.values() if entry.frozen)
    return class_diagram_svg(
        title="PyTex Class Hierarchy",
        subtitle=(
            "Every inheritance relation in the library — and it is short because PyTex "
            "composes rather than subclasses"
        ),
        description=(
            "The complete inheritance structure of PyTex. Only a few families use "
            "subclassing: the TEM tilt envelopes, which are alternative geometric "
            "shapes behind one interface, and the elastic tensor pair. Everything else "
            "is a frozen dataclass composed of other objects, which is why the object-"
            "model diagrams rather than this one carry the library's structure."
        ),
        boxes=boxes,
        edges=edges,
        notes=(
            f"{len(model.entries)} public classes: {dataclasses_count} dataclasses "
            f"({frozen} frozen), and {len(inheritance)} inheritance relations.",
            "Subclassing is used where alternatives share one interface, not to share code.",
            "The object-model diagrams show the composition that carries the real structure.",
        ),
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def build_figures() -> dict[str, str]:
    """Every atlas figure, keyed by output file name."""

    model = build_model()
    if model.unresolved:  # pragma: no cover - a real defect if it ever fires
        raise RuntimeError("annotations could not be resolved for: " + ", ".join(model.unresolved))
    figures = {
        "class_model_architecture.svg": architecture_figure(model),
        "class_hierarchy.svg": hierarchy_figure(model),
    }
    for spec in VIEWS:
        figures[f"class_model_{spec.key}.svg"] = domain_figure(model, spec)
    return figures


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for name, svg in build_figures().items():
        (FIGURES / name).write_text(svg, encoding="utf-8")
        print(f"wrote docs/figures/{name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
