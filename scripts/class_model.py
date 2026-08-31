"""Introspect PyTex and build its class and object model.

This module is the *source of truth extractor* behind the Class & Object Model
Atlas in the Sphinx site. It imports the public `pytex` packages and reads what
is actually declared — base classes, dataclass fields, resolved type hints — so
that every diagram in the atlas is a statement about the code as it exists, not
a hand-drawn approximation that decays the moment a field is renamed.

Why introspection rather than pyreverse or ``sphinx.ext.inheritance_diagram``:
both render through the Graphviz ``dot`` binary, which is a system dependency
the docs build does not otherwise carry, and neither reads dataclass field types
— the very relation that carries PyTex's structure. The renderer used instead is
`pytex.plotting.class_diagrams`, on the same canonical SVG stack as the
reference-frame and algorithm figures.

What the model records
----------------------
``ClassEntry``
    One public class: its module, its kind (dataclass, enum, protocol, plain
    class), its declared fields with shortened type names, and its docstring
    summary line.
``Relation``
    A typed edge. ``inheritance`` comes from ``__bases__``. ``composition`` and
    ``association`` come from dataclass field annotations: a required field is
    read as composition (the owner cannot exist without it), an optional or
    defaulted field as association. Multiplicity is read from the annotation —
    ``*`` for a sequence, ``0..1`` for an optional, ``1`` otherwise.

A deliberate finding, recorded here because it governs the atlas
----------------------------------------------------------------
PyTex is composition-first. Of its public classes, the overwhelming majority are
frozen dataclasses and there are only a handful of internal inheritance edges.
An inheritance-led atlas in the style of a C++ project would be almost empty and
would misdescribe the library. The atlas therefore leads with the object model,
and shows the class hierarchy in full precisely because it is small.

Usage::

    from class_model import build_model
    model = build_model()
    print(len(model.entries), len(model.relations))
"""

from __future__ import annotations

import dataclasses
import enum
import functools
import importlib
import inspect
import pkgutil
import sys
import types
import typing
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy.typing

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # pragma: no cover - import-path setup
    sys.path.insert(0, str(REPO_ROOT / "src"))

__all__ = [
    "ClassEntry",
    "ClassModel",
    "Relation",
    "build_model",
    "short_type_name",
]

#: Modules excluded from the model: private helpers, vendored data, the
#: plotting validation-case registry, and the application layer.
#:
#: ``app`` is excluded deliberately. The atlas documents the *scientific* object
#: model — phases, orientations, patterns, reports — and `pytex.app` contributes
#: only user-interface plumbing (parameter descriptors, result envelopes, wire
#: forms) that would triple the class count without describing any
#: crystallography. The application's own architecture is documented in
#: ``docs/architecture/application_platform.md``.
_EXCLUDED_MODULE_PARTS = ("_data", "themes", "fixtures", "app")


@dataclass(frozen=True)
class FieldEntry:
    """One declared dataclass field, with its type rendered for display."""

    name: str
    type_name: str
    required: bool
    """False when the field has a default or a default factory."""


@dataclass(frozen=True)
class ClassEntry:
    """One public class as declared in the source."""

    qualname: str
    """Fully qualified name, e.g. ``pytex.core.lattice.Phase``."""

    name: str
    module: str
    kind: str
    """One of ``dataclass``, ``enum``, ``protocol``, ``class``."""

    summary: str
    """First sentence of the docstring, or an empty string."""

    fields: tuple[FieldEntry, ...] = ()
    frozen: bool = False

    @property
    def package(self) -> str:
        """The subpackage the class belongs to, e.g. ``core`` or ``diffraction``."""

        parts = self.module.split(".")
        return parts[1] if len(parts) > 1 else ""


@dataclass(frozen=True)
class Relation:
    """A typed edge between two classes in the model."""

    source: str
    target: str
    kind: str
    """``inheritance``, ``composition`` or ``association``."""

    label: str = ""
    """The field name carrying the relation; empty for inheritance."""

    multiplicity: str = "1"
    """``1``, ``0..1`` or ``*``, read from the annotation."""


@dataclass(frozen=True)
class ClassModel:
    """Every public class and every relation between them."""

    entries: Mapping[str, ClassEntry]
    relations: tuple[Relation, ...]
    unresolved: tuple[str, ...] = ()
    """Classes whose annotations could not be resolved; empty in a healthy tree."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "entries", dict(self.entries))

    def of_kind(self, kind: str) -> tuple[Relation, ...]:
        """Every relation of one kind, in declaration order."""

        return tuple(rel for rel in self.relations if rel.kind == kind)

    def induced(self, qualnames: Iterable[str]) -> tuple[Relation, ...]:
        """Relations whose endpoints both lie in ``qualnames``.

        Self-edges are dropped: a class that references its own type (a tree node
        holding children, for instance) adds a loop that no layered layout can
        draw legibly, and the field list already states it.
        """

        selected = set(qualnames)
        seen: set[tuple[str, str, str, str]] = set()
        kept: list[Relation] = []
        for rel in self.relations:
            if rel.source not in selected or rel.target not in selected:
                continue
            if rel.source == rel.target:
                continue
            key = (rel.source, rel.target, rel.kind, rel.label)
            if key in seen:
                continue
            seen.add(key)
            kept.append(rel)
        return tuple(kept)

    def degree(self, qualname: str, within: Iterable[str] | None = None) -> int:
        """How many relations touch ``qualname`` (optionally within a subset)."""

        scope = set(within) if within is not None else None
        total = 0
        for rel in self.relations:
            if scope is not None and (rel.source not in scope or rel.target not in scope):
                continue
            if qualname in (rel.source, rel.target):
                total += 1
        return total


# ---------------------------------------------------------------------------
# Type-annotation reading
# ---------------------------------------------------------------------------


def _iter_annotation(annotation: object) -> Iterator[object]:
    """Every type mentioned anywhere inside an annotation."""

    stack = [annotation]
    while stack:
        current = stack.pop()
        yield current
        stack.extend(typing.get_args(current))


def _is_optional(annotation: object) -> bool:
    """True for ``X | None`` in either the ``Optional`` or PEP 604 spelling."""

    origin = typing.get_origin(annotation)
    if origin is typing.Union or origin is types.UnionType:
        return type(None) in typing.get_args(annotation)
    return False


def _is_many(annotation: object) -> bool:
    """True when the annotation holds a variable number of the referenced type.

    ``tuple[X, Y]`` is *not* many — it is a pair of distinct roles — while
    ``tuple[X, ...]``, ``Sequence[X]``, and ``dict[str, X]`` all are.
    """

    for current in _iter_annotation(annotation):
        origin = typing.get_origin(current)
        if origin is None:
            continue
        args = typing.get_args(current)
        if origin is tuple:
            if len(args) == 2 and args[1] is Ellipsis:
                return True
            continue
        if isinstance(origin, type) and issubclass(origin, (list, set, frozenset, dict, tuple)):
            return True
        name = getattr(origin, "_name", "") or getattr(origin, "__name__", "")
        if name in {"Sequence", "Iterable", "Mapping", "Collection", "List", "Tuple"}:
            return True
    return False


@functools.lru_cache(maxsize=1)
def _expandable_aliases() -> tuple[tuple[str, frozenset[object]], ...]:
    """Named unions the type system may flatten, and what they flatten into.

    NumPy's ``ArrayLike`` is a PEP 695 ``TypeAliasType`` on Python 3.12 and
    later and a plain ``Union`` alias on 3.11, so ``get_type_hints`` keeps the
    name on one interpreter and expands it into a dozen implementation types on
    the other. A card that reads ``ArrayLike`` here and
    ``_Buffer | _SupportsArray | ...`` on a 3.11 runner is not a canonical
    asset: the atlas SVGs are compared byte for byte, and this alone made them
    fail on half the CI matrix.

    The alias name is what the source wrote, so the alias name is what the card
    shows, whichever way the interpreter resolved it. Longest first, so a
    larger alias is recognised before any alias contained inside it.
    """

    aliases: list[tuple[str, frozenset[object]]] = []
    for name in ("ArrayLike", "DTypeLike"):
        alias = getattr(numpy.typing, name, None)
        if alias is None:  # pragma: no cover - a NumPy without the alias
            continue
        # `__value__` is the TypeAliasType's right-hand side; on 3.11 the alias
        # *is* its right-hand side already.
        members = typing.get_args(getattr(alias, "__value__", alias))
        if members:
            aliases.append((name, frozenset(members)))
    aliases.sort(key=lambda item: len(item[1]), reverse=True)
    return tuple(aliases)


def _collapse_expanded_aliases(members: list[object]) -> list[object | str]:
    """Rewrite flattened alias members back to the alias name. See above."""

    rendered: list[object | str] = list(members)
    for name, expansion in _expandable_aliases():
        present = [item for item in rendered if item in expansion]
        if len(present) != len(expansion):
            continue
        first = rendered.index(present[0])
        rendered = [item for item in rendered if item not in expansion]
        rendered.insert(first, name)
    return rendered


def short_type_name(annotation: object) -> str:
    """Render an annotation compactly enough to sit inside a diagram card.

    Module paths are stripped (``pytex.core.lattice.Phase`` becomes ``Phase``),
    ``X | None`` becomes ``X?``, homogeneous tuples become ``X[]``, and NumPy
    arrays collapse to ``ndarray`` — the element type of an array is carried by
    the docstring, not by a name a reader can act on.
    """

    if annotation is None or annotation is type(None):
        return "None"
    if annotation is Ellipsis:
        return "..."
    if isinstance(annotation, str):
        return annotation.rsplit(".", 1)[-1]

    origin = typing.get_origin(annotation)
    args = typing.get_args(annotation)

    if origin is typing.Union or origin is types.UnionType:
        inner = [arg for arg in args if arg is not type(None)]
        parts = _collapse_expanded_aliases(inner)
        rendered = " | ".join(
            part if isinstance(part, str) else short_type_name(part) for part in parts
        )
        return f"{rendered}?" if len(args) != len(inner) else rendered

    if origin is not None:
        base = getattr(origin, "_name", None) or getattr(origin, "__name__", str(origin))
        if base in {"ndarray", "NDArray"}:
            return "ndarray"
        if origin is tuple and len(args) == 2 and args[1] is Ellipsis:
            return f"{short_type_name(args[0])}[]"
        if base in {"list", "List", "Sequence", "frozenset", "set"} and args:
            return f"{short_type_name(args[0])}[]"
        if base in {"dict", "Dict", "Mapping"} and len(args) == 2:
            return f"{{{short_type_name(args[0])}: {short_type_name(args[1])}}}"
        if args:
            rendered = ", ".join(short_type_name(arg) for arg in args)
            return f"{base}[{rendered}]"
        return str(base)

    name = getattr(annotation, "__name__", None)
    if name:
        return "ndarray" if name in {"ndarray", "NDArray"} else name
    return str(annotation).rsplit(".", 1)[-1]


def _summary(obj: type) -> str:
    """First sentence of a docstring, normalized to one line."""

    doc = inspect.getdoc(obj) or ""
    first = doc.strip().split("\n\n", 1)[0].replace("\n", " ").strip()
    if not first:
        return ""
    sentence = first.split(". ")[0].strip()
    return sentence.rstrip(".")


def _kind(obj: type) -> str:
    if isinstance(obj, type) and issubclass(obj, enum.Enum):
        return "enum"
    if getattr(obj, "_is_protocol", False):
        return "protocol"
    if dataclasses.is_dataclass(obj):
        return "dataclass"
    return "class"


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------


def _public_modules(package: types.ModuleType) -> list[str]:
    names: list[str] = [package.__name__]
    for info in pkgutil.walk_packages(package.__path__, f"{package.__name__}."):
        parts = info.name.split(".")
        if any(part.startswith("_") for part in parts[1:]):
            continue
        if any(part in _EXCLUDED_MODULE_PARTS for part in parts):
            continue
        names.append(info.name)
    return sorted(names)


def build_model(package_name: str = "pytex") -> ClassModel:
    """Import ``package_name`` and read its public class and object model.

    Returns every public class declared in a public module — classes are
    attributed to the module that *defines* them, so a symbol re-exported by a
    package ``__init__`` appears once, at its definition site.
    """

    package = importlib.import_module(package_name)
    modules: dict[str, types.ModuleType] = {}
    for name in _public_modules(package):
        modules[name] = importlib.import_module(name)

    classes: dict[type, str] = {}
    entries: dict[str, ClassEntry] = {}
    unresolved: list[str] = []

    for module_name, module in modules.items():
        for attr, obj in sorted(vars(module).items()):
            if not inspect.isclass(obj) or attr.startswith("_"):
                continue
            if obj.__module__ != module_name:
                continue
            classes[obj] = f"{module_name}.{attr}"

    bare_index = _bare_name_index(classes)
    hints_by_class: dict[str, dict[str, object]] = {}

    for obj, qualname in classes.items():
        hints: dict[str, object] = {}
        if dataclasses.is_dataclass(obj):
            hints, ok = _resolve_hints(obj, bare_index)
            if not ok:
                unresolved.append(qualname)
        hints_by_class[qualname] = hints

        fields: list[FieldEntry] = []
        if dataclasses.is_dataclass(obj):
            for f in dataclasses.fields(obj):
                annotation = hints.get(f.name, f.type)
                required = (
                    f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING
                )
                fields.append(
                    FieldEntry(
                        name=f.name,
                        type_name=short_type_name(annotation),
                        required=required,
                    )
                )
        elif issubclass(obj, enum.Enum):
            fields = [
                FieldEntry(name=member.name.lower(), type_name="", required=True) for member in obj
            ]

        entries[qualname] = ClassEntry(
            qualname=qualname,
            name=obj.__name__,
            module=obj.__module__,
            kind=_kind(obj),
            summary=_summary(obj),
            fields=tuple(fields),
            frozen=bool(getattr(obj, "__dataclass_params__", None))
            and bool(obj.__dataclass_params__.frozen),  # type: ignore[attr-defined]
        )

    relations = _relations(classes, hints_by_class)
    return ClassModel(entries=entries, relations=relations, unresolved=tuple(sorted(unresolved)))


def _bare_name_index(classes: Mapping[type, str]) -> dict[str, type]:
    """Bare class name -> class, for names that are unambiguous library-wide.

    Several modules import the types they reference under ``TYPE_CHECKING`` to
    break an import cycle (``ebsd.models`` and ``texture``, for instance). Those
    names are absent from the module globals at runtime, so ``get_type_hints``
    cannot resolve the annotation that mentions them — and those annotations
    carry exactly the cross-domain relations the atlas exists to show. Supplying
    the library's own class names as a fallback namespace resolves them. Only
    names that are unique across the library are offered, and only where the
    module does not already define the name itself, so the fallback can never
    silently redirect a relation to a different class.
    """

    counts: dict[str, list[type]] = {}
    for obj in classes:
        counts.setdefault(obj.__name__, []).append(obj)
    return {name: found[0] for name, found in counts.items() if len(found) == 1}


def _resolve_hints(obj: type, bare_index: Mapping[str, type]) -> tuple[dict[str, object], bool]:
    """Resolved annotations for ``obj``, and whether resolution fully succeeded."""

    module_globals = vars(sys.modules[obj.__module__])
    fallback = {name: cls for name, cls in bare_index.items() if name not in module_globals}
    try:
        return typing.get_type_hints(obj, module_globals, fallback), True
    except Exception:  # pragma: no cover - defensive; empty in a healthy tree
        return {}, False


def _relations(
    classes: Mapping[type, str], hints_by_class: Mapping[str, dict[str, object]]
) -> tuple[Relation, ...]:
    """Every inheritance and field-typed relation between modelled classes."""

    relations: list[Relation] = []

    for obj, qualname in classes.items():
        for base in obj.__bases__:
            if base in classes:
                relations.append(
                    Relation(source=classes[base], target=qualname, kind="inheritance")
                )

    for obj, qualname in classes.items():
        if not dataclasses.is_dataclass(obj):
            continue
        hints = hints_by_class.get(qualname, {})
        for f in dataclasses.fields(obj):
            annotation = hints.get(f.name)
            if annotation is None:
                continue
            required = f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING
            optional = _is_optional(annotation)
            many = _is_many(annotation)
            multiplicity = "*" if many else ("0..1" if optional else "1")
            kind = "composition" if required and not optional else "association"
            seen: set[str] = set()
            for referenced in _iter_annotation(annotation):
                if not isinstance(referenced, type):
                    continue
                target = classes.get(referenced)
                if target is None or target in seen:
                    continue
                seen.add(target)
                relations.append(
                    Relation(
                        source=qualname,
                        target=target,
                        kind=kind,
                        label=f.name,
                        multiplicity=multiplicity,
                    )
                )

    relations.sort(key=lambda rel: (rel.kind, rel.source, rel.target, rel.label))
    return tuple(relations)


def package_flow(model: ClassModel) -> tuple[tuple[str, str, int], ...]:
    """Cross-package reference counts, source package -> target package.

    This is what the architecture overview draws: not an intention about how the
    layers *should* depend on one another, but a count of the typed references
    that actually cross each boundary.
    """

    counts: dict[tuple[str, str], int] = {}
    for rel in model.relations:
        source = model.entries[rel.source].package
        target = model.entries[rel.target].package
        if not source or not target or source == target:
            continue
        counts[(source, target)] = counts.get((source, target), 0) + 1
    return tuple(
        (source, target, count)
        for (source, target), count in sorted(counts.items(), key=lambda item: item[0])
    )


def most_connected(model: ClassModel, qualnames: Sequence[str], *, limit: int) -> tuple[str, ...]:
    """The ``limit`` most connected classes within ``qualnames``.

    Ties break on the qualified name so the selection is reproducible, which
    matters: the figures are committed assets checked byte-for-byte.
    """

    scope = list(qualnames)
    ranked = sorted(scope, key=lambda name: (-model.degree(name, within=scope), name))
    return tuple(ranked[:limit])


def classes_in(model: ClassModel, module_prefix: str) -> tuple[str, ...]:
    """Every modelled class whose module starts with ``module_prefix``."""

    return tuple(
        sorted(
            name
            for name, entry in model.entries.items()
            if entry.module == module_prefix or entry.module.startswith(f"{module_prefix}.")
        )
    )


_ = field  # re-exported dataclasses helper kept importable for view specs
