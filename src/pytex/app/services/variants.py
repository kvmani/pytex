"""Transformation variants: where they point, and how they differ from each other.

What it does
    Takes an orientation relationship and shows the *set* of child orientations
    one parent grain produces. Two views, because the two questions a variant
    set is asked are different questions.

    The **pole figure** answers "where do the child poles land?". Every variant
    contributes its whole symmetry family of a chosen child plane, projected
    into the parent frame, so the figure is directly comparable with a measured
    pole figure from the same parent grain. Variants can be coloured by packet,
    which is what makes the structure visible: under Kurdjumov-Sachs the 24
    variants fall into 4 packets of 6, one per member of the parent {111}
    family, and the packet is what a lath martensite micrograph shows as a
    block.

    The **intervariant misorientation spectrum** answers "what boundaries do
    these variants make with each other?". Two child grains descended from one
    parent can only meet at certain misorientations, and that discrete spectrum
    is what an EBSD misorientation histogram is compared against when deciding
    whether a region shares a parent.

When to use it
    Before interpreting a variant-selection argument, to know what no selection
    looks like. After indexing an EBSD map, to compare the measured boundary
    spectrum against the admissible one. When teaching why one parent grain
    gives 24 child orientations and only 4 apparent plate directions.

What every pole carries
    Its variant, its packet, its indices in the child crystal, its position on
    the projection and its polar and azimuthal angles — the same row the CSV
    export writes, so a hover and a file cannot disagree.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import DependencyMissingError, InvalidInputError
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import (
    _RELATIONSHIP_CONSTRUCTORS,
    _RELATIONSHIPS,
    direction_label,
    phase_parameter,
    plane_label,
    relationship_name,
)
from pytex.app.services.crystal import (
    _EULER_CONVENTIONS,
    _euler_convention,
    scene_payload,
)
from pytex.core.miller import canonicalize_sign

__all__: tuple[str, ...] = ()

_CITATION_MORITO = (
    "Morito, Tanaka, Konishi, Furuhara & Maki, Acta Mater. 51 (2003) 1789 "
    "(packet and block structure of lath martensite)."
)
_CITATION_BUNGE = "Bunge, Texture Analysis in Materials Science (1982), chapter 2."
_CITATION_RANDLE = (
    "Randle & Engler, Introduction to Texture Analysis, 2nd ed., chapter 2 (projections)."
)

#: Columns of the pole table, shared by the hover card and every export.
_POLE_COLUMNS: tuple[Column, ...] = (
    Column("variant", "Variant", numeric=True, help_text="Index in `generate_variants()` order."),
    Column(
        "packet",
        "Packet",
        numeric=True,
        help_text=(
            "Variants sharing the parent plane they carry into exact parallelism. Under "
            "Kurdjumov-Sachs and the {111} family this is the packet of lath martensite."
        ),
    ),
    Column("pole", "Child plane", help_text="The family member of the plotted child plane."),
    Column(
        "x",
        "x",
        numeric=True,
        digits=5,
        help_text="Projected coordinate in the unit disc, in the parent frame.",
    ),
    Column("y", "y", numeric=True, digits=5),
    Column(
        "polar_deg",
        "Polar angle",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle from the projection pole (the parent z axis). 0° is the centre.",
    ),
    Column(
        "azimuth_deg",
        "Azimuth",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle around the projection, measured from the parent x axis.",
    ),
)

#: Columns of the intervariant misorientation table.
_PAIR_COLUMNS: tuple[Column, ...] = (
    Column("variant_a", "Variant A", numeric=True),
    Column("variant_b", "Variant B", numeric=True),
    Column(
        "angle_deg",
        "Disorientation",
        units="°",
        numeric=True,
        digits=4,
        help_text="Minimal angle over the child point group, which is what EBSD reports.",
    ),
    Column(
        "axis",
        "Axis",
        help_text=(
            "The nearest low-index direction of the child crystal to the rotation axis of the "
            "minimal representative. These axes are not all rational, so the deviation column "
            "beside it says how far the label is from the exact axis."
        ),
    ),
    Column(
        "axis_deviation_deg",
        "Axis deviation",
        units="°",
        numeric=True,
        digits=3,
        help_text="How far the exact axis lies from the low-index label. 0 means it is exact.",
    ),
    Column(
        "axis_x",
        "Axis x",
        numeric=True,
        digits=5,
        help_text="The exact axis, in Cartesian child-crystal components, for re-plotting.",
    ),
    Column("axis_y", "Axis y", numeric=True, digits=5),
    Column("axis_z", "Axis z", numeric=True, digits=5),
    Column(
        "same_packet",
        "Same packet",
        help_text="Whether the two variants share a parent habit plane.",
    ),
)


#: Packet colours, as hex, matching `PACKET_COLORS` in `js/panels/variants.js`.
#:
#: The two lists are the same colours in two languages, which is a duplication
#: worth naming: the browser cannot read a Python constant and matplotlib cannot
#: read a CSS one. `test_app_variants.py` compares them, so a change to either
#: without the other is a test failure rather than a figure that quietly stops
#: matching the screen it was published from.
_PACKET_COLORS: tuple[str, ...] = (
    "#2c7fdd",
    "#e87917",
    "#259d71",
    "#ce3b98",
    "#8c56d2",
    "#b8a41e",
)

#: The parent's own poles: achromatic, because they are the reference.
_PARENT_COLOR = "#333333"


def _unit_from_angles(polar_deg: float, azimuth_deg: float) -> list[float]:
    """Rebuild a unit vector from the polar and azimuthal angles in a pole row.

    The figure is drawn from the table rather than from the vectors that
    produced it, so that the published figure and the exported numbers are the
    same numbers. The round trip is exact to floating point: the table stores
    the angles of the folded vector, and the projection depends on nothing else.
    """

    polar = math.radians(polar_deg)
    azimuth = math.radians(azimuth_deg)
    return [
        math.sin(polar) * math.cos(azimuth),
        math.sin(polar) * math.sin(azimuth),
        math.cos(polar),
    ]


#: Largest |u|, |v| or |w| considered when naming a disorientation axis.
#:
#: Eight rather than six because the martensite intervariant axes are quoted in
#: the literature up to ⟨5 5 7⟩ (Morito et al., Table 2), and a limit of six
#: cannot reach them: it labelled the 14.88° pair as ⟨2 5 0⟩ nearly four degrees
#: away, which is a label worse than none. Beyond eight an axis is irrational
#: for practical purposes, and the deviation column says so rather than the
#: label growing longer.
_AXIS_MAX_INDEX = 8


def _canonical_sign(triple: tuple[int, int, int]) -> tuple[int, int, int]:
    """Return the representative of an axis with its first nonzero index positive.

    A rotation axis is a line, not a ray: [1 -1 1] and [-1 1 -1] name the same
    axis, and the sign that comes out of a symmetry reduction is arbitrary.
    Without this the same physical axis appears under two labels in one table,
    and the reader is left to work out that they are the same thing.
    """

    for value in triple:
        if value != 0:
            return triple if value > 0 else (-triple[0], -triple[1], -triple[2])
    return triple


def _nearest_axis(axis_cartesian: np.ndarray, phase: Any, spec: PhaseSpec) -> tuple[str, float]:
    """Name the nearest low-index lattice direction to a Cartesian axis.

    A disorientation axis comes out of the symmetry reduction as three Cartesian
    components, and three decimals is not how anyone quotes one: the interesting
    fact about the 60° Kurdjumov-Sachs pairs is that they are about ⟨111⟩, and
    ``+0.577 +0.577 +0.577`` is that fact in a form that has to be decoded.

    The residual is returned with the label rather than hidden, because these
    axes are *not* all rational — an orientation relationship carrying a
    rational parent axis produces irrational child axes in general — and a
    label with no error attached would silently claim otherwise.

    Returns the label and the angular residual in degrees. An axis further than
    a few degrees from every candidate is reported as irrational by the caller.
    """

    from itertools import product
    from math import gcd

    # Columns of `matrix` are the Cartesian images of a, b, c, so `matrix @ uvw`
    # is the Cartesian vector of the direction [uvw].
    basis = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    target = np.asarray(axis_cartesian, dtype=float)
    norm = float(np.linalg.norm(target))
    if norm < 1e-12:
        return "—", 180.0
    target = target / norm

    # Candidate triples, primitive and sign-canonical: an axis and its reverse
    # describe the same rotation axis once the angle's sign is fixed, and a
    # non-primitive triple names the same direction as its reduced form.
    candidates: list[tuple[int, int, int]] = []
    span = range(-_AXIS_MAX_INDEX, _AXIS_MAX_INDEX + 1)
    for triple in product(span, span, span):
        if triple == (0, 0, 0):
            continue
        if gcd(gcd(abs(triple[0]), abs(triple[1])), abs(triple[2])) != 1:
            continue
        candidates.append(triple)

    images = np.asarray(candidates, dtype=float) @ basis.T
    units = images / np.linalg.norm(images, axis=1)[:, None]
    cosines = units @ target
    best = int(np.argmax(np.abs(cosines)))
    residual = float(np.degrees(np.arccos(min(abs(float(cosines[best])), 1.0))))
    return direction_label(_canonical_sign(candidates[best]), spec=spec), residual


def _relationship(name: str, parent: Any, child: Any) -> Any:
    """Build a named relationship, or explain why these phases cannot carry it."""

    from pytex.core.transformation import OrientationRelationship

    constructor = getattr(OrientationRelationship, _RELATIONSHIP_CONSTRUCTORS[name])
    try:
        return constructor(parent_phase=parent, child_phase=child)
    except (ValueError, TypeError) as error:
        raise InvalidInputError(
            f"The {relationship_name(name)} relationship does not apply to these phases: {error}",
            field="relationship",
            hint=(
                "The fcc-to-bcc relationships need a cubic parent and a cubic child; Burgers "
                "needs a cubic parent and a hexagonal child."
            ),
        ) from error


def _packet_labels(relationship: Any, parent_spec: PhaseSpec, indices: tuple[int, ...]) -> Any:
    """Group the variants by the parent plane each carries into parallelism.

    Raises an :class:`InvalidInputError` naming the packet-plane field rather
    than letting a library ``ValueError`` reach the user as an internal message:
    choosing a plane that is not the relationship's defining family is an easy
    and reasonable thing for a user to try.
    """

    from pytex.core.lattice import CrystalPlane, MillerIndex
    from pytex.core.transformation import variant_close_packed_groups

    parent_phase = relationship.parent_phase
    plane = CrystalPlane(
        MillerIndex(np.asarray(indices, dtype=int), phase=parent_phase),
        phase=parent_phase,
    )
    try:
        return variant_close_packed_groups(relationship, plane)
    except (ValueError, TypeError) as error:
        raise InvalidInputError(
            f"Variants cannot be grouped by {plane_label(indices, spec=parent_spec)} of "
            f"{parent_spec.name}: {error}",
            field="packet_plane",
            hint=(
                "Use the family the relationship is defined on: {111} for the fcc-to-bcc "
                "relationships, {110} for Burgers."
            ),
        ) from error


def _project(normals: np.ndarray, method: str) -> np.ndarray:
    """Project unit normals onto the unit disc, whichever projection is asked for.

    :func:`pytex.texture.projections.project_directions` returns each projection
    in its own natural radius: the stereographic net reaches 1 at the equator,
    the equal-area net reaches √2. Both are correct and the difference is
    invisible until something draws them in the same circle — at which point the
    equal-area figure spills over the rim by 41%.

    Dividing by the equatorial radius makes the rim mean "90° from the
    projection axis" in both, which is what a reader assumes and what makes the
    exported x and y comparable between the two. It changes no angle and no
    relative spacing: it is a choice of unit, fixed here so that the panel, the
    CSV and any figure built from the same rows cannot disagree about it.
    """

    from pytex.texture.projections import project_directions

    projected = np.asarray(project_directions(normals, method=method), dtype=float)
    equatorial_radius = math.sqrt(2.0) if method == "equal_area" else 1.0
    return projected / equatorial_radius


def _child_family(child_phase: Any, indices: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Return the symmetry family of a child plane, as indices and unit normals.

    Both are needed and they must correspond row for row: the normals are what
    is projected, the indices are what the hover card and the CSV name.
    """

    from pytex.core.miller import MillerPlaneSet

    single = MillerPlaneSet.from_hkl(np.asarray([indices], dtype=int), phase=child_phase)
    members, mask = single.symmetry_equivalent_indices()
    valid = np.asarray(members[0])[np.asarray(mask[0], dtype=bool)]
    if valid.size == 0:
        raise InvalidInputError(
            "That plane has no symmetry family; check the indices.", field="pole"
        )
    family = MillerPlaneSet.from_hkl(np.asarray(valid, dtype=int), phase=child_phase)
    normals = np.asarray(family.normals_cartesian(), dtype=float)
    return np.asarray(valid, dtype=int), normals


@REGISTRY.operation(
    "variants.pole_figure",
    title="Variant pole figure",
    summary="Where the child poles of every transformation variant land, in the parent frame.",
    help_text=(
        "Plots the chosen child plane family of every transformation variant, projected into "
        "the parent crystal frame. This is the figure that makes a variant set legible: one "
        "parent grain under Kurdjumov-Sachs produces 24 child orientations, and asking where "
        "each one puts its {100} is the question a measured pole figure of that grain answers "
        "too, so the two are directly comparable.\n\n"
        "**Colour by packet to see the structure.** The 24 Kurdjumov-Sachs variants fall into 4 "
        "groups of 6, one per member of the parent {111} family, because each variant carries "
        "exactly one {111} into exact parallelism with a child {110}. That group is the *packet* "
        "of lath martensite, and it is what a micrograph shows as a block of parallel laths. "
        "Burgers gives 6 groups of 2 on the parent {110} family instead.\n\n"
        "**Reading the projection.** The centre of the disc is the parent z axis and the rim is "
        "the equator; poles in the lower hemisphere are folded up, which is standard for a plane "
        "normal because a plane and its opposite normal are the same plane. Equal-area is the "
        "projection to use when density matters, because it does not concentrate area toward the "
        "rim; stereographic is the one that preserves angles, and is what a Wulff net measures.\n\n"
        "**What this does not claim.** Every variant is drawn with equal weight. A real "
        "microstructure shows variant *selection* — some variants far more often than others — "
        "and that is a property of the transformation conditions, not of the crystallography. "
        "This figure is the no-selection baseline a selection argument is measured against."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms, and the frame the figure is drawn in.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase — ferrite or martensite, or alpha for Burgers.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="pole",
            label="Child plane to plot",
            help_text=(
                "The child plane whose symmetry family is projected. (100) is the usual choice "
                "for a cubic product: three poles per variant, so the figure stays readable."
            ),
            default=(1, 0, 0),
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text=(
                "The parent family whose members the variants are grouped by. Use the family "
                "the relationship is built on: (111) for the fcc-to-bcc relationships, (110) "
                "for Burgers."
            ),
            default=(1, 1, 1),
        ),
        ChoiceParameter(
            name="projection",
            label="Projection",
            help_text="How the sphere is flattened onto the disc.",
            options=(
                (
                    "stereographic",
                    "Stereographic",
                    "Preserves angles; what a Wulff net measures.",
                ),
                (
                    "equal_area",
                    "Equal area (Schmidt)",
                    "Preserves area; use when pole density is the point.",
                ),
            ),
            default="stereographic",
        ),
        BooleanParameter(
            name="include_parent",
            label="Include the parent's own poles",
            help_text=(
                "Draws the same indices in the parent crystal, unrotated. It gives the figure a "
                "reference frame: the parent poles are where a pole would sit with no "
                "transformation at all."
            ),
            default=True,
        ),
    ),
    returns="One row per plotted pole; per-variant and per-packet summaries under `data`.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_RANDLE, _CITATION_BUNGE),
    tags=(
        "variant",
        "pole figure",
        "packet",
        "martensite",
        "orientation relationship",
        "OR",
        "stereographic",
        "projection",
    ),
)
def _variant_pole_figure(request: dict[str, Any]) -> dict[str, Any]:
    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    name = str(request["relationship"])
    relationship = _relationship(name, parent_phase, child_phase)
    variants = relationship.generate_variants()
    packets = _packet_labels(relationship, parent_spec, tuple(request["packet_plane"]))
    method = str(request["projection"])

    pole_indices, pole_normals = _child_family(child_phase, tuple(request["pole"]))

    # Every variant contributes the whole family at once. `parent_to_child`
    # carries a parent vector into the child frame, so its transpose is what
    # brings a child pole back into the parent frame the figure is drawn in.
    matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in variants], axis=0
    )
    in_parent = np.einsum("vji,pj->vpi", matrices, pole_normals, optimize=True)

    flat = in_parent.reshape(-1, 3)
    projected = _project(flat, method)
    # `project_directions` folds to the upper hemisphere before projecting, and
    # the polar angle must be read from the folded vector, not the original, or
    # a pole at 175° would be reported as such while being drawn at 5°.
    folded = np.where(flat[:, 2:3] < 0.0, -flat, flat)
    polar_deg = np.degrees(np.arccos(np.clip(folded[:, 2], -1.0, 1.0)))
    azimuth_deg = np.degrees(np.arctan2(folded[:, 1], folded[:, 0])) % 360.0

    rows: list[dict[str, Any]] = []
    position = 0
    for variant_position, variant in enumerate(variants):
        for member in pole_indices:
            rows.append(
                {
                    "variant": int(variant.variant_index),
                    "packet": int(packets[variant_position]) + 1,
                    "pole": plane_label(tuple(int(value) for value in member), spec=child_spec),
                    "x": float(projected[position, 0]),
                    "y": float(projected[position, 1]),
                    "polar_deg": float(polar_deg[position]),
                    "azimuth_deg": float(azimuth_deg[position]),
                }
            )
            position += 1

    parent_rows: list[dict[str, Any]] = []
    if bool(request["include_parent"]):
        parent_indices, parent_normals = _child_family(parent_phase, tuple(request["pole"]))
        parent_projected = _project(parent_normals, method)
        parent_folded = np.where(parent_normals[:, 2:3] < 0.0, -parent_normals, parent_normals)
        for member, point, vector in zip(
            parent_indices, parent_projected, parent_folded, strict=True
        ):
            parent_rows.append(
                {
                    "variant": 0,
                    "packet": 0,
                    "pole": plane_label(tuple(int(value) for value in member), spec=parent_spec),
                    "x": float(point[0]),
                    "y": float(point[1]),
                    "polar_deg": float(np.degrees(np.arccos(np.clip(vector[2], -1.0, 1.0)))),
                    "azimuth_deg": float(np.degrees(np.arctan2(vector[1], vector[0])) % 360.0),
                }
            )

    packet_sizes: dict[int, int] = {}
    for label in packets:
        packet_sizes[int(label) + 1] = packet_sizes.get(int(label) + 1, 0) + 1
    # Uneven packets mean the chosen parent family is not the one this
    # relationship is built on — grouping Pitsch, which is defined on {100}, by
    # {111} gives 4/3/4/1 rather than four equal groups. The grouping is still
    # computed, because it is a real property of a real question; but a figure
    # whose colours mean nothing is worse than one that says so.
    uneven = len(set(packet_sizes.values())) > 1
    variant_rows = [
        {
            "variant": int(variant.variant_index),
            "packet": int(packets[position]) + 1,
            "angle_deg": float(variant.parent_to_child_rotation.angle_deg),
        }
        for position, variant in enumerate(variants)
    ]

    display = relationship_name(name)
    pole_text = plane_label(tuple(int(v) for v in request["pole"]), spec=child_spec)
    packet_text = plane_label(tuple(int(v) for v in request["packet_plane"]), spec=parent_spec)
    packet_count = len(packet_sizes)
    sizes = sorted(set(packet_sizes.values()))
    size_text = f"{sizes[0]}" if len(sizes) == 1 else "/".join(str(value) for value in sizes)
    result = AppResult(
        title=f"{display} variant pole figure: {pole_text} of {child_spec.name}",
        summary=(
            f"{len(variants)} variants of {child_spec.name} in {parent_spec.name}, each "
            f"contributing its {len(pole_indices)}-member {pole_text} family, so "
            f"{len(rows)} poles on a "
            f"{'stereographic' if method == 'stereographic' else 'equal-area'} projection. "
            f"Grouping by {packet_text} of the parent gives {packet_count} packets of "
            f"{size_text} variants."
        ),
        table=ResultTable(
            columns=_POLE_COLUMNS,
            rows=tuple(rows + parent_rows),
            caption=(
                f"{pole_text} poles of every {display} variant, projected into the "
                f"{parent_spec.name} frame."
            ),
        ),
        data={
            "poles": rows,
            "parent_poles": parent_rows,
            "variants": variant_rows,
            "packet_sizes": packet_sizes,
            "packet_count": packet_count,
            "variant_count": len(variants),
            "family_size": len(pole_indices),
            "projection": method,
            "pole_label": pole_text,
            "packet_plane_label": packet_text,
            "columns": [column.to_json() for column in _POLE_COLUMNS],
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": name,
            "pole": [int(value) for value in request["pole"]],
            "packet_plane": [int(value) for value in request["packet_plane"]],
            "projection": method,
            "include_parent": bool(request["include_parent"]),
        },
        notes=(
            (
                f"The packets are uneven ({size_text} variants), which means "
                f"{packet_text} is not the family the {display} relationship is built on. The "
                "grouping is still the nearest-parallel parent plane of each variant, but it "
                "does not carry the packet meaning it does for the defining family."
            )
            if uneven
            else (
                f"Every packet holds {size_text} variants, as it must when the grouping plane "
                "is the family the relationship is defined on."
            ),
            "Every variant is drawn with equal weight. A real microstructure shows variant "
            "selection, which is a property of the transformation conditions rather than of the "
            "crystallography; this figure is the no-selection baseline.",
            "Poles in the lower hemisphere are folded onto the upper one, because a plane and "
            "its opposite normal describe the same plane.",
        ),
        citations=(_CITATION_MORITO, _CITATION_RANDLE),
    )
    return result.to_json()


@REGISTRY.operation(
    "variants.intervariant_misorientations",
    title="Intervariant misorientation spectrum",
    summary="Every boundary two child grains of one parent can make with each other.",
    help_text=(
        "Lists the disorientation of every pair of transformation variants. Two child grains "
        "that grew from the same parent cannot meet at an arbitrary misorientation: the set of "
        "possible boundaries is fixed by the relationship and by the two point groups, and it "
        "is discrete.\n\n"
        "That is what makes this table useful rather than merely descriptive. A measured "
        "misorientation histogram from an EBSD map of prior-austenite grains should show peaks "
        "at these angles and nowhere else; peaks away from them are boundaries between "
        "*different* parent grains. It is the same reasoning parent-grain reconstruction rests "
        "on.\n\n"
        "**Same-packet pairs are the low-angle ones.** Variants sharing a parent habit plane "
        "differ by small rotations about that plane normal, so they make the sub-block "
        "boundaries a micrograph shows as faint lines within a packet. Cross-packet pairs give "
        "the high-angle boundaries that delimit blocks.\n\n"
        "**Disorientation, not misorientation.** The angle reported is the minimum over the "
        "child point group, which is the convention EBSD software reports and the only one under "
        "which a histogram is comparable between tools."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase, whose symmetry reduces each misorientation.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text=(
                "Used only to label each pair as within or across a packet. (111) for the "
                "fcc-to-bcc relationships, (110) for Burgers."
            ),
            default=(1, 1, 1),
        ),
        BooleanParameter(
            name="merge_equal_angles",
            label="Report the spectrum, not every pair",
            help_text=(
                "On, pairs at the same angle are collapsed to one row carrying the count, which "
                "is the discrete spectrum to compare a histogram against. Off lists all pairs, "
                "which is what is wanted when a particular pair is the question."
            ),
            default=False,
            advanced=True,
        ),
    ),
    returns="One row per variant pair, or per distinct angle; the spectrum under `data`.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_BUNGE),
    tags=(
        "variant",
        "misorientation",
        "disorientation",
        "boundary",
        "packet",
        "block",
        "EBSD",
        "parent reconstruction",
    ),
)
def _intervariant_misorientations(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.transformation import intervariant_misorientations

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    name = str(request["relationship"])
    relationship = _relationship(name, parent_phase, child_phase)
    variants = relationship.generate_variants()
    if len(variants) < 2:
        raise InvalidInputError(
            f"The {relationship_name(name)} relationship gives only {len(variants)} variant, so "
            "there are no variant pairs to report.",
            field="relationship",
            hint="Bain gives 3 variants; Kurdjumov-Sachs gives 24.",
        )
    packets = _packet_labels(relationship, parent_spec, tuple(request["packet_plane"]))
    by_variant = {
        variant.variant_index: int(packets[position]) + 1
        for position, variant in enumerate(variants)
    }

    pairs = intervariant_misorientations(relationship, variants=variants)
    rows: list[dict[str, Any]] = []
    for pair in pairs:
        axis = np.asarray(pair.axis_child_frame, dtype=float)
        same = by_variant[pair.variant_a] == by_variant[pair.variant_b]
        label, deviation = _nearest_axis(axis, child_phase, child_spec)
        rows.append(
            {
                "variant_a": int(pair.variant_a),
                "variant_b": int(pair.variant_b),
                "angle_deg": float(pair.angle_deg),
                "axis": label,
                "axis_deviation_deg": deviation,
                "axis_x": float(axis[0]),
                "axis_y": float(axis[1]),
                "axis_z": float(axis[2]),
                "same_packet": "yes" if same else "no",
            }
        )
    rows.sort(key=lambda row: float(row["angle_deg"]))

    # The spectrum: distinct angles with their multiplicity. Rounding to three
    # decimals merges representatives that differ only at the floating-point
    # floor of chained matrix products, which is the whole point of asking for
    # a spectrum rather than for a list of pairs.
    spectrum: dict[float, dict[str, Any]] = {}
    for row in rows:
        key = round(float(row["angle_deg"]), 3)
        entry = spectrum.setdefault(
            key, {"angle_deg": key, "pairs": 0, "same_packet": 0, "cross_packet": 0}
        )
        entry["pairs"] += 1
        if row["same_packet"] == "yes":
            entry["same_packet"] += 1
        else:
            entry["cross_packet"] += 1
    spectrum_rows = sorted(spectrum.values(), key=lambda entry: float(entry["angle_deg"]))

    if bool(request["merge_equal_angles"]):
        columns: tuple[Column, ...] = (
            Column("angle_deg", "Disorientation", units="°", numeric=True, digits=3),
            Column("pairs", "Pairs", numeric=True),
            Column("same_packet", "Within a packet", numeric=True),
            Column("cross_packet", "Across packets", numeric=True),
        )
        table_rows: tuple[dict[str, Any], ...] = tuple(spectrum_rows)
        caption = f"Distinct disorientations among the {len(variants)} variants."
    else:
        columns = _PAIR_COLUMNS
        table_rows = tuple(rows)
        caption = f"Disorientation of every pair of the {len(variants)} variants."

    angles = np.asarray([row["angle_deg"] for row in rows], dtype=float)
    low_angle = int(np.count_nonzero(angles < 15.0))
    low_angle_within = sum(
        1 for row in rows if float(row["angle_deg"]) < 15.0 and row["same_packet"] == "yes"
    )
    # The within-packet spectrum is the interesting subset and is smaller than
    # the whole: for Kurdjumov-Sachs it is three angles out of ten. Saying "the
    # low-angle pairs are the within-packet ones" would be the tidy sentence and
    # the wrong one — 24 of the 48 pairs below 15 degrees cross a packet.
    within = sorted(
        {entry["angle_deg"] for entry in spectrum_rows if entry["same_packet"]},
    )
    within_pairs = sum(int(entry["same_packet"]) for entry in spectrum_rows)
    within_text = (
        "no pair shares a packet under this grouping"
        if not within
        else (
            f"the {within_pairs} pairs that share a packet take only "
            + ("one disorientation, " if len(within) == 1 else f"{len(within)} of them, ")
            + ", ".join(f"{value:g}°" for value in within)
        )
    )
    display = relationship_name(name)
    result = AppResult(
        title=f"{display}: misorientations among {len(variants)} variants",
        summary=(
            f"{len(rows)} variant pairs of {child_spec.name} from one {parent_spec.name} grain, "
            f"falling at {len(spectrum_rows)} distinct disorientations between "
            f"{angles.min():.2f}° and {angles.max():.2f}° — a discrete spectrum, not a spread. "
            f"Of those, {within_text}. "
            f"{low_angle} pairs sit below 15°, of which {low_angle_within} share a packet — "
            "sharing a packet and being low-angle are related but not the same thing."
        ),
        table=ResultTable(columns=columns, rows=table_rows, caption=caption),
        data={
            "pairs": rows,
            "spectrum": spectrum_rows,
            "variant_count": len(variants),
            "pair_count": len(rows),
            "distinct_angles": len(spectrum_rows),
            "min_angle_deg": float(angles.min()),
            "max_angle_deg": float(angles.max()),
            "low_angle_pairs": low_angle,
            "packets": by_variant,
            "columns": [column.to_json() for column in columns],
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": name,
            "packet_plane": [int(value) for value in request["packet_plane"]],
            "merge_equal_angles": bool(request["merge_equal_angles"]),
        },
        notes=(
            "Angles are disorientations: the minimum over the child point group, which is the "
            "convention EBSD software reports.",
            "This is the set of boundaries variants of one parent *can* make. Which of them a "
            "microstructure actually shows depends on variant selection and on which variants "
            "happen to be neighbours.",
        ),
        citations=(_CITATION_MORITO, _CITATION_BUNGE),
    )
    return result.to_json()


@REGISTRY.operation(
    "variants.render",
    title="Publication figure of the pole figure",
    summary="Render the variant pole figure through the publication renderer.",
    help_text=(
        "Draws the same poles as the interactive figure through "
        "`pytex.plotting.spherical` in the journal style, with the Wulff net, the packet "
        "colours, and the parent poles marked.\n\n"
        "The screen figure and this one are the same numbers seen through two renderers: the "
        "browser draws from the pole table, and so does this. Nothing is recomputed, so the "
        "figure cannot disagree with the CSV.\n\n"
        "**Choosing a format.** A pole figure is line art and a few dozen markers, so SVG is "
        "small, stays editable, and is the better artefact — the opposite of the crystal "
        "viewer, where a sphere mesh makes SVG the expensive choice. PNG at 600 dpi is offered "
        "for submission systems that will not take vector art."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms, and the frame the figure is drawn in.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="pole",
            label="Child plane to plot",
            help_text="The child plane whose symmetry family is projected.",
            default=(1, 0, 0),
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text="The parent family the variants are grouped and coloured by.",
            default=(1, 1, 1),
        ),
        ChoiceParameter(
            name="projection",
            label="Projection",
            help_text="How the sphere is flattened onto the disc.",
            options=(
                ("stereographic", "Stereographic", "Preserves angles."),
                ("equal_area", "Equal area (Schmidt)", "Preserves area."),
            ),
            default="stereographic",
        ),
        BooleanParameter(
            name="include_parent",
            label="Include the parent's own poles",
            help_text="Draws the same indices in the parent crystal, unrotated.",
            default=True,
        ),
        ChoiceParameter(
            name="format",
            label="Format",
            help_text="SVG for an editable vector figure; PNG for a raster one.",
            options=(
                ("svg", "SVG", "Vector, editable, and small for line art like this."),
                ("png", "PNG", "Raster, for submission systems that refuse vector art."),
            ),
            default="svg",
        ),
        IntegerParameter(
            name="dpi",
            label="Resolution",
            help_text="Pixels per inch, for PNG only.",
            default=600,
            minimum=72,
            maximum=1200,
            advanced=True,
        ),
    ),
    returns="The encoded image under `data.image`, with its format and encoding.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_RANDLE),
    tags=("variant", "pole figure", "figure", "export", "publication"),
)
def _variant_render(request: dict[str, Any]) -> dict[str, Any]:
    try:
        import matplotlib
    except ImportError as error:
        raise DependencyMissingError(
            "matplotlib", purpose="Rendering a publication figure", extra="plotting"
        ) from error

    matplotlib.use("Agg", force=False)
    import base64
    import io

    import matplotlib.pyplot as plt

    from pytex.plotting.spherical import plot_stereographic_vectors

    # The figure is drawn from the same computation the table reports, called
    # here rather than reimplemented, so the published figure and the exported
    # numbers cannot describe different poles.
    computed = _variant_pole_figure(
        {
            "phase": request["phase"],
            "child_phase": request["child_phase"],
            "relationship": request["relationship"],
            "pole": request["pole"],
            "packet_plane": request["packet_plane"],
            "projection": request["projection"],
            "include_parent": request["include_parent"],
        }
    )
    data = computed["data"]
    image_format = str(request["format"])
    dpi = int(request["dpi"])

    vectors: list[list[float]] = []
    colors: list[str] = []
    for pole in data["poles"]:
        vectors.append(_unit_from_angles(pole["polar_deg"], pole["azimuth_deg"]))
        colors.append(_PACKET_COLORS[(int(pole["packet"]) - 1) % len(_PACKET_COLORS)])
    for pole in data["parent_poles"]:
        vectors.append(_unit_from_angles(pole["polar_deg"], pole["azimuth_deg"]))
        colors.append(_PARENT_COLOR)

    axes = plot_stereographic_vectors(
        np.asarray(vectors, dtype=float),
        colors=colors,
        method=str(request["projection"]),
        render="pole",
        include_wulff_net=True,
        title=computed["title"],
    )
    figure = axes.get_figure()
    buffer = io.BytesIO()
    try:
        figure.savefig(buffer, format=image_format, dpi=dpi, bbox_inches="tight")
    finally:
        # A leaked figure is a defect here and a memory leak in a long-running
        # server, which is the same rule the crystal renderer follows.
        plt.close(figure)
    payload = buffer.getvalue()
    image = (
        payload.decode("utf-8")
        if image_format == "svg"
        else base64.b64encode(payload).decode("ascii")
    )

    result = AppResult(
        title=f"{computed['title']}: publication figure",
        summary=(
            f"{image_format.upper()} figure of {len(vectors)} poles "
            f"({len(data['poles'])} from {data['variant_count']} variants in "
            f"{data['packet_count']} packets, {len(data['parent_poles'])} from the parent), "
            f"rendered through the journal style"
            + (f" at {dpi} dpi" if image_format == "png" else "")
            + f" ({len(payload) / 1024:.0f} kB). The poles are the same rows the table exports."
        ),
        data={
            "image": image,
            "format": image_format,
            "encoding": "text" if image_format == "svg" else "base64",
            "bytes": len(payload),
            "dpi": dpi,
            "pole_count": len(vectors),
        },
        inputs=dict(computed["inputs"], format=image_format, dpi=dpi),
        citations=(_CITATION_MORITO, _CITATION_RANDLE),
    )
    return result.to_json()


#: Placement modes for the two crystals of a composite scene.
_PLACEMENTS: tuple[tuple[str, str, str], ...] = (
    (
        "interpenetrating",
        "Interpenetrating",
        "Both crystals share one origin. This is the placement in which the parallelism is "
        "visible, because the parallel plane and direction physically coincide.",
    ),
    (
        "side_by_side",
        "Side by side",
        "The child is translated clear of the parent along the world x axis. Easier to read "
        "each structure, at the cost of the coincidence being implied rather than seen.",
    ),
)


def _canonical_plane(values: Any) -> tuple[int, int, int]:
    """A plane's index triple under the repository's one sign rule.

    A plane has no sign, and which of ``(1 1 -1)`` and ``(-1 -1 1)`` a symmetry
    image comes back as is an artefact. ``canonicalize_sign`` is that rule, and
    it is deliberately *not* applied to directions, where the two spellings are
    opposite directions.
    """

    canonical = canonicalize_sign(_index_triple(values))[0]
    return (int(canonical[0]), int(canonical[1]), int(canonical[2]))


def _index_triple(values: Any) -> tuple[int, int, int]:
    """An index triple as three plain ints, rounded from float coordinates."""

    rounded = [round(float(value)) for value in np.asarray(values).reshape(-1)]
    if len(rounded) != 3:
        raise InvalidInputError("Expected a three-index triple.", field="relationship")
    return (rounded[0], rounded[1], rounded[2])


def _composite_relationship(request: dict[str, Any]) -> tuple[Any, Any, Any, Any, Any]:
    """Resolve both phases and the relationship, refusing phases with no atoms."""

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    for spec in (parent_spec, child_spec):
        if not spec.has_structure:
            raise InvalidInputError(
                f"{spec.name} carries no atomic basis, so there is nothing to draw.",
                field="phase" if spec is parent_spec else "child_phase",
                hint=(
                    "Choose a built-in phase, or add atomic sites to the phase description. "
                    "Lattice geometry alone supports the calculator but not the viewer."
                ),
            )
    relationship = _relationship(str(request["relationship"]), parent_phase, child_phase)
    return parent_spec, parent_phase, child_spec, child_phase, relationship


def _resolved_variants(relationship: Any) -> tuple[Any, ...]:
    variants = relationship.generate_variants()
    if not variants:  # pragma: no cover - a relationship always has at least one
        raise InvalidInputError(
            "This relationship generates no variants.", field="relationship"
        )
    return tuple(variants)


def _variant_at(relationship: Any, index: int) -> Any:
    """The one-based variant, or an error naming the field the user can fix."""

    variants = _resolved_variants(relationship)
    if not 1 <= index <= len(variants):
        raise InvalidInputError(
            f"This relationship has {len(variants)} variants, so variant {index} "
            "does not exist.",
            field="variant",
            hint=(
                f"Choose a variant between 1 and {len(variants)}. Kurdjumov-Sachs and "
                "Greninger-Troiano have 24; Nishiyama-Wassermann, Pitsch and Burgers have 12; "
                "Bain has 3."
            ),
        )
    return variants[index - 1]


def _child_translation(placement: str, parent_scene: Any, child_scene: Any) -> list[float]:
    """Where the child sits relative to the parent, for the chosen placement."""

    if placement != "side_by_side":
        return [0.0, 0.0, 0.0]
    parent_bounds = np.asarray(parent_scene.bounds(), dtype=float)
    child_bounds = np.asarray(child_scene.bounds(), dtype=float)
    span = float(parent_bounds[1][0] - parent_bounds[0][0])
    span += float(child_bounds[1][0] - child_bounds[0][0])
    return [0.65 * span, 0.0, 0.0]


def _primitive_payload(primitives: Any) -> dict[str, Any]:
    """The world-frame OR primitives as JSON the browser can draw directly.

    Everything is in Cartesian angstrom in the world frame — which is the parent
    crystal frame — so the browser applies the camera rotation and nothing else.
    """

    return {
        "arrows": [
            {
                "tail": [float(value) for value in arrow.tail],
                "head": [float(value) for value in arrow.head],
                "color": str(arrow.color),
                "label": arrow.label,
            }
            for arrow in primitives.arrows
        ],
        "patches": [
            {
                "vertices": [[float(value) for value in vertex] for vertex in patch.vertices],
                "normal": [float(value) for value in patch.normal],
                "color": str(patch.color),
                "alpha": float(patch.alpha),
                "label": patch.label,
            }
            for patch in primitives.patches
        ],
    }


def _parallelism_rows(
    variant: Any, *, parent_spec: PhaseSpec, child_spec: PhaseSpec
) -> list[dict[str, Any]]:
    """What *this* variant holds parallel, labelled with its own indices.

    Read from ``TransformationVariant.parallel_planes`` / ``.parallel_directions``
    rather than from the relationship: under variant k the parent-side objects
    are the symmetry images under that variant's operator, so quoting the
    relationship's nominal pair would label the figure with another variant's
    indices.
    """

    rows: list[dict[str, Any]] = []
    for parent_plane, child_plane in variant.parallel_planes:
        # A plane has no sign: (111) and (-1-1-1) name the same plane. Left
        # unsigned, the 24 Kurdjumov-Sachs variants appear to name eight parent
        # planes where they name four, and the packet column stops agreeing
        # with the plane column beside it.
        parent_indices = _canonical_plane(parent_plane.miller.indices)
        child_indices = _canonical_plane(child_plane.miller.indices)
        rows.append(
            {
                "kind": "plane",
                "parent": plane_label(parent_indices, spec=parent_spec),
                "child": plane_label(child_indices, spec=child_spec),
                "parent_indices": list(parent_indices),
                "child_indices": list(child_indices),
            }
        )
    for parent_direction, child_direction in variant.parallel_directions:
        parent_indices = _index_triple(parent_direction.coordinates)
        child_indices = _index_triple(child_direction.coordinates)
        rows.append(
            {
                "kind": "direction",
                "parent": direction_label(parent_indices, spec=parent_spec),
                "child": direction_label(child_indices, spec=child_spec),
                "parent_indices": list(parent_indices),
                "child_indices": list(child_indices),
            }
        )
    return rows


def _world_extent(points: np.ndarray) -> dict[str, Any]:
    """Centre and radius of a point cloud, for a camera that must frame both crystals."""

    lower = np.min(points, axis=0)
    upper = np.max(points, axis=0)
    centre = 0.5 * (lower + upper)
    radius = float(np.linalg.norm(upper - lower) / 2.0) or 1.0
    return {
        "centre": [float(value) for value in centre],
        "radius": radius,
        "bounds": [[float(v) for v in lower], [float(v) for v in upper]],
    }


_COMPOSITE_PARAMETERS = (
    phase_parameter(
        label="Parent phase",
        help_text="The phase that transforms. The world frame is its crystal frame.",
        builtin="austenite_fcc",
    ),
    phase_parameter(
        name="child_phase",
        label="Child phase",
        help_text="The product phase, placed by the variant's rotation.",
        builtin="fe_bcc",
    ),
    ChoiceParameter(
        name="relationship",
        label="Orientation relationship",
        help_text="Which relationship the two crystals are held in.",
        options=_RELATIONSHIPS,
        default="kurdjumov_sachs",
    ),
    IntegerParameter(
        name="repeats",
        label="Cells along each axis",
        help_text=(
            "How many unit cells of each crystal to build, in every direction. Two crystals "
            "are drawn at once, so this costs roughly twice what the single-crystal viewer does."
        ),
        default=1,
        minimum=1,
        maximum=3,
        group="Extent",
    ),
    ChoiceParameter(
        name="placement",
        label="Placement",
        help_text="Whether the two crystals share an origin or stand apart.",
        options=_PLACEMENTS,
        default="interpenetrating",
        group="Extent",
    ),
    BooleanParameter(
        name="show_parallel_planes",
        label="Draw the parallel planes",
        help_text=(
            "The variant's own parallel planes, as translucent patches in the world frame. "
            "Under the relationship they coincide exactly, which is the statement being made."
        ),
        default=True,
        group="Overlays",
    ),
    BooleanParameter(
        name="show_parallel_directions",
        label="Draw the parallel directions",
        help_text="The variant's own parallel directions, as arrows from the world origin.",
        default=True,
        group="Overlays",
    ),
    BooleanParameter(
        name="show_bonds",
        label="Draw bonds",
        help_text=(
            "Bonds are inferred from covalent radii plus a tolerance, not read from a file, "
            "so they aid reading rather than assert chemistry."
        ),
        default=True,
        advanced=True,
        group="Overlays",
    ),
    BooleanParameter(
        name="show_unit_cells",
        label="Outline every cell",
        help_text="Draw the edges of each repeated cell, not only the outer box.",
        default=True,
        advanced=True,
        group="Overlays",
    ),
)


@REGISTRY.operation(
    "variants.composite_scene",
    title="Both crystals of one variant",
    summary="Parent and product structures in one world frame, with the parallelism drawn on them.",
    help_text=(
        "Builds the two-crystal scene of a single transformation variant: the parent crystal in "
        "the world frame, the child placed by that variant's rotation, and the planes and "
        "directions the relationship holds parallel drawn across both. It is the picture behind "
        "the parallelism statement — instead of reading that (111) of austenite is parallel to "
        "(011) of ferrite, you see one plane through two lattices.\n\n"
        "**The variant matters, and so do its indices.** A relationship is realized by a family "
        "of variants — 24 for Kurdjumov-Sachs, 12 for Nishiyama-Wassermann and Burgers — and "
        "each one holds a *different* member of the parent family parallel. The overlay is "
        "labelled with the chosen variant's own indices, not with variant 1's, because drawing "
        "the nominal pair on variant 17 gives a picture that looks right and is wrong.\n\n"
        "**Everything crossing the wire is already placed.** Both structures come back in "
        "Cartesian angstrom in one world frame, so a viewer applies a camera rotation and "
        "nothing else; no crystallography happens in the browser, and one camera cannot drift "
        "between the two crystals."
    ),
    parameters=(
        *_COMPOSITE_PARAMETERS[:3],
        IntegerParameter(
            name="variant",
            label="Variant",
            help_text=(
                "Which variant to place, numbered as `generate_variants()` orders them. "
                "Variant 1 is the relationship exactly as it is written."
            ),
            default=1,
            minimum=1,
            maximum=24,
        ),
        *_COMPOSITE_PARAMETERS[3:],
    ),
    returns=(
        "Both placed scenes under `data.parent.scene` and `data.child.scene`, the world-frame "
        "overlays under `data.primitives`, and the variant's rotation under `data.variant`; "
        "the parallelisms it realizes as the table."
    ),
    panel="variants",
    citations=(_CITATION_MORITO,),
    tags=(
        "variant",
        "composite",
        "3D",
        "viewer",
        "orientation relationship",
        "OR",
        "parallel plane",
        "martensite",
    ),
)
def _variant_composite_scene(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.plotting.primitives import Transform3D
    from pytex.plotting.scene3d import WorldScene3D

    parent_spec, parent_phase, child_spec, child_phase, relationship = _composite_relationship(
        request
    )
    variants = _resolved_variants(relationship)
    variant = _variant_at(relationship, int(request["variant"]))
    repeats = int(request["repeats"])
    build_kwargs: dict[str, Any] = {
        "show_bonds": bool(request["show_bonds"]),
        "show_unit_cells": bool(request["show_unit_cells"]),
    }

    # The translation is measured from the two unplaced scenes, so a
    # side-by-side offset does not depend on which variant is showing and the
    # crystals do not jump as the user steps through the family.
    from pytex.plotting.crystal3d import build_crystal_scene

    parent_reference = build_crystal_scene(parent_phase, repeats=(repeats,) * 3, **build_kwargs)
    child_reference = build_crystal_scene(child_phase, repeats=(repeats,) * 3, **build_kwargs)
    translation = _child_translation(
        str(request["placement"]), parent_reference, child_reference
    )

    world = WorldScene3D.from_orientation_relationship(
        relationship,
        variant=variant,
        repeats=(repeats,) * 3,
        child_translation=translation,
        show_parallel_planes=bool(request["show_parallel_planes"]),
        show_parallel_directions=bool(request["show_parallel_directions"]),
        parent_build_kwargs=build_kwargs,
        child_build_kwargs=build_kwargs,
    )
    parent_placed, child_placed = world.placed_scenes()
    parent_payload = scene_payload(parent_placed, spec=parent_spec)
    child_payload = scene_payload(child_placed, spec=child_spec)
    bounds = np.asarray(world.bounds(), dtype=float)
    extent = _world_extent(bounds)

    rows = _parallelism_rows(variant, parent_spec=parent_spec, child_spec=child_spec)
    child_transform: Transform3D = world.crystals[1].transform
    misorientation = relationship.misorientation()

    result = AppResult(
        title=f"{parent_spec.name} and {child_spec.name}: variant {variant.variant_index}",
        summary=(
            f"Variant {variant.variant_index} of {len(variants)} under "
            f"{relationship_name(str(request['relationship']))}: "
            f"{len(parent_payload['atoms'])} parent atoms and "
            f"{len(child_payload['atoms'])} child atoms in one world frame, with "
            f"{len(rows)} parallelism(s) drawn across both. Both structures are already placed, "
            "in Cartesian angstrom in the parent crystal frame, so one camera drives them both. "
            "The overlay carries this variant's own indices, which are not variant 1's."
        ),
        table=ResultTable(
            columns=(
                Column("kind", "Kind", help_text="Whether the pair is a plane or a direction."),
                Column("parent", f"{parent_spec.name}"),
                Column("child", f"{child_spec.name}"),
            ),
            rows=tuple(
                {"kind": row["kind"], "parent": row["parent"], "child": row["child"]}
                for row in rows
            ),
            caption=(
                f"What variant {variant.variant_index} holds parallel. These are the variant's "
                "own symmetry images of the defining pair, not the relationship's nominal pair."
            ),
        ),
        data={
            "world": extent,
            "parent": {"label": parent_spec.name, "scene": parent_payload},
            "child": {"label": child_spec.name, "scene": child_payload},
            "primitives": _primitive_payload(world.primitives),
            "parallelisms": rows,
            "variant": {
                "index": int(variant.variant_index),
                "count": len(variants),
                "child_matrix": [
                    [float(value) for value in row] for row in child_transform.matrix
                ],
                "translation": [float(value) for value in translation],
                "parent_to_child_matrix": [
                    [float(value) for value in row]
                    for row in variant.parent_to_child_rotation.as_matrix()
                ],
            },
            "relationship": {
                "name": relationship_name(str(request["relationship"])),
                "angle_deg": float(misorientation.angle_deg),
            },
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": request["relationship"],
            "variant": int(variant.variant_index),
            "repeats": repeats,
            "placement": request["placement"],
            "show_parallel_planes": bool(request["show_parallel_planes"]),
            "show_parallel_directions": bool(request["show_parallel_directions"]),
            "show_bonds": bool(request["show_bonds"]),
            "show_unit_cells": bool(request["show_unit_cells"]),
        },
        citations=(_CITATION_MORITO,),
    )
    return result.to_json()


@REGISTRY.operation(
    "variants.contact_sheet",
    title="Every variant at once",
    summary="One placement matrix and one parallelism statement per variant, for an N-up grid.",
    help_text=(
        "The contact sheet behind the composite viewer: the same two structures, and the "
        "placement of the child for every variant of the relationship. Seeing 24 panels of one "
        "parent with 24 differently oriented children is what makes 'twenty-four variants' mean "
        "something.\n\n"
        "**Two scenes, N matrices.** The structures are sent once, each in its own crystal "
        "frame, together with one 3x3 placement matrix per variant. Sending 24 fully placed "
        "copies of both crystals would be tens of megabytes for information a matrix multiply "
        "reproduces exactly. Every matrix is computed here in Python; applying one is the same "
        "arithmetic the camera already does, so no crystallography moves into the browser.\n\n"
        "**Each panel carries its own parallelism.** Every variant holds a different member of "
        "the parent family parallel, and its row says which. Variants sharing that member form "
        "a packet — 4 packets of 6 under Kurdjumov-Sachs, 6 of 2 under Burgers — which is the "
        "structure a lath martensite micrograph shows as a block."
    ),
    parameters=(
        *_COMPOSITE_PARAMETERS,
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text=(
                "The parent family whose members the variants are grouped by. Use the family "
                "the relationship is built on: (111) for the fcc-to-bcc relationships, (110) "
                "for Burgers."
            ),
            default=(1, 1, 1),
            advanced=True,
        ),
    ),
    returns=(
        "The two structures under `data.parent.scene` and `data.child.scene`, each in its own "
        "crystal frame, and one entry per variant under `data.variants` carrying its placement "
        "matrix, its overlays and its packet; the same variants as the table."
    ),
    panel="variants",
    citations=(_CITATION_MORITO,),
    tags=(
        "variant",
        "contact sheet",
        "composite",
        "3D",
        "orientation relationship",
        "OR",
        "packet",
        "martensite",
    ),
)
def _variant_contact_sheet(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.plotting.crystal3d import build_crystal_scene
    from pytex.plotting.scene3d import WorldScene3D

    parent_spec, parent_phase, child_spec, child_phase, relationship = _composite_relationship(
        request
    )
    variants = _resolved_variants(relationship)
    repeats = int(request["repeats"])
    build_kwargs: dict[str, Any] = {
        "show_bonds": bool(request["show_bonds"]),
        "show_unit_cells": bool(request["show_unit_cells"]),
    }
    parent_scene = build_crystal_scene(parent_phase, repeats=(repeats,) * 3, **build_kwargs)
    child_scene = build_crystal_scene(child_phase, repeats=(repeats,) * 3, **build_kwargs)
    translation = _child_translation(str(request["placement"]), parent_scene, child_scene)
    packets = _packet_labels(
        relationship, parent_spec, tuple(int(v) for v in request["packet_plane"])
    )

    entries: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    cloud = [np.asarray(parent_scene.bounds(), dtype=float)]
    # All eight corners of the child box, not the two of its bounds array: a
    # rotated box is not bounded by the images of its min and max corners, so
    # taking only those would frame the sheet too tightly and clip variants.
    child_bounds = np.asarray(child_scene.bounds(), dtype=float)
    child_corners = np.array(
        [
            [child_bounds[i, 0], child_bounds[j, 1], child_bounds[k, 2]]
            for i in (0, 1)
            for j in (0, 1)
            for k in (0, 1)
        ],
        dtype=float,
    )
    for variant, packet in zip(variants, packets, strict=True):
        world = WorldScene3D.from_orientation_relationship(
            relationship,
            variant=variant,
            repeats=(repeats,) * 3,
            child_translation=translation,
            show_parallel_planes=bool(request["show_parallel_planes"]),
            show_parallel_directions=bool(request["show_parallel_directions"]),
            parent_build_kwargs=build_kwargs,
            child_build_kwargs=build_kwargs,
        )
        matrix = np.asarray(world.crystals[1].transform.matrix, dtype=float)
        parallelisms = _parallelism_rows(
            variant, parent_spec=parent_spec, child_spec=child_spec
        )
        entries.append(
            {
                "index": int(variant.variant_index),
                "packet": int(packet) + 1,
                "child_matrix": [[float(value) for value in row] for row in matrix],
                "translation": [float(value) for value in translation],
                "primitives": _primitive_payload(world.primitives),
                "parallelisms": parallelisms,
            }
        )
        plane_pairs = [row for row in parallelisms if row["kind"] == "plane"]
        direction_pairs = [row for row in parallelisms if row["kind"] == "direction"]
        rows.append(
            {
                "variant": int(variant.variant_index),
                "packet": int(packet) + 1,
                "planes": " ; ".join(f"{row['parent']} || {row['child']}" for row in plane_pairs),
                "directions": " ; ".join(
                    f"{row['parent']} || {row['child']}" for row in direction_pairs
                ),
            }
        )
        cloud.append((matrix @ child_corners.T).T + np.asarray(translation, dtype=float))

    extent = _world_extent(np.vstack(cloud))
    packet_count = len({entry["packet"] for entry in entries})
    result = AppResult(
        title=(
            f"{parent_spec.name} to {child_spec.name}: all "
            f"{len(variants)} {relationship_name(str(request['relationship']))} variants"
        ),
        summary=(
            f"{len(variants)} variants in {packet_count} packets. Both structures are sent once "
            f"({len(parent_scene.atoms)} parent atoms, {len(child_scene.atoms)} child atoms) in "
            "their own crystal frames, with one 3x3 placement matrix per variant; applying a "
            "matrix is the same arithmetic the camera does, so no crystallography moves into "
            "the viewer. Each variant carries its own parallel plane and direction, which are "
            "not the same indices from panel to panel."
        ),
        table=ResultTable(
            columns=(
                Column("variant", "Variant", numeric=True),
                Column(
                    "packet",
                    "Packet",
                    numeric=True,
                    help_text=(
                        "Variants sharing the parent plane they carry into exact parallelism."
                    ),
                ),
                Column("planes", "Parallel planes"),
                Column("directions", "Parallel directions"),
            ),
            rows=tuple(rows),
            caption=(
                "One row per variant. The indices are that variant's own symmetry images of the "
                "defining pair, which is why they differ down the column."
            ),
        ),
        data={
            "world": extent,
            "parent": {
                "label": parent_spec.name,
                "scene": scene_payload(parent_scene, spec=parent_spec),
            },
            "child": {
                "label": child_spec.name,
                "scene": scene_payload(child_scene, spec=child_spec),
            },
            "variants": entries,
            "variant_count": len(variants),
            "packet_count": packet_count,
            "frames": "own_crystal_frame",
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": request["relationship"],
            "repeats": repeats,
            "placement": request["placement"],
            "packet_plane": [int(v) for v in request["packet_plane"]],
            "show_parallel_planes": bool(request["show_parallel_planes"]),
            "show_parallel_directions": bool(request["show_parallel_directions"]),
            "show_bonds": bool(request["show_bonds"]),
            "show_unit_cells": bool(request["show_unit_cells"]),
        },
        citations=(_CITATION_MORITO,),
    )
    return result.to_json()


#: The four angles this panel puts on one screen, and what each one measures.
#:
#: They are four different quantities, and calling them all "deviation" would be
#: stating something false about three of them. The table's help text and the
#: summary both use these words.
_ANGLE_MEANINGS: dict[str, str] = {
    "residual": (
        "Scatter of the measurement: how far each measured pair sits from the single rotation "
        "fitted to all of them. Zero for one pair, by construction."
    ),
    "catalog": (
        "Distance from the fitted rotation to a named relationship, symmetry-reduced. This is "
        "what identifies the relationship."
    ),
    "rationalization": (
        "Price of writing the fitted rotation with integer indices: the symmetry-reduced angle "
        "between the integer statement and the fit. Compare it against the scatter."
    ),
    "clause": (
        "Per-clause rationalization residual: how far the exact image of the parent object sits "
        "from the integer child indices reported."
    ),
}


def _euler_angles(request: dict[str, Any], prefix: str) -> tuple[float, float, float]:
    values = [float(request[f"{prefix}_angle{index}"]) for index in (1, 2, 3)]
    if not all(math.isfinite(value) for value in values):
        raise InvalidInputError(
            "Euler angles must be finite numbers in degrees.",
            field=f"{prefix}_angle1",
        )
    return (values[0], values[1], values[2])


def _orientation_from_euler(
    angles: tuple[float, float, float],
    *,
    phase: Any,
    convention: str,
    frame: Any,
) -> Any:
    from pytex.core.orientation import Orientation

    return Orientation.from_euler(
        angles[0],
        angles[1],
        angles[2],
        specimen_frame=frame,
        symmetry=phase.symmetry,
        phase=phase,
        convention=convention,
        degrees=True,
    )


@REGISTRY.operation(
    "variants.or_from_grains",
    title="Relationship between two measured grains",
    summary="Two orientations in, a named relationship and its integer statement out.",
    help_text=(
        "The everyday EBSD question, answered end to end: two grains were indexed, one of each "
        "phase, and what is wanted is the orientation relationship between them — named if it is "
        "a known one, and written the way a paper writes it.\n\n"
        "Enter the two orientations as Euler angles in the convention your software exports "
        "(Bunge unless you know otherwise) and pick the two phases. What comes back is the "
        "rotation fitted to the pair, its distance from every relationship in the catalogue, a "
        "conclusive-or-not verdict, and the relationship restated in integers.\n\n"
        "**Four different angles appear here, and they are not interchangeable.** The *scatter* "
        "is how far the measured pairs sit from one fitted rotation — zero for a single pair, by "
        "construction, which is why one pair can never be contradicted by its own residual. The "
        "*catalogue distance* is how far the fit sits from a named relationship, and it is what "
        "identifies it. The *rationalization cost* is what writing the fit in integers costs. "
        "The *clause deviation* is how far one index pair sits from the exact image. Each is "
        "labelled where it appears.\n\n"
        "**The integer statement is an idealization, and its price is reported beside it.** "
        "Greninger-Troiano held to low indices comes out as the Kurdjumov-Sachs statement at a "
        "cost of 2.40 degrees, which is the separation between them; returned without that "
        "number it would read as a measurement of Kurdjumov-Sachs. Raise the index bound to buy "
        "a closer statement with untidier indices, and watch the cost fall."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase of the first grain — the one that transformed.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The phase of the second grain — the product.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="euler_convention",
            label="Euler convention",
            help_text=(
                "Which axis sequence the six angles below name. Both grains are read in the "
                "same convention, because they came from one indexing run."
            ),
            options=_EULER_CONVENTIONS,
            default="bunge",
        ),
        *(
            NumberParameter(
                name=f"parent_angle{index}",
                label=f"Parent {label}",
                help_text=f"{ordinal} Euler angle of the parent grain, in degrees.",
                units="deg",
                default=default,
                minimum=-360.0,
                maximum=720.0,
                group="Parent grain",
            )
            for index, label, ordinal, default in (
                (1, "phi1 / alpha", "First", 30.0),
                (2, "Phi / beta", "Second", 40.0),
                (3, "phi2 / gamma", "Third", 10.0),
            )
        ),
        *(
            NumberParameter(
                name=f"child_angle{index}",
                label=f"Child {label}",
                help_text=f"{ordinal} Euler angle of the child grain, in degrees.",
                units="deg",
                default=default,
                minimum=-360.0,
                maximum=720.0,
                group="Child grain",
            )
            for index, label, ordinal, default in (
                (1, "phi1 / alpha", "First", 45.2774),
                (2, "Phi / beta", "Second", 34.9979),
                (3, "phi2 / gamma", "Third", 316.2482),
            )
        ),
        NumberParameter(
            name="catalog_tolerance_deg",
            label="Naming tolerance",
            help_text=(
                "How close the fit must sit to a catalogued relationship before it is named. "
                "Three degrees is the working figure: above the orientation noise of a "
                "well-calibrated map, below the 5.26 degrees separating Kurdjumov-Sachs from "
                "Nishiyama-Wassermann."
            ),
            units="deg",
            default=3.0,
            minimum=0.1,
            maximum=15.0,
            advanced=True,
        ),
        IntegerParameter(
            name="max_index",
            label="Largest index in the statement",
            help_text=(
                "Bound on the integers the statement may use. Two gives the tidiest statement "
                "and the largest cost; raising it buys a closer one with untidier indices."
            ),
            default=3,
            minimum=1,
            maximum=6,
            advanced=True,
        ),
    ),
    returns=(
        "The catalogue ranking as the table; the fit, the verdict and the integer statement "
        "with its cost under `data`."
    ),
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_BUNGE),
    tags=(
        "orientation relationship",
        "OR",
        "EBSD",
        "measured",
        "grain",
        "Euler",
        "rationalization",
        "variant",
    ),
)
def _or_from_grains(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.frame_catalog import specimen_frame
    from pytex.core.orientation import OrientationSet
    from pytex.core.transformation import characterize_orientation_relationship

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    convention = _euler_convention(request["euler_convention"])
    frame = specimen_frame()
    parent_angles = _euler_angles(request, "parent")
    child_angles = _euler_angles(request, "child")
    parent_orientation = _orientation_from_euler(
        parent_angles, phase=parent_phase, convention=convention, frame=frame
    )
    child_orientation = _orientation_from_euler(
        child_angles, phase=child_phase, convention=convention, frame=frame
    )
    tolerance = float(request["catalog_tolerance_deg"])
    try:
        report = characterize_orientation_relationship(
            OrientationSet.from_orientations([parent_orientation]),
            OrientationSet.from_orientations([child_orientation]),
            catalog_tolerance_deg=tolerance,
        )
    except ValueError as error:
        raise InvalidInputError(
            f"These two grains cannot be characterized: {error}",
            field="phase",
            hint=(
                "The relationship is defined between two distinct phases; choose a parent and a "
                "child that are not the same phase."
            ),
        ) from error

    order = np.argsort(np.asarray(report.catalog_deviations_deg, dtype=float))
    rows = [
        {
            "relationship": relationship_name(report.catalog_names[int(index)]),
            "deviation_deg": float(report.catalog_deviations_deg[int(index)]),
            "within_tolerance": (
                "yes" if report.catalog_deviations_deg[int(index)] <= tolerance else "no"
            ),
        }
        for index in order
    ]

    max_index = int(request["max_index"])
    statement: dict[str, Any] | None = None
    statement_note: str | None = None
    try:
        rationalized = report.as_rational_relationship(
            max_index=max_index, tolerance_deg=tolerance
        )
    except ValueError as error:
        # Not an input error: a rotation with no low-index statement is a real
        # answer, and the panel says so rather than failing the whole run.
        statement_note = str(error)
    else:
        statement = {
            "text": rationalized.statement,
            "plane": {
                "parent": rationalized.plane_statement.parent_label,
                "child": rationalized.plane_statement.child_label,
                "deviation_deg": float(rationalized.plane_statement.deviation_deg),
            },
            "direction": {
                "parent": rationalized.direction_statement.parent_label,
                "child": rationalized.direction_statement.child_label,
                "deviation_deg": float(rationalized.direction_statement.deviation_deg),
            },
            "rationalization_cost_deg": float(rationalized.residual_rotation_deg),
            "zone_law_deviation_deg": float(rationalized.zone_law_deviation_deg),
            "max_index": rationalized.max_index,
            "describe": rationalized.describe(),
        }

    misorientation = report.relationship.misorientation()
    verdict = (
        f"{relationship_name(report.best_catalog_name)} within "
        f"{report.best_catalog_deviation_deg:.2f} deg"
        if report.is_conclusive and report.best_catalog_name is not None
        else "no conclusive match"
    )
    statement_text = (
        f" The integer statement is {statement['text']}, which costs "
        f"{statement['rationalization_cost_deg']:.2f} deg to write."
        if statement is not None
        else " No integer statement was found within the index bound; the rotation stands alone."
    )
    result = AppResult(
        title=f"{parent_spec.name} to {child_spec.name}: relationship between two grains",
        summary=(
            f"One measured pair. The relationship is a {misorientation.angle_deg:.2f} deg "
            "disorientation; the nearest catalogued relationship is "
            f"{relationship_name(report.catalog_names[int(order[0])])} at "
            f"{float(report.catalog_deviations_deg[int(order[0])]):.2f} deg, and the verdict is "
            f"{verdict}." + statement_text + " A single pair has no scatter to contradict it, so "
            "the residual column is identically zero and says nothing about the measurement."
        ),
        table=ResultTable(
            columns=(
                Column("relationship", "Relationship"),
                Column(
                    "deviation_deg",
                    "Catalogue distance",
                    units="°",
                    numeric=True,
                    digits=3,
                    help_text=_ANGLE_MEANINGS["catalog"],
                ),
                Column(
                    "within_tolerance",
                    "Within tolerance",
                    help_text="Whether this relationship is close enough to be named.",
                ),
            ),
            rows=tuple(rows),
            caption=(
                "Every catalogued relationship, ordered by how far the fitted rotation sits from "
                "it. A conclusive naming needs the winner to lead the runner-up by more than the "
                "scatter and its own misfit."
            ),
        ),
        data={
            "fit": {
                "angle_deg": float(misorientation.angle_deg),
                "axis": [float(value) for value in misorientation.rotation.axis],
                "matrix": [
                    [float(value) for value in row]
                    for row in report.relationship.parent_to_child_rotation.as_matrix()
                ],
                "mean_residual_deg": float(report.mean_residual_deg),
                "pair_count": int(report.pair_count),
                "converged": bool(report.converged),
            },
            "naming": {
                "best": report.best_catalog_name,
                "best_label": (
                    None
                    if report.best_catalog_name is None
                    else relationship_name(report.best_catalog_name)
                ),
                "best_deviation_deg": float(report.best_catalog_deviation_deg),
                "margin_deg": float(report.margin_deg),
                "is_conclusive": bool(report.is_conclusive),
                "tolerance_deg": tolerance,
            },
            "statement": statement,
            "statement_note": statement_note,
            "angle_meanings": dict(_ANGLE_MEANINGS),
            "euler_convention": convention,
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "euler_convention": convention,
            "parent_angle1": parent_angles[0],
            "parent_angle2": parent_angles[1],
            "parent_angle3": parent_angles[2],
            "child_angle1": child_angles[0],
            "child_angle2": child_angles[1],
            "child_angle3": child_angles[2],
            "catalog_tolerance_deg": tolerance,
            "max_index": max_index,
        },
        citations=(_CITATION_MORITO, _CITATION_BUNGE),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="variants.example.ks_packets",
            title="The 24 variants, and the 4 packets they fall into",
            panel="variants",
            summary="Kurdjumov-Sachs {100} poles of austenite-to-ferrite, coloured by packet.",
            teaches=(
                "Twenty-four variants put seventy-two poles on the figure, and they are not "
                "scattered: colouring by packet shows four groups of six, one per member of the "
                "parent {111} family. That grouping is what a lath martensite micrograph shows "
                "as a block, and it is why one parent grain gives twenty-four orientations but "
                "only four apparent plate directions."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.nw_twelve",
            title="Nishiyama-Wassermann: half as many variants",
            panel="variants",
            summary="The same figure under the 12-variant relationship, for comparison.",
            teaches=(
                "Run this straight after the Kurdjumov-Sachs figure. The same parent, the same "
                "poles plotted, half the variants — and the packets still number four, because "
                "packets are counted by the parent {111} family, which has four members "
                "whatever relationship sits on it. Nishiyama-Wassermann simply puts three "
                "variants in each rather than six."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "nishiyama_wassermann",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.bain_three",
            title="Bain: three variants and nothing to choose between",
            panel="variants",
            summary="The pure-strain path, whose three variants are the cube axes themselves.",
            teaches=(
                "The Bain path has no shear and only three variants, and its {100} poles land "
                "on the parent's own cube poles. That is why Bain is the reference the "
                "shear-carrying relationships are read against: it shows what the lattice "
                "correspondence alone requires, before any rotation is added to make a plane "
                "and a direction match."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "bain",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.burgers_poles",
            title="Burgers: six packets of two in zirconium",
            panel="variants",
            summary="The bcc-to-hcp path, with the basal pole of every variant plotted.",
            teaches=(
                "Burgers is grouped on the parent {110} family, which has six members, so the "
                "twelve variants fall into six packets of two rather than four of six. Plotting "
                "the child (0001) shows it directly: each packet contributes one basal pole "
                "position, because both of its variants share the parent plane the basal plane "
                "lies on."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "zr_bcc_beta"},
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
                "pole": [0, 0, 1],
                "packet_plane": [1, 1, 0],
            },
        ),
        ExampleScenario(
            id="variants.example.or_from_grains_ks",
            title="Two grains in, Kurdjumov-Sachs out",
            panel="variants",
            summary="An austenite and a ferrite orientation, and the relationship between them.",
            teaches=(
                "This is the everyday EBSD question run end to end. The two Euler triples are an "
                "exact Kurdjumov-Sachs pair, so pressing the button without touching anything "
                "recovers a known answer rather than producing an unverifiable number: "
                "Kurdjumov-Sachs at zero, Greninger-Troiano 2.40 degrees behind it and "
                "Nishiyama-Wassermann 5.26 degrees behind that. Those spacings are the "
                "literature's, and seeing the whole ladder is what makes the naming a judgement "
                "rather than an assertion."
            ),
            operation="variants.or_from_grains",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "euler_convention": "bunge",
                "parent_angle1": 30.0,
                "parent_angle2": 40.0,
                "parent_angle3": 10.0,
                "child_angle1": 45.2774,
                "child_angle2": 34.9979,
                "child_angle3": 316.2482,
                "max_index": 3,
            },
        ),
        ExampleScenario(
            id="variants.example.or_from_grains_gt_cost",
            title="What the tidy statement costs",
            panel="variants",
            summary=(
                "A Greninger-Troiano pair written with low indices, and the price of doing so."
            ),
            teaches=(
                "Run this after the Kurdjumov-Sachs example. The pair is exact "
                "Greninger-Troiano and the panel names it correctly at zero — but the integer "
                "statement, held to index two, comes out as the *Kurdjumov-Sachs* one, because "
                "Greninger-Troiano has no low-index direction pair. The cost of writing it that "
                "way is 2.40 degrees, exactly the separation between the two relationships, and "
                "it is reported beside the statement. Raise the index bound and watch the cost "
                "fall as the indices get untidier: that trade is the user's to make, so it is "
                "shown rather than decided."
            ),
            operation="variants.or_from_grains",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "euler_convention": "bunge",
                "parent_angle1": 30.0,
                "parent_angle2": 40.0,
                "parent_angle3": 10.0,
                "child_angle1": 126.5536,
                "child_angle2": 73.6774,
                "child_angle3": 352.9432,
                "max_index": 2,
            },
        ),
        ExampleScenario(
            id="variants.example.composite_variant_one",
            title="One plane through two lattices",
            panel="variants",
            summary="Austenite and ferrite in the Kurdjumov-Sachs relationship, variant 1.",
            teaches=(
                "The parallelism statement stops being a line of notation here. The parent and "
                "the product are drawn in one frame, and the plane the relationship holds "
                "parallel is a single translucent sheet cutting through both lattices with the "
                "shared direction as an arrow lying in it. Rotate it until you are looking down "
                "that arrow: the two structures are edge-on to the same plane, which is what "
                "(111) austenite parallel to (011) ferrite actually means."
            ),
            operation="variants.composite_scene",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "variant": 1,
                "repeats": 1,
                "placement": "interpenetrating",
            },
        ),
        ExampleScenario(
            id="variants.example.composite_variant_seventeen",
            title="The same relationship, a different variant",
            panel="variants",
            summary="Variant 17 of the same Kurdjumov-Sachs pair, for comparison with variant 1.",
            teaches=(
                "Run this straight after variant 1. The relationship has not changed and the "
                "two crystals are the same, but the product sits somewhere else entirely -- and "
                "the overlay is labelled with a different parent plane. That is the point: a "
                "variant is not a redrawing of the same picture, it is a different member of "
                "the parent family being carried into parallelism, and the label has to move "
                "with it. Drawing variant 1's indices here would give a figure that looks right "
                "and is wrong."
            ),
            operation="variants.composite_scene",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "variant": 17,
                "repeats": 1,
                "placement": "interpenetrating",
            },
        ),
        ExampleScenario(
            id="variants.example.contact_sheet_ks",
            title="All 24 variants, and the 4 packets in the table",
            panel="variants",
            summary="The Kurdjumov-Sachs contact sheet: one placement per variant.",
            teaches=(
                "Twenty-four rows, four values in the parallel-plane column, six rows each. "
                "That is the packet structure of lath martensite as a table rather than as a "
                "claim, and it comes out of the variants themselves: each carries exactly one "
                "member of the parent {111} family into parallelism, and there are four "
                "members. Compare the same run under Burgers, where the parent {110} family "
                "has six members and the twelve variants fall into six packets of two."
            ),
            operation="variants.contact_sheet",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "repeats": 1,
                "placement": "interpenetrating",
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.ks_spectrum",
            title="The boundaries 24 variants can make with each other",
            panel="variants",
            summary="The Kurdjumov-Sachs intervariant disorientation spectrum.",
            teaches=(
                "The 276 variant pairs do not spread over the angle range: they collapse onto a "
                "handful of discrete disorientations. That discreteness is the test — a "
                "measured misorientation histogram from prior-austenite grains should show "
                "peaks here and nowhere else, and a peak somewhere else is a boundary between "
                "two different parent grains."
            ),
            operation="variants.intervariant_misorientations",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "packet_plane": [1, 1, 1],
                "merge_equal_angles": True,
            },
        ),
        ExampleScenario(
            id="variants.example.ks_pairs",
            title="Which pairs are the low-angle ones",
            panel="variants",
            summary="Every Kurdjumov-Sachs variant pair, sorted by disorientation.",
            teaches=(
                "Sorted by angle, the top of the table is almost entirely same-packet pairs: "
                "variants sharing a parent {111} differ by small rotations about it, and those "
                "are the sub-block boundaries within a packet. Scroll down and the same-packet "
                "column turns to 'no' — those are the block boundaries a micrograph shows as "
                "sharp lines."
            ),
            operation="variants.intervariant_misorientations",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "packet_plane": [1, 1, 1],
                "merge_equal_angles": False,
            },
        ),
    )
)
