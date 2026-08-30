"""Composite 3D world scenes: several placed crystals plus loose primitives.

`build_crystal_scene` / `plot_crystal_structure_3d` model **one** crystal in its
own frame. This module lifts that to a shared **world frame** so a figure can
hold multiple crystals, each placed by a `Transform3D`, together with the
renderer-independent primitives from `pytex.plotting.primitives`:

- `PlacedCrystal` — a `CrystalScene` positioned by a `Transform3D` (from a
  `Rotation`/`Orientation`).
- `WorldScene3D` — an immutable bag of placed crystals and a `PrimitiveScene3D`;
  `add_crystal` / `add_primitives` grow it functionally, and
  `from_orientation_relationship` assembles the canonical two-crystal figure.
- `render_world_scene_3d` — draws every placed crystal into **one** globally
  depth-sorted `Poly3DCollection` (so atoms and bonds of different crystals
  occlude each other correctly) and overlays the primitives.

The orientation-relationship constructor places the child crystal by the
inverse of the parent-to-child rotation in the parent frame, which makes the
OR's parallel planes and directions coincide in world coordinates — the
geometric statement of the relationship, shown directly. Passing ``variant=k``
uses that variant's rotation *and* that variant's own parallel plane and
direction (its symmetry images, not the nominal pair), and
`WorldScene3D.variant_scenes` plus `render_variant_contact_sheet` put the whole
variant family on one sheet.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pytex.core.lattice import Phase
from pytex.core.miller import canonicalize_sign
from pytex.core.notation import format_miller_indices
from pytex.plotting.crystal3d import (
    CrystalScene,
    _accumulate_crystal_mesh,
    _apply_depth_cue,
    _draw_crystal_frame,
    _draw_crystal_planes_and_directions,
    _normalize_light_direction,
    _view_angles_from_direction,
    _view_vector_from_angles,
    build_crystal_scene,
)
from pytex.plotting.primitives import (
    Arrow3D,
    PlanePatch3D,
    PrimitiveScene3D,
    Transform3D,
    _draw_primitive_scene,
    crystal_plane_patch,
    plane_normal_endpoints,
    scene_span,
    segment_in_cell,
    segment_in_polygon,
)
from pytex.plotting.styles import resolve_style


def _to_hex(color: Any) -> str:
    """``matplotlib.colors.to_hex``, imported on demand.

    Matplotlib is an *optional* dependency behind the ``pytex[plotting]`` extra,
    and the repository forbids import-time coupling to optional stacks. A
    module-level ``from matplotlib.colors import to_hex`` made matplotlib
    mandatory for ``import pytex``; importing inside the call restores the
    declared contract at the cost of one cached dict lookup.
    """

    try:
        from matplotlib.colors import to_hex
    except ImportError as exc:  # pragma: no cover - exercised only without matplotlib
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return str(to_hex(color))


def _to_rgb(color: Any) -> tuple[float, float, float]:
    """``matplotlib.colors.to_rgb``, imported on demand. See :func:`_to_hex`."""

    try:
        from matplotlib.colors import to_rgb
    except ImportError as exc:  # pragma: no cover - exercised only without matplotlib
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    red, green, blue = to_rgb(color)
    return (float(red), float(green), float(blue))


# Distinct parent/child accent colors for orientation-relationship figures,
# drawn from the categorical palette so the two crystals stay separable.
_PARENT_ACCENT = "#2563eb"
_CHILD_ACCENT = "#dc2626"


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:  # pragma: no cover - environment-dependent branch
        raise ImportError(
            "PyTex plotting requires matplotlib. Install the 'pytex[plotting]' extra."
        ) from exc
    return plt, Poly3DCollection


@dataclass(frozen=True, slots=True)
class PlacedCrystal:
    """A `CrystalScene` positioned in the world frame by a `Transform3D`."""

    scene: CrystalScene
    transform: Transform3D = field(default_factory=Transform3D.identity)
    label: str | None = None

    def placed_scene(self) -> CrystalScene:
        """Return the scene with `transform` baked into every coordinate."""

        if np.allclose(self.transform.matrix, np.eye(3)) and np.allclose(
            self.transform.translation, 0.0
        ):
            return self.scene
        return self.scene.transformed(self.transform)


@dataclass(frozen=True, slots=True)
class WorldScene3D:
    """An immutable composite of placed crystals and world-frame primitives."""

    crystals: tuple[PlacedCrystal, ...] = ()
    primitives: PrimitiveScene3D = field(default_factory=PrimitiveScene3D)

    def add_crystal(
        self,
        scene_or_phase: CrystalScene | Phase,
        *,
        transform: Transform3D | None = None,
        label: str | None = None,
        **build_kwargs: Any,
    ) -> WorldScene3D:
        """Return a new world scene with one more placed crystal.

        ``scene_or_phase`` may be a prebuilt `CrystalScene` or a `Phase` (built
        via `build_crystal_scene` with ``build_kwargs``). ``transform`` places it
        in the world frame (default: identity).
        """

        scene = (
            scene_or_phase
            if isinstance(scene_or_phase, CrystalScene)
            else build_crystal_scene(scene_or_phase, **build_kwargs)
        )
        placed = PlacedCrystal(
            scene=scene,
            transform=transform or Transform3D.identity(),
            label=label,
        )
        return WorldScene3D(crystals=(*self.crystals, placed), primitives=self.primitives)

    def add_primitives(self, primitives: PrimitiveScene3D) -> WorldScene3D:
        """Return a new world scene with ``primitives`` merged in."""

        return WorldScene3D(
            crystals=self.crystals,
            primitives=self.primitives.merge(primitives),
        )

    def placed_scenes(self) -> tuple[CrystalScene, ...]:
        """The child scenes with their placement transforms already applied.

        Each sub-scene is defined in its own local frame; this returns them in
        world coordinates, which is what a renderer consumes.
        """

        return tuple(placed.placed_scene() for placed in self.crystals)

    def bounds(self) -> np.ndarray:
        """Axis-aligned ``(2, 3)`` [min; max] bounds over crystals and primitives."""

        mins: list[np.ndarray] = []
        maxs: list[np.ndarray] = []
        for scene in self.placed_scenes():
            crystal_bounds = scene.bounds()
            mins.append(crystal_bounds[0])
            maxs.append(crystal_bounds[1])
        if not self.primitives.is_empty():
            primitive_bounds = self.primitives.bounds()
            mins.append(primitive_bounds[0])
            maxs.append(primitive_bounds[1])
        if not mins:
            return np.zeros((2, 3), dtype=np.float64)
        return np.vstack([np.min(np.vstack(mins), axis=0), np.max(np.vstack(maxs), axis=0)])

    @classmethod
    def from_orientation_relationship(
        cls,
        relationship: Any,
        *,
        variant: int | Any | None = None,
        repeats: tuple[int, int, int] = (1, 1, 1),
        parent_transform: Transform3D | None = None,
        child_translation: Any = (0.0, 0.0, 0.0),
        show_parallel_directions: bool = True,
        show_parallel_planes: bool = True,
        show_plane_normals: bool = False,
        parent_build_kwargs: dict[str, Any] | None = None,
        child_build_kwargs: dict[str, Any] | None = None,
    ) -> WorldScene3D:
        """Assemble the two-crystal figure of an `OrientationRelationship`.

        The parent crystal is placed at ``parent_transform`` (default identity,
        i.e. the world frame *is* the parent crystal frame); the child is placed
        by the inverse of the parent-to-child rotation plus ``child_translation``,
        so the OR's parallel planes and directions coincide in world coordinates.
        With the defaults the parallel directions are drawn as arrows and the
        parallel planes as translucent patches, which visually verify the
        alignment.

        Parameters
        ----------
        variant : int or TransformationVariant, optional
            Which transformation variant to draw. ``None`` (the default) draws
            the relationship as stated. An ``int`` is a **one-based** index into
            ``relationship.generate_variants()``. The variant's own rotation
            places the child, and — this is the part that is easy to get wrong —
            the parallel planes and directions drawn are that variant's, taken
            from ``TransformationVariant.parallel_planes`` /
            ``.parallel_directions``, not the nominal pair the relationship was
            defined by. Drawing the nominal pair on variant 17 yields a figure
            that looks right and is wrong.

        Use it as the minimal-code entry point for OR schematics; pass
        ``child_translation`` to separate the two crystals side by side, or
        ``parent/child_build_kwargs`` to control repeats, bonds, or overlays.

        See Also
        --------
        variant_scenes : one scene per variant, for a contact sheet.
        render_variant_contact_sheet : draw such a tuple as a grid.
        """

        resolved = resolve_transformation_variant(relationship, variant)
        source = relationship if resolved is None else resolved
        parent_transform = parent_transform or Transform3D.identity()
        child_rotation = source.parent_to_child_rotation.inverse().as_matrix()
        child_transform = parent_transform.compose(
            Transform3D.from_matrix(child_rotation, translation=child_translation)
        )
        parent_kwargs = {"repeats": repeats, "show_unit_cells": True, **(parent_build_kwargs or {})}
        child_kwargs = {"repeats": repeats, "show_unit_cells": True, **(child_build_kwargs or {})}
        world = (
            cls()
            .add_crystal(
                relationship.parent_phase,
                transform=parent_transform,
                label=relationship.parent_phase.name,
                **parent_kwargs,
            )
            .add_crystal(
                relationship.child_phase,
                transform=child_transform,
                label=relationship.child_phase.name,
                **child_kwargs,
            )
        )
        primitives = _orientation_relationship_primitives(
            parallel_directions=source.parallel_directions if show_parallel_directions else (),
            parallel_planes=source.parallel_planes if show_parallel_planes else (),
            parent_transform=parent_transform,
            child_transform=child_transform,
            parent_phase=relationship.parent_phase,
            child_phase=relationship.child_phase,
            repeats=repeats,
            length=_relationship_reference_length(relationship, repeats),
            show_plane_normals=show_plane_normals,
        )
        if not primitives.is_empty():
            world = world.add_primitives(primitives)
        return world

    @classmethod
    def variant_scenes(
        cls,
        relationship: Any,
        *,
        variants: Sequence[Any] | None = None,
        **scene_kwargs: Any,
    ) -> tuple[WorldScene3D, ...]:
        """One composite scene per transformation variant, in variant order.

        Purpose
        -------
        The figure that shows what "24 variants" actually means: the same parent
        crystal with 24 differently oriented children, each carrying its *own*
        parallel plane and direction.

        Parameters
        ----------
        relationship : OrientationRelationship
        variants : sequence of TransformationVariant, optional
            Defaults to ``relationship.generate_variants()``.
        **scene_kwargs
            Forwarded to :meth:`from_orientation_relationship`; ``variant`` is
            supplied per scene and must not appear here.

        Returns
        -------
        tuple of WorldScene3D
            Parallel to ``variants``, so ``scenes[i]`` is the scene of
            ``variants[i]``. Feed it to :func:`render_variant_contact_sheet`.
        """

        if "variant" in scene_kwargs:
            raise ValueError("variant_scenes supplies 'variant' itself; pass 'variants' instead.")
        resolved = (
            relationship.generate_variants() if variants is None else tuple(variants)
        )
        return tuple(
            cls.from_orientation_relationship(relationship, variant=item, **scene_kwargs)
            for item in resolved
        )

    def render(self, **kwargs: Any) -> Any:
        """Render this world scene; forwards to `render_world_scene_3d`."""

        return render_world_scene_3d(self, **kwargs)


def resolve_transformation_variant(relationship: Any, variant: int | Any | None) -> Any | None:
    """Coerce a ``variant`` argument to a `TransformationVariant` (or ``None``).

    ``None`` passes through, meaning "the relationship as stated". An ``int`` is
    a **one-based** index into ``relationship.generate_variants()``, matching the
    ``variant_index`` the variants carry; anything else is returned unchanged and
    is expected to be a `TransformationVariant`.
    """

    if variant is None or not isinstance(variant, int):
        return variant
    variants = relationship.generate_variants()
    if not 1 <= variant <= len(variants):
        raise ValueError(
            f"variant must be a one-based index in 1..{len(variants)} for "
            f"'{relationship.name}'; got {variant}."
        )
    return variants[variant - 1]


def _relationship_reference_length(relationship: Any, repeats: tuple[int, int, int]) -> float:
    basis = relationship.parent_phase.lattice.direct_basis().matrix
    edge = float(np.max(np.linalg.norm(basis, axis=0)))
    return edge * float(max(repeats)) * 1.05


def _integer_indices(values: Any) -> list[int] | None:
    """``values`` rounded to integers, or ``None`` if they are not integral."""

    array = np.asarray(values, dtype=np.float64)
    rounded = np.rint(array)
    if not np.allclose(array, rounded, atol=1e-8):
        return None
    return [int(value) for value in rounded]


def _parallelism_label(parent_indices: Any, child_indices: Any, *, family: str) -> str:
    """``(111) ∥ (011)``-style label for a parallel pair, or a generic fallback."""

    parent = _integer_indices(parent_indices)
    child = _integer_indices(child_indices)
    kind = "plane" if family == "plane" else "direction"
    if parent is None or child is None:
        return f"∥ {kind}"
    if family == "plane":
        # A plane has no sign, and the sign a symmetry image comes back with is
        # an artefact. Spelling it one way here and another in a table beside it
        # leaves the reader to work out that they are the same plane.
        parent = [int(value) for value in canonicalize_sign(parent)[0]]
        child = [int(value) for value in canonicalize_sign(child)[0]]
    return (
        f"{format_miller_indices(parent, family=family, style='plain')}"
        f" ∥ {format_miller_indices(child, family=family, style='plain')}"
    )


def _orientation_relationship_primitives(
    *,
    parallel_directions: Any,
    parallel_planes: Any,
    parent_transform: Transform3D,
    child_transform: Transform3D | None = None,
    parent_phase: Any = None,
    child_phase: Any = None,
    repeats: tuple[int, int, int] = (1, 1, 1),
    length: float,
    show_plane_normals: bool = False,
) -> PrimitiveScene3D:
    """The parallel objects, drawn on **both** crystals and clipped to each cell.

    The pairs are passed in rather than read off the relationship, because under
    a transformation variant the objects that are actually parallel are that
    variant's symmetry images, not the relationship's nominal pair.

    Three properties, and each was previously false:

    **Every overlay is clipped to the cell it belongs to.** A plane is a
    statement about a lattice, so it is drawn where that lattice is — entering
    the box through one edge and leaving through another. The patches used to
    be squares centred on the world origin and sized by a scene-scale
    heuristic, which straddled the origin corner and hung outside the crystal.

    **Both crystals carry the overlay, each in its own cell.** The parent's
    ``(110)`` clipped to a cubic box and the child's ``(0001)`` clipped to a
    hexagonal one are the same physical plane with different outlines, and that
    is the truthful picture: the claim is that the planes are parallel, not that
    the cells are alike. A single patch translated onto the second crystal
    showed the parent's outline in the child's cell.

    **A direction lies in its plane, so its arrow is a chord of the patch.**
    Clipped to the polygon, the arrow is inside the cell by construction and is
    visibly *in* the plane it belongs to — which is the whole claim a
    plane-and-direction pair makes. It used to start at the origin and run on a
    scene-scale length, so it lay in nothing and left the crystal.
    """

    arrows: list[Arrow3D] = []
    patches: list[PlanePatch3D] = []

    def side(
        plane_index: int, is_parent: bool
    ) -> tuple[Any, Transform3D, Any] | None:
        """The plane, placement and phase of one side of the pair."""

        pairs = list(parallel_planes)
        if plane_index >= len(pairs):
            return None
        parent_plane, child_plane = pairs[plane_index]
        if is_parent:
            return parent_plane, parent_transform, parent_phase
        if child_transform is None:
            return None
        return child_plane, child_transform, child_phase

    # The polygons first: the arrows are chords of them, so they have to exist
    # before a direction can be placed in one.
    polygons: list[tuple[np.ndarray, Transform3D]] = []
    for index, (parent_plane, child_plane) in enumerate(parallel_planes):
        label = _parallelism_label(
            parent_plane.miller.indices, child_plane.miller.indices, family="plane"
        )
        for is_parent, plane in ((True, parent_plane), (False, child_plane)):
            resolved = side(index, is_parent)
            if resolved is None:
                continue
            _, side_transform, phase = resolved
            # Built in the crystal's own frame and placed afterwards, so the
            # cell the overlay is clipped against is the cell it belongs to
            # rather than an axis-aligned box in the world.
            local = crystal_plane_patch(
                plane,
                cell_repeats=repeats,
                color="#7c3aed",
                alpha=0.16,
                # The pair is named once, on the parent, so the two copies of
                # one statement do not print it twice over one figure.
                label=label if is_parent else "",
            )
            patch = local.transformed(side_transform)
            patches.append(patch)
            polygons.append((patch.vertices, side_transform))
            if show_plane_normals and phase is not None:
                spacing = float(plane.d_spacing_angstrom)
                tail, head = plane_normal_endpoints(
                    phase, local.vertices, local.normal, length=spacing, repeats=repeats
                )
                arrows.append(
                    Arrow3D(
                        tail=side_transform.apply_points(tail),
                        head=side_transform.apply_points(head),
                        color="#0f766e",
                        role="normal",
                        # The arrow is one interplanar spacing long, so it
                        # carries a measurement rather than a scene-scale
                        # guess: the gap it spans is the gap between this
                        # plane and the next of its family.
                        label=_normal_label(plane, spacing) if is_parent else "",
                    )
                )

    for parent_direction, child_direction in parallel_directions:
        label = _parallelism_label(
            parent_direction.coordinates, child_direction.coordinates, family="direction"
        )
        for is_parent, direction, phase, placement in (
            (True, parent_direction, parent_phase, parent_transform),
            (False, child_direction, child_phase, child_transform),
        ):
            if placement is None:
                continue
            transform: Transform3D = placement
            world_direction = transform.apply_vector(direction.unit_vector)
            world_direction = world_direction / np.linalg.norm(world_direction)
            tail, head = _direction_endpoints(
                world_direction,
                polygons=[points for points, owner in polygons if owner is transform],
                phase=phase,
                transform=transform,
                repeats=repeats,
                length=length,
            )
            arrows.append(
                Arrow3D(
                    tail=tail,
                    head=head,
                    color="#be123c",
                    role="direction",
                    label=label if is_parent else "",
                )
            )
    return PrimitiveScene3D(arrows=tuple(arrows), patches=tuple(patches))


def _normal_label(plane: Any, spacing: float) -> str:
    """``n (110) - d = 2.552 A``: which normal, and how far apart the planes are."""

    indices = _integer_indices(plane.miller.indices)
    if indices is None:
        return f"normal - d = {spacing:.3f} A"
    text = format_miller_indices(
        tuple(indices), family="plane", style="plain", scope="specific"
    )
    return f"n {text} - d = {spacing:.3f} A"


def _direction_endpoints(
    world_direction: np.ndarray,
    *,
    polygons: Sequence[np.ndarray],
    phase: Any,
    transform: Transform3D,
    repeats: tuple[int, int, int],
    length: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Where a direction arrow starts and ends: in its plane, or across the cell.

    A direction that lies in one of the drawn planes is drawn as a chord of
    that plane, because that is the relationship being asserted. One that lies
    in none of them — the planes switched off, or a pair whose direction is not
    in its plane — is clipped to the cell instead, so the rule that nothing
    leaves the crystal still holds without claiming a parallelism it has not
    got.
    """

    for polygon in polygons:
        chord = segment_in_polygon(polygon, world_direction)
        if chord is not None:
            return chord
    if phase is not None:
        # The crystal's own box, in the crystal frame, then placed.
        local_direction = np.linalg.solve(transform.matrix, world_direction)
        tail, head = segment_in_cell(phase, local_direction, repeats=repeats)
        return transform.apply_points(tail), transform.apply_points(head)
    origin = transform.translation
    return origin, origin + length * world_direction


def render_world_scene_3d(
    world: WorldScene3D,
    *,
    ax: Any | None = None,
    theme: str = "journal",
    style_path: str | None = None,
    style_overrides: dict[str, Any] | None = None,
    elev_deg: float = 22.0,
    azim_deg: float = 34.0,
    projection: str = "persp",
    view_direction: Any | None = None,
    show_legend: bool = False,
    title: str | None = None,
) -> Any:
    """Render a `WorldScene3D` to a single VESTA-class composite 3D figure.

    Every placed crystal's atom, bond, and polyhedron faces accumulate into one
    depth-sorted `Poly3DCollection`, so crystals occlude one another correctly;
    lattice frames, plane patches, direction arrows, and the loose primitives are
    drawn on top. Pass ``ax`` to compose into an existing 3D axes, ``view_direction``
    to look along a world vector, and ``show_legend`` for a merged species key.
    """

    plt, poly3d_collection = _require_matplotlib()
    style = resolve_style(theme=theme, style_path=style_path, overrides=style_overrides)
    common = style["common"]
    crystal_style = style["crystal"]
    background = crystal_style["background"]
    if ax is None:
        fig = plt.figure(
            figsize=tuple(common["figure"]["figsize"]),
            dpi=int(common["figure"]["dpi"]),
            facecolor=background,
        )
        axes = fig.add_subplot(111, projection="3d", proj_type=projection)
    else:
        axes = ax
        fig = axes.figure
    axes.set_facecolor(background)
    light_direction = _normalize_light_direction(crystal_style["light_direction"])
    # resolve the viewing direction FIRST so depth cueing can fade along it
    if view_direction is not None:
        elev_deg, azim_deg = _view_angles_from_direction(
            np.asarray(view_direction, dtype=np.float64)
        )
    bounds = world.bounds()
    span = scene_span(bounds)
    placed_scenes = world.placed_scenes()
    mesh_faces: list[np.ndarray] = []
    mesh_colors: list[np.ndarray] = []
    camera_direction = _view_vector_from_angles(elev_deg, azim_deg)
    for scene in placed_scenes:
        _draw_crystal_frame(axes, scene, crystal_style)
        faces, colors = _accumulate_crystal_mesh(
            axes,
            scene,
            crystal_style,
            light_direction=light_direction,
            view_direction=camera_direction,
        )
        mesh_faces.extend(faces)
        mesh_colors.extend(colors)
    if mesh_faces:
        all_faces = np.concatenate(mesh_faces, axis=0)
        all_colors = _apply_depth_cue(
            all_faces,
            np.concatenate(mesh_colors, axis=0),
            elev_deg=elev_deg,
            azim_deg=azim_deg,
            strength=float(crystal_style.get("depth_cue_strength", 0.0)),
            background=str(background),
        )
        mesh = poly3d_collection(
            all_faces,
            facecolors=all_colors,
            edgecolors="none",
            linewidths=0.0,
        )
        mesh.set_zsort("average")
        axes.add_collection3d(mesh)
    for scene in placed_scenes:
        _draw_crystal_planes_and_directions(axes, scene, crystal_style, scene_span=span)
    _draw_primitive_scene(axes, world.primitives, scene_span=span)
    center = 0.5 * (bounds[0] + bounds[1])
    radius = 0.55 * span
    axes.set_xlim(center[0] - radius, center[0] + radius)
    axes.set_ylim(center[1] - radius, center[1] + radius)
    axes.set_zlim(center[2] - radius, center[2] + radius)
    axes.set_box_aspect((1.0, 1.0, 1.0))
    axes.view_init(elev=elev_deg, azim=azim_deg)
    if show_legend:
        _draw_species_legend(axes, placed_scenes)
    if bool(crystal_style.get("hide_grid", True)):
        axes.grid(False)
    pane_rgba = (*_to_rgb(background), float(crystal_style["pane_alpha"]))
    axes.xaxis.set_pane_color(pane_rgba)
    axes.yaxis.set_pane_color(pane_rgba)
    axes.zaxis.set_pane_color(pane_rgba)
    if not bool(crystal_style.get("show_axes", False)):
        axes.set_axis_off()
    if title is not None:
        axes.set_title(title)
    fig.tight_layout()
    return fig


def _draw_species_legend(axes: Any, placed_scenes: tuple[CrystalScene, ...]) -> None:
    from matplotlib.lines import Line2D

    species_colors: dict[str, str] = {}
    for scene in placed_scenes:
        for atom in scene.atoms:
            species_colors.setdefault(atom.species, atom.color)
    if not species_colors:
        return
    handles = [
        Line2D(
            [0.0],
            [0.0],
            marker="o",
            linestyle="none",
            markersize=9,
            markerfacecolor=color,
            markeredgecolor="#334155",
            markeredgewidth=0.4,
            label=species,
        )
        for species, color in sorted(species_colors.items())
    ]
    axes.legend(handles=handles, loc="upper right", framealpha=0.85)


def render_variant_contact_sheet(
    scenes: Sequence[WorldScene3D],
    *,
    variants: Sequence[Any] | None = None,
    titles: Sequence[str] | None = None,
    columns: int = 4,
    panel_size_inches: float = 2.4,
    suptitle: str | None = None,
    **render_kwargs: Any,
) -> Any:
    """Draw one composite scene per panel of a grid — the variant contact sheet.

    Purpose
    -------
    Twenty-four separate figures do not show variant selection; one sheet does.
    Each panel is a full `render_world_scene_3d` of the corresponding scene, so
    every panel carries its own parallel plane and direction.

    Parameters
    ----------
    scenes : sequence of WorldScene3D
        Typically :meth:`WorldScene3D.variant_scenes` output.
    variants : sequence of TransformationVariant, optional
        Used only to title the panels ``V1``, ``V2``, ... by
        ``variant_index``; must be the same length as ``scenes``.
    titles : sequence of str, optional
        Explicit panel titles, overriding ``variants``.
    columns : int
        Panels per row; the row count follows. Must be strictly positive.
    panel_size_inches : float
        Edge length of one square panel; the figure size follows from the grid.
    suptitle : str, optional
    **render_kwargs
        Forwarded to `render_world_scene_3d` for every panel (``ax`` and
        ``title`` are supplied here and must not appear).

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If ``scenes`` is empty, ``columns`` is not positive, a label sequence
        has the wrong length, or ``ax``/``title`` is passed through.
    """

    plt, _ = _require_matplotlib()
    panels = tuple(scenes)
    if not panels:
        raise ValueError("render_variant_contact_sheet requires at least one scene.")
    if columns <= 0:
        raise ValueError("columns must be strictly positive.")
    for reserved in ("ax", "title"):
        if reserved in render_kwargs:
            raise ValueError(f"render_variant_contact_sheet controls '{reserved}' itself.")
    if titles is not None:
        panel_titles: tuple[str | None, ...] = tuple(titles)
    elif variants is not None:
        panel_titles = tuple(f"V{variant.variant_index}" for variant in variants)
    else:
        panel_titles = tuple(f"V{index + 1}" for index in range(len(panels)))
    if len(panel_titles) != len(panels):
        raise ValueError("titles/variants must have the same length as scenes.")
    rows = -(-len(panels) // columns)
    figure = plt.figure(figsize=(columns * panel_size_inches, rows * panel_size_inches))
    for index, (scene, panel_title) in enumerate(zip(panels, panel_titles, strict=True)):
        axes = figure.add_subplot(rows, columns, index + 1, projection="3d")
        render_world_scene_3d(scene, ax=axes, title=panel_title, **render_kwargs)
    if suptitle is not None:
        figure.suptitle(suptitle)
    figure.tight_layout()
    return figure


__all__ = [
    "PlacedCrystal",
    "WorldScene3D",
    "render_variant_contact_sheet",
    "render_world_scene_3d",
    "resolve_transformation_variant",
]
