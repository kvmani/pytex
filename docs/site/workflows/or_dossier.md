# The Orientation-Relationship Dossier

An orientation relationship is reported in the literature as a scatter of quantities across a paper
— a cell table here, an axis/angle there, a variant count in a caption, a parallelism in the
abstract. Reproducing such a report means recomputing all of it from a sentence.

`pytex.or_dossier(relationship, ...)` turns the declaration into one object that carries every
number the declaration implies, states the convention each was computed under, explains itself in
prose, serialises against a published schema, and writes itself to disk as a bundle another person
can check.

```python
from pytex import or_dossier

dossier = or_dossier(ks, variant=17)
print(dossier.describe())
dossier.export("kurdjumov_sachs_v17")
```

## What Is In It

Five blocks, following F17 of the
{doc}`OR-analysis foundation <../architecture/orientation_relationship_analysis_foundation>`:

| Block | Contents |
| --- | --- |
| **Lattice** (one per phase) | cell parameters, cell volume, direct and reciprocal structure matrices, both metric tensors, point group |
| **Transformation** | the parent-to-child rotation, the direction correspondence carrying parent `[uvw]` to child `[uvw]`, the plane correspondence carrying parent `(hkl)` to child `(hkl)`, and the lattice-correspondence deformation with its principal strains, volume change and residual polar rotation |
| **Misorientation** | the symmetry-reduced disorientation representative, the variant count, the packet grouping, and the distinct intervariant spectrum |
| **Parallelism** | the chosen variant's own defining parallelisms, plus any near-parallelisms found among nominated families, in publication notation |
| **Figures** | written by `export`, not held in memory: the OR stereogram and the variant contact sheet |

`describe()` produces the prose; `to_json()` produces a dict against
`schemas/or_dossier.schema.json`; `export(directory)` writes the whole bundle — `or_dossier.json`,
`describe.md`, `parallelisms.csv` and `.md`, `intervariant_angles.csv`, and the two SVG figures.

## The Rule It Is Built On

**The dossier calls the existing functions and never reimplements them.** The cell volume comes from
`Lattice.volume_angstrom3`, the correspondence matrices from
`OrientationRelationship.correspondence_direct` / `correspondence_reciprocal`, the strain from
`deformation_gradient`, the spectrum from `intervariant_misorientation_angles_deg`, the packets from
`variant_close_packed_groups`, the parallelisms from the variant itself and from
`find_parallel_planes` / `find_parallel_directions`, and the figures from `plot_or_stereogram` and
`render_variant_contact_sheet`.

That rule is what makes the bundle worth anything. A dossier number that disagreed with the function
a reader would check it against is the exact class of defect this repository exists to prevent, so
the tests in `tests/unit/test_or_dossier.py` assert the agreement value by value rather than
comparing against recorded output.

## Two Things To Read Precisely

**The defining parallelisms are the variant's own.** Each variant holds a *different* member of the
parent family parallel, so `or_dossier(ks, variant=17)` reports variant 17's plane and not variant
1's. Over the 24 Kurdjumov-Sachs variants the parent plane takes four values, one per packet.

**The discovered deviations are rationalization residuals.** Nominating a family through `planes=`
or `directions=` runs the parallelism search over it. The exact child image of a parent plane is
parallel to it *by construction*, so what that search reports — and what the dossier's
`deviation_deg` column carries for a `discovered` row — is the angle by which the nearest low-index
child index misses that exact image. A small tolerance therefore keeps the pairs for which a
low-index child object really is parallel, and drops the parent members for which none is. The
`origin` column exists so the two kinds of row are never read the same way: a `defining` row's
deviation is zero by construction and measures nothing.

## What Is Not In It

The **interface block** is not implemented: no habit plane, no misfit, no terrace decomposition.
`to_json()` carries `"interface": null` and `describe()` says so in a sentence, rather than omitting
a section — a reader must be able to tell "not analysed" from "analysed and empty".

The figure bundle omits the variant pole figure and the per-variant SAED patterns. Both are
available elsewhere in PyTex, and both need something a relationship does not supply: a measured
parent orientation for the pole figure, a diffraction geometry for the patterns.

## Verification

The worked example {doc}`../examples/generated/transformation` computes six values from one dossier:
the difference between the volume the dossier reports and the volume the lattice reports (an
identity, exactly zero), the difference between that volume and the cube of the cubic edge (an
identity), and then the published Kurdjumov-Sachs figures — 24 variants, 4 packets, 10 distinct
intervariant angles, the largest being the 60° Σ3 twin relation.

## Related Material

- {doc}`../concepts/visualization_primitives` — the composite scenes and the variant-aware overlays
- {doc}`stereographic_projections` — the OR stereogram the bundle draws
- {doc}`phase_transformation_manifests_and_scoring` — the measured-data side of the same subsystem
- {doc}`workbench_application` — the composite views in the application

## References

### Normative

- International Tables for Crystallography, Volume A (cell conventions and point groups).
- `docs/standards/notation_and_conventions.md`
- `docs/standards/data_contracts_and_manifests.md`

### Informative

- Morito, Tanaka, Konishi, Furuhara and Maki, *Acta Materialia* **51** (2003) 1789 — packet
  structure and the intervariant table.
- Kurdjumov and Sachs, *Z. Phys.* **64** (1930) 325.
