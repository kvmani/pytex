# Parent-Grain Reconstruction And Variant Selection

**Surface:** `pytex.core.parent_reconstruction.reconstruct_parent_orientation`,
`select_variants`, `ParentReconstructionConfig`, `ParentReconstructionReport`,
`VariantSelectionReport`, `OrientationRelationshipCatalog` and the
`standard_*_relationships` constructors, resting on
`pytex.experimental.phase_transformation.score_parent_orientations`; exposed as
the workbench operation `ebsd.or_from_grains`.

```{admonition} Stability
:class: warning

Parent reconstruction is **experimental at map scale**. The scoring primitive is
bounded and tested; what is not settled is clustering a whole map into prior
grains. This page documents what exists — scoring candidates you supply — and is
explicit about what it does not do.
```

A steel transformed from austenite retains no austenite to measure, yet the
prior austenite grain structure governs toughness. The martensite it left behind
carries the information, because each martensite variant is related to its
parent by a known orientation relationship. Recovering the parent is inverting
that relationship — and the inversion is not unique, which is the whole
difficulty.

## 1. The forward relation

An orientation relationship $R$ plus the parent symmetry generates a set of
**variants**: symmetry-equivalent ways the child can sit on the parent,

$$
\mathbf{V}_i = R\,S_{p,i},\qquad S_{p,i} \in G_p,
$$

so a child orientation is $\mathbf{C}_i = \mathbf{P}\,\mathbf{V}_i^{\mathsf{T}}$
in the canonical crystal-to-specimen convention. Kurdjumov-Sachs on fcc→bcc
gives 24 variants, Nishiyama-Wassermann 12, Burgers on bcc→hcp 12.

**Reconstruction is the inverse**: given several child orientations known to
share a parent, find $\mathbf{P}$.

## 2. Why it is ambiguous, and why that is reported rather than hidden

One child orientation with 24 possible variants gives **24 candidate parents**,
all exactly consistent with the observation. Nothing in a single measurement
distinguishes them.

Adding children helps only when they came from *different* variants: two
children of the same variant are as ambiguous as one. So the resolving power
depends not on how many children were measured but on **how many distinct
variants they sample** — which is precisely what variant selection reduces.
Strong variant selection therefore makes reconstruction harder, and it is common
in exactly the materials where reconstruction is wanted.

`ParentReconstructionReport` carries `is_ambiguous` and `ambiguous_indices`, not
just a best answer. **`ambiguity_tolerance_deg` set to zero does not remove the
ambiguity; it hides it**, and the parameter's documentation says so at the point
of use.

## 3. The scoring algorithm

The implemented surface **scores candidates you supply**. It does not search
orientation space, so the candidate set bounds what can be found.

```text
input : record (child orientations + the OR), candidate parents P_m,
        reduction, symmetry_aware, ambiguity_tolerance

1  validate the candidates against the record's parent frame, phase and symmetry
2  for each candidate P_m, predict every child it implies under the OR's variants
3  per (candidate m, child n): the crystal-frame relative rotation
       M = C_observed^T C_predicted
4  angle of M, symmetry-reduced to a disorientation when symmetry_aware
5  reduce the per-child residuals to one score per candidate:
       "mean"   -- the default
       "median" -- robust to a few misindexed children
       "max"    -- worst case
6  best = argmin score
7  ambiguous = every candidate scoring within ambiguity_tolerance of the best
```

Steps 2-4 are one `einsum` over candidates and children together, so scoring a
few hundred candidates against a few hundred children is a single vectorised
operation.

### 3.1 The reduction is a modelling choice

| Reduction | Use when |
| --- | --- |
| `mean` | children are all trustworthy; the default |
| `median` | a few children may be misindexed — one bad child cannot dominate |
| `max` | every child must be explained; conservative |

`median` is the one to reach for on real EBSD data, where a small fraction of
misindexed points is normal and a mean residual is pulled by them.

### 3.2 Symmetry awareness

`symmetry_aware=True` compares observed and predicted children by
**disorientation** rather than by raw rotation angle — the same reduction as in
{doc}`misorientation_and_disorientation`. Without it, a correct parent scores
badly whenever the comparison lands on a non-minimal representative, and the
ranking is noise.

## 4. Variant selection

`select_variants` answers the complementary question: given a parent and an OR,
**which variant does each child correspond to**?

```text
per child: the variant whose predicted orientation is nearest, by
           symmetry-reduced disorientation
```

The distribution of the resulting indices is the observable of interest:

- **Uniform** across variants — no selection; the transformation sampled the
  variants freely.
- **Skewed** — selection, pointing to stress during transformation, to
  boundary-nucleation effects, or to prior deformation.

`VariantSelectionReport` returns **one-based, strictly positive** variant
indices, matching the literature's numbering rather than a zero-based array
convention that would silently disagree with every published variant table.

Variant selection is where reconstruction and the OR literature meet: the
packet and block structure of lath martensite is a statement about which
variants group together, and {doc}`variant_correspondence` carries the
plane-and-direction side of it.

## 5. The relationship catalogue

`OrientationRelationshipCatalog` and the `standard_*_relationships` constructors
supply the named relationships per transformation system:

| Constructor | System | Contains |
| --- | --- | --- |
| `standard_fcc_bcc_relationships` | austenite → ferrite/martensite | Kurdjumov-Sachs, Nishiyama-Wassermann, Bain, Pitsch |
| `standard_bcc_hcp_relationships` | β → α (Ti, Zr) | Burgers |
| `standard_fcc_hcp_relationships` | fcc → hcp | Shoji-Nishiyama |
| `standard_cubic_cubic_relationships` | precipitates, twins | cube-on-cube, the coherent twin |
| `standard_ferrite_cementite_relationships` | steel carbides | Bagaryatskii, Pitsch-Petch, Isaichev |

Supplying a catalogue rather than a single relationship lets a reconstruction be
scored against several hypotheses, which is the honest procedure when the
operative OR is itself uncertain. To *determine* the relationship from the data
instead, see {doc}`orientation_relationship_determination`.

## 6. What is not implemented

Stated plainly, because the gap between this and a published reconstruction
workflow is the interesting part:

- **No orientation-space search.** Candidates are scored, not found. In
  practice candidates are generated by applying the inverse variants to measured
  child orientations, which is a small set, not a search.
- **No map-scale clustering.** Grouping the points of a map into prior grains —
  the graph-clustering step that makes reconstruction a microstructural tool —
  is not implemented. This is the substantive missing piece and the reason for
  the experimental marking.
- **No sub-block or packet grouping** from the reconstruction itself.
- **No confidence beyond the ambiguity set.** The report says which candidates
  are indistinguishable; it does not attach a probability to the winner.

## Verification

- Round trip: generate children from a known parent through an OR's variants,
  reconstruct, and recover the parent — with the ambiguity set reported, in
  {doc}`../examples/generated/transformation`.

## See also

- {doc}`orientation_relationship_determination` — determining $R$ from measured
  pairs, when it is not known in advance.
- {doc}`variant_correspondence` — the plane and direction correspondence of each
  variant.
- {doc}`misorientation_and_disorientation` — the reduction the scoring depends on.
- {doc}`../theory/experimental_parent_candidate_scoring` — the scoring
  primitive's derivation and limits.
- {doc}`../theory/phase_transformation_relationship_construction` — how a
  relationship is built.

## References

### Normative

- Kurdjumov, G. & Sachs, G. (1930). Über den Mechanismus der Stahlhärtung.
  *Zeitschrift für Physik* **64**, 325-343.
  <https://doi.org/10.1007/BF01397346>
- Burgers, W. G. (1934). On the process of transition of the cubic-body-centred
  modification into the hexagonal-close-packed modification of zirconium.
  *Physica* **1**, 561-586.
  <https://doi.org/10.1016/S0031-8914(34)80244-3>

### Informative

- Miyamoto, G., Takayama, N. & Furuhara, T. (2009). Accurate measurement of the
  orientation relationship of lath martensite and bainite by electron
  backscatter diffraction analysis. *Scripta Materialia* **60**, 1113-1116.
  <https://doi.org/10.1016/j.scriptamat.2009.02.053>
- Nyyssönen, T., Isakov, M., Peura, P. & Kuokkala, V.-T. (2016). Iterative
  determination of the orientation relationship between austenite and
  martensite from a large amount of grain pair misorientations. *Metallurgical
  and Materials Transactions A* **47**, 2587-2590.
  <https://doi.org/10.1007/s11661-016-3462-2>
