# Stereographic Projections and XRDML Texture Import in PyTex

## Scope

This note records the current scientific assumptions behind two newly linked PyTex surfaces:

- stereographic or equal-area projection of crystallographic directions, plane traces, and rotational symmetry axes,
- XRDML pole-figure import and its handoff into the discrete kernel ODF inversion model.

## Spherical Projection Convention

PyTex treats specimen directions as right-handed Cartesian unit vectors with the positive specimen $z$ axis at the center of the projection. For a unit vector $\mathbf{v} = (x, y, z)$ in the upper hemisphere, the current projection formulas are

$$
\mathbf{p}_{\text{stereo}} = \left(\frac{x}{1+z}, \frac{y}{1+z}\right),
  \qquad
  \mathbf{p}_{\text{equal-area}} = \left(x\sqrt{\frac{2}{1+z}}, y\sqrt{\frac{2}{1+z}}\right)
$$

For antipodal surfaces, lower-hemisphere directions are first folded onto the upper hemisphere before projection. Great-circle plane traces are generated explicitly from the plane normal and then projected as sampled curves rather than retrieved from static image templates.

## Rotational Symmetry Elements

The current symmetry-element plotting surface covers proper rotational axes only. Each non-identity operator contributes an axis and an order inferred from its rotation angle. Duplicate axes are consolidated by keeping the highest order present on that axis. The rendered symbols are therefore semantic summaries of rotational order, not exhaustive lists of every generator or power of a generator.

## XRDML Pole-Figure Import

The implemented XRDML adapter targets pole-figure style measurements carrying `Phi` and `Psi` goniometer positions. These angles are interpreted as specimen spherical coordinates, with $\Psi$ taken as the polar angle from the specimen $z$ axis and $\Phi$ as the azimuth around that axis.

The adapter normalizes the XML payload into a PyTex `PoleFigure` through three steps:

1. read scan intensities or count-time normalized counts,
2. reconstruct explicit specimen directions from the `Phi`/`Psi` scan geometry,
3. require the caller to provide the crystallographic pole explicitly because the reflection indices are not reliably encoded in a stable machine-readable field across practical XRDML examples.

## ODF Reconstruction Path

Once imported, XRDML-backed pole figures use the same forward model as the in-memory discrete inversion surface:

$$
A_{mj} = \frac{1}{|\mathcal{H}|}\sum_{h \in \mathcal{H}} K\!\left(\angle(\mathbf{s}_m, g_j h)\right)
$$

with measurement directions $\mathbf{s}_m$, orientation dictionary $\{g_j\}$, pole family $\mathcal{H}$, and kernel $K$. PyTex therefore validates the XML import boundary separately from the inversion doctrine: a pinned real XRDML file checks parser semantics, while deterministic synthetic fixtures check inversion behavior.

## Current Limits

The current surface does not yet claim:

- full MTEX visual parity for all stereographic render details,
- mirror-plane, inversion-center, or nonsymmorphic symmetry-element symbol coverage,
- full harmonic or experimentally corrected X-ray texture inversion,
- robust inference of reflection metadata from free-text XRDML comments.

## Normative References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, 1969. DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- MTEX documentation, “Spherical Projections” and related pole-figure import notes. <https://mtex-toolbox.github.io/>
- D. Kriegner and contributors, *xrayutilities*, <https://github.com/dkriegner/xrayutilities>.
