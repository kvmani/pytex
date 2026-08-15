# Multiphase EBSD Graph Workflows

This note records the current PyTex algorithmic contract for graph-backed and multiphase EBSD map
workflows.

## Scope

The current implementation supports:

- regular-grid adjacency through deterministic four- and eight-neighbor construction,
- staggered hexagonal-grid adjacency through immutable logical row lengths and six-neighbor
  construction,
- irregular coordinate graphs through a distance-based neighborhood radius inferred from step size or the minimum positive coordinate spacing,
- phase-resolved crystal maps with explicit per-point phase identifiers, and
- phase-aware KAM and segmentation that preserve topology but exclude cross-phase misorientation pairs from within-phase metrics.

## Neighbor Doctrine

For a crystal map with point set $\{x_i\}$, PyTex constructs an undirected neighbor graph
$\mathcal{G} = (\mathcal{V}, \mathcal{E})$.

- On rectangular grids, $\mathcal{E}$ is defined by four- or eight-connectivity and neighbor order.
- On hexagonal grids, each point joins its horizontal neighbors and the nearest points in the
  staggered rows above and below. Interior points therefore have six first-shell neighbors.
  Higher orders are cumulative graph distance, not a Euclidean-radius approximation.
- On irregular coordinates, $\mathcal{E}$ is defined by a radius threshold derived from the base spacing and the requested order.

This keeps neighborhood semantics explicit and allows KAM, grain segmentation, and boundary
extraction to share the same adjacency substrate.

## Ragged Hexagonal Rows

An EDAX/TSL `HexGrid` scan alternates `NCOLS_ODD` and `NCOLS_EVEN`. PyTex stores those counts as
`CrystalMap.row_lengths` with `grid_kind="hexagonal"`; it does not place the scan into a
rectangular `grid_shape`, because the padded positions are not measurements. Numerical local
metrics consequently return one value per measured point. Display helpers may return a padded
array, using `NaN` for scalar properties and `-1` for grain labels so the absence cannot be
mistaken for data.

For the analytic 3/2/3-row fixture, the first shell has five within-row edges and four edges
across each row boundary, hence $5+4+4=13$ unique pairs. That count and a planted six-degree
orientation perturbation are pinned by the executable EBSD example.

## Phase Masking

For phase-resolved maps, each point carries a phase identifier $p_i$. PyTex preserves the full
neighbor graph but restricts within-phase angular metrics to pairs satisfying

$$
p_i = p_j
$$

This means the geometry remains visible while phase-incompatible misorientation calculations are not
silently mixed into KAM or grain-union logic.

## Texture Extraction

Texture outputs from a multiphase crystal map require an explicit phase selector. PyTex does not
collapse a multiphase EBSD map into a single texture object without an explicit phase decision by
the user.

## Current Limits

- Irregular-coordinate graph construction is geometric rather than confidence-weighted.
- Majority-vote smoothing remains a regular-grid label operation.
- Curvature/GND finite differences and pixel-face grain perimeters remain rectangular-only. A
  hexagonal center lattice does not by itself declare the cell-boundary model those quantities
  require.
- Multiphase EBSD import normalization is now supported, but detector-pattern semantics are still outside the current normalized dataset contract.

## Informative Reference

- [MTEX Gridded EBSD Data](https://mtex-toolbox.github.io/EBSDGrid.html) distinguishes `EBSDhex`
  from square grids and shows why resampling a hexagonal scan onto squares introduces distortion.
