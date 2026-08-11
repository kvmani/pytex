# Multiphase EBSD Graph Workflows

This note records the current PyTex algorithmic contract for graph-backed and multiphase EBSD map
workflows.

## Scope

The current implementation supports:

- regular-grid adjacency through deterministic four- and eight-neighbor construction,
- irregular coordinate graphs through a distance-based neighborhood radius inferred from step size or the minimum positive coordinate spacing,
- phase-resolved crystal maps with explicit per-point phase identifiers, and
- phase-aware KAM and segmentation that preserve topology but exclude cross-phase misorientation pairs from within-phase metrics.

## Neighbor Doctrine

For a crystal map with point set $\{x_i\}$, PyTex constructs an undirected neighbor graph
$\mathcal{G} = (\mathcal{V}, \mathcal{E})$.

- On regular grids, $\mathcal{E}$ is defined by the requested connectivity and neighbor order.
- On irregular coordinates, $\mathcal{E}$ is defined by a radius threshold derived from the base spacing and the requested order.

This keeps neighborhood semantics explicit and allows KAM, grain segmentation, and boundary
extraction to share the same adjacency substrate.

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
- Multiphase EBSD import normalization is now supported, but detector-pattern semantics are still outside the current normalized dataset contract.
