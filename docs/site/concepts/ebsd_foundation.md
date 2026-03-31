# EBSD Foundation

This page drills into the EBSD layer as the map-and-neighborhood specialization built on the canonical orientation model.

## What The EBSD Layer Owns

- `CrystalMap` as the canonical orientation-map container
- map-frame metadata and regular-grid semantics
- neighborhood topology on 2D grids
- KAM, segmentation, GROD, boundaries, and cleanup
- grain-graph aggregation and reproducible normalization contracts
- import-manifest normalization for vendor or third-party inputs

## EBSD Flow

```{mermaid}
flowchart LR
    core["Canonical core<br/>Frames, symmetry, orientation, provenance"]
    map["CrystalMap<br/>Coordinates, orientations, map frame"]
    neighbors["Neighborhood semantics<br/>Regular-grid adjacency"]
    grains["Grains and boundaries<br/>Segmentation, GROD, cleanup"]
    manifest["Import manifests<br/>Normalization and provenance"]
    adapters["Adapters<br/>Vendor and third-party bridges"]

    core --> map
    map --> neighbors
    neighbors --> grains
    adapters --> manifest
    manifest --> map
    core --> grains

    classDef core fill:#f0f4f8,stroke:#334e68,color:#102a43,stroke-width:1.5px;
    classDef domain fill:#d9f2e6,stroke:#2d6a4f,color:#1b4332,stroke-width:1.5px;
    classDef boundary fill:#fff3d6,stroke:#b7791f,color:#744210,stroke-width:1.5px;
    class core core;
    class map,neighbors,grains domain;
    class manifest,adapters boundary;
```

## Why This Layer Matters

EBSD is where PyTex must combine scientific semantics with measurement topology.

- the map frame must stay distinct from the specimen frame unless a workflow explicitly links them
- neighborhood logic must be deterministic
- grain cleanup must remain reproducible
- normalization must preserve source-system meaning rather than flattening it away

## Current State

- `CrystalMap` exists as the canonical map container
- KAM, segmentation, GROD, boundary extraction, cleanup, and grain graphs are implemented
- stable import manifests exist
- richer vendor detector/pattern metadata normalization remains ahead

## Related Material

- {doc}`../workflows/ebsd_grains`
- {doc}`../workflows/ebsd_kam`
- {doc}`../architecture/ebsd_foundation`
- {doc}`../architecture/multimodal_characterization_foundation`

## References

### Normative

- {doc}`../architecture/ebsd_foundation`
- {doc}`../standards/reference_canon`

### Informative

- {doc}`../workflows/ebsd_import_normalization`
