# Core Model

PyTex treats scientific meaning as part of the API surface. Stable operations should not require a user to remember hidden frame, basis, or symmetry assumptions.

## Design Goal

The core model exists to make incorrect combinations difficult and correct combinations obvious. A stable scientific object should answer three questions immediately:

- what coordinate frame does this live in?
- what symmetry assumptions are active?
- what provenance or conversion context should stay attached?

## Key Rules

- `ReferenceFrame` and `FrameTransform` make coordinate meaning explicit.
- `SymmetrySpec` owns reusable symmetry operators and sector or reduction rules.
- `Lattice`, `Phase`, `CrystalPlane`, and related types prevent direct or reciprocal ambiguity.
- `Rotation`, `Orientation`, and `OrientationSet` carry mathematically precise orientation semantics instead of exposing anonymous arrays.
- `ProvenanceRecord` keeps import and conversion context attached to scientific objects.

## Why This Matters

Texture and diffraction workflows frequently fail at tool boundaries. The common failure mode is not linear algebra itself; it is semantic drift. The same numeric array can mean:

- crystal to specimen
- specimen to crystal
- active rotation
- passive coordinate change
- Bunge angles in degrees
- another Euler convention in radians

PyTex normalizes data once at the boundary and keeps later computations on PyTex-native domain objects.

![Reference Frames](../../figures/reference_frames_vectors.svg)

## Reading The Figure

The figure above fixes the repository-wide frame vocabulary.

- The crystal frame is phase-attached and uses lattice-linked crystal axes.
- The specimen frame is the target frame for orientations in texture and EBSD workflows.
- The map frame carries scan coordinates for regular or irregular measurement layouts.
- The detector frame carries image-plane or diffraction-geometry semantics rather than sample-position semantics.

That separation is deliberate. It prevents later algorithms from inventing private coordinate systems that silently disagree.

## Practical Consequences For Users

- Construct `ReferenceFrame` objects early and reuse them.
- Construct `SymmetrySpec` once per phase or orientation family instead of repeating ad hoc symmetry logic in each workflow.
- Prefer `Orientation` or `OrientationSet` over bare matrices when a value crosses module boundaries.
- Treat provenance as retained scientific metadata, not merely import bookkeeping.

```{note}
The core model is not “extra structure around NumPy arrays.” It is the mechanism that lets later EBSD, texture, and diffraction algorithms remain interpretable, testable, and interoperable.
```

## Related Material

- `docs/architecture/canonical_data_model.md`
- [../../tex/theory/canonical_data_model.tex](../../tex/theory/canonical_data_model.tex)
- [../../tex/theory/reference_frames.tex](../../tex/theory/reference_frames.tex)
- [../../figures/reference_frames_vectors.svg](../../figures/reference_frames_vectors.svg)

## References

### Normative

- `docs/standards/notation_and_conventions.md`
- `docs/architecture/canonical_data_model.md`

### Informative

- [IUCr: Crystallography overview](https://www.iucr.org/publ/50yearsofxraydiffraction/full-text/crystallography)
- Bunge, *Texture Analysis in Materials Science* (1982)
