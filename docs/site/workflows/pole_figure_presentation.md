# Presenting Pole Figures: Contours, Scales, And Comparison Plates

A pole figure that will be shown to somebody is a different object from one being examined. The
examining view is a scatter or a histogram: it shows where the data actually are, including the gaps
a contour would smooth over. The presenting view is a contoured density with levels a reader can
name, drawn on a scale shared with every figure it will be compared against, and carrying enough
identification to survive being cut out and pasted into a montage.

This page covers the second. The inversion and reconstruction that produce the figures are in
{doc}`texture_odf_inversion` and {doc}`harmonic_odf_reconstruction`.

## What The Drawing Is Made Of

```python
from pytex.plotting import (
    ContourSpec,
    PoleFigureSet,
    PoleFigureStyle,
    plot_pole_figure_comparison,
    plot_pole_figure_contours,
)
```

Three settings objects, so that every choice is declared rather than defaulted silently:

| Object | Governs | The failure it prevents |
| --- | --- | --- |
| `ContourSpec` | levels, level scale, range, filled bands, line labels | levels that change with the data, so no two figures can be compared and none can be quoted |
| `PoleFigureStyle` | projection, raster resolution, smoothing halfwidth, annotation, rotation | a figure whose smoothing was chosen until it looked right |
| `PoleFigureSet` | several figures on one scale, each identified | a plate on which a weak texture and a strong one look alike |

## The Density Is Estimated On The Sphere

`pole_figure_density_grid` builds the raster of the drawing, inverse-projects each raster point to
the direction it stands for with
{func}`~pytex.core.sphere.unproject_plane_points`, and evaluates the density there with
{meth}`~pytex.texture.PoleFigure.density_on_directions`. Nothing is binned into drawing pixels and
nothing is smoothed in the projection plane, where the distortion is largest exactly where pole
figures are most crowded.

The calibration is an identity rather than a tolerance: an isotropic specimen contours flat at
1 multiple of a random distribution, which `tests/unit/test_pole_figure_plotting.py` checks
directly. Points outside the projection boundary are NaN, so the contour stops at the rim.

## Choosing Levels

```python
ContourSpec(scale="linear", count=8)                          # equal steps
ContourSpec(scale="geometric", count=6, vmin=0.25)            # equal ratios
ContourSpec(scale="explicit", values=(0.5, 1.0, 2.0, 4.0, 8.0))
```

- **linear** — equal steps across the range. The default, and the right choice for a weak or
  moderate texture.
- **geometric** — equal ratios. A sharp texture spans a large dynamic range, and a linear ladder
  spends most of its levels on the empty part of it. A ratio ladder needs a strictly positive lower
  bound; asking for one over a range that reaches zero raises rather than silently shifting the
  range, because a shifted range misstates every level on the figure.
- **explicit** — the levels verbatim, which is how a paper quotes them and how a plate is matched
  to one published earlier.

`include_random_level` (on by default) puts a contour at exactly 1 m.r.d. It is the only level on a
pole figure with an absolute meaning — the boundary between orientations that are over-represented
and under-represented — and a generated ladder that stepped over it would hide that boundary. When
the generated ladder already has a level close to 1, that level is *replaced* rather than joined,
because a ladder carrying both 1.00 and 1.07 draws a band nobody can see.

`vmin` and `vmax` fix the range. Leaving them unset makes the figure self-scaled, which is correct
for a figure standing alone and wrong for one in a set.

## Comparing Samples On One Scale

The comparison case is a plate: samples down the rows, poles across the columns, one scale
throughout.

```python
plate = PoleFigureSet(
    figures=(cold_rolled_figures, annealed_figures, recrystallized_figures),
    sample_labels=("cold rolled", "annealed 400 C", "recrystallized"),
    style=PoleFigureStyle(halfwidth_deg=12.0, rotation_deg=90.0),
)
figure = plot_pole_figure_comparison(plate, suptitle="Rolling texture evolution")
```

Every panel is contoured at the same levels, computed from the pooled densities of the whole plate;
each panel is labelled with its sample identifier on the drawing itself, so a reader can attribute
a panel without counting rows; and the plate carries **one** colour bar, because a bar per panel
would suggest each panel had its own scale.

`shared_scale=False` gives each panel its own ladder. That is occasionally what a reader wants, and
it is never a comparison — the plate then also reverts to a colour bar per panel, so the drawing
says which it is.

Rows must be equal in length, so that a column of the plate means one pole; a ragged set raises at
construction.

## What Is Annotated, And Why

The publication defaults draw no Cartesian axes and no grid — they carry no information on a disc
— and instead put on the figure the things a reader needs to interpret it:

- the specimen axes where they meet the rim, taken from the specimen frame's own axis names;
- the ``{hkl}`` family, in the notation that matches what the figure actually plots (a family in
  braces, or a single plane in parentheses when family expansion was switched off);
- the maximum and minimum density, because a contoured figure without them cannot be read
  quantitatively at all;
- a marker at the strongest point;
- the sample identifier, when one is given.

## The Rotation, And What PyTex Will Not Do Silently

`rotation_deg` rotates the drawing in its own plane. Zero — the default — draws the specimen
frame's first axis to the right, which is what the projection actually does. Setting it to 90 puts
that axis at the top, which is the usual rolling-plane presentation with RD up.

PyTex rotates; it does not mirror. The familiar RD-up/TD-right layout is a *reflection* of the
projection of a right-handed specimen frame, not a rotation of it, and applying one silently would
reverse the sense of every asymmetric feature on the figure — a shear texture would lean the wrong
way. The rim labels rotate with the data, so the drawing states where the axes are whatever
rotation is chosen.

## See Also

- {doc}`texture_odf_inversion` — where the pole figures being drawn come from.
- {doc}`harmonic_odf_reconstruction` — the harmonic route, and the ghost correction that should
  precede quoting a maximum density from a PF-derived ODF.
- {doc}`stereographic_projections` — the projection itself, and the angle-preserving alternative.
- {doc}`style_customization` — themes and rcParams, which apply to these figures like any other.
- {doc}`plotting_primitives` — the figure-spec layer these builders emit.
