# YAML Style Customization

PyTex now keeps runtime plotting style in a shared YAML-driven system rather than scattering figure constants across plotting routines and notebooks.

## Scope

- built-in themes such as `journal`, `presentation`, and `dark`
- YAML loading from user files
- deterministic merge order from base defaults to runtime overrides
- shared application across texture, diffraction, and crystal-visualization plots

## Merge Order

PyTex resolves style in the following order:

1. built-in `base` theme
2. selected named theme such as `journal`
3. optional user YAML file
4. optional runtime override mapping

This lets tutorials, notebooks, and user scripts share the same plotting contract while still allowing targeted customization.

## Example

```python
from pathlib import Path

from pytex import list_style_themes, load_style_theme, resolve_style

print(list_style_themes())
print(load_style_theme("journal")["xrd"]["line_color"])

Path("custom_style.yaml").write_text(
    """
common:
  figure:
    facecolor: "#fff7ed"
xrd:
  line_color: "#9a3412"
  fill_color: "#fdba74"
""".strip()
    + "\n",
    encoding="utf-8",
)

style = resolve_style(
    theme="journal",
    style_path="custom_style.yaml",
    overrides={"xrd": {"annotate_peaks": False}},
)
print(style["common"]["figure"]["facecolor"])
print(style["xrd"]["annotate_peaks"])
```

## Design Notes

- style validation happens at load time, not inside each plotting routine
- computation code should not carry hard-coded visual policy
- runtime plots remain ordinary Matplotlib figures even though repo-tracked documentation figures remain SVG assets

## Related Material

- {doc}`plotting_primitives`
- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`xrd_generation`
- {doc}`saed_generation`
- {doc}`crystal_visualization`

## References

### Normative

- `../../standards/documentation_architecture.md`
- `../../standards/latex_and_figures.md`

### Informative

- `../../testing/plotting_validation_matrix.md`
