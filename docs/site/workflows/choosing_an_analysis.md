# Choose An Analysis From The Measurement

Start with the question and the data you actually measured. A diffraction pattern, an
orientation map and a list of grain means support different inferences, even when they
refer to the same specimen. This guide connects those inputs to a useful first result
and the next diagnostic to inspect.

For installation and a first calculation, use {doc}`../tutorials/quickstart`. For the
graphical interface, use {doc}`workbench_application`. Numerical examples linked below
are executed and compared with independent reference values during testing.

## Find The Relevant Workflow

| Your question | Required input | Useful first output | What to inspect before interpreting it |
| --- | --- | --- | --- |
| Which crystallographic plane or direction is this? | Phase lattice, crystal frame and indices | {doc}`vectorized_miller_workflows`: spacing, vectors, angles and families | Distinguish a direct-lattice direction from a reciprocal-lattice plane normal. In a non-cubic cell, equal index triples need not be parallel. |
| What relationship connects the parent and product grains? | Row-matched parent/product orientations, both phases, a shared specimen frame and Euler convention | {doc}`workbench_application`: **OR from grains**, fitted rotation, catalog distances, parallelisms and pair residuals | Check scatter and convergence, inspect suspect pairs, and compare the catalog margin with the residuals. A single pair can match a name but cannot establish population consistency. |
| What variants and diffraction coincidences should that relationship produce? | A specified OR, parent orientation and both phases | {doc}`or_dossier` and {doc}`composite_or_diffraction`: variant correspondence and composite patterns | Distinguish a catalog relationship from a measured fit and a rationalized approximation. Preserve the cost of replacing the fit by integer parallelisms. |
| Where does local orientation variation occur? | A normalized EBSD map with spatial calibration and phase assignments | {doc}`ebsd_kam` and {doc}`ebsd_grains`: local misorientation, segmentation and grain metrics | Compare the neighborhood size with map resolution; document boundary exclusion and cleanup choices. Missing or unindexed pixels are not zero misorientation. |
| What texture does this map represent? | Indexed orientations, phase selection and a sampling choice | {doc}`ebsd_to_texture_outputs`: PF, IPF and ODF | Decide whether pixels or grains are the sampling units. Repeating many pixels from one grain changes its influence; it does not add independent grains. |
| What texture explains measured pole figures? | Reflection identities, measurement directions and intensities, with background and defocusing information | {doc}`xrdml_texture_import` or {doc}`labotex_open_pole_figures`, then {doc}`texture_odf_inversion` | Inspect coverage and reconstructed pole-figure residuals. An ODF fitting incomplete data is not uniquely established by that fit. |
| Which supplied phase best explains a powder scan? | Measured angular scan, radiation and candidate structures | {doc}`xrd_generation`: ranked phase-identification evidence | Inspect unexplained peaks, missing predicted lines and the fitted cell scale. A ranking among supplied candidates is not a database search or a multiphase quantitative refinement. |
| What zone axis or orientation explains a TEM pattern? | Measured spots, candidate phase and detector calibration | {doc}`tem_pattern_indexing` and {doc}`saed_pattern_solving`: indexed spots and geometric residuals | Inspect calibration, unassigned spots and possible multiple patterns. Geometric agreement does not establish a dynamical intensity model. |
| How will the microscope optics change HREM contrast? | Atomic snapshot, beam energy, aberrations, coherence and sampling | {doc}`../theory/hrem_multislice_and_ctf`: transfer function and simulated image | Record whether the phase-object or abTEM multislice backend ran. Compare sampling and model assumptions before comparing simulated and measured contrast. |
| How does orientation affect a material response? | Orientations and the relevant single-crystal stiffness or slip-system assumptions | {doc}`../examples/generated/elastic-anisotropy` and {doc}`../examples/generated/schmid-and-taylor` | State the imposed load/strain and homogenization assumptions; orientation geometry alone does not supply constitutive parameters. |

## Match The Statistic To The Question

**Residuals measure agreement with a model.** They are not automatically measurement
uncertainties. Inspect individual residuals as well as their mean: a small mean can conceal
one mismatched grain pair or one unexplained diffraction peak. Retain rejected observations
and the reason for rejecting them in the experimental record.

**Weights describe influence.** Equal weighting of pixels, equal weighting of grains and
weighting by independently assessed precision answer different questions. In measured-pair
OR fitting, use the optional weight column deliberately: zero excludes a pair from the fit
while preserving its residual. The reported effective pair count describes concentration
of weights; it does not correct correlations between measurements. The
{doc}`weighted OR examples <../examples/generated/weighted-or-fitting>` show exactly what
changes when a pair is excluded.

**Directional summaries discard information.** A Kearns parameter summarizes a second
moment of the basal-pole distribution; distinct textures can share that summary. Read the
{doc}`Kearns examples <../examples/generated/kearns-parameter>` together with the pole
figure, rather than using one scalar as a complete texture description. Likewise,
{doc}`GND estimates <../theory/lattice_curvature_and_gnd_density>` from a planar EBSD map
must retain the distinction between observable components and unmeasured derivatives.

## Keep A Result Reproducible

Save the input data and phase definition, frame/convention declarations, units, selection
or exclusion decisions, numerical settings and PyTex version with the result. Prefer the
result's `describe()` prose and supported JSON/CSV exports over manually transcribing the
screen. Keep the image and numerical evidence together: a plot without the calculation
settings is insufficient for reproducing an analysis.

When two tools disagree, compare the declared frames, orientation mapping direction,
symmetry groups and reciprocal normalization before comparing arrays element by element.
Use {doc}`../concepts/reference_frames_and_conventions` for that reconciliation. Consult
{doc}`../validation/mtex_parity_matrix` for the scope of external comparison; an analytic
identity test, a synthetic recovery test and a measured external baseline provide different
evidence and should be cited as such.

## References

- {doc}`../standards/notation_and_conventions`
- {doc}`../standards/data_contracts_and_manifests`
- {doc}`../theory/orientation_relationship_determination`
- {doc}`../theory/hrem_multislice_and_ctf`
