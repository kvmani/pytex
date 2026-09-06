"""Generate the canonical algorithm flow-sheet figures.

These are *generated assets*, not hand-authored ones: they are produced by
`pytex.plotting.algorithm_diagrams` from a declarative description of each
algorithm, so a documentation figure shares one visual language with the
reference-frame diagrams and cannot drift into a hand-drawn approximation.
Re-run this script whenever an algorithm's stages or constraints change.

Usage::

    python scripts/generate_algorithm_figures.py

Outputs (tracked as canonical documentation assets):

- ``docs/figures/or_determination_algorithm.svg`` — determining the orientation
  relationship from measured parent/child orientations.
- ``docs/figures/variant_correspondence_algorithm.svg`` — mapping a plane or
  direction through every transformation variant.
- ``docs/figures/composite_saed_algorithm.svg`` — assembling a composite
  parent + variant zone-axis pattern on one shared detector.
- ``docs/figures/saed_indexing_algorithm.svg`` — solving a measured pattern by
  ratio/angle indexing.
- ``docs/figures/pole_figure_inversion_algorithm.svg`` — recovering an ODF from
  measured pole figures.
- ``docs/figures/rietveld_refinement_algorithm.svg`` — whole-profile powder
  refinement.
- ``docs/figures/phase_identification_algorithm.svg`` — ranking candidate
  structures against a measured powder pattern.
- ``docs/figures/ebsd_grain_metrics_algorithm.svg`` — grains, the local
  misorientation family, and GND density.
- ``docs/figures/kikuchi_geometry_algorithm.svg`` — the Kikuchi forward model.

Each figure states its stages, the constraint governing each stage, and — in the
footer — what the algorithm deliberately does not do.
"""

from __future__ import annotations

from pathlib import Path

from pytex.plotting.algorithm_diagrams import AlgorithmStage, SideNote, algorithm_flow_svg

FIGURES = Path(__file__).resolve().parents[1] / "docs" / "figures"


def or_determination_figure() -> str:
    """Flow sheet for `characterize_orientation_relationship`."""

    return algorithm_flow_svg(
        [
            (
                "1 - measurement",
                [
                    AlgorithmStage(
                        label="Paired orientations",
                        role="input",
                        formula="P_i (parent), C_i (child)",
                        detail=[
                            "n grain-mean orientations per phase,",
                            "row-matched, one specimen frame.",
                            "Euler angles accepted directly.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Per-pair rotation",
                        formula="V_i = C_i^T P_i",
                        detail=[
                            "The operative parent-to-child rotation",
                            "of each pair, in the canonical",
                            "crystal-to-specimen convention C = P V^T.",
                        ],
                    ),
                ],
            ),
            (
                "2 - fit, with no nominal relationship required",
                [
                    AlgorithmStage(
                        label="Double-coset seed",
                        formula="max-trace element of G_c V_0 G_p",
                        detail=[
                            "One pair reduced to its minimum-angle",
                            "representative. The coset absorbs the",
                            "parent symmetry operation that",
                            "distinguishes one variant from another.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Align to the estimate",
                        formula="argmax_S,S' tr(S_c V_i S_p R^T)",
                        detail=[
                            "Each measurement is replaced by its",
                            "symmetry-equivalent description",
                            "nearest the current estimate.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Quaternion eigen-mean",
                        formula="R <- top eigenvector of sum q q^T",
                        detail=[
                            "Markley's rotation average of the",
                            "aligned set. Align and average iterate",
                            "until the assignment is stable.",
                        ],
                    ),
                ],
            ),
            (
                "3 - name it, state it, judge it",
                [
                    AlgorithmStage(
                        label="Rank the catalog",
                        role="decision",
                        formula="min over G_c, G_p of angle(R, R_cand)",
                        detail=[
                            "Symmetry-reduced distance to each named",
                            "relationship for the two crystal systems.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Extract parallelisms",
                        detail=[
                            "Low-index parent planes and directions",
                            "carried through R, kept when the image",
                            "is near a low-index child object.",
                            "Preference: the winner's own families.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Report with a verdict",
                        role="output",
                        detail=[
                            "Fitted rotation, pair scatter, ranking,",
                            "the (hkl)||(hkl) and [uvw]||[uvw]",
                            "statement, and is_conclusive.",
                        ],
                    ),
                ],
            ),
        ],
        title="Determining an Orientation Relationship From Measured Orientations",
        subtitle=(
            "characterize_orientation_relationship: fit without a nominal, name it against "
            "the catalog, state it crystallographically, and say whether to trust it."
        ),
        description=(
            "Three-lane flow sheet. Lane 1 turns paired parent and child orientations into "
            "per-pair rotations V = C^T P. Lane 2 seeds the fit from the double-coset "
            "reduction of one pair and refines it by alternating symmetry alignment with the "
            "quaternion eigen-mean. Lane 3 ranks the standard catalog under both symmetry "
            "groups, extracts the parallel-plane and parallel-direction statement, and emits "
            "a report carrying an explicit conclusiveness verdict. Side notes give the "
            "constraint governing each stage."
        ),
        notes=[
            SideNote(
                stage_index=1,
                title="Constraint: one convention",
                lines=[
                    "V = C^T P has exactly one definition in",
                    "the library. Both orientation sets must",
                    "share a specimen frame and carry phases;",
                    "row order is the pairing.",
                ],
            ),
            SideNote(
                stage_index=2,
                title="Constraint: seed ONE pair",
                lines=[
                    "The maximum-trace coset element is not",
                    "unique when the relationship's rotation is",
                    "itself symmetric. Averaging every pair's",
                    "reduced representative turns Bain (45 deg,",
                    "<100>, 3 variants) into 26.9 deg, which",
                    "then reads as Kurdjumov-Sachs.",
                ],
            ),
            SideNote(
                stage_index=5,
                title="Constraint: catalog scope",
                lines=[
                    "Only the supplied candidates can win.",
                    "Cubic-to-cubic assumes the fcc->bcc class:",
                    "point-group symmetry cannot tell fcc from",
                    "bcc. Supply a catalog when that is wrong.",
                ],
            ),
            SideNote(
                stage_index=6,
                title="Constraint: several are true at once",
                lines=[
                    "A rotation satisfies many exact low-index",
                    "parallelisms. Which one the literature",
                    "quotes depends on the structures, not the",
                    "rotation, so a preference is required.",
                ],
            ),
            SideNote(
                stage_index=7,
                title="Failure mode: admitted, not hidden",
                lines=[
                    "is_conclusive needs the winner to fit within",
                    "tolerance AND to lead the runner-up by more",
                    "than the scatter and its own misfit.",
                    "Measured: conclusive to 2 deg of scatter,",
                    "inconclusive at 5 deg (KS-GT spacing 2.4).",
                ],
            ),
        ],
        footer=[
            "Validation is synthetic: planted variants of a known relationship, recovered "
            "with the relationship withheld. Measured-EBSD fixtures remain outstanding.",
            "Cubic example: Kurdjumov-Sachs, fcc -> bcc.  Hexagonal example: Burgers, "
            "bcc -> hcp, whose statement is reported in four-index Miller-Bravais form.",
        ],
    )


def variant_correspondence_figure() -> str:
    """Flow sheet for `variant_correspondence_table`."""

    return algorithm_flow_svg(
        [
            (
                "1 - the object and the variants",
                [
                    AlgorithmStage(
                        label="Nominated object",
                        role="input",
                        formula="(hkl) or [uvw] in one phase",
                        detail=[
                            "One or many, all of one kind.",
                            "Sense chooses parent->child or",
                            "child->parent.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Variant set",
                        role="input",
                        formula="V_k = R S_p,k, reduced by G_c",
                        detail=[
                            "Parent operators generate candidates;",
                            "the child-symmetry orbit defines",
                            "distinctness. KS 24, Burgers 12.",
                        ],
                    ),
                ],
            ),
            (
                "2 - map each object through each variant",
                [
                    AlgorithmStage(
                        label="Exact image",
                        formula="h_c = A_c*^-1 V_k A_p* h_p",
                        detail=[
                            "Planes ride the reciprocal basis,",
                            "directions the direct basis, so the",
                            "zone law h.u is preserved.",
                            "The image is generally irrational.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Rationalize",
                        role="decision",
                        formula="argmin_t angle(A t, A h_exact)",
                        detail=[
                            "Nearest primitive integer triple, by",
                            "the true angle between Cartesian",
                            "images. Sign-sensitive: the triple",
                            "pointing along the image wins.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Group by family",
                        formula="orbit of the image under G",
                        detail=[
                            "Variants whose images are symmetry-",
                            "equivalent share an equivalence-group",
                            "label, which is what turns 24 rows",
                            "into four distinct answers.",
                        ],
                    ),
                ],
            ),
            (
                "3 - the table",
                [
                    AlgorithmStage(
                        label="Rows with residuals",
                        role="output",
                        detail=[
                            "Per (object, variant): exact",
                            "components, integer indices, the",
                            "angle between them, labels, group.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Exact subset",
                        role="output",
                        formula="residual <= 1e-6 deg",
                        detail=[
                            "The variants that really do carry the",
                            "object onto that low-index image -",
                            "the packet, physically.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Export",
                        role="output",
                        detail=[
                            "CSV, Markdown, records, JSON, and a",
                            "describe() that states how many",
                            "distinct answers there are.",
                        ],
                    ),
                ],
            ),
        ],
        title="Variant-Resolved Plane And Direction Correspondence",
        subtitle=(
            "variant_correspondence_table: what one parent (hkl) or [uvw] becomes in every "
            "product variant, grouped so the answer is readable."
        ),
        description=(
            "Three-lane flow sheet. Lane 1 takes the nominated crystallographic object and "
            "the variant set. Lane 2 maps each object through each variant on the correct "
            "basis, rationalizes the generally irrational image to the nearest primitive "
            "integer triple by true angle, and groups variants whose images are "
            "symmetry-equivalent. Lane 3 emits the table, the exactly-parallel subset that "
            "carries the physics, and the export surfaces. Side notes give the constraint "
            "governing each stage."
        ),
        notes=[
            SideNote(
                stage_index=2,
                title="Constraint: basis by object kind",
                lines=[
                    "Using the direct basis for a plane is the",
                    "classic error. Plane indices are reciprocal-",
                    "basis components already; the library routes",
                    "every mapping through one pair of helpers so",
                    "the choice cannot be made per call site.",
                ],
            ),
            SideNote(
                stage_index=3,
                title="Constraint: the index bound",
                lines=[
                    "max_index bounds the candidate triples.",
                    "Raising it never worsens a residual, and",
                    "never changes WHICH variants are exact -",
                    "only how the irrational images are labelled,",
                    "so 'four distinct images' is partly a",
                    "bookkeeping choice. describe() says so.",
                ],
            ),
            SideNote(
                stage_index=5,
                title="Result: the packet structure",
                lines=[
                    "Cubic: KS (111) gives 4 distinct images",
                    "across 24 variants, 6 each; the 6 exact ones",
                    "are {011} at zero residual - Morito's packet.",
                    "Hexagonal: Burgers (011)beta is basal in 2",
                    "of 12 - six packets of two.",
                ],
            ),
            SideNote(
                stage_index=6,
                title="Asymmetry worth knowing",
                lines=[
                    "The reverse map is not selective: the child",
                    "(0001) maps back onto a {110} parent plane",
                    "in ALL 12 variants, because every variant's",
                    "basal plane came from some {110}.",
                ],
            ),
        ],
        footer=[
            "Rationalization is sign-sensitive by design, so the hexagonal basal image may be "
            "reported as (0001) or its antiparallel (000-1); both name the same plane.",
            "Hexagonal phases are labelled in four-index Miller-Bravais form throughout.",
        ],
    )


def composite_saed_figure() -> str:
    """Flow sheet for the composite SAED engine and its two anchoring modes."""

    return algorithm_flow_svg(
        [
            (
                "1 - choose the viewing direction",
                [
                    AlgorithmStage(
                        label="Parent zone axis",
                        role="input",
                        formula="z_p, rational [uvw]",
                        detail=[
                            "The derivation's natural choice.",
                        ],
                    ),
                    AlgorithmStage(
                        label="or a product zone axis",
                        role="input",
                        formula="z_p = R_k^T z_c",
                        detail=[
                            "The microscope's natural choice: tilt",
                            "the product on zone. The implied",
                            "parent direction is irrational and is",
                            "reported exactly plus nearest rational.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Shared detector basis",
                        formula="(u, v, z), u x v = z",
                        detail=[
                            "Built once from z_p by one function,",
                            "so both anchoring routes give",
                            "identical geometry.",
                        ],
                    ),
                ],
            ),
            (
                "2 - simulate each sub-pattern on that basis",
                [
                    AlgorithmStage(
                        label="Enumerate and filter",
                        role="decision",
                        formula="centring allowed, |g| <= g_max",
                        detail=[
                            "Reflection cube to max_index, then the",
                            "phase's lattice-centring condition.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Excitation error",
                        role="decision",
                        formula="s_g = g_z - lambda |g|^2 / 2",
                        detail=[
                            "Keep |s_g| <= tolerance. This is the",
                            "Ewald-sphere proximity test; exact",
                            "zone-axis spots have s_g = -lambda|g|^2/2.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Project and sort",
                        formula="r = (L lambda) |g_in-plane|",
                        detail=[
                            "Structure factors, max-normalized",
                            "intensity, then a sort whose",
                            "continuous keys are quantized so",
                            "ties fall to the exact hkl order.",
                        ],
                    ),
                ],
            ),
            (
                "3 - assemble and export",
                [
                    AlgorithmStage(
                        label="Child bases",
                        formula="B_k = V_k B_parent",
                        detail=[
                            "Each variant's basis is the parent's",
                            "rotated into that child frame, so all",
                            "sub-patterns overlay one detector.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Coincidences",
                        detail=[
                            "Parent/child reflection pairs within a",
                            "detector tolerance - the measurable",
                            "content of an OR in TEM practice.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Figure, table, manifest",
                        role="output",
                        detail=[
                            "One row per rendered spot, the",
                            "rendered figure, and a JSON manifest",
                            "recording how both were produced.",
                        ],
                    ),
                ],
            ),
        ],
        title="Composite Parent + Variant SAED Assembly",
        subtitle=(
            "simulate_composite_saed and simulate_composite_saed_from_child_zone: one shared "
            "detector geometry, whichever crystal's zone axis was chosen."
        ),
        description=(
            "Three-lane flow sheet. Lane 1 resolves the viewing direction from either a "
            "parent zone axis or a product-variant zone axis mapped back through the variant "
            "rotation, then builds the shared detector triad once. Lane 2 simulates each "
            "sub-pattern: enumerate reflections, apply the lattice-centring condition, select "
            "by excitation error, compute structure factors, project onto the detector and "
            "sort deterministically. Lane 3 assembles the variant bases, finds parent/child "
            "spot coincidences and exports the figure, reflection table and manifest. Side "
            "notes give the constraint governing each stage."
        ),
        notes=[
            SideNote(
                stage_index=2,
                title="Identity: the two routes agree",
                lines=[
                    "Anchoring on variant k's image of a parent",
                    "zone reproduces the parent-anchored pattern",
                    "for that zone exactly - measured to 1e-13 mm",
                    "- because both build the basis from the same",
                    "parent direction through the same function.",
                ],
            ),
            SideNote(
                stage_index=3,
                title="Constraint: declare the space group",
                lines=[
                    "Centring is read from the first letter of the",
                    "space-group symbol. A phase without one is",
                    "simulated as PRIMITIVE, keeping reflections",
                    "the real structure forbids. centering_audit()",
                    "reports declared vs assumed; describe() warns.",
                ],
            ),
            SideNote(
                stage_index=4,
                title="Constraint: zero-order Laue zone",
                lines=[
                    "The excitation-error window selects the ZOLZ",
                    "at the default tolerance. Higher-order Laue",
                    "zones and double diffraction are out of scope",
                    "for this engine.",
                ],
            ),
            SideNote(
                stage_index=5,
                title="Constraint: sort keys are quantized",
                lines=[
                    "Symmetry-equivalent reflections have equal",
                    "intensity and radius differing by ~1e-14, so",
                    "raw keys let noise decide the order and the",
                    "same pattern came out permuted. Quantized to",
                    "1 pm of radius and 1e-12 of full-scale I.",
                ],
            ),
            SideNote(
                stage_index=8,
                title="Constraint: intensities are per-pattern",
                lines=[
                    "Each sub-pattern is normalized to its own",
                    "maximum. Kinematic theory defines no intensity",
                    "ratio between two phases, so comparing across",
                    "sources is meaningless - a shared scale would",
                    "manufacture a number the theory lacks.",
                ],
            ),
        ],
        footer=[
            "Kinematic only: no dynamical (Bloch-wave) intensities, no double diffraction, no "
            "HOLZ rings. Intensities rank reflections; they do not predict a measured plate.",
            "Cubic example: Kurdjumov-Sachs down a parent <001>.  Hexagonal example: Burgers "
            "down beta [110], where the alpha basal zone appears exactly.",
        ],
    )


def saed_indexing_figure() -> str:
    """Flow sheet for `solve_saed_pattern`."""

    return algorithm_flow_svg(
        [
            (
                "1 - calibrate the picked spots",
                [
                    AlgorithmStage(
                        label="Picked positions",
                        role="input",
                        detail=[
                            "Clicked, or listed in the measured-",
                            "pattern YAML. The transmitted beam is",
                            "the centre, not a spot.",
                        ],
                    ),
                    AlgorithmStage(
                        label="To reciprocal space",
                        formula="|g| = r / (L lambda)",
                        detail=[
                            "Pixels scale by the pixel size first.",
                            "The camera constant may be given, or",
                            "derived from camera length and kV.",
                        ],
                    ),
                ],
            ),
            (
                "2 - seed from two reflections",
                [
                    AlgorithmStage(
                        label="Seed pairs",
                        detail=[
                            "The shortest non-collinear measured",
                            "vectors: best determined relative to",
                            "picking error, and two of them fix",
                            "the zone.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Ratio and angle test",
                        role="decision",
                        formula="| |g_calc| - |g_obs| | / |g_obs| <= tol",
                        detail=[
                            "Both lengths must match a calculated",
                            "allowed reflection, and the",
                            "interplanar angle must match within",
                            "the angular tolerance.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Zone and orientation",
                        formula="z = g_1 x g_2,  R = F E^T",
                        detail=[
                            "Right-handed triads from the",
                            "calculated pair and the observed",
                            "pair; their product is the",
                            "crystal-to-pattern rotation.",
                        ],
                    ),
                ],
            ),
            (
                "3 - verify, rank, and judge",
                [
                    AlgorithmStage(
                        label="Index every spot",
                        detail=[
                            "Project all allowed zone reflections",
                            "and give each measured spot its",
                            "nearest free prediction within its",
                            "match radius.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Deduplicate and canonicalize",
                        formula="R ~ R S for S in G",
                        detail=[
                            "Symmetry-equivalent solutions are one",
                            "answer; the survivor is rewritten as",
                            "the conventional description.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Ranked report",
                        role="output",
                        detail=[
                            "Matched fraction first, then residual.",
                            "Phase, zone axis, orientation, every",
                            "spot's indices, and is_conclusive.",
                        ],
                    ),
                ],
            ),
        ],
        title="Solving A Measured SAED Pattern",
        subtitle=(
            "solve_saed_pattern: classical ratio/angle indexing from a calibrated spot list. "
            "Geometry alone decides - intensities are never used."
        ),
        description=(
            "Three-lane flow sheet. Lane 1 converts picked spot positions into reciprocal-"
            "space vectors through the camera constant. Lane 2 seeds a solution from the two "
            "shortest non-collinear vectors, admits a calculated reflection pair only when "
            "both lengths and the interplanar angle match, and builds the zone axis and the "
            "crystal-to-pattern rotation from that pair. Lane 3 indexes every remaining spot "
            "by projection, deduplicates symmetry-equivalent solutions into the conventional "
            "description, and ranks the survivors with an explicit conclusiveness verdict. "
            "Side notes give the constraint governing each stage."
        ),
        notes=[
            SideNote(
                stage_index=1,
                title="Constraint: calibrate or fail early",
                lines=[
                    "Coordinates in pixels or millimetres without",
                    "a camera constant are rejected at",
                    "construction, not at the first spot: an",
                    "uncalibrated length is not recoverable.",
                ],
            ),
            SideNote(
                stage_index=3,
                title="Constraint: absences come from the phase",
                lines=[
                    "Only reflections the phase's centring allows",
                    "are offered. A phase without a space group is",
                    "treated as primitive - the same trap as in",
                    "simulation, with the same consequence.",
                ],
            ),
            SideNote(
                stage_index=4,
                title="Constraint: zone-axis pattern assumed",
                lines=[
                    "The spots must lie in one zero-order Laue",
                    "zone. A crystal tilted off zone - a variant",
                    "seen from a PARENT zone axis - is only partly",
                    "indexed, and that partial match is the honest",
                    "outcome, not a bug.",
                ],
            ),
            SideNote(
                stage_index=6,
                title="Ambiguity: intrinsic, not a failure",
                lines=[
                    "A single pattern cannot distinguish a zone",
                    "axis from its reverse for a centrosymmetric",
                    "reflection set: inverting the crystal leaves",
                    "the pattern unchanged. The report names this",
                    "rather than presenting one sense as the answer.",
                ],
            ),
            SideNote(
                stage_index=7,
                title="Failure mode: no answer is an answer",
                lines=[
                    "best() raises rather than guessing when",
                    "nothing was solved. is_conclusive is False",
                    "whenever spots are left unindexed or a",
                    "genuinely different candidate explains the",
                    "pattern equally well.",
                ],
            ),
        ],
        footer=[
            "Intensities are carried through for plotting but never used to index: a kinematic "
            "intensity model is not reliable enough, and a printed pattern rarely carries them.",
            "Cubic example: an fcc pattern down [001], where the {220}/{200} radius ratio is "
            "exactly sqrt(2) at 45 degrees.  Hexagonal example: alpha-Ti down [0001].",
        ],
    )


def pole_figure_inversion_figure() -> str:
    """Flow sheet for `ODF.invert_pole_figures` and the harmonic route."""

    return algorithm_flow_svg(
        [
            (
                "1 - the measurement, corrected before it is inverted",
                [
                    AlgorithmStage(
                        label="Measured pole figures",
                        role="input",
                        formula="P_hkl(y), in m.r.d.",
                        detail=[
                            "Several independent {hkl}; one is never",
                            "enough. Intensities on the",
                            "multiples-of-random scale.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Defocus correction",
                        detail=[
                            "Calibrated from a texture-free standard.",
                            "Uncorrected, the inversion reproduces an",
                            "instrumental rim as a texture feature.",
                        ],
                    ),
                ],
            ),
            (
                "2 - build the operator on the observations' scale",
                [
                    AlgorithmStage(
                        label="Kernel response",
                        formula="A[i,j] = sum_family psi(angle(g_j h, y_i))",
                        detail=[
                            "Pole density of dictionary orientation j",
                            "at measured direction i, summed over the",
                            "{hkl} symmetry family.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Rescale to m.r.d.",
                        role="decision",
                        formula="A <- A / random_pole_density(psi)",
                        detail=[
                            "The raw kernel sum for a random texture",
                            "is the kernel mean, not 1 - a factor of",
                            "order 64 at a 12 degree halfwidth.",
                        ],
                    ),
                ],
            ),
            (
                "3 - solve, constrained by what an ODF is",
                [
                    AlgorithmStage(
                        label="Non-negative simplex fit",
                        formula="min |Aw - b|^2 + lam|w|^2,  w>=0, sum w = 1",
                        detail=[
                            "Both constraints are physics: a density",
                            "cannot be negative and integrates to one.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Projected gradient",
                        formula="w <- proj_simplex(w - grad / L)",
                        detail=[
                            "Step 1/L with L the Lipschitz constant.",
                            "Stationarity measured scale-free.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Harmonic alternative",
                        detail=[
                            "Unknowns are symmetry-projected",
                            "coefficients to degree L instead.",
                            "Even degrees only: pole figures carry",
                            "no odd information.",
                        ],
                    ),
                ],
            ),
            (
                "4 - judge the fit, do not trust it",
                [
                    AlgorithmStage(
                        label="Residual report",
                        role="output",
                        formula="relative residual, MAE, max error, coverage",
                        detail=[
                            "Errors in m.r.d., so they read directly",
                            "against the texture strength.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Recalculate unfitted poles",
                        role="decision",
                        detail=[
                            "The only check that carries information:",
                            "the map is non-injective, so a low",
                            "residual on the fitted poles proves",
                            "nothing about the ODF.",
                        ],
                    ),
                ],
            ),
        ],
        title="Pole-figure inversion to an orientation distribution",
        subtitle="ODF.invert_pole_figures - the discrete route, with the harmonic route alongside",
        description=(
            "Four-lane flow sheet. Lane 1 takes measured pole figures on the m.r.d. scale "
            "and applies the defocus correction. Lane 2 builds the kernel response operator "
            "and rescales it to the observations' scale. Lane 3 solves the non-negative, "
            "simplex-constrained, Tikhonov-regularised least-squares problem by projected "
            "gradient, with the harmonic series expansion as the alternative unknown. Lane 4 "
            "reports residuals and recalculates poles that were not fitted, which is the only "
            "check that carries information because the forward map is not injective."
        ),
        notes=[
            SideNote(
                stage_index=3,
                title="Failure: a fit that reports success",
                lines=[
                    "Skip the rescale and the weights cannot",
                    "absorb the kernel mean, because they are",
                    "constrained to sum to one. The solver stalls",
                    "near a relative residual of 1, its step stops",
                    "moving, and it reports convergence.",
                ],
            ),
            SideNote(
                stage_index=5,
                title="Constraint: scale-free stationarity",
                lines=[
                    "The raw step is proportional to 1/L, so on a",
                    "large-operator system the first step is tiny",
                    "for that reason alone. Multiplying by L and",
                    "dividing by |A^T b| makes one tolerance mean",
                    "the same thing in any units.",
                ],
            ),
            SideNote(
                stage_index=6,
                title="Structural: the odd part is absent",
                lines=[
                    "Friedel's law makes pole figures",
                    "centrosymmetric, so they determine only the",
                    "even-order ODF. No regularisation recovers",
                    "the odd part; ghost correction infers one",
                    "from positivity and states its cost.",
                ],
            ),
        ],
        footer=[
            "Not modelled here: ghost correction (see the ghost-correction page), grain",
            "statistics, and any claim of uniqueness - several ODFs reproduce one pole-figure",
            "set exactly.",
        ],
    )


def rietveld_refinement_figure() -> str:
    """Flow sheet for `refine_rietveld`."""

    return algorithm_flow_svg(
        [
            (
                "1 - the model, assembled once",
                [
                    AlgorithmStage(
                        label="Measured profile",
                        role="input",
                        formula="y_obs(2theta), background NOT subtracted",
                        detail=[
                            "Subtracting first removes the",
                            "background-scale correlation the",
                            "uncertainties depend on.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Enumerate reflections",
                        formula="over a PADDED angular window",
                        detail=[
                            "Once, from the starting cell. Padding",
                            "stops a dilating cell moving a",
                            "reflection in or out mid-refinement.",
                        ],
                    ),
                ],
            ),
            (
                "2 - the forward calculation, per evaluation",
                [
                    AlgorithmStage(
                        label="Positions",
                        formula="from the dilated cell + zero shift",
                    ),
                    AlgorithmStage(
                        label="Intensities",
                        formula="m |F|^2 L(theta) P_march",
                        detail=[
                            "Multiplicity, structure factor with",
                            "B_iso, Lorentz-polarisation, and the",
                            "March-Dollase texture factor.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Profile + background",
                        formula="Thompson-Cox-Hastings pV + Chebyshev",
                        detail=[
                            "Widths from the Caglioti form",
                            "U tan^2 t + V tan t + W, plus Y.",
                        ],
                    ),
                ],
            ),
            (
                "3 - refine, bounded",
                [
                    AlgorithmStage(
                        label="Trust-region least squares",
                        formula="min sum w (y_obs - y_calc)^2",
                        detail=[
                            "scipy least_squares. Every parameter",
                            "bounded, because an unbounded fit",
                            "reaches a lower R_wp at an impossible",
                            "cell.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Uncertainties",
                        detail=[
                            "From the Jacobian at the solution,",
                            "with the background still in the model.",
                        ],
                    ),
                ],
            ),
            (
                "4 - judge it, and not by R_wp alone",
                [
                    AlgorithmStage(
                        label="R factors",
                        role="output",
                        formula="R_p, R_wp, R_exp, GoF = R_wp/R_exp",
                        detail=["Fractions, not percentages."],
                    ),
                    AlgorithmStage(
                        label="Durbin-Watson",
                        role="decision",
                        detail=[
                            "Serial correlation of the weighted",
                            "residuals. Near 2 uncorrelated; well",
                            "under 1 means systematic misfit that",
                            "R_wp does not show.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Residual curve",
                        role="output",
                        detail=[
                            "The most informative single output.",
                            "Derivative-shaped at peaks: positions.",
                            "Symmetric at peaks: widths.",
                            "High-angle only: B_iso or absorption.",
                        ],
                    ),
                ],
            ),
        ],
        title="Rietveld whole-profile refinement",
        subtitle="refine_rietveld - fitting the pattern point by point, not extracted intensities",
        description=(
            "Four-lane flow sheet. Lane 1 takes the raw measured profile with its background "
            "intact and enumerates reflection families once over a padded angular window. "
            "Lane 2 recomputes positions, intensities, profiles and background at every "
            "evaluation. Lane 3 minimises the weighted residual by bounded trust-region least "
            "squares and takes uncertainties from the Jacobian. Lane 4 reports the R factors, "
            "the Durbin-Watson statistic and the residual curve, which diagnoses what to "
            "refine next in a way no scalar does."
        ),
        notes=[
            SideNote(
                stage_index=5,
                title="Constraint: bounds encode physics",
                lines=[
                    "Zero shift within +/- 1 degree: beyond that",
                    "it is a broken diffractometer, and left free",
                    "it swaps places with the cell. Cell dilation",
                    "within 10 percent: a refinement wanting more",
                    "has misidentified the phase.",
                ],
            ),
            SideNote(
                stage_index=6,
                title="Practice: refine incrementally",
                lines=[
                    "The default set is scale, zero, cell, one",
                    "width, and background. Turning everything on",
                    "at once reaches a low R_wp at a meaningless",
                    "minimum, because zero-against-cell and",
                    "background-against-scale trade freely.",
                ],
            ),
        ],
        footer=[
            "Use this when phase fractions, structure or profile parameters are the target.",
            "For the cell alone, the extrapolation route makes fewer assumptions and gives a",
            "better lattice parameter.",
        ],
    )


def ebsd_grain_metrics_figure() -> str:
    """Flow sheet for segmentation and the local-misorientation family."""

    return algorithm_flow_svg(
        [
            (
                "1 - the map, as a graph",
                [
                    AlgorithmStage(
                        label="Orientations on a grid",
                        role="input",
                        detail=[
                            "Square (4/8) or hexagonal (6)",
                            "topology, with phases per point.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Neighbour pairs",
                        formula="first shell; phase boundaries dropped",
                        detail=[
                            "A misorientation between different",
                            "phases is not defined here.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Pair disorientation",
                        formula="min over the symmetry orbit",
                        detail=[
                            "576 candidates per cubic pair,",
                            "evaluated vectorised over all pairs.",
                        ],
                    ),
                ],
            ),
            (
                "2 - grains",
                [
                    AlgorithmStage(
                        label="Threshold the edges",
                        role="decision",
                        formula="keep pairs with angle <= theta_c",
                        detail=["Conventionally 5 to 15 degrees."],
                    ),
                    AlgorithmStage(
                        label="Connected components",
                        formula="flood fill = components of the edge set",
                        detail=[
                            "Grain ids numbered by each",
                            "component's lowest member, so they are",
                            "a function of the data.",
                        ],
                    ),
                ],
            ),
            (
                "3 - the four local metrics, by what each compares",
                [
                    AlgorithmStage(
                        label="KAM",
                        formula="point vs its neighbours",
                        detail=[
                            "Short-wavelength gradient.",
                            "Exclude boundaries or it becomes a",
                            "boundary map.",
                        ],
                    ),
                    AlgorithmStage(
                        label="GROD",
                        formula="point vs its grain reference",
                        detail=["Long-wavelength rotation."],
                    ),
                    AlgorithmStage(
                        label="GOS / GAM",
                        formula="grain averages of the two",
                        detail=["Per-grain deformation measures."],
                    ),
                ],
            ),
            (
                "4 - dislocation content",
                [
                    AlgorithmStage(
                        label="Curvature tensor",
                        formula="kappa_ij = d omega_i / d x_j",
                    ),
                    AlgorithmStage(
                        label="Nye tensor",
                        formula="alpha = kappa^T - tr(kappa) I",
                    ),
                    AlgorithmStage(
                        label="GND density",
                        role="output",
                        formula="rho = sum |alpha_measurable| / b",
                        detail=[
                            "A LOWER BOUND: a surface scan cannot",
                            "measure the gradient normal to itself,",
                            "and statistically stored dislocations",
                            "produce no curvature at all.",
                        ],
                    ),
                ],
            ),
        ],
        title="Grains, local misorientation, and dislocation density",
        subtitle="segment_grains, the KAM/GROD/GOS/GAM family, and the Nye route to GND",
        description=(
            "Four-lane flow sheet. Lane 1 turns the orientation grid into a neighbour graph "
            "and computes symmetry-reduced pair disorientations. Lane 2 thresholds the edges "
            "and takes connected components, giving grains. Lane 3 derives the four local "
            "misorientation metrics, distinguished by what each compares a point with. Lane 4 "
            "forms the lattice curvature and Nye tensors and reports a lower-bound "
            "geometrically necessary dislocation density."
        ),
        notes=[
            SideNote(
                stage_index=3,
                title="Inherent: flood fill merges gradients",
                lines=[
                    "Points connected by a chain of small steps",
                    "join, so a grain with a continuous gradient",
                    "can exceed theta_c end to end and remain one",
                    "grain. That is the definition, not a defect -",
                    "and it is why GROD exists.",
                ],
            ),
            SideNote(
                stage_index=5,
                title="Constraint: KAM is grid-dependent",
                lines=[
                    "A larger kernel smooths and lowers KAM, so",
                    "values are comparable only at equal step size",
                    "and equal order. Report both.",
                ],
            ),
        ],
        footer=[
            "Every number here depends on a choice: theta_c and connectivity for grains, step",
            "and kernel for KAM, method and Burgers vector for GND. A result quoted without",
            "them is not reproducible.",
        ],
    )


def kikuchi_geometry_figure() -> str:
    """Flow sheet for `simulate_kikuchi_pattern`."""

    return algorithm_flow_svg(
        [
            (
                "1 - which planes can give a band",
                [
                    AlgorithmStage(
                        label="Enumerate planes",
                        role="input",
                        formula="(hkl) to max_index",
                        detail=["Cost is cubic in the index limit."],
                    ),
                    AlgorithmStage(
                        label="Antipodal reduction",
                        detail=[
                            "(hkl) and its opposite are the same",
                            "plane and would draw one band twice.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Centring and intensity",
                        role="decision",
                        formula="reflection condition, then |F|^2",
                        detail=[
                            "A systematically absent reflection",
                            "produces no band at all.",
                        ],
                    ),
                ],
            ),
            (
                "2 - into the laboratory frame",
                [
                    AlgorithmStage(
                        label="Carry the normal through",
                        formula="n_lab = T_spec->lab . R_cry->spec . n_cry",
                        detail=[
                            "Frames are checked, not assumed: a",
                            "mismatch gives a plausible wrong",
                            "pattern.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Bragg angle",
                        formula="sin(theta_B) = lambda / 2d",
                        detail=[
                            "Relativistic wavelength. Planes that",
                            "cannot satisfy it are dropped.",
                        ],
                    ),
                ],
            ),
            (
                "3 - onto the detector",
                [
                    AlgorithmStage(
                        label="Gnomonic projection",
                        formula="great circles -> straight lines",
                        detail=[
                            "Origin at the pattern centre; one unit",
                            "is one detector distance, so the frame",
                            "is detector-size independent.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Band edges",
                        role="output",
                        formula="angular width = 2 theta_B",
                        detail=[
                            "Width measures d INVERSELY: wide bands",
                            "are low-d planes.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Zone axes",
                        role="output",
                        formula="h u + k v + l w = 0",
                        detail=[
                            "The bright intersections a human",
                            "indexes by.",
                        ],
                    ),
                ],
            ),
        ],
        title="Kikuchi band geometry and pattern simulation",
        subtitle="simulate_kikuchi_pattern - the forward model that indexing inverts",
        description=(
            "Three-lane flow sheet. Lane 1 enumerates candidate lattice planes, reduces "
            "antipodal pairs, and filters by the centring reflection condition and a kinematic "
            "intensity threshold. Lane 2 carries each plane normal into the laboratory frame "
            "through the crystal orientation and computes the Bragg angle from the "
            "relativistic electron wavelength. Lane 3 projects gnomonically, where great "
            "circles become straight lines, and reports band centre lines, edges at twice the "
            "Bragg angle, and the zone axes where bands intersect."
        ),
        notes=[
            SideNote(
                stage_index=6,
                title="Constraint: the forward hemisphere only",
                lines=[
                    "Directions travelling away from the detector",
                    "have no intersection and are reported invalid.",
                    "A silent wrap would place a band on the",
                    "opposite side of the pattern.",
                ],
            ),
            SideNote(
                stage_index=7,
                title="Reading: min_d_spacing drops WIDE bands",
                lines=[
                    "Band width grows as spacing falls, so a",
                    "minimum-spacing cut removes the widest,",
                    "weakest, high-order bands - not the",
                    "narrowest.",
                ],
            ),
        ],
        footer=[
            "Geometric and kinematic: band positions and widths are exact, intensities are a",
            "|F|^2 proxy. No excess/deficiency asymmetry, no dynamical contrast, no inelastic",
            "background, no detector response.",
        ],
    )


def phase_identification_figure() -> str:
    """Flow sheet for `identify_phase_from_pattern`."""

    return algorithm_flow_svg(
        [
            (
                "1 - the measurement, once, shared by every candidate",
                [
                    AlgorithmStage(
                        label="Measured diffractogram",
                        role="input",
                        formula="I(2*theta), lambda",
                        detail=[
                            "A raw scan with its background intact,",
                            "and the radiation it was measured with.",
                            "A position without its wavelength is not",
                            "a measurement of anything.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Detect and fit",
                        formula="detect_and_fit_peaks",
                        detail=[
                            "Ricker filter on the variance-stabilised",
                            "profile, then a pseudo-Voigt per peak.",
                            "Yields positions with ESDs and",
                            "integrated intensities.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Candidate structures",
                        role="input",
                        detail=[
                            "CIF files, or catalogue entries.",
                            "Two or more make it a choice; one is",
                            "a check on that phase, and the report",
                            "says so.",
                        ],
                    ),
                ],
            ),
            (
                "2 - per candidate: place its lines, then assign them",
                [
                    AlgorithmStage(
                        label="Refine one cell dilation",
                        formula="min_s sum_p min_j clip(|2th_p - 2th_j(s)|, eps)",
                        detail=[
                            "Grid search on s in [1-d, 1+d], d = 0.02.",
                            "A grid, not a gradient: the objective is",
                            "piecewise linear with a minimum at every",
                            "near-coincidence.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Enumerate reflections",
                        formula="generate_powder_reflections",
                        detail=[
                            "The candidate's own symmetry and",
                            "systematic absences, not a generic",
                            "(hkl) list. Families below the intensity",
                            "floor are never offered for matching.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Global assignment",
                        formula="Hungarian on |2th_obs - 2th_calc|",
                        detail=[
                            "One-to-one over all pairings at once.",
                            "A greedy nearest-line pass can assign two",
                            "peaks to one reflection and strand the",
                            "true partner, invisibly.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Cannot be indexed",
                        role="reject",
                        detail=[
                            "Predicts no line in range at all:",
                            "scored zero with the reason stated,",
                            "never raised. One impossible CIF among",
                            "five must not cost the other four.",
                        ],
                    ),
                ],
            ),
            (
                "3 - four criteria, each bounded, each failing differently",
                [
                    AlgorithmStage(
                        label="Explained intensity",
                        role="decision",
                        formula="E = sum_indexed A_p / sum_all A_p",
                        detail=[
                            "Intensity-weighted, not counted:",
                            "a strong unindexed peak is a second",
                            "phase; a weak one is a trace.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Completeness",
                        role="decision",
                        formula="C = |observed strong| / |predicted strong|",
                        detail=[
                            "Inside the measured span only.",
                            "This is what separates two cells",
                            "differing by a centring: a centring is",
                            "a claim about absent lines.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Position",
                        role="decision",
                        formula="P = max(0, 1 - <|d2th|> / eps)",
                        detail=[
                            "How far inside the window the lines",
                            "landed, not whether they landed inside",
                            "it. Widening eps judges every match",
                            "against the laxer standard.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Intensity agreement",
                        role="decision",
                        formula="S = 1 - (1/2) sum |o_i - c_i|, unit sum",
                        detail=[
                            "One minus Bray-Curtis: scale-free and",
                            "bounded. Undefined below two indexed",
                            "lines, and then renormalised away rather",
                            "than counted as a failure.",
                        ],
                    ),
                ],
            ),
            (
                "4 - rank, then qualify the winner twice",
                [
                    AlgorithmStage(
                        label="Weighted mean",
                        formula="score = sum_D w_k x_k / sum_D w_k",
                        detail=[
                            "Over the criteria that are defined.",
                            "Weights declared and overridable:",
                            "they encode a judgement about evidence,",
                            "not a law of diffraction.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Conclusive?",
                        role="decision",
                        formula="best score >= minimum_score",
                        detail=[
                            "If not: none of the candidates offered",
                            "accounts for this pattern. Widen the",
                            "list, suspect a mixture, or check the",
                            "tolerance against the aberrations.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Decisive?",
                        role="decision",
                        formula="best - runner_up >= decisive_margin",
                        detail=[
                            "If not: this scan does not tell the top",
                            "two apart. The remedy is measurement,",
                            "not computation - high-angle counts,",
                            "another wavelength, or chemistry.",
                        ],
                    ),
                    AlgorithmStage(
                        label="Ranked report",
                        role="output",
                        detail=[
                            "Every candidate with its four criteria,",
                            "its refined dilation, M_N and F_N, its",
                            "unindexed peaks, and describe().",
                        ],
                    ),
                ],
            ),
        ],
        title="Phase identification among candidate structures",
        subtitle=(
            "identify_phase_from_pattern - ranking supplied candidates, "
            "and declining to choose"
        ),
        description=(
            "Four-lane flow sheet. Lane 1 detects and fits the peaks of the measured scan once, "
            "and takes the candidate structures as CIF files or catalogue entries. Lane 2 runs "
            "per candidate: one uniform cell dilation is refined by grid search, the candidate's "
            "own reflections are enumerated under its symmetry, and peaks are assigned to lines "
            "by the Hungarian algorithm; a candidate that predicts no line at all is scored zero "
            "with a stated reason rather than aborting the comparison. Lane 3 scores four "
            "bounded criteria that fail for different physical reasons: explained intensity, "
            "completeness, position agreement and intensity agreement. Lane 4 ranks by their "
            "weighted mean and then qualifies the winner twice - whether it explains the pattern "
            "in absolute terms, and whether it is distinguished from the runner-up."
        ),
        notes=[
            SideNote(
                stage_index=3,
                title="Why a dilation cannot rescue a wrong phase",
                lines=[
                    "d -> s d leaves every ratio d_hkl / d_h'k'l'",
                    "unchanged exactly, and the ratios are what",
                    "indexing tests. A candidate a scale factor",
                    "rescues is the right structure with the wrong",
                    "cell size. The factor is reported, not hidden:",
                    "pinned at the edge means stretched as far as",
                    "allowed and still not fitting.",
                ],
            ),
            SideNote(
                stage_index=10,
                title="Why intensity is weighted least",
                lines=[
                    "Preferred orientation, microabsorption,",
                    "extinction and a coarse powder all move a",
                    "measured intensity by factors while moving no",
                    "peak position at all. Weighting S heavily would",
                    "reject the correct phase of any textured",
                    "specimen, which is most engineering specimens.",
                ],
            ),
            SideNote(
                stage_index=13,
                title="Not retrieval, and not quantification",
                lines=[
                    "The candidates must be supplied: nothing here",
                    "searches a database for a structure nobody",
                    "proposed, so a low best score may mean the",
                    "right phase was never offered. And when several",
                    "candidates each explain part of the pattern the",
                    "next step is a multi-phase Rietveld refinement,",
                    "not a larger score.",
                ],
            ),
        ],
    )


def main() -> int:
    """Write every algorithm figure into ``docs/figures/``."""

    figures = {
        "or_determination_algorithm.svg": or_determination_figure(),
        "variant_correspondence_algorithm.svg": variant_correspondence_figure(),
        "composite_saed_algorithm.svg": composite_saed_figure(),
        "saed_indexing_algorithm.svg": saed_indexing_figure(),
        "pole_figure_inversion_algorithm.svg": pole_figure_inversion_figure(),
        "rietveld_refinement_algorithm.svg": rietveld_refinement_figure(),
        "phase_identification_algorithm.svg": phase_identification_figure(),
        "ebsd_grain_metrics_algorithm.svg": ebsd_grain_metrics_figure(),
        "kikuchi_geometry_algorithm.svg": kikuchi_geometry_figure(),
    }
    FIGURES.mkdir(parents=True, exist_ok=True)
    for name, svg in figures.items():
        path = FIGURES / name
        path.write_text(svg, encoding="utf-8")
        print(f"wrote {path.relative_to(FIGURES.parents[1])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
