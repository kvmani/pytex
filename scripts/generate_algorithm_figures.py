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


def main() -> int:
    """Write every algorithm figure into ``docs/figures/``."""

    figures = {
        "or_determination_algorithm.svg": or_determination_figure(),
        "variant_correspondence_algorithm.svg": variant_correspondence_figure(),
        "composite_saed_algorithm.svg": composite_saed_figure(),
        "saed_indexing_algorithm.svg": saed_indexing_figure(),
    }
    FIGURES.mkdir(parents=True, exist_ok=True)
    for name, svg in figures.items():
        path = FIGURES / name
        path.write_text(svg, encoding="utf-8")
        print(f"wrote {path.relative_to(FIGURES.parents[1])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
