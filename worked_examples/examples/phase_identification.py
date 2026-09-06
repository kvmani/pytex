"""Worked examples: identifying a phase among several candidate structures.

Every expected value here has independent provenance, and none is a copied
output of this code.

The generating phase of a synthetic pattern is known before the calculation
starts, so "the ranking returns it first" is checkable against the fixture that
produced the pattern rather than against a previous run. The cell dilation of
the second example is *imposed* by the example itself, so recovering it is a
measurement against a value the example set. The d-spacing-ratio identity is
algebra: a uniform dilation multiplies every spacing by the same factor and so
leaves every ratio unchanged, exactly, whatever the structure. And a comparison
must return one row per candidate offered whether or not that candidate could be
indexed, which is a contract the multi-candidate case depends on.

See :doc:`../../theory/phase_identification_from_powder_patterns`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

_THEORY = SeeAlso(
    "Phase identification from powder patterns",
    "../../theory/phase_identification_from_powder_patterns",
)
_ALGORITHM = SeeAlso("Phase identification algorithm", "../../algorithms/phase_identification")
_XRD_THEORY = SeeAlso("Powder XRD and SAED theory", "../../theory/powder_xrd_and_saed")

_THETA = SymbolUse(r"\theta", "Bragg half-angle.")
_D_SPACING = SymbolUse(r"d_{hkl}", "Interplanar spacing of the (hkl) family.")

#: The synthetic-scan preamble every end-to-end example below shares. Written
#: once because the examples differ in the question they ask, not in the
#: specimen they ask it about, and a reader comparing two of them should see
#: only the difference.
_SCAN_SETUP = (
    "import numpy as np\n"
    "from pytex.app.phases import builtin_phase\n"
    "from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern\n"
    "from pytex.diffraction.xrd_measurement import MeasuredPowderPattern\n"
    "from pytex.diffraction.xrd_phase_identification import identify_phase_from_pattern\n"
    "\n"
    "def scan(phase, seed=7):\n"
    "    radiation = RadiationSpec.cu_ka()\n"
    "    pattern = generate_xrd_pattern(\n"
    "        phase,\n"
    "        radiation=radiation,\n"
    "        two_theta_range_deg=(25.0, 140.0),\n"
    "        resolution_deg=0.01,\n"
    "        broadening_fwhm_deg=0.12,\n"
    "        profile='pseudo_voigt',\n"
    "    )\n"
    "    profile = np.asarray(pattern.intensity_grid, dtype=float)\n"
    "    counts = profile / profile.max() * 30000.0 + 150.0\n"
    "    return MeasuredPowderPattern(\n"
    "        name='synthetic specimen',\n"
    "        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),\n"
    "        intensity=np.random.default_rng(seed).poisson(counts).astype(float),\n"
    "        radiation=radiation,\n"
    "        synthetic=True,\n"
    "    )\n"
)

#: Dilating a phase means dilating its cell and the unit cell's copy of it
#: together; `Phase` requires the two to agree.
_DILATE_SETUP = (
    "from dataclasses import replace\n"
    "\n"
    "def dilated(phase, scale, name):\n"
    "    lattice = replace(\n"
    "        phase.lattice,\n"
    "        a=phase.lattice.a * scale,\n"
    "        b=phase.lattice.b * scale,\n"
    "        c=phase.lattice.c * scale,\n"
    "    )\n"
    "    return replace(\n"
    "        phase,\n"
    "        lattice=lattice,\n"
    "        unit_cell=replace(phase.unit_cell, lattice=lattice),\n"
    "        name=name,\n"
    "    )\n"
)


RANKING_RETURNS_THE_GENERATING_PHASE = WorkedExample(
    id="phase-id-ranking-returns-the-generating-phase",
    title="A pattern generated from nickel ranks nickel first",
    domain="phase-identification",
    scenario=(
        "The elementary check on any identification: generate a powder pattern from a known "
        "structure, hand the ranking that structure along with three plausible competitors, and "
        "require that it comes back first. The answer is fixed by the fixture the pattern was "
        "generated from, before the calculation starts. The competitors are not chosen to be "
        "easy: copper is face-centred cubic like nickel and differs only in cell size, ferrite "
        "is cubic of a similar size but body-centred, and halite is face-centred cubic with a "
        "two-species basis. Between them they exercise all three ways a candidate can be wrong "
        "- the wrong cell dimension, the wrong centring, and the wrong basis. The reported "
        "value is the rank of the true phase, which must be 1."
    ),
    setup=_SCAN_SETUP,
    code=(
        "truth = builtin_phase('ni_fcc').to_phase()\n"
        "candidates = {\n"
        "    key: builtin_phase(key).to_phase()\n"
        "    for key in ('ni_fcc', 'cu_fcc', 'fe_bcc', 'nacl')\n"
        "}\n"
        "report, _ = identify_phase_from_pattern(scan(truth), candidates)\n"
        "result = float(\n"
        "    1 + [item.phase_name for item in report].index('ni_fcc')\n"
        ")"
    ),
    expected=1.0,
    unit="",
    tolerance=0.0,
    reference=(
        "The pattern was generated from the pinned ni_fcc fixture, so the identity of the "
        "specimen is known independently of the ranking. A correct ranking places it first."
    ),
    citation=(
        "Hanawalt, Rinn & Frevel, Ind. Eng. Chem. Anal. Ed. 10 (1938) 457, "
        "doi:10.1021/ac50125a001."
    ),
    symbols=(_THETA, _D_SPACING),
    see_also=(_THEORY, _ALGORITHM, _XRD_THEORY),
    result_format="{:.0f}",
)


REFINEMENT_RECOVERS_AN_IMPOSED_DILATION = WorkedExample(
    id="phase-id-refinement-recovers-an-imposed-cell-dilation",
    title="The refined cell dilation recovers a dilation the example imposed",
    domain="phase-identification",
    scenario=(
        "A CIF records the cell of somebody else's specimen. Yours is a solid solution, or at "
        "another temperature, or stressed, and its cell differs by a fraction of a per cent - "
        "which by Delta(2*theta) = 2 e tan(theta) displaces a back-reflection line by far more "
        "than any sensible matching tolerance. The ranking therefore refines one uniform cell "
        "dilation per candidate before indexing, and reports it.\n\n"
        "Here the specimen's cell is dilated by exactly 1.0040 relative to the tabulated nickel "
        "fixture, and the tabulated fixture is offered as the candidate. The refinement must "
        "recover the factor that was imposed. Note what this also demonstrates: the *candidate* "
        "is the undilated tabulated cell, exactly as it would arrive from a CIF."
    ),
    setup=_SCAN_SETUP + _DILATE_SETUP,
    code=(
        "tabulated = builtin_phase('ni_fcc').to_phase()\n"
        "specimen = dilated(tabulated, 1.0040, 'nickel solid solution')\n"
        "report, _ = identify_phase_from_pattern(\n"
        "    scan(specimen), {'tabulated nickel': tabulated}\n"
        ")\n"
        "result = float(report.best.cell_scale)"
    ),
    expected=1.0040,
    unit="",
    tolerance=5.0e-4,
    reference=(
        "The dilation is imposed by this example, so the value to recover is set before the "
        "calculation runs. The tolerance is the resolution of the scale grid searched over the "
        "default two per cent range, which is 1.0e-4, widened to 5.0e-4 to absorb the counting "
        "noise on the fitted peak positions."
    ),
    citation=(
        "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11 "
        "- Delta d / d = -cot(theta) Delta(theta), the relation that makes a small cell error a "
        "large high-angle displacement."
    ),
    symbols=(_THETA, _D_SPACING),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.4f}",
)


DILATION_PRESERVES_SPACING_RATIOS = WorkedExample(
    id="phase-id-uniform-dilation-preserves-spacing-ratios",
    title="A uniform dilation leaves every ratio of d spacings unchanged",
    domain="phase-identification",
    scenario=(
        "This is the algebra that makes the cell-scale refinement safe rather than a way of "
        "flattering any candidate at all. For a lattice scaled uniformly by s, every "
        "interplanar spacing becomes s*d_hkl, so the ratio d_hkl / d_h'k'l' is unchanged "
        "exactly, for every pair, whatever the structure. Those ratios are precisely what "
        "indexing tests: a candidate whose relative line positions are wrong is wrong at every "
        "scale, and one a scale factor can rescue is the right structure with the wrong cell "
        "size.\n\n"
        "The check dilates nickel by 1.05 - far beyond anything the refinement would search - "
        "and reports the largest absolute deviation between the two sets of spacing ratios. "
        "The families are matched by their Miller indices rather than by position in the list, "
        "because a dilated cell moves every line to lower angle and so pulls an extra one into "
        "a fixed angular window; comparing the two lists elementwise would be comparing "
        "different reflections. The deviation must be zero to machine precision."
    ),
    setup=(
        "import numpy as np\n"
        "from pytex.app.phases import builtin_phase\n"
        "from pytex.diffraction.xrd import RadiationSpec, generate_powder_reflections\n"
        + _DILATE_SETUP
    ),
    code=(
        "radiation = RadiationSpec.cu_ka()\n"
        "phase = builtin_phase('ni_fcc').to_phase()\n"
        "\n"
        "def spacings(target):\n"
        "    lines = generate_powder_reflections(\n"
        "        target, radiation=radiation, two_theta_range_deg=(20.0, 150.0)\n"
        "    )\n"
        "    return {\n"
        "        tuple(int(value) for value in item.miller_indices): item.d_spacing_angstrom\n"
        "        for item in lines\n"
        "    }\n"
        "\n"
        "plain = spacings(phase)\n"
        "stretched = spacings(dilated(phase, 1.05, 'stretched'))\n"
        "shared = sorted(set(plain) & set(stretched))\n"
        "reference = shared[0]\n"
        "result = float(\n"
        "    max(\n"
        "        abs(plain[hkl] / plain[reference] - stretched[hkl] / stretched[reference])\n"
        "        for hkl in shared\n"
        "    )\n"
        ")"
    ),
    expected=0.0,
    unit="",
    tolerance=1.0e-12,
    reference=(
        "Algebra: d_hkl -> s d_hkl under a uniform dilation, so every ratio d_hkl / d_h'k'l' is "
        "identically unchanged. The tolerance is floating-point round-off, not a physical "
        "allowance."
    ),
    citation=(
        "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 2 "
        "- the plane-spacing equations, in which the cell edges enter as an overall scale."
    ),
    symbols=(_D_SPACING,),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.1e}",
)


EVERY_CANDIDATE_IS_RANKED = WorkedExample(
    id="phase-id-every-candidate-offered-is-ranked",
    title="A candidate that cannot be indexed is scored, not dropped",
    domain="phase-identification",
    scenario=(
        "A user who opens five CIF files and gets back four rows cannot tell which one was "
        "discarded or why, and a comparison that aborts because one candidate is impossible "
        "wastes the other four. So a candidate whose lines cannot be matched at all - here a "
        "cell shrunk until every d spacing falls below lambda/2, for which Bragg's law has no "
        "solution at any angle - is recorded with a stated reason and a score of zero, and the "
        "ranking of the rest proceeds.\n\n"
        "Five candidates are offered, one of them impossible. The reported value is the number "
        "of rows returned, which must be five."
    ),
    setup=_SCAN_SETUP + _DILATE_SETUP,
    code=(
        "truth = builtin_phase('ni_fcc').to_phase()\n"
        "candidates = {\n"
        "    key: builtin_phase(key).to_phase()\n"
        "    for key in ('ni_fcc', 'cu_fcc', 'fe_bcc', 'nacl')\n"
        "}\n"
        "candidates['impossible'] = dilated(truth, 0.30, 'impossible')\n"
        "report, _ = identify_phase_from_pattern(scan(truth), candidates)\n"
        "rejected = [item for item in report if item.indexing is None]\n"
        "result = float(len(report) if len(rejected) == 1 else -1)"
    ),
    expected=5.0,
    unit="",
    tolerance=0.0,
    reference=(
        "Five candidates are offered and exactly one of them - the cell shrunk by 0.30, whose "
        "largest d spacing is below lambda/2 for Cu K-alpha - can produce no reflection at any "
        "angle. A ranking that returned four rows, or that raised, would fail this check; the "
        "result is set to -1 unless exactly one candidate was rejected, so a run that indexed "
        "the impossible cell would fail too."
    ),
    citation=(
        "Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 3 "
        "- lambda = 2 d sin(theta) has no solution for d < lambda/2."
    ),
    symbols=(_THETA, _D_SPACING),
    see_also=(_THEORY, _ALGORITHM),
    result_format="{:.0f}",
)


GROUP = ExampleGroup(
    slug="phase-identification",
    title="Phase identification from a powder pattern",
    summary=(
        "A pattern generated from a known fixture ranked against three candidates chosen to be "
        "wrong in three different ways; a cell dilation imposed by the example and recovered by "
        "the refinement; the algebraic identity that makes that refinement safe rather than a "
        "way of flattering any candidate; and the contract a comparison of several uploaded "
        "structures depends on, that an impossible candidate is scored rather than dropped."
    ),
    examples=(
        RANKING_RETURNS_THE_GENERATING_PHASE,
        REFINEMENT_RECOVERS_AN_IMPOSED_DILATION,
        DILATION_PRESERVES_SPACING_RATIOS,
        EVERY_CANDIDATE_IS_RANKED,
    ),
)
