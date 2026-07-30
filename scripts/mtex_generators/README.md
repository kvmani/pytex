# MTEX Generator Scripts

These MATLAB/MTEX scripts generate external parity artifacts from the shared JSON campaign files in
`fixtures/mtex_parity/campaigns/`.

The target baseline for the current runner is MTEX `6.0.0`, released in November 2024.

## External MTEX Workflow

### 1. Open MATLAB At The Repository Root

```matlab
cd("C:/path/to/pytex")
```

### 2. Ensure MTEX Is Available In The Session

Start MTEX using your normal local startup flow before running the parity commands.

### 3. Add The Generator Directory

```matlab
addpath("scripts/mtex_generators")
```

### 4. Run One Campaign

```matlab
run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
                         "fixtures/mtex_parity/results/mtex")
```

### 5. Run The Full Current Campaign Set

```matlab
campaigns = {
    "fixtures/mtex_parity/campaigns/orientation_core_cases.json"
    "fixtures/mtex_parity/campaigns/ipf_color_cases.json"
    "fixtures/mtex_parity/campaigns/miller_geometry_cases.json"
    "fixtures/mtex_parity/campaigns/odf_discrete_cases.json"
    "fixtures/mtex_parity/campaigns/or_transformation_cases.json"
    "fixtures/mtex_parity/campaigns/xrdml_pole_figure_cases.json"
    "fixtures/mtex_parity/campaigns/xrdml_odf_reconstruction_cases.json"
};

for k = 1:numel(campaigns)
    run_mtex_parity_campaign(campaigns{k}, "fixtures/mtex_parity/results/mtex");
end
```

### 6. Bring The Results Back To The PyTex Machine

Copy the populated `fixtures/mtex_parity/results/mtex/` tree back into this repository.

### 7. Generate PyTex Result JSON

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
```

### 8. Compare MTEX And PyTex Results

```powershell
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

## Notes

- **`mtex_parity_transformation.m` has never been run against a real MTEX installation.** The
  `or_transformation_v1` campaign (OR-as-misorientation, variant count, and the
  `calcParent2Child` fit) was authored on a machine with no MATLAB/MTEX available. Its PyTex
  side is generated and verified — the Kurdjumov-Sachs and Nishiyama-Wassermann representatives
  come out at the literature 42.85 deg and 45.99 deg, and the fit recovers Greninger-Troiano
  exactly from a Kurdjumov-Sachs nominal — but the MATLAB handler encodes only the *intended*
  comparison and should be expected to need small API corrections on first run.

  Until it has run cleanly once, **a mismatch is not evidence of a PyTex/MTEX disagreement.**
  Please record any corrections in this file so the next person does not rediscover them. Points
  most likely to need attention: the exact spelling of the named-relationship constructors
  (`orientation.GreningerTrojano` in particular), whether `variants(p2c)` is indexed in the same
  order as PyTex's `generate_variants()`, and the argument form accepted by `calcParent2Child`
  in the installed MTEX version.

- No PyTex document may claim MTEX parity for the orientation-relationship stack until this
  campaign's MTEX-side results exist and `scripts/compare_parity_results.py` has been run over
  them. The validation ledger states this limitation explicitly.

- Pending XRDML cases emit skipped JSON records until the cubic and hexagonal pole-figure fixture
  files are added and their case status is changed to `active`.
- The output layout must remain `fixtures/mtex_parity/results/mtex/<campaign_id>/<case_id>.json`
  so the comparator can resolve the shared campaign and case ids.
- The public Sphinx-facing version of these instructions is in
  `docs/site/validation/mtex_regeneration.md`.
