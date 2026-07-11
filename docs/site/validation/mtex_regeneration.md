# MTEX Parity Regeneration

This page describes the full workflow for generating MTEX parity result artifacts on a separate
MATLAB + MTEX machine, then bringing those artifacts back into the PyTex repository for comparison.

The current target baseline is MTEX `6.0.0`, released in November 2024. That version is recent
enough for the current validation program without depending on the newest MTEX release line.

## What This Produces

The MTEX runner reads the shared campaign JSON files under:

- `fixtures/mtex_parity/campaigns/`

It writes one machine-readable result JSON file per case under:

- `fixtures/mtex_parity/results/mtex/<campaign_id>/`

Each output file uses the shared result schema:

- `schemas/parity_result.schema.json`

## Before You Start

On the MATLAB machine, make sure all of the following are true:

1. the PyTex repository is available locally
2. MATLAB can open that repository directory
3. MTEX is installed
4. MTEX is initialized in the MATLAB session
5. you are running from the repository root, not from some unrelated working directory

If MTEX is not yet initialized in the current MATLAB session, start it first using your normal
local MTEX startup workflow.

## Step By Step

### 1. Open MATLAB At The Repository Root

Change MATLAB into the repository root:

```matlab
cd("C:/path/to/pytex")
```

Use the actual repository path on the MTEX machine.

### 2. Add The MTEX Generator Scripts

```matlab
addpath("scripts/mtex_generators")
```

This exposes the parity runner and helper functions used by the campaign files.

### 3. Run One Campaign

Example for the orientation campaign:

```matlab
run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
                         "fixtures/mtex_parity/results/mtex")
```

That command reads the shared input JSON and writes result JSON files into:

```text
fixtures/mtex_parity/results/mtex/orientation_core_v1/
```

### 4. Run All Active Campaigns

Run the same command once per campaign file:

```matlab
campaigns = {
    "fixtures/mtex_parity/campaigns/orientation_core_cases.json"
    "fixtures/mtex_parity/campaigns/ipf_color_cases.json"
    "fixtures/mtex_parity/campaigns/miller_geometry_cases.json"
    "fixtures/mtex_parity/campaigns/odf_discrete_cases.json"
    "fixtures/mtex_parity/campaigns/xrdml_pole_figure_cases.json"
    "fixtures/mtex_parity/campaigns/xrdml_odf_reconstruction_cases.json"
};

for k = 1:numel(campaigns)
    run_mtex_parity_campaign(campaigns{k}, "fixtures/mtex_parity/results/mtex");
end
```

The XRDML campaigns are currently marked `pending`. The runner still writes skipped result JSON for
them so the artifact set remains provenance-complete.

### 5. Inspect The Output Tree

After the run, you should see folders like:

```text
fixtures/mtex_parity/results/mtex/orientation_core_v1/
fixtures/mtex_parity/results/mtex/ipf_color_v1/
fixtures/mtex_parity/results/mtex/miller_geometry_v1/
fixtures/mtex_parity/results/mtex/odf_discrete_v1/
fixtures/mtex_parity/results/mtex/xrdml_pole_figure_v1/
fixtures/mtex_parity/results/mtex/xrdml_odf_reconstruction_v1/
```

Each folder should contain one `*.json` file per case.

### 6. Bring The MTEX Results Back To The PyTex Machine

Copy the populated `fixtures/mtex_parity/results/mtex/` tree back into this repository on the
PyTex machine.

Do not rename the campaign directories or individual case JSON files. The comparator expects the
shared `campaign_id` and `case_id` layout.

### 7. Generate Matching PyTex Result JSON

On the PyTex machine, from the repository root:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
```

That writes PyTex result files into:

```text
fixtures/mtex_parity/results/pytex/
```

### 8. Compare MTEX And PyTex Results

```powershell
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

If everything matches within the configured tolerances, the script prints:

```text
Parity result comparison passed.
```

If not, it prints field-by-field discrepancies using the shared case ids and field paths.

## How Pending Cases Behave

Some campaigns are deliberately present before their real scientific fixtures are available:

- `xrdml_pole_figure_cases.json`
- `xrdml_odf_reconstruction_cases.json`

Those cases currently have `status: "pending"`.

Behavior:

- MTEX writes skipped result JSON with the pending reason
- PyTex writes skipped result JSON with the pending reason
- the comparator accepts those as matching metadata records

Once real cubic and hexagonal XRDML fixtures are available, update those case files to:

1. replace the placeholder input paths
2. change `status` from `pending` to `active`
3. rerun MTEX and PyTex generation

## Troubleshooting

If MATLAB says a parity helper function is missing:

- verify you ran `addpath("scripts/mtex_generators")`
- verify your current working directory is the repository root

If the runner writes only skipped outputs:

- check the campaign JSON file and confirm the relevant cases are `active`

If the comparator reports missing files:

- confirm the copied MTEX result directory still has the expected `<campaign_id>/<case_id>.json`
  layout

If the comparator reports numeric mismatches:

- first confirm the MTEX machine is using the intended version line
- then confirm the campaign input JSON has not been modified on only one side
- then inspect the specific result field named in the comparator output

## Related Material

- {doc}`mtex_parity_matrix`
- Repository README for the MATLAB runner: `scripts/mtex_generators/README.md`
- Repository README for exchanged result artifacts: `fixtures/mtex_parity/results/README.md`
