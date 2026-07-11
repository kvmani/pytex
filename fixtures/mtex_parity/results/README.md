# MTEX And PyTex Parity Results

This directory is the exchange point for generated parity artifacts.

- `mtex/`: JSON and figure artifacts generated on a MATLAB system with MTEX.
- `pytex/`: JSON and figure artifacts generated locally by PyTex.

The result files are intentionally machine-readable and use `schemas/parity_result.schema.json`.
They are produced from the shared input campaigns in `fixtures/mtex_parity/campaigns/`.

## Directory Layout

- `mtex/<campaign_id>/`: result JSON generated on the MATLAB + MTEX machine
- `pytex/<campaign_id>/`: result JSON generated on the PyTex machine

Do not rename the campaign or case files after generation. The comparator matches results using the
shared campaign ids and case ids.

## Typical Workflow

1. Run the MTEX generator scripts on the MATLAB machine so they populate `mtex/`.
2. Bring that populated `mtex/` tree back into this repository.
3. Generate matching PyTex result JSON into `pytex/`.
4. Compare the two roots field by field.

Typical local comparison after bringing MTEX results back to this machine:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```
