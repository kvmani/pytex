# MTEX And PyTex Parity Results

This directory is the exchange point for generated parity artifacts.

- `mtex/`: JSON and figure artifacts generated on a MATLAB system with MTEX.
- `pytex/`: JSON and figure artifacts generated locally by PyTex.

The result files are intentionally machine-readable and use `schemas/parity_result.schema.json`.
They are produced from the shared input campaigns in `fixtures/mtex_parity/campaigns/`.

Typical local comparison after bringing MTEX results back to this machine:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

