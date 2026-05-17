# MTEX Generator Scripts

These MATLAB/MTEX scripts generate external parity artifacts from the shared JSON campaign files in
`fixtures/mtex_parity/campaigns/`.

The target baseline for the new campaign runner is MTEX `6.0.0`, released in November 2024. This is
recent enough for the current validation program without depending on the latest MTEX release line.

Run from the repository root on a MATLAB machine where MTEX has already been started:

```matlab
addpath("scripts/mtex_generators")
run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
                         "fixtures/mtex_parity/results/mtex")
```

Repeat for each campaign file. Bring the resulting `fixtures/mtex_parity/results/mtex/` tree back
to the PyTex machine, then generate PyTex results and compare:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

Pending XRDML cases emit skipped JSON records until the cubic and hexagonal pole-figure fixture
files are added and their case status is changed to `active`.
