# ASAI: reproducible article analysis

Release **v2026.10.26-r3**, derived from the author-supplied Roman v9 audit and the r2 run. This is the single production code route for Tables 1–3, Supplementary Tables S1–S3 and Figures 1–4. Original manuscript files and older alternative R renderers are not required.

## Files to obtain

1. This code archive, also supplied as the GitHub package.
2. `ASAI_Zenodo_data_v2026.10.26-r3.zip`, deposited separately. It contains frozen inputs, train/test ranks, annotations, classification audit, reference tables, numerical outputs and figure source data. It contains no executable analysis code or rendered manuscript figures.

The GitHub code ZIP and the Zenodo code ZIP contain exactly the same bytes. Use either one, not both. Copy the code folder contents to the root of the existing GitHub repository; retain the repository's existing `LICENSE`. See `LICENSING.md`.

## Windows: reproduce the six tables

Install Python 3.12 and Node.js 24. Reopen the terminal after installation. Open a terminal in this code folder, then run (replace the data ZIP path):

```powershell
py -3.12 prepare_data.py "C:\Downloads\ASAI_Zenodo_data_v2026.10.26-r3.zip"
py -3.12 run_tables.py
```

No Python packages need to be installed for the table-only run. The six UTF-8 CSV tables are written to `outputs/tables/`. Full-precision values are in `source_data/`. `outputs/table_verification.json` records the arithmetic checks. Do not run older table generators alongside these commands.

If Node is not found, set its full path before running, for example:

```powershell
$env:ASAI_NODE = "C:\Program Files\nodejs\node.exe"
```

## Windows: reproduce every analysis and all four figures

After restoring the data as above:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run.py --full
```

The full run includes 20,000-permutation intrinsic-AA CCA analyses, parallel analysis, PC3 checks and 200,000 whole-chain allocations. Figures are exported to `outputs/figures/` as PDF, SVG and 600-dpi PNG, plus `Figures_1-4.pdf`. Figure 2 is `Figure2_PC2_dispersion.pdf`. No browser or interactive HTML is required for PDF export.

On macOS/Linux use `python3` instead of `py -3.12`, and `.venv/bin/python` instead of the Windows virtual-environment path. Node 24 must be on PATH. The release was validated with Python 3.12.14 and Node 24.22.1, using the exact dependency versions in `requirements.txt` and one BLAS thread.

## Run scope and references

`run_tables.py` reproduces only the six article tables. `run.py` also performs the remaining JavaScript analyses, PC3 checks and figures. **Use `run.py --full` to reproduce every packaged analysis.** `--no-figures` skips rendering.

`prepare_data.py` checks all archive hashes and pinned scientific input hashes. It restores previously calculated outputs only to `reference/published_r3/`, never to the live calculation folders. The production pipeline regenerates its inputs from the frozen residue CSV and annotations. Printed Final13 tables under `reference/printed_Final13/` are used only for a difference report, after the new values have been calculated.

The archived residue CSV has 1,500,829 residues in 6,342 chains and 120,822 chain–amino-acid-type pairs. The selected-chain list has 6,908 entries; 6,287 chains meet the context-model eligibility rules. The kingdom comparison uses 1,778 eukaryotic and 2,877 prokaryotic chains. Different counts describe different subsets.

See `NUMERICAL_METHODS.md` for calculation conventions, `CHANGELOG.md` for r3 corrections, and `FILE_MAP.md` for output locations. Offline reproduction starts from the processed data; complete original DSSP downloads and the exact historical upstream extraction environment are not recovered by this release.
