# ASAI analysis code

Version 2026.10.26. Authors: Pauli Tikka and Roman A. Zubarev.

This repository contains the analysis and figure-generation code accompanying
the ASAI study. It includes the original numerical algorithms, their local
runtime adapter, Python plotting, optional R plotting and verification programs.

## Get the frozen inputs

Download `ASAI_Zenodo_data_and_code_2026-10-26.zip` from the Zenodo record cited in
the article (https://zenodo.org/records/23662053). Keep the ZIP
outside the repository, then run:

```bash
python prepare_data.py /path/to/ASAI_Zenodo_data_and_code_2026-10-26.zip
```

On Windows, a quoted path such as
`python prepare_data.py "C:/Users/yourname/Downloads/ASAI_Zenodo_data_and_code_2026-10-26.zip"`
can be used. This copies the frozen data and numerical reference files into
their expected local folders, preserves the repository's code, and verifies
the residue dataset checksum.

The frozen dataset contains 1,500,829 residue observations from 6,342 chains.
The selected-chain annotation input contains 6,908 chains. The data are held
in Zenodo; the .gitignore file excludes restored data and generated results
from normal Git commits.

## Run

Install Python 3.12 or a compatible environment and Node.js. From this folder:

```bash
python -m pip install -r requirements.txt
python run.py --full
python verify.py
python verify_release.py
```

If Node is not on PATH, set `ASAI_NODE` to its executable. The complete Python/
Node analysis was tested with Python 3.12.14, Node.js v24.19.0 and the pinned
Python dependencies. `run.py` without `--full` skips the additional whole-chain
randomisations and intrinsic-index dimensionality/CCA rerun.

To regenerate just the figures and their source tables after data restoration:

```bash
python plot_figures.py
```

Figure files are written to `outputs/figures/` as PDF, SVG and 600-dpi PNG.
Source tables are written to `source_data/`.

## Optional R plotting

```r
install.packages('ggplot2')
source('plot_figures.R')
plot_asai_figures('.')
```

The R renderer uses the same exported source CSVs. It was not executed in the
packaging environment because R was unavailable. Python-generated figures
were visually checked. The upstream R script is retained as additional source
material and is not required by the frozen-data workflow.

## Citation and reuse

Cite the versioned Zenodo record linked from this repository and the associated
article. The release tag for this source bundle is `v2026.10.26`.

Select a licence for the authors' code before the public repository release.
This preparation step does not assign a new licence. Source databases retain
their respective attribution and reuse terms.

## Figures based on the current manuscript tables

After restoring the data release, run `python plot_final12_figures.py` to write
the four review figures to `manuscript_figures/`. These use the unchanged printed
Table 1, Supplementary Table S1 and Supplementary Table S2 values as appropriate.
The residue histogram and SD inputs come from the frozen data release.

The three PDF attachments selected by the authors and their recovered original
R source have separate provenance in `selected_figure_sources/`. Their earlier
input values and figure numbering differ from Final12; see that folder's README.
The optional R renderers have not been executed here. The Python figure route
was executed and checked. Packaging version 2026.10.26 adds the table extracts,
figure-source mapping and manuscript-based plotting route; the statistical
analysis algorithms and frozen residue observations are unchanged.
