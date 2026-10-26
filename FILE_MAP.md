# Source and output map

| Article item | Calculation | Numerical output |
|---|---|---|
| Table 1 | pool.js, decompose.js | source_data/Table1_full_precision.csv |
| Table 2 | admissible_parameters.js | outputs/tables/Table_2.csv |
| Table 3 | admissible_parameters.js, frozen design | source_data/Table3_full_precision.csv; Table3_scores_by_split_type.csv; Table3_fit_diagnostics.csv |
| Table S1 | decompose.js, numeric_policy.cjs | source_data/TableS1_full_precision.csv; TableS1_fractions_by_type.csv |
| Table S2 | decompose.js | source_data/TableS2_full_precision.csv |
| Table S3 | s6_windows.js | source_data/TableS3_full_precision.csv |
| Figure 1 | plot_figures.py, residue CSV | source_data/Figure1_histogram_bins.csv |
| Figure 2 | Table 1 PCA and chain means | source_data/Figure2_PC2_dispersion.csv; Figure2_fit_line.csv |
| Figure 3 | Same S1 calculation | source_data/Figure3_divergence.csv; Figure3_interval_fractions_by_type.csv |
| Figure 4 | Same joint PCA as S2 | source_data/Figure4_kingdom_coordinates.csv |
| Organisms and classification | classify.js | outputs/spec.json; chain_selection_organism_classification.csv |
| PC3 checks | pc3_stability.js, rsa_pc3/*.py | outputs/pc3_stability.json; rsa_pc3/validation_results.json |
| Length/context analysis | parameter_correlations.js, length_slope.js | outputs/parameter_correlations.json; length_slope.json |
| Kingdom controls | kingdom_composition.js, kingdom_joint_control.js, kingdom_null/test.py | corresponding outputs/*.json; kingdom_null/results.json |
| Intrinsic-AA selection/PCA/CCA | intrinsic/*.py | intrinsic/AAindex_553_audit.csv; results.json; family_sensitivity.json |
| AAindex rank enrichment | aaindex_ranks.js | outputs/aaindex_ranks.json |

The six formatted CSV tables are under `outputs/tables/`. Reproduction of the table-only route requires Python's standard library and Node; every other analysis and figure uses pinned Python requirements. Full results live under `numerical_results/` in the separate data ZIP. On restoration, that snapshot is placed under `reference/published_r3/`.

No separate supplementary figure is designated by this production release. Additional PC3 and intrinsic-AA diagnostic calculations are archived numerically; they are not silently assigned manuscript supplementary-figure numbers.
