# Numerical conventions, v2026.10.26-r3

## RSA and principal components

The frozen corrected residue file is unchanged from the Roman v9/r2 source (SHA-256 in `release.json`). Chain/type mean RSA is the sequential sum divided by the number of residues of that type. Overall Table 1 mean RSA gives equal weight to chains containing the type, not to individual residues. Percentiles 1–99 of those chain means form the 20 × 99 matrix; columns are centred, and PCA coordinates are unit eigenvectors (equivalently, left singular vectors), not singular-value-scaled scores. The joint kingdom PCA uses the 40 eukaryotic/prokaryotic type distributions in one common space. PC signs are set consistently. Full-precision coordinates drive tables, figures and downstream comparisons.

Every displayed number is rounded once from the full-precision result. Printed Word values are reference evidence, never target values for fitting or forced rounding. Independent compensated aggregation and NumPy SVD verify the production JavaScript PCA.

## Supplementary Table S1 and Figure 3

Exactly zero has its own interval. Positive intervals are (0, 0.1], (0.1, 0.2], …, (0.8, 0.9], followed by >0.9. Only for interval assignment, values within an absolute RSA tolerance of 10^-12 of a positive 0.1 boundary are treated as exactly on that boundary. Unmodified chain means are retained for PCA, regressions and other numerical analyses. This resolves binary floating-point boundary artefacts without decimal-rounding the analytical dataset.

For each amino-acid type, the fraction is the number of chains in the interval divided by the number containing that type. CV is the sample SD (ddof = 1) of the twenty fractions divided by their mean. The pooled percentage is 100 times the interval's chain/type count divided by 120,822. Figure 3 displays 100 × CV; the source CSV also retains raw CV. Zero, lower endpoints and upper endpoints follow the same convention in the table and figure. An independent integer-ASA reconstruction validates all interval counts without tolerance-based classification.

## Context models and frozen splits

Five parameter groups comprise log10 resolved chain length, 20 resolved-residue composition fractions, net charge, signed charge balance, and three sequence-encoded hydrophobicity/complexity descriptors. The eligible resolved chain length is at least 40 and the original annotated sequence must be available. Charge, signed charge balance, composition and length use the processed resolved residues. Composition uses a 0.5-residue pseudocount per amino-acid type. Hydrophobicity and complexity use the annotated sequence, including X as zero hydropathy and a distinct symbol for block complexity. The formulas are implemented in `scripts/admissible_parameters.js` and independently checked in `verify_independent.py`.

All predictor sets for a given amino-acid type and replicate share the same train/test split. The twenty original Roman LCG/Fisher–Yates permutations are archived as explicit ranks in `design/train_test_splits.csv`. The first floor(N/2) ranks train the model, and the remaining ranks test it. Rank order is also preserved within each half. The seed alone is insufficient for reproducing these splits in a different RNG implementation.

Regression retains the original ridge regularisation: 10^-9 times the mean diagonal of XᵀX is added to every diagonal entry, including the intercept. Held-out R² is 1 − test SSE / test SST, where test SST is centred on the test-set response mean. Table 3 averages R² over 20 types × 20 splits (400 values per model); unexplained percentage is 100 × (1 − mean R²). All 2,000 fit diagnostics are exported. Within-split paired t statistics compare twenty type-specific model scores; their mean across splits is a descriptive summary, not a single independent-sample significance test.

The length permutation statistic is the **mean of the absolute Pearson correlations across the twenty types**, both in the observed data and in every shuffled dataset. Five hundred permutations use P = (b + 1) / 501. This corrects the original observed/null statistic mismatch; all six reported marginal P values remain 0.002.

## Kingdom comparisons and supplementary tables

All differences and t statistics use eukaryotes minus prokaryotes. S2 deltas are calculated before display rounding; subtracting two printed rounded coordinates may differ in the last digit. S3 t statistics use the twenty paired type coordinates, not chains as independent replicates. Restricted windows are 150 <= L < 350 and 200 <= L < 300. These are common length restrictions, not distribution matching.

## Other article analyses

AAindex rank enrichment now uses recomputed, full-precision RSA coordinates rather than hard-coded printed coordinates. It retains the original protocol of ranking absolute correlations after rounding to three decimals; ties retain AAindex input order. The original normal approximation for the rank-sum statistic is retained and explicitly identified in the JSON output. This protocol choice is distinct from the data/PCA precision correction.

The 51-property selection audit, intrinsic-AA PCA/parallel analysis, CCA closed testing, family sensitivity, PC3 stability/property checks and whole-chain null experiments retain their documented seeds and permutation counts. Their RSA coordinates now come from the same canonical chain means. Tiny floating-point changes may alter a few permutation ties; regenerated exact P values are archived rather than forced to their earlier values. `intrinsic/AAindex_553_audit.csv`, regenerated from `inputs/aaindex_records.json`, documents inclusion decisions. Structural provenance and the limitations of recovering upstream DSSP processing are in `provenance/STRUCTURAL_PROVENANCE.md`.
