# v2026.10.26-r3

- Publication base: Roman v9-derived r2, with the same frozen corrected residue data and organism annotations.
- Use one interval policy for S1 and Figure 3. Correct the 0.1-boundary floating-point artefacts and the unused residue-bucket diagnostic's fractional-RSA scale. Record each of 28 boundary assignments in `outputs/interval_boundary_adjustments.csv`.
- Export full-precision numerical sources and apply one final display rounding. Table 1 mean RSA is A 0.2009, R 0.3395 and Y 0.1687; printed Final13 uses 0.2008, 0.3396 and 0.1686. Do not double-round to force the printed values.
- Freeze all twenty train/test permutations, including within-half order. The five Table 3 model summaries reproduce Roman's run (0.1154, 0.1371, 0.1176, 0.1168, 0.1358); Final13's length-only value 0.1164 is replaced by 0.1154.
- Correct the length permutation test to use mean absolute r in both observed and null calculations. All six marginal P values remain 0.002.
- Use full-precision recomputed PCA coordinates for AAindex enrichment and downstream PC3/intrinsic-AA checks. Use canonical chain means for the whole-chain null test.
- Parse quoted CSV fields correctly and export classification rules and annotation decisions. Current kingdom labels remain unchanged.
- Use a common eukaryote-minus-prokaryote direction in S2/S3. S3 all-chain PC1 t is +10.14 and PC2 t is −5.86.
- Rebuild Figures 1–4 from the same full-precision outputs as the tables; export PDF/SVG/PNG and numerical source data.
- Add arithmetic verification and independent integer-ASA, SVD, predictor and 2,000-fit checks.
- Separate data and code deposits. Saved results are restored as reference material, not live calculation inputs. No old alternative table generator is included in the production route.

All numeric differences from printed Final13, including pre-existing PCA-coordinate differences, are listed cell by cell in `outputs/differences_from_printed_Final13.csv`. Table 2's calculation structure is retained; its charge term is labelled "Signed charge balance" and the logarithm is explicitly log10.
