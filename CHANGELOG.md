# Changelog

## 0.1.0 — first validated portfolio release (November 2026)

- Nextflow DSL2 inspection, exploratory general analysis, and matched-dose RA case study.
- Synthetic input generator, sparse/dense numerical tests and GitHub Actions Python tests.
- Dataset-specific metadata audit of the Gehring–Pachter scPerturb dataset.
- RA 0 vs 10 µM comparison with fixed background doses, within-batch contrasts, and equal batch weighting.
- Independent validation of output tables, expression-detection fractions and descriptive batch-concordance figure.
- Documented biological limits: no independent biological replicates established; no inferential differential expression.
- Fixed the standard-library `inspect` naming collision by renaming the script to `inspect_dataset.py`.
- Bounded floating-point gene-detection fractions to the mathematically valid interval [0, 1].

**User-verified execution on Hypatia (Ubuntu):** 18/18 Python tests passed. Synthetic Nextflow workflow and the public Gehring–Pachter RA + quality-control workflow both completed. Java 21 and Nextflow 26.10.6 were used. The final small bounding change should be rerun before GitHub publication.
