# Scope and possible future work

## Implemented and tested in version 0.1

- AnnData inspection and biological metadata audit
- General exploratory analysis on synthetic counts (QC, PCA and descriptive comparisons)
- Matched RA-dose case study on public scPerturb data
- Separate within-batch expression summaries and equal-weight descriptive average
- Independent re-computation of saved results, detection fractions and batch-concordance plot
- Python numerical and integration tests, plus Python GitHub Actions CI
- Manual end-to-end Nextflow execution on an Ubuntu server

## Not implemented (deliberately out of scope)

- General sample-aware pseudobulk differential expression or inferential statistics
- Arbitrary treatment-combination modelling or causal drug-effect estimation
- Production-level containers, cloud/HPC schedulers or large-cohort scaling
- Automated Nextflow end-to-end CI and per-platform environment lockfiles

These are potential future enhancements, **not** claimed capabilities of this release.
