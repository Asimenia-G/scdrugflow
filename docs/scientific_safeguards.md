# Scientific safeguards

## Matrix provenance

The general `analyse` mode requires users to declare whether the selected matrix represents `counts` or `log1p` values. Integer-like values alone cannot prove raw-count provenance.

For the Gehring–Pachter `ra` mode, `X` contains non-integer values and no raw-count layer; an `expm1(X)` library-size diagnostic strongly suggests log1p-normalised expression. The case study deliberately does **not** normalise again, perform count-level QC or call differential expression.

## Treatment definitions and matching

Drug names in Gehring–Pachter metadata are constant. Treatment conditions are encoded in four dose columns. The RA case study compares RA 0 and 10 µM **while holding three other dose factors fixed**, stratified within batches 0 and 1. RA 0 is not a completely untreated culture.

The general `analyse` mode does **not** automatically match doses or adjust for batch; it must only be used with a scientifically meaningful condition column and reference.

## What is reported

For each gene and batch, `delta = mean(stored log1p X in RA-treated cells) - mean(stored log1p X in RA-zero cells)`. The reported overall value is the arithmetic mean of these two batch-specific deltas (equal batch weighting). It is **not** a log2 fold change or an estimate of drug efficacy.

Gene detection fraction is defined as the proportion of cells with **stored** `X > 0`. It is provided as a QC diagnostic, not used to selectively rerank effects. Roundoff outside [0, 1] is clipped.

## Replication and claims

Cells and batch labels are **not automatically biological replicates**. Independent replication has not been established in this release. The pipeline computes no p-values, FDR, confidence intervals, pseudobulk differential expression or causal estimates. Correlation and sign concordance between two batches are descriptive observations, not inferential replication.

## QC, figures and scope

QC thresholds in the synthetic count demonstration illustrate validation logic; they have not been optimised for all experiments. PCA is exploratory. Any gene annotations or pathway interpretations require separate biological review, and the output is not suitable for clinical or pharmacological decisions.
