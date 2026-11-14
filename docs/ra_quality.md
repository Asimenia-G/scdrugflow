# RA descriptive quality control

This step runs **after** the Gehring-Pachter matched RA contrast.
It does not perform differential-expression testing or make causal claims.

It reopens the original, unchanged `.h5ad`, selects the same fixed-dose
comparison, checks per-batch group sizes, and independently recomputes the
saved gene-wise differences to verify the result table.

It also calculates the fraction of cells with non-zero expression (`X > 0`)
within each of the four batch-by-treatment groups. Sparse and dense matrices
are handled consistently. A gene's minimum detection fraction is **reported,
not used as a new filter**, so this QC does not silently change rankings.

The scatter plot compares descriptive gene differences across batches.
Correlation and directional concordance are **descriptive**, not statistical
replication or significance evidence. The group counts are numbers of cells,
not independent biological samples.

## Run as part of Nextflow

Run the same RA command as before, but use a new results directory, e.g.:

```bash
nextflow run main.nf -profile local \
  --input data/GehringPachter2019.h5ad \
  --mode ra --ra_bmp4 8 --ra_egf_bfgf 200 --ra_mix 0 \
  --ra_dose 10 --ra_batches 0,1 --min_group_cells 20 \
  --outdir results/gehring_ra_qc
```

Additional results will appear at `results/gehring_ra_qc/ra_quality/`:

- `quality_summary.json`: counts and result verification, descriptive concordance
- `gene_quality.tsv`: per-gene detection fractions across all four groups
- `top20_gene_quality.tsv`: detection fractions for the original top 20 genes
- `batch_concordance.png`: scatterplot of effects in batches 0 and 1

This quality-control step does not establish biological replication.
