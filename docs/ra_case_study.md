# Gehring–Pachter retinoic acid case study

The public `GehringPachter2019.h5ad` contains 20,382 mouse neural stem cells
with four perturbation factors. The drug-name columns are constant: **dose
columns**, not drug names, specify experimental treatments. The data have 60
unique four-factor dose combinations.

## Input assumptions

- Public scPerturb dataset, checksum-verified against the Zenodo release.
- `X` is non-integer, has no count layer, and the median over sampled cells of
  `sum(expm1(X))` is about 9,995. This strongly suggests prior library-size
  normalisation to 10,000 followed by `log1p`, but is not proof of the exact
  upstream preprocessing protocol. The source should be cited.
- We do not run raw-count filtering, renormalisation or differential-expression
  tests on this dataset.

## Selected comparison

| Batch | RA 0 µM | RA 10 µM |
| --- | ---: | ---: |
| 0 | 661 | 510 |
| 1 | 591 | 537 |

Hold BMP4 = 8 ng/mL, EGF/bFGF = 200 ng/mL, and the
Scriptaid/decitabine mixture = 0 µM. RA = 0 is the *comparison reference*,
not a completely untreated culture; EGF/bFGF and BMP4 are present.

Batch-specific difference for gene *g*:

`delta_g,b = mean(X_g | RA10, batch b) - mean(X_g | RA0, batch b)`

Overall exploratory contrast is the unweighted average across the two
batch-specific differences. This prevents unequal cell counts across batches
from directly determining the pooled contrast. The expression differences
are **not** log2 fold changes. They are not statistical effect estimates
with independent biological replicates.

**Batch IDs have not been verified as independent biological replicates.**
There are no p-values, adjusted p-values, or causality claims. Sample size
is the number of **cells**, not the number of independent experiments.

## Run on Hypatia

```bash
source /data/scdrugflow_setup/activate_paths.sh
conda activate /data/scdrugflow_setup/envs/scdrugflow
cd /data/scdrugflow

python -m pytest -q

nextflow run main.nf -profile local \
  --input data/GehringPachter2019.h5ad \
  --mode ra \
  --ra_bmp4 8 \
  --ra_egf_bfgf 200 \
  --ra_mix 0 \
  --ra_dose 10 \
  --ra_batches 0,1 \
  --min_group_cells 20 \
  --outdir results/gehring_ra
```

The outputs are in `results/gehring_ra/ra_analysis/`:
`contrast_summary.json`, `group_counts.tsv`, `per_batch_gene_differences.tsv`,
`gene_differences.tsv`, `top_genes.png`, and `README.md`.

This pipeline is not an nf-core pipeline; it is a reproducible teaching-sized
case study that demonstrates a defensible workflow choice.
