# Gehring–Pachter metadata audit

**Purpose:** find comparisons that differ in retinoic-acid dose, while holding
BMP4, EGF/bFGF, Scriptaid/decitabine doses *and* the batch label fixed.

Run from the scDrugFlow repository root, with its Conda environment activated:

```bash
python scripts/audit_gehring.py \
  --input data/GehringPachter2019.h5ad \
  --outdir results/gehring_audit
```

Outputs (all under `/data/scdrugflow/results/gehring_audit` on Hypatia):

- `dose_conditions.tsv`: cell counts per **joint four-factor** dose combination.
- `dose_conditions_by_batch.tsv`: counts for each combination within each batch.
- `ra_matched_pairs.tsv`: RA dose contrasts that occur within the same batch
  and have identical other three dose factors. These are *candidate* matched
  comparisons, not evidence of independent replicates or valid causal inference.
- `matrix_diagnostics.json`: basic expression-scale clues from the first 200
  cells. A constant `sum(expm1(X))` is suggestive of prior log1p library-size
  normalisation, **not proof**; confirm the original processing source.

This diagnostic does not modify or normalise the input file and does not
compute differential-expression statistics. In this dataset no all-zero-dose
condition exists because EGF/bFGF concentration is always positive. An
RA-zero group is a reference for a specifically matched RA contrast, not a
fully untreated control.

The `batch` column cannot be assumed to represent biological replication
without original experimental methods.
