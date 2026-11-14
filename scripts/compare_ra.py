"""Describe a matched retinoic-acid contrast in GehringPachter2019.

Input X is previously log1p-library-normalised expression. No new normalisation,
raw-count QC, cell-level differential-expression testing, or p-values.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

DOSE_COLS = ('dose_value', 'dose_value_2', 'dose_value_3', 'dose_value_4')


def select_cells(obs, bmp4, egf_bfgf, mix, ra, batches, min_cells):
    """Select only exact-dose RA=0 and RA=ra cells, stratified by batch."""
    needed = list(DOSE_COLS) + ['batch']
    missing = [column for column in needed if column not in obs.columns]
    if missing:
        raise ValueError(f'Missing metadata columns: {missing}')
    if obs[needed].isna().any().any():
        raise ValueError('Missing dose or batch annotations; cannot match groups')
    if ra <= 0 or min_cells < 1 or not batches or len(set(batches)) != len(batches):
        raise ValueError('Use positive treated RA dose, positive cell threshold, and unique batch IDs')

    numeric = {column: pd.to_numeric(obs[column], errors='raise').to_numpy()
               for column in DOSE_COLS}
    match = (np.isclose(numeric['dose_value'], bmp4, atol=1e-7, rtol=0) &
             np.isclose(numeric['dose_value_2'], egf_bfgf, atol=1e-7, rtol=0) &
             np.isclose(numeric['dose_value_3'], mix, atol=1e-7, rtol=0))
    reference = np.isclose(numeric['dose_value_4'], 0, atol=1e-7, rtol=0)
    treated = np.isclose(numeric['dose_value_4'], ra, atol=1e-7, rtol=0)
    batch_labels = obs['batch'].astype(str).to_numpy()

    rows = []
    for batch in batches:
        in_batch = batch_labels == batch
        n_ref = int(np.sum(match & reference & in_batch))
        n_treated = int(np.sum(match & treated & in_batch))
        if min(n_ref, n_treated) < min_cells:
            raise ValueError(f'Batch {batch}: RA 0 has {n_ref} cells and RA {ra:g} has '
                             f'{n_treated}; require at least {min_cells} each')
        rows.append({'batch': batch, 'n_ra_0': n_ref, 'n_ra_treated': n_treated})

    keep = match & (reference | treated) & np.isin(batch_labels, batches)
    return keep, pd.DataFrame(rows)


def means_and_differences(X, obs, batches, ra):
    """Batch-specific differences of mean(log1p(expression)), not log2FC."""
    dose = pd.to_numeric(obs['dose_value_4']).to_numpy()
    batch_labels = obs['batch'].astype(str).to_numpy()
    effects = []
    ref_means = []
    treated_means = []
    for batch in batches:
        b = batch_labels == batch
        ref = X[b & np.isclose(dose, 0, atol=1e-7, rtol=0)]
        treated = X[b & np.isclose(dose, ra, atol=1e-7, rtol=0)]
        if ref.shape[0] == 0 or treated.shape[0] == 0:
            raise ValueError(f'No comparable cells in batch {batch}')
        r = np.asarray(ref.mean(axis=0)).ravel()
        t = np.asarray(treated.mean(axis=0)).ravel()
        ref_means.append(r)
        treated_means.append(t)
        effects.append(t-r)
    return np.stack(ref_means), np.stack(treated_means), np.stack(effects)


def check_expression_scale(X):
    """Guard against obvious wrong inputs; matching 10k is suggestive, not proof."""
    # Validate every stored value, then sample rows for the scale diagnostic.
    all_values = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    if (not np.isfinite(all_values).all()) or (all_values < 0).any():
        raise ValueError('X contains negative or non-finite entries')
    x = X[:min(X.shape[0], 100), :]
    if sparse.issparse(x):
        x = x.tocsr(copy=True)
    else:
        x = np.asarray(x)
    if sparse.issparse(x):
        x.data = np.expm1(x.data.astype(np.float64))
        totals = np.asarray(x.sum(axis=1)).ravel()
    else:
        totals = np.expm1(x.astype(np.float64)).sum(axis=1)
    median = float(np.median(totals))
    if not 9500 <= median <= 10500:
        raise ValueError('X does not resemble log1p(expression normalised to 10,000): '
                         f'median sum(expm1(X))={median:.1f}. Stop and inspect provenance.')
    return median


def run(args):
    import anndata as ad

    batches = [b.strip() for b in args.batches.split(',') if b.strip()]
    if len(batches) != 2:
        raise ValueError('For this case study, supply exactly two batch IDs')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    source = ad.read_h5ad(args.input, backed='r')
    try:
        if source.layers or source.raw is not None:
            # Such a file may still work, but it is not the audited data release.
            raise ValueError('Unexpected layers/raw found; please re-audit the input')
        if not source.var_names.is_unique:
            raise ValueError('Duplicate gene IDs prevent reliable gene reporting')
        keep, counts = select_cells(source.obs, args.bmp4, args.egf_bfgf,
                                    args.mix, args.ra, batches, args.min_cells)
        # Read only the selected 2,299 cells rather than the entire ~20k-cell object.
        cells = source[keep].to_memory()
    finally:
        source.file.close()

    median_total = check_expression_scale(cells.X)
    means_ref, means_treated, delta = means_and_differences(cells.X, cells.obs,
                                                            batches, args.ra)
    genes = cells.var_names.astype(str).to_numpy()
    per_batch = []
    for i, batch in enumerate(batches):
        per_batch.append(pd.DataFrame({
            'batch': batch, 'gene': genes,
            'mean_log1p_ra0': means_ref[i],
            'mean_log1p_ra_treated': means_treated[i],
            'delta_mean_log1p': delta[i],
        }))
    pd.concat(per_batch, ignore_index=True).to_csv(
        outdir / 'per_batch_gene_differences.tsv', sep='\t', index=False)

    mean_delta = delta.mean(axis=0)  # Equal batch weights, not cell-count weights.
    summary = pd.DataFrame({
        'gene': genes,
        'mean_delta_log1p_equal_batch_weight': mean_delta,
        f'delta_batch_{batches[0]}': delta[0],
        f'delta_batch_{batches[1]}': delta[1],
        'same_direction_in_both_batches': np.sign(delta[0]) == np.sign(delta[1]),
    })
    summary = summary.assign(abs_delta=np.abs(mean_delta)).sort_values(
        'abs_delta', ascending=False, kind='stable').drop(columns='abs_delta')
    summary.to_csv(outdir / 'gene_differences.tsv', sep='\t', index=False)
    counts.to_csv(outdir / 'group_counts.tsv', sep='\t', index=False)

    best = summary.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(best['gene'].astype(str), best['mean_delta_log1p_equal_batch_weight'])
    ax.axvline(0, color='gray', linewidth=0.8)
    ax.set_xlabel('RA-treated − RA-zero: equal-weight mean of batch-specific\n'
                  'differences in mean log1p-normalised expression')
    ax.set_title('Largest descriptive gene-expression contrasts (not differential expression)')
    fig.tight_layout()
    fig.savefig(outdir / 'top_genes.png', dpi=160)
    plt.close(fig)

    details = {
        'input_file': Path(args.input).name,
        'cells_selected': int(cells.n_obs),
        'genes': int(cells.n_vars),
        'fixed_doses': {'BMP4_ng_per_ml': args.bmp4,
                        'EGF_and_bFGF_ng_per_ml': args.egf_bfgf,
                        'Scriptaid_decitabine_mixture_uM': args.mix},
        'ra_reference_uM': 0,
        'ra_treated_uM': args.ra,
        'batches': batches,
        'counts_per_batch': counts.to_dict(orient='records'),
        'median_sum_expm1_X_selected_cells_sample': median_total,
        'expression_assumption': 'Input X appears log1p library-size-normalised to 10000',
        'analysis': 'Per-batch difference in mean(log1p(expression)); equal-weight average over batches',
        'independent_biological_replicates_verified': False,
        'statistical_tests': 'none',
        'warning': 'Batch labels are strata, NOT established biological replicates. '
                   'No p-values, log2 fold changes, or causal-effect claims.',
    }
    (outdir / 'contrast_summary.json').write_text(json.dumps(details, indent=2) + '\n')
    (outdir / 'README.md').write_text(
        '# Gehring–Pachter exploratory RA contrast\n\n'
        'Comparison: retinoic acid at 0 vs the specified treated dose, holding '
        'BMP4, EGF/bFGF, and Scriptaid/decitabine doses constant.\n\n'
        'Calculated within each specified batch, then averaged with equal batch '
        'weight. Gene differences are treated minus reference in '
        '**mean log1p-normalised expression**, NOT log2 fold changes.\n\n'
        'The source has no raw count layer, so count-based QC, re-normalisation, '
        'pseudobulk DE and p-values were NOT performed. Batch labels are not '
        'verified independent biological replicates. This is exploratory and '
        'must not be interpreted as a causal drug effect.\n', encoding='utf-8'
    )
    print(f'RA comparison finished: {cells.n_obs} cells, {cells.n_vars} genes; batches {batches}')
    print(counts.to_string(index=False))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', default='ra_analysis')
    p.add_argument('--bmp4', type=float, required=True)
    p.add_argument('--egf-bfgf', type=float, required=True)
    p.add_argument('--mix', type=float, required=True)
    p.add_argument('--ra', type=float, required=True)
    p.add_argument('--batches', required=True, help='Exactly two comma-separated batch labels')
    p.add_argument('--min-cells', type=int, default=20)
    return p.parse_args()


if __name__ == '__main__':
    run(parse_args())
