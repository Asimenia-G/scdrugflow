"""Audit Gehring–Pachter dose combinations, batches, and expression scale.

This does not perform differential expression or modify the input file.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

DOSES = ['dose_value', 'dose_value_2', 'dose_value_3', 'dose_value_4']
BACKGROUND = DOSES[:3]


def make_tables(obs):
    required = DOSES + ['batch']
    missing = [column for column in required if column not in obs]
    if missing:
        raise ValueError(f'Missing metadata columns: {missing}')

    table = obs[required].copy()
    if table.isna().any().any():
        raise ValueError('Missing values in dose or batch metadata')
    table['batch'] = table['batch'].astype(str)

    conditions = (table.groupby(DOSES, observed=True)
                  .size().reset_index(name='n_cells')
                  .sort_values('n_cells', ascending=False))
    by_batch = (table.groupby(DOSES + ['batch'], observed=True)
                .size().reset_index(name='n_cells'))

    pairs = []
    for keys, group in by_batch.groupby(BACKGROUND + ['batch'], observed=True):
        baseline = group.loc[group['dose_value_4'] == 0]
        if len(baseline) != 1:
            continue
        n_reference = int(baseline['n_cells'].iloc[0])
        for _, treated in group.loc[group['dose_value_4'] > 0].iterrows():
            pairs.append(dict(zip(BACKGROUND + ['batch'], keys)),
                         )
            pairs[-1].update({
                'ra_reference_uM': 0,
                'ra_treated_uM': int(treated['dose_value_4']),
                'n_reference': n_reference,
                'n_treated': int(treated['n_cells']),
                'min_group_size': min(n_reference, int(treated['n_cells'])),
            })

    pairs = pd.DataFrame(pairs, columns=BACKGROUND + ['batch',
        'ra_reference_uM', 'ra_treated_uM', 'n_reference', 'n_treated', 'min_group_size'])
    if not pairs.empty:
        pairs = pairs.sort_values('min_group_size', ascending=False)
    return conditions, by_batch, pairs


def expression_diagnostics(adata, n=200):
    """Sample X only; expm1 totals are suggestive, not proof of log1p input."""
    x = adata.X[: min(n, adata.n_obs), :]
    if sparse.issparse(x):
        x = x.tocsr(copy=True)
        vals = x.data
        logged = x.copy()
        logged.data = np.expm1(logged.data.astype(np.float64))
        totals = np.asarray(logged.sum(axis=1)).ravel()
    else:
        x = np.asarray(x)
        vals = x.ravel()
        totals = np.expm1(x.astype(np.float64)).sum(axis=1)
    nonzero = vals[vals != 0]
    result = {
        'n_cells_sampled': int(len(totals)),
        'min_nonzero_X': float(nonzero.min()) if len(nonzero) else None,
        'max_nonzero_X': float(nonzero.max()) if len(nonzero) else None,
        'nonnegative': bool(np.all(vals >= 0)),
        'integer_like': bool(np.allclose(vals, np.rint(vals), atol=1e-6, rtol=0)),
        'median_row_sum_expm1_X': float(np.median(totals)),
        'min_row_sum_expm1_X': float(np.min(totals)),
        'max_row_sum_expm1_X': float(np.max(totals)),
        'uns_keys': list(adata.uns.keys()),
        'raw_present': adata.raw is not None,
        'caution': 'exp(X)-1 row sums near a constant suggest log1p-normalised values, but do not establish provenance.',
    }
    return result


def main():
    import anndata as ad

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--outdir', default='results/gehring_audit')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(args.input, backed='r')
    try:
        conditions, by_batch, pairs = make_tables(a.obs)
        diagnostic = expression_diagnostics(a)
    finally:
        a.file.close()

    conditions.to_csv(outdir / 'dose_conditions.tsv', sep='\t', index=False)
    by_batch.to_csv(outdir / 'dose_conditions_by_batch.tsv', sep='\t', index=False)
    pairs.to_csv(outdir / 'ra_matched_pairs.tsv', sep='\t', index=False)
    (outdir / 'matrix_diagnostics.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
    print(f'Unique four-factor combinations: {len(conditions)}')
    print(f'RA contrasts matched on the other doses AND batch: {len(pairs)}')
    print('Matrix diagnostics:')
    print(json.dumps(diagnostic, indent=2))
    print('Top potentially comparable RA pairs (not yet validated replicates):')
    print(pairs.head(12).to_string(index=False) if not pairs.empty else 'None found')
    print(f'Full reports: {outdir.resolve()}')


if __name__ == '__main__':
    main()
