"""Recheck the descriptive Gehring–Pachter RA comparison and expression prevalence.

This is quality control, not a differential-expression or significance test.
The original .h5ad is never modified.
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

from compare_ra import select_cells, means_and_differences, check_expression_scale


def detected_fractions(X, obs, batches, treated_dose):
    """Proportion of cells with non-zero stored expression in each matched group."""
    batch_values = obs['batch'].astype(str).to_numpy()
    dose_values = pd.to_numeric(obs['dose_value_4'], errors='raise').to_numpy()
    results = {}
    for batch in batches:
        for dose, name in [(0, 'reference'), (treated_dose, 'treated')]:
            keep = ((batch_values == str(batch)) &
                    np.isclose(dose_values, dose, rtol=0, atol=1e-7))
            if not keep.any():
                raise ValueError(f'No cells for {name} in batch {batch}')
            block = X[keep, :]
            # The comparison > 0 excludes stored explicit zeros in sparse matrices.
            fractions = np.asarray((block > 0).mean(axis=0)).ravel()
            # Sparse reductions can exceed 1.0 by tiny floating-point error.
            fractions = np.clip(fractions, 0.0, 1.0)
            results[f'detected_{name}_batch_{batch}'] = fractions
    return pd.DataFrame(results)


def check_saved_differences(saved, gene_names, delta, batches, tolerance=2e-5):
    """Recalculate effects and ensure published TSV matches the source data."""
    if not gene_names.is_unique or not saved['gene'].is_unique:
        raise ValueError('Duplicate gene names in input or saved results')
    if set(saved['gene']) != set(gene_names):
        raise ValueError('Saved results contain a different set of genes')
    ordered = saved.set_index('gene').loc[gene_names]
    for i, batch in enumerate(batches):
        column = f'delta_batch_{batch}'
        if column not in ordered.columns:
            raise ValueError(f'Missing required column: {column}')
        if not np.allclose(ordered[column].to_numpy(dtype=float),
                           delta[i], rtol=0, atol=tolerance):
            raise ValueError(f'Saved effects for batch {batch} do not match recomputed data')
    averaged = delta.mean(axis=0)
    column = 'mean_delta_log1p_equal_batch_weight'
    if not np.allclose(ordered[column].to_numpy(dtype=float), averaged,
                       rtol=0, atol=tolerance):
        raise ValueError('Saved mean effects do not match recomputed data')
    return ordered


def summarize_directions(delta, min_magnitude=0.1):
    """Simple concordance descriptives; no hypothesis testing."""
    if delta.shape[0] != 2:
        raise ValueError('Expected exactly two batches')
    d0, d1 = delta
    strong = (np.abs(d0) >= min_magnitude) & (np.abs(d1) >= min_magnitude)
    same = np.sign(d0) == np.sign(d1)
    if np.std(d0) > 0 and np.std(d1) > 0:
        corr = float(np.corrcoef(d0, d1)[0, 1])
    else:
        corr = None
    return {
        'all_gene_batch_effect_pearson_r_descriptive': corr,
        'n_genes_with_abs_delta_at_least_0_1_in_each_batch': int(strong.sum()),
        'n_of_those_with_same_direction': int((strong & same).sum()),
        'fraction_of_those_with_same_direction':
            float((strong & same).sum() / strong.sum()) if strong.any() else None,
    }


def run(args):
    import anndata as ad

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    details = json.loads(Path(args.summary).read_text())
    fixed = details['fixed_doses']
    batches = [str(x) for x in details['batches']]
    if len(batches) != 2 or len(set(batches)) != 2:
        raise ValueError('Two distinct batches are required')
    ra = float(details['ra_treated_uM'])
    saved = pd.read_csv(args.genes, sep='\t')
    required = {'gene', 'mean_delta_log1p_equal_batch_weight'}
    if not required.issubset(saved.columns):
        raise ValueError('Saved results lack required gene columns')

    original = ad.read_h5ad(args.input, backed='r')
    try:
        keep, counts = select_cells(
            original.obs,
            fixed['BMP4_ng_per_ml'],
            fixed['EGF_and_bFGF_ng_per_ml'],
            fixed['Scriptaid_decitabine_mixture_uM'],
            ra, batches, min_cells=1)
        data = original[keep].to_memory()
    finally:
        original.file.close()

    expected_counts = pd.DataFrame(details['counts_per_batch'])
    if not counts.equals(expected_counts):
        raise ValueError('Group counts do not match saved contrast summary')
    if data.n_obs != details['cells_selected'] or data.n_vars != details['genes']:
        raise ValueError('Dataset dimensions do not match saved contrast summary')
    median_total = check_expression_scale(data.X)
    _, _, delta = means_and_differences(data.X, data.obs, batches, ra)
    genes = pd.Index(data.var_names.astype(str), name='gene')
    ordered = check_saved_differences(saved, genes, delta, batches)

    fractions = detected_fractions(data.X, data.obs, batches, ra)
    fractions.index = genes
    quality = ordered.join(fractions)
    fraction_cols = list(fractions.columns)
    quality['min_detection_fraction_across_groups'] = quality[fraction_cols].min(axis=1)
    quality.reset_index().to_csv(outdir / 'gene_quality.tsv', sep='\t', index=False)

    # Display detection prevalence for top genes, without dropping genes or
    # retrospectively changing the selection or ranking used in the contrast.
    top_genes = saved.head(20)['gene'].tolist()
    quality.loc[top_genes].reset_index().to_csv(
        outdir / 'top20_gene_quality.tsv', sep='\t', index=False)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(delta[0], delta[1], s=5, alpha=0.2)
    limits = max(float(np.max(np.abs(delta))), 0.1) * 1.05
    ax.plot([-limits, limits], [-limits, limits], linestyle='--', linewidth=1)
    ax.axhline(0, linewidth=0.6)
    ax.axvline(0, linewidth=0.6)
    ax.set_xlim(-limits, limits)
    ax.set_ylim(-limits, limits)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(f'Batch {batches[0]}: RA-treated − RA-zero (mean log1p)')
    ax.set_ylabel(f'Batch {batches[1]}: RA-treated − RA-zero (mean log1p)')
    ax.set_title('Descriptive concordance across batches')
    for gene in top_genes[:6]:
        pos = genes.get_loc(gene)
        ax.annotate(gene, (delta[0, pos], delta[1, pos]), fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / 'batch_concordance.png', dpi=160)
    plt.close(fig)

    summary = {
        'input_file': Path(args.input).name,
        'source_contrast_summary': Path(args.summary).name,
        'cells_verified': int(data.n_obs),
        'genes_verified': int(data.n_vars),
        'batches': batches,
        'group_counts_verified': True,
        'saved_gene_effects_recomputed_and_verified': True,
        'median_sum_expm1_X_selected_cells_sample': median_total,
        'detection_definition': 'fraction of cells with stored X > 0; X is presumed log1p-normalised',
        'min_detection_fraction_is_not_a_filter': True,
        'biological_replicates_verified': False,
        'inference': 'descriptive only; no independent biological replicates established',
        **summarize_directions(delta),
    }
    (outdir / 'quality_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Quality check passed: {data.n_obs} cells, {data.n_vars} genes; '
          f'{len(batches)} batches, output: {outdir}')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, help='Original Gehring .h5ad')
    p.add_argument('--summary', required=True, help='RA contrast_summary.json')
    p.add_argument('--genes', required=True, help='RA gene_differences.tsv')
    p.add_argument('--outdir', default='ra_quality')
    return p.parse_args()


if __name__ == '__main__':
    run(parse_args())
