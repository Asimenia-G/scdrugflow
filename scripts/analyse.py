"""Conservative, descriptive analysis of a single-cell perturbation dataset.

No p-values, cell-level DE tests, or claims of causal drug response are produced.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from numerics import cell_metrics, top_differences, validate_matrix


def selected_matrix(adata, name):
    if name == 'X':
        return adata.X
    if name not in adata.layers:
        raise ValueError(f'Layer {name!r} not present; available: {list(adata.layers)}')
    return adata.layers[name]


def filter_by_metadata(adata, condition_column, subset_column, subset_value):
    if condition_column not in adata.obs:
        raise ValueError(
            f'Condition column {condition_column!r} not found. '
            'Run the inspect mode and review inspection.json.'
        )
    if bool(subset_column) != bool(subset_value):
        raise ValueError('Supply both --subset-column and --subset-value, or neither')
    if subset_column:
        if subset_column not in adata.obs:
            raise ValueError(f'Subset column {subset_column!r} not present')
        subset = adata.obs[subset_column].astype('string').str.strip()
        adata = adata[subset.eq(subset_value).fillna(False).to_numpy(bool)].copy()

    labels = adata.obs[condition_column].astype('string').str.strip()
    adata = adata[labels.notna().to_numpy(bool) & labels.ne('').fillna(False).to_numpy(bool)].copy()
    # Convert to ordinary Python strings for portable .h5ad writing.
    adata.obs[condition_column] = (
        adata.obs[condition_column].astype('string').str.strip().astype(str)
    )
    if adata.n_obs == 0:
        raise ValueError('No cells left after metadata filtering')
    return adata


def draw_condition_counts(counts, output):
    counts = counts.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(3, 0.42 * len(counts) + 1)))
    ax.barh(range(len(counts)), counts.to_numpy())
    ax.set_yticks(range(len(counts)), [str(c) for c in counts.index])
    ax.set_xlabel('Cells retained after filtering')
    ax.set_title('Most frequent conditions (up to 12)')
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def draw_pca(adata, condition_column, output):
    # Sparse-friendly Scanpy PCA on all retained genes, descriptive only.
    if min(adata.n_obs, adata.n_vars) <= 3:
        return False
    sc.pp.pca(adata, n_comps=2, svd_solver='arpack', random_state=42)
    coords = adata.obsm['X_pca']
    labels = adata.obs[condition_column].astype(str)
    most_common = labels.value_counts().head(7).index
    display_labels = labels.where(labels.isin(most_common), 'Other')

    fig, ax = plt.subplots(figsize=(8, 6))
    for label in sorted(display_labels.unique()):
        mask = (display_labels == label).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=8, alpha=0.7, label=label)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('PCA of log1p-normalised expression')
    ax.legend(loc='best', fontsize=8, markerscale=2)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)
    return True


def analyse(args):
    if args.min_group_cells < 1 or args.min_cells_per_gene < 1:
        raise ValueError('Minimum cell thresholds must be positive integers')
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input)
    if not adata.var_names.is_unique:
        raise ValueError('Gene IDs are not unique; resolve them before analysis')

    original_shape = [int(adata.n_obs), int(adata.n_vars)]
    adata = filter_by_metadata(
        adata, args.condition_column, args.subset_column, args.subset_value
    )
    metadata_filtered_shape = [int(adata.n_obs), int(adata.n_vars)]

    matrix = selected_matrix(adata, args.matrix)
    validate_matrix(matrix, args.input_kind)
    # Work on a fresh copy so that X is always the explicitly selected expression matrix.
    adata.X = matrix.copy()
    if sparse.issparse(adata.X):
        adata.X = adata.X.tocsr().astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)

    if args.input_kind == 'counts':
        if args.min_genes < 1 or args.min_counts < 1:
            raise ValueError('For raw counts, min_genes and min_counts must be >= 1')
        totals, detected = cell_metrics(adata.X)
        keep = (totals >= args.min_counts) & (detected >= args.min_genes)
        adata = adata[keep].copy()
        if adata.n_obs == 0:
            raise ValueError('All cells failed count-based QC; review thresholds')
        sc.pp.filter_genes(adata, min_cells=args.min_cells_per_gene)
        if adata.n_vars < 3:
            raise ValueError('Too few genes retained to analyse')
        # Preserve filtered, unnormalised counts for auditability.
        adata.layers['counts'] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        preprocessing = 'Validated integer-like counts; cell/gene QC; library-size normalisation to 10,000; log1p'
    else:
        # Already log1p-normalised: do NOT infer raw-count QC metrics or renormalise.
        if adata.n_vars < 3:
            raise ValueError('Too few genes to analyse')
        preprocessing = 'Input declared already log1p-normalised; count-based QC and normalisation skipped'

    counts = adata.obs[args.condition_column].astype(str).value_counts()
    if args.control_label not in counts:
        raise ValueError(
            f'Control label {args.control_label!r} not found after filtering. '
            f'Available conditions include: {list(counts.head(15).index)}'
        )
    if counts[args.control_label] < args.min_group_cells:
        raise ValueError('Not enough control cells for the requested --min-group-cells')

    (outdir / 'inspection.json').write_text(
        Path(args.inspection).read_text(encoding='utf-8'), encoding='utf-8'
    )
    pd.DataFrame({
        'condition': counts.index,
        'n_cells': counts.to_numpy(),
    }).to_csv(outdir / 'condition_counts.tsv', sep='\t', index=False)
    draw_condition_counts(counts, outdir / 'condition_counts.png')
    pca_created = draw_pca(adata, args.condition_column, outdir / 'pca.png')

    labels = adata.obs[args.condition_column].astype(str).to_numpy()
    control_mask = labels == args.control_label
    control_matrix = adata.X[control_mask]
    rows = []
    skipped = []
    for group, n_cells in counts.items():
        if group == args.control_label:
            continue
        if n_cells < args.min_group_cells:
            skipped.append(str(group))
            continue
        group_matrix = adata.X[labels == group]
        for gene, mean_group, mean_control, delta in top_differences(
            group_matrix, control_matrix, adata.var_names, top_n=20
        ):
            rows.append({
                'condition': str(group), 'n_group_cells': int(n_cells),
                'n_control_cells': int(counts[args.control_label]),
                'gene': gene, 'mean_log1p_group': mean_group,
                'mean_log1p_control': mean_control, 'delta_mean_log1p': delta,
            })

    result_cols = [
        'condition', 'n_group_cells', 'n_control_cells', 'gene',
        'mean_log1p_group', 'mean_log1p_control', 'delta_mean_log1p',
    ]
    pd.DataFrame(rows, columns=result_cols).to_csv(
        outdir / 'top_expression_differences.tsv', sep='\t', index=False
    )

    summary = {
        'input_file': Path(args.input).name,
        'original_shape_cells_genes': original_shape,
        'after_metadata_filter_cells_genes': metadata_filtered_shape,
        'after_qc_cells_genes': [int(adata.n_obs), int(adata.n_vars)],
        'condition_column': args.condition_column,
        'control_label': args.control_label,
        'selected_matrix': args.matrix,
        'input_kind_declared_by_user': args.input_kind,
        'preprocessing': preprocessing,
        'min_group_cells': args.min_group_cells,
        'n_conditions_including_control': len(counts),
        'n_conditions_compared': len(set(row['condition'] for row in rows)),
        'groups_excluded_for_small_cell_count': skipped,
        'pca_generated': pca_created,
        'synthetic_data': bool(adata.uns.get('synthetic_example', False)),
    }
    (outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (outdir / 'README.md').write_text(
        '# scDrugFlow results\n\n'
        f'- Input: `{Path(args.input).name}`\n'
        f'- Cells: {original_shape[0]} original, {adata.n_obs} retained\n'
        f'- Genes: {original_shape[1]} original, {adata.n_vars} retained\n'
        f'- Control: `{args.control_label}`\n'
        f'- Preprocessing: {preprocessing}\n'
        f'- Matched/stratified analysis: **not performed automatically**\n'
        '- Expression differences: mean(log1p-normalised expression) for the '
        'condition minus the control. These are **not** log2 fold changes, '
        'differential-expression tests or causal drug effects.\n'
        '- No p-values, false-discovery rates, or pathway-enrichment claims are produced.\n'
        '- Interpret comparisons only after reviewing dose, batch, cell identity '
        'and independent biological replication in the original study.\n'
        + ('- **Toy synthetic data: results have no biological meaning.**\n' if summary['synthetic_data'] else ''),
        encoding='utf-8',
    )
    adata.write_h5ad(outdir / 'processed.h5ad', compression='gzip')
    print(f'Finished exploratory analysis: {adata.n_obs} cells, {adata.n_vars} genes')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--inspection', required=True)
    parser.add_argument('--outdir', default='analysis')
    parser.add_argument('--condition-column', required=True)
    parser.add_argument('--control-label', required=True)
    parser.add_argument('--matrix', default='X')
    parser.add_argument('--input-kind', choices=['counts', 'log1p'], required=True)
    parser.add_argument('--subset-column', default='')
    parser.add_argument('--subset-value', default='')
    parser.add_argument('--min-genes', type=int, default=100)
    parser.add_argument('--min-counts', type=int, default=300)
    parser.add_argument('--min-cells-per-gene', type=int, default=3)
    parser.add_argument('--min-group-cells', type=int, default=20)
    return parser.parse_args()


if __name__ == '__main__':
    analyse(parse_args())
