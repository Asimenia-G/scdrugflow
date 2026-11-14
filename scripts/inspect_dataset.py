"""Examine an .h5ad file without assuming its metadata or count scale."""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


def inspect_dataset(path):
    # backed mode avoids loading the full expression matrix just for inspection.
    adata = ad.read_h5ad(path, backed='r')
    try:
        matrix = adata.X
        sample = matrix[: min(100, adata.n_obs), :]
        values = sample.data if sparse.issparse(sample) else np.asarray(sample).ravel()
        present = values[np.abs(values) > 0]

        obs_info = {}
        for column in adata.obs.columns:
            series = adata.obs[column]
            detail = {
                'dtype': str(series.dtype),
                'missing': int(series.isna().sum()),
                'unique_values': int(series.nunique(dropna=True)),
            }
            if series.nunique(dropna=True) <= 1000:
                detail['most_common_values'] = {
                    str(k): int(v)
                    for k, v in series.astype('string').value_counts().head(12).items()
                }
            obs_info[column] = detail

        return {
            'file': Path(path).name,
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars),
            'obs_columns': obs_info,
            'var_columns': list(adata.var.columns),
            'layers': list(adata.layers.keys()),
            'X_dtype': str(matrix.dtype),
            'X_sample_first_100_cells': {
                'n_nonzero_sampled': int(present.size),
                'min_nonzero': float(present.min()) if present.size else None,
                'max_nonzero': float(present.max()) if present.size else None,
                'nonnegative': bool((values >= 0).all()),
                'integer_like': bool(
                    np.isfinite(values).all()
                    and np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
                ),
            },
            'important': (
                'Integer-like values do NOT prove that X contains raw counts. '
                'Inspect the source documentation and any layers before analysing.'
            ),
        }
    finally:
        adata.file.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', default='inspection.json')
    args = parser.parse_args()
    report = inspect_dataset(args.input)
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'obs_columns'}, indent=2))
    print('Metadata inspection saved to', args.output)


if __name__ == '__main__':
    main()
