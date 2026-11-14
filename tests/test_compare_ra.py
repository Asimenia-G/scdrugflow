"""Tests with simple synthetic expression and treatment metadata."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_ra import select_cells, means_and_differences, check_expression_scale


def make_obs():
    rows = []
    for batch in ('0', '1'):
        for ra in (0, 10):
            for _ in range(4):
                rows.append([8, 200, 0, ra, batch])
    rows.append([8, 40, 0, 10, '0'])  # wrong background
    rows.append([8, 200, 0.2, 10, '1'])  # wrong mixture
    return pd.DataFrame(rows, columns=['dose_value', 'dose_value_2',
                                       'dose_value_3', 'dose_value_4', 'batch'])


def test_exact_matched_selection_and_counts():
    obs = make_obs()
    keep, counts = select_cells(obs, 8, 200, 0, 10, ['0', '1'], min_cells=3)
    assert int(keep.sum()) == 16
    assert counts[['n_ra_0', 'n_ra_treated']].values.tolist() == [[4, 4], [4, 4]]


def test_sparse_and_dense_contrasts_agree():
    obs = make_obs().iloc[:16].copy()
    matrix = np.ones((16, 3), dtype=float)
    matrix[obs['dose_value_4'].to_numpy() == 10, 0] = 3.0
    for X in (matrix, sparse.csr_matrix(matrix)):
        reference, treated, difference = means_and_differences(X, obs, ['0', '1'], 10)
        assert difference.shape == (2, 3)
        np.testing.assert_allclose(difference[:, 0], [2, 2])
        np.testing.assert_allclose(difference[:, 1:], 0)
        np.testing.assert_allclose(reference[:, 0], 1)
        np.testing.assert_allclose(treated[:, 0], 3)


def test_logged_expression_guard():
    X = np.log1p(np.array([[2500, 7500], [4000, 6000]], dtype=float))
    assert 9999.99 < check_expression_scale(sparse.csr_matrix(X)) < 10000.01


def test_missing_comparison_fails():
    obs = make_obs()
    try:
        select_cells(obs, 8, 200, 0, 2, ['0', '1'], min_cells=3)
    except ValueError as exc:
        assert 'require at least' in str(exc)
    else:
        raise AssertionError('Missing comparison must fail')


def test_realistic_h5ad_command(tmp_path):
    """Exercise the actual command when AnnData is installed (as on Hypatia)."""
    import subprocess
    import pytest

    ad = pytest.importorskip('anndata')
    obs = []
    X = []
    for batch in ('0', '1'):
        for ra in (0, 10):
            for i in range(5):
                obs.append([8, 200, 0, ra, batch])
                counts = [2500, 5500, 2000] if ra == 0 else [6500, 1500, 2000]
                X.append(np.log1p(counts))
    metadata = pd.DataFrame(obs, columns=['dose_value', 'dose_value_2',
                                          'dose_value_3', 'dose_value_4', 'batch'])
    metadata.index = [f'cell_{i}' for i in range(len(metadata))]
    a = ad.AnnData(np.array(X, dtype=np.float32), obs=metadata,
                   var=pd.DataFrame(index=['GeneA', 'GeneB', 'GeneC']))
    input_path = tmp_path / 'toy_ra.h5ad'
    outdir = tmp_path / 'results'
    a.write_h5ad(input_path)
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'compare_ra.py'
    subprocess.run([sys.executable, str(script), '--input', str(input_path),
                    '--outdir', str(outdir), '--bmp4', '8', '--egf-bfgf', '200',
                    '--mix', '0', '--ra', '10', '--batches', '0,1',
                    '--min-cells', '3'], check=True)
    result = pd.read_csv(outdir / 'gene_differences.tsv', sep='\t')
    assert result.loc[result['gene'] == 'GeneA',
                      'mean_delta_log1p_equal_batch_weight'].iloc[0] > 0
    assert (outdir / 'top_genes.png').is_file()
    assert json_from_path(outdir / 'contrast_summary.json')['cells_selected'] == 20


def json_from_path(path):
    import json
    return json.loads(path.read_text())
