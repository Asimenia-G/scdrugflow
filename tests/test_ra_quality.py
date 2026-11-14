"""Pure numerical checks for RA quality-control calculations."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from ra_quality import detected_fractions, check_saved_differences, summarize_directions


def test_detected_fractions_sparse_and_dense_agree():
    obs = pd.DataFrame({'batch': ['0', '0', '0', '0', '1', '1', '1', '1'],
                        'dose_value_4': [0, 0, 10, 10] * 2})
    X = np.array([[0, 1], [2, 0], [1, 0], [3, 1],
                  [0, 0], [2, 0], [1, 1], [0, 1]], dtype=float)
    dense = detected_fractions(X, obs, ['0', '1'], 10)
    sparse_result = detected_fractions(sparse.csr_matrix(X), obs, ['0', '1'], 10)
    pd.testing.assert_frame_equal(dense, sparse_result)
    assert dense.to_numpy().min() >= 0.0
    assert dense.to_numpy().max() <= 1.0
    assert sparse_result.to_numpy().min() >= 0.0
    assert sparse_result.to_numpy().max() <= 1.0
    assert dense.loc[0, 'detected_reference_batch_0'] == 0.5
    assert dense.loc[1, 'detected_treated_batch_1'] == 1.0


def test_saved_effects_are_recomputed():
    genes = pd.Index(['A', 'B'], name='gene')
    delta = np.array([[0.4, -0.2], [0.2, -0.4]])
    saved = pd.DataFrame({'gene': ['B', 'A'],
                          'delta_batch_0': [-0.2, 0.4],
                          'delta_batch_1': [-0.4, 0.2],
                          'mean_delta_log1p_equal_batch_weight': [-0.3, 0.3]})
    ordered = check_saved_differences(saved, genes, delta, ['0', '1'])
    assert ordered.index.tolist() == ['A', 'B']
    saved.loc[0, 'delta_batch_0'] = 99
    with pytest.raises(ValueError, match='do not match'):
        check_saved_differences(saved, genes, delta, ['0', '1'])


def test_direction_summary_is_descriptive():
    delta = np.array([[0.4, -0.3, 0.01], [0.5, 0.4, -0.2]])
    info = summarize_directions(delta)
    assert info['n_genes_with_abs_delta_at_least_0_1_in_each_batch'] == 2
    assert info['n_of_those_with_same_direction'] == 1
    assert info['fraction_of_those_with_same_direction'] == 0.5


def test_missing_group_is_rejected():
    obs = pd.DataFrame({'batch': ['0'], 'dose_value_4': [0]})
    with pytest.raises(ValueError, match='No cells'):
        detected_fractions(np.array([[1]]), obs, ['0'], 10)


def test_quality_command_with_small_h5ad(tmp_path):
    """Runs only when AnnData is installed; Hypatia has that dependency."""
    import json
    import subprocess

    ad = pytest.importorskip('anndata')
    root = Path(__file__).resolve().parents[1]
    obs_rows, X_rows = [], []
    for batch in ('0', '1'):
        for ra in (0, 10):
            for _ in range(5):
                obs_rows.append([8, 200, 0, ra, batch])
                counts = [2500, 5500, 2000] if ra == 0 else [6500, 1500, 2000]
                X_rows.append(np.log1p(counts))
    obs = pd.DataFrame(obs_rows, columns=[
        'dose_value', 'dose_value_2', 'dose_value_3', 'dose_value_4', 'batch'])
    obs.index = [f'cell_{i}' for i in range(len(obs))]
    a = ad.AnnData(np.array(X_rows, dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=['GeneA', 'GeneB', 'GeneC']))
    input_file = tmp_path / 'test.h5ad'
    a.write_h5ad(input_file)

    # compare_ra.py is an existing file in the user's scDrugFlow project.
    out_compare = tmp_path / 'comparison'
    subprocess.run([
        sys.executable, str(root / 'scripts/compare_ra.py'), '--input', str(input_file),
        '--outdir', str(out_compare), '--bmp4', '8', '--egf-bfgf', '200',
        '--mix', '0', '--ra', '10', '--batches', '0,1', '--min-cells', '3',
    ], check=True)
    out_quality = tmp_path / 'qc'
    subprocess.run([
        sys.executable, str(root / 'scripts/ra_quality.py'), '--input', str(input_file),
        '--summary', str(out_compare / 'contrast_summary.json'),
        '--genes', str(out_compare / 'gene_differences.tsv'),
        '--outdir', str(out_quality),
    ], check=True)
    summary = json.loads((out_quality / 'quality_summary.json').read_text())
    assert summary['saved_gene_effects_recomputed_and_verified'] is True
    assert summary['cells_verified'] == 20
    assert (out_quality / 'batch_concordance.png').exists()
    quality = pd.read_csv(out_quality / 'gene_quality.tsv', sep='\t')
    assert quality.shape[0] == 3
    assert quality['min_detection_fraction_across_groups'].between(0, 1).all()
