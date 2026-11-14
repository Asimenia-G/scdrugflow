"""Optional end-to-end Python test; Nextflow orchestration is a separate check."""

from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

pytest.importorskip("anndata")
pytest.importorskip("scanpy")


def test_toy_example(tmp_path):
    root = Path(__file__).resolve().parents[1]
    dataset = tmp_path / 'toy.h5ad'
    inspection = tmp_path / 'inspection.json'
    out = tmp_path / 'analysis'
    subprocess.run([
        sys.executable, str(root / 'scripts/make_toy.py'), '--output', str(dataset)
    ], check=True)
    subprocess.run([
        sys.executable, str(root / 'scripts/inspect_dataset.py'), '--input', str(dataset),
        '--output', str(inspection)
    ], check=True)
    subprocess.run([
        sys.executable, str(root / 'scripts/analyse.py'), '--input', str(dataset),
        '--inspection', str(inspection), '--outdir', str(out),
        '--condition-column', 'perturbation', '--control-label', 'control',
        '--matrix', 'X', '--input-kind', 'counts', '--min-genes', '10',
        '--min-counts', '30', '--min-cells-per-gene', '3', '--min-group-cells', '10',
    ], check=True)
    for name in ['summary.json', 'pca.png', 'condition_counts.tsv',
                 'top_expression_differences.tsv', 'processed.h5ad']:
        assert (out / name).is_file(), name

    summary = pd.read_csv(out / 'condition_counts.tsv', sep='\t')
    assert set(summary['condition']) == {'control', 'drug_A', 'drug_B'}
    effects = pd.read_csv(out / 'top_expression_differences.tsv', sep='\t')
    assert not effects.empty
    drug_a = effects[effects['condition'] == 'drug_A']
    drug_b = effects[effects['condition'] == 'drug_B']
    assert 'Gene000' in set(drug_a['gene'])
    assert 'Gene005' in set(drug_b['gene'])
    assert float(drug_a.loc[drug_a['gene'] == 'Gene000', 'delta_mean_log1p'].iloc[0]) > 0
    assert float(drug_b.loc[drug_b['gene'] == 'Gene005', 'delta_mean_log1p'].iloc[0]) > 0
