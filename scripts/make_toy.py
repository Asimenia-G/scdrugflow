"""Create a tiny, synthetic example for tests. It is NOT biological evidence."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def create_toy(path):
    rng = np.random.default_rng(42)
    groups = np.repeat(['control', 'drug_A', 'drug_B'], 40)
    n_cells, n_genes = len(groups), 60
    counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.int32)
    # Synthetic signal to check code; not actual pharmacological response.
    counts[groups == 'drug_A', 0:5] += 8
    counts[groups == 'drug_B', 5:10] += 8
    obs = pd.DataFrame(
        {'perturbation': groups, 'batch': np.tile(np.repeat(['B1', 'B2', 'B3', 'B4'], 10), 3)},
        index=[f'cell_{i:03}' for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f'Gene{i:03}' for i in range(n_genes)])
    data = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    data.uns['synthetic_example'] = True
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(path)
    print('Created SYNTHETIC test data:', path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='data/toy.h5ad')
    args = parser.parse_args()
    create_toy(args.output)
