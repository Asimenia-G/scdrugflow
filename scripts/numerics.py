"""Pure numerical helpers, tested for dense and sparse matrices."""

import numpy as np
from scipy import sparse


def stored_values(matrix):
    """Only stored entries matter for sparse matrices (implicit zeros are valid)."""
    return matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()


def validate_matrix(matrix, kind):
    """Fail rather than silently renormalise an unsuitable expression matrix.

    Integer-like values are necessary but NOT sufficient proof of raw UMI counts.
    Users must confirm the source's preprocessing/provenance themselves.
    """
    values = stored_values(matrix)
    if values.size and (not np.isfinite(values).all() or (values < 0).any()):
        raise ValueError('Expression matrix contains negative or non-finite values')
    if kind == 'counts' and values.size and not np.allclose(
        values, np.rint(values), rtol=0, atol=1e-6
    ):
        raise ValueError(
            'Matrix is not integer-like. Check inspection.json and the data source. '
            'Choose a raw-count layer, or set --input_kind log1p ONLY if the '
            'matrix is already log1p-normalised expression.'
        )


def cell_metrics(matrix):
    """Return per-cell library sizes and number of detected genes."""
    total = np.asarray(matrix.sum(axis=1)).ravel()
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
    return total, detected


def mean_expression(matrix):
    """Per-gene arithmetic mean, for a matrix of cells by genes."""
    if matrix.shape[0] == 0:
        raise ValueError('Cannot compute expression mean for zero cells')
    return np.asarray(matrix.mean(axis=0)).ravel()


def top_differences(group_matrix, control_matrix, genes, top_n=20):
    """Rank differences of mean log1p-normalised expression. Not log-fold changes."""
    group_mean = mean_expression(group_matrix)
    control_mean = mean_expression(control_matrix)
    differences = group_mean - control_mean
    idx = np.argsort(-np.abs(differences), kind='stable')[:top_n]
    return [
        (str(genes[i]), float(group_mean[i]), float(control_mean[i]), float(differences[i]))
        for i in idx
    ]
