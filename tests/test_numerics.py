import numpy as np
import pytest
from scipy import sparse

from scripts.numerics import cell_metrics, mean_expression, top_differences, validate_matrix


@pytest.mark.parametrize('as_sparse', [False, True])
def test_count_validation_and_metrics(as_sparse):
    matrix = np.array([[0, 2, 4], [1, 0, 0]], dtype=float)
    if as_sparse:
        matrix = sparse.csr_matrix(matrix)
    validate_matrix(matrix, 'counts')
    totals, detected = cell_metrics(matrix)
    np.testing.assert_array_equal(totals, [6, 1])
    np.testing.assert_array_equal(detected, [2, 1])
    np.testing.assert_allclose(mean_expression(matrix), [0.5, 1., 2.])


@pytest.mark.parametrize('as_sparse', [False, True])
def test_negative_and_fractional_counts_are_rejected(as_sparse):
    for values in [np.array([[1.2, 3.0]]), np.array([[0., -1.]]), np.array([[np.nan, 1.]])]:
        values = sparse.csr_matrix(values) if as_sparse else values
        with pytest.raises(ValueError):
            validate_matrix(values, 'counts')


def test_log1p_input_can_be_noninteger_but_not_negative():
    validate_matrix(np.array([[0.0, 1.3]]), 'log1p')
    with pytest.raises(ValueError):
        validate_matrix(np.array([[-1.0]]), 'log1p')


def test_descriptive_differences_are_not_fold_changes():
    treatment = sparse.csr_matrix([[2., 1.], [4., 3.]])
    control = sparse.csr_matrix([[1., 2.], [1., 2.]])
    rows = top_differences(treatment, control, ['A', 'B'], top_n=2)
    assert rows[0][0] == 'A'
    assert rows[0][3] == pytest.approx(2.0)
    assert rows[1][3] == pytest.approx(0.0)
