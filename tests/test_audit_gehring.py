"""Unit tests for exact-factor and batch matching."""

import importlib.util
from pathlib import Path

import pandas as pd

SPEC = importlib.util.spec_from_file_location(
    'audit_gehring', Path(__file__).resolve().parents[1] / 'scripts' / 'audit_gehring.py'
)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_same_background_and_batch_only():
    rows = [
        # matched RA doses in batch 0
        (0, 40, 0, 0, '0'), (0, 40, 0, 0, '0'), (0, 40, 0, 2, '0'),
        # not a matched pair: different batch
        (0, 40, 0, 10, '1'),
        # not a matched pair: different BMP4 dose
        (8, 40, 0, 2, '0'),
    ]
    obs = pd.DataFrame(rows, columns=module.DOSES + ['batch'])
    conditions, by_batch, pairs = module.make_tables(obs)
    assert len(conditions) == 4
    assert len(by_batch) == 4
    assert len(pairs) == 1
    assert pairs.iloc[0]['ra_treated_uM'] == 2
    assert pairs.iloc[0]['n_reference'] == 2
    assert pairs.iloc[0]['n_treated'] == 1
