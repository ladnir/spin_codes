"""Tiny exact/toy/mock tests; no GF256/a3 catalog or actual gate run."""
from collections import Counter
from itertools import product
from math import log2
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_categories as cats
import hd_lambda_shear_exact3_control as gate
from hd_pair_identity_kernel import mask_table_from_entries


MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


def output_mask(forms, values):
    mask = 0
    for row, form in enumerate(forms):
        value = 0
        for i, x in enumerate(values): value ^= MUL[form >> (2*i) & 3][x]
        mask |= bool(value) << row
    return mask


class ShearExactThreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.toy = cats.build_catalog(3, 4)

    def test_universal_80_forms_and_high16bits_match_exhaustive_GF4(self):
        forms = cats.universal_forms(3)
        self.assertEqual(len(forms), 80)
        self.assertEqual(len(set(forms)), 80)
        expected = Counter()
        for values in product(range(4), repeat=6):
            pattern = ((1 << 80)-1) ^ output_mask(forms, values)
            expected[sum(bool(x) for x in values), pattern] += 1
        actual = {(int(w), int(lo) | (int(hi) << 64)): int(n) for w, lo, hi, n in
                  zip(self.toy.input_weights, self.toy.zero_low, self.toy.zero_high, self.toy.counts)}
        self.assertEqual(actual, expected)
        self.assertEqual(int(self.toy.zero_low[0]), (1 << 64)-1)
        self.assertEqual(int(self.toy.zero_high[0]), 65535)
        self.assertFalse(self.toy.zero_high.flags.writeable)

    def test_actual_joint_a3_pair_matches_GF4_brute_and_has_mass4096(self):
        for source, column in (((0, 1, 2), 0), ((0, 4, 8), 10)):
            forms = sum(cats.shear.paired_rows(source, column, 1), ())
            expected = Counter()
            for values in product(range(4), repeat=6):
                expected[sum(bool(x) for x in values), output_mask(forms, values)] += 1
            actual = cats.pair_counts(source, column, self.toy)
            self.assertEqual(actual, expected)
            self.assertEqual(sum(actual.values()), 4**6)
            self.assertEqual(actual[0, 0], 1)
            self.assertTrue(any(mask >> 16 and mask & 65535 for _, mask in actual))

    def test_raw_a1_a2_GF256_categories_replay_frozen_counts(self):
        for source in ((5,), (0, 7)):
            catalog = cats.build_catalog(len(source), 256)
            for column in (0, 6, 14):
                old = gate.frozen.pair_data(source, column)
                expected = Counter()
                for w, lo, hi, n in zip(old['input_weights'], old['first'], old['second'], old['raw_counts']):
                    if n: expected[w, lo | (hi << 16)] += n
                self.assertEqual(cats.pair_counts(source, column, catalog), expected)

    def test_tiny_noncommuting_a1_a2_tail_replay_frozen_engine(self):
        zero = np.array(((.7, .2), (.1, .9)))
        nonzero = np.array(((.25, .35), (.4, .2)))
        table = mask_table_from_entries(np.log2(zero), np.log2(nonzero))
        for source in ((5,), (0, 7)):
            a, r = len(source), 1
            catalog = cats.build_catalog(a, 4)
            new_columns, old_blocks = [], []
            for j in range(0, 16, 2):
                weights, masks, counts = cats.pair_arrays(source, j, catalog)
                new_columns.append(gate.pair_coefficients(weights, masks, counts, table, a, groups=r))
                old = gate.frozen.pair_data(source, j, q=4)
                old_blocks.append(gate.frozen.branch_log_matrices(np.array(old['first']), np.array(old['second']), table))
            old_counts = np.array([log2(n) if n else -np.inf for n in old['raw_counts']])
            expected = gate.frozen.word_tail(old_blocks, old_counts, np.array(old['input_weights']), r, 9*a*r)
            observed = gate.tail_backend.positive_deficit_tail(np.stack(new_columns), r, 9*a*r)
            self.assertAlmostEqual(float(observed), float(expected), places=10)

    def test_mock_selected_scope_and_RS_normalization(self):
        saved = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        fake_catalog = SimpleNamespace(arity=3, q=256, counts=np.array([1]), positive_mass_kernels=1)
        with patch.object(gate.frozen, 'load_measure', return_value=(saved, {}, 17.)), \
             patch.object(gate.frozen, 'mask_log_matrices', return_value=None), \
             patch.object(gate.catalog_api, 'build_catalog', return_value=fake_catalog), \
             patch.object(gate, 'prepare_source', return_value=(np.zeros((8, 7, 2, 2)), {})), \
             patch.object(gate.tail_backend, 'positive_deficit_tail', return_value=18*(24*log2(255)+2)) as tail:
            data = gate.run()
        self.assertEqual(tail.call_args.args[1:], (18, 486))
        self.assertTrue(data['complete'])
        self.assertEqual(data['completed_sources'], 3)
        self.assertAlmostEqual(data['selected_rowset_sum_bound_log2'], 17+18*(2+log2(3)))
        self.assertFalse(data['complete_all560_original_a3'])
        self.assertFalse(data['independent_column_support_assumption'])
        self.assertEqual(data['regional_provider_runs'], 0)
        self.assertEqual(data['rows'][0]['RS_factor_per_group'], '255^(-24)')

    def test_frozen_replay_validator_and_exclusive_catalog(self):
        saved, rows = gate.load_replay()
        self.assertEqual(len(rows), 136)
        self.assertEqual(saved['all_labelled_arities'], [1, 2])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'toy.npz'
            report = cats.save_catalog(self.toy, path)
            loaded = cats.load_catalog(path, report['sha256'])
            self.assertTrue(np.array_equal(loaded.zero_high, self.toy.zero_high))
            self.assertTrue(np.array_equal(loaded.counts, self.toy.counts))
            with self.assertRaises(FileExistsError): cats.save_catalog(self.toy, path)

    def test_guards_keep_catalog_small_and_forms_exact(self):
        with self.assertRaises(ValueError): cats.build_catalog(4, 4)
        with self.assertRaises(ValueError): cats.joint_arrays((0,), self.toy)
        with self.assertRaises(ValueError): cats.pair_arrays((0, 1, 2), 1, self.toy)
        with self.assertRaises(ValueError): gate.run(sources=((0, 1, 2),)*2)
        with self.assertRaises(ValueError): gate.run(sources=tuple((0, 1, i) for i in range(2, 16))+( (0, 2, 3), (0, 2, 4), (0, 2, 5)))


if __name__ == '__main__':
    unittest.main()
