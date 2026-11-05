"""Tiny catalogs/tests only; never construct a GF256 a5/a6 catalog or gate."""
from collections import Counter
from itertools import product
from math import comb
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

import hd_lambda_shear_shared_ratio_extended as extended
import hd_lambda_shear_shared_ratio as frozen
from hd_lambda_scheduled_relations_tests import direct_scheduled


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, dtype=float))


class ExtendedRatioTests(unittest.TestCase):
    def test_low_arity_GF256_replay_against_frozen_single_columns(self):
        for source in ((9,), (0, 11), (0, 4, 13), (0, 1, 4, 8)):
            catalog = extended.build_catalog(len(source), q=256)
            for column in (0, 7, 15):
                self.assertEqual(extended.column_histogram(source, column, catalog),
                                 frozen.column_histogram(source, column))
            self.assertEqual(sum(map(int, catalog.counts)), 256**len(source))

    def test_a5_a6_GF4_exhaustive_actual_Hlambda_columns(self):
        for source, column in (((0, 1, 2, 3, 4), 0), ((0, 1, 2, 3, 4, 5), 9)):
            a = len(source)
            catalog = extended.build_catalog(a, q=4)
            actual = extended.column_histogram(source, column, catalog)
            expected = Counter()
            for values in product(range(4), repeat=a):
                output = direct_scheduled(values, source, frozen.SCHEDULE, column, 2)
                expected[sum(bool(x) for x in values), output] += 1
            self.assertEqual(actual, tuple((w, support, n) for (w, support), n in sorted(expected.items())))
            self.assertEqual(sum(row[2] for row in actual), 4**a)
            for weight in range(a+1):
                self.assertEqual(sum(n for w, _, n in actual if w == weight), comb(a, weight)*3**weight)

    def test_dense_dictionary_64th_bit_and_unique_zero(self):
        catalog = extended.build_catalog(6, q=4)
        self.assertEqual(len(extended.universal_forms(6)), 64)
        self.assertEqual(int(catalog.zero_masks[0]), (1 << 64)-1)
        self.assertEqual(int(catalog.counts[0]), 1)
        self.assertEqual(int(catalog.input_weights[0]), 0)
        self.assertTrue(np.any(catalog.zero_masks[1:] & np.uint64(1 << 63)))
        self.assertFalse(catalog.zero_masks.flags.writeable)

    def test_GF256_a5_a6_chart_mass_without_catalog_build(self):
        for arity in (5, 6):
            plan = extended.catalog_plan(arity, q=256)
            self.assertEqual(plan['exact_raw_mass'], 256**arity)
            self.assertEqual(plan['dense_forms'], 1 << arity)
            self.assertGreater(plan['positive_mass_kernels'], 10000)
        self.assertEqual(extended.catalog_plan(5)['positive_mass_kernels'], 12277)

    def test_pair_coefficients_exact_frozen_replay(self):
        prepared = frozen.prepare_atoms(logs([.2, .8]),
            logs([[[.7, .2], [.1, .6]], [[.2, .5], [.3, .1]]]), 3, rows=2)
        first = ((0, 0, 1), (1, 1, 3), (2, 3, 9), (4, 2, 7))
        second = ((0, 0, 1), (1, 2, 3), (3, 3, 27))
        expected = frozen.pair_coefficients(first, second, prepared, q=4, chunk_size=3)
        actual = extended.pair_coefficients(first, second, prepared, q=4, chunk_size=3)
        np.testing.assert_array_equal(actual, expected)

    def test_weights_through_twelve_and_raw_product_above_uint64(self):
        prepared = frozen.prepare_atoms([0.], logs([[[1., 2.], [3., 1.]]]), 2, rows=2)
        count = 256**6-1
        hist = ((0, 0, 1), (6, 3, count))
        got = extended.pair_coefficients(hist, hist, prepared, chunk_size=1)
        self.assertEqual(got.shape, (13, 2, 2))
        base = frozen.mixture_log_norms(np.array([3]), np.array([3]), prepared)[0]
        np.testing.assert_allclose(got[12], base+2*np.log2(float(count)), rtol=0, atol=1e-13)
        self.assertGreater(count*count, 2**64)
        self.assertTrue(np.all(np.isneginf(got[1:6])))

    def test_catalog_archive_provenance_and_no_overwrite(self):
        catalog = extended.build_catalog(2, q=4)
        with TemporaryDirectory() as folder:
            path = Path(folder)/'tiny.npz'
            receipt = extended.save_catalog(catalog, path)
            restored = extended.load_catalog(path, receipt['sha256'])
            np.testing.assert_array_equal(restored.zero_masks, catalog.zero_masks)
            np.testing.assert_array_equal(restored.counts, catalog.counts)
            self.assertEqual(restored.source_pins, catalog.source_pins)
            with self.assertRaises(FileExistsError):
                extended.save_catalog(catalog, path)
            with self.assertRaises(ValueError):
                extended.load_catalog(path, '0'*64)

    def test_budget_guards_and_scope(self):
        for arity, q in ((7, 256), (0, 4), (5, 16)):
            with self.assertRaises(ValueError):
                extended.catalog_plan(arity, q)
        prepared = frozen.prepare_atoms([0.], np.zeros((1, 2, 2)), 2, rows=2)
        hist = ((0, 0, 1), (6, 3, 3))
        with self.assertRaises(ValueError):
            extended.pair_coefficients(hist, hist, prepared, max_products=3)
        with self.assertRaises(ValueError):
            extended.pair_coefficients(hist, hist, prepared, max_products=1_000_001)
        with self.assertRaises(ValueError):
            extended.pair_coefficients(((7, 3, 1),), hist, prepared)
        with self.assertRaises(ValueError):
            extended.pair_coefficients(hist, hist, prepared, chunk_size=8193)
        metadata = extended.metadata()
        self.assertEqual(metadata['maximum_pair_input_weight'], 12)
        for name in ('old_double_matrix_used', 'whole_code_certificate', 'raw_pair_counts_fit_uint64_claimed'):
            self.assertFalse(metadata[name])


if __name__ == '__main__':
    unittest.main()
