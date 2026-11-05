"""Small catalogs and GF4 checks; never build the production GF256/a6 catalog."""
from collections import Counter
from itertools import product
from pathlib import Path
import tempfile
import unittest

import numpy as np

import hd_lambda_double_categories as target


class DenseCategoryTests(unittest.TestCase):
    def test_small_catalogs_match_existing_exact_relations(self):
        for arity in range(1, 5):
            catalog = target.build_catalog(arity, 256)
            forms = target.dense_forms(arity)
            selected = tuple(forms[(i*7+3) % len(forms)] for i in range(16))
            expected = target.exact.base.joint_counts(selected, arity)
            self.assertEqual(target.joint_counts(selected, catalog), expected)
            self.assertFalse(catalog.zero_patterns.flags.writeable)
            self.assertFalse(catalog.counts.flags.writeable)
            self.assertEqual(sum(map(int, catalog.counts)), 256**arity)

    def test_a6_GF4_complete_zero_patterns_and_uint64_top_bit(self):
        catalog = target.build_catalog(6, 4)
        forms = target.dense_forms(6)
        expected = Counter()
        mul = target.exact.base.multiply_gf4
        for values in product(range(4), repeat=6):
            pattern = 0
            for index, form in enumerate(forms):
                value = 0
                for i, x in enumerate(values): value ^= mul(form >> (2*i) & 3, x)
                pattern |= int(value == 0) << index
            expected[sum(bool(x) for x in values), pattern] += 1
        actual = {(int(w), int(p)): int(n) for w, p, n in zip(catalog.input_weights, catalog.zero_patterns, catalog.counts)}
        self.assertEqual(actual, expected)
        self.assertEqual(int(catalog.zero_patterns[0]), (1 << 64)-1)
        self.assertEqual(int(catalog.counts[0]), 1)

    def test_scheduled_columns_match_frozen_provider_and_source_reordering(self):
        for source, q in (((0, 3, 7), 256), (tuple(range(6)), 4)):
            catalog = target.build_catalog(len(source), q)
            cache = target.MappedHistogramCache(catalog, max_entries=16, max_terms=100000)
            observed, info = cache.scheduled_columns(source, columns=(0, 5, 15))
            expected, _ = target.exact.exact_columns(source, q=q, columns=(0, 5, 15))
            self.assertEqual(observed, expected)
            self.assertEqual([list(column.items()) for column in observed],
                             [list(column.items()) for column in expected])
            reversed_columns, _ = cache.scheduled_columns(tuple(reversed(source)), columns=(0, 5, 15))
            self.assertEqual(reversed_columns, expected)
            self.assertEqual(info['original_source'], source)
            self.assertTrue(info['total_input_weight_preserved'])
            self.assertFalse(info['individual_input_row_identity_preserved'])

    def test_cache_is_bounded_and_dictionaries_are_fresh(self):
        catalog = target.build_catalog(2, 4)
        cache = target.MappedHistogramCache(catalog, max_entries=1, max_terms=1024)
        first = cache.mapped_counts((0, 1))
        first[(0, 0)] = 99
        self.assertEqual(cache.mapped_counts((1, 0))[(0, 0)], 1)
        self.assertEqual(cache.hits, 1)
        cache.mapped_counts((1, 7))
        self.assertEqual(len(cache._cache), 1)
        self.assertLessEqual(cache._terms, 1024)
        for indexes, counts in cache._cache.values():
            self.assertFalse(indexes.flags.writeable)
            self.assertFalse(counts.flags.writeable)

    def test_dense_form_index_and_rejection(self):
        for arity in range(1, 7):
            for index, form in enumerate(target.dense_forms(arity)):
                self.assertEqual(target.dense_form_index(form, arity), index)
        for form in (0, 1, 4**6):
            with self.assertRaises(ValueError): target.dense_form_index(form, 6)
        with self.assertRaises(ValueError): target.build_catalog(7, 4)

    def test_exclusive_pinned_small_catalog_roundtrip(self):
        catalog = target.build_catalog(2, 4)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'toy.npz'
            report = target.save_catalog(catalog, path)
            read = target.load_catalog(path, report['sha256'])
            self.assertTrue(np.array_equal(read.counts, catalog.counts))
            self.assertTrue(np.array_equal(read.zero_patterns, catalog.zero_patterns))
            with self.assertRaises(FileExistsError): target.save_catalog(catalog, path)
            with self.assertRaises(ValueError): target.load_catalog(path, '0'*64)


if __name__ == '__main__':
    unittest.main()
