"""Tiny independent checks; no selected bank or production replay is built."""
from copy import deepcopy
from fractions import Fraction as F
from io import BytesIO
from itertools import combinations
import json
import unittest

import numpy as np
from flint import arb, ctx

import hd_bit_syndrome_directed_recheck as check


def brute_counts(size, weight):
    counts = dict(F0=0, F1=0, F_le1=0, F_ge2=0, all_sources=0)
    for positions in combinations(range(size), weight):
        support = set(positions)
        first = sum(i in support and i+1 in support for i in range(0, size, 2))
        second = sum(i in support and (i+1) % size in support for i in range(1, size, 2))
        value = max(first, second)
        counts['all_sources'] += 1
        counts['F0' if value == 0 else 'F1' if value == 1 else 'F_ge2'] += 1
        counts['F_le1'] += value <= 1
    return counts


def exact_word(matrices):
    left, right = F(1), F(0)
    for a, b, c, d in matrices:
        left, right = left*a+right*c, left*b+right*d
    return left+right


def cache_fixture(**changes):
    arrays = dict(quiet=np.zeros((check.ROWS, 2, 2)), quiet_scale=0,
                  nonquiet=np.zeros((check.ROWS, 2, 2)), nonquiet_scale=0,
                  combined=np.zeros((check.ROWS, 2, 2)), combined_scale=0,
                  words=np.zeros(check.ROWS), words_scale=0,
                  rowsets=np.array(list(combinations(range(16), 8)), dtype=np.int8),
                  mark='993/1000', precision=256)
    arrays.update(changes)
    buffer = BytesIO()
    np.savez(buffer, **arrays)
    return buffer.getvalue()


class IndependentRecheckTests(unittest.TestCase):
    def test_cycle_transfer_against_indexed_exhaustion(self):
        for size in (4, 6, 8, 10):
            for weight in range(size+1):
                self.assertEqual(check.independent_source_counts(size, weight),
                                 brute_counts(size, weight))

    def test_noncommuting_order_and_terminal_vector(self):
        a = tuple(map(F, (1, 2, 0, 3)))
        b = tuple(map(F, (4, 0, 5, 1)))
        sequence = (a, b, a, a, b)
        with ctx.workprec(512):
            value = check.row_word(sequence)
            self.assertEqual(check._endpoint(value), exact_word(sequence))
            self.assertNotEqual(exact_word(sequence), exact_word(tuple(reversed(sequence))))
            # Initial e0 and terminal ones: the off-diagonal entry matters.
            self.assertEqual(check._endpoint(check.row_word((a,))), 3)

    def test_sixteen_products_extreme_scales(self):
        for exponent in (-1200, 1200):
            matrix = tuple(check.exact_scaled(v, exponent) for v in (0.5, 0.25, 0.0, 0.5))
            expected = exact_word((matrix,)*16)
            with ctx.workprec(512):
                value = check.row_word((matrix,)*16)
                self.assertLessEqual(check._endpoint(value, False), expected)
                self.assertGreaterEqual(check._endpoint(value), expected)
                check.assert_dominated(value, check._endpoint(value), 'exact upper')
                with self.assertRaises(ArithmeticError):
                    check.assert_dominated(value, expected/2, 'too small')

    def test_upper_not_lower_comparison(self):
        with ctx.workprec(512):
            uncertain = arb('1 +/- 0.125')
            with self.assertRaises(ArithmeticError):
                check.assert_dominated(uncertain, F(1), 'overlap is insufficient')
            check.assert_dominated(uncertain, F(2), 'upper covered')

    def test_root_sum_and_independent_prefix(self):
        with ctx.workprec(512):
            value = check.minkowski([arb(2)**18, arb(3)**18], 18)
            self.assertLessEqual(check._endpoint(value, False), 5**18)
            self.assertGreaterEqual(check._endpoint(value), 5**18)
            # The exact preimage cancels the independently implemented prefix.
            preimage = check._arb(check.MARK)**check.LIMIT*arb(255)**check.DENSITY_EXPONENT
            result = check.independent_prefix(preimage, 7)
            self.assertLessEqual(check._endpoint(result, False), 7)
            self.assertGreaterEqual(check._endpoint(result), 7)

    def test_scope_identity_and_boolean_types(self):
        receipt = dict(check.EXPECTED_SCOPE, precision_bits=256,
                       source_sha256=check.DRIVER_SHA, reviewed_sources=dict(check.SOURCE_PINS))
        check._checked_scope(receipt)
        for key, bad in (('construction', 'old'), ('directed', 1), ('groups', 17),
                         ('post_refresh_bank_used', True), ('young_weights', ['1/2', '1/2'])):
            broken = deepcopy(receipt)
            broken[key] = bad
            with self.assertRaises(ValueError):
                check._checked_scope(broken)

    def test_exact_labels_and_cache_rejections(self):
        receipt = dict(precision_bits=256)
        arrays, scales = check._cache_arrays(cache_fixture(), receipt)
        self.assertEqual(arrays['combined'].shape, (12870, 2, 2))
        self.assertEqual(scales['words'], 0)
        rows = np.array(list(combinations(range(16), 8)), dtype=np.int8)
        rows[[0, 1]] = rows[[1, 0]]
        for changes in (dict(rowsets=rows), dict(mark='99/100'), dict(precision=384),
                        dict(combined_scale=1.5), dict(words=np.full(12870, np.nan)),
                        dict(words=np.zeros(12870, dtype=np.float32))):
            with self.assertRaises(ValueError):
                check._cache_arrays(cache_fixture(**changes), receipt)

    def test_duplicate_json_and_invalid_inputs(self):
        with self.assertRaises(ValueError):
            json.loads('{"x":1,"x":2}', object_pairs_hook=check._unique_object)
        for value, scale in ((-1.0, 0), (float('inf'), 0), (1.0, True), (1.0, 100001)):
            with self.assertRaises(ValueError):
                check.exact_scaled(value, scale)
        for size, weight in ((3, 1), (8, 9), (8192, 33)):
            with self.assertRaises(ValueError):
                check.independent_source_counts(size, weight)
        with self.assertRaises(ValueError):
            check._sha('0'*63)


if __name__ == '__main__':
    unittest.main()
