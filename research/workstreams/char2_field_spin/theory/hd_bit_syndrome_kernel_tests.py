"""Tiny exact tests; no production distance contraction or provider."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from gf8_packet_maps import rank
from hd_bit_syndrome_kernel import (
    _fiber_projection_terms, affine_eight_orbits, affine_image_counts,
    build_kernel, evaluate, outside_weight_envelopes,
    evaluate_kernel, exact_sparse_terms, information_cover, information_sets,
)
from hd_bit_syndrome_outer import affine_rank


class BitSyndromeKernelTests(unittest.TestCase):
    def test_affine_evaluation_code_counts_by_gf4_enumeration(self):
        positions = (0, 1, 2, 3)
        expected = Counter()
        # Four-point affine square, dimension3 over GF4.
        for coeff in product(range(4), repeat=3):
            y = [coeff[0] ^ (coeff[1] if i & 1 else 0) ^ (coeff[2] if i & 2 else 0)
                 for i in positions]
            expected[sum(bool(value) << i for i, value in zip(positions, y))] += 1
        self.assertEqual(affine_image_counts(positions, q=4, bits=2), dict(expected))

    def test_exact_translated_fibers_and_projection_dominance(self):
        # Independent eight-coordinate binary H=I+VV^T, over GF4. The
        # selected four-coordinate source has cross rank3 and quiet dim1.
        n, source = 8, (0, 1, 2, 3)
        h = [[int(i == j) ^ ((1 + (i & j).bit_count()) & 1) for j in range(n)] for i in range(n)]
        self.assertEqual(rank(tuple(tuple(h[i][j] for j in source) for i in range(4, 8))), 3)
        quiet, nonquiet, info, _ = _fiber_projection_terms(source, q=4, bits=3)
        projection = sum(1 << i for i in (*range(4, 8), *info))
        expected_quiet, expected_nonquiet = Counter(), Counter()
        actual_words = []
        for values in product(range(4), repeat=4):
            y = []
            for i in range(n):
                value = 0
                for j, x in zip(source, values):
                    if h[i][j]:
                        value ^= x
                y.append(value)
            support = sum(bool(value) << i for i, value in enumerate(y))
            weight = sum(bool(value) for value in values)
            if not support & 0xf0:
                expected_quiet[(weight, support)] += 1
            else:
                expected_nonquiet[(4, support & projection)] += 1
            actual_words.append((weight, support))
        self.assertEqual(quiet, dict(expected_quiet))
        self.assertEqual(nonquiet, dict(expected_nonquiet))
        for t in ((F(1),) * 8, (F(1, 2),) * 8,
                  (F(0), F(1, 3), F(1), F(1, 5), F(1, 7), F(1), F(1, 2), F(0))):
            for x in (F(1), F(3, 2)):
                true_value = evaluate(Counter(actual_words), t, x)
                upper = evaluate(quiet, t, x) + evaluate(nonquiet, t, x)
                self.assertLessEqual(true_value, upper)
        # The original polynomial is not claimed above1. The explicit
        # constrained-coordinate factors safely restore that domain.
        model = dict(quiet_terms=quiet, nonquiet_projection_terms=nonquiet,
                     exact_sparse_terms={}, source=source, information=info)
        markers = (F(2), F(3), F(1, 2), F(4), F(2), F(1, 3), F(1), F(3, 2))
        truth = 0
        for weight, mask in actual_words:
            term = F(3, 2) ** weight
            for i, marker in enumerate(markers):
                if mask >> i & 1:
                    term *= marker
            truth += term
        self.assertLessEqual(truth, evaluate_kernel(model, markers, F(3, 2)))

    def test_information_sets_and_fractional_cover(self):
        sets = information_sets(tuple(range(8)))
        self.assertTrue(all(len(row) == 4 for row in sets))
        self.assertEqual(information_cover(tuple(range(8))), {i: F(1, 2) for i in range(8)})
        # Exact Holder bound for a translated repetition fiber over GF4:
        # (sum_x product_i t_i^[x+offset_i !=0])^4 <= product_i(1+3*t_i^4).
        markers = (F(1, 2), F(1, 3), F(1), F(1, 5))
        for offsets in product(range(4), repeat=4):
            lhs = sum(evaluate({(0, sum(bool(x ^ off) << i for i, off in enumerate(offsets))): 1}, markers)
                      for x in range(4))
            rhs = 1
            for t in markers:
                rhs *= 1 + 3 * t ** 4
            self.assertLessEqual(lhs ** 4, rhs)

    def test_sparse_counts_and_positive_restoration(self):
        for source in (tuple(range(8)), (0, 1, 2, 3, 4, 5, 6, 8)):
            sparse = exact_sparse_terms(source)
            by_degree = Counter()
            for (degree, _), count in sparse.items():
                by_degree[degree] += count
            self.assertEqual(by_degree, {0: 1, 1: 8 * 255, 2: comb(8, 2) * 255 ** 2})
            base = build_kernel(source, exact_sparse=False)
            restored = build_kernel(source)
            self.assertEqual(evaluate(base['terms'], (1,) * 16), 256 ** 8)
            self.assertEqual(evaluate(restored['terms'], (1,) * 16), 256 ** 8)
            self.assertTrue(all(value > 0 for value in restored['terms'].values()))
            for markers in ((F(1, 2),) * 16, tuple(F(i % 5, 5) for i in range(16))):
                for x in (F(1), F(2)):
                    self.assertEqual(evaluate_kernel(restored, markers, x),
                                     evaluate(restored['terms'], markers, x))
                    self.assertLessEqual(evaluate(restored['terms'], markers, x),
                                         evaluate(base['terms'], markers, x))
            self.assertFalse(restored['metadata']['whole_code_certificate'])
            self.assertFalse(restored['metadata']['rs_normalization_included'])

    def test_marker_and_scope_guards(self):
        with self.assertRaises(ValueError):
            build_kernel((0, 1))
        with self.assertRaises(ValueError):
            build_kernel(tuple(range(8)), information=(0,))
        with self.assertRaises(ValueError):
            evaluate({(0, 0): 1}, (F(2),))
        with self.assertRaises(ValueError):
            evaluate({(0, 0): 1}, (F(1),), F(1, 2))
        with self.assertRaises(ValueError):
            affine_image_counts((0, 1), q=3)

    def test_count_orbits_and_coefficient_envelopes(self):
        orbits = affine_eight_orbits()
        self.assertEqual(sum(count for _, count in orbits), comb(16, 8))
        self.assertEqual(sum(count for mask, count in orbits if affine_rank(mask) == 4), 30)
        caps = dict(outside_weight_envelopes())
        self.assertEqual(caps[4][0], 0)
        self.assertEqual(caps[5][0], 0)
        # Independent sample includes labelled sets away from representatives.
        samples = [tuple(range(8)), (1, 2, 4, 7, 8, 9, 12, 15), (0, 2, 4, 6, 8, 10, 12, 14)]
        samples += [tuple(i for i in range(16) if mask >> i & 1) for mask, _ in orbits]
        for positions in samples:
            r = affine_rank(sum(1 << i for i in positions))
            for mask, count in affine_image_counts(positions).items():
                if mask:
                    self.assertLessEqual(count, caps[r][mask.bit_count()])
        # Rank4 outside sets are affine cubes: nonzero evaluation weights
        # are exactly4,6,7,8. These zeros must survive the envelope.
        self.assertEqual(tuple(w for w, count in enumerate(caps[4]) if count), (4, 6, 7, 8))


if __name__ == '__main__':
    unittest.main()
