"""Independent exact checks for sparse swap transfer and scope."""
from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

import hd_shifted_swap_sparse as s
from hd_shifted_general_screen import general_coefficients
from hd_shifted_good_pair_screen import good_pair_coefficients
from hd_shifted_pairs_count import exact_profiles, profile_count_upper


def reference_product(left, right, maximum):
    result = {}
    for (a, b), u in left.items():
        for (c, d), v in right.items():
            if a+c+2*(b+d)>maximum:
                continue
            dst = result.setdefault((a+c, b+d), [0]*4)
            for i in range(2):
                for j in range(2):
                    dst[2*i+j] += u[2*i]*v[j]+u[2*i+1]*v[2+j]
    return result


class SparseSwapTests(unittest.TestCase):
    def test_geometry_counts_before_and_after_fair_coins(self):
        for maximum in (1, 2, 4, 8, 16):
            counts = {}
            for a, b, count, occupancy in s.sparse_geometry(maximum):
                self.assertEqual(sum(occupancy), a+2*b)
                counts[a, b] = counts.get((a, b), 0)+count
            for a, b in s.valid_degrees(maximum):
                expected = comb(32, a)*comb(32-a, b)*2**a
                self.assertEqual(counts[a, b], expected)

    def test_sparse_coefficients_match_independent_oriented_sum(self):
        z, alpha, maximum = .73, .37, 6
        full = general_coefficients(z, alpha, (0., 0., 0.))
        sparse = s.swap_single_screen(z, alpha, maximum)
        for a, b in s.valid_degrees(maximum):
            expected = sum(full[:, u, a-u, b] for u in range(a+1))/2**a
            np.testing.assert_allclose(sparse[:, a, b], expected, rtol=3e-13, atol=0.)

    def test_good_pair_truncation_and_unit_exponent(self):
        for alpha in (.37, 1.):
            z, maximum = .73, 6
            full = good_pair_coefficients(z, alpha).transpose(2, 3, 0, 1).reshape(4, 3, 65)
            actual = s.good_pair_screen(z, alpha, maximum)
            np.testing.assert_allclose(actual, s.truncate_screen(full, maximum), rtol=3e-13, atol=0.)
            if alpha==1:
                one = s.swap_single_screen(z, alpha, maximum).copy()
                one[:, 2:, :] = 0.
                np.testing.assert_allclose(actual, s.product_screen(one, one, maximum), rtol=3e-13, atol=0.)

    def test_carry_free_noncommuting_product(self):
        maximum = 4
        left = {(0, 0): (1, 2, 3, 4), (4, 0): (4, 1, 2, 3), (0, 1): (3, 2, 0, 1)}
        right = {(0, 0): (5, 6, 7, 8), (2, 0): (2, 4, 1, 3), (0, 2): (3, 1, 0, 2)}
        expected = reference_product(left, right, maximum)
        self.assertNotEqual(expected, reference_product(right, left, maximum))
        with ctx.workprec(256):
            actual = s.product_arb(s.encode_arb(left, maximum), s.encode_arb(right, maximum), maximum)
            for a, b in s.valid_degrees(maximum):
                for i in range(4):
                    self.assertEqual(actual[i][a+(2*maximum+1)*b], expected.get((a, b), (0,)*4)[i])
            overflow = s.product_arb(s.encode_arb({(4, 0): (1, 0, 0, 1)}, maximum),
                                     s.encode_arb({(2, 0): (1, 0, 0, 1)}, maximum), maximum)
            self.assertTrue(all(len(poly)==0 for poly in overflow))

    def test_multinomial_normalizer_and_marker_cancellation(self):
        maximum, marker = 6, F(3, 4096)
        coefficients = {(a, b): (comb(32, a)*comb(32-a, b), 0, 0, comb(32, a)*comb(32-a, b))
                        for a, b in s.valid_degrees(maximum)}
        with ctx.workprec(256):
            run = s.power_arb(s.encode_arb(coefficients, maximum, marker), 128, maximum)
            for a, b in s.valid_degrees(maximum):
                expected = comb(4096, a)*comb(4096-a, b)*s.rational(marker)**(a+b)
                for i in (0, 3):
                    ratio = run[i][a+(2*maximum+1)*b]/expected
                    self.assertLess(abs(ratio-1), arb(2)**-200)
                for i in (1, 2):
                    self.assertEqual(len(run[i]), 0)

    def test_indexed_support_counts_and_nonempty_scope(self):
        for pairs in range(2, 7):
            counts = exact_profiles(pairs)
            for (r, a, b), count in counts.items():
                self.assertGreaterEqual(profile_count_upper(pairs, r, a, b), count)
            for r in range(1, 2*pairs):
                self.assertEqual(sum(count for (weight, _, _), count in counts.items() if weight==r), comb(2*pairs, r))
        for bad in (0, 65, -1):
            with self.assertRaises(ValueError):
                s.weight_receipt(bad, F(1199, 1200), F(1))

    def test_directed_precision_replay(self):
        receipts = [s.weight_receipt(1, F(1199, 1200), F(1), precision=p) for p in (256, 384)]
        self.assertEqual(receipts[0]['log2_upper_rational'], receipts[1]['log2_upper_rational'])
        self.assertLess(F(receipts[0]['log2_upper_rational']), -107)
        self.assertTrue(receipts[0]['covers_all_supports_at_weight'])
        self.assertTrue(receipts[0]['covers_all_even_odd_splits'])
        self.assertEqual(receipts[0]['construction'], s.SWAP_CONSTRUCTION)
        self.assertEqual(receipts[0]['ranges'], ((1, 1), (0, 4096), (0, 4096)))
        self.assertEqual(receipts[0]['integer_profile_triples'], 1)
        self.assertFalse(receipts[0]['whole_code_certificate'])

    def test_fixed_family_scope_and_union(self):
        expected = {1: (F(1199, 1200), F(1)), 8: (F(149, 150), F(1)),
                    9: (F(397, 400), F(2, 5)), 17: (F(1183, 1200), F(2, 5)),
                    18: (F(691, 700), F(1, 5)), 32: (F(171, 175), F(1, 5)),
                    33: (F(1367, 1400), F(1, 10)), 64: (F(167, 175), F(1, 10))}
        for r, witness in expected.items():
            self.assertEqual(s.family_witness(r), witness)
        result = s.family_receipt(2)
        self.assertEqual(result['integer_profile_triples'], 4)
        self.assertEqual(result['log2_upper_rational'], '-107393558356229/1000000000000')
        self.assertEqual(result['ranges'], ((1, 2), (0, 4096), (0, 4096)))
        self.assertTrue(result['pair_swaps'])
        self.assertFalse(result['whole_code_certificate'])

    def test_batch_reuses_regions_and_matches_single_weights(self):
        z, alpha = F(1199, 1200), F(2, 5)
        for precision in (256, 384):
            with patch.object(s, 'regional_directed', wraps=s.regional_directed) as regional:
                result = s.batch_receipt(1, 3, z, alpha, maximum=6, precision=precision)
            regional.assert_called_once_with(z, alpha, 6, 60)
            for row in result['weight_receipts']:
                single = s.weight_receipt(row['active_groups'], z, alpha, precision=precision)
                self.assertEqual(row['log2_upper_rational'], single['log2_upper_rational'])
                self.assertEqual(row['profiles'], single['profiles'])
                self.assertEqual(row['integer_profile_triples'], single['integer_profile_triples'])
                self.assertEqual(row['witness']['coefficient_maximum'], 6)
                self.assertEqual(row['witness']['coefficient_marker'], '3/2048')
            self.assertEqual(result['ranges'], ((1, 3), (0, 4096), (0, 4096)))
            self.assertEqual(result['integer_profile_triples'],
                             sum(row['integer_profile_triples'] for row in result['weight_receipts']))
            self.assertTrue(result['covers_all_supports_in_range'])
            self.assertTrue(result['covers_all_even_odd_splits'])
            self.assertEqual(result['construction'], s.SWAP_CONSTRUCTION)
            self.assertTrue(result['pair_swaps'])
            self.assertFalse(result['whole_code_certificate'])

    def test_batch_scope_excludes_unrequested_weights(self):
        result = s.batch_receipt(2, 3, F(1199, 1200), F(1), maximum=6)
        self.assertEqual([row['active_groups'] for row in result['weight_receipts']], [2, 3])
        self.assertEqual(result['ranges'], ((2, 3), (0, 4096), (0, 4096)))
        self.assertEqual(result['integer_profile_triples'], s.profile_count(result['ranges']))
        for first, last, maximum in ((0, 3, 6), (3, 2, 6), (2, 7, 6), (1, 257, 257), (1, 3, 0)):
            with self.assertRaises(ValueError):
                s.batch_receipt(first, last, F(1199, 1200), F(1), maximum=maximum)
        for call in (s.family_witness, s.family_receipt):
            with self.assertRaises(ValueError):
                call(65)

    def test_fixed_batch_schedule_is_disjoint_and_complete(self):
        self.assertEqual(s.BATCH_SCHEDULE,
                         ((65, 72, F(333, 350), F(1, 10)),
                          (73, 90, F(659, 700), F(1, 10)),
                          (91, 110, F(13, 14), F(1, 10)),
                          (111, 128, F(32, 35), F(1, 10))))
        weights = [r for first, last, _, _ in s.BATCH_SCHEDULE for r in range(first, last+1)]
        self.assertEqual(weights, list(range(65, 129)))
        self.assertEqual(sum(s.profile_count(((first, last), (0, 4096), (0, 4096)))
                             for first, last, _, _ in s.BATCH_SCHEDULE), 159104)

    def test_cutoff_256_degree_and_receipt_scope(self):
        maximum = 256
        degrees = s.valid_degrees(maximum)
        indices = [a+(2*maximum+1)*b for a, b in degrees]
        self.assertEqual(len(indices), 129**2)
        self.assertEqual(len(set(indices)), len(indices))
        self.assertEqual(max(indices), 65664)
        self.assertLess(2*maximum, 2*maximum+1)
        # Check the public boundary without replaying the large polynomial.
        # The all-profile helper receives the requested r, not the cutoff.
        with patch.object(s, 'regional_directed', return_value={(1, 0): (arb(1), arb(0), arb(0), arb(1))}) as regional:
            result = s.batch_receipt(1, 1, F(1199, 1200), F(1), maximum=maximum)
        regional.assert_called_once_with(F(1199, 1200), F(1), maximum, 60)
        self.assertEqual(result['ranges'], ((1, 1), (0, 4096), (0, 4096)))
        self.assertEqual(result['witness']['coefficient_marker'], '1/16')
        self.assertEqual(result['integer_profile_triples'], 1)
        self.assertFalse(result['whole_code_certificate'])

    def test_repaired_256_schedule_is_disjoint_and_complete(self):
        self.assertEqual(s.BATCH_SCHEDULE_256,
                         ((129, 160, F(157, 175), F(1, 10)),
                          (161, 192, F(153, 175), F(1, 10)),
                          (193, 225, F(149, 175), F(1, 10)),
                          (226, 256, F(29, 35), F(1, 10))))
        weights = [r for first, last, _, _ in s.BATCH_SCHEDULE_256 for r in range(first, last+1)]
        self.assertEqual(weights, list(range(129, 257)))
        self.assertEqual([s.profile_count(((first, last), (0, 4096), (0, 4096)))
                          for first, last, _, _ in s.BATCH_SCHEDULE_256],
                         [171200, 254144, 366241, 456415])
        self.assertEqual(s.profile_count(((1, 256), (0, 4096), (0, 4096))), 1431040)


if __name__ == '__main__':
    unittest.main()
