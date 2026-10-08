"""Independent coefficient counts, scaling identities, and directed receipts."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
from unittest.mock import patch
import unittest

from hd_randomized_outer_proof import Parameters
from hd_randomized_outer_joint_quiet import joint_moment_enclosures
import hd_randomized_outer_quiet_exact as exact


class QuietExactTests(unittest.TestCase):
    def setUp(self):
        from flint import ctx
        self.previous = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        from flint import ctx
        ctx.prec = self.previous

    def assert_contains_fraction(self, value, expected):
        from flint import arb
        self.assertTrue((value-arb(expected.numerator)/expected.denominator).contains(0))

    def test_coefficients_match_exhaustive_fixed_subsets_and_free_last(self):
        quiet, free = (F(1), F(1, 4), F(1, 4)), (F(1),)*3
        width, packets, rows = 2, 3, 6
        for alpha in (F(1), F(1, 2)):
            roots = quiet if alpha == 1 else (F(1), F(1, 2), F(1, 2))
            for active in range(rows+1):
                quiet_sum = last_sum = F(0)
                for support in combinations(range(rows), active):
                    counts = tuple(sum(x//width == packet for x in support) for packet in range(packets))
                    weight = roots[counts[0]]*roots[counts[1]]
                    last_sum += weight
                    quiet_sum += weight*roots[counts[2]]
                quiet_sum /= comb(rows, active)
                last_sum /= comb(rows, active)
                outputs = []
                for marker in (F(1, 3), F(4, 5)):
                    result = exact.regional_coefficients(quiet, free, alpha, marker, packets, active)
                    self.assert_contains_fraction(result['quiet'], quiet_sum)
                    self.assert_contains_fraction(result['last'], last_sum)
                    self.assertTrue(result['quiet_probability'] > 0)
                    self.assertTrue(result['quiet_probability'] < 1)
                    self.assertTrue(result['last_probability'] < 1)
                    outputs.append(result)
                self.assertTrue((outputs[0]['quiet']-outputs[1]['quiet']).contains(0))
                self.assertTrue((outputs[0]['last']-outputs[1]['last']).contains(0))

    def test_actual_z_one_coefficients_by_direct_packet_count_sum(self):
        q, alpha = 65536, F(1, 2)
        quiet, free = joint_moment_enclosures(F(1), q)
        for active in range(7):
            expected_quiet = expected_last = F(0)
            for j in range(active+1):
                coefficient = F(comb(64, j)*comb(64, active-j), comb(128, active))
                expected_quiet += coefficient/F(256)**(min(j, 4)+min(active-j, 4))
                expected_last += coefficient/F(256)**min(j, 4)
            result = exact.regional_coefficients(quiet, free, alpha, F(1, 3), 2, active)
            self.assert_contains_fraction(result['quiet'], expected_quiet)
            self.assert_contains_fraction(result['last'], expected_last)

    def test_optimized_marker_matches_closed_form_binomial_tilt(self):
        p = Parameters(field=16, rows=128, outer_n=4, outer_k=2, limit=24)
        local = tuple(-4*j for j in range(65))
        with patch.object(exact, 'moment_logs', return_value=(local, local)):
            witness = exact.optimized_marker(17, F(1, 2), F(1, 2), p)
        self.assertAlmostEqual(witness['log_odds'], log2(F(68, 111)), places=12)
        self.assertAlmostEqual(float(F(witness['marker'])), float(F(68, 179)), places=13)

    def test_marker_gain_is_exactly_the_missing_coefficient_masses(self):
        from flint import arb
        quiet, free = (F(1), F(1, 4), F(1, 4)), (F(1),)*3
        p = exact.regional_coefficients(quiet, free, F(1, 2), F(1, 3), 3, 2)
        regions = 4
        actual = p['quiet']**(regions-1)*p['last']
        marker = p['ordinary']**(3*regions-1)*p['terminal']/p['binomial_mass']**regions
        ratio = p['quiet_probability']**(regions-1)*p['last_probability']
        self.assertTrue((actual-marker*ratio).contains(0))
        self.assertTrue(actual < marker)
        self.assertTrue((marker/actual-1/ratio).contains(0))

    def test_rank_refinement_gap_at_most_two_for_all_residual_profiles(self):
        largest_gap = 0
        profiles = 0
        for counts in product(range(9), repeat=5):
            if sum(counts) != 8:
                continue
            dimensions = tuple(d for d, number in enumerate(counts) for unused in range(number))
            positive_tiles = sum(d > 0 for d in dimensions)
            lower_rank = max(max(dimensions), min(2, positive_tiles))
            upper_rank = min(4, sum(dimensions))
            gap = upper_rank-lower_rank
            self.assertGreaterEqual(gap, 0)
            self.assertLessEqual(gap, 2)
            if lower_rank <= 1:
                self.assertEqual(gap, 0)
            largest_gap = max(largest_gap, gap)
            profiles += 1
        self.assertEqual(profiles, comb(12, 4))
        self.assertEqual(largest_gap, 2)
        receipt = dict(active_rows=1536, mark='1/4', alpha='1/32',
                       exact_coefficient_log2_upper='33495')
        from flint import arb
        floor = exact.rank_refinement_floor(receipt)
        self.assertTrue((arb(floor['maximum_rank_only_gain_bits'])-32767).contains(0))
        self.assertTrue(arb(floor['best_rank_only_log2_lower']) > 727)
        self.assertTrue(floor['same_witness_cannot_certify'])
        receipt['mark'] = '1'
        floor = exact.rank_refinement_floor(receipt)
        self.assertEqual(arb(floor['maximum_rank_only_gain_bits']), 0)

    def test_nested_subset_monotonicity_and_prefactor_envelope(self):
        # Uniformly adding a missing position couples consecutive occupancies.
        roots = (F(1), F(1, 2), F(1, 2))
        values = []
        for active in range(7):
            ordinary = terminal = F(0)
            for support in combinations(range(6), active):
                counts = tuple(sum(i//2 == tile for i in support) for tile in range(3))
                partial = roots[counts[0]]*roots[counts[1]]
                ordinary += partial*roots[counts[2]]
                terminal += partial
            ordinary /= comb(6, active)
            terminal /= comb(6, active)
            values.append(ordinary*terminal)
        self.assertTrue(all(a >= b for a, b in zip(values, values[1:])))
        first = 1
        for last in (1, 2, 3):
            actual = comb(6, last)*2**last*values[last]
            envelope = comb(6, last)*2**last*values[first]
            self.assertLessEqual(actual, envelope)

    def test_sparse_whole_interval_receipt(self):
        from flint import arb
        receipt = exact.directed_sparse_extension()
        self.assertEqual((receipt['active_rows_first'], receipt['active_rows_last']), (954, 962))
        self.assertTrue(arb(receipt['union_log2_upper']) < -226)
        self.assertTrue(receipt['below_2_to_minus_40'])
        self.assertEqual(receipt['exact_regional_occupancies_evaluated'], 1)

    def test_dense_whole_interval_receipt_and_boundary(self):
        from flint import arb
        receipt = exact.directed_dense_extension()
        self.assertEqual((receipt['active_rows_first'], receipt['active_rows_last']), (5865, 8192))
        self.assertTrue(arb(receipt['union_log2_upper']) < -355)
        self.assertTrue(receipt['below_2_to_minus_40'])
        self.assertEqual(tuple(point['active_rows'] for point in receipt['exact_points']),
                         (5865, 5866, 5867, 5868))
        self.assertEqual(receipt['marker_interval']['active_rows_first'], 5869)
        self.assertEqual(receipt['marker_interval']['occupancy_endpoints_evaluated'], 2)
        previous = exact.directed_endpoint(5864, F(1, 16384), F(3, 8))
        self.assertTrue(arb(previous['exact_coefficient_log2_upper']) > 251)
        self.assertFalse(previous['below_2_to_minus_40'])


if __name__ == '__main__':
    unittest.main()
