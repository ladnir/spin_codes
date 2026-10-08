"""Exact all-route occupancy averages and integer quiet-subcode caps."""

from fractions import Fraction
from itertools import combinations
import unittest

from hd_blockwise_shared_budget import covered_mean, tiled_quiet_cap


class SharedBudgetTests(unittest.TestCase):
    def test_exhaustive_row_set_average_for_fixed_routes(self):
        for rows in (4, 6, 8):
            for tile in (g for g in (1, 2, 4) if rows % g == 0):
                # Deliberately correlated, nonrandom regional permutations.
                routes = (tuple(range(rows)), tuple(reversed(range(rows))),
                          tuple((i+1) % rows for i in range(rows)))
                for active in range(rows+1):
                    totals = []
                    for subset in combinations(range(rows), active):
                        selected = set(subset)
                        totals.append(sum(tile for region in routes
                                          for start in range(0, rows, tile)
                                          if selected.intersection(region[start:start+tile])))
                    self.assertEqual(Fraction(sum(totals), len(totals)),
                                     covered_mean(rows, active, tile, 3*rows))

    def test_concavity_makes_endpoint_gate_sufficient(self):
        for rows, tile in ((8, 2), (8, 4), (16, 8)):
            values = [covered_mean(rows, a, tile, 4*rows)-2*a
                      for a in range(rows+1)]
            for a in range(1, rows):
                self.assertLessEqual(values[a+1]-values[a], values[a]-values[a-1])
            for first in range(rows+1):
                self.assertEqual(min(values[first:]), min(values[first], values[-1]))

    def test_k20_structural_gates(self):
        self.assertTrue(tiled_quiet_cap(state=1)['rejects_40_percent'])
        self.assertTrue(tiled_quiet_cap(width=256)['rejects_40_percent'])
        for options in (dict(state=2), dict(tile=4), {}, dict(width=128)):
            result = tiled_quiet_cap(**options)
            self.assertFalse(result['rejects_40_percent'])
            self.assertFalse(result['whole_certificate'])
            self.assertGreater(result['quiet_dimension'], 0)
        two = tiled_quiet_cap(state=2)
        self.assertLess(Fraction(two['distance_upper']), Fraction(41, 100))
        wide = tiled_quiet_cap(width=128)
        self.assertEqual(two['selected_rows'], wide['selected_rows'])
        self.assertEqual(two['quiet_weight_cap']-wide['quiet_weight_cap'], 2)

    def test_tile_one_recovers_identity_gate(self):
        result = tiled_quiet_cap(tile=1)
        self.assertEqual(result['selected_rows'], 1024)
        self.assertEqual(result['quiet_dimension'], 4)
        self.assertEqual(result['quiet_weight_cap'], 1024*256-3)

    def test_one_packet_and_unavailable_quiet_witness(self):
        result = tiled_quiet_cap(8, 4, 2, 16, 8, 2)
        self.assertEqual(result['quiet_dimension'], 2)
        self.assertEqual(result['selected_rows'], 1)
        self.assertEqual(result['quiet_weight_cap'], 7)
        self.assertFalse(tiled_quiet_cap(8, 4, 2, 2, 8, 2)['quiet_witness'])
        rate_one=tiled_quiet_cap(8, 4, 4, 8, 1, 2)
        self.assertEqual(rate_one['distance_upper'], '1/8')
        self.assertTrue(rate_one['rejects_40_percent'])

    def test_validation(self):
        for options in (dict(width=63), dict(tile=0), dict(state=-1),
                        dict(message_symbols=7), dict(message_symbols=128*7)):
            with self.assertRaises(ValueError):
                tiled_quiet_cap(**options)


if __name__ == '__main__':
    unittest.main()
