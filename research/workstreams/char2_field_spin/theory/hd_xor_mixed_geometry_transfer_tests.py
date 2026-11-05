"""Small exact tests for the separate, geometry-independent count transfer."""
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_mixed_geometry_transfer import BASE_CEILING, transfer


class GeometryTransferTests(unittest.TestCase):
    def test_arbitrary_labelled_pairs_satisfy_resource_constraints(self):
        # Every pair assignment for two labelled groups, with all empty,
        # singleton, uncancelled-double, and cancelled-double column states.
        pairs = tuple(combinations(range(8), 2))
        for edges in product(pairs, repeat=2):
            for states in product(range(5), repeat=2):
                occupancy = [0]*8
                s = d = k = 0
                for (i, j), state in zip(edges, states):
                    if state == 0:
                        mask = 0
                    elif state in (1, 2):
                        s += 1
                        mask = 255 ^ (1 << (i if state == 1 else j))
                    elif state == 3:
                        d += 1
                        mask = 255
                    else:
                        d += 1
                        k += 1
                        mask = (1 << i) | (1 << j)
                    for region in range(8):
                        occupancy[region] += (mask >> region) & 1
                self.assertEqual(sum(occupancy), 7*s+8*d-6*k)
                self.assertTrue(all(d-k <= m <= s+d for m in occupancy))
                self.assertTrue(all(occupancy[i]+occupancy[j] <= 2*(s+d)
                                    for i, j in pairs))

    def test_independent_higher_precision_arithmetic_and_scope(self):
        expected = F('-7550725335987/200000000000')
        with ctx.workprec(1024):
            base = arb(BASE_CEILING.numerator)/BASE_CEILING.denominator
            self.assertEqual(rational_ceiling(base+15*arb(28).log()/arb(2).log()), expected)
        result = transfer()
        self.assertEqual(F(result['failure_log2_upper']), expected)
        self.assertEqual(result['exact_mass_multiplier'], 28**15)
        self.assertEqual(result['previous_sparse_pair_assignments'], 28)
        self.assertEqual(result['derived_sparse_pair_assignments'], 28**16)
        self.assertEqual(result['source_matching_counts'], [8, 0])
        self.assertEqual(result['cohort_counts'], [0, 16, 1])
        self.assertEqual(result['operator_replays'], 0)
        self.assertFalse(result['reaches_40_bits'])
        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
