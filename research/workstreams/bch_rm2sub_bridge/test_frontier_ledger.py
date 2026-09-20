"""Coverage regressions, including equality with the frozen K24 ledger."""
from fractions import Fraction as F
import unittest
import bridge as base
import frontier_ledger as ledger
from audit_frontier_k24 import FILES


class LedgerTest(unittest.TestCase):
    def test_exact_interval_partition(self):
        ledger.validate_coverage([(5,8),(1,1),(2,4)],8)
        for intervals in ([(1,4),(6,8)],[(1,4),(4,8)],[(1,7)],[(2,8)],[(1,9)]):
            with self.assertRaises(AssertionError):ledger.validate_coverage(intervals,8)

    def test_complete_density_and_integer_tree(self):
        self.assertEqual(ledger.validate_tree(['q','p',-12,-13,-14],2,9),3)
        for tree,lo,hi in ((['q',-1],2,9),([-1,-2],2,9),(['q',-1,-2],2,2),([True],2,9)):
            with self.assertRaises(AssertionError):ledger.validate_tree(tree,lo,hi)

    def test_matches_frozen_per_occupancy_ledger(self):
        expected=base.read(base.HERE/'generated/frontier_k24_full_v1.json')
        got=ledger.build(24,[base.HERE/'generated'/name for name in FILES])
        for key in ('failure_upper','higher_occupancy_upper','cutoff','covered_occupancy_range'):
            self.assertEqual(got[key],expected[key])
        self.assertTrue(got['full_distance_proved'])
        self.assertEqual(got['strict_certified_bits'],46)
        self.assertLess(base.decode(got['failure_upper']),F(1,1<<46))


if __name__=='__main__':unittest.main()
