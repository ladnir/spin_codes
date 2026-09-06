import unittest

from parameter_cost_frontiers import frontier


class CostFrontierTest(unittest.TestCase):
    def test_tradeoff_preserves_incomparable_candidates(self):
        rows=[dict(name='a',xor_per_bit=10,updates_per_bit=1/64),
              dict(name='b',xor_per_bit=11,updates_per_bit=1/256),
              dict(name='dominated',xor_per_bit=12,updates_per_bit=1/64),
              dict(name='tie',xor_per_bit=10,updates_per_bit=1/64)]
        self.assertEqual({r['name'] for r in frontier(rows)},{'a','b','tie'})


if __name__=='__main__':unittest.main()
