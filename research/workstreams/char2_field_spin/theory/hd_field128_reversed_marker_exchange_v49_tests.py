"""Small normalization and witness-order tests; no production search."""
import unittest
from fractions import Fraction as F
from math import log2
from types import SimpleNamespace

import numpy as np
import hd_field128_reversed_marker_exchange_v49 as subject
import hd_field128_lift_alias_screen as alias


class Tests(unittest.TestCase):
    def test_marker_normalization_and_fallback(self):
        rng = np.random.default_rng(49)
        tables = [SimpleNamespace(logs=rng.normal(size=(4, 128, 128)),
                  reliable=rng.random((4, 128, 128)) > .4,
                  positive_numerator=rng.normal(size=4)) for _ in range(2)]
        odds = [(.13, -.4), (-.2, .8)]
        bank = object.__new__(alias.AliasBankScreen)
        bank.tables, bank.odds, bank.grid = tables, odds, 128
        bank.floor_policy = 'positive-fallback'
        actual = np.minimum(subject.marker_regional(tables[0], odds[0]),
                            subject.marker_regional(tables[1], odds[1]))
        np.testing.assert_array_equal(actual, bank.regional(4096))

    def test_literals(self):
        values, odds = subject.literal_marker((-.347, 1.243))
        self.assertTrue(all(F(v) > 0 for v in values))
        self.assertEqual(odds, tuple(log2(float(F(v))) for v in values))
        for bad in [(float('nan'), 0), (33., 0), (float('inf'), 0)]:
            with self.assertRaises(ValueError):
                subject.literal_marker(bad)

    def test_complete_transfer_min_and_monotonicity(self):
        transfer = subject.previous.previous.counts_module.near.provider.profiles.transfer
        rng = np.random.default_rng(50)
        regions = [rng.normal(size=(3, 2, 2)), rng.normal(size=(3, 2, 2))]
        prefix = [.1, -.3]
        budget = SimpleNamespace(check=lambda: None)
        counts = SimpleNamespace(log_counts=lambda a, b, **kw: np.zeros(np.broadcast_shapes(a.shape, b.shape)))
        result = subject.envelope(regions, prefix, counts, budget, repetitions=2)
        terms = []
        for a in range(3):
            for b in range(3):
                witnesses = []
                for table, p in zip(regions, prefix):
                    product = transfer.power(transfer.multiply(table[a], table[b]), 2)
                    witnesses.append(p+np.logaddexp2(product[0, 0], product[0, 1]))
                terms.append(min(witnesses))
        self.assertAlmostEqual(result['score'], float(np.logaddexp2.reduce(terms)), places=11)
        self.assertEqual(result['profiles_evaluated'], 9)
        improved = [np.minimum(regions[0], regions[0]-.1), regions[1]]
        self.assertLessEqual(subject.envelope(improved, prefix, counts, budget, repetitions=2)['score'], result['score'])


if __name__ == '__main__':
    unittest.main()
