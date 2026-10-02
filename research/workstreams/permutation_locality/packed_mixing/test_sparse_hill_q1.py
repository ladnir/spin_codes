from fractions import Fraction as Q
from itertools import product
import unittest

from flint import arb, ctx
from canonical_counts import transport_cdf
import sparse_hill_q1 as hill


class SparseHillQ1Tests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        ctx.prec = 256

    def test_monotone_weights_commute_with_transport(self):
        local = (Q(0), Q(1, 4), Q(3, 4))
        caps = [0, 2, 5, 9]
        laws = hill.support_laws(local, 3)
        weights = [hill.aq(Q(1, 3)**u) for u in range(7)]
        value, terms = hill.canonical_terms(caps, laws, weights, 1)
        transported = transport_cdf(caps, local)
        exact = sum((transported[u]-transported[u-1])*Q(1, 3)**u for u in range(1, 7))
        mantissa, exponent = hill.endpoint(value)
        upper = Q(mantissa)*Q(2)**exponent
        self.assertGreaterEqual(upper, exact)
        self.assertLess(upper-exact, Q(1, 1 << 240))
        self.assertEqual(len(terms), 3)

    def test_nonmonotone_weights_valid_for_every_compatible_count_profile(self):
        caps = [0, 1, 2, 3]
        laws = hill.support_laws((Q(0), Q(1, 2), Q(1, 2)), 3)
        exact_weights = [0, 1, 0, 2, 0, 1, 0]
        weights = [arb(v) for v in exact_weights]
        bound, _ = hill.canonical_terms(caps, laws, weights, 1)
        for counts in product(range(4), repeat=3):
            if any(sum(counts[:h]) > caps[h] for h in range(1, 4)):
                continue
            exact = sum(counts[h-1]*sum(p*exact_weights[u] for u, p in enumerate(laws[h])) for h in range(1, 4))
            self.assertTrue(hill.aq(exact) <= bound)

    def test_bad_caps_laws_and_weights_rejected(self):
        laws = hill.support_laws((Q(0), Q(1)), 3)
        for caps, weights in (([1, 2, 3, 4], [arb(1)]*4), ([0, 2, 1, 3], [arb(1)]*4),
                ([0, 1, 2], [arb(1)]*4), ([0, 1, 2, 3], [arb(1)]*3),
                ([0, 1, 2, 3], [arb(1), arb(-1), arb(1), arb(1)])):
            with self.assertRaises(ValueError):
                hill.canonical_terms(caps, laws, weights)
        with self.assertRaises(ValueError):
            hill.support_laws((Q(1, 2), Q(1, 2)), 3)


if __name__ == '__main__':
    unittest.main()
