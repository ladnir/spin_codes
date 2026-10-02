import unittest
from fractions import Fraction as Q
from itertools import combinations
from math import comb

import pair_mixer as mixer


class PairMixerTests(unittest.TestCase):
    def test_pair_support_and_iid_conditional_labels(self):
        expected = {
            (0,0): {(False,False): Q(1)},
            (1,0): {(True,True): Q(1)},
            (0,7): {(True,True): Q(1)},
            (3,9): {(False,True): Q(1,15), (True,False): Q(1,15), (True,True): Q(13,15)},
        }
        for alpha in (2, 7, 15):
            for pair, supports in expected.items():
                law = mixer.pair_law(*pair, alpha=alpha)
                self.assertEqual(sum(law.values()), 1)
                actual = {s: sum(p for (x,y),p in law.items() if (bool(x),bool(y)) == s) for s in supports}
                self.assertEqual(actual, supports)
                for (x,y), probability in law.items():
                    support = bool(x), bool(y)
                    self.assertEqual(probability, supports[support]/15**sum(support))

    def test_without_postmultiplier_labels_are_correlated(self):
        law = mixer.pair_law(1, 0, postmultiply=False)
        self.assertEqual(len(law), 15)
        self.assertTrue(all(x == y for x,y in law))

    def test_small_kernel_exhaustive_support_placements(self):
        for n in (2, 4, 6, 8):
            kernel = mixer.support_kernel(n)
            for u in range(n+1):
                expected = [Q(0)]*(n+1)
                for positions in combinations(range(n),u):
                    occupied = set(positions)
                    j = sum(2*i in occupied and 2*i+1 in occupied for i in range(n//2))
                    for loss in range(j+1):
                        v = 2*u-2*j-loss
                        expected[v] += Q(comb(j,loss)*2**loss*13**(j-loss), comb(n,u)*15**j)
                self.assertEqual(kernel[u],expected)
                self.assertEqual(sum(v*p for v,p in enumerate(kernel[u])), mixer.mean_support(u,n))

    def test_kernel_is_not_stochastically_monotone(self):
        kernel = mixer.support_kernel(16)
        self.assertGreater(sum(v*p for v,p in enumerate(kernel[15])), sum(v*p for v,p in enumerate(kernel[16])))

    def test_gf256_pair_law_and_monotonicity(self):
        # Uniform nonzero eight-bit values give independent nonzero nibble
        # labels conditional on which nibbles are active.
        supports = {}
        for value in range(1,256):
            support = bool(value & 15), bool(value >> 4)
            supports[support] = supports.get(support,0)+1
        self.assertEqual(supports, {(True,False):15, (False,True):15, (True,True):225})
        for n in (2,4,8,16):
            kernel = mixer.gf256_support_kernel(n)
            for u in range(n+1):
                self.assertEqual(sum(kernel[u]),1)
                self.assertEqual(sum(v*p for v,p in enumerate(kernel[u])),
                    Q(32,17)*u-Q(16*u*(u-1),17*(n-1)))
                if u:
                    self.assertTrue(all(sum(kernel[u][:v+1]) <= sum(kernel[u-1][:v+1]) for v in range(n+1)))

    def test_cdf_transport_dominates_every_single_atom_placement(self):
        n = 8
        kernel = mixer.support_kernel(n)
        # An upper CDF with a unit jump at u allows the real atom anywhere
        # in u..n. The suffix envelope must dominate all such placements.
        for u in range(1,n+1):
            cdf = [int(v >= u) for v in range(n+1)]
            bound = mixer.transport_cdf(cdf,kernel)
            for location in range(u,n+1):
                for v in range(n+1):
                    self.assertGreaterEqual(bound[v], sum(kernel[location][:v+1]))

    def test_transport_is_safe_but_naive_cdf_difference_can_fail(self):
        n = 16
        kernel = mixer.support_kernel(n)
        cdf = [int(v >= n-1) for v in range(n+1)]
        actual = [sum(kernel[n][:v+1]) for v in range(n+1)]
        naive = [sum(kernel[n-1][:v+1]) for v in range(n+1)]
        self.assertTrue(any(a>b for a,b in zip(actual,naive)))
        self.assertTrue(all(a<=b for a,b in zip(actual,mixer.transport_cdf(cdf,kernel))))

    def test_fixed_witness_cutoff_retarget(self):
        from flint import arb, ctx
        from pair_mixer_retarget import retarget
        ctx.prec = 128
        point = dict(upper=[3,-10], witness=dict(parameters=['1/16']))
        same = retarget(point,100,100)
        self.assertEqual(Q(same['upper'][0])*Q(2)**same['upper'][1],Q(3,1024))
        lower = retarget(point,100,99)
        self.assertLess(Q(lower['upper'][0])*Q(2)**lower['upper'][1],Q(3,1024))
        with self.assertRaises(ValueError):
            retarget(point,100,-1)


if __name__ == '__main__':
    unittest.main()
