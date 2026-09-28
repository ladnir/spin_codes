"""Exact small-array checks for total input-weight band domination."""
from fractions import Fraction as Q
from itertools import product
from math import comb, log, prod
import unittest

from flint import ctx

from aggregate_weight import AggregateWeights


def rational(point):
    value = point.fmpq()
    return Q(int(value.numerator),int(value.denominator))


class AggregateWeightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 192

    def test_max_product_counts_individual_arrays(self):
        for length in range(1,5):
            caps = [Q(w+1,3) for w in range(length+1)]
            for rows in range(1,4):
                model = AggregateWeights(caps,0,length,rows,exclude_zero=False)
                expected = {}
                for weights in product(range(length+1),repeat=rows):
                    w = sum(weights)
                    value = prod(caps[a]/comb(length,a) for a in weights)
                    expected[w] = max(expected.get(w,Q(0)),value)
                self.assertEqual(model.exact_products(),expected)

    def test_every_toy_support_and_band(self):
        caps, length, rows = [Q(1),Q(2),Q(1,3),Q(4)],3,2
        model = AggregateWeights(caps,1,length,rows)
        for lower,upper in ((0,6),(2,2),(2,3),(4,6),(0,1)):
            for p,s in product((Q(1,5),Q(1,2),Q(4,5)),repeat=2):
                exact = model.exact_band(p,lower,upper)
                bound = model.majorant_exact(p,lower,upper,s)
                self.assertLessEqual(exact,bound)
                for words in product(range(1<<length),repeat=rows):
                    weights = tuple(word.bit_count() for word in words)
                    w = sum(weights)
                    original = (prod(caps[a]/comb(length,a) for a in weights)
                                if min(weights)>=1 and lower<=w<=upper else Q(0))
                    reference = p**w*(1-p)**(length*rows-w)
                    self.assertLessEqual(original,bound*reference)

    def test_outward_replay_and_float_proposals(self):
        model = AggregateWeights([1,0,Q(7,3),2,0,Q(11,5)],1,5,3)
        for lower,upper in ((0,15),(3,6),(7,12),(14,15)):
            for p,s in product((Q(1,1000),Q(1,2),Q(123,127)),repeat=2):
                exact = model.majorant_exact(p,lower,upper,s)
                outward = model.majorant_arb(p,lower,upper,s)
                self.assertGreaterEqual(rational(outward),exact)
                if exact:
                    self.assertAlmostEqual(model.majorant_log(p,lower,upper,s),
                                           log(exact.numerator)-log(exact.denominator),places=9)
                else:
                    self.assertEqual(model.majorant_log(p,lower,upper,s),float('-inf'))
                best,witness = model.best_arb(p,lower,upper,(s,p))
                self.assertIn(witness,(s,p))
                self.assertLessEqual(best,outward)
                self.assertGreaterEqual(rational(best),model.exact_band(p,lower,upper))

    def test_s_equals_p_recovers_original(self):
        model = AggregateWeights([1,0,6,0,1],0,4,3,exclude_zero=False)
        for p in (Q(1,5),Q(1,2),Q(4,5)):
            self.assertEqual(model.majorant_exact(p,0,12,p),model.gamma(p)**3)

    def test_band_can_strictly_improve_scalar(self):
        model = AggregateWeights([1,0,6,0,1],0,4,2,exclude_zero=False)
        p,s = Q(1,4),Q(1,2)
        self.assertLess(model.majorant_exact(p,0,2,s),model.gamma(p)**2)

    def test_disjoint_bands_with_different_reference_moments(self):
        caps, length, rows = [1,2,3,1],3,2
        model = AggregateWeights(caps,1,3,rows)
        bands = ((2,3,Q(1,3)),(4,4,Q(1,2)),(5,6,Q(3,4)))
        actual = Q(0)
        moments = [Q(0) for _ in bands]
        for words in product(range(1<<length),repeat=rows):
            weights = tuple(word.bit_count() for word in words)
            total = sum(weights)
            function = Q(1+(words[0] ^ (3*words[1]))**2,7)
            if min(weights)>=1:
                actual += function*prod(Q(caps[w],comb(length,w)) for w in weights)
            for i,(_,_,p) in enumerate(bands):
                moments[i] += function*p**total*(1-p)**(length*rows-total)
        bound = sum(model.majorant_exact(p,lo,hi,Q(1,2))*moment
                    for (lo,hi,p),moment in zip(bands,moments))
        self.assertLessEqual(actual,bound)

    def test_empty_support_gaps_and_operation_guard(self):
        empty = AggregateWeights([1,0,0],1,2,3)
        self.assertEqual(empty.exact_products(),{})
        self.assertEqual(empty.majorant_exact(Q(1,2),0,6,Q(1,2)),0)
        gaps = AggregateWeights([1,0,1],0,2,2,exclude_zero=False)
        self.assertEqual(gaps.exact_band(Q(1,2),1,1),0)
        self.assertGreater(gaps.majorant_exact(Q(1,2),1,1,Q(1,2)),0)
        with self.assertRaises(ValueError):
            gaps.exact_products(max_operations=1)
        for p in (0,1,-1):
            with self.assertRaises(ValueError):
                gaps.majorant_arb(p,0,4,Q(1,2))
        with self.assertRaises(ValueError):
            gaps.majorant_arb(Q(1,2),3,2,Q(1,2))


if __name__ == '__main__':
    unittest.main()
