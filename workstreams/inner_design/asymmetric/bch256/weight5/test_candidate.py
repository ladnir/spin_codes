"""Exact candidate binding, cancellation, and transfer regression checks."""
from collections import Counter
from fractions import Fraction as F
import itertools
import unittest
from flint import arb, ctx
import candidate
model = candidate.model


class CandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 256
        cls.engine = candidate.Engine(20)

    def test_matches_retained_candidate(self):
        saved = model.base.read(model.HERE/'overlap/KERNEL_DIAGNOSIS.json')
        selected = next(r for r in saved['results'] if r['name'] == 'weight5_seed0')
        self.assertEqual(self.engine.columns, selected['columns'])
        self.assertEqual(self.engine.kernel[4], selected['zero_kernel_weight4'])
        self.assertEqual(min(self.engine.b_spectrum), selected['dual_minimum_weight'])
        self.assertEqual(set(c.bit_count() for c in self.engine.columns), {5})

    def test_exact_small_weight_fibers(self):
        engine = self.engine
        for j in (1,2,3):
            fibers = Counter()
            for support in itertools.combinations(range(128),j):
                value = 0
                for i in support:
                    value ^= engine.columns[i]
                fibers[value] += 1
            self.assertEqual(fibers.pop(0,0), engine.kernel[j])
            self.assertLessEqual(max(fibers.values()), int(engine.caps[j]['cap']))

    def test_independent_cancellation_histograms(self):
        engine = self.engine
        search = model.independent.search
        for j in (1,2):
            by_weight = {}
            for support in itertools.combinations(range(128),j):
                x = sum(1 << i for i in support)
                q = search.inject(engine.columns,x)
                aq = search.image(engine.a_columns,q)
                by_weight.setdefault(aq.bit_count(),Counter())[(x^aq).bit_count()] += 1
            self.assertEqual(by_weight, engine.low[j]['by_weight'])

    def test_zero_row_and_linear_q1(self):
        engine = candidate.Engine(16)
        rows = engine.epoch(F(-3))
        z = (-model.number(F(-3)).exp()).exp()
        import math
        for j in (0,1,2,4,64,128):
            exact = arb(engine.kernel[j])/math.comb(128,j)*z**j
            self.assertLessEqual(exact, rows[j][0])
        first = engine.q1(F(-5))
        second = engine.q1(F(-5),linear=True)
        for a,b in zip(first,second):
            self.assertLess(abs(float(a/b)-1),1e-12)

    def test_checker_uses_candidate_data(self):
        checker = candidate.Checker(20)
        self.assertEqual(checker.engine.columns, self.engine.columns)
        self.assertEqual(checker.engine.caps, self.engine.caps)
        self.assertEqual(checker.engine.low, self.engine.low)
        gap = model.base.read(model.HERE/'DENSE_GAP_POINT.json')
        q,v = gap['occupation'],F(gap['coordinate'])
        self.assertLess(checker.bound(q,q,v,v,gap['witness']), -700*arb(2).log())


if __name__ == '__main__':
    unittest.main()
