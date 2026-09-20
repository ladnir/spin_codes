from fractions import Fraction as F
import itertools
import math
import unittest

from flint import arb, ctx

import ladder_dense_tilt as tilted


class TiltTests(unittest.TestCase):
    def test_change_of_measure_exact(self):
        # Includes zero and all-one rows, so endpoint handling is exercised.
        qs = [F(0),F(1,4),F(2,3),F(1)]
        for x in (F(1,5),F(1),F(7,2)):
            ps = [q/(x+(1-x)*q) for q in qs]
            normalizer = math.prod(1-p+p*x for p in ps)
            for bits in itertools.product((0,1),repeat=len(ps)):
                old = math.prod(p if bit else 1-p for p,bit in zip(ps,bits))
                new = math.prod(q if bit else 1-q for q,bit in zip(qs,bits))
                self.assertEqual(old,normalizer*x**(-sum(bits))*new)

    def test_tilted_row_cost_exact_identity(self):
        for q in (F(1,4),F(3,5)):
            for x in (F(1,3),F(2)):
                p = q/(x+(1-x)*q)
                z = 1-p+p*x
                for w in range(9):
                    self.assertEqual(z**8/(p**w*(1-p)**(8-w)),
                                     x**w/(q**w*(1-q)**(8-w)))

    def test_zero_input_tilt_matches_old_bound(self):
        ctx.prec = 256
        ps = tilted.labels.row_probabilities(.5)
        checker = tilted.Checker(tilted.labels.identity.instance(20),ps)
        old = tilted.labels.Checker(checker.spec,ps)
        witness = dict(tilt=-100,input_tilt=0,eta=tilted.base.encode(F(10)))
        a = checker.bound(512,512,F(1,3),F(1,3),witness)
        b = old.bound(512,512,F(1,3),F(1,3),witness)
        self.assertLess(abs(float(a-b)),1e-8)

    def test_all_one_normalizer_cancels(self):
        ctx.prec = 256
        for k in (-80,0,80):
            log_x = arb(k)/40
            costs = tilted.tilted_costs([],[],{},log_x)
            self.assertEqual(len(costs),1)
            self.assertTrue((costs[0]+256*((-log_x).exp()).log()).contains(0))

    def test_coupled_scalar_upper(self):
        ctx.prec = 256
        rows, lo, hi = 128,30,75
        nlo,nhi = F(1,4),F(7,8)
        for x in (F(1,7),F(1),F(5)):
            for eta in (F(-500),F(0),F(320)):
                log_s = arb(40)
                upper = tilted.scalar_bound(rows,lo,hi,nlo,nhi,log_s,eta,tilted.aa(x))
                for q in range(lo,hi+1):
                    for i in range(11):
                        nu = nlo+(nhi-nlo)*F(i,10)
                        theta = F(q,rows)*nu
                        actual = q*(log_s+tilted.aa(eta*nu))+256*rows*(1+tilted.aa((1/x-1)*theta)).log()
                        self.assertGreaterEqual(upper,actual.lower())


if __name__ == '__main__':
    unittest.main()
