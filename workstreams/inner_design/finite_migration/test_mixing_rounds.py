"""Exact mixer-law and represented-measure checks, plus a screen cross-check."""
from collections import Counter
from fractions import Fraction as F
import math
import unittest

from flint import ctx
import mixing_rounds as study
import replay_mixing_rounds as replay


class MixingTests(unittest.TestCase):
    def test_native_recurrence_matches_scalar_and_exact_toy(self):
        from flint import arb
        ctx.prec = 128
        zero = (arb(1),arb(0),arb(0),arb(1)/2)
        one = (arb(0),arb(1)/3,arb(1)/5,arb(1)/7)
        native = replay.moments(zero,one,2,7)
        scalar = study.model.independent.single.moments(zero,one,2,7)
        for left,right in zip(native,scalar):
            self.assertTrue(left.overlaps(right))
        self.assertTrue(native[0].contains(arb(1)))

    def test_round_count_and_exact_transvection_product(self):
        with self.assertRaises(ValueError):
            study.epsilon('0')
        m = 7
        pairs = [(u,v) for u in range(1,m+1) for v in range(m+1) if (u&v).bit_count()%2 == 0]
        p = [[F() for _ in range(m)] for _ in range(m)]
        for q in range(1,m+1):
            for u,v in pairs:
                dest = q ^ (u if (q&v).bit_count()%2 else 0)
                p[q-1][dest-1] += F(1,len(pairs))
        current = [[F(i==j) for j in range(m)] for i in range(m)]
        for rounds in range(1,5):
            current = [[sum((current[i][k]*p[k][j] for k in range(m)),F())
                        for j in range(m)] for i in range(m)]
            eps = study.epsilon(str(rounds))
            for i in range(m):
                for j in range(m):
                    self.assertEqual(current[i][j],eps*(i==j)+(1-eps)/m)

    def test_all_transfer_rows_dominate_exact_laws(self):
        a,b = [1,2,3,4,5],[3,5,6,7,1]
        t,m = len(a),7
        image = lambda q: sum(((q&col).bit_count()&1)<<i for i,col in enumerate(a))
        weights = {q:image(q).bit_count() for q in range(1,m+1)}
        spectrum = Counter(weights.values())
        levels = sorted(spectrum)
        iw = [weights[q] for q in b]
        cw = [(image(q)^(1<<i)).bit_count() for i,q in enumerate(b)]
        laws = [(0,[F(q==0) for q in range(m+1)])]
        for q in range(1,m+1):
            atom = [F(y==q) for y in range(m+1)]
            laws.extend([(1,atom),(len(levels)+2+levels.index(weights[q]),atom)])
        for i,v in enumerate(levels):
            laws.append((i+2,[F(int(q!=0 and weights.get(q)==v),spectrum[v]) for q in range(m+1)]))
        for rounds in ('1','2','3','8','refresh'):
            eps = study.epsilon(rounds)
            for z in (F(1,2),F(9,10),F(1)):
                zero,one,n = study.transfers(spectrum,iw,cw,z,eps)
                for j,transfer in enumerate((zero,one)):
                    for index,law in laws:
                        actual = [F()]*(m+1)
                        inputs = [(0,0)] if j==0 else [(1<<i,b[i]) for i in range(t)]
                        for q,pq in enumerate(law):
                            for x,c in inputs:
                                mass = pq*z**((x^image(q)).bit_count())/len(inputs)
                                if q==0:
                                    actual[c] += mass
                                else:
                                    actual[q^c] += eps*mass
                                    for u in range(1,m+1):
                                        actual[u^c] += (1-eps)*mass/m
                        bound = transfer[index*n:(index+1)*n]
                        self.assertLessEqual(actual[0],bound[0])
                        residual = {q:max(F(),actual[q]-bound[2+levels.index(weights[q])]/spectrum[weights[q]])
                                    for q in range(1,m+1)}
                        uncovered = sum((max(F(),sum(residual[q] for q in residual if weights[q]==v)
                                             -bound[len(levels)+2+i]) for i,v in enumerate(levels)),F())
                        self.assertLessEqual(uncovered,bound[1],(rounds,z,j,index))

    def test_binary64_recurrence_matches_outward_evaluation(self):
        ctx.prec = 128
        engine = study.Engine(16,'2')
        tilt = F(-5)
        screen = engine.screen_coefficients(float(tilt))
        exact = study.audit.coefficients(engine,tilt)
        for w in (38,40,42,44,48,64):
            self.assertAlmostEqual(screen[w],math.log(float(exact[w])),places=9)


if __name__ == '__main__':
    unittest.main()
