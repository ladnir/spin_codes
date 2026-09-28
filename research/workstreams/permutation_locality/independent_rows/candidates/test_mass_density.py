"""Exact local-model checks for the unintegrated mass-density candidate."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
from pathlib import Path
import sys
import unittest

from flint import arb, arb_mat, ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from test_column_density import linear, packets
from mass_density import coefficients, nested_coefficients, blend, conditioned_peaks, conditioned_coefficients
from occupancy_memory import Z, M, C
from mature_tail import L48, L56


def census(maximum=3):
    return {shape:(1,sum(shape)+shape.count(4),128,{})
            for j in range(1,maximum+1)
            for shape in combinations_with_replacement(range(1,5),j)}


def operators(maximum=4):
    matrices = [arb_mat([[arb(1)/8]*11 for _ in range(11)]) for _ in range(maximum+1)]
    for matrix in matrices:
        for source in (M,L56,L48):
            matrix[source,C] = 0
    return matrices


class MassDensity(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.precision

    def test_nested_domination_including_nonmonotone_moments(self):
        for moments in product((Q(0),Q(1,8),Q(1,2),Q(1)),repeat=5):
            levels = (48,56,64,72,80)
            a,b,c = nested_coefficients(dict(zip(levels,moments)))
            self.assertTrue(min(a,b,c) >= 0)
            for v,f in zip(levels,moments):
                self.assertLessEqual(f,a+b*(v<=56)+c*(v<=48))

    def test_exact_toy_density_for_all_targets_and_shapes(self):
        cases = [(3,(0b010101,0b001111),(1,2,3,1,3,2)),
                 (4,(0b01010101,0b00110011),(1,2,3,1,2,3,3,1)),
                 (3,(0b001,0b011),(0,)*6)]
        for windows,expansion,feedback in cases:
            weights = [linear(expansion,s).bit_count() for s in range(4)]
            for shape in [h for j in range(1,min(3,windows)+1)
                          for h in combinations_with_replacement((1,2),j)]:
                inputs = list(packets(shape,windows,feedback))
                counts = Counter(syndrome for _,syndrome in inputs)
                beta = Q(max(counts.values()),len(inputs))
                for z in (Q(1,2),Q(4,5),Q(1)):
                    a,b,c = nested_coefficients({v:z**abs(v-sum(shape)) for v in set(weights[1:])},1,2)
                    measures = [(Q(0),Q(1,7),Q(2,7),Q(4,7))]
                    measures += [tuple(Q(int(s==k)) for s in range(4)) for k in range(1,4)]
                    for mu in measures:
                        mass = sum(mu)
                        low = sum(mu[s] for s in range(1,4) if weights[s]<=1)
                        mid = sum(mu[s] for s in range(1,4) if weights[s]<=2)
                        upper = beta*(a*mass+b*mid+c*low)/4
                        for target in range(4):
                            actual = sum((mu[target^syndrome]
                                          *z**(linear(expansion,target^syndrome)^word).bit_count()
                                          for word,syndrome in inputs),Q(0))/(4*len(inputs))
                            self.assertLessEqual(actual,upper)

    def test_zero_feedback_atom_cannot_be_omitted(self):
        data = census(1)
        data[(4,)] = (128,0,128,{})
        result = coefficients(data,1,'0')[1]
        self.assertEqual(result,(arb(1)/4,arb(0),arb(0)))

    def test_exact_weighted_maxima_and_single_lazy_factor(self):
        data = census()
        for rho,a in ((Q(1),Q(1)),(Q(3,4),Q(9,8)),(Q(9,10),Q(7,8))):
            result = coefficients(data,3,0,rho,a)
            for j in range(1,4):
                expected = max(Q(max(z,p),den)*rho**shape.count(4)*a**sum(shape)/4
                               for shape,(z,p,den,_) in data.items() if len(shape)==j)
                ctx.prec=512
                self.assertGreaterEqual(result[j][0],arb(expected.numerator)/expected.denominator)
                self.assertAlmostEqual(float(result[j][0]),float(expected),places=14)
                self.assertEqual(result[j][1:],(arb(0),arb(0)))
                ctx.prec=192
        self.assertTrue(all(x.is_exact() and x>=0
                            for row in coefficients(data,3,'.056').values() for x in row))

    def test_convex_column_endpoints_scope_and_no_mutation(self):
        before = operators()
        saved = [x*1 for x in before]
        candidates = coefficients(census(),3,'.056')
        self.assertEqual(blend(before,candidates,{2:Q(0)}),before)
        for fraction in (Q(1,4),Q(1,2),Q(1)):
            after = blend(before,candidates,{2:fraction})
            for j in range(5):
                for s,t in product(range(11),repeat=2):
                    if j==2 and t==C and s in (M,L56,L48,C):
                        if s==C:
                            expected = before[j][C,C]*(fraction.denominator-fraction.numerator)/fraction.denominator
                        else:
                            expected = candidates[2][(M,L56,L48).index(s)]*fraction.numerator/fraction.denominator
                        self.assertTrue(after[j][s,t] >= expected)
                        self.assertAlmostEqual(float(after[j][s,t]),float(expected),places=14)
                    else:
                        self.assertEqual(after[j][s,t],before[j][s,t])
            self.assertEqual(before,saved)

    def test_incomplete_census_coupled_column_and_bad_parameters(self):
        data = census()
        del data[(4,4,4)]
        with self.assertRaises(ValueError): coefficients(data,3,'.056')
        data = census()
        data[(1,)] = (0,129,128,{})
        with self.assertRaises(ValueError): coefficients(data,3,'.056')
        for value in (-1,2,True,.5):
            with self.assertRaises(ValueError): blend(operators(),coefficients(census(),3,'.056'),{2:value})
        before=operators(); before[2][M,C]=1
        with self.assertRaises(ValueError): blend(before,coefficients(census(),3,'.056'),{2:Q(1,2)})
        for j in (0,5):
            with self.assertRaises(ValueError): blend(operators(),coefficients(census(),3,'.056'),{j:Q(1,2)})

    def test_conditioned_peaks_dominate_exact_small_models(self):
        for windows,feedback in ((3,(1,2,3,1,3,2)),(4,(1,2,3,1,2,3,3,1)),(4,(0,)*8)):
            exact = {}
            for j in range(1,windows+1):
                for shape in combinations_with_replacement((1,2),j):
                    counts = Counter(syndrome for _,syndrome in packets(shape,windows,feedback))
                    exact[shape] = (counts[0],max((c for s,c in counts.items() if s),default=0),sum(counts.values()),{})
            for through in range(1,windows+1):
                partial = {shape:record for shape,record in exact.items() if len(shape)<=through}
                bounds = conditioned_peaks(partial,through,windows,windows=windows,packet_bits=2)
                self.assertEqual(set(bounds),set(exact))
                for shape,(zero,peak,den,_) in exact.items():
                    self.assertGreaterEqual(bounds[shape],Q(max(zero,peak),den))
                    self.assertLessEqual(bounds[shape],1)
                    if len(shape)<=through:
                        self.assertEqual(bounds[shape],Q(max(zero,peak),den))

    def test_zero_target_preserves_refresh_and_other_components(self):
        before=operators()
        candidates=coefficients(census(),3,'.056')
        for fraction in (Q(0),Q(1,4),Q(1)):
            after=blend(before,candidates,{2:fraction},target=Z)
            for j in range(5):
                for s,t in product(range(11),repeat=2):
                    if j==2 and t==Z and s in (M,L56,L48,C):
                        if s==C:
                            expected=before[j][C,Z]*(fraction.denominator-fraction.numerator)/fraction.denominator
                        else:
                            expected=before[j][s,Z]+candidates[j][(M,L56,L48).index(s)]*fraction.numerator/fraction.denominator
                        self.assertGreaterEqual(after[j][s,t],expected)
                        self.assertAlmostEqual(float(after[j][s,t]),float(expected),places=14)
                    else:
                        self.assertEqual(after[j][s,t],before[j][s,t])

    def test_conditioned_peaks_all_subshapes_and_zero_atom(self):
        data=census()
        data[(4,4)] = (128,0,128,{})
        actual=conditioned_peaks(data,3,7,windows=8)
        for shape,value in actual.items():
            if len(shape)<=3:
                z,p,d,_=data[shape]
                self.assertEqual(value,Q(max(z,p),d))
                continue
            candidates=[Q(1)]
            for sub,(z,p,d,_) in data.items():
                if all(sub.count(b)<=shape.count(b) for b in range(1,5)):
                    from math import comb
                    candidates.append(Q(max(z,p)*comb(8,len(sub)),d*comb(8-len(shape)+len(sub),len(sub))))
            self.assertEqual(value,min(candidates))
        low=coefficients(data,3,'.056')
        extended=conditioned_coefficients(data,3,7,'.056')
        self.assertEqual({j:extended[j] for j in low},low)
        for through,maximum in ((0,2),(4,3),(3,33)):
            with self.assertRaises(ValueError): conditioned_peaks(data,through,maximum)

    def test_combined_zero_and_density_retain_their_separate_columns(self):
        before=operators()
        candidate=coefficients(census(),3,'.056')
        fractions={2:Q(1),3:Q(1,2)}
        zero=blend(before,candidate,fractions,target=Z)
        density=blend(before,candidate,fractions,target=C)
        both=blend(zero,candidate,fractions,target=C)
        for j in range(len(before)):
            for s,t in product(range(11),repeat=2):
                expected=zero[j][s,t] if t==Z else density[j][s,t] if t==C else before[j][s,t]
                self.assertEqual(both[j][s,t],expected)
        self.assertEqual(before,operators())


if __name__ == '__main__':
    unittest.main()
