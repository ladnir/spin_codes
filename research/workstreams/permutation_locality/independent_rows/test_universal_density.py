"""Complete maxima, coefficient weights, and envelope-scope tests."""
from fractions import Fraction as Q
from itertools import product
import unittest

from flint import arb, arb_mat, ctx

from column_density import bound
from universal_density import shapes, coefficients, refine
from occupancy_memory import M, C, U
from mature_tail import L48, L56
from group_rank_one_verify import up


SPECTRUM = {48:5, 56:7, 64:11, 72:13, 80:17}
WINDOWS = {b:(1000, b*5, b*7, 1000) for b in range(1, 5)}


def feedback(maximum=3, rounds=2):
    records = {}
    for j in range(1, maximum+1):
        for weights in shapes(j):
            value = arb(sum(weights)+weights.count(4)+1)/1024
            records[weights] = dict(density=up(value), uniform={v:up(value/h) for v,h in SPECTRUM.items()})
    return dict(bounds={'.052':records}, provenance=dict(maximum=maximum, rounds=rounds,
                spectrum=SPECTRUM.copy(), checked_shapes=len(records)))


def matrices(maximum=4):
    result = [arb_mat([[1]*11 for _ in range(11)]) for _ in range(maximum+1)]
    for matrix in result:
        for source in (M, L48, L56):
            matrix[source, C] = 0
    return result


class UniversalDensity(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.precision

    def test_independent_per_shape_oracle_and_weights(self):
        data = feedback()
        for rho, a in ((Q(1),Q(1)), (Q(3,4),Q(9,8)), (Q(9,10),Q(7,8))):
            for windows, exact in ((WINDOWS,None), (None,data), (WINDOWS,data)):
                actual = coefficients(4,SPECTRUM,'.052',rho,a,window_averages=windows,feedback=exact)
                for j, row in actual.items():
                    oracle = [[] for _ in range(6)]
                    # High precision and a separate raw bound implementation
                    # avoid testing the adapter against its own power cache.
                    ctx.prec = 512
                    for weights in shapes(j):
                        candidate = None
                        if windows is not None and j >= 2:
                            col = bound(weights,windows,'.052')/4
                            candidate = [col, *(col/h for h in SPECTRUM.values())]
                        if exact is not None and j <= 3:
                            record = exact['bounds']['.052'][weights]
                            values = [record['density'], *record['uniform'].values()]
                            candidate = values if candidate is None else [min(x,y) for x,y in zip(candidate,values)]
                        scale = rho**weights.count(4)*a**sum(weights)
                        for i,value in enumerate(candidate):
                            oracle[i].append(value*scale.numerator/scale.denominator)
                    for value, candidates in zip(row,oracle):
                        expected = max(candidates)
                        self.assertTrue(value >= expected.upper(),(j,rho,a,windows is not None,exact is not None,value,expected))
                        self.assertAlmostEqual(float(value),float(expected),places=14)
                        self.assertTrue(value.is_exact())
                    ctx.prec = 192

    def test_penalty_is_applied_before_shape_maximum(self):
        data = feedback(1)
        for weights,record in data['bounds']['.052'].items():
            record['density'] = arb(8 if weights == (4,) else 4)/16
        row = coefficients(1,SPECTRUM,'.052','1/4',feedback=data)[1]
        self.assertEqual(row[0],arb(1)/4)
        # Feedback records include alpha: do not apply another factor 1/4.
        self.assertEqual(coefficients(1,SPECTRUM,'.052',feedback=data)[1][0],arb(1)/2)

    def test_intersection_before_maximum(self):
        data = feedback(2)
        # Make the two competing bounds favor different shapes.
        for weights,record in data['bounds']['.052'].items():
            record['density'] = arb(1)/1024 if weights == (4,4) else arb(1)
        combined = coefficients(2,SPECTRUM,'.052',window_averages=WINDOWS,feedback=data)[2][0]
        column = coefficients(2,SPECTRUM,'.052',window_averages=WINDOWS)[2][0]
        census = coefficients(2,SPECTRUM,'.052',feedback=data)[2][0]
        self.assertLess(combined,min(column,census))

    def test_only_requested_entries_and_occupancies_change(self):
        before = matrices()
        preserved = [t*1 for t in before]
        after = refine(before,SPECTRUM,'.052',feedback=feedback(2))
        allowed = {(C,C)} | {(U+i,C) for i in range(5)}
        for j in range(5):
            for source,target in product(range(11),repeat=2):
                if 1 <= j <= 2 and (source,target) in allowed:
                    self.assertLess(after[j][source,target],before[j][source,target])
                else:
                    self.assertEqual(after[j][source,target],before[j][source,target])
            self.assertEqual(before[j],preserved[j])
        column_only = refine(before,SPECTRUM,'.052',window_averages=WINDOWS)
        self.assertEqual(column_only[0],before[0])
        self.assertEqual(column_only[1],before[1])
        self.assertEqual(refine(before,SPECTRUM,'.052'),before)

    def test_reject_incomplete_or_mismatched_census(self):
        for mutation in ('missing','extra','count','rounds','spectrum','uniform'):
            data = feedback()
            records = data['bounds']['.052']
            if mutation == 'missing':
                del records[(4,4,4)]
            elif mutation == 'extra':
                records[(1,1,1,1)] = records[(1,)]
            elif mutation == 'uniform':
                del records[(1,)]['uniform'][48]
            elif mutation == 'count':
                data['provenance']['checked_shapes'] -= 1
            elif mutation == 'rounds':
                data['provenance']['rounds'] = 1
            else:
                data['provenance']['spectrum'][48] += 1
            with self.assertRaises(ValueError):
                coefficients(4,SPECTRUM,'.052',feedback=data)
        with self.assertRaises(ValueError):
            coefficients(4,SPECTRUM,'.052',window_averages={1:WINDOWS[1]})

    def test_reject_coupled_density_and_invalid_parameters(self):
        for source in (M,L48,L56):
            before = matrices()
            before[2][source,C] = 1
            with self.assertRaises(ValueError):
                refine(before,SPECTRUM,'.052')
        for maximum,tilt,rho,a in ((0,0,1,1),(33,0,1,1),(4,-1,1,1),(4,0,0,1),(4,0,2,1),(4,0,1,0)):
            with self.assertRaises(ValueError):
                coefficients(maximum,SPECTRUM,tilt,rho,a)
        with self.assertRaises(ValueError):
            refine([arb_mat(9,9)]*3,SPECTRUM,'.052')


if __name__ == '__main__':
    unittest.main()
