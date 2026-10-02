import unittest
from itertools import accumulate, product
from fractions import Fraction as Q
import pairwise_macwilliams as mw
from joint_support import span


class PairwiseMacWilliamsTests(unittest.TestCase):
    def test_exact_small_code_rows_and_verified_duals(self):
        for basis,n in (([3,5],4),([3,5,9],4),([15,51,85],7),([1,2,4],3)):
            words = span(basis)
            dual = [x for x in range(1<<n) if all((x&r).bit_count()%2==0 for r in basis)]
            enumerators=[]
            for code in (words,dual):
                counts=[0]*(n+1)
                for x,y in product(code,repeat=2):
                    counts[(x|y).bit_count()]+=1
                enumerators.append(counts)
            # Loose caps: the equations, rather than exact shell priors,
            # must retain every true enumerator in the feasible region.
            cdfs=[[1]+[sum(row)]*n for row in enumerators]
            caps=[[1]+[sum(row)]*n for row in enumerators]
            system=mw.Constraints(cdfs,caps,len(basis))
            system.check_values(enumerators)
            exact_cdf_system=mw.Constraints([list(accumulate(row)) for row in enumerators],caps,len(basis))
            for u in range(1,n+1):
                self.assertGreaterEqual(system.transform_cap(u),enumerators[0][u])
                self.assertGreaterEqual(exact_cdf_system.transform_cap(u),enumerators[0][u])
                result=system.solve(u)
                self.assertGreaterEqual(result['cap'],enumerators[0][u])
                self.assertLessEqual(result['cap'],result['prior'])
                if 'dual' in result:
                    self.assertEqual(system.verify(u,result['dual'])[0],result['cap'])
                    fake=dict(result,cap=-10**100)  # Numerical caps are not proof inputs.
                    replay,checked=mw.replay_shell_witnesses(system,
                        dict(schema='pairwise-quaternary-shell-lp-1',rows=[fake]),system.caps[0])
                    self.assertEqual(checked,1)
                    self.assertEqual(replay[u],result['cap'])
            with self.assertRaises(ValueError):system.verify(1,[[0,'-1']])
            with self.assertRaises(ValueError):system.verify(1,[[0,'1'],[0,'1']])

    def test_forward_macwilliams_identity(self):
        for basis,n in (([3,5],4),([15,51,85],7)):
            code=span(basis)
            dual=[x for x in range(1<<n) if all((x&r).bit_count()%2==0 for r in basis)]
            laws=[]
            for words in (code,dual):
                laws.append([sum((x|y).bit_count()==u for x,y in product(words,repeat=2))
                             for u in range(n+1)])
            kw=[mw.pairwise_moments.krawtchouk4(n,u,n) for u in range(n+1)]
            for j in range(n+1):
                self.assertEqual(sum(c*kw[u][j] for u,c in enumerate(laws[0])),
                                 len(code)**2*laws[1][j])


if __name__ == '__main__':
    unittest.main()
