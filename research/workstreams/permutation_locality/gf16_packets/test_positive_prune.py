import unittest
from fractions import Fraction as Q
from positive_prune import repair,prune,pool,limiting_shells,union_gradient_logs
from math import exp
from shared_mixture import verify


class PositivePruneTests(unittest.TestCase):
    def test_exact_repair_and_componentwise_decrease(self):
        caps=[0,2,3,2,1]
        old=[(Q(64),Q(1,2)),(Q(32),Q(3,4))]
        for factors in ([0,0],[Q(1,100),Q(1,100)],[1,0],[1,1]):
            new,delta=repair(caps,old,factors)
            verify(caps,new)
            self.assertTrue(0<=delta<=1)
            by_p={p:c for c,p in new}
            self.assertTrue(all(by_p.get(p,0)<=c for c,p in old))
        with self.assertRaises(ValueError):repair(caps,old,[2,0])

    def test_lp_is_only_a_proposal(self):
        caps=[0,2,3,2,1]
        old=[(Q(64),Q(1,2)),(Q(32),Q(3,4))]
        new,info=prune(caps,old)
        verify(caps,new)
        self.assertEqual(info['status'],0)
        self.assertLess(sum(c for c,p in new),sum(c for c,p in old))

    def test_empty_required_measure(self):
        self.assertEqual(prune([0,0,0],[(Q(8),Q(1,2))])[0],[])

    def test_limiting_shell_diagnostic_is_sorted(self):
        rows=limiting_shells([0,1,2],[(Q(8),Q(1,2))])
        self.assertEqual(rows,[dict(support=2,ratio=1.),dict(support=1,ratio=4.)])
        self.assertEqual(limiting_shells([0,1,2],[(Q(8),Q(1,2))],1),rows[:1])
        with self.assertRaises(ValueError):limiting_shells([0,1],[(Q(2),Q(1,2))],0)

    def test_union_objective_matches_exact_gradient(self):
        mixture=[(Q(5),Q(1,3)),(Q(7),Q(3,4))];n=4;t=Q(3,16)
        values=union_gradient_logs(mixture,n,t)
        for (c,p),value in zip(mixture,values):
            expected=c*((1-p+p*t)**n+sum(d*(1-(p+q-p*q)+(p+q-p*q)*t)**n for d,q in mixture))
            self.assertAlmostEqual(exp(value),float(expected),places=12)
        with self.assertRaises(ValueError):union_gradient_logs(mixture,n,0)
        with self.assertRaises(ValueError):union_gradient_logs(mixture,0,t)

    def test_union_objective_retains_exact_shell_guarantee(self):
        caps=[0,2,3,2,1];old=[(Q(64),Q(1,2)),(Q(32),Q(3,4))]
        new,info=prune(caps,old,union_tilt=Q(3,16))
        verify(caps,new)
        self.assertEqual(info['union_cost_tilt'],'3/16')
        self.assertTrue(all(c<=dict((p,c) for c,p in old)[p] for c,p in new))

    def test_pool_preserves_each_majorant_and_merges_duplicates(self):
        caps=[0,1,1]
        result=pool(caps,[[(Q(4),Q(1,2)),(Q(4),Q(1,2))],[(Q(8),Q(3,4))]])
        self.assertEqual(result,[(Q(8),Q(1,2)),(Q(8),Q(3,4))])
        verify(caps,prune(caps,result)[0])


if __name__=='__main__':unittest.main()
