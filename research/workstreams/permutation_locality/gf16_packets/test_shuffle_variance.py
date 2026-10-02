from fractions import Fraction as Q
from itertools import combinations_with_replacement,product
from math import comb
import unittest
from flint import ctx
import shuffle_variance as sv
from test_scalar_cover import endpoint


class ShuffleVarianceTests(unittest.TestCase):
    def test_ratio_mode_and_outward_bound_on_exact_grid(self):
        ctx.prec=192
        grid=[Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)]
        for n in range(1,7):
            for probabilities in combinations_with_replacement(grid,n):
                mean=sum(probabilities);variance=sum(p*(1-p) for p in probabilities)
                if mean in (0,n):
                    self.assertEqual(sv.density_upper(n,mean,0),1)
                    continue
                masses=[Q(1)]
                for p in probabilities:
                    updated=[Q(0)]*(len(masses)+1)
                    for k,m in enumerate(masses):updated[k]+=m*(1-p);updated[k+1]+=m*p
                    masses=updated
                p=mean/n
                ratios=[mass/(comb(n,k)*p**k*(1-p)**(n-k)) for k,mass in enumerate(masses)]
                lo=mean.numerator//mean.denominator;hi=-(-mean.numerator//mean.denominator)
                self.assertEqual(max(ratios),max(ratios[lo],ratios[hi]))
                bound=endpoint(sv.density_upper(n,mean,variance))
                self.assertGreaterEqual(bound,max(ratios))
                self.assertGreaterEqual(endpoint(sv.density_upper(n,mean,variance/2)),bound)

    def test_validation(self):
        for args in ((0,0,0),(4,-1,0),(4,5,0),(4,2,-1),(4,2,2)):
            with self.assertRaises(ValueError):sv.density_upper(*args)

    def test_interval_bound_at_means_and_integer_boundaries(self):
        ctx.prec=192
        for lo,hi in ((Q(1,5),Q(7,5)),(Q(2),Q(4)),(Q(9,4),Q(11,4))):
            bound=endpoint(sv.density_interval_upper(8,(lo,hi),Q(1,10)))
            for i in range(21):
                mean=lo+(hi-lo)*i/20
                self.assertGreaterEqual(bound,endpoint(sv.density_upper(8,mean,Q(1,10))))
        for interval in ((0,1),(1,8),(2,1)):
            with self.assertRaises(ValueError):sv.density_interval_upper(8,interval,0)

    def test_variance_dual_covers_every_small_composition(self):
        features=[Q(0),Q(1,4),Q(2,3),Q(1)];active=[0,1,1,1];n=5
        for cell in ((Q(1,10),Q(2,5)),(Q(1,3),Q(2,3)),(Q(1,2),Q(1,2))):
            lower,dual=sv.variance_dual(features,active,cell,Q(2,n))
            self.assertEqual((lower,dual),sv.variance_dual(features,active,cell,Q(2,n),dual))
            for counts in product(range(n+1),repeat=4):
                if sum(counts)!=n or sum(c*a for c,a in zip(counts,active))<2:continue
                mean=sum(c*f for c,f in zip(counts,features))/n
                if not cell[0]<=mean<=cell[1]:continue
                variance=sum(c*f*(1-f) for c,f in zip(counts,features))/n
                self.assertLessEqual(lower,variance)
            # Any rational slopes are safe after checking the intercept;
            # they need not be numerically optimal.
            self.assertEqual(sv.variance_dual(features,active,cell,Q(2,n),[0,0])[0],0)
        with self.assertRaises(ValueError):sv.variance_dual(features,active,(0,1),Q(2,n),[0,-1])


if __name__=='__main__':unittest.main()
