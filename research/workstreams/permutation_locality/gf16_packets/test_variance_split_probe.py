from fractions import Fraction as Q
from itertools import product
from math import factorial
import unittest
import numpy as np
import variance_split_probe as probe


class VarianceSplitProbeTests(unittest.TestCase):
    def test_rounded_outer_dual_dominates_exact_composition_counts(self):
        features=[Q(0),Q(1,4),Q(3,4),Q(1)];active=[0,1,1,1]
        coefficients=[Q(1),Q(2),Q(3),Q(1)];n=4;mean=Q(1,2);occupancy=Q(1,2)
        logs=np.log(list(map(float,coefficients)))
        for interval in ((Q(0),Q(1,16)),(Q(1,16),Q(1,8)),(Q(1,8),Q(1,4))):
            value,dual=probe.outer_witness(logs,features,active,float(mean),float(occupancy),interval)
            eta,mu,gamma=map(float,dual);self.assertGreaterEqual(mu,0)
            total=Q(0)
            for counts in product(range(n+1),repeat=4):
                if sum(counts)!=n or sum(a*c for a,c in zip(active,counts))<n*occupancy:continue
                if sum(f*c for f,c in zip(features,counts))!=n*mean:continue
                variance=sum(f*(1-f)*c for f,c in zip(features,counts))/n
                if not interval[0]<=variance<=interval[1]:continue
                term=Q(factorial(n))
                for count,coefficient in zip(counts,coefficients):term*=coefficient**count/factorial(count)
                total+=term
            if total:self.assertGreaterEqual(n*value+1e-9,np.log(float(total)))
            # The witness formula remains valid regardless of optimizer
            # quality; these are proposals, not rounded certificate data.
            exponents=[eta*float(f)+mu*a+gamma*float(f*(1-f)) for f,a in zip(features,active)]
            direct=probe.logsumexp(logs+np.array(exponents))-eta*float(mean)-mu*float(occupancy)-min(gamma*float(v) for v in interval)
            self.assertAlmostEqual(value,direct,places=10)


if __name__=='__main__':unittest.main()
