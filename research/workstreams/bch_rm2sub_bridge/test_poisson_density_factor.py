from fractions import Fraction as F
import itertools
import math
import unittest
from poisson_density_factor import density_factor


class DensityTest(unittest.TestCase):
    def test_mode_mass(self):
        for n in range(1,65):
            for j in range(n+1):
                p=F(j,n);mass=math.comb(n,j)*p**j*(1-p)**(n-j)
                self.assertGreaterEqual(mass,F(1,density_factor(n)))

    def test_small_poisson_binomial_distributions(self):
        for n in range(1,6):
            for ps in itertools.product((F(0),F(1,4),F(3,4),F(1)),repeat=n):
                dist=[F(1)]
                for p in ps:
                    new=[F(0)]*(len(dist)+1)
                    for j,value in enumerate(dist):new[j]+=value*(1-p);new[j+1]+=value*p
                    dist=new
                r=sum(ps,F(0))/n
                for j,value in enumerate(dist):
                    self.assertLessEqual(value,density_factor(n)*math.comb(n,j)*r**j*(1-r)**(n-j))

    def test_large_length_factor(self):
        self.assertEqual(density_factor(1<<19),2900)
        self.assertEqual(density_factor(1<<21),5796)
        for n in (1,2,3,4,17,1<<19,1<<21):
            self.assertLessEqual(density_factor(n),n+1)


if __name__=='__main__':unittest.main()
