import math
import unittest
from fractions import Fraction as F

import activation_occupation as baseline
import random_spectrum_variance as variance


class RandomSpectrumVarianceTest(unittest.TestCase):
    def test_exact_moments_over_all_binary_four_by_two_subspaces(self):
        codes={frozenset((0,a,b,a^b)) for a in range(1,16) for b in range(1,16) if a!=b}
        self.assertEqual(len(codes),35)
        for w in range(1,5):
            samples=[sum(x.bit_count()==w for x in code) for code in codes]
            mean=F(sum(samples),len(samples))
            var=sum((F(x)-mean)**2 for x in samples)/len(samples)
            self.assertEqual((mean,var),variance.shell_moments(4,2,w))

    def test_simultaneous_failure_charge_and_markov_domination(self):
        for b,d,h in ((16,8,0),(128,64,20),(512,256,60),(1024,512,60)):
            caps=variance.caps(b,d,h)
            old=baseline.random_spectrum_caps(b,d,h)
            charge=F(0)
            for w in range(1,b+1):
                cap=caps.get(w,0)
                self.assertLessEqual(cap,old.get(w,0))
                if cap>=min((1<<d)-1,math.comb(b,w)):continue
                mean,var=variance.shell_moments(b,d,w)
                probability=mean/(cap+1)
                if cap+1>mean:probability=min(probability,var/(cap+1-mean)**2)
                charge+=probability
            self.assertLessEqual(charge,F(1,1<<h))

    def test_large_central_shell_is_close_to_its_mean(self):
        mean,_=variance.shell_moments(1024,512,512)
        cap=variance.caps(1024,512,60)[512]
        self.assertLess(cap,mean*F(1000001,1000000))
        self.assertGreaterEqual(cap,mean)


if __name__=='__main__':unittest.main()
