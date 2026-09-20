"""Exact rational domination checks for the independent-map outward code."""
from collections import Counter
from fractions import Fraction as F
import math
import unittest
from unittest.mock import patch

from flint import arb,ctx
import certify as cert
import search_feedback as search


def fraction(value):
    mantissa,exponent=map(int,value.upper().man_exp())
    return F(mantissa)*(F(2)**exponent)


class OutwardTests(unittest.TestCase):
    def setUp(self):
        ctx.prec=256
        self.a=[1,2,3,4,5];self.b=[3,5,6,7,1];self.t=5;self.s=3
        self.weights=[search.image(self.a,q).bit_count() for q in range(8)]
        spectrum=Counter(self.weights);del spectrum[0]
        dual=Counter(search.image(self.b,q).bit_count() for q in range(8));del dual[0]
        kernel=[0]*6
        for x in range(32):kernel[x.bit_count()]+=int(search.inject(self.b,x)==0)
        self.engine=object.__new__(cert.Engine)
        self.engine.a_columns=self.a;self.engine.columns=self.b
        self.engine.spectrum=dict(spectrum);self.engine.b_spectrum=dict(dual)
        self.engine.kernel=kernel;self.engine.levels=sorted(spectrum)
        self.engine.n=len(spectrum)+2;self.engine.m=7
        self.engine.caps=search.g.fiber_caps(5,3,dict(dual),kernel)
        self.engine.low=search.low_cancellation(self.a,self.b,3)

    def check(self,matrix,probabilities,z):
        engine=self.engine;n=engine.n
        bound=[[fraction(matrix[i*n+j]) for j in range(n)] for i in range(n)]
        exact=[[F(0) for _ in range(8)] for _ in range(8)]
        for x,p in enumerate(probabilities):
            syndrome=search.inject(self.b,x)
            for q in range(8):
                mass=p*z**((x^search.image(self.a,q)).bit_count())
                if q==0:exact[q][syndrome]+=mass
                else:
                    exact[q][q^syndrome]+=mass/2
                    for mixed in range(1,8):exact[q][mixed^syndrome]+=mass/14
        laws=[([F(q==source) for q in range(8)],int(source!=0)) for source in range(8)]
        laws += [([F(self.weights[q]==v,engine.spectrum[v]) for q in range(8)],i+2) for i,v in enumerate(engine.levels)]
        for law,index in laws:
            actual=[sum(law[q]*exact[q][y] for q in range(8)) for y in range(8)]
            self.assertLessEqual(actual[0],bound[index][0])
            caps={v:bound[index][i+2]/engine.spectrum[v] for i,v in enumerate(engine.levels)}
            residual=sum(max(F(0),actual[q]-caps[self.weights[q]]) for q in range(1,8))
            self.assertLessEqual(residual,bound[index][1])

    def test_fixed_weight_arb_envelope(self):
        with patch.object(cert.base,'T',5):
            for log_lam in (F(-1),F(-3)):
                # The exact test uses the rational lower endpoint of z.
                # The producer's z enclosure lies above that endpoint.
                rows=self.engine.epoch(log_lam)
                z=fraction((-cert.old.number(log_lam).exp()).exp().lower())
                for j in range(6):
                    probabilities=[F(x.bit_count()==j,math.comb(5,j)) for x in range(32)]
                    self.check(rows[j],probabilities,z)

    def test_bernoulli_arb_envelope(self):
        with patch.object(cert,'T',5):
            for theta in (F(1,10),F(1,2),F(9,10)):
                for z in (F(1,2),F(9,10)):
                    lam=-(arb(z.numerator)/z.denominator).log()
                    matrix=self.engine.bernoulli(arb(theta.numerator)/theta.denominator,lam)
                    probabilities=[theta**x.bit_count()*(1-theta)**(5-x.bit_count()) for x in range(32)]
                    self.check(matrix,probabilities,z)


if __name__=='__main__':unittest.main()
