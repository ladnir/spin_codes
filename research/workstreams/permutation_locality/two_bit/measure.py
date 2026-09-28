"""Pointwise domination of the shuffled pair's support measure.

If S_u bounds the number of outer pairs with union size u, the mass on
each specific support is at most S_u/C(n,u). Against iid Bernoulli(p)
supports its density is at most max_u S_u/Pr[Bin(n,p)=u]. This bounds
the *whole* interval simultaneously, not one conditional event at a time.
CDF caps may cap an individual shell; their differences may not.
"""
from math import comb,log

import numpy as np
from flint import arb

import model
from shell_cover import IntervalFolds


class DensityFolds(IntervalFolds):
    def __init__(self,cdf,shells,**kwargs):
        super().__init__(cdf,shells,**kwargs)
        self.caps=tuple(min(a,b) for a,b in zip(self.cdf,self.shells))
        self.density_functions={}

    def function(self,lo,hi):
        key=lo,hi
        if key not in self.density_functions:
            old=super().function(lo,hi)
            supports=np.arange(lo,hi+1)
            values=np.array([log(self.caps[u].numerator)-log(self.caps[u].denominator)-log(comb(self.n,u))
                             if self.caps[u] else -np.inf for u in supports])
            def evaluate(p):
                if not 0<p<1:raise ValueError('interior support witness required')
                return min(old(p),float(np.max(values-supports*np.log(p)-(self.n-supports)*np.log1p(-p))))
            self.density_functions[key]=evaluate
        return self.density_functions[key]

    def outward(self,lo,hi,p):
        if not 0<=lo<=hi<=self.n or not 0<p<1:raise ValueError('bounded support interval and interior witness required')
        def term(u):
            c=self.caps[u]
            return model.up((arb(c.numerator)/c.denominator)/(comb(self.n,u)*p**u*(1-p)**(self.n-u)))
        return min(super().outward(lo,hi,p),max(term(u) for u in range(lo,hi+1)))
