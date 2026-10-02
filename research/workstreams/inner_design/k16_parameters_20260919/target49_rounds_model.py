"""Two-round fixed-weight envelopes for the selected t64 subspace.

On a fixed nonzero state, P2 = P1/2 + J/2, where J is uniform nonzero
refresh. The zero source is unchanged. Each matrix is a positive majorant,
so combine whole source rows, not minima of incompatible source measures.
"""
from fractions import Fraction as F
from functools import lru_cache
import math
from flint import arb
import model
import subspace
import subspace_cover


class Engine(subspace.Engine):
    def identity(self):
        result=super().identity()
        result['inner']['transvection_rounds']=2
        result['setup']='independent row/region permutations; two independent transvections before feedback per update; zero start; output before update; persistent state; no flush'
        return result

    @lru_cache(maxsize=24)
    def epoch(self,tilt):
        old=super().epoch(tilt)
        z=(-model.number(tilt).exp()).exp()
        powers=[z**j for j in range(2*self.t+1)]
        result=[]
        for j,original in enumerate(old):
            total=math.comb(self.t,j)
            nk=arb(total-self.kernel[j])/total
            moments=[model.up(sum((math.comb(v,h)*math.comb(self.t-v,j-h)*powers[v+j-2*h]
                for h in range(max(0,j-self.t+v),min(v,j)+1)),arb(0))/total) for v in self.levels]
            new=list(original)
            for source,moment in [(1,max(moments))]+[(i+2,v) for i,v in enumerate(moments)]:
                fresh=[arb(0)]*self.n
                fresh[0]=min(model.up(moment),model.up(nk))/self.m
                for k,w in enumerate(self.levels): fresh[k+2]=moment*self.spectrum[w]/self.m
                for target in range(self.n):
                    index=source*self.n+target
                    new[index]=model.up((original[index]+fresh[target])/2)
            result.append(tuple(new))
        return result

    def bernoulli(self,theta,lam):
        # Conservative whole-row domination proved in target49_dense_bridge.py.
        original=super().bernoulli(theta,lam)
        factor=model.number(F(4*self.m+3,4*self.m+2))
        return original[:self.n]+tuple(model.up(factor*v) for v in original[self.n:])


class Checker(subspace_cover.Checker):
    def __init__(self,record):
        super().__init__(record)
        self.engine=Engine(record)
