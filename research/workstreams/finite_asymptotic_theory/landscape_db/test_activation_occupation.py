"""Exact finite checks including nonzero kernel inputs and setup caps."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np

import activation_occupation as general
from test_activation_q1 import mul, identity


def exact_epochs(z):
    result=[]
    for j in range(9):
        moments=[]
        for w,count in ((4,14),(8,1)):
            moment=sum((F(math.comb(w,v)*math.comb(8-w,j-v),math.comb(8,j))*z**(w+j-2*v)
                        for v in range(max(0,j-8+w),min(w,j)+1)),F(0))
            moments.append((moment,count))
        arbitrary=max(v for v,c in moments)
        uniform=sum(v*c for v,c in moments)/15
        beta=F({0:1,4:14,8:1}.get(j,0),math.comb(8,j))
        result.append([[beta*z**j,(1-beta)*z**j,F(0)],
                       [min(1-beta,arbitrary)/15,F(0),arbitrary],
                       [min(1-beta,uniform)/14,F(0),F(15,14)*uniform]])
    return result


class GeneralOccupationTest(unittest.TestCase):
    def test_all_epoch_weights_and_region_coefficients(self):
        z=F(3,4)
        exact=exact_epochs(z)
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],-math.log(float(z)),8)
        np.testing.assert_allclose(np.exp(epoch),np.array(exact,dtype=float),rtol=2e-14,atol=2e-15)
        regions=general.region_logs(epoch,8,16,16)
        for j in range(17):
            total=[[F(0)]*3 for _ in range(3)]
            for a in range(max(0,j-8),min(8,j)+1):
                product=mul(exact[a],exact[j-a])
                count=F(math.comb(8,a)*math.comb(8,j-a),math.comb(16,j))
                total=[[x+count*y for x,y in zip(row,p)] for row,p in zip(total,product)]
            np.testing.assert_allclose(np.exp(regions[j]),np.array(total,dtype=float),rtol=2e-13,atol=2e-14)
        # Truncation must leave low coefficients unchanged.
        short=general.region_logs(epoch[:5],8,16,4)
        np.testing.assert_allclose(short,regions[:5],rtol=1e-13,atol=1e-13)

    def test_actual_gf16_paths_include_kernel_weights_four_and_eight(self):
        generators=[255,170,204,240]
        words=[0]*16
        for state in range(16):
            for i,g in enumerate(generators):
                if state>>i&1: words[state]^=g
        def field_mul(a,b):
            out=0
            while b:
                if b&1: out^=a
                b>>=1; a<<=1
                if a&16: a^=19
            return out
        z=F(3,4)
        actual=[]
        for weight in range(9):
            matrix=[[F(0)]*16 for _ in range(16)]
            for selected in itertools.combinations(range(8),weight):
                x=sum(1<<j for j in selected)
                syndrome=sum(((x&g).bit_count()%2)<<i for i,g in enumerate(generators))
                for state in range(16):
                    factor=z**((x^words[state]).bit_count())/(15*math.comb(8,weight))
                    for alpha in range(1,16):
                        matrix[state][field_mul(alpha,state)^syndrome]+=factor
            actual.append(matrix)
        envelope=exact_epochs(z)
        for a,b in itertools.product(range(9),repeat=2):
            tail=[sum(row) for row in actual[b]]
            exact=[sum(x*y for x,y in zip(row,tail)) for row in actual[a]]
            bound=[sum(row) for row in mul(envelope[a],envelope[b])]
            self.assertLessEqual(exact[0],bound[0])
            self.assertLessEqual(max(exact[1:]),bound[1])
            # Extremal density-bounded laws put mass 1/14 on any 14 states.
            self.assertLessEqual((sum(exact[1:])-min(exact[1:]))/14,bound[2])

    def test_adaptive_recurrence_dominates_all_band_assignments(self):
        rng=np.random.default_rng(89123)
        regions=rng.uniform(.001,.9,size=(5,3,3))
        ps=[F(1,4),F(3,4),F(1)]
        roots=[2,3,1]
        active=np.log(np.array(ps,dtype=float))
        inactive=np.array([math.log(float(1-p)) if p!=1 else -np.inf for p in ps])
        results=list(general.adaptive_logs(np.log(regions),np.log(roots),active,inactive))
        for q,actual in results:
            for assignment in itertools.product(range(3),repeat=q):
                probabilities=[F(1)]
                factor=1
                for g in assignment:
                    p=ps[g]; factor*=roots[g]
                    updated=[F(0)]*(len(probabilities)+1)
                    for j,v in enumerate(probabilities):
                        updated[j]+=v*(1-p); updated[j+1]+=v*p
                    probabilities=updated
                expected=factor*np.sum(np.array(probabilities,dtype=float)[:,None,None]*regions[:q+1],axis=0)
                self.assertTrue(np.all(np.exp(actual)+1e-11 >= expected))

    def test_density_domination_at_every_supported_weight(self):
        counts={1:3,2:4,4:1}
        bands=general.bands_for(counts,4,2)
        roots,p,not_p=general.density_roots(counts,4,bands,shift=.4)
        for index,band in enumerate(bands):
            for w in band:
                auxiliary=roots[index]*4+w*p[index]
                if w<4: auxiliary+=(4-w)*not_p[index]
                self.assertGreaterEqual(auxiliary+1e-13,math.log(counts[w])-math.log(math.comb(4,w)))

    def test_random_caps_charge_one_simultaneous_failure_event(self):
        for b,d,h in ((16,8,0),(64,32,8),(1024,512,60)):
            caps=general.random_spectrum_caps(b,d,h)
            charge=F(0)
            for w in range(1,b+1):
                cap=caps.get(w,0)
                if cap < min((1<<d)-1,math.comb(b,w)):
                    expectation=F(math.comb(b,w)*((1<<d)-1),(1<<b)-1)
                    charge+=expectation/(cap+1)
            self.assertLessEqual(charge,F(1,1<<h))
        self.assertNotIn(1,general.random_spectrum_caps(1024,512,60))


if __name__ == '__main__':
    unittest.main()
