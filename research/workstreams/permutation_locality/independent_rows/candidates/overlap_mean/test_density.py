"""Compare integer translated-chord evaluation with complete small laws."""
from collections import Counter
from fractions import Fraction as Q
import unittest

import numpy as np
from flint import arb,ctx

from mean import moments,chord_histogram
from density import shape_bounds,powers
from feedback_density import class_transforms
from test_mean import linear


class OverlapDensity(unittest.TestCase):
    def setUp(self):
        self.old=ctx.prec;ctx.prec=192
        self.E=(0b000000001111,0b111111000000,0b010101010101)
        self.B=(1,2,3,4,5,6,7,1,3,5,7,2)
        self.images=[linear(self.E,s) for s in range(8)]
        self.expansion=np.array([x.bit_count() for x in self.images],dtype=np.int64)
        self.prepared=class_transforms(self.expansion)
        self.overlap=moments(self.E,self.B,2,width=2)
        self.words={s:[] for s in self.overlap}
        for x in range(1,1<<len(self.B)):
            shape=tuple(sorted(b for b in ((x>>(2*w)&3).bit_count() for w in range(len(self.B)//2)) if b))
            if shape in self.words:self.words[shape].append(x)

    def tearDown(self):ctx.prec=self.old

    def test_every_shape_target_and_uniform_class_against_direct_sum(self):
        for shape,words in self.words.items():
            counts=np.zeros(len(self.images),dtype=np.int64)
            for x in words: counts[linear(self.B,x)]+=1
            D,_,mu=self.overlap[shape]
            self.assertEqual(D,len(words))
            for bits in (2,40):
                ctx.prec=192
                values=shape_bounds(counts,D,shape,self.expansion,self.prepared,mu,['.17'],census_bits=bits)['.17']
                ctx.prec=512
                for t in range(1,len(self.images)):
                    actual=arb(0)
                    by_class={v:arb(0) for v in self.prepared[0]}
                    class_counts=Counter()
                    for x in words:
                        s=t^linear(self.B,x)
                        if not s: continue
                        v=self.expansion[s].item(); class_counts[v]+=1
                        term=(-arb('.17')*(x^self.images[s]).bit_count()).exp()/(4*D)
                        actual+=term;by_class[v]+=term
                    self.assertLessEqual(actual,values['density'])
                    for v,mass in by_class.items():
                        self.assertLessEqual(mass/self.prepared[1][v],values['uniform'][v])
                    if bits==40:
                        # With no count rounding, integer budget rounding
                        # adds at most one unit to the exact overlap budget.
                        budget=mu*D; raised=Q((budget.numerator+budget.denominator-1)//budget.denominator,D)
                        histogram=chord_histogram({v:Q(n,D) for v,n in class_counts.items()},sum(shape),raised)
                        exact=sum((arb(p.numerator)/p.denominator*(-arb('.17')*w).exp()
                                   for w,p in histogram.items()),arb(0))/4
                        self.assertLessEqual(exact,values['density'])

    def test_upward_counts_keep_the_unnormalized_factor(self):
        shape=(1,1);words=self.words[shape]
        counts=np.array([sum(linear(self.B,x)==s for x in words) for s in range(len(self.images))],dtype=np.int64)
        D=len(words);mu=self.overlap[shape][2]
        # Scaling the counting measure leaves its law and overlap unchanged.
        multiplier=(1<<45)+1
        values=shape_bounds(counts*multiplier,D*multiplier,shape,self.expansion,self.prepared,mu,
                            ['.01','.17'],census_bits=2)
        self.assertTrue(all(row['rounding_factor']>1 for row in values.values()))
        for tilt,row in values.items():
            for t in range(1,len(self.images)):
                actual=sum(((-arb(tilt)*(x^self.images[t^linear(self.B,x)]).bit_count()).exp()
                            for x in words if t^linear(self.B,x)),arb(0))/(4*D)
                self.assertLessEqual(actual,row['density'])

    def test_rounding_powers_is_upward_and_ordered(self):
        values=powers((48,56,64,72,80),28,'.072',17)
        for v,(low,high) in values.items():
            self.assertLessEqual(low,high)
            self.assertGreaterEqual(arb(low)/2**17,(-arb('.072')*(v+28)).exp())
            self.assertGreaterEqual(arb(high)/2**17,(-arb('.072')*(v-28)).exp())

    def test_reject_unsupported_weight_and_invalid_cap(self):
        counts=np.ones(len(self.images),dtype=np.int64)
        for weights,mean in (((4,4),Q(4)),((1,),.5),((1,),Q(2))):
            with self.assertRaises(ValueError):
                shape_bounds(counts,len(self.images),weights,self.expansion,self.prepared,mean,['.17'])


if __name__=='__main__':unittest.main()
