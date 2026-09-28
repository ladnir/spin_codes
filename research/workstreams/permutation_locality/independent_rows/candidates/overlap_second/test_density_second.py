"""Direct small-state comparisons for every translated quadratic bound."""
from fractions import Fraction as Q
import unittest

import numpy as np
from flint import arb,ctx

from second import moments
from density_second import shape_bounds,dyadic_intercepts
from feedback_density import class_transforms
from test_mean import linear


class QuadraticDensity(unittest.TestCase):
    def setUp(self):
        self.old=ctx.prec;ctx.prec=192
        self.E=(0b000000001111,0b111111000000,0b010101010101)
        self.B=(1,2,3,4,5,6,7,1,3,5,7,2)
        self.images=[linear(self.E,s) for s in range(8)]
        self.expansion=np.array([x.bit_count() for x in self.images],dtype=np.int64)
        self.prepared=class_transforms(self.expansion)
        self.second=moments(self.E,self.B,2,width=2)
        self.words={s:[] for s in self.second}
        for x in range(1,1<<len(self.B)):
            shape=tuple(sorted(b for b in ((x>>(2*w)&3).bit_count() for w in range(len(self.B)//2)) if b))
            if shape in self.words:self.words[shape].append(x)

    def tearDown(self):ctx.prec=self.old

    def test_all_shapes_targets_classes_and_upward_count_rounding(self):
        for shape,words in self.words.items():
            counts=np.array([sum(linear(self.B,x)==s for x in words) for s in range(8)],dtype=np.int64)
            D,_,second=self.second[shape];self.assertEqual(D,len(words))
            for census_bits,multiplier in ((40,1),(2,1),(2,(1<<45)+1)):
                ctx.prec=192
                bounds=shape_bounds(counts*multiplier,D*multiplier,shape,self.expansion,self.prepared,
                                    second,['.01','.17'],census_bits=census_bits)
                ctx.prec=512
                for tilt,row in bounds.items():
                    for t in range(1,8):
                        by_class={v:arb(0) for v in self.prepared[0]}
                        for x in words:
                            s=t^linear(self.B,x)
                            if not s:continue
                            v=int(self.expansion[s])
                            by_class[v]+=(-arb(tilt)*(x^self.images[s]).bit_count()).exp()/(4*D)
                        self.assertLessEqual(sum(by_class.values(),arb(0)),row['density'])
                        for v,mass in by_class.items():
                            self.assertLessEqual(mass/self.prepared[1][v],row['uniform'][v])
                    if census_bits==2:self.assertGreater(row['rounding_factor'],1)

    def test_dyadic_constraints_for_every_overlap_and_slope(self):
        for W in (3,4,20,40):
            for c in (Q(0),Q(1,1000),Q(1,2)):
                ctx.prec=192
                coefficients=dyadic_intercepts((48,56,64),W,'.072',c,19)
                ctx.prec=512
                for v,a in coefficients.items():
                    for h in range(W+1):
                        d=2*h-W
                        self.assertGreaterEqual(arb(a)/(1<<19)+arb(c.numerator)/c.denominator*d*d,
                                                (-arb('.072')*(v-d)).exp())

    def test_reject_unsupported_parameters(self):
        counts=np.ones(8,dtype=np.int64)
        for shape,second in (((4,4),Q(8)),((1,),.5),((1,),Q(2)),((1,),Q(-1))):
            with self.assertRaises(ValueError):
                shape_bounds(counts,8,shape,self.expansion,self.prepared,second,['.17'])


if __name__=='__main__':unittest.main()
