"""Exhaustive short-word checks of the exact overlap and chord bound."""
from collections import Counter
from fractions import Fraction as Q
import unittest
from itertools import product
from unittest.mock import patch
from flint import arb,arb_mat,ctx

from mean import moments,chord_histogram,zero_refine,C,Z


def linear(columns,x):
    y=0
    for i,c in enumerate(columns):
        if x>>i&1:y^=c
    return y


class OverlapMean(unittest.TestCase):
    def test_all_short_shapes_offsets_and_output_tilts(self):
        checked=0; strict=0
        cases=[((0b101101,0b011110,0b110001),(1,2,3,4,5,6)),
               ((0b10101010,0b11001100,0b11110000),(1,3,5,7,2,4,6,1)),
               ((0b011001,0b101010,0b110100),(0,0,0,1,2,3))]
        for expansion,feedback in cases:
            n=len(feedback)//2
            records=moments(expansion,feedback,3,width=2)
            groups={}
            for x in range(1<<(2*n)):
                shape=tuple(sorted((x>>(2*w)&3).bit_count() for w in range(n) if x>>(2*w)&3))
                if 1<=len(shape)<=3: groups.setdefault(shape,[]).append(x)
            self.assertEqual(set(records),set(groups))
            for shape,words in groups.items():
                D,mean,universal=records[shape]; W=sum(shape)
                self.assertEqual(D,len(words))
                self.assertEqual(mean,Q(sum((x&linear(expansion,linear(feedback,x))).bit_count() for x in words),D))
                for target in range(1<<len(expansion)):
                    translated=[(x,linear(expansion,target^linear(feedback,x))) for x in words]
                    actual_mean=Q(sum((x&y).bit_count() for x,y in translated),D)
                    self.assertLessEqual(actual_mean,universal)
                    counts=Counter(y.bit_count() for x,y in translated if y)
                    classes={v:Q(c,D) for v,c in counts.items()}
                    histogram=chord_histogram(classes,W,mean if target==0 else universal)
                    for z in (Q(1,2),Q(4,5),Q(1)):
                        actual=sum((z**(y^x).bit_count() for x,y in translated if y),Q(0))/D
                        bound=sum((p*z**w for w,p in histogram.items()),Q(0))
                        triangle=sum((p*z**(v-W) for v,p in classes.items()),Q(0))
                        self.assertLessEqual(actual,bound)
                        self.assertLessEqual(bound,triangle)
                        checked+=1;strict+=bound<triangle
        self.assertEqual(checked,648)
        self.assertEqual(strict,290)

    def test_allocation_endpoints_and_scope(self):
        classes={48:Q(1,8),56:Q(3,8),64:Q(1,2)}
        self.assertEqual(chord_histogram(classes,20,Q(0)),{68:Q(1,8),76:Q(3,8),84:Q(1,2)})
        self.assertEqual(chord_histogram(classes,20,Q(20)),{28:Q(1,8),36:Q(3,8),44:Q(1,2)})
        self.assertEqual(chord_histogram(classes,20,Q(10)),{28:Q(1,8),36:Q(3,8),84:Q(1,2)})
        self.assertEqual(chord_histogram({},20,Q(10)),{})

    def test_bad_geometry_and_nonexact_parameters(self):
        for width,maximum in ((0,1),(2,0),(2,3)):
            with self.assertRaises(ValueError): moments((1,2),(1,2,1,2),maximum,width=width)
        for classes,W,mean in (({48:Q(2)},20,Q(10)),({0:Q(1)},20,Q(10)),
                               ({48:Q(1)},0,Q(0)),({48:Q(1)},20,.5),({48:Q(1)},20,Q(21))):
            with self.assertRaises(ValueError): chord_histogram(classes,W,mean)

    def test_shape_maximum_penalty_and_outward_refinement_scope(self):
        saved_precision=ctx.prec
        self.addCleanup(setattr,ctx,'prec',saved_precision)
        ctx.prec=192
        expansion=(0b10110101,0b01111010,0b11000110)
        feedback=(1,2,3,4,5,6,7,2)
        overlap=moments(expansion,feedback,2,width=4)
        words={s:[] for s in overlap}
        for x in range(1,256):
            shape=tuple(sorted(b for b in ((x&15).bit_count(),(x>>4).bit_count()) if b))
            words[shape].append(x)
        census={}
        for shape,group in words.items():
            counts=Counter(linear(expansion,linear(feedback,x)).bit_count() for x in group)
            zero=counts.pop(0,0)
            census[shape]=zero,0,len(group),counts
        base=[arb_mat([[1]*11 for _ in range(11)]) for _ in range(3)]
        before=[t*1 for t in base]
        with patch('builtins.print'):
            after=zero_refine(base,overlap,census,'.1','3/4',minimum=1,maximum=2)
        ctx.prec=512
        for j in range(1,3):
            expected=arb(0)
            actual=arb(0)
            for shape,group in words.items():
                if len(shape)!=j: continue
                z,p,D,classes=census[shape]
                hist=chord_histogram({v:Q(n,D) for v,n in classes.items()},sum(shape),overlap[shape][1])
                scale=Q(3,4)**shape.count(4)/4
                value=sum((arb(p.numerator)/p.denominator*(-arb('.1')*w).exp() for w,p in hist.items()),arb(0))
                expected=max(expected,value*scale.numerator/scale.denominator)
                true=sum(((-arb('.1')*(x^linear(expansion,linear(feedback,x))).bit_count()).exp()
                          for x in group if linear(feedback,x)),arb(0))/D
                actual=max(actual,true*scale.numerator/scale.denominator)
            self.assertGreaterEqual(after[j][C,Z],expected)
            self.assertGreaterEqual(after[j][C,Z],actual)
            self.assertAlmostEqual(float(after[j][C,Z]),float(expected),places=13)
        self.assertEqual(base,before)
        for j,s,t in product(range(3),range(11),range(11)):
            if j==0 or (s,t)!=(C,Z):self.assertEqual(after[j][s,t],before[j][s,t])
        del census[(4,4)]
        with self.assertRaisesRegex(ValueError,'complete'):
            zero_refine(base,overlap,census,'.1','3/4',minimum=1,maximum=2)


if __name__=='__main__':
    unittest.main()
