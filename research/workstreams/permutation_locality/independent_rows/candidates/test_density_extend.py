"""Positive-measure domination and scaled convolution bounds."""
import unittest
import numpy as np
from flint import arb,arb_mat,ctx

from density_extend import rounded_counts,bounds,refine
from feedback_density import class_transforms
from occupancy_memory import C,U


class ExtendedDensity(unittest.TestCase):
    def test_count_rounding_is_pointwise_upward(self):
        for values in ([0,1,2,3],[(1<<60)+1,0,(1<<58)+123,1], [1,0,0,0]):
            counts=np.array(values,dtype=np.int64); denominator=sum(values)
            for bits in (2,8,40):
                small,total,unit=rounded_counts(counts,denominator,census_bits=bits)
                self.assertTrue(all(int(a)<=int(b)*unit for a,b in zip(counts,small)))
                self.assertGreaterEqual(total*unit,denominator)
                self.assertTrue(all(not a or b for a,b in zip(counts,small)))

    def test_small_convolution_encloses_exact_weighted_value(self):
        ctx.prec=192
        expansion=np.array([0,2,2,4,2,4,4,6],dtype=np.int64)
        prepared=class_transforms(expansion)
        counts=np.array([0,11,1,0,7,3,13,2],dtype=np.int64)
        denominator=sum(map(int,counts)); tilt='.17'; weights=(1,2)
        actual=bounds(counts,denominator,weights,expansion,prepared,[tilt],census_bits=2)[tilt]
        for target in range(1,len(counts)):
            density=arb(0)
            for v,size in prepared[1].items():
                mass=sum(int(counts[target^s]) for s in range(1,len(counts)) if expansion[s]==v)
                value=arb(mass)/(denominator*4)*(-arb(tilt)*abs(v-sum(weights))).exp()
                density+=value
                self.assertTrue(value/size<=actual['uniform'][v])
            self.assertTrue(density<=actual['density'])

    def test_rejects_nonmeasures_and_incomplete_shapes(self):
        for counts,total in (([1.,2.],3),([-1,3],2),([1,2],4),([1,1,1],3)):
            with self.assertRaises(ValueError): rounded_counts(np.array(counts),total)
        with self.assertRaises(ValueError): refine([None]*11,{}, {48:1},'.9')

    def test_shape_penalty_precedes_maximum_and_scope_is_scalar(self):
        spectrum={v:1 for v in (48,56,64,72,80)}
        base=[arb_mat([[10]*11 for _ in range(11)]) for _ in range(2)]
        records={(b,):{'density':arb(8 if b==4 else b),
                       'uniform':{v:arb(2*b) for v in spectrum}} for b in range(1,5)}
        result=refine(base,records,spectrum,'1/2',minimum=1,maximum=1)
        self.assertEqual(result[1][C,C],4)
        for source in range(U,U+5): self.assertEqual(result[1][source,C],6)
        self.assertEqual(base[1][C,C],10)
        self.assertEqual(result[1][C,0],10)
        self.assertEqual(result[0],base[0])


if __name__=='__main__':
    unittest.main()
