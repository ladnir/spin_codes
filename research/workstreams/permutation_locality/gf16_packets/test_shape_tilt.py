from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest
from unittest.mock import patch
from flint import arb,arb_mat,ctx
import sparse_cover
from shape_tilt import normalizer,build_operators


def endpoint(value):
    m,e=value.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


class ShapeTiltTests(unittest.TestCase):
    def test_exact_change_of_measure(self):
        for rho,a in product((Q(1),Q(3,4),Q(1,2)),(Q(1),Q(3,4),Q(5,4))):
            c=normalizer(rho,a)
            weights=[a**x.bit_count()*rho**int(x==15) for x in range(1,16)]
            reference=[1/(15*c*w) for w in weights]
            self.assertEqual(sum(reference),1)
            self.assertEqual(c,sum(Q(comb(4,w),15)/(a**w*rho**int(w==4)) for w in range(1,5)))
            for q in (0,1,2):
                total=Q(0);tilted=Q(0)
                for xs in product(range(1,16),repeat=q):
                    value=Q(3,5)**(sum(xs)&15).bit_count()
                    mass=Q(1);penalty=Q(1)
                    for x in xs:mass*=reference[x-1];penalty*=weights[x-1]
                    total+=value/Q(15)**q;tilted+=mass*penalty*value
                self.assertEqual(total,c**q*tilted)

    def test_weight_shells_stay_uniform(self):
        for rho in (Q(1),Q(3,4),Q(1,2)):
            c=normalizer(rho)
            self.assertEqual(c,(14+1/rho)/15)
            for w in range(1,5):
                shell=[1/(15*c*rho**int(x==15)) for x in range(1,16) if x.bit_count()==w]
                self.assertTrue(all(p==shell[0] for p in shell))
                self.assertEqual(sum(shell),Q(comb(4,w),15)/(c*rho**int(w==4)))

    def test_invalid_weights(self):
        for rho,a in ((0,1),(-1,1),(2,1),(1,0),(1,-1)):
            with self.assertRaises(ValueError):normalizer(rho,a)

    def test_region_compensation_and_unchanged_outer_arguments(self):
        ctx.prec=192
        args=sparse_cover.build_args(3,['.064'],192,2,52,None)
        matrices=[arb_mat([[arb(1)/8,0],[arb(1)/4,arb(1)/2]]) for _ in range(4)]
        with patch.object(sparse_cover.sparse,'build_operators',return_value={('.064','.75'):(matrices,None)}) as build:
            exact,floating=build_operators(args,['.75'],'1')[('.064','.75')]
        self.assertEqual(args.penalties,['1'])
        self.assertEqual(args.weight_tilt,'1')
        self.assertEqual(build.call_args.args[0].penalties,['.75'])
        for j,(matrix,array) in enumerate(zip(exact,floating)):
            for row,col,value in ((0,0,Q(1,8)),(0,1,Q(0)),(1,0,Q(1,4)),(1,1,Q(1,2))):
                target=value*Q(46,45)**j;upper=endpoint(matrix[row,col])
                self.assertGreaterEqual(upper,target)
                self.assertLess(upper-target,Q(1,2)**180)
                self.assertEqual(array[row,col],float(matrix[row,col]))


if __name__=='__main__':unittest.main()
