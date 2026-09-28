from fractions import Fraction as Q
from itertools import product
import unittest

import numpy as np
from flint import arb,arb_mat,ctx

import iid_kernel
from test_model import linear


class IidKernelTests(unittest.TestCase):
    def test_tilted_feedback_atom_for_every_small_state_word(self):
        ctx.prec=192
        feedback=[1,2,3,5,4,7];images=[linear([0b010101,0b101010,0b110001],s) for s in range(8)]
        data=iid_kernel.prepare(images,feedback,3)
        v,r,tilt=Q(3,5),Q(2,7),Q(3,2)
        atom=iid_kernel.outward_tilted_atom(data,v,r,tilt)
        self.assertAlmostEqual(float(atom),iid_kernel.floating_tilted_atom(data,float(v),float(r),float(tilt)),places=13)
        probabilities=[1-v,v*(1-r)/2,v*(1-r)/2,v*r];z=(-arb(3)/2).exp()
        for y in range(64):
            bins=[arb(0) for _ in range(8)]
            for packets in product(range(4),repeat=3):
                x=sum(p<<(2*i) for i,p in enumerate(packets));mass=Q(1)
                for p in packets:mass*=probabilities[p]
                bins[linear(feedback,x)]+=(arb(mass.numerator)/mass.denominator)*z**(x^y).bit_count()
            total=sum(bins,arb(0))
            for value in bins:self.assertTrue(value.lower()<=(atom*total).upper())

    def test_exact_chain_is_dominated(self):
        ctx.prec=192
        expansion=[0b010101,0b101010,0b110001];feedback=[1,2,3,5,4,7]
        images=[linear(expansion,s) for s in range(8)];data=iid_kernel.prepare(images,feedback,3)
        cases=((2,Q(3,4),Q(1,3),Q(2)),(2,Q(1,10),Q(4,5),Q(1,10)),(2,Q(9,10),Q(9,10),Q(3,2)),
               (1,Q(1,10),Q(4,5),Q(1,10)),(3,Q(1,10),Q(4,5),Q(1,10)),(4,Q(3,4),Q(1,3),Q(2)))
        for updates,v,r,tilt in cases:
            data={**data,'updates':updates};alpha=arb(2)**-updates
            envelope=iid_kernel.outward(data,v,r,tilt)
            floating=iid_kernel.floating(data,float(v),float(r),float(tilt))
            np.testing.assert_allclose(floating,[[float(envelope[i,j]) for j in range(2)] for i in range(2)],rtol=2e-12,atol=1e-15)
            probabilities=[1-v,v*(1-r)/2,v*(1-r)/2,v*r]
            z=(-arb(tilt.numerator)/tilt.denominator).exp();actual=arb_mat(8,8)
            for s in range(8):
                for packets in product(range(4),repeat=3):
                    x=sum(p<<(2*i) for i,p in enumerate(packets));mass=Q(1)
                    for p in packets:mass*=probabilities[p]
                    weight=(arb(mass.numerator)/mass.denominator)*z**(x^images[s]).bit_count()
                    c=linear(feedback,x)
                    if s==0:actual[0,c]+=weight
                    else:
                        actual[s,s^c]+=alpha*weight
                        for fresh in range(1,8):actual[s,fresh^c]+=(1-alpha)*weight/7
            for s in range(8):
                row=int(s!=0)
                self.assertTrue(actual[s,0].lower()<=envelope[row,0])
                self.assertTrue(sum((actual[s,t] for t in range(1,8)),arb(0)).lower()<=envelope[row,1])
            for length in (1,2,3,7):
                exact=actual**length;upper=envelope**length
                self.assertTrue(sum((exact[0,t] for t in range(8)),arb(0)).lower()<=sum((upper[0,t] for t in range(2)),arb(0)).upper())


if __name__=='__main__':unittest.main()
