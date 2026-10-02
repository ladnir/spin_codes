from fractions import Fraction as Q
from math import comb
import unittest
import numpy as np
from flint import arb,ctx
import kernel


def endpoint(x):
    m,e=x.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


class KernelTests(unittest.TestCase):
    def test_every_state_against_exhaustive_eight_bit_input(self):
        ctx.prec=192
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        columns=[1,2,4,3,5,7,6,1]
        for updates in (1,2,3):
            data=kernel.prepare(images,columns,3,updates)
            for p in ([Q(1,5)]*5,[Q(1,2),Q(1,4),Q(1,4),0,0],[0,0,0,0,1],[1,0,0,0,0]):
                z=Q(3,4);matrix=kernel.outward_at_z(data,p,arb(3)/4)
                for state in range(8):
                    actual=[Q(0),Q(0)]
                    for x in range(256):
                        probability=Q(1)
                        for offset in (0,4):
                            weight=((x>>offset)&15).bit_count()
                            probability*=p[weight]/Q(comb(4,weight))
                        feedback=0
                        for b,c in enumerate(columns):
                            if x>>b&1:feedback^=c
                        weight=probability*z**(images[state]^x).bit_count()
                        if state==0:actual[int(feedback!=0)]+=weight
                        else:
                            for refreshed in range(1,8):
                                chance=Q(int(refreshed==state),2**updates)+(1-Q(1,2**updates))/7
                                actual[int(refreshed!=feedback)]+=weight*chance
                    for j in range(2):
                        self.assertGreaterEqual(endpoint(matrix[int(state!=0),j]),actual[j])

    def test_float_matches_outward_kernel(self):
        ctx.prec=192
        images=[s for s in range(8)];columns=[1,2,4,3]
        data=kernel.prepare(images,columns,3)
        probabilities=[Q(1,5)]*5
        exact=kernel.outward(data,probabilities,Q(1,8))
        proposed=kernel.floating(data,probabilities,.125)
        np.testing.assert_allclose(proposed,np.array([[float(exact[i,j]) for j in range(2)] for i in range(2)]),rtol=2e-14)


if __name__=='__main__':unittest.main()
