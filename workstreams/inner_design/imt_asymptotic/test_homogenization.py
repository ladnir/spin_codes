"""Small-state checks of the analytic IMT region homogenization identities."""
from fractions import Fraction as F
import math
import unittest

import numpy as np


def exact_impulse(a,b):
    size=1<<max(c.bit_length() for c in a+b)
    matrix=[[F(0) for _ in range(size)] for _ in range(size)]
    for q in range(size):
        for column in b:
            if q==0:
                matrix[q][column]+=F(1,len(b))
            else:
                matrix[q][q^column]+=F(1,2*len(b))
                for fresh in range(1,size):
                    matrix[q][fresh^column]+=F(1,2*(size-1)*len(b))
    return matrix


def kernels(a,b,z):
    size=1<<max(c.bit_length() for c in a+b)
    zero=np.zeros((size,size));one=zero.copy()
    for q in range(size):
        image=sum(((column&q).bit_count()%2)<<j for j,column in enumerate(a))
        law=np.zeros(size)
        if q==0:law[0]=1
        else:
            law[1:]=.5/(size-1);law[q]+=.5
        zero[q]=z**image.bit_count()*law
        for j,column in enumerate(b):
            moment=z**(image^(1<<j)).bit_count()/len(b)
            for r,p in enumerate(law):one[q,r^column]+=moment*p
    return zero,one


class HomogenizationTests(unittest.TestCase):
    a=[1,2,4,7]
    b=[1,3,5,7]

    def test_impulse_projection_exact(self):
        matrix=exact_impulse(self.a,self.b)
        size=len(matrix);m=size-1
        self.assertEqual(matrix[0][0],0)
        self.assertEqual(sum(matrix[0]),1)
        self.assertEqual(sum(matrix[q][0] for q in range(1,size))/m,F(1,m))
        self.assertEqual(sum(sum(matrix[q][1:]) for q in range(1,size))/m,1-F(1,m))

    def test_empty_transition_and_mean(self):
        zero,_=kernels(self.a,self.b,1.)
        p=zero[1:,1:]
        m=len(p)
        uniform=np.ones((m,m))/m
        for r in (1,2,7):
            np.testing.assert_allclose(np.linalg.matrix_power(p,r),2.**(-r)*np.eye(m)+(1-2.**(-r))*uniform,atol=1e-14)
        weights=[sum((c&q).bit_count()%2 for c in self.a) for q in range(1,m+1)]
        self.assertEqual(F(sum(weights),m*len(self.a)),F(4,7))

    def test_full_state_region_limits(self):
        size=8;m=7;t=4;theta=3.;p0=4/7;gamma=p0*theta
        a=math.exp(-gamma);f=-math.expm1(-gamma)/gamma
        gap=2*(gamma-1+a)/gamma**2
        edge=2*(1-(1+gamma)*a)/gamma**2
        r=1/m;s=1-r
        limits=[np.diag([1.,a]),np.array([[0,f],[r*f,s*a]]),
                np.array([[r*gap,s*edge],[s*r*edge,r*edge+s*s*a]])]
        lift=np.zeros((size,2));lift[0,0]=1;lift[1:,1]=1
        project=np.zeros((2,size));project[0,0]=1;project[1,1:]=1/m
        errors=[]
        for length in (1<<12,1<<16,1<<20):
            zero,one=kernels(self.a,self.b,math.exp(-theta/length))
            block=np.zeros((3*size,3*size))
            for i in range(3):block[i*size:(i+1)*size,i*size:(i+1)*size]=zero
            for i in range(2):block[i*size:(i+1)*size,(i+1)*size:(i+2)*size]=one
            count=length//t
            result=np.linalg.matrix_power(block,count)
            error=0.
            for q in range(3):
                actual=result[:size,q*size:(q+1)*size]/math.comb(count,q)
                expected=lift @ limits[q] @ project
                error=max(error,float(np.max(abs(actual-expected))))
            errors.append(error)
        self.assertGreater(errors[0],errors[1])
        self.assertGreater(errors[1],errors[2])
        self.assertLess(errors[2],1e-5)


if __name__=='__main__':unittest.main()
