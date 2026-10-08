"""Exact checks of the block-emission screen's conditional moments."""
from fractions import Fraction as F
from itertools import product
import unittest
import numpy as np
from hd_blockwise_shared_transfer import dot
from hd_bucket_feedback import profiles
from hd_block_emission_screen import exact_moments,coefficients


class BlockEmissionTests(unittest.TestCase):
    def test_active_class_jensen_pair_and_alpha_one(self):
        from hd_block_emission_classes import local_matrices,coefficients as classes
        z=.75;degrees,counts,values=local_matrices(z,2,4,1)
        expected=np.zeros((9,3,3))
        for j,n,a in zip(degrees,counts,values):
            for k,m,b in zip(degrees,counts,values):expected[j+k]+=n*m*np.sqrt(a@b)
        majorant=np.exp2(classes(z,.5,2,4,1,True))
        self.assertTrue(np.all(majorant>=expected*(1-1e-12)))
        np.testing.assert_allclose(np.exp2(classes(z,1,2,4,1,True)),np.exp2(classes(z,1,2,4,1)),rtol=1e-12,atol=1e-20)

    def test_active_state_class_moments(self):
        from hd_block_emission_classes import local_matrices
        q=4;z=F(3,4);a=((1,0),(1,0),(0,1),(0,1));c=((1,1,0,0),(0,0,1,1))
        degrees,counts,actual=local_matrices(float(z),2,4,1,float(q))
        for index,(occ,count) in enumerate(profiles(2,4,1)):
            support=list(range(occ[0]))+list(range(2,2+occ[1]));j=len(support)
            k=t=F(0);vs=[F(0)]*3;quiet={}
            for word in product(range(q),repeat=j):
                x=[0]*4
                for i,value in zip(support,word):x[i]=value
                zero=not any(dot(row,x,q) for row in c)
                mark=z**sum(value!=0 for value in x);t+=mark
                if zero:k+=mark
                for state in product(range(q),repeat=2):
                    u=sum(value!=0 for value in state)
                    if not u:continue
                    mark=z**sum(x[i]^dot(row,state,q)!=0 for i,row in enumerate(a));vs[u]+=mark
                    if zero:quiet[state]=quiet.get(state,F(0))+mark
            scale=q**j;k/=scale;t/=scale;sizes=[1,2*(q-1),(q-1)**2];d=q*q-1
            expected=np.zeros((3,3));expected[0,0]=float(k)
            for v in (1,2):expected[0,v]=float((t-k)*sizes[v]/d)
            for u in (1,2):
                v=vs[u]/(scale*sizes[u]);rho=sum(min(1,l) for l in occ)
                floor=sum(max(0,2-l) for l in occ[-u:]);gamma=min(k,z**floor/q**rho,sizes[u]*v)
                worst=max((mass for state,mass in quiet.items() if sum(x!=0 for x in state)==u),default=F(0))/scale
                self.assertLessEqual(worst,gamma)
                expected[u,0]=float(v/d) if j else 0
                for dest in (1,2):expected[u,dest]=float((v if j else F(q-1,q)*v)*sizes[dest]/d)
                expected[u,u]+=float(gamma/q)
            np.testing.assert_allclose(actual[index],expected,rtol=1e-13,atol=1e-15)
            self.assertEqual(degrees[index],j);self.assertEqual(counts[index],count)

    def test_native_block_emission_digests(self):
        from hd_bucket_feedback import maps
        for h,expected in ((4,5069833586872545664),(8,9075152494480514616)):
            _,c=maps(local_rank=h);a=tuple(zip(*c));f=tuple(tuple(int(i==j) for j in range(64)) for i in range(64));value=0
            for matrix in (f,a,c):
                for row in matrix:
                    for x in row:value=((value^x)*0x100000001b3)&((1<<64)-1)
            self.assertEqual(value,expected)

    def test_gf4_every_input_support_and_nonzero_state(self):
        q=4;z=F(3,4);a=((1,0),(1,0),(0,1),(0,1))
        for local in (False,True):
            c=((1,1,0,0),(0,0,1,1)) if local else ((1,1,1,1),(0,1,2,3))
            for mask in range(16):
                u=[i for i in range(4) if mask>>i&1];occ=(sum(i<2 for i in u),sum(i>=2 for i in u))
                k=t=v=F(0);quiet={}
                for word in product(range(q),repeat=len(u)):
                    x=[0]*4
                    for i,value in zip(u,word):x[i]=value
                    zero=not any(dot(row,x,q) for row in c)
                    mark=z**sum(value!=0 for value in x);t+=mark
                    if zero:k+=mark
                    for state in product(range(q),repeat=2):
                        if not any(state):continue
                        mark=z**sum(x[i]^dot(row,state,q)!=0 for i,row in enumerate(a));v+=mark
                        if zero:quiet[state]=quiet.get(state,F(0))+mark
                scale=q**len(u);k/=scale;t/=scale;v/=scale*(q*q-1)
                expected=exact_moments(occ,z,2,4,1,q,local)
                self.assertEqual((k,t,v),expected[:3])
                self.assertLessEqual(max(quiet.values(),default=F(0))/scale,expected[3])

    def test_pair_screen_keeps_distinct_emission_profiles(self):
        z=F(3,4);q=65536;state=2;width=4
        for local in (False,True):
            exact=np.zeros((9,2,2));items=profiles(state,width,1)
            def matrix(row):
                k,t,v,gamma=exact_moments(row,z,state,width,1,q,local)
                j=sum(row)
                return np.array([[float(k),float(t-k)],[float(v/(q*q-1)) if j else 0,float((v if j else F(q-1,q)*v)+gamma/q)]])
            for left,n in items:
                for right,m in items:exact[sum(left)+sum(right)]+=n*m*np.sqrt(matrix(left)@matrix(right))
            actual=np.exp2(coefficients(float(z),.5,state,width,1,local))
            self.assertTrue(np.all(actual>=exact*(1-2e-12)))
            self.assertTrue(np.all(actual<=exact*(1+1/2000)))


if __name__=='__main__':unittest.main()
