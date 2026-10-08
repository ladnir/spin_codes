"""Independent rank, enumeration, and grouping tests for block feedback."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from math import comb
import random
import unittest
import numpy as np
from hd_blockwise_shared_transfer import rank,dot
from hd_bucket_feedback import maps,moments,matrix,profiles,grouped_packets,digest


class BucketFeedbackTests(unittest.TestCase):
    def test_native_map_digests(self):
        self.assertEqual(digest(4),16990099864024496792)
        self.assertEqual(digest(8),5242589667199767848)
        self.assertEqual(digest(4,True),12097737318272255584)
        self.assertEqual(digest(8,True),1636729993102572008)

    def test_actual_mask_ranks(self):
        rng=random.Random(71936)
        for h,global_nodes in product((4,8),(False,True)):
            _,c=maps(local_rank=h,global_nodes=global_nodes);size=4*h
            for _ in range(90):
                u=sorted(rng.sample(range(64),rng.randrange(65)))
                z=sorted(rng.sample(range(64),rng.randrange(65)))
                overlap=set(u)&set(z)
                f=tuple(tuple(int(i==j) for j in u) for i in z)
                quiet=f+tuple(tuple(row[j] for j in u) for row in c)
                residual=[sum(i//size==b and i not in overlap for i in u) for b in range(16//h)]
                self.assertEqual(rank(quiet,65536),len(overlap)+sum(min(h,l) for l in residual))

    def test_small_field_exact_moments_and_quiet_cap(self):
        q=4;z=F(3,4);a=tuple((1,i) for i in range(4));c=((1,1,0,0),(0,0,1,1))
        for mask in range(16):
            u=tuple(i for i in range(4) if mask>>i&1);occ=(sum(i<2 for i in u),sum(i>=2 for i in u))
            k=t=v=F(0);quiet=Counter()
            for values in product(range(q),repeat=len(u)):
                x=[0]*4
                for i,value in zip(u,values):x[i]=value
                r=tuple(dot(row,x,q) for row in c)
                mark=z**sum(value!=0 for value in x);t+=mark
                if not any(r):k+=mark
                for state in product(range(q),repeat=2):
                    if not any(state):continue
                    weight=z**sum(x[i]^dot(row,state,q)!=0 for i,row in enumerate(a));v+=weight
                    if not any(r):quiet[state]+=weight
            scale=q**len(u);k/=scale;t/=scale;v/=(q*q-1)*scale
            self.assertEqual((k,t,v),moments(occ,z,2,4,1,q))
            rank_c=sum(min(1,j) for j in occ)
            cap=min(k,z**max(0,3-len(u))/q**rank_c,(q*q-1)*v)
            self.assertLessEqual(max(quiet.values(),default=F(0))/scale,cap)

    def test_profile_multiplicities(self):
        for state,width,h in ((2,8,1),(4,8,2),(16,64,4),(16,64,8)):
            rows=profiles(state,width,h)
            for j in range(width+1):self.assertEqual(sum(n for row,n in rows if sum(row)==j),comb(width,j))
        observed=Counter()
        for mask in range(1<<8):
            observed[tuple(sorted((sum(mask>>i&1 for i in range(4)),sum(mask>>i&1 for i in range(4,8)))))]+=1
        self.assertEqual(dict(profiles(2,8,1)),observed)

    def test_grouping_and_pair_screen(self):
        from hd_bucket_feedback_screen import paired_coefficients
        z=F(3,4);state=4;width=8;h=2
        exact=np.zeros((17,2,2))
        items=profiles(state,width,h)
        for row,n in items:
            left=np.array([float(x) for x in matrix(row,z,state,width,h)]).reshape((2,2))
            for other,m in items:
                right=np.array([float(x) for x in matrix(other,z,state,width,h)]).reshape((2,2))
                exact[sum(row)+sum(other)]+=n*m*np.sqrt(left@right)
        screened=np.exp2(paired_coefficients(float(z),.5,state,width,cutoff=width,local_rank=h))
        self.assertTrue(np.all(screened>=exact*(1-2e-13)))
        self.assertTrue(np.all(screened<=exact*(1+1/2000)))
        for cutoff in (0,3,8):
            majorant=np.exp2(paired_coefficients(float(z),.5,state,width,cutoff=cutoff,local_rank=h))
            self.assertTrue(np.all(majorant>=exact*(1-2e-13)))
        grouped=grouped_packets(z,state,width,h)
        for j in range(width+1):self.assertEqual(sum(n for degree,n,_ in grouped if degree==j),comb(width,j))

    def test_alternating_pair_screen(self):
        from hd_bucket_feedback_screen import paired_coefficients
        z=F(3,4);state=4;width=8;exact=np.zeros((17,2,2))
        for row,n in profiles(state,width,2):
            left=np.array([float(x) for x in matrix(row,z,state,width,2)]).reshape((2,2))
            for other,m in profiles(state,width,1):
                right=np.array([float(x) for x in matrix(other,z,state,width,1)]).reshape((2,2))
                exact[sum(row)+sum(other)]+=n*m*np.sqrt(left@right)
        actual=np.exp2(paired_coefficients(float(z),.5,state,width,width,2,1))
        self.assertTrue(np.all(actual>=exact*(1-2e-13)))
        self.assertTrue(np.all(actual<=exact*(1+1/2000)))

    def test_directed_alternating_pair_coefficients(self):
        from flint import arb,ctx
        from hd_randomized_outer_quiet_exact import _arb_fraction
        from hd_sumrank_paired_certificate import pair_coefficients,grouped_packets as selected_groups
        with ctx.workprec(256):
            z=F(3,4);state=8;width=16
            actual,_=pair_coefficients(z,1,state,cutoff=width,feedthrough='bucket4',width=width,right_feedthrough='bucket8')
            expected=[[F(0) for _ in range(33)] for _ in range(4)]
            left=selected_groups(z,state,cutoff=width,feedthrough='bucket4',width=width)
            right=selected_groups(z,state,cutoff=width,feedthrough='bucket8',width=width)
            for j,n,a in left:
                for k,m,b in right:
                    for r in range(2):
                        for s in range(2):expected[2*r+s][j+k]+=n*m*(a[2*r]*b[s]+a[2*r+1]*b[2+s])
            for row,want in zip(actual,expected):
                for value,exact in zip(row,want):self.assertTrue(value.overlaps(_arb_fraction(exact)))


if __name__=='__main__':unittest.main()
