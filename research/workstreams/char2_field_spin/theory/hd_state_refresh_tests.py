"""Exact local validation of the proposed state-only refresh envelope."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import unittest
from hd_blockwise_shared_transfer import multiply,dot,mm
from hd_block_emission_screen import exact_moments


class StateRefreshTests(unittest.TestCase):
    def test_all_small_packet_input_supports_and_incoming_states(self):
        q=4;s=2;d=15;z=F(3,4);h=((3,2),(2,3));zero=(0,0)
        a=((1,0),(1,0),(0,1),(0,1));c=mm(h,tuple(zip(*a)),q)
        vectors=tuple(product(range(q),repeat=s));diags=tuple(product(range(1,q),repeat=s))
        b=F(d,(q-1)**s);den=len(diags)**2
        @lru_cache(None)
        def refresh(incoming):
            counts=Counter()
            for d1 in diags:
                mid=tuple(dot(row,tuple(multiply(x,y,q) for x,y in zip(d1,incoming)),q) for row in h)
                for d2 in diags:counts[tuple(multiply(x,y,q) for x,y in zip(d2,mid))]+=1
            return counts
        for mask in range(16):
            u=[i for i in range(4) if mask>>i&1];j=len(u)
            occ=(sum(i<2 for i in u),sum(i>=2 for i in u));rho=sum(min(1,l) for l in occ)
            k,t,v,_=exact_moments(occ,z,s,4,1,q,True)
            live=t*z**max(0,2-max(occ));density=F(d,q**rho)*(z+(1-z)/q)**(j-rho)
            flat=Counter()
            for incoming in vectors:
                output=Counter()
                for word in product(range(q),repeat=j):
                    x=[0]*4
                    for i,value in zip(u,word):x[i]=value
                    cx=tuple(dot(row,x,q) for row in c)
                    weight=z**sum(x[i]^dot(row,incoming,q)!=0 for i,row in enumerate(a))
                    for state,count in refresh(incoming).items():
                        following=tuple(x^y for x,y in zip(state,cx))
                        output[following]+=weight*count/F(den*q**j)
                if incoming==zero:
                    self.assertEqual(output[zero],k)
                    self.assertEqual(sum(value for w,value in output.items() if w!=zero),t-k)
                    self.assertLessEqual(d*max((value for w,value in output.items() if w!=zero),default=0),density)
                else:
                    self.assertLessEqual(output[zero],b*live/d if j else 0)
                    self.assertLessEqual(d*max((value for w,value in output.items() if w!=zero),default=0),b*live)
                    for w,value in output.items():flat[w]+=value/d
            self.assertLessEqual(flat[zero],b*v/d if j else 0)
            self.assertLessEqual(d*max((value for w,value in flat.items() if w!=zero),default=0),b*v)


if __name__=='__main__':unittest.main()
