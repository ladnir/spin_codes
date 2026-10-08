"""Exact local validation of the common-map recurrence and return bound."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,permutations,product
from math import comb,prod
import unittest
import numpy as np

from hd_blockwise_shared_transfer import multiply,dot,mm
from hd_block_emission_screen import exact_moments


def return_cap(occ,z,h,length,q):
    quiet=lambda l:sum(F(comb(l,c),q**(c+min(h,l-c)))*z**(l-c)*(1-z)**c for c in range(l+1))
    k=[quiet(l) for l in occ]
    r=[min((z+(1-z)/q)**l-kb,(1-F(1,q**min(h,l)))*z**max(0,length-h+1-l)) for l,kb in zip(occ,k)]
    return prod(a+b for a,b in zip(k,r))-prod(k)


class CommonRefreshTests(unittest.TestCase):
    def test_scattered_map_digest_and_check_ranks(self):
        from hd_scattered_feedback import maps,check_ranks,small_return_counts
        from hd_common_refresh_enumerator import extend
        from hd_systematic_mds import novel
        import random
        a,c=maps();f=tuple(tuple(int(i==j) for j in range(64)) for i in range(64));value=0
        for matrix in (f,a,c):
            for row in matrix:
                for x in row:value=((value^x)*0x100000001b3)&((1<<64)-1)
        self.assertEqual(value,1225673016892679168)
        self.assertFalse(any(x for row in mm(c,a,65536) for x in row))
        rng=random.Random(75191)
        for mask in range(16):
            ranks=check_ranks(mask)
            for _ in range(100):
                support=rng.sample(range(16),rng.randrange(17));basis=()
                for j in support:basis=extend(basis,tuple(novel(j^16,i) for i in range(4) if mask>>i&1))
                self.assertLessEqual(ranks[len(support)],len(basis))
        counts,denominator=small_return_counts()[(0,0,0,1)]
        self.assertEqual(min(counts),56)
        self.assertEqual(sum(counts.values()),64*65535)

    def test_pinned_small_support_enumerators(self):
        from hd_common_refresh_enumerator import average_counts,pinned_counts
        for l in (1,2):self.assertEqual(average_counts(l),pinned_counts(l))
        for l in (3,4):
            values,den=pinned_counts(l)
            self.assertEqual(sum(values),den)
            self.assertEqual(values[0],comb(16,l))

    def test_zero_set_lattice_against_exhaustive_field_words(self):
        from hd_common_refresh_enumerator import enumerator
        import random
        rng=random.Random(81931)
        # Binary coefficients have the same ranks in GF4 and GF256.
        for dimension in (1,2,3,4):
            rows=tuple(tuple(int(i==j) for j in range(dimension)) for i in range(dimension))
            rows+=tuple(tuple(rng.randrange(2) for _ in range(dimension)) for _ in range(7))
            direct=[0]*(len(rows)+1)
            for word in product(range(4),repeat=dimension):
                direct[sum(dot(row,word,4)!=0 for row in rows)]+=1
            self.assertEqual(tuple(direct),enumerator(rows,dimension,4))
        # Also exercise genuinely nonbinary coefficients over the AES field.
        rows=((1,0),(0,1),(1,1),(1,2),(0,0),(2,3));direct=[0]*7
        table=tuple(tuple(multiply(i,j,65536) for j in range(256)) for i in range(4))
        for x,y in product(range(256),repeat=2):
            direct[sum(table[a][x]^table[b][y]!=0 for a,b in rows)]+=1
        self.assertEqual(tuple(direct),enumerator(rows,2,256))

    def test_return_rank_caps_on_concrete_supports(self):
        from hd_bucket_feedback import maps
        from hd_common_refresh_enumerator import extend
        from hd_common_refresh_screen import principal_defects
        import random
        _,full=maps(local_rank=4);c=tuple(row[:16] for row in full[:4]);a=tuple(zip(*c));ac=mm(a,c,65536)
        p=tuple(tuple(x^int(i==j) for j,x in enumerate(row)) for i,row in enumerate(ac))
        def rank(rows):
            basis=()
            for row in rows:basis=extend(basis,row)
            return len(basis)
        rng=random.Random(59137);nu=principal_defects(4);q=65536
        for _ in range(1500):
            u=rng.sample(range(16),rng.randrange(17));zeros=rng.sample(range(16),rng.randrange(17))
            l=len(u);cc=len(set(u)&set(zeros));e=len(zeros)-cc
            pp=tuple(tuple(p[i][j] for j in u) for i in zeros)
            rr=tuple(tuple(row[j] for j in u) for row in c)
            actual=F(1,q**rank(pp))-F(1,q**rank(pp+rr))
            if e>=4:self.assertEqual(actual,0);continue
            exponent=max(0,cc-min(nu[cc],max(0,4-l+cc)))
            count_cap=F(q**min(min(4,l),4-e)-1,q**(l-max(0,l-cc-4)))
            bound=min(1-F(1,q**min(4,l)),F(1,q**exponent),count_cap)
            self.assertLessEqual(actual,bound)

    def test_all_actual_principal_nullities_by_independent_minors(self):
        from hd_bucket_feedback import maps
        from hd_common_refresh_principal import nullities
        _,full=maps(local_rank=4);c=np.array([row[:16] for row in full[:4]],dtype=np.uint16)
        table=np.array([[multiply(i,j,65536) for j in range(256)] for i in range(256)],dtype=np.uint16)
        # Build subsets by an explicit bit-matrix reduction, independently
        # of the scan's doubling construction and Gaussian elimination.
        masks=np.arange(1<<16);gram=np.broadcast_to(np.eye(4,dtype=np.uint16),(1<<16,4,4)).copy()
        for i in range(16):gram^=((masks>>i&1).astype(np.uint16)[:,None,None]*table[c[:,i,None],c[None,:,i]])
        ranks=np.zeros(1<<16,dtype=int)
        for size in (1,2,3,4):
            for rows in combinations(range(4),size):
                for columns in combinations(range(4),size):
                    det=np.zeros(1<<16,dtype=np.uint16)
                    for perm in permutations(columns):
                        term=np.ones(1<<16,dtype=np.uint16)
                        for i,j in zip(rows,perm):term=table[term,gram[:,i,j]]
                        det^=term
                    ranks[det!=0]=size
        np.testing.assert_array_equal(4-ranks,nullities())
        self.assertEqual(int(np.count_nonzero(4-ranks)),704)
        self.assertEqual(int((4-ranks).max()),2)

    def test_actual_maps_self_orthogonal(self):
        from hd_bucket_feedback import maps
        for h in (4,8):
            _,c=maps(local_rank=h);a=tuple(zip(*c))
            self.assertFalse(any(x for row in mm(c,a,65536) for x in row))

    def test_all_small_packet_operators(self):
        self.check_small_packet(False)

    def test_cross_block_small_packet_operators(self):
        self.check_small_packet(True)

    def check_small_packet(self,cross):
        q=4;s=2;d=q**s-1;z=F(3,4);h=((3,2),(2,3));zero=(0,0)
        a=((1,0),(1,0),(0,1),(0,1));local_c=tuple(zip(*a))
        c=mm(h,local_c,q) if cross else local_c
        vectors=tuple(product(range(q),repeat=s));diags=tuple(product(range(1,q),repeat=s))
        cap=F(d,(q-1)**s);den=(q-1)**(2*s)
        @lru_cache(None)
        def transitions(w):
            counts=Counter()
            for d1 in diags:
                mid=tuple(dot(row,tuple(multiply(x,y,q) for x,y in zip(d1,w)),q) for row in h)
                for d2 in diags:counts[tuple(multiply(x,y,q) for x,y in zip(d2,mid))]+=1
            if w==zero:self.assertEqual(counts,{zero:den})
            else:
                self.assertEqual(counts[zero],0)
                self.assertLessEqual(F(max(counts.values()),den),cap/d)
            return tuple(counts.items())
        for mask in range(16):
            u=[i for i in range(4) if mask>>i&1];j=len(u);occ=(sum(i<2 for i in u),sum(i>=2 for i in u))
            outputs=(Counter(),Counter());ret=F(0)
            for word in product(range(q),repeat=j):
                x=[0]*4
                for i,v in zip(u,word):x[i]=v
                r=tuple(dot(row,x,q) for row in c)
                for incoming in vectors:
                    weight=z**sum(x[i]^dot(row,incoming,q)!=0 for i,row in enumerate(a))
                    if incoming==r and incoming!=zero:ret+=weight
                    for following,count in transitions(tuple(i^v for i,v in zip(incoming,r))):
                        outputs[int(incoming!=zero)][following]+=weight*count
            scale=q**j;ret/=scale
            k,t,v,_=exact_moments(occ,z,s,4,1,q,True)
            if cross:
                p0=[F(1,q**min(1,l)) for l in occ]
                probabilities=(p0[0]*p0[1],p0[0]*(1-p0[1])+p0[1]*(1-p0[0]),(1-p0[0])*(1-p0[1]))
                jcap=min(t-k,sum(probabilities[w]*z**max(0,2*(3-w)-j) for w in (1,2)),d*v)
            else:jcap=min(return_cap(occ,z,1,2,q),d*v)
            self.assertLessEqual(ret,jcap)
            if not cross:
                from hd_common_refresh_screen import local_return_cap
                # This small construction has one state coordinate per
                # local block. Enumerate each block return independently.
                for local_mask in range(4):
                    positions=[i for i in range(2) if local_mask>>i&1];local_ret=F(0)
                    for entries in product(range(q),repeat=len(positions)):
                        xx=[0,0]
                        for i,value in zip(positions,entries):xx[i]=value
                        rr=xx[0]^xx[1]
                        if rr:local_ret+=z**sum(v^rr!=0 for v in xx)
                    local_ret/=q**len(positions)
                    self.assertLessEqual(float(local_ret),local_return_cap(len(positions),float(z),1,2,q)+1e-14)
            # Test the sharper exact-J operator, as well as the computable cap.
            bounds=((k,t-k),(ret/d,v-ret/d))
            for values,normalizer,bound in zip(outputs,(den*scale,den*scale*d),bounds):
                self.assertEqual(values[zero]/normalizer,bound[0])
                self.assertLessEqual(d*max((x for w,x in values.items() if w!=zero),default=0)/normalizer,cap*bound[1])


if __name__=='__main__':unittest.main()
