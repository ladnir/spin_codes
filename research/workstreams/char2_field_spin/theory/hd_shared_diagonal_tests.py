"""Exact small-field checks of the correlated refresh and its transfer bound."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
import random
import unittest

from hd_blockwise_shared_transfer import multiply,power,dot,mm,rank
from hd_shared_diagonal import transfer
from hd_sumrank_maps import inverse,Q


class SharedDiagonalTests(unittest.TestCase):
    def test_augmented_cauchy_actual_map(self):
        a=tuple(tuple(power(i,j,Q) for j in range(16)) for i in range(16))
        b=tuple(tuple(power(i+16,j,Q) for j in range(16)) for i in range(16))
        h=mm(b,inverse(a),Q)
        constant=multiply(h[0][0],16,Q)
        self.assertNotEqual(constant,0)
        for i,row in enumerate(h):
            for j,x in enumerate(row):self.assertEqual(multiply(x,16^i^j,Q),constant)
            self.assertEqual(dot(row,(1,)*16,Q),1)
        rng=random.Random(525831)
        # Sample checks supplement the rational-function proof; they do not
        # establish the universal minor property by themselves.
        for w in range(16):
            for _ in range(8):
                rows=rng.sample(range(16),w+1);cols=rng.sample(range(16),w)
                self.assertEqual(rank(tuple(tuple(h[i][j] for j in cols)+(1,) for i in rows),Q),w+1)

    def test_all_small_field_joint_distributions(self):
        q=4;h=((3,2),(2,3));zero=(0,0)
        table=tuple(tuple(multiply(a,b,q) for b in range(q)) for a in range(q))
        vectors=tuple(product(range(q),repeat=2));diags=tuple(product(range(1,q),repeat=2))
        dot4=lambda a,b:table[a[0]][b[0]]^table[a[1]][b[1]]
        for cols in combinations(range(3),2):
            self.assertEqual(rank(tuple(tuple((h[i]+(1,))[j] for j in cols) for i in range(2)),q),2)
        for incoming in vectors:
            for r in vectors:
                counts=Counter();unconditioned=Counter()
                for d1 in diags:
                    dr=tuple(table[a][x] for a,x in zip(d1,r))
                    hr=tuple(dot4(row,dr) for row in h)
                    for d2 in diags:
                        for v in vectors:
                            tau=dot4(v,incoming)
                            out=tuple(x^table[a][y^tau] for x,a,y in zip(incoming,d2,hr))
                            unconditioned[out]+=1
                            if dot4(v,d2)!=1:counts[out]+=1
                total=(q-1)**4*q*q;accepted=total*(q-1)//q
                self.assertEqual(sum(unconditioned.values()),total)
                self.assertEqual(sum(counts.values()),accepted)
                if incoming!=zero and r!=zero:
                    self.assertLessEqual(F(max(unconditioned.values()),total),F(1,(q-1)**2))
                    self.assertLessEqual(F(max(counts.values()),accepted),F(q,(q-1)**3))
                if r==zero and incoming!=zero:
                    self.assertEqual(counts[zero],0)
                    self.assertEqual(F(unconditioned[incoming],total),F(1,q))
                    for target,count in unconditioned.items():
                        if target!=incoming:self.assertLessEqual(F(count,total),F(1,q*(q-1)))
                if incoming==zero:
                    if r==zero:self.assertEqual(counts[zero],accepted)
                    else:
                        self.assertEqual(counts[zero],0)
                        self.assertLessEqual(F(max(counts.values()),accepted),F(1,(q-1)**2))

    def test_directed_transfer(self):
        from flint import arb,ctx
        with ctx.workprec(384):
            result=transfer(precision=384)
            self.assertTrue(result['whole_code_certificate'])
            self.assertTrue(arb(result['additional_loss_bits'])<arb('.722'))
            self.assertTrue(arb(result['union_log2_upper'])<arb('-49.338'))
            self.assertFalse(result['underlying_certificate_replayed'])

    def test_all_small_packet_weighted_operators(self):
        self.check_small_packet(False)

    def test_block_feedback_weighted_operators(self):
        self.check_small_packet(True)

    def check_small_packet(self,block_feedback):
        # Independent direct packet enumeration checks both incoming measure
        # components, including the quiet-state atom and both zero entries.
        from functools import lru_cache
        q=4;s=2;dim=15;h=((3,2),(2,3));zero=(0,0);z=F(3,4)
        vectors=tuple(product(range(q),repeat=s));diags=tuple(product(range(1,q),repeat=s))
        a=tuple((1,i) for i in range(4))
        c=((1,1,0,0),(0,0,1,1)) if block_feedback else tuple(zip(*a))
        cap=F(q,q-1)*F(dim,(q-1)**s)
        table=tuple(tuple(multiply(i,j,q) for j in range(q)) for i in range(q))
        dot4=lambda x,y:table[x[0]][y[0]]^table[x[1]][y[1]]
        @lru_cache(None)
        def transitions(incoming,r):
            counts=Counter()
            for d1 in diags:
                dr=tuple(table[v][x] for v,x in zip(d1,r))
                hr=tuple(dot4(row,dr) for row in h)
                for d2 in diags:
                    for v in vectors:
                        if dot4(v,d2)==1:continue
                        tau=dot4(v,incoming)
                        out=tuple(x^table[t][y^tau] for x,t,y in zip(incoming,d2,hr))
                        counts[out]+=1
            return tuple(counts.items())
        accepted=(q-1)**4*q*(q-1)
        for mask in range(16):
            u=tuple(i for i in range(4) if mask>>i&1);j=len(u)
            zero_out=Counter();flat_out=Counter();quiet=Counter();k=F(0);t=F(0);vtotal=F(0)
            for coefficients in product(range(q),repeat=j):
                x=[0]*4
                for i,value in zip(u,coefficients):x[i]=value
                r=tuple(dot(row,x,q) for row in c)
                mark=z**sum(value!=0 for value in x);t+=mark
                if r==zero:k+=mark
                for incoming in vectors:
                    output=tuple(x[i]^dot4(row,incoming) for i,row in enumerate(a))
                    weight=z**sum(value!=0 for value in output)
                    target=zero_out if incoming==zero else flat_out
                    for following,count in transitions(incoming,r):target[following]+=weight*count
                    if incoming!=zero:
                        vtotal+=weight
                        if r==zero:quiet[incoming]+=weight
            scale=q**j;k/=scale;t/=scale;v=vtotal/(dim*scale)
            gamma=max(quiet.values(),default=F(0))/scale
            bounds=((k,t-k),(v/dim if j else 0,(v if j else F(q-1,q)*v)+gamma/q))
            for values,denominator,bound in ((zero_out,accepted*scale,bounds[0]),(flat_out,accepted*scale*dim,bounds[1])):
                self.assertLessEqual(values[zero]/denominator,cap*bound[0])
                self.assertLessEqual(dim*max((value for key,value in values.items() if key!=zero),default=0)/denominator,cap*bound[1])


if __name__=='__main__':unittest.main()
