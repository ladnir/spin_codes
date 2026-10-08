"""Exact local checks for the proposed rank-two shared refresh.

These do not constitute a whole-code certificate or a native cost claim.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import random
import unittest
from hd_blockwise_shared_transfer import multiply,power,rank,mm
from hd_sumrank_maps import inverse,Q


class RankTwoRefreshTests(unittest.TestCase):
    def test_augmented_actual_cauchy_map(self):
        a=tuple(tuple(power(i,j,Q) for j in range(16)) for i in range(16))
        b=tuple(tuple(power(i+16,j,Q) for j in range(16)) for i in range(16))
        h=mm(b,inverse(a),Q);rng=random.Random(72718)
        for w in range(15):
            for _ in range(12):
                rows=rng.sample(range(16),w+2);cols=rng.sample(range(16),w)
                augmented=tuple(tuple(h[i][j] for j in cols)+(1,16+i) for i in rows)
                self.assertEqual(rank(augmented,Q),w+2)

    def test_gf8_all_supports_pointwise_law(self):
        q=8;s=3;nodes=(3,4,5)
        table=tuple(tuple(multiply(a,b,q) for b in range(q)) for a in range(q))
        h=tuple(tuple(power(x^y,q-2,q) for y in range(s)) for x in nodes)
        # D1 times any fixed R is uniform nonzero on its support. D2
        # makes every target with the same support equiprobable, so exact
        # counts of zero masks suffice for the full output pointwise law.
        for support in range(1<<s):
            indices=[i for i in range(s) if support>>i&1];counts=Counter();total=0
            for entries in product(range(1,q),repeat=len(indices)):
                image=[0]*s
                for i in range(s):
                    for j,value in zip(indices,entries):image[i]^=table[h[i][j]][value]
                for t0,t1 in product(range(q),repeat=2):
                    value=[image[i]^t0^table[nodes[i]][t1] for i in range(s)]
                    counts[sum(1<<i for i,v in enumerate(value) if not v)]+=1;total+=1
            self.assertEqual(total,(q-1)**len(indices)*q*q)
            for mask,count in counts.items():
                zeros=mask.bit_count()
                if support:self.assertLessEqual(count*(q-1)**zeros,total)
                elif mask==(1<<s)-1:self.assertEqual(F(count,total),F(1,q*q))
                else:self.assertLessEqual(F(count*(q-1)**zeros,total),1-F(1,q*q))

    def test_gf4_acceptance_for_every_diagonal(self):
        q=4;table=tuple(tuple(multiply(a,b,q) for b in range(q)) for a in range(q))
        for d0,d1 in product(range(1,q),repeat=2):
            u=((d0,table[d0][2]),(d1,table[d1][3]));accepted=0
            for a,b,c,d in product(range(q),repeat=4):
                vu=((table[a][u[0][0]]^table[b][u[1][0]],table[a][u[0][1]]^table[b][u[1][1]]),
                    (table[c][u[0][0]]^table[d][u[1][0]],table[c][u[0][1]]^table[d][u[1][1]]))
                determinant=table[1^vu[0][0]][1^vu[1][1]]^table[vu[0][1]][vu[1][0]]
                accepted+=determinant!=0
            self.assertEqual(F(accepted,q**4),(1-F(1,q))*(1-F(1,q*q)))


if __name__=='__main__':unittest.main()
