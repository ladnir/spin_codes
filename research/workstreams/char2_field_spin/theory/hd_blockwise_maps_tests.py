"""Exact small-field validation of the fixed maps and local measure cone."""

from fractions import Fraction as F
from itertools import combinations,product
from random import Random
import unittest

import hd_blockwise_inner as inner
import hd_blockwise_maps as proof


def rank(matrix,q):
    work=[list(row) for row in matrix]
    found=0
    for col in range(len(work[0]) if work else 0):
        pivot=next((i for i in range(found,len(work)) if work[i][col]),None)
        if pivot is None:
            continue
        work[found],work[pivot]=work[pivot],work[found]
        inverse=proof.power(work[found][col],q-2,q)
        work[found]=[proof.mul(inverse,x,q) for x in work[found]]
        for i in range(found+1,len(work)):
            factor=work[i][col]
            work[i]=[x^proof.mul(factor,y,q) for x,y in zip(work[i],work[found])]
        found+=1
    return found


class BlockwiseMapsTests(unittest.TestCase):
    def test_native_geometry_and_ca(self):
        receipt=proof.native_receipt()
        self.assertEqual(receipt['A_sha256'],
                         '4a8399e103c2fb2ca5824883d7da1f6d3e235c4e5873c8467701a6cde071586a')
        a,c=proof.maps()
        self.assertEqual(c,tuple(zip(*a)))
        self.assertEqual(rank(a,256),16)
        self.assertEqual(rank(c,256),16)
        # Interpolation of the constant polynomial, original state coordinates.
        self.assertEqual(proof.matvec(a,(1,)*16),(1,)*64)
        # Four coset translations, each self-inverse in characteristic two.
        identity=tuple(tuple(int(i==j) for j in range(16)) for i in range(16))
        for offset in (0,16,32,48):
            block=a[offset:offset+16]
            self.assertEqual(proof.matmul(block,block),identity)

    def test_small_mds_and_nilpotent_map(self):
        q,t,s=8,4,2
        a,c=proof.maps(t,s,4,q)
        self.assertFalse(any(x for row in proof.matmul(c,a,q) for x in row))
        for size in (1,2):
            for rows in combinations(range(t),size):
                for cols in combinations(range(s),size):
                    self.assertEqual(rank([[a[i][j] for j in cols] for i in rows],q),size)
        emission=[sum(x!=0 for x in proof.matvec(a,z,q))
                  for z in product(range(q),repeat=s) if any(z)]
        kernel=[sum(v!=0 for v in x) for x in product(range(q),repeat=t)
                if any(x) and not any(proof.matvec(c,x,q))]
        self.assertEqual(min(emission),3)
        self.assertEqual(min(kernel),3)

    def test_smaller_packet_candidate_maps(self):
        for width,state in ((32,8),(16,4)):
            a,c=proof.maps(width,state)
            self.assertEqual(rank(a,256),state)
            self.assertEqual(rank(c,256),state)
            self.assertFalse(any(x for row in proof.matmul(c,a) for x in row))
            self.assertEqual(proof.matvec(a,(1,)*state),(1,)*width)

    def test_uniform_emission_bound_all_small_inputs(self):
        q,t,s=8,4,2
        a,c=proof.maps(t,s,4,q)
        states=tuple(z for z in product(range(q),repeat=s) if any(z))
        emissions=[proof.matvec(a,z,q) for z in states]
        theta=F(2,3)
        for x in product(range(q),repeat=t):
            weight=sum(v!=0 for v in x)
            flag=any(proof.matvec(c,x,q))
            floor=proof.emission_minimum(weight,t,s,flag)
            marks=[theta**sum(u!=v for u,v in zip(x,az)) for az in emissions]
            self.assertTrue(all(mark<=theta**floor for mark in marks))
            self.assertLessEqual(sum(marks)/len(states),
                                 proof.emission_upper(weight,theta,t,s,q))

    def test_local_cone_pointwise_and_mass_domination(self):
        checked=0
        for q,t,s,base in ((4,2,1,2),(8,4,2,4)):
            a,c=proof.maps(t,s,base,q)
            states=tuple(product(range(q),repeat=s))
            live=states[1:]
            index={z:i for i,z in enumerate(states)}
            emissions={z:proof.matvec(a,z,q) for z in live}
            kernel={z:tuple(inner.theoretical_image(z,target,q) for target in states) for z in live}
            rng=Random(8416)
            inputs=(tuple(product(range(q),repeat=t)) if q==4 else
                    ((0,)*t,)+tuple(tuple(rng.randrange(q) for _ in range(t)) for _ in range(14))
                    +tuple(emissions[z] for z in live[:8]))
            for theta in (F(1,3),F(1)):
                for x in inputs:
                    weight=sum(v!=0 for v in x)
                    syndrome=proof.matvec(c,x,q)
                    flag=any(syndrome)
                    actual_uniform=[F(0)]*len(states)
                    matrix=proof.local_transfer(weight,theta,t,s,q,syndrome_nonzero=flag)
                    for z in live:
                        mark=theta**sum(u!=v for u,v in zip(x,emissions[z]))
                        actual=[mark*kernel[z][index[tuple(u^v for u,v in zip(target,syndrome))]]
                                for target in states]
                        for i,value in enumerate(actual):
                            actual_uniform[i]+=value/len(live)
                        self.assertLessEqual(actual[0],matrix[1][0])
                        cap=matrix[1][2]/len(live)
                        self.assertLessEqual(sum(max(F(0),v-cap) for v in actual[1:]),matrix[1][1])
                        checked+=1
                    for mode in ('density','mass'):
                        m=proof.local_transfer(weight,theta,t,s,q,syndrome_nonzero=flag,uniform_lazy=mode)
                        self.assertLessEqual(actual_uniform[0],m[2][0])
                        cap=m[2][2]/len(live)
                        self.assertLessEqual(sum(max(F(0),v-cap) for v in actual_uniform[1:]),m[2][1])
                    # Unknown syndrome is a simultaneous entrywise upper bound.
                    unknown=proof.local_transfer(weight,theta,t,s,q)
                    self.assertTrue(all(unknown[i][j]>=matrix[i][j] for i in range(3) for j in range(3)))
        self.assertEqual(checked,2994)

    def test_multi_packet_weighted_envelope(self):
        q,t,s=8,4,2
        a,c=proof.maps(t,s,4,q)
        states=tuple(product(range(q),repeat=s))
        index={z:i for i,z in enumerate(states)}
        kernel={z:tuple(inner.theoretical_image(z,target,q) for target in states) for z in states}
        rng=Random(6491)
        theta=F(1,2)
        actual=[F(int(not any(z))) for z in states]
        envelope=(F(1),F(0),F(0))
        for step in range(10):
            x=((0,)*t if step in (2,5,6) else tuple(rng.randrange(q) for _ in range(t)))
            syndrome=proof.matvec(c,x,q)
            following=[F(0)]*len(states)
            for i,z in enumerate(states):
                mark=theta**sum(u!=v for u,v in zip(x,proof.matvec(a,z,q)))
                for j,target in enumerate(states):
                    shifted=tuple(u^v for u,v in zip(target,syndrome))
                    following[j]+=actual[i]*mark*kernel[z][index[shifted]]
            actual=following
            matrix=proof.local_transfer(sum(v!=0 for v in x),theta,t,s,q,
                                        syndrome_nonzero=any(syndrome),
                                        uniform_lazy='mass' if step%2 else 'density')
            envelope=proof.apply_local(envelope,matrix)
            self.assertLessEqual(actual[0],envelope[0])
            cap=envelope[2]/(q**s-1)
            self.assertLessEqual(sum(max(F(0),v-cap) for v in actual[1:]),envelope[1])
            self.assertLessEqual(sum(actual),sum(envelope))


if __name__=='__main__':
    unittest.main()
