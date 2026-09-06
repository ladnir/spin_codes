"""Exact domination checks for adaptive group choice and root weighting."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np
from flint import arb,ctx
import activation_bridge as q1
import general_occupancy as general
import bridge as base
import adaptive_range as screen
from certify_adaptive_range import adaptive_kernel
from audit_bch_q1_full_arb import rational


def exact_adaptive(region,ps,roots):
    current=region
    for _ in range(len(region)-1):
        current=[tuple(max(r*((1-p)*a[k]+p*b[k]) for p,r in zip(ps,roots)) for k in range(9))
                 for a,b in zip(current,current[1:])]
    return current[0]


class AdaptiveTests(unittest.TestCase):
    def test_every_fixed_group_sequence(self):
        region=general.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},F(3,4),F,4,length=16)
        ps=(F(1,4),F(1,2),F(3,4));roots=(F(3,2),F(5,4),F(2))
        checked=0
        for q in range(1,5):
            upper=exact_adaptive(region[:q+1],ps,roots)
            for sequence in itertools.product(range(3),repeat=q):
                exact=[F(0)]*9
                for bits in itertools.product((0,1),repeat=q):
                    mass=math.prod(roots[g]*(ps[g] if bit else 1-ps[g]) for g,bit in zip(sequence,bits))
                    exact=[a+mass*b for a,b in zip(exact,region[sum(bits)])]
                self.assertTrue(all(a<=b for a,b in zip(exact,upper)));checked+=1
            ctx.prec=256
            rounded=adaptive_kernel([tuple(arb(v.numerator)/v.denominator for v in m) for m in region[:q+1]],
                                    [arb(p.numerator)/p.denominator for p in ps],
                                    [arb(r.numerator)/r.denominator for r in roots])
            self.assertTrue(all(exact<=rational(v.upper()) for exact,v in zip(upper,rounded)))
        self.assertEqual(checked,120)

    def test_counting_cost_distributed_across_regions(self):
        # Four outer coordinates: choose exact roots whose fourth powers bound
        # each group's entire unnormalized counting density.
        positions=4;caps={1:2,2:3,3:1,4:1};bands=((1,2),(3,4));ps=(F(1,3),F(3,4))
        gamma=[max(F(caps[w],math.comb(positions,w))/(p**w*(1-p)**(positions-w)) for w in band)
               for band,p in zip(bands,ps)]
        roots=tuple(F(math.ceil(float(g)**(1/positions))+1) for g in gamma)
        self.assertTrue(all(r**positions>=g for r,g in zip(roots,gamma)))
        region=general.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},F(3,4),F,3,length=16)
        upper=exact_adaptive(region,ps,roots)
        upper_power=upper
        for _ in range(positions-1):upper_power=q1.positive_mul(upper_power,upper)
        total=F(0)
        for sequence in itertools.product(range(2),repeat=3):
            mixed=[F(0)]*9
            for bits in itertools.product((0,1),repeat=3):
                mass=math.prod(ps[g] if bit else 1-ps[g] for g,bit in zip(sequence,bits))
                mixed=[a+mass*b for a,b in zip(mixed,region[sum(bits)])]
            value=tuple(mixed)
            for _ in range(positions-1):value=q1.positive_mul(value,mixed)
            bound=math.prod(gamma[g] for g in sequence)*sum(value[:3])
            self.assertLessEqual(bound,sum(upper_power[:3]));total+=bound
        self.assertLessEqual(total,2**3*sum(upper_power[:3]))

    def test_normalized_float_regions(self):
        epoch=general.epoch_matrices(8,4,{4:14,8:1},{0:1,4:14,8:1},F(3,4),F,8)
        expected=general.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},F(3,4),F,8,length=16)
        actual=screen.normalized_regions(np.array([[float(v) for v in row] for row in epoch]).reshape(-1,3,3),8,t=8,length=16)
        np.testing.assert_allclose(actual,np.array([[float(v) for v in row] for row in expected]).reshape(-1,3,3),rtol=2e-13,atol=1e-15)

    def test_saved_range_sum(self):
        saved=base.read(base.HERE/'generated'/'adaptive_q17_q64_outward.json')
        self.assertEqual([row['occupation'] for row in saved['rows']],list(range(17,65)))
        total=sum((base.decode(row['upper']) for row in saved['rows']),F(0))
        self.assertEqual(total,base.decode(saved['Q17_through_Q64_upper']))
        self.assertLess(total,F(1,1<<142))
        prior=base.read(base.HERE/'generated'/'general_q4_q16_outward.json')
        partial=total+base.decode(prior['Q1_through_Q16_upper'])
        self.assertEqual(partial,base.decode(saved['Q1_through_Q64_upper']))
        self.assertLess(partial,F(1,1<<49))


if __name__=='__main__':unittest.main()
