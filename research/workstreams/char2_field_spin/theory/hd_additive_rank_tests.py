"""Small exact/support checks for the unconditioned rank envelope."""
from itertools import product
from fractions import Fraction
import math
import unittest

import numpy as np

from hd_additive_rank_screen import (double_prefix_rank_transfers, product_first_bound,
                                     rank_transfers, row_log_transfer)
from hd_rational_post_tests import mul


def plus(a,b):
    return a[0]^b[0],a[1]^b[1]


def scale(c,a):
    return mul(c,a[0]),mul(c,a[1])


def det(a,b):
    return mul(a[0],b[1]) ^ mul(a[1],b[0])


def alive_rank(zeros,live):
    nonzero = [z for z in zeros if z != (0,0)]
    if any(v == (0,0) for v in live):
        return None
    if not nonzero:
        return 0
    pivot = nonzero[0]
    if any(det(pivot,z) for z in nonzero[1:]):
        return None
    if any(not det(pivot,v) for v in live):
        return None
    return 1


class AdditiveRankTests(unittest.TestCase):
    def test_sharp_pivot_keeps_required_live_scalar(self):
        forms = tuple(product(range(4), repeat=2))
        for tau in (1, 2, 3, 8):
            for live in forms[1:]:
                for fixed in forms:
                    mass = Fraction(0)
                    for coefficient in range(1, 4):
                        equation = plus(fixed, scale(coefficient, live))
                        rank = alive_rank([equation], [live])
                        if rank is not None:
                            mass += Fraction(1, 3 * tau**rank)
                    self.assertLessEqual(mass, max(Fraction(1, 3), Fraction(1, tau)))

    def test_local_envelope_on_dependent_gf4_forms(self):
        forms = ((0,0),(1,0),(0,1),(1,1),(2,1))
        q,tau,y = 4,3.,.4
        g = max(1/(q-1),1/tau)
        matrices = rank_transfers(y,g)
        for x in ((0,0),(1,0)):
            active = int(x != (0,0))
            for s1,s2 in product(forms,repeat=2):
                old = 2*int(s1 != (0,0))+int(s2 != (0,0))
                history = [v for v in (x,s1,s2,(1,2)) if v != (0,0)]
                mass = np.zeros(4)
                for c,d1,d2,f in product(range(1,4),repeat=4):
                    if not(c ^ d1 ^ mul(f,d2)):
                        continue
                    output = plus(plus(s1,scale(f,s2)),scale(c,x))
                    new1,new2 = plus(s1,scale(d1,x)),plus(s2,scale(d2,x))
                    for new,live_output in product(range(4),range(2)):
                        zeros,live = [],list(history)
                        for value,flag in ((new1,new&2),(new2,new&1),(output,live_output)):
                            (live if flag else zeros).append(value)
                        rank = alive_rank(zeros,live)
                        if rank is not None:
                            mass[new] += tau**(-rank)*y**live_output/81
                # Division by81 means this is P[valid AND path] in the
                # independent experiment, before the separate1/Z factor.
                self.assertTrue(np.all(mass <= matrices[active][old]+1e-12))

    def test_one_prefix_exact_support_row_count(self):
        b,k,tau,g = 4,2,3.,.4
        actual = np.exp(row_log_transfer(b,k,math.log(tau),g)+k*math.log(tau)/2)
        for old_weight in range(b+1):
            old = (1,)*old_weight+(0,)*(b-old_weight)
            expected = np.zeros(b+1)
            for source in product(range(2),repeat=b):
                w = sum(source)
                if 0 < w <= b-k:
                    continue
                source_weight = 1. if not w else tau**(w-b+k)
                for new in product(range(2),repeat=b):
                    value = source_weight
                    for before,active,after in zip(old,source,new):
                        if (not active and before != after) or (active and not before and not after):
                            value = 0
                            break
                        if active and before and not after:
                            value *= g
                    expected[sum(new)] += value
            np.testing.assert_allclose(actual[old_weight],expected,rtol=1e-12,atol=1e-12)

    def test_double_prefix_product_envelope_all_small_histories(self):
        b,k,rows,tau,g = 4,2,1,5.,.3
        for mark in ((-.8,-1.,.2,.4),(-2.,.5,-.2,-1.)):
            logs = np.array([0.,*mark[1:]])
            v = np.exp(logs)
            t0,t1 = double_prefix_rank_transfers(math.exp(mark[0]),g)
            logupper = product_first_bound(mark,b,k,math.log(tau),rows,g)
            logupper -= b*(math.log(4)+float(logs.max()-logs.min()))
            upper = math.exp(logupper)
            ratio = ((t0 @ v)/v,(t1 @ v)/v)
            for old in product(range(4),repeat=b):
                actual = 0.
                for active in product(range(2),repeat=b):
                    w = sum(active)
                    if 0 < w <= b-k:
                        continue
                    coefficient = 1. if not w else tau**(w-b+k)
                    actual += coefficient*math.prod(ratio[a][o] for a,o in zip(active,old))
                self.assertLessEqual(actual,upper+1e-10)


if __name__ == '__main__':
    unittest.main()
