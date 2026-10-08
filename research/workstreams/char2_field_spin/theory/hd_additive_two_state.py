"""Exact four-state support law for a two-prefix observation inner.

At each coordinate sample (c,d1,d2,f) uniformly among nonzero field tuples
with c+d1+f*d2 !=0. Independently across coordinates, compute
    output = s1+f*s2+c*x; s1 += d1*x; s2 += d2*x.
Initial state is (0,0), final state is free. The diagonal c makes the
complete input/output map injective even without the extra determinant
condition; the condition makes the local input/state map invertible too.

The state index is 2*[s1!=0]+[s2!=0]. Values within a support class have
the same transition law by nonzero coordinate scaling of the coefficient
tuple. This helper establishes that local law, not a full code certificate.
"""
from fractions import Fraction as F
import math

import numpy as np


def support_law(old, active, q):
    if old not in range(4) or active not in (0,1) or q < 4:
        raise ValueError('invalid support or field size')
    if not active:
        if old == 0:
            return {(0,0): F(1)}
        if old == 3:
            return {(3,0): F(1,q-1), (3,1): F(q-2,q-1)}
        return {(old,1): F(1)}
    if old == 0:
        return {(3,1): F(1)}
    choices = q*q-3*q+3
    denominator = (q-1)**2*choices
    if old in (1,2):
        # The previously zero coordinate is necessarily live afterward.
        only = 3-old
        counts = {
            only: ((q-1)*choices, (q-1)**2),
            3: ((q-2)*(q-1)*choices, (q-2)**2*(q-1)),
        }
    else:
        counts = {
            0: (choices, 0),
            1: ((q-2)*choices, (q-2)**2),
            2: ((q-2)*choices, (q-2)**2),
            3: ((q-2)**2*choices, (q-2)*(q*q-5*q+7)),
        }
    result = {}
    for new, (total, zeros) in counts.items():
        if zeros:
            result[new,0] = F(zeros,denominator)
        if total != zeros:
            result[new,1] = F(total-zeros,denominator)
    if sum(result.values()) != 1:
        raise ArithmeticError('transition mass is not one')
    return result


def local_transfer(output_mark, q):
    return tuple(np.array([
        [sum(float(p)*output_mark**out for (new,out),p
             in support_law(old,active,q).items() if new == target)
         for target in range(4)] for old in range(4)]) for active in (0,1))


def final_log_moment(logx, logy, n, q, steps=96):
    """Floating positive-vector bound, not an interval certificate."""
    t0,t1 = local_transfer(math.exp(logy),q)
    matrix = t0+math.exp(logx)*t1
    vector = np.ones(4)
    for _ in range(steps):
        following = matrix @ vector
        vector = following/float(following.max())
    upper = float(np.max((matrix @ vector)/vector))
    return n*math.log(upper)+math.log(vector[0]/float(vector.min()))


def exact_probability(n, weight, cutoff, q):
    if not 0 <= weight <= n or not 0 <= cutoff <= n:
        raise ValueError('invalid finite class')
    states = {(0,0,0): F(1)}
    for _ in range(n):
        following = {}
        for (used,old,emitted),mass in states.items():
            for active in (0,1):
                if used+active > weight:
                    continue
                for (new,out),p in support_law(old,active,q).items():
                    if emitted+out > cutoff:
                        continue
                    key = used+active,new,emitted+out
                    following[key] = following.get(key,F(0))+mass*p
        states = following
    return sum((mass for (used,_,_),mass in states.items() if used == weight),
               F(0))/math.comb(n,weight)
