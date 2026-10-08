"""Exact support law and floating moment screen for a rational scalar inner.

For each position, independently sample distinct c,d from F_q minus zero.
Set s'=s+d*x and y=s+c*x. Initial state is zero; final state is free.
The two products are independent of the serial XOR state update. This is
not an accumulator: its inverse generally has dense causal support.

The optional middle-class screen composes a coordinate-major MDS outer,
one weighted prefix, ONE UNIFORM GLOBAL permutation, and this inner.
It reuses a positive-transfer upper bound for the first prefix. It is not
valid for independent regional permutations, and is not a full certificate.
"""
import argparse
from fractions import Fraction as F
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from hd_baa_moment import first_log_moment


def support_law(old, active, q):
    """Return {(new state active, output active): exact probability}."""
    if old not in (0, 1) or active not in (0, 1) or q < 3:
        raise ValueError('invalid support or field size')
    if not active:
        return {(old, old): F(1)}
    if not old:
        return {(1, 1): F(1)}
    return {(0, 1): F(1, q-1), (1, 0): F(1, q-1),
            (1, 1): F(q-3, q-1)}


def local_transfer(output_mark, q):
    return tuple(np.array([
        [sum(float(p)*output_mark**out for (new, out), p
             in support_law(old, active, q).items() if new == target)
         for target in (0, 1)] for old in (0, 1)])
        for active in (0, 1))


def final_log_moment(logx, logy, n, q):
    """Positive-vector upper bound for e0*(T0+x*T1)^n*1.

    The scaling proof is elementary: Mv<=lambda*v implies the stated
    bound for each positive v. The numerical eigenvector proposes v;
    the row-ratio maximum is used instead of trusting an eigenvalue.
    """
    x, y = math.exp(logx), math.exp(logy)
    b = x*y
    c = x*y/(q-1)
    d = y+x*(1+(q-3)*y)/(q-1)
    radical = math.hypot(1-d, 2*math.sqrt(b*c))
    # Solve b*v^2+(1-d)*v-c=0 without subtracting nearly equal numbers.
    v = (2*c/(1-d+radical) if d <= 1 else (d-1+radical)/(2*b))
    upper = max(1+b*v, c/v+d)
    return n*math.log(upper)+max(0., -math.log(v))


def exact_probability(n, weight, cutoff, q):
    """Exact bad-output probability for uniform support of fixed weight.

    Sum over binary input support strings, with each fresh coefficient pair
    marginalized using the exact support law. No independent-output model.
    """
    if not 0 <= weight <= n or not 0 <= cutoff <= n:
        raise ValueError('invalid finite class')
    states = {(0, 0, 0): F(1)}
    for _ in range(n):
        following = {}
        for (used, old, emitted), mass in states.items():
            for active in (0, 1):
                if used+active > weight:
                    continue
                for (new, out), p in support_law(old, active, q).items():
                    if emitted+out > cutoff:
                        continue
                    key = used+active, new, emitted+out
                    following[key] = following.get(key, F(0))+mass*p
        states = following
    return sum((mass for (used, _, _), mass in states.items() if used == weight),
               F(0))/math.comb(n, weight)


def middle_class(alpha, b=128, k=64, q=65536, k_log=20, distance=.4, states=1):
    if states == 1:
        post_moment = final_log_moment
    elif states == 2:
        from hd_additive_two_state import final_log_moment as post_moment
    else:
        raise ValueError('supported scalar state counts are one and two')
    n = (2**k_log//k)*b
    weight = round(alpha*n)
    cutoff = math.ceil(distance*n)-1
    # Unlike an accumulator, there is NO weight<=2*cutoff inverse cutoff.
    if not 0 < weight <= n:
        raise ValueError('middle weight outside full nonzero range')
    first = lambda z: first_log_moment(z, b, k, q, k_log)-weight*z
    first_fit = minimize_scalar(first, bounds=(-math.log(q)-3, 0),
                                method='bounded', options={'maxiter':50})
    logchoose = math.lgamma(n+1)-math.lgamma(weight+1)-math.lgamma(n-weight+1)
    second = lambda marks: (post_moment(*marks, n, q)-weight*marks[0]
                            -cutoff*marks[1]-logchoose)
    starts = ((0., -1.), (math.log(q)/2, -math.log(q)/2),
              (math.log(q), -math.log(q)))
    fits = [minimize(second, start, method='L-BFGS-B', bounds=((-30,30),(-30,0)),
                     options={'maxiter':60, 'ftol':1e-12}) for start in starts]
    best = min(fits, key=lambda r:r.fun)
    logbound = first_fit.fun+best.fun-math.log(q-1)
    return dict(outer=[b,k], q=q, k_log=k_log, n=n, post_state_count=states,
                intermediate_weight=weight,
                cutoff=cutoff, margin_bits=-float(logbound)/math.log(2),
                first_log_mark=float(first_fit.x), final_log_marks=best.x.tolist(),
                first_bits=float(first_fit.fun)/math.log(2),
                final_bits=float(best.fun)/math.log(2),
                first_converged=bool(first_fit.success),
                final_converged=bool(best.success), rigorous=False, whole_code=False,
                route='one uniform global permutation, NOT regional')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outer', type=int, default=128)
    p.add_argument('--states', type=int, choices=(1,2), default=1)
    p.add_argument('--density', type=float, nargs='+', default=[.2,.4,.6,.8,.95])
    args = p.parse_args()
    for alpha in args.density:
        print(json.dumps(middle_class(alpha, args.outer, args.outer//2,
                                      states=args.states)), flush=True)
