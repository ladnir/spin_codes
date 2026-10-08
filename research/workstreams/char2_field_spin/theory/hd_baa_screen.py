"""Exact q-ary accumulator counts and one-permutation striped-band gates.

No native timing, asymptotic-to-finite inference, or whole-code certificate.
The band gate is a genuine failure-probability lower bound, independent of
the field and both accumulators' fixed nonzero coefficients.
"""
import argparse
from fractions import Fraction as F
import json
import math

import numpy as np
from scipy.stats import hypergeom


def choose(n,k):
    return math.comb(n,k) if 0 <= k <= n else 0


def accumulator_probability(n,w,h,q):
    """Exact P[output weight h | uniform input support w], zero initial state.

    Live inputs are independently uniform nonzero, or the accumulator has
    independent uniform nonzero feedback multipliers. Final state is free.
    """
    if not (n >= 0 and 0 <= w <= n and 0 <= h <= n and q >= 2):
        raise ValueError('invalid accumulator parameters')
    if not w:
        return F(int(h == 0))
    if not h:
        return F(0)
    total=F(0)
    cancel=F(1,q-1)
    survive=F(q-2,q-1)
    for c in range(min(h,n-h,w//2)+1):
        # Closed last run; then open last run. Binomial conventions handle c=0.
        j=w-2*c
        closed=choose(n-h,c)*choose(h-1,c-1)*choose(h-c,j)
        if closed:
            total += closed*cancel**c*survive**j
        j=w-2*c-1
        opened=choose(n-h,c)*choose(h-1,c)*choose(h-c-1,j)
        if opened:
            total += opened*cancel**c*survive**j
    return total/choose(n,w)


def accumulator_iowe(n,w,h,q):
    """Number of q-ary input words of weight w producing output weight h."""
    value=accumulator_probability(n,w,h,q)*choose(n,w)*(q-1)**w
    if value.denominator != 1:
        raise ArithmeticError('IOWE is not integral')
    return value.numerator


def local_transfer(y,q):
    return (np.array([[1.,0.],[0.,y]]),
            np.array([[0.,y],[1/(q-1),(q-2)/(q-1)*y]]))


def double_prefix_transfer(y,q):
    """Exact four-state support law for two consecutive weighted prefixes.

    State is (first-prefix value !=0, second-prefix value !=0). Each stage
    uses its own independent uniform nonzero multiplier at every position.
    Only the second stage's output weight is marked by y.
    """
    scalar=local_transfer(1.,q)
    result=[]
    for active in range(2):
        matrix=np.zeros((4,4))
        for old_first in range(2):
            for old_second in range(2):
                for new_first in range(2):
                    for new_second in range(2):
                        matrix[2*old_first+old_second,2*new_first+new_second] = (
                            scalar[active][old_first,new_first]
                            *scalar[new_first][old_second,new_second]*y**new_second)
        result.append(matrix)
    return tuple(result)


def group_ceiling(b,k_log=20,profile=(1,1)):
    """Exact deterministic parity-check group ceiling at final rate one half.

    Profile records consecutive scalar scans on each side of a permutation.
    An unpermuted block of t scans has inverse column weight at most t+1.
    No distributional assumption or MDS property is needed for this lemma.
    """
    if b < 2 or b%2 or not profile or any(t < 1 for t in profile):
        raise ValueError('invalid rate-half local block or scan profile')
    parity=b//2
    if 2**k_log % parity:
        raise ValueError('complete unpadded rate-half blocks required')
    groups=2**k_log//parity
    degree=math.prod(t+1 for t in profile)
    if degree > groups:
        raise ValueError('degree exceeds number of check groups')
    denominator=choose(groups-1,degree-1)
    low,high=degree,groups
    while low < high:
        middle=(low+high)//2
        if 2*choose(middle-1,degree-1) > denominator:
            high=middle
        else:
            low=middle+1
    cutoff=parity*low+1
    return dict(outer=[b,parity],k_log=k_log,n=2**(k_log+1),
                profile=list(profile),inverse_group_degree_at_most=degree,
                check_groups=groups,checks_per_group=parity,selected_groups=low,
                distance_upper_bound=cutoff,relative_upper_bound=cutoff/2**(k_log+1),
                exact=True,deterministic=True)


def distant_repeat_moment(k,y,q,source_mark=1.):
    """E sum_x source_mark^wt(x) y^wt(A D [x;x]), including x=0.

    Four-state transfer retains the same message support in both halves.
    This bounded evaluator is for small k; a log-scaled power is needed at K20.
    """
    t0,t1=local_transfer(y,q)
    matrix=np.kron(t0,t0)+(q-1)*source_mark*np.kron(t1,t1)
    product=np.linalg.matrix_power(matrix,k)
    return sum(product[a,2*a+b] for a in range(2) for b in range(2))


def band_gate(b,a,k_log=20,distance=F(2,5),terms=8):
    """Rigorous finite lower bound from disjoint bands of a adjacent rows.

    Outer [b,b/2], coordinate-major first scan, scalar nonzero weighted
    accumulator, ONE uniform global permutation, scalar final accumulator.
    Each of the b band exits is zero, deleting b known output coordinates.
    The first `terms` hypergeometric tail terms already give a lower bound.
    """
    if b < 2 or b%2 or a < 3 or terms < 1:
        raise ValueError('even outer length, a>=3, and positive term count required')
    dimension=b//2
    if 2**k_log % dimension:
        raise ValueError('this gate assumes complete source rows without padding')
    rows=2**k_log//dimension
    n=b*rows
    e=(F(distance).numerator*n-1)//F(distance).denominator
    size=b*(a-1)
    threshold=(b-dimension)*a+1
    groups=rows//a
    if not groups or size > n:
        raise ValueError('band exceeds source geometry')
    denominator=choose(n,size)
    p=sum((F(choose(e,t)*choose(n-e,size-t),denominator)
           for t in range(threshold,min(size,e,threshold+terms-1)+1)),F(0))
    failure=groups*p/(1+(groups-1)*p)
    bits=(math.log2(failure.denominator)-math.log2(failure.numerator)) if failure else math.inf
    return dict(outer=[b,dimension],adjacent_rows=a,rows=rows,n=n,
                target_distance=str(F(distance)),suffix=e,container_size=size,
                preacc_subspace_dimension_at_least=a*dimension-b,
                threshold=threshold,disjoint_groups=groups,
                exact_tail_terms=terms,single_group_probability_lower=float(p),
                failure_probability_lower=float(failure),reliability_upper_bits=bits,
                exact_failure_fraction=[str(failure.numerator),str(failure.denominator)],
                exact=True,whole_code_upper_certificate=False)


def select_band(b,k_log=20,distance=F(2,5),last=48):
    """Floating selection only; the returned gate is evaluated exactly."""
    rows=2**k_log//(b//2)
    n=b*rows
    e=(F(distance).numerator*n-1)//F(distance).denominator
    def log_expected(a):
        size=b*(a-1)
        threshold=b//2*a+1
        return math.log(rows//a)+float(hypergeom.logsf(threshold-1,n,e,size))
    a=max(range(3,min(last,rows)+1),key=log_expected)
    return band_gate(b,a,k_log,distance)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outer',type=int,nargs='+',default=[16,32,64,128])
    p.add_argument('--last',type=int,default=48)
    p.add_argument('--fractions',action='store_true')
    p.add_argument('--ceilings',action='store_true')
    args=p.parse_args()
    if args.ceilings:
        for b in args.outer:
            for profile in ((1,1),(2,1),(3,1),(2,2)):
                print(json.dumps(group_ceiling(b,profile=profile)),flush=True)
        raise SystemExit
    for b in args.outer:
        result=select_band(b,last=args.last)
        if not args.fractions:
            del result['exact_failure_fraction']
        print(json.dumps(result),flush=True)
