"""Positive reference bounds for the fixed BCH/RandomStepConv occupation tail.

Binary64 routines are diagnostic only. Prefix envelopes are exact rationals.
"""
from __future__ import annotations
import math
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
sys.path.insert(0,str(GEN/'bch256_q2_envelope_diagnostic_sources'))
import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_bch_q2_envelopes import envelopes
from probe_cumulative_reference_q3 import prefix_envelope

L=8192
B=256
D=209716


def normalized_regions(zero,candidate,length,maximum):
    current=np.zeros((maximum+1,2,2))
    scratch=np.zeros_like(current)
    current[0]=np.eye(2)
    degrees=np.arange(maximum+1,dtype=float)
    for n in range(1,length+1):
        size=min(n,maximum)+1
        old_size=min(n-1,maximum)+1
        updated=scratch[:size]
        updated.fill(0)
        updated[:old_size]=np.matmul(current[:old_size],zero)*((n-degrees[:old_size])/n)[:,None,None]
        updated[1:]+=np.matmul(current[:size-1],candidate)*(degrees[1:size]/n)[:,None,None]
        current,scratch=scratch,current
    return current


def log_choose_array(n,values):
    return np.array([math.lgamma(n+1)-math.lgamma(int(q)+1)-math.lgamma(n-int(q)+1) for q in values])


def dense_moment_log(s,rho,q,p,log_choose):
    zero,active=base.step_matrices(math.exp(-s),22)
    eta=p*rho
    mixed=(1-eta)[:,None,None]*zero+eta[:,None,None]*active
    full=base.log_power(base.log_entries(mixed),L*B)
    log_moment=np.logaddexp(full[:,0,0],full[:,0,1])
    log_point=log_choose+q*np.log(p)
    interior=q<L
    log_point[interior]+=(L-q[interior])*np.log1p(-p[interior])
    return log_moment-B*log_point


def contiguous(values):
    result=[]
    for value in values:
        q=int(value)
        if result and q==result[-1][1]+1:
            result[-1][1]=q
        else:
            result.append([q,q])
    return result


def toy_check():
    import itertools
    zero=np.array(((1.,0.),(.02,.8)))
    active=np.array(((.02,.8),(.02,.8)))
    candidate=.6*zero+.4*active
    actual=normalized_regions(zero,candidate,5,5)
    for q in range(6):
        direct=np.zeros((2,2))
        for selected in itertools.combinations(range(5),q):
            value=np.eye(2)
            for i in range(5):
                value=value@(candidate if i in selected else zero)
            direct+=value/math.comb(5,q)
        assert np.max(np.abs(direct-actual[q]))<1e-14
        p=.37 if q<5 else 1.
        upper=np.linalg.matrix_power((1-p)*zero+p*candidate,5)/(math.comb(5,q)*p**q*(1-p)**(5-q))
        assert np.all(upper+1e-14>=direct)
    # Summation-by-parts domination on an exact finite weight space.
    mass=[0,1,2,1]
    reference=[Fraction(1,8),Fraction(3,8),Fraction(3,8),Fraction(1,8)]
    factor=max(Fraction(sum(mass[:i+1]),sum(reference[:i+1])) for i in range(4))
    for drops in itertools.product(range(3),repeat=4):
        decreasing=[sum(drops[i:]) for i in range(4)]
        assert sum(a*g for a,g in zip(mass,decreasing))<=factor*sum(p*g for p,g in zip(reference,decreasing))
    return dict(all_32_supports_checked=True,all_81_decreasing_toy_functions_checked=True,
                coefficient_conditioning_checked=True)
