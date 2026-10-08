"""Bounded source-row diagnostics for the preserved vector regional envelope.

Track one row degree exactly under the positive relaxed support measure.
The all-row minimum-weight ceiling uses a union bound, not row independence.
No optimization, native code, timing, or certificate generation occurs here.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import json
import math

import numpy as np
from scipy.special import logsumexp

from hd_vector_region import (SCREEN_WITNESSES,IDENTITY,ZERO,mm,mv,
                             exact_coefficients,log_coefficients,screen_bound,local)
from cauchy_firstmoment_screen import log_choose
from reduced_distance_certificate import log_polynomial_product


def exact_row_degree(active,length,rows,y,q,state,tau):
    """Small exact row-degree polynomial, before the common source prefactor."""
    if not 1<=active<=rows or not 1<=length<=24:
        raise ValueError('small positive source geometry required')
    raw=exact_coefficients(rows,y,q,state,tau)
    kernels=[tuple(tuple(v/comb(rows,h) for v in row) for row in matrix)
             for h,matrix in enumerate(raw)]
    matrices=[]
    for x in (0,1):
        matrices.append(tuple(tuple(sum((comb(active-1,h)*F(tau)**(h+x)*kernels[h+x][i][j]
                                              for h in range(active)),F(0))
                                    for j in range(2)) for i in range(2)))
    coefficients=[IDENTITY]
    for _ in range(length):
        following=[ZERO]*(len(coefficients)+1)
        for degree,old in enumerate(coefficients):
            for x,factor in enumerate(matrices):
                term=mm(old,factor)
                following[degree+x]=tuple(tuple(following[degree+x][i][j]+term[i][j]
                                                for j in range(2)) for i in range(2))
        coefficients=following
    return tuple(sum(matrix[0]) for matrix in coefficients)


def log_row_degree(active,length,rows,y,q,state,tau):
    raw=log_coefficients(rows,active,y,q,state,tau)
    kernels=raw-np.asarray([log_choose(rows,h) for h in range(active+1)])[:,None,None]
    marks=np.asarray([log_choose(active-1,h)+h*math.log(tau) for h in range(active)])
    factors=np.stack((logsumexp(kernels[:active]+marks[:,None,None],axis=0),
                      logsumexp(kernels[1:active+1]+marks[:,None,None]+math.log(tau),axis=0)))
    answer=np.asarray([[[0.,-math.inf],[-math.inf,0.]]])
    exponent=length
    while exponent:
        if exponent&1:answer=log_polynomial_product(answer,factors,length)
        exponent>>=1
        if exponent:factors=log_polynomial_product(factors,factors,length)
    return logsumexp(answer[:,0,:],axis=1)


def diagnose(active,point,length=128,dimension=64,k_log=20,q=65536,state=8):
    rows=(2**k_log+dimension-1)//dimension
    n=rows*length
    y=math.exp(-math.exp(point[0])*active/n)
    tau=math.exp(point[1])
    degrees=log_row_degree(active,length,rows,y,q,state,tau)
    total=float(logsumexp(degrees))
    probabilities=np.exp(degrees-total)
    parity=length-dimension
    bad=float(sum(probabilities[:parity+1]))
    mean=float(sum(w*p for w,p in enumerate(probabilities)))
    variance=float(sum((w-mean)**2*p for w,p in enumerate(probabilities)))
    base=-screen_bound(active,point,length,dimension,k_log,q,state)
    one_gain=-math.log1p(-bad)/math.log(2) if bad<1 else math.inf
    lower=max(0.,1-active*bad)
    all_ceiling=-math.log(lower)/math.log(2) if lower else None
    return dict(outer=[length,dimension],active=active,q=q,state=state,k_log=k_log,
                log_parameters=list(point),base_margin_bits=base,
                row_degree_mean=mean,row_degree_variance=variance,
                one_row_bad_probability=bad,one_row_filter_gain_bits=one_gain,
                all_rows_good_probability_lower=lower,
                all_row_filter_gain_ceiling_bits=all_ceiling,
                best_possible_margin_at_same_marks=base+all_ceiling if lower else None,
                rigorous=False,whole_code=False,scope='fixed-mark floating diagnostic of source relaxation')


def exact_all_row_filter(active,length,rows,y,q,state,tau,minimum):
    """Tiny independent enumeration for tests, not a scalable source algorithm."""
    if active*length>16:raise ValueError('at most16 source support bits')
    raw=exact_coefficients(rows,y,q,state,tau)
    kernels=[tuple(tuple(v/comb(rows,h) for v in row) for row in matrix)
             for h,matrix in enumerate(raw)]
    row_degrees=[F(0)]*(length+1)
    total=good=F(0)
    for pattern in product((0,1),repeat=active*length):
        v=(F(1),F(0))
        for j in range(length):
            h=sum(pattern[i*length+j] for i in range(active))
            v=mv(v,kernels[h])
        weight=sum(v)*F(tau)**sum(pattern)
        degrees=[sum(pattern[i*length:(i+1)*length]) for i in range(active)]
        row_degrees[degrees[0]]+=weight
        total+=weight
        if min(degrees)>=minimum:good+=weight
    return total,good,tuple(row_degrees)


def return_rank_charge(rank_mod_input,q,state,tau):
    """Exact return moment retaining only a required-live scalar input.

    h=dim(span(S_1,...,S_r,X)/span(X)). Additional live histories can only
    reduce this value. A full rank-state transition system is not supplied.
    """
    if not 0<=rank_mod_input<=state or q<2 or not 1<=tau<=q:
        raise ValueError('invalid conditional return-rank parameters')
    return F(q)**(rank_mod_input-state)*F(tau)**(-rank_mod_input)


def directed_ceiling(active,point,length=128,dimension=64,k_log=20,q=65536,state=8,precision=4096):
    """Certify that exact per-row truncation cannot close these fixed marks.

    This bounds a proof expression, not the encoder's actual failure rate.
    The outward-rounded computation only handles the two saved diagnostics.
    """
    from flint import arb,arb_poly,ctx
    from state8_sparse_certificate import ball,poly_product
    if active not in (16,64) or (length,dimension,k_log,q,state)!=(128,64,20,65536,8):
        raise ValueError('directed diagnostic is fixed to the two saved classes')
    rows=2**k_log//dimension;n=rows*length
    threshold=(2*n+4)//5-1
    y=F(math.exp(-math.exp(point[0])*active/n)).limit_denominator(10**12)
    tau=F(math.exp(point[1])).limit_denominator(10**12)
    old=ctx.prec
    try:
        ctx.prec=precision
        def power(factor,exponent,degree):
            result=(arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
            while exponent:
                if exponent&1:result=poly_product(result,factor,degree)
                exponent>>=1
                if exponent:factor=poly_product(factor,factor,degree)
            return result
        transfer=local(y,q,state,tau)
        raw=power(tuple(arb_poly([ball(transfer[x][i][j]) for x in (0,1)])
                        for i,j in ((0,0),(0,1),(1,0),(1,1))),rows,active)
        factors=[]
        for polynomial in raw:
            coefficients=[]
            for x in (0,1):
                coefficients.append(sum((polynomial[h+x]*ball(F(comb(active-1,h),comb(rows,h+x))*tau**(h+x))
                                         for h in range(active)),arb(0)))
            factors.append(arb_poly(coefficients))
        degrees=power(tuple(factors),length,length)
        mass=[degrees[0][w]+degrees[1][w] for w in range(length+1)]
        total=sum(mass,arb(0));bad=sum(mass[:length-dimension+1],arb(0))
        bad_upper=(bad/total).upper()
        good_lower=(arb(1)-active*bad_upper).lower()
        if not good_lower>0:raise ArithmeticError('precision or union bound insufficient')
        gain_upper=(-good_lower.log()/arb(2).log()).upper()
        gain_limit=F(1,1<<(42 if active==16 else 16))
        if not gain_upper<=ball(gain_limit):raise ArithmeticError('claimed truncation ceiling not proved')
        base=arb(comb(rows,active))*total/ball(tau)**((length-dimension)*active+1)/ball(y)**threshold
        filtered_lower=(base.lower()*good_lower).lower()
        lower_bits=math.floor(float(filtered_lower.log()/arb(2).log()))
        while not filtered_lower>=arb(2)**lower_bits:lower_bits-=1
        if lower_bits<0:raise ArithmeticError('fixed witness might close after truncation')
        return dict(outer=[length,dimension],active=active,q=q,state=state,k_log=k_log,
                    y=str(y),tau=str(tau),precision_bits=precision,
                    row_bad_probability_upper=bad_upper.str(24),
                    all_rows_good_probability_lower=good_lower.str(24),
                    all_row_filter_gain_bits_upper=gain_upper.str(24),
                    certified_gain_ceiling_bits=str(gain_limit),
                    filtered_proof_expression_lower_log2=lower_bits,
                    rigorous=True,whole_code=False,actual_failure_probability_claim=False,
                    scope='fixed rational marks: exact per-row source truncation cannot close this proof expression')
    finally:
        ctx.prec=old


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directed',action='store_true')
    args=parser.parse_args()
    for length,dimension,active,is_rank,point in SCREEN_WITNESSES:
        if is_rank and length==128 and active in (16,64) and len(point)==2:
            function=directed_ceiling if args.directed else diagnose
            print(json.dumps(function(active,point,length,dimension)),flush=True)
