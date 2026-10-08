"""Directed sparse bounds for a fixed periodic packet-output restriction.

The chronological local polynomial product is formed before exponentiation.
No averaging over puncturing phases or additional routing is used.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from flint import arb, arb_poly, ctx
from field128_sparse_certificate import complete_cover, region
from gf8_periodic_climb_profile import metadata
from gf8_punctured_screen import local_matrices
from gf8_rate_climb_sparse import witnesses
from sparse_active_interval_screen import anchor_zones
from state8_sparse_certificate import (
    ball, integer_bits, poly_product, rational_marks, source_maximum,
    terminal_power)


def periodic_polynomial_power(y, tau, pattern, packets, degree):
    """Return coefficients of the ordered width-eight packet product."""
    pattern=tuple(pattern)
    if (not pattern or any(not isinstance(h,int) or not 1<=h<=8 for h in pattern)
            or packets<0 or packets%len(pattern) or degree<0):
        raise ValueError('a valid period must divide the packet count')
    identity=(arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
    factor=identity
    for outputs in pattern:
        local=local_matrices(y,tau,8,outputs,8)
        packet=tuple(arb_poly([ball(local[w][i][j])*math.comb(8,w)
                              for w in range(min(8,degree)+1)])
                     for i,j in ((0,0),(0,1),(1,0),(1,1)))
        factor=poly_product(factor,packet,degree)
    result=identity
    periods=packets//len(pattern)
    while periods:
        if periods&1:
            result=poly_product(result,factor,degree)
        periods>>=1
        if periods:
            factor=poly_product(factor,factor,degree)
    return result


def certify_zone(lower,upper,anchor,point,dimension=84,pattern=(6,6,6,7),
                 k_log=20,cell_bits=80,precision=4096,split=True):
    g=metadata(dimension,pattern,k_log)
    if not 1<=lower<=upper<=g['rows'] or cell_bits<1:
        raise ValueError('invalid active interval')
    y,tau=rational_marks(anchor,point,g['output_n'])
    previous=ctx.prec
    try:
        ctx.prec=precision
        coefficients=periodic_polynomial_power(y,tau,pattern,g['rows']//8,upper)
        output_factor=ball(y)**(-g['output_threshold'])
        target=arb(2)**(-cell_bits)
        cache={}
        def bound(lo,hi):
            if hi not in cache:
                cache[hi]=terminal_power(region(coefficients,hi,g['rows'],tau),128)
            return (source_maximum(lo,hi,g['rows'],128-dimension,tau)
                    *output_factor*cache[hi]*(hi-lo+1))
        stack,entries,total=[(lower,upper)],[],arb(0)
        while stack:
            lo,hi=stack.pop()
            endpoint=bound(lo,hi).upper()
            if not endpoint.is_finite() or not endpoint>0:
                raise ArithmeticError('increase Arb precision')
            if endpoint<=target or not split:
                entries.append(dict(lower=lo,upper=hi,anchor=anchor,
                                    certified_integer_bits=integer_bits(endpoint)))
                total+=endpoint
            elif lo==hi:
                raise ArithmeticError(f'active={lo}: only {integer_bits(endpoint)} bits at anchor={anchor}')
            else:
                mid=(lo+hi)//2
                stack.extend(((mid+1,hi),(lo,mid)))
        complete_cover(entries,lower,upper)
        endpoint=total.upper()
        exact=endpoint.fmpq()
        return dict(anchor=anchor,y=str(y),tau=str(tau),precision_bits=precision,
                    entries=entries,certified_integer_bits=integer_bits(endpoint),
                    failure_upper_dyadic_hex=[hex(int(exact.p)),hex(int(exact.q))])
    finally:
        ctx.prec=previous


def certify(last=1070,dimension=84,pattern=(6,6,6,7),k_log=20,
            cell_bits=80,precision=4096,progress=None):
    g=metadata(dimension,pattern,k_log)
    if not 1<=last<=g['rows']:
        raise ValueError('active endpoint outside row range')
    anchors=witnesses(last,dimension)
    receipts,entries,total=[],[],F(0)
    for lo,hi,anchor in anchor_zones(1,last,anchors):
        receipt=certify_zone(lo,hi,anchor,anchors[anchor],dimension,pattern,k_log,
                             cell_bits,precision)
        receipts.append(receipt)
        entries.extend(receipt['entries'])
        p,q=(int(x,16) for x in receipt['failure_upper_dyadic_hex'])
        total+=F(p,q)
        if progress:
            progress(dict(zone=[lo,hi],anchor=anchor,
                          certified_integer_bits=receipt['certified_integer_bits'],
                          intervals=len(receipt['entries'])))
    complete_cover(entries,1,last)
    previous=ctx.prec
    try:
        ctx.prec=256
        bits=integer_bits(ball(total).upper())
    finally:
        ctx.prec=previous
    digest=sha256(''.join(f"{e['lower']}:{e['upper']}:{e['anchor']}\n" for e in entries).encode('ascii')).hexdigest()
    return dict(**g,active_range=[1,last],complete=True,rigorous=True,whole_code=False,
                margin_bits_lower=bits,interval_count=len(entries),partition_sha256=digest,
                receipts=receipts,
                failure_upper_dyadic_hex=[hex(total.numerator),hex(total.denominator)],
                scope='complete stated active-row interval only')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,default=84)
    parser.add_argument('--pattern',type=int,nargs='+',default=[6,6,6,7])
    parser.add_argument('--last',type=int,default=1070)
    parser.add_argument('--precision',type=int,default=4096)
    parser.add_argument('--cell-bits',type=int,default=80)
    parser.add_argument('--quiet',action='store_true')
    args=parser.parse_args()
    result=certify(args.last,args.dimension,tuple(args.pattern),cell_bits=args.cell_bits,
                   precision=args.precision,
                   progress=None if args.quiet else lambda x:print(json.dumps(x),flush=True))
    del result['receipts'],result['failure_upper_dyadic_hex']
    print(json.dumps(result),flush=True)
