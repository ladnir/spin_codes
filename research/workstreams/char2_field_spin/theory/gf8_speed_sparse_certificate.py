"""Directed sparse cover for GRS64/40 and the punctured width4/state4 map.

The receipt covers only its stated active-row interval. Local coefficients
are exact rationals; all polynomial and terminal powers use outward Arb
arithmetic. The fixed profile supplies the matrix identity and geometry.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from flint import arb,arb_poly,ctx
from field128_sparse_certificate import complete_cover,region
from gf8_punctured_screen import local_matrices
from gf8_speed_profile import metadata
from sparse_active_interval_screen import anchor_zones,interpolated_anchors
from state8_sparse_certificate import ball,integer_bits,poly_product,rational_marks,source_maximum,terminal_power


ANCHORS={
    1:(5.902081713326922,5.541263545158426),
    2:(5.889930278515905,5.541263545158426),
    4:(6.178089401638577,5.541263545158426),
    8:(5.960047465228139,3.1782902769390806),
    16:(5.876953454328387,2.5527284957990455),
    32:(5.809891536014135,2.3228594248421537),
    64:(5.774315897049173,2.2263565643290972),
    96:(5.767137870007008,2.201106692133474),
    192:(5.766363296752179,2.1861070240560623),
    256:(5.768725049531589,2.18626193956568),
    512:(5.781122910266329,2.2008322525372517),
}


def witnesses(last):
    raw=dict(ANCHORS)
    for anchor in (768,1024,1536,last):
        if anchor>max(ANCHORS):
            raw[anchor]=ANCHORS[max(ANCHORS)]
    # Fixed y is sensitive to its anchor near the shallow sparse dip. Use
    # a separate witness for every class there, without changing the bound.
    ordered=sorted(raw)
    for anchor in range(32,min(last,192)+1):
        if anchor in raw:
            continue
        left=max(a for a in ordered if a<anchor)
        right=min(a for a in ordered if a>anchor)
        t=math.log(anchor/left)/math.log(right/left)
        raw[anchor]=tuple(a+(b-a)*t for a,b in zip(raw[left],raw[right]))
    return interpolated_anchors(interpolated_anchors(raw))


def polynomial_power(local,packets,degree):
    if len(local)!=5 or packets<0 or degree<0:
        raise ValueError('width-four local matrices and nonnegative exponent/degree required')
    result=(arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
    factor=tuple(arb_poly([ball(local[w][i][j])*math.comb(4,w)
                          for w in range(min(4,degree)+1)])
                 for i,j in ((0,0),(0,1),(1,0),(1,1)))
    while packets:
        if packets&1:
            result=poly_product(result,factor,degree)
        packets>>=1
        if packets:
            factor=poly_product(factor,factor,degree)
    return result


def certify_zone(lower,upper,anchor,point,dimension=40,k_log=20,
                 cell_bits=50,precision=4096,split=True):
    g=metadata(dimension,k_log)
    length,dimension=g['outer']
    if not 1<=lower<=upper<=g['rows'] or cell_bits<1:
        raise ValueError('invalid active interval')
    y,tau=rational_marks(anchor,point,g['output_n'])
    previous=ctx.prec
    try:
        ctx.prec=precision
        local=local_matrices(y,tau,g['width'],g['outputs'],g['state'])
        coefficients=polynomial_power(local,g['rows']//g['width'],upper)
        output_factor=ball(y)**(-g['output_threshold'])
        target=arb(2)**(-cell_bits)
        cache={}
        def bound(lo,hi):
            if hi not in cache:
                cache[hi]=terminal_power(region(coefficients,hi,g['rows'],tau),length)
            return source_maximum(lo,hi,g['rows'],length-dimension,tau)*output_factor*cache[hi]*(hi-lo+1)
        stack,entries,total=[(lower,upper)],[],arb(0)
        while stack:
            lo,hi=stack.pop()
            endpoint=bound(lo,hi).upper()
            if not endpoint.is_finite() or not endpoint>0:
                raise ArithmeticError('increase Arb precision')
            if endpoint<=target or not split:
                entries.append(dict(lower=lo,upper=hi,anchor=anchor,certified_integer_bits=integer_bits(endpoint)))
                total+=endpoint
            elif lo==hi:
                raise ArithmeticError(f'active={lo}: only {integer_bits(endpoint)} bits at anchor={anchor}')
            else:
                mid=(lo+hi)//2
                stack.extend(((mid+1,hi),(lo,mid)))
        complete_cover(entries,lower,upper)
        endpoint=total.upper()
        exact=endpoint.fmpq()
        return dict(anchor=anchor,y=str(y),tau=str(tau),precision_bits=precision,entries=entries,
                    certified_integer_bits=integer_bits(endpoint),
                    failure_upper_dyadic_hex=[hex(int(exact.p)),hex(int(exact.q))])
    finally:
        ctx.prec=previous


def certify(last=2015,dimension=40,k_log=20,cell_bits=50,precision=4096,progress=None):
    g=metadata(dimension,k_log)
    if not 1<=last<=g['rows']:
        raise ValueError('active endpoint outside row range')
    anchors=witnesses(last)
    receipts,entries,total=[],[],F(0)
    for lo,hi,anchor in anchor_zones(1,last,anchors):
        receipt=certify_zone(lo,hi,anchor,anchors[anchor],dimension,k_log,cell_bits,precision)
        receipts.append(receipt)
        entries.extend(receipt['entries'])
        p,q=(int(x,16) for x in receipt['failure_upper_dyadic_hex'])
        total+=F(p,q)
        if progress:
            progress(dict(zone=[lo,hi],anchor=anchor,certified_integer_bits=receipt['certified_integer_bits'],
                          intervals=len(receipt['entries'])))
    complete_cover(entries,1,last)
    previous=ctx.prec
    try:
        ctx.prec=256
        bits=integer_bits(ball(total).upper())
    finally:
        ctx.prec=previous
    digest=sha256(''.join(f"{e['lower']}:{e['upper']}:{e['anchor']}\n" for e in entries).encode('ascii')).hexdigest()
    return dict(**g,active_range=[1,last],complete=True,rigorous=True,whole_code=False,margin_bits_lower=bits,
                interval_count=len(entries),partition_sha256=digest,receipts=receipts,
                failure_upper_dyadic_hex=[hex(total.numerator),hex(total.denominator)],
                scope='complete stated active-row interval only')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,choices=(39,40),default=40)
    parser.add_argument('--k-log',type=int,choices=(16,18,20),default=20)
    parser.add_argument('--last',type=int,default=2015)
    parser.add_argument('--cell-bits',type=int,default=50)
    parser.add_argument('--precision',type=int,default=4096)
    parser.add_argument('--quiet',action='store_true')
    args=parser.parse_args()
    result=certify(args.last,args.dimension,args.k_log,args.cell_bits,args.precision,
                   progress=None if args.quiet else lambda x:print(json.dumps(x),flush=True))
    del result['receipts'],result['failure_upper_dyadic_hex']
    print(json.dumps(result),flush=True)
