"""Directed sparse cover for the fixed7of8 full-Cauchy [128,84] candidate.

Every accepted interval uses exact rational local matrices and outward
Arb arithmetic. Routes consume input_n coordinates; distance and the
output mark use output_n after fixed local puncturing. This component
certifies only its explicit active-row range, not the whole code.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

from flint import arb,ctx

from field128_sparse_certificate import complete_cover,region
from gf8_packet_sparse_certificate import polynomial_power
from gf8_punctured_screen import local_matrices
from gf8_punctured_profile import metadata,CONSTRUCTION
from sparse_active_interval_screen import anchor_zones,interpolated_anchors
from state8_sparse_certificate import ball,integer_bits,rational_marks,source_maximum,terminal_power


ANCHORS={
    1:(6.582490666673327,5.541263545158426),
    2:(6.575740080595776,5.541263545158426),
    4:(6.572338631885011,5.541263545158426),
    8:(6.813835865834685,5.517059244781045),
    16:(6.56182951888675,3.0885713532164702),
    32:(6.488077110299875,2.372657676497856),
    64:(6.422781063050977,2.094529384846494),
    128:(6.377552482340238,1.9923267723029066),
    256:(6.3655320014320464,1.9654282320856027),
    512:(6.383178863602478,2.0032797880840616),
}


def certify_zone(lower,upper,anchor,point,cell_bits=64,precision=4096,split=True):
    g=metadata()
    if not 1<=lower<=upper<=g['rows'] or cell_bits<1:
        raise ValueError('invalid active interval')
    y,tau=rational_marks(anchor,point,g['output_n'])
    previous=ctx.prec
    try:
        ctx.prec=precision
        local=local_matrices(y,tau,8,7,8)
        coefficients=polynomial_power(local,g['rows']//8,upper)
        output_factor=ball(y)**(-g['output_threshold'])
        target=arb(2)**(-cell_bits)
        cache={}
        def bound(lo,hi):
            if hi not in cache:
                cache[hi]=terminal_power(region(coefficients,hi,g['rows'],tau),128)
            return source_maximum(lo,hi,g['rows'],44,tau)*output_factor*cache[hi]*(hi-lo+1)
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


def certify(last=713,cell_bits=80,precision=4096,progress=None):
    g=metadata()
    if not 1<=last<=g['rows']:
        raise ValueError('active endpoint outside row range')
    raw=dict(ANCHORS)
    for anchor in (512,768,1024,last):
        if anchor>max(ANCHORS):
            raw[anchor]=ANCHORS[max(ANCHORS)]
    anchors=interpolated_anchors(interpolated_anchors(raw))
    receipts,entries,total=[],[],F(0)
    for lo,hi,anchor in anchor_zones(1,last,anchors):
        receipt=certify_zone(lo,hi,anchor,anchors[anchor],cell_bits,precision)
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
    parser.add_argument('--last',type=int,default=713)
    parser.add_argument('--precision',type=int,default=4096)
    parser.add_argument('--cell-bits',type=int,default=80)
    parser.add_argument('--quiet',action='store_true')
    args=parser.parse_args()
    result=certify(args.last,args.cell_bits,args.precision,
                   None if args.quiet else lambda x:print(json.dumps(x),flush=True))
    del result['receipts'],result['failure_upper_dyadic_hex']
    print(json.dumps(result),flush=True)
