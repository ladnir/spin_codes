"""Directed dense-class receipts for the GF256 GRS[128,d] rate sweep.

The optimizer only proposes rational tilts. Exact endpoint brackets and Arb
outward rounding decide acceptance. Existing certificates remain untouched.
This is one component certificate, not a whole-code distance guarantee.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from flint import arb, ctx


CONSTRUCTION = 'symmetric-cauchy-w16-s8-eta2-regional-nonzero-diagonals'
BALL_EXPONENT = F(127,100)


def geometry(dimension, k_log=20):
    if not 1<=dimension<128 or not 1<=k_log<=30:
        raise ValueError('invalid GRS dimension or requested power of two')
    rows=64*((2**k_log+64*dimension-1)//(64*dimension))
    n=128*rows
    return dict(outer=[128,dimension],q=256,width=16,state=8,k_log=k_log,
                requested_k=2**k_log,padded_k=dimension*rows,rows=rows,n=n,
                output_threshold=(n+9)//10-1,construction=CONSTRUCTION)


def rational(value):
    value=F(value)
    return arb(value.numerator)/value.denominator


def bracket(dimension, lower, upper, tilt):
    """Uniform endpoint bracket for every real density in [lower,upper]."""
    lower,upper,tilt=map(F,(lower,upper,tilt))
    if not 1<=dimension<128 or not 0<=lower<=upper<=1 or tilt<1:
        raise ValueError('invalid endpoint witness')
    return (1-lower)**128+F(1,256**(128-dimension))*(1+(F(256,255)*tilt-1)*upper)**128


def propose_tilt(dimension, lower, upper):
    """Floating maximizer of the endpoint formula; never an acceptance test."""
    lower,upper=float(lower),float(upper)
    if not 0<lower<=upper<=1:
        raise ValueError('invalid proposal interval')
    parity=128-dimension
    def evaluate(log_tilt):
        tilt=math.exp(log_tilt)
        value=1+(256/255*tilt-1)*upper
        x=128*math.log1p(-lower) if lower<1 else -math.inf
        y=-parity*math.log(256)+128*math.log(value)
        maximum=max(x,y)
        log_bracket=maximum+math.log(math.exp(x-maximum)+math.exp(y-maximum))
        derivative=lower-upper*(256/255)*tilt*math.exp(
            -parity*math.log(256)+127*math.log(value)-log_bracket)
        return (lower*log_tilt-log_bracket/128)/math.log(2),derivative
    if evaluate(0)[1]<=0:
        return F(1)
    left,right=0.,40.
    for _ in range(80):
        middle=(left+right)/2
        if evaluate(middle)[1]>0:left=middle
        else:right=middle
    return F(math.exp((left+right)/2)).limit_denominator(10000)


def interval_margin(g, first, last, tilt):
    """Arb interval enclosing a uniform negative log2 failure bound."""
    if not 1<=first<=last<=g['n']:
        raise ValueError('invalid integer dense interval')
    base=bracket(g['outer'][1],F(first,g['n']),F(last,g['n']),tilt)
    log2=arb(2).log()
    return (first*rational(tilt).log()-g['rows']*rational(base).log()
            -128*arb(g['rows']+1).log())/log2-rational(BALL_EXPONENT)*g['n']


def lower_integer(value):
    bound=value.lower().fmpq()
    return int(bound.p//bound.q)


def certify(dimension, k_log=20, first_density=F(6,25), cell_bits=64,
            initial_grid=200, precision=192, include_witnesses=False):
    first_density=F(first_density)
    if not 0<first_density<=1 or cell_bits<40 or initial_grid<1 or precision<80:
        raise ValueError('invalid dense replay parameters')
    if not F(10,9)**100*2295**10<2**127:
        raise ArithmeticError('failed exact radius-.1 output-ball comparison')
    g=geometry(dimension,k_log)
    first=(first_density.numerator*g['n']+first_density.denominator-1)//first_density.denominator
    old_precision=ctx.prec
    try:
        ctx.prec=precision
        span=(g['n']+initial_grid-1)//initial_grid
        pending=[(lo,min(g['n'],lo+span-1)) for lo in range(first,g['n']+1,span)][::-1]
        entries=[]
        while pending:
            lo,hi=pending.pop()
            tilt=propose_tilt(dimension,F(lo,g['n']),F(hi,g['n']))
            bits=lower_integer(interval_margin(g,lo,hi,tilt))
            if bits>=cell_bits:
                entries.append((lo,hi,tilt,bits))
            elif lo==hi:
                raise ArithmeticError(f'dense point weight {lo} has only {bits} directed bits')
            else:
                middle=(lo+hi)//2
                pending.extend(((middle+1,hi),(lo,middle)))
        if (entries[0][0]!=first or entries[-1][1]!=g['n'] or
            any(a[1]+1!=b[0] for a,b in zip(entries,entries[1:]))):
            raise ArithmeticError('gap in dense integer coverage')
        # Each cell certifies a bound for EVERY fixed forbidden output whose
        # inverse input has that density. One uniform output union is used;
        # the number of cells is not an additional probability multiplier.
        bits=min(entry[3] for entry in entries)
        text=''.join(f'{lo}:{hi}:{tilt}:{margin}\n' for lo,hi,tilt,margin in entries)
        result=dict(**g,weights=[first,g['n']],first_density=str(first_density),
                    covered_weights=g['n']-first+1,checked_cells=len(entries),
                    precision_bits=precision,initial_grid=initial_grid,target_cell_bits=cell_bits,
                    ball_exponent=str(BALL_EXPONENT),conditioning='(rows+1)^128',
                    margin_bits_lower=bits,failure_upper_dyadic_exponent=bits,
                    witness_sha256=sha256(text.encode('ascii')).hexdigest(),
                    complete=True,rigorous=True,whole_code_certified=False)
        if include_witnesses:
            result['witnesses']=[[lo,hi,str(tilt),margin] for lo,hi,tilt,margin in entries]
        return result
    finally:
        ctx.prec=old_precision


def point_screen(dimension,k_log=20,target_bits=40):
    """Approximate frontier of the same bound, explicitly not a certificate."""
    g=geometry(dimension,k_log)
    previous=ctx.prec
    try:
        ctx.prec=128
        left,right=F(1,100),F(99,100)
        right_weight=(right.numerator*g['n'])//right.denominator
        right_tilt=propose_tilt(dimension,F(right_weight,g['n']),F(right_weight,g['n']))
        if not interval_margin(g,right_weight,right_weight,right_tilt)>=target_bits:
            return dict(**g,approximate_first_density=None,approximate_first_weight=None,
                        rigorous=False,whole_code_certified=False,
                        scope='no passing point at the upper search endpoint .99')
        for _ in range(40):
            midpoint=(left+right)/2
            weight=max(1,(midpoint.numerator*g['n'])//midpoint.denominator)
            tilt=propose_tilt(dimension,F(weight,g['n']),F(weight,g['n']))
            margin=interval_margin(g,weight,weight,tilt)
            if margin>=target_bits:right=midpoint
            else:left=midpoint
        return dict(**g,approximate_first_density=float(right),
                    approximate_first_weight=math.ceil(float(right)*g['n']),
                    rigorous=False,whole_code_certified=False)
    finally:
        ctx.prec=previous


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,nargs='+',default=[88,90,92,96])
    parser.add_argument('--k-log',type=int,nargs='+',default=[20])
    parser.add_argument('--first',type=F,default=F(6,25))
    parser.add_argument('--cell-bits',type=int,default=64)
    parser.add_argument('--point-screen',action='store_true')
    parser.add_argument('--show-witnesses',action='store_true')
    args=parser.parse_args()
    for k_log in args.k_log:
        for dimension in args.dimension:
            if args.point_screen:result=point_screen(dimension,k_log)
            else:result=certify(dimension,k_log,args.first,args.cell_bits,
                                include_witnesses=args.show_witnesses)
            print(json.dumps(result),flush=True)
