"""Directed dense certificates for higher-rate fixed local puncturing.

The full causal inverse and omitted-output completion argument are the
same as the frozen [128,84] checkpoint. This helper changes outer
dimension and retained output count with new, explicitly bound metadata.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

from flint import arb,ctx

from gf8_outer_rate_dense import BALL_EXPONENT,bracket,propose_tilt,rational,lower_integer
from gf8_rate_climb_profile import metadata


def interval_margin(g,first,last,tilt):
    if not 1<=first<=last<=g['input_n']:
        raise ValueError('invalid dense integer cell')
    base=bracket(g['outer'][1],F(first,g['input_n']),F(last,g['input_n']),tilt)
    return ((first*rational(tilt).log()-g['rows']*rational(base).log()
             -128*arb(g['rows']+1).log())/arb(2).log()
            -rational(BALL_EXPONENT)*g['output_n']-8*(g['input_n']-g['output_n']))


def point_frontier(dimension=80,outputs=6,k_log=20,target_bits=40):
    """Point diagnostic only; no monotonicity or interval coverage claim."""
    g=metadata(dimension,outputs,k_log)
    previous=ctx.prec
    try:
        ctx.prec=128
        def score(weight):
            tilt=propose_tilt(dimension,F(weight,g['input_n']),F(weight,g['input_n']))
            return interval_margin(g,weight,weight,tilt)
        if not score(g['input_n'])>=target_bits:
            return dict(dimension=dimension,outputs=outputs,first_point_weight=None,
                        first_point_density=None,rigorous=False,complete=False)
        lo,hi=1,g['input_n']
        while lo<hi:
            middle=(lo+hi)//2
            if score(middle)>=target_bits:
                hi=middle
            else:
                lo=middle+1
        return dict(dimension=dimension,outputs=outputs,first_point_weight=lo,
                    first_point_density=str(F(lo,g['input_n'])),
                    approximate_first_density=lo/g['input_n'],rigorous=False,complete=False)
    finally:
        ctx.prec=previous


def certify(dimension=80,outputs=6,first_density=F(4,5),k_log=20,cell_bits=128,
            initial_grid=200,precision=192,include_witnesses=False):
    return certify_profile(metadata(dimension,outputs,k_log),first_density,cell_bits,
                           initial_grid,precision,include_witnesses)


def certify_profile(g,first_density,cell_bits=128,initial_grid=200,precision=192,
                    include_witnesses=False):
    """The same dense theorem for a supplied fixed or periodic geometry."""
    if not isinstance(first_density,F) or not 0<first_density<=1:
        raise ValueError('exact positive first density required')
    if cell_bits<40 or initial_grid<1 or precision<80:
        raise ValueError('invalid directed replay settings')
    if not F(10,9)**100*2295**10<2**127:
        raise ArithmeticError('exact radius-.1 output-ball inequality failed')
    dimension=g['outer'][1]
    if g['q']!=256 or g['outer'][0]!=128 or not 0<g['output_n']<=g['input_n']:
        raise ValueError('unsupported dense geometry')
    n=g['input_n']
    first=(first_density.numerator*n+first_density.denominator-1)//first_density.denominator
    previous=ctx.prec
    try:
        ctx.prec=precision
        span=(n+initial_grid-1)//initial_grid
        pending=[(lo,min(n,lo+span-1)) for lo in range(first,n+1,span)][::-1]
        entries=[]
        while pending:
            lo,hi=pending.pop()
            tilt=propose_tilt(dimension,F(lo,n),F(hi,n))
            bits=lower_integer(interval_margin(g,lo,hi,tilt))
            if bits>=cell_bits:
                entries.append((lo,hi,tilt,bits))
            elif lo==hi:
                raise ArithmeticError(f'dense weight{lo} has only{bits} directed bits')
            else:
                middle=(lo+hi)//2
                pending.extend(((middle+1,hi),(lo,middle)))
        if (entries[0][0]!=first or entries[-1][1]!=n or
                any(a[1]+1!=b[0] for a,b in zip(entries,entries[1:]))):
            raise ArithmeticError('dense coverage gap')
        bits=min(entry[3] for entry in entries)
        digest=sha256(''.join(f'{lo}:{hi}:{tilt}:{margin}\n'
                             for lo,hi,tilt,margin in entries).encode('ascii')).hexdigest()
        result=dict(**g,weights=[first,n],covered_weights=n-first+1,first_density=str(first_density),
            checked_cells=len(entries),precision_bits=precision,initial_grid=initial_grid,
            target_cell_bits=cell_bits,ball_exponent=str(BALL_EXPONENT),conditioning='(rows+1)^128',
            omitted_symbols=n-g['output_n'],completion_factor_log2=8*(n-g['output_n']),
            margin_bits_lower=bits,failure_upper_dyadic_exponent=bits,
            failure_upper_dyadic_hex=['0x1',hex(1<<bits)],witness_sha256=digest,
            complete=True,rigorous=True,whole_code=False,
            scope='complete dense inverse-input interval, including omitted-output completions',
            dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals')
        if include_witnesses:
            result['witnesses']=[[lo,hi,str(tilt),margin] for lo,hi,tilt,margin in entries]
        return result
    finally:
        ctx.prec=previous


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,nargs='+',default=[80])
    parser.add_argument('--outputs',type=int,default=6)
    parser.add_argument('--first',type=F,default=F(4,5))
    parser.add_argument('--point-frontier',action='store_true')
    args=parser.parse_args()
    for dimension in args.dimension:
        result=(point_frontier(dimension,args.outputs) if args.point_frontier
                else certify(dimension,args.outputs,args.first))
        print(json.dumps(result),flush=True)
