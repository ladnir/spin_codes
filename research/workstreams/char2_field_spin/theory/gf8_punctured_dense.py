"""Directed dense-input certificate after fixed seven-of-eight puncturing.

The full causal inner is invertible because its current-input feedthrough
F is invertible. Each retained output vector has exactly 256^lost full
output completions. The old fixed-full-output membership theorem applies
to each completion, after conditioning on state diagonals; input
diagonals and the single regional route remain independent and fresh.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

from flint import arb, ctx

from gf8_outer_rate_dense import BALL_EXPONENT, bracket, propose_tilt, rational, lower_integer
from gf8_punctured_profile import metadata


def interval_margin(g, first, last, tilt):
    if not 1 <= first <= last <= g['input_n']:
        raise ValueError('invalid dense inverse-input interval')
    base = bracket(84, F(first, g['input_n']), F(last, g['input_n']), tilt)
    return ((first*rational(tilt).log()-g['rows']*rational(base).log()
             -128*arb(g['rows']+1).log())/arb(2).log()
            -rational(BALL_EXPONENT)*g['output_n']
            -8*(g['input_n']-g['output_n']))


def certify(first_density=F(9,20), cell_bits=128, initial_grid=200, precision=192,
            include_witnesses=False, k_log=20):
    if not isinstance(first_density, F) or not 0 < first_density <= 1:
        raise ValueError('exact positive first density required')
    if cell_bits < 40 or initial_grid < 1 or precision < 80:
        raise ValueError('invalid directed replay settings')
    if not F(10,9)**100*2295**10 < 2**127:
        raise ArithmeticError('exact radius-.1 output-ball inequality failed')
    g = metadata(k_log)
    n = g['input_n']
    first = (first_density.numerator*n+first_density.denominator-1)//first_density.denominator
    previous = ctx.prec
    try:
        ctx.prec = precision
        span = (n+initial_grid-1)//initial_grid
        pending = [(lo,min(n,lo+span-1)) for lo in range(first,n+1,span)][::-1]
        entries = []
        while pending:
            lo,hi = pending.pop()
            tilt = propose_tilt(84,F(lo,n),F(hi,n))
            bits = lower_integer(interval_margin(g,lo,hi,tilt))
            if bits >= cell_bits:
                entries.append((lo,hi,tilt,bits))
            elif lo == hi:
                raise ArithmeticError(f'dense weight {lo} has only {bits} directed bits')
            else:
                middle = (lo+hi)//2
                pending.extend(((middle+1,hi),(lo,middle)))
        if (entries[0][0] != first or entries[-1][1] != n or
                any(a[1]+1 != b[0] for a,b in zip(entries,entries[1:]))):
            raise ArithmeticError('dense coverage has a gap')
        bits = min(item[3] for item in entries)
        digest = sha256(''.join(f'{lo}:{hi}:{tilt}:{margin}\n'
                               for lo,hi,tilt,margin in entries).encode('ascii')).hexdigest()
        result = dict(**g,weights=[first,n],covered_weights=n-first+1,
            first_density=str(first_density),checked_cells=len(entries),
            precision_bits=precision,initial_grid=initial_grid,target_cell_bits=cell_bits,
            ball_exponent=str(BALL_EXPONENT),conditioning='(rows+1)^128',
            omitted_symbols=g['input_n']-g['output_n'],
            completion_factor_log2=8*(g['input_n']-g['output_n']),
            margin_bits_lower=bits,failure_upper_dyadic_exponent=bits,
            failure_upper_dyadic_hex=['0x1',hex(1<<bits)],witness_sha256=digest,
            complete=True,rigorous=True,whole_code=False,
            scope='complete dense inverse-input interval, with all output completions unioned',
            dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals')
        if include_witnesses:
            result['witnesses'] = [[lo,hi,str(tilt),margin] for lo,hi,tilt,margin in entries]
        return result
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=F,default=F(9,20))
    parser.add_argument('--cell-bits',type=int,default=128)
    parser.add_argument('--show-witnesses',action='store_true')
    args = parser.parse_args()
    print(json.dumps(certify(args.first,args.cell_bits,include_witnesses=args.show_witnesses)),flush=True)
