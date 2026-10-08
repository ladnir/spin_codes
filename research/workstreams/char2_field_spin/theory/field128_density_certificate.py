"""Directed K20 middle-input certificate for GRS[128,87], state8/width16.

One regional route, the unchanged symmetric Cauchy maps, eta2, independent
nonzero input/state diagonals. Covers every integer W in [.03N,.24N), not
the whole code. Replay uses frozen rational marks and never an optimizer.
The directed local transfer is reused from the certified identical inner.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_density_certificate import (integer_span, refined_cells,
                                        retarget_witness, transfer)
from subfield_certificate import HIGH, LOW, ZERO, integer_margin


LOW.Emax = HIGH.Emax = 1_000_000_000
LOW.Emin = HIGH.Emin = -1_000_000_000

LENGTH, DIMENSION, PARITY = 128, 87, 41
WIDTH, STATE, FIELD = 16, 8, 256
WITNESSES = tuple((F(density), Witness(*(F(value) for value in values)))
    for density, values in (
        ('3/100', ('2674/8821', '5272/7271', '39826/5971', '22721/5865')),
        ('3/50', ('6204/9799', '4843/8725', '56642/6625', '5950/3389')),
        ('2/25', ('7039/8486', '263/540', '40619/4011', '8792/7225')),
        ('1/10', ('3098/2787', '2978/7217', '89303/7064', '1053/1205')),
        ('3/25', ('8118/5429', '25/74', '96787/6149', '2159/3328')),
        ('7/50', ('22163/9020', '58/269', '185051/8846', '3435/7207')),
        ('9/50', ('93301/5502', '59/3229', '406462/9959', '1869/7429')),
        ('6/25', ('1356591/8030', '8/5339', '255', '298/4939'))))


def weight_geometry(k_log=20):
    period = 64
    rows = ((2**k_log+DIMENSION*period-1)//(DIMENSION*period))*period
    n = rows*LENGTH
    return dict(k_log=k_log, length=LENGTH, dimension=DIMENSION, parity=PARITY,
                rows=rows, k=DIMENSION*rows, n=n, threshold=(n+9)//10-1)


def cell_upper(lo, hi, n, rows, witness, prefactor):
    """Sum bound for every integer input weight in the inclusive cell."""
    if not 0 < lo <= hi < n:
        raise ValueError('invalid integer input-weight cell')
    source = row_upper(F(lo,n), F(hi,n), witness, length=LENGTH, parity=PARITY)
    base = witness.x*witness.tilt
    endpoint = lo if base >= 1 else hi
    return (hi-lo+1)*prefactor*source**rows/interval(base)**endpoint


def prefactor(witness, geometry):
    return (interval(geometry['rows']+1)**LENGTH
            *terminal_power(transfer(witness),geometry['n']//WIDTH)
            /(interval(witness.tau)*interval(witness.y)**geometry['threshold']))


def certify(k_log=20, first=F(3,100), last=F(6,25), initial_grid=1000,
            cell_bits=256, witnesses=WITNESSES):
    if k_log != 20:
        raise ValueError('this frozen checkpoint certifies K20 only')
    if initial_grid < 1 or cell_bits < 40 or not witnesses:
        raise ValueError('invalid directed replay parameters')
    if any(not isinstance(anchor,F) or not 0<anchor<1 or
           any(not isinstance(getattr(w,name),F) for name in ('x','y','tau','tilt'))
           for anchor,w in witnesses):
        raise ValueError('directed replay requires exact Fraction witnesses and anchors')
    g = weight_geometry(k_log)
    first_weight, last_weight = integer_span(g['n'],first,last)
    prefactors = tuple(prefactor(w,g) for _,w in witnesses)

    def evaluate(lo,hi):
        midpoint = F(lo+hi,2*g['n'])
        candidates = sorted(range(len(witnesses)),
                            key=lambda i: abs(witnesses[i][0]-midpoint))[:3]
        options = [(i,cell_upper(lo,hi,g['n'],g['rows'],
                    retarget_witness(*witnesses[i],midpoint),prefactors[i]))
                   for i in candidates]
        return min(options,key=lambda item:item[1].hi)

    maximum = interval(F(1,2**cell_bits)).lo
    span = (g['n']+initial_grid-1)//initial_grid
    total, next_weight, cells, largest = ZERO, first_weight, 0, ZERO.hi
    uses, smallest_width, largest_width = [0]*len(witnesses), None, 0
    partition = sha256()
    for lo,hi,index,bound in refined_cells(first_weight,last_weight,span,evaluate,maximum):
        if lo != next_weight or hi < lo:
            raise ArithmeticError('density cells have a gap, overlap, or reversed endpoints')
        next_weight = hi+1
        total += bound
        largest = max(largest,bound.hi)
        cells += 1
        uses[index] += 1
        width = hi-lo+1
        smallest_width = width if smallest_width is None else min(smallest_width,width)
        largest_width = max(largest_width,width)
        partition.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight != last_weight+1:
        raise ArithmeticError('density cells do not reach the required endpoint')
    margin = integer_margin(total.hi)
    if margin < 40:
        raise ArithmeticError(f'complete class margin below40: {margin}')
    return dict(construction='symmetric-Cauchy width16/state8 fresh diagonals',
                outer=[LENGTH,DIMENSION],parity=PARITY,width=WIDTH,state=STATE,q=FIELD,
                rate=str(F(DIMENSION,LENGTH)),distance='1/10',k_log=k_log,
                requested_k=2**k_log,k=g['k'],n=g['n'],rows=g['rows'],
                output_threshold=g['threshold'],minimum_input_weight=first_weight,
                maximum_input_weight=last_weight,covered_weights=last_weight-first_weight+1,
                density_min=str(first),density_max_exclusive=str(last),checked_cells=cells,
                smallest_cell_weights=smallest_width,largest_cell_weights=largest_width,
                partition_sha256=partition.hexdigest(),witness_cell_counts=uses,
                failure_upper=str(total.hi),class_margin_bits_lower=margin,
                failure_upper_dyadic_exponent=margin,
                minimum_cell_margin_bits_lower=integer_margin(largest),
                precision_digits=80,rigorous=True,whole_code=False,
                scope='complete inverse-input-weight interval only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-grid',type=int,default=1000)
    parser.add_argument('--cell-bits',type=int,default=256)
    args = parser.parse_args()
    print(json.dumps(certify(initial_grid=args.initial_grid,cell_bits=args.cell_bits)),flush=True)
