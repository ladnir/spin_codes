"""Outer-dimension sweep and directed middle-class replay for GF256 SPIN.

The length128 outer, width16/state8 symmetric inner, and independent
regional route are fixed. Only the outer dimension and natural padding
change. Floating density samples are not certificates. Directed replay
uses exact frozen witnesses and covers every weight in its interval.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from cauchy_firstmoment_screen import geometry
from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_density_certificate import integer_span, refined_cells, retarget_witness, transfer
from subfield_certificate import ZERO, integer_margin
from symmetric_density_screen import search as control_search


CONSTRUCTION = 'symmetric-cauchy-w16-s8-eta2-regional-nonzero-diagonals'

WITNESSES = tuple((F(density), Witness(*(F(value) for value in values)))
    for density, values in (
        ('3/100', ('569/1882', '6709/9250', '57961/8745', '28327/7546')),
        ('1/25', ('2219/5322', '5421/8273', '64024/8987', '573/209')),
        ('3/50', ('1253/1990', '2180/3921', '9247/1091', '16729/9830')),
        ('2/25', ('3504/4255', '4289/8783', '18661/1863', '629/533')),
        ('1/10', ('10543/9592', '4105/9896', '40117/3221', '8359/9853')),
        ('3/25', ('9005/6118', '1211/3553', '14513/940', '4847/7685')),
        ('7/50', ('2296/953', '860/3937', '132917/6497', '1749/3772')),
        ('9/50', ('151091/9369', '158/8267', '5099/131', '1342/5421')),
        ('6/25', ('1052497/6230', '8/5339', '255', '125/2183'))))

WITNESSES_K18 = tuple(sorted(WITNESSES+(
    (F(13, 100), Witness(F(341, 183), F(953, 3379),
                        F(138157, 7632), F(3659, 6837))),), key=lambda item: item[0]))


def profile(k_log, dimension):
    if k_log not in (16, 18, 20) or not 1 <= dimension < 128:
        raise ValueError('supported sizes are K16/K18/K20 and dimensions1..127')
    g = geometry(k_log, 128, dimension)
    return dict(g, k_log=k_log, requested_k=2**k_log,
                padded_k=g['k'], output_threshold=g['threshold'],
                nominal_rate=str(F(dimension, 128)),
                requested_rate=str(F(2**k_log, g['n'])))


def metadata(k_log, dimension):
    g = profile(k_log, dimension)
    return dict(construction=CONSTRUCTION, outer=[128, dimension], q=256,
                width=16, state=8, k_log=k_log, requested_k=g['requested_k'],
                padded_k=g['padded_k'], rows=g['rows'], n=g['n'],
                output_threshold=g['output_threshold'], distance='1/10',
                nominal_rate=g['nominal_rate'], requested_rate=g['requested_rate'])


def conditioning(rows, modal=True):
    if rows <= 0 or rows % 2:
        raise ValueError('positive even regional length required')
    return F(2**rows, math.comb(rows, rows//2)) if modal else F(rows+1)


def cell_upper(lo, hi, g, witness, prefactor):
    if not 0 < lo <= hi < g['n']:
        raise ValueError('invalid integer input-weight cell')
    source = row_upper(F(lo, g['n']), F(hi, g['n']), witness,
                       length=128, parity=128-g['dimension'])
    base = witness.x*witness.tilt
    endpoint = lo if base >= 1 else hi
    return (hi-lo+1)*prefactor*source**g['rows']/interval(base)**endpoint


def certify(k_log=20, dimension=88, first=F(3, 100), last=F(6, 25),
            initial_grid=1000, cell_bits=128, witnesses=None, modal=True):
    """Directed replay with unchanged two-state cover and source envelope.

    Modal conditioning is an exact uniform replacement of rows+1. The
    feasible-source truncation and rank-cap screen are deliberately not
    imported into this certificate. Existing replay witnesses remain valid
    at any outer dimension; only their numerical usefulness may change.
    """
    if witnesses is None:
        witnesses = WITNESSES_K18 if k_log == 18 else WITNESSES
    if initial_grid < 1 or cell_bits < 40 or not witnesses:
        raise ValueError('invalid directed replay parameters')
    if any(not isinstance(anchor, F) or not 0 < anchor < 1 or
           any(not isinstance(getattr(w, name), F) for name in ('x', 'y', 'tau', 'tilt'))
           for anchor, w in witnesses):
        raise ValueError('directed replay requires exact Fraction witnesses')
    g = profile(k_log, dimension)
    lo_weight, hi_weight = integer_span(g['n'], first, last)
    route = interval(conditioning(g['rows'], modal))**128
    prefactors = tuple(route*terminal_power(transfer(w), g['n']//16)
                       /(interval(w.tau)*interval(w.y)**g['threshold'])
                       for _, w in witnesses)

    def evaluate(lo, hi):
        middle = F(lo+hi, 2*g['n'])
        candidates = sorted(range(len(witnesses)),
                            key=lambda i: abs(witnesses[i][0]-middle))[:3]
        return min(((i, cell_upper(lo, hi, g,
                                 retarget_witness(*witnesses[i], middle), prefactors[i]))
                    for i in candidates), key=lambda item: item[1].hi)

    maximum = interval(F(1, 2**cell_bits)).lo
    span = (g['n']+initial_grid-1)//initial_grid
    total, next_weight, count, largest = ZERO, lo_weight, 0, ZERO.hi
    digest, uses = sha256(), [0]*len(witnesses)
    for lo, hi, index, value in refined_cells(lo_weight, hi_weight, span, evaluate, maximum):
        if lo != next_weight or hi < lo:
            raise ArithmeticError('middle-class coverage is not contiguous')
        next_weight = hi+1
        total += value
        largest = max(largest, value.hi)
        count += 1
        uses[index] += 1
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight != hi_weight+1:
        raise ArithmeticError('middle-class coverage ends early')
    margin = integer_margin(total.hi)
    if margin < 40:
        raise ArithmeticError(f'middle-class margin below40: {margin}')
    return dict(**metadata(k_log, dimension), complete=True, rigorous=True,
                whole_code=False, weights=[lo_weight, hi_weight],
                covered_weights=hi_weight-lo_weight+1,
                density_min=str(first), density_max_exclusive=str(last),
                failure_upper=str(total.hi), class_margin_bits_lower=margin,
                minimum_cell_margin_bits_lower=integer_margin(largest),
                checked_cells=count, witness_cell_counts=uses,
                partition_sha256=digest.hexdigest(), modal_conditioning=modal,
                precision_digits=80, scope='complete middle-input interval only')


def screen(k_log, dimension, densities, iterations=55, rank_history=False):
    """Yield bounded screens; each retains actual requested/padded geometry."""
    previous = None
    for density in densities:
        if rank_history:
            from gf8_rank_history_screen import search
            value = search(density, k_log, 128, dimension, 8, iterations=iterations,
                           start=previous, joint=True, cap_joint=True,
                           return_minima=(9, 9), feasible_source=True, modal=True)
        else:
            value = control_search(density, k_log, 128, dimension, 8,
                                   iterations=iterations, start=previous, joint=True)
            rows = profile(k_log, dimension)['rows']
            gain = 128*(math.log2(rows+1)-math.log2(conditioning(rows)))
            value['margin_bits'] += gain
            value['modal_conditioning'] = True
            value['conditioning_gain_bits'] = gain
        previous = value['log_parameters']
        value.update(metadata(k_log, dimension))
        value['complete'] = False
        yield value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--dimension', type=int, nargs='+', default=[88, 89, 90, 92, 96])
    parser.add_argument('--density', type=float, nargs='+',
                        default=[.03, .04, .06, .08, .10, .12, .13, .14, .16, .20, .24])
    parser.add_argument('--iterations', type=int, default=55)
    parser.add_argument('--rank-history', action='store_true')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--first', type=F, default=F(3, 100))
    parser.add_argument('--last', type=F, default=F(6, 25))
    parser.add_argument('--cell-bits', type=int, default=128)
    parser.add_argument('--legacy-conditioning', action='store_true')
    args = parser.parse_args()
    for dimension in args.dimension:
        if args.certify:
            print(json.dumps(certify(args.k_log, dimension, args.first, args.last,
                                     cell_bits=args.cell_bits,
                                     modal=not args.legacy_conditioning)), flush=True)
        else:
            for value in screen(args.k_log, dimension, args.density,
                                args.iterations, args.rank_history):
                print(json.dumps(value), flush=True)
