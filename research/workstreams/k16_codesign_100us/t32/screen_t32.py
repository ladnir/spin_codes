"""Fresh partial screens for physical t32/s8 inners and the K16 RS16 outer.

The expansion is RM1(5) plus x0*x1 and one other quadratic monomial.
Feedback is its transpose. Every32-bit physical step has an independently
sampled uniform GL8 update. State starts at zero, persists across every
step and region, and has no final flush. No old numerical endpoint enters.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb, log
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'k16_design'))

import numpy as np
from flint import arb, ctx
from scipy.special import logsumexp
import packet_q1 as q1
import packet_uniform_tail as tail
import packet_regional_log as logarithmic

aq, up = q1.kernel_t64.aq, q1.kernel_t64.up
PAIRS = tuple(pair for pair in combinations(range(5), 2) if pair != (0, 1))


def prepare(pair):
    if pair != 'rank4' and tuple(pair) not in PAIRS:
        raise ValueError('a second quadratic pair other than(0,1) required')
    rows = [2**32-1]
    rows.extend(sum(1 << x for x in range(32) if (x >> i) & 1) for i in range(5))
    forms = (((0, 1), (2, 3)), ((0, 2), (1, 4))) if pair == 'rank4' else (((0, 1),), (tuple(pair),))
    for terms in forms:
        row = 0
        for i, j in terms:
            row ^= sum(1 << x for x in range(32) if ((x >> i) & (x >> j)) & 1)
        rows.append(row)
    columns = [sum(((row >> x) & 1) << i for i, row in enumerate(rows)) for x in range(32)]
    rank = q1.kernel_t64.s16_maps.binary_rank
    if rank(rows) != 8 or any((a & b).bit_count() & 1 for a in rows for b in rows):
        raise ArithmeticError('rank or CA=0 failure')
    if [rank(columns[i:i+4]) for i in range(0, 32, 4)] != [4]*8:
        raise ArithmeticError('single-packet feedback is not injective')
    images = tuple(q1.kernel_t64.s16_maps.images_from_rows(rows))
    spectrum = Counter(x.bit_count() for x in images)
    data = q1.kernel_t64.kernel_maps.prepare_maps(images, columns, bits=8,
               distribution='uniform_gl', birth_density='capped')
    record = dict(physical_t=32, state_bits=8, packet_bits=4, physical_windows=8,
                  quadratic_forms=forms, expansion_rows_hex=list(map(hex, rows)),
                  feedback_columns=columns, expansion_spectrum=dict(sorted(spectrum.items())),
                  expansion_rank=8, packet_ranks=[4]*8, feedback_times_expansion_zero=True,
                  feedback_definition='C=A^T', physical_steps=4096,
                  physical_steps_per_region=64, regions=64, group_count=512,
                  independent_uniform_gl_per_physical_step=True,
                  zero_initial_state=True, final_flush=False,
                  state_continuity='retained between every physical step and region',
                  whole_code_certificate=False)
    return data, record


def local_operators(data, tilt):
    q1.kernel_t64.kernel_maps.authenticate(data)
    if data['windows'] != 8 or data['bits'] != 8:
        raise ValueError('explicit physical t32/s8 maps required')
    z = (-aq(Q(tilt))).exp()
    local = q1.kernel_t64.sparse_kernel.outward_at_z(data, z)
    return q1.kernel_t64.kernel_birth_density.refine_local(data, local, z, Q(1, 2))


def macro_operators(data, tilt):
    """Four chronological physical steps, with exact hypergeometric splits.

    First compose two8-window steps, then compose the two16-window halves.
    The second composition retains the intermediate state. Every uniform
    32-window subset has exactly its hypergeometric weight in both splits.
    No t64 wrapper or t128 physical map is used.
    """
    physical = local_operators(data, tilt)
    two = q1.kernel_t64.convolve(physical)
    return q1.kernel_t64.convolve(two)


def screen_q1(pair, tilts):
    ctx.prec = 256
    data, record = prepare(pair)
    _, counts = tail.uniform_envelope(16, 8, 4)
    best = [arb(1)]*65
    for tilt in tilts:
        # Direct physical-step placement avoids any assumption about macros.
        local = local_operators(data, tilt)
        regional = q1.placement(local, epochs=64, windows=8,
                      maximum_groups=1, rounding=q1.rounded)
        moments = q1.support_moments(regional[0], regional[1], 64)
        factor = (aq(Q(tilt))*13107).exp()
        best = [min(old, up(factor*moment)) for old, moment in zip(best, moments)]
    terms = [up(512*aq(count)*value) for count, value in zip(counts, best)]
    upper = up(sum(terms, arb(0)))
    margin = -upper.log()/arb(2).log()
    print(f'OUTWARD q1 pair={pair} margin={margin} '
          f'dominant={sorted(range(65),key=lambda v:float(terms[v]),reverse=True)[:5]} '
          f'spectrum={record["expansion_spectrum"]}; whole_code_certificate=false', flush=True)
    return data, upper


def screen_tail(pair, tilts, q_min, q_max):
    ctx.prec = 256
    data, _ = prepare(pair)
    beta, _ = tail.uniform_envelope(16, 8, 4)
    beta_log = float(aq(beta).log())
    best = np.full(q_max+1, -np.inf)
    choices = {}
    for tilt in tilts:
        local = macro_operators(data, tilt)
        regional = logarithmic.log_placement(logarithmic.upper_logs(local), q_max,
                                             epochs=16, windows=32)
        for q in range(q_min, q_max+1):
            mixed = logarithmic.regional_uniform_log(regional, q)
            moment = logarithmic.log_power_matrix(mixed, 64)
            margin = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+moment)/log(2)
            if margin > best[q]:
                best[q], choices[q] = margin, tilt
        worst = sorted(range(q_min, q_max+1), key=lambda q: best[q])[:8]
        print(f'PROPOSAL pair={pair} tilt={tilt} worst={[(q,float(best[q]),choices[q]) for q in worst]}', flush=True)
    print(f'PROPOSAL pair={pair} union={-logsumexp(-best[q_min:]*log(2))/log(2)}; whole_code_certificate=false', flush=True)
    return best, choices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair', nargs=2, type=int)
    parser.add_argument('--rank4', action='store_true')
    parser.add_argument('--tilts', nargs='+', default=['.00128', '.00256', '.00512', '.01024', '.0256', '.0512'])
    parser.add_argument('--tail', action='store_true')
    parser.add_argument('--q-min', type=int, default=3)
    parser.add_argument('--q-max', type=int, default=128)
    args = parser.parse_args()
    pairs = ['rank4'] if args.rank4 else [tuple(args.pair)] if args.pair else PAIRS
    for pair in pairs:
        if args.tail:
            screen_tail(pair, args.tilts, args.q_min, args.q_max)
        else:
            screen_q1(pair, args.tilts)


if __name__ == '__main__':
    main()
