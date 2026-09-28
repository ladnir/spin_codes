"""Outward replay for selected mixtures of two row-interval group types.

Type i has r_i positive rows whose weights lie in I_i, and 4-r_i zero
rows. Its witness probability is fixed across all256 regions. Counts fix
the number of groups of each type; all group and row locations are counted.
This driver is not a complete type or occupancy cover.
"""
import argparse
from fractions import Fraction as Q
from math import comb

from flint import arb, ctx

from row_verify import interval, packet_law
from row_counts import row_gamma_exact
import shape_inner
from typed_inner import build_typed_operators
from typed_placement import typed_placement
from bch_joint_support import authenticated_caps
from group_rank_one_verify import up


def packet_weights(p, rows):
    active, conditional = packet_law(p, rows)
    return [1-active]+[active*value for value in conditional]


def replay(model, caps, counts, rows, domains, probabilities, tilt):
    q = sum(counts)
    laws = [packet_weights(p, r) for p, r in zip(probabilities, rows)]
    operators = build_typed_operators(model, *laws, max_counts=counts)
    print('Typed epoch operators:', len(operators), 'counts', counts, flush=True)
    region = typed_placement(operators, counts)
    result = region**256
    moment = sum((result[0, j] for j, flag in enumerate(model.terminal) if flag), arb(0))
    scalar = arb(comb(2048, q)*comb(q, counts[0]))
    for count, r, domain, p in zip(counts, rows, domains, probabilities):
        if count:
            gamma = row_gamma_exact(caps, *domain, p, exclude_zero=True)
            scalar *= comb(4, r)**count*(arb(gamma.numerator)/gamma.denominator)**(r*count)
    return up(moment*scalar*(arb(tilt)*209715).exp())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=64)
    parser.add_argument('--second-counts', type=int, nargs='+', default=[1])
    parser.add_argument('--rows', type=int, nargs=2, choices=(1, 2, 3, 4), default=[4, 3])
    parser.add_argument('--intervals', type=interval, nargs=2, default=[(64,192), (64,192)])
    parser.add_argument('--probabilities', type=Q, nargs=2, default=[Q(1,2), Q(1,2)])
    parser.add_argument('--tilt', default='.052')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--cut', type=int, default=8)
    parser.add_argument('--full-feedback', type=int, default=6)
    parser.add_argument('--window-histogram', type=int, default=8)
    parser.add_argument('--joint-cancellation', action='store_true')
    parser.add_argument('--target-bits', type=int, default=52)
    args = parser.parse_args()
    if (not 1 <= args.groups <= 2048
            or any(not 0 <= n <= args.groups for n in args.second_counts)
            or any(not 0 < p <= 1 for p in args.probabilities)
            or Q(args.tilt) <= 0):
        parser.error('invalid group counts, probabilities, or tilt')
    model = shape_inner.build(args)
    caps = authenticated_caps()
    ctx.prec = args.precision
    for second in args.second_counts:
        counts = args.groups-second, second
        upper = replay(model, caps, counts, args.rows, args.intervals, args.probabilities, args.tilt)
        print('OUTWARD TWO-TYPE counts', counts, 'rows', args.rows, 'intervals', args.intervals,
              'probabilities', [str(p) for p in args.probabilities], 'tilt', args.tilt,
              'upper', upper, 'log2 upper', upper.log()/arb(2).log(),
              'meets class budget', bool(upper < arb(2)**(-args.target_bits)), flush=True)
    print('Only the stated type-count classes; no full-occupancy or full-code certificate.', flush=True)


if __name__ == '__main__':
    main()
