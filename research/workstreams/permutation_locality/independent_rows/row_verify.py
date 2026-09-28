"""Outward bounds for selected homogeneous row-interval classes.

Each selected four-row group has r active rows with weights in one interval.
All group locations and choices of row positions are counted. Different rows
may have different weights inside that interval. This is not a full cover.
Binary64 chooses witnesses only; every reported bound is replayed in Arb.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log
from pathlib import Path
import sys

import numpy as np
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
from scipy.optimize import minimize_scalar
from row_counts import row_gamma_exact, row_gamma_function

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bch_joint_support import authenticated_caps
from group_rank_one_verify import up
from mature_tail import TAIL_TERMINAL
from occupancy_memory import rounded
from occupancy_model import placement
from occupancy_screen import matrix_for_probabilities
from occupancy_sensitivity import float_placement
from two_group_screen import log_power_moment


def packet_law(p, active_rows):
    """Nonzero-packet probability and its conditional weight law."""
    assert 1 <= active_rows <= 4 and 0 < p <= 1
    active = 1-(1-p)**active_rows
    theta = [comb(active_rows, b)*p**b*(1-p)**(active_rows-b)/active
             if b <= active_rows else p*0 for b in range(1, 5)]
    return active, theta


def interval(value):
    parts = value.split(':')
    if len(parts) == 1:
        lo = hi = int(parts[0])
    elif len(parts) == 2:
        lo, hi = map(int, parts)
    else:
        raise argparse.ArgumentTypeError('use a weight or inclusive lo:hi')
    if not 1 <= lo <= hi <= 256:
        raise argparse.ArgumentTypeError('require 1 <= lo <= hi <= 256')
    return lo, hi


def float_inner(inner, q, tilt, p, active_rows, clamp=False):
    active, theta = packet_law(p, active_rows)
    local = inner.float_mix(theta)
    region = float_placement(local, q)
    matrix = matrix_for_probabilities(region, [active]*q)
    if clamp:
        from cone_moment import log_moment
        value = log_moment(matrix)
    else:
        value = log_power_moment(matrix, 256, TAIL_TERMINAL)
    return value+float(tilt)*209715


def replay(inner, caps, q, domain, active_rows, p, tilt, clamp=False):
    active, theta = packet_law(p, active_rows)
    local = inner.mix(theta)
    region = placement(local, rounding=rounded, maximum_groups=q)
    size = len(TAIL_TERMINAL)
    a = arb(active.numerator)/active.denominator
    matrix = sum((arb(comb(q, j))*a**j*(1-a)**(q-j)*region[j]
                  for j in range(q+1)), arb_mat(size, size))
    if clamp:
        from cone_moment import moment as clamped_moment
        moment = clamped_moment(rounded(matrix))
    else:
        matrix = rounded(matrix)**256
        moment = sum((matrix[0, j] for j, flag in enumerate(TAIL_TERMINAL) if flag), arb(0))
    gamma = row_gamma_exact(caps, *domain, p, exclude_zero=True)
    assert gamma > 0 and moment > 0
    scale = (arb(gamma.numerator)/gamma.denominator)**(q*active_rows)
    scale *= comb(2048, q)*comb(4, active_rows)**q
    return up(moment*scale*(arb(tilt)*209715).exp())


def self_test():
    # Thinning selected packets before placement equals averaging zero packets
    # inside each epoch. Noncommuting rational matrices test ordering too.
    operators = [fmpq_mat([[1, 1], [0, 1]]),
                 fmpq_mat([[1, 0], [1, 1]]),
                 fmpq_mat([[2, 1], [0, 1]]),
                 fmpq_mat([[1, 2], [1, 1]])]
    checks = 0
    for epochs in (1, 2, 3):
        for windows in (2, 3):
            region = placement(operators[:windows+1], epochs, windows,
                               fmpq_mat, lambda x: x, epochs*windows)
            for a in (Q(0), Q(1, 3), Q(2, 3), Q(1)):
                mixed = [sum((operators[k]*fmpq(str(comb(j, k)*a**k*(1-a)**(j-k)))
                              for k in range(j+1)), fmpq_mat(2, 2))
                         for j in range(windows+1)]
                reference = placement(mixed, epochs, windows,
                                      fmpq_mat, lambda x: x, epochs*windows)
                for q in range(epochs*windows+1):
                    expected = sum((region[k]*fmpq(str(comb(q, k)*a**k*(1-a)**(q-k)))
                                    for k in range(q+1)), fmpq_mat(2, 2))
                    assert expected == reference[q]
                    checks += 1
    for r in range(1, 5):
        for p in (Q(1, 5), Q(1, 2), Q(1)):
            a, theta = packet_law(p, r)
            assert sum(theta) == 1 and 0 < a <= 1
            assert sum((b+1)*theta[b] for b in range(4))*a == r*p
            checks += 1
    print('Exact thinning/placement and packet-law checks:', checks, 'passed', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=64)
    parser.add_argument('--occupancies', type=int, nargs='+', help='reuse local operators across group counts')
    parser.add_argument('--active-rows', type=int, nargs='+', choices=(1, 2, 3, 4), default=[4])
    parser.add_argument('--intervals', type=interval, nargs='+', default=[(w, w) for w in (64, 80, 96, 112, 128)])
    parser.add_argument('--tilt', default='.048')
    parser.add_argument('--tilts', nargs='+', help='reuse exact censuses across several positive tilts')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--cut', type=int, default=8)
    parser.add_argument('--full-feedback', type=int, default=6)
    parser.add_argument('--window-histogram', type=int, default=8)
    parser.add_argument('--joint-cancellation', action='store_true')
    parser.add_argument('--averaged-high', action='store_true', help='average reference packet weights in high-occupancy output moments')
    parser.add_argument('--class-tail', action='store_true', help='use coupled class-specific mature-tail arrivals')
    parser.add_argument('--column-density', action='store_true', help='bound multi-packet density columns by conditioned one-window averages')
    parser.add_argument('--feedback-density', type=int, choices=range(0,7), default=0,
                        help='exact feedback/expansion-class convolution through this local occupancy')
    parser.add_argument('--clamp', action='store_true', help='enforce true mature-subset inequalities after each region')
    parser.add_argument('--max-evaluations', type=int, default=35)
    parser.add_argument('--p', type=Q, help='explicit rational witness instead of optimization')
    parser.add_argument('--screen-only', action='store_true')
    parser.add_argument('--replay-below', type=float, help='also replay screened classes at or below this log2 proposal')
    parser.add_argument('--self-test-only', action='store_true')
    parser.add_argument('--target-bits', type=int, default=52)
    args = parser.parse_args()
    if any(not 1 <= q <= 2048 for q in (args.occupancies or [args.groups])) or any(Q(x) <= 0 for x in (args.tilts or [args.tilt])):
        parser.error('require groups in [1,2048] and positive tilt')
    if args.p is not None and not 0 < args.p <= 1:
        parser.error('p must be in (0,1]')
    self_test()
    if args.self_test_only:
        return
    import shape_inner
    args.penalty = '1'
    shared = shape_inner.prepare_shared(args)
    for tilt in args.tilts or [args.tilt]:
        args.tilt = tilt
        inner = shape_inner.build(args, shared)
        if args.averaged_high:
            from averaged_windows import AveragedHighInner
            inner = AveragedHighInner(inner)
            print('Using outward categorical window moments above the shape cutoff', flush=True)
        for q in args.occupancies or [args.groups]:
            args.groups = q
            run_classes(args, inner)


def run_classes(args, inner):
    ctx.prec = args.precision
    caps = authenticated_caps()
    q = args.groups
    cache = {}
    def moment(p, r):
        key = p, r
        if key not in cache:
            cache[key] = float_inner(inner, q, args.tilt, p, r, args.clamp)
        return cache[key]
    for r in args.active_rows:
        for domain in args.intervals:
            if not any(caps[w] for w in range(domain[0], domain[1]+1)):
                print('Excluded interval:', domain, flush=True)
                continue
            constant = log(comb(2048, q))+q*log(comb(4, r))
            log_gamma = row_gamma_function(caps, *domain, exclude_zero=True)
            def score(p):
                return (moment(p, r)+constant+q*r*log_gamma(p))/log(2)
            if args.p is not None:
                p = args.p
            elif domain == (256, 256):
                p = Q(1)
            else:
                fit = minimize_scalar(lambda x: score(1/(1+np.exp(-x))),
                                      bounds=(-7, 7), method='bounded',
                                      options={'maxiter': args.max_evaluations, 'xatol': 1e-5})
                witness = 1/(1+np.exp(-fit.x))
                p = Q(max(1, min(10**9-1, round(witness*10**9))), 10**9)
            proposal = score(float(p))
            print('BINARY64 proposal q/r/interval', q, r, domain, 'p', str(p),
                  'log2 upper', proposal, flush=True)
            if not args.screen_only or (args.replay_below is not None and proposal <= args.replay_below):
                result = replay(inner, caps, q, domain, r, p, args.tilt, args.clamp)
                print('OUTWARD ROW-INTERVAL q/r/interval', q, r, domain, 'p', str(p),
                      'tilt', args.tilt, 'precision', ctx.prec, 'cut', args.cut,
                      'refinements', {'class_tail':args.class_tail,'averaged_high':args.averaged_high,'clamp':args.clamp,
                                      'column_density':args.column_density,'feedback_density':args.feedback_density},
                      'upper', result, 'log2 upper', result.log()/arb(2).log(),
                      'meets class budget', bool(result < arb(2)**(-args.target_bits)), flush=True)
    print('Selected homogeneous interval classes only; unequal weights inside each interval included.', flush=True)
    print('Not an occupancy-wide cover or full-code certificate.', flush=True)


if __name__ == '__main__':
    main()
