"""Outward common-potential replay for one homogeneous support event.

This is NOT a cover of all supports at an occupancy, nor a full-code cert.
Local bounds are rebuilt or read from the current source-bound attack memo.
No historical local-family input is accepted here.
"""
import argparse
from fractions import Fraction as Q
from math import comb

from flint import arb, arb_mat, ctx

import attack_cache
import mass_density_screen as screen
from mixing_attack import retarget
from shape_potential import matrix_potential
from group_rank_one_verify import up


def point(x):
    x = Q(x)
    return arb(x.numerator) / x.denominator


def region_action(operators, potential, degree, *, epochs=64, windows=32):
    """Exact-coefficient polynomial recurrence, outward at each epoch."""
    n = len(potential)
    if (any(type(x) is not int for x in (degree, epochs, windows))
            or min(epochs, windows) < 1 or not 0 <= degree <= epochs*windows
            or len(operators) <= min(degree, windows)
            or any(t.nrows() != n or t.ncols() != n for t in operators)
            or any(not t[i,j].is_finite() or not t[i,j] >= 0
                   for t in operators for i in range(n) for j in range(n))
            or any(not v.is_finite() or not v > 0 for v in potential)):
        raise ValueError('complete operators, positive potential and valid geometry required')
    coefficients = [arb_mat([[v] for v in potential])]
    for epoch in range(1, epochs+1):
        following = []
        for r in range(min(degree, epoch*windows)+1):
            column = arb_mat(n, 1)
            for k in range(max(0, r-len(coefficients)+1), min(r, windows)+1):
                column += operators[k] * coefficients[r-k] * comb(windows, k)
            following.append(arb_mat([[up(column[i, 0])] for i in range(n)]))
        coefficients = following
    return [[up(column[i, 0] / comb(epochs*windows, r)) for i in range(n)]
            for r, column in enumerate(coefficients)]


def replay(operators, potential, p, q, support, count, tilt, cutoff):
    """Certified finite bound for the explicitly named homogeneous event."""
    potential = [point(x) for x in potential]
    if potential[0] != 1:
        raise ValueError('normalize the zero-state potential to exactly one')
    p = point(p)
    if not 0 < p < 1:
        raise ValueError('interior reference probability required')
    action = region_action(operators, potential, q)
    averaged = [arb(0) for _ in potential]
    for j, column in enumerate(action):
        probability = comb(q, j)*p**j*(1-p)**(q-j)
        for i, value in enumerate(column):
            averaged[i] += probability*value
    lam = max(up(x/v) for x, v in zip(averaged, potential))
    constant = max(up(arb(int(t))/v) for t, v in zip(screen.baseline.TAIL_TERMINAL, potential))
    mass = comb(256, support)*p**support*(1-p)**(256-support)
    if not mass > 0:
        raise ArithmeticError('conditioning probability must be positive')
    upper = up(constant * lam**256 * (point(tilt)*cutoff).exp()
               * comb(2048, q) * (point(count)/mass)**q)
    return upper, lam, constant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt', default='.072')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--groups', type=int, default=96)
    parser.add_argument('--support', type=int, default=200)
    parser.add_argument('--cutoff', type=int, default=193986)
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--local-precision', type=int, default=192)
    parser.add_argument('--directory', default='tmp/four-bit-attack-cache')
    args = parser.parse_args()
    if (args.local_precision < 128 or args.precision < args.local_precision or not 2 <= args.rounds <= 32
            or not 1 <= args.groups <= 2048 or not 38 <= args.support < 256
            or not 0 <= args.cutoff < 1 << 21 or Q(args.tilt) <= 0):
        parser.error('invalid parameters')
    ctx.prec = args.precision
    base, _, details = attack_cache.build([args.tilt], precision=args.local_precision,
                                         directory=args.directory)[args.tilt]
    ctx.prec = args.precision
    operators = retarget(base, details['spectrum'], args.tilt, 2, args.rounds)
    caps = screen.baseline.authenticated_caps()
    cdf = screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(
        caps, 1 << 128, full_weight=Q(10, 9)))
    shells = screen.baseline.weighted_union_shells(caps, full_weight=Q(10, 9))
    count = min(Q(cdf[args.support]), shells[args.support])
    # Legacy BCH verifier imports can reset the process-wide Flint context.
    # Restore the requested global precision AFTER loading those helpers.
    ctx.prec = args.precision
    region = screen.float_placement(screen.as_array(operators), args.groups)
    proposal, probability = screen.score(region, args.groups, args.support, count,
                                        args.tilt, cutoff=args.cutoff)
    p = Q(round(probability*10**12), 10**12)
    matrix = screen.baseline.matrix_for_probabilities(region, [float(p)]*args.groups)
    v = matrix_potential(matrix)
    # Decimal witnesses are exact rationals; binary64 affects only the search.
    potential = [Q(str(float(x))) for x in v]
    potential[0] = Q(1)
    assert ctx.prec == args.precision
    print('SELECTED REPLAY q/u/R', args.groups, args.support, args.rounds,
          'tilt', args.tilt, 'cutoff', args.cutoff, 'precision', ctx.prec,
          'local precision', args.local_precision,
          'binary64 proposal', proposal, 'p', p, 'v', list(map(str, potential)), flush=True)
    upper, lam, constant = replay(operators, potential, p, args.groups, args.support,
                                  count, args.tilt, args.cutoff)
    assert ctx.prec == args.precision
    print('OUTWARD SELECTED-EVENT upper', upper, 'log2 upper', upper.log()/arb(2).log(),
          'lambda', lam, 'terminal factor', constant, 'below 2^-40', upper < arb(2)**-40,
          flush=True)
    print('Only the stated homogeneous union-support event; NOT an occupancy cover or full-code certificate.')


if __name__ == '__main__':
    main()
