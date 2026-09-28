"""Measure remaining transition losses on the current unrestricted envelope.

Deleted entries are counterfactuals, NOT upper bounds. The two mass-based
choices use proved local inequalities, but this binary64 selected-point
screen is not an outward certificate or complete support cover.
"""
import argparse
from fractions import Fraction as Q

from mass_density_screen import epochs, as_array, score, baseline, float_placement
from mass_density import blend, conditioned_coefficients
from occupancy_memory import Z, F, M, C, U


def choices(base, census, tilt, penalty):
    yield 'baseline', base
    candidate = conditioned_coefficients(census, 6, 12, tilt, penalty)
    fractions = {j: Q(1) for j in range(6, 13)}
    zero = blend(base, candidate, fractions, target=Z)
    yield 'valid-local/mass-zero', zero
    # Apply zero first: blend checks that its input density column is unsplit.
    yield 'valid-local/mass-zero-and-density', blend(zero, candidate, fractions, target=C)
    for name, sources, targets in (
        ('fresh-return', (F,), (Z,)),
        ('fresh-density', (F,), (C,)),
        ('fresh-return-and-density', (F,), (Z, C)),
        ('mature-return', (C,), (Z,)),
        ('mature-density', (C,), (C,)),
        ('uniform-return', tuple(range(U, U+5)), (Z,)),
    ):
        changed = [t*1 for t in base]
        for matrix in changed[1:]:
            for source in sources:
                for target in targets:
                    matrix[source, target] = 0
        yield 'INVALID-ABLATION/'+name, changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt', default='.056')
    parser.add_argument('--penalty', default='.9')
    parser.add_argument('--groups', type=int, nargs='+', default=[64, 80])
    parser.add_argument('--supports', type=int, nargs='+', default=[192, 200])
    args = parser.parse_args()
    if (Q(args.tilt) <= 0 or not 0 < Q(args.penalty) <= 1
            or any(not 1 <= q <= 2048 for q in args.groups)
            or any(not 38 <= u <= 256 for u in args.supports)):
        parser.error('positive tilt, penalty in (0,1], q in 1..2048 and support in 38..256 required')
    print('DIAGNOSTIC ONLY: deleted transitions are not valid probability bounds', flush=True)
    base, census = epochs(args.tilt, args.penalty)
    caps = baseline.authenticated_caps()
    cdf = baseline.integer_cdf(baseline.weighted_cdf_upper(caps, 1 << 128, full_weight=1/Q(args.penalty)))
    shells = baseline.weighted_union_shells(caps, full_weight=1/Q(args.penalty))
    counts = {u: min(Q(cdf[u]), shells[u]) for u in args.supports}
    for name, matrices in choices(base, census, args.tilt, args.penalty):
        region = float_placement(as_array(matrices), max(args.groups))
        for q in args.groups:
            for u in args.supports:
                value, p = score(region, q, u, counts[u], args.tilt)
                print(name, 'q/u', q, u, 'tilt/rho', args.tilt, args.penalty,
                      'log2 score', value, 'p', p, flush=True)


if __name__ == '__main__':
    main()
