"""Outward selected-family replay with persistent group types.

One family is the all-middle row category; the other is selected explicitly.
Every original group type inside each family is covered by an exact-input
Holder envelope. This is not a cover of all family-count assignments.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb, factorial, prod

from flint import arb, ctx

from holder_exact import holder_envelope
from holder_pilot import INTERVALS, families
from row_counts import row_gamma_exact
import shape_inner
from typed_inner import build_typed_operators
from typed_placement import typed_placement
from bch_joint_support import authenticated_caps
from group_rank_one_verify import up


def exact_types(caps, p):
    """Construct all69 nonzero four-row types from exact witness probabilities."""
    if not isinstance(p, Q) or not 0 < p < Q(1,2):
        raise ValueError('an exact rational witness in (0,1/2) is required')
    covered = [w for lo, hi in INTERVALS for w in range(lo, hi+1)]
    assert len(covered) == len(set(covered))
    assert all(w in covered for w, count in enumerate(caps) if count)
    probabilities = (Q(0), p, Q(1,2), 1-p, Q(1))
    gamma = [row_gamma_exact(caps, lo, hi, x)
             for (lo,hi), x in zip(INTERVALS, probabilities)]
    result = []
    for labels in combinations_with_replacement(range(5), 4):
        if labels == (0,0,0,0):
            continue
        multiplicity = factorial(4)//prod(factorial(n) for n in Counter(labels).values())
        h = Q(multiplicity)*prod(gamma[label] for label in labels)
        law = [Q(1)]
        for label in labels:
            x = probabilities[label]
            updated = [Q(0)]*(len(law)+1)
            for b, mass in enumerate(law):
                updated[b] += mass*(1-x)
                updated[b+1] += mass*x
            law = updated
        assert sum(law) == 1
        result.append((labels, h, tuple(law)))
    assert len(result) == 69
    return result


def replay(model, envelopes, counts, tilt):
    (z1, law1), (z2, law2) = envelopes
    assert z1 > 0 and z2 > 0
    local = build_typed_operators(model, law1, law2, max_counts=counts)
    region = typed_placement(local, counts)
    matrix = region**256
    moment = sum((matrix[0,j] for j, flag in enumerate(model.terminal) if flag), arb(0))
    q = sum(counts)
    scalar = arb(comb(2048,q)*comb(q,counts[1]))
    for z, n in zip((z1,z2), counts):
        scalar *= (arb(z.numerator)/z.denominator)**(256*n)
    return up(scalar*moment*(arb(tilt)*209715).exp())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=64)
    parser.add_argument('--exception-counts', type=int, nargs='+', default=[1])
    parser.add_argument('--family', default='exactly_1_low_strict')
    parser.add_argument('--p', type=Q, default=Q(1,4))
    parser.add_argument('--tilt', default='.052')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--cut', type=int, default=8)
    parser.add_argument('--full-feedback', type=int, default=6)
    parser.add_argument('--window-histogram', type=int, default=8)
    parser.add_argument('--joint-cancellation', action='store_true')
    parser.add_argument('--target-bits', type=int, default=52)
    parser.add_argument('--envelope-only', action='store_true')
    args = parser.parse_args()
    if (not 1 <= args.groups <= 2048
            or any(not 0 <= n <= args.groups for n in args.exception_counts)
            or not 0 < args.p < Q(1,2) or Q(args.tilt) <= 0):
        parser.error('invalid occupancy, rational witness, or tilt')
    ctx.prec = args.precision
    caps = authenticated_caps()
    family = families(exact_types(caps, args.p))
    if args.family not in family or args.family == 'middle':
        parser.error('choose a listed exceptional family from holder_pilot.py')
    selected = family['middle'], family[args.family]
    if set(row[0] for row in selected[0]) & set(row[0] for row in selected[1]):
        parser.error('exception family must be disjoint from the middle family')
    envelopes = []
    for name, rows in zip(('middle',args.family), selected):
        z, law = holder_envelope(((h,theta) for _,h,theta in rows), precision=args.precision)
        envelopes.append((z,law))
        z_arb = arb(z.numerator)/z.denominator
        print('OUTWARD FAMILY', name, 'types', len(rows), 'witness', str(args.p),
              'scalar log2 cost', 256*z_arb.log()/arb(2).log(),
              'exact rational packet law established', flush=True)
    if args.envelope_only:
        print('Scalar product-envelope costs, not distance margins.', flush=True)
        return
    model = shape_inner.build(args)
    ctx.prec = args.precision
    for exceptions in args.exception_counts:
        counts = args.groups-exceptions, exceptions
        upper = replay(model, envelopes, counts, args.tilt)
        print('OUTWARD FAMILY MIXTURE', args.family, 'counts', counts,
              'p', str(args.p), 'tilt', args.tilt, 'upper', upper,
              'log2 upper', upper.log()/arb(2).log(),
              'meets class budget', bool(upper < arb(2)**(-args.target_bits)), flush=True)
    print('Selected family-count classes only; not a complete occupancy or code certificate.', flush=True)


if __name__ == '__main__':
    main()
