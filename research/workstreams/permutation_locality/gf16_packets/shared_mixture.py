"""Exact positive-mixture majorant for the SHARED coordinate-shuffle ensemble.

This supplies a new outer interface, not an independent-row certificate.
An upper CDF B(u) is an upper bound on the shell count S(u) separately;
differences B(u)-B(u-1) are never treated as shell bounds.
The optional Fourier-monotone comparison instead constructs a different
measure that dominates the whole linear output-weight moment; its shell
masses are not bounds on the original code's shell counts.
"""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

from flint import ctx
import scalar_cover
import birth_classes
from shared_support import shared_counts

SCHEMA = 'shared-gf16-positive-mixture-1'


def envelope(caps, centers, zero_bits=None, cost_tilt=Q(1), mass_bits=None, allow_empty=False):
    """Dominate every nonzero shell by an exactly checked Bernoulli mixture.

    Assign each shell to a center minimizing its required tilted mass.
    With cost_tilt=1 this maximizes its binomial probability. The optional
    empty-input constraint excludes centers with a large artificial atom.
    A separate optional cap limits the total mass of each component.
    The coefficient of that component is the maximum required mass over
    its assigned shells. Maximization and verification are rational.
    Other positive components only strengthen the checked majorant.
    """
    n = len(caps) - 1
    centers = tuple(map(Q, centers))
    cost_tilt = Q(cost_tilt)
    if (n < 1 or type(allow_empty) is not bool or (caps[0] != 0 and not allow_empty)
            or any(type(v) not in (int, Q) or v < 0 for v in caps)
            or not centers or any(not 0 < p <= 1 for p in centers)
            or len(set(centers)) != len(centers)
            or not 0 < cost_tilt <= 1
            or (mass_bits is not None and (type(mass_bits) is not int or not 0 <= mass_bits <= 2048))
            or (zero_bits is not None and (type(zero_bits) is not int or not 0 <= zero_bits <= 512))):
        raise ValueError('nonnegative exact shell caps and distinct positive centers required')
    masses = [Q(0)] * len(centers)
    kernels = [[comb(n, u) * p**u * (1-p)**(n-u) for u in range(n+1)] for p in centers]
    costs = [(1-p+p*cost_tilt)**n for p in centers]
    for u in range(n+1):
        if not caps[u]:
            continue
        eligible = [j for j, p in enumerate(centers) if kernels[j][u] and
                    (mass_bits is None or caps[u] <= (1 << mass_bits)*kernels[j][u]) and
                    (zero_bits is None or Q(caps[u])*(1-p)**n <= Q(2)**-zero_bits*kernels[j][u])]
        if not eligible:
            raise ValueError('centers cannot cover a nonzero shell')
        j = max(eligible, key=lambda j: kernels[j][u]/costs[j])
        masses[j] = max(masses[j], Q(caps[u]) / kernels[j][u])
    result = [(c, p) for c, p in zip(masses, centers) if c]
    verify(caps, result, allow_empty)
    if zero_bits is not None:
        assert all(c*(1-p)**n <= Q(2)**-zero_bits for c, p in result)
    if mass_bits is not None:
        assert all(c <= 1 << mass_bits for c, _ in result)
    return result


def verify(caps, components, allow_empty=False):
    n = len(caps) - 1
    if (n < 1 or type(allow_empty) is not bool or (caps[0] != 0 and not allow_empty)
            or any(type(v) not in (int, Q) or v < 0 for v in caps)
            or any(Q(c) <= 0 or not 0 < Q(p) <= 1 for c, p in components)):
        raise ValueError('valid positive components and nonzero-message shell caps required')
    for u in range(n+1):
        upper = comb(n, u) * sum((Q(c)*Q(p)**u*(1-Q(p))**(n-u) for c, p in components), Q(0))
        if upper < caps[u]:
            raise ArithmeticError(f'shell {u} is not dominated')


def as_components(mixture):
    # Keep the real zero message separate. Positive comparison components
    # retain the active label even when their Bernoulli draw is all zero.
    return [('zero', Q(1), (Q(1), Q(0), Q(0), Q(0), Q(0)), 0)] + [
        (f'shared-{i}', c, tuple(scalar_cover.probabilities(p)), 1)
        for i, (c, p) in enumerate(mixture)]


def actual_components(step=8, refined=True, coupled=False, zero_bits=None, cost_tilt=Q(1),
                      mass_bits=None, joint=False, wide=False, monotone=False):
    if type(step) is not int or not 1 <= step <= 64:
        raise ValueError('integer support-grid step in 1..64 required')
    cdf, _ = shared_counts(refined, coupled, joint, wide)
    if type(monotone) is not bool:
        raise ValueError('boolean monotone-comparison flag required')
    caps = cdf
    if monotone:
        from monotone_comparison import thinned_shells
        caps = thinned_shells(cdf)
    support = max(1, next(u for u, value in enumerate(caps) if value))
    centers = sorted({Q(u, 256) for u in range(support, 257, step)} | {Q(1)})
    mixture = envelope(caps, centers, zero_bits, cost_tilt, mass_bits, allow_empty=monotone)
    return as_components(mixture), caps, mixture


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=int, default=8)
    parser.add_argument('--updates', type=int, choices=(2, 3, 4), default=4)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--distance', default='.1')
    parser.add_argument('--minimum-groups', type=int, default=17)
    parser.add_argument('--base-tilt', default='1/32')
    parser.add_argument('--coupled-counts', action='store_true')
    parser.add_argument('--joint-counts', action='store_true')
    parser.add_argument('--wide-counts', action='store_true')
    parser.add_argument('--monotone-comparison', action='store_true', help='Use Fourier-monotone CDF comparison; not actual shell bounds')
    parser.add_argument('--zero-bits', type=int, help='Cap empty input mass of each active comparison component by 2^-bits')
    parser.add_argument('--cost-tilt', default='1', help='Optimize tilted mixture mass rather than raw mass; shell checks remain exact')
    parser.add_argument('--mass-bits', type=int, help='Cap the total mass of each comparison component by 2^bits')
    parser.add_argument('--variance-bins', type=int, default=0)
    parser.add_argument('--regional-count', action='store_true')
    parser.add_argument('--probe', nargs='*', default=[], help='Selected mean fractions; diagnostics, not a cover')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.precision < 128 or not 0 < Q(args.distance) < 1:
        parser.error('precision >=128 and distance in (0,1) required')
    if any(not 0 <= Q(x) <= 1 for x in args.probe):
        parser.error('probe means must be in [0,1]')
    if args.variance_bins < 0 or args.variance_bins > 64 or (args.regional_count and not args.variance_bins):
        parser.error('regional counts require positive variance bins (at most 64)')
    if (args.joint_counts and not args.coupled_counts) or (args.wide_counts and not args.joint_counts):
        parser.error('joint counts require coupled counts; wide counts require joint counts')
    components, cdf, mixture = actual_components(args.step, coupled=args.coupled_counts,
        zero_bits=args.zero_bits, cost_tilt=Q(args.cost_tilt), mass_bits=args.mass_bits,
        joint=args.joint_counts, wide=args.wide_counts, monotone=args.monotone_comparison)
    print('EXACT COMPARISON SHELL MAJORANT checked all 257 shells including zero;'
          if args.monotone_comparison else 'EXACT SHARED SHELL MAJORANT checked all 256 nonzero shells;',
          len(mixture), 'positive components', flush=True)
    for c, p in mixture:
        print('  activity', p, 'log2 mass', scalar_cover.logq(c)/scalar_cover.np.log(2), flush=True)
    empty = sum((c*(1-p)**256 for c, p in mixture), Q(0))
    print('ACTIVE COMPARISON empty-input log2 mass',
          scalar_cover.logq(empty)/scalar_cover.np.log(2) if empty else '-inf', flush=True)
    record = dict(schema=SCHEMA, step=args.step, refined_counts=True, coupled_counts=args.coupled_counts,
                  joint_counts=args.joint_counts, wide_counts=args.wide_counts,
                  monotone_comparison=args.monotone_comparison,
                  zero_bits=args.zero_bits, cost_tilt=args.cost_tilt, mass_bits=args.mass_bits,
                  empty_active_mass=str(empty),
                  shell_cap_source=('Fourier-monotone thinned CDF comparison, not actual shell bounds'
                      if args.monotone_comparison else 'shared union-support CDF itself, not its differences'),
                  mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
                  shell_caps=list(map(str, cdf)), proof_status='shell majorant only; not a distance certificate')
    if args.probe:
        data = birth_classes.actual(args.updates)
        ctx.prec = args.precision
        model = scalar_cover.Model(components, data, int(Q(args.distance)*scalar_cover.N),
                                   args.minimum_groups, Q(args.base_tilt), inner=birth_classes,
                                   variance_shuffle=bool(args.variance_bins), variance_bins=args.variance_bins,
                                   regional_count=args.regional_count)
        print('SHARED DIAGNOSTIC ROOT', tuple(map(str, model.root)), flush=True)
        rows = []
        for value in map(Q, args.probe):
            if model.empty((value, value)):
                print('Empty mean', value, flush=True)
                continue
            score, witness = model.proposal((value, value))
            print('SHARED DIAGNOSTIC mean', value, 'log2 proposal', score, flush=True)
            upper = model.outward((value, value), witness)
            print('SHARED DIAGNOSTIC outward log2', upper.log()/scalar_cover.arb(2).log(), flush=True)
            rows.append(dict(mean=str(value), proposal=score, witness=witness,
                             upper=[int(x) for x in upper.upper().man_exp()]))
        record.update(probes=rows, updates=args.updates, minimum_groups=args.minimum_groups,
                      distance=args.distance, base_tilt=args.base_tilt,
                      variance_bins=args.variance_bins, regional_count=args.regional_count)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    main()
