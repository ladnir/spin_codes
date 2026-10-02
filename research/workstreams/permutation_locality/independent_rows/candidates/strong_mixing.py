"""Screen stronger mixing against the ten-window local refinements.

Historical memo input is explicitly diagnostic only: its stored source hash
is reported, not authenticated against today's generators. Rebuild mode
regenerates local bounds, but the global binary64 search is still not a cert.
"""
import argparse
import json
from pathlib import Path

from flint import arb, ctx

import local_family
import mass_density_screen as screen
from mass_density import blend
from mixing_attack import retarget
from occupancy_memory import Z, C


def variants(base, coefficients, spectrum, tilt, rounds):
    # Retarget BEFORE blending: the unsplit M/L -> Z entries contain refresh
    # only. Blending inserts lazy contributions and destroys that invariant.
    changed = retarget(base, spectrum, tilt, 2, rounds)
    candidate = {j: tuple(x * arb(2) ** (2-rounds) for x in row)
                 for j, row in coefficients.items()}
    yield 'unsplit', changed
    for start in (5, 6, 7):
        fractions = {j: screen.Q(1) for j in range(start, 17)}
        zero = blend(changed, candidate, fractions, target=Z)
        yield f'zero-{start}-16', zero
        yield f'both-{start}-16', blend(zero, candidate, fractions, target=C)


def historical(directory):
    """Checksum-check saved data without claiming source authentication."""
    for path in sorted(Path(directory).glob('local-family-*.json')):
        key = json.loads(path.read_bytes())['body']['key']
        expected = dict(rounds=2, dimension=11, windows=32, maximum=16,
                        penalty='9/10', tilt='9/125', exact_feedback=10,
                        density_through=10, joint_four=True, input_weight='1')
        if any(key.get(k) != v for k, v in expected.items()):
            continue
        if local_family.cache_path(directory, key) != path:
            raise ValueError('historical filename does not match its key')
        ctx.prec = max(192, key['precision'])
        family = local_family.load(directory, key)
        if family is None:
            raise ValueError('historical record vanished')
        print('HISTORICAL DIAGNOSTIC ONLY; source NOT authenticated:', path,
              'sources', key['sources'], 'stage', key['stage'],
              'density_chord', key.get('density_chord'), flush=True)
        yield path.stem, family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--historical-directory')
    source.add_argument('--rebuild', action='store_true')
    parser.add_argument('--groups', type=int, nargs='+', default=[96, 128])
    parser.add_argument('--supports', type=int, nargs='+', default=[200])
    parser.add_argument('--rounds', type=int, nargs='+', default=[2, 3, 4, 8])
    args = parser.parse_args()
    if (any(not 1 <= q <= 2048 for q in args.groups)
            or any(not 38 <= u <= 256 for u in args.supports)
            or any(not 2 <= r <= 32 for r in args.rounds)):
        parser.error('invalid groups, supports, or updates')
    print('SELECTED-POINT BINARY64 SEARCH, NOT A CERTIFICATE', flush=True)
    if args.rebuild:
        families = [('regenerated', local_family.build(
            ['.072'], '.9', exact_feedback=10, joint_four=True,
            density_through=10)['.072'])]
    else:
        families = historical(args.historical_directory)
    spectrum = screen.baseline.prepare_inputs()[3]
    caps = screen.baseline.authenticated_caps()
    cdf = screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(
        caps, 1 << 128, full_weight=screen.Q(10, 9)))
    shells = screen.baseline.weighted_union_shells(caps, full_weight=screen.Q(10, 9))
    counts = {u: min(screen.Q(cdf[u]), shells[u]) for u in args.supports}
    for label, (base, coefficients) in families:
        for rounds in args.rounds:
            for variant, operators in variants(base, coefficients, spectrum, '.072', rounds):
                region = screen.float_placement(screen.as_array(operators), max(args.groups))
                for q in args.groups:
                    for u in args.supports:
                        value, p = screen.score(region, q, u, counts[u], '.072', cutoff=193986)
                        print('STRONG MIXING', label, variant, 'R/q/u', rounds, q, u,
                              'log2', value, 'p', p, flush=True)


if __name__ == '__main__':
    main()
