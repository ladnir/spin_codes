"""FLOAT regional transfer with the existing sharpened joint-profile count.

The sparse proof already uses profile_count_upper. This diagnostic applies
the same bound to the alias matrices: minimum of both matching marginals
and the cyclic-run count, with exact contracted counts when a matching has
no singleton active groups. No construction or coefficient law changes.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaln

import hd_field128_gamma_transfer_diagnostic as transfer
from hd_shifted_pairs_count import profile_count_upper


def run_logs(length, weight):
    result = np.full(weight+1, -np.inf)
    if weight in (0, length):
        result[0] = 0.
    else:
        c = np.arange(1, min(weight, length-weight)+1)
        result[c] = (np.log2(length/c)+(gammaln(weight)-gammaln(c)-gammaln(weight-c+1)
            +gammaln(length-weight)-gammaln(c)-gammaln(length-weight-c+1))/np.log(2.))
    return result


def profile_log_counts(pairs, weight, first, second):
    first, second = np.broadcast_arrays(np.asarray(first), np.asarray(second))
    runs = weight-first-second
    raw = run_logs(2*pairs, weight)[runs]
    def marginal(full):
        singles, empty = weight-2*full, pairs-weight+full
        return (gammaln(pairs+1)-gammaln(full+1)-gammaln(singles+1)-gammaln(empty+1))/np.log(2.)+singles
    result = np.minimum(raw, np.minimum(marginal(first), marginal(second)))
    if weight % 2 == 0:
        contracted = np.full(weight+1, -np.inf)
        small = run_logs(pairs, weight//2)
        contracted[:len(small)] = small
        result = np.where((2*first == weight) | (2*second == weight), contracted[runs], result)
    return result


def support_sum(pairs, weight, regions, *, repetitions=256, chunk=32, budget=None):
    full = transfer.tail.full_counts(pairs, weight)
    regions = np.asarray(regions)
    if (regions.shape != (len(full), 2, 2) or type(repetitions) is not int
            or repetitions < 2 or repetitions % 2 or not 1 <= chunk <= 128):
        raise ValueError('Invalid regional matrices, repetitions, or chunk size')
    budget = budget or transfer.gamma.old.Budget(60.)
    total, largest = -np.inf, None
    for start in range(0, len(full), chunk):
        budget.check()
        left, right = full[start:start+chunk, None], full[None, :]
        counts = profile_log_counts(pairs, weight, left, right)
        product = transfer.multiply(regions[start:start+chunk, None], regions[None])
        repeated = transfer.power(product, repetitions//2)
        output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
        terms = counts+output
        index = np.unravel_index(np.argmax(terms), terms.shape)
        candidate = float(terms[index])
        if largest is None or candidate > largest['log2_term']:
            even, odd = int(full[start+index[0]]), int(full[index[1]])
            largest = dict(even_full=even, odd_full=odd, runs=weight-even-odd,
                           log2_term=candidate, log2_count_upper=float(counts[index]),
                           log2_transfer=float(output[index]))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(terms.ravel())))
    budget.check()
    return dict(log2_support_sum_screen=total, largest_profile=largest,
                profile_pairs_evaluated=len(full)**2, common_potential_used=False,
                support_count='existing profile_count_upper: marginals, runs, singleton-free contraction')


def source_pins():
    pins = transfer.source_pins()
    directory = Path(__file__).resolve().parent
    for name in (Path(__file__).name, 'hd_shifted_pairs_count.py'):
        pins[name] = sha256((directory/name).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--state', type=int, choices=transfer.gamma.STATES, required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be unused')
    raw, pins = args.source.read_bytes(), source_pins()
    artifact = json.loads(raw)
    row = artifact.get('best', artifact)
    if row['active_groups'] != args.r:
        parser.error('Source selected weight must match')
    w = row['witness']
    budget = transfer.gamma.old.Budget(args.max_seconds)
    bank = transfer.gamma.AliasBankScreen(float(F(w['z'])), float(F(w['alpha'])), w['marker_log_odds'],
             row['log2_potential_ratio'], state=args.state, distance=F(args.distance), grid=w['grid'], _budget=budget)
    regions = bank.regional(args.r)
    result = support_sum(transfer.gamma.old.PAIRS, args.r, regions, budget=budget)
    prefix = bank.alpha*bank.density*args.r+bank.constant
    result.update(bank.metadata(), active_groups=args.r, target_distance=str(F(args.distance)),
                  log2_outer_prefix=prefix, log2_upper_screen=prefix+result['log2_support_sum_screen'],
                  log2_potential_ratio=row['log2_potential_ratio'], regional_alias_diagnostics=bank.regional_diagnostics,
                  witness=bank.witness(), source=dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest()),
                  source_sha256=pins, old_numerical_bounds_reused=False,
                  scope='one selected weight with sharpened counts and chronological transfer; FLOAT only')
    if source_pins() != pins or args.source.read_bytes() != raw:
        raise RuntimeError('Pinned source changed')
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('active_groups', 'state_dimension', 'log2_upper_screen',
                                                  'largest_profile', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
