"""FLOAT diagnosis of the common-potential relaxation, not a certificate.

For every even/odd matching profile, retain the complete chronological
2-by-2 regional transfer product. The old cyclic-run count still upper
bounds the number of indexed supports with that profile. Only the final
common-potential contraction is removed; local and coefficient caps stay.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaln

import hd_field128_full_rs_gamma_screen as gamma
import hd_field128_lift_tail_screen as tail


def multiply(left, right):
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def power(matrix, exponent):
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Expected a nonnegative integer exponent')
    matrix = np.asarray(matrix, dtype=float)
    result = np.full_like(matrix, -np.inf)
    result[..., 0, 0] = result[..., 1, 1] = 0.
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = multiply(matrix, matrix)
    return result


def support_sum(pairs, weight, regions, *, repetitions=256, chunk=32, budget=None):
    full = tail.full_counts(pairs, weight)
    regions = np.asarray(regions)
    if (regions.shape != (len(full), 2, 2) or type(repetitions) is not int
            or repetitions < 2 or repetitions % 2 or not 1 <= chunk <= 128):
        raise ValueError('Invalid regional matrices, repetitions, or row chunk')
    budget = budget or gamma.old.Budget(60.)
    length = 2*pairs
    runs = np.arange(weight+1)
    counts = np.full(weight+1, -np.inf)
    if weight in (0, length):
        counts[0] = 0.
    else:
        valid = (runs >= 1) & (runs <= min(weight, length-weight))
        c = runs[valid]
        counts[valid] = (np.log2(length/c)
            +(gammaln(weight)-gammaln(c)-gammaln(weight-c+1)
              +gammaln(length-weight)-gammaln(c)-gammaln(length-weight-c+1))/np.log(2.))
    total, largest = -np.inf, None
    for start in range(0, len(full), chunk):
        budget.check()
        product = multiply(regions[start:start+chunk, None, :, :], regions[None, :, :, :])
        repeated = power(product, repetitions//2)
        output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
        run_count = weight-full[start:start+chunk, None]-full[None, :]
        terms = counts[run_count]+output
        index = np.unravel_index(np.argmax(terms), terms.shape)
        candidate = float(terms[index])
        if largest is None or candidate > largest['log2_term']:
            largest = dict(even_full=int(full[start+index[0]]), odd_full=int(full[index[1]]),
                           runs=int(run_count[index]), log2_term=candidate,
                           log2_count_upper=float(counts[run_count[index]]),
                           log2_transfer=float(output[index]))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(terms.ravel())))
    budget.check()
    return dict(log2_support_sum_screen=total, largest_profile=largest,
                profile_pairs_evaluated=len(full)**2, common_potential_used=False,
                support_count='inherited cyclic-run count upper bound per matching profile')


def source_pins():
    directory, pending, seen = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--state', type=int, choices=gamma.STATES, required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    raw = args.source.read_bytes()
    artifact = json.loads(raw)
    row = artifact.get('best', artifact)
    if row['active_groups'] != args.r:
        parser.error('Source must be a selected-weight row at the requested weight')
    w = row['witness']
    pins = source_pins()
    budget = gamma.old.Budget(args.max_seconds)
    bank = gamma.AliasBankScreen(float(F(w['z'])), float(F(w['alpha'])), w['marker_log_odds'],
                                row['log2_potential_ratio'], state=args.state, distance=F(args.distance),
                                grid=w['grid'], _budget=budget)
    regions = bank.regional(args.r)
    result = support_sum(gamma.old.PAIRS, args.r, regions, budget=budget)
    prefix = bank.alpha*bank.density*args.r+bank.constant
    common = tail.support_sum_logs(gamma.old.PAIRS, args.r, regions,
                    row['log2_potential_ratio'], method='runs', _budget=budget)
    result.update(bank.metadata(), active_groups=args.r, target_distance=str(F(args.distance)),
                  log2_outer_prefix=prefix, log2_upper_screen=prefix+result['log2_support_sum_screen'],
                  common_potential_control_log2_upper=prefix+common,
                  log2_potential_ratio=row['log2_potential_ratio'], regional_alias_diagnostics=bank.regional_diagnostics,
                  witness=bank.witness(), source=dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest()),
                  source_sha256=pins, old_numerical_bounds_reused=False,
                  scope='all matching profiles at one weight; FLOAT transfer diagnostic, not a directed certificate')
    if source_pins() != pins or args.source.read_bytes() != raw:
        raise RuntimeError('Source changed during diagnosis')
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('active_groups', 'state_dimension', 'log2_upper_screen',
                    'common_potential_control_log2_upper', 'largest_profile', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
