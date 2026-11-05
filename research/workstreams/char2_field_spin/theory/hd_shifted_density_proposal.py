"""Recompute only the outer prefix of a saved FLOAT fixed-atlas proposal.

The normalized support sums and exact witness are preserved. No atlas is
built, no witness is fitted, and no saved FLOAT value becomes a certificate.
Only the established outer-density modes are accepted.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import isfinite, log2
from pathlib import Path

import numpy as np

from hd_grouped_packet_screen import LIMIT
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
    checked_outer_measure, log_density_float, bound_description)
from hd_shifted_swap_screen import CONSTRUCTION


def reweight_proposal(screen, outer_measure, *, threshold=None):
    """Update every recorded weight from its saved support sum and exact z/alpha."""
    outer_measure = checked_outer_measure(outer_measure)
    witness = screen.get('fixed_witness', {})
    if (screen.get('construction') != CONSTRUCTION
            or witness.get('construction') != CONSTRUCTION
            or witness.get('pair_swaps') is not True
            or witness.get('orientation_aggregation') != 'independent-fair-swaps'
            or screen.get('directed') is not False
            or screen.get('every_weight_evaluated') is not True
            or screen.get('method') != 'fixed-atlas-common-potential-run-convolution'):
        raise ValueError('A complete FLOAT fixed-atlas proposal for the swapped construction is required')
    saved_measure = checked_outer_measure(witness.get('outer_measure', DEFAULT_OUTER_MEASURE))
    if ('outer_measure' in screen
            and checked_outer_measure(screen['outer_measure']) != saved_measure):
        raise ValueError('Top-level and witness outer-bound identities must agree')
    try:
        first, last = screen['requested_weight_range']
        rows = screen['results']
        exact = witness['source_exact_witness']
        # Exact literals, rather than the approximate convenience fields, define the prefix.
        if any(not isinstance(exact[name], str) for name in ('z', 'alpha', 'potential')):
            raise ValueError('Exact rational witness strings are required')
        z, alpha, potential = (F(exact[name]) for name in ('z', 'alpha', 'potential'))
        markers = tuple(tuple(F(value) for value in row.split()) for row in exact['markers'])
        if (not 0 < z < 1 or not 0 < alpha <= 1 or potential <= 0
                or not markers or len(markers) > 100
                or any(len(row) != 2 or min(row) <= 0 for row in markers)):
            raise ValueError('Admissible exact fixed-witness literals are required')
        threshold = screen['threshold'] if threshold is None else threshold
        if (type(first) is not int or type(last) is not int or not 1 <= first <= last <= 8192
                or not isinstance(rows, list) or screen['weights_evaluated'] != last-first+1
                or [row['active_groups'] for row in rows] != list(range(first, last+1))
                or any(type(row['active_groups']) is not int for row in rows)
                or not isfinite(threshold)
                or any(not isfinite(row['log2_support_sum_screen']) for row in rows)):
            raise ValueError('Finite support sums for every integer weight exactly once are required')
        z, alpha = float(z), float(alpha)
        if not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Exact z/alpha must have admissible finite FLOAT representations')
    except (KeyError, TypeError, ZeroDivisionError, OverflowError, AttributeError) as error:
        raise ValueError('Malformed exact witness or complete per-weight FLOAT rows') from error

    result = deepcopy(screen)
    density = log_density_float(outer_measure)
    constant = -alpha*(log2(255)+LIMIT*log2(z))
    for row in result['results']:
        prefix = alpha*row['active_groups']*density+constant
        row.update(log2_outer_prefix=prefix,
                   log2_upper_screen=row['log2_support_sum_screen']+prefix,
                   directed=False, whole_code_certificate=False)
    values = [row['log2_upper_screen'] for row in result['results']]
    if not all(map(isfinite, values)):
        raise ArithmeticError('Nonfinite reweighted FLOAT bound')
    union = float(np.logaddexp2.reduce(values))
    failed = [row['active_groups'] for row in result['results'] if row['log2_upper_screen'] > threshold]
    result.update(threshold=threshold, failed_weights=failed,
                  all_weights_pass_threshold=not failed, union_log2_upper_screen=union,
                  union_passes_threshold=union <= threshold, directed=False, whole_code_certificate=False,
                  density_reweighted_from=saved_measure, prefix_recomputed_from_exact_witness=True,
                  new_fft_tables_built=0, new_witness_searches=0, new_central_marker_fits=0)
    for record in (result, result['fixed_witness']):
        if outer_measure == DEFAULT_OUTER_MEASURE:
            record.pop('outer_measure', None)
            record.pop('outer_measure_bound', None)
        else:
            record.update(outer_measure=outer_measure,
                          outer_measure_bound=bound_description(outer_measure))
    return result


if __name__ == '__main__':
    import argparse, json
    from hashlib import sha256
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--outer-measure', required=True,
                        choices=(DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE))
    parser.add_argument('--threshold', type=float, help='Default: preserve the source threshold')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    raw = args.source.read_bytes()
    result = reweight_proposal(json.loads(raw), args.outer_measure, threshold=args.threshold)
    result['source_screen_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, separators=(',', ':'), allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()),
                          source_screen_sha256=result['source_screen_sha256'],
                          requested_weight_range=result['requested_weight_range'],
                          weights_evaluated=result['weights_evaluated'],
                          union_log2_upper_screen=result['union_log2_upper_screen'],
                          union_passes_threshold=result['union_passes_threshold'],
                          directed=False)), flush=True)
