"""One fresh Fmax1 regional bank and two frozen late-alpha quiet controls.

Only r18/arity8/z691/700/alpha1/5. For a fixed source set G with the maximum
of its even/odd full-pair counts equal to1, every regional subset has f<=1.
Source counts come from the pinned exact matching partition, once outside
alpha. The local envelope does NOT assume fresh source geometry per mask.

Both controls cover every FIXED arity8 row-set vector in this F1 source
class, but only the zero-moment/all-regional-state0 majorant contribution.
They do not cover all paths/syndromes or give a failure lower bound.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_factorized import forced_nonzero_regions
from hd_two_sum_late_alpha_quiet import (ALPHA, ARITY, COLUMNS, DIMENSION,
    GROUPS, MARK, MINIMUM, compare_contributions, zero_intercept_envelope)
from hd_two_sum_late_alpha_moment import quiet_mixture_coefficients, two_atom_envelope, word_tail
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
PARTITION = 'hd_two_sum_gf16_source_strata_anchors_float.json'
PARTITION_SHA = '5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf'
FULL_MAX = 1


def checked_source_count(data):
    if (data.get('source_size') != 8192 or data.get('field_order') != 256
            or data.get('forbidden_weight') != LIMIT
            or data.get('complete_requested_anchors') is not True):
        raise ValueError('The complete pinned source partition is required')
    banks = [bank for bank in data.get('banks', []) if bank.get('active_groups') == GROUPS]
    if len(banks) != 1 or (banks[0].get('mark'), banks[0].get('alpha')) != (str(MARK), str(ALPHA)):
        raise ValueError('Exactly the frozen r18 witness partition is required')
    rows = banks[0].get('source_strata', [])
    if (len(rows) != GROUPS//2+1 or any(type(row.get('full_max')) is not int
            or type(row.get('count')) is not int or row['count'] <= 0 for row in rows)):
        raise ValueError('Unique positive integer counts for every F are required')
    counts = {row['full_max']: row['count'] for row in rows}
    if (set(counts) != set(range(GROUPS//2+1))
            or sum(counts.values()) != comb(8192, GROUPS)):
        raise ValueError('Source strata must partition every size18 source set')
    return counts[FULL_MAX]


def load_source_count():
    raw = (DIRECTORY/PARTITION).read_bytes()
    if sha256(raw).hexdigest() != PARTITION_SHA:
        raise ValueError('Pinned source partition changed')
    return checked_source_count(json.loads(raw))


def restrict_regions(regions, maximum=GROUPS, full_max=FULL_MAX):
    """Replicate each F-restricted matrix maximum for the line helper.

    The complete returned layout is an ENVELOPE, not a claim that the
    omitted f histograms occur in the fixed source class. All four entries
    are enveloped, although this diagnostic subsequently uses only00.
    """
    expected = {(m, f) for m in range(maximum+1) for f in range(m//2+1)}
    if (type(maximum) is not int or maximum < 1 or type(full_max) is not int
            or not 0 <= full_max <= maximum//2 or set(regions) != expected):
        raise ValueError('Complete regional histograms and a valid fixed source maximum required')
    answer = {}
    for m in range(maximum+1):
        selected = np.asarray([regions[m, f] for f in range(min(full_max, m//2)+1)])
        if (selected.shape != (min(full_max, m//2)+1, 2, 2)
                or np.any(np.isnan(selected)) or np.any(np.isposinf(selected))):
            raise ValueError('Finite-or-zero regional log matrices required')
        envelope = np.max(selected, axis=0)
        for f in range(m//2+1):
            answer[m, f] = envelope
    return answer


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    count = load_source_count()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Restricted F1 late-alpha gate exhausted its budget')
    prefix = (log2(count)+GROUPS*log2(comb(16, ARITY))
              +float(ALPHA)*(-LIMIT*log2(float(MARK))-log2(255)))
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=FULL_MAX,
        mark=str(MARK), alpha=str(ALPHA), field_order=256, forbidden_weight=LIMIT,
        columns=COLUMNS, constituent_dimension=DIMENSION,
        scope='Fmax1 sources; zero moments; all regional boundary states zero',
        all_sources=False, all_sources_in_declared_stratum=True,
        all_fixed_arity_rowset_vectors=True, rowsets_resampled=False,
        source_count=str(count), rowset_count=f'{comb(16, ARITY)}^{GROUPS}',
        source_count_once=True, source_and_rowset_counts_outside_alpha=True,
        partition_artifact=PARTITION, partition_sha256=PARTITION_SHA,
        partition_FLOAT_bounds_used=False, fresh_regional_alpha=1,
        regional_full_pair_restriction='f<=min(1,floor(m/2))',
        regional_geometry_resampled=False, raw_regional_banks=0,
        requested_controls=2, moment_secant_candidates=1, mark_fits=0,
        aggregate_input_weight_floor=MINIMUM*ARITY*GROUPS,
        individual_row_floors=False, all_rows_seen_constraint=False,
        directed=False, whole_code_certificate=False,
        full_quiet_event_bounded=False, failure_probability_lower_bound=False,
        compute_budget_seconds=max_seconds, complete_requested_controls=False)
    try:
        check_time()
        raw = forced_nonzero_regions(float(MARK), GROUPS)
        result['raw_regional_banks'] = 1
        check_time()
        regions = restrict_regions(raw)
        line = zero_intercept_envelope(regions, GROUPS)
        result['line_envelope'] = line
        if progress:
            progress(dict(stage='F1-R1-built', elapsed_seconds=perf_counter()-started,
                log2_slope=line['log2_slope'], contacts=line['contacts'],
                saved_F1_log2_curves=line['log2_curves']))
        control = compare_contributions(GROUPS, ARITY, line['log2_slope'], float(ALPHA), check_time=check_time)
        result['line_control'] = dict(control, log2_global_prefix=prefix,
            late_majorant_term_log2=prefix+control['late_column_log2'],
            early_majorant_term_log2=prefix+control['early_mask_log2'],
            late_minus_early_bits=control['late_column_log2']-control['early_mask_log2'])
        curve = np.asarray([-np.inf if value is None else value for value in line['log2_curves']])
        moments, envelope = two_atom_envelope(curve)
        result['two_atom_envelope'] = envelope
        coefficients = quiet_mixture_coefficients(GROUPS, ARITY, moments, check_time=check_time)
        tail = word_tail(coefficients, check_time=check_time)
        result['moment_control'] = dict(normalized_word_log2=tail,
            majorant_term_log2=prefix+tail,
            old_matched_F1_late_majorant_log2=prefix+control['late_column_log2'],
            moment_minus_line_bits=tail-control['late_column_log2'],
            coefficientwise_moment_balancing=True,
            synthetic_weights_one_and_two_admitted=True,
            raw_column_coefficient_degree81_log2=float(coefficients[81]))
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['dependency_sha256'] = {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in (
        'hd_two_sum_late_alpha_quiet.py', 'hd_two_sum_late_alpha_moment.py',
        'hd_two_sum_gf16_factorized.py', 'hd_two_sum_closed_kernel.py')}
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', default=30, type=float)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    result = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(line=result.get('line_control'), moment=result.get('moment_control'),
        complete=result['complete_requested_controls'], elapsed_seconds=result['elapsed_seconds'])))
