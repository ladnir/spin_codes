"""Three prescribed directed controls for the scheduled-Hlambda dense bridge.

Fix r=18 and total constituent arity A>=216. For each of the three marks,
rebuild the uniform-amplitude regional bank at alpha=1/5. The reviewed
Hlambda adapter supplies the construction binding and dependency checks.
All 99 indexed even/odd source profiles are summed in chronological order.
No old receipt or forced-nonzero regional measure supplies arithmetic.

For this fixed source size, changing A0 changes only the density factor.
If U216 is the newly computed log bound, the same witness bounds A>=A0 by
U216 - 8*alpha*(A0-216)*log2(255). Directed evaluation at every integer
A0 in216..288 identifies the first cutoff whose ceiling is at most-60.
This is the smallest cutoff certified by this computed witness, not a
claim that the true probability fails below that cutoff.

The marks are alternatives for the same class, not disjoint events.
Their bounds may be minimized, but must not be summed as a partition.
This exploratory control does not cover lower total arity, other r, or
the whole code. Importing it performs no provider or production run.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_high_arity_certificate as binding


DIRECTORY = Path(__file__).resolve().parent
ARITHMETIC = binding.arithmetic
R, MINIMUM_ARITY, MAXIMUM_ARITY = 18, 216, 288
ALPHA, TARGET = F(1, 5), F(-60)
MARKS = (F(691, 700), F(172, 175), F(171, 175))
PACKET_BITS = 60
ADAPTER_FILE = 'hd_lambda_high_arity_certificate.py'
ADAPTER_SHA256 = '156ca16d6f18370eb2484e6a1bf74ef83cbede670b63b2ba03db27f465d2e036'


def checked_dependencies():
    digest = sha256((DIRECTORY/ADAPTER_FILE).read_bytes()).hexdigest()
    if digest != ADAPTER_SHA256:
        raise ValueError('Reviewed Hlambda construction adapter changed')
    # Includes exact Hlambda involution/permutation checks and all current
    # numerical dependencies, including the paired-certificate dependency.
    return dict(binding.checked_dependencies(), **{ADAPTER_FILE: digest,
        Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def directed_cutoffs(base_log_bound):
    """Scalar density adjustment only; call inside the replay precision context."""
    if not base_log_bound.is_finite():
        raise ValueError('A finite fresh directed log-bound is required')
    step = ARITHMETIC.uniform.rational(8*ALPHA)*ARITHMETIC.log2(arb(255))
    rows = []
    for minimum in range(MINIMUM_ARITY, MAXIMUM_ARITY+1):
        value = base_log_bound-(minimum-MINIMUM_ARITY)*step
        ceiling = ARITHMETIC.uniform.rational_ceiling(value)
        rows.append(dict(total_arity_minimum=minimum, failure_log2_upper=str(ceiling),
                         reaches_60_bits=bool(ceiling <= TARGET)))
    selected = next((row for row in rows if row['reaches_60_bits']), None)
    return dict(density_improvement_per_added_row_log2_enclosure=str(step),
        cutoff_ceilings=rows,
        minimal_certified_total_arity=(selected['total_arity_minimum'] if selected else None),
        minimal_cutoff_failure_log2_upper=(selected['failure_log2_upper'] if selected else None),
        no_cutoff_in_requested_range=selected is None,
        minimality_scope='first directed ceiling <=-60 under this witness, not true-probability optimality')


def summarize_terms(mark, terms):
    """Validate full indexed-profile coverage before the directed source union."""
    if mark not in MARKS:
        raise ValueError('Only the three prescribed marks are permitted')
    counts, seen = ARITHMETIC.profile_counts(R), set()
    if len(counts) != 99:
        raise AssertionError('The fixed r18 source inventory must have99 profiles')
    prefix = ARITHMETIC.high_prefix(R, MINIMUM_ARITY, mark, ALPHA)
    rows, union = [], arb(0)
    for even, odd, count, contribution in terms:
        key = even, odd
        if (type(even) is not int or type(odd) is not int
                or key not in counts or key in seen or type(count) is not int
                or count != counts[key]):
            raise ValueError('Missing, duplicate, or incorrect indexed-source profile')
        if not contribution.is_finite() or not contribution > 0:
            raise ValueError('A positive finite fresh profile contribution is required')
        seen.add(key)
        union += contribution
        rows.append(dict(even_full=even, odd_full=odd, indexed_source_count_upper=count,
            failure_log2_upper=str(ARITHMETIC.uniform.rational_ceiling(
                prefix+ARITHMETIC.log2(contribution)))))
    if seen != set(counts):
        raise ValueError('Incomplete r18 source-profile union')
    value = prefix+ARITHMETIC.log2(union)
    ceiling = ARITHMETIC.uniform.rational_ceiling(value)
    return dict(mark=str(mark), alpha=str(ALPHA), total_arity_minimum=MINIMUM_ARITY,
        total_arity_maximum=MAXIMUM_ARITY, profiles=rows, profile_count=len(rows),
        failure_log2_upper=str(ceiling), failure_log2_enclosure=str(value),
        cutoff216_reaches_60_bits=bool(ceiling <= TARGET),
        **directed_cutoffs(value))


def receipt(precision=256, *, max_seconds_per_mark=60, progress=None):
    """Three fresh banks, serially; budget checked between numerical stages."""
    if (type(precision) is not int or precision not in (256, 384)
            or type(max_seconds_per_mark) not in (int, float)
            or not isfinite(max_seconds_per_mark) or not 0 < max_seconds_per_mark <= 60):
        raise ValueError('Use256/384 bits and a positive per-mark budget at most60 seconds')
    pins = checked_dependencies()
    started, controls, stop = perf_counter(), [], None
    attempted_banks = 0
    with ctx.workprec(precision):
        for mark in MARKS:
            mark_started = perf_counter()
            def check_time():
                if perf_counter()-mark_started >= max_seconds_per_mark:
                    raise TimeoutError('Dense-bridge mark exceeded its wall-time budget')
            try:
                if progress:
                    progress(dict(stage='fresh-regional-bank', mark=str(mark),
                                  precision_bits=precision))
                attempted_banks += 1
                regions = ARITHMETIC.uniform.regional_directed(mark, ALPHA, R, PACKET_BITS)
                check_time()
                terms = ARITHMETIC.source_profile_terms(R, regions, check_time=check_time)
                result = summarize_terms(mark, terms)
                check_time()
                result.update(fresh_regional_banks=1, elapsed_seconds=perf_counter()-mark_started)
                controls.append(result)
                if progress:
                    progress(dict(stage='mark-complete', mark=str(mark),
                        cutoff216_failure_log2_upper=result['failure_log2_upper'],
                        minimal_certified_total_arity=result['minimal_certified_total_arity'],
                        elapsed_seconds=result['elapsed_seconds']))
            except TimeoutError as error:
                stop = dict(mark=str(mark), type=type(error).__name__, message=str(error))
                break
    pins_after = checked_dependencies()
    if pins != pins_after:
        raise ValueError('Reviewed sources changed during the control')
    complete = len(controls) == len(MARKS) and stop is None
    best = min(controls, key=lambda row: F(row['failure_log2_upper'])) if controls else None
    scope = dict(binding.SCOPE)
    scope.pop('mark')
    scope.update(total_arity_minimum=MINIMUM_ARITY, total_arity_maximum=MAXIMUM_ARITY,
        prescribed_marks=[str(mark) for mark in MARKS], alpha=str(ALPHA),
        exploratory_dense_bridge_control=True, fixed_witness_search=True,
        threshold_log2=str(TARGET), cutoff_search=[MINIMUM_ARITY, MAXIMUM_ARITY],
        cutoff_change_uses_density_only=True, witness_choice_not_a_partition=True,
        no_extra_union_over_marks=True, minimal_cutoff_claim='certified by fixed-witness upper bound',
        precision_bits=precision, controls=controls, completed_controls=len(controls),
        attempted_fresh_regional_banks=attempted_banks, complete=complete, stop_reason=stop,
        max_seconds_per_mark=max_seconds_per_mark,
        best_completed_mark=best['mark'] if best else None,
        best_completed_cutoff216_log2_upper=best['failure_log2_upper'] if best else None,
        cutoff216_reaches_60_bits=bool(best and best['cutoff216_reaches_60_bits']),
        all_three_marks_tested=complete, old_source_receipts_used=False,
        reviewed_source_sha256=pins, reviewed_source_sha256_after=pins_after,
        source_sha256=pins[Path(__file__).name], source_pins_unchanged=True,
        elapsed_seconds=perf_counter()-started)
    return scope


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    parser.add_argument('--max-seconds-per-mark', type=float, default=60)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.precision, max_seconds_per_mark=args.max_seconds_per_mark,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
