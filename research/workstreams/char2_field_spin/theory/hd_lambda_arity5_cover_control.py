"""Bounded FLOAT coverage of all original arity-five row sets at r18/F<=1.

There are 4368 labelled original sets A. For each A, the marked-fibre
method retains quiet input weight and inflates nonquiet weight to five.
Positive deficit truncation imposes aggregate RS input activity >=9*5*18.
At most 128 weak fixed-A bounds receive an exact GF4-category fallback.
Each method completes its entire ordered word tail before their roots are
minimized. Their activity polynomials are never minimized coefficientwise.

The all-A sum is prefix + 18*logsum_A(2**root_A), with
root_A=logtail_A/18-40*log2(255). No additional row-set count is multiplied.
The source/output-mark/projective prefix appears once outside that sum.

The full root sum permits a different labelled five-row set in every
active group. It does not cover arbitrary arity mixtures of total arity90.

This controller records partial roots on timeout. A partial sum covers only
its recorded subset, not all arity-five sources. A full selected-class bound
is exposed only after exact labelled coverage and unchanged source pins are
checked. Even a complete result is FLOAT, not a directed certificate or a
whole-code claim. No regional provider or native work is performed.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_column_deficit_tail import column_deficit_tail
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_lambda_fiber_control import moment_table
from hd_lambda_marked_fiber_control import prepare_source, UnsupportedQuietFiber
from hd_lambda_middle_control import exact_columns
from hd_pair_scheduled_control import load_measure


DIRECTORY = Path(__file__).resolve().parent
ARITY, GROUPS, ROWS = 5, 18, 16
SOURCE_COUNT = comb(ROWS, ARITY)
TARGET_LOG2 = -60.
FALLBACK_THRESHOLD = TARGET_LOG2-GROUPS*log2(SOURCE_COUNT)
RS_ROOT_CHARGE = 8*ARITY*log2(255)
SOURCE_FILES = (
    'hd_lambda_arity5_cover_control.py', 'hd_column_deficit_tail.py',
    'hd_lambda_middle_control.py', 'hd_lambda_arity5_relations.py',
    'hd_quaternary_outer_relations.py', 'hd_lambda_scheduled_relations.py',
    'hd_lambda_marked_fiber_control.py', 'hd_lambda_fiber_control.py',
    'hd_lambda_fiber_outer.py', 'hd_lambda_outer.py', 'hd_outer_pair_layers.py',
    'hd_bit_syndrome_outer.py', 'hd_bit_syndrome_flat_geometry.py',
    'gf8_packet_maps.py', 'hd_pair_scheduled_outer.py', 'hd_pair_scheduled_control.py',
    'hd_pair_scalar_control.py', 'hd_inverse_column_control.py',
    'hd_pair_identity_kernel.py', 'hd_pair_state_kernel.py',
    'hd_post_refresh_quiet_rows.py', 'hd_scheduled_column_control.py',
    'hd_pair_state_tail_general.py', 'hd_pair_state_tail.py')


def original_sources():
    """The complete deterministic inventory, with no orbit compression."""
    return tuple(combinations(range(ROWS), ARITY))


def source_mask(source):
    source = tuple(source)
    if (source != tuple(sorted(set(source))) or len(source) != ARITY
            or any(type(i) is not int or not 0 <= i < ROWS for i in source)):
        raise ValueError('A sorted original labelled set of exactly five rows is required')
    return sum(1 << i for i in source)


def _log_value(value):
    value = float(value)
    if not isfinite(value) and value != -np.inf:
        raise ArithmeticError('Invalid FLOAT bound: NaN or positive infinity')
    return value


def coverage_summary(records, prefix):
    """Check exact full coverage before distinguishing full and subset sums.

    This does not authorize incomplete sums as bounds for omitted row sets.
    It also rechecks each chosen root as the minimum of complete word roots.
    """
    expected = {source_mask(source) for source in original_sources()}
    seen, roots = set(), []
    duplicate_masks, invalid_records = [], []
    for position, row in enumerate(records):
        try:
            mask = source_mask(row['original_source'])
            if row['original_mask'] != mask or mask not in expected:
                raise ValueError('Wrong original mask')
            marked = _log_value(row['marked_normalized_root_log2'])
            chosen = _log_value(row['normalized_fixed_A_root_log2'])
            exact = row.get('exact_normalized_root_log2')
            wanted = min(marked, _log_value(exact)) if exact is not None else marked
            if chosen != wanted:
                raise ValueError('Chosen root is not a complete-method minimum')
            if mask in seen:
                duplicate_masks.append(mask)
            seen.add(mask)
            roots.append(chosen)
        except (KeyError, TypeError, ValueError, ArithmeticError) as error:
            invalid_records.append(dict(position=position, error=str(error)))
    exact_coverage = (len(records) == SOURCE_COUNT and seen == expected
                      and not duplicate_masks and not invalid_records)
    subset_root = float(np.logaddexp2.reduce(roots)) if roots else -np.inf
    # On malformed/duplicated input, withhold even the subset interpretation.
    subset_value = (prefix+GROUPS*subset_root
                    if not duplicate_masks and not invalid_records else None)
    return dict(expected_original_A=SOURCE_COUNT, recorded_rows=len(records),
        unique_original_A=len(seen), all_original_A_covered=exact_coverage,
        missing_original_masks=sorted(expected-seen), duplicate_original_masks=duplicate_masks,
        invalid_records=invalid_records, covered_subset_root_sum_log2=subset_root,
        covered_subset_bound_log2=subset_value,
        covered_subset_bound_scope='recorded labelled A only; omitted A are NOT bounded',
        full_selected_arity_bound_log2=subset_value if exact_coverage else None)


def _pin_files():
    return {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def _input_hash():
    return sha256((DIRECTORY/INPUT).read_bytes()).hexdigest()


def _process_sources(sources, table, table2, prefix, *, batch_size,
                     fallback_limit, check_time, progress=None):
    """Internal bounded loop; tests inject tiny preparations, not real banks."""
    records, fallbacks = [], []
    fallback_skipped = 0
    stop_reason = None
    pending_preparations = 0
    try:
        for offset in range(0, len(sources), batch_size):
            chunk = sources[offset:offset+batch_size]
            preparations = []
            for source in chunk:
                check_time()
                columns, _ = prepare_source(source, table, table2, groups=GROUPS,
                                             check_time=check_time)
                preparations.append(columns)
                pending_preparations = len(preparations)
            check_time()
            # Axes: column, fixed original A, activity degree, state, state.
            columns = np.stack(preparations, axis=1)
            tails = np.asarray(column_deficit_tail(columns, ARITY, GROUPS,
                                                   check_time=check_time), dtype=float)
            if tails.shape != (len(chunk),):
                raise ArithmeticError('Batched deficit tail returned an unexpected shape')
            new_records = []
            for source, tail in zip(chunk, tails):
                tail = _log_value(tail)
                root = tail/GROUPS-RS_ROOT_CHARGE
                row = dict(original_source=list(source), original_mask=source_mask(source),
                    marked_tail_log2=tail, marked_normalized_root_log2=root,
                    normalized_fixed_A_root_log2=root,
                    fixed_A_repeated_r_bound_log2=prefix+GROUPS*root,
                    selected_method='marked-fibre', exact_fallback_status='not_needed')
                new_records.append(row)
            # Commit the whole valid marked batch before any costly fallback.
            records.extend(new_records)
            pending_preparations = 0
            for row in new_records:
                if row['fixed_A_repeated_r_bound_log2'] <= FALLBACK_THRESHOLD:
                    continue
                if len(fallbacks) >= fallback_limit:
                    row['exact_fallback_status'] = 'attempt_cap_skipped'
                    fallback_skipped += 1
                    continue
                attempt = dict(original_source=row['original_source'],
                    original_mask=row['original_mask'], status='started',
                    marked_fixed_A_bound_log2=row['fixed_A_repeated_r_bound_log2'])
                fallbacks.append(attempt)
                row['exact_fallback_status'] = 'started'
                try:
                    check_time()
                    exact, diagnostic = exact_columns(tuple(row['original_source']), table, GROUPS,
                                                       check_time=check_time)
                    tail_array = np.asarray(column_deficit_tail(exact, ARITY, GROUPS,
                                                                check_time=check_time), dtype=float)
                    if tail_array.shape != ():
                        raise ArithmeticError('Exact fixed-A tail must be scalar')
                    tail = _log_value(tail_array)
                    root = tail/GROUPS-RS_ROOT_CHARGE
                    row.update(exact_tail_log2=tail, exact_normalized_root_log2=root,
                               exact_fallback_status='complete')
                    if root < row['normalized_fixed_A_root_log2']:
                        row['normalized_fixed_A_root_log2'] = root
                        row['selected_method'] = 'exact-GF4'
                    row['fixed_A_repeated_r_bound_log2'] = prefix+GROUPS*row['normalized_fixed_A_root_log2']
                    attempt.update(status='complete', exact_tail_log2=tail,
                        exact_normalized_root_log2=root, selected_method=row['selected_method'],
                        selected_fixed_A_bound_log2=row['fixed_A_repeated_r_bound_log2'],
                        diagnostics=diagnostic)
                except TimeoutError:
                    row['exact_fallback_status'] = 'timeout_marked_retained'
                    attempt['status'] = 'timeout_marked_retained'
                    raise
            if progress and len(records) % 256 == 0:
                progress(dict(stage='covered_original_A', covered=len(records), total=SOURCE_COUNT,
                    exact_attempts=len(fallbacks), exact_cap_skipped=fallback_skipped,
                    largest_selected_fixed_A_bound_log2=max(
                        row['fixed_A_repeated_r_bound_log2'] for row in records)))
    except (TimeoutError, UnsupportedQuietFiber) as error:
        stop_reason = dict(type=type(error).__name__, message=str(error))
    return dict(records=records, exact_fallbacks=fallbacks,
                exact_fallback_cap_skipped=fallback_skipped,
                pending_preparations_without_tail=pending_preparations,
                stop_reason=stop_reason)


def _jsonable(value):
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and not isfinite(value):
        if value == -np.inf:
            return None
        raise ArithmeticError('NaN/+infinity cannot be exported as proof evidence')
    return value


def run(*, batch_size=16, max_seconds=360, exact_fallback_limit=128, progress=None):
    """Run the complete inventory with strict production caps; no provider."""
    if (type(batch_size) is not int or batch_size not in (16, 32)
            or not isfinite(max_seconds) or not 0 < max_seconds <= 360
            or type(exact_fallback_limit) is not int or not 0 <= exact_fallback_limit <= 128):
        raise ValueError('Batch16/32, at most360 seconds, and at most128 exact fallbacks required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('All-A5 control exceeded its wall-time budget')
    pins_before, raw_before = _pin_files(), _input_hash()
    if raw_before != INPUT_SHA:
        raise ValueError('Pinned regional artifact changed before the run')
    saved, measure, prefix = load_measure()
    work = dict(records=[], exact_fallbacks=[], exact_fallback_cap_skipped=0,
                pending_preparations_without_tail=0, stop_reason=None)
    try:
        check_time()
        table, table2 = moment_table(measure, GROUPS), moment_table(measure, 2*GROUPS)
        check_time()
        work = _process_sources(original_sources(), table, table2, prefix,
            batch_size=batch_size, fallback_limit=exact_fallback_limit,
            check_time=check_time, progress=progress)
    except TimeoutError as error:
        work['stop_reason'] = dict(type=type(error).__name__, message=str(error))
    pins_after, raw_after = _pin_files(), _input_hash()
    unchanged = pins_before == pins_after and raw_before == raw_after == INPUT_SHA
    coverage = coverage_summary(work['records'], prefix)
    complete = coverage['all_original_A_covered'] and work['stop_reason'] is None and unchanged
    full_bound = coverage['full_selected_arity_bound_log2'] if complete else None
    # Withhold the full result even if an early-ended run happened to record
    # every marked root. The partial roots remain available for explicit reuse.
    coverage['full_selected_arity_bound_log2'] = full_bound
    result = dict(construction='scheduled-Hlambda-bc-no-cross', lambda_aes_byte=0xbc,
        schedule='input-inverse-output-xor', groups=GROUPS, arity=ARITY,
        source_class='Fmax<=1', mark=saved['mark'], source_count=saved['source_count'],
        common_measure=measure, common_prefix_log2=prefix,
        original_row_arity_per_active_group=ARITY,
        heterogeneous_labelled_five_row_choices_permitted=True,
        arbitrary_arity_mixtures_with_total_arity90_covered=False,
        original_A_fixed_across_columns=True, full_independent_D1=True,
        raw_field_order=256, rs_root_normalization='255**(-40)',
        aggregate_input_activity_floor=9*ARITY*GROUPS, retained_deficit=7*ARITY*GROUPS,
        original_per_row_activity_floors_relaxed=True,
        minimum_only_after_complete_fixed_A_word_tail=True,
        coefficientwise_minimum_between_exact_and_surrogate=False,
        extra_rowset_multiplier=False, source_projective_prefix_included_once=True,
        batch_size=batch_size, max_seconds=max_seconds,
        exact_fallback_limit=exact_fallback_limit, exact_fallback_trigger_log2=FALLBACK_THRESHOLD,
        fallback_budget_exhausted=bool(work['exact_fallback_cap_skipped']),
        refinement_cap_does_not_stop_marked_coverage=True,
        target_selected_arity_log2=TARGET_LOG2,
        target_met=bool(complete and full_bound <= TARGET_LOG2), complete=bool(complete),
        coverage=coverage, **work, input_artifact=INPUT, input_sha256=INPUT_SHA,
        source_pins_before=pins_before, source_pins_after=pins_after,
        source_pins_unchanged=unchanged, raw_artifact_sha256_after=raw_after,
        arithmetic='positive-FLOAT-log2-not-directed', null_log2_means_exact_zero=True,
        directed=False, whole_code_certificate=False, regional_provider_runs=0,
        elapsed_seconds=perf_counter()-start)
    if progress:
        progress(dict(stage='finished', covered=coverage['unique_original_A'], total=SOURCE_COUNT,
            complete=bool(complete), target_met=result['target_met'],
            full_selected_arity_bound_log2=full_bound,
            exact_attempts=len(work['exact_fallbacks']), elapsed_seconds=result['elapsed_seconds']))
    return _jsonable(result)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--batch-size', type=int, choices=(16, 32), default=16)
    parser.add_argument('--max-seconds', type=float, default=360)
    parser.add_argument('--exact-fallback-limit', type=int, default=128)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(batch_size=args.batch_size, max_seconds=args.max_seconds,
        exact_fallback_limit=args.exact_fallback_limit,
        progress=lambda row: print(json.dumps(_jsonable(row), allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
