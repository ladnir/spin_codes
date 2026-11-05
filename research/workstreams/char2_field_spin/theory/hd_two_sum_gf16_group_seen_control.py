"""Six saved-bank group-activity controls, restricted to source stratum F=1.

Use exactly the weights18/24/32, two stored type templates and frozen
regional slopes from state_cover_control. Reproduce every old empty-cut
baseline before evaluating the seen-group and aggregate9r constraints.
No provider, fit, new witness or directed arithmetic is invoked.
"""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_group_seen as seen
import hd_two_sum_gf16_state_cover_control as reference


DIRECTORY = Path(__file__).resolve().parent
PINS = {
    'hd_two_sum_gf16_group_seen.py': '9f049719e50f04a663ef3c4415fd34e43f81e28cd062312aac847d87dcb0e733',
    'hd_two_sum_gf16_state_cover_control.py': '0549c92786c5d2ea7a635014991d279ab3a73dd0490bea2f7809acd491739a06',
}


def checked_source():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed saved-bank control changed: {name}')
    return reference.load_source()


def run(max_seconds=30, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('An explicit at-most30s saved-bank budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Six saved-bank group-activity controls exhausted their budget')
    result = dict(construction=reference.low.fractional.CONSTRUCTION, field_order=256,
        forbidden_weight=reference.low.fractional.LIMIT, source_size=8192,
        source_bank=reference.BANK_NAME, source_bank_sha256=reference.BANK_SHA,
        dependency_pins=dict(PINS), inherited_reference_dependency_pins=dict(reference.PINS),
        source_full_pair_max=reference.FULL_MAX, source_stratum_only=True,
        all_source_supports=False, all_source_supports_in_declared_stratum=True,
        source_partition_complete=False, source_count_outside_alpha=True,
        fixed_types_across_all_columns=True, all_cohorts_per_anchor=True,
        all_constituent_row_sets_per_anchor=True, unchanged_literal_markers=True,
        unchanged_numeric_slopes=True, regional_banks_rebuilt=0, fitted_parameters=0,
        templates=list(reference.TEMPLATES), requested_controls=6,
        independent_seen_and_aggregate_controls=True, constraints_imposed_jointly=False,
        individual_nine_column_floors_enforced=False, best_with_old_baseline_reported=True,
        exact_total_arity_restriction_used=False, continuous_weight_band_covered=False,
        directed=False, whole_code_certificate=False, requires_fresh_directed_replay=True,
        compute_budget_seconds=max_seconds, witnesses=[], complete_requested_controls=False)
    try:
        check_time()
        source = checked_source()
        for bank, (r, _, mark) in zip(source['banks'], reference.low.ANCHORS):
            check_time()
            curves = reference.restricted_curves(bank)
            count = next(row['count'] for row in bank['source_strata']
                         if row['full_max'] == reference.FULL_MAX)
            slopes = np.asarray(bank['log2_slopes'])
            for metadata, weights in reference.fixed_templates(bank):
                check_time()
                result['in_progress'] = dict(active_groups=r, template=metadata['name'])
                prepared = reference.cover.emptycut.prepared_envelope(curves, slopes, weights)
                old = reference.restrict_count(reference.cover.emptycut.evaluate_prepared(
                    prepared, r, mark, reference.low.ALPHA), r, count)
                previous = [row for row in source['witnesses'] if row['active_groups'] == r
                    and row['source_full_pair_max'] == reference.FULL_MAX
                    and row['template']['name'] == metadata['name']]
                if (len(previous) != 1 or abs(old['failure_log2_screen']
                        -previous[0]['failure_log2_screen']) > 1e-7):
                    raise ArithmeticError('Frozen F1 empty-cut baseline failed its replay identity')
                row = seen.evaluate(prepared, weights, r, mark, reference.low.ALPHA,
                    source_count=count, check_time=check_time)
                best = min(old['failure_log2_screen'], row['failure_log2_screen'])
                row.update(template=metadata, source_full_pair_max=reference.FULL_MAX,
                    source_stratum_only=True, all_source_supports=False,
                    all_source_supports_in_declared_stratum=True, source_supports=count,
                    previous_emptycut_failure_log2_screen=old['failure_log2_screen'],
                    best_with_old_failure_log2_screen=best,
                    improvement_over_old_bits=old['failure_log2_screen']-best,
                    log2_slopes=slopes.tolist())
                result['witnesses'].append(row)
                result.pop('in_progress', None)
                if progress:
                    progress(dict(active_groups=r, source_full_pair_max=reference.FULL_MAX,
                        template=metadata['name'], seen=row['seen_all_failure_log2_screen'],
                        aggregate=row['aggregate_9r_failure_log2_screen'],
                        previous=old['failure_log2_screen'], best_with_old=best,
                        elapsed_seconds=perf_counter()-started))
        check_time()
        selected = [min((row for row in result['witnesses'] if row['active_groups'] == r),
                        key=lambda row: row['best_with_old_failure_log2_screen'])
                    for r, _, _ in reference.low.ANCHORS]
        result.update(status='completed', complete_requested_controls=True, selected=selected,
            selected_three_strata_log2_screen=float(reference.low.state.logsum(
                [row['best_with_old_failure_log2_screen'] for row in selected])))
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen saved-bank control')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('witnesses', 'selected')}))
