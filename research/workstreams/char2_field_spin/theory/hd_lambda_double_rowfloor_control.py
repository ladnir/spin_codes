"""One original-row RS floor for the failing double-transvection A=0..5.

The pinned baseline already marks every quiet input weight exactly. Its
aggregate six-row floor nevertheless permits activity to migrate between
original rows. Here original row0 must be active in at least nine columns.
For r moment replicas, retain total row0 activity at least9r.

For each fixed scheduled column, the zero-outside fibre is pointwise quiet.
Thus row0 activity is exactly its physical output coordinate's activity.
Every nonquiet input is assigned activity1 before the existing unmarked
coset comparison. This is upper-tail domination, not coefficientwise equality.
The zero vector remains activity0; a nonzero quiet vector may also score0.

Original rows0 and1 occur together in all quiet supports. The nonquiet
surrogate activates both. This score is therefore exactly half the proposed
two-row score, with the same retained event. Its degree is1, floor162 and
deficit126. The actual source arity remains6: the RS factor is255^(-48)
per group, not255^(-8). Positive entrywise rth powers precede the ordered
sixteen-column matrix product. Only complete same-A word bounds are minimized.

This is a bounded FLOAT diagnostic using a pinned, outer-independent regional
measure. It performs no new regional calculation, marker fit, or native test.
It does not cover other row sets or certify the whole code.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_middle_control as base


DIRECTORY = Path(__file__).resolve().parent
SOURCE = tuple(range(6))
GROUPS = 18
BASELINE = 'hd_lambda_double_middle_representatives_float.json'
BASELINE_SHA = '5f6cf543bc57646a6aba273c8f8ab8888d27106b3b489018edb23447f10e28f3'
BASE_SHA = '78655f83f08ae01118beb6d9f6258daaddd7a5035bc5837694b6243d762e9820'


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_double_middle_control.py').read_bytes()).hexdigest() != BASE_SHA:
        raise ValueError('Frozen double-transvection controller changed')
    return {**base.dependencies(),
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def load_baseline():
    raw = (DIRECTORY/BASELINE).read_bytes()
    if sha256(raw).hexdigest() != BASELINE_SHA:
        raise ValueError('Pinned representative receipt changed')
    saved = json.loads(raw)
    if (saved.get('construction') != base.fiber.outer.CONSTRUCTION
            or saved.get('schedule') != base.fiber.SCHEDULE or saved.get('groups') != GROUPS
            or saved.get('source_class') != 'Fmax<=1' or saved.get('mark') != '9927/10000'
            or saved.get('mode') != 'representatives' or saved.get('directed') is not False
            or saved.get('complete') is not True or saved.get('source_pins_unchanged') is not True
            or saved.get('reviewed_source_sha256') != base.dependencies()
            or saved.get('input_artifact') != base.INPUT or saved.get('input_sha256') != base.INPUT_SHA):
        raise ValueError('Baseline is not the frozen complete same-construction FLOAT gate')
    matches = [row for row in saved['controls'] if tuple(row['original_source']) == SOURCE]
    if len(matches) != 1 or matches[0]['status'] != 'complete' or matches[0]['arity'] != 6:
        raise ValueError('Baseline needs one complete original A=0..5 bound')
    row = matches[0]
    root = row['normalized_fixed_A_root_log2']
    if (not isfinite(root) or abs(saved['common_prefix_log2']+GROUPS*root
            -row['fixed_A_repeated18_bound_log2']) > 1e-9):
        raise ValueError('Baseline fixed-A normalization is inconsistent')
    return saved, row


def original_score(mask, data, label=0):
    """Exact selected original-coordinate activity on this pointwise quiet fibre."""
    if (label not in SOURCE or tuple(data['original_source']) != SOURCE
            or not data['all_zero_fiber_coordinates_protected']
            or data['physical_protected_mask'] != data['physical_inside_mask']):
        raise ValueError('The fixed A=0..5 pointwise-quiet fibre is required')
    if type(mask) is not int or mask < 0 or mask & ~data['physical_inside_mask']:
        raise ValueError('Quiet output mask must lie in the physical inside set')
    coordinate = data['output_permutation'][data['input_permutation'][label]]
    return int(bool(mask & (1 << coordinate)))


def prepare_columns(table_r, table_2r, *, check_time=None):
    columns, nonquiet, diagnostics = [], [], []
    for j in range(16):
        if check_time: check_time()
        data = base.fiber.column_fiber_counts(SOURCE, j)
        bounds = base.nonquiet_norm(data, table_r, table_2r, GROUPS, check_time=check_time)
        coefficients = np.full((2, 2, 2), -np.inf)
        scores = []
        for mask, count in data['zero_fiber_counts'].items():
            score = original_score(mask, data, 0)
            if original_score(mask, data, 1) != score:
                raise ValueError('The two-row score no longer admits exact degree halving')
            coefficients[score] = np.logaddexp2(coefficients[score], log2(count)+table_r[mask]/GROUPS)
            scores.append(dict(physical_output_mask=mask, exact_count=count, original_row0_activity=score))
        # Replace every nonquiet input's row0 activity by1 BEFORE coset domination.
        coefficients[1] = np.logaddexp2(coefficients[1], bounds['selected'])
        columns.append(coefficients)
        nonquiet.append(bounds['selected'])
        diagnostics.append(dict(column=j, mapped_source=data['mapped_source'],
            cross_rank=data['cross_rank'], syndrome_rank=data['syndrome_rank'],
            physical_row0=data['output_permutation'][data['input_permutation'][0]],
            physical_row1=data['output_permutation'][data['input_permutation'][1]],
            all_zero_fiber_coordinates_protected=True, quiet_scores=scores,
            direct_union_terms=bounds['direct_terms'], factor_norm_log2=bounds['factor_2r'].tolist(),
            direct_norm_log2=bounds['direct_union'].tolist()))
    return np.stack(columns), np.stack(nonquiet), diagnostics


def run(*, max_seconds=120, progress=None):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('A positive wall budget of at most120 seconds is required')
    pins, started = dependencies(), perf_counter()
    baseline, old = load_baseline()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded original-row floor control exhausted its wall budget')
    saved, measure, prefix = base.load_measure()
    if (saved['source_count'] != baseline['source_count'] or saved['mark'] != baseline['mark']
            or measure != baseline['common_measure'] or prefix != baseline['common_prefix_log2']):
        raise ValueError('New and baseline word bounds need the same common measure and prefix')
    complete, error, row, diagnostic, details = False, None, None, None, []
    try:
        table_r, table_2r = base.moment_table(measure, 18), base.moment_table(measure, 36)
        check_time()
        columns, nonquiet, details = prepare_columns(table_r, table_2r, check_time=check_time)
        diagnostic = base.all_nonquiet_diagnostic(nonquiet, 6, prefix)
        diagnostic.update(surrogate_degree=16*18, survives_floor=9*18,
                          score_original_rows=[0], actual_source_arity=6)
        if progress:
            progress(base._jsonable(dict(stage='all-nonquiet-before-tail', diagnostic=diagnostic,
                elapsed_seconds=perf_counter()-started)))
        tail = float(base.positive_deficit_tail(columns, GROUPS, 9*GROUPS, check_time=check_time))
        root = tail/GROUPS-8*6*log2(255)
        best = min(root, old['normalized_fixed_A_root_log2'])
        row = dict(original_source=SOURCE, actual_source_arity=6, score_original_rows=[0],
            word_tail_log2=tail, normalized_fixed_A_root_log2=root,
            fixed_A_repeated18_bound_log2=prefix+GROUPS*root,
            selected_complete_fixed_A_root_log2=best,
            selected_complete_fixed_A_bound_log2=prefix+GROUPS*best,
            selected_method='row0-floor' if root < old['normalized_fixed_A_root_log2'] else 'baseline-all-row-floor',
            min_taken_after_complete_chronological_word_tail=True)
        complete = True
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins:
        raise ValueError('Reviewed sources changed during the control')
    if sha256((DIRECTORY/base.INPUT).read_bytes()).hexdigest() != base.INPUT_SHA:
        raise ValueError('Pinned common regional artifact changed during the control')
    if sha256((DIRECTORY/BASELINE).read_bytes()).hexdigest() != BASELINE_SHA:
        raise ValueError('Baseline receipt changed during the control')
    return base._jsonable(dict(construction=base.fiber.outer.CONSTRUCTION,
        supports=base.fiber.outer.SUPPORTS, lambda_aes_byte=188, schedule=base.fiber.SCHEDULE,
        original_source=SOURCE, original_source_fixed_across16columns=True,
        actual_source_arity=6, groups=GROUPS, source_class='Fmax<=1',
        source_count=saved['source_count'], mark=saved['mark'], common_measure=measure,
        common_prefix_log2=prefix, input_artifact=base.INPUT, input_sha256=base.INPUT_SHA,
        baseline_artifact=BASELINE, baseline_sha256=BASELINE_SHA,
        baseline_fixed_A_bound_log2=old['fixed_A_repeated18_bound_log2'],
        score_original_rows=[0], equivalent_two_row_scores=[1, 1, 0, 0, 0, 0],
        exact_degree_halving_checked=True, retained_score_floor=162, retained_score_deficit=126,
        RS_factor_per_group='255^(-48)', RS_factor_whole_word='255^(-864)',
        nonquiet_score_inflated_before_unmarked_comparison=True,
        coefficientwise_true_input_activity_equality_claimed=False,
        complete=complete, error=error, result=row, column_diagnostics=details,
        all_nonquiet_diagnostic=diagnostic,
        all_nonquiet_diagnostic_is_method_lower_bound_not_failure=True,
        reviewed_source_sha256=pins, source_pins_unchanged=True,
        regional_provider_runs=0, input_marker_searches=0, directed=False,
        old_outer_receipts_transferred=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
