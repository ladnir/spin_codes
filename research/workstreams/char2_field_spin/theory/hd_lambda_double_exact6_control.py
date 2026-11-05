"""Bounded selected a6 gate using exact input-weight/output-support counts.

All eighteen active groups use the fixed original row set A=0..5.
The new enumerator replaces the coset/nonquiet majorant by exact raw
GF256 counts in each scheduled column. Minkowski is applied within each
input degree; coefficientwise Holder then retains total input weight>=972.
The pinned r18/Fmax<=1 regional measure and RS factor255^-864 are unchanged.
This is a FLOAT diagnostic for one class, not a whole-code certificate.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_exact6 as exact
import hd_lambda_double_middle_control as tail_backend
import hd_lambda_double_rowfloor_control as baseline


DIRECTORY = Path(__file__).resolve().parent
PINS = {
    'hd_lambda_double_exact6.py': 'c6065ca7f6b9605aec096283d8ff2db94630692f99a4d2deed9e20b30d4dc960',
    'hd_lambda_double_middle_control.py': '78655f83f08ae01118beb6d9f6258daaddd7a5035bc5837694b6243d762e9820',
    'hd_lambda_double_rowfloor_control.py': 'de46c03f2633a35de3fd12ed2f3c608d5b5001d77c3c2c2fb2401bed00c96da4',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**exact.dependencies(), **tail_backend.dependencies(), **baseline.dependencies(),
            **PINS, Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def coefficients(counts, table, *, groups=18):
    """Raw multiplicities stay outside the norm root."""
    result = np.full((7, 2, 2), -np.inf)
    for (weight, mask), count in counts.items():
        if not 0 <= weight <= 6 or not 0 <= mask < 65536 or type(count) is not int or count <= 0:
            raise ValueError('Positive exact raw joint counts are required')
        result[weight] = np.logaddexp2(result[weight], log2(count)+table[mask]/groups)
    return result


def run(*, max_seconds=120, progress=None):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Bounded wall budget at most120 seconds required')
    pins, started = dependencies(), perf_counter()
    saved, old = baseline.load_baseline()
    bank, measure, prefix = tail_backend.load_measure()
    if measure != saved['common_measure'] or prefix != saved['common_prefix_log2']:
        raise ValueError('Baseline and new bound require the same common measure')
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Exact a6 gate exhausted its wall budget')
    complete, error, diagnostics, result = False, None, None, None
    try:
        counts, diagnostics = exact.exact_columns(tuple(range(6)), check_time=check_time)
        diagnostics['ordered_joint_counts_sha256'] = sha256(json.dumps(
            [[[w, mask, count] for (w, mask), count in sorted(column.items())] for column in counts],
            separators=(',', ':')).encode()).hexdigest()
        diagnostics['joint_terms_per_column'] = [len(column) for column in counts]
        if progress: progress(dict(stage='exact-columns-complete', **diagnostics))
        table = tail_backend.moment_table(measure, 18)
        columns = np.stack([coefficients(column, table) for column in counts])
        tail = float(tail_backend.positive_deficit_tail(columns, 18, 972, check_time=check_time))
        root = tail/18-48*log2(255)
        result = dict(word_tail_log2=tail, normalized_fixed_A_root_log2=root,
            fixed_A_repeated18_bound_log2=prefix+18*root,
            baseline_fixed_A_bound_log2=old['fixed_A_repeated18_bound_log2'])
        complete = True
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins:
        raise ValueError('Sources changed during the exact gate')
    baseline.load_baseline()
    return tail_backend._jsonable(dict(construction=exact.outer.CONSTRUCTION,
        original_source=list(range(6)), actual_source_arity=6, groups=18,
        original_source_fixed_across16columns=True, source_class='Fmax<=1',
        source_count=bank['source_count'], mark=bank['mark'], common_prefix_log2=prefix,
        exact_joint_counts=True, nonquiet_coset_majorant_used=False,
        aggregate_original_input_floor=972, RS_factor_whole_word='255^(-864)',
        raw_counts_outside_norm_roots=True, complete=complete, error=error,
        enumeration=diagnostics, result=result, source_pins=pins,
        common_measure_receipt=tail_backend.INPUT, common_measure_sha256=tail_backend.INPUT_SHA,
        baseline_receipt=baseline.BASELINE, baseline_sha256=baseline.BASELINE_SHA,
        directed=False, other_original_rowsets_covered=False, whole_code_certificate=False,
        native_operations=0, regional_provider_runs=0, elapsed_seconds=perf_counter()-started))


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
    print(json.dumps(dict(complete=result['complete'], result=result['result'],
        error=result['error'], elapsed_seconds=result['elapsed_seconds']), allow_nan=False))
