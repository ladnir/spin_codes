"""Read the frozen post-refresh gate and locate its largest positive terms.

This diagnoses an upper bound, not actual codewords or a failure lower bound.
There is no provider call or fitted parameter. All chronology is retained.
"""
from collections import Counter
from hashlib import sha256
import json
from math import comb, expm1, log, log2
from pathlib import Path

import numpy as np

import hd_coupled_state_moment_gate as engine
from hd_two_sum_matrix_moment_diagnostic import largest_column_word


DIRECTORY = Path(__file__).resolve().parent


def path_breakdown(weights, groups, measure, arity, initial, final, degree):
    p = measure['log2_probability']
    other = -np.inf if p == 0 else log2(-expm1(p*log(2)))
    cache, rows = {}, []
    for source, target, counts, multiplicity in engine.path_classes(len(weights)-1):
        if (source, target) != (initial, final):
            continue
        n00, n01, n10, n11 = counts
        scale = (log2(multiplicity)+n01*measure['log2_scale01']
                 +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active in range(n00+1):
            if p == 0 and active != n00:
                continue
            key = active, n01, n10, n11
            if key not in cache:
                poly = engine.group_polynomial(weights, counts, active, measure, arity)
                cache[key] = engine.log_power(poly, groups)
            coefficient = cache[key][degree]
            if not np.isfinite(coefficient):
                continue
            mass = (log2(comb(n00, active))+active*p
                    +((n00-active)*other if active < n00 else 0.))
            rows.append(dict(counts=list(counts), active00=active,
                chronological_paths=multiplicity, log2_term=float(coefficient+mass+scale)))
    total = float(np.logaddexp2.reduce([row['log2_term'] for row in rows])) if rows else -np.inf
    rows.sort(key=lambda row: row['log2_term'], reverse=True)
    return dict(initial=initial, final=final, degree=degree, coefficient_log2=total,
                positive_terms=len(rows), top_terms=rows[:5])


def run(path, expected_hash):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != expected_hash:
        raise ValueError('Frozen completed gate hash required')
    data = json.loads(raw)
    if (not data.get('complete_requested_control') or data.get('construction') !=
            'proposed-g16-post-refresh-coupled-inner' or data.get('whole_code_certificate')):
        raise ValueError('Completed post-refresh FLOAT gate required')
    for name, pin in data['reviewed_source_sha256'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != pin:
            raise ValueError('Gate dependency changed: '+name)
    weights, measure = engine.raw_holder_weights(), data['state_measure']
    matrix, _ = engine.column_matrix(weights, 18, measure, arity=8)
    word = engine.matrix_poly_power(matrix, 16, lambda: None)
    degrees = np.logaddexp2.reduce(word[0], axis=0)
    tail = float(np.logaddexp2.reduce(degrees[1296:]))
    if abs(tail-data['positive_tail_log2']) > 1e-8:
        raise ArithmeticError('Saved gate tail did not reproduce')
    peak = 1296+int(np.argmax(degrees[1296:]))
    dominant = largest_column_word(matrix, 16, peak)
    counts = Counter((row['initial'], row['final'], row['degree']) for row in dominant['columns'])
    breakdown = []
    for (initial, final, degree), count in sorted(counts.items()):
        row = path_breakdown(weights, 18, measure, 8, initial, final, degree)
        if abs(row['coefficient_log2']-matrix[initial, final, degree]) > 1e-8:
            raise ArithmeticError('Path terms do not reproduce the column coefficient')
        row['columns_in_largest_word'] = count
        breakdown.append(row)
    return dict(gate_sha256=expected_hash, raw_regional_banks_built=0,
        parameter_fits=0, whole_code_certificate=False, failure_probability_lower_bound=False,
        positive_tail_log2=tail, peak_tail_degree=peak,
        peak_degree_fraction_log2=float(degrees[peak]-tail),
        dominant_column_word=dominant, column_path_breakdowns=breakdown,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    parser.add_argument('--gate-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen diagnostic')
    result = run(args.gate, args.gate_sha256)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, allow_nan=False))
