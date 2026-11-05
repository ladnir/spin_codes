"""Positive-term diagnosis of the pinned full-a8 four-atom comparator.

No provider, new witness, optimizer, or certificate. Rebuild the pinned Q,
find its whole-word tail's peak degree, then max-plus backtrack one largest
column-coefficient word at that degree. Decompose selected Q coefficients
by regional path class and00 atom count. A comparator term is NOT an
actual failure probability lower bound.
"""
from collections import Counter
from hashlib import sha256
import json
from math import comb, expm1, log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_matrix_moment as frozen


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_matrix_moment_gate_float.json'
ARTIFACT_SHA = '0465d355fd56c5c0cdfabac206a377ed74f62bf638b78af354bb4a147a9b5783'
HELPER_SHA = 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67'


def largest_column_word(matrix, columns, degree, check_time=lambda: None):
    """Max-plus single positive word, retaining state/degree backpointers."""
    matrix = np.asarray(matrix)
    if (matrix.ndim != 3 or matrix.shape[:2] != (2, 2)
            or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix))
            or type(columns) is not int or columns < 1 or type(degree) is not int
            or not 0 <= degree <= columns*(matrix.shape[2]-1)):
        raise ValueError('Valid finite-or-zero matrix polynomial and reachable degree required')
    values = np.full((2, degree+1), -np.inf)
    values[0, 0] = 0.
    back = []
    local_degree = matrix.shape[2]-1
    for column in range(columns):
        check_time()
        new = np.full_like(values, -np.inf)
        prior = np.full((2, degree+1), -1, dtype=np.int8)
        last = np.full((2, degree+1), -1, dtype=np.int16)
        max_old = min(degree, column*local_degree)
        for initial in range(2):
            for final in range(2):
                for d in range(min(local_degree, degree)+1):
                    if not np.isfinite(matrix[initial, final, d]):
                        continue
                    length = min(max_old, degree-d)+1
                    candidate = values[initial, :length]+matrix[initial, final, d]
                    target = new[final, d:d+length]
                    take = candidate > target
                    target[take] = candidate[take]
                    prior[final, d:d+length][take] = initial
                    last[final, d:d+length][take] = d
        back.append((prior, last))
        values = new
    final = int(np.argmax(values[:, degree]))
    value = float(values[final, degree])
    if not np.isfinite(value):
        raise ValueError('No positive word reaches the requested degree')
    state, remaining, path = final, degree, []
    for column in range(columns-1, -1, -1):
        prior, last = back[column]
        initial, d = int(prior[state, remaining]), int(last[state, remaining])
        if initial < 0 or d < 0:
            raise ArithmeticError('Missing max-plus path backpointer')
        path.append(dict(column=column, initial=initial, final=state,
                         degree=d, coefficient_log2=float(matrix[initial, state, d])))
        state, remaining = initial, remaining-d
    if state != 0 or remaining != 0:
        raise ArithmeticError('Max-plus path does not start at zero state and degree')
    path.reverse()
    return dict(log2_term=value, final_state=final, total_degree=degree, columns=path)


def path_breakdown(weights, groups, measure, arity, initial, final, degree,
                   check_time=lambda: None):
    """Positive class terms of ONE Q_initial,final coefficient."""
    if initial not in (0, 1) or final not in (0, 1) or not 0 <= degree <= groups*arity:
        raise ValueError('Valid state endpoints and column degree required')
    p = measure['log2_probability']
    other = -np.inf if p == 0 else log2(-expm1(p*log(2)))
    cache, rows = {}, []
    for source, target, counts, multiplicity in frozen.path_classes(len(weights)-1):
        if source != initial or target != final:
            continue
        n00, n01, n10, n11 = counts
        state_mass = (log2(multiplicity)+n01*measure['log2_scale01']
            +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active in range(n00+1):
            check_time()
            if p == 0 and active != n00:
                continue
            key = (active, n01+n10, n11)
            if key not in cache:
                polynomial = frozen.group_polynomial(weights, counts, active, measure, arity)
                cache[key] = frozen.log_power(polynomial, groups, check_time)
            coefficient = cache[key][degree]
            if not np.isfinite(coefficient):
                continue
            atom = log2(comb(n00, active))+active*p+((n00-active)*other if active < n00 else 0.)
            rows.append(dict(counts=list(counts), active00=active,
                chronological_paths=multiplicity, log2_term=float(coefficient+atom+state_mass)))
    total = float(np.logaddexp2.reduce([row['log2_term'] for row in rows])) if rows else -np.inf
    rows.sort(key=lambda row: row['log2_term'], reverse=True)
    for row in rows:
        row['log2_fraction_of_coefficient'] = row['log2_term']-total
    return dict(initial=initial, final=final, degree=degree, coefficient_log2=total,
                positive_terms=len(rows), top_terms=rows[:10])


def run(*, max_seconds=30):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit budget at most30 seconds required')
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA or sha256(Path(frozen.__file__).read_bytes()).hexdigest() != HELPER_SHA:
        raise ValueError('Frozen full-class comparator changed')
    source = json.loads(raw)
    if (source.get('scope') != 'complete homogeneous arity8/F1 class'
            or source.get('complete_requested_controls') is not True):
        raise ValueError('Completed full homogeneous comparator required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Saved-bank dominant-term diagnostic exhausted its budget')
    result = dict(scope=source['scope'], source=ARTIFACT, source_sha256=ARTIFACT_SHA,
        helper_sha256=HELPER_SHA, raw_regional_banks_built=0, new_parameters=0,
        fitting=False, failure_probability_lower_bound=False, directed=False,
        whole_code_certificate=False, compute_budget_seconds=max_seconds,
        complete_requested_diagnostic=False)
    try:
        data, curves = frozen.load_saved()
        old = data['two_atom_envelope']
        measure = frozen.entry_measure(curves, old['log2_probability'], old['log2_nonzero_atom'])
        weights = frozen.raw_holder_weights()
        matrix, costs = frozen.column_matrix(weights, frozen.GROUPS, measure,
            arity=frozen.ARITY, check_time=check_time)
        word = frozen.matrix_poly_power(matrix, frozen.COLUMNS, check_time)
        degrees = np.logaddexp2.reduce(word[0], axis=0)
        tail = float(np.logaddexp2.reduce(degrees[frozen.FLOOR:]))
        baseline = next(row for row in source['controls'] if row['mode'] == 'polynomial')
        if abs(tail-baseline['positive_tail_log2']) > 1e-8:
            raise ArithmeticError('Frozen whole-word tail did not reproduce')
        peak = frozen.FLOOR+int(np.argmax(degrees[frozen.FLOOR:]))
        dominant = largest_column_word(matrix, frozen.COLUMNS, peak, check_time)
        counts = Counter((row['initial'], row['final'], row['degree']) for row in dominant['columns'])
        q00 = [(count, d) for (i, j, d), count in counts.items() if i == j == 0]
        targets = {(0, 0, frozen.FLOOR//frozen.COLUMNS)}
        if q00:
            targets.add((0, 0, max(q00)[1]))
        else:
            targets.add(max(counts, key=counts.get))
        breakdown = []
        for initial, final, degree in sorted(targets):
            row = path_breakdown(weights, frozen.GROUPS, measure, frozen.ARITY,
                initial, final, degree, check_time)
            if abs(row['coefficient_log2']-matrix[initial, final, degree]) > 1e-8:
                raise ArithmeticError('Regional path decomposition did not reproduce Q')
            breakdown.append(row)
        result.update(positive_tail_log2=tail, peak_tail_degree=peak,
            peak_degree_coefficient_log2=float(degrees[peak]),
            peak_degree_fraction_log2=float(degrees[peak]-tail),
            dominant_column_word=dominant,
            dominant_word_fraction_of_peak_log2=dominant['log2_term']-float(degrees[peak]),
            column_entry_degree_counts=[dict(initial=i, final=j, degree=d, count=count)
                for (i, j, d), count in sorted(counts.items())],
            column_path_breakdowns=breakdown, complete_requested_diagnostic=True, **costs)
        check_time()
    except TimeoutError as error:
        result['complete_requested_diagnostic'] = False
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['implementation_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen diagnosis')
    answer = run(max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(answer))
