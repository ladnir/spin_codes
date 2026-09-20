"""Audit coverage, aggregation and provenance of all 130 BCH evidence rows."""
import csv
import json
import math
from pathlib import Path

import numpy as np

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
KEYS = ('block_bits', 'step_bits', 'state_bits', 'message_exponent')
LN2 = math.log(2)


def main():
    hashes = {}

    def sha(path):
        path = path.resolve()
        if path not in hashes:
            hashes[path] = study.sha(path)
        return hashes[path]

    def authenticate(data):
        for name, digest in data['source_sha256'].items():
            if sha(study.ROOT/name) != digest:
                raise ValueError(f'Stale evidence dependency: {name}')

    def key(row):
        return tuple(int(row[k]) for k in KEYS)

    receipt_path = HERE/'bch_engineering_evidence_v8.json'
    receipt = json.loads(receipt_path.read_text()); authenticate(receipt)
    csv_path = HERE/'bch_engineering_evidence_v8.csv'
    if sha(csv_path) != receipt['csv_sha256']:
        raise ValueError('Evidence CSV changed')
    rows = list(csv.DictReader(csv_path.open(newline='')))
    original_receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    authenticate(original_receipt)
    if sha(HERE/'bch_dominance_sparse.csv') != original_receipt['csv_sha256']:
        raise ValueError('Original grid changed')
    baseline = list(csv.DictReader((HERE/'bch_dominance_sparse.csv').open(newline='')))
    if len(rows) != 130 or len({key(r) for r in rows}) != len(rows) or {key(r) for r in rows} != {key(r) for r in baseline}:
        raise ValueError('The report does not cover the complete original BCH grid')
    baseline_index = {key(r): r for r in baseline}
    lower_checks = {}
    for version in (2, 3):
        path = HERE/f'bch_zero_state_verification_v{version}.json'
        data = json.loads(path.read_text()); authenticate(data)
        if data['status'] != 'ALL_POSITIVE_GRID_WITNESSES_REPLAYED' or data['decimal_digits'] < 90:
            raise ValueError('Expected high-precision positive lower-witness replay')
        for check in data['checks']:
            lower_checks.setdefault(key(check), []).append(check['lower_bits'])
    full_count = 0; lower_count = 0; useful = []; weak = []; errors = []
    for row in rows:
        b, t, s, e = geometry = key(row)
        if abs(float(row['q1_margin_bits'])-float(baseline_index[geometry]['q1_margin_bits'])) > 1e-8:
            raise ValueError('The report changed the original Q1 grid')
        if row['evidence'] == 'FIRST_MOMENT_OBSTRUCTION':
            lower = float(row['first_moment_lower_bits'])
            if lower <= 0 or not any(abs(lower-value) < 1e-6 for value in lower_checks.get(geometry, [])):
                raise ValueError('Obstruction lacks a matching positive replay')
            gap = lower+float(row['q1_margin_bits'])
            expected = gap+math.log(-math.expm1(-gap*LN2))/LN2
            if abs(expected-float(row['higher_to_q1_log2_lower_bound'])) > 1e-7:
                raise ArithmeticError('Obstruction ratio lower bound changed')
            if float(row['full_margin_upper_from_first_moment_lower_bits']) != -lower:
                raise ArithmeticError('First-moment margin obstruction changed')
            lower_count += 1
            continue
        path = HERE/f'bch_full_reference_b{b}_t{t}_s{s}_e{e}.json'
        if not path.exists():
            raise ValueError(f'Uncovered geometry: {geometry}')
        full = json.loads(path.read_text()); authenticate(full)
        if full['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE' or key(full) != geometry:
            raise ValueError('Incompatible full reference')
        components = full['components']; length = (1 << e)//(b//2)
        intervals = [(c['occupation_min'], c['occupation_max']) for c in components]
        if (intervals[:4] != [(1, 1), (2, 2), (3, 3), (4, 4)] or intervals[-1][1] != length
                or any(lo > hi for lo, hi in intervals)
                or any(a[1]+1 != z[0] for a, z in zip(intervals, intervals[1:]))):
            raise ValueError('Full reference has an incomplete occupation interval cover')
        q1 = components[0]['log_upper']
        tail = float(np.logaddexp.reduce([c['log_upper'] for c in components[1:]]))
        total = float(np.logaddexp(q1, tail)); margin = -total/LN2
        penalty = float(np.logaddexp(0., tail-q1))/LN2
        discrepancies = [abs(margin-float(row['full_margin_bits'])),
                         abs(penalty-float(row['full_aggregation_penalty_bits'])),
                         abs((tail-q1)/LN2-float(row['full_higher_to_q1_log2_ratio'])),
                         abs(-q1/LN2-float(row['q1_margin_bits']))]
        if max(discrepancies) > 1e-7:
            raise ArithmeticError('Grid aggregation or ratio does not match the complete reference')
        expected = ('FULL_BOUND_UNINFORMATIVE' if margin <= 0 else
                    'FULL_BOUND_Q1_DOMINANT' if full['margin_penalty_bits'] <= math.log2(1.1) else
                    'FULL_BOUND_LARGER_PENALTY')
        if row['evidence'] != expected:
            raise ValueError('Incorrect evidence classification')
        errors.extend(discrepancies); full_count += 1
        result = dict(zip(KEYS, geometry), full_margin_bits=margin, margin_penalty_bits=penalty,
                      higher_to_q1_log2_ratio=(tail-q1)/LN2)
        if margin <= 0:
            result['q1_q4_margin_bits'] = -float(np.logaddexp.reduce(
                [c['log_upper'] for c in components[:4]]))/LN2
            result['nonpositive_component_margins'] = [
                dict(occupation_min=c['occupation_min'], occupation_max=c['occupation_max'],
                     margin_bits=-c['log_upper']/LN2)
                for c in components if c['log_upper'] >= 0]
        (useful if margin > 0 else weak).append(result)
    model_path = HERE/'bch_k_counting_model_verification_v2.json'
    model = json.loads(model_path.read_text()); authenticate(model)
    if model['status'] != 'VERIFIED_90_DIGIT_COUNT_SCALING':
        raise ValueError('Missing count-scaling model replay')
    authenticate(json.loads((HERE/'bch_k_counting_model_v4.json').read_text()))
    if {(c['block_bits'], c['comparison_exponent']) for c in model['checks']} != {(b, e) for b in (64, 128) for e in range(21, 27)}:
        raise ValueError('The large-K model is not checked at every grid exponent')
    counts = {label: sum(r['evidence'] == label for r in rows) for label in sorted({r['evidence'] for r in rows})}
    if counts != receipt['evidence_counts'] or full_count+lower_count != len(rows):
        raise ValueError('Coverage totals do not reconcile')
    result = dict(status='VERIFIED_COMPLETE_BCH_EVIDENCE_GRID', geometries=len(rows),
                  full_reference_geometries=full_count, first_moment_obstructions=lower_count,
                  useful_full_bounds=len(useful), weak_full_bounds=len(weak), missing_geometries=0,
                  evidence_counts=counts, maximum_aggregation_error_bits=max(errors),
                  largest_useful_margin_penalty_bits=max(r['margin_penalty_bits'] for r in useful),
                  q1_majority_of_selected_bound=sum(r['higher_to_q1_log2_ratio'] < 0 for r in useful),
                  weak_geometries=weak, authenticated_files=len(hashes),
                  source_sha256={str(p): sha(p) for p in (Path(__file__), receipt_path, csv_path, model_path)},
                  limitations=['Complete evidence coverage does not imply a useful full bound at every point.',
                               'Weak upper bounds do not prove failure or Q1 dominance.',
                               'Obstructions concern first moments, not failure probabilities.',
                               'Numerical bound replays are authenticated; no outward certificate is asserted.'])
    (HERE/'bch_evidence_grid_audit_v1.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('weak_geometries', 'source_sha256', 'limitations')}, indent=2))


if __name__ == '__main__':
    main()
