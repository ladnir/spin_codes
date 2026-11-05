"""Independent cached-FLOAT recombination; no provider or contraction DP.

Arb treats every stored binary64 logarithm as its exact dyadic value.
It reconstructs the common constant functional and all distinct activity
residuals, then checks sampled complete words against a separate batched
chronological row-vector calculation. These checks do not turn the input
FLOAT column bounds into directed certificates.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np
from flint import arb, ctx

DIRECTORY = Path(__file__).resolve().parent
RECEIPT = 'hd_mds_exact_empty_mixed789_9927_pinned_float.json'
RECEIPT_SHA = 'c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8'
SCALAR = 'hd_mds_empty_weight_9927_exact.json'
SCALAR_SHA = '83df60adddbd2979372db1a9f74ace1542ce6febbd5dfdb318c8148632014b5d'
OLD = 'hd_bit_syndrome_mixed89_9927_x1_float.json'
OLD_SHA = 'd94ccf2b9a9bc073b67c90256608fa67ef7be02f44efc47eede301f627c4e2c3'
SOURCE_COUNT = 4230999418401723055235639193772735965409871930858905600


def checked_json(name, expected):
    raw = (DIRECTORY/name).read_bytes()
    if sha256(raw).hexdigest() != expected:
        raise ValueError('Frozen input changed: '+name)
    return json.loads(raw)


def rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def exp2(value):
    return arb(0) if value == -np.inf else (rational(float(value))*arb(2).log()).exp()


def mm(left, right):
    return [[sum((left[i][k]*right[k][j] for k in range(2)), arb(0))
             for j in range(2)] for i in range(2)]


def empty_functional(measure):
    # The two 00 masses sum to exactly one, not a rounded p+(1-p).
    step = [[arb(1), exp2(measure['log2_scale01'])],
            [exp2(measure['log2_scale10']), exp2(measure['log2_scale11'])]]
    result = [[arb(1), arb(0)], [arb(0), arb(1)]]
    for _ in range(16):
        result = mm(result, step)
    return result


def word_logs(column):
    zero, live = np.zeros(len(column)), np.full(len(column), -np.inf)
    for _ in range(16):
        zero, live = (np.logaddexp2(zero+column[:, 0, 0], live+column[:, 1, 0]),
                      np.logaddexp2(zero+column[:, 0, 1], live+column[:, 1, 1]))
    return np.logaddexp2(zero, live)


def arb_word(column, empty=None, marker=F(1)):
    matrix = [[exp2(column[i, j]) for j in range(2)] for i in range(2)]
    if marker != 1:
        factor = rational(marker)**18
        matrix = [[empty[i][j]+factor*(matrix[i][j]-empty[i][j])
                   for j in range(2)] for i in range(2)]
    row = [arb(1), arb(0)]
    for _ in range(16):
        row = [row[0]*matrix[0][j]+row[1]*matrix[1][j] for j in range(2)]
    return (row[0]+row[1]).log()/arb(2).log()


def run(max_seconds=5):
    started = perf_counter()
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Cached recheck budget exceeded')
    receipt = checked_json(RECEIPT, RECEIPT_SHA)
    scalar = checked_json(SCALAR, SCALAR_SHA)
    old = checked_json(OLD, OLD_SHA)
    expected = dict(construction='G16-RS16x8-D2-bit-syndrome-D1-research', groups=18,
                    arities=[7, 8, 9], source_class='Fmax<=1', source_count=str(SOURCE_COUNT),
                    mark='9927/10000', input_marker=1, common_measure=True,
                    directed=False, whole_code_certificate=False)
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError('Wrong selected-class scope')
    pins = dict(receipt['reviewed_sources'], **{'hd_mds_exact_empty_control.py': receipt['source_sha256']})
    for name, expected_sha in pins.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected_sha:
            raise ValueError('Dependency changed: '+name)
    provider = receipt['empty_packet_provider']
    if (provider['exact_cap'] != scalar['exact_empty11'] or scalar['mark'] != receipt['mark']
            or sha256(provider['exact_cap'].encode()).hexdigest() != provider['exact_cap_sha256']
            or not provider['reference_cone_unchanged'] or not provider['nonempty_packet_caps_unchanged']
            or scalar['source_sha256'] != pins['hd_mds_empty_weight.py']):
        raise ValueError('Exact empty-cap provenance mismatch')
    measure, curves = receipt['state_measure'], receipt['regional_log2_curves']
    if any(len(row) != 19 for row in curves) or len(curves) != 4:
        raise ValueError('Four complete regional curves required')
    p, b0, b1, b2, c1, c2, c3 = map(F, (measure[key] for key in
        ('log2_probability', 'log2_atom00', 'log2_atom01', 'log2_atom10',
         'log2_scale01', 'log2_scale10', 'log2_scale11')))
    moment_slacks, monotonic_slacks = [], []
    for state in range(4):
        for m in range(19):
            value, previous = curves[state][m], old['regional_log2_curves'][state][m]
            bound = (F(0) if m == 0 else p+b0*m, c1+b1*m, c2+b2*m, c3)[state]
            if value is not None:
                moment_slacks.append(bound-F(value))
                if previous is None:
                    raise ArithmeticError('New finite entry exceeds a structural zero')
                monotonic_slacks.append(F(previous)-F(value))
            elif previous is not None:
                monotonic_slacks.append(F(0))
    if min(moment_slacks) < -F(1, 10**9) or min(monotonic_slacks) < -F(1, 10**9):
        raise ArithmeticError('FLOAT moment or monotonicity defect exceeds roundoff scale')
    cache_path = Path(receipt['cache'])
    if not cache_path.is_absolute():
        cache_path = DIRECTORY/cache_path
    if sha256(cache_path.read_bytes()).hexdigest() != receipt['cache_sha256']:
        raise ValueError('Column cache changed')
    with np.load(cache_path, allow_pickle=False) as cache:
        columns = {a: cache['column'+str(a)].copy() for a in (7, 8, 9)}
        for a in columns:
            if (columns[a].shape != (comb(16, a), 2, 2)
                    or not np.array_equal(cache['rowsets'+str(a)], list(combinations(range(16), a)))
                    or not np.isfinite(columns[a]).all()):
                raise ValueError('Incomplete labelled column cache')
        # Independent reconstruction of the fixed a8 Young combination.
        combined8 = np.logaddexp2(cache['quiet8']-17*np.log2(1/8),
                                 cache['nonquiet8']-17*np.log2(7/8))
        branch_error = float(np.max(np.abs(combined8-columns[8])))
        if branch_error > 1e-10:
            raise ArithmeticError('Arity8 Young branches disagree')
    check_time()
    with ctx.workprec(384):
        ln2 = arb(2).log()
        empty = empty_functional(measure)
        marker, corrected = F(5, 2), columns[7].copy()
        min_relative = arb(1)
        for i in range(2):
            for j in range(2):
                values, inverse = np.unique(columns[7][:, i, j], return_inverse=True)
                replacements = []
                for value in values:
                    q = exp2(value)
                    residual = q-empty[i][j]
                    if not residual > 0:
                        raise ArithmeticError('Activity residual not rigorously positive')
                    relative = residual/q
                    if relative < min_relative:
                        min_relative = relative
                    replacements.append(float((empty[i][j]+rational(marker)**18*residual).log()/ln2))
                corrected[:, i, j] = np.asarray(replacements)[inverse]
        words = {a: word_logs(corrected if a == 7 else columns[a]) for a in columns}
        sample_errors = []
        for a in columns:
            for index in sorted({0, len(words[a])//2, int(np.argmax(words[a]))}):
                exact = arb_word(columns[a][index], empty, marker if a == 7 else F(1))
                sample_errors.append(abs(float(exact)-words[a][index]))
        if max(sample_errors) > 1e-9:
            raise ArithmeticError('Chronological word controls disagree')
        prefix = (arb(SOURCE_COUNT).log()-734003*rational(F(9927, 10000)).log()-arb(255).log())/ln2
        root_masses = {}
        for a, values in words.items():
            unique, multiplicities = np.unique(values, return_counts=True)
            total = sum((int(n)*(rational(float(v))*ln2/18).exp()
                         for v, n in zip(unique, multiplicities)), arb(0))
            root_masses[a] = total/(arb(255)**(8*a))
            if a == 7:
                root_masses[a] /= rational(marker)**9
        per_arity = {str(a): float(prefix+18*mass.log()/ln2) for a, mass in root_masses.items()}
        mixed = float(prefix+18*sum(root_masses.values(), arb(0)).log()/ln2)
        target = next(row for row in receipt['controls'] if row['arity7_marker'] == '5/2')
        if max([abs(mixed-target['mixed_log2'])]+[abs(per_arity[a]-target['per_arity_log2'][a])
                                                  for a in per_arity]) > 1e-8:
            raise ArithmeticError('Saved selected result does not reproduce')
        result = dict(scope=expected, precision_bits=384, dependency_pins_checked=len(pins),
            regional_moment_slots_checked=76, moment_minimum_exact_dyadic_log_slack=str(min(moment_slacks)),
            moment_dyadic_roundoff_violations=sum(x < 0 for x in moment_slacks),
            old_minus_new_minimum_exact_dyadic_log_slack=str(min(monotonic_slacks)),
            regional_monotonicity_roundoff_violations=sum(x < 0 for x in monotonic_slacks),
            exact_scalar_receipt_matches=True, all_labelled_rows_checked={str(a): len(columns[a]) for a in columns},
            constant_functional_log2=[[float(x.log()/ln2) for x in row] for row in empty],
            minimum_activity_relative_residual_log2=float(min_relative.log()/ln2),
            sampled_arb_word_max_difference_bits=float(max(sample_errors)),
            arity8_branch_recombination_difference_bits=branch_error,
            prefix_log2=str(prefix), per_arity_log2=per_arity, mixed_log2=mixed,
            saved_mixed_difference_bits=mixed-target['mixed_log2'])
    check_time()
    return dict(result, receipt_sha256=RECEIPT_SHA, cache_sha256=receipt['cache_sha256'],
        exact_scalar_receipt_sha256=SCALAR_SHA, original_curve_receipt_sha256=OLD_SHA,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), elapsed_seconds=perf_counter()-started,
        regional_providers_run=0, column_contractions_run=0, directed=False, whole_code_certificate=False,
        input_interpretation='Stored binary64 logs are exact dyadics; source FLOAT columns are not certified uppers')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen recheck')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))
