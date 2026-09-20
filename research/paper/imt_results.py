"""Pinned selected finite IMT evidence; no search, replay, or benchmark implicit.

The pins identify accepted local receipts, not bundled data. Authentication and
exact union arithmetic do not independently prove the producer's inequalities.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'workstreams/inner_design'
MIGRATION = BASE / 'finite_migration'
HALF_PINS = {
    16: ('finite_migration/FULL_M16_VERIFIED.json', '03a70bbc631806205faff70a8402248225e3ae7bc6119e9fd4c26e131931f8f1'),
    18: ('finite_migration/FULL_M18_VERIFIED.json', '46fef112dbd9d6effcf402679fc17ba02d02569406bf610deafc8c0071b394ff'),
    20: ('asymmetric/bch256/weight5/FULL_M20_VERIFIED.json', 'f5cc1aa3ae775768e25685f94433fa4b51cb70b9dba775b2674b617d12e89a42'),
    22: ('finite_migration/FULL_M22_VERIFIED.json', '18d270c7d4f3b3470629c4fe1728f2e326a822e5628969c25724ca9df018e1d5'),
    24: ('finite_migration/FULL_M24_VERIFIED.json', '438d85071567063269f19fa9b09e422802d37a3f870882b6ef998e07e9f13108'),
}
TIMING_PIN = ('finite_migration/HALF_LADDER_PERFORMANCE_VERIFIED.json',
              'fe63621127c5f0821da65409ce8b32ac0e4ffc50403907a69099eb65a2e3ae13')
QUARTER_PIN = ('finite_migration/QUARTER_DEPLOYMENT_VERIFIED.json',
               'c399584893c6def0d0fa8e6c5725d2a03edaa4fddf6920d917fbc78e91995084')
QUARTER_PROOF_PIN = ('asymmetric/ASYMMETRIC_MARGIN_CERTIFICATE.json',
                    '1d136f7a74cc093a803a925544fb43aa0f5e733df5d3ea6c1c95f2e71dabfd3c')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fraction(value):
    if 'mantissa' in value:
        return int(value['mantissa']) * Fraction(2) ** int(value['exponent'])
    return Fraction(int(value['numerator']), int(value['denominator']))


def margin(value):
    require(value > 0, 'Expected positive probability bound')
    return math.log2(value.denominator) - math.log2(value.numerator)


def exact_union(parts, total, target, *, rounded=False):
    require(parts, 'Empty union')
    if 'mantissa' in total and all('mantissa' in p for p in parts):
        exponent = min([int(p['exponent']) for p in parts] + [int(total['exponent']), -target])
        require(all(int(p['mantissa']) >= 0 for p in parts), 'Negative union component')
        summed = sum(int(p['mantissa']) << (int(p['exponent']) - exponent) for p in parts)
        ceiling = int(total['mantissa']) << (int(total['exponent']) - exponent)
        require(summed <= ceiling if rounded else summed == ceiling, 'Union arithmetic differs')
        require(0 < ceiling < 1 << (-target-exponent), 'Target is not established')
        return margin(fraction(total))
    if all('mantissa' in p for p in parts):
        # Dense receipts can have million-bit exponents. Align dyadics once;
        # repeated Fraction addition would repeatedly normalize huge integers.
        exponent = min(int(p['exponent']) for p in parts)
        require(all(int(p['mantissa']) >= 0 for p in parts), 'Negative union component')
        numerator = sum(int(p['mantissa']) << (int(p['exponent']) - exponent) for p in parts)
        summed = numerator * Fraction(2) ** exponent
    else:
        values = [fraction(p) for p in parts]
        require(all(v >= 0 for v in values), 'Invalid union components')
        summed = sum(values, Fraction())
    bound = fraction(total)
    require(summed <= bound if rounded else summed == bound, 'Union arithmetic differs')
    require(0 < bound < Fraction(2) ** -target, 'Target is not established')
    return margin(bound)


def authenticate(pin, seen):
    name, expected = pin
    path = BASE / name
    require(path.is_file(), f'Missing IMT evidence: {path}; see finite_migration/README.md')
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, f'Receipt changed: {name}')
    data = json.loads(raw)
    for relative, digest in data['source_sha256'].items():
        target = (ROOT / relative).resolve()
        require(target.is_relative_to(ROOT), 'Dependency path escapes repository')
        if target in seen:
            require(seen[target] == digest, f'Conflicting source binding: {relative}')
            continue
        require(target.is_file(), f'Missing pinned dependency: {target}')
        with target.open('rb') as stream:
            actual = hashlib.file_digest(stream, 'sha256').hexdigest()
        require(actual == digest, f'Changed pinned dependency: {relative}')
        seen[target] = digest
    return data


def load():
    seen, half = {}, {}
    for m, pin in HALF_PINS.items():
        record = authenticate(pin, seen)
        instance = record['instance']
        require(record['full_distance_proved'] and record['target_bits'] == 40,
                'Expected complete finite certificate')
        require(record['covered_occupancies'] == [1, 2 ** (m-7)], 'Wrong coverage')
        require((instance['message_bits'], instance['output_bits'], instance['cutoff'])
                == (2**m, 2**(m+1), 2**(m+1)//10), 'Wrong half-rate geometry')
        require((instance['outer_length'], instance['outer_dimension'], instance['distance'])
                == (256, 128, '1/10'), 'Wrong outer or distance')
        inner = instance['inner']
        require((inner['t'], inner['s'], inner['transvection_rounds'], inner['feedback_name'])
                == (128, 19, 1, 'weight5_seed0'), 'Wrong finite IMT map family')
        if half:
            require(inner == half[16]['instance']['inner'], 'Ladder changes the inner')
        parts = record.get('component_upper')
        if parts is None:
            parts = [record[k] for k in ('q1_upper', 'sparse_upper', 'dense_upper')]
        bits = exact_union(parts, record['union_upper'], 40)
        require(abs(bits-record['margin_bits']) < 1e-8, 'Margin differs from exact union')
        half[m] = record

    timing = authenticate(TIMING_PIN, seen)
    require(timing['status'] == 'VERIFIED_IMT_HALF_RATE_LADDER_TIMINGS', 'Wrong timing binding')
    require(timing['full_certificate_exponents'] == [16, 18, 20], 'Missing timing certificate')
    require([r['message_exponent'] for r in timing['cells']] == [16, 18, 20], 'Wrong timing cells')
    for cell in timing['cells']:
        require(cell['distance_certificate_available'], 'Uncertified timing')
        require(cell['certified_margin_bits'] == half[cell['message_exponent']]['margin_bits'],
                'Timing/proof mismatch')
        for row in cell['summaries'].values():
            require(len(row['process_medians_ms']) == 3 and
                    statistics.median(row['process_medians_ms']) == row['median_ms'],
                    'Wrong process-median summary')

    quarter = authenticate(QUARTER_PIN, seen)
    proof = authenticate(QUARTER_PROOF_PIN, seen)
    require(quarter['status'] == 'VERIFIED_QUARTER_IMT_DEPLOYMENT_BINDING', 'Wrong quarter binding')
    require(proof['inner']['feedback_name'] == 'greedy3_2' and
            (proof['message_bits'], proof['output_bits']) == (2**20, 2**22), 'Wrong quarter instance')
    require(len(proof['results']) == 2, 'Expected two quarter thresholds')
    for row, delta, target in zip(proof['results'], ('33/200', '19/100'), (40, 30), strict=True):
        require(row['distance_target'] == delta and row['target_bits'] == target, 'Wrong cutoff')
        require(row['bad_weight'] == int(Fraction(delta) * 2**22), 'Wrong integer cutoff')
        parts = [row['q1_upper']] + row['q2_terms'] + row['adaptive_terms'] + row['dense_terms']
        bits = exact_union(parts, row['union_upper'], target, rounded=True)
        require(abs(bits-row['margin_bits_diagnostic']) < 1e-8, 'Wrong quarter margin')
    selected = [r for r in quarter['performance'] if r['distance_certificate_available']]
    require(len(selected) == 1 and selected[0]['message_exponent'] == 20, 'Wrong quarter timing scope')
    qtime = selected[0]['on']
    require(len(qtime['medians_ms']) == 3 and statistics.median(qtime['medians_ms']) == qtime['median_ms'],
            'Wrong quarter process median')
    return dict(half=half, timing=timing, quarter=quarter, quarter_proof=proof,
                quarter_timing=qtime, authenticated_files=len(seen))
