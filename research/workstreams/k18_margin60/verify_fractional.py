"""Independent authentication and exact-union audit of the K18 fractional receipt.

This verifier does not import the producer or its arithmetic helpers. It
authenticates the retained components and the fresh physical export, checks
every old/new selection using rational arithmetic, and checks the exact
2048-occupancy union. It does not regenerate the state census or independently
recompute the fractional numerical endpoints.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isclose, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MAP_SHA256 = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
OLD_SHA256 = '20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b'
LOCAL_SHA256 = 'fb5d68b115c8010726a7ff4cbf47da199e5d9cf217a411e1a970327d2a28d478'
DEFAULT_RECEIPT = ROOT / 'tmp/k18-margin60-fractional-whole-p256-v2.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def dyadic(pair):
    require(isinstance(pair, list) and len(pair) == 2, 'dyadic pair required')
    numerator, exponent = pair
    require(type(numerator) is int and numerator >= 0 and type(exponent) is int,
            'nonnegative integer dyadic required')
    return Fraction(numerator << exponent) if exponent >= 0 else Fraction(numerator, 1 << -exponent)


def exact_union(value):
    require(isinstance(value, dict) and set(value) == {'numerator_hex', 'exponent'},
            'hexadecimal exact-union record required')
    encoded, exponent = value['numerator_hex'], value['exponent']
    require(type(encoded) is str and encoded.startswith('0x') and type(exponent) is int,
            'hexadecimal numerator and integer exponent required')
    numerator = int(encoded, 16)
    require(numerator > 0 and hex(numerator) == encoded, 'positive canonical hexadecimal numerator required')
    return dyadic([numerator, exponent])


def upward_round(value, precision):
    """Round a positive dyadic upward to at most precision significant bits."""
    require(value > 0 and type(precision) is int and precision > 0, 'positive dyadic and precision required')
    require(value.denominator & (value.denominator - 1) == 0, 'power-of-two denominator required')
    excess = value.numerator.bit_length() - precision
    if excess <= 0:
        return value
    quantum = 1 << excess
    return Fraction(((value.numerator + quantum - 1) // quantum) * quantum, value.denominator)


def margin(value):
    require(value > 0, 'positive union required')
    return log2(value.denominator) - log2(value.numerator)


class SourceAudit:
    def __init__(self):
        self.observed = {}

    def digest(self, path):
        path = str(Path(path).resolve())
        if path not in self.observed:
            self.observed[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        return self.observed[path]

    def pins(self, pins):
        require(isinstance(pins, dict) and pins, 'nonempty source pins required')
        for path, expected in pins.items():
            require(type(expected) is str and len(expected) == 64, 'SHA256 source pin required')
            require(self.digest(path) == expected, f'source pin mismatch: {path}')


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def audit_old(path, sources):
    require(sources.digest(path) == OLD_SHA256, 'retained whole receipt SHA256 mismatch')
    old = load_json(path)
    expected = dict(schema='rs16-state-ladder-complete-replay-1', K=262144, N=524288,
        state_bits=22, groups=2048, regions=64, threshold=52428, physical_t=64,
        macro_t=128, macro_steps=4096, macros_per_region=64, physical_steps_per_macro=2,
        zero_initial_state=True, final_flush=False, whole_code_certificate=True,
        all_occupancies_covered=True, state_continuity='retained_between_every_physical_step_and_region')
    require(all(old.get(key) == value for key, value in expected.items()), 'retained geometry or scope mismatch')
    require(old['map_record']['map_sha256'] == MAP_SHA256, 'retained map mismatch')
    sources.pins(old['source_sha256'])
    require(set(old['occupancy_uppers']) == {str(q) for q in range(1, 2049)}, 'retained coverage mismatch')
    components = []
    require(len(old['components']) == 4, 'four retained components required')
    for item in old['components']:
        component = load_json(item['path'])
        canonical = json.dumps(component, sort_keys=True, separators=(',', ':')).encode()
        require(hashlib.sha256(canonical).hexdigest() == item['sha256'], 'retained component hash mismatch')
        sources.pins(component['source_sha256'])
        require(component['map_record']['map_sha256'] == MAP_SHA256, 'component map mismatch')
        components.append(component)
    for q, pair in old['occupancy_uppers'].items():
        dyadic(pair)
        require(any(component.get('occupancy_uppers', {}).get(q) == pair
                    or (q == '1' and component.get('q1_upper') == pair)
                    or (q == '2' and component.get('q2_upper') == pair)
                    for component in components), f'retained q{q} lacks an authenticated component')
    total = sum((dyadic(pair) for pair in old['occupancy_uppers'].values()), Fraction())
    require(upward_round(total, old['precision']) == dyadic(old['union_upper']),
            'retained rounded union differs from exact occupancy sum')
    return old


def audit_local(path, sources, old):
    require(sources.digest(path) == LOCAL_SHA256, 'fresh physical export SHA256 mismatch')
    local = load_json(path)
    require(local.get('schema') == 'k18-certified22-outward-physical-export-1', 'physical export schema mismatch')
    require(local.get('precision') == 192 and local.get('whole_code_certificate') is False,
            'physical export arithmetic scope mismatch')
    require(local['map_record'] == old['map_record'], 'physical export and retained map records differ')
    require(Fraction(local['tilt']) == Fraction(17, 200)
            and Fraction(local['local_selection_activity']) == Fraction(1, 2),
            'physical export witness mismatch')
    sources.pins(local['source_sha256'])
    metadata = local['physical_metadata']
    expected = dict(bits=22, windows=16, distance=16, map_sha256=MAP_SHA256,
        distribution='uniform_gl', birth_density='capped', birth_degree=8,
        birth_class_levels=[16, 24, 28, 32, 36, 40, 48, 64],
        map_images={'length': 4194304}, columns={'length': 64},
        records={'shape': [25, 5], 'dtype': 'int64'},
        birth_class_census={'shape': [8, 25], 'dtype': 'int64'})
    require(all(metadata.get(key) == value for key, value in expected.items()), 'physical metadata mismatch')
    raw = local['matrices']
    require(len(raw) == 17 and all(len(matrix) == 11 and all(len(row) == 11 for row in matrix)
                                  for matrix in raw), 'physical family dimensions mismatch')
    for matrix in raw:
        for row in matrix:
            for pair in row:
                dyadic(pair)
    require(dyadic(raw[0][0][0]) == 1 and all(dyadic(pair) == 0 for pair in raw[0][0][1:]),
            'zero-input zero-state transition mismatch')
    require(len(local['birth_row_choices']) == 17
            and all(type(index) is int and index >= 0 for index in local['birth_row_choices']),
            'complete birth-row choices required')
    return local


def audit_record(record):
    """Check a decoded receipt; deliberately independent of producer imports."""
    sources = SourceAudit()
    expected = dict(schema='k18-routing-conditioned-fractional-whole-1', whole_code_certificate=True,
        all_occupancies_covered=True, K=262144, N=524288, distance='1/10', threshold=52428,
        minimum_distance=52429, map_sha256=MAP_SHA256, outer_group_dimension=128, packet_bits=4,
        state_bits=22, regions=64, slots_per_region=2048, physical_step_bits=64,
        physical_steps_per_macro=2, macro_slots=32, macros_per_region=64,
        continuous_state=True, zero_initial_state=True, final_flush=False,
        precision=256, exact_geometry=True, thinning_before_fractional_power=True,
        chronological_product_before_power=True, placement_weights_outside_power=True,
        group_union_outside_power=True, local_census_recomputed_by_this_replay=False,
        retained_receipt_sha256=OLD_SHA256, target_bits=60, target_met=True)
    require(all(record.get(key) == value for key, value in expected.items()), 'whole receipt geometry or scope mismatch')
    require(Fraction(record['alpha']) == Fraction(4, 5) and Fraction(record['tilt']) == Fraction(17, 200),
            'fractional witness mismatch')
    old = audit_old(record['retained_receipt'], sources)
    local = audit_local(record['local_receipt'], sources, old)
    pins = record['source_sha256']
    sources.pins(pins)
    for origin in (old['source_sha256'], local['source_sha256']):
        require(all(pins.get(path) == value for path, value in origin.items()), 'whole receipt omits inherited source pins')
    for path in (Path(record['local_receipt']).resolve(), Path(record['retained_receipt']).resolve(),
                 HERE / 'replay_fractional.py', HERE / 'SLACK_REVIEW.md',
                 HERE.parent / 'k16_design/packet_regional_power.py'):
        require(pins.get(str(path)) == sources.digest(path), f'required receipt or source pin missing: {path}')
    candidates = record['new_occupancy_bounds']
    chosen = record['selected_occupancy_bounds']
    require(set(candidates) == {str(q) for q in range(62, 144)}, 'new occupancy coverage must be q62..143')
    require(set(chosen) == {str(q) for q in range(1, 2049)}, 'whole occupancy coverage must be q1..2048')
    changed = []
    for q, retained in old['occupancy_uppers'].items():
        retained_value = dyadic(retained)
        candidate = candidates.get(q)
        if candidate is not None:
            candidate_value = dyadic(candidate)
            require(candidate_value > 0, f'positive outward endpoint required for q{q}')
        use_new = candidate is not None and candidate_value < retained_value
        expected_pair = candidate if use_new else retained
        require(chosen[q] == expected_pair, f'wrong exact old/new selection at q{q}')
        if use_new:
            changed.append(int(q))
    require(record['improved_occupancies'] == sorted(changed), 'improved occupancy list mismatch')
    total = sum((dyadic(pair) for pair in chosen.values()), Fraction())
    tail = sum((dyadic(pair) for q, pair in chosen.items() if int(q) >= 3), Fraction())
    require(exact_union(record['exact_union']) == total, 'whole exact union mismatch')
    require(exact_union(record['exact_tail_union']) == tail, 'tail exact union mismatch')
    require(0 < total < Fraction(1, 1 << 60), 'exact union does not meet the 60-bit target')
    require(isclose(record['margin_bits_display'], margin(total), rel_tol=0, abs_tol=1e-8)
            and isclose(record['tail_margin_bits_display'], margin(tail), rel_tol=0, abs_tol=1e-8),
            'displayed margin mismatch')
    # Re-read each source at completion, avoiding a stale cache across the audit.
    for path, expected in sources.observed.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, f'source changed during audit: {path}')
    return dict(schema='k18-fractional-independent-union-audit-1', verified=True,
        authenticated_file_count=len(sources.observed), authenticated_old_components=4,
        authenticated_occupancies=2048, new_occupancies=82, improved_occupancies=len(changed),
        exact_selections_match=True, exact_union_matches=True, exact_tail_matches=True,
        target_bits=60, target_met=True, margin_bits=margin(total), tail_margin_bits=margin(tail),
        retained_receipt_sha256=OLD_SHA256, physical_receipt_sha256=LOCAL_SHA256,
        local_census_recomputed=False, new_fractional_endpoints_independently_recomputed=False,
        scope='Independent hash, geometry, component provenance, exact selection, and exact-union audit. '
              'The state census and new fractional numerical endpoints are authenticated, not independently recomputed.')


def verify(path):
    raw = Path(path).read_bytes()
    result = audit_record(json.loads(raw))
    result['receipt_sha256'] = hashlib.sha256(raw).hexdigest()
    require(Path(path).read_bytes() == raw, 'receipt changed during audit')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        raise ValueError('fresh audit output required')
    result = verify(args.receipt)
    print(json.dumps(result, indent=2))
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')


if __name__ == '__main__':
    main()
