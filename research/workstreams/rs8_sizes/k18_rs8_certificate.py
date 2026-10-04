"""Complete bounded K18 RS8/GF32 certificate with the exact retained s22 inner.

Authenticate every retained endpoint, refold all q1 and q2 support shells,
and transfer each q>=3 endpoint through a smaller uniform-input envelope.
Only standard-library exact arithmetic is needed; no census or benchmark runs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WORKSTREAMS = HERE.parent
GROUPS = 2048
MAP = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
REFERENCE = ROOT / 'tmp/k18-margin60-fractional-whole-p256-v2.json'
REFERENCE_SHA = '44fcdb0389d5cc7be870bdfcbeb0ec5500193706b7e690313054971d052a01d0'
Q1 = ROOT / 'tmp/rs16-s22-k262144-whole-fresh-p256-q1.json'
Q1_SHA = 'bc55924c4c4708b293a9bb1e1baacb81344e2a7fd7a25a0546ad85761c3f880c'
Q2 = ROOT / 'tmp/rs16-s22-k262144-whole-fresh-p256-q2.json'
Q2_SHA = '5e411ebefb9c5c8d48c22e20fc070fa98a4cfa01213231e7cc75045aeaed37d9'
HELPER = WORKSTREAMS / 'k20_two_ms/field128_transfer.py'
HELPER_SHA = '45d1fd4dbc4ecec20ea8ffac2fff26038fe044160582d6b0aa746e5ccc6b8153'
AUDIT = WORKSTREAMS / 'k18_margin60/verify_fractional.py'
AUDIT_SHA = '303cd2a0918637e44c854fd82792eae20c9dc903033274724a4bb5a0583919b0'
OUTER = WORKSTREAMS / 'k16_design/rs_outer.py'
OUTER_SHA = 'f0fda40a61536fdb672ffd474d45b9f6f6c58015c8cb4114366096d7749d7d2a'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pinned_module(name, path, expected):
    require(digest(path) == expected, f'frozen source changed: {path}')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


helpers = pinned_module('rs8_sizes_frozen_arithmetic', HELPER, HELPER_SHA)


def encode_exact(value):
    """Store an exact nonnegative dyadic without decimal digit limits."""
    value = F(value)
    require(value >= 0 and value.denominator & (value.denominator - 1) == 0,
            'nonnegative dyadic required')
    return dict(numerator_hex=hex(value.numerator),
                exponent=-(value.denominator.bit_length() - 1))


def probabilities(record, occupancy):
    """Decode the complete support interface, including zero-mass old shells."""
    if occupancy == 1:
        raw = record['support_probability_uppers']
        require(isinstance(raw, list) and len(raw) == 65, 'all 65 q1 supports required')
        values = list(map(helpers.dyadic, raw))
        flattened = values
    elif occupancy == 2:
        raw = record['support_pair_probability_uppers']
        require(isinstance(raw, list) and len(raw) == 65
                and all(isinstance(row, list) and len(row) == v + 1
                        for v, row in enumerate(raw)), 'all 2145 triangular q2 supports required')
        values = [list(map(helpers.dyadic, row)) for row in raw]
        flattened = [value for row in values for value in row]
    else:
        raise ValueError('only q1/q2 support endpoints exist')
    require(all(0 <= value <= 1 for value in flattened), 'support probabilities must lie in [0,1]')
    return values


def fold_one(counts, bounds, groups):
    require(len(counts) == len(bounds) == 65, 'complete q1 shell interface required')
    return groups * sum((count * value for count, value in zip(counts, bounds)), F())


def fold_two(counts, bounds, groups, *, low_support_only=False):
    """Unordered active groups, ordered support assignments within each pair."""
    require(len(counts) == len(bounds) == 65
            and all(len(row) == v + 1 for v, row in enumerate(bounds)),
            'complete triangular q2 shell interface required')
    result = F()
    for v, row in enumerate(bounds):
        for u, value in enumerate(row):
            if low_support_only and not (5 <= u <= 8 or 5 <= v <= 8):
                continue
            result += counts[v] * counts[u] * value * (1 if v == u else 2)
    return comb(groups, 2) * result


def validate_sparse(record, occupancy, old, outer):
    expected = dict(schema=f'finite-packet-rs-state-q{occupancy}-1',
        K=262144, N=524288, groups=GROUPS, regions=64, group_dimension=128,
        group_output_bits=256, distance='1/10', threshold=52428,
        physical_t=64, state_bits=22, physical_steps=8192, macro_t=128,
        macro_steps=4096, macros_per_region=64, physical_steps_per_macro=2,
        inner_distribution='uniform_gl', inner_group='GL(22,2)',
        physical_updates_independent=True, zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        occupancy_covered=[occupancy], all_supports_covered=True,
        all_two_group_support_pairs_covered=occupancy == 2, precision=256,
        count_kind='exact_expected_shells',
        regional_placement='exact without-replacement average of local upper envelopes')
    require(all(record.get(key) == value for key, value in expected.items()),
            f'q{occupancy} geometry or support scope mismatch')
    require(record['geometry'] == dict(group_count=GROUPS, regions=64,
            group_dimension=128, macro_windows=32, packet_bits=4), 'sparse packet geometry mismatch')
    require(record['map_record'] == old['map_record']
            and record['map_record']['map_sha256'] == MAP, 'sparse exact map mismatch')
    require(record['count_premises'] == dict(
        outer='four parallel GF16 RS[16,8] words',
        label_mixing='independent uniform GL16 per aligned four-nibble symbol',
        independent_setups_between_groups=True, exact_expected_shell_counts=True,
        independent_uniform_group_packet_shuffles=True,
        independent_uniform_regional_shuffles=True), 'sparse routing premises mismatch')
    if occupancy == 2:
        require(record['support_pair_storage'] ==
                'row v contains u=0..v; off-diagonal count products doubled', 'q2 multiplicity mismatch')
    old_counts = outer.expected_group_support_counts(16, 8, 4, 4)
    require(hashlib.sha256(json.dumps(list(map(str, old_counts))).encode()).hexdigest()
            == record['count_sha256'], 'old shell identity mismatch')
    values = probabilities(record, occupancy)
    reconstructed = (fold_one(old_counts, values, GROUPS) if occupancy == 1
                     else fold_two(old_counts, values, GROUPS))
    require(reconstructed <= helpers.dyadic(record[f'q{occupancy}_upper']),
            'old sparse fold exceeds authenticated endpoint')
    return values


def calculate(counts, one, two, selected):
    require(len(counts) == 65 and sum(counts) == (1 << 128) - 1
            and [v for v, count in enumerate(counts) if count] == list(range(5, 65)),
            'complete RS8 nonzero expected measure required')
    require(set(selected) == {str(q) for q in range(1, GROUPS + 1)},
            'retained coverage must include every q1..2048')
    beta_new = F(1 << 256, ((1 << 32) - 1)**4)
    beta_old = F(1 << 256, 65535**8)
    require(beta_new < beta_old, 'new uniform-input envelope must improve')
    q1, q2 = fold_one(counts, one, GROUPS), fold_two(counts, two, GROUPS)
    # Every saved tail endpoint keeps its full occupancy multiplicity.
    occupancy = dict(selected)
    occupancy['1'], occupancy['2'] = helpers.ceil_dyadic(q1), helpers.ceil_dyadic(q2)
    total = sum(map(helpers.dyadic, occupancy.values()), F())
    tail = sum((helpers.dyadic(pair) for q, pair in occupancy.items() if int(q) >= 3), F())
    upper = helpers.ceil_dyadic(total)
    sparse_q1 = GROUPS * sum((counts[v] * one[v] for v in range(5, 9)), F())
    sparse_q2 = fold_two(counts, two, GROUPS, low_support_only=True)
    return dict(exact_expected_shells=list(map(str, counts)),
        beta_new=str(beta_new), beta_old=str(beta_old),
        occupancy_uppers=occupancy, exact_union=encode_exact(total),
        exact_tail_union=encode_exact(tail), union_upper=upper,
        component_uppers=dict(q1=occupancy['1'], q2=occupancy['2'], tail=helpers.ceil_dyadic(tail)),
        component_margin_bits_display=dict(q1=helpers.margin(q1), q2=helpers.margin(q2),
                                           tail=helpers.margin(tail)),
        low_support_accounting=dict(packet_supports=[5, 6, 7, 8],
            q1_charge_upper=helpers.ceil_dyadic(sparse_q1),
            q2_charge_upper=helpers.ceil_dyadic(sparse_q2),
            q1_each_support_upper={str(v):helpers.ceil_dyadic(GROUPS * counts[v] * one[v])
                                   for v in range(5, 9)},
            q2_definition='at least one active group has packet support five through eight',
            included_in_full_folds=True, global_exclusion_used=False,
            margin_bits_display=dict(q1=helpers.margin(sparse_q1), q2=helpers.margin(sparse_q2))),
        margin_bits_display=helpers.margin(helpers.dyadic(upper)),
        target_bits=50, target_met=helpers.dyadic(upper) < F(1, 1 << 50),
        target_58_met=helpers.dyadic(upper) < F(1, 1 << 58),
        target_60_met=helpers.dyadic(upper) < F(1, 1 << 60))


def run():
    for path, expected in ((REFERENCE, REFERENCE_SHA), (Q1, Q1_SHA), (Q2, Q2_SHA)):
        require(digest(path) == expected, f'retained receipt changed: {path}')
    audit = pinned_module('rs8_sizes_k18_audit', AUDIT, AUDIT_SHA)
    authentication = audit.verify(REFERENCE)
    reference = json.loads(REFERENCE.read_text(encoding='utf-8'))
    old_path = Path(reference['retained_receipt'])
    old = json.loads(old_path.read_text(encoding='utf-8'))
    outer = pinned_module('rs8_sizes_outer', OUTER, OUTER_SHA)
    require(reference['source_sha256'].get(str(OUTER)) == OUTER_SHA,
            'reference does not authenticate exact RS enumerator')
    components = {Path(item['path']).resolve():item for item in old['components']}
    require(Q1.resolve() in components and Q2.resolve() in components,
            'sparse receipts must belong to authenticated original whole receipt')
    sparse1, sparse2 = (json.loads(path.read_text(encoding='utf-8')) for path in (Q1, Q2))
    one = validate_sparse(sparse1, 1, old, outer)
    two = validate_sparse(sparse2, 2, old, outer)
    counts = outer.expected_group_support_counts(8, 4, 4, 8)
    result = calculate(counts, one, two, reference['selected_occupancy_bounds'])
    require(audit.exact_union(result['exact_tail_union']) == audit.exact_union(reference['exact_tail_union']),
            'authenticated full q>=3 tail changed')
    pins = dict(reference['source_sha256'])
    for path in (REFERENCE, Q1, Q2, HELPER, AUDIT, OUTER, Path(__file__).resolve(),
                 HERE / 'K18_CERTIFICATE.md', HERE / 'test_k18_rs8_certificate.py'):
        actual = digest(path)
        require(str(path) not in pins or pins[str(path)] == actual, 'conflicting source pin')
        pins[str(path)] = actual
    for path, expected in pins.items():
        require(digest(path) == expected, f'source changed during calculation: {path}')
    result.update(schema='rs8-gf32-k18-s22-complete-1', whole_code_certificate=result['target_met'],
        all_occupancies_covered=True, K=262144, N=524288, distance='1/10',
        bad_cutoff=52428, minimum_distance=52429, groups=GROUPS,
        group_dimension=128, group_output_bits=256, packet_bits=4, regions=64,
        slots_per_region=GROUPS, physical_t=64, state_bits=22,
        outer='eight parallel GF16 RS[8,4] words; eight aligned 32-bit columns',
        randomizer='independent uniform nonzero GF(2^32) multipliers or their binary adjoints',
        route='independent uniform 64-packet shuffle per group and group permutation per region',
        inner_distribution='independent uniform GL(22,2) per physical step',
        map_sha256=MAP, map_record=old['map_record'], continuous_state=True,
        zero_initial_state=True, final_flush=False,
        reference=dict(path=str(REFERENCE), sha256=REFERENCE_SHA),
        sparse_references=[dict(path=str(Q1), sha256=Q1_SHA), dict(path=str(Q2), sha256=Q2_SHA)],
        reference_audit=authentication, source_sha256=pins,
        q1_supports_evaluated=list(range(65)), q2_triangular_pairs_evaluated=2145,
        q1_group_multiplier=GROUPS, q2_group_multiplier=comb(GROUPS, 2),
        q2_off_diagonal_multiplier=2, q2_diagonal_multiplier=1,
        inherited_occupancies=list(range(3, GROUPS + 1)),
        local_operators_recomputed=False, numerical_endpoints_recomputed=False,
        sparse_folds_recomputed=True,
        proof_method='complete exact q1/q2 folds; no-worse uniform-envelope transfer for every q>=3',
        arithmetic='exact rational shells and folds; outward integer dyadics; exact occupancy sum',
        scope='Ideal independent-setup probability that some nonzero message has output weight at most 52428. '
              'Not a distance test of a concrete seed. Local operator endpoints are authenticated, not regenerated.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'fresh output required')
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in
          ('margin_bits_display', 'component_margin_bits_display', 'target_met',
           'target_58_met', 'target_60_met')}, indent=2))
