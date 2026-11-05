"""Verify the observed directed r18 dense-group partial cover.

This checks the pinned fresh replay receipts, not merely their names.
Hashes bind observed arithmetic; they do not substitute for that arithmetic
or its reviewed comparison proof. The cover is not the whole construction.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_lambda_renyi_mixed as mixed

DIRECTORY = Path(__file__).resolve().parent
SOURCE_PIN = '9dc5d2cfba096bab6e3fc3b2ceb39425f62450ea0364a53407b7eff3405c1723'
INPUTS = {
    'hd_lambda_renyi_a12_inventory.json': 'a5560281ebc57e0e1782aa5e40f653ac4547c9ec498d19edf16a5bbcb1ee56a9',
    'hd_lambda_renyi_orbit_inventory34.json': '227fd5c05f5946698d6fac8a703abe4657ec154f817f48352d0ce346f818e673',
    'hd_lambda_renyi_mixed_256.json': '534de633bc59f905463e5260253a0b947d7c4ca70747539692a3d32b78579401',
    'hd_lambda_renyi_mixed_384.json': 'bf5f2d5dcfd0d12393350779dab811b1054783f57d9f6447bde8feec21fdac74',
    'hd_lambda_dense_bridge_256.json': '6b93a7e7b194ca932118bdfca891efd1f6324149e99b18d580bebb33cb43597d',
    'hd_lambda_dense_bridge_384.json': '062138e3c1e22a9c8a1c837ca1e44ab70deb5b11160028961a00f3d6dec0bba4',
}


def verify():
    pins = mixed.dependencies()
    if pins['hd_lambda_renyi_mixed.py'] != SOURCE_PIN:
        raise ValueError('Reviewed mixed controller changed')
    data = {}
    for name, digest in INPUTS.items():
        raw = (DIRECTORY/name).read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Observed replay input changed: '+name)
        data[name] = json.loads(raw)
    direct = data['hd_lambda_renyi_a12_inventory.json']
    orbit = data['hd_lambda_renyi_orbit_inventory34.json']
    moments = mixed.checked_inventory(orbit)
    if {tuple(row['complement']): F(row['moment']) for row in direct['rows']} != {
            key: value[3] for key, value in moments.items() if len(key) == 4}:
        raise ValueError('Orbit moments disagree with independent full labelled enumeration')
    receipts = [data[f'hd_lambda_renyi_mixed_{precision}.json'] for precision in (256, 384)]
    for precision, receipt in zip((256, 384), receipts):
        expected = dict(precision_bits=precision, directed=True, whole_code_certificate=False,
            active_groups=18, every_group_arity_range=[12, 16], total_arity_range=[216, 246],
            original_rowsets_fixed_across_columns=True, scheduled_moment_products=True,
            source_profiles=99, fresh_regional_banks=1, reviewed_source_sha256=pins,
            inventory_sha256=INPUTS['hd_lambda_renyi_orbit_inventory34.json'])
        mixed.A.require_fields(receipt, expected)
    if receipts[0]['selected_totals'] != receipts[1]['selected_totals'] or [
            row['totals'] for row in receipts[0]['controls']] != [
            row['totals'] for row in receipts[1]['controls']]:
        raise ValueError('Fresh directed mixed replays disagree')
    low_rows = receipts[0]['selected_totals'][:30]
    if [row['total_arity'] for row in low_rows] != list(range(216, 246)):
        raise ValueError('Wrong mixed total-arity partition')
    low = mixed.A.rounded_union([F(row['failure_log2_upper']) for row in low_rows], precision=512)
    high_rows = [data[f'hd_lambda_dense_bridge_{p}.json']['controls'][0]['cutoff_ceilings'][30]
                 for p in (256, 384)]
    if high_rows[0] != high_rows[1] or high_rows[0]['total_arity_minimum'] != 246:
        raise ValueError('Fresh high-side cutoff replays disagree')
    high = F(high_rows[0]['failure_log2_upper'])
    union = mixed.A.rounded_union([low, high], precision=512)
    return dict(construction=mixed.binding.SCOPE['construction'], field_order=256,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        active_groups=18, every_group_arity_at_least12_covered=True,
        additional_class_all_total_arity_at_least246_covered=True,
        all99_indexed_source_profiles=True,
        mixed_low_class=dict(every_group_arity_range=[12, 16], total_arity_range=[216, 245],
            failure_log2_upper=str(low)),
        high_class=dict(total_arity_range=[246, 288], any_individual_arity=True,
            failure_log2_upper=str(high)),
        failure_log2_upper=str(union), margin_bits_float=-float(union),
        directed=True, whole_code_certificate=False,
        remaining_r18_class='total arity <=245 with at least one group of arity <=11',
        other_group_counts_covered=False, replay_precisions=[256, 384], sum_precision=512,
        all62_method_ceilings_and31_selected_ceilings_agree=True,
        all1820_third_moments_match_direct_labelled_enumeration=True,
        probability_union_does_not_assume_disjoint_bad_setup_events=True,
        input_sha256=INPUTS, reviewed_source_sha256=pins,
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = verify()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('input_sha256', 'reviewed_source_sha256')}))
