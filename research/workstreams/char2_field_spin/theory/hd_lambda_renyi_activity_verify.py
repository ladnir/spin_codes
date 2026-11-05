"""Verify the expanded, still partial r18 cover and paired diagnostics."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_lambda_renyi_cover_verify as old
import hd_lambda_renyi_activity as activity
import hd_lambda_renyi_marked_control as marked
import hd_lambda_renyi_interpolate as interpolate
import hd_lambda_renyi_limits as limits

DIRECTORY = Path(__file__).resolve().parent
INPUTS = {
    'hd_lambda_renyi_activity_a11_256.json': '77edcd1ea1f70b3557e249198912a267721b37dca82a3bba4d8abac1fe841464',
    'hd_lambda_renyi_activity_a11_384.json': 'af2b5fb312136bb37ebabc41c9523122cb8a1c8d30d8d7e8773fddbf ee97d766'.replace(' ', ''),
    'hd_lambda_renyi_marked_control_256.json': '5f2dd8a84852d7760f62a57b14479aa1e0b8457364cd4fdd5da17c1b09c97ad2',
    'hd_lambda_renyi_marked_control_384.json': 'ceb3c9ec647280e997cdd16aa2015c3b63bec9bdbd79c11733a9cb7236bea563',
    'hd_lambda_renyi_interpolate_256.json': 'a22ecfd2a976c3c8c6d53a1486597d484ce9a1ade6104b4141ed0f26f6066fa6',
    'hd_lambda_renyi_interpolate_384.json': '30fc85e78bbbce3cd314a362a5e7bee70dadb1a69a78605dc0af207d0c2c7e6c',
    'hd_lambda_renyi_limits_256.json': 'f3d433f3e55d369741f93e53370d57e1423ab5b12ff033328c168e817ddf71be',
    'hd_lambda_renyi_limits_384.json': '0a425a1acb6466880226c08abcd131ab910b4f76f3f9ed88a939a2b99644f31e',
}


def verify():
    previous = old.verify()
    families = {}
    for stem, module in (('activity_a11', activity), ('marked_control', marked),
                         ('interpolate', interpolate), ('limits', limits)):
        receipts = []
        for precision in (256, 384):
            name = f'hd_lambda_renyi_{stem}_{precision}.json'
            raw = (DIRECTORY/name).read_bytes()
            if sha256(raw).hexdigest() != INPUTS[name]:
                raise ValueError('Observed receipt changed: '+name)
            data = json.loads(raw)
            activity.A.require_fields(data, dict(precision_bits=precision, active_groups=18,
                construction=activity.binding.SCOPE['construction'], whole_code_certificate=False,
                source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
                reviewed_source_sha256=module.dependencies(), fresh_regional_banks=1))
            receipts.append(data)
        left, right = receipts
        projections = lambda d: [(c.get('totals'), c.get('homogeneous'), c.get('mixed')) for c in d['controls']]
        if projections(left) != projections(right) or left.get('selected_totals') != right.get('selected_totals'):
            raise ValueError('Fresh directed ceilings disagree: '+stem)
        families[stem] = left
    new = families['activity_a11']
    if (new['every_group_arity_range'] != [11, 16] or new['source_profiles'] != 99
            or new['minimum_active_columns_per_group'] != 9):
        raise ValueError('Wrong activity scope')
    lower = [r for r in new['selected_totals'] if r['total_arity'] <= 225]
    if [r['total_arity'] for r in lower] != list(range(198, 226)):
        raise ValueError('Wrong new class partition')
    low_bound = activity.A.rounded_union([F(r['failure_log2_upper']) for r in lower], precision=512)
    sparse = families['marked_control']
    if sparse['profile_count'] != 99 or sparse['every_group_arity_range'] != [1, 5]:
        raise ValueError('Wrong marked comparison scope')
    singleton = min(F(c['homogeneous'][0]['failure_log2_upper']) for c in sparse['controls'])
    if any(c['homogeneous'][0]['every_group_arity'] != 1 for c in sparse['controls']):
        raise ValueError('Wrong homogeneous arity')
    union = activity.A.rounded_union([F(previous['failure_log2_upper']), low_bound, singleton], precision=512)
    return dict(construction=activity.binding.SCOPE['construction'], active_groups=18,
        field_order=256, source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        classes=[dict(every_group_arity_range=[11, 16], total_arity_range=[198, 225], failure_log2_upper=str(low_bound)),
                 dict(every_group_arity_range=[12, 16], total_arity_range=[216, 245],
                      failure_log2_upper=previous['mixed_low_class']['failure_log2_upper']),
                 previous['high_class'],
                 dict(every_group_arity=1, failure_log2_upper=str(singleton))],
        failure_log2_upper=str(union), margin_bits_float=-float(union),
        all99_indexed_source_profiles=True, directed=True, whole_code_certificate=False,
        other_group_counts_covered=False, overlap_handled_by_probability_union=True,
        diagnostics=dict(sparse_all1to5_failure_log2_upper=sparse['selected_mixed']['failure_log2_upper'],
            ideal_method_not_failure_bound=True,
            ideal_fixed_method_log2_lower=families['limits']['ideal_all_totals_method_expression_log2_lower'],
            intermediate_orders_do_not_extend_cutoff225=True),
        replay_precisions=[256, 384], sum_precision=512, paired_method_ceilings_agree=True,
        input_sha256=INPUTS, prior_cover_verifier_sha256=sha256((DIRECTORY/'hd_lambda_renyi_cover_verify.py').read_bytes()).hexdigest(),
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
    print(json.dumps({k: v for k, v in result.items() if k != 'input_sha256'}))
