"""Replay the expanded PARTIAL r18 cover; the whole-code goal remains open."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_lambda_renyi_a10_control as control
import hd_lambda_renyi_activity_verify as prior

DIRECTORY = Path(__file__).resolve().parent
PINS = {
    'hd_lambda_renyi_a10_control.py': '1aa06f252a2f68ff0090281d9b3fbc8d9e7b9bd233de2e913e6cf60703208911',
    'hd_lambda_renyi_activity_verify.py': 'd64f82392d72a19a4c0cbb89ae5e1bee8d6b4c3190a832cb48184c1b18aada15',
    'hd_lambda_renyi_a10_control_256.json': 'c6b56c7a25b1bf9ea88d38e6fb75b60bb6495d820ff40e2107aecda2944cfd6d',
    'hd_lambda_renyi_a10_control_384.json': '3f12920d9a93f1982ac4c07d9f942d0aabc2eb4b3aa584a585923639b9ff11e3',
}


def verify():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pinned proof input changed: '+name)
    control.dependencies()
    control.inputs()
    old = prior.verify()
    observed = []
    for precision in (256, 384):
        data = json.loads((DIRECTORY/f'hd_lambda_renyi_a10_control_{precision}.json').read_text())
        expected = dict(construction=control.activity.binding.SCOPE['construction'],
            active_groups=18, every_group_arity_range=[10, 16], total_arity_range=[180, 245],
            source_profiles=99, precision_bits=precision, directed=True, whole_code_certificate=False,
            reviewed_source_sha256=control.dependencies(), inventory_sha256=control.INPUTS,
            all14893_labelled_rowsets=True, minimum_active_columns_per_group=9,
            original_rowsets_fixed_across_columns=True, source_count_outside_theta=True,
            rowset_sum_inside_theta=True, projective_normalizations=1)
        for key, value in expected.items():
            if type(data.get(key)) is not type(value) or data[key] != value:
                raise ValueError('Wrong observed scope: '+key)
        if [c['moment_order'] for c in data['controls']] != [2, 3, 4]:
            raise ValueError('All three method controls required')
        for c in data['controls']:
            if [row['total_arity'] for row in c['totals']] != list(range(180, 246)):
                raise ValueError('Wrong method total partition')
        for index, selected in enumerate(data['selected_totals']):
            best = min(data['controls'], key=lambda c: F(c['totals'][index]['failure_log2_upper']))
            if (selected['total_arity'] != index+180
                    or selected['failure_log2_upper'] != best['totals'][index]['failure_log2_upper']
                    or selected['selected_moment_order'] != best['moment_order']):
                raise ValueError('Selected method does not match completed bounds')
        observed.append(data)
    if ([c['totals'] for c in observed[0]['controls']] != [c['totals'] for c in observed[1]['controls']]
            or observed[0]['selected_totals'] != observed[1]['selected_totals']):
        raise ValueError('Fresh precision replays disagree')
    selected = [row for row in observed[0]['selected_totals'] if row['total_arity'] <= 202]
    if [row['total_arity'] for row in selected] != list(range(180, 203)):
        raise ValueError('Only the passing complete lower-total range is permitted')
    local = control.A.rounded_union([F(row['failure_log2_upper']) for row in selected], precision=512)
    union = control.A.rounded_union([local, F(old['failure_log2_upper'])], precision=512)
    return dict(construction=old['construction'], active_groups=18, field_order=256,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        classes=old['classes']+[dict(every_group_arity_range=[10, 16],
            total_arity_range=[180, 202], failure_log2_upper=str(local))],
        new_class_failure_log2_upper=str(local), new_class_margin_bits_float=-float(local),
        failure_log2_upper=str(union), margin_bits_float=-float(union),
        all99_indexed_source_profiles=True, directed=True, whole_code_certificate=False,
        other_group_counts_covered=False, overlap_handled_by_probability_union=True,
        replay_precisions=[256, 384], sum_precision=512, paired_method_ceilings_agree=True,
        input_sha256=PINS, verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


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
    print(json.dumps({k: v for k, v in result.items() if k not in ('classes', 'input_sha256')}))
