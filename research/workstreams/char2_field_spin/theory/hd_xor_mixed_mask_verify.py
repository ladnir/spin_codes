"""Independent scope/assignment check and512-bit sum of fresh cell ceilings.

This checks receipts from the observed fresh directed replays. Hashes bind
their bytes and reviewed sources; hashes do not independently prove their
arithmetic. The only arithmetic recomputed here is the positive sum of all
219 rational logarithmic ceilings, in a new512-bit Arb precision context.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_mixed_mask_certificate import EXPECTED_SCOPE
from hd_xor_mixed_mask_witnesses import checked_classes, PROPOSAL_SHA256


DIRECTORY = Path(__file__).resolve().parent
CERTIFICATE_SOURCE_SHA256 = 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b'
WITNESS_SOURCE_SHA256 = '6b6c4ceb63a9624ccdcbbf7d5ab78beabd19caa3224dbb193eb68509c4e09a61'
SOURCE_COUNT = 15708305994709246428578365440


def checked_rows(receipt, precision):
    """Validate the exact fixed assignment; infer arity from each mask."""
    expected = {key: list(value) if isinstance(value, tuple) else value for key, value in EXPECTED_SCOPE.items()}
    expected.update(precision_bits=precision, directed=True, fitted_new_witnesses=0,
        frozen_proposal_sha256=PROPOSAL_SHA256, source_sha256=CERTIFICATE_SOURCE_SHA256,
        frozen_witness_source_sha256=WITNESS_SOURCE_SHA256, source_supports=SOURCE_COUNT,
        dense_source_colorings=17, two_row_geometries=28, fixed_row_sets=219,
        scoped_class_complete=True, complete_mixed_class=False, whole_code_certificate=False,
        dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
        dense_input_includes_zero=True, exact_uniform_dense_mask_mixture=True,
        masks_and_cancellations_mixed_before_potential=True, same_potential_both_dense_branches=True,
        dense_beta_formula='q^(32a)/(q-1)^(16a)', exact_positive_bivariate_tail=True,
        chernoff_dense_floor=False, product_density_comparison_used=False,
        projective_normalizations=1, endpoint_potential_factors=1)
    for key, value in expected.items():
        if type(receipt.get(key)) is not type(value) or receipt[key] != value:
            raise ValueError(f'Wrong directed receipt field: {key}')
    assignment = {mask: (label, str(z), str(logp))
                  for label, z, logp, masks in checked_classes() for mask in masks}
    entries = receipt.get('entries', [])
    if len(entries) != 219:
        raise ValueError('Exactly219 fixed row-set contributions are required')
    rows = {}
    for row in entries:
        mask = row.get('active_row_mask')
        if type(mask) is not int or mask not in assignment or mask in rows:
            raise ValueError('Unexpected or duplicated fixed dense row set')
        if type(row.get('active_constituents')) is not int or row['active_constituents'] != mask.bit_count():
            raise ValueError('Arity must agree with the fixed mask')
        if (row.get('witness_label'), row.get('mark'), row.get('log2_potential')) != assignment[mask]:
            raise ValueError('The per-class witness must match the frozen assignment exactly')
        if not isinstance(row.get('failure_log2_upper'), str):
            raise ValueError('An exact rational logarithmic ceiling is required')
        ceiling = F(row['failure_log2_upper'])
        if ceiling >= -40:
            raise ValueError('This frozen successful replay requires each cell below the40-bit threshold')
        rows[mask] = ceiling
    if set(rows) != set(assignment):
        raise ValueError('The219-class partition is incomplete')
    return rows


def rounded_union_512(rows):
    """No cached aggregate or FLOAT arithmetic participates in this sum."""
    with ctx.workprec(512):
        total = sum((arb(2)**(arb(value.numerator)/value.denominator)
                     for _, value in sorted(rows.items())), arb(0))
        return rational_ceiling(total.log()/arb(2).log())


def verify(directory=DIRECTORY):
    directory = Path(directory)
    for filename, pin in (('hd_xor_mixed_mask_certificate.py', CERTIFICATE_SOURCE_SHA256),
                          ('hd_xor_mixed_mask_witnesses.py', WITNESS_SOURCE_SHA256)):
        if sha256((DIRECTORY/filename).read_bytes()).hexdigest() != pin:
            raise ValueError('Reviewed source pin changed before dual-receipt verification')
    inputs, receipts, rows = [], [], []
    for precision in (256, 384):
        path = directory/f'hd_xor_mixed_mask_8_0_mu16_receipt{precision}.json'
        data = path.read_bytes()
        receipt = json.loads(data)
        rows.append(checked_rows(receipt, precision))
        receipts.append(receipt)
        inputs.append(dict(filename=path.name, sha256=sha256(data).hexdigest(),
                           precision_bits=precision, elapsed_seconds=receipt['elapsed_seconds']))
    if rows[0] != rows[1]:
        raise ValueError('The per-class rational ceilings differ between256 and384 bits')
    if receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('The raw aggregate rational ceilings differ between precisions')
    union = rounded_union_512(rows[0])
    direct = F(receipts[0]['failure_log2_upper'])
    if union < direct:
        raise ArithmeticError('Rounded-cell union unexpectedly lies below the raw aggregate ceiling')
    arities = [dict(active_constituents=active, fixed_row_sets=sum(mask.bit_count() == active for mask in rows[0]),
                   rounded_cell_union_log2_upper=str(rounded_union_512(
                       {mask: ceiling for mask, ceiling in rows[0].items() if mask.bit_count() == active})))
               for active in range(3, 9)]
    return dict(construction=EXPECTED_SCOPE['construction'], cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], mu_range=[16, 16], fixed_row_sets=219,
        exact_partition_checked=True, all_per_class_ceilings_agree=True,
        aggregate_ceilings_agree=True, receipt_precision_bits=[256, 384],
        independent_sum_precision_bits=512, raw_aggregate_log2_upper=str(direct),
        rounded_cell_union_log2_upper=str(union), arity_summaries=arities,
        selected_classes_reach_40_bits=union < -40, scoped_class_complete=True,
        complete_mixed_class=False, whole_code_certificate=False, inputs=inputs,
        receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes',
        certificate_source_sha256=CERTIFICATE_SOURCE_SHA256,
        frozen_witness_source_sha256=WITNESS_SOURCE_SHA256,
        frozen_proposal_sha256=PROPOSAL_SHA256)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a saved dual verification')
    result = verify()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)
