"""Separate count transfer from28 common pairs to all28^16 assignments.

No regional/operator computation is performed. The reviewed mixed resource
DP used only geometry-independent box and total-sum constraints. Its caps
therefore apply to every fixed assignment of sixteen constituent pairs.
The previous receipts remain unchanged and still count28 common-pair cases.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling


DIRECTORY = Path(__file__).resolve().parent
VERIFIED_NAME = 'hd_xor_mixed_mask_8_0_mu16_dual_verified.json'
VERIFIED_SHA256 = '1db1308fc80eb9787ffbcd0a6eb694113b363dc17ce2d1cb596ecc18bb1d1e3b'
CERTIFICATE_SHA256 = 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b'
CONSTRUCTION = 'xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half'
BASE_CEILING = F('-274659876277/2500000000')


def transfer():
    """Read pinned observed evidence and recompute only the positive mass multiplier."""
    data = (DIRECTORY/VERIFIED_NAME).read_bytes()
    if sha256(data).hexdigest() != VERIFIED_SHA256:
        raise ValueError('The pinned verified common-pair evidence changed')
    if sha256((DIRECTORY/'hd_xor_mixed_mask_certificate.py').read_bytes()).hexdigest() != CERTIFICATE_SHA256:
        raise ValueError('The reviewed geometry-independent operator source changed')
    source = json.loads(data)
    expected = dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], mu_range=[16, 16], fixed_row_sets=219,
        exact_partition_checked=True, all_per_class_ceilings_agree=True,
        aggregate_ceilings_agree=True, independent_sum_precision_bits=512,
        rounded_cell_union_log2_upper=str(BASE_CEILING), whole_code_certificate=False)
    for key, value in expected.items():
        if type(source.get(key)) is not type(value) or source[key] != value:
            raise ValueError(f'Wrong source evidence field: {key}')
    factor = 28**15
    with ctx.workprec(512):
        base = arb(BASE_CEILING.numerator)/BASE_CEILING.denominator
        # This multiplication is the derived bound. It uses neither the
        # unrounded aggregate nor a reinterpreted previous receipt.
        mass = factor*arb(2)**base
        result = rational_ceiling(mass.log()/arb(2).log())
        independent_log_form = rational_ceiling(base+15*arb(28).log()/arb(2).log())
        if result != independent_log_form:
            raise ArithmeticError('Independent mass and logarithmic evaluations disagree')
        penalty = rational_ceiling(15*arb(28).log()/arb(2).log())
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], active_groups=17,
        source_matching_counts=[8, 0], source_classes=1, field_order=256,
        forbidden_weight=734003, total_columns=32, fixed_dense_row_sets=219,
        sparse_support_floor=544, dense_nonempty_column_floor=17,
        previous_mu_range=[16, 16], derived_mu_range=[1, 16],
        previous_sparse_pair_assignments=28, derived_sparse_pair_assignments=28**16,
        exact_mass_multiplier=factor, multiplier_log2_upper=str(penalty),
        previous_rounded_cell_union_log2_upper=str(BASE_CEILING),
        failure_log2_upper=str(result), precision_bits=512, directed=True,
        arithmetic='positive-mass-multiplication-cross-checked-by-logarithmic-addition',
        reviewed_resource_constraints='d-k<=m_i<=s+d; sum(m_i)=7s+8d-6k',
        pairwise_geometry_constraints_used=False,
        constituent_geometry_fixed_across_columns=True,
        geometry_resampled_across_columns=False, setup_dependent_witness_selection=False,
        all_fixed_dense_sets_covered=True, all_sparse_pair_assignments_covered=True,
        reaches_40_bits=result < -40, whole_code_certificate=False,
        old_receipts_unchanged=True, operator_replays=0, fitted_new_witnesses=0,
        source_evidence_filename=VERIFIED_NAME, source_evidence_sha256=VERIFIED_SHA256,
        reviewed_certificate_source_sha256=CERTIFICATE_SHA256,
        evidence_trust='observed-directed-replays-and-reviewed-universal-resource-inequality')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a count-transfer record')
    result = transfer()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)
