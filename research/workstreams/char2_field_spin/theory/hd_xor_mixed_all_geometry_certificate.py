"""Fresh directed replay of one frozen all-geometry mixed source class.

The five-bank219-mask assignment is pinned independently of setup. Existing
reviewed primitives are reused without changing their sources or globals.
The universal resource DP's mu=16 parameter imposes no multiplicity bound;
the one-time28-to-28^16 count transfer is explicit in every receipt.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_mask_certificate import (mixed_vectors_arb,
    fixed_branch_caps_arb, column_coefficients_arb, positive_powers_arb,
    bivariate_tail_arb, _rational, _upper_max)
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import arb_pair_regions
from hd_xor_paired_joint_gate import ordered_tables


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_xor_mixed_all_geometry_8_0_refine64.json'
PROPOSAL_SHA256 = '44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6'
PRIMITIVES_SHA256 = 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b'
SOURCE_COUNT = 15708305994709246428578365440
ALL_MASKS = tuple(mask for mask in range(256) if mask.bit_count() >= 3)
EXPECTED_SCOPE = dict(
    construction='xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half',
    cohort_counts=[0, 16, 1], source_matching_counts=[8, 0], mu_range=[1, 16],
    source_classes=1, active_groups=17, sparse_active_rows=32, field_order=256,
    forbidden_weight=734003, total_columns=32, sparse_support_floor=544,
    dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
    pair_swaps=True, alpha=1, source_supports=SOURCE_COUNT,
    dense_source_colorings=17, base_sparse_pair_assignments=28,
    sparse_pair_assignments=28**16, exact_mass_multiplier=28**15,
    fixed_row_sets=219, scoped_class_complete=True, complete_mixed_class=False,
    whole_code_certificate=False, dense_row_set_fixed_across_columns=True,
    sparse_geometry_fixed_across_columns=True, geometry_resampled_across_columns=False,
    setup_dependent_selection=False, dense_input_includes_zero=True,
    exact_uniform_dense_mask_mixture=True, masks_and_cancellations_mixed_before_potential=True,
    same_potential_both_dense_branches=True, exact_positive_bivariate_tail=True,
    dense_beta_formula='q^(32a)/(q-1)^(16a)', projective_normalizations=1,
    endpoint_potential_factors=1, pairwise_geometry_constraints_used=False,
    directed=True, fitted_new_witnesses=0)


def require_fields(data, expected):
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError(f'Wrong frozen scope field: {key}')


def validate_scope(scope=None):
    require_fields(EXPECTED_SCOPE if scope is None else scope, EXPECTED_SCOPE)
    if (Q, REGION_LENGTH, LIMIT) != (256, 8192, 734003):
        raise ValueError('Imported construction constants changed')
    if sha256((DIRECTORY/'hd_xor_mixed_mask_certificate.py').read_bytes()).hexdigest() != PRIMITIVES_SHA256:
        raise ValueError('Reviewed primitive source changed')
    return dict(EXPECTED_SCOPE)


def checked_classes():
    data = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(data).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen all-geometry proposal changed')
    proposal = json.loads(data)['combined']
    require_fields(proposal, {key: EXPECTED_SCOPE[key] for key in (
        'construction', 'cohort_counts', 'source_matching_counts', 'mu_range',
        'field_order', 'forbidden_weight', 'total_columns', 'sparse_pair_assignments',
        'base_sparse_pair_assignments', 'exact_mass_multiplier', 'fixed_row_sets',
        'sparse_support_floor', 'dense_nonempty_column_floor', 'pair_swaps', 'alpha',
        'dense_row_set_fixed_across_columns', 'sparse_geometry_fixed_across_columns',
        'geometry_resampled_across_columns', 'setup_dependent_selection')})
    require_fields(proposal, dict(directed=False, proposal_requires_fresh_directed_replay=True))
    assigned = {}
    seen = set()
    for row in proposal['entries']:
        mask = row.get('active_row_mask')
        if type(mask) is not int or mask not in ALL_MASKS or mask in seen:
            raise ValueError('Each of219 fixed dense masks must occur exactly once')
        if type(row.get('active_constituents')) is not int or row['active_constituents'] != mask.bit_count():
            raise ValueError('Arity must be derived from the fixed dense mask')
        label = row['witness_label']
        witness = proposal['witnesses'][label]
        z, logp = F(witness['mark']), F(witness['log2_potential'])
        if not 0 < z < 1 or logp > 0:
            raise ValueError('Invalid frozen rational witness')
        assigned.setdefault((label, z, logp), []).append(mask)
        seen.add(mask)
    if seen != set(ALL_MASKS):
        raise ValueError('Incomplete219-class assignment')
    return tuple((label, z, logp, tuple(masks)) for (label, z, logp), masks in sorted(assigned.items()))


def receipt(precision=256, *, scope=None, progress=None):
    metadata = validate_scope(scope)
    if type(precision) is not int or precision < 128:
        raise ValueError('At least128 precision bits are required')
    assigned = checked_classes()
    if matching_joint_counts(8192, 17)[8, 0] != SOURCE_COUNT:
        raise ArithmeticError('Exact source-support count changed')
    started, entries, bank_costs = perf_counter(), [], []
    with ctx.workprec(precision):
        common_total = arb(0)
        for label, z, logp, masks in assigned:
            bank_started = perf_counter()
            potential = arb(2)**_rational(logp)
            raw = arb_pair_regions(z, 17)
            regions = ordered_tables(raw, 8, 0, directed=True)
            dense_masks = sorted({mask for active in masks
                                  for mask, count in enumerate(comparison_counts(active)) if count})
            # mu=groups disables pair restrictions. The resulting box/sum
            # resource constraints hold for every fixed pair assignment.
            vectors = mixed_vectors_arb(regions, potential, dense_masks, groups=16, mu=16, progress=progress)
            empty_caps = {kind: _upper_max(table[0][0], table[0][1]/potential)
                          for kind, table in vectors.items()}
            unscaled_empty = column_coefficients_arb(empty_caps)
            empty_powers = positive_powers_arb(unscaled_empty, 15)
            prefix = SOURCE_COUNT*17*28*_rational(z)**(-734003)/(255*potential*255**512)
            for active_mask in masks:
                arity = active_mask.bit_count()
                empty, nonempty = fixed_branch_caps_arb(vectors, active_mask, potential)
                tail, by_degree = bivariate_tail_arb(
                    column_coefficients_arb(empty), column_coefficients_arb(nonempty),
                    empty_powers=empty_powers, empty_scale=arb(1)/256**arity)
                common = prefix*(arb(256**(32*arity))/255**(16*arity))*tail
                if not common.is_finite() or not common > 0:
                    raise ArithmeticError('The positive polynomial tail is not finite and positive')
                all_geometry = common*28**15
                entry = dict(active_row_mask=active_mask, active_constituents=arity,
                    witness_label=label, mark=str(z), log2_potential=str(logp),
                    common_pair_failure_log2_upper=str(rational_ceiling(common.log()/arb(2).log())),
                    failure_log2_upper=str(rational_ceiling(all_geometry.log()/arb(2).log())),
                    largest_upper_dense_degree=17+max(range(len(by_degree)), key=lambda i: by_degree[i].upper()))
                entries.append(entry)
                common_total += common
                if progress:
                    progress(dict(stage='directed-fixed-A', **entry))
            cost = dict(witness_label=label, assigned_masks=len(masks), dense_masks=len(dense_masks),
                        elapsed_seconds=perf_counter()-bank_started)
            bank_costs.append(cost)
            if progress:
                progress(dict(stage='directed-bank-complete', **cost))
        entries.sort(key=lambda row: row['active_row_mask'])
        if [row['active_row_mask'] for row in entries] != list(ALL_MASKS):
            raise ArithmeticError('Directed replay did not visit every fixedA exactly once')
        metadata.update(precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
            reviewed_primitives_sha256=PRIMITIVES_SHA256,
            common_pair_failure_log2_upper=str(rational_ceiling(common_total.log()/arb(2).log())),
            failure_log2_upper=str(rational_ceiling((common_total*28**15).log()/arb(2).log())),
            scoped_class_reaches_60_bits=bool(common_total*28**15 < arb(2)**-60),
            entries=entries, bank_costs=bank_costs, elapsed_seconds=perf_counter()-started)
        return metadata


def checked_rows(data, precision):
    expected = dict(EXPECTED_SCOPE, precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_primitives_sha256=PRIMITIVES_SHA256,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    require_fields(data, expected)
    assignment = {mask: (label, str(z), str(logp))
                  for label, z, logp, masks in checked_classes() for mask in masks}
    rows = {}
    for row in data.get('entries', []):
        mask = row.get('active_row_mask')
        if type(mask) is not int or mask not in assignment or mask in rows:
            raise ValueError('Unexpected or duplicate fixed dense mask')
        if type(row.get('active_constituents')) is not int or row['active_constituents'] != mask.bit_count():
            raise ValueError('Arity must match the fixed mask')
        if (row.get('witness_label'), row.get('mark'), row.get('log2_potential')) != assignment[mask]:
            raise ValueError('Exact per-class witness identity changed')
        if any(type(row.get(key)) is not str for key in ('common_pair_failure_log2_upper', 'failure_log2_upper')):
            raise ValueError('Exact rational cell ceilings required')
        rows[mask] = (F(row['common_pair_failure_log2_upper']), F(row['failure_log2_upper']))
    if set(rows) != set(ALL_MASKS):
        raise ValueError('Incomplete219-class receipt')
    return rows


def rounded_union_512(common_ceilings):
    """Sum rounded common-pair cell masses, then multiply28^15 exactly once."""
    with ctx.workprec(512):
        total = sum((arb(2)**_rational(value) for value in common_ceilings), arb(0))
        common = rational_ceiling(total.log()/arb(2).log())
        transferred = rational_ceiling((total*28**15).log()/arb(2).log())
    return common, transferred


def verify(directory=DIRECTORY):
    validate_scope()
    directory = Path(directory)
    rows, receipts, inputs = [], [], []
    for precision in (256, 384):
        path = directory/f'hd_xor_mixed_all_geometry_8_0_receipt{precision}.json'
        data = path.read_bytes()
        result = json.loads(data)
        rows.append(checked_rows(result, precision))
        receipts.append(result)
        inputs.append(dict(filename=path.name, sha256=sha256(data).hexdigest(),
                           precision_bits=precision, elapsed_seconds=result['elapsed_seconds']))
    if rows[0] != rows[1]:
        raise ValueError('Per-class rational ceilings differ between256 and384 bits')
    for key in ('common_pair_failure_log2_upper', 'failure_log2_upper'):
        if receipts[0][key] != receipts[1][key]:
            raise ValueError('Raw aggregate rational ceilings differ between precisions')
    common, transferred = rounded_union_512([rows[0][mask][0] for mask in ALL_MASKS])
    if common < F(receipts[0]['common_pair_failure_log2_upper']) or transferred < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded-cell union lies below the raw aggregate ceiling')
    return dict(EXPECTED_SCOPE, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_primitives_sha256=PRIMITIVES_SHA256, receipt_precision_bits=[256, 384],
        independent_sum_precision_bits=512, exact_partition_checked=True,
        all_per_class_ceilings_agree=True, aggregate_ceilings_agree=True,
        rounded_common_pair_union_log2_upper=str(common),
        rounded_cell_union_log2_upper=str(transferred),
        scoped_class_reaches_60_bits=transferred < -60,
        raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'], inputs=inputs,
        receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--precision', type=int, choices=(256, 384))
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen directed artifact')
    result = verify(args.directory) if args.verify else receipt(args.precision,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
