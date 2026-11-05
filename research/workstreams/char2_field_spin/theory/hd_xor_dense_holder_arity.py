"""Last bounded FLOAT control: split each frozen witness class by arity.

Exactly five existing regional/vector banks are rebuilt once. Each nonempty
class/arity intersection uses that class's unchanged rational mark and
potential. The output is a diagnostic, not a directed certificate.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_dense_holder import SOURCE, SOURCE_SHA256
from hd_xor_dense_holder_classes import (class_moments, class_log_weights,
    mask_branch_caps, frozen_classes)
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import (mixed_mask_vectors, column_log_coefficients,
                                    bivariate_log_tail)
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION


def arity_classes(masks):
    return tuple((arity, tuple(mask for mask in masks if mask.bit_count() == arity))
                 for arity in sorted({mask.bit_count() for mask in masks}))


def run(progress=None):
    controls, comparator = frozen_classes()
    source_count = matching_joint_counts(REGION_LENGTH, 17)[8, 0]
    original = json.loads(SOURCE.read_bytes())['combined']
    rows_by_mask = {row['active_row_mask']: row for row in original['entries']}
    entries, costs = [], []
    for label, z, logp, masks, _ in controls:
        started = perf_counter()
        if progress:
            progress(dict(stage='bank-start', witness_label=label))
        raw = pair_regions_screen(float(z), 17)
        tables = ordered_region_tables(raw, 17, 8, 0)
        vectors = mixed_mask_vectors(tables, float(logp))
        prefix = (log2(source_count)+log2(17)+16*log2(28)-log2(Q-1)
                  -LIMIT*log2(float(z))-float(logp)-512*log2(Q-1))
        for arity, selected in arity_classes(masks):
            nums, denominator = class_moments(selected)
            weights = class_log_weights(nums, denominator)
            empty, nonempty = mask_branch_caps(vectors, float(logp), weights)
            tail, _ = bivariate_log_tail(column_log_coefficients(empty, 16),
                                         column_log_coefficients(nonempty, 16))
            exact = float(np.logaddexp2.reduce(
                [rows_by_mask[mask]['failure_log2_screen'] for mask in selected]))
            row = dict(witness_label=label, mark=str(z), log2_potential=str(logp),
                active_constituents=arity, active_row_masks=list(selected), fixed_row_sets=len(selected),
                failure_log2_screen=prefix+tail, exact_fixed_A_subclass_log2_screen=exact,
                holder_loss_bits=prefix+tail-exact)
            entries.append(row)
            if progress:
                progress(dict(stage='arity-class-complete', **row))
        costs.append(dict(witness_label=label, elapsed_seconds=perf_counter()-started))
    union = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries]))
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], active_groups=17,
        source_matching_counts=[8, 0], mu_range=[1, 16], sparse_pair_assignments=28**16,
        source_colorings=17, field_order=Q, forbidden_weight=LIMIT, total_columns=32,
        sparse_support_floor=544, dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
        holder_exponent=32, fixed_row_sets=219, deterministic_classes=len(entries),
        comparator_proposal_sha256=SOURCE_SHA256, regional_bank_evaluations=len(controls),
        maximum_regional_bank_evaluations=5, fitted_new_witnesses=0,
        witness_classes_refined_only_by_arity=True, same_frozen_witness_each_intersection=True,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, iid_rowset_resampling_assumed=False,
        setup_dependent_selection=False, labelled_dense_masks_retained=True,
        common_potential_within_each_class=True, projective_normalizations=1,
        endpoint_potential_factors=1, alpha=1, entries=entries, bank_costs=costs,
        failure_log2_screen=union, exact_selected_fixed_A_union_log2_screen=comparator,
        holder_loss_bits=union-comparator, directed=False, whole_code_certificate=False,
        complete_mixed_class=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded diagnostic artifact')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}), flush=True)
