"""Five frozen fixed-A classes: a bounded FLOAT Holder control, not a proof.

For each deterministic class B, W_B(T)^32=sum_(A in B) beta_A*nu_A(T)^32.
Apply Holder within B, bound its 32-column word, then sum the five classes.
Each class uses its already-frozen common mark and potential. No per-column
row-set or witness resampling is assumed. Labelled masks remain intact.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_dense_holder import uniform_mask_count, SOURCE, SOURCE_SHA256
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_all_geometry_certificate import checked_classes
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import (mixed_mask_vectors, column_log_coefficients,
                                    bivariate_log_tail)
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION


def class_moments(active_masks, q=256, width=8, columns=32, dimension=16):
    """Return exact integer numerators and one common positive denominator."""
    masks = tuple(active_masks)
    if (not masks or len(set(masks)) != len(masks)
            or any(type(mask) is not int or not 0 < mask < 1 << width for mask in masks)
            or not 0 < dimension < columns):
        raise ValueError('A nonempty class of distinct fixed row sets is required')
    denominator = (q-1)**((columns-dimension)*width)
    result = [0]*(1 << width)
    for active in masks:
        factor = (q-1)**((columns-dimension)*(width-active.bit_count()))
        for target in range(1 << width):
            count = uniform_mask_count(q, width, active, target)
            if count:
                result[target] += count**columns*factor
    return tuple(result), denominator


def class_log_weights(numerators, denominator, columns=32):
    if denominator <= 0 or columns < 1 or any(n < 0 for n in numerators):
        raise ValueError('Nonnegative moments with a positive denominator required')
    return np.array([(log2(n)-log2(denominator))/columns if n else -np.inf
                     for n in numerators])


def mask_branch_caps(vectors, log_potential, log_weights):
    """Mix labelled masks before taking one shared-potential quotient."""
    vectors = np.asarray(vectors, dtype=float)
    weights = np.asarray(log_weights, dtype=float)
    if (vectors.ndim != 3 or vectors.shape[1:] != (len(weights), 2)
            or np.any(np.isnan(vectors)) or np.any(np.isposinf(vectors))
            or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))
            or not np.isfinite(log_potential)):
        raise ValueError('Complete nonnegative mask vectors and weights required')
    quotient = np.array([0., log_potential])
    empty = np.max(vectors[:, 0]+weights[0]-quotient, axis=1)
    nonempty = np.max(np.logaddexp2.reduce(
        vectors[:, 1:]+weights[None, 1:, None], axis=1)-quotient, axis=1)
    return empty, nonempty


def frozen_classes():
    raw = SOURCE.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError('Frozen proposal changed')
    combined = json.loads(raw)['combined']
    entries = {row['active_row_mask']: row for row in combined['entries']}
    assigned = checked_classes()
    if len(assigned) != 5:
        raise ValueError('This bounded control has exactly five frozen classes')
    result = []
    for label, z, logp, masks in assigned:
        comparator = float(np.logaddexp2.reduce(
            [entries[mask]['failure_log2_screen'] for mask in masks]))
        result.append((label, z, logp, masks, comparator))
    return tuple(result), combined['failure_log2_screen']


def run(progress=None):
    controls, exact_union = frozen_classes()
    source_count = matching_joint_counts(REGION_LENGTH, 17)[8, 0]
    results = []
    for label, z, logp, masks, comparator in controls:
        started = perf_counter()
        if progress:
            progress(dict(stage='class-start', witness_label=label, fixed_row_sets=len(masks)))
        numerators, denominator = class_moments(masks)
        weights = class_log_weights(numerators, denominator)
        raw = pair_regions_screen(float(z), 17)
        tables = ordered_region_tables(raw, 17, 8, 0)
        vectors = mixed_mask_vectors(tables, float(logp))
        empty, nonempty = mask_branch_caps(vectors, float(logp), weights)
        tail, by_degree = bivariate_log_tail(column_log_coefficients(empty, 16),
                                            column_log_coefficients(nonempty, 16))
        prefix = (log2(source_count)+log2(17)+16*log2(28)-log2(Q-1)
                  -LIMIT*log2(float(z))-float(logp)-512*log2(Q-1))
        bound = prefix+tail
        row = dict(witness_label=label, mark=str(z), log2_potential=str(logp),
            active_row_masks=list(masks), fixed_row_sets=len(masks),
            failure_log2_screen=bound, exact_fixed_A_class_log2_screen=comparator,
            holder_loss_bits=bound-comparator,
            log2_mask_holder_weights=[float(w) if np.isfinite(w) else None for w in weights],
            dense_degree_log2_masses=by_degree, elapsed_seconds=perf_counter()-started)
        results.append(row)
        if progress:
            progress({key: value for key, value in row.items()
                      if key not in ('log2_mask_holder_weights', 'dense_degree_log2_masses', 'active_row_masks')})
    union = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in results]))
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], active_groups=17,
        source_matching_counts=[8, 0], mu_range=[1, 16], sparse_pair_assignments=28**16,
        source_colorings=17, field_order=Q, forbidden_weight=LIMIT, total_columns=32,
        sparse_support_floor=544, dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
        holder_exponent=32, fixed_row_sets=219, deterministic_classes=5,
        comparator_proposal_sha256=SOURCE_SHA256, evaluations=len(results), maximum_evaluations=5,
        fitted_new_witnesses=0, dense_rowsets_charged_once_inside_holder=True,
        extra_dense_beta_factor=False, extra_dense_rowset_factor=False,
        iid_rowset_resampling_assumed=False, setup_dependent_selection=False,
        labelled_dense_masks_retained=True, common_potential_within_each_class=True,
        projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
        entries=results, failure_log2_screen=union,
        exact_selected_fixed_A_union_log2_screen=exact_union, holder_loss_bits=union-exact_union,
        directed=False, whole_code_certificate=False, complete_mixed_class=False)


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
