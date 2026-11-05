"""Bounded FLOAT full-arity g16 dense gates, without labelled-mask enumeration.

Fix the complete source support and each group's constituent row set. For
one full-arity group, invertibility gives N(T)=(q-1)^|T| for each column
mask T. Its comparison weight is N(T)/(q-1)^8. Across16 columns this pays
the16 RS[16,8] row normalizers once. No row set is sampled anew.

For h<=2 such groups, a column DP records how many groups have appeared in
earlier regions. Their identities are exchangeable under this full-arity
counting measure after applying the universal source-conditioned occupancy
envelopes. Actual subsets of the fixed source support need not be
exchangeable. The matrices still act in labelled regional order.
Sparse holes and dense masks are mixed before one potential quotient.
The final seen count marks nonempty group-columns; the retained aggregate
floor9h relaxes the individual per-group floors9. This is FLOAT evidence,
not a directed receipt, and covers only source(8,0) at total support17.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS, ordered_regions
from hd_xor_dense_two_feature import positive_bivariate_tail
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


# Frozen winning g8 h1 literal; no fit and no claim that its receipt transfers.
Z = F('8970988817957597/9007199254740992')
LOGP = F('-2266954815/33554432')
GROUPS, COLUMNS, DIMENSION = 17, 16, 8
SOURCE = (8, 0)


def seen_transitions(groups, q=256):
    """(old_seen,new_seen,occupancy,integer_mass), for any fixed seen set.

    Choose k nonzeros among its s members and j among the h-s unseen
    members. The next seen count is s+j, not s+k+j. No quotient is taken.
    """
    if type(groups) is not int or not 1 <= groups <= 2 or type(q) is not int or q < 3:
        raise ValueError('One or two full-arity dense groups and integer q>=3 required')
    return tuple((s, s+j, k+j, comb(s, k)*comb(groups-s, j)*(q-1)**(k+j))
                 for s in range(groups+1) for k in range(s+1)
                 for j in range(groups-s+1))


def _check_deadline(deadline):
    if deadline is not None and perf_counter() >= deadline:
        raise TimeoutError('The bounded FLOAT gate exhausted its compute budget')


def column_vectors(regions, log_potential, sparse_groups, dense_groups, q=256, *, deadline=None):
    """Unnormalized positive vectors indexed by sparse category and final seen count.

    A suffix entry V[s,e,k] sums all assignments of k labelled sparse holes
    to that suffix, starting with s dense groups seen and ending with e.
    Assigning j holes to this region contributes C(k,j). Its occupancy is
    t-j+o, where t is the number of occupied sparse groups and o is the
    dense occupancy. The common terminal vector is (1,2^log_potential).
    """
    transitions = seen_transitions(dense_groups, q)
    tables = np.asarray(regions, dtype=float)
    if (type(sparse_groups) is not int or not 0 <= sparse_groups <= 17
            or tables.ndim != 4 or tables.shape[2:] != (2, 2)
            or not 1 <= len(tables) < q-1
            or tables.shape[1] <= sparse_groups+dense_groups
            or not np.isfinite(log_potential)
            or np.any(np.isnan(tables)) or np.any(np.isposinf(tables))):
        raise ValueError('Complete finite-or-zero ordered log matrices and bounded sparse count required')
    transition_logs = tuple((s, next_s, occupancy, log2(mass))
                            for s, next_s, occupancy, mass in transitions)

    @lru_cache(None)
    def suffix(occupied):
        _check_deadline(deadline)
        states = dense_groups+1
        values = np.full((states, states, occupied+1, 2), -np.inf)
        for s in range(states):
            values[s, s, 0] = [0., log_potential]
        binomials = [np.array([log2(comb(k, holes)) for k in range(holes, occupied+1)])
                     for holes in range(occupied+1)]
        for table in reversed(tables):
            following = np.full_like(values, -np.inf)
            for holes in range(occupied+1):
                length = occupied+1-holes
                for s, next_s, occupancy, mass in transition_logs:
                    matrix = table[occupied-holes+occupancy]
                    previous = values[next_s, :, :length]
                    action = np.logaddexp2(matrix[:, 0]+previous[..., 0, None],
                                          matrix[:, 1]+previous[..., 1, None])
                    target = following[s, :, holes:]
                    np.logaddexp2(target, action+mass+binomials[holes][None, :, None], out=target)
            values = following
        return values[0]

    result = []
    for single, double in categories(sparse_groups):
        values = suffix(single+double)
        probabilities = np.array([log2(comb(double, k))
            +(double-k)*log2(q-1-len(tables))-double*log2(q-1) for k in range(double+1)])
        result.append(np.logaddexp2.reduce(values[:, :double+1]+probabilities[None, :, None], axis=1))
    _check_deadline(deadline)
    return np.array(result)


def run(*, max_seconds=120, progress=None):
    """Exactly h1 and h2, one source geometry and one inherited witness."""
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 120):
        raise ValueError('An explicit positive budget at most120 seconds is required')
    started = perf_counter()
    deadline = started+max_seconds
    result = dict(construction=CONSTRUCTION, field_order=256, field_polynomial='0x11b',
        labels=list(LABELS), active_groups=GROUPS, source_matching_counts=list(SOURCE),
        source_supports=matching_joint_counts(8192, GROUPS)[SOURCE],
        forbidden_weight=LIMIT, total_columns=COLUMNS, constituent_dimension=DIMENSION,
        requested_dense_groups=[1, 2], full_arity=16, fixed_row_sets_per_dense_group=1,
        dense_rowsets_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
        iid_rowset_resampling_assumed=False, local_full_arity_mask_count='255^popcount(T)',
        local_dense_comparison_weight='255^(popcount(T)-8)',
        dense_identity_exchangeability_only=True, labelled_region_order_preserved=True,
        hole_and_dense_mixtures_before_potential=True, pair_swaps=True,
        one_global_mark_and_potential=True, mark=str(Z), log2_potential=str(LOGP),
        witness_label='inherited-old-dense-a8', witness_banks=1, fitted_new_witnesses=0,
        raw_regional_banks=0, projective_normalizations=1, endpoint_potential_factors=1,
        alpha=1, directed=False, whole_code_certificate=False,
        complete_source_coverage=False, complete_dense_rowset_coverage=False,
        old_certificates_transferred=False, max_seconds=max_seconds, entries=[])
    try:
        _check_deadline(deadline)
        raw = pair_regions_screen(float(Z), GROUPS)
        result['raw_regional_banks'] = 1
        regions = ordered_regions(raw, GROUPS, *SOURCE)
        _check_deadline(deadline)
        for dense_groups in (1, 2):
            entry_started = perf_counter()
            sparse_groups = GROUPS-dense_groups
            vectors = column_vectors(regions, float(LOGP), sparse_groups, dense_groups, deadline=deadline)
            caps = np.max(vectors-np.array([0., float(LOGP)]), axis=2)-DIMENSION*dense_groups*log2(255)
            branches = np.array([column_log_coefficients(caps[:, e], sparse_groups)
                                 for e in range(dense_groups+1)])
            tail, peak = positive_bivariate_tail(branches, COLUMNS,
                                                 18*sparse_groups, 9*dense_groups)
            _check_deadline(deadline)
            prefix = (log2(result['source_supports'])+log2(comb(GROUPS, dense_groups))
                +sparse_groups*log2(120)-log2(255)-LIMIT*log2(float(Z))-float(LOGP)
                -16*sparse_groups*log2(255))
            entry = dict(cohort_counts=[0, sparse_groups, dense_groups],
                source_colorings=comb(GROUPS, dense_groups), sparse_pair_assignments=120**sparse_groups,
                sparse_support_floor=18*sparse_groups, dense_nonempty_column_floor=9*dense_groups,
                sparse_normalizer_exponent=16*sparse_groups,
                dense_per_column_normalizer_exponent=DIMENSION*dense_groups,
                dense_nonempty_floor_is_aggregate_relaxation=dense_groups > 1,
                seen_count_states=dense_groups+1, sparse_categories=len(caps),
                dense_labelled_masks_enumerated=0, failure_log2_screen=prefix+tail,
                prefix_log2=prefix, positive_tail_log2_screen=tail, **peak,
                elapsed_seconds=perf_counter()-entry_started)
            result['entries'].append(entry)
            if progress:
                progress(entry)
        result.update(status='complete-selected-gates', complete_selected_gates=True,
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in result['entries']])))
    except TimeoutError:
        result.update(status='compute-budget-exceeded', complete_selected_gates=False)
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded FLOAT artifact')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, 'hd_two_sum_gf16_sparse.py', 'hd_two_sum_affine_outer.py',
                     'hd_xor_dense_two_feature.py', 'hd_xor_mixed_mask_floor.py',
                     'hd_xor_pair_route.py', 'hd_xor_group_geometry.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False))
