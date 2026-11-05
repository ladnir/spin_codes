"""Bounded FLOAT control of a fixed-row-set Holder envelope.

Fix a dense row set A through all n=32 columns. Write N_A(T) for the
exact uniform-input count of output mask T and a=|A|. Holder gives

 sum_A prod_j N_A(T_j)/(q-1)^((n-k)*a)
 <= prod_j (sum_A N_A(T_j)^n/(q-1)^((n-k)*a))^(1/n).

This dominates the sum of the fixed-A comparison measures. It is not an
equality and does not assume row sets are independently resampled. The
right side depends only on |T| after summing every A of arity3..8.
Dense beta and row-set counts are already included; do not add them again.

The executable control has one dense group and sixteen two-row groups,
source(8,0), all28^16 fixed sparse pair assignments. It retains each labelled
dense mask through the existing ordered resource DP, so this first control
isolates Holder loss before any further dense-count resource relaxation.
Both X>=544 and Y>=17 are retained. It is FLOAT only, with five frozen
witnesses and no optimizer. No native implementation or certificate claim.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_mixed_all_geometry_refine import load_banks, PENALTY
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import (mixed_mask_vectors, column_log_coefficients,
                                    bivariate_log_tail)
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


DIRECTORY = Path(__file__).resolve().parent
SOURCE = DIRECTORY/'hd_xor_mixed_all_geometry_8_0_refine64.json'
SOURCE_SHA256 = '44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6'
LABELS = ('inherited', 'a8_fit', 'a4_mask240_fit',
          'all_geometry_mask195_fit1', 'all_geometry_mask224_fit2')


def _choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def _sum_count(q, size, nonzero_sum):
    """Nonzero field entries summing to zero or one fixed nonzero value."""
    if nonzero_sum:
        return ((q-1)**size-(-1)**size)//q
    return ((q-1)**size+(q-1)*(-1)**size)//q


def uniform_mask_count(q, width, active_mask, target_mask):
    """Exact count for H=I+11^T, uniform input on the fixed coordinate A."""
    if (q < 2 or q & (q-1) or width < 2 or width % 2
            or not 0 < active_mask < 1 << width or not 0 <= target_mask < 1 << width):
        raise ValueError('An even width, characteristic-two field, and fixed masks are required')
    active, weight = active_mask.bit_count(), target_mask.bit_count()
    if active == width:
        return (q-1)**weight
    outside = ((1 << width)-1)^active_mask
    if not target_mask & outside:
        return _sum_count(q, weight, False)
    if target_mask & outside == outside:
        inside = (target_mask & active_mask).bit_count()
        return (q-1)*_sum_count(q, inside, not (outside.bit_count() & 1))
    return 0


def holder_moments(q=256, width=8, columns=32, dimension=16, minimum_active=3):
    """Exact rational W_w**columns, with row sets and dense beta included.

    For a<width, positive counts require T subset A or A-complement subset T.
    These two families are disjoint. Their cardinalities give the two terms
    below; the full-A case has (q-1)**w unconstrained output vectors.
    """
    if (q < 2 or q & (q-1) or width < 2 or width % 2
            or not 0 < dimension < columns or not 1 <= minimum_active <= width):
        raise ValueError('Invalid field, width, code dimensions, or arity cutoff')
    result = []
    for weight in range(width+1):
        value = F((q-1)**(columns*weight), (q-1)**((columns-dimension)*width))
        zero_sum = _sum_count(q, weight, False)
        for active in range(minimum_active, width):
            overlap = active+weight-width
            second = ((q-1)*_sum_count(q, overlap, not ((width-active) & 1))
                      if overlap >= 0 else 0)
            numerator = (_choose(width-weight, active-weight)*zero_sum**columns
                         +_choose(weight, overlap)*second**columns)
            value += F(numerator, (q-1)**((columns-dimension)*active))
        result.append(value)
    return tuple(result)


def holder_log_weights(moments, columns=32):
    if columns < 1 or any(value <= 0 for value in moments):
        raise ValueError('Positive exact moments and a positive column count required')
    return np.array([(log2(value.numerator)-log2(value.denominator))/columns for value in moments])


def dense_feature_logs(log_weights, groups):
    """D(y,x,u)^groups by (nonempty groups,total occupancy,full groups).

    D=W0+sum_(w=1..width-1) C(width,w) Ww*y*x^w
      +Wwidth*y*x^width*u. This counts complete labelled masks, not iid
    regional bits. The control below retains the masks instead of using
    a further box/sum envelope, but this count polynomial is reusable.
    """
    weights = tuple(float(value) for value in log_weights)
    width = len(weights)-1
    if width < 1 or type(groups) is not int or groups < 0 or not all(np.isfinite(weights)):
        raise ValueError('Finite per-mask weights and a nonnegative group count required')
    terms = [((int(w > 0), w, int(w == width)), weights[w]+log2(comb(width, w)))
             for w in range(width+1)]
    result = {(0, 0, 0): 0.}
    for _ in range(groups):
        following = {}
        for key, value in result.items():
            for degree, coefficient in terms:
                target = tuple(x+y for x, y in zip(key, degree))
                following[target] = float(np.logaddexp2(following.get(target, -np.inf), value+coefficient))
        result = following
    return result


def holder_branch_caps(vectors, log_potential, log_weights):
    """Mix whole labelled dense masks before the shared-potential quotient."""
    vectors = np.asarray(vectors, dtype=float)
    width = len(log_weights)-1
    if (vectors.ndim != 3 or vectors.shape[1:] != (1 << width, 2)
            or np.any(np.isnan(vectors)) or np.any(np.isposinf(vectors))
            or not np.isfinite(log_potential) or not np.all(np.isfinite(log_weights))):
        raise ValueError('Complete finite-or-zero mask vectors and finite witnesses required')
    quotient = np.array([0., log_potential])
    empty = np.max(vectors[:, 0]+log_weights[0]-quotient, axis=1)
    weights = np.array([log_weights[mask.bit_count()] for mask in range(1, 1 << width)])
    nonempty = np.max(np.logaddexp2.reduce(vectors[:, 1:]+weights[None, :, None], axis=1)-quotient, axis=1)
    return empty, nonempty


def evaluate(mark, log_potential, log_weights, progress=None):
    z, logp = F(mark), F(log_potential)
    if not 0 < z < 1 or logp > 0:
        raise ValueError('A fixed rational mark and nonpositive log potential are required')
    started = perf_counter()
    raw = pair_regions_screen(float(z), 17)
    tables = ordered_region_tables(raw, 17, 8, 0)
    vectors = mixed_mask_vectors(tables, float(logp), progress=progress)
    empty, nonempty = holder_branch_caps(vectors, float(logp), log_weights)
    tail, by_degree = bivariate_log_tail(column_log_coefficients(empty, 16),
                                        column_log_coefficients(nonempty, 16))
    source_count = matching_joint_counts(REGION_LENGTH, 17)[8, 0]
    prefix = (log2(source_count)+log2(17)+16*log2(28)-log2(Q-1)
              -LIMIT*log2(float(z))-float(logp)-512*log2(Q-1))
    peak = max(by_degree, key=lambda row: row['log2_mass'])
    return dict(mark=str(z), log2_potential=str(logp), failure_log2_screen=prefix+tail,
        largest_dense_degree=peak['nonempty_dense_columns'],
        largest_sparse_degree=peak['largest_sparse_degree'],
        dense_degree_log2_masses=by_degree, elapsed_seconds=perf_counter()-started)


def frozen_controls():
    data = SOURCE.read_bytes()
    if sha256(data).hexdigest() != SOURCE_SHA256:
        raise ValueError('The frozen comparison proposal changed')
    proposal = json.loads(data)
    banks = load_banks()
    banks.update(proposal['added_banks'])
    controls = []
    for label in LABELS:
        bank = banks[label]
        if (bank['cohort_counts'] != [0, 16, 1] or bank['source_matching_counts'] != [8, 0]
                or bank['mu_range'] != [16, 16] or bank['two_row_geometries'] != 28
                or bank['fixed_row_sets'] != 219 or bank['directed']):
            raise ValueError('The comparison bank has the wrong scope')
        controls.append((label, F(bank['mark']), F(bank['log2_potential']), bank['failure_log2_screen']+PENALTY))
    return tuple(controls), proposal['combined']['failure_log2_screen']


def run(progress=None):
    moments = holder_moments()
    weights = holder_log_weights(moments)
    controls, combined = frozen_controls()
    entries = []
    for label, z, logp, comparator in controls:
        if progress:
            progress(dict(stage='witness-start', label=label, mark=str(z), log2_potential=str(logp)))
        row = evaluate(z, logp, weights, progress)
        row.update(witness_label=label, exact_fixed_A_same_witness_log2_screen=comparator,
                   holder_loss_bits=row['failure_log2_screen']-comparator)
        entries.append(row)
        if progress:
            progress({key: value for key, value in row.items() if key != 'dense_degree_log2_masses'})
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], active_groups=17,
        source_matching_counts=[8, 0], mu_range=[1, 16], sparse_pair_assignments=28**16,
        source_colorings=17, field_order=Q, forbidden_weight=LIMIT, total_columns=32,
        sparse_support_floor=544, dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
        holder_exponent=32, exact_moment_literals=[str(value) for value in moments],
        log2_per_mask_holder_weights=weights.tolist(),
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, fixed_A_replaced_by_proved_word_domination=True,
        iid_rowset_resampling_assumed=False, dense_masks_kept_through_ordered_resource_dp=True,
        additional_dense_feature_resource_relaxation_used=False, same_potential_both_branches=True,
        projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
        selected_fixed_A_comparator_log2_screen=combined,
        comparator_proposal_sha256=SOURCE_SHA256, maximum_evaluations=5, evaluations=len(entries),
        fitted_new_witnesses=0, entries=entries,
        best_log2_screen=min(row['failure_log2_screen'] for row in entries),
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
    print(json.dumps(dict(best_log2_screen=result['best_log2_screen'], evaluations=result['evaluations']), allow_nan=False))
