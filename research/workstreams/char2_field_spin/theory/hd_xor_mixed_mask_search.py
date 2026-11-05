"""At most36 FLOAT witness proposals for the single A255 mixed gate.

No directed or whole-code claim is made. Source(8,0), mu16, the exact
X/Y floors, fixed A, and all normalizations are inherited from the isolated
mixed-mask helper. Only z and the common endpoint potential are proposed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

from hd_xor_fixed_mask_floor import beta_log2
from hd_xor_geometry_strata import geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import inherited_witnesses, ordered_region_tables, SOURCE
from hd_xor_mixed_mask_floor import (mixed_mask_vectors, fixed_mask_branch_caps,
    column_log_coefficients, log_powers, bivariate_log_tail)
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


BOUNDS = ((.99, .9998), (-100., -20.))


def evaluate_at(mark, log_potential, active_masks=(255,), progress=None):
    """One exact-literal FLOAT witness and an explicitly listed fixed-A scope."""
    z, logp = F(mark), F(log_potential)
    masks = tuple(active_masks)
    if not 0 < z < 1 or logp > 0 or not masks or len(set(masks)) != len(masks):
        raise ValueError('One valid witness and distinct nonempty fixed row sets required')
    if any(type(mask) is not int or not 0 < mask < 256 or mask.bit_count() < 3 for mask in masks):
        raise ValueError('Every dense row mask must have arity3..8')
    started = perf_counter()
    raw = pair_regions_screen(float(z), 17)
    tables = ordered_region_tables(raw, 17, 8, 0)
    vectors = mixed_mask_vectors(tables, float(logp))
    empty_caps = np.max(vectors[:, 0]-[0., float(logp)], axis=1)
    empty_poly = column_log_coefficients(empty_caps, 16)
    empty_powers = log_powers(empty_poly, 15)
    source_count = matching_joint_counts(REGION_LENGTH, 17)[8, 0]
    geometries = geometry_stratum_counts(16)[16]
    prefix = (log2(source_count)+log2(17)+log2(geometries)-log2(Q-1)
              -LIMIT*log2(float(z))-float(logp)-512*log2(Q-1))
    rows = []
    for mask in masks:
        empty, nonempty = fixed_mask_branch_caps(vectors, mask, float(logp))
        empty_branch = column_log_coefficients(empty, 16)
        nonempty_branch = column_log_coefficients(nonempty, 16)
        shift = -8*mask.bit_count()
        if not np.allclose(empty_branch, empty_poly+shift, rtol=0, atol=1e-11):
            raise ArithmeticError('Cached empty-branch powers have the wrong scale')
        tail, degrees = bivariate_log_tail(empty_branch, nonempty_branch,
            empty_powers=empty_powers, empty_shift=shift)
        peak = max(degrees, key=lambda row: row['log2_mass'])
        row = dict(active_row_mask=mask, active_constituents=mask.bit_count(),
            failure_log2_screen=prefix+beta_log2(mask.bit_count())+tail,
            largest_dense_degree=peak['nonempty_dense_columns'],
            largest_sparse_degree_at_peak=peak['largest_sparse_degree'],
            dense_degree_log2_masses=degrees)
        rows.append(row)
        if progress and (len(rows) == 1 or len(rows) % 20 == 0 or len(rows) == len(masks)):
            progress(dict(stage='fixed-row-set-tail', completed=len(rows), total=len(masks),
                          **{key: value for key, value in row.items() if key != 'dense_degree_log2_masses'}))
    arities = []
    for active in range(3, 9):
        chosen = [row for row in rows if row['active_constituents'] == active]
        if chosen:
            arities.append(dict(active_constituents=active, fixed_row_sets=len(chosen),
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in chosen]))))
    return dict(construction=CONSTRUCTION, active_groups=17, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], source_classes=1, source_supports=source_count,
        dense_source_colorings=17, mu_range=[16, 16], two_row_geometries=geometries,
        fixed_row_sets=len(masks), all_dense_row_sets=set(masks) == {m for m in range(256) if m.bit_count() >= 3},
        mark=str(z), log2_potential=str(logp), sparse_support_floor=544,
        dense_nonempty_column_floor=17, total_columns=32, sparse_support_normalizer_exponent=512,
        dense_beta_formula='q^(32a)/(q-1)^(16a)', field_order=Q, forbidden_weight=LIMIT,
        pair_swaps=True, alpha=1, dense_row_set_fixed_across_columns=True,
        sparse_geometry_fixed_across_columns=True, exact_uniform_dense_mask_mixture=True,
        dense_input_includes_zero=True, masks_and_cancellations_mixed_before_potential=True,
        same_potential_both_dense_branches=True, exact_positive_bivariate_tail=True,
        chernoff_dense_floor=False, projective_normalizations=1, endpoint_potential_factors=1,
        product_density_comparison_used=False, complete_mixed_class=False,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])),
        arity_summaries=arities, entries=rows, elapsed_seconds=perf_counter()-started)


def bounded_search(progress=None):
    inherited, inherited_hash = inherited_witnesses()
    initial = next((z, logp) for lo, hi, z, logp in inherited[8, 0] if lo <= 16 <= hi)
    proposals, cache = [], {}

    def objective(point):
        key = tuple(float(value) for value in point)
        if key in cache:
            return cache[key]
        if len(proposals) >= 36:
            raise ArithmeticError('The fixed36-evaluation budget was exceeded')
        literals = tuple(F(value) for value in key)
        result = evaluate_at(*literals)
        trial = dict(evaluation=len(proposals)+1, mark=str(literals[0]),
            log2_potential=str(literals[1]), failure_log2_screen=result['failure_log2_screen'],
            largest_dense_degree=result['entries'][0]['largest_dense_degree'],
            largest_sparse_degree=result['entries'][0]['largest_sparse_degree_at_peak'],
            elapsed_seconds=result['elapsed_seconds'])
        proposals.append(trial)
        cache[key] = trial['failure_log2_screen']
        if progress:
            progress(trial)
        return cache[key]

    start = np.array([float(value) for value in initial])
    simplex = np.array([start, start+[-.0015, 0.], start+[0., 15.]])
    fit = minimize(objective, start, method='Nelder-Mead', bounds=BOUNDS,
        options=dict(maxfev=36, maxiter=36, initial_simplex=simplex,
                     xatol=1e-10, fatol=.01, adaptive=False))
    best = min(proposals, key=lambda row: row['failure_log2_screen'])
    point = (float(F(best['mark'])), float(F(best['log2_potential'])))
    boundary = [name for name, value, (low, high) in zip(('mark', 'log2_potential'), point, BOUNDS)
                if abs(value-low) <= 1e-10 or abs(value-high) <= 1e-10]
    result = dict(construction=CONSTRUCTION, source_matching_counts=[8, 0], mu_range=[16, 16],
        active_row_mask=255, active_constituents=8, evaluated_dense_row_sets=1,
        inherited_witness_source_sha256=inherited_hash, bounds=BOUNDS,
        maximum_evaluations=36, evaluations=len(proposals), optimizer='bounded Nelder-Mead',
        optimizer_message=str(fit.message), boundary_coordinates=boundary,
        best=best, proposals=proposals, directed=False, whole_code_certificate=False,
        proposal_requires_fresh_directed_replay=True, complete_mixed_class=False,
        exact_positive_bivariate_tail=True, sparse_support_floor=544,
        dense_nonempty_column_floor=17, fitted_new_witnesses=len(proposals))
    improvement = proposals[0]['failure_log2_screen']-best['failure_log2_screen']
    # One all-A follow-up is permitted only after a substantial improvement.
    if improvement >= 50:
        result['selected_all_dense_row_sets'] = evaluate_at(best['mark'], best['log2_potential'],
            [mask for mask in range(1, 256) if mask.bit_count() >= 3], progress)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite the bounded-search record')
    result = bounded_search(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(evaluations=result['evaluations'], best=result['best'],
        boundary=result['boundary_coordinates'], all_A_bound=result.get('selected_all_dense_row_sets', {}).get('failure_log2_screen'))), flush=True)
