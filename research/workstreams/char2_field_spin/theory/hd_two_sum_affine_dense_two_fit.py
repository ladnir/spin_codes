"""Bounded global witness fit for only the affine h2 a8 x a8 subclass.

At most23 optimizer calls plus one fresh exact-literal FLOAT replay.
No source/type expansion or directed certificate. Pinned helpers unchanged.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

from hd_two_sum_affine_dense import affine_class_moments
from hd_two_sum_affine_dense_frozen import Z, LOGP
from hd_two_sum_affine_dense_two import (two_mask_masses, occupancy_vectors, branch_caps, CONSTRUCTION)
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_dense_two_feature import positive_bivariate_tail
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


DIRECTORY = Path(__file__).resolve().parent
INITIAL_NAME = 'hd_two_sum_affine_dense_h2_a8_8_0_frozen.json'
INITIAL_SHA256 = 'a2d482c9c2168a1c1649e70842a0685726ddc6be6c1b8a015c4ceef71190d518'
Z_BOUNDS, P_BOUNDS = (.992, .998), (-100., -40.)


def evaluate(mark, log_potential, mass_logs):
    """Fresh arithmetic at one global witness; no setup-dependent selection."""
    z, logp = F(mark), F(log_potential)
    if not Z_BOUNDS[0] <= float(z) <= Z_BOUNDS[1] or not P_BOUNDS[0] <= float(logp) <= P_BOUNDS[1]:
        raise ValueError('Witness outside the bounded local-fit rectangle')
    started = perf_counter()
    raw = pair_regions_screen(float(z), 17)
    regions = ordered_region_tables(raw, 17, 8, 0)
    vectors = occupancy_vectors(regions, float(logp))
    caps = branch_caps(vectors, mass_logs, float(logp))
    branches = np.array([column_log_coefficients(caps[:, e], 15) for e in range(3)])
    tail, peak = positive_bivariate_tail(branches, 32, 510, 34)
    source_count = matching_joint_counts(8192, 17)[8, 0]
    prefix = (log2(source_count)+log2(comb(17, 2))+15*log2(28)-log2(255)
              -LIMIT*log2(float(z))-float(logp)-480*log2(255))
    return dict(mark=str(z), log2_potential=str(logp), failure_log2_screen=prefix+tail,
                elapsed_seconds=perf_counter()-started, **peak)


def fit(progress=None):
    raw = (DIRECTORY/INITIAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != INITIAL_SHA256:
        raise ValueError('Frozen a8 x a8 initial control changed')
    initial = json.loads(raw)
    if (initial['construction'] != CONSTRUCTION or initial['dense_arities'] != [8, 8]
            or initial['cohort_counts'] != [0, 15, 2] or initial['source_matching_counts'] != [8, 0]
            or initial['mark'] != str(Z) or initial['log2_potential'] != str(LOGP)
            or initial['directed'] is not False):
        raise ValueError('Wrong immutable initial control scope')
    nums, denominator = affine_class_moments((255,))
    weights = np.exp2(class_log_weights(nums, denominator))
    masses = np.asarray(two_mask_masses(weights, weights), dtype=float)
    with np.errstate(divide='ignore'):
        mass_logs = np.log2(masses)
    trials = []
    started = perf_counter()
    def objective(point):
        if len(trials) >= 23:
            raise RuntimeError('No more than23 optimizer evaluations authorized')
        z = float(Z)+1e-4*float(point[0])
        logp = float(LOGP)+float(point[1])
        row = evaluate(F(z), F(logp), mass_logs)
        row['evaluation'] = len(trials)+1
        trials.append(row)
        if progress:
            progress(dict(stage='bounded-global-fit', **row))
        return row['failure_log2_screen']
    bounds = [((Z_BOUNDS[0]-float(Z))/1e-4, (Z_BOUNDS[1]-float(Z))/1e-4),
              (P_BOUNDS[0]-float(LOGP), P_BOUNDS[1]-float(LOGP))]
    result = minimize(objective, [0., 0.], method='Nelder-Mead', bounds=bounds,
        options=dict(maxfev=23, xatol=1e-6, fatol=1e-6,
                     initial_simplex=np.array([[0., 0.], [-3., 0.], [0., -4.]])))
    best = min(trials, key=lambda row: row['failure_log2_screen'])
    frozen = evaluate(F(best['mark']), F(best['log2_potential']), mass_logs)
    if abs(frozen['failure_log2_screen']-best['failure_log2_screen']) > 1e-8:
        raise ArithmeticError('Fresh exact-literal FLOAT replay disagreed with winning trial')
    boundary = (min(abs(float(F(best['mark']))-bound) for bound in Z_BOUNDS) < 1e-10
                or min(abs(float(F(best['log2_potential']))-bound) for bound in P_BOUNDS) < 1e-8)
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 15, 2], active_groups=17,
        source_matching_counts=[8, 0], dense_arities=[8, 8], fixed_dense_rowsets=[[255], [255]],
        ordered_fixed_dense_rowset_assignments=1, complete_dense_rowset_coverage=False,
        source_colorings=comb(17, 2), sparse_pair_assignments=28**15,
        field_order=256, labels=list(range(8, 16)), forbidden_weight=LIMIT, total_columns=32,
        sparse_support_floor=510, dense_nonempty_column_floor=34, sparse_normalizer_exponent=480,
        dense_floor_is_aggregate_relaxation=True, one_global_mark_and_potential=True,
        old_xor_mask_counts_used=False, old_certificates_transferred=False,
        initial_filename=INITIAL_NAME, initial_sha256=INITIAL_SHA256,
        initial_failure_log2_screen=initial['failure_log2_screen'], bounds=dict(mark=Z_BOUNDS, log2_potential=P_BOUNDS),
        optimizer='bounded-Nelder-Mead', optimizer_evaluations=len(trials),
        fresh_literal_replays=1, total_evaluations=len(trials)+1, maximum_evaluations=24,
        optimizer_success=bool(result.success), optimizer_message=str(result.message),
        winning_witness_hits_boundary=boundary, trials=trials, winning_literal=frozen,
        failure_log2_screen=frozen['failure_log2_screen'], elapsed_seconds=perf_counter()-started,
        directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite bounded fit evidence')
    result = fit(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'trials'}, allow_nan=False))
