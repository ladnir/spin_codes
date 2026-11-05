"""Two bounded fixed-A FLOAT refinements; no additional source or mu class."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from hd_xor_mixed_mask_search import BOUNDS, evaluate_at
from hd_xor_paired_gate import CONSTRUCTION


DIRECTORY = Path(__file__).resolve().parent
INHERITED = DIRECTORY/'hd_xor_mixed_mask_floor_8_0_mu16_float.json'
A8_SEARCH = DIRECTORY/'hd_xor_mixed_mask_8_0_mu16_search36.json'
ALL_MASKS = tuple(mask for mask in range(256) if mask.bit_count() >= 3)


def combine_fixed_banks(banks):
    """Choose a proposal separately for each deterministic message class A.

    This validates proposal scope, not its arithmetic. All inputs remain
    FLOAT proposals requiring fresh directed replay of each selected literal.
    """
    assigned = {mask: [] for mask in ALL_MASKS}
    expected = dict(construction=CONSTRUCTION, active_groups=17, source_matching_counts=[8, 0],
        mu_range=[16, 16], fixed_row_sets=219, all_dense_row_sets=True,
        sparse_support_floor=544, dense_nonempty_column_floor=17,
        sparse_support_normalizer_exponent=512, projective_normalizations=1,
        endpoint_potential_factors=1, exact_positive_bivariate_tail=True,
        dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
        directed=False, whole_code_certificate=False)
    witnesses = {}
    for label, bank in banks.items():
        if any(bank.get(key) != value for key, value in expected.items()):
            raise ValueError('All fixed banks must have the same bounded mixed-class scope')
        if sorted(row['active_row_mask'] for row in bank['entries']) != list(ALL_MASKS):
            raise ValueError('Each bank must enumerate each of the219 fixed row sets exactly once')
        z, logp = F(bank['mark']), F(bank['log2_potential'])
        if not 0 < z < 1 or logp > 0:
            raise ValueError('Invalid saved witness')
        witnesses[label] = dict(mark=str(z), log2_potential=str(logp))
        for row in bank['entries']:
            if not np.isfinite(row['failure_log2_screen']):
                raise ValueError('A finite FLOAT proposal value is required')
            assigned[row['active_row_mask']].append(dict(
                active_row_mask=row['active_row_mask'], active_constituents=row['active_constituents'],
                witness_label=label, failure_log2_screen=row['failure_log2_screen']))
    if not banks:
        raise ValueError('At least one complete proposal bank is required')
    rows = [min(candidates, key=lambda row: row['failure_log2_screen']) for candidates in assigned.values()]
    arities = [dict(active_constituents=active, fixed_row_sets=sum(row['active_constituents'] == active for row in rows),
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows
                                                      if row['active_constituents'] == active])))
        for active in range(3, 9)]
    return dict(expected, witnesses=witnesses, entries=rows, arity_summaries=arities,
        witness_choice_is_fixed_by_active_row_mask=True, setup_dependent_selection=False,
        proposal_requires_fresh_directed_replay=True, complete_mixed_class=False,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])))


def search_one(active_mask, initial, direction, progress=None):
    """At most32 evaluations for one fixed row set, with the original bounds."""
    proposals, cache = [], {}

    def objective(point):
        key = tuple(map(float, point))
        if key in cache:
            return cache[key]
        if len(proposals) >= 32:
            raise ArithmeticError('The32-evaluation per-mask budget was exceeded')
        z, logp = (F(value) for value in key)
        result = evaluate_at(z, logp, (active_mask,))
        row = dict(active_row_mask=active_mask, evaluation=len(proposals)+1,
            mark=str(z), log2_potential=str(logp), failure_log2_screen=result['failure_log2_screen'],
            largest_dense_degree=result['entries'][0]['largest_dense_degree'],
            largest_sparse_degree=result['entries'][0]['largest_sparse_degree_at_peak'],
            elapsed_seconds=result['elapsed_seconds'])
        proposals.append(row)
        cache[key] = row['failure_log2_screen']
        if progress:
            progress(dict(stage='fixed-A-search', **row))
        return cache[key]

    start = np.array([float(F(initial['mark'])), float(F(initial['log2_potential']))])
    simplex = np.array([start, start+[direction*.00025, 0.], start+[0., 4.]])
    fit = minimize(objective, start, method='Nelder-Mead', bounds=BOUNDS,
        options=dict(maxfev=32, maxiter=32, initial_simplex=simplex,
                     xatol=1e-10, fatol=.01, adaptive=False))
    best = min(proposals, key=lambda row: row['failure_log2_screen'])
    point = float(F(best['mark'])), float(F(best['log2_potential']))
    boundary = [name for name, value, (low, high) in zip(('mark', 'log2_potential'), point, BOUNDS)
                if abs(value-low) <= 1e-10 or abs(value-high) <= 1e-10]
    return dict(active_row_mask=active_mask, active_constituents=active_mask.bit_count(),
        maximum_evaluations=32, evaluations=len(proposals), proposals=proposals, best=best,
        optimizer='bounded Nelder-Mead', optimizer_message=str(fit.message), bounds=BOUNDS,
        boundary_coordinates=boundary, directed=False, whole_code_certificate=False)


def run(progress=None):
    inherited = json.loads(INHERITED.read_text(encoding='utf-8'))
    a8 = json.loads(A8_SEARCH.read_text(encoding='utf-8'))['selected_all_dense_row_sets']
    banks = dict(inherited=inherited, a8_fit=a8)
    baseline = combine_fixed_banks(banks)
    selected = [max((row for row in baseline['entries'] if row['active_constituents'] == active),
                    key=lambda row: row['failure_log2_screen']) for active in (3, 4)]
    searches = []
    for selected_row, direction in zip(selected, (-1., 1.)):
        mask = selected_row['active_row_mask']
        initial = baseline['witnesses'][selected_row['witness_label']]
        search = search_one(mask, initial, direction, progress)
        improvement = selected_row['failure_log2_screen']-search['best']['failure_log2_screen']
        search['improvement_bits'] = improvement
        searches.append(search)
        if improvement >= 10:
            witness = search['best']
            bank = evaluate_at(witness['mark'], witness['log2_potential'], ALL_MASKS, progress)
            banks[f'a{mask.bit_count()}_mask{mask}_fit'] = bank
    return dict(construction=CONSTRUCTION, source_matching_counts=[8, 0], mu_range=[16, 16],
        source_classes=1, selected_masks=selected, maximum_search_evaluations=64,
        search_evaluations=sum(row['evaluations'] for row in searches), searches=searches,
        proposal_sources=[dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
                          for path in (INHERITED, A8_SEARCH)],
        added_banks={key: value for key, value in banks.items() if key not in ('inherited', 'a8_fit')},
        combined=combine_fixed_banks(banks), directed=False, whole_code_certificate=False,
        proposal_requires_fresh_directed_replay=True, complete_mixed_class=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite the bounded refinement record')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(search_evaluations=result['search_evaluations'],
        selected_masks=[row['active_row_mask'] for row in result['selected_masks']],
        boundary=[row['boundary_coordinates'] for row in result['searches']],
        combined=result['combined']['failure_log2_screen'],
        arities=result['combined']['arity_summaries'])), flush=True)
