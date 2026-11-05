"""At most two32-trial FLOAT refinements for one all-geometry mixed class.

The old evaluator counts28 common-pair choices. Its reviewed resource
envelope is universal over fixed pair assignments; this adapter adds
15log2(28) only after each complete deterministic fixed-A bank is evaluated.
No saved helper, literal, or directed receipt is changed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_mixed_mask_refine import (ALL_MASKS, INHERITED, A8_SEARCH,
    combine_fixed_banks, search_one)
from hd_xor_mixed_mask_search import evaluate_at
from hd_xor_mixed_mask_witnesses import PROPOSAL_SHA256, checked_classes
from hd_xor_paired_gate import CONSTRUCTION


DIRECTORY = Path(__file__).resolve().parent
REFINEMENT = DIRECTORY/'hd_xor_mixed_mask_8_0_mu16_refine64.json'
PINS = {
    INHERITED.name: '2b3fe0af20331ca123c5370d729c4ca031e3357aab98be3ea09436d24fd3546c',
    A8_SEARCH.name: '88faefb4973e9a10599da7b6997158a8f596886b81b5fcb0e38ce38f0d025cc2',
    REFINEMENT.name: PROPOSAL_SHA256,
}
PENALTY = 15*log2(28)
TARGET = -60.
MINIMUM_IMPROVEMENT = 5.


def load_banks():
    sources = {}
    for name, expected in PINS.items():
        data = (DIRECTORY/name).read_bytes()
        if sha256(data).hexdigest() != expected:
            raise ValueError(f'Frozen input changed: {name}')
        sources[name] = json.loads(data)
    banks = dict(inherited=sources[INHERITED.name],
                 a8_fit=sources[A8_SEARCH.name]['selected_all_dense_row_sets'])
    banks.update(sources[REFINEMENT.name]['added_banks'])
    combined = combine_fixed_banks(banks)
    frozen = {mask: (label, z, logp) for label, z, logp, masks in checked_classes() for mask in masks}
    for row in combined['entries']:
        label, z, logp = frozen[row['active_row_mask']]
        if (row['witness_label'] != label or
                F(combined['witnesses'][label]['mark']) != z or
                F(combined['witnesses'][label]['log2_potential']) != logp):
            raise ValueError('Initial assignment differs from the reviewed frozen assignment')
    if max(combined['entries'], key=lambda row: row['failure_log2_screen'])['active_row_mask'] != 195:
        raise ValueError('The initial worst fixed row set must be195')
    return banks


def all_geometry_summary(banks):
    """Apply only the reviewed one-time geometry count change, then sum A."""
    for bank in banks.values():
        if (bank.get('cohort_counts') != [0, 16, 1] or bank.get('two_row_geometries') != 28
                or bank.get('field_order') != 256 or bank.get('forbidden_weight') != 734003
                or bank.get('total_columns') != 32):
            raise ValueError('Only the pinned common-pair mixed scope can be transferred')
    common = combine_fixed_banks(banks)
    entries = []
    for row in common['entries']:
        if row['active_constituents'] != row['active_row_mask'].bit_count():
            raise ValueError('Arity must agree with the fixed dense row mask')
        entries.append(dict(row, common_pair_failure_log2_screen=row['failure_log2_screen'],
                            failure_log2_screen=row['failure_log2_screen']+PENALTY))
    arities = [dict(active_constituents=a, fixed_row_sets=sum(row['active_constituents'] == a for row in entries),
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries
                                                      if row['active_constituents'] == a])))
        for a in range(3, 9)]
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], active_groups=17,
        source_matching_counts=[8, 0], source_classes=1, field_order=256,
        forbidden_weight=734003, total_columns=32, pair_swaps=True, alpha=1,
        mu_range=[1, 16], sparse_pair_assignments=28**16,
        base_sparse_pair_assignments=28, exact_mass_multiplier=28**15,
        fixed_row_sets=219, sparse_support_floor=544, dense_nonempty_column_floor=17,
        sparse_support_normalizer_exponent=512, exact_positive_bivariate_tail=True,
        dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
        geometry_resampled_across_columns=False, setup_dependent_selection=False,
        witness_choice_is_fixed_by_active_row_mask=True,
        witnesses=common['witnesses'], entries=entries, arity_summaries=arities,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
        directed=False, whole_code_certificate=False, complete_mixed_class=False,
        proposal_requires_fresh_directed_replay=True)


def run(progress=None):
    started = perf_counter()
    banks = load_banks()
    baseline = all_geometry_summary(banks)
    combined = baseline
    searches, added = [], {}
    for index in range(2):
        if combined['failure_log2_screen'] <= TARGET:
            break
        selected = max(combined['entries'], key=lambda row: row['failure_log2_screen'])
        mask = selected['active_row_mask']
        initial = combined['witnesses'][selected['witness_label']]
        if progress:
            progress(dict(stage='selected-worst-A', search=index+1, **selected, initial=initial))
        search = search_one(mask, initial, 1., progress)
        improvement = selected['common_pair_failure_log2_screen']-search['best']['failure_log2_screen']
        search.update(improvement_bits=improvement, initial_witness=initial,
                      trial_values_scope='common-pair-count28; add15log2(28) for all assignments')
        searches.append(search)
        if improvement < MINIMUM_IMPROVEMENT:
            break
        best = search['best']
        label = f'all_geometry_mask{mask}_fit{index+1}'
        bank = evaluate_at(best['mark'], best['log2_potential'], ALL_MASKS, progress)
        added[label] = bank
        banks[label] = bank
        combined = all_geometry_summary(banks)
        if progress:
            progress(dict(stage='combined-all-geometry', search=index+1,
                failure_log2_screen=combined['failure_log2_screen'],
                worst=max(combined['entries'], key=lambda row: row['failure_log2_screen'])))
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], mu_range=[1, 16], sparse_pair_assignments=28**16,
        maximum_search_evaluations=64, maximum_evaluations_per_mask=32,
        search_evaluations=sum(search['evaluations'] for search in searches),
        maximum_full_bank_evaluations=2, full_bank_evaluations=len(added),
        target_log2_screen=TARGET, minimum_material_improvement_bits=MINIMUM_IMPROVEMENT,
        input_pins=PINS, baseline=baseline, searches=searches, added_banks=added, combined=combined,
        elapsed_seconds=perf_counter()-started, directed=False, whole_code_certificate=False,
        proposal_requires_fresh_directed_replay=True, old_receipts_unchanged=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen improvement proposal')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(search_evaluations=result['search_evaluations'],
        full_bank_evaluations=result['full_bank_evaluations'],
        failure_log2_screen=result['combined']['failure_log2_screen'],
        elapsed_seconds=result['elapsed_seconds'])), flush=True)
