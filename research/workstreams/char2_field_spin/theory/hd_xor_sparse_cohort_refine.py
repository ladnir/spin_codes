"""One bounded joint fit and one conditional fixed-bank sparse-cohort screen.

Only two source(0,7) cells enter the32-evaluation objective.  If their
log-sum improves by at least10 bits, the resulting ONE literal witness
is screened over all81 sources/all18 sparse r17 cohorts.  Existing
controls are read by exact hash and never overwritten.  All new numerical
values are FLOAT proposals, including unions using some directed inputs.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

from hd_xor_geometry_strata import MultiplicityCapsFloat
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_pair_route import pair_regions_screen
from hd_xor_sparse_cohorts import (GROUPS, CANONICAL_BINS, CONSTRUCTION, Q,
    REGION_LENGTH, LIMIT, class_row, geometry_bins)


MAX_EVALUATIONS = 32
MINIMUM_GAIN_BITS = 10.
BOUNDS = ((.99, .9998), (-100., -20.))
INITIAL = (F('4487948476156805/4503599627370496'), F('-8814329/131072'))
FIT_SOURCE = (0, 7)
FIT_CELLS = ((2, (15, 15)), (3, (9, 14)))
HOMOGENEOUS_DIRECTED_CEILING = F('-50891481584655/1000000000000')
INPUT_PINS = {
    'hd_xor_sparse_cohorts17_all_sources_two_banks_float.json':
        'ff786dc9a2c4a2791b759e7b244eab09206cc0213df43f7100eb71795def63ab',
    'hd_xor_sparse_mixed_refine32_all_sources_float.json':
        '10a7116e0376f203fee0256f524ce1c3614cbc358f872dba0afb863d9d831b75',
    'hd_xor_paired_joint_gate17_verified256.json':
        'f1a30c19232b117b116545f66ef5cf7885e0203b7f6285684334fa501079ff35',
    'hd_xor_paired_joint_gate17_verified384.json':
        'd6ecffaf84c570c263fdadb7e20b47ad6349085f7b6faca41d17effe1073e94c',
}


@lru_cache(None)
def source_counts():
    counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    if (set(counts) != {(a, b) for a in range(9) for b in range(9)}
            or sum(counts.values()) != comb(REGION_LENGTH, GROUPS)):
        raise ArithmeticError('Expected the complete81-source partition')
    return counts


@lru_cache(None)
def partitions():
    return {n2: geometry_bins(n2) for n2 in range(18)}


def canonical_mu(high):
    """Exact cap reuse on used categories d<=n2, not a loosened envelope."""
    if type(high) is not int or not 0 <= high <= 17:
        raise ValueError('Expected a truncated multiplicity endpoint0..17')
    return next(upper for low, upper in CANONICAL_BINS if low <= max(1, high) <= upper)


def evaluate_joint(z, logp):
    z, logp = F(z), F(logp)
    if not BOUNDS[0][0] <= float(z) <= BOUNDS[0][1] or not BOUNDS[1][0] <= float(logp) <= BOUNDS[1][1]:
        raise ValueError('Witness is outside the explicitly bounded fit domain')
    raw = pair_regions_screen(float(z), GROUPS)
    evaluator = MultiplicityCapsFloat(ordered_region_tables(raw, GROUPS, *FIT_SOURCE), float(logp), GROUPS)
    caps = {mu: evaluator(mu) for mu in (14, 15)}
    entries = []
    for n1, bounds in FIT_CELLS:
        geometries = next(count for low, high, count in partitions()[17-n1] if (low, high) == bounds)
        row = class_row(n1, bounds, geometries, source_counts()[FIT_SOURCE],
            caps[bounds[1]], 'joint-cohort-refit', z, logp)
        row['source_matching_counts'] = list(FIT_SOURCE)
        entries.append(row)
    return dict(mark=str(z), log2_potential=str(logp), entries=entries,
        objective_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])))


def fit_joint(*, progress=None):
    started, proposals, cache = perf_counter(), [], {}
    def objective(point):
        key = tuple(map(float, point))
        if key in cache:
            return cache[key]
        if len(proposals) >= MAX_EVALUATIONS:
            raise ArithmeticError('The joint fit exceeded its32-trial budget')
        row = evaluate_joint(*map(F, key))
        row['evaluation'] = len(proposals)+1
        proposals.append(row)
        cache[key] = row['objective_log2_screen']
        if progress:
            progress(dict(stage='joint-fit', **row))
        return cache[key]
    start = np.array(list(map(float, INITIAL)))
    fit = minimize(objective, start, method='Nelder-Mead', bounds=BOUNDS,
        options=dict(maxfev=MAX_EVALUATIONS, maxiter=MAX_EVALUATIONS,
            initial_simplex=np.array([start, start+[-.001, 0], start+[0, 8.]]),
            xatol=1e-10, fatol=.01, adaptive=False))
    best = min(proposals, key=lambda row: row['objective_log2_screen'])
    gain = proposals[0]['objective_log2_screen']-best['objective_log2_screen']
    return dict(best=best, proposals=proposals, evaluations=len(proposals),
        maximum_evaluations=MAX_EVALUATIONS, gain_bits=gain,
        minimum_gain_for_full_screen=MINIMUM_GAIN_BITS,
        full_screen_authorized_by_gain=gain >= MINIMUM_GAIN_BITS,
        optimizer_message=str(fit.message), elapsed_seconds=perf_counter()-started)


def row_key(row):
    source, cohort, bounds = row['source_matching_counts'], row['cohort_counts'], row['mu_range']
    if (len(source) != 2 or len(cohort) != 3 or len(bounds) != 2
            or any(type(value) is not int for value in list(source)+list(cohort)+list(bounds))
            or tuple(source) not in source_counts() or cohort[2] != 0 or sum(cohort) != 17
            or not 0 <= cohort[0] <= 17
            or tuple(bounds) not in [(low, high) for low, high, _ in partitions()[cohort[1]]]):
        raise ValueError('Invalid sparse source/cohort/multiplicity cell')
    return tuple(source), cohort[0], tuple(bounds)


def expected_keys():
    return {(source, n1, (low, high)) for source in source_counts() for n1 in range(18)
            for low, high, _ in partitions()[17-n1]}


def load_controls(directory):
    data = {}
    for name, expected in INPUT_PINS.items():
        raw = (Path(directory)/name).read_bytes()
        if sha256(raw).hexdigest() != expected:
            raise ValueError(f'Pinned control changed: {name}')
        data[name] = json.loads(raw)
    baseline, sibling, hom256, hom384 = [data[name] for name in INPUT_PINS]
    if any(item['construction'] != CONSTRUCTION or item['active_groups'] != 17 for item in data.values()):
        raise ValueError('Control construction or fixed-r scope differs')
    if hom256['entries'] != hom384['entries']:
        raise ValueError('Directed homogeneous cell ceilings do not agree')
    selected = {}
    for source in baseline['sources']:
        for cohort in source['cohorts']:
            for old in cohort['entries']:
                row = dict(old, input_evidence='frozen-two-bank-FLOAT')
                key = row_key(row)
                if key in selected:
                    raise ValueError('Duplicate baseline cell')
                selected[key] = row
    if set(selected) != expected_keys():
        raise ValueError('Incomplete baseline sparse-cohort partition')
    added = []
    for old in sibling['screen']['entries']:
        added.append(dict(old, cohort_counts=[1, 16, 0], input_evidence='sibling-FLOAT-schedule'))
    for old in hom256['entries']:
        added.append(dict(old, cohort_counts=[0, 17, 0], input_evidence='directed-homogeneous-cell',
                          failure_log2_screen=float(F(old['failure_log2_upper']))))
    for row in added:
        key = row_key(row)
        if row['failure_log2_screen'] < selected[key]['failure_log2_screen']:
            selected[key] = row
    return selected


def summarize_rows(rows):
    mapped = {row_key(row): row for row in rows}
    if len(mapped) != len(rows) or set(mapped) != expected_keys():
        raise ValueError('Complete2673-cell partition required')
    cohorts = []
    for n1 in range(18):
        entries = [row for key, row in mapped.items() if key[1] == n1]
        value = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries]))
        worst = max(entries, key=lambda row: row['failure_log2_screen'])
        cohorts.append(dict(cohort_counts=[n1, 17-n1, 0], failure_log2_screen=value,
            worst_cell=worst, selected_cohort_reaches_50_bits=value <= -50))
    return dict(cohorts=cohorts,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts])),
        nonhomogeneous_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts[1:]])),
        preserved_homogeneous_mixed_log2_screen=float(np.logaddexp2.reduce(
            [float(HOMOGENEOUS_DIRECTED_CEILING)]+[row['failure_log2_screen'] for row in cohorts[1:]])))


def fixed_bank_screen(best, controls, *, progress=None):
    """One raw bank,81 source-conditioned cap banks, no witness search."""
    if set(controls) != expected_keys():
        raise ValueError('Complete frozen controls are required before any new table build')
    started = perf_counter()
    z, logp = F(best['mark']), F(best['log2_potential'])
    raw = pair_regions_screen(float(z), GROUPS)
    selected, new_rows = dict(controls), []
    for completed, (source, count) in enumerate(sorted(source_counts().items()), 1):
        evaluator = MultiplicityCapsFloat(ordered_region_tables(raw, GROUPS, *source), float(logp), GROUPS)
        caps = {upper: evaluator(upper) for _, upper in CANONICAL_BINS}
        source_rows = []
        for n1 in range(18):
            for low, high, geometries in partitions()[17-n1]:
                row = class_row(n1, (low, high), geometries, count, caps[canonical_mu(high)],
                                'joint-cohort-refit', z, logp)
                row.update(source_matching_counts=list(source), evaluated_cap_mu=canonical_mu(high),
                           input_evidence='new-fixed-bank-FLOAT')
                new_rows.append(row)
                source_rows.append(row)
                key = row_key(row)
                if row['failure_log2_screen'] < selected[key]['failure_log2_screen']:
                    selected[key] = row
        if progress:
            progress(dict(stage='fixed-bank-source', completed=completed, total=81,
                source_matching_counts=list(source),
                raw_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in source_rows])),
                combined_log2_screen=float(np.logaddexp2.reduce([
                    row['failure_log2_screen'] for key, row in selected.items() if key[0] == source])),
                elapsed_seconds=perf_counter()-started))
    rows = [selected[key] for key in sorted(selected)]
    return dict(raw_new_bank=dict(entries=new_rows, **summarize_rows(new_rows)),
        combined=dict(entries=rows, **summarize_rows(rows)),
        new_fixed_witness=dict(mark=str(z), log2_potential=str(logp)),
        input_sha256=dict(INPUT_PINS), complete_source_coverage=True,
        complete_sparse_cohort_coverage=True, source_classes=81, cohort_count=18,
        selected_cells=len(rows), new_raw_regional_banks=1, new_source_cap_banks=81,
        cap_mu_values=[upper for _, upper in CANONICAL_BINS],
        all_new_values_directed=False, combined_evidence='MIXED-PRELIMINARY',
        preserved_homogeneous_directed_ceiling=str(HOMOGENEOUS_DIRECTED_CEILING),
        fitted_new_witnesses_during_screen=0, elapsed_seconds=perf_counter()-started)


def metadata():
    return dict(construction=CONSTRUCTION, active_groups=17, field_order=Q, forbidden_weight=LIMIT,
        pair_swaps=True, fixed_geometry_across_columns=True, projective_normalizations=1,
        endpoint_potential_factors=1, directed=False, proposal_requires_fresh_directed_replay=True,
        whole_code_certificate=False)


def save(path, result):
    result = dict(metadata(), **result, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    with Path(path).open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-output', type=Path, required=True)
    parser.add_argument('--screen-output', type=Path, required=True)
    parser.add_argument('--controls-directory', type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    if args.fit_output.resolve() == args.screen_output.resolve() or args.fit_output.exists() or args.screen_output.exists():
        parser.error('Two distinct fresh output paths are required')
    controls = load_controls(args.controls_directory)
    progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    fit = fit_joint(progress=progress)
    save(args.fit_output, dict(fit=fit, source_matching_counts=list(FIT_SOURCE),
                              objective_cells=[dict(n1=n1, mu_range=list(bounds)) for n1, bounds in FIT_CELLS]))
    progress(dict(stage='fit-complete', **{key: value for key, value in fit.items() if key not in ('best', 'proposals')},
                  best=fit['best']))
    if fit['full_screen_authorized_by_gain']:
        result = fixed_bank_screen(fit['best'], controls, progress=progress)
        result['fit_artifact_sha256'] = sha256(args.fit_output.read_bytes()).hexdigest()
        save(args.screen_output, result)
        progress(dict(stage='screen-complete', elapsed_seconds=result['elapsed_seconds'],
            combined_log2_screen=result['combined']['failure_log2_screen'],
            preserved_homogeneous_mixed_log2_screen=result['combined']['preserved_homogeneous_mixed_log2_screen']))
